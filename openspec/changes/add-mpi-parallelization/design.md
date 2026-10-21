# Design

## Context

See proposal.md for motivation. The current code (branch `sequential`, `032ee87`) holds
all orbitals in one column-major block `es.psi` (`ld x norb`, orbital i at
`psi + i * ld`), and `runDft` in `src/dft.c` runs five phases per repetition:

```
  phase     operation                           access pattern of Psi
  gto       putGtoSymOrthoOrb per orbital       whole orbital (column)
  overlap   cblas_dsyrk  S = dV Psi^T Psi       sum over grid rows
  eigen     LAPACKE_dsyevd + dgemm -> X         norb x norb only
  rotate    blocked dgemm Psi <- Psi X          row-local
  kinetic   FFTW r2c/c2r per orbital            whole orbital (column)
```

Constraints that shape the design:

- The grid spacing is fixed, so N = 2916 · norb at the default spacing for every
  POSCAR. All big phases and the eigen phase therefore scale as norb³, and the serial
  eigen share is about the same for every input size (measured: 0.4 % at 216 atoms).
- The overlap spectrum is fixed by the crystal structure (min 2.28e-2, max 1.98,
  κ ≈ 87 for all inputs), so all S-based orthogonalizations lose orthogonality
  only at about 1e-14.
- The largest input (8000 atoms) needs 5.4 TiB of orbitals. The replicated
  norb x norb matrices (S, X, work), about 5.7 GiB per rank, would take
  about 730 GiB per node at 128 ranks per node.
- The code is C99, formatted with `.clang-format`, and uses the `include_*.mk` pattern
  for build options. There is no test suite; correctness is checked with `-c`.

## Goals / Non-Goals

**Goals:**
- One code base for the sequential and the MPI build, with the numerics free of
  `#ifdef _MPI`.
- Multi-node runs that serve both strong scaling (64 to 1728 atoms) and the largest
  inputs, without special-casing either.
- One copy per node of the orbital data, S and X.
- The eigen phase stays a dense Löwdin eigensolve and is allowed to be the scaling
  bottleneck. The benchmark measures that bottleneck.
- The multi-node code paths can be tested on one machine.

**Non-Goals:**
- Overlapping communication with computation (pipelined transposes). Chunking makes
  this possible later, but the first version uses blocking collectives so phases
  can be attributed cleanly.
- Splitting a single orbital's FFT across ranks (a distributed FFT). `gto` and
  `kinetic` stop scaling once there are more ranks than orbitals. This is accepted
  and reported.
- Threading (OpenMP) or GPU offload.
- Cholesky or Newton–Schulz orthogonalization, and a Rayleigh–Ritz phase. These were
  discussed as possible later variants.
- Non-orthogonal cells, which the sequential code does not support either.

## Decisions

### D1. Communication layer with serial stubs

New `src/comm.h` and `src/comm.c` encapsulate all MPI calls. `comm.c` contains both
implementations under one `#ifdef _MPI`. In the serial version every collective is a
no-op or an identity, `commAllocateShared` calls `allocate()`, and the transposed form
is an alias of the orbital form. `dft.c` calls the layer unconditionally.

- *Alternative:* `#ifdef _MPI` inside `dft.c`. Rejected because it would interleave two
  programs in the kernels and duplicate timing and LIKWID code.
- *Alternative:* a separate MPI source tree or branch. Rejected by the user in favor
  of a build option.

### D2. Hierarchical decomposition: groups, ranks, two forms

```
  world
   +-- group g (MPI_Comm_split_type SHARED, optionally split into -s n ranks)
   |     ranks 0..L_g-1 share MPI-3 windows
   +-- lane l  (MPI_Comm_split by group-local rank, one rank per group)
   +-- leaders (lane 0)

  orbital form   (gto, kinetic)       transposed form  (overlap, rotate)
  group g: ld x norb_g column-major   group g: ldt_g x norb column-major
  rank r : contiguous orbital range   rank r : contiguous grid-row range
```

- **Orbitals** are assigned to groups in proportion to each group's rank count, then
  in contiguous ranges to ranks. Remainders go to the first ranks.
- **Grid rows** (the flat `ngtot` index) are split the same way. Row counts per rank
  are rounded to multiples of 8 so blocks stay 64-byte aligned. The split ignores
  the 3D geometry, because overlap and rotate do not use it.
- **Both forms are one contiguous column-major matrix per group.** Rank r's orbital
  range is contiguous in memory. Its row range is a strided block of the transposed
  matrix, so `rotate` and the tile products use the existing BLAS calls with
  `lda = ldt_g`.
- **With a single group**, the transposed form is the orbital form, with `ldt = ld`,
  all rows and all orbitals. No transpose and no second copy are needed. The
  sequential build is the single-group case, which is why it stays unchanged.
- **The `-s n` option** splits the shared communicator further, to emulate several
  nodes on one machine (needed for testing) or to experiment with groups per socket
  or NUMA domain.
- *Alternative:* a pure grid distribution with a distributed FFT (FFTW-MPI slabs).
  Rejected because slabs limit the rank count to `nz` (216 to 360), and pencil
  decompositions need a new library.
- *Alternative:* a pure orbital distribution with ring exchanges for overlap and
  rotate. Rejected because it moves all of Psi through every rank twice per
  repetition.

### D3. Shared windows and NUMA placement

The windows are allocated with `MPI_Win_allocate_shared`, with contiguous layout,
called on the group communicator. Each rank's segment is its own orbital columns,
padded to whole pages. The transposed window is allocated in one segment by group
rank 0. Base pointers come from `MPI_Win_shared_query`.

First touch decides NUMA placement:
- the orbital form is touched by its owners in `gto`;
- each rank zeroes its own row range of the transposed form right after allocation.

Synchronization: `MPI_Win_lock_all(MPI_MODE_NOCHECK)` once, then
`MPI_Win_sync` + `MPI_Barrier(group)` at every phase boundary where ranks read
other ranks' data. The phases are already bulk-synchronous.

S, X and the eigenvalues are small group windows allocated by group rank 0.

For FFTW, `kinetic` plans on rank-private buffers. Before each new-array execute it
checks `fftw_alignment_of` of the orbital column, and falls back to a copy into the
private buffer if the alignment differs. Columns are 64-byte aligned by construction
(`ld` padding plus page-aligned segments), so the copy should never happen.

### D4. Transposes between groups

Only data that crosses groups goes through MPI:
- **Rows within the own group:** each rank copies its own orbitals' rows into the
  group's transposed window (`memcpy` from the shared pointers).
- **Rows for other groups:** one `MPI_Alltoallv` over `world`, with zero counts
  inside the group. Rank p (group a) sends to rank q (group b ≠ a) the rows of q's
  range of p's orbitals, packed column by column. Rank q unpacks them into columns
  `o_p...` of its rows.
- The untranspose does the reverse.

Pack and unpack buffers are sized per chunk. MPI derived datatypes, which would avoid
the packing, are a later optimization.

### D5. Memory budget: stored and streamed mode

```
  stored   : orbital form + transposed form + chunk buffers        (2 transposes)
  streamed : orbital form + transposed chunk + chunk buffers       (3 transposes)

  streamed, per repetition:
    pass 1: for chunk k: transpose-in(k)  -> accumulate S tiles (beta = 1)
    reduce, eigen, gather
    pass 2: for chunk k: transpose-in(k)  -> rotate(k) -> untranspose(k)
            (overwrites rows of chunk k in the orbital form: rotate of chunk k
             reads only old rows of chunk k, so the update is safe in place)
```

A chunk covers the same fraction of every rank's row range. The default budget is
80 % of the node's physical memory (`sysconf(_SC_PHYS_PAGES)`) divided by the groups
on the node; `-m` overrides it. The program picks stored mode if it fits, otherwise
the smallest chunk count K that fits, and stops with an error if the orbital form
alone exceeds the budget.

- *Alternative:* always streamed. Rejected because it costs 1.5x the communication in
  strong-scaling runs, where memory is plentiful.
- *Alternative:* always stored. Rejected because it needs about twice the nodes for
  the 8000-atom input.

### D6. Overlap as tiles within a group, reduced per lane

```
  group-shared transposed Psi ----> upper-triangle tiles of S (OVERLAP_TILE, e.g. 512)
                                    assigned round-robin to the group's ranks
                                    diagonal tile: dsyrk, off-diagonal: dgemm,
                                    over all rows of the group (lda = ldt_g)
  barrier(group)
  S (group window) split into segments l = 0..Lmin-1 (Lmin = min ranks per group)
  lane l: MPI_Allreduce(MPI_IN_PLACE, segment l, SUM) across groups
  barrier(group)
```

- No rank holds a private full-size partial S.
- The inter-group reduction runs over Lmin independent lanes in parallel, so it
  does not serialize through one leader per node.
- Separating tile computation (all group ranks) from the reduction segments
  (Lmin lanes) allows groups with different rank counts.
- In streamed mode the tiles accumulate over the chunks and are reduced once.
- The lower triangle is filled in by symmetry wherever a full matrix is needed (the
  eigen backends and the validation check).
- *Alternative:* a per-rank `dsyrk` followed by `MPI_Allreduce` of the full S.
  Rejected because it brings back the norb² memory per rank (2 GB at norb = 16000)
  and sends the full matrix from every rank.

### D7. Eigen backends behind one interface

New `src/eigen.h` with one entry point: it takes the group-shared S and returns the
eigenvalues and the group-shared X. The backend is chosen at compile time
(`EIGEN_SOLVER`) and lives in `eigenLapack.c`, `eigenScalapack.c` and `eigenElpa.c`.
Only the selected file is compiled.

- **LAPACK:** global rank 0 runs the existing `loewdinMatrix` code on its group's S.
  X and the eigenvalues are broadcast over the leader communicator into every
  group's window, then a group barrier follows. The sequential build uses the same
  file with no broadcast. All ranks get a bitwise-identical X.
- **SCALAPACK:**
  - Grid: square by default, q = clamp(floor(norb / 512), 1, floor(sqrt(P))), so each
    rank holds a local block of about 512 x 512. `-e pr x pc` overrides it.
    Placement follows `-p` (packed: the first pr·pc ranks in world order; spread:
    round-robin over groups). A BLACS context is built from the selected ranks
    with `blacs_gridmap`.
  - Input: grid ranks copy their block-cyclic blocks (block size `EIGEN_NB`, e.g. 64)
    straight from the group-shared S. No communication is needed, because every
    group holds the full S after D6.
  - Solve: `pdsyevd` gives U and the eigenvalues. Each rank scales its local columns
    to form W = U · diag(eval^-1/4), then `pdsyrk` computes X = W Wᵀ. This is n³
    flops; the current code uses a 2 n³ `dgemm`.
  - Output: the result is redistributed with `pdgemr2d` onto the 1x1 grid of global
    rank 0 and then broadcast over the leaders, the same path as the LAPACK backend.
    A direct all-gather into every group window is a later optimization once timings
    show the gather matters.
- **ELPA:** the same grid, blocks and data path as SCALAPACK, with the eigensolve
  done by `elpa_eigenvectors` (2-stage solver by default) instead of `pdsyevd`.
  ELPA uses the ScaLAPACK descriptors, so it shares all the code except the solver
  call.
- **Reference check (`-c`, norb ≤ 4096):** global rank 0 copies S before the solve,
  runs the LAPACK path outside the timers, and compares the eigenvalues and X.
- *Alternatives considered:* Newton–Schulz iteration (all `dgemm`, no dependency,
  scales best, but has no eigensolver character) and Cholesky (cheapest, but not
  Löwdin). The user chose to keep a real eigensolver. SLATE, DLA-Future and EigenExa
  were rejected as heavier dependencies.

### D8. Timing and phases

- `phaseStart` and `phaseStop` in `dft.c` are kept. In the MPI build `phaseStart`
  first calls `commBarrierWorld`. At the end of each repetition the per-rank phase
  times are combined with one `MPI_Reduce(MAX)` into rank 0.
- Streamed mode: within the passes, transposes and computation alternate per chunk.
  Their times are accumulated into the `transpose`, `overlap`, `rotate` and
  `untranspose` slots with a barrier only before each pass.
- The phase enum grows to nine entries. The sequential build prints only the
  original five, in the old order, so its output is unchanged.
- Communication volume (bytes sent between groups) is counted per phase by the comm
  layer and turned into GB/s.
- LIKWID regions are registered for every phase that is printed.

### D9. Startup, input, output and errors

- `main` calls `commInit` before option parsing and `commFinalize` at exit.
- Every rank parses the options and reads the POSCAR file, which is small and gives
  identical results everywhere. Only rank 0 prints.
- Errors go through a new `dftAbort(msg)`: in the sequential build it prints and
  calls `exit(EXIT_FAILURE)`; in the MPI build it prints and calls
  `MPI_Abort(world, 1)`. The `exit(EXIT_FAILURE)` calls in `dft.c` are replaced.
  `poscarReader.c` and `gtoOnGrid.c` keep `exit()`, since a failure there happens on
  every rank, and `mpirun` ends the job. Hangs are avoided because none of those
  paths is inside a collective.
- The MPI-only options `-s`, `-m`, `-e` and `-p` are parsed only in the MPI build, so
  the sequential usage text does not change.

### D10. Build integration

- **`mk/include_MPI.mk`** (new), included after the tool chain file. If
  `ENABLE_MPI=true`:
  - `CC = LD = $(MPICC)`, where `MPICC ?= mpicc`, or `mpiicx` for ICX;
  - `DEFINES += -D_MPI`;
  - `MPI_SUFFIX = -MPI`.
- **`Makefile`:** `TARGET = dftbench-$(TOOLCHAIN)$(MPI_SUFFIX)` and
  `BUILD_DIR = ./build/$(TOOLCHAIN)$(MPI_SUFFIX)`. Both are recursive variables, so
  they pick up the suffix after the includes. The object rule also depends on the
  new `.mk` files.
- **`mk/include_EIGEN.mk`** (new):
  - `EIGEN_SOLVER ?= LAPACK` adds `-DEIGEN_LAPACK`, the default.
  - `SCALAPACK` with `BLAS=MKL` links `-lmkl_scalapack_lp64` and
    `-lmkl_blacs_$(MKL_BLACS)_lp64`, where `MKL_BLACS ?= openmpi` (alternative
    `intelmpi`).
  - `SCALAPACK` otherwise links `-L$(SCALAPACK_ROOT)/lib -lscalapack`.
  - `ELPA` adds `ELPA_ROOT` includes and `-lelpa` on top of the ScaLAPACK flags.
  - Any distributed backend without `ENABLE_MPI=true` stops with `$(error ...)`.
- **The eigen source files** are selected with `EIGEN_SRC` in `include_EIGEN.mk`. The
  Makefile's `wildcard` source list excludes the backends that are not chosen.
- **`mk/config-default.mk`** documents `ENABLE_MPI`, `MPICC`, `EIGEN_SOLVER`,
  `MKL_BLACS`, `SCALAPACK_ROOT` and `ELPA_ROOT`. All variables are `?=`, so an older
  `config.mk` keeps working.

### D11. Validation under MPI

Each check is computed locally and combined:
- **GTO residual and normalization:** computed by the orbital owners in the orbital
  form, combined with `MPI_MAX`.
- **Kinetic energy:** a per-rank partial sum, combined with `MPI_SUM`.
- **Orthonormality:** the distributed overlap is recomputed (transpose, tiles,
  reduce) and checked on rank 0 against the full S.
- **Pass flag:** combined with a logical AND reduction.

The tolerances in `constants.h` stay as they are.

## Risks / Trade-offs

- **[Many small Alltoallv messages in strong scaling]** For 1728 atoms on 3456 ranks a
  message pair is about 23 KB, less with chunks. → Measure first. A later step can
  aggregate per group (leaders, or a 2-level all-to-all). The comm layer hides the
  change.
- **[The serial eigen phase limits scaling with the LAPACK backend]** This is the
  documented and intended behavior of the fallback. With ScaLAPACK or ELPA it is
  still the phase that runs into its limit first. → Report it, and use `-e` and `-p`
  to explore it.
- **[Gathering X through rank 0 costs O(norb²) through one node]** At norb = 16000
  that's 2 GB, about 0.1 to 0.5 s against minutes for overlap and rotate. → Accept
  for now. A direct all-gather is a contained follow-up.
- **[Heuristic defaults for the eigen grid and tile sizes are untuned]** → An early
  task measures `dsyevd` against `pdsyevd` and ELPA at norb 1024 and 3456, to
  calibrate the default grid and `EIGEN_NB`.
- **[Different results across groups if X differed]** → X is computed once and
  broadcast with every backend.
- **[Shared windows on systems with restricted shared memory]** Some containers or
  MPI setups disallow it. → `-s 1` gives one rank per group, which works anywhere, at
  the cost of per-rank copies.
- **[NUMA effects on large nodes]** → First touch by the owning rank (D3), and `-s`
  to form groups per socket or NUMA domain if whole-node windows turn out slower.
- **[`exit()` in shared kernels under MPI]** → Only reached on every rank alike,
  outside collectives (D9).
- **[Toolchain and MPI combinations]** MKL BLACS must match the MPI library. →
  `MKL_BLACS`, documented in the README.

## Migration Plan

- All work is on branch `mpi`, created from `sequential` at `032ee87`. `sequential`
  stays untouched.
- The default build (`ENABLE_MPI=false`) is the regression reference. At every step
  `dftbench-<TOOLCHAIN> -c` must still pass and print the same results.
- A rollback is the branch itself: nothing lands on `sequential` or `main` as part of
  this change.

## Open Questions

- The exact default values of `OVERLAP_TILE` and `EIGEN_NB`, and the constant in
  the eigen grid heuristic. They are tuned by the spike in the tasks and don't affect
  the specs.
- Whether `spread` or `packed` should become the default placement once
  measurements exist. Specs fix `packed` as the default for now.
