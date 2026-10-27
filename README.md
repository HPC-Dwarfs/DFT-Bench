# DFT-Bench

A planewave quantum chemistry inspired benchmark application.

DFT-Bench mimics the core operations of a planewave / real-space DFT code on a
set of trial orbitals for an atomistic structure read from a VASP POSCAR file:

1. **gto** – place two normalized Gaussian orbitals per atom on a periodic real-space
   grid (separable 1D Gaussians, folded back into the periodic cell).
2. **overlap** – compute the overlap matrix S = ΔV · Ψᵀ Ψ (BLAS `dsyrk`).
3. **eigen** – diagonalize S and build S^{-1/2} (LAPACK `dsyevd` + `dgemm`).
4. **rotate** – Löwdin orthogonalization Ψ ← Ψ · S^{-1/2} (blocked BLAS `dgemm`).
5. **kinetic** – kinetic energy ⟨ψ|−½∇²|ψ⟩ of all orbitals via a 3D real-to-complex
   FFT, multiplication by ½|k|² and a complex-to-real FFT (FFTW).

The grid spacing is fixed, so the number of grid points N grows linearly with the
number of atoms. The cost of the overlap and rotation phases is therefore O(N_atoms³),
as in production DFT codes.

The default build is sequential: it runs on a single core and contains no threading.
The MPI build (`ENABLE_MPI=true`) distributes the orbitals over many nodes, see
[MPI version](#mpi-version).

## Getting Started

To build and run **DFT-Bench**, you need a C compiler, GNU make, FFTW3, and a
BLAS/LAPACK library with CBLAS and LAPACKE interfaces. Intel oneMKL provides both
BLAS/LAPACK and the FFTW3 interface.

1. **Install a supported compiler**
   - GCC
   - Clang
   - Intel ICX

2. **Clone the repository**

   ```sh
   git clone https://github.com/HPC-Dwarfs/DFT-Bench.git

   cd DFT-Bench
   ```

3. **Adjust configuration**

   The first `make` creates `config.mk` from `mk/config-default.mk`. Edit it to select
   the tool chain, the BLAS backend, etc.

4. **Build**

   ```sh
   make
   ```

5. **Run**

   ```sh
   ./dftbench-<TOOLCHAIN> -f poscars/POSCAR-0064-mp-149 -c
   ```

## Build

Configure the tool chain and additional options in `config.mk`:

```make
# Supported: GCC, CLANG, ICX
TOOLCHAIN ?= GCC
ENABLE_LIKWID ?= false

# BLAS/LAPACK and FFT backend. Supported: MKL, OPENBLAS
BLAS ?= MKL
# Installation prefixes for BLAS=OPENBLAS (leave empty for system paths)
OPENBLAS_ROOT ?=
FFTW_ROOT ?=

#Feature options
OPTIONS +=  -DARRAY_ALIGNMENT=64
#OPTIONS +=  -DVERBOSE_TIMER
```

- `BLAS=MKL` requires `MKLROOT` to be set (e.g. by sourcing oneAPI `setvars.sh`) and
  uses MKL's FFTW3 interface and its sequential layer.
- `BLAS=OPENBLAS` links `-lopenblas -lfftw3`. Any BLAS library that provides `cblas.h`
  and `lapacke.h` can be used this way, e.g. with a symlink or by adjusting
  `mk/include_BLAS.mk`. If the OpenBLAS library is multi-threaded, run with
  `OPENBLAS_NUM_THREADS=1`. Example for Homebrew on macOS:

  ```sh
  make TOOLCHAIN=CLANG BLAS=OPENBLAS \
       OPENBLAS_ROOT=/opt/homebrew/opt/openblas FFTW_ROOT=/opt/homebrew/opt/fftw
  ```

- `ENABLE_LIKWID=true` enables LIKWID marker regions named after the phases above
  (`gto`, `overlap`, `eigen`, `rotate`, `kinetic`):

  ```sh
  likwid-perfctr -C 0 -g FLOPS_DP -m ./dftbench-GCC -f poscars/POSCAR-0216-mp-149
  ```

- `-DVERBOSE_TIMER` prints the timer resolution.

You can build multiple tool chains in the same directory, but notice that the
Makefile is only acting on the one currently set. Intermediate build results are
located in the `./build/<TOOLCHAIN>` directory.

Clean up intermediate build results for the active tool chain with `make clean`, and
for all tool chains with `make distclean`. `make asm` generates assembler files in
`./build/<TOOLCHAIN>`, `make format` reformats all sources with `clang-format`.

## Usage

```
Usage: ./dftbench-<TOOLCHAIN> [options]
  -f <file>     POSCAR input file (default: poscars/POSCAR-0008-mp-149)
  -g <h>        target grid spacing in bohr (default: 0.30)
  -n <x,y,z>    explicit grid size, overrides -g
  -r <reps>     number of timed repetitions (default: 3)
  -c            run validation checks (exit code 1 on failure)
  -v            verbose output
  -h            print this help
```

The grid size per direction is the smallest even, FFT friendly (2, 3, 5 smooth) number
of points that achieves the requested spacing. Only orthogonal cells are supported.

The program prints the minimum and average time of every phase over all repetitions
and a nominal GFLOP/s rate. The rate is based on n(n+1)·N flops for the overlap,
2·N·n² for the rotation and 5·N·log₂N + 4·N per orbital for the kinetic energy, where
n is the number of orbitals and N the number of grid points.

### Validation

With `-c` the first repetition additionally checks:

- the placed orbitals against an independently computed periodic Gaussian,
- the orbital normalization,
- the numerical kinetic energy against the analytical value 3/(2·gw²) per orbital
  (fails if the grid is too coarse for the Gaussians),
- the orthonormality max|S − I| after the Löwdin orthogonalization.

The exit code is 1 if any check fails.

### Inputs

`poscars/` contains silicon diamond supercells (Materials Project mp-149). Grid size
and orbital memory at the default spacing of 0.3 bohr:

| Atoms | Orbitals | Grid  | Orbital memory |
| ----: | -------: | ----- | -------------: |
| 8     | 16       | 36³   | 0.01 GB        |
| 64    | 128      | 72³   | 0.36 GB        |
| 216   | 432      | 108³  | 4.05 GB        |
| 512   | 1024     | 144³  | 22.8 GB        |
| 1000  | 2000     | 180³  | 86.9 GB        |
| 1728  | 3456     | 216³  | 259 GB         |

The larger inputs (2744 to 8000 atoms) need several TB at this spacing. They are
intended for the MPI version or can be run with a coarser grid (`-g`), which fails
the kinetic-energy validation.

## MPI version

### Build

```make
# Distributed-memory (MPI) build, binary dftbench-<TOOLCHAIN>-MPI
ENABLE_MPI ?= false
# MPI compiler wrapper (default: mpicc, mpiicx for TOOLCHAIN=ICX)
#MPICC ?= mpicc

# Eigensolver of the Loewdin eigen phase. Supported: LAPACK, SCALAPACK, ELPA
# SCALAPACK and ELPA require ENABLE_MPI=true.
EIGEN_SOLVER ?= LAPACK
# BLACS layer of MKL's ScaLAPACK for BLAS=MKL. Supported: openmpi, intelmpi
MKL_BLACS ?= openmpi
# Installation prefixes of ScaLAPACK (for BLAS=OPENBLAS) and ELPA (leave empty for system paths)
SCALAPACK_ROOT ?=
ELPA_ROOT ?=
```

- `ENABLE_MPI=true` compiles and links with the MPI compiler wrapper of the tool chain
  (`MPICC` overrides it) and needs MPI-3 shared-memory windows. The binary is
  `dftbench-<TOOLCHAIN>-MPI`, its objects go to `./build/<TOOLCHAIN>-MPI`, so the
  sequential and the MPI build of a tool chain coexist.
- `EIGEN_SOLVER` selects the backend of the eigen phase:
  - `LAPACK` (default): `dsyevd` on rank 0, X is broadcast to all nodes.
  - `SCALAPACK`: `pdsyevd` on a 2D BLACS process grid. `BLAS=MKL` links MKL's
    ScaLAPACK and the BLACS layer `MKL_BLACS`, which must match the MPI library.
    Otherwise `-lscalapack` is linked from `SCALAPACK_ROOT`.
  - `ELPA`: the ELPA 2-stage solver on the same grid, linked from `ELPA_ROOT` in
    addition to ScaLAPACK. An OpenMP build of ELPA (`include/elpa_openmp-*`) is
    linked as `-lelpa_openmp`; `ELPA_LIB` overrides the library name.

Examples (command line variables override `config.mk`), the second one for Homebrew
on macOS, whose ScaLAPACK is built with Homebrew's Open MPI:

```sh
make TOOLCHAIN=CLANG BLAS=OPENBLAS ENABLE_MPI=true
make TOOLCHAIN=CLANG BLAS=OPENBLAS ENABLE_MPI=true EIGEN_SOLVER=SCALAPACK \
     MPICC=/opt/homebrew/bin/mpicc SCALAPACK_ROOT=/opt/homebrew/opt/scalapack
```

On an Intel MPI system with oneAPI and an ELPA module that sets `ELPA_ROOT`, e.g. Fritz
at NHR@FAU:

```sh
module load intel mkl intelmpi elpa
make TOOLCHAIN=ICX BLAS=MKL ENABLE_MPI=true EIGEN_SOLVER=ELPA MKL_BLACS=intelmpi
```

The ScaLAPACK library must be built for the same MPI library as `MPICC`, and the
program must be started with that library's `mpirun`. With `BLAS=MKL`, set
`MKL_BLACS=intelmpi` for Intel MPI.

### Usage

The MPI build accepts the options of the sequential build and additionally:

```
MPI options:
  -s <n>        at most n ranks per shared-memory group (default: node)
  -m <GB>       memory budget per group (default: 80% of node memory)
  -e <pr>x<pc>  eigensolver process grid (default: about 160^2 elements per rank)
  -p <mode>     eigensolver rank placement: packed or spread (default: packed)
```

```sh
mpirun -np 4 ./dftbench-CLANG-MPI -f poscars/POSCAR-0064-mp-149 -c
mpirun -np 8 ./dftbench-CLANG-MPI -f poscars/POSCAR-0064-mp-149 -c -s 4
mpirun -np 8 ./dftbench-CLANG-MPI -f poscars/POSCAR-0064-mp-149 -c -s 4 -m 0.25
```

Only rank 0 prints. The header additionally reports the ranks, the shared-memory
groups, the transpose mode, the memory budget and use per group and the eigensolver.

### Decomposition

The ranks of a node form a shared-memory group, which holds one copy of the orbitals,
of S and of X in MPI-3 shared windows. `-s <n>` splits the node into groups of at most
n ranks, e.g. one group per socket, or several groups on one machine to test the
multi-node code paths: `-s 2` on 6 ranks behaves like 3 nodes.

Orbitals and grid points are split over the groups in proportion to their ranks. The
orbitals are held in two forms:

- **orbital form**, whole orbitals per group and rank, used by `gto` and `kinetic`,
- **transposed form**, all orbitals at the group's grid points, used by `overlap` and
  `rotate`.

With a single group both forms are the same memory and no data is moved between
ranks. With several groups the `transpose` and `untranspose` phases move the rows of
other groups with `MPI_Alltoallv`. The overlap tiles of a group are summed over the
groups (`reduce`), the eigen phase computes X, and `gather` distributes it to all
groups. Each group needs the orbitals of its share (orbital form), the transposed form
and S, X:

- **stored mode**: the transposed form is kept from `overlap` to `rotate`, two
  transposes per repetition.
- **streamed mode**: if the stored mode does not fit into the memory budget, the
  transposed form is processed in chunks of grid rows. It needs about half the memory
  and three transposes per repetition.

The budget defaults to 80 % of the node memory divided by the groups on the node,
`-m` sets it per group in GB. The program picks the stored mode if it fits, otherwise
the streamed mode with the fewest chunks that fit, and stops if the orbitals of a
group alone exceed the budget.

### Eigensolver grid

With `SCALAPACK` or `ELPA` the eigen phase runs on a `pr x pc` subset of the ranks.
By default the grid is square, with about 160 x 160 matrix elements per grid rank and
at most all ranks, e.g. 6 x 6 for 1024 orbitals on 36 or more ranks. `-e pr x pc` sets
it, `-p packed` takes the grid ranks from as few groups as possible and `-p spread`
round-robin from all groups. The other ranks wait for the result. With `LAPACK`, `-e`
and `-p` have no effect. ELPA needs at least one `EIGEN_NB` x `EIGEN_NB` block (32 x 32)
per grid row and column, so with `ELPA` a grid with 8 rows or columns needs at least
225 orbitals; the program stops with an error for a grid that is too large.

```sh
mpirun -np 8 ./dftbench-CLANG-MPI -f poscars/POSCAR-0064-mp-149 -c -e 2x3
mpirun -np 8 ./dftbench-CLANG-MPI -f poscars/POSCAR-0064-mp-149 -c -s 4 -e 2x2 -p spread
```

With `-c` and at most 4096 orbitals, the distributed backends are additionally
checked against the LAPACK path (eigenvalues and X, tolerance 1e-10), outside the
timed phases.

### Timing

The MPI build reports nine phases:

| Phase         | Operation                                                    |
| ------------- | ------------------------------------------------------------ |
| `gto`         | place the orbitals of each rank                              |
| `transpose`   | orbital form to transposed form, between groups              |
| `overlap`     | tiles of S over the group's grid points (`dsyrk`, `dgemm`)   |
| `reduce`      | sum S over the groups (`MPI_Allreduce`, in parallel lanes)   |
| `eigen`       | eigendecomposition and X = S^{-1/2}                          |
| `gather`      | distribute X to all groups                                   |
| `rotate`      | Ψ ← Ψ · X on the rows of each rank                           |
| `untranspose` | transposed form back to orbital form                         |
| `kinetic`     | kinetic energy of the orbitals of each rank                  |

Every phase starts after a barrier, and its time is the maximum over the ranks. In
the streamed mode the transposes alternate with `overlap` and `rotate` per chunk; the
time per phase is accumulated. Phases that do not occur are reported with zero time,
e.g. `transpose` with a single group. GFLOP/s is computed from the same global flop
counts as in the sequential build; GB/s from the bytes moved between groups. With
`ENABLE_LIKWID=true` every phase is a marker region on every rank:

```sh
likwid-mpirun -np 4 -g FLOPS_DP -m ./dftbench-GCC-MPI -f poscars/POSCAR-0216-mp-149
```

### Scaling limits

- `gto` and `kinetic` distribute whole orbitals. Ranks beyond the number of orbitals
  idle in these phases; the program prints how many.
- `overlap` distributes tiles of 128 x 128 elements of S over the ranks of a group,
  n/128 (n/128 + 1) / 2 tiles for n orbitals. Ranks beyond the number of tiles idle
  in this phase, e.g. 92 of 128 ranks for 1024 orbitals.
- The eigen phase is a dense eigensolve of n x n with n orbitals. With `LAPACK` it
  runs on one rank, with `SCALAPACK` or `ELPA` on the eigen grid. It is the phase
  that stops scaling first, which is intended: real DFT codes face the same limit.
