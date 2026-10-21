# Tasks

## 1. Baseline and build options

- [x] 1.1 Record sequential reference outputs (`-c -r 1` for POSCAR-0008, -0064, -0216, and `POSCAR-0008 -g 0.8 -c`) from the unchanged build into the scratchpad. Verify: the files exist and contain `Validation PASSED`, or `FAILED` for `-g 0.8`.
- [x] 1.2 Add `mk/include_MPI.mk` (`ENABLE_MPI`, `MPICC` default per tool chain, `-D_MPI`, `MPI_SUFFIX`) and make `TARGET`/`BUILD_DIR` in `Makefile` use the suffix. Verify: `make` and `make ENABLE_MPI=true` produce `dftbench-CLANG` and `dftbench-CLANG-MPI` in separate build dirs.
- [x] 1.3 Add `mk/include_EIGEN.mk` with `EIGEN_SOLVER` (LAPACK default, SCALAPACK, ELPA), the library flags, the backend source selection and the `$(error)` checks. Verify: `make EIGEN_SOLVER=SCALAPACK` without MPI and `make EIGEN_SOLVER=FOO` both stop with the specified messages.
- [x] 1.4 Document the new variables in `mk/config-default.mk` using `?=`. Verify: a build with the existing local `config.mk`, which lacks the variables, still produces the sequential binary.

## 2. Communication layer skeleton

- [x] 2.1 Create `src/comm.h`/`src/comm.c` with serial stubs and the MPI implementation of init/finalize, rank and size queries, rank-0 printing, barriers, reductions and `dftAbort`. Call them from `main.c` and replace `exit(EXIT_FAILURE)` in `dft.c` with `dftAbort`. Verify: sequential output matches 1.1 apart from timings, and `mpirun -np 4 dftbench-CLANG-MPI -c` prints the header once.
- [x] 2.2 Build the group, lane and leader communicators (`MPI_Comm_split_type` plus the `-s` split) and print groups and ranks per group in the header. Verify: `-np 4` reports 1 group of 4, `-np 6 -s 2` reports 3 groups of 2, and `-s 0` exits with an error.
- [x] 2.3 Implement shared window allocation (contiguous, page-padded segments, `MPI_Win_shared_query`, `lock_all`) and the group sync helper (`MPI_Win_sync` + barrier), with the serial fallback to `allocate()`. Verify: an MPI run of POSCAR-0008 on 4 ranks allocates the orbital window and passes `-c` once task 3.3 is done. Until then, check that the window base pointers on all group ranks see each other's writes, using a temporary debug print removed before completion.
- [x] 2.4 Parse the MPI-only options `-s`, `-m`, `-e`, `-p` in the MPI build only, with usage text and error handling. Verify: the sequential `-h` output is unchanged, the MPI `-h` lists the options, and invalid values exit non-zero with one message.

## 3. Decomposition within a single group

- [x] 3.1 Add the distribution data (group and rank orbital ranges, grid-row ranges rounded to multiples of 8) to `dft.h` and compute it in `runDft`. Print the notice for ranks without orbitals. Verify: POSCAR-0008 on 20 ranks prints that 4 ranks own no orbitals.
- [x] 3.2 Let `gto` and `kinetic` work on the rank's orbital range in the group orbital window, with the FFTW alignment check (D3), and reduce `validateOrbitals` results with MAX/SUM. Verify: on 1 to 4 ranks with one group, the GTO, normalization and Ekin check values match 1.1.
- [x] 3.3 Implement the overlap tiles (dsyrk/dgemm over the group's rows, round-robin tile assignment, symmetric fill), the row-range `rotate` on the aliased transposed form, and the LAPACK eigen backend in `eigenLapack.c` behind `eigen.h`. Verify: `-c` passes with 1, 2, 3, 4 and 7 ranks in one group, and Ekin and the eigenvalue range match 1.1 to 1e-10 relative.

## 4. Multiple groups

- [x] 4.1 Implement the stored-mode transpose and untranspose: intra-group copy plus `MPI_Alltoallv` with packing, with inter-group byte counters. Verify: `-np 4 -s 2 -c` and `-np 6 -s 2 -c` pass, and Ekin matches 1.1 to 1e-10.
- [x] 4.2 Implement the lane-wise `MPI_Allreduce` of S segments over `Lmin` lanes. Verify: `-np 7 -s 3 -c` (groups of 3, 3, 1) passes, and so do `-np 5 -s 2 -c` and `-np 8 -s 4 -c`.
- [x] 4.3 Broadcast X and the eigenvalues from global rank 0 over the leaders into all group windows (LAPACK backend). Verify: multi-group runs from 4.1 and 4.2 still pass, and the printed eigenvalue range matches 1.1.
- [x] 4.4 Implement the memory budget: default from physical memory and groups per node, `-m`, mode and chunk selection, error before allocation, and header lines. Verify: the header shows stored mode by default for `-np 4 -s 2` on POSCAR-0216, and a too-small `-m` exits with the required memory per group.
- [x] 4.5 Implement the streamed mode: pass 1 transposes and accumulates S, pass 2 transposes, rotates and untransposes in place. Verify: POSCAR-0216 with `-np 4 -s 2` and a `-m` value that forces streamed mode with at least 3 chunks passes `-c`, and Ekin matches 1.1 to 1e-10.

## 5. Reporting

- [x] 5.1 Extend the phases to the nine MPI phases, with a barrier in `phaseStart`, accumulated timing in streamed mode, and a per-repetition MAX reduction. Keep the original five phases in the sequential build. Verify: the MPI table lists nine phases plus `total`, single-group runs show zero transpose times, and sequential output matches 1.1 in structure.
- [x] 5.2 Report GFLOP/s from global flop counts for `overlap`, `rotate` and `kinetic`, and GB/s from the inter-group byte counters for the communication phases. Verify: the `overlap` GFLOP/s equals n(n+1)·N / min time for a 4-rank run, and the communication phases show non-zero GB/s with `-s 2`.
- [ ] 5.3 Register LIKWID marker regions for all printed phases. Verify: `make ENABLE_MPI=true ENABLE_LIKWID=true` compiles on a system with LIKWID, and `likwid-mpirun -m` lists the nine regions per rank.
- [x] 5.4 Check error termination under MPI. Verify: `mpirun -np 4 dftbench-CLANG-MPI -f does-not-exist` exits non-zero without hanging, and the same holds for an invalid `-e` grid.

## 6. ScaLAPACK backend

- [x] 6.1 Provide ScaLAPACK for development: MKL on the cluster, or netlib ScaLAPACK for OpenBLAS (e.g. `brew install scalapack`, after user approval). Wire up `SCALAPACK_ROOT`/`MKL_BLACS`. Verify: `make ENABLE_MPI=true EIGEN_SOLVER=SCALAPACK` links and the header reports the backend `SCALAPACK`.
- [x] 6.2 Implement eigen grid selection (default heuristic, `-e`, `-p packed|spread`) and the BLACS context via `blacs_gridmap`. Print grid and placement, and print the notice for `-e` with LAPACK. Verify: `-np 8 -e 2x3` reports a 2 x 3 grid, `-np 8 -s 4 -e 2x2 -p spread` uses 2 ranks per group, and `-e 4x4` on 8 ranks exits with an error.
- [x] 6.3 Implement `eigenScalapack.c`: block-cyclic copy from the group-shared S, `pdsyevd`, column scaling, `pdsyrk`, `pdgemr2d` to rank 0, and the leader broadcast of X. Verify: POSCAR-0064 `-c` passes on 4 and 8 ranks, with and without `-s`, and the eigenvalue range matches 1.1 to 1e-10.
- [x] 6.4 Implement the LAPACK reference check under `-c` for norb ≤ 4096, run outside the timers, and the skip message above that. Verify: POSCAR-0064 shows PASS lines for the eigenvalue and X deviations, and POSCAR-2744 with `-g 1.0 -c -r 1` prints that the reference check is skipped.

## 7. ELPA backend

- [ ] 7.1 Implement `eigenElpa.c` (ELPA setup with the ScaLAPACK descriptors, 2-stage solver, the rest shared with 6.3) and the `ELPA_ROOT` build flags. Verify: on a system with ELPA, POSCAR-0064 `-c` passes on 4 and 8 ranks, including the reference check.

## 8. Calibration

- [ ] 8.1 On a cluster node, time the eigen phase for the LAPACK, SCALAPACK and ELPA backends at norb 1024 and 3456 for several grid sizes. Use the results to set the defaults of `OVERLAP_TILE`, `EIGEN_NB` and the grid heuristic constant, and record them in design.md. Verify: the table of measurements is in design.md, and the defaults in `constants.h` match it.

## 9. Documentation and final verification

- [x] 9.1 Update README.md: MPI build, `EIGEN_SOLVER` and library variables, the new options `-s -m -e -p`, the phase table under MPI, stored and streamed modes, and the scaling limits (ranks beyond the orbital count, eigen phase). Verify: README sections exist, and every example command in them runs.
- [x] 9.2 Run the full regression matrix: the sequential build (GCC and CLANG) against 1.1; MPI with 1, 2, 3, 4 and 7 ranks; `-s 1`, `-s 2`, `-s 3`; stored and streamed modes; the LAPACK and SCALAPACK backends; `-g 0.8 -c` exit code 1. Run `make format`. Verify: all runs pass, or fail exactly as specified, and `git diff` after `make format` is empty.
