# Proposal

## Why

DFT-Bench is sequential: it runs on one core and its working set is bounded by one
node's memory. The README already reserves the 2744 to 8000 atom inputs (several TB
of orbitals) for a distributed-memory version, and the mid-sized inputs (64 to 1728
atoms) are natural strong-scaling cases. A multi-node MPI version serves both regimes
with one code path, and it keeps the dense eigensolver as the scaling bottleneck that
real DFT codes face, without hiding it.

## What Changes

- New build option `ENABLE_MPI` (default `false`). The existing sequential build keeps
  its dependencies, options, phases and results, and serves as the reference for the
  MPI build.
- New build option `EIGEN_SOLVER` = `LAPACK` (default) | `SCALAPACK` | `ELPA` that
  selects the backend of the Löwdin eigen phase. `SCALAPACK` and `ELPA` require
  `ENABLE_MPI=true`.
- Multi-node distribution of the orbital block with a hierarchical layout:
  - orbital distribution (whole orbitals per node and rank) for `gto` and `kinetic`,
  - transposed grid-row distribution (all orbitals at a range of grid points) for
    `overlap` and `rotate`,
  - node-level MPI-3 shared-memory windows, so that the orbital block, the overlap
    matrix S and the Löwdin matrix X exist once per node instead of once per rank,
  - inter-node transposes only, done in chunks under a per-node memory budget
    (stored mode: 2 transposes per repetition; streamed mode: 3 transposes and about
    half the orbital memory).
- The overlap matrix is computed as node-local output tiles and reduced tile-wise
  across nodes; no rank holds a private full-size partial S.
- The Löwdin orthogonalization is unchanged in meaning (S^{-1/2} via an
  eigendecomposition). With `SCALAPACK` or `ELPA` it runs on a sub-communicator sized
  to the matrix on a 2D BLACS process grid. The `LAPACK` backend computes on one rank
  and broadcasts X; it is also the validation reference.
- New command line options: memory budget per node, eigen process grid, and eigen
  rank placement.
- Reporting: new timed phases for the communication steps (`transpose`,
  `untranspose`, `reduce`, `gather`), phase times as the maximum over ranks behind a
  barrier, rank and node counts in the header, and a notice when ranks exceed the
  number of orbitals (idle ranks in `gto`/`kinetic`).
- Validation (`-c`) works unchanged under MPI with reductions across ranks. On request
  it also compares the distributed eigensolver against the LAPACK reference.
- README documents the MPI build, the new options and the scaling limits.

## Capabilities

### New Capabilities
- `mpi-build`: build configuration for the MPI version (`ENABLE_MPI`, MPI compiler
  wrappers per tool chain, `EIGEN_SOLVER` backends and their libraries) and the
  guarantee that the sequential build is preserved.
- `distributed-orbitals`: how the orbital block is distributed over nodes and ranks,
  the node-shared windows, the chunked transposes under a memory budget, and the
  distributed overlap and rotation.
- `distributed-eigensolver`: the Löwdin eigen phase on a sub-communicator with the
  LAPACK, ScaLAPACK and ELPA backends, its process grid and placement options, and the
  distribution of X back to the nodes.
- `parallel-reporting`: phase timing and output under MPI, the new communication
  phases, and validation with reductions across ranks.

### Modified Capabilities
<!-- None: openspec/specs/ is empty, the sequential behavior is not yet captured in
     specs. The sequential build's behavior is required to stay unchanged
     (see mpi-build). -->

## Impact

- **Code**: `src/dft.c` (phases operate on local ranges, call the comm layer),
  `src/main.c` (MPI init/finalize, new options, rank-0 output), `src/dft.h`
  (configuration and distribution data). New `src/comm.[ch]` (MPI and serial stub
  implementations) and `src/eigen*.[ch]` (eigensolver backends).
  `src/gtoOnGrid.c`, `src/poscarReader.c` and the FFT code stay sequential kernels.
- **Build**: `Makefile`, `mk/config-default.mk`, new `mk/include_MPI.mk` and
  `mk/include_EIGEN.mk`, tool chain files for MPI compiler wrappers.
- **Dependencies**: MPI-3 (shared-memory windows) when `ENABLE_MPI=true`; ScaLAPACK
  plus BLACS (MKL or netlib) for `SCALAPACK`; ELPA plus ScaLAPACK for `ELPA`. No new
  dependency for the default build.
- **Output**: additional phases and header lines. Scripts that parse the output need to
  accept the new rows. The sequential build prints the same phases as before.
- **Branch**: all work lands on branch `mpi`. Branch `sequential` stays untouched.
