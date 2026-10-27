# parallel-reporting Specification

## Purpose

Defines what the MPI build prints and how it times phases, validates results and
terminates, so that scaling runs produce comparable and honest measurements.

## Requirements

### Requirement: Single output stream
In the MPI build only rank 0 SHALL print the header, validation results and timing
table. The header SHALL additionally report the number of ranks, the number of
shared-memory groups and ranks per group, the eigensolver backend and grid, the
transpose mode, and the memory per group.

#### Scenario: Output on several ranks
- **WHEN** the MPI build runs on 4 ranks
- **THEN** the header, validation lines and timing table appear exactly once

### Requirement: Phase timing under MPI
The MPI build SHALL report the phases `gto`, `transpose`, `overlap`, `reduce`, `eigen`,
`gather`, `rotate`, `untranspose` and `kinetic` in this order, and a `total` line. Each
phase time of a repetition SHALL be the maximum over all ranks. In the stored mode and
with a single group, every phase SHALL start after a barrier over all ranks. In the
streamed mode, transposes interleave with `overlap` and `rotate`; the per-rank time
spent in each of these phases SHALL be accumulated and its maximum over ranks reported.
As in the sequential build, the table SHALL list the minimum and average over
repetitions. Phases that do not occur in a run SHALL be reported with zero time.

#### Scenario: Timing table in an MPI run
- **WHEN** the MPI build runs with 2 groups
- **THEN** the timing table lists the nine phases and `total` with min and avg time
  columns

#### Scenario: Single group
- **WHEN** the MPI build runs with one shared-memory group
- **THEN** `transpose` and `untranspose` are reported with zero time

### Requirement: Rates for computation and communication phases
The timing table SHALL report a GFLOP/s rate for `overlap`, `rotate` and `kinetic`, based
on the global flop counts of the sequential build divided by the phase time. It SHALL
report an effective bandwidth in GB/s for `transpose`, `untranspose`, `reduce` and
`gather`, based on the bytes sent between shared-memory groups divided by the phase
time.

#### Scenario: Aggregate rate
- **WHEN** a run with 4 ranks reports the `overlap` phase
- **THEN** its GFLOP/s value equals n(n+1)·N divided by the reported minimum time,
  where n is the number of orbitals and N the number of grid points

### Requirement: Validation with reductions across ranks
With `-c`, the MPI build SHALL perform the same checks with the same tolerances as the
sequential build, combining per-rank results so that each check covers all orbitals and
grid points. The exit code SHALL be 1 on every rank if any check fails and 0 otherwise.

#### Scenario: Validation passes on several ranks
- **WHEN** the MPI build runs `POSCAR-0064-mp-149 -c` on 4 ranks
- **THEN** the output contains the GTO residual, normalization, kinetic energy and
  orthonormality checks, each with PASS
- **AND** `mpirun` returns exit code 0

#### Scenario: Validation fails on a coarse grid
- **WHEN** the MPI build runs `POSCAR-0008-mp-149 -g 0.8 -c` on 4 ranks
- **THEN** the normalization and kinetic energy checks report FAIL with the same values
  as the sequential build, and the job exits with code 1

### Requirement: Errors terminate the whole job
In the MPI build, any error that terminates the sequential program (invalid options,
unreadable input, unsupported cell, failed eigensolver, insufficient memory) SHALL
terminate all ranks with a non-zero exit code and print its message once.

#### Scenario: Missing input file
- **WHEN** the MPI build runs on 4 ranks with `-f does-not-exist`
- **THEN** the job terminates with a non-zero exit code and does not hang

### Requirement: LIKWID regions for all phases
With `ENABLE_LIKWID=true`, the MPI build SHALL provide a marker region for every
reported phase, named as in the timing table, on every rank.

#### Scenario: Marker regions in an MPI run
- **WHEN** an MPI build with LIKWID markers is run under `likwid-mpirun -m`
- **THEN** marker regions `gto`, `transpose`, `overlap`, `reduce`, `eigen`, `gather`,
  `rotate`, `untranspose` and `kinetic` are reported per rank
