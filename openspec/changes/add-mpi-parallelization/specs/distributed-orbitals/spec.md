# Spec Delta

## Purpose

Defines how the MPI build distributes the orbital block across nodes and ranks, how much
memory it uses per node, and that results do not depend on the decomposition.

## ADDED Requirements

### Requirement: Results are independent of the decomposition
The MPI build SHALL produce, for any number of ranks and shared-memory groups, the same
validation outcome as the sequential build, and the printed kinetic energy of the
orthonormalized orbitals and the overlap eigenvalue range SHALL agree with the
sequential build to a relative deviation of at most 1e-10.

#### Scenario: Validation across rank counts
- **WHEN** `dftbench-<TOOLCHAIN>-MPI -f poscars/POSCAR-0064-mp-149 -c` is run with 1, 2, 3,
  4 and 7 ranks
- **THEN** every run reports `Validation PASSED`
- **AND** the printed `Ekin (orthonormalized orbitals)` agrees with the sequential build
  to a relative deviation of at most 1e-10

#### Scenario: Uneven distribution
- **WHEN** the number of orbitals is not divisible by the number of ranks or groups
- **THEN** the run completes and passes validation

### Requirement: More ranks than orbitals
The MPI build SHALL run when the number of ranks exceeds the number of orbitals. It
SHALL print a notice stating how many ranks own no orbitals and therefore idle in the
`gto` and `kinetic` phases.

#### Scenario: Small input on many ranks
- **WHEN** `POSCAR-0008-mp-149` (16 orbitals) is run on 20 ranks with `-c`
- **THEN** the run passes validation
- **AND** the output contains a notice that 4 ranks own no orbitals

### Requirement: Shared-memory groups
Ranks SHALL be organized in shared-memory groups that hold one copy of the orbital data,
the overlap matrix and the Löwdin matrix per group. By default a group SHALL be all ranks
of one compute node. The option `-s <n>` SHALL limit groups to at most `n` ranks; groups
SHALL never span compute nodes. The number of groups and ranks per group SHALL be
printed in the header.

#### Scenario: Default grouping
- **WHEN** the MPI build runs with 4 ranks on one compute node without `-s`
- **THEN** the header reports 1 shared-memory group with 4 ranks

#### Scenario: Emulating several nodes on one machine
- **WHEN** the MPI build runs with 6 ranks on one compute node with `-s 2`
- **THEN** the header reports 3 shared-memory groups with 2 ranks each
- **AND** the run behaves as on 3 nodes, including inter-group transposes and
  reductions, and passes validation with `-c`

#### Scenario: Invalid group size
- **WHEN** `-s` is given a value smaller than 1
- **THEN** the program terminates with an error message and a non-zero exit code

### Requirement: Memory budget and transpose mode
The MPI build SHALL keep the per-group orbital memory within a budget. The default
budget SHALL be 80 percent of the physical memory of the node divided by the number of
groups on that node; the option `-m <GB>` SHALL set the budget per group in GB. If the
orbital data and a full transposed copy fit within the budget, the program SHALL use the
stored mode (transposed copy kept between overlap and rotation). Otherwise it SHALL use
the streamed mode with the smallest number of chunks that fits. If the orbital data
alone does not fit, the program SHALL terminate with an error stating the required
memory per group. With a single group, no transposed copy SHALL be allocated. The
chosen mode, the number of chunks and the memory per group SHALL be printed in the
header.

#### Scenario: Stored mode
- **WHEN** a run with 2 groups has a budget that holds the orbital data and the
  transposed copy
- **THEN** the header reports the stored transpose mode

#### Scenario: Streamed mode forced by a small budget
- **WHEN** a run with 2 groups is given `-m` with a budget below the stored-mode
  requirement but above the orbital data alone
- **THEN** the header reports the streamed transpose mode and its number of chunks
- **AND** the run passes validation with `-c`

#### Scenario: Budget too small
- **WHEN** `-m` is below the memory needed for the orbital data of one group
- **THEN** the program terminates before allocating the orbitals with an error stating
  the required memory per group and a non-zero exit code

#### Scenario: Single group
- **WHEN** all ranks form one shared-memory group
- **THEN** the run performs no transposes and allocates no transposed copy
