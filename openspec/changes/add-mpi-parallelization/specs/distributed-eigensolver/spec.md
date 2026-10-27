# Spec Delta

## Purpose

Defines the Löwdin eigen phase of the MPI build: which eigensolver backends compute
S^{-1/2}, on which ranks they run, and how their results are checked.

## ADDED Requirements

### Requirement: Löwdin orthogonalization semantics are kept
For every backend, the eigen phase SHALL compute all eigenpairs of the overlap matrix S
and the symmetric Löwdin matrix X = S^{-1/2}, which rotates the orbitals. The program
SHALL terminate with an error message and non-zero exit code if the eigensolver fails or
S is not positive definite, as the sequential build does.

#### Scenario: Orthonormal orbitals with a distributed backend
- **WHEN** an MPI build with `EIGEN_SOLVER=SCALAPACK` or `ELPA` runs
  `POSCAR-0064-mp-149` with `-c`
- **THEN** the orthonormality check max|S-I| passes its tolerance
- **AND** the printed overlap eigenvalue range agrees with the sequential build to a
  relative deviation of at most 1e-10

### Requirement: LAPACK backend under MPI
With `EIGEN_SOLVER=LAPACK` the MPI build SHALL compute the eigendecomposition on a
single rank and distribute X to all shared-memory groups, so that every rank uses an
identical X.

#### Scenario: LAPACK backend on several groups
- **WHEN** an MPI build with `EIGEN_SOLVER=LAPACK` runs with 4 ranks and `-s 2` and `-c`
- **THEN** the run passes validation

### Requirement: Eigen process grid
With a distributed backend, the eigen phase SHALL run on a 2D process grid that is a
subset of all ranks. The option `-e <pr>x<pc>` SHALL set the grid. Without `-e`, the
grid SHALL be square and sized so that each grid rank holds a local block of about
160 x 160 matrix elements, limited by the number of ranks. The grid SHALL be printed in
the header. Ranks outside the grid SHALL wait for the result.

#### Scenario: Explicit grid
- **WHEN** a run with 8 ranks uses `-e 2x3`
- **THEN** the header reports the eigen grid 2 x 3 and the run passes validation with `-c`

#### Scenario: Grid larger than the number of ranks
- **WHEN** `-e` requests more grid ranks than the job has
- **THEN** the program terminates with an error message and a non-zero exit code

#### Scenario: Grid option with the LAPACK backend
- **WHEN** `-e` is given in a build with `EIGEN_SOLVER=LAPACK`
- **THEN** the program prints a notice that the option has no effect and continues

### Requirement: Eigen rank placement
The option `-p packed|spread` SHALL select whether the eigen grid ranks are taken from as
few shared-memory groups as possible (`packed`, default) or distributed round-robin over
all groups (`spread`). The placement SHALL be printed in the header.

#### Scenario: Spread placement
- **WHEN** a run with 8 ranks, `-s 4` and `-e 2x2 -p spread` is started
- **THEN** the eigen grid uses 2 ranks of each of the 2 groups and the run passes
  validation with `-c`

### Requirement: Reference check of distributed backends
When validation (`-c`) is requested with a distributed backend and the number of
orbitals is at most 4096, the program SHALL additionally compute the eigendecomposition
with the LAPACK backend outside the timed region and SHALL report as checks the maximum
relative eigenvalue deviation and the maximum absolute deviation of X, each with a
tolerance of 1e-10. For larger matrices it SHALL print that the reference check is
skipped.

#### Scenario: Reference check passes
- **WHEN** an MPI build with `EIGEN_SOLVER=SCALAPACK` runs `POSCAR-0064-mp-149` with `-c`
- **THEN** the validation output contains PASS lines for the eigenvalue and X
  deviations against the LAPACK reference

#### Scenario: Reference check skipped for large matrices
- **WHEN** a distributed backend runs an input with more than 4096 orbitals with `-c`
- **THEN** the output states that the LAPACK reference check is skipped
