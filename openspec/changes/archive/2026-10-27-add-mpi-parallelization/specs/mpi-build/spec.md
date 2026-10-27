# Spec Delta

## Purpose

Defines how DFT-Bench is built with or without MPI and with which eigensolver backend,
so that the sequential benchmark stays available unchanged next to the distributed one.

## ADDED Requirements

### Requirement: MPI is an opt-in build option
The build SHALL provide an `ENABLE_MPI` option in `config.mk` with default `false`.
With `ENABLE_MPI=false` the build SHALL NOT require MPI headers, libraries or compiler
wrappers, and the resulting binary SHALL accept the same command line options and print
the same phases as the sequential version.

#### Scenario: Default build stays sequential
- **WHEN** the user builds with an unmodified `mk/config-default.mk`
- **THEN** the binary `dftbench-<TOOLCHAIN>` is built without any MPI dependency
- **AND** `./dftbench-<TOOLCHAIN> -f poscars/POSCAR-0064-mp-149 -c` passes validation and
  reports the phases `gto`, `overlap`, `eigen`, `rotate` and `kinetic`

#### Scenario: Existing config.mk without the new options
- **WHEN** a `config.mk` created before this change, which does not define `ENABLE_MPI`
  or `EIGEN_SOLVER`, is used
- **THEN** the build behaves as with `ENABLE_MPI=false` and `EIGEN_SOLVER=LAPACK`

### Requirement: MPI build uses the tool chain's MPI compiler wrapper
With `ENABLE_MPI=true` the build SHALL compile and link with the MPI compiler wrapper
that matches the selected tool chain (`mpicc` for GCC and CLANG, `mpiicx` for ICX by
default). The wrapper SHALL be overridable with an `MPICC` variable.

#### Scenario: MPI build with GCC
- **WHEN** the user builds with `TOOLCHAIN=GCC ENABLE_MPI=true`
- **THEN** all sources are compiled and linked with `mpicc` and the MPI code paths are
  enabled

#### Scenario: Custom wrapper
- **WHEN** the user builds with `ENABLE_MPI=true MPICC=/opt/mpi/bin/mpicc`
- **THEN** that wrapper is used for compiling and linking

### Requirement: MPI and sequential builds coexist
The MPI build SHALL use a binary name `dftbench-<TOOLCHAIN>-MPI` and a separate build
directory `build/<TOOLCHAIN>-MPI`, so that sequential and MPI builds of the same tool
chain can coexist in one checkout without sharing object files.

#### Scenario: Both builds in one checkout
- **WHEN** the user runs `make` with `ENABLE_MPI=false` and then with `ENABLE_MPI=true`
- **THEN** both `dftbench-GCC` and `dftbench-GCC-MPI` exist
- **AND** the second build recompiles all objects rather than reusing sequential ones

### Requirement: Eigensolver backend selection
The build SHALL provide an `EIGEN_SOLVER` option with the values `LAPACK` (default),
`SCALAPACK` and `ELPA`. `SCALAPACK` SHALL link ScaLAPACK and BLACS, from MKL when
`BLAS=MKL` and from a configurable installation prefix otherwise. `ELPA` SHALL link ELPA
from a configurable installation prefix in addition to ScaLAPACK. The selected backend
SHALL be printed in the program header.

#### Scenario: ScaLAPACK with MKL
- **WHEN** the user builds with `BLAS=MKL ENABLE_MPI=true EIGEN_SOLVER=SCALAPACK`
- **THEN** the binary links MKL's ScaLAPACK and the BLACS layer matching the MPI library
- **AND** the program header reports the eigensolver backend `SCALAPACK`

#### Scenario: Distributed backend without MPI
- **WHEN** the user builds with `ENABLE_MPI=false EIGEN_SOLVER=SCALAPACK` or
  `EIGEN_SOLVER=ELPA`
- **THEN** `make` stops with an error stating that the backend requires `ENABLE_MPI=true`

#### Scenario: Unknown backend
- **WHEN** the user builds with an `EIGEN_SOLVER` value other than `LAPACK`,
  `SCALAPACK` or `ELPA`
- **THEN** `make` stops with an error listing the supported values
