# Supported: GCC, CLANG, ICX
TOOLCHAIN ?= GCC
ENABLE_LIKWID ?= false
# Distributed-memory (MPI) build, binary dftbench-<TOOLCHAIN>-MPI
ENABLE_MPI ?= false
# MPI compiler wrapper (default: mpicc, mpiicx for TOOLCHAIN=ICX)
#MPICC ?= mpicc

# BLAS/LAPACK and FFT backend. Supported: MKL, OPENBLAS
BLAS ?= MKL
# Installation prefixes for BLAS=OPENBLAS (leave empty for system paths)
OPENBLAS_ROOT ?=
FFTW_ROOT ?=

# Eigensolver of the Loewdin eigen phase. Supported: LAPACK, SCALAPACK, ELPA
# SCALAPACK and ELPA require ENABLE_MPI=true.
EIGEN_SOLVER ?= LAPACK
# BLACS layer of MKL's ScaLAPACK for BLAS=MKL. Supported: openmpi, intelmpi
MKL_BLACS ?= openmpi
# Installation prefixes of ScaLAPACK (for BLAS=OPENBLAS) and ELPA (leave empty for system paths)
SCALAPACK_ROOT ?=
ELPA_ROOT ?=
# ELPA library name (default: elpa_openmp for an OpenMP build of ELPA, otherwise elpa)
#ELPA_LIB ?= elpa

#Feature options
OPTIONS +=  -DARRAY_ALIGNMENT=64
#OPTIONS +=  -DVERBOSE_TIMER
