# Supported: GCC, CLANG, ICX
TOOLCHAIN ?= GCC
ENABLE_OPENMP ?= false
ENABLE_LIKWID ?= false

# BLAS/LAPACK and FFT backend. Supported: MKL, OPENBLAS
BLAS ?= MKL
# Installation prefixes for BLAS=OPENBLAS (leave empty for system paths)
OPENBLAS_ROOT ?=
FFTW_ROOT ?=

#Feature options
OPTIONS +=  -DARRAY_ALIGNMENT=64
#OPTIONS +=  -DVERBOSE_AFFINITY
#OPTIONS +=  -DVERBOSE_TIMER
