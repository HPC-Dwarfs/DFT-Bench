# Selects the BLAS/LAPACK(E) and FFTW libraries, see BLAS in config.mk
#   MKL      - Intel oneMKL, also provides the FFTW3 interface. Requires MKLROOT.
#   OPENBLAS - any library providing cblas.h/lapacke.h linked with -lopenblas,
#              plus FFTW3. Set OPENBLAS_ROOT/FFTW_ROOT if not in default paths.

BLAS ?= MKL

ifeq ($(strip $(BLAS)),MKL)
ifeq ($(strip $(MKLROOT)),)
$(error BLAS=MKL but MKLROOT is not set. Source the oneAPI setvars.sh or set BLAS=OPENBLAS in config.mk)
endif
DEFINES  += -DUSE_MKL
INCLUDES += -I$(MKLROOT)/include -I$(MKLROOT)/include/fftw

ifeq ($(TOOLCHAIN),ICX)
ifeq ($(ENABLE_OPENMP),true)
LIBS     += -qmkl=parallel
else
LIBS     += -qmkl=sequential
endif
else
ifeq ($(ENABLE_OPENMP),true)
ifeq ($(TOOLCHAIN),GCC)
MKL_THREAD_LIB = -lmkl_gnu_thread -lgomp
else
MKL_THREAD_LIB = -lmkl_intel_thread -liomp5
endif
else
MKL_THREAD_LIB = -lmkl_sequential
endif
LIBS     += -L$(MKLROOT)/lib -Wl,-rpath,$(MKLROOT)/lib -Wl,--no-as-needed \
            -lmkl_intel_lp64 $(MKL_THREAD_LIB) -lmkl_core -lpthread -lm -ldl
endif

else ifeq ($(strip $(BLAS)),OPENBLAS)
ifneq ($(strip $(OPENBLAS_ROOT)),)
INCLUDES += -I$(OPENBLAS_ROOT)/include
LIBS     += -L$(OPENBLAS_ROOT)/lib -Wl,-rpath,$(OPENBLAS_ROOT)/lib
endif
ifneq ($(strip $(FFTW_ROOT)),)
INCLUDES += -I$(FFTW_ROOT)/include
LIBS     += -L$(FFTW_ROOT)/lib -Wl,-rpath,$(FFTW_ROOT)/lib
endif
LIBS     += -lfftw3 -lopenblas -lm

else
$(error Unsupported BLAS=$(BLAS), use MKL or OPENBLAS)
endif
