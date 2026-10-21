# Selects the BLAS/LAPACK(E) and FFTW libraries, see BLAS in config.mk
#   MKL      - Intel oneMKL (sequential layer), also provides the FFTW3 interface.
#              Requires MKLROOT.
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
LIBS     += -qmkl=sequential
else
LIBS     += -L$(MKLROOT)/lib -Wl,-rpath,$(MKLROOT)/lib -Wl,--no-as-needed \
            -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
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
