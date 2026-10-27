# Backend of the Loewdin eigen phase, see EIGEN_SOLVER in config.mk
#   LAPACK    - dsyevd on one rank, X is broadcast to all ranks (default)
#   SCALAPACK - pdsyevd on a 2D BLACS grid, requires ENABLE_MPI=true.
#               BLAS=MKL links MKL's ScaLAPACK and BLACS (MKL_BLACS: openmpi, intelmpi),
#               otherwise -lscalapack from SCALAPACK_ROOT (empty: system paths).
#   ELPA      - ELPA 2-stage solver on the ScaLAPACK grid, requires ENABLE_MPI=true.
#               Links ELPA from ELPA_ROOT in addition to ScaLAPACK.
# eigenLapack.c is always compiled, it is the reference of the distributed backends.

EIGEN_SOLVER   ?= LAPACK
MKL_BLACS      ?= openmpi
SCALAPACK_ROOT ?=
ELPA_ROOT      ?=

EIGEN_SOURCES = eigenLapack.c eigenScalapack.c eigenElpa.c

ifeq ($(strip $(EIGEN_SOLVER)),LAPACK)
DEFINES  += -DEIGEN_LAPACK
EIGEN_SRC = eigenLapack.c

else ifneq ($(filter $(strip $(EIGEN_SOLVER)),SCALAPACK ELPA),)
ifneq ($(strip $(ENABLE_MPI)),true)
$(error EIGEN_SOLVER=$(EIGEN_SOLVER) requires ENABLE_MPI=true)
endif

ifeq ($(strip $(BLAS)),MKL)
SCALAPACK_LIBS = -lmkl_scalapack_lp64 -lmkl_blacs_$(MKL_BLACS)_lp64
else ifneq ($(strip $(SCALAPACK_ROOT)),)
SCALAPACK_LIBS = -L$(SCALAPACK_ROOT)/lib -Wl,-rpath,$(SCALAPACK_ROOT)/lib -lscalapack
else
SCALAPACK_LIBS = -lscalapack
endif

ifeq ($(strip $(EIGEN_SOLVER)),SCALAPACK)
DEFINES  += -DEIGEN_SCALAPACK
EIGEN_SRC = eigenLapack.c eigenScalapack.c
# ScaLAPACK and BLACS must precede the BLAS/LAPACK libraries
LIBS     := $(SCALAPACK_LIBS) $(LIBS)
else
DEFINES  += -DEIGEN_SCALAPACK -DEIGEN_ELPA
EIGEN_SRC = eigenLapack.c eigenScalapack.c eigenElpa.c
ifneq ($(strip $(ELPA_ROOT)),)
ELPA_INC  = $(firstword $(wildcard $(ELPA_ROOT)/include/elpa*))
INCLUDES += $(addprefix -I,$(ELPA_INC))
# An OpenMP build of ELPA installs include/elpa_openmp-<version> and libelpa_openmp
ELPA_LIB ?= $(if $(findstring elpa_openmp,$(ELPA_INC)),elpa_openmp,elpa)
ELPA_LIBS = -L$(ELPA_ROOT)/lib -Wl,-rpath,$(ELPA_ROOT)/lib -l$(ELPA_LIB)
else
ELPA_LIB ?= elpa
ELPA_LIBS = -l$(ELPA_LIB)
endif
LIBS     := $(ELPA_LIBS) $(SCALAPACK_LIBS) $(LIBS)
endif

else
$(error Unsupported EIGEN_SOLVER=$(EIGEN_SOLVER), use LAPACK, SCALAPACK or ELPA)
endif

EIGEN_EXCLUDE = $(filter-out $(EIGEN_SRC),$(EIGEN_SOURCES))
