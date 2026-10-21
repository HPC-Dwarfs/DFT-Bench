# Distributed-memory build, see ENABLE_MPI in config.mk
#   MPICC - MPI compiler wrapper (default: mpicc, mpiicx for TOOLCHAIN=ICX)

ENABLE_MPI ?= false

ifeq ($(strip $(ENABLE_MPI)),true)
ifeq ($(TOOLCHAIN),ICX)
MPICC ?= mpiicx
else
MPICC ?= mpicc
endif
CC         = $(MPICC)
LD         = $(CC)
DEFINES   += -D_MPI
MPI_SUFFIX = -MPI
else ifneq ($(strip $(ENABLE_MPI)),false)
$(error Unsupported ENABLE_MPI=$(ENABLE_MPI), use true or false)
endif
