CC   = clang
LD = $(CC)

ifeq ($(ENABLE_OPENMP),true)
ifeq ($(shell uname -s),Darwin)
# Apple clang needs the homebrew libomp
LIBOMP_ROOT ?= /opt/homebrew/opt/libomp
OPENMP   = -Xpreprocessor -fopenmp
OMP_INC  = -I$(LIBOMP_ROOT)/include
OMP_LIB  = -L$(LIBOMP_ROOT)/lib -Wl,-rpath,$(LIBOMP_ROOT)/lib -lomp
else
OPENMP   = -fopenmp
endif
endif

VERSION  = --version
CFLAGS   = -O3 -ffast-math -std=c99 $(OPENMP)
#CFLAGS   = -Ofast -fnt-store=aggressive  -std=c99 $(OPENMP) #AMD CLANG
LFLAGS   = $(OPENMP)
DEFINES  = -D_GNU_SOURCE
INCLUDES = $(OMP_INC)
LIBS     = $(OMP_LIB)
