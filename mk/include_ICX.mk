CC   = icx
LD = $(CC)

ifeq ($(ENABLE_LTO),true)
FAST_WORKAROUND = -ipo -O3 -fp-model=fast
else
FAST_WORKAROUND = -O3 -fp-model=fast
endif

VERSION  = --version
CFLAGS   = $(FAST_WORKAROUND) -xHost -std=c99 -Wno-unused-command-line-argument
LFLAGS   =
DEFINES  = -D_GNU_SOURCE
INCLUDES =
LIBS     =
