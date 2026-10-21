CC   = clang
LD = $(CC)

VERSION  = --version
CFLAGS   = -O3 -ffast-math -std=c99
#CFLAGS   = -Ofast -fnt-store=aggressive  -std=c99 #AMD CLANG
LFLAGS   =
DEFINES  = -D_GNU_SOURCE
INCLUDES =
LIBS     =
