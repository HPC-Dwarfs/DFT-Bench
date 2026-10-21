# Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
# All rights reserved. This file is part of DFT-Bench.
# Use of this source code is governed by a MIT-style
# license that can be found in the LICENSE file.

#CONFIGURE BUILD SYSTEM
TARGET	   = dftbench-$(TOOLCHAIN)$(MPI_SUFFIX)
BUILD_DIR  = ./build/$(TOOLCHAIN)$(MPI_SUFFIX)
SRC_DIR    = ./src
MAKE_DIR   = ./mk
Q         ?= @

#DO NOT EDIT BELOW
ifeq (,$(wildcard config.mk))
$(info )
$(info ====================================================================)
$(info config.mk does not exist!)
$(info Creating config.mk from ./mk/config-default.mk)
$(info Please adapt config.mk to your needs and run make again.)
$(info ====================================================================)
$(info )
$(shell cp ./mk/config-default.mk config.mk)
$(error Stopping after creating config.mk - please review and run make again)
endif
include config.mk
include $(MAKE_DIR)/include_$(TOOLCHAIN).mk
include $(MAKE_DIR)/include_MPI.mk
include $(MAKE_DIR)/include_LIKWID.mk
include $(MAKE_DIR)/include_BLAS.mk
include $(MAKE_DIR)/include_EIGEN.mk
CFLAGS    += -Wall -Wextra

VPATH     = $(SRC_DIR)
CSRC      = $(filter-out $(addprefix $(SRC_DIR)/,$(EIGEN_EXCLUDE)),$(wildcard $(SRC_DIR)/*.c))
ASM       = $(patsubst $(SRC_DIR)/%.c, $(BUILD_DIR)/%.s,$(CSRC))
OBJ       = $(patsubst $(SRC_DIR)/%.c, $(BUILD_DIR)/%.o,$(CSRC))
SRC       =  $(wildcard $(SRC_DIR)/*.h $(SRC_DIR)/*.c)
CPPFLAGS := $(CPPFLAGS) $(DEFINES) $(OPTIONS) $(INCLUDES)
c := ,
clist = $(subst $(eval) ,$c,$(strip $1))

define CLANGD_TEMPLATE
CompileFlags:
  Add: [$(call clist,$(CPPFLAGS)), $(call clist,$(CFLAGS)), -xc]
  Compiler: clang
endef

${TARGET}: $(BUILD_DIR) .clangd $(OBJ)
	$(info ===>  LINKING  $(TARGET))
	$(Q)${LD} ${LFLAGS} -o $(TARGET) $(OBJ) $(LIBS)

MK_DEPS   = $(addprefix $(MAKE_DIR)/,include_$(TOOLCHAIN).mk include_MPI.mk include_BLAS.mk include_EIGEN.mk) config.mk

# Recompile when the flags change, e.g. EIGEN_SOLVER given on the command line
FLAGS_STAMP = $(BUILD_DIR)/.flags
FLAGS_LINE  = $(CC) $(CPPFLAGS) $(CFLAGS)
ifneq ($(FLAGS_LINE),$(shell cat $(FLAGS_STAMP) 2>/dev/null))
$(shell mkdir -p $(BUILD_DIR) && echo '$(FLAGS_LINE)' > $(FLAGS_STAMP))
endif
MK_DEPS  += $(FLAGS_STAMP)

$(BUILD_DIR)/%.o:  %.c $(MK_DEPS)
	$(info ===>  COMPILE  $@)
	$(Q)$(CC) -c $(CPPFLAGS) $(CFLAGS) $< -o $@
	$(Q)$(CC) $(CPPFLAGS) -MT $(@:.d=.o) -MM  $< > $(BUILD_DIR)/$*.d

$(BUILD_DIR)/%.s:  %.c
	$(info ===>  GENERATE ASM  $@)
	$(Q)$(CC) -S $(CPPFLAGS) $(CFLAGS) $< -o $@

.PHONY: clean distclean info asm format

clean:
	$(info ===>  CLEAN)
	@rm -rf $(BUILD_DIR)

distclean:
	$(info ===>  DIST CLEAN)
	@rm -rf build
	@rm -f dftbench-*
	@rm -f .clangd compile_commands.json

info:
	$(info $(CFLAGS))
	$(Q)$(CC) $(VERSION)

asm:  $(BUILD_DIR) $(ASM)

format:
	@for src in $(SRC) ; do \
		echo "Formatting $$src" ; \
		clang-format -i $$src ; \
	done
	@echo "Done"

$(BUILD_DIR):
	@mkdir -p $(BUILD_DIR)

.clangd:
	$(file > .clangd,$(CLANGD_TEMPLATE))

-include $(OBJ:.o=.d)
