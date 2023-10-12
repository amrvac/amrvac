ifndef GMUNU_DIR
$(error GMUNU_DIR is not set)
endif

ARCH ?= default
NDIM := 1

# By exporting these can be used when building libgmunu
export ARCH NDIM

SRC_DIR := $(GMUNU_DIR)/src
LIB_DIR := $(GMUNU_DIR)/lib/$(NDIM)d_$(ARCH)
LIB_MAKE := $(GMUNU_DIR)/arch/lib.make
LIB_GMUNU := $(LIB_DIR)/libgmunu.a

# These are used for compilation
INC_DIRS := $(LIB_DIR)
LIB_DIRS := $(LIB_DIR)
LIBS := gmunu

.PHONY: all clean allclean force

all: gmunu

# Include architecture and compilation rules
include $(GMUNU_DIR)/arch/$(ARCH).defs
include $(GMUNU_DIR)/arch/rules.make

# Optionally include a local user makefile
-include local.make

# Where to find gmunu.t
vpath %.t $(SRC_DIR)

# Keep mod_usr.f for inspection
.PRECIOUS: mod_usr.f

# Intermediate files are removed
.INTERMEDIATE: gmunu.o mod_usr.o mod_usr.mod

# Always try to build/update the gmunu library
$(LIB_GMUNU): force
	@mkdir -p $(LIB_DIR)
	$(MAKE) -C $(LIB_DIR) -f $(LIB_MAKE)

clean:
	@echo 'Cleaning local objects ("make allclean" cleans libgmunu)'
	$(RM) gmunu *.o *.o *.mod

# Also clean the library
allclean: clean
	@echo 'Cleaning libgmunu'
	@mkdir -p $(LIB_DIR)	# Prevent error message
	$(MAKE) -C $(LIB_DIR) -f $(LIB_MAKE) clean

# Dependencies
gmunu: mod_usr.o gmunu.o
gmunu.o mod_usr.o: $(LIB_GMUNU)
gmunu.o: mod_usr.o
