ifndef AMRVAC_DIR
$(error AMRVAC_DIR is not set)
endif

ARCH ?= default
NDIM := 2

# By exporting these can be used when building libamrvac
export ARCH NDIM

SRC_DIR := $(AMRVAC_DIR)/src
LIB_DIR := $(AMRVAC_DIR)/lib/$(NDIM)d_$(ARCH)
LIB_MAKE := $(AMRVAC_DIR)/arch/lib.make
LIB_AMRVAC := $(LIB_DIR)/libamrvac.a

# These are used for compilation
INC_DIRS := $(LIB_DIR)
LIB_DIRS := $(LIB_DIR)
LIBS := amrvac



.PHONY: all clean allclean force hdr



all: | hdr amrvac

# copy amrvac.h (in order to use the std preprocessor in the main code files, e.g. twofl); create the file if it does not exist
hdr:
ifeq ("$(wildcard amrvac.h)","")
	touch amrvac.h
endif
	@mkdir -p $(LIB_DIR)	# Prevent error message
	rsync -c amrvac.h $(LIB_DIR)/amrvac.h


# Include architecture and compilation rules
include $(AMRVAC_DIR)/arch/$(ARCH).defs
include $(AMRVAC_DIR)/arch/rules.make

# Optionally include a local user makefile
-include local.make

# Where to find amrvac.t
vpath %.t $(SRC_DIR)

# Keep mod_usr.f for inspection
.PRECIOUS: mod_usr.f

# Intermediate files are removed
.INTERMEDIATE: amrvac.o mod_usr.o mod_usr.mod

# Always try to build/update the amrvac library
$(LIB_AMRVAC): force
	@mkdir -p $(LIB_DIR)
	$(MAKE) -C $(LIB_DIR) -f $(LIB_MAKE)

clean:
	@echo 'Cleaning local objects ("make allclean" cleans libamrvac)'
	$(RM) amrvac amrvac.o mod_usr.o mod_usr.mod

# Also clean the library
allclean: clean
	@echo 'Cleaning libamrvac'
	@mkdir -p $(LIB_DIR)	# Prevent error message
	$(MAKE) -C $(LIB_DIR) -f $(LIB_MAKE) clean

# Dependencies
amrvac: mod_usr.o amrvac.o


amrvac.o mod_usr.o: $(LIB_AMRVAC)
amrvac.o: mod_usr.o
mod_usr.o: mod_usr.f amrvac.h
