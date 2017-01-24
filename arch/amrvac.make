ifndef AMRVAC_DIR
$(error AMRVAC_DIR is not set)
endif

ARCH ?= default
NDIM := 1

# By exporting these can be used when building libamrvac
export ARCH NDIM

LIB_DIR := $(AMRVAC_DIR)/lib/$(NDIM)d_$(ARCH)
LIB_MAKE := $(AMRVAC_DIR)/arch/lib.make
LIB_AMRVAC := $(LIB_DIR)/libamrvac.a

# These are used for compilation
INC_DIRS := $(LIB_DIR)
LIB_DIRS := $(LIB_DIR)
LIBS := amrvac

.PHONY: all clean force

all: amrvac

# Include architecture and compilation rules
include $(AMRVAC_DIR)/arch/$(ARCH).defs
include $(AMRVAC_DIR)/arch/rules.make

# To get amrvac.f
vpath %.f $(LIB_DIR)

# Keep .f files for inspection
.PRECIOUS: %.f

# Intermediate files are removed
.INTERMEDIATE: amrvac.o mod_usr.o mod_usr.mod

# Always try to build/update the amrvac library
$(LIB_AMRVAC): force
	@mkdir -p $(LIB_DIR)
	$(MAKE) -C $(LIB_DIR) -f $(LIB_MAKE)

clean:
	$(RM) amrvac

# Dependencies
amrvac: mod_usr.o amrvac.o
$(LIB_DIR)/amrvac.f amrvac.o mod_usr.o: $(LIB_AMRVAC)
amrvac.o: mod_usr.o
