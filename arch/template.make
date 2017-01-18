ifndef AMRVAC_DIR
$(error AMRVAC_DIR is not set)
endif

ARCH ?= default
NDIM := 1

# By exporting these can be used when building libamrvac
export ARCH NDIM

BUILD_DIR := $(AMRVAC_DIR)/builds/$(NDIM)d_$(ARCH)
BUILD_MAKEFILE := $(AMRVAC_DIR)/arch/build.make
LIB_AMRVAC := $(BUILD_DIR)/libamrvac.a

INC_DIRS := $(BUILD_DIR)
LIB_DIRS := $(BUILD_DIR)
LIBS := amrvac

.PHONY: all clean force

all: amrvac

# Include architecture and compilation rules
include $(AMRVAC_DIR)/arch/$(ARCH).defs
include $(AMRVAC_DIR)/arch/rules.make

# To get amrvac.f from BUILD_DIR
vpath %.f $(BUILD_DIR)

# Always try to build/update the amrvac library
$(LIB_AMRVAC): force
	@mkdir -p $(BUILD_DIR)
	$(MAKE) -C $(BUILD_DIR) -f $(BUILD_MAKEFILE)

clean:
	$(RM) mod_usr.o mod_usr.mod amrvac.o amrvac

# Dependencies
amrvac: mod_usr.o amrvac.o
amrvac.o: mod_usr.o $(LIB_AMRVAC)
mod_usr.o: $(LIB_AMRVAC)
