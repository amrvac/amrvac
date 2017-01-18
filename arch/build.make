ifndef ARCH
$(error build.make: ARCH is not set)
endif

ifndef NDIM
$(error build.make: NDIM is not set)
endif

SRC_DIRS := . modules amrvacio physics rho hd #mhd
SRC_DIRS := $(addprefix $(AMRVAC_DIR)/src/, $(SRC_DIRS))
LIB_AMRVAC := libamrvac.a
PPFLAGS := -z=-2 -phi=-1	# Remove in future

.PHONY: libamrvac clean

# Ensure that LIB_AMRVAC is the default target
libamrvac: $(LIB_AMRVAC)

# Include makefiles, which define FOBJECTS and dependencies
include $(addsuffix /makefile, $(SRC_DIRS))

# Include architecture and rules
include $(AMRVAC_DIR)/arch/$(ARCH).defs
include $(AMRVAC_DIR)/arch/rules.make

# Get .t files from SRC_DIRS
vpath %.t $(SRC_DIRS)

OBJECTS := $(FOBJECTS:.t=.o) $(INCLUDES:.t=.o)

# amrvac.f depends on user code so it will be compiled later
$(LIB_AMRVAC): $(OBJECTS) amrvac.f
	$(RM) $@
	$(AR) rcs $@ $^

clean:
	$(RM) $(OBJECTS) *.mod $(LIB_AMRVAC)

# INCLUDES are always compiled before FOBJECTS
$(FOBJECTS:.t=.o): $(INCLUDES:.t=.o)
