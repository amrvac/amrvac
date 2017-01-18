ifndef ARCH
	$(error build.make: ARCH is not set)
endif

ifndef NDIM
	$(error build.make: NDIM is not set)
endif

SRC_DIRS := . modules amrvacio physics rho hd #mhd
SRC_DIRS := $(addprefix $(AMRVAC_DIR)/src/, $(SRC_DIRS))

LIB := libamrvac_$(ARCH).a

.PHONY: all clean force
all: $(LIB)

.SUFFIXES:			# Disable built-in rules
.PRECIOUS: %.f			# Keep .f files

vpath %.t $(SRC_DIRS)		# Get .t files from SRC_DIRS

# Include makefiles, which define FOBJECTS and dependencies
include $(addsuffix /makefile, $(SRC_DIRS))

# Include architecture and rules
include $(AMRVAC_DIR)/arch/$(ARCH).defs
include $(AMRVAC_DIR)/arch/rules.make

OBJECTS := $(FOBJECTS:.t=.o) $(INCLUDES:.t=.o)
PPFLAGS := -z=-2 -phi=-1

# amrvac.f will be compiled later
$(LIB): $(OBJECTS) amrvac.f
	$(RM) $@
	$(AR) rcs $@ $^

# Used to rebuild objects when ARCH changes
current_arch: force
	@echo '$(ARCH)' | cmp -s - $@ || echo '$(ARCH)' > $@

clean:
	$(RM) $(OBJECTS) *.mod $(LIB)

# Dependencies
$(FOBJECTS:.t=.o): $(INCLUDES:.t=.o)
$(OBJECTS): current_arch
