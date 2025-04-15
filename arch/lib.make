ifndef AMRVAC_DIR
$(error AMRVAC_DIR is not set)
endif

ifndef ARCH
$(error build.make: ARCH is not set)
endif

ifndef NDIM
$(error build.make: NDIM is not set)
endif

SRC_DIRS := . modules amr io limiter physics hd particle
SRC_DIRS := $(addprefix $(AMRVAC_DIR)/src/, $(SRC_DIRS))
LIB_AMRVAC := libamrvac.a
PPFLAGS :=

.PHONY: libamrvac clean
.PRECIOUS: %.F90			# Don't remove intermediate .F90 files

libamrvac: $(LIB_AMRVAC) amrvac.F90

# Include makefiles, which define FOBJECTS and dependencies
include $(addsuffix /makefile, $(SRC_DIRS))

# Include architecture and rules
include $(AMRVAC_DIR)/arch/$(ARCH).defs
include $(AMRVAC_DIR)/arch/rules.make

# Get .t files from SRC_DIRS
vpath %.t $(SRC_DIRS)

OBJECTS := $(FOBJECTS:.t=.o) 

$(LIB_AMRVAC): $(OBJECTS)
	$(RM) $@
	$(AR) rcs $@ $^

clean:
	$(RM) *.o *.mod *.F90 $(LIB_AMRVAC)
