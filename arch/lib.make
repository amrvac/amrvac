ifndef GMUNU_DIR
$(error GMUNU_DIR is not set)
endif

ifndef ARCH
$(error build.make: ARCH is not set)
endif

ifndef NDIM
$(error build.make: NDIM is not set)
endif

SRC_DIRS := . modules amr io physics limiter hd

SRC_DIRS := $(addprefix $(GMUNU_DIR)/src/, $(SRC_DIRS))
LIB_GMUNU := libgmunu.a
PPFLAGS :=

.PHONY: libgmunu clean
.PRECIOUS: %.f			# Don't remove intermediate .f files

libgmunu: $(LIB_GMUNU) gmunu.f

# Include makefiles, which define FOBJECTS and dependencies
include $(addsuffix /makefile, $(SRC_DIRS))

# Include architecture and rules
include $(GMUNU_DIR)/arch/$(ARCH).defs
include $(GMUNU_DIR)/arch/rules.make

# Get .t files from SRC_DIRS
vpath %.t $(SRC_DIRS)

OBJECTS := $(FOBJECTS:.t=.o) $(INCLUDES:.t=.o)

$(LIB_GMUNU): $(OBJECTS)
	$(RM) $@
	$(AR) rcs $@ $^

clean:
	$(RM) *.o *.mod *.f $(LIB_GMUNU)

# INCLUDES are always compiled before FOBJECTS
$(FOBJECTS:.t=.o): $(INCLUDES:.t=.mod)
