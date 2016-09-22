# Location of setup script
SETUP_SCRIPT := $(AMRVAC_DIR)/setup.pl

# The files that define a setup
SETUP_FILES := amrvacusr.t amrvacusrpar.t amrvacsettings.t definitions.h	\
mod_indices.t makefile

SETUP_FLAGS := -d=11 -g=22 -p=rho

.PHONY: all

all: amrvac

$(SETUP_FILES):
	$(SETUP_SCRIPT) $(SETUP_FLAGS) > setup.log

# How to compile amrvac
amrvac: $(SETUP_FILES) $(AMRVAC_DIR)/src/*
	$(MAKE)

