# This is a template Makefile to simplify writing tests

# Disable built-in rules
.SUFFIXES:

# Location of setup script
SETUP_SCRIPT := $(AMRVAC_DIR)/setup.pl

# Tool to compare test output numerically
LOG_CMP := $(AMRVAC_DIR)/tools/python/compare_logs.py

# Number of MPI processes to use
NUM_PROCS ?= 4

# The files that define a setup
SETUP_FILES := makefile

.PHONY: all always_rebuild

all: $(TESTS)

%.log: amrvac always_rebuild
	@$(RM) $@		# Remove log to prevent pass when aborted
	mpirun -np $(NUM_PROCS) ./amrvac -i $(PARS) > run.log
	@if $(LOG_CMP) $@ correct_output/$@ ; \
	then echo "PASSED $@" ; \
	else echo "** FAILED ** $@" ; \
	fi

$(SETUP_FILES):
	$(SETUP_SCRIPT) $(SETUP_FLAGS) > setup.log

# How to compile amrvac
amrvac: $(SETUP_FILES) $(AMRVAC_DIR)/src/*
	$(MAKE) > compilation.log


