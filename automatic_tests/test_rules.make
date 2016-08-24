# This is a template Makefile to simplify writing tests

# Location of setup script
SETUP_SCRIPT := $(AMRVAC_DIR)/setup.pl

# Tool to compare test output numerically
LOG_CMP := $(AMRVAC_DIR)/tools/python/compare_logs.py

# Number of MPI processes to use
NUM_PROCS ?= 4

# Colors to highlight the test results, see
# https://en.wikipedia.org/wiki/ANSI_escape_code
RED := \033[0;31m
GREEN := \033[0;32m
NOCOLOR := \033[0m

# The files that define a setup
SETUP_FILES := amrvacusr.t amrvacusrpar.t amrvacsettings.t definitions.h	\
mod_indices.t makefile

# The default targets are the $(TESTS)
all: $(TESTS)

# Declare amrvac phony so that it is always recompiled
.PHONY: all $(TESTS)

$(SETUP_FILES):
	$(SETUP_SCRIPT) $(SETUP_FLAGS) > setup.log

# How to compile amrvac
amrvac: $(SETUP_FILES)
	$(MAKE) &> compilation.log

# How to run the tests
$(TESTS): amrvac
	mpirun -np $(NUM_PROCS) ./amrvac -i $(@:%=%.par) > run.log
	@if $(LOG_CMP) ${@:%=%.par.log} correct_output/${@:%=%.par.log} ; \
	then echo -e "$(GREEN)PASSED $@ $(NOCOLOR)" ; \
	else echo -e "$(RED)FAILED $@ $(NOCOLOR)" ; \
	fi
