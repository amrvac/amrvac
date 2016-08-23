# This is a template Makefile to simplify writing tests

# Location of setup script
SETUP_SCRIPT := $(AMRVAC_DIR)/setup.pl

# Number of MPI processes to use
NUM_PROCS ?= 4

# The default targets are the $(TESTS)
all: $(TESTS)

# Declare amrvac phony so that it is always recompiled
.PHONY: all amrvac $(TESTS)

# How to compile amrvac
amrvac:
	@$(SETUP_SCRIPT) $(SETUP_FLAGS) &> setup.log || "FAILED setup, see setup.log"
	@$(MAKE) &> compilation.log || "FAILED compilation, see compilation.log"

# How to run the tests
$(TESTS): amrvac
	@$(RM) amrvac.log	# Clear old log file
	@mpirun -np $(NUM_PROCS) ./amrvac -i $(@:%=%.par) &> run.log || "FAILED run, see run.log"
	@if diff ${@:%=%.par.log} correct_output/${@:%=%.par.log} ; \
	then echo "PASSED $@" ; \
	else echo "FAILED $@" ; \
	fi
