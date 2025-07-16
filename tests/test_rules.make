_red := "\\e[31m"
_green := "\\e[32m"
_reset := "\\e[m"



# This is a template Makefile to simplify writing tests
ifndef AMRVAC_DIR
$(error AMRVAC_DIR is not set)
endif

ARCH ?= gnu

# Disable built-in make rules
.SUFFIXES:

# Tool to compare test output numerically
LOG_CMP := $(AMRVAC_DIR)/tools/fortran/compare_logs

# Number of MPI processes to use
NUM_PROCS ?= 4

# force is a dummy to force re-running tests
.PHONY: all clean force

all: $(TESTS)

clean:
	$(RM) $(TESTS) amrvac *.vtu *.dat *.log *.f90 *.mod

# Include architecture and compilation rules for the compare_log utility
include $(AMRVAC_DIR)/arch/$(ARCH).mk

F90 := $(compile)
F90FLAGS := $(f90_flags)

%.log: $(LOG_CMP) amrvac force
	@$(RM) $@		# Remove log to prevent pass when aborted
# for Intel same machine
# @mpirun -genv I_MPI_FABRICS shm  -np $(NUM_PROCS) ./amrvac -i $(filter %.par,$^) > run.log
	@mpirun   -np $(NUM_PROCS) ./amrvac -i $(filter %.par,$^) > run.log
	@if $(LOG_CMP) 1.0e-5 1.0e-8 $@ correct_output/$@ ; \
	then echo -e "$(_green)PASSED$(_reset) $@" ; \
	else echo -e "$(_red)** FAILED **$(_reset) $@" ; \
	fi

amrvac: force		# Always try to build
	@$(MAKE) arch=$(ARCH)

# To make sure the comparison utility can be build
$(LOG_CMP).o: $(LOG_CMP).f
	$(F90) $(F90FLAGS) -c $< -o $@

$(LOG_CMP): $(LOG_CMP).o
	$(F90) $< -o $@
