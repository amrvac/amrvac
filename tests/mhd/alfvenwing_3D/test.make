# Name of .par file you want to use
PAR_FILE := AW_3D.par

# Base file name as set in the .par file
BASE_NAME := my_test

# Change to -d=2 or -d=3 for 2D/3D problems
SETUP_FLAGS := -d=3

# This directory contains combinations of numerical schemes that you can use
SCHEME_DIR := ../../schemes

# Which of the schemes you want to use for your test
SCHEMES := 2step_tvdlf_mm 3step_hll_cada 3step_hll_vl 3step_hll_ppm

# Define the targets for the makefile, which are the log files produced by your
# program.
TESTS := $(SCHEMES:%=$(BASE_NAME)_%.log)

# Include rules on how to compile and run the tests
include ../../test_rules.make

# Generate dependency rules for the tests, which are used to run them
$(foreach s, $(SCHEMES),\
    $(eval $(s:%=$(BASE_NAME)_%.log): $(PAR_FILE) $(SCHEME_DIR)/$(s).par))
