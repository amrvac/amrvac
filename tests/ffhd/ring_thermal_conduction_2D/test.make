PAR_FILE := tc_ring_ffhd.par
BASE_NAME := tc_ring_ffhd
SETUP_FLAGS := -d=2
SCHEME_DIR := ../../schemes
SCHEMES := 2step_tvdlf_mm 3step_hll_cada 3step_hll_vl
TESTS := $(SCHEMES:%=$(BASE_NAME)_%.log)
include ../../test_rules.make
$(foreach s, $(SCHEMES),\
    $(eval $(s:%=$(BASE_NAME)_%.log): $(PAR_FILE) $(SCHEME_DIR)/$(s).par))
