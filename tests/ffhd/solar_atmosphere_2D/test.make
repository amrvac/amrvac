PAR_FILE := solar_atm_ffhdver.par
BASE_NAME := solar_atm_ffhdver
SETUP_FLAGS := -d=2
SCHEME_DIR := ../../schemes
SCHEMES := 2step_tvdlf_mm 3step_hll_cada 3step_hll_vl
TESTS := $(SCHEMES:%=$(BASE_NAME)_%.log)
include ../../test_rules.make
$(foreach s, $(SCHEMES),\
    $(eval $(s:%=$(BASE_NAME)_%.log): $(PAR_FILE) $(SCHEME_DIR)/$(s).par))
