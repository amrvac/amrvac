SETUP_FLAGS := -d=1
SCHEME_DIR := ../../schemes

SCHEMES := 2step_tvdlf_vl 2step_tvdmu_al 3step_hll_vl 3step_hll_ko

TESTS := $(SCHEMES:%=cak1d_%.log)

include ../../test_rules.make

# Generate dependency rules for the tests, which are used to run them
$(foreach s, $(SCHEMES),\
    $(eval $(s:%=cak1d_%.log): cak_1d.par $(SCHEME_DIR)/$(s).par))