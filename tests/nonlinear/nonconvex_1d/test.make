SETUP_FLAGS := -d=1 -v=1
SCHEME_DIR := ../../schemes

SCHEMES := 3step_hll_cada 3step_hll_vl 3step_hll_ppm 4step_hll_mc \
3step_tvdlf_woodward 3step_tvdlf_superbee 3step_tvdlf_cada1

TESTS := $(SCHEMES:%=compound_1d_%.log)

include ../../test_rules.make

# Generate dependency rules for the tests
$(foreach s, $(SCHEMES),\
	$(eval $(s:%=compound_1d_%.log): compound_1d.par $(SCHEME_DIR)/$(s).par))
