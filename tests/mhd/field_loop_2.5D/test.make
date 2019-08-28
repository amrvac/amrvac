SETUP_FLAGS := -d=2
SCHEME_DIR := ../../schemes
SCHEMES := 3step_hll_cada 3step_hlld_cada 3step_hll_cada_ct 3step_hll_cada_ct_hll 3step_hll_cada_ct_average

TESTS := $(SCHEMES:%=fl_2d_%.log)

include ../../test_rules.make

# Generate dependency rules for the tests
$(foreach s, $(SCHEMES),\
	$(eval $(s:%=fl_2d_%.log): fl_2d.par $(SCHEME_DIR)/$(s).par))
