SETUP_FLAGS := -d=2
SCHEME_DIR := ../../schemes
SCHEMES := 3step_hll_cada_ct 3step_hll_cada_ct_hll 3step_hll_cada_ct_average 3step_hlld_cada_ct 3step_hll_cada_ambi 3step_hll_cada_ct_ambi 3step_hll_cada_ct_ambi_ex 3step_hll_cada_ambi_ex

TESTS := $(SCHEMES:%=ot_2d_%.log)

include ../../test_rules.make

# Generate dependency rules for the tests
$(foreach s, $(SCHEMES),\
	$(eval $(s:%=ot_2d_%.log): ot_2d.par $(SCHEME_DIR)/$(s).par))
