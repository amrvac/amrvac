SETUP_FLAGS := -d=2
SCHEME_DIR := ../../schemes
SCHEMES := 3step_hll_cada 3step_hll_cada_ct mhd_internal_e

TESTS := $(SCHEMES:%=tc_2d_%.log)

include ../../test_rules.make

# Generate dependency rules for the tests
$(foreach s, $(SCHEMES),\
	$(eval $(s:%=tc_2d_%.log): tc_2d.par $(SCHEME_DIR)/$(s).par))
