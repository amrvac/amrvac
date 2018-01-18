SETUP_FLAGS := -d=2
SCHEME_DIR := ../../schemes
SCHEMES := 3step_hll_cada 4step_hll_mc 4step_hllc_ko rk4_tvdlf_cada \
3step_hlld_cada

TESTS := $(SCHEMES:%=ls_2d_%.log)

include ../../test_rules.make

# Generate dependency rules for the tests
$(foreach s, $(SCHEMES),\
	$(eval $(s:%=ls_2d_%.log): ls_2d.par $(SCHEME_DIR)/$(s).par))
