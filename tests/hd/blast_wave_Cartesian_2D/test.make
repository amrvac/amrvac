SETUP_FLAGS := -d=2
SCHEME_DIR := ../../schemes
SCHEMES := 2step_tvdlf_mm 2step_tvdmu_al 3step_hll_cada 4step_hll_mc	\
4step_hllc_ko rk4_tvdlf_cada

TESTS := $(SCHEMES:%=bw_2d_%.log)

include ../../test_rules.make

# Generate dependency rules for the tests
$(foreach s, $(SCHEMES),\
	$(eval $(s:%=bw_2d_%.log): bw_2d.par $(SCHEME_DIR)/$(s).par))
