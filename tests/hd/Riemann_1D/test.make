SETUP_FLAGS := -d=1
SCHEME_DIR := ../../schemes

SCHEMES := 2step_tvdlf_mm 2step_tvdmu_al 3step_hll_cada 4step_hll_mc	\
4step_hllc_ko rk4_tvdlf_cada 3step_hll_ppm

TESTS := $(SCHEMES:%=rm_1d_%.log)

include ../../test_rules.make

# Generate dependency rules for the tests
$(foreach s, $(SCHEMES),\
	$(eval $(s:%=rm_1d_%.log): rm_1d.par $(SCHEME_DIR)/$(s).par))

