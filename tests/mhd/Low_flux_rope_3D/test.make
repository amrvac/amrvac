SETUP_FLAGS := -d=3
SCHEME_DIR := ../../schemes
SCHEMES := 2step_tvdlf_mm 2step_tvdmu_al 3step_hll_cada 4step_hll_mc rk4_tvdlf_cada

TESTS := $(SCHEMES:%=lfr_3d_%.log)

include ../../test_rules.make

# Generate dependency rules for the tests
$(foreach s, $(SCHEMES),\
	$(eval $(s:%=lfr_3d_%.log): lfr_3d.par $(SCHEME_DIR)/$(s).par))
