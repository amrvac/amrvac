SETUP_FLAGS := -d=2 -v=3
SCHEME_DIR := ../../schemes
SCHEMES := 2step_tvdlf_mm 3step_hll_cada 3step_hll_cada 4step_hll_mc \
4step_hllc_ko rk4_tvdlf_cada 3step_hlld_cada

TESTS := $(SCHEMES:%=rip_2.5d_%.log)

include ../../test_rules.make

# Generate dependency rules for the tests
$(foreach s, $(SCHEMES),\
	$(eval $(s:%=rip_2.5d_%.log): rip_2.5d.par $(SCHEME_DIR)/$(s).par))
