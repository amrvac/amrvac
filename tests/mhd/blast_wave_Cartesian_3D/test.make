SETUP_FLAGS := -d=3
SCHEME_DIR := ../../schemes
SCHEMES := 3step_hll_cada_ct 2step_tvdlf_mm 3step_hll_cada 3step_hll_ko 4step_hll_mc \
4step_hllc_ko rk4_tvdlf_cada 3step_hlld_cada ssprk54_hlld_mp5

TESTS := $(SCHEMES:%=bw_3d_%.log)

include ../../test_rules.make

# Generate dependency rules for the tests
$(foreach s, $(SCHEMES),\
	$(eval $(s:%=bw_3d_%.log): bw_3d.par $(SCHEME_DIR)/$(s).par))
