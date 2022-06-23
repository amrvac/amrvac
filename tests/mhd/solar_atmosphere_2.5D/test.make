SETUP_FLAGS := -d=2
SCHEME_DIR := ../../schemes
SCHEMES := 3step_hll_cada_ct 2step_tvdlf_mm 3step_hll_cada 3step_hlld_cada 4step_hll_mc \
4step_hllc_ko rk4_tvdlf_cada ssprk54_hlld_mp5 B0split mhd_internal_e mhd_hydrodynamic_e

TESTS := $(SCHEMES:%=solar_atm_25D_%.log)

include ../../test_rules.make

# Generate dependency rules for the tests
$(foreach s, $(SCHEMES),\
	$(eval $(s:%=solar_atm_25D_%.log): solar_atm_25D.par $(SCHEME_DIR)/$(s).par))
