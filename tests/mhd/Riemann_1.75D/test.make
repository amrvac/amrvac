SETUP_FLAGS := -d=1 -v=3
SCHEME_DIR := ../../schemes

SCHEMES := 2step_tvdlf_mm 2step_tvdmu_al 3step_hll_cada 3step_hlld_cada		\
4step_hll_mc 4step_hllc_ko rk4_tvdlf_cada ssprk54_fd_mp5 ssprk54_hll_mp5	\
ssprk54_hlld_mp5 3step_hll_ppm

TESTS := $(SCHEMES:%=R_1d_%.log)

include ../../test_rules.make

# Generate dependency rules for the tests
$(foreach s, $(SCHEMES),\
	$(eval $(s:%=R_1d_%.log): R_1d.par $(SCHEME_DIR)/$(s).par))
