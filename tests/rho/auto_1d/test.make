SETUP_FLAGS := -d=1
SCHEME_DIR := ../../schemes

SCHEMES := 1step_tvd 2step_tvdlf_mm 2step_tvdmu_al trap_hll_ko 3step_hll_cada	\
3step_hll_vl 3step_hll_ppm 4step_hll_mc 4step_hllc_ko rk4_tvdlf_cada		\
ssprk43_fd_mp5 ssprk54_fd_mp5 3step_tvdlf_woodward 3step_tvdlf_superbee		\
3step_tvdlf_cada1

TESTS := $(SCHEMES:%=ball_1d_%.log) $(SCHEMES:%=ball_1d_rev_%.log)

include ../../test_rules.make

# Generate dependency rules for the tests
$(foreach s, $(SCHEMES),\
	$(eval $(s:%=ball_1d_%.log): ball_1d.par $(SCHEME_DIR)/$(s).par))

$(foreach s, $(SCHEMES),\
	$(eval $(s:%=ball_1d_rev_%.log): ball_1d.par rev.par $(SCHEME_DIR)/$(s).par))
