SETUP_FLAGS := -d=2
SCHEME_DIR := ../../schemes

SCHEMES := 1step_tvd 2step_tvdlf_mm 2step_tvdmu_al trap_hll_ko 3step_hll_cada	\
3step_hll_ppm rk4_tvdlf_cada ssprk43_fd_mp5		\
ssprk54_fd_mp5
# SCHEMES := 1step_tvd 2step_tvdlf_mm 2step_tvdmu_al trap_hll_ko 3step_hll_cada	\
# 3step_hll_ppm 4step_hll_mc 4step_hllc_ko rk4_tvdlf_cada ssprk43_fd_mp5		\
# ssprk54_fd_mp5

TESTS := $(SCHEMES:%=ball_2d_%.log) $(SCHEMES:%=ball_2d_rev_%.log) $(SCHEMES:%=vac_2d_%.log)

include ../../test_rules.make

# Generate dependency rules for the tests
$(foreach s, $(SCHEMES),\
	$(eval $(s:%=ball_2d_%.log): ball_2d.par $(SCHEME_DIR)/$(s).par))

$(foreach s, $(SCHEMES),\
	$(eval $(s:%=ball_2d_rev_%.log): ball_2d.par rev.par $(SCHEME_DIR)/$(s).par))

$(foreach s, $(SCHEMES),\
	$(eval $(s:%=vac_2d_%.log): vac_2d.par $(SCHEME_DIR)/$(s).par))
