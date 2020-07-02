        limiter= 20*'mp52
SETUP_FLAGS := -d=3
SCHEME_DIR := ../../schemes

SCHEMES := 1step_tvd 2step_tvdlf_mm 2step_tvdmu_al 3step_hll_cada rk4_tvdlf_cada ssprk43_fd_mp5 ssprk54_fd_mp5
# SCHEMES := 1step_tvd 2step_tvdlf_mm 2step_tvdmu_al 3step_hll_cada 4step_hll_mc	\
# 4step_hllc_ko rk4_tvdlf_cada ssprk43_fd_mp5 ssprk54_fd_mp5

TESTS := $(SCHEMES:%=ball_3d_%.log) $(SCHEMES:%=ball_3d_rev_%.log)

include ../../test_rules.make

# Generate dependency rules for the tests
$(foreach s, $(SCHEMES),\
	$(eval $(s:%=ball_3d_%.log): ball_3d.par $(SCHEME_DIR)/$(s).par))

$(foreach s, $(SCHEMES),\
	$(eval $(s:%=ball_3d_rev_%.log): ball_3d.par rev.par $(SCHEME_DIR)/$(s).par))
