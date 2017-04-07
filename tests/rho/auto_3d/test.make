SETUP_FLAGS := -d=3
SCHEME_DIR := ../../schemes

TESTS := ball_3d_1step_tvd.log ball_3d_2step_tvdlf_mm.log			\
ball_3d_2step_tvdmu_al.log ball_3d_3step_hll_cada.log ball_3d_4step_hll_mc.log	\
ball_3d_4step_hllc_ko.log ball_3d_rk4_tvdlf_cada.log ball_3d_ssprk43_fd_mp5.log	\
ball_3d_ssprk54_fd_mp5.log ball_3d_rev_1step_tvd.log				\
ball_3d_rev_2step_tvdlf_mm.log ball_3d_rev_2step_tvdmu_al.log			\
ball_3d_rev_3step_hll_cada.log ball_3d_rev_4step_hll_mc.log			\
ball_3d_rev_4step_hllc_ko.log ball_3d_rev_rk4_tvdlf_cada.log			\
ball_3d_rev_ssprk43_fd_mp5.log ball_3d_rev_ssprk54_fd_mp5.log

ball_3d_1step_tvd.log: PARS = ball_3d.par $(SCHEME_DIR)/1step_tvd.par
ball_3d_2step_tvdlf_mm.log: PARS = ball_3d.par $(SCHEME_DIR)/2step_tvdlf_mm.par
ball_3d_2step_tvdmu_al.log: PARS = ball_3d.par $(SCHEME_DIR)/2step_tvdmu_al.par
ball_3d_3step_hll_cada.log: PARS = ball_3d.par $(SCHEME_DIR)/3step_hll_cada.par
ball_3d_4step_hll_mc.log: PARS = ball_3d.par $(SCHEME_DIR)/4step_hll_mc.par
ball_3d_4step_hllc_ko.log: PARS = ball_3d.par $(SCHEME_DIR)/4step_hllc_ko.par
ball_3d_rk4_tvdlf_cada.log: PARS = ball_3d.par $(SCHEME_DIR)/rk4_tvdlf_cada.par
ball_3d_ssprk43_fd_mp5.log: PARS = ball_3d.par $(SCHEME_DIR)/ssprk43_fd_mp5.par
ball_3d_ssprk54_fd_mp5.log: PARS = ball_3d.par $(SCHEME_DIR)/ssprk54_fd_mp5.par

ball_3d_rev_1step_tvd.log: PARS = ball_3d.par rev.par $(SCHEME_DIR)/1step_tvd.par
ball_3d_rev_2step_tvdlf_mm.log: PARS = ball_3d.par rev.par $(SCHEME_DIR)/2step_tvdlf_mm.par
ball_3d_rev_2step_tvdmu_al.log: PARS = ball_3d.par rev.par $(SCHEME_DIR)/2step_tvdmu_al.par
ball_3d_rev_3step_hll_cada.log: PARS = ball_3d.par rev.par $(SCHEME_DIR)/3step_hll_cada.par
ball_3d_rev_4step_hll_mc.log: PARS = ball_3d.par rev.par $(SCHEME_DIR)/4step_hll_mc.par
ball_3d_rev_4step_hllc_ko.log: PARS = ball_3d.par rev.par $(SCHEME_DIR)/4step_hllc_ko.par
ball_3d_rev_rk4_tvdlf_cada.log: PARS = ball_3d.par rev.par $(SCHEME_DIR)/rk4_tvdlf_cada.par
ball_3d_rev_ssprk43_fd_mp5.log: PARS = ball_3d.par rev.par $(SCHEME_DIR)/ssprk43_fd_mp5.par
ball_3d_rev_ssprk54_fd_mp5.log: PARS = ball_3d.par rev.par $(SCHEME_DIR)/ssprk54_fd_mp5.par

include ../../test_rules.make
