SETUP_FLAGS := -d=2
SCHEME_DIR := ../../schemes
TESTS := ti_2d_2step_tvdlf_mm.log ti_2d_2step_tvdmu_al.log		\
ti_2d_3step_hll_cada.log ti_2d_4step_hll_mc.log ti_2d_4step_hllc_ko.log	\
ti_2d_rk4_tvdlf_cada.log ti_2d_ssprk54_fd_mp5.log ti_2d_ssprk54_hll_mp5.log

ti_2d_2step_tvdlf_mm.log: PARS = ti_2d.par $(SCHEME_DIR)/2step_tvdlf_mm.par
ti_2d_2step_tvdmu_al.log: PARS = ti_2d.par $(SCHEME_DIR)/2step_tvdmu_al.par
ti_2d_3step_hll_cada.log: PARS = ti_2d.par $(SCHEME_DIR)/3step_hll_cada.par
ti_2d_4step_hll_mc.log: PARS = ti_2d.par $(SCHEME_DIR)/4step_hll_mc.par
ti_2d_4step_hllc_ko.log: PARS = ti_2d.par $(SCHEME_DIR)/4step_hllc_ko.par
ti_2d_rk4_tvdlf_cada.log: PARS = ti_2d.par $(SCHEME_DIR)/rk4_tvdlf_cada.par
ti_2d_ssprk54_fd_mp5.log: PARS = ti_2d.par $(SCHEME_DIR)/ssprk54_fd_mp5.par
ti_2d_ssprk54_hll_mp5.log: PARS = ti_2d.par $(SCHEME_DIR)/ssprk54_hll_mp5.par

include ../../test_rules.make
