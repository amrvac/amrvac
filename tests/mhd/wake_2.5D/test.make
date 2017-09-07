SETUP_FLAGS := -d=2
SCHEME_DIR := ../../schemes
TESTS := wake_2d_2step_tvdlf_mm.log wake_2d_2step_tvdmu_al.log		\
wake_2d_3step_hll_cada.log wake_2d_4step_hll_mc.log wake_2d_4step_hllc_ko.log	\
wake_2d_rk4_tvdlf_cada.log wake_2d_3step_hlld_cada.log

wake_2d_2step_tvdlf_mm.log: PARS = wake_2d.par $(SCHEME_DIR)/2step_tvdlf_mm.par
wake_2d_2step_tvdmu_al.log: PARS = wake_2d.par $(SCHEME_DIR)/2step_tvdmu_al.par
wake_2d_3step_hll_cada.log: PARS = wake_2d.par $(SCHEME_DIR)/3step_hll_cada.par
wake_2d_3step_hlld_cada.log: PARS = wake_2d.par $(SCHEME_DIR)/3step_hlld_cada.par
wake_2d_4step_hll_mc.log: PARS = wake_2d.par $(SCHEME_DIR)/4step_hll_mc.par
wake_2d_4step_hllc_ko.log: PARS = wake_2d.par $(SCHEME_DIR)/4step_hllc_ko.par
wake_2d_rk4_tvdlf_cada.log: PARS = wake_2d.par $(SCHEME_DIR)/rk4_tvdlf_cada.par

include ../../test_rules.make
