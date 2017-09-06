SETUP_FLAGS := -d=33
SCHEME_DIR := ../../schemes
TESTS := rc_3d_2step_tvdlf_mm.log rc_3d_2step_tvdmu_al.log		\
rc_3d_3step_hll_cada.log rc_3d_4step_hll_mc.log rc_3d_4step_hllc_ko.log	\
rc_3d_rk4_tvdlf_cada.log rc_3d_3step_hlld_cada.log

rc_3d_2step_tvdlf_mm.log: PARS = rc_3d.par $(SCHEME_DIR)/2step_tvdlf_mm.par
rc_3d_2step_tvdmu_al.log: PARS = rc_3d.par $(SCHEME_DIR)/2step_tvdmu_al.par
rc_3d_3step_hll_cada.log: PARS = rc_3d.par $(SCHEME_DIR)/3step_hll_cada.par
rc_3d_3step_hlld_cada.log: PARS = rc_3d.par $(SCHEME_DIR)/3step_hlld_cada.par
rc_3d_4step_hll_mc.log: PARS = rc_3d.par $(SCHEME_DIR)/4step_hll_mc.par
rc_3d_4step_hllc_ko.log: PARS = rc_3d.par $(SCHEME_DIR)/4step_hllc_ko.par
rc_3d_rk4_tvdlf_cada.log: PARS = rc_3d.par $(SCHEME_DIR)/rk4_tvdlf_cada.par

include ../../test_rules.make
