SETUP_FLAGS := -d=33
SCHEME_DIR := ../../schemes
TESTS := rt_3d_2step_tvdlf_mm.log rt_3d_2step_tvdmu_al.log		\
rt_3d_3step_hll_cada.log rt_3d_4step_hll_mc.log rt_3d_4step_hllc_ko.log	\
rt_3d_rk4_tvdlf_cada.log

rt_3d_2step_tvdlf_mm.log: PARS = test.par $(SCHEME_DIR)/2step_tvdlf_mm.par
rt_3d_2step_tvdmu_al.log: PARS = test.par $(SCHEME_DIR)/2step_tvdmu_al.par
rt_3d_3step_hll_cada.log: PARS = test.par $(SCHEME_DIR)/3step_hll_cada.par
rt_3d_4step_hll_mc.log: PARS = test.par $(SCHEME_DIR)/4step_hll_mc.par
rt_3d_4step_hllc_ko.log: PARS = test.par $(SCHEME_DIR)/4step_hllc_ko.par
rt_3d_rk4_tvdlf_cada.log: PARS = test.par $(SCHEME_DIR)/rk4_tvdlf_cada.par

include ../../test_rules.make
