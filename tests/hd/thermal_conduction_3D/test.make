SETUP_FLAGS := -d=33
SCHEME_DIR := ../../schemes
TESTS := tc_3d_2step_tvdlf_mm.log tc_3d_2step_tvdmu_al.log		\
tc_3d_3step_hll_cada.log tc_3d_4step_hll_mc.log tc_3d_4step_hllc_ko.log	\
tc_3d_rk4_tvdlf_cada.log

tc_3d_2step_tvdlf_mm.log: tc_3d.par $(SCHEME_DIR)/2step_tvdlf_mm.par
tc_3d_2step_tvdmu_al.log: tc_3d.par $(SCHEME_DIR)/2step_tvdmu_al.par
tc_3d_3step_hll_cada.log: tc_3d.par $(SCHEME_DIR)/3step_hll_cada.par
tc_3d_4step_hll_mc.log: tc_3d.par $(SCHEME_DIR)/4step_hll_mc.par
tc_3d_4step_hllc_ko.log: tc_3d.par $(SCHEME_DIR)/4step_hllc_ko.par
tc_3d_rk4_tvdlf_cada.log: tc_3d.par $(SCHEME_DIR)/rk4_tvdlf_cada.par

include ../../test_rules.make
