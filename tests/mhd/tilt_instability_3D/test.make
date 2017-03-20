SETUP_FLAGS := -d=3
SCHEME_DIR := ../../schemes
TESTS := ti_3d_2step_tvdlf_mm.log ti_3d_2step_tvdmu_al.log		\
ti_3d_3step_hll_cada.log ti_3d_4step_hll_mc.log ti_3d_4step_hllc_ko.log	\
ti_3d_rk4_tvdlf_cada.log

ti_3d_2step_tvdlf_mm.log: PARS = ti_3d.par $(SCHEME_DIR)/2step_tvdlf_mm.par
ti_3d_2step_tvdmu_al.log: PARS = ti_3d.par $(SCHEME_DIR)/2step_tvdmu_al.par
ti_3d_3step_hll_cada.log: PARS = ti_3d.par $(SCHEME_DIR)/3step_hll_cada.par
ti_3d_4step_hll_mc.log: PARS = ti_3d.par $(SCHEME_DIR)/4step_hll_mc.par
ti_3d_4step_hllc_ko.log: PARS = ti_3d.par $(SCHEME_DIR)/4step_hllc_ko.par
ti_3d_rk4_tvdlf_cada.log: PARS = ti_3d.par $(SCHEME_DIR)/rk4_tvdlf_cada.par

include ../../test_rules.make
