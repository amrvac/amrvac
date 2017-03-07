SETUP_FLAGS := -d=2
SCHEME_DIR := ../../schemes
TESTS := ls_2d_2step_tvdlf_mm.log ls_2d_2step_tvdmu_al.log		\
ls_2d_3step_hll_cada.log ls_2d_4step_hll_mc.log ls_2d_4step_hllc_ko.log	\
ls_2d_rk4_tvdlf_cada.log

ls_2d_2step_tvdlf_mm.log: PARS = ls_2d.par $(SCHEME_DIR)/2step_tvdlf_mm.par
ls_2d_2step_tvdmu_al.log: PARS = ls_2d.par $(SCHEME_DIR)/2step_tvdmu_al.par
ls_2d_3step_hll_cada.log: PARS = ls_2d.par $(SCHEME_DIR)/3step_hll_cada.par
ls_2d_4step_hll_mc.log: PARS = ls_2d.par $(SCHEME_DIR)/4step_hll_mc.par
ls_2d_4step_hllc_ko.log: PARS = ls_2d.par $(SCHEME_DIR)/4step_hllc_ko.par
ls_2d_rk4_tvdlf_cada.log: PARS = ls_2d.par $(SCHEME_DIR)/rk4_tvdlf_cada.par

include ../../test_rules.make
