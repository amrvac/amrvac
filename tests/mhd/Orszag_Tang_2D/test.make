SETUP_FLAGS := -d=2
SCHEME_DIR := ../../schemes
TESTS := ot_2d_2step_tvdlf_mm.log ot_2d_2step_tvdmu_al.log		\
ot_2d_3step_hll_cada.log ot_2d_4step_hll_mc.log ot_2d_4step_hllc_ko.log	\
ot_2d_rk4_tvdlf_cada.log ot_2d_3step_hlld_cada.log

ot_2d_2step_tvdlf_mm.log: PARS = ot_2d.par $(SCHEME_DIR)/2step_tvdlf_mm.par
ot_2d_2step_tvdmu_al.log: PARS = ot_2d.par $(SCHEME_DIR)/2step_tvdmu_al.par
ot_2d_3step_hll_cada.log: PARS = ot_2d.par $(SCHEME_DIR)/3step_hll_cada.par
ot_2d_4step_hll_mc.log: PARS = ot_2d.par $(SCHEME_DIR)/4step_hll_mc.par
ot_2d_4step_hllc_ko.log: PARS = ot_2d.par $(SCHEME_DIR)/4step_hllc_ko.par
ot_2d_rk4_tvdlf_cada.log: PARS = ot_2d.par $(SCHEME_DIR)/rk4_tvdlf_cada.par
ot_2d_3step_hlld_cada.log: PARS = ot_2d.par $(SCHEME_DIR)/3step_hlld_cada.par

include ../../test_rules.make
