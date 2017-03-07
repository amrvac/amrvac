SETUP_FLAGS := -d=2
SCHEME_DIR := ../../schemes
TESTS := sc_2d_2step_tvdlf_mm.log sc_2d_2step_tvdmu_al.log		\
sc_2d_3step_hll_cada.log sc_2d_4step_hll_mc.log sc_2d_4step_hllc_ko.log	\
sc_2d_rk4_tvdlf_cada.log

sc_2d_2step_tvdlf_mm.log: PARS = sc_2d.par $(SCHEME_DIR)/2step_tvdlf_mm.par
sc_2d_2step_tvdmu_al.log: PARS = sc_2d.par $(SCHEME_DIR)/2step_tvdmu_al.par
sc_2d_3step_hll_cada.log: PARS = sc_2d.par $(SCHEME_DIR)/3step_hll_cada.par
sc_2d_4step_hll_mc.log: PARS = sc_2d.par $(SCHEME_DIR)/4step_hll_mc.par
sc_2d_4step_hllc_ko.log: PARS = sc_2d.par $(SCHEME_DIR)/4step_hllc_ko.par
sc_2d_rk4_tvdlf_cada.log: PARS = sc_2d.par $(SCHEME_DIR)/rk4_tvdlf_cada.par

include ../../test_rules.make
