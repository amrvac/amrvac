SETUP_FLAGS := -d=1
SCHEME_DIR := ../../schemes
TESTS := rm_1d_2step_tvdlf_mm.log rm_1d_2step_tvdmu_al.log		\
rm_1d_3step_hll_cada.log rm_1d_4step_hll_mc.log rm_1d_4step_hllc_ko.log	\
rm_1d_rk4_tvdlf_cada.log

rm_1d_2step_tvdlf_mm.log: PARS = rm_1d.par $(SCHEME_DIR)/2step_tvdlf_mm.par
rm_1d_2step_tvdmu_al.log: PARS = rm_1d.par $(SCHEME_DIR)/2step_tvdmu_al.par
rm_1d_3step_hll_cada.log: PARS = rm_1d.par $(SCHEME_DIR)/3step_hll_cada.par
rm_1d_4step_hll_mc.log: PARS = rm_1d.par $(SCHEME_DIR)/4step_hll_mc.par
rm_1d_4step_hllc_ko.log: PARS = rm_1d.par $(SCHEME_DIR)/4step_hllc_ko.par
rm_1d_rk4_tvdlf_cada.log: PARS = rm_1d.par $(SCHEME_DIR)/rk4_tvdlf_cada.par

include ../../test_rules.make
