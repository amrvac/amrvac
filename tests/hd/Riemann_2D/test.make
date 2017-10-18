SETUP_FLAGS := -d=22
SCHEME_DIR := ../../schemes
TESTS := rm_2d_2step_tvdlf_mm.log rm_2d_2step_tvdmu_al.log		\
rm_2d_3step_hll_cada.log rm_2d_4step_hll_mc.log rm_2d_4step_hllc_ko.log	\
rm_2d_rk4_tvdlf_cada.log

rm_2d_2step_tvdlf_mm.log: PARS = rm_2d.par $(SCHEME_DIR)/2step_tvdlf_mm.par
rm_2d_2step_tvdmu_al.log: PARS = rm_2d.par $(SCHEME_DIR)/2step_tvdmu_al.par
rm_2d_3step_hll_cada.log: PARS = rm_2d.par $(SCHEME_DIR)/3step_hll_cada.par
rm_2d_4step_hll_mc.log: PARS = rm_2d.par $(SCHEME_DIR)/4step_hll_mc.par
rm_2d_4step_hllc_ko.log: PARS = rm_2d.par $(SCHEME_DIR)/4step_hllc_ko.par
rm_2d_rk4_tvdlf_cada.log: PARS = rm_2d.par $(SCHEME_DIR)/rk4_tvdlf_cada.par

include ../../test_rules.make
