SETUP_FLAGS := -d=2
SCHEME_DIR := ../../schemes
TESTS := gem_2d_2step_tvdlf_mm.log gem_2d_2step_tvdmu_al.log		\
gem_2d_3step_hll_cada.log gem_2d_4step_hll_mc.log gem_2d_4step_hllc_ko.log	\
gem_2d_rk4_tvdlf_cada.log

gem_2d_2step_tvdlf_mm.log: PARS = gem_2d.par $(SCHEME_DIR)/2step_tvdlf_mm.par
gem_2d_2step_tvdmu_al.log: PARS = gem_2d.par $(SCHEME_DIR)/2step_tvdmu_al.par
gem_2d_3step_hll_cada.log: PARS = gem_2d.par $(SCHEME_DIR)/3step_hll_cada.par
gem_2d_4step_hll_mc.log: PARS = gem_2d.par $(SCHEME_DIR)/4step_hll_mc.par
gem_2d_4step_hllc_ko.log: PARS = gem_2d.par $(SCHEME_DIR)/4step_hllc_ko.par
gem_2d_rk4_tvdlf_cada.log: PARS = gem_2d.par $(SCHEME_DIR)/rk4_tvdlf_cada.par

include ../../test_rules.make
