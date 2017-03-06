SETUP_FLAGS := -d=2
SCHEME_DIR := ../../schemes
TESTS := kh_3d_2step_tvdlf_mm.log kh_3d_2step_tvdmu_al.log		\
kh_3d_3step_hll_cada.log kh_3d_4step_hll_mc.log kh_3d_4step_hllc_ko.log	\
kh_3d_rk4_tvdlf_cada.log

kh_3d_2step_tvdlf_mm.log: PARS = kh_3d.par $(SCHEME_DIR)/2step_tvdlf_mm.par
kh_3d_2step_tvdmu_al.log: PARS = kh_3d.par $(SCHEME_DIR)/2step_tvdmu_al.par
kh_3d_3step_hll_cada.log: PARS = kh_3d.par $(SCHEME_DIR)/3step_hll_cada.par
kh_3d_4step_hll_mc.log: PARS = kh_3d.par $(SCHEME_DIR)/4step_hll_mc.par
kh_3d_4step_hllc_ko.log: PARS = kh_3d.par $(SCHEME_DIR)/4step_hllc_ko.par
kh_3d_rk4_tvdlf_cada.log: PARS = kh_3d.par $(SCHEME_DIR)/rk4_tvdlf_cada.par

include ../../test_rules.make
