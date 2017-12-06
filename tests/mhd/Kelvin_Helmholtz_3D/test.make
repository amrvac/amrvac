SETUP_FLAGS := -d=3
SCHEME_DIR := ../../schemes
TESTS := kh_3D_2step_tvdlf_mm.log kh_3D_2step_tvdmu_al.log		\
kh_3D_3step_hll_cada.log kh_3D_4step_hll_mc.log kh_3D_4step_hllc_ko.log	\
kh_3D_rk4_tvdlf_cada.log kh_3D_3step_hlld_cada.log

kh_3D_2step_tvdlf_mm.log: kh_3d.par $(SCHEME_DIR)/2step_tvdlf_mm.par
kh_3D_2step_tvdmu_al.log: kh_3d.par $(SCHEME_DIR)/2step_tvdmu_al.par
kh_3D_3step_hll_cada.log: kh_3d.par $(SCHEME_DIR)/3step_hll_cada.par
kh_3D_4step_hll_mc.log: kh_3d.par $(SCHEME_DIR)/4step_hll_mc.par
kh_3D_4step_hllc_ko.log: kh_3d.par $(SCHEME_DIR)/4step_hllc_ko.par
kh_3D_rk4_tvdlf_cada.log: kh_3d.par $(SCHEME_DIR)/rk4_tvdlf_cada.par
kh_3D_3step_hlld_cada.log: kh_3d.par $(SCHEME_DIR)/3step_hlld_cada.par

include ../../test_rules.make
