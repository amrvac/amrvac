SETUP_FLAGS := -d=3
SCHEME_DIR := ../../schemes
TESTS := bw_cylindrical_3D_2step_tvdlf_mm.log \
bw_cylindrical_3D_3step_hll_cada.log bw_cylindrical_3D_4step_hll_mc.log bw_cylindrical_3D_4step_hllc_ko.log	\
bw_cylindrical_3D_rk4_tvdlf_cada.log

bw_cylindrical_3D_2step_tvdlf_mm.log: PARS = bw_cylindrical_3D.par $(SCHEME_DIR)/2step_tvdlf_mm.par
bw_cylindrical_3D_3step_hll_cada.log: PARS = bw_cylindrical_3D.par $(SCHEME_DIR)/3step_hll_cada.par
bw_cylindrical_3D_4step_hll_mc.log: PARS = bw_cylindrical_3D.par $(SCHEME_DIR)/4step_hll_mc.par
bw_cylindrical_3D_4step_hllc_ko.log: PARS = bw_cylindrical_3D.par $(SCHEME_DIR)/4step_hllc_ko.par
bw_cylindrical_3D_rk4_tvdlf_cada.log: PARS = bw_cylindrical_3D.par $(SCHEME_DIR)/rk4_tvdlf_cada.par

include ../../test_rules.make
