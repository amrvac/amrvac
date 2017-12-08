SETUP_FLAGS := -d=3
SCHEME_DIR := ../../schemes
TESTS := bw_3d_2step_tvdlf_mm.log bw_3d_3step_hll_cada.log \
bw_3d_3step_hllc_cada.log bw_3d_4step_hll_mc.log bw_3d_4step_hllc_ko.log	\
bw_3d_rk4_tvdlf_cada.log

bw_3d_2step_tvdlf_mm.log: bw_3d.par $(SCHEME_DIR)/2step_tvdlf_mm.par
bw_3d_3step_hll_cada.log: bw_3d.par $(SCHEME_DIR)/3step_hll_cada.par
bw_3d_3step_hllc_cada.log: bw_3d.par $(SCHEME_DIR)/3step_hllc_cada.par
bw_3d_4step_hll_mc.log: bw_3d.par $(SCHEME_DIR)/4step_hll_mc.par
bw_3d_4step_hllc_ko.log: bw_3d.par $(SCHEME_DIR)/4step_hllc_ko.par
bw_3d_rk4_tvdlf_cada.log: bw_3d.par $(SCHEME_DIR)/rk4_tvdlf_cada.par

include ../../test_rules.make
