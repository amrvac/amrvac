SETUP_FLAGS := -d=2
SCHEME_DIR := ../../schemes
TESTS := lfff_2d_2step_tvdlf_mm.log lfff_2d_2step_tvdmu_al.log		\
lfff_2d_3step_hll_cada.log lfff_2d_4step_hll_mc.log lfff_2d_4step_hllc_ko.log	\
lfff_2d_rk4_tvdlf_cada.log

lfff_2d_2step_tvdlf_mm.log: lfff_2d.par $(SCHEME_DIR)/2step_tvdlf_mm.par
lfff_2d_2step_tvdmu_al.log: lfff_2d.par $(SCHEME_DIR)/2step_tvdmu_al.par
lfff_2d_3step_hll_cada.log: lfff_2d.par $(SCHEME_DIR)/3step_hll_cada.par
lfff_2d_4step_hll_mc.log: lfff_2d.par $(SCHEME_DIR)/4step_hll_mc.par
lfff_2d_4step_hllc_ko.log: lfff_2d.par $(SCHEME_DIR)/4step_hllc_ko.par
lfff_2d_rk4_tvdlf_cada.log: lfff_2d.par $(SCHEME_DIR)/rk4_tvdlf_cada.par

include ../../test_rules.make
