SETUP_FLAGS := -d=2
SCHEME_DIR := ../../schemes
TESTS := rip_2.5d_2step_tvdlf_mm.log rip_2.5d_2step_tvdmu_al.log		\
rip_2.5d_3step_hll_cada.log rip_2.5d_4step_hll_mc.log rip_2.5d_4step_hllc_ko.log	\
rip_2.5d_rk4_tvdlf_cada.log rip_2.5d_3step_hlld_cada.log

rip_2.5d_2step_tvdlf_mm.log: rip_2.5d.par $(SCHEME_DIR)/2step_tvdlf_mm.par
rip_2.5d_2step_tvdmu_al.log: rip_2.5d.par $(SCHEME_DIR)/2step_tvdmu_al.par
rip_2.5d_3step_hll_cada.log: rip_2.5d.par $(SCHEME_DIR)/3step_hll_cada.par
rip_2.5d_3step_hlld_cada.log: rip_2.5d.par $(SCHEME_DIR)/3step_hlld_cada.par
rip_2.5d_4step_hll_mc.log: rip_2.5d.par $(SCHEME_DIR)/4step_hll_mc.par
rip_2.5d_4step_hllc_ko.log: rip_2.5d.par $(SCHEME_DIR)/4step_hllc_ko.par
rip_2.5d_rk4_tvdlf_cada.log: rip_2.5d.par $(SCHEME_DIR)/rk4_tvdlf_cada.par

include ../../test_rules.make
