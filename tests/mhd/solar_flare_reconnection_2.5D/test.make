SETUP_FLAGS := -d=2
SCHEME_DIR := ../../schemes
TESTS := sfr_2.5d_2step_tvdlf_mm.log sfr_2.5d_2step_tvdmu_al.log		\
sfr_2.5d_3step_hll_cada.log sfr_2.5d_4step_hll_mc.log sfr_2.5d_4step_hllc_ko.log	\
sfr_2.5d_rk4_tvdlf_cada.log sfr_2.5d_3step_hlld_cada.log sfr_2.5d_2step_hll_vl.log

sfr_2.5d_2step_tvdlf_mm.log: PARS = sfr_2.5d.par $(SCHEME_DIR)/2step_tvdlf_mm.par
sfr_2.5d_2step_tvdmu_al.log: PARS = sfr_2.5d.par $(SCHEME_DIR)/2step_tvdmu_al.par
sfr_2.5d_3step_hll_cada.log: PARS = sfr_2.5d.par $(SCHEME_DIR)/3step_hll_cada.par
sfr_2.5d_3step_hlld_cada.log: PARS = sfr_2.5d.par $(SCHEME_DIR)/3step_hlld_cada.par
sfr_2.5d_4step_hll_mc.log: PARS = sfr_2.5d.par $(SCHEME_DIR)/4step_hll_mc.par
sfr_2.5d_4step_hllc_ko.log: PARS = sfr_2.5d.par $(SCHEME_DIR)/4step_hllc_ko.par
sfr_2.5d_rk4_tvdlf_cada.log: PARS = sfr_2.5d.par $(SCHEME_DIR)/rk4_tvdlf_cada.par
sfr_2.5d_2step_hll_vl.log: PARS = sfr_2.5d.par $(SCHEME_DIR)/2step_hll_vl.par

include ../../test_rules.make
