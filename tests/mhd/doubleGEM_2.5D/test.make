SETUP_FLAGS := -d=2
SCHEME_DIR := ../../schemes
TESTS := dgem_2d_2step_tvdlf_mm.log dgem_2d_2step_tvdmu_al.log		\
dgem_2d_3step_hll_cada.log dgem_2d_4step_hll_mc.log dgem_2d_4step_hllc_ko.log	\
dgem_2d_rk4_tvdlf_cada.log

dgem_2d_2step_tvdlf_mm.log: PARS = dgem_2d.par $(SCHEME_DIR)/2step_tvdlf_mm.par
dgem_2d_2step_tvdmu_al.log: PARS = dgem_2d.par $(SCHEME_DIR)/2step_tvdmu_al.par
dgem_2d_3step_hll_cada.log: PARS = dgem_2d.par $(SCHEME_DIR)/3step_hll_cada.par
dgem_2d_4step_hll_mc.log: PARS = dgem_2d.par $(SCHEME_DIR)/4step_hll_mc.par
dgem_2d_4step_hllc_ko.log: PARS = dgem_2d.par $(SCHEME_DIR)/4step_hllc_ko.par
dgem_2d_rk4_tvdlf_cada.log: PARS = dgem_2d.par $(SCHEME_DIR)/rk4_tvdlf_cada.par

include ../../test_rules.make
