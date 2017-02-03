SETUP_FLAGS := -d=22
SCHEME_DIR := ../../schemes
TESTS := kh_2d_3step_hll_cada.log kh_2d_4step_hll_mc.log	\
kh_2d_4step_hllc_ko.log kh_2d_rk4_tvdlf_cada.log

kh_2d_3step_hll_cada.log: PARS = test.par $(SCHEME_DIR)/3step_hll_cada.par
kh_2d_4step_hll_mc.log: PARS = test.par $(SCHEME_DIR)/4step_hll_mc.par
kh_2d_4step_hllc_ko.log: PARS = test.par $(SCHEME_DIR)/4step_hllc_ko.par
kh_2d_rk4_tvdlf_cada.log: PARS = test.par $(SCHEME_DIR)/rk4_tvdlf_cada.par

include ../../test_rules.make
