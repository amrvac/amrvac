SETUP_FLAGS := -d=2
SCHEME_DIR := ../../schemes
TESTS := kh_2d_3step_hll_cada.log kh_2d_4step_hll_mc.log	\
kh_2d_4step_hllc_ko.log kh_2d_rk4_tvdlf_cada.log kh_2d_ssprk54_hll_mp5.log

kh_2d_3step_hll_cada.log: PARS = kh_2d.par $(SCHEME_DIR)/3step_hll_cada.par
kh_2d_4step_hll_mc.log: PARS = kh_2d.par $(SCHEME_DIR)/4step_hll_mc.par
kh_2d_4step_hllc_ko.log: PARS = kh_2d.par $(SCHEME_DIR)/4step_hllc_ko.par
kh_2d_rk4_tvdlf_cada.log: PARS = kh_2d.par $(SCHEME_DIR)/rk4_tvdlf_cada.par
kh_2d_ssprk54_hll_mp5.log: PARS = kh_2d.par $(SCHEME_DIR)/ssprk54_hll_mp5.par

include ../../test_rules.make
