SETUP_FLAGS := -d=2
SCHEME_DIR := ../../schemes
TESTS := rotor_2d_2step_tvdlf_mm.log \
rotor_2d_3step_hll_cada.log rotor_2d_4step_hll_mc.log rotor_2d_4step_hllc_ko.log	\
rotor_2d_rk4_tvdlf_cada.log rotor_2d_ssprk54_fd_mp5.log rotor_2d_ssprk54_hll_mp5.log \
rotor_2d_3step_hlld_cada.log rotor_2d_ssprk54_hlld_mp5.log

rotor_2d_2step_tvdlf_mm.log: PARS = rotor_2d.par $(SCHEME_DIR)/2step_tvdlf_mm.par
rotor_2d_3step_hll_cada.log: PARS = rotor_2d.par $(SCHEME_DIR)/3step_hll_cada.par
rotor_2d_3step_hlld_cada.log: PARS = rotor_2d.par $(SCHEME_DIR)/3step_hlld_cada.par
rotor_2d_4step_hll_mc.log: PARS = rotor_2d.par $(SCHEME_DIR)/4step_hll_mc.par
rotor_2d_4step_hllc_ko.log: PARS = rotor_2d.par $(SCHEME_DIR)/4step_hllc_ko.par
rotor_2d_rk4_tvdlf_cada.log: PARS = rotor_2d.par $(SCHEME_DIR)/rk4_tvdlf_cada.par
rotor_2d_ssprk54_fd_mp5.log: PARS = rotor_2d.par $(SCHEME_DIR)/ssprk54_fd_mp5.par
rotor_2d_ssprk54_hll_mp5.log: PARS = rotor_2d.par $(SCHEME_DIR)/ssprk54_hll_mp5.par
rotor_2d_ssprk54_hlld_mp5.log: PARS = rotor_2d.par $(SCHEME_DIR)/ssprk54_hlld_mp5.par

include ../../test_rules.make
