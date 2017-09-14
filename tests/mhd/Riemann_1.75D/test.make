SETUP_FLAGS := -d=1
SCHEME_DIR := ../../schemes
TESTS := R_1d_2step_tvdlf_mm.log R_1d_2step_tvdmu_al.log \
R_1d_3step_hll_cada.log R_1d_4step_hll_mc.log R_1d_4step_hllc_ko.log	\
R_1d_rk4_tvdlf_cada.log R_1d_ssprk54_fd_mp5.log R_1d_ssprk54_hll_mp5.log \
R_1d_3step_hlld_cada.log R_1d_ssprk54_hlld_mp5.log

R_1d_2step_tvdlf_mm.log: PARS = R_1d.par $(SCHEME_DIR)/2step_tvdlf_mm.par
R_1d_2step_tvdmu_al.log: PARS = R_1d.par $(SCHEME_DIR)/2step_tvdmu_al.par
R_1d_3step_hll_cada.log: PARS = R_1d.par $(SCHEME_DIR)/3step_hll_cada.par
R_1d_3step_hlld_cada.log: PARS = R_1d.par $(SCHEME_DIR)/3step_hlld_cada.par
R_1d_4step_hll_mc.log: PARS = R_1d.par $(SCHEME_DIR)/4step_hll_mc.par
R_1d_4step_hllc_ko.log: PARS = R_1d.par $(SCHEME_DIR)/4step_hllc_ko.par
R_1d_rk4_tvdlf_cada.log: PARS = R_1d.par $(SCHEME_DIR)/rk4_tvdlf_cada.par
R_1d_ssprk54_fd_mp5.log: PARS = R_1d.par $(SCHEME_DIR)/ssprk54_fd_mp5.par
R_1d_ssprk54_hll_mp5.log: PARS = R_1d.par $(SCHEME_DIR)/ssprk54_hll_mp5.par
R_1d_ssprk54_hlld_mp5.log: PARS = R_1d.par $(SCHEME_DIR)/ssprk54_hlld_mp5.par

include ../../test_rules.make
