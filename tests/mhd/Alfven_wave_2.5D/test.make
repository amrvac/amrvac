SETUP_FLAGS := -d=2
SCHEME_DIR := ../../schemes
TESTS := aw_2d_2step_tvdlf_mm.log aw_2d_2step_tvdmu_al.log		\
aw_2d_3step_hll_cada.log aw_2d_4step_hll_mc.log aw_2d_4step_hllc_ko.log	\
aw_2d_rk4_tvdlf_cada.log aw_2d_ssprk54_fd_mp5.log aw_2d_ssprk54_hll_mp5.log \
aw_2d_3step_hlld_cada.log aw_2d_ssprk54_hlld_mp5.log

aw_2d_2step_tvdlf_mm.log: aw_2d.par $(SCHEME_DIR)/2step_tvdlf_mm.par
aw_2d_2step_tvdmu_al.log: aw_2d.par $(SCHEME_DIR)/2step_tvdmu_al.par
aw_2d_3step_hll_cada.log: aw_2d.par $(SCHEME_DIR)/3step_hll_cada.par
aw_2d_4step_hll_mc.log: aw_2d.par $(SCHEME_DIR)/4step_hll_mc.par
aw_2d_4step_hllc_ko.log: aw_2d.par $(SCHEME_DIR)/4step_hllc_ko.par
aw_2d_rk4_tvdlf_cada.log: aw_2d.par $(SCHEME_DIR)/rk4_tvdlf_cada.par
aw_2d_ssprk54_fd_mp5.log: aw_2d.par $(SCHEME_DIR)/ssprk54_fd_mp5.par
aw_2d_ssprk54_hll_mp5.log: aw_2d.par $(SCHEME_DIR)/ssprk54_hll_mp5.par
aw_2d_3step_hlld_cada.log: aw_2d.par $(SCHEME_DIR)/3step_hlld_cada.par
aw_2d_ssprk54_hlld_mp5.log: aw_2d.par $(SCHEME_DIR)/ssprk54_hlld_mp5.par

include ../../test_rules.make
