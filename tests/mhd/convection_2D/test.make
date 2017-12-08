SETUP_FLAGS := -d=2
SCHEME_DIR := ../../schemes
TESTS := cv_2d_2step_tvdlf_mm.log cv_2d_2step_tvdmu_al.log		\
cv_2d_3step_hll_cada.log cv_2d_4step_hll_mc.log cv_2d_4step_hllc_ko.log	\
cv_2d_rk4_tvdlf_cada.log cv_2d_3step_hlld_cada.log

cv_2d_2step_tvdlf_mm.log: cv_2d.par $(SCHEME_DIR)/2step_tvdlf_mm.par
cv_2d_2step_tvdmu_al.log: cv_2d.par $(SCHEME_DIR)/2step_tvdmu_al.par
cv_2d_3step_hll_cada.log: cv_2d.par $(SCHEME_DIR)/3step_hll_cada.par
cv_2d_3step_hlld_cada.log: cv_2d.par $(SCHEME_DIR)/3step_hlld_cada.par
cv_2d_4step_hll_mc.log: cv_2d.par $(SCHEME_DIR)/4step_hll_mc.par
cv_2d_4step_hllc_ko.log: cv_2d.par $(SCHEME_DIR)/4step_hllc_ko.par
cv_2d_rk4_tvdlf_cada.log: cv_2d.par $(SCHEME_DIR)/rk4_tvdlf_cada.par

include ../../test_rules.make
