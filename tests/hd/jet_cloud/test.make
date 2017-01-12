SETUP_FLAGS := -d=22
SCHEME_DIR := ../../schemes
TESTS := jetcl_2d_ssprk54_hllc_mp5.log

jetcl_2d_ssprk54_hllc_mp5.log: PARS = jetcl_2d_ssprk54_hllc_mp5.par

include ../../test_rules.make
