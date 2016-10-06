SETUP_FLAGS := -d=22 -phi=0 -z=0 -g=16,16 -p=hd -eos=default -nf=1
SCHEME_DIR := ../../schemes
TESTS := jetcl_2d_ssprk54_hllc_mp5.log

jetcl_2d_ssprk54_hllc_mp5.log: PARS = jetcl_2d_ssprk54_hllc_mp5.par

include ../../test_rules.make
