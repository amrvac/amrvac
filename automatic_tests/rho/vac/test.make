TESTS := vac2d vac2d_mxnest6
SETUP_FLAGS := -d=22 -g=16,16 -p=rho -u=testrho

include ../../test_rules.make
