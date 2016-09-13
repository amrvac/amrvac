TESTS := vac2d_mxnest6 testvac22_ssprk43_fd_mp5 testvac22_ssprk54_fd_mp5
SETUP_FLAGS := -d=22 -g=16,16 -p=rho

include ../../test_rules.make
