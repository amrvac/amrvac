TESTS := testvac22_1step_tvd testvac22_4step_hll_mc testvac22_2step_tvdlf_mm	\
testvac22_4step_tvdlf_cada testvac22_2step_tvdmu_al testvac22_rk4_tvdlf_cada	\
testvac22_3step_hll_cada testvac22_4step_hllc_ko

SETUP_FLAGS := -d=22 -g=14,14 -p=rho

include ../../test_rules.make
