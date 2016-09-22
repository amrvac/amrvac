SETUP_FLAGS := -d=23 -phi=0 -z=0 -g=18,18 -p=mhd -eos=default
SCHEME_DIR := ../../schemes
TESTS := doubleGEM_25D_A.log doubleGEM_25D_B.log

doubleGEM_25D_A.log: PARS = doubleGEM_25D_A.par
doubleGEM_25D_B.log: PARS = doubleGEM_25D_B.par

include ../../test_rules.make
