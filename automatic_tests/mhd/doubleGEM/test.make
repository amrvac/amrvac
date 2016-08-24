TESTS := doubleGEM_A doubleGEM_B
SETUP_FLAGS := -d=23 -phi=0 -z=0 -g=18,18 -p=mhd -eos=default

include ../../test_rules.make
