TESTS := jetcloud_2d
SETUP_FLAGS := -d=22 -phi=0 -z=0 -g=16,16 -p=hd -eos=default -nf=1 -u=nul

include ../../test_rules.make
