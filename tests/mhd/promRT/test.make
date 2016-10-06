SETUP_FLAGS := -d=23 -phi=0 -z=0 -g=14,14 -p=mhd -eos=default -nf=1
SCHEME_DIR := ../../schemes
TESTS := promRT_25D.log

promRT_25D.log: PARS = promRT_25D.par

include ../../test_rules.make
