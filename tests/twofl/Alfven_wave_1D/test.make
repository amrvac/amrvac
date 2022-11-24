SETUP_FLAGS := -d=1
TESTS := implicit_coupled.log explicit_uncoupled.log 

include ../../test_rules.make

explicit_uncoupled.log: explicit_tvdlf_minmod_uncoupled.par
implicit_coupled.log: implicit_hll_cada3_coupled.par

