SETUP_FLAGS := -d=1
TESTS := implicit-coupled.log explicit-uncoupled.log 

include ../../test_rules.make

explicit-uncoupled.log: explicit-tvdlf-minmod-uncoupled.par
implicit-coupled.log: implicit-hll-cada3-coupled.par

