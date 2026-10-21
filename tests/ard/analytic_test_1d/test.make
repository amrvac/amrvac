SETUP_FLAGS := -d=1 -v=1

TESTS := test_ana.log

include ../../test_rules.make

test_ana.log: test_ana.par
