SETUP_FLAGS := -d=2
TESTS := test_schnakenberg_imex.log test_schnakenberg_split.log

include ../../test_rules.make

# Dependency rules for the tests
test_schnakenberg_imex.log: test_schnakenberg_imex.par
test_schnakenberg_split.log: test_schnakenberg_split.par
