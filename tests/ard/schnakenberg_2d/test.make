SETUP_FLAGS := -d=2
TESTS := imex_nonlin_adv.log

include ../../test_rules.make

# Dependency rules for the tests
imex_nonlin_adv.log: imex_nonlin_adv.par
