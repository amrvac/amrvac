SETUP_FLAGS := -d=2 -v=2

TESTS := test_heatcool_2D_high.log test_heatcool_2D_low.log

include ../../test_rules.make

# Dependency rules for the tests
test_heatcool_2D_high.log: high_test.par
test_heatcool_2D_low.log: low_test.par
