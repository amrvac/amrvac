SETUP_FLAGS := -d=2
TESTS := test_mg_conv_2d.log

include ../../test_rules.make

# Dependency rules for the tests
test_mg_conv_2d.log: test_mg_conv_2d.par
