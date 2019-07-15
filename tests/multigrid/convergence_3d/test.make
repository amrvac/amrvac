SETUP_FLAGS := -d=3
TESTS := test_mg_conv_3d.log

include ../../test_rules.make

# Dependency rules for the tests
test_mg_conv_3d.log: test_mg_conv_3d.par
