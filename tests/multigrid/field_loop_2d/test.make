SETUP_FLAGS := -d=2
TESTS := test_field_loop_2d.log

include ../../test_rules.make

# Dependency rules for the tests
test_field_loop_2d.log: test_field_loop_2d.par
