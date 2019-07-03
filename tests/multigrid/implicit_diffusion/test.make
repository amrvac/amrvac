SETUP_FLAGS := -d=2
TESTS := test_impl_diff_2d.log

include ../../test_rules.make

# Dependency rules for the tests
test_impl_diff_2d.log: test_impl_diff_2d.par
