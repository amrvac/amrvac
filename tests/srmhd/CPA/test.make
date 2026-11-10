
TESTS := test.log

include ../../test_rules.make

# Generate dependency rules for the tests
test.log: test.par
