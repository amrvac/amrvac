TESTS := uflow.log

include ../../test_rules.make

# Generate dependency rules for the tests
uflow.log: uflow.par
