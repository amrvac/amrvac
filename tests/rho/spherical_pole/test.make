TESTS := uflow.log blob.log

include ../../test_rules.make

# Generate dependency rules for the tests
uflow.log: uflow.par
blob.log: blob.par
