TESTS := uflow.log blast.log

include ../../test_rules.make

# Generate dependency rules for the tests. Both share one build of ./agile:
# they agree on every compile-time parameter and differ only in their par file.
uflow.log: uflow.par
blast.log: blast.par
