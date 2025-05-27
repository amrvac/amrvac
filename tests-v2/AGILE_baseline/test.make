SETUP_FLAGS := -d=2
TESTS := khi.log \

khi.log: amrvac_test.par

include ../../test_rules.make
