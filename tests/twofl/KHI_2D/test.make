SETUP_FLAGS := -d=2
TESTS := khi.log 

include ../../test_rules.make

khi.log: amrvac.par

