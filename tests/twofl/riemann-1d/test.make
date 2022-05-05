SETUP_FLAGS := -d=1
TESTS := shock.log 

include ../../test_rules.make

shock.log: amrvac.par

