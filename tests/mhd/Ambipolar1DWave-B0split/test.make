SETUP_FLAGS := -d=1
TESTS := ambi.log \

ambi.log: amrvac.par

include ../../test_rules.make
