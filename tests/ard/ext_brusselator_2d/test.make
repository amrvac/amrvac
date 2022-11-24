SETUP_FLAGS := -d=1
TESTS := ext_bruselator_2d.log 

include ../../test_rules.make

ext_bruselator_2d.log: spirals.par

