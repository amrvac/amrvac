SETUP_FLAGS := -d=2 -v=2

# The transfer solve on the exact gray radiative-equilibrium source: mod_usr stops the run if the
# flux or the heating departs from the exact solution in any cell
TESTS := gray_atm.log

include ../../test_rules.make

gray_atm.log: gray_atm.par
