SETUP_FLAGS := -d=2
TESTS := solar_atm_25D.log

solar_atm_25D.log: PARS = solar_atm_25D.par

include ../../test_rules.make
