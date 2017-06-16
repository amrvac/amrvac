SETUP_FLAGS := -d=2
TESTS := solar_atm_25D.log solar_atm_B0split.log

solar_atm_25D.log: PARS = solar_atm_25D.par
solar_atm_B0split.log: PARS = solar_atm_B0split.par

include ../../test_rules.make
