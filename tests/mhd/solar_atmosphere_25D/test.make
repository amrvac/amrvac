SETUP_FLAGS := -d=23 -phi=0 -z=0 -p=mhd
TESTS := solar_atm_25D.log

solar_atm_25D.log: PARS = solar_atm_25D.par

include ../../test_rules.make
