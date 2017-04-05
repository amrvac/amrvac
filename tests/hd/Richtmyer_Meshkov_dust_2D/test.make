SETUP_FLAGS := -d=2
TESTS := RM2D_dust_sticking.log RM2D_dust_Kwok.log

RM2D_dust_sticking.log: PARS = RM2D_dust_sticking.par
RM2D_dust_Kwok.log: PARS = RM2D_dust_Kwok.par

include ../../test_rules.make
