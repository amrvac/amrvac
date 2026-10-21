SETUP_FLAGS := -d=2 -v=3
TESTS := uawsom_solar_short.log

include ../../test_rules.make

uawsom_solar_short.log: amrvac.par short.par
