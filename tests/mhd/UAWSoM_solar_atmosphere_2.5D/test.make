SETUP_FLAGS := -d=2 -v=3
TESTS := whiski_solar_short.log

include ../../test_rules.make

whiski_solar_short.log: amrvac.par short.par
