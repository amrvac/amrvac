SETUP_FLAGS := -d=1 -v=1

TESTS := whiski_1d_whiski.log whiski_1dwhiski_reflection_whiski.log \
	whiski_1dwhiski_damping_whiski.log whiski_1dwhiski_negative_whiski.log

include ../../test_rules.make

whiski_1d_whiski.log: amrvac.par
whiski_1dwhiski_reflection_whiski.log: amrvac.par reflection.par
whiski_1dwhiski_damping_whiski.log: amrvac.par damping.par
whiski_1dwhiski_negative_whiski.log: amrvac.par negative.par
