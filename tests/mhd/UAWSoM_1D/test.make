SETUP_FLAGS := -d=1 -v=1

TESTS := uawsom_1d_uawsom.log uawsom_1duawsom_reflection_uawsom.log \
	uawsom_1duawsom_damping_uawsom.log uawsom_1duawsom_negative_uawsom.log

include ../../test_rules.make

uawsom_1d_uawsom.log: amrvac.par
uawsom_1duawsom_reflection_uawsom.log: amrvac.par reflection.par
uawsom_1duawsom_damping_uawsom.log: amrvac.par damping.par
uawsom_1duawsom_negative_uawsom.log: amrvac.par negative.par
