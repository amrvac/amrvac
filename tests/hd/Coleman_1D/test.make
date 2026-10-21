SETUP_FLAGS := -d=1 -v=1

# 36-run matrix: 6 Coleman (2020) Riemann tests x {FI, LTE(entropy), state} x {HHe, pureH}.
# state = the 'state' LTE method (uniform256 direct-state tables); validated to
# agree with LTE(entropy) on the conserved dynamics, differing only on the
# ionisation/temperature diagnostics at the table/interpolation-accuracy level.
# Each parfile is standalone (no scheme overlay), so the per-test rule
# just lists its own parfile as the dependency.
TEST_NUMS := 1 2 3 4 5 6
EOS_VARIANTS := FI_HHe FI_pureH LTE_HHe LTE_pureH state_HHe state_pureH

TESTS := $(foreach n,$(TEST_NUMS), \
           $(foreach v,$(EOS_VARIANTS), coleman$(n)_$(v).log))

include ../../test_rules.make

# Each {name}.log is produced by AMRVAC reading {name}.par.
$(foreach t,$(TESTS),$(eval $(t): $(t:%.log=%.par)))
