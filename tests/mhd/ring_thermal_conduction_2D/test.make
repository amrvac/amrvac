SETUP_FLAGS := -d=2
SCHEME_DIR := ../../schemes
SCHEMES := 2step_tvdlf_mm 3step_hll_cada \
 mhd_internal_e mhd_hyperbolic_thermal_conduction

TESTS := $(SCHEMES:%=rtc_2d_%.log)

include ../../test_rules.make

# Generate dependency rules for the tests
$(foreach s, $(SCHEMES),\
	$(eval $(s:%=rtc_2d_%.log): rtc_2d.par $(SCHEME_DIR)/$(s).par))
