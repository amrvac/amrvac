SETUP_FLAGS := -d=2
SCHEME_DIR := ../../schemes
SCHEMES := 3step_hll_cada_ct 3step_hll_cada 4step_hll_mc

TESTS := $(SCHEMES:%=solar_wind_25D_%.log)

include ../../test_rules.make

# Generate dependency rules for the tests
$(foreach s, $(SCHEMES),\
	$(eval $(s:%=solar_wind_25D_%.log): sw_25d.par $(SCHEME_DIR)/$(s).par))
