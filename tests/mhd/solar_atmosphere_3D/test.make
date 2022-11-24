SETUP_FLAGS := -d=3
SCHEME_DIR := ../../schemes
SCHEMES := 3step_hll_cada_ct 3step_hll_cada 3step_hlld_cada_ct semirelati_mhd

TESTS := $(SCHEMES:%=solar_atm_3D_%.log)

include ../../test_rules.make

# Generate dependency rules for the tests
$(foreach s, $(SCHEMES),\
	$(eval $(s:%=solar_atm_3D_%.log): solar_atm_3D.par $(SCHEME_DIR)/$(s).par))
