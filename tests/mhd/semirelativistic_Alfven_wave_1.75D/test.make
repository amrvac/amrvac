SETUP_FLAGS := -d=1
SCHEME_DIR := ../../schemes
SCHEMES := 2step_tvdlf_mm 3step_hll_cada 4step_hll_mc

TESTS := $(SCHEMES:%=semirelati_Alfven_wave_1d_%.log)

include ../../test_rules.make

# Generate dependency rules for the tests
$(foreach s, $(SCHEMES),\
	$(eval $(s:%=semirelati_Alfven_wave_1d_%.log): sraw_1.75d.par $(SCHEME_DIR)/$(s).par))
