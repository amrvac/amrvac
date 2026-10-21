SETUP_FLAGS := -d=1 -v=1
SCHEME_DIR := ../../schemes

SCHEMES := dummy

SCHEMES_FIXED := 

SCHEMES_BLOW := 

TESTS := $(SCHEMES:%=rhd_wave_%.log) $(SCHEMES_FIXED:%=rhd_wave_%_fixes_on.log)

include ../../test_rules.make

# Generate dependency rules for the tests
$(foreach s, $(SCHEMES),\
	$(eval $(s:%=rhd_wave_%.log): rhd_wave.par $(SCHEME_DIR)/$(s).par))

$(foreach s, $(SCHEMES_FIXED),\
	$(eval $(s:%=rhd_wave_%_fixes_on.log): rhd_wave.par $(SCHEME_DIR)/$(s).par $(SCHEME_DIR)/fixes_on.par))
