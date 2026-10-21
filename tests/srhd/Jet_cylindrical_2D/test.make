SETUP_FLAGS := -d=2 -v=2
SCHEME_DIR := ../../schemes

SCHEMES := dummy

TESTS := $(SCHEMES:%=SRHD_Jet_%.log)

include ../../test_rules.make

# Generate dependency rules for the tests
$(foreach s, $(SCHEMES),\
	$(eval $(s:%=SRHD_Jet_%.log): test.par $(SCHEME_DIR)/$(s).par))
