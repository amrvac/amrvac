SETUP_FLAGS := -d=1 -v=1
SCHEME_DIR := ../../schemes

SCHEMES := dummy

TESTS := $(SCHEMES:%=ShuO_1D_weno7_%.log)

include ../../test_rules.make

# Generate dependency rules for the tests
$(foreach s, $(SCHEMES),\
	$(eval $(s:%=ShuO_1D_weno7_%.log): test.par $(SCHEME_DIR)/$(s).par))
