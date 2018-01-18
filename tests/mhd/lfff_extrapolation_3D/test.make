SETUP_FLAGS := -d=3
SCHEME_DIR := ../../schemes
SCHEMES := 3step_hll_cada

TESTS := $(SCHEMES:%=lfff_3d_%.log)

include ../../test_rules.make

# Generate dependency rules for the tests
$(foreach s, $(SCHEMES),\
	$(eval $(s:%=lfff_3d_%.log): lfff_3d.par $(SCHEME_DIR)/$(s).par))
