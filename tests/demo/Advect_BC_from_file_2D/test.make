SETUP_FLAGS := -d=2 -v=2
SCHEME_DIR := ../../schemes

SCHEMES := 3step_tvdlf_woodward 4step_hll_mc

TESTS := $(SCHEMES:%=advect_alfven_test_%.log)

include ../../test_rules.make

# Generate dependency rules for the tests
$(foreach s, $(SCHEMES),\
	$(eval $(s:%=advect_alfven_test_%.log): test.par $(SCHEME_DIR)/$(s).par))
