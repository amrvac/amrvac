SETUP_FLAGS := -d=2 -v=2
SCHEME_DIR := ../../schemes

SCHEMES := 3step_hll_cada 4step_hll_mc 3step_tvdlf_superbee 5step_tvdlf_vl

TESTS := $(SCHEMES:%=burgers_2d_%.log)

include ../../test_rules.make

# Generate dependency rules for the tests
$(foreach s, $(SCHEMES),\
	$(eval $(s:%=burgers_2d_%.log): burgers_2d.par $(SCHEME_DIR)/$(s).par))
