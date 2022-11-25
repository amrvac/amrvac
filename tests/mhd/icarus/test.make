SETUP_FLAGS := -d=3
SCHEME_DIR := ../../schemes
SCHEMES := 2step_tvdlf_mm

TESTS := $(SCHEMES:%=icarus_3d_%.log)

include ../../test_rules.make

# Generate dependency rules for the tests
$(foreach s, $(SCHEMES),\
	$(eval $(s:%=icarus_3d_%.log): amrvac.par $(SCHEME_DIR)/$(s).par))
