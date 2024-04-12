SETUP_FLAGS := -d=3
SCHEME_DIR := ../../schemes
SCHEMES := 1step_cd4_mp5 3step_hll_cada 3step_fd_cada

TESTS := $(SCHEMES:%=saa_3d_%.log)

include ../../test_rules.make

# Generate dependency rules for the tests
$(foreach s, $(SCHEMES),\
	$(eval $(s:%=saa_3d_%.log): saa_3d.par $(SCHEME_DIR)/$(s).par))
