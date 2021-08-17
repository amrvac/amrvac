SETUP_FLAGS := -d=3
SCHEME_DIR := ../../schemes
SCHEMES := 3step_hll_cada 3step_hll_cada_ct 4step_cd4_cada 3step_fd_cada 1step_cd4_mp5

TESTS := $(SCHEMES:%=mb_3d_%.log)

include ../../test_rules.make

# Generate dependency rules for the tests
$(foreach s, $(SCHEMES),\
	$(eval $(s:%=mb_3d_%.log): mb_3d.par $(SCHEME_DIR)/$(s).par))
