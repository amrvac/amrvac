SETUP_FLAGS := -d=2
SCHEME_DIR := ../../schemes

SCHEMES := 3step_tvdlf_mm 3step_hll_cada 3step_hlld_cada 3step_fd_mp5

TESTS := $(SCHEMES:%=twofl_khi_2d_%.log)

include ../../test_rules.make

# Generate dependency rules for the tests
$(foreach s, $(SCHEMES),\
	$(eval $(s:%=twofl_khi_2d_%.log): twofl_khi_2d.par $(SCHEME_DIR)/$(s).par))

