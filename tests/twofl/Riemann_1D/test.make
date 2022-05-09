SETUP_FLAGS := -d=1
SCHEME_DIR := ../../schemes

SCHEMES := 3step_tvdlf_mm 3step_hll_cada 3step_hlld_cada 3step_fd_mp5

TESTS := $(SCHEMES:%=twofl_shock_1d_%.log)

include ../../test_rules.make

# Generate dependency rules for the tests
$(foreach s, $(SCHEMES),\
	$(eval $(s:%=twofl_shock_1d_%.log): twofl_R1d.par $(SCHEME_DIR)/$(s).par))
