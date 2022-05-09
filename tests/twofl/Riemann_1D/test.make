SETUP_FLAGS := -d=1
SCHEME_DIR := ../../schemes

SCHEMES := 2step_tvdlf_mm 3step_hll_cada \
4step_hll_mc rk4_tvdlf_cada ssprk54_fd_mp5 ssprk54_hll_mp5

TESTS := $(SCHEMES:%=twofl_shock_1d_%.log)

include ../../test_rules.make

# Generate dependency rules for the tests
$(foreach s, $(SCHEMES),\
	$(eval $(s:%=twofl_shock_1d_%.log): twofl_R1d.par $(SCHEME_DIR)/$(s).par))
