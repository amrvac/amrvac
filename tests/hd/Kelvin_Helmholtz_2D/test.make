SETUP_FLAGS := -d=2
SCHEME_DIR := ../../schemes
SCHEMES := 3step_hll_cada 4step_hll_mc 4step_hllc_ko rk4_tvdlf_cada	\
ssprk54_hll_mp5

TESTS := $(SCHEMES:%=kh_2d_%.log)

include ../../test_rules.make

# Generate dependency rules for the tests
$(foreach s, $(SCHEMES),\
	$(eval $(s:%=kh_2d_%.log): kh_2d.par $(SCHEME_DIR)/$(s).par))
