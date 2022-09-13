SETUP_FLAGS := -d=1
SCHEME_DIR := ../../schemes

SCHEMES := IMEX_SP_tvdlf_ko IMEX_Trap_hll_w5
# SCHEMES := IMEX_Euler_hll_w5 IMEX_Midp_hll_w5 IMEX_SP_hll_ko \
# IMEX_SP_hll_mm IMEX_SP_hll_w5 IMEX_SP_tvdlf_ko IMEX_Trap_hll_w5

TESTS := $(SCHEMES:%=rhd_shock_%.log)

include ../../test_rules.make

# Generate dependency rules for the tests
$(foreach s, $(SCHEMES),\
	$(eval $(s:%=rhd_shock_%.log): rhd_shock.par $(SCHEME_DIR)/$(s).par))
