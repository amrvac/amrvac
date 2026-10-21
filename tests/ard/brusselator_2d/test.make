SETUP_FLAGS := -d=2 -v=2
SCHEME_DIR := ../../schemes

SCHEMES := IMEX_SP_hll_ko \
 IMEX_SP_hll_mm IMEX_SP_tvdlf_ko IMEX_Midp_hll_w5 \
 IMEX_222_hll_w5 IMEX_ARS3_hll_cada3 IMEX_SP_hll_w5  IMEX_Euler_hll_w5

TESTS := $(SCHEMES:%=test_brus_%.log) 

include ../../test_rules.make

# Generate dependency rules for the tests
$(foreach s, $(SCHEMES),\
	$(eval $(s:%=test_brus_%.log): test_brus.par $(SCHEME_DIR)/$(s).par))
