SETUP_FLAGS := -d=1 -v=3
SCHEME_DIR := ../../schemes

SCHEMES := IMEX_SP_tvdlf_ko IMEX_ARS3_hll_cada3 IMEX_Midp_hll_w5

TESTS := $(SCHEMES:%=Torrilhon_%.log)

include ../../test_rules.make

# Generate dependency rules for the tests
$(foreach s, $(SCHEMES),\
	$(eval $(s:%=Torrilhon_%.log): Torrilhon.par $(SCHEME_DIR)/$(s).par))
