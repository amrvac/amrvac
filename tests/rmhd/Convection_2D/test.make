SETUP_FLAGS := -d=2 -v=2
SCHEME_DIR := ../../schemes

SCHEMES := IMEX_ARS3_hll_cada3 IMEX_Midp_hll_w5 IMEX_SP_tvdlf_ko

TESTS := $(SCHEMES:%=Convection2d_%.log)

include ../../test_rules.make

# Generate dependency rules for the tests
$(foreach s, $(SCHEMES),\
	$(eval $(s:%=Convection2d_%.log): Convection2d.par $(SCHEME_DIR)/$(s).par))
