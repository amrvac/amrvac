SETUP_FLAGS := -d=1 -v=1
SCHEME_DIR := ../../schemes

SCHEMES := IMEX_232_hll_cada3 

SCHEMES_FIXED := IMEX_222_hll_w5 IMEX_Midp_hll_w5 IMEX_CB3a_hll_ko IMEX_ARS3_hll_cada3 \
         IMEX_SP_hll_mm IMEX_SP_tvdlf_ko IMEX_Trap_hll_w5 amr_3level

SCHEMES_FIXED_SLOW := IMEX_Euler_hll_w5 IMEX_SP_hll_w5

TESTS := $(SCHEMES:%=rhd_shock_%.log) $(SCHEMES_FIXED:%=rhd_shock_%_fixes_on.log) $(SCHEMES_FIXED_SLOW:%=rhd_shock_%_fixes_on_slow_fld.log)

include ../../test_rules.make

# Generate dependency rules for the tests
$(foreach s, $(SCHEMES),\
	$(eval $(s:%=rhd_shock_%.log): rhd_shock.par $(SCHEME_DIR)/$(s).par))

$(foreach s, $(SCHEMES_FIXED),\
	$(eval $(s:%=rhd_shock_%_fixes_on.log): rhd_shock.par $(SCHEME_DIR)/$(s).par $(SCHEME_DIR)/fixes_on.par))

$(foreach s, $(SCHEMES_FIXED_SLOW),\
	$(eval $(s:%=rhd_shock_%_fixes_on_slow_fld.log): rhd_shock.par $(SCHEME_DIR)/$(s).par $(SCHEME_DIR)/fixes_on_slow_fld.par))
