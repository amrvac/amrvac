SETUP_FLAGS := -d=1 -v=1
SCHEME_DIR := ../../schemes

SCHEMES := IMEX_222_hll_w5 IMEX_232_hll_cada3 IMEX_Euler_hll_w5 \
 IMEX_SP_hll_mm IMEX_SP_tvdlf_ko IMEX_Midp_hll_w5 \
 IMEX_SP_hll_ko IMEX_Trap_hll_w5 IMEX_SP_hll_w5 IMEX_CB3a_hll_ko \
 amr_3level fld_minerbo fld_pomraning

SCHEMES_FIXED := IMEX_ARS3_hll_cada3 IMEX_ARS3_hll_ko

SCHEMES_BLOW := 

TESTS := $(SCHEMES:%=rhd_pulse_%.log) $(SCHEMES_FIXED:%=rhd_pulse_%_fixes_on.log)

include ../../test_rules.make

# Generate dependency rules for the tests
$(foreach s, $(SCHEMES),\
	$(eval $(s:%=rhd_pulse_%.log): rhd_pulse.par $(SCHEME_DIR)/$(s).par))

$(foreach s, $(SCHEMES_FIXED),\
	$(eval $(s:%=rhd_pulse_%_fixes_on.log): rhd_pulse.par $(SCHEME_DIR)/$(s).par $(SCHEME_DIR)/fixes_on.par))
