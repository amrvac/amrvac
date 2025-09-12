SETUP_FLAGS := -d=1 -v=3
SCHEME_DIR := ../../schemes

SCHEMES := IMEX_SP_tvdlf_ko

TESTS := $(SCHEMES:%=Brio_Wu_%.log)

include ../../test_rules.make

# Generate dependency rules for the tests
$(foreach s, $(SCHEMES),\
	$(eval $(s:%=Brio_Wu_%.log): Brio_Wu.par $(SCHEME_DIR)/$(s).par))
