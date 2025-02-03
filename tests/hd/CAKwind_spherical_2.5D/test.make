SETUP_FLAGS := -d=2 -v=3
SCHEME_DIR := ../../schemes

SCHEMES := 3step_hll_ko 4step_hll_mc

TESTS := $(SCHEMES:%=cak2.5d_%.log)

include ../../test_rules.make

# Generate dependency rules for the tests, which are used to run them
$(foreach s, $(SCHEMES),\
    $(eval $(s:%=cak2.5d_%.log): cak_2.5d.par $(SCHEME_DIR)/$(s).par))
