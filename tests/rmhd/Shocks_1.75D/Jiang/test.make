SETUP_FLAGS := -d=1 -v=3
SCHEME_DIR := ../../../schemes

SCHEMES_FIXED := dummy

TESTS := $(SCHEMES_FIXED:%=rhd_jiang_%_fixes_on.log)

include ../../../test_rules.make

# Generate dependency rules for the tests
$(foreach s, $(SCHEMES_FIXED),\
	$(eval $(s:%=rhd_jiang_%_fixes_on.log): rhd_jiang.par $(SCHEME_DIR)/$(s).par $(SCHEME_DIR)/fixes_on.par))
