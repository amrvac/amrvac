TEST_DIRS := rho/vac_16 rho/vac_14x14 rho/vac_16x16 rho/vac_16x16x16	\
hd/jetCloudInteractions hd/RM2D hd/RT3D mhd/doubleGEM mhd/promRT

# Directories with altered names (useful for cleaning)
CLEAN_TEST := $(TEST_DIRS:%=clean-%)

.PHONY: all clean $(TEST_DIRS) $(CLEAN_TEST)

all: $(TEST_DIRS)

$(TEST_DIRS):
	@$(MAKE) -C $@ -f test.make

clean:	$(CLEAN_TEST)

$(CLEAN_TEST):
	$(MAKE) -C $(@:clean-%=%) clean
