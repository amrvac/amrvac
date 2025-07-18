# Makefile for automatically testing AGILE
#
# 1. You can run the hd tests with: make hd. The same applies to other physics
# modules
#
# 2. You can run them in parallel with for example: -j 4
#
# 3. You can hide the Entering/leaving directory messages with: -s

HD_DIRS := hd/KH3D hd/RT3D

TEST_DIRS := $(HD_DIRS)

# Directories with altered names (useful for cleaning)
CLEAN_TEST := $(TEST_DIRS:%=clean-%)

.PHONY: $(TEST_DIRS) $(CLEAN_TEST)

all: hd
hd: $(HD_DIRS)

$(TEST_DIRS):
	@$(MAKE) -C $@ -f test.make

clean:	$(CLEAN_TEST)

$(CLEAN_TEST):
	$(MAKE) -C $(@:clean-%=%) -f test.make clean
