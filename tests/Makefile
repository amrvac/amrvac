# Makefile for automatically testing AGILE
#
# 1. You can run the hd tests with: make hd. The same applies to other physics
# modules
#
# 2. You can run them in parallel with for example: -j 4
#
# 3. You can hide the Entering/leaving directory messages with: -s

HD_DIRS := hd/KH3D hd/RT3D hd/KH3D-SMR
MHD_DIRS := mhd/blast_wave_Cartesian_3D
FFHD_DIRS := ffhd/TI3D

TEST_DIRS := $(HD_DIRS) $(MHD_DIRS) $(FFHD_DIRS)

# Directories with altered names (useful for cleaning)
CLEAN_TEST := $(TEST_DIRS:%=clean-%)

.PHONY: $(TEST_DIRS) $(CLEAN_TEST)

all: hd mhd ffhd
hd: $(HD_DIRS)
mhd: $(MHD_DIRS)
ffhd: $(FFHD_DIRS)

$(TEST_DIRS):
	@$(MAKE) -C $@ -f test.make

clean:	$(CLEAN_TEST)

$(CLEAN_TEST):
	$(MAKE) -C $(@:clean-%=%) -f test.make clean
