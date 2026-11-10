# Makefile for automatically testing AGILE
#
# 1. You can run the hd tests with: make hd. The same applies to other physics
# modules
#
# 2. You can run them in parallel with for example: -j 4
#
# 3. You can hide the Entering/leaving directory messages with: -s

# 4. Each physics module's directories are split into a basic suite, run by
# default and on every push/PR in CI, and an extended suite. Select with
# SUITE=basic (default), SUITE=extended or SUITE=all, e.g.
#    make SUITE=all -j4 hd
# Individual directories (make hd/spherical) run regardless of SUITE.
# make list prints the directories the current SUITE selects.
#
# A directory is one build, and builds dominate the cost, so the split is per
# directory. The basic suite keeps one Cartesian case per physics plus the
# curvilinear cases whose checks nothing else provides.

SUITE ?= basic

HD_BASIC := hd/KH3D-SMR hd/spherical hd/cylindrical_pole hd/log_spherical \
            hd/log_cylindrical_pole
HD_EXTENDED := hd/KH3D hd/KH3D-SMR-integrators hd/RT3D hd/KH3D_tracer \
               hd/sphere_advection hd/TRML hd/cloud_crushing hd/cylindrical \
               hd/spherical_pole hd/log_cylindrical hd/log_spherical_pole
MHD_BASIC := mhd/Orszag-Tang_3D mhd/cylindrical_pole
MHD_EXTENDED := mhd/blast_wave_Cartesian_3D mhd/BrioWu mhd/spherical mhd/cylindrical \
                mhd/spherical_pole mhd/log_spherical mhd/log_cylindrical
FFHD_BASIC := ffhd/TI3D ffhd/spherical_pole
FFHD_EXTENDED := ffhd/spherical ffhd/cylindrical ffhd/cylindrical_pole \
                 ffhd/log_spherical ffhd/log_cylindrical
SRHD_BASIC := srhd/shockTube srhd/cylindrical
SRHD_EXTENDED := srhd/jet_3D srhd/spherical srhd/spherical_pole \
                 srhd/cylindrical_pole srhd/log_spherical srhd/log_cylindrical
SRMHD_BASIC := srmhd/shock_tube srmhd/spherical_pole
SRMHD_EXTENDED := srmhd/CPA srmhd/spherical srmhd/cylindrical \
                  srmhd/cylindrical_pole srmhd/log_spherical srmhd/log_cylindrical

HD_DIRS := $(HD_BASIC) $(HD_EXTENDED)
MHD_DIRS := $(MHD_BASIC) $(MHD_EXTENDED)
FFHD_DIRS := $(FFHD_BASIC) $(FFHD_EXTENDED)
SRHD_DIRS := $(SRHD_BASIC) $(SRHD_EXTENDED)
SRMHD_DIRS := $(SRMHD_BASIC) $(SRMHD_EXTENDED)


TEST_DIRS := $(HD_DIRS) $(MHD_DIRS) $(FFHD_DIRS) $(SRHD_DIRS) $(SRMHD_DIRS)

# $(call suite,BASIC,EXTENDED) gives the directories selected by SUITE
ifeq ($(SUITE),basic)
suite = $(1)
else ifeq ($(SUITE),extended)
suite = $(2)
else ifeq ($(SUITE),all)
suite = $(1) $(2)
else
$(error SUITE must be basic, extended or all, not '$(SUITE)')
endif

HD_SEL := $(call suite,$(HD_BASIC),$(HD_EXTENDED))
MHD_SEL := $(call suite,$(MHD_BASIC),$(MHD_EXTENDED))
FFHD_SEL := $(call suite,$(FFHD_BASIC),$(FFHD_EXTENDED))
SRHD_SEL := $(call suite,$(SRHD_BASIC),$(SRHD_EXTENDED))
SRMHD_SEL := $(call suite,$(SRMHD_BASIC),$(SRMHD_EXTENDED))

# Directories with altered names (useful for cleaning)
CLEAN_TEST := $(TEST_DIRS:%=clean-%)

.PHONY: all hd mhd ffhd srhd srmhd list clean $(TEST_DIRS) $(CLEAN_TEST)

all: hd mhd ffhd srhd srmhd
hd: $(HD_SEL)
mhd: $(MHD_SEL)
ffhd: $(FFHD_SEL)
srhd: $(SRHD_SEL)
srmhd: $(SRMHD_SEL)

list:
	@echo "SUITE=$(SUITE)"
	@for d in $(HD_SEL) $(MHD_SEL) $(FFHD_SEL) $(SRHD_SEL) $(SRMHD_SEL); do echo "  $$d"; done

$(TEST_DIRS):
	@$(MAKE) -C $@ -f test.make

clean:	$(CLEAN_TEST)

$(CLEAN_TEST):
	$(MAKE) -C $(@:clean-%=%) -f test.make clean
