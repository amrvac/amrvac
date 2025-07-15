# When make is called with command-line overrides, we want to store
# those for provenance. Also those overrides become defaults for
# recurrent invocations to make. Currently this is done for
# DEBUG and OPENMP flags. In general, any flag that doesn't modify
# the output of AMRVAC (so no physics, or anything meaningful),
# but change the efficiency or amount of analysis done should go here:
# OPENACC, PROFILE etc. How these flags change compilation under
# any specific architecture should be defined in the architecture file.

all: $(build_dir)/flags.mk $(build_dir)/config.mk

$(build_dir)/flags.mk:
	@mkdir -p $(@D)
	@echo "" > $@
ifdef DEBUG
	@echo "DEBUG ?= $(DEBUG)" >> $@
endif
ifdef OPENMP
	@echo "OPENMP ?= $(OPENMP)" >> $@
endif
ifdef OPENACC
	@echo "OPENACC ?= $(OPENACC)" >> $@
endif
ifdef NOGPUDIRECT
	@echo "NOGPUDIRECT ?= $(NOGPUDIRECT)" >> $@
endif
ifdef NVTX
	@echo "NVTX ?= $(NVTX)" >> $@
endif

$(build_dir)/config.mk: config.mk
	@mkdir -p $(@D)
	@cp $< $@
