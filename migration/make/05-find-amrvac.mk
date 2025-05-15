$(info AMRVAC Path: $(AMRVAC_DIR))

#FIXME remove migration subdir once done
amrvac := $(AMRVAC_DIR)/migration
build := $(amrvac)/build

ifeq ($(findstring help, $(MAKECMDGOALS)), help)
disable_precompile := 1
endif

ifeq ($(findstring clean, $(MAKECMDGOALS)), clean)
disable_precompile := 1
endif

.PHONY: help

help:
> less $(amrvac)/make/README.md

