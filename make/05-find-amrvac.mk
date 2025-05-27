$(info AMRVAC Path: $(AMRVAC_DIR))

amrvac := $(AMRVAC_DIR)
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

