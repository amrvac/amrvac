# Everything that DOES have meaning in terms of model output,
# should be configured in `amrvac.par`. The part  of this file
# that is relevant at  compile time is read by a Python
# script and converted into Makefile.

# TODO: maybe enabling a component only touches a small part
# of  the  code, so we don't need to recompile everything.
# To fix this we would need something of a tree based system
# where we can fetch module files from lower parts of a tree
# if they're unaffected. This implies a lot of added complexity
# to the build system though.

# TODO: Provenance. Fypp could insert a print statement into the
# Fortran code so that we can obtain git id and compile flags from
# simulation logs.

config.mk: amrvac.par
	@echo -e "Generating $(_magenta)amrvac.par$(_reset) -> $(_blue)config.mk$(_reset)"
	@python $(amrvac)/make/config_reader.py < $< > $@

include config.mk

