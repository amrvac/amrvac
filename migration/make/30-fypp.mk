ifdef OPENMP
fypp_flags += -DOPENACC
endif

ifdef DEBUG
fypp_flags += -DDEBUG
endif

fypp_flags += -n

$(info Fypp flags: $(fypp_flags))

source_files := $(shell find $(amrvac)/src -name '*.fpp') mod_usr.fpp
f90_files := $(patsubst %.fpp, $(build_dir)/f90/%.f90, \
	     $(patsubst $(amrvac)/src/%.fpp, $(build_dir)/f90/%.f90, \
	     $(source_files)))

.PRECIOUS: $(f90_files)

$(build_dir)/f90/%.f90: $(amrvac)/src/%.fpp
> @mkdir -p $(@D)
> @echo -e "Preprocessing $(_magenta)$(^:$(amrvac)/%=%)$(_reset)"
> @fypp $(fypp_flags) $^ $@

$(build_dir)/f90/mod_usr.f90: mod_usr.fpp | $(build_dir)
> @mkdir -p $(@D)
> @echo -e "Preprocessing $(_magenta)$^$(_reset)"
> @fypp $(fypp_flags) $^ $@

