ifdef OPENMP
fypp_flags += -DOPENACC
endif

ifdef DEBUG
fypp_flags += -DDEBUG
endif

fypp_flags += -nn

$(info Fypp flags: $(fypp_flags))

.PRECIOUS: $(build_dir)/f90/%.f90
$(build_dir)/f90/%.f90: src/%.fypp
> @mkdir -p $(@D)
> @echo -e "Preprocessing $(_green)$<$(_reset)"
> @fypp $(fypp_flags) $< $@

