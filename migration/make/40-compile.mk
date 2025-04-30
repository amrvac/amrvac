# obj_files = $(patsubst %.f90, $(build_dir)/obj/%, $(notdir $(f90_files)))
# TODO: get rid of amrvac.h?
# TODO: set location of mod files, this location add to -I

compile_flags += -I$(build_dir)/f90

$(build_dir)/f90/amrvac.h:
> @touch $@

$(build_dir)/obj/%.o: $(build_dir)/f90/%.f90 $(build_dir)/f90/amrvac.h
> @mkdir -p $(@D)
> @echo -e "Compiling $(_magenta)$(notdir $<)$(_reset)"
> @$(compile) $(compile_flags) -J$(build_dir)/f90 $< -o $@

