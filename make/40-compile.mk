# obj_files = $(patsubst %.f90, $(build_dir)/obj/%, $(notdir $(f90_files)))
# TODO: get rid of amrvac.h? YES!

compile_flags += -I$(build_dir)/f90

$(build_dir)/f90/amrvac.h:
> @touch $@

# The GNU Fortran compiler has a -J flag to place the .mod files,
# but this is not supported by other compilers. Running the compile
# command in the f90 directory should fix this.

$(build_dir)/obj/%.o: $(build_dir)/f90/%.f90 $(build_dir)/f90/amrvac.h
> @mkdir -p $(@D)
> @echo -e "Compiling $(_magenta)$(notdir $<)$(_reset)"
> @cd $(build_dir)/f90; $(compile) $(compile_flags) $< -o $@

