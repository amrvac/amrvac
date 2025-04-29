# obj_files = $(patsubst %.f90, $(build_dir)/obj/%, $(notdir $(f90_files)))
# TODO: get rid of amrvac.h?
# TODO: set location of mod files, this location add to -I

compile_flags += -I$(build_dir)/f90

$(build_dir)/f90/amrvac.h:
> @touch $@

define compile_rule
$(patsubst %.f90, $(build_dir)/obj/%.o, $(notdir $(1))): $(1) $(build_dir)/f90/amrvac.h
> @mkdir -p $$(@D)
> @echo -e "Compiling $(_magenta)$$<$(_reset)"
> @$(compile) $(compile_flags) $$< -o $$@
endef

$(foreach f, $(f90_files), $(eval $(call compile_rule, $(f))))

