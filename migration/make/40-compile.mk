# obj_files = $(patsubst %.f90, $(build_dir)/obj/%, $(notdir $(f90_files)))

define compile_rule
$(patsubst %.f90, $(build_dir)/obj/%.o, $(notdir $(1))): $(1)
> @mkdir -p $$(@D)
> @echo -e "Compiling \\e[1;32m$$<\\e[m"
> @$(compile) $(compile_flags) $$< -o $$@
endef

$(foreach f, $(f90_files), $(eval $(call compile_rule, $(f))))

