
$(build_dir)/obj/%.o: $(build_dir)/f90/%.f90
> @mkdir -p $(@D)
> @echo -e "Compiling \\e[1;32m$<\\e[m"
> @$(compile) $(compile_flags) $< -o $@

