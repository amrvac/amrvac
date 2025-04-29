$(info Generating dependcy files)

# dep_files := $(f90_files:$(build_dir)/f90/%.f90=$(build_dir)/dep/%.dep)

# TODO: is there alternative to  fortdepend without dependencies?
# build check that  fortdepend is installed and inform user.

$(build_dir)/dependencies.mk: $(f90_files) | $(build_dir)
> @echo "Regenerating depencies"
> @fortdepend -f $(f90_files) -i mpi -b $(build_dir)/obj -w -o $@

include $(build_dir)/dependencies.mk

