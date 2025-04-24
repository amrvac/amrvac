source_files := $(shell find $(amrvac)/migration/src -name '*.fpp') mod_usr.f

$(build_dir)/dependencies.mk:
> fortdepend -f $(source_files) -i mpi -b $(build_dir)/obj -w -o $(build_dir)/dependencies.mk

-include $(build_dir)/dependencies.mk

