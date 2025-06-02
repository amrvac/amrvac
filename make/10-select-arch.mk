# if there was a previous build, default to the same arch,
# or use the default `gnu`

# This line reads the flags from latest build, but this turns
# out to be confusing
# -include $(build)/latest/flags.mk

# This selects the architecture from latest build if not given
ifndef arch
-include $(amrvac)/build/latest/arch.mk
arch ?= gnu
endif

# validate architecture selection
ifneq ("$(wildcard $(amrvac)/arch/$(arch).mk)","")
include $(amrvac)/arch/$(arch).mk
else
$(error Unknown architecture: $(arch))
endif

$(info Selected architecture: $(arch))

compile_flags ?= -c $(f90_flags)
link ?= $(compile)
link_flags ?= $(f90_flags)

