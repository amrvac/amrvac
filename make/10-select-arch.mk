# if there was a previous build, default to the same arch,
# or use the default `gnu`
-include $(build)/latest/flags.mk

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

