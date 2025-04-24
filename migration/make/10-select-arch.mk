# if there was a previous build, default to the same arch,
# or use the default `gnu`
ifndef arch
-include build/latest/arch.mk
arch ?= gnu
endif

# validate architecture selection
ifneq ("$(wildcard arch/$(arch).mk)","")
include $(amrvac)/arch/$(arch).mk
else
$(error Unknown architecture: $(arch))
endif

$(info Selected architecture: $(arch))

