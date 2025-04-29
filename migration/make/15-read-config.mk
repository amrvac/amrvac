# Compile time options are read from `./config.mk`
ifneq ($(wildcard config.mk),)
$(info Reading config.mk)
include config.mk
else
$(info Warning: no config.mk found)
endif

# TODO add rule to write config.mk to $(build_dir)/config.mk
# also have a extra makefile with command line variables stored

