# Compile time options are read from `./config.mk`
ifneq ($(wildcard config.mk),)
$(info Reading config.mk)
include config.mk
else
$(info Warning: no config.mk found)
endif

