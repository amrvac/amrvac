# TODO: check AMRVAC_DIR
ARCH ?= default
BUILD_DIRS := build_1d build_2d build_3d
BUILD_MAKEFILE := $(AMRVAC_DIR)/arch/build.make

all .PHONY: $(BUILD_DIRS)

build_1d: NDIM=1
build_2d: NDIM=2
build_3d: NDIM=3

export ARCH NDIM

$(BUILD_DIRS):
	@mkdir -p $@
	$(MAKE) -C $@ -f $(BUILD_MAKEFILE)

clean .PHONY:
	$(RM) -r $(BUILD_DIRS)

