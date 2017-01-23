# TODO: check AMRVAC_DIR
ARCH ?= default
BUILD_DIRS := 1d_$(ARCH) 2d_$(ARCH) 3d_$(ARCH)
BUILD_MAKEFILE := $(AMRVAC_DIR)/arch/build.make

.PHONY: all clean 1d 2d 3d $(BUILD_DIRS)

all: 1d 2d 3d

# To allow user to type e.g. make 1d
1d: 1d_$(ARCH)
2d: 2d_$(ARCH)
3d: 3d_$(ARCH)

# NDIM will be used during the build
1d_$(ARCH): NDIM=1
2d_$(ARCH): NDIM=2
3d_$(ARCH): NDIM=3

export ARCH NDIM

# How to make the builds
$(BUILD_DIRS):
	@mkdir -p $@
	$(MAKE) -C $@ -f $(BUILD_MAKEFILE)

# Remove everything
clean:
	$(RM) -r 1d_* 2d_* 3d_*

