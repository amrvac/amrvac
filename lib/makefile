# TODO: check AMRVAC_DIR
ARCH ?= default
DIRS := 1d_$(ARCH) 2d_$(ARCH) 3d_$(ARCH)
LIB_MAKE := $(AMRVAC_DIR)/arch/lib.make

.PHONY: all clean 1d 2d 3d $(DIRS)

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
$(DIRS):
	@mkdir -p $@
ifeq ("$(wildcard amrvac.h)","")
	touch amrvac.h
endif
	rsync -c amrvac.h $@/amrvac.h
	$(MAKE) -C $@ -f $(LIB_MAKE)

# Remove everything
clean:
	$(RM) -r 1d_* 2d_* 3d_*

