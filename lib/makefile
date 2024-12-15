# TODO: check AMRVAC_DIR
ARCH ?= default
DIRS := 1d1_$(ARCH) 1d2_$(ARCH) 1d3_$(ARCH) 2d2_$(ARCH) 2d3_$(ARCH) 3d3_$(ARCH)
LIB_MAKE := $(AMRVAC_DIR)/arch/lib.make

.PHONY: all clean 1d1 1d2 1d3 2d2 2d3 3d3 $(DIRS)

all: 1d1 1d2 1d3 2d2 2d3 3d3

# To allow user to type e.g. make 1d1
1d1: 1d1_$(ARCH)
1d2: 1d2_$(ARCH)
1d3: 1d3_$(ARCH)
2d2: 2d2_$(ARCH)
2d3: 2d3_$(ARCH)
3d3: 3d3_$(ARCH)

# NDIM will be used during the build
1d1_$(ARCH): NDIM=1
1d2_$(ARCH): NDIM=1
1d3_$(ARCH): NDIM=1
2d2_$(ARCH): NDIM=2
2d3_$(ARCH): NDIM=2
3d3_$(ARCH): NDIM=3
# NDIR will be used during the build
1d1_$(ARCH): NDIR=1
1d2_$(ARCH): NDIR=2
1d3_$(ARCH): NDIR=3
2d2_$(ARCH): NDIR=2
2d3_$(ARCH): NDIR=3
3d3_$(ARCH): NDIR=3

export ARCH NDIM NDIR

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
	$(RM) -r 1d1_* 1d2_* 1d3_* 2d2_* 2d3_* 3d3_* 1d_* 2d_* 3d_*

