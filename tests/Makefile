TEST_DIRS := rho/auto_1d rho/auto_2d rho/auto_3d hd/Riemann_2D		\
hd/Rayleigh_Taylor hd/Kelvin_Helmholtz hd/jet_cloud mhd/doubleGEM_25D	\
mhd/solar_atmosphere_25D mhd/promRT

# Directories with altered names (useful for cleaning)
CLEAN_TEST := $(TEST_DIRS:%=clean-%)

.PHONY: all clean $(TEST_DIRS) $(CLEAN_TEST)

all: $(TEST_DIRS)

$(TEST_DIRS):
	@$(MAKE) -C $@ -f test.make

clean:	$(CLEAN_TEST)

$(CLEAN_TEST):
	$(MAKE) -C $(@:clean-%=%) clean
