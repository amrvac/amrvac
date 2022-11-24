TEST_DIRS := Riemann_1D Riemann_2D	\
Rayleigh_Taylor_2D Rayleigh_Taylor_3D Kelvin_Helmholtz_2D		\
blast_wave_Cartesian_2D blast_wave_polar_2D blast_wave_polar_stretched_2D \
blast_wave_spherical_3D Richtmyer_Meshkov_dust_2D thermal_conduction_2D \
blast_wave_spherical_stretched_3D blast_wave_Cartesian_stretched_2D \
blast_wave_cylindrical_2D blast_wave_cylindrical_3D blast_wave_polar_3D \
CAKwind_spherical_1D CAKwind_spherical_2.5D

# Directories with altered names (useful for cleaning)
CLEAN_TEST := $(TEST_DIRS:%=clean-%)

.PHONY: all clean $(TEST_DIRS) $(CLEAN_TEST)

all: $(TEST_DIRS)

$(TEST_DIRS):
	@$(MAKE) -C $@ -f test.make

clean:	$(CLEAN_TEST)

$(CLEAN_TEST):
	$(MAKE) -C $(@:clean-%=%) clean
