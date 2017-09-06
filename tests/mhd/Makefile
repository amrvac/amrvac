TEST_DIRS := Riemann_1.75D rotor_2D wake_2.5D shock_cloud_2D Kelvin_Helmholtz_2D		\
Orszag_Tang_2D blast_wave_Cartesian_2D blast_wave_polar_2D blast_wave_spherical_2D \
blast_wave_spherical_3D ripple_2.5D prominence_Rayleigh_Taylor_2.5D	\
solar_atmosphere_2.5D radiative_cooling_3D thermal_conduction_3D	\
tilt_instability_2D doubleGEM_2.5D Longcope_Strauss_2D blast_wave_cylindrical_3D \
blast_wave_spherical_stretched_3D convection_2D Kelvin_Helmholtz_3D Kelvin_Helmholtz_double_2D

# Directories with altered names (useful for cleaning)
CLEAN_TEST := $(TEST_DIRS:%=clean-%)

.PHONY: all clean $(TEST_DIRS) $(CLEAN_TEST)

all: $(TEST_DIRS)

$(TEST_DIRS):
	@$(MAKE) -C $@ -f test.make

clean:	$(CLEAN_TEST)

$(CLEAN_TEST):
	$(MAKE) -C $(@:clean-%=%) clean
