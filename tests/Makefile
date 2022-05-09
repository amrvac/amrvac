# Makefile for automatically testing MPI-AMRVAC
#
# 1. You can run the rho tests with: make rho The same applies to other physics
# modules
#
# 2. You can run them in parallel with for example: -j 4
#
# 3. You can hide the Entering/leaving directory messages with: -s

RHO_DIRS := rho/auto_1d rho/auto_2d rho/auto_3d

HD_DIRS := hd/Riemann_1D hd/Riemann_2D hd/Rayleigh_Taylor_2D		\
hd/Rayleigh_Taylor_3D hd/Kelvin_Helmholtz_2D hd/blast_wave_Cartesian_2D	\
hd/blast_wave_polar_2D hd/blast_wave_polar_stretched_2D			\
hd/blast_wave_spherical_3D hd/Richtmyer_Meshkov_dust_2D hd/thermal_conduction_2D \
hd/blast_wave_spherical_stretched_3D hd/blast_wave_Cartesian_stretched_2D \
hd/blast_wave_cylindrical_2D hd/blast_wave_cylindrical_3D hd/blast_wave_polar_3D

MHD_DIRS := mhd/Riemann_1.75D mhd/rotor_2D mhd/wake_2.5D mhd/shock_cloud_2D	\
mhd/Kelvin_Helmholtz_2D mhd/Orszag_Tang_2D mhd/blast_wave_Cartesian_2D		\
mhd/blast_wave_polar_2D mhd/blast_wave_spherical_3D mhd/ripple_2.5D		\
mhd/solar_atmosphere_2.5D  mhd/solar_atmosphere_3D  mhd/blast_wave_Cartesian_3D \
mhd/radiative_cooling_3D mhd/blast_wave_polar_stretched_2D			\
mhd/tilt_instability_2D mhd/doubleGEM_2.5D mhd/Longcope_Strauss_2D		\
mhd/blast_wave_cylindrical_3D mhd/blast_wave_spherical_stretched_3D		\
mhd/convection_2D mhd/Kelvin_Helmholtz_3D mhd/ring_thermal_conduction_2D	\
mhd/ring_thermal_conduction_3D mhd/solar_flare_reconnection_2.5D                \
mhd/lfff_extrapolation_3D mhd/tilt_instability_3D mhd/Low_flux_rope_3D          \
mhd/field_loop_2.5D mhd/solar_wind_2.5D mhd/Orszag_Tang_2.5D mhd/Ambipolar1DWave mhd/Ambipolar1DWave-B0split

MG_DIRS := multigrid/convergence_2d multigrid/convergence_3d	\
multigrid/field_loop_2d multigrid/implicit_diffusion

RD_DIRS := rd/schnakenberg_2d rd/test_suite_2d

MF_DIRS := mf/Mackay_bipole_Cartesian_3D

TWOFL_DIRS := twofl/Alfven_wave_1D twofl/Riemann_1D twofl/KHI_2D

TEST_DIRS := $(RHO_DIRS) $(HD_DIRS) $(MHD_DIRS) $(MG_DIRS) $(RD_DIRS) $(MF_DIRS) $(TWOFL_DIRS)

# Directories with altered names (useful for cleaning)
CLEAN_TEST := $(TEST_DIRS:%=clean-%)

.PHONY: all rho hd mhd mg rd mf twofl clean $(TEST_DIRS) $(CLEAN_TEST)

all: rho hd mhd mg rd mf twofl
rho: $(RHO_DIRS)
hd: $(HD_DIRS)
mhd: $(MHD_DIRS)
mg: $(MG_DIRS)
rd: $(RD_DIRS)
mf: $(MF_DIRS)
twofl: $(TWOFL_DIRS)

$(TEST_DIRS):
	@$(MAKE) -C $@ -f test.make

clean:	$(CLEAN_TEST)

$(CLEAN_TEST):
	$(MAKE) -C $(@:clean-%=%) -f test.make clean
