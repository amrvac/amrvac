# Makefile for automatically testing MPI-AMRVAC
#
# 1. You can run the rho tests with: make rho The same applies to other physics
# modules
#
# 2. You can run them in parallel with for example: -j 4
#
# 3. You can hide the Entering/leaving directory messages with: -s

RHO_DIRS := rho/auto_1d rho/auto_2d rho/auto_3d

NL_DIRS := nonlinear/burgers2d nonlinear/nonconvex_1d

DEMO_DIRS := demo/Advect_BC_from_file_2D demo/Shu_Osher_1D_HD demo/AlfvenShock_MHD2D

HD_DIRS := hd/Riemann_1D hd/Riemann_2D hd/Rayleigh_Taylor_2D		\
hd/Rayleigh_Taylor_3D hd/Kelvin_Helmholtz_2D hd/blast_wave_Cartesian_2D	\
hd/blast_wave_polar_2D hd/blast_wave_polar_stretched_2D			\
hd/blast_wave_spherical_3D hd/Richtmyer_Meshkov_dust_2D hd/thermal_conduction_2D \
hd/blast_wave_spherical_stretched_3D hd/blast_wave_Cartesian_stretched_2D \
hd/blast_wave_cylindrical_2D hd/blast_wave_cylindrical_3D hd/blast_wave_polar_3D \
hd/CAKwind_spherical_1D hd/CAKwind_spherical_2.5D hd/thermal_instability_1D hd/thermal_instability_2D \
hd/Coleman_1D

MHD_DIRS := mhd/Riemann_1.75D mhd/rotor_2D mhd/wake_2.5D mhd/shock_cloud_2D mhd/Kelvin_Helmholtz_2D \
mhd/Orszag_Tang_2D mhd/blast_wave_Cartesian_2D mhd/blast_wave_polar_2D mhd/blast_wave_polar_stretched_2D \
mhd/blast_wave_spherical_3D mhd/ripple_2.5D	\
mhd/Alfven_wave_2.5D  \
mhd/solar_atmosphere_2.5D  mhd/radiative_cooling_3D \
mhd/tilt_instability_2D mhd/tilt_instability_3D mhd/doubleGEM_2.5D mhd/Longcope_Strauss_2D mhd/blast_wave_cylindrical_3D \
mhd/blast_wave_spherical_stretched_3D mhd/convection_2D mhd/Kelvin_Helmholtz_3D mhd/Kelvin_Helmholtz_double_2D \
mhd/ring_thermal_conduction_2D mhd/ring_thermal_conduction_3D mhd/solar_flare_reconnection_2.5D mhd/lfff_extrapolation_3D  \
mhd/Low_flux_rope_3D mhd/field_loop_2.5D mhd/solar_wind_2.5D mhd/Orszag_Tang_2.5D mhd/solar_atmosphere_3D  mhd/alfvenwing_3D \
mhd/semirelativistic_Alfven_wave_1.75D \
mhd/blast_wave_Cartesian_3D mhd/Ambipolar1DWave-B0split mhd/Ambipolar1DWave mhd/icarus

RHD_DIRS := rhd/RHD_heatcool rhd/RHD_pulse_1D rhd/RHD_shock_1D rhd/RHD_wave_1D

SRHD_DIRS := srhd/Jet_cylindrical_2D

RMHD_DIRS := rmhd/Shocks_1.75D/Jiang rmhd/Shocks_1.75D/Fast rmhd/Shocks_1.75D/Slow rmhd/Shocks_1.75D/Intermediate rmhd/Waves_fine_1.75D/Alfven rmhd/Waves_fine_1.75D/Fast
# rmhd/Brio_Wu_1.75D rmhd/Torrilhon_1.75D rmhd/Convection_2D

MG_DIRS := multigrid/convergence_2d multigrid/convergence_3d	\
multigrid/field_loop_2d multigrid/implicit_diffusion

RD_DIRS := rd/schnakenberg_2d rd/test_suite_2d

FFHD_DIRS := ffhd/solar_atmosphere_2D ffhd/ring_thermal_conduction_2D

ARD_DIRS := ard/schnakenberg_2d

MF_DIRS := mf/Mackay_bipole_Cartesian_3D

TWOFL_DIRS := twofl/Alfven_wave_1D twofl/Riemann_1D twofl/KHI_2D

TEST_DIRS := $(RHO_DIRS) $(NL_DIRS) $(DEMO_DIRS) $(HD_DIRS) $(SRHD_DIRS) $(RHD_DIRS) $(RMHD_DIRS) $(MHD_DIRS) $(MG_DIRS) $(RD_DIRS) $(ARD_DIRS)  $(MF_DIRS) $(TWOFL_DIRS) $(FFHD_DIRS)

# Directories with altered names (useful for cleaning)
CLEAN_TEST := $(TEST_DIRS:%=clean-%)

.PHONY: all rho nl demo hd srhd rhd rmhd mhd mg rd ard mf twofl ffhd clean $(TEST_DIRS) $(CLEAN_TEST)

all: rho nl demo hd srhd rhd rmhd mhd mg rd ard mf twofl ffhd
rho: $(RHO_DIRS)
nl: $(NL_DIRS)
demo: $(DEMO_DIRS)
hd: $(HD_DIRS)
srhd: $(SRHD_DIRS)
rhd: $(RHD_DIRS)
rmhd: $(RMHD_DIRS)
mhd: $(MHD_DIRS)
mg: $(MG_DIRS)
rd: $(RD_DIRS)
ard: $(ARD_DIRS)
mf: $(MF_DIRS)
twofl: $(TWOFL_DIRS)
ffhd: $(FFHD_DIRS)

$(TEST_DIRS):
	@$(MAKE) -C $@ -f test.make

clean:	$(CLEAN_TEST)

$(CLEAN_TEST):
	$(MAKE) -C $(@:clean-%=%) -f test.make clean
