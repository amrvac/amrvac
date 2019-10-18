#!/usr/bin/env bash

cd $AMRVAC_DIR/src
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/src/amrvacio
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/src/physics
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/src/rho
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/src/hd
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/src/mhd
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/src/nonlinear
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/src/particle
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/src/modules
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/src/includes
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/schemes
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/rho/auto_1d
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/rho/auto_2d
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/rho/auto_3d
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/rho/convergence_new
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/hd/Kelvin_Helmholtz_2D
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/hd/Riemann_1D
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/hd/Riemann_2D
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/hd/Rayleigh_Taylor_2D
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/hd/Rayleigh_Taylor_3D
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/hd/Rayleigh_Taylor_particles_advect_2D
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/hd/Richtmyer_Meshkov
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/hd/Richtmyer_Meshkov_dust_2D
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/hd/Riemann_1D
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/hd/Riemann_2D
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/hd/Woodward_Collela_2D
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/hd/blast_wave_Cartesian_2D
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/hd/blast_wave_Cartesian_3D
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/hd/blast_wave_Cartesian_stretched_2D
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/hd/blast_wave_cylindrical_2D
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/hd/blast_wave_cylindrical_3D
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/hd/blast_wave_polar_2D
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/hd/blast_wave_polar_3D
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/hd/blast_wave_polar_stretched_2D
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/hd/blast_wave_spherical_3D
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/hd/blast_wave_spherical_stretched_3D
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/hd/radiative_cooling_3D
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/hd/thermal_conduction_3D
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/hd/jet_cloud
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/mhd/Alfven_wave_2.5D
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/mhd/GEM_2D
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/mhd/Kelvin_Helmholtz_2D
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/mhd/Kelvin_Helmholtz_3D
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/mhd/Kelvin_Helmholtz_double_2D
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/mhd/Longcope_Strauss_2D
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/mhd/Orszag_Tang_2D
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/mhd/Orszag_Tang_particles_advect_2D
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/mhd/Orszag_Tang_particles_gca_2.5D
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/mhd/Riemann_1.75D
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/mhd/blast_wave_Cartesian_2D
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/mhd/blast_wave_Cartesian_3D
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/mhd/blast_wave_cylindrical_3D
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/mhd/blast_wave_polar_2D
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/mhd/blast_wave_polar_stretched_2D
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/mhd/blast_wave_spherical_2D
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/mhd/blast_wave_spherical_3D
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/mhd/blast_wave_spherical_stretched_3D
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/mhd/convection_2D
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/mhd/doubleGEM_2.5D
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/mhd/Low_flux_rope_3D
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/mhd/prominence_Rayleigh_Taylor_2.5D
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/mhd/radiative_cooling_3D
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/mhd/ripple_2.5D
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/mhd/rotor_2D
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/mhd/shock_cloud_2D
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/mhd/solar_atmosphere_2.5D
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/mhd/thermal_conduction_3D
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/mhd/tilt_instability_2D
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/mhd/tilt_instability_3D
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/mhd/wake_2.5D
$AMRVAC_DIR/tools/upgrade.pl
