#!/usr/bin/env bash

cd $AMRVAC_DIR/src
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/src/amrvacio
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/src/physics
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/src/modules
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/src/rho
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/src/hd
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/src/mhd
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/src/nonlinear
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/schemes
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/rho/auto_1d
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/rho/auto_2d
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/rho/auto_3d
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/rho/auto_3d
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/hd/Riemann_2D
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/hd/Rayleigh_Taylor_2D
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/hd/Rayleigh_Taylor_3D
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/hd/Kelvin_Helmholtz
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/hd/jet_cloud
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/mhd/ripple_2.5D
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/mhd/Kelvin_Helmholtz_2D
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/mhd/radiative_cooling_3D
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/mhd/thermal_conduction_3D
$AMRVAC_DIR/tools/upgrade.pl

cd $AMRVAC_DIR/tests/mhd/solar_atmosphere_2.5D
$AMRVAC_DIR/tools/upgrade.pl
