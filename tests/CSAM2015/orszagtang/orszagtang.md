# Test orszagtang

\test orszagtang

# Setup instructions

Setup this test case in 2D:

    $AMRVAC_DIR/setup.pl -d=22 -phi=0 -z=0 -g=14,14 -p=mhd -eos=default -nf=0 -ndust=0 -u=nul -arch=default

Or in 3D:

    $AMRVAC_DIR/setup.pl -d=33 -phi=0 -z=0 -g=14,14,14 -p=mhd -eos=default -nf=0 -ndust=0 -u=nul -arch=default

# Description

This directory contains a setup for the Orszag-Tang problem, a typical MHD
testcase. The code is set up to handle both 2D and 3D variants, under either an
isothermal closure (no energy equation) or in fully adiabatic ideal MHD (with
energy evolution). This example serves to show the use of the `-eos=` flag
at compile time. The test lets a Mach 1 vortical flow distort a series of
magnetic islands.

Examples are given for 2D and 3D, for isothermal and full MHD
(`amrvac_otmhd2D.par`, `amrvac_otmhd3D.par`, `amrvac_otmhdiso2D.par`,
`amrvac_otmhdiso3D.par`). Use `diff` to spot obvious differences in these
parfiles.

Boundaries are (double to triple) periodic. Note the settings for resolution (up
to 5 grid levels in the 2D runs, domain decomposition at $100^3$ in the 3D
runs). When simulated far enough, and at high resolution, these ideal MHD 2D
runs will show nice tearing instabilities on the developing current sheets. Note
that there is no explicit resistivity here. The monopole condition is handled
through a diffusive approach (`typedivbfix='linde'`).


