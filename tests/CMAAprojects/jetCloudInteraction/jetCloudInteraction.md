# Test jetCloudinteraction

\test jetCloudinteraction
\todo Merge with CSAM2015 version of this test

# Setup instructions

Setup this test case in 2D with

    $AMRVAC_DIR/setup.pl -d=22 -phi=0 -z=0 -g=16,16 -p=hd -eos=default -nf=1 -ndust=4 -u=nul -arch=default

Or in 3D with

    $AMRVAC_DIR/setup.pl -d=22 -phi=0 -z=0 -g=16,16 -p=hd -eos=default -nf=1 -ndust=4 -u=nul -arch=default

Then create an output directory:

    mkdir output

# Description

In this Euler (gas dynamical) test, we simulate the impact of a Mach 12 jet flow
as it passes through lighter external medium (density contrast 3), and finally
impacts on a spherical (density stratified) cloud which is denser than the jet.
It is inspired from a 3D SPH study done by De Gouveia del Pino, E., ApJ 526,
862-873 (1999).

The user file can be used for either a 2D (purely planar, Cartesian) or a 3D
setup. It is meant to demonstrate the use of having an additional tracer
quantity added to the Euler equations.

The amrvacusr.t file shows how to set up the initial condition, and how to
control the left ($x=0$) inlet boundary by fixing all (primitive) quantities.

The parfiles demonstrate the use of combining boundary prescriptions that are
already available (free-outflow with `cont' or `noinflow' conditions) with a
special (here, fixed) treatment of the inlet.

Use

    diff amrvac2D.par amrvac3D.par

to spot simple differences between a 2D and a 3D variant. Check for `mxnest` to
see how many grid levels you use, together with the (AMR level one) grid
settings in `nxlone1, nxlone2, nxlone3`.

You can try to experiment with different resolutions, different numerical
schemes (it now uses a strong stability preserving Runge-Kutta scheme, and a MP5
limiter in an HLLC discretization).


