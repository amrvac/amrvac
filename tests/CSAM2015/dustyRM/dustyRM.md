# Test dustyRM

\test dustyRM
\todo The test crashed (using gfortran compiler):

    too small pressure =  -0.11018446360550127       with limit=   0.0000000000000000 ...
    ERROR for processor           0 :
    === primitive pressure problem===

# Setup instructions

Setup this test case with no dust species:

    $AMRVAC_DIR/setup.pl -d=22 -phi=0 -z=0 -g=14,14 -p=hd -eos=default -nf=0 -ndust=0 -u=nul -arch=default

Or with four dust species:

    $AMRVAC_DIR/setup.pl -d=22 -phi=0 -z=0 -g=14,14 -p=hd -eos=default -nf=0 -ndust=4 -u=nul -arch=default

# Description

This is a 2D hydro setup, meant to demonstrate the use of gas-dust evolutions. The setup can be done without any dust species (`parfile_nodust`), or with an arbitrary number of dust bins respresenting a dust grain size distribution. An example for 4 dust species (`parfile_4dust`) is given.

The setup studies a Mach 2 shock, impinging on an inclined density discontinuity (density contrast of 3). Dust is only present beyond the density discontinuity. The pure hydro setup will show you double mach reflection (twice, both bottom and top are handled reflective), along with Richtmyer-Meshkov (i.e. Kelvin-Helmholtz on the shocked contact) instability development. If dust is included, the dust grain size will determine in how far the dust is coupled to the gas, and possibily evacuated from the vortices.

The `amrvacusr.t` file shows you how to add additional output variables to the vtu files which are generated during runtime. In this example, the vtu files contain the primitive variables, together with a Schlieren plot of the density (i.e. the stretched density gradient).

Check out how we run at Courant parameter of 1.5 (thanks to the `ssprk54` time stepping).

Note: the `normvar(*)` array is used only in the conversion to vtu files, and here also in the drag force computation when dust is present. Note that Paraview may think that the ranges for (dust) densities are too small to render, but this is obviously not true (use the calculator to multiply to a range that paraview does recognize...).

Related research is in *Effect of dust on Kelvin-Helmholtz instabilities*, T.
Hendrix \& R. Keppens, 2014, Astron. \& Astrophys. 526, A114,
`doi:10.1051/0004-6361/201322322`


