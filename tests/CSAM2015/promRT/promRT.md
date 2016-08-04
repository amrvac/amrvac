# Test promRT

\test promRT

# Setup instructions

Setup this test case in 2D:

    $AMRVAC_DIR/setup.pl -d=23 -phi=0 -z=0 -g=14,14 -p=mhd -eos=default -nf=1 -ndust=0 -u=nul -arch=default

Or in 3D:

    $AMRVAC_DIR/setup.pl -d=33 -phi=0 -z=0 -g=14,14,14 -p=mhd -eos=default -nf=1 -ndust=0 -u=nul -arch=default

# Description

This directory allows you to explore (recent) research into Rayleigh-Taylor
dynamics in macroscopic prominence segments, using the ideal MHD assumption. In
a stratified atmosphere, connecting chromosphere to corona, we have a prominence
which is unstable to Rayleigh-Taylor instabilities. The `amrvacusr.t` contains
many aspects for you to explore, namely: how we setup a stratified atmosphere
for a given temperature distribution, and how we use the initial stratification
also in the boundary conditions (top and bottom boundaries are treated special
here, where the $y$-direction has gravity). The code is ok for both 2.5D or 3D
runs, but obviously the 3D run, for simulations over long timescale, will need
to be done overnight.

Check out how we activate a tracer to identify different parts (chromosphere,
corona, prominence), and how in post-processing, we can use the user convert
option to collect integrals over domain parts, to quantify information like the
average temperature in some region, the volume and mass, etc.

The initial condition uses velocity pertubations with random phases, check out
how this is to be done when running the parallel, block-adaptive code.

The boundary conditions ate top and bottom are instructive to see how
extrapolations can be easily done on primitive variables. We have external
gravity as a main physical ingedient here. Also, the log file is coded up very
generally to store integrals over the domain.

The 3D setup also shows you how to generate on-the-fly collapsed views, along
the three coordinate axes.

Related work is described in: *Solar prominences: "double, double ... boil and
bubble"*, R. Keppens, X. Cia, \& O. Porth, 2015, ApJ Letters 806, L13 (7pp).
`doi:10.1088/2041-8205/806/1/L13`


