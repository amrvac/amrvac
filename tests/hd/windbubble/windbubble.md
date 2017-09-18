# Test windbubble

\test windbubble
\todo There seems to be no (working) par file

# Setup instructions

Setup this test case with

    $AMRVAC_DIR/setup.pl -d=22 -phi=0 -z=2 -g=14,14 -p=hd -eos=default -nf=0 -ndust=0 -u=nul -arch=default

# Description

This is also a 2D hydro evolution, but in an axisymmetric setup, where there is
a user-defined source term, as well as (local) radiative loss effects
incorporated. The setup models a stellar wind bubble, as it moves through the
ISM (the simulation is performed in the frame fixed on the star, the ISM blows
past the star from above).

The `amrvacusr.t` and `amrvacusrpar.t` shows how to incorporate the radiative
loss module (see e.g. `INCLUDE:amrvacmodules/cooling.t` and other statements
accessing its functionality), where we have several means to handle the energy
loss term, and many of the frequently used cooling tables pre-coded.

Check out how we enforce AMR in the wind bubble region, using the
`specialrefine_grid` subroutine.


