# Test 2D wind bubble

\test 2D wind bubble

# Description

This is a 2D hydro evolution, but in an axisymmetric setup, where there is
a user-defined source term that acts to set the stellar wind injection. We can activate 
(local) radiative loss effects. The setup models a stellar wind bubble, as it moves through the
ISM (the simulation is performed in the frame fixed on the star, the ISM blows
past the star from above).

We have several means to handle the energy 
loss term, and many of the frequently used cooling tables pre-coded.

Check out how we enforce AMR in the wind bubble region.

