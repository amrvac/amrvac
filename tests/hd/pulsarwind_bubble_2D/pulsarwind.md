# Test 2D pulsar wind bubble

\test 2D pulsar wind bubble

# Description

This is a 2D hydro evolution, but in an axisymmetric setup, where there is
a user-defined source term that acts to set a pulsar wind injection. This reproduces the
setting in Swaluw et al., A&A 397, 913 (2003) doi:10.1051/0004-6361:20021488.
The setup reproduces conditions for a pulsar which achieved a kick velocity large enough to meet up with the
shocked shell of the supernova (that led to the pulsar) in the late sedov phase. The simulation just models the
pulsar wind bubble, as it moves with a mach number M_psr=3.13, as predicted in Swaluw et al, 2003.
The setup realizes a steady state bow shock, and illustrates the usage of stretched (AMR) grids.


Check out how we enforce AMR in the pulsar wind injection region. A fix_small_values may be needed for numerical stability.

