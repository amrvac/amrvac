# Test pfss

\test pfss
\todo Does not compile:

    eigensystem.f:430.68:
    call gradient(w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,p_),ixImin1,&
                                                                    1
    Error: Symbol 'p_' at (1) has no IMPLICIT type

# Setup instructions

Setup this test case with

    $AMRVAC_DIR/setup.pl -d=33 -phi=3 -z=0 -g=14,14,14 -p=mhd -eos=iso

# Description

This is a very similar exercise as the one above, but this time generating a
Potential Field Source Surface model from a global magnetogram. It also has an
IDL routine `convertsynoptic.pro` provided for converting downloadable synoptic
maps to more easily handled maps of the radial field component, on a spherical
grid, uniform in both spherical angles. This file is read in (and an example
*dat for a specific date is given) by the user file (again with isothermal MHD).

Check the user file to see how we use the PFSS module (see the
`INCLUDE::amrvacmodules/pfss.t` statement). The par file is set to stop after
the field is generated. You can visualize again with Paraview. Note how we can
use AMR to zoom in on specific active regions.

A description of this functionality and how the PFSS extrapolations compares to
local cartesian extrapolations, can be found in *MPI-AMRVAC for solar and
astrophysics*, O. Porth, C. Xia, T. Hendrix, S.P. Moschou, \& R. Keppens, 2014,
ApJS 214, 4 (26pp) `doi:10.1088/0067-0049/214/1/4`.


