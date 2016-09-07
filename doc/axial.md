# Using axial and translational symmetry

This document describes how to use properly axial symmetry and polar grids. It
is relevant when doing simulations with **typeaxial** different from 'slab' in
the **methodlist** of your **par/PROBLEM** parameter file.

## What does Axial Symmetry Mean?

In 1 and 2D simulations, there is an inherent assumption about the behavior of
quantities in the ignored direction(s). The most usual assumption is 'slab'
symmetry, i.e. all the quantities are invariant to a translation in the
ignored direction.

In 2D, one may assume a rotational invariance around some axis instead, which
gives a cylindrical symmetry, thus **typeaxial='cylindrical'** should be
selected. MPI-AMRVAC always assumes that the radial distance from the symmetry
axis is in the first coordinate, which can therefore be referred to as **r_**,
where **r_=1** is an integer parameter set in **src/amrvacdef.t**. When
setting

    setup.pl -d=23 -z=2 -phi=3

the second coordinate is interpreted as the one parallel to the axis,
**z_=2**, while the ignored direction is **phi_=3** in 2.5D. The PHI
components of vector variables are then stored in the third component in 2.5D
simulations.

Actually, for 2D and 'cylindrical' (which can be 2D or 2.5D) the grid and the
symmetry depend on the settings for the -phi and -z flags. When -d=22 -z=2, a
cartesian grid is used in a poloidal plane, with axial symmetry for the r- and
z- vector components. The same is true for -d=23 -z=2 -phi=3, when all three
vector components are then depending on (r,z) coordinates only. The vector
components then appear in the r,z,phi order. One can use 2.5D on a circular
grid also with translational symmetry in the third, i.e. axial, direction by
the use of _setup.pl -d=23 -phi=2 -z=3_. The vector components then appear in
the r,phi,z order and they are depending on (r,phi) coordinates.

In spherical symmetry, **typeaxial='spherical'**, the invariant transformation
is a rotation around the origin. Again, the first coordinate is **r_=1**. For
2D and 'spherical', the coordinates denote radial and polar angle in a
poloidal plane, and similarly for 2.5D in combination with spherical where the
third vector component is then the phi-component, depending on (r,theta)
alone. This means that 2D and 'spherical' implies the use of a 2D spherical 
grid in a poloidal cross-section (i.e. one containing the symmetry axis) of a 
3D sphercial coordinate and axial symmetry. Note that the _-phi -z_ flags  
have no influence on the 2D spherical grid.

## What does Polar Grid Mean?

In case of a polar grid (**typeaxial='cylindrical'**) the PHI direction is not
ignored.

In 2D, the data in the **x** coordinates is interpreted as radial distance
**r** and the angular coordinate **phi** is measured in radians. The **dx**
array contains **dr** and **r*dphi** and MPI-AMRVAC treats everything the same
way as in a Cartesian grid. Meaningfull 2D settings are

    setup.pl -d=22 -phi=2 -z=0
    setup.pl -d=23 -phi=2 -z=3

The main difference is that the geometric source terms are added due to
**typeaxial='cylindrical'**. The former setting is purely 2D, the latter
setting is 2.5D, and assumes translational invariance with the z direction.

In 3D, there are two possibilities: the coordinates can be ordered as **r, z,
phi** or **r, phi, z**, where **z** is the vertical position.

    setup.pl -d=33 -phi=2 -z=3
    setup.pl -d=33 -phi=3 -z=2

The vector variables also have their components in the **r, phi**, and **z**
directions. These facts should be kept in mind when the initial conditions are
defined, when source terms are applied, or when the results are visualized.

## Implementation

In slab symmetry all fluxes entering the cells from the ignored direction
cancel exactly. In cylindrical or spherical symmetry, however, there are non-
vanishing contributions. These extra fluxes are always functions of the
quantities in the local cell and of the geometry only, since due to the
symmetry assumption, the quantities in the virtual neighbouring cells in the
ignored direction(s) are the same as in the local cell except for the
geometrical transformation, e.g. rotation of vector quantities around the
symmetry axis. Consequently the contributions are added like source terms in
the **addgeometry** subroutine in **src/amrvacphys**. To maintain equilibrium
better, the 'geometrical sources' are added at the same time as fluxes in the
radial direction.

On a polar grid the extra terms are exactly the same as for cylindrical
symmetry.

## Boundary at the Symmetry Axis

The symmetry axis is actually a special case of periodicity. For 3D
cylindrical and spherical grid computations, the singular polar axis is
trivially handled using a so-called pi-periodic boundary treatment, where
periodicity across the pole comes from the grid cell diagonally across the
pole, i.e. displaced over pi instead of 2 pi. These are automatically
recognized from the typeaxial setting, if the radial range starts at zero. The
corresponding range in angle phi must span 2 pi for cylindrical, and theta
must then start at zero (to include the north pole) and/or end at pi (for the
south pole) for spherical grids. The user just needs to set the typeB as if
the singular axis is a symmetry boundary (using symm and asymm combinations).
**There is one important restriction: the number of grid blocks must be even
around the pole.**

