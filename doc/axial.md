# Using axial and translational symmetry

A number of different grid geometries are *in principle* supported in AMRVAC:

* Cartesian (or 'slab'): x, y, z
* Cylindrical: r, z, phi
* Polar: r, phi, z
* Spherical: r, phi, theta

For each grid type, there are 1D, 2D, and 3D versions. For some, there are also
special versions:

* 1.5D: 1D grid, but at each grid cell vectors have 2 components
* 2.5D: 2D grid, but at each grid cell vectors have 3 components

A user can specify the coordinate system by calling `set_coordinate_system` in his/her `usr_init()` routine. Here are some examples:

    call set_coordinate_system("Cartesian_2D")
    call set_coordinate_system("Cartesian_2.5D")
    call set_coordinate_system("Cartesian_3D")
    call set_coordinate_system("cylindrical_2D")
    call set_coordinate_system("polar_2D")
    call set_coordinate_system("spherical_2.5D")

The code should print an error when you select an option that is unavailable.

# Geometric source terms

In slab symmetry all fluxes entering the cells from the ignored direction cancel
exactly. In cylindrical or spherical symmetry, however, there are non- vanishing
contributions. These extra fluxes are always functions of the quantities in the
local cell and of the geometry only, since due to the symmetry assumption, the
quantities in the virtual neighbouring cells in the ignored direction(s) are the
same as in the local cell except for the geometrical transformation, e.g.
rotation of vector quantities around the symmetry axis. Consequently the
contributions are added like source terms. To maintain equilibrium better, the
'geometrical sources' are added at the same time as fluxes in the radial
direction.

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

