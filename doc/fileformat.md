# MPI-AMRVAC's data file format

All data files consist of a single snapshot, and they can be used for restart
and/or further conversion to other data formats directly usable for
visualization. Note that restart is possible on a differing number of CPUs,
and may suddenly allow more refinement levels. Also, note that the individual
snapshots will typically have different lengths, as the number of grid blocks
will vary dynamically. The data is saved in binary format (double precision).
You can find the exact implementation in _amrio.t_, in the subroutine
_saveamrfile_. In fact, there are three options to write out the data files,
differing in their means to use parallel I/O from MPI-2. They are selected
through **typeparIO** which is an integer taking values **1,0,-1**. **Each
snapshot is in a single _*.dat_ file, and contains the following:**

The first part of each _*.dat_ file stores all the conservative variables _w_
for all grids which happen to be present at that time in the tree hierarchy.
The Morton-ordered space filling curve through the grid hierarchy, together
with the fixed block size, allows us to perform fully parallel I/O read and
write operations. Each processor writes all grids stored in its local memory
in the order defined by the Morton number of the grid. All processors write
simultaneously to a single file when **typeparIO=1**. Each grid block writes
the _nw_ variables out over the full mesh (without the ghost cell layers).

This part is then followed by all extra info needed to reconstruct/restart the
simulation, which is written at the end _*.dat_ file by the master CPU. The
extra info consists of

    (1)  the boolean representation of the grid structure;
    (2)  info on grid block size per dimension, actually the mesh size _nx^D_
    (3)  equation-specific variable values (i.e. the _eqpar_ array)
    (4)  number of active tree leafs _nleafs_,
    (5)  maximal refinement level present _levmax_,
    (6)  dimensionality NDIM,
    (7)  number of vector components NDIR,
    (8)  number of variables _nw_,
    (9)  number of equation-specific variables _neqpar+nspecialpar_,
    (10) integer time counter _it_,
    (11) global time _t_.

