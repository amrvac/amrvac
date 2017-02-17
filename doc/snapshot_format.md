# Snapshot format

# Header

    integer :: Version number
    integer :: Byte offset where tree information starts
    integer :: Byte offset where block data starts
    double precision :: global_time
    integer :: it
    integer :: nw
    integer :: ndir
    integer :: ndim
    integer :: levmax
    integer :: nleafs
    integer :: nparents
    double precision :: xprobmin(ndim)
    double precision :: xprobmax(ndim)
    integer :: domain_nx(ndim)
    integer :: block_nx(ndim)
    character(len=10) :: w_names(nw)

# Tree information

    logical :: leaf(nleafs+nparents)
    integer :: offset_block(nleafs)
    integer :: refinement_level(nleafs)
    integer :: spatial_index(ndim, nleafs)

# Block

    integer :: ixO^L (2 * ndim indices)
    double precision :: w(ixO^L, nw)

# Version history

## Version 1

Version 1 contained the following information

    1. block data (nw variables)
    2. leaf/parent logical array
    3. header:
        nx^D
        domain_nx^D
        xprobmin^D
        xprobmax^D\}
        nleafs
        levmax
        ndim
        ndir
        nw
        it
        global_time

The idea is that you can reconstruct the full grid when you know the Morton
order used for the leaf/parent logical array.

## Version 2

Version 2 had the same information as version 1, but changes were made to the
Morton order on the coarse grid, causing incompatibility.
