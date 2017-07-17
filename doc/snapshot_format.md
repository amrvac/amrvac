# Snapshot format

# Header

    integer           :: Version number
    integer           :: Byte offset where tree information starts
    integer           :: Byte offset where block data starts
    integer           :: nw
    integer           :: ndir
    integer           :: ndim
    integer           :: levmax
    integer           :: nleafs
    integer           :: nparents
    integer           :: it
    double precision  :: global_time
    double precision  :: xprobmin(ndim)
    double precision  :: xprobmax(ndim)
    integer           :: domain_nx(ndim)
    integer           :: block_nx(ndim)
    character(len=16) :: w_names(nw)
    character(len=16) :: physics_type
    ! The physics parameters, such as gamma
    integer           :: n_params
    double precision  :: parameters(n_params)
    character(len=16) :: parameter_names(n_params)
    ! Future implementations will put geometry info below

# Tree information

    logical :: leaf(nleafs+nparents)
    integer :: refinement_level(nleafs)
    integer :: spatial_index(ndim, nleafs)
    integer :: offset_block(nleafs)

# Block

    integer :: n_ghost_lo(ndim) ! number of ghost cells on lower boundaries
    integer :: n_ghost_hi(ndim) ! number of ghost cells on upper boundaries
    ! block_shape = 1-n_ghost_lo:block_nx+n_ghost_hi
    double precision :: w(block_shape, nw)

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
