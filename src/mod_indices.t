!> This module contains global integers and indices
!> @todo Why are there logicals/reals here?
!> @todo Clean this list up
module mod_indices

  implicit none

  !> The maximum number of grid blocks
  !> @todo Don't use a fixed upper bound
  integer, parameter :: ngridshi  = 4000

  !> The maximum number of levels in the grid refinement
  !> @todo Don't use a fixed upper bound
  integer, parameter :: nlevelshi = 13

  !> The number of interleaving sending buffers for ghost cells
  integer, parameter :: npwbuf=2

  integer, save :: ixM^LL

  integer, dimension(nlevelshi), save :: ng^D
  double precision, dimension(nlevelshi), save :: dg^D

  logical, save :: slab

  !> The number of MPI tasks
  integer, save :: npe

  !> The rank of the current MPI task
  integer, save :: mype

  !> The MPI communicator
  integer, save :: icomm

  !> A global MPI error return code
  !> @todo Make local
  integer, save :: ierrmpi

  integer, save :: log_fh
  integer, save :: type_block, type_coarse_block, type_sub_block(2^D&)
  integer, save :: type_block_io, size_block_io, size_block
  {#IFDEF TRANSFORMW
  integer, save :: type_block_io_tf, size_block_io_tf}
  integer, save :: type_subblock_io, type_subblock_x_io
  integer, save :: type_block_xc_io,type_block_xcc_io
  integer, save :: type_block_wc_io,type_block_wcc_io
  integer, save :: itag
  integer, save :: irecv, isend
  integer, dimension(:), allocatable, save :: recvrequest, sendrequest
  integer, dimension(:,:), allocatable, save :: recvstatus, sendstatus

  integer, save :: snapshot, snapshotnext, slice, slicenext, collapseNext, icollapse

  logical, allocatable, dimension(:^D&), save :: patchfalse

  !> For MHD: split B=B0+B1 (time-independent potential field B0 in set_B0.t)
  !>   In spherical coordinates: Dipole/Quadrupole/Octopole strengths
  !>   In general: Busr different from zero interfaces to specialset_B0
  !> @todo Make local to the MHD module only
  logical, save :: B0field
  double precision, save :: Bdip, Bquad, Boct, Busr
end module mod_indices
