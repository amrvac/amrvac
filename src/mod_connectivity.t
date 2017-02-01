!> This module contains variables that describe the connectivity of the mesh and
!> also data structures for connectivity-related communication.
module mod_connectivity
   implicit none
   save

   integer, dimension(:,:^D&,:), allocatable :: neighbor
   integer, dimension(:,:^D&,:), allocatable :: neighbor_child
   integer, dimension(:^D&,:), allocatable :: neighbor_type
   logical, dimension(:^D&,:), allocatable :: neighbor_active
   logical, dimension(-1:1^D&) :: leveljump
   integer, dimension(:^D&,:), allocatable :: neighbor_pole

   integer, dimension(:), allocatable :: igrids
   integer, dimension(:), allocatable :: igrids_active
   integer, dimension(:), allocatable :: igrids_passive
   integer :: igridstail
   integer :: igridstail_active
   integer :: igridstail_passive

   integer, dimension(^ND) :: nrecv_fc, nsend_fc

   integer :: nrecv_bc_srl, nsend_bc_srl, &
                    nrecv_bc_r, nsend_bc_r, &
                    nrecv_bc_p, nsend_bc_p
!$OMP THREADPRIVATE(leveljump)

end module mod_connectivity
