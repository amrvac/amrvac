module mod_connectivity
   use mod_indices, only: ngridshi
   implicit none

   integer, dimension(2,-1:1^D&,ngridshi), save :: neighbor
   integer, dimension(2,0:3^D&,ngridshi), save :: neighbor_child
   integer, dimension(-1:1^D&,ngridshi), save :: neighbor_type
   logical, dimension(-1:1^D&,ngridshi), save :: neighbor_active
   logical, dimension(-1:1^D&), save :: leveljump
   {^IFPHI integer, dimension(-1:1^D&,ngridshi), save :: neighbor_pole}

   integer, dimension(ngridshi), save :: igrids
   integer, save :: igridstail

   integer, dimension(ngridshi), save :: igrids_active
   integer, save :: igridstail_active
   integer, dimension(ngridshi), save :: igrids_passive
   integer, save :: igridstail_passive

   integer, dimension(^ND), save :: nrecv_fc, nsend_fc

   integer, save :: nrecv_bc_srl, nsend_bc_srl, &
                    nrecv_bc_r, nsend_bc_r, &
                    nrecv_bc_p, nsend_bc_p
!$OMP THREADPRIVATE(leveljump,nrecv_bc_srl,nsend_bc_srl,nrecv_bc_r,&
!$OMP& nsend_bc_r,nrecv_bc_p,nsend_bc_p,nrecv_fc,nsend_fc)

end module mod_connectivity
