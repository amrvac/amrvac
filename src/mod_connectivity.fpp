!> This module contains variables that describe the connectivity of the mesh and
!> also data structures for connectivity-related communication.
module mod_connectivity
   implicit none

   integer, parameter :: neighbor_boundary = 1
   integer, parameter :: neighbor_coarse = 2
   integer, parameter :: neighbor_sibling = 3
   integer, parameter :: neighbor_fine = 4

   integer, dimension(:,:,:,:,:), allocatable :: neighbor
   integer, dimension(:,:,:,:,:), allocatable :: neighbor_child
   integer, dimension(:,:,:,:), allocatable :: neighbor_type
   logical, dimension(:,:,:,:), allocatable :: neighbor_active
   integer, dimension(:,:,:,:), allocatable :: neighbor_pole
   !$acc declare create(neighbor, neighbor_type, neighbor_pole, neighbor_child)

   ! grid number array per processor
   integer, dimension(:), allocatable :: igrids
   integer, dimension(:), allocatable :: igrids_active
   integer, dimension(:), allocatable :: igrids_passive
   !$acc declare create(igrids, igrids_active, igrids_passive)
   
   ! number of grids on current processor
   integer :: igridstail
   integer :: igridstail_active
   integer :: igridstail_passive
   !$acc declare create(igridstail, igridstail_active, igridstail_passive)

   integer, dimension(3) :: nrecv_fc, nsend_fc
   ! cc for corner coarse
   integer, dimension(3) :: nrecv_cc, nsend_cc

   ! srl neighbor info 
   type nbinfo_srl_t
      integer                            :: nigrids=0
      integer                            :: iexpand=4   ! realloc with iexpand bigger arrays
      integer, allocatable, dimension(:) :: igrid, i1, i2, i3
    contains
      procedure, non_overridable         :: init_srl
      procedure, non_overridable         :: expand_srl
   end type nbinfo_srl_t

   ! neighbor cpu info structure
   type nbprocs_info_t
      integer              :: nbprocs=0            ! number of neighboring processes overall
      integer              :: nbprocs_srl=0        ! number of neighboring processes at srl
      integer, allocatable :: nbprocs_srl_list(:)  ! list of neighboring ipe at srl
      integer, allocatable :: ipe_to_inbpe_srl(:)  ! inverse to nbprocs_srl_list
      type(nbinfo_srl_t), allocatable   :: srl(:)  ! list of the ipelist for each nbproc
    contains
      procedure, non_overridable :: add_ipe_to_srl_list
      procedure, non_overridable :: init
      procedure, non_overridable :: add_igrid_to_srl
      procedure, non_overridable :: reset
   end type nbprocs_info_t

   type(nbprocs_info_t) :: nbprocs_info
   !$acc declare create(nbprocs_info)

   public :: nbprocs_info_t, nbprocs_info
   
   
 contains

   subroutine expand_srl(self)
     class(nbinfo_srl_t) :: self
     type(nbinfo_srl_t)  :: tmp
     integer             :: size_old

     ! make a copy:
     size_old = size( self%igrid )
     call tmp%init_srl( size_old )
     tmp%nigrids  = self%nigrids
     tmp%igrid    = self%igrid
     tmp%i1       = self%i1
     tmp%i2       = self%i2
     tmp%i3       = self%i3

     ! reallocate and copy back ( cumbersome in fortran :-( )
     call self%init_srl( size(self%igrid) * self%iexpand )
     self%nigrids             = tmp%nigrids
     self%igrid(1:size_old)   = tmp%igrid
     self%i1(1:size_old)      = tmp%i1
     self%i2(1:size_old)      = tmp%i2
     self%i3(1:size_old)      = tmp%i3

   end subroutine expand_srl

   subroutine init_srl(self, nigrids)
     class(nbinfo_srl_t)   :: self
     integer, intent(in)   :: nigrids

     if ( allocated(self%igrid) ) then
        deallocate(self%igrid, self%i1, self%i2, self%i3)
     end if
     allocate(self%igrid(nigrids), self%i1(nigrids), self%i2(nigrids), self%i3(nigrids))

   end subroutine init_srl

   subroutine reset(self)
     class(nbprocs_info_t) :: self
     integer               :: i

     do i=1, self%nbprocs_srl
        self%srl(i)%nigrids=0
     end do

     self%nbprocs             = 0
     self%nbprocs_srl         = 0
     self%ipe_to_inbpe_srl(:) = -1
     
   end subroutine reset
   
   subroutine init(self, npe, nigrids)
     class(nbprocs_info_t) :: self
     integer, intent(in)   :: npe, nigrids
     integer               :: i

     self%nbprocs = npe-1
     
     allocate(self%nbprocs_srl_list(npe-1), &
          self%srl(npe-1))

     allocate(self%ipe_to_inbpe_srl(0:npe-1))
     self%ipe_to_inbpe_srl(:) = -1
     
     do i = 1, npe-1
        call self%srl(i)%init_srl(nigrids)
     end do

   end subroutine init     

   subroutine add_ipe_to_srl_list(self, ipe)
     class(nbprocs_info_t) :: self
     integer, intent(in)   :: ipe

     if (self%ipe_to_inbpe_srl(ipe) == -1) then
     
        ! enlarge counter
        self%nbprocs_srl = self%nbprocs_srl + 1
        ! add the process
        self%nbprocs_srl_list(self%nbprocs_srl) = ipe
        ! add to inverse list
        self%ipe_to_inbpe_srl(ipe) = self%nbprocs_srl
        
     end if
     
   end subroutine add_ipe_to_srl_list

   subroutine add_igrid_to_srl(self, ipe, igrid, i1, i2, i3)
     class(nbprocs_info_t) :: self
     integer, intent(in)   :: ipe, igrid, i1, i2, i3
     integer               :: inbpe

     ! translate to neighbor processor index for srl
     inbpe = self%ipe_to_inbpe_srl(ipe)
     
     ! enlarge counter
     self%srl(inbpe)%nigrids = self%srl(inbpe)%nigrids + 1

     ! need to enlarge storage
     if ( self%srl(inbpe)%nigrids > size( self%srl(inbpe)%igrid ) ) call self%srl(inbpe)%expand_srl

     ! add the data at the counter
     self%srl(inbpe)%igrid( self%srl(inbpe)%nigrids ) = igrid
     self%srl(inbpe)%i1( self%srl(inbpe)%nigrids )    = i1
     self%srl(inbpe)%i2( self%srl(inbpe)%nigrids )    = i2
     self%srl(inbpe)%i3( self%srl(inbpe)%nigrids )    = i3

   end subroutine add_igrid_to_srl
   
end module mod_connectivity
