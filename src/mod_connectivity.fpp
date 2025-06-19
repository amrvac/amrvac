!> This module contains variables that describe the connectivity of the mesh and
!> also data structures for connectivity-related communication.
module mod_connectivity
   implicit none

   integer, parameter :: neighbor_boundary = 1
   integer, parameter :: neighbor_coarse = 2
   integer, parameter :: neighbor_sibling = 3
   integer, parameter :: neighbor_fine = 4

   integer, private, parameter  :: nb_max=100 ! maximum number of neighbor processes
   integer, private, parameter  :: igrids_nb_max=1000 ! maximum number of igrids bordering with particular neighbor ipe
      
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
      integer, dimension(igrids_nb_max) :: igrid, i1, i2, i3
   end type nbinfo_srl_t

   ! neighbor cpu info structure
   type nbprocs_info_t
      ! max ranges, can be replaced by dynamic allocation if bothersome
      integer              :: nbprocs_srl=0  ! number of neighboring processes at srl
      integer              :: nbprocs_srl_list(nb_max)=-1  ! list of neighboring ipe at srl
      type(nbinfo_srl_t)   :: srl(nb_max)
    contains
      procedure, non_overridable :: add_ipe_to_srl_list
   end type nbprocs_info_t

   type(nbprocs_info_t) :: nbprocs_info

   public :: nbprocs_info_t
   
 contains

   subroutine add_ipe_to_srl_list(self, ipe)
     class(nbprocs_info_t) :: self
     integer, intent(in)   :: ipe
     integer               :: i
     logical               :: already_there

     ! ok as long as small loop
     already_there = .false.
     do i = 1, self%nbprocs_srl
        if (self%nbprocs_srl_list(i) == ipe) then
           already_there = .true.
           exit ! already in list, get out
        end if
     end do

     if (.not. already_there) then

        if (self%nbprocs_srl == nb_max) then
           print *, 'reached nb_max neighbors, enlarge limit'
           stop
        else        
           ! enlarge counter
           self%nbprocs_srl = self%nbprocs_srl + 1
           ! add the process
           self%nbprocs_srl_list(self%nbprocs_srl) = ipe
        end if

     end if
     
   end subroutine add_ipe_to_srl_list
   
end module mod_connectivity
