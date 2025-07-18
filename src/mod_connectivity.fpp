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

   ! phys boundary indices
   integer, dimension(:,:), allocatable :: idphyb
   !$acc declare create(idphyb)
   
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
      integer, allocatable, dimension(:) :: igrid, iencode, isize, ibuf_start
    contains
      procedure, non_overridable         :: init => init_srl
      procedure, non_overridable         :: expand => expand_srl
   end type nbinfo_srl_t

   ! buffer
   type nbinfo_buffer_t
      double precision, allocatable, dimension(:) :: buffer
    contains
      procedure, non_overridable :: alloc => alloc_buffer
   end type nbinfo_buffer_t
   
   type nbinfo_buffer_i_t
      integer, allocatable, dimension(:) :: buffer
    contains
      procedure, non_overridable :: alloc => alloc_buffer_info
   end type nbinfo_buffer_i_t

   ! neighbor cpu info structure
   type nbprocs_info_t
      integer                              :: nbprocs_srl=0        ! number of neighboring processes at srl
      integer, allocatable                 :: nbprocs_srl_list(:)  ! list of neighboring ipe at srl
      integer, allocatable                 :: ipe_to_inbpe_srl(:)  ! inverse to nbprocs_srl_list
      type(nbinfo_srl_t), allocatable      :: srl(:)               ! list of the ipelist for each nbproc
      type(nbinfo_buffer_t), allocatable   :: srl_send(:)          ! double precision send data. One for each nb proc
      type(nbinfo_buffer_t), allocatable   :: srl_rcv(:)           ! double precision receive data
      type(nbinfo_buffer_i_t), allocatable :: srl_info_send(:)     ! info package send
      type(nbinfo_buffer_i_t), allocatable :: srl_info_rcv(:)      ! info package receive
    contains
      procedure, non_overridable :: add_ipe_to_srl_list
      procedure, non_overridable :: init
      procedure, non_overridable :: add_igrid_to_srl
      procedure, non_overridable :: add_to_srl      
      procedure, non_overridable :: reset
      procedure, non_overridable :: alloc_buffers_srl
      procedure, non_overridable :: iencode
      procedure, non_overridable :: idecode
   end type nbprocs_info_t

   type(nbprocs_info_t) :: nbprocs_info
   !$acc declare create(nbprocs_info)

   public  :: nbprocs_info_t, nbprocs_info
   private :: alloc_buffer, alloc_buffer_info, expand_srl, init_srl, reset, init, add_to_srl
   private :: add_ipe_to_srl_list, add_igrid_to_srl, alloc_buffers_srl, iencode, idecode
     
 contains

   subroutine alloc_buffer(self, isize)
     class(nbinfo_buffer_t) :: self
     integer, intent(in)    :: isize

     if (allocated(self%buffer)) deallocate(self%buffer)

     allocate(self%buffer(1:isize))
     
   end subroutine alloc_buffer
   
   subroutine alloc_buffer_info(self, isize)
     class(nbinfo_buffer_i_t) :: self
     integer, intent(in)    :: isize

     if (allocated(self%buffer)) deallocate(self%buffer)

     allocate(self%buffer(1:isize))
     
   end subroutine alloc_buffer_info

   subroutine expand_srl(self)
     class(nbinfo_srl_t) :: self
     type(nbinfo_srl_t)  :: tmp
     ! .. local ..
     integer             :: size_old

     ! make a copy:
     size_old = size( self%igrid )
     call tmp%init( size_old )
     tmp = self

     ! reallocate and copy back ( cumbersome in fortran :-( )
     call self%init( size(self%igrid) * self%iexpand )
     self%nigrids                 = tmp%nigrids
     self%igrid(1:size_old)       = tmp%igrid
     self%iencode(1:size_old)     = tmp%iencode
     self%ibuf_start(1:size_old)  = tmp%ibuf_start
     self%isize(1:size_old)       = tmp%isize

   end subroutine expand_srl

   subroutine init_srl(self, nigrids)
     class(nbinfo_srl_t)   :: self
     integer, intent(in)   :: nigrids

     if ( allocated(self%igrid) ) then
        deallocate(self%igrid, self%iencode, self%ibuf_start, self%isize)
     end if
     allocate(self%igrid(nigrids), self%iencode(nigrids), self%ibuf_start(nigrids), self%isize(nigrids))

   end subroutine init_srl

   subroutine reset(self)
     class(nbprocs_info_t) :: self
     ! .. local ..
     integer               :: i

     do i=1, self%nbprocs_srl
        self%srl(i)%nigrids=0
     end do

     self%nbprocs_srl         = 0
     self%ipe_to_inbpe_srl(:) = -1
     
   end subroutine reset
   
   subroutine init(self, npe, nigrids)
     class(nbprocs_info_t) :: self
     integer, intent(in)   :: npe, nigrids
     ! .. local ..
     integer               :: i

     allocate(self%nbprocs_srl_list(npe-1), &
          self%srl(npe-1))

     allocate(self%ipe_to_inbpe_srl(0:npe-1))
     self%ipe_to_inbpe_srl(:) = -1
     
     do i = 1, npe-1
        call self%srl(i)%init(nigrids)
     end do

   end subroutine init   

   subroutine add_to_srl(self, ipe, igrid, i1, i2, i3)
     class(nbprocs_info_t) :: self
     integer, intent(in)   :: ipe, igrid, i1, i2, i3
     
     call self%add_ipe_to_srl_list(ipe)
     call self%add_igrid_to_srl(ipe, igrid, i1, i2, i3)
     
   end subroutine add_to_srl

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
     ! .. local ..
     integer               :: inbpe, i

     ! translate to neighbor processor index for srl
     inbpe = self%ipe_to_inbpe_srl(ipe)
     
     ! enlarge counter
     self%srl(inbpe)%nigrids = self%srl(inbpe)%nigrids + 1

     ! need to enlarge storage
     if ( self%srl(inbpe)%nigrids > size( self%srl(inbpe)%igrid ) ) call self%srl(inbpe)%expand

     call self%iencode(i1,i2,i3,i)
     ! add the data at the counter
     self%srl(inbpe)%igrid( self%srl(inbpe)%nigrids )    = igrid
     self%srl(inbpe)%iencode( self%srl(inbpe)%nigrids )  = i
     
   end subroutine add_igrid_to_srl

   subroutine alloc_buffers_srl(self, nwgc, &
         ixS_srl_min1, ixS_srl_max1, &
         ixS_srl_min2, ixS_srl_max2, &
         ixS_srl_min3, ixS_srl_max3, &
         ixR_srl_min1, ixR_srl_max1, &
         ixR_srl_min2, ixR_srl_max2, &
         ixR_srl_min3, ixR_srl_max3)
     
     class(nbprocs_info_t)                     :: self
     integer, intent(in)                       :: nwgc
     integer, dimension(-1:2,-1:1), intent(in) :: &
         ixS_srl_min1, ixS_srl_max1, &
         ixS_srl_min2, ixS_srl_max2, &
         ixS_srl_min3, ixS_srl_max3, &
         ixR_srl_min1, ixR_srl_max1, &
         ixR_srl_min2, ixR_srl_max2, &
         ixR_srl_min3, ixR_srl_max3
     ! .. local ..
     integer         :: inb, igrid, isize_S, isize_R
     integer         :: n_i1, n_i2, n_i3, i1, i2, i3
     integer         :: ixSmin1, ixSmin2, ixSmin3
     integer         :: ixRmin1, ixRmin2, ixRmin3
     integer         :: ixSmax1, ixSmax2, ixSmax3
     integer         :: ixRmax1, ixRmax2, ixRmax3
     integer         :: iib1, iib2, iib3, ibuf_start

     if (allocated(self%srl_send)) then
        deallocate( self%srl_send, self%srl_rcv, self%srl_info_send, self%srl_info_rcv )
     end if
     
     allocate( &
          self%srl_send(1:self%nbprocs_srl), &
          self%srl_rcv(1:self%nbprocs_srl), &
          self%srl_info_send(1:self%nbprocs_srl), &
          self%srl_info_rcv(1:self%nbprocs_srl) &
          )

     do inb = 1, self%nbprocs_srl
        call self%srl_info_send(inb)%alloc( 3 * self%srl(inb)%nigrids )
        call self%srl_info_rcv(inb)%alloc( 3 * self%srl(inb)%nigrids ) 
        
        isize_S = 0; isize_R = 0; ibuf_start = 1
        do igrid = 1, self%srl(inb)%nigrids
           iib1 = idphyb(1,self%srl(inb)%igrid(igrid))
           iib2 = idphyb(2,self%srl(inb)%igrid(igrid))
           iib3 = idphyb(3,self%srl(inb)%igrid(igrid))
           call self%idecode( i1, i2, i3, self%srl(inb)%iencode(igrid) )
           
           n_i1=-i1; n_i2=-i2; n_i3=-i3
           
           ixSmin1=ixS_srl_min1(iib1,i1);ixSmin2=ixS_srl_min2(iib2,i2)
           ixSmin3=ixS_srl_min3(iib3,i3);ixSmax1=ixS_srl_max1(iib1,i1)
           ixSmax2=ixS_srl_max2(iib2,i2);ixSmax3=ixS_srl_max3(iib3,i3)
           ixRmin1=ixR_srl_min1(iib1,n_i1);ixRmin2=ixR_srl_min2(iib2,n_i2)
           ixRmin3=ixR_srl_min3(iib3,n_i3);ixRmax1=ixR_srl_max1(iib1,n_i1)
           ixRmax2=ixR_srl_max2(iib2,n_i2);ixRmax3=ixR_srl_max3(iib3,n_i3)

           isize_S = isize_S + (ixSmax1-ixSmin1+1) * (ixSmax2-ixSmin2+1) * (ixSmax3-ixSmin3+1)
           isize_R = isize_R + (ixRmax1-ixRmin1+1) * (ixRmax2-ixRmin2+1) * (ixRmax3-ixRmin3+1)

           self%srl(inb)%ibuf_start(igrid) = ibuf_start
           self%srl(inb)%isize(igrid) = (ixSmax1-ixSmin1+1) * (ixSmax2-ixSmin2+1) * (ixSmax3-ixSmin3+1) * nwgc
           
           ibuf_start = ibuf_start + self%srl(inb)%isize(igrid)

        end do
        call self%srl_send(inb)%alloc(isize_S*nwgc)
        call self%srl_rcv(inb)%alloc(isize_R*nwgc)
     end do

   end subroutine alloc_buffers_srl

   subroutine iencode(self, i1, i2, i3, i)
     class(nbprocs_info_t)                     :: self
     integer, intent(in)                       :: i1, i2, i3
     integer, intent(out)                      :: i
     ! .. local ..
     integer, parameter                        :: i1min=-1, i2min=-1, i3min=-1

     i = 1 + (i1-i1min) + 3 * (i2-i2min) + 9 * (i3-i3min)

   end subroutine iencode
   
   subroutine idecode(self, i1, i2, i3, i)
     class(nbprocs_info_t)                     :: self
     integer, intent(in)                       :: i
     integer, intent(out)                      :: i1, i2, i3
     ! .. local ..
     integer, parameter                        :: i1min=-1, i2min=-1, i3min=-1
     integer                                   :: id, idd
     
     i3  = ceiling(dble(i)/9.0d0) - 1 + i3min
     id  = i - 9 * (i3-i3min)
     i2  = ceiling(dble(id)/3.0d0) -1 + i2min
     idd = id - 3 * (i2-i2min)
     i1  = idd + i1min - 1
     
   end subroutine idecode

 end module mod_connectivity
