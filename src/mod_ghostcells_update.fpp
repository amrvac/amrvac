!> update ghost cells of all blocks including physical boundaries
module mod_ghostcells_update

  implicit none
  ! Special buffer for pole copy
  type wbuffer
     double precision, dimension(:,:,:,:), allocatable :: w
  end type wbuffer

  ! A switch of update physical boundary or not
  logical, public :: bcphys=.true.

  integer :: ixMmin1,ixMmin2,ixMmin3,ixMmax1,ixMmax2,ixMmax3, ixCoGmin1,&
       ixCoGmin2,ixCoGmin3,ixCoGmax1,ixCoGmax2,ixCoGmax3, ixCoMmin1,ixCoMmin2,&
       ixCoMmin3,ixCoMmax1,ixCoMmax2,ixCoMmax3, ixCoGsmin1,ixCoGsmin2,ixCoGsmin3,&
       ixCoGsmax1,ixCoGsmax2,ixCoGsmax3

  ! The number of interleaving sending buffers for ghost cells
  integer, parameter :: npwbuf=2

  ! The first index goes from -1:2, where -1 is used when a block touches the
  ! lower boundary, 1 when a block touches an upper boundary, and 0 a situation
  ! away from boundary conditions, 2 when a block touched both lower and upper
  ! boundary

  ! index ranges to send (S) to sibling blocks, receive (R) from sibling blocks
  integer, dimension(-1:2,-1:1) :: ixS_srl_min1,ixS_srl_min2,ixS_srl_min3,&
       ixS_srl_max1,ixS_srl_max2,ixS_srl_max3, ixR_srl_min1,ixR_srl_min2,&
       ixR_srl_min3,ixR_srl_max1,ixR_srl_max2,ixR_srl_max3
  !$acc declare create(ixS_srl_min1,ixS_srl_min2,ixS_srl_min3)
  !$acc declare create(ixS_srl_max1,ixS_srl_max2,ixS_srl_max3, ixR_srl_min1,ixR_srl_min2)
  !$acc declare create(ixR_srl_min3,ixR_srl_max1,ixR_srl_max2,ixR_srl_max3)

  ! index ranges of staggered variables to send (S) to sibling blocks, receive (R) from sibling blocks
  integer, dimension(3,-1:1) :: ixS_srl_stg_min1,ixS_srl_stg_min2,&
       ixS_srl_stg_min3,ixS_srl_stg_max1,ixS_srl_stg_max2,ixS_srl_stg_max3,&
       ixR_srl_stg_min1,ixR_srl_stg_min2,ixR_srl_stg_min3,ixR_srl_stg_max1,&
       ixR_srl_stg_max2,ixR_srl_stg_max3

  ! index ranges to send (S) restricted (r) ghost cells to coarser blocks
  integer, dimension(-1:1,-1:1) :: ixS_r_min1,ixS_r_min2,ixS_r_min3,ixS_r_max1,&
       ixS_r_max2,ixS_r_max3

  ! index ranges of staggered variables to send (S) restricted (r) ghost cells to coarser blocks
  integer, dimension(3,-1:1) :: ixS_r_stg_min1,ixS_r_stg_min2,ixS_r_stg_min3,&
       ixS_r_stg_max1,ixS_r_stg_max2,ixS_r_stg_max3

  ! index ranges to receive restriced ghost cells from finer blocks
  integer, dimension(-1:1, 0:3) :: ixR_r_min1,ixR_r_min2,ixR_r_min3,ixR_r_max1,&
       ixR_r_max2,ixR_r_max3

  ! index ranges of staggered variables to receive restriced ghost cells from finer blocks
  integer, dimension(3,0:3)  :: ixR_r_stg_min1,ixR_r_stg_min2,ixR_r_stg_min3,&
       ixR_r_stg_max1,ixR_r_stg_max2,ixR_r_stg_max3

  ! send prolongated (p) ghost cells to finer blocks, receive prolongated from coarser blocks
  integer, dimension(-1:1, 0:3) :: ixS_p_min1,ixS_p_min2,ixS_p_min3,ixS_p_max1,&
       ixS_p_max2,ixS_p_max3, ixR_p_min1,ixR_p_min2,ixR_p_min3,ixR_p_max1,&
       ixR_p_max2,ixR_p_max3

  ! send prolongated (p) staggered ghost cells to finer blocks, receive prolongated from coarser blocks
  integer, dimension(3,0:3)  :: ixS_p_stg_min1,ixS_p_stg_min2,ixS_p_stg_min3,&
       ixS_p_stg_max1,ixS_p_stg_max2,ixS_p_stg_max3, ixR_p_stg_min1,&
       ixR_p_stg_min2,ixR_p_stg_min3,ixR_p_stg_max1,ixR_p_stg_max2,&
       ixR_p_stg_max3

  ! number of MPI receive-send pairs, srl: same refinement level; r: restrict; p: prolong
  integer :: nrecv_bc_srl, nsend_bc_srl, nrecv_bc_r, nsend_bc_r, nrecv_bc_p,&
       nsend_bc_p

  ! record index position of buffer arrays
  integer :: ibuf_send_srl, ibuf_recv_srl, ibuf_send_r, ibuf_recv_r,&
       ibuf_send_p, ibuf_recv_p

  ! count of times of send and receive
  integer :: isend_srl, irecv_srl, isend_r, irecv_r, isend_p, irecv_p

  ! count of times of send and receive for cell center ghost cells
  integer :: isend_c, irecv_c

  ! tag of MPI send and recv
  integer, private :: itag

  ! total sizes = cell-center normal flux + stagger-grid flux of send and receive
  integer, dimension(-1:1,-1:1,-1:1) :: sizes_srl_send_total,&
       sizes_srl_recv_total

  ! sizes of buffer arrays for center-grid variable for siblings and restrict
  integer, dimension(:), allocatable :: recvrequest_c_sr, sendrequest_c_sr
  integer, dimension(:,:), allocatable :: recvstatus_c_sr, sendstatus_c_sr

  ! sizes of buffer arrays for center-grid variable for prolongation
  integer, dimension(:), allocatable :: recvrequest_c_p, sendrequest_c_p
  integer, dimension(:,:), allocatable :: recvstatus_c_p, sendstatus_c_p

  ! sizes of buffer arrays for stagger-grid variable
  integer, dimension(3,-1:1,-1:1,-1:1) :: sizes_srl_send_stg,&
       sizes_srl_recv_stg

  integer, dimension(:), allocatable :: recvrequest_srl, sendrequest_srl
  integer, dimension(:,:), allocatable :: recvstatus_srl, sendstatus_srl

  ! buffer arrays for send and receive of siblings, allocate in build_connectivity
  double precision, dimension(:), allocatable :: recvbuffer_srl,&
       sendbuffer_srl

  integer, dimension(:), allocatable :: recvrequest_r, sendrequest_r
  integer, dimension(:,:), allocatable :: recvstatus_r, sendstatus_r

  ! buffer arrays for send and receive in restriction
  double precision, dimension(:), allocatable :: recvbuffer_r, sendbuffer_r

  integer, dimension(:), allocatable :: recvrequest_p, sendrequest_p
  integer, dimension(:,:), allocatable :: recvstatus_p, sendstatus_p

  ! buffer arrays for send and receive in prolongation
  double precision, dimension(:), allocatable :: recvbuffer_p, sendbuffer_p

  ! sizes to allocate buffer arrays for send and receive for restriction
  integer, dimension(-1:1,-1:1,-1:1)     :: sizes_r_send_total
  integer, dimension(0:3,0:3,0:3)      :: sizes_r_recv_total
  integer, dimension(3,-1:1,-1:1,-1:1) :: sizes_r_send_stg
  integer, dimension(3,0:3,0:3,0:3)  :: sizes_r_recv_stg

  ! sizes to allocate buffer arrays for send and receive for restriction
  integer, dimension(0:3,0:3,0:3)      :: sizes_p_send_total,&
       sizes_p_recv_total
  integer, dimension(3,0:3,0:3,0:3)  :: sizes_p_send_stg, sizes_p_recv_stg

  ! There are two variants, _f indicates that all flux variables are filled,
  ! whereas _p means that part of the variables is filled
  ! Furthermore _r_ stands for restrict, _p_ for prolongation.
  integer, dimension(-1:2,-1:2,-1:2,-1:1,-1:1,-1:1), target :: type_send_srl_f,&
       type_recv_srl_f
  integer, dimension(-1:1,-1:1,-1:1,-1:1,-1:1,-1:1), target :: type_send_r_f
  integer, dimension(-1:1,-1:1,-1:1, 0:3,0:3,0:3), target :: type_recv_r_f,&
       type_send_p_f, type_recv_p_f
  integer, dimension(-1:2,-1:2,-1:2,-1:1,-1:1,-1:1),&
       target :: type_send_srl_p1, type_recv_srl_p1
  integer, dimension(-1:1,-1:1,-1:1,-1:1,-1:1,-1:1), target :: type_send_r_p1
  integer, dimension(-1:1,-1:1,-1:1, 0:3,0:3,0:3), target :: type_recv_r_p1,&
       type_send_p_p1, type_recv_p_p1
  integer, dimension(-1:2,-1:2,-1:2,-1:1,-1:1,-1:1),&
       target :: type_send_srl_p2, type_recv_srl_p2
  integer, dimension(-1:1,-1:1,-1:1,-1:1,-1:1,-1:1), target :: type_send_r_p2
  integer, dimension(-1:1,-1:1,-1:1, 0:3,0:3,0:3), target :: type_recv_r_p2,&
       type_send_p_p2, type_recv_p_p2
  integer, dimension(:,:,:,:,:,:), pointer :: type_send_srl, type_recv_srl,&
       type_send_r
  integer, dimension(:,:,:,:,:,:), pointer :: type_recv_r, type_send_p,&
       type_recv_p

contains

  subroutine init_bc()
    use mod_global_parameters
    use mod_physics, only: phys_req_diagonal, physics_type
    use mod_comm_lib, only: mpistop

    integer :: nghostcellsCo, interpolation_order
    integer :: nx1,nx2,nx3, nxCo1,nxCo2,nxCo3, ixGmin1,ixGmin2,ixGmin3,ixGmax1,&
         ixGmax2,ixGmax3, i1,i2,i3, ic1,ic2,ic3, inc1,inc2,inc3, idir

    ixGmin1=ixGlo1;ixGmin2=ixGlo2;ixGmin3=ixGlo3;ixGmax1=ixGhi1
    ixGmax2=ixGhi2;ixGmax3=ixGhi3;
    ixMmin1=ixGmin1+nghostcells;ixMmin2=ixGmin2+nghostcells
    ixMmin3=ixGmin3+nghostcells;ixMmax1=ixGmax1-nghostcells
    ixMmax2=ixGmax2-nghostcells;ixMmax3=ixGmax3-nghostcells;
    ixCoGmin1=1;ixCoGmin2=1;ixCoGmin3=1;
    ixCoGmax1=(ixGhi1-2*nghostcells)/2+2*nghostcells
    ixCoGmax2=(ixGhi2-2*nghostcells)/2+2*nghostcells
    ixCoGmax3=(ixGhi3-2*nghostcells)/2+2*nghostcells;
    ixCoGsmin1=0;ixCoGsmin2=0;ixCoGsmin3=0;
    ixCoGsmax1=ixCoGmax1;ixCoGsmax2=ixCoGmax2;ixCoGsmax3=ixCoGmax3;

    ixCoMmin1=ixCoGmin1+nghostcells;ixCoMmin2=ixCoGmin2+nghostcells
    ixCoMmin3=ixCoGmin3+nghostcells;ixCoMmax1=ixCoGmax1-nghostcells
    ixCoMmax2=ixCoGmax2-nghostcells;ixCoMmax3=ixCoGmax3-nghostcells;

    nx1=ixMmax1-ixMmin1+1;nx2=ixMmax2-ixMmin2+1;nx3=ixMmax3-ixMmin3+1;
    nxCo1=nx1/2;nxCo2=nx2/2;nxCo3=nx3/2;

    if(ghost_copy) then
       interpolation_order=1
    else
       interpolation_order=2
    end if
    nghostcellsCo=int((nghostcells+1)/2)

    if (nghostcellsCo+interpolation_order-1>nghostcells) then
       call mpistop("interpolation order for prolongation in getbc too high")
    end if

    ! (iib,i) index has following meanings: iib = 0 means it is not at any physical boundary
    ! iib=-1 means it is at the minimum side of a physical boundary
    ! iib= 1 means it is at the maximum side of a physical boundary
    ! i=-1 means subregion prepared for the neighbor at its minimum side
    ! i= 1 means subregion prepared for the neighbor at its maximum side

    ixS_srl_min1(:,-1)=ixMmin1
    ixS_srl_min1(:, 0)=ixMmin1
    ixS_srl_min1(:, 1)=ixMmax1+1-nghostcells
    ixS_srl_max1(:,-1)=ixMmin1-1+nghostcells
    ixS_srl_max1(:, 0)=ixMmax1
    ixS_srl_max1(:, 1)=ixMmax1

    ixR_srl_min1(:,-1)=1
    ixR_srl_min1(:, 0)=ixMmin1
    ixR_srl_min1(:, 1)=ixMmax1+1
    ixR_srl_max1(:,-1)=nghostcells
    ixR_srl_max1(:, 0)=ixMmax1
    ixR_srl_max1(:, 1)=ixGmax1

    ixS_r_min1(:,-1)=ixCoMmin1
    ixS_r_min1(:, 0)=ixCoMmin1
    ixS_r_min1(:, 1)=ixCoMmax1+1-nghostcells
    ixS_r_max1(:,-1)=ixCoMmin1-1+nghostcells
    ixS_r_max1(:, 0)=ixCoMmax1
    ixS_r_max1(:, 1)=ixCoMmax1

    ixR_r_min1(:, 0)=1
    ixR_r_min1(:, 1)=ixMmin1
    ixR_r_min1(:, 2)=ixMmin1+nxCo1
    ixR_r_min1(:, 3)=ixMmax1+1
    ixR_r_max1(:, 0)=nghostcells
    ixR_r_max1(:, 1)=ixMmin1-1+nxCo1
    ixR_r_max1(:, 2)=ixMmax1
    ixR_r_max1(:, 3)=ixGmax1

    ixS_p_min1(:, 0)=ixMmin1-(interpolation_order-1)
    ixS_p_min1(:, 1)=ixMmin1-(interpolation_order-1)
    ixS_p_min1(:, 2)=ixMmin1+nxCo1-nghostcellsCo-(interpolation_order-1)
    ixS_p_min1(:, 3)=ixMmax1+1-nghostcellsCo-(interpolation_order-1)
    ixS_p_max1(:, 0)=ixMmin1-1+nghostcellsCo+(interpolation_order-1)
    ixS_p_max1(:, 1)=ixMmin1-1+nxCo1+nghostcellsCo+(interpolation_order-1)
    ixS_p_max1(:, 2)=ixMmax1+(interpolation_order-1)
    ixS_p_max1(:, 3)=ixMmax1+(interpolation_order-1)

    if(.not.phys_req_diagonal) then
       ! exclude ghost-cell region when diagonal cells are unknown
       ixS_p_min1(:, 0)=ixMmin1
       ixS_p_max1(:, 3)=ixMmax1
       ixS_p_max1(:, 1)=ixMmin1-1+nxCo1+(interpolation_order-1)
       ixS_p_min1(:, 2)=ixMmin1+nxCo1-(interpolation_order-1)
    end if

    ixR_p_min1(:, 0)=ixCoMmin1-nghostcellsCo-(interpolation_order-1)
    ixR_p_min1(:, 1)=ixCoMmin1-(interpolation_order-1)
    ixR_p_min1(:, 2)=ixCoMmin1-nghostcellsCo-(interpolation_order-1)
    ixR_p_min1(:, 3)=ixCoMmax1+1-(interpolation_order-1)
    ixR_p_max1(:, 0)=nghostcells+(interpolation_order-1)
    ixR_p_max1(:, 1)=ixCoMmax1+nghostcellsCo+(interpolation_order-1)
    ixR_p_max1(:, 2)=ixCoMmax1+(interpolation_order-1)
    ixR_p_max1(:, 3)=ixCoMmax1+nghostcellsCo+(interpolation_order-1)

    if(.not.phys_req_diagonal) then
       ixR_p_max1(:, 0)=nghostcells
       ixR_p_min1(:, 3)=ixCoMmax1+1
       ixR_p_max1(:, 1)=ixCoMmax1+(interpolation_order-1)
       ixR_p_min1(:, 2)=ixCoMmin1-(interpolation_order-1)
    end if



    ixS_srl_min2(:,-1)=ixMmin2
    ixS_srl_min2(:, 0)=ixMmin2
    ixS_srl_min2(:, 1)=ixMmax2+1-nghostcells
    ixS_srl_max2(:,-1)=ixMmin2-1+nghostcells
    ixS_srl_max2(:, 0)=ixMmax2
    ixS_srl_max2(:, 1)=ixMmax2

    ixR_srl_min2(:,-1)=1
    ixR_srl_min2(:, 0)=ixMmin2
    ixR_srl_min2(:, 1)=ixMmax2+1
    ixR_srl_max2(:,-1)=nghostcells
    ixR_srl_max2(:, 0)=ixMmax2
    ixR_srl_max2(:, 1)=ixGmax2

    ixS_r_min2(:,-1)=ixCoMmin2
    ixS_r_min2(:, 0)=ixCoMmin2
    ixS_r_min2(:, 1)=ixCoMmax2+1-nghostcells
    ixS_r_max2(:,-1)=ixCoMmin2-1+nghostcells
    ixS_r_max2(:, 0)=ixCoMmax2
    ixS_r_max2(:, 1)=ixCoMmax2

    ixR_r_min2(:, 0)=1
    ixR_r_min2(:, 1)=ixMmin2
    ixR_r_min2(:, 2)=ixMmin2+nxCo2
    ixR_r_min2(:, 3)=ixMmax2+1
    ixR_r_max2(:, 0)=nghostcells
    ixR_r_max2(:, 1)=ixMmin2-1+nxCo2
    ixR_r_max2(:, 2)=ixMmax2
    ixR_r_max2(:, 3)=ixGmax2

    ixS_p_min2(:, 0)=ixMmin2-(interpolation_order-1)
    ixS_p_min2(:, 1)=ixMmin2-(interpolation_order-1)
    ixS_p_min2(:, 2)=ixMmin2+nxCo2-nghostcellsCo-(interpolation_order-1)
    ixS_p_min2(:, 3)=ixMmax2+1-nghostcellsCo-(interpolation_order-1)
    ixS_p_max2(:, 0)=ixMmin2-1+nghostcellsCo+(interpolation_order-1)
    ixS_p_max2(:, 1)=ixMmin2-1+nxCo2+nghostcellsCo+(interpolation_order-1)
    ixS_p_max2(:, 2)=ixMmax2+(interpolation_order-1)
    ixS_p_max2(:, 3)=ixMmax2+(interpolation_order-1)

    if(.not.phys_req_diagonal) then
       ! exclude ghost-cell region when diagonal cells are unknown
       ixS_p_min2(:, 0)=ixMmin2
       ixS_p_max2(:, 3)=ixMmax2
       ixS_p_max2(:, 1)=ixMmin2-1+nxCo2+(interpolation_order-1)
       ixS_p_min2(:, 2)=ixMmin2+nxCo2-(interpolation_order-1)
    end if

    ixR_p_min2(:, 0)=ixCoMmin2-nghostcellsCo-(interpolation_order-1)
    ixR_p_min2(:, 1)=ixCoMmin2-(interpolation_order-1)
    ixR_p_min2(:, 2)=ixCoMmin2-nghostcellsCo-(interpolation_order-1)
    ixR_p_min2(:, 3)=ixCoMmax2+1-(interpolation_order-1)
    ixR_p_max2(:, 0)=nghostcells+(interpolation_order-1)
    ixR_p_max2(:, 1)=ixCoMmax2+nghostcellsCo+(interpolation_order-1)
    ixR_p_max2(:, 2)=ixCoMmax2+(interpolation_order-1)
    ixR_p_max2(:, 3)=ixCoMmax2+nghostcellsCo+(interpolation_order-1)

    if(.not.phys_req_diagonal) then
       ixR_p_max2(:, 0)=nghostcells
       ixR_p_min2(:, 3)=ixCoMmax2+1
       ixR_p_max2(:, 1)=ixCoMmax2+(interpolation_order-1)
       ixR_p_min2(:, 2)=ixCoMmin2-(interpolation_order-1)
    end if



    ixS_srl_min3(:,-1)=ixMmin3
    ixS_srl_min3(:, 0)=ixMmin3
    ixS_srl_min3(:, 1)=ixMmax3+1-nghostcells
    ixS_srl_max3(:,-1)=ixMmin3-1+nghostcells
    ixS_srl_max3(:, 0)=ixMmax3
    ixS_srl_max3(:, 1)=ixMmax3

    ixR_srl_min3(:,-1)=1
    ixR_srl_min3(:, 0)=ixMmin3
    ixR_srl_min3(:, 1)=ixMmax3+1
    ixR_srl_max3(:,-1)=nghostcells
    ixR_srl_max3(:, 0)=ixMmax3
    ixR_srl_max3(:, 1)=ixGmax3

    ixS_r_min3(:,-1)=ixCoMmin3
    ixS_r_min3(:, 0)=ixCoMmin3
    ixS_r_min3(:, 1)=ixCoMmax3+1-nghostcells
    ixS_r_max3(:,-1)=ixCoMmin3-1+nghostcells
    ixS_r_max3(:, 0)=ixCoMmax3
    ixS_r_max3(:, 1)=ixCoMmax3

    ixR_r_min3(:, 0)=1
    ixR_r_min3(:, 1)=ixMmin3
    ixR_r_min3(:, 2)=ixMmin3+nxCo3
    ixR_r_min3(:, 3)=ixMmax3+1
    ixR_r_max3(:, 0)=nghostcells
    ixR_r_max3(:, 1)=ixMmin3-1+nxCo3
    ixR_r_max3(:, 2)=ixMmax3
    ixR_r_max3(:, 3)=ixGmax3

    ixS_p_min3(:, 0)=ixMmin3-(interpolation_order-1)
    ixS_p_min3(:, 1)=ixMmin3-(interpolation_order-1)
    ixS_p_min3(:, 2)=ixMmin3+nxCo3-nghostcellsCo-(interpolation_order-1)
    ixS_p_min3(:, 3)=ixMmax3+1-nghostcellsCo-(interpolation_order-1)
    ixS_p_max3(:, 0)=ixMmin3-1+nghostcellsCo+(interpolation_order-1)
    ixS_p_max3(:, 1)=ixMmin3-1+nxCo3+nghostcellsCo+(interpolation_order-1)
    ixS_p_max3(:, 2)=ixMmax3+(interpolation_order-1)
    ixS_p_max3(:, 3)=ixMmax3+(interpolation_order-1)

    if(.not.phys_req_diagonal) then
       ! exclude ghost-cell region when diagonal cells are unknown
       ixS_p_min3(:, 0)=ixMmin3
       ixS_p_max3(:, 3)=ixMmax3
       ixS_p_max3(:, 1)=ixMmin3-1+nxCo3+(interpolation_order-1)
       ixS_p_min3(:, 2)=ixMmin3+nxCo3-(interpolation_order-1)
    end if

    ixR_p_min3(:, 0)=ixCoMmin3-nghostcellsCo-(interpolation_order-1)
    ixR_p_min3(:, 1)=ixCoMmin3-(interpolation_order-1)
    ixR_p_min3(:, 2)=ixCoMmin3-nghostcellsCo-(interpolation_order-1)
    ixR_p_min3(:, 3)=ixCoMmax3+1-(interpolation_order-1)
    ixR_p_max3(:, 0)=nghostcells+(interpolation_order-1)
    ixR_p_max3(:, 1)=ixCoMmax3+nghostcellsCo+(interpolation_order-1)
    ixR_p_max3(:, 2)=ixCoMmax3+(interpolation_order-1)
    ixR_p_max3(:, 3)=ixCoMmax3+nghostcellsCo+(interpolation_order-1)

    if(.not.phys_req_diagonal) then
       ixR_p_max3(:, 0)=nghostcells
       ixR_p_min3(:, 3)=ixCoMmax3+1
       ixR_p_max3(:, 1)=ixCoMmax3+(interpolation_order-1)
       ixR_p_min3(:, 2)=ixCoMmin3-(interpolation_order-1)
    end if



    if (stagger_grid) then
       allocate(pole_buf%ws(ixGslo1:ixGshi1,ixGslo2:ixGshi2,ixGslo3:ixGshi3,&
            nws))
       ! Staggered (face-allocated) variables
       do idir=1,ndim
          ixS_srl_stg_min1(idir,-1)=ixMmin1-kr(idir,1)
          ixS_srl_stg_max1(idir,-1)=ixMmin1-1+nghostcells
          ixS_srl_stg_min1(idir,0) =ixMmin1-kr(idir,1)
          ixS_srl_stg_max1(idir,0) =ixMmax1
          ixS_srl_stg_min1(idir,1) =ixMmax1-nghostcells+1-kr(idir,1)
          ixS_srl_stg_max1(idir,1) =ixMmax1

          ixR_srl_stg_min1(idir,-1)=1-kr(idir,1)
          ixR_srl_stg_max1(idir,-1)=nghostcells
          ixR_srl_stg_min1(idir,0) =ixMmin1-kr(idir,1)
          ixR_srl_stg_max1(idir,0) =ixMmax1
          ixR_srl_stg_min1(idir,1) =ixMmax1+1-kr(idir,1)
          ixR_srl_stg_max1(idir,1) =ixGmax1

          ixS_r_stg_min1(idir,-1)=ixCoMmin1-kr(idir,1)
          ixS_r_stg_max1(idir,-1)=ixCoMmin1-1+nghostcells
          ixS_r_stg_min1(idir,0) =ixCoMmin1-kr(idir,1)
          ixS_r_stg_max1(idir,0) =ixCoMmax1
          ixS_r_stg_min1(idir,1) =ixCoMmax1+1-nghostcells-kr(idir,1)
          ixS_r_stg_max1(idir,1) =ixCoMmax1

          ixR_r_stg_min1(idir,0)=1-kr(idir,1)
          ixR_r_stg_max1(idir,0)=nghostcells
          ixR_r_stg_min1(idir,1)=ixMmin1-kr(idir,1)
          ixR_r_stg_max1(idir,1)=ixMmin1-1+nxCo1
          ixR_r_stg_min1(idir,2)=ixMmin1+nxCo1-kr(idir,1)
          ixR_r_stg_max1(idir,2)=ixMmax1
          ixR_r_stg_min1(idir,3)=ixMmax1+1-kr(idir,1)
          ixR_r_stg_max1(idir,3)=ixGmax1

          ixS_srl_stg_min2(idir,-1)=ixMmin2-kr(idir,2)
          ixS_srl_stg_max2(idir,-1)=ixMmin2-1+nghostcells
          ixS_srl_stg_min2(idir,0) =ixMmin2-kr(idir,2)
          ixS_srl_stg_max2(idir,0) =ixMmax2
          ixS_srl_stg_min2(idir,1) =ixMmax2-nghostcells+1-kr(idir,2)
          ixS_srl_stg_max2(idir,1) =ixMmax2

          ixR_srl_stg_min2(idir,-1)=1-kr(idir,2)
          ixR_srl_stg_max2(idir,-1)=nghostcells
          ixR_srl_stg_min2(idir,0) =ixMmin2-kr(idir,2)
          ixR_srl_stg_max2(idir,0) =ixMmax2
          ixR_srl_stg_min2(idir,1) =ixMmax2+1-kr(idir,2)
          ixR_srl_stg_max2(idir,1) =ixGmax2

          ixS_r_stg_min2(idir,-1)=ixCoMmin2-kr(idir,2)
          ixS_r_stg_max2(idir,-1)=ixCoMmin2-1+nghostcells
          ixS_r_stg_min2(idir,0) =ixCoMmin2-kr(idir,2)
          ixS_r_stg_max2(idir,0) =ixCoMmax2
          ixS_r_stg_min2(idir,1) =ixCoMmax2+1-nghostcells-kr(idir,2)
          ixS_r_stg_max2(idir,1) =ixCoMmax2

          ixR_r_stg_min2(idir,0)=1-kr(idir,2)
          ixR_r_stg_max2(idir,0)=nghostcells
          ixR_r_stg_min2(idir,1)=ixMmin2-kr(idir,2)
          ixR_r_stg_max2(idir,1)=ixMmin2-1+nxCo2
          ixR_r_stg_min2(idir,2)=ixMmin2+nxCo2-kr(idir,2)
          ixR_r_stg_max2(idir,2)=ixMmax2
          ixR_r_stg_min2(idir,3)=ixMmax2+1-kr(idir,2)
          ixR_r_stg_max2(idir,3)=ixGmax2

          ixS_srl_stg_min3(idir,-1)=ixMmin3-kr(idir,3)
          ixS_srl_stg_max3(idir,-1)=ixMmin3-1+nghostcells
          ixS_srl_stg_min3(idir,0) =ixMmin3-kr(idir,3)
          ixS_srl_stg_max3(idir,0) =ixMmax3
          ixS_srl_stg_min3(idir,1) =ixMmax3-nghostcells+1-kr(idir,3)
          ixS_srl_stg_max3(idir,1) =ixMmax3

          ixR_srl_stg_min3(idir,-1)=1-kr(idir,3)
          ixR_srl_stg_max3(idir,-1)=nghostcells
          ixR_srl_stg_min3(idir,0) =ixMmin3-kr(idir,3)
          ixR_srl_stg_max3(idir,0) =ixMmax3
          ixR_srl_stg_min3(idir,1) =ixMmax3+1-kr(idir,3)
          ixR_srl_stg_max3(idir,1) =ixGmax3

          ixS_r_stg_min3(idir,-1)=ixCoMmin3-kr(idir,3)
          ixS_r_stg_max3(idir,-1)=ixCoMmin3-1+nghostcells
          ixS_r_stg_min3(idir,0) =ixCoMmin3-kr(idir,3)
          ixS_r_stg_max3(idir,0) =ixCoMmax3
          ixS_r_stg_min3(idir,1) =ixCoMmax3+1-nghostcells-kr(idir,3)
          ixS_r_stg_max3(idir,1) =ixCoMmax3

          ixR_r_stg_min3(idir,0)=1-kr(idir,3)
          ixR_r_stg_max3(idir,0)=nghostcells
          ixR_r_stg_min3(idir,1)=ixMmin3-kr(idir,3)
          ixR_r_stg_max3(idir,1)=ixMmin3-1+nxCo3
          ixR_r_stg_min3(idir,2)=ixMmin3+nxCo3-kr(idir,3)
          ixR_r_stg_max3(idir,2)=ixMmax3
          ixR_r_stg_min3(idir,3)=ixMmax3+1-kr(idir,3)
          ixR_r_stg_max3(idir,3)=ixGmax3

          if (idir==1) then
             ! Parallel components

             ixS_p_stg_min1(idir,0)=ixMmin1-1 ! -1 to make redundant
             ixS_p_stg_max1(idir,0)=ixMmin1-1+nghostcellsCo
             ixS_p_stg_min1(idir,1)=ixMmin1-1 ! -1 to make redundant
             ixS_p_stg_max1(idir,1)=ixMmin1-1+nxCo1+nghostcellsCo
             ixS_p_stg_min1(idir,2)=ixMmax1-nxCo1-nghostcellsCo
             ixS_p_stg_max1(idir,2)=ixMmax1
             ixS_p_stg_min1(idir,3)=ixMmax1-nghostcellsCo
             ixS_p_stg_max1(idir,3)=ixMmax1

             ixR_p_stg_min1(idir,0)=ixCoMmin1-1-nghostcellsCo
             ixR_p_stg_max1(idir,0)=ixCoMmin1-1
             ixR_p_stg_min1(idir,1)=ixCoMmin1-1 ! -1 to make redundant
             ixR_p_stg_max1(idir,1)=ixCoMmax1+nghostcellsCo
             ixR_p_stg_min1(idir,2)=ixCoMmin1-1-nghostcellsCo
             ixR_p_stg_max1(idir,2)=ixCoMmax1
             ixR_p_stg_min1(idir,3)=ixCoMmax1+1-1 ! -1 to make redundant
             ixR_p_stg_max1(idir,3)=ixCoMmax1+nghostcellsCo


             ixS_p_stg_min1(idir,0)=ixMmin1-1 ! -1 to make redundant
             ixS_p_stg_max1(idir,0)=ixMmin1-1+nghostcellsCo
             ixS_p_stg_min1(idir,1)=ixMmin1-1 ! -1 to make redundant
             ixS_p_stg_max1(idir,1)=ixMmin1-1+nxCo1+nghostcellsCo
             ixS_p_stg_min1(idir,2)=ixMmax1-nxCo1-nghostcellsCo
             ixS_p_stg_max1(idir,2)=ixMmax1
             ixS_p_stg_min1(idir,3)=ixMmax1-nghostcellsCo
             ixS_p_stg_max1(idir,3)=ixMmax1

             ixR_p_stg_min1(idir,0)=ixCoMmin1-1-nghostcellsCo
             ixR_p_stg_max1(idir,0)=ixCoMmin1-1
             ixR_p_stg_min1(idir,1)=ixCoMmin1-1 ! -1 to make redundant
             ixR_p_stg_max1(idir,1)=ixCoMmax1+nghostcellsCo
             ixR_p_stg_min1(idir,2)=ixCoMmin1-1-nghostcellsCo
             ixR_p_stg_max1(idir,2)=ixCoMmax1
             ixR_p_stg_min1(idir,3)=ixCoMmax1+1-1 ! -1 to make redundant
             ixR_p_stg_max1(idir,3)=ixCoMmax1+nghostcellsCo


             ixS_p_stg_min1(idir,0)=ixMmin1-1 ! -1 to make redundant
             ixS_p_stg_max1(idir,0)=ixMmin1-1+nghostcellsCo
             ixS_p_stg_min1(idir,1)=ixMmin1-1 ! -1 to make redundant
             ixS_p_stg_max1(idir,1)=ixMmin1-1+nxCo1+nghostcellsCo
             ixS_p_stg_min1(idir,2)=ixMmax1-nxCo1-nghostcellsCo
             ixS_p_stg_max1(idir,2)=ixMmax1
             ixS_p_stg_min1(idir,3)=ixMmax1-nghostcellsCo
             ixS_p_stg_max1(idir,3)=ixMmax1

             ixR_p_stg_min1(idir,0)=ixCoMmin1-1-nghostcellsCo
             ixR_p_stg_max1(idir,0)=ixCoMmin1-1
             ixR_p_stg_min1(idir,1)=ixCoMmin1-1 ! -1 to make redundant
             ixR_p_stg_max1(idir,1)=ixCoMmax1+nghostcellsCo
             ixR_p_stg_min1(idir,2)=ixCoMmin1-1-nghostcellsCo
             ixR_p_stg_max1(idir,2)=ixCoMmax1
             ixR_p_stg_min1(idir,3)=ixCoMmax1+1-1 ! -1 to make redundant
             ixR_p_stg_max1(idir,3)=ixCoMmax1+nghostcellsCo

          else

             ! Perpendicular component
             ixS_p_stg_min1(idir,0)=ixMmin1
             ixS_p_stg_max1(idir,0)=ixMmin1-1+nghostcellsCo+&
                  (interpolation_order-1)
             ixS_p_stg_min1(idir,1)=ixMmin1
             ixS_p_stg_max1(idir,1)=ixMmin1-1+nxCo1+nghostcellsCo+&
                  (interpolation_order-1)
             ixS_p_stg_min1(idir,2)=ixMmax1+&
                  1-nxCo1-nghostcellsCo-(interpolation_order-1)
             ixS_p_stg_max1(idir,2)=ixMmax1
             ixS_p_stg_min1(idir,3)=ixMmax1+&
                  1-nghostcellsCo-(interpolation_order-1)
             ixS_p_stg_max1(idir,3)=ixMmax1

             ixR_p_stg_min1(idir,0)=ixCoMmin1-nghostcellsCo-(interpolation_order-&
                  1)
             ixR_p_stg_max1(idir,0)=ixCoMmin1-1
             ixR_p_stg_min1(idir,1)=ixCoMmin1
             ixR_p_stg_max1(idir,1)=ixCoMmax1+nghostcellsCo+&
                  (interpolation_order-1)
             ixR_p_stg_min1(idir,2)=ixCoMmin1-nghostcellsCo-(interpolation_order-&
                  1)
             ixR_p_stg_max1(idir,2)=ixCoMmax1
             ixR_p_stg_min1(idir,3)=ixCoMmax1+1
             ixR_p_stg_max1(idir,3)=ixCoMmax1+nghostcellsCo+&
                  (interpolation_order-1)


             ! Perpendicular component
             ixS_p_stg_min1(idir,0)=ixMmin1
             ixS_p_stg_max1(idir,0)=ixMmin1-1+nghostcellsCo+&
                  (interpolation_order-1)
             ixS_p_stg_min1(idir,1)=ixMmin1
             ixS_p_stg_max1(idir,1)=ixMmin1-1+nxCo1+nghostcellsCo+&
                  (interpolation_order-1)
             ixS_p_stg_min1(idir,2)=ixMmax1+&
                  1-nxCo1-nghostcellsCo-(interpolation_order-1)
             ixS_p_stg_max1(idir,2)=ixMmax1
             ixS_p_stg_min1(idir,3)=ixMmax1+&
                  1-nghostcellsCo-(interpolation_order-1)
             ixS_p_stg_max1(idir,3)=ixMmax1

             ixR_p_stg_min1(idir,0)=ixCoMmin1-nghostcellsCo-(interpolation_order-&
                  1)
             ixR_p_stg_max1(idir,0)=ixCoMmin1-1
             ixR_p_stg_min1(idir,1)=ixCoMmin1
             ixR_p_stg_max1(idir,1)=ixCoMmax1+nghostcellsCo+&
                  (interpolation_order-1)
             ixR_p_stg_min1(idir,2)=ixCoMmin1-nghostcellsCo-(interpolation_order-&
                  1)
             ixR_p_stg_max1(idir,2)=ixCoMmax1
             ixR_p_stg_min1(idir,3)=ixCoMmax1+1
             ixR_p_stg_max1(idir,3)=ixCoMmax1+nghostcellsCo+&
                  (interpolation_order-1)


             ! Perpendicular component
             ixS_p_stg_min1(idir,0)=ixMmin1
             ixS_p_stg_max1(idir,0)=ixMmin1-1+nghostcellsCo+&
                  (interpolation_order-1)
             ixS_p_stg_min1(idir,1)=ixMmin1
             ixS_p_stg_max1(idir,1)=ixMmin1-1+nxCo1+nghostcellsCo+&
                  (interpolation_order-1)
             ixS_p_stg_min1(idir,2)=ixMmax1+&
                  1-nxCo1-nghostcellsCo-(interpolation_order-1)
             ixS_p_stg_max1(idir,2)=ixMmax1
             ixS_p_stg_min1(idir,3)=ixMmax1+&
                  1-nghostcellsCo-(interpolation_order-1)
             ixS_p_stg_max1(idir,3)=ixMmax1

             ixR_p_stg_min1(idir,0)=ixCoMmin1-nghostcellsCo-(interpolation_order-&
                  1)
             ixR_p_stg_max1(idir,0)=ixCoMmin1-1
             ixR_p_stg_min1(idir,1)=ixCoMmin1
             ixR_p_stg_max1(idir,1)=ixCoMmax1+nghostcellsCo+&
                  (interpolation_order-1)
             ixR_p_stg_min1(idir,2)=ixCoMmin1-nghostcellsCo-(interpolation_order-&
                  1)
             ixR_p_stg_max1(idir,2)=ixCoMmax1
             ixR_p_stg_min1(idir,3)=ixCoMmax1+1
             ixR_p_stg_max1(idir,3)=ixCoMmax1+nghostcellsCo+&
                  (interpolation_order-1)

          end if
          if (idir==2) then
             ! Parallel components

             ixS_p_stg_min2(idir,0)=ixMmin2-1 ! -1 to make redundant
             ixS_p_stg_max2(idir,0)=ixMmin2-1+nghostcellsCo
             ixS_p_stg_min2(idir,1)=ixMmin2-1 ! -1 to make redundant
             ixS_p_stg_max2(idir,1)=ixMmin2-1+nxCo2+nghostcellsCo
             ixS_p_stg_min2(idir,2)=ixMmax2-nxCo2-nghostcellsCo
             ixS_p_stg_max2(idir,2)=ixMmax2
             ixS_p_stg_min2(idir,3)=ixMmax2-nghostcellsCo
             ixS_p_stg_max2(idir,3)=ixMmax2

             ixR_p_stg_min2(idir,0)=ixCoMmin2-1-nghostcellsCo
             ixR_p_stg_max2(idir,0)=ixCoMmin2-1
             ixR_p_stg_min2(idir,1)=ixCoMmin2-1 ! -1 to make redundant
             ixR_p_stg_max2(idir,1)=ixCoMmax2+nghostcellsCo
             ixR_p_stg_min2(idir,2)=ixCoMmin2-1-nghostcellsCo
             ixR_p_stg_max2(idir,2)=ixCoMmax2
             ixR_p_stg_min2(idir,3)=ixCoMmax2+1-1 ! -1 to make redundant
             ixR_p_stg_max2(idir,3)=ixCoMmax2+nghostcellsCo


             ixS_p_stg_min2(idir,0)=ixMmin2-1 ! -1 to make redundant
             ixS_p_stg_max2(idir,0)=ixMmin2-1+nghostcellsCo
             ixS_p_stg_min2(idir,1)=ixMmin2-1 ! -1 to make redundant
             ixS_p_stg_max2(idir,1)=ixMmin2-1+nxCo2+nghostcellsCo
             ixS_p_stg_min2(idir,2)=ixMmax2-nxCo2-nghostcellsCo
             ixS_p_stg_max2(idir,2)=ixMmax2
             ixS_p_stg_min2(idir,3)=ixMmax2-nghostcellsCo
             ixS_p_stg_max2(idir,3)=ixMmax2

             ixR_p_stg_min2(idir,0)=ixCoMmin2-1-nghostcellsCo
             ixR_p_stg_max2(idir,0)=ixCoMmin2-1
             ixR_p_stg_min2(idir,1)=ixCoMmin2-1 ! -1 to make redundant
             ixR_p_stg_max2(idir,1)=ixCoMmax2+nghostcellsCo
             ixR_p_stg_min2(idir,2)=ixCoMmin2-1-nghostcellsCo
             ixR_p_stg_max2(idir,2)=ixCoMmax2
             ixR_p_stg_min2(idir,3)=ixCoMmax2+1-1 ! -1 to make redundant
             ixR_p_stg_max2(idir,3)=ixCoMmax2+nghostcellsCo


             ixS_p_stg_min2(idir,0)=ixMmin2-1 ! -1 to make redundant
             ixS_p_stg_max2(idir,0)=ixMmin2-1+nghostcellsCo
             ixS_p_stg_min2(idir,1)=ixMmin2-1 ! -1 to make redundant
             ixS_p_stg_max2(idir,1)=ixMmin2-1+nxCo2+nghostcellsCo
             ixS_p_stg_min2(idir,2)=ixMmax2-nxCo2-nghostcellsCo
             ixS_p_stg_max2(idir,2)=ixMmax2
             ixS_p_stg_min2(idir,3)=ixMmax2-nghostcellsCo
             ixS_p_stg_max2(idir,3)=ixMmax2

             ixR_p_stg_min2(idir,0)=ixCoMmin2-1-nghostcellsCo
             ixR_p_stg_max2(idir,0)=ixCoMmin2-1
             ixR_p_stg_min2(idir,1)=ixCoMmin2-1 ! -1 to make redundant
             ixR_p_stg_max2(idir,1)=ixCoMmax2+nghostcellsCo
             ixR_p_stg_min2(idir,2)=ixCoMmin2-1-nghostcellsCo
             ixR_p_stg_max2(idir,2)=ixCoMmax2
             ixR_p_stg_min2(idir,3)=ixCoMmax2+1-1 ! -1 to make redundant
             ixR_p_stg_max2(idir,3)=ixCoMmax2+nghostcellsCo

          else

             ! Perpendicular component
             ixS_p_stg_min2(idir,0)=ixMmin2
             ixS_p_stg_max2(idir,0)=ixMmin2-1+nghostcellsCo+&
                  (interpolation_order-1)
             ixS_p_stg_min2(idir,1)=ixMmin2
             ixS_p_stg_max2(idir,1)=ixMmin2-1+nxCo2+nghostcellsCo+&
                  (interpolation_order-1)
             ixS_p_stg_min2(idir,2)=ixMmax2+&
                  1-nxCo2-nghostcellsCo-(interpolation_order-1)
             ixS_p_stg_max2(idir,2)=ixMmax2
             ixS_p_stg_min2(idir,3)=ixMmax2+&
                  1-nghostcellsCo-(interpolation_order-1)
             ixS_p_stg_max2(idir,3)=ixMmax2

             ixR_p_stg_min2(idir,0)=ixCoMmin2-nghostcellsCo-(interpolation_order-&
                  1)
             ixR_p_stg_max2(idir,0)=ixCoMmin2-1
             ixR_p_stg_min2(idir,1)=ixCoMmin2
             ixR_p_stg_max2(idir,1)=ixCoMmax2+nghostcellsCo+&
                  (interpolation_order-1)
             ixR_p_stg_min2(idir,2)=ixCoMmin2-nghostcellsCo-(interpolation_order-&
                  1)
             ixR_p_stg_max2(idir,2)=ixCoMmax2
             ixR_p_stg_min2(idir,3)=ixCoMmax2+1
             ixR_p_stg_max2(idir,3)=ixCoMmax2+nghostcellsCo+&
                  (interpolation_order-1)


             ! Perpendicular component
             ixS_p_stg_min2(idir,0)=ixMmin2
             ixS_p_stg_max2(idir,0)=ixMmin2-1+nghostcellsCo+&
                  (interpolation_order-1)
             ixS_p_stg_min2(idir,1)=ixMmin2
             ixS_p_stg_max2(idir,1)=ixMmin2-1+nxCo2+nghostcellsCo+&
                  (interpolation_order-1)
             ixS_p_stg_min2(idir,2)=ixMmax2+&
                  1-nxCo2-nghostcellsCo-(interpolation_order-1)
             ixS_p_stg_max2(idir,2)=ixMmax2
             ixS_p_stg_min2(idir,3)=ixMmax2+&
                  1-nghostcellsCo-(interpolation_order-1)
             ixS_p_stg_max2(idir,3)=ixMmax2

             ixR_p_stg_min2(idir,0)=ixCoMmin2-nghostcellsCo-(interpolation_order-&
                  1)
             ixR_p_stg_max2(idir,0)=ixCoMmin2-1
             ixR_p_stg_min2(idir,1)=ixCoMmin2
             ixR_p_stg_max2(idir,1)=ixCoMmax2+nghostcellsCo+&
                  (interpolation_order-1)
             ixR_p_stg_min2(idir,2)=ixCoMmin2-nghostcellsCo-(interpolation_order-&
                  1)
             ixR_p_stg_max2(idir,2)=ixCoMmax2
             ixR_p_stg_min2(idir,3)=ixCoMmax2+1
             ixR_p_stg_max2(idir,3)=ixCoMmax2+nghostcellsCo+&
                  (interpolation_order-1)


             ! Perpendicular component
             ixS_p_stg_min2(idir,0)=ixMmin2
             ixS_p_stg_max2(idir,0)=ixMmin2-1+nghostcellsCo+&
                  (interpolation_order-1)
             ixS_p_stg_min2(idir,1)=ixMmin2
             ixS_p_stg_max2(idir,1)=ixMmin2-1+nxCo2+nghostcellsCo+&
                  (interpolation_order-1)
             ixS_p_stg_min2(idir,2)=ixMmax2+&
                  1-nxCo2-nghostcellsCo-(interpolation_order-1)
             ixS_p_stg_max2(idir,2)=ixMmax2
             ixS_p_stg_min2(idir,3)=ixMmax2+&
                  1-nghostcellsCo-(interpolation_order-1)
             ixS_p_stg_max2(idir,3)=ixMmax2

             ixR_p_stg_min2(idir,0)=ixCoMmin2-nghostcellsCo-(interpolation_order-&
                  1)
             ixR_p_stg_max2(idir,0)=ixCoMmin2-1
             ixR_p_stg_min2(idir,1)=ixCoMmin2
             ixR_p_stg_max2(idir,1)=ixCoMmax2+nghostcellsCo+&
                  (interpolation_order-1)
             ixR_p_stg_min2(idir,2)=ixCoMmin2-nghostcellsCo-(interpolation_order-&
                  1)
             ixR_p_stg_max2(idir,2)=ixCoMmax2
             ixR_p_stg_min2(idir,3)=ixCoMmax2+1
             ixR_p_stg_max2(idir,3)=ixCoMmax2+nghostcellsCo+&
                  (interpolation_order-1)

          end if
          if (idir==3) then
             ! Parallel components

             ixS_p_stg_min3(idir,0)=ixMmin3-1 ! -1 to make redundant
             ixS_p_stg_max3(idir,0)=ixMmin3-1+nghostcellsCo
             ixS_p_stg_min3(idir,1)=ixMmin3-1 ! -1 to make redundant
             ixS_p_stg_max3(idir,1)=ixMmin3-1+nxCo3+nghostcellsCo
             ixS_p_stg_min3(idir,2)=ixMmax3-nxCo3-nghostcellsCo
             ixS_p_stg_max3(idir,2)=ixMmax3
             ixS_p_stg_min3(idir,3)=ixMmax3-nghostcellsCo
             ixS_p_stg_max3(idir,3)=ixMmax3

             ixR_p_stg_min3(idir,0)=ixCoMmin3-1-nghostcellsCo
             ixR_p_stg_max3(idir,0)=ixCoMmin3-1
             ixR_p_stg_min3(idir,1)=ixCoMmin3-1 ! -1 to make redundant
             ixR_p_stg_max3(idir,1)=ixCoMmax3+nghostcellsCo
             ixR_p_stg_min3(idir,2)=ixCoMmin3-1-nghostcellsCo
             ixR_p_stg_max3(idir,2)=ixCoMmax3
             ixR_p_stg_min3(idir,3)=ixCoMmax3+1-1 ! -1 to make redundant
             ixR_p_stg_max3(idir,3)=ixCoMmax3+nghostcellsCo


             ixS_p_stg_min3(idir,0)=ixMmin3-1 ! -1 to make redundant
             ixS_p_stg_max3(idir,0)=ixMmin3-1+nghostcellsCo
             ixS_p_stg_min3(idir,1)=ixMmin3-1 ! -1 to make redundant
             ixS_p_stg_max3(idir,1)=ixMmin3-1+nxCo3+nghostcellsCo
             ixS_p_stg_min3(idir,2)=ixMmax3-nxCo3-nghostcellsCo
             ixS_p_stg_max3(idir,2)=ixMmax3
             ixS_p_stg_min3(idir,3)=ixMmax3-nghostcellsCo
             ixS_p_stg_max3(idir,3)=ixMmax3

             ixR_p_stg_min3(idir,0)=ixCoMmin3-1-nghostcellsCo
             ixR_p_stg_max3(idir,0)=ixCoMmin3-1
             ixR_p_stg_min3(idir,1)=ixCoMmin3-1 ! -1 to make redundant
             ixR_p_stg_max3(idir,1)=ixCoMmax3+nghostcellsCo
             ixR_p_stg_min3(idir,2)=ixCoMmin3-1-nghostcellsCo
             ixR_p_stg_max3(idir,2)=ixCoMmax3
             ixR_p_stg_min3(idir,3)=ixCoMmax3+1-1 ! -1 to make redundant
             ixR_p_stg_max3(idir,3)=ixCoMmax3+nghostcellsCo


             ixS_p_stg_min3(idir,0)=ixMmin3-1 ! -1 to make redundant
             ixS_p_stg_max3(idir,0)=ixMmin3-1+nghostcellsCo
             ixS_p_stg_min3(idir,1)=ixMmin3-1 ! -1 to make redundant
             ixS_p_stg_max3(idir,1)=ixMmin3-1+nxCo3+nghostcellsCo
             ixS_p_stg_min3(idir,2)=ixMmax3-nxCo3-nghostcellsCo
             ixS_p_stg_max3(idir,2)=ixMmax3
             ixS_p_stg_min3(idir,3)=ixMmax3-nghostcellsCo
             ixS_p_stg_max3(idir,3)=ixMmax3

             ixR_p_stg_min3(idir,0)=ixCoMmin3-1-nghostcellsCo
             ixR_p_stg_max3(idir,0)=ixCoMmin3-1
             ixR_p_stg_min3(idir,1)=ixCoMmin3-1 ! -1 to make redundant
             ixR_p_stg_max3(idir,1)=ixCoMmax3+nghostcellsCo
             ixR_p_stg_min3(idir,2)=ixCoMmin3-1-nghostcellsCo
             ixR_p_stg_max3(idir,2)=ixCoMmax3
             ixR_p_stg_min3(idir,3)=ixCoMmax3+1-1 ! -1 to make redundant
             ixR_p_stg_max3(idir,3)=ixCoMmax3+nghostcellsCo

          else

             ! Perpendicular component
             ixS_p_stg_min3(idir,0)=ixMmin3
             ixS_p_stg_max3(idir,0)=ixMmin3-1+nghostcellsCo+&
                  (interpolation_order-1)
             ixS_p_stg_min3(idir,1)=ixMmin3
             ixS_p_stg_max3(idir,1)=ixMmin3-1+nxCo3+nghostcellsCo+&
                  (interpolation_order-1)
             ixS_p_stg_min3(idir,2)=ixMmax3+&
                  1-nxCo3-nghostcellsCo-(interpolation_order-1)
             ixS_p_stg_max3(idir,2)=ixMmax3
             ixS_p_stg_min3(idir,3)=ixMmax3+&
                  1-nghostcellsCo-(interpolation_order-1)
             ixS_p_stg_max3(idir,3)=ixMmax3

             ixR_p_stg_min3(idir,0)=ixCoMmin3-nghostcellsCo-(interpolation_order-&
                  1)
             ixR_p_stg_max3(idir,0)=ixCoMmin3-1
             ixR_p_stg_min3(idir,1)=ixCoMmin3
             ixR_p_stg_max3(idir,1)=ixCoMmax3+nghostcellsCo+&
                  (interpolation_order-1)
             ixR_p_stg_min3(idir,2)=ixCoMmin3-nghostcellsCo-(interpolation_order-&
                  1)
             ixR_p_stg_max3(idir,2)=ixCoMmax3
             ixR_p_stg_min3(idir,3)=ixCoMmax3+1
             ixR_p_stg_max3(idir,3)=ixCoMmax3+nghostcellsCo+&
                  (interpolation_order-1)


             ! Perpendicular component
             ixS_p_stg_min3(idir,0)=ixMmin3
             ixS_p_stg_max3(idir,0)=ixMmin3-1+nghostcellsCo+&
                  (interpolation_order-1)
             ixS_p_stg_min3(idir,1)=ixMmin3
             ixS_p_stg_max3(idir,1)=ixMmin3-1+nxCo3+nghostcellsCo+&
                  (interpolation_order-1)
             ixS_p_stg_min3(idir,2)=ixMmax3+&
                  1-nxCo3-nghostcellsCo-(interpolation_order-1)
             ixS_p_stg_max3(idir,2)=ixMmax3
             ixS_p_stg_min3(idir,3)=ixMmax3+&
                  1-nghostcellsCo-(interpolation_order-1)
             ixS_p_stg_max3(idir,3)=ixMmax3

             ixR_p_stg_min3(idir,0)=ixCoMmin3-nghostcellsCo-(interpolation_order-&
                  1)
             ixR_p_stg_max3(idir,0)=ixCoMmin3-1
             ixR_p_stg_min3(idir,1)=ixCoMmin3
             ixR_p_stg_max3(idir,1)=ixCoMmax3+nghostcellsCo+&
                  (interpolation_order-1)
             ixR_p_stg_min3(idir,2)=ixCoMmin3-nghostcellsCo-(interpolation_order-&
                  1)
             ixR_p_stg_max3(idir,2)=ixCoMmax3
             ixR_p_stg_min3(idir,3)=ixCoMmax3+1
             ixR_p_stg_max3(idir,3)=ixCoMmax3+nghostcellsCo+&
                  (interpolation_order-1)


             ! Perpendicular component
             ixS_p_stg_min3(idir,0)=ixMmin3
             ixS_p_stg_max3(idir,0)=ixMmin3-1+nghostcellsCo+&
                  (interpolation_order-1)
             ixS_p_stg_min3(idir,1)=ixMmin3
             ixS_p_stg_max3(idir,1)=ixMmin3-1+nxCo3+nghostcellsCo+&
                  (interpolation_order-1)
             ixS_p_stg_min3(idir,2)=ixMmax3+&
                  1-nxCo3-nghostcellsCo-(interpolation_order-1)
             ixS_p_stg_max3(idir,2)=ixMmax3
             ixS_p_stg_min3(idir,3)=ixMmax3+&
                  1-nghostcellsCo-(interpolation_order-1)
             ixS_p_stg_max3(idir,3)=ixMmax3

             ixR_p_stg_min3(idir,0)=ixCoMmin3-nghostcellsCo-(interpolation_order-&
                  1)
             ixR_p_stg_max3(idir,0)=ixCoMmin3-1
             ixR_p_stg_min3(idir,1)=ixCoMmin3
             ixR_p_stg_max3(idir,1)=ixCoMmax3+nghostcellsCo+&
                  (interpolation_order-1)
             ixR_p_stg_min3(idir,2)=ixCoMmin3-nghostcellsCo-(interpolation_order-&
                  1)
             ixR_p_stg_max3(idir,2)=ixCoMmax3
             ixR_p_stg_min3(idir,3)=ixCoMmax3+1
             ixR_p_stg_max3(idir,3)=ixCoMmax3+nghostcellsCo+&
                  (interpolation_order-1)

          end if

       end do
       ! calculate sizes for buffer arrays for siblings
       do i3=-1,1
          do i2=-1,1
             do i1=-1,1
                ! Staggered (face-allocated) variables
                do idir=1,ndim
                   sizes_srl_send_stg(idir,i1,i2,i3)=(ixS_srl_stg_max1(idir,&
                        i1)-ixS_srl_stg_min1(idir,i1)+1)*(ixS_srl_stg_max2(idir,&
                        i2)-ixS_srl_stg_min2(idir,i2)+1)*(ixS_srl_stg_max3(idir,&
                        i3)-ixS_srl_stg_min3(idir,i3)+1)
                   sizes_srl_recv_stg(idir,i1,i2,i3)=(ixR_srl_stg_max1(idir,&
                        i1)-ixR_srl_stg_min1(idir,i1)+1)*(ixR_srl_stg_max2(idir,&
                        i2)-ixR_srl_stg_min2(idir,i2)+1)*(ixR_srl_stg_max3(idir,&
                        i3)-ixR_srl_stg_min3(idir,i3)+1)
                   sizes_r_send_stg(idir,i1,i2,i3)=(ixS_r_stg_max1(idir,&
                        i1)-ixS_r_stg_min1(idir,i1)+1)*(ixS_r_stg_max2(idir,&
                        i2)-ixS_r_stg_min2(idir,i2)+1)*(ixS_r_stg_max3(idir,&
                        i3)-ixS_r_stg_min3(idir,i3)+1)
                end do
                sizes_srl_send_total(i1,i2,i3)=sum(sizes_srl_send_stg(:,i1,i2,i3))
                sizes_srl_recv_total(i1,i2,i3)=sum(sizes_srl_recv_stg(:,i1,i2,i3))
                sizes_r_send_total(i1,i2,i3)=sum(sizes_r_send_stg(:,i1,i2,i3))
             end do
          end do
       end do

       do i3=0,3
          do i2=0,3
             do i1=0,3
                ! Staggered (face-allocated) variables
                do idir=1,ndim
                   sizes_r_recv_stg(idir,i1,i2,i3)=(ixR_r_stg_max1(idir,&
                        i1)-ixR_r_stg_min1(idir,i1)+1)*(ixR_r_stg_max2(idir,&
                        i2)-ixR_r_stg_min2(idir,i2)+1)*(ixR_r_stg_max3(idir,&
                        i3)-ixR_r_stg_min3(idir,i3)+1)
                   sizes_p_send_stg(idir,i1,i2,i3)=(ixS_p_stg_max1(idir,&
                        i1)-ixS_p_stg_min1(idir,i1)+1)*(ixS_p_stg_max2(idir,&
                        i2)-ixS_p_stg_min2(idir,i2)+1)*(ixS_p_stg_max3(idir,&
                        i3)-ixS_p_stg_min3(idir,i3)+1)
                   sizes_p_recv_stg(idir,i1,i2,i3)=(ixR_p_stg_max1(idir,&
                        i1)-ixR_p_stg_min1(idir,i1)+1)*(ixR_p_stg_max2(idir,&
                        i2)-ixR_p_stg_min2(idir,i2)+1)*(ixR_p_stg_max3(idir,&
                        i3)-ixR_p_stg_min3(idir,i3)+1)
                end do
                sizes_r_recv_total(i1,i2,i3)=sum(sizes_r_recv_stg(:,i1,i2,i3))
                sizes_p_send_total(i1,i2,i3)=sum(sizes_p_send_stg(:,i1,i2,i3))
                sizes_p_recv_total(i1,i2,i3)=sum(sizes_p_recv_stg(:,i1,i2,i3))
             end do
          end do
       end do
    end if
    if(.not.stagger_grid .or. physics_type=='mf') then
       ! extend index range to physical boundary

       ixS_srl_min1(-1,0)=1
       ixS_srl_min1( 1,0)=ixMmin1
       ixS_srl_min1( 2,0)=1
       ixS_srl_max1(-1,0)=ixMmax1
       ixS_srl_max1( 1,0)=ixGmax1
       ixS_srl_max1( 2,0)=ixGmax1

       ixR_srl_min1(-1,0)=1
       ixR_srl_min1( 1,0)=ixMmin1
       ixR_srl_min1( 2,0)=1
       ixR_srl_max1(-1,0)=ixMmax1
       ixR_srl_max1( 1,0)=ixGmax1
       ixR_srl_max1( 2,0)=ixGmax1

       ixS_r_min1(-1,0)=1
       ixS_r_min1( 1,0)=ixCoMmin1
       ixS_r_max1(-1,0)=ixCoMmax1
       ixS_r_max1( 1,0)=ixCoGmax1

       ixR_r_min1(-1,1)=1
       ixR_r_max1(-1,1)=ixMmin1-1+nxCo1
       ixR_r_min1( 1,2)=ixMmin1+nxCo1
       ixR_r_max1( 1,2)=ixGmax1

       ixS_p_min1(-1,1)=1
       ixS_p_max1( 1,2)=ixGmax1

       ixR_p_min1(-1,1)=1
       ixR_p_max1( 1,2)=ixCoGmax1


       ixS_srl_min2(-1,0)=1
       ixS_srl_min2( 1,0)=ixMmin2
       ixS_srl_min2( 2,0)=1
       ixS_srl_max2(-1,0)=ixMmax2
       ixS_srl_max2( 1,0)=ixGmax2
       ixS_srl_max2( 2,0)=ixGmax2

       ixR_srl_min2(-1,0)=1
       ixR_srl_min2( 1,0)=ixMmin2
       ixR_srl_min2( 2,0)=1
       ixR_srl_max2(-1,0)=ixMmax2
       ixR_srl_max2( 1,0)=ixGmax2
       ixR_srl_max2( 2,0)=ixGmax2

       ixS_r_min2(-1,0)=1
       ixS_r_min2( 1,0)=ixCoMmin2
       ixS_r_max2(-1,0)=ixCoMmax2
       ixS_r_max2( 1,0)=ixCoGmax2

       ixR_r_min2(-1,1)=1
       ixR_r_max2(-1,1)=ixMmin2-1+nxCo2
       ixR_r_min2( 1,2)=ixMmin2+nxCo2
       ixR_r_max2( 1,2)=ixGmax2

       ixS_p_min2(-1,1)=1
       ixS_p_max2( 1,2)=ixGmax2

       ixR_p_min2(-1,1)=1
       ixR_p_max2( 1,2)=ixCoGmax2


       ixS_srl_min3(-1,0)=1
       ixS_srl_min3( 1,0)=ixMmin3
       ixS_srl_min3( 2,0)=1
       ixS_srl_max3(-1,0)=ixMmax3
       ixS_srl_max3( 1,0)=ixGmax3
       ixS_srl_max3( 2,0)=ixGmax3

       ixR_srl_min3(-1,0)=1
       ixR_srl_min3( 1,0)=ixMmin3
       ixR_srl_min3( 2,0)=1
       ixR_srl_max3(-1,0)=ixMmax3
       ixR_srl_max3( 1,0)=ixGmax3
       ixR_srl_max3( 2,0)=ixGmax3

       ixS_r_min3(-1,0)=1
       ixS_r_min3( 1,0)=ixCoMmin3
       ixS_r_max3(-1,0)=ixCoMmax3
       ixS_r_max3( 1,0)=ixCoGmax3

       ixR_r_min3(-1,1)=1
       ixR_r_max3(-1,1)=ixMmin3-1+nxCo3
       ixR_r_min3( 1,2)=ixMmin3+nxCo3
       ixR_r_max3( 1,2)=ixGmax3

       ixS_p_min3(-1,1)=1
       ixS_p_max3( 1,2)=ixGmax3

       ixR_p_min3(-1,1)=1
       ixR_p_max3( 1,2)=ixCoGmax3

    end if

    !$acc update device(ixS_srl_min1,ixS_srl_min2,ixS_srl_min3)
    !$acc update device(ixS_srl_max1,ixS_srl_max2,ixS_srl_max3, ixR_srl_min1,ixR_srl_min2)
    !$acc update device(ixR_srl_min3,ixR_srl_max1,ixR_srl_max2,ixR_srl_max3)

  end subroutine init_bc

  subroutine create_bc_mpi_datatype(nwstart,nwbc)
    use mod_global_parameters

    integer, intent(in) :: nwstart, nwbc
    integer :: i1,i2,i3, ic1,ic2,ic3, inc1,inc2,inc3, iib1,iib2,iib3

    do i3=-1,1
       do i2=-1,1
          do i1=-1,1
             if (i1==0.and.i2==0.and.i3==0) cycle
             do iib3=-1,2
                do iib2=-1,2
                   do iib1=-1,2
                      call get_bc_comm_type(type_send_srl(iib1,iib2,iib3,i1,i2,i3),&
                           ixS_srl_min1(iib1,i1),ixS_srl_min2(iib2,i2),ixS_srl_min3(iib3,i3),&
                           ixS_srl_max1(iib1,i1),ixS_srl_max2(iib2,i2),ixS_srl_max3(iib3,i3),&
                           ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,nwstart,nwbc)
                      call get_bc_comm_type(type_recv_srl(iib1,iib2,iib3,i1,i2,i3),&
                           ixR_srl_min1(iib1,i1),ixR_srl_min2(iib2,i2),ixR_srl_min3(iib3,i3),&
                           ixR_srl_max1(iib1,i1),ixR_srl_max2(iib2,i2),ixR_srl_max3(iib3,i3),&
                           ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,nwstart,nwbc)
                      if (iib1==2.or.iib2==2.or.iib3==2) cycle
                      call get_bc_comm_type(type_send_r(iib1,iib2,iib3,i1,i2,i3),&
                           ixS_r_min1(iib1,i1),ixS_r_min2(iib2,i2),ixS_r_min3(iib3,i3),&
                           ixS_r_max1(iib1,i1),ixS_r_max2(iib2,i2),ixS_r_max3(iib3,i3),&
                           ixCoGmin1,ixCoGmin2,ixCoGmin3,ixCoGmax1,ixCoGmax2,ixCoGmax3,&
                           nwstart,nwbc)
                      do ic3=1+int((1-i3)/2),2-int((1+i3)/2)
                         inc3=2*i3+ic3
                         do ic2=1+int((1-i2)/2),2-int((1+i2)/2)
                            inc2=2*i2+ic2
                            do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
                               inc1=2*i1+ic1
                               call get_bc_comm_type(type_recv_r(iib1,iib2,iib3,inc1,inc2,inc3),&
                                    ixR_r_min1(iib1,inc1),ixR_r_min2(iib2,inc2),ixR_r_min3(iib3,&
                                    inc3),ixR_r_max1(iib1,inc1),ixR_r_max2(iib2,inc2),&
                                    ixR_r_max3(iib3,inc3), ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,&
                                    ixGhi3,nwstart,nwbc)
                               call get_bc_comm_type(type_send_p(iib1,iib2,iib3,inc1,inc2,inc3),&
                                    ixS_p_min1(iib1,inc1),ixS_p_min2(iib2,inc2),ixS_p_min3(iib3,&
                                    inc3),ixS_p_max1(iib1,inc1),ixS_p_max2(iib2,inc2),&
                                    ixS_p_max3(iib3,inc3), ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,&
                                    ixGhi3,nwstart,nwbc)
                               call get_bc_comm_type(type_recv_p(iib1,iib2,iib3,inc1,inc2,inc3),&
                                    ixR_p_min1(iib1,inc1),ixR_p_min2(iib2,inc2),ixR_p_min3(iib3,&
                                    inc3),ixR_p_max1(iib1,inc1),ixR_p_max2(iib2,inc2),&
                                    ixR_p_max3(iib3,inc3),ixCoGmin1,ixCoGmin2,ixCoGmin3,ixCoGmax1,&
                                    ixCoGmax2,ixCoGmax3,nwstart,nwbc)
                            end do
                         end do
                      end do
                   end do
                end do
             end do
          end do
       end do
    end do

  end subroutine create_bc_mpi_datatype

  subroutine get_bc_comm_type(comm_type,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,&
       ixmax3,ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,nwstart,nwbc)
    use mod_global_parameters

    integer, intent(inout) :: comm_type
    integer, intent(in) :: ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3, ixGmin1,&
         ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3, nwstart, nwbc

    integer, dimension(ndim+1) :: fullsize, subsize, start

    fullsize(1)=ixGmax1;fullsize(2)=ixGmax2;fullsize(3)=ixGmax3;
    fullsize(ndim+1)=nw
    subsize(1)=ixmax1-ixmin1+1;subsize(2)=ixmax2-ixmin2+1
    subsize(3)=ixmax3-ixmin3+1;
    subsize(ndim+1)=nwbc
    start(1)=ixmin1-1;start(2)=ixmin2-1;start(3)=ixmin3-1;
    start(ndim+1)=nwstart-1

    call MPI_TYPE_CREATE_SUBARRAY(ndim+1,fullsize,subsize,start,&
         MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION,comm_type,ierrmpi)
    call MPI_TYPE_COMMIT(comm_type,ierrmpi)

  end subroutine get_bc_comm_type

  subroutine put_bc_comm_types()
    use mod_global_parameters

    integer :: i1,i2,i3, ic1,ic2,ic3, inc1,inc2,inc3, iib1,iib2,iib3

    do i3=-1,1
       do i2=-1,1
          do i1=-1,1
             if (i1==0.and.i2==0.and.i3==0) cycle
             do iib3=-1,2
                do iib2=-1,2
                   do iib1=-1,2
                      call MPI_TYPE_FREE(type_send_srl(iib1,iib2,iib3,i1,i2,i3),ierrmpi)
                      call MPI_TYPE_FREE(type_recv_srl(iib1,iib2,iib3,i1,i2,i3),ierrmpi)
                      if (levmin==levmax) cycle
                      if (iib1==2.or.iib2==2.or.iib3==2) cycle
                      call MPI_TYPE_FREE(type_send_r(iib1,iib2,iib3,i1,i2,i3),ierrmpi)
                      do ic3=1+int((1-i3)/2),2-int((1+i3)/2)
                         inc3=2*i3+ic3
                         do ic2=1+int((1-i2)/2),2-int((1+i2)/2)
                            inc2=2*i2+ic2
                            do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
                               inc1=2*i1+ic1
                               call MPI_TYPE_FREE(type_recv_r(iib1,iib2,iib3,inc1,inc2,inc3),&
                                    ierrmpi)
                               call MPI_TYPE_FREE(type_send_p(iib1,iib2,iib3,inc1,inc2,inc3),&
                                    ierrmpi)
                               call MPI_TYPE_FREE(type_recv_p(iib1,iib2,iib3,inc1,inc2,inc3),&
                                    ierrmpi)
                            end do
                         end do
                      end do
                   end do
                end do
             end do
          end do
       end do
    end do

  end subroutine put_bc_comm_types

  !> do update ghost cells of all blocks including physical boundaries
  subroutine getbc(time,qdt,psb,nwstart,nwbc,req_diag)
    use mod_nvtx
    use mod_global_parameters
    use mod_physics
    use mod_coarsen, only: coarsen_grid
    use mod_boundary_conditions, only: getintbc, bc_phys
    use mod_comm_lib, only: mpistop

    double precision, intent(in)      :: time, qdt
    type(state), target               :: psb(max_blocks)
    integer, intent(in)               :: nwstart ! Fill from nwstart
    integer, intent(in)               :: nwbc    ! Number of variables to fill
    logical, intent(in), optional     :: req_diag !If false, skip diagonal ghost cells

    double precision :: time_bcin
    integer :: ipole, nwhead, nwtail
    integer :: iigrid, igrid, ineighbor, ipe_neighbor, isizes
    integer :: ixRmin1,ixRmin2,ixRmin3,ixRmax1,ixRmax2,ixRmax3, ixSmin1,&
         ixSmin2,ixSmin3,ixSmax1,ixSmax2,ixSmax3
    integer :: i1,i2,i3, n_i1,n_i2,n_i3, ic1,ic2,ic3, inc1,inc2,inc3, n_inc1,&
         n_inc2,n_inc3, iib1,iib2,iib3, idir, istage
    ! store physical boundary indicating index
    integer :: idphyb(ndim,max_blocks)
    integer :: isend_buf(npwbuf), ipwbuf, nghostcellsco
    ! index pointer for buffer arrays as a start for a segment
    integer :: ibuf_start, ibuf_next
    ! shapes of reshape
    integer, dimension(1) :: shapes
    logical  :: req_diagonal, update
    type(wbuffer) :: pwbuf(npwbuf)

    integer :: ix1,ix2,ix3, iw

    time_bcin=MPI_WTIME()

    call nvtxStartRange("getbc",2)

    nwhead=nwstart
    nwtail=nwstart+nwbc-1

    req_diagonal = .true.
    if (present(req_diag)) req_diagonal = req_diag

    ! fill internal physical boundary
    if (internalboundary) then
       call getintbc(time,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3)
    end if

    ! fill physical-boundary ghost cells before internal ghost-cell values exchange
    if(bcphys.and. .not.stagger_grid) then
       !$OMP PARALLEL DO SCHEDULE(dynamic) PRIVATE(igrid)
       do iigrid=1,igridstail; igrid=igrids(iigrid);
          if(.not.phyboundblock(igrid)) cycle
          call fill_boundary_before_gc(igrid)
       end do
       !$OMP END PARALLEL DO
    end if

    ! prepare coarse values to send to coarser neighbors
    !$OMP PARALLEL DO SCHEDULE(dynamic) PRIVATE(igrid)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       if(any(neighbor_type(:,:,:,igrid)==neighbor_coarse)) then
          call coarsen_grid(psb(igrid),ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,&
               ixMmin1,ixMmin2,ixMmin3,ixMmax1,ixMmax2,ixMmax3,psc(igrid),&
               ixCoGmin1,ixCoGmin2,ixCoGmin3,ixCoGmax1,ixCoGmax2,ixCoGmax3,&
               ixCoMmin1,ixCoMmin2,ixCoMmin3,ixCoMmax1,ixCoMmax2,ixCoMmax3)
          do i3=-1,1
             do i2=-1,1
                do i1=-1,1
                   if(skip_direction([ i1,i2,i3 ])) cycle
                   if(neighbor_type(i1,i2,i3,igrid)==neighbor_coarse) call &
                        fill_coarse_boundary(igrid,i1,i2,i3)
                end do
             end do
          end do
       end if
    end do
    !$OMP END PARALLEL DO

    ! default : no singular axis
    ipole=0
    irecv_c=0
    isend_c=0
    isend_buf=0
    ipwbuf=1

    if(stagger_grid) then
       ibuf_recv_srl=1
       ibuf_recv_r=1
       ibuf_recv_p=1
       ibuf_send_srl=1
       ibuf_send_r=1
       ibuf_send_p=1
       irecv_srl=0
       irecv_r=0
       irecv_p=0
       isend_srl=0
       isend_r=0
       isend_p=0
    end if

    ! MPI receive ghost-cell values from sibling blocks and finer neighbors in different processors
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       call identifyphysbound(ps(igrid),iib1,iib2,iib3)
       idphyb(1,igrid)=iib1;idphyb(2,igrid)=iib2;idphyb(3,igrid)=iib3;
       do i3=-1,1
          do i2=-1,1
             do i1=-1,1
                if (skip_direction([ i1,i2,i3 ])) cycle
                select case (neighbor_type(i1,i2,i3,igrid))
                case (neighbor_sibling)
                   call bc_recv_srl
                case (neighbor_fine)
                   call bc_recv_restrict
                end select
             end do
          end do
       end do
    end do

#ifdef _OPENACC
    ! update blocks on host prior to MPI send
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       call identifyphysbound(ps(igrid),iib1,iib2,iib3)
       idphyb(1,igrid)=iib1;idphyb(2,igrid)=iib2;idphyb(3,igrid)=iib3;
       update = .false.
       do i3=-1,1
          do i2=-1,1
             do i1=-1,1
                if (skip_direction([ i1,i2,i3 ])) cycle
                if (neighbor(2,i1,i2,i3,igrid) /= mype) then
                   update = .true.
                end if
             end do
          end do
       end do
       if (update) then
          !$acc update host(psb(igrid)%w)
       end if
    end do
#endif
    
    ! MPI send ghost-cell values to sibling blocks and coarser neighbors in different processors
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       iib1=idphyb(1,igrid);iib2=idphyb(2,igrid);iib3=idphyb(3,igrid);
       do i3=-1,1
          do i2=-1,1
             do i1=-1,1
                if(skip_direction([ i1,i2,i3 ])) cycle
                select case (neighbor_type(i1,i2,i3,igrid))
                case (neighbor_sibling)
                   call bc_send_srl
                case (neighbor_coarse)
                   call bc_send_restrict
                end select
             end do
          end do
       end do
    end do

    call MPI_WAITALL(irecv_c,recvrequest_c_sr,recvstatus_c_sr,ierrmpi)
    call MPI_WAITALL(isend_c,sendrequest_c_sr,sendstatus_c_sr,ierrmpi)

#ifdef _OPENACC    
    ! update blocks on device after MPI comm    
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       call identifyphysbound(ps(igrid),iib1,iib2,iib3)
       idphyb(1,igrid)=iib1;idphyb(2,igrid)=iib2;idphyb(3,igrid)=iib3;
       update = .false.
       do i3=-1,1
          do i2=-1,1
             do i1=-1,1
                if (skip_direction([ i1,i2,i3 ])) cycle
                if (neighbor(2,i1,i2,i3,igrid) /= mype) then
                   update = .true.
                end if
             end do
          end do
       end do
       if (update) then
          !$acc update device(psb(igrid)%w)
       end if
    end do
#endif
    
    ! fill ghost-cell values of sibling blocks and coarser neighbors in the same processor
    !$OMP PARALLEL DO SCHEDULE(dynamic) PRIVATE(igrid,iib1,iib2,iib3)
     !$acc parallel loop default(present) copyin(idphyb,ixS_srl_min1,ixS_srl_min2,ixS_srl_min3,ixS_srl_max1,ixS_srl_max2,ixS_srl_max3,ixR_srl_min1,ixR_srl_min2,ixR_srl_min3,ixR_srl_max1,ixR_srl_max2,ixR_srl_max3) private(igrid,iib1,iib2,iib3,ineighbor,n_i1,n_i2,n_i3,ixSmin1,ixSmin2,ixSmin3,ixSmax1,ixSmax2,ixSmax3,ixRmin1,ixRmin2,ixRmin3,ixRmax1,ixRmax2,ixRmax3,iw,ix1,ix2,ix3) firstprivate(nwhead,nwtail)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       iib1=idphyb(1,igrid);iib2=idphyb(2,igrid);iib3=idphyb(3,igrid);
       !$acc loop seq
       do i3=-1,1
          do i2=-1,1
             do i1=-1,1
                ! next line is inlined version of skip_direction
                if ((all([ i1,i2,i3 ] == 0)) .or. (.not. req_diagonal .and. count([ &
                     i1,i2,i3 ] /= 0) > 1)) cycle
                select case (neighbor_type(i1,i2,i3,igrid))
                case(neighbor_sibling)
                   call bc_fill_srl(psb,igrid,nwhead,nwtail,i1,i2,i3,iib1,iib2,iib3)
                   !         case(neighbor_coarse)
                   ! call bc_fill_restrict(igrid,i^D,iib^D)
                end select
             end do
          end do
       end do
    end do
    !$OMP END PARALLEL DO
    
    if(stagger_grid) then
       call MPI_WAITALL(nrecv_bc_srl,recvrequest_srl,recvstatus_srl,ierrmpi)
       call MPI_WAITALL(nsend_bc_srl,sendrequest_srl,sendstatus_srl,ierrmpi)
       call MPI_WAITALL(nrecv_bc_r,recvrequest_r,recvstatus_r,ierrmpi)
       call MPI_WAITALL(nsend_bc_r,sendrequest_r,sendstatus_r,ierrmpi)
       ! unpack the received data from sibling blocks and finer neighbors to fill ghost-cell staggered values
       ibuf_recv_srl=1
       ibuf_recv_r=1
       do iigrid=1,igridstail; igrid=igrids(iigrid);
          iib1=idphyb(1,igrid);iib2=idphyb(2,igrid);iib3=idphyb(3,igrid);
          do i3=-1,1
             do i2=-1,1
                do i1=-1,1
                   if (skip_direction([ i1,i2,i3 ])) cycle
                   select case (neighbor_type(i1,i2,i3,igrid))
                   case (neighbor_sibling)
                      call bc_fill_srl_stg
                   case (neighbor_fine)
                      call bc_fill_restrict_stg
                   end select
                end do
             end do
          end do
       end do
    end if

    do ipwbuf=1,npwbuf
       if (isend_buf(ipwbuf)/=0) deallocate(pwbuf(ipwbuf)%w)
    end do

    irecv_c=0
    isend_c=0
    isend_buf=0
    ipwbuf=1

    ! MPI receive ghost-cell values from coarser neighbors in different processors
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       iib1=idphyb(1,igrid);iib2=idphyb(2,igrid);iib3=idphyb(3,igrid);
       do i3=-1,1
          do i2=-1,1
             do i1=-1,1
                if (skip_direction([ i1,i2,i3 ])) cycle
                if (neighbor_type(i1,i2,i3,igrid)==neighbor_coarse) call &
                     bc_recv_prolong
             end do
          end do
       end do
    end do
    ! MPI send ghost-cell values to finer neighbors in different processors
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       iib1=idphyb(1,igrid);iib2=idphyb(2,igrid);iib3=idphyb(3,igrid);
       do i3=-1,1
          do i2=-1,1
             do i1=-1,1
                if (skip_direction([ i1,i2,i3 ])) cycle
                if (neighbor_type(i1,i2,i3,igrid)==neighbor_fine) call &
                     bc_send_prolong
             end do
          end do
       end do
    end do

    ! fill coarse ghost-cell values of finer neighbors in the same processor
    !$OMP PARALLEL DO SCHEDULE(dynamic) PRIVATE(igrid,iib1,iib2,iib3)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       iib1=idphyb(1,igrid);iib2=idphyb(2,igrid);iib3=idphyb(3,igrid);
       do i3=-1,1
          do i2=-1,1
             do i1=-1,1
                if (skip_direction([ i1,i2,i3 ])) cycle
                if (neighbor_type(i1,i2,i3,igrid)==neighbor_fine) call &
                     bc_fill_prolong(igrid,i1,i2,i3,iib1,iib2,iib3)
             end do
          end do
       end do
    end do
    !$OMP END PARALLEL DO

    call MPI_WAITALL(irecv_c,recvrequest_c_p,recvstatus_c_p,ierrmpi)
    call MPI_WAITALL(isend_c,sendrequest_c_p,sendstatus_c_p,ierrmpi)

    if(stagger_grid) then
       call MPI_WAITALL(nrecv_bc_p,recvrequest_p,recvstatus_p,ierrmpi)
       call MPI_WAITALL(nsend_bc_p,sendrequest_p,sendstatus_p,ierrmpi)

       ! fill coarser representative ghost cells after receipt
       ibuf_recv_p=1
       do iigrid=1,igridstail; igrid=igrids(iigrid);
          iib1=idphyb(1,igrid);iib2=idphyb(2,igrid);iib3=idphyb(3,igrid);
          do i3=-1,1
             do i2=-1,1
                do i1=-1,1
                   if (skip_direction([ i1,i2,i3 ])) cycle
                   if(neighbor_type(i1,i2,i3,igrid)==neighbor_coarse) call &
                        bc_fill_prolong_stg
                end do
             end do
          end do
       end do
    end if
    ! do prolongation on the ghost-cell values based on the received coarse values from coarser neighbors
    !$OMP PARALLEL DO SCHEDULE(dynamic) PRIVATE(igrid)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       call gc_prolong(igrid)
    end do
    !$OMP END PARALLEL DO

    do ipwbuf=1,npwbuf
       if (isend_buf(ipwbuf)/=0) deallocate(pwbuf(ipwbuf)%w)
    end do

    ! fill physical boundary ghost cells after internal ghost-cell values exchange
    if(bcphys.and.stagger_grid) then
       !$OMP PARALLEL DO SCHEDULE(dynamic) PRIVATE(igrid)
       do iigrid=1,igridstail; igrid=igrids(iigrid);
          if(.not.phyboundblock(igrid)) cycle
          call fill_boundary_after_gc(igrid)
       end do
       !$OMP END PARALLEL DO
    end if

    ! modify normal component of magnetic field to fix divB=0
    if(bcphys.and.associated(phys_boundary_adjust)) then
       !$OMP PARALLEL DO SCHEDULE(dynamic) PRIVATE(igrid)
       do iigrid=1,igridstail; igrid=igrids(iigrid);
          if(.not.phyboundblock(igrid)) cycle
          call phys_boundary_adjust(igrid,psb)
       end do
       !$OMP END PARALLEL DO
    end if

    time_bc=time_bc+(MPI_WTIME()-time_bcin)
    call nvtxEndRange

  contains

    logical function skip_direction(dir)
      integer, intent(in) :: dir(3)

      if (all(dir == 0)) then
         skip_direction = .true.
      else if (.not. req_diagonal .and. count(dir /= 0) > 1) then
         skip_direction = .true.
      else
         skip_direction = .false.
      end if
    end function skip_direction

    !> Physical boundary conditions
    subroutine fill_boundary_before_gc(igrid)

      integer, intent(in) :: igrid

      integer :: idims,iside,i1,i2,i3,kmin1,kmin2,kmin3,kmax1,kmax2,kmax3,&
           ixBmin1,ixBmin2,ixBmin3,ixBmax1,ixBmax2,ixBmax3

      block=>psb(igrid)
        dxlevel(1)=rnode(rpdx1_,igrid);dxlevel(2)=rnode(rpdx2_,igrid)
        dxlevel(3)=rnode(rpdx3_,igrid);
        do idims=1,ndim
           ! to avoid using as yet unknown corner info in more than 1D, we
           ! fill only interior mesh ranges of the ghost cell ranges at first,
           ! and progressively enlarge the ranges to include corners later

           kmin1=merge(0, 1, idims==1)
           kmax1=merge(0, 1, idims==1)
           ixBmin1=ixGlo1+kmin1*nghostcells
           ixBmax1=ixGhi1-kmax1*nghostcells


           kmin2=merge(0, 1, idims==2)
           kmax2=merge(0, 1, idims==2)
           ixBmin2=ixGlo2+kmin2*nghostcells
           ixBmax2=ixGhi2-kmax2*nghostcells


           kmin3=merge(0, 1, idims==3)
           kmax3=merge(0, 1, idims==3)
           ixBmin3=ixGlo3+kmin3*nghostcells
           ixBmax3=ixGhi3-kmax3*nghostcells



           if(idims > 1 .and. neighbor_type(-1,0,0,&
                igrid)==neighbor_boundary) ixBmin1=ixGlo1
           if(idims > 1 .and. neighbor_type( 1,0,0,&
                igrid)==neighbor_boundary) ixBmax1=ixGhi1
           if(idims > 2 .and. neighbor_type(0,-1,0,&
                igrid)==neighbor_boundary) ixBmin2=ixGlo2
           if(idims > 2 .and. neighbor_type(0, 1,0,&
                igrid)==neighbor_boundary) ixBmax2=ixGhi2
           do iside=1,2
              i1=kr(1,idims)*(2*iside-3);i2=kr(2,idims)*(2*iside-3)
              i3=kr(3,idims)*(2*iside-3);
              if (aperiodB(idims)) then
                 if (neighbor_type(i1,i2,i3,&
                      igrid) /= neighbor_boundary .and. .not. &
                      psb(igrid)%is_physical_boundary(2*idims-2+iside)) cycle
              else
                 if (neighbor_type(i1,i2,i3,igrid) /= neighbor_boundary) cycle
              end if
              call bc_phys(iside,idims,time,qdt,psb(igrid),ixGlo1,ixGlo2,ixGlo3,&
                   ixGhi1,ixGhi2,ixGhi3,ixBmin1,ixBmin2,ixBmin3,ixBmax1,ixBmax2,&
                   ixBmax3)
           end do
        end do

      end subroutine fill_boundary_before_gc

      !> Physical boundary conditions
      subroutine fill_boundary_after_gc(igrid)

        integer, intent(in) :: igrid

        integer :: idims,iside,i1,i2,i3,kmin1,kmin2,kmin3,kmax1,kmax2,kmax3,&
             ixBmin1,ixBmin2,ixBmin3,ixBmax1,ixBmax2,ixBmax3

        block=>psb(igrid)
          dxlevel(1)=rnode(rpdx1_,igrid);dxlevel(2)=rnode(rpdx2_,igrid)
          dxlevel(3)=rnode(rpdx3_,igrid);
          do idims=1,ndim
             ! to avoid using as yet unknown corner info in more than 1D, we
             ! fill only interior mesh ranges of the ghost cell ranges at first,
             ! and progressively enlarge the ranges to include corners later
             kmin1=0; kmax1=0;


             kmin2=merge(1, 0, idims .lt. 2 .and. neighbor_type(0,-1,0,&
                  igrid)==1)
             kmax2=merge(1, 0, idims .lt. 2 .and. neighbor_type(0, 1,0,&
                  igrid)==1)
             kmin3=merge(1, 0, idims .lt. 3 .and. neighbor_type(0,0,-1,&
                  igrid)==1)
             kmax3=merge(1, 0, idims .lt. 3 .and. neighbor_type(0,0, 1,&
                  igrid)==1)
             ixBmin1=ixGlo1+kmin1*nghostcells;ixBmin2=ixGlo2+kmin2*nghostcells
             ixBmin3=ixGlo3+kmin3*nghostcells;
             ixBmax1=ixGhi1-kmax1*nghostcells;ixBmax2=ixGhi2-kmax2*nghostcells
             ixBmax3=ixGhi3-kmax3*nghostcells;
             do iside=1,2
                i1=kr(1,idims)*(2*iside-3);i2=kr(2,idims)*(2*iside-3)
                i3=kr(3,idims)*(2*iside-3);
                if (aperiodB(idims)) then
                   if (neighbor_type(i1,i2,i3,&
                        igrid) /= neighbor_boundary .and. .not. &
                        psb(igrid)%is_physical_boundary(2*idims-2+iside)) cycle
                else
                   if (neighbor_type(i1,i2,i3,igrid) /= neighbor_boundary) cycle
                end if
                call bc_phys(iside,idims,time,qdt,psb(igrid),ixGlo1,ixGlo2,ixGlo3,&
                     ixGhi1,ixGhi2,ixGhi3,ixBmin1,ixBmin2,ixBmin3,ixBmax1,ixBmax2,&
                     ixBmax3)
             end do
          end do

        end subroutine fill_boundary_after_gc

        !> Receive from sibling at same refinement level
        subroutine bc_recv_srl
          integer                     :: istep

          ipe_neighbor=neighbor(2,i1,i2,i3,igrid)
          if (ipe_neighbor/=mype) then
             irecv_c=irecv_c+1
             itag=(3**3+4**3)*(igrid-1)+(i1+1)*3**(1-1)+(i2+1)*3**(2-1)+(i3+&
                  1)*3**(3-1)
             istep = psb(igrid)%istep
!             !$acc host_data use_device(bg(istep)%w)
             call MPI_IRECV(bg(istep)%w(:,:,:,:,igrid),1,type_recv_srl(iib1,iib2,iib3,i1,i2,&
                  i3), ipe_neighbor,itag,icomm,recvrequest_c_sr(irecv_c),ierrmpi)
!             !$acc end host_data
             if(stagger_grid) then
                irecv_srl=irecv_srl+1
                call MPI_IRECV(recvbuffer_srl(ibuf_recv_srl),&
                     sizes_srl_recv_total(i1,i2,i3),MPI_DOUBLE_PRECISION,&
                     ipe_neighbor,itag,icomm,recvrequest_srl(irecv_srl),ierrmpi)
                ibuf_recv_srl=ibuf_recv_srl+sizes_srl_recv_total(i1,i2,i3)
             end if
          end if

        end subroutine bc_recv_srl

        !> Send to sibling at same refinement level
        subroutine bc_send_srl
          integer                     :: istep

          ipe_neighbor=neighbor(2,i1,i2,i3,igrid)

          if(ipe_neighbor/=mype) then
             ineighbor=neighbor(1,i1,i2,i3,igrid)
             ipole=neighbor_pole(i1,i2,i3,igrid)
             if(ipole==0) then
                n_i1=-i1;n_i2=-i2;n_i3=-i3;
                isend_c=isend_c+1
                itag=(3**3+4**3)*(ineighbor-1)+(n_i1+1)*3**(1-1)+(n_i2+1)*3**(2-1)+&
                     (n_i3+1)*3**(3-1)
                istep = psb(igrid)%istep
!                !$acc host_data use_device(bg(istep)%w)
                call MPI_ISEND(bg(istep)%w(:,:,:,:,igrid),1,type_send_srl(iib1,iib2,iib3,i1,i2,&
                     i3), ipe_neighbor,itag,icomm,sendrequest_c_sr(isend_c),ierrmpi)
!                !$acc end host_data
                if(stagger_grid) then
                   ibuf_start=ibuf_send_srl
                   do idir=1,ndim
                      ixSmin1=ixS_srl_stg_min1(idir,i1)
                      ixSmin2=ixS_srl_stg_min2(idir,i2)
                      ixSmin3=ixS_srl_stg_min3(idir,i3)
                      ixSmax1=ixS_srl_stg_max1(idir,i1)
                      ixSmax2=ixS_srl_stg_max2(idir,i2)
                      ixSmax3=ixS_srl_stg_max3(idir,i3);
                      ibuf_next=ibuf_start+sizes_srl_send_stg(idir,i1,i2,i3)
                      shapes=(/sizes_srl_send_stg(idir,i1,i2,i3)/)
                      sendbuffer_srl(ibuf_start:ibuf_next-&
                           1)=reshape(psb(igrid)%ws(ixSmin1:ixSmax1,ixSmin2:ixSmax2,&
                           ixSmin3:ixSmax3,idir),shapes)
                      ibuf_start=ibuf_next
                   end do
                   isend_srl=isend_srl+1
                   call MPI_ISEND(sendbuffer_srl(ibuf_send_srl),&
                        sizes_srl_send_total(i1,i2,i3),MPI_DOUBLE_PRECISION,&
                        ipe_neighbor,itag,icomm,sendrequest_srl(isend_srl),ierrmpi)
                   ibuf_send_srl=ibuf_next
                end if
             else
                ixSmin1=ixS_srl_min1(iib1,i1);ixSmin2=ixS_srl_min2(iib2,i2)
                ixSmin3=ixS_srl_min3(iib3,i3);ixSmax1=ixS_srl_max1(iib1,i1)
                ixSmax2=ixS_srl_max2(iib2,i2);ixSmax3=ixS_srl_max3(iib3,i3);
                select case (ipole)
                case (1)
                   n_i1=i1;n_i2=-i2;n_i3=-i3;
                case (2)
                   n_i1=-i1;n_i2=i2;n_i3=-i3;
                case (3)
                   n_i1=-i1;n_i2=-i2;n_i3=i3;
                end select
                if (isend_buf(ipwbuf)/=0) then
                   call MPI_WAIT(sendrequest_c_sr(isend_buf(ipwbuf)),&
                        sendstatus_c_sr(:,isend_buf(ipwbuf)),ierrmpi)
                   deallocate(pwbuf(ipwbuf)%w)
                end if
                allocate(pwbuf(ipwbuf)%w(ixSmin1:ixSmax1,ixSmin2:ixSmax2,&
                     ixSmin3:ixSmax3,nwhead:nwtail))
                call pole_buffer(pwbuf(ipwbuf)%w,ixSmin1,ixSmin2,ixSmin3,ixSmax1,&
                     ixSmax2,ixSmax3,ixSmin1,ixSmin2,ixSmin3,ixSmax1,ixSmax2,ixSmax3,&
                     psb(igrid)%w,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,ixSmin1,&
                     ixSmin2,ixSmin3,ixSmax1,ixSmax2,ixSmax3)
                isend_c=isend_c+1
                isend_buf(ipwbuf)=isend_c
                itag=(3**3+4**3)*(ineighbor-1)+(n_i1+1)*3**(1-1)+(n_i2+1)*3**(2-1)+&
                     (n_i3+1)*3**(3-1)
                isizes=(ixSmax1-ixSmin1+1)*(ixSmax2-ixSmin2+1)*(ixSmax3-ixSmin3+&
                     1)*nwbc
                call MPI_ISEND(pwbuf(ipwbuf)%w,isizes,MPI_DOUBLE_PRECISION,&
                     ipe_neighbor,itag,icomm,sendrequest_c_sr(isend_c),ierrmpi)
                ipwbuf=1+modulo(ipwbuf,npwbuf)
                if(stagger_grid) then
                   ibuf_start=ibuf_send_srl
                   do idir=1,ndim
                      ixSmin1=ixS_srl_stg_min1(idir,i1)
                      ixSmin2=ixS_srl_stg_min2(idir,i2)
                      ixSmin3=ixS_srl_stg_min3(idir,i3)
                      ixSmax1=ixS_srl_stg_max1(idir,i1)
                      ixSmax2=ixS_srl_stg_max2(idir,i2)
                      ixSmax3=ixS_srl_stg_max3(idir,i3);
                      ibuf_next=ibuf_start+sizes_srl_send_stg(idir,i1,i2,i3)
                      shapes=(/sizes_srl_send_stg(idir,i1,i2,i3)/)
                      sendbuffer_srl(ibuf_start:ibuf_next-&
                           1)=reshape(psb(igrid)%ws(ixSmin1:ixSmax1,ixSmin2:ixSmax2,&
                           ixSmin3:ixSmax3,idir),shapes)
                      ibuf_start=ibuf_next
                   end do
                   isend_srl=isend_srl+1
                   call MPI_ISEND(sendbuffer_srl(ibuf_send_srl),&
                        sizes_srl_send_total(i1,i2,i3),MPI_DOUBLE_PRECISION,&
                        ipe_neighbor,itag,icomm,sendrequest_srl(isend_srl),ierrmpi)
                   ibuf_send_srl=ibuf_next
                end if
             end if
          end if

        end subroutine bc_send_srl

        !> Receive from fine neighbor
        subroutine bc_recv_restrict

          do ic3=1+int((1-i3)/2),2-int((1+i3)/2)
             inc3=2*i3+ic3
             do ic2=1+int((1-i2)/2),2-int((1+i2)/2)
                inc2=2*i2+ic2
                do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
                   inc1=2*i1+ic1
                   ipe_neighbor=neighbor_child(2,inc1,inc2,inc3,igrid)
                   if (ipe_neighbor/=mype) then
                      irecv_c=irecv_c+1
                      itag=(3**3+4**3)*(igrid-1)+3**3+inc1*4**(1-1)+inc2*4**(2-1)+&
                           inc3*4**(3-1)
                      call MPI_IRECV(psb(igrid)%w,1,type_recv_r(iib1,iib2,iib3,inc1,&
                           inc2,inc3), ipe_neighbor,itag,icomm,recvrequest_c_sr(irecv_c),&
                           ierrmpi)
                      if(stagger_grid) then
                         irecv_r=irecv_r+1
                         call MPI_IRECV(recvbuffer_r(ibuf_recv_r),&
                              sizes_r_recv_total(inc1,inc2,inc3), MPI_DOUBLE_PRECISION,&
                              ipe_neighbor,itag, icomm,recvrequest_r(irecv_r),ierrmpi)
                         ibuf_recv_r=ibuf_recv_r+sizes_r_recv_total(inc1,inc2,inc3)
                      end if
                   end if
                end do
             end do
          end do

        end subroutine bc_recv_restrict

        subroutine fill_coarse_boundary(igrid,i1,i2,i3)
          integer, intent(in) :: igrid,i1,i2,i3

          integer :: idims,iside,kmin1,kmin2,kmin3,kmax1,kmax2,kmax3,ixBmin1,&
               ixBmin2,ixBmin3,ixBmax1,ixBmax2,ixBmax3,ii1,ii2,ii3

          if(phyboundblock(igrid).and..not.stagger_grid.and.bcphys) then
             ! to use block in physical boundary setup for coarse representative
             block=>psc(igrid)
               ! filling physical boundary ghost cells of a coarser representative block for
               ! sending swap region with width of nghostcells to its coarser neighbor
               do idims=1,ndim
                  ! to avoid using as yet unknown corner info in more than 1D, we
                  ! fill only interior mesh ranges of the ghost cell ranges at first,
                  ! and progressively enlarge the ranges to include corners later
                  kmin1=merge(0, 1, idims==1)
                  kmax1=merge(0, 1, idims==1)
                  ixBmin1=ixCoGmin1+kmin1*nghostcells
                  ixBmax1=ixCoGmax1-kmax1*nghostcells
                  kmin2=merge(0, 1, idims==2)
                  kmax2=merge(0, 1, idims==2)
                  ixBmin2=ixCoGmin2+kmin2*nghostcells
                  ixBmax2=ixCoGmax2-kmax2*nghostcells
                  kmin3=merge(0, 1, idims==3)
                  kmax3=merge(0, 1, idims==3)
                  ixBmin3=ixCoGmin3+kmin3*nghostcells
                  ixBmax3=ixCoGmax3-kmax3*nghostcells


                  if(idims > 1 .and. neighbor_type(-1,0,0,&
                       igrid)==neighbor_boundary) ixBmin1=ixCoGmin1
                  if(idims > 1 .and. neighbor_type( 1,0,0,&
                       igrid)==neighbor_boundary) ixBmax1=ixCoGmax1
                  if(idims > 2 .and. neighbor_type(0,-1,0,&
                       igrid)==neighbor_boundary) ixBmin2=ixCoGmin2
                  if(idims > 2 .and. neighbor_type(0, 1,0,&
                       igrid)==neighbor_boundary) ixBmax2=ixCoGmax2
                  if(i1==-1) then
                     ixBmin1=ixCoGmin1+nghostcells
                     ixBmax1=ixCoGmin1+2*nghostcells-1
                  else if(i1==1) then
                     ixBmin1=ixCoGmax1-2*nghostcells+1
                     ixBmax1=ixCoGmax1-nghostcells
                  end if
                  if(i2==-1) then
                     ixBmin2=ixCoGmin2+nghostcells
                     ixBmax2=ixCoGmin2+2*nghostcells-1
                  else if(i2==1) then
                     ixBmin2=ixCoGmax2-2*nghostcells+1
                     ixBmax2=ixCoGmax2-nghostcells
                  end if
                  if(i3==-1) then
                     ixBmin3=ixCoGmin3+nghostcells
                     ixBmax3=ixCoGmin3+2*nghostcells-1
                  else if(i3==1) then
                     ixBmin3=ixCoGmax3-2*nghostcells+1
                     ixBmax3=ixCoGmax3-nghostcells
                  end if
                  do iside=1,2
                     ii1=kr(1,idims)*(2*iside-3);ii2=kr(2,idims)*(2*iside-3)
                     ii3=kr(3,idims)*(2*iside-3);
                     if (abs(i1)==1.and.abs(ii1)==1.or.abs(i2)==1.and.abs(ii2)==&
                          1.or.abs(i3)==1.and.abs(ii3)==1) cycle
                     if (neighbor_type(ii1,ii2,ii3,igrid)/=neighbor_boundary) cycle
                     call bc_phys(iside,idims,time,0.d0,psc(igrid),ixCoGmin1,&
                          ixCoGmin2,ixCoGmin3,ixCoGmax1,ixCoGmax2,ixCoGmax3,ixBmin1,&
                          ixBmin2,ixBmin3,ixBmax1,ixBmax2,ixBmax3)
                  end do
               end do
            end if

          end subroutine fill_coarse_boundary

          !> Send to coarser neighbor
          subroutine bc_send_restrict

            ipe_neighbor=neighbor(2,i1,i2,i3,igrid)
            if(ipe_neighbor/=mype) then
               ic1=1+modulo(node(pig1_,igrid)-1,2)
               ic2=1+modulo(node(pig2_,igrid)-1,2)
               ic3=1+modulo(node(pig3_,igrid)-1,2);
               if(.not.(i1==0.or.i1==2*ic1-3).or..not.(i2==0.or.i2==&
                    2*ic2-3).or..not.(i3==0.or.i3==2*ic3-3)) return
               ineighbor=neighbor(1,i1,i2,i3,igrid)
               ipole=neighbor_pole(i1,i2,i3,igrid)
               if(ipole==0) then
                  n_inc1=-2*i1+ic1;n_inc2=-2*i2+ic2;n_inc3=-2*i3+ic3;
                  isend_c=isend_c+1
                  itag=(3**3+4**3)*(ineighbor-1)+3**3+n_inc1*4**(1-1)+&
                       n_inc2*4**(2-1)+n_inc3*4**(3-1)
                  call MPI_ISEND(psc(igrid)%w,1,type_send_r(iib1,iib2,iib3,i1,i2,i3),&
                       ipe_neighbor,itag,icomm,sendrequest_c_sr(isend_c),ierrmpi)
                  if(stagger_grid) then
                     ibuf_start=ibuf_send_r
                     do idir=1,ndim
                        ixSmin1=ixS_r_stg_min1(idir,i1)
                        ixSmin2=ixS_r_stg_min2(idir,i2)
                        ixSmin3=ixS_r_stg_min3(idir,i3)
                        ixSmax1=ixS_r_stg_max1(idir,i1)
                        ixSmax2=ixS_r_stg_max2(idir,i2)
                        ixSmax3=ixS_r_stg_max3(idir,i3);
                        ibuf_next=ibuf_start+sizes_r_send_stg(idir,i1,i2,i3)
                        shapes=(/sizes_r_send_stg(idir,i1,i2,i3)/)
                        sendbuffer_r(ibuf_start:ibuf_next-&
                             1)=reshape(psc(igrid)%ws(ixSmin1:ixSmax1,ixSmin2:ixSmax2,&
                             ixSmin3:ixSmax3,idir),shapes)
                        ibuf_start=ibuf_next
                     end do
                     isend_r=isend_r+1
                     call MPI_ISEND(sendbuffer_r(ibuf_send_r),sizes_r_send_total(i1,&
                          i2,i3),MPI_DOUBLE_PRECISION,ipe_neighbor,itag, icomm,&
                          sendrequest_r(isend_r),ierrmpi)
                     ibuf_send_r=ibuf_next
                  end if
               else
                  ixSmin1=ixS_r_min1(iib1,i1);ixSmin2=ixS_r_min2(iib2,i2)
                  ixSmin3=ixS_r_min3(iib3,i3);ixSmax1=ixS_r_max1(iib1,i1)
                  ixSmax2=ixS_r_max2(iib2,i2);ixSmax3=ixS_r_max3(iib3,i3);
                  select case (ipole)
                  case (1)
                     n_inc1=2*i1+(3-ic1);n_inc2=-2*i2+ic2;n_inc3=-2*i3+ic3;
                  case (2)
                     n_inc1=-2*i1+ic1;n_inc2=2*i2+(3-ic2);n_inc3=-2*i3+ic3;
                  case (3)
                     n_inc1=-2*i1+ic1;n_inc2=-2*i2+ic2;n_inc3=2*i3+(3-ic3);
                  end select
                  if(isend_buf(ipwbuf)/=0) then
                     call MPI_WAIT(sendrequest_c_sr(isend_buf(ipwbuf)),&
                          sendstatus_c_sr(:,isend_buf(ipwbuf)),ierrmpi)
                     deallocate(pwbuf(ipwbuf)%w)
                  end if
                  allocate(pwbuf(ipwbuf)%w(ixSmin1:ixSmax1,ixSmin2:ixSmax2,&
                       ixSmin3:ixSmax3,nwhead:nwtail))
                  call pole_buffer(pwbuf(ipwbuf)%w,ixSmin1,ixSmin2,ixSmin3,ixSmax1,&
                       ixSmax2,ixSmax3,ixSmin1,ixSmin2,ixSmin3,ixSmax1,ixSmax2,ixSmax3,&
                       psc(igrid)%w,ixCoGmin1,ixCoGmin2,ixCoGmin3,ixCoGmax1,ixCoGmax2,&
                       ixCoGmax3,ixSmin1,ixSmin2,ixSmin3,ixSmax1,ixSmax2,ixSmax3)
                  isend_c=isend_c+1
                  isend_buf(ipwbuf)=isend_c
                  itag=(3**3+4**3)*(ineighbor-1)+3**3+n_inc1*4**(1-1)+&
                       n_inc2*4**(2-1)+n_inc3*4**(3-1)
                  isizes=(ixSmax1-ixSmin1+1)*(ixSmax2-ixSmin2+1)*(ixSmax3-ixSmin3+&
                       1)*nwbc
                  call MPI_ISEND(pwbuf(ipwbuf)%w,isizes,MPI_DOUBLE_PRECISION,&
                       ipe_neighbor,itag,icomm,sendrequest_c_sr(isend_c),ierrmpi)
                  ipwbuf=1+modulo(ipwbuf,npwbuf)
                  if(stagger_grid) then
                     ibuf_start=ibuf_send_r
                     do idir=1,ndim
                        ixSmin1=ixS_r_stg_min1(idir,i1)
                        ixSmin2=ixS_r_stg_min2(idir,i2)
                        ixSmin3=ixS_r_stg_min3(idir,i3)
                        ixSmax1=ixS_r_stg_max1(idir,i1)
                        ixSmax2=ixS_r_stg_max2(idir,i2)
                        ixSmax3=ixS_r_stg_max3(idir,i3);
                        ibuf_next=ibuf_start+sizes_r_send_stg(idir,i1,i2,i3)
                        shapes=(/sizes_r_send_stg(idir,i1,i2,i3)/)
                        sendbuffer_r(ibuf_start:ibuf_next-&
                             1)=reshape(psc(igrid)%ws(ixSmin1:ixSmax1,ixSmin2:ixSmax2,&
                             ixSmin3:ixSmax3,idir),shapes)
                        ibuf_start=ibuf_next
                     end do
                     isend_r=isend_r+1
                     call MPI_ISEND(sendbuffer_r(ibuf_send_r),sizes_r_send_total(i1,&
                          i2,i3),MPI_DOUBLE_PRECISION,ipe_neighbor,itag, icomm,&
                          sendrequest_r(isend_r),ierrmpi)
                     ibuf_send_r=ibuf_next
                  end if
               end if
            end if

          end subroutine bc_send_restrict

          !> fill coarser neighbor's ghost cells
          subroutine bc_fill_restrict(igrid,i1,i2,i3,iib1,iib2,iib3)
            integer, intent(in) :: igrid,i1,i2,i3,iib1,iib2,iib3

            integer :: ic1,ic2,ic3,n_inc1,n_inc2,n_inc3,ixSmin1,ixSmin2,ixSmin3,&
                 ixSmax1,ixSmax2,ixSmax3,ixRmin1,ixRmin2,ixRmin3,ixRmax1,ixRmax2,&
                 ixRmax3,ipe_neighbor,ineighbor,ipole,idir

            ipe_neighbor=neighbor(2,i1,i2,i3,igrid)
            if(ipe_neighbor==mype) then
               ic1=1+modulo(node(pig1_,igrid)-1,2)
               ic2=1+modulo(node(pig2_,igrid)-1,2)
               ic3=1+modulo(node(pig3_,igrid)-1,2);
               if(.not.(i1==0.or.i1==2*ic1-3).or..not.(i2==0.or.i2==&
                    2*ic2-3).or..not.(i3==0.or.i3==2*ic3-3)) return
               ineighbor=neighbor(1,i1,i2,i3,igrid)
               ipole=neighbor_pole(i1,i2,i3,igrid)
               if(ipole==0) then
                  n_inc1=-2*i1+ic1;n_inc2=-2*i2+ic2;n_inc3=-2*i3+ic3;
                  ixSmin1=ixS_r_min1(iib1,i1);ixSmin2=ixS_r_min2(iib2,i2)
                  ixSmin3=ixS_r_min3(iib3,i3);ixSmax1=ixS_r_max1(iib1,i1)
                  ixSmax2=ixS_r_max2(iib2,i2);ixSmax3=ixS_r_max3(iib3,i3);
                  ixRmin1=ixR_r_min1(iib1,n_inc1);ixRmin2=ixR_r_min2(iib2,n_inc2)
                  ixRmin3=ixR_r_min3(iib3,n_inc3);ixRmax1=ixR_r_max1(iib1,n_inc1)
                  ixRmax2=ixR_r_max2(iib2,n_inc2);ixRmax3=ixR_r_max3(iib3,n_inc3);
                  psb(ineighbor)%w(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,&
                       nwhead:nwtail)=psc(igrid)%w(ixSmin1:ixSmax1,ixSmin2:ixSmax2,&
                       ixSmin3:ixSmax3,nwhead:nwtail)
                  if(stagger_grid) then
                     do idir=1,ndim
                        ixSmin1=ixS_r_stg_min1(idir,i1)
                        ixSmin2=ixS_r_stg_min2(idir,i2)
                        ixSmin3=ixS_r_stg_min3(idir,i3)
                        ixSmax1=ixS_r_stg_max1(idir,i1)
                        ixSmax2=ixS_r_stg_max2(idir,i2)
                        ixSmax3=ixS_r_stg_max3(idir,i3);
                        ixRmin1=ixR_r_stg_min1(idir,n_inc1)
                        ixRmin2=ixR_r_stg_min2(idir,n_inc2)
                        ixRmin3=ixR_r_stg_min3(idir,n_inc3)
                        ixRmax1=ixR_r_stg_max1(idir,n_inc1)
                        ixRmax2=ixR_r_stg_max2(idir,n_inc2)
                        ixRmax3=ixR_r_stg_max3(idir,n_inc3);
                        psb(ineighbor)%ws(ixRmin1:ixRmax1,ixRmin2:ixRmax2,&
                             ixRmin3:ixRmax3,idir)=psc(igrid)%ws(ixSmin1:ixSmax1,&
                             ixSmin2:ixSmax2,ixSmin3:ixSmax3,idir)
                     end do
                  end if
               else
                  ixSmin1=ixS_r_min1(iib1,i1);ixSmin2=ixS_r_min2(iib2,i2)
                  ixSmin3=ixS_r_min3(iib3,i3);ixSmax1=ixS_r_max1(iib1,i1)
                  ixSmax2=ixS_r_max2(iib2,i2);ixSmax3=ixS_r_max3(iib3,i3);
                  select case (ipole)
                  case (1)
                     n_inc1=2*i1+(3-ic1);n_inc2=-2*i2+ic2;n_inc3=-2*i3+ic3;
                  case (2)
                     n_inc1=-2*i1+ic1;n_inc2=2*i2+(3-ic2);n_inc3=-2*i3+ic3;
                  case (3)
                     n_inc1=-2*i1+ic1;n_inc2=-2*i2+ic2;n_inc3=2*i3+(3-ic3);
                  end select
                  ixRmin1=ixR_r_min1(iib1,n_inc1);ixRmin2=ixR_r_min2(iib2,n_inc2)
                  ixRmin3=ixR_r_min3(iib3,n_inc3);ixRmax1=ixR_r_max1(iib1,n_inc1)
                  ixRmax2=ixR_r_max2(iib2,n_inc2);ixRmax3=ixR_r_max3(iib3,n_inc3);
                  call pole_copy(psb(ineighbor)%w,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,&
                       ixGhi3,ixRmin1,ixRmin2,ixRmin3,ixRmax1,ixRmax2,ixRmax3,&
                       psc(igrid)%w,ixCoGmin1,ixCoGmin2,ixCoGmin3,ixCoGmax1,ixCoGmax2,&
                       ixCoGmax3,ixSmin1,ixSmin2,ixSmin3,ixSmax1,ixSmax2,ixSmax3,&
                       ipole)
                  if(stagger_grid) then
                     do idir=1,ndim
                        ixSmin1=ixS_r_stg_min1(idir,i1)
                        ixSmin2=ixS_r_stg_min2(idir,i2)
                        ixSmin3=ixS_r_stg_min3(idir,i3)
                        ixSmax1=ixS_r_stg_max1(idir,i1)
                        ixSmax2=ixS_r_stg_max2(idir,i2)
                        ixSmax3=ixS_r_stg_max3(idir,i3);
                        ixRmin1=ixR_r_stg_min1(idir,n_inc1)
                        ixRmin2=ixR_r_stg_min2(idir,n_inc2)
                        ixRmin3=ixR_r_stg_min3(idir,n_inc3)
                        ixRmax1=ixR_r_stg_max1(idir,n_inc1)
                        ixRmax2=ixR_r_stg_max2(idir,n_inc2)
                        ixRmax3=ixR_r_stg_max3(idir,n_inc3);
                        !! Fill ghost cells
                        call pole_copy_stg(psb(ineighbor)%ws,ixGslo1,ixGslo2,ixGslo3,&
                             ixGshi1,ixGshi2,ixGshi3,ixRmin1,ixRmin2,ixRmin3,ixRmax1,&
                             ixRmax2,ixRmax3,psc(igrid)%ws,ixCoGsmin1,ixCoGsmin2,&
                             ixCoGsmin3,ixCoGsmax1,ixCoGsmax2,ixCoGsmax3,ixSmin1,ixSmin2,&
                             ixSmin3,ixSmax1,ixSmax2,ixSmax3,idir,ipole)
                     end do
                  end if
               end if
            end if

          end subroutine bc_fill_restrict

          !> fill siblings ghost cells with received data
          subroutine bc_fill_srl_stg
            double precision :: tmp(ixGslo1:ixGshi1,ixGslo2:ixGshi2,&
                 ixGslo3:ixGshi3)
            integer :: ixSmin1,ixSmin2,ixSmin3,ixSmax1,ixSmax2,ixSmax3,ixRmin1,&
                 ixRmin2,ixRmin3,ixRmax1,ixRmax2,ixRmax3,n_i1,n_i2,n_i3,ixSsyncmin1,&
                 ixSsyncmin2,ixSsyncmin3,ixSsyncmax1,ixSsyncmax2,ixSsyncmax3,&
                 ixRsyncmin1,ixRsyncmin2,ixRsyncmin3,ixRsyncmax1,ixRsyncmax2,&
                 ixRsyncmax3
            integer :: idir

            ipe_neighbor=neighbor(2,i1,i2,i3,igrid)
            if(ipe_neighbor/=mype) then
               ineighbor=neighbor(1,i1,i2,i3,igrid)
               ipole=neighbor_pole(i1,i2,i3,igrid)

               !! Now the special treatment of the pole is done here, at the receive step
               if (ipole==0) then
                  ixRmin1=ixR_srl_min1(iib1,i1);ixRmin2=ixR_srl_min2(iib2,i2)
                  ixRmin3=ixR_srl_min3(iib3,i3);ixRmax1=ixR_srl_max1(iib1,i1)
                  ixRmax2=ixR_srl_max2(iib2,i2);ixRmax3=ixR_srl_max3(iib3,i3);
                  !! Unpack the buffer and fill the ghost cells
                  n_i1=-i1;n_i2=-i2;n_i3=-i3;
                  do idir=1,ndim
                     ixSmin1=ixS_srl_stg_min1(idir,n_i1)
                     ixSmin2=ixS_srl_stg_min2(idir,n_i2)
                     ixSmin3=ixS_srl_stg_min3(idir,n_i3)
                     ixSmax1=ixS_srl_stg_max1(idir,n_i1)
                     ixSmax2=ixS_srl_stg_max2(idir,n_i2)
                     ixSmax3=ixS_srl_stg_max3(idir,n_i3);
                     ixRmin1=ixR_srl_stg_min1(idir,i1)
                     ixRmin2=ixR_srl_stg_min2(idir,i2)
                     ixRmin3=ixR_srl_stg_min3(idir,i3)
                     ixRmax1=ixR_srl_stg_max1(idir,i1)
                     ixRmax2=ixR_srl_stg_max2(idir,i2)
                     ixRmax3=ixR_srl_stg_max3(idir,i3);
                     ibuf_next=ibuf_recv_srl+sizes_srl_recv_stg(idir,i1,i2,i3)
                     tmp(ixSmin1:ixSmax1,ixSmin2:ixSmax2,&
                          ixSmin3:ixSmax3) = reshape(source=recvbuffer_srl(&
                          ibuf_recv_srl:ibuf_next-1),&
                          shape=shape(psb(igrid)%ws(ixSmin1:ixSmax1,ixSmin2:ixSmax2,&
                          ixSmin3:ixSmax3,idir)))
                     psb(igrid)%ws(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,&
                          idir) = tmp(ixSmin1:ixSmax1,ixSmin2:ixSmax2,ixSmin3:ixSmax3)
                     ibuf_recv_srl=ibuf_next
                  end do
               else ! There is a pole
                  select case (ipole)
                  case (1)
                     n_i1=i1;n_i2=-i2;n_i3=-i3;
                  case (2)
                     n_i1=-i1;n_i2=i2;n_i3=-i3;
                  case (3)
                     n_i1=-i1;n_i2=-i2;n_i3=i3;
                  end select
                  pole_buf%ws=zero
                  do idir=1,ndim
                     ixRmin1=ixR_srl_stg_min1(idir,i1)
                     ixRmin2=ixR_srl_stg_min2(idir,i2)
                     ixRmin3=ixR_srl_stg_min3(idir,i3)
                     ixRmax1=ixR_srl_stg_max1(idir,i1)
                     ixRmax2=ixR_srl_stg_max2(idir,i2)
                     ixRmax3=ixR_srl_stg_max3(idir,i3);
                     ixSmin1=ixS_srl_stg_min1(idir,n_i1)
                     ixSmin2=ixS_srl_stg_min2(idir,n_i2)
                     ixSmin3=ixS_srl_stg_min3(idir,n_i3)
                     ixSmax1=ixS_srl_stg_max1(idir,n_i1)
                     ixSmax2=ixS_srl_stg_max2(idir,n_i2)
                     ixSmax3=ixS_srl_stg_max3(idir,n_i3);
                     ibuf_next=ibuf_recv_srl+sizes_srl_recv_stg(idir,i1,i2,i3)
                     pole_buf%ws(ixSmin1:ixSmax1,ixSmin2:ixSmax2,ixSmin3:ixSmax3,&
                          idir)=reshape(source=recvbuffer_srl(ibuf_recv_srl:ibuf_next-1),&
                          shape=shape(psb(igrid)%ws(ixSmin1:ixSmax1,ixSmin2:ixSmax2,&
                          ixSmin3:ixSmax3,idir)))
                     ibuf_recv_srl=ibuf_next
                     call pole_copy_stg(psb(igrid)%ws,ixGslo1,ixGslo2,ixGslo3,ixGshi1,&
                          ixGshi2,ixGshi3,ixRmin1,ixRmin2,ixRmin3,ixRmax1,ixRmax2,&
                          ixRmax3,pole_buf%ws,ixGslo1,ixGslo2,ixGslo3,ixGshi1,ixGshi2,&
                          ixGshi3,ixSmin1,ixSmin2,ixSmin3,ixSmax1,ixSmax2,ixSmax3,idir,&
                          ipole)
                  end do
               end if
            end if

          end subroutine bc_fill_srl_stg

          subroutine indices_for_syncing(idir,i1,i2,i3,ixRmin1,ixRmin2,ixRmin3,&
               ixRmax1,ixRmax2,ixRmax3,ixSmin1,ixSmin2,ixSmin3,ixSmax1,ixSmax2,&
               ixSmax3,ixRsyncmin1,ixRsyncmin2,ixRsyncmin3,ixRsyncmax1,ixRsyncmax2,&
               ixRsyncmax3,ixSsyncmin1,ixSsyncmin2,ixSsyncmin3,ixSsyncmax1,&
               ixSsyncmax2,ixSsyncmax3)
            integer, intent(in)       :: i1,i2,i3,idir
            integer, intent(inout)    :: ixRmin1,ixRmin2,ixRmin3,ixRmax1,ixRmax2,&
                 ixRmax3,ixSmin1,ixSmin2,ixSmin3,ixSmax1,ixSmax2,ixSmax3
            integer, intent(out)      :: ixRsyncmin1,ixRsyncmin2,ixRsyncmin3,&
                 ixRsyncmax1,ixRsyncmax2,ixRsyncmax3,ixSsyncmin1,ixSsyncmin2,&
                 ixSsyncmin3,ixSsyncmax1,ixSsyncmax2,ixSsyncmax3

            ixRsyncmin1=ixRmin1;ixRsyncmin2=ixRmin2;ixRsyncmin3=ixRmin3
            ixRsyncmax1=ixRmax1;ixRsyncmax2=ixRmax2;ixRsyncmax3=ixRmax3;
            ixSsyncmin1=ixSmin1;ixSsyncmin2=ixSmin2;ixSsyncmin3=ixSmin3
            ixSsyncmax1=ixSmax1;ixSsyncmax2=ixSmax2;ixSsyncmax3=ixSmax3;


            if (i1 == -1 .and. idir == 1) then
               ixRsyncmin1 = ixRmax1
               ixRsyncmax1 = ixRmax1
               ixSsyncmin1 = ixSmax1
               ixSsyncmax1 = ixSmax1
               ixRmax1 = ixRmax1 - 1
               ixSmax1 = ixSmax1 - 1
            else if (i1 == 1 .and. idir == 1) then
               ixRsyncmin1 = ixRmin1
               ixRsyncmax1 = ixRmin1
               ixSsyncmin1 = ixSmin1
               ixSsyncmax1 = ixSmin1
               ixRmin1 = ixRmin1 + 1
               ixSmin1 = ixSmin1 + 1
            end if


            if (i2 == -1 .and. idir == 2) then
               ixRsyncmin2 = ixRmax2
               ixRsyncmax2 = ixRmax2
               ixSsyncmin2 = ixSmax2
               ixSsyncmax2 = ixSmax2
               ixRmax2 = ixRmax2 - 1
               ixSmax2 = ixSmax2 - 1
            else if (i2 == 1 .and. idir == 2) then
               ixRsyncmin2 = ixRmin2
               ixRsyncmax2 = ixRmin2
               ixSsyncmin2 = ixSmin2
               ixSsyncmax2 = ixSmin2
               ixRmin2 = ixRmin2 + 1
               ixSmin2 = ixSmin2 + 1
            end if


            if (i3 == -1 .and. idir == 3) then
               ixRsyncmin3 = ixRmax3
               ixRsyncmax3 = ixRmax3
               ixSsyncmin3 = ixSmax3
               ixSsyncmax3 = ixSmax3
               ixRmax3 = ixRmax3 - 1
               ixSmax3 = ixSmax3 - 1
            else if (i3 == 1 .and. idir == 3) then
               ixRsyncmin3 = ixRmin3
               ixRsyncmax3 = ixRmin3
               ixSsyncmin3 = ixSmin3
               ixSsyncmax3 = ixSmin3
               ixRmin3 = ixRmin3 + 1
               ixSmin3 = ixSmin3 + 1
            end if


          end subroutine indices_for_syncing

          !> fill restricted ghost cells after receipt
          subroutine bc_fill_restrict_stg

            ipole=neighbor_pole(i1,i2,i3,igrid)
            if (ipole==0) then
               ! Loop over the children ic^D to and their neighbors inc^D
               do ic3=1+int((1-i3)/2),2-int((1+i3)/2)
                  inc3=2*i3+ic3
                  do ic2=1+int((1-i2)/2),2-int((1+i2)/2)
                     inc2=2*i2+ic2
                     do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
                        inc1=2*i1+ic1
                        ipe_neighbor=neighbor_child(2,inc1,inc2,inc3,igrid)
                        if(ipe_neighbor/=mype) then
                           ineighbor=neighbor_child(1,inc1,inc2,inc3,igrid)
                           n_i1=-i1;n_i2=-i2;n_i3=-i3;
                           !! Unpack the buffer and fill the ghost cells
                           do idir=1,ndim
                              ixRmin1=ixR_r_stg_min1(idir,inc1)
                              ixRmin2=ixR_r_stg_min2(idir,inc2)
                              ixRmin3=ixR_r_stg_min3(idir,inc3)
                              ixRmax1=ixR_r_stg_max1(idir,inc1)
                              ixRmax2=ixR_r_stg_max2(idir,inc2)
                              ixRmax3=ixR_r_stg_max3(idir,inc3);
                              ibuf_next=ibuf_recv_r+sizes_r_recv_stg(idir,inc1,inc2,inc3)
                              psb(igrid)%ws(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,&
                                   idir)=reshape(source=recvbuffer_r(ibuf_recv_r:ibuf_next-1),&
                                   shape=shape(psb(igrid)%ws(ixRmin1:ixRmax1,ixRmin2:ixRmax2,&
                                   ixRmin3:ixRmax3,idir)))
                              ibuf_recv_r=ibuf_next
                           end do
                        end if
                     end do
                  end do
               end do
            else !! There is a pole
               do ic3=1+int((1-i3)/2),2-int((1+i3)/2)
                  inc3=2*i3+ic3
                  do ic2=1+int((1-i2)/2),2-int((1+i2)/2)
                     inc2=2*i2+ic2
                     do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
                        inc1=2*i1+ic1
                        ipe_neighbor=neighbor_child(2,inc1,inc2,inc3,igrid)
                        if(ipe_neighbor/=mype) then
                           ineighbor=neighbor_child(1,inc1,inc2,inc3,igrid)
                           select case(ipole)
                           case (1)
                              n_i1=i1;n_i2=-i2;n_i3=-i3;
                           case (2)
                              n_i1=-i1;n_i2=i2;n_i3=-i3;
                           case (3)
                              n_i1=-i1;n_i2=-i2;n_i3=i3;
                           end select
                           ixRmin1=ixR_r_min1(iib1,inc1);ixRmin2=ixR_r_min2(iib2,inc2)
                           ixRmin3=ixR_r_min3(iib3,inc3);ixRmax1=ixR_r_max1(iib1,inc1)
                           ixRmax2=ixR_r_max2(iib2,inc2);ixRmax3=ixR_r_max3(iib3,inc3);
                           !! Unpack the buffer and fill an auxiliary array
                           pole_buf%ws=zero
                           do idir=1,ndim
                              ixSmin1=ixS_r_stg_min1(idir,n_i1)
                              ixSmin2=ixS_r_stg_min2(idir,n_i2)
                              ixSmin3=ixS_r_stg_min3(idir,n_i3)
                              ixSmax1=ixS_r_stg_max1(idir,n_i1)
                              ixSmax2=ixS_r_stg_max2(idir,n_i2)
                              ixSmax3=ixS_r_stg_max3(idir,n_i3);
                              ixRmin1=ixR_r_stg_min1(idir,inc1)
                              ixRmin2=ixR_r_stg_min2(idir,inc2)
                              ixRmin3=ixR_r_stg_min3(idir,inc3)
                              ixRmax1=ixR_r_stg_max1(idir,inc1)
                              ixRmax2=ixR_r_stg_max2(idir,inc2)
                              ixRmax3=ixR_r_stg_max3(idir,inc3);
                              ibuf_next=ibuf_recv_r+sizes_r_recv_stg(idir,inc1,inc2,inc3)
                              pole_buf%ws(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,&
                                   idir)=reshape(source=recvbuffer_r(ibuf_recv_r:ibuf_next-1),&
                                   shape=shape(psb(igrid)%ws(ixRmin1:ixRmax1,ixRmin2:ixRmax2,&
                                   ixRmin3:ixRmax3,idir)))
                              call pole_copy_stg(psb(igrid)%ws,ixGslo1,ixGslo2,ixGslo3,&
                                   ixGshi1,ixGshi2,ixGshi3,ixRmin1,ixRmin2,ixRmin3,ixRmax1,&
                                   ixRmax2,ixRmax3,pole_buf%ws,ixGslo1,ixGslo2,ixGslo3,&
                                   ixGshi1,ixGshi2,ixGshi3,ixRmin1,ixRmin2,ixRmin3,ixRmax1,&
                                   ixRmax2,ixRmax3,idir,ipole)
                              ibuf_recv_r=ibuf_next
                           end do
                        end if
                     end do
                  end do
               end do
            end if

          end subroutine bc_fill_restrict_stg

          !> Receive from coarse neighbor
          subroutine bc_recv_prolong

            ic1=1+modulo(node(pig1_,igrid)-1,2)
            ic2=1+modulo(node(pig2_,igrid)-1,2)
            ic3=1+modulo(node(pig3_,igrid)-1,2);
            if (.not.(i1==0.or.i1==2*ic1-3).or..not.(i2==0.or.i2==&
                 2*ic2-3).or..not.(i3==0.or.i3==2*ic3-3)) return

            ipe_neighbor=neighbor(2,i1,i2,i3,igrid)
            if (ipe_neighbor/=mype) then
               irecv_c=irecv_c+1
               inc1=ic1+i1;inc2=ic2+i2;inc3=ic3+i3;
               itag=(3**3+4**3)*(igrid-1)+3**3+inc1*4**(1-1)+inc2*4**(2-1)+&
                    inc3*4**(3-1)
               call MPI_IRECV(psc(igrid)%w,1,type_recv_p(iib1,iib2,iib3,inc1,inc2,&
                    inc3), ipe_neighbor,itag,icomm,recvrequest_c_p(irecv_c),ierrmpi)
               if(stagger_grid) then
                  irecv_p=irecv_p+1
                  call MPI_IRECV(recvbuffer_p(ibuf_recv_p),sizes_p_recv_total(inc1,&
                       inc2,inc3),MPI_DOUBLE_PRECISION,ipe_neighbor,itag,icomm,&
                       recvrequest_p(irecv_p),ierrmpi)
                  ibuf_recv_p=ibuf_recv_p+sizes_p_recv_total(inc1,inc2,inc3)
               end if
            end if

          end subroutine bc_recv_prolong

          !> Send to finer neighbor
          subroutine bc_send_prolong
            integer :: ii1,ii2,ii3

            ipole=neighbor_pole(i1,i2,i3,igrid)

            do ic3=1+int((1-i3)/2),2-int((1+i3)/2)
               inc3=2*i3+ic3
               do ic2=1+int((1-i2)/2),2-int((1+i2)/2)
                  inc2=2*i2+ic2
                  do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
                     inc1=2*i1+ic1
                     ipe_neighbor=neighbor_child(2,inc1,inc2,inc3,igrid)
                     if(ipe_neighbor/=mype) then
                        ixSmin1=ixS_p_min1(iib1,inc1);ixSmin2=ixS_p_min2(iib2,inc2)
                        ixSmin3=ixS_p_min3(iib3,inc3);ixSmax1=ixS_p_max1(iib1,inc1)
                        ixSmax2=ixS_p_max2(iib2,inc2);ixSmax3=ixS_p_max3(iib3,inc3);
                        ineighbor=neighbor_child(1,inc1,inc2,inc3,igrid)
                        if(ipole==0) then
                           n_i1=-i1;n_i2=-i2;n_i3=-i3;
                           n_inc1=ic1+n_i1;n_inc2=ic2+n_i2;n_inc3=ic3+n_i3;
                           isend_c=isend_c+1
                           itag=(3**3+4**3)*(ineighbor-1)+3**3+n_inc1*4**(1-1)+&
                                n_inc2*4**(2-1)+n_inc3*4**(3-1)
                           call MPI_ISEND(psb(igrid)%w,1,type_send_p(iib1,iib2,iib3,inc1,&
                                inc2,inc3), ipe_neighbor,itag,icomm,sendrequest_c_p(isend_c),&
                                ierrmpi)
                           if(stagger_grid) then
                              ibuf_start=ibuf_send_p
                              do idir=1,ndim
                                 ixSmin1=ixS_p_stg_min1(idir,inc1)
                                 ixSmin2=ixS_p_stg_min2(idir,inc2)
                                 ixSmin3=ixS_p_stg_min3(idir,inc3)
                                 ixSmax1=ixS_p_stg_max1(idir,inc1)
                                 ixSmax2=ixS_p_stg_max2(idir,inc2)
                                 ixSmax3=ixS_p_stg_max3(idir,inc3);
                                 ibuf_next=ibuf_start+sizes_p_send_stg(idir,inc1,inc2,inc3)
                                 shapes=(/sizes_p_send_stg(idir,inc1,inc2,inc3)/)
                                 sendbuffer_p(ibuf_start:ibuf_next-&
                                      1)=reshape(psb(igrid)%ws(ixSmin1:ixSmax1,ixSmin2:ixSmax2,&
                                      ixSmin3:ixSmax3,idir),shapes)
                                 ibuf_start=ibuf_next
                              end do
                              isend_p=isend_p+1
                              call MPI_ISEND(sendbuffer_p(ibuf_send_p),&
                                   sizes_p_send_total(inc1,inc2,inc3),MPI_DOUBLE_PRECISION,&
                                   ipe_neighbor,itag, icomm,sendrequest_p(isend_p),ierrmpi)
                              ibuf_send_p=ibuf_next
                           end if
                        else
                           select case (ipole)
                           case (1)
                              n_inc1=inc1;n_inc2=ic2-i2;n_inc3=ic3-i3;
                           case (2)
                              n_inc1=ic1-i1;n_inc2=inc2;n_inc3=ic3-i3;
                           case (3)
                              n_inc1=ic1-i1;n_inc2=ic2-i2;n_inc3=inc3;
                           end select
                           if(isend_buf(ipwbuf)/=0) then
                              call MPI_WAIT(sendrequest_c_p(isend_buf(ipwbuf)),&
                                   sendstatus_c_p(:,isend_buf(ipwbuf)),ierrmpi)
                              deallocate(pwbuf(ipwbuf)%w)
                           end if
                           allocate(pwbuf(ipwbuf)%w(ixSmin1:ixSmax1,ixSmin2:ixSmax2,&
                                ixSmin3:ixSmax3,nwhead:nwtail))
                           call pole_buffer(pwbuf(ipwbuf)%w,ixSmin1,ixSmin2,ixSmin3,&
                                ixSmax1,ixSmax2,ixSmax3,ixSmin1,ixSmin2,ixSmin3,ixSmax1,&
                                ixSmax2,ixSmax3,psb(igrid)%w,ixGlo1,ixGlo2,ixGlo3,ixGhi1,&
                                ixGhi2,ixGhi3,ixSmin1,ixSmin2,ixSmin3,ixSmax1,ixSmax2,&
                                ixSmax3)
                           isend_c=isend_c+1
                           isend_buf(ipwbuf)=isend_c
                           itag=(3**3+4**3)*(ineighbor-1)+3**3+n_inc1*4**(1-1)+&
                                n_inc2*4**(2-1)+n_inc3*4**(3-1)
                           isizes=(ixSmax1-ixSmin1+1)*(ixSmax2-ixSmin2+1)*(ixSmax3-ixSmin3+&
                                1)*nwbc
                           call MPI_ISEND(pwbuf(ipwbuf)%w,isizes,MPI_DOUBLE_PRECISION,&
                                ipe_neighbor,itag,icomm,sendrequest_c_p(isend_c),ierrmpi)
                           ipwbuf=1+modulo(ipwbuf,npwbuf)
                           if(stagger_grid) then
                              ibuf_start=ibuf_send_p
                              do idir=1,ndim
                                 ixSmin1=ixS_p_stg_min1(idir,inc1)
                                 ixSmin2=ixS_p_stg_min2(idir,inc2)
                                 ixSmin3=ixS_p_stg_min3(idir,inc3)
                                 ixSmax1=ixS_p_stg_max1(idir,inc1)
                                 ixSmax2=ixS_p_stg_max2(idir,inc2)
                                 ixSmax3=ixS_p_stg_max3(idir,inc3);
                                 ibuf_next=ibuf_start+sizes_p_send_stg(idir,inc1,inc2,inc3)
                                 shapes=(/sizes_p_send_stg(idir,inc1,inc2,inc3)/)
                                 sendbuffer_p(ibuf_start:ibuf_next-&
                                      1)=reshape(psb(igrid)%ws(ixSmin1:ixSmax1,ixSmin2:ixSmax2,&
                                      ixSmin3:ixSmax3,idir),shapes)
                                 ibuf_start=ibuf_next
                              end do
                              isend_p=isend_p+1
                              call MPI_ISEND(sendbuffer_p(ibuf_send_p),&
                                   sizes_p_send_total(inc1,inc2,inc3),MPI_DOUBLE_PRECISION,&
                                   ipe_neighbor,itag, icomm,sendrequest_p(isend_p),ierrmpi)
                              ibuf_send_p=ibuf_next
                           end if
                        end if
                     end if
                  end do
               end do
            end do

          end subroutine bc_send_prolong

          !> Send to finer neighbor
          subroutine bc_fill_prolong(igrid,i1,i2,i3,iib1,iib2,iib3)
            integer, intent(in) :: igrid,i1,i2,i3,iib1,iib2,iib3

            integer :: ipe_neighbor,ineighbor,ixSmin1,ixSmin2,ixSmin3,ixSmax1,&
                 ixSmax2,ixSmax3,ixRmin1,ixRmin2,ixRmin3,ixRmax1,ixRmax2,ixRmax3,ic1,&
                 ic2,ic3,inc1,inc2,inc3,ipole,idir

            ipole=neighbor_pole(i1,i2,i3,igrid)

            if(ipole==0) then
               do ic3=1+int((1-i3)/2),2-int((1+i3)/2)
                  inc3=2*i3+ic3
                  do ic2=1+int((1-i2)/2),2-int((1+i2)/2)
                     inc2=2*i2+ic2
                     do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
                        inc1=2*i1+ic1
                        ipe_neighbor=neighbor_child(2,inc1,inc2,inc3,igrid)
                        if(ipe_neighbor==mype) then
                           ixSmin1=ixS_p_min1(iib1,inc1);ixSmin2=ixS_p_min2(iib2,inc2)
                           ixSmin3=ixS_p_min3(iib3,inc3);ixSmax1=ixS_p_max1(iib1,inc1)
                           ixSmax2=ixS_p_max2(iib2,inc2);ixSmax3=ixS_p_max3(iib3,inc3);
                           ineighbor=neighbor_child(1,inc1,inc2,inc3,igrid)
                           ipole=neighbor_pole(i1,i2,i3,igrid)
                           n_i1=-i1;n_i2=-i2;n_i3=-i3;
                           n_inc1=ic1+n_i1;n_inc2=ic2+n_i2;n_inc3=ic3+n_i3;
                           ixRmin1=ixR_p_min1(iib1,n_inc1)
                           ixRmin2=ixR_p_min2(iib2,n_inc2)
                           ixRmin3=ixR_p_min3(iib3,n_inc3)
                           ixRmax1=ixR_p_max1(iib1,n_inc1)
                           ixRmax2=ixR_p_max2(iib2,n_inc2)
                           ixRmax3=ixR_p_max3(iib3,n_inc3);
                           psc(ineighbor)%w(ixRmin1:ixRmax1,ixRmin2:ixRmax2,&
                                ixRmin3:ixRmax3,nwhead:nwtail) =psb(igrid)%w(ixSmin1:ixSmax1,&
                                ixSmin2:ixSmax2,ixSmin3:ixSmax3,nwhead:nwtail)
                           if(stagger_grid) then
                              do idir=1,ndim
                                 ixSmin1=ixS_p_stg_min1(idir,inc1)
                                 ixSmin2=ixS_p_stg_min2(idir,inc2)
                                 ixSmin3=ixS_p_stg_min3(idir,inc3)
                                 ixSmax1=ixS_p_stg_max1(idir,inc1)
                                 ixSmax2=ixS_p_stg_max2(idir,inc2)
                                 ixSmax3=ixS_p_stg_max3(idir,inc3);
                                 ixRmin1=ixR_p_stg_min1(idir,n_inc1)
                                 ixRmin2=ixR_p_stg_min2(idir,n_inc2)
                                 ixRmin3=ixR_p_stg_min3(idir,n_inc3)
                                 ixRmax1=ixR_p_stg_max1(idir,n_inc1)
                                 ixRmax2=ixR_p_stg_max2(idir,n_inc2)
                                 ixRmax3=ixR_p_stg_max3(idir,n_inc3);
                                 psc(ineighbor)%ws(ixRmin1:ixRmax1,ixRmin2:ixRmax2,&
                                      ixRmin3:ixRmax3,idir)=psb(igrid)%ws(ixSmin1:ixSmax1,&
                                      ixSmin2:ixSmax2,ixSmin3:ixSmax3,idir)
                              end do
                           end if
                        end if
                     end do
                  end do
               end do
            else
               do ic3=1+int((1-i3)/2),2-int((1+i3)/2)
                  inc3=2*i3+ic3
                  do ic2=1+int((1-i2)/2),2-int((1+i2)/2)
                     inc2=2*i2+ic2
                     do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
                        inc1=2*i1+ic1
                        ipe_neighbor=neighbor_child(2,inc1,inc2,inc3,igrid)
                        if(ipe_neighbor==mype) then
                           ixSmin1=ixS_p_min1(iib1,inc1);ixSmin2=ixS_p_min2(iib2,inc2)
                           ixSmin3=ixS_p_min3(iib3,inc3);ixSmax1=ixS_p_max1(iib1,inc1)
                           ixSmax2=ixS_p_max2(iib2,inc2);ixSmax3=ixS_p_max3(iib3,inc3);
                           ineighbor=neighbor_child(1,inc1,inc2,inc3,igrid)
                           ipole=neighbor_pole(i1,i2,i3,igrid)
                           select case (ipole)
                           case (1)
                              n_inc1=inc1;n_inc2=ic2-i2;n_inc3=ic3-i3;
                           case (2)
                              n_inc1=ic1-i1;n_inc2=inc2;n_inc3=ic3-i3;
                           case (3)
                              n_inc1=ic1-i1;n_inc2=ic2-i2;n_inc3=inc3;
                           end select
                           ixRmin1=ixR_p_min1(iib1,n_inc1)
                           ixRmin2=ixR_p_min2(iib2,n_inc2)
                           ixRmin3=ixR_p_min3(iib3,n_inc3)
                           ixRmax1=ixR_p_max1(iib1,n_inc1)
                           ixRmax2=ixR_p_max2(iib2,n_inc2)
                           ixRmax3=ixR_p_max3(iib3,n_inc3);
                           call pole_copy(psc(ineighbor)%w,ixCoGmin1,ixCoGmin2,ixCoGmin3,&
                                ixCoGmax1,ixCoGmax2,ixCoGmax3,ixRmin1,ixRmin2,ixRmin3,&
                                ixRmax1,ixRmax2,ixRmax3,psb(igrid)%w,ixGlo1,ixGlo2,ixGlo3,&
                                ixGhi1,ixGhi2,ixGhi3,ixSmin1,ixSmin2,ixSmin3,ixSmax1,ixSmax2,&
                                ixSmax3,ipole)
                           if(stagger_grid) then
                              do idir=1,ndim
                                 ixSmin1=ixS_p_stg_min1(idir,inc1)
                                 ixSmin2=ixS_p_stg_min2(idir,inc2)
                                 ixSmin3=ixS_p_stg_min3(idir,inc3)
                                 ixSmax1=ixS_p_stg_max1(idir,inc1)
                                 ixSmax2=ixS_p_stg_max2(idir,inc2)
                                 ixSmax3=ixS_p_stg_max3(idir,inc3);
                                 ixRmin1=ixR_p_stg_min1(idir,n_inc1)
                                 ixRmin2=ixR_p_stg_min2(idir,n_inc2)
                                 ixRmin3=ixR_p_stg_min3(idir,n_inc3)
                                 ixRmax1=ixR_p_stg_max1(idir,n_inc1)
                                 ixRmax2=ixR_p_stg_max2(idir,n_inc2)
                                 ixRmax3=ixR_p_stg_max3(idir,n_inc3);
                                 call pole_copy_stg(psc(ineighbor)%ws,ixCoGsmin1,ixCoGsmin2,&
                                      ixCoGsmin3,ixCoGsmax1,ixCoGsmax2,ixCoGsmax3,ixRmin1,&
                                      ixRmin2,ixRmin3,ixRmax1,ixRmax2,ixRmax3,psb(igrid)%ws,&
                                      ixGslo1,ixGslo2,ixGslo3,ixGshi1,ixGshi2,ixGshi3,ixSmin1,&
                                      ixSmin2,ixSmin3,ixSmax1,ixSmax2,ixSmax3,idir,ipole)
                              end do
                           end if
                        end if
                     end do
                  end do
               end do
            end if
          end subroutine bc_fill_prolong

          subroutine gc_prolong(igrid)
            integer, intent(in) :: igrid

            integer :: iib1,iib2,iib3,i1,i2,i3,idims,iside
            logical,dimension(-1:1,-1:1,-1:1) :: NeedProlong

            iib1=idphyb(1,igrid);iib2=idphyb(2,igrid);iib3=idphyb(3,igrid);
            NeedProlong=.false.
            do i3=-1,1
               do i2=-1,1
                  do i1=-1,1
                     if (skip_direction([ i1,i2,i3 ])) cycle
                     if (neighbor_type(i1,i2,i3,igrid)==neighbor_coarse) then
                        call bc_prolong(igrid,i1,i2,i3,iib1,iib2,iib3)
                        NeedProlong(i1,i2,i3)=.true.
                     end if
                  end do
               end do
            end do
            if(stagger_grid) then
               ! Ghost cell prolongation for staggered variables
               ! must be done in a specific order.
               ! First the first neighbours, which have 2 indices=0 in 3D
               ! or one index=0 in 2D
               block=>psb(igrid)
                 do idims=1,ndim
                    i1=0;i2=0;i3=0;
                    select case(idims)
                    case(1)
                       do i1=-1,1,2
                          if (NeedProlong(i1,i2,i3)) call bc_prolong_stg(igrid,i1,i2,i3,&
                               iib1,iib2,iib3,NeedProlong)
                       end do

                    case(2)
                       do i2=-1,1,2
                          if (NeedProlong(i1,i2,i3)) call bc_prolong_stg(igrid,i1,i2,i3,&
                               iib1,iib2,iib3,NeedProlong)
                       end do

                    case(3)
                       do i3=-1,1,2
                          if (NeedProlong(i1,i2,i3)) call bc_prolong_stg(igrid,i1,i2,i3,&
                               iib1,iib2,iib3,NeedProlong)
                       end do

                    end select
                 end do
                 ! Then the second neighbours which have 1 index=0 in 3D
                 ! (Only in 3D)

                 i1=0;
                 do i2=-1,1,2
                    do i3=-1,1,2
                       if (NeedProlong(i1,i2,i3)) call bc_prolong_stg(igrid,i1,i2,i3,&
                            iib1,iib2,iib3,NeedProlong)
                    end do
                 end do
                 i2=0;
                 do i3=-1,1,2
                    do i1=-1,1,2
                       if (NeedProlong(i1,i2,i3)) call bc_prolong_stg(igrid,i1,i2,i3,&
                            iib1,iib2,iib3,NeedProlong)
                    end do
                 end do
                 i3=0;
                 do i1=-1,1,2
                    do i2=-1,1,2
                       if (NeedProlong(i1,i2,i3)) call bc_prolong_stg(igrid,i1,i2,i3,&
                            iib1,iib2,iib3,NeedProlong)
                    end do
                 end do

                 ! Finally, the corners, that have no index=0
                 do i1=-1,1,2
                    do i2=-1,1,2
                       do i3=-1,1,2
                          if (NeedProlong(i1,i2,i3)) call bc_prolong_stg(igrid,i1,i2,i3,iib1,&
                               iib2,iib3,NeedProlong)
                       end do
                    end do
                 end do
              end if
            end subroutine gc_prolong

            !> fill coarser representative with data from coarser neighbors
            subroutine bc_fill_prolong_stg
              ic1=1+modulo(node(pig1_,igrid)-1,2)
              ic2=1+modulo(node(pig2_,igrid)-1,2)
              ic3=1+modulo(node(pig3_,igrid)-1,2);
              if (.not.(i1==0.or.i1==2*ic1-3).or..not.(i2==0.or.i2==&
                   2*ic2-3).or..not.(i3==0.or.i3==2*ic3-3)) return

              ipe_neighbor=neighbor(2,i1,i2,i3,igrid)
              if(ipe_neighbor/=mype) then
                 ineighbor=neighbor(1,i1,i2,i3,igrid)
                 ipole=neighbor_pole(i1,i2,i3,igrid)

                 if (ipole==0) then   !! There is no pole
                    inc1=ic1+i1;inc2=ic2+i2;inc3=ic3+i3;
                    ixRmin1=ixR_p_min1(iib1,inc1);ixRmin2=ixR_p_min2(iib2,inc2)
                    ixRmin3=ixR_p_min3(iib3,inc3);ixRmax1=ixR_p_max1(iib1,inc1)
                    ixRmax2=ixR_p_max2(iib2,inc2);ixRmax3=ixR_p_max3(iib3,inc3);
                    do idir=1,ndim
                       ixRmin1=ixR_p_stg_min1(idir,inc1)
                       ixRmin2=ixR_p_stg_min2(idir,inc2)
                       ixRmin3=ixR_p_stg_min3(idir,inc3)
                       ixRmax1=ixR_p_stg_max1(idir,inc1)
                       ixRmax2=ixR_p_stg_max2(idir,inc2)
                       ixRmax3=ixR_p_stg_max3(idir,inc3);
                       ibuf_next=ibuf_recv_p+sizes_p_recv_stg(idir,inc1,inc2,inc3)
                       psc(igrid)%ws(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,&
                            idir)=reshape(source=recvbuffer_p(ibuf_recv_p:ibuf_next-1),&
                            shape=shape(psc(igrid)%ws(ixRmin1:ixRmax1,ixRmin2:ixRmax2,&
                            ixRmin3:ixRmax3,idir)))
                       ibuf_recv_p=ibuf_next
                    end do
                 else !! There is a pole
                    inc1=ic1+i1;inc2=ic2+i2;inc3=ic3+i3;
                    select case (ipole)
                    case (1)
                       n_inc1=2*i1+(3-ic1);n_inc2=-2*i2+ic2;n_inc3=-2*i3+ic3;
                    case (2)
                       n_inc1=-2*i1+ic1;n_inc2=2*i2+(3-ic2);n_inc3=-2*i3+ic3;
                    case (3)
                       n_inc1=-2*i1+ic1;n_inc2=-2*i2+ic2;n_inc3=2*i3+(3-ic3);
                    end select
                    !! Unpack the buffer and fill an auxiliary array
                    pole_buf%ws=zero
                    do idir=1,ndim
                       ixRmin1=ixR_p_stg_min1(idir,inc1)
                       ixRmin2=ixR_p_stg_min2(idir,inc2)
                       ixRmin3=ixR_p_stg_min3(idir,inc3)
                       ixRmax1=ixR_p_stg_max1(idir,inc1)
                       ixRmax2=ixR_p_stg_max2(idir,inc2)
                       ixRmax3=ixR_p_stg_max3(idir,inc3);
                       ibuf_next=ibuf_recv_p+sizes_p_recv_stg(idir,inc1,inc2,inc3)
                       pole_buf%ws(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,&
                            idir)=reshape(source=recvbuffer_p(ibuf_recv_p:ibuf_next-1),&
                            shape=shape(psc(igrid)%ws(ixRmin1:ixRmax1,ixRmin2:ixRmax2,&
                            ixRmin3:ixRmax3,idir)))
                       call pole_copy_stg(psc(igrid)%ws,ixCoGsmin1,ixCoGsmin2,&
                            ixCoGsmin3,ixCoGsmax1,ixCoGsmax2,ixCoGsmax3,ixRmin1,ixRmin2,&
                            ixRmin3,ixRmax1,ixRmax2,ixRmax3,pole_buf%ws,ixGslo1,ixGslo2,&
                            ixGslo3,ixGshi1,ixGshi2,ixGshi3,ixRmin1,ixRmin2,ixRmin3,&
                            ixRmax1,ixRmax2,ixRmax3,idir,ipole)
                       ibuf_recv_p=ibuf_next
                    end do
                 end if
              end if

            end subroutine bc_fill_prolong_stg

            !> do prolongation for fine blocks after receipt data from coarse neighbors
            subroutine bc_prolong(igrid,i1,i2,i3,iib1,iib2,iib3)
              use mod_physics, only: phys_to_primitive, phys_to_conserved

              integer :: i1,i2,i3,iib1,iib2,iib3,igrid
              integer :: ixFimin1,ixFimin2,ixFimin3,ixFimax1,ixFimax2,ixFimax3,&
                   ixComin1,ixComin2,ixComin3,ixComax1,ixComax2,ixComax3,ii1,ii2,ii3,&
                   idims,iside,ixBmin1,ixBmin2,ixBmin3,ixBmax1,ixBmax2,ixBmax3
              double precision :: dxFi1,dxFi2,dxFi3, dxCo1,dxCo2,dxCo3, xFimin1,&
                   xFimin2,xFimin3, xComin1,xComin2,xComin3, invdxCo1,invdxCo2,&
                   invdxCo3

              ixFimin1=ixR_srl_min1(iib1,i1);ixFimin2=ixR_srl_min2(iib2,i2)
              ixFimin3=ixR_srl_min3(iib3,i3);ixFimax1=ixR_srl_max1(iib1,i1)
              ixFimax2=ixR_srl_max2(iib2,i2);ixFimax3=ixR_srl_max3(iib3,i3);
              dxFi1=rnode(rpdx1_,igrid);dxFi2=rnode(rpdx2_,igrid)
              dxFi3=rnode(rpdx3_,igrid);
              dxCo1=two*dxFi1;dxCo2=two*dxFi2;dxCo3=two*dxFi3;
              invdxCo1=1.d0/dxCo1;invdxCo2=1.d0/dxCo2;invdxCo3=1.d0/dxCo3;

              ! compute the enlarged grid lower left corner coordinates
              ! these are true coordinates for an equidistant grid,
              ! but we can temporarily also use them for getting indices
              ! in stretched grids
              xFimin1=rnode(rpxmin1_,igrid)-dble(nghostcells)*dxFi1
              xFimin2=rnode(rpxmin2_,igrid)-dble(nghostcells)*dxFi2
              xFimin3=rnode(rpxmin3_,igrid)-dble(nghostcells)*dxFi3;
              xComin1=rnode(rpxmin1_,igrid)-dble(nghostcells)*dxCo1
              xComin2=rnode(rpxmin2_,igrid)-dble(nghostcells)*dxCo2
              xComin3=rnode(rpxmin3_,igrid)-dble(nghostcells)*dxCo3;

              if(stagger_grid.and.phyboundblock(igrid).and.bcphys) then
                 block=>psc(igrid)
                   do idims=1,ndim
                      ixComin1=int((xFimin1+(dble(ixFimin1)-half)*dxFi1-&
                           xComin1)*invdxCo1)+1-1
                      ixComin2=int((xFimin2+(dble(ixFimin2)-half)*dxFi2-&
                           xComin2)*invdxCo2)+1-1
                      ixComin3=int((xFimin3+(dble(ixFimin3)-half)*dxFi3-&
                           xComin3)*invdxCo3)+1-1;
                      ixComax1=int((xFimin1+(dble(ixFimax1)-half)*dxFi1-&
                           xComin1)*invdxCo1)+1+1
                      ixComax2=int((xFimin2+(dble(ixFimax2)-half)*dxFi2-&
                           xComin2)*invdxCo2)+1+1
                      ixComax3=int((xFimin3+(dble(ixFimax3)-half)*dxFi3-&
                           xComin3)*invdxCo3)+1+1;

                      ! avoid using undetermined ghost cells at physical boundary edges
                      if(idims == 1) then
                         if(neighbor_type(-1,0,0,igrid)==neighbor_boundary .or. &
                              neighbor_type(1,0,0,igrid)==neighbor_boundary) then
                            if(neighbor_type(0,-1,0,igrid)==neighbor_boundary) &
                                 ixComin2=ixCoMmin2
                            if(neighbor_type(0,0,-1,igrid)==neighbor_boundary) &
                                 ixComin3=ixCoMmin3
                            if(neighbor_type(0,1,0,igrid)==neighbor_boundary) &
                                 ixComax2=ixCoMmax2
                            if(neighbor_type(0,0,1,igrid)==neighbor_boundary) &
                                 ixComax3=ixCoMmax3
                         end if
                      else if(idims == 2) then
                         if(neighbor_type(0,-1,0,igrid)==neighbor_boundary .or. &
                              neighbor_type(0,1,0,igrid)==neighbor_boundary) then
                            if(neighbor_type(0,0,-1,igrid)==neighbor_boundary) &
                                 ixComin3=ixCoMmin3
                            if(neighbor_type(0,0,1,igrid)==neighbor_boundary) &
                                 ixComax3=ixCoMmax3
                         end if
                      end if

                      do iside=1,2
                         ii1=kr(1,idims)*(2*iside-3);ii2=kr(2,idims)*(2*iside-3)
                         ii3=kr(3,idims)*(2*iside-3);
                         if(neighbor_type(ii1,ii2,ii3,igrid)/=neighbor_boundary) cycle
                         if(( (iside==1.and.idims==1.and.ixComin1<ixCoGmin1+&
                              nghostcells).or.(iside==1.and.idims==&
                              2.and.ixComin2<ixCoGmin2+nghostcells).or.(iside==1.and.idims==&
                              3.and.ixComin3<ixCoGmin3+nghostcells) ) .or.( &
                              (iside==2.and.idims==1.and.ixComax1>ixCoGmax1-&
                              nghostcells).or. (iside==2.and.idims==&
                              2.and.ixComax2>ixCoGmax2-nghostcells).or. &
                              (iside==2.and.idims==3.and.ixComax3>ixCoGmax3-nghostcells))) &
                              then
                            ixBmin1=merge(ixCoGmin1,ixComin1,idims==1)
                            ixBmin2=merge(ixCoGmin2,ixComin2,idims==2)
                            ixBmin3=merge(ixCoGmin3,ixComin3,idims==3);
                            ixBmax1=merge(ixCoGmax1,ixComax1,idims==1)
                            ixBmax2=merge(ixCoGmax2,ixComax2,idims==2)
                            ixBmax3=merge(ixCoGmax3,ixComax3,idims==3);
                            call bc_phys(iside,idims,time,0.d0,psc(igrid),ixCoGmin1,&
                                 ixCoGmin2,ixCoGmin3,ixCoGmax1,ixCoGmax2,ixCoGmax3,ixBmin1,&
                                 ixBmin2,ixBmin3,ixBmax1,ixBmax2,ixBmax3)
                         end if
                      end do
                   end do
                end if

                if(prolongprimitive) then
                   ! following line again assumes equidistant grid, but
                   ! just computes indices, so also ok for stretched case
                   ! reason for +1-1 and +1+1: the coarse representation has
                   ! also nghostcells at each side. During
                   ! prolongation, we need cells to left and right, hence -1/+1
                   block=>psc(igrid)
                     ixComin1=int((xFimin1+(dble(ixFimin1)-half)*dxFi1-xComin1)*invdxCo1)+&
                          1-1
                     ixComin2=int((xFimin2+(dble(ixFimin2)-half)*dxFi2-xComin2)*invdxCo2)+&
                          1-1
                     ixComin3=int((xFimin3+(dble(ixFimin3)-half)*dxFi3-xComin3)*invdxCo3)+&
                          1-1;
                     ixComax1=int((xFimin1+(dble(ixFimax1)-half)*dxFi1-xComin1)*invdxCo1)+&
                          1+1
                     ixComax2=int((xFimin2+(dble(ixFimax2)-half)*dxFi2-xComin2)*invdxCo2)+&
                          1+1
                     ixComax3=int((xFimin3+(dble(ixFimax3)-half)*dxFi3-xComin3)*invdxCo3)+&
                          1+1;
                     call phys_to_primitive(ixCoGmin1,ixCoGmin2,ixCoGmin3,ixCoGmax1,&
                          ixCoGmax2,ixCoGmax3,ixComin1,ixComin2,ixComin3,ixComax1,ixComax2,&
                          ixComax3,psc(igrid)%w,psc(igrid)%x)
                  end if

                  if(ghost_copy) then
                     call interpolation_copy(igrid,ixFimin1,ixFimin2,ixFimin3,ixFimax1,&
                          ixFimax2,ixFimax3,dxFi1,dxFi2,dxFi3,xFimin1,xFimin2,xFimin3,dxCo1,&
                          dxCo2,dxCo3,invdxCo1,invdxCo2,invdxCo3,xComin1,xComin2,xComin3)
                  else
                     call interpolation_linear(igrid,ixFimin1,ixFimin2,ixFimin3,ixFimax1,&
                          ixFimax2,ixFimax3,dxFi1,dxFi2,dxFi3,xFimin1,xFimin2,xFimin3,dxCo1,&
                          dxCo2,dxCo3,invdxCo1,invdxCo2,invdxCo3,xComin1,xComin2,xComin3)
                  end if

                  if(prolongprimitive) then
                     block=>psc(igrid)
                       call phys_to_conserved(ixCoGmin1,ixCoGmin2,ixCoGmin3,ixCoGmax1,&
                            ixCoGmax2,ixCoGmax3,ixComin1,ixComin2,ixComin3,ixComax1,ixComax2,&
                            ixComax3,psc(igrid)%w,psc(igrid)%x)
                    end if

                  end subroutine bc_prolong

                  subroutine bc_prolong_stg(igrid,i1,i2,i3,iib1,iib2,iib3,NeedProlong)
                    use mod_amr_fct
                    integer                    :: igrid,i1,i2,i3,iib1,iib2,iib3
                    logical,dimension(-1:1,-1:1,-1:1) :: NeedProlong
                    logical                    :: fine_min1in,fine_min2in,fine_min3in,&
                         fine_max1in,fine_max2in,fine_max3in
                    integer                    :: ixFimin1,ixFimin2,ixFimin3,ixFimax1,&
                         ixFimax2,ixFimax3,ixComin1,ixComin2,ixComin3,ixComax1,ixComax2,&
                         ixComax3
                    double precision           :: dxFi1,dxFi2,dxFi3,dxCo1,dxCo2,dxCo3,&
                         xFimin1,xFimin2,xFimin3,xComin1,xComin2,xComin3,invdxCo1,invdxCo2,&
                         invdxCo3
                    ! Check what is already at the desired level
                    fine_min1in=.false.;fine_min2in=.false.;fine_min3in=.false.
                    fine_max1in=.false.;fine_max2in=.false.;fine_max3in=.false.;

                    if(i1>-1) fine_min1in=(.not.NeedProlong(i1-kr(1,1),i2-kr(1,2),i3-kr(1,&
                         3)).and.neighbor_type(i1-kr(1,1),i2-kr(1,2),i3-kr(1,3),igrid)/=1)
                    if(i1<1)  fine_max1in=(.not.NeedProlong(i1+kr(1,1),i2+kr(1,2),i3+kr(1,&
                         3)).and.neighbor_type(i1+kr(1,1),i2+kr(1,2),i3+kr(1,3),igrid)/=1)


                    if(i2>-1) fine_min2in=(.not.NeedProlong(i1-kr(2,1),i2-kr(2,2),i3-kr(2,&
                         3)).and.neighbor_type(i1-kr(2,1),i2-kr(2,2),i3-kr(2,3),igrid)/=1)
                    if(i2<1)  fine_max2in=(.not.NeedProlong(i1+kr(2,1),i2+kr(2,2),i3+kr(2,&
                         3)).and.neighbor_type(i1+kr(2,1),i2+kr(2,2),i3+kr(2,3),igrid)/=1)


                    if(i3>-1) fine_min3in=(.not.NeedProlong(i1-kr(3,1),i2-kr(3,2),i3-kr(3,&
                         3)).and.neighbor_type(i1-kr(3,1),i2-kr(3,2),i3-kr(3,3),igrid)/=1)
                    if(i3<1)  fine_max3in=(.not.NeedProlong(i1+kr(3,1),i2+kr(3,2),i3+kr(3,&
                         3)).and.neighbor_type(i1+kr(3,1),i2+kr(3,2),i3+kr(3,3),igrid)/=1)


                    ixFimin1=ixR_srl_min1(iib1,i1);ixFimin2=ixR_srl_min2(iib2,i2)
                    ixFimin3=ixR_srl_min3(iib3,i3);ixFimax1=ixR_srl_max1(iib1,i1)
                    ixFimax2=ixR_srl_max2(iib2,i2);ixFimax3=ixR_srl_max3(iib3,i3);

                    dxFi1=rnode(rpdx1_,igrid);dxFi2=rnode(rpdx2_,igrid)
                    dxFi3=rnode(rpdx3_,igrid);
                    dxCo1=two*dxFi1;dxCo2=two*dxFi2;dxCo3=two*dxFi3;
                    invdxCo1=1.d0/dxCo1;invdxCo2=1.d0/dxCo2;invdxCo3=1.d0/dxCo3;

                    xFimin1=rnode(rpxmin1_,igrid)-dble(nghostcells)*dxFi1
                    xFimin2=rnode(rpxmin2_,igrid)-dble(nghostcells)*dxFi2
                    xFimin3=rnode(rpxmin3_,igrid)-dble(nghostcells)*dxFi3;
                    xComin1=rnode(rpxmin1_,igrid)-dble(nghostcells)*dxCo1
                    xComin2=rnode(rpxmin2_,igrid)-dble(nghostcells)*dxCo2
                    xComin3=rnode(rpxmin3_,igrid)-dble(nghostcells)*dxCo3;

                    ! moved the physical boundary filling here, to only fill the
                    ! part needed

                    ixComin1=int((xFimin1+(dble(ixFimin1)-half)*dxFi1-xComin1)*invdxCo1)+&
                         1-1
                    ixComin2=int((xFimin2+(dble(ixFimin2)-half)*dxFi2-xComin2)*invdxCo2)+&
                         1-1
                    ixComin3=int((xFimin3+(dble(ixFimin3)-half)*dxFi3-xComin3)*invdxCo3)+&
                         1-1;
                    ixComax1=int((xFimin1+(dble(ixFimax1)-half)*dxFi1-xComin1)*invdxCo1)+1+&
                         1
                    ixComax2=int((xFimin2+(dble(ixFimax2)-half)*dxFi2-xComin2)*invdxCo2)+1+&
                         1
                    ixComax3=int((xFimin3+(dble(ixFimax3)-half)*dxFi3-xComin3)*invdxCo3)+1+&
                         1;

                    if(prolongprimitive) call phys_to_primitive(ixGlo1,ixGlo2,ixGlo3,&
                         ixGhi1,ixGhi2,ixGhi3,ixFimin1,ixFimin2,ixFimin3,ixFimax1,ixFimax2,&
                         ixFimax3,psb(igrid)%w,psb(igrid)%x)

                    call prolong_2nd_stg(psc(igrid),psb(igrid),ixComin1,ixComin2,ixComin3,&
                         ixComax1,ixComax2,ixComax3,ixFimin1,ixFimin2,ixFimin3,ixFimax1,&
                         ixFimax2,ixFimax3,dxCo1,dxCo2,dxCo3,xComin1,xComin2,xComin3,dxFi1,&
                         dxFi2,dxFi3,xFimin1,xFimin2,xFimin3,.true.,fine_min1in,fine_min2in,&
                         fine_min3in,fine_max1in,fine_max2in,fine_max3in)

                    if(prolongprimitive) call phys_to_conserved(ixGlo1,ixGlo2,ixGlo3,&
                         ixGhi1,ixGhi2,ixGhi3,ixFimin1,ixFimin2,ixFimin3,ixFimax1,ixFimax2,&
                         ixFimax3,psb(igrid)%w,psb(igrid)%x)

                    ! The current region has already been refined, so it does not need to be prolonged again
                    NeedProlong(i1,i2,i3)=.false.

                  end subroutine bc_prolong_stg

                  subroutine interpolation_linear(igrid,ixFimin1,ixFimin2,ixFimin3,&
                       ixFimax1,ixFimax2,ixFimax3,dxFi1,dxFi2,dxFi3,xFimin1,xFimin2,xFimin3,&
                       dxCo1,dxCo2,dxCo3,invdxCo1,invdxCo2,invdxCo3,xComin1,xComin2,&
                       xComin3)
                    use mod_physics, only: phys_to_conserved
                    integer, intent(in) :: igrid, ixFimin1,ixFimin2,ixFimin3,ixFimax1,&
                         ixFimax2,ixFimax3
                    double precision, intent(in) :: dxFi1,dxFi2,dxFi3, xFimin1,xFimin2,&
                         xFimin3,dxCo1,dxCo2,dxCo3, invdxCo1,invdxCo2,invdxCo3, xComin1,&
                         xComin2,xComin3

                    integer :: ixCo1,ixCo2,ixCo3, jxCo1,jxCo2,jxCo3, hxCo1,hxCo2,hxCo3,&
                         ixFi1,ixFi2,ixFi3, ix1,ix2,ix3, iw, idims, nwmin,nwmax
                    double precision :: xCo1,xCo2,xCo3, xFi1,xFi2,xFi3, eta1,eta2,eta3
                    double precision :: slopeL, slopeR, slopeC, signC, signR
                    double precision :: slope(1:nw,ndim)
                    !!double precision :: local_invdxCo^D
                    double precision :: signedfactorhalf1,signedfactorhalf2,&
                         signedfactorhalf3
                    !integer :: ixshift^D, icase

                    !icase=mod(nghostcells,2)

                    if(prolongprimitive) then
                       nwmin=1
                       nwmax=nw
                    else
                       nwmin=nwhead
                       nwmax=nwtail
                    end if

                    do ixFi3 = ixFimin3,ixFimax3
                       ! cell-centered coordinates of fine grid point
                       ! here we temporarily use an equidistant grid
                       xFi3=xFimin3+(dble(ixFi3)-half)*dxFi3

                       ! indices of coarse cell which contains the fine cell
                       ! since we computed lower left corner earlier
                       ! in equidistant fashion: also ok for stretched case
                       ixCo3=int((xFi3-xComin3)*invdxCo3)+1

                       ! cell-centered coordinates of coarse grid point
                       ! here we temporarily use an equidistant grid
                       xCo3=xComin3+(dble(ixCo3)-half)*dxCo3 
                       do ixFi2 = ixFimin2,ixFimax2
                          ! cell-centered coordinates of fine grid point
                          ! here we temporarily use an equidistant grid
                          xFi2=xFimin2+(dble(ixFi2)-half)*dxFi2

                          ! indices of coarse cell which contains the fine cell
                          ! since we computed lower left corner earlier
                          ! in equidistant fashion: also ok for stretched case
                          ixCo2=int((xFi2-xComin2)*invdxCo2)+1

                          ! cell-centered coordinates of coarse grid point
                          ! here we temporarily use an equidistant grid
                          xCo2=xComin2+(dble(ixCo2)-half)*dxCo2 
                          do ixFi1 = ixFimin1,ixFimax1
                             ! cell-centered coordinates of fine grid point
                             ! here we temporarily use an equidistant grid
                             xFi1=xFimin1+(dble(ixFi1)-half)*dxFi1

                             ! indices of coarse cell which contains the fine cell
                             ! since we computed lower left corner earlier
                             ! in equidistant fashion: also ok for stretched case
                             ixCo1=int((xFi1-xComin1)*invdxCo1)+1

                             ! cell-centered coordinates of coarse grid point
                             ! here we temporarily use an equidistant grid
                             xCo1=xComin1+(dble(ixCo1)-half)*dxCo1 

                             !if(.not.slab) then
                             !   ^D&local_invdxCo^D=1.d0/psc(igrid)%dx({ixCo^DD},^D)\
                             !endif

                             if(slab_uniform) then
                                ! actual cell-centered coordinates of fine grid point
                                !!^D&xFi^D=block%x({ixFi^DD},^D)\
                                ! actual cell-centered coordinates of coarse grid point
                                !!^D&xCo^D=psc(igrid)%x({ixCo^DD},^D)\
                                ! normalized distance between fine/coarse cell center
                                ! in coarse cell: ranges from -0.5 to 0.5 in each direction
                                ! (origin is coarse cell center)
                                ! this is essentially +1/4 or -1/4 on cartesian mesh
                                eta1=(xFi1-xCo1)*invdxCo1;eta2=(xFi2-xCo2)*invdxCo2
                                eta3=(xFi3-xCo3)*invdxCo3;
                             else
                                !select case(icase)
                                ! case(0)
                                !{! here we assume an even number of ghostcells!!!
                                !ixshift^D=2*(mod(ixFi^D,2)-1)+1
                                !if(ixshift^D>0.0d0)then
                                !   ! oneven fine grid points
                                !   eta^D=-0.5d0*(one-block%dvolume(ixFi^DD) &
                                !     /sum(block%dvolume(ixFi^D:ixFi^D+1^D%ixFi^DD)))
                                !else
                                !   ! even fine grid points
                                !   eta^D=+0.5d0*(one-block%dvolume(ixFi^DD) &
                                !     /sum(block%dvolume(ixFi^D-1:ixFi^D^D%ixFi^DD)))
                                !endif\}
                                ! case(1)
                                !{! here we assume an odd number of ghostcells!!!
                                !ixshift^D=2*(mod(ixFi^D,2)-1)+1
                                !if(ixshift^D>0.0d0)then
                                !   ! oneven fine grid points
                                !   eta^D=+0.5d0*(one-block%dvolume(ixFi^DD) &
                                !     /sum(block%dvolume(ixFi^D-1:ixFi^D^D%ixFi^DD)))
                                !else
                                !   ! even fine grid points
                                !   eta^D=-0.5d0*(one-block%dvolume(ixFi^DD) &
                                !     /sum(block%dvolume(ixFi^D:ixFi^D+1^D%ixFi^DD)))
                                !endif\}
                                ! case default
                                !  call mpistop("no such case")
                                !end select
                                ! the different cases for even/uneven number of ghost cells
                                ! are automatically handled using the relative index to ixMlo
                                ! as well as the pseudo-coordinates xFi and xCo
                                ! these latter differ from actual cell centers when stretching is used
                                ix1=2*int((ixFi1+ixMlo1)/2)-ixMlo1
                                ix2=2*int((ixFi2+ixMlo2)/2)-ixMlo2
                                ix3=2*int((ixFi3+ixMlo3)/2)-ixMlo3;
                                if(xFi1>xCo1) then
                                   signedfactorhalf1=0.5d0
                                else
                                   signedfactorhalf1=-0.5d0
                                end if
                                eta1=signedfactorhalf1*(one-psb(igrid)%dvolume(ixFi1,ixFi2,&
                                     ixFi3) /sum(psb(igrid)%dvolume(ix1:ix1+1,ixFi2,ixFi3))) 
                                if(xFi2>xCo2) then
                                   signedfactorhalf2=0.5d0
                                else
                                   signedfactorhalf2=-0.5d0
                                end if
                                eta2=signedfactorhalf2*(one-psb(igrid)%dvolume(ixFi1,ixFi2,&
                                     ixFi3) /sum(psb(igrid)%dvolume(ixFi1,ix2:ix2+1,ixFi3))) 
                                if(xFi3>xCo3) then
                                   signedfactorhalf3=0.5d0
                                else
                                   signedfactorhalf3=-0.5d0
                                end if
                                eta3=signedfactorhalf3*(one-psb(igrid)%dvolume(ixFi1,ixFi2,&
                                     ixFi3) /sum(psb(igrid)%dvolume(ixFi1,ixFi2,ix3:ix3+1))) 
                                !{eta^D=(xFi^D-xCo^D)*invdxCo^D &
                                !      *two*(one-block%dvolume(ixFi^DD) &
                                !      /sum(block%dvolume(ix^D:ix^D+1^D%ixFi^DD))) \}
                             end if

                             do idims=1,ndim
                                hxCo1=ixCo1-kr(1,idims)
                                hxCo2=ixCo2-kr(2,idims)
                                hxCo3=ixCo3-kr(3,idims)
                                jxCo1=ixCo1+kr(1,idims)
                                jxCo2=ixCo2+kr(2,idims)
                                jxCo3=ixCo3+kr(3,idims)

                                do iw=nwmin,nwmax
                                   slopeL=psc(igrid)%w(ixCo1,ixCo2,ixCo3,iw)-psc(igrid)%w(hxCo1,&
                                        hxCo2,hxCo3,iw)
                                   slopeR=psc(igrid)%w(jxCo1,jxCo2,jxCo3,iw)-psc(igrid)%w(ixCo1,&
                                        ixCo2,ixCo3,iw)
                                   slopeC=half*(slopeR+slopeL)

                                   ! get limited slope
                                   signR=sign(one,slopeR)
                                   signC=sign(one,slopeC)
                                   !select case(prolong_limiter)
                                   !case(1)
                                   !  ! unlimit
                                   !  slope(iw,idims)=slopeC
                                   !case(2)
                                   !  ! minmod
                                   !  slope(iw,idims)=signR*max(zero,min(dabs(slopeR), &
                                   !                                    signR*slopeL))
                                   !case(3)
                                   !  ! woodward
                                   !  slope(iw,idims)=two*signR*max(zero,min(dabs(slopeR), &
                                   !                     signR*slopeL,signR*half*slopeC))
                                   !case(4)
                                   !  ! koren
                                   !  slope(iw,idims)=signR*max(zero,min(two*signR*slopeL, &
                                   !   (dabs(slopeR)+two*slopeL*signR)*third,two*dabs(slopeR)))
                                   !case default
                                   slope(iw,idims)=signC*max(zero,min(dabs(slopeC),&
                                        signC*slopeL,signC*slopeR))
                                   !end select
                                end do
                             end do

                             ! Interpolate from coarse cell using limited slopes
                             psb(igrid)%w(ixFi1,ixFi2,ixFi3,nwmin:nwmax)=psc(igrid)%w(ixCo1,&
                                  ixCo2,ixCo3,nwmin:nwmax)+(slope(nwmin:nwmax,&
                                  1)*eta1)+(slope(nwmin:nwmax,2)*eta2)+(slope(nwmin:nwmax,3)*eta3)

                          end do
                       end do
                    end do

                    if(prolongprimitive) then
                       block=>psb(igrid)
                         call phys_to_conserved(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,&
                              ixFimin1,ixFimin2,ixFimin3,ixFimax1,ixFimax2,ixFimax3,&
                              psb(igrid)%w,psb(igrid)%x)
                      end if

                    end subroutine interpolation_linear

                    subroutine interpolation_copy(igrid, ixFimin1,ixFimin2,ixFimin3,ixFimax1,&
                         ixFimax2,ixFimax3,dxFi1,dxFi2,dxFi3,xFimin1,xFimin2,xFimin3, dxCo1,&
                         dxCo2,dxCo3,invdxCo1,invdxCo2,invdxCo3,xComin1,xComin2,xComin3)
                      use mod_physics, only: phys_to_conserved
                      integer, intent(in) :: igrid, ixFimin1,ixFimin2,ixFimin3,ixFimax1,&
                           ixFimax2,ixFimax3
                      double precision, intent(in) :: dxFi1,dxFi2,dxFi3, xFimin1,xFimin2,&
                           xFimin3,dxCo1,dxCo2,dxCo3, invdxCo1,invdxCo2,invdxCo3, xComin1,&
                           xComin2,xComin3

                      integer :: ixCo1,ixCo2,ixCo3, ixFi1,ixFi2,ixFi3, nwmin,nwmax
                      double precision :: xFi1,xFi2,xFi3

                      if(prolongprimitive) then
                         nwmin=1
                         nwmax=nw
                      else
                         nwmin=nwhead
                         nwmax=nwtail
                      end if

                      do ixFi3 = ixFimin3,ixFimax3
                         ! cell-centered coordinates of fine grid point
                         xFi3=xFimin3+(dble(ixFi3)-half)*dxFi3

                         ! indices of coarse cell which contains the fine cell
                         ! note: this also works for stretched grids
                         ixCo3=int((xFi3-xComin3)*invdxCo3)+1
                         do ixFi2 = ixFimin2,ixFimax2
                            ! cell-centered coordinates of fine grid point
                            xFi2=xFimin2+(dble(ixFi2)-half)*dxFi2

                            ! indices of coarse cell which contains the fine cell
                            ! note: this also works for stretched grids
                            ixCo2=int((xFi2-xComin2)*invdxCo2)+1
                            do ixFi1 = ixFimin1,ixFimax1
                               ! cell-centered coordinates of fine grid point
                               xFi1=xFimin1+(dble(ixFi1)-half)*dxFi1

                               ! indices of coarse cell which contains the fine cell
                               ! note: this also works for stretched grids
                               ixCo1=int((xFi1-xComin1)*invdxCo1)+1

                               ! Copy from coarse cell
                               psb(igrid)%w(ixFi1,ixFi2,ixFi3,nwmin:nwmax)=psc(igrid)%w(ixCo1,&
                                    ixCo2,ixCo3,nwmin:nwmax)

                            end do
                         end do
                      end do

                      if(prolongprimitive) call phys_to_conserved(ixGlo1,ixGlo2,ixGlo3,&
                           ixGhi1,ixGhi2,ixGhi3,ixFimin1,ixFimin2,ixFimin3,ixFimax1,ixFimax2,&
                           ixFimax3,psb(igrid)%w,psb(igrid)%x)

                    end subroutine interpolation_copy

                    subroutine pole_copy(wrecv,ixIRmin1,ixIRmin2,ixIRmin3,ixIRmax1,ixIRmax2,&
                         ixIRmax3,ixRmin1,ixRmin2,ixRmin3,ixRmax1,ixRmax2,ixRmax3,wsend,&
                         ixISmin1,ixISmin2,ixISmin3,ixISmax1,ixISmax2,ixISmax3,ixSmin1,ixSmin2,&
                         ixSmin3,ixSmax1,ixSmax2,ixSmax3,ipole)

                      integer, intent(in) :: ixIRmin1,ixIRmin2,ixIRmin3,ixIRmax1,ixIRmax2,&
                           ixIRmax3,ixRmin1,ixRmin2,ixRmin3,ixRmax1,ixRmax2,ixRmax3,ixISmin1,&
                           ixISmin2,ixISmin3,ixISmax1,ixISmax2,ixISmax3,ixSmin1,ixSmin2,&
                           ixSmin3,ixSmax1,ixSmax2,ixSmax3,ipole
                      double precision :: wrecv(ixIRmin1:ixIRmax1,ixIRmin2:ixIRmax2,&
                           ixIRmin3:ixIRmax3,1:nw), wsend(ixISmin1:ixISmax1,ixISmin2:ixISmax2,&
                           ixISmin3:ixISmax3,1:nw)

                      integer :: iw, iside, iB

                      select case (ipole)
                      case (1)
                         iside=int((i1+3)/2)
                         iB=2*(1-1)+iside
                         do iw=nwhead,nwtail
                            select case (typeboundary(iw,iB))
                            case (bc_symm)
                               wrecv(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,&
                                    iw) = wsend(ixSmax1:ixSmin1:-1,ixSmin2:ixSmax2,&
                                    ixSmin3:ixSmax3,iw)
                            case (bc_asymm)
                               wrecv(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,&
                                    iw) =-wsend(ixSmax1:ixSmin1:-1,ixSmin2:ixSmax2,&
                                    ixSmin3:ixSmax3,iw)
                            case default
                               call mpistop("Pole boundary condition should be symm or asymm")
                            end select
                         end do
                      case (2)
                         iside=int((i2+3)/2)
                         iB=2*(2-1)+iside
                         do iw=nwhead,nwtail
                            select case (typeboundary(iw,iB))
                            case (bc_symm)
                               wrecv(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,&
                                    iw) = wsend(ixSmin1:ixSmax1,ixSmax2:ixSmin2:-1,&
                                    ixSmin3:ixSmax3,iw)
                            case (bc_asymm)
                               wrecv(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,&
                                    iw) =-wsend(ixSmin1:ixSmax1,ixSmax2:ixSmin2:-1,&
                                    ixSmin3:ixSmax3,iw)
                            case default
                               call mpistop("Pole boundary condition should be symm or asymm")
                            end select
                         end do
                      case (3)
                         iside=int((i3+3)/2)
                         iB=2*(3-1)+iside
                         do iw=nwhead,nwtail
                            select case (typeboundary(iw,iB))
                            case (bc_symm)
                               wrecv(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,&
                                    iw) = wsend(ixSmin1:ixSmax1,ixSmin2:ixSmax2,&
                                    ixSmax3:ixSmin3:-1,iw)
                            case (bc_asymm)
                               wrecv(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,&
                                    iw) =-wsend(ixSmin1:ixSmax1,ixSmin2:ixSmax2,&
                                    ixSmax3:ixSmin3:-1,iw)
                            case default
                               call mpistop("Pole boundary condition should be symm or asymm")
                            end select
                         end do
                      end select

                    end subroutine pole_copy

                    subroutine pole_copy_stg(wrecv,ixIRmin1,ixIRmin2,ixIRmin3,ixIRmax1,&
                         ixIRmax2,ixIRmax3,ixRmin1,ixRmin2,ixRmin3,ixRmax1,ixRmax2,ixRmax3,&
                         wsend,ixISmin1,ixISmin2,ixISmin3,ixISmax1,ixISmax2,ixISmax3,ixSmin1,&
                         ixSmin2,ixSmin3,ixSmax1,ixSmax2,ixSmax3,idirs,ipole)

                      integer, intent(in) :: ixIRmin1,ixIRmin2,ixIRmin3,ixIRmax1,ixIRmax2,&
                           ixIRmax3,ixRmin1,ixRmin2,ixRmin3,ixRmax1,ixRmax2,ixRmax3,ixISmin1,&
                           ixISmin2,ixISmin3,ixISmax1,ixISmax2,ixISmax3,ixSmin1,ixSmin2,&
                           ixSmin3,ixSmax1,ixSmax2,ixSmax3,idirs,ipole

                      double precision :: wrecv(ixIRmin1:ixIRmax1,ixIRmin2:ixIRmax2,&
                           ixIRmin3:ixIRmax3,1:nws), wsend(ixISmin1:ixISmax1,ixISmin2:ixISmax2,&
                           ixISmin3:ixISmax3,1:nws)
                      integer :: iB, iside

                      select case (ipole)
                      case (1)
                         iside=int((i1+3)/2)
                         iB=2*(1-1)+iside
                         select case (typeboundary(iw_mag(idirs),iB))
                         case (bc_symm)
                            wrecv(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,&
                                 idirs) = wsend(ixSmax1:ixSmin1:-1,ixSmin2:ixSmax2,&
                                 ixSmin3:ixSmax3,idirs)
                         case (bc_asymm)
                            wrecv(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,&
                                 idirs) =-wsend(ixSmax1:ixSmin1:-1,ixSmin2:ixSmax2,&
                                 ixSmin3:ixSmax3,idirs)
                         case default
                            call mpistop("Pole boundary condition should be symm or asymm")
                         end select

                      case (2)
                         iside=int((i2+3)/2)
                         iB=2*(2-1)+iside
                         select case (typeboundary(iw_mag(idirs),iB))
                         case (bc_symm)
                            wrecv(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,&
                                 idirs) = wsend(ixSmin1:ixSmax1,ixSmax2:ixSmin2:-1,&
                                 ixSmin3:ixSmax3,idirs)
                         case (bc_asymm)
                            wrecv(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,&
                                 idirs) =-wsend(ixSmin1:ixSmax1,ixSmax2:ixSmin2:-1,&
                                 ixSmin3:ixSmax3,idirs)
                         case default
                            call mpistop("Pole boundary condition should be symm or asymm")
                         end select

                      case (3)
                         iside=int((i3+3)/2)
                         iB=2*(3-1)+iside
                         select case (typeboundary(iw_mag(idirs),iB))
                         case (bc_symm)
                            wrecv(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,&
                                 idirs) = wsend(ixSmin1:ixSmax1,ixSmin2:ixSmax2,&
                                 ixSmax3:ixSmin3:-1,idirs)
                         case (bc_asymm)
                            wrecv(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,&
                                 idirs) =-wsend(ixSmin1:ixSmax1,ixSmin2:ixSmax2,&
                                 ixSmax3:ixSmin3:-1,idirs)
                         case default
                            call mpistop("Pole boundary condition should be symm or asymm")
                         end select

                      end select

                    end subroutine pole_copy_stg

                    subroutine pole_buffer(wrecv,ixIRmin1,ixIRmin2,ixIRmin3,ixIRmax1,&
                         ixIRmax2,ixIRmax3,ixRmin1,ixRmin2,ixRmin3,ixRmax1,ixRmax2,ixRmax3,&
                         wsend,ixISmin1,ixISmin2,ixISmin3,ixISmax1,ixISmax2,ixISmax3,ixSmin1,&
                         ixSmin2,ixSmin3,ixSmax1,ixSmax2,ixSmax3)

                      integer, intent(in) :: ixIRmin1,ixIRmin2,ixIRmin3,ixIRmax1,ixIRmax2,&
                           ixIRmax3,ixRmin1,ixRmin2,ixRmin3,ixRmax1,ixRmax2,ixRmax3,ixISmin1,&
                           ixISmin2,ixISmin3,ixISmax1,ixISmax2,ixISmax3,ixSmin1,ixSmin2,&
                           ixSmin3,ixSmax1,ixSmax2,ixSmax3
                      double precision :: wrecv(ixIRmin1:ixIRmax1,ixIRmin2:ixIRmax2,&
                           ixIRmin3:ixIRmax3,nwhead:nwtail), wsend(ixISmin1:ixISmax1,&
                           ixISmin2:ixISmax2,ixISmin3:ixISmax3,1:nw)

                      integer :: iw, iside, iB

                      select case (ipole)
                      case (1)
                         iside=int((i1+3)/2)
                         iB=2*(1-1)+iside
                         do iw=nwhead,nwtail
                            select case (typeboundary(iw,iB))
                            case (bc_symm)
                               wrecv(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,&
                                    iw) = wsend(ixSmax1:ixSmin1:-1,ixSmin2:ixSmax2,&
                                    ixSmin3:ixSmax3,iw)
                            case (bc_asymm)
                               wrecv(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,&
                                    iw) =-wsend(ixSmax1:ixSmin1:-1,ixSmin2:ixSmax2,&
                                    ixSmin3:ixSmax3,iw)
                            case default
                               call mpistop("Pole boundary condition should be symm or asymm")
                            end select
                         end do
                      case (2)
                         iside=int((i2+3)/2)
                         iB=2*(2-1)+iside
                         do iw=nwhead,nwtail
                            select case (typeboundary(iw,iB))
                            case (bc_symm)
                               wrecv(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,&
                                    iw) = wsend(ixSmin1:ixSmax1,ixSmax2:ixSmin2:-1,&
                                    ixSmin3:ixSmax3,iw)
                            case (bc_asymm)
                               wrecv(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,&
                                    iw) =-wsend(ixSmin1:ixSmax1,ixSmax2:ixSmin2:-1,&
                                    ixSmin3:ixSmax3,iw)
                            case default
                               call mpistop("Pole boundary condition should be symm or asymm")
                            end select
                         end do
                      case (3)
                         iside=int((i3+3)/2)
                         iB=2*(3-1)+iside
                         do iw=nwhead,nwtail
                            select case (typeboundary(iw,iB))
                            case (bc_symm)
                               wrecv(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,&
                                    iw) = wsend(ixSmin1:ixSmax1,ixSmin2:ixSmax2,&
                                    ixSmax3:ixSmin3:-1,iw)
                            case (bc_asymm)
                               wrecv(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,&
                                    iw) =-wsend(ixSmin1:ixSmax1,ixSmin2:ixSmax2,&
                                    ixSmax3:ixSmin3:-1,iw)
                            case default
                               call mpistop("Pole boundary condition should be symm or asymm")
                            end select
                         end do
                      end select

                    end subroutine pole_buffer


                  end subroutine getbc

                  subroutine bc_fill_srl(psb,igrid,nwhead,nwtail,i1,i2,i3,iib1,iib2,iib3)
                    !$acc routine vector
                    use mod_physicaldata, only: state
                    use mod_global_parameters, only: max_blocks, mype, ndim, npe
                    use mod_connectivity
                    integer, intent(in) :: igrid,i1,i2,i3,iib1,iib2,iib3,nwhead,nwtail
                    type(state), target :: psb(max_blocks)
                    integer :: ineighbor,ipe_neighbor,ipole,ixSmin1,ixSmin2,ixSmin3,&
                         ixSmax1,ixSmax2,ixSmax3,ixRmin1,ixRmin2,ixRmin3,ixRmax1,ixRmax2,&
                         ixRmax3,n_i1,n_i2,n_i3,idir,iw, ix1,ix2,ix3

                    ipe_neighbor=neighbor(2,i1,i2,i3,igrid)
                    if(ipe_neighbor==mype) then
                       ineighbor=neighbor(1,i1,i2,i3,igrid)
                       ipole=neighbor_pole(i1,i2,i3,igrid)
                       if(ipole==0) then
                          n_i1=-i1;n_i2=-i2;n_i3=-i3;
                          ixSmin1=ixS_srl_min1(iib1,i1);ixSmin2=ixS_srl_min2(iib2,i2)
                          ixSmin3=ixS_srl_min3(iib3,i3);ixSmax1=ixS_srl_max1(iib1,i1)
                          ixSmax2=ixS_srl_max2(iib2,i2);ixSmax3=ixS_srl_max3(iib3,i3);
                          ixRmin1=ixR_srl_min1(iib1,n_i1);ixRmin2=ixR_srl_min2(iib2,n_i2)
                          ixRmin3=ixR_srl_min3(iib3,n_i3);ixRmax1=ixR_srl_max1(iib1,n_i1)
                          ixRmax2=ixR_srl_max2(iib2,n_i2);ixRmax3=ixR_srl_max3(iib3,n_i3);

                          !$acc loop collapse(ndim+1) independent vector
                          do iw = nwhead, nwtail
                             do ix3=1,ixSmax3-ixSmin3+1
                                do ix2=1,ixSmax2-ixSmin2+1
                                   do ix1=1,ixSmax1-ixSmin1+1
                                      psb(ineighbor)%w(ixRmin1+ix1-1,ixRmin2+ix2-1,ixRmin3+ix3-1,&
                                           iw) = psb(igrid)%w(ixSmin1+ix1-1,ixSmin2+ix2-1,&
                                           ixSmin3+ix3-1,iw)
                                   end do
                                end do
                             end do
                          end do

                          !    if(stagger_grid) then
                          !       do idir=1,ndim
                          !          ixSmin1=ixS_srl_stg_min1(idir,i1)
                          !          ixSmin2=ixS_srl_stg_min2(idir,i2)
                          !          ixSmin3=ixS_srl_stg_min3(idir,i3)
                          !          ixSmax1=ixS_srl_stg_max1(idir,i1)
                          !          ixSmax2=ixS_srl_stg_max2(idir,i2)
                          !          ixSmax3=ixS_srl_stg_max3(idir,i3);
                          !          ixRmin1=ixR_srl_stg_min1(idir,n_i1)
                          !          ixRmin2=ixR_srl_stg_min2(idir,n_i2)
                          !          ixRmin3=ixR_srl_stg_min3(idir,n_i3)
                          !          ixRmax1=ixR_srl_stg_max1(idir,n_i1)
                          !          ixRmax2=ixR_srl_stg_max2(idir,n_i2)
                          !          ixRmax3=ixR_srl_stg_max3(idir,n_i3);
                          !          psb(ineighbor)%ws(ixRmin1:ixRmax1,ixRmin2:ixRmax2,&
                          !               ixRmin3:ixRmax3,idir)=psb(igrid)%ws(ixSmin1:ixSmax1,&
                          !               ixSmin2:ixSmax2,ixSmin3:ixSmax3,idir)
                          !       end do
                          !    end if
                          ! else
                          !    ixSmin1=ixS_srl_min1(iib1,i1);ixSmin2=ixS_srl_min2(iib2,i2)
                          !    ixSmin3=ixS_srl_min3(iib3,i3);ixSmax1=ixS_srl_max1(iib1,i1)
                          !    ixSmax2=ixS_srl_max2(iib2,i2);ixSmax3=ixS_srl_max3(iib3,i3);
                          !    select case (ipole)
                          !    case (1)
                          !       n_i1=i1;n_i2=-i2;n_i3=-i3;
                          !    case (2)
                          !       n_i1=-i1;n_i2=i2;n_i3=-i3;
                          !    case (3)
                          !       n_i1=-i1;n_i2=-i2;n_i3=i3;
                          !    end select
                          !    ixRmin1=ixR_srl_min1(iib1,n_i1);ixRmin2=ixR_srl_min2(iib2,n_i2)
                          !    ixRmin3=ixR_srl_min3(iib3,n_i3);ixRmax1=ixR_srl_max1(iib1,n_i1)
                          !    ixRmax2=ixR_srl_max2(iib2,n_i2);ixRmax3=ixR_srl_max3(iib3,n_i3);
                          !    call pole_copy(psb(ineighbor)%w,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,&
                          !         ixGhi3,ixRmin1,ixRmin2,ixRmin3,ixRmax1,ixRmax2,ixRmax3,&
                          !         psb(igrid)%w,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,ixSmin1,&
                          !         ixSmin2,ixSmin3,ixSmax1,ixSmax2,ixSmax3,ipole)
                          !    if(stagger_grid) then
                          !       do idir=1,ndim
                          !          ixSmin1=ixS_srl_stg_min1(idir,i1)
                          !          ixSmin2=ixS_srl_stg_min2(idir,i2)
                          !          ixSmin3=ixS_srl_stg_min3(idir,i3)
                          !          ixSmax1=ixS_srl_stg_max1(idir,i1)
                          !          ixSmax2=ixS_srl_stg_max2(idir,i2)
                          !          ixSmax3=ixS_srl_stg_max3(idir,i3);
                          !          ixRmin1=ixR_srl_stg_min1(idir,n_i1)
                          !          ixRmin2=ixR_srl_stg_min2(idir,n_i2)
                          !          ixRmin3=ixR_srl_stg_min3(idir,n_i3)
                          !          ixRmax1=ixR_srl_stg_max1(idir,n_i1)
                          !          ixRmax2=ixR_srl_stg_max2(idir,n_i2)
                          !          ixRmax3=ixR_srl_stg_max3(idir,n_i3);
                          !          call pole_copy_stg(psb(ineighbor)%ws,ixGslo1,ixGslo2,ixGslo3,&
                          !               ixGshi1,ixGshi2,ixGshi3,ixRmin1,ixRmin2,ixRmin3,ixRmax1,&
                          !               ixRmax2,ixRmax3,psb(igrid)%ws,ixGslo1,ixGslo2,ixGslo3,&
                          !               ixGshi1,ixGshi2,ixGshi3,ixSmin1,ixSmin2,ixSmin3,ixSmax1,&
                          !               ixSmax2,ixSmax3,idir,ipole)
                          !       end do
                          !    end if
                       end if
                    end if

                  end subroutine bc_fill_srl

                  subroutine identifyphysbound(s,iib1,iib2,iib3)
                    use mod_global_parameters

                    type(state)          :: s
                    integer, intent(out) :: iib1,iib2,iib3


                    if(s%is_physical_boundary(2*1) .and. s%is_physical_boundary(2*1-1)) then
                       iib1=2
                    else if(s%is_physical_boundary(2*1-1)) then
                       iib1=-1
                    else if(s%is_physical_boundary(2*1)) then
                       iib1=1
                    else
                       iib1=0
                    end if


                    if(s%is_physical_boundary(2*2) .and. s%is_physical_boundary(2*2-1)) then
                       iib2=2
                    else if(s%is_physical_boundary(2*2-1)) then
                       iib2=-1
                    else if(s%is_physical_boundary(2*2)) then
                       iib2=1
                    else
                       iib2=0
                    end if


                    if(s%is_physical_boundary(2*3) .and. s%is_physical_boundary(2*3-1)) then
                       iib3=2
                    else if(s%is_physical_boundary(2*3-1)) then
                       iib3=-1
                    else if(s%is_physical_boundary(2*3)) then
                       iib3=1
                    else
                       iib3=0
                    end if


                  end subroutine identifyphysbound

                end module mod_ghostcells_update
