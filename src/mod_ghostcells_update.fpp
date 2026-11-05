#:mute
#:include "mod_gpu_directives.fpp"
#:endmute
#ifdef _CRAYFTN
! locally switch off gpudirect since it does not work reliably with the nested datastructures on LUMI
#ifndef NOGPUDIRECT
#define NOGPUDIRECT 1
#endif
#endif

!> update ghost cells of all blocks including physical boundaries
module mod_ghostcells_update
#ifdef USE_MPIWRAPPERS
  use mod_mpi_wrapper
#else
#define mpi_irecv_wrapper MPI_IRECV
#define mpi_isend_wrapper MPI_ISEND
#endif
  implicit none
  ! A switch of update physical boundary or not
  logical, public :: bcphys=.true.
  ${GPU_DECLARE_COPYIN('bcphys')}$

  integer :: ixMmin1,ixMmin2,ixMmin3,ixMmax1,ixMmax2,ixMmax3, ixCoGmin1,&
       ixCoGmin2,ixCoGmin3,ixCoGmax1,ixCoGmax2,ixCoGmax3, ixCoMmin1,ixCoMmin2,&
       ixCoMmin3,ixCoMmax1,ixCoMmax2,ixCoMmax3, ixCoGsmin1,ixCoGsmin2,ixCoGsmin3,&
       ixCoGsmax1,ixCoGsmax2,ixCoGsmax3
  ${GPU_DECLARE_CREATE('ixCoGmin1,ixCoGmin2,ixCoGmin3,ixCoGmax1,ixCoGmax2,ixCoGmax3')}$
  ${GPU_DECLARE_CREATE('ixMmin1,ixMmin2,ixMmin3,ixMmax1,ixMmax2,ixMmax3')}$
  ${GPU_DECLARE_CREATE('ixCoMmin1,ixCoMmin2,ixCoMmin3,ixCoMmax1,ixCoMmax2,ixCoMmax3')}$
  ${GPU_DECLARE_CREATE('ixCoGsmin1,ixCoGsmin2,ixCoGsmin3,ixCoGsmax1,ixCoGsmax2,ixCoGsmax3')}$

  logical  :: req_diagonal = .true.
  ${GPU_DECLARE_COPYIN('req_diagonal')}$


  ! The first index goes from -1:2, where -1 is used when a block touches the
  ! lower boundary, 1 when a block touches an upper boundary, and 0 a situation
  ! away from boundary conditions, 2 when a block touched both lower and upper
  ! boundary

  ! index ranges to send (S) to sibling blocks, receive (R) from sibling blocks
  integer, dimension(-1:2,-1:1) :: ixS_srl_min1,ixS_srl_min2,ixS_srl_min3,&
       ixS_srl_max1,ixS_srl_max2,ixS_srl_max3, ixR_srl_min1,ixR_srl_min2,&
       ixR_srl_min3,ixR_srl_max1,ixR_srl_max2,ixR_srl_max3
  ${GPU_DECLARE_CREATE('ixS_srl_min1,ixS_srl_min2,ixS_srl_min3')}$
  ${GPU_DECLARE_CREATE('ixS_srl_max1,ixS_srl_max2,ixS_srl_max3, ixR_srl_min1,ixR_srl_min2')}$
  ${GPU_DECLARE_CREATE('ixR_srl_min3,ixR_srl_max1,ixR_srl_max2,ixR_srl_max3')}$

  ! index ranges of staggered variables to send (S) to sibling blocks, receive (R) from sibling blocks
  integer, dimension(3,-1:1) :: ixS_srl_stg_min1,ixS_srl_stg_min2,&
       ixS_srl_stg_min3,ixS_srl_stg_max1,ixS_srl_stg_max2,ixS_srl_stg_max3,&
       ixR_srl_stg_min1,ixR_srl_stg_min2,ixR_srl_stg_min3,ixR_srl_stg_max1,&
       ixR_srl_stg_max2,ixR_srl_stg_max3

  ! index ranges to send (S) restricted (r) ghost cells to coarser blocks
  integer, dimension(-1:1,-1:1) :: ixS_r_min1,ixS_r_min2,ixS_r_min3,ixS_r_max1,&
       ixS_r_max2,ixS_r_max3
  ${GPU_DECLARE_CREATE('ixS_r_min1,ixS_r_min2,ixS_r_min3,ixS_r_max1,ixS_r_max2,ixS_r_max3')}$

  ! index ranges of staggered variables to send (S) restricted (r) ghost cells to coarser blocks
  integer, dimension(3,-1:1) :: ixS_r_stg_min1,ixS_r_stg_min2,ixS_r_stg_min3,&
       ixS_r_stg_max1,ixS_r_stg_max2,ixS_r_stg_max3

  ! index ranges to receive restriced ghost cells from finer blocks
  integer, dimension(-1:1, 0:3) :: ixR_r_min1,ixR_r_min2,ixR_r_min3,ixR_r_max1,&
       ixR_r_max2,ixR_r_max3
  ${GPU_DECLARE_CREATE('ixR_r_min1,ixR_r_min2,ixR_r_min3,ixR_r_max1,ixR_r_max2,ixR_r_max3')}$

  ! index ranges of staggered variables to receive restriced ghost cells from finer blocks
  integer, dimension(3,0:3)  :: ixR_r_stg_min1,ixR_r_stg_min2,ixR_r_stg_min3,&
       ixR_r_stg_max1,ixR_r_stg_max2,ixR_r_stg_max3

  ! send prolongated (p) ghost cells to finer blocks, receive prolongated from coarser blocks
  integer, dimension(-1:1, 0:3) :: ixS_p_min1,ixS_p_min2,ixS_p_min3,ixS_p_max1,&
       ixS_p_max2,ixS_p_max3, ixR_p_min1,ixR_p_min2,ixR_p_min3,ixR_p_max1,&
       ixR_p_max2,ixR_p_max3
  ${GPU_DECLARE_CREATE('ixS_p_min1,ixS_p_min2,ixS_p_min3,ixS_p_max1,ixS_p_max2,ixS_p_max3')}$
  ${GPU_DECLARE_CREATE('ixR_p_min1,ixR_p_min2,ixR_p_min3,ixR_p_max1,ixR_p_max2,ixR_p_max3')}$

  !> One set of the prolongation ranges above, kept on the host. init_bc fills
  !> a full and a narrow set; select_prolong_ranges copies one of them into
  !> ixS_p_*/ixR_p_*, which are what the getbc kernels read.
  type prolong_ranges_t
     integer, dimension(-1:1, 0:3) :: S_min1, S_min2, S_min3, S_max1, S_max2, &
          S_max3, R_min1, R_min2, R_min3, R_max1, R_max2, R_max3
  end type prolong_ranges_t
  type(prolong_ranges_t) :: prolong_full, prolong_narrow
  ! which set ixS_p_*/ixR_p_* currently hold
  logical :: prolong_ranges_full = .true.

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

  ! total sizes = cell-center normal flux + stagger-grid flux of send and receive
  integer, dimension(-1:1,-1:1,-1:1) :: sizes_srl_send_total,&
       sizes_srl_recv_total

  ! size of srl buffers (non-staggered)
  integer, dimension(-1:1,-1:1,-1:1) :: sizes_srl_send, sizes_srl_recv

  ! sizes of buffer arrays for center-grid variable for siblings and restrict
  integer, dimension(:), allocatable :: recvrequest_c_sr, sendrequest_c_sr
  integer, dimension(:,:), allocatable :: recvstatus_c_sr, sendstatus_c_sr

  ! MPI requests and status for srl neighbor exchange, used with nprocs_info
  integer, dimension(:), allocatable :: recv_srl_nb, send_srl_nb, recv_f_nb, send_f_nb, send_c_nb, recv_c_nb
  integer, dimension(:,:), allocatable :: recvstatus_srl_nb, sendstatus_srl_nb
  integer, dimension(:,:), allocatable :: recvstatus_f_nb, sendstatus_f_nb
  integer, dimension(:,:), allocatable :: recvstatus_c_nb, sendstatus_c_nb

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


contains

#! Decode a neighbour code i (1..27, from iencode) into the direction offsets
#! i1, i2, i3 in {-1,0,1}. Written inline rather than as a device routine: a
#! call passed the kernel's private i1, i2, i3 by reference, and NVHPC could
#! hand the routine an invalid address for them, so the unpack kernel wrote
#! to a wild, misaligned location (compute-sanitizer: invalid __global__
#! write). Whether it happened depended on unrelated code generation.
#:def IDECODE(i1, i2, i3, i)
          ${i3}$ = ceiling(dble(${i}$)/9.0d0) - 2
          ${i2}$ = ceiling(dble(${i}$ - 9*(${i3}$+1))/3.0d0) - 2
          ${i1}$ = ${i}$ - 9*(${i3}$+1) - 3*(${i2}$+1) - 2
#:enddef

#! True when direction (i1,i2,i3) needs no ghost-cell exchange: the block
#! itself (0,0,0), and edge/corner directions when diagonal ghost cells are
#! not required. A logical expression rather than a device function, for the
#! same reason as IDECODE: no temporaries passed by reference into a call.
#:def SKIP_DIRECTION(i1, i2, i3)
((${i1}$ == 0 .and. ${i2}$ == 0 .and. ${i3}$ == 0) .or. (.not. req_diagonal .and. (merge(1,0,${i1}$ /= 0) + merge(1,0,${i2}$ /= 0) + merge(1,0,${i3}$ /= 0)) > 1))
#:enddef

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

    ixR_p_min1(:, 0)=ixCoMmin1-nghostcellsCo-(interpolation_order-1)
    ixR_p_min1(:, 1)=ixCoMmin1-(interpolation_order-1)
    ixR_p_min1(:, 2)=ixCoMmin1-nghostcellsCo-(interpolation_order-1)
    ixR_p_min1(:, 3)=ixCoMmax1+1-(interpolation_order-1)
    ixR_p_max1(:, 0)=nghostcells+(interpolation_order-1)
    ixR_p_max1(:, 1)=ixCoMmax1+nghostcellsCo+(interpolation_order-1)
    ixR_p_max1(:, 2)=ixCoMmax1+(interpolation_order-1)
    ixR_p_max1(:, 3)=ixCoMmax1+nghostcellsCo+(interpolation_order-1)



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

    ixR_p_min2(:, 0)=ixCoMmin2-nghostcellsCo-(interpolation_order-1)
    ixR_p_min2(:, 1)=ixCoMmin2-(interpolation_order-1)
    ixR_p_min2(:, 2)=ixCoMmin2-nghostcellsCo-(interpolation_order-1)
    ixR_p_min2(:, 3)=ixCoMmax2+1-(interpolation_order-1)
    ixR_p_max2(:, 0)=nghostcells+(interpolation_order-1)
    ixR_p_max2(:, 1)=ixCoMmax2+nghostcellsCo+(interpolation_order-1)
    ixR_p_max2(:, 2)=ixCoMmax2+(interpolation_order-1)
    ixR_p_max2(:, 3)=ixCoMmax2+nghostcellsCo+(interpolation_order-1)



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

    ixR_p_min3(:, 0)=ixCoMmin3-nghostcellsCo-(interpolation_order-1)
    ixR_p_min3(:, 1)=ixCoMmin3-(interpolation_order-1)
    ixR_p_min3(:, 2)=ixCoMmin3-nghostcellsCo-(interpolation_order-1)
    ixR_p_min3(:, 3)=ixCoMmax3+1-(interpolation_order-1)
    ixR_p_max3(:, 0)=nghostcells+(interpolation_order-1)
    ixR_p_max3(:, 1)=ixCoMmax3+nghostcellsCo+(interpolation_order-1)
    ixR_p_max3(:, 2)=ixCoMmax3+(interpolation_order-1)
    ixR_p_max3(:, 3)=ixCoMmax3+nghostcellsCo+(interpolation_order-1)

    if (stagger_grid) then
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

    ! The prolongation ranges come in two sets, and getbc picks one per call
    ! (select_prolong_ranges). The full set reaches into the ghost cells that
    ! the limited slopes of an edge or corner region read; the narrow set,
    ! used when diagonal ghost cells are skipped, keeps the sender's ghost
    ! cells out of the coarse buffer, because such a getbc leaves the
    ! sender's diagonal ghost cells stale. The physical-boundary adjustments
    ! above touch none of the entries the narrowing changes, so they hold
    ! for both sets. The MPI buffers in nbprocs_info are laid out for the
    ! full set, and a narrow call packs each chunk into the front of its slot.
    call store_prolong_ranges(prolong_full)
    if(.not.phys_req_diagonal) then
       ixS_p_min1(:, 0)=ixMmin1
       ixS_p_max1(:, 3)=ixMmax1
       ixS_p_max1(:, 1)=ixMmin1-1+nxCo1+(interpolation_order-1)
       ixS_p_min1(:, 2)=ixMmin1+nxCo1-(interpolation_order-1)
       ixR_p_max1(:, 0)=nghostcells
       ixR_p_min1(:, 3)=ixCoMmax1+1
       ixR_p_max1(:, 1)=ixCoMmax1+(interpolation_order-1)
       ixR_p_min1(:, 2)=ixCoMmin1-(interpolation_order-1)
       ixS_p_min2(:, 0)=ixMmin2
       ixS_p_max2(:, 3)=ixMmax2
       ixS_p_max2(:, 1)=ixMmin2-1+nxCo2+(interpolation_order-1)
       ixS_p_min2(:, 2)=ixMmin2+nxCo2-(interpolation_order-1)
       ixR_p_max2(:, 0)=nghostcells
       ixR_p_min2(:, 3)=ixCoMmax2+1
       ixR_p_max2(:, 1)=ixCoMmax2+(interpolation_order-1)
       ixR_p_min2(:, 2)=ixCoMmin2-(interpolation_order-1)
       ixS_p_min3(:, 0)=ixMmin3
       ixS_p_max3(:, 3)=ixMmax3
       ixS_p_max3(:, 1)=ixMmin3-1+nxCo3+(interpolation_order-1)
       ixS_p_min3(:, 2)=ixMmin3+nxCo3-(interpolation_order-1)
       ixR_p_max3(:, 0)=nghostcells
       ixR_p_min3(:, 3)=ixCoMmax3+1
       ixR_p_max3(:, 1)=ixCoMmax3+(interpolation_order-1)
       ixR_p_min3(:, 2)=ixCoMmin3-(interpolation_order-1)
    end if
    call store_prolong_ranges(prolong_narrow)
    call load_prolong_ranges(prolong_full)
    prolong_ranges_full = .true.

    ${GPU_UPDATE_DEVICE('ixS_srl_min1,ixS_srl_min2,ixS_srl_min3')}$
    ${GPU_UPDATE_DEVICE('ixS_srl_max1,ixS_srl_max2,ixS_srl_max3, ixR_srl_min1,ixR_srl_min2')}$
    ${GPU_UPDATE_DEVICE('ixR_srl_min3,ixR_srl_max1,ixR_srl_max2,ixR_srl_max3')}$
    ${GPU_UPDATE_DEVICE('ixCoGmin1,ixCoGmin2,ixCoGmin3,ixCoGmax1,ixCoGmax2,ixCoGmax3')}$
    ${GPU_UPDATE_DEVICE('ixCoMmin1,ixCoMmin2,ixCoMmin3,ixCoMmax1,ixCoMmax2,ixCoMmax3')}$
    ${GPU_UPDATE_DEVICE('ixCoGsmin1,ixCoGsmin2,ixCoGsmin3,ixCoGsmax1,ixCoGsmax2,ixCoGsmax3')}$
    ${GPU_UPDATE_DEVICE('ixMmin1,ixMmin2,ixMmin3,ixMmax1,ixMmax2,ixMmax3')}$
    ${GPU_UPDATE_DEVICE('ixS_r_min1,ixS_r_min2,ixS_r_min3,ixS_r_max1,ixS_r_max2,ixS_r_max3')}$
    ${GPU_UPDATE_DEVICE('ixR_r_min1,ixR_r_min2,ixR_r_min3,ixR_r_max1,ixR_r_max2,ixR_r_max3')}$
    ${GPU_UPDATE_DEVICE('ixS_p_min1,ixS_p_min2,ixS_p_min3,ixS_p_max1,ixS_p_max2,ixS_p_max3')}$
    ${GPU_UPDATE_DEVICE('ixR_p_min1,ixR_p_min2,ixR_p_min3,ixR_p_max1,ixR_p_max2,ixR_p_max3')}$

  end subroutine init_bc

  !> copy the current prolongation ranges into a set
  subroutine store_prolong_ranges(r)
    type(prolong_ranges_t), intent(out) :: r

    r%S_min1 = ixS_p_min1; r%S_min2 = ixS_p_min2; r%S_min3 = ixS_p_min3
    r%S_max1 = ixS_p_max1; r%S_max2 = ixS_p_max2; r%S_max3 = ixS_p_max3
    r%R_min1 = ixR_p_min1; r%R_min2 = ixR_p_min2; r%R_min3 = ixR_p_min3
    r%R_max1 = ixR_p_max1; r%R_max2 = ixR_p_max2; r%R_max3 = ixR_p_max3

  end subroutine store_prolong_ranges

  !> make a set the current prolongation ranges (host only)
  subroutine load_prolong_ranges(r)
    type(prolong_ranges_t), intent(in) :: r

    ixS_p_min1 = r%S_min1; ixS_p_min2 = r%S_min2; ixS_p_min3 = r%S_min3
    ixS_p_max1 = r%S_max1; ixS_p_max2 = r%S_max2; ixS_p_max3 = r%S_max3
    ixR_p_min1 = r%R_min1; ixR_p_min2 = r%R_min2; ixR_p_min3 = r%R_min3
    ixR_p_max1 = r%R_max1; ixR_p_max2 = r%R_max2; ixR_p_max3 = r%R_max3

  end subroutine load_prolong_ranges

  !> Give getbc the prolongation ranges for a fill with (full) or without
  !> diagonal ghost cells, on the host and the device. A no-op when that set
  !> is already current, which it is for every getbc of an advance in a row.
  subroutine select_prolong_ranges(full)
    logical, intent(in) :: full

    if (full .eqv. prolong_ranges_full) return

    if (full) then
       call load_prolong_ranges(prolong_full)
    else
       call load_prolong_ranges(prolong_narrow)
    end if
    prolong_ranges_full = full
    ${GPU_UPDATE_DEVICE('ixS_p_min1,ixS_p_min2,ixS_p_min3,ixS_p_max1,ixS_p_max2,ixS_p_max3')}$
    ${GPU_UPDATE_DEVICE('ixR_p_min1,ixR_p_min2,ixR_p_min3,ixR_p_max1,ixR_p_max2,ixR_p_max3')}$

  end subroutine select_prolong_ranges


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
    integer :: nwhead, nwtail, bgstep
    integer :: iigrid, igrid, ineighbor, ipe_neighbor
    integer :: ixRmin1,ixRmin2,ixRmin3,ixRmax1,ixRmax2,ixRmax3, ixSmin1,&
         ixSmin2,ixSmin3,ixSmax1,ixSmax2,ixSmax3
    integer :: i1,i2,i3, n_i1,n_i2,n_i3, ic1,ic2,ic3, inc1,inc2,inc3, n_inc1,&
         n_inc2, n_inc3, iib1, iib2, iib3
    real(dp) :: tempval

    integer :: ix1,ix2,ix3, iw, inb, i, Nx1, Nx2, Nx3, ienc, ibuf_start
    double precision :: CoFiratio
    ! for prolongation:
    integer  :: ixCo1, ixCo2, ixCo3, ixFi1, ixFi2, ixFi3, idims
    integer  :: ixFimin1, ixFimin2, ixFimin3, ixFimax1, ixFimax2, ixFimax3
    integer  :: hxCo1, hxCo2, hxCo3, jxCo1, jxCo2, jxCo3, nigrids
    real(dp) :: slope(ndim), dxFi1, dxFi2, dxFi3, dxCo1, dxCo2, dxCo3
    real(dp) :: invdxCo1, invdxCo2, invdxCo3
    real(dp) :: xFimin1, xFimin2, xFimin3, xComin1, xComin2, xComin3
    real(dp) :: xFi1, xFi2, xFi3, xCo3, xCo2, xCo1
    real(dp) :: eta1, eta2, eta3, slopeL, slopeR, slopeC, signR, signC
    ! Polar axis.  ipole is the dimension the pole lies across (0 for an
    ! ordinary neighbour, so every expression below collapses to the plain
    ! copy), iB_pole is the physical-boundary index whose typeboundary column
    ! holds the per-variable sign, and sbase/sstep walk the source backwards
    ! along the pole dimension so the row nearest the axis lands nearest the
    ! axis.  See "The polar axis" in CLAUDE.md for the derivation.
    integer  :: ipole, iB_pole
    integer  :: sbase1, sbase2, sbase3, sstep1, sstep2, sstep3

    time_bcin=MPI_WTIME()
    call nvtxStartRange("getbc",2)

    nwhead = nwstart
    nwtail = nwstart+nwbc-1
    ! igrids(1) is undefined when this rank owns no blocks
    if (igridstail > 0) then
       bgstep = psb(igrids(1))%istep
    else
       bgstep = 1
    end if

    req_diagonal = .true.
    if (present(req_diag)) req_diagonal = req_diag
    ${GPU_UPDATE_DEVICE('req_diagonal')}$
    ! the slopes of a diagonal prolongation need the full coarse buffer
    call select_prolong_ranges(req_diagonal)

    ! fill internal physical boundary
    if (internalboundary) then
       call getintbc(time,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3)
    end if

    ! fill physical-boundary ghost cells before internal ghost-cell values exchange
    if(bcphys.and. .not.stagger_grid) then
       ${GPU_PARALLEL_LOOP_GANG("private(igrid)")}$
       do iigrid = 1, igridstail; igrid=igrids(iigrid);
          if (.not.phyboundblock(igrid)) cycle
          call fill_boundary_before_gc(psb(igrid),igrid,time,qdt)
       end do
    end if


    ! prepare coarse values to send to coarser neighbors
    ${GPU_PARALLEL_LOOP_GANG("private(igrid, CoFiratio, i1, i2, i3)")}$
    do iigrid = 1, igridstail; igrid=igrids(iigrid);
       if (any(neighbor_type(:,:,:,igrid)==neighbor_coarse)) then

          CoFiratio=one/dble(2**ndim)
          ${GPU_LOOP_VECTOR("collapse(4) private(ixFi1, ixFi2, ixFi3, ix1, ix2, ix3)")}$
          do iw = nwhead, nwtail
             do ixCo3 = ixCoMmin3,ixCoMmax3
                do ixCo2 = ixCoMmin2,ixCoMmax2
                   do ixCo1 = ixCoMmin1,ixCoMmax1
                      ixFi3=2*(ixCo3-ixCoMmin3)+ixMmin3
                      ixFi2=2*(ixCo2-ixCoMmin2)+ixMmin2
                      ixFi1=2*(ixCo1-ixCoMmin1)+ixMmin1

                      bgc(1)%w(ixCo1,ixCo2,ixCo3,iw,igrid) = 0.0d0
                      do ix3 = ixFi3,ixFi3+1
                         do ix2 = ixFi2,ixFi2+1
                            do ix1 = ixFi1,ixFi1+1
                               bgc(1)%w(ixCo1,ixCo2,ixCo3,iw,igrid) = bgc(1)%w(ixCo1,ixCo2,ixCo3,iw,igrid) &
                                    + bg(bgstep)%w(ix1,ix2,ix3,iw,igrid)
                            end do
                         end do
                      end do
                      bgc(1)%w(ixCo1,ixCo2,ixCo3,iw,igrid) = bgc(1)%w(ixCo1,ixCo2,ixCo3,iw,igrid)*coFiRatio

                   end do
                end do
             end do
          end do

          do i3 = -1, 1
             do i2 = -1, 1
                do i1 = -1, 1
                   if (${SKIP_DIRECTION('i1','i2','i3')}$) cycle
                   if (neighbor_type(i1,i2,i3,igrid)==neighbor_coarse) call &
                        fill_coarse_boundary(time,igrid,i1,i2,i3)
                end do
             end do
          end do
       end if
    end do

    ! MPI receive SRL
    do inb = 1, nbprocs_info%nbprocs_srl
#ifndef NOGPUDIRECT
#ifdef _CRAYFTN
      ${GPU_HOST_DATA_USE_DEVICE('nbprocs_info')}$
#else
      ${GPU_HOST_DATA_USE_DEVICE('nbprocs_info%srl_nb(inb)%rcv%buffer, nbprocs_info%srl_nb(inb)%info_rcv%buffer')}$
#endif
#endif
       call mpi_irecv_wrapper(nbprocs_info%srl_nb(inb)%rcv%buffer, &
            nbprocs_info%srl_nb(inb)%rcv%size, &
            MPI_DOUBLE_PRECISION, nbprocs_info%nbprocs_srl_list(inb), 1, icomm, recv_srl_nb(inb), ierrmpi)
       call mpi_irecv_wrapper(nbprocs_info%srl_nb(inb)%info_rcv%buffer, &
            nbprocs_info%srl_nb(inb)%info_rcv%size, &
            MPI_INTEGER, nbprocs_info%nbprocs_srl_list(inb), 2, icomm, recv_srl_nb(nbprocs_info%nbprocs_srl + inb), ierrmpi)

#ifndef NOGPUDIRECT
      ${GPU_END_HOST_DATA()}$
#endif
    end do

    ! MPI receive restrict (neighbor is finer)
    do inb = 1, nbprocs_info%nbprocs_f
#ifndef NOGPUDIRECT
#ifdef _CRAYFTN
      ${GPU_HOST_DATA_USE_DEVICE('nbprocs_info')}$
#else
      ${GPU_HOST_DATA_USE_DEVICE('nbprocs_info%fine_nb(inb)%rcv%buffer, nbprocs_info%fine_nb(inb)%info_rcv%buffer')}$
#endif
#endif
       call mpi_irecv_wrapper(nbprocs_info%fine_nb(inb)%rcv%buffer, &
            nbprocs_info%fine_nb(inb)%rcv%size, &
            MPI_DOUBLE_PRECISION, nbprocs_info%nbprocs_f_list(inb), 1, icomm, recv_f_nb(inb), ierrmpi)
       call mpi_irecv_wrapper(nbprocs_info%fine_nb(inb)%info_rcv%buffer, &
            nbprocs_info%fine_nb(inb)%info_rcv%size, &
            MPI_INTEGER, nbprocs_info%nbprocs_f_list(inb), 2, icomm, recv_f_nb(nbprocs_info%nbprocs_f + inb), ierrmpi)
#ifndef NOGPUDIRECT
      ${GPU_END_HOST_DATA()}$
#endif
    end do

    ! fill the SRL send buffers on GPU
    do inb = 1, nbprocs_info%nbprocs_srl
       ${GPU_PARALLEL_LOOP_GANG("private(igrid, ienc, ibuf_start, i1, i2, i3, iib1, iib2, iib3, ixSmin1, ixSmin2, ixSmin3, ixSmax1, ixSmax2, ixSmax3, Nx1, Nx2, Nx3, ipole, iB_pole, n_i1, n_i2, n_i3, sbase1, sbase2, sbase3, sstep1, sstep2, sstep3)")}$
       do i = 1, nbprocs_info%srl_nb(inb)%info%nigrids
             igrid = nbprocs_info%srl_nb(inb)%info%igrid(i)
             ienc = nbprocs_info%srl_nb(inb)%info%iencode(i)
             ibuf_start = nbprocs_info%srl_nb(inb)%info%ibuf_start(i)
             @:IDECODE(i1, i2, i3, ienc)
             iib1=idphyb(1,igrid); iib2=idphyb(2,igrid); iib3=idphyb(3,igrid)

             ! now fill the data and info buffers
             ixSmin1=ixS_srl_min1(iib1,i1); ixSmax1=ixS_srl_max1(iib1,i1)
             ixSmin2=ixS_srl_min2(iib2,i2); ixSmax2=ixS_srl_max2(iib2,i2)
             ixSmin3=ixS_srl_min3(iib3,i3); ixSmax3=ixS_srl_max3(iib3,i3)
             Nx1=ixSmax1-ixSmin1+1; Nx2=ixSmax2-ixSmin2+1; Nx3=ixSmax3-ixSmin3+1

             ! The mirror and the sign are applied here, on the send side, so
             ! that the buffer already holds what the receiver has to store and
             ! the unpack stays a plain strided copy - which is what lets the
             ! receiver get away with knowing nothing about the pole beyond the
             ! destination offset shipped in the info record below.
             ipole=neighbor_pole(i1,i2,i3,igrid)
             n_i1=-i1; n_i2=-i2; n_i3=-i3
             sbase1=ixSmin1; sstep1=1; sbase2=ixSmin2; sstep2=1
             sbase3=ixSmin3; sstep3=1; iB_pole=1
             if (ipole==1) then
                n_i1=i1; sbase1=ixSmax1; sstep1=-1; iB_pole=(i1+3)/2
             else if (ipole==2) then
                n_i2=i2; sbase2=ixSmax2; sstep2=-1; iB_pole=2+(i2+3)/2
             else if (ipole==3) then
                n_i3=i3; sbase3=ixSmax3; sstep3=-1; iB_pole=4+(i3+3)/2
             end if

             ${GPU_LOOP_VECTOR("collapse(4)")}$ ${GPU_INDEPENDENT()}$
             do iw = nwhead, nwtail
                do ix3 = ixSmin3, ixSmax3
                   do ix2 = ixSmin2, ixSmax2
                      do ix1 = ixSmin1, ixSmax1
                         nbprocs_info%srl_nb(inb)%send%buffer( &
                              ibuf_start &
                              + sstep1*(ix1-sbase1) &
                              + Nx1 * (sstep2*(ix2-sbase2)) &
                              + Nx1*Nx2 * (sstep3*(ix3-sbase3)) &
                              + Nx1*Nx2*Nx3 * (iw-nwhead) &
                              ) = merge(-one,one,ipole/=0 .and. &
                              typeboundary(iw,iB_pole)==bc_asymm) &
                              * bg(bgstep)%w( ix1, ix2, ix3, iw, igrid)
                      end do
                   end do
                end do
             end do
             ! ship the destination offset rather than our own, so that the
             ! receiver does not have to undo the mirror it cannot see
             nbprocs_info%srl_nb(inb)%info_send%buffer( 1 + 3 * (i - 1) : 3 * i ) = &
                  [neighbor(1,i1,i2,i3,igrid), &
                   1 + (n_i1+1) + 3*(n_i2+1) + 9*(n_i3+1), ibuf_start]
       end do
    end do

    ! fill the C send buffers on GPU (send_restrict)
    do inb = 1, nbprocs_info%nbprocs_c
       ${GPU_PARALLEL_LOOP_GANG("private(Nx1,Nx2,Nx3,i1,i2,i3,inc1,inc2,inc3,igrid,ibuf_start,iib1,iib2,iib3,ixSmin1,ixSmin2,ixSmin3,ixSmax1,ixSmax2,ixSmax3,ipole,iB_pole,sbase1,sbase2,sbase3,sstep1,sstep2,sstep3)")}$
       do i = 1, nbprocs_info%course_nb(inb)%info%nigrids

          igrid = nbprocs_info%course_nb(inb)%info%igrid(i)
          i1 = nbprocs_info%course_nb(inb)%info%i1(i)
          i2 = nbprocs_info%course_nb(inb)%info%i2(i)
          i3 = nbprocs_info%course_nb(inb)%info%i3(i)
          ipole=neighbor_pole(i1,i2,i3,igrid)

          ! how they will be used at the receiving end.  The plain swap is
          ! n_inc = -2*i+ic; across a pole it becomes n_inc = 2*i+(3-ic),
          ! which on this 0..3 index leaves 0 and 3 alone and swaps 1 with 2.
          inc1 = nbprocs_info%course_nb(inb)%info%inc1(i)
          if (ipole == 1) then
             if (inc1 == 1) then; inc1 = 2; else if (inc1 == 2) then; inc1 = 1; end if
          else
             if (inc1 == 0) then; inc1 = 3; else if (inc1 == 3) then; inc1 = 0; end if
          end if

          inc2 = nbprocs_info%course_nb(inb)%info%inc2(i)
          if (ipole == 2) then
             if (inc2 == 1) then; inc2 = 2; else if (inc2 == 2) then; inc2 = 1; end if
          else
             if (inc2 == 0) then; inc2 = 3; else if (inc2 == 3) then; inc2 = 0; end if
          end if

          inc3 = nbprocs_info%course_nb(inb)%info%inc3(i)
          if (ipole == 3) then
             if (inc3 == 1) then; inc3 = 2; else if (inc3 == 2) then; inc3 = 1; end if
          else
             if (inc3 == 0) then; inc3 = 3; else if (inc3 == 3) then; inc3 = 0; end if
          end if

          ibuf_start = nbprocs_info%course_nb(inb)%info%ibuf_start(i)
          iib1=idphyb(1,igrid); iib2=idphyb(2,igrid); iib3=idphyb(3,igrid)

          ! now fill the data and info buffers
          ixSmin1=ixS_r_min1(iib1,i1);ixSmin2=ixS_r_min2(iib2,i2)
          ixSmin3=ixS_r_min3(iib3,i3);ixSmax1=ixS_r_max1(iib1,i1)
          ixSmax2=ixS_r_max2(iib2,i2);ixSmax3=ixS_r_max3(iib3,i3)
          Nx1=ixSmax1-ixSmin1+1; Nx2=ixSmax2-ixSmin2+1; Nx3=ixSmax3-ixSmin3+1

          ! mirror and sign on the send side, as for srl above
          sbase1=ixSmin1; sstep1=1; sbase2=ixSmin2; sstep2=1
          sbase3=ixSmin3; sstep3=1; iB_pole=1
          if (ipole==1) then
             sbase1=ixSmax1; sstep1=-1; iB_pole=(i1+3)/2
          else if (ipole==2) then
             sbase2=ixSmax2; sstep2=-1; iB_pole=2+(i2+3)/2
          else if (ipole==3) then
             sbase3=ixSmax3; sstep3=-1; iB_pole=4+(i3+3)/2
          end if

          ${GPU_LOOP_VECTOR("collapse(4)")}$ ${GPU_INDEPENDENT()}$
          do iw = nwhead, nwtail
             do ix3 = ixSmin3, ixSmax3
                do ix2 = ixSmin2, ixSmax2
                   do ix1 = ixSmin1, ixSmax1
                      nbprocs_info%course_nb(inb)%send%buffer( &
                           ibuf_start &
                           + sstep1*(ix1-sbase1) &
                           + Nx1 * (sstep2*(ix2-sbase2)) &
                           + Nx1*Nx2 * (sstep3*(ix3-sbase3)) &
                           + Nx1*Nx2*Nx3 * (iw-nwhead) &
                           ) = merge(-one,one,ipole/=0 .and. &
                           typeboundary(iw,iB_pole)==bc_asymm) &
                           * bgc(1)%w( ix1, ix2, ix3, iw, igrid )
                   end do
                end do
             end do
          end do

          nbprocs_info%course_nb(inb)%info_send%buffer( 1 + 5 * (i - 1) : 5 * i ) = &
               [neighbor(1,i1,i2,i3,igrid), inc1, inc2, inc3, ibuf_start]
       end do
    end do

#ifdef NOGPUDIRECT
    do inb = 1, nbprocs_info%nbprocs_srl
       ${GPU_UPDATE_HOST('nbprocs_info%srl_nb(inb)%info_send%buffer')}$
       ${GPU_UPDATE_HOST('nbprocs_info%srl_nb(inb)%send%buffer')}$
    end do
    do inb = 1, nbprocs_info%nbprocs_c
       ${GPU_UPDATE_HOST('nbprocs_info%course_nb(inb)%info_send%buffer')}$
       ${GPU_UPDATE_HOST('nbprocs_info%course_nb(inb)%send%buffer')}$
    end do
#endif

    ! MPI send SRL
    do inb = 1, nbprocs_info%nbprocs_srl
#ifndef NOGPUDIRECT
#ifdef _CRAYFTN
      ${GPU_HOST_DATA_USE_DEVICE('nbprocs_info')}$
#else
      ${GPU_HOST_DATA_USE_DEVICE('nbprocs_info%srl_nb(inb)%send%buffer, nbprocs_info%srl_nb(inb)%info_send%buffer')}$
#endif
#endif
       call mpi_isend_wrapper(nbprocs_info%srl_nb(inb)%send%buffer, &
            nbprocs_info%srl_nb(inb)%send%size, &
            MPI_DOUBLE_PRECISION, nbprocs_info%nbprocs_srl_list(inb), 1, icomm, send_srl_nb(inb), ierrmpi)
       call mpi_isend_wrapper(nbprocs_info%srl_nb(inb)%info_send%buffer, &
            nbprocs_info%srl_nb(inb)%info_send%size, &
            MPI_INTEGER, nbprocs_info%nbprocs_srl_list(inb), 2, icomm, send_srl_nb(nbprocs_info%nbprocs_srl + inb), ierrmpi)
#ifndef NOGPUDIRECT
      ${GPU_END_HOST_DATA()}$
#endif
    end do

    ! MPI send C (send_restrict)
    do inb = 1, nbprocs_info%nbprocs_c
#ifndef NOGPUDIRECT
#ifdef _CRAYFTN
      ${GPU_HOST_DATA_USE_DEVICE('nbprocs_info')}$
#else
      ${GPU_HOST_DATA_USE_DEVICE('nbprocs_info%course_nb(inb)%send%buffer, nbprocs_info%course_nb(inb)%info_send%buffer')}$
#endif
#endif
       call mpi_isend_wrapper(nbprocs_info%course_nb(inb)%send%buffer, &
            nbprocs_info%course_nb(inb)%send%size, &
            MPI_DOUBLE_PRECISION, nbprocs_info%nbprocs_c_list(inb), 1, icomm, send_c_nb(inb), ierrmpi)
       call mpi_isend_wrapper(nbprocs_info%course_nb(inb)%info_send%buffer, &
            nbprocs_info%course_nb(inb)%info_send%size, &
            MPI_INTEGER, nbprocs_info%nbprocs_c_list(inb), 2, icomm, send_c_nb(nbprocs_info%nbprocs_c + inb), ierrmpi)
#ifndef NOGPUDIRECT
      ${GPU_END_HOST_DATA()}$
#endif
    end do

    ! fill ghost-cell values of sibling blocks and if neighbor is coarser (f2c)
    ! same process case
    ${GPU_PARALLEL_LOOP_GANG("collapse(2) private(i1,i2,i3, igrid, iib1,iib2,iib3, ipe_neighbor, ineighbor, n_i1,n_i2,n_i3, ixSmin1,ixSmin2,ixSmin3,ixSmax1,ixSmax2,ixSmax3, ixRmin1,ixRmin2,ixRmin3,ixRmax1,ixRmax2,ixRmax3, ic1,ic2,ic3, n_inc1,n_inc2,n_inc3, ipole,iB_pole, sbase1,sbase2,sbase3, sstep1,sstep2,sstep3)")}$
    do iigrid = 1, igridstail
       do i = 1, 27
          @:IDECODE(i1, i2, i3, i)
          if (${SKIP_DIRECTION('i1','i2','i3')}$) cycle
          igrid=igrids(iigrid)
          iib1=idphyb(1,igrid); iib2=idphyb(2,igrid); iib3=idphyb(3,igrid)
          ipe_neighbor=neighbor(2,i1,i2,i3,igrid)

          if (ipe_neighbor==mype) then
             ineighbor=neighbor(1,i1,i2,i3,igrid)

             select case (neighbor_type(i1,i2,i3,igrid))
             case(neighbor_sibling)

                   ipole=neighbor_pole(i1,i2,i3,igrid)
                   ! across a pole the partner faces us from the same side, so
                   ! that component of the offset keeps its sign
                   n_i1=-i1; n_i2=-i2; n_i3=-i3
                   if (ipole==1) n_i1=i1
                   if (ipole==2) n_i2=i2
                   if (ipole==3) n_i3=i3

                   ixSmin1=ixS_srl_min1(iib1,i1);   ixSmin2=ixS_srl_min2(iib2,i2)
                   ixSmin3=ixS_srl_min3(iib3,i3);   ixSmax1=ixS_srl_max1(iib1,i1)
                   ixSmax2=ixS_srl_max2(iib2,i2);   ixSmax3=ixS_srl_max3(iib3,i3)
                   ixRmin1=ixR_srl_min1(iib1,n_i1); ixRmin2=ixR_srl_min2(iib2,n_i2)
                   ixRmin3=ixR_srl_min3(iib3,n_i3); ixRmax1=ixR_srl_max1(iib1,n_i1)
                   ixRmax2=ixR_srl_max2(iib2,n_i2); ixRmax3=ixR_srl_max3(iib3,n_i3)

                   sbase1=ixSmin1; sstep1=1; sbase2=ixSmin2; sstep2=1
                   sbase3=ixSmin3; sstep3=1; iB_pole=1
                   if (ipole==1) then
                      sbase1=ixSmax1; sstep1=-1; iB_pole=(i1+3)/2
                   else if (ipole==2) then
                      sbase2=ixSmax2; sstep2=-1; iB_pole=2+(i2+3)/2
                   else if (ipole==3) then
                      sbase3=ixSmax3; sstep3=-1; iB_pole=4+(i3+3)/2
                   end if

                   ${GPU_LOOP_VECTOR("collapse(ndim+1)")}$ ${GPU_INDEPENDENT()}$
                   do iw = nwhead, nwtail
                      do ix3=1,ixSmax3-ixSmin3+1
                         do ix2=1,ixSmax2-ixSmin2+1
                            do ix1=1,ixSmax1-ixSmin1+1
                               bg(bgstep)%w(ixRmin1+ix1-1,ixRmin2+ix2-1,ixRmin3+ix3-1,&
                                    iw,ineighbor) = merge(-one,one,ipole/=0 .and. &
                                    typeboundary(iw,iB_pole)==bc_asymm) &
                                    * bg(bgstep)%w(sbase1+sstep1*(ix1-1),sbase2+sstep2*(ix2-1),&
                                    sbase3+sstep3*(ix3-1),iw,igrid)
                            end do
                         end do
                      end do
                   end do

              case(neighbor_coarse)

                 ic1=1+modulo(node(pig1_,igrid)-1,2)
                 ic2=1+modulo(node(pig2_,igrid)-1,2)
                 ic3=1+modulo(node(pig3_,igrid)-1,2)

                 if(.not.(i1==0.or.i1==2*ic1-3).or..not.(i2==0.or.i2==2*ic2-3)&
                      .or..not.(i3==0.or.i3==2*ic3-3)) cycle

                   ipole=neighbor_pole(i1,i2,i3,igrid)
                   ! same rule as srl, expressed in the child-offset index:
                   ! the pole component keeps 2*i and mirrors the child index
                   n_inc1=-2*i1+ic1;n_inc2=-2*i2+ic2;n_inc3=-2*i3+ic3;
                   if (ipole==1) n_inc1=2*i1+(3-ic1)
                   if (ipole==2) n_inc2=2*i2+(3-ic2)
                   if (ipole==3) n_inc3=2*i3+(3-ic3)

                   ixSmin1=ixS_r_min1(iib1,i1);ixSmin2=ixS_r_min2(iib2,i2)
                   ixSmin3=ixS_r_min3(iib3,i3);ixSmax1=ixS_r_max1(iib1,i1)
                   ixSmax2=ixS_r_max2(iib2,i2);ixSmax3=ixS_r_max3(iib3,i3);
                   ixRmin1=ixR_r_min1(iib1,n_inc1);ixRmin2=ixR_r_min2(iib2,n_inc2)
                   ixRmin3=ixR_r_min3(iib3,n_inc3);ixRmax1=ixR_r_max1(iib1,n_inc1)
                   ixRmax2=ixR_r_max2(iib2,n_inc2);ixRmax3=ixR_r_max3(iib3,n_inc3);

                   sbase1=ixSmin1; sstep1=1; sbase2=ixSmin2; sstep2=1
                   sbase3=ixSmin3; sstep3=1; iB_pole=1
                   if (ipole==1) then
                      sbase1=ixSmax1; sstep1=-1; iB_pole=(i1+3)/2
                   else if (ipole==2) then
                      sbase2=ixSmax2; sstep2=-1; iB_pole=2+(i2+3)/2
                   else if (ipole==3) then
                      sbase3=ixSmax3; sstep3=-1; iB_pole=4+(i3+3)/2
                   end if

                   ${GPU_LOOP_VECTOR("collapse(ndim+1)")}$ ${GPU_INDEPENDENT()}$
                   do iw = nwhead, nwtail
                      do ix3=1,ixSmax3-ixSmin3+1
                         do ix2=1,ixSmax2-ixSmin2+1
                            do ix1=1,ixSmax1-ixSmin1+1
                               bg(bgstep)%w(ixRmin1+ix1-1,ixRmin2+ix2-1,&
                                    ixRmin3+ix3-1,iw,ineighbor) = merge(-one,one,&
                                    ipole/=0 .and. typeboundary(iw,iB_pole)==bc_asymm) &
                                    * bgc(1)%w(sbase1+sstep1*(ix1-1),&
                                    sbase2+sstep2*(ix2-1),sbase3+sstep3*(ix3-1),iw, igrid)
                            end do
                         end do
                      end do
                   end do

                end select

            end if

         end do
    end do


    call MPI_WAITALL(nbprocs_info%nbprocs_srl*2, recv_srl_nb, recvstatus_srl_nb, ierrmpi)
    call MPI_WAITALL(nbprocs_info%nbprocs_srl*2, send_srl_nb, sendstatus_srl_nb, ierrmpi)
#ifdef NOGPUDIRECT
    do inb = 1, nbprocs_info%nbprocs_srl
      ${GPU_UPDATE_DEVICE('nbprocs_info%srl_nb(inb)%info_rcv%buffer')}$
      ${GPU_UPDATE_DEVICE('nbprocs_info%srl_nb(inb)%rcv%buffer')}$
    end do
#endif

    ! unpack the MPI buffers
    do inb = 1, nbprocs_info%nbprocs_srl
      ${GPU_PARALLEL_LOOP_GANG("private(igrid, ienc, ibuf_start, i1, i2, i3, iib1, iib2, iib3, ixRmin1, ixRmin2, ixRmin3, ixRmax1, ixRmax2, ixRmax3, Nx1, Nx2, Nx3)")}$ ${GPU_INDEPENDENT()}$
       do i = 1, nbprocs_info%srl_nb(inb)%info%nigrids

          igrid       = nbprocs_info%srl_nb(inb)%info_rcv%buffer( 3 * (i - 1) + 1 )
          ienc        = nbprocs_info%srl_nb(inb)%info_rcv%buffer( 3 * (i - 1) + 2 )
          ibuf_start  = nbprocs_info%srl_nb(inb)%info_rcv%buffer( 3 * (i - 1) + 3 )

          ! ienc already encodes our own ghost-region offset: the sender
          ! computed it, mirrored across a pole where there is one
          @:IDECODE(i1, i2, i3, ienc)

          iib1 = idphyb(1,igrid); iib2 = idphyb(2,igrid); iib3 = idphyb(3,igrid)

          ixRmin1=ixR_srl_min1(iib1,i1); ixRmin2=ixR_srl_min2(iib2,i2)
          ixRmin3=ixR_srl_min3(iib3,i3); ixRmax1=ixR_srl_max1(iib1,i1)
          ixRmax2=ixR_srl_max2(iib2,i2); ixRmax3=ixR_srl_max3(iib3,i3)
          Nx1=ixRmax1-ixRmin1+1; Nx2=ixRmax2-ixRmin2+1; Nx3=ixRmax3-ixRmin3+1

          ${GPU_LOOP_VECTOR("collapse(4) private(tempval)")}$ ${GPU_INDEPENDENT()}$
          do iw = nwhead, nwtail
             do ix3 = ixRmin3, ixRmax3
                do ix2 = ixRmin2, ixRmax2
                   do ix1 = ixRmin1, ixRmax1
                      ! going through a tempval is a workaround for Cray, which gives
                      ! a memory access fault on the GPUs otherwise
                      tempval = nbprocs_info%srl_nb(inb)%rcv%buffer( &
                           ibuf_start &
                           + (ix1-ixRmin1) &
                           + Nx1 * (ix2-ixRmin2) &
                           + Nx1*Nx2 * (ix3-ixRmin3) &
                           + Nx1*Nx2*Nx3 * (iw-nwhead) &
                           )
                      bg(bgstep)%w( ix1, ix2, ix3, iw , igrid) = tempval
                   end do
                end do
             end do
          end do

       end do
    end do

    call MPI_WAITALL(nbprocs_info%nbprocs_f*2, recv_f_nb, recvstatus_f_nb, ierrmpi)
    call MPI_WAITALL(nbprocs_info%nbprocs_c*2, send_c_nb, sendstatus_c_nb, ierrmpi)

#ifdef NOGPUDIRECT
    do inb = 1, nbprocs_info%nbprocs_f
       ${GPU_UPDATE_DEVICE('nbprocs_info%fine_nb(inb)%info_rcv%buffer(1:nbprocs_info%fine_nb(inb)%info_rcv%size)')}$
       ${GPU_UPDATE_DEVICE('nbprocs_info%fine_nb(inb)%rcv%buffer(1:nbprocs_info%fine_nb(inb)%rcv%size)')}$
    end do
#endif

    ! unpack the MPI buffers, fine neighbor, (f_recv), recv_restrict
    do inb = 1, nbprocs_info%nbprocs_f
       ${GPU_PARALLEL_LOOP_GANG("private(igrid, inc1, inc2, inc3, ibuf_start, iib1, iib2, iib3, ixRmin1, ixRmin2, ixRmin3, ixRmax1, ixRmax2, ixRmax3, Nx1, Nx2, Nx3)")}$
       do i = 1,nbprocs_info%fine_nb(inb)%info%nigrids

          igrid       = nbprocs_info%fine_nb(inb)%info_rcv%buffer( 5 * (i - 1) + 1 )
          inc1        = nbprocs_info%fine_nb(inb)%info_rcv%buffer( 5 * (i - 1) + 2 )
          inc2        = nbprocs_info%fine_nb(inb)%info_rcv%buffer( 5 * (i - 1) + 3 )
          inc3        = nbprocs_info%fine_nb(inb)%info_rcv%buffer( 5 * (i - 1) + 4 )
          ibuf_start  = nbprocs_info%fine_nb(inb)%info_rcv%buffer( 5 * (i - 1) + 5 )

          iib1 = idphyb(1,igrid); iib2 = idphyb(2,igrid); iib3 = idphyb(3,igrid)

          ixRmin1=ixR_r_min1(iib1,inc1); ixRmin2=ixR_r_min2(iib2,inc2)
          ixRmin3=ixR_r_min3(iib3,inc3); ixRmax1=ixR_r_max1(iib1,inc1)
          ixRmax2=ixR_r_max2(iib2,inc2); ixRmax3=ixR_r_max3(iib3,inc3)
          Nx1=ixRmax1-ixRmin1+1; Nx2=ixRmax2-ixRmin2+1; Nx3=ixRmax3-ixRmin3+1

          ${GPU_LOOP_VECTOR("collapse(4) private(tempval)")}$ ${GPU_INDEPENDENT()}$
          do iw = nwhead, nwtail
             do ix3 = ixRmin3, ixRmax3
                do ix2 = ixRmin2, ixRmax2
                   do ix1 = ixRmin1, ixRmax1
                      ! going through a tempval is a workaround for Cray, which gives
                      ! a memory access fault on the GPUs otherwise
                      tempval = nbprocs_info%fine_nb(inb)%rcv%buffer( &
                           ibuf_start &
                           + (ix1-ixRmin1) &
                           + Nx1 * (ix2-ixRmin2) &
                           + Nx1*Nx2 * (ix3-ixRmin3) &
                           + Nx1*Nx2*Nx3 * (iw-nwhead) &
                           )
                      bg(bgstep)%w( ix1, ix2, ix3, iw, igrid ) = tempval
                   end do
                end do
             end do
          end do

       end do
    end do

    ! MPI receive prolong (neighbor is coarser)
    do inb = 1, nbprocs_info%nbprocs_c
#ifndef NOGPUDIRECT
#ifdef _CRAYFTN
      ${GPU_HOST_DATA_USE_DEVICE('nbprocs_info')}$
#else
      ${GPU_HOST_DATA_USE_DEVICE('nbprocs_info%course_nb(inb)%rcv%buffer, nbprocs_info%course_nb(inb)%info_rcv%buffer')}$
#endif
#endif
       call mpi_irecv_wrapper(nbprocs_info%course_nb(inb)%rcv%buffer, &
            nbprocs_info%course_nb(inb)%rcv%size, &
            MPI_DOUBLE_PRECISION, nbprocs_info%nbprocs_c_list(inb), 1, icomm, recv_c_nb(inb), ierrmpi)
       call mpi_irecv_wrapper(nbprocs_info%course_nb(inb)%info_rcv%buffer, &
            nbprocs_info%course_nb(inb)%info_rcv%size, &
            MPI_INTEGER, nbprocs_info%nbprocs_c_list(inb), 2, icomm, recv_c_nb(nbprocs_info%nbprocs_c + inb), ierrmpi)
#ifndef NOGPUDIRECT
      ${GPU_END_HOST_DATA()}$
#endif
    end do

    ! fill the F (neighbor is finer) send buffer on GPU (send_prolong)
    do inb = 1, nbprocs_info%nbprocs_f
       ${GPU_PARALLEL_LOOP_GANG("private(Nx1,Nx2,Nx3,inc1,inc2,inc3,n_inc1,n_inc2,n_inc3,i1,i2,i3,igrid,ibuf_start,iib1,iib2,iib3,ixSmin1,ixSmin2,ixSmin3,ixSmax1,ixSmax2,ixSmax3,ipole,iB_pole,sbase1,sbase2,sbase3,sstep1,sstep2,sstep3)")}$ ${GPU_INDEPENDENT()}$
       do i = 1,nbprocs_info%fine_nb(inb)%info%nigrids

          igrid = nbprocs_info%fine_nb(inb)%info%igrid(i)
          inc1 = nbprocs_info%fine_nb(inb)%info%inc1(i)
          inc2 = nbprocs_info%fine_nb(inb)%info%inc2(i)
          inc3 = nbprocs_info%fine_nb(inb)%info%inc3(i)

          ! only inc = 2*i+ic is stored, and neighbor_pole is indexed by i;
          ! ic is 1 or 2, so inc pins i down (0 -> -1, 1 or 2 -> 0, 3 -> +1)
          if (inc1 == 0) then; i1 = -1; else if (inc1 == 3) then; i1 = 1; else; i1 = 0; end if
          if (inc2 == 0) then; i2 = -1; else if (inc2 == 3) then; i2 = 1; else; i2 = 0; end if
          if (inc3 == 0) then; i3 = -1; else if (inc3 == 3) then; i3 = 1; else; i3 = 0; end if
          ipole=neighbor_pole(i1,i2,i3,igrid)

          ibuf_start = nbprocs_info%fine_nb(inb)%info%ibuf_start(i)
          iib1=idphyb(1,igrid); iib2=idphyb(2,igrid); iib3=idphyb(3,igrid)

          ! now fill the data and info buffers
          ixSmin1=ixS_p_min1(iib1,inc1);ixSmin2=ixS_p_min2(iib2,inc2)
          ixSmin3=ixS_p_min3(iib3,inc3);ixSmax1=ixS_p_max1(iib1,inc1)
          ixSmax2=ixS_p_max2(iib2,inc2);ixSmax3=ixS_p_max3(iib3,inc3)
          Nx1=ixSmax1-ixSmin1+1; Nx2=ixSmax2-ixSmin2+1; Nx3=ixSmax3-ixSmin3+1

          ! mirror and sign on the send side, as for srl and restrict above
          sbase1=ixSmin1; sstep1=1; sbase2=ixSmin2; sstep2=1
          sbase3=ixSmin3; sstep3=1; iB_pole=1
          if (ipole==1) then
             sbase1=ixSmax1; sstep1=-1; iB_pole=(i1+3)/2
          else if (ipole==2) then
             sbase2=ixSmax2; sstep2=-1; iB_pole=2+(i2+3)/2
          else if (ipole==3) then
             sbase3=ixSmax3; sstep3=-1; iB_pole=4+(i3+3)/2
          end if

          ${GPU_LOOP_VECTOR("collapse(4)")}$ ${GPU_INDEPENDENT()}$
          do iw = nwhead, nwtail
             do ix3 = ixSmin3, ixSmax3
                do ix2 = ixSmin2, ixSmax2
                   do ix1 = ixSmin1, ixSmax1
                      nbprocs_info%fine_nb(inb)%send%buffer( &
                           ibuf_start &
                           + sstep1*(ix1-sbase1) &
                           + Nx1 * (sstep2*(ix2-sbase2)) &
                           + Nx1*Nx2 * (sstep3*(ix3-sbase3)) &
                           + Nx1*Nx2*Nx3 * (iw-nwhead) &
                           ) = merge(-one,one,ipole/=0 .and. &
                           typeboundary(iw,iB_pole)==bc_asymm) &
                           * bg(bgstep)%w( ix1, ix2, ix3, iw, igrid )
                   end do
                end do
             end do
          end do

          ! how they will be used at the receiving end.  The plain rule is
          ! n_inc = ic-i, the 0<->3 swap; across a pole n_inc = inc, no swap.
          n_inc1 = inc1; n_inc2 = inc2; n_inc3 = inc3
          if (ipole /= 1) then
             if (n_inc1 == 0) then; n_inc1 = 3; else if (n_inc1 == 3) then; n_inc1 = 0; end if
          end if
          if (ipole /= 2) then
             if (n_inc2 == 0) then; n_inc2 = 3; else if (n_inc2 == 3) then; n_inc2 = 0; end if
          end if
          if (ipole /= 3) then
             if (n_inc3 == 0) then; n_inc3 = 3; else if (n_inc3 == 3) then; n_inc3 = 0; end if
          end if

          nbprocs_info%fine_nb(inb)%info_send%buffer( 1 + 5 * (i - 1) : 5 * i ) = &
               [neighbor_child(1,inc1,inc2,inc3,igrid), n_inc1, n_inc2, n_inc3, ibuf_start]
       end do
    end do

#ifdef NOGPUDIRECT
    do inb = 1, nbprocs_info%nbprocs_f
      ${GPU_UPDATE_HOST('nbprocs_info%fine_nb(inb)%info_send%buffer')}$
      ${GPU_UPDATE_HOST('nbprocs_info%fine_nb(inb)%send%buffer')}$
    end do
#endif

    ! MPI send F (send_prolong)
    do inb = 1, nbprocs_info%nbprocs_f
#ifndef NOGPUDIRECT
#ifdef _CRAYFTN
      ${GPU_HOST_DATA_USE_DEVICE('nbprocs_info')}$
#else
      ${GPU_HOST_DATA_USE_DEVICE('nbprocs_info%fine_nb(inb)%send%buffer, nbprocs_info%fine_nb(inb)%info_send%buffer')}$
#endif
#endif
       call mpi_isend_wrapper(nbprocs_info%fine_nb(inb)%send%buffer, &
            nbprocs_info%fine_nb(inb)%send%size, &
            MPI_DOUBLE_PRECISION, nbprocs_info%nbprocs_f_list(inb), 1, icomm, send_f_nb(inb), ierrmpi)
       call mpi_isend_wrapper(nbprocs_info%fine_nb(inb)%info_send%buffer, &
            nbprocs_info%fine_nb(inb)%info_send%size, &
            MPI_INTEGER, nbprocs_info%nbprocs_f_list(inb), 2, icomm, send_f_nb(nbprocs_info%nbprocs_f + inb), ierrmpi)
#ifndef NOGPUDIRECT
      ${GPU_END_HOST_DATA()}$
#endif
    end do


    ! fill coarse ghost-cell values of finer neighbors in the same processor
    ${GPU_PARALLEL_LOOP_GANG("collapse(4) private(iib1,iib2,iib3,igrid,ipole,iB_pole,sbase1,sbase2,sbase3,sstep1,sstep2,sstep3, inc1,inc2,inc3, ineighbor,ipe_neighbor, n_i1,n_i2,n_i3, n_inc1,n_inc2,n_inc3, ixSmin1,ixSmin2,ixSmin3,ixSmax1,ixSmax2,ixSmax3, ixRmin1,ixRmin2,ixRmin3,ixRmax1,ixRmax2,ixRmax3)")}$
    do iigrid=1,igridstail
       do i3=-1,1
          do i2=-1,1
             do i1=-1,1
                if (${SKIP_DIRECTION('i1','i2','i3')}$) cycle
                igrid = igrids(iigrid)
                if (neighbor_type(i1,i2,i3,igrid)==neighbor_fine) then
                   iib1=idphyb(1,igrid);iib2=idphyb(2,igrid);iib3=idphyb(3,igrid);
                   ipole=neighbor_pole(i1,i2,i3,igrid)
                   iB_pole=1
                   if (ipole==1) then
                      iB_pole=(i1+3)/2
                   else if (ipole==2) then
                      iB_pole=2+(i2+3)/2
                   else if (ipole==3) then
                      iB_pole=4+(i3+3)/2
                   end if
                   !  inline of call bc_fill_prolong(igrid,i1,i2,i3,iib1,iib2,iib3) :

                   do ic3 = 1+int((1-i3)/2), 2-int((1+i3)/2)
                      inc3 = 2*i3+ic3
                      do ic2 = 1+int((1-i2)/2), 2-int((1+i2)/2)
                         inc2 = 2*i2+ic2
                         do ic1 = 1+int((1-i1)/2), 2-int((1+i1)/2)
                            inc1 = 2*i1+ic1
                            ipe_neighbor = neighbor_child(2,inc1,inc2,inc3,igrid)
                            if(ipe_neighbor==mype) then
                               ixSmin1=ixS_p_min1(iib1,inc1);ixSmin2=ixS_p_min2(iib2,inc2)
                               ixSmin3=ixS_p_min3(iib3,inc3);ixSmax1=ixS_p_max1(iib1,inc1)
                               ixSmax2=ixS_p_max2(iib2,inc2);ixSmax3=ixS_p_max3(iib3,inc3)

                               ineighbor = neighbor_child(1,inc1,inc2,inc3,igrid)
                               ! plain rule n_inc = ic-i; across a pole the
                               ! child sits at the same offset, n_inc = inc
                               n_i1=-i1; n_i2=-i2; n_i3=-i3
                               n_inc1=ic1+n_i1; n_inc2=ic2+n_i2; n_inc3=ic3+n_i3
                               if (ipole==1) n_inc1=inc1
                               if (ipole==2) n_inc2=inc2
                               if (ipole==3) n_inc3=inc3

                               ixRmin1=ixR_p_min1(iib1,n_inc1)
                               ixRmin2=ixR_p_min2(iib2,n_inc2)
                               ixRmin3=ixR_p_min3(iib3,n_inc3)
                               ixRmax1=ixR_p_max1(iib1,n_inc1)
                               ixRmax2=ixR_p_max2(iib2,n_inc2)
                               ixRmax3=ixR_p_max3(iib3,n_inc3)

                               sbase1=ixSmin1; sstep1=1; sbase2=ixSmin2; sstep2=1
                               sbase3=ixSmin3; sstep3=1
                               if (ipole==1) then
                                  sbase1=ixSmax1; sstep1=-1
                               else if (ipole==2) then
                                  sbase2=ixSmax2; sstep2=-1
                               else if (ipole==3) then
                                  sbase3=ixSmax3; sstep3=-1
                               end if

                               ${GPU_LOOP_VECTOR("collapse(4)")}$ ${GPU_INDEPENDENT()}$
                               do iw = nwhead, nwtail
                                  do ix3 =0, ixRmax3-ixRmin3
                                     do ix2 = 0, ixRmax2-ixRmin2
                                        do ix1 = 0, ixRmax1-ixRmin1
                                           bgc(1)%w(ixRmin1+ix1,ixRmin2+ix2,&
                                                ixRmin3+ix3,iw,ineighbor) = merge(-one,one,&
                                                ipole/=0 .and. typeboundary(iw,iB_pole)==bc_asymm) &
                                                * bg(bgstep)%w(sbase1+sstep1*ix1,&
                                                sbase2+sstep2*ix2,sbase3+sstep3*ix3,iw, igrid)
                                        end do
                                     end do
                                  end do
                               end do

                            end if

                         end do
                      end do
                   end do

                end if

             end do
          end do
       end do

    end do

    call MPI_WAITALL(nbprocs_info%nbprocs_c*2, recv_c_nb, recvstatus_c_nb, ierrmpi)
    call MPI_WAITALL(nbprocs_info%nbprocs_f*2, send_f_nb, sendstatus_f_nb, ierrmpi)


#ifdef NOGPUDIRECT
    do inb = 1, nbprocs_info%nbprocs_c
      ${GPU_UPDATE_DEVICE('nbprocs_info%course_nb(inb)%info_rcv%buffer(1:nbprocs_info%course_nb(inb)%info_rcv%size)')}$
      ${GPU_UPDATE_DEVICE('nbprocs_info%course_nb(inb)%rcv%buffer(1:nbprocs_info%course_nb(inb)%rcv%size)')}$
    end do
#endif

    ! unpack the MPI buffers, coarse neighbor, (c_recv), recv_prolong
    do inb = 1, nbprocs_info%nbprocs_c
       ${GPU_PARALLEL_LOOP_GANG("private(Nx1,Nx2,Nx3,inc1,inc2,inc3,igrid,ibuf_start,iib1,iib2,iib3,ixRmin1,ixRmin2,ixRmin3,ixRmax1,ixRmax2,ixRmax3)")}$ ${GPU_INDEPENDENT()}$
       do i = 1, nbprocs_info%course_nb(inb)%info%nigrids

          igrid       = nbprocs_info%course_nb(inb)%info_rcv%buffer( 5 * (i - 1) + 1 )
          inc1        = nbprocs_info%course_nb(inb)%info_rcv%buffer( 5 * (i - 1) + 2 )
          inc2        = nbprocs_info%course_nb(inb)%info_rcv%buffer( 5 * (i - 1) + 3 )
          inc3        = nbprocs_info%course_nb(inb)%info_rcv%buffer( 5 * (i - 1) + 4 )
          ibuf_start  = nbprocs_info%course_nb(inb)%info_rcv%buffer( 5 * (i - 1) + 5 )

          iib1 = idphyb(1,igrid); iib2 = idphyb(2,igrid); iib3 = idphyb(3,igrid)

          ixRmin1=ixR_p_min1(iib1,inc1); ixRmin2=ixR_p_min2(iib2,inc2)
          ixRmin3=ixR_p_min3(iib3,inc3); ixRmax1=ixR_p_max1(iib1,inc1)
          ixRmax2=ixR_p_max2(iib2,inc2); ixRmax3=ixR_p_max3(iib3,inc3)
          Nx1=ixRmax1-ixRmin1+1; Nx2=ixRmax2-ixRmin2+1; Nx3=ixRmax3-ixRmin3+1

          ${GPU_LOOP_VECTOR("collapse(4) private(tempval)")}$ ${GPU_INDEPENDENT()}$
          do iw = nwhead, nwtail
             do ix3 = ixRmin3, ixRmax3
                do ix2 = ixRmin2, ixRmax2
                   do ix1 = ixRmin1, ixRmax1
                      ! going through a tempval is a workaround for Cray, which gives
                      ! a memory access fault on the GPUs otherwise
                      tempval = nbprocs_info%course_nb(inb)%rcv%buffer( &
                           ibuf_start &
                           + (ix1-ixRmin1) &
                           + Nx1 * (ix2-ixRmin2) &
                           + Nx1*Nx2 * (ix3-ixRmin3) &
                           + Nx1*Nx2*Nx3 * (iw-nwhead) &
                           )
                      bgc(1)%w( ix1, ix2, ix3, iw, igrid ) = tempval
                   end do
                end do
             end do
          end do

       end do
    end do

    ! do prolongation on the ghost-cell values based on the received coarse values from coarser neighbors (f2c)
    ${GPU_PARALLEL_LOOP_GANG("collapse(4) private(igrid, iib1,iib2,iib3, ixFimin1,ixFimin2,ixFimin3,ixFimax1,ixFimax2,ixFimax3, dxFi1,dxFi2,dxFi3, dxCo1,dxCo2,dxCo3, invdxCo1,invdxCo2,invdxCo3, xFimin1,xFimin2,xFimin3, xComin1,xComin2,xComin3)")}$
    do iigrid=1, igridstail
       !      inline variant of call gc_prolong(igrid)
       do i3 = -1, 1
          do i2 = -1, 1
             do i1 = -1, 1
                if (${SKIP_DIRECTION('i1','i2','i3')}$) cycle
                igrid = igrids(iigrid)
                if (neighbor_type(i1,i2,i3,igrid)==neighbor_coarse) then
                   iib1=idphyb(1,igrid);iib2=idphyb(2,igrid);iib3=idphyb(3,igrid)
                   !     inline variant of call bc_prolong(igrid,i1,i2,i3,iib1,iib2,iib3)

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

                   ${GPU_LOOP_VECTOR("collapse(3) private(slope, xFi1,xFi2,xFi3, ixCo1,ixCo2,ixCo3, xCo1,xCo2,xCo3, eta1,eta2,eta3, iw, idims, hxCo1,hxCo2,hxCo3, jxCo1,jxCo2,jxCo3, slopeL,slopeR,slopeC,signR,signC)")}$ ${GPU_INDEPENDENT()}$
                   do ixFi3 = ixFimin3,ixFimax3
                      do ixFi2 = ixFimin2,ixFimax2
                         do ixFi1 = ixFimin1,ixFimax1
                            xFi3=xFimin3+(dble(ixFi3)-half)*dxFi3
                            ixCo3=int((xFi3-xComin3)*invdxCo3)+1
                            xCo3=xComin3+(dble(ixCo3)-half)*dxCo3

                            xFi2=xFimin2+(dble(ixFi2)-half)*dxFi2
                            ixCo2=int((xFi2-xComin2)*invdxCo2)+1
                            xCo2=xComin2+(dble(ixCo2)-half)*dxCo2

                            xFi1=xFimin1+(dble(ixFi1)-half)*dxFi1
                            ixCo1=int((xFi1-xComin1)*invdxCo1)+1
                            xCo1=xComin1+(dble(ixCo1)-half)*dxCo1

                            eta1 = (xFi1-xCo1) * invdxCo1
                            eta2 = (xFi2-xCo2) * invdxCo2
                            eta3 = (xFi3-xCo3) * invdxCo3

                            do iw = nwhead, nwtail
                               do idims = 1, ndim
                                  hxCo1=ixCo1-kr(1,idims)
                                  hxCo2=ixCo2-kr(2,idims)
                                  hxCo3=ixCo3-kr(3,idims)
                                  jxCo1=ixCo1+kr(1,idims)
                                  jxCo2=ixCo2+kr(2,idims)
                                  jxCo3=ixCo3+kr(3,idims)

                                  slopeL = bgc(1)%w(ixCo1,ixCo2,ixCo3,iw,igrid) &
                                       - bgc(1)%w(hxCo1,hxCo2,hxCo3,iw,igrid)

                                  slopeR = bgc(1)%w(jxCo1,jxCo2,jxCo3,iw,igrid) &
                                       - bgc(1)%w(ixCo1,ixCo2,ixCo3,iw,igrid)

                                  slopeC = half * ( slopeR + slopeL )

                                  ! get limited slope
                                  signR = sign(one,slopeR)
                                  signC = sign(one,slopeC)

                                  slope(idims) = signC * max(zero,min(dabs(slopeC),&
                                       signC*slopeL,signC*slopeR))
                               end do

                               ! Interpolate from coarse cell using limited slopes
                               bg(bgstep)%w(ixFi1,ixFi2,ixFi3,iw,igrid) = &
                                    bgc(1)%w(ixCo1,ixCo2,ixCo3,iw, igrid) &
                                    + (slope(1)*eta1) + (slope(2)*eta2) + (slope(3)*eta3)
                            end do

                         end do
                      end do
                   end do

                end if
             end do
          end do
       end do
    end do

    ! Re-fill the physical-boundary ghost cells now that the internal
    ! ghost-cell exchange has run, mirroring upstream MPI-AMRVAC's
    ! post-exchange bc_phys pass. fill_boundary_before_gc fills only the
    ! interior transverse range of each boundary (corner info is not yet
    ! known); fill_boundary_after_gc enlarges it to the corners, which is
    ! where a physical boundary meets a pole face. No current test needs
    ! this second pass, but it keeps corner-ghost handling complete and
    ! matches upstream for setups this fork does not yet exercise.
    if(bcphys.and. .not.stagger_grid) then
       ${GPU_PARALLEL_LOOP_GANG("private(igrid)")}$
       do iigrid = 1, igridstail; igrid=igrids(iigrid);
          if (.not.phyboundblock(igrid)) cycle
          call fill_boundary_after_gc(psb(igrid),igrid,time,qdt)
       end do
    end if

    ! modify normal component of magnetic field to fix divB=0
    if(bcphys.and.associated(phys_boundary_adjust)) then
       do iigrid=1,igridstail; igrid=igrids(iigrid);
          if(.not.phyboundblock(igrid)) cycle
          call phys_boundary_adjust(igrid,psb)
       end do
    end if

    time_bc=time_bc+(MPI_WTIME()-time_bcin)

    call nvtxEndRange

  end subroutine getbc


  subroutine fill_coarse_boundary(time,igrid,i1,i2,i3)
    use mod_global_parameters
    use mod_boundary_conditions, only: bc_phys
    ${GPU_ROUTINE_VECTOR()}$
    integer, intent(in) :: igrid,i1,i2,i3
    double precision, intent(in) :: time

    integer :: idims,iside,kmin1,kmin2,kmin3,kmax1,kmax2,kmax3,ixBmin1,&
         ixBmin2,ixBmin3,ixBmax1,ixBmax2,ixBmax3,ii1,ii2,ii3

    if(phyboundblock(igrid).and..not.stagger_grid.and.bcphys) then
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
             call bc_phys(iside,idims,time,0.d0,psc(igrid),&
                  bgeoc%x(:,:,:,:,igrid),ixCoGmin1,&
                  ixCoGmin2,ixCoGmin3,ixCoGmax1,ixCoGmax2,ixCoGmax3,ixBmin1,&
                  ixBmin2,ixBmin3,ixBmax1,ixBmax2,ixBmax3)
          end do
       end do
    end if

  end subroutine fill_coarse_boundary

  !> Physical boundary conditions
  subroutine fill_boundary_before_gc(s,igrid,time,qdt)
    use mod_global_parameters
    use mod_boundary_conditions, only: bc_phys
    ${GPU_ROUTINE_VECTOR()}$

    type(state), intent(inout) :: s
    integer, intent(in) :: igrid
    double precision, intent(in) :: time, qdt

    integer :: idims,iside,i1,i2,i3,kmin1,kmin2,kmin3,kmax1,kmax2,kmax3,&
         ixBmin1,ixBmin2,ixBmin3,ixBmax1,ixBmax2,ixBmax3

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
                  s%is_physical_boundary(2*idims-2+iside)) cycle
          else
             if (neighbor_type(i1,i2,i3,igrid) /= neighbor_boundary) cycle
          end if
          call bc_phys(iside,idims,time,qdt,s,bgeo%x(:,:,:,:,igrid),&
               ixGlo1,ixGlo2,ixGlo3,&
               ixGhi1,ixGhi2,ixGhi3,ixBmin1,ixBmin2,ixBmin3,ixBmax1,ixBmax2,&
               ixBmax3)
       end do
    end do

  end subroutine fill_boundary_before_gc

  !> Physical boundary conditions
  subroutine fill_boundary_after_gc(s,igrid,time,qdt)
    use mod_global_parameters
    use mod_boundary_conditions, only: bc_phys
    ${GPU_ROUTINE_VECTOR()}$

    type(state), intent(inout) :: s
    integer, intent(in) :: igrid
    double precision, intent(in) :: time, qdt

    integer :: idims,iside,i1,i2,i3,kmin1,kmin2,kmin3,kmax1,kmax2,kmax3,&
         ixBmin1,ixBmin2,ixBmin3,ixBmax1,ixBmax2,ixBmax3

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
                  s%is_physical_boundary(2*idims-2+iside)) cycle
          else
             if (neighbor_type(i1,i2,i3,igrid) /= neighbor_boundary) cycle
          end if
          call bc_phys(iside,idims,time,qdt,s,bgeo%x(:,:,:,:,igrid),&
               ixGlo1,ixGlo2,ixGlo3,&
               ixGhi1,ixGhi2,ixGhi3,ixBmin1,ixBmin2,ixBmin3,ixBmax1,ixBmax2,&
               ixBmax3)
       end do
    end do

  end subroutine fill_boundary_after_gc

end module mod_ghostcells_update
