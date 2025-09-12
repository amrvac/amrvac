!> update ghost cells of all blocks including physical boundaries 
module mod_ghostcells_update

  implicit none

  ! The number of interleaving sending buffers for ghost cells
  integer, parameter :: npwbuf=2

  integer :: ixM^L, ixCoG^L, ixCoM^L, ixCoGs^L

  ! index ranges to send (S) to sibling blocks, receive (R) from sibling blocks
  integer, dimension(-1:1) :: ixS_srl_^L, ixR_srl_^L

  ! index ranges of staggered variables to send (S) to sibling blocks, receive (R) from sibling blocks
  integer, dimension(^ND,-1:1) :: ixS_srl_stg_^L, ixR_srl_stg_^L

  ! index ranges to send (S) restricted (r) ghost cells to coarser blocks 
  integer, dimension(-1:1) :: ixS_r_^L

  ! index ranges of staggered variables to send (S) restricted (r) ghost cells to coarser blocks 
  integer, dimension(^ND,-1:1) :: ixS_r_stg_^L

  ! index ranges to receive restriced ghost cells from finer blocks 
  integer, dimension(0:3) :: ixR_r_^L

  ! index ranges of staggered variables to receive restriced ghost cells from finer blocks 
  integer, dimension(^ND,0:3)  :: ixR_r_stg_^L

  ! send prolongated (p) ghost cells to finer blocks, receive prolongated from coarser blocks
  integer, dimension(0:3) :: ixS_p_^L, ixR_p_^L

  ! send prolongated (p) staggered ghost cells to finer blocks, receive prolongated from coarser blocks
  integer, dimension(^ND,0:3)  :: ixS_p_stg_^L, ixR_p_stg_^L

  ! number of MPI receive-send pairs, srl: same refinement level; r: restrict; p: prolong
  integer :: nrecv_bc_srl, nsend_bc_srl, nrecv_bc_r, nsend_bc_r, nrecv_bc_p, nsend_bc_p

  ! record index position of buffer arrays
  integer :: ibuf_send_srl, ibuf_recv_srl, ibuf_send_r, ibuf_recv_r, ibuf_send_p, ibuf_recv_p

  ! count of times of send and receive
  integer :: isend_srl, irecv_srl, isend_r, irecv_r, isend_p, irecv_p

  ! count of times of send and receive for cell center ghost cells
  integer :: isend_c, irecv_c

  ! tag of MPI send and recv
  integer, private :: itag

  ! total sizes = cell-center normal flux + stagger-grid flux of send and receive
  integer, dimension(-1:1^D&) :: sizes_srl_send_total, sizes_srl_recv_total

  ! sizes of buffer arrays for center-grid variable for siblings and restrict
  integer, dimension(:), allocatable :: recvrequest_c_sr, sendrequest_c_sr
  integer, dimension(:,:), allocatable :: recvstatus_c_sr, sendstatus_c_sr

  ! sizes of buffer arrays for center-grid variable for prolongation
  integer, dimension(:), allocatable :: recvrequest_c_p, sendrequest_c_p
  integer, dimension(:,:), allocatable :: recvstatus_c_p, sendstatus_c_p

  ! sizes of buffer arrays for stagger-grid variable
  integer, dimension(^ND,-1:1^D&) :: sizes_srl_send_stg, sizes_srl_recv_stg

  integer, dimension(:), allocatable :: recvrequest_srl, sendrequest_srl
  integer, dimension(:,:), allocatable :: recvstatus_srl, sendstatus_srl
 
  ! buffer arrays for send and receive of siblings, allocate in build_connectivity
  double precision, dimension(:), allocatable :: recvbuffer_srl, sendbuffer_srl

  integer, dimension(:), allocatable :: recvrequest_r, sendrequest_r
  integer, dimension(:,:), allocatable :: recvstatus_r, sendstatus_r

  ! buffer arrays for send and receive in restriction
  double precision, dimension(:), allocatable :: recvbuffer_r, sendbuffer_r

  integer, dimension(:), allocatable :: recvrequest_p, sendrequest_p
  integer, dimension(:,:), allocatable :: recvstatus_p, sendstatus_p

  ! buffer arrays for send and receive in prolongation
  double precision, dimension(:), allocatable :: recvbuffer_p, sendbuffer_p

  ! sizes to allocate buffer arrays for send and receive for restriction
  integer, dimension(-1:1^D&)     :: sizes_r_send_total
  integer, dimension(0:3^D&)      :: sizes_r_recv_total
  integer, dimension(^ND,-1:1^D&) :: sizes_r_send_stg
  integer, dimension(^ND,0:3^D&)  :: sizes_r_recv_stg

  ! sizes to allocate buffer arrays for send and receive for restriction
  integer, dimension(0:3^D&)      :: sizes_p_send_total, sizes_p_recv_total
  integer, dimension(^ND,0:3^D&)  :: sizes_p_send_stg, sizes_p_recv_stg

  ! There are two variants, _f indicates that all flux variables are filled,
  ! whereas _p means that part of the variables is filled 
  ! Furthermore _r_ stands for restrict, _p_ for prolongation.
  integer, dimension(-1:1^D&), target :: type_send_srl_f, type_recv_srl_f
  integer, dimension(-1:1^D&), target :: type_send_r_f
  integer, dimension( 0:3^D&), target :: type_recv_r_f, type_send_p_f, type_recv_p_f
  integer, dimension(-1:1^D&), target :: type_send_srl_p1, type_recv_srl_p1
  integer, dimension(-1:1^D&), target :: type_send_r_p1
  integer, dimension( 0:3^D&), target :: type_recv_r_p1, type_send_p_p1, type_recv_p_p1
  integer, dimension(-1:1^D&), target :: type_send_srl_p2, type_recv_srl_p2
  integer, dimension(-1:1^D&), target :: type_send_r_p2
  integer, dimension( 0:3^D&), target :: type_recv_r_p2, type_send_p_p2, type_recv_p_p2
  integer, dimension(  :^D&), pointer :: type_send_srl, type_recv_srl, type_send_r
  integer, dimension(  :^D&), pointer :: type_recv_r, type_send_p, type_recv_p

  ! A switch of update physical boundary or not
  logical, public :: bcphys=.true.
  ! Special buffer for pole copy
  type wbuffer
    double precision, dimension(:^D&,:), allocatable :: w
  end type wbuffer

contains

  subroutine init_bc()
    use mod_global_parameters
    use mod_physics, only: physics_type
    use mod_comm_lib, only: mpistop

    integer :: nghostcellsCo, interpolation_order
    integer :: nx^D, nxCo^D, ixG^L, i^D, idir

    ixG^L=ixG^LL;
    ixM^L=ixG^L^LSUBnghostcells;
    ixCoGmin^D=1;
    ixCoGmax^D=(ixGhi^D-2*nghostcells)/2+2*nghostcells;
    ixCoGsmin^D=0;
    ixCoGsmax^D=ixCoGmax^D;

    ixCoM^L=ixCoG^L^LSUBnghostcells;

    nx^D=ixMmax^D-ixMmin^D+1;
    nxCo^D=nx^D/2;

    if(ghost_copy) then
       interpolation_order=1
    else
       interpolation_order=2
    end if
    nghostcellsCo=int((nghostcells+1)/2)
    
    if (nghostcellsCo+interpolation_order-1>nghostcells) then
       call mpistop("interpolation order for prolongation in getbc too high")
    end if

    ! i=-1 means subregion prepared for the neighbor at its minimum side 
    ! i= 1 means subregion prepared for the neighbor at its maximum side 
    {
    ixS_srl_min^D(-1)=ixMmin^D
    ixS_srl_min^D( 0)=ixMmin^D
    ixS_srl_min^D( 1)=ixMmax^D+1-nghostcells
    ixS_srl_max^D(-1)=ixMmin^D-1+nghostcells
    ixS_srl_max^D( 0)=ixMmax^D
    ixS_srl_max^D( 1)=ixMmax^D
     
    ixR_srl_min^D(-1)=1
    ixR_srl_min^D( 0)=ixMmin^D
    ixR_srl_min^D( 1)=ixMmax^D+1
    ixR_srl_max^D(-1)=nghostcells
    ixR_srl_max^D( 0)=ixMmax^D
    ixR_srl_max^D( 1)=ixGmax^D
    
    ixS_r_min^D(-1)=ixCoMmin^D
    ixS_r_min^D( 0)=ixCoMmin^D
    ixS_r_min^D( 1)=ixCoMmax^D+1-nghostcells
    ixS_r_max^D(-1)=ixCoMmin^D-1+nghostcells
    ixS_r_max^D( 0)=ixCoMmax^D
    ixS_r_max^D( 1)=ixCoMmax^D
    
    ixR_r_min^D(0)=1
    ixR_r_min^D(1)=ixMmin^D
    ixR_r_min^D(2)=ixMmin^D+nxCo^D
    ixR_r_min^D(3)=ixMmax^D+1
    ixR_r_max^D(0)=nghostcells
    ixR_r_max^D(1)=ixMmin^D-1+nxCo^D
    ixR_r_max^D(2)=ixMmax^D
    ixR_r_max^D(3)=ixGmax^D

    ixS_p_min^D(0)=ixMmin^D-(interpolation_order-1)
    ixS_p_min^D(1)=ixMmin^D-(interpolation_order-1)
    ixS_p_min^D(2)=ixMmin^D+nxCo^D-nghostcellsCo-(interpolation_order-1)
    ixS_p_min^D(3)=ixMmax^D+1-nghostcellsCo-(interpolation_order-1)
    ixS_p_max^D(0)=ixMmin^D-1+nghostcellsCo+(interpolation_order-1)
    ixS_p_max^D(1)=ixMmin^D-1+nxCo^D+nghostcellsCo+(interpolation_order-1)
    ixS_p_max^D(2)=ixMmax^D+(interpolation_order-1)
    ixS_p_max^D(3)=ixMmax^D+(interpolation_order-1)

    ixR_p_min^D(0)=ixCoMmin^D-nghostcellsCo-(interpolation_order-1)
    ixR_p_min^D(1)=ixCoMmin^D-(interpolation_order-1)
    ixR_p_min^D(2)=ixCoMmin^D-nghostcellsCo-(interpolation_order-1)
    ixR_p_min^D(3)=ixCoMmax^D+1-(interpolation_order-1)
    ixR_p_max^D(0)=nghostcells+(interpolation_order-1)
    ixR_p_max^D(1)=ixCoMmax^D+nghostcellsCo+(interpolation_order-1)
    ixR_p_max^D(2)=ixCoMmax^D+(interpolation_order-1)
    ixR_p_max^D(3)=ixCoMmax^D+nghostcellsCo+(interpolation_order-1)

    \}

    if (stagger_grid) then
      allocate(pole_buf%ws(ixGs^T,nws))
      ! Staggered (face-allocated) variables
      do idir=1,ndim
      { ixS_srl_stg_min^D(idir,-1)=ixMmin^D-kr(idir,^D)
        ixS_srl_stg_max^D(idir,-1)=ixMmin^D-1+nghostcells
        ixS_srl_stg_min^D(idir,0) =ixMmin^D-kr(idir,^D)
        ixS_srl_stg_max^D(idir,0) =ixMmax^D
        ixS_srl_stg_min^D(idir,1) =ixMmax^D-nghostcells+1-kr(idir,^D)
        ixS_srl_stg_max^D(idir,1) =ixMmax^D
        
        ixR_srl_stg_min^D(idir,-1)=1-kr(idir,^D)
        ixR_srl_stg_max^D(idir,-1)=nghostcells
        ixR_srl_stg_min^D(idir,0) =ixMmin^D-kr(idir,^D)
        ixR_srl_stg_max^D(idir,0) =ixMmax^D
        ixR_srl_stg_min^D(idir,1) =ixMmax^D+1-kr(idir,^D)
        ixR_srl_stg_max^D(idir,1) =ixGmax^D

        ixS_r_stg_min^D(idir,-1)=ixCoMmin^D-kr(idir,^D)
        ixS_r_stg_max^D(idir,-1)=ixCoMmin^D-1+nghostcells
        ixS_r_stg_min^D(idir,0) =ixCoMmin^D-kr(idir,^D)
        ixS_r_stg_max^D(idir,0) =ixCoMmax^D
        ixS_r_stg_min^D(idir,1) =ixCoMmax^D+1-nghostcells-kr(idir,^D)
        ixS_r_stg_max^D(idir,1) =ixCoMmax^D
 
        ixR_r_stg_min^D(idir,0)=1-kr(idir,^D)
        ixR_r_stg_max^D(idir,0)=nghostcells
        ixR_r_stg_min^D(idir,1)=ixMmin^D-kr(idir,^D)
        ixR_r_stg_max^D(idir,1)=ixMmin^D-1+nxCo^D
        ixR_r_stg_min^D(idir,2)=ixMmin^D+nxCo^D-kr(idir,^D)
        ixR_r_stg_max^D(idir,2)=ixMmax^D
        ixR_r_stg_min^D(idir,3)=ixMmax^D+1-kr(idir,^D)
        ixR_r_stg_max^D(idir,3)=ixGmax^D
        \}
       {if (idir==^D) then
          ! Parallel components
          {
          ixS_p_stg_min^D(idir,0)=ixMmin^D-1 ! -1 to make redundant 
          ixS_p_stg_max^D(idir,0)=ixMmin^D-1+nghostcellsCo
          ixS_p_stg_min^D(idir,1)=ixMmin^D-1 ! -1 to make redundant 
          ixS_p_stg_max^D(idir,1)=ixMmin^D-1+nxCo^D+nghostcellsCo
          ixS_p_stg_min^D(idir,2)=ixMmax^D-nxCo^D-nghostcellsCo
          ixS_p_stg_max^D(idir,2)=ixMmax^D
          ixS_p_stg_min^D(idir,3)=ixMmax^D-nghostcellsCo
          ixS_p_stg_max^D(idir,3)=ixMmax^D

          ixR_p_stg_min^D(idir,0)=ixCoMmin^D-1-nghostcellsCo
          ixR_p_stg_max^D(idir,0)=ixCoMmin^D-1
          ixR_p_stg_min^D(idir,1)=ixCoMmin^D-1 ! -1 to make redundant 
          ixR_p_stg_max^D(idir,1)=ixCoMmax^D+nghostcellsCo
          ixR_p_stg_min^D(idir,2)=ixCoMmin^D-1-nghostcellsCo
          ixR_p_stg_max^D(idir,2)=ixCoMmax^D
          ixR_p_stg_min^D(idir,3)=ixCoMmax^D+1-1 ! -1 to make redundant 
          ixR_p_stg_max^D(idir,3)=ixCoMmax^D+nghostcellsCo
          \}
        else
          {
          ! Perpendicular component
          ixS_p_stg_min^D(idir,0)=ixMmin^D
          ixS_p_stg_max^D(idir,0)=ixMmin^D-1+nghostcellsCo+(interpolation_order-1)
          ixS_p_stg_min^D(idir,1)=ixMmin^D
          ixS_p_stg_max^D(idir,1)=ixMmin^D-1+nxCo^D+nghostcellsCo+(interpolation_order-1)
          ixS_p_stg_min^D(idir,2)=ixMmax^D+1-nxCo^D-nghostcellsCo-(interpolation_order-1)
          ixS_p_stg_max^D(idir,2)=ixMmax^D
          ixS_p_stg_min^D(idir,3)=ixMmax^D+1-nghostcellsCo-(interpolation_order-1)
          ixS_p_stg_max^D(idir,3)=ixMmax^D
 
          ixR_p_stg_min^D(idir,0)=ixCoMmin^D-nghostcellsCo-(interpolation_order-1)
          ixR_p_stg_max^D(idir,0)=ixCoMmin^D-1
          ixR_p_stg_min^D(idir,1)=ixCoMmin^D
          ixR_p_stg_max^D(idir,1)=ixCoMmax^D+nghostcellsCo+(interpolation_order-1)
          ixR_p_stg_min^D(idir,2)=ixCoMmin^D-nghostcellsCo-(interpolation_order-1)
          ixR_p_stg_max^D(idir,2)=ixCoMmax^D
          ixR_p_stg_min^D(idir,3)=ixCoMmax^D+1
          ixR_p_stg_max^D(idir,3)=ixCoMmax^D+nghostcellsCo+(interpolation_order-1)
          \}
        end if
        }
      end do
      ! calculate sizes for buffer arrays for siblings
      {do i^DB=-1,1\}
         ! Staggered (face-allocated) variables
         do idir=1,ndim
           sizes_srl_send_stg(idir,i^D)={(ixS_srl_stg_max^D(idir,i^D)-ixS_srl_stg_min^D(idir,i^D)+1)|*}
           sizes_srl_recv_stg(idir,i^D)={(ixR_srl_stg_max^D(idir,i^D)-ixR_srl_stg_min^D(idir,i^D)+1)|*}
           sizes_r_send_stg(idir,i^D)={(ixS_r_stg_max^D(idir,i^D)-ixS_r_stg_min^D(idir,i^D)+1)|*}
         end do
         sizes_srl_send_total(i^D)=sum(sizes_srl_send_stg(:,i^D))
         sizes_srl_recv_total(i^D)=sum(sizes_srl_recv_stg(:,i^D))
         sizes_r_send_total(i^D)=sum(sizes_r_send_stg(:,i^D))
      {end do\}

      {do i^DB=0,3\}
         ! Staggered (face-allocated) variables
           do idir=1,ndim
             sizes_r_recv_stg(idir,i^D)={(ixR_r_stg_max^D(idir,i^D)-ixR_r_stg_min^D(idir,i^D)+1)|*}
             sizes_p_send_stg(idir,i^D)={(ixS_p_stg_max^D(idir,i^D)-ixS_p_stg_min^D(idir,i^D)+1)|*}
             sizes_p_recv_stg(idir,i^D)={(ixR_p_stg_max^D(idir,i^D)-ixR_p_stg_min^D(idir,i^D)+1)|*}
           end do
           sizes_r_recv_total(i^D)=sum(sizes_r_recv_stg(:,i^D))
           sizes_p_send_total(i^D)=sum(sizes_p_send_stg(:,i^D))
           sizes_p_recv_total(i^D)=sum(sizes_p_recv_stg(:,i^D))
      {end do\}
    end if

  end subroutine init_bc

  subroutine create_bc_mpi_datatype(nwstart,nwbc) 
    use mod_global_parameters

    integer, intent(in) :: nwstart, nwbc
    integer :: i^D, ic^D, inc^D

   {do i^DB=-1,1\}
      if(i^D==0|.and.) cycle
      call get_bc_comm_type(type_send_srl(i^D),ixS_srl_^L(i^D),ixG^LL,nwstart,nwbc)
      call get_bc_comm_type(type_recv_srl(i^D),ixR_srl_^L(i^D),ixG^LL,nwstart,nwbc)
      call get_bc_comm_type(type_send_r(i^D),   ixS_r_^L(i^D),ixCoG^L,nwstart,nwbc)
      {do ic^DB=1+int((1-i^DB)/2),2-int((1+i^DB)/2)
         inc^DB=2*i^DB+ic^DB\}
         call get_bc_comm_type(type_recv_r(inc^D),ixR_r_^L(inc^D), ixG^LL,nwstart,nwbc)
         call get_bc_comm_type(type_send_p(inc^D),ixS_p_^L(inc^D), ixG^LL,nwstart,nwbc)
         call get_bc_comm_type(type_recv_p(inc^D),ixR_p_^L(inc^D),ixCoG^L,nwstart,nwbc)
      {end do\}
   {end do\}
  
  end subroutine create_bc_mpi_datatype

  subroutine get_bc_comm_type(comm_type,ix^L,ixG^L,nwstart,nwbc)
    use mod_global_parameters
  
    integer, intent(inout) :: comm_type
    integer, intent(in) :: ix^L, ixG^L, nwstart, nwbc
    
    integer, dimension(ndim+1) :: fullsize, subsize, start

    ^D&fullsize(^D)=ixGmax^D;
    fullsize(ndim+1)=nw
    ^D&subsize(^D)=ixmax^D-ixmin^D+1;
    subsize(ndim+1)=nwbc
    ^D&start(^D)=ixmin^D-1;
    start(ndim+1)=nwstart-1
    
    call MPI_TYPE_CREATE_SUBARRAY(ndim+1,fullsize,subsize,start,MPI_ORDER_FORTRAN, &
                                  MPI_DOUBLE_PRECISION,comm_type,ierrmpi)
    call MPI_TYPE_COMMIT(comm_type,ierrmpi)
    
  end subroutine get_bc_comm_type

  !> do update ghost cells of all blocks including physical boundaries
  subroutine getbc(time,qdt,psb,nwstart,nwbc)
    use mod_global_parameters
    use mod_physics
    use mod_coarsen, only: coarsen_grid
    use mod_boundary_conditions, only: getintbc, bc_phys
    use mod_comm_lib, only: mpistop

    double precision, intent(in)      :: time, qdt
    type(state), target               :: psb(max_blocks)
    integer, intent(in)               :: nwstart ! Fill from nwstart
    integer, intent(in)               :: nwbc    ! Number of variables to fill

    double precision :: time_bcin
    integer :: nwhead, nwtail
    integer :: iigrid, igrid, isizes, i^D
    integer :: isend_buf(npwbuf), ipwbuf, nghostcellsco
    ! index pointer for buffer arrays as a start for a segment
    integer :: ibuf_start, ibuf_next
    ! shapes of reshape
    integer, dimension(1) :: shapes
    type(wbuffer) :: pwbuf(npwbuf)

    time_bcin=MPI_WTIME()

    nwhead=nwstart
    nwtail=nwstart+nwbc-1

    ! fill internal physical boundary
    if (internalboundary) then 
       call getintbc(time,ixG^LL)
    end if

    ! prepare coarse values to send to coarser neighbors
    !$OMP PARALLEL DO SCHEDULE(dynamic) PRIVATE(igrid)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
      if(any(neighbor_type(:^D&,igrid)==neighbor_coarse)) then
        call coarsen_grid(psb(igrid),ixG^LL,ixM^L,psc(igrid),ixCoG^L,ixCoM^L)
      end if
    end do
    !$OMP END PARALLEL DO

    ! default : no singular axis
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

    !$OMP PARALLEL SECTIONS PRIVATE(iigrid,igrid,i^D)
    !$OMP SECTION
    ! MPI receive ghost-cell values from sibling blocks and finer neighbors in different processors
    do iigrid=1,igridstail; igrid=igrids(iigrid);
      {do i^DB=-1,1\}
         if (skip_direction([ i^D ])) cycle
         select case (neighbor_type(i^D,igrid))
         case (neighbor_sibling)
            call bc_recv_srl(igrid,i^D)
         case (neighbor_fine)
            call bc_recv_restrict(igrid,i^D)
         end select
      {end do\}
    end do

    ! MPI send ghost-cell values to sibling blocks and coarser neighbors in different processors
    do iigrid=1,igridstail; igrid=igrids(iigrid);
      {do i^DB=-1,1\}
         if(skip_direction([ i^D ])) cycle
         select case (neighbor_type(i^D,igrid))
         case (neighbor_sibling)
            call bc_send_srl(igrid,i^D)
         case (neighbor_coarse)
            call bc_send_restrict(igrid,i^D)
         end select
      {end do\}
    end do

    call MPI_WAITALL(irecv_c,recvrequest_c_sr,recvstatus_c_sr,ierrmpi)
    call MPI_WAITALL(isend_c,sendrequest_c_sr,sendstatus_c_sr,ierrmpi)
    if(stagger_grid) then
      call MPI_WAITALL(nrecv_bc_srl,recvrequest_srl,recvstatus_srl,ierrmpi)
      call MPI_WAITALL(nsend_bc_srl,sendrequest_srl,sendstatus_srl,ierrmpi)
      call MPI_WAITALL(nrecv_bc_r,recvrequest_r,recvstatus_r,ierrmpi)
      call MPI_WAITALL(nsend_bc_r,sendrequest_r,sendstatus_r,ierrmpi)
    end if

    do ipwbuf=1,npwbuf
       if (isend_buf(ipwbuf)/=0) deallocate(pwbuf(ipwbuf)%w)
    end do
    !$OMP SECTION
    ! fill ghost-cell values of sibling blocks and coarser neighbors in the same processor
    do iigrid=1,igridstail; igrid=igrids(iigrid);
      {do i^DB=-1,1\}
        if(skip_direction([ i^D ])) cycle
         select case (neighbor_type(i^D,igrid))
         case(neighbor_sibling)
           call bc_fill_srl(igrid,i^D)
         case(neighbor_coarse)
           call bc_fill_restrict(igrid,i^D)
         end select
      {end do\}
    end do
    !$OMP END PARALLEL SECTIONS

    if(stagger_grid) then
      ! unpack the received data from sibling blocks and finer neighbors to fill ghost-cell staggered values
      ibuf_recv_srl=1
      ibuf_recv_r=1
      do iigrid=1,igridstail; igrid=igrids(iigrid);
       {do i^DB=-1,1\}
          if (skip_direction([ i^D ])) cycle
          select case (neighbor_type(i^D,igrid))
          case (neighbor_sibling)
             call bc_fill_srl_stg(igrid,i^D)
          case (neighbor_fine)
             call bc_fill_restrict_stg(igrid,i^D)
          end select
       {end do\}
      end do
    end if

    irecv_c=0
    isend_c=0
    isend_buf=0
    ipwbuf=1

    !$OMP PARALLEL SECTIONS PRIVATE(iigrid,igrid,i^D)
    !$OMP SECTION
    ! MPI receive ghost-cell values from coarser neighbors in different processors
    do iigrid=1,igridstail; igrid=igrids(iigrid);
      {do i^DB=-1,1\}
         if (skip_direction([ i^D ])) cycle
         if (neighbor_type(i^D,igrid)==neighbor_coarse) call bc_recv_prolong(igrid,i^D)
      {end do\}
    end do
    ! MPI send ghost-cell values to finer neighbors in different processors
    do iigrid=1,igridstail; igrid=igrids(iigrid);
      {do i^DB=-1,1\}
         if (skip_direction([ i^D ])) cycle
         if (neighbor_type(i^D,igrid)==neighbor_fine) call bc_send_prolong(igrid,i^D)
      {end do\}
    end do

    call MPI_WAITALL(irecv_c,recvrequest_c_p,recvstatus_c_p,ierrmpi)
    call MPI_WAITALL(isend_c,sendrequest_c_p,sendstatus_c_p,ierrmpi)

    if(stagger_grid) then
      call MPI_WAITALL(nrecv_bc_p,recvrequest_p,recvstatus_p,ierrmpi)
      call MPI_WAITALL(nsend_bc_p,sendrequest_p,sendstatus_p,ierrmpi)
    end if

    do ipwbuf=1,npwbuf
       if (isend_buf(ipwbuf)/=0) deallocate(pwbuf(ipwbuf)%w)
    end do

    !$OMP SECTION
    ! fill coarse ghost-cell values of finer neighbors in the same processor
    do iigrid=1,igridstail; igrid=igrids(iigrid);
      {do i^DB=-1,1\}
         if (skip_direction([ i^D ])) cycle
         if (neighbor_type(i^D,igrid)==neighbor_fine) call bc_fill_prolong(igrid,i^D)
      {end do\}
    end do
    !$OMP END PARALLEL SECTIONS

    if(stagger_grid) then
      ! fill coarser representative ghost cells after receipt
      ibuf_recv_p=1
      do iigrid=1,igridstail; igrid=igrids(iigrid);
        {do i^DB=-1,1\}
           if (skip_direction([ i^D ])) cycle
           if(neighbor_type(i^D,igrid)==neighbor_coarse) call bc_fill_prolong_stg(igrid,i^D)
        {end do\}
      end do
    end if
    ! do prolongation on the ghost-cell values based on the received coarse values from coarser neighbors
    !$OMP PARALLEL DO SCHEDULE(dynamic) PRIVATE(igrid)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
      call gc_prolong(igrid)
    end do
    !$OMP END PARALLEL DO

    ! fill physical boundary ghost cells after internal ghost-cell values exchange
    if(bcphys) then
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
    
    contains

      logical function skip_direction(dir)
        integer, intent(in) :: dir(^ND)

        if (all(dir == 0)) then
           skip_direction = .true.
        else
           skip_direction = .false.
        end if
      end function skip_direction

      !> Physical boundary conditions
      subroutine fill_boundary_after_gc(igrid)

        integer, intent(in) :: igrid

        integer :: idims,iside,i^D,k^L,ixB^L

        block=>psb(igrid)
        ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
        do idims=1,ndim
          ! to avoid using as yet unknown corner info in more than 1D, we
          ! fill only interior mesh ranges of the ghost cell ranges at first,
          ! and progressively enlarge the ranges to include corners later
          kmin1=0; kmax1=0;
          {^IFTWOD
           kmin2=merge(1, 0,  idims .lt. 2 .and. neighbor_type(0,-1,igrid)==1)
           kmax2=merge(1, 0,  idims .lt. 2 .and. neighbor_type(0, 1,igrid)==1)}
          {^IFTHREED
           kmin2=merge(1, 0, idims .lt. 2 .and. neighbor_type(0,-1,0,igrid)==1)
           kmax2=merge(1, 0, idims .lt. 2 .and. neighbor_type(0, 1,0,igrid)==1)
           kmin3=merge(1, 0, idims .lt. 3 .and. neighbor_type(0,0,-1,igrid)==1)
           kmax3=merge(1, 0, idims .lt. 3 .and. neighbor_type(0,0, 1,igrid)==1)}
          ixBmin^D=ixGlo^D+kmin^D*nghostcells;
          ixBmax^D=ixGhi^D-kmax^D*nghostcells;
          do iside=1,2
            i^D=kr(^D,idims)*(2*iside-3);
            if (aperiodB(idims)) then 
              if (neighbor_type(i^D,igrid) /= neighbor_boundary .and. &
                 .not. psb(igrid)%is_physical_boundary(2*idims-2+iside)) cycle
            else 
              if (neighbor_type(i^D,igrid) /= neighbor_boundary) cycle
            end if
            call bc_phys(iside,idims,time,qdt,psb(igrid),ixG^LL,ixB^L)
          end do
        end do

      end subroutine fill_boundary_after_gc

      !> MPI receive from sibling at same refinement level
      subroutine bc_recv_srl(igrid,i^D)
        integer, intent(in) :: igrid,i^D

        integer :: ipe_neighbor

        ipe_neighbor=neighbor(2,i^D,igrid)
        if (ipe_neighbor/=mype) then
           irecv_c=irecv_c+1
           itag=(3**^ND+4**^ND)*(igrid-1)+{(i^D+1)*3**(^D-1)+}
           call MPI_IRECV(psb(igrid)%w,1,type_recv_srl(i^D), &
                          ipe_neighbor,itag,icomm,recvrequest_c_sr(irecv_c),ierrmpi)
           if(stagger_grid) then
             irecv_srl=irecv_srl+1
             call MPI_IRECV(recvbuffer_srl(ibuf_recv_srl),sizes_srl_recv_total(i^D),MPI_DOUBLE_PRECISION, &
                            ipe_neighbor,itag,icomm,recvrequest_srl(irecv_srl),ierrmpi)
             ibuf_recv_srl=ibuf_recv_srl+sizes_srl_recv_total(i^D)
           end if
        end if

      end subroutine bc_recv_srl

      !> MPI receive from fine neighbor
      subroutine bc_recv_restrict(igrid,i^D)
        integer, intent(in) :: igrid,i^D

        integer :: ic^D,inc^D,ipe_neighbor

        {do ic^DB=1+int((1-i^DB)/2),2-int((1+i^DB)/2)
           inc^DB=2*i^DB+ic^DB\}
           ipe_neighbor=neighbor_child(2,inc^D,igrid)
           if (ipe_neighbor/=mype) then
              irecv_c=irecv_c+1
              itag=(3**^ND+4**^ND)*(igrid-1)+3**^ND+{inc^D*4**(^D-1)+}
              call MPI_IRECV(psb(igrid)%w,1,type_recv_r(inc^D), &
                             ipe_neighbor,itag,icomm,recvrequest_c_sr(irecv_c),ierrmpi)
              if(stagger_grid) then
                irecv_r=irecv_r+1
                call MPI_IRECV(recvbuffer_r(ibuf_recv_r),sizes_r_recv_total(inc^D), &
                               MPI_DOUBLE_PRECISION,ipe_neighbor,itag, &
                               icomm,recvrequest_r(irecv_r),ierrmpi)
                ibuf_recv_r=ibuf_recv_r+sizes_r_recv_total(inc^D)
              end if
           end if
        {end do\}

      end subroutine bc_recv_restrict

      !> MPI send to sibling at same refinement level
      subroutine bc_send_srl(igrid,i^D)
        integer, intent(in) :: igrid,i^D

        integer :: n_i^D,ipole,idir,ineighbor,ipe_neighbor,ixS^L

        ipe_neighbor=neighbor(2,i^D,igrid)
        if(ipe_neighbor/=mype) then
          ineighbor=neighbor(1,i^D,igrid)
          ipole=neighbor_pole(i^D,igrid)
          if(ipole==0) then
            n_i^D=-i^D;
            isend_c=isend_c+1
            itag=(3**^ND+4**^ND)*(ineighbor-1)+{(n_i^D+1)*3**(^D-1)+}
            call MPI_ISEND(psb(igrid)%w,1,type_send_srl(i^D), &
                           ipe_neighbor,itag,icomm,sendrequest_c_sr(isend_c),ierrmpi)
            if(stagger_grid) then
              ibuf_start=ibuf_send_srl
              do idir=1,ndim
                ixS^L=ixS_srl_stg_^L(idir,i^D);
                ibuf_next=ibuf_start+sizes_srl_send_stg(idir,i^D)
                shapes=(/sizes_srl_send_stg(idir,i^D)/)
                sendbuffer_srl(ibuf_start:ibuf_next-1)=&
                  reshape(psb(igrid)%ws(ixS^S,idir),shapes)
                ibuf_start=ibuf_next
              end do
              isend_srl=isend_srl+1
              call MPI_ISEND(sendbuffer_srl(ibuf_send_srl),sizes_srl_send_total(i^D),&
                MPI_DOUBLE_PRECISION, ipe_neighbor,itag,icomm,sendrequest_srl(isend_srl),ierrmpi)
              ibuf_send_srl=ibuf_next
            end if
          else
            ixS^L=ixS_srl_^L(i^D);
            select case (ipole)
            {case (^D)
               n_i^D=i^D^D%n_i^DD=-i^DD;\}
            end select
            if (isend_buf(ipwbuf)/=0) then
               call MPI_WAIT(sendrequest_c_sr(isend_buf(ipwbuf)), &
                             sendstatus_c_sr(:,isend_buf(ipwbuf)),ierrmpi)
               deallocate(pwbuf(ipwbuf)%w)
            end if
            allocate(pwbuf(ipwbuf)%w(ixS^S,nwhead:nwtail))
            call pole_buffer(pwbuf(ipwbuf)%w,ixS^L,ixS^L,psb(igrid)%w,ixG^LL,ixS^L,ipole)
            isend_c=isend_c+1
            isend_buf(ipwbuf)=isend_c
            itag=(3**^ND+4**^ND)*(ineighbor-1)+{(n_i^D+1)*3**(^D-1)+}
            isizes={(ixSmax^D-ixSmin^D+1)*}*nwbc
            call MPI_ISEND(pwbuf(ipwbuf)%w,isizes,MPI_DOUBLE_PRECISION, &
                           ipe_neighbor,itag,icomm,sendrequest_c_sr(isend_c),ierrmpi)
            ipwbuf=1+modulo(ipwbuf,npwbuf)
            if(stagger_grid) then
              ibuf_start=ibuf_send_srl
              do idir=1,ndim
                ixS^L=ixS_srl_stg_^L(idir,i^D);
                ibuf_next=ibuf_start+sizes_srl_send_stg(idir,i^D)
                shapes=(/sizes_srl_send_stg(idir,i^D)/)
                sendbuffer_srl(ibuf_start:ibuf_next-1)=&
                  reshape(psb(igrid)%ws(ixS^S,idir),shapes)
                ibuf_start=ibuf_next
              end do
              isend_srl=isend_srl+1
              call MPI_ISEND(sendbuffer_srl(ibuf_send_srl),sizes_srl_send_total(i^D),&
                MPI_DOUBLE_PRECISION, ipe_neighbor,itag,icomm,sendrequest_srl(isend_srl),ierrmpi)
              ibuf_send_srl=ibuf_next
            end if
          end if
        end if

      end subroutine bc_send_srl

      !> MPI send to coarser neighbor's ghost cells
      subroutine bc_send_restrict(igrid,i^D)
        integer, intent(in) :: igrid,i^D

        integer :: ic^D,n_inc^D,ipole,idir,ineighbor,ipe_neighbor,ixS^L

        ipe_neighbor=neighbor(2,i^D,igrid)
        if(ipe_neighbor/=mype) then
          ic^D=1+modulo(node(pig^D_,igrid)-1,2);
          if({.not.(i^D==0.or.i^D==2*ic^D-3)|.or.}) return
          ineighbor=neighbor(1,i^D,igrid)
          ipole=neighbor_pole(i^D,igrid)
          if(ipole==0) then
            n_inc^D=-2*i^D+ic^D;
            isend_c=isend_c+1
            itag=(3**^ND+4**^ND)*(ineighbor-1)+3**^ND+{n_inc^D*4**(^D-1)+}
            call MPI_ISEND(psc(igrid)%w,1,type_send_r(i^D), &
                           ipe_neighbor,itag,icomm,sendrequest_c_sr(isend_c),ierrmpi)
            if(stagger_grid) then 
              ibuf_start=ibuf_send_r
              do idir=1,ndim
                ixS^L=ixS_r_stg_^L(idir,i^D);
                ibuf_next=ibuf_start+sizes_r_send_stg(idir,i^D)
                shapes=(/sizes_r_send_stg(idir,i^D)/)
                sendbuffer_r(ibuf_start:ibuf_next-1)=&
                  reshape(psc(igrid)%ws(ixS^S,idir),shapes)
                ibuf_start=ibuf_next
              end do
              isend_r=isend_r+1
              call MPI_ISEND(sendbuffer_r(ibuf_send_r),sizes_r_send_total(i^D),&
                             MPI_DOUBLE_PRECISION,ipe_neighbor,itag, &
                             icomm,sendrequest_r(isend_r),ierrmpi)
              ibuf_send_r=ibuf_next
            end if
          else
            ixS^L=ixS_r_^L(i^D);
            select case (ipole)
            {case (^D)
               n_inc^D=2*i^D+(3-ic^D)^D%n_inc^DD=-2*i^DD+ic^DD;\}
            end select
            if(isend_buf(ipwbuf)/=0) then
              call MPI_WAIT(sendrequest_c_sr(isend_buf(ipwbuf)), &
                            sendstatus_c_sr(:,isend_buf(ipwbuf)),ierrmpi)
              deallocate(pwbuf(ipwbuf)%w)
            end if
            allocate(pwbuf(ipwbuf)%w(ixS^S,nwhead:nwtail))
            call pole_buffer(pwbuf(ipwbuf)%w,ixS^L,ixS^L,psc(igrid)%w,ixCoG^L,ixS^L,ipole)
            isend_c=isend_c+1
            isend_buf(ipwbuf)=isend_c
            itag=(3**^ND+4**^ND)*(ineighbor-1)+3**^ND+{n_inc^D*4**(^D-1)+}
            isizes={(ixSmax^D-ixSmin^D+1)*}*nwbc
            call MPI_ISEND(pwbuf(ipwbuf)%w,isizes,MPI_DOUBLE_PRECISION, &
                           ipe_neighbor,itag,icomm,sendrequest_c_sr(isend_c),ierrmpi)
            ipwbuf=1+modulo(ipwbuf,npwbuf)
            if(stagger_grid) then 
              ibuf_start=ibuf_send_r
              do idir=1,ndim
                ixS^L=ixS_r_stg_^L(idir,i^D);
                ibuf_next=ibuf_start+sizes_r_send_stg(idir,i^D)
                shapes=(/sizes_r_send_stg(idir,i^D)/)
                sendbuffer_r(ibuf_start:ibuf_next-1)=&
                  reshape(psc(igrid)%ws(ixS^S,idir),shapes)
                ibuf_start=ibuf_next
              end do
              isend_r=isend_r+1
              call MPI_ISEND(sendbuffer_r(ibuf_send_r),sizes_r_send_total(i^D),&
                             MPI_DOUBLE_PRECISION,ipe_neighbor,itag, &
                             icomm,sendrequest_r(isend_r),ierrmpi)
              ibuf_send_r=ibuf_next
            end if
          end if
        end if

      end subroutine bc_send_restrict

      !> fill same-level neighbor's ghost cells in the same processor
      subroutine bc_fill_srl(igrid,i^D)
        integer, intent(in) :: igrid,i^D

        integer :: ineighbor,ipe_neighbor,ipole,ixS^L,ixR^L,n_i^D,idir

        ipe_neighbor=neighbor(2,i^D,igrid)
        if(ipe_neighbor==mype) then
          ineighbor=neighbor(1,i^D,igrid)
          ipole=neighbor_pole(i^D,igrid)
          if(ipole==0) then
            n_i^D=-i^D;
            ixS^L=ixS_srl_^L(i^D);
            ixR^L=ixR_srl_^L(n_i^D);
            psb(ineighbor)%w(ixR^S,nwhead:nwtail)=&
                psb(igrid)%w(ixS^S,nwhead:nwtail)
            if(stagger_grid) then
              do idir=1,ndim
                ixS^L=ixS_srl_stg_^L(idir,i^D);
                ixR^L=ixR_srl_stg_^L(idir,n_i^D);
                psb(ineighbor)%ws(ixR^S,idir)=psb(igrid)%ws(ixS^S,idir)
              end do
            end if
          else
            ixS^L=ixS_srl_^L(i^D);
            select case (ipole)
            {case (^D)
              n_i^D=i^D^D%n_i^DD=-i^DD;\}
            end select
            ixR^L=ixR_srl_^L(n_i^D);
            call pole_copy(psb(ineighbor)%w,ixG^LL,ixR^L,psb(igrid)%w,ixG^LL,ixS^L,ipole)
            if(stagger_grid) then
              do idir=1,ndim
                ixS^L=ixS_srl_stg_^L(idir,i^D);
                ixR^L=ixR_srl_stg_^L(idir,n_i^D);
                call pole_copy_stg(psb(ineighbor)%ws,ixGs^LL,ixR^L,psb(igrid)%ws,ixGs^LL,ixS^L,idir,ipole)
              end do
            end if
          end if
        end if

      end subroutine bc_fill_srl

      !> fill coarser neighbor's ghost cells in the same processor
      subroutine bc_fill_restrict(igrid,i^D)
        integer, intent(in) :: igrid,i^D

        integer :: ic^D,n_inc^D,ixS^L,ixR^L,ipe_neighbor,ineighbor,ipole,idir

        ipe_neighbor=neighbor(2,i^D,igrid)
        if(ipe_neighbor==mype) then
          ic^D=1+modulo(node(pig^D_,igrid)-1,2);
          if({.not.(i^D==0.or.i^D==2*ic^D-3)|.or.}) return
          ineighbor=neighbor(1,i^D,igrid)
          ipole=neighbor_pole(i^D,igrid)
          if(ipole==0) then
            n_inc^D=-2*i^D+ic^D;
            ixS^L=ixS_r_^L(i^D);
            ixR^L=ixR_r_^L(n_inc^D);
            psb(ineighbor)%w(ixR^S,nwhead:nwtail)=&
                psc(igrid)%w(ixS^S,nwhead:nwtail)
            if(stagger_grid) then
              do idir=1,ndim
                 ixS^L=ixS_r_stg_^L(idir,i^D);
                 ixR^L=ixR_r_stg_^L(idir,n_inc^D);
                 psb(ineighbor)%ws(ixR^S,idir)=psc(igrid)%ws(ixS^S,idir)
              end do
            end if
          else
            ixS^L=ixS_r_^L(i^D);
            select case (ipole)
            {case (^D)
              n_inc^D=2*i^D+(3-ic^D)^D%n_inc^DD=-2*i^DD+ic^DD;\}
            end select
            ixR^L=ixR_r_^L(n_inc^D);
            call pole_copy(psb(ineighbor)%w,ixG^LL,ixR^L,psc(igrid)%w,ixCoG^L,ixS^L,ipole)
            if(stagger_grid) then
              do idir=1,ndim
                ixS^L=ixS_r_stg_^L(idir,i^D);
                ixR^L=ixR_r_stg_^L(idir,n_inc^D);
                !! Fill ghost cells
                call pole_copy_stg(psb(ineighbor)%ws,ixGs^LL,ixR^L,psc(igrid)%ws,ixCoGs^L,ixS^L,idir,ipole)
              end do
            end if
          end if
        end if

      end subroutine bc_fill_restrict

      !> fill siblings ghost cells with received data
      subroutine bc_fill_srl_stg(igrid,i^D)
        integer, intent(in) :: igrid,i^D

        integer :: ixS^L,ixR^L,n_i^D,idir,ineighbor,ipe_neighbor,ipole

        ipe_neighbor=neighbor(2,i^D,igrid)
        if(ipe_neighbor/=mype) then
          ineighbor=neighbor(1,i^D,igrid)
          ipole=neighbor_pole(i^D,igrid)

        !! Now the special treatment of the pole is done here, at the receive step
          if (ipole==0) then
            ixR^L=ixR_srl_^L(i^D);
            !! Unpack the buffer and fill the ghost cells
            n_i^D=-i^D;
            do idir=1,ndim
              ixS^L=ixS_srl_stg_^L(idir,n_i^D);
              ixR^L=ixR_srl_stg_^L(idir,i^D);
              ibuf_next=ibuf_recv_srl+sizes_srl_recv_stg(idir,i^D)
              psb(igrid)%ws(ixR^S,idir)=reshape(source=recvbuffer_srl(ibuf_recv_srl:ibuf_next-1),&
                shape=shape(psb(igrid)%ws(ixS^S,idir)))       
              ibuf_recv_srl=ibuf_next
            end do
          else ! There is a pole
            select case (ipole)
            {case (^D)
               n_i^D=i^D^D%n_i^DD=-i^DD;\}
            end select
            pole_buf%ws=zero
            do idir=1,ndim
             ixR^L=ixR_srl_stg_^L(idir,i^D);
             ixS^L=ixS_srl_stg_^L(idir,n_i^D);
             ibuf_next=ibuf_recv_srl+sizes_srl_recv_stg(idir,i^D)
             pole_buf%ws(ixS^S,idir)=reshape(source=recvbuffer_srl(ibuf_recv_srl:ibuf_next-1),&
               shape=shape(psb(igrid)%ws(ixS^S,idir)))
             ibuf_recv_srl=ibuf_next
             call pole_copy_stg(psb(igrid)%ws,ixGs^LL,ixR^L,pole_buf%ws,ixGs^LL,ixS^L,idir,ipole)
            end do
          end if
        end if

      end subroutine bc_fill_srl_stg

      !> fill restricted ghost cells after receipt
      subroutine bc_fill_restrict_stg(igrid,i^D)
        integer, intent(in) :: igrid,i^D

        integer :: ipole,ic^D,inc^D,ineighbor,ipe_neighbor,ixS^L,ixR^L,n_i^D,idir

        ipole=neighbor_pole(i^D,igrid)
        if (ipole==0) then
          ! Loop over the children ic^D to and their neighbors inc^D
          {do ic^DB=1+int((1-i^DB)/2),2-int((1+i^DB)/2)
             inc^DB=2*i^DB+ic^DB\}
             ipe_neighbor=neighbor_child(2,inc^D,igrid)
             if(ipe_neighbor/=mype) then
               ineighbor=neighbor_child(1,inc^D,igrid)
               n_i^D=-i^D;
               !! Unpack the buffer and fill the ghost cells
               do idir=1,ndim
                 ixR^L=ixR_r_stg_^L(idir,inc^D);
                 ibuf_next=ibuf_recv_r+sizes_r_recv_stg(idir,inc^D)
                 psb(igrid)%ws(ixR^S,idir)=reshape(source=recvbuffer_r(ibuf_recv_r:ibuf_next-1),&
                       shape=shape(psb(igrid)%ws(ixR^S,idir)))
                 ibuf_recv_r=ibuf_next
               end do
             end if
          {end do\}
        else !! There is a pole
          {do ic^DB=1+int((1-i^DB)/2),2-int((1+i^DB)/2)
             inc^DB=2*i^DB+ic^DB\}
             ipe_neighbor=neighbor_child(2,inc^D,igrid)
             if(ipe_neighbor/=mype) then
               ineighbor=neighbor_child(1,inc^D,igrid)
               select case(ipole)
              {case (^D)
                 n_i^D=i^D^D%n_i^DD=-i^DD;\}
               end select
               ixR^L=ixR_r_^L(inc^D);
               !! Unpack the buffer and fill an auxiliary array
               pole_buf%ws=zero
               do idir=1,ndim
                 ixS^L=ixS_r_stg_^L(idir,n_i^D);
                 ixR^L=ixR_r_stg_^L(idir,inc^D);
                 ibuf_next=ibuf_recv_r+sizes_r_recv_stg(idir,inc^D)
                 pole_buf%ws(ixR^S,idir)=reshape(source=recvbuffer_r(ibuf_recv_r:ibuf_next-1),&
                   shape=shape(psb(igrid)%ws(ixR^S,idir)))
                 call pole_copy_stg(psb(igrid)%ws,ixGs^LL,ixR^L,pole_buf%ws,ixGs^LL,ixR^L,idir,ipole)
                 ibuf_recv_r=ibuf_next
               end do
             end if
          {end do\}
        end if

      end subroutine bc_fill_restrict_stg

      !> Receive from coarse neighbor
      subroutine bc_recv_prolong(igrid,i^D)
        integer, intent(in) :: igrid,i^D

        integer :: ic^D,ipe_neighbor,inc^D

        ic^D=1+modulo(node(pig^D_,igrid)-1,2);
        if ({.not.(i^D==0.or.i^D==2*ic^D-3)|.or.}) return

        ipe_neighbor=neighbor(2,i^D,igrid)
        if (ipe_neighbor/=mype) then
           irecv_c=irecv_c+1
           inc^D=ic^D+i^D;
           itag=(3**^ND+4**^ND)*(igrid-1)+3**^ND+{inc^D*4**(^D-1)+}
           call MPI_IRECV(psc(igrid)%w,1,type_recv_p(inc^D), &
                          ipe_neighbor,itag,icomm,recvrequest_c_p(irecv_c),ierrmpi)  
           if(stagger_grid) then
             irecv_p=irecv_p+1
             call MPI_IRECV(recvbuffer_p(ibuf_recv_p),sizes_p_recv_total(inc^D),& 
                            MPI_DOUBLE_PRECISION,ipe_neighbor,itag,&
                            icomm,recvrequest_p(irecv_p),ierrmpi)
             ibuf_recv_p=ibuf_recv_p+sizes_p_recv_total(inc^D)
           end if
        end if

      end subroutine bc_recv_prolong

      !> Send to finer neighbor
      subroutine bc_send_prolong(igrid,i^D)
        integer, intent(in) :: igrid,i^D

        integer :: ic^D,inc^D,n_i^D,n_inc^D,ineighbor,ipe_neighbor,ixS^L,ipole,idir

        ipole=neighbor_pole(i^D,igrid)

        {do ic^DB=1+int((1-i^DB)/2),2-int((1+i^DB)/2)
           inc^DB=2*i^DB+ic^DB\}
           ipe_neighbor=neighbor_child(2,inc^D,igrid)
           if(ipe_neighbor/=mype) then
             ineighbor=neighbor_child(1,inc^D,igrid)
             if(ipole==0) then
               n_i^D=-i^D;
               n_inc^D=ic^D+n_i^D;
               isend_c=isend_c+1
               itag=(3**^ND+4**^ND)*(ineighbor-1)+3**^ND+{n_inc^D*4**(^D-1)+}
               call MPI_ISEND(psb(igrid)%w,1,type_send_p(inc^D), &
                              ipe_neighbor,itag,icomm,sendrequest_c_p(isend_c),ierrmpi)
               if(stagger_grid) then 
                 ibuf_start=ibuf_send_p
                 do idir=1,ndim
                   ixS^L=ixS_p_stg_^L(idir,inc^D);
                   ibuf_next=ibuf_start+sizes_p_send_stg(idir,inc^D)
                   shapes=(/sizes_p_send_stg(idir,inc^D)/)
                   sendbuffer_p(ibuf_start:ibuf_next-1)=&
                     reshape(psb(igrid)%ws(ixS^S,idir),shapes)   
                   ibuf_start=ibuf_next
                 end do
                 isend_p=isend_p+1
                 call MPI_ISEND(sendbuffer_p(ibuf_send_p),sizes_p_send_total(inc^D),&
                                MPI_DOUBLE_PRECISION,ipe_neighbor,itag, &
                                icomm,sendrequest_p(isend_p),ierrmpi)
                 ibuf_send_p=ibuf_next
               end if
             else
               ixS^L=ixS_p_^L(inc^D);
               select case (ipole)
               {case (^D)
                 n_inc^D=inc^D^D%n_inc^DD=ic^DD-i^DD;\}
               end select
               if(isend_buf(ipwbuf)/=0) then
                 call MPI_WAIT(sendrequest_c_p(isend_buf(ipwbuf)), &
                               sendstatus_c_p(:,isend_buf(ipwbuf)),ierrmpi)
                 deallocate(pwbuf(ipwbuf)%w)
               end if
               allocate(pwbuf(ipwbuf)%w(ixS^S,nwhead:nwtail))
               call pole_buffer(pwbuf(ipwbuf)%w,ixS^L,ixS^L,psb(igrid)%w,ixG^LL,ixS^L,ipole)
               isend_c=isend_c+1
               isend_buf(ipwbuf)=isend_c
               itag=(3**^ND+4**^ND)*(ineighbor-1)+3**^ND+{n_inc^D*4**(^D-1)+}
               isizes={(ixSmax^D-ixSmin^D+1)*}*nwbc
               call MPI_ISEND(pwbuf(ipwbuf)%w,isizes,MPI_DOUBLE_PRECISION, &
                              ipe_neighbor,itag,icomm,sendrequest_c_p(isend_c),ierrmpi)
               ipwbuf=1+modulo(ipwbuf,npwbuf)
               if(stagger_grid) then 
                 ibuf_start=ibuf_send_p
                 do idir=1,ndim
                   ixS^L=ixS_p_stg_^L(idir,inc^D);
                   ibuf_next=ibuf_start+sizes_p_send_stg(idir,inc^D)
                   shapes=(/sizes_p_send_stg(idir,inc^D)/)
                   sendbuffer_p(ibuf_start:ibuf_next-1)=&
                     reshape(psb(igrid)%ws(ixS^S,idir),shapes)   
                   ibuf_start=ibuf_next
                 end do
                 isend_p=isend_p+1
                 call MPI_ISEND(sendbuffer_p(ibuf_send_p),sizes_p_send_total(inc^D),&
                                MPI_DOUBLE_PRECISION,ipe_neighbor,itag, &
                                icomm,sendrequest_p(isend_p),ierrmpi)
                 ibuf_send_p=ibuf_next
               end if
             end if
           end if
        {end do\}

      end subroutine bc_send_prolong

      !> Send to finer neighbor
      subroutine bc_fill_prolong(igrid,i^D)
        integer, intent(in) :: igrid,i^D

        integer :: ipe_neighbor,ineighbor,ixS^L,ixR^L,ic^D,inc^D,n_i^D,n_inc^D,ipole,idir

        ipole=neighbor_pole(i^D,igrid)

        if(ipole==0) then
          {do ic^DB=1+int((1-i^DB)/2),2-int((1+i^DB)/2)
             inc^DB=2*i^DB+ic^DB\}
             ipe_neighbor=neighbor_child(2,inc^D,igrid)
             if(ipe_neighbor==mype) then
               ixS^L=ixS_p_^L(inc^D);
               ineighbor=neighbor_child(1,inc^D,igrid)
               n_i^D=-i^D;
               n_inc^D=ic^D+n_i^D;
               ixR^L=ixR_p_^L(n_inc^D);
               psc(ineighbor)%w(ixR^S,nwhead:nwtail) &
                  =psb(igrid)%w(ixS^S,nwhead:nwtail)
               if(stagger_grid) then
                 do idir=1,ndim
                   ixS^L=ixS_p_stg_^L(idir,inc^D);
                   ixR^L=ixR_p_stg_^L(idir,n_inc^D);
                   psc(ineighbor)%ws(ixR^S,idir)=psb(igrid)%ws(ixS^S,idir)
                 end do
               end if
             end if
          {end do\}
        else
          {do ic^DB=1+int((1-i^DB)/2),2-int((1+i^DB)/2)
             inc^DB=2*i^DB+ic^DB\}
             ipe_neighbor=neighbor_child(2,inc^D,igrid)
             if(ipe_neighbor==mype) then
               ixS^L=ixS_p_^L(inc^D);
               ineighbor=neighbor_child(1,inc^D,igrid)
               select case (ipole)
               {case (^D)
                  n_inc^D=inc^D^D%n_inc^DD=ic^DD-i^DD;\}
               end select
               ixR^L=ixR_p_^L(n_inc^D);
               call pole_copy(psc(ineighbor)%w,ixCoG^L,ixR^L,psb(igrid)%w,ixG^LL,ixS^L,ipole)
               if(stagger_grid) then
                 do idir=1,ndim
                   ixS^L=ixS_p_stg_^L(idir,inc^D);
                   ixR^L=ixR_p_stg_^L(idir,n_inc^D);
                   call pole_copy_stg(psc(ineighbor)%ws,ixCoGs^L,ixR^L,psb(igrid)%ws,ixGs^LL,ixS^L,idir,ipole)
                 end do
               end if
             end if
          {end do\}
        end if
      end subroutine bc_fill_prolong

      subroutine gc_prolong(igrid)
        integer, intent(in) :: igrid

        integer :: i^D,idims,iside
        logical,dimension(-1:1^D&) :: NeedProlong

        NeedProlong=.false.
        {do i^DB=-1,1\}
           if (skip_direction([ i^D ])) cycle
           if (neighbor_type(i^D,igrid)==neighbor_coarse) then
             call bc_prolong(igrid,i^D)
             NeedProlong(i^D)=.true.
           end if
        {end do\}
        if(stagger_grid) then
          ! Ghost cell prolongation for staggered variables
          ! must be done in a specific order.
          ! First the first neighbours, which have 2 indices=0 in 3D
          ! or one index=0 in 2D
          block=>psb(igrid)
          do idims=1,ndim
            i^D=0;
            select case(idims)
           {case(^D)
              do i^D=-1,1,2
                if (NeedProlong(i^DD)) call bc_prolong_stg(igrid,i^DD,NeedProlong)
              end do
            \}
            end select
          end do
          ! Then the second neighbours which have 1 index=0 in 3D
          ! (Only in 3D)
          {^IFTHREED
          i1=0;
          do i2=-1,1,2
            do i3=-1,1,2
              if (NeedProlong(i^D)) call bc_prolong_stg(igrid,i^D,NeedProlong)
            end do
          end do
          i2=0;
          do i3=-1,1,2
            do i1=-1,1,2
              if (NeedProlong(i^D)) call bc_prolong_stg(igrid,i^D,NeedProlong)
            end do
          end do
          i3=0;
          do i1=-1,1,2
            do i2=-1,1,2
              if (NeedProlong(i^D)) call bc_prolong_stg(igrid,i^D,NeedProlong)
            end do
          end do
          }
          ! Finally, the corners, that have no index=0
         {do i^D=-1,1,2\}
            if (NeedProlong(i^D)) call bc_prolong_stg(igrid,i^D,NeedProlong)
         {end do\}
        end if
      end subroutine gc_prolong

      !> fill coarser representative with data from coarser neighbors
      subroutine bc_fill_prolong_stg(igrid,i^D)
        integer, intent(in) :: igrid,i^D

        integer :: ipe_neighbor,ineighbor,ipole,ixR^L,ic^D,inc^D,n_inc^D,idir

        ic^D=1+modulo(node(pig^D_,igrid)-1,2);
        if ({.not.(i^D==0.or.i^D==2*ic^D-3)|.or.}) return

        ipe_neighbor=neighbor(2,i^D,igrid)
        if(ipe_neighbor/=mype) then
          ineighbor=neighbor(1,i^D,igrid)
          ipole=neighbor_pole(i^D,igrid)

          if (ipole==0) then   !! There is no pole 
            inc^D=ic^D+i^D;
            ixR^L=ixR_p_^L(inc^D);
            do idir=1,ndim
              ixR^L=ixR_p_stg_^L(idir,inc^D);
              ibuf_next=ibuf_recv_p+sizes_p_recv_stg(idir,inc^D)
              psc(igrid)%ws(ixR^S,idir)=reshape(source=recvbuffer_p(ibuf_recv_p:ibuf_next-1),&
                    shape=shape(psc(igrid)%ws(ixR^S,idir)))
              ibuf_recv_p=ibuf_next
            end do
          else !! There is a pole
            inc^D=ic^D+i^D;
            select case (ipole)
            {case (^D)
               n_inc^D=2*i^D+(3-ic^D)^D%n_inc^DD=-2*i^DD+ic^DD;\}
            end select
            !! Unpack the buffer and fill an auxiliary array
            pole_buf%ws=zero
            do idir=1,ndim
              ixR^L=ixR_p_stg_^L(idir,inc^D);
              ibuf_next=ibuf_recv_p+sizes_p_recv_stg(idir,inc^D)
              pole_buf%ws(ixR^S,idir)=reshape(source=recvbuffer_p(ibuf_recv_p:ibuf_next-1),&
                shape=shape(psc(igrid)%ws(ixR^S,idir)))
              call pole_copy_stg(psc(igrid)%ws,ixCoGs^L,ixR^L,pole_buf%ws,ixGs^LL,ixR^L,idir,ipole)
              ibuf_recv_p=ibuf_next
            end do
          end if
        end if

      end subroutine bc_fill_prolong_stg

      !> do prolongation for fine blocks after receipt data from coarse neighbors
      subroutine bc_prolong(igrid,i^D)
        use mod_physics, only: phys_to_primitive, phys_to_conserved

        double precision :: dxFi^D, dxCo^D, xFimin^D, xComin^D, invdxCo^D
        integer :: i^D,igrid
        integer :: ixFi^L,ixCo^L,ii^D, idims,iside,ixB^L

        ixFi^L=ixR_srl_^L(i^D);
        dxFi^D=rnode(rpdx^D_,igrid);
        dxCo^D=two*dxFi^D;
        invdxCo^D=1.d0/dxCo^D;

        ! compute the enlarged grid lower left corner coordinates
        ! these are true coordinates for an equidistant grid, 
        ! but we can temporarily also use them for getting indices 
        ! in stretched grids
        xFimin^D=rnode(rpxmin^D_,igrid)-dble(nghostcells)*dxFi^D;
        xComin^D=rnode(rpxmin^D_,igrid)-dble(nghostcells)*dxCo^D;

        if(phyboundblock(igrid).and.bcphys) then
          block=>psc(igrid)
          do idims=1,ndim
            ixComin^D=int((xFimin^D+(dble(ixFimin^D)-half)*dxFi^D-xComin^D)*invdxCo^D)+1-1;
            ixComax^D=int((xFimin^D+(dble(ixFimax^D)-half)*dxFi^D-xComin^D)*invdxCo^D)+1+1;
            {^IFTHREED
            ! avoid using undetermined ghost cells at physical boundary edges
            if(idims == 1) then
              if(neighbor_type(-1,0,0,igrid)==neighbor_boundary .or. &
                 neighbor_type(1,0,0,igrid)==neighbor_boundary) then
                if(neighbor_type(0,-1,0,igrid)==neighbor_boundary) ixComin2=ixCoMmin2
                if(neighbor_type(0,0,-1,igrid)==neighbor_boundary) ixComin3=ixCoMmin3
                if(neighbor_type(0,1,0,igrid)==neighbor_boundary) ixComax2=ixCoMmax2
                if(neighbor_type(0,0,1,igrid)==neighbor_boundary) ixComax3=ixCoMmax3
              end if
            else if(idims == 2) then
              if(neighbor_type(0,-1,0,igrid)==neighbor_boundary .or. &
                 neighbor_type(0,1,0,igrid)==neighbor_boundary) then
                if(neighbor_type(0,0,-1,igrid)==neighbor_boundary) ixComin3=ixCoMmin3
                if(neighbor_type(0,0,1,igrid)==neighbor_boundary) ixComax3=ixCoMmax3
              end if
            end if
            }
            do iside=1,2
              ii^D=kr(^D,idims)*(2*iside-3);
              if(neighbor_type(ii^D,igrid)/=neighbor_boundary) cycle
              if(( {(iside==1.and.idims==^D.and.ixComin^D<ixCoGmin^D+nghostcells)|.or.} ) &
               .or.( {(iside==2.and.idims==^D.and.ixComax^D>ixCoGmax^D-nghostcells)|.or. })) then
                {ixBmin^D=merge(ixCoGmin^D,ixComin^D,idims==^D);}
                {ixBmax^D=merge(ixCoGmax^D,ixComax^D,idims==^D);}
                call bc_phys(iside,idims,time,0.d0,psc(igrid),ixCoG^L,ixB^L)
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
          ixComin^D=int((xFimin^D+(dble(ixFimin^D)-half)*dxFi^D-xComin^D)*invdxCo^D)+1-1;
          ixComax^D=int((xFimin^D+(dble(ixFimax^D)-half)*dxFi^D-xComin^D)*invdxCo^D)+1+1;
          call phys_to_primitive(ixCoG^L,ixCo^L,psc(igrid)%w,psc(igrid)%x)
        end if

        if(ghost_copy) then
          call interpolation_copy(igrid,ixFi^L,dxFi^D,xFimin^D,dxCo^D,invdxCo^D,xComin^D)
        else
          call interpolation_linear(igrid,ixFi^L,dxFi^D,xFimin^D,dxCo^D,invdxCo^D,xComin^D)
        end if

        if(prolongprimitive) then
          block=>psc(igrid)
          call phys_to_conserved(ixCoG^L,ixCo^L,psc(igrid)%w,psc(igrid)%x)
        end if

      end subroutine bc_prolong

      subroutine bc_prolong_stg(igrid,i^D,NeedProlong)
        use mod_amr_fct
        double precision           :: dxFi^D,dxCo^D,xFimin^D,xComin^D,invdxCo^D
        integer                    :: igrid,i^D
        integer                    :: ixFi^L,ixCo^L
        logical,dimension(-1:1^D&) :: NeedProlong
        logical                    :: fine_^Lin
        ! Check what is already at the desired level
        fine_^Lin=.false.;
        {
        if(i^D>-1) fine_min^Din=(.not.NeedProlong(i^DD-kr(^D,^DD)).and.neighbor_type(i^DD-kr(^D,^DD),igrid)/=1)
        if(i^D<1)  fine_max^Din=(.not.NeedProlong(i^DD+kr(^D,^DD)).and.neighbor_type(i^DD+kr(^D,^DD),igrid)/=1)
        \}

        ixFi^L=ixR_srl_^L(i^D);

        dxFi^D=rnode(rpdx^D_,igrid);
        dxCo^D=two*dxFi^D;
        invdxCo^D=1.d0/dxCo^D;

        xFimin^D=rnode(rpxmin^D_,igrid)-dble(nghostcells)*dxFi^D;
        xComin^D=rnode(rpxmin^D_,igrid)-dble(nghostcells)*dxCo^D;

        ! moved the physical boundary filling here, to only fill the
        ! part needed

        ixComin^D=int((xFimin^D+(dble(ixFimin^D)-half)*dxFi^D-xComin^D)*invdxCo^D)+1-1;
        ixComax^D=int((xFimin^D+(dble(ixFimax^D)-half)*dxFi^D-xComin^D)*invdxCo^D)+1+1;

        if(prolongprimitive) call phys_to_primitive(ixG^LL,ixFi^L,psb(igrid)%w,psb(igrid)%x)

        call prolong_2nd_stg(psc(igrid),psb(igrid),ixCo^L,ixFi^L,dxCo^D,xComin^D,dxFi^D,xFimin^D,.true.,fine_^Lin)

        if(prolongprimitive) call phys_to_conserved(ixG^LL,ixFi^L,psb(igrid)%w,psb(igrid)%x)

        ! The current region has already been refined, so it does not need to be prolonged again
        NeedProlong(i^D)=.false. 

      end subroutine bc_prolong_stg

      subroutine interpolation_linear(igrid,ixFi^L,dxFi^D,xFimin^D, &
                                      dxCo^D,invdxCo^D,xComin^D)
        use mod_physics, only: phys_to_conserved
        integer, intent(in) :: igrid, ixFi^L
        double precision, intent(in) :: dxFi^D, xFimin^D,dxCo^D, invdxCo^D, xComin^D

        double precision :: xCo^D, xFi^D, eta^D
        double precision :: slopeL, slopeR, slopeC, signC, signR
        double precision :: slope(1:nw,ndim)
        !!double precision :: local_invdxCo^D
        double precision :: signedfactorhalf^D
        integer :: ixCo^D, jxCo^D, hxCo^D, ixFi^D, ix^D, iw, idims, nwmin,nwmax
        !integer :: ixshift^D, icase

        !icase=mod(nghostcells,2)

        if(prolongprimitive) then
          nwmin=1
          nwmax=nw
        else
          nwmin=nwhead
          nwmax=nwtail
        end if

        {do ixFi^DB = ixFi^LIM^DB
           ! cell-centered coordinates of fine grid point
           ! here we temporarily use an equidistant grid
           xFi^DB=xFimin^DB+(dble(ixFi^DB)-half)*dxFi^DB

           ! indices of coarse cell which contains the fine cell
           ! since we computed lower left corner earlier 
           ! in equidistant fashion: also ok for stretched case
           ixCo^DB=int((xFi^DB-xComin^DB)*invdxCo^DB)+1

           ! cell-centered coordinates of coarse grid point
           ! here we temporarily use an equidistant grid
           xCo^DB=xComin^DB+(dble(ixCo^DB)-half)*dxCo^DB \}

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
             eta^D=(xFi^D-xCo^D)*invdxCo^D;
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
             ix^D=2*int((ixFi^D+ixMlo^D)/2)-ixMlo^D;
             {if(xFi^D>xCo^D) then
                signedfactorhalf^D=0.5d0
              else
                signedfactorhalf^D=-0.5d0
              end if
              eta^D=signedfactorhalf^D*(one-psb(igrid)%dvolume(ixFi^DD) &
                   /sum(psb(igrid)%dvolume(ix^D:ix^D+1^D%ixFi^DD))) \}
             !{eta^D=(xFi^D-xCo^D)*invdxCo^D &
             !      *two*(one-block%dvolume(ixFi^DD) &
             !      /sum(block%dvolume(ix^D:ix^D+1^D%ixFi^DD))) \}
           end if

           do idims=1,ndim
              hxCo^D=ixCo^D-kr(^D,idims)\
              jxCo^D=ixCo^D+kr(^D,idims)\

              do iw=nwmin,nwmax
                 slopeL=psc(igrid)%w(ixCo^D,iw)-psc(igrid)%w(hxCo^D,iw)
                 slopeR=psc(igrid)%w(jxCo^D,iw)-psc(igrid)%w(ixCo^D,iw)
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
                   slope(iw,idims)=signC*max(zero,min(dabs(slopeC), &
                                     signC*slopeL,signC*slopeR))
                 !end select
              end do
           end do

           ! Interpolate from coarse cell using limited slopes
           psb(igrid)%w(ixFi^D,nwmin:nwmax)=psc(igrid)%w(ixCo^D,nwmin:nwmax)+&
             {(slope(nwmin:nwmax,^D)*eta^D)+}

        {end do\}

        if(prolongprimitive) then
          block=>psb(igrid)
          call phys_to_conserved(ixG^LL,ixFi^L,psb(igrid)%w,psb(igrid)%x)
        end if

      end subroutine interpolation_linear

      subroutine interpolation_copy(igrid, ixFi^L,dxFi^D,xFimin^D, &
                                    dxCo^D,invdxCo^D,xComin^D)
        use mod_physics, only: phys_to_conserved
        integer, intent(in) :: igrid, ixFi^L
        double precision, intent(in) :: dxFi^D, xFimin^D,dxCo^D, invdxCo^D, xComin^D

        double precision :: xFi^D
        integer :: ixCo^D, ixFi^D, nwmin,nwmax

        if(prolongprimitive) then
          nwmin=1
          nwmax=nw
        else
          nwmin=nwhead
          nwmax=nwtail
        end if

        {do ixFi^DB = ixFi^LIM^DB
           ! cell-centered coordinates of fine grid point
           xFi^DB=xFimin^DB+(dble(ixFi^DB)-half)*dxFi^DB
        
           ! indices of coarse cell which contains the fine cell
           ! note: this also works for stretched grids
           ixCo^DB=int((xFi^DB-xComin^DB)*invdxCo^DB)+1\}
        
           ! Copy from coarse cell
           psb(igrid)%w(ixFi^D,nwmin:nwmax)=psc(igrid)%w(ixCo^D,nwmin:nwmax)
        
        {end do\}
        
        if(prolongprimitive) call phys_to_conserved(ixG^LL,ixFi^L,psb(igrid)%w,psb(igrid)%x)
      
      end subroutine interpolation_copy

      subroutine pole_copy(wrecv,ixIR^L,ixR^L,wsend,ixIS^L,ixS^L,ipole)
      
        integer, intent(in) :: ixIR^L,ixR^L,ixIS^L,ixS^L,ipole
        double precision :: wrecv(ixIR^S,1:nw), wsend(ixIS^S,1:nw)

        integer :: iw, iside, iB

        select case (ipole)
        {case (^D)
           iside=int((i^D+3)/2)
           iB=2*(^D-1)+iside
           do iw=nwhead,nwtail
             select case (typeboundary(iw,iB))
             case (bc_symm)
               wrecv(ixR^S,iw) = wsend(ixSmax^D:ixSmin^D:-1^D%ixS^S,iw)
             case (bc_asymm)
               wrecv(ixR^S,iw) =-wsend(ixSmax^D:ixSmin^D:-1^D%ixS^S,iw)
             case default
               call mpistop("Pole boundary condition should be symm or asymm")
             end select
           end do \}
        end select
      
      end subroutine pole_copy

      subroutine pole_copy_stg(wrecv,ixIR^L,ixR^L,wsend,ixIS^L,ixS^L,idirs,ipole)
      
        integer, intent(in) :: ixIR^L,ixR^L,ixIS^L,ixS^L,idirs,ipole

        double precision :: wrecv(ixIR^S,1:nws), wsend(ixIS^S,1:nws)
        integer :: iB, iside

        select case (ipole)
        {case (^D)
           iside=int((i^D+3)/2)
           iB=2*(^D-1)+iside
           select case (typeboundary(iw_mag(idirs),iB))
           case (bc_symm)
             wrecv(ixR^S,idirs) = wsend(ixSmax^D:ixSmin^D:-1^D%ixS^S,idirs)
           case (bc_asymm)
             wrecv(ixR^S,idirs) =-wsend(ixSmax^D:ixSmin^D:-1^D%ixS^S,idirs)
           case default
             call mpistop("Pole boundary condition should be symm or asymm")
           end select
         \}
        end select

      end subroutine pole_copy_stg

      subroutine pole_buffer(wrecv,ixIR^L,ixR^L,wsend,ixIS^L,ixS^L,ipole)
      
        integer, intent(in) :: ixIR^L,ixR^L,ixIS^L,ixS^L,ipole
        double precision :: wrecv(ixIR^S,nwhead:nwtail), wsend(ixIS^S,1:nw)

        integer :: iw, iside, iB

        select case (ipole)
        {case (^D)
           iside=int((i^D+3)/2)
           iB=2*(^D-1)+iside
           do iw=nwhead,nwtail
             select case (typeboundary(iw,iB))
             case (bc_symm)
               wrecv(ixR^S,iw) = wsend(ixSmax^D:ixSmin^D:-1^D%ixS^S,iw)
             case (bc_asymm)
               wrecv(ixR^S,iw) =-wsend(ixSmax^D:ixSmin^D:-1^D%ixS^S,iw)
             case default
               call mpistop("Pole boundary condition should be symm or asymm")
             end select
           end do \}
        end select
      
      end subroutine pole_buffer

  end subroutine getbc

end module mod_ghostcells_update
