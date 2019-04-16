!> update ghost cells of all blocks including physical boundaries 
module mod_ghostcells_update

  implicit none

  ! A switch of update physical boundary or not
  logical, public :: bcphys=.true.

  integer :: ixM^L, ixCoG^L, ixCoM^L

  ! The first index goes from -1:2, where -1 is used when a block touches the
  ! lower boundary, 1 when a block touches an upper boundary, and 0 a situation
  ! away from boundary conditions, 2 when a block touched both lower and upper
  ! boundary

  ! index ranges to send (S) to sibling blocks, receive (R) from sibling blocks
  integer, dimension(-1:2,-1:1) :: ixS_srl_^L, ixR_srl_^L

  ! index ranges to send (S) restricted (r) ghost cells to coarser blocks 
  integer, dimension(-1:1,-1:1) :: ixS_r_^L

  ! index ranges to receive restriced ghost cells from finer blocks 
  integer, dimension(-1:1, 0:3) :: ixR_r_^L

  ! send prolongated (p) ghost cells to finer blocks, receive prolongated 
  ! ghost from coarser blocks
  integer, dimension(-1:1, 0:3) :: ixS_p_^L, ixR_p_^L

  integer, dimension(^ND,-1:1) :: ixS_srl_stg_^L, ixR_srl_stg_^L
  integer, dimension(^ND,-1:1) :: ixS_r_stg_^L
  integer, dimension(^ND,0:3)  :: ixS_p_stg_^L, ixR_p_stg_^L
  integer, dimension(^ND,0:3)  :: ixR_r_stg_^L

  ! number of MPI receive-send pairs, srl: same refinement level; r: restrict; p: prolong
  integer :: nrecv_bc_srl, nsend_bc_srl, nrecv_bc_r, nsend_bc_r, nrecv_bc_p, nsend_bc_p

  ! total size of buffer arrays
  integer :: nbuff_bc_recv_srl, nbuff_bc_send_srl, nbuff_bc_recv_r, nbuff_bc_send_r, nbuff_bc_recv_p, nbuff_bc_send_p

  ! record index position of buffer arrays
  integer :: ibuf_send_srl, ibuf_recv_srl, ibuf_send_r, ibuf_recv_r, ibuf_send_p, ibuf_recv_p

  ! count of times of send and receive
  integer :: isend_srl, irecv_srl, isend_r, irecv_r, isend_p, irecv_p

  ! sizes to send and receive for siblings
  integer, dimension(-1:1^D&) :: sizes_srl_send, sizes_srl_recv

  ! total sizes = cell-center normal flux + stagger-grid flux of send and receive
  integer, dimension(-1:1^D&) :: sizes_srl_send_total, sizes_srl_recv_total

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
  integer, dimension(-1:1^D&)     :: sizes_r_send, sizes_r_send_total
  integer, dimension(0:3^D&)      :: sizes_r_recv, sizes_r_recv_total
  integer, dimension(^ND,-1:1^D&) :: sizes_r_send_stg
  integer, dimension(^ND,0:3^D&)  :: sizes_r_recv_stg

  ! sizes to allocate buffer arrays for send and receive for restriction
  integer, dimension(0:3^D&)      :: sizes_p_send, sizes_p_recv
  integer, dimension(0:3^D&)      :: sizes_p_send_total, sizes_p_recv_total
  integer, dimension(^ND,0:3^D&)  :: sizes_p_send_stg, sizes_p_recv_stg

  ! index pointer for buffer arrays as a start for a segment
  integer :: ibuf_start, ibuf_next

  ! shapes of reshape
  integer, dimension(1) :: shapes

contains

  subroutine init_bc()
    use mod_global_parameters 
    use mod_physics, only: phys_req_diagonal

    integer :: nghostcellsCo, interpolation_order
    integer :: nx^D, nxCo^D, ixG^L, i^D, ic^D, inc^D, idir

    ixG^L=ixG^LL;
    ixM^L=ixG^L^LSUBnghostcells;
    ixCoGmin^D=1;
    !ixCoGmax^D=ixGmax^D/2+nghostcells;
    ixCoGmax^D=(ixGhi^D-2*nghostcells)/2+2*nghostcells;

    ixCoM^L=ixCoG^L^LSUBnghostcells;

    nx^D=ixMmax^D-ixMmin^D+1;
    nxCo^D=nx^D/2;

    select case (typeghostfill)
    case ("copy")
       interpolation_order=1
    case ("linear")
       interpolation_order=2
    case default
       write (unitterm,*) "Undefined typeghostfill ",typeghostfill
       call mpistop("Undefined typeghostfill")
    end select
    nghostcellsCo=int((nghostcells+1)/2)

    if (nghostcellsCo+interpolation_order-1>nghostcells) then
       call mpistop("interpolation order for prolongation in getbc too high")
    end if
    
    allocate(pole_buf%w(ixG^T,1:nwflux+nwaux))
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
    end if

    ! (iib,i) index has following meanings: iib = 0 means it is not at any physical boundary
    ! iib=-1 means it is at the minimum side of a physical boundary  
    ! iib= 1 means it is at the maximum side of a physical boundary  
    ! i=-1 means subregion prepared for the neighbor at its minimum side 
    ! i= 1 means subregion prepared for the neighbor at its maximum side 

    ! index limits for siblings
    {
    ixS_srl_min^D(:,-1)=ixMmin^D
    ixS_srl_min^D(:, 1)=ixMmax^D+1-nghostcells
    ixS_srl_max^D(:,-1)=ixMmin^D-1+nghostcells
    ixS_srl_max^D(:, 1)=ixMmax^D
    
    ixS_srl_min^D(-1,0)=1
    ixS_srl_min^D( 0,0)=ixMmin^D
    ixS_srl_min^D( 1,0)=ixMmin^D
    ixS_srl_min^D( 2,0)=1
    ixS_srl_max^D(-1,0)=ixMmax^D
    ixS_srl_max^D( 0,0)=ixMmax^D
    ixS_srl_max^D( 1,0)=ixGmax^D
    ixS_srl_max^D( 2,0)=ixGmax^D
     
    ixR_srl_min^D(:,-1)=1
    ixR_srl_min^D(:, 1)=ixMmax^D+1
    ixR_srl_max^D(:,-1)=nghostcells
    ixR_srl_max^D(:, 1)=ixGmax^D
    
    ixR_srl_min^D(-1,0)=1
    ixR_srl_min^D( 0,0)=ixMmin^D
    ixR_srl_min^D( 1,0)=ixMmin^D
    ixR_srl_min^D( 2,0)=1
    ixR_srl_max^D(-1,0)=ixMmax^D
    ixR_srl_max^D( 0,0)=ixMmax^D
    ixR_srl_max^D( 1,0)=ixGmax^D
    ixR_srl_max^D( 2,0)=ixGmax^D
    \}

    ! index limits for restrict
    { 
    ixS_r_min^D(:,-1)=ixCoMmin^D
    ixS_r_min^D(:, 1)=ixCoMmax^D+1-nghostcells
    ixS_r_max^D(:,-1)=ixCoMmin^D-1+nghostcells
    ixS_r_max^D(:, 1)=ixCoMmax^D

    ixS_r_min^D(-1,0)=1
    ixS_r_min^D( 0,0)=ixCoMmin^D
    ixS_r_min^D( 1,0)=ixCoMmin^D
    ixS_r_max^D(-1,0)=ixCoMmax^D
    ixS_r_max^D( 0,0)=ixCoMmax^D
    ixS_r_max^D( 1,0)=ixCoGmax^D

    ixR_r_min^D(:, 0)=1
    ixR_r_min^D(:, 1)=ixMmin^D
    ixR_r_min^D(:, 2)=ixMmin^D+nxCo^D
    ixR_r_min^D(:, 3)=ixMmax^D+1
    ixR_r_max^D(:, 0)=nghostcells
    ixR_r_max^D(:, 1)=ixMmin^D-1+nxCo^D
    ixR_r_max^D(:, 2)=ixMmax^D
    ixR_r_max^D(:, 3)=ixGmax^D

    ixR_r_min^D(-1,1)=1
    ixR_r_max^D(-1,1)=ixMmin^D-1+nxCo^D
    ixR_r_min^D( 1,2)=ixMmin^D+nxCo^D
    ixR_r_max^D( 1,2)=ixGmax^D
    \}

    ! index limits for prolong
    {
    ixS_p_min^D(:, 0)=ixMmin^D-(interpolation_order-1)
    ixS_p_min^D(:, 1)=ixMmin^D-(interpolation_order-1)
    ixS_p_min^D(:, 2)=ixMmin^D+nxCo^D-nghostcellsCo-(interpolation_order-1)
    ixS_p_min^D(:, 3)=ixMmax^D+1-nghostcellsCo-(interpolation_order-1)
    ixS_p_max^D(:, 0)=ixMmin^D-1+nghostcellsCo+(interpolation_order-1)
    ixS_p_max^D(:, 1)=ixMmin^D-1+nxCo^D+nghostcellsCo+(interpolation_order-1)
    ixS_p_max^D(:, 2)=ixMmax^D+(interpolation_order-1)
    ixS_p_max^D(:, 3)=ixMmax^D+(interpolation_order-1)

    if(.not.phys_req_diagonal) then
      ! exclude ghost-cell region when diagonal cells are unknown
      ixS_p_min^D(:, 0)=ixMmin^D
      ixS_p_max^D(:, 3)=ixMmax^D
      ixS_p_max^D(:, 1)=ixMmin^D-1+nxCo^D+(interpolation_order-1)
      ixS_p_min^D(:, 2)=ixMmin^D+nxCo^D-(interpolation_order-1)
    end if

    ! extend index range to physical boundary
    ixS_p_min^D(-1,1)=1
    ixS_p_max^D( 1,2)=ixGmax^D

    ixR_p_min^D(:, 0)=ixCoMmin^D-nghostcellsCo-(interpolation_order-1)
    ixR_p_min^D(:, 1)=ixCoMmin^D-(interpolation_order-1)
    ixR_p_min^D(:, 2)=ixCoMmin^D-nghostcellsCo-(interpolation_order-1)
    ixR_p_min^D(:, 3)=ixCoMmax^D+1-(interpolation_order-1)
    ixR_p_max^D(:, 0)=nghostcells+(interpolation_order-1)
    ixR_p_max^D(:, 1)=ixCoMmax^D+nghostcellsCo+(interpolation_order-1)
    ixR_p_max^D(:, 2)=ixCoMmax^D+(interpolation_order-1)
    ixR_p_max^D(:, 3)=ixCoMmax^D+nghostcellsCo+(interpolation_order-1)

    if(.not.phys_req_diagonal) then
      ! exclude ghost-cell region when diagonal cells are unknown
      ixR_p_max^D(:, 0)=nghostcells
      ixR_p_min^D(:, 3)=ixCoMmax^D+1
      ixR_p_max^D(:, 1)=ixCoMmax^D+(interpolation_order-1)
      ixR_p_min^D(:, 2)=ixCoMmin^D-(interpolation_order-1)
    end if

    ! extend index range to physical boundary
    ixR_p_min^D(-1,1)=1
    ixR_p_max^D( 1,2)=ixCoGmax^D
    \}

    ! calculate sizes for buffer arrays for siblings
    {do i^DB=-1,1\}
      sizes_srl_send(i^D)=(nwflux+nwaux)*{(ixS_srl_max^D(2,i^D)-ixS_srl_min^D(2,i^D)+1)|*}
      sizes_srl_recv(i^D)=(nwflux+nwaux)*{(ixR_srl_max^D(2,i^D)-ixR_srl_min^D(2,i^D)+1)|*}
      sizes_srl_send_total(i^D)=sizes_srl_send(i^D)
      sizes_srl_recv_total(i^D)=sizes_srl_recv(i^D)
      if(stagger_grid) then
        ! Staggered (face-allocated) variables
        do idir=1,ndim
          sizes_srl_send_stg(idir,i^D)={(ixS_srl_stg_max^D(idir,i^D)-ixS_srl_stg_min^D(idir,i^D)+1)|*}
          sizes_srl_recv_stg(idir,i^D)={(ixR_srl_stg_max^D(idir,i^D)-ixR_srl_stg_min^D(idir,i^D)+1)|*}
          sizes_srl_send_total(i^D)=sizes_srl_send_total(i^D)+sizes_srl_send_stg(idir,i^D)
          sizes_srl_recv_total(i^D)=sizes_srl_recv_total(i^D)+sizes_srl_recv_stg(idir,i^D)
        end do
      end if
    {end do\}

    ! Sizes for multi-resolution communications
    {do i^DB=-1,1\}
       ! Cell-centred variables
       sizes_r_send(i^D)=(nwflux+nwaux)*{(ixS_r_max^D(1,i^D)-ixS_r_min^D(1,i^D)+1)|*}
       sizes_r_send_total(i^D)=sizes_r_send(i^D)
       if(stagger_grid) then
         ! Staggered (face-allocated) variables
         do idir=1,ndim
           sizes_r_send_stg(idir,i^D)={(ixS_r_stg_max^D(idir,i^D)-ixS_r_stg_min^D(idir,i^D)+1)|*}
           sizes_r_send_total(i^D)=sizes_r_send_total(i^D)+sizes_r_send_stg(idir,i^D)
         end do
       end if
    {end do\}

    {do i^DB=0,3\}
       ! Cell-centred variables
       sizes_r_recv(i^D)=(nwflux+nwaux)*{(ixR_r_max^D(1,i^D)-ixR_r_min^D(-1,i^D)+1)|*}
       sizes_p_send(i^D)=(nwflux+nwaux)*{(ixS_p_max^D(1,i^D)-ixS_p_min^D(-1,i^D)+1)|*}
       sizes_p_recv(i^D)=(nwflux+nwaux)*{(ixR_p_max^D(1,i^D)-ixR_p_min^D(-1,i^D)+1)|*}

       sizes_r_recv_total(i^D)=sizes_r_recv(i^D)
       sizes_p_send_total(i^D)=sizes_p_send(i^D)
       sizes_p_recv_total(i^D)=sizes_p_recv(i^D)
    
       if (stagger_grid) then
       ! Staggered (face-allocated) variables
         do idir=1,^ND
           sizes_r_recv_stg(idir,i^D)={(ixR_r_stg_max^D(idir,i^D)-ixR_r_stg_min^D(idir,i^D)+1)|*}
           sizes_p_send_stg(idir,i^D)={(ixS_p_stg_max^D(idir,i^D)-ixS_p_stg_min^D(idir,i^D)+1)|*}
           sizes_p_recv_stg(idir,i^D)={(ixR_p_stg_max^D(idir,i^D)-ixR_p_stg_min^D(idir,i^D)+1)|*}
           sizes_r_recv_total(i^D)=sizes_r_recv_total(i^D)+sizes_r_recv_stg(idir,i^D)
           sizes_p_send_total(i^D)=sizes_p_send_total(i^D)+sizes_p_send_stg(idir,i^D)
           sizes_p_recv_total(i^D)=sizes_p_recv_total(i^D)+sizes_p_recv_stg(idir,i^D)
         end do
       end if

    {end do\}

  end subroutine init_bc

  subroutine getbc(time,qdt,psb,nwstart,nwbc,req_diag)
    use mod_global_parameters
    use mod_physics

    double precision, intent(in)      :: time, qdt
    type(state), target               :: psb(max_blocks)
    integer, intent(in)               :: nwstart ! Fill from nwstart+1
    integer, intent(in)               :: nwbc    ! Number of variables to fill
    logical, intent(in), optional     :: req_diag ! If false, skip diagonal ghost cells

    double precision :: time_bcin
    integer :: my_neighbor_type, ipole, idims, iside, nwhead, nwtail
    integer :: iigrid, igrid, ineighbor, ipe_neighbor
    integer :: nrecvs, nsends, isizes
    integer :: ixG^L, ixR^L, ixS^L, ixB^L, ixI^L, k^L
    integer :: i^D, n_i^D, ic^D, inc^D, n_inc^D, iib^D, isb^D, idir
    ! store physical boundary indicating index
    integer :: idphyb(ndim,max_blocks),bindex(ndim)
    integer :: nghostcellsco,iB
    logical  :: req_diagonal, NeedProlong(-1:1^D&)

    nwhead=nwstart+1
    nwtail=nwstart+nwbc

    req_diagonal = .true.
    if (present(req_diag)) req_diagonal = req_diag

    time_bcin=MPI_WTIME()
    ixG^L=ixG^LL;
    
    if (internalboundary) then 
       call getintbc(time,ixG^L)
    end if
    ! fill ghost cells in physical boundaries
    if(bcphys) then
      do iigrid=1,igridstail; igrid=igrids(iigrid);
         if(.not.phyboundblock(igrid)) cycle
         saveigrid=igrid
         block=>ps(igrid)
         ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
         do idims=1,ndim
            ! to avoid using as yet unknown corner info in more than 1D, we
            ! fill only interior mesh ranges of the ghost cell ranges at first,
            ! and progressively enlarge the ranges to include corners later
            {
             kmin^D=merge(0, 1, idims==^D)
             kmax^D=merge(0, 1, idims==^D)
             ixBmin^D=ixGmin^D+kmin^D*nghostcells
             ixBmax^D=ixGmax^D-kmax^D*nghostcells
            \}
            {^IFTWOD
             if(idims > 1 .and. neighbor_type(-1,0,igrid)==neighbor_boundary) ixBmin1=ixGmin1
             if(idims > 1 .and. neighbor_type( 1,0,igrid)==neighbor_boundary) ixBmax1=ixGmax1}
            {^IFTHREED
             if(idims > 1 .and. neighbor_type(-1,0,0,igrid)==neighbor_boundary) ixBmin1=ixGmin1
             if(idims > 1 .and. neighbor_type( 1,0,0,igrid)==neighbor_boundary) ixBmax1=ixGmax1
             if(idims > 2 .and. neighbor_type(0,-1,0,igrid)==neighbor_boundary) ixBmin2=ixGmin2
             if(idims > 2 .and. neighbor_type(0, 1,0,igrid)==neighbor_boundary) ixBmax2=ixGmax2}
            do iside=1,2
               i^D=kr(^D,idims)*(2*iside-3);
               if (aperiodB(idims)) then 
                  if (neighbor_type(i^D,igrid) /= neighbor_boundary .and. &
                       .not. psb(igrid)%is_physical_boundary(2*idims-2+iside)) cycle
               else 
                  if (neighbor_type(i^D,igrid) /= neighbor_boundary) cycle
               end if
               call bc_phys(iside,idims,time,qdt,psb(igrid)%w,ps(igrid)%x,ixG^L,ixB^L)
            end do
         end do
      end do
    end if

    ! receiving ghost-cell values from sibling blocks at another cpu
    irecv_srl=0
    ibuf_recv_srl=1
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       call identifyphysbound(ps(igrid),iib^D)
       ^D&idphyb(^D,igrid)=iib^D;
       {do i^DB=-1,1\}
          if (skip_direction([ i^D ])) cycle
          if(neighbor_type(i^D,igrid)==neighbor_sibling) call bc_recv_srl
       {end do\}
    end do

    ! sending ghost-cell values to sibling blocks at another cpu
    isend_srl=0
    ibuf_send_srl=1
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       ! Used stored data to identify physical boundaries
       ^D&iib^D=idphyb(^D,igrid);
       {do i^DB=-1,1\}
          if (skip_direction([ i^D ])) cycle
          if(neighbor_type(i^D,igrid)==neighbor_sibling) call bc_send_srl
       {end do\}
    end do

    if (nrecv_bc_srl>0) then
      call MPI_WAITALL(nrecv_bc_srl,recvrequest_srl,recvstatus_srl,ierrmpi)
    end if
    if (nsend_bc_srl>0) then
      call MPI_WAITALL(nsend_bc_srl,sendrequest_srl,sendstatus_srl,ierrmpi)
    end if

    ! unpack the received data to fill ghost cells
    ibuf_recv_srl=1
    do iigrid=1,igridstail; igrid=igrids(iigrid);
      ^D&iib^D=idphyb(^D,igrid);
      {do i^DB=-1,1\}
         if (skip_direction([ i^D ])) cycle
         if(neighbor_type(i^D,igrid)==neighbor_sibling) call bc_fill_srl
      {end do\}
    end do

    ! receiving ghost-cell values from finer neighbors
    irecv_r=0
    ibuf_recv_r=1
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       ^D&iib^D=idphyb(^D,igrid);
       {do i^DB=-1,1\}
          if (skip_direction([ i^D ])) cycle
          if(neighbor_type(i^D,igrid)==neighbor_fine) call bc_recv_restrict
       {end do\}
    end do

    ! sending ghost-cell values coarser neighbors
    isend_r=0
    ibuf_send_r=1
    nghostcellsco=ceiling(nghostcells*0.5d0)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       ! Used stored data to identify physical boundaries
       ^D&iib^D=idphyb(^D,igrid);
       if (any(neighbor_type(:^D&,igrid)==neighbor_coarse)) then
          call coarsen_grid(psb(igrid),ixG^L,ixM^L,psc(igrid),ixCoG^L,ixCoM^L)
         {do i^DB=-1,1\}
            if (skip_direction([ i^D ])) cycle
            if(neighbor_type(i^D,igrid)==neighbor_coarse) call bc_send_restrict
         {end do\}
       end if
    end do

    if (nrecv_bc_r>0) then
      call MPI_WAITALL(nrecv_bc_r,recvrequest_r,recvstatus_r,ierrmpi)
    end if
    if (nsend_bc_r>0) then
      call MPI_WAITALL(nsend_bc_r,sendrequest_r,sendstatus_r,ierrmpi)
    end if

    ! unpack the received data to fill ghost cells
    ibuf_recv_r=1
    do iigrid=1,igridstail; igrid=igrids(iigrid);
      ^D&iib^D=idphyb(^D,igrid);
      {do i^DB=-1,1\}
         if (skip_direction([ i^D ])) cycle
         if(neighbor_type(i^D,igrid)==neighbor_fine) call bc_fill_r
      {end do\}
    end do

    ! receiving ghost-cell values from coarser neighbors
    isend_p=0
    ibuf_recv_p=1
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       ^D&iib^D=idphyb(^D,igrid);
       {do i^DB=-1,1\}
          if (skip_direction([ i^D ])) cycle
          if(neighbor_type(i^D,igrid)==neighbor_coarse) call bc_recv_prolong
       {end do\}
    end do
    ! sending ghost-cell values to finer neighbors 
    irecv_p=0
    ibuf_send_p=1
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       ^D&iib^D=idphyb(^D,igrid);
       {do i^DB=-1,1\}
          if (skip_direction([ i^D ])) cycle
          if(neighbor_type(i^D,igrid)==neighbor_fine) call bc_send_prolong
       {end do\}
    end do

    if (nrecv_bc_p>0) then
      call MPI_WAITALL(nrecv_bc_p,recvrequest_p,recvstatus_p,ierrmpi)
    end if
    if (nsend_bc_p>0) then
      call MPI_WAITALL(nsend_bc_p,sendrequest_p,sendstatus_p,ierrmpi)
    end if

    ! fill coarser representative after receipt
    ibuf_recv_p=1
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       ^D&iib^D=idphyb(^D,igrid);
       {do i^DB=-1,1\}
          if (skip_direction([ i^D ])) cycle
          if(neighbor_type(i^D,igrid)==neighbor_coarse) call bc_fill_p
       {end do\}
    end do

    ! do prolongation on the ghost-cell values received from coarser neighbors 
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       block=>ps(igrid)
       ^D&iib^D=idphyb(^D,igrid);
       if (any(neighbor_type(:^D&,igrid)==neighbor_coarse)) then
         NeedProlong=.false.
         {do i^DB=-1,1\}
            if(skip_direction([ i^D ])) cycle
            if(neighbor_type(i^D,igrid)==neighbor_coarse) then
              call bc_prolong
              NeedProlong(i^D)=.true.
            end if
         {end do\}
         if(stagger_grid) then
           ! Ghost cell prolongation for staggered variables
           ! must be done in a specific order.
           ! First the first neighbours, which have 2 indices=0 in 3D
           ! or one index=0 in 2D
           do idims=1,ndim
             i^D=0;
             select case(idims)
            {case(^D)
               do i^D=-1,1,2
                 if (NeedProlong(i^DD)) call bc_prolong_stg(NeedProlong)
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
               if (NeedProlong(i^D)) call bc_prolong_stg(NeedProlong)
             end do
           end do
           i2=0;
           do i3=-1,1,2
             do i1=-1,1,2
               if (NeedProlong(i^D)) call bc_prolong_stg(NeedProlong)
             end do
           end do
           i3=0;
           do i1=-1,1,2
             do i2=-1,1,2
               if (NeedProlong(i^D)) call bc_prolong_stg(NeedProlong)
             end do
           end do
           }
           ! Finally, the corners, that have no index=0
          {do i^D=-1,1,2\}
             if (NeedProlong(i^D)) call bc_prolong_stg(NeedProlong)
          {end do\}
         end if
       end if
    end do

     ! modify normal component of magnetic field to fix divB=0 
    if(bcphys .and. physics_type=='mhd' .and. ndim>1) call phys_boundary_adjust()
    
    if (nwaux>0) call fix_auxiliary
    
    time_bc=time_bc+(MPI_WTIME()-time_bcin)
    
    contains

      !> skip corner if corner ghost cells are not used
      logical function skip_direction(dir)
        integer, intent(in) :: dir(^ND)

        if (all(dir == 0)) then
           skip_direction = .true.
        else if (.not. req_diagonal .and. count(dir /= 0) > 1) then
           skip_direction = .true.
        else
           skip_direction = .false.
        end if
      end function skip_direction

      !> Receive from sibling at same refinement level
      subroutine bc_recv_srl

        ipe_neighbor=neighbor(2,i^D,igrid)
        if (ipe_neighbor/=mype) then
           sizes_srl_recv(i^D)=nwbc*{(ixR_srl_max^D(iib^D,i^D)-ixR_srl_min^D(iib^D,i^D)+1)|*}
           if(stagger_grid) then
             sizes_srl_recv_total(i^D)=sizes_srl_recv(i^D)+sum(sizes_srl_recv_stg(:,i^D))
           else
             sizes_srl_recv_total(i^D)=sizes_srl_recv(i^D)
           end if
           irecv_srl=irecv_srl+1
           itag=(3**^ND+4**^ND)*(igrid-1)+{(i^D+1)*3**(^D-1)+}
           call MPI_IRECV(recvbuffer_srl(ibuf_recv_srl),sizes_srl_recv_total(i^D),MPI_DOUBLE_PRECISION, &
                          ipe_neighbor,itag,icomm,recvrequest_srl(irecv_srl),ierrmpi)
           ibuf_recv_srl=ibuf_recv_srl+sizes_srl_recv_total(i^D)
        end if

      end subroutine bc_recv_srl

      !> Send data to sibling blocks at another cpu
      subroutine bc_send_srl

        ineighbor=neighbor(1,i^D,igrid)
        ipe_neighbor=neighbor(2,i^D,igrid)
        ipole=neighbor_pole(i^D,igrid)

        if (ipe_neighbor/=mype) then
          select case (ipole)
           case (0) ! No pole
             n_i^D=-i^D;
          {case (^D)
             n_i^D=i^D^D%n_i^DD=-i^DD;\}
          end select
          ixS^L=ixS_srl_^L(iib^D,i^D);
          sizes_srl_send(i^D)=nwbc*{(ixS_srl_max^D(iib^D,i^D)-ixS_srl_min^D(iib^D,i^D)+1)|*}
          ibuf_next=ibuf_send_srl+sizes_srl_send(i^D)
          shapes=(/sizes_srl_send(i^D)/)
          sendbuffer_srl(ibuf_send_srl:ibuf_next-1)=reshape(psb(igrid)%w(ixS^S,nwhead:nwtail),shapes)
          ibuf_start=ibuf_next
          if(stagger_grid) then
            do idir=1,ndim
              ixS^L=ixS_srl_stg_^L(idir,i^D);
              ibuf_next=ibuf_start+sizes_srl_send_stg(idir,i^D)
              shapes=(/sizes_srl_send_stg(idir,i^D)/)
              sendbuffer_srl(ibuf_start:ibuf_next-1)=&
                reshape(psb(igrid)%ws(ixS^S,idir),shapes)
              ibuf_start=ibuf_next
            end do
          end if
          itag=(3**^ND+4**^ND)*(ineighbor-1)+{(n_i^D+1)*3**(^D-1)+}
          isend_srl=isend_srl+1
          call MPI_ISEND(sendbuffer_srl(ibuf_send_srl),ibuf_next-ibuf_send_srl,MPI_DOUBLE_PRECISION, &
                         ipe_neighbor,itag,icomm,sendrequest_srl(isend_srl),ierrmpi)
          ibuf_send_srl=ibuf_next
        end if

      end subroutine bc_send_srl

      !> fill siblings ghost cells with received data
      subroutine bc_fill_srl
        double precision :: tmp(ixGs^T)
        integer :: ixS^L,ixR^L,n_i^D,ixSsync^L,ixRsync^L
        integer :: idir, idirect

        ineighbor=neighbor(1,i^D,igrid)
        ipe_neighbor=neighbor(2,i^D,igrid)
        ipole=neighbor_pole(i^D,igrid)
        idirect={abs(i^D)|+}

        !! Now the special treatment of the pole is done here, at the receive step
        if (ipole==0) then    
          ixR^L=ixR_srl_^L(iib^D,i^D);
          if (ipe_neighbor==mype) then
            n_i^D=-i^D;
            ^D&isb^D=idphyb(^D,ineighbor);
            ixS^L=ixS_srl_^L(isb^D,n_i^D);
            psb(igrid)%w(ixR^S,nwhead:nwtail)=psb(ineighbor)%w(ixS^S,nwhead:nwtail)
            if(stagger_grid) then
              do idir=1,ndim
                ixS^L=ixS_srl_stg_^L(idir,n_i^D);
                ixR^L=ixR_srl_stg_^L(idir,i^D);
                if (idirect == 1) then
                  call indices_for_syncing(idir,i^D,ixR^L,ixS^L,ixRsync^L,ixSsync^L) ! Overwrites ixR, ixS
                  psb(igrid)%ws(ixRsync^S,idir) = half*(psb(igrid)%ws(ixRsync^S,idir)+psb(ineighbor)%ws(ixSsync^S,idir))
                end if
                psb(igrid)%ws(ixR^S,idir) = psb(ineighbor)%ws(ixS^S,idir)
              end do
            end if
          else
            !! Unpack the buffer and fill the ghost cells
            sizes_srl_recv(i^D)=nwbc*{(ixR_srl_max^D(iib^D,i^D)-ixR_srl_min^D(iib^D,i^D)+1)|*}
            ibuf_next=ibuf_recv_srl+sizes_srl_recv(i^D)
            psb(igrid)%w(ixR^S,nwhead:nwtail)=reshape(source=recvbuffer_srl(ibuf_recv_srl:ibuf_next-1),shape=shape(psb(igrid)%w(ixR^S,nwhead:nwtail)))
            ibuf_recv_srl=ibuf_next
            if(stagger_grid) then
              n_i^D=-i^D;
              do idir=1,ndim
                ixS^L=ixS_srl_stg_^L(idir,n_i^D);
                ixR^L=ixR_srl_stg_^L(idir,i^D);
                ibuf_next=ibuf_recv_srl+sizes_srl_recv_stg(idir,i^D)
                tmp(ixS^S) = reshape(source=recvbuffer_srl(ibuf_recv_srl:ibuf_next-1),shape=shape(psb(igrid)%ws(ixS^S,idir)))       
                if (idirect == 1) then
                   ! ixR ixS maybe changed
                   call indices_for_syncing(idir,i^D,ixR^L,ixS^L,ixRsync^L,ixSsync^L) ! Overwrites ixR, ixS
                   psb(igrid)%ws(ixRsync^S,idir) = half*(tmp(ixSsync^S) + psb(igrid)%ws(ixRsync^S,idir))
                end if
                psb(igrid)%ws(ixR^S,idir) = tmp(ixS^S)
                ibuf_recv_srl=ibuf_next
              end do
            end if
          end if

        else ! There is a pole
          select case (ipole)
          {case (^D)
             n_i^D=i^D^D%n_i^DD=-i^DD;\}
          end select
          ixR^L=ixR_srl_^L(iib^D,i^D);
        
          if (ipe_neighbor==mype) then
            ^D&isb^D=idphyb(^D,ineighbor);
            ixS^L=ixS_srl_^L(isb^D,n_i^D);
            !! Fill ghost cells
            call pole_copy(psb(igrid),ixR^L,psb(ineighbor),ixS^L)
            if(stagger_grid) then
              do idir=1,ndim
                ixR^L=ixR_srl_stg_^L(idir,i^D);
                ixS^L=ixS_srl_stg_^L(idir,n_i^D);
                !! Fill ghost cells
                call pole_copy_stg(psb(igrid)%ws,ixR^L,psb(ineighbor)%ws,ixS^L,idir)
              end do
            end if
          else
            ixS^L=ixS_srl_^L(iib^D,n_i^D);
            !! Unpack the buffer and fill an auxiliary array
            sizes_srl_recv(i^D)=nwbc*{(ixR_srl_max^D(iib^D,i^D)-ixR_srl_min^D(iib^D,i^D)+1)|*}
            ibuf_next=ibuf_recv_srl+sizes_srl_recv(i^D)
            pole_buf%w=zero
            pole_buf%w(ixS^S,nwhead:nwtail)=reshape(source=recvbuffer_srl(ibuf_recv_srl:ibuf_next-1),&
                     shape=shape(psb(igrid)%w(ixS^S,nwhead:nwtail)))
            ibuf_recv_srl=ibuf_next
         
            !! Fill ghost cells
            call pole_copy(psb(igrid),ixR^L,pole_buf,ixS^L)
        
            if(stagger_grid) then
              pole_buf%ws=zero
              do idir=1,ndim
                ixR^L=ixR_srl_stg_^L(idir,i^D);
                ixS^L=ixS_srl_stg_^L(idir,n_i^D);
                ibuf_next=ibuf_recv_srl+sizes_srl_recv_stg(idir,i^D)
                pole_buf%ws(ixS^S,idir)=reshape(source=recvbuffer_srl(ibuf_recv_srl:ibuf_next-1),&
                  shape=shape(psb(igrid)%ws(ixS^S,idir)))
                ibuf_recv_srl=ibuf_next
                call pole_copy_stg(psb(igrid)%ws,ixR^L,pole_buf%ws,ixS^L,idir)
              end do
            end if
          end if
        end if
      
      end subroutine bc_fill_srl

      subroutine indices_for_syncing(idir,i^D,ixR^L,ixS^L,ixRsync^L,ixSsync^L)
        integer, intent(in)       :: i^D,idir
        integer, intent(inout)    :: ixR^L,ixS^L
        integer, intent(out)      :: ixRsync^L,ixSsync^L
      
        ixRsync^L=ixR^L;
        ixSsync^L=ixS^L;
        
        {
        if (i^D == -1 .and. idir == ^D) then
           ixRsyncmin^D = ixRmax^D
           ixRsyncmax^D = ixRmax^D
           ixSsyncmin^D = ixSmax^D
           ixSsyncmax^D = ixSmax^D
           ixRmax^D = ixRmax^D - 1
           ixSmax^D = ixSmax^D - 1
        else if (i^D == 1 .and. idir == ^D) then
           ixRsyncmin^D = ixRmin^D
           ixRsyncmax^D = ixRmin^D
           ixSsyncmin^D = ixSmin^D
           ixSsyncmax^D = ixSmin^D
           ixRmin^D = ixRmin^D + 1
           ixSmin^D = ixSmin^D + 1
        end if
        \}

      end subroutine indices_for_syncing

      !> Receive from fine neighbor
      subroutine bc_recv_restrict

        {do ic^DB=1+int((1-i^DB)/2),2-int((1+i^DB)/2)
           inc^DB=2*i^DB+ic^DB\}
           ipe_neighbor=neighbor_child(2,inc^D,igrid)
           if (ipe_neighbor/=mype) then
              itag=(3**^ND+4**^ND)*(igrid-1)+3**^ND+{inc^D*4**(^D-1)+}
              sizes_r_recv(inc^D)=nwbc*{(ixR_r_max^D(iib^D,inc^D)-ixR_r_min^D(iib^D,inc^D)+1)|*}
              if(stagger_grid) then
                sizes_r_recv_total(inc^D)=sizes_r_recv(inc^D)+sum(sizes_r_recv_stg(:,inc^D))
              else
                sizes_r_recv_total(inc^D)=sizes_r_recv(inc^D)
              end if
              irecv_r=irecv_r+1
              call MPI_IRECV(recvbuffer_r(ibuf_recv_r),sizes_r_recv_total(inc^D), &
                             MPI_DOUBLE_PRECISION,ipe_neighbor,itag, &
                             icomm,recvrequest_r(irecv_r),ierrmpi)
              ibuf_recv_r=ibuf_recv_r+sizes_r_recv_total(inc^D)
           end if
        {end do\}

      end subroutine bc_recv_restrict

      !> Send to coarser neighbor
      subroutine bc_send_restrict
        integer :: ii^D

        ic^D=1+modulo(node(pig^D_,igrid)-1,2);
        if ({.not.(i^D==0.or.i^D==2*ic^D-3)|.or.}) return
        if(phyboundblock(igrid)) then
          ! filling physical boundary ghost cells of a coarser representative block for
          ! sending swap region with width of nghostcells to its coarser neighbor
          do idims=1,ndim
             ! to avoid using as yet unknown corner info in more than 1D, we
             ! fill only interior mesh ranges of the ghost cell ranges at first,
             ! and progressively enlarge the ranges to include corners later
             {kmin^D=merge(0, 1, idims==^D)
             kmax^D=merge(0, 1, idims==^D)
             ixBmin^D=ixCoGmin^D+kmin^D*nghostcells
             ixBmax^D=ixCoGmax^D-kmax^D*nghostcells\}
             {^IFTWOD
             if(idims > 1 .and. neighbor_type(-1,0,igrid)==neighbor_boundary) ixBmin1=ixCoGmin1
             if(idims > 1 .and. neighbor_type( 1,0,igrid)==neighbor_boundary) ixBmax1=ixCoGmax1}
             {^IFTHREED
             if(idims > 1 .and. neighbor_type(-1,0,0,igrid)==neighbor_boundary) ixBmin1=ixCoGmin1
             if(idims > 1 .and. neighbor_type( 1,0,0,igrid)==neighbor_boundary) ixBmax1=ixCoGmax1
             if(idims > 2 .and. neighbor_type(0,-1,0,igrid)==neighbor_boundary) ixBmin2=ixCoGmin2
             if(idims > 2 .and. neighbor_type(0, 1,0,igrid)==neighbor_boundary) ixBmax2=ixCoGmax2}
             {if(i^D==-1) then
               ixBmin^D=ixCoGmin^D+nghostcells
               ixBmax^D=ixCoGmin^D+2*nghostcells-1
             else if(i^D==1) then
               ixBmin^D=ixCoGmax^D-2*nghostcells+1
               ixBmax^D=ixCoGmax^D-nghostcells
             end if\}
             do iside=1,2
                ii^D=kr(^D,idims)*(2*iside-3);
                if ({abs(i^D)==1.and.abs(ii^D)==1|.or.}) cycle
                if (neighbor_type(ii^D,igrid)/=neighbor_boundary) cycle
                call bc_phys(iside,idims,time,0.d0,psc(igrid)%w,&
                       psc(igrid)%x,ixCoG^L,ixB^L)
             end do
          end do
        end if

        ineighbor=neighbor(1,i^D,igrid)
        ipe_neighbor=neighbor(2,i^D,igrid)
        ipole=neighbor_pole(i^D,igrid)

        if (ipe_neighbor/=mype) then
          inc^D=i^D+ic^D;
          ! The ghost region of the neighbor changes
          ! if there is a pole
          select case (ipole)
          case(0) ! No pole
            n_inc^D=-2*i^D+ic^D;
         {case (^D) ! Pole in this direction
            n_inc^D=2*i^D+(3-ic^D)^D%n_inc^DD=-2*i^DD+ic^DD;\}
          end select
        
          ! fill corresponding part of the send buffer...
          ixS^L=ixS_r_^L(iib^D,i^D);
          sizes_r_send(i^D)=nwbc*{(ixS_r_max^D(iib^D,i^D)-ixS_r_min^D(iib^D,i^D)+1)|*}
          ibuf_next=ibuf_send_r+sizes_r_send(i^D)
          shapes=(/sizes_r_send(i^D)/)
          sendbuffer_r(ibuf_send_r:ibuf_next-1)=reshape(psc(igrid)%w(ixS^S,nwhead:nwtail),shapes)

          if(stagger_grid) then 
            ibuf_start=ibuf_next
            do idir=1,ndim
              ixS^L=ixS_r_stg_^L(idir,i^D);
              ibuf_next=ibuf_start+sizes_r_send_stg(idir,i^D)
              shapes=(/sizes_r_send_stg(idir,i^D)/)
              sendbuffer_r(ibuf_start:ibuf_next-1)=&
                reshape(psc(igrid)%ws(ixS^S,idir),shapes)
              ibuf_start=ibuf_next
            end do
          end if

          itag=(3**^ND+4**^ND)*(ineighbor-1)+3**^ND+{n_inc^D*4**(^D-1)+}
          isend_r=isend_r+1
          call MPI_ISEND(sendbuffer_r(ibuf_send_r),ibuf_next-ibuf_send_r, &
                         MPI_DOUBLE_PRECISION,ipe_neighbor,itag, &
                         icomm,sendrequest_r(isend_r),ierrmpi)
        
          ibuf_send_r=ibuf_next
        
        end if

      end subroutine bc_send_restrict

      !> fill restricted ghost cells after receipt
      subroutine bc_fill_r

        ipole=neighbor_pole(i^D,igrid)
        if (ipole==0) then
          ! Loop over the children ic^D to and their neighbors inc^D
          {do ic^DB=1+int((1-i^DB)/2),2-int((1+i^DB)/2)
             inc^DB=2*i^DB+ic^DB\}
             n_i^D=-i^D;
             ineighbor=neighbor_child(1,inc^D,igrid)
             ipe_neighbor=neighbor_child(2,inc^D,igrid)
             ixR^L=ixR_r_^L(iib^D,inc^D);
             if (ipe_neighbor==mype) then ! Same processor
               ^D&isb^D=idphyb(^D,ineighbor);
               ixS^L=ixS_r_^L(isb^D,n_i^D);
               psb(igrid)%w(ixR^S,nwhead:nwtail)=psc(ineighbor)%w(ixS^S,nwhead:nwtail)
               if(stagger_grid) then
                 do idir=1,ndim
                    ixS^L=ixS_r_stg_^L(idir,n_i^D);
                    ixR^L=ixR_r_stg_^L(idir,inc^D);
                    psb(igrid)%ws(ixR^S,idir)=psc(ineighbor)%ws(ixS^S,idir)
                 end do
               end if
             else ! Different processor
               !! Unpack the buffer and fill the ghost cells
               sizes_r_recv(inc^D)=nwbc*{(ixR_r_max^D(iib^D,inc^D)-ixR_r_min^D(iib^D,inc^D)+1)|*}
               ibuf_next=ibuf_recv_r+sizes_r_recv(inc^D)
               psb(igrid)%w(ixR^S,nwhead:nwtail)=reshape(source=recvbuffer_r(ibuf_recv_r:ibuf_next-1),&
                     shape=shape(psb(igrid)%w(ixR^S,nwhead:nwtail)))
               ibuf_recv_r=ibuf_next

               if(stagger_grid) then
                 do idir=1,ndim
                   ixR^L=ixR_r_stg_^L(idir,inc^D);
                   ibuf_next=ibuf_recv_r+sizes_r_recv_stg(idir,inc^D)
                   psb(igrid)%ws(ixR^S,idir)=reshape(source=recvbuffer_r(ibuf_recv_r:ibuf_next-1),&
                         shape=shape(psb(igrid)%ws(ixR^S,idir)))
                   ibuf_recv_r=ibuf_next
                 end do
               end if
             end if
          {end do\}
        
        else !! There is a pole
          {do ic^DB=1+int((1-i^DB)/2),2-int((1+i^DB)/2)
             inc^DB=2*i^DB+ic^DB\}
             select case(ipole)
            {case (^D)
               n_i^D=i^D^D%n_i^DD=-i^DD;\}
             end select
             ineighbor=neighbor_child(1,inc^D,igrid)
             ipe_neighbor=neighbor_child(2,inc^D,igrid)
             ixR^L=ixR_r_^L(iib^D,inc^D);
             if (ipe_neighbor==mype) then ! Same processor
               ^D&isb^D=idphyb(^D,ineighbor);
               ixS^L=ixS_r_^L(isb^D,n_i^D);
               !! Fill ghost cells
               call pole_copy(psb(igrid),ixR^L,psc(ineighbor),ixS^L)
               if(stagger_grid) then
                 do idir=1,ndim
                   ixS^L=ixS_r_stg_^L(idir,n_i^D);
                   ixR^L=ixR_r_stg_^L(idir,inc^D);
                   !! Fill ghost cells
                   call pole_copy_stg(psb(igrid)%ws,ixR^L,psc(ineighbor)%ws,ixS^L,idir)
                 end do
               end if

             else ! Different processor
               !! Unpack the buffer and fill an auxiliary array
               sizes_r_recv(inc^D)=nwbc*{(ixR_r_max^D(iib^D,inc^D)-ixR_r_min^D(iib^D,inc^D)+1)|*}
               ibuf_next=ibuf_recv_r+sizes_r_recv(inc^D)
               pole_buf%w=zero
               pole_buf%w(ixR^S,nwhead:nwtail)=&
                 reshape(source=recvbuffer_r(ibuf_recv_r:ibuf_next-1),&
                        shape=shape(psb(igrid)%w(ixR^S,nwhead:nwtail)))
               ibuf_recv_r=ibuf_next
           
               !! Fill ghost cells
               call pole_copy(psb(igrid),ixR^L,pole_buf,ixR^L)

               if(stagger_grid) then
                 pole_buf%ws=zero
                 do idir=1,ndim
                    ixS^L=ixS_r_stg_^L(idir,n_i^D);
                    ixR^L=ixR_r_stg_^L(idir,inc^D);
                    ibuf_next=ibuf_recv_r+sizes_r_recv_stg(idir,inc^D)
                    pole_buf%ws(ixR^S,idir)=reshape(source=recvbuffer_r(ibuf_recv_r:ibuf_next-1),&
                      shape=shape(psb(igrid)%ws(ixR^S,idir)))
                    call pole_copy_stg(psb(igrid)%ws,ixR^L,pole_buf%ws,ixR^L,idir)
                    ibuf_recv_r=ibuf_next
                 end do
               end if

             end if
          {end do\}
        
        end if

      end subroutine bc_fill_r

      !> Receive from coarse neighbor
      subroutine bc_recv_prolong

        ic^D=1+modulo(node(pig^D_,igrid)-1,2);
        if ({.not.(i^D==0.or.i^D==2*ic^D-3)|.or.}) return

        ipe_neighbor=neighbor(2,i^D,igrid)
        if (ipe_neighbor/=mype) then
           inc^D=ic^D+i^D;
           itag=(3**^ND+4**^ND)*(igrid-1)+3**^ND+{inc^D*4**(^D-1)+}
           sizes_p_recv(inc^D)=nwbc*{(ixR_p_max^D(iib^D,inc^D)-ixR_p_min^D(iib^D,inc^D)+1)|*}
           if(stagger_grid) then
             sizes_p_recv_total(inc^D)=sizes_p_recv(inc^D)+sum(sizes_p_recv_stg(:,inc^D))
           else
             sizes_p_recv_total(inc^D)=sizes_p_recv(inc^D)
           end if
           irecv_p=irecv_p+1
           call MPI_IRECV(recvbuffer_p(ibuf_recv_p),sizes_p_recv_total(inc^D),&                                                                            
                          MPI_DOUBLE_PRECISION,ipe_neighbor,itag,&
                          icomm,recvrequest_p(irecv_p),ierrmpi)
           ibuf_recv_p=ibuf_recv_p+sizes_p_recv_total(inc^D)
        end if

      end subroutine bc_recv_prolong

      !> Send to finer neighbor
      subroutine bc_send_prolong
        integer :: ii^D

        ipole=neighbor_pole(i^D,igrid)
        {do ic^DB=1+int((1-i^DB)/2),2-int((1+i^DB)/2)
           inc^DB=2*i^DB+ic^DB\}

           ineighbor=neighbor_child(1,inc^D,igrid)
           ipe_neighbor=neighbor_child(2,inc^D,igrid)
        
           if (ipe_neighbor/=mype) then
             select case (ipole)
             case(0) ! No pole
                n_i^D=-i^D;
                n_inc^D=ic^D+n_i^D;
            {case (^D) ! Pole in this direction
                n_inc^D=inc^D^D%n_inc^DD=ic^DD-i^DD;
                \}
             end select
             ! fill corresponding part of the send buffer...
             ixS^L=ixS_p_^L(iib^D,inc^D);
             sizes_p_send(inc^D)=nwbc*{(ixS_p_max^D(iib^D,inc^D)-ixS_p_min^D(iib^D,inc^D)+1)|*}
             ibuf_next=ibuf_send_p+sizes_p_send(inc^D)
             shapes=(/sizes_p_send(inc^D)/)
             sendbuffer_p(ibuf_send_p:ibuf_next-1)=reshape(psb(igrid)%w(ixS^S,nwhead:nwtail),shapes)
             if(stagger_grid) then 
               ibuf_start=ibuf_next
               do idir=1,ndim
                 ixS^L=ixS_p_stg_^L(idir,inc^D);
                 ibuf_next=ibuf_start+sizes_p_send_stg(idir,inc^D)
                 shapes=(/sizes_p_send_stg(idir,inc^D)/)
                 sendbuffer_p(ibuf_start:ibuf_next-1)=&
                   reshape(psb(igrid)%ws(ixS^S,idir),shapes)   
                 ibuf_start=ibuf_next
               end do
             end if
             ! ...and send
             itag=(3**^ND+4**^ND)*(ineighbor-1)+3**^ND+{n_inc^D*4**(^D-1)+}
             isend_p=isend_p+1
             call MPI_ISEND(sendbuffer_p(ibuf_send_p),ibuf_next-ibuf_send_p, &
                            MPI_DOUBLE_PRECISION,ipe_neighbor,itag, &
                            icomm,sendrequest_p(isend_p),ierrmpi)
             ibuf_send_p=ibuf_next
           end if

        {end do\}

      end subroutine bc_send_prolong

      !> fill coarser representative with data from coarser neighbors
      subroutine bc_fill_p
        ic^D=1+modulo(node(pig^D_,igrid)-1,2);
        if ({.not.(i^D==0.or.i^D==2*ic^D-3)|.or.}) return

        ineighbor=neighbor(1,i^D,igrid)
        ipe_neighbor=neighbor(2,i^D,igrid)
        ipole=neighbor_pole(i^D,igrid)

        if (ipole==0) then   !! There is no pole 
          inc^D=ic^D+i^D;
          ixR^L=ixR_p_^L(iib^D,inc^D);
          if(ipe_neighbor==mype) then !! Same processor
            n_inc^D=-2*i^D+ic^D;
            ^D&isb^D=idphyb(^D,ineighbor);
            ixS^L=ixS_p_^L(isb^D,n_inc^D);
            psc(igrid)%w(ixR^S,nwhead:nwtail) &
                    =psb(ineighbor)%w(ixS^S,nwhead:nwtail)
            if(stagger_grid) then
              do idir=1,ndim
                ixS^L=ixS_p_stg_^L(idir,n_inc^D);
                ixR^L=ixR_p_stg_^L(idir,inc^D);
                psc(igrid)%ws(ixR^S,idir)=psb(ineighbor)%ws(ixS^S,idir)
              end do
            end if
          else !! Different processor
            !! Unpack the buffer and fill the ghost cells
            sizes_p_recv(inc^D)=nwbc*{(ixR_p_max^D(iib^D,inc^D)-ixR_p_min^D(iib^D,inc^D)+1)|*}
            ibuf_next=ibuf_recv_p+sizes_p_recv(inc^D)
            psc(igrid)%w(ixR^S,nwhead:nwtail)=reshape(source=recvbuffer_p(ibuf_recv_p:ibuf_next-1),&
                  shape=shape(psc(igrid)%w(ixR^S,nwhead:nwtail)))
            ibuf_recv_p=ibuf_next
            if(stagger_grid) then
              do idir=1,^ND
                ixR^L=ixR_p_stg_^L(idir,inc^D);
                ibuf_next=ibuf_recv_p+sizes_p_recv_stg(idir,inc^D)
                psc(igrid)%ws(ixR^S,idir)=reshape(source=recvbuffer_p(ibuf_recv_p:ibuf_next-1),&
                      shape=shape(psc(igrid)%ws(ixR^S,idir)))
                ibuf_recv_p=ibuf_next
              end do
            end if
          end if

        else !! There is a pole
          inc^D=ic^D+i^D;
          select case (ipole)
          {case (^D)
             n_inc^D=2*i^D+(3-ic^D)^D%n_inc^DD=-2*i^DD+ic^DD;\}
          end select

          ixR^L=ixR_p_^L(iib^D,inc^D);
          if (ipe_neighbor==mype) then
            ^D&isb^D=idphyb(^D,ineighbor);
            ixS^L=ixS_p_^L(isb^D,n_inc^D);
            call pole_copy(psc(igrid),ixR^L,psb(ineighbor),ixS^L)
            if(stagger_grid) then
              do idir=1,ndim
                ixS^L=ixS_p_stg_^L(idir,n_inc^D);
                ixR^L=ixR_p_stg_^L(idir,inc^D);
                call pole_copy_stg(psc(igrid)%ws,ixR^L,psb(ineighbor)%ws,ixS^L,idir)
              end do
            end if
          else
            !! Unpack the buffer and fill an auxiliary array
            sizes_p_recv(inc^D)=nwbc*{(ixR_p_max^D(iib^D,inc^D)-ixR_p_min^D(iib^D,inc^D)+1)|*}
            ibuf_next=ibuf_recv_p+sizes_p_recv(inc^D)
            pole_buf%w=zero
            pole_buf%w(ixR^S,nwhead:nwtail)=reshape(source=recvbuffer_p(ibuf_recv_p:ibuf_next-1),&
                     shape=shape(psc(igrid)%w(ixR^S,nwhead:nwtail)))
            ibuf_recv_p=ibuf_next
            !! Fill ghost cells
            call pole_copy(psc(igrid),ixR^L,pole_buf,ixR^L)

            if(stagger_grid) then
              pole_buf%ws=zero
              do idir=1,ndim
                ixR^L=ixR_p_stg_^L(idir,inc^D);
                ibuf_next=ibuf_recv_p+sizes_p_recv_stg(idir,inc^D)
                pole_buf%ws(ixR^S,idir)=reshape(source=recvbuffer_p(ibuf_recv_p:ibuf_next-1),&
                  shape=shape(psc(igrid)%ws(ixR^S,idir)))
                ibuf_recv_p=ibuf_next
                call pole_copy_stg(psc(igrid)%ws,ixR^L,pole_buf%ws,ixR^L,idir)
              end do
            end if
          end if

        end if

      end subroutine bc_fill_p

      subroutine bc_prolong
        use mod_physics, only: phys_to_primitive, phys_to_conserved

        integer :: ixFi^L,ixCo^L,ii^D
        double precision :: dxFi^D, dxCo^D, xFimin^D, xComin^D, invdxCo^D

        ixFi^L=ixR_srl_^L(iib^D,i^D);
        dxFi^D=rnode(rpdx^D_,igrid);
        dxCo^D=two*dxFi^D;
        invdxCo^D=1.d0/dxCo^D;

        ! compute the enlarged grid lower left corner coordinates
        ! these are true coordinates for an equidistant grid, 
        ! but we can temporarily also use them for getting indices 
        ! in stretched grids
        xFimin^D=rnode(rpxmin^D_,igrid)-dble(nghostcells)*dxFi^D;
        xComin^D=rnode(rpxmin^D_,igrid)-dble(nghostcells)*dxCo^D;

        if(prolongprimitive) then
           ! following line again assumes equidistant grid, but 
           ! just computes indices, so also ok for stretched case
           ! reason for +1-1 and +1+1: the coarse representation has 
           ! also nghostcells at each side. During
           ! prolongation, we need cells to left and right, hence -1/+1
           ixComin^D=int((xFimin^D+(dble(ixFimin^D)-half)*dxFi^D-xComin^D)*invdxCo^D)+1-1;
           ixComax^D=int((xFimin^D+(dble(ixFimax^D)-half)*dxFi^D-xComin^D)*invdxCo^D)+1+1;
           call phys_to_primitive(ixCoG^L,ixCo^L,&
             psc(igrid)%w,psc(igrid)%x)
        endif

        select case (typeghostfill)
        case ("linear")
           call interpolation_linear(ixFi^L,dxFi^D,xFimin^D,dxCo^D,invdxCo^D,xComin^D)
        case ("copy")
           call interpolation_copy(ixFi^L,dxFi^D,xFimin^D,dxCo^D,invdxCo^D,xComin^D)
        case default
           write (unitterm,*) "Undefined typeghostfill ",typeghostfill
           call mpistop("Undefined typeghostfill")
        end select

        if(prolongprimitive) call phys_to_conserved(ixCoG^L,ixCo^L,&
             psc(igrid)%w,psc(igrid)%x)

      end subroutine bc_prolong

      subroutine bc_prolong_stg(NeedProlong)
        use mod_amr_fct
        logical,dimension(-1:1^D&) :: NeedProlong
        logical                    :: fine_^Lin
        integer                    :: ixFi^L,ixCo^L
        double precision           :: dxFi^D,dxCo^D,xFimin^D,xComin^D,invdxCo^D
        ! Check what is already at the desired level
        fine_^Lin=.false.;
        {
        if(i^D>-1) fine_min^Din=.not.NeedProlong(i^DD-kr(^D,^DD))
        if(i^D<1)  fine_max^Din=.not.NeedProlong(i^DD+kr(^D,^DD))
        \}

        ixFi^L=ixR_srl_^L(iib^D,i^D);

        dxFi^D=rnode(rpdx^D_,igrid);
        dxCo^D=two*dxFi^D;
        invdxCo^D=1.d0/dxCo^D;

        xFimin^D=rnode(rpxmin^D_,igrid)-dble(nghostcells)*dxFi^D;
        xComin^D=rnode(rpxmin^D_,igrid)-dble(nghostcells)*dxCo^D;

        ! moved the physical boundary filling here, to only fill the
        ! part needed

        ixComin^D=int((xFimin^D+(dble(ixFimin^D)-half)*dxFi^D-xComin^D)*invdxCo^D)+1-1;
        ixComax^D=int((xFimin^D+(dble(ixFimax^D)-half)*dxFi^D-xComin^D)*invdxCo^D)+1+1;

        call prolong_2nd_stg(psc(igrid),ps(igrid),ixCo^L,ixFi^L,dxCo^D,xComin^D,dxFi^D,xFimin^D,.true.,fine_^Lin)

        ! The current region has already been refined, so it doesn t need to be prolonged again
         NeedProlong(i^D)=.false. 

      end subroutine bc_prolong_stg

      subroutine interpolation_linear(ixFi^L,dxFi^D,xFimin^D, &
                                      dxCo^D,invdxCo^D,xComin^D)
        use mod_physics, only: phys_to_conserved
        integer, intent(in) :: ixFi^L
        double precision, intent(in) :: dxFi^D, xFimin^D,dxCo^D, invdxCo^D, xComin^D

        integer :: ixCo^D, jxCo^D, hxCo^D, ixFi^D, ix^D, iw, idims, nwmin,nwmax
        double precision :: xCo^D, xFi^D, eta^D
        double precision :: slopeL, slopeR, slopeC, signC, signR
        double precision :: slope(1:nw,ndim)
        !!double precision :: local_invdxCo^D
        double precision :: signedfactorhalf^D
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

           if(slab) then
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
             {signedfactorhalf^D=(xFi^D-xCo^D)*invdxCo^D*two
              if(dabs(signedfactorhalf^D**2-1.0d0/4.0d0)>smalldouble) call mpistop("error in bc_prolong")
              eta^D=signedfactorhalf^D*(one-block%dvolume(ixFi^DD) &
                   /sum(block%dvolume(ix^D:ix^D+1^D%ixFi^DD))) \}
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
                 select case(typeprolonglimit)
                 case('unlimit')
                   slope(iw,idims)=slopeC
                 case('minmod')
                   slope(iw,idims)=signR*max(zero,min(dabs(slopeR), &
                                                     signR*slopeL))
                 case('woodward')
                   slope(iw,idims)=two*signR*max(zero,min(dabs(slopeR), &
                                      signR*slopeL,signR*half*slopeC))
                 case('koren')
                   slope(iw,idims)=signR*max(zero,min(two*signR*slopeL, &
                    (dabs(slopeR)+two*slopeL*signR)*third,two*dabs(slopeR)))
                 case default
                   slope(iw,idims)=signC*max(zero,min(dabs(slopeC), &
                                     signC*slopeL,signC*slopeR))
                 end select
              end do
           end do
        
           ! Interpolate from coarse cell using limited slopes
           psb(igrid)%w(ixFi^D,nwmin:nwmax)=psc(igrid)%w(ixCo^D,nwmin:nwmax)+&
             {(slope(nwmin:nwmax,^D)*eta^D)+}
        
        {end do\}
        
        if(prolongprimitive) call phys_to_conserved(ixG^LL,ixFi^L,psb(igrid)%w,psb(igrid)%x)
      
      end subroutine interpolation_linear

      subroutine interpolation_copy(ixFi^L,dxFi^D,xFimin^D, &
                                    dxCo^D,invdxCo^D,xComin^D)
        use mod_physics, only: phys_to_conserved
        integer, intent(in) :: ixFi^L
        double precision, intent(in) :: dxFi^D, xFimin^D,dxCo^D, invdxCo^D, xComin^D

        integer :: ixCo^D, ixFi^D, nwmin,nwmax
        double precision :: xFi^D

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

      subroutine pole_copy(psrecv,ixR^L,pssend,ixS^L)
      
        integer, intent(in) :: ixR^L,ixS^L
        type(state) :: psrecv, pssend

        integer :: iw

        select case (ipole)
        {case (^D)
           iside=int((i^D+3)/2)
           iB=2*(^D-1)+iside
           do iw=nwhead,nwtail
             select case (typeboundary(iw,iB))
             case ("symm")
               psrecv%w(ixR^S,iw) = pssend%w(ixSmax^D:ixSmin^D:-1^D%ixS^S,iw)
             case ("asymm")
               psrecv%w(ixR^S,iw) =-pssend%w(ixSmax^D:ixSmin^D:-1^D%ixS^S,iw)
             case default
               call mpistop("Pole boundary condition should be symm or asymm")
             end select
           end do \}
        end select
      
      end subroutine pole_copy

      subroutine pole_copy_stg(wrecv,ixR^L,wsend,ixS^L,idir)
      
        integer, intent(in) :: ixR^L,ixS^L,idir
        double precision :: wrecv(ixGs^T,1:nws), wsend(ixGs^T,1:nws)

        select case (ipole)
        {case (^D)
           iside=int((i^D+3)/2)
           iB=2*(^D-1)+iside
           select case (typeboundary(iw_s0+idir,iB))
           case ("symm")
             wrecv(ixR^S,idir) = wsend(ixSmax^D:ixSmin^D:-1^D%ixS^S,idir)
           case ("asymm")
             wrecv(ixR^S,idir) =-wsend(ixSmax^D:ixSmin^D:-1^D%ixS^S,idir)
           case default
             call mpistop("Pole boundary condition should be symm or asymm")
           end select
         \}
        end select

      end subroutine pole_copy_stg

      subroutine pole_buffer(wrecv,ixIR^L,ixR^L,wsend,ixIS^L,ixS^L)
      
        integer, intent(in) :: ixIR^L,ixR^L,ixIS^L,ixS^L
        double precision :: wrecv(ixIR^S,nwhead:nwtail), wsend(ixIS^S,1:nw)

        integer :: iw

        select case (ipole)
        {case (^D)
           iside=int((i^D+3)/2)
           iB=2*(^D-1)+iside
           do iw=nwhead,nwtail
             select case (typeboundary(iw,iB))
             case ("symm")
               wrecv(ixR^S,iw) = wsend(ixSmax^D:ixSmin^D:-1^D%ixS^S,iw)
             case ("asymm")
               wrecv(ixR^S,iw) =-wsend(ixSmax^D:ixSmin^D:-1^D%ixS^S,iw)
             case default
               call mpistop("Pole boundary condition should be symm or asymm")
             end select
           end do \}
        end select
      
      end subroutine pole_buffer

      subroutine fix_auxiliary
        use mod_physics, only: phys_get_aux
      
        integer :: ix^L

        do iigrid=1,igridstail; igrid=igrids(iigrid);
          saveigrid=igrid
          block=>psb(igrid)
          call identifyphysbound(psb(igrid),iib^D)   
          {do i^DB=-1,1\}
             if (skip_direction([ i^D ])) cycle
             ix^L=ixR_srl_^L(iib^D,i^D);
             call phys_get_aux(.true.,psb(igrid)%w,ps(igrid)%x,ixG^L,ix^L,"bc")
          {end do\}
        end do
      
      end subroutine fix_auxiliary

  end subroutine getbc

  subroutine identifyphysbound(s,iib^D)
    use mod_global_parameters

    type(state)          :: s
    integer, intent(out) :: iib^D

    {
    if(s%is_physical_boundary(2*^D) .and. &
       s%is_physical_boundary(2*^D-1)) then
      iib^D=2
    else if(s%is_physical_boundary(2*^D-1)) then
      iib^D=-1
    else if(s%is_physical_boundary(2*^D)) then
      iib^D=1
    else
      iib^D=0
    end if
    \}

  end subroutine identifyphysbound

end module mod_ghostcells_update
