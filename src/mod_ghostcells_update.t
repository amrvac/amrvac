!> update ghost cells of all blocks including physical boundaries 
module mod_ghostcells_update

  implicit none
  ! Special buffer for pole copy
  type wbuffer
    double precision, dimension(:^D&,:), allocatable :: w
  end type wbuffer

  ! A switch of update physical boundary or not
  logical :: bcphys=.true.
  integer :: ixM^L, ixCoG^L, ixCoM^L

  !> The number of interleaving sending buffers for ghost cells
  integer, parameter :: npwbuf=2

  ! index ranges to send (S) to sibling blocks, receive (R) from 
  ! sibling blocks, send restricted (r) ghost cells to coarser blocks 
  integer, dimension(-1:1,-1:1) :: ixS_srl_^L, ixR_srl_^L, ixS_r_^L

  ! index ranges to receive restriced ghost cells from finer blocks, 
  ! send prolongated (p) ghost cells to finer blocks, receive prolongated 
  ! ghost from coarser blocks
  integer, dimension(-1:1, 0:3) :: ixR_r_^L, ixS_p_^L, ixR_p_^L

  ! MPI derived datatype to send and receive subarrays of ghost cells to
  ! neighbor blocks in a different processor.
  !
  ! The first index goes from -1:1, where -1 is used when a block touches the
  ! lower boundary, 1 when a block touches an upper boundary, and 0 a situation
  ! away from boundary conditions.
  !
  ! There are two variants, _f indicates that all flux variables are filled,
  ! whereas _p means that part of the variables is filled (currently used for
  ! energy). Furthermore _r_ stands for restrict, _p_ for prolongation.
  integer, dimension(-1:1^D&,-1:1^D&), target :: type_send_srl_f, type_recv_srl_f, type_send_r_f
  integer, dimension(-1:1^D&, 0:3^D&), target :: type_recv_r_f, type_send_p_f, type_recv_p_f
  integer, dimension(-1:1^D&,-1:1^D&), target :: type_send_srl_p, type_recv_srl_p, type_send_r_p
  integer, dimension(-1:1^D&, 0:3^D&), target :: type_recv_r_p, type_send_p_p, type_recv_p_p
  integer, dimension(:^D&,:^D&), pointer :: type_send_srl, type_recv_srl, type_send_r
  integer, dimension(:^D&,:^D&), pointer :: type_recv_r, type_send_p, type_recv_p

contains

  subroutine init_bc()
    use mod_global_parameters 

    integer :: nghostcellsCo, interpolation_order
    integer :: nx^D, nxCo^D, ixG^L, i^D, ic^D, inc^D, iib^D

    ixG^L=ixG^LL;
    ixM^L=ixG^L^LSUBnghostcells;
    ixCoGmin^D=1;
    ixCoGmax^D=ixGmax^D/2+nghostcells;
    ixCoM^L=ixCoG^L^LSUBnghostcells;
    
    nx^D=ixMmax^D-ixMmin^D+1;
    nxCo^D=nx^D/2;
    
    select case (typeghostfill)
    case ("copy")
       interpolation_order=1
    case ("linear","unlimit")
       interpolation_order=2
    case default
       interpolation_order=2
       write (unitterm,*) "Undefined typeghostfill ",typeghostfill
       call mpistop("")
    end select
    nghostcellsCo=int((nghostcells+1)/2)
    
    if (nghostcellsCo+interpolation_order-1>nghostcells) then
       call mpistop("interpolation order for prolongation in getbc to high")
    end if
    
    ! (iib,i) index has following meanings: iib = 0 means it is not at any physical boundary
    ! iib=-1 means it is at the minimum side of a physical boundary  
    ! iib= 1 means it is at the maximum side of a physical boundary  
    ! i=-1 means subregion prepared for the neighbor at its minimum side 
    ! i= 1 means subregion prepared for the neighbor at its maximum side 
    {
    ixS_srl_min^D(:,-1)=ixMmin^D
    ixS_srl_min^D(:, 0)=ixMmin^D
    ixS_srl_min^D(:, 1)=ixMmax^D+1-nghostcells
    ixS_srl_max^D(:,-1)=ixMmin^D-1+nghostcells
    ixS_srl_max^D(:, 0)=ixMmax^D
    ixS_srl_max^D(:, 1)=ixMmax^D
    
    ixS_srl_min^D(-1,0)=1
    ixS_srl_min^D( 1,0)=ixMmin^D
    ixS_srl_max^D(-1,0)=ixMmax^D
    ixS_srl_max^D( 1,0)=ixGmax^D
     
    ixR_srl_min^D(:,-1)=1
    ixR_srl_min^D(:, 0)=ixMmin^D
    ixR_srl_min^D(:, 1)=ixMmax^D+1
    ixR_srl_max^D(:,-1)=nghostcells
    ixR_srl_max^D(:, 0)=ixMmax^D
    ixR_srl_max^D(:, 1)=ixGmax^D
    
    ixR_srl_min^D(-1,0)=1
    ixR_srl_min^D( 1,0)=ixMmin^D
    ixR_srl_max^D(-1,0)=ixMmax^D
    ixR_srl_max^D( 1,0)=ixGmax^D
    
    ixS_r_min^D(:,-1)=ixCoMmin^D
    ixS_r_min^D(:, 0)=ixCoMmin^D
    ixS_r_min^D(:, 1)=ixCoMmax^D+1-nghostcells
    ixS_r_max^D(:,-1)=ixCoMmin^D-1+nghostcells
    ixS_r_max^D(:, 0)=ixCoMmax^D
    ixS_r_max^D(:, 1)=ixCoMmax^D
    
    ixS_r_min^D(-1,0)=1
    ixS_r_min^D( 1,0)=ixCoMmin^D
    ixS_r_max^D(-1,0)=ixCoMmax^D
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
    
    ixS_p_min^D(:, 0)=ixMmin^D-(interpolation_order-1)
    ixS_p_min^D(:, 1)=ixMmin^D-(interpolation_order-1)
    ixS_p_min^D(:, 2)=ixMmin^D+nxCo^D-nghostcellsCo-(interpolation_order-1)
    ixS_p_min^D(:, 3)=ixMmax^D+1-nghostcellsCo-(interpolation_order-1)
    ixS_p_max^D(:, 0)=ixMmin^D-1+nghostcellsCo+(interpolation_order-1)
    ixS_p_max^D(:, 1)=ixMmin^D-1+nxCo^D+nghostcellsCo+(interpolation_order-1)
    ixS_p_max^D(:, 2)=ixMmax^D+(interpolation_order-1)
    ixS_p_max^D(:, 3)=ixMmax^D+(interpolation_order-1)
    
    ixS_p_min^D(-1,1)=1
    ixS_p_max^D(-1,1)=ixMmin^D-1+nxCo^D+nghostcellsCo+(interpolation_order-1)
    ixS_p_min^D( 1,2)=ixMmin^D+nxCo^D-nghostcellsCo-(interpolation_order-1)
    ixS_p_max^D( 1,2)=ixGmax^D
    
    ixR_p_min^D(:, 0)=ixCoMmin^D-nghostcellsCo-(interpolation_order-1)
    ixR_p_min^D(:, 1)=ixCoMmin^D-(interpolation_order-1)
    ixR_p_min^D(:, 2)=ixCoMmin^D-nghostcellsCo-(interpolation_order-1)
    ixR_p_min^D(:, 3)=ixCoMmax^D+1-(interpolation_order-1)
    ixR_p_max^D(:, 0)=nghostcells+(interpolation_order-1)
    ixR_p_max^D(:, 1)=ixCoMmax^D+nghostcellsCo+(interpolation_order-1)
    ixR_p_max^D(:, 2)=ixCoMmax^D+(interpolation_order-1)
    ixR_p_max^D(:, 3)=ixCoMmax^D+nghostcellsCo+(interpolation_order-1)
    
    ixR_p_min^D(-1,1)=1
    ixR_p_max^D(-1,1)=ixCoMmax^D+nghostcellsCo+(interpolation_order-1)
    ixR_p_min^D( 1,2)=ixCoMmin^D-nghostcellsCo-(interpolation_order-1)
    ixR_p_max^D( 1,2)=ixCoGmax^D
    \}

  end subroutine init_bc

  subroutine create_bc_mpi_datatype(nwstart,nwbc) 
    use mod_global_parameters 

    integer, intent(in) :: nwstart, nwbc
    integer :: i^D, ic^D, inc^D, iib^D

    {do i^DB=-1,1\}
      {do iib^DB=-1,1\}
         if (i^D==0|.and.) cycle
         call get_bc_comm_type(type_send_srl(iib^D,i^D),ixS_srl_^L(iib^D,i^D),ixG^LL,nwstart,nwbc)
         call get_bc_comm_type(type_recv_srl(iib^D,i^D),ixR_srl_^L(iib^D,i^D),ixG^LL,nwstart,nwbc)
         call get_bc_comm_type(type_send_r(iib^D,i^D),  ixS_r_^L(iib^D,i^D),ixCoG^L,nwstart,nwbc)
         {do ic^DB=1+int((1-i^DB)/2),2-int((1+i^DB)/2)
            inc^DB=2*i^DB+ic^DB\}
            call get_bc_comm_type(type_recv_r(iib^D,inc^D),ixR_r_^L(iib^D,inc^D), ixG^LL,nwstart,nwbc)
            call get_bc_comm_type(type_send_p(iib^D,inc^D),ixS_p_^L(iib^D,inc^D), ixG^LL,nwstart,nwbc)
            call get_bc_comm_type(type_recv_p(iib^D,inc^D),ixR_p_^L(iib^D,inc^D),ixCoG^L,nwstart,nwbc)
         {end do\}
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
    start(ndim+1)=nwstart
    
    call MPI_TYPE_CREATE_SUBARRAY(ndim+1,fullsize,subsize,start,MPI_ORDER_FORTRAN, &
                                  MPI_DOUBLE_PRECISION,comm_type,ierrmpi)
    call MPI_TYPE_COMMIT(comm_type,ierrmpi)
    
  end subroutine get_bc_comm_type

  subroutine put_bc_comm_types()
    use mod_global_parameters 
 
    integer :: i^D, ic^D, inc^D, iib^D

    {do i^DB=-1,1\}
       {do iib^DB=-1,1\}
           if (i^D==0|.and.) cycle
           call MPI_TYPE_FREE(type_send_srl(iib^D,i^D),ierrmpi)
           call MPI_TYPE_FREE(type_recv_srl(iib^D,i^D),ierrmpi)
           if (levmin==levmax) cycle
           call MPI_TYPE_FREE(type_send_r(iib^D,i^D),ierrmpi)
           {do ic^DB=1+int((1-i^DB)/2),2-int((1+i^DB)/2)
              inc^DB=2*i^DB+ic^DB\}
              call MPI_TYPE_FREE(type_recv_r(iib^D,inc^D),ierrmpi)
              call MPI_TYPE_FREE(type_send_p(iib^D,inc^D),ierrmpi)
              call MPI_TYPE_FREE(type_recv_p(iib^D,inc^D),ierrmpi)
           {end do\}
       {end do\}
    {end do\}
  
  end subroutine put_bc_comm_types

  subroutine getbc(time,qdt,nwstart,nwbc)
    use mod_global_parameters
    use mod_physics
    
    double precision, intent(in)      :: time, qdt
    integer, intent(in)               :: nwstart ! Fill from nw = nwstart+1
    integer, intent(in)               :: nwbc    ! Number of variables to fill
    
    integer :: my_neighbor_type, ipole, idims, iside
    integer :: iigrid, igrid, ineighbor, ipe_neighbor
    integer :: nrecvs, nsends, isizes
    integer :: ixG^L, ixR^L, ixS^L, ixB^L, ixI^L, k^L
    integer :: i^D, n_i^D, ic^D, inc^D, n_inc^D, iib^D
    ! store physical boundary indicating index
    integer :: idphyb(max_blocks,ndim),bindex(ndim)
    integer :: isend_buf(npwbuf), ipwbuf, nghostcellsco,iB
    logical  :: isphysbound
    type(wbuffer) :: pwbuf(npwbuf)

    double precision :: time_bcin
    ! Stretching grid parameters for coarsened block of the current block
    double precision :: logGl,qstl

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
         block=>pw(igrid)
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
             if(idims > 1 .and. neighbor_type(-1,0,igrid)==1) ixBmin1=ixGmin1
             if(idims > 1 .and. neighbor_type( 1,0,igrid)==1) ixBmax1=ixGmax1}
            {^IFTHREED
             if(idims > 1 .and. neighbor_type(-1,0,0,igrid)==1) ixBmin1=ixGmin1
             if(idims > 1 .and. neighbor_type( 1,0,0,igrid)==1) ixBmax1=ixGmax1
             if(idims > 2 .and. neighbor_type(0,-1,0,igrid)==1) ixBmin2=ixGmin2
             if(idims > 2 .and. neighbor_type(0, 1,0,igrid)==1) ixBmax2=ixGmax2}
            do iside=1,2
               i^D=kr(^D,idims)*(2*iside-3);
               if (aperiodB(idims)) then 
                  call physbound(i^D,igrid,isphysbound)
                  if (neighbor_type(i^D,igrid)/=1 .and. .not. isphysbound) cycle
               else 
                  if (neighbor_type(i^D,igrid)/=1) cycle
               end if
               call bc_phys(iside,idims,time,qdt,pw(igrid)%wb,pw(igrid)%x,ixG^L,ixB^L)
            end do
         end do
      end do
    end if
    
    ! default : no singular axis
    ipole=0
    
    irecv=0
    isend=0
    isend_buf=0
    ipwbuf=1
    ! total number of times to call MPI_IRECV in each processor between sibling blocks or from finer neighbors
    nrecvs=nrecv_bc_srl+nrecv_bc_r
    ! total number of times to call MPI_ISEND in each processor between sibling blocks or to coarser neighors
    nsends=nsend_bc_srl+nsend_bc_r
    if (nrecvs>0) then
       allocate(recvstatus(MPI_STATUS_SIZE,nrecvs),recvrequest(nrecvs))
       recvrequest=MPI_REQUEST_NULL
    end if
    if (nsends>0) then
       allocate(sendstatus(MPI_STATUS_SIZE,nsends),sendrequest(nsends))
       sendrequest=MPI_REQUEST_NULL
    end if
    
    ! receiving ghost-cell values from sibling blocks and finer neighbors
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       saveigrid=igrid
       call identifyphysbound(igrid,isphysbound,iib^D)   
       ^D&idphyb(igrid,^D)=iib^D;
       {do i^DB=-1,1\}
          if (i^D==0|.and.) cycle
          my_neighbor_type=neighbor_type(i^D,igrid)
          select case (my_neighbor_type)
          case (3)
             call bc_recv_srl
          case (4)
             call bc_recv_restrict
          end select
       {end do\}
    end do
    
    ! sending ghost-cell values to sibling blocks and coarser neighbors
    nghostcellsco=ceiling(nghostcells*0.5d0)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       saveigrid=igrid
       block=>pw(igrid)
       call identifyphysbound(igrid,isphysbound,iib^D)   
       if (any(neighbor_type(:^D&,igrid)==2)) then
          ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
    {#IFDEF EVOLVINGBOUNDARY
          if(isphysbound) then
            ! coarsen finer ghost cells at physical boundaries
            ixCoMmin^D=ixCoGmin^D+nghostcellsco;
            ixCoMmax^D=ixCoGmax^D-nghostcellsco;
            ixMmin^D=ixGmin^D+(nghostcellsco-1);
            ixMmax^D=ixGmax^D-(nghostcellsco-1);
          else
            ixCoM^L=ixCoG^L^LSUBnghostcells;
            ixM^L=ixG^L^LSUBnghostcells;
          end if
    }
          call coarsen_grid(pw(igrid)%wb,pw(igrid)%x,ixG^L,ixM^L,pw(igrid)%wcoarse,pw(igrid)%xcoarse,&
                            ixCoG^L,ixCoM^L,igrid,igrid)
       end if
    
       {do i^DB=-1,1\}
          if (i^D==0|.and.) cycle
          if (phi_ > 0) ipole=neighbor_pole(i^D,igrid)
          my_neighbor_type=neighbor_type(i^D,igrid)
          select case (my_neighbor_type)
          case (2)
             call bc_send_restrict
          case (3)
             call bc_send_srl
          end select
       {end do\}
    end do
    
    !if (irecv/=nrecvs) then
    !   call mpistop("number of recvs in phase1 in amr_ghostcells is incorrect")
    !end if
    !if (isend/=nsends) then
    !   call mpistop("number of sends in phase1 in amr_ghostcells is incorrect")
    !end if
    
    if (irecv>0) then
       call MPI_WAITALL(irecv,recvrequest,recvstatus,ierrmpi)
       deallocate(recvstatus,recvrequest)
    end if
    if (isend>0) then
       call MPI_WAITALL(isend,sendrequest,sendstatus,ierrmpi)
       deallocate(sendstatus,sendrequest)
       do ipwbuf=1,npwbuf
          if (isend_buf(ipwbuf)/=0) deallocate(pwbuf(ipwbuf)%w)
       end do
    end if
    
    irecv=0
    isend=0
    isend_buf=0
    ipwbuf=1
    ! total number of times to call MPI_IRECV in each processor from coarser neighbors
    nrecvs=nrecv_bc_p
    ! total number of times to call MPI_ISEND in each processor to finer neighbors
    nsends=nsend_bc_p
    if (nrecvs>0) then
       allocate(recvstatus(MPI_STATUS_SIZE,nrecvs),recvrequest(nrecvs))
       recvrequest=MPI_REQUEST_NULL
    end if
    if (nsends>0) then
       allocate(sendstatus(MPI_STATUS_SIZE,nsends),sendrequest(nsends))
       sendrequest=MPI_REQUEST_NULL
    end if
    
    ! receiving ghost-cell values from coarser neighbors
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       saveigrid=igrid
       ^D&iib^D=idphyb(igrid,^D);
       {do i^DB=-1,1\}
          if (i^D==0|.and.) cycle
          my_neighbor_type=neighbor_type(i^D,igrid)
          if (my_neighbor_type==2) call bc_recv_prolong
       {end do\}
    end do
    ! sending ghost-cell values to finer neighbors 
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       saveigrid=igrid
       block=>pw(igrid)
       ^D&iib^D=idphyb(igrid,^D);
       ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
       if (any(neighbor_type(:^D&,igrid)==4)) then
          {do i^DB=-1,1\}
             if (i^D==0|.and.) cycle
             if (phi_ > 0) ipole=neighbor_pole(i^D,igrid)
             my_neighbor_type=neighbor_type(i^D,igrid)
             if (my_neighbor_type==4) call bc_send_prolong
          {end do\}
       end if
    end do
    
    !if (irecv/=nrecvs) then
    !   call mpistop("number of recvs in phase2 in amr_ghostcells is incorrect")
    !end if
    !if (isend/=nsends) then
    !   call mpistop("number of sends in phase2 in amr_ghostcells is incorrect")
    !end if
    
    if (irecv>0) then
       call MPI_WAITALL(irecv,recvrequest,recvstatus,ierrmpi)
       deallocate(recvstatus,recvrequest)
    end if
    
    ! do prolongation on the ghost-cell values received from coarser neighbors 
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       saveigrid=igrid
       block=>pw(igrid)
       ^D&iib^D=idphyb(igrid,^D);
       ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
       if (any(neighbor_type(:^D&,igrid)==2)) then
          {do i^DB=-1,1\}
             if (i^D==0|.and.) cycle
             my_neighbor_type=neighbor_type(i^D,igrid)
             if (my_neighbor_type==2) call bc_prolong
          {end do\}
       end if
    end do
    
    if (isend>0) then
       call MPI_WAITALL(isend,sendrequest,sendstatus,ierrmpi)
       deallocate(sendstatus,sendrequest)
       do ipwbuf=1,npwbuf
          if (isend_buf(ipwbuf)/=0) deallocate(pwbuf(ipwbuf)%w)
       end do
    end if

    !> @todo Implement generic interface for modifying ghost cells
     ! modify normal component of magnetic field to fix divB=0 
    if(bcphys .and. physics_type=='mhd' .and. ndim>1) call phys_boundary_adjust()
    
    if (nwaux>0) call fix_auxiliary
    
    time_bc=time_bc+(MPI_WTIME()-time_bcin)
    
    contains

      subroutine bc_send_srl

        ineighbor=neighbor(1,i^D,igrid)
        ipe_neighbor=neighbor(2,i^D,igrid)

        if (ipole==0) then
           n_i^D=-i^D;
           if (ipe_neighbor==mype) then
              ixS^L=ixS_srl_^L(iib^D,i^D);
              ixR^L=ixR_srl_^L(iib^D,n_i^D);
              pw(ineighbor)%wb(ixR^S,nwstart+1:nwstart+nwbc)=&
                  pw(igrid)%wb(ixS^S,nwstart+1:nwstart+nwbc)
           else
              isend=isend+1
              itag=(3**^ND+4**^ND)*(ineighbor-1)+{(n_i^D+1)*3**(^D-1)+}
              call MPI_ISEND(pw(igrid)%wb,1,type_send_srl(iib^D,i^D), &
                             ipe_neighbor,itag,icomm,sendrequest(isend),ierrmpi)
           end if
        else
           ixS^L=ixS_srl_^L(iib^D,i^D);
           select case (ipole)
           {case (^D)
              n_i^D=i^D^D%n_i^DD=-i^DD;\}
           end select
           if (ipe_neighbor==mype) then
              ixR^L=ixR_srl_^L(iib^D,n_i^D);
              call pole_copy(pw(ineighbor)%wb,ixG^L,ixR^L,pw(igrid)%wb,ixG^L,ixS^L)
           else
              if (isend_buf(ipwbuf)/=0) then
                 call MPI_WAIT(sendrequest(isend_buf(ipwbuf)), &
                               sendstatus(:,isend_buf(ipwbuf)),ierrmpi)
                 deallocate(pwbuf(ipwbuf)%w)
              end if
              allocate(pwbuf(ipwbuf)%w(ixS^S,nwstart+1:nwstart+nwbc))
              call pole_buf(pwbuf(ipwbuf)%w,ixS^L,ixS^L,pw(igrid)%wb,ixG^L,ixS^L)
              isend=isend+1
              isend_buf(ipwbuf)=isend
              itag=(3**^ND+4**^ND)*(ineighbor-1)+{(n_i^D+1)*3**(^D-1)+}
              isizes={(ixSmax^D-ixSmin^D+1)*}*nwbc
              call MPI_ISEND(pwbuf(ipwbuf)%w,isizes,MPI_DOUBLE_PRECISION, &
                             ipe_neighbor,itag,icomm,sendrequest(isend),ierrmpi)
              ipwbuf=1+modulo(ipwbuf,npwbuf)
           end if
        end if

      end subroutine bc_send_srl

      subroutine bc_send_restrict
        integer :: ii^D

        ic^D=1+modulo(node(pig^D_,igrid)-1,2);
        if ({.not.(i^D==0.or.i^D==2*ic^D-3)|.or.}) return
        if(isphysbound) then
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
             if(idims > 1 .and. neighbor_type(-1,0,igrid)==1) ixBmin1=ixCoGmin1
             if(idims > 1 .and. neighbor_type( 1,0,igrid)==1) ixBmax1=ixCoGmax1}
             {^IFTHREED
             if(idims > 1 .and. neighbor_type(-1,0,0,igrid)==1) ixBmin1=ixCoGmin1
             if(idims > 1 .and. neighbor_type( 1,0,0,igrid)==1) ixBmax1=ixCoGmax1
             if(idims > 2 .and. neighbor_type(0,-1,0,igrid)==1) ixBmin2=ixCoGmin2
             if(idims > 2 .and. neighbor_type(0, 1,0,igrid)==1) ixBmax2=ixCoGmax2}
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
                if (neighbor_type(ii^D,igrid)/=1) cycle
                call bc_phys(iside,idims,time,0.d0,pw(igrid)%wcoarse,&
                       pw(igrid)%xcoarse,ixCoG^L,ixB^L)
             end do
          end do
        end if

        ineighbor=neighbor(1,i^D,igrid)
        ipe_neighbor=neighbor(2,i^D,igrid)

        if (ipole==0) then
           n_inc^D=-2*i^D+ic^D;
           if (ipe_neighbor==mype) then
              ixS^L=ixS_r_^L(iib^D,i^D);
              ixR^L=ixR_r_^L(iib^D,n_inc^D);
              pw(ineighbor)%wb(ixR^S,nwstart+1:nwstart+nwbc)=&
               pw(igrid)%wcoarse(ixS^S,nwstart+1:nwstart+nwbc)
           else
              isend=isend+1
              itag=(3**^ND+4**^ND)*(ineighbor-1)+3**^ND+{n_inc^D*4**(^D-1)+}
              call MPI_ISEND(pw(igrid)%wcoarse,1,type_send_r(iib^D,i^D), &
                             ipe_neighbor,itag,icomm,sendrequest(isend),ierrmpi)
           end if
        else
           ixS^L=ixS_r_^L(iib^D,i^D);
           select case (ipole)
           {case (^D)
              n_inc^D=2*i^D+(3-ic^D)^D%n_inc^DD=-2*i^DD+ic^DD;\}
           end select
           if (ipe_neighbor==mype) then
              ixR^L=ixR_r_^L(iib^D,n_inc^D);
              call pole_copy(pw(ineighbor)%wb,ixG^L,ixR^L,pw(igrid)%wcoarse,ixCoG^L,ixS^L)
           else
              if (isend_buf(ipwbuf)/=0) then
                 call MPI_WAIT(sendrequest(isend_buf(ipwbuf)), &
                               sendstatus(:,isend_buf(ipwbuf)),ierrmpi)
                 deallocate(pwbuf(ipwbuf)%w)
              end if
              allocate(pwbuf(ipwbuf)%w(ixS^S,nwstart+1:nwstart+nwbc))
              call pole_buf(pwbuf(ipwbuf)%w,ixS^L,ixS^L,pw(igrid)%wcoarse,ixCoG^L,ixS^L)
              isend=isend+1
              isend_buf(ipwbuf)=isend
              itag=(3**^ND+4**^ND)*(ineighbor-1)+3**^ND+{n_inc^D*4**(^D-1)+}
              isizes={(ixSmax^D-ixSmin^D+1)*}*nwbc
              call MPI_ISEND(pwbuf(ipwbuf)%w,isizes,MPI_DOUBLE_PRECISION, &
                             ipe_neighbor,itag,icomm,sendrequest(isend),ierrmpi)
              ipwbuf=1+modulo(ipwbuf,npwbuf)
           end if
        end if

      end subroutine bc_send_restrict

      subroutine bc_send_prolong
        integer :: ii^D

        {do ic^DB=1+int((1-i^DB)/2),2-int((1+i^DB)/2)
           inc^DB=2*i^DB+ic^DB\}
           ixS^L=ixS_p_^L(iib^D,inc^D);

           ineighbor=neighbor_child(1,inc^D,igrid)
           ipe_neighbor=neighbor_child(2,inc^D,igrid)
        
           if (ipole==0) then
              n_i^D=-i^D;
              n_inc^D=ic^D+n_i^D;
              if (ipe_neighbor==mype) then
                 ixR^L=ixR_p_^L(iib^D,n_inc^D);
                 pw(ineighbor)%wcoarse(ixR^S,nwstart+1:nwstart+nwbc) &
                    =pw(igrid)%wb(ixS^S,nwstart+1:nwstart+nwbc)
              else
                 isend=isend+1
                 itag=(3**^ND+4**^ND)*(ineighbor-1)+3**^ND+{n_inc^D*4**(^D-1)+}
                 call MPI_ISEND(pw(igrid)%wb,1,type_send_p(iib^D,inc^D), &
                                ipe_neighbor,itag,icomm,sendrequest(isend),ierrmpi)
              end if
           else
              select case (ipole)
              {case (^D)
                 n_inc^D=inc^D^D%n_inc^DD=ic^DD-i^DD;\}
              end select
              if (ipe_neighbor==mype) then
                 ixR^L=ixR_p_^L(iib^D,n_inc^D);
                 call pole_copy(pw(ineighbor)%wcoarse,ixCoG^L,ixR^L,pw(igrid)%wb,ixG^L,ixS^L)
              else
                 if (isend_buf(ipwbuf)/=0) then
                    call MPI_WAIT(sendrequest(isend_buf(ipwbuf)), &
                                  sendstatus(:,isend_buf(ipwbuf)),ierrmpi)
                    deallocate(pwbuf(ipwbuf)%w)
                 end if
                 allocate(pwbuf(ipwbuf)%w(ixS^S,nwstart+1:nwstart+nwbc))
                 call pole_buf(pwbuf(ipwbuf)%w,ixS^L,ixS^L,pw(igrid)%wb,ixG^L,ixS^L)
                 isend=isend+1
                 isend_buf(ipwbuf)=isend
                 itag=(3**^ND+4**^ND)*(ineighbor-1)+3**^ND+{n_inc^D*4**(^D-1)+}
                 isizes={(ixSmax^D-ixSmin^D+1)*}*nwbc
                 call MPI_ISEND(pwbuf(ipwbuf)%w,isizes,MPI_DOUBLE_PRECISION, &
                                ipe_neighbor,itag,icomm,sendrequest(isend),ierrmpi)
                 ipwbuf=1+modulo(ipwbuf,npwbuf)
              end if
           end if
        {end do\}

      end subroutine bc_send_prolong

      subroutine bc_recv_srl

        ipe_neighbor=neighbor(2,i^D,igrid)
        if (ipe_neighbor/=mype) then
           irecv=irecv+1
           itag=(3**^ND+4**^ND)*(igrid-1)+{(i^D+1)*3**(^D-1)+}
           call MPI_IRECV(pw(igrid)%wb,1,type_recv_srl(iib^D,i^D), &
                          ipe_neighbor,itag,icomm,recvrequest(irecv),ierrmpi)
        end if

      end subroutine bc_recv_srl

      subroutine bc_recv_restrict

        {do ic^DB=1+int((1-i^DB)/2),2-int((1+i^DB)/2)
           inc^DB=2*i^DB+ic^DB\}
           ipe_neighbor=neighbor_child(2,inc^D,igrid)
           if (ipe_neighbor/=mype) then
              irecv=irecv+1
              itag=(3**^ND+4**^ND)*(igrid-1)+3**^ND+{inc^D*4**(^D-1)+}
              call MPI_IRECV(pw(igrid)%wb,1,type_recv_r(iib^D,inc^D), &
                             ipe_neighbor,itag,icomm,recvrequest(irecv),ierrmpi)
           end if
        {end do\}

      end subroutine bc_recv_restrict

      subroutine bc_recv_prolong

        ic^D=1+modulo(node(pig^D_,igrid)-1,2);
        if ({.not.(i^D==0.or.i^D==2*ic^D-3)|.or.}) return

        ipe_neighbor=neighbor(2,i^D,igrid)
        if (ipe_neighbor/=mype) then
           irecv=irecv+1
           inc^D=ic^D+i^D;
           itag=(3**^ND+4**^ND)*(igrid-1)+3**^ND+{inc^D*4**(^D-1)+}
           call MPI_IRECV(pw(igrid)%wcoarse,1,type_recv_p(iib^D,inc^D), &
                          ipe_neighbor,itag,icomm,recvrequest(irecv),ierrmpi)  
        end if

      end subroutine bc_recv_prolong

      subroutine bc_prolong
        use mod_physics, only: phys_convert_before_prolong, &
             phys_convert_after_prolong

        integer :: ixFi^L,ixCo^L,ii^D
        double precision :: dxFi^D, dxCo^D, xFimin^D, xComin^D, invdxCo^D

        ixFi^L=ixR_srl_^L(iib^D,i^D);
        dxFi^D=rnode(rpdx^D_,igrid);
        dxCo^D=two*dxFi^D;
        invdxCo^D=1.d0/dxCo^D;

        xFimin^D=rnode(rpxmin^D_,igrid)-dble(nghostcells)*dxFi^D;
        xComin^D=rnode(rpxmin^D_,igrid)-dble(nghostcells)*dxCo^D;
        if(stretched_grid) then
          qst=qsts(node(plevel_,igrid))
          logG=logGs(node(plevel_,igrid))
          qstl=qsts(node(plevel_,igrid)-1)
          logGl=logGs(node(plevel_,igrid)-1)
          xFimin1=rnode(rpxmin1_,igrid)*qst**(-nghostcells)
          xComin1=rnode(rpxmin1_,igrid)*qstl**(-nghostcells)
        end if

        ixComin^D=int((xFimin^D+(dble(ixFimin^D)-half)*dxFi^D-xComin^D)*invdxCo^D)+1-1;
        ixComax^D=int((xFimin^D+(dble(ixFimax^D)-half)*dxFi^D-xComin^D)*invdxCo^D)+1+1;

        call phys_convert_before_prolong(ixCoG^L,ixCo^L,&
             pw(igrid)%wcoarse,pw(igrid)%xcoarse)

        select case (typeghostfill)
        case ("linear")
           call interpolation_linear(ixFi^L,dxFi^D,xFimin^D,dxCo^D,invdxCo^D,xComin^D)
        case ("copy")
           call interpolation_copy(ixFi^L,dxFi^D,xFimin^D,dxCo^D,invdxCo^D,xComin^D)
        case ("unlimit")
           call interpolation_unlimit(ixFi^L,dxFi^D,xFimin^D,dxCo^D,invdxCo^D,xComin^D)
        case default
           write (unitterm,*) "Undefined typeghostfill ",typeghostfill
           call mpistop("")
        end select

        call phys_convert_after_prolong(ixCoG^L,ixCo^L,&
             pw(igrid)%wcoarse,pw(igrid)%xcoarse)

      end subroutine bc_prolong

      subroutine interpolation_linear(ixFi^L,dxFi^D,xFimin^D, &
                                      dxCo^D,invdxCo^D,xComin^D)
        use mod_physics, only: phys_convert_after_prolong
        integer, intent(in) :: ixFi^L
        double precision, intent(in) :: dxFi^D, xFimin^D,dxCo^D, invdxCo^D, xComin^D

        integer :: ixCo^D, jxCo^D, hxCo^D, ixFi^D, ix^D, iw, idims
        double precision :: xCo^D, xFi^D, eta^D
        double precision :: slopeL, slopeR, slopeC, signC, signR
        double precision :: slope(nwstart+1:nwstart+nwbc,ndim)

        {do ixFi^DB = ixFi^LIM^DB
           ! cell-centered coordinates of fine grid point
           xFi^DB=xFimin^DB+(dble(ixFi^DB)-half)*dxFi^DB
        
           ! indices of coarse cell which contains the fine cell
           ixCo^DB=int((xFi^DB-xComin^DB)*invdxCo^DB)+1
        
           ! cell-centered coordinate for coarse cell
           xCo^DB=xComin^DB+(dble(ixCo^DB)-half)*dxCo^DB\}
           if(stretched_grid) then
             xFi1=xFimin1/(one-half*logG)*qst**(ixFi1-1)
             do ixCo1=1,ixCoGmax1
               xCo1=xComin1/(one-half*logGl)*qstl**(ixCo1-1)
               if(dabs(xFi1-xCo1)<half*logGl*xCo1) exit
             end do
           end if
           ! normalized distance between fine/coarse cell center
           ! in coarse cell: ranges from -0.5 to 0.5 in each direction
           ! (origin is coarse cell center)
           if(slab) then
             eta^D=(xFi^D-xCo^D)*invdxCo^D;
           else
             ix^D=2*int((ixFi^D+ixMlo^D)/2)-ixMlo^D;
             {eta^D=(xFi^D-xCo^D)*invdxCo^D &
                   *two*(one-block%dvolume(ixFi^DD) &
                   /sum(block%dvolume(ix^D:ix^D+1^D%ixFi^DD))) \}
             if(stretched_grid) then
               eta1=(xFi1-xCo1)/(logGl*xCo1)*two*(one-block%dvolume(ixFi^D) &
                   /sum(block%dvolume(ix1:ix1+1^%1ixFi^D))) 
             end if
           end if
        
           do idims=1,ndim
              hxCo^D=ixCo^D-kr(^D,idims)\
              jxCo^D=ixCo^D+kr(^D,idims)\
        
              do iw=nwstart+1,nwstart+nwbc
                 slopeL=pw(igrid)%wcoarse(ixCo^D,iw)-pw(igrid)%wcoarse(hxCo^D,iw)
                 slopeR=pw(igrid)%wcoarse(jxCo^D,iw)-pw(igrid)%wcoarse(ixCo^D,iw)
                 slopeC=half*(slopeR+slopeL)
        
                 ! get limited slope
                 signR=sign(one,slopeR)
                 signC=sign(one,slopeC)
                 select case(typeprolonglimit)
                 case('minmod')
                   slope(iw,idims)=signR*max(zero,min(dabs(slopeR), &
                                                     signR*slopeL))
                 case('woodward')
                   slope(iw,idims)=two*signR*max(zero,min(dabs(slopeR), &
                                      signR*slopeL,signR*half*slopeC))
                 case('mcbeta')
                   slope(iw,idims)=signR*max(zero,min(mcbeta*dabs(slopeR), &
                                      mcbeta*signR*slopeL,signR*slopeC))
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
           pw(igrid)%wb(ixFi^D,nwstart+1:nwstart+nwbc)=pw(igrid)%wcoarse(ixCo^D,nwstart+1:&
             nwstart+nwbc)+{(slope(nwstart+1:nwstart+nwbc,^D)*eta^D)+}
        
        {end do\}
        
        call phys_convert_after_prolong(ixG^LL,ixFi^L,pw(igrid)%wb,pw(igrid)%x)
      
      end subroutine interpolation_linear

      subroutine interpolation_copy(ixFi^L,dxFi^D,xFimin^D, &
                                    dxCo^D,invdxCo^D,xComin^D)
        use mod_physics, only: phys_convert_after_prolong
        integer, intent(in) :: ixFi^L
        double precision, intent(in) :: dxFi^D, xFimin^D,dxCo^D, invdxCo^D, xComin^D

        integer :: ixCo^D, ixFi^D
        double precision :: xFi^D

        {do ixFi^DB = ixFi^LIM^DB
           ! cell-centered coordinates of fine grid point
           xFi^DB=xFimin^DB+(dble(ixFi^DB)-half)*dxFi^DB
        
           ! indices of coarse cell which contains the fine cell
           ixCo^DB=int((xFi^DB-xComin^DB)*invdxCo^DB)+1\}
        
           ! Copy from coarse cell
           pw(igrid)%wb(ixFi^D,nwstart+1:nwstart+nwbc)=pw(igrid)%wcoarse(ixCo^D,nwstart+1:nwstart+nwbc)
        
        {end do\}
        
        call phys_convert_after_prolong(ixG^LL,ixFi^L,pw(igrid)%wb,pw(igrid)%x)
      
      end subroutine interpolation_copy

      subroutine interpolation_unlimit(ixFi^L,dxFi^D,xFimin^D, &
                                       dxCo^D,invdxCo^D,xComin^D)
        use mod_physics, only: phys_convert_after_prolong
        integer, intent(in) :: ixFi^L
        double precision, intent(in) :: dxFi^D, xFimin^D, dxCo^D,invdxCo^D, xComin^D

        integer :: ixCo^D, jxCo^D, hxCo^D, ixFi^D, ix^D, idims
        double precision :: xCo^D, xFi^D, eta^D
        double precision :: slope(nwstart+1:nwstart+nwbc,ndim)

        {do ixFi^DB = ixFi^LIM^DB
           ! cell-centered coordinates of fine grid point
           xFi^DB=xFimin^DB+(dble(ixFi^DB)-half)*dxFi^DB
        
           ! indices of coarse cell which contains the fine cell
           ixCo^DB=int((xFi^DB-xComin^DB)*invdxCo^DB)+1
        
           ! cell-centered coordinate for coarse cell
           xCo^DB=xComin^DB+(dble(ixCo^DB)-half)*dxCo^DB\}
        
           ! normalized distance between fine/coarse cell center
           ! in coarse cell: ranges from -0.5 to 0.5 in each direction
           ! (origin is coarse cell center)
           if (slab) then
              eta^D=(xFi^D-xCo^D)*invdxCo^D;
           else
              ix^D=2*int((ixFi^D+ixMlo^D)/2)-ixMlo^D;
              {eta^D=(xFi^D-xCo^D)*invdxCo^D &
                    *two*(one-block%dvolume(ixFi^DD) &
                    /sum(block%dvolume(ix^D:ix^D+1^D%ixFi^DD))) \}
           end if
        
           do idims=1,ndim
              hxCo^D=ixCo^D-kr(^D,idims)\
              jxCo^D=ixCo^D+kr(^D,idims)\
        
              ! get centered slope
              slope(nwstart+1:nwstart+nwbc,idims)=half*(pw(igrid)%wcoarse(jxCo^D,&
                nwstart+1:nwstart+nwbc)-pw(igrid)%wcoarse(hxCo^D,nwstart+1:nwstart+nwbc))
           end do
        
           ! Interpolate from coarse cell using centered slopes
           pw(igrid)%wb(ixFi^D,nwstart+1:nwstart+nwbc)=pw(igrid)%wcoarse(ixCo^D,nwstart+1:&
             nwstart+nwbc)+{(slope(nwstart+1:nwstart+nwbc,^D)*eta^D)+}
        {end do\}
        
        call phys_convert_after_prolong(ixG^LL,ixFi^L,pw(igrid)%wb,pw(igrid)%x)
      
      end subroutine interpolation_unlimit

      subroutine pole_copy(wrecv,ixIR^L,ixR^L,wsend,ixIS^L,ixS^L)
      
        integer, intent(in) :: ixIR^L,ixR^L,ixIS^L,ixS^L
        double precision :: wrecv(ixIR^S,1:nw), wsend(ixIS^S,1:nw)

        integer :: iw, iB

        select case (ipole)
        {case (^D)
           iside=int((i^D+3)/2)
           iB=2*(^D-1)+iside
           do iw=nwstart+1,nwstart+nwbc
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
      
      end subroutine pole_copy

      subroutine pole_buf(wrecv,ixIR^L,ixR^L,wsend,ixIS^L,ixS^L)
      
        integer, intent(in) :: ixIR^L,ixR^L,ixIS^L,ixS^L
        double precision :: wrecv(ixIR^S,nwstart+1:nwstart+nwbc), wsend(ixIS^S,1:nw)

        integer :: iw, iB

        select case (ipole)
        {case (^D)
           iside=int((i^D+3)/2)
           iB=2*(^D-1)+iside
           do iw=nwstart+1,nwstart+nwbc
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
      
      end subroutine pole_buf

      subroutine fix_auxiliary
        use mod_physics, only: phys_get_aux
      
        integer :: ix^L

        do iigrid=1,igridstail; igrid=igrids(iigrid);
          saveigrid=igrid
          block=>pw(igrid)
          call identifyphysbound(igrid,isphysbound,iib^D)   
             
          {do i^DB=-1,1\}
             if (i^D==0|.and.) cycle
        
             ix^L=ixR_srl_^L(iib^D,i^D);
             call phys_get_aux(.true.,pw(igrid)%wb,pw(igrid)%x,ixG^L,ix^L,"bc")
          {end do\}
        end do
      
      end subroutine fix_auxiliary

  end subroutine getbc

  subroutine physbound(i^D,igrid,isphysbound)
    use mod_forest
    use mod_global_parameters
    
    integer, intent(in)  :: i^D, igrid
    logical, intent(out) :: isphysbound
    type(tree_node_ptr)  :: tree
    integer              :: level, ig^D, ign^D

    isphysbound = .false.
    
    tree%node => igrid_to_node(igrid,mype)%node
    level = tree%node%level
    {ig^D = tree%node%ig^D; }
    
    {ign^D = ig^D + i^D; }
    if ({ign^D .gt. ng^D(level) .or. ign^D .lt. 1|.or.}) isphysbound = .true.
  
  end subroutine physbound

  subroutine identifyphysbound(igrid,isphysbound,iib^D)
    use mod_forest
    use mod_global_parameters

    integer, intent(in)  :: igrid
    logical, intent(out) :: isphysbound
    integer              :: i^D,level, ig^D, ign^D, iib^D,ipole

    isphysbound = .false.
    ipole=0

    level = node(plevel_,igrid)
    {ig^D = node(pig^D_,igrid);}
    iib^D=0;
    {do i^DB=-1,1\}
       if (i^D==0|.and.) cycle
       if (phi_ > 0) then
         if (neighbor_pole(i^D,igrid)>0) then
           ipole=neighbor_pole(i^D,igrid)
         end if
       end if
    {end do\}
    {
    do i^D=-1,1
      ign^D=ig^D+i^D
      ! blocks at periodic boundary have neighbors in the physical domain
      ! thus threated at internal blocks with no physical boundary
      if (periodB(^D)) ign^D=1+modulo(ign^D-1,ng^D(level))
      if (ign^D .gt. ng^D(level)) then
         iib^D=1
      else if (ign^D .lt. 1) then
         iib^D=-1
      end if
      if(ipole==^D) iib^D=0
    end do
    \}
    if ({iib^D/=0|.or.}) isphysbound = .true.

  end subroutine identifyphysbound

end module mod_ghostcells_update
