!> update ghost cells of all blocks including physical boundaries 
module mod_getbc
implicit none
integer :: ixM^L, ixCoG^L, ixCoM^L
! index ranges to send (S) to sibling blocks, receive (R) from 
! sibling blocks, send restricted (r) ghost cells to coarser blocks 
integer, dimension(-1:1,-1:1) :: ixS_srl_^L, ixR_srl_^L, ixS_r_^L
! index ranges to receive restriced ghost cells from finer blocks, 
! send prolongated (p) ghost cells to finer blocks, receive prolongated 
! ghost from coarser blocks
integer, dimension(-1:1, 0:3) :: ixR_r_^L, ixS_p_^L, ixR_p_^L
! MPI derived datatype to send and receive subarrays of ghost cells to 
! neighbor blocks in a different processor
integer, dimension(-1:1^D&,-1:1^D&) :: type_send_srl, type_recv_srl, type_send_r
integer, dimension(-1:1^D&,0:3^D&) :: type_recv_r, type_send_p, type_recv_p

end module mod_getbc
!=============================================================================
subroutine init_bc
use mod_getbc
use mod_global_parameters

integer :: dixBCo, interpolation_order
integer :: nx^D, nxCo^D, ixG^L, i^D, ic^D, inc^D, iib^D
!-----------------------------------------------------------------------------
ixG^L=ixG^LL;
ixM^L=ixG^L^LSUBdixB;
ixCoGmin^D=1;
ixCoGmax^D=ixGmax^D/2+dixB;
ixCoM^L=ixCoG^L^LSUBdixB;

nx^D=ixMmax^D-ixMmin^D+1;
nxCo^D=nx^D/2;

select case (typeghostfill)
case ("copy")
   interpolation_order=1
case ("linear","unlimit")
   interpolation_order=2
case default
   write (unitterm,*) "Undefined typeghostfill ",typeghostfill
   call mpistop("")
end select
dixBCo=int((dixB+1)/2)

if (dixBCo+interpolation_order-1>dixB) then
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
ixS_srl_min^D(:, 1)=ixMmax^D+1-dixB
ixS_srl_max^D(:,-1)=ixMmin^D-1+dixB
ixS_srl_max^D(:, 0)=ixMmax^D
ixS_srl_max^D(:, 1)=ixMmax^D

ixS_srl_min^D(-1,0)=1
ixS_srl_min^D( 1,0)=ixMmin^D
ixS_srl_max^D(-1,0)=ixMmax^D
ixS_srl_max^D( 1,0)=ixGmax^D
 
ixR_srl_min^D(:,-1)=1
ixR_srl_min^D(:, 0)=ixMmin^D
ixR_srl_min^D(:, 1)=ixMmax^D+1
ixR_srl_max^D(:,-1)=dixB
ixR_srl_max^D(:, 0)=ixMmax^D
ixR_srl_max^D(:, 1)=ixGmax^D

ixR_srl_min^D(-1,0)=1
ixR_srl_min^D( 1,0)=ixMmin^D
ixR_srl_max^D(-1,0)=ixMmax^D
ixR_srl_max^D( 1,0)=ixGmax^D

ixS_r_min^D(:,-1)=ixCoMmin^D
ixS_r_min^D(:, 0)=ixCoMmin^D
ixS_r_min^D(:, 1)=ixCoMmax^D+1-dixB
ixS_r_max^D(:,-1)=ixCoMmin^D-1+dixB
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
ixR_r_max^D(:, 0)=dixB
ixR_r_max^D(:, 1)=ixMmin^D-1+nxCo^D
ixR_r_max^D(:, 2)=ixMmax^D
ixR_r_max^D(:, 3)=ixGmax^D

ixR_r_min^D(-1,1)=1
ixR_r_max^D(-1,1)=ixMmin^D-1+nxCo^D
ixR_r_min^D( 1,2)=ixMmin^D+nxCo^D
ixR_r_max^D( 1,2)=ixGmax^D

ixS_p_min^D(:, 0)=ixMmin^D-(interpolation_order-1)
ixS_p_min^D(:, 1)=ixMmin^D-(interpolation_order-1)
ixS_p_min^D(:, 2)=ixMmin^D+nxCo^D-dixBCo-(interpolation_order-1)
ixS_p_min^D(:, 3)=ixMmax^D+1-dixBCo-(interpolation_order-1)
ixS_p_max^D(:, 0)=ixMmin^D-1+dixBCo+(interpolation_order-1)
ixS_p_max^D(:, 1)=ixMmin^D-1+nxCo^D+dixBCo+(interpolation_order-1)
ixS_p_max^D(:, 2)=ixMmax^D+(interpolation_order-1)
ixS_p_max^D(:, 3)=ixMmax^D+(interpolation_order-1)

ixS_p_min^D(-1,1)=1
ixS_p_max^D(-1,1)=ixMmin^D-1+nxCo^D+dixBCo+(interpolation_order-1)
ixS_p_min^D( 1,2)=ixMmin^D+nxCo^D-dixBCo-(interpolation_order-1)
ixS_p_max^D( 1,2)=ixGmax^D

ixR_p_min^D(:, 0)=ixCoMmin^D-dixBCo-(interpolation_order-1)
ixR_p_min^D(:, 1)=ixCoMmin^D-(interpolation_order-1)
ixR_p_min^D(:, 2)=ixCoMmin^D-dixBCo-(interpolation_order-1)
ixR_p_min^D(:, 3)=ixCoMmax^D+1-(interpolation_order-1)
ixR_p_max^D(:, 0)=dixB+(interpolation_order-1)
ixR_p_max^D(:, 1)=ixCoMmax^D+dixBCo+(interpolation_order-1)
ixR_p_max^D(:, 2)=ixCoMmax^D+(interpolation_order-1)
ixR_p_max^D(:, 3)=ixCoMmax^D+dixBCo+(interpolation_order-1)

ixR_p_min^D(-1,1)=1
ixR_p_max^D(-1,1)=ixCoMmax^D+dixBCo+(interpolation_order-1)
ixR_p_min^D( 1,2)=ixCoMmin^D-dixBCo-(interpolation_order-1)
ixR_p_max^D( 1,2)=ixCoGmax^D
\}

{do i^DB=-1,1\}
   {do iib^DB=-1,1\}
       if (i^D==0|.and.) cycle
       call get_bc_comm_type(type_send_srl(iib^D,i^D),ixS_srl_^L(iib^D,i^D),ixG^L)
       call get_bc_comm_type(type_recv_srl(iib^D,i^D),ixR_srl_^L(iib^D,i^D),ixG^L)
       call get_bc_comm_type(type_send_r(iib^D,i^D),ixS_r_^L(iib^D,i^D),ixCoG^L)
       {do ic^DB=1+int((1-i^DB)/2),2-int((1+i^DB)/2)
          inc^DB=2*i^DB+ic^DB\}
          call get_bc_comm_type(type_recv_r(iib^D,inc^D),ixR_r_^L(iib^D,inc^D),ixG^L)
          call get_bc_comm_type(type_send_p(iib^D,inc^D),ixS_p_^L(iib^D,inc^D),ixG^L)
          call get_bc_comm_type(type_recv_p(iib^D,inc^D),ixR_p_^L(iib^D,inc^D),ixCoG^L)
       {end do\}
   {end do\}
{end do\}

contains
!=============================================================================
subroutine get_bc_comm_type(comm_type,ix^L,ixG^L)

integer, intent(inout) :: comm_type
integer, intent(in) :: ix^L, ixG^L

integer, dimension(ndim+1) :: size, subsize, start
!-----------------------------------------------------------------------------
^D&size(^D)=ixGmax^D;
size(ndim+1)=nw
^D&subsize(^D)=ixmax^D-ixmin^D+1;
subsize(ndim+1)=nwflux+nwaux
^D&start(^D)=ixmin^D-1;
start(ndim+1)=0

call MPI_TYPE_CREATE_SUBARRAY(ndim+1,size,subsize,start,MPI_ORDER_FORTRAN, &
                              MPI_DOUBLE_PRECISION,comm_type,ierrmpi)
call MPI_TYPE_COMMIT(comm_type,ierrmpi)

end subroutine get_bc_comm_type
end subroutine init_bc
!=============================================================================
subroutine put_bc_comm_types
use mod_getbc
use mod_global_parameters
integer :: i^D, ic^D, inc^D, iib^D
!-----------------------------------------------------------------------------
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
!=============================================================================
subroutine getbc(time,qdt,pwuse,nwstart,nwbc)
use mod_getbc
use mod_global_parameters

double precision, intent(in)               :: time, qdt
integer, intent(in)                        :: nwstart,nwbc
type(walloc), dimension(ngridshi)          :: pwuse

integer :: my_neighbor_type, ipole, idims, iside
integer :: iigrid, igrid, ineighbor, ipe_neighbor
integer :: nrecvs, nsends, isizes
integer :: ixG^L, ixR^L, ixS^L, ixB^L, ixI^L, k^L
integer :: i^D, n_i^D, ic^D, inc^D, n_inc^D, iib^D
! store physical boundary indicating index
integer :: idphyb(ngridshi,ndim),bindex(ndim)
integer :: isend_buf(npwbuf), ipwbuf, dixBco,iB
type(walloc) :: pwbuf(npwbuf)
logical  :: isphysbound

double precision :: time_bcin
{#IFDEF STRETCHGRID
! Stretching grid parameters for coarsened block of the current block
double precision :: logGl,qstl
}
!-----------------------------------------------------------------------------
time_bcin=MPI_WTIME()
ixG^L=ixG^LL;

if (internalboundary) then 
   call getintbc(time,ixG^L,pwuse)
end if
! fill ghost cells in physical boundaries
if(bcphys) then
  do iigrid=1,igridstail; igrid=igrids(iigrid);
     if(.not.phyboundblock(igrid)) cycle
     saveigrid=igrid
     ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
     do idims=1,ndim
        ! to avoid using as yet unknown corner info in more than 1D, we
        ! fill only interior mesh ranges of the ghost cell ranges at first,
        ! and progressively enlarge the ranges to include corners later
        {
         kmin^D=merge(0, 1, idims==^D)
         kmax^D=merge(0, 1, idims==^D)
         ixBmin^D=ixGmin^D+kmin^D*dixB
         ixBmax^D=ixGmax^D-kmax^D*dixB
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
           if(.not.slab)mygeo=>pgeo(igrid)
           if (B0field) then
              myB0_cell => pB0_cell(igrid)
              {^D&myB0_face^D => pB0_face^D(igrid)\}
           end if
           call bc_phys(iside,idims,time,qdt,pwuse(igrid)%w,px(igrid)%x,ixG^L,ixB^L)
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
dixBco=ceiling(dixB*0.5d0)
do iigrid=1,igridstail; igrid=igrids(iigrid);
   saveigrid=igrid
   call identifyphysbound(igrid,isphysbound,iib^D)   
   if (any(neighbor_type(:^D&,igrid)==2)) then
      ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
{#IFDEF EVOLVINGBOUNDARY
      if(isphysbound) then
        ! coarsen finer ghost cells at physical boundaries
        ixCoMmin^D=ixCoGmin^D+dixBco;
        ixCoMmax^D=ixCoGmax^D-dixBco;
        ixMmin^D=ixGmin^D+(dixBco-1);
        ixMmax^D=ixGmax^D-(dixBco-1);
      else
        ixCoM^L=ixCoG^L^LSUBdixB;
        ixM^L=ixG^L^LSUBdixB;
      end if
}
      call coarsen_grid(pwuse(igrid)%w,px(igrid)%x,ixG^L,ixM^L,pwCoarse(igrid)%w,pxCoarse(igrid)%x,&
                        ixCoG^L,ixCoM^L,pgeo(igrid),pgeoCoarse(igrid),coarsenprimitive,.true.)
   end if

   {do i^DB=-1,1\}
      if (i^D==0|.and.) cycle
      {^IFPHI ipole=neighbor_pole(i^D,igrid)}
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
   ^D&iib^D=idphyb(igrid,^D);
   ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
   if (any(neighbor_type(:^D&,igrid)==4)) then
      {do i^DB=-1,1\}
         if (i^D==0|.and.) cycle
         {^IFPHI ipole=neighbor_pole(i^D,igrid)}
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
{^IFMHD
! modify normal component of magnetic field to fix divB=0 
if(bcphys .and. b0_>0) then
  do iigrid=1,igridstail; igrid=igrids(iigrid);
     if(.not.phyboundblock(igrid)) cycle
     saveigrid=igrid
     ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
     do idims=1,ndim
        ! to avoid using as yet unknown corner info in more than 1D, we
        ! fill only interior mesh ranges of the ghost cell ranges at first,
        ! and progressively enlarge the ranges to include corners later
        do iside=1,2
           i^D=kr(^D,idims)*(2*iside-3);
           if (neighbor_type(i^D,igrid)/=1) cycle
           iB=(idims-1)*2+iside
           if(any(typeB(:,iB)=="special")) then
             if(.not.slab)mygeo=>pgeo(igrid)
             select case (idims)
             {case (^D)
                if (iside==2) then
                   ! maximal boundary
                   ixImin^DD=ixGmax^D+1-dixB^D%ixImin^DD=ixGmin^DD;
                   ixImax^DD=ixGmax^DD;
                else
                   ! minimal boundary
                   ixImin^DD=ixGmin^DD;
                   ixImax^DD=ixGmin^D-1+dixB^D%ixImax^DD=ixGmax^DD;
                end if \}
             end select
             call fixdivB_boundary(ixG^L,ixI^L,pwuse(igrid)%w,px(igrid)%x,iB)
           end if
        end do
     end do
  end do
end if
}

if (nwaux>0) call fix_auxiliary

time_bc=time_bc+(MPI_WTIME()-time_bcin)

contains
!=============================================================================
! internal procedures
!=============================================================================
subroutine bc_send_srl
!-----------------------------------------------------------------------------
ineighbor=neighbor(1,i^D,igrid)
ipe_neighbor=neighbor(2,i^D,igrid)


if (ipole==0) then
   n_i^D=-i^D;
   if (ipe_neighbor==mype) then
      ixS^L=ixS_srl_^L(iib^D,i^D);
      ixR^L=ixR_srl_^L(iib^D,n_i^D);
      pwuse(ineighbor)%w(ixR^S,nwstart+1:nwstart+nwbc)=pwuse(igrid)%w(ixS^S,nwstart+1:nwstart+nwbc)
   else
      isend=isend+1
      itag=(3**^ND+4**^ND)*(ineighbor-1)+{(n_i^D+1)*3**(^D-1)+}
      call MPI_ISEND(pwuse(igrid)%w,1,type_send_srl(iib^D,i^D), &
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
      call pole_copy(pwuse(ineighbor),ixR^L,pwuse(igrid),ixS^L)
   else
      if (isend_buf(ipwbuf)/=0) then
         call MPI_WAIT(sendrequest(isend_buf(ipwbuf)), &
                       sendstatus(:,isend_buf(ipwbuf)),ierrmpi)
         deallocate(pwbuf(ipwbuf)%w)
      end if
      allocate(pwbuf(ipwbuf)%w(ixS^S,nwstart+1:nwstart+nwbc))
      call pole_copy(pwbuf(ipwbuf),ixS^L,pwuse(igrid),ixS^L)
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
!=============================================================================
subroutine bc_send_restrict
integer :: ii^D
!-----------------------------------------------------------------------------
ic^D=1+modulo(node(pig^D_,igrid)-1,2);
if ({.not.(i^D==0.or.i^D==2*ic^D-3)|.or.}) return
if(isphysbound) then
  ! filling physical boundary ghost cells of a coarser representative block for
  ! sending swap region with width of dixB to its coarser neighbor
  do idims=1,ndim
     ! to avoid using as yet unknown corner info in more than 1D, we
     ! fill only interior mesh ranges of the ghost cell ranges at first,
     ! and progressively enlarge the ranges to include corners later
     {kmin^D=merge(0, 1, idims==^D)
     kmax^D=merge(0, 1, idims==^D)
     ixBmin^D=ixCoGmin^D+kmin^D*dixB
     ixBmax^D=ixCoGmax^D-kmax^D*dixB\}
     {^IFTWOD
     if(idims > 1 .and. neighbor_type(-1,0,igrid)==1) ixBmin1=ixCoGmin1
     if(idims > 1 .and. neighbor_type( 1,0,igrid)==1) ixBmax1=ixCoGmax1}
     {^IFTHREED
     if(idims > 1 .and. neighbor_type(-1,0,0,igrid)==1) ixBmin1=ixCoGmin1
     if(idims > 1 .and. neighbor_type( 1,0,0,igrid)==1) ixBmax1=ixCoGmax1
     if(idims > 2 .and. neighbor_type(0,-1,0,igrid)==1) ixBmin2=ixCoGmin2
     if(idims > 2 .and. neighbor_type(0, 1,0,igrid)==1) ixBmax2=ixCoGmax2}
     {if(i^D==-1) then
       ixBmin^D=ixCoGmin^D+dixB
       ixBmax^D=ixCoGmin^D+2*dixB-1
     else if(i^D==1) then
       ixBmin^D=ixCoGmax^D-2*dixB+1
       ixBmax^D=ixCoGmax^D-dixB
     end if\}
     do iside=1,2
        ii^D=kr(^D,idims)*(2*iside-3);
        if ({abs(i^D)==1.and.abs(ii^D)==1|.or.}) cycle
        if (neighbor_type(ii^D,igrid)/=1) cycle
        if (.not.slab) mygeo=>pgeoCoarse(igrid)
        call bc_phys(iside,idims,time,0.d0,pwCoarse(igrid)%w,pxCoarse(igrid)%x,ixCoG^L,ixB^L)
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
      pwuse(ineighbor)%w(ixR^S,nwstart+1:nwstart+nwbc)=pwCoarse(igrid)%w(ixS^S,nwstart+1:nwstart+nwbc)
   else
      isend=isend+1
      itag=(3**^ND+4**^ND)*(ineighbor-1)+3**^ND+{n_inc^D*4**(^D-1)+}
      call MPI_ISEND(pwCoarse(igrid)%w,1,type_send_r(iib^D,i^D), &
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
      call pole_copy(pwuse(ineighbor),ixR^L,pwCoarse(igrid),ixS^L)
   else
      if (isend_buf(ipwbuf)/=0) then
         call MPI_WAIT(sendrequest(isend_buf(ipwbuf)), &
                       sendstatus(:,isend_buf(ipwbuf)),ierrmpi)
         deallocate(pwbuf(ipwbuf)%w)
      end if
      allocate(pwbuf(ipwbuf)%w(ixS^S,nwstart+1:nwstart+nwbc))
      call pole_copy(pwbuf(ipwbuf),ixS^L,pwCoarse(igrid),ixS^L)
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
!=============================================================================
subroutine bc_send_prolong
integer :: ii^D
!-----------------------------------------------------------------------------
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
         pwCoarse(ineighbor)%w(ixR^S,nwstart+1:nwstart+nwbc) &
            =pwuse(igrid)%w(ixS^S,nwstart+1:nwstart+nwbc)
      else
         isend=isend+1
         itag=(3**^ND+4**^ND)*(ineighbor-1)+3**^ND+{n_inc^D*4**(^D-1)+}
         call MPI_ISEND(pwuse(igrid)%w,1,type_send_p(iib^D,inc^D), &
                        ipe_neighbor,itag,icomm,sendrequest(isend),ierrmpi)
      end if
   else
      select case (ipole)
      {case (^D)
         n_inc^D=inc^D^D%n_inc^DD=ic^DD-i^DD;\}
      end select
      if (ipe_neighbor==mype) then
         ixR^L=ixR_p_^L(iib^D,n_inc^D);
         call pole_copy(pwCoarse(ineighbor),ixR^L,pwuse(igrid),ixS^L)
      else
         if (isend_buf(ipwbuf)/=0) then
            call MPI_WAIT(sendrequest(isend_buf(ipwbuf)), &
                          sendstatus(:,isend_buf(ipwbuf)),ierrmpi)
            deallocate(pwbuf(ipwbuf)%w)
         end if
         allocate(pwbuf(ipwbuf)%w(ixS^S,nwstart+1:nwstart+nwbc))
         call pole_copy(pwbuf(ipwbuf),ixS^L,pwuse(igrid),ixS^L)
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
!=============================================================================
subroutine bc_recv_srl
!-----------------------------------------------------------------------------
ipe_neighbor=neighbor(2,i^D,igrid)
if (ipe_neighbor/=mype) then
   irecv=irecv+1
   itag=(3**^ND+4**^ND)*(igrid-1)+{(i^D+1)*3**(^D-1)+}
   call MPI_IRECV(pwuse(igrid)%w,1,type_recv_srl(iib^D,i^D), &
                  ipe_neighbor,itag,icomm,recvrequest(irecv),ierrmpi)
end if

end subroutine bc_recv_srl
!=============================================================================
subroutine bc_recv_restrict
!-----------------------------------------------------------------------------
{do ic^DB=1+int((1-i^DB)/2),2-int((1+i^DB)/2)
   inc^DB=2*i^DB+ic^DB\}
   ipe_neighbor=neighbor_child(2,inc^D,igrid)
   if (ipe_neighbor/=mype) then
      irecv=irecv+1
      itag=(3**^ND+4**^ND)*(igrid-1)+3**^ND+{inc^D*4**(^D-1)+}
      call MPI_IRECV(pwuse(igrid)%w,1,type_recv_r(iib^D,inc^D), &
                     ipe_neighbor,itag,icomm,recvrequest(irecv),ierrmpi)
   end if
{end do\}

end subroutine bc_recv_restrict
!=============================================================================
subroutine bc_recv_prolong
!-----------------------------------------------------------------------------
ic^D=1+modulo(node(pig^D_,igrid)-1,2);
if ({.not.(i^D==0.or.i^D==2*ic^D-3)|.or.}) return

ipe_neighbor=neighbor(2,i^D,igrid)
if (ipe_neighbor/=mype) then
   irecv=irecv+1
   inc^D=ic^D+i^D;
   itag=(3**^ND+4**^ND)*(igrid-1)+3**^ND+{inc^D*4**(^D-1)+}
   call MPI_IRECV(pwCoarse(igrid)%w,1,type_recv_p(iib^D,inc^D), &
                  ipe_neighbor,itag,icomm,recvrequest(irecv),ierrmpi)  
end if

end subroutine bc_recv_prolong
!=============================================================================
subroutine bc_prolong

integer :: ixFi^L,ixCo^L,ii^D
double precision :: dxFi^D, dxCo^D, xFimin^D, xComin^D, invdxCo^D
!-----------------------------------------------------------------------------
ixFi^L=ixR_srl_^L(iib^D,i^D);

dxFi^D=rnode(rpdx^D_,igrid);
dxCo^D=two*dxFi^D;
invdxCo^D=1.d0/dxCo^D;

xFimin^D=rnode(rpxmin^D_,igrid)-dble(dixB)*dxFi^D;
xComin^D=rnode(rpxmin^D_,igrid)-dble(dixB)*dxCo^D;
{#IFDEF STRETCHGRID
qst=qsts(node(plevel_,igrid))
logG=logGs(node(plevel_,igrid))
qstl=qsts(node(plevel_,igrid)-1)
logGl=logGs(node(plevel_,igrid)-1)
xFimin1=rnode(rpxmin1_,igrid)*qst**(-dixB)
xComin1=rnode(rpxmin1_,igrid)*qstl**(-dixB)
}

ixComin^D=int((xFimin^D+(dble(ixFimin^D)-half)*dxFi^D-xComin^D)*invdxCo^D)+1-1;
ixComax^D=int((xFimin^D+(dble(ixFimax^D)-half)*dxFi^D-xComin^D)*invdxCo^D)+1+1;

if (amrentropy) then
   call e_to_rhos(ixCoG^L,ixCo^L,pwCoarse(igrid)%w,pxCoarse(igrid)%x)
else if (prolongprimitive) then
   call primitive(ixCoG^L,ixCo^L,pwCoarse(igrid)%w,pxCoarse(igrid)%x)
end if

select case (typeghostfill)
case ("linear")
   call interpolation_linear(pwuse(igrid),ixFi^L,dxFi^D,xFimin^D, &
                           pwCoarse(igrid),dxCo^D,invdxCo^D,xComin^D)
case ("copy")
   call interpolation_copy(pwuse(igrid),ixFi^L,dxFi^D,xFimin^D, &
                           pwCoarse(igrid),dxCo^D,invdxCo^D,xComin^D)
case ("unlimit")
   call interpolation_unlimit(pwuse(igrid),ixFi^L,dxFi^D,xFimin^D, &
                           pwCoarse(igrid),dxCo^D,invdxCo^D,xComin^D)
case default
   write (unitterm,*) "Undefined typeghostfill ",typeghostfill
   call mpistop("")
end select

if (amrentropy) then
    call rhos_to_e(ixCoG^L,ixCo^L,pwCoarse(igrid)%w,pxCoarse(igrid)%x)
else if (prolongprimitive) then
    call conserve(ixCoG^L,ixCo^L,pwCoarse(igrid)%w,pxCoarse(igrid)%x,patchfalse)
end if

end subroutine bc_prolong
!=============================================================================
subroutine interpolation_linear(pwFi,ixFi^L,dxFi^D,xFimin^D, &
                                pwCo,dxCo^D,invdxCo^D,xComin^D)

integer, intent(in) :: ixFi^L
double precision, intent(in) :: dxFi^D, xFimin^D,dxCo^D, invdxCo^D, xComin^D
type(walloc) :: pwCo, pwFi

integer :: ixCo^D, jxCo^D, hxCo^D, ixFi^D, ix^D, iw, idims
double precision :: xCo^D, xFi^D, eta^D
double precision :: slopeL, slopeR, slopeC, signC, signR
double precision :: slope(nwstart+1:nwstart+nwbc,ndim)
!-----------------------------------------------------------------------------
{do ixFi^DB = ixFi^LIM^DB
   ! cell-centered coordinates of fine grid point
   xFi^DB=xFimin^DB+(dble(ixFi^DB)-half)*dxFi^DB

   ! indices of coarse cell which contains the fine cell
   ixCo^DB=int((xFi^DB-xComin^DB)*invdxCo^DB)+1

   ! cell-centered coordinate for coarse cell
   xCo^DB=xComin^DB+(dble(ixCo^DB)-half)*dxCo^DB\}
{#IFDEF STRETCHGRID
   xFi1=xFimin1/(one-half*logG)*qst**(ixFi1-1)
   do ixCo1=1,ixCoGmax1
     xCo1=xComin1/(one-half*logGl)*qstl**(ixCo1-1)
     if(dabs(xFi1-xCo1)<half*logGl*xCo1) exit
   end do
}
   ! normalized distance between fine/coarse cell center
   ! in coarse cell: ranges from -0.5 to 0.5 in each direction
   ! (origin is coarse cell center)
   if (slab) then
      eta^D=(xFi^D-xCo^D)*invdxCo^D;
   else
      ix^D=2*int((ixFi^D+ixMlo^D)/2)-ixMlo^D;
      {eta^D=(xFi^D-xCo^D)*invdxCo^D &
            *two*(one-pgeo(igrid)%dvolume(ixFi^DD) &
            /sum(pgeo(igrid)%dvolume(ix^D:ix^D+1^D%ixFi^DD))) \}
{#IFDEF STRETCHGRID
      eta1=(xFi1-xCo1)/(logGl*xCo1)*two*(one-pgeo(igrid)%dvolume(ixFi^D) &
            /sum(pgeo(igrid)%dvolume(ix1:ix1+1^%1ixFi^D))) 
}
   end if

   do idims=1,ndim
      hxCo^D=ixCo^D-kr(^D,idims)\
      jxCo^D=ixCo^D+kr(^D,idims)\

      do iw=nwstart+1,nwstart+nwbc
         slopeL=pwCo%w(ixCo^D,iw)-pwCo%w(hxCo^D,iw)
         slopeR=pwCo%w(jxCo^D,iw)-pwCo%w(ixCo^D,iw)
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
   pwFi%w(ixFi^D,nwstart+1:nwstart+nwbc)=pwCo%w(ixCo^D,nwstart+1:nwstart+nwbc)+{(slope(nwstart+1:nwstart+nwbc,^D)*eta^D)+}

{end do\}

if (amrentropy) then
   call rhos_to_e(ixG^LL,ixFi^L,pwFi%w,px(igrid)%x)
else if (prolongprimitive) then
   call conserve(ixG^LL,ixFi^L,pwFi%w,px(igrid)%x,patchfalse)
end if

end subroutine interpolation_linear
!=============================================================================
subroutine interpolation_copy(pwFi,ixFi^L,dxFi^D,xFimin^D, &
                              pwCo,dxCo^D,invdxCo^D,xComin^D)

integer, intent(in) :: ixFi^L
double precision, intent(in) :: dxFi^D, xFimin^D,dxCo^D, invdxCo^D, xComin^D
type(walloc) :: pwCo, pwFi

integer :: ixCo^D, ixFi^D
double precision :: xFi^D
!-----------------------------------------------------------------------------
{do ixFi^DB = ixFi^LIM^DB
   ! cell-centered coordinates of fine grid point
   xFi^DB=xFimin^DB+(dble(ixFi^DB)-half)*dxFi^DB

   ! indices of coarse cell which contains the fine cell
   ixCo^DB=int((xFi^DB-xComin^DB)*invdxCo^DB)+1\}

   ! Copy from coarse cell
   pwFi%w(ixFi^D,nwstart+1:nwstart+nwbc)=pwCo%w(ixCo^D,nwstart+1:nwstart+nwbc)

{end do\}

if (amrentropy) then
   call rhos_to_e(ixG^LL,ixFi^L,pwFi%w,px(igrid)%x)
else if (prolongprimitive) then
   call conserve(ixG^LL,ixFi^L,pwFi%w,px(igrid)%x,patchfalse)
end if

end subroutine interpolation_copy
!=============================================================================
subroutine interpolation_unlimit(pwFi,ixFi^L,dxFi^D,xFimin^D, &
                                 pwCo,dxCo^D,invdxCo^D,xComin^D)

integer, intent(in) :: ixFi^L
double precision, intent(in) :: dxFi^D, xFimin^D, dxCo^D,invdxCo^D, xComin^D
type(walloc) :: pwCo, pwFi

integer :: ixCo^D, jxCo^D, hxCo^D, ixFi^D, ix^D, idims
double precision :: xCo^D, xFi^D, eta^D
double precision :: slope(nwstart+1:nwstart+nwbc,ndim)
!-----------------------------------------------------------------------------
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
            *two*(one-pgeo(igrid)%dvolume(ixFi^DD) &
            /sum(pgeo(igrid)%dvolume(ix^D:ix^D+1^D%ixFi^DD))) \}
   end if

   do idims=1,ndim
      hxCo^D=ixCo^D-kr(^D,idims)\
      jxCo^D=ixCo^D+kr(^D,idims)\

      ! get centered slope
      slope(nwstart+1:nwstart+nwbc,idims)=half*(pwCo%w(jxCo^D,nwstart+1:nwstart+nwbc)-pwCo%w(hxCo^D,nwstart+1:nwstart+nwbc))
   end do

   ! Interpolate from coarse cell using centered slopes
   pwFi%w(ixFi^D,nwstart+1:nwstart+nwbc)=pwCo%w(ixCo^D,nwstart+1:nwstart+nwbc)+{(slope(nwstart+1:nwstart+nwbc,^D)*eta^D)+}
{end do\}

if (amrentropy) then
   call rhos_to_e(ixG^LL,ixFi^L,pwFi%w,px(igrid)%x)
else if (prolongprimitive) then
   call conserve(ixG^LL,ixFi^L,pwFi%w,px(igrid)%x,patchfalse)
end if

end subroutine interpolation_unlimit
!=============================================================================
subroutine pole_copy(pwrecv,ixR^L,pwsend,ixS^L)

integer, intent(in) :: ixR^L, ixS^L
type(walloc) :: pwrecv, pwsend

integer :: iw, iB
!-----------------------------------------------------------------------------
select case (ipole)
{case (^D)
   iside=int((i^D+3)/2)
   iB=2*(^D-1)+iside
   do iw=nwstart+1,nwstart+nwbc
      select case (typeB(iw,iB))
      case ("symm")
         pwrecv%w(ixR^S,iw) = pwsend%w(ixSmax^D:ixSmin^D:-1^D%ixS^S,iw)
      case ("asymm")
         pwrecv%w(ixR^S,iw) =-pwsend%w(ixSmax^D:ixSmin^D:-1^D%ixS^S,iw)
      case default
         call mpistop("Boundary condition at pole should be symm or asymm")
      end select
   end do \}
end select

end subroutine pole_copy
!=============================================================================
subroutine fix_auxiliary

integer :: ix^L
!-----------------------------------------------------------------------------
do iigrid=1,igridstail; igrid=igrids(iigrid);
   saveigrid=igrid
   call identifyphysbound(igrid,isphysbound,iib^D)   
      
   {do i^DB=-1,1\}
      if (i^D==0|.and.) cycle

      ix^L=ixR_srl_^L(iib^D,i^D);
      if(.not.slab)mygeo=>pgeo(igrid)
      call getaux(.true.,pwuse(igrid)%w,px(igrid)%x,ixG^L,ix^L,"bc")
   {end do\}
end do

end subroutine fix_auxiliary
!=============================================================================
! end of internal procedures
!=============================================================================
end subroutine getbc
!=============================================================================
subroutine physbound(i^D,igrid,isphysbound)
use mod_forest
use mod_global_parameters

integer, intent(in)  :: i^D, igrid
logical, intent(out) :: isphysbound
type(tree_node_ptr)  :: tree
integer              :: level, ig^D, ign^D
!-----------------------------------------------------------------------------
isphysbound = .false.

tree%node => igrid_to_node(igrid,mype)%node
level = tree%node%level
{ig^D = tree%node%ig^D; }

{ign^D = ig^D + i^D; }
if ({ign^D .gt. ng^D(level) .or. ign^D .lt. 1|.or.}) isphysbound = .true.

end subroutine physbound
!=============================================================================
subroutine identifyphysbound(igrid,isphysbound,iib^D)
use mod_forest
use mod_global_parameters

integer, intent(in)  :: igrid
logical, intent(out) :: isphysbound
type(tree_node_ptr)  :: tree
integer              :: i^D,level, ig^D, ign^D, iib^D
!-----------------------------------------------------------------------------
isphysbound = .false.

tree%node => igrid_to_node(igrid,mype)%node
level = tree%node%level
{ig^D = tree%node%ig^D; }
iib^D=0;
{do i^DB=-1,1\}
    if (i^D==0|.and.) cycle
   {ign^D = ig^D + i^D; }
   ! blocks at periodic boundary have neighbors in the physical domain
   ! thus threated at internal blocks with no physical boundary 
   {if (periodB(^D)) ign^D=1+modulo(ign^D-1,ng^D(level))\}
   {
   if (ign^D .gt. ng^D(level)) then
      iib^D=1
      isphysbound = .true.
   else if (ign^D .lt. 1) then
      iib^D=-1
      isphysbound = .true.
   end if
   \}
{end do\}
end subroutine identifyphysbound
!=============================================================================
subroutine fixdivB_boundary(ixG^L,ixO^L,w,x,iB)
use mod_global_parameters

integer, intent(in) :: ixG^L,ixO^L,iB
double precision, intent(inout) :: w(ixG^S,1:nw)
double precision, intent(in) :: x(ixG^S,1:ndim)

double precision :: dx1x2,dx1x3,dx2x1,dx2dx3,dx3x1,dx3x2
integer :: ix^D
!-----------------------------------------------------------------------------
select case(iB)
 case(1)
   ! 2nd order CD for divB=0 to set normal B component better
   call primitive(ixG^L,ixO^L,w,x)
   {^IFTWODMHD
   dx1x2=dxlevel(1)/dxlevel(2)
   do ix2=ixOmin2+1,ixOmax2-1
     do ix1=ixOmax1,ixOmin1,-1
       w(ix1,ix2,b1_)=w(ix1+2,ix2,b1_) &
        +dx1x2*(w(ix1+1,ix2+1,b2_)-w(ix1+1,ix2-1,b2_))
     enddo
   enddo
   }
   {^IFTHREEDMHD
   dx1x2=dxlevel(1)/dxlevel(2)
   dx1x3=dxlevel(1)/dxlevel(3)
   do ix3=ixOmin3+1,ixOmax3-1
     do ix2=ixOmin2+1,ixOmax2-1
       do ix1=ixOmax1,ixOmin1,-1
         w(ix1,ix2,ix3,b1_)=w(ix1+2,ix2,ix3,b1_) &
          +dx1x2*(w(ix1+1,ix2+1,ix3,b2_)-w(ix1+1,ix2-1,ix3,b2_))&
          +dx1x3*(w(ix1+1,ix2,ix3+1,b3_)-w(ix1+1,ix2,ix3-1,b3_))
       enddo
     enddo
   enddo
   }
   call conserve(ixG^L,ixO^L,w,x,patchfalse)
 case(2)
   call primitive(ixG^L,ixO^L,w,x)
   {^IFTWODMHD
   dx1x2=dxlevel(1)/dxlevel(2)
   do ix2=ixOmin2+1,ixOmax2-1
     do ix1=ixOmin1,ixOmax1
       w(ix1,ix2,b1_)=w(ix1-2,ix2,b1_) &
        -dx1x2*(w(ix1-1,ix2+1,b2_)-w(ix1-1,ix2-1,b2_))
     enddo
   enddo
   }
   {^IFTHREEDMHD
   dx1x2=dxlevel(1)/dxlevel(2)
   dx1x3=dxlevel(1)/dxlevel(3)
   do ix3=ixOmin3+1,ixOmax3-1
     do ix2=ixOmin2+1,ixOmax2-1
       do ix1=ixOmin1,ixOmax1
         w(ix1,ix2,ix3,b1_)=w(ix1-2,ix2,ix3,b1_) &
          -dx1x2*(w(ix1-1,ix2+1,ix3,b2_)-w(ix1-1,ix2-1,ix3,b2_))&
          -dx1x3*(w(ix1-1,ix2,ix3+1,b3_)-w(ix1-1,ix2,ix3-1,b3_))
       enddo
     enddo
   enddo
   }
   call conserve(ixG^L,ixO^L,w,x,patchfalse)
 case(3)
   call primitive(ixG^L,ixO^L,w,x)
   {^IFTWODMHD
   dx2x1=dxlevel(2)/dxlevel(1)
   do ix2=ixOmax2,ixOmin2,-1
     do ix1=ixOmin1+1,ixOmax1-1
       w(ix1,ix2,b2_)=w(ix1,ix2+2,b2_) &
        +dx2x1*(w(ix1+1,ix2+1,b1_)-w(ix1-1,ix2+1,b1_))
     enddo
   enddo
   }
   {^IFTHREEDMHD
   dx2x1=dxlevel(2)/dxlevel(1)
   dx2x3=dxlevel(2)/dxlevel(3)
   do ix3=ixOmin3+1,ixOmax3-1
     do ix2=ixOmax2,ixOmin2,-1
       do ix1=ixOmin1+1,ixOmax1-1
         w(ix1,ix2,ix3,b2_)=w(ix1,ix2+2,ix3,b2_) &
          +dx2x1*(w(ix1+1,ix2+1,ix3,b1_)-w(ix1-1,ix2+1,ix3,b1_))&
          +dx2x3*(w(ix1,ix2+1,ix3+1,b3_)-w(ix1,ix2+1,ix3-1,b3_))
       enddo
     enddo
   enddo
   }
   call conserve(ixG^L,ixO^L,w,x,patchfalse)
 case(4)
   call primitive(ixG^L,ixO^L,w,x)
   {^IFTWODMHD
   dx2x1=dxlevel(2)/dxlevel(1)
   do ix2=ixOmin2,ixOmax2
     do ix1=ixOmin1+1,ixOmax1-1
         w(ix1,ix2,b2_)=w(ix1,ix2-2,b2_) &
          -dx2x1*(w(ix1+1,ix2-1,b1_)-w(ix1-1,ix2-1,b1_))
     enddo
   enddo
   }
   {^IFTHREEDMHD
   dx2x1=dxlevel(2)/dxlevel(1)
   dx2x3=dxlevel(2)/dxlevel(3)
   do ix3=ixOmin3+1,ixOmax3-1
     do ix2=ixOmin2,ixOmax2
       do ix1=ixOmin1+1,ixOmax1-1
         w(ix1,ix2,ix3,b2_)=w(ix1,ix2-2,ix3,b2_) &
          -dx2x1*(w(ix1+1,ix2-1,ix3,b1_)-w(ix1-1,ix2-1,ix3,b1_))&
          -dx2x3*(w(ix1,ix2-1,ix3+1,b3_)-w(ix1,ix2-1,ix3-1,b3_))
       enddo
     enddo
   enddo
   }
   call conserve(ixG^L,ixO^L,w,x,patchfalse)
 case(5)
   call primitive(ixG^L,ixO^L,w,x)
   {^IFTHREEDMHD
   dx3x1=dxlevel(3)/dxlevel(1)
   dx3x2=dxlevel(3)/dxlevel(2)
   do ix3=ixOmax3,ixOmin3,-1
     do ix2=ixOmin2+1,ixOmax2-1
       do ix1=ixOmin1+1,ixOmax1-1
         w(ix1,ix2,ix3,b3_)=w(ix1,ix2,ix3+2,b3_) &
         +dx3x1*(w(ix1+1,ix2,ix3+1,b1_)-w(ix1-1,ix2,ix3+1,b1_))&
         +dx3x2*(w(ix1,ix2+1,ix3+1,b2_)-w(ix1,ix2-1,ix3+1,b2_))
       enddo
     enddo
   enddo
   }
   call conserve(ixG^L,ixO^L,w,x,patchfalse)
 case(6)
   call primitive(ixG^L,ixO^L,w,x)
   {^IFTHREEDMHD
   dx3x1=dxlevel(3)/dxlevel(1)
   dx3x2=dxlevel(3)/dxlevel(2)
   do ix3=ixOmin3,ixOmax3
     do ix2=ixOmin2+1,ixOmax2-1
       do ix1=ixOmin1+1,ixOmax1-1
         w(ix1,ix2,ix3,b3_)=w(ix1,ix2,ix3-2,b3_) &
         -dx3x1*(w(ix1+1,ix2,ix3-1,b1_)-w(ix1-1,ix2,ix3-1,b1_))&
         -dx3x2*(w(ix1,ix2+1,ix3-1,b2_)-w(ix1,ix2-1,ix3-1,b2_))
       enddo
     enddo
   enddo
   }
   call conserve(ixG^L,ixO^L,w,x,patchfalse)
 case default
   call mpistop("Special boundary is not defined for this region")
end select

end subroutine fixdivB_boundary
!=============================================================================
