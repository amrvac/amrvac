!=============================================================================
subroutine getbc(time,ixGmin1,ixGmax1,pwuse,pwuseCo,pgeoFi,pgeoCo,richardson)

include 'amrvacdef.f'

double precision, intent(in)               :: time
integer, intent(in)                        :: ixGmin1,ixGmax1
type(walloc), dimension(ngridshi)          :: pwuse, pwuseCo
type(geoalloc), target,dimension(ngridshi) :: pgeoFi, pgeoCo
logical, intent(in)                        :: richardson

integer :: ixMmin1,ixMmax1, ixCoGmin1,ixCoGmax1, ixCoMmin1,ixCoMmax1, idims,&
    iside
integer :: my_neighbor_type, ipole
integer :: iigrid, igrid, ineighbor, ipe_neighbor
integer :: nrecvs, nsends, isizes
integer :: ixRmin1,ixRmax1, ixSmin1,ixSmax1
integer :: ixBmin1,ixBmax1
integer :: kmin1,kmax1
integer :: i1, n_i1, ic1, inc1, n_inc1
integer, dimension(-1:1) :: ixS_srl_min1,ixS_srl_max1, ixR_srl_min1,&
   ixR_srl_max1, ixS_r_min1,ixS_r_max1
integer, dimension(0:3) :: ixR_r_min1,ixR_r_max1, ixS_p_min1,ixS_p_max1,&
    ixR_p_min1,ixR_p_max1, ixS_old_min1,ixS_old_max1, ixR_old_min1,&
   ixR_old_max1
integer, dimension(-1:1) :: type_send_srl, type_recv_srl, type_send_r
integer, dimension(0:3) :: type_recv_r, type_send_p, type_recv_p,&
    type_send_old, type_recv_old
integer :: isend_buf(npwbuf), ipwbuf
type(walloc) :: pwbuf(npwbuf)
logical  :: isphysbound

double precision :: time_bcin
!-----------------------------------------------------------------------------
time_bcin=MPI_WTIME()

call init_bc
if (internalboundary) then 
   call getintbc(time,ixGmin1,ixGmax1,pwuse)
end if

! default : no singular axis
ipole=0

irecv=0
isend=0
isend_buf=0
ipwbuf=1
nrecvs=nrecv_bc_srl+nrecv_bc_r
nsends=nsend_bc_srl+nsend_bc_r
if (richardson) then
   nrecvs=nrecvs+nrecv_bc_p
   nsends=nsends+nsend_bc_p
end if
if (nrecvs>0) then
   allocate(recvstatus(MPI_STATUS_SIZE,nrecvs),recvrequest(nrecvs))
   recvrequest=MPI_REQUEST_NULL
end if
if (nsends>0) then
   allocate(sendstatus(MPI_STATUS_SIZE,nsends),sendrequest(nsends))
   sendrequest=MPI_REQUEST_NULL
end if

do iigrid=1,igridstail; igrid=igrids(iigrid);
      saveigrid=igrid
      
   do i1=-1,1
      if (i1==0) cycle

      my_neighbor_type=neighbor_type(i1,igrid)
      select case (my_neighbor_type)
      case (2)
         if (richardson) call bc_recv_old
      case (3)
         call bc_recv_srl
      case (4)
         call bc_recv_restrict
      end select
   end do
end do

do iigrid=1,igridstail; igrid=igrids(iigrid);
      saveigrid=igrid
   if (any(neighbor_type(:,igrid)==2)) then
      if (richardson) then
         dxlevel(1)=two*rnode(rpdx1_,igrid);
      else
         dxlevel(1)=rnode(rpdx1_,igrid);
      end if
      call coarsen_grid(pwuse(igrid)%w,px(igrid)%x,ixGmin1,ixGmax1,ixMmin1,&
         ixMmax1,pwuseCo(igrid)%w,pxCoarse(igrid)%x, ixCoGmin1,ixCoGmax1,&
         ixCoMmin1,ixCoMmax1,pgeoFi(igrid),pgeoCo(igrid), coarsenprimitive,&
         .true.)
   end if
   do i1=-1,1
      if (i1==0) cycle

      
      my_neighbor_type=neighbor_type(i1,igrid)
      select case (my_neighbor_type)
      case (2)
         call bc_send_restrict
      case (3)
         call bc_send_srl
      case (4)
         if (richardson) call bc_send_old
      end select
   end do
end do

if (irecv/=nrecvs) then
   call mpistop("number of recvs in phase1 in amr_ghostcells is incorrect")
end if
if (isend/=nsends) then
   call mpistop("number of sends in phase1 in amr_ghostcells is incorrect")
end if

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


if (.not.richardson) then
   irecv=0
   isend=0
   isend_buf=0
   ipwbuf=1
   nrecvs=nrecv_bc_p
   nsends=nsend_bc_p
   if (nrecvs>0) then
      allocate(recvstatus(MPI_STATUS_SIZE,nrecvs),recvrequest(nrecvs))
      recvrequest=MPI_REQUEST_NULL
   end if
   if (nsends>0) then
      allocate(sendstatus(MPI_STATUS_SIZE,nsends),sendrequest(nsends))
      sendrequest=MPI_REQUEST_NULL
   end if

   do iigrid=1,igridstail; igrid=igrids(iigrid);
      saveigrid=igrid
      
      do i1=-1,1
         if (i1==0) cycle

         my_neighbor_type=neighbor_type(i1,igrid)
         if (my_neighbor_type==2) call bc_recv_prolong
      end do
   end do
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      saveigrid=igrid
      dxlevel(1)=rnode(rpdx1_,igrid);
     if (any(neighbor_type(:,igrid)==4)) then
      do i1=-1,1
         if (i1==0) cycle

         
         my_neighbor_type=neighbor_type(i1,igrid)
         if (my_neighbor_type==4) call bc_send_prolong
      end do
     end if
   end do


   if (irecv/=nrecvs) then
      call mpistop("number of recvs in phase2 in amr_ghostcells is incorrect")
   end if
   if (isend/=nsends) then
      call mpistop("number of sends in phase2 in amr_ghostcells is incorrect")
   end if

   if (irecv>0) then
      call MPI_WAITALL(irecv,recvrequest,recvstatus,ierrmpi)
      deallocate(recvstatus,recvrequest)
   end if

   do iigrid=1,igridstail; igrid=igrids(iigrid);
      saveigrid=igrid
      dxlevel(1)=rnode(rpdx1_,igrid);
      if (any(neighbor_type(:,igrid)==2)) then
         do i1=-1,1
            if (i1==0) cycle
            my_neighbor_type=neighbor_type(i1,igrid)
            if (my_neighbor_type==2) call bc_prolong
         end do
      end if
   end do

   if (isend>0) then
      call MPI_WAITALL(isend,sendrequest,sendstatus,ierrmpi)
      deallocate(sendstatus,sendrequest)
      do ipwbuf=1,npwbuf
         if (isend_buf(ipwbuf)/=0) deallocate(pwbuf(ipwbuf)%w)
      end do
   end if
end if

if(.not.energyonly) then
do iigrid=1,igridstail; igrid=igrids(iigrid);
   dxlevel(1)=rnode(rpdx1_,igrid);
   do idims=1,ndim
      ! to avoid using as yet unknown corner info in more than 1D, we
      ! fill only interior mesh ranges of the ghost cell ranges at first,
      ! and progressively enlarge the ranges to include corners later
      kmin1=0; kmax1=0;
      
      
      ixBmin1=ixGmin1+kmin1*dixB;
      ixBmax1=ixGmax1-kmax1*dixB;
      do iside=1,2
         i1=kr(1,idims)*(2*iside-3);
         if (aperiodB(idims)) then 
            call physbound(i1,igrid,isphysbound)
            if (neighbor_type(i1,igrid)/=1 .and. .not. isphysbound) cycle
         else 
            if (neighbor_type(i1,igrid)/=1) cycle
         end if
         if (richardson) then
            if(.not.slab)mygeo=>pgeoCo(igrid)
            call bc_phys(iside,idims,time,pwuse(igrid)%w,pxCoarse(igrid)%x,&
               ixGmin1,ixGmax1,ixBmin1,ixBmax1)
         else
            if(.not.slab)mygeo=>pgeoFi(igrid)
            if (B0field) then
               myB0_cell => pB0_cell(igrid)
               myB0_face1 => pB0_face1(igrid)
            end if
            call bc_phys(iside,idims,time,pwuse(igrid)%w,px(igrid)%x,ixGmin1,&
               ixGmax1,ixBmin1,ixBmax1)
         end if
      end do
   end do
end do
end if

if (npe>1) call put_bc_comm_types

if (nwaux>0) call fix_auxiliary

time_bc=time_bc+(MPI_WTIME()-time_bcin)

contains
!=============================================================================
! internal procedures
!=============================================================================
subroutine bc_send_srl
!-----------------------------------------------------------------------------
ineighbor=neighbor(1,i1,igrid)
ipe_neighbor=neighbor(2,i1,igrid)


if (ipole==0) then
   n_i1=-i1;
   if (ipe_neighbor==mype) then
      ixSmin1=ixS_srl_min1(i1);ixSmax1=ixS_srl_max1(i1);
      ixRmin1=ixR_srl_min1(n_i1);ixRmax1=ixR_srl_max1(n_i1);
      pwuse(ineighbor)%w(ixRmin1:ixRmax1,1:nwflux+nwaux)=pwuse(igrid)%w&
         (ixSmin1:ixSmax1,1:nwflux+nwaux)
   else
      isend=isend+1
      itag=(3**1+4**1)*(ineighbor-1)+(n_i1+1)*3**(1-1)
      call MPI_ISEND(pwuse(igrid)%w,1,type_send_srl(i1), ipe_neighbor,itag,&
         icomm,sendrequest(isend),ierrmpi)
   end if
else
   ixSmin1=ixS_srl_min1(i1);ixSmax1=ixS_srl_max1(i1);
   select case (ipole)
   case (1)
      n_i1=i1;
   end select
   if (ipe_neighbor==mype) then
      ixRmin1=ixR_srl_min1(n_i1);ixRmax1=ixR_srl_max1(n_i1);
      call pole_copy(pwuse(ineighbor),ixRmin1,ixRmax1,pwuse(igrid),ixSmin1,&
         ixSmax1)
   else
      if (isend_buf(ipwbuf)/=0) then
         call MPI_WAIT(sendrequest(isend_buf(ipwbuf)), sendstatus(1,&
            isend_buf(ipwbuf)),ierrmpi)
         deallocate(pwbuf(ipwbuf)%w)
      end if
      allocate(pwbuf(ipwbuf)%w(ixSmin1:ixSmax1,1:nwflux+nwaux))
      call pole_copy(pwbuf(ipwbuf),ixSmin1,ixSmax1,pwuse(igrid),ixSmin1,&
         ixSmax1)
      isend=isend+1
      isend_buf(ipwbuf)=isend
      itag=(3**1+4**1)*(ineighbor-1)+(n_i1+1)*3**(1-1)
      isizes=(ixSmax1-ixSmin1+1)*(nwflux+nwaux)
      call MPI_ISEND(pwbuf(ipwbuf)%w,isizes,MPI_DOUBLE_PRECISION,&
          ipe_neighbor,itag,icomm,sendrequest(isend),ierrmpi)
      ipwbuf=1+modulo(ipwbuf,npwbuf)
   end if
end if

end subroutine bc_send_srl
!=============================================================================
subroutine bc_send_restrict
!-----------------------------------------------------------------------------
ic1=1+modulo(node(pig1_,igrid)-1,2);
if (.not.(i1==0.or.i1==2*ic1-3)) return

ineighbor=neighbor(1,i1,igrid)
ipe_neighbor=neighbor(2,i1,igrid)

if (ipole==0) then
   n_inc1=-2*i1+ic1;
   if (ipe_neighbor==mype) then
      ixSmin1=ixS_r_min1(i1);ixSmax1=ixS_r_max1(i1);
      ixRmin1=ixR_r_min1(n_inc1);ixRmax1=ixR_r_max1(n_inc1);
      pwuse(ineighbor)%w(ixRmin1:ixRmax1,1:nwflux+nwaux)=pwuseCo(igrid)%w&
         (ixSmin1:ixSmax1,1:nwflux+nwaux)
   else
      isend=isend+1
      itag=(3**1+4**1)*(ineighbor-1)+3**1+n_inc1*4**(1-1)
      call MPI_ISEND(pwuseCo(igrid)%w,1,type_send_r(i1), ipe_neighbor,itag,&
         icomm,sendrequest(isend),ierrmpi)
   end if
else
   ixSmin1=ixS_r_min1(i1);ixSmax1=ixS_r_max1(i1);
   select case (ipole)
   case (1)
      n_inc1=2*i1+(3-ic1);
   end select
   if (ipe_neighbor==mype) then
      ixRmin1=ixR_r_min1(n_inc1);ixRmax1=ixR_r_max1(n_inc1);
      call pole_copy(pwuse(ineighbor),ixRmin1,ixRmax1,pwuseCo(igrid),ixSmin1,&
         ixSmax1)
   else
      if (isend_buf(ipwbuf)/=0) then
         call MPI_WAIT(sendrequest(isend_buf(ipwbuf)), sendstatus(1,&
            isend_buf(ipwbuf)),ierrmpi)
         deallocate(pwbuf(ipwbuf)%w)
      end if
      allocate(pwbuf(ipwbuf)%w(ixSmin1:ixSmax1,1:nwflux+nwaux))
      call pole_copy(pwbuf(ipwbuf),ixSmin1,ixSmax1,pwuseCo(igrid),ixSmin1,&
         ixSmax1)
      isend=isend+1
      isend_buf(ipwbuf)=isend
      itag=(3**1+4**1)*(ineighbor-1)+3**1+n_inc1*4**(1-1)
      isizes=(ixSmax1-ixSmin1+1)*(nwflux+nwaux)
      call MPI_ISEND(pwbuf(ipwbuf)%w,isizes,MPI_DOUBLE_PRECISION,&
          ipe_neighbor,itag,icomm,sendrequest(isend),ierrmpi)
      ipwbuf=1+modulo(ipwbuf,npwbuf)
   end if
end if

end subroutine bc_send_restrict
!=============================================================================
subroutine bc_send_prolong
integer :: ii1
!-----------------------------------------------------------------------------
do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
   inc1=2*i1+ic1

   ixSmin1=ixS_p_min1(inc1);ixSmax1=ixS_p_max1(inc1);

   do idims=1,ndim
      do iside=1,2
         ii1=kr(1,idims)*(2*iside-3);

         if (neighbor_type(ii1,igrid)/=1) cycle

         if ((  (iside==1.and.idims==1.and.ixSmin1<ixMlo1)) &
            .or.( (iside==2.and.idims==1.and.ixSmax1>ixMhi1)))then
          ixBmin1=merge(ixGmin1,ixSmin1,idims==1);
          ixBmax1=merge(ixGmax1,ixSmax1,idims==1);
         ! to avoid using as yet unknown corner info in more than 1D, we
         ! fill only interior mesh ranges of the ghost cell ranges at first,
         ! and progressively enlarge the ranges to include corners later
          kmin1=0; kmax1=0;
         
         
          ixBmin1=ixBmin1+kmin1;
          ixBmax1=ixBmax1-kmax1;

          if(.not.slab)mygeo=>pgeoFi(igrid)
          if (B0field) then
            myB0_cell => pB0_cell(igrid)
            myB0_face1 => pB0_face1(igrid)
          end if

          call bc_phys(iside,idims,time,pwuse(igrid)%w, px(igrid)%x,ixGmin1,&
             ixGmax1,ixBmin1,ixBmax1)
         end if
      end do
   end do

   ineighbor=neighbor_child(1,inc1,igrid)
   ipe_neighbor=neighbor_child(2,inc1,igrid)

   if (ipole==0) then
      n_i1=-i1;
      n_inc1=ic1+n_i1;
      if (ipe_neighbor==mype) then
         ixRmin1=ixR_p_min1(n_inc1);ixRmax1=ixR_p_max1(n_inc1);
         pwuseCo(ineighbor)%w(ixRmin1:ixRmax1,1:nwflux+nwaux) &
            =pwuse(igrid)%w(ixSmin1:ixSmax1,1:nwflux+nwaux)
      else
         isend=isend+1
         itag=(3**1+4**1)*(ineighbor-1)+3**1+n_inc1*4**(1-1)
         call MPI_ISEND(pwuse(igrid)%w,1,type_send_p(inc1), ipe_neighbor,itag,&
            icomm,sendrequest(isend),ierrmpi)
      end if
   else
      select case (ipole)
      case (1)
         n_inc1=inc1;
      end select
      if (ipe_neighbor==mype) then
         ixRmin1=ixR_p_min1(n_inc1);ixRmax1=ixR_p_max1(n_inc1);
         call pole_copy(pwuseCo(ineighbor),ixRmin1,ixRmax1,pwuse(igrid),&
            ixSmin1,ixSmax1)
      else
         if (isend_buf(ipwbuf)/=0) then
            call MPI_WAIT(sendrequest(isend_buf(ipwbuf)), sendstatus(1,&
               isend_buf(ipwbuf)),ierrmpi)
            deallocate(pwbuf(ipwbuf)%w)
         end if
         allocate(pwbuf(ipwbuf)%w(ixSmin1:ixSmax1,1:nwflux+nwaux))
         call pole_copy(pwbuf(ipwbuf),ixSmin1,ixSmax1,pwuse(igrid),ixSmin1,&
            ixSmax1)
         isend=isend+1
         isend_buf(ipwbuf)=isend
         itag=(3**1+4**1)*(ineighbor-1)+3**1+n_inc1*4**(1-1)
         isizes=(ixSmax1-ixSmin1+1)*(nwflux+nwaux)
         call MPI_ISEND(pwbuf(ipwbuf)%w,isizes,MPI_DOUBLE_PRECISION,&
             ipe_neighbor,itag,icomm,sendrequest(isend),ierrmpi)
         ipwbuf=1+modulo(ipwbuf,npwbuf)
      end if
   end if
end do

end subroutine bc_send_prolong
!=============================================================================
subroutine bc_send_old
!-----------------------------------------------------------------------------
do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
   inc1=2*i1+ic1

   ineighbor=neighbor_child(1,inc1,igrid)
   ipe_neighbor=neighbor_child(2,inc1,igrid)

   if (ipole==0) then
      n_i1=-i1;
      n_inc1=ic1+n_i1;
      if (ipe_neighbor==mype) then
         ixSmin1=ixS_old_min1(inc1);ixSmax1=ixS_old_max1(inc1);
         ixRmin1=ixR_old_min1(n_inc1);ixRmax1=ixR_old_max1(n_inc1);
         pwuse(ineighbor)%w(ixRmin1:ixRmax1,1:nwflux+nwaux) &
            =pwold(igrid)%w(ixSmin1:ixSmax1,1:nwflux+nwaux)
      else
         isend=isend+1
         itag=(3**1+4**1)*(ineighbor-1)+3**1+n_inc1*4**(1-1)
         call MPI_ISEND(pwold(igrid)%w,1,type_send_old(inc1), ipe_neighbor,&
            itag,icomm,sendrequest(isend),ierrmpi)
      end if
   else
      ixSmin1=ixS_old_min1(inc1);ixSmax1=ixS_old_max1(inc1);
      select case (ipole)
      case (1)
         n_inc1=inc1;
      end select
      if (ipe_neighbor==mype) then
         ixRmin1=ixR_old_min1(n_inc1);ixRmax1=ixR_old_max1(n_inc1);
         call pole_copy(pwuse(ineighbor),ixRmin1,ixRmax1,pwold(igrid),ixSmin1,&
            ixSmax1)
      else
         if (isend_buf(ipwbuf)/=0) then
            call MPI_WAIT(sendrequest(isend_buf(ipwbuf)), sendstatus(1,&
               isend_buf(ipwbuf)),ierrmpi)
            deallocate(pwbuf(ipwbuf)%w)
         end if
         allocate(pwbuf(ipwbuf)%w(ixSmin1:ixSmax1,1:nwflux+nwaux))
         call pole_copy(pwbuf(ipwbuf),ixSmin1,ixSmax1,pwold(igrid),ixSmin1,&
            ixSmax1)
         isend=isend+1
         isend_buf(ipwbuf)=isend
         itag=(3**1+4**1)*(ineighbor-1)+3**1+n_inc1*4**(1-1)
         isizes=(ixSmax1-ixSmin1+1)*(nwflux+nwaux)
         call MPI_ISEND(pwbuf(ipwbuf)%w,isizes,MPI_DOUBLE_PRECISION,&
             ipe_neighbor,itag,icomm,sendrequest(isend),ierrmpi)
         ipwbuf=1+modulo(ipwbuf,npwbuf)
      end if
   end if

end do

end subroutine bc_send_old
!=============================================================================
subroutine bc_recv_srl
!-----------------------------------------------------------------------------
ipe_neighbor=neighbor(2,i1,igrid)
if (ipe_neighbor/=mype) then
   irecv=irecv+1
   itag=(3**1+4**1)*(igrid-1)+(i1+1)*3**(1-1)
   call MPI_IRECV(pwuse(igrid)%w,1,type_recv_srl(i1), ipe_neighbor,itag,icomm,&
      recvrequest(irecv),ierrmpi)
end if

end subroutine bc_recv_srl
!=============================================================================
subroutine bc_recv_restrict
!-----------------------------------------------------------------------------
do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
   inc1=2*i1+ic1
   ipe_neighbor=neighbor_child(2,inc1,igrid)
   if (ipe_neighbor/=mype) then
      irecv=irecv+1
      itag=(3**1+4**1)*(igrid-1)+3**1+inc1*4**(1-1)
      call MPI_IRECV(pwuse(igrid)%w,1,type_recv_r(inc1), ipe_neighbor,itag,&
         icomm,recvrequest(irecv),ierrmpi)
   end if
end do

end subroutine bc_recv_restrict
!=============================================================================
subroutine bc_recv_prolong
!-----------------------------------------------------------------------------
ic1=1+modulo(node(pig1_,igrid)-1,2);
if (.not.(i1==0.or.i1==2*ic1-3)) return

ipe_neighbor=neighbor(2,i1,igrid)
if (ipe_neighbor/=mype) then
   irecv=irecv+1
   inc1=ic1+i1;
   itag=(3**1+4**1)*(igrid-1)+3**1+inc1*4**(1-1)
   call MPI_IRECV(pwuseCo(igrid)%w,1,type_recv_p(inc1), ipe_neighbor,itag,&
      icomm,recvrequest(irecv),ierrmpi)  
end if

end subroutine bc_recv_prolong
!=============================================================================
subroutine bc_recv_old
!-----------------------------------------------------------------------------
ic1=1+modulo(node(pig1_,igrid)-1,2);
if (.not.(i1==0.or.i1==2*ic1-3)) return

ipe_neighbor=neighbor(2,i1,igrid)
if (ipe_neighbor/=mype) then
   irecv=irecv+1
   inc1=ic1+i1;
   itag=(3**1+4**1)*(igrid-1)+3**1+inc1*4**(1-1)
   call MPI_IRECV(pwuse(igrid)%w,1,type_recv_old(inc1), ipe_neighbor,itag,&
      icomm,recvrequest(irecv),ierrmpi)
end if

end subroutine bc_recv_old
!=============================================================================
subroutine bc_prolong

integer :: ixFimin1,ixFimax1,ixComin1,ixComax1,ii1
double precision :: dxFi1, dxCo1, xFimin1, xComin1, invdxCo1
!-----------------------------------------------------------------------------
ixFimin1=ixR_srl_min1(i1);ixFimax1=ixR_srl_max1(i1);

dxFi1=rnode(rpdx1_,igrid);
dxCo1=two*dxFi1;
invdxCo1=1.d0/dxCo1;

xFimin1=rnode(rpxmin1_,igrid)-dble(dixB)*dxFi1;
xComin1=rnode(rpxmin1_,igrid)-dble(dixB)*dxCo1;

! moved the physical boundary filling here, to only fill the
! part needed

ixComin1=int((xFimin1+(dble(ixFimin1)-half)*dxFi1-xComin1)*invdxCo1)+1-1;
ixComax1=int((xFimin1+(dble(ixFimax1)-half)*dxFi1-xComin1)*invdxCo1)+1+1;

do idims=1,ndim
   do iside=1,2
      ii1=kr(1,idims)*(2*iside-3);

      if (neighbor_type(ii1,igrid)/=1) cycle

      if  (( (iside==1.and.idims==1.and.ixComin1<ixCoGmin1+dixB) ) &
         .or.( (iside==2.and.idims==1.and.ixComax1>ixCoGmax1-dixB)))then
        ixBmin1=merge(ixCoGmin1,ixComin1,idims==1);
        ixBmax1=merge(ixCoGmax1,ixComax1,idims==1);
        if(.not.slab)mygeo=>pgeoCo(igrid)

        call bc_phys(iside,idims,time,pwuseCo(igrid)%w, pxCoarse(igrid)%x,&
           ixCoGmin1,ixCoGmax1,ixBmin1,ixBmax1)
      end if
   end do
end do

if (amrentropy) then
   call e_to_rhos(ixCoGmin1,ixCoGmax1,ixComin1,ixComax1,pwuseCo(igrid)%w,&
      pxCoarse(igrid)%x)
else if (prolongprimitive) then
   call primitive(ixCoGmin1,ixCoGmax1,ixComin1,ixComax1,pwuseCo(igrid)%w,&
      pxCoarse(igrid)%x)
end if

select case (typeghostfill)
case ("linear")
   call interpolation_linear(pwuse(igrid),ixFimin1,ixFimax1,dxFi1,xFimin1,&
       pwuseCo(igrid),dxCo1,invdxCo1,xComin1)
case ("copy")
   call interpolation_copy(pwuse(igrid),ixFimin1,ixFimax1,dxFi1,xFimin1,&
       pwuseCo(igrid),dxCo1,invdxCo1,xComin1)
case ("unlimit")
   call interpolation_unlimit(pwuse(igrid),ixFimin1,ixFimax1,dxFi1,xFimin1,&
       pwuseCo(igrid),dxCo1,invdxCo1,xComin1)
case default
   write (unitterm,*) "Undefined typeghostfill ",typeghostfill
   call mpistop("")
end select

if (amrentropy) then
    call rhos_to_e(ixCoGmin1,ixCoGmax1,ixComin1,ixComax1,pwuseCo(igrid)%w,&
       pxCoarse(igrid)%x)
else if (prolongprimitive) then
    call conserve(ixCoGmin1,ixCoGmax1,ixComin1,ixComax1,pwuseCo(igrid)%w,&
       pxCoarse(igrid)%x,patchfalse)
end if

end subroutine bc_prolong
!=============================================================================
subroutine interpolation_linear(pwFi,ixFimin1,ixFimax1,dxFi1,xFimin1, pwCo,&
   dxCo1,invdxCo1,xComin1)

integer, intent(in) :: ixFimin1,ixFimax1
double precision, intent(in) :: dxFi1, xFimin1,dxCo1, invdxCo1, xComin1
type(walloc) :: pwCo, pwFi

integer :: ixCo1, jxCo1, hxCo1, ixFi1, ix1, iw, idims
double precision :: xCo1, xFi1, eta1
double precision :: slopeL, slopeR, slopeC, signC, signR, slope(nwflux&
   +nwaux,ndim)
!-----------------------------------------------------------------------------
do ixFi1 = ixFimin1,ixFimax1
   ! cell-centered coordinates of fine grid point
   xFi1=xFimin1+(dble(ixFi1)-half)*dxFi1

   ! indices of coarse cell which contains the fine cell
   ixCo1=int((xFi1-xComin1)*invdxCo1)+1

   ! cell-centered coordinate for coarse cell
   xCo1=xComin1+(dble(ixCo1)-half)*dxCo1

   ! normalized distance between fine/coarse cell center
   ! in coarse cell: ranges from -0.5 to 0.5 in each direction
   ! (origin is coarse cell center)
   if (slab) then
      eta1=(xFi1-xCo1)*invdxCo1;
   else
      ix1=2*int((ixFi1+ixMlo1)/2)-ixMlo1;
      eta1=(xFi1-xCo1)*invdxCo1 *two*(one-pgeoFi(igrid)%dvolume(ixFi1) &
         /sum(pgeoFi(igrid)%dvolume(ix1:ix1+1))) 
   end if

   do idims=1,ndim
      hxCo1=ixCo1-kr(1,idims)
      jxCo1=ixCo1+kr(1,idims)

      do iw=1,nwflux+nwaux
         slopeL=pwCo%w(ixCo1,iw)-pwCo%w(hxCo1,iw)
         slopeR=pwCo%w(jxCo1,iw)-pwCo%w(ixCo1,iw)
         slopeC=half*(slopeR+slopeL)

         ! get limited slope
         signR=sign(one,slopeR)
         signC=sign(one,slopeC)
         select case(typeprolonglimit)
         case('minmod')
           slope(iw,idims)=signR*max(zero,min(dabs(slopeR), signR*slopeL))
         case('woodward')
           slope(iw,idims)=two*signR*max(zero,min(dabs(slopeR), signR&
              *slopeL,signR*half*slopeC))
         case('mcbeta')
           slope(iw,idims)=signR*max(zero,min(mcbeta*dabs(slopeR), mcbeta&
              *signR*slopeL,signR*slopeC))
         case('koren')
           slope(iw,idims)=signR*max(zero,min(two*signR*slopeL, (dabs(slopeR)&
              +two*slopeL*signR)*third,two*dabs(slopeR)))
         case default
           slope(iw,idims)=signC*max(zero,min(dabs(slopeC), signC&
              *slopeL,signC*slopeR))
         end select
      end do
   end do

   ! Interpolate from coarse cell using limited slopes
   pwFi%w(ixFi1,1:nwflux+nwaux)=pwCo%w(ixCo1,1:nwflux+nwaux)+(slope(1:nwflux&
      +nwaux,1)*eta1)

end do

if (amrentropy) then
   call rhos_to_e(ixGlo1,ixGhi1,ixFimin1,ixFimax1,pwFi%w,px(igrid)%x)
else if (prolongprimitive) then
   call conserve(ixGlo1,ixGhi1,ixFimin1,ixFimax1,pwFi%w,px(igrid)%x,&
      patchfalse)
end if

end subroutine interpolation_linear
!=============================================================================
subroutine interpolation_copy(pwFi,ixFimin1,ixFimax1,dxFi1,xFimin1, pwCo,&
   dxCo1,invdxCo1,xComin1)

integer, intent(in) :: ixFimin1,ixFimax1
double precision, intent(in) :: dxFi1, xFimin1,dxCo1, invdxCo1, xComin1
type(walloc) :: pwCo, pwFi

integer :: ixCo1, ixFi1
double precision :: xFi1
!-----------------------------------------------------------------------------
do ixFi1 = ixFimin1,ixFimax1
   ! cell-centered coordinates of fine grid point
   xFi1=xFimin1+(dble(ixFi1)-half)*dxFi1

   ! indices of coarse cell which contains the fine cell
   ixCo1=int((xFi1-xComin1)*invdxCo1)+1

   ! Copy from coarse cell
   pwFi%w(ixFi1,1:nwflux+nwaux)=pwCo%w(ixCo1,1:nwflux+nwaux)

end do

if (amrentropy) then
   call rhos_to_e(ixGlo1,ixGhi1,ixFimin1,ixFimax1,pwFi%w,px(igrid)%x)
else if (prolongprimitive) then
   call conserve(ixGlo1,ixGhi1,ixFimin1,ixFimax1,pwFi%w,px(igrid)%x,&
      patchfalse)
end if

end subroutine interpolation_copy
!=============================================================================
subroutine interpolation_unlimit(pwFi,ixFimin1,ixFimax1,dxFi1,xFimin1, pwCo,&
   dxCo1,invdxCo1,xComin1)

integer, intent(in) :: ixFimin1,ixFimax1
double precision, intent(in) :: dxFi1, xFimin1, dxCo1,invdxCo1, xComin1
type(walloc) :: pwCo, pwFi

integer :: ixCo1, jxCo1, hxCo1, ixFi1, ix1, idims
double precision :: xCo1, xFi1, eta1
double precision :: slope(nwflux+nwaux,ndim)
!-----------------------------------------------------------------------------
do ixFi1 = ixFimin1,ixFimax1
   ! cell-centered coordinates of fine grid point
   xFi1=xFimin1+(dble(ixFi1)-half)*dxFi1

   ! indices of coarse cell which contains the fine cell
   ixCo1=int((xFi1-xComin1)*invdxCo1)+1

   ! cell-centered coordinate for coarse cell
   xCo1=xComin1+(dble(ixCo1)-half)*dxCo1

   ! normalized distance between fine/coarse cell center
   ! in coarse cell: ranges from -0.5 to 0.5 in each direction
   ! (origin is coarse cell center)
   if (slab) then
      eta1=(xFi1-xCo1)*invdxCo1;
   else
      ix1=2*int((ixFi1+ixMlo1)/2)-ixMlo1;
      eta1=(xFi1-xCo1)*invdxCo1 *two*(one-pgeoFi(igrid)%dvolume(ixFi1) &
         /sum(pgeoFi(igrid)%dvolume(ix1:ix1+1))) 
   end if

   do idims=1,ndim
      hxCo1=ixCo1-kr(1,idims)
      jxCo1=ixCo1+kr(1,idims)

      ! get centered slope
      slope(1:nwflux+nwaux,idims)=half*(pwCo%w(jxCo1,1:nwflux&
         +nwaux)-pwCo%w(hxCo1,1:nwflux+nwaux))
   end do

   ! Interpolate from coarse cell using centered slopes
   pwFi%w(ixFi1,1:nwflux+nwaux)=pwCo%w(ixCo1,1:nwflux+nwaux)+(slope(1:nwflux&
      +nwaux,1)*eta1)
end do

if (amrentropy) then
   call rhos_to_e(ixGlo1,ixGhi1,ixFimin1,ixFimax1,pwFi%w,px(igrid)%x)
else if (prolongprimitive) then
   call conserve(ixGlo1,ixGhi1,ixFimin1,ixFimax1,pwFi%w,px(igrid)%x,&
      patchfalse)
end if

end subroutine interpolation_unlimit
!=============================================================================
subroutine init_bc

integer :: dixBCo, interpolation_order
integer :: ixoldGmin1,ixoldGmax1, ixoldMmin1,ixoldMmax1, nx1, nxCo1
!-----------------------------------------------------------------------------
ixMmin1=ixGmin1+dixB;ixMmax1=ixGmax1-dixB;
ixCoGmin1=1;
ixCoGmax1=ixGmax1/2+dixB;
ixCoMmin1=ixCoGmin1+dixB;ixCoMmax1=ixCoGmax1-dixB;

if (richardson) then
   ixoldGmin1=1; ixoldGmax1=2*ixMmax1;
   ixoldMmin1=ixoldGmin1+dixB;ixoldMmax1=ixoldGmax1-dixB;
end if

nx1=ixMmax1-ixMmin1+1;
nxCo1=nx1/2;

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


ixS_srl_min1(-1)=ixMmin1
ixS_srl_min1(0) =ixMmin1
ixS_srl_min1(1) =ixMmax1+1-dixB
ixS_srl_max1(-1)=ixMmin1-1+dixB
ixS_srl_max1(0) =ixMmax1
ixS_srl_max1(1) =ixMmax1

ixR_srl_min1(-1)=1
ixR_srl_min1(0) =ixMmin1
ixR_srl_min1(1) =ixMmax1+1
ixR_srl_max1(-1)=dixB
ixR_srl_max1(0) =ixMmax1
ixR_srl_max1(1) =ixGmax1


if (levmin/=levmax) then

   ixS_r_min1(-1)=ixCoMmin1
   ixS_r_min1(0) =ixCoMmin1
   ixS_r_min1(1) =ixCoMmax1+1-dixB
   ixS_r_max1(-1)=ixCoMmin1-1+dixB
   ixS_r_max1(0) =ixCoMmax1
   ixS_r_max1(1) =ixCoMmax1

   ixR_r_min1(0)=1
   ixR_r_min1(1)=ixMmin1
   ixR_r_min1(2)=ixMmin1+nxCo1
   ixR_r_min1(3)=ixMmax1+1
   ixR_r_max1(0)=dixB
   ixR_r_max1(1)=ixMmin1-1+nxCo1
   ixR_r_max1(2)=ixMmax1
   ixR_r_max1(3)=ixGmax1

   if (richardson) then
      ixS_old_min1(0)=ixoldMmin1
      ixS_old_min1(1)=ixoldMmin1
      ixS_old_min1(2)=ixoldMmin1+nx1-dixB
      ixS_old_min1(3)=ixoldMmax1+1-dixB
      ixS_old_max1(0)=ixoldMmin1-1+dixB
      ixS_old_max1(1)=ixoldMmin1-1+nx1+dixB
      ixS_old_max1(2)=ixoldMmax1
      ixS_old_max1(3)=ixoldMmax1

      ixR_old_min1(0)=1
      ixR_old_min1(1)=ixMmin1
      ixR_old_min1(2)=1
      ixR_old_min1(3)=ixMmax1+1
      ixR_old_max1(0)=dixB
      ixR_old_max1(1)=ixGmax1
      ixR_old_max1(2)=ixMmax1
      ixR_old_max1(3)=ixGmax1
   else
      ixS_p_min1(0)=ixMmin1-(interpolation_order-1)
      ixS_p_min1(1)=ixMmin1-(interpolation_order-1)
      ixS_p_min1(2)=ixMmin1+nxCo1-dixBCo-(interpolation_order-1)
      ixS_p_min1(3)=ixMmax1+1-dixBCo-(interpolation_order-1)
      ixS_p_max1(0)=ixMmin1-1+dixBCo+(interpolation_order-1)
      ixS_p_max1(1)=ixMmin1-1+nxCo1+dixBCo+(interpolation_order-1)
      ixS_p_max1(2)=ixMmax1+(interpolation_order-1)
      ixS_p_max1(3)=ixMmax1+(interpolation_order-1)

      ixR_p_min1(0)=ixCoMmin1-dixBCo-(interpolation_order-1)
      ixR_p_min1(1)=ixCoMmin1-(interpolation_order-1)
      ixR_p_min1(2)=ixCoMmin1-dixBCo-(interpolation_order-1)
      ixR_p_min1(3)=ixCoMmax1+1-(interpolation_order-1)
      ixR_p_max1(0)=dixB+(interpolation_order-1)
      ixR_p_max1(1)=ixCoMmax1+dixBCo+(interpolation_order-1)
      ixR_p_max1(2)=ixCoMmax1+(interpolation_order-1)
      ixR_p_max1(3)=ixCoMmax1+dixBCo+(interpolation_order-1)
   end if

end if

if (npe>1) then
   do i1=-1,1
      if (i1==0) cycle

      call get_bc_comm_type(type_send_srl(i1),ixS_srl_min1(i1),&
         ixS_srl_max1(i1),ixGmin1,ixGmax1)
      call get_bc_comm_type(type_recv_srl(i1),ixR_srl_min1(i1),&
         ixR_srl_max1(i1),ixGmin1,ixGmax1)

      if (levmin==levmax) cycle

      call get_bc_comm_type(type_send_r(i1),ixS_r_min1(i1),ixS_r_max1(i1),&
         ixCoGmin1,ixCoGmax1)
      do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
         inc1=2*i1+ic1
         call get_bc_comm_type(type_recv_r(inc1),ixR_r_min1(inc1),&
            ixR_r_max1(inc1),ixGmin1,ixGmax1)
         if (richardson) then
            call get_bc_comm_type(type_send_old(inc1), ixS_old_min1(inc1),&
               ixS_old_max1(inc1),ixoldGmin1,ixoldGmax1)
            call get_bc_comm_type(type_recv_old(inc1), ixR_old_min1(inc1),&
               ixR_old_max1(inc1),ixGmin1,ixGmax1)
         else
            call get_bc_comm_type(type_send_p(inc1),ixS_p_min1(inc1),&
               ixS_p_max1(inc1),ixGmin1,ixGmax1)
            call get_bc_comm_type(type_recv_p(inc1),ixR_p_min1(inc1),&
               ixR_p_max1(inc1),ixCoGmin1,ixCoGmax1)
         end if
      end do
   end do
end if

end subroutine init_bc
!=============================================================================
subroutine get_bc_comm_type(comm_type,ixmin1,ixmax1,ixGmin1,ixGmax1)

integer, intent(inout) :: comm_type
integer, intent(in) :: ixmin1,ixmax1, ixGmin1,ixGmax1

integer, dimension(ndim+1) :: size, subsize, start
!-----------------------------------------------------------------------------
size(1)=ixGmax1;
size(ndim+1)=nw
subsize(1)=ixmax1-ixmin1+1;
!subsize(ndim+1)=nwflux
!sorry i need to communicate auxilary lfac for initial guess in inversion:
subsize(ndim+1)=nwflux+nwaux
start(1)=ixmin1-1;
start(ndim+1)=0

if(energyonly) then
  subsize(ndim+1)=1
  start(ndim+1)=e_-1
end if

call MPI_TYPE_CREATE_SUBARRAY(ndim+1,size,subsize,start,MPI_ORDER_FORTRAN,&
    MPI_DOUBLE_PRECISION,comm_type,ierrmpi)
call MPI_TYPE_COMMIT(comm_type,ierrmpi)

end subroutine get_bc_comm_type
!=============================================================================
subroutine put_bc_comm_types
!-----------------------------------------------------------------------------
do i1=-1,1
   if (i1==0) cycle

   call MPI_TYPE_FREE(type_send_srl(i1),ierrmpi)
   call MPI_TYPE_FREE(type_recv_srl(i1),ierrmpi)

   if (levmin==levmax) cycle

   call MPI_TYPE_FREE(type_send_r(i1),ierrmpi)
   do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
      inc1=2*i1+ic1
      call MPI_TYPE_FREE(type_recv_r(inc1),ierrmpi)
      if (richardson) then
         call MPI_TYPE_FREE(type_send_old(inc1),ierrmpi)
         call MPI_TYPE_FREE(type_recv_old(inc1),ierrmpi)
      else
         call MPI_TYPE_FREE(type_send_p(inc1),ierrmpi)
         call MPI_TYPE_FREE(type_recv_p(inc1),ierrmpi)
      end if
   end do
end do

end subroutine put_bc_comm_types
!=============================================================================
subroutine pole_copy(pwrecv,ixRmin1,ixRmax1,pwsend,ixSmin1,ixSmax1)

integer, intent(in) :: ixRmin1,ixRmax1, ixSmin1,ixSmax1
type(walloc) :: pwrecv, pwsend

integer :: iw, iB
!-----------------------------------------------------------------------------
select case (ipole)
case (1)
   iside=int((i1+3)/2)
   iB=2*(1-1)+iside
   do iw=1,nwflux+nwaux
      select case (typeB(iw,iB))
      case ("symm")
         pwrecv%w(ixRmin1:ixRmax1,iw) = pwsend%w(ixSmax1:ixSmin1:-1,iw)
      case ("asymm")
         pwrecv%w(ixRmin1:ixRmax1,iw) =-pwsend%w(ixSmax1:ixSmin1:-1,iw)
      case default
         call mpistop("Boundary condition at pole should be symm or asymm")
      end select
   end do 
end select

end subroutine pole_copy
!=============================================================================
subroutine fix_auxiliary

integer :: ixmin1,ixmax1
!-----------------------------------------------------------------------------
do iigrid=1,igridstail; igrid=igrids(iigrid);
      saveigrid=igrid
      
   do i1=-1,1
      if (i1==0) cycle

      ixmin1=ixR_srl_min1(i1);ixmax1=ixR_srl_max1(i1);
      if(.not.slab)mygeo=>pgeoFi(igrid)
      call getaux(.true.,pwuse(igrid)%w,px(igrid)%x,ixGmin1,ixGmax1,ixmin1,&
         ixmax1,"bc")
   end do
end do

end subroutine fix_auxiliary
!=============================================================================
! end of internal procedures
!=============================================================================
end subroutine getbc
!=============================================================================
subroutine physbound(i1,igrid,isphysbound)
use mod_forest
include 'amrvacdef.f'

integer, intent(in)  :: i1, igrid
logical, intent(out) :: isphysbound
type(tree_node_ptr)  :: tree
integer              :: level, ig1, ign1
!-----------------------------------------------------------------------------
isphysbound = .false.

tree%node => igrid_to_node(igrid,mype)%node
level = tree%node%level
ig1 = tree%node%ig1;

ign1 = ig1 + i1;
if (ign1 .gt. ng1(level) .or. ign1 .lt. 1) isphysbound = .true.

end subroutine physbound
!=============================================================================
