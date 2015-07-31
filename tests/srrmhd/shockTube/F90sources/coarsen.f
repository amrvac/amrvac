!=============================================================================
subroutine coarsen_grid_siblings(igrid,ipe,child_igrid,child_ipe,active)

include 'amrvacdef.f'

integer, intent(in) :: igrid, ipe
integer, dimension(2), intent(in) :: child_igrid, child_ipe
logical, intent(in) :: active

integer :: igridFi, ipeFi, ixComin1,ixComax1, ixCoGmin1,ixCoGmax1, ixCoMmin1,&
   ixCoMmax1, ic1
!-----------------------------------------------------------------------------
if(addmpibarrier) call MPI_BARRIER(icomm,ierrmpi)

if (ipe==mype) call alloc_node(igrid)

! New passive cell, coarsen from initial condition:
if (.not. active) then

   if (ipe == mype) then
      call initial_condition(igrid)
      
      do ic1=1,2
      igridFi=child_igrid(ic1)
      ipeFi=child_ipe(ic1)
      if (ipeFi==mype) then
         ! remove solution space of child      
         call dealloc_node(igridFi)
      end if
      end do
      
   end if

   return
end if



do ic1=1,2
   igridFi=child_igrid(ic1)
   ipeFi=child_ipe(ic1)

   if (ipeFi==mype) then
      dxlevel(1)=rnode(rpdx1_,igridFi);
      if (ipe==mype) then
         ixComin1=ixMlo1+(ic1-1)*(ixMhi1-ixMlo1+1)/2;
         ixComax1=ixMhi1+(ic1-2)*(ixMhi1-ixMlo1+1)/2;

         call coarsen_grid(pw(igridFi)%w,px(igridFi)%x,ixGlo1,ixGhi1,ixMlo1,&
            ixMhi1,pw(igrid)%w,px(igrid)%x,ixGlo1,ixGhi1, ixComin1,ixComax1,&
            pgeo(igridFi),pgeo(igrid),restrictprimitive,.false.)

         ! remove solution space of child
         call dealloc_node(igridFi)
      else
         ixCoGmin1=1;
         ixCoGmax1=ixGhi1/2+dixB;
         ixCoMmin1=ixCoGmin1+dixB;ixCoMmax1=ixCoGmax1-dixB;
         call coarsen_grid(pw(igridFi)%w,px(igridFi)%x,ixGlo1,ixGhi1,ixMlo1,&
            ixMhi1,pwCoarse(igridFi)%w,pxCoarse(igridFi)%x, ixCoGmin1,&
            ixCoGmax1,ixCoMmin1,ixCoMmax1,pgeo(igridFi),pgeoCoarse(igridFi),&
            restrictprimitive,.false.)

         itag=ipeFi*ngridshi+igridFi
         isend=isend+1
         call MPI_ISEND(pwCoarse(igridFi)%w,1,type_coarse_block,ipe,itag,&
             icomm,sendrequest(isend),ierrmpi)
      end if
   else
      if (ipe==mype) then
         itag=ipeFi*ngridshi+igridFi
         irecv=irecv+1
         call MPI_IRECV(pw(igrid)%w,1,type_sub_block(ic1),ipeFi,itag, icomm,&
            recvrequest(irecv),ierrmpi)
      end if
   end if
end do

if(addmpibarrier) call MPI_BARRIER(icomm,ierrmpi)
end subroutine coarsen_grid_siblings
!=============================================================================
subroutine coarsen_grid(wFi,xFi,ixFiGmin1,ixFiGmax1,ixFimin1,ixFimax1,wCo,xCo,&
   ixCoGmin1,ixCoGmax1,ixComin1,ixComax1,pgeogrid,pgeoCoarsegrid,coarsenprim,&
   keepFi)

include 'amrvacdef.f'

integer, intent(in) :: ixFiGmin1,ixFiGmax1, ixFimin1,ixFimax1, ixCoGmin1,&
   ixCoGmax1, ixComin1,ixComax1
double precision, intent(inout) :: wFi(ixFiGmin1:ixFiGmax1,1:nw),&
    xFi(ixFiGmin1:ixFiGmax1,1:ndim)
double precision,intent(inout) :: wCo(ixCoGmin1:ixCoGmax1,1:nw),&
    xCo(ixCoGmin1:ixCoGmax1,1:ndim)
type(geoalloc) :: pgeogrid, pgeoCoarsegrid
logical, intent(in) :: coarsenprim, keepFi

integer :: ixCo1, ixFi1, iw
double precision :: CoFiratio
!-----------------------------------------------------------------------------
! coarsen by 2 in every direction - conservatively

if (amrentropy) then
   call e_to_rhos(ixFiGmin1,ixFiGmax1,ixFimin1,ixFimax1,wFi,xFi)
else if (coarsenprim) then
   call primitive(ixFiGmin1,ixFiGmax1,ixFimin1,ixFimax1,wFi,xFi)
end if

if (slab) then
   CoFiratio=one/dble(2**ndim)
   do iw=1,nw
      do ixCo1 = ixComin1,ixComax1
         ixFi1=2*(ixCo1-ixComin1)+ixFimin1
         wCo(ixCo1,iw)=sum(wFi(ixFi1:ixFi1+1,iw))*CoFiratio
      end do
   end do
else
   do iw=1,nw
      do ixCo1 = ixComin1,ixComax1
         ixFi1=2*(ixCo1-ixComin1)+ixFimin1
         wCo(ixCo1,iw)= sum(pgeogrid%dvolume(ixFi1:ixFi1+1)*wFi(ixFi1:ixFi1&
            +1,iw)) /pgeoCoarsegrid%dvolume(ixCo1)
      end do
   end do
end if

if (amrentropy) then
   if (keepFi) call rhos_to_e(ixFiGmin1,ixFiGmax1,ixFimin1,ixFimax1,wFi,xFi)
   call rhos_to_e(ixCoGmin1,ixCoGmax1,ixComin1,ixComax1,wCo,xCo)
else if (coarsenprim) then
   if (keepFi) call conserve(ixFiGmin1,ixFiGmax1,ixFimin1,ixFimax1,wFi,xFi,&
      patchfalse)
   call conserve(ixCoGmin1,ixCoGmax1,ixComin1,ixComax1,wCo,xCo,patchfalse)
end if

end subroutine coarsen_grid
!=============================================================================
