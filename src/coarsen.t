!=============================================================================
subroutine coarsen_grid_siblings(igrid,ipe,child_igrid,child_ipe,active)

include 'amrvacdef.f'

integer, intent(in) :: igrid, ipe
integer, dimension(2^D&), intent(in) :: child_igrid, child_ipe
logical, intent(in) :: active

integer :: igridFi, ipeFi, ixCo^L, ixCoG^L, ixCoM^L, ic^D
!-----------------------------------------------------------------------------
if(addmpibarrier) call MPI_BARRIER(icomm,ierrmpi)

if (ipe==mype) call alloc_node(igrid)

! New passive cell, coarsen from initial condition:
if (.not. active) then

   if (ipe == mype) then
      call initial_condition(igrid)
      
      {do ic^DB=1,2\}
      igridFi=child_igrid(ic^D)
      ipeFi=child_ipe(ic^D)
      if (ipeFi==mype) then
         ! remove solution space of child      
         call dealloc_node(igridFi)
      end if
      {end do\}
      
   end if

   return
end if



{do ic^DB=1,2\}
   igridFi=child_igrid(ic^D)
   ipeFi=child_ipe(ic^D)

   if (ipeFi==mype) then
      ^D&dxlevel(^D)=rnode(rpdx^D_,igridFi);
      if (ipe==mype) then
         ixComin^D=ixMlo^D+(ic^D-1)*(ixMhi^D-ixMlo^D+1)/2;
         ixComax^D=ixMhi^D+(ic^D-2)*(ixMhi^D-ixMlo^D+1)/2;

         call coarsen_grid(pw(igridFi)%w,px(igridFi)%x,ixG^LL,ixM^LL,pw(igrid)%w,px(igrid)%x,ixG^LL, &
                     ixCo^L,pgeo(igridFi),pgeo(igrid),restrictprimitive,.false.)

         ! remove solution space of child
         call dealloc_node(igridFi)
      else
         ixCoGmin^D=1;
         ixCoGmax^D=ixGhi^D/2+dixB;
         ixCoM^L=ixCoG^L^LSUBdixB;
         call coarsen_grid(pw(igridFi)%w,px(igridFi)%x,ixG^LL,ixM^LL,pwCoarse(igridFi)%w,pxCoarse(igridFi)%x, &
                           ixCoG^L,ixCoM^L,pgeo(igridFi),pgeoCoarse(igridFi),&
                           restrictprimitive,.false.)

         itag=ipeFi*ngridshi+igridFi
         isend=isend+1
         call MPI_ISEND(pwCoarse(igridFi)%w,1,type_coarse_block,ipe,itag, &
                        icomm,sendrequest(isend),ierrmpi)
      end if
   else
      if (ipe==mype) then
         itag=ipeFi*ngridshi+igridFi
         irecv=irecv+1
         call MPI_IRECV(pw(igrid)%w,1,type_sub_block(ic^D),ipeFi,itag, &
                        icomm,recvrequest(irecv),ierrmpi)
      end if
   end if
{end do\}

if(addmpibarrier) call MPI_BARRIER(icomm,ierrmpi)
end subroutine coarsen_grid_siblings
!=============================================================================
subroutine coarsen_grid(wFi,xFi,ixFiG^L,ixFi^L,wCo,xCo,ixCoG^L,ixCo^L,&
                        pgeogrid,pgeoCoarsegrid,coarsenprim,keepFi)

include 'amrvacdef.f'

integer, intent(in) :: ixFiG^L, ixFi^L, ixCoG^L, ixCo^L
double precision, intent(inout) :: wFi(ixFiG^S,1:nw), xFi(ixFiG^S,1:ndim)
double precision,intent(inout) :: wCo(ixCoG^S,1:nw), xCo(ixCoG^S,1:ndim)
type(geoalloc) :: pgeogrid, pgeoCoarsegrid
logical, intent(in) :: coarsenprim, keepFi

integer :: ixCo^D, ixFi^D, iw
double precision :: CoFiratio
!-----------------------------------------------------------------------------
! coarsen by 2 in every direction - conservatively

if (amrentropy) then
   call e_to_rhos(ixFiG^L,ixFi^L,wFi,xFi)
else if (coarsenprim) then
   call primitive(ixFiG^L,ixFi^L,wFi,xFi)
end if

if (slab) then
   CoFiratio=one/dble(2**ndim)
   do iw=1,nw
      {do ixCo^DB = ixCo^LIM^DB
         ixFi^DB=2*(ixCo^DB-ixComin^DB)+ixFimin^DB\}
         wCo(ixCo^D,iw)=sum(wFi(ixFi^D:ixFi^D+1,iw))*CoFiratio
      {end do\}
   end do
else
   do iw=1,nw
      {do ixCo^DB = ixCo^LIM^DB
         ixFi^DB=2*(ixCo^DB-ixComin^DB)+ixFimin^DB\}
         wCo(ixCo^D,iw)= &
             sum(pgeogrid%dvolume(ixFi^D:ixFi^D+1)*wFi(ixFi^D:ixFi^D+1,iw)) &
            /pgeoCoarsegrid%dvolume(ixCo^D)
      {end do\}
   end do
end if

if (amrentropy) then
   if (keepFi) call rhos_to_e(ixFiG^L,ixFi^L,wFi,xFi)
   call rhos_to_e(ixCoG^L,ixCo^L,wCo,xCo)
else if (coarsenprim) then
   if (keepFi) call conserve(ixFiG^L,ixFi^L,wFi,xFi,patchfalse)
   call conserve(ixCoG^L,ixCo^L,wCo,xCo,patchfalse)
end if

end subroutine coarsen_grid
!=============================================================================
