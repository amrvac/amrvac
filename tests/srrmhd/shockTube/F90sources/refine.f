!=============================================================================
subroutine refine_grid(child_igrid,child_ipe,igrid,ipe,active)

include 'amrvacdef.f'

integer, dimension(2), intent(in) :: child_igrid, child_ipe
integer, intent(in) :: igrid, ipe
logical, intent(in) :: active

integer :: ic1
!-----------------------------------------------------------------------------

! allocate solution space for new children
do ic1=1,2
   call alloc_node(child_igrid(ic1))
end do

if ((time_advance .and. active).or.convert.or.firstprocess) then
   ! prolong igrid to new children
   call prolong_grid(child_igrid,child_ipe,igrid,ipe)
else
   ! Fill new created children with initial condition
   do ic1=1,2
      call initial_condition(child_igrid(ic1))
   end do
end if

! remove solution space of igrid
call dealloc_node(igrid)

end subroutine refine_grid
!=============================================================================
subroutine prolong_grid(child_igrid,child_ipe,igrid,ipe)

include 'amrvacdef.f'

integer, dimension(2), intent(in) :: child_igrid, child_ipe
integer, intent(in) :: igrid, ipe

integer :: ixmin1,ixmax1, ichild, ixComin1,ixComax1, ic1
double precision :: dxCo1, xComin1, dxFi1, xFimin1
!-----------------------------------------------------------------------------
if (typegridfill=="linear") then
   dxlevel(1)=rnode(rpdx1_,igrid);
   if (amrentropy) then
      ixmin1=ixMlo1-1;ixmax1=ixMhi1+1;
      call e_to_rhos(ixGlo1,ixGhi1,ixmin1,ixmax1,pw(igrid)%w,px(igrid)%x)
   else if (prolongprimitive) then
      ixmin1=ixMlo1-1;ixmax1=ixMhi1+1;
      call primitive(ixGlo1,ixGhi1,ixmin1,ixmax1,pw(igrid)%w,px(igrid)%x)
   end if

   xComin1=rnode(rpxmin1_,igrid)
   dxCo1=rnode(rpdx1_,igrid)
end if

do ic1=1,2
   ichild=child_igrid(ic1)

   ixComin1=ixMlo1+(ic1-1)*(ixMhi1-ixMlo1+1)/2
   ixComax1=ixMhi1+(ic1-2)*(ixMhi1-ixMlo1+1)/2

   if (typegridfill=="linear") then
      xFimin1=rnode(rpxmin1_,ichild)
      dxFi1=rnode(rpdx1_,ichild)

      call prolong_2nd(pw(igrid)%w,px(igrid)%x,ixComin1,ixComax1,pw(ichild)%w,&
         px(ichild)%x, dxCo1,xComin1,dxFi1,xFimin1,ichild)
   else
      call prolong_1st(pw(igrid)%w,ixComin1,ixComax1,pw(ichild)%w,&
         px(ichild)%x)
   end if
end do

if (typegridfill=="linear") then
   if (amrentropy) then
      call rhos_to_e(ixGlo1,ixGhi1,ixmin1,ixmax1,pw(igrid)%w,px(igrid)%x)
   else if (prolongprimitive) then
      call conserve(ixGlo1,ixGhi1,ixmin1,ixmax1,pw(igrid)%w,px(igrid)%x,&
         patchfalse)
   end if
end if

end subroutine prolong_grid
!=============================================================================
subroutine prolong_2nd(wCo,xCo,ixComin1,ixComax1,wFi,xFi,dxCo1,xComin1,dxFi1,&
   xFimin1,igridFi)

include 'amrvacdef.f'

integer, intent(in) :: ixComin1,ixComax1, igridFi
double precision, intent(in) :: dxCo1, xComin1, dxFi1, xFimin1
double precision, intent(in) :: wCo(ixGlo1:ixGhi1,nw), xCo(ixGlo1:ixGhi1,&
   1:ndim), xFi(ixGlo1:ixGhi1,1:ndim)
double precision, intent(inout) :: wFi(ixGlo1:ixGhi1,nw)

integer :: ixCo1, jxCo1, hxCo1, ixFi1, ix1, idim, iw
integer :: ixFimin1,ixFimax1
double precision :: slopeL, slopeR, slopeC, signC, signR
double precision :: slope(nw,ndim)
double precision :: xCo1, xFi1, eta1, invdxCo1
!-----------------------------------------------------------------------------
invdxCo1=dxCo1;
do ixCo1 = ixComin1,ixComax1
   ! cell-centered coordinates of coarse grid point
   xCo1=xComin1+(dble(ixCo1-dixB)-half)*dxCo1

   ixFi1=2*(ixCo1-ixComin1)+ixMlo1

   do idim=1,ndim
      hxCo1=ixCo1-kr(1,idim)
      jxCo1=ixCo1+kr(1,idim)

      do iw=1,nw
         slopeL=wCo(ixCo1,iw)-wCo(hxCo1,iw)
         slopeR=wCo(jxCo1,iw)-wCo(ixCo1,iw)
         slopeC=half*(slopeR+slopeL)

         ! get limited slope
         signR=sign(one,slopeR)
         signC=sign(one,slopeC)
         select case(typeprolonglimit)
         case('minmod')
           slope(iw,idim)=signR*max(zero,min(dabs(slopeR), signR*slopeL))
         case('woodward')
           slope(iw,idim)=two*signR*max(zero,min(dabs(slopeR), signR&
              *slopeL,signR*half*slopeC))
         case('mcbeta')
           slope(iw,idim)=signR*max(zero,min(mcbeta*dabs(slopeR), mcbeta&
              *signR*slopeL,signR*slopeC))
         case('koren')
           slope(iw,idim)=signR*max(zero,min(two*signR*slopeL, (dabs(slopeR)&
              +two*slopeL*signR)*third,two*dabs(slopeR)))
         case default
           slope(iw,idim)=signC*max(zero,min(dabs(slopeC), signC*slopeL,signC&
              *slopeR))
         end select
      end do
   end do
   do ix1=ixFi1,ixFi1+1
      ! cell-centered coordinates of fine grid point
      xFi1=xFimin1+(dble(ix1-dixB)-half)*dxFi1

      ! normalized distance between fine/coarse cell center
      ! in coarse cell: ranges from -0.5 to 0.5 in each direction
      ! (origin is coarse cell center)
      if (slab) then
         eta1=(xFi1-xCo1)*invdxCo1;
      else
         eta1=(xFi1-xCo1)*invdxCo1 *two*(one-pgeo(igridFi)%dvolume(ix1) &
            /sum(pgeo(igridFi)%dvolume(ixFi1:ixFi1+1))) 
      end if

      wFi(ix1,1:nw) = wCo(ixCo1,1:nw) + (slope(1:nw,1)*eta1)
   end do
end do

if (amrentropy) then
   call rhos_to_e(ixGlo1,ixGhi1,ixMlo1,ixMhi1,wFi,xFi)
else if (prolongprimitive) then
   call conserve(ixGlo1,ixGhi1,ixMlo1,ixMhi1,wFi,xFi,patchfalse)
end if

end subroutine prolong_2nd
!=============================================================================
subroutine prolong_1st(wCo,ixComin1,ixComax1,wFi,xFi)

include 'amrvacdef.f'

integer, intent(in) :: ixComin1,ixComax1
double precision, intent(in) :: wCo(ixGlo1:ixGhi1,nw), xFi(ixGlo1:ixGhi1,&
   1:ndim)
double precision, intent(out) :: wFi(ixGlo1:ixGhi1,nw)

integer :: ixCo1, ixFi1, iw
integer :: ixFimin1,ixFimax1
!-----------------------------------------------------------------------------
do ixCo1 = ixComin1,ixComax1
   ixFi1=2*(ixCo1-ixComin1)+ixMlo1
   forall(iw=1:nw) wFi(ixFi1:ixFi1+1,iw)=wCo(ixCo1,iw)
end do

end subroutine prolong_1st
!=============================================================================
