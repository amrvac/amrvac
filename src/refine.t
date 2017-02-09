!=============================================================================
subroutine refine_grid(child_igrid,child_ipe,igrid,ipe,active)

use mod_global_parameters

integer, dimension(2^D&), intent(in) :: child_igrid, child_ipe
integer, intent(in) :: igrid, ipe
logical, intent(in) :: active

integer :: ic^D
!-----------------------------------------------------------------------------

! allocate solution space for new children
{do ic^DB=1,2\}
   call alloc_node(child_igrid(ic^D))
{end do\}

if ((time_advance .and. active).or.convert.or.firstprocess&
   {#IFDEF MAGNETOFRICTION .or. mf_advance}) then
   ! prolong igrid to new children
   call prolong_grid(child_igrid,child_ipe,igrid,ipe)
else
   ! Fill new created children with initial condition
   {do ic^DB=1,2\}
      call initial_condition(child_igrid(ic^D))
   {end do\}
end if

! remove solution space of igrid
call dealloc_node(igrid)

end subroutine refine_grid
!=============================================================================
subroutine prolong_grid(child_igrid,child_ipe,igrid,ipe)
  use mod_physics, only: phys_convert_before_prolong, &
       phys_convert_after_prolong
  use mod_global_parameters

  integer, dimension(2^D&), intent(in) :: child_igrid, child_ipe
  integer, intent(in) :: igrid, ipe

  integer :: ix^L, ichild, ixCo^L, ic^D
  double precision :: dxCo^D, xComin^D, dxFi^D, xFimin^D
  !-----------------------------------------------------------------------------
  ! TODO: discuss/clean commented code below
  if (prolongation_method=="linear") then
     ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
     {#IFDEF EVOLVINGBOUNDARY
     if (phyboundblock(igrid)) then
        ix^L=ixG^LL;
     else
        ix^L=ixM^LL^LADD1;
     end if
     }{#IFNDEF EVOLVINGBOUNDARY
     ix^L=ixM^LL^LADD1;
     }

     call phys_convert_before_prolong(ixG^LL,ix^L,pw(igrid)%w,px(igrid)%x)

     xComin^D=rnode(rpxmin^D_,igrid)\
     dxCo^D=rnode(rpdx^D_,igrid)\
     {#IFDEF STRETCHGRID
     logG=logGs(node(plevel_,igrid))
     qst=qsts(node(plevel_,igrid))
     }
  end if

  {do ic^DB=1,2\}
  ichild=child_igrid(ic^D)

  ixComin^D=ixMlo^D+(ic^D-1)*(ixMhi^D-ixMlo^D+1)/2\
  ixComax^D=ixMhi^D+(ic^D-2)*(ixMhi^D-ixMlo^D+1)/2\

  if (prolongation_method=="linear") then
     xFimin^D=rnode(rpxmin^D_,ichild)\
     dxFi^D=rnode(rpdx^D_,ichild)\
     {#IFDEF EVOLVINGBOUNDARY
     if (phyboundblock(ichild)) then
        call prolong_2ab(pw(igrid)%w,px(igrid)%x,ixCo^L,pw(ichild)%w,px(ichild)%x, &
             dxCo^D,xComin^D,dxFi^D,xFimin^D,ichild)
     else
        call prolong_2nd(pw(igrid)%w,px(igrid)%x,ixCo^L,pw(ichild)%w,px(ichild)%x, &
             dxCo^D,xComin^D,dxFi^D,xFimin^D,ichild)
     end if
     }{#IFNDEF EVOLVINGBOUNDARY
     call prolong_2nd(pw(igrid)%w,px(igrid)%x,ixCo^L,pw(ichild)%w,px(ichild)%x, &
          dxCo^D,xComin^D,dxFi^D,xFimin^D,ichild)
     }
  else
     call prolong_1st(pw(igrid)%w,ixCo^L,pw(ichild)%w,px(ichild)%x)
  end if
  {end do\}

  if (prolongation_method=="linear") then
     call phys_convert_after_prolong(ixG^LL,ix^L,pw(igrid)%w,px(igrid)%x)
     ! TODO: clean this up
     !    if (amrentropy) then
     !       call rhos_to_e(ixG^LL,ix^L,pw(igrid)%w,px(igrid)%x)
     !    else if (prolongprimitive) then
     !       call conserve(ixG^LL,ix^L,pw(igrid)%w,px(igrid)%x,patchfalse)
     !    end if
  end if

end subroutine prolong_grid
!=============================================================================
subroutine prolong_2ab(wCo,xCo,ixCo^L,wFi,xFi,dxCo^D,xComin^D,dxFi^D,xFimin^D,igridFi)
! interpolate children blocks including ghost cells

use mod_physics, only: phys_convert_after_prolong
use mod_global_parameters

integer, intent(in) :: ixCo^L, igridFi
double precision, intent(in) :: dxCo^D, xComin^D, dxFi^D, xFimin^D
double precision, intent(in) :: wCo(ixG^T,nw), xCo(ixG^T,1:ndim), xFi(ixG^T,1:ndim)
double precision, intent(inout) :: wFi(ixG^T,nw)

integer :: ixCo^D, jxCo^D, hxCo^D, ixFi^D, ix^D, idim, iw
integer :: ixFi^L, ixCg^L, el
double precision :: slopeL, slopeR, slopeC, signC, signR
double precision :: slope(nw,ndim)
double precision :: xCo^D, xFi^D, eta^D, invdxCo^D
{#IFDEF STRETCHGRID
double precision :: logGh,qsth
}
!-----------------------------------------------------------------------------
{#IFDEF STRETCHGRID
qsth=dsqrt(qst)
logGh=2.d0*(qsth-1.d0)/(qsth+1.d0)
}
invdxCo^D=1.d0/dxCo^D;
el=ceiling(real(nghostcells)/2.)
ixCgmin^D=ixComin^D-el\
ixCgmax^D=ixComax^D+el\
{do ixCo^DB = ixCg^LIM^DB
   ! cell-centered coordinates of coarse grid point
   xCo^DB=xComin^DB+(dble(ixCo^DB-nghostcells)-half)*dxCo^DB

   ixFi^DB=2*(ixCo^DB-ixComin^DB)+ixMlo^DB\}
{#IFDEF STRETCHGRID
   xCo1=xComin1/(one-half*logG)*qst**(ixCo1-nghostcells-1)
}

   do idim=1,ndim
      hxCo^D=ixCo^D-kr(^D,idim)\
      jxCo^D=ixCo^D+kr(^D,idim)\

      do iw=1,nw
         slopeL=wCo(ixCo^D,iw)-wCo(hxCo^D,iw)
         slopeR=wCo(jxCo^D,iw)-wCo(ixCo^D,iw)
         slopeC=half*(slopeR+slopeL)

         ! get limited slope
         signR=sign(one,slopeR)
         signC=sign(one,slopeC)
         select case(typeprolonglimit)
         case('minmod')
           slope(iw,idim)=signR*max(zero,min(dabs(slopeR), &
                                             signR*slopeL))
         case('woodward')
           slope(iw,idim)=two*signR*max(zero,min(dabs(slopeR), &
                              signR*slopeL,signR*half*slopeC))
         case('mcbeta')
           slope(iw,idim)=signR*max(zero,min(mcbeta*dabs(slopeR), &
                              mcbeta*signR*slopeL,signR*slopeC))
         case('koren')
           slope(iw,idim)=signR*max(zero,min(two*signR*slopeL, &
            (dabs(slopeR)+two*slopeL*signR)*third,two*dabs(slopeR)))
         case default
           slope(iw,idim)=signC*max(zero,min(dabs(slopeC), &
                             signC*slopeL,signC*slopeR))
         end select
      end do
   end do
   {do ix^DB=ixFi^DB,ixFi^DB+1
      if (ixFi^DB==0) cycle
      ! cell-centered coordinates of fine grid point
      xFi^DB=xFimin^DB+(dble(ix^DB-nghostcells)-half)*dxFi^DB\}
{#IFDEF STRETCHGRID
      xFi1=xFimin1/(one-half*logGh)*qsth**(ixFi1-nghostcells-1)
}

      ! normalized distance between fine/coarse cell center
      ! in coarse cell: ranges from -0.5 to 0.5 in each direction
      ! (origin is coarse cell center)
      if (slab) then
         eta^D=(xFi^D-xCo^D)*invdxCo^D;
      else
         {eta^D=(xFi^D-xCo^D)*invdxCo^D &
               *two*(one-pgeo(igridFi)%dvolume(ix^DD) &
               /sum(pgeo(igridFi)%dvolume(ixFi^D:ixFi^D+1^D%ix^DD))) \}
{#IFDEF STRETCHGRID
         eta1=(xFi1-xCo1)/(logG*xCo1) &
               *two*(one-pgeo(igridFi)%dvolume(ix^D) &
               /sum(pgeo(igridFi)%dvolume(ixFi1:ixFi1+1^%1ix^D)))
}
      end if

      wFi(ix^D,1:nw) = wCo(ixCo^D,1:nw) &
                            + {(slope(1:nw,^D)*eta^D)+}
   {end do\}
{end do\}

call phys_convert_after_prolong(ixG^LL,ixM^LL,wFi,xFi)

end subroutine prolong_2ab
!=============================================================================
subroutine prolong_2nd(wCo,xCo,ixCo^L,wFi,xFi,dxCo^D,xComin^D,dxFi^D,xFimin^D,igridFi)

use mod_physics, only: phys_convert_after_prolong
use mod_global_parameters

integer, intent(in) :: ixCo^L, igridFi
double precision, intent(in) :: dxCo^D, xComin^D, dxFi^D, xFimin^D
double precision, intent(in) :: wCo(ixG^T,nw), xCo(ixG^T,1:ndim), xFi(ixG^T,1:ndim)
double precision, intent(inout) :: wFi(ixG^T,nw)

integer :: ixCo^D, jxCo^D, hxCo^D, ixFi^D, ix^D, idim, iw
integer :: ixFi^L
double precision :: slopeL, slopeR, slopeC, signC, signR
double precision :: slope(nw,ndim)
double precision :: xCo^D, xFi^D, eta^D, invdxCo^D
{#IFDEF STRETCHGRID
double precision :: logGh,qsth
}
!-----------------------------------------------------------------------------
{#IFDEF STRETCHGRID
qsth=dsqrt(qst)
logGh=2.d0*(qsth-1.d0)/(qsth+1.d0)
}
invdxCo^D=1.d0/dxCo^D;
{do ixCo^DB = ixCo^LIM^DB
   ! cell-centered coordinates of coarse grid point
   xCo^DB=xComin^DB+(dble(ixCo^DB-nghostcells)-half)*dxCo^DB

   ixFi^DB=2*(ixCo^DB-ixComin^DB)+ixMlo^DB\}
{#IFDEF STRETCHGRID
   xCo1=xComin1/(one-half*logG)*qst**(ixCo1-nghostcells-1)
}

   do idim=1,ndim
      hxCo^D=ixCo^D-kr(^D,idim)\
      jxCo^D=ixCo^D+kr(^D,idim)\

      do iw=1,nw
         slopeL=wCo(ixCo^D,iw)-wCo(hxCo^D,iw)
         slopeR=wCo(jxCo^D,iw)-wCo(ixCo^D,iw)
         slopeC=half*(slopeR+slopeL)

         ! get limited slope
         signR=sign(one,slopeR)
         signC=sign(one,slopeC)
         select case(typeprolonglimit)
         case('minmod')
           slope(iw,idim)=signR*max(zero,min(dabs(slopeR), &
                                             signR*slopeL))
         case('woodward')
           slope(iw,idim)=two*signR*max(zero,min(dabs(slopeR), &
                              signR*slopeL,signR*half*slopeC))
         case('mcbeta')
           slope(iw,idim)=signR*max(zero,min(mcbeta*dabs(slopeR), &
                              mcbeta*signR*slopeL,signR*slopeC))
         case('koren')
           slope(iw,idim)=signR*max(zero,min(two*signR*slopeL, &
            (dabs(slopeR)+two*slopeL*signR)*third,two*dabs(slopeR)))
         case default
           slope(iw,idim)=signC*max(zero,min(dabs(slopeC), &
                             signC*slopeL,signC*slopeR))
         end select
      end do
   end do
   {do ix^DB=ixFi^DB,ixFi^DB+1
      ! cell-centered coordinates of fine grid point
      xFi^DB=xFimin^DB+(dble(ix^DB-nghostcells)-half)*dxFi^DB\}
{#IFDEF STRETCHGRID
      xFi1=xFimin1/(one-half*logGh)*qsth**(ixFi1-nghostcells-1)
}

      ! normalized distance between fine/coarse cell center
      ! in coarse cell: ranges from -0.5 to 0.5 in each direction
      ! (origin is coarse cell center)
      if (slab) then
         eta^D=(xFi^D-xCo^D)*invdxCo^D;
      else
         {eta^D=(xFi^D-xCo^D)*invdxCo^D &
               *two*(one-pgeo(igridFi)%dvolume(ix^DD) &
               /sum(pgeo(igridFi)%dvolume(ixFi^D:ixFi^D+1^D%ix^DD))) \}
{#IFDEF STRETCHGRID
         eta1=(xFi1-xCo1)/(logG*xCo1) &
               *two*(one-pgeo(igridFi)%dvolume(ix^D) &
               /sum(pgeo(igridFi)%dvolume(ixFi1:ixFi1+1^%1ix^D)))
}
      end if

      wFi(ix^D,1:nw) = wCo(ixCo^D,1:nw) &
                            + {(slope(1:nw,^D)*eta^D)+}
   {end do\}
{end do\}

call phys_convert_after_prolong(ixG^LL,ixM^LL,wFi,xFi)

end subroutine prolong_2nd
!=============================================================================
subroutine prolong_1st(wCo,ixCo^L,wFi,xFi)

use mod_global_parameters

integer, intent(in) :: ixCo^L
double precision, intent(in) :: wCo(ixG^T,nw), xFi(ixG^T,1:ndim)
double precision, intent(out) :: wFi(ixG^T,nw)

integer :: ixCo^D, ixFi^D, iw
integer :: ixFi^L
!-----------------------------------------------------------------------------
{do ixCo^DB = ixCo^LIM^DB
   ixFi^DB=2*(ixCo^DB-ixComin^DB)+ixMlo^DB\}
   forall(iw=1:nw) wFi(ixFi^D:ixFi^D+1,iw)=wCo(ixCo^D,iw)
{end do\}

end subroutine prolong_1st
!=============================================================================
