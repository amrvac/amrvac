!> coarsen one grid to its coarser representative
subroutine coarsen_grid(sFi,ixFiG^L,ixFi^L,sCo,ixCoG^L,ixCo^L,do_m,do_w)
  use mod_global_parameters
  use mod_physics

  type(state), intent(inout) :: sFi, sCo
  integer, intent(in)        :: ixFiG^L, ixFi^L, ixCoG^L, ixCo^L
  logical, intent(in)        :: do_w, do_m

  integer                    :: ixCo^D, ixFi^D, iw, nwmin, nwmax
  double precision           :: CoFiratio

  ! if coarsen metric variables
  if (do_m) then
    associate(wFi=>sFi%metric%vars(ixFiG^S,1:nmetric), &
              wCo=>sCo%metric%vars(ixCoG^S,1:nmetric))
    ! coarsen by 2 in every direction 
    nwmin = 1
    nwmax = nmetric
    if(slab_uniform) then
      CoFiratio=1.0d0/dble(2**ndim)
      do iw=nwmin,nwmax
         {do ixCo^DB = ixCo^LIM^DB
            ixFi^DB=2*(ixCo^DB-ixComin^DB)+ixFimin^DB\}
            wCo(ixCo^D,iw)=sum(wFi(ixFi^D:ixFi^D+1,iw))*CoFiratio
         {end do\}
      end do
    else
      do iw=nwmin,nwmax
        {do ixCo^DB = ixCo^LIM^DB
           ixFi^DB=2*(ixCo^DB-ixComin^DB)+ixFimin^DB\}
           wCo(ixCo^D,iw)= &
               sum(sFi%mesh%dvolume(ixFi^D:ixFi^D+1)*wFi(ixFi^D:ixFi^D+1,iw)) &
              /sCo%mesh%dvolume(ixCo^D)
        {end do\}
      end do
    end if
    end associate 
  end if

  if (.not. do_w) return

  associate(wFi=>sFi%w(ixFiG^S,1:nw), wCo=>sCo%w(ixCoG^S,1:nw))
  staggered: associate(wFis=>sFi%ws,wCos=>sCo%ws)
  ! coarsen by 2 in every direction - conservatively

  nwmin = 1
  nwmax = nwflux

  if(slab_uniform) then
    CoFiratio=1.0d0/dble(2**ndim)
    do iw=nwmin,nwmax
       {do ixCo^DB = ixCo^LIM^DB
          ixFi^DB=2*(ixCo^DB-ixComin^DB)+ixFimin^DB\}
          wCo(ixCo^D,iw)=sum(wFi(ixFi^D:ixFi^D+1,iw))*CoFiratio
       {end do\}
    end do
  else
    do iw=nwmin,nwmax
      {do ixCo^DB = ixCo^LIM^DB
         ixFi^DB=2*(ixCo^DB-ixComin^DB)+ixFimin^DB\}
         wCo(ixCo^D,iw)= &
             sum(sFi%mesh%dvolume(ixFi^D:ixFi^D+1)*wFi(ixFi^D:ixFi^D+1,iw)) &
            /sCo%mesh%dvolume(ixCo^D)
      {end do\}
    end do
  end if

  if(stagger_grid) then
    do iw=1,nws
      ! Start one layer before
      {do ixCo^DB = ixComin^DB-kr(^DB,iw),ixComax^DB
         ixFi^DB=2*(ixCo^DB-ixComin^DB+kr(^DB,iw))+ixFimin^DB-kr(^DB,iw)\}
         ! This if statement catches the axis where surface is zero:
         if (sCo%mesh%surfaceC(ixCo^D,iw)>1.0d-9*sCo%mesh%dvolume(ixCo^D)) then ! Normal case
           wCos(ixCo^D,iw)=sum(sFi%mesh%surfaceC(ixFi^D:ixFi^D+1-kr(iw,^D),iw)*wFis(ixFi^D:ixFi^D+1-kr(iw,^D),iw)) &
                /sCo%mesh%surfaceC(ixCo^D,iw)
         else ! On axis
           wCos(ixCo^D,iw)=zero
         end if
      {end do\}
    end do
    ! average to fill cell-centred values
    call phys_face_to_center(ixCo^L,sCo)
  end if

  end associate staggered
  end associate
end subroutine coarsen_grid
