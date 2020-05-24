!> coarsen one grid to its coarser representative
subroutine coarsen_grid(sFi,ixFiG^L,ixFi^L,sCo,ixCoG^L,ixCo^L)
  use mod_global_parameters
  use mod_physics

  type(state), intent(inout)      :: sFi, sCo
  integer, intent(in) :: ixFiG^L, ixFi^L, ixCoG^L, ixCo^L

  integer :: ixCo^D, ixFi^D, iw
  double precision :: CoFiratio
  double precision :: B_energy_change(ixCoG^S)

  associate(wFi=>sFi%w(ixFiG^S,1:nw), wCo=>sCo%w(ixCoG^S,1:nw))
  staggered: associate(wFis=>sFi%ws,wCos=>sCo%ws)
  ! coarsen by 2 in every direction - conservatively

  if(coarsenprimitive) call phys_to_primitive(ixFiG^L,ixFi^L,wFi,sFi%x)

  if(slab_uniform) then
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
             sum(sFi%dvolume(ixFi^D:ixFi^D+1)*wFi(ixFi^D:ixFi^D+1,iw)) &
            /sCo%dvolume(ixCo^D)
      {end do\}
    end do
  end if

  if(stagger_grid) then
    do iw=1,nws
      ! Start one layer before
      {do ixCo^DB = ixComin^DB-kr(^DB,iw),ixComax^DB
         ixFi^DB=2*(ixCo^DB-ixComin^DB+kr(^DB,iw))+ixFimin^DB-kr(^DB,iw)\}
         ! This if statement catches the axis where surface is zero:
         if (sCo%surfaceC(ixCo^D,iw)>1.0d-9*sCo%dvolume(ixCo^D)) then ! Normal case
           wCos(ixCo^D,iw)=sum(sFi%surfaceC(ixFi^D:ixFi^D+1-kr(iw,^D),iw)*wFis(ixFi^D:ixFi^D+1-kr(iw,^D),iw)) &
                /sCo%surfaceC(ixCo^D,iw)
         else ! On axis
           wCos(ixCo^D,iw)=zero
         end if
      {end do\}
    end do
    if(phys_energy.and. .not.coarsenprimitive) then
      B_energy_change(ixCo^S)=0.5d0*sum(wCo(ixCo^S,iw_mag(:))**2,dim=ndim+1)
    end if
    ! average to fill cell-centred values
    call phys_face_to_center(ixCo^L,sCo)
    if(phys_energy.and. .not.coarsenprimitive) then
      wCo(ixCo^S,iw_e)=wCo(ixCo^S,iw_e)-B_energy_change(ixCo^S)+&
         0.5d0*sum(wCo(ixCo^S,iw_mag(:))**2,dim=ndim+1)
    end if
  end if

  if(coarsenprimitive) then
    call phys_to_conserved(ixFiG^L,ixFi^L,wFi,sFi%x)
    call phys_to_conserved(ixCoG^L,ixCo^L,wCo,sCo%x)
  end if
  end associate staggered
  end associate
end subroutine coarsen_grid
