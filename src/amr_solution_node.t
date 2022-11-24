!> Get first available igrid on processor ipe
integer function getnode(ipe)
  use mod_forest, only: igrid_inuse
  use mod_global_parameters

  integer, intent(in) :: ipe
  integer :: igrid, igrid_available

  igrid_available=0

  do igrid=1,max_blocks
     if (igrid_inuse(igrid,ipe)) cycle

     igrid_available=igrid
     exit
  end do

  if (igrid_available == 0) then
     getnode = -1
     print *, "Current maximum number of grid blocks:", max_blocks
     call mpistop("Insufficient grid blocks; increase max_blocks in meshlist")
  else
     getnode=igrid_available
     igrid_inuse(igrid,ipe)=.true.
  end if

  if (ipe==mype) then
     ! initialize nodal block
     node(1:nodehi,getnode) = 0
     rnode(1:rnodehi,getnode) = zero
  end if

end function getnode

! put igrid on processor ipe to be not in use
subroutine putnode(igrid,ipe)
  use mod_forest
  implicit none

  integer, intent(in) :: igrid, ipe

  igrid_inuse(igrid,ipe)=.false.

end subroutine putnode

!> allocate arrays on igrid node
subroutine alloc_node(igrid)
  use mod_forest
  use mod_global_parameters
  use mod_geometry
  use mod_usr_methods, only: usr_set_surface
  use mod_physics, only: phys_set_equi_vars

  integer, intent(in) :: igrid

  integer :: level, ig^D, ign^D, ixCoG^L, ix, i^D
  integer :: imin, imax, index, igCo^D, ixshift, offset, ifirst
  integer :: icase, ixGext^L
  double precision :: dx^D, summeddx, sizeuniformpart^D
  double precision :: xext(ixGlo^D-1:ixGhi^D+1,1:ndim)
  double precision :: delx_ext(ixGlo1-1:ixGhi1+1)
  double precision :: exp_factor_ext(ixGlo1-1:ixGhi1+1),del_exp_factor_ext(ixGlo1-1:ixGhi1+1),exp_factor_primitive_ext(ixGlo1-1:ixGhi1+1)
  double precision :: xc(ixGlo1:ixGhi1),delxc(ixGlo1:ixGhi1)
  double precision :: exp_factor_coarse(ixGlo1:ixGhi1),del_exp_factor_coarse(ixGlo1:ixGhi1),exp_factor_primitive_coarse(ixGlo1:ixGhi1)

  ixCoGmin^D=1;
  ixCoGmax^D=(ixGhi^D-2*nghostcells)/2+2*nghostcells;

  icase=mod(nghostcells,2)
  if(stagger_grid) icase=1
  select case(icase)
    case(0)
      ixGext^L=ixG^LL;
    case(1)
      ! for ghost cell related prolongations, we need
      ! an extra layer with known volumes and dx-intervals
      ! in case the number of ghost cells is odd
      ixGext^L=ixG^LL^LADD1;
    case default
      call mpistop("no such case")
  end select

  ! set level information
  level=igrid_to_node(igrid,mype)%node%level

  if(.not. allocated(ps(igrid)%w)) then

    ! allocate arrays for solution and space
    call alloc_state(igrid, ps(igrid), ixG^LL, ixGext^L, .true.)
    ! allocate arrays for one level coarser solution
    call alloc_state_coarse(igrid, psc(igrid), ixCoG^L, ixCoG^L)
    if(.not.convert) then
      ! allocate arrays for old solution
      call alloc_state(igrid, pso(igrid), ixG^LL, ixGext^L, .false.)
      ! allocate arrays for temp solution 1
      call alloc_state(igrid, ps1(igrid), ixG^LL, ixGext^L, .false.)

      ! allocate temporary solution space
      select case (t_integrator)
      case(ssprk3,ssprk4,IMEX_Midpoint,IMEX_Trapezoidal,IMEX_222)
        call alloc_state(igrid, ps2(igrid), ixG^LL, ixGext^L, .false.)
      case(RK3_BT,rk4,ssprk5,IMEX_CB3a)
        call alloc_state(igrid, ps2(igrid), ixG^LL, ixGext^L, .false.)
        call alloc_state(igrid, ps3(igrid), ixG^LL, ixGext^L, .false.)
      case(IMEX_ARS3,IMEX_232)
        call alloc_state(igrid, ps2(igrid), ixG^LL, ixGext^L, .false.)
        call alloc_state(igrid, ps3(igrid), ixG^LL, ixGext^L, .false.)
        call alloc_state(igrid, ps4(igrid), ixG^LL, ixGext^L, .false.)
      end select
    end if

  end if

  ! avoid dividing by zero rho in skipped corner ghostcells when phys_req_diagonal=F
  ps(igrid)%w(:^D&,1)=1.d0
  ps(igrid)%level=level
  psc(igrid)%level=level-1
  ! avoid dividing by zero rho in skipped corner ghostcells when phys_req_diagonal=F
  psc(igrid)%w(:^D&,1)=1.d0
  if(phys_trac) ps(igrid)%special_values=0.d0
  if(.not.convert) then
    pso(igrid)%level=level
    ps1(igrid)%level=level
    select case (t_integrator)
    case(ssprk3,ssprk4,IMEX_Midpoint,IMEX_Trapezoidal,IMEX_222)
      ps2(igrid)%level=level
    case(RK3_BT,rk4,ssprk5,IMEX_CB3a)
      ps2(igrid)%level=level
      ps3(igrid)%level=level
    case(IMEX_ARS3,IMEX_232)
      ps2(igrid)%level=level
      ps3(igrid)%level=level
      ps4(igrid)%level=level
    end select
  end if

  ! block pointer to current block
  block=>ps(igrid)

  ig^D=igrid_to_node(igrid,mype)%node%ig^D;

  node(plevel_,igrid)=level
  ^D&node(pig^D_,igrid)=ig^D\

  ! set dx information
  ^D&rnode(rpdx^D_,igrid)=dx(^D,level)\
  dxlevel(:)=dx(:,level)

  ! uniform cartesian case as well as all unstretched coordinates
  ! determine the minimal and maximal corners
  ^D&rnode(rpxmin^D_,igrid)=xprobmin^D+dble(ig^D-1)*dg^D(level)\
  ^D&rnode(rpxmax^D_,igrid)=xprobmin^D+dble(ig^D)*dg^D(level)\
!  ^D&rnode(rpxmax^D_,igrid)=xprobmax^D-dble(ng^D(level)-ig^D)*dg^D(level)\

  ^D&dx^D=rnode(rpdx^D_,igrid)\
 {do ix=ixGlo^D,ixMhi^D-nghostcells
    ps(igrid)%x(ix^D%ixG^T,^D)=rnode(rpxmin^D_,igrid)+(dble(ix-nghostcells)-half)*dx^D
  end do\}
 ! update overlap cells of neighboring blocks in the same way to get the same values
 {do ix=ixMhi^D-nghostcells+1,ixGhi^D
    ps(igrid)%x(ix^D%ixG^T,^D)=rnode(rpxmax^D_,igrid)+(dble(ix-ixMhi^D)-half)*dx^D
  end do\}

  ^D&dx^D=2.0d0*rnode(rpdx^D_,igrid)\
 {do ix=ixCoGmin^D,ixCoGmax^D
    psc(igrid)%x(ix^D%ixCoG^S,^D)=rnode(rpxmin^D_,igrid)+(dble(ix-nghostcells)-half)*dx^D
  end do\}

  ^D&ps(igrid)%dx(ixGext^S,^D)=rnode(rpdx^D_,igrid);
  ^D&psc(igrid)%dx(ixCoG^S,^D)=2.0d0*rnode(rpdx^D_,igrid);
  ^D&dx^D=rnode(rpdx^D_,igrid)\
 {do ix=ixGextmin^D,ixGextmax^D
    xext(ix^D%ixGext^S,^D)=rnode(rpxmin^D_,igrid)+(dble(ix-nghostcells)-half)*dx^D
  end do\}

  if(any(stretched_dim)) then
   {if(stretch_type(^D) == stretch_uni)then
      imin=(ig^D-1)*block_nx^D
      imax=ig^D*block_nx^D
      rnode(rpxmin^D_,igrid)=xprobmin^D+dxfirst_1mq(level,^D) &
                   *(1.0d0-qstretch(level,^D)**imin)
      rnode(rpxmax^D_,igrid)=xprobmin^D+dxfirst_1mq(level,^D) &
                   *(1.0d0-qstretch(level,^D)**imax)
      ! fix possible out of bound due to precision
      if(rnode(rpxmax^D_,igrid)>xprobmax^D) rnode(rpxmax^D_,igrid)=xprobmax^D
      ixshift=(ig^D-1)*block_nx^D-nghostcells
      do ix=ixGextmin^D,ixGextmax^D
        index=ixshift+ix
        ps(igrid)%dx(ix^D%ixGext^S,^D)=dxfirst(level,^D)*qstretch(level,^D)**(index-1)
      enddo
      igCo^D=(ig^D-1)/2
      ixshift=igCo^D*block_nx^D+(1-modulo(ig^D,2))*block_nx^D/2-nghostcells
      do ix=ixCoGmin^D,ixCoGmax^D
        index=ixshift+ix
        psc(igrid)%dx(ix^D%ixCoG^S,^D)=dxfirst(level-1,^D)*qstretch(level-1,^D)**(index-1)
        psc(igrid)%x(ix^D%ixCoG^S,^D)=xprobmin^D+dxfirst_1mq(level-1,^D)&
                                *(1.0d0-qstretch(level-1,^D)**(index-1))&
                    + 0.5d0*dxfirst(level-1,^D)*qstretch(level-1,^D)**(index-1)
      end do
      ! now that dx and grid boundaries are known: fill cell centers
      ifirst=nghostcells+1
      ! first fill the mesh
      summeddx=0.0d0
      do ix=ixMlo^D,ixMhi^D
        ps(igrid)%x(ix^D%ixG^T,^D)=rnode(rpxmin^D_,igrid)+summeddx+0.5d0*ps(igrid)%dx(ix^D%ixG^T,^D)
        summeddx=summeddx+ps(igrid)%dx(ix^D%ifirst,^D)
      enddo
      ! then ghost cells to left
      summeddx=0.0d0
      do ix=nghostcells,1,-1
        ps(igrid)%x(ix^D%ixG^T,^D)=rnode(rpxmin^D_,igrid)-summeddx-0.5d0*ps(igrid)%dx(ix^D%ixG^T,^D)
        summeddx=summeddx+ps(igrid)%dx(ix^D%ifirst,^D)
      enddo
      ! then ghost cells to right
      summeddx=0.0d0
      do ix=ixGhi^D-nghostcells+1,ixGhi^D
        ps(igrid)%x(ix^D%ixG^T,^D)=rnode(rpxmax^D_,igrid)+summeddx+0.5d0*ps(igrid)%dx(ix^D%ixG^T,^D)
        summeddx=summeddx+ps(igrid)%dx(ix^D%ifirst,^D)
      enddo
      select case(icase)
        case(0)
          ! if even number of ghost cells: xext is just copy of local x
          xext(ixGext^S,^D)=ps(igrid)%x(ixGext^S,^D)
        case(1)
          ! if uneven number of ghost cells: extra layer left/right
          summeddx=0.0d0
          do ix=ixMlo^D,ixMhi^D
            xext(ix^D%ixGext^S,^D)=rnode(rpxmin^D_,igrid)+summeddx+0.5d0*ps(igrid)%dx(ix^D%ixGext^S,^D)
           summeddx=summeddx+ps(igrid)%dx(ix^D%ifirst,^D)
          enddo
          ! then ghost cells to left
          summeddx=0.0d0
          do ix=nghostcells,ixGextmin^D,-1
            xext(ix^D%ixGext^S,^D)=rnode(rpxmin^D_,igrid)-summeddx-0.5d0*ps(igrid)%dx(ix^D%ixGext^S,^D)
            summeddx=summeddx+ps(igrid)%dx(ix^D%ifirst,^D)
          enddo
          ! then ghost cells to right
          summeddx=0.0d0
          do ix=ixGhi^D-nghostcells+1,ixGextmax^D
             xext(ix^D%ixGext^S,^D)=rnode(rpxmax^D_,igrid)+summeddx+0.5d0*ps(igrid)%dx(ix^D%ixGext^S,^D)
             summeddx=summeddx+ps(igrid)%dx(ix^D%ifirst,^D)
          enddo
        case default
          call mpistop("no such case")
      end select
     endif\}
    {if(stretch_type(^D) == stretch_symm)then
       ! here we distinguish three kinds of grid blocks
       ! depending on their ig-index, set per level
       !      the first n_stretchedblocks/2  will stretch to the left
       !      the middle ntotal-n_stretchedblocks will be uniform
       !      the last  n_stretchedblocks/2  will stretch to the right
       if(ig^D<=nstretchedblocks(level,^D)/2)then
         ! stretch to the left
         offset=block_nx^D*nstretchedblocks(level,^D)/2
         imin=(ig^D-1)*block_nx^D
         imax=ig^D*block_nx^D
         rnode(rpxmin^D_,igrid)=xprobmin^D+xstretch^D-dxfirst_1mq(level,^D) &
                                    *(1.0d0-qstretch(level,^D)**(offset-imin))
         rnode(rpxmax^D_,igrid)=xprobmin^D+xstretch^D-dxfirst_1mq(level,^D) &
                                    *(1.0d0-qstretch(level,^D)**(offset-imax))
         ! fix possible out of bound due to precision
         if(rnode(rpxmin^D_,igrid)<xprobmin^D) rnode(rpxmin^D_,igrid)=xprobmin^D
         ixshift=(ig^D-1)*block_nx^D-nghostcells
         do ix=ixGextmin^D,ixGextmax^D
           index=ixshift+ix
           ps(igrid)%dx(ix^D%ixGext^S,^D)=dxfirst(level,^D)*qstretch(level,^D)**(offset-index)
         enddo
         ixshift=(nstretchedblocks(level,^D)/2-ig^D)*(block_nx^D/2)+block_nx^D/2+nghostcells
         do ix=ixCoGmin^D,ixCoGmax^D
           index=ixshift-ix
           psc(igrid)%dx(ix^D%ixCoG^S,^D)=dxfirst(level-1,^D)*qstretch(level-1,^D)**index
         enddo
         ! last block: to modify ghost cells!!!
         if(ig^D==nstretchedblocks(level,^D)/2)then
           if(ng^D(level)==nstretchedblocks(level,^D))then
             ! if middle blocks do not exist then use symmetry
             do ix=ixGhi^D-nghostcells+1,ixGextmax^D
                ps(igrid)%dx(ix^D%ixGext^S,^D)= &
                ps(igrid)%dx(2*(ixGhi^D-nghostcells)+1-ix^D%ixGext^S,^D)
             enddo
             do ix=ixCoGmax^D-nghostcells+1,ixCoGmax^D
                psc(igrid)%dx(ix^D%ixCoG^S,^D)= &
                psc(igrid)%dx(2*(ixCoGmax^D-nghostcells)+1-ix^D%ixCoG^S,^D)
             enddo
           else
             ! if middle blocks exist then use same as middle blocks:
             do ix=ixGhi^D-nghostcells+1,ixGextmax^D
                ps(igrid)%dx(ix^D%ixGext^S,^D)=dxmid(level,^D)
             enddo
             do ix=ixCoGmax^D-nghostcells+1,ixCoGmax^D
                psc(igrid)%dx(ix^D%ixCoG^S,^D)=dxmid(level-1,^D)
             enddo
           endif
         endif
         ! first block: make ghost cells symmetric (to allow periodicity)
         if(ig^D==1)then
           do ix=ixGextmin^D,nghostcells
             ps(igrid)%dx(ix^D%ixGext^S,^D)=ps(igrid)%dx(2*nghostcells+1-ix^D%ixGext^S,^D)
           enddo
           do ix=1,nghostcells
             psc(igrid)%dx(ix^D%ixCoG^S,^D)=psc(igrid)%dx(2*nghostcells+1-ix^D%ixCoG^S,^D)
           enddo
         endif
       else
         if(ig^D<=ng^D(level)-nstretchedblocks(level,^D)/2) then
           ! keep uniform
           ps(igrid)%dx(ixGext^S,^D)=dxmid(level,^D)
           psc(igrid)%dx(ixCoG^S,^D)=dxmid(level-1,^D)
           rnode(rpxmin^D_,igrid)=xprobmin^D+xstretch^D+(ig^D-nstretchedblocks(level,^D)/2-1)*block_nx^D*dxmid(level,^D)
           rnode(rpxmax^D_,igrid)=xprobmin^D+xstretch^D+(ig^D-nstretchedblocks(level,^D)/2)  *block_nx^D*dxmid(level,^D)
           ! first and last block: to modify the ghost cells!!!
           if(ig^D==nstretchedblocks(level,^D)/2+1)then
             do ix=ixGextmin^D,nghostcells
               ps(igrid)%dx(ix^D%ixGext^S,^D)=dxfirst(level,^D)*qstretch(level,^D)**(nghostcells-ix)
             enddo
             do ix=1,nghostcells
               psc(igrid)%dx(ix^D%ixCoG^S,^D)=dxfirst(level-1,^D)*qstretch(level-1,^D)**(nghostcells-ix)
             enddo
           endif
           if(ig^D==ng^D(level)-nstretchedblocks(level,^D))then
             do ix=ixGhi^D-nghostcells+1,ixGextmax^D
               ps(igrid)%dx(ix^D%ixGext^S,^D)=dxfirst(level,^D)*qstretch(level,^D)**(ix-block_nx^D-nghostcells-1)
             enddo
             do ix=ixCoGmax^D-nghostcells+1,ixCoGmax^D
               psc(igrid)%dx(ix^D%ixCoG^S,^D)=dxfirst(level-1,^D)*qstretch(level-1,^D)**(ix-ixCoGmax^D+nghostcells-1)
             enddo
           endif
         else
           ! stretch to the right
           offset=block_nx^D*(ng^D(level)-nstretchedblocks(level,^D)/2)
           sizeuniformpart^D=dxmid(1,^D)*(domain_nx^D-nstretchedblocks_baselevel(^D)*block_nx^D)
           imin=(ig^D-1)*block_nx^D-offset
           imax=ig^D*block_nx^D-offset
           rnode(rpxmin^D_,igrid)=xprobmin^D+xstretch^D+sizeuniformpart^D+dxfirst_1mq(level,^D) &
                                   *(1.0d0-qstretch(level,^D)**imin)
           rnode(rpxmax^D_,igrid)=xprobmin^D+xstretch^D+sizeuniformpart^D+dxfirst_1mq(level,^D) &
                                   *(1.0d0-qstretch(level,^D)**imax)
           ! fix possible out of bound due to precision
           if(rnode(rpxmax^D_,igrid)>xprobmax^D) rnode(rpxmax^D_,igrid)=xprobmax^D
           ixshift=(ig^D-1)*block_nx^D-nghostcells-offset
           do ix=ixGextmin^D,ixGextmax^D
             index=ixshift+ix
             ps(igrid)%dx(ix^D%ixGext^S,^D)=dxfirst(level,^D)*qstretch(level,^D)**(index-1)
           enddo
           ixshift=(ig^D+nstretchedblocks(level,^D)/2-ng^D(level)-1)*(block_nx^D/2)-nghostcells
           do ix=ixCoGmin^D,ixCoGmax^D
             index=ixshift+ix
             psc(igrid)%dx(ix^D%ixCoG^S,^D)=dxfirst(level-1,^D)*qstretch(level-1,^D)**(index-1)
           enddo
           ! first block: modify ghost cells!!!
           if(ig^D==ng^D(level)-nstretchedblocks(level,^D)+1)then
             if(ng^D(level)==nstretchedblocks(level,^D))then
               ! if middle blocks do not exist then use symmetry
               do ix=ixGextmin^D,nghostcells
                 ps(igrid)%dx(ix^D%ixGext^S,^D)=ps(igrid)%dx(2*nghostcells+1-ix^D%ixGext^S,^D)
               enddo
               do ix=1,nghostcells
                 psc(igrid)%dx(ix^D%ixCoG^S,^D)=psc(igrid)%dx(2*nghostcells+1-ix^D%ixCoG^S,^D)
               enddo
             else
               ! if middle blocks exist then use same as middle blocks:
               do ix=ixGextmin^D,nghostcells
                 ps(igrid)%dx(ix^D%ixGext^S,^D)=dxmid(level,^D)
               enddo
               do ix=1,nghostcells
                 psc(igrid)%dx(ix^D%ixCoG^S,^D)=dxmid(level-1,^D)
               enddo
             endif
           endif
           ! last block: make ghost cells symmetric (to allow periodicity)
           if(ig^D==ng^D(level))then
             do ix=ixGhi^D-nghostcells+1,ixGextmax^D
               ps(igrid)%dx(ix^D%ixGext^S,^D)=ps(igrid)%dx(2*(ixGhi^D-nghostcells)+1-ix^D%ixGext^S,^D)
             enddo
             do ix=ixCoGmax^D-nghostcells+1,ixCoGmax^D
               psc(igrid)%dx(ix^D%ixCoG^S,^D)=psc(igrid)%dx(2*(ixCoGmax^D-nghostcells)+1-ix^D%ixCoG^S,^D)
             enddo
           endif
         endif
       endif
       ! now that dx and grid boundaries are known: fill cell centers
       ifirst=nghostcells+1
       ! first fill the mesh
       summeddx=0.0d0
       do ix=ixMlo^D,ixMhi^D
         ps(igrid)%x(ix^D%ixG^T,^D)=rnode(rpxmin^D_,igrid)+summeddx+0.5d0*ps(igrid)%dx(ix^D%ixG^T,^D)
         summeddx=summeddx+ps(igrid)%dx(ix^D%ifirst,^D)
       enddo
       ! then ghost cells to left
       summeddx=0.0d0
       do ix=nghostcells,1,-1
         ps(igrid)%x(ix^D%ixG^T,^D)=rnode(rpxmin^D_,igrid)-summeddx-0.5d0*ps(igrid)%dx(ix^D%ixG^T,^D)
         summeddx=summeddx+ps(igrid)%dx(ix^D%ifirst,^D)
       enddo
       ! then ghost cells to right
       summeddx=0.0d0
       do ix=ixGhi^D-nghostcells+1,ixGhi^D
         ps(igrid)%x(ix^D%ixG^T,^D)=rnode(rpxmax^D_,igrid)+summeddx+0.5d0*ps(igrid)%dx(ix^D%ixG^T,^D)
         summeddx=summeddx+ps(igrid)%dx(ix^D%ifirst,^D)
       enddo
       ! and next for the coarse representation
       ! first fill the mesh
       summeddx=0.0d0
       do ix=nghostcells+1,ixCoGmax^D-nghostcells
         psc(igrid)%x(ix^D%ixCoG^S,^D)=rnode(rpxmin^D_,igrid)+summeddx+0.5d0*psc(igrid)%dx(ix^D%ixCoG^S,^D)
         summeddx=summeddx+psc(igrid)%dx(ix^D%ifirst,^D)
       enddo
       ! then ghost cells to left
       summeddx=0.0d0
       do ix=nghostcells,1,-1
         psc(igrid)%x(ix^D%ixCoG^S,^D)=rnode(rpxmin^D_,igrid)-summeddx-0.5d0*psc(igrid)%dx(ix^D%ixCoG^S,^D)
         summeddx=summeddx+psc(igrid)%dx(ix^D%ifirst,^D)
       enddo
       ! then ghost cells to right
       summeddx=0.0d0
       do ix=ixCoGmax^D-nghostcells+1,ixCoGmax^D
         psc(igrid)%x(ix^D%ixCoG^S,^D)=rnode(rpxmax^D_,igrid)+summeddx+0.5d0*psc(igrid)%dx(ix^D%ixCoG^S,^D)
         summeddx=summeddx+psc(igrid)%dx(ix^D%ifirst,^D)
       enddo
       select case(icase)
         case(0)
           ! if even number of ghost cells: xext is just copy of local x
           xext(ixGext^S,^D)=ps(igrid)%x(ixGext^S,^D)
         case(1)
           ! if uneven number of ghost cells: extra layer left/right
           summeddx=0.0d0
           do ix=ixMlo^D,ixMhi^D
             xext(ix^D%ixGext^S,^D)=rnode(rpxmin^D_,igrid)+summeddx+0.5d0*ps(igrid)%dx(ix^D%ixGext^S,^D)
            summeddx=summeddx+ps(igrid)%dx(ix^D%ifirst,^D)
           enddo
           ! then ghost cells to left
           summeddx=0.0d0
           do ix=nghostcells,ixGextmin^D,-1
             xext(ix^D%ixGext^S,^D)=rnode(rpxmin^D_,igrid)-summeddx-0.5d0*ps(igrid)%dx(ix^D%ixGext^S,^D)
             summeddx=summeddx+ps(igrid)%dx(ix^D%ifirst,^D)
           enddo
          ! then ghost cells to right
          summeddx=0.0d0
          do ix=ixGhi^D-nghostcells+1,ixGextmax^D
             xext(ix^D%ixGext^S,^D)=rnode(rpxmax^D_,igrid)+summeddx+0.5d0*ps(igrid)%dx(ix^D%ixGext^S,^D)
             summeddx=summeddx+ps(igrid)%dx(ix^D%ifirst,^D)
          enddo
         case default
           call mpistop("no such case")
       end select
     endif\}
  endif

  ! calculate area of cell surfaces for standard block
  call get_surface_area(ps(igrid),ixG^LL)
  ! calculate area of cell surfaces for coarser representative block
  call get_surface_area(psc(igrid),ixCoG^L)
  ! calculate volume and distance of cells
  ps(igrid)%dsC=1.d0
  select case (coordinate)
    case (Cartesian)
      ps(igrid)%dvolume(ixGext^S)= {^D&rnode(rpdx^D_,igrid)|*}
      ps(igrid)%ds(ixGext^S,1:ndim)=ps(igrid)%dx(ixGext^S,1:ndim)
      ps(igrid)%dsC(ixGext^S,1:ndim)=ps(igrid)%dx(ixGext^S,1:ndim)
      psc(igrid)%dvolume(ixCoG^S)= {^D&2.d0*rnode(rpdx^D_,igrid)|*}
      psc(igrid)%ds(ixCoG^S,1:ndim)=psc(igrid)%dx(ixCoG^S,1:ndim)
    case (Cartesian_stretched)
      ps(igrid)%dvolume(ixGext^S)= {^D&ps(igrid)%dx(ixGext^S,^D)|*}
      ps(igrid)%ds(ixGext^S,1:ndim)=ps(igrid)%dx(ixGext^S,1:ndim)
      ps(igrid)%dsC(ixGext^S,1:ndim)=ps(igrid)%dx(ixGext^S,1:ndim)
      psc(igrid)%dvolume(ixCoG^S)= {^D&psc(igrid)%dx(ixCoG^S,^D)|*}
      psc(igrid)%ds(ixCoG^S,1:ndim)=psc(igrid)%dx(ixCoG^S,1:ndim)
    case (Cartesian_expansion)
      {^IFONED
      delx_ext(ixGext^S) = ps(igrid)%dx(ixGext^S,1)
      if(associated(usr_set_surface)) call usr_set_surface(ixGext^L,xext(ixGext^S,1),delx_ext(ixGext^S),exp_factor_ext(ixGext^S),del_exp_factor_ext(ixGext^S),exp_factor_primitive_ext(ixGext^S))
      ps(igrid)%dvolume(ixGext^S)= exp_factor_primitive_ext(ixGext^S)
      ps(igrid)%ds(ixGext^S,1)=ps(igrid)%dx(ixGext^S,1)
      ps(igrid)%dsC(ixGext^S,1)=ps(igrid)%dx(ixGext^S,1)
      xc(ixCoG^S) = psc(igrid)%x(ixCoG^S,1)
      delxc(ixCoG^S) = psc(igrid)%dx(ixCoG^S,1)
      if(associated(usr_set_surface)) call usr_set_surface(ixCoG^L,xc(ixCoG^S),delxc(ixCoG^S),exp_factor_coarse(ixCoG^S),del_exp_factor_coarse(ixCoG^S),exp_factor_primitive_coarse(ixCoG^S))
      psc(igrid)%dvolume(ixCoG^S)= exp_factor_primitive_coarse(ixCoG^S)
      psc(igrid)%ds(ixCoG^S,1)=psc(igrid)%dx(ixCoG^S,1)
      }
    case (spherical)
      ps(igrid)%dvolume(ixGext^S)=(xext(ixGext^S,1)**2 &
                                +ps(igrid)%dx(ixGext^S,1)**2/12.0d0)*&
              ps(igrid)%dx(ixGext^S,1){^NOONED &
             *two*dabs(dsin(xext(ixGext^S,2))) &
             *dsin(half*ps(igrid)%dx(ixGext^S,2))}{^IFTHREED*ps(igrid)%dx(ixGext^S,3)}
      psc(igrid)%dvolume(ixCoG^S)=(psc(igrid)%x(ixCoG^S,1)**2 &
                                       +psc(igrid)%dx(ixCoG^S,1)**2/12.0d0)*&
              psc(igrid)%dx(ixCoG^S,1){^NOONED &
             *two*dabs(dsin(psc(igrid)%x(ixCoG^S,2))) &
             *dsin(half*psc(igrid)%dx(ixCoG^S,2))}{^IFTHREED*psc(igrid)%dx(ixCoG^S,3)}
      ps(igrid)%ds(ixGext^S,1)=ps(igrid)%dx(ixGext^S,1)
      {^NOONED   ps(igrid)%ds(ixGext^S,2)=xext(ixGext^S,1)*ps(igrid)%dx(ixGext^S,2)}
      {^IFTHREED ps(igrid)%ds(ixGext^S,3)=xext(ixGext^S,1)*dsin(xext(ixGext^S,2))*&
                                          ps(igrid)%dx(ixGext^S,3)}
      ps(igrid)%dsC(ixGext^S,1)=ps(igrid)%dx(ixGext^S,1)
      {^NOONED   ps(igrid)%dsC(ixGext^S,2)=(xext(ixGext^S,1)+half*ps(igrid)%dx(ixGext^S,1))*&
                                          ps(igrid)%dx(ixGext^S,2)
      if(ndir>ndim) then
        ps(igrid)%dsC(ixGext^S,3)=(xext(ixGext^S,1)+half*ps(igrid)%dx(ixGext^S,1))*&
                                       dsin(xext(ixGext^S,2)+half*ps(igrid)%dx(ixGext^S,2))
      end if
      }
      {^IFTHREED ps(igrid)%dsC(ixGext^S,3)=(xext(ixGext^S,1)+half*ps(igrid)%dx(ixGext^S,1))*&
                                       dsin(xext(ixGext^S,2)+half*ps(igrid)%dx(ixGext^S,2))*&
                                       ps(igrid)%dx(ixGext^S,3)}
    case (cylindrical)
      ps(igrid)%dvolume(ixGext^S)=dabs(xext(ixGext^S,1)) &
           *ps(igrid)%dx(ixGext^S,1){^DE&*ps(igrid)%dx(ixGext^S,^DE) }
      psc(igrid)%dvolume(ixCoG^S)=dabs(psc(igrid)%x(ixCoG^S,1)) &
           *psc(igrid)%dx(ixCoG^S,1){^DE&*psc(igrid)%dx(ixCoG^S,^DE) }
      ps(igrid)%ds(ixGext^S,r_)=ps(igrid)%dx(ixGext^S,r_)
      ps(igrid)%dsC(ixGext^S,r_)=ps(igrid)%dx(ixGext^S,r_)
      if(z_>0.and.z_<=ndim) then
        ps(igrid)%ds(ixGext^S,z_)=ps(igrid)%dx(ixGext^S,z_)
        ps(igrid)%dsC(ixGext^S,z_)=ps(igrid)%dx(ixGext^S,z_)
        if(phi_>z_.and.ndir>ndim) then
          ps(igrid)%dsC(ixGext^S,phi_)=xext(ixGext^S,1)+half*ps(igrid)%dx(ixGext^S,1)
        end if
      end if
      if(phi_>0.and.phi_<=ndim) then
        ps(igrid)%ds(ixGext^S,phi_)=xext(ixGext^S,1)*ps(igrid)%dx(ixGext^S,phi_)
        ps(igrid)%dsC(ixGext^S,phi_)=(xext(ixGext^S,1)+&
                   half*ps(igrid)%dx(ixGext^S,1))*ps(igrid)%dx(ixGext^S,phi_)
        if(z_>phi_.and.ndir>ndim) ps(igrid)%dsC(ixGext^S,z_)=1.d0
      end if
    case default
      call mpistop("Sorry, coordinate unknown")
  end select

  ! initialize background non-evolving solution
  if (B0field) call set_B0_grid(igrid)
  if (number_equi_vars>0) call phys_set_equi_vars(igrid)

  ! find the blocks on the boundaries
  ps(igrid)%is_physical_boundary=.false.
  {
  do i^D=-1,1
    if(i^D==0) cycle
    ign^D=ig^D+i^D
    ! blocks at periodic boundary have neighbors in the physical domain
    ! thus threated at internal blocks with no physical boundary
    if (periodB(^D)) ign^D=1+modulo(ign^D-1,ng^D(level))
    if (ign^D > ng^D(level)) then
       if(phi_ > 0 .and. poleB(2,^D)) then
         ! if at a pole, the boundary is not physical boundary
         ps(igrid)%is_physical_boundary(2*^D)=.false.
       else
         ps(igrid)%is_physical_boundary(2*^D)=.true.
       end if
    else if (ign^D < 1) then
       if(phi_ > 0 .and. poleB(1,^D)) then
         ! if at a pole, the boundary is not physical boundary
         ps(igrid)%is_physical_boundary(2*^D-1)=.false.
       else
         ps(igrid)%is_physical_boundary(2*^D-1)=.true.
       end if
    end if
  end do
  \}
  if(any(ps(igrid)%is_physical_boundary)) then
    phyboundblock(igrid)=.true.
  else
    phyboundblock(igrid)=.false.
  end if

end subroutine alloc_node

!> allocate memory to physical state of igrid node
subroutine alloc_state(igrid, s, ixG^L, ixGext^L, alloc_once_for_ps)
  use mod_global_parameters
  type(state) :: s
  integer, intent(in) :: igrid, ixG^L, ixGext^L
  logical, intent(in) :: alloc_once_for_ps
  integer             :: ixGs^L

  allocate(s%w(ixG^S,1:nw))
  s%igrid=igrid
  s%w=0.d0
  s%ixG^L=ixG^L;
  {^D& ixGsmin^D = ixGmin^D-1; ixGsmax^D = ixGmax^D|;}
  if(stagger_grid) then
    allocate(s%ws(ixGs^S,1:nws))
    s%ws=0.d0
    if(record_electric_field) then
      allocate(s%we(ixGs^S,1:nws))
      s%we=0.d0
    end if
    s%ixGs^L=ixGs^L;
  end if
  if(alloc_once_for_ps) then
    ! allocate extra variables for ps state
    if(nw_extra>0) allocate(s%wextra(ixG^S,1:nw_extra))
    ! allocate coordinates
    allocate(s%x(ixG^S,1:ndim))
    allocate(s%dx(ixGext^S,1:ndim), &
               s%ds(ixGext^S,1:ndim),s%dsC(ixGext^S,1:3))
    allocate(s%dvolume(ixGext^S))
    allocate(s%surfaceC(ixGs^S,1:ndim), &
             s%surface(ixG^S,1:ndim))
    ! allocate physical boundary flag
    allocate(s%is_physical_boundary(2*ndim))
    if(B0field) then
      allocate(s%B0(ixG^S,1:ndir,0:ndim))
      allocate(s%J0(ixG^S,7-2*ndir:3))
    end if
    if(number_equi_vars > 0) then
      allocate(s%equi_vars(ixG^S,1:number_equi_vars,0:ndim))
    endif

    ! allocate space for special values for each block state
    if(phys_trac) then
      ! special_values(1) Tcoff local
      ! special_values(2) Tmax local
      ! special_values(3:2+ndim) Bdir local
      allocate(s%special_values(ndim+2))
    end if
  else
    ! share common info from ps states to save memory
    if(nw_extra>0) s%wextra=>ps(igrid)%wextra
    s%x=>ps(igrid)%x
    s%dx=>ps(igrid)%dx
    s%ds=>ps(igrid)%ds
    s%dsC=>ps(igrid)%dsC
    s%dvolume=>ps(igrid)%dvolume
    s%surfaceC=>ps(igrid)%surfaceC
    s%surface=>ps(igrid)%surface
    s%is_physical_boundary=>ps(igrid)%is_physical_boundary
    if(B0field) then
      s%B0=>ps(igrid)%B0
      s%J0=>ps(igrid)%J0
    end if
    if(number_equi_vars > 0) then
      s%equi_vars=>ps(igrid)%equi_vars
    endif
    if(phys_trac) s%special_values=>ps(igrid)%special_values
  end if
end subroutine alloc_state

!> allocate memory to one-level coarser physical state of igrid node
subroutine alloc_state_coarse(igrid, s, ixG^L, ixGext^L)
  use mod_global_parameters 
  use mod_mhd_phys, only: mhd_semirelativistic
  type(state) :: s
  integer, intent(in) :: igrid, ixG^L, ixGext^L
  integer             :: ixGs^L

  allocate(s%w(ixG^S,1:nw))
  s%igrid=igrid
  s%w=0.d0
  s%ixG^L=ixG^L;
  {^D& ixGsmin^D = ixGmin^D-1; ixGsmax^D = ixGmax^D|;}
  if(stagger_grid) then
    allocate(s%ws(ixGs^S,1:nws))
    s%ws=0.d0
    s%ixGs^L=ixGs^L;
  end if
  if(B0field.and.mhd_semirelativistic) then
    allocate(s%B0(ixG^S,1:ndir,0:ndim))
  end if
  ! allocate coordinates
  allocate(s%x(ixG^S,1:ndim))
  allocate(s%dx(ixGext^S,1:ndim), &
             s%ds(ixGext^S,1:ndim),s%dsC(ixGext^S,1:3))
  allocate(s%dvolume(ixGext^S))
  allocate(s%surfaceC(ixGs^S,1:ndim), &
           s%surface(ixG^S,1:ndim))
  ! allocate physical boundary flag
  allocate(s%is_physical_boundary(2*ndim))
end subroutine alloc_state_coarse

subroutine dealloc_state(igrid, s,dealloc_x)
  use mod_global_parameters
  integer, intent(in) :: igrid
  type(state) :: s
  logical, intent(in) :: dealloc_x

  deallocate(s%w)
  if(stagger_grid) then
    deallocate(s%ws)
  end if
  if(dealloc_x) then
    if(nw_extra>0) deallocate(s%wextra)
    ! deallocate coordinates
    deallocate(s%x)
    deallocate(s%dx,s%ds,s%dsC)
    deallocate(s%dvolume)
    deallocate(s%surfaceC,s%surface)
    deallocate(s%is_physical_boundary)
    if(B0field) then
      deallocate(s%B0)
      deallocate(s%J0)
    end if
    if(number_equi_vars > 0) then
      deallocate(s%equi_vars)
    end if
  else
    nullify(s%x,s%dx,s%ds,s%dsC,s%dvolume,s%surfaceC,s%surface)
    nullify(s%is_physical_boundary)
    if(B0field) nullify(s%B0,s%J0)
    if(number_equi_vars > 0) then
      nullify(s%equi_vars)
    end if
    if(nw_extra>0) nullify(s%wextra)
  end if
end subroutine dealloc_state

subroutine dealloc_state_coarse(igrid, s)
  use mod_global_parameters
  use mod_mhd_phys, only: mhd_semirelativistic
  integer, intent(in) :: igrid
  type(state) :: s

  deallocate(s%w)
  if(stagger_grid) then
    deallocate(s%ws)
  end if
  if(B0field.and.mhd_semirelativistic) then
    deallocate(s%B0)
  end if
  ! deallocate coordinates
  deallocate(s%x)
  deallocate(s%dx,s%ds,s%dsC)
  deallocate(s%dvolume)
  deallocate(s%surfaceC,s%surface)
  deallocate(s%is_physical_boundary)
end subroutine dealloc_state_coarse

subroutine dealloc_node(igrid)
  use mod_global_parameters

  integer, intent(in) :: igrid

  if (igrid==0) then
     call mpistop("trying to delete a non-existing grid in dealloc_node")
  end if

  call dealloc_state(igrid, ps(igrid),.true.)
  call dealloc_state_coarse(igrid, psc(igrid))
  if(.not.convert) then
    call dealloc_state(igrid, ps1(igrid),.false.)
    call dealloc_state(igrid, pso(igrid),.false.)
    ! deallocate temporary solution space
    select case (t_integrator)
    case(ssprk3,ssprk4,IMEX_Midpoint,IMEX_Trapezoidal,IMEX_222)
      call dealloc_state(igrid, ps2(igrid),.false.)
    case(RK3_BT,rk4,ssprk5,IMEX_CB3a)
      call dealloc_state(igrid, ps2(igrid),.false.)
      call dealloc_state(igrid, ps3(igrid),.false.)
    case(IMEX_ARS3,IMEX_232)
      call dealloc_state(igrid, ps2(igrid),.false.)
      call dealloc_state(igrid, ps3(igrid),.false.)
      call dealloc_state(igrid, ps4(igrid),.false.)
    end select
  end if

end subroutine dealloc_node
