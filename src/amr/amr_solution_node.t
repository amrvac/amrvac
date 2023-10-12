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

  integer, intent(in) :: igrid

  integer :: level, ig^D, ign^D, ixCoG^L, ix, i^D
  integer :: imin, imax, index, igCo^D, ixshift, offset, ifirst
  integer :: icase, ixGext^L
  double precision :: dx^D, summeddx, sizeuniformpart^D
  double precision :: xext(ixGlo^D-1:ixGhi^D+1,1:ndim)

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

  if(.not. allocated(mesh(igrid)%x)) then
     call alloc_mesh(igrid, level, mesh(igrid), ixG^LL, ixGext^L)
     call alloc_mesh(igrid, level-1, mesh_co(igrid), ixCoG^L, ixCoG^L)
  end if

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
  {if(rnode(rpxmax^D_,igrid)>xprobmax^D) rnode(rpxmax^D_,igrid)=xprobmax^D\}

  ^D&dx^D=rnode(rpdx^D_,igrid)\
 {do ix=ixGlo^D,ixMhi^D-nghostcells
    mesh(igrid)%x(ix^D%ixG^T,^D)=rnode(rpxmin^D_,igrid)+(dble(ix-nghostcells)-0.5d0)*dx^D
  end do\}
 ! update overlap cells of neighboring blocks in the same way to get the same values
 {do ix=ixMhi^D-nghostcells+1,ixGhi^D
    mesh(igrid)%x(ix^D%ixG^T,^D)=rnode(rpxmax^D_,igrid)+(dble(ix-ixMhi^D)-0.5d0)*dx^D
  end do\}

  ^D&dx^D=2.0d0*rnode(rpdx^D_,igrid)\
 {do ix=ixCoGmin^D,ixCoGmax^D
    mesh_co(igrid)%x(ix^D%ixCoG^S,^D)=rnode(rpxmin^D_,igrid)+(dble(ix-nghostcells)-0.5d0)*dx^D
  end do\}

  ^D&mesh(igrid)%dx(ixGext^S,^D)=rnode(rpdx^D_,igrid);
  ^D&mesh_co(igrid)%dx(ixCoG^S,^D)=2.0d0*rnode(rpdx^D_,igrid);
  ^D&dx^D=rnode(rpdx^D_,igrid)\
 {do ix=ixGextmin^D,ixGextmax^D
    xext(ix^D%ixGext^S,^D)=rnode(rpxmin^D_,igrid)+(dble(ix-nghostcells)-0.5d0)*dx^D
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
        mesh(igrid)%dx(ix^D%ixGext^S,^D)=dxfirst(level,^D)*qstretch(level,^D)**(index-1)
      enddo
      igCo^D=(ig^D-1)/2
      ixshift=igCo^D*block_nx^D+(1-modulo(ig^D,2))*block_nx^D/2-nghostcells
      do ix=ixCoGmin^D,ixCoGmax^D
        index=ixshift+ix
        mesh_co(igrid)%dx(ix^D%ixCoG^S,^D)=dxfirst(level-1,^D)*qstretch(level-1,^D)**(index-1)
        mesh_co(igrid)%x(ix^D%ixCoG^S,^D)=xprobmin^D+dxfirst_1mq(level-1,^D)&
                                *(1.0d0-qstretch(level-1,^D)**(index-1))&
                    + 0.5d0*dxfirst(level-1,^D)*qstretch(level-1,^D)**(index-1)
      end do
      ! now that dx and grid boundaries are known: fill cell centers
      ifirst=nghostcells+1
      ! first fill the mesh
      summeddx=0.0d0
      do ix=ixMlo^D,ixMhi^D
        mesh(igrid)%x(ix^D%ixG^T,^D)=rnode(rpxmin^D_,igrid)+summeddx+0.5d0*mesh(igrid)%dx(ix^D%ixG^T,^D)
        summeddx=summeddx+mesh(igrid)%dx(ix^D%ifirst,^D)
      enddo
      ! then ghost cells to left
      summeddx=0.0d0
      do ix=nghostcells,1,-1
        mesh(igrid)%x(ix^D%ixG^T,^D)=rnode(rpxmin^D_,igrid)-summeddx-0.5d0*mesh(igrid)%dx(ix^D%ixG^T,^D)
        summeddx=summeddx+mesh(igrid)%dx(ix^D%ifirst,^D)
      enddo
      ! then ghost cells to right
      summeddx=0.0d0
      do ix=ixGhi^D-nghostcells+1,ixGhi^D
        mesh(igrid)%x(ix^D%ixG^T,^D)=rnode(rpxmax^D_,igrid)+summeddx+0.5d0*mesh(igrid)%dx(ix^D%ixG^T,^D)
        summeddx=summeddx+mesh(igrid)%dx(ix^D%ifirst,^D)
      enddo
      select case(icase)
        case(0)
          ! if even number of ghost cells: xext is just copy of local x
          xext(ixGext^S,^D)=mesh(igrid)%x(ixGext^S,^D)
        case(1)
          ! if uneven number of ghost cells: extra layer left/right
          summeddx=0.0d0
          do ix=ixMlo^D,ixMhi^D
            xext(ix^D%ixGext^S,^D)=rnode(rpxmin^D_,igrid)+summeddx+0.5d0*mesh(igrid)%dx(ix^D%ixGext^S,^D)
           summeddx=summeddx+mesh(igrid)%dx(ix^D%ifirst,^D)
          enddo
          ! then ghost cells to left
          summeddx=0.0d0
          do ix=nghostcells,ixGextmin^D,-1
            xext(ix^D%ixGext^S,^D)=rnode(rpxmin^D_,igrid)-summeddx-0.5d0*mesh(igrid)%dx(ix^D%ixGext^S,^D)
            summeddx=summeddx+mesh(igrid)%dx(ix^D%ifirst,^D)
          enddo
          ! then ghost cells to right
          summeddx=0.0d0
          do ix=ixGhi^D-nghostcells+1,ixGextmax^D
             xext(ix^D%ixGext^S,^D)=rnode(rpxmax^D_,igrid)+summeddx+0.5d0*mesh(igrid)%dx(ix^D%ixGext^S,^D)
             summeddx=summeddx+mesh(igrid)%dx(ix^D%ifirst,^D)
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
           mesh(igrid)%dx(ix^D%ixGext^S,^D)=dxfirst(level,^D)*qstretch(level,^D)**(offset-index)
         enddo
         ixshift=(nstretchedblocks(level,^D)/2-ig^D)*(block_nx^D/2)+block_nx^D/2+nghostcells
         do ix=ixCoGmin^D,ixCoGmax^D
           index=ixshift-ix
           mesh_co(igrid)%dx(ix^D%ixCoG^S,^D)=dxfirst(level-1,^D)*qstretch(level-1,^D)**index
         enddo
         ! last block: to modify ghost cells!!!
         if(ig^D==nstretchedblocks(level,^D)/2)then
           if(ng^D(level)==nstretchedblocks(level,^D))then
             ! if middle blocks do not exist then use symmetry
             do ix=ixGhi^D-nghostcells+1,ixGextmax^D
                mesh(igrid)%dx(ix^D%ixGext^S,^D)= &
                mesh(igrid)%dx(2*(ixGhi^D-nghostcells)+1-ix^D%ixGext^S,^D)
             enddo
             do ix=ixCoGmax^D-nghostcells+1,ixCoGmax^D
                mesh_co(igrid)%dx(ix^D%ixCoG^S,^D)= &
                mesh_co(igrid)%dx(2*(ixCoGmax^D-nghostcells)+1-ix^D%ixCoG^S,^D)
             enddo
           else
             ! if middle blocks exist then use same as middle blocks:
             do ix=ixGhi^D-nghostcells+1,ixGextmax^D
                mesh(igrid)%dx(ix^D%ixGext^S,^D)=dxmid(level,^D)
             enddo
             do ix=ixCoGmax^D-nghostcells+1,ixCoGmax^D
                mesh_co(igrid)%dx(ix^D%ixCoG^S,^D)=dxmid(level-1,^D)
             enddo
           endif
         endif
         ! first block: make ghost cells symmetric (to allow periodicity)
         if(ig^D==1)then
           do ix=ixGextmin^D,nghostcells
             mesh(igrid)%dx(ix^D%ixGext^S,^D)=mesh(igrid)%dx(2*nghostcells+1-ix^D%ixGext^S,^D)
           enddo
           do ix=1,nghostcells
             mesh_co(igrid)%dx(ix^D%ixCoG^S,^D)=mesh_co(igrid)%dx(2*nghostcells+1-ix^D%ixCoG^S,^D)
           enddo
         endif
       else
         if(ig^D<=ng^D(level)-nstretchedblocks(level,^D)/2) then
           ! keep uniform
           mesh(igrid)%dx(ixGext^S,^D)=dxmid(level,^D)
           mesh_co(igrid)%dx(ixCoG^S,^D)=dxmid(level-1,^D)
           rnode(rpxmin^D_,igrid)=xprobmin^D+xstretch^D+(ig^D-nstretchedblocks(level,^D)/2-1)*block_nx^D*dxmid(level,^D)
           rnode(rpxmax^D_,igrid)=xprobmin^D+xstretch^D+(ig^D-nstretchedblocks(level,^D)/2)  *block_nx^D*dxmid(level,^D)
           ! first and last block: to modify the ghost cells!!!
           if(ig^D==nstretchedblocks(level,^D)/2+1)then
             do ix=ixGextmin^D,nghostcells
               mesh(igrid)%dx(ix^D%ixGext^S,^D)=dxfirst(level,^D)*qstretch(level,^D)**(nghostcells-ix)
             enddo
             do ix=1,nghostcells
               mesh_co(igrid)%dx(ix^D%ixCoG^S,^D)=dxfirst(level-1,^D)*qstretch(level-1,^D)**(nghostcells-ix)
             enddo
           endif
           if(ig^D==ng^D(level)-nstretchedblocks(level,^D))then
             do ix=ixGhi^D-nghostcells+1,ixGextmax^D
               mesh(igrid)%dx(ix^D%ixGext^S,^D)=dxfirst(level,^D)*qstretch(level,^D)**(ix-block_nx^D-nghostcells-1)
             enddo
             do ix=ixCoGmax^D-nghostcells+1,ixCoGmax^D
               mesh_co(igrid)%dx(ix^D%ixCoG^S,^D)=dxfirst(level-1,^D)*qstretch(level-1,^D)**(ix-ixCoGmax^D+nghostcells-1)
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
             mesh(igrid)%dx(ix^D%ixGext^S,^D)=dxfirst(level,^D)*qstretch(level,^D)**(index-1)
           enddo
           ixshift=(ig^D+nstretchedblocks(level,^D)/2-ng^D(level)-1)*(block_nx^D/2)-nghostcells
           do ix=ixCoGmin^D,ixCoGmax^D
             index=ixshift+ix
             mesh_co(igrid)%dx(ix^D%ixCoG^S,^D)=dxfirst(level-1,^D)*qstretch(level-1,^D)**(index-1)
           enddo
           ! first block: modify ghost cells!!!
           if(ig^D==ng^D(level)-nstretchedblocks(level,^D)+1)then
             if(ng^D(level)==nstretchedblocks(level,^D))then
               ! if middle blocks do not exist then use symmetry
               do ix=ixGextmin^D,nghostcells
                 mesh(igrid)%dx(ix^D%ixGext^S,^D)=mesh(igrid)%dx(2*nghostcells+1-ix^D%ixGext^S,^D)
               enddo
               do ix=1,nghostcells
                 mesh_co(igrid)%dx(ix^D%ixCoG^S,^D)=mesh_co(igrid)%dx(2*nghostcells+1-ix^D%ixCoG^S,^D)
               enddo
             else
               ! if middle blocks exist then use same as middle blocks:
               do ix=ixGextmin^D,nghostcells
                 mesh(igrid)%dx(ix^D%ixGext^S,^D)=dxmid(level,^D)
               enddo
               do ix=1,nghostcells
                 mesh_co(igrid)%dx(ix^D%ixCoG^S,^D)=dxmid(level-1,^D)
               enddo
             endif
           endif
           ! last block: make ghost cells symmetric (to allow periodicity)
           if(ig^D==ng^D(level))then
             do ix=ixGhi^D-nghostcells+1,ixGextmax^D
               mesh(igrid)%dx(ix^D%ixGext^S,^D)=mesh(igrid)%dx(2*(ixGhi^D-nghostcells)+1-ix^D%ixGext^S,^D)
             enddo
             do ix=ixCoGmax^D-nghostcells+1,ixCoGmax^D
               mesh_co(igrid)%dx(ix^D%ixCoG^S,^D)=mesh_co(igrid)%dx(2*(ixCoGmax^D-nghostcells)+1-ix^D%ixCoG^S,^D)
             enddo
           endif
         endif
       endif
       ! now that dx and grid boundaries are known: fill cell centers
       ifirst=nghostcells+1
       ! first fill the mesh
       summeddx=0.0d0
       do ix=ixMlo^D,ixMhi^D
         mesh(igrid)%x(ix^D%ixG^T,^D)=rnode(rpxmin^D_,igrid)+summeddx+0.5d0*mesh(igrid)%dx(ix^D%ixG^T,^D)
         summeddx=summeddx+mesh(igrid)%dx(ix^D%ifirst,^D)
       enddo
       ! then ghost cells to left
       summeddx=0.0d0
       do ix=nghostcells,1,-1
         mesh(igrid)%x(ix^D%ixG^T,^D)=rnode(rpxmin^D_,igrid)-summeddx-0.5d0*mesh(igrid)%dx(ix^D%ixG^T,^D)
         summeddx=summeddx+mesh(igrid)%dx(ix^D%ifirst,^D)
       enddo
       ! then ghost cells to right
       summeddx=0.0d0
       do ix=ixGhi^D-nghostcells+1,ixGhi^D
         mesh(igrid)%x(ix^D%ixG^T,^D)=rnode(rpxmax^D_,igrid)+summeddx+0.5d0*mesh(igrid)%dx(ix^D%ixG^T,^D)
         summeddx=summeddx+mesh(igrid)%dx(ix^D%ifirst,^D)
       enddo
       ! and next for the coarse representation
       ! first fill the mesh
       summeddx=0.0d0
       do ix=nghostcells+1,ixCoGmax^D-nghostcells
         mesh_co(igrid)%x(ix^D%ixCoG^S,^D)=rnode(rpxmin^D_,igrid)+summeddx+0.5d0*mesh_co(igrid)%dx(ix^D%ixCoG^S,^D)
         summeddx=summeddx+mesh_co(igrid)%dx(ix^D%ifirst,^D)
       enddo
       ! then ghost cells to left
       summeddx=0.0d0
       do ix=nghostcells,1,-1
         mesh_co(igrid)%x(ix^D%ixCoG^S,^D)=rnode(rpxmin^D_,igrid)-summeddx-0.5d0*mesh_co(igrid)%dx(ix^D%ixCoG^S,^D)
         summeddx=summeddx+mesh_co(igrid)%dx(ix^D%ifirst,^D)
       enddo
       ! then ghost cells to right
       summeddx=0.0d0
       do ix=ixCoGmax^D-nghostcells+1,ixCoGmax^D
         mesh_co(igrid)%x(ix^D%ixCoG^S,^D)=rnode(rpxmax^D_,igrid)+summeddx+0.5d0*mesh_co(igrid)%dx(ix^D%ixCoG^S,^D)
         summeddx=summeddx+mesh_co(igrid)%dx(ix^D%ifirst,^D)
       enddo
       select case(icase)
         case(0)
           ! if even number of ghost cells: xext is just copy of local x
           xext(ixGext^S,^D)=mesh(igrid)%x(ixGext^S,^D)
         case(1)
           ! if uneven number of ghost cells: extra layer left/right
           summeddx=0.0d0
           do ix=ixMlo^D,ixMhi^D
             xext(ix^D%ixGext^S,^D)=rnode(rpxmin^D_,igrid)+summeddx+0.5d0*mesh(igrid)%dx(ix^D%ixGext^S,^D)
            summeddx=summeddx+mesh(igrid)%dx(ix^D%ifirst,^D)
           enddo
           ! then ghost cells to left
           summeddx=0.0d0
           do ix=nghostcells,ixGextmin^D,-1
             xext(ix^D%ixGext^S,^D)=rnode(rpxmin^D_,igrid)-summeddx-0.5d0*mesh(igrid)%dx(ix^D%ixGext^S,^D)
             summeddx=summeddx+mesh(igrid)%dx(ix^D%ifirst,^D)
           enddo
          ! then ghost cells to right
          summeddx=0.0d0
          do ix=ixGhi^D-nghostcells+1,ixGextmax^D
             xext(ix^D%ixGext^S,^D)=rnode(rpxmax^D_,igrid)+summeddx+0.5d0*mesh(igrid)%dx(ix^D%ixGext^S,^D)
             summeddx=summeddx+mesh(igrid)%dx(ix^D%ifirst,^D)
          enddo
         case default
           call mpistop("no such case")
       end select
     endif\}
  endif

  ! calculate barycenter for standard block
  call get_xbar(mesh(igrid),ixG^LL)
  ! calculate area of cell surfaces for standard block
  call get_surface_area(mesh(igrid),ixG^LL)
  ! calculate area of cell surfaces for coarser representative block
  call get_surface_area(mesh_co(igrid),ixCoG^L)
  ! calculate volume and distance of cells
  mesh(igrid)%dsC=1.d0
  mesh(igrid)%dlE=1.d0
  select case (coordinate)
    case (Cartesian)
      mesh(igrid)%dvolume(ixGext^S)= {^D&rnode(rpdx^D_,igrid)|*}
      mesh(igrid)%ds(ixGext^S,1:ndim)=mesh(igrid)%dx(ixGext^S,1:ndim)
      mesh(igrid)%dsC(ixGext^S,1:ndim)=mesh(igrid)%dx(ixGext^S,1:ndim)
      mesh(igrid)%dlE(ixGext^S,1:ndim)=mesh(igrid)%dx(ixGext^S,1:ndim)
      mesh_co(igrid)%dvolume(ixCoG^S)= {^D&2.d0*rnode(rpdx^D_,igrid)|*}
      mesh_co(igrid)%ds(ixCoG^S,1:ndim)=mesh_co(igrid)%dx(ixCoG^S,1:ndim)
    case (Cartesian_stretched)
      mesh(igrid)%dvolume(ixGext^S)= {^D&mesh(igrid)%dx(ixGext^S,^D)|*}
      mesh(igrid)%ds(ixGext^S,1:ndim)=mesh(igrid)%dx(ixGext^S,1:ndim)
      mesh(igrid)%dsC(ixGext^S,1:ndim)=mesh(igrid)%dx(ixGext^S,1:ndim)
      mesh(igrid)%dlE(ixGext^S,1:ndim)=mesh(igrid)%dx(ixGext^S,1:ndim)
      mesh_co(igrid)%dvolume(ixCoG^S)= {^D&mesh_co(igrid)%dx(ixCoG^S,^D)|*}
      mesh_co(igrid)%ds(ixCoG^S,1:ndim)=mesh_co(igrid)%dx(ixCoG^S,1:ndim)
    case (spherical)
      mesh(igrid)%dvolume(ixGext^S)=(xext(ixGext^S,1)**2 &
                                +mesh(igrid)%dx(ixGext^S,1)**2/12.0d0)*&
              mesh(igrid)%dx(ixGext^S,1){^NOONED &
             *2.0d0*dabs(dsin(xext(ixGext^S,2))) &
             *dsin(0.5d0*mesh(igrid)%dx(ixGext^S,2))}{^IFTHREED*mesh(igrid)%dx(ixGext^S,3)}
      mesh_co(igrid)%dvolume(ixCoG^S)=(mesh_co(igrid)%x(ixCoG^S,1)**2 &
                                       +mesh_co(igrid)%dx(ixCoG^S,1)**2/12.0d0)*&
              mesh_co(igrid)%dx(ixCoG^S,1){^NOONED &
             *2.0d0*dabs(dsin(mesh_co(igrid)%x(ixCoG^S,2))) &
             *dsin(0.5d0*mesh_co(igrid)%dx(ixCoG^S,2))}{^IFTHREED*mesh_co(igrid)%dx(ixCoG^S,3)}
      mesh(igrid)%ds(ixGext^S,1)=mesh(igrid)%dx(ixGext^S,1)
      {^NOONED   mesh(igrid)%ds(ixGext^S,2)=xext(ixGext^S,1)*mesh(igrid)%dx(ixGext^S,2)}
      {^IFTHREED mesh(igrid)%ds(ixGext^S,3)=xext(ixGext^S,1)*dsin(xext(ixGext^S,2))*&
                                          mesh(igrid)%dx(ixGext^S,3)}
      mesh(igrid)%dsC(ixGext^S,1)=mesh(igrid)%dx(ixGext^S,1)
      {^NOONED   mesh(igrid)%dsC(ixGext^S,2)=(xext(ixGext^S,1)+0.5d0*mesh(igrid)%dx(ixGext^S,1))*&
                                          mesh(igrid)%dx(ixGext^S,2)
      if(ndir>ndim) then
        mesh(igrid)%dsC(ixGext^S,3)=(xext(ixGext^S,1)+0.5d0*mesh(igrid)%dx(ixGext^S,1))*&
                                       dsin(xext(ixGext^S,2)+0.5d0*mesh(igrid)%dx(ixGext^S,2))
      end if
      }
      {^IFTHREED mesh(igrid)%dsC(ixGext^S,3)=(xext(ixGext^S,1)+0.5d0*mesh(igrid)%dx(ixGext^S,1))*&
                                       dsin(xext(ixGext^S,2)+0.5d0*mesh(igrid)%dx(ixGext^S,2))*&
                                       mesh(igrid)%dx(ixGext^S,3)}
      mesh(igrid)%dlE(ixGext^S,1)=(xext(ixGext^S,1)**2+mesh(igrid)%dx(ixGext^S,1)**2/12.0d0)*&
              mesh(igrid)%dx(ixGext^S,1) {^NOONED * dabs(dsin(xext(ixGext^S,2)+0.5d0*mesh(igrid)%dx(ixGext^S,2))) }
      {^NOONED   mesh(igrid)%dlE(ixGext^S,2)=(xext(ixGext^S,1)+0.5d0*mesh(igrid)%dx(ixGext^S,1))**2 * &
                      2.0d0*dabs(dsin(xext(ixGext^S,2)))*dsin(0.5d0*mesh(igrid)%dx(ixGext^S,2)) }
      {^IFTHREED mesh(igrid)%dlE(ixGext^S,3)=(xext(ixGext^S,1)+0.5d0*mesh(igrid)%dx(ixGext^S,1))**2 *&
                                       dsin(xext(ixGext^S,2)+0.5d0*mesh(igrid)%dx(ixGext^S,2))*&
                                       mesh(igrid)%dx(ixGext^S,3)}
    case (cylindrical)
      mesh(igrid)%dvolume(ixGext^S)=dabs(xext(ixGext^S,1)) &
           *mesh(igrid)%dx(ixGext^S,1){^DE&*mesh(igrid)%dx(ixGext^S,^DE) }
      mesh_co(igrid)%dvolume(ixCoG^S)=dabs(mesh_co(igrid)%x(ixCoG^S,1)) &
           *mesh_co(igrid)%dx(ixCoG^S,1){^DE&*mesh_co(igrid)%dx(ixCoG^S,^DE) }
      mesh(igrid)%ds(ixGext^S,r_)=mesh(igrid)%dx(ixGext^S,r_)
      mesh(igrid)%dsC(ixGext^S,r_)=mesh(igrid)%dx(ixGext^S,r_)
      mesh(igrid)%dlE(ixGext^S,r_)=(xext(ixGext^S,1)+&
                   0.5d0*mesh(igrid)%dx(ixGext^S,1))*mesh(igrid)%dx(ixGext^S,r_)
      if(z_>0.and.z_<=ndim) then
        mesh(igrid)%ds(ixGext^S,z_)=mesh(igrid)%dx(ixGext^S,z_)
        mesh(igrid)%dsC(ixGext^S,z_)=mesh(igrid)%dx(ixGext^S,z_)
        mesh(igrid)%dlE(ixGext^S,z_)=(xext(ixGext^S,1)+&
                   0.5d0*mesh(igrid)%dx(ixGext^S,1))*mesh(igrid)%dx(ixGext^S,z_)
        if(phi_>z_.and.ndir>ndim) then
          mesh(igrid)%dsC(ixGext^S,phi_)= 1.0d0 !xext(ixGext^S,1)+0.5d0*mesh(igrid)%dx(ixGext^S,1)
        end if
      end if
      if(phi_>0.and.phi_<=ndim) then
        mesh(igrid)%ds(ixGext^S,phi_)=xext(ixGext^S,1)*mesh(igrid)%dx(ixGext^S,phi_)
        mesh(igrid)%dsC(ixGext^S,phi_)=(xext(ixGext^S,1)+&
                   0.5d0*mesh(igrid)%dx(ixGext^S,1))*mesh(igrid)%dx(ixGext^S,phi_)
        mesh(igrid)%dlE(ixGext^S,phi_)=(xext(ixGext^S,1)+&
                   0.5d0*mesh(igrid)%dx(ixGext^S,1))*mesh(igrid)%dx(ixGext^S,phi_)
        if(z_>phi_.and.ndir>ndim) then
           mesh(igrid)%dsC(ixGext^S,z_)=1.d0
           mesh(igrid)%dlE(ixGext^S,z_)=xext(ixGext^S,1)+0.5d0*mesh(igrid)%dx(ixGext^S,1)
        end if
      end if
    case default
      call mpistop("Sorry, coordinate unknown")
  end select

  if ( coordinate /= cartesian ) &
     call get_christoffel(mesh(igrid),ixG^LL)

  ! find the blocks on the boundaries
  mesh(igrid)%is_physical_boundary=.false.
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
         mesh(igrid)%is_physical_boundary(2*^D)=.false.
       else
         mesh(igrid)%is_physical_boundary(2*^D)=.true.
       end if
    else if (ign^D < 1) then
       if(phi_ > 0 .and. poleB(1,^D)) then
         ! if at a pole, the boundary is not physical boundary
         mesh(igrid)%is_physical_boundary(2*^D-1)=.false.
       else
         mesh(igrid)%is_physical_boundary(2*^D-1)=.true.
       end if
    end if
  end do
  \}
  if(any(mesh(igrid)%is_physical_boundary)) then
    phyboundblock(igrid)=.true.
  else
    phyboundblock(igrid)=.false.
  end if

  if(.not. allocated(metric(igrid)%vars)) then
    ! allocate arrays for solution and space
    call alloc_metric(igrid, metric(igrid), ixG^LL, ixGext^L)
    ! allocate arrays for one level coarser solution
    call alloc_metric(igrid, metric_co(igrid), ixCoG^L, ixCoG^L)
  end if

  if(.not. allocated(ps(igrid)%w)) then

    ! allocate arrays for solution and space
    call alloc_state(igrid, ps(igrid), ixG^LL, ixGext^L, .false.)
    ! allocate arrays for one level coarser solution
    call alloc_state(igrid, psc(igrid), ixCoG^L, ixCoG^L, .true.)
    if(.not.convert) then
      ! allocate arrays for old solution
      call alloc_state(igrid, pso(igrid), ixG^LL, ixGext^L, .false.)
      ! allocate arrays for temp solution 1
      call alloc_state(igrid, ps1(igrid), ixG^LL, ixGext^L, .false.)

      ! allocate temperary solution space
      select case (t_integrator)
      case(ssprk3,ssprk4,jameson,IMEX_Midpoint,IMEX_Trapezoidal,IMEX_222)
        call alloc_state(igrid, ps2(igrid), ixG^LL, ixGext^L, .false.)
      case(RK3_BT,rk4,ssprk5,IMEX_CB3a,IMEX_RK4)
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
  ps(igrid)%w  = smalldouble
  ps1(igrid)%w = smalldouble
  psc(igrid)%w = smalldouble
  ps(igrid)%w  = smalldouble
  ps1(igrid)%w = smalldouble
  psc(igrid)%w = smalldouble

!  ps(igrid)%level=level
!  psc(igrid)%level=level-1
!  if(.not.convert) then
!    pso(igrid)%level=level
!    ps1(igrid)%level=level
!    select case (t_integrator)
!    case(ssprk3,ssprk4,jameson,IMEX_Midpoint,IMEX_Trapezoidal,IMEX_222)
!      ps2(igrid)%level=level
!    case(RK3_BT,rk4,ssprk5,IMEX_CB3a)
!      ps2(igrid)%level=level
!      ps3(igrid)%level=level
!    case(IMEX_ARS3,IMEX_232)
!      ps2(igrid)%level=level
!      ps3(igrid)%level=level
!      ps4(igrid)%level=level
!    end select
!  end if

  ! block pointer to current block
  block=>ps(igrid)

end subroutine alloc_node

!> allocate memory to physical state of igrid node
subroutine alloc_mesh(igrid, level, s, ixG^L, ixGext^L)
  use mod_global_parameters
  use mod_geometry
  type(mesh_t)        :: s
  integer, intent(in) :: igrid, level, ixG^L, ixGext^L
  integer             :: ixGs^L

  s%igrid=igrid
  s%level=level
  s%ixG^L=ixG^L;
  {^D& ixGsmin^D = ixGmin^D-1; ixGsmax^D = ixGmax^D|;}
  s%ixGs^L=ixGs^L;

  ! allocate coordinates
  allocate(s%x(ixG^S,1:ndim))
  allocate(s%xbar(ixG^S,1:ndim))
  allocate(s%dx(ixGext^S,1:ndim), &
             s%ds(ixGext^S,1:ndim),s%dsC(ixGext^S,1:3),s%dlE(ixGext^S,1:3))
  allocate(s%dvolume(ixGext^S))
  allocate(s%surfaceC(ixGs^S,1:ndim), &
           s%surface(ixG^S,1:ndim))
  select case (coordinate)
    case (spherical)
      allocate(s%christoffel(ixGext^S,0:3 {^NOONED +2 }))
    case (cylindrical)
      allocate(s%christoffel(ixGext^S,0:2))
  end select
  ! allocate physical boundary flag
  allocate(s%is_physical_boundary(2*ndim))
end subroutine alloc_mesh

!> allocate memory to metric of igrid node
subroutine alloc_metric(igrid, s, ixG^L, ixGext^L)
  use mod_global_parameters
  use mod_geometry
  type(metric_t)       :: s
  integer, intent(in)  :: igrid, ixG^L, ixGext^L

  allocate(s%vars(ixG^S,1:nmetric))
  s%vars = 0.d0
  s%vars(ixG^S,1:2) = 1.d0
end subroutine alloc_metric

!> allocate memory to physical state of igrid node
subroutine alloc_state(igrid, s, ixG^L, ixGext^L, is_coarse)
  use mod_global_parameters
  use mod_geometry
  type(state)          :: s
  integer, intent(in)  :: igrid, ixG^L, ixGext^L
  logical, intent(in)  :: is_coarse
  integer              :: ixGs^L

  s%is_prim = .False.
  allocate(s%w(ixG^S,1:nw))
  s%w=0.d0
  {^D& ixGsmin^D = ixGmin^D-1; ixGsmax^D = ixGmax^D|;}
  if(stagger_grid) then
    allocate(s%ws(ixGs^S,1:nws))
    s%ws=0.d0
  end if
  if (nw_extra>0) then
    allocate(s%wextra(ixG^S,1:nw_extra))
    s%wextra=0.d0
  end if

  if (is_coarse) then
    s%mesh   => mesh_co(igrid)
    s%metric => metric_co(igrid)
  else
    s%mesh   => mesh(igrid)
    s%metric => metric(igrid)
  end if
end subroutine alloc_state

subroutine dealloc_mesh(s)
  use mod_global_parameters
  type(mesh_t)        :: s
  ! deallocate coordinates
  deallocate(s%x)
  deallocate(s%xbar)
  deallocate(s%dx,s%ds,s%dsC,s%dlE)
  deallocate(s%dvolume)
  deallocate(s%surfaceC,s%surface)
  if ( coordinate /= cartesian ) &
     deallocate(s%christoffel)
  deallocate(s%is_physical_boundary)
end subroutine dealloc_mesh

subroutine dealloc_metric(s)
  use mod_global_parameters
  type(metric_t) :: s
  deallocate(s%vars)
end subroutine dealloc_metric

subroutine dealloc_state(s)
  use mod_global_parameters
  type(state) :: s
  deallocate(s%w)
  if(stagger_grid) then
    deallocate(s%ws)
  end if
  if (nw_extra>0) then
    deallocate(s%wextra)
  end if
  nullify(s%metric)
  nullify(s%mesh)
end subroutine dealloc_state

subroutine dealloc_node(igrid)
  use mod_global_parameters

  integer, intent(in) :: igrid

  ! fixme: no psc?
  if (igrid==0) then
     call mpistop("trying to delete a non-existing grid in dealloc_node")
  end if

  call dealloc_state(ps(igrid))
  call dealloc_state(ps1(igrid))
  call dealloc_state(pso(igrid))
  ! deallocate temporary solution space
  select case (t_integrator)
  case(ssprk3,ssprk4,jameson,IMEX_Midpoint,IMEX_Trapezoidal,IMEX_222)
    call dealloc_state(ps2(igrid))
  case(RK3_BT,rk4,ssprk5,IMEX_CB3a)
    call dealloc_state(ps2(igrid))
    call dealloc_state(ps3(igrid))
  case(IMEX_ARS3,IMEX_232)
    call dealloc_state(ps2(igrid))
    call dealloc_state(ps3(igrid))
    call dealloc_state(ps4(igrid))
  end select

  call dealloc_metric(metric(igrid))

  call dealloc_mesh(mesh(igrid))

end subroutine dealloc_node
