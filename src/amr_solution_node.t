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
!=============================================================================
subroutine putnode(igrid,ipe)
use mod_forest
implicit none

! putnode = return igrid node on processor ipe
 
integer, intent(in) :: igrid, ipe
!----------------------------------------------------------------------------
igrid_inuse(igrid,ipe)=.false.

end subroutine putnode
!=============================================================================
subroutine alloc_node(igrid)
use mod_forest
use mod_global_parameters

integer, intent(in) :: igrid

integer :: level, ig^D, ign^D, ixCoG^L, ix, i^D
integer :: imin, imax, index, igCo^D, ixshift, offset, ifirst
integer:: icase, ixGext^L
double precision :: rXmin^D, dx^D, summeddx, sizeuniformpart^D
double precision :: xext(ixGlo^D-1:ixGhi^D+1,1:ndim)
!-----------------------------------------------------------------------------
ixCoGmin^D=1;
!ixCoGmax^D=ixGhi^D/2+nghostcells;
ixCoGmax^D=(ixGhi^D-2*nghostcells)/2+2*nghostcells;

icase=mod(nghostcells,2)
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

if(.not. allocated(pw(igrid)%w)) then
  
  ! initialize solution space
  allocate(pw(igrid)%w(ixG^T,1:nw), &
           pw(igrid)%wold(ixG^T,1:nw), &
           pw(igrid)%w1(ixG^T,1:nw), &
           pw(igrid)%wcoarse(ixCoG^S,1:nw))
  
  ! wio for visualization data
  allocate(pw(igrid)%wio(ixG^T,1:nw+nwauxio))
  
  ! allocate temperary solution space
  select case (time_integrator)
  case("threestep","fourstep","jameson","twostep_trapezoidal")
    allocate(pw(igrid)%w2(ixG^T,1:nw))
  case("rk4","ssprk43")
    allocate(pw(igrid)%w2(ixG^T,1:nw))
    allocate(pw(igrid)%w3(ixG^T,1:nw))
  case("ssprk54")
    allocate(pw(igrid)%w2(ixG^T,1:nw))
    allocate(pw(igrid)%w3(ixG^T,1:nw))
    allocate(pw(igrid)%w4(ixG^T,1:nw))
  end select
  
  ! allocate coordinates
  allocate(pw(igrid)%x(ixG^T,1:ndim), &
           pw(igrid)%xcoarse(ixCoG^S,1:ndim))
  if(.not.slab)then
      allocate(pw(igrid)%dx(ixGext^S,1:ndim), &
               pw(igrid)%dxcoarse(ixCoG^S,1:ndim),&
               pw(igrid)%ds(ixGext^S,1:ndim))
      allocate(pw(igrid)%dvolume(ixGext^S), &
               pw(igrid)%dvolumecoarse(ixCoG^S))
      allocate(pw(igrid)%surfaceC(ixG^T,1:ndim), &
               pw(igrid)%surface(ixG^T,1:ndim))
  endif
end if

pw(igrid)%w=0.d0
pw(igrid)%wold=0.d0
pw(igrid)%w1=0.d0
pw(igrid)%wcoarse=0.d0
pw(igrid)%igrid=igrid
pw(igrid)%wio=0.d0

! wb is w by default
pw(igrid)%wb=>pw(igrid)%w


! block pointer to current block
block=>pw(igrid)

if(phys_energy .and. solve_internal_e) then
  block%e_is_internal=.true.
endif

! set level information
level=igrid_to_node(igrid,mype)%node%level
ig^D=igrid_to_node(igrid,mype)%node%ig^D;

node(plevel_,igrid)=level
^D&node(pig^D_,igrid)=ig^D\

! set dx information
^D&rnode(rpdx^D_,igrid)=dx(^D,level)\
dxlevel(:)=dx(:,level)

! uniform cartesian case as well as all unstretched coordinates
! determine the minimal and maximal corners
^D&rnode(rpxmin^D_,igrid)=xprobmin^D+dble(ig^D-1)*dg^D(level)\
^D&rnode(rpxmax^D_,igrid)=xprobmax^D-dble(ng^D(level)-ig^D)*dg^D(level)\

^D&dx^D=rnode(rpdx^D_,igrid)\
^D&rXmin^D=rnode(rpxmin^D_,igrid)-nghostcells*dx^D\
{do ix=ixGlo^D,ixGhi^D
   pw(igrid)%x(ix^D%ixG^T,^D)=rXmin^D+(dble(ix)-half)*dx^D
end do\}
^D&dx^D=2.0d0*rnode(rpdx^D_,igrid)\
^D&rXmin^D=rnode(rpxmin^D_,igrid)-nghostcells*dx^D\
{do ix=ixCoGmin^D,ixCoGmax^D
   pw(igrid)%xcoarse(ix^D%ixCoG^S,^D)=rXmin^D+(dble(ix)-half)*dx^D
end do\}

if (.not.slab) then
  ^D&pw(igrid)%dx(ixGext^S,^D)=rnode(rpdx^D_,igrid);
  ^D&pw(igrid)%dxcoarse(ixCoG^S,^D)=2.0d0*rnode(rpdx^D_,igrid);
  ^D&dx^D=rnode(rpdx^D_,igrid)\
  ^D&rXmin^D=rnode(rpxmin^D_,igrid)-nghostcells*dx^D\
  {do ix=ixGextmin^D,ixGextmax^D
      xext(ix^D%ixGext^S,^D)=rXmin^D+(dble(ix)-half)*dx^D
  end do\}
endif

if(stretched_grid) then
  {if(stretched_dim(^D))then
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
      pw(igrid)%dx(ix^D%ixG^T,^D)=dxfirst(level,^D)*qstretch(level,^D)**(index-1)
    enddo
    igCo^D=(ig^D-1)/2
    ixshift=igCo^D*block_nx^D+(1-modulo(ig^D,2))*block_nx^D/2-nghostcells
    do ix=ixCoGmin^D,ixCoGmax^D
      index=ixshift+ix
      pw(igrid)%dxcoarse(ix^D%ixCoG^S,^D)=dxfirst(level-1,^D)*qstretch(level-1,^D)**(index-1)
      pw(igrid)%xcoarse(ix^D%ixCoG^S,^D)=xprobmin^D+dxfirst_1mq(level-1,^D)&
                                         *(1.0d0-qstretch(level-1,^D)**(index-1)) &
                  + 0.5d0*dxfirst(level-1,^D)*qstretch(level-1,^D)**(index-1)
    end do
    ! now that dx and grid boundaries are known: fill cell centers
    ifirst=nghostcells+1
    ! first fill the mesh
    summeddx=0.0d0
    do ix=ixMlo^D,ixMhi^D
       pw(igrid)%x(ix^D%ixG^T,^D)=rnode(rpxmin^D_,igrid)+summeddx+0.5d0*pw(igrid)%dx(ix^D%ixG^T,^D)
       summeddx=summeddx+pw(igrid)%dx(ix^D%ifirst,^D)
    enddo
    ! then ghost cells to left
    summeddx=0.0d0
    do ix=nghostcells,1,-1
       pw(igrid)%x(ix^D%ixG^T,^D)=rnode(rpxmin^D_,igrid)-summeddx-0.5d0*pw(igrid)%dx(ix^D%ixG^T,^D)
       summeddx=summeddx+pw(igrid)%dx(ix^D%ifirst,^D)
    enddo
    ! then ghost cells to right
    summeddx=0.0d0
    do ix=ixGhi^D-nghostcells+1,ixGhi^D
       pw(igrid)%x(ix^D%ixG^T,^D)=rnode(rpxmax^D_,igrid)+summeddx+0.5d0*pw(igrid)%dx(ix^D%ixG^T,^D)
       summeddx=summeddx+pw(igrid)%dx(ix^D%ifirst,^D)
    enddo
    select case(icase)
      case(0)
        ! if even number of ghost cells: xext is just copy of local x
        xext(ixGext^S,^D)=pw(igrid)%x(ixGext^S,^D)
      case(1)
        ! if uneven number of ghost cells: extra layer left/right
        summeddx=0.0d0
        do ix=ixMlo^D,ixMhi^D
          xext(ix^D%ixGext^S,^D)=rnode(rpxmin^D_,igrid)+summeddx+0.5d0*pw(igrid)%dx(ix^D%ixGext^S,^D)
         summeddx=summeddx+pw(igrid)%dx(ix^D%ifirst,^D)
        enddo
        ! then ghost cells to left
        summeddx=0.0d0
        do ix=nghostcells,ixGextmin^D,-1
          xext(ix^D%ixGext^S,^D)=rnode(rpxmin^D_,igrid)-summeddx-0.5d0*pw(igrid)%dx(ix^D%ixGext^S,^D)
          summeddx=summeddx+pw(igrid)%dx(ix^D%ifirst,^D)
        enddo
       ! then ghost cells to right
       summeddx=0.0d0
       do ix=ixGhi^D-nghostcells+1,ixGextmax^D
          xext(ix^D%ixGext^S,^D)=rnode(rpxmax^D_,igrid)+summeddx+0.5d0*pw(igrid)%dx(ix^D%ixGext^S,^D)
          summeddx=summeddx+pw(igrid)%dx(ix^D%ifirst,^D)
       enddo
      case default
        call mpistop("no such case")
    end select
   endif\}
  {if(stretched_symm_dim(^D))then
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
         pw(igrid)%dx(ix^D%ixG^T,^D)=dxfirst(level,^D)*qstretch(level,^D)**(offset-index)
      enddo
      ixshift=(nstretchedblocks(level,^D)/2-ig^D)*(block_nx^D/2)+block_nx^D/2+nghostcells
      do ix=ixCoGmin^D,ixCoGmax^D
         index=ixshift-ix
         pw(igrid)%dxcoarse(ix^D%ixCoG^S,^D)=dxfirst(level-1,^D)*qstretch(level-1,^D)**index
      enddo
      ! last block: to modify ghost cells!!!
      if(ig^D==nstretchedblocks(level,^D)/2)then
        if(ng^D(level)==nstretchedblocks(level,^D))then
           ! if middle blocks do not exist then use symmetry
           do ix=ixGhi^D-nghostcells+1,ixGextmax^D
              pw(igrid)%dx(ix^D%ixG^T,^D)= &
              pw(igrid)%dx(2*(ixGhi^D-nghostcells)+1-ix^D%ixG^T,^D)
           enddo
           do ix=ixCoGmax^D-nghostcells+1,ixCoGmax^D
              pw(igrid)%dxcoarse(ix^D%ixCoG^S,^D)= &
              pw(igrid)%dxcoarse(2*(ixCoGmax^D-nghostcells)+1-ix^D%ixCoG^S,^D)
           enddo
        else
           ! if middle blocks exist then use same as middle blocks: 
           do ix=ixGhi^D-nghostcells+1,ixGextmax^D
              pw(igrid)%dx(ix^D%ixG^T,^D)=dxmid(level,^D)
           enddo
           do ix=ixCoGmax^D-nghostcells+1,ixCoGmax^D
              pw(igrid)%dxcoarse(ix^D%ixCoG^S,^D)=dxmid(level-1,^D)
           enddo
        endif
      endif
      ! first block: make ghost cells symmetric (to allow periodicity)
      if(ig^D==1)then
         do ix=ixGextmin^D,nghostcells
            pw(igrid)%dx(ix^D%ixGext^S,^D)=pw(igrid)%dx(2*nghostcells+1-ix^D%ixGext^S,^D)
         enddo
         do ix=1,nghostcells
            pw(igrid)%dxcoarse(ix^D%ixCoG^S,^D)=pw(igrid)%dxcoarse(2*nghostcells+1-ix^D%ixCoG^S,^D)
         enddo
      endif
    else 
      if (ig^D<=ng^D(level)-nstretchedblocks(level,^D)/2) then
         ! keep uniform
         pw(igrid)%dx(ixGext^S,^D)=dxmid(level,^D)
         pw(igrid)%dxcoarse(ixCoG^S,^D)=dxmid(level-1,^D)
         rnode(rpxmin^D_,igrid)=xprobmin^D+xstretch^D+(ig^D-nstretchedblocks(level,^D)/2-1)*block_nx^D*dxmid(level,^D)
         rnode(rpxmax^D_,igrid)=xprobmin^D+xstretch^D+(ig^D-nstretchedblocks(level,^D)/2)  *block_nx^D*dxmid(level,^D)
         ! first and last block: to modify the ghost cells!!!
         if(ig^D==nstretchedblocks(level,^D)/2+1)then
            do ix=ixGextmin^D,nghostcells
               pw(igrid)%dx(ix^D%ixGext^S,^D)=dxfirst(level,^D)*qstretch(level,^D)**(nghostcells-ix)
            enddo
            do ix=1,nghostcells
               pw(igrid)%dxcoarse(ix^D%ixCoG^S,^D)=dxfirst(level-1,^D)*qstretch(level-1,^D)**(nghostcells-ix)
            enddo
         endif
         if(ig^D==ng^D(level)-nstretchedblocks(level,^D))then
            do ix=ixGhi^D-nghostcells+1,ixGextmax^D
              pw(igrid)%dx(ix^D%ixG^T,^D)=dxfirst(level,^D)*qstretch(level,^D)**(ix-block_nx^D-nghostcells-1)
            enddo
            do ix=ixCoGmax^D-nghostcells+1,ixCoGmax^D
              pw(igrid)%dxcoarse(ix^D%ixCoG^S,^D)=dxfirst(level-1,^D)*qstretch(level-1,^D)**(ix-ixCoGmax^D+nghostcells-1)
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
            pw(igrid)%dx(ix^D%ixGext^S,^D)=dxfirst(level,^D)*qstretch(level,^D)**(index-1)
         enddo
         ixshift=(ig^D+nstretchedblocks(level,^D)/2-ng^D(level)-1)*(block_nx^D/2)-nghostcells
         do ix=ixCoGmin^D,ixCoGmax^D
            index=ixshift+ix
            pw(igrid)%dxcoarse(ix^D%ixCoG^S,^D)=dxfirst(level-1,^D)*qstretch(level-1,^D)**(index-1)
         enddo
         ! first block: modify ghost cells!!!
         if(ig^D==ng^D(level)-nstretchedblocks(level,^D)+1)then
            if(ng^D(level)==nstretchedblocks(level,^D))then
               ! if middle blocks do not exist then use symmetry
               do ix=ixGextmin^D,nghostcells
                  pw(igrid)%dx(ix^D%ixGext^S,^D)=pw(igrid)%dx(2*nghostcells+1-ix^D%ixGext^S,^D)
               enddo
               do ix=1,nghostcells
                  pw(igrid)%dxcoarse(ix^D%ixCoG^S,^D)=pw(igrid)%dxcoarse(2*nghostcells+1-ix^D%ixCoG^S,^D)
               enddo
            else
               ! if middle blocks exist then use same as middle blocks: 
               do ix=ixGextmin^D,nghostcells
                  pw(igrid)%dx(ix^D%ixGext^S,^D)=dxmid(level,^D)
               enddo
               do ix=1,nghostcells
                  pw(igrid)%dxcoarse(ix^D%ixCoG^S,^D)=dxmid(level-1,^D)
               enddo
            endif
         endif
         ! last block: make ghost cells symmetric (to allow periodicity)
         if(ig^D==ng^D(level))then
            do ix=ixGhi^D-nghostcells+1,ixGextmax^D
               pw(igrid)%dx(ix^D%ixGext^S,^D)=pw(igrid)%dx(2*(ixGhi^D-nghostcells)+1-ix^D%ixGext^S,^D)
            enddo
           do ix=ixCoGmax^D-nghostcells+1,ixCoGmax^D
              pw(igrid)%dxcoarse(ix^D%ixCoG^S,^D)=pw(igrid)%dxcoarse(2*(ixCoGmax^D-nghostcells)+1-ix^D%ixCoG^S,^D)
           enddo
         endif
      endif
    endif
    ! now that dx and grid boundaries are known: fill cell centers
    ifirst=nghostcells+1
    ! first fill the mesh
    summeddx=0.0d0
    do ix=ixMlo^D,ixMhi^D
       pw(igrid)%x(ix^D%ixG^T,^D)=rnode(rpxmin^D_,igrid)+summeddx+0.5d0*pw(igrid)%dx(ix^D%ixG^T,^D)
       summeddx=summeddx+pw(igrid)%dx(ix^D%ifirst,^D)
    enddo
    ! then ghost cells to left
    summeddx=0.0d0
    do ix=nghostcells,1,-1
       pw(igrid)%x(ix^D%ixG^T,^D)=rnode(rpxmin^D_,igrid)-summeddx-0.5d0*pw(igrid)%dx(ix^D%ixG^T,^D)
       summeddx=summeddx+pw(igrid)%dx(ix^D%ifirst,^D)
    enddo
    ! then ghost cells to right
    summeddx=0.0d0
    do ix=ixGhi^D-nghostcells+1,ixGhi^D
       pw(igrid)%x(ix^D%ixG^T,^D)=rnode(rpxmax^D_,igrid)+summeddx+0.5d0*pw(igrid)%dx(ix^D%ixG^T,^D)
       summeddx=summeddx+pw(igrid)%dx(ix^D%ifirst,^D)
    enddo
    ! and next for the coarse representation
    ! first fill the mesh
    summeddx=0.0d0
    do ix=nghostcells+1,ixCoGmax^D-nghostcells
       pw(igrid)%xcoarse(ix^D%ixCoG^S,^D)=rnode(rpxmin^D_,igrid)+summeddx+0.5d0*pw(igrid)%dxcoarse(ix^D%ixCoG^S,^D)
       summeddx=summeddx+pw(igrid)%dxcoarse(ix^D%ifirst,^D)
    enddo
    ! then ghost cells to left
    summeddx=0.0d0
    do ix=nghostcells,1,-1
       pw(igrid)%xcoarse(ix^D%ixCoG^S,^D)=rnode(rpxmin^D_,igrid)-summeddx-0.5d0*pw(igrid)%dxcoarse(ix^D%ixCoG^S,^D)
       summeddx=summeddx+pw(igrid)%dxcoarse(ix^D%ifirst,^D)
    enddo
    ! then ghost cells to right
    summeddx=0.0d0
    do ix=ixCoGmax^D-nghostcells+1,ixCoGmax^D
       pw(igrid)%xcoarse(ix^D%ixCoG^S,^D)=rnode(rpxmax^D_,igrid)+summeddx+0.5d0*pw(igrid)%dxcoarse(ix^D%ixCoG^S,^D)
       summeddx=summeddx+pw(igrid)%dxcoarse(ix^D%ifirst,^D)
    enddo
    select case(icase)
      case(0)
        ! if even number of ghost cells: xext is just copy of local x
        xext(ixGext^S,^D)=pw(igrid)%x(ixGext^S,^D)
      case(1)
        ! if uneven number of ghost cells: extra layer left/right
        summeddx=0.0d0
        do ix=ixMlo^D,ixMhi^D
          xext(ix^D%ixGext^S,^D)=rnode(rpxmin^D_,igrid)+summeddx+0.5d0*pw(igrid)%dx(ix^D%ixGext^S,^D)
         summeddx=summeddx+pw(igrid)%dx(ix^D%ifirst,^D)
        enddo
        ! then ghost cells to left
        summeddx=0.0d0
        do ix=nghostcells,ixGextmin^D,-1
          xext(ix^D%ixGext^S,^D)=rnode(rpxmin^D_,igrid)-summeddx-0.5d0*pw(igrid)%dx(ix^D%ixGext^S,^D)
          summeddx=summeddx+pw(igrid)%dx(ix^D%ifirst,^D)
        enddo
       ! then ghost cells to right
       summeddx=0.0d0
       do ix=ixGhi^D-nghostcells+1,ixGextmax^D
          xext(ix^D%ixGext^S,^D)=rnode(rpxmax^D_,igrid)+summeddx+0.5d0*pw(igrid)%dx(ix^D%ixGext^S,^D)
          summeddx=summeddx+pw(igrid)%dx(ix^D%ifirst,^D)
       enddo
      case default
        call mpistop("no such case")
    end select
   endif\}
endif

if (.not.slab) then
   call fillgeo(igrid,ixG^LL)
   select case (typeaxial)
      case ("slabstretch")
         pw(igrid)%dvolume(ixGext^S)= {^D&pw(igrid)%dx(ixGext^S,^D)|*}
         pw(igrid)%dvolumecoarse(ixCoG^S)= {^D&pw(igrid)%dxcoarse(ixCoG^S,^D)|*}
         pw(igrid)%ds(ixGext^S,1:ndim)=pw(igrid)%dx(ixGext^S,1:ndim)
      case ("spherical")
         pw(igrid)%dvolume(ixGext^S)=(xext(ixGext^S,1)**2 &
                                   +pw(igrid)%dx(ixGext^S,1)**2/12.0d0)*&
                 pw(igrid)%dx(ixGext^S,1){^NOONED &
                *two*dabs(dsin(xext(ixGext^S,2))) &
                *dsin(half*pw(igrid)%dx(ixGext^S,2))}{^IFTHREED*pw(igrid)%dx(ixGext^S,3)}
         pw(igrid)%dvolumecoarse(ixCoG^S)=(pw(igrid)%xcoarse(ixCoG^S,1)**2 &
                                          +pw(igrid)%dxcoarse(ixCoG^S,1)**2/12.0d0)*&
                 pw(igrid)%dxcoarse(ixCoG^S,1){^NOONED &
                *two*dabs(dsin(pw(igrid)%xcoarse(ixCoG^S,2))) &
                *dsin(half*pw(igrid)%dxcoarse(ixCoG^S,2))}{^IFTHREED*pw(igrid)%dxcoarse(ixCoG^S,3)}
         pw(igrid)%ds(ixGext^S,1)=pw(igrid)%dx(ixGext^S,1)
         {^NOONED pw(igrid)%ds(ixGext^S,2)=xext(ixGext^S,1)*pw(igrid)%dx(ixGext^S,2)}
         {^IFTHREED pw(igrid)%ds(ixGext^S,3)= &
                  xext(ixGext^S,1)*dsin(xext(ixGext^S,2))*pw(igrid)%dx(ixGext^S,3)}
      case ("cylindrical")
         pw(igrid)%dvolume(ixGext^S)=dabs(xext(ixGext^S,1)) &
              *pw(igrid)%dx(ixGext^S,1){^DE&*pw(igrid)%dx(ixGext^S,^DE) }
         pw(igrid)%dvolumecoarse(ixCoG^S)=dabs(pw(igrid)%xcoarse(ixCoG^S,1)) &
              *pw(igrid)%dxcoarse(ixCoG^S,1){^DE&*pw(igrid)%dxcoarse(ixCoG^S,^DE) }
         pw(igrid)%ds(ixGext^S,r_)=pw(igrid)%dx(ixGext^S,r_)
         if(z_>0) pw(igrid)%ds(ixGext^S,z_)=pw(igrid)%dx(ixGext^S,z_)
         if (phi_ > 0) then
           {if (^DE==phi_) pw(igrid)%ds(ixGext^S,^DE)= &
                    xext(ixGext^S,1)*pw(igrid)%dx(ixGext^S,^DE)\}
         end if
      case default
         call mpistop("Sorry, typeaxial unknown")
      end select
endif

if (B0field) then
   ! initialize background non-evolving solution
   call alloc_B0_grid(igrid)
   call set_B0_grid(igrid)
end if

! find the blocks on the boundaries
pw(igrid)%is_physical_boundary=.false.
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
       pw(igrid)%is_physical_boundary(2*^D)=.false.
     else
       pw(igrid)%is_physical_boundary(2*^D)=.true.
     end if
  else if (ign^D < 1) then
     if(phi_ > 0 .and. poleB(1,^D)) then
       ! if at a pole, the boundary is not physical boundary
       pw(igrid)%is_physical_boundary(2*^D-1)=.false.
     else
       pw(igrid)%is_physical_boundary(2*^D-1)=.true.
     end if
  end if
end do
\}
if(any(pw(igrid)%is_physical_boundary)) then
  phyboundblock(igrid)=.true.
else
  phyboundblock(igrid)=.false.
end if

end subroutine alloc_node
!=============================================================================
subroutine dealloc_node(igrid)

use mod_global_parameters

integer, intent(in) :: igrid
!-----------------------------------------------------------------------------
if (igrid==0) then
   call mpistop("trying to delete a non-existing grid in dealloc_node")
end if

deallocate(pw(igrid)%w,pw(igrid)%w1,pw(igrid)%x)
deallocate(pw(igrid)%wcoarse,pw(igrid)%xcoarse)
deallocate(pw(igrid)%wold,pw(igrid)%wio)
! deallocate temperary solution space
select case (time_integrator)
case("threestep","fourstep","jameson","twostep_trapezoidal")
  deallocate(pw(igrid)%w2)
case("rk4","ssprk43")
  deallocate(pw(igrid)%w2)
  deallocate(pw(igrid)%w3)
case("ssprk54")
  deallocate(pw(igrid)%w2)
  deallocate(pw(igrid)%w3)
  deallocate(pw(igrid)%w4)
end select
if(allocated(pw(igrid)%w2)) deallocate(pw(igrid)%w2)
if(allocated(pw(igrid)%w3)) deallocate(pw(igrid)%w3)

if (.not.slab) call putgridgeo(igrid)

if (B0field) call dealloc_B0_grid(igrid)

end subroutine dealloc_node
!=============================================================================
