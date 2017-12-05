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

integer :: level, ig^D, ign^D, ixCoG^L, ixCoCoG^L, ix, i^D
integer :: imin, imax, index, igCo^D, ixshift
double precision :: rXmin^D, dx^D
logical, save:: first(1:ndim)=.true.

integer :: itheta
double precision :: aval, cval
double precision, allocatable :: dtheta(:), theta(:), sumdtheta(:), &
                                dthetacoarse(:), sumdthetacoarse(:)
!-----------------------------------------------------------------------------
ixCoGmin^D=1;
!ixCoGmax^D=ixGhi^D/2+nghostcells;
ixCoGmax^D=(ixGhi^D-2*nghostcells)/2+2*nghostcells;

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

if(stretched_grid) then
  {if(stretched_dim(^D))then
    imin=(ig^D-1)*block_nx^D+1
    imax=ig^D*block_nx^D
    rnode(rpxmin^D_,igrid)=xprobmin^D+dxfirst_1mq(level,^D) &
                 *(1.0d0-qstretch(level,^D)**(imin-1))
    rnode(rpxmax^D_,igrid)=xprobmin^D+dxfirst_1mq(level,^D) &
                 *(1.0d0-qstretch(level,^D)**imax)
    ! fix possible out of bound due to precision
    if(rnode(rpxmax^D_,igrid)>xprobmax^D) then
       if(first(^D)) then
          write(*,*) 'Warning: edge beyond domain?', ^D,igrid,imax,rnode(rpxmax^D_,igrid)
          first(^D)=.false.
       endif
       rnode(rpxmax^D_,igrid)=xprobmax^D
    endif
    ixshift=(ig^D-1)*block_nx^D-nghostcells
    do ix=ixGlo^D,ixGhi^D
      index=ixshift+ix
      pw(igrid)%x(ix^D%ixG^T,^D)=xprobmin^D+dxfirst_1mq(level,^D)&
                                *(1.0d0-qstretch(level,^D)**(index-1)) &
                   + 0.5d0*dxfirst(level,^D)*qstretch(level,^D)**(index-1)
    enddo
    igCo^D=(ig^D-1)/2
    ixshift=igCo^D*block_nx^D+(1-mod(ig^D,2))*block_nx^D/2-nghostcells
    do ix=ixCoGmin^D,ixCoGmax^D
      index=ixshift+ix
      pw(igrid)%xcoarse(ix^D%ixCoG^S,^D)=xprobmin^D+dxfirst_1mq(level-1,^D)&
                                         *(1.0d0-qstretch(level-1,^D)**(index-1)) &
                  + 0.5d0*dxfirst(level-1,^D)*qstretch(level-1,^D)**(index-1)
    end do
   endif\}
  {if(stretched_symm_dim(^D))then
    allocate(theta(1:block_nx^D*ng^D(level)))
    allocate(dtheta(1:block_nx^D*ng^D(level)))
    allocate(sumdtheta(1:block_nx^D*ng^D(level)))
    allocate(dthetacoarse(1:block_nx^D*ng^D(level)/2))
    allocate(sumdthetacoarse(1:block_nx^D*ng^D(level)/2))
    dtheta(1)=dxfirst(level,^D)
    theta(1)=dxfirst(level,^D)*0.5d0
    sumdtheta(1)=0.0d0
    do itheta=2,block_nx^D*ng^D(level)
       aval=3.0d0*dsin(theta(itheta-1))*dsin(dtheta(itheta-1)*0.5d0) &
             +dcos(theta(itheta-1))*dcos(dtheta(itheta-1)*0.5d0)
       cval=dsin(theta(itheta-1)+dtheta(itheta-1)*0.5d0)
       dtheta(itheta)=acos(-aval*dsqrt(1.0d0-cval**2)+cval*dsqrt(1.0d0-aval**2))
       theta(itheta)=theta(itheta-1)+dtheta(itheta-1)*0.5d0+dtheta(itheta)*0.5d0
       sumdtheta(itheta)=sumdtheta(itheta-1)+dtheta(itheta-1)
    enddo
    dthetacoarse(1)=dtheta(1)+dtheta(2)
    sumdthetacoarse(1)=0.0d0
    do itheta=2,block_nx^D*ng^D(level)/2
       dthetacoarse(itheta)=dtheta(2*itheta-1)+dtheta(2*itheta)
       sumdthetacoarse(itheta)=sumdthetacoarse(itheta-1)+dthetacoarse(itheta-1)
    enddo
    imin=(ig^D-1)*block_nx^D+1
    imax=ig^D*block_nx^D
    rnode(rpxmin^D_,igrid)=xprobmin^D+(xprobmax^D-xprobmin^D)*sumdtheta(imin)/dpi
    rnode(rpxmax^D_,igrid)=xprobmin^D+(xprobmax^D-xprobmin^D)*(sumdtheta(imax)+dtheta(imax))/dpi
    ! fix possible out of bound due to precision
    if(rnode(rpxmax^D_,igrid)>xprobmax^D) then
       if(first(^D)) then
          write(*,*) 'Warning: edge beyond domain?', ^D,igrid,imax,rnode(rpxmax^D_,igrid)
          first(^D)=.false.
       endif
       rnode(rpxmax^D_,igrid)=xprobmax^D
    endif
    ixshift=(ig^D-1)*block_nx^D-nghostcells
    do ix=ixGlo^D,ixGhi^D
      index=ixshift+ix
      pw(igrid)%x(ix^D%ixG^T,^D)=xprobmin^D+(xprobmax^D-xprobmin^D)*sumdtheta(index)/dpi &
                   + 0.5d0*dtheta(index)
    enddo
    igCo^D=(ig^D-1)/2
    ixshift=igCo^D*block_nx^D+(1-mod(ig^D,2))*block_nx^D/2-nghostcells
    do ix=ixCoGmin^D,ixCoGmax^D
      index=ixshift+ix
      pw(igrid)%xcoarse(ix^D%ixCoG^S,^D)=xprobmin^D+(xprobmax^D-xprobmin^D)*sumdthetacoarse(index)/dpi &
                   + 0.5d0*dthetacoarse(index)
    end do
    deallocate(theta)
    deallocate(dtheta)
    deallocate(sumdtheta)
    deallocate(dthetacoarse)
    deallocate(sumdthetacoarse)
   endif\}
endif

if (.not.slab) call getgridgeo(igrid)

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
