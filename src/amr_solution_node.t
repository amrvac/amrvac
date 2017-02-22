!=============================================================================
integer function getnode(ipe)
use mod_forest, only: igrid_inuse
use mod_global_parameters

! getnode = get first available igrid on processor ipe

integer, intent(in) :: ipe

integer :: igrid, igrid_available
!----------------------------------------------------------------------------
igrid_available=0

do igrid=1,max_blocks
   if (igrid_inuse(igrid,ipe)) cycle

   igrid_available=igrid
   exit
end do

if (igrid_available==0) then
   write(unitterm,*) " out of nodal space - allowed ",max_blocks," grids"
   call mpistop("")
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

integer :: level, ig^D, ixCoG^L, ixCoCoG^L, ix
double precision :: rXmin^D, dx^D
!-----------------------------------------------------------------------------
ixCoGmin^D=1;
ixCoGmax^D=ixGhi^D/2+nghostcells;

! initialize solution space
allocate(pwold(igrid)%w(ixG^T,1:nw), &
         pw(igrid)%w(ixG^T,1:nw), &
         pwCoarse(igrid)%w(ixCoG^S,1:nw))
pwold(igrid)%w(ixG^T,1:nw)=0.d0
pw(igrid)%w(ixG^T,1:nw)=0.d0
pwCoarse(igrid)%w(ixCoG^S,1:nw)=0.d0

! initialize background non-evolving solution
if(B0field) then
  allocate(pw(igrid)%w0(ixG^T,1:ndir), &
           pw(igrid)%w0face(ixG^T,1:ndir,1:ndim))
  pw(igrid)%w0=0.d0
  pw(igrid)%w0face=0.d0
end if

! initialize coordinates
allocate(pw(igrid)%x(ixG^T,1:ndim), &
         pwCoarse(igrid)%x(ixCoG^S,1:ndim))


! set level information
level=igrid_to_node(igrid,mype)%node%level
ig^D=igrid_to_node(igrid,mype)%node%ig^D;

node(plevel_,igrid)=level
^D&node(pig^D_,igrid)=ig^D\

! set dx information
^D&rnode(rpdx^D_,igrid)=dx(^D,level)\

! determine the minimal and maximal corners
^D&rnode(rpxmin^D_,igrid)=xprobmin^D+dble(ig^D-1)*dg^D(level)\
^D&rnode(rpxmax^D_,igrid)=xprobmax^D-dble(ng^D(level)-ig^D)*dg^D(level)\


allocate(px(igrid)%x(ixG^T,1:ndim),pxCoarse(igrid)%x(ixCoG^S,1:ndim))
^D&dx^D=rnode(rpdx^D_,igrid)\
^D&rXmin^D=rnode(rpxmin^D_,igrid)-nghostcells*dx^D\
{do ix=ixGlo^D,ixGhi^D
   px(igrid)%x(ix^D%ixG^T,^D)=rXmin^D+(dble(ix)-half)*dx^D
end do\}
^D&dx^D=2.0d0*rnode(rpdx^D_,igrid)\
^D&rXmin^D=rnode(rpxmin^D_,igrid)-nghostcells*dx^D\
{do ix=ixCoGmin^D,ixCoGmax^D
   pxCoarse(igrid)%x(ix^D%ixCoG^S,^D)=rXmin^D+(dble(ix)-half)*dx^D
end do\}

{#IFDEF STRETCHGRID
logG=logGs(level)
qst=qsts(level)
rnode(rpxmin1_,igrid)=xprobmin1*qst**((ixMhi1-ixMlo1+1)*(ig1-1))
rnode(rpxmax1_,igrid)=xprobmin1*qst**((ixMhi1-ixMlo1+1)*ig1)
! fix possible out of bound due to precision
if(rnode(rpxmax1_,igrid)>xprobmax1) rnode(rpxmax1_,igrid)=xprobmax1
rXmin1=rnode(rpxmin1_,igrid)*qst**(-nghostcells)
do ix=ixGlo1,ixGhi1
  px(igrid)%x(ix^%1ixG^T,1)=rXmin1/(one-half*logG)*qst**(ix-1)
end do
rnode(rpdx1_,igrid)=minval(px(igrid)%x(ixG^T,1)*logG)
logG=logGs(level-1)
qst=qsts(level-1)
rXmin1=rnode(rpxmin1_,igrid)*qst**(-nghostcells)
do ix=ixCoGmin1,ixCoGmax1
  pxCoarse(igrid)%x(ix^%1ixCoG^S,1)=rXmin1/(one-half*logG)*qst**(ix-1)
end do
}

! find the blocks on the boundaries
if ({rnode(rpxmin^D_,igrid)-dx(^D,level)<xprobmin^D|.or.}.or.&
    {rnode(rpxmax^D_,igrid)+dx(^D,level)>xprobmax^D|.or.}) then
   phyboundblock(igrid)=.true.
else
   phyboundblock(igrid)=.false.
end if

if (.not.slab) call getgridgeo(igrid)

if (B0field) then
   call alloc_B0_grid(igrid)
   call set_B0_grid(igrid)
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

deallocate(pw(igrid)%w,pwold(igrid)%w,pwCoarse(igrid)%w)
deallocate(px(igrid)%x,pxCoarse(igrid)%x)

if (.not.slab) call putgridgeo(igrid)

if (B0field) call dealloc_B0_grid(igrid)

end subroutine dealloc_node
!=============================================================================
