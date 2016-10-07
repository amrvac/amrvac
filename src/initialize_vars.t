!=============================================================================
!> Initialize (and allocate) simulation and grid variables
!> @todo Explain which ones are not initialized here
subroutine initialize_vars
use mod_forest
use mod_global_parameters

integer :: igrid, level, ipe, ig^D
logical :: ok
!-----------------------------------------------------------------------------

! set time, time counter
if(.not.treset)t=zero
if(.not.itreset)it=0
dt=zero
dtimpl=zero
itmin=0

if (residmin > smalldouble) then
   residual = one
else
   residual = zero
endif

! set all dt to zero
dt_grid(1:ngridshi)=zero

! check resolution
if ({mod(ixGhi^D,2)/=0|.or.}) then
   call mpistop("mesh widths must give even number grid points")
end if
ixM^LL=ixG^LL^LSUBdixB;

if (nbufferx^D>(ixMhi^D-ixMlo^D+1)|.or.) then
   write(unitterm,*) "nbufferx^D bigger than mesh size makes no sense."
   write(unitterm,*) "Decrease nbufferx or increase mesh size"
   call mpistop("")
end if

! initialize dx arrays on finer (>1) levels
do level=2,mxnest
   {dx(^D,level) = dx(^D,level-1) * half\}  ! refine ratio 2
end do

! domain decomposition
! physical extent of a grid block at level 1, per dimension
^D&dg^D(1)=dx(^D,1)*dble(ixGhi^D-2*dixB)\
! number of grid blocks at level 1 in simulation domain, per dimension
^D&ng^D(1)=nint((xprobmax^D-xprobmin^D)/dg^D(1))\
! total number of grid blocks at level 1
nglev1={ng^D(1)*}

do level=2,mxnest
   dg^D(level)=half*dg^D(level-1);
   ng^D(level)=ng^D(level-1)*2;
end do

! check that specified stepsize correctly divides domain
ok=({(abs(dble(ng^D(1))*dg^D(1)-(xprobmax^D-xprobmin^D))<=smalldouble)|.and.})
if (.not.ok) then
   write(unitterm,*)"domain cannot be divided by meshes of given gridsize"
   call mpistop("domain cannot be divided by meshes of given gridsize")
end if


poleB=.false.
if (.not.slab) call set_pole

do igrid=1,ngridshi
   nullify(pwold(igrid)%w,pw(igrid)%w,pw1(igrid)%w, &
           pwCoarse(igrid)%w,pwCoCo(igrid)%w)
   nullify(px(igrid)%x,pxCoarse(igrid)%x)
   nullify(pgeo(igrid)%surfaceC^D,pgeo(igrid)%surface^D, &
           pgeo(igrid)%dvolume,pgeo(igrid)%dx)
   nullify(pgeoCoarse(igrid)%surfaceC^D,pgeoCoarse(igrid)%surface^D, &
           pgeoCoarse(igrid)%dvolume,pgeoCoarse(igrid)%dx)
   if (B0field) then
    nullify(pB0_cell(igrid)%w,pB0_face^D(igrid)%w)
   end if
   if (nstep>2) then
     nullify(pw2(igrid)%w)
   end if
   if (nstep>3) then
     nullify(pw3(igrid)%w)
   end if
   if (nstep>4) then
      nullify(pw4(igrid)%w)
   end if
   if (residmin>smalldouble) then
     nullify(pwres(igrid)%w)
  end if
end do

! on each processor, create for later use a default patch array
allocate(patchfalse(ixG^T))
patchfalse(ixG^T)=.false.

! initialize connectivity data
igridstail=0

! allocate memory for forest data structures
allocate(level_head(mxnest),level_tail(mxnest))
do level=1,mxnest
   nullify(level_head(level)%node,level_tail(level)%node)
end do

allocate(igrid_to_node(ngridshi,0:npe-1))
do ipe=0,npe-1
   do igrid=1,ngridshi
      nullify(igrid_to_node(igrid,ipe)%node)
   end do
end do

allocate(sfc(1:3,ngridshi*npe))

allocate(igrid_to_sfc(ngridshi))

sfc=0
allocate(Morton_start(0:npe-1),Morton_stop(0:npe-1))
allocate(Morton_sub_start(0:npe-1),Morton_sub_stop(0:npe-1))

allocate(nleafs_level(1:nlevelshi))

allocate(coarsen(ngridshi,0:npe-1),refine(ngridshi,0:npe-1))
coarsen=.false.
refine=.false.
if (nbufferx^D/=0|.or.) then
   allocate(buffer(ngridshi,0:npe-1))
   buffer=.false.
end if
allocate(igrid_inuse(ngridshi,0:npe-1))
igrid_inuse=.false.

allocate(tree_root(1:ng^D(1)))
{do ig^DB=1,ng^DB(1)\}
   nullify(tree_root(ig^D)%node)
{end do\}

{#IFDEF STRETCHGRID
logGs(1)=logG
qsts(1)=qst
qsts(0)=qst**2
logGs(0)=2.d0*(qsts(0)-1.d0)/(qsts(0)+1.d0)
if(mxnest>1) then
  do level=2,mxnest
    qsts(level)=dsqrt(qsts(level-1))
    logGs(level)=2.d0*(qsts(level)-1.d0)/(qsts(level)+1.d0) 
  end do
end if
}

! default the physical scaling parameters:
UNIT_LENGTH   = ONE
UNIT_DENSITY  = ONE
UNIT_VELOCITY = ONE

end subroutine initialize_vars
!=============================================================================
