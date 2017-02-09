!=============================================================================
!> Generate and initialize all grids at the coarsest level (level one)
subroutine initlevelone
{#IFDEF EVOLVINGBOUNDARY
use mod_forest
}

use mod_global_parameters

integer :: iigrid, igrid{#IFDEF EVOLVINGBOUNDARY , Morton_no}
!-----------------------------------------------------------------------------
levmin=1
levmax=1

call init_forest_root

call getigrids
call build_connectivity

! fill solution space of all root grids
do iigrid=1,igridstail; igrid=igrids(iigrid);
   saveigrid=igrid
   call alloc_node(igrid)
   ! in case gradient routine used in initial condition, ensure geometry known
   ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
   if (.not.slab) mygeo => pgeo(igrid)
   if (B0field) then
      myB0_cell => pB0_cell(igrid)
      {^D&myB0_face^D => pB0_face^D(igrid)\}
   end if
   typelimiter=limiter(node(plevel_,igrid))
   typegradlimiter=gradient_limiter(node(plevel_,igrid))
   call initial_condition(igrid)
end do
{#IFDEF EVOLVINGBOUNDARY
! mark physical-boundary blocks on space-filling curve
do Morton_no=Morton_start(mype),Morton_stop(mype)
   igrid=sfc_to_igrid(Morton_no)
   if (phyboundblock(igrid)) sfc_phybound(Morton_no)=1
end do
call MPI_ALLREDUCE(MPI_IN_PLACE,sfc_phybound,nleafs,MPI_INTEGER,&
                   MPI_SUM,icomm,ierrmpi)
}

call selectgrids

end subroutine initlevelone
!=============================================================================
subroutine initial_condition(igrid)

! Need only to set the mesh values (can leave ghost cells untouched)
use mod_usr_methods, only: usr_init_one_grid
use mod_global_parameters

integer, intent(in) :: igrid

!----------------------------------------------------------------------------
pw(igrid)%w(ixG^T,1:nw)=zero

saveigrid=igrid
! in case gradient routine used in initial condition, ensure geometry known
^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
if (.not.slab) mygeo => pgeo(igrid)
if (B0field) then
   myB0_cell => pB0_cell(igrid)
   {^D&myB0_face^D => pB0_face^D(igrid)\}
end if
typelimiter=limiter(node(plevel_,igrid))
typegradlimiter=gradient_limiter(node(plevel_,igrid))

if (.not. associated(usr_init_one_grid)) then
   call mpistop("usr_init_one_grid not defined")
else
   call usr_init_one_grid(ixG^LL,ixM^LL,pw(igrid)%w,px(igrid)%x)
end if

end subroutine initial_condition
!=============================================================================
subroutine modify_IC
use mod_usr_methods, only: usr_init_one_grid
use mod_global_parameters

integer :: iigrid, igrid

!-----------------------------------------------------------------------------
do iigrid=1,igridstail; igrid=igrids(iigrid);
   saveigrid=igrid

   if (.not. associated(usr_init_one_grid)) then
      call mpistop("usr_init_one_grid not defined")
   else
      call usr_init_one_grid(ixG^LL,ixM^LL,pw(igrid)%w,px(igrid)%x)
   end if
end do

end subroutine modify_IC
!=============================================================================
