!=============================================================================
subroutine initlevelone

include 'amrvacdef.f'

integer :: iigrid, igrid
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
   dxlevel(1)=rnode(rpdx1_,igrid);
   if (.not.slab) mygeo => pgeo(igrid)
   if (B0field) then
      myB0_cell => pB0_cell(igrid)
      myB0_face1 => pB0_face1(igrid)
   end if
   typelimiter=typelimiter1(node(plevel_,igrid))
   typegradlimiter=typegradlimiter1(node(plevel_,igrid))
   call initial_condition(igrid)
end do

call selectgrids

end subroutine initlevelone
!=============================================================================
subroutine initial_condition(igrid)

! Need only to set the mesh values (can leave ghost cells untouched)

include 'amrvacdef.f'

integer, intent(in) :: igrid

external initonegrid_usr
!----------------------------------------------------------------------------
pw(igrid)%w(ixGlo1:ixGhi1,1:nw)=zero

saveigrid=igrid
! in case gradient routine used in initial condition, ensure geometry known
dxlevel(1)=rnode(rpdx1_,igrid);
if (.not.slab) mygeo => pgeo(igrid)
if (B0field) then
   myB0_cell => pB0_cell(igrid)
   myB0_face1 => pB0_face1(igrid)
end if
typelimiter=typelimiter1(node(plevel_,igrid))
typegradlimiter=typegradlimiter1(node(plevel_,igrid))

call initonegrid_usr(ixGlo1,ixGhi1,ixMlo1,ixMhi1,pw(igrid)%w,px(igrid)%x)

end subroutine initial_condition
!=============================================================================
subroutine modify_IC

include 'amrvacdef.f'

integer :: iigrid, igrid

external initonegrid_usr
!-----------------------------------------------------------------------------
do iigrid=1,igridstail; igrid=igrids(iigrid);
   saveigrid=igrid
   call initonegrid_usr(ixGlo1,ixGhi1,ixMlo1,ixMhi1,pw(igrid)%w,px(igrid)%x)
end do

end subroutine modify_IC
!=============================================================================
