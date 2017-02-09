!=============================================================================
subroutine amr_coarsen_refine
use mod_forest
use mod_global_parameters
use mod_ghostcells_update

integer :: iigrid, igrid, ipe, igridCo, ipeCo, level, ic^D
integer, dimension(2^D&) :: igridFi, ipeFi
type(tree_node_ptr) :: tree, sibling
logical             :: active
integer, external :: getnode
!-----------------------------------------------------------------------------
if(addmpibarrier) call MPI_BARRIER(icomm,ierrmpi)
call proper_nesting

! to save memory: first coarsen then refine
irecv=0
isend=0
allocate(recvstatus(MPI_STATUS_SIZE,max_blocks),recvrequest(max_blocks), &
         sendstatus(MPI_STATUS_SIZE,max_blocks),sendrequest(max_blocks))
recvrequest=MPI_REQUEST_NULL
sendrequest=MPI_REQUEST_NULL

do ipe=0,npe-1
   do igrid=1,max_blocks
      if (coarsen(igrid,ipe)) then
         if (.not.associated(igrid_to_node(igrid,ipe)%node)) cycle

         tree%node => igrid_to_node(igrid,ipe)%node%parent%node
         {do ic^DB=1,2\}
            sibling%node => tree%node%child(ic^D)%node
            ipeFi(ic^D)=sibling%node%ipe
            igridFi(ic^D)=sibling%node%igrid
         {end do\}

         ipeCo=ipeFi(1^D&)
         igridCo=getnode(ipeCo)

         call coarsen_tree_leaf(igridCo,ipeCo,igridFi,ipeFi,active)

         call coarsen_grid_siblings(igridCo,ipeCo,igridFi,ipeFi,active)

         ! local coarsening done
         {do ic^DB=1,2\}
            if (ipeFi(ic^D)==ipeCo) then
               call putnode(igridFi(ic^D),ipeFi(ic^D))
               coarsen(igridFi(ic^D),ipeFi(ic^D))=.false.
            end if
         {end do\}
      end if
   end do
end do

if (irecv>0) call MPI_WAITALL(irecv,recvrequest,recvstatus,ierrmpi)
if (isend>0) call MPI_WAITALL(isend,sendrequest,sendstatus,ierrmpi)

! non-local coarsening done
do ipe=0,npe-1
   do igrid=1,max_blocks
      if (coarsen(igrid,ipe)) then
         if (ipe==mype) call dealloc_node(igrid)
         call putnode(igrid,ipe)
         coarsen(igrid,ipe)=.false.
      end if
   end do
end do

do ipe=0,npe-1
   do igrid=1,max_blocks
      if (refine(igrid,ipe)) then

         {do ic^DB=1,2\}
            igridFi(ic^D)=getnode(ipe)
            ipeFi(ic^D)=ipe
         {end do\}

         call refine_tree_leaf(igridFi,ipeFi,igrid,ipe,active)

         if (ipe==mype) call refine_grid(igridFi,ipeFi,igrid,ipe,active)

         ! refinement done
         call putnode(igrid,ipe)
         refine(igrid,ipe)=.false.
      end if
   end do
end do

call get_level_range

call load_balance
deallocate(recvstatus,recvrequest,sendstatus,sendrequest)

! Rebuild tree connectivity
call getigrids
call build_connectivity

! Update the list of active grids
call selectgrids
!  grid structure now complete again.

! since we only filled mesh values, and advance assumes filled
! ghost cells, do boundary filling for the new levels
if (time_advance) then
   call getbc(global_time+dt,0.d0,pw,0,nwflux+nwaux)
else
   call getbc(global_time,0.d0,pw,0,nwflux+nwaux)
end if



if(addmpibarrier) call MPI_BARRIER(icomm,ierrmpi)
end subroutine amr_coarsen_refine
!=============================================================================
subroutine proper_nesting
use mod_forest
use mod_global_parameters

logical, dimension(:,:), allocatable :: refine2
integer :: iigrid, igrid, level, ic^D, inp^D, i^D, my_neighbor_type,ipe
logical :: coarsening, pole(ndim), sendbuf(max_blocks)
type(tree_node_ptr) :: tree, p_neighbor, my_parent, sibling, my_neighbor, &
                       neighborchild
!-----------------------------------------------------------------------------

! For all grids on all processors, do a check on refinement w_for_refine. Make
! sure that neighbors will not differ more than one level of refinement.

if (nbufferx^D/=0|.or.) then
   allocate(refine2(max_blocks,npe))
   call MPI_ALLREDUCE(refine,refine2,max_blocks*npe,MPI_LOGICAL,MPI_LOR, &
                      icomm,ierrmpi)
   refine=refine2
else
   sendbuf(:)=refine(:,mype)
   call MPI_ALLGATHER(sendbuf,max_blocks,MPI_LOGICAL,refine,max_blocks, &
                      MPI_LOGICAL,icomm,ierrmpi)
end if

do level=min(levmax,refine_max_level-1),levmin+1,-1
   tree%node => level_head(level)%node
   do
      if (.not.associated(tree%node)) exit

      if (refine(tree%node%igrid,tree%node%ipe)) then
         ic^D=1+modulo(tree%node%ig^D-1,2);
         {do inp^DB=ic^DB-2,ic^DB-1\}
            if (inp^D==0|.and.) cycle
            p_neighbor%node => tree%node%parent%node
            {if (inp^D/=0) then
               p_neighbor%node => p_neighbor%node%neighbor(ic^D,^D)%node
               if (.not.associated(p_neighbor%node)) cycle
            end if\}
            if (p_neighbor%node%leaf) then
               refine(p_neighbor%node%igrid,p_neighbor%node%ipe)=.true.
            end if
         {end do\}
      end if

      tree%node => tree%node%next%node
   end do
end do

! On each processor locally, check if grids set for coarsening are already
! set for refinement.

do iigrid=1,igridstail; igrid=igrids(iigrid);
   if (refine(igrid,mype).and.coarsen(igrid,mype)) coarsen(igrid,mype)=.false.
end do

! For all grids on all processors, do a check on coarsen flags.
sendbuf(:)=coarsen(:,mype)
call MPI_ALLGATHER(sendbuf,max_blocks,MPI_LOGICAL,coarsen,max_blocks, &
                   MPI_LOGICAL,icomm,ierrmpi)

do level=levmax,max(2,levmin),-1
   tree%node => level_head(level)%node
   do
      if (.not.associated(tree%node)) exit

      if (coarsen(tree%node%igrid,tree%node%ipe)) then
         coarsening=.true.
         my_parent%node => tree%node%parent%node

         ! are all siblings flagged for coarsen ?
check1:  {do ic^DB=1,2\}
            sibling%node => my_parent%node%child(ic^D)%node
            if (sibling%node%leaf) then
               if (coarsen(sibling%node%igrid,sibling%node%ipe)) cycle
            end if
            call unflag_coarsen_siblings
            exit check1
         {end do\} check1

         ! Make sure that neighbors will not differ more than one level of
         ! refinement, otherwise unflag all siblings
         if (coarsening) then
check2:     {do ic^DB=1,2\}
               sibling%node => my_parent%node%child(ic^D)%node
               {do i^DB=ic^DB-2,ic^DB-1\}
                  if (i^D==0|.and.) cycle
                  call find_neighbor(my_neighbor,my_neighbor_type, &
                                     sibling,i^D,pole)
                  select case (my_neighbor_type)
                  case (3)
                     if (refine(my_neighbor%node%igrid, &
                                my_neighbor%node%ipe)) then
                        call unflag_coarsen_siblings
                        exit check2
                     else
                        cycle
                     end if
                  case (4)
                     neighborchild%node=>my_neighbor%node%child(1^D&)%node
                     if (neighborchild%node%leaf) then
                        if (coarsen(neighborchild%node%igrid, &
                                    neighborchild%node%ipe)) then
                           cycle
                        end if
                     end if
                     call unflag_coarsen_siblings
                     exit check2
                  end select
               {end do\}
            {end do\} check2
         end if

      end if

      tree%node => tree%node%next%node
   end do
end do

contains
!=============================================================================
! internal procedures
!=============================================================================
subroutine unflag_coarsen_siblings

integer :: ic^D
type(tree_node_ptr) :: sibling
!-----------------------------------------------------------------------------

{do ic^DB=1,2\}
   sibling%node => my_parent%node%child(ic^D)%node
   if (sibling%node%leaf) then
      coarsen(sibling%node%igrid,sibling%node%ipe)=.false.
   end if
{end do\}
coarsening=.false.

end subroutine unflag_coarsen_siblings
!=============================================================================
! end of internal procedures
!=============================================================================
end subroutine proper_nesting
!=============================================================================
