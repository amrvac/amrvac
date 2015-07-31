!=============================================================================
subroutine amr_coarsen_refine
use mod_forest
include 'amrvacdef.f'

integer :: iigrid, igrid, ipe, igridCo, ipeCo, level, ic1
integer, dimension(2) :: igridFi, ipeFi
type(tree_node_ptr) :: tree, sibling
logical             :: active
integer, external :: getnode
!-----------------------------------------------------------------------------
if(addmpibarrier) call MPI_BARRIER(icomm,ierrmpi)
call proper_nesting

! to save memory: first coarsen then refine
irecv=0
isend=0
allocate(recvstatus(MPI_STATUS_SIZE,ngridshi),recvrequest(ngridshi),&
    sendstatus(MPI_STATUS_SIZE,ngridshi),sendrequest(ngridshi))
recvrequest=MPI_REQUEST_NULL
sendrequest=MPI_REQUEST_NULL

do ipe=0,npe-1
   do igrid=1,ngridshi
      if (coarsen(igrid,ipe)) then
         if (.not.associated(igrid_to_node(igrid,ipe)%node)) cycle

         tree%node => igrid_to_node(igrid,ipe)%node%parent%node
         do ic1=1,2
            sibling%node => tree%node%child(ic1)%node
            ipeFi(ic1)=sibling%node%ipe
            igridFi(ic1)=sibling%node%igrid
         end do

         ipeCo=ipeFi(1)
         igridCo=getnode(ipeCo)

         call coarsen_tree_leaf(igridCo,ipeCo,igridFi,ipeFi,active)

         call coarsen_grid_siblings(igridCo,ipeCo,igridFi,ipeFi,active)

         ! local coarsening done
         do ic1=1,2
            if (ipeFi(ic1)==ipeCo) then
               call putnode(igridFi(ic1),ipeFi(ic1))
               coarsen(igridFi(ic1),ipeFi(ic1))=.false.
            end if
         end do
      end if
   end do
end do

if (irecv>0) call MPI_WAITALL(irecv,recvrequest,recvstatus,ierrmpi)
if (isend>0) call MPI_WAITALL(isend,sendrequest,sendstatus,ierrmpi)

! non-local coarsening done
do ipe=0,npe-1
   do igrid=1,ngridshi
      if (coarsen(igrid,ipe)) then
         if (ipe==mype) call dealloc_node(igrid)
         call putnode(igrid,ipe)
         coarsen(igrid,ipe)=.false.
      end if
   end do
end do

do ipe=0,npe-1
   do igrid=1,ngridshi
      if (refine(igrid,ipe)) then

         do ic1=1,2
            igridFi(ic1)=getnode(ipe)
            ipeFi(ic1)=ipe
         end do

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
   call getbc(t+dt,ixGlo1,ixGhi1,pw,pwCoarse,pgeo,pgeoCoarse,.false.)
else
   call getbc(t,ixGlo1,ixGhi1,pw,pwCoarse,pgeo,pgeoCoarse,.false.)
end if



if(addmpibarrier) call MPI_BARRIER(icomm,ierrmpi)
end subroutine amr_coarsen_refine
!=============================================================================
subroutine proper_nesting
use mod_forest
include 'amrvacdef.f'

!! following alternative
!!logical, dimension(1:ngridshi,0:npe-1):: coarsen2

integer :: iigrid, igrid, level, ic1, inp1, i1, my_neighbor_type,ipe
logical :: coarsening, pole(ndim)
type(tree_node_ptr) :: tree, p_neighbor, my_parent, sibling, my_neighbor,&
    neighborchild
!-----------------------------------------------------------------------------

! For all grids on all processors, do a check on refinement flags. Make
! sure that neighbors will not differ more than one level of refinement.

if (nbufferx1/=0) then
   call MPI_ALLREDUCE(MPI_IN_PLACE,refine,ngridshi*npe,MPI_LOGICAL,MPI_LOR,&
       icomm,ierrmpi)
else
   call MPI_ALLGATHER(MPI_IN_PLACE,ngridshi,MPI_LOGICAL,refine,ngridshi,&
       MPI_LOGICAL,icomm,ierrmpi)
end if

do level=min(levmax,mxnest-1),levmin+1,-1
   tree%node => level_head(level)%node
   do
      if (.not.associated(tree%node)) exit

      if (refine(tree%node%igrid,tree%node%ipe)) then
         ic1=1+modulo(tree%node%ig1-1,2);
         do inp1=ic1-2,ic1-1
            if (inp1==0) cycle
            p_neighbor%node => tree%node%parent%node
            if (inp1/=0) then
               p_neighbor%node => p_neighbor%node%neighbor(ic1,1)%node
               if (.not.associated(p_neighbor%node)) cycle
            end if
            if (p_neighbor%node%leaf) then
               refine(p_neighbor%node%igrid,p_neighbor%node%ipe)=.true.
            end if
         end do
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

!!coarsen2(1:ngridshi,:)=coarsen(1:ngridshi,:)

!!do ipe=0,npe-1
!! if(ipe/=mype) coarsen2(1:ngridshi,ipe)=.false.
!!enddo

call MPI_ALLGATHER(MPI_IN_PLACE,ngridshi,MPI_LOGICAL,coarsen,ngridshi,&
    MPI_LOGICAL,icomm,ierrmpi)
!!call MPI_ALLREDUCE(MPI_IN_PLACE,coarsen2,ngridshi*npe,MPI_LOGICAL,MPI_LOR, &
!!                      icomm,ierrmpi)

!!coarsen(1:ngridshi,:)=coarsen2(1:ngridshi,:)

do level=levmax,max(2,levmin),-1
   tree%node => level_head(level)%node
   do
      if (.not.associated(tree%node)) exit

      if (coarsen(tree%node%igrid,tree%node%ipe)) then
         coarsening=.true.
         my_parent%node => tree%node%parent%node

         ! are all siblings flagged for coarsen ?
check1:  do ic1=1,2
            sibling%node => my_parent%node%child(ic1)%node
            if (sibling%node%leaf) then
               if (coarsen(sibling%node%igrid,sibling%node%ipe)) cycle
            end if
            call unflag_coarsen_siblings
            exit check1
         end do check1

         ! Make sure that neighbors will not differ more than one level of
         ! refinement, otherwise unflag all siblings
         if (coarsening) then
check2:     do ic1=1,2
               sibling%node => my_parent%node%child(ic1)%node
               do i1=ic1-2,ic1-1
                  if (i1==0) cycle
                  call find_neighbor(my_neighbor,my_neighbor_type, sibling,i1,&
                     pole)
                  select case (my_neighbor_type)
                  case (3)
                     if (refine(my_neighbor%node%igrid, my_neighbor%node%ipe))&
                         then
                        call unflag_coarsen_siblings
                        exit check2
                     else
                        cycle
                     end if
                  case (4)
                     neighborchild%node=>my_neighbor%node%child(1)%node
                     if (neighborchild%node%leaf) then
                        if (coarsen(neighborchild%node%igrid,&
                            neighborchild%node%ipe)) then
                           cycle
                        end if
                     end if
                     call unflag_coarsen_siblings
                     exit check2
                  end select
               end do
            end do check2
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

integer :: ic1
type(tree_node_ptr) :: sibling
!-----------------------------------------------------------------------------

do ic1=1,2
   sibling%node => my_parent%node%child(ic1)%node
   if (sibling%node%leaf) then
      coarsen(sibling%node%igrid,sibling%node%ipe)=.false.
   end if
end do
coarsening=.false.

end subroutine unflag_coarsen_siblings
!=============================================================================
! end of internal procedures
!=============================================================================
end subroutine proper_nesting
!=============================================================================
