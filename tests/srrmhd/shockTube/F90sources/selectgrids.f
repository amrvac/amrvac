!=============================================================================
subroutine selectgrids

use mod_forest
include 'amrvacdef.f'
integer :: iigrid, igrid, jgrid, kgrid, isave, my_isafety
integer, allocatable,  dimension(:,:)  :: isafety

! Set the number of safety-blocks (additional blocks after 
! flag_grid_usr): 
integer, parameter :: nsafety = 1
integer             :: ixOmin1,ixOmax1, ipe
type(tree_node_ptr) :: tree
integer             :: igrid_active, userflag
!-----------------------------------------------------------------------------
if (.not. allocated(isafety)) allocate(isafety(ngridshi,0:npe-1))

! reset all grids to active:
neighbor_active = .true.

      jgrid=0
      kgrid=0
      isafety = -1

!     Check the user flag:
      do iigrid=1,igridstail; igrid=igrids(iigrid);
         userflag = igrid_active(igrid)
         if (userflag <= 0) then
            jgrid=jgrid+1
            igrids_active(jgrid)=igrid
         else
            kgrid=kgrid+1
            igrids_passive(kgrid)=igrid
         end if
         isafety(igrid,mype) = userflag
      end do

      igridstail_active  = jgrid
      igridstail_passive = kgrid

!     Check if user wants to deactivate grids at all and return if not:
      if (userflag == -1) return

!     Got the passive grids. 
!     Now, we re-activate a safety belt of radius nsafety blocks.

!     First communicate the current isafety buffer:
      call MPI_ALLGATHER(MPI_IN_PLACE,ngridshi,MPI_INTEGER,isafety,ngridshi,&
         MPI_INTEGER,icomm,ierrmpi)

!     Now check the distance of neighbors to the active zone:
      do isave = 1, nsafety
         do iigrid=1,igridstail_passive; igrid=igrids_passive(iigrid);
            if ( isafety(igrid,mype) /= isave) cycle
!     Get the minimum neighbor isafety:
            if(min_isafety_neighbor(igrid) >= isafety(igrid,mype)) then
!     Increment the buffer:
               isafety(igrid,mype)=isafety(igrid,mype)+1
            end if
         end do
         
!     Communicate the incremented buffers:
         call MPI_ALLGATHER(MPI_IN_PLACE,ngridshi,MPI_INTEGER,isafety,&
            ngridshi,MPI_INTEGER,icomm,ierrmpi)
      end do

!     Update the active and passive arrays:
      jgrid=0
      kgrid=0
      do iigrid=1,igridstail; igrid=igrids(iigrid);
         if (isafety(igrid,mype) <= nsafety) then
            jgrid=jgrid+1
            igrids_active(jgrid)=igrid
         else 
            kgrid=kgrid+1
            igrids_passive(kgrid)=igrid
         end if
!     Create the neighbor flags:
         call set_neighbor_state(igrid)
      end do
      igridstail_active  = jgrid
      igridstail_passive = kgrid

!     Update the tree:
      nleafs_active = 0
      do ipe=0,npe-1
         do igrid=1,ngridshi
            if (isafety(igrid,ipe) == -1) cycle
            if (.not.associated(igrid_to_node(igrid,ipe)%node)) cycle
            tree%node => igrid_to_node(igrid,ipe)%node
            if (isafety(igrid,ipe) > nsafety) then
               tree%node%active=.false.
            else
               tree%node%active=.true.
              nleafs_active = nleafs_active + 1
            end if
         end do
      end do

!     Just for output and testing: 
!      ixO^L=ixG^LL^LSUBdixB;      
!      do iigrid=1,igridstail; igrid=igrids(iigrid);        
!         pw(igrid)%w(ixO^S,flg_) = dble(isafety(igrid,mype))
!         pw(igrid)%w(ixO^S,cpu_) = dble(mype)
!      end do

      contains
!=============================================================================
subroutine set_neighbor_state(igrid)

integer, intent(in)  :: igrid
integer :: my_neighbor_type, i1, isafety_neighbor
!-----------------------------------------------------------------------------


   do i1=-1,1
      if (i1==0) then 
         if (isafety(igrid,mype) > nsafety) neighbor_active(i1,igrid) &
            = .false.
      end if
      my_neighbor_type=neighbor_type(i1,igrid)

      select case (my_neighbor_type)
      case (1) ! boundary
         isafety_neighbor = nsafety+1
      case (2) ! fine-coarse
         isafety_neighbor = isafety_fc(i1,igrid)
      case (3) ! same level
         isafety_neighbor = isafety_srl(i1,igrid)
      case (4) ! coarse-fine
         isafety_neighbor = isafety_cf_max(i1,igrid)
      end select

      if (isafety_neighbor > nsafety) neighbor_active(i1,igrid) = .false.

   end do

end subroutine set_neighbor_state
!=============================================================================
integer function min_isafety_neighbor(igrid)

integer, intent(in) :: igrid
integer :: my_neighbor_type, i1
!-----------------------------------------------------------------------------

min_isafety_neighbor = biginteger

   do i1=-1,1
      if (i1==0) cycle
      my_neighbor_type=neighbor_type(i1,igrid)

      select case (my_neighbor_type)
      case (2) ! fine-coarse
         min_isafety_neighbor = min(isafety_fc(i1,igrid),min_isafety_neighbor)
      case (3) ! same level
         min_isafety_neighbor = min(isafety_srl(i1,igrid),&
            min_isafety_neighbor)
      case (4) ! coarse-fine
         min_isafety_neighbor = min(isafety_cf_min(i1,igrid),&
            min_isafety_neighbor)
      end select

   end do

end function min_isafety_neighbor
!=============================================================================
integer function isafety_fc(i1,igrid)

integer, intent(in) :: i1, igrid
integer            :: ineighbor, ipe_neighbor
!-----------------------------------------------------------------------------
ineighbor=neighbor(1,i1,igrid)
ipe_neighbor=neighbor(2,i1,igrid)

      isafety_fc = isafety(ineighbor,ipe_neighbor)

end function isafety_fc
!=============================================================================
integer function isafety_srl(i1,igrid)

integer, intent(in) :: i1, igrid
integer            :: ineighbor, ipe_neighbor
!-----------------------------------------------------------------------------
ineighbor=neighbor(1,i1,igrid)
ipe_neighbor=neighbor(2,i1,igrid)

      isafety_srl = isafety(ineighbor,ipe_neighbor)

end function isafety_srl
!=============================================================================
integer function isafety_cf_min(i1,igrid)

integer, intent(in) :: i1, igrid
integer            :: ic1, inc1
integer            :: ineighbor, ipe_neighbor
!-----------------------------------------------------------------------------

isafety_cf_min = biginteger

      do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
      inc1=2*i1+ic1
      
        ineighbor    = neighbor_child(1,inc1,igrid)
        ipe_neighbor = neighbor_child(2,inc1,igrid)

        isafety_cf_min = min(isafety_cf_min,isafety(ineighbor,ipe_neighbor))

      end do

end function isafety_cf_min
!=============================================================================
integer function isafety_cf_max(i1,igrid)

integer, intent(in) :: i1, igrid
integer            :: ic1, inc1
integer            :: ineighbor, ipe_neighbor
!-----------------------------------------------------------------------------

isafety_cf_max = - biginteger

      do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
      inc1=2*i1+ic1
      
        ineighbor    = neighbor_child(1,inc1,igrid)
        ipe_neighbor = neighbor_child(2,inc1,igrid)

        isafety_cf_max = max(isafety_cf_max,isafety(ineighbor,ipe_neighbor))

      end do

end function isafety_cf_max
!=============================================================================
end subroutine selectgrids
!=============================================================================
      integer function igrid_active(igrid)
      include 'amrvacdef.f'
      
      integer, intent(in) :: igrid
      integer             :: ixOmin1,ixOmax1, flag
!-----------------------------------------------------------------------------
      ixOmin1=ixGlo1+dixB;ixOmax1=ixGhi1-dixB;

      igrid_active = -1      
      
      call flag_grid_usr(t,ixGlo1,ixGhi1,ixOmin1,ixOmax1,pw(igrid)%w,&
         px(igrid)%x,igrid_active)
      
      end function igrid_active
!=============================================================================

