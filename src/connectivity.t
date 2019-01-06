!=============================================================================
subroutine get_level_range
use mod_forest
use mod_global_parameters

integer :: level
!-----------------------------------------------------------------------------

! determine new finest level
do level=refine_max_level,1,-1
   if (associated(level_tail(level)%node)) then
      levmax=level
      exit
   end if
end do

! determine coarsest level
do level=1,levmax
   if (associated(level_tail(level)%node)) then
      levmin=level
      exit
   end if
end do

end subroutine get_level_range
!=============================================================================
subroutine getigrids
use mod_forest
use mod_global_parameters
implicit none

integer :: iigrid, igrid
!-----------------------------------------------------------------------------
iigrid=0
do igrid=1,max_blocks
   if (igrid_inuse(igrid,mype)) then
      iigrid=iigrid+1
      igrids(iigrid)=igrid
   end if
end do

igridstail=iigrid

end subroutine getigrids
!=============================================================================
subroutine build_connectivity
use mod_forest
use mod_global_parameters

integer :: iigrid, igrid, i^D, my_neighbor_type
integer :: iside, idim, ic^D, inc^D, ih^D, icdim
type(tree_node_ptr) :: tree, my_neighbor, child
logical, dimension(^ND) :: pole
logical :: nopole
! Variables to detect special corners for stagger grid
integer :: idir,pi^D, mi^D, ph^D, mh^D, ipe_neighbor
!-----------------------------------------------------------------------------
nrecv_bc_srl=0; nsend_bc_srl=0
nrecv_bc_r=0; nsend_bc_r=0
nrecv_bc_p=0; nsend_bc_p=0
nrecv_fc=0; nsend_fc=0
if(stagger_grid) nrecv_cc=0; nsend_cc=0

do iigrid=1,igridstail; igrid=igrids(iigrid);
   tree%node => igrid_to_node(igrid,mype)%node

   {do i^DB=-1,1\}
      ! skip the grid itself
      if (i^D==0|.and.) then
         neighbor_type(0^D&,igrid)=0
         neighbor(1,0^D&,igrid)=igrid
         neighbor(2,0^D&,igrid)=mype
      else
         call find_neighbor(my_neighbor,my_neighbor_type,tree,i^D,pole)
         nopole=.not.any(pole)

         select case (my_neighbor_type)
         ! adjacent to physical boundary
         case (neighbor_boundary)
            neighbor(1,i^D,igrid)=0
            neighbor(2,i^D,igrid)=-1
         ! fine-coarse transition
         case (neighbor_coarse)
            neighbor(1,i^D,igrid)=my_neighbor%node%igrid
            neighbor(2,i^D,igrid)=my_neighbor%node%ipe
            if (my_neighbor%node%ipe/=mype) then
               ic^D=1+modulo(tree%node%ig^D-1,2);
               if ({(i^D==0.or.i^D==2*ic^D-3)|.and.}) then
                  nrecv_bc_p=nrecv_bc_p+1
                  nsend_bc_r=nsend_bc_r+1
               end if
            end if
         ! same refinement level
         case (neighbor_sibling)
            neighbor(1,i^D,igrid)=my_neighbor%node%igrid
            neighbor(2,i^D,igrid)=my_neighbor%node%ipe
            if (my_neighbor%node%ipe/=mype) then
               nrecv_bc_srl=nrecv_bc_srl+1
               nsend_bc_srl=nsend_bc_srl+1
            end if
         ! coarse-fine transition
         case (neighbor_fine)
            neighbor(1,i^D,igrid)=0
            neighbor(2,i^D,igrid)=-1
            {do ic^DB=1+int((1-i^DB)/2),2-int((1+i^DB)/2)
               inc^DB=2*i^DB+ic^DB
               if (pole(^DB)) then
                  ih^DB=3-ic^DB
               else
                  ih^DB=ic^DB
               end if\}
               child%node => my_neighbor%node%child(ih^D)%node
               neighbor_child(1,inc^D,igrid)=child%node%igrid
               neighbor_child(2,inc^D,igrid)=child%node%ipe
               if (child%node%ipe/=mype) then
                  nrecv_bc_r=nrecv_bc_r+1
                  nsend_bc_p=nsend_bc_p+1
               end if
            {end do\}
         end select

         ! flux fix for conservation only for pure directional shifts
         if ({abs(i^D)+}==1) then
            {if (i^D/=0) then
               idim=^D
               iside=int((i^D+3)/2)
            end if\}
            select case (my_neighbor_type)
            ! only across fine-coarse or coarse-fine boundaries
            case (neighbor_coarse)
               if (my_neighbor%node%ipe/=mype) then
                  if (.not.pole(idim)) nsend_fc(idim)=nsend_fc(idim)+1
               end if
            case (neighbor_fine)
               if (pole(idim)) then
                  icdim=iside
               else
                  icdim=3-iside
               end if
               select case (idim)
               {case (^D)
                  {do ic^D=icdim,icdim^D%do ic^DD=1,2\}
                     child%node => my_neighbor%node%child(ic^DD)%node
                     if (child%node%ipe/=mype) then
                        if (.not.pole(^D)) nrecv_fc(^D)=nrecv_fc(^D)+1
                     end if
                  {end do\} \}
               end select
            end select
         end if

         if (phi_ > 0) then
           neighbor_pole(i^D,igrid)=0
           if (my_neighbor_type>1) then
             do idim=1,^ND
               if (pole(idim)) then
                 neighbor_pole(i^D,igrid)=idim
                 exit ! there can only be one pole between two meshes
               end if
             end do
           end if
         end if
         neighbor_type(i^D,igrid)=my_neighbor_type

      end if
   {end do\}

   if(stagger_grid) then
   !Now all the neighbour information is known.
   !Check if there are special corners that need to be communicated
   !To determine whether to send/receive, we must check three neighbours
    {do i^DB=-1,1\}
       if ({abs(i^D)+}==1) then
         if (neighbor_pole(i^D,igrid)/=0) cycle
        ! Assign value to idim and iside
        {if (i^D/=0) then
            idim=^D
            iside=int((i^D+3)/2)
         end if\}
         ! Fine block surrounded by coarse blocks
         if (neighbor_type(i^D,igrid)==2) then
           do idir=idim+1,ndim
             pi^D=i^D+kr(idir,^D);
             mi^D=i^D-kr(idir,^D);
             ph^D=pi^D-kr(idim,^D)*(2*iside-3);
             mh^D=mi^D-kr(idim,^D)*(2*iside-3);

             if (neighbor_type(pi^D,igrid)==2.and.&
                 neighbor_type(ph^D,igrid)==2.and.&
                 mype/=neighbor(2,pi^D,igrid).and.&
                 neighbor_pole(pi^D,igrid)==0) then
                nsend_cc(idim) = nsend_cc(idim) + 1
             end if

             if (neighbor_type(mi^D,igrid)==2.and.&
                 neighbor_type(mh^D,igrid)==2.and.&
                 mype/=neighbor(2,mi^D,igrid).and.&
                 neighbor_pole(mi^D,igrid)==0) then
                nsend_cc(idim) = nsend_cc(idim) + 1
             end if
           end do
         end if
         ! Coarse block diagonal to fine block(s)
         if (neighbor_type(i^D,igrid)==3) then
           do idir=idim+1,ndim
             pi^D=i^D+kr(idir,^D);
             mi^D=i^D-kr(idir,^D);
             ph^D=pi^D-kr(idim,^D)*(2*iside-3);
             mh^D=mi^D-kr(idim,^D)*(2*iside-3);

             if (neighbor_type(pi^D,igrid)==4.and.&
                 neighbor_type(ph^D,igrid)==3.and.&
                 neighbor_pole(pi^D,igrid)==0) then
              ! Loop on children (several in 3D)
              {do ic^DB=1+int((1-pi^DB)/2),2-int((1+pi^DB)/2)
                  inc^DB=2*pi^DB+ic^DB\}
                 if (mype.ne.neighbor_child(2,inc^D,igrid)) then
                   nrecv_cc(idim) = nrecv_cc(idim) + 1
                 end if
              {end do\}
             end if

             if (neighbor_type(mi^D,igrid)==4.and.&
                 neighbor_type(mh^D,igrid)==3.and.&
                 neighbor_pole(mi^D,igrid)==0) then
              ! Loop on children (several in 3D)
              {do ic^DB=1+int((1-mi^DB)/2),2-int((1+mi^DB)/2)
                  inc^DB=2*mi^DB+ic^DB\}
                 if (mype.ne.neighbor_child(2,inc^D,igrid)) then
                   nrecv_cc(idim) = nrecv_cc(idim) + 1
                 end if
              {end do\}
             end if
           end do
         end if
       end if
    {end do\}
   end if

end do

end subroutine build_connectivity
!=============================================================================
