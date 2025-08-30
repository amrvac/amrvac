module mod_boundary_conditions

  implicit none
  private

  public :: bc_phys
  public :: getintbc

contains

  !> fill ghost cells at a physical boundary
  subroutine bc_phys(iside,idims,time,qdt,s,ixGmin1,ixGmin2,ixGmin3,ixGmax1,&
       ixGmax2,ixGmax3,ixBmin1,ixBmin2,ixBmin3,ixBmax1,ixBmax2,ixBmax3)
    !$acc routine vector
#:if defined('SPECIALBOUNDARY')    
    use mod_usr, only: specialbound_usr
#:endif
    use mod_global_parameters

    integer, intent(in) :: iside, idims, ixGmin1,ixGmin2,ixGmin3,ixGmax1,&
       ixGmax2,ixGmax3,ixBmin1,ixBmin2,ixBmin3,ixBmax1,ixBmax2,ixBmax3
    double precision, intent(in) :: time,qdt
    type(state), intent(inout) :: s

    integer :: idir, is
    integer :: ixOsmin1,ixOsmin2,ixOsmin3,ixOsmax1,ixOsmax2,ixOsmax3,hxOmin1,&
       hxOmin2,hxOmin3,hxOmax1,hxOmax2,hxOmax3,jxOmin1,jxOmin2,jxOmin3,jxOmax1,&
       jxOmax2,jxOmax3
    integer :: iw, iB, ix1,ix2,ix3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
       ixOmax3, ixMmin1,ixMmin2,ixMmin3,ixMmax1,ixMmax2,ixMmax3, nghostcellsi,&
       iib1,iib2,iib3
    logical  :: isphysbound

    associate(x=>s%x,w=>s%w,ws=>s%ws)
    select case (idims)
    case (1)
       if (iside==2) then
          ! maximal boundary
          iB=2*1
          ixOmin1=ixBmax1+1-nghostcells;ixOmin2=ixBmin2;ixOmin3=ixBmin3;
          ixOmax1=ixBmax1;ixOmax2=ixBmax2;ixOmax3=ixBmax3;
          ! cont/symm/asymm types
          do iw=1,nwflux+nwaux
             select case (typeboundary(iw,iB))
             case (bc_special)
                ! skip it here, do AFTER all normal type boundaries are set
             case (bc_cont)
                do ix1=ixOmin1,ixOmax1
                   w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw) = w(ixOmin1-1,&
                      ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw)
                end do
             case (bc_symm)
                w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                   iw) = w(ixOmin1-1:ixOmin1-nghostcells:-1,ixOmin2:ixOmax2,&
                   ixOmin3:ixOmax3,iw)
             case (bc_asymm)
                w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                   iw) =-w(ixOmin1-1:ixOmin1-nghostcells:-1,ixOmin2:ixOmax2,&
                   ixOmin3:ixOmax3,iw)
             case (bc_periodic)
                ! skip it here, periodic bc info should come from neighbors
             case(bc_noinflow)
                if (iw==1+1)then
                  do ix1=ixOmin1,ixOmax1
                      w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                         iw) = max(w(ixOmin1-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                         iw),zero)
                  end do
                else
                  do ix1=ixOmin1,ixOmax1
                      w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw) = w(ixOmin1-1,&
                         ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw)
                  end do
                end if
             case (bc_aperiodic)
 !this just multiplies the variables with (-), they have been set from neighbors just like periodic.
                w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                   iw) = - w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                   iw)
             case (bc_data)
                ! skip it here, do AFTER all normal type boundaries are set
           case (bc_icarus)
              ! skip it here, do AFTER all normal type boundaries are set
             case (bc_character)
                ! skip it here, do AFTER all normal type boundaries are set
             case default
                write (unitterm,*) "Undefined boundarytype found in bc_phys",&
                    "for variable iw=",iw," and side iB=",iB
             end select
          end do
          
       else
          ! minimal boundary
          iB=2*1-1
          ixOmin1=ixBmin1;ixOmin2=ixBmin2;ixOmin3=ixBmin3;
          ixOmax1=ixBmin1-1+nghostcells;ixOmax2=ixBmax2;ixOmax3=ixBmax3;
          ! cont/symm/asymm types
          do iw=1,nwflux+nwaux
             select case (typeboundary(iw,iB))
             case (bc_special)
                ! skip it here, do AFTER all normal type boundaries are set
             case (bc_cont)
                do ix1=ixOmin1,ixOmax1
                   w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw) = w(ixOmax1+1,&
                      ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw)
                end do
             case (bc_symm)
                w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                   iw) = w(ixOmax1+nghostcells:ixOmax1+1:-1,ixOmin2:ixOmax2,&
                   ixOmin3:ixOmax3,iw)
             case (bc_asymm)
                w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                   iw) =-w(ixOmax1+nghostcells:ixOmax1+1:-1,ixOmin2:ixOmax2,&
                   ixOmin3:ixOmax3,iw)
             case (bc_periodic)
                ! skip it here, periodic bc info should come from neighbors
             case(bc_noinflow)
                if (iw==1+1)then
                   do ix1=ixOmin1,ixOmax1
                     w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                        iw) = min(w(ixOmax1+1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                        iw),zero)
                   end do
                else
                   do ix1=ixOmin1,ixOmax1
                     w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw) = w(ixOmax1+1,&
                        ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw)
                   end do
                end if
             case (bc_aperiodic)
 !this just multiplies the variables with (-), they have been set from neighbors just like periodic.
                w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                   iw) = - w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                   iw)
             case (bc_data)
                ! skip it here, do AFTER all normal type boundaries are set
           case (bc_icarus)
              ! skip it here, do AFTER all normal type boundaries are set
             case (bc_character)
                ! skip it here, do AFTER all normal type boundaries are set
             case default
                write (unitterm,*) "Undefined boundarytype found in bc_phys",&
                    "for variable iw=",iw," and side iB=",iB
             end select
          end do
          
       end if
    case (2)
       if (iside==2) then
          ! maximal boundary
          iB=2*2
          ixOmin1=ixBmin1;ixOmin2=ixBmax2+1-nghostcells;ixOmin3=ixBmin3;
          ixOmax1=ixBmax1;ixOmax2=ixBmax2;ixOmax3=ixBmax3;
          ! cont/symm/asymm types
          do iw=1,nwflux+nwaux
             select case (typeboundary(iw,iB))
             case (bc_special)
                ! skip it here, do AFTER all normal type boundaries are set
             case (bc_cont)
                do ix2=ixOmin2,ixOmax2
                   w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,&
                      iw) = w(ixOmin1:ixOmax1,ixOmin2-1,ixOmin3:ixOmax3,iw)
                end do
             case (bc_symm)
                w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                   iw) = w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-nghostcells:-1,&
                   ixOmin3:ixOmax3,iw)
             case (bc_asymm)
                w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                   iw) =-w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-nghostcells:-1,&
                   ixOmin3:ixOmax3,iw)
             case (bc_periodic)
                ! skip it here, periodic bc info should come from neighbors
             case(bc_noinflow)
                if (iw==1+2)then
                  do ix2=ixOmin2,ixOmax2
                      w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,&
                         iw) = max(w(ixOmin1:ixOmax1,ixOmin2-1,ixOmin3:ixOmax3,&
                         iw),zero)
                  end do
                else
                  do ix2=ixOmin2,ixOmax2
                      w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,&
                         iw) = w(ixOmin1:ixOmax1,ixOmin2-1,ixOmin3:ixOmax3,iw)
                  end do
                end if
             case (bc_aperiodic)
 !this just multiplies the variables with (-), they have been set from neighbors just like periodic.
                w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                   iw) = - w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                   iw)
             case (bc_data)
                ! skip it here, do AFTER all normal type boundaries are set
           case (bc_icarus)
              ! skip it here, do AFTER all normal type boundaries are set
             case (bc_character)
                ! skip it here, do AFTER all normal type boundaries are set
             case default
                write (unitterm,*) "Undefined boundarytype found in bc_phys",&
                    "for variable iw=",iw," and side iB=",iB
             end select
          end do
          
       else
          ! minimal boundary
          iB=2*2-1
          ixOmin1=ixBmin1;ixOmin2=ixBmin2;ixOmin3=ixBmin3;
          ixOmax1=ixBmax1;ixOmax2=ixBmin2-1+nghostcells;ixOmax3=ixBmax3;
          ! cont/symm/asymm types
          do iw=1,nwflux+nwaux
             select case (typeboundary(iw,iB))
             case (bc_special)
                ! skip it here, do AFTER all normal type boundaries are set
             case (bc_cont)
                do ix2=ixOmin2,ixOmax2
                   w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,&
                      iw) = w(ixOmin1:ixOmax1,ixOmax2+1,ixOmin3:ixOmax3,iw)
                end do
             case (bc_symm)
                w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                   iw) = w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,&
                   ixOmin3:ixOmax3,iw)
             case (bc_asymm)
                w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                   iw) =-w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,&
                   ixOmin3:ixOmax3,iw)
             case (bc_periodic)
                ! skip it here, periodic bc info should come from neighbors
             case(bc_noinflow)
                if (iw==1+2)then
                   do ix2=ixOmin2,ixOmax2
                     w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,&
                        iw) = min(w(ixOmin1:ixOmax1,ixOmax2+1,ixOmin3:ixOmax3,&
                        iw),zero)
                   end do
                else
                   do ix2=ixOmin2,ixOmax2
                     w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,&
                        iw) = w(ixOmin1:ixOmax1,ixOmax2+1,ixOmin3:ixOmax3,iw)
                   end do
                end if
             case (bc_aperiodic)
 !this just multiplies the variables with (-), they have been set from neighbors just like periodic.
                w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                   iw) = - w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                   iw)
             case (bc_data)
                ! skip it here, do AFTER all normal type boundaries are set
           case (bc_icarus)
              ! skip it here, do AFTER all normal type boundaries are set
             case (bc_character)
                ! skip it here, do AFTER all normal type boundaries are set
             case default
                write (unitterm,*) "Undefined boundarytype found in bc_phys",&
                    "for variable iw=",iw," and side iB=",iB
             end select
          end do
          
       end if 
    case (3)
       if (iside==2) then
          ! maximal boundary
          iB=2*3
          ixOmin1=ixBmin1;ixOmin2=ixBmin2;ixOmin3=ixBmax3+1-nghostcells;
          ixOmax1=ixBmax1;ixOmax2=ixBmax2;ixOmax3=ixBmax3;
          ! cont/symm/asymm types
          do iw=1,nwflux+nwaux
             select case (typeboundary(iw,iB))
             case (bc_special)
                ! skip it here, do AFTER all normal type boundaries are set
             case (bc_cont)
                do ix3=ixOmin3,ixOmax3
                   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3,&
                      iw) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3-1,iw)
                end do
             case (bc_symm)
                w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                   iw) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                   ixOmin3-1:ixOmin3-nghostcells:-1,iw)
             case (bc_asymm)
                w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                   iw) =-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                   ixOmin3-1:ixOmin3-nghostcells:-1,iw)
             case (bc_periodic)
                ! skip it here, periodic bc info should come from neighbors
             case(bc_noinflow)
                if (iw==1+3)then
                  do ix3=ixOmin3,ixOmax3
                      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3,&
                         iw) = max(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3-1,&
                         iw),zero)
                  end do
                else
                  do ix3=ixOmin3,ixOmax3
                      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3,&
                         iw) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3-1,iw)
                  end do
                end if
             case (bc_aperiodic)
 !this just multiplies the variables with (-), they have been set from neighbors just like periodic.
                w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                   iw) = - w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                   iw)
             case (bc_data)
                ! skip it here, do AFTER all normal type boundaries are set
           case (bc_icarus)
              ! skip it here, do AFTER all normal type boundaries are set
             case (bc_character)
                ! skip it here, do AFTER all normal type boundaries are set
             case default
                write (unitterm,*) "Undefined boundarytype found in bc_phys",&
                    "for variable iw=",iw," and side iB=",iB
             end select
          end do
          
       else
          ! minimal boundary
          iB=2*3-1
          ixOmin1=ixBmin1;ixOmin2=ixBmin2;ixOmin3=ixBmin3;
          ixOmax1=ixBmax1;ixOmax2=ixBmax2;ixOmax3=ixBmin3-1+nghostcells;
          ! cont/symm/asymm types
          do iw=1,nwflux+nwaux
             select case (typeboundary(iw,iB))
             case (bc_special)
                ! skip it here, do AFTER all normal type boundaries are set
             case (bc_cont)
                do ix3=ixOmin3,ixOmax3
                   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3,&
                      iw) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmax3+1,iw)
                end do
             case (bc_symm)
                w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                   iw) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                   ixOmax3+nghostcells:ixOmax3+1:-1,iw)
             case (bc_asymm)
                w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                   iw) =-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                   ixOmax3+nghostcells:ixOmax3+1:-1,iw)
             case (bc_periodic)
                ! skip it here, periodic bc info should come from neighbors
             case(bc_noinflow)
                if (iw==1+3)then
                   do ix3=ixOmin3,ixOmax3
                     w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3,&
                        iw) = min(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmax3+1,&
                        iw),zero)
                   end do
                else
                   do ix3=ixOmin3,ixOmax3
                     w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3,&
                        iw) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmax3+1,iw)
                   end do
                end if
             case (bc_aperiodic)
 !this just multiplies the variables with (-), they have been set from neighbors just like periodic.
                w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                   iw) = - w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                   iw)
             case (bc_data)
                ! skip it here, do AFTER all normal type boundaries are set
           case (bc_icarus)
              ! skip it here, do AFTER all normal type boundaries are set
             case (bc_character)
                ! skip it here, do AFTER all normal type boundaries are set
             case default
                write (unitterm,*) "Undefined boundarytype found in bc_phys",&
                    "for variable iw=",iw," and side iB=",iB
             end select
          end do
          
       end if 
    end select

#:if defined('SPECIALBOUNDARY')    
    ! do user defined special boundary conditions
    if (any(typeboundary(1:nwflux+nwaux,iB)==bc_special)) then

       call specialbound_usr(time,ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,&
            ixGmax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,iB,w,x)

    end if
#:endif
    
  end associate
  end subroutine bc_phys

  !> fill inner boundary values
  subroutine getintbc(time,ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3)
    use mod_usr_methods, only: usr_internal_bc
    use mod_global_parameters

    double precision, intent(in)   :: time
    integer, intent(in)            :: ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,&
       ixGmax3
  
    integer :: iigrid, igrid, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
  
    ixOmin1=ixGmin1+nghostcells;ixOmin2=ixGmin2+nghostcells
    ixOmin3=ixGmin3+nghostcells;ixOmax1=ixGmax1-nghostcells
    ixOmax2=ixGmax2-nghostcells;ixOmax3=ixGmax3-nghostcells;

    !$OMP PARALLEL DO SCHEDULE(dynamic) PRIVATE(igrid)
    do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
       block=>ps(igrid)
       dxlevel(1)=rnode(rpdx1_,igrid);dxlevel(2)=rnode(rpdx2_,igrid)
       dxlevel(3)=rnode(rpdx3_,igrid);

       if (associated(usr_internal_bc)) then
          call usr_internal_bc(node(plevel_,igrid),time,ixGmin1,ixGmin2,&
             ixGmin3,ixGmax1,ixGmax2,ixGmax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
             ixOmax2,ixOmax3,ps(igrid)%w,ps(igrid)%x)
       end if
    end do
    !$OMP END PARALLEL DO

  end subroutine getintbc

end module mod_boundary_conditions
