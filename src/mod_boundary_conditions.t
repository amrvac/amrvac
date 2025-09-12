module mod_boundary_conditions

  implicit none
  private

  public :: bc_phys
  public :: getintbc

contains

  !> fill ghost cells at a physical boundary
  subroutine bc_phys(iside,idims,time,qdt,s,ixG^L,ixB^L)
    use mod_usr_methods, only: usr_special_bc
    use mod_bc_data, only: bc_data_set
    use mod_global_parameters
    use mod_physics
    use mod_functions_bfield, only: get_divb
    use mod_comm_lib, only: mpistop

    integer, intent(in) :: iside, idims, ixG^L,ixB^L
    double precision, intent(in) :: time,qdt
    type(state), intent(inout) :: s

    double precision :: Q(ixG^S),Qp(ixG^S) 
    integer :: idir, ixOs^L,hxO^L,jxO^L
    integer :: iw, iB, ix^D, ixO^L

    associate(x=>s%x,w=>s%w,ws=>s%ws)
    select case (idims)
    {case (^D)
       if (iside==2) then
          ! maximal boundary
          iB=2*^D
          ixOmin^DD=ixBmax^D+1-nghostcells^D%ixOmin^DD=ixBmin^DD;
          ixOmax^DD=ixBmax^DD;
          ! cont/symm/asymm types
          do iw=1,nwflux+nwaux
             select case (typeboundary(iw,iB))
             case (bc_special)
                ! skip it here, do AFTER all normal type boundaries are set
             case (bc_cont)
                do ix^D=ixOmin^D,ixOmax^D
                   w(ix^D^D%ixO^S,iw) = w(ixOmin^D-1^D%ixO^S,iw)
                end do
             case (bc_symm)
                w(ixO^S,iw) = w(ixOmin^D-1:ixOmin^D-nghostcells:-1^D%ixO^S,iw)
             case (bc_asymm)
                w(ixO^S,iw) =-w(ixOmin^D-1:ixOmin^D-nghostcells:-1^D%ixO^S,iw)
             case (bc_periodic)
                ! skip it here, periodic bc info should come from neighbors
             case(bc_noinflow)
                if (iw==1+^D)then
                  do ix^D=ixOmin^D,ixOmax^D
                      w(ix^D^D%ixO^S,iw) = max(w(ixOmin^D-1^D%ixO^S,iw),zero)
                  end do
                else
                  do ix^D=ixOmin^D,ixOmax^D
                      w(ix^D^D%ixO^S,iw) = w(ixOmin^D-1^D%ixO^S,iw)
                  end do
                end if
             case (bc_aperiodic)
                !this just multiplies the variables with (-), they have been set from neighbors just like periodic.
                w(ixO^S,iw) = - w(ixO^S,iw)
             case (bc_data)
                ! skip it here, do AFTER all normal type boundaries are set
           case (bc_icarus)
              ! skip it here, do AFTER all normal type boundaries are set
             case (bc_character)
                ! skip it here, do AFTER all normal type boundaries are set
             case default
                write (unitterm,*) "Undefined boundarytype found in bc_phys", &
                   "for variable iw=",iw," and side iB=",iB
             end select
          end do
          if(stagger_grid) then
            do idir=1,nws
            ! At this stage, extrapolation is applied only to the tangential components
              if(idir==^D) cycle 
              ixOsmax^DD=ixOmax^DD;
              ixOsmin^DD=ixOmin^DD-kr(^DD,idir);
              select case(typeboundary(iw_mag(idir),iB))
              case (bc_special)
                 ! skip it here, do AFTER all normal type boundaries are set
              case (bc_cont)
                do ix^D=ixOsmin^D,ixOsmax^D
                   ws(ix^D^D%ixOs^S,idir) = ws(ixOsmin^D-1^D%ixOs^S,idir)
                end do
              case (bc_symm)
                ws(ixOs^S,idir) = ws(ixOsmin^D-1:ixOsmin^D-nghostcells:-1^D%ixOs^S,idir)
              case (bc_asymm)
                ws(ixOs^S,idir) =-ws(ixOsmin^D-1:ixOsmin^D-nghostcells:-1^D%ixOs^S,idir)
              case (bc_periodic)
              case default
                write (unitterm,*) "Undefined boundarytype found in bc_phys", &
                   "for variable iw=",iw," and side iB=",iB
              end select
            end do
            ! Now that the tangential components are set,
            ! fill the normal components using a prescription for the divergence.
            ! This prescription is given by the typeB for the normal component.
            do idir=1,nws
              ! Consider only normal direction
              if (idir/=^D) cycle
              ixOs^L=ixO^L;
              select case(typeboundary(iw_mag(idir),iB))
              case(bc_cont)
                hxO^L=ixO^L-nghostcells*kr(^DD,^D);
                ! Calculate divergence and partial divergence
                call get_divb(s%w,ixG^L,hxO^L,Q)
                ws(ixOs^S,idir)=zero
                do ix^D=0,nghostcells-1
                  call get_divb(s%w,ixG^L,ixO^L,Qp)
                  ws(ixOsmin^D+ix^D^D%ixOs^S,idir)=&
                    (Q(hxOmax^D^D%hxO^S)*s%dvolume(hxOmax^D^D%hxO^S)&
                   -Qp(ixOmin^D+ix^D^D%ixO^S)*s%dvolume(ixOmin^D+ix^D^D%ixO^S))&
                    /s%surfaceC(ixOsmin^D+ix^D^D%ixOs^S,^D)
                end do
                ! Fill cell averages
                call phys_face_to_center(ixO^L,s)
              case(bc_symm)
                ! (a)symmetric normal B ensures symmetric divB
                ws(ixOs^S,idir)= ws(ixOsmin^D-2:ixOsmin^D-nghostcells-1:-1^D%ixOs^S,idir)
                ! Fill cell averages
                call phys_face_to_center(ixO^L,s)
              case(bc_asymm)
                ! (a)symmetric normal B ensures symmetric divB
                ws(ixOs^S,idir)=-ws(ixOsmin^D-2:ixOsmin^D-nghostcells-1:-1^D%ixOs^S,idir)
                ! Fill cell averages
                call phys_face_to_center(ixO^L,s)
              case(bc_periodic)
                ! Fill cell averages
                call phys_face_to_center(ixO^L,s)
              end select
            end do
          end if
       else
          ! minimal boundary
          iB=2*^D-1
          ixOmin^DD=ixBmin^DD;
          ixOmax^DD=ixBmin^D-1+nghostcells^D%ixOmax^DD=ixBmax^DD;
          ! cont/symm/asymm types
          do iw=1,nwflux+nwaux
             select case (typeboundary(iw,iB))
             case (bc_special)
                ! skip it here, do AFTER all normal type boundaries are set
             case (bc_cont)
                do ix^D=ixOmin^D,ixOmax^D
                   w(ix^D^D%ixO^S,iw) = w(ixOmax^D+1^D%ixO^S,iw)
                end do
             case (bc_symm)
                w(ixO^S,iw) = w(ixOmax^D+nghostcells:ixOmax^D+1:-1^D%ixO^S,iw)
             case (bc_asymm)
                w(ixO^S,iw) =-w(ixOmax^D+nghostcells:ixOmax^D+1:-1^D%ixO^S,iw)
             case (bc_periodic)
                ! skip it here, periodic bc info should come from neighbors
             case(bc_noinflow)
                if (iw==1+^D)then
                   do ix^D=ixOmin^D,ixOmax^D
                     w(ix^D^D%ixO^S,iw) = min(w(ixOmax^D+1^D%ixO^S,iw),zero)
                   end do
                else
                   do ix^D=ixOmin^D,ixOmax^D
                     w(ix^D^D%ixO^S,iw) = w(ixOmax^D+1^D%ixO^S,iw)
                   end do
                end if
             case (bc_aperiodic)
                !this just multiplies the variables with (-), they have been set from neighbors just like periodic.
                w(ixO^S,iw) = - w(ixO^S,iw)
             case (bc_data)
                ! skip it here, do AFTER all normal type boundaries are set
           case (bc_icarus)
              ! skip it here, do AFTER all normal type boundaries are set
             case (bc_character)
                ! skip it here, do AFTER all normal type boundaries are set
             case default
                write (unitterm,*) "Undefined boundarytype found in bc_phys", &
                   "for variable iw=",iw," and side iB=",iB
             end select
          end do
          if(stagger_grid) then
            do idir=1,nws
            ! At this stage, extrapolation is applied only to the tangential components
              if(idir==^D) cycle 
              ixOsmax^DD=ixOmax^DD;
              ixOsmin^DD=ixOmin^DD-kr(^DD,idir);
              select case(typeboundary(iw_mag(idir),iB))
              case (bc_special)
                ! skip it here, periodic bc info should come from neighbors
              case (bc_cont)
                do ix^D=ixOsmin^D,ixOsmax^D
                   ws(ix^D^D%ixOs^S,idir) = ws(ixOsmax^D+1^D%ixOs^S,idir)
                end do
              case (bc_symm)
                ws(ixOs^S,idir) = ws(ixOsmax^D+nghostcells:ixOsmax^D+1:-1^D%ixOs^S,idir)
              case (bc_asymm)
                ws(ixOs^S,idir) =-ws(ixOsmax^D+nghostcells:ixOsmax^D+1:-1^D%ixOs^S,idir)
              case (bc_periodic)
              case default
                write (unitterm,*) "Undefined boundarytype in bc_phys", &
                   "for variable iw=",iw," and side iB=",iB
              end select
            end do
            ! Now that the tangential components are set,
            ! fill the normal components using a prescription for the divergence.
            ! This prescription is given by the typeB for the normal component.
            do idir=1,nws
              ! Consider only normal direction
              if (idir/=^D) cycle
              ixOs^L=ixO^L-kr(^DD,^D);
              select case(typeboundary(iw_mag(idir),iB))
              case(bc_cont)
                jxO^L=ixO^L+nghostcells*kr(^DD,^D);
                ! Calculate divergence and partial divergence
                call get_divb(s%w,ixG^L,jxO^L,Q)
                ws(ixOs^S,idir)=zero
                do ix^D=0,nghostcells-1
                  call get_divb(s%w,ixG^L,ixO^L,Qp)
                  ws(ixOsmax^D-ix^D^D%ixOs^S,idir)=&
                   -(Q(jxOmin^D^D%jxO^S)*s%dvolume(jxOmin^D^D%jxO^S)&
                   -Qp(ixOmax^D-ix^D^D%ixO^S)*s%dvolume(ixOmax^D-ix^D^D%ixO^S))&
                    /s%surfaceC(ixOsmax^D-ix^D^D%ixOs^S,^D)
                end do
                ! Fill cell averages
                call phys_face_to_center(ixO^L,s)
              case(bc_symm)
                ! (a)symmetric normal B ensures symmetric divB
                ws(ixOs^S,idir)= ws(ixOsmax^D+nghostcells+1:ixOsmax^D+2:-1^D%ixOs^S,idir)
                ! Fill cell averages
                call phys_face_to_center(ixO^L,s)
              case(bc_asymm)
                ! (a)symmetric normal B ensures symmetric divB
                ws(ixOs^S,idir)=-ws(ixOsmax^D+nghostcells+1:ixOsmax^D+2:-1^D%ixOs^S,idir)
                ! Fill cell averages
                call phys_face_to_center(ixO^L,s)
              case(bc_periodic)
                ! Fill cell averages
                call phys_face_to_center(ixO^L,s)
              end select
            end do
          end if
       end if \}
    end select

    ! do user defined special boundary conditions
    if (any(typeboundary(1:nwflux+nwaux,iB)==bc_special)) then
       if (.not. associated(usr_special_bc)) &
            call mpistop("usr_special_bc not defined")
       call usr_special_bc(time,ixG^L,ixO^L,iB,w,x)
    end if

    ! fill boundary conditions from external data vtk files
    if (any(typeboundary(1:nwflux+nwaux,iB)==bc_data)) then
       call bc_data_set(time,ixG^L,ixO^L,iB,w,x)
    end if

    ! fill boundary conditions from external data vtk files and do user defined special boundary conditions
    if (any(typeboundary(1:nwflux+nwaux,iB)==bc_icarus)) then
       call bc_data_set(time,ixG^L,ixO^L,iB,w,x)
       if (.not. associated(usr_special_bc)) &
            call mpistop("usr_special_bc not defined")
       call usr_special_bc(time,ixG^L,ixO^L,iB,w,x)
    end if

    end associate
  end subroutine bc_phys

  !> fill inner boundary values
  subroutine getintbc(time,ixG^L)
    use mod_usr_methods, only: usr_internal_bc
    use mod_global_parameters

    double precision, intent(in)   :: time
    integer, intent(in)            :: ixG^L
  
    integer :: iigrid, igrid, ixO^L
  
    ixO^L=ixG^L^LSUBnghostcells;

    !$OMP PARALLEL DO SCHEDULE(dynamic) PRIVATE(igrid)
    do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
       block=>ps(igrid)
       ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);

       if (associated(usr_internal_bc)) then
          call usr_internal_bc(node(plevel_,igrid),time,ixG^L,ixO^L,ps(igrid)%w,ps(igrid)%x)
       end if
    end do
    !$OMP END PARALLEL DO

  end subroutine getintbc

end module mod_boundary_conditions
