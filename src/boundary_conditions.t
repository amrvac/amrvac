!> fill ghost cells at a physical boundary
subroutine bc_init()
  use mod_global_parameters
  use mod_physics
  use mod_geometry
  use mod_usr_methods, only: usr_set_bc

  integer :: iw, ib

  if (nwflux > 0) then
     allocate(bc_type%w(nwflux, 2 * ndim))
     bc_type%w = 0
  end if
  if (nmetric > 0) then
     allocate(bc_type%m(nmetric, 2 * ndim))
     bc_type%m = 0
  end if

  ! user specify boundary type
  if (.not. associated(usr_set_bc)) then
     call mpistop("usr_set_bc not defined")
  else
     call usr_set_bc()
  end if

  if (any(bc_type%m == 0)) then
    if (mype==0) then
      do ib = 1, 2*ndim; do iw = 1, nmetric
         if (bc_type%m(iw,ib) == 0) then
            write(*,*) "Undefined boundary type for metric: iw=",iw,&
               " and side iB=", iB
         end if
      end do; end do
    end if
    call mpistop("Not all boundary conditions for metric have been defined")
  end if

  do idim=1,ndim
     periodB(idim)=(any(bc_type%m(:,2*idim-1:2*idim)==bc_periodic))
     aperiodB(idim)=(any(bc_type%m(:,2*idim-1:2*idim)==bc_aperiodic))
     if (periodB(idim).or.aperiodB(idim)) then
        do iw=1,nmetric
           if (bc_type%m(iw,2*idim-1) .ne. bc_type%m(iw,2*idim)) &
                call mpistop("Wrong counterpart in periodic boundary")

           if (bc_type%m(iw,2*idim-1) /= bc_periodic .and. &
                bc_type%m(iw,2*idim-1) /= bc_aperiodic) then
             call mpistop("Each dimension should either have all "//&
                  "or no variables periodic, some can be aperiodic")
           end if
        end do
     end if
  end do
  ! fixme: refine this part
  {^NOONED
  do idim=1,ndim
    if(any(bc_type%m(:,2*idim-1)==bc_pole)) then
      if(any(bc_type%m(:,2*idim-1)/=bc_pole)) bc_type%m(:,2*idim-1)=bc_pole
        windex=2
      bc_type%m(:,2*idim-1)=bc_symm
      select case(coordinate)
      case(cylindrical)
        bc_type%m(phi_+1,2*idim-1)=bc_asymm
        !if(physics_type=='mhd') bc_type%m(ndir+windex+phi_,2*idim-1)=bc_asymm
      case(spherical)
        bc_type%m(3:ndir+1,2*idim-1)=bc_asymm
        !if(physics_type=='mhd') bc_type%m(ndir+windex+2:ndir+windex+ndir,2*idim-1)=bc_asymm
      case default
        call mpistop('Pole is in cylindrical, polar, spherical coordinates!')
      end select
    end if
    if(any(bc_type%m(:,2*idim)==bc_pole)) then
      if(any(bc_type%m(:,2*idim)/=bc_pole)) bc_type%m(:,2*idim)=bc_pole
        windex=2
      bc_type%m(:,2*idim)=bc_symm
      select case(coordinate)
      case(cylindrical)
        bc_type%m(phi_+1,2*idim)=bc_asymm
        !if(physics_type=='mhd') bc_type%m(ndir+windex+phi_,2*idim)=bc_asymm
      case(spherical)
        bc_type%m(3:ndir+1,2*idim)=bc_asymm
        !if(physics_type=='mhd') bc_type%m(ndir+windex+2:ndir+windex+ndir,2*idim)=bc_asymm
      case default
        call mpistop('Pole is in cylindrical, polar, spherical coordinates!')
      end select
    end if
  end do
  }

  if (any(bc_type%w == 0)) then
    if (mype==0) then
      do ib = 1, 2*ndim; do iw = 1, nwflux
         if (bc_type%w(iw,ib) == 0) then
            write(*,*) "Undefined boundary type for w: iw=",iw,&
               " and side iB=", iB
         end if
      end do; end do
    end if
    call mpistop("Not all boundary conditions for w have been defined")
  end if

  do idim=1,ndim
     periodB(idim)=(any(bc_type%w(:,2*idim-1:2*idim)==bc_periodic))
     aperiodB(idim)=(any(bc_type%w(:,2*idim-1:2*idim)==bc_aperiodic))
     if (periodB(idim).or.aperiodB(idim)) then
        do iw=1,nwflux
           if (bc_type%w(iw,2*idim-1) .ne. bc_type%w(iw,2*idim)) &
                call mpistop("Wrong counterpart in periodic boundary")

           if (bc_type%w(iw,2*idim-1) /= bc_periodic .and. &
                bc_type%w(iw,2*idim-1) /= bc_aperiodic) then
             call mpistop("Each dimension should either have all "//&
                  "or no variables periodic, some can be aperiodic")
           end if
        end do
     end if
  end do
  ! fixme: refine this part
  {^NOONED
  do idim=1,ndim
    if(any(bc_type%w(:,2*idim-1)==bc_pole)) then
      if(any(bc_type%w(:,2*idim-1)/=bc_pole)) bc_type%w(:,2*idim-1)=bc_pole
        windex=2
      bc_type%w(:,2*idim-1)=bc_symm
      select case(coordinate)
      case(cylindrical)
        bc_type%w(phi_+1,2*idim-1)=bc_asymm
        !if(physics_type=='mhd') bc_type%w(ndir+windex+phi_,2*idim-1)=bc_asymm
      case(spherical)
        bc_type%w(3:ndir+1,2*idim-1)=bc_asymm
        !if(physics_type=='mhd') bc_type%w(ndir+windex+2:ndir+windex+ndir,2*idim-1)=bc_asymm
      case default
        call mpistop('Pole is in cylindrical, polar, spherical coordinates!')
      end select
    end if
    if(any(bc_type%w(:,2*idim)==bc_pole)) then
      if(any(bc_type%w(:,2*idim)/=bc_pole)) bc_type%w(:,2*idim)=bc_pole
        windex=2
      bc_type%w(:,2*idim)=bc_symm
      select case(coordinate)
      case(cylindrical)
        bc_type%w(phi_+1,2*idim)=bc_asymm
        !if(physics_type=='mhd') bc_type%w(ndir+windex+phi_,2*idim)=bc_asymm
      case(spherical)
        bc_type%w(3:ndir+1,2*idim)=bc_asymm
        !if(physics_type=='mhd') bc_type%w(ndir+windex+2:ndir+windex+ndir,2*idim)=bc_asymm
      case default
        call mpistop('Pole is in cylindrical, polar, spherical coordinates!')
      end select
    end if
  end do
  }


end subroutine bc_init

!> fill ghost cells at a physical boundary
subroutine bc_metric(iside,idims,time,qdt,s,ixG^L,ixB^L)
  use mod_global_parameters
  use mod_physics

  integer, intent(in)          :: iside, idims, ixG^L, ixB^L
  double precision, intent(in) :: time,qdt
  type(state), intent(inout)   :: s

  integer          :: iw, iB, ix^D, ixO^L

  associate(x=>s%mesh%x,w=>s%metric%vars)
  select case (idims)
  {case (^D)
     if (iside==2) then
        ! maximal boundary
        iB=2*^D
        ixOmin^DD=ixBmax^D+1-nghostcells^D%ixOmin^DD=ixBmin^DD;
        ixOmax^DD=ixBmax^DD;
        ! cont/symm/asymm types
        do iw=1,nmetric
           select case (bc_type%m(iw,iB))
           case (bc_symm)
              w(ixO^S,iw) = w(ixOmin^D-1:ixOmin^D-nghostcells:-1^D%ixO^S,iw)
           case (bc_asymm)
              w(ixO^S,iw) =-w(ixOmin^D-1:ixOmin^D-nghostcells:-1^D%ixO^S,iw)
           case (bc_cont, bc_noinflow)
              do ix^D=ixOmin^D,ixOmax^D
                 w(ix^D^D%ixO^S,iw) = w(ixOmin^D-1^D%ixO^S,iw)
              end do
           case default
              write (unitterm,*) "Undefined boundarytype found in bc_metric ", &
                 "for variable iw=",iw," and side iB=",iB
           end select
        end do
     else
        ! minimal boundary
        iB=2*^D-1
        ixOmin^DD=ixBmin^DD;
        ixOmax^DD=ixBmin^D-1+nghostcells^D%ixOmax^DD=ixBmax^DD;
        ! cont/symm/asymm types
        do iw=1,nmetric
           select case (bc_type%m(iw,iB))
           case (bc_symm)
              w(ixO^S,iw) = w(ixOmax^D+nghostcells:ixOmax^D+1:-1^D%ixO^S,iw)
           case (bc_asymm)
              w(ixO^S,iw) =-w(ixOmax^D+nghostcells:ixOmax^D+1:-1^D%ixO^S,iw)
           case (bc_cont, bc_noinflow)
              do ix^D=ixOmin^D,ixOmax^D
                 w(ix^D^D%ixO^S,iw) = w(ixOmax^D+1^D%ixO^S,iw)
              end do
           case default
              write (unitterm,*) "Undefined boundarytype found in bc_metric ", &
                 "for variable iw=",iw," and side iB=",iB
           end select
        end do
     end if \}
  end select

  end associate
end subroutine bc_metric

!> fill ghost cells at a physical boundary
subroutine bc_phys(iside,idims,time,qdt,s,ixG^L,ixB^L)
  use mod_usr_methods, only: usr_special_bc
  use mod_bc_data, only: bc_data_set
  use mod_global_parameters
  use mod_physics

  integer, intent(in)          :: iside, idims, ixG^L,ixB^L
  double precision, intent(in) :: time,qdt
  type(state), intent(inout)   :: s

  double precision :: wtmp(ixG^S,1:nwflux)
  integer          :: idir, is
  integer          :: ixOs^L,hxO^L,jxO^L
  double precision :: Q(ixG^S),Qp(ixG^S)
  integer          :: iw, iB, ix^D, ixO^L, ixM^L, nghostcellsi,iib^D
  logical          :: isphysbound
  integer          :: vel_index

  associate(x=>s%mesh%x,w=>s%w,ws=>s%ws)
  select case (idims)
  {case (^D)
     if (iside==2) then
        ! maximal boundary
        iB=2*^D
        ixOmin^DD=ixBmax^D+1-nghostcells^D%ixOmin^DD=ixBmin^DD;
        ixOmax^DD=ixBmax^DD;
        ! cont/symm/asymm types
        do iw=1,nwflux
           select case (bc_type%w(iw,iB))
           case (bc_symm)
              w(ixO^S,iw) = w(ixOmin^D-1:ixOmin^D-nghostcells:-1^D%ixO^S,iw)
           case (bc_asymm)
              w(ixO^S,iw) =-w(ixOmin^D-1:ixOmin^D-nghostcells:-1^D%ixO^S,iw)
           case (bc_cont)
              do ix^D=ixOmin^D,ixOmax^D
                 w(ix^D^D%ixO^S,iw) = w(ixOmin^D-1^D%ixO^S,iw)
              end do
           case(bc_noinflow)
              ! noinflow applies only on velocities, for other vectors, use cont
              ! Note that although ^D = ndim is not ndir, in case of ndir > ndim, 
              ! the lastest velocity has no corresponding boundary
              ! so apply cont should be enough.
              vel_index = W_vel(^D)

              if ( iw==vel_index )then
                do ix^D=ixOmin^D,ixOmax^D
                    w(ix^D^D%ixO^S,iw) = max(w(ixOmin^D-1^D%ixO^S,iw),zero)
                end do
              else
                do ix^D=ixOmin^D,ixOmax^D
                    w(ix^D^D%ixO^S,iw) = w(ixOmin^D-1^D%ixO^S,iw)
                end do
              end if
           case (bc_periodic)
              ! skip it here, periodic bc info should come from neighbors
           case (bc_aperiodic)
              !this just multiplies the variables with (-), they have been set from neighbors just like periodic.
              w(ixO^S,iw) = - w(ixO^S,iw)
           case (bc_special, bc_data, bc_character)
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
            select case(bc_type%w(Bvec(idir),iB))
            case (bc_symm)
              ws(ixOs^S,idir) = ws(ixOsmin^D-1:ixOsmin^D-nghostcells:-1^D%ixOs^S,idir)
            case (bc_asymm)
              ws(ixOs^S,idir) =-ws(ixOsmin^D-1:ixOsmin^D-nghostcells:-1^D%ixOs^S,idir)
            case (bc_cont,bc_noinflow)
              do ix^D=ixOsmin^D,ixOsmax^D
                 ws(ix^D^D%ixOs^S,idir) = ws(ixOsmin^D-1^D%ixOs^S,idir)
              end do
            case (bc_periodic, bc_special)
               ! skip it here, do AFTER all normal type boundaries are set
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
            select case(bc_type%w(Bvec(idir),iB))
            case(bc_symm)
              ! (a)symmetric normal B ensures symmetric divB
              ws(ixOs^S,idir)= ws(ixOsmin^D-2:ixOsmin^D-nghostcells-1:-1^D%ixOs^S,idir)
            case(bc_asymm)
              ! (a)symmetric normal B ensures symmetric divB
              ws(ixOs^S,idir)=-ws(ixOsmin^D-2:ixOsmin^D-nghostcells-1:-1^D%ixOs^S,idir)
            case(bc_cont, bc_noinflow)
              hxO^L=ixO^L-nghostcells*kr(^DD,^D);
              ! Calculate divergence and partial divergence
              call phys_get_divb(ixG^L,hxO^L,s,Q)
              ws(ixOs^S,idir)=zero
              do ix^D=0,nghostcells-1
                call phys_get_divb(ixG^L,ixO^L,s,Qp)
                ws(ixOsmin^D+ix^D^D%ixOs^S,idir)=&
                  (Q(hxOmax^D^D%hxO^S)*s%mesh%dvolume(hxOmax^D^D%hxO^S)&
                 -Qp(ixOmin^D+ix^D^D%ixO^S)*s%mesh%dvolume(ixOmin^D+ix^D^D%ixO^S))&
                  /s%mesh%surfaceC(ixOsmin^D+ix^D^D%ixOs^S,^D)
              end do
            case(bc_periodic)
            end select
          end do
          ! Fill cell averages
          !call phys_face_to_center(ixO^L,s)
        end if
     else
        ! minimal boundary
        iB=2*^D-1
        ixOmin^DD=ixBmin^DD;
        ixOmax^DD=ixBmin^D-1+nghostcells^D%ixOmax^DD=ixBmax^DD;
        ! cont/symm/asymm types
        do iw=1,nwflux
           select case (bc_type%w(iw,iB))
           case (bc_symm)
              w(ixO^S,iw) = w(ixOmax^D+nghostcells:ixOmax^D+1:-1^D%ixO^S,iw)
           case (bc_asymm)
              w(ixO^S,iw) =-w(ixOmax^D+nghostcells:ixOmax^D+1:-1^D%ixO^S,iw)
           case (bc_cont)
              do ix^D=ixOmin^D,ixOmax^D
                 w(ix^D^D%ixO^S,iw) = w(ixOmax^D+1^D%ixO^S,iw)
              end do
           case(bc_noinflow)
              vel_index = W_vel(^D)

              if ( iw==vel_index )then
                 do ix^D=ixOmin^D,ixOmax^D
                   w(ix^D^D%ixO^S,iw) = min(w(ixOmax^D+1^D%ixO^S,iw),zero)
                 end do
              else
                 do ix^D=ixOmin^D,ixOmax^D
                   w(ix^D^D%ixO^S,iw) = w(ixOmax^D+1^D%ixO^S,iw)
                 end do
              end if
           case (bc_special)
              ! skip it here, do AFTER all normal type boundaries are set
           case (bc_periodic)
              ! skip it here, periodic bc info should come from neighbors
           case (bc_aperiodic)
              !this just multiplies the variables with (-), they have been set from neighbors just like periodic.
              w(ixO^S,iw) = - w(ixO^S,iw)
           case (bc_data)
              ! skip it here, do AFTER all normal type boundaries are set
           case (bc_character)
              ! skip it here, do AFTER all normal type boundaries are set
           case default
              write (unitterm,*) "Undefined boundarytype found in bc_phys ", &
                 "for variable iw=",iw," and side iB=",iB
           end select
        end do
        if(stagger_grid) then
          do idir=1,nws
          ! At this stage, extrapolation is applied only to the tangential components
            if(idir==^D) cycle 
            ixOsmax^DD=ixOmax^DD;
            ixOsmin^DD=ixOmin^DD-kr(^DD,idir);
            select case(bc_type%w(Bvec(idir),iB))
            case (bc_symm)
              ws(ixOs^S,idir) = ws(ixOsmax^D+nghostcells:ixOsmax^D+1:-1^D%ixOs^S,idir)
            case (bc_asymm)
              ws(ixOs^S,idir) =-ws(ixOsmax^D+nghostcells:ixOsmax^D+1:-1^D%ixOs^S,idir)
            case (bc_cont, bc_noinflow)
              do ix^D=ixOsmin^D,ixOsmax^D
                 ws(ix^D^D%ixOs^S,idir) = ws(ixOsmax^D+1^D%ixOs^S,idir)
              end do
            case (bc_periodic)
            case (bc_special)
              ! skip it here, periodic bc info should come from neighbors
            case default
              write (unitterm,*) "Undefined boundarytype found in bc_phys ", &
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
            select case(bc_type%w(Bvec(idir),iB))
            case(bc_symm)
              ! (a)symmetric normal B ensures symmetric divB
              ws(ixOs^S,idir)= ws(ixOsmax^D+nghostcells+1:ixOsmax^D+2:-1^D%ixOs^S,idir)
            case(bc_asymm)
              ! (a)symmetric normal B ensures symmetric divB
              ws(ixOs^S,idir)=-ws(ixOsmax^D+nghostcells+1:ixOsmax^D+2:-1^D%ixOs^S,idir)
            case(bc_cont, bc_noinflow)
              jxO^L=ixO^L+nghostcells*kr(^DD,^D);
              ! Calculate divergence and partial divergence
              call phys_get_divb(ixG^L,jxO^L,s,Q)
              ws(ixOs^S,idir)=zero
              do ix^D=0,nghostcells-1
                call phys_get_divb(ixG^L,ixO^L,s,Qp)
                ws(ixOsmax^D-ix^D^D%ixOs^S,idir)=&
                 -(Q(jxOmin^D^D%jxO^S)*s%mesh%dvolume(jxOmin^D^D%jxO^S)&
                 -Qp(ixOmax^D-ix^D^D%ixO^S)*s%mesh%dvolume(ixOmax^D-ix^D^D%ixO^S))&
                  /s%mesh%surfaceC(ixOsmax^D-ix^D^D%ixOs^S,^D)
              end do
            case(bc_periodic)
            end select
          end do
          ! Fill cell averages
          !call phys_face_to_center(ixO^L,s)
        end if
     end if \}
  end select

  ! do user defined special boundary conditions
  if (any(bc_type%w(1:nwflux,iB)==bc_special)) then
     if (.not. associated(usr_special_bc)) &
          call mpistop("usr_special_bc not defined")
     call usr_special_bc(time,ixG^L,ixO^L,iB,s)
  end if

  ! fill boundary conditions from external data vtk files
  if (any(bc_type%w(1:nwflux,iB)==bc_data)) then
     call bc_data_set(time,ixG^L,ixO^L,iB,w,x)
  end if

  {#IFDEF EVOLVINGBOUNDARY
  if (any(bc_type%w(1:nwflux,iB)==bc_character)) then
    ixM^L=ixM^LL;
    if(ixGmax1==ixGhi1) then
      nghostcellsi=nghostcells
    else
      nghostcellsi=ceiling(nghostcells*0.5d0)
    end if
    select case (idims)
    {case (^D)
       if (iside==2) then
          ! maximal boundary
          ixOmin^DD=ixGmax^D+1-nghostcellsi^D%ixOmin^DD=ixBmin^DD;
          ixOmax^DD=ixBmax^DD;
          if(all(w(ixO^S,1:nwflux)==0.d0)) then
            do ix^D=ixOmin^D,ixOmax^D
               w(ix^D^D%ixO^S,1:nwflux) = w(ixOmin^D-1^D%ixO^S,1:nwflux)
            end do
          end if
          if(qdt>0.d0.and.ixGmax^D==ixGhi^D) then
            ixOmin^DD=ixOmin^D^D%ixOmin^DD=ixMmin^DD;
            ixOmax^DD=ixOmax^D^D%ixOmax^DD=ixMmax^DD;
            wtmp(ixG^S,1:nwflux)=pso(block%igrid)%w(ixG^S,1:nwflux)
            call characteristic_project(idims,iside,ixG^L,ixO^L,wtmp,x,dxlevel,qdt)
            w(ixO^S,1:nwflux)=wtmp(ixO^S,1:nwflux)
          end if
       else
          ! minimal boundary
          ixOmin^DD=ixBmin^DD;
          ixOmax^DD=ixGmin^D-1+nghostcellsi^D%ixOmax^DD=ixBmax^DD;
          if(all(w(ixO^S,1:nwflux)==0.d0)) then
            do ix^D=ixOmin^D,ixOmax^D
               w(ix^D^D%ixO^S,1:nwflux) = w(ixOmax^D+1^D%ixO^S,1:nwflux)
            end do
          end if
          if(qdt>0.d0.and.ixGmax^D==ixGhi^D) then
            ixOmin^DD=ixOmin^D^D%ixOmin^DD=ixMmin^DD;
            ixOmax^DD=ixOmax^D^D%ixOmax^DD=ixMmax^DD;
            wtmp(ixG^S,1:nwflux)=pso(block%igrid)%w(ixG^S,1:nwflux)
            call characteristic_project(idims,iside,ixG^L,ixO^L,wtmp,x,dxlevel,qdt)
            w(ixO^S,1:nwflux)=wtmp(ixO^S,1:nwflux)
          end if
       end if \}
    end select
    if(ixGmax1==ixGhi1) then
      call identifyphysbound(block%igrid,isphysbound,iib^D)   
      if(iib1==-1.and.iib2==-1) then
        do ix2=nghostcells,1,-1 
          do ix1=nghostcells,1,-1 
            w(ix^D,1:nwflux)=(w(ix1+1,ix2+1,1:nwflux)+w(ix1+1,ix2,1:nwflux)+w(ix1,ix2+1,1:nwflux))/3.d0
          end do
        end do
      end if
      if(iib1== 1.and.iib2==-1) then
        do ix2=nghostcells,1,-1 
          do ix1=ixMmax1+1,ixGmax1
            w(ix^D,1:nwflux)=(w(ix1-1,ix2+1,1:nwflux)+w(ix1-1,ix2,1:nwflux)+w(ix1,ix2+1,1:nwflux))/3.d0
          end do
        end do
      end if
      if(iib1==-1.and.iib2== 1) then
        do ix2=ixMmax2+1,ixGmax2
          do ix1=nghostcells,1,-1 
            w(ix^D,1:nwflux)=(w(ix1+1,ix2-1,1:nwflux)+w(ix1+1,ix2,1:nwflux)+w(ix1,ix2-1,1:nwflux))/3.d0
          end do
        end do
      end if
      if(iib1== 1.and.iib2== 1) then
        do ix2=ixMmax2+1,ixGmax2
          do ix1=ixMmax1+1,ixGmax1
            w(ix^D,1:nwflux)=(w(ix1-1,ix2-1,1:nwflux)+w(ix1-1,ix2,1:nwflux)+w(ix1,ix2-1,1:nwflux))/3.d0
          end do
        end do
      end if
    end if
  end if
  }
  !end do

  ! Fill cell averages
  if (stagger_grid) call phys_face_to_center(ixO^L,s)

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
        call usr_internal_bc(node(plevel_,igrid),time,ixG^L,ixO^L,ps(igrid))
     end if
  end do
  !$OMP END PARALLEL DO

end subroutine getintbc
