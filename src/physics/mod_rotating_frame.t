!> Module for including rotating frame in hydrodynamics simulations
!> The rotation vector is assumed to be along z direction
!>(both in cylindrical and spherical)

module mod_rotating_frame
  implicit none

  !> Rotation frequency of the frame
  double precision :: omega_frame

  !> Index of the density (in the w array)
  integer, private, parameter              :: rho_ = 1
  
contains
  !> Read this module's parameters from a file
  subroutine rotating_frame_params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n
    
    namelist /rotating_frame_list/ omega_frame
    
    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, rotating_frame_list, end=111)
111    close(unitpar)
    end do
    
  end subroutine rotating_frame_params_read
  
  !> Initialize the module
  subroutine rotating_frame_init()
    use mod_global_parameters
    integer :: nwx,idir
    
    call rotating_frame_params_read(par_files)
    
  end subroutine rotating_frame_init
  
  !> w[iw]=w[iw]+qdt*S[wCT,qtC,x] where S is the source based on wCT within ixO
  subroutine rotating_frame_add_source(qdt,ixI^L,ixO^L,wCT,w,x)
    use mod_global_parameters
    use mod_usr_methods
    use mod_geometry
    
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, x(ixI^S,1:ndim)
    double precision, intent(in)    :: wCT(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    integer                         :: idir

    ! .. local ..

    double precision                :: rotating_terms(ixI^S),frame_omega(ixI^S)

    select case (coordinate)
    case (cylindrical)
       rotating_terms(ixO^S) =  omega_frame**2 * x(ixO^S,r_) * wCT(ixO^S,iw_rho)
       if (phi_ > 0) then
          rotating_terms(ixO^S) = rotating_terms(ixO^S) + 2.d0 * omega_frame *wCT(ixO^S,iw_mom(phi_))
       end if
       w(ixO^S, iw_mom(r_)) = w(ixO^S, iw_mom(r_)) + qdt * rotating_terms(ixO^S)
       
       if(phi_>0 .and. .not. angmomfix) then
          rotating_terms(ixO^S)   = - two * omega_frame * wCT(ixO^S,iw_mom(r_))
          w(ixO^S, iw_mom(phi_)) = w(ixO^S, iw_mom(phi_)) + qdt * rotating_terms(ixO^S)
       end if
    case (spherical)
       !call mpistop("Rotating frame not implemented yet with spherical geometries")
       ! source[mrad] = 2mphi*vframe/r + rho*vframe**2/r

       call rotating_frame_omega(x,ixI^L,ixO^L,frame_omega)
       rotating_terms(ixO^S) = frame_omega(ixO^S)**2 * x(ixO^S,r_) * wCT(ixO^S,iw_rho)
        if (phi_ > 0) then
           rotating_terms(ixO^S) = rotating_terms(ixO^S) + &
                2.d0 * frame_omega(ixO^S) * wCT(ixO^S,iw_mom(phi_))
       end if
       w(ixO^S, iw_mom(r_)) = w(ixO^S, iw_mom(r_)) + qdt * rotating_terms(ixO^S)

       {^NOONED
       ! source[mtheta] = cot(theta)*(2mphi*vframe/r + rho*vframe**2/r )
       rotating_terms(ixO^S) =  frame_omega(ixO^S)**2 * x(ixO^S,r_) * wCT(ixO^S,iw_rho)
       if (phi_>0 ) then
          rotating_terms(ixO^S) = rotating_terms(ixO^S) + &
               2.d0*wCT(ixO^S, iw_mom(phi_))* frame_omega(ixO^S)
       end if
       w(ixO^S, iw_mom(2)) = w(ixO^S, iw_mom(2)) + qdt * rotating_terms(ixO^S)/ tan(x(ixO^S, 2))

       if (phi_>0  .and. .not. angmomfix) then
          ! source[mphi]=-2*mr*vframe/r-2*cot(theta)*mtheta*vframe/r
          rotating_terms(ixO^S) = -2.d0*frame_omega(ixO^S)* wCT(ixO^S, iw_mom(r_))&
               - 2.d0*wCT(ixO^S, iw_mom(2)) * frame_omega(ixO^S)/ tan(x(ixO^S, 2))
          w(ixO^S, iw_mom(3)) = w(ixO^S, iw_mom(3)) + qdt * rotating_terms(ixO^S)
       end if
       }
    case default
       call mpistop("Rotating frame not implemented in this geometries")
    end select
    
  end subroutine rotating_frame_add_source


  subroutine rotating_frame_velocity(x,ixI^L,ixO^L,frame_vel)
        use mod_global_parameters
    use mod_usr_methods
    use mod_geometry
    
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out)   :: frame_vel(ixI^S)

    frame_vel(ixO^S) = x(ixO^S,r_) * omega_frame
    if (coordinate == spherical .and. ndim > 1) then
       frame_vel(ixO^S) = frame_vel(ixO^S) * dsin(x(ixO^S,2))
    end if

       
  end subroutine rotating_frame_velocity

    subroutine rotating_frame_omega(x,ixI^L,ixO^L,frame_omega)
        use mod_global_parameters
    use mod_usr_methods
    use mod_geometry
    
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out)   :: frame_omega(ixI^S)

    frame_omega(ixO^S) =  omega_frame
    if (coordinate == spherical .and. ndim > 1) then
       frame_omega(ixO^S) = frame_omega(ixO^S)* dsin(x(ixO^S,2))
    end if

       
  end subroutine rotating_frame_omega
end module mod_rotating_frame
