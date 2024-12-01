!> Module for including rotating frame in (magneto)hydrodynamics simulations
!> The rotation vector is assumed to be along z direction
!> (both in cylindrical and spherical)

module mod_rotating_frame
  implicit none

  !> Rotation frequency of the frame
  double precision :: omega_frame

  !> Index of the density (in the w array)
  integer, private, parameter :: rho_ = 1
  
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
  subroutine rotating_frame_add_source(qdt,dtfactor,ixI^L,ixO^L,wCT,w,x)
    use mod_global_parameters
    use mod_usr_methods
    use mod_geometry
    use mod_physics, only: phys_energy, phys_internal_e
    use mod_comm_lib, only: mpistop
    
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, dtfactor,x(ixI^S,1:ndim)
    double precision, intent(in)    :: wCT(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    integer                         :: idir

    ! .. local ..
    double precision :: rotating_terms(ixI^S), frame_omega(ixI^S)
    double precision :: work(ixI^S)

    select case (coordinate)
    case (cylindrical)

      ! S[mrad] = 2*mphi*Omegaframe + rho*r*Omegaframe**2
      rotating_terms(ixO^S) = omega_frame**2 * x(ixO^S,r_) * wCT(ixO^S,iw_rho)

      if (phi_ > 0) then
        rotating_terms(ixO^S) = rotating_terms(ixO^S) + 2.d0 * omega_frame *wCT(ixO^S,iw_mom(phi_))*wCT(ixO^S,iw_rho)
      end if

      if(local_timestep) then
        w(ixO^S, iw_mom(r_)) = w(ixO^S, iw_mom(r_)) + block%dt(ixO^S)*dtfactor * rotating_terms(ixO^S)
      else
        w(ixO^S, iw_mom(r_)) = w(ixO^S, iw_mom(r_)) + qdt * rotating_terms(ixO^S)
      endif
      ! S[mphi] = -2*mrad*Omegaframe
      if (phi_ > 0) then
        rotating_terms(ixO^S) = - 2.0d0*omega_frame * wCT(ixO^S,iw_mom(r_))*wCT(ixO^S,iw_rho)
        if(local_timestep) then
          w(ixO^S, iw_mom(phi_)) = w(ixO^S, iw_mom(phi_)) + block%dt(ixO^S)*dtfactor * rotating_terms(ixO^S)
        else
          w(ixO^S, iw_mom(phi_)) = w(ixO^S, iw_mom(phi_)) + qdt * rotating_terms(ixO^S)
        endif
      end if

      ! S[etot] = mrad*r*Omegaframe**2
      if (phys_energy .and. (.not.phys_internal_e)) then
        if(local_timestep) then
          w(ixO^S, iw_e) = w(ixO^S, iw_e) + block%dt(ixO^S) *dtfactor * omega_frame**2 * x(ixO^S,r_) * wCT(ixO^S,iw_mom(r_))*wCT(ixO^S,iw_rho)
        else
          w(ixO^S, iw_e) = w(ixO^S, iw_e) + qdt * omega_frame**2 * x(ixO^S,r_) * wCT(ixO^S,iw_mom(r_))*wCT(ixO^S,iw_rho)
        endif
      endif

    case (spherical)
       frame_omega(ixO^S) = omega_frame{^NOONED * dsin(x(ixO^S,2))}

       ! S[mrad] = 2*mphi*Omegaframe + rho*r*Omegaframe**2
       rotating_terms(ixO^S) = frame_omega(ixO^S)**2 * x(ixO^S,r_) * wCT(ixO^S,iw_rho)

       if (phi_ > 0) then
         rotating_terms(ixO^S) = rotating_terms(ixO^S) + &
                2.d0 * frame_omega(ixO^S) * wCT(ixO^S,iw_mom(phi_))*wCT(ixO^S,iw_rho)
       end if
        if(local_timestep) then
          w(ixO^S, iw_mom(r_)) = w(ixO^S, iw_mom(r_)) + block%dt(ixO^S)*dtfactor * rotating_terms(ixO^S)
        else
          w(ixO^S, iw_mom(r_)) = w(ixO^S, iw_mom(r_)) + qdt * rotating_terms(ixO^S)
        endif
       {^NOONED
       ! S[mtheta] = cot(theta) * S[mrad], reuse above rotating_terms
        if(local_timestep) then
          w(ixO^S, iw_mom(2)) = w(ixO^S, iw_mom(2)) + block%dt(ixO^S)*dtfactor * rotating_terms(ixO^S)/ tan(x(ixO^S, 2))
        else
          w(ixO^S, iw_mom(2)) = w(ixO^S, iw_mom(2)) + qdt * rotating_terms(ixO^S)/ tan(x(ixO^S, 2))
        endif
       ! S[mphi] = -2*Omegaframe * (mrad + cot(theta)*mtheta)
       if (phi_ > 0) then
         rotating_terms(ixO^S) = -2.d0*frame_omega(ixO^S)* wCT(ixO^S, iw_mom(r_))*wCT(ixO^S,iw_rho)&
               - 2.d0*wCT(ixO^S, iw_mom(2))*wCT(ixO^S,iw_rho) * frame_omega(ixO^S)/ tan(x(ixO^S, 2))
        if(local_timestep) then
          w(ixO^S, iw_mom(3)) = w(ixO^S, iw_mom(3)) + block%dt(ixO^S)*dtfactor * rotating_terms(ixO^S)
        else
         w(ixO^S, iw_mom(3)) = w(ixO^S, iw_mom(3)) + qdt * rotating_terms(ixO^S)
        endif
       end if
       }

       ! S[etot] = r*Omegaframe**2 * (mrad + cot(theta)*mtheta)
       if (phys_energy .and. (.not.phys_internal_e)) then
         work(ixO^S) = frame_omega(ixO^S)**2 * x(ixO^S,r_) * wCT(ixO^S,iw_mom(r_))*wCT(ixO^S,iw_rho)
         {^NOONED
         work(ixO^S) = work(ixO^S) + frame_omega(ixO^S)**2 * x(ixO^S,r_) * wCT(ixO^S, iw_mom(2))*wCT(ixO^S,iw_rho)/ tan(x(ixO^S, 2))
         }
          if(local_timestep) then
            w(ixO^S, iw_e) = w(ixO^S, iw_e) + block%dt(ixO^S)*dtfactor* work(ixO^S)
          else
            w(ixO^S, iw_e) = w(ixO^S, iw_e) + qdt * work(ixO^S)
        endif
       endif

    case default
       call mpistop("Rotating frame not implemented in this geometry")
    end select
    
  end subroutine rotating_frame_add_source

end module mod_rotating_frame
