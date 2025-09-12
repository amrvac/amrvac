!> Module to include Coriolis and centrifugal forces in (radiation-)(magneto)
!> hydrodynamics simulations when going into the rotating frame.
!> The rotation vector is assumed to be along z direction in all geometries
!> (Cartesian, cylindrical, and spherical).
!> 11.2019 developed for cylindrical and spherical geometry by Heloise Meheut
!>         and Clement Robert
!> 12.2022 update for fictitious work by Rony Keppens
!> 08.2025 update to Cartesian geometry by Florian Driessen
module mod_rotating_frame
  implicit none

  !> Rotation frequency of the frame
  double precision :: omega_frame

  !> Index of the density (in the w array)
  integer, private, protected :: rho_

  !> Indices of the momentum density
  integer, allocatable, private, protected :: mom(:)

  !> Index of the energy density (-1 if not present)
  integer, private, protected :: e_

  !> Indices of Cartesian coordinate directions
  integer, private, parameter :: x_ = 1, y_ = 2
  
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
    
    rho_ = iw_rho
    if (.not.allocated(mom)) allocate(mom(ndir))
    mom = iw_mom
    e_ = iw_e

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
    double precision, intent(in)    :: qdt, dtfactor, x(ixI^S,1:ndim)
    double precision, intent(in)    :: wCT(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    ! .. local ..
    double precision :: rotating_terms(ixI^S), frame_omega(ixI^S), work(ixI^S)
    double precision :: v(ixI^S,ndir)
    integer          :: idir

    ! NOTE: wCT is passed into subroutine as array of primitive variables
    do idir = 1,ndir
      v(ixI^S,idir) = wCT(ixI^S,mom(idir))
    end do

    select case (coordinate)
    case (Cartesian)

      ! S[mx] = rho*x*Omegaframe^2 + 2*my*Omegaframe
      rotating_terms(ixO^S) = omega_frame**2.0d0 * x(ixO^S,x_) * wCT(ixO^S,rho_)

      if (ndir > 1) then
        rotating_terms(ixO^S) = rotating_terms(ixO^S) &
             + 2.0d0*omega_frame * v(ixO^S,y_)*wCT(ixO^S,rho_)
      end if

      if (local_timestep) then
        w(ixO^S,mom(x_)) = w(ixO^S,mom(x_)) &
             + block%dt(ixO^S)*dtfactor * rotating_terms(ixO^S)
      else
        w(ixO^S,mom(x_)) = w(ixO^S,mom(x_)) + qdt * rotating_terms(ixO^S)
      end if

      ! S[my] = rho*y*Omegaframe^2 - 2*mx*Omegaframe
      if (ndir > 1) then
        rotating_terms(ixO^S) = -2.0d0*omega_frame * v(ixO^S,x_)*wCT(ixO^S,rho_)

        {^NOONED
        rotating_terms(ixO^S) = rotating_terms(ixO^S) &
             + omega_frame**2.0d0 * x(ixO^S,y_) * wCT(ixO^S,rho_)
        }

        if (local_timestep) then
          w(ixO^S,mom(y_)) = w(ixO^S,mom(y_)) &
               + block%dt(ixO^S)*dtfactor * rotating_terms(ixO^S)
        else
          w(ixO^S,mom(y_)) = w(ixO^S,mom(y_)) + qdt * rotating_terms(ixO^S)
        end if
      end if

      ! S[etot] = Omegaframe^2 * (x*mx + y*my)
      if (phys_energy .and. (.not.phys_internal_e)) then
        work(ixO^S) = &
             omega_frame**2.0d0 * x(ixO^S,x_) * v(ixO^S,x_)*wCT(ixO^S,rho_)

        {^NOONED
        work(ixO^S) = work(ixO^S) + omega_frame**2.0d0 * x(ixO^S,y_) &
             * v(ixO^S,y_)*wCT(ixO^S,rho_)
        }

        if (local_timestep) then
          w(ixO^S,e_) = w(ixO^S,e_) + block%dt(ixO^S)*dtfactor * work(ixO^S)
        else
          w(ixO^S,e_) = w(ixO^S,e_) + qdt * work(ixO^S)
        end if
      end if

    case (cylindrical)

      ! S[mrad] = 2*mphi*Omegaframe + rho*r*Omegaframe^2
      rotating_terms(ixO^S) = omega_frame**2.0d0 * x(ixO^S,r_) * wCT(ixO^S,rho_)

      if (phi_ > 0) then
        rotating_terms(ixO^S) = rotating_terms(ixO^S) &
             + 2.0d0*omega_frame * v(ixO^S,phi_)*wCT(ixO^S,rho_)
      end if

      if (local_timestep) then
        w(ixO^S,mom(r_)) = w(ixO^S,mom(r_)) &
             + block%dt(ixO^S)*dtfactor * rotating_terms(ixO^S)
      else
        w(ixO^S,mom(r_)) = w(ixO^S,mom(r_)) + qdt * rotating_terms(ixO^S)
      end if

      ! S[mphi] = -2*mrad*Omegaframe
      if (phi_ > 0) then
        rotating_terms(ixO^S) = -2.0d0*omega_frame * v(ixO^S,r_)*wCT(ixO^S,rho_)

        if (local_timestep) then
          w(ixO^S,mom(phi_)) = w(ixO^S,mom(phi_)) &
               + block%dt(ixO^S)*dtfactor * rotating_terms(ixO^S)
        else
          w(ixO^S,mom(phi_)) = w(ixO^S,mom(phi_)) &
               + qdt * rotating_terms(ixO^S)
        end if
      end if

      ! S[etot] = mrad*r*Omegaframe^2
      if (phys_energy .and. (.not.phys_internal_e)) then
        work(ixO^S) = &
             omega_frame**2.0d0 * x(ixO^S,r_) * v(ixO^S,r_)*wCT(ixO^S,rho_)

        if (local_timestep) then
          w(ixO^S,e_) = w(ixO^S,e_) + block%dt(ixO^S)*dtfactor * work(ixO^S)
        else
          w(ixO^S,e_) = w(ixO^S,e_) + qdt * work(ixO^S)
        end if
      end if

    case (spherical)
      frame_omega(ixO^S) = omega_frame{^NOONED * sin(x(ixO^S,2))}

      ! S[mrad] = 2*mphi*Omegaframe + rho*r*Omegaframe^2
      rotating_terms(ixO^S) = &
           frame_omega(ixO^S)**2.0d0 * x(ixO^S,r_) * wCT(ixO^S,rho_)

      if (phi_ > 0) then
        rotating_terms(ixO^S) = rotating_terms(ixO^S) &
             + 2.0d0*frame_omega(ixO^S) * v(ixO^S,phi_)*wCT(ixO^S,rho_)
      end if

      if (local_timestep) then
        w(ixO^S,mom(r_)) = w(ixO^S,mom(r_)) &
             + block%dt(ixO^S)*dtfactor * rotating_terms(ixO^S)
      else
        w(ixO^S,mom(r_)) = w(ixO^S,mom(r_)) + qdt * rotating_terms(ixO^S)
      end if

      {^NOONED
      ! S[mtheta] = cot(theta) * S[mrad], reuse above rotating_terms
      if (local_timestep) then
        w(ixO^S,mom(2)) = w(ixO^S,mom(2)) &
             + block%dt(ixO^S)*dtfactor * rotating_terms(ixO^S)/tan(x(ixO^S,2))
      else
        w(ixO^S,mom(2)) = w(ixO^S,mom(2)) &
             + qdt * rotating_terms(ixO^S)/tan(x(ixO^S,2))
      end if

      ! S[mphi] = -2*Omegaframe * (mrad + cot(theta)*mtheta)
      if (phi_ > 0) then
        rotating_terms(ixO^S) = -2.0d0*frame_omega(ixO^S) * wCT(ixO^S,rho_) &
             * ( v(ixO^S,r_) + v(ixO^S,2)/tan(x(ixO^S,2)) )

        if (local_timestep) then
          w(ixO^S,mom(phi_)) = w(ixO^S,mom(phi_)) &
               + block%dt(ixO^S)*dtfactor * rotating_terms(ixO^S)
        else
          w(ixO^S,mom(phi_)) = w(ixO^S,mom(phi_)) &
               + qdt * rotating_terms(ixO^S)
        end if
      end if
      }

      ! S[etot] = r*Omegaframe^2 * (mrad + cot(theta)*mtheta)
      if (phys_energy .and. (.not.phys_internal_e)) then
        work(ixO^S) = frame_omega(ixO^S)**2.0d0 * x(ixO^S,r_) &
             * v(ixO^S,r_)*wCT(ixO^S,rho_)

        {^NOONED
        work(ixO^S) = work(ixO^S) + frame_omega(ixO^S)**2.0d0 * x(ixO^S,r_) &
             * v(ixO^S,2)*wCT(ixO^S,rho_)/tan(x(ixO^S,2))
        }

        if (local_timestep) then
          w(ixO^S,e_) = w(ixO^S,e_) + block%dt(ixO^S)*dtfactor * work(ixO^S)
        else
          w(ixO^S,e_) = w(ixO^S,e_) + qdt * work(ixO^S)
        end if
      end if

    case default
      call mpistop("Rotating frame not implemented in this geometry")
    end select
    
  end subroutine rotating_frame_add_source

end module mod_rotating_frame
