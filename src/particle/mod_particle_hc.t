!> Higuera-Cary particle mover (improvement of Boris method)
module mod_particle_hc
  use mod_particle_base

  implicit none

  private

  public :: HC_init

contains

  subroutine HC_init()
    use mod_global_parameters
    integer :: idir, nwx

    if (physics_type/='mhd') call mpistop("HC mover need magnetic field!")
    if (ndir/=3) call mpistop("HC mover need ndir=3!")

    dtsave_particles = dtsave_particles*unit_time
    ngridvars        = ndir*2
    nwx              = 0

    allocate(bp(ndir))
    do idir = 1, ndir
      nwx = nwx + 1
      bp(idir) = nwx
    end do

    allocate(ep(ndir))
    do idir = 1, ndir
      nwx = nwx + 1
      ep(idir) = nwx
    end do

    integrator_velocity_factor = const_c

    particles_fill_gridvars => fill_gridvars_default
    particles_integrate     => HC_integrate_particles
  end subroutine HC_init

  !> Relativistic Boris scheme
  subroutine HC_integrate_particles(end_time)
    use mod_global_parameters
    use mod_geometry
    double precision, intent(in)      :: end_time
    integer                           :: ipart, iipart
    double precision                  :: lfac, q, m, dt_p, cosphi, sinphi, phi1, phi2, r, re
    double precision, dimension(ndir) :: b, e, emom, uminus, t_geom, s, udash, tmp, uplus, xcart1, xcart2, ucart2, radmom
    double precision, dimension(ndir) :: uhalf, tau, uprime, ustar, time
    double precision                  :: lfacprime, sscal, sigma,lfachalf
    do iipart=1,nparticles_active_on_mype
      ipart = particles_active_on_mype(iipart);
      q     = particle(ipart)%self%q
      m     = particle(ipart)%self%m

      dt_p = HC_get_particle_dt(particle(ipart), end_time)
      particle(ipart)%self%dt = dt_p

      ! Push particle over half time step
      call get_lfac(particle(ipart)%self%u,lfac)
      particle(ipart)%self%x(1:ndir) = particle(ipart)%self%x(1:ndir) &
           + 0.5d0 * dt_p * particle(ipart)%self%u(1:ndir)/lfac &
           * const_c / unit_length

      ! Get E, B at new position
      call get_vec(bp, particle(ipart)%igrid, &
           particle(ipart)%self%x,particle(ipart)%self%time,b)
      call get_vec(ep, particle(ipart)%igrid, &
           particle(ipart)%self%x,particle(ipart)%self%time,e)

      ! 'Kick' particle (update velocity)
      select case(coordinate)

        ! CARTESIAN COORDINATES
      case(Cartesian,Cartesian_stretched)
        ! HC mover
        uhalf(1:ndir) = particle(ipart)%self%u(1:ndir) + &
             q * dt_p /(2.0d0 * m * const_c) * e(1:ndir)

        tau(1:ndir) = q * dt_p / (2.d0 * m * const_c) * b(1:ndir)
        call get_lfac(uhalf,lfachalf)
        sigma = lfachalf**2 - sum(tau(:)*tau(:))
        ustar = sum(uhalf(:)*tau(:))
        lfacprime = sqrt((sigma + sqrt(sigma**2 &
             + 4.d0 * (sum(tau(:)*tau(:)) + sum(ustar(:)*ustar(:))))) / 2.d0)

        time(1:ndir) = tau(1:ndir) / lfacprime
        sscal = 1.d0 / (1.d0 + sum(time(:)*time(:)))
        call cross(uhalf,time,tmp)
        uprime(1:ndir) = sscal * (uhalf(1:ndir) + sum(uhalf(:)*time(:)) * time(1:ndir) + tmp(1:ndir))
        call cross(uprime,time,tmp)
        particle(ipart)%self%u(1:ndir) = uprime(1:ndir) &
             + q * dt_p /(2.0d0 * m * const_c) * e(1:ndir) + tmp(1:ndir)
      case default
        call mpistop("This coordinate is not supported in mod_particle_HC")
      end select

      call get_lfac(particle(ipart)%self%u,lfac)

      ! Push particle over half time step at end
      particle(ipart)%self%x(1:ndir) = particle(ipart)%self%x(1:ndir) &
           + 0.5d0 * dt_p * particle(ipart)%self%u(1:ndir)/lfac &
           * const_c / unit_length

      ! Time update
      particle(ipart)%self%time = particle(ipart)%self%time + dt_p

      ! Payload update
      if (npayload > 0) then
        ! current gyroradius
        call get_vec(bp, particle(ipart)%igrid, &
           particle(ipart)%self%x,particle(ipart)%self%time,b)
        call cross(particle(ipart)%self%u,b,tmp)
        tmp = tmp / sqrt(sum(b(:)**2))
        particle(ipart)%payload(1) = sqrt(sum(tmp(:)**2)) / sqrt(sum(b(:)**2)) * &
             m / abs(q) * const_c**2
      end if

      ! magnetic moment
      if (npayload>1) then
        particle(ipart)%payload(2) = &
             sqrt(sum(tmp(:)**2))**2/(2.0d0*sqrt(sum(b(:)**2)))*&
             m * const_c**2
      end if

      ! e.b payload
      if (npayload>2) then
        particle(ipart)%payload(3) = &
             sum(e(:)*b(:))/sqrt(sum(b(:)**2)*sum(e(:)**2))
      end if

    end do ! ipart loop

  end subroutine HC_integrate_particles

  function HC_get_particle_dt(partp, end_time) result(dt_p)
    use mod_global_parameters
    use mod_geometry
    type(particle_ptr), intent(in)   :: partp
    double precision, intent(in)     :: end_time
    double precision                 :: dt_p
    integer                          :: ipart, iipart, nout
    double precision,dimension(ndir) :: b,v
    double precision                 :: lfac,absb,dt_cfl
    double precision                 :: tout
    double precision, parameter      :: cfl = 0.5d0

    if (const_dt_particles > 0) then
      dt_p = const_dt_particles
      return
    end if

    call get_vec(bp, partp%igrid,partp%self%x,partp%self%time,b)
    absb = sqrt(sum(b(:)**2))
    call get_lfac(partp%self%u,lfac)

    ! CFL timestep
    ! make sure we step only one cell at a time:
    v(:) = abs(const_c * partp%self%u(:) / lfac)

    ! convert to angular velocity:
    if(coordinate ==cylindrical.and.phi_>0) then
      v(phi_) = abs(v(phi_)/partp%self%x(r_))
    end if

    dt_cfl = min(bigdouble, &
         {rnode(rpdx^D_,partp%igrid)/max(v(^D), smalldouble)})

    if(coordinate ==cylindrical.and.phi_>0) then
      ! phi-momentum leads to radial velocity:
      if(phi_ .gt. ndim) dt_cfl = min(dt_cfl, &
           sqrt(rnode(rpdx1_,partp%igrid)/partp%self%x(r_)) &
           / v(phi_))
      ! limit the delta phi of the orbit (just for aesthetic reasons):
      dt_cfl = min(dt_cfl,0.1d0/max(v(phi_), smalldouble))
      ! take some care at the axis:
      dt_cfl = min(dt_cfl,(partp%self%x(r_)+smalldouble)/max(v(r_), smalldouble))
    end if

    dt_cfl = dt_cfl * cfl

    ! bound by gyro-rotation:
    dt_p = abs( dtheta * const_c/unit_length * partp%self%m * lfac &
         / (partp%self%q * absb) )

    dt_p = min(dt_p, dt_cfl)*unit_length

    ! Make sure we don't advance beyond end_time
    call limit_dt_endtime(end_time - partp%self%time, dt_p)

  end function HC_get_particle_dt

end module mod_particle_hc
