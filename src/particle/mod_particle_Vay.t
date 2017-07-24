!> Vay particle mover (improvement of Boris method)
module mod_particle_Vay
  use mod_particle_base

  implicit none

  private

  public :: Vay_init

contains

  subroutine Vay_init()
    use mod_global_parameters
    integer :: idir, nwx

    if (physics_type/='mhd') call mpistop("Vay mover need magnetic field!")
    if (ndir/=3) call mpistop("Vay mover need ndir=3!")

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
    particles_integrate     => Vay_integrate_particles
  end subroutine Vay_init

  !> Relativistic Boris scheme
  subroutine Vay_integrate_particles(end_time)
    use mod_global_parameters
    double precision, intent(in)      :: end_time
    integer                           :: ipart, iipart
    double precision                  :: lfac, q, m, dt_p, cosphi, sinphi, phi1, phi2, r, re
    double precision, dimension(ndir) :: b, e, emom, uminus, t_geom, s, udash, tmp, uplus, xcart1, xcart2, ucart2, radmom
    double precision, dimension(ndir) :: uhalf, tau, uprime, ustar, t
    double precision                  :: lfacprime, sscal, sigma

    do iipart=1,nparticles_active_on_mype
      ipart = particles_active_on_mype(iipart);
      q     = particle(ipart)%self%q
      m     = particle(ipart)%self%m

      dt_p = Vay_get_particle_dt(particle(ipart), end_time)
      particle(ipart)%self%dt = dt_p

      ! Push particle over half time step
      call get_lfac(particle(ipart)%self%u,lfac)
      particle(ipart)%self%x(1:ndir) = particle(ipart)%self%x(1:ndir) &
           + 0.5d0 * dt_p * particle(ipart)%self%u(1:ndir)/lfac &
           * const_c / unit_length

      ! Get E, B at new position
      call get_vec(bp, particle(ipart)%igrid, &
           particle(ipart)%self%x,particle(ipart)%self%t,b)
      call get_vec(ep, particle(ipart)%igrid, &
           particle(ipart)%self%x,particle(ipart)%self%t,e)

      ! 'Kick' particle (update velocity)
      select case(typeaxial)

        ! CARTESIAN COORDINATES
      case('slab')
        ! Vay mover
        call cross(particle(ipart)%self%u,b,tmp)
        uhalf(1:ndir) = particle(ipart)%self%u(1:ndir) + &
             q * dt_p /(2.0d0 * m * const_c) * (e(1:ndir) + tmp(1:ndir)/lfac)

        tau = q * dt_p / (2.d0 * m * const_c) * b
        uprime = uhalf + q * dt_p / (2.d0 * m * const_c) * e
        call get_lfac(uprime,lfacprime)
        sigma = lfacprime**2 - dot_product(tau,tau)
        ustar = dot_product(uprime,tau) / const_c
        lfac = sqrt((sigma + sqrt(sigma**2 + 4.d0 * &
             dot_product(tau,tau) + dot_product(ustar,ustar))) / 2.d0)

        t = tau / lfac
        sscal = 1.d0 / (1.d0 + dot_product(t,t))
        call cross(uprime,t,tmp)
        particle(ipart)%self%u = sscal * (uprime + dot_product(uprime,t) * t + tmp)

      case default
        call mpistop("This geometry is not supported in mod_particle_Vay")
      end select

      call get_lfac(particle(ipart)%self%u,lfac)

      ! Push particle over half time step at end
      particle(ipart)%self%x(1:ndir) = particle(ipart)%self%x(1:ndir) &
           + 0.5d0 * dt_p * particle(ipart)%self%u(1:ndir)/lfac &
           * const_c / unit_length

      ! Time update
      particle(ipart)%self%t = particle(ipart)%self%t + dt_p

      ! Payload update
      if (npayload > 0) then
        ! current gyroradius
        call cross(particle(ipart)%self%u,b,tmp)
        tmp = tmp / sqrt(sum(b(:)**2))
        particle(ipart)%payload(1) = sqrt(sum(tmp(:)**2)) / sqrt(sum(b(:)**2)) * &
             m / abs(q) * const_c**2
      end if

      ! e.b payload
      if (npayload>1) then
        particle(ipart)%payload(2) = &
             sum(e(:)*b(:))/sqrt(sum(b(:)**2)*sum(e(:)**2))
      end if

    end do ! ipart loop

  end subroutine Vay_integrate_particles

  function Vay_get_particle_dt(partp, end_time) result(dt_p)
    use mod_global_parameters
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

    call get_vec(bp, partp%igrid,partp%self%x,partp%self%t,b)
    absb = sqrt(sum(b(:)**2))
    call get_lfac(partp%self%u,lfac)

    ! CFL timestep
    ! make sure we step only one cell at a time:
    v(:) = abs(const_c * partp%self%u(:) / lfac)

    ! convert to angular velocity:
    if(typeaxial =='cylindrical'.and.phi_>0) then
      v(phi_) = abs(v(phi_)/partp%self%x(r_))
    end if

    dt_cfl = min(bigdouble, &
         {rnode(rpdx^D_,partp%igrid)/max(v(^D), smalldouble)})

    if(typeaxial =='cylindrical'.and.phi_>0) then
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
    call limit_dt_endtime(end_time - partp%self%t, dt_p)

  end function Vay_get_particle_dt

end module mod_particle_Vay
