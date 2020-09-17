!> Particle mover with Newtonian/relativistic Boris scheme for Lorentz dynamics
!> By Jannis Teunissen, Bart Ripperda, Oliver Porth, and Fabio Bacchini (2016-2020)
module mod_particle_lorentz
  use mod_particle_base

  private

  public :: Lorentz_init
  public :: Lorentz_create_particles
  integer, parameter :: Boris=1, Vay=2, HC=3, LM=4

  ! Variables
  public :: bp, ep, vp, jp

contains

  subroutine Lorentz_init()
    use mod_global_parameters
    integer :: idir, nwx

    if(physics_type/='mhd') call mpistop("Lorentz particles need magnetic field!")
    if(ndir/=3) call mpistop("Lorentz particles need ndir=3!")

    ! The first 6 gridvars are always B and E
    ngridvars = ndir*2
    nwx = 0

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
    allocate(vp(ndir))
    do idir = 1, ndir
      nwx = nwx + 1
      vp(idir) = nwx
    end do
    allocate(jp(ndir))
    do idir = 1, ndir
      nwx = nwx + 1
      jp(idir) = nwx
    end do 
    ngridvars=nwx

!    particles_fill_gridvars => fill_gridvars_default
    particles_fill_gridvars => lorentz_fill_gridvars

    if (associated(particles_define_additional_gridvars)) then
      call particles_define_additional_gridvars(ngridvars)
    end if

    select case(integrator_type_particles)
    case('Boris','boris')
      integrator = Boris
    case('Vay','vay')
      integrator = Vay
    case('HC','hc','higueracary')
      integrator = HC
    case('LM','lm','lapentamarkidis')
      integrator = LM
    end select

    particles_integrate => Lorentz_integrate_particles

  end subroutine Lorentz_init

  subroutine Lorentz_create_particles()

    use mod_global_parameters
    use mod_usr_methods, only: usr_create_particles, usr_update_payload, usr_check_particle

    integer          :: n, idir, igrid, ipe_particle
    double precision :: lfac
    double precision :: x(3, num_particles)
    double precision :: v(3, num_particles)
    double precision :: q(num_particles)
    double precision :: m(num_particles)
    double precision :: rrd(num_particles,ndir)
    double precision :: defpayload(ndefpayload)
    double precision :: usrpayload(nusrpayload)
    logical          :: follow(num_particles), check

    follow = .false.
    x      = 0.0d0

    if (mype==0) then
      if (.not. associated(usr_create_particles)) then
        ! Randomly distributed
        do idir=1,ndir
          do n = 1, num_particles
            rrd(n,idir) = rng%unif_01()
          end do
        end do
        do n=1, num_particles
          {^D&x(^D,n) = xprobmin^D + rrd(n+1,^D) * (xprobmax^D - xprobmin^D)\}
        end do
      else
        call usr_create_particles(num_particles, x, v, q, m, follow)
      end if
    end if

    call MPI_BCAST(x,3*num_particles,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
    call MPI_BCAST(v,3*num_particles,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
    call MPI_BCAST(q,num_particles,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
    call MPI_BCAST(m,num_particles,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
    call MPI_BCAST(follow,num_particles,MPI_LOGICAL,0,icomm,ierrmpi)

    nparticles = num_particles

    ! Find ipe and igrid responsible for particle
    do n = 1, num_particles
      call find_particle_ipe(x(:, n),igrid,ipe_particle)

      particle(n)%igrid = igrid
      particle(n)%ipe   = ipe_particle

      if (ipe_particle == mype) then
        check = .true.

        ! Check for user-defined modifications or rejection conditions
        if (associated(usr_check_particle)) call usr_check_particle(igrid, x(:,n), v(:,n), q(n), m(n), follow(n), check)
        if (check) then
          call push_particle_into_particles_on_mype(n)
        else
          cycle
        end if

        call get_lfac_from_velocity(v(:,n), lfac)

        allocate(particle(n)%self)
        particle(n)%self%x      = x(:,n)
        particle(n)%self%u      = v(:,n) * lfac
        particle(n)%self%q      = q(n)
        particle(n)%self%m      = m(n)
        particle(n)%self%follow = follow(n)
        particle(n)%self%index  = n
        particle(n)%self%time   = global_time
        particle(n)%self%dt     = 0.0d0

        ! initialise payloads for Lorentz module
        allocate(particle(n)%payload(npayload))
        call Lorentz_update_payload(igrid,ps(igrid)%w,pso(igrid)%w,ps(igrid)%x,x(:,n),v(:,n)*lfac,q(n),m(n),defpayload,ndefpayload,0.d0)
        particle(n)%payload(1:ndefpayload) = defpayload
        if (associated(usr_update_payload)) then
          call usr_update_payload(igrid,ps(igrid)%w,pso(igrid)%w,ps(igrid)%x,x(:,n),v(:,n)*lfac,q(n),m(n),usrpayload,nusrpayload,0.d0)
          particle(n)%payload(ndefpayload+1:npayload) = usrpayload
        end if
      end if
    end do

  end subroutine Lorentz_create_particles

  subroutine lorentz_fill_gridvars
    use mod_global_parameters
    use mod_usr_methods, only: usr_particle_fields
    use mod_geometry

    integer                                   :: igrid, iigrid, idir, idim
    double precision, dimension(ixG^T,1:ndir) :: beta
    double precision, dimension(ixG^T,1:nw)   :: w,wold
    double precision                          :: current(ixG^T,7-2*ndir:3)
    integer                                   :: idirmin
    double precision, dimension(ixG^T,1:ndir) :: vE, bhat
    double precision, dimension(ixG^T)        :: kappa, kappa_B, absB, tmp

    call fill_gridvars_default()

    do iigrid=1,igridstail; igrid=igrids(iigrid);
      w(ixG^T,1:nw) = ps(igrid)%w(ixG^T,1:nw)
      call phys_to_primitive(ixG^LL,ixG^LL,w,ps(igrid)%x)
      ! fill with velocity:
      gridvars(igrid)%w(ixG^T,vp(:)) = w(ixG^T,iw_mom(:))

      if(time_advance) then
        ! Fluid velocity
        w(ixG^T,1:nw) = pso(igrid)%w(ixG^T,1:nw)
        call phys_to_primitive(ixG^LL,ixG^LL,w,ps(igrid)%x)
        gridvars(igrid)%wold(ixG^T,vp(:)) = w(ixG^T,iw_mom(:))
      end if
    end do

  end subroutine lorentz_fill_gridvars

  !> Relativistic particle integrator
  subroutine Lorentz_integrate_particles(end_time)
    use mod_global_parameters
    use mod_geometry
    use mod_usr_methods, only: usr_create_particles, usr_update_payload
    double precision, intent(in)      :: end_time
    double precision                  :: defpayload(ndefpayload)
    double precision                  :: usrpayload(nusrpayload)
    integer                           :: ipart, iipart
    double precision                  :: lfac, q, m, dt_p, cosphi, sinphi, phi1, phi2, r, re
    double precision, dimension(ndir) :: b, e, emom, uminus, t_geom, s, udash, tmp, uplus, xcart1, xcart2, ucart2, radmom, vfluid, current

    do iipart=1,nparticles_active_on_mype
      ipart = particles_active_on_mype(iipart);
      q     = particle(ipart)%self%q
      m     = particle(ipart)%self%m

      dt_p = Lorentz_get_particle_dt(particle(ipart), end_time)
      particle(ipart)%self%dt = dt_p

      ! Push particle over half time step
      call get_lfac(particle(ipart)%self%u,lfac)
      particle(ipart)%self%x(1:ndir) = particle(ipart)%self%x(1:ndir) &
           + 0.5d0 * dt_p * particle(ipart)%self%u(1:ndir)/lfac

      ! Get E, B at new position
      call get_vec(bp, particle(ipart)%igrid, &
           particle(ipart)%self%x,particle(ipart)%self%time,b)
!      if (particles_eta > 0.d0) then
!        call get_vec(ep, particle(ipart)%igrid, &
!             particle(ipart)%self%x,particle(ipart)%self%time,e)
!      else
        call get_vec(vp, particle(ipart)%igrid, &
             particle(ipart)%self%x,particle(ipart)%self%time,vfluid)
        call get_vec(jp, particle(ipart)%igrid, &
             particle(ipart)%self%x,particle(ipart)%self%time,current)
        e(1) = -vfluid(2)*b(3)+vfluid(3)*b(2) + particles_eta*current(1)
        e(2) = vfluid(1)*b(3)-vfluid(3)*b(1) + particles_eta*current(2)
        e(3) = -vfluid(1)*b(2)+vfluid(2)*b(1) + particles_eta*current(3)
!      end if

      ! 'Kick' particle (update velocity) based on the chosen integrator
      call Lorentz_kick(particle(ipart)%self%x,particle(ipart)%self%u,e,b,q,m,dt_p)

      ! Push particle over half time step at end
      call get_lfac(particle(ipart)%self%u,lfac)
      particle(ipart)%self%x(1:ndir) = particle(ipart)%self%x(1:ndir) &
           + 0.5d0 * dt_p * particle(ipart)%self%u(1:ndir)/lfac

      ! Time update
      particle(ipart)%self%time = particle(ipart)%self%time + dt_p

      ! Update payload
      call Lorentz_update_payload(particle(ipart)%igrid,ps(particle(ipart)%igrid)%w,pso(particle(ipart)%igrid)%w,ps(particle(ipart)%igrid)%x, &
           particle(ipart)%self%x,particle(ipart)%self%u,q,m,defpayload,ndefpayload,particle(ipart)%self%time)
      particle(ipart)%payload(1:ndefpayload) = defpayload
      if (associated(usr_update_payload)) then
        call usr_update_payload(particle(ipart)%igrid,ps(particle(ipart)%igrid)%w,pso(particle(ipart)%igrid)%w,ps(particle(ipart)%igrid)%x,&
             particle(ipart)%self%x,particle(ipart)%self%u,q,m,usrpayload,nusrpayload,particle(ipart)%self%time)
        particle(ipart)%payload(ndefpayload+1:npayload) = usrpayload
      end if

    end do ! ipart loop

  end subroutine Lorentz_integrate_particles

  !> Momentum update subroutine for full Lorentz dynamics
  subroutine Lorentz_kick(xpart,upart,e,b,q,m,dtp)
    use mod_global_parameters
    use mod_geometry
    double precision, intent(in)      :: e(ndir), b(ndir), q, m, dtp
    double precision, intent(inout)   :: xpart(ndim), upart(ndir)
    double precision                  :: lfac, cosphi, sinphi, phi1, phi2, r, re, sigma
    double precision, dimension(ndir) :: emom, uprime, tau, s, tmp, uplus, xcart1, xcart2, ucart2, radmom
    double precision, dimension(ndir) :: upartk, vbar, Fk, C1, C2, dupartk
    double precision                  :: abserr, tol, lfack, J11, J12, J13, J21, J22, J23, J31, J32, J33
    double precision                  :: iJ11, iJ12, iJ13, iJ21, iJ22, iJ23, iJ31, iJ32, iJ33, Det
    integer                           :: nk, nkmax

    ! Perform momentum update based on the chosen integrator
    select case(integrator)

    ! Boris integrator (works in Cartesian and cylindrical)
    case(Boris)

      select case(coordinate)

      ! CARTESIAN COORDINATES
      case(Cartesian)!,Cartesian_stretched)
        ! TODO: Adjust and document Cartesian_stretched

        ! Momentum update
        emom = q * e * dtp /(2.0d0 * m)
        !!!!!! TODO: Adjust and document losses
!        if(losses) then
!          call get_lfac(particle(ipart)%self%u,lfac)
!          re = abs(q)**2 / (m * const_c**2)
!          call cross(particle(ipart)%self%u,b,tmp)
!          radmom = - third * re**2 * lfac &
!               * ( sum((e(:)+tmp(:)/lfac)**2)  &
!               -  (sum(e(:)*particle(ipart)%self%u(:))/lfac)**2 ) &
!               * particle(ipart)%self%u / m / const_c * dt_p
!        else
          radmom = 0.0d0
!        end if

        uprime = upart + emom + radmom

        call get_lfac(uprime,lfac)
        tau = q * b * dtp / (2.0d0 * lfac * m)
        s = 2.0d0 * tau / (1.0d0+sum(tau(:)**2))

        call cross(uprime,tau,tmp)
        call cross(uprime+tmp,s,tmp)
        uplus = uprime + tmp

        !!!!!! TODO: Adjust and document losses
!        if(losses) then
!          call cross(uplus,b,tmp)
!          radmom = - third * re**2 * lfac &
!               * ( sum((e(:)+tmp(:)/lfac)**2)  &
!               -  (sum(e(:)*uplus(:))/lfac)**2 ) &
!               * uplus / m / const_c * dt_p
!        else
          radmom = 0.0d0
!        end if

        ! Update momentum
        upart = uplus + emom + radmom

      ! CYLINDRICAL COORDINATES
      case (cylindrical)

        !  Momentum update
        emom = q * e * dtp /(2.0d0 * m)

        !!!!!! TODO: Adjust and document losses
!        if(losses) then
!          call get_lfac(particle(ipart)%self%u,lfac)
!          re = abs(q)**2 / (m * const_c**2)
!          call cross(particle(ipart)%self%u,b,tmp)
!          radmom = - third * re**2 * lfac &
!               * ( sum((e(:)+tmp(:)/lfac)**2)  &
!               -  (sum(e(:)*particle(ipart)%self%u(:))/lfac)**2 ) &
!               * particle(ipart)%self%u / m / const_c * dt_p
!        else
          radmom = 0.0d0
!        end if

        uprime = upart + emom + radmom

        call get_lfac(uprime,lfac)
        tau = q * b * dtp / (2.0d0 * lfac * m)
        s = 2.0d0 * tau / (1.0d0+sum(tau(:)**2))

        call cross(uprime,tau,tmp)
        call cross(uprime+tmp,s,tmp)
        uplus = uprime + tmp

        !!!!!! TODO: Adjust and document losses
!        if(losses) then
!          call cross(uplus,b,tmp)
!          radmom = - third * re**2 * lfac &
!               * ( sum((e(:)+tmp(:)/lfac)**2)  &
!               -  (sum(e(:)*uplus(:))/lfac)**2 ) &
!               * uplus / m / const_c * dt_p
!        else
          radmom = 0.0d0
!        end if

        upart = uplus + emom + radmom
        ! Position update
        ! Get cartesian coordinates:
        phi1       = xpart(phi_)
        cosphi     = cos(phi1)
        sinphi     = sin(phi1)

        xcart1(1)  = xpart(r_) * cosphi
        xcart1(2)  = xpart(r_) * sinphi
        xcart1(3)  = xpart(z_)

        ucart2(1)   = cosphi * upart(r_) - sinphi * upart(phi_)
        ucart2(2)   = cosphi * upart(phi_) + sinphi * upart(r_)
        ucart2(3)   = upart(z_)

        ! update position
        xcart2(1:ndir) = xcart1(1:ndir) &
             + dtp * ucart2(1:ndir)/lfac

        ! back to cylindrical coordinates
        phi2 = atan2(xcart2(2),xcart2(1))
        if(phi2 .lt. 0.0d0) phi2 = 2.0d0*dpi + phi2
        r = sqrt(xcart2(1)**2 + xcart2(2)**2)
        xpart(r_)   = r
        xpart(phi_) = phi2
        xpart(z_)   = xcart2(3)

        ! Rotate the momentum to the new cooridnates
        ! rotate velocities
        cosphi     = cos(phi2-phi1)
        sinphi     = sin(phi2-phi1)

        tmp = upart
        upart(r_)   = cosphi * tmp(r_)   + sinphi * tmp(phi_)
        upart(phi_) = cosphi * tmp(phi_) - sinphi * tmp(r_)
        upart(z_)   = tmp(z_)

      case default
        call mpistop("Boris particle pusher incompatible with the chosen geometry")
      end select

    ! Vay integrator (works in Cartesian)
    case(Vay)

      select case(coordinate)

      ! CARTESIAN COORDINATES
      case(Cartesian)!,Cartesian_stretched)
        ! TODO: Adjust and document Cartesian_stretched

        call get_lfac(upart,lfac)
        emom = q * e * dtp /(2.0d0 * m)
        tau = q * b * dtp / (2.0d0 * m)

        call cross(upart,tau,tmp)
        uprime = upart + 2.d0*emom + tmp/lfac

        call get_lfac(uprime,lfac)
        sigma = lfac**2 - sum(tau(:)*tau(:))
        lfac = sqrt((sigma + sqrt(sigma**2 &
             + 4.d0 * (sum(tau(:)*tau(:)) + sum(uprime(:)*tau(:)/c_norm)**2))) / 2.d0)
        
        call cross(uprime,tau,tmp)
        upart = (uprime + sum(uprime(:)*tau(:))*tau/lfac**2 + tmp/lfac) / (1.d0+sum(tau(:)*tau(:))/lfac**2)

      case default
        call mpistop("Vay particle pusher incompatible with the chosen geometry")
      end select

    ! Higuera-Cary integrator (works in Cartesian)
    case(HC)

      select case(coordinate)

      ! CARTESIAN COORDINATES
      case(Cartesian)!,Cartesian_stretched)
        ! TODO: Adjust and document Cartesian_stretched

        call get_lfac(upart,lfac)
        emom = q * e * dtp /(2.0d0 * m)
        tau = q * b * dtp / (2.0d0 * m)
        uprime = upart + emom

        call get_lfac(uprime,lfac)
        sigma = lfac**2 - sum(tau(:)*tau(:))
        lfac = sqrt((sigma + sqrt(sigma**2 &
             + 4.d0 * (sum(tau(:)*tau(:)) + sum(uprime(:)*tau(:)/c_norm)**2))) / 2.d0)
        
        call cross(uprime,tau,tmp)
        upart = (uprime + sum(uprime(:)*tau(:))*tau/lfac**2 + tmp/lfac) / (1.d0+sum(tau(:)*tau(:))/lfac**2) &
                + emom + tmp/lfac

      case default
        call mpistop("Higuera-Cary particle pusher incompatible with the chosen geometry")
      end select

    ! Lapenta-Markidis integrator (works in Cartesian)
    case(LM)

      select case(coordinate)

      ! CARTESIAN COORDINATES
      case(Cartesian)!,Cartesian_stretched)
        ! TODO: Adjust and document Cartesian_stretched

        ! Initialise iteration quantities
        call get_lfac(upart,lfac)
        upartk = upart
  
        ! START OF THE NONLINEAR CYCLE
        abserr = 1.d0
        tol=1.d-14
        nkmax=10
        nk=0
        do while(abserr > tol .and. nk < nkmax)

          nk=nk+1

          call get_lfac(upartk,lfack)
          vbar = (upart + upartk) / (lfac + lfack)
          call cross(vbar,b,tmp)

          ! Compute residual vector
          Fk = upartk - upart - q*dtp/m * (e + tmp)

          ! Compute auxiliary coefficients
          C1 = (lfack + lfac - upartk(1:ndim) / lfack / c_norm**2 * (upartk + upart)) / (lfack + lfac)**2
          C2 = - upartk / lfack / c_norm**2 / (lfack + lfac)**2

          ! Compute Jacobian
          J11 = 1. - q*dtp/m * (C2(1) * (upartk(2) + upart(2)) * b(3) - C2(1) * (upartk(3) + upart(3)) * b(2))
          J12 = - q*dtp/m * (C1(2) * b(3) - C2(2) * (upartk(3) + upart(3)) * b(2))
          J13 = - q*dtp/m * (C2(3) * (upartk(2) + upart(2)) * b(3) - C1(3) * b(2))
          J21 = - q*dtp/m * (- C1(1) * b(3) + C2(1) * (upartk(3) + upart(3)) * b(1))
          J22 = 1. - q*dtp/m * (- C2(2) * (upartk(1) + upart(1)) * b(3) + C2(2) * (upartk(3) + upart(3)) * b(1))
          J23 = - q*dtp/m * (- C2(3) * (upartk(1) + upart(1)) * b(3) + C1(3) * b(1))
          J31 = - q*dtp/m * (C1(1) * b(2) - C2(1) * (upartk(2) + upart(2)) * b(1))
          J32 = - q*dtp/m * (C2(2) * (upartk(1) + upart(1)) * b(2) - C1(2) * b(1))
          J33 = 1. - q*dtp/m * (C2(3) * (upartk(1) + upart(1)) * b(2) - C2(3) * (upartk(2) + upart(2)) * b(1))

          ! Compute inverse Jacobian
          Det = J11*J22*J33 + J21*J32*J13 + J31*J12*J23 - J11*J32*J23 - J31*J22*J13 - J21*J12*J33
          iJ11 = (J22*J33 - J23*J32) / Det
          iJ12 = (J13*J32 - J12*J33) / Det
          iJ13 = (J12*J23 - J13*J22) / Det
          iJ21 = (J23*J31 - J21*J33) / Det
          iJ22 = (J11*J33 - J13*J31) / Det
          iJ23 = (J13*J21 - J11*J23) / Det
          iJ31 = (J21*J32 - J22*J31) / Det
          iJ32 = (J12*J31 - J11*J32) / Det
          iJ33 = (J11*J22 - J12*J21) / Det

          ! Compute new upartk = upartk - J^(-1) * F(upartk)
          dupartk(1) = - (iJ11 * Fk(1) + iJ12 * Fk(2) + iJ13 * Fk(3))
          dupartk(2) = - (iJ21 * Fk(1) + iJ22 * Fk(2) + iJ23 * Fk(3))
          dupartk(3) = - (iJ31 * Fk(1) + iJ32 * Fk(2) + iJ33 * Fk(3))

          ! Check convergence
          upartk=upartk+dupartk
          abserr=sqrt(sum(dupartk(:)*dupartk(:)))

        end do
        ! END OF THE NONLINEAR CYCLE

        ! Update velocity
        upart = upartk

      case default
        call mpistop("Lapenta-Markidis particle pusher incompatible with the chosen geometry")
      end select

    end select

  end subroutine Lorentz_kick

  !> Update payload subroutine
  subroutine Lorentz_update_payload(igrid,w,wold,xgrid,xpart,upart,qpart,mpart,mypayload,mynpayload,particle_time)
    use mod_global_parameters
    integer, intent(in)           :: igrid,mynpayload
    double precision, intent(in)  :: w(ixG^T,1:nw),wold(ixG^T,1:nw)
    double precision, intent(in)  :: xgrid(ixG^T,1:ndim),xpart(1:ndir),upart(1:ndir),qpart,mpart,particle_time
    double precision, intent(out) :: mypayload(mynpayload)
    double precision              :: b(3), e(3), tmp(3), lfac, vfluid(3), current(3)

    call get_vec(bp, igrid, xpart,particle_time,b)
!    if (particles_eta > 0.d0) then
!      call get_vec(ep, igrid, xpart,particle_time,e)
!    else
      call get_vec(vp, igrid, xpart,particle_time,vfluid)
      call get_vec(jp, igrid, xpart,particle_time,current)
      e(1) = -vfluid(2)*b(3)+vfluid(3)*b(2) + particles_eta*current(1)
      e(2) = vfluid(1)*b(3)-vfluid(3)*b(1) + particles_eta*current(2)
      e(3) = -vfluid(1)*b(2)+vfluid(2)*b(1) + particles_eta*current(3)
!    end if 

    ! Payload update
    ! Lorentz factor
    if (mynpayload > 0) then
      call get_lfac(upart,lfac)     
      mypayload(1) = lfac
    end if

    ! current gyroradius
    if (mynpayload > 1) then
      call cross(upart,b,tmp)
      tmp = tmp / sqrt(sum(b(:)**2))
      mypayload(2) = mpart/abs(qpart)*sqrt(sum(tmp(:)**2)) / sqrt(sum(b(:)**2))
    end if

    ! magnetic moment
    if (mynpayload > 2) then
      mypayload(3) = mpart*sum(tmp(:)**2)/(2.0d0*sqrt(sum(b(:)**2)))
    end if

    ! e.b
    if (mynpayload > 3) then
      mypayload(4) = sum(e(:)*b(:))
    end if

  end subroutine Lorentz_update_payload

  function Lorentz_get_particle_dt(partp, end_time) result(dt_p)
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
    v(:) = abs(partp%self%u(:) / lfac)

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
    dt_p = abs( dtheta * partp%self%m * lfac &
         / (partp%self%q * absb) )

    dt_p = min(dt_p, dt_cfl)

    ! Make sure we don't advance beyond end_time
    call limit_dt_endtime(end_time - partp%self%time, dt_p)

  end function Lorentz_get_particle_dt

end module mod_particle_lorentz
