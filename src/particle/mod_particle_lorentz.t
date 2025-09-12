!> Particle mover with Newtonian/relativistic Boris scheme for Lorentz dynamics
!> By Jannis Teunissen, Bart Ripperda, Oliver Porth, and Fabio Bacchini (2016-2020)
module mod_particle_lorentz
  use mod_particle_base

  private

  public :: Lorentz_init
  public :: Lorentz_create_particles
  integer, parameter :: Boris=1, Vay=2, HC=3, LM=4

  ! Variables
  public :: bp, ep, vp, jp, rhop

contains

  subroutine Lorentz_init()
    use mod_global_parameters
    integer :: idir, nwx

    if(physics_type/='mhd') call mpistop("Lorentz particles need magnetic field!")
    if(ndir/=3) call mpistop("Lorentz particles need ndir=3!")

    ! The first 6 gridvars are always B and E
    ngridvars = ndir*2
    nwx=0
    ! density
    if(particles_etah>0) nwx = 1

    allocate(bp(ndir))
    do idir = 1, ndir
      nwx = nwx + 1
      bp(idir) = nwx
    end do
    allocate(ep(ndir)) ! electric field
    do idir = 1, ndir
      nwx = nwx + 1
      ep(idir) = nwx
    end do
    allocate(vp(ndir)) ! fluid velocity
    do idir = 1, ndir
      nwx = nwx + 1
      vp(idir) = nwx
    end do
    allocate(jp(ndir)) ! current
    do idir = 1, ndir
      nwx = nwx + 1
      jp(idir) = nwx
    end do 
    nwx = nwx + 1 ! density
    rhop = nwx

    ngridvars=nwx

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
    case default
      integrator = Boris
    end select

    particles_integrate => Lorentz_integrate_particles

  end subroutine Lorentz_init

  subroutine Lorentz_create_particles()

    use mod_global_parameters
    use mod_usr_methods, only: usr_create_particles, usr_update_payload, usr_check_particle

    double precision :: lfac
    double precision :: x(3, num_particles)
    double precision :: v(3, num_particles)
    double precision :: q(num_particles)
    double precision :: m(num_particles)
    double precision :: rrd(num_particles,ndir)
    double precision :: defpayload(ndefpayload)
    double precision :: usrpayload(nusrpayload)
    integer          :: n, idir, igrid, ipe_particle, nparticles_local
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

    nparticles_local = 0

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

        nparticles_local = nparticles_local + 1

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
        call Lorentz_update_payload(igrid,x(:,n),v(:,n)*lfac,q(n),m(n),defpayload,ndefpayload,0.d0)
        particle(n)%payload(1:ndefpayload) = defpayload
        if (associated(usr_update_payload)) then
          call usr_update_payload(igrid,x(:,n),v(:,n)*lfac,q(n),m(n),usrpayload,nusrpayload,0.d0)
          particle(n)%payload(ndefpayload+1:npayload) = usrpayload
        end if
      end if
    end do

    call MPI_ALLREDUCE(nparticles_local,nparticles,1,MPI_INTEGER,MPI_SUM,icomm,ierrmpi)

    ! write the first csv file of particles
    t_next_output=global_time
    call write_particle_output()
    t_next_output=t_next_output+dtsave_particles

  end subroutine Lorentz_create_particles

  subroutine lorentz_fill_gridvars
    use mod_global_parameters
    use mod_usr_methods, only: usr_particle_fields
    use mod_geometry

    double precision, dimension(ixG^T,1:ndir) :: beta
    double precision, dimension(ixG^T,1:nw)   :: w,wold
    double precision                          :: current(ixG^T,7-2*ndir:3)
    double precision, dimension(ixG^T,1:ndir) :: vE, bhat
    double precision, dimension(ixG^T)        :: kappa, kappa_B, absB, tmp
    integer                                   :: igrid, iigrid, idir, idim
    integer                                   :: idirmin

    ! Fill electromagnetic quantities
    call fill_gridvars_default()

  end subroutine lorentz_fill_gridvars

  !> Relativistic particle integrator
  subroutine Lorentz_integrate_particles(end_time)
    use mod_global_parameters
    use mod_geometry
    use mod_usr_methods, only: usr_create_particles, usr_update_payload
    double precision, intent(in)      :: end_time
    double precision                  :: defpayload(ndefpayload)
    double precision                  :: usrpayload(nusrpayload)
    double precision                  :: lfac, q, m, dt_p
    double precision                  :: xp(ndir), xpm(ndir), xpc(ndir), xpcm(ndir)
    double precision                  :: up(ndir), upc(ndir), tp
    double precision, dimension(ndir) :: b, e, bc, ec, vfluid, current
    double precision                  :: rho, rhoold, td
    integer                           :: ipart, iipart, igrid

    do iipart=1,nparticles_active_on_mype
      ipart = particles_active_on_mype(iipart);
      q     = particle(ipart)%self%q
      m     = particle(ipart)%self%m
      igrid = particle(ipart)%igrid
      xp    = particle(ipart)%self%x
      up    = particle(ipart)%self%u
      tp    = particle(ipart)%self%time   

      dt_p = Lorentz_get_particle_dt(particle(ipart), end_time)
      particle(ipart)%self%dt = dt_p

      ! Push particle over first half time step
      ! all pushes done in Cartesian coordinates
      call partcoord_to_cartesian(xp,xpc)
      call partvec_to_cartesian(xp,up,upc)
      call get_lfac(upc,lfac)
      xpcm = xpc + 0.5d0 * dt_p * upc/lfac
      call partcoord_from_cartesian(xpm,xpcm)
      ! Fix xp if the 0,2*pi boundary was crossed in cylindrical/spherical coords
      call fix_phi_crossing(xpm,igrid)
      
      ! Get E, B at n+1/2 position
      call get_bfield(igrid, xpm, tp+dt_p/2.d0, b)
      call get_efield(igrid, xpm, tp+dt_p/2.d0, e)

!      call get_vec(vp, igrid, xpm, tp+dt_p/2.d0, vfluid)
!      call get_vec(jp, igrid, xpm, tp+dt_p/2.d0, current)
!      select case (coordinate)
!      case (Cartesian,Cartesian_stretched,spherical)
!        e(1) = -vfluid(2)*b(3)+vfluid(3)*b(2) + particles_eta*current(1)
!        e(2) = vfluid(1)*b(3)-vfluid(3)*b(1) + particles_eta*current(2)
!        e(3) = -vfluid(1)*b(2)+vfluid(2)*b(1) + particles_eta*current(3)
!      case (cylindrical)
!        e(r_) = -vfluid(phi_)*b(z_)+vfluid(z_)*b(phi_) + particles_eta*current(r_)
!        e(phi_) = vfluid(r_)*b(z_)-vfluid(z_)*b(r_) + particles_eta*current(phi_)
!        e(z_) = -vfluid(r_)*b(phi_)+vfluid(phi_)*b(r_) + particles_eta*current(z_)
!      end select
!      if (particles_etah > zero) then
!        call interpolate_var(igrid,ixG^LL,ixM^LL,ps(igrid)%w(ixG^T,1),ps(igrid)%x,xpm,rho)
!        if (time_advance) then
!          td = (tp+dt_p/2.d0 - global_time) / dt
!          call interpolate_var(igrid,ixG^LL,ixM^LL,pso(igrid)%w(ixG^T,1),ps(igrid)%x,xpm,rhoold)
!          rho = rhoold * (1.d0-td) + rho * td
!        end if
!        select case (coordinate)
!        case (Cartesian,Cartesian_stretched,spherical)
!          e(1) = e(1) + particles_etah/rho * (current(2)*b(3) - current(3)*b(2))
!          e(2) = e(2) + particles_etah/rho * (-current(1)*b(3) + current(3)*b(1))
!          e(3) = e(3) + particles_etah/rho * (current(1)*b(2) - current(2)*b(1))
!        case (cylindrical)
!          e(r_) = e(r_) + particles_etah/rho * (current(phi_)*b(z_) - current(z_)*b(phi_))
!          e(phi_) = e(phi_) + particles_etah/rho * (-current(r_)*b(z_) + current(z_)*b(r_))
!          e(z_) = e(z_) + particles_etah/rho * (current(r_)*b(phi_) - current(phi_)*b(r_))
!        end select
!      end if

      ! Convert fields to Cartesian frame
      call partvec_to_cartesian(xpm,b,bc)
      call partvec_to_cartesian(xpm,e,ec)

      ! 'Kick' particle (update velocity) based on the chosen integrator
      call Lorentz_kick(upc,ec,bc,q,m,dt_p)

      ! Push particle over second half time step
      ! all pushes done in Cartesian coordinates
      call get_lfac(upc,lfac)
      xpc = xpcm + 0.5d0 * dt_p * upc/lfac
      call partcoord_from_cartesian(xp,xpc)
      ! Fix xp if the 0,2*pi boundary was crossed in cylindrical/spherical coords
      call fix_phi_crossing(xp,igrid)
      call partvec_from_cartesian(xp,up,upc)

      ! Store updated x,u
      particle(ipart)%self%x = xp
      particle(ipart)%self%u = up

      ! Time update
      tp = tp + dt_p
      particle(ipart)%self%time = tp

      ! Update payload
      call Lorentz_update_payload(igrid,xp,up,q,m,defpayload,ndefpayload,tp)
      particle(ipart)%payload(1:ndefpayload) = defpayload
      if (associated(usr_update_payload)) then
        call usr_update_payload(igrid,xp,up,q,m,usrpayload,nusrpayload,tp)
        particle(ipart)%payload(ndefpayload+1:npayload) = usrpayload
      end if

    end do ! ipart loop

  end subroutine Lorentz_integrate_particles

  !> Momentum update subroutine for full Lorentz dynamics
  subroutine Lorentz_kick(upart,e,b,q,m,dtp)
    use mod_global_parameters
    use mod_geometry
    double precision, intent(in)      :: e(ndir), b(ndir), q, m, dtp
    double precision, intent(inout)   :: upart(ndir)
    double precision                  :: lfac, cosphi, sinphi, phi1, phi2, r, re, sigma, td
    double precision, dimension(ndir) :: emom, uprime, tau, s, tmp, uplus, xcart1, xcart2, ucart2, radmom
    double precision, dimension(ndir) :: upartk, vbar, Fk, C1, C2, dupartk
    double precision                  :: abserr, tol, lfack, J11, J12, J13, J21, J22, J23, J31, J32, J33
    double precision                  :: iJ11, iJ12, iJ13, iJ21, iJ22, iJ23, iJ31, iJ32, iJ33, Det
    integer                           :: nk, nkmax

    ! Perform momentum update based on the chosen integrator
    select case(integrator)

    ! Boris integrator (works in Cartesian and cylindrical)
    case(Boris)
      ! Momentum update
      emom = q * e * dtp /(2.0d0 * m)
      uprime = upart + emom
      !!!!!! TODO: Adjust and document losses
!      if (losses) then
!        call get_lfac(particle(ipart)%self%u,lfac)
!        re = abs(q)**2 / (m * const_c**2)
!        call cross(upart,b,tmp)
!        radmom = - third * re**2 * lfac &
!             * ( sum((e(:)+tmp(:)/lfac)**2)  &
!             -  (sum(e(:)*upart(:))/lfac)**2 ) &
!             * particle(ipart)%self%u / m / const_c * dt_p
!        uprime = uprime + radmom
!      end if

      call get_lfac(uprime,lfac)
      tau = q * b * dtp / (2.0d0 * lfac * m)
      s = 2.0d0 * tau / (1.0d0+sum(tau(:)**2))

      call cross(uprime,tau,tmp)
      call cross(uprime+tmp,s,tmp)
      uplus = uprime + tmp

      upart = uplus + emom
      !!!!!! TODO: Adjust and document losses
!      if(losses) then
!        call cross(uplus,b,tmp)
!        radmom = - third * re**2 * lfac &
!             * ( sum((e(:)+tmp(:)/lfac)**2)  &
!             -  (sum(e(:)*uplus(:))/lfac)**2 ) &
!             * uplus / m / const_c * dt_p
!        upart = upart + radmom
!      end if

    ! Vay integrator
    case(Vay)
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

    ! Higuera-Cary integrator
    case(HC)
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

    ! Lapenta-Markidis integrator
    case(LM)
      ! Initialise iteration quantities
      call get_lfac(upart,lfac)
      upartk = upart
      tau(:) = b(:)

      ! START OF THE NONLINEAR CYCLE
      abserr = 1.d0
      tol=1.d-14
      nkmax=10
      nk=0
      do while(abserr > tol .and. nk < nkmax)

        nk=nk+1

        call get_lfac(upartk,lfack)
        vbar = (upart + upartk) / (lfac + lfack)
        call cross(vbar,tau,tmp)

        ! Compute residual vector
        Fk = upartk - upart - q*dtp/m * (e + tmp)

        ! Compute auxiliary coefficients
        C1 = (lfack + lfac - upartk(1:ndim) / lfack / c_norm**2 * (upartk + upart)) / (lfack + lfac)**2
        C2 = - upartk / lfack / c_norm**2 / (lfack + lfac)**2

        ! Compute Jacobian
        J11 = 1. - q*dtp/m * (C2(1) * (upartk(2) + upart(2)) * tau(3) - C2(1) * (upartk(3) + upart(3)) * tau(2))
        J12 = - q*dtp/m * (C1(2) * tau(3) - C2(2) * (upartk(3) + upart(3)) * tau(2))
        J13 = - q*dtp/m * (C2(3) * (upartk(2) + upart(2)) * tau(3) - C1(3) * tau(2))
        J21 = - q*dtp/m * (- C1(1) * tau(3) + C2(1) * (upartk(3) + upart(3)) * tau(1))
        J22 = 1. - q*dtp/m * (- C2(2) * (upartk(1) + upart(1)) * tau(3) + C2(2) * (upartk(3) + upart(3)) * tau(1))
        J23 = - q*dtp/m * (- C2(3) * (upartk(1) + upart(1)) * tau(3) + C1(3) * tau(1))
        J31 = - q*dtp/m * (C1(1) * tau(2) - C2(1) * (upartk(2) + upart(2)) * tau(1))
        J32 = - q*dtp/m * (C2(2) * (upartk(1) + upart(1)) * tau(2) - C1(2) * tau(1))
        J33 = 1. - q*dtp/m * (C2(3) * (upartk(1) + upart(1)) * tau(2) - C2(3) * (upartk(2) + upart(2)) * tau(1))

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

    end select

  end subroutine Lorentz_kick

  !> Update payload subroutine
  subroutine Lorentz_update_payload(igrid,xpart,upart,qpart,mpart,mypayload,mynpayload,particle_time)
    use mod_global_parameters
    use mod_geometry
    integer, intent(in)           :: igrid,mynpayload
    double precision, intent(in)  :: xpart(1:ndir),upart(1:ndir),qpart,mpart,particle_time
    double precision, intent(out) :: mypayload(mynpayload)
    double precision              :: b(3), e(3), tmp(3), lfac, vfluid(3), current(3), rho, rhoold, td

    call get_bfield(igrid, xpart, particle_time, b)
    call get_efield(igrid, xpart, particle_time, e)

!    call get_vec(bp, igrid, xpart,particle_time,b)
!    call get_vec(vp, igrid, xpart,particle_time,vfluid)
!    call get_vec(jp, igrid, xpart,particle_time,current)
!    select case (coordinate)
!    case (Cartesian,Cartesian_stretched,spherical)
!      e(1) = -vfluid(2)*b(3)+vfluid(3)*b(2) + particles_eta*current(1)
!      e(2) = vfluid(1)*b(3)-vfluid(3)*b(1) + particles_eta*current(2)
!      e(3) = -vfluid(1)*b(2)+vfluid(2)*b(1) + particles_eta*current(3)
!    case (cylindrical)
!      e(r_) = -vfluid(phi_)*b(z_)+vfluid(z_)*b(phi_) + particles_eta*current(r_)
!      e(phi_) = vfluid(r_)*b(z_)-vfluid(z_)*b(r_) + particles_eta*current(phi_)
!      e(z_) = -vfluid(r_)*b(phi_)+vfluid(phi_)*b(r_) + particles_eta*current(z_)
!    end select
!    if (particles_etah > zero) then
!      call interpolate_var(igrid,ixG^LL,ixM^LL,ps(igrid)%w(ixG^T,1),ps(igrid)%x,xpart,rho)
!      if (time_advance) then
!        td = (particle_time - global_time) / dt
!        call interpolate_var(igrid,ixG^LL,ixM^LL,pso(igrid)%w(ixG^T,1),ps(igrid)%x,xpart,rhoold)
!        rho = rhoold * (1.d0-td) + rho * td
!      end if
!      select case (coordinate)
!      case (Cartesian,Cartesian_stretched,spherical)
!        e(1) = e(1) + particles_etah/rho * (current(2)*b(3) - current(3)*b(2))
!        e(2) = e(2) + particles_etah/rho * (-current(1)*b(3) + current(3)*b(1))
!        e(3) = e(3) + particles_etah/rho * (current(1)*b(2) - current(2)*b(1))
!      case (cylindrical)
!        e(r_) = e(r_) + particles_etah/rho * (current(phi_)*b(z_) - current(z_)*b(phi_))
!        e(phi_) = e(phi_) + particles_etah/rho * (-current(r_)*b(z_) + current(z_)*b(r_))
!        e(z_) = e(z_) + particles_etah/rho * (current(r_)*b(phi_) - current(phi_)*b(r_))
!      end select
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
    double precision,dimension(ndir) :: b,v
    double precision                 :: lfac,absb,dt_cfl
    double precision                 :: tout
    double precision, parameter      :: cfl = 0.5d0
    integer                          :: ipart, iipart, nout

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

!    ! convert to angular velocity:
!    if(coordinate ==cylindrical.and.phi_>0) then
!      v(phi_) = abs(v(phi_)/partp%self%x(r_))
!    end if

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
