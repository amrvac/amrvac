!> Leap-frog fully implicit particle mover (implemented by Fabio Bacchini)
module mod_particle_lfimp
  use mod_particle_base

  private

  public :: lfimp_init
  public :: lfimp_create_particles

contains

  subroutine lfimp_init()
    use mod_global_parameters
    integer :: idir, nwx

    if(physics_type/='mhd') call mpistop("lfimp particles need magnetic field!")
    if(ndir/=3) call mpistop("lfimp particles need ndir=3!")
    dtsave_particles=dtsave_particles*unit_time
    ngridvars=ndir*2
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

    integrator_velocity_factor = const_c

    particles_fill_gridvars => fill_gridvars_default
    particles_integrate     => lfimp_integrate_particles
  end subroutine lfimp_init

  subroutine lfimp_create_particles()

    use mod_global_parameters
    use mod_usr_methods, only: usr_create_particles

    integer          :: n, igrid_particle, ipe_particle
    double precision :: lfac
    double precision :: x(3, num_particles)
    double precision :: v(3, num_particles)
    double precision :: q(num_particles)
    double precision :: m(num_particles)
    logical          :: follow(num_particles)

    if (.not. associated(usr_create_particles)) then
      call mpistop("Error: no usr_create_particles method specified")
    else if (mype == 0) then
      call usr_create_particles(num_particles, x, v, q, m, follow)
    end if

    call MPI_BCAST(x,3*num_particles,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
    call MPI_BCAST(v,3*num_particles,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
    call MPI_BCAST(q,num_particles,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
    call MPI_BCAST(m,num_particles,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
    call MPI_BCAST(follow,num_particles,MPI_LOGICAL,0,icomm,ierrmpi)

    nparticles = num_particles

    ! Find ipe and igrid responsible for particle
    do n = 1, num_particles
      call find_particle_ipe(x(:, n),igrid_particle,ipe_particle)

      particle(n)%igrid = igrid_particle
      particle(n)%ipe   = ipe_particle

      if (ipe_particle == mype) then
        call push_particle_into_particles_on_mype(n)

        call get_lfac_from_velocity(v(:, n), lfac)

        allocate(particle(n)%self)
        particle(n)%self%x      = x(:, n)
        particle(n)%self%u      = v(:, n) * lfac / const_c
        particle(n)%self%q      = q(n)
        particle(n)%self%m      = m(n)
        particle(n)%self%follow = follow(n)
        particle(n)%self%index  = n
        particle(n)%self%time      = 0.0d0
        particle(n)%self%dt     = 0.0d0

        ! initialise payloads for lfimp module
        allocate(particle(n)%payload(npayload))
        particle(n)%payload(:) = 0.0d0
      end if
    end do

  end subroutine lfimp_create_particles

  !> Relativistic Fully Implicit scheme
  subroutine lfimp_integrate_particles(end_time)
    use mod_global_parameters
    double precision, intent(in)      :: end_time
    integer                           :: ipart, iipart,nk
    double precision                  :: lfac, q, m, dt_p
    double precision, dimension(ndir) :: b, e, vbar, pnp, pkp, xbar, dpkp, tmp
    double precision, dimension(ndir) :: dxb, dyb, dzb, dxe, dye, dze, C1, C2, Fk
    double precision                  :: abserrx, abserry, abserrz, tolera, J11, J12, J13, J21, J22, J23, J31, J32, J33, Det, lfack
    double precision                  :: iJ11, iJ12, iJ13, iJ21, iJ22, iJ23, iJ31, iJ32, iJ33

    do iipart=1,nparticles_active_on_mype
      ipart = particles_active_on_mype(iipart);
      q     = particle(ipart)%self%q
      m     = particle(ipart)%self%m

      dt_p = lfimp_get_particle_dt(particle(ipart), end_time)
      particle(ipart)%self%dt = dt_p

      ! Initialise pkrylov
      call get_lfac(particle(ipart)%self%u,lfac)
      pkp(1:ndir) = particle(ipart)%self%u(1:ndir)
      pnp(1:ndir) = particle(ipart)%self%u(1:ndir)

      abserrx=1.
      abserrz=1.
      abserry=1.
      tolera=1.d-14

      nk=0
      do while((abserrx>tolera .or. abserry>tolera .or. abserrz>tolera) .and. nk <10)

      ! START OF THE NONLINEAR CYCLE
      ! Push particle over half time step

      call get_lfac(pkp,lfack)

      vbar(1:ndir) = (pnp(1:ndir) + pkp(1:ndir)) / (lfac + lfack)

      xbar(1:ndir) = particle(ipart)%self%x(1:ndir) &
           + 0.5d0 * dt_p * vbar(1:ndir) &
           * const_c / unit_length

      ! Get E, B at new position
      call get_vec(bp, particle(ipart)%igrid, &
           xbar,particle(ipart)%self%time,b)
      call get_vec(ep, particle(ipart)%igrid, &
           xbar,particle(ipart)%self%time,e)

      ! Compute residual vector
      Fk(1) = pkp(1) - pnp(1) - q * dt_p /(m * const_c) * (e(1) + vbar(2) * b(3) - vbar(3) * b(2))
      Fk(2) = pkp(2) - pnp(2) - q * dt_p /(m * const_c) * (e(2) - vbar(1) * b(3) + vbar(3) * b(1))
      Fk(3) = pkp(3) - pnp(3) - q * dt_p /(m * const_c) * (e(3) + vbar(1) * b(2) - vbar(2) * b(1))

      ! Get specially derived E, B at new position
      call get_dervec(ep, particle(ipart)%igrid, &
           xbar,particle(ipart)%self%time,dxe,1)
      call get_dervec(ep, particle(ipart)%igrid, &
           xbar,particle(ipart)%self%time,dye,2)
      call get_dervec(ep, particle(ipart)%igrid, &
           xbar,particle(ipart)%self%time,dze,3)
      call get_dervec(bp, particle(ipart)%igrid, &
           xbar,particle(ipart)%self%time,dxb,1)
      call get_dervec(bp, particle(ipart)%igrid, &
           xbar,particle(ipart)%self%time,dyb,2)
      call get_dervec(bp, particle(ipart)%igrid, &
           xbar,particle(ipart)%self%time,dzb,3)

      ! Compute auxiliary coefficients
      C1(1:ndim) = (lfack + lfac - pkp(1:ndim) / lfack * (pkp(1:ndim) + pnp(1:ndim))) / (lfack + lfac)**2
      C2(1:ndim) = - pkp(1:ndim) / lfack / (lfack + lfac)**2

      ! Compute Jacobian
      J11 = 1. - q * dt_p**2 /(2 * m * const_c) &
                 * (C1(1) * dxe(1) + (pkp(2) + pnp(2)) * C2(1) * dye(1) + (pkp(3) + pnp(3)) * C2(1) * dze(1) &
                 + vbar(2) * (C1(1) * dxb(3) + (pkp(2) + pnp(2)) * C2(1) * dyb(3) + (pkp(3) + pnp(3)) * C2(1) * dzb(3)) &
                 - vbar(3) * (C1(1) * dxb(2) + (pkp(2) + pnp(2)) * C2(1) * dyb(2) + (pkp(3) + pnp(3)) * C2(1) * dzb(2))) &
               - q * dt_p /(m * const_c) * (C2(1) * (pkp(2) + pnp(2)) * b(3) - C2(1) * (pkp(3) + pnp(3)) * b(2))

      J12 = - q * dt_p**2 /(2 * m * const_c) &
              * ((pkp(1) + pnp(1)) * C2(2) * dxe(1) + C1(2) * dye(1) + (pkp(3) + pnp(3)) * C2(2) * dze(1) &
              + vbar(2) * ((pkp(1) + pnp(1)) * C2(2) * dxb(3) + C1(2) * dyb(3) + (pkp(3) + pnp(3)) * C2(2) * dzb(3)) &
              - vbar(3) * ((pkp(1) + pnp(1)) * C2(2) * dxb(2) + C1(2) * dyb(2) + (pkp(3) + pnp(3)) * C2(2) * dzb(2))) &
            - q * dt_p /(m * const_c) * (C1(2) * b(3) - C2(2) * (pkp(3) + pnp(3)) * b(2))

      J13 = - q * dt_p**2 /(2 * m * const_c) &
              * ((pkp(1) + pnp(1)) * C2(3) * dxe(1) + (pkp(2) + pnp(2)) * C2(3) * dye(1) + C1(3) * dze(1) &
              + vbar(2) * ((pkp(1) + pnp(1)) * C2(3) * dxb(3) + (pkp(2) + pnp(2)) * C2(3) * dyb(3) + C1(3) * dzb(3)) &
              - vbar(3) * ((pkp(1) + pnp(1)) * C2(3) * dxb(2) + (pkp(2) + pnp(2)) * C2(3) * dyb(2) + C1(3) * dzb(2))) &
            - q * dt_p /(m * const_c) * (C2(3) * (pkp(2) + pnp(2)) * b(3) - C1(3) * b(2))

      J21 = - q * dt_p**2 /(2 * m * const_c) &
              * (C1(1) * dxe(2) + (pkp(2) + pnp(2)) * C2(1) * dye(2) + (pkp(3) + pnp(3)) * C2(1) * dze(2) &
              - vbar(1) * (C1(1) * dxb(3) + (pkp(2) + pnp(2)) * C2(1) * dyb(3) + (pkp(3) + pnp(3)) * C2(1) * dzb(3)) &
              + vbar(3) * (C1(1) * dxb(1) + (pkp(2) + pnp(2)) * C2(1) * dyb(1) + (pkp(3) + pnp(3)) * C2(1) * dzb(1))) &
            - q * dt_p /(m * const_c) * (- C1(1) * b(3) + C2(1) * (pkp(3) + pnp(3)) * b(1))

      J22 = 1. - q * dt_p**2 /(2 * m * const_c) &
                 * ((pkp(1) + pnp(1)) * C2(2) * dxe(2) + C1(2) * dye(2) + (pkp(3) + pnp(3)) * C2(2) * dze(2) &
                 - vbar(1) * ((pkp(1) + pnp(1)) * C2(2) * dxb(3) + C1(2) * dyb(3) + (pkp(3) + pnp(3)) * C2(2) * dzb(3)) &
                 + vbar(3) * ((pkp(1) + pnp(1)) * C2(2) * dxb(1) + C1(2) * dyb(1) + (pkp(3) + pnp(3)) * C2(2) * dzb(1))) &
               - q * dt_p /(m * const_c) * (- C2(2) * (pkp(1) + pnp(1)) * b(3) + C2(2) * (pkp(3) + pnp(3)) * b(1))

      J23 = - q * dt_p**2 /(2 * m * const_c) &
              * ((pkp(1) + pnp(1)) * C2(3) * dxe(2) + (pkp(2) + pnp(2)) * C2(3) * dye(2) + C1(3) * dze(2) &
              - vbar(1) * ((pkp(1) + pnp(1)) * C2(3) * dxb(3) + (pkp(2) + pnp(2)) * C2(3) * dyb(3) + C1(3) * dzb(3)) &
              + vbar(3) * ((pkp(1) + pnp(1)) * C2(3) * dxb(1) + (pkp(2) + pnp(2)) * C2(3) * dyb(1) + C1(3) * dzb(1))) &
            - q * dt_p /(m * const_c) * (- C2(3) * (pkp(1) + pnp(1)) * b(3) + C1(3) * b(1))

      J31 = - q * dt_p**2 /(2 * m * const_c) &
              * (C1(1) * dxe(3) + (pkp(2) + pnp(2)) * C2(1) * dye(3) + (pkp(3) + pnp(3)) * C2(1) * dze(3) &
              + vbar(1) * (C1(1) * dxb(2) + (pkp(2) + pnp(2)) * C2(1) * dyb(2) + (pkp(3) + pnp(3)) * C2(1) * dzb(2)) &
              - vbar(2) * (C1(1) * dxb(1) + (pkp(2) + pnp(2)) * C2(1) * dyb(1) + (pkp(3) + pnp(3)) * C2(1) * dzb(1))) &
            - q * dt_p /(m * const_c) * (C1(1) * b(2) - C2(1) * (pkp(2) + pnp(2)) * b(1))

      J32 = - q * dt_p**2 /(2 * m * const_c) &
              * ((pkp(1) + pnp(1)) * C2(2) * dxe(3) + C1(2) * dye(3) + (pkp(3) + pnp(3)) * C2(2) * dze(3) &
              + vbar(1) * ((pkp(1) + pnp(1)) * C2(2) * dxb(2) + C1(2) * dyb(2) + (pkp(3) + pnp(3)) * C2(2) * dzb(2)) &
              - vbar(2) * ((pkp(1) + pnp(1)) * C2(2) * dxb(1) + C1(2) * dyb(1) + (pkp(3) + pnp(3)) * C2(2) * dzb(1))) &
            - q * dt_p /(m * const_c) * (C2(2) * (pkp(1) + pnp(1)) * b(2) - C1(2) * b(1))

      J33 = 1. - q * dt_p**2 /(2 * m * const_c) &
                 * ((pkp(1) + pnp(1)) * C2(3) * dxe(3) + (pkp(2) + pnp(2)) * C2(3) * dye(3) + C1(3) * dze(3) &
                 + vbar(1) * ((pkp(1) + pnp(1)) * C2(3) * dxb(2) + (pkp(2) + pnp(2)) * C2(3) * dyb(2) + C1(3) * dzb(2)) &
                 - vbar(2) * ((pkp(1) + pnp(1)) * C2(3) * dxb(1) + (pkp(2) + pnp(2)) * C2(3) * dyb(1) + C1(3) * dzb(1))) &
               - q * dt_p /(m * const_c) * (C2(3) * (pkp(1) + pnp(1)) * b(2) - C2(3) * (pkp(2) + pnp(2)) * b(1))

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

      ! Compute new pkrylov pkp = pkp - J^(-1) * F(pkp)
      dpkp(1) = - (iJ11 * Fk(1) + iJ12 * Fk(2) + iJ13 * Fk(3))
      dpkp(2) = - (iJ21 * Fk(1) + iJ22 * Fk(2) + iJ23 * Fk(3))
      dpkp(3) = - (iJ31 * Fk(1) + iJ32 * Fk(2) + iJ33 * Fk(3))

      pkp=pkp+dpkp

      abserrx=abs(dpkp(1))
      abserry=abs(dpkp(2))
      abserrz=abs(dpkp(3))

      nk=nk+1
      !write(*,*) "Newton iteration #", nk, "Residuals:", abserrx, abserry, abserrz
      end do
      ! END OF THE NONLINEAR CYCLE


      ! Recompute final vbar
      call get_lfac(pkp, lfack)

      vbar(1:ndir) = (pnp(1:ndir) + pkp(1:ndir)) / (lfac + lfack)

      ! Push particles to n+1 position
      particle(ipart)%self%x(1:ndir) = particle(ipart)%self%x(1:ndir) &
           + dt_p * vbar(1:ndir) &
           * const_c / unit_length

      ! Update velocity
      particle(ipart)%self%u = pkp

      call get_lfac(particle(ipart)%self%u,lfac)

      ! Time update
      particle(ipart)%self%time = particle(ipart)%self%time + dt_p
      !write(*,*) particle(ipart)%self%time

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
             sqrt(sum(tmp(:)**2))**2 / (2.0d0*sqrt(sum(b(:)**2))) * &
             m * const_c**2
      end if

      ! e.b payload
      if (npayload>2) then
        particle(ipart)%payload(3) = &
             sum(e(:)*b(:))/sqrt(sum(b(:)**2)*sum(e(:)**2))
      end if

    end do ! ipart loop

  end subroutine lfimp_integrate_particles

  function lfimp_get_particle_dt(partp, end_time) result(dt_p)
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
    else
       call mpistop("Fully implicit particle mover requires constant time step")
    end if

    ! Make sure we don't advance beyond end_time
    call limit_dt_endtime(end_time - partp%self%time, dt_p)

  end function lfimp_get_particle_dt

  ! Get the specially interpolated field in the grid at postion x.
  ! For ideal SRMHD, we first interpolate b and u=lfac*v/c
  ! The electric field then follows from e = b x beta, where beta=u/lfac.
  ! This ensures for the resulting e that e<b and e.b=0. Interpolating on u
  ! avoids interpolation-errors leading to v>c.
  ! For (non-ideal) MHD, we directly interpolate the electric field as
  ! there is no such constraint.
  subroutine get_dervec(ix,igrid,x,tloc,vec,dir)
    use mod_global_parameters

    integer,intent(in)                                 :: ix(ndir) !< Indices in gridvars
    integer,intent(in)                                 :: igrid
    double precision,dimension(ndir), intent(in)       :: x
    double precision, intent(in)                       :: tloc
    double precision,dimension(ndir), intent(out)      :: vec
    double precision,dimension(ndir)                   :: vec1, vec2
    double precision                                   :: td
    integer                                            :: ic^D,idir
    integer,intent(in)                                 :: dir

    if (.not.time_advance) then
      do idir=1,ndir
        call derinterpolate_var(igrid,ixG^LL,ixM^LL,gridvars(igrid)%w(ixG^T,ix(idir)), &
             ps(igrid)%x(ixG^T,1:ndim),x,vec(idir),dir)
      end do
    else
      do idir=1,ndir
        call derinterpolate_var(igrid,ixG^LL,ixM^LL,gridvars(igrid)%wold(ixG^T,ix(idir)), &
             ps(igrid)%x(ixG^T,1:ndim),x,vec1(idir),dir)
        call derinterpolate_var(igrid,ixG^LL,ixM^LL,gridvars(igrid)%w(ixG^T,ix(idir)), &
             ps(igrid)%x(ixG^T,1:ndim),x,vec2(idir),dir)
      end do
      td = (tloc/unit_time - global_time) / dt
      vec(:) = vec1(:) * (1.0d0 - td) + vec2(:) * td
    end if

  end subroutine get_dervec

  subroutine derinterpolate_var(igrid,ixI^L,ixO^L,gf,x,xloc,gfloc,dir)
    use mod_global_parameters

    integer, intent(in)                   :: igrid,ixI^L, ixO^L
    double precision, intent(in)          :: gf(ixI^S)
    double precision, intent(in)          :: x(ixI^S,1:ndim)
    double precision, intent(in)          :: xloc(1:ndir)
    double precision, intent(out)         :: gfloc
    double precision                      :: xd^D
    integer, intent(in)                   :: dir
    {^IFTWOD
    double precision                      :: c00, c10
    }
    {^IFTHREED
    double precision                      :: c0, c1, c00, c10, c01, c11
    }
    integer                               :: ic^D, ic1^D, ic2^D, idir

    ! flat interpolation:
    {ic^D = int((xloc(^D)-rnode(rpxmin^D_,igrid))/rnode(rpdx^D_,igrid)) + 1 + nghostcells \}
    !gfloc = gf(ic^D)

    ! linear interpolation:
    {
    if (x({ic^DD},^D) .lt. xloc(^D)) then
      ic1^D = ic^D
    else
      ic1^D = ic^D -1
    end if
    ic2^D = ic1^D + 1
    \}

    {^D&
    if(ic1^D.lt.ixGlo^D+1 .or. ic2^D.gt.ixGhi^D-1) then
      print *, 'direction: ',^D
      print *, 'position: ',xloc(1:ndim)
      print *, 'indices:', ic1^D,ic2^D
      call mpistop('Trying to interpolate from out of grid!')
    end if
    \}

    {^IFONED
    if(dir .eq. 1) then
      xd1 = sign(1.0 / (x(ic21,1) - x(ic11,1)),xloc(1)-x(ic11,1))
      gfloc  = gf(ic11) * (- xd1) + gf(ic21) * xd1
    else
      gfloc = 0.0
    end if
    }
    {^IFTWOD
    if(dir .eq. 1) then
      xd1 = sign(1.0 / (x(ic21,ic12,1) - x(ic11,ic12,1)), xloc(1)-x(ic11,ic12,1))
      xd2 = (xloc(2)-x(ic11,ic12,2)) / (x(ic11,ic22,2) - x(ic11,ic12,2))
      c00 = gf(ic11,ic12) * (- xd1) + gf(ic21,ic12) * xd1
      c10 = gf(ic11,ic22) * (1.0d0 - xd1) + gf(ic21,ic22) * xd1
      gfloc  = c00 * (1.0d0 - xd2) + c10 * xd2
    else if( dir .eq. 2) then
      xd1 = (xloc(1)-x(ic11,ic12,1)) / (x(ic21,ic12,1) - x(ic11,ic12,1))
      xd2 = sign(1.0 / (x(ic11,ic22,2) - x(ic11,ic12,2)), xloc(2)-x(ic11,ic12,2))
      c00 = gf(ic11,ic12) * (1.0d0 - xd1) + gf(ic21,ic12) * xd1
      c10 = gf(ic11,ic22) * (1.0d0 - xd1) + gf(ic21,ic22) * xd1
      gfloc  = c00 * (- xd2) + c10 * xd2
    else
      gfloc = 0.0
    end if
    }
    {^IFTHREED
    if(dir .eq. 1) then
      xd1 = sign(1.0 / (x(ic21,ic12,ic13,1) - x(ic11,ic12,ic13,1)), xloc(1)-x(ic11,ic12,ic13,1))
      xd2 = (xloc(2)-x(ic11,ic12,ic13,2)) / (x(ic11,ic22,ic13,2) - x(ic11,ic12,ic13,2))
      xd3 = (xloc(3)-x(ic11,ic12,ic13,3)) / (x(ic11,ic12,ic23,3) - x(ic11,ic12,ic13,3))
      c00 = gf(ic11,ic12,ic13) * (- xd1) + gf(ic21,ic12,ic13) * xd1
      c10 = gf(ic11,ic22,ic13) * (- xd1) + gf(ic21,ic22,ic13) * xd1
      c01 = gf(ic11,ic12,ic23) * (- xd1) + gf(ic21,ic12,ic23) * xd1
      c11 = gf(ic11,ic22,ic23) * (- xd1) + gf(ic21,ic22,ic23) * xd1

      c0  = c00 * (1.0d0 - xd2) + c10 * xd2
      c1  = c01 * (1.0d0 - xd2) + c11 * xd2

      gfloc = c0 * (1.0d0 - xd3) + c1 * xd3
    else if(dir .eq. 2) then
      xd1 = (xloc(1)-x(ic11,ic12,ic13,1)) / (x(ic21,ic12,ic13,1) - x(ic11,ic12,ic13,1))
      xd2 = sign(1.0 / (x(ic11,ic22,ic13,2) - x(ic11,ic12,ic13,2)), xloc(2)-x(ic11,ic12,ic13,2))
      xd3 = (xloc(3)-x(ic11,ic12,ic13,3)) / (x(ic11,ic12,ic23,3) - x(ic11,ic12,ic13,3))
      c00 = gf(ic11,ic12,ic13) * (1.0d0 - xd1) + gf(ic21,ic12,ic13) * xd1
      c10 = gf(ic11,ic22,ic13) * (1.0d0 - xd1) + gf(ic21,ic22,ic13) * xd1
      c01 = gf(ic11,ic12,ic23) * (1.0d0 - xd1) + gf(ic21,ic12,ic23) * xd1
      c11 = gf(ic11,ic22,ic23) * (1.0d0 - xd1) + gf(ic21,ic22,ic23) * xd1

      c0  = c00 * (- xd2) + c10 * xd2
      c1  = c01 * (- xd2) + c11 * xd2

      gfloc = c0 * (1.0d0 - xd3) + c1 * xd3
    else
      xd1 = (xloc(1)-x(ic11,ic12,ic13,1)) / (x(ic21,ic12,ic13,1) - x(ic11,ic12,ic13,1))
      xd2 = (xloc(2)-x(ic11,ic12,ic13,2)) / (x(ic11,ic22,ic13,2) - x(ic11,ic12,ic13,2))
      xd3 = sign(1.0 / (x(ic11,ic12,ic23,3) - x(ic11,ic12,ic13,3)), xloc(3)-x(ic11,ic12,ic13,3))

      c00 = gf(ic11,ic12,ic13) * (1.0d0 - xd1) + gf(ic21,ic12,ic13) * xd1
      c10 = gf(ic11,ic22,ic13) * (1.0d0 - xd1) + gf(ic21,ic22,ic13) * xd1
      c01 = gf(ic11,ic12,ic23) * (1.0d0 - xd1) + gf(ic21,ic12,ic23) * xd1
      c11 = gf(ic11,ic22,ic23) * (1.0d0 - xd1) + gf(ic21,ic22,ic23) * xd1

      c0  = c00 * (1.0d0 - xd2) + c10 * xd2
      c1  = c01 * (1.0d0 - xd2) + c11 * xd2

      gfloc = c0 * (- xd3) + c1 * xd3
    end if
    }

  end subroutine derinterpolate_var

end module mod_particle_lfimp
