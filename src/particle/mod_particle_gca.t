module mod_particle_gca
  use mod_particle_base

  private

  !> Variable index for magnetic field
  integer, protected, allocatable      :: bp(:)
  !> Variable index for electric field
  integer, protected, allocatable      :: ep(:)

  !> Variable index for gradient B (but with relativistic correction kappa,
  !> which is given by kappa = sqrt(1 - E_perp^2/B^2)
  integer, protected, allocatable      :: grad_kappa_B(:)
  !> Variable index for B . grad(B) (curvature B drift)
  integer, protected, allocatable      :: b_dot_grad_b(:)

  ! ExB related drifts (ue = ExB/B^2)
  !> Variable index for curvature drift
  integer, protected, allocatable      :: ue_dot_grad_b(:)
  !> Variable index for polarization drift
  integer, protected, allocatable      :: b_dot_grad_ue(:)
  !> Variable index for polarization drift
  integer, protected, allocatable      :: ue_dot_grad_ue(:)

  public :: gca_init
  public :: gca_create_particles

  ! Variables
  public :: bp, ep, grad_kappa_B, b_dot_grad_b
  public :: ue_dot_grad_b, b_dot_grad_ue, ue_dot_grad_ue

contains

  subroutine gca_init()
    use mod_global_parameters
    integer :: idir, nwx

    if(physics_type/='mhd') call mpistop("GCA particles need magnetic field!")
    if(ndir/=3) call mpistop("GCA particles need ndir=3!")
    dtsave_particles=dtsave_particles*unit_time
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
    allocate(grad_kappa_B(ndir))
    do idir = 1, ndir
      nwx = nwx + 1
      grad_kappa_B(idir) = nwx
    end do
    allocate(b_dot_grad_b(ndir))
    do idir = 1, ndir
      nwx = nwx + 1
      b_dot_grad_b(idir) = nwx
    end do
    allocate(ue_dot_grad_b(ndir))
    do idir = 1, ndir
      nwx = nwx + 1
      ue_dot_grad_b(idir) = nwx
    end do
    allocate(b_dot_grad_ue(ndir))
    do idir = 1, ndir
      nwx = nwx + 1
      b_dot_grad_ue(idir) = nwx
    end do
    allocate(ue_dot_grad_ue(ndir))
    do idir = 1, ndir
      nwx = nwx + 1
      ue_dot_grad_ue(idir) = nwx
    end do
    ngridvars=nwx

    particles_fill_gridvars => gca_fill_gridvars
    particles_integrate     => gca_integrate_particles
  end subroutine gca_init

  subroutine gca_create_particles()
    ! initialise the particles
    use mod_global_parameters
    use mod_usr_methods, only: usr_create_particles

    double precision :: B(ndir), u(ndir), magmom
    double precision :: Bnorm, lfac, vnorm, vperp, vpar
    integer          :: igrid_particle, ipe_particle
    integer          :: n, idir
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

    ! first find ipe and igrid responsible for particle
    do n = 1, num_particles
      call find_particle_ipe(x(:, n),igrid_particle,ipe_particle)

      particle(n)%igrid = igrid_particle
      particle(n)%ipe   = ipe_particle

      if(ipe_particle == mype) then
        call push_particle_into_particles_on_mype(n)
        call get_lfac_from_velocity(v(:, n), lfac)

        allocate(particle(n)%self)
        particle(n)%self%x      = x(:, n)
        particle(n)%self%q      = q(n)
        particle(n)%self%m      = m(n)
        particle(n)%self%follow = follow(n)
        particle(n)%self%index  = n
        particle(n)%self%t      = 0.0d0
        particle(n)%self%dt     = 0.0d0

        call get_vec(bp, igrid_particle, x(:, n), particle(n)%self%t, B)

        Bnorm = norm2(B(:))
        vnorm = norm2(v(:, n))
        vpar  = sum(v(:, n) * B/Bnorm)
        vperp = sqrt(vnorm**2 - vpar**2)

        ! The momentum vector u(1:3) is filled with the following components

        ! parallel momentum component (gamma v||)
        particle(n)%self%u(1) = lfac * vpar

        ! Mr: the conserved magnetic moment
        magmom = m(n) * (vperp * lfac)**2 / (2.0d0 * Bnorm)
        particle(n)%self%u(2) = magmom

        ! Lorentz factor
        particle(n)%self%u(3) = lfac

        ! initialise payloads for guiding centre module
        allocate(particle(n)%payload(npayload))
        particle(n)%payload(:) = 0.0d0

        if (npayload >= 3) then
          ! gyroradius
          particle(n)%payload(1)=sqrt(2.0d0*m(n)*&
               magmom*Bnorm)/abs(q(n)*Bnorm)*const_c
          ! pitch angle
          particle(n)%payload(2)=atan(vperp, vpar)
          ! perpendicular velocity
          particle(n)%payload(3)=vperp
        end if

      end if

    end do

  end subroutine gca_create_particles

  subroutine gca_fill_gridvars()
    use mod_global_parameters
    use mod_usr_methods, only: usr_particle_fields

    integer                                   :: igrid, iigrid, idir, idim
    double precision, dimension(ixG^T,1:ndir) :: beta
    double precision, dimension(ixG^T,1:nw)   :: w,wold
    double precision                          :: current(ixG^T,7-2*ndir:3)
    integer                                   :: idirmin
    double precision, dimension(ixG^T,1:ndir) :: ue, bhat
    double precision, dimension(ixG^T)        :: kappa, kappa_B, absB, tmp
    double precision                          :: E_field(ixG^T, 1:ndir)
    double precision                          :: B_field(ixG^T, 1:ndir)

    do iigrid=1,igridstail; igrid=igrids(iigrid);

       ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
       gridvars(igrid)%w(ixG^T,1:ngridvars) = 0.0d0

       if (associated(usr_particle_fields)) then
         call usr_particle_fields(pw(igrid)%x, E_field, B_field)
         gridvars(igrid)%w(ixG^T,ep(:)) = E_field
         gridvars(igrid)%w(ixG^T,bp(:)) = B_field
       else
         ! Determine fields from MHD variables
         w(ixG^T,1:nw) = pw(igrid)%w(ixG^T,1:nw)
         call phys_to_primitive(ixG^LL,ixG^LL,w,pw(igrid)%x)

         ! fill with magnetic field:
         gridvars(igrid)%w(ixG^T,bp(:)) = w(ixG^T,iw_mag(:))

         current = zero
         call get_current(w,ixG^LL,ixG^LLIM^D^LSUB1,idirmin,current)
         ! fill electric field
         gridvars(igrid)%w(ixG^T,ep(1)) = gridvars(igrid)%w(ixG^T,bp(2)) * w(ixG^T,iw_mom(3)) &
              - gridvars(igrid)%w(ixG^T,bp(3)) * w(ixG^T,iw_mom(2)) + particles_eta * current(ixG^T,1)
         gridvars(igrid)%w(ixG^T,ep(2)) = gridvars(igrid)%w(ixG^T,bp(3)) * w(ixG^T,iw_mom(1)) &
              - gridvars(igrid)%w(ixG^T,bp(1)) * w(ixG^T,iw_mom(3)) + particles_eta * current(ixG^T,2)
         gridvars(igrid)%w(ixG^T,ep(3)) = gridvars(igrid)%w(ixG^T,bp(1)) * w(ixG^T,iw_mom(2)) &
              - gridvars(igrid)%w(ixG^T,bp(2)) * w(ixG^T,iw_mom(1)) + particles_eta * current(ixG^T,3)

         ! scale to cgs units:
         ! TODO: here unit_magneticfield is used, below something else
         gridvars(igrid)%w(ixG^T,bp(:)) = &
              gridvars(igrid)%w(ixG^T,bp(:)) * unit_magneticfield
         gridvars(igrid)%w(ixG^T,ep(:)) = &
              gridvars(igrid)%w(ixG^T,ep(:)) * unit_magneticfield * unit_velocity / const_c
       end if

       ! grad(kappa B)
       absB(ixG^T) = sqrt(sum(gridvars(igrid)%w(ixG^T,bp(:))**2,dim=ndim+1))
       ue(ixG^T,1) = gridvars(igrid)%w(ixG^T,ep(2)) * gridvars(igrid)%w(ixG^T,bp(3)) &
            - gridvars(igrid)%w(ixG^T,ep(3)) * gridvars(igrid)%w(ixG^T,bp(2))
       ue(ixG^T,2) = gridvars(igrid)%w(ixG^T,ep(3)) * gridvars(igrid)%w(ixG^T,bp(1)) &
            - gridvars(igrid)%w(ixG^T,ep(1)) * gridvars(igrid)%w(ixG^T,bp(3))
       ue(ixG^T,3) = gridvars(igrid)%w(ixG^T,ep(1)) * gridvars(igrid)%w(ixG^T,bp(2)) &
            - gridvars(igrid)%w(ixG^T,ep(2)) * gridvars(igrid)%w(ixG^T,bp(1))
       do idir=1,ndir
         ue(ixG^T,idir) = ue(ixG^T,idir) * const_c / absB(ixG^T)**2
       end do

       kappa(ixG^T) = sqrt(1.0d0 - sum(ue(ixG^T,:)**2,dim=ndim+1)/const_c**2)
       kappa_B(ixG^T) = kappa(ixG^T) * absB(ixG^T)

       do idim=1,ndim
         call gradient(kappa_B,ixG^LL,ixG^LL^LSUB1,idim,tmp)
         gridvars(igrid)%w(ixG^T,grad_kappa_B(idim)) = tmp(ixG^T)/unit_length
       end do

       do idir=1,ndir
         bhat(ixG^T,idir) = gridvars(igrid)%w(ixG^T,bp(idir)) / absB(ixG^T)
       end do

       do idir=1,ndir
         ! (b dot grad) b and the other directional derivatives
         do idim=1,ndim
           call gradient(bhat(ixG^T,idir),ixG^LL,ixG^LL^LSUB1,idim,tmp)
           gridvars(igrid)%w(ixG^T,b_dot_grad_b(idir)) = gridvars(igrid)%w(ixG^T,b_dot_grad_b(idir)) &
                + bhat(ixG^T,idim) * tmp(ixG^T)/unit_length
           gridvars(igrid)%w(ixG^T,ue_dot_grad_b(idir)) = gridvars(igrid)%w(ixG^T,ue_dot_grad_b(idir)) &
                + ue(ixG^T,idim) * tmp(ixG^T)/unit_length
           call gradient(ue(ixG^T,idir),ixG^LL,ixG^LL^LSUB1,idim,tmp)
           gridvars(igrid)%w(ixG^T,b_dot_grad_ue(idir)) = gridvars(igrid)%w(ixG^T,b_dot_grad_ue(idir)) &
                + bhat(ixG^T,idim) * tmp(ixG^T)/unit_length
           gridvars(igrid)%w(ixG^T,ue_dot_grad_ue(idir)) = gridvars(igrid)%w(ixG^T,ue_dot_grad_b(idir)) &
                + ue(ixG^T,idim) * tmp(ixG^T)/unit_length
         end do
       end do

       if(time_advance) then
         gridvars(igrid)%wold(ixG^T,1:ngridvars) = 0.0d0

         if (associated(usr_particle_fields)) then
           !> @todo Compute these fields at a different time?
           call usr_particle_fields(pw(igrid)%x, E_field, B_field)
           gridvars(igrid)%wold(ixG^T,ep(:)) = E_field
           gridvars(igrid)%wold(ixG^T,bp(:)) = B_field
         else
           wold(ixG^T,1:nw) = pw(igrid)%wold(ixG^T,1:nw)
           call phys_to_primitive(ixG^LL,ixG^LL,wold,pw(igrid)%x)
           gridvars(igrid)%wold(ixG^T,bp(:)) = wold(ixG^T,iw_mag(:))
           current = zero
           call get_current(wold,ixG^LL,ixG^LLIM^D^LSUB1,idirmin,current)
           ! fill electric field
           gridvars(igrid)%wold(ixG^T,ep(1)) = gridvars(igrid)%wold(ixG^T,bp(2)) * wold(ixG^T,iw_mom(3)) &
                - gridvars(igrid)%wold(ixG^T,bp(3)) * wold(ixG^T,iw_mom(2)) + particles_eta * current(ixG^T,1)
           gridvars(igrid)%wold(ixG^T,ep(2)) = gridvars(igrid)%wold(ixG^T,bp(3)) * wold(ixG^T,iw_mom(1)) &
                - gridvars(igrid)%wold(ixG^T,bp(1)) * wold(ixG^T,iw_mom(3)) + particles_eta * current(ixG^T,2)
           gridvars(igrid)%wold(ixG^T,ep(3)) = gridvars(igrid)%wold(ixG^T,bp(1)) * wold(ixG^T,iw_mom(2)) &
                - gridvars(igrid)%wold(ixG^T,bp(2)) * wold(ixG^T,iw_mom(1)) + particles_eta * current(ixG^T,3)

           ! scale to cgs units:
           gridvars(igrid)%wold(ixG^T,bp(:)) = &
                gridvars(igrid)%wold(ixG^T,bp(:)) * sqrt(4.0d0*dpi*unit_velocity**2 * unit_density)
           gridvars(igrid)%wold(ixG^T,ep(:)) = &
                gridvars(igrid)%wold(ixG^T,ep(:)) * sqrt(4.0d0*dpi*unit_velocity**2 * unit_density) * unit_velocity / const_c
         end if

         ! grad(kappa B)
         absB(ixG^T) = sqrt(sum(gridvars(igrid)%wold(ixG^T,bp(:))**2,dim=ndim+1))
         ue(ixG^T,1) = gridvars(igrid)%wold(ixG^T,ep(2)) * gridvars(igrid)%wold(ixG^T,bp(3)) &
              - gridvars(igrid)%wold(ixG^T,ep(3)) * gridvars(igrid)%wold(ixG^T,bp(2))
         ue(ixG^T,2) = gridvars(igrid)%wold(ixG^T,ep(3)) * gridvars(igrid)%wold(ixG^T,bp(1)) &
              - gridvars(igrid)%wold(ixG^T,ep(1)) * gridvars(igrid)%wold(ixG^T,bp(3))
         ue(ixG^T,3) = gridvars(igrid)%wold(ixG^T,ep(1)) * gridvars(igrid)%wold(ixG^T,bp(2)) &
              - gridvars(igrid)%wold(ixG^T,ep(2)) * gridvars(igrid)%wold(ixG^T,bp(1))
         do idir=1,ndir
           ue(ixG^T,idir) = ue(ixG^T,idir) * const_c / absB(ixG^T)**2
         end do

         kappa(ixG^T) = sqrt(1.0d0 - sum(ue(ixG^T,:)**2,dim=ndim+1)/const_c**2)
         kappa_B(ixG^T) = kappa(ixG^T) * absB(ixG^T)

         do idim=1,ndim
           call gradient(kappa_B,ixG^LL,ixG^LL^LSUB1,idim,tmp)
           gridvars(igrid)%wold(ixG^T,grad_kappa_B(idim)) = tmp(ixG^T)/unit_length
         end do

         do idir=1,ndir
           bhat(ixG^T,idir) = gridvars(igrid)%wold(ixG^T,bp(idir)) / absB(ixG^T)
         end do

         do idir=1,ndir
           ! (b dot grad) b and the other directional derivatives
           do idim=1,ndim
             call gradient(bhat(ixG^T,idir),ixG^LL,ixG^LL^LSUB1,idim,tmp)
             gridvars(igrid)%wold(ixG^T,b_dot_grad_b(idir)) = gridvars(igrid)%wold(ixG^T,b_dot_grad_b(idir)) &
                  + bhat(ixG^T,idim) * tmp(ixG^T)/unit_length
             gridvars(igrid)%wold(ixG^T,ue_dot_grad_b(idir)) = gridvars(igrid)%wold(ixG^T,ue_dot_grad_b(idir)) &
                  + ue(ixG^T,idim) * tmp(ixG^T)/unit_length
             call gradient(ue(ixG^T,idir),ixG^LL,ixG^LL^LSUB1,idim,tmp)
             gridvars(igrid)%wold(ixG^T,b_dot_grad_ue(idir)) = gridvars(igrid)%wold(ixG^T,b_dot_grad_ue(idir)) &
                  + bhat(ixG^T,idim) * tmp(ixG^T)/unit_length
             gridvars(igrid)%wold(ixG^T,ue_dot_grad_ue(idir)) = gridvars(igrid)%wold(ixG^T,ue_dot_grad_b(idir)) &
                  + ue(ixG^T,idim) * tmp(ixG^T)/unit_length
           end do
         end do
       end if
    end do

  end subroutine gca_fill_gridvars

    subroutine gca_integrate_particles(end_time)
    use mod_odeint
    use mod_global_parameters
    double precision, intent(in)        :: end_time

    double precision                    :: lfac, absS
    double precision                    :: dt_p, tloc, y(ndir+2),dydt(ndir+2),ytmp(ndir+2), euler_cfl, int_factor
    double precision, dimension(1:ndir) :: x, ue, e, b, bhat, x_new
    double precision, dimension(1:ndir) :: drift1, drift2
    double precision, dimension(1:ndir) :: drift3, drift4, drift5, drift6, drift7
    double precision, dimension(1:ndir) :: bdotgradb, uedotgradb, gradkappaB
    double precision, dimension(1:ndir) :: bdotgradue, uedotgradue
    double precision, dimension(1:ndir) :: gradBdrift, reldrift, bdotgradbdrift
    double precision, dimension(1:ndir) :: uedotgradbdrift, bdotgraduedrift
    double precision, dimension(1:ndir) :: uedotgraduedrift
    double precision                    :: kappa, Mr, upar, m, absb, gamma, q, mompar, vpar, ueabs
    double precision                    :: gradBdrift_abs, reldrift_abs, epar
    double precision                    :: bdotgradbdrift_abs, uedotgradbdrift_abs
    double precision                    :: bdotgraduedrift_abs, uedotgraduedrift_abs
    double precision                    :: momentumpar1, momentumpar2, momentumpar3, momentumpar4
    ! Precision of time-integration:
    double precision,parameter          :: eps=1.0d-6
    ! for odeint:
    double precision                    :: h1, hmin, h_old
    integer                             :: nok, nbad, ic1^D, ic2^D, ierror, nvar
    integer                             :: ipart, iipart, seed, ic^D,igrid_particle, ipe_particle, ipe_working
    logical                             :: int_choice
    logical                             :: BC_applied

    nvar=ndir+2

    do iipart=1,nparticles_active_on_mype
      ipart                   = particles_active_on_mype(iipart)
      int_choice              = .false.
      dt_p                    = gca_get_particle_dt(particle(ipart), end_time)
      particle(ipart)%self%dt = dt_p

      igrid_working           = particle(ipart)%igrid
      ipart_working           = particle(ipart)%self%index
      tloc                    = particle(ipart)%self%t
      x(1:ndir)               = particle(ipart)%self%x(1:ndir)

      ! Adaptive stepwidth RK4:
      ! initial solution vector:
      y(1:ndir) = x(1:ndir) ! position of guiding center
      y(ndir+1) = particle(ipart)%self%u(1) ! parallel momentum component (gamma v||)
      y(ndir+2) = particle(ipart)%self%u(2) ! conserved magnetic moment Mr
    ! y(ndir+3) = particle(ipart)%self%u(3) ! Lorentz factor of particle

      ! we temporarily save the solution vector, to replace the one from the euler
      ! timestep after euler integration
      ytmp=y

      call derivs_gca(particle(ipart)%self%t,y,dydt)

      ! make an Euler step with the proposed timestep:
      ! factor to ensure we capture all particles near the internal ghost cells.
      ! Can be adjusted during a run, after an interpolation error.
      euler_cfl=2.5d0

      ! new solution vector:
      y(1:ndir+2) = y(1:ndir+2) + euler_cfl * dt_p * dydt(1:ndir+2)
      particle(ipart)%self%x(1:ndir) = y(1:ndir) ! position of guiding center
      particle(ipart)%self%u(1)      = y(ndir+1) ! parallel momentum component(gamma v||)
      particle(ipart)%self%u(2)      = y(ndir+2) ! conserved magnetic moment

      ! check if the particle is in the internal ghost cells
      int_factor =1.0d0

      if(.not. particle_in_igrid(ipart_working,igrid_working)) then
        ! if particle is not in the grid an euler timestep is taken instead of a RK4
        ! timestep. Then based on that we do an interpolation and check how much further
        ! the timestep for the RK4 has to be restricted.
        ! factor to make integration more accurate for particles near the internal
        ! ghost cells. This factor can be changed during integration after an
        ! interpolation error. But one should be careful with timesteps for i/o

        ! flat interpolation:
        {ic^D = int((y(^D)-rnode(rpxmin^D_,igrid_working))/rnode(rpdx^D_,igrid_working)) + 1 + nghostcells\}

        ! linear interpolation:
        {
        if (pw(igrid_working)%x({ic^DD},^D) .lt. y(^D)) then
           ic1^D = ic^D
        else
           ic1^D = ic^D -1
        end if
        ic2^D = ic1^D + 1
        \}

        int_factor =0.5d0

        {^D&
        if (ic1^D .le. ixGlo^D-2 .or. ic2^D .ge. ixGhi^D+2) then
          int_factor = 0.05d0
        end if
        \}

        {^D&
        if (ic1^D .eq. ixGlo^D-1 .or. ic2^D .eq. ixGhi^D+1) then
          int_factor = 0.1d0
        end if
        \}

        dt_p=int_factor*dt_p
      end if

      ! replace the solution vector with the original as it was before the Euler timestep
      y(1:ndir+2) = ytmp(1:ndir+2)

      particle(ipart)%self%x(1:ndir) = ytmp(1:ndir) ! position of guiding center
      particle(ipart)%self%u(1)      = ytmp(ndir+1) ! parallel momentum component (gamma v||)
      particle(ipart)%self%u(2)      = ytmp(ndir+2) ! conserved magnetic moment

      ! specify a minimum step hmin. If the timestep reaches this minimum, multiply by
      ! a factor 100 to make sure the RK integration doesn't crash
      h1 = dt_p/2.0d0; hmin=1.0d-9; h_old=dt_p/2.0d0

      if(h1 .lt. hmin)then
        h1=hmin
        dt_p=2.0d0*h1
      endif

      ! RK4 integration with adaptive stepwidth
      call odeint(y,nvar,tloc,tloc+dt_p,eps,h1,hmin,nok,nbad,derivs_gca_rk,rkqs,ierror)

      if (ierror /= 0) then
         print *, "odeint returned error code", ierror
         print *, "1 means hmin too small, 2 means MAXSTP exceeded"
         print *, "Having a problem with particle", iipart
      end if

      ! original RK integration without interpolation in ghost cells
      ! call odeint(y,nvar,tloc,tloc+dt_p,eps,h1,hmin,nok,nbad,derivs_gca,rkqs)

      ! final solution vector after rk integration
      particle(ipart)%self%x(1:ndir) = y(1:ndir)
      particle(ipart)%self%u(1)      = y(ndir+1)
      particle(ipart)%self%u(2)      = y(ndir+2)
      !particle(ipart)%self%u(3)      = y(ndir+3)

      ! now calculate other quantities, mean Lorentz factor, drifts, perpendicular velocity:
      call get_vec(bp, igrid_working,y(1:ndir),tloc+dt_p,b)
      call get_vec(ep, igrid_working,y(1:ndir),tloc+dt_p,e)

      absb         = sqrt(sum(b(:)**2))
      bhat(1:ndir) = b(1:ndir) / absb

      epar         = sum(e(:)*bhat(:))

      call cross(e,bhat,ue)

      ue(1:ndir)   = ue(1:ndir)*const_c / absb
      ueabs = sqrt(sum(ue(:)**2))
      kappa = sqrt(1.0d0 - sum(ue(:)**2)/const_c**2)
      Mr = y(ndir+2); upar = y(ndir+1); m=particle(ipart)%self%m; q=particle(ipart)%self%q
      gamma = sqrt(1.0d0+upar**2/const_c**2+2.0d0*Mr*absb/m/const_c**2)/kappa

      particle(ipart)%self%u(3)      = gamma

      vpar = particle(ipart)%self%u(1)/particle(ipart)%self%u(3)
      mompar = particle(ipart)%self%u(1)

      call get_vec(b_dot_grad_b, igrid_working,y(1:ndir),tloc+dt_p,bdotgradb)
      call get_vec(ue_dot_grad_b, igrid_working,y(1:ndir),tloc+dt_p,uedotgradb)
      call get_vec(grad_kappa_B, igrid_working,y(1:ndir),tloc+dt_p,gradkappaB)
      call get_vec(b_dot_grad_ue, igrid_working,y(1:ndir),tloc+dt_p,bdotgradue)
      call get_vec(ue_dot_grad_ue, igrid_working,y(1:ndir),tloc+dt_p,uedotgradue)

      drift1(1:ndir) = bhat(1:ndir)/(absb*kappa**2)
      drift2(1:ndir) = Mr*const_c/(gamma*q)*gradkappaB(1:ndir)

      call cross(drift1,drift2,gradBdrift)
      gradBdrift_abs = sqrt(sum(gradBdrift(:)**2))

      drift3(1:ndir) = upar*epar/(gamma*const_c)*ue(1:ndir)
      call cross(drift1,drift3,reldrift)
      reldrift_abs = sqrt(sum(reldrift(:)**2))

      drift4(1:ndir) = m*const_c/q* ( upar**2/gamma*bdotgradb(1:ndir))
      call cross(drift1,drift4,bdotgradbdrift)
      bdotgradbdrift_abs = sqrt(sum(bdotgradbdrift(:)**2))

      drift5(1:ndir) = m*const_c/q* ( upar*uedotgradb(1:ndir))
      call cross(drift1,drift5,uedotgradbdrift)
      uedotgradbdrift_abs = sqrt(sum(uedotgradbdrift(:)**2))

      drift6(1:ndir) = m*const_c/q* ( upar*bdotgradue(1:ndir))
      call cross(drift1,drift6,bdotgraduedrift)
      bdotgraduedrift_abs = sqrt(sum(bdotgraduedrift(:)**2))

      drift7(1:ndir) = m*const_c/q* (gamma*uedotgradue(1:ndir))
      call cross(drift1,drift7,uedotgraduedrift)
      uedotgraduedrift_abs = sqrt(sum(uedotgraduedrift(:)**2))

      momentumpar1 = m*gamma*vpar**2*sum(ue(:)*bdotgradb(:))
      momentumpar2 = m*gamma*vpar*sum(ue(:)*uedotgradb(:))
      momentumpar3 = q*epar
      momentumpar4 = -(Mr/gamma)*sum(bhat(:)*gradkappaB(:))

      ! Payload update
      if (npayload >= 14) then
        ! current gyroradius
        particle(ipart)%payload(1) = sqrt(2.0d0*m*Mr*absb)/abs(q*absb)*const_c
        ! pitch angle
        particle(ipart)%payload(2) = atan(sqrt((2.0d0*Mr*absb)/(m*gamma**2)), vpar)
        ! particle v_perp
        particle(ipart)%payload(3) = sqrt((2.0d0*Mr*absb)/(m*gamma**2))
        ! particle parallel momentum term 1
        particle(ipart)%payload(4) = momentumpar1
        ! particle parallel momentum term 2
        particle(ipart)%payload(5) = momentumpar2
        ! particle parallel momentum term 3
        particle(ipart)%payload(6) = momentumpar3
        ! particle parallel momentum term 4
        particle(ipart)%payload(7) = momentumpar4
        ! particle ExB drift
        particle(ipart)%payload(8) = ueabs
        ! relativistic drift
        particle(ipart)%payload(9) = reldrift_abs
        ! gradB drift
        particle(ipart)%payload(10) = gradBdrift_abs
        ! bdotgradb drift
        particle(ipart)%payload(11) = bdotgradbdrift_abs
        ! uedotgradb drift
        particle(ipart)%payload(12) = uedotgradbdrift_abs
        ! bdotgradue drift
        particle(ipart)%payload(13) = bdotgraduedrift_abs
        ! uedotgradue drift
        particle(ipart)%payload(14) = uedotgraduedrift_abs
      end if

      ! Time update
      particle(ipart)%self%t = particle(ipart)%self%t + dt_p

    end do

  end subroutine gca_integrate_particles

  subroutine derivs_gca_rk(t_s,y,dydt)
    use mod_global_parameters

    double precision                :: t_s, y(ndir+2)
    double precision                :: dydt(ndir+2)

    double precision,dimension(ndir):: ue, b, e, x, bhat, bdotgradb, uedotgradb, gradkappaB
    double precision,dimension(ndir):: bdotgradue, uedotgradue, u, utmp1, utmp2, utmp3
    double precision                :: upar, Mr, gamma, absb, q, m, epar, kappa
    integer                         :: ic^D

    ! Here the terms in the guiding centre equations of motion are interpolated for
    ! the RK integration. The interpolation is also done in the ghost cells such
    ! that the RK integration does not give an error

    q = particle(ipart_working)%self%q
    m = particle(ipart_working)%self%m

    x(1:ndir) = y(1:ndir)
    upar      = y(ndir+1) ! gamma v||
    Mr        = y(ndir+2)
    !gamma     = y(ndir+3)

    call get_vec(bp, igrid_working,x,t_s,b)
    call get_vec(ep, igrid_working,x,t_s,e)
    call get_vec(b_dot_grad_b, igrid_working,x,t_s,bdotgradb)
    call get_vec(ue_dot_grad_b, igrid_working,x,t_s,uedotgradb)
    call get_vec(grad_kappa_B, igrid_working,x,t_s,gradkappaB)
    call get_vec(b_dot_grad_ue, igrid_working,x,t_s,bdotgradue)
    call get_vec(ue_dot_grad_ue, igrid_working,x,t_s,uedotgradue)

    absb         = sqrt(sum(b(:)**2))
    bhat(1:ndir) = b(1:ndir) / absb
    epar         = sum(e(:)*bhat(:))

    call cross(e,bhat,ue)
    ue(1:ndir)   = ue(1:ndir)*const_c / absb

    kappa = sqrt(1.0d0 - sum(ue(:)**2)/const_c**2)
    gamma = sqrt(1.0d0+upar**2/const_c**2+2.0d0*Mr*absb/m/const_c**2)/kappa

    utmp1(1:ndir) = bhat(1:ndir)/(absb*kappa**2)
    utmp2(1:ndir) = Mr*const_c/(gamma*q)*gradkappaB(1:ndir) &
     + upar*epar/(gamma*const_c)*ue(1:ndir) &
         + m*const_c/q* ( upar**2/gamma*bdotgradb(1:ndir) + upar*uedotgradb(1:ndir) &
         + upar*bdotgradue(1:ndir) + gamma*uedotgradue(1:ndir))

    call cross(utmp1,utmp2,utmp3)
    u(1:ndir) = ue(1:ndir) + utmp3(1:ndir)

    ! done assembling the terms, now write rhs:
    dydt(1:ndir) = ( u(1:ndir) + upar/gamma * bhat(1:ndir) )/ unit_length
    dydt(ndir+1) = sum(ue(:)*(upar*bdotgradb(:)+gamma*uedotgradb(:))) &
         + q/m*epar - Mr/(m*gamma) * sum(bhat(:)*gradkappaB(:))

    dydt(ndir+2) = 0.0d0 ! magnetic moment is conserved
    !dydt(ndir+3) = q/(m*const_c**2) * ({^C& dydt(^C)*e(^C)|+}) * unit_length

  end subroutine derivs_gca_rk

  subroutine derivs_gca(t_s,y,dydt)
    use mod_global_parameters

    double precision                :: t_s, y(ndir+2)
    double precision                :: dydt(ndir+2)

    double precision,dimension(ndir):: ue, b, e, x, bhat, bdotgradb, uedotgradb, gradkappaB
    double precision,dimension(ndir):: bdotgradue, uedotgradue, u, utmp1, utmp2, utmp3
    double precision                :: upar, Mr, gamma, absb, q, m, epar, kappa

    ! Here the normal interpolation is done for the terms in the GCA equations of motion

    q = particle(ipart_working)%self%q
    m = particle(ipart_working)%self%m

    x(1:ndir) = y(1:ndir)
    upar      = y(ndir+1) ! gamma v||
    Mr        = y(ndir+2)
    !gamma     = y(ndir+3)

    call get_vec(bp, igrid_working,x,t_s,b)
    call get_vec(ep, igrid_working,x,t_s,e)
    call get_vec(b_dot_grad_b, igrid_working,x,t_s,bdotgradb)
    call get_vec(ue_dot_grad_b, igrid_working,x,t_s,uedotgradb)
    call get_vec(grad_kappa_B, igrid_working,x,t_s,gradkappaB)
    call get_vec(b_dot_grad_ue, igrid_working,x,t_s,bdotgradue)
    call get_vec(ue_dot_grad_ue, igrid_working,x,t_s,uedotgradue)

    absb         = sqrt(sum(b(:)**2))
    bhat(1:ndir) = b(1:ndir) / absb

    epar         = sum(e(:)*bhat(:))
    call cross(e,bhat,ue)
    ue(1:ndir)   = ue(1:ndir)*const_c / absb

    kappa = sqrt(1.0d0 - sum(ue(:)**2)/const_c**2)
    gamma = sqrt(1.0d0+upar**2/const_c**2+2.0d0*Mr*absb/m/const_c**2)/kappa
    utmp1(1:ndir) = bhat(1:ndir)/(absb*kappa**2)
    utmp2(1:ndir) = Mr*const_c/(gamma*q)*gradkappaB(1:ndir) &
     + upar*epar/(gamma*const_c)*ue(1:ndir) &
         + m*const_c/q* ( upar**2/gamma*bdotgradb(1:ndir) + upar*uedotgradb(1:ndir) &
         + upar*bdotgradue(1:ndir) + gamma*uedotgradue(1:ndir))

    call cross(utmp1,utmp2,utmp3)
    u(1:ndir) = ue(1:ndir) + utmp3(1:ndir)

    ! done assembling the terms, now write rhs:
    dydt(1:ndir) = ( u(1:ndir) + upar/gamma * bhat(1:ndir) )/ unit_length
    dydt(ndir+1) = sum(ue(:)*(upar*bdotgradb(:)+gamma*uedotgradb(:))) &
         + q/m*epar - Mr/(m*gamma) * sum(bhat(:)*gradkappaB(:))
    dydt(ndir+2) = 0.0d0 ! magnetic moment is conserved
    !dydt(ndir+3) = q/(m*const_c**2) * ({^C& dydt(^C)*e(^C)|+}) * unit_length

  end subroutine derivs_gca

  function gca_get_particle_dt(partp, end_time) result(dt_p)
    use mod_odeint
    use mod_global_parameters
    type(particle_ptr), intent(in) :: partp
    double precision, intent(in)   :: end_time
    double precision               :: dt_p

    double precision            :: tout, dt_particles_mype, dt_cfl0, dt_cfl1, dt_a
    double precision            :: dxmin, vp, a, gammap
    double precision            :: v(ndir), y(ndir+2),ytmp(ndir+2), dydt(ndir+2), v0(ndir), v1(ndir), dydt1(ndir+2)
    double precision            :: ap0, ap1, dt_cfl_ap0, dt_cfl_ap1, dt_cfl_ap
    double precision            :: dt_euler, dt_tmp
    ! make these particle cfl conditions more restrictive if you are interpolating out of the grid
    double precision, parameter :: cfl=0.8d0, uparcfl=0.8d0
    double precision, parameter :: uparmin=1.0d-6*const_c
    integer                     :: ipart, iipart, nout, ic^D, igrid_particle, ipe_particle, ipe
    logical                     :: BC_applied

    if (const_dt_particles > 0) then
      dt_p = const_dt_particles
      return
    end if

    igrid_working = partp%igrid
    ipart_working = partp%self%index
    dt_tmp = (end_time - partp%self%t)
    if(dt_tmp .le. 0.0d0) dt_tmp = smalldouble
    ! make sure we step only one cell at a time, first check CFL at current location
    ! then we make an Euler step to the new location and check the new CFL
    ! we simply take the minimum of the two timesteps.
    ! added safety factor cfl:
    dxmin  = min({rnode(rpdx^D_,partp%igrid)},bigdouble)*cfl
    ! initial solution vector:
    y(1:ndir) = partp%self%x(1:ndir) ! position of guiding center
    y(ndir+1) = partp%self%u(1) ! parallel momentum component (gamma v||)
    y(ndir+2) = partp%self%u(2) ! conserved magnetic moment
    ytmp=y
    !y(ndir+3) = partp%self%u(3) ! Lorentz factor of guiding centre

    call derivs_gca(partp%self%t,y,dydt)
    v0(1:ndir) = dydt(1:ndir)
    ap0        = dydt(ndir+1)

    ! guiding center velocity:
    v(1:ndir) = abs(dydt(1:ndir))
    vp = sqrt(sum(v(:)**2))

    dt_cfl0    = dxmin / max(vp, smalldouble)
    dt_cfl_ap0 = uparcfl * abs(max(abs(y(ndir+1)),uparmin) / max(ap0, smalldouble))
    !dt_cfl_ap0 = min(dt_cfl_ap0, uparcfl * sqrt(abs(unit_length*dxmin/(ap0+smalldouble))) )

    ! make an Euler step with the proposed timestep:
    ! new solution vector:
    dt_euler = min(dt_tmp,dt_cfl0,dt_cfl_ap0)
    y(1:ndir+2) = y(1:ndir+2) + dt_euler * dydt(1:ndir+2)

    partp%self%x(1:ndir) = y(1:ndir) ! position of guiding center
    partp%self%u(1)      = y(ndir+1) ! parallel momentum component (gamma v||)
    partp%self%u(2)      = y(ndir+2) ! conserved magnetic moment

    ! first check if the particle is outside the physical domain or in the ghost cells
    if(.not. particle_in_igrid(ipart_working,igrid_working)) then
      y(1:ndir+2) = ytmp(1:ndir+2)
    end if

    call derivs_gca_rk(partp%self%t+dt_euler,y,dydt)
    !call derivs_gca(partp%self%t+dt_euler,y,dydt)

    v1(1:ndir) = dydt(1:ndir)
    ap1        = dydt(ndir+1)

    ! guiding center velocity:
    v(1:ndir) = abs(dydt(1:ndir))
    vp = sqrt(sum(v(:)**2))

    dt_cfl1    = dxmin / max(vp, smalldouble)
    dt_cfl_ap1 = uparcfl * abs(max(abs(y(ndir+1)),uparmin) / max(ap1, smalldouble))
    !dt_cfl_ap1 = min(dt_cfl_ap1, uparcfl * sqrt(abs(unit_length*dxmin/(ap1+smalldouble))) )

    dt_tmp = min(dt_euler, dt_cfl1, dt_cfl_ap1)

    partp%self%x(1:ndir) = ytmp(1:ndir) ! position of guiding center
    partp%self%u(1)      = ytmp(ndir+1) ! parallel momentum component (gamma v||)
    partp%self%u(2)      = ytmp(ndir+2) ! conserved magnetic moment
    !dt_tmp = min(dt_cfl1, dt_cfl_ap1)

    ! time step due to parallel acceleration:
    ! The standart thing, dt=sqrt(dx/a) where we comupte a from d(gamma v||)/dt and d(gamma)/dt
    ! dt_ap = sqrt(abs(dxmin*unit_length*y(ndir+3)/( dydt(ndir+1) - y(ndir+1)/y(ndir+3)*dydt(ndir+3) ) ) )
    ! vp = sqrt({^C& (v(^C)*unit_length)**2|+})
    ! gammap = sqrt(1.0d0/(1.0d0-(vp/const_c)**2))
    ! ap = const_c**2/vp*gammap**(-3)*dydt(ndir+3)
    ! dt_ap = sqrt(dxmin*unit_length/ap)

    !dt_a = bigdouble
    !if (dt_euler .gt. smalldouble) then
    !   a = sqrt({^C& (v1(^C)-v0(^C))**2 |+})/dt_euler
    !   dt_a = min(sqrt(dxmin/a),bigdouble)
    !end if

    !dt_p = min(dt_tmp , dt_a)
    dt_p = dt_tmp

    ! Make sure we don't advance beyond end_time
    call limit_dt_endtime(end_time - partp%self%t, dt_p)

  end function gca_get_particle_dt


end module mod_particle_gca
