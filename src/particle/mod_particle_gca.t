!> Particle mover with Newtonian/relativistic Guiding Center Approximation (GCA)
!> By Jannis Teunissen, Bart Ripperda, Oliver Porth, and Fabio Bacchini (2016-2020)
module mod_particle_gca
  use mod_particle_base

  private

  !> Variable index for gradient B, with relativistic correction 1/kappa
  !> where kappa = 1/sqrt(1 - E_perp^2/B^2)
  integer, protected, allocatable      :: grad_kappa_B(:)
  !> Variable index for (B . grad)B (curvature B drift)
  integer, protected, allocatable      :: b_dot_grad_b(:)

  ! ExB related drifts (vE = ExB/B^2)
  !> Variable index for curvature drift
  integer, protected, allocatable      :: vE_dot_grad_b(:)
  !> Variable index for polarization drift
  integer, protected, allocatable      :: b_dot_grad_vE(:)
  !> Variable index for polarization drift
  integer, protected, allocatable      :: vE_dot_grad_vE(:)

  public :: gca_init
  public :: gca_create_particles

  ! Variables
  public :: bp, ep, grad_kappa_B, b_dot_grad_b
  public :: vE_dot_grad_b, b_dot_grad_vE, vE_dot_grad_vE

contains

  subroutine gca_init()
    use mod_global_parameters
    integer :: idir, nwx

    if (physics_type/='mhd') call mpistop("GCA particles need magnetic field!")
    if (ndir/=3) call mpistop("GCA particles need ndir=3!")

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
    allocate(vE_dot_grad_b(ndir))
    do idir = 1, ndir
      nwx = nwx + 1
      vE_dot_grad_b(idir) = nwx
    end do
    allocate(b_dot_grad_vE(ndir))
    do idir = 1, ndir
      nwx = nwx + 1
      b_dot_grad_vE(idir) = nwx
    end do
    allocate(vE_dot_grad_vE(ndir))
    do idir = 1, ndir
      nwx = nwx + 1
      vE_dot_grad_vE(idir) = nwx
    end do
    ngridvars=nwx

    particles_fill_gridvars => gca_fill_gridvars
    particles_integrate     => gca_integrate_particles
  end subroutine gca_init

  subroutine gca_create_particles()
    ! initialise the particles
    use mod_global_parameters
    use mod_usr_methods, only: usr_create_particles, usr_update_payload

    double precision :: b(ndir), u(ndir), magmom
    double precision :: bnorm, lfac, vnorm, vperp, vpar
    integer          :: igrid, ipe_particle
    integer          :: n, idir
    double precision :: x(3, num_particles)
    double precision :: v(3, num_particles)
    double precision :: q(num_particles)
    double precision :: m(num_particles)
    double precision :: rrd(num_particles,ndir)
    double precision :: payload(npayload)
    logical          :: follow(num_particles)

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

    ! first find ipe and igrid responsible for particle
    do n = 1, num_particles
      call find_particle_ipe(x(:, n),igrid,ipe_particle)

      particle(n)%igrid = igrid
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
        particle(n)%self%time   = global_time
        particle(n)%self%dt     = 0.0d0

        call get_vec(bp, igrid, x(:, n), particle(n)%self%time, b)

        bnorm = norm2(b(:))
        vnorm = norm2(v(:, n))
        vpar  = sum(v(:, n) * b/bnorm)
        vperp = sqrt(vnorm**2 - vpar**2)

        ! The momentum vector u(1:3) is filled with the following components

        ! parallel momentum component (gamma v||)
        particle(n)%self%u(1) = lfac * vpar

        ! Mr: the conserved magnetic moment
        magmom = m(n) * (vperp * lfac)**2 / (2.0d0 * bnorm)
        particle(n)%self%u(2) = magmom

        ! Lorentz factor
        particle(n)%self%u(3) = lfac

        ! initialise payloads for GCA module
        allocate(particle(n)%payload(npayload))
        if (.not. associated(usr_update_payload)) then
          call gca_update_payload(igrid,ps(igrid)%w,pso(igrid)%w,ps(igrid)%x,x(:,n),v(:,n)*lfac,q(n),m(n),payload,npayload,0.d0)
        else
          call usr_update_payload(igrid,ps(igrid)%w,pso(igrid)%w,ps(igrid)%x,x(:,n),v(:,n)*lfac,q(n),m(n),payload,npayload,0.d0)
        end if
        particle(n)%payload(:) = payload
      end if
    end do
  end subroutine gca_create_particles

  subroutine gca_fill_gridvars
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
       ! grad(kappa B)
       absB(ixG^T) = sqrt(sum(gridvars(igrid)%w(ixG^T,bp(:))**2,dim=ndim+1))
       vE(ixG^T,1) = gridvars(igrid)%w(ixG^T,ep(2)) * gridvars(igrid)%w(ixG^T,bp(3)) &
            - gridvars(igrid)%w(ixG^T,ep(3)) * gridvars(igrid)%w(ixG^T,bp(2))
       vE(ixG^T,2) = gridvars(igrid)%w(ixG^T,ep(3)) * gridvars(igrid)%w(ixG^T,bp(1)) &
            - gridvars(igrid)%w(ixG^T,ep(1)) * gridvars(igrid)%w(ixG^T,bp(3))
       vE(ixG^T,3) = gridvars(igrid)%w(ixG^T,ep(1)) * gridvars(igrid)%w(ixG^T,bp(2)) &
            - gridvars(igrid)%w(ixG^T,ep(2)) * gridvars(igrid)%w(ixG^T,bp(1))
       do idir=1,ndir
         vE(ixG^T,idir) = vE(ixG^T,idir) / absB(ixG^T)**2
       end do

       if (relativistic) then
         kappa(ixG^T) = 1.d0/sqrt(1.0d0 - sum(vE(ixG^T,:)**2,dim=ndim+1)/c_norm**2)
       else
         kappa(ixG^T) = 1.d0
       end if
       kappa_B(ixG^T) = absB(ixG^T) / kappa(ixG^T)

       do idim=1,ndim
         call gradient(kappa_B,ixG^LL,ixG^LL^LSUB1,idim,tmp)
         gridvars(igrid)%w(ixG^T,grad_kappa_B(idim)) = tmp(ixG^T)
       end do

       do idir=1,ndir
         bhat(ixG^T,idir) = gridvars(igrid)%w(ixG^T,bp(idir)) / absB(ixG^T)
       end do

       do idir=1,ndir
         ! (b dot grad) b and the other directional derivatives
         do idim=1,ndim
           call gradient(bhat(ixG^T,idir),ixG^LL,ixG^LL^LSUB1,idim,tmp)
           gridvars(igrid)%w(ixG^T,b_dot_grad_b(idir)) = gridvars(igrid)%w(ixG^T,b_dot_grad_b(idir)) &
                + bhat(ixG^T,idim) * tmp(ixG^T)
           gridvars(igrid)%w(ixG^T,vE_dot_grad_b(idir)) = gridvars(igrid)%w(ixG^T,vE_dot_grad_b(idir)) &
                + vE(ixG^T,idim) * tmp(ixG^T)
           call gradient(vE(ixG^T,idir),ixG^LL,ixG^LL^LSUB1,idim,tmp)
           gridvars(igrid)%w(ixG^T,b_dot_grad_vE(idir)) = gridvars(igrid)%w(ixG^T,b_dot_grad_vE(idir)) &
                + bhat(ixG^T,idim) * tmp(ixG^T)
           gridvars(igrid)%w(ixG^T,vE_dot_grad_vE(idir)) = gridvars(igrid)%w(ixG^T,vE_dot_grad_vE(idir)) &
                + vE(ixG^T,idim) * tmp(ixG^T)
         end do
       end do

       if(time_advance) then
         ! grad(kappa B)
         absB(ixG^T) = sqrt(sum(gridvars(igrid)%wold(ixG^T,bp(:))**2,dim=ndim+1))
         vE(ixG^T,1) = gridvars(igrid)%wold(ixG^T,ep(2)) * gridvars(igrid)%wold(ixG^T,bp(3)) &
              - gridvars(igrid)%wold(ixG^T,ep(3)) * gridvars(igrid)%wold(ixG^T,bp(2))
         vE(ixG^T,2) = gridvars(igrid)%wold(ixG^T,ep(3)) * gridvars(igrid)%wold(ixG^T,bp(1)) &
              - gridvars(igrid)%wold(ixG^T,ep(1)) * gridvars(igrid)%wold(ixG^T,bp(3))
         vE(ixG^T,3) = gridvars(igrid)%wold(ixG^T,ep(1)) * gridvars(igrid)%wold(ixG^T,bp(2)) &
              - gridvars(igrid)%wold(ixG^T,ep(2)) * gridvars(igrid)%wold(ixG^T,bp(1))
         do idir=1,ndir
           vE(ixG^T,idir) = vE(ixG^T,idir) / absB(ixG^T)**2
         end do

         
         if (relativistic) then
           kappa(ixG^T) = 1.d0/sqrt(1.0d0 - sum(vE(ixG^T,:)**2,dim=ndim+1)/c_norm**2)
         else
           kappa(ixG^T) = 1.d0
         end if
         kappa_B(ixG^T) = absB(ixG^T) / kappa(ixG^T)

         do idim=1,ndim
           call gradient(kappa_B,ixG^LL,ixG^LL^LSUB1,idim,tmp)
           gridvars(igrid)%wold(ixG^T,grad_kappa_B(idim)) = tmp(ixG^T)
         end do

         do idir=1,ndir
           bhat(ixG^T,idir) = gridvars(igrid)%wold(ixG^T,bp(idir)) / absB(ixG^T)
         end do

         do idir=1,ndir
           ! (b dot grad) b and the other directional derivatives
           do idim=1,ndim
             call gradient(bhat(ixG^T,idir),ixG^LL,ixG^LL^LSUB1,idim,tmp)
             gridvars(igrid)%wold(ixG^T,b_dot_grad_b(idir)) = gridvars(igrid)%wold(ixG^T,b_dot_grad_b(idir)) &
                  + bhat(ixG^T,idim) * tmp(ixG^T)
             gridvars(igrid)%wold(ixG^T,vE_dot_grad_b(idir)) = gridvars(igrid)%wold(ixG^T,vE_dot_grad_b(idir)) &
                  + vE(ixG^T,idim) * tmp(ixG^T)
             call gradient(vE(ixG^T,idir),ixG^LL,ixG^LL^LSUB1,idim,tmp)
             gridvars(igrid)%wold(ixG^T,b_dot_grad_vE(idir)) = gridvars(igrid)%wold(ixG^T,b_dot_grad_vE(idir)) &
                  + bhat(ixG^T,idim) * tmp(ixG^T)
             gridvars(igrid)%wold(ixG^T,vE_dot_grad_vE(idir)) = gridvars(igrid)%wold(ixG^T,vE_dot_grad_vE(idir)) &
                  + vE(ixG^T,idim) * tmp(ixG^T)
           end do
         end do
       end if
    end do

  end subroutine gca_fill_gridvars

  subroutine gca_integrate_particles(end_time)
    use mod_odeint
    use mod_global_parameters
    use mod_usr_methods, only: usr_create_particles, usr_update_payload
    double precision, intent(in)        :: end_time

    double precision                    :: lfac, absS
    double precision                    :: payload(npayload)
    double precision                    :: dt_p, tloc, y(ndir+2),dydt(ndir+2),ytmp(ndir+2), euler_cfl, int_factor
    double precision, dimension(1:ndir) :: x, vE, e, b, bhat, x_new
    double precision, dimension(1:ndir) :: drift1, drift2
    double precision, dimension(1:ndir) :: drift3, drift4, drift5, drift6, drift7
    double precision, dimension(1:ndir) :: bdotgradb, vEdotgradb, gradkappaB
    double precision, dimension(1:ndir) :: bdotgradvE, vEdotgradvE
    double precision, dimension(1:ndir) :: gradBdrift, reldrift, bdotgradbdrift
    double precision, dimension(1:ndir) :: vEdotgradbdrift, bdotgradvEdrift
    double precision, dimension(1:ndir) :: vEdotgradvEdrift
    double precision                    :: kappa, Mr, upar, m, absb, gamma, q, mompar, vpar, vEabs
    double precision                    :: gradBdrift_abs, reldrift_abs, epar
    double precision                    :: bdotgradbdrift_abs, vEdotgradbdrift_abs
    double precision                    :: bdotgradvEdrift_abs, vEdotgradvEdrift_abs
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
      tloc                    = particle(ipart)%self%time
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

      call derivs_gca(particle(ipart)%self%time,y,dydt)

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
        if (ps(igrid_working)%x({ic^DD},^D) .lt. y(^D)) then
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

      call cross(e,bhat,vE)

      vE(1:ndir)   = vE(1:ndir) / absb
      vEabs = sqrt(sum(vE(:)**2))
      if (relativistic) then
        kappa = 1.d0/sqrt(1.0d0 - sum(vE(:)**2)/c_norm**2)
      else
        kappa = 1.d0
      end if
      Mr = y(ndir+2); upar = y(ndir+1); m=particle(ipart)%self%m; q=particle(ipart)%self%q
      if (relativistic) then
        gamma = sqrt(1.0d0+upar**2/c_norm**2+2.0d0*Mr*absb/m/c_norm**2)*kappa
      else
        gamma = 1.d0
      end if

      particle(ipart)%self%u(3)      = gamma

      ! Time update
      particle(ipart)%self%time = particle(ipart)%self%time + dt_p

      ! Update payload
      if (.not. associated(usr_update_payload)) then
        call gca_update_payload(particle(ipart)%igrid,ps(particle(ipart)%igrid)%w,pso(particle(ipart)%igrid)%w,ps(particle(ipart)%igrid)%x, &
             particle(ipart)%self%x,particle(ipart)%self%u,q,m,payload,npayload,particle(ipart)%self%time)
      else
        call usr_update_payload(particle(ipart)%igrid,ps(particle(ipart)%igrid)%w,pso(particle(ipart)%igrid)%w,ps(particle(ipart)%igrid)%x,&
             particle(ipart)%self%x,particle(ipart)%self%u,q,m,payload,npayload,particle(ipart)%self%time)
      end if
      particle(ipart)%payload = payload

    end do

  end subroutine gca_integrate_particles

  subroutine derivs_gca_rk(t_s,y,dydt)
    use mod_global_parameters

    double precision                :: t_s, y(ndir+2)
    double precision                :: dydt(ndir+2)

    double precision,dimension(ndir):: vE, b, e, x, bhat, bdotgradb, vEdotgradb, gradkappaB
    double precision,dimension(ndir):: bdotgradvE, vEdotgradvE, u, utmp1, utmp2, utmp3
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
    call get_vec(vE_dot_grad_b, igrid_working,x,t_s,vEdotgradb)
    call get_vec(grad_kappa_B, igrid_working,x,t_s,gradkappaB)
    call get_vec(b_dot_grad_vE, igrid_working,x,t_s,bdotgradvE)
    call get_vec(vE_dot_grad_vE, igrid_working,x,t_s,vEdotgradvE)

    absb         = sqrt(sum(b(:)**2))
    bhat(1:ndir) = b(1:ndir) / absb
    epar         = sum(e(:)*bhat(:))

    call cross(e,bhat,vE)
    vE(1:ndir)   = vE(1:ndir) / absb

    if (relativistic) then
      kappa = 1.d0/sqrt(1.0d0 - sum(vE(:)**2)/c_norm**2)
      gamma = sqrt(1.0d0+upar**2/c_norm**2+2.0d0*Mr*absb/m/c_norm**2)*kappa
    else
      kappa = 1.d0
      gamma = 1.d0
    end if

    utmp1(1:ndir) = bhat(1:ndir)/(absb/kappa**2)
    utmp2(1:ndir) = Mr/(gamma*q)*gradkappaB(1:ndir) &
         + m/q* (upar**2/gamma*bdotgradb(1:ndir) + upar*vEdotgradb(1:ndir) &
                 + upar*bdotgradvE(1:ndir) + gamma*vEdotgradvE(1:ndir))
    if (relativistic) then
      utmp2(1:ndir) = utmp2(1:ndir) + upar*epar/(gamma)*vE(1:ndir)
    end if

    call cross(utmp1,utmp2,utmp3)
    u(1:ndir) = vE(1:ndir) + utmp3(1:ndir)

    ! done assembling the terms, now write rhs:
    dydt(1:ndir) = ( u(1:ndir) + upar/gamma * bhat(1:ndir) )
    dydt(ndir+1) = q/m*epar - Mr/(m*gamma) * sum(bhat(:)*gradkappaB(:)) &
                   + sum(vE(:)*(upar*bdotgradb(:)+gamma*vEdotgradb(:)))
    dydt(ndir+2) = 0.0d0 ! magnetic moment is conserved

  end subroutine derivs_gca_rk

  subroutine derivs_gca(t_s,y,dydt)
    use mod_global_parameters

    double precision                :: t_s, y(ndir+2)
    double precision                :: dydt(ndir+2)

    double precision,dimension(ndir):: vE, b, e, x, bhat, bdotgradb, vEdotgradb, gradkappaB
    double precision,dimension(ndir):: bdotgradvE, vEdotgradvE, u, utmp1, utmp2, utmp3
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
    call get_vec(vE_dot_grad_b, igrid_working,x,t_s,vEdotgradb)
    call get_vec(grad_kappa_B, igrid_working,x,t_s,gradkappaB)
    call get_vec(b_dot_grad_vE, igrid_working,x,t_s,bdotgradvE)
    call get_vec(vE_dot_grad_vE, igrid_working,x,t_s,vEdotgradvE)

    absb         = sqrt(sum(b(:)**2))
    bhat(1:ndir) = b(1:ndir) / absb

    epar         = sum(e(:)*bhat(:))
    call cross(e,bhat,vE)
    vE(1:ndir)   = vE(1:ndir) / absb

    if (relativistic) then
      kappa = sqrt(1.0d0 - sum(vE(:)**2)/c_norm**2)
      gamma = sqrt(1.0d0+upar**2/c_norm**2+2.0d0*Mr*absb/m/c_norm**2)*kappa
    else
      kappa = 1.d0
      gamma = 1.d0
    end if
    utmp1(1:ndir) = bhat(1:ndir)/(absb/kappa**2)
    utmp2(1:ndir) = Mr/(gamma*q)*gradkappaB(1:ndir) &
         + m/q* (upar**2/gamma*bdotgradb(1:ndir) + upar*vEdotgradb(1:ndir) &
                 + upar*bdotgradvE(1:ndir) + gamma*vEdotgradvE(1:ndir))
    if (relativistic) then
      utmp2(1:ndir) = utmp2(1:ndir) + upar*epar/(gamma)*vE(1:ndir)
    end if

    call cross(utmp1,utmp2,utmp3)
    u(1:ndir) = vE(1:ndir) + utmp3(1:ndir)

    ! done assembling the terms, now write rhs:
    dydt(1:ndir) = ( u(1:ndir) + upar/gamma * bhat(1:ndir) )
    dydt(ndir+1) = q/m*epar - Mr/(m*gamma) * sum(bhat(:)*gradkappaB(:)) &
                   + sum(vE(:)*(upar*bdotgradb(:)+gamma*vEdotgradb(:)))
    dydt(ndir+2) = 0.0d0 ! magnetic moment is conserved

  end subroutine derivs_gca

  !> Update payload subroutine
  subroutine gca_update_payload(igrid,w,wold,xgrid,xpart,upart,qpart,mpart,payload,npayload,particle_time)
    use mod_global_parameters
    integer, intent(in)           :: igrid,npayload
    double precision, intent(in)  :: w(ixG^T,1:nw),wold(ixG^T,1:nw)
    double precision, intent(in)  :: xgrid(ixG^T,1:ndim),xpart(1:ndir),upart(1:ndir),qpart,mpart,particle_time
    double precision, intent(out) :: payload(npayload)
    double precision, dimension(1:ndir) :: vE, e, b, bhat
    double precision, dimension(1:ndir) :: drift1, drift2
    double precision, dimension(1:ndir) :: drift3, drift4, drift5, drift6, drift7
    double precision, dimension(1:ndir) :: bdotgradb, vEdotgradb, gradkappaB
    double precision, dimension(1:ndir) :: bdotgradvE, vEdotgradvE
    double precision, dimension(1:ndir) :: gradBdrift, reldrift, bdotgradbdrift
    double precision, dimension(1:ndir) :: vEdotgradbdrift, bdotgradvEdrift
    double precision, dimension(1:ndir) :: vEdotgradvEdrift
    double precision                    :: kappa, upar, absb, gamma, vpar, vEabs
    double precision                    :: gradBdrift_abs, reldrift_abs, epar
    double precision                    :: bdotgradbdrift_abs, vEdotgradbdrift_abs
    double precision                    :: bdotgradvEdrift_abs, vEdotgradvEdrift_abs
    double precision                    :: momentumpar1, momentumpar2, momentumpar3, momentumpar4

    call get_vec(bp, igrid,xpart(1:ndir),particle_time,b)
    call get_vec(ep, igrid,xpart(1:ndir),particle_time,e)

    absb         = sqrt(sum(b(:)**2))
    bhat(1:ndir) = b(1:ndir) / absb
    epar         = sum(e(:)*bhat(:))
    call cross(e,bhat,vE)
    vE(1:ndir)   = vE(1:ndir) / absb
    vEabs = sqrt(sum(vE(:)**2))
    if (relativistic) then
      kappa = 1.d0/sqrt(1.0d0 - sum(vE(:)**2)/c_norm**2)
    else
      kappa = 1.d0
    end if
    vpar = upart(1)/upart(3)
    upar = upart(1)

    call get_vec(b_dot_grad_b, igrid,xpart(1:ndir),particle_time,bdotgradb)
    call get_vec(vE_dot_grad_b, igrid,xpart(1:ndir),particle_time,vEdotgradb)
    call get_vec(grad_kappa_B, igrid,xpart(1:ndir),particle_time,gradkappaB)
    call get_vec(b_dot_grad_vE, igrid,xpart(1:ndir),particle_time,bdotgradvE)
    call get_vec(vE_dot_grad_vE, igrid,xpart(1:ndir),particle_time,vEdotgradvE)

    drift1(1:ndir) = bhat(1:ndir)/(absb/kappa**2)
    drift2(1:ndir) = upart(2)/(upart(3)*q)*gradkappaB(1:ndir)

    call cross(drift1,drift2,gradBdrift)
    gradBdrift_abs = sqrt(sum(gradBdrift(:)**2))

    drift3(1:ndir) = upar*epar/upart(3)*vE(1:ndir)
    call cross(drift1,drift3,reldrift)
    reldrift_abs = sqrt(sum(reldrift(:)**2))

    drift4(1:ndir) = mpart/qpart* ( upar**2/upart(3)*bdotgradb(1:ndir))
    call cross(drift1,drift4,bdotgradbdrift)
    bdotgradbdrift_abs = sqrt(sum(bdotgradbdrift(:)**2))

    drift5(1:ndir) = mpart/qpart* ( upar*vEdotgradb(1:ndir))
    call cross(drift1,drift5,vEdotgradbdrift)
    vEdotgradbdrift_abs = sqrt(sum(vEdotgradbdrift(:)**2))

    drift6(1:ndir) = mpart/qpart* ( upar*bdotgradvE(1:ndir))
    call cross(drift1,drift6,bdotgradvEdrift)
    bdotgradvEdrift_abs = sqrt(sum(bdotgradvEdrift(:)**2))

    drift7(1:ndir) = mpart/qpart* (upart(3)*vEdotgradvE(1:ndir))
    call cross(drift1,drift7,vEdotgradvEdrift)
    vEdotgradvEdrift_abs = sqrt(sum(vEdotgradvEdrift(:)**2))

    momentumpar1 = qpart/mpart*epar
    momentumpar2 = -(upart(2)/m/upart(3))*sum(bhat(:)*gradkappaB(:))
    momentumpar3 = upar*sum(vE(:)*bdotgradb(:))
    momentumpar4 = upart(3)*sum(vE(:)*vEdotgradb(:))

    ! Payload update
    if (npayload > 0) then
      ! current gyroradius
      payload(1) = sqrt(2.0d0*m*upart(2)*absb)/abs(q*absb)
    end if
    if (npayload > 1) then
      ! pitch angle
      payload(2) = atan2(sqrt((2.0d0*upart(2)*absb)/(m*upart(3)**2)),vpar)
    end if
    if (npayload > 2) then
      ! particle v_perp
      payload(3) = sqrt((2.0d0*upart(2)*absb)/(m*upart(3)**2))
    end if
    if (npayload > 3) then
      ! particle parallel momentum term 1
      payload(4) = momentumpar1
    end if
    if (npayload > 4) then
      ! particle parallel momentum term 2
      payload(5) = momentumpar2
    end if
    if (npayload > 5) then
      ! particle parallel momentum term 3
      payload(6) = momentumpar3
    end if
    if (npayload > 6) then
      ! particle parallel momentum term 4
      payload(7) = momentumpar4
    end if
    if (npayload > 7) then
      ! particle ExB drift
      payload(8) = vEabs
    end if
    if (npayload > 8) then
      ! relativistic drift
      payload(9) = gradBdrift_abs
    end if
    if (npayload > 9) then
      ! gradB drift
      payload(10) = reldrift_abs
    end if
    if (npayload > 10) then
      ! bdotgradb drift
      payload(11) = bdotgradbdrift_abs
    end if
    if (npayload > 11) then
      ! vEdotgradb drift
      payload(12) = vEdotgradbdrift_abs
    end if
    if (npayload > 12) then
      ! bdotgradvE drift
      payload(13) = bdotgradvEdrift_abs
    end if
    if (npayload > 13) then
      ! vEdotgradvE drift
      payload(14) = vEdotgradvEdrift_abs
    end if

  end subroutine gca_update_payload

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
    dt_tmp = (end_time - partp%self%time)
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

    call derivs_gca(partp%self%time,y,dydt)
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

    call derivs_gca_rk(partp%self%time+dt_euler,y,dydt)
    !call derivs_gca(partp%self%time+dt_euler,y,dydt)

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
    ! The standard thing, dt=sqrt(dx/a) where we compute a from d(gamma v||)/dt and d(gamma)/dt
    ! dt_ap = sqrt(abs(dxmin*unit_length*y(ndir+3)/( dydt(ndir+1) - y(ndir+1)/y(ndir+3)*dydt(ndir+3) ) ) )
    ! vp = sqrt(sum(v(1:ndir)**))
    ! gammap = sqrt(1.0d0/(1.0d0-(vp/const_c)**2))
    ! ap = const_c**2/vp*gammap**(-3)*dydt(ndir+3)
    ! dt_ap = sqrt(dxmin*unit_length/ap)

    !dt_a = bigdouble
    !if (dt_euler .gt. smalldouble) then
    !   a = sqrt(sum((v1(1:ndir)-v0(1:ndir))**2))/dt_euler
    !   dt_a = min(sqrt(dxmin/a),bigdouble)
    !end if

    !dt_p = min(dt_tmp , dt_a)
    dt_p = dt_tmp

    ! Make sure we don't advance beyond end_time
    call limit_dt_endtime(end_time - partp%self%time, dt_p)

  end function gca_get_particle_dt


end module mod_particle_gca
