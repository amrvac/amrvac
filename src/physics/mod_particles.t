module mod_particles
  use mod_global_parameters, only: name_len, std_len
  use mod_usr_methods, only: usr_init_particles, usr_update_payload
  use mod_constants
  use mod_physics
  implicit none

  !> String describing the particle physics type
  character(len=name_len) :: physics_type_particles = ""
  !> Maximum number of particles
  integer                                :: nparticleshi
  !> Maximum number of particles in one processor
  integer                                :: nparticles_per_cpu_hi
  !> Number of additional variables for a particle
  integer                                :: npayload
  !> Number of variables for grid field
  integer                                :: ngridvars
  !> Number of particles
  integer                                :: num_particles
  !> Time of particles
  double precision                       :: t_particles
  !> Time step of particles
  double precision                       :: dt_particles
  !> Time limit of particles
  double precision                       :: tmax_particles
  !> Iteration limit of particles
  integer                                :: itmax_particles
  !> Iteration interval to save info of traced particles
  integer                                :: ditsave_particles
  !> Time interval to save snapshots of all particles
  double precision                       :: dtsave_ensemble
  !> Resistivity
  double precision                       :: particles_eta
  double precision                       :: dtheta
  logical                                :: losses
  !> Identity number and total number of particles
  integer                                :: nparticles
  !> Iteration number of paritcles
  integer                                :: it_particles
  integer                                :: itsavelast_particles

  ! these two save the list of neighboring cpus:
  integer, dimension(:), allocatable,save :: ipe_neighbor
  integer                                 :: npe_neighbors

  integer                                 :: type_particle
  !> set the current igrid for the particle integrator:
  integer                                :: igrid_working
  !> set the current ipart for the particle integrator:
  integer                                :: ipart_working

  integer, parameter                      :: unitparticles=15

  !> Array of identity numbers of particles in current processor
  integer, dimension(:), allocatable      :: particles_on_mype
  !> Array of identity numbers of active particles in current processor
  integer, dimension(:), allocatable      :: particles_active_on_mype
  !> Number of particles in current processor
  integer                                 :: nparticles_on_mype
  !> Number of active particles in current processor
  integer                                 :: nparticles_active_on_mype

  !> Variable index for velocity
  integer, dimension(:), allocatable      :: vp(:)
  !> Variable index for magnetic field
  integer, dimension(:), allocatable      :: bp(:)
  !> Variable index for electric field
  integer, dimension(:), allocatable      :: ep(:)
  !> Variable index for
  integer, dimension(:), allocatable      :: grad_kappa_B(:)
  !> Variable index for
  integer, dimension(:), allocatable      :: b_dot_grad_b(:)
  !> Variable index for
  integer, dimension(:), allocatable      :: ue_dot_grad_b(:)
  !> Variable index for
  integer, dimension(:), allocatable      :: b_dot_grad_ue(:)
  !> Variable index for
  integer, dimension(:), allocatable      :: ue_dot_grad_ue(:)

  type particle_ptr
     type(particle_t), pointer         :: self
     !> extra information carried by the particle
     double precision, allocatable        :: payload(:)
     integer                              :: igrid, ipe
  end type particle_ptr

  type particle_t
     !> follow the history of the particle
     logical                        :: follow
     !> identity number
     integer                        :: index
     !> charge
     double precision               :: q
     !> mass
     double precision               :: m
     !> time
     double precision               :: t
     !> time step
     double precision               :: dt
     !> coordinates
     double precision, dimension(3) :: x
     !> velocity, momentum, or special ones
     double precision, dimension(3) :: u
  end type particle_t

  ! Array containing all particles
  type(particle_ptr), dimension(:), allocatable  :: particle

  procedure(init_particles), pointer  :: phys_init_particles => null()
  procedure(fill_gridvars), pointer   :: phys_fill_gridvars => null()
  procedure(integrate_particles), pointer   :: phys_integrate_particles => null()
  procedure(set_particles_dt), pointer   :: phys_set_particles_dt => null()

  abstract interface

    subroutine init_particles()
      use mod_Knuth_random
      use mod_global_parameters
    end subroutine init_particles

    subroutine fill_gridvars()
      use mod_global_parameters
    end subroutine fill_gridvars

    subroutine integrate_particles()
      use mod_global_parameters
    end subroutine integrate_particles

    subroutine set_particles_dt()
      use mod_global_parameters
    end subroutine set_particles_dt

  end interface

  contains

  !> Initialize the particles module
  subroutine particles_init()

    call init_particles_vars
    call init_particles_com

  end subroutine particles_init

  !> Finalize the particles module
  subroutine particles_finish()

    call finish_particles
    call finish_particles_com
    call finish_particles_vars

  end subroutine particles_finish

  subroutine finish_particles()

    call destroy_particles(nparticles_on_mype,particles_on_mype(1:nparticles_on_mype))

  end subroutine finish_particles

  subroutine finish_particles_com()
    use mod_global_parameters

    call MPI_TYPE_FREE(type_particle,ierrmpi)

  end subroutine finish_particles_com

  subroutine finish_particles_vars()
    use mod_global_parameters

    deallocate(particle)
    deallocate(ipe_neighbor)
    deallocate(particles_on_mype)
    deallocate(particles_active_on_mype)
    deallocate(gridvars)

  end subroutine finish_particles_vars

  !> Give initial values to paramters
  subroutine init_particles_vars()
    use mod_global_parameters
    integer :: nwx, idir

    physics_type_particles='advect'
    nparticleshi = 100000
    nparticles_per_cpu_hi = 100000
    num_particles     = 1000
    npayload          = 1
    dt_particles      = bigdouble
    t_particles       = 0.0d0
    tmax_particles    = bigdouble
    itmax_particles   = biginteger
    ditsave_particles = 8
    dtsave_ensemble   = bigdouble
    dtheta            = 2.0d0*dpi / 60.0d0
    particles_eta     = 0.d0
    losses            = .false.
    nparticles = 0
    it_particles = 0
    itsavelast_particles = 0
    nparticles_on_mype   = 0
    nparticles_active_on_mype   = 0

    call particles_params_read(par_files)

    allocate(particle(1:nparticleshi))
    allocate(ipe_neighbor(0:npe-1))
    allocate(particles_on_mype(nparticles_per_cpu_hi))
    allocate(particles_active_on_mype(nparticles_per_cpu_hi))

    particles_on_mype(:) = 0
    particles_active_on_mype(:) = 0

    select case(physics_type_particles)
    case('advect')
      ngridvars=ndir
      allocate(vp(ndir))
      do idir = 1, ndir
        vp(idir) = idir
      end do
      phys_init_particles => advect_init_particles
      phys_fill_gridvars => advect_fill_gridvars
      phys_integrate_particles => advect_integrate_particles
      phys_set_particles_dt => advect_set_particles_dt
      if(.not. associated(usr_update_payload)) usr_update_payload => advect_update_payload
    case('Lorentz')
      if(physics_type/='mhd') call mpistop("Lorentz particles need magnetic field!")
      if(ndir/=3) call mpistop("Lorentz particles need ndir=3!")
      dtsave_ensemble=dtsave_ensemble*unit_time
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
      phys_init_particles => Lorentz_init_particles
      phys_fill_gridvars => Lorentz_fill_gridvars
      phys_integrate_particles => Lorentz_integrate_particles
      phys_set_particles_dt => Lorentz_set_particles_dt
    case('gca')
      if(physics_type/='mhd') call mpistop("GCA particles need magnetic field!")
      if(ndir/=3) call mpistop("GCA particles need ndir=3!")
      dtsave_ensemble=dtsave_ensemble*unit_time
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
      phys_init_particles => gca_init_particles
      phys_fill_gridvars => gca_fill_gridvars
      phys_integrate_particles => gca_integrate_particles
      phys_set_particles_dt => gca_set_particles_dt
    case default
      call mpistop("unknown physics_type_particles (advect,gca,Lorentz)")
    end select

    if(associated(usr_init_particles)) phys_init_particles => usr_init_particles

  end subroutine init_particles_vars

  !> Read this module's parameters from a file
  subroutine particles_params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /particles_list/ physics_type_particles,nparticleshi, &
      nparticles_per_cpu_hi,ditsave_particles, particles_eta, &
      dtsave_ensemble,num_particles,npayload,itmax_particles,tmax_particles,dtheta,losses

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, particles_list, end=111)
111    close(unitpar)
    end do

  end subroutine particles_params_read

  !> Initialise communicators for particles
  subroutine init_particles_com()
    use mod_global_parameters

    integer :: oldtypes(0:7), offsets(0:7), blocklengths(0:7)

    oldtypes(0) = MPI_LOGICAL
    oldtypes(1) = MPI_INTEGER
    oldtypes(2:7) = MPI_DOUBLE_PRECISION

    blocklengths(0:5)=1
    blocklengths(6)=3
    blocklengths(7)=3

    offsets(0) = 0
    offsets(1) = size_logical * blocklengths(0)
    offsets(2) = offsets(1) + size_int * blocklengths(1)
    offsets(3) = offsets(2) + size_double * blocklengths(2)
    offsets(4) = offsets(3) + size_double * blocklengths(3)
    offsets(5) = offsets(4) + size_double * blocklengths(4)
    offsets(6) = offsets(5) + size_double * blocklengths(5)
    offsets(7) = offsets(6) + size_double * blocklengths(6)

    call MPI_TYPE_STRUCT(8,blocklengths,offsets,oldtypes,type_particle,ierrmpi)
    call MPI_TYPE_COMMIT(type_particle,ierrmpi)

  end subroutine init_particles_com

  subroutine advect_init_particles()
    ! initialise the particles
    use mod_Knuth_random
    use mod_global_parameters

    double precision, dimension(ndir)    :: x, v
    double precision, dimension(num_particles,ndir) :: rrd
    double precision                     :: w(ixG^T,1:nw)
    integer                              :: igrid_particle, ipe_particle
    integer                              :: idir,seed
    logical, dimension(1:num_particles)  :: follow

    follow=.false.
    x(:)=0.0d0
    ! initialise the random number generator
    seed = 310952
    call rnstrt(seed)
    do idir=1,ndir
      call rand(rrd(:,idir),num_particles)
    end do
    follow(num_particles/2)=.true.

    do while (nparticles .lt. num_particles)

   {^D&x(^D) = xprobmin^D + rrd(nparticles+1,^D) * (xprobmax^D - xprobmin^D)\}

      call find_particle_ipe(x,igrid_particle,ipe_particle)

      nparticles=nparticles+1
      particle(nparticles)%igrid  = igrid_particle
      particle(nparticles)%ipe    = ipe_particle

      if(ipe_particle == mype) then
        call push_particle_into_particles_on_mype(nparticles)
        allocate(particle(nparticles)%self)
        particle(nparticles)%self%follow = follow(nparticles)
        particle(nparticles)%self%index  = nparticles
        particle(nparticles)%self%t      = 0.0d0
        particle(nparticles)%self%dt     = 0.0d0
        particle(nparticles)%self%x = 0.d0
        particle(nparticles)%self%x(1:ndir) = x(1:ndir)

        w=pw(igrid_particle)%w
        call phys_to_primitive(ixG^LL,ixG^LL,w,pw(igrid_particle)%x)
        do idir=1,ndir
          call interpolate_var(igrid_particle,ixG^LL,ixM^LL,&
            w(ixG^T,iw_mom(idir)),pw(igrid_particle)%x,x,v(idir))
        end do
        particle(nparticles)%self%u(:) = 0.d0
        particle(nparticles)%self%u(1:ndir) = v(1:ndir)
        allocate(particle(nparticles)%payload(npayload))
        particle(nparticles)%payload=0.d0

        ! Payload update
        call usr_update_payload(igrid_particle,pw(igrid_particle)%w,pw(igrid_particle)%wold,&
          pw(igrid_particle)%x,x,particle(nparticles)%payload,npayload,0.d0)

      end if

    end do

  end subroutine advect_init_particles

  subroutine Lorentz_init_particles()
    ! initialise the particles
    use mod_Knuth_random
    use mod_global_parameters

    double precision, dimension(ndir)   :: x,B,u
    double precision, dimension(num_particles,ndir) :: rrd, srd, trd
    double precision                    :: absB, absS, theta, prob, theta2, prob2
    integer                             :: igrid_particle, ipe_particle
    integer                             :: seed, idir
    logical, dimension(1:num_particles) :: follow

    follow=.false.
    ! initialise the random number generator
    seed = 310952
    call rnstrt(seed)
    do idir=1,ndir
      call rand(rrd(:,idir),num_particles)
      call rand(srd(:,idir),num_particles)
      call rand(trd(:,idir),num_particles)
    end do

    ! first find ipe and igrid responsible for particle

    x(:)=0.0d0

    do while (nparticles .lt. num_particles)

   {^D&x(^D) = xprobmin^D + rrd(nparticles+1,^D) * (xprobmax^D - xprobmin^D)\}

      call find_particle_ipe(x,igrid_particle,ipe_particle)

      nparticles=nparticles+1
      particle(nparticles)%igrid  = igrid_particle
      particle(nparticles)%ipe    = ipe_particle

      if(ipe_particle == mype) then
        call push_particle_into_particles_on_mype(nparticles)
        allocate(particle(nparticles)%self)
        particle(nparticles)%self%follow = follow(nparticles)
        particle(nparticles)%self%index  = nparticles
        particle(nparticles)%self%q      = - const_e
        particle(nparticles)%self%m      =   const_me

        particle(nparticles)%self%t      = 0.0d0
        particle(nparticles)%self%dt     = 0.0d0
        particle(nparticles)%self%x = 0.d0
        particle(nparticles)%self%x(1:ndir) = x(1:ndir)

        ! Maxwellian velocity distribution
        absS = sqrt(sum(srd(nparticles,:)**2))

        ! Maxwellian velocity distribution assigned here
        prob  = sqrt(-2.0d0*log(1.0-0.999999*srd(nparticles,1)));
        prob2 = sqrt(-2.0d0*log(1.0-0.999999*srd(nparticles,2)));

        ! random pitch angle given to each particle
        theta  = 2.0d0*dpi*trd(nparticles,1)
        theta2 = 2.0d0*dpi*trd(nparticles,2)

        ! random maxwellian velocity
        ! momentum gamma*v/c normalised by speed of light
        u(1) =  unit_velocity *prob*dcos(theta)/const_c
        u(2) =  unit_velocity *prob*dsin(theta)/const_c
        u(3) =  unit_velocity *prob2*dcos(theta2)/const_c
        particle(nparticles)%self%u(:) = u(:)
        ! initialise payloads for Lorentz module
        allocate(particle(nparticles)%payload(npayload))
        particle(nparticles)%payload(:) = 0.0d0

      end if

    end do

  end subroutine Lorentz_init_particles

  subroutine gca_init_particles()
    ! initialise the particles
    use mod_Knuth_random
    use mod_global_parameters

    double precision, dimension(ndir)   :: x,B,u
    double precision, dimension(num_particles,ndir) :: rrd, srd, trd
    double precision                    :: absB, absS, theta, prob, lfac, gamma
    integer                             :: igrid_particle, ipe_particle
    integer                             :: seed, idir
    logical, dimension(1:num_particles) :: follow

    follow=.false.
    ! initialise the random number generator
    seed = 310952
    call rnstrt(seed)
    do idir=1,ndir
      call rand(rrd(:,idir),num_particles)
      call rand(srd(:,idir),num_particles)
      call rand(trd(:,idir),num_particles)
    end do

    ! first find ipe and igrid responsible for particle

    x(:)=0.0d0

    do while (nparticles .lt. num_particles)

   {^D&x(^D) = xprobmin^D + rrd(nparticles+1,^D) * (xprobmax^D - xprobmin^D)\}

      call find_particle_ipe(x,igrid_particle,ipe_particle)

      nparticles=nparticles+1
      particle(nparticles)%igrid  = igrid_particle
      particle(nparticles)%ipe    = ipe_particle

      if(ipe_particle == mype) then
        call push_particle_into_particles_on_mype(nparticles)
        allocate(particle(nparticles)%self)
        particle(nparticles)%self%follow = follow(nparticles)
        particle(nparticles)%self%index  = nparticles
        particle(nparticles)%self%q      = - const_e
        particle(nparticles)%self%m      =   const_me

        particle(nparticles)%self%t      = 0.0d0
        particle(nparticles)%self%dt     = 0.0d0
        particle(nparticles)%self%x = 0.d0
        particle(nparticles)%self%x(1:ndir) = x(1:ndir)

        ! Maxwellian velocity distribution
        absS = sqrt(sum(srd(nparticles,:)**2))

        ! Maxwellian velocity distribution assigned here
        prob  = sqrt(-2.0d0*log(1.0-0.999999*srd(nparticles,1)))

        ! random pitch angle given to each particle
        theta = 2.0d0*dpi*trd(nparticles,1)

        !> random Maxwellian velocity. In this way v_perp = u(2)^2+u(3)^2, v// =
        !> sqrt(u(1)^2) and vthermal = sqrt(v_perp^2 + v//^2)

        !*sqrt(const_mp/const_me)
        u(1) = sqrt(2.0d0) * unit_velocity *prob*dsin(theta)
        !*sqrt(const_mp/const_me)
        u(2) =  unit_velocity *prob*dcos(theta)
        !*sqrt(const_mp/const_me)
        u(3) =  unit_velocity *prob*dcos(theta)

        lfac = one/sqrt(one-sum(u(:)**2)/const_c**2)

        ! particles Lorentz factor
        gamma = sqrt(1.0d0 + lfac**2*sum(u(:)**2)/const_c**2)

        do idir=1,ndir
          call interpolate_var(igrid_particle,ixG^LL,ixM^LL,&
            pw(igrid_particle)%w(ixG^T,iw_mag(idir)),pw(igrid_particle)%x,x,B(idir))
        end do
        B=B*unit_magneticfield
        absB = sqrt(sum(B(:)**2))

        ! parallel momentum component (gamma v||)
        particle(nparticles)%self%u(1) = lfac * u(1)

        ! Mr: the conserved magnetic moment
        particle(nparticles)%self%u(2) = (sqrt(u(2)**2+u(3)**2)* &
               particle(nparticles)%self%m )**2 / (2.0d0* &
               particle(nparticles)%self%m * absB)

        ! Lorentz factor
        particle(nparticles)%self%u(3) = lfac

        ! initialise payloads for guiding centre module
        if(npayload>2) then
          allocate(particle(nparticles)%payload(npayload))
          particle(nparticles)%payload(:) = 0.0d0
          ! gyroradius
          particle(nparticles)%payload(1)=sqrt(2.0d0*particle(nparticles)%self%m*&
            particle(nparticles)%self%u(2)*absB)/abs(particle(nparticles)%self%q*absB)*const_c
          ! pitch angle
          particle(nparticles)%payload(2)=datan(sqrt((2.0d0*particle(nparticles)%self%u(2)*&
            absB)/(particle(nparticles)%self%m*particle(nparticles)%self%u(3)**2))/(u(1)))
          ! perpendicular velocity
          particle(nparticles)%payload(3)=sqrt((2.0d0*particle(nparticles)%self%u(2)*&
            absB)/(particle(nparticles)%self%m*particle(nparticles)%self%u(3)**2))
        end if

      end if

    end do

  end subroutine gca_init_particles

  subroutine init_gridvars()
    use mod_global_parameters

    integer :: igrid, iigrid

    do iigrid=1,igridstail; igrid=igrids(iigrid);
      allocate(gridvars(igrid)%w(ixG^T,1:ngridvars))
      if(time_advance) allocate(gridvars(igrid)%wold(ixG^T,1:ngridvars))
    end do

    call phys_fill_gridvars()

  end subroutine init_gridvars

  subroutine finish_gridvars()
    use mod_global_parameters

    integer :: iigrid, igrid

    do iigrid=1,igridstail; igrid=igrids(iigrid);
      deallocate(gridvars(igrid)%w)
      if(time_advance) deallocate(gridvars(igrid)%wold)
    end do

  end subroutine finish_gridvars

  subroutine advect_fill_gridvars()
    use mod_global_parameters

    integer                                   :: igrid, iigrid, idir
    double precision, dimension(ixG^T,1:nw)   :: w

    do iigrid=1,igridstail; igrid=igrids(iigrid);

       ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);

       gridvars(igrid)%w(ixG^T,1:ngridvars) = 0.0d0
       w(ixG^T,1:nw) = pw(igrid)%w(ixG^T,1:nw)
       call phys_to_primitive(ixG^LL,ixG^LL,w,pw(igrid)%x)
       ! fill with velocity:
       gridvars(igrid)%w(ixG^T,vp(:)) = w(ixG^T,iw_mom(:))

       if(time_advance) then
         gridvars(igrid)%wold(ixG^T,1:ngridvars) = 0.0d0
         w(ixG^T,1:nw) = pw(igrid)%wold(ixG^T,1:nw)
         call phys_to_primitive(ixG^LL,ixG^LL,w,pw(igrid)%x)
         gridvars(igrid)%wold(ixG^T,vp(:)) = w(ixG^T,iw_mom(:))
       end if

    end do

  end subroutine advect_fill_gridvars

  subroutine gca_fill_gridvars()
    use mod_global_parameters

    integer                                   :: igrid, iigrid, idir, idim
    double precision, dimension(ixG^T,1:ndir) :: beta
    double precision, dimension(ixG^T,1:nw)   :: w,wold
    double precision                          :: current(ixG^T,7-2*ndir:3)
    integer                                   :: idirmin
    double precision, dimension(ixG^T,1:ndir) :: ue, bhat
    double precision, dimension(ixG^T)        :: kappa, kappa_B, absB, tmp

    do iigrid=1,igridstail; igrid=igrids(iigrid);

       ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
       gridvars(igrid)%w(ixG^T,1:ngridvars) = 0.0d0
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
       gridvars(igrid)%w(ixG^T,bp(:)) = &
            gridvars(igrid)%w(ixG^T,bp(:)) * unit_magneticfield
       gridvars(igrid)%w(ixG^T,ep(:)) = &
            gridvars(igrid)%w(ixG^T,ep(:)) * unit_magneticfield * unit_velocity / const_c

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

  subroutine Lorentz_fill_gridvars()
    use mod_global_parameters

    integer                                   :: igrid, iigrid, idir
    double precision, dimension(ixG^T,1:ndir) :: beta
    double precision, dimension(ixG^T,1:nw)   :: w,wold
    double precision                          :: current(ixG^T,7-2*ndir:3)
    integer                                   :: idirmin

    do iigrid=1,igridstail; igrid=igrids(iigrid);

       gridvars(igrid)%w(ixG^T,1:ngridvars) = 0.0d0
       ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
       w(ixG^T,1:nw) = pw(igrid)%w(ixG^T,1:nw)
       call phys_to_primitive(ixG^LL,ixG^LL,w,pw(igrid)%x)

       ! fill with magnetic field:
       gridvars(igrid)%w(ixG^T,bp(:)) = w(ixG^T,iw_mag(:))

       ! fill with electric field
       current = zero
       call get_current(w,ixG^LL,ixG^LLIM^D^LSUB1,idirmin,current)
       gridvars(igrid)%w(ixG^T,ep(1)) = gridvars(igrid)%w(ixG^T,bp(2)) * w(ixG^T,iw_mom(3)) &
            - gridvars(igrid)%w(ixG^T,bp(3)) * w(ixG^T,iw_mom(2)) + particles_eta * current(ixG^T,1)
       gridvars(igrid)%w(ixG^T,ep(2)) = gridvars(igrid)%w(ixG^T,bp(3)) * w(ixG^T,iw_mom(1)) &
            - gridvars(igrid)%w(ixG^T,bp(1)) * w(ixG^T,iw_mom(3)) + particles_eta * current(ixG^T,2)
       gridvars(igrid)%w(ixG^T,ep(3)) = gridvars(igrid)%w(ixG^T,bp(1)) * w(ixG^T,iw_mom(2)) &
            - gridvars(igrid)%w(ixG^T,bp(2)) * w(ixG^T,iw_mom(1)) + particles_eta * current(ixG^T,3)

       ! scale to cgs units:
       gridvars(igrid)%w(ixG^T,bp(:)) = &
            gridvars(igrid)%w(ixG^T,bp(:)) * sqrt(4.0d0*dpi*unit_velocity**2 * unit_density)
       gridvars(igrid)%w(ixG^T,ep(:)) = &
            gridvars(igrid)%w(ixG^T,ep(:)) * sqrt(4.0d0*dpi*unit_velocity**2 * unit_density) * unit_velocity / const_c

       if(time_advance) then
         gridvars(igrid)%wold(ixG^T,1:ngridvars) = 0.0d0
         wold(ixG^T,1:nw) = pw(igrid)%wold(ixG^T,1:nw)
         call phys_to_primitive(ixG^LL,ixG^LL,wold,pw(igrid)%x)
         ! fill with magnetic field:
         gridvars(igrid)%wold(ixG^T,bp(:)) = wold(ixG^T,iw_mag(:))
         ! fill with electric field
         current = zero
         call get_current(wold,ixG^LL,ixG^LLIM^D^LSUB1,idirmin,current)
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

    end do

  end subroutine Lorentz_fill_gridvars

  !> Calculate idirmin and the idirmin:3 components of the common current array
  !> make sure that dxlevel(^D) is set correctly.
  subroutine get_current(w,ixI^L,ixO^L,idirmin,current)
    use mod_global_parameters

    integer :: idirmin0
    integer :: ixO^L, idirmin, ixI^L
    double precision :: w(ixI^S,1:nw)
    integer :: idir
    ! For ndir=2 only 3rd component of J can exist, ndir=1 is impossible for MHD
    double precision :: current(ixI^S,7-2*ndir:3),bvec(ixI^S,1:ndir)

    idirmin0 = 7-2*ndir

    if (B0field) then
       do idir = 1, ndir
          bvec(ixI^S,idir)=w(ixI^S,iw_mag(idir))+block%B0(ixI^S,idir,0)
       end do
    else
       do idir = 1, ndir
          bvec(ixI^S,idir)=w(ixI^S,iw_mag(idir))
       end do
    end if

    call curlvector(bvec,ixI^L,ixO^L,current,idirmin,idirmin0,ndir)

  end subroutine get_current

  subroutine handle_particles()
    use mod_timing
    use mod_global_parameters

    if(time_advance) tmax_particles = global_time + dt
    if(physics_type_particles/='advect') tmax_particles=tmax_particles*unit_time

    tpartc0 = MPI_WTIME()

    call set_neighbor_ipe

    tpartc_grid_0=MPI_WTIME()
    call init_gridvars()
    tpartc_grid = tpartc_grid + (MPI_WTIME()-tpartc_grid_0)

    tpartc_com0=MPI_WTIME()
    call comm_particles_global
    tpartc_com=tpartc_com + (MPI_WTIME()-tpartc_com0)

    ! main integration loop:
    particle_evol: do

       call select_active_particles

       call phys_set_particles_dt

       tpartc_io_0=MPI_WTIME()
       call check_particles_output
       timeio_tot=timeio_tot+(MPI_WTIME()-tpartc_io_0)
       tpartc_io=tpartc_io+(MPI_WTIME()-tpartc_io_0)

       if (exit_condition() .eqv. .true.) exit particle_evol

       tpartc_int_0=MPI_WTIME()
       call phys_integrate_particles
       tpartc_int=tpartc_int+(MPI_WTIME()-tpartc_int_0)

       tpartc_com0=MPI_WTIME()
       call comm_particles
       tpartc_com=tpartc_com + (MPI_WTIME()-tpartc_com0)

       it_particles = it_particles + 1

    end do particle_evol

    tpartc_com0=MPI_WTIME()
    call comm_particles
    tpartc_com=tpartc_com + (MPI_WTIME()-tpartc_com0)

    call finish_gridvars()

    tpartc = tpartc + (MPI_WTIME() - tpartc0)

  end subroutine handle_particles

  subroutine advect_integrate_particles()
    ! this solves dx/dt=v for particles
    use mod_odeint
    use mod_global_parameters

    double precision, dimension(1:ndir) :: v, x
    double precision                 :: tloc, tlocnew, dt_p, h1
    double precision,parameter       :: eps=1.0d-6, hmin=1.0d-8
    integer                          :: ipart, iipart, igrid
    integer                          :: nok, nbad, ierror

    do iipart=1,nparticles_active_on_mype;ipart=particles_active_on_mype(iipart);

      dt_p = particle(ipart)%self%dt
      igrid = particle(ipart)%igrid
      igrid_working = igrid
      tloc = particle(ipart)%self%t
      x(1:ndir) = particle(ipart)%self%x(1:ndir)
      tlocnew=tloc+dt_p

      ! Position update
      ! Simple forward Euler start
      !call get_vec_advect(igrid,x,tloc,v,vp(1),vp(ndir))
      !particle(ipart)%self%u(1:ndir) = v(1:ndir)
      !particle(ipart)%self%x(1:ndir) = particle(ipart)%self%x(1:ndir) &
      !     + dt_p * v(1:ndir)
      ! Simple forward Euler end

      ! Adaptive stepwidth RK4:
      h1 = dt_p/2.0d0
      call odeint(x,ndir,tloc,tlocnew,eps,h1,hmin,nok,nbad,derivs_advect,rkqs,ierror)

      if (ierror /= 0) then
        print *, "odeint returned error code", ierror
        print *, "1 means hmin too small, 2 means MAXSTP exceeded"
        print *, "Having a problem with particle", iipart
      end if

      particle(ipart)%self%x(1:ndir) = x(1:ndir)

      ! Velocity update
      call get_vec_advect(igrid,x,tlocnew,v,vp(1),vp(ndir))
      particle(ipart)%self%u(1:ndir) = v(1:ndir)

      ! Payload update
      call usr_update_payload(igrid,pw(igrid)%w,pw(igrid)%wold,pw(igrid)%x,&
        x,particle(ipart)%payload,npayload,tlocnew)

      ! Time update
      particle(ipart)%self%t = tlocnew

    end do

  end subroutine advect_integrate_particles

  !> Example of update payload with local density
  subroutine advect_update_payload(igrid,w,wold,xgrid,xpart,payload,npayload,particle_time)
    use mod_global_parameters
    integer, intent(in)           :: igrid,npayload
    double precision, intent(in)  :: w(ixG^T,1:nw),wold(ixG^T,1:nw)
    double precision, intent(in)  :: xgrid(ixG^T,1:ndim),xpart(1:ndir),particle_time
    double precision, intent(out) :: payload(npayload)
    double precision              :: rho, rho1, rho2, td

    if (.not.time_advance) then
       call interpolate_var(igrid,ixG^LL,ixM^LL,w(ixG^T,iw_rho),xgrid,xpart,rho)
    else
       call interpolate_var(igrid,ixG^LL,ixM^LL,wold(ixG^T,iw_rho),xgrid,xpart,rho1)
       call interpolate_var(igrid,ixG^LL,ixM^LL,w(ixG^T,iw_rho),xgrid,xpart,rho2)
       td = (particle_time - global_time) / dt
       rho = rho1 * (1.0d0 - td) + rho2 * td
    end if
    payload(1) = rho * w_convert_factor(1)

  end subroutine advect_update_payload

  subroutine derivs_advect(t_s,x,dxdt)
    use mod_global_parameters
    double precision :: t_s, x(ndir)
    double precision :: dxdt(ndir)

    double precision :: v(ndir)

    call get_vec_advect(igrid_working,x,t_s,v,vp(1),vp(ndir))
    dxdt(:) = v(:)

  end subroutine derivs_advect

  subroutine advect_set_particles_dt()
    use mod_global_parameters

    integer                         :: ipart, iipart, nout
    double precision                :: t_min_mype, tout, dt_particles_mype, dt_cfl
    double precision                :: v(1:ndir)

    dt_particles      = bigdouble
    dt_particles_mype = bigdouble
    t_min_mype        = bigdouble

    do iipart=1,nparticles_active_on_mype;ipart=particles_active_on_mype(iipart);

      ! make sure we step only one cell at a time:
      v(1:ndir)=abs(particle(ipart)%self%u(1:ndir))

      ! convert to angular velocity:
      if(typeaxial =='cylindrical'.and.phi_>0) v(phi_) = abs(v(phi_)/particle(ipart)%self%x(r_))

      dt_cfl = min({rnode(rpdx^D_,particle(ipart)%igrid)/v(^D)},bigdouble)

      if(typeaxial =='cylindrical'.and.phi_>0) then
        ! phi-momentum leads to radial velocity:
        if(phi_ .gt. ndim) dt_cfl = min(dt_cfl, &
          sqrt(rnode(rpdx1_,particle(ipart)%igrid)/particle(ipart)%self%x(r_)) &
           / v(phi_))
        ! limit the delta phi of the orbit (just for aesthetic reasons):
        dt_cfl = min(dt_cfl,0.1d0/v(phi_))
        ! take some care at the axis:
        dt_cfl = min(dt_cfl,(particle(ipart)%self%x(r_)+smalldouble)/v(r_))
      end if

      particle(ipart)%self%dt = dt_cfl

      ! Make sure we don't miss an output or tmax_particles:
      ! corresponding output slot:
      nout = int(particle(ipart)%self%t/dtsave_ensemble) + 1
      tout = dble(nout) * dtsave_ensemble
      if(particle(ipart)%self%t+particle(ipart)%self%dt .gt. tout) &
           particle(ipart)%self%dt = max(tout - particle(ipart)%self%t, smalldouble * tout)
      ! bring to tmax_particles:
      if (particle(ipart)%self%t+particle(ipart)%self%dt .gt. tmax_particles) &
           particle(ipart)%self%dt = max(tmax_particles - particle(ipart)%self%t, smalldouble * tmax_particles)

      dt_particles_mype = min(particle(ipart)%self%dt,dt_particles_mype)
      t_min_mype = min(t_min_mype,particle(ipart)%self%t)

    end do ! ipart loop

    ! keep track of the global minimum:
    call MPI_ALLREDUCE(dt_particles_mype,dt_particles,1,MPI_DOUBLE_PRECISION,MPI_MIN,&
         icomm,ierrmpi)

    ! keep track of the minimum particle time:
    call MPI_ALLREDUCE(t_min_mype,t_particles,1,MPI_DOUBLE_PRECISION,MPI_MIN,&
         icomm,ierrmpi)

  end subroutine advect_set_particles_dt

  subroutine gca_integrate_particles()
    use mod_odeint
    use mod_global_parameters

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

    do iipart=1,nparticles_active_on_mype;ipart=particles_active_on_mype(iipart);
      int_choice=.false.
      dt_p = particle(ipart)%self%dt
      igrid_working = particle(ipart)%igrid
      ipart_working = particle(ipart)%self%index
      tloc = particle(ipart)%self%t
      x(1:ndir) = particle(ipart)%self%x(1:ndir)

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
      call get_vec_rk(igrid_working,y(1:ndir),tloc+dt_p,b,bp(1),bp(ndir))
      call get_vec_rk(igrid_working,y(1:ndir),tloc+dt_p,e,ep(1),ep(ndir))

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

      call get_vec_rk(igrid_working,y(1:ndir),tloc+dt_p,bdotgradb,b_dot_grad_b(1),      b_dot_grad_b(ndir))
      call get_vec_rk(igrid_working,y(1:ndir),tloc+dt_p,uedotgradb,ue_dot_grad_b(1),   ue_dot_grad_b(ndir))
      call get_vec_rk(igrid_working,y(1:ndir),tloc+dt_p,gradkappaB,grad_kappa_B(1),     grad_kappa_B(ndir))
      call get_vec_rk(igrid_working,y(1:ndir),tloc+dt_p,bdotgradue,b_dot_grad_ue(1),   b_dot_grad_ue(ndir))
      call get_vec_rk(igrid_working,y(1:ndir),tloc+dt_p,uedotgradue,ue_dot_grad_ue(1),ue_dot_grad_ue(ndir))

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
      ! current gyroradius
      particle(ipart)%payload(1) = sqrt(2.0d0*m*Mr*absb)/abs(q*absb)*const_c
      ! pitch angle
      particle(ipart)%payload(2) = datan(sqrt((2.0d0*Mr*absb)/(m*gamma**2))/vpar)
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

    call get_vec_rk(igrid_working,x,t_s,b,bp(1),bp(ndir))
    call get_vec_rk(igrid_working,x,t_s,e,ep(1),ep(ndir))
    call get_vec_rk(igrid_working,x,t_s,bdotgradb,b_dot_grad_b(1),      b_dot_grad_b(ndir))
    call get_vec_rk(igrid_working,x,t_s,uedotgradb,ue_dot_grad_b(1),   ue_dot_grad_b(ndir))
    call get_vec_rk(igrid_working,x,t_s,gradkappaB,grad_kappa_B(1),     grad_kappa_B(ndir))
    call get_vec_rk(igrid_working,x,t_s,bdotgradue,b_dot_grad_ue(1),   b_dot_grad_ue(ndir))
    call get_vec_rk(igrid_working,x,t_s,uedotgradue,ue_dot_grad_ue(1),ue_dot_grad_ue(ndir))

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

    call get_vec(igrid_working,x,t_s,b,bp(1),bp(ndir))
    call get_vec(igrid_working,x,t_s,e,ep(1),ep(ndir))
    call get_vec(igrid_working,x,t_s,bdotgradb,b_dot_grad_b(1),      b_dot_grad_b(ndir))
    call get_vec(igrid_working,x,t_s,uedotgradb,ue_dot_grad_b(1),   ue_dot_grad_b(ndir))
    call get_vec(igrid_working,x,t_s,gradkappaB,grad_kappa_B(1),     grad_kappa_B(ndir))
    call get_vec(igrid_working,x,t_s,bdotgradue,b_dot_grad_ue(1),   b_dot_grad_ue(ndir))
    call get_vec(igrid_working,x,t_s,uedotgradue,ue_dot_grad_ue(1),ue_dot_grad_ue(ndir))

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

  subroutine gca_set_particles_dt()
    use mod_odeint
    use mod_global_parameters

    double precision            :: t_min_mype, tout, dt_particles_mype, dt_cfl0, dt_cfl1, dt_a
    double precision            :: dxmin, vp, a, gammap
    double precision            :: v(ndir), y(ndir+2),ytmp(ndir+2), dydt(ndir+2), v0(ndir), v1(ndir), dydt1(ndir+2)
    double precision            :: ap0, ap1, dt_cfl_ap0, dt_cfl_ap1, dt_cfl_ap
    double precision            :: dt_max_output, dt_max_time, dt_euler, dt_tmp
    ! make these particle cfl conditions more restrictive if you are interpolating out of the grid
    double precision, parameter :: cfl=0.8d0, uparcfl=0.8d0
    double precision, parameter :: uparmin=1.0d-6*const_c
    integer                     :: ipart, iipart, nout, ic^D, igrid_particle, ipe_particle, ipe
    logical                     :: BC_applied

    ! Here the timestep for the guiding centre integration is chosen
    dt_particles      = bigdouble
    dt_particles_mype = bigdouble
    t_min_mype        = bigdouble
    dt_max_output     = bigdouble
    dt_max_time       = bigdouble

    do iipart=1,nparticles_active_on_mype;ipart=particles_active_on_mype(iipart);

      igrid_working = particle(ipart)%igrid
      ipart_working = particle(ipart)%self%index
      dt_tmp = (tmax_particles - particle(ipart)%self%t)
      if(dt_tmp .le. 0.0d0) dt_tmp = smalldouble
      ! make sure we step only one cell at a time, first check CFL at current location
      ! then we make an Euler step to the new location and check the new CFL
      ! we simply take the minimum of the two timesteps.
      ! added safety factor cfl:
      dxmin  = min({rnode(rpdx^D_,particle(ipart)%igrid)},bigdouble)*cfl
      ! initial solution vector:
      y(1:ndir) = particle(ipart)%self%x(1:ndir) ! position of guiding center
      y(ndir+1) = particle(ipart)%self%u(1) ! parallel momentum component (gamma v||)
      y(ndir+2) = particle(ipart)%self%u(2) ! conserved magnetic moment
      ytmp=y
      !y(ndir+3) = particle(ipart)%self%u(3) ! Lorentz factor of guiding centre

      call derivs_gca(particle(ipart)%self%t,y,dydt)
      v0(1:ndir) = dydt(1:ndir)
      ap0        = dydt(ndir+1)

      ! guiding center velocity:
      v(1:ndir) = abs(dydt(1:ndir))
      vp = sqrt(sum(v(:)**2))

      dt_cfl0    = dxmin/vp
      dt_cfl_ap0 = uparcfl * abs(max(abs(y(ndir+1)),uparmin) / ap0)
      !dt_cfl_ap0 = min(dt_cfl_ap0, uparcfl * sqrt(abs(unit_length*dxmin/(ap0+smalldouble))) )

      ! make an Euler step with the proposed timestep:
      ! new solution vector:
      dt_euler = min(dt_tmp,dt_cfl0,dt_cfl_ap0)
      y(1:ndir+2) = y(1:ndir+2) + dt_euler * dydt(1:ndir+2)

      particle(ipart)%self%x(1:ndir) = y(1:ndir) ! position of guiding center
      particle(ipart)%self%u(1)      = y(ndir+1) ! parallel momentum component (gamma v||)
      particle(ipart)%self%u(2)      = y(ndir+2) ! conserved magnetic moment

      ! first check if the particle is outside the physical domain or in the ghost cells
      if(.not. particle_in_igrid(ipart_working,igrid_working)) then
        y(1:ndir+2) = ytmp(1:ndir+2)
      end if

      call derivs_gca_rk(particle(ipart)%self%t+dt_euler,y,dydt)
      !call derivs_gca(particle(ipart)%self%t+dt_euler,y,dydt)

      v1(1:ndir) = dydt(1:ndir)
      ap1        = dydt(ndir+1)

      ! guiding center velocity:
      v(1:ndir) = abs(dydt(1:ndir))
      vp = sqrt(sum(v(:)**2))

      dt_cfl1    = dxmin/vp
      dt_cfl_ap1 = uparcfl * abs(max(abs(y(ndir+1)),uparmin) / ap1)
      !dt_cfl_ap1 = min(dt_cfl_ap1, uparcfl * sqrt(abs(unit_length*dxmin/(ap1+smalldouble))) )

      dt_tmp = min(dt_euler, dt_cfl1, dt_cfl_ap1)

      particle(ipart)%self%x(1:ndir) = ytmp(1:ndir) ! position of guiding center
      particle(ipart)%self%u(1)      = ytmp(ndir+1) ! parallel momentum component (gamma v||)
      particle(ipart)%self%u(2)      = ytmp(ndir+2) ! conserved magnetic moment
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

      !particle(ipart)%self%dt = min(dt_tmp , dt_a)
      particle(ipart)%self%dt = dt_tmp

      ! Make sure we don't miss an output or tmax_particles:
      ! corresponding output slot:
      nout = int(particle(ipart)%self%t/dtsave_ensemble) + 1
      tout = dble(nout) * dtsave_ensemble
      if(particle(ipart)%self%t+particle(ipart)%self%dt .gt. tout) then
        dt_max_output = tout - particle(ipart)%self%t
        if(dt_max_output .le. 0.0d0) dt_max_output = smalldouble * tout
      end if

      ! bring to tmax_particles:
      if(particle(ipart)%self%t+particle(ipart)%self%dt .gt. tmax_particles) then
        dt_max_time = (tmax_particles - particle(ipart)%self%t)
        if (dt_max_time .le. 0.0d0) dt_max_time = smalldouble * tmax_particles
      end if

      particle(ipart)%self%dt = min(particle(ipart)%self%dt,dt_max_time,dt_max_output)
      dt_particles_mype = min(particle(ipart)%self%dt,dt_particles_mype)
      t_min_mype = min(t_min_mype,particle(ipart)%self%t)

    end do ! ipart loop

    ! keep track of the global minimum:
    call MPI_ALLREDUCE(dt_particles_mype,dt_particles,1,MPI_DOUBLE_PRECISION,MPI_MIN,&
         icomm,ierrmpi)

    ! keep track of the minimum particle time:
    call MPI_ALLREDUCE(t_min_mype,t_particles,1,MPI_DOUBLE_PRECISION,MPI_MIN,&
         icomm,ierrmpi)

  end subroutine gca_set_particles_dt

  !> this is the relativistic Boris scheme, a leapfrog integrator
  subroutine Lorentz_integrate_particles()
    use mod_global_parameters

    integer                           :: ipart, iipart
    double precision                  :: lfac, q, m, dt_p, cosphi, sinphi, phi1, phi2, r, re
    double precision, dimension(ndir) :: b, e, emom, uminus, t_geom, s, udash, tmp, uplus, xcart1, xcart2, ucart2, radmom

    do iipart=1,nparticles_active_on_mype;ipart=particles_active_on_mype(iipart);

      q  = particle(ipart)%self%q
      m  = particle(ipart)%self%m
      dt_p = particle(ipart)%self%dt
      call get_b(particle(ipart)%igrid,particle(ipart)%self%x,particle(ipart)%self%t,b)
      call get_e(particle(ipart)%igrid,particle(ipart)%self%x,particle(ipart)%self%t,e)

      select case(typeaxial)

      ! CARTESIAN COORDINATES
      case('slab')
        ! Momentum update
        emom = q * e * dt_p /(2.0d0 * m * const_c)
        if(losses) then
          call get_lfac(particle(ipart)%self%u,lfac)
          re = abs(q)**2 / (m * const_c**2)
          call cross(particle(ipart)%self%u,b,tmp)
          radmom = - third * re**2 * lfac &
               * ( sum((e(:)+tmp(:)/lfac)**2)  &
               -  (sum(e(:)*particle(ipart)%self%u(:))/lfac)**2 ) &
               * particle(ipart)%self%u / m / const_c * dt_p
        else
          radmom = 0.0d0
        end if

        uminus = particle(ipart)%self%u + emom + radmom

        call get_lfac(uminus,lfac)
        call get_t(b,lfac,dt_p,q,m,t_geom)
        call get_s(t_geom,s)

        call cross(uminus,t_geom,tmp)
        udash = uminus + tmp

        call cross(udash,s,tmp)
        uplus = uminus + tmp

        if(losses) then
          call cross(uplus,b,tmp)
          radmom = - third * re**2 * lfac &
               * ( sum((e(:)+tmp(:)/lfac)**2)  &
               -  (sum(e(:)*uplus(:))/lfac)**2 ) &
               * uplus / m / const_c * dt_p
        else
          radmom = 0.0d0
        end if

        particle(ipart)%self%u = uplus + emom + radmom

        ! Position update
        ! update position
        particle(ipart)%self%x(1:ndir) = particle(ipart)%self%x(1:ndir) &
             + dt_p * particle(ipart)%self%u(1:ndir)/lfac &
             * const_c / unit_length

        ! Payload update
        ! current gyroradius
        call cross(particle(ipart)%self%u,b,tmp)
        tmp = tmp / sqrt(sum(b(:)**2))
        particle(ipart)%payload(1) = sqrt(sum(tmp(:)**2)) / sqrt(sum(b(:)**2)) * m / abs(q) * 8.9875d+20

        ! e.b
        if(npayload>1) particle(ipart)%payload(2) = sum(e(:)*b(:))/sqrt(sum(b(:)**2)*sum(e(:)**2))

      case ('cylindrical')

        !  Momentum update
        emom = q * e * dt_p /(2.0d0 * m * const_c)

        if(losses) then
          call get_lfac(particle(ipart)%self%u,lfac)
          re = abs(q)**2 / (m * const_c**2)
          call cross(particle(ipart)%self%u,b,tmp)
          radmom = - third * re**2 * lfac &
               * ( sum((e(:)+tmp(:)/lfac)**2)  &
               -  (sum(e(:)*particle(ipart)%self%u(:))/lfac)**2 ) &
               * particle(ipart)%self%u / m / const_c * dt_p
        else
          radmom = 0.0d0
        end if

        uminus = particle(ipart)%self%u + emom + radmom

        call get_lfac(uminus,lfac)
        call get_t(b,lfac,dt_p,q,m,t_geom)
        call get_s(t_geom,s)

        call cross(uminus,t_geom,tmp)
        udash = uminus + tmp

        call cross(udash,s,tmp)
        uplus = uminus + tmp

        if(losses) then
          call cross(uplus,b,tmp)
          radmom = - third * re**2 * lfac &
               * ( sum((e(:)+tmp(:)/lfac)**2)  &
               -  (sum(e(:)*uplus(:))/lfac)**2 ) &
               * uplus / m / const_c * dt_p
        else
           radmom = 0.0d0
        end if

        particle(ipart)%self%u = uplus + emom + radmom
        ! Position update
        ! Get cartesian coordinates:
        phi1       = particle(ipart)%self%x(phi_)
        cosphi     = cos(phi1)
        sinphi     = sin(phi1)

        xcart1(1)  = particle(ipart)%self%x(r_) * cosphi
        xcart1(2)  = particle(ipart)%self%x(r_) * sinphi
        xcart1(3)  = particle(ipart)%self%x(z_)

        ucart2(1)   = cosphi * particle(ipart)%self%u(r_) - sinphi * particle(ipart)%self%u(phi_)
        ucart2(2)   = cosphi * particle(ipart)%self%u(phi_) + sinphi * particle(ipart)%self%u(r_)
        ucart2(3)   = particle(ipart)%self%u(z_)

        ! update position
        xcart2(1:ndir) = xcart1(1:ndir) &
             + dt_p * ucart2(1:ndir)/lfac &
             * const_c/unit_length

        ! back to cylindrical coordinates
        phi2 = atan2(xcart2(2),xcart2(1))
        if(phi2 .lt. 0.0d0) phi2 = 2.0d0*dpi + phi2
        r = sqrt(xcart2(1)**2 + xcart2(2)**2)
        particle(ipart)%self%x(r_)   = r
        particle(ipart)%self%x(phi_) = phi2
        particle(ipart)%self%x(z_)   = xcart2(3)

        ! Rotate the momentum to the new cooridnates
        ! rotate velocities
        cosphi     = cos(phi2-phi1)
        sinphi     = sin(phi2-phi1)

        tmp = particle(ipart)%self%u
        particle(ipart)%self%u(r_)   = cosphi * tmp(r_)   + sinphi * tmp(phi_)
        particle(ipart)%self%u(phi_) = cosphi * tmp(phi_) - sinphi * tmp(r_)
        particle(ipart)%self%u(z_)   = tmp(z_)

        ! Payload update
        ! current gyroradius
        call cross(particle(ipart)%self%u,b,tmp)
        tmp = tmp / sqrt(sum(b(:)**2))
        particle(ipart)%payload(1) = sqrt(sum(tmp(:)**2)) / sqrt(sum(b(:)**2)) * m / abs(q) * 8.9875d+20

      end select

      ! Time update
      particle(ipart)%self%t = particle(ipart)%self%t + dt_p

    end do ! ipart loop

    contains

    subroutine get_t(b,lfac,dt,q,m,t)
      implicit none
      double precision, dimension(ndir), intent(in)      :: b
      double precision, intent(in)                      :: lfac, dt, q, m
      double precision, dimension(ndir), intent(out)     :: t

      t = q * b * dt / (2.0d0 * lfac * m * const_c)

    end subroutine get_t

    subroutine get_s(t,s)
      implicit none
      double precision, dimension(ndir), intent(in)   :: t
      double precision, dimension(ndir), intent(out)  :: s

      s = 2.0d0 * t / (1.0d0+sum(t(:)**2))

    end subroutine get_s

  end subroutine Lorentz_integrate_particles

  subroutine Lorentz_set_particles_dt()
    use mod_global_parameters

    integer                         :: ipart, iipart, nout
    double precision,dimension(ndir):: b,v
    double precision                :: lfac,absb,dt_particles_mype,dt_cfl
    double precision                :: t_min_mype, tout
    double precision, parameter     :: cfl=0.5d0

    dt_particles      = bigdouble
    dt_particles_mype = bigdouble
    t_min_mype        = bigdouble

    do iipart=1,nparticles_active_on_mype;ipart=particles_active_on_mype(iipart);

      call get_b(particle(ipart)%igrid,particle(ipart)%self%x,particle(ipart)%self%t,b)
      absb = sqrt(sum(b(:)**2))
      call get_lfac(particle(ipart)%self%u,lfac)

      ! CFL timestep
      ! make sure we step only one cell at a time:
      v(:) = abs(const_c * particle(ipart)%self%u(:) / lfac)

      ! convert to angular velocity:
      if(typeaxial =='cylindrical'.and.phi_>0) v(phi_) = abs(v(phi_)/particle(ipart)%self%x(r_))

      dt_cfl = min({rnode(rpdx^D_,particle(ipart)%igrid)/v(^D)},bigdouble)

      if(typeaxial =='cylindrical'.and.phi_>0) then
        ! phi-momentum leads to radial velocity:
        if(phi_ .gt. ndim) dt_cfl = min(dt_cfl, &
          sqrt(rnode(rpdx1_,particle(ipart)%igrid)/particle(ipart)%self%x(r_)) &
           / v(phi_))
        ! limit the delta phi of the orbit (just for aesthetic reasons):
        dt_cfl = min(dt_cfl,0.1d0/v(phi_))
        ! take some care at the axis:
        dt_cfl = min(dt_cfl,(particle(ipart)%self%x(r_)+smalldouble)/v(r_))
      end if

      dt_cfl = dt_cfl * cfl

      ! bound by gyro-rotation:
      particle(ipart)%self%dt = &
           abs( dtheta * const_c/unit_length * particle(ipart)%self%m * lfac &
           / (particle(ipart)%self%q * absb) )

      particle(ipart)%self%dt = min(particle(ipart)%self%dt,dt_cfl)*unit_length

      ! Make sure we don't miss an output or tmax_particles:
      ! corresponding output slot:
      nout = int(particle(ipart)%self%t/dtsave_ensemble) + 1
      tout = dble(nout) * dtsave_ensemble
      if(particle(ipart)%self%t+particle(ipart)%self%dt .gt. tout) &
           particle(ipart)%self%dt = max(tout - particle(ipart)%self%t , smalldouble * tout)

      ! bring to tmax_particles:
      if(particle(ipart)%self%t+particle(ipart)%self%dt .gt. tmax_particles) &
           particle(ipart)%self%dt = max(tmax_particles - particle(ipart)%self%t , smalldouble * tmax_particles)

      dt_particles_mype = min(particle(ipart)%self%dt,dt_particles_mype)

      t_min_mype = min(t_min_mype,particle(ipart)%self%t)

    end do ! ipart loop

    ! keep track of the global minimum:
    call MPI_ALLREDUCE(dt_particles_mype,dt_particles,1,MPI_DOUBLE_PRECISION,MPI_MIN,&
         icomm,ierrmpi)

    ! keep track of the minimum particle time:
    call MPI_ALLREDUCE(t_min_mype,t_particles,1,MPI_DOUBLE_PRECISION,MPI_MIN,&
         icomm,ierrmpi)

  end subroutine Lorentz_set_particles_dt

  subroutine get_e(igrid,x,tloc,e)
    ! Get the electric field in the grid at postion x.
    ! For ideal SRMHD, we first interpolate b and u=lfac*v/c
    ! The electric field then follows from e = b x beta, where beta=u/lfac.
    ! This ensures for the resulting e that e<b and e.b=0. Interpolating on u
    ! avoids interpolation-errors leading to v>c.
    ! For (non-ideal) MHD, we directly interpolate the electric field as
    ! there is no such constraint.
    use mod_global_parameters

    integer,intent(in)                                 :: igrid
    double precision,dimension(ndir), intent(in)       :: x
    double precision, intent(in)                       :: tloc
    double precision,dimension(ndir), intent(out)      :: e
    double precision,dimension(ndir)                   :: e1, e2
    double precision                                   :: td
    integer                                            :: ic^D,idir

    if(.not.time_advance) then
      do idir=1,ndir
        call interpolate_var(igrid,ixG^LL,ixM^LL,gridvars(igrid)%w(ixG^T,ep(idir)),pw(igrid)%x(ixG^T,1:ndim),x,e(idir))
      end do
    else
      do idir=1,ndir
        call interpolate_var(igrid,ixG^LL,ixM^LL,gridvars(igrid)%wold(ixG^T,ep(idir)),pw(igrid)%x(ixG^T,1:ndim),x,e1(idir))
        call interpolate_var(igrid,ixG^LL,ixM^LL,gridvars(igrid)%w(ixG^T,ep(idir)),pw(igrid)%x(ixG^T,1:ndim),x,e2(idir))
      end do
      td = (tloc/unit_time - global_time) / dt
      e(:) = e1(:) * (1.0d0 - td) + e2(:) * td
    end if

  end subroutine get_e

  subroutine get_b(igrid,x,tloc,b)
    use mod_global_parameters

    integer,intent(in)                                 :: igrid
    double precision,dimension(ndir), intent(in)       :: x
    double precision, intent(in)                       :: tloc
    double precision,dimension(ndir), intent(out)      :: b
    double precision,dimension(ndir)                   :: b1, b2
    double precision                                   :: td
    integer                                            :: ic^D, idir

    if(.not.time_advance) then
      do idir=1,ndir
        call interpolate_var(igrid,ixG^LL,ixM^LL,gridvars(igrid)%w(ixG^T,bp(idir)),pw(igrid)%x(ixG^T,1:ndim),x,b(idir))
      end do
    else
      do idir=1,ndir
        call interpolate_var(igrid,ixG^LL,ixM^LL,gridvars(igrid)%wold(ixG^T,bp(idir)),pw(igrid)%x(ixG^T,1:ndim),x,b1(idir))
        call interpolate_var(igrid,ixG^LL,ixM^LL,gridvars(igrid)%w(ixG^T,bp(idir)),pw(igrid)%x(ixG^T,1:ndim),x,b2(idir))
      end do
      td = (tloc/unit_time - global_time) / dt
      b(:) = b1(:) * (1.0d0 - td) + b2(:) * td
    end if

  end subroutine get_b

  subroutine get_lfac(u,lfac)
    use mod_global_parameters, only: ndir
    double precision,dimension(ndir), intent(in)        :: u
    double precision, intent(out)                      :: lfac

    lfac = sqrt(1.0d0 + sum(u(:)**2))

  end subroutine get_lfac

  subroutine cross(a,b,c)
    use mod_global_parameters
    double precision, dimension(ndir), intent(in)   :: a,b
    double precision, dimension(ndir), intent(out)  :: c

    ! ndir needs to be three for this to work!!!
    if(ndir==3) then
      select case(typeaxial)
      case('slab')
        c(1) = a(2)*b(3) - a(3)*b(2)
        c(2) = a(3)*b(1) - a(1)*b(3)
        c(3) = a(1)*b(2) - a(2)*b(1)
      case ('cylindrical')
        c(r_) = a(phi_)*b(z_) - a(z_)*b(phi_)
        c(phi_) = a(z_)*b(r_) - a(r_)*b(z_)
        c(z_) = a(r_)*b(phi_) - a(phi_)*b(r_)
      case default
         call mpistop('geometry not implemented in cross(a,b,c)')
      end select
    else
      call mpistop("cross in particles needs to be run with three components!")
    end if

  end subroutine cross

  subroutine interpolate_var_rk(igrid,ixI^L,ixO^L,gf,x,xloc,gfloc)
    use mod_global_parameters

    integer, intent(in)                   :: igrid,ixI^L, ixO^L
    double precision, intent(in)          :: gf(ixI^S)
    double precision, intent(in)          :: x(ixI^S,1:ndim)
    double precision, intent(in)          :: xloc(1:ndir)
    double precision, intent(out)         :: gfloc
    double precision                      :: xd^D
    {^IFTWOD
    double precision                      :: c00, c10
    }
    {^IFTHREED
    double precision                      :: c0, c1, c00, c10, c01, c11
    }
    integer                               :: ic^D, ic1^D, ic2^D, idir
    character(len=1024)                   :: line

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
    ! for the RK4 integration in integrate_particles.t we allow interpolation in
    ! to enter the ghost cells. After his integration the particle is communicated
    ! to the neighbouring grid block.
    !if (ic1^D.lt.ixGlo^D+1 .or. ic2^D.gt.ixGhi^D-1) then
    if(ic1^D.lt.ixGlo^D .or. ic2^D.gt.ixGhi^D) then
      line = ''
      write(line,"(a)")'Trying to interpolate from out of grid!'
      write(line,"(a,a,i3.2)")trim(line),' direction: ',^D
      write(line,"(a,a,^NDes14.6)")trim(line),' position: ',xloc(1:ndim)
      write(line,"(a,a,^NDi3.2,^NDi3.2)"),trim(line),' indices:', ic1^D,ic2^D
      call mpistop(line)
    end if
    \}

    {^IFONED
    xd1 = (xloc(1)-x(ic11,1)) / (x(ic21,1) - x(ic11,1))
    gfloc  = gf(ic11) * (1.0d0 - xd1) + gf(ic21) * xd1
    }
    {^IFTWOD
    xd1 = (xloc(1)-x(ic11,ic12,1)) / (x(ic21,ic12,1) - x(ic11,ic12,1))
    xd2 = (xloc(2)-x(ic11,ic12,2)) / (x(ic11,ic22,2) - x(ic11,ic12,2))
    c00 = gf(ic11,ic12) * (1.0d0 - xd1) + gf(ic21,ic12) * xd1
    c10 = gf(ic11,ic22) * (1.0d0 - xd1) + gf(ic21,ic22) * xd1
    gfloc  = c00 * (1.0d0 - xd2) + c10 * xd2
    }
    {^IFTHREED
    xd1 = (xloc(1)-x(ic11,ic12,ic13,1)) / (x(ic21,ic12,ic13,1) - x(ic11,ic12,ic13,1))
    xd2 = (xloc(2)-x(ic11,ic12,ic13,2)) / (x(ic11,ic22,ic13,2) - x(ic11,ic12,ic13,2))
    xd3 = (xloc(3)-x(ic11,ic12,ic13,3)) / (x(ic11,ic12,ic23,3) - x(ic11,ic12,ic13,3))

    c00 = gf(ic11,ic12,ic13) * (1.0d0 - xd1) + gf(ic21,ic12,ic13) * xd1
    c10 = gf(ic11,ic22,ic13) * (1.0d0 - xd1) + gf(ic21,ic22,ic13) * xd1
    c01 = gf(ic11,ic12,ic23) * (1.0d0 - xd1) + gf(ic21,ic12,ic23) * xd1
    c11 = gf(ic11,ic22,ic23) * (1.0d0 - xd1) + gf(ic21,ic22,ic23) * xd1

    c0  = c00 * (1.0d0 - xd2) + c10 * xd2
    c1  = c01 * (1.0d0 - xd2) + c11 * xd2

    gfloc = c0 * (1.0d0 - xd3) + c1 * xd3
    }

  end subroutine interpolate_var_rk

  subroutine interpolate_var(igrid,ixI^L,ixO^L,gf,x,xloc,gfloc)
    use mod_global_parameters

    integer, intent(in)                   :: igrid,ixI^L, ixO^L
    double precision, intent(in)          :: gf(ixI^S)
    double precision, intent(in)          :: x(ixI^S,1:ndim)
    double precision, intent(in)          :: xloc(1:ndir)
    double precision, intent(out)         :: gfloc
    double precision                      :: xd^D
    {^IFTWOD
    double precision                      :: c00, c10
    }
    {^IFTHREED
    double precision                      :: c0, c1, c00, c10, c01, c11
    }
    integer                               :: ic^D, ic1^D, ic2^D, idir
    character(len=1024)                   :: line

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
      line = ''
      write(line,"(a)")'Trying to interpolate from out of grid!'
      write(line,"(a,a,i3.2)")trim(line),' direction: ',^D
      write(line,"(a,a,^NDes14.6)")trim(line),' position: ',xloc(1:ndim)
      write(line,"(a,a,^NDi3.2,^NDi3.2)"),trim(line),' indices:', ic1^D,ic2^D
      call mpistop(line)
    end if
    \}

    {^IFONED
    xd1 = (xloc(1)-x(ic11,1)) / (x(ic21,1) - x(ic11,1))
    gfloc  = gf(ic11) * (1.0d0 - xd1) + gf(ic21) * xd1
    }
    {^IFTWOD
    xd1 = (xloc(1)-x(ic11,ic12,1)) / (x(ic21,ic12,1) - x(ic11,ic12,1))
    xd2 = (xloc(2)-x(ic11,ic12,2)) / (x(ic11,ic22,2) - x(ic11,ic12,2))
    c00 = gf(ic11,ic12) * (1.0d0 - xd1) + gf(ic21,ic12) * xd1
    c10 = gf(ic11,ic22) * (1.0d0 - xd1) + gf(ic21,ic22) * xd1
    gfloc  = c00 * (1.0d0 - xd2) + c10 * xd2
    }
    {^IFTHREED
    xd1 = (xloc(1)-x(ic11,ic12,ic13,1)) / (x(ic21,ic12,ic13,1) - x(ic11,ic12,ic13,1))
    xd2 = (xloc(2)-x(ic11,ic12,ic13,2)) / (x(ic11,ic22,ic13,2) - x(ic11,ic12,ic13,2))
    xd3 = (xloc(3)-x(ic11,ic12,ic13,3)) / (x(ic11,ic12,ic23,3) - x(ic11,ic12,ic13,3))

    c00 = gf(ic11,ic12,ic13) * (1.0d0 - xd1) + gf(ic21,ic12,ic13) * xd1
    c10 = gf(ic11,ic22,ic13) * (1.0d0 - xd1) + gf(ic21,ic22,ic13) * xd1
    c01 = gf(ic11,ic12,ic23) * (1.0d0 - xd1) + gf(ic21,ic12,ic23) * xd1
    c11 = gf(ic11,ic22,ic23) * (1.0d0 - xd1) + gf(ic21,ic22,ic23) * xd1

    c0  = c00 * (1.0d0 - xd2) + c10 * xd2
    c1  = c01 * (1.0d0 - xd2) + c11 * xd2

    gfloc = c0 * (1.0d0 - xd3) + c1 * xd3
    }

  end subroutine interpolate_var

  subroutine get_vec_rk(igrid,x,tloc,var,ibeg,iend)
    use mod_global_parameters

    integer,intent(in)                                   :: igrid, ibeg, iend
    double precision,dimension(ndir), intent(in)         :: x
    double precision, intent(in)                         :: tloc
    double precision,dimension(iend-ibeg+1), intent(out) :: var
    double precision,dimension(iend-ibeg+1)              :: e1, e2
    integer                                              :: ivar, iloc
    double precision                                     :: td

    if(.not.time_advance) then
      do ivar=ibeg,iend
        iloc = ivar-ibeg+1
        call interpolate_var_rk(igrid,ixG^LL,ixM^LL,gridvars(igrid)%w(ixG^T,ivar),pw(igrid)%x(ixG^T,1:ndim),x,var(iloc))
      end do
    else
      td = (tloc/unit_time - global_time) / dt
      do ivar=ibeg,iend
        iloc = ivar-ibeg+1
        call interpolate_var_rk(igrid,ixG^LL,ixM^LL,gridvars(igrid)%wold(ixG^T,ivar),pw(igrid)%x(ixG^T,1:ndim),x,e1(iloc))
        call interpolate_var_rk(igrid,ixG^LL,ixM^LL,gridvars(igrid)%w(ixG^T,ivar),pw(igrid)%x(ixG^T,1:ndim),x,e2(iloc))
        var(iloc) = e1(iloc) * (1.0d0 - td) + e2(iloc) * td
      end do
    end if

  end subroutine get_vec_rk

  subroutine get_vec(igrid,x,tloc,var,ibeg,iend)
    use mod_global_parameters

    integer,intent(in)                                   :: igrid, ibeg, iend
    double precision,dimension(ndir), intent(in)         :: x
    double precision, intent(in)                         :: tloc
    double precision,dimension(iend-ibeg+1), intent(out) :: var
    double precision,dimension(iend-ibeg+1)              :: e1, e2
    integer                                              :: ivar, iloc
    double precision                                     :: td

    if(.not.time_advance) then
      do ivar=ibeg,iend
        iloc = ivar-ibeg+1
        call interpolate_var(igrid,ixG^LL,ixM^LL,gridvars(igrid)%w(ixG^T,ivar),pw(igrid)%x(ixG^T,1:ndim),x,var(iloc))
      end do
    else
      td = (tloc/unit_time - global_time) / dt
      do ivar=ibeg,iend
        iloc = ivar-ibeg+1
        call interpolate_var(igrid,ixG^LL,ixM^LL,gridvars(igrid)%wold(ixG^T,ivar),pw(igrid)%x(ixG^T,1:ndim),x,e1(iloc))
        call interpolate_var(igrid,ixG^LL,ixM^LL,gridvars(igrid)%w(ixG^T,ivar),pw(igrid)%x(ixG^T,1:ndim),x,e2(iloc))
        var(iloc) = e1(iloc) * (1.0d0 - td) + e2(iloc) * td
      end do
    end if

  end subroutine get_vec

  subroutine get_vec_advect(igrid,x,tloc,var,ibeg,iend)
    use mod_global_parameters

    integer,intent(in)                                   :: igrid, ibeg, iend
    double precision,dimension(ndir), intent(in)         :: x
    double precision, intent(in)                         :: tloc
    double precision,dimension(iend-ibeg+1), intent(out) :: var
    double precision,dimension(iend-ibeg+1)              :: e1, e2
    integer                                              :: ivar, iloc
    double precision                                     :: td

    if(.not.time_advance) then
      do ivar=ibeg,iend
        iloc = ivar-ibeg+1
        call interpolate_var(igrid,ixG^LL,ixM^LL,gridvars(igrid)%w(ixG^T,ivar),pw(igrid)%x(ixG^T,1:ndim),x,var(iloc))
      end do
    else
      td = (tloc - global_time) / dt
      do ivar=ibeg,iend
        iloc = ivar-ibeg+1
        call interpolate_var(igrid,ixG^LL,ixM^LL,gridvars(igrid)%wold(ixG^T,ivar),pw(igrid)%x(ixG^T,1:ndim),x,e1(iloc))
        call interpolate_var(igrid,ixG^LL,ixM^LL,gridvars(igrid)%w(ixG^T,ivar),pw(igrid)%x(ixG^T,1:ndim),x,e2(iloc))
        var(iloc) = e1(iloc) * (1.0d0 - td) + e2(iloc) * td
      end do
    end if

  end subroutine get_vec_advect

  !> Check if we should go out of the integration loop
  logical function exit_condition()

    exit_condition = (&
         t_particles >= tmax_particles &
         .or. nparticles == 0 &
         .or. it_particles >= itmax_particles &
         )

  end function exit_condition

  subroutine time_spent_on_particles()
    use mod_timing
    use mod_global_parameters

    double precision         :: tpartc_avg, tpartc_int_avg, tpartc_io_avg, tpartc_com_avg, tpartc_grid_avg

    call MPI_REDUCE(tpartc,tpartc_avg,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,icomm,ierrmpi)
    call MPI_REDUCE(tpartc_int,tpartc_int_avg,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,icomm,ierrmpi)
    call MPI_REDUCE(tpartc_io,tpartc_io_avg,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,icomm,ierrmpi)
    call MPI_REDUCE(tpartc_com,tpartc_com_avg,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,icomm,ierrmpi)
    call MPI_REDUCE(tpartc_grid,tpartc_grid_avg,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,icomm,ierrmpi)

    if( mype ==0 ) then
       tpartc_avg     = tpartc_avg/dble(npe)
       tpartc_int_avg = tpartc_int_avg/dble(npe)
       tpartc_io_avg  = tpartc_io_avg/dble(npe)
       tpartc_com_avg  = tpartc_com_avg/dble(npe)
       tpartc_grid_avg  = tpartc_grid_avg/dble(npe)
       write(*,'(a,f12.3,a)')' Particle handling took     : ',tpartc,' sec'
       write(*,'(a,f12.2,a)')'                  Percentage: ',100.0d0*tpartc/timeloop,' %'
       write(*,'(a,f12.3,a)')' Particle IO took           : ',tpartc_io_avg,' sec'
       write(*,'(a,f12.2,a)')'                  Percentage: ',100.0d0*tpartc_io_avg/timeloop,' %'
       write(*,'(a,f12.3,a)')' Particle COM took          : ',tpartc_com_avg,' sec'
       write(*,'(a,f12.2,a)')'                  Percentage: ',100.0d0*tpartc_com_avg/timeloop,' %'
       write(*,'(a,f12.3,a)')' Particle integration took  : ',tpartc_int_avg,' sec'
       write(*,'(a,f12.2,a)')'                  Percentage: ',100.0d0*tpartc_int_avg/timeloop,' %'
       write(*,'(a,f12.3,a)')' Particle init grid took    : ',tpartc_grid_avg,' sec'
       write(*,'(a,f12.2,a)')'                  Percentage: ',100.0d0*tpartc_grid_avg/timeloop,' %'
    end if

  end subroutine time_spent_on_particles

  subroutine read_particles_snapshot()
    use mod_global_parameters

    logical,save                    :: file_exists=.false.
    character(len=std_len)          :: filename
    integer                         :: mynpayload, mynparticles
    integer, dimension(0:1)         :: buff

    ! some initialisations:
    nparticles_on_mype = 0
    mynparticles       = 0
    nparticles         = 0
    it_particles       = 0

    ! open the snapshot file on the headnode
    if (mype .eq. 0) then
       write(filename,"(a,a,i4.4,a)") trim(base_filename),'_particles',snapshotini,'.dat'
       INQUIRE(FILE=filename, EXIST=file_exists)
       if (.not. file_exists) then
          write(*,*) 'File '//trim(filename)//' with particle data does not exist, calling phys_init_particles instead'
          buff(0) = -1
          buff(1) = -1
       else
          open(unit=unitparticles,file=filename,form='unformatted',action='read',status='unknown',access='stream')
          read(unitparticles) nparticles,it_particles,mynpayload
          if (mynpayload .ne. npayload) &
               call mpistop('npayload in restart file does not match npayload in mod_particles')
          buff(0) = nparticles
          buff(1) = it_particles
       end if
    end if

    if (npe>0) call MPI_BCAST(buff,2,MPI_INTEGER,0,icomm,ierrmpi)

    ! check if the particle data was found:
    if (buff(1) .eq. -1) then
       call phys_init_particles
       return
    end if

    ! particle data is there, fill variables:
    nparticles   = buff(0)
    it_particles = buff(1)

    do while (mynparticles .lt. nparticles)
       if (mype .eq. 0) then
          do while (nparticles_on_mype .lt. nparticles_per_cpu_hi &
               .and. mynparticles .lt. nparticles)
             call read_from_snapshot
             mynparticles = mynparticles + 1
          end do
       end if ! mype==0
       if (npe>0) call MPI_BCAST(mynparticles,1,MPI_INTEGER,0,icomm,ierrmpi)
       call comm_particles_global
    end do

    if (mype .eq. 0) close(unit=unitparticles)

    itsavelast_particles = it_particles

  end subroutine read_particles_snapshot

  subroutine read_from_snapshot()

    integer :: index,igrid_particle,ipe_particle

    read(unitparticles) index
    allocate(particle(index)%self)
    allocate(particle(index)%payload(1:npayload))
    particle(index)%self%index = index

    read(unitparticles) particle(index)%self%follow
    read(unitparticles) particle(index)%self%q
    read(unitparticles) particle(index)%self%m
    read(unitparticles) particle(index)%self%t
    read(unitparticles) particle(index)%self%dt
    read(unitparticles) particle(index)%self%x
    read(unitparticles) particle(index)%self%u
    read(unitparticles) particle(index)%payload(1:npayload)

    if (particle(index)%self%follow) print*, 'follow index:', index

    call find_particle_ipe(particle(index)%self%x,igrid_particle,ipe_particle)
    particle(index)%igrid = igrid_particle
    particle(index)%ipe   = ipe_particle

    call push_particle_into_particles_on_mype(index)

  end subroutine read_from_snapshot

  subroutine write_particles_snapshot()
    use mod_global_parameters

    character(len=std_len)          :: filename
    type(particle_t), dimension(nparticles_per_cpu_hi)  :: send_particles
    type(particle_t), dimension(nparticles_per_cpu_hi)  :: receive_particles
    double precision, dimension(npayload,nparticles_per_cpu_hi)  :: send_payload
    double precision, dimension(npayload,nparticles_per_cpu_hi)  :: receive_payload
    integer                         :: status(MPI_STATUS_SIZE)
    integer,dimension(0:npe-1)      :: receive_n_particles_for_output_from_ipe
    integer                         :: ipe, ipart, iipart, send_n_particles_for_output
    logical,save                    :: file_exists=.false.

    receive_n_particles_for_output_from_ipe(:) = 0

    ! open the snapshot file on the headnode
    if (mype .eq. 0) then
       write(filename,"(a,a,i4.4,a)") trim(base_filename),'_particles',snapshotnext,'.dat'
       INQUIRE(FILE=filename, EXIST=file_exists)
       if (.not. file_exists) then
          open(unit=unitparticles,file=filename,form='unformatted',status='new',access='stream')
       else
          open(unit=unitparticles,file=filename,form='unformatted',status='replace',access='stream')
       end if
       write(unitparticles) nparticles,it_particles,npayload
    end if

    if (npe==1) then
       do iipart=1,nparticles_on_mype;ipart=particles_on_mype(iipart);
          call append_to_snapshot(particle(ipart)%self,particle(ipart)%payload)
       end do
       return
    end if

    if (mype .ne. 0) then
       call MPI_SEND(nparticles_on_mype,1,MPI_INTEGER,0,mype,icomm,ierrmpi)
       ! fill the send_buffer
       send_n_particles_for_output = nparticles_on_mype
       do iipart=1,nparticles_on_mype;ipart=particles_on_mype(iipart);
          send_particles(iipart) = particle(ipart)%self
          send_payload(1:npayload,iipart) = particle(ipart)%payload(1:npayload)
       end do
    else
       ! get number of particles on other ipes
       do ipe=1,npe-1
          call MPI_RECV(receive_n_particles_for_output_from_ipe(ipe),1,MPI_INTEGER,ipe,ipe,icomm,status,ierrmpi)
       end do
    end if

    if (mype .ne. 0) then
       call MPI_SEND(send_particles,send_n_particles_for_output,type_particle,0,mype,icomm,ierrmpi)
       call MPI_SEND(send_payload,npayload*send_n_particles_for_output,MPI_DOUBLE_PRECISION,0,mype,icomm,ierrmpi)
    end if

    if (mype==0) then
       ! first output own particles (already in receive buffer)
       do iipart=1,nparticles_on_mype;ipart=particles_on_mype(iipart);
          call append_to_snapshot(particle(ipart)%self,particle(ipart)%payload)
       end do
       ! now output the particles sent from the other ipes
       do ipe=1,npe-1
          call MPI_RECV(receive_particles,receive_n_particles_for_output_from_ipe(ipe),type_particle,ipe,ipe,icomm,status,ierrmpi)
          call MPI_RECV(receive_payload,npayload*receive_n_particles_for_output_from_ipe(ipe),MPI_DOUBLE_PRECISION,ipe,ipe,icomm,status,ierrmpi)
          do ipart=1,receive_n_particles_for_output_from_ipe(ipe)
             call append_to_snapshot(receive_particles(ipart),receive_payload(1:npayload,ipart))
          end do ! ipart
       end do ! ipe
       close(unit=unitparticles)
    end if ! mype == 0

  end subroutine write_particles_snapshot

  subroutine append_to_snapshot(myparticle,mypayload)

    type(particle_t),intent(in) :: myparticle
    double precision :: mypayload(1:npayload)

    write(unitparticles) myparticle%index
    write(unitparticles) myparticle%follow
    write(unitparticles) myparticle%q
    write(unitparticles) myparticle%m
    write(unitparticles) myparticle%t
    write(unitparticles) myparticle%dt
    write(unitparticles) myparticle%x
    write(unitparticles) myparticle%u
    write(unitparticles) mypayload(1:npayload)

  end subroutine append_to_snapshot

  subroutine init_particles_output()
    use mod_global_parameters

    character(len=std_len)            :: filename
    character(len=1024)               :: line, strdata
    integer                           :: iipart, ipart, icomp

    do iipart=1,nparticles_on_mype;ipart=particles_on_mype(iipart);
      if (particle(ipart)%self%follow) then
        filename = make_particle_filename(particle(ipart)%self%t,particle(ipart)%self%index,'individual')
        open(unit=unitparticles,file=filename,status='replace')
        line=''
        do icomp=1, npayload
          write(strdata,"(a,i2.2,a)") 'payload',icomp,','
          line = trim(line)//trim(strdata)
        end do
        write(unitparticles,"(a,a,a)") 'time,dt,x1,x2,x3,u1,u2,u3,',trim(line),'ipe,iteration,index'
        close(unit=unitparticles)
      end if
    end do

  end subroutine init_particles_output

  subroutine select_active_particles
    use mod_global_parameters

    integer                                         :: ipart, iipart
    logical                                         :: activate

    nparticles_active_on_mype = 0
    do iipart=1,nparticles_on_mype;ipart=particles_on_mype(iipart);
      activate = particle(ipart)%self%t .le. tmax_particles
      if(activate) then
        nparticles_active_on_mype = nparticles_active_on_mype + 1
        particles_active_on_mype(nparticles_active_on_mype) = ipart
      end if
    end do

  end subroutine select_active_particles

  subroutine locate_particle(index,igrid_particle,ipe_particle)
    ! given the particles unique index, tell me on which cpu and igrid it is
    ! returns -1,-1 if particle was not found
    use mod_global_parameters

    integer, intent(in)                            :: index
    integer, intent(out)                           :: igrid_particle, ipe_particle
    integer                                        :: iipart,ipart,ipe_has_particle,ipe
    logical                                        :: has_particle(0:npe-1)
    integer,dimension(0:1)                         :: buff

    has_particle(:) = .false.
    do iipart=1,nparticles_on_mype;ipart=particles_on_mype(iipart);
       if (particle(ipart)%self%index == index) then
          has_particle(mype) = .true.
          exit
       end if
    end do

    if (has_particle(mype)) then
       buff(0) = particle(ipart)%igrid
       buff(1) = mype
    end if

    if (npe>0) call MPI_ALLGATHER(has_particle(mype),1,MPI_LOGICAL,has_particle,1,MPI_LOGICAL,icomm,ierrmpi)

    ipe_has_particle = -1
    do ipe=0,npe-1
       if (has_particle(ipe) .eqv. .true.) then
          ipe_has_particle = ipe
          exit
       end if
    end do

    if (ipe_has_particle .ne. -1) then
       if (npe>0) call MPI_BCAST(buff,2,MPI_INTEGER,ipe_has_particle,icomm,ierrmpi)
       igrid_particle=buff(0)
       ipe_particle = buff(1)
    else
       igrid_particle=-1
       ipe_particle = -1
    end if

  end subroutine locate_particle

  subroutine find_particle_ipe(x,igrid_particle,ipe_particle)
    use mod_forest, only: tree_node_ptr, tree_root
    use mod_global_parameters

    double precision, dimension(ndir), intent(in)   :: x
    integer, intent(out)                            :: igrid_particle, ipe_particle

    integer, dimension(ndir,nlevelshi)              :: ig
    integer                                         :: idim, ic(ndim)
    type(tree_node_ptr)                             :: branch

    ! first check if the particle is in the domain
    if (.not. particle_in_domain(x)) then
       igrid_particle = -1
       ipe_particle   = -1
       return
    end if

    ! get the index on each level
    do idim = 1, ndim
       call get_igslice(idim,x(idim),ig(idim,:))
    end do

    ! traverse the tree until leaf is found
    branch=tree_root({ig(^D,1)})
    do while (.not.branch%node%leaf)
       {ic(^D)=ig(^D,branch%node%level+1) - 2 * branch%node%ig^D +2\}
       branch%node => branch%node%child({ic(^D)})%node
    end do

    igrid_particle = branch%node%igrid
    ipe_particle   = branch%node%ipe

  end subroutine find_particle_ipe

  logical function particle_in_domain(x)
    use mod_global_parameters

    double precision, dimension(ndim), intent(in)  :: x
    integer                                        :: idim

    particle_in_domain = .true.

    do idim=1,ndim
       select case(idim)
          {case (^D)
          if (x(^D) .lt. xprobmin^D) then
             particle_in_domain = .false.
             exit
          end if
          if (x(^D) .ge. xprobmax^D) then
             particle_in_domain = .false.
             exit
          end if
          \}
       end select
    end do

  end function particle_in_domain

  !> Quick check if particle is still in igrid
  logical function particle_in_igrid(ipart,igrid)
    use mod_global_parameters
    integer, intent(in)                            :: igrid,ipart
    integer                                        :: idim

    ! first check if the igrid is still there:
    if (.not. allocated(pw(igrid)%w)) then
       particle_in_igrid = .false.
       return
    end if

    particle_in_igrid = .true.

    do idim=1,ndim
       select case(idim)
          {case (^D)
          if (particle(ipart)%self%x(^D) .lt. rnode(rpxmin^D_,igrid)) then
             particle_in_igrid = .false.
             exit
          end if
          if (particle(ipart)%self%x(^D) .ge. rnode(rpxmax^D_,igrid)) then
             particle_in_igrid = .false.
             exit
          end if
          \}
       end select
    end do

    return
  end function particle_in_igrid


  subroutine set_neighbor_ipe()
    use mod_global_parameters

    integer              :: igrid, iigrid,ipe
    integer              :: my_neighbor_type, i^D
    logical              :: ipe_is_neighbor(0:npe-1)

    ipe_is_neighbor(:) = .false.
    do iigrid=1,igridstail; igrid=igrids(iigrid);

       {do i^DB=-1,1\}
       if (i^D==0|.and.) then
          cycle
       end if
       my_neighbor_type=neighbor_type(i^D,igrid)

       select case (my_neighbor_type)
        case (1) ! boundary
           ! do nothing
        case (2) ! fine-coarse
           call ipe_fc(i^D,igrid,ipe_is_neighbor)
        case (3) ! same level
           call ipe_srl(i^D,igrid,ipe_is_neighbor)
        case (4) ! coarse-fine
           call ipe_cf(i^D,igrid,ipe_is_neighbor)
        end select

        {end do\}

    end do

    ! remove self as neighbor
    ipe_is_neighbor(mype) = .false.

    ! Now make the list of neighbors
    npe_neighbors = 0
    do ipe=0,npe-1
       if (ipe_is_neighbor(ipe)) then
          npe_neighbors = npe_neighbors + 1
          ipe_neighbor(npe_neighbors) = ipe
       end if
    end do ! ipe

  end subroutine set_neighbor_ipe

  subroutine ipe_fc(i^D,igrid,ipe_is_neighbor)
    use mod_global_parameters
    integer, intent(in) :: i^D, igrid
    logical, intent(inout) :: ipe_is_neighbor(0:npe-1)

    ipe_is_neighbor(neighbor(2,i^D,igrid)) = .true.

  end subroutine ipe_fc

  subroutine ipe_srl(i^D,igrid,ipe_is_neighbor)
    use mod_global_parameters
    integer, intent(in) :: i^D, igrid
    logical, intent(inout) :: ipe_is_neighbor(0:npe-1)

    ipe_is_neighbor(neighbor(2,i^D,igrid)) = .true.

  end subroutine ipe_srl

  subroutine ipe_cf(i^D,igrid,ipe_is_neighbor)
    use mod_global_parameters
    integer, intent(in)    :: i^D, igrid
    logical, intent(inout) :: ipe_is_neighbor(0:npe-1)
    integer                :: ic^D, inc^D

    {do ic^DB=1+int((1-i^DB)/2),2-int((1+i^DB)/2)
        inc^DB=2*i^DB+ic^DB
    \}
    ipe_is_neighbor( neighbor_child(2,inc^D,igrid) ) = .true.

    {end do\}

  end subroutine ipe_cf

  subroutine check_particles_output()
    use mod_global_parameters

    integer                         :: ipart,iipart
    type(particle_t), dimension(nparticles_per_cpu_hi)  :: send_particles
    double precision, dimension(npayload,nparticles_per_cpu_hi)  :: send_payload
    integer                         :: send_n_particles_for_output
    integer                         :: nout
    double precision                :: tout

    ! is it time for output of individual particle?
    if ((it_particles == itsavelast_particles + ditsave_particles) &
         .or.  (it_particles == 0)) call output_individual()

    ! append to ensemble files if its time to do so
    send_n_particles_for_output = 0

    do iipart=1,nparticles_active_on_mype;ipart=particles_active_on_mype(iipart);

       ! corresponding output slot:
       nout = nint(particle(ipart)%self%t/dtsave_ensemble)
       tout = dble(nout) * dtsave_ensemble
       ! should the particle be dumped?
       if (particle(ipart)%self%t .le. tout &
            .and. particle(ipart)%self%t+particle(ipart)%self%dt .gt. tout) then
          ! have to send particle to rank zero for output
          send_n_particles_for_output = send_n_particles_for_output + 1
          send_particles(send_n_particles_for_output) = particle(ipart)%self
          send_payload(1:npayload,send_n_particles_for_output) = particle(ipart)%payload(1:npayload)
       end if ! time to output?

    end do

    call output_ensemble(send_n_particles_for_output,send_particles,send_payload,'ensemble')

  end subroutine check_particles_output

  character(len=128) function make_particle_filename(tout,index,type)
    use mod_global_parameters

    character(len=*), intent(in)    :: type
    double precision, intent(in)    :: tout
    integer, intent(in)             :: index
    integer                         :: nout, mysnapshot

    if (snapshotini .ne. -1) then
       mysnapshot = snapshotini
    else
       mysnapshot = 0
    end if

    if (convert) then
      select case(type)
      case ('ensemble')
         nout = nint(tout / dtsave_ensemble)
         write(make_particle_filename,"(a,i4.4,a,i6.6,a)") trim(base_filename),mysnapshot,'_ensemble',nout,'.csv'
      case ('destroy')
         write(make_particle_filename,"(a,i4.4,a)") trim(base_filename),mysnapshot,'_destroyed.csv'
      case ('individual')
         write(make_particle_filename,"(a,i4.4,a,i6.6,a)") trim(base_filename),mysnapshot,'_particle',index,'.csv'
      case ('followed')
         nout = nint(tout / dtsave_ensemble)
         write(make_particle_filename,"(a,i4.4,a,i6.6,a)") trim(base_filename),mysnapshot,'_followed',nout,'.csv'
      end select
    else
      select case(type)
      case ('ensemble')
         nout = nint(tout / dtsave_ensemble)
         write(make_particle_filename,"(a,a,i6.6,a)") trim(base_filename),'_ensemble',nout,'.csv'
      case ('destroy')
         write(make_particle_filename,"(a,a)") trim(base_filename),'_destroyed.csv'
      case ('individual')
         write(make_particle_filename,"(a,a,i6.6,a)") trim(base_filename),'_particle',index,'.csv'
      case ('followed')
         nout = nint(tout / dtsave_ensemble)
         write(make_particle_filename,"(a,a,i6.6,a)") trim(base_filename),'_followed',nout,'.csv'
      end select
    end if

  end function make_particle_filename

  subroutine output_ensemble(send_n_particles_for_output,send_particles,send_payload,type)
    use mod_global_parameters

    integer, intent(in)             :: send_n_particles_for_output
    type(particle_t), dimension(send_n_particles_for_output), intent(in)  :: send_particles
    double precision, dimension(npayload,send_n_particles_for_output), intent(in)  :: send_payload
    character(len=*), intent(in)    :: type
    character(len=std_len)              :: filename, filename2
    type(particle_t), dimension(nparticles_per_cpu_hi)  :: receive_particles
    double precision, dimension(npayload,nparticles_per_cpu_hi) :: receive_payload
    integer                         :: status(MPI_STATUS_SIZE)
    integer,dimension(0:npe-1)      :: receive_n_particles_for_output_from_ipe
    integer                         :: ipe, ipart

    receive_n_particles_for_output_from_ipe(:) = 0

    if (npe==1) then
       do ipart=1,send_n_particles_for_output
          filename = make_particle_filename(send_particles(ipart)%t,send_particles(ipart)%index,type)
          call output_particle(send_particles(ipart),send_payload(1:npayload,ipart),mype,filename)
          if(type=='ensemble' .and. send_particles(ipart)%follow) then
            filename2 = make_particle_filename(send_particles(ipart)%t,send_particles(ipart)%index,'followed')
            call output_particle(send_particles(ipart),send_payload(1:npayload,ipart),mype,filename2)
          end if
       end do
       return
    end if

    if (mype .ne. 0) then
       call MPI_SEND(send_n_particles_for_output,1,MPI_INTEGER,0,mype,icomm,ierrmpi)
    else
       do ipe=1,npe-1
          call MPI_RECV(receive_n_particles_for_output_from_ipe(ipe),1,MPI_INTEGER,ipe,ipe,icomm,status,ierrmpi)
       end do
    end if

    if (mype .ne. 0) then
       call MPI_SEND(send_particles,send_n_particles_for_output,type_particle,0,mype,icomm,ierrmpi)
       call MPI_SEND(send_payload,npayload*send_n_particles_for_output,MPI_DOUBLE_PRECISION,0,mype,icomm,ierrmpi)
    end if

    if (mype==0) then
       do ipart=1,send_n_particles_for_output
          filename = make_particle_filename(send_particles(ipart)%t,send_particles(ipart)%index,type)
          call output_particle(send_particles(ipart),send_payload(1:npayload,ipart),0,filename)
          if(type=='ensemble' .and. send_particles(ipart)%follow) then
            filename2 = make_particle_filename(send_particles(ipart)%t,send_particles(ipart)%index,'followed')
            call output_particle(send_particles(ipart),send_payload(1:npayload,ipart),0,filename2)
          end if
       end do
       do ipe=1,npe-1
          call MPI_RECV(receive_particles,receive_n_particles_for_output_from_ipe(ipe),type_particle,ipe,ipe,icomm,status,ierrmpi)
          call MPI_RECV(receive_payload,npayload*receive_n_particles_for_output_from_ipe(ipe),MPI_DOUBLE_PRECISION,ipe,ipe,icomm,status,ierrmpi)
          do ipart=1,receive_n_particles_for_output_from_ipe(ipe)
             filename = make_particle_filename(receive_particles(ipart)%t,receive_particles(ipart)%index,type)
             call output_particle(receive_particles(ipart),receive_payload(1:npayload,ipart),ipe,filename)
             if(type=='ensemble' .and. receive_particles(ipart)%follow) then
               filename2 = make_particle_filename(receive_particles(ipart)%t,receive_particles(ipart)%index,'followed')
               call output_particle(receive_particles(ipart),receive_payload(1:npayload,ipart),ipe,filename2)
             end if
          end do ! ipart
       end do ! ipe
    end if ! mype == 0

  end subroutine output_ensemble

  subroutine output_individual()
    use mod_global_parameters
    logical,parameter               :: output_from_root=.false.
    character(len=std_len)          :: filename
    integer                         :: ipart,iipart,ipe
    integer                         :: send_n_particles_for_output
    type(particle_t), dimension(nparticles_per_cpu_hi)  :: send_particles
    type(particle_t), dimension(nparticles_per_cpu_hi)  :: receive_particles
    double precision, dimension(npayload,nparticles_per_cpu_hi)  :: send_payload
    double precision, dimension(npayload,nparticles_per_cpu_hi)  :: receive_payload
    integer                         :: status(MPI_STATUS_SIZE)
    integer,dimension(0:npe-1)      :: receive_n_particles_for_output_from_ipe

    send_n_particles_for_output = 0
    receive_n_particles_for_output_from_ipe(:) = 0

    do iipart=1,nparticles_on_mype;ipart=particles_on_mype(iipart);
       ! should the particle be dumped?
       if (particle(ipart)%self%follow) then
          ! have to send particle to rank zero for output
          send_n_particles_for_output = send_n_particles_for_output + 1
          send_particles(send_n_particles_for_output) = particle(ipart)%self
          send_payload(1:npayload,send_n_particles_for_output) = particle(ipart)%payload(1:npayload)
       end if ! follow
    end do ! ipart

    if (npe==1) then
       do ipart=1,send_n_particles_for_output
          filename = make_particle_filename(send_particles(ipart)%t,send_particles(ipart)%index,'individual')
          call output_particle(send_particles(ipart),send_payload(1:npayload,ipart),mype,filename)
       end do
       itsavelast_particles = it_particles
       return
    end if

    if (output_from_root) then
       if (mype .ne. 0) then
          call MPI_SEND(send_n_particles_for_output,1,MPI_INTEGER,0,mype,icomm,ierrmpi)
       else
          do ipe=1,npe-1
             call MPI_RECV(receive_n_particles_for_output_from_ipe(ipe),1,MPI_INTEGER,ipe,ipe,icomm,status,ierrmpi)
          end do
       end if

       if (mype .ne. 0) then
          call MPI_SEND(send_particles,send_n_particles_for_output,type_particle,0,mype,icomm,ierrmpi)
          call MPI_SEND(send_payload,npayload*send_n_particles_for_output,MPI_DOUBLE_PRECISION,0,mype,icomm,ierrmpi)
       end if

       if (mype==0) then
          do ipart=1,send_n_particles_for_output
            filename = make_particle_filename(send_particles(ipart)%t,send_particles(ipart)%index,'individual')
            call output_particle(send_particles(ipart),send_payload(1:npayload,ipart),0,filename)
          end do
          do ipe=1,npe-1
             call MPI_RECV(receive_particles,receive_n_particles_for_output_from_ipe(ipe),type_particle,ipe,ipe,icomm,status,ierrmpi)
             call MPI_RECV(receive_payload,npayload*receive_n_particles_for_output_from_ipe(ipe),MPI_DOUBLE_PRECISION,ipe,ipe,icomm,status,ierrmpi)
             do ipart=1,receive_n_particles_for_output_from_ipe(ipe)
                filename = make_particle_filename(receive_particles(ipart)%t,receive_particles(ipart)%index,'individual')
                call output_particle(receive_particles(ipart),receive_payload(1:npayload,ipart),ipe,filename)
             end do
          end do
       end if
    else
       do ipart=1,send_n_particles_for_output
          ! generate filename
          filename = make_particle_filename(send_particles(ipart)%t,send_particles(ipart)%index,'individual')
          call output_particle(send_particles(ipart),send_payload(1:npayload,ipart),mype,filename)
       end do
    end if

    itsavelast_particles = it_particles

  end subroutine output_individual

  subroutine output_particle(myparticle,payload,ipe,filename)
    use mod_global_parameters

    type(particle_t),intent(in)                 :: myparticle
    double precision, intent(in)                   :: payload(npayload)
    integer, intent(in)                            :: ipe
    character(len=std_len),intent(in)              :: filename
    logical,save                                   :: file_exists=.false.
    character(len=20)                              :: formatstring
    double precision                               :: x(3)
    character(len=1024)                            :: line, strdata
    integer                                        :: icomp
    double precision, parameter                    :: minvalue = 1.0d-99
    double precision                               :: roundoff

    ! normalize the position
    if (typeaxial == 'slab') x = myparticle%x*length_convert_factor
    if (typeaxial == 'cylindrical') then
       x(r_)   = myparticle%x(r_)*length_convert_factor
       x(z_)   = myparticle%x(z_)*length_convert_factor
       x(phi_) = myparticle%x(phi_)
    end if
    if (typeaxial == 'spherical') then
       x(:) = myparticle%x(:)
       x(1) = x(1)*length_convert_factor
    end if

    INQUIRE(FILE=filename, EXIST=file_exists)

    if (.not. file_exists) then
       open(unit=unitparticles,file=filename,status='unknown',access='append')
       line=''
       do icomp=1, npayload
          write(strdata,"(a,i2.2,a)") 'payload',icomp,','
          line = trim(line)//trim(strdata)
       end do
       write(unitparticles,"(a,a,a)") 'time,dt,x1,x2,x3,u1,u2,u3,',trim(line),'ipe,iteration,index'
    else
       open(unit=unitparticles,file=filename,status='unknown',access='append')
    end if

    ! create the formatstring:
    line = ''
    write(strdata,"(es13.6, a)")roundoff(myparticle%t,minvalue), ','
    line = trim(line)//trim(strdata)
    write(strdata,"(es13.6, a)")roundoff(myparticle%dt,minvalue), ','
    line = trim(line)//trim(strdata)
    do icomp = 1, 3
       write(strdata,"(es13.6, a)")roundoff(x(icomp),minvalue), ','
       line = trim(line)//trim(strdata)
    end do
    do icomp = 1, 3
       write(strdata,"(es13.6, a)")roundoff(myparticle%u(icomp),minvalue), ','
       line = trim(line)//trim(strdata)
    end do
    do icomp = 1, npayload
       write(strdata,"(es13.6, a)")roundoff(payload(icomp),minvalue), ','
       line = trim(line)//trim(strdata)
    end do
    write(strdata,"(i8.7, a)")ipe, ','
    line = trim(line)//trim(strdata)
    write(strdata,"(i11.10, a)")it_particles, ','
    line = trim(line)//trim(strdata)
    write(strdata,"(i8.7)")myparticle%index
    line = trim(line)//trim(strdata)
    write(unitparticles,"(a)") trim(line)
    close(unit=unitparticles)

  end subroutine output_particle


  subroutine comm_particles()
    use mod_global_parameters

    integer                         :: ipart, iipart, igrid_particle, ipe_particle, ipe, iipe
    integer                         :: index
    integer                         :: tag_send, tag_receive, send_buff, rcv_buff
    integer                         :: status(MPI_STATUS_SIZE)
    integer, dimension(0:npe-1)     :: send_n_particles_to_ipe
    integer, dimension(0:npe-1)    :: receive_n_particles_from_ipe
    type(particle_t), dimension(nparticles_per_cpu_hi)  :: send_particles
    type(particle_t), dimension(nparticles_per_cpu_hi)  :: receive_particles
    double precision, dimension(npayload,nparticles_per_cpu_hi)  :: send_payload
    double precision, dimension(npayload,nparticles_per_cpu_hi)  :: receive_payload
    integer, dimension(nparticles_per_cpu_hi,0:npe-1)  :: particle_index_to_be_sent_to_ipe
    integer, dimension(nparticles_per_cpu_hi) :: particle_index_to_be_destroyed
    integer                                   :: destroy_n_particles_mype
    logical                                   :: BC_applied

    send_n_particles_to_ipe(:)      = 0
    receive_n_particles_from_ipe(:) = 0
    destroy_n_particles_mype        = 0

    ! check if and where to send each particle, destroy if necessary
    do iipart=1,nparticles_on_mype;ipart=particles_on_mype(iipart);

       ! first check if the particle should be destroyed (user defined criterion)
       if ( .not.time_advance .and. particle(ipart)%self%t .gt. tmax_particles ) then
          destroy_n_particles_mype  = destroy_n_particles_mype + 1
          particle_index_to_be_destroyed(destroy_n_particles_mype) = ipart
          cycle
       end if

       ! is my particle still in the same igrid?
       if (.not.particle_in_igrid(ipart,particle(ipart)%igrid)) then
          call find_particle_ipe(particle(ipart)%self%x,igrid_particle,ipe_particle)

          ! destroy particle if out of domain (signalled by return value of -1)
          if (igrid_particle == -1 )then
             call apply_periodB(particle(ipart)%self,igrid_particle,ipe_particle,BC_applied)
             if (.not. BC_applied .or. igrid_particle == -1) then
                destroy_n_particles_mype  = destroy_n_particles_mype + 1
                particle_index_to_be_destroyed(destroy_n_particles_mype) = ipart
                cycle
             end if
          end if

          ! particle still present
          particle(ipart)%igrid = igrid_particle
          particle(ipart)%ipe = ipe_particle

          ! if we have more than one core, is it on another cpu?
          if (npe .gt. 1 .and. particle(ipart)%ipe .ne. mype) then
             send_n_particles_to_ipe(ipe_particle) = &
                  send_n_particles_to_ipe(ipe_particle) + 1
             particle_index_to_be_sent_to_ipe(send_n_particles_to_ipe(ipe_particle),ipe_particle) = ipart
          end if ! ipe_particle

       end if ! particle_in_grid

    end do ! ipart

    call destroy_particles(destroy_n_particles_mype,particle_index_to_be_destroyed(1:destroy_n_particles_mype))

    ! get out when only one core:
    if (npe == 1) return

    ! communicate amount of particles to be sent/received
    do iipe=1,npe_neighbors;ipe=ipe_neighbor(iipe);
       tag_send    = mype * npe + ipe
       tag_receive = ipe * npe + mype
       call MPI_SEND(send_n_particles_to_ipe(ipe),1,MPI_INTEGER,ipe,tag_send,icomm,ierrmpi)
       call MPI_RECV(receive_n_particles_from_ipe(ipe),1,MPI_INTEGER,ipe,tag_receive,icomm,status,ierrmpi)
    end do

    ! send and receive the data of the particles
    do iipe=1,npe_neighbors;ipe=ipe_neighbor(iipe);
       tag_send    = mype * npe + ipe
       tag_receive = ipe * npe + mype

       ! should i send some particles to ipe?
       if (send_n_particles_to_ipe(ipe) .gt. 0) then

          ! create the send buffer
          do ipart = 1, send_n_particles_to_ipe(ipe)
             send_particles(ipart) = particle(particle_index_to_be_sent_to_ipe(ipart,ipe))%self
             send_payload(1:npayload,ipart) = particle(particle_index_to_be_sent_to_ipe(ipart,ipe))%payload(1:npayload)
          end do ! ipart
          call MPI_SEND(send_particles,send_n_particles_to_ipe(ipe),type_particle,ipe,tag_send,icomm,ierrmpi)
          call MPI_SEND(send_payload,npayload*send_n_particles_to_ipe(ipe),MPI_DOUBLE_PRECISION,ipe,tag_send,icomm,ierrmpi)
          do ipart = 1, send_n_particles_to_ipe(ipe)
             deallocate(particle(particle_index_to_be_sent_to_ipe(ipart,ipe))%self)
             deallocate(particle(particle_index_to_be_sent_to_ipe(ipart,ipe))%payload)
             call pull_particle_from_particles_on_mype(particle_index_to_be_sent_to_ipe(ipart,ipe))
          end do ! ipart

       end if ! send .gt. 0

       ! should i receive some particles from ipe?
       if (receive_n_particles_from_ipe(ipe) .gt. 0) then

          call MPI_RECV(receive_particles,receive_n_particles_from_ipe(ipe),type_particle,ipe,tag_receive,icomm,status,ierrmpi)
          call MPI_RECV(receive_payload,npayload*receive_n_particles_from_ipe(ipe),MPI_DOUBLE_PRECISION,ipe,tag_receive,icomm,status,ierrmpi)

          do ipart = 1, receive_n_particles_from_ipe(ipe)
             index = receive_particles(ipart)%index
             allocate(particle(index)%self)
             particle(index)%self = receive_particles(ipart)
             allocate(particle(index)%payload(npayload))
             particle(index)%payload(1:npayload) = receive_payload(1:npayload,ipart)
             call push_particle_into_particles_on_mype(index)
             ! since we don't send the igrid, need to re-locate it
             call find_particle_ipe(particle(index)%self%x,igrid_particle,ipe_particle)
             particle(index)%igrid = igrid_particle
             particle(index)%ipe = ipe_particle
          end do ! ipart

       end if ! receive .gt. 0
    end do ! ipe

  end subroutine comm_particles

  subroutine comm_particles_global()
    use mod_global_parameters

    integer                         :: ipart, iipart, igrid_particle, ipe_particle, ipe, iipe
    integer                         :: index
    integer                         :: tag_send, tag_receive, send_buff, rcv_buff
    integer                         :: status(MPI_STATUS_SIZE)
    integer, dimension(0:npe-1)     :: send_n_particles_to_ipe
    integer, dimension(0:npe-1)    :: receive_n_particles_from_ipe
    type(particle_t), dimension(nparticles_per_cpu_hi)  :: send_particles
    type(particle_t), dimension(nparticles_per_cpu_hi)  :: receive_particles
    double precision, dimension(npayload,nparticles_per_cpu_hi)  :: send_payload
    double precision, dimension(npayload,nparticles_per_cpu_hi)  :: receive_payload
    integer, dimension(nparticles_per_cpu_hi,0:npe-1)  :: particle_index_to_be_sent_to_ipe
    integer, dimension(nparticles_per_cpu_hi) :: particle_index_to_be_destroyed
    integer                                   :: destroy_n_particles_mype
    logical                                   :: BC_applied

    send_n_particles_to_ipe(:)      = 0
    receive_n_particles_from_ipe(:) = 0
    destroy_n_particles_mype        = 0

    ! check if and where to send each particle, relocate all of them (in case grid has changed)
    do iipart=1,nparticles_on_mype;ipart=particles_on_mype(iipart);

       call find_particle_ipe(particle(ipart)%self%x,igrid_particle,ipe_particle)

       particle(ipart)%igrid = igrid_particle
       particle(ipart)%ipe = ipe_particle

          ! if we have more than one core, is it on another cpu?
       if (particle(ipart)%ipe .ne. mype) then
          send_n_particles_to_ipe(ipe_particle) = &
               send_n_particles_to_ipe(ipe_particle) + 1
          particle_index_to_be_sent_to_ipe(send_n_particles_to_ipe(ipe_particle),ipe_particle) = ipart
       end if ! ipe_particle

    end do ! ipart

    ! get out when only one core:
    if (npe == 1) return

    ! communicate amount of particles to be sent/received
    do ipe=0,npe-1;if (ipe .eq. mype) cycle;
       tag_send    = mype * npe + ipe
       tag_receive = ipe * npe + mype
       call MPI_SEND(send_n_particles_to_ipe(ipe),1,MPI_INTEGER,ipe,tag_send,icomm,ierrmpi)
       call MPI_RECV(receive_n_particles_from_ipe(ipe),1,MPI_INTEGER,ipe,tag_receive,icomm,status,ierrmpi)
    end do

    ! send and receive the data of the particles
    do ipe=0,npe-1;if (ipe .eq. mype) cycle;
       tag_send    = mype * npe + ipe
       tag_receive = ipe * npe + mype

       ! should i send some particles to ipe?
       if (send_n_particles_to_ipe(ipe) .gt. 0) then

          ! create the send buffer
          do ipart = 1, send_n_particles_to_ipe(ipe)
             send_particles(ipart) = particle(particle_index_to_be_sent_to_ipe(ipart,ipe))%self
             send_payload(1:npayload,ipart) = particle(particle_index_to_be_sent_to_ipe(ipart,ipe))%payload(1:npayload)
          end do ! ipart
          call MPI_SEND(send_particles,send_n_particles_to_ipe(ipe),type_particle,ipe,tag_send,icomm,ierrmpi)
          call MPI_SEND(send_payload,npayload*send_n_particles_to_ipe(ipe),MPI_DOUBLE_PRECISION,ipe,tag_send,icomm,ierrmpi)
          do ipart = 1, send_n_particles_to_ipe(ipe)
             deallocate(particle(particle_index_to_be_sent_to_ipe(ipart,ipe))%self)
             deallocate(particle(particle_index_to_be_sent_to_ipe(ipart,ipe))%payload)
             call pull_particle_from_particles_on_mype(particle_index_to_be_sent_to_ipe(ipart,ipe))
          end do ! ipart

       end if ! send .gt. 0

       ! should i receive some particles from ipe?
       if (receive_n_particles_from_ipe(ipe) .gt. 0) then

          call MPI_RECV(receive_particles,receive_n_particles_from_ipe(ipe),type_particle,ipe,tag_receive,icomm,status,ierrmpi)
          call MPI_RECV(receive_payload,npayload*receive_n_particles_from_ipe(ipe),MPI_DOUBLE_PRECISION,ipe,tag_receive,icomm,status,ierrmpi)

          do ipart = 1, receive_n_particles_from_ipe(ipe)

             index = receive_particles(ipart)%index
             allocate(particle(index)%self)
             particle(index)%self = receive_particles(ipart)
             allocate(particle(index)%payload(npayload))
             particle(index)%payload(1:npayload) = receive_payload(1:npayload,ipart)
             call push_particle_into_particles_on_mype(index)

             ! since we don't send the igrid, need to re-locate it
             call find_particle_ipe(particle(index)%self%x,igrid_particle,ipe_particle)
             particle(index)%igrid = igrid_particle
             particle(index)%ipe = ipe_particle

          end do ! ipart

       end if ! receive .gt. 0
    end do ! ipe

  end subroutine comm_particles_global

  subroutine apply_periodB(particle,igrid_particle,ipe_particle,BC_applied)
    use mod_global_parameters

    type(particle_t), intent(inout)        :: particle
    integer, intent(inout)                    :: igrid_particle, ipe_particle
    logical,intent(out)                       :: BC_applied
    integer                                   :: idim

    BC_applied = .false.
    ! get out if we don't have any periodic BC:
    if (.not. any(periodB(1:ndim))) return

    ! go through dimensions and try re-inject the particle at the other side
    do idim=1,ndim

       if (.not. periodB(idim)) cycle

       select case(idim)
          {case (^D)
          if (particle%x(^D) .lt. xprobmin^D) then
             particle%x(^D) = particle%x(^D) + (xprobmax^D - xprobmin^D)
             BC_applied = .true.
          end if
          if (particle%x(^D) .ge. xprobmax^D) then
             particle%x(^D) = particle%x(^D) - (xprobmax^D - xprobmin^D)
             BC_applied = .true.
          end if
          \}
       end select

    end do

    call find_particle_ipe(particle%x,igrid_particle,ipe_particle)

  end subroutine apply_periodB

  !> clean up destroyed particles on all cores
  subroutine destroy_particles(destroy_n_particles_mype,particle_index_to_be_destroyed)
    use mod_global_parameters

    integer, intent(in)                                   :: destroy_n_particles_mype
    integer, dimension(1:destroy_n_particles_mype), intent(in) :: particle_index_to_be_destroyed
    type(particle_t), dimension(destroy_n_particles_mype):: destroy_particles_mype
    double precision, dimension(npayload,destroy_n_particles_mype):: destroy_payload_mype
    integer                                               :: iipart,ipart,destroy_n_particles

    destroy_n_particles             = 0

    ! append the particle to list of destroyed particles
    do iipart=1,destroy_n_particles_mype;ipart=particle_index_to_be_destroyed(iipart);
       destroy_particles_mype(iipart) = particle(ipart)%self
       destroy_payload_mype(1:npayload,iipart) = particle(ipart)%payload(1:npayload)
    end do

    call output_ensemble(destroy_n_particles_mype,destroy_particles_mype,destroy_payload_mype,'destroy')

    if (npe > 1) then
       call MPI_ALLREDUCE(destroy_n_particles_mype,destroy_n_particles,1,MPI_INTEGER,MPI_SUM,icomm,ierrmpi)
    else
       destroy_n_particles = destroy_n_particles_mype
    end if

    nparticles = nparticles - destroy_n_particles

    do iipart=1,destroy_n_particles_mype;ipart=particle_index_to_be_destroyed(iipart);
       particle(ipart)%igrid = -1
       particle(ipart)%ipe   = -1
       if(time_advance) then
         write(*,*), particle(ipart)%self%t, ': particle',ipart,' has left at it ',it_particles,' on pe', mype
         write(*,*), particle(ipart)%self%t, ': particles remaining:',nparticles, '(total)', nparticles_on_mype-1, 'on pe', mype
       endif
       deallocate(particle(ipart)%self)
       deallocate(particle(ipart)%payload)
       call pull_particle_from_particles_on_mype(ipart)
    end do

  end subroutine destroy_particles

  subroutine push_particle_into_particles_on_mype(ipart)
    implicit none

    integer, intent(in)            :: ipart

    nparticles_on_mype = nparticles_on_mype + 1
    particles_on_mype(nparticles_on_mype) = ipart

  end subroutine push_particle_into_particles_on_mype

  subroutine pull_particle_from_particles_on_mype(ipart)
    implicit none

    integer, intent(in)            :: ipart
    integer                        :: i

    do i=1,nparticles_on_mype
       if (particles_on_mype(i) == ipart) then
          particles_on_mype(i) = particles_on_mype(nparticles_on_mype)
          exit
       end if
    end do

    nparticles_on_mype = nparticles_on_mype - 1

  end subroutine pull_particle_from_particles_on_mype

end module mod_particles
