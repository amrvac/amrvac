!> Module with shared functionality for all the particle movers
module mod_particle_base
  use mod_global_parameters, only: name_len, std_len
  use mod_physics
  use mod_random
  use mod_constants

  !> String describing the particle physics type
  character(len=name_len) :: physics_type_particles = ""
  !> String describing the particle integrator type
  character(len=name_len) :: integrator_type_particles = ""
  !> Header string used in CSV files
  character(len=400)      :: csv_header
  !> Format string used in CSV files
  character(len=60)       :: csv_format
  !> Maximum number of particles
  integer                 :: nparticleshi
  !> Maximum number of particles in one processor
  integer                 :: nparticles_per_cpu_hi
  !> Number of default payload variables for a particle
  integer                 :: ndefpayload
  !> Number of user-defined payload variables for a particle
  integer                 :: nusrpayload
  !> Number of total payload variables for a particle
  integer                 :: npayload
  !> Number of variables for grid field
  integer                 :: ngridvars
  !> Number of particles
  integer                 :: num_particles
  !> Time limit of particles
  double precision        :: tmax_particles
  !> Minimum time of all particles
  double precision        :: min_particle_time
  !> Time interval to save particles
  double precision        :: dtsave_particles
  !> If positive, a constant time step for the particles
  double precision        :: const_dt_particles
  !> Time to write next particle output
  double precision        :: t_next_output
  !> Whether to write individual particle output (followed particles)
  logical                 :: write_individual
  !> Whether to write ensemble output
  logical                 :: write_ensemble
  !> Whether to write particle snapshots
  logical                 :: write_snapshot
  !> Use a relativistic particle mover?
  logical                 :: relativistic
  !> Resistivity
  double precision        :: particles_eta, particles_etah
  double precision        :: dtheta
  logical                 :: losses
  !> Identity number and total number of particles
  integer                 :: nparticles
  !> Iteration number of paritcles
  integer                 :: it_particles
  !> Output count for ensembles
  integer                 :: n_output_ensemble
  !> Output count for ensembles of destroyed particles
  integer                 :: n_output_destroyed

  ! these two save the list of neighboring cpus:
  integer, allocatable :: ipe_neighbor(:)
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
  !> Normalization factor for velocity in the integrator
  double precision                        :: integrator_velocity_factor(3)
  !> Integrator to be used for particles
  integer                                 :: integrator

  !> Variable index for magnetic field
  integer, allocatable :: bp(:)
  !> Variable index for electric field
  integer, allocatable :: ep(:)
  !> Variable index for fluid velocity
  integer, allocatable :: vp(:)
  !> Variable index for current
  integer, allocatable :: jp(:)

  type particle_ptr
    type(particle_t), pointer         :: self
    !> extra information carried by the particle
    double precision, allocatable        :: payload(:)
    integer                              :: igrid, ipe
  end type particle_ptr

  type particle_t
    !> follow the history of the particle
    logical          :: follow
    !> identity number
    integer          :: index
    !> charge
    double precision :: q
    !> mass
    double precision :: m
    !> time
    double precision :: time
    !> time step
    double precision :: dt
    !> coordinates
    double precision :: x(3)
    !> velocity, momentum, or special ones
    double precision :: u(3)
  end type particle_t

  ! Array containing all particles
  type(particle_ptr), dimension(:), allocatable  :: particle

  ! The pseudo-random number generator
  type(rng_t) :: rng

  ! Pointers for user-defined stuff
  procedure(sub_fill_gridvars), pointer              :: particles_fill_gridvars => null()
  procedure(sub_define_additional_gridvars), pointer :: particles_define_additional_gridvars => null()
  procedure(sub_fill_additional_gridvars), pointer   :: particles_fill_additional_gridvars => null()
  procedure(sub_integrate), pointer                  :: particles_integrate => null()
  procedure(fun_destroy), pointer                    :: usr_destroy_particle => null()

  abstract interface

    subroutine sub_fill_gridvars
    end subroutine sub_fill_gridvars

    subroutine sub_define_additional_gridvars(ngridvars)
      integer, intent(inout) :: ngridvars
    end subroutine sub_define_additional_gridvars

    subroutine sub_fill_additional_gridvars
    end subroutine sub_fill_additional_gridvars

    subroutine sub_integrate(end_time)
      double precision, intent(in) :: end_time
    end subroutine sub_integrate

    function fun_destroy(myparticle)
      import particle_ptr
      logical                         :: fun_destroy
      type(particle_ptr), intent(in)    :: myparticle
    end function fun_destroy

  end interface

contains

  !> Finalize the particles module
  subroutine particles_finish()
    use mod_global_parameters

    call destroy_particles(nparticles_on_mype,particles_on_mype(1:nparticles_on_mype))

    ! Clean up particle type
    call MPI_TYPE_FREE(type_particle,ierrmpi)

    ! Clean up variables
    deallocate(particle)
    deallocate(ipe_neighbor)
    deallocate(particles_on_mype)
    deallocate(particles_active_on_mype)
    deallocate(gridvars)

  end subroutine particles_finish

  !> Read this module's parameters from a file
  subroutine particles_params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /particles_list/ physics_type_particles,nparticleshi, &
         nparticles_per_cpu_hi, particles_eta, particles_etah, write_individual, write_ensemble, &
         write_snapshot, dtsave_particles,num_particles,ndefpayload,nusrpayload,tmax_particles, &
         dtheta,losses, const_dt_particles, relativistic, integrator_type_particles

    do n = 1, size(files)
      open(unitpar, file=trim(files(n)), status="old")
      read(unitpar, particles_list, end=111)
111   close(unitpar)
    end do

  end subroutine particles_params_read

  !> Give initial values to parameters
  subroutine particle_base_init()
    use mod_global_parameters
    integer            :: n, idir
    integer, parameter :: i8 = selected_int_kind(18)
    integer(i8)        :: seed(2)
    character(len=20)  :: strdata

    physics_type_particles    = 'advect'
    nparticleshi              = 10000000
    nparticles_per_cpu_hi     = 1000000
    num_particles             = 1000
    ndefpayload               = 1
    nusrpayload               = 0
    tmax_particles            = bigdouble
    dtsave_particles          = bigdouble
    const_dt_particles        = -bigdouble
    write_individual          = .true.
    write_ensemble            = .true.
    write_snapshot            = .true.
    relativistic              = .true.
    t_next_output             = 0.0d0
    dtheta                    = 2.0d0*dpi / 60.0d0
    particles_eta             = mhd_eta
    particles_etah            = mhd_etah
    losses                    = .false.
    nparticles                = 0
    it_particles              = 0
    n_output_ensemble         = 0
    n_output_destroyed        = 0
    nparticles_on_mype        = 0
    nparticles_active_on_mype = 0
    integrator_velocity_factor(:) = 1.0d0
    integrator_type_particles = 'Boris'

    call particles_params_read(par_files)

    ! If sampling, ndefpayload = nw:
    if (physics_type_particles == 'sample') ndefpayload = nw
    ! Total number of payloads
    npayload = ndefpayload + nusrpayload

    ! initialise the random number generator
    seed = [310952_i8, 8948923749821_i8]
    call rng%set_seed(seed)

    allocate(particle(1:nparticleshi))
    allocate(ipe_neighbor(0:npe-1))
    allocate(particles_on_mype(nparticles_per_cpu_hi))
    allocate(particles_active_on_mype(nparticles_per_cpu_hi))

    particles_on_mype(:) = 0
    particles_active_on_mype(:) = 0

    ! Generate header for CSV files
    csv_header = ' time, dt, x1, x2, x3, u1, u2, u3,'
    ! If sampling, name the payloads like the fluid quantities
    if (physics_type_particles == 'sample') then
      do n = 1, nw
        write(strdata,"(a,a,a)") ' ', trim(prim_wnames(n)), ','
        csv_header = trim(csv_header) // trim(strdata)
      end do
    else ! Otherwise, payloads are called pl01, pl02, ...
      do n = 1, ndefpayload
        write(strdata,"(a,i2.2,a)") ' pl', n, ','
        csv_header = trim(csv_header) // trim(strdata)
      end do
    end if
    ! User-defined payloads are called usrpl01, usrpl02, ...
    if (nusrpayload > 0) then
      do n = 1, nusrpayload
        write(strdata,"(a,i2.2,a)") ' usrpl', n, ','
        csv_header = trim(csv_header) // trim(strdata)
      end do
    end if
    csv_header = trim(csv_header) // ' ipe, iteration, index'

    ! Generate format string for CSV files
    write(csv_format, '(A,I0,A,A)') '(', 8 + npayload, '(es14.6,", ")', &
         'i8.7,", ",i11.10,", ",i8.7)'

  end subroutine particle_base_init

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

  subroutine get_maxwellian_velocity(v, velocity)
    double precision, intent(out) :: v(3)
    double precision, intent(in)  :: velocity
    double precision              :: vtmp(3), vnorm

    vtmp(1) = velocity * rng%normal()
    vtmp(2) = velocity * rng%normal()
    vtmp(3) = velocity * rng%normal()
    vnorm   = norm2(vtmp)
    v       = rng%sphere(vnorm)
  end subroutine get_maxwellian_velocity

  subroutine get_uniform_pos(x)
    use mod_global_parameters
    double precision, intent(out) :: x(3)

    call rng%unif_01_vec(x)

    {^D&x(^D) = xprobmin^D + x(^D) * (xprobmax^D - xprobmin^D)\}
    x(ndim+1:) = 0.0d0
  end subroutine get_uniform_pos

  !> Initialize grid variables for particles
  subroutine init_gridvars()
    use mod_global_parameters

    integer :: igrid, iigrid

    do iigrid=1,igridstail; igrid=igrids(iigrid);
      allocate(gridvars(igrid)%w(ixG^T,1:ngridvars))
      gridvars(igrid)%w = 0.0d0

      if (time_advance) then
        allocate(gridvars(igrid)%wold(ixG^T,1:ngridvars))
        gridvars(igrid)%wold = 0.0d0
      end if
    end do

    call particles_fill_gridvars()
    if (associated(particles_define_additional_gridvars)) then
      call particles_fill_additional_gridvars()
    end if

  end subroutine init_gridvars

  !> Deallocate grid variables for particles
  subroutine finish_gridvars()
    use mod_global_parameters

    integer :: iigrid, igrid

    do iigrid=1,igridstail; igrid=igrids(iigrid);
      deallocate(gridvars(igrid)%w)
      if(time_advance) deallocate(gridvars(igrid)%wold)
    end do

  end subroutine finish_gridvars

  !> This routine fills the particle grid variables with the default for mhd, i.e. only E and B
  subroutine fill_gridvars_default
    use mod_global_parameters
    use mod_usr_methods, only: usr_particle_fields

    integer :: igrid, iigrid
    double precision :: E(ixG^T, ndir)
    double precision :: B(ixG^T, ndir)

    do iigrid=1,igridstail; igrid=igrids(iigrid);
      if (associated(usr_particle_fields)) then
        call usr_particle_fields(ps(igrid)%w, ps(igrid)%x, E, B)
        gridvars(igrid)%w(ixG^T,ep(:)) = E
        gridvars(igrid)%w(ixG^T,bp(:)) = B
      else
        call fields_from_mhd(igrid, ps(igrid)%w, gridvars(igrid)%w)
      end if

      if (time_advance) then
        if (associated(usr_particle_fields)) then
          call usr_particle_fields(pso(igrid)%w, ps(igrid)%x, E, B)
          gridvars(igrid)%wold(ixG^T,ep(:)) = E
          gridvars(igrid)%wold(ixG^T,bp(:)) = B
        else
          call fields_from_mhd(igrid, pso(igrid)%w, gridvars(igrid)%wold)
        end if
      end if
    end do
  end subroutine fill_gridvars_default

  !> Determine fields from MHD variables
  subroutine fields_from_mhd(igrid, w_mhd, w_part)
!    use mod_mhd
    use mod_global_parameters
    integer, intent(in)             :: igrid
    double precision, intent(in)    :: w_mhd(ixG^T,nw)
    double precision, intent(inout) :: w_part(ixG^T,ngridvars)
    integer                         :: idirmin
    double precision                :: current(ixG^T,7-2*ndir:3)
    double precision                :: w(ixG^T,1:nw)

    ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);

    w(ixG^T,1:nw) = w_mhd(ixG^T,1:nw)
    w_part(ixG^T,bp(1):bp(3)) = 0.0d0
    w_part(ixG^T,ep(1):ep(3)) = 0.0d0

    call phys_to_primitive(ixG^LL,ixG^LL,w,ps(igrid)%x)

    ! fill with magnetic field:
    w_part(ixG^T,bp(:)) = w(ixG^T,iw_mag(:))

    ! fill with current
    current = zero
    call particle_get_current(w,ixG^LL,ixG^LLIM^D^LSUB1,idirmin,current)
    w_part(ixG^T,jp(:)) = current(ixG^T,:)

    ! fill with electric field
    w_part(ixG^T,ep(1)) = w_part(ixG^T,bp(2)) * &
         w(ixG^T,iw_mom(3)) - w_part(ixG^T,bp(3)) * &
         w(ixG^T,iw_mom(2)) + particles_eta * current(ixG^T,1)
    w_part(ixG^T,ep(2)) = w_part(ixG^T,bp(3)) * &
         w(ixG^T,iw_mom(1)) - w_part(ixG^T,bp(1)) * &
         w(ixG^T,iw_mom(3)) + particles_eta * current(ixG^T,2)
    w_part(ixG^T,ep(3)) = w_part(ixG^T,bp(1)) * &
         w(ixG^T,iw_mom(2)) - w_part(ixG^T,bp(2)) * &
         w(ixG^T,iw_mom(1)) + particles_eta * current(ixG^T,3)

    ! Hall term
    if (particles_etah > zero) then
      w_part(ixG^T,ep(1)) = w_part(ixG^T,ep(1)) + particles_etah/w(ixG^T,iw_rho) * &
           (current(ixG^T,2) * w_part(ixG^T,bp(3)) - current(ixG^T,3) * w_part(ixG^T,bp(2)))
      w_part(ixG^T,ep(2)) = w_part(ixG^T,ep(2)) + particles_etah/w(ixG^T,iw_rho) * &
           (-current(ixG^T,1) * w_part(ixG^T,bp(3)) + current(ixG^T,3) * w_part(ixG^T,bp(1)))
      w_part(ixG^T,ep(3)) = w_part(ixG^T,ep(3)) + particles_etah/w(ixG^T,iw_rho) * &
           (current(ixG^T,1) * w_part(ixG^T,bp(2)) - current(ixG^T,2) * w_part(ixG^T,bp(1)))
    end if

  end subroutine fields_from_mhd

  !> Calculate idirmin and the idirmin:3 components of the common current array
  !> make sure that dxlevel(^D) is set correctly.
  subroutine particle_get_current(w,ixI^L,ixO^L,idirmin,current)
    use mod_global_parameters
    use mod_geometry

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

  end subroutine particle_get_current

  !> Let particles evolve in time. The routine also handles grid variables and
  !> output.
  subroutine handle_particles()
    use mod_timing
    use mod_global_parameters

    integer :: steps_taken, nparticles_left

    tpartc0 = MPI_WTIME()

    call set_neighbor_ipe

    tpartc_grid_0=MPI_WTIME()
    call init_gridvars()
    tpartc_grid = tpartc_grid + (MPI_WTIME()-tpartc_grid_0)

    tpartc_com0=MPI_WTIME()
    call comm_particles_global
    tpartc_com=tpartc_com + (MPI_WTIME()-tpartc_com0)

    if(time_advance) then
      tmax_particles = global_time + dt
    else
      tmax_particles = global_time + (time_max-global_time)
    end if

    ! main integration loop
    do
      if (tmax_particles >= t_next_output) then
        call advance_particles(t_next_output, steps_taken)
        tpartc_io_0 = MPI_WTIME()
        call write_particle_output()
        timeio_tot  = timeio_tot+(MPI_WTIME()-tpartc_io_0)
        tpartc_io   = tpartc_io+(MPI_WTIME()-tpartc_io_0)

        ! If we are using a constant time step, this prevents multiple outputs
        ! at the same time
        t_next_output = max(t_next_output, min_particle_time) + &
             dtsave_particles
      else
        call advance_particles(tmax_particles, steps_taken)
!        call write_particle_output()
        exit
      end if

      call MPI_ALLREDUCE(nparticles_on_mype, nparticles_left, 1, MPI_INTEGER, &
           MPI_SUM, icomm, ierrmpi)
      if (nparticles_left == 0 .and. convert) call mpistop('No particles left')

    end do

    call finish_gridvars()

    tpartc = tpartc + (MPI_WTIME() - tpartc0)

  end subroutine handle_particles

  !> Routine handles the particle evolution
  subroutine advance_particles(end_time, steps_taken)
    use mod_timing
    use mod_global_parameters
    double precision, intent(in) :: end_time !< Advance at most up to this time
    integer, intent(out)         :: steps_taken
    integer                      :: n_active

    steps_taken = 0

    do
      call select_active_particles(end_time)

      ! Determine total number of active particles
      call MPI_ALLREDUCE(nparticles_active_on_mype, n_active, 1, MPI_INTEGER, &
           MPI_SUM, icomm, ierrmpi)

      if (n_active == 0) exit

      tpartc_int_0=MPI_WTIME()
      call particles_integrate(end_time)
      steps_taken = steps_taken + 1
      tpartc_int=tpartc_int+(MPI_WTIME()-tpartc_int_0)

      tpartc_com0=MPI_WTIME()
      call comm_particles()
      tpartc_com=tpartc_com + (MPI_WTIME()-tpartc_com0)

      it_particles = it_particles + 1
    end do

  end subroutine advance_particles

  pure subroutine limit_dt_endtime(t_left, dt_p)
    double precision, intent(in)    :: t_left
    double precision, intent(inout) :: dt_p
    double precision                :: n_steps
    double precision, parameter     :: eps = 1d-10

    n_steps = t_left / dt_p

    if (n_steps < 1+eps) then
      ! Last step
      dt_p = t_left
    else if (n_steps < 2-eps) then
      ! Divide time left over two steps, to prevent one tiny time step
      dt_p = 0.5d0 * t_left
    end if

  end subroutine limit_dt_endtime

  ! Get the electric field in the grid at postion x.
  ! For ideal SRMHD, we first interpolate b and u=lfac*v/c
  ! The electric field then follows from e = b x beta, where beta=u/lfac.
  ! This ensures for the resulting e that e<b and e.b=0. Interpolating on u
  ! avoids interpolation-errors leading to v>c.
  ! For (non-ideal) MHD, we directly interpolate the electric field as
  ! there is no such constraint.

  subroutine get_vec(ix,igrid,x,tloc,vec)
    use mod_global_parameters
    use mod_usr_methods, only: usr_particle_analytic

    integer,intent(in)                                 :: ix(ndir) !< Indices in gridvars
    integer,intent(in)                                 :: igrid
    double precision,dimension(3), intent(in)          :: x
    double precision, intent(in)                       :: tloc
    double precision,dimension(ndir), intent(out)      :: vec
    double precision,dimension(ndir)                   :: vec1, vec2
    double precision                                   :: td
    integer                                            :: ic^D,idir

    if (associated(usr_particle_analytic)) then
      call usr_particle_analytic(ix, x, tloc, vec)
    else if (.not.time_advance) then
      do idir=1,ndir
        call interpolate_var(igrid,ixG^LL,ixM^LL,gridvars(igrid)%w(ixG^T,ix(idir)), &
             ps(igrid)%x(ixG^T,1:ndim),x,vec(idir))
      end do
    else
      do idir=1,ndir
        call interpolate_var(igrid,ixG^LL,ixM^LL,gridvars(igrid)%wold(ixG^T,ix(idir)), &
             ps(igrid)%x(ixG^T,1:ndim),x,vec1(idir))
        call interpolate_var(igrid,ixG^LL,ixM^LL,gridvars(igrid)%w(ixG^T,ix(idir)), &
             ps(igrid)%x(ixG^T,1:ndim),x,vec2(idir))
      end do
      td = (tloc - global_time) / dt
      vec(:) = vec1(:) * (1.0d0 - td) + vec2(:) * td
    end if

  end subroutine get_vec

  !> Get Lorentz factor from relativistic momentum
  pure subroutine get_lfac(u,lfac)
    use mod_global_parameters, only: ndir, c_norm
    double precision,dimension(3), intent(in)       :: u
    double precision, intent(out)                      :: lfac

    if (relativistic) then
       lfac = sqrt(1.0d0 + sum(u(:)**2)/c_norm**2)
    else
       lfac = 1.0d0
    end if
  end subroutine get_lfac

  !> Get Lorentz factor from velocity
  pure subroutine get_lfac_from_velocity(v,lfac)
    use mod_global_parameters, only: ndir, c_norm
    double precision,dimension(3), intent(in)       :: v
    double precision, intent(out)                      :: lfac

    if (relativistic) then
       lfac = 1.0d0 / sqrt(1.0d0 - sum(v(:)**2)/c_norm**2)
    else
       lfac = 1.0d0
    end if

  end subroutine get_lfac_from_velocity

  subroutine cross(a,b,c)
    use mod_global_parameters
    use mod_geometry
    double precision, dimension(ndir), intent(in)   :: a,b
    double precision, dimension(ndir), intent(out)  :: c

    ! ndir needs to be three for this to work!!!
    if(ndir==3) then
      select case(coordinate)
      case(Cartesian,Cartesian_stretched)
        c(1) = a(2)*b(3) - a(3)*b(2)
        c(2) = a(3)*b(1) - a(1)*b(3)
        c(3) = a(1)*b(2) - a(2)*b(1)
      case (cylindrical)
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

  subroutine interpolate_var(igrid,ixI^L,ixO^L,gf,x,xloc,gfloc)
    use mod_global_parameters

    integer, intent(in)                   :: igrid,ixI^L, ixO^L
    double precision, intent(in)          :: gf(ixI^S)
    double precision, intent(in)          :: x(ixI^S,1:ndim)
    double precision, intent(in)          :: xloc(1:3)
    double precision, intent(out)         :: gfloc
    double precision                      :: xd^D
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

  subroutine read_particles_snapshot(file_exists)
    use mod_global_parameters

    logical,intent(out)             :: file_exists
    character(len=std_len)          :: filename
    integer                         :: mynpayload, mynparticles, pos

    ! some initialisations:
    nparticles_on_mype = 0
    mynparticles       = 0
    nparticles         = 0
    it_particles       = 0

    ! open the snapshot file on the headnode
    file_exists=.false.
    if (mype == 0) then
!      write(filename,"(a,a,i4.4,a)") trim(base_filename),'_particles',snapshotini,'.dat'
      ! Strip restart_from_filename of the ending 
      pos = scan(restart_from_file, '.dat', back=.true.)
      write(filename,"(a,a,i4.4,a)") trim(restart_from_file(1:pos-8)),'_particles',snapshotini,'.dat'
      INQUIRE(FILE=filename, EXIST=file_exists)
      if (.not. file_exists) then
        write(*,*) 'WARNING: File '//trim(filename)//' with particle data does not exist.'
        write(*,*) 'Initialising particles from user or default routines'
      else
        open(unit=unitparticles,file=filename,form='unformatted',action='read',status='unknown',access='stream')
        read(unitparticles) nparticles,it_particles,mynpayload
        if (mynpayload .ne. npayload) &
             call mpistop('npayload in restart file does not match npayload in mod_particles')
      end if
    end if

    call MPI_BCAST(file_exists,1,MPI_LOGICAL,0,icomm,ierrmpi)
    if (.not. file_exists) return

    call MPI_BCAST(nparticles,1,MPI_INTEGER,0,icomm,ierrmpi)
    call MPI_BCAST(it_particles,1,MPI_INTEGER,0,icomm,ierrmpi)

    do while (mynparticles .lt. nparticles)
      if (mype == 0) then
        do while (nparticles_on_mype .lt. nparticles_per_cpu_hi &
             .and. mynparticles .lt. nparticles)
          call read_from_snapshot
          mynparticles = mynparticles + 1
        end do
      end if ! mype==0
      call MPI_BCAST(mynparticles,1,MPI_INTEGER,0,icomm,ierrmpi)
      call comm_particles_global
    end do

    if (mype == 0) close(unit=unitparticles)

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
    read(unitparticles) particle(index)%self%time
    read(unitparticles) particle(index)%self%dt
    read(unitparticles) particle(index)%self%x
    read(unitparticles) particle(index)%self%u
    read(unitparticles) particle(index)%payload(1:npayload)

!    if (particle(index)%self%follow) print*, 'follow index:', index

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

    if (.not. write_snapshot) return

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
    write(unitparticles) myparticle%time
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
        filename = make_particle_filename(particle(ipart)%self%time,particle(ipart)%self%index,'individual')
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

  subroutine select_active_particles(end_time)
    use mod_global_parameters
    double precision, intent(in) :: end_time

    integer          :: ipart, iipart
    logical          :: activate
    double precision :: t_min_mype

    t_min_mype = bigdouble
    nparticles_active_on_mype = 0

    do iipart = 1, nparticles_on_mype
      ipart      = particles_on_mype(iipart);
      activate   = (particle(ipart)%self%time < end_time)
      t_min_mype = min(t_min_mype, particle(ipart)%self%time)

      if (activate) then
        nparticles_active_on_mype = nparticles_active_on_mype + 1
        particles_active_on_mype(nparticles_active_on_mype) = ipart
      end if
    end do

    call MPI_ALLREDUCE(t_min_mype, min_particle_time, 1, MPI_DOUBLE_PRECISION, &
         MPI_MIN, icomm, ierrmpi)

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
    use mod_slice, only: get_igslice
    use mod_global_parameters

    double precision, intent(in) :: x(3)
    integer, intent(out)         :: igrid_particle, ipe_particle
    integer                      :: ig(ndir,nlevelshi), ig_lvl(nlevelshi)
    integer                      :: idim, ic(ndim)
    type(tree_node_ptr)          :: branch

    ! first check if the particle is in the domain
    if (.not. particle_in_domain(x)) then
      igrid_particle = -1
      ipe_particle   = -1
      return
    end if

    ! get the index on each level
    do idim = 1, ndim
      call get_igslice(idim,x(idim), ig_lvl)
      ig(idim,:) = ig_lvl
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

  !> Check if particle is inside computational domain
  logical function particle_in_domain(x)
    use mod_global_parameters
    double precision, dimension(ndim), intent(in)  :: x

    particle_in_domain = all(x >= [ xprobmin^D ]) .and. &
         all(x < [ xprobmax^D ]) ! Jannis: as before < instead of <= here, necessary?
  end function particle_in_domain

  !> Quick check if particle is still in igrid
  logical function particle_in_igrid(ipart, igrid)
    use mod_global_parameters
    integer, intent(in) :: igrid, ipart
    double precision    :: x(ndim), grid_rmin(ndim), grid_rmax(ndim)

    ! First check if the igrid is still there
    if (.not. allocated(ps(igrid)%w)) then
      particle_in_igrid = .false.
    else
      grid_rmin         = [ {rnode(rpxmin^D_,igrid)} ]
      grid_rmax         = [ {rnode(rpxmax^D_,igrid)} ]
      x                 = particle(ipart)%self%x(1:ndim)
      particle_in_igrid = all(x >= grid_rmin) .and. all(x < grid_rmax)
    end if
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
      case (neighbor_boundary) ! boundary
        ! do nothing
      case (neighbor_coarse) ! fine-coarse
        call ipe_fc(i^D,igrid,ipe_is_neighbor)
      case (neighbor_sibling) ! same level
        call ipe_srl(i^D,igrid,ipe_is_neighbor)
      case (neighbor_fine) ! coarse-fine
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

  subroutine write_particle_output()
    use mod_global_parameters

    character(len=std_len) :: filename
    integer                         :: ipart,iipart
    type(particle_t), dimension(nparticles_per_cpu_hi)  :: send_particles
    double precision, dimension(npayload,nparticles_per_cpu_hi)  :: send_payload
    integer                         :: send_n_particles_for_output
    integer                         :: nout
    double precision                :: tout

    if (write_individual) then
      call output_individual()
    end if

    if (write_ensemble) then
      send_n_particles_for_output = 0

      do iipart=1,nparticles_on_mype
        ipart=particles_on_mype(iipart);

        ! have to send particle to rank zero for output
        send_n_particles_for_output = send_n_particles_for_output + 1
        send_particles(send_n_particles_for_output) = particle(ipart)%self
        send_payload(1:npayload,send_n_particles_for_output) = particle(ipart)%payload(1:npayload)
      end do

      call output_ensemble(send_n_particles_for_output,send_particles, &
           send_payload,'ensemble')
    end if

  end subroutine write_particle_output

  character(len=128) function make_particle_filename(tout,index,type)
    use mod_global_parameters

    character(len=*), intent(in)    :: type
    double precision, intent(in)    :: tout
    integer, intent(in)             :: index
    integer                         :: nout, mysnapshot

    select case(type)
    case ('ensemble')
      nout = nint(tout / dtsave_particles)
      write(make_particle_filename,"(a,a,i6.6,a)") trim(base_filename),'_ensemble',nout,'.csv'
    case ('destroy')
      write(make_particle_filename,"(a,a)") trim(base_filename),'_destroyed.csv'
    case ('individual')
      write(make_particle_filename,"(a,a,i6.6,a)") trim(base_filename),'_particle',index,'.csv'
    end select

  end function make_particle_filename

  subroutine output_ensemble(send_n_particles_for_output,send_particles,send_payload,typefile)
    use mod_global_parameters

    integer, intent(in)             :: send_n_particles_for_output
    type(particle_t), dimension(send_n_particles_for_output), intent(in)  :: send_particles
    double precision, dimension(npayload,send_n_particles_for_output), intent(in)  :: send_payload
    character(len=*), intent(in)    :: typefile
    character(len=std_len)              :: filename
    type(particle_t), dimension(nparticles_per_cpu_hi)  :: receive_particles
    double precision, dimension(npayload,nparticles_per_cpu_hi) :: receive_payload
    integer                         :: status(MPI_STATUS_SIZE)
    integer,dimension(0:npe-1)      :: receive_n_particles_for_output_from_ipe
    integer                         :: ipe, ipart, nout
    logical                         :: file_exists

    receive_n_particles_for_output_from_ipe(:) = 0

    call MPI_ALLGATHER(send_n_particles_for_output, 1, MPI_INTEGER, &
         receive_n_particles_for_output_from_ipe, 1, MPI_INTEGER, icomm, ierrmpi)

    ! If there are no particles to be written, skip the output
    if (sum(receive_n_particles_for_output_from_ipe(:)) == 0) return

    if (mype > 0) then
      call MPI_SEND(send_particles,send_n_particles_for_output,type_particle,0,mype,icomm,ierrmpi)
      call MPI_SEND(send_payload,npayload*send_n_particles_for_output,MPI_DOUBLE_PRECISION,0,mype,icomm,ierrmpi)
    else
      ! Create file and write header
      if(typefile=='destroy') then ! Destroyed file
        write(filename,"(a,a,i6.6,a)") trim(base_filename) // '_', &
             trim(typefile) // '.csv'
        n_output_destroyed= n_output_destroyed + 1
        inquire(file=filename, exist=file_exists)

        if (.not. file_exists) then
          open(unit=unitparticles, file=filename)
          write(unitparticles,"(a)") trim(csv_header)
        else
          open(unit=unitparticles, file=filename, access='append')
        end if

      else ! Ensemble file
        write(filename,"(a,a,i6.6,a)") trim(base_filename) // '_', &
             trim(typefile) // '_', n_output_ensemble,'.csv'
        n_output_ensemble = n_output_ensemble + 1
        open(unit=unitparticles,file=filename)
        write(unitparticles,"(a)") trim(csv_header)
      end if

      ! Write own particles
      do ipart=1,send_n_particles_for_output
        call output_particle(send_particles(ipart),send_payload(1:npayload,ipart),0,unitparticles)
      end do

      ! Write particles from other tasks
      do ipe=1,npe-1
        call MPI_RECV(receive_particles,receive_n_particles_for_output_from_ipe(ipe),type_particle,ipe,ipe,icomm,status,ierrmpi)
        call MPI_RECV(receive_payload,npayload*receive_n_particles_for_output_from_ipe(ipe),MPI_DOUBLE_PRECISION,ipe,ipe,icomm,status,ierrmpi)
        do ipart=1,receive_n_particles_for_output_from_ipe(ipe)
          call output_particle(receive_particles(ipart),receive_payload(1:npayload,ipart),ipe,unitparticles)
        end do ! ipart
      end do ! ipe

      close(unitparticles)
    end if

  end subroutine output_ensemble

  subroutine output_individual()
    use mod_global_parameters
    character(len=std_len) :: filename
    integer                :: ipart,iipart,ipe
    logical                :: file_exists

    do iipart=1,nparticles_on_mype
      ipart=particles_on_mype(iipart)

      ! should the particle be dumped?
      if (particle(ipart)%self%follow) then
        write(filename,"(a,a,i6.6,a)") trim(base_filename), '_particle_', &
             particle(ipart)%self%index, '.csv'
        inquire(file=filename, exist=file_exists)

        ! Create empty file and write header on first iteration, or when the
        ! file does not exist yet
        if (it_particles == 0 .or. .not. file_exists) then
          open(unit=unitparticles, file=filename)
          write(unitparticles,"(a)") trim(csv_header)
        else
          open(unit=unitparticles, file=filename, access='append')
        end if

        call output_particle(particle(ipart)%self,&
             particle(ipart)%payload(1:npayload),mype,unitparticles)

        close(unitparticles)
      end if
    end do

  end subroutine output_individual

  subroutine output_particle(myparticle,payload,ipe,file_unit)
    use mod_global_parameters
    use mod_geometry

    type(particle_t),intent(in)  :: myparticle
    double precision, intent(in) :: payload(npayload)
    integer, intent(in)          :: ipe
    integer, intent(in)          :: file_unit
    double precision             :: x(3), u(3)
    integer                      :: icomp

    ! convert and normalize the position
    if (nocartesian) then
      select case(coordinate)
        case(Cartesian,Cartesian_stretched)
          x = myparticle%x*length_convert_factor
        case(cylindrical)
          x(r_)   = myparticle%x(r_)*length_convert_factor
          x(z_)   = myparticle%x(z_)*length_convert_factor
          x(phi_) = myparticle%x(phi_)
        case(spherical)
          x(:) = myparticle%x(:)
          x(1) = x(1)*length_convert_factor
      end select
    else
      call partcoord_to_cartesian(myparticle%x,x)
      x(:) = x(:)*length_convert_factor
    end if

    u = myparticle%u(1:3) * integrator_velocity_factor

    write(file_unit, csv_format) myparticle%time, myparticle%dt, x(1:3), &
         u, payload(1:npayload), ipe, it_particles, &
         myparticle%index

  end subroutine output_particle

  !> convert to cartesian coordinates
  subroutine partcoord_to_cartesian(xp,xpcart)
    ! conversion of particle coordinate from cylindrical/spherical to cartesian coordinates done here
    ! NOT converting velocity components: TODO?
    ! Also: nullifying values lower than smalldouble
    use mod_global_parameters
    use mod_geometry

    double precision, intent(in)  :: xp(1:3)
    double precision, intent(out) :: xpcart(1:3)

    select case (coordinate)
       case (Cartesian,Cartesian_stretched,Cartesian_expansion)
          xpcart(1:3)=xp(1:3)
       case (cylindrical)
          xpcart(1)=xp(1)*cos(xp(phi_))
          xpcart(2)=xp(1)*sin(xp(phi_))
          xpcart(3)=xp(z_)
       case (spherical)
          xpcart(1)=xp(1)*sin(xp(2))*cos(xp(3))
          {^IFTWOD
          xpcart(2)=xp(1)*cos(xp(2))
          xpcart(3)=xp(1)*sin(xp(2))*sin(xp(3))}
          {^IFTHREED
          xpcart(2)=xp(1)*sin(xp(2))*sin(xp(3))
          xpcart(3)=xp(1)*cos(xp(2))}
       case default
          write(*,*) "No converter for coordinate=",coordinate
       end select

  end subroutine partcoord_to_cartesian

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

      ! first check if the particle should be destroyed (TODO: user-defined criterion)
      if ( .not.time_advance .and. particle(ipart)%self%time .gt. tmax_particles ) then
        destroy_n_particles_mype  = destroy_n_particles_mype + 1
        particle_index_to_be_destroyed(destroy_n_particles_mype) = ipart
        cycle
      end if

      ! Then check user-defined condition
      if (associated(usr_destroy_particle)) then
        if (usr_destroy_particle(particle(ipart))) then
          destroy_n_particles_mype  = destroy_n_particles_mype + 1
          particle_index_to_be_destroyed(destroy_n_particles_mype) = ipart
          cycle
        end if
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
        write(*,*) particle(ipart)%self%time, ': particle',ipart,' has left at it ',it_particles,' on pe', mype
        write(*,*) particle(ipart)%self%time, ': particles remaining:',nparticles, '(total)', nparticles_on_mype-1, 'on pe', mype
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

end module mod_particle_base
