module mod_usr
  use mod_mhd
  use mod_eos, only: eos
  use mod_lookup_table
  use mod_bc_data
  !SI units, constants
  use mod_constants, only: mp_SI, kB_SI, miu0_SI
  use, intrinsic :: ieee_arithmetic
  use mod_datacube
  use mod_particles
  implicit none

  character(len=20)                              :: printsettingformat
  double precision                               :: omega_frame
  double precision, parameter                    :: omega_HEEQ = 2.d0*dpi/(365.25d0*24.d0) ! rad/h
  character(len=500)                             :: amr_criterion, cme_parameter_file
  integer, parameter                             :: AMR_NONE=0, AMR_SHOCK=1, AMR_LONWINDOW=2, AMR_TRACING=3, AMR_SPIRAL=4
  integer                                        :: amr_mode = AMR_NONE

  double precision                               :: spiral_lon_min
  double precision                               :: spiral_lon_max
  double precision                               :: spiral_speed

  double precision                               :: amr_start_hour

  integer                                        ::  num_cmes, relaxation, cme_insertion
  type satellite_pos
    real(kind=8), dimension(:,:), allocatable    :: positions
  end type satellite_pos

  !Shared over subroutines

  real(kind=8), allocatable                      :: coord_grid_init(:,:,:),variables_init(:,:,:)
  type(satellite_pos), dimension(:), allocatable :: positions_list
  character(len=250), dimension(10)              :: trajectory_list
  character(len=11), dimension(10)               :: sat_name = (/'earth      ','mars       ','mercury    ','venus      ', &
                                                                 'sta        ','stb        ','psp        ','SolO       ', &
                                                                 'bepi       ','juno       '/)

  integer, dimension(10)                          :: which_satellite = (/0, 0, 0, 0, 0, 0, 0, 0, 0, 0/)     ! intended order: earth, mars, mercury, venus, sta, stb, psp, solo, bepi, juno
  integer, dimension(10)                          :: sat_indx = (/0, 0, 0, 0, 0, 0, 0, 0, 0, 0/)     ! intended order: earth, mars, mercury, venus, sta, stb, psp, solo, bepi, juno
  integer                                         :: sat_count=0, zero_count=0

  integer, dimension(10)                          :: last_index = (/0, 0, 0, 0, 0, 0, 0, 0, 0, 0/)  ! intended order: earth, mars, mercury, venus, sta, stb, psp, solo, bepi, juno
  integer, dimension(10)                          :: last_index_s = (/0, 0, 0, 0, 0, 0, 0, 0, 0, 0/)

  integer, dimension(:,:), allocatable           :: starting_index, cme_index    ! first coordinate: satellite index; second coordinate: cme index
  integer, dimension(10)                          :: magnetogram_index = (/0, 0, 0, 0, 0, 0, 0, 0,0,0/)

  ! CME parameters and simulation details from the parameter file
  ! define my cme parameters here
  character(len=100), dimension(:), allocatable   :: cme_type, cme_date
  integer, dimension(:), allocatable             :: cme_year, cme_month, cme_day, cme_hour, cme_minute, cme_second
  double precision, dimension(:), allocatable    ::  vr_cme, w_half, clt_cme, lon_cme, rho_cme, temperature_cme
  double precision, dimension(:), allocatable    :: timestamp, longitudes_fix
  double precision, dimension(:,:), allocatable  :: time_difference_cme_magn
  double precision, dimension(:), allocatable    :: lon_updated, lon_original
  integer             :: magnetogram_timestamp(6)
  character(len=20)   :: magnetogram_time

  integer             :: cme_exists
      public :: bc_data_get_3d



contains

  subroutine usr_params_read(files)
    implicit none
    character(len=*), intent(in) :: files(:)
    integer :: n, ios

    namelist /rotating_frame_list/ omega_frame
    namelist /icarus_list/ amr_criterion,  num_cmes, relaxation, cme_insertion, &
    cme_parameter_file, magnetogram_time,  amr_start_hour

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old", action="read")
       ! Read in the order they appear;
      ios = 0
      read(unitpar, nml=rotating_frame_list, iostat=ios)
      ios = 0
      read(unitpar, nml=icarus_list,        iostat=ios)
      close(unitpar)
    end do
    
    if (num_cmes  == 0) cme_insertion = 0
    ! --- Map AMR mode (case-sensitive) ---
    select case (trim(adjustl(amr_criterion)))
      case ('shock');     amr_mode = AMR_SHOCK
      case ('lonwindow'); amr_mode = AMR_LONWINDOW
      case ('tracing');   amr_mode = AMR_TRACING
      case default;       amr_mode = AMR_NONE
    end select
    
   
    ! --- Parse magnetogram_time 'YYYY_MM_DD_HH_MM_SS'  ---
      read (magnetogram_time(1:4),*) magnetogram_timestamp(1)
      read (magnetogram_time(6:7),*) magnetogram_timestamp(2)
      read (magnetogram_time(9:10),*) magnetogram_timestamp(3)
      read (magnetogram_time(12:13),*) magnetogram_timestamp(4)
      read (magnetogram_time(15:16),*) magnetogram_timestamp(5)
      read (magnetogram_time(18:19),*) magnetogram_timestamp(6)
      
  end subroutine usr_params_read

  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods


    call usr_params_read(par_files)
    usr_set_parameters  => initglobaldata_usr
    usr_init_one_grid   => initonegrid_usr
    usr_aux_output      => specialvar_output
    usr_add_aux_names   => specialvarnames_output
    usr_special_bc      => specialbound_usr
    usr_refine_grid     => specialrefine_grid
    usr_source          => specialsource
    usr_create_particles => generate_particles
    usr_particle_position => move_particle
  


    call set_coordinate_system('spherical_3D')


    call mhd_activate()

   



    !  Note: mhd_activate sets the physical units used by MPI-AMRVAC as governed
    ! in subroutine mhd_phys_init (in mod_mhd_phys.t) which in turn calls
    ! subroutine mhd_physical_units (also in mod_mhd_phys.t)
    !  There, the parameters SI_unit, eq_state_units, mhd_partial_ionization enter
    !  Sometime we use eos%He_abundance, H_ion_fr, He_ion_fr, He_ion_fr2
    !  Moreover, we use 3 out of
    !      unit_density, unit_numberdensity, unit_length, unit_time, unit_velocity,
    !      unit_pressure, unit_magneticfield, unit_temperature,
    !      unit_mass, unit_charge
    !  Note we have factor RR (p=RR rho T)
    printsettingformat='(1x,A50,ES15.7)'
    if(mype==0) then
      write(*,*)'----------------PARAMETERS--   ----------------------'
      write(*,printsettingformat) "eos%gamma ",eos%gamma
      write(*,printsettingformat) "mhd_eta ",mhd_eta
      write(*,*)'----------------BEGIN UNITS  ------------------------'
      write(*,*)'----------------UNIT CONTROLS------------------------'
      write(*,*) "SI_unit",SI_unit
      write(*,*) "eos_type",eos%eos_type
      write(*,printsettingformat) "eos%He_abundance",eos%He_abundance
      write(*,printsettingformat) "H_ion_fr",H_ion_fr
      write(*,printsettingformat) "He_ion_fr",He_ion_fr
      write(*,printsettingformat) "He_ion_fr2",He_ion_fr2
      write(*,*)'----------------UNIT CONTROLS------------------------'
      write(*,*)'----------------UNITS---------        ---------------'
      write(*,printsettingformat) "unit density in g/cm-3: ",unit_density
      write(*,printsettingformat) "unit number density cm-3: ",unit_numberdensity
      write(*,printsettingformat) "unit length in cm: ",unit_length
      write(*,printsettingformat) "unit time  in seconds: ",unit_time
      write(*,printsettingformat) "unit velocity in cm/s: ",unit_velocity
      write(*,printsettingformat) "unit pressure in cgs: ",unit_pressure
      write(*,printsettingformat) "unit magnetic field in gauss: ",unit_magneticfield
      write(*,printsettingformat) "unit temperature in K: ",unit_temperature
      write(*,printsettingformat) "unit mass in g: ",unit_mass
      write(*,printsettingformat) "unit charge: ",unit_charge
      write(*,*)'----------------UNITS--------------------------------'
      write(*,printsettingformat) "p=RR rho T with RR: ",RR
      write(*,*)'----------------END UNITS----------------------------'
    end if

  end subroutine usr_init

  subroutine initglobaldata_usr
    use mod_global_parameters
    logical, save       :: firstglobalusr=.true.
    integer             :: AllocateStatus, DeAllocateStatus
    !Parameters related to the unit conversions
    double precision    :: Lunit_in, Tunit_in, Rhounit_in, Vunit_in, Bunit_in, Eunit_in, Punit_in
    !Parameters related to the read-in of the boundary file
    character(len=50)   :: earth_trajectory, mars_trajectory, venus_trajectory
    character(len=50)   :: sta_trajectory, stb_trajectory, mercury_trajectory

    character(len=200)  :: path_satellite_trajectories
    integer             :: nr_colat, nr_lon, k, n, i


    if(firstglobalusr) then
      call print_initial_information()
      !Read-in coronal model
      !Coronal model at boundary has 2 coordinates: colat and lon
      !Coronal model has data for 4 parameters: vr, n, T, Br
      path_satellite_trajectories = './orbit/'


      nr_colat = lt_3d(1)%n_points(1)
      nr_lon = lt_3d(1)%n_points(2)
      ! WARNING: ASSUMES STENCIL IS USING 2 GHOSTCELLS: may need to GENERALIZE!
      !We have 4 extra points for longitude because the boundary is periodic
      ALLOCATE(coord_grid_init(nr_colat, nr_lon+4, ndim-1), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) call mpistop('*** Not enough memory ***')
      ALLOCATE(variables_init(nr_colat, nr_lon+4, 4), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) call mpistop('*** Not enough memory ***')

     


      ! read in cme parameters
      if (num_cmes == 0) then
        ALLOCATE(timestamp(1))
        ALLOCATE(cme_index(10,1))
        ALLOCATE(starting_index(10,1))
        ALLOCATE(time_difference_cme_magn(10,1))
        timestamp(:) = 0.0
        cme_index(:,:) = 1   ! ensure not zero (valid Fortran index)
        starting_index(:,:) = 1
        time_difference_cme_magn(:,:) = 0.0
      else
        call read_cme_parameters(cme_parameter_file)
      end if
      ! Initialize cme starting index in the trajectory file, cme index in the trajectory file and the time difference between the start and cme indexes
      if (num_cmes > 0) then
        do n = 1, num_cmes
          do i = 1, 10
            cme_index(i, n) = 0
            starting_index(i, n) = 0
            time_difference_cme_magn(i, n) = 0
          end do
        end do
        else
          do i = 1, 10
            cme_index(i, 1) = 0
            starting_index(i, 1) = 0
            time_difference_cme_magn(i, 1) = 0
          end do
       end if


        do i=1, 10
            call find_trajectory_file(i, path_satellite_trajectories)
        end do

        ALLOCATE(positions_list(10), STAT=AllocateStatus)


      ! for each satellite, read the trajectory data and save in the arrays of time and locations
      do i = 1, 10
        if (which_satellite(i)==1) then
          sat_indx(i-zero_count) = i
          sat_count = sat_count+1
          call read_satellite_trajectory(trajectory_list(i), i)
          
        end if
         if (which_satellite(i) == 0) then
          zero_count = zero_count+1
         end if
      end do

      if (mype == 0) then
        write(*,*) 'Particle -> spacecraft mapping:'
        do i = 1, sat_count
           write(*,'(I6,2X,I2,2X,A)') i, sat_indx(i), trim(sat_name(sat_indx(i)))
        end do
      end if


      ! calculate timestamp for cme insertion
      timestamp(:) = relaxation*24.0+cme_insertion*24.0
      if (num_cmes >0) then
        do k=1, num_cmes
          timestamp(k) = timestamp(k) + time_difference_cme_magn(1, k)
        end do
        call cme_insertion_longitudes_fix()
      end if

      firstglobalusr = .false.

    end if

    !> gamma is set in the parfile (&eos_list gamma=1.5d0), read once at
    !> eos_init so the derived gamma constants stay consistent. Do not set
    !> eos%gamma from user code.
    call set_units(Lunit_in, Tunit_in, Rhounit_in, Vunit_in, Bunit_in, Eunit_in, Punit_in)

    w_convert_factor(rho_) = Rhounit_in ! in km/m^3
    w_convert_factor(mom(1)) =Vunit_in*1d-3  ! in km/s
    w_convert_factor(mom(2)) = w_convert_factor(mom(1))
    w_convert_factor(mom(3)) = w_convert_factor(mom(1))
    w_convert_factor(p_) = Punit_in  ! in kg/s/m^2
    w_convert_factor(mag(1)) = Bunit_in*1d9 ! in nT
    w_convert_factor(mag(2)) = w_convert_factor(mag(1))
    w_convert_factor(mag(3)) = w_convert_factor(mag(1))

  end subroutine initglobaldata_usr



  subroutine generate_particles(n_particles, x, v, q, m, follow)
    use mod_particles
    integer, intent(in)           :: n_particles
    double precision, intent(out) :: x(3, n_particles), v(3, n_particles)
    double precision, intent(out) :: q(n_particles), m(n_particles)
    logical, intent(out)          :: follow(n_particles)
    integer                       :: i, satellite_index, delta_sat

    ! only create as many particles as we have satellites
    sat_count = min(sat_count, n_particles)
    do i = 1, sat_count
       v(:, i)       = 0.d0
       q(i)          = 0.d0
       m(i)          = 0.d0
       call initialize_particle(x(:, i), sat_indx(i) )
       follow(i)     = .true.
    end do

    ! explicitly disable the rest
    do i = sat_count+1, n_particles
       x(:, i)    = 0.d0
       v(:, i)    = 0.d0
       q(i)       = 0.d0
       m(i)       = 0.d0
       follow(i)  = .false.
    end do

  end subroutine generate_particles

  subroutine initialize_particle(x, satellite_index)
    use mod_global_parameters
    implicit none
    double precision, intent(out) :: x(3)            ! (r, theta, phi)
    integer,          intent(in)  :: satellite_index

    integer :: npos, idx0
    double precision :: t0
    double precision :: r0, lat0, phi_heeq, phi_rot

    ! guard
    npos = size(positions_list(satellite_index)%positions, 2)
    if (npos <= 0) then
      x(:) = 0.d0
      return
    end if

    ! Start index for this spacecraft (clamped)
    idx0 = max(1, min(npos, starting_index(satellite_index, 1)))

    ! Hours between magnetogram epoch and simulation start:
    t0 = (relaxation + cme_insertion)*24.d0

    ! Position at simulation start
    r0   = positions_list(satellite_index)%positions(7, idx0)
    lat0 = positions_list(satellite_index)%positions(8, idx0)   ! latitude [rad]
    x(1) = r0
    x(2) = dpi/2.d0 - lat0                                      ! theta = pi/2 - lat

    ! HEEQ longitude at start index
    phi_heeq = positions_list(satellite_index)%positions(9, idx0)

    ! Convert to simulation rotating frame: subtract relative rotation over t0
    phi_rot = phi_heeq +  (omega_frame - omega_HEEQ)*t0 

    ! wrap to [0, 2π)
    phi_rot = modulo(phi_rot, 2.d0*dpi)
    if (phi_rot < 0.d0) phi_rot = phi_rot + 2.d0*dpi
    x(3) = phi_rot
  end subroutine initialize_particle



    subroutine move_particle(x, particle_id, told, tnew)
    use mod_global_parameters
    implicit none
    ! Inputs/outputs
    double precision, intent(inout) :: x(3)     ! (r, theta, phi) in rotating frame
    integer,          intent(in)    :: particle_id
    double precision, intent(in)    :: told, tnew  ! hours since simulation start

    ! Locals
    integer :: sidx                 ! mapped spacecraft index (1..10)
    integer :: npos, base           ! number of samples, start index
    double precision :: idx0, idx1  ! fractional minute indices at told/tnew
    integer :: i0, j0, i1, j1       ! bracketing minute samples [i, i+1]
    double precision :: a0, a1      ! linear weights within minute [0,1)

    double precision :: lon0f, lon0c, lon1f, lon1c
    double precision :: d0, d1, lon0, lon1, dlon, dframe
    double precision :: r_f, r_c, lat_f, lat_c, lat_new


    !--- Map particle -> spacecraft; exit if invalid
    sidx = sat_indx(particle_id)
    if (sidx <= 0) return

    !--- Ephemeris size
    npos = size(positions_list(sidx)%positions, 2)
    if (npos <= 1) return     ! need at least 2 samples to interpolate

    !--- Starting index on minute grid (1-based, clamped)
    base = max(1, starting_index(sidx, 1))
    if (base > npos) base = npos

    !--- Hours -> fractional minute indices (relative to 'base')
    idx0 = dble(base) + told*60.0d0
    idx1 = dble(base) + tnew*60.0d0

    ! Low/high neighbors and weights at told
    i0 = max(1, min(npos-1, int(idx0)))     ! floor for positive idx
    j0 = i0 + 1
    a0 = idx0 - dble(i0)                    ! in [0,1)

    ! Low/high neighbors and weights at tnew
    i1 = max(1, min(npos-1, int(idx1)))
    j1 = i1 + 1
    a1 = idx1 - dble(i1)


    !--- Wrap-aware interpolation of inertial longitude at told
    lon0f = positions_list(sidx)%positions(9, i0)
    lon0c = positions_list(sidx)%positions(9, j0)
    d0    = modulo((lon0c - lon0f) + dpi, 2.d0*dpi) - dpi
    lon0  = lon0f + a0*d0

    ! and at tnew
    lon1f = positions_list(sidx)%positions(9, i1)
    lon1c = positions_list(sidx)%positions(9, j1)
    d1    = modulo((lon1c - lon1f) + dpi, 2.d0*dpi) - dpi
    lon1  = lon1f + a1*d1

    ! Inertial change over the step (wrap-aware)
    dlon  = modulo((lon1 - lon0) + dpi, 2.d0*dpi) - dpi

    ! subtract rotation of the simulation frame relative to HEEQ frame of spacecraft
    dframe = (omega_frame - omega_HEEQ) * (tnew - told)
    ! Advance azimuth in rotating frame and wrap to [0, 2π)
    x(3) = modulo(x(3) + dlon - dframe, 2.d0*dpi)
    if (x(3) < 0.d0) x(3) = x(3) + 2.d0*dpi

    !--- Interpolate radius at tnew
    r_f  = positions_list(sidx)%positions(7, i1)
    r_c  = positions_list(sidx)%positions(7, j1)
    x(1) = r_f + a1*(r_c - r_f)

    !--- Interpolate latitude at tnew and convert to co-latitude
    lat_f   = positions_list(sidx)%positions(8, i1)   ! latitude [rad]
    lat_c   = positions_list(sidx)%positions(8, j1)
    lat_new = lat_f + a1*(lat_c - lat_f)
    x(2)    = dpi/2.d0 - lat_new                      ! theta = π/2 − latitude

    ! Optional bookkeeping
    last_index_s(sidx) = i1
  end subroutine move_particle

  subroutine initonegrid_usr(ixG^L,ix^L,w,x)
    use mod_global_parameters

    logical, save:: first=.true.
    integer, intent(in)             :: ixG^L, ix^L
    double precision, intent(in)    :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)
    integer                         :: ix^D
    double precision    :: xloc(1:ndim)
    double precision    :: r_boundary
    integer             :: point11_clt, point11_lon, point22_clt, point22_lon
    double precision :: phi(ixmin2:ixmax2,ixmin3:ixmax3)
    double precision :: velocity2d(ixmin2:ixmax2,ixmin3:ixmax3)
    double precision :: rho2d(ixmin2:ixmax2,ixmin3:ixmax3)
    double precision :: p2d(ixmin2:ixmax2,ixmin3:ixmax3)
    double precision :: br2d(ixmin2:ixmax2,ixmin3:ixmax3)
    integer ::  idir, i

    double precision :: ur, rho, p, br, bphi, r, theta, sin_theta, u_phi_corot
    w(ix^S,1:nw) = zero

    r_boundary   = xprobmin1 !in R_sun
    phi = x(ixmin1,ixmin2:ixmax2,ixmin3:ixmax3,3) 

    velocity2d(ixmin2:ixmax2, ixmin3:ixmax3) = bc_data_get_3d(bc_data_ix(mom(1), 1), &
           x(ixmin1, ixmin2:ixmax2, ixmin3:ixmax3, 2), &
           phi, 0d0)
    rho2d(ixmin2:ixmax2, ixmin3:ixmax3) = bc_data_get_3d(bc_data_ix(rho_, 1), &
           x(ixmin1, ixmin2:ixmax2, ixmin3:ixmax3, 2), &
           phi, 0d0)
    p2d(ixmin2:ixmax2, ixmin3:ixmax3) = bc_data_get_3d(bc_data_ix(p_, 1), &
           x(ixmin1, ixmin2:ixmax2, ixmin3:ixmax3, 2), &
           phi, 0d0)
    br2d(ixmin2:ixmax2, ixmin3:ixmax3) = bc_data_get_3d(bc_data_ix(mag(1), 1), &
           x(ixmin1, ixmin2:ixmax2, ixmin3:ixmax3, 2), &
           phi, 0d0)




    do ix1 = ixmin1, ixmax1
    do ix2 = ixmin2, ixmax2
      do ix3 = ixmin3, ixmax3

        ur   = velocity2d(ix2, ix3)
        r    = x(ix1, ix2, ix3, 1)
        theta = x(ix1, ix2, ix3, 2)
        sin_theta = sin(theta)

        rho  = rho2d(ix2, ix3) * (r_boundary / r)**2
        p    = p2d(ix2, ix3)   * (r_boundary / r)**2
        br   = br2d(ix2, ix3)  * (r_boundary / r)**2
        
        u_phi_corot = -omega_frame * r * sin_theta ! if  radial flow as inner BC in the inertial frame
        
        !u_phi_corot = -omega_frame * (r - r_boundary) * sin_theta ! if radial flow as inner BC in the corotating frame
        bphi = 0.d0                                         ! (u_phi_corot / ur) * br (creates divB)


        ! Fill primitive variables
        w(ix1, ix2, ix3, rho_) = rho
        w(ix1, ix2, ix3, p_)   = p
        w(ix1, ix2, ix3, mom(1)) = ur
        w(ix1, ix2, ix3, mom(2)) = 0.d0
        w(ix1, ix2, ix3, mom(3)) = u_phi_corot

        w(ix1, ix2, ix3, mag(1)) = br
        w(ix1, ix2, ix3, mag(2)) = 0.d0
        w(ix1, ix2, ix3, mag(3)) = bphi

      end do
    end do
  end do

    !Convert to conserved values
    call eos%to_conserved(ixG^L,ix^L,w,x)

    if(mhd_n_tracer ==  1) then
       w(ix^S, tracer(1)) = 0.0d0
    end if

  end subroutine initonegrid_usr

  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
    use mod_global_parameters
    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)

    double precision :: divb(ixI^S), divmom(ixI^S)
    double precision :: v(ixI^S,ndir), divV(ixI^S), momentum(ixI^S, ndir)
    integer :: i
    double precision :: r_boundary
    double precision ::  r(ixO^S), theta(ixO^S), sin_theta(ixO^S)

    ! output divB1
    call get_divb(w,ixI^L,ixO^L,divb)
    w(ixO^S,nw+1)=divb(ixO^S)
   

    do i=1,ndir
      v(ixI^S,i)=w(ixI^S,mom(i))/w(ixI^S,rho_)
    end do

    call divvector(v,ixI^L,ixO^L,divV)
    !w(ixO^S,nw+2)=divV(ixO^S)*step_size(ixI^S)
    w(ixO^S,nw+2)=divV(ixO^S)
    do i=1,ndir
      momentum(ixI^S,i)=w(ixI^S,mom(i))
    end do

    call divvector(momentum,ixI^L,ixO^L,divmom)
    w(ixO^S,nw+3)=divmom(ixO^S)


   ! w(ixO^S,nw+4) = (v(ixO^S,3) + &
   ! omega_frame*(r)* sin_theta)*unit_velocity*1d-3 ! the unit in km/s

  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
    character(len=*) :: varnames

    varnames='divB divV div_mom'
  end subroutine specialvarnames_output


  

  subroutine specialsource(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, iw^LIM
    integer, intent(in)            :: ixO^L
    double precision, intent(in)    :: qdt, qtC, qt
    double precision, intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: Gravconst_Msun_normalized
    double precision :: omega_normalized, omega = 2.97d-6
    double precision   :: xloc(1:ndim)

    omega_normalized = omega * unit_length/unit_velocity

    !Gravity
    Gravconst_Msun_normalized = const_G*1d-3*const_MSun/1.0d3 / (unit_velocity**2 * unit_length )
    !momentum equation -> source = F = rho.g
    !energy equation   -> source  = v . F
    w(ixO^S,e_)  = w(ixO^S,e_)  - qdt*wCT(ixO^S,mom(1))*Gravconst_Msun_normalized/(x(ixO^S,1)**2)
    w(ixO^S,mom(1)) = w(ixO^S,mom(1)) - qdt*wCT(ixO^S,rho_)*Gravconst_Msun_normalized/(x(ixO^S,1)**2)

  end subroutine specialsource



  subroutine specialrefine_grid(igrid,level,ixI^L,ixO^L,qt,w,x,refine,coarsen)
    use mod_global_parameters
    implicit none
    integer, intent(inout)          :: refine, coarsen

    ! locals
    integer, intent(in)             :: igrid, level, ixI^L, ixO^L
    double precision, intent(in)    :: qt, w(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision                :: v(ixI^S,ndir), divV(ixI^S)
    double precision                :: momvec(ixI^S,ndir), divRhoV(ixI^S)
    double precision                :: rcell
    integer                         :: ix1, ix2, ix3, i, cme_hit
    double precision, parameter     :: thrV = -1.0d0   ! threshold for normalized div(V)
    double precision, parameter     :: thrR = -5.0d0   ! threshold for normalized div(rho*V)
    double precision, parameter     :: thrTR = 5.0d-3   ! threshold for normalized div(rho*V)
    double precision, parameter     :: fudge_cap  = 1.0d0   ! 
    double precision, parameter     :: eps_rho = 1.0d-30
    double precision, parameter     :: Vchar   = 2.6d0     ! 500 ~ 2.6*193 km/s (= unit_velocity)
    double precision, parameter     :: rhochar = 1.0d0     ! unit_density is characteristic density at inner boundary 
    double precision, parameter     :: rchar   = 21.5d0   ! inner boundary in unit_length


    double precision                :: phi_satellite, before_cme
    integer                         :: npos, i_now 
      
    double precision, parameter     :: halfw = 15.d0*dpi/180.d0 
    double precision, parameter     :: phi0  = 0.5d0*dpi 
    double precision, parameter     :: thr   = 1.0d0
    logical                         :: has_tracer


    ! default: coarsen
      refine  = -1
      coarsen = 1
        ! common precomputes used by several modes
      if (num_cmes > 0) then
        before_cme = (cme_index(1,1) - magnetogram_index(1))/60.0d0
      else
        before_cme = 0.0d0
      end if

      npos  = size(positions_list(1)%positions, 2)
      i_now = max(1, min(npos, starting_index(1,1) + floor(qt*60.0d0)))

      phi_satellite = positions_list(1)%positions(9, i_now) &
                    - omega_frame * (qt - (timestamp(1) - before_cme))
      phi_satellite = modulo(phi_satellite, 2.d0*dpi)
      if (phi_satellite < 0.d0) phi_satellite = phi_satellite + 2.d0*dpi

        select case (amr_mode)

      case (AMR_SHOCK)
        if (qt <= timestamp(1)) return
    
       ! ---- Phase 1: cheap angular mask pass 
        cme_hit = 0
        do ix2 = ixOmin2, ixOmax2
          do ix3 = ixOmin3, ixOmax3
            do i = 1, num_cmes
              if (qt < timestamp(i)) cycle
              call mask_cap( x(ixOmin1,ix2,ix3,2), &
                             x(ixOmin1,ix2,ix3,3), &
                             qt, fudge_cap, cme_hit, i )
              if (cme_hit == 1) exit
            end do
            if (cme_hit == 1) exit
          end do
          if (cme_hit == 1) exit
        end do

        if (cme_hit == 0) return    

        ! ---- Phase 2: compute divergences once
        do i = 1, ndir
          v(ixI^S,i)      = w(ixI^S,mom(i)) / max(w(ixI^S,rho_), eps_rho)
          momvec(ixI^S,i) = w(ixI^S,mom(i))
        end do
        call divvector(v,      ixI^L, ixO^L, divV)
        call divvector(momvec, ixI^L, ixO^L, divRhoV)

        ! ---- Phase 3: scan only masked lines; refine on first match
        do ix2 = ixOmin2, ixOmax2
          do ix3 = ixOmin3, ixOmax3
            do ix1 = ixOmin1, ixOmax1
              rcell = x(ix1,ix2,ix3,1) 
              if ( (rcell * divV(ix1,ix2,ix3)/Vchar)   < thrV .and. &
                   (rcell**3 *  divRhoV(ix1,ix2,ix3) / (rhochar * rchar**2 * Vchar))  < thrR) then
                refine  =  1
                coarsen = -1
                return
              end if
            end do
          end do
        end do

      case (AMR_LONWINDOW)
        ! refine if any cell is within fixed window around phi0 (wrap-safe)
        if (any( abs( modulo((x(ixI^S,3)-phi0)+dpi, 2.d0*dpi)-dpi ) <= halfw )) then
          refine=1; coarsen=-1
        end if

      case (AMR_TRACING)

        if (qt > timestamp(1)) then
         

         has_tracer = any( w(ixI^S, tracer(1)) > thrTR )

         !has_tracer = any( w(ixI^S, tracer(1)) > thr * w(ixI^S, rho_) )

         ! Old satellite/longitude window (DISABLED):
         ! double precision, parameter :: halfw = 30.d0*dpi/180.d0
         ! logical :: in_window
         ! in_window = any( abs( modulo((x(ixI^S,3) - phi_satellite) + dpi, 2.d0*dpi) - dpi ) <= halfw )

         if (has_tracer) then
           refine=1; coarsen=-1
         else
           refine=-1; coarsen=1
          end if
        end if

      case default
        ! leave defaults
      end select
  end subroutine specialrefine_grid

  subroutine mask_cap(clt_point, lon_point, qt, fudge, is_inside_cme, i)
    use mod_global_parameters
    implicit none
    double precision, intent(in)  :: clt_point, lon_point, qt, fudge
    integer,          intent(out) :: is_inside_cme
    integer,          intent(in)  :: i
    double precision :: dphi, cosd, cth, phi_cme


    phi_cme = modulo(lon_cme(i) - omega_frame * (qt - timestamp(i)) + 2.d0*dpi , 2.d0*dpi)
    
    ! wrap Δφ to (-π, π]
    dphi = modulo((lon_point - phi_cme) + dpi, 2.d0*dpi) - dpi

    ! cos(d) via spherical law of cosines (θ is co-lat)
    cosd = sin(clt_point)*sin(clt_cme(i))*cos(dphi) + cos(clt_point)*cos(clt_cme(i))
    cosd = max(-1.d0, min(1.d0, cosd))              ! clamp for safety

    ! inside if d ≤ w_half  <=>  cos d ≥ cos w_half
    cth  = cos(fudge * w_half(i))
    if (cosd >= cth) then
      is_inside_cme = 1
    else
      is_inside_cme = 0
    end if
  end subroutine mask_cap

  subroutine specialbound_usr(qdt,qt,ixI^L,ixO^L,iB,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L, iB
    double precision, intent(in) :: qdt,qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    integer             :: ix3,ix2
    double precision    :: vr_bc, vt_bc, vp_bc, br_bc, bt_bc, bp_bc, md_bc,rho_bc, temp_bc, p_bc, theta, phi
    integer             :: i,j,k, valuej, valuek, i_in
    integer             :: nr_r, nr_colat, nr_lon
    integer             :: point11_clt, point11_lon, point22_clt, point22_lon
    double precision    :: xloc(1:ndim)

    double precision    :: clt_zero, lon_zero, ts_magnetogram
    double precision    :: r_ref, r_g
    integer             :: mask_cme, n, local_check

    real, allocatable             :: mask(:)
    real, allocatable             :: vr(:), vp(:), vt(:), br(:), bt(:), bp(:), md(:), temp(:)

    double precision, allocatable :: clts(:), lons(:)
    double precision :: time
    integer :: idx, ndata, time_found

    select case(iB)
      case(1)! Lower radial boundary
        nr_colat = lt_3d(1)%n_points(1)
        nr_lon = lt_3d(1)%n_points(2)
! Should we want a radial velocity in the inertial frame, we need to set v3 = -omega_frame * r  and B3 = -omega_frame * r * sin\theta / v1
! we cannot use asymm BC for v3 and B3 in the par file, since that will kill any component in the corot-frame
      
      call eos%to_primitive(ixI^L,ixO^L,w,x)
      do ix3 = ixOmin3, ixOmax3
          do ix2 = ixOmin2,ixOmax2
            ! bc_data_set() has already filled ghost cells; keep its value at the
            ! ghost cell adjacent to the domain and enforce r^2 Br = const inward.
            br_bc  = w(ixOmax1, ix2, ix3, mag(1))
            rho_bc = w(ixOmax1, ix2, ix3, rho_)
            p_bc   = w(ixOmax1, ix2, ix3, p_)
            r_ref  = x(ixOmax1, ix2, ix3, 1)

            do i = ixOmin1, ixOmax1-1   ! fill all other ghost cells
              r_g = x(i, ix2, ix3, 1)
              w(i,ix2,ix3,mag(1)) = br_bc  * (r_ref / r_g)**2
              w(i,ix2,ix3,rho_)   = rho_bc * (r_ref / r_g)**2
              w(i,ix2,ix3,mom(3)) = -omega_frame * (r_g)*sin(x(i, ix2, ix3, 2))
              w(i,ix2,ix3,mag(3)) = -w(i,ix2,ix3,mag(1))*omega_frame * (r_g)&
              *sin(x(i, ix2, ix3, 2))/w(i,ix2,ix3,mom(1))

              ! (A) isothermal:
              ! w(i,ix2,ix3,p_) = p_bc * ( w(i,ix2,ix3,rho_) / rho_bc )

              ! (B) polytropic:
              w(i,ix2,ix3,p_)  =  p_bc * ( w(i,ix2,ix3,rho_) / rho_bc )**eos%gamma

            end do

            ! leave Btheta, Bphi as set by bc_data_set
 
            ! default: no tracer dye in the inner ghosts
            if (mhd_n_tracer >= 1) then
              do i = ixOmin1, ixOmax1
                w(i, ix2, ix3, tracer(1)) = 0.d0
              end do
            end if
           end do 
         end do

        do n = 1, num_cmes
            time = qt - timestamp(n)
            ts_magnetogram = timestamp(n) - relaxation*24 - cme_insertion*24

            if (time > 0.0d0) then
                call get_datacube_arrays(time, n, base_filename, time_found, clts, lons, mask, vr, vt, vp, br, bt, bp, md, temp)
            else
                time_found = 0
            end if

            if (time_found == 1) then
                ndata   = size(mask)
            end if

            do ix3 = ixOmin3, ixOmax3
                do ix2 = ixOmin2, ixOmax2
                    theta = x(ixOmin1+2, ix2, ix3, 2)
                    phi   = x(ixOmin1+2, ix2, ix3, 3)

                    mask_cme = 0
                   if (time_found == 1) then
                        call get_value(theta, phi, ts_magnetogram, &
                                       clts, lons, mask, vr, vp, vt, br, bt, bp, temp, md, ndata, mask_cme, &
                                       vr_bc, vp_bc, vt_bc, br_bc, bt_bc, bp_bc, md_bc, temp_bc)

                            if (mask_cme == 1) then
                            if (.not. ieee_is_nan(vr_bc)) then
                                w(ixOmin1, ix2, ix3, mom(1)) = vr_bc/unit_velocity
                                w(ixOmax1, ix2, ix3, mom(1)) = vr_bc/unit_velocity
                            end if

                            if ((.not. ieee_is_nan(vt_bc))) then
                                w(ixOmin1, ix2, ix3, mom(2)) = vt_bc/unit_velocity
                                w(ixOmax1, ix2, ix3, mom(2)) = vt_bc/unit_velocity
                            end if

                            if ((.not. ieee_is_nan(vp_bc))) then
                                w(ixOmin1, ix2, ix3, mom(3)) = vp_bc/unit_velocity
                                w(ixOmax1, ix2, ix3, mom(3)) = vp_bc/unit_velocity
                            end if

                            if (.not. ieee_is_nan(br_bc)) then
                                w(ixOmin1, ix2, ix3, mag(1)) = br_bc/unit_magneticfield
                                w(ixOmax1, ix2, ix3, mag(1)) = br_bc/unit_magneticfield
                            end if

                            if (.not. ieee_is_nan(bt_bc)) then
                                w(ixOmin1, ix2, ix3, mag(2)) = bt_bc/unit_magneticfield
                                w(ixOmax1, ix2, ix3, mag(2)) = bt_bc/unit_magneticfield
                            end if

                            if (.not. ieee_is_nan(bp_bc)) then
                                w(ixOmin1, ix2, ix3, mag(3)) = bp_bc/unit_magneticfield
                                w(ixOmax1, ix2, ix3, mag(3)) = bp_bc/unit_magneticfield
                            end if

                            if (.not. ieee_is_nan(md_bc)) then
                                w(ixOmin1, ix2, ix3, rho_) = md_bc/unit_density
                                w(ixOmax1, ix2, ix3, rho_) = md_bc/unit_density
                                w(ixOmax1, ix2, ix3, tracer(1)) = md_bc/unit_density
                                w(ixOmin1, ix2, ix3, tracer(1)) = md_bc/unit_density
                            end if

                            p_bc = md_bc/(0.5*mp_SI) * kB_SI * temp_bc
                            if (.not. ieee_is_nan(p_bc)) then
                                w(ixOmin1, ix2, ix3, p_) = p_bc/unit_pressure
                                w(ixOmax1, ix2, ix3, p_) = p_bc/unit_pressure
                            end if
                        else
                            w(ixOmax1, ix2, ix3, tracer(1)) = - w(ixOmin1+2, ix2, ix3, tracer(1))
                            w(ixOmin1, ix2, ix3, tracer(1)) = - w(ixOmin1+3, ix2, ix3, tracer(1))
                        end if
                    end if
                end do
            end do
        end do

        ! convert back to conserved values
        call eos%to_conserved(ixI^L, ixO^L, w, x)
    end select
  end subroutine specialbound_usr


  subroutine set_units(Lunit, Tunit, Rhounit, Vunit, Bunit, Eunit, Punit)
    use mod_global_parameters
    double precision, intent(out)    :: Lunit, Tunit, Rhounit, Vunit, Bunit, Eunit, Punit

    !Unit Length: [m] = 1 solar radii
    Lunit = const_RSun*1d-2
    !Unit Time : [s] = 1 hour
    Tunit = 6.0d1*6.0d1
    !Unit Mass density: [kg/m^3] = 1.6726d-13 - scaled so that rho is +- 1 in dimensionless form
    Rhounit = 1.6726d-19
    !Velocity unit
    Vunit = Lunit/Tunit
    !Magnetic filed unit
    Bunit = dsqrt( miu0_SI * Vunit**2 * Rhounit)
    !Energy density unit
    Eunit = Vunit*Vunit * Rhounit
    !Pressure unit: [kg/s/m^2]
    Punit = Vunit*Vunit * Rhounit

    if (mype ==0) then
      print *, ''//NEW_LINE('A'),' Length Unit: 1 Rs in m: ', Lunit
      print *, "Time Unit: 1 hour in s: ", Tunit
      print *, "Mass Density Unit: 1.6726d-19 so that rho = +- 1 in dimensionless form: ", Rhounit
      print *, "Velcoity Unit: Lunit/Tunit: ", Vunit
      print *, "Magnetic field Unit: sqrt(mu0 Vunit**2 Rhounit): ", Bunit
      print *, "Pressure Unit: Vunit*Vunit*rhounit: ",  Punit
      print *, ''//NEW_LINE('A')
      print *, "CME characteristic parameters"
      print *, ''//NEW_LINE('A')
      if (num_cmes > 0) then
        print *, "CME Timestamp[Y/M/D H/M/S]: ", cme_year, cme_month, cme_day, cme_hour, cme_minute, cme_second
        print *,  "Vr [m/s] = ", vr_cme
        print *, "Half width [deg]= ", w_half*180.0/dpi
        print *, "Co-latitude [deg] = ", clt_cme*180.0/dpi
        print *, "Longitude [deg]= ",  lon_cme*180.0/dpi
        print *, "Density [kg/m^3] = ", rho_cme
        print *, "Temperature [K] = ", temperature_cme
        print *, '=================================================================='
        print *, '=================================================================='//NEW_LINE('A')
      else
        print *, " No CME injected"//NEW_LINE('A')
      end if
    end if

    !Set conversion units
    unit_length=Lunit
    unit_density = Rhounit
    unit_time    = Tunit
    unit_velocity  = Vunit
    unit_magneticfield  = Bunit
    unit_pressure = Punit
    unit_numberdensity = unit_density / half / mp_SI
  end subroutine set_units

  subroutine print_initial_information()
    use mod_global_parameters

    if(mype==0) then
      print *, ''//NEW_LINE('A')
      print *, '=================================================================='
      print *, '=================================================================='//NEW_LINE('A')
      print *, 'Icarus'//NEW_LINE('A')
    end if
  end subroutine print_initial_information

  subroutine find_indices_coord_grid(coordinate, index_clt_1, index_lon_1, index_clt_2, index_lon_2)
    use mod_global_parameters
    double precision, dimension(3), intent(in)    :: coordinate
    integer, intent(out)    :: index_clt_1, index_lon_1, index_clt_2, index_lon_2
    integer                 :: counter_clt, counter_lon, nr_colat, nr_lon, j,k
    double precision        :: minimum

    nr_colat     = lt_3d(1)%n_points(1)
    nr_lon       = lt_3d(1)%n_points(2)

    counter_clt= 1
    minimum = abs(coord_grid_init(1,1,1)-coordinate(2))
    do j = 1,nr_colat
      if (abs(coord_grid_init(j,1,1)-coordinate(2)) < minimum ) then
        counter_clt = j
        minimum = abs(coord_grid_init(j,1,1)-coordinate(2))
      end if
    end do
    counter_lon= 1
    minimum = abs(coord_grid_init(1,1,2)-coordinate(3))
    do k = 1, nr_lon
      if (abs(coord_grid_init(1,k,2)-coordinate(3)) < minimum) then
        minimum = abs(coord_grid_init(1,k,2)-coordinate(3))
        counter_lon = k
      end if
    end do

    !SET THE COUNTERS FOR THE FOUR POINTS CONNECTED TO THE interpolation
    !see Q11,Q12,Q21,Q22 from wiki
    if (coordinate(2)>=coord_grid_init(counter_clt,counter_lon,1) .and. coordinate(3)>=coord_grid_init(counter_clt,counter_lon,2)) then
      index_clt_1 = counter_clt
      index_lon_1 = counter_lon
      index_clt_2 = counter_clt+1
      index_lon_2 = counter_lon+1
    end if
    if (coordinate(2)>=coord_grid_init(counter_clt,counter_lon,1) .and. coordinate(3)<coord_grid_init(counter_clt,counter_lon,2)) then
      index_clt_1 = counter_clt
      index_lon_1 = counter_lon-1
      index_clt_2 = counter_clt+1
      index_lon_2 = counter_lon
    end if
    if (coordinate(2)<coord_grid_init(counter_clt,counter_lon,1) .and. coordinate(3)>=coord_grid_init(counter_clt,counter_lon,2)) then
      index_clt_1 = counter_clt-1
      index_lon_1 = counter_lon
      index_clt_2 = counter_clt
      index_lon_2 = counter_lon+1
    end if
    if (coordinate(2)<coord_grid_init(counter_clt,counter_lon,1) .and. coordinate(3)<coord_grid_init(counter_clt,counter_lon,2)) then
      index_clt_1 = counter_clt-1
      index_lon_1 = counter_lon-1
      index_clt_2 = counter_clt
      index_lon_2 = counter_lon
    end if

  end subroutine find_indices_coord_grid

  subroutine find_trajectory_file(satellite_index, path_satellite_trajectories)
    use mod_global_parameters
    integer, intent(in)             :: satellite_index
    character(len=200) , intent(in) :: path_satellite_trajectories

    character(len=10), dimension(20) :: satellite_begin_dates = (/'1975_01_01', '1977_10_03', '1980_10_03', '1983_10_04', '1986_10_04', '1989_10_04', '1992_10_04', '1995_10_05', '1998_10_05', '2001_10_05', '2004_10_05', '2007_10_06', '2010_10_06', '2013_10_06', '2016_10_06', '2019_10_07', '2022_10_07', '2025_10_07', '2028_10_07', '2031_10_08'/)
    character(len=10), dimension(20) :: satellite_end_dates = (/'1978_04_01', '1981_04_01', '1984_04_01', '1987_04_02', '1990_04_02', '1993_04_02', '1996_04_02', '1999_04_03', '2002_04_03', '2005_04_03', '2008_04_03', '2011_04_04', '2014_04_04', '2017_04_04', '2020_04_04', '2023_04_05', '2026_04_05', '2029_04_05', '2032_04_05', '2034_12_31'/)
    character(len=10), dimension(9) :: sta_begin_dates = (/'2006_10_10', '2007_10_06', '2010_10_06', '2013_10_06', '2016_10_06', '2019_10_07', '2022_10_07', '2025_10_07', '2028_10_07'/)
    character(len=10), dimension(9) :: sta_end_dates = (/'2008_04_03', '2011_04_04', '2014_04_04', '2017_04_04', '2020_04_04', '2023_04_05', '2026_04_05', '2029_04_05', '2030_10_09'/)
    character(len=10), dimension(4) :: stb_begin_dates = (/'2006_10_10', '2007_10_06', '2010_10_06', '2013_10_06'/)
    character(len=10), dimension(4) :: stb_end_dates = (/'2008_04_03', '2011_04_04', '2014_04_04', '2016_09_12'/)
    character(len=10), dimension(3) :: psp_begin_dates = (/'2018_08_13', '2019_10_07', '2022_10_07'/)
    character(len=10), dimension(3) :: psp_end_dates = (/'2020_04_04', '2023_04_05', '2025_08_30'/)
    character(len=10), dimension(4) :: solo_begin_dates = (/'2020_02_11', '2022_10_07', '2025_10_07', '2028_10_07'/)
    character(len=10), dimension(4) :: solo_end_dates = (/'2023_04_05', '2026_04_05', '2029_04_05', '2030_11_17'/)
    character(len=10), dimension(3) :: bepi_begin_dates = (/'2018_10_21', '2019_10_03', '2022_10_02'/)
    character(len=10), dimension(3) :: bepi_end_dates = (/'2020_03_31', '2023_03_31', '2025_11_01'/)
    character(len=10), dimension(3) :: juno_begin_dates = (/'2011_08_06', '2013_10_03', '2016_10_03'/)
    character(len=10), dimension(3) :: juno_end_dates = (/'2014_04_01', '2017_04_01', '2019_06_19'/)
    character(len=10), dimension(:), allocatable :: begin_dates, end_dates

    character(len=11), dimension(10) :: satellite_list = (/'earth      ', 'mars       ', 'mercury    ', 'venus      ', 'sta        ', 'stb        ', 'psp_nom_R02', 'SolO       ', 'mpo        ', 'juno       '/)

    integer, dimension(10)           :: dates_lengths = (/20, 20, 20, 20, 9, 4, 3, 4, 3, 3/)
    integer :: begin_year, begin_month, begin_day, begin_year_previous, begin_month_previous
    integer :: first_year, first_month, first_day
    integer :: last_year, last_month, last_day

    integer :: AllocateStatus, DeAllocateStatus
    integer :: i, j, length

    length = dates_lengths(satellite_index)
    ALLOCATE(begin_dates(length), STAT = AllocateStatus)
    ALLOCATE(end_dates(length), STAT = AllocateStatus)

    if (satellite_index <= 4) then
      begin_dates = satellite_begin_dates
      end_dates = satellite_end_dates
    else if (satellite_index == 5) then
      begin_dates = sta_begin_dates
      end_dates = sta_end_dates
    else if (satellite_index == 6) then
      begin_dates = stb_begin_dates
      end_dates = stb_end_dates
    else if (satellite_index == 7) then
      begin_dates = psp_begin_dates
      end_dates = psp_end_dates
    else if (satellite_index == 8) then
      begin_dates = solo_begin_dates
      end_dates = solo_end_dates
    else if (satellite_index == 9) then
      begin_dates = bepi_begin_dates
      end_dates = bepi_end_dates
    else
      begin_dates = juno_begin_dates
      end_dates = juno_end_dates
    end if




    read(begin_dates(1)(1:4), '(i4)') first_year
    read(begin_dates(1)(6:7), '(i2)') first_month
    read(begin_dates(1)(9:10), '(i2)') first_day
    read(end_dates(length)(1:4), '(i4)') last_year
    read(end_dates(length)(6:7), '(i2)') last_month
    read(end_dates(length)(9:10), '(i2)') last_day

    if ((magnetogram_timestamp(1)>first_year .and. .not.(magnetogram_timestamp(1)==first_year+1 .and. magnetogram_timestamp(2)==1 .and. first_month==12) .or. magnetogram_timestamp(1)==first_year .and. magnetogram_timestamp(2)>first_month+1) .and. &
        (magnetogram_timestamp(1)<last_year .and. .not.(magnetogram_timestamp(1)==last_year-1 .and. magnetogram_timestamp(2)==12 .and. first_month==1) .or. magnetogram_timestamp(1)==last_year .and. magnetogram_timestamp(2)<last_month-1)) then

      do i=2, length
        read(begin_dates(i)(1:4), '(i4)') begin_year
        read(begin_dates(i)(6:7), '(i2)') begin_month
        read(begin_dates(i)(9:10), '(i2)') begin_day
        if (magnetogram_timestamp(1)<begin_year .or. magnetogram_timestamp(1)==begin_year .and. magnetogram_timestamp(2)<begin_month .or. &
            magnetogram_timestamp(1)==begin_year .and. magnetogram_timestamp(2)==begin_month .and. magnetogram_timestamp(3)<begin_day) then
          j = i-1
          if (i > 2) then
            if ((magnetogram_timestamp(1)==begin_year_previous .and. (magnetogram_timestamp(2)==begin_month_previous .or. magnetogram_timestamp(2)-1==begin_month_previous)) .or. &
                (magnetogram_timestamp(1)-1==begin_year_previous .and. magnetogram_timestamp(2)==1 .and. begin_month_previous==12)) then    ! if magnetogram time is too close to the begin date of the file (at most 1 month)
                j = i-2     ! change to the one file before
            end if
          end if
          trajectory_list(satellite_index) = trim(path_satellite_trajectories)//trim(satellite_list(satellite_index))//'__'//begin_dates(j)//'__'//end_dates(j)//'.unf'
          which_satellite(satellite_index) = 1
          exit
        end if
        begin_year_previous = begin_year
        begin_month_previous = begin_month
      end do

      if (which_satellite(satellite_index)==0) then
        if ((magnetogram_timestamp(1)==begin_year_previous .and. (magnetogram_timestamp(2)==begin_month_previous .or. magnetogram_timestamp(2)-1==begin_month_previous)) .or. &
            (magnetogram_timestamp(1)-1==begin_year_previous .and. magnetogram_timestamp(2)==1 .and. begin_month_previous==12)) then    ! if magnetogram time is too close to the begin date of the file (at most 1 month)
            length = length-1       ! change to the one file before
        end if
        trajectory_list(satellite_index) = trim(path_satellite_trajectories)//trim(satellite_list(satellite_index))//'__'//begin_dates(length)//'__'//end_dates(length)//'.unf'
        which_satellite(satellite_index) = 1
      end if
    end if

    DEALLOCATE(begin_dates, STAT = DeAllocateStatus)
    DEALLOCATE(end_dates, STAT = DeAllocateStatus)

  end subroutine find_trajectory_file

  subroutine read_satellite_trajectory(trajectory_file, index)
    use mod_global_parameters
    character(len=50), intent(in)   :: trajectory_file
    integer, intent(in)             :: index
    integer                         :: iUnit=16, iError, i, j
    integer                         :: arr_size(2), nr_positions, nr_coordinates
    double precision, allocatable   :: radii(:), latitudes(:), longitudes(:)
    integer, allocatable            :: year(:), month(:), day(:)
    integer, allocatable            :: hour(:), minute(:), second(:)
    integer                         :: AllocateStatus, DeAllocateStatus
    double precision                :: delta_time
    integer                         :: delta_steps, i_date, j_date, n

    open(iUnit, file=trajectory_file, action="read", form='unformatted', iostat=iError)
    if(iError /= 0) call mpistop('Importdata could not open real4 file = '//trim(trajectory_file))
    read(iUnit) arr_size
    nr_positions = arr_size(1)
    nr_coordinates = arr_size(2)

    ALLOCATE(year(nr_positions), STAT = AllocateStatus)
    ALLOCATE(month(nr_positions), STAT = AllocateStatus)
    ALLOCATE(day(nr_positions), STAT = AllocateStatus)
    ALLOCATE(hour(nr_positions), STAT = AllocateStatus)
    ALLOCATE(minute(nr_positions), STAT = AllocateStatus)
    ALLOCATE(second(nr_positions), STAT = AllocateStatus)
    ALLOCATE(radii(nr_positions), STAT = AllocateStatus)
    ALLOCATE(latitudes(nr_positions), STAT = AllocateStatus)
    ALLOCATE(longitudes(nr_positions), STAT = AllocateStatus)

    read(iUnit) year
    read(iUnit) month
    read(iUnit) day
    read(iUnit) hour
    read(iUnit) minute
    read(iUnit) second
    read(iUnit) radii
    read(iUnit) latitudes
    read(iUnit) longitudes
    close(iUnit)

    ALLOCATE(positions_list(index)%positions(nr_coordinates, nr_positions), STAT = AllocateStatus)

    positions_list(index)%positions(1,:) = year
    positions_list(index)%positions(2,:) = month
    positions_list(index)%positions(3,:) = day
    positions_list(index)%positions(4,:) = hour
    positions_list(index)%positions(5,:) = minute
    positions_list(index)%positions(6,:) = second
    positions_list(index)%positions(7,:) = radii
    positions_list(index)%positions(8,:) = latitudes
    positions_list(index)%positions(9,:) = longitudes

    delta_time = 0.25d0
    if (magnetogram_index(index) .eq. 0) then
      do j_date = 1, size(year)
        if ((year(j_date) == magnetogram_timestamp(1)) .and.  (month(j_date) == magnetogram_timestamp(2))) then
          if ((day(j_date) == magnetogram_timestamp(3)) .and. (hour(j_date) == magnetogram_timestamp(4))) then
            if (minute(j_date) == magnetogram_timestamp(5)) then
              magnetogram_index(index) = j_date
              exit
            end if
          end if
        end if
      end do
    end if

    if (magnetogram_index(index) == 0) then
       ! pick nearest minute (simple fallback: first index)
       magnetogram_index(index) = 1
    end if

    if (starting_index(index, 1) == 0 .and. (num_cmes == 0)) then
      starting_index(index, 1) = max(1,magnetogram_index(index))
      cme_index(index,1) = magnetogram_index(index)
    end if


    if (starting_index(index, 1) .eq. 0 .and. (num_cmes > 0)) then
      do n = 1, num_cmes
        do i_date = 1, size(year)
          if ((year(i_date) == cme_year(n)) .and.  (month(i_date) == cme_month(n))) then
            if ((day(i_date) == cme_day(n)) .and. (hour(i_date) == cme_hour(n))) then
              if (minute(i_date) == cme_minute(n)) then
                cme_index(index, n) = i_date
                time_difference_cme_magn(index, n) = (cme_index(index, n) - magnetogram_index(index))/60.0 !hours

                delta_steps = int((relaxation*24.0+cme_insertion*24.0+time_difference_cme_magn(index, n))*60)
                starting_index(index, n) = i_date - delta_steps
                starting_index(index, n) = max(1, min(nr_positions, starting_index(index, n))) 

                exit
              end if
            end if
          end if
        end do
      end do
    end if

    DEALLOCATE(year, STAT = DEAllocateStatus)
    DEALLOCATE(month, STAT = DEAllocateStatus)
    DEALLOCATE(day, STAT = DEAllocateStatus)
    DEALLOCATE(hour, STAT = DEAllocateStatus)
    DEALLOCATE(minute, STAT = DEAllocateStatus)
    DEALLOCATE(second, STAT = DEAllocateStatus)
    DEALLOCATE(radii, STAT = DEAllocateStatus)
    DEALLOCATE(latitudes, STAT = DEAllocateStatus)
    DEALLOCATE(longitudes, STAT = DEAllocateStatus)


  end subroutine read_satellite_trajectory

  subroutine read_cme_parameters(filename)
    use mod_global_parameters
    character(len=50), intent(in)   :: filename
    integer                         :: iUnit=30, iError, i, DEAllocateStatus, j
    character(len=50)               :: cme_parameter_file
    character(len=500)              :: commented_line

    open(iUnit, file = filename, status = 'old', action="read", iostat=iError)
    if(iError /= 0) call mpistop('Importdata could not open real4 file = '//trim(filename))

    ALLOCATE(cme_type(num_cmes))
    ALLOCATE(cme_date(num_cmes))
    ALLOCATE(cme_year(num_cmes))
    ALLOCATE(cme_month(num_cmes))
    ALLOCATE(cme_day(num_cmes))
    ALLOCATE(cme_hour(num_cmes))
    ALLOCATE(cme_minute(num_cmes))
    ALLOCATE(cme_second(num_cmes))

    ALLOCATE(timestamp(num_cmes))
    ALLOCATE(vr_cme(num_cmes))
    ALLOCATE(w_half(num_cmes))
    ALLOCATE(clt_cme(num_cmes))
    ALLOCATE(lon_cme(num_cmes))
    ALLOCATE(rho_cme(num_cmes))
    ALLOCATE(temperature_cme(num_cmes))
    ALLOCATE(cme_index(10,num_cmes))
    ALLOCATE(starting_index(10, num_cmes))
    ALLOCATE(time_difference_cme_magn(10, num_cmes))
    ALLOCATE(longitudes_fix(num_cmes))

    do i=1, num_cmes
      read(iUnit,*) cme_type(i), cme_date(i), clt_cme(i), lon_cme(i), w_half(i),  vr_cme(i), rho_cme(i), temperature_cme(i)
      w_half(i) = w_half(i) * dpi/180.0
      clt_cme(i) = (-clt_cme(i) + 90.0) * dpi/180.0
      lon_cme(i) = lon_cme(i)*dpi/180.0
      vr_cme(i) = vr_cme(i)*1000.0
    end do

    close(iUnit)
    do i=1, num_cmes
      read (cme_date(i)(1:4),*) cme_year(i)
      read (cme_date(i)(6:7),*) cme_month(i)
      read (cme_date(i)(9:10),*) cme_day(i)
      read (cme_date(i)(12:13),*) cme_hour(i)
      read (cme_date(i)(15:16),*) cme_minute(i)
      read (cme_date(i)(18:19),*) cme_second(i)
    end do
  end subroutine read_cme_parameters

! apply frame correction to CME longitudes using the simulation frame rate
  subroutine cme_insertion_longitudes_fix()
    use mod_global_parameters
    implicit none
    integer :: i
    double precision :: time_shift  ! hours between magnetogram epoch and CME timestamp

    if (num_cmes <= 0) return

    do i = 1, num_cmes
      ! Prefer the precomputed difference if available (hours):
      if (allocated(time_difference_cme_magn)) then
        time_shift = time_difference_cme_magn(1, i)
      else
        ! Fallback: recover the offset from the fields you already have (hours)
        time_shift  = timestamp(i) - (relaxation + cme_insertion) * 24.0d0
      end if

! phi_rot = phi_heeq +  (omega_frame - omega_HEEQ)*t0

      ! Frame-rotation correction (rad): how much the corotating frame turned over time_shift
      longitudes_fix(i) = (omega_frame-omega_HEEQ) * time_shift

      ! Express CME longitude at the simulation’s reference epoch (subtract forward rotation)
      lon_cme(i) = lon_cme(i) - longitudes_fix(i)

      ! Wrap to [0, 2π)
      lon_cme(i) = modulo(lon_cme(i), 2.0d0*dpi)
      if (lon_cme(i) < 0.d0) lon_cme(i) = lon_cme(i) + 2.0d0*dpi
    end do
  end subroutine cme_insertion_longitudes_fix

! half time (scalar) for CME i
  subroutine find_half_time(t_half, i)
    use mod_global_parameters
    implicit none
    double precision, intent(out) :: t_half
    integer,          intent(in)  :: i
    double precision              :: au

    au      = 1.49d11
    ! Eq. 2 (Shapes paper), evaluated at 0.1 AU; return in simulation time units
    t_half  = 0.1d0 * sin(w_half(i)) * au / vr_cme(i) / unit_time
  end subroutine find_half_time

! opening half-angle (scalar) for CME i; spheroidal (planar-b) model
  subroutine find_opening_angle_spherical(t_half, theta_half, i)
    use mod_global_parameters
    implicit none
    double precision, intent(in)  :: t_half
    double precision, intent(out) :: theta_half
    integer,          intent(in)  :: i
    double precision              :: tau, arg

    theta_half = -1.0d0
    if (timestamp(i) <= global_time .and. global_time <= timestamp(i) + 2.0d0*t_half) then
      tau = (global_time - timestamp(i)) / t_half - 1.0d0
      arg = 1.0d0 - tau*tau
      if (arg > 0.0d0) then
        theta_half = w_half(i) * sqrt(arg)
      end if
    end if
  end subroutine find_opening_angle_spherical


! opening half-angle (scalar) using Eq. 5 (spheroidal in planar-b)
  subroutine find_opening_angle_spherical3(t_half, theta_half, i)
    use mod_global_parameters
    implicit none
    double precision, intent(in)  :: t_half
    double precision, intent(out) :: theta_half
    integer,          intent(in)  :: i
    double precision              :: au, r_half, r_new, num, den, cval

    theta_half = -1.0d0
    if (timestamp(i) <= global_time .and. global_time <= timestamp(i) + 2.0d0*t_half) then
      au     = 1.49d11
      r_half = 0.1d0 * au * sin(w_half(i))
      r_new  = (global_time - timestamp(i)) * vr_cme(i) * unit_time + 0.1d0*au - r_half
      if (r_new > 0.0d0) then
        num  = r_new**2 + (0.1d0*au)**2 - r_half**2
        den  = 2.0d0 * r_new * 0.1d0 * au
        cval = max(-1.0d0, min(1.0d0, num/den))    ! clamp for acos safety
        theta_half = acos(cval)
      end if
    end if
  end subroutine find_opening_angle_spherical3

! point-in-CME test (great-circle distance with wrap-around)
  subroutine mask(clt_point, lon_point, is_inside_cme, i)
    use mod_global_parameters
    implicit none
    double precision, intent(in)  :: clt_point, lon_point    ! co-lat (θ), lon (φ)
    integer,          intent(out) :: is_inside_cme
    integer,          intent(in)  :: i

    double precision :: half_time, half_angle
    double precision :: dphi, cosd, d

    is_inside_cme = 0
    half_angle    = -1.0d0

    call find_half_time(half_time, i)
    call find_opening_angle_spherical3(half_time, half_angle, i)

    if (half_angle <= 0.d0) return

    ! Wrap-aware longitude difference Δφ ∈ (-π, π]
    dphi = modulo((lon_point - lon_cme(i)) + dpi, 2.d0*dpi) - dpi

    ! Spherical law of cosines for great-circle distance:
    ! cos d = sinθ sinθc cosΔφ + cosθ cosθc
    cosd = sin(clt_point)*sin(clt_cme(i))*cos(dphi) + cos(clt_point)*cos(clt_cme(i))
    cosd = max(-1.d0, min(1.d0, cosd))
    d    = acos(cosd)

    if (d < half_angle) is_inside_cme = 1
  end subroutine mask
end module mod_usr
