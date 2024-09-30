module mod_usr
  use mod_mhd
  use mod_lookup_table
  use mod_bc_data
  !SI units, constants
  use mod_constants, only: mp_SI, kB_SI, miu0_SI
  implicit none

  character(len=20)                              :: printsettingformat
  double precision                               :: omega_frame
  character(len=500)                             :: amr_criterion, cme_parameter_file

  integer                                        :: cme_flag, num_cmes, relaxation, cme_insertion
  type satellite_pos
    real(kind=8), dimension(:,:), allocatable    :: positions
  end type satellite_pos

  !Shared over subroutines
  real(kind=8), allocatable                      :: coord_grid_init(:,:,:),variables_init(:,:,:)
  type(satellite_pos), dimension(:), allocatable :: positions_list
  character(len=250), dimension(8)               :: trajectory_list
  integer, dimension(8)                          :: which_satellite = (/1, 1,1, 1, 1, 1, 0, 0/)     ! intended order: earth, mars, mercury, venus, sta, stb, psp, solo
  integer, dimension(8)                          :: sat_indx = (/0, 0, 0, 0, 0, 0, 0, 0/)     ! intended order: earth, mars, mercury, venus, sta, stb, psp, solo
  integer                                        :: sat_count=0, zero_count=0
  double precision, dimension(8)                 :: last = (/0, 0, 0, 0, 0, 0, 0, 0/)        ! intended order: earth, mars, mercury, venus, sta, stb, psp, solo
  integer, dimension(8)                          :: last_index = (/0, 0, 0, 0, 0, 0, 0, 0/)  ! intended order: earth, mars, mercury, venus, sta, stb, psp, solo
  integer, dimension(8)                          :: last_index_s = (/0, 0, 0, 0, 0, 0, 0, 0/)
  double precision                               :: delta_phi
  integer, dimension(:,:), allocatable           :: starting_index, cme_index    ! first coordinate: satellite index; second coordinate: cme index
  integer, dimension(8)                          :: magnetogram_index = (/0, 0, 0, 0, 0, 0, 0, 0/)

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
    character(len=*), intent(in) :: files(:)
    integer :: n

    namelist /rotating_frame_list/ omega_frame
    namelist /icarus_list/ amr_criterion, cme_flag, num_cmes, relaxation, cme_insertion, &
    cme_parameter_file, magnetogram_time

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, rotating_frame_list, end=111)

       read(unitpar, icarus_list, end=111)
111    close(unitpar)
    end do
   if (num_cmes .eq.0) then
      cme_flag = 0
   end if
   if (cme_flag .eq. 0) then
     cme_insertion = 0
   end if
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
    !  Sometime we use He_abundance, H_ion_fr, He_ion_fr, He_ion_fr2
    !  Moreover, we use 3 out of
    !      unit_density, unit_numberdensity, unit_length, unit_time, unit_velocity,
    !      unit_pressure, unit_magneticfield, unit_temperature,
    !      unit_mass, unit_charge
    !  Note we have factor RR (p=RR rho T)
    printsettingformat='(1x,A50,ES15.7)'
    if(mype==0) then
      write(*,*)'----------------PARAMETERS--   ----------------------'
      write(*,printsettingformat) "mhd_gamma ",mhd_gamma
      write(*,printsettingformat) "mhd_eta ",mhd_eta
      write(*,*)'----------------BEGIN UNITS  ------------------------'
      write(*,*)'----------------UNIT CONTROLS------------------------'
      write(*,*) "SI_unit",SI_unit
      write(*,*) "eq_state_units",eq_state_units
      write(*,*) "mhd_partial_ionization",mhd_partial_ionization
      write(*,printsettingformat) "He_abundance",He_abundance
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



      ! WARNING: ASSUMES STENCIL IS USING 2 GHOSTCELLS: may need to GENERALIZE!
      !We have 4 extra points for longitude because the boundary is periodic
      ALLOCATE(coord_grid_init(nr_colat, nr_lon+4, ndim-1), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) call mpistop('*** Not enough memory ***')
      ALLOCATE(variables_init(nr_colat, nr_lon+4, 4), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) call mpistop('*** Not enough memory ***')



      ! read in cme parameters
      if (num_cmes == 0) then
        ALLOCATE(timestamp(1))
        ALLOCATE(cme_index(8,1))
        ALLOCATE(starting_index(8, 1))
        ALLOCATE(time_difference_cme_magn(8, 1))
      else
        call read_cme_parameters(cme_parameter_file)
      end if
      ! Initialize cme starting index in the trajectory file, cme index in the trajectory file and the time difference between the start and cme indexes
      if (num_cmes > 0) then
        do n = 1, num_cmes
          do i = 1, 8
            cme_index(i, n) = 0
            starting_index(i, n) = 0
            time_difference_cme_magn(i, n) = 0
          end do
        end do
        else
          do i = 1, 8
            cme_index(i, 1) = 0
            starting_index(i, 1) = 0
            time_difference_cme_magn(i, 1) = 0
          end do
       end if



        earth_trajectory = "orbit/2015_june_earth_ext.unf"
        mars_trajectory = "orbit/2015_june_mars_ext.unf"
        venus_trajectory = "orbit/2015_june_venus_ext.unf"
        mercury_trajectory = "orbit/2015_june_mercury_ext.unf"
        sta_trajectory = "orbit/2015_june_sta_ext.unf"
        stb_trajectory = "orbit/2015_june_stb_ext.unf"


        trajectory_list(1) = earth_trajectory
        trajectory_list(2) = mars_trajectory
        trajectory_list(4) = venus_trajectory
        trajectory_list(3) = mercury_trajectory
        trajectory_list(5) = sta_trajectory
        trajectory_list(6) = stb_trajectory

        ALLOCATE(positions_list(8), STAT=AllocateStatus)

      ! for each satellite, read the trajectory data and save in the arrays of time and locations
      do i = 1, 8
        if (which_satellite(i)==1) then
          sat_indx(i-zero_count) = i
          sat_count = sat_count+1
          call read_satellite_trajectory(trajectory_list(i), i)
        end if
         if (which_satellite(i) == 0) then
          zero_count = zero_count+1
         end if
      end do

      ! calculate timestamp for cme insertion
      timestamp(:) = relaxation*24.0+cme_insertion*24.0
      if (num_cmes >0) then
        do k=1, num_cmes
          timestamp(k) = timestamp(k) + time_difference_cme_magn(1, k)
        end do
        call cme_insertion_longitudes_fix()
      end if

    end if

    mhd_gamma= 3.0d0/2.0d0
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
    double precision, intent(out) :: x(3, n_particles)
    double precision, intent(out) :: v(3, n_particles)
    double precision, intent(out) :: q(n_particles)
    double precision, intent(out) :: m(n_particles)
    logical, intent(out)          :: follow(n_particles)
    integer                       :: satellite_index, delta_sat

    if (sat_count < n_particles) then
      delta_sat = n_particles-sat_count
    end if
    do satellite_index = 1, n_particles-delta_sat
      v(:, satellite_index) = 0.d0
      q(satellite_index) = 0.d0
      m(satellite_index) = 0.d0
      call get_particle(x(:, sat_indx(satellite_index)), sat_indx(satellite_index), n_particles-delta_sat)
    end do
    follow(:) = .true.

  end subroutine generate_particles

  subroutine get_particle(x, satellite_index, n_particles)
    double precision, intent(out)      :: x(3)
    integer, intent(in)                :: satellite_index
    integer                            :: n_particles
    double precision, dimension(8)     :: orbital_period = (/365.24, 686.98, 87.969, 224.7, 346.0, 388.0, 88.0, 168.0/)    ! earth, mars, mercury, venus, sta, stb, psp, solo
    double precision                   :: phi_satellite, before_cme


    before_cme = (cme_index(1,1) - magnetogram_index(1))/60.0
    x(1) = positions_list(satellite_index)%positions(7, starting_index(satellite_index,1))
    x(2) = (dpi/2.0 - positions_list(satellite_index)%positions(8, starting_index(satellite_index,1)))

    ! phi_satellite here is at qt = 0, so at the simulation start
    phi_satellite = delta_phi+positions_list(satellite_index)%positions(9, starting_index(satellite_index,1))&
     + ((timestamp(1)-before_cme))*(2.0*dpi)/24.0*(1/2.447d1-1/orbital_period(1))

    if (phi_satellite < 0) then
      phi_satellite = 2 * dpi + phi_satellite
    else if (phi_satellite > 2*dpi) then
      phi_satellite = mod(phi_satellite, 2.0*dpi)
    end if
    x(3) = phi_satellite
  end subroutine get_particle

  subroutine move_particle(x, satellite_index, told, tnew)
    double precision, intent(inout) :: x(3)
    double precision, intent(in)    :: told,tnew
    double precision                :: xf(3), xc(3), x_test(3)
    integer, intent(in)             :: satellite_index

    double precision                :: phi_satellite, before_cme, delta_lon, lon_old, lon_new
    double precision, dimension(8)     :: orbital_period = (/365.24, 686.98, 87.969, 224.7, 346.0, 388.0, 88.0, 168.0/)    ! earth, mars, mercury, venus, sta, stb, psp, solo
    double precision                :: curr_lon, final_fix

    last_index_s(satellite_index) = starting_index(satellite_index, 1) + floor(tnew*60.0)
    before_cme = (cme_index(1,1) - magnetogram_index(1))/60.0

    xf(1) = positions_list(satellite_index)%positions(7, last_index_s(satellite_index))
    xf(2) = dpi/2.0 - positions_list(satellite_index)%positions(8, last_index_s(satellite_index))
    xf(3) = x(3) - (tnew-told)*(2.0*dpi)/24.0*(1/2.447d1-1/orbital_period(1))
    last_index_s(satellite_index) = starting_index(satellite_index, 1) + ceiling(tnew*60.0)

    xc(1) = positions_list(satellite_index)%positions(7, last_index_s(satellite_index))
    xc(2) = dpi/2.0 - positions_list(satellite_index)%positions(8, last_index_s(satellite_index))
    xc(3) = x(3) - (tnew-told)*(2.0*dpi)/24.0*(1/2.447d1-1/orbital_period(1))

    if ((ceiling(tnew) - floor(tnew)) .gt. 0.0) then
      x(1) = xf(1) + (tnew - floor(tnew))*(xc(1)-xf(1))/(ceiling(tnew) - floor(tnew))
      x(2) = xf(2) + (tnew - floor(tnew))*(xc(2)-xf(2))/(ceiling(tnew) - floor(tnew))
      x(3) = xf(3) + (tnew - floor(tnew))*(xc(3)-xf(3))/(ceiling(tnew) - floor(tnew))
    else
      x(1) = xf(1)
      x(2) = xf(2)
      x(3) = xf(3)
    end if


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

    double precision :: velocity2d(ixmin2:ixmax2,ixmin3:ixmax3)
    double precision :: rho2d(ixmin2:ixmax2,ixmin3:ixmax3)
    double precision :: p2d(ixmin2:ixmax2,ixmin3:ixmax3)
    double precision :: br2d(ixmin2:ixmax2,ixmin3:ixmax3)
    integer ::  idir, i


    w(ix^S,1:nw) = zero


    r_boundary   = xprobmin1 !in R_sun

    velocity2d(ixmin2:ixmax2, ixmin3:ixmax3) = bc_data_get_3d(bc_data_ix(mom(1), 1), &
           x(ixmin1, ixmin2:ixmax2, ixmin3:ixmax3, 2), &
           x(ixmin1, ixmin2:ixmax2, ixmin3:ixmax3, 3), 0d0)
    rho2d(ixmin2:ixmax2, ixmin3:ixmax3) = bc_data_get_3d(bc_data_ix(rho_, 1), &
           x(ixmin1, ixmin2:ixmax2, ixmin3:ixmax3, 2), &
           x(ixmin1, ixmin2:ixmax2, ixmin3:ixmax3, 3), 0d0)
    p2d(ixmin2:ixmax2, ixmin3:ixmax3) = bc_data_get_3d(bc_data_ix(p_, 1), &
           x(ixmin1, ixmin2:ixmax2, ixmin3:ixmax3, 2), &
           x(ixmin1, ixmin2:ixmax2, ixmin3:ixmax3, 3), 0d0)
    br2d(ixmin2:ixmax2, ixmin3:ixmax3) = bc_data_get_3d(bc_data_ix(mag(1), 1), &
           x(ixmin1, ixmin2:ixmax2, ixmin3:ixmax3, 2), &
           x(ixmin1, ixmin2:ixmax2, ixmin3:ixmax3, 3), 0d0)



   do ix1=ixmin1,ixmax1
      w(ix1,ixmin2:ixmax2,ixmin3:ixmax3,mom(1))=velocity2d(ixmin2:ixmax2, ixmin3:ixmax3)
      w(ix1,ixmin2:ixmax2,ixmin3:ixmax3,rho_)=rho2d(ixmin2:ixmax2, ixmin3:ixmax3)&
      *(r_boundary/x(ix1,ixmin2:ixmax2, ixmin3:ixmax3, 1))**2
      w(ix1,ixmin2:ixmax2,ixmin3:ixmax3,p_)=p2d(ixmin2:ixmax2, ixmin3:ixmax3)&
      *(r_boundary/x(ix1,ixmin2:ixmax2, ixmin3:ixmax3, 1))**2
      w(ix1,ixmin2:ixmax2,ixmin3:ixmax3,mag(1))=br2d(ixmin2:ixmax2, ixmin3:ixmax3)&
      *(r_boundary/x(ix1,ixmin2:ixmax2, ixmin3:ixmax3, 1))**2
    enddo

    !Convert to conserved values
    call mhd_to_conserved(ixG^L,ix^L,w,x)

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

    ! output divB1
    call get_divb(w,ixI^L,ixO^L,divb)
    w(ixO^S,nw+1)=divb(ixO^S)

    do i=1,ndir
      v(ixI^S,i)=w(ixI^S,mom(i))/w(ixI^S,rho_)
    end do

    call divvector(v,ixI^L,ixO^L,divV)
    !w(ixO^S,nw+2)=divV(ixO^S)*step_size(ixI^S)
    w(ixO^S,nw+2)=divV(ixO^S)*1.37
    do i=1,ndir
      momentum(ixI^S,i)=w(ixI^S,mom(i))
    end do

    call divvector(momentum,ixI^L,ixO^L,divmom)
    w(ixO^S,nw+3)=divmom(ixO^S)

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
    ! refine = -1 enforce to not refine
    ! refine =  0 doesn't enforce anything
    ! refine =  1 enforce refinement
    ! coarsen = -1 enforce to not coarsen
    !  coarsen =  0 doesn't enforce anything
    ! coarsen =  1 enforce coarsen
    integer, intent(inout)          :: refine, coarsen

    ! added by me
    integer, intent(in)             :: igrid, level, ixI^L, ixO^L
    double precision, intent(in)    :: qt, w(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision                :: v(ixI^S,ndir), divV(ixI^S)
    integer                         :: i, block_num
    double precision                :: threshold
    integer                         :: ix1, ix2, ix3
    double precision                :: xloc(1:ndim)
    double precision                :: r, theta, phi, r_this, theta_this, phi_this
    double precision                :: lon_cir, u_artificial
    double precision                :: before_cme
    double precision                :: phi_satellite

    ! To Follow Earth location in the domain

    before_cme = (cme_index(1,1) - magnetogram_index(1))/60.0
    phi_satellite = delta_phi + positions_list(1)%positions(9, last_index(1))&
      - (qt-(timestamp(1)-before_cme))*(2.0*dpi)/24.0*(1/2.447d1-1/365.24)

    ! Threshold for negative nabla V
    threshold = -0.005

    ! Refinement criterion for div(V)
    if (amr_criterion == "shock") then
      if (qt > timestamp(1)) then
        do i=1,ndir
          v(ixI^S,i)=w(ixI^S,mom(i))/w(ixI^S, rho_)
        end do
        call divvector(v,ixI^L,ixO^L,divV)

        if (any(divV(ixO^S) < threshold)) then
          refine = 1
          coarsen = -1
        else
          refine = -1
          coarsen = 1
        end if
      end if
    end if

   if (amr_criterion == "lonwindow" .and. qt>360.) then
    if (any(x(ixI^S,3) >= dpi/2) .and. any(x(ixI^S,3) <= dpi/2)) then
      refine = 1
      coarsen = -1
      else
      refine = -1
      coarsen = 1
    end if
   end if

    ! Refinement criterion for tracing function
    if (amr_criterion == "tracing") then
      if (qt > timestamp(1)) then
      !  if (any(w(ixI^S,tracer(1)) > 0.001) .and. any(x(ixI^S,1) >= 50.0)) then
          if (any(w(ixI^S,tracer(1)) > 0.001)) then
            if (any(x(ixI^S,3) < phi_satellite + 30*dpi/180) .and. any(x(ixI^S,3) > max(phi_satellite -30*dpi/180,0.0)) ) then
              refine = 1
              coarsen = -1
            else
              refine = -1
              coarsen = 1
            end if
        else
          refine = -1
          coarsen = 1
        end if

      end if
    end if

  end subroutine specialrefine_grid

  subroutine specialbound_usr(qt,ixI^L,ixO^L,iB,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L, iB
    double precision, intent(in)    :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    integer             :: ix3,ix2
    double precision    :: r, theta, phi, minimum
    integer             :: i,j,k, valuej, valuek
    integer             :: nr_r, nr_colat, nr_lon
    integer             :: point11_clt, point11_lon, point22_clt, point22_lon
    double precision    :: xloc(1:ndim)

    double precision    :: clt_zero, lon_zero
    integer             :: mask_cme, n, local_check

    select case(iB)
      case(1)! Lower radial boundary
        nr_colat = SIZE(coord_grid_init(:,1,1))
        nr_lon = SIZE(coord_grid_init(1,:,1))
!        w(ixO^S,:) = 0.d0

        !W IS USING CONSERVATIVE VALUES
        !SO TO GET v1 WE NEED TO USE m1/rho!!

        !rghost1 = x(ixOmax1,1,1,r_): radial coordinate at first ghost cell (closest to in-domain)
        !rghost2 = x(ixOmin1,1,1,r_): radial coordinate at second ghost cell
        !Inner boundary has a linear relation: f(r) = a + b*r
        !We have f(r) boundary 21.5 and r = 21.5 and first physical cell f(r)
        !then b. = f(r)_bound - f(r)_phys / (R_bound - R_phys)
        !and a = f(r)_bound - b*21.5R_sun
        !FOR NOW WE ASSUME NON-STRETCHED GRID!
        call mhd_to_primitive(ixI^L,ixO^L,w,x)
        do ix3 = ixOmin3, ixOmax3
          do ix2 = ixOmin2,ixOmax2
            !Get data from the first phsyical cell going in the r-direction
            r       = x(ixOmin1 +2 , ix2, ix3, 1)
            theta   = x(ixOmin1, ix2, ix3, 2)
            phi     = x(ixOmin1, ix2, ix3, 3)
            xloc(:) = x(ixOmin1 +2, ix2, ix3,:)
            !IMPORTANT NOTE:
!            !AMRVAC has ghostcells that are largers than 2*pi
            !We made sure that coord_grid_init also has data >2*pi

            !Find the corresponding data in the initial grid
!            !For now we take the values at 21.5+delta_r/2 as our grid coord_grid_init does not have 21.5R_sun
            !TODO: add the 21.5 R_sun grid!!

            call find_indices_coord_grid(xloc, point11_clt, point11_lon, point22_clt, point22_lon)
            mask_cme = 0
            if (num_cmes == 0) then
              local_check = 1
            else
              local_check = num_cmes
            end if
            do n=1, local_check
              if (mask_cme .eq. 0) then
                if (cme_flag == 1) then
                  call mask(xloc(2), xloc(3), mask_cme, n)
                end if
                if (mask_cme .eq. 1) then
                  ! Mass density
                  w(ixOmax1, ix2, ix3, rho_) = rho_cme(n)/unit_density
                  w(ixOmin1, ix2, ix3, rho_) = rho_cme(n)/unit_density
                  ! speed
                  w(ixOmax1, ix2, ix3, mom(1)) =  vr_cme(n)/unit_velocity
                  w(ixOmin1, ix2, ix3, mom(1)) =  vr_cme(n)/unit_velocity
                  !pressure: p
                  w(ixOmax1, ix2, ix3, p_) = rho_cme(n) / (0.5 * mp_SI) * kB_SI * temperature_cme(n)/unit_pressure
                  w(ixOmin1, ix2, ix3, p_) = rho_cme(n) / (0.5 * mp_SI) * kB_SI * temperature_cme(n)/unit_pressure
                  ! Setting tracer function to the value of densicy inside CME
                  w(ixOmax1, ix2, ix3,  tracer(1)) = rho_cme(n)/unit_density
                  w(ixOmin1, ix2, ix3,  tracer(1)) = rho_cme(n)/unit_density
                 else
                  w(ixOmax1, ix2, ix3,  tracer(1)) = - w(ixOmin1+2, ix2, ix3,  tracer(1))
                  w(ixOmin1, ix2, ix3,  tracer(1)) = - w(ixOmin1+3, ix2, ix3,  tracer(1))
                end if
              end if
             end do
            !Mass density
      end do
    end do
        !convert back to conserved values
        call mhd_to_conserved(ixI^L,ixO^L,w,x)
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
      print *, 'EUHFORIA'//NEW_LINE('A')
    end if
  end subroutine print_initial_information

  subroutine find_indices_coord_grid(coordinate, index_clt_1, index_lon_1, index_clt_2, index_lon_2)
    use mod_global_parameters
    double precision, dimension(3), intent(in)    :: coordinate
    integer, intent(out)    :: index_clt_1, index_lon_1, index_clt_2, index_lon_2
    integer                 :: counter_clt, counter_lon, nr_colat, nr_lon, j,k
    double precision        :: minimum

    nr_colat     = SIZE(coord_grid_init(:,1,1))
    nr_lon       = SIZE(coord_grid_init(1,:,1))

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
    character(len=10), dimension(:), allocatable :: begin_dates, end_dates

    character(len=11), dimension(8) :: satellite_list = (/'earth      ', 'mars       ', 'mercury    ', 'venus      ', 'sta        ', 'stb        ', 'psp_nom_R02', 'SolO       '/)

    integer, dimension(8)           :: dates_lengths = (/20, 20, 20, 20, 9, 4, 3, 4/)
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
    else
      begin_dates = solo_begin_dates
      end_dates = solo_end_dates
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
    ! delta_steps = int(timestamp*60)

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

    if (starting_index(index, 1) .eq. 0 .and. (num_cmes == 0)) then
      starting_index(index, 1) = magnetogram_index(index)
      !time_difference_cme_magn(index, 1) = 0.0
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
    ALLOCATE(cme_index(8,num_cmes))
    ALLOCATE(starting_index(8, num_cmes))
    ALLOCATE(time_difference_cme_magn(8, num_cmes))
    ALLOCATE(longitudes_fix(num_cmes))

    do i=1, num_cmes
      read(iUnit,*) cme_type(i), cme_date(i), clt_cme(i), lon_cme(i), w_half(i),  vr_cme(i), rho_cme(i), temperature_cme(i)
      w_half(i) = w_half(i) * dpi/180.0
      clt_cme(i) = (-clt_cme(i) + 90.0) * dpi/180.0
      lon_cme(i) = delta_phi + lon_cme(i) * dpi/180.0
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

  subroutine cme_insertion_longitudes_fix()
    integer       :: i

    do i=1, num_cmes
       longitudes_fix(i) = (timestamp(i)-relaxation*24.0 - cme_insertion*24.0)*(2.0*dpi)/24.0*(1/2.447d1-1/365.24)
       lon_cme(i) = lon_cme(i) - longitudes_fix(i)
    end do

  end subroutine cme_insertion_longitudes_fix

  subroutine find_half_time(t_half, i)
    use mod_global_parameters
    double precision, dimension(num_cmes)          :: t_half
    double precision                               :: au = 1.49d11
    integer                                        :: i
    !	 	we need to calculate at 0.1 AU

    ! t_half defined as equation 2 in Shapes paper
     t_half(i) = 0.1 * sin(w_half(i)) * au / vr_cme(i)/unit_time
    ! t_half defined as equation 1 in Shapes paper
    !  t_half(i) = 0.1 * tan(w_half(i)) * au / vr_cme(i) /unit_time
     ! t_half(i) = t_half(i)/unit_time

  end subroutine find_half_time

  subroutine find_opening_angle_spherical(t_half, theta_opening_angle, i)
    use mod_global_parameters
    double precision, dimension(num_cmes)        :: t_half,  theta_opening_angle
    integer                                      :: i

    if (timestamp(i) <= global_time) then
      if (global_time <= (timestamp(i) + 2*t_half(i))) then
        ! opening angle defined as eq 4 in Shapes paper. spheroidal in planar b
        theta_opening_angle(i) = w_half(i) * sqrt(1-((global_time-timestamp(i)) / t_half(i) - 1)**2)
      end if
    end if

  end subroutine find_opening_angle_spherical

  subroutine find_opening_angle_spherical3(t_half, theta_opening_angle, i)
    use mod_global_parameters
    double precision, dimension(num_cmes)    :: t_half, r_half,  theta_opening_angle, r_new
    double precision           :: au = 1.49d11
    integer  :: i

    if (timestamp(i) <= global_time) then
      if (global_time <= (timestamp(i) + 2*t_half(i))) then
        ! opening angle defined as eq 5 in Shapes paper. spheroidal in planar b
        r_half(i) = 0.1 * au * sin(w_half(i))
        r_new(i) = (global_time - timestamp(i)) * vr_cme(i)*unit_time + 0.1*au - r_half(i)

        theta_opening_angle(i)=acos((r_new(i)**2+(0.1*au)**2-r_half(i)**2)/(2.0*r_new(i)*0.1*au))
      end if
    end if

  end subroutine find_opening_angle_spherical3

  subroutine mask(clt_p, lon_p, mask_value, i)
    use mod_global_parameters
    double precision          :: clt_p, lon_p
    double precision, dimension(num_cmes)           :: theta_opening_angle, t_half, distance, distance1
    integer, intent(out)      :: mask_value
    integer                   :: i

    theta_opening_angle(i) = -1

    call find_half_time(t_half, i)
    !call find_opening_angle(t_half, theta_opening_angle,i)
    !call find_opening_angle_spherical(t_half, theta_opening_angle,i)
    call find_opening_angle_spherical3(t_half, theta_opening_angle, i)

    if (theta_opening_angle(i) > 0) then

      !distance(i) = (clt_p - clt_cme(i))**2 + (lon_p - lon_cme(i))**2
      distance(i) = acos(sin(clt_p)*sin(clt_cme(i))*cos(lon_p-lon_cme(i)) + cos(clt_p)*cos(clt_cme(i)))
      if (lon_cme(i) - theta_opening_angle(i) < 0) then
        !distance1(i) = (clt_p - clt_cme(i))**2 + (-(2*dpi - lon_p) - lon_cme(i))**2
        distance1(i) = acos(sin(clt_p)*sin(clt_cme(i))*cos(-(2*dpi - lon_p)-lon_cme(i)) + cos(clt_p)*cos(clt_cme(i)))
      else if (lon_cme(i) + theta_opening_angle(i) > 2*dpi) then
        !distance1(i) = (clt_p - clt_cme(i))**2 + (2*dpi + lon_p - lon_cme(i))**2
        distance1(i) = acos(sin(clt_p)*sin(clt_cme(i))*cos((lon_p + 2*dpi)-lon_cme(i)) + cos(clt_p)*cos(clt_cme(i)))
      else
        distance1(i) = 1000
      end if

      if (distance(i) < theta_opening_angle(i)**2) then
        mask_value = 1
      else if (distance1(i) < theta_opening_angle(i)**2) then
        mask_value = 1
      else
        mask_value = 0
      end if
    else
      mask_value = 0
    end if

  end subroutine mask

end module mod_usr
