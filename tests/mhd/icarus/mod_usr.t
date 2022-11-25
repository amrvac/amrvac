module mod_usr
  use mod_mhd
  implicit none
  !Shared over subroutines
  real(kind=8), allocatable :: coord_grid_init(:,:,:),variables_init(:,:,:)
  real(kind=8), allocatable :: satellite_positions(:,:)
  real(kind=8), allocatable :: mars_positions(:,:)
  real(kind=8), allocatable :: venus_positions(:,:)
  real(kind=8), allocatable :: sta_positions(:,:)
  real(kind=8), allocatable :: stb_positions(:,:)
  real(kind=8), allocatable :: mercury_positions(:,:)

  !SI units
  ! define my cme parameters here
  double precision    :: k_b    = 1.38064853d-23
  double precision    :: mu_0   = 4.d0*dpi*1.0d-7
  double precision    :: m_p    = 1.67262158d-27
  double precision    :: last_earth=0, last_mercury=0, last_mars=0, last_sta=0, last_stb=0
  double precision    :: last_venus=0
  double precision    :: delta_phi
  !integer             :: starting_index = 0, magnetogram_index = 0, cme_index = 0
  integer, dimension(:), allocatable  :: starting_index(:), cme_index(:)
  integer              :: magnetogram_index = 0, last_index = 0

  ! CME parameters
  integer, dimension(:), allocatable             :: cme_year, cme_month, cme_day, cme_hour, cme_minute, cme_second
  double precision, dimension(:), allocatable    :: relaxation, cme_insertion, vr_cme, w_half, clt_cme, lon_cme, rho_cme, temperature_cme
  double precision, dimension(:), allocatable    :: timestamp, time_difference_cme_magn, longitudes_fix
  double precision, dimension(:), allocatable    :: lon_updated, lon_original
  integer             :: magnetogram_timestamp(6)

  integer             :: num_cmes, cme_exists

contains

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods

    usr_set_parameters  => initglobaldata_usr
    usr_init_one_grid   => initonegrid_usr
    usr_aux_output      => specialvar_output
    usr_add_aux_names   => specialvarnames_output
    usr_special_bc      => specialbound_usr
    usr_refine_grid     => specialrefine_grid
    usr_source          => specialsource

    call set_coordinate_system('spherical_3D')
    call mhd_activate()

  end subroutine usr_init
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------





!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
  subroutine initglobaldata_usr
    use mod_global_parameters
    logical, save       :: firstglobalusr=.true.
    integer             :: AllocateStatus, DeAllocateStatus
    !Parameters related to the unit conversions
    double precision    :: Lunit, Tunit, Rhounit, Vunit, Bunit, Eunit, Punit
    !Parameters related to the read-in of the boundary file
    character(len=50)   :: boundary_file
    character(len=50)   :: cme_parameter_file
    character(len=50)   :: satelite_trajectory, mars_trajectory, venus_trajectory
    character(len=50)   :: sta_trajectory, stb_trajectory, mercury_trajectory
    integer             :: nr_colat, nr_lon, k, n


    double precision    :: time_difference
    !-----------------------------------------------------------------------------

    if(firstglobalusr) then
        call print_initial_information()
        !Read-in coronal model
        !Coronal model at boundary has 2 coordinates: colat and lon
        !Coronal model has data for 4 parameters: vr, n, T, Br
        cme_parameter_file = "cme_input_parameters.dat"   !"cme_cone_model_old.par"
        boundary_file = "solar_wind_bc_used_in_paper.dat"

        ! 2015 event june corresponding satellite data
        satelite_trajectory = "2015_june_earth_ext.unf"
        mars_trajectory = "2015_june_mars_ext.unf"
        venus_trajectory = "2015_june_venus_ext.unf"
        mercury_trajectory = "2015_june_mercury_ext.unf"
        sta_trajectory = "2015_june_sta_ext.unf"
        stb_trajectory = "2015_june_stb_ext.unf"


        call grid_info_coronal_model(boundary_file, nr_colat, nr_lon)

        !We have 4 extra points for longitude because the boundary is periodic
        ALLOCATE(coord_grid_init(nr_colat, nr_lon+4, ndim-1), STAT = AllocateStatus)
            IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
        ALLOCATE(variables_init(nr_colat, nr_lon+4, 4), STAT = AllocateStatus)
            IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

        call read_boundary_coronal_model(boundary_file, coord_grid_init, variables_init, delta_phi)

        !read in cme parameters
        call read_cme_parameters(cme_parameter_file)

        do n = 1, num_cmes
          cme_index(n) = 0
          starting_index(n) = 0
          time_difference_cme_magn(n) = 0
        end do
        ! read in trajectory coordinates
        call read_earth_trajectory(satelite_trajectory)
        call read_mars_trajectory(mars_trajectory)
        call read_venus_trajectory(venus_trajectory)
        call read_mercury_trajectory(mercury_trajectory)
        call read_sta_trajectory(sta_trajectory)
        call read_stb_trajectory(stb_trajectory)

        ! calculate timestamp for cme insertion
        do k=1, num_cmes
          timestamp(k) = relaxation(k)*24.0 + cme_insertion(k)*24.0 + time_difference_cme_magn(k)
        end do
        call cme_insertion_longitudes_fix()
    end if

    mhd_gamma= 3.0d0/2.0d0
    call set_units(Lunit, Tunit, Rhounit, Vunit, Bunit, Eunit, Punit)

  end subroutine initglobaldata_usr
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------





!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
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
    !-----------------------------------------------------------------------------

    w(ix^S,1:nw) = zero
    r_boundary   = 2.15d1 !in R_sun

    !Loop non-ghost cells of cell-block
    {do ix^DB=ixmin^DB,ixmax^DB\}
      !Information of grid coordinate
      xloc(:)=x(ix^D,:)
      !Find the corresponding grid numbers of clt and lon that are closest to the value of xloc
      !point11 and point22 similar to wiki of linear interpolation in 2D
      call find_indices_coord_grid(xloc, point11_clt, point11_lon, point22_clt, point22_lon)

      ! INTERPOLATION
      w(ix^D, mom(1)) = linear_interpolation(xloc(2), xloc(3), point11_clt, point11_lon, point22_clt, point22_lon, 1)/unit_velocity
      w(ix^D, rho_)   = linear_interpolation(xloc(2),xloc(3),point11_clt, point11_lon, point22_clt, point22_lon, 2)/unit_density*(r_boundary/xloc(1))**2
      w(ix^D, p_)     = linear_interpolation(xloc(2),xloc(3),point11_clt, point11_lon, point22_clt, point22_lon,3)/unit_pressure*(r_boundary/xloc(1))**2
      w(ix^D, mag(1)) = linear_interpolation(xloc(2),xloc(3),point11_clt, point11_lon, point22_clt, point22_lon,4)/unit_magneticfield*(r_boundary/xloc(1))**2


      if(w(ix^D, p_)<0.0d0) then
          print*, "NEGATIVE PRESSURE when setting the initial grid: ", w(ix^D, p_)
      end if

    {end do\}

    !Convert to conserved values
    call mhd_to_conserved(ixG^L,ix^L,w,x)

    if(mhd_n_tracer ==  1) then
       w(ix^S, tracer(1)) = 0.0d0
    end if

  end subroutine initonegrid_usr
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------



!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------


  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
    use mod_global_parameters
    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)

    double precision :: divb(ixI^S), divmom(ixI^S)
    double precision :: v(ixI^S,ndir), divV(ixI^S), momentum(ixI^S, ndir)
    integer :: i
    double precision :: r_out, r_in, ratio, block_num, degree, step_size(ixI^S)


    block_num = 61.0
    r_out = 432.5
    r_in = 21.5
    ratio = r_out/r_in
    degree = (block_num-1)/block_num
    step_size(ixI^S) = x(ixI^S,1)*(ratio**degree-1)

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


!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------



!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------


  subroutine specialvarnames_output(varnames)
      character(len=*) :: varnames

      varnames='divB divV div_mom'
end subroutine specialvarnames_output



!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

subroutine readout_satellite(satellite_name, file_name, ixI^L,ixO^L, before_cme, qt,w,x)
  double precision, intent(inout) :: w(ixI^S,1:nw)
  integer, intent(in)             :: ixI^L, ixO^L
  double precision, intent(in)    :: x(ixI^S,1:ndim), qt
  double precision   :: rho_e, P_e, vr_e, vclt_e, vlon_e, Br_e, Bclt_e, Blon_e
  character(len=250) :: file_name
  character(len=20)  :: satellite_name

  !integer            :: last_index
  double precision   :: before_cme
  logical, save      :: exist
  double precision   :: xloc(1:ndim)
  double precision   :: r_satellite, theta_satellite, phi_satellite
  integer            :: rmin, thmin, phimin
  double precision   :: xd, yd, zd
  integer            :: ir1, ith2, iphi3, i_date
  double precision   :: p_r, p_th, p_phi, p_rme, p_thme, p_phime
  integer            :: i, check
  real(8)            :: yr, mnth, day, hr, mn, sc


  if (satellite_name == "earth") then
  r_satellite = satellite_positions(7,last_index)
  theta_satellite = dpi/2.0 - satellite_positions(8, last_index)
  phi_satellite = delta_phi + &
            satellite_positions(9, last_index)&
           - (qt-(timestamp(1)-before_cme))*(2.0*dpi)/24.0*(1/2.447d1-1/365.24)

  else if (satellite_name == "mercury") then
  r_satellite = mercury_positions(7,last_index)
  theta_satellite = dpi/2.0 - mercury_positions(8, last_index)
  phi_satellite = delta_phi + mercury_positions(9, last_index)&
           - (qt-(timestamp(1)-before_cme))*(2.0*dpi)/24.0*(1/2.447d1-1/87.969)

  else if (satellite_name == "mars") then
  r_satellite = mars_positions(7,last_index)
  theta_satellite = dpi/2.0 - mars_positions(8, last_index)
  phi_satellite = delta_phi + mars_positions(9, last_index)&
           - (qt-(timestamp(1)-before_cme))*(2.0*dpi)/24.0*(1/2.447d1-1/686.98)

  else if (satellite_name == "sta") then
  r_satellite = sta_positions(7,last_index)
  theta_satellite = dpi/2.0 - sta_positions(8, last_index)
  phi_satellite = delta_phi + sta_positions(9, last_index)&
           - (qt-(timestamp(1)-before_cme))*(2.0*dpi)/24.0*(1/2.447d1-1/346.0)

  else if (satellite_name == "stb") then
  r_satellite = stb_positions(7,last_index)
  theta_satellite = dpi/2.0 - stb_positions(8, last_index)
  phi_satellite = delta_phi + stb_positions(9, last_index)&
           - (qt-(timestamp(1)-before_cme))*(2.0*dpi)/24.0*(1/2.447d1-1/388.0)

  else if (satellite_name == "venus") then
  r_satellite = venus_positions(7,last_index)
  theta_satellite = dpi/2.0 - venus_positions(8, last_index)
  phi_satellite = delta_phi + venus_positions(9, last_index)&
           - (qt-(timestamp(1)-before_cme))*(2.0*dpi)/24.0*(1/2.447d1-1/224.7)
  end if

  if (phi_satellite < 0) then
     phi_satellite = 2 * dpi + phi_satellite
  else if (phi_satellite > 2*dpi) then
     phi_satellite = mod(phi_satellite, 2.0*dpi)
  end if



     if((x(ixImin1+1,ixImin2,ixImin3,1)<= r_satellite) .and. &
        (x(ixImax1-2,ixImin2,ixImin3,1) >= r_satellite)) then

     if((x(ixImin1,ixImin2+1,ixImin3,2)<= theta_satellite) .and. &
       (x(ixImin1,ixImax2-2,ixImin3,2) >= theta_satellite)) then

     if((x(ixImin1,ixImin2,ixImin3+1,3)<= phi_satellite) .and. &
        (x(ixImin1,ixImin2,ixImax3-2,3) >= phi_satellite)) then

  if (satellite_name == "earth") then
    last_earth= qt
  else if (satellite_name == "mercury") then
    last_mercury= qt
  else if (satellite_name == "mars") then
    last_mars= qt
  else if (satellite_name == "sta") then
    last_sta= qt
  else if (satellite_name == "stb") then
    last_stb= qt
  else if (satellite_name == "venus") then
    last_venus= qt
  end if



   rmin = 0
   thmin = 0
   phimin = 0

   do ir1=ixImin1+1,ixImax1-2
       if( ( x(ir1,ixImin2,ixImin3,1) <= r_satellite) .and. &
           ( x(ir1+1,ixImin2,ixImin3,1) > r_satellite)) then
         rmin = ir1
       end if
   end do
   do ith2=ixImin2+1,ixImax2-2
       if( ( x(ixImin1,ith2,ixImin3,2) <= theta_satellite) .and. &
           ( x(ixImin1,ith2+1,ixImin3,2) > theta_satellite)) then
         thmin=ith2
       end if
   end do
   do iphi3=ixImin3+1,ixImax3-2
       if( ( x(ixImin1,ixImin2,iphi3,3) <= phi_satellite) .and. &
           ( x(ixImin1,ixImin2,iphi3+1,3) > phi_satellite)) then
         phimin=iphi3
       end if
   end do


xd = (r_satellite - x(rmin, thmin, phimin,1))/(x(rmin+1, thmin+1, phimin+1,1)-x(rmin, thmin, phimin,1))
yd = (theta_satellite - x(rmin, thmin, phimin,2))/(x(rmin+1, thmin+1, phimin+1,2)-x(rmin, thmin, phimin,2))
zd = (phi_satellite - x(rmin, thmin, phimin,3))/(x(rmin+1, thmin+1, phimin+1,3)-x(rmin, thmin, phimin,3))


rho_e = linear_interpolation_3D(xd, yd, zd, &
       w(rmin,thmin,phimin,rho_), &
       w(rmin+1,thmin,phimin,rho_), &
       w(rmin,thmin,phimin+1,rho_), &
       w(rmin+1,thmin,phimin+1,rho_), &
       w(rmin,thmin+1,phimin,rho_), &
       w(rmin+1,thmin+1,phimin,rho_), &
       w(rmin,thmin+1,phimin+1,rho_), &
       w(rmin+1,thmin+1,phimin+1,rho_))

P_e = linear_interpolation_3D(xd, yd, zd, &
       w(rmin,thmin,phimin,p_), &
       w(rmin+1,thmin,phimin,p_), &
       w(rmin,thmin,phimin+1,p_), &
       w(rmin+1,thmin,phimin+1,p_), &
       w(rmin,thmin+1,phimin,p_), &
       w(rmin+1,thmin+1,phimin,p_), &
       w(rmin,thmin+1,phimin+1,p_), &
       w(rmin+1,thmin+1,phimin+1,p_))

vr_e = linear_interpolation_3D(xd, yd, zd, &
       w(rmin,thmin,phimin,mom(1))/w(rmin,thmin,phimin,rho_), &
       w(rmin+1,thmin,phimin,mom(1))/w(rmin+1,thmin,phimin,rho_), &
       w(rmin,thmin,phimin+1,mom(1))/w(rmin,thmin,phimin+1,rho_), &
       w(rmin+1,thmin,phimin+1,mom(1))/w(rmin+1,thmin,phimin+1,rho_)&
       ,w(rmin,thmin+1,phimin,mom(1))/w(rmin,thmin+1,phimin,rho_), &
       w(rmin+1,thmin+1,phimin,mom(1))/w(rmin+1,thmin+1,phimin,rho_)&
       ,w(rmin,thmin+1,phimin+1,mom(1))/w(rmin,thmin+1,phimin+1,rho_)&
       ,w(rmin+1,thmin+1,phimin+1,mom(1))/w(rmin+1,thmin+1,phimin+1,rho_))

vclt_e = linear_interpolation_3D(xd, yd, zd, &
       w(rmin,thmin,phimin,mom(2))/w(rmin,thmin,phimin,rho_), &
       w(rmin+1,thmin,phimin,mom(2))/w(rmin+1,thmin,phimin,rho_), &
       w(rmin,thmin,phimin+1,mom(2))/w(rmin,thmin,phimin+1,rho_), &
       w(rmin+1,thmin,phimin+1,mom(2))/w(rmin+1,thmin,phimin+1,rho_)&
       ,w(rmin,thmin+1,phimin,mom(2))/w(rmin,thmin+1,phimin,rho_), &
       w(rmin+1,thmin+1,phimin,mom(2))/w(rmin+1,thmin+1,phimin,rho_)&
       ,w(rmin,thmin+1,phimin+1,mom(2))/w(rmin,thmin+1,phimin+1,rho_)&
       ,w(rmin+1,thmin+1,phimin+1,mom(2))/w(rmin+1,thmin+1,phimin+1,rho_))

vlon_e = linear_interpolation_3D(xd, yd, zd, &
       w(rmin,thmin,phimin,mom(3))/w(rmin,thmin,phimin,rho_), &
       w(rmin+1,thmin,phimin,mom(3))/w(rmin+1,thmin,phimin,rho_), &
       w(rmin,thmin,phimin+1,mom(3))/w(rmin,thmin,phimin+1,rho_), &
       w(rmin+1,thmin,phimin+1,mom(3))/w(rmin+1,thmin,phimin+1,rho_)&
       ,w(rmin,thmin+1,phimin,mom(3))/w(rmin,thmin+1,phimin,rho_), &
       w(rmin+1,thmin+1,phimin,mom(3))/w(rmin+1,thmin+1,phimin,rho_)&
       ,w(rmin,thmin+1,phimin+1,mom(3))/w(rmin,thmin+1,phimin+1,rho_)&
       ,w(rmin+1,thmin+1,phimin+1,mom(3))/w(rmin+1,thmin+1,phimin+1,rho_))

Br_e = linear_interpolation_3D(xd, yd, zd, &
       w(rmin,thmin,phimin,mag(1)), &
       w(rmin+1,thmin,phimin,mag(1)), &
       w(rmin,thmin,phimin+1,mag(1)), &
       w(rmin+1,thmin,phimin+1,mag(1)), &
       w(rmin,thmin+1,phimin,mag(1)), &
       w(rmin+1,thmin+1,phimin,mag(1)), &
       w(rmin,thmin+1,phimin+1,mag(1)), &
       w(rmin+1,thmin+1,phimin+1,mag(1)))

Bclt_e = linear_interpolation_3D(xd, yd, zd, &
       w(rmin,thmin,phimin,mag(2)), &
       w(rmin+1,thmin,phimin,mag(2)), &
       w(rmin,thmin,phimin+1,mag(2)), &
       w(rmin+1,thmin,phimin+1,mag(2)), &
       w(rmin,thmin+1,phimin,mag(2)), &
       w(rmin+1,thmin+1,phimin,mag(2)), &
       w(rmin,thmin+1,phimin+1,mag(2)), &
       w(rmin+1,thmin+1,phimin+1,mag(2)))

Blon_e = linear_interpolation_3D(xd, yd, zd, &
       w(rmin,thmin,phimin,mag(3)), &
       w(rmin+1,thmin,phimin,mag(3)), &
       w(rmin,thmin,phimin+1,mag(3)), &
       w(rmin+1,thmin,phimin+1,mag(3)), &
       w(rmin,thmin+1,phimin,mag(3)), &
       w(rmin+1,thmin+1,phimin,mag(3)), &
       w(rmin,thmin+1,phimin+1,mag(3)), &
       w(rmin+1,thmin+1,phimin+1,mag(3)))

p_r  = r_satellite/214.93946938      ![AU] from solar radii
p_th = theta_satellite               ![rad]
p_phi= phi_satellite                 ![rad]

rho_e = rho_e * unit_density / 0.5d0 / m_p /1.0d6
!rho_e = rho_e * unit_density / 0.5d0 / m_p
P_e = P_e*unit_pressure ![ks/s/m^2] this is energy
vr_e = vr_e *unit_velocity / 1.0d3 !km/s
vclt_e = vclt_e *unit_velocity / 1.0d3 !km/s
vlon_e = vlon_e *unit_velocity / 1.0d3 !km/s
Br_e =Br_e *unit_magneticfield*1.0d9 ![nT], from T to nT
Bclt_e=Bclt_e *unit_magneticfield*1.0d9 ![nT], from T to nT
Blon_e = Blon_e *unit_magneticfield*1.0d9 ![nT], from T to nT
yr = satellite_positions(1,last_index)
mnth = satellite_positions(2,last_index)
day = satellite_positions(3,last_index)
hr = satellite_positions(4,last_index)
mn = satellite_positions(5,last_index)
sc = satellite_positions(6,last_index)

inquire(file=file_name, exist=exist)
if (exist) then
open(12, file=file_name, status="old", position="append", action="write")
write(12,*) yr, mnth, day, hr, mn, sc, qt, p_r, p_th, p_phi, rho_e, P_e, vr_e, vclt_e, vlon_e, Br_e, Bclt_e, Blon_e
!write(12,*) qt, p_r, p_th, p_phi, rho_e, P_e, vr_e, vclt_e, vlon_e, Br_e, Bclt_e, Blon_e
close(12)
else
open(12, file=file_name, status="new", action="write")
write(12, *) "yy mm dd hh MM ss time r clt lon n P vr vclt vlon Br Bclt Blon"
write(12,*) yr, mnth, day, hr, mn, sc, qt, p_r, p_th, p_phi, rho_e, P_e, vr_e, vclt_e, vlon_e, Br_e, Bclt_e, Blon_e
!write(12,*) qt, p_r, p_th, p_phi, rho_e, P_e, vr_e, vclt_e, vlon_e, Br_e, Bclt_e, Blon_e
close(12)
end if



end if
end if
end if



end subroutine readout_satellite




!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
  subroutine specialsource(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, iw^LIM
    integer, intent(in)            :: ixO^L
    double precision, intent(in)    :: qdt, qtC, qt
    double precision, intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision   :: test(ixI^S, 1:ndim)
    double precision :: omega_normalized, omega = 2.97d-6     !Rotation rate: rad/s
    double precision :: MSun = 1.9891d30   !kg
    double precision :: Gravconst_Msun_normalized , GravConst = 6.67384d-11 !SI units

    logical, save      :: exist
    double precision   :: xloc(1:ndim)
    double precision   :: r_earth, theta_earth, phi_earth
    double precision   :: r_mercury, theta_mercury, phi_mercury
    integer            :: ir1, ith2, iphi3, i_date
    integer            :: ir1_me, ith2_me, iphi3_me, i_date_me
    integer            :: rmin, thmin, phimin
    integer            :: rmin_me, thmin_me, phimin_me
    !double precision   :: x1, x2, y1, y2, z1, z
    double precision   :: xd, yd, zd, yd_me, xd_me, zd_me
    double precision   :: p_r, p_th, p_phi, p_rme, p_thme, p_phime
    double precision   :: rho_e, P_e, vr_e, vclt_e, vlon_e, Br_e, Bclt_e, Blon_e
    double precision   :: rho_me, P_me, vr_me, vclt_me, vlon_me, Br_me, Bclt_me, Blon_me
    double precision   :: rho_earth
    character(len=250) :: earth_file, mars_file, sta_file, stb_file
    character(len=250) :: mercury_file, venus_file
    character(len=50)  :: satellite_name

    double precision   :: delta_time
    integer            :: delta_steps
    !integer            :: last_index
    double precision   :: elapsed_time !from the start of the simulation until magnetogram time
    double precision   :: before_cme ! Time from magnetogram time to cme starting time
    double precision   :: current_time
    integer           :: r_min, r_max, th_min, th_max, phi_min, phi_max
    integer            :: r_min1, r_max1, th_min1, th_max1, phi_min1, phi_max1
    !integer            :: i
    integer            :: ix1, ix2, ix3
    integer            :: check
    ! ---------------------------------------------------------------------------------
    earth_file = "./multi_cme_2015june_master_lc_earth.txt"
    mars_file = "./multi_cme_2015june_master_lc_mars.txt"
    mercury_file = "./multi_cme_2015june_master_lc_mercury.txt"
    sta_file = "./multi_cme_2015june_master_lc_sta.txt"
    stb_file = "./multi_cme_2015june_master_lc_stb.txt"
    venus_file = "./multi_cme_2015june_master_lc_venus.txt"

    delta_time = 0.25d0
    delta_steps = int(timestamp(1)*60) !each step is one minute there fore 1 : 60 h
    elapsed_time = (magnetogram_index-starting_index(1))*1/60.0
    before_cme = (cme_index(1) - magnetogram_index)/60.0

    if (mod(qt,0.25d0)<0.25d0 .and. mod(qt+qdt,0.25d0)<mod(qt,0.25d0)) then
        last_index = starting_index(1) + int(qt*60.0)
        if (qt>last_earth) then
          satellite_name = "earth"
          !call readout_satellite(satellite_name, earth_file, last_index, ixI^L, ixO^L, before_cme, qt,w,x)
           call readout_satellite(satellite_name, earth_file, ixI^L, ixO^L, before_cme, qt,w,x)
        end if
        if (qt>last_mercury) then
            satellite_name = "mercury"
            !call readout_satellite(satellite_name, mercury_file, last_index, ixI^L, ixO^L, before_cme, qt,w,x)
            call readout_satellite(satellite_name, mercury_file, ixI^L, ixO^L, before_cme, qt,w,x)
        end if
        if (qt>last_mars) then
            satellite_name = "mars"
            !call readout_satellite(satellite_name, mars_file, last_index, ixI^L, ixO^L, before_cme, qt,w,x)
            call readout_satellite(satellite_name, mars_file, ixI^L, ixO^L, before_cme, qt,w,x)
        end if
        if (qt>last_sta) then
            satellite_name = "sta"
            !call readout_satellite(satellite_name, sta_file, last_index, ixI^L, ixO^L, before_cme, qt,w,x)
            call readout_satellite(satellite_name, sta_file, ixI^L, ixO^L, before_cme, qt,w,x)
        end if
        if (qt>last_stb) then
            satellite_name = "stb"
            !call readout_satellite(satellite_name, stb_file, last_index, ixI^L, ixO^L, before_cme, qt,w,x)
            call readout_satellite(satellite_name, stb_file, ixI^L, ixO^L, before_cme, qt,w,x)
        end if
        if (qt>last_venus) then
            satellite_name = "venus"
            !call readout_satellite(satellite_name, venus_file, last_index, ixI^L, ixO^L, before_cme, qt,w,x)
            call readout_satellite(satellite_name, venus_file, ixI^L, ixO^L, before_cme, qt,w,x)
        end if

  end if




    omega_normalized = omega * unit_length/unit_velocity
    Gravconst_Msun_normalized = GravConst*MSun / (unit_velocity**2 * unit_length )
    !momentum equation -> source = F = rho.g + rho (omega x r)x omega + 2(v x omega)
    !energy equation   -> source  = v . F
    !Gravity
    w(ixO^S,e_)  = w(ixO^S,e_)  - qdt*wCT(ixO^S,mom(1))*Gravconst_Msun_normalized/(x(ixO^S,1)**2)
    w(ixO^S,mom(1)) = w(ixO^S,mom(1)) - qdt*wCT(ixO^S,rho_)*Gravconst_Msun_normalized/(x(ixO^S,1)**2)
    !Centrifugal
    w(ixO^S,mom(1)) = w(ixO^S,mom(1)) + qdt*( omega_normalized*omega_normalized*x(ixO^S,1)*wCT(ixO^S,rho_)*sin(x(ixO^S,2))*sin(x(ixO^S,2)) )
    w(ixO^S,mom(2)) = w(ixO^S,mom(2)) + qdt*( omega_normalized*omega_normalized*x(ixO^S,1)*wCT(ixO^S,rho_)*cos(x(ixO^S,2))*sin(x(ixO^S,2)) )
    !No phi component for the centrifugal force.
    w(ixO^S,e_)  = w(ixO^S,e_)  + qdt*( omega_normalized*omega_normalized*x(ixO^S,1)*sin(x(ixO^S,2))*( cos(x(ixO^S,2))*wCT(ixO^S,mom(2))&
       +sin(x(ixO^S,2))*wCT(ixO^S,mom(1)) ) )
    !Coriolis
    w(ixO^S,mom(1)) = w(ixO^S,mom(1)) + qdt*2.0d0*(omega_normalized*sin(x(ixO^S,2)))*wCT(ixO^S,mom(3))
    w(ixO^S,mom(2)) = w(ixO^S,mom(2)) + qdt*2.0d0*(omega_normalized*cos(x(ixO^S,2)))*wCT(ixO^S,mom(3))
    w(ixO^S,mom(3)) = w(ixO^S,mom(3)) + qdt*2.0d0*omega_normalized*( -cos(x(ixO^S,2))*wCT(ixO^S,mom(2))-sin(x(ixO^S,2))*wCT(ixO^S,mom(1)) )
    !v.F for coriolis is zero.

  end subroutine specialsource
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
! LINEAR INTERPOLATION FUNCTION IN 3D
  double precision function linear_interpolation_3D(xd, yd, zd, C000, C100, C001, C101, C010, C110, C011, C111)
    use mod_global_parameters
    double precision :: xd, yd, zd, C000, C100, C001, C101, C010, C110, C011, C111
    double precision :: C00, C01, C10, C11, C0, C1

    !-----------------------------------------------------------------------------
    ! variables CXXX are taken as on wikipedia page on linear interpolation in 3D

    C00 = C000*(1.0d0-xd)+C100*xd
    C01 = C001*(1.0d0-xd)+C101*xd
    C10 = C010*(1.0d0-xd)+C110*xd
    C11 = C011*(1.0d0-xd)+C111*xd

    C0 = C00*(1.0d0-yd)+C10*yd
    C1 = C01*(1.0d0-yd)+C11*yd

    linear_interpolation_3D = C0*(1.0d0-zd)+C1*zd

    return
  end function linear_interpolation_3D
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------



!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
  !subroutine specialrefine_grid(igrid,level,ixG^L,ix^L,qt,w,x,refine,coarsen)
  subroutine specialrefine_grid(igrid,level,ixI^L,ixO^L,qt,w,x,refine,coarsen)
    ! refine = -1 enforce to not refine
    ! refine =  0 doesn't enforce anything
    ! refine =  1 enforce refinement
    ! coarsen = -1 enforce to not coarsen
    !  coarsen =  0 doesn't enforce anything
    ! coarsen =  1 enforce coarsen
    !integer, intent(in)             :: igrid, level, ixG^L, ix^L

    !double precision, intent(in)    :: qt, w(ixG^S,1:nw+nwauxio), x(ixG^S,1:ndim)
    integer, intent(inout)          :: refine, coarsen

    ! added by me
    integer, intent(in)             :: igrid, level, ixI^L, ixO^L
    double precision, intent(in)    :: qt, w(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision                :: v(ixI^S,ndir), divV(ixI^S)
    integer                         :: i, block_num
    double precision                :: threshold, step_size(ixI^S), r_out, r_in, ratio
    integer                         :: ix1, ix2, ix3
    double precision                :: xloc(1:ndim)
    double precision                :: r, theta, phi, r_this, theta_this, phi_this
    double precision                :: lon_cir, u_artificial, omega
    double precision                :: before_cme
    double precision                :: phi_satellite
    !-----------------------------------------------------------------------------
    !Currently no refining!
    block_num = 61
    r_out = 432.5
    r_in = 21.5
    ratio = r_out/r_in
    step_size = x(ixI^S,1)*(ratio**((block_num-1)/block_num)-1)


    !Implementation of a specific cir, values are for 2012 wind

    !if (qt > timestamp) then
    !u_artificial = 500000
    !omega = (2.0*dpi)/24.0/3600.0*(1/2.447d1-1/365.24)
    !do ix3 = ixImin3, ixImax3
    !    do ix2 = ixImin2,ixImax2
    !        do ix1 = ixImin1, ixImax1
    !            lon_cir = x(ix1, ix2, ix3,3) + (x(ix1, ix2, ix3, 1) - 21.5)*6.955d8/u_artificial*omega
    !            if (x(ix1, ix2, ix3, 2) > 85*dpi/180.0) then
    !            if (x(ix1, ix2, ix3, 2) < 95*dpi/180.0) then
    !
    !            if (lon_cir > 5.0) then
    !                if (lon_cir < 7.3) then
    !                    refine = 1
    !                    coarsen = -1
    !                else
    !                    refine = -1
    !                    coarsen = 1
    !                end if
    !            else
    !                refine = -1
    !                coarsen = 1
    !            end if
    !            end if
    !            end if
    !       end do
    !    end do
   !end do
   !end if


    ! To Follow Earth location in the domain

    before_cme = (cme_index(1) - magnetogram_index)/60.0
    phi_satellite = 1.1982114855026409 + satellite_positions(9, last_index)&
           - (qt-(timestamp(1)-before_cme))*(2.0*dpi)/24.0*(1/2.447d1-1/365.24)

    ! Threshold for negative nabla V
    threshold = -0.005

    ! Refinement criterion for div(V)

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

    ! Refinement criterion for tracing function

    !if (qt > timestamp) then
    !    if (any(w(ixI^S,tracer(1)) > 0.001)) then
    !        refine = 1
    !        coarsen = -1
    !    else
    !        refine = -1
    !        coarsen = 1
    !    end if
    !end if



  end subroutine specialrefine_grid
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------



!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
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
    integer             :: mask_cme, n
    double precision    :: omega = 2.97d-6
    !-----------------------------------------------------------------------------

    select case(iB)
        case(1)! Lower radial boundary
            nr_colat = SIZE(coord_grid_init(:,1,1))
            nr_lon = SIZE(coord_grid_init(1,:,1))
            w(ixO^S,:) = 0.d0

            !W IS USING CONSERVATIVE VALUES
            !SO TO GET v1 WE NEED TO USE m1/rho!!

            !rghost1 = x(ixOmax1,1,1,r_): radial coordinate at first ghost cell (closest to in-domain)
            !rghost2 = x(ixOmin1,1,1,r_): radial coordinate at second ghost cell
            !Inner boundary has a linear relation: f(r) = a + b*r
            !We have f(r) boundary 21.5 and r = 21.5 and first physical cell f(r)
            !then b. = f(r)_bound - f(r)_phys / (R_bound - R_phys)
            !and a = f(r)_bound - b*21.5R_sun
            !FOR NOW WE ASSUME NON-STRETCHED GRID!


            do ix3 = ixOmin3, ixOmax3
                do ix2 = ixOmin2,ixOmax2
                    !Get data from the first phsyical cell going in the r-direction
                    r       = x(ixOmin1 +2 , ix2, ix3, 1)
                    theta   = x(ixOmin1, ix2, ix3, 2)
                    phi     = x(ixOmin1, ix2, ix3, 3)
                    xloc(:) = x(ixOmin1 +2, ix2, ix3,:)
                    !IMPORTANT NOTE:
                    !AMRVAC has ghostcells that are largers than 2*pi
                    !We made sure that coord_grid_init also has data >2*pi

                    !Find the corresponding data in the initial grid
                    !For now we take the values at 21.5+delta_r/2 as our grid coord_grid_init does not have 21.5R_sun
                    !TODO: add the 21.5 R_sun grid!!

                    call find_indices_coord_grid(xloc, point11_clt, point11_lon, point22_clt, point22_lon)


                    mask_cme = 0

                    do n=1, num_cmes
                    if (mask_cme .eq. 0) then
                    call mask(xloc(2), xloc(3), mask_cme, n)


                    if (mask_cme .eq. 1) then

                        ! Mass density
                        w(ixOmax1, ix2, ix3, rho_) = rho_cme(n)/unit_density
                        w(ixOmin1, ix2, ix3, rho_) = rho_cme(n)/unit_density

                        ! speed
                        w(ixOmax1, ix2, ix3, mom(1)) =  vr_cme(n)/unit_velocity
                        w(ixOmin1, ix2, ix3, mom(1)) =  vr_cme(n)/unit_velocity


                        w(ixOmax1, ix2, ix3, mom(2)) = - w(ixOmin1+2,ix2, ix3,mom(2))/w(ixOmin1+2,ix2, ix3,rho_)
                        w(ixOmin1, ix2, ix3, mom(2)) = - w(ixOmin1+3,ix2, ix3,mom(2))/w(ixOmin1+2,ix2, ix3,rho_)

                        w(ixOmax1, ix2, ix3, mom(3)) = - w(ixOmin1+2,ix2, ix3,mom(3))/w(ixOmin1+2,ix2, ix3,rho_)
                        w(ixOmin1, ix2, ix3, mom(3)) = - w(ixOmin1+3,ix2, ix3,mom(3))/w(ixOmin1+2,ix2, ix3,rho_)

                        ! magnetic field (keep the way it is)
                        w(ixOmax1, ix2, ix3, mag(1)) =  2.0d0* linear_interpolation(theta, phi, point11_clt, point11_lon, point22_clt, point22_lon, 4)/unit_magneticfield&
                                                    - w(ixOmin1+2,ix2, ix3,mag(1))
                        w(ixOmin1, ix2, ix3, mag(1)) =  4.0d0* linear_interpolation(theta, phi, point11_clt, point11_lon, point22_clt, point22_lon, 4)/unit_magneticfield&
                                                    - 3.0d0*w(ixOmin1+2,ix2, ix3,mag(1))

                        ! B_theta B_phi are 0
                        w(ixOmax1, ix2, ix3, mag(2)) = - w(ixOmin1+2,ix2, ix3,mag(2))
                        w(ixOmax1, ix2, ix3, mag(3)) = - w(ixOmin1+2,ix2, ix3,mag(3))
                        w(ixOmin1, ix2, ix3, mag(2)) = - w(ixOmin1+3,ix2, ix3,mag(2))
                        w(ixOmin1, ix2, ix3, mag(3)) = - w(ixOmin1+3,ix2, ix3,mag(3))

                        !pressure: p
                        w(ixOmax1, ix2, ix3, p_) = rho_cme(n) / (0.5 * M_P) * k_b * temperature_cme(n)/unit_pressure
                        w(ixOmin1, ix2, ix3, p_) = rho_cme(n) / (0.5 * M_P) * k_b * temperature_cme(n)/unit_pressure


                        ! Setting tracer function to the value of densicy inside CME
                        w(ixOmax1, ix2, ix3,  tracer(1)) = rho_cme(n)/unit_density
                        w(ixOmin1, ix2, ix3,  tracer(1)) = rho_cme(n)/unit_density

                    else

                       !Mass density
                       w(ixOmax1, ix2, ix3, rho_) =  2.0d0* linear_interpolation(theta, phi, point11_clt, point11_lon, point22_clt, point22_lon, 2)/unit_density&
                                                    - w(ixOmin1+2,ix2, ix3,rho_)
                       w(ixOmin1, ix2, ix3, rho_) =  4.0d0* linear_interpolation(theta, phi, point11_clt, point11_lon, point22_clt, point22_lon, 2)/unit_density &
                                                    - 3.0d0*w(ixOmin1+2,ix2, ix3,rho_)

                       !Change momentum to speed
                       w(ixOmax1, ix2, ix3, mom(1)) =  2.0d0* linear_interpolation(theta, phi, point11_clt, point11_lon, point22_clt, point22_lon, 1)/unit_velocity&
                                                    - w(ixOmin1+2,ix2, ix3,mom(1))/w(ixOmin1+2,ix2, ix3,rho_)
                       w(ixOmin1, ix2, ix3, mom(1)) =  4.0d0* linear_interpolation(theta, phi, point11_clt, point11_lon, point22_clt, point22_lon, 1)/unit_velocity &
                                                    - 3.0d0*w(ixOmin1+2,ix2, ix3,mom(1))/w(ixOmin1+2,ix2, ix3,rho_)


                       !v2 and v3 are zero at boundary- so we can just copy the value and change the sign
                       w(ixOmax1, ix2, ix3, mom(2)) = - w(ixOmin1+2,ix2, ix3,mom(2))/w(ixOmin1+2,ix2, ix3,rho_)
                       w(ixOmin1, ix2, ix3, mom(2)) = - w(ixOmin1+3,ix2, ix3,mom(2))/w(ixOmin1+2,ix2, ix3,rho_)

                       w(ixOmax1, ix2, ix3, mom(3)) = - w(ixOmin1+2,ix2, ix3,mom(3))/w(ixOmin1+2,ix2, ix3,rho_)
                       w(ixOmin1, ix2, ix3, mom(3)) = - w(ixOmin1+3,ix2, ix3,mom(3))/w(ixOmin1+2,ix2, ix3,rho_)

                       !magnetic field: b1
                       w(ixOmax1, ix2, ix3, mag(1)) =  2.0d0* linear_interpolation(theta, phi, point11_clt, point11_lon, point22_clt, point22_lon, 4)/unit_magneticfield&
                                                    - w(ixOmin1+2,ix2, ix3,mag(1))
                       w(ixOmin1, ix2, ix3, mag(1)) =  4.0d0* linear_interpolation(theta, phi, point11_clt, point11_lon, point22_clt, point22_lon, 4)/unit_magneticfield&
                                                    - 3.0d0*w(ixOmin1+2,ix2, ix3,mag(1))

                       !b2 and b3 are zero at boundary- so we can just copy the value and change the sign
                       w(ixOmax1, ix2, ix3, mag(2)) = - w(ixOmin1+2,ix2, ix3,mag(2))
                       w(ixOmax1, ix2, ix3, mag(3)) = - w(ixOmin1+2,ix2, ix3,mag(3))
                       w(ixOmin1, ix2, ix3, mag(2)) = - w(ixOmin1+3,ix2, ix3,mag(2))
                       w(ixOmin1, ix2, ix3, mag(3)) = - w(ixOmin1+3,ix2, ix3,mag(3))

                       !pressure: p
                       w(ixOmax1, ix2, ix3, p_) = linear_interpolation(theta, phi, point11_clt, point11_lon, point22_clt, point22_lon, 3)/unit_pressure
                       w(ixOmin1, ix2, ix3, p_) = linear_interpolation(theta, phi, point11_clt, point11_lon, point22_clt, point22_lon, 3)/unit_pressure


                       w(ixOmax1, ix2, ix3,  tracer(1)) = - w(ixOmin1+2, ix2, ix3,  tracer(1))
                       w(ixOmin1, ix2, ix3,  tracer(1)) = - w(ixOmin1+3, ix2, ix3,  tracer(1))

                    end if
                    end if

                end do
                end do
            end do


            !convert back to conserved values
            call mhd_to_conserved(ixI^L,ixO^L,w,x)
    end select
  end subroutine specialbound_usr
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------








!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
! LINEAR INTERPOLATION FUNCTION
  double precision function linear_interpolation(x, y, counter_x1, counter_y1, counter_x2, counter_y2, variable)
    use mod_global_parameters

    integer          :: variable, counter_x1, counter_x2, counter_y1, counter_y2
    double precision :: x1, x2, y1, y2, x, y, interpolation_x1, interpolation_x2

    !-----------------------------------------------------------------------------

    x1 = coord_grid_init(counter_x1, counter_y1, 1)
    x2 = coord_grid_init(counter_x2, counter_y2, 1)
    y1 = coord_grid_init(counter_x1, counter_y1, 2)
    y2 = coord_grid_init(counter_x2, counter_y2, 2)
    interpolation_x1 = (x2-x)/(x2-x1)*variables_init(counter_x1, counter_y1, variable)
    interpolation_x1 = interpolation_x1 + (x-x1)/(x2-x1)*variables_init(counter_x2, counter_y1,variable)
    interpolation_x2 = (x2-x)/(x2-x1)*variables_init(counter_x1, counter_y2, variable)
    interpolation_x2 = interpolation_x2 + (x-x1)/(x2-x1)*variables_init(counter_x2, counter_y2, variable)

    linear_interpolation = (y2-y)/(y2-y1)*interpolation_x1 + (y-y1)/(y2-y1)*interpolation_x2

    return
  end function linear_interpolation
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------



!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!Get grid info on coronal model
  subroutine grid_info_coronal_model(filename, nr_colat, nr_lon)
    use mod_global_parameters
    character(len=50), intent(in)   :: filename
    integer, intent(out)            :: nr_colat, nr_lon
    integer                         :: iUnit=20, iError, i
    !-----------------------------------------------------------------------------

    open(iUnit, file=filename, status="old", action="read", iostat=iError)
    if(iError /= 0) call mpistop('Importdata could not open real4 file = '//trim(filename))
    !Note: second line is the time, in case this is needed
    do i = 1, 5
        read(iUnit,*)
    end do
    read(iUnit,*) nr_colat
    do i = 1, nr_colat+2
        read(iUnit,*)
    end do
    read(iUnit,*) nr_lon
    close(iUnit)

  end subroutine grid_info_coronal_model
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------



!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
! read in boundary file
  subroutine read_boundary_coronal_model(filename, coord_grid, variables, delta_phi)
    use mod_global_parameters
    character(len=50), intent(in)   :: filename
    double precision, intent(out)   :: coord_grid(:,:,:), variables(:,:,:)
    double precision, intent(out)   :: delta_phi
    integer                         :: AllocateStatus, DeAllocateStatus
    integer                         :: iUnit=20, iError
    integer                         :: nVar =4, nr_colat, nr_lon, counter, i,j,k
    integer                         :: total_boundary_points
    double precision, dimension(:), allocatable     :: colat_points, lon_points
    double precision, dimension(:,:), allocatable   :: variables_boundary_data

    character(len=20)               :: magnetogram_time
    !-----------------------------------------------------------------------------
    open(iUnit, file=filename, status="old", action="read", iostat=iError)
    if(iError /= 0) call mpistop('Importdata could not open real4 file = '//trim(filename))
    !Note: second line is the time, in case this would ever be needed
    do i = 1, 5
        if (i == 2) then
            read(iUnit, *) magnetogram_time
        else
            read(iUnit,*)
        end if
    end do

    read (magnetogram_time(1:4),*) magnetogram_timestamp(1)
    read (magnetogram_time(6:7),*) magnetogram_timestamp(2)
    read (magnetogram_time(9:10),*) magnetogram_timestamp(3)
    read (magnetogram_time(12:13),*) magnetogram_timestamp(4)
    read (magnetogram_time(15:16),*) magnetogram_timestamp(5)
    read (magnetogram_time(18:19),*) magnetogram_timestamp(6)

    read(iUnit,*) nr_colat
    read(iUnit,*)
    ALLOCATE(colat_points(nr_colat), STAT = AllocateStatus)
            IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
    do i = 1, nr_colat
        read(iUnit,*) colat_points(i)
    end do
    read(iUnit,*)
    read(iUnit,*) nr_lon
    read(iUnit,*)
    ALLOCATE(lon_points(nr_lon), STAT = AllocateStatus)
        IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
    do i = 1, nr_lon
        read(iUnit,*) lon_points(i)
    end do
    read(iUnit,*)
    total_boundary_points = nr_colat*nr_lon
    ALLOCATE(variables_boundary_data(total_boundary_points, nVar), STAT = AllocateStatus)
        IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
    do j = 1,nVar
        do i = 1, total_boundary_points
            read(iUnit,*) variables_boundary_data(i,j)
        end do
        if (j<4) then
            read(iUnit,*)
        end if
    end do
    close(iUnit)

    !Set longitude points to be between 0 and 2pi by shifting smallest value to zero
    if (lon_points(1)<=0) then
        delta_phi =  abs(lon_points(1))
    else if (lon_points(1)>0) then
        delta_phi = - lon_points(1)
    end if



    !Note: Boundary points loop over colat (inner loop) then over longitude (outer loop)
    !k+2 for periodic reasons (longitude)
    !4 parameters: vr, n, T, Br
    counter = 1
    do k = 1, nr_lon
        do j = 1, nr_colat
            coord_grid(j,k+2,1) = colat_points(j)
            coord_grid(j,k+2,2) = lon_points(k) + delta_phi
            !vr stays vr
            variables(j,k+2,1)  = variables_boundary_data(counter,1)
            !number density -> density rho = n_bound * 0.5 * m_p * (r_0/r)**2
            variables(j,k+2,2)  = variables_boundary_data(counter,2) * 0.5 * m_p
            !temperature --> change to pressure = n k_b T
            variables(j,k+2,3) = variables_boundary_data(counter,3) * k_b * variables_boundary_data(counter, 2)
            !magnetic field: Br = Br_bound * sqrt(r_0/r)
            variables(j,k+2,4)  = variables_boundary_data(counter,4)

            if(variables(j,k+1,3)<0) call mpistop('NEGATIVE PRESSURE in boundary file.')

            counter = counter + 1
        end do
    end do

    !k+2 for periodic reasons longitude (copy last two to front and first two to back)
    do j = 1, nr_colat
        coord_grid(j,1,1) = colat_points(j)
        coord_grid(j,1,2) = lon_points(nr_lon-1) + delta_phi - 2*dpi
        variables(j,1,:)  = variables(j,nr_lon-1,:)

        coord_grid(j,2,1) = colat_points(j)
        coord_grid(j,2,2) = lon_points(nr_lon-2) + delta_phi - 2*dpi
        variables(j,2,:)  = variables(j,nr_lon-2,:)

        coord_grid(j,nr_lon+3,1) = colat_points(j)
        coord_grid(j,nr_lon+3,2) = lon_points(1) +delta_phi + 2*dpi
        variables(j,nr_lon+3,:)  = variables(j,3,:)

        coord_grid(j,nr_lon+4,1) = colat_points(j)
        coord_grid(j,nr_lon+4,2) = lon_points(2) + delta_phi + 2*dpi
        variables(j,nr_lon+4,:)  = variables(j,4,:)
    end do

    DEALLOCATE (variables_boundary_data, STAT = DeAllocateStatus)
    DEALLOCATE (lon_points, STAT = DeAllocateStatus)
    DEALLOCATE (colat_points, STAT = DeAllocateStatus)



  end subroutine read_boundary_coronal_model
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

  subroutine set_units(Lunit, Tunit, Rhounit, Vunit, Bunit, Eunit, Punit)
    use mod_global_parameters
    double precision, intent(out)    :: Lunit, Tunit, Rhounit, Vunit, Bunit, Eunit, Punit
    !-----------------------------------------------------------------------------
    !Unit Length: [m] = 1 solar radii
    Lunit = 6.955d8
    !Unit Time : [s] = 1 hour
    Tunit = 6.0d1*6.0d1
    !Unit Mass density: [kg/m^3] = 1.6726d-13 - scaled so that rho is +- 1 in dimensionless form
    Rhounit = 1.6726d-19
    !Velocity unit
    Vunit = Lunit/Tunit
    !Magnetic filed unit
    Bunit = dsqrt( mu_0 * Vunit**2 * Rhounit)
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
        print *, "CME Timestamp[Y/M/D H/M/S]: ", cme_year, cme_month, cme_day, cme_hour, cme_minute, cme_second
        print *,  "Vr [m/s] = ", vr_cme
        print *, "Half width [deg]= ", w_half*180.0/dpi
        print *, "Co-latitude [deg] = ", clt_cme*180.0/dpi
        print *, "Longitude [deg]= ",  lon_cme*180.0/dpi
        print *, "Density [kg/m^3] = ", rho_cme
        print *, "Temperature [K] = ", temperature_cme
        print *, '=================================================================='
        print *, '=================================================================='//NEW_LINE('A')
    end if

    !Set conversion units
    unit_length=Lunit
    unit_density = Rhounit
    unit_time    = Tunit
    unit_velocity  = Vunit
    unit_magneticfield  = Bunit
    unit_pressure = Punit
    unit_numberdensity = unit_density / 0.5 / m_p
  end subroutine set_units

  subroutine print_initial_information()
    use mod_global_parameters
    !-----------------------------------------------------------------------------
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
    !-----------------------------------------------------------------------------
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

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

  subroutine read_earth_trajectory(trajectory_file)

    use mod_global_parameters
    character(len=50), intent(in)   :: trajectory_file
    integer                         :: iUnit=15, iError, i, j
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

    ALLOCATE(satellite_positions(nr_coordinates, nr_positions), STAT = AllocateStatus)


    satellite_positions(1,:) = year
    satellite_positions(2,:) = month
    satellite_positions(3,:) = day
    satellite_positions(4,:) = hour
    satellite_positions(5,:) = minute
    satellite_positions(6,:) = second
    satellite_positions(7,:) = radii
    satellite_positions(8,:) = latitudes
    satellite_positions(9,:) = longitudes

    delta_time = 0.25d0
    ! delta_steps = int(timestamp*60)

    if (magnetogram_index .eq. 0) then
      do j_date = 1, size(year)
        if ((year(j_date) == magnetogram_timestamp(1)) .and.  (month(j_date) == magnetogram_timestamp(2))) then
          if ((day(j_date) == magnetogram_timestamp(3)) .and. (hour(j_date) == magnetogram_timestamp(4))) then
            if (minute(j_date) == magnetogram_timestamp(5)) then
              magnetogram_index = j_date
              exit
            end if
          end if
        end if
      end do
    end if


    if (starting_index(1) .eq. 0) then
      do n = 1, num_cmes
      do i_date = 1, size(year)
        if ((year(i_date) == cme_year(n)) .and.  (month(i_date) == cme_month(n))) then
          if ((day(i_date) == cme_day(n)) .and. (hour(i_date) == cme_hour(n))) then
            if (minute(i_date) == cme_minute(n)) then
              cme_index(n) = i_date
              time_difference_cme_magn(n) = (cme_index(n) - magnetogram_index)/60.0 !hours
              delta_steps = int((relaxation(1)+cme_insertion(1)+time_difference_cme_magn(n))*60)
              starting_index(n) = i_date - delta_steps
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



  end subroutine read_earth_trajectory
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

  subroutine read_mars_trajectory(trajectory_file)

    use mod_global_parameters
    character(len=50), intent(in)   :: trajectory_file
    integer                         :: iUnit=16, iError, i, j
    integer                         :: arr_size(2), nr_positions, nr_coordinates
    double precision, allocatable   :: radii(:), latitudes(:), longitudes(:)
    integer, allocatable            :: year(:), month(:), day(:)
    integer, allocatable            :: hour(:), minute(:), second(:)
    integer                         :: AllocateStatus, DeAllocateStatus
    double precision                :: delta_time
    integer                         :: delta_steps, i_date, j_date



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

    ALLOCATE(mars_positions(nr_coordinates, nr_positions), STAT = AllocateStatus)


    mars_positions(1,:) = year
    mars_positions(2,:) = month
    mars_positions(3,:) = day
    mars_positions(4,:) = hour
    mars_positions(5,:) = minute
    mars_positions(6,:) = second
    mars_positions(7,:) = radii
    mars_positions(8,:) = latitudes
    mars_positions(9,:) = longitudes


    DEALLOCATE(year, STAT = DEAllocateStatus)
    DEALLOCATE(month, STAT = DEAllocateStatus)
    DEALLOCATE(day, STAT = DEAllocateStatus)
    DEALLOCATE(hour, STAT = DEAllocateStatus)
    DEALLOCATE(minute, STAT = DEAllocateStatus)
    DEALLOCATE(second, STAT = DEAllocateStatus)
    DEALLOCATE(radii, STAT = DEAllocateStatus)
    DEALLOCATE(latitudes, STAT = DEAllocateStatus)
    DEALLOCATE(longitudes, STAT = DEAllocateStatus)



  end subroutine read_mars_trajectory
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

  subroutine read_venus_trajectory(trajectory_file)

    use mod_global_parameters
    character(len=50), intent(in)   :: trajectory_file
    integer                         :: iUnit=17, iError, i, j
    integer                         :: arr_size(2), nr_positions, nr_coordinates
    double precision, allocatable   :: radii(:), latitudes(:), longitudes(:)
    integer, allocatable            :: year(:), month(:), day(:)
    integer, allocatable            :: hour(:), minute(:), second(:)
    integer                         :: AllocateStatus, DeAllocateStatus
    double precision                :: delta_time
    integer                         :: delta_steps, i_date, j_date



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

    ALLOCATE(venus_positions(nr_coordinates, nr_positions), STAT = AllocateStatus)


    venus_positions(1,:) = year
    venus_positions(2,:) = month
    venus_positions(3,:) = day
    venus_positions(4,:) = hour
    venus_positions(5,:) = minute
    venus_positions(6,:) = second
    venus_positions(7,:) = radii
    venus_positions(8,:) = latitudes
    venus_positions(9,:) = longitudes


    DEALLOCATE(year, STAT = DEAllocateStatus)
    DEALLOCATE(month, STAT = DEAllocateStatus)
    DEALLOCATE(day, STAT = DEAllocateStatus)
    DEALLOCATE(hour, STAT = DEAllocateStatus)
    DEALLOCATE(minute, STAT = DEAllocateStatus)
    DEALLOCATE(second, STAT = DEAllocateStatus)
    DEALLOCATE(radii, STAT = DEAllocateStatus)
    DEALLOCATE(latitudes, STAT = DEAllocateStatus)
    DEALLOCATE(longitudes, STAT = DEAllocateStatus)



  end subroutine read_venus_trajectory
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

  subroutine read_mercury_trajectory(trajectory_file)

    use mod_global_parameters
    character(len=50), intent(in)   :: trajectory_file
    integer                         :: iUnit=18, iError, i, j
    integer                         :: arr_size(2), nr_positions, nr_coordinates
    double precision, allocatable   :: radii(:), latitudes(:), longitudes(:)
    integer, allocatable            :: year(:), month(:), day(:)
    integer, allocatable            :: hour(:), minute(:), second(:)
    integer                         :: AllocateStatus, DeAllocateStatus
    double precision                :: delta_time
    integer                         :: delta_steps, i_date, j_date



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

    ALLOCATE(mercury_positions(nr_coordinates, nr_positions), STAT = AllocateStatus)


    mercury_positions(1,:) = year
    mercury_positions(2,:) = month
    mercury_positions(3,:) = day
    mercury_positions(4,:) = hour
    mercury_positions(5,:) = minute
    mercury_positions(6,:) = second
    mercury_positions(7,:) = radii
    mercury_positions(8,:) = latitudes
    mercury_positions(9,:) = longitudes


    DEALLOCATE(year, STAT = DEAllocateStatus)
    DEALLOCATE(month, STAT = DEAllocateStatus)
    DEALLOCATE(day, STAT = DEAllocateStatus)
    DEALLOCATE(hour, STAT = DEAllocateStatus)
    DEALLOCATE(minute, STAT = DEAllocateStatus)
    DEALLOCATE(second, STAT = DEAllocateStatus)
    DEALLOCATE(radii, STAT = DEAllocateStatus)
    DEALLOCATE(latitudes, STAT = DEAllocateStatus)
    DEALLOCATE(longitudes, STAT = DEAllocateStatus)



  end subroutine read_mercury_trajectory
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

  subroutine read_sta_trajectory(trajectory_file)

    use mod_global_parameters
    character(len=50), intent(in)   :: trajectory_file
    integer                         :: iUnit=19, iError, i, j
    integer                         :: arr_size(2), nr_positions, nr_coordinates
    double precision, allocatable   :: radii(:), latitudes(:), longitudes(:)
    integer, allocatable            :: year(:), month(:), day(:)
    integer, allocatable            :: hour(:), minute(:), second(:)
    integer                         :: AllocateStatus, DeAllocateStatus
    double precision                :: delta_time
    integer                         :: delta_steps, i_date, j_date



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

    ALLOCATE(sta_positions(nr_coordinates, nr_positions), STAT = AllocateStatus)


    sta_positions(1,:) = year
    sta_positions(2,:) = month
    sta_positions(3,:) = day
    sta_positions(4,:) = hour
    sta_positions(5,:) = minute
    sta_positions(6,:) = second
    sta_positions(7,:) = radii
    sta_positions(8,:) = latitudes
    sta_positions(9,:) = longitudes


    DEALLOCATE(year, STAT = DEAllocateStatus)
    DEALLOCATE(month, STAT = DEAllocateStatus)
    DEALLOCATE(day, STAT = DEAllocateStatus)
    DEALLOCATE(hour, STAT = DEAllocateStatus)
    DEALLOCATE(minute, STAT = DEAllocateStatus)
    DEALLOCATE(second, STAT = DEAllocateStatus)
    DEALLOCATE(radii, STAT = DEAllocateStatus)
    DEALLOCATE(latitudes, STAT = DEAllocateStatus)
    DEALLOCATE(longitudes, STAT = DEAllocateStatus)



  end subroutine read_sta_trajectory
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

  subroutine read_stb_trajectory(trajectory_file)

    use mod_global_parameters
    character(len=50), intent(in)   :: trajectory_file
    integer                         :: iUnit=14, iError, i, j
    integer                         :: arr_size(2), nr_positions, nr_coordinates
    double precision, allocatable   :: radii(:), latitudes(:), longitudes(:)
    integer, allocatable            :: year(:), month(:), day(:)
    integer, allocatable            :: hour(:), minute(:), second(:)
    integer                         :: AllocateStatus, DeAllocateStatus
    double precision                :: delta_time
    integer                         :: delta_steps, i_date, j_date



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

    ALLOCATE(stb_positions(nr_coordinates, nr_positions), STAT = AllocateStatus)


    stb_positions(1,:) = year
    stb_positions(2,:) = month
    stb_positions(3,:) = day
    stb_positions(4,:) = hour
    stb_positions(5,:) = minute
    stb_positions(6,:) = second
    stb_positions(7,:) = radii
    stb_positions(8,:) = latitudes
    stb_positions(9,:) = longitudes


    DEALLOCATE(year, STAT = DEAllocateStatus)
    DEALLOCATE(month, STAT = DEAllocateStatus)
    DEALLOCATE(day, STAT = DEAllocateStatus)
    DEALLOCATE(hour, STAT = DEAllocateStatus)
    DEALLOCATE(minute, STAT = DEAllocateStatus)
    DEALLOCATE(second, STAT = DEAllocateStatus)
    DEALLOCATE(radii, STAT = DEAllocateStatus)
    DEALLOCATE(latitudes, STAT = DEAllocateStatus)
    DEALLOCATE(longitudes, STAT = DEAllocateStatus)



  end subroutine read_stb_trajectory
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------



  subroutine read_cme_parameters(filename)
    use mod_global_parameters
    character(len=50), intent(in)   :: filename
    integer                         :: iUnit=30, iError, i, DEAllocateStatus, j
    character(len=50)               :: cme_parameter_file


    open(iUnit, file = filename, status = 'old', action="read", iostat=iError)
    if(iError /= 0) call mpistop('Importdata could not open real4 file = '//trim(filename))
    read(iUnit,*) num_cmes
    ALLOCATE(cme_year(num_cmes))
    ALLOCATE(cme_month(num_cmes))
    ALLOCATE(cme_day(num_cmes))
    ALLOCATE(cme_hour(num_cmes))
    ALLOCATE(cme_minute(num_cmes))
    ALLOCATE(cme_second(num_cmes))

    ALLOCATE(timestamp(num_cmes))
    ALLOCATE(relaxation(num_cmes))
    ALLOCATE(cme_insertion(num_cmes))
    ALLOCATE(vr_cme(num_cmes))
    ALLOCATE(w_half(num_cmes))
    ALLOCATE(clt_cme(num_cmes))
    ALLOCATE(lon_cme(num_cmes))
    ALLOCATE(rho_cme(num_cmes))
    ALLOCATE(temperature_cme(num_cmes))
    ALLOCATE(cme_index(num_cmes))
    ALLOCATE(starting_index(num_cmes))
    ALLOCATE(time_difference_cme_magn(num_cmes))
    ALLOCATE(longitudes_fix(num_cmes))

    do i=1, num_cmes
      read(iUnit,*) cme_year(i), cme_month(i), cme_day(i), cme_hour(i), cme_minute(i), cme_second(i)
      read(iUnit,*) relaxation(i), cme_insertion(i), vr_cme(i), w_half(i), clt_cme(i), lon_cme(i), rho_cme(i), temperature_cme(i)
      w_half(i) = w_half(i) * dpi/180.0
      clt_cme(i) = (-clt_cme(i) + 90.0) * dpi/180.0
      lon_cme(i) = delta_phi + lon_cme(i) * dpi/180.0

    end do



    close(iUnit)

  end subroutine read_cme_parameters


  subroutine cme_insertion_longitudes_fix()
     integer       :: i

     do i=1, num_cmes
        longitudes_fix(i) = (timestamp(i)-relaxation(i)*24.0 - cme_insertion(i)*24.0)*(2.0*dpi)/24.0*(1/2.447d1-1/365.24)
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
        r_new(i) = (global_time - timestamp(i)) * vr_cme(i)*3600.0 + 0.1*au - r_half(i)

        theta_opening_angle(i)=acos((r_new(i)**2+(0.1*au)**2-r_half(i)**2)/(2.0*r_new(i)*0.1*au))
      end if
    end if

  end subroutine find_opening_angle_spherical3





  subroutine mask(clt_p, lon_p, mask_value, i)
    use mod_global_parameters
    double precision          :: clt_p, lon_p, omega = 2.97d-6
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
