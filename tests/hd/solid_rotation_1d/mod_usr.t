module mod_usr

  use mod_hd

  implicit none
  double precision :: au2cm  = 1.49597870691d13 ! a.u. to cm         conversion factor
  double precision :: msun2g = 1.988d33         ! solar mass to gram conversion factor
  double precision :: yr2s = 3.1536d7

  ! &usr_list
  double precision :: rho0  = 1d0
  double precision :: vphi0 = 1d0
  double precision :: usr_g = -1d0
  double precision :: dust2gas_ratio = 1d-2

contains

  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods, only: usr_init_one_grid, usr_set_parameters, usr_special_bc
  
    {^IFONED call set_coordinate_system("polar_1.5D")}
    {^IFTWOD call set_coordinate_system("polar_2D")}
    {^IFTHREED call mpistop("3D case not implemented")}

    usr_init_one_grid   => ic
    usr_set_parameters  => set_params
    usr_gravity         => rgravity

    ! Choose independent normalization units if using dimensionless variables.
    unit_length  = au2cm                     ! 1au            (cm)
    unit_density = msun2g / au2cm**2         ! 1M_sun / au^2  (g/cm^2)
    unit_time    = yr2s

    call hd_activate()
    call read_usr_parameters(par_files)
  end subroutine usr_init

  subroutine read_usr_parameters(files)
    !> Read parameters from .par files
    use mod_global_parameters, only: unitpar

    character(len=*), intent(in) :: files(:)
    integer ifile

    namelist /usr_list/ usr_g, dust2gas_ratio, rho0, vphi0

    do ifile = 1, size(files)
       open(unitpar, file=trim(files(ifile)), status="old")
       read(unitpar, usr_list, end=111)
111    rewind(unitpar)
    end do
  end subroutine read_usr_parameters


  subroutine set_params()
    ! Overwrite some default parameters
    use mod_dust, only: dust_n_species, dust_density
    use mod_global_parameters
    ! .. local ..
    double precision :: norm_density
    integer i

    if (hd_energy) call mpistop("usr case not implemented (hd_energy)")

  end subroutine set_params


  subroutine ic(ixI^L, ixO^L, w, x)
    ! Set up initial conditions
    use mod_global_parameters
    use mod_dust, only: dust_n_species, dust_rho, dust_mom, dust_size
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: gradp_r(ixI^S)
    double precision :: dust2gas_frac0, partial_dust2gas_fracs(dust_n_species)
    integer idust

    w(ixO^S, mom(r_)) = 0.0d0

    w(ixO^S, rho_) = rho0
    w(ixO^S, mom(phi_)) = rho0*vphi0

    ! dust ----------------------------------
    ! we compute partial dust to gas fractions assuming the total
    ! fraction = 1/gas2dust_ratio and a size distribution of
    ! n(s) \proto s**-3.5 (standard collisional equilibriuum
    ! assumption from debris disks)

    if (hd_dust) then
       do idust = 1, dust_n_species
          w(ixO^S, dust_rho(idust)) = w(ixO^S, rho_) * dust2gas_ratio
          w(ixO^S, dust_mom(r_, idust)) = 0.0d0
          w(ixO^S, dust_mom(phi_, idust)) = w(ixO^S, mom(phi_)) * dust2gas_ratio
       end do
    end if
  end subroutine ic


  subroutine rgravity(ixI^L, ixO^L, wCT, x, gravity_field)
    ! Set up a time-indep, radial gravity field g \propto r^-1
    ! so as to support solid disk rotation
    use mod_global_parameters
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: x(ixI^S, 1:ndim), wCT(ixI^S, 1:nw)
    double precision, intent(out) :: gravity_field(ixI^S, 1:ndim)

    ! usr_g has dimensionality_of(Gm/r^2) = LT-2
    gravity_field(ixO^S, r_) = usr_g / x(ixO^S, r_)
  end subroutine rgravity

end module mod_usr
