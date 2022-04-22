!> Module for including dust species, which interact with the gas through a drag
!> force
module mod_dust
  use mod_global_parameters, only: std_len
  use mod_physics

  implicit none
  private

  !> The number of dust species
  integer, public, protected      :: dust_n_species = 0

  integer, protected              :: gas_rho_ = -1
  integer, allocatable, protected :: gas_mom(:)
  integer, protected              :: gas_e_   = -1

  !> Indices of the dust densities
  integer, allocatable, public, protected :: dust_rho(:)

  !> Indices of the dust momentum densities
  integer, allocatable, public, protected :: dust_mom(:, :)

  !> Size of each dust species, dimensionless expression
  double precision, allocatable, public :: dust_size(:)

  !> Internal density of each dust species, dimensionless expression
  double precision, allocatable, public :: dust_density(:)

  !> Reduction of stopping time timestep limit
  double precision :: dust_dtpar = 0.5d0

  !> Factor used in squared thermal velocity
  double precision :: gas_vtherm_factor = 3.0d0

  !> Dust temperature in K (if dust_temperature_type is constant)
  double precision :: dust_temperature = -1.0d0

  !> Dust drag coefficient for linear drag (for testing dust_method=linear)
  double precision :: dust_K_lineardrag = -1.0d0

  !> If dust_temperature_type is stellar, it will be calculated according to Tielens (2005),
  !> eqn. 5.44 using an input stellar luminosity in solar luminosities
  double precision :: dust_stellar_luminosity = -1.0d0

  !> Set small dust densities to zero to avoid numerical problems
  logical, public, protected :: dust_small_to_zero = .false.

  !> Minimum dust density as used when dust_small_to_zero=T
  double precision, public, protected :: dust_min_rho = -1.0d0

  !> Adding dust in sourcesplit manner or not
  logical :: dust_source_split = .false.

  !> This can be turned off for testing purposes, if F then gas uncouples from dust
  logical :: dust_backreaction = .true.

  !> What type of dust drag force to use. Can be 'Kwok', 'sticking', 'linear', 'usr' or 'none'.
  character(len=std_len), public, protected :: dust_method = 'Kwok'

  !> Can be 'graphite' or 'silicate', affects the dust temperature
  character(len=std_len) :: dust_species = 'graphite'

  !> Determines the dust temperature, can be 'constant', 'ism', or 'stellar'
  character(len=std_len) :: dust_temperature_type = 'constant'

  ! Public methods
  public :: dust_init
  public :: dust_get_dt
  public :: dust_get_flux
  public :: dust_get_cmax
  public :: dust_get_flux_prim
  public :: dust_get_cmax_prim
  public :: dust_add_source
  public :: dust_to_conserved
  public :: dust_to_primitive
  public :: dust_check_params
  public :: dust_check_w
  public :: set_dusttozero

contains

  subroutine dust_init(g_rho, g_mom, g_energy)
    use mod_global_parameters

    integer, intent(in) :: g_rho
    integer, intent(in) :: g_mom(ndir)
    integer, intent(in) :: g_energy ! Negative value if not present
    integer             :: n, idir
    character(len=2)    :: dim

    call dust_read_params(par_files)

    allocate(gas_mom(ndir))
    gas_rho_ = g_rho
    gas_mom  = g_mom
    gas_e_   = g_energy

    allocate(dust_size(dust_n_species))
    allocate(dust_density(dust_n_species))
    dust_size(:) = -1.0d0
    dust_density(:) = -1.0d0

    allocate(dust_rho(dust_n_species))
    allocate(dust_mom(ndir, dust_n_species))

    ! Set index of dust densities
    do n = 1, dust_n_species
      dust_rho(n) = var_set_fluxvar("rhod", "rhod", n)
    end do

    ! Dust momentum
    do idir = 1, ndir
      write(dim, "(I0,A)") idir, "d"
      do n = 1, dust_n_species
        dust_mom(idir, n) = var_set_fluxvar("m"//dim, "v"//dim, n)
      end do
    end do

  end subroutine dust_init

  !> Read dust_list module parameters from a file
  subroutine dust_read_params(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /dust_list/ dust_n_species, dust_min_rho, dust_method, &
         dust_K_lineardrag, dust_small_to_zero, dust_source_split, dust_temperature, &
         dust_temperature_type, dust_backreaction, dust_dtpar, gas_vtherm_factor, dust_stellar_luminosity

    do n = 1, size(files)
      open(unitpar, file=trim(files(n)), status="old")
      read(unitpar, dust_list, end=111)
111   close(unitpar)
    end do

  end subroutine dust_read_params

  subroutine dust_check_params()
    use mod_usr_methods, only: usr_get_3d_dragforce,usr_dust_get_dt
    use mod_global_parameters, only : mype, SI_unit

    if (dust_method == 'sticking') then
       if (SI_unit) call mpistop("Dust error: sticking assumes cgs units")
       if (dust_temperature_type == "constant") then
          if (dust_temperature < 0.0d0) then
             call mpistop("Dust error: dust_temperature (in K) < 0 or not set")
          end if
       else if (dust_temperature_type == "stellar") then
          if (dust_stellar_luminosity < 0.0d0) then
             call mpistop("Dust error: dust_stellar_luminosity (in solar) < 0 or not set")
          end if
       end if
    end if

    if (dust_method == 'linear') then
        if(dust_K_lineardrag<0.0d0) then
          call mpistop("With dust_method=='linear', you must set a positive dust_K_lineardrag")
       end if
    end if

    if (any(dust_size < 0.0d0)) &
            call mpistop("Dust error: any(dust_size < 0) or not set")
    if (any(dust_density < 0.0d0)) &
            call mpistop("Dust error: any(dust_density < 0) or not set")

    if (dust_method == 'usr') then
       if (.not. associated(usr_get_3d_dragforce) .or. .not. associated(usr_dust_get_dt)) &
            call mpistop("Dust error:usr_get_3d_dragforce and usr_dust_get_dt not defined")
    end if

  end subroutine dust_check_params

  subroutine dust_check_w(ixI^L,ixO^L,w,flag)
    use mod_global_parameters
    
    integer, intent(in)         :: ixI^L,ixO^L
    double precision, intent(in):: w(ixI^S,1:nw)
    logical, intent(inout)      :: flag(ixI^S,1:nw)
    integer                     :: n

    do n = 1, dust_n_species
       flag(ixO^S,dust_rho(n))=(w(ixO^S,dust_rho(n))<0.0d0)
    enddo

  end subroutine dust_check_w

  subroutine dust_to_conserved(ixI^L, ixO^L, w, x)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, 1:nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    integer                         :: n, idir

    if(fix_small_values .and. dust_small_to_zero) call set_dusttozero(ixI^L, ixO^L, w, x)

    do n = 1, dust_n_species
      ! Convert velocity to momentum
      do idir = 1, ndir
        w(ixO^S, dust_mom(idir, n)) = w(ixO^S, dust_rho(n)) * &
             w(ixO^S, dust_mom(idir, n))
      end do
    end do

  end subroutine dust_to_conserved

  subroutine dust_to_primitive(ixI^L, ixO^L, w, x)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, 1:nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    integer                         :: n, idir

    do n = 1, dust_n_species
      ! Convert momentum to velocity
      do idir = 1, ndir
        where (w(ixO^S, dust_rho(n)) > 0.0d0)
          w(ixO^S, dust_mom(idir, n)) = w(ixO^S, dust_mom(idir, n)) / &
               w(ixO^S, dust_rho(n))
        elsewhere
          w(ixO^S, dust_mom(idir, n)) = 0.0d0
        end where
      end do
    end do

    if(fix_small_values .and. dust_small_to_zero) call set_dusttozero(ixI^L, ixO^L, w, x)

  end subroutine dust_to_primitive

  subroutine dust_get_flux(w, x, ixI^L, ixO^L, idim, f)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L, idim
    double precision, intent(in)    :: w(ixI^S, 1:nw), x(ixI^S, 1:ndim)
    double precision, intent(inout) :: f(ixI^S, nwflux)
    integer                         :: n, idir

    do n = 1, dust_n_species
       where (w(ixO^S, dust_rho(n)) > 0.0d0)
          f(ixO^S, dust_rho(n)) = w(ixO^S, dust_mom(idim, n))
       elsewhere
          f(ixO^S, dust_rho(n)) = 0.0d0
       end where

       do idir = 1, ndir
        f(ixO^S, dust_mom(idir, n)) = w(ixO^S, dust_mom(idir, n)) * &
             get_vdust(w, ixI^L, ixO^L, idim, n)
       end do
    end do

  end subroutine dust_get_flux

  subroutine dust_get_flux_prim(w, x, ixI^L, ixO^L, idim, f)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L, idim
    double precision, intent(in)    :: w(ixI^S, 1:nw), x(ixI^S, 1:ndim)
    double precision, intent(inout) :: f(ixI^S, nwflux)
    integer                         :: n, idir

    do n = 1, dust_n_species
       where (w(ixO^S, dust_rho(n)) > 0.0d0)
          f(ixO^S, dust_rho(n)) = w(ixO^S, dust_mom(idim, n))*w(ixO^S, dust_rho(n))
       elsewhere
          f(ixO^S, dust_rho(n)) = 0.0d0
       end where

       do idir = 1, ndir
        f(ixO^S, dust_mom(idir, n)) = w(ixO^S, dust_mom(idir, n)) * &
         w(ixO^S, dust_rho(n)) * get_vdust_prim(w, ixI^L, ixO^L, idim, n)
       end do
    end do

  end subroutine dust_get_flux_prim

  function get_vdust(w, ixI^L, ixO^L, idim, n) result(vdust)
    use mod_global_parameters, only: nw
    integer, intent(in)           :: ixI^L, ixO^L, idim, n
    double precision, intent(in)  :: w(ixI^S, nw)
    double precision              :: vdust(ixO^S)

    where (w(ixO^S, dust_rho(n)) > 0.0d0)
      vdust(ixO^S) = w(ixO^S, dust_mom(idim, n)) / w(ixO^S, dust_rho(n))
    elsewhere
      vdust(ixO^S) = 0.0d0
    end where

  end function get_vdust

  function get_vdust_prim(w, ixI^L, ixO^L, idim, n) result(vdust)
    use mod_global_parameters, only: nw
    integer, intent(in)           :: ixI^L, ixO^L, idim, n
    double precision, intent(in)  :: w(ixI^S, nw)
    double precision              :: vdust(ixO^S)

    where (w(ixO^S, dust_rho(n)) > 0.0d0)
      vdust(ixO^S) = w(ixO^S, dust_mom(idim, n))
    elsewhere
      vdust(ixO^S) = 0.0d0
    end where

  end function get_vdust_prim

  ! Force dust density to zero if dust_rho <= dust_min_rho
  subroutine set_dusttozero(ixI^L, ixO^L, w, x)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    double precision, intent(inout) :: w(ixI^S, 1:nw)
    logical                         :: flag(ixI^S)
    integer                         :: n, idir

    do n = 1, dust_n_species
      flag(ixO^S)=(w(ixO^S, dust_rho(n)) <= dust_min_rho)
      where (flag(ixO^S))
        w(ixO^S, dust_rho(n)) = 0.0d0
      end where
      do idir = 1, ndir
        where (flag(ixO^S))
          w(ixO^S, dust_mom(idir, n)) = 0.0d0
        end where
      end do
    end do

  end subroutine set_dusttozero

  ! Calculate drag force based on Epstein's law
  ! From Kwok 1975, page 584 (between eqn 8 and 9)
  subroutine get_3d_dragforce(ixI^L, ixO^L, w, x, fdrag, ptherm, vgas)
    use mod_global_parameters
    use mod_usr_methods, only: usr_get_3d_dragforce
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    double precision, intent(in)    :: w(ixI^S, 1:nw)
    double precision, intent(out)   :: &
         fdrag(ixI^S, 1:ndir, 1:dust_n_species)
    double precision, intent(in)    :: ptherm(ixI^S), vgas(ixI^S, 1:ndir)

    double precision, dimension(ixI^S) :: vt2, deltav, fd, vdust
    double precision                   :: alpha_T(ixI^S, 1:dust_n_species)
    integer                            :: n, idir

    vt2(ixO^S) = gas_vtherm_factor*ptherm(ixO^S)/w(ixO^S, gas_rho_)

    select case( TRIM(dust_method) )
    case ('Kwok') ! assume sticking coefficient equals 0.25

      do idir = 1, ndir
        do n = 1, dust_n_species
          where(w(ixO^S, dust_rho(n)) > 0.0d0)
            vdust(ixO^S)  = w(ixO^S, dust_mom(idir, n)) / w(ixO^S, dust_rho(n))
            deltav(ixO^S) = vgas(ixO^S, idir)-vdust(ixO^S)

            ! 0.75 from sticking coefficient
            fd(ixO^S)     = 0.75d0*w(ixO^S, dust_rho(n))*w(ixO^S, gas_rho_)*deltav(ixO^S) &
                 / (dust_density(n) * dust_size(n))

            ! 0.75 from spherical grainvolume
            fd(ixO^S)     = -fd(ixO^S)*0.75d0*dsqrt(vt2(ixO^S) + deltav(ixO^S)**2)
          elsewhere
            fd(ixO^S) = 0.0d0
          end where
          fdrag(ixO^S, idir, n) = fd(ixO^S)
        end do
      end do

    case ('sticking') ! Calculate sticking coefficient based on the gas and dust temperatures

      call get_sticking(w, x, ixI^L, ixO^L, alpha_T, ptherm)

      do idir = 1, ndir
        do n = 1, dust_n_species
          where(w(ixO^S, dust_rho(n))>0.0d0)
            vdust(ixO^S)  = w(ixO^S,dust_mom(idir, n)) / w(ixO^S, dust_rho(n))
            deltav(ixO^S) = vgas(ixO^S, idir)-vdust(ixO^S)
            fd(ixO^S)     = (one-alpha_T(ixO^S,n)) * w(ixO^S, dust_rho(n))*w(ixO^S, gas_rho_) * &
                 deltav(ixO^S) / (dust_density(n)*dust_size(n))
            fd(ixO^S)     = -fd(ixO^S)*0.75d0*dsqrt(vt2(ixO^S) + deltav(ixO^S)**2)
          else where
            fd(ixO^S) = 0.0d0
          end where
          fdrag(ixO^S, idir,n) = fd(ixO^S)
        end do
      end do

    case('linear') !linear with Deltav, for testing (see Laibe & Price 2011)
      do idir = 1, ndir
        do n = 1, dust_n_species
          where(w(ixO^S, dust_rho(n))>0.0d0)
            vdust(ixO^S)  = w(ixO^S,dust_mom(idir, n))/w(ixO^S, dust_rho(n))
            deltav(ixO^S) = vgas(ixO^S, idir)-vdust(ixO^S)

            fd(ixO^S)     = -dust_K_lineardrag*deltav(ixO^S)
          else where
            fd(ixO^S) = 0.0d0
          end where
          fdrag(ixO^S, idir,n) = fd(ixO^S)
        end do
      end do

    case('usr')
      call usr_get_3d_dragforce(ixI^L, ixO^L, w, x, fdrag, ptherm, vgas, dust_n_species)
    case('none')
      fdrag(ixO^S, :, :) = 0.0d0
    case default
      call mpistop( "=== This dust method has not been implemented===" )
    end select

  end subroutine get_3d_dragforce

  !> Get sticking coefficient alpha_T (always between 0 and 1)
  !>
  !> Uses Temperatures in K
  !> Equation from Decin et al. 2006
  subroutine get_sticking(w, x, ixI^L, ixO^L, alpha_T, ptherm)
    use mod_global_parameters
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: x(ixI^S, 1:ndim)
    double precision, intent(in)  :: w(ixI^S, 1:nw)
    double precision, intent(out) :: alpha_T(ixI^S, 1:dust_n_species)
    double precision, intent(in)  :: ptherm(ixI^S)
    double precision              :: Tgas(ixI^S)
    integer                       :: n

    ! get the dust species T in K
    call get_tdust(w, x, ixI^L, ixO^L, alpha_T)

    ! convert dimensionless gas T to K
    Tgas(ixO^S) = (ptherm(ixO^S)/w(ixO^S, gas_rho_))*unit_temperature

    do n = 1, dust_n_species
      alpha_T(ixO^S,n) =  0.35d0 * dexp(-dsqrt((Tgas(ixO^S) + &
           alpha_T(ixO^S,n))/5.0d2))+0.1d0
    end do

  end subroutine get_sticking

  !> Returns dust temperature (in K), either as constant or based on equ. 5.41,
  !> 5.42 and 5.44 from Tielens (2005)
  !>
  !> Note that this calculation assumes cgs!!!! 
  !>
  !> It takes as input the stellar luminosity in solar units in 'stellar' case
  !> or a fixed dust input temperature in Kelvin when 'constant' or does case 'ism'
  subroutine get_tdust(w, x, ixI^L, ixO^L, Td)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: x(ixI^S, 1:ndim)
    double precision, intent(in)  :: w(ixI^S, 1:nw)
    double precision, intent(out) :: Td(ixI^S, 1:dust_n_species)
    double precision              :: G0(ixO^S)
    integer                       :: n

    select case( trim(dust_temperature_type) )
    case( 'constant' )
      Td(ixO^S, :) = dust_temperature
    case( 'ism' )
      select case( trim(dust_species) )
      case( 'graphite' )
        do n = 1, dust_n_species
          Td(ixO^S, n) = 15.8d0*((0.0001d0/(dust_size(n)*unit_length))**0.06d0)
        end do
      case( 'silicate' )
        do n = 1, dust_n_species
          Td(ixO^S, n) = 13.6d0*((0.0001d0/(dust_size(n)*unit_length))**0.06d0)
        end do
      case default
        call mpistop( "=== Dust species undetermined===" )
      end select
    case( 'stellar' )
      select case(coordinate)
      case(spherical)
        G0(ixO^S) = max(x(ixO^S, 1)*unit_length, smalldouble)
      !!!case(cylindrical) convert R,Z to spherical radial coordinate r here
      !!! but only ok for 2D (R,Z) or 2.5D (R,Z) case
      !!! G0(ixO^S) = max(dsqrt(sum(x(ixO^S,:)**2,dim=ndim+1))*unit_length, smalldouble)
      case default
        call mpistop('stellar case not available in this coordinate system')
      end select

      G0(ixO^S) = 2.1d4*(dust_stellar_luminosity/1.0d8)*((3.0857d17/G0(ixO^S))**2)

      select case( trim(dust_species) )
      case( 'graphite' )
        do n = 1, dust_n_species
          Td(ixO^S, n) = 61.0d0*((0.0001d0/(dust_size(n)*unit_length))**0.06d0) &
               *(G0(ixO^S)**(one/5.8d0))
        end do
      case( 'silicate' )
        do n = 1, dust_n_species
          Td(ixO^S, n) = 50.0d0*((0.0001d0/(dust_size(n)*unit_length))**0.06d0) &
               *(G0(ixO^S)**(one/6.0d0))
        end do
      case default
        call mpistop( "=== Dust species undetermined===" )
      end select
    case default
      call mpistop( "=== Dust temperature undetermined===" )
    end select

  end subroutine get_tdust

  !> w[iw]= w[iw]+qdt*S[wCT,  x] where S is the source based on wCT within ixO
  subroutine dust_add_source(qdt, ixI^L, ixO^L, wCT, w, x, qsourcesplit, active)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixI^S, 1:nw), x(ixI^S, 1:ndim)
    double precision, intent(inout) :: w(ixI^S, 1:nw)
    logical, intent(in)             :: qsourcesplit
    logical, intent(inout)          :: active

    double precision :: ptherm(ixI^S), vgas(ixI^S, 1:ndir)
    double precision :: fdrag(ixI^S, 1:ndir, 1:dust_n_species)
    integer          :: n, idir

    select case( TRIM(dust_method) )
    case( 'none' )
      !do nothing here
    case default !all regular dust methods here
      if (qsourcesplit .eqv. dust_source_split) then
        active = .true.

        call phys_get_pthermal(wCT, x, ixI^L, ixO^L, ptherm)
        do idir=1,ndir
          vgas(ixO^S,idir)=wCT(ixO^S,gas_mom(idir))/wCT(ixO^S,gas_rho_)
        end do

        call get_3d_dragforce(ixI^L, ixO^L, wCT, x, fdrag, ptherm, vgas)
        fdrag(ixO^S, 1:ndir, 1:dust_n_species) = fdrag(ixO^S, 1:ndir, 1:dust_n_species) * qdt

        do idir = 1, ndir

          do n = 1, dust_n_species
            if (dust_backreaction) then
               w(ixO^S, gas_mom(idir))  = w(ixO^S, gas_mom(idir)) + &
                    fdrag(ixO^S, idir, n)
            end if

            if (gas_e_ > 0) then
              w(ixO^S, gas_e_) = w(ixO^S, gas_e_) + vgas(ixO^S, idir)  &
                   * fdrag(ixO^S, idir, n)
            end if

            w(ixO^S, dust_mom(idir, n)) = w(ixO^S, dust_mom(idir, n)) - &
                 fdrag(ixO^S, idir, n)
          end do
        end do

      endif
    end select

  end subroutine dust_add_source

  !> Get dt related to dust and gas stopping time (Laibe 2011)
  subroutine dust_get_dt(w, ixI^L, ixO^L, dtnew, dx^D, x)
    use mod_global_parameters
    use mod_usr_methods, only: usr_dust_get_dt

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: dx^D, x(ixI^S, 1:ndim)
    double precision, intent(in)    :: w(ixI^S, 1:nw)
    double precision, intent(inout) :: dtnew

    double precision                :: ptherm(ixI^S), vgas(ixI^S, 1:ndir)
    double precision, dimension(1:dust_n_species):: dtdust
    double precision, dimension(ixI^S)           :: vt2, deltav, tstop, vdust
    double precision, dimension(ixI^S, 1:dust_n_species) :: alpha_T
    integer                                    :: n, idir

    call phys_get_pthermal(w, x, ixI^L, ixO^L, ptherm)
    do idir = 1, ndir
      vgas(ixO^S,idir)=w(ixO^S,gas_mom(idir))/w(ixO^S,gas_rho_)
    end do

    select case( TRIM(dust_method) )

    case( 'Kwok' ) ! assume sticking coefficient equals 0.25
      dtdust(:) = bigdouble

      vt2(ixO^S) = gas_vtherm_factor*ptherm(ixO^S)/w(ixO^S, gas_rho_)

      do idir = 1, ndir
        do n = 1, dust_n_species
          where(w(ixO^S, dust_rho(n))>0.0d0)
            vdust(ixO^S)  = w(ixO^S,dust_mom(idir, n))/w(ixO^S, dust_rho(n))
            deltav(ixO^S) = vgas(ixO^S, idir)-vdust(ixO^S)
            tstop(ixO^S)  = 4.0d0*(dust_density(n)*dust_size(n))/ &
                 (3.0d0*(0.75d0)*dsqrt(vt2(ixO^S) + &
                 deltav(ixO^S)**2)*(w(ixO^S, dust_rho(n)) + &
                 w(ixO^S, gas_rho_)))
          else where
            tstop(ixO^S) = bigdouble
          end where

          dtdust(n) = min(minval(tstop(ixO^S)), dtdust(n))
        end do
      end do

      dtnew = min(minval(dust_dtpar*dtdust(:)), dtnew)

    case( 'sticking' ) ! Calculate sticking coefficient based on the gas temperature
      dtdust(:) = bigdouble

      vt2(ixO^S) = gas_vtherm_factor*ptherm(ixO^S)/w(ixO^S, gas_rho_)

      ! Sticking coefficient
      call get_sticking(w, x, ixI^L, ixO^L, alpha_T, ptherm)

      do idir = 1, ndir
        do n = 1, dust_n_species
          where(w(ixO^S, dust_rho(n))>0.0d0)
            vdust(ixO^S)  = w(ixO^S,dust_mom(idir, n))/w(ixO^S, dust_rho(n))
            deltav(ixO^S) = vgas(ixO^S, idir)-vdust(ixO^S)
            tstop(ixO^S)  = 4.0d0*(dust_density(n)*dust_size(n))/ &
                 (3.0d0*(one-alpha_T(ixO^S,n))*dsqrt(vt2(ixO^S) + &
                 deltav(ixO^S)**2)*(w(ixO^S, dust_rho(n)) + &
                 w(ixO^S, gas_rho_)))
          else where
            tstop(ixO^S) = bigdouble
          end where

          dtdust(n) = min(minval(tstop(ixO^S)), dtdust(n))
        end do
      end do

      dtnew = min(minval(dust_dtpar*dtdust(:)), dtnew)

    case('linear') !linear with Deltav, for testing (see Laibe & Price 2011)
      dtdust(:) = bigdouble

      do n = 1, dust_n_species
        where(w(ixO^S, dust_rho(n))>0.0d0)
          tstop(ixO^S)  = (w(ixO^S, dust_rho(n))*w(ixO^S, gas_rho_))/ &
               (dust_K_lineardrag*(w(ixO^S, dust_rho(n)) + w(ixO^S, gas_rho_)))
        else where
          tstop(ixO^S) = bigdouble
        end where

        dtdust(n) = min(minval(tstop(ixO^S)), dtdust(n))
      end do

      dtnew = min(minval(dust_dtpar*dtdust(:)), dtnew)
    case('usr')
      dtdust(:) = bigdouble
      call usr_dust_get_dt(w, ixI^L, ixO^L, dtdust, dx^D, x, dust_n_species)
      dtnew = min(minval(dust_dtpar*dtdust(:)), dtnew)
    case('none')
      ! no dust timestep
    case default
      call mpistop( "=== This dust method has not been implemented===" )
    end select

    if (dtnew < dtmin) then
      write(unitterm,*)"-------------------------------------"
      write(unitterm,*)"Warning: found DUST related time step too small! dtnew=", dtnew
      write(unitterm,*)"on grid with index:", block%igrid," grid level=", block%level
      write(unitterm,*)"grid corners are=",{^D&rnode(rpxmin^D_, block%igrid), rnode(rpxmax^D_, block%igrid)}
      write(unitterm,*)" dtdust =", dtdust(:)
      write(unitterm,*)"on processor:", mype
      write(unitterm,*)"-------------------------------------"
    endif

  end subroutine dust_get_dt

  ! Note that cmax and cmin are assumed to be initialized
  subroutine dust_get_cmax(w, x, ixI^L, ixO^L, idim, cmax, cmin)
    use mod_global_parameters
    use mod_variables

    integer, intent(in)                       :: ixI^L, ixO^L, idim
    double precision, intent(in)              :: w(ixI^S, 1:nw), x(ixI^S, 1:ndim)
    double precision, intent(inout)           :: cmax(ixI^S,1:number_species)
    double precision, intent(inout), optional :: cmin(ixI^S,1:number_species)
    double precision                          :: vdust(ixO^S)
    integer                                   :: n

    do n = 1, dust_n_species
      vdust(ixO^S) = get_vdust(w, ixI^L, ixO^L, idim, n)

      if (present(cmin)) then
        cmin(ixO^S,1) = min(cmin(ixO^S,1), vdust(ixO^S))
        cmax(ixO^S,1) = max(cmax(ixO^S,1), vdust(ixO^S))
      else
        cmax(ixO^S,1) = max(cmax(ixO^S,1), abs(vdust(ixO^S)))
      end if
    end do
  end subroutine dust_get_cmax

  ! Note that cmax and cmin are assumed to be initialized
  subroutine dust_get_cmax_prim(w, x, ixI^L, ixO^L, idim, cmax, cmin)
    use mod_global_parameters

    integer, intent(in)                       :: ixI^L, ixO^L, idim
    double precision, intent(in)              :: w(ixI^S, 1:nw), x(ixI^S, 1:ndim)
    double precision, intent(inout)           :: cmax(ixI^S)
    double precision, intent(inout), optional :: cmin(ixI^S)
    double precision                          :: vdust(ixO^S)
    integer                                   :: n

    do n = 1, dust_n_species
      vdust(ixO^S) = get_vdust_prim(w, ixI^L, ixO^L, idim, n)

      if (present(cmin)) then
        cmin(ixO^S) = min(cmin(ixO^S), vdust(ixO^S))
        cmax(ixO^S) = max(cmax(ixO^S), vdust(ixO^S))
      else
        cmax(ixO^S) = max(cmax(ixO^S), abs(vdust(ixO^S)))
      end if
    end do
  end subroutine dust_get_cmax_prim

end module mod_dust
