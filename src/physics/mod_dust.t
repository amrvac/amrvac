module mod_dust
  use mod_global_parameters, only: std_len
  use mod_physics

  implicit none
  private

  !> The number of dust species
  integer, public, protected :: dust_n_species

  integer, protected              :: gas_rho_
  integer, allocatable, protected :: gas_mom(:)
  integer, protected              :: gas_e_
  double precision, protected     :: gas_mu = -huge(1.0d0)

  !> Indices of the dust densities
  integer, allocatable, public, protected :: dust_rho(:)

  !> Indices of the dust momentum densities
  integer, allocatable, public, protected :: dust_mom(:, :)

  !> Size of each dust species
  double precision, allocatable, public :: dust_size(:)

  !> Internal density of each dust species
  double precision, allocatable, public :: dust_density(:)

  !> Dust temperature (if dust_temperature_type is constant)
  double precision :: dust_temperature

  !> If dust_temperature_type is stellar, it will be calculated according to Tielens (2005),
  !> eqn. 5.44 using a stellar luminosity in solar luminosities
  double precision :: dust_stellar_luminosity

  ! ???
  double precision :: mhcgspar = 1.6733D-24
  double precision :: kbcgspar = 1.38065D-16

  !> Set small dust densities to zero to avoid numerical problems
  logical :: dust_small_to_zero = .false.

  !> Minimum dust density
  double precision :: dust_min_rho = -1.0d0

  !> TODO: 1. Introduce this generically in advance, 2: document
  logical :: dust_source_split = .false.

  !> What type of dust drag force to use. Can be 'Kwok', 'sticking', 'linear',or 'none'.
  character(len=std_len) :: dust_method = 'Kwok'

  !> Can be 'graphite' or 'silicate', affects the dust temperature
  character(len=std_len) :: dust_species = 'graphite'

  !> Determines the dust temperature, can be 'constant', 'ism', or 'stellar'
  character(len=std_len) :: dust_temperature_type = 'constant'

  ! Public methods
  public :: dust_init
  public :: dust_get_dt
  public :: dust_get_flux
  public :: dust_get_cmax
  public :: dust_add_source

contains

  subroutine dust_init(g_rho, g_mom, g_energy)
    use mod_global_parameters
    integer, intent(in) :: g_rho
    integer, intent(in) :: g_mom(ndir)
    integer, intent(in) :: g_energy ! Negative value if not present
    integer             :: n, idir

    dust_n_species=1

    call dust_read_params(par_files)

    allocate(gas_mom(ndir))
    gas_rho_ = g_rho
    gas_mom  = g_mom
    gas_e_   = g_energy

    ! TODO: how will dust read settings? Use m_config?

    allocate(dust_size(dust_n_species))
    allocate(dust_density(dust_n_species))
    allocate(dust_rho(dust_n_species))
    allocate(dust_mom(ndir, dust_n_species))

    ! Set starting index of dust species
    do n = 1, dust_n_species
       nwflux = nwflux + 1
       dust_rho(n) = nwflux     ! Dust density

       do idir = 1, ndir
          nwflux = nwflux + 1
          dust_mom(idir, n) = nwflux ! Dust momentum
       end do
    end do

  end subroutine dust_init

  !> Read this module"s parameters from a file
  subroutine dust_read_params(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /dust_list/ dust_n_species, dust_min_rho

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, dust_list, end=111)
111    close(unitpar)
    end do

  end subroutine dust_read_params

  subroutine dust_check_params()
    if (gas_mu <= 0.0d0) call mpistop ("Error: mu (molecular weight) negative")

  end subroutine dust_check_params

  subroutine dust_to_conserved(ixI^L, ixO^L, w, x)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    integer                         :: n, idir

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
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    integer                         :: n, idir

    do n = 1, dust_n_species
       ! Convert momentum to velocity
       do idir = 1, ndir
          where (w(ixO^S, dust_rho(n)) > dust_min_rho)
             w(ixO^S, dust_mom(idir, n)) = w(ixO^S, dust_rho(n)) / &
                  w(ixO^S, dust_mom(idir, n))
          elsewhere
             w(ixO^S, dust_mom(idir, n)) = 0
          end where
       end do
    end do
  end subroutine dust_to_primitive

  subroutine dust_get_flux(w, x, ixI^L, ixO^L, iw, idim, f, transport)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L, iw, idim
    double precision, intent(in)    :: w(ixI^S, 1:nw), x(ixI^S, 1:^ND)
    double precision, intent(inout) :: f(ixI^S)
    double precision                :: vdust(ixI^S)
    logical, intent(out)            :: transport
    integer                         :: n

    transport = .false.
    f(ixO^S) = 0.0d0

    do n = 1, dust_n_species
       if (iw == dust_rho(n)) then ! A dust density
          where (w(ixO^S, dust_rho(n)) > dust_min_rho)
             f(ixO^S) = w(ixO^S, dust_mom(idim, n))
          elsewhere             ! TODO: remove?
             f(ixO^S) = 0.0d0
          end where
          exit
       else if (iw == dust_mom(idim, n)) then ! A dust momentum component
          f(ixO^S) = w(ixO^S, iw) * get_vdust(w, ixI^L, ixO^L, idim, n)
          exit
       end if
    end do
  end subroutine dust_get_flux

  function get_vdust(w, ixI^L, ixO^L, idim, n) result(vdust)
    use mod_global_parameters, only: nw
    integer, intent(in)           :: ixI^L, ixO^L, idim, n
    double precision, intent(in)  :: w(ixI^S, nw)
    double precision              :: vdust(ixO^S)

    where (w(ixO^S, dust_rho(n)) > dust_min_rho)
       vdust = w(ixO^S, dust_mom(idim, n)) / w(ixO^S, dust_rho(n))
    elsewhere
       vdust = 0.0d0;
    end where
  end function get_vdust

  ! Force dust density to zero if dust_rho <= dust_min_rho
  subroutine set_dusttozero(qdt, ixI^L, ixO^L,  wCT,  w, x)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixI^S, 1:nw), x(ixI^S, 1:ndim)
    double precision, intent(inout) :: w(ixI^S, 1:nw)
    integer                         :: n, idir

    do n = 1, dust_n_species
       where (w(ixO^S, dust_rho(n)) <= dust_min_rho)
          w(ixO^S, dust_rho(n)) = 0.0d0
       end where

       do idir = 1, ndir
          where (w(ixO^S, dust_rho(n)) <= dust_min_rho)
             w(ixO^S, dust_mom(idir, n)) = 0.0d0
          end where
       end do
    end do
  end subroutine set_dusttozero

  ! Calculate drag force based on Epstein's or Stokes' law
  ! From Kwok 1975, page 584 (between eqn 8 and 9)
  subroutine get_3d_dragforce(ixI^L, ixO^L, w, x, fdrag, ptherm, vgas)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    double precision, intent(in)    :: w(ixI^S, 1:nw)
    double precision, intent(out)   :: &
         fdrag(ixI^S, 1:ndir, 1:dust_n_species)
    double precision, intent(in)    :: ptherm(ixI^S), vgas(ixI^S, ndir)

    double precision, dimension(ixI^S) :: vt2, deltav, fd, vdust, Tgas
    double precision                   :: alpha_T(ixI^S, 1:dust_n_species)
    integer                            :: n, idir
    double precision                   :: K

    ! call hd_get_pthermal(w, x, ixI^L, ixO^L, ptherm)

    vt2(ixO^S) = 3.0d0*ptherm(ixO^S)/w(ixO^S, gas_rho_)

    select case( TRIM(dust_method) )
    case ('Kwok') ! assume sticking coefficient equals 0.25

       do idir = 1, ndir
          ! call hd_get_v(w, x, ixI^L, ixO^L, idir, vgas)

          do n = 1, dust_n_species
             ! TODO: simplify
             where(w(ixO^S, dust_rho(n)) > dust_min_rho)
                vdust(ixO^S)  = w(ixO^S, dust_mom(idir, n)) / w(ixO^S, dust_rho(n))
                deltav(ixO^S) = (vgas(ixO^S, idir)-vdust(ixO^S))

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
       !  Equation from Decin et al. 2006
       if (gas_e_ < 0) call mpistop("dust sticking requires gas energy")

       Tgas(ixO^S) = ( ptherm(ixO^S)*normvar(gas_e_)*mhcgspar) / &
            (w(ixO^S, gas_rho_)*normvar(gas_rho_)*kbcgspar)
       call get_sticking(w, x, ixI^L, ixO^L, alpha_T, ptherm)

       do idir = 1, ndir
          ! call hd_get_v(w, x, ixI^L, ixO^L, idir, vgas)

          do n = 1, dust_n_species
             where(w(ixO^S, dust_rho(n))>dust_min_rho)
                vdust(ixO^S)  = w(ixO^S,dust_mom(idir, n)) / w(ixO^S, dust_rho(n))
                deltav(ixO^S) = (vgas(ixO^S, idir)-vdust(ixO^S))
                fd(ixO^S)     = (one-alpha_T(ixO^S,n)) * w(ixO^S, dust_rho(n))*w(ixO^S, gas_rho_) * &
                     deltav(ixO^S) / (dust_density(n)*dust_size(n))
                fd(ixO^S)     = -fd(ixO^S) * 0.75d0 * dsqrt(vt2(ixO^S) + deltav(ixO^S)**2)
             else where
                fd(ixO^S) = 0.0d0
             end where
             fdrag(ixO^S, idir,n) = fd(ixO^S)
          end do
       end do
    case('linear') !linear with Deltav, for testing (see Laibe & Price 2011)
       K = 3.4d5 / dust_n_species
       do idir = 1, ndir
          ! call hd_get_v(w, x, ixI^L, ixO^L, idir, vgas)

          do n = 1, dust_n_species
             where(w(ixO^S, dust_rho(n))>dust_min_rho)
                vdust(ixO^S)  = w(ixO^S,dust_mom(idir, n))/w(ixO^S, dust_rho(n))
                deltav(ixO^S) = (vgas(ixO^S, idir)-vdust(ixO^S))

                fd(ixO^S)     = -K*deltav(ixO^S)
             else where
                fd(ixO^S) = 0.0d0
             end where
             fdrag(ixO^S, idir,n) = fd(ixO^S)
          end do
       end do
    case('none')
       fdrag(ixO^S, idir, :) = 0.0d0
    case default
       call mpistop( "=== This dust method has not been implemented===" )
    end select

  end subroutine get_3d_dragforce

  !> get sticking coefficient
  !>
  !> Assume cgs units, and use of normvar(0:nw) array for conversion
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

    ! call getpthermal(w, x, ixI^L, ixO^L, Tgas)
    call get_tdust(w, x, ixI^L, ixO^L, alpha_T)

    Tgas(ixO^S) = (ptherm(ixO^S)*normvar(gas_e_)*mhcgspar) / &
         (w(ixO^S, gas_rho_) * normvar(gas_rho_) * kbcgspar)

    do n = 1, dust_n_species
       alpha_T(ixO^S,n) =  max(0.35d0 * exp(-dsqrt((Tgas(ixO^S) + &
            alpha_T(ixO^S,n))/5.0d2))+0.1d0, smalldouble)
    end do
  end subroutine get_sticking

  !> Returns dust temperature (in K), either as constant or based on equ. 5.41,
  !> 5.42 and 5.44 from Tielens (2005)
  !>
  !> Note that this calculation assumes cgs!!!! with conversion between physical
  !> and scaled quantities done through the normvar(0:nw) array!!!!
  !>
  !> It takes as input the stellar luminosoity in solar units and/or a fixed
  !> dust temperature in Kelvin
  subroutine get_tdust(w, x, ixI^L, ixO^L, Td)
    use mod_global_parameters

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
             Td(ixO^S, n) = 15.8d0*((0.0001d0/(dust_size(n)*normvar(0)))**0.06d0)
          end do
       case( 'silicate' )
          do n = 1, dust_n_species
             Td(ixO^S, n) = 13.6d0*((0.0001d0/(dust_size(n)*normvar(0)))**0.06d0)
          end do
       case default
          call mpistop( "=== Dust species undetermined===" )
       end select
    case( 'stellar' )
       select case( trim(typeaxial) )
       case( 'spherical' )
          G0(ixO^S) = max(x(ixO^S, 1)*normvar(0), smalldouble)
       case( 'cylindrical' )
          G0(ixO^S) = max(dsqrt(sum(x(ixO^S,:)**2,dim=ndim+1))*normvar(0), smalldouble)
       case( 'slab' )
          {^IFTHREED
          G0(ixO^S) = max(dsqrt((x(ixO^S, 1)-x1ptms)**2 + (x(ixO^S, 2)-x2ptms)**2  &
               + (x(ixO^S, 3)-x3ptms)**2)*normvar(0), smalldouble)
          }
       end select

       G0(ixO^S) = 2.1d4*(dust_stellar_luminosity/1.0d8)*((3.0857d17/G0(ixO^S))**2)

       select case( trim(dust_species) )
       case( 'graphite' )
          do n = 1, dust_n_species
             Td(ixO^S, n) = 61.0d0*((0.0001d0/(dust_size(n)*normvar(0)))**0.06d0) &
                  *(G0(ixO^S)**(one/5.8d0))
          end do
       case( 'silicate' )
          do n = 1, dust_n_species
             Td(ixO^S, n) = 50.0d0*((0.0001d0/(dust_size(n)*normvar(0)))**0.06d0) &
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
  subroutine dust_add_source(qdt, ixI^L, ixO^L, wCT,w, x, qsourcesplit)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixI^S, 1:nw), x(ixI^S, 1:ndim)
    double precision, intent(inout) :: w(ixI^S, 1:nw)
    logical, intent(in)             :: qsourcesplit

    double precision                :: ptherm(ixI^S), vgas(ixI^S, ndir)
    double precision, dimension(ixI^S, 1:ndir, 1:dust_n_species) :: fdrag
    integer                                             :: n, idir, sum_dim

    sum_dim = ndim + 2          ! Temporary variable, used below

    select case( TRIM(dust_method) )
    case( 'none' )
       !do nothing here
    case default !all regular dust methods here
       if (qsourcesplit .eqv. dust_source_split) then
          call phys_get_pthermal(wCT, x, ixI^L, ixO^L, ptherm)
          do idir=1,ndir
            vgas(ixO^S,idir)=wCT(ixO^S,gas_mom(idir))/wCT(ixO^S,gas_rho_)
          end do
          call get_3d_dragforce(ixI^L, ixO^L, wCT, x, fdrag, ptherm, vgas)

          do idir = 1, ndir
             fdrag(ixO^S, idir, :) = fdrag(ixO^S, idir, :) * qdt

             w(ixO^S, gas_mom(idir))  = w(ixO^S, gas_mom(idir)) + &
                  sum(fdrag(ixO^S, idir, :), dim=sum_dim)

             if (gas_e_ > 0) then
                w(ixO^S, gas_e_) = w(ixO^S, gas_e_) + (wCT(ixO^S, gas_mom(idir)) / &
                     wCT(ixO^S, gas_rho_)) * sum(fdrag(ixO^S, idir, :), dim=sum_dim)
             end if
             w(ixO^S,dust_mom(idir, :)) = w(ixO^S,dust_mom(idir, :)) - fdrag(ixO^S, idir, :)
          end do

          if ( dust_small_to_zero ) call set_dusttozero(qdt, ixI^L, ixO^L,  wCT,  w, x)
       endif
    end select

  end subroutine dust_add_source

  !> Get dt related to dust and gas stopping time (Laibe 2011)
  subroutine dust_get_dt(w, ixI^L, ixO^L, dtnew, dx^D, x)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: dx^D, x(ixI^S, 1:^ND)
    double precision, intent(in)    :: w(ixI^S, 1:nw)
    double precision, intent(inout) :: dtnew

    double precision                :: ptherm(ixI^S), vgas(ixI^S, ndir)
    double precision, dimension(1:dust_n_species)        :: dtdust
    double precision, dimension(ixI^S)         :: vt2, deltav, tstop, vdust
    double precision, dimension(ixI^S, 1:dust_n_species) :: alpha_T
    double precision                           :: K
    integer                                    :: n, idims, idust, idir

    call phys_get_pthermal(w, x, ixI^L, ixO^L, ptherm)
    do idir = 1, ndir
      vgas(ixO^S,idir)=w(ixO^S,gas_mom(idir))/w(ixO^S,gas_rho_)
    end do
    select case( TRIM(dust_method) )

    case( 'Kwok' ) ! assume sticking coefficient equals 0.25
       dtdust(:) = bigdouble

       vt2(ixO^S) = 3.0d0*ptherm(ixO^S)/w(ixO^S, gas_rho_)

       ! Tgas, mu = mean molecular weight
       ! ptherm(ixO^S) = ( ptherm(ixO^S)*normvar(gas_e_)*mhcgspar*gas_mu)/(w(ixO^S, gas_rho_)*normvar(gas_rho_)*kbcgspar)

       do idir = 1, ndir
          ! call hd_get_v(w, x, ixI^L, ixO^L, idir, vgas)

          do n = 1, dust_n_species
             where(w(ixO^S, dust_rho(n))>dust_min_rho)
                vdust(ixO^S)  = w(ixO^S,dust_mom(idir, n))/w(ixO^S, dust_rho(n))
                deltav(ixO^S) = (vgas(ixO^S, idir)-vdust(ixO^S))
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

       dtnew = min(minval(dtdiffpar*dtdust(:)), dtnew)

    case( 'sticking' ) ! Calculate sticking coefficient based on the gas temperature
       dtdust(:) = bigdouble

       vt2(ixO^S) = 3.0d0*ptherm(ixO^S)/w(ixO^S, gas_rho_)

       ! Sticking coefficient
       call get_sticking(w, x, ixI^L, ixO^L, alpha_T, ptherm)

       ! Tgas, mu = mean molecular weight
       ! ptherm(ixO^S) = ( ptherm(ixO^S)*normvar(gas_e_) * mhcgspar*gas_mu) / &
            ! (w(ixO^S, gas_rho_)*normvar(gas_rho_)*kbcgspar)

       do idir = 1, ndir
          ! call hd_get_v(w, x, ixI^L, ixO^L, idir, vgas)

          do n = 1, dust_n_species
             where(w(ixO^S, dust_rho(n))>dust_min_rho)
                vdust(ixO^S)  = w(ixO^S,dust_mom(idir, n))/w(ixO^S, dust_rho(n))
                deltav(ixO^S) = (vgas(ixO^S, idir)-vdust(ixO^S))
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

       dtnew = min(minval(dtdiffpar*dtdust(:)), dtnew)

    case('linear') !linear with Deltav, for testing (see Laibe & Price 2011)
       K = 3.4d5/dust_n_species
       dtdust(:) = bigdouble

       do n = 1, dust_n_species
          where(w(ixO^S, dust_rho(n))>dust_min_rho)
             tstop(ixO^S)  = (w(ixO^S, dust_rho(n))*w(ixO^S, gas_rho_))/ &
                  (K*(w(ixO^S, dust_rho(n)) + w(ixO^S, gas_rho_)))
          else where
             tstop(ixO^S) = bigdouble
          end where

          dtdust(n) = min(minval(tstop(ixO^S)), dtdust(n))
       end do

       dtnew = min(minval(dtdiffpar*dtdust(:)), dtnew)
    case('none')
       ! no dust timestep
    case default
       call mpistop( "=== This dust method has not been implemented===" )
    end select

    if (dtnew < dtmin) then
       write(unitterm,*)"-------------------------------------"
       write(unitterm,*)"Warning: found DUST related time step too small! dtnew=", dtnew
       write(unitterm,*)"on grid with index:", saveigrid," grid level=", node(plevel_, saveigrid)
       write(unitterm,*)"grid corners are=",{^D&rnode(rpxmin^D_, saveigrid), rnode(rpxmax^D_, saveigrid)}
       write(unitterm,*)" dtdust =", dtdust(:)
       write(unitterm,*)"on processor:", mype
       write(unitterm,*)"-------------------------------------"
    endif

  end subroutine dust_get_dt

  ! Note that cmax and cmin are assumed to be initialized
  subroutine dust_get_cmax(w, x, ixI^L, ixO^L, idim, cmax, cmin)
    use mod_global_parameters

    integer, intent(in)                       :: ixI^L, ixO^L, idim
    double precision, intent(in)              :: w(ixI^S, nw), x(ixI^S, 1:^ND)
    double precision, intent(inout)           :: cmax(ixI^S)
    double precision, intent(inout), optional :: cmin(ixI^S)
    double precision                          :: vdust(ixI^S, dust_n_species)
    integer                                   :: n

    do n = 1, dust_n_species
       vdust(ixO^S, n) = get_vdust(w, ixI^L, ixO^L, idim, n)

       if (present(cmin)) then
          cmin(ixO^S) = min(cmin(ixO^S), vdust(ixO^S, n))
          cmax(ixO^S) = max(cmax(ixO^S), vdust(ixO^S, n))
       else
          cmax(ixO^S) = max(cmax(ixO^S), abs(vdust(ixO^S, n)))
       end if
    end do
  end subroutine dust_get_cmax

end module mod_dust
