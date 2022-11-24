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

  !> whether second order terms (relevant only when dust_n_species >=2) are included
  !> there are the terms  n2, ni2, d2 in Eqs 6,7,8 in amrvac 3.0 paper
  logical :: dust_implicit_second_order = .true.  

  !> whether fh is added for gas energy:  is only added in the impliict implementation, the explicit one was left as before
  logical :: dust_backreaction_fh = .false.  


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
  public :: dust_implicit_update
  public :: dust_evaluate_implicit


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
         dust_temperature_type, dust_backreaction, dust_dtpar, gas_vtherm_factor, dust_stellar_luminosity,&
         dust_implicit_second_order, dust_backreaction_fh 

    do n = 1, size(files)
      open(unitpar, file=trim(files(n)), status="old")
      read(unitpar, dust_list, end=111)
111   close(unitpar)
    end do

  end subroutine dust_read_params

  subroutine dust_check_params()
    use mod_usr_methods, only: usr_get_3d_dragforce,usr_dust_get_dt
    use mod_global_parameters, only : mype, SI_unit, use_imex_scheme

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

    if(.not. use_imex_scheme .and. ((dust_dtpar .ge. 1d0).or.(dust_dtpar.le.0))) then
      if(mype .eq. 0) print*, "EXPLICIT source for dust requires 0<dt_dustpar < 1, set to 0.8"
      dust_dtpar = 0.8
    endif

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
               if (gas_e_ > 0) then
                 w(ixO^S, gas_e_) = w(ixO^S, gas_e_) + vgas(ixO^S, idir)  &
                     * fdrag(ixO^S, idir, n)
               end if
            end if


            w(ixO^S, dust_mom(idir, n)) = w(ixO^S, dust_mom(idir, n)) - &
                 fdrag(ixO^S, idir, n)
          end do
        end do

      endif
    end select

  end subroutine dust_add_source

  !> inplace update of psa==>F_im(psa)
  subroutine dust_evaluate_implicit(qtC,psa)
    use mod_global_parameters
    type(state), target :: psa(max_blocks)
    double precision, intent(in) :: qtC

    integer :: iigrid, igrid, level

    !dust_method = 'none' not used

    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
      block=>psa(igrid)
       call dust_terms(ixG^LL,ixM^LL,psa(igrid)%w,psa(igrid)%x)
    end do
    !$OMP END PARALLEL DO

  end subroutine dust_evaluate_implicit




  subroutine dust_terms(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, 1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)

    double precision :: tmp(ixI^S), vgas(ixI^S, 1:ndir)
    double precision :: alpha(ixI^S, 1:ndir, 1:dust_n_species)
    integer          :: n, idir

    do idir=1,ndir
      vgas(ixO^S,idir)=w(ixO^S,gas_mom(idir))/w(ixO^S,gas_rho_)
    end do
    call get_alpha_dust(ixI^L, ixO^L, w, vgas,x, alpha)
    w(ixO^S, gas_e_)=0d0
    do idir = 1, ndir

      w(ixO^S, gas_mom(idir))=0d0
      do n = 1, dust_n_species
        ! contribution for gas momentum
        tmp(ixO^S) = alpha(ixO^S, idir,n) * (  &
            w(ixO^S,dust_rho(n)) * w(ixO^S, gas_mom(idir)) - &
            w(ixO^S,gas_rho_) * w(ixO^S, dust_mom(idir, n)))
        w(ixO^S, dust_mom(idir, n)) = -tmp(ixO^S)
        if (dust_backreaction) then
          w(ixO^S, gas_mom(idir)) = w(ixO^S, gas_mom(idir)) + tmp(ixO^S)
          if (gas_e_ > 0) then
            if(dust_backreaction_fh) then
              where(w(ixO^S,dust_rho(n)) > 0d0)
                w(ixO^S, gas_e_) = w(ixO^S, gas_e_) + alpha(ixO^S, idir,n) * &
                (w(ixO^S,gas_rho_) * (w(ixO^S, dust_mom(idir,n))**2/w(ixO^S,dust_rho(n))) - &
                w(ixO^S,dust_rho(n)) * (w(ixO^S, gas_mom(idir))**2/w(ixO^S,gas_rho_)))
              elsewhere
                w(ixO^S, gas_e_) = w(ixO^S, gas_e_) + alpha(ixO^S, idir,n) * ( - &
                w(ixO^S,dust_rho(n)) * (w(ixO^S, gas_mom(idir))**2/w(ixO^S,gas_rho_)))
              endwhere
            else  
              w(ixO^S, gas_e_) = w(ixO^S, gas_e_) + vgas(ixO^S, idir)  &
                 * tmp(ixO^S)
            end if
          end if
        end if
      end do
    end do
  end subroutine dust_terms

    !> Implicit solve of psb=psa+dtfactor*dt*F_im(psb)
  subroutine dust_implicit_update(dtfactor,qdt,qtC,psb,psa)
    use mod_global_parameters
    !use mod_ghostcells_update

    type(state), target :: psa(max_blocks)
    type(state), target :: psb(max_blocks)
    double precision, intent(in) :: qdt
    double precision, intent(in) :: qtC
    double precision, intent(in) :: dtfactor

    integer :: iigrid, igrid

    !call getbc(global_time,0.d0,psa,1,nw)
    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
      block=>psa(igrid)
      call dust_advance_implicit_grid(ixG^LL, ixG^LL, psa(igrid)%w, psb(igrid)%w, psa(igrid)%x, dtfactor,qdt)
    end do
    !$OMP END PARALLEL DO

   end subroutine dust_implicit_update 

  subroutine dust_advance_implicit_grid(ixI^L, ixO^L, w, wout, x, dtfactor,qdt)
    use mod_global_parameters
    integer, intent(in)                :: ixI^L, ixO^L
    double precision, intent(in) :: qdt
    double precision, intent(in) :: dtfactor
    double precision, intent(in)       :: w(ixI^S,1:nw)
    double precision, intent(in)    ::  x(ixI^S,1:ndim)
    double precision, intent(out)       :: wout(ixI^S,1:nw)

    integer                            :: n, m, idir
    double precision :: alpha(ixI^S, 1:ndir, 1:dust_n_species)
    double precision :: tmp(ixI^S),tmp2(ixI^S)
    double precision :: tmp3(ixI^S)
    double precision    :: vgas(ixI^S, 1:ndir)


    do idir = 1, ndir
      vgas(ixO^S,idir)=w(ixO^S,gas_mom(idir))/w(ixO^S,gas_rho_)
    end do
    call get_alpha_dust(ixI^L, ixO^L, w, vgas, x, alpha)
    !TODO this is still neeed?
    wout(ixO^S,1:nw) = w(ixO^S,1:nw)

    do idir = 1, ndir
      ! d1 from Eq 6
      tmp2(ixO^S) = 0d0
      do n = 1, dust_n_species
        tmp2(ixO^S) = tmp2(ixO^S) +  alpha(ixO^S, idir,n) * &
          (w(ixO^S,gas_rho_) + w(ixO^S,dust_rho(n))) 

      enddo
      ! store D in tmp
      tmp(ixO^S) = 1d0 + tmp2(ixO^S) * qdt 
      if(dust_implicit_second_order) then
        ! d2 from Eq 6
        tmp2(ixO^S) = 0d0
        do n = 1, dust_n_species
          do m = n+1, dust_n_species
            tmp2(ixO^S) = tmp3(ixO^S) +  alpha(ixO^S, idir,n) * alpha(ixO^S, idir,m) *&
                 (w(ixO^S,gas_rho_) + w(ixO^S,dust_rho(n))+w(ixO^S,dust_rho(m)))
          enddo
        enddo
        ! multiplied at the end by rho_gas 
        tmp(ixO^S) = tmp(ixO^S) + w(ixO^S,gas_rho_)*tmp2(ixO^S) * (qdt**2)
      endif



      do n = 1, dust_n_species
        ! ni1 from eq 7
        tmp2(ixO^S) = alpha(ixO^S, idir,n) * (  &
            w(ixO^S,dust_rho(n)) * w(ixO^S, gas_mom(idir)) - &
            w(ixO^S,gas_rho_) * w(ixO^S, dust_mom(idir, n))) * qdt

        if(dust_implicit_second_order) then 
          ! ni2 from eq 7 
          tmp3(ixO^S) = 0d0
          do m = n+1, dust_n_species
              tmp3(ixO^S) = tmp3(ixO^S) +  alpha(ixO^S, idir,n) * alpha(ixO^S, idir,m) * &
              ( w(ixO^S,dust_rho(n)) * (w(ixO^S, dust_mom(idir, n)) +  w(ixO^S, gas_mom(idir))) - &
                (w(ixO^S,gas_rho_) + w(ixO^S,dust_rho(m))) * w(ixO^S, dust_mom(idir, n)) )  
          enddo
          ! tmp3 multiplied at the end by rho_gas 
          tmp2(ixO^S) = tmp2(ixO^S) + tmp3(ixO^S) * w(ixO^S,gas_rho_)* (qdt**2) 
        endif
        tmp2(ixO^S) = tmp2(ixO^S)/tmp(ixO^S)
        wout(ixO^S, dust_mom(idir,n)) = w(ixO^S, dust_mom(idir,n)) + tmp2(ixO^S)
      enddo

      if (dust_backreaction) then
        tmp2(ixO^S) = 0d0
        !n1 from eq 8 
        do n = 1, dust_n_species
          tmp2(ixO^S) = tmp2(ixO^S) + alpha(ixO^S, idir,n) * &
              (w(ixO^S,gas_rho_) * w(ixO^S, dust_mom(idir,n)) - &
              w(ixO^S,dust_rho(n)) * w(ixO^S, gas_mom(idir)))

        enddo
        tmp2(ixO^S) = qdt *  tmp2(ixO^S) 
        if(dust_implicit_second_order) then 
          !n2 from eq 8 
          tmp3(ixO^S) = 0d0
          do n = 1, dust_n_species
            do m = n+1, dust_n_species
               tmp3(ixO^S) = tmp3(ixO^S) + alpha(ixO^S, idir,n) * alpha(ixO^S, idir,m) * &
                    (w(ixO^S,gas_rho_) * (w(ixO^S, dust_mom(idir, n)) + w(ixO^S, dust_mom(idir, m))) - &
                    (w(ixO^S,dust_rho(n)) + w(ixO^S,dust_rho(m)))* w(ixO^S, gas_mom(idir)))  
            enddo
          enddo
          ! tmp3 multiplied at the end by rho_gas 
          tmp2(ixO^S) = tmp2(ixO^S) + (qdt**2)*tmp3(ixO^S)* w(ixO^S,gas_rho_)
        endif
        ! store in tmp2 contribution to momentum
        ! so that it is used when dust_backreaction_fh = .false.
        tmp2(ixO^S) = tmp2(ixO^S) / tmp(ixO^S)
        wout(ixO^S, gas_mom(idir)) = w(ixO^S, gas_mom(idir)) + tmp2(ixO^S)

        ! kinetic energy update
         if (gas_e_ > 0) then
          if(dust_backreaction_fh) then 
            ! add work done by coll terms + FrictionalHeating
            tmp2(ixO^S) = 0d0
            do n = 1, dust_n_species
              ! 2*dust kinetic energy: dust rho can be 0
              where(w(ixO^S,dust_rho(n)) > 0d0)
                tmp3(ixO^S)= w(ixO^S, dust_mom(idir,n))**2/w(ixO^S,dust_rho(n))
              elsewhere
                tmp3(ixO^S) = 0d0
              endwhere
              tmp2(ixO^S) = tmp2(ixO^S) + alpha(ixO^S, idir,n) * &
                (w(ixO^S,gas_rho_) * tmp3(ixO^S) - &
                w(ixO^S,dust_rho(n)) * (w(ixO^S, gas_mom(idir))**2/w(ixO^S,gas_rho_)))
  
            enddo
            tmp2(ixO^S) = qdt *  tmp2(ixO^S) 
            if(dust_implicit_second_order) then
              tmp3(ixO^S) = 0d0
              do n = 1, dust_n_species
                do m = n+1, dust_n_species
                    tmp3(ixO^S) = tmp3(ixO^S) + alpha(ixO^S, idir,n) * alpha(ixO^S, idir,m) * &
                      (w(ixO^S,gas_rho_) * (w(ixO^S, dust_mom(idir, n))**2/w(ixO^S,dust_rho(n)) + w(ixO^S, dust_mom(idir,m))**2/w(ixO^S,dust_rho(m))) - &
                      (w(ixO^S,dust_rho(n)) + w(ixO^S,dust_rho(m)))* w(ixO^S, gas_mom(idir))**2/w(ixO^S,gas_rho_))  
                enddo
              enddo
              ! tmp3 multiplied at the end by rho_gas 
              tmp2(ixO^S) = tmp2(ixO^S) + (qdt**2)*tmp3(ixO^S)* w(ixO^S,gas_rho_)
            endif
            wout(ixO^S, gas_e_) = wout(ixO^S, gas_e_) + 0.5d0 * tmp2(ixO^S) / tmp(ixO^S)
          else
            ! dust_backreaction_fh = .false.
            ! add only work done by coll term by multiplyting the contribution in mom eq. by velocity
            wout(ixO^S, gas_e_) = wout(ixO^S, gas_e_) + vgas(ixO^S, idir) * tmp2(ixO^S)
          endif
         end if
      end if
    end do !1..ndir
        

  end subroutine dust_advance_implicit_grid

  ! copied from  get_3d_dragforce subroutine
  subroutine get_alpha_dust(ixI^L, ixO^L, w, vgas,x, alpha)
    use mod_global_parameters
    use mod_usr_methods, only: usr_get_3d_dragforce
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    double precision, intent(in)    :: w(ixI^S, 1:nw)
    double precision,intent(in)    :: vgas(ixI^S, 1:ndir)
    double precision, intent(out)   :: &
         alpha(ixI^S, 1:ndir, 1:dust_n_species)

    double precision    :: ptherm(ixI^S)
    double precision, dimension(ixI^S) :: vt2, deltav, fd, vdust
    double precision                   :: alpha_T(ixI^S, 1:dust_n_species)
    integer                            :: n, idir

    call phys_get_pthermal(w, x, ixI^L, ixO^L, ptherm)

    vt2(ixO^S) = gas_vtherm_factor*ptherm(ixO^S)/w(ixO^S, gas_rho_)

    select case( TRIM(dust_method) )
    case ('Kwok') ! assume sticking coefficient equals 0.25

      do idir = 1, ndir
        do n = 1, dust_n_species
          where(w(ixO^S, dust_rho(n)) > 0.0d0)

            ! 0.75 from sticking coefficient
            fd(ixO^S)     = 0.75d0 / (dust_density(n) * dust_size(n))

            ! 0.75 from spherical grainvolume
            vdust(ixO^S)  = w(ixO^S, dust_mom(idir, n)) / w(ixO^S, dust_rho(n))
            deltav(ixO^S) = vgas(ixO^S, idir)-vdust(ixO^S)
            fd(ixO^S)     = fd(ixO^S)*0.75d0*dsqrt(vt2(ixO^S) + deltav(ixO^S)**2)
          elsewhere
            fd(ixO^S) = 0.0d0
          end where
          alpha(ixO^S, idir, n) = fd(ixO^S)
        end do
      end do

    case ('sticking') ! Calculate sticking coefficient based on the gas and dust temperatures

      call get_sticking(w, x, ixI^L, ixO^L, alpha_T, ptherm)

      do idir = 1, ndir
        do n = 1, dust_n_species
          where(w(ixO^S, dust_rho(n))>0.0d0)
            ! sticking
            fd(ixO^S)     = (one-alpha_T(ixO^S,n)) / (dust_density(n)*dust_size(n))
            ! 0.75 from spherical grainvolume
            vdust(ixO^S)  = w(ixO^S,dust_mom(idir, n)) / w(ixO^S, dust_rho(n))
            deltav(ixO^S) = vgas(ixO^S, idir)-vdust(ixO^S)
            fd(ixO^S)     = fd(ixO^S)*0.75d0*dsqrt(vt2(ixO^S) + deltav(ixO^S)**2)
          else where
            fd(ixO^S) = 0.0d0
          end where
          alpha(ixO^S, idir,n) = fd(ixO^S)
        end do
      end do

    case('linear') !linear with Deltav, for testing (see Laibe & Price 2011)
      do idir = 1, ndir
        do n = 1, dust_n_species
          where(w(ixO^S, dust_rho(n))>0.0d0)
            fd(ixO^S)     = dust_K_lineardrag/(w(ixO^S,gas_rho_)*w(ixO^S, dust_rho(n)))
          else where
            fd(ixO^S) = 0.0d0
          end where
          alpha(ixO^S, idir,n) = fd(ixO^S)
        end do
      end do

    case('none')
      alpha(ixO^S, :, :) = 0.0d0
    case default
      call mpistop( "=== This dust method has not been implemented===" )
    end select

  end subroutine get_alpha_dust


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

    if(dust_dtpar .le. 0) return

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
