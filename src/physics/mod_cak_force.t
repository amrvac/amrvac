!> Module to include CAK radiation line force in (magneto)hydrodynamic models
!> Computes both the force from free electrons and the force from an ensemble of
!> lines (various possibilities for the latter).
!> There is an option to only simulate the pure radial CAK force (with various
!> corrections applied) as well as the full vector CAK force. Depending on the
!> chosen option additional output are the CAK line force component(s) and,
!> when doing a 1-D radial force, the finite disc factor.
!>
!> USAGE:
!>
!>  1. Include a cak_list in the .par file and activate (m)hd_cak_force in the
!>     (m)hd_list
!>  2. Create a mod_usr.t file for the problem with appropriate initial and
!>     boundary conditions
!>  3. In the mod_usr.t header call the mod_cak_force module to have access to
!>     global variables from mod_cak_force, which may be handy for printing or
!>     the computation of other variables inside mod_usr.t
!>  4. In usr_init of mod_usr.t call the set_cak_force_norm routine and pass
!>     along the stellar radius and wind temperature---this is needed to
!>     correctly compute the (initial) force normalisation inside mod_cak_force
!>  5. Always ensure that the set_cak_force_norm routine is called in mod_usr.t
!      either in usr_init or in usr_set_parameters and after the hydro unit
!      variables have been set/computed
!>
!> Developed by Florian Driessen (2022, 2026)
module mod_cak_force

  use mod_physics,     only: phys_get_pthermal, physics_type
  use mod_cak_opacity, only: init_cak_table, set_cak_opacity

  implicit none
  private

  !> Line-ensemble parameters in the Gayley (1995) formalism
  real(8), public :: cak_alpha, gayley_qbar, gayley_q0

  !> Free-electron scattering opacity
  real(8), public :: kappae_cgs = 0.0d0

  !> Ray positions + weights for impact parameter and azimuthal radiation angle
  real(8), allocatable :: ay(:), wy(:), aphi(:), wphi(:)

  !> The adiabatic index
  real(8) :: cak_gamma

  !> Variables needed to compute force normalisation fnorm in initialisation
  real(8) :: lstar_cgs, lstar, rstar, kappae, clight

  !> To enforce a floor temperature (~wind temperature) for adiabatic (M)HD
  real(8) :: tfloor

  ! Telectron/Teff = 0.8 from Puls+ (2000), A&AS 141; assume Twind = Telectron
  real(8), parameter :: ratio_twind_teff = 0.8d0

  !> Method and option for CAK force
  integer :: method_cakforce, type_cak_1d

  !> Type CAK line force method
  integer, parameter :: radialforce=0
  integer, parameter :: vectorforce=1

  !> Type 1-D CAK line force options
  integer, parameter :: pointstar=0
  integer, parameter :: fdisc=1
  integer, parameter :: fdisc_cutoff=2

  !> Amount of rays in radiation polar and radiation azimuthal direction
  integer :: nthetaray, nphiray

  !> Extra slots to store quantities in w-array
  integer, public :: gcak1_, gcak2_, gcak3_, fdf_
  integer, public :: alpha_, qbar_, q0_, kappae_

  !> To treat source term in split or unsplit (default) fashion
  logical :: cak_split=.false.

  !> To activate the pure radial vector CAK line force computation
  logical :: fix_vector_force_1d=.false.

  !> To compute CAK line-force parameters from opacity tables
  logical :: use_cak_table=.false.

  !> Allow reading different opacity table than default src/tables/CAK_tables
  !> If true, name_cak_table requires absolute/relative path to file location
  logical :: use_custom_cak_table=.false.

  !> String to specify the CAK force method: 'radial' or 'vector'
  character(len=6) :: cak_force_method

  !> String to choose between the 1-D CAK line force options
  !> Can be 'pointstar', 'finitedisc', or 'finitedisc_cutoff'
  character(len=256) :: cak_1d_type

  !> String of the opacity table to use from src/tables/CAK_tables
  character(len=256) :: name_cak_table=''

  !> Public methods for mod_hd_phys or mod_mhd_phys
  public :: cak_init
  public :: cak_add_source
  public :: cak_get_dt

  !> Public method for mod_usr
  public :: set_cak_force_norm
  
contains

  !> Read this module's parameters from a file
  subroutine cak_params_read(files)
    use mod_global_parameters, only: unitpar, unitterm
    use mod_comm_lib, only: mpistop

    character(len=*), intent(in) :: files(:)

    ! Local variable
    integer :: n

    namelist /cak_list/ cak_alpha, gayley_qbar, gayley_q0, kappae_cgs, &
         cak_force_method, cak_1d_type, cak_split, &
         nphiray, nthetaray, fix_vector_force_1d, &
         use_cak_table, name_cak_table, use_custom_cak_table

    do n = 1,size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, cak_list, end=111)
       111 close(unitpar)
    enddo

    select case(trim(cak_force_method))
    case('radial')
      method_cakforce = radialforce
    case('vector')
      method_cakforce = vectorforce
    case default
      write(unitterm,*) 'cak_force_method = ', trim(cak_force_method)
      call mpistop('cak_params_read: unknown CAK force in cak_list')
    end select

    ! Set CAK method for pure radial force
    select case(trim(cak_1d_type))
    case('pointstar')
      type_cak_1d = pointstar
    case('finitedisc')
      type_cak_1d = fdisc
    case('finitedisc_cutoff')
      type_cak_1d = fdisc_cutoff
    case default
      write(unitterm,*) 'cak_1d_type = ', trim(cak_1d_type)
      call mpistop('cak_params_read: unknown CAK wind method in cak_list')
    end select

  end subroutine cak_params_read

  !> Initialize the module
  subroutine cak_init(phys_gamma)
    use mod_global_parameters
    use mod_comm_lib, only: mpistop

    real(8), intent(in) :: phys_gamma

    cak_gamma = phys_gamma

    ! Set some defaults when user does not
    cak_alpha   = 0.65d0
    gayley_qbar = 2000.0d0
    gayley_q0   = 2000.0d0

    cak_force_method = 'radial'
    cak_1d_type = 'finitedisc'

    nthetaray = 6
    nphiray   = 6

    call cak_params_read(par_files)

    if (method_cakforce == radialforce) then
      gcak1_ = var_set_extravar("gcak1", "gcak1")
      fdf_   = var_set_extravar("fdfac", "fdfac")
    endif

    if (method_cakforce == vectorforce) then
      gcak1_ = var_set_extravar("gcak1", "gcak1")
      gcak2_ = var_set_extravar("gcak2", "gcak2")
      gcak3_ = var_set_extravar("gcak3", "gcak3")
      call rays_init(nthetaray,nphiray)
    endif

    if (cak_split) any_source_split = .true.

    if (use_cak_table) then
      call init_cak_table(trim(name_cak_table), use_custom_cak_table)

      select case(trim(name_cak_table))
      ! Electron opacity ~= 0.2*(1+X) [cgs] for a fully-ionised plasma with X
      ! the hydrogen mass fraction
      case('Y02400')
        kappae_cgs = 0.34d0
      case('Y09800')
        kappae_cgs = 0.2d0
      end select

      alpha_  = var_set_extravar("alpha", "alpha")
      qbar_   = var_set_extravar("Qbar", "Qbar")
      q0_     = var_set_extravar("Q0", "Q0")
      kappae_ = var_set_extravar("kappae", "kappae")
    else
      alpha_  = var_set_wextra()
      qbar_   = var_set_wextra()
      q0_     = var_set_wextra()
      kappae_ = var_set_wextra()
    endif

    ! Some sanity checks
    if (SI_unit .and. use_cak_table) then
      call mpistop('cak_init: SI_unit=T but LTE tables assume cgs units')
    endif

    if (slab) then
      call mpistop('cak_init: Cartesian geometry not supported')
    endif

    if ((cak_alpha < 0.0d0) .or. (cak_alpha >= 1.0d0)) then
      call mpistop('cak_init: input alpha in [0,1[')
    endif

    if ((gayley_qbar < 0.0d0) .or. (gayley_q0 < 0.0d0)) then
      call mpistop('cak_init: input Qbar or Q0 is < 0')
    endif

    if (method_cakforce == vectorforce .and. ndir < 3) then
      call mpistop('cak_init: vector CAK force only for 2.5D and 3D')
    endif

    if (kappae_cgs < smalldouble) then
      call mpistop('cak_init: set input kappae to a reasonable constant')
    endif

  end subroutine cak_init

  !> Compute some (unitless) variables for CAK force normalisation
  subroutine set_cak_force_norm(rstar_cgs,twind_cgs)
    use mod_global_parameters
    use mod_constants

    real(8), intent(in) :: rstar_cgs, twind_cgs

    lstar_cgs = 4.0d0*dpi * rstar_cgs**2.0d0 * sigma_SB_cgs * twind_cgs**4.0d0

    ! Dimensionless quantities used in this module computations
    kappae = kappae_cgs * unit_density * unit_length
    clight = const_c/unit_velocity
    lstar  = lstar_cgs/(unit_density * unit_length**5.0d0 / unit_time**3.0d0)
    rstar  = rstar_cgs/unit_length
    tfloor = ratio_twind_teff * twind_cgs/unit_temperature

  end subroutine set_cak_force_norm
  
  !> w[iw]=w[iw]+qdt*S[wCT,qtC,x] where S is the source based on wCT within ixO
  subroutine cak_add_source(qdt,ixI^L,ixO^L,wCT,w,x,energy,qsourcesplit,active)
    use mod_global_parameters
    use mod_comm_lib, only: mpistop

    integer, intent(in)    :: ixI^L, ixO^L
    real(8), intent(in)    :: qdt, x(ixI^S,1:ndim), wCT(ixI^S,1:nw)
    real(8), intent(inout) :: w(ixI^S,1:nw)
    logical, intent(in)    :: energy, qsourcesplit
    logical, intent(inout) :: active

    ! Local variables
    integer :: idir, ix^D
    real(8) :: gcak(ixO^S,1:3), ge(ixO^S), ptherm(ixI^S), pmin(ixI^S)
    real(8) :: local_rho_cgs, local_twind_cgs
    real(8) :: table_alpha, table_qbar, table_q0, table_kappae_cgs

    ! By default add source in unsplit fashion together with the fluxes
    if (qsourcesplit .eqv. cak_split) then

      active = .true.

      if (use_cak_table) then
        ! Set line-statistic parameters from local density and wind temperature
        {do ix^DB = ixO^LIM^DB\}
          local_rho_cgs = wCT(ix^D,iw_rho)*unit_density
          local_twind_cgs = tfloor*unit_temperature

          call set_cak_opacity(local_rho_cgs,local_twind_cgs, &
               table_alpha,table_qbar,table_q0,table_kappae_cgs)

          w(ix^D,alpha_)  = table_alpha
          w(ix^D,qbar_)   = table_qbar
          w(ix^D,q0_)     = table_q0
          w(ix^D,kappae_) = table_kappae_cgs
        {enddo^D&\}

        ! Convert electron opacity from table to unitless
        w(ixO^S,kappae_) = w(ixO^S,kappae_) * unit_density * unit_length

      else
        ! Constant line-statistic parameters
        block%wextra(ixO^S,alpha_)  = cak_alpha
        block%wextra(ixO^S,qbar_)   = gayley_qbar
        block%wextra(ixO^S,q0_)     = gayley_q0
        block%wextra(ixO^S,kappae_) = kappae
      endif

      ! Thomson force
      call get_gelectron(ixI^L,ixO^L,w,x,ge)

      ! CAK line force
      select case(method_cakforce)
      case(radialforce)
        call get_cak_force_radial(ixI^L,ixO^L,wCT,w,x,gcak)
      case(vectorforce)
        call get_cak_force_vector(ixI^L,ixO^L,wCT,w,x,gcak)
      end select

      ! Update conservative vars: w = w + qdt*gsource
      do idir = 1,ndir
        if (idir == 1) gcak(ixO^S,idir) = gcak(ixO^S,idir) + ge(ixO^S)
        
        w(ixO^S,iw_mom(idir)) = w(ixO^S,iw_mom(idir)) &
             + qdt * gcak(ixO^S,idir) * wCT(ixO^S,iw_rho)
      enddo

      if (energy) then
        w(ixO^S,iw_e) = w(ixO^S,iw_e) &
             + qdt * sum(gcak(ixO^S,1:ndir) * wCT(ixO^S,iw_mom(1:ndir)))

        ! Impose fixed floor temperature to mimic stellar heating
        call phys_get_pthermal(w,x,ixI^L,ixO^L,ptherm)
        pmin(ixO^S) = w(ixO^S,iw_rho) * tfloor

        where (ptherm(ixO^S) < pmin(ixO^S))
          w(ixO^S,iw_e) = w(ixO^S,iw_e) &
               + (pmin(ixO^S) - ptherm(ixO^S)) / (cak_gamma - 1.0d0)
        endwhere
      endif
    endif

  end subroutine cak_add_source

  !> 1-D CAK line force in the Gayley line-ensemble distribution parametrisation
  subroutine get_cak_force_radial(ixI^L,ixO^L,wCT,w,x,gcak)
    use mod_global_parameters
    use mod_comm_lib, only: mpistop

    integer, intent(in)    :: ixI^L, ixO^L
    real(8), intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    real(8), intent(inout) :: w(ixI^S,1:nw)
    real(8), intent(out)   :: gcak(ixO^S,1:3)
  
    ! Local variables
    real(8) :: vr(ixI^S), dvrdr(ixO^S)
    real(8) :: beta_fd(ixO^S), fdfac(ixO^S), taus(ixO^S), ge(ixO^S)
    real(8) :: alpha(ixO^S), qbar(ixO^S), q0(ixO^S), kappae(ixO^S)

    if (use_cak_table) then
      alpha(ixO^S)  = w(ixO^S,alpha_)
      qbar(ixO^S)   = w(ixO^S,qbar_)
      q0(ixO^S)     = w(ixO^S,q0_)
      kappae(ixO^S) = w(ixO^S,kappae_)
    else
      alpha(ixO^S)  = block%wextra(ixO^S,alpha_)
      qbar(ixO^S)   = block%wextra(ixO^S,qbar_)
      q0(ixO^S)     = block%wextra(ixO^S,q0_)
      kappae(ixO^S) = block%wextra(ixO^S,kappae_)
    endif

    vr(ixI^S) = wCT(ixI^S,iw_mom(1)) / wCT(ixI^S,iw_rho)
    call get_velocity_gradient(ixI^L,ixO^L,vr,x,1,dvrdr)

    if (physics_type == 'hd') then
      ! Monotonic flow to avoid multiple resonances and radiative coupling
      dvrdr(ixO^S) = max(abs(dvrdr(ixO^S)), smalldouble)
    elseif (physics_type == 'mhd') then
      ! Allow material to fallback to the star in a magnetosphere model
      dvrdr(ixO^S) = max(dvrdr(ixO^S), smalldouble)
    endif
  
    ! Thomson force
    call get_gelectron(ixI^L,ixO^L,w,x,ge)

    ! Sobolev optical depth for line ensemble (tau = Qbar * t_r) and the force
    select case (type_cak_1d)
    case(pointstar, fdisc)
      taus(ixO^S) = qbar(ixO^S) * kappae(ixO^S) * clight &
           * wCT(ixO^S,iw_rho) / dvrdr(ixO^S)
      gcak(ixO^S,1) = qbar(ixO^S) / (1.0d0 - alpha(ixO^S)) * ge(ixO^S) &
           / taus(ixO^S)**alpha(ixO^S)

    case(fdisc_cutoff)
      taus(ixO^S) = q0(ixO^S) * kappae(ixO^S) * clight &
           * wCT(ixO^S,iw_rho) / dvrdr(ixO^S)
      gcak(ixO^S,1) = qbar(ixO^S) * ge(ixO^S) / (1.0d0 - alpha(ixO^S)) &
           * ( (1.0d0 + taus(ixO^S))**(1.0d0 - alpha(ixO^S)) - 1.0d0 ) &
           / taus(ixO^S)
    end select

    ! Finite disk factor parameterisation (Owocki & Puls 1996)
    beta_fd(ixO^S) = ( 1.0d0 - vr(ixO^S)/(x(ixO^S,1) * dvrdr(ixO^S)) ) &
         * (rstar/x(ixO^S,1))**2.0d0

    select case (type_cak_1d)
    case(pointstar)
      fdfac(ixO^S) = 1.0d0
    case(fdisc, fdisc_cutoff)
      where (beta_fd(ixO^S) >= 1.0d0)
        fdfac(ixO^S) = 1.0d0/(1.0d0 + alpha(ixO^S))
      elsewhere (beta_fd(ixO^S) < -1.0d10)
        fdfac(ixO^S) = abs(beta_fd(ixO^S))**alpha(ixO^S) &
             / (1.0d0 + alpha(ixO^S))
      elsewhere (abs(beta_fd(ixO^S)) > 1.0d-3)
        fdfac(ixO^S) = &
             (1.0d0 - (1.0d0 - beta_fd(ixO^S))**(1.0d0 + alpha(ixO^S))) &
             / (beta_fd(ixO^S)*(1.0d0 + alpha(ixO^S)))
      elsewhere
        fdfac(ixO^S) = 1.0d0 - 0.5d0*alpha(ixO^S)*beta_fd(ixO^S) &
             * ( 1.0d0 + 1.0d0/3.0d0 * (1.0d0 - alpha(ixO^S))*beta_fd(ixO^S) )
      endwhere
    end select

    ! Correct radial line force for finite disc (if applicable)
    gcak(ixO^S,1) = gcak(ixO^S,1) * fdfac(ixO^S)
    gcak(ixO^S,2) = 0.0d0
    gcak(ixO^S,3) = 0.0d0
      
    ! Fill the nwextra slots for output
    w(ixO^S,gcak1_) = gcak(ixO^S,1)
    w(ixO^S,fdf_)   = fdfac(ixO^S)
    
  end subroutine get_cak_force_radial

  !> Vector CAK line force in the Gayley line-ensemble distribution parametrisation
  subroutine get_cak_force_vector(ixI^L,ixO^L,wCT,w,x,gcak)
    use mod_global_parameters
    use mod_usr_methods

    ! Subroutine arguments
    integer, intent(in)    :: ixI^L, ixO^L
    real(8), intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    real(8), intent(inout) :: w(ixI^S,1:nw)
    real(8), intent(out)   :: gcak(ixO^S,1:3)

    ! Local variables
    integer :: ix^D, itray, ipray
    real(8) :: a1, a2, a3, wyray, y, wpray, phiray, wtot, mustar, dvndn
    real(8) :: costp, sintp, cospp, sinpp, cott0
    real(8) :: vr(ixI^S), vt(ixI^S), vp(ixI^S), inv_rho(ixI^S), inv_r(ixI^S)
    real(8) :: vrr(ixI^S), vtr(ixI^S), vpr(ixI^S)
    real(8) :: dvrdr(ixO^S), dvtdr(ixO^S), dvpdr(ixO^S)
    real(8) :: dvrdt(ixO^S), dvtdt(ixO^S), dvpdt(ixO^S)
    real(8) :: dvrdp(ixO^S), dvtdp(ixO^S), dvpdp(ixO^S)
    real(8) :: gcaktmp1, gcaktmp2, gcaktmp3, taus, integrand
    real(8) :: alpha(ixO^S), qbar(ixO^S), q0(ixO^S), kappae(ixO^S)

    if (use_cak_table) then
      alpha(ixO^S)  = w(ixO^S,alpha_)
      qbar(ixO^S)   = w(ixO^S,qbar_)
      q0(ixO^S)     = w(ixO^S,q0_)
      kappae(ixO^S) = w(ixO^S,kappae_)
    else
      alpha(ixO^S)  = block%wextra(ixO^S,alpha_)
      qbar(ixO^S)   = block%wextra(ixO^S,qbar_)
      q0(ixO^S)     = block%wextra(ixO^S,q0_)
      kappae(ixO^S) = block%wextra(ixO^S,kappae_)
    endif

    inv_rho(ixI^S) = 1.0d0/wCT(ixI^S,iw_rho)
    inv_r(ixI^S)   = 1.0d0/x(ixI^S,1)

    ! Initialisation to have full velocity strain tensor expression at all times
    vt(ixO^S) = 0.0d0; vtr(ixO^S) = 0.0d0
    vp(ixO^S) = 0.0d0; vpr(ixO^S) = 0.0d0
    cott0 = 0.0d0
    dvrdr(ixO^S) = 0.0d0; dvtdr(ixO^S) = 0.0d0; dvpdr(ixO^S) = 0.0d0
    dvrdt(ixO^S) = 0.0d0; dvtdt(ixO^S) = 0.0d0; dvpdt(ixO^S) = 0.0d0
    dvrdp(ixO^S) = 0.0d0; dvtdp(ixO^S) = 0.0d0; dvpdp(ixO^S) = 0.0d0

    ! Populate velocity field(s) depending on dimensions and directions
    vr(ixI^S)  = wCT(ixI^S,iw_mom(1)) * inv_rho(ixI^S)
    vrr(ixI^S) = vr(ixI^S) * inv_r(ixI^S)

    {^NOONED
    vt(ixI^S)  = wCT(ixI^S,iw_mom(2)) * inv_rho(ixI^S)
    vtr(ixI^S) = vt(ixI^S) * inv_r(ixI^S)
    
    if (ndir > 2) then
      vp(ixI^S)  = wCT(ixI^S,iw_mom(3)) * inv_rho(ixI^S)
      vpr(ixI^S) = vp(ixI^S) * inv_r(ixI^S)
    endif
    }
    
    ! Derivatives of velocity field in each coordinate direction (r=1,t=2,p=3)
    call get_velocity_gradient(ixI^L,ixO^L,vr,x,1,dvrdr)
    
    {^NOONED
    call get_velocity_gradient(ixI^L,ixO^L,vr,x,2,dvrdt)
    call get_velocity_gradient(ixI^L,ixO^L,vt,x,1,dvtdr)
    call get_velocity_gradient(ixI^L,ixO^L,vt,x,2,dvtdt)

    if (ndir > 2) then
      call get_velocity_gradient(ixI^L,ixO^L,vp,x,1,dvpdr)
      call get_velocity_gradient(ixI^L,ixO^L,vp,x,2,dvpdt)
    endif
    }
    {^IFTHREED
    call get_velocity_gradient(ixI^L,ixO^L,vr,x,3,dvrdp)
    call get_velocity_gradient(ixI^L,ixO^L,vt,x,3,dvtdp)
    call get_velocity_gradient(ixI^L,ixO^L,vp,x,3,dvpdp)
    }

    ! Get total acceleration from all rays at a certain grid point
    {do ix^DB=ixOmin^DB,ixOmax^DB\}
      gcaktmp1 = 0.0d0
      gcaktmp2 = 0.0d0
      gcaktmp3 = 0.0d0

      ! Loop over the rays; first theta then phi radiation angle
      ! Get weights from current ray and their position
      do itray = 1,nthetaray
        wyray = wy(itray)
        y = ay(itray)

        do ipray = 1,nphiray
          wpray = wphi(ipray)
          phiray = aphi(ipray)

          ! Redistribute the phi rays by a small offset
          ! if (mod(itp,3) == 1) then
          !   phip = phip + dphi/3.0d0
          ! elseif (mod(itp,3) == 2) then
          !   phip = phip - dphi/3.0d0
          ! endif
          
          ! === Geometrical factors ===
          ! Make y quadrature linear in mu, not mu**2; better for gtheta,gphi
          ! y -> mu quadrature is preserved; y=0 <=> mu=1; y=1 <=> mu=mustar
          mustar = sqrt(max(1.0d0 - (rstar*inv_r(ix^D))**2.0d0, 0.0d0))
          costp  = 1.0d0 - y*(1.0d0 - mustar)
          sintp  = sqrt(max(1.0d0 - costp*costp, 0.0d0))
          sinpp  = sin(phiray)
          cospp  = cos(phiray)
          {^NOONED cott0  = cos(x(ix^D,2))/max(sin(x(ix^D,2)), smalldouble)}

          ! More weight close to star, less farther away
          wtot  = wyray * wpray * (1.0d0 - mustar)

          ! Convenients a la Cranmer & Owocki (1995), ApJ 440, eq. 42
          a1 = costp
          a2 = sintp * cospp
          a3 = sintp * sinpp

          ! Get total velocity gradient for one ray with given (theta', phi')
          dvndn = a1*a1 * dvrdr(ix^D) + a2*a2 * (dvtdt(ix^D) + vrr(ix^D))  &
                 + a3*a3 * (dvpdp(ix^D) + cott0 * vtr(ix^D) + vrr(ix^D))   &
                 + a1*a2 * (dvtdr(ix^D) + dvrdt(ix^D) - vtr(ix^D))         &
                 + a1*a3 * (dvpdr(ix^D) + dvrdp(ix^D) - vpr(ix^D))         &
                 + a2*a3 * (dvpdt(ix^D) + dvtdp(ix^D) - cott0 * vpr(ix^D))

          ! No multiple resonances in CAK
          dvndn = abs(dvndn)

          taus = q0(ix^D) * kappae(ix^D) * clight * wCT(ix^D,iw_rho) / dvndn
          integrand = ((1.0d0 + taus)**(1.0d0 - alpha(ix^D)) - 1.0d0) / taus

          ! Convert gradient back from wind coordinates (r',theta',phi') to
          ! stellar coordinates (r,theta,phi)
          gcaktmp1 = gcaktmp1 + integrand * a1 * wtot
          gcaktmp2 = gcaktmp2 + integrand * a2 * wtot
          gcaktmp3 = gcaktmp3 + integrand * a3 * wtot
        enddo
      enddo

      gcak(ix^D,1:3) = [gcaktmp1, gcaktmp2, gcaktmp3] &
           * kappae(ix^D) * qbar(ix^D) / (1.0d0 - alpha(ix^D))
    {enddo\}

    ! Normalisation for line force array
    ! NOTE: extra 1/pi factor comes from integration in radiation Phi angle
    gcak = gcak/dpi * lstar/(4.0d0*dpi*rstar**2.0d0 * clight)

    if (fix_vector_force_1d) then
      gcak(ixO^S,2) = 0.0d0
      gcak(ixO^S,3) = 0.0d0
    endif
          
    ! Fill the nwextra slots for output
    w(ixO^S,gcak1_) = gcak(ixO^S,1)
    w(ixO^S,gcak2_) = gcak(ixO^S,2)
    w(ixO^S,gcak3_) = gcak(ixO^S,3)

  end subroutine get_cak_force_vector
  
  !> Compute continuum radiation force from Thomson scattering
  subroutine get_gelectron(ixI^L,ixO^L,w,x,ge)
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    real(8), intent(in) :: w(ixI^S,1:nw), x(ixI^S,1:ndim)
    real(8), intent(out):: ge(ixO^S)

    real(8) :: kappae(ixO^S)

    if (use_cak_table) then
      kappae(ixO^S) = w(ixO^S,kappae_)
    else
      kappae(ixO^S) = block%wextra(ixO^S,kappae_)
    endif

    ge(ixO^S) = kappae(ixO^S) * lstar/(4.0d0*dpi * clight * x(ixO^S,1)**2.0d0)

  end subroutine get_gelectron

  !> Check time step for total radiation contribution
  subroutine cak_get_dt(wprim,ixI^L,ixO^L,dtnew,dx^D,x)
    use mod_global_parameters

    integer, intent(in)    :: ixI^L, ixO^L
    real(8), intent(in)    :: dx^D, x(ixI^S,1:ndim)
    real(8), intent(in)    :: wprim(ixI^S,1:nw)
    real(8), intent(inout) :: dtnew
    
    ! Local variables
    real(8) :: ge(ixO^S), max_gr, dt_cak

    call get_gelectron(ixI^L,ixO^L,wprim,x,ge)

    dtnew = bigdouble

    ! Get dt from line force that is saved in the w-array in nwextra slot
    max_gr = max( maxval(abs(ge(ixO^S) + wprim(ixO^S,gcak1_))), epsilon(1.0d0) )
    dt_cak = minval( sqrt(block%dx(ixO^S,1)/max_gr) )
    dtnew  = min(dtnew, courantpar*dt_cak)

    {^NOONED
    if (method_cakforce == vectorforce) then
      max_gr = max( maxval(abs(wprim(ixO^S,gcak2_))), epsilon(1.0d0) )
      dt_cak = minval( sqrt(block%dx(ixO^S,1) * block%dx(ixO^S,2)/max_gr) )
      dtnew  = min(dtnew, courantpar*dt_cak)

      {^IFTHREED
      max_gr = max( maxval(abs(wprim(ixO^S,gcak3_))), epsilon(1.0d0) )
      dt_cak = minval( sqrt(block%dx(ixO^S,1) * sin(block%dx(ixO^S,3))/max_gr) )
      dtnew  = min(dtnew, courantpar*dt_cak)
      }
    endif
    }

  end subroutine cak_get_dt

  !> Compute velocity gradient in direction 'idir' on a non-uniform grid
  subroutine get_velocity_gradient(ixI^L,ixO^L,v,x,idir,grad_vn)
    use mod_global_parameters

    integer, intent(in)  :: ixI^L, ixO^L, idir
    real(8), intent(in)  :: v(ixI^S), x(ixI^S,1:ndim)
    real(8), intent(out) :: grad_vn(ixO^S)

    ! Local variables
    real(8) :: forw(ixO^S), backw(ixO^S), cent(ixO^S)
    integer :: jrx^L, hrx^L{^NOONED,jtx^L, htx^L}{^IFTHREED,jpx^L, hpx^L}

    ! Index +1 (j) and index -1 (h) in radial direction; kr(dir,dim)=1, dir=dim
    jrx^L=ixO^L+kr(1,^D);
    hrx^L=ixO^L-kr(1,^D);

    {^NOONED
    ! Index +1 (j) and index -1 (h) in polar direction
    jtx^L=ixO^L+kr(2,^D);
    htx^L=ixO^L-kr(2,^D);
    }

    {^IFTHREED
    ! Index +1 (j) and index -1 (h) in azimuthal direction
    jpx^L=ixO^L+kr(3,^D);
    hpx^L=ixO^L-kr(3,^D);
    }

    ! grad(v.n) on non-uniform grid according to Sundqvist & Veronis (1970)
    select case (idir)
    case(1) ! Radial forward, backward, and central derivatives
      forw(ixO^S) = (x(ixO^S,1) - x(hrx^S,1)) * v(jrx^S) &
           / ((x(jrx^S,1) - x(ixO^S,1)) * (x(jrx^S,1) - x(hrx^S,1)))

      backw(ixO^S) = -(x(jrx^S,1) - x(ixO^S,1)) * v(hrx^S) &
           / ((x(ixO^S,1) - x(hrx^S,1)) * (x(jrx^S,1) - x(hrx^S,1)))

      cent(ixO^S) = (x(jrx^S,1) + x(hrx^S,1) - 2.0d0*x(ixO^S,1)) * v(ixO^S) &
           / ((x(ixO^S,1) - x(hrx^S,1)) * (x(jrx^S,1) - x(ixO^S,1)))
    {^NOONED
    case(2) ! Polar forward, backward, and central derivatives
      forw(ixO^S) = (x(ixO^S,2) - x(htx^S,2)) * v(jtx^S) &
           / (x(ixO^S,1) * (x(jtx^S,2) - x(ixO^S,2)) * (x(jtx^S,2) - x(htx^S,2)))

      backw(ixO^S) = -(x(jtx^S,2) - x(ixO^S,2)) * v(htx^S) &
           / ( x(ixO^S,1) * (x(ixO^S,2) - x(htx^S,2)) * (x(jtx^S,2) - x(htx^S,2)))

      cent(ixO^S) = (x(jtx^S,2) + x(htx^S,2) - 2.0d0*x(ixO^S,2)) * v(ixO^S) &
           / ( x(ixO^S,1) * (x(ixO^S,2) - x(htx^S,2)) * (x(jtx^S,2) - x(ixO^S,2)))
    }
    {^IFTHREED
    case(3) ! Azimuthal forward, backward, and central derivatives
      forw(ixO^S) = (x(ixO^S,3) - x(hpx^S,3)) * v(jpx^S) &
           / ( x(ixO^S,1)*sin(x(ixO^S,2)) * (x(jpx^S,3) - x(ixO^S,3)) * (x(jpx^S,3) - x(hpx^S,3)))

      backw(ixO^S) = -(x(jpx^S,3) - x(ixO^S,3)) * v(hpx^S) &
           / ( x(ixO^S,1)*sin(x(ixO^S,2)) * (x(ixO^S,3) - x(hpx^S,3)) * (x(jpx^S,3) - x(hpx^S,3)))

      cent(ixO^S) = (x(jpx^S,3) + x(hpx^S,3) - 2.0d0*x(ixO^S,3)) * v(ixO^S) &
           / ( x(ixO^S,1)*sin(x(ixO^S,2)) * (x(ixO^S,3) - x(hpx^S,3)) * (x(jpx^S,3) - x(ixO^S,3)))
    }
    end select

    ! Total gradient for given velocity field
    grad_vn(ixO^S) = backw(ixO^S) + cent(ixO^S) + forw(ixO^S)

  end subroutine get_velocity_gradient

  !> Initialise (theta',phi') radiation angles coming from stellar disc
  subroutine rays_init(ntheta_point,nphi_point)
    use mod_global_parameters

    ! Subroutine arguments
    integer, intent(in) :: ntheta_point, nphi_point

    ! Local variables
    real(8) :: ymin, ymax, phipmin, phipmax, adum
    integer :: ii

    ! Minimum and maximum range of theta and phi rays
    ! NOTE: theta points are cast into y-space
    ymin    = 0.0d0
    ymax    = 1.0d0
    phipmin = -dpi !0.0d0
    phipmax = dpi !2.0d0*dpi
    ! dphi    = (phipmax - phipmin) / nphi_point

    if (mype == 0) then
      allocate(ay(ntheta_point))
      allocate(wy(ntheta_point))
      allocate(aphi(nphi_point))
      allocate(wphi(nphi_point))

      ! theta and phi ray positions and weights: Gauss-Legendre
      call gauss_legendre_quadrature(ymin,ymax,ntheta_point,ay,wy)
      call gauss_legendre_quadrature(phipmin,phipmax,nphi_point,aphi,wphi)

      ! theta rays and weights: uniform
      ! dth = 1.0d0 / nthetap
      ! adum = ymin + 0.5d0*dth
      ! do ip = 1,nthetap
      !   ay(ip) = adum
      !   wy(ip) = 1.0d0/nthetap
      !   adum = adum + dth
      !   !print*,'phipoints'
      !   !print*,ip,aphi(ip),wphi(ip),dphi
      ! enddo

      ! phi ray position and weights: uniform
      ! adum = phipmin + 0.5d0*dphi
      ! do ii = 1,nphi_point
      !   aphi(ii) = adum
      !   wphi(ii) = 1.0d0/nphi_point
      !   adum     = adum + dphi
      ! enddo

      print*, '==========================='
      print*, '    Radiation ray setup    '
      print*, '==========================='
      print*, 'Theta ray points + weights '
      do ii = 1,ntheta_point
        print*,ii,ay(ii),wy(ii)
      enddo
      print*
      print*, 'Phi ray points + weights   '
      do ii = 1,nphi_point
        print*,ii,aphi(ii),wphi(ii)
      enddo
      print*
    endif

    call MPI_BARRIER(icomm,ierrmpi)

    !===========================
    ! Broadcast what mype=0 read
    !===========================
    if (npe > 1) then
      if (mype /= 0) then
        allocate(ay(ntheta_point))
        allocate(wy(ntheta_point))
        allocate(aphi(nphi_point))
        allocate(wphi(nphi_point))
      endif

      call MPI_BCAST(ay,ntheta_point,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
      call MPI_BCAST(wy,ntheta_point,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
      call MPI_BCAST(aphi,nphi_point,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
      call MPI_BCAST(wphi,nphi_point,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
    endif

  end subroutine rays_init
  
  !> Fast Gauss-Legendre N-point quadrature algorithm by G. Rybicki
  subroutine gauss_legendre_quadrature(xlow,xhi,n,x,w)
    ! Given the lower and upper limits of integration xlow and xhi, and given n,
    ! this routine returns arrays x and w of length n, containing the abscissas
    ! and weights of the Gauss-Legendre N-point quadrature
    use mod_global_parameters

    ! Subroutine arguments
    real(8), intent(in)  :: xlow, xhi
    integer, intent(in)  :: n
    real(8), intent(out) :: x(n), w(n)

    ! Local variables
    real(8) :: p1, p2, p3, pp, xl, xm, z, z1
    real(8), parameter :: error=3.0d-14
    integer :: i, j, m

    m = (n + 1)/2
    xm = 0.5d0*(xhi + xlow)
    xl = 0.5d0*(xhi - xlow)

    do i = 1,m
      z = cos( dpi * (i - 0.25d0)/(n + 0.5d0) )
      z1 = 2.0d0 * z

      do while (abs(z1 - z) > error)
        p1 = 1.0d0
        p2 = 0.0d0

        do j = 1,n
          p3 = p2
          p2 = p1
          p1 = ( (2.0d0*j - 1.0d0)*z*p2 - (j - 1.0d0)*p3 )/j
        enddo

        pp = n*(z*p1 - p2) / (z*z - 1.0d0)
        z1 = z
        z = z1 - p1/pp
      enddo

      x(i)     = xm - xl*z
      x(n+1-i) = xm + xl*z
      w(i)     = 2.0d0*xl / ((1.0d0 - z*z) * pp*pp)
      w(n+1-i) = w(i)
    enddo

  end subroutine gauss_legendre_quadrature

end module mod_cak_force
