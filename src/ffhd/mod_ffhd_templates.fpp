#:mute
#:include "../mod_gpu_directives.fpp"
#:endmute
#:if PHYS == 'ffhd'

#:def phys_vars()

  integer, parameter :: dp = kind(0.0d0)
  !> Base number of physical variables (density, momentum, energy, (q), b1, b2, b3)
  integer, parameter, public              :: nw_phys_base=6
  !> Total number of physical variables (base + conditional additions)
  !> This will be set to nw_phys_base + 1 if HYPERTC is defined
  integer, parameter, public              :: nw_phys=nw_phys_base &
#:if defined('HYPERTC')
    + 1
#:else
    + 0
#:endif

  integer, parameter, public              :: nw_flux_base=3
  integer, parameter, public              :: nw_flux=nw_flux_base &
#:if defined('HYPERTC')
    + 1
#:else
    + 0
#:endif

  !> Number of leading variables the MUSCL reconstruction limits (1:nw_recon).
  !> Slots nw_recon+1:nw_phys are left to fix_prim_state, which must rebuild
  !> them before the Riemann solver reads the face state.
  !> ffhd reconstructs everything: get_flux reads the frozen field b1..b3 at
  !> the faces, and it sits behind the flux variables.
  integer, parameter, public              :: nw_recon=nw_phys

  !> Whether an energy equation is used
  logical, public                         :: ffhd_energy = .true.
  ${GPU_DECLARE_COPYIN('ffhd_energy')}$

  !> Index of the density (in the w array)
  integer, public                         :: rho_
  ${GPU_DECLARE_CREATE('rho_')}$

  !> Indices of the momentum density
  integer, allocatable, public            :: mom(:)
  ${GPU_DECLARE_CREATE('mom')}$

  !> Index of the energy density (-1 if not present)
  integer, public                         :: e_
  ${GPU_DECLARE_CREATE('e_')}$

  !> Index of the gas pressure (-1 if not present) should equal e_
  integer, public                         :: p_
  ${GPU_DECLARE_CREATE('p_')}$

  !> Index of the hyperbolic flux variable
  integer, public                         :: q_
  ${GPU_DECLARE_CREATE('q_')}$

  ! !> Index of the frozen magnetic field
  integer, public                         :: iw_b1
  integer, public                         :: iw_b2
  integer, public                         :: iw_b3
  ${GPU_DECLARE_CREATE('iw_b1, iw_b2, iw_b3')}$

  !> The adiabatic index
  double precision, public                :: ffhd_gamma = 5.d0/3.0d0
  double precision, public                :: gamma_1, inv_gamma_1
  ${GPU_DECLARE_COPYIN('ffhd_gamma, gamma_1, inv_gamma_1')}$
  ${GPU_DECLARE_CREATE('gamma_1, inv_gamma_1')}$

  !> The adiabatic constant
  double precision, public                :: ffhd_adiab = 1.0d0
  ${GPU_DECLARE_COPYIN('ffhd_adiab')}$

  !> The helium abundance
  double precision, public                :: He_abundance=0.1d0
  ${GPU_DECLARE_COPYIN('He_abundance')}$

  !> The thermal conductivity kappa in hyperbolic thermal conduction
  double precision, public                :: hypertc_kappa=-1.0d0
  ${GPU_DECLARE_COPYIN('hypertc_kappa')}$

  !> Whether p*divb source term is not zero
  logical, public                         :: ffhd_pdivb = .false.
  ${GPU_DECLARE_COPYIN('ffhd_pdivb')}$

  !> switch for hyperbolic thermal conduction
  logical, public                         :: ffhd_hyperbolic_thermal_conduction = .false.
  ${GPU_DECLARE_COPYIN('ffhd_hyperbolic_thermal_conduction')}$

  !> Whether plasma is partially ionized
  logical, public                         :: ffhd_partial_ionization = .false.
  ${GPU_DECLARE_COPYIN('ffhd_partial_ionization')}$

  !> switch for gravity
  logical, public                         :: ffhd_gravity = .false.
  ${GPU_DECLARE_COPYIN('ffhd_gravity')}$

  !> switch for radiative cooling
  logical, public                         :: ffhd_radiative_cooling = .false.
  ${GPU_DECLARE_COPYIN('ffhd_radiative_cooling')}$

  !> switch for source user
  logical, public                         :: ffhd_source_usr = .false.
  ${GPU_DECLARE_COPYIN('ffhd_source_usr')}$

  !> Smallest gas pressure allowed in a reconstructed face state, applied by
  !> fix_prim_state. Zero, the default, clips only negative values.
  double precision, public                :: ffhd_small_pressure = 0.0d0
  !> Smallest density allowed in a reconstructed face state
  double precision, public                :: ffhd_small_density  = 0.0d0
  ${GPU_DECLARE_COPYIN('ffhd_small_pressure,ffhd_small_density')}$

#:enddef

#:def read_params()
    !> Read this module's parameters from a file
  subroutine read_params(files)
    use mod_global_parameters
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /ffhd_list/ ffhd_energy, ffhd_gamma, ffhd_partial_ionization, ffhd_gravity, &
          ffhd_radiative_cooling, ffhd_hyperbolic_thermal_conduction, ffhd_source_usr, ffhd_pdivb, He_abundance, &
      ffhd_small_pressure, ffhd_small_density

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, ffhd_list, end=111)
111    close(unitpar)
    end do

    if (ffhd_small_pressure < 0.0d0) call mpistop(&
       "ffhd_small_pressure should be positive.")
    if (ffhd_small_density < 0.0d0) call mpistop(&
       "ffhd_small_density should be positive.")

    ${GPU_UPDATE_DEVICE('ffhd_energy, ffhd_gamma, ffhd_partial_ionization, ffhd_gravity, ffhd_radiative_cooling, ffhd_hyperbolic_thermal_conduction, ffhd_source_usr, ffhd_pdivb, He_abundance')}$
    ${GPU_UPDATE_DEVICE('ffhd_small_pressure, ffhd_small_density')}$


  end subroutine read_params
#:enddef

#:def phys_activate() 
  subroutine phys_activate()
    call phys_init()
  end subroutine phys_activate
#:enddef

#:def phys_units()
  subroutine phys_units()
    use mod_global_parameters
    double precision :: mp, kB
    double precision :: a,b

    !> here no SI_UNIT used by default, to be implemented
    mp = mp_cgs
    kB = kB_cgs
    !> eq_state_units by default, to be implemented
    a = 1.d0+4.d0*He_abundance
    b = 2.d0+3.d0*He_abundance

    if(unit_density/=1.d0 .or. unit_numberdensity/=1.d0) then
      if(unit_density/=1.d0) then
        unit_numberdensity=unit_density/(a*mp)
      else if(unit_numberdensity/=1.d0) then
        unit_density=a*mp*unit_numberdensity
      end if
      if(unit_temperature/=1.d0) then
        unit_pressure=b*unit_numberdensity*kB*unit_temperature
        unit_velocity=dsqrt(unit_pressure/unit_density)
        if(unit_length/=1.d0) then
          unit_time=unit_length/unit_velocity
        else if(unit_time/=1.d0) then
          unit_length=unit_velocity*unit_time
        end if
      else if(unit_pressure/=1.d0) then
        unit_temperature=unit_pressure/(b*unit_numberdensity*kB)
        unit_velocity=dsqrt(unit_pressure/unit_density)
        if(unit_length/=1.d0) then
          unit_time=unit_length/unit_velocity
        else if(unit_time/=1.d0) then
          unit_length=unit_velocity*unit_time
        end if
      else if(unit_velocity/=1.d0) then
        unit_pressure=unit_density*unit_velocity**2
        unit_temperature=unit_pressure/(b*unit_numberdensity*kB)
        if(unit_length/=1.d0) then
          unit_time=unit_length/unit_velocity
        else if(unit_time/=1.d0) then
          unit_length=unit_velocity*unit_time
        end if
      else if(unit_time/=1.d0) then
        unit_velocity=unit_length/unit_time
        unit_pressure=unit_density*unit_velocity**2
        unit_temperature=unit_pressure/(b*unit_numberdensity*kB)
      end if
    else if(unit_temperature/=1.d0) then
      ! units of temperature and velocity are dependent
      if(unit_pressure/=1.d0) then
        unit_numberdensity=unit_pressure/(b*unit_temperature*kB)
        unit_density=a*mp*unit_numberdensity
        unit_velocity=dsqrt(unit_pressure/unit_density)
        if(unit_length/=1.d0) then
          unit_time=unit_length/unit_velocity
        else if(unit_time/=1.d0) then
          unit_length=unit_velocity*unit_time
        end if
      end if
    else if(unit_pressure/=1.d0) then
      if(unit_velocity/=1.d0) then
        unit_density=unit_pressure/unit_velocity**2
        unit_numberdensity=unit_density/(a*mp)
        unit_temperature=unit_pressure/(b*unit_numberdensity*kB)
        if(unit_length/=1.d0) then
          unit_time=unit_length/unit_velocity
        else if(unit_time/=1.d0) then
          unit_length=unit_velocity*unit_time
        end if
      else if(unit_time/=0.d0) then
        unit_velocity=unit_length/unit_time
        unit_density=unit_pressure/unit_velocity**2
        unit_numberdensity=unit_density/(a*mp)
        unit_temperature=unit_pressure/(b*unit_numberdensity*kB)
      end if
    end if
    unit_mass=unit_density*unit_length**3

    ${GPU_UPDATE_DEVICE('unit_density, unit_numberdensity, unit_temperature, unit_pressure, unit_velocity, unit_length, unit_time, unit_mass')}$
  end subroutine phys_units
#:enddef
  
#:def phys_init()
    !> Initialize the module
  subroutine phys_init()

    use mod_global_parameters
    #:if defined('COOLING')
    use mod_radiative_cooling, only: rc_fl, radiative_cooling_init_params, radiative_cooling_init
    #:endif

    call read_params(par_files)
    call phys_units()

    phys_energy  = ffhd_energy
    phys_total_energy  = ffhd_energy
    phys_internal_e = .false.
    phys_gamma = ffhd_gamma
    phys_partial_ionization=ffhd_partial_ionization
    gamma_1 = ffhd_gamma - 1.0d0
    inv_gamma_1 = 1.0d0/gamma_1
    ${GPU_UPDATE_DEVICE('physics_type, phys_energy, phys_total_energy, phys_internal_e, phys_gamma, phys_partial_ionization, gamma_1, inv_gamma_1')}$

    ! Determine flux variables
    rho_ = var_set_rho()
    ${GPU_UPDATE_DEVICE('rho_')}$

    allocate(mom(1))
    mom(:) = var_set_momentum(1)
    ${GPU_UPDATE_DEVICE('mom')}$

    ! Set index of energy variable
    if (ffhd_energy) then
       e_ = var_set_energy()
       p_ = e_
    else
       e_ = -1
       p_ = -1
    end if
    ${GPU_UPDATE_DEVICE('e_,p_')}$

    ! Set index for heat flux
#:if defined('HYPERTC')
    q_ = var_set_q()
    need_global_cmax = .true.
    hypertc_kappa = 8.d-7*unit_temperature**3.5_dp/unit_length/unit_density/unit_velocity**3.0_dp
    ${GPU_UPDATE_DEVICE('q_')}$
    ${GPU_UPDATE_DEVICE('need_global_cmax')}$
    ${GPU_UPDATE_DEVICE('hypertc_kappa')}$
#:endif

    ! Set index of frozen magnetic field
    iw_b1 = var_set_extravar('b1', 'b1')
    iw_b2 = var_set_extravar('b2', 'b2')
    iw_b3 = var_set_extravar('b3', 'b3')
    ${GPU_UPDATE_DEVICE('iw_b1, iw_b2, iw_b3')}$

    ! set number of variables which need update ghostcells
    ! The frozen field (b1,b2,b3) is registered via var_set_extravar (it is
    ! not advected), but it is NOT exchanged through getbc: it is a pure
    ! function of position, so with the FILL_NWEXTRA_ANALYTIC flag (which
    ! config_schema.toml's `implies` sets for phys='ffhd') fill_nwextra_device
    ! re-derives it from usr_set_nwextra in every cell of every block -
    ! interior, inter-block ghosts, physical-boundary ghosts and polar-axis
    ! ghosts alike - after each tree change. Keeping it out of getbc is what
    ! lets the polar axis work: evaluating the user's analytic field at a ghost
    ! cell's own (mirrored) coordinates reproduces the pole sign flips for
    ! free, with no entry in typeboundary. So nwgc stops at nwflux, exactly as
    ! for hd and mhd.
    nwgc=nwflux
    ${GPU_UPDATE_DEVICE('nwgc')}$

    ! Define custom flux types:
    if (.not. allocated(flux_type)) then
       allocate(flux_type(ndir, nw_flux))
       flux_type = flux_default
    else if (any(shape(flux_type) /= [ndir, nw_flux])) then
       call mpistop("phys_check error: flux_type has wrong shape")
    end if
    ${GPU_UPDATE_DEVICE('flux_type')}$

    ! Only the frozen field (b1,b2,b3) is a genuine vector here: mom(1) is a
    ! scalar (the field-aligned momentum m_par), not 3 consecutive slots, so
    ! it must not be registered as one.
    nvector      = 1 ! No. vector vars
    allocate(iw_vector(nvector))
    iw_vector(1) = iw_b1 - 1
    ${GPU_UPDATE_DEVICE('nvector, iw_vector')}$

#:if defined('COOLING')
    call radiative_cooling_init_params(phys_gamma,He_abundance)
    call radiative_cooling_init(rc_fl)
    ${GPU_UPDATE_DEVICE('rc_fl')}$
    ${GPU_ENTER_DATA_COPYIN('rc_fl%tcool,rc_fl%Lcool, rc_fl%Yc')}$
#:endif

  end subroutine phys_init
#:enddef

#:def phys_get_dt()
  subroutine phys_get_dt(w, x, dx, dtnew)
#:if defined('GRAVITY')
  use mod_usr, only: gravity_field
#:endif    
  ${GPU_ROUTINE_SEQ()}$
    real(dp), intent(in)   :: w(nw_phys), x(1:ndim), dx(1:ndim)
    real(dp), intent(out)  :: dtnew
    ! .. local ..
    integer                :: idim
    real(dp)               :: field

    dtnew = huge(1.0d0)
    
#:if defined('GRAVITY')
    do idim = 1, ndim
       field = gravity_field(w, x, idim)
       field = max( abs(field), epsilon(1.0d0) )
       dtnew = min( dtnew, 1_dp / sqrt( field/dx(idim) ) )
    end do
#:endif    
    
  end subroutine phys_get_dt
#:enddef  

#:def addsource_local()
subroutine addsource_local(qdt, dtfactor, qtC, wCT, wCTprim, qt, wnew, x, dr, &
    qsourcesplit)
#:if defined('SOURCE_USR')
  use mod_usr, only: addsource_usr
#:endif
#:if defined('GRAVITY')
  use mod_usr, only: gravity_field
#:endif    
#:if defined('COOLING')
  use mod_radiative_cooling, only: radiative_cooling_add_source
#:endif
  ${GPU_ROUTINE_SEQ()}$

  real(dp), intent(in)     :: qdt, dtfactor, qtC, qt
  real(dp), intent(in)     :: wCT(nw_phys), wCTprim(nw_phys)
  real(dp), intent(in)     :: x(1:ndim), dr(ndim)
  real(dp), intent(inout)  :: wnew(nw_phys)
  logical, intent(in)      :: qsourcesplit
  ! .. local ..
  integer                  :: idim
  real(dp)                 :: field, mag, divb

  if (.not. qsourcesplit) then 
     !---------------------------------
     ! unsplit sources
     !---------------------------------
     
#:if defined('GRAVITY')
     do idim = 1, ndim
        field = gravity_field(wCT, x, idim)
        mag = wCTprim(iw_b1-1+idim)
        wnew(iw_mom(1)) = wnew(iw_mom(1)) + qdt * field * wCT(iw_rho) * mag
        wnew(iw_e)      = wnew(iw_e) + qdt * field * wCT(iw_mom(1)) * mag
     end do
#:endif  

#:if defined('COOLING')
     call radiative_cooling_add_source(qdt,wCT,wCTprim,wnew,x)
#:endif

#:if defined('SOURCE_USR')
     call addsource_usr(qdt, qt, wCT, wCTprim, wnew, x, .false.)
#:endif

  else
     !---------------------------------
     ! split sources     
     !---------------------------------

     ! Not yet implemented

  end if
  
end subroutine addsource_local
#:enddef

#:def addsource_nonlocal()
subroutine addsource_nonlocal(qdt, dtfactor, qtC, wCTprim, qt, wnew, x, dx, idir, &
     qsourcesplit)
  use mod_global_parameters, only: dt, cmax_global, courantpar, third
  ${GPU_ROUTINE_SEQ()}$

  real(dp), intent(in)     :: qdt, dtfactor, qtC, qt
  real(dp), intent(in)     :: wCTprim(nw_phys,5)
  real(dp), intent(in)     :: x(1:ndim), dx(1:ndim)
  real(dp), intent(inout)  :: wnew(nw_phys)
  integer, intent(in)      :: idir
  logical, intent(in)      :: qsourcesplit
  ! .. local ..
  real(dp)                 :: tau, htc_qrsc, sigT
  real(dp)                 :: Te(1:5), gradT
  real(dp)                 :: mag5(1:5), divb, mag

  if (.not. qsourcesplit) then 
     !---------------------------------
     ! unsplit sources
     !---------------------------------

#:if defined('PDIVB')
     ! > p*divb 
     mag5(1:5) = wCTprim(iw_b1-1+idir,1:5)
     divb = (8*mag5(4)-8*mag5(2)-mag5(5)+mag5(1))/12.0_dp/dx(idir)
     wnew(iw_mom(1)) = wnew(iw_mom(1)) + qdt*wCTprim(iw_e,3)*divb
#:endif

#:if defined('HYPERTC')
     !> gradient of temperature:
     Te(1:5) = wCTprim(iw_e,1:5) / wCTprim(iw_rho,1:5)
     mag = wCTprim(iw_b1-1+idir,3)

     gradT = (8.d0*(Te(4)-Te(2))-Te(5)+Te(1))/(12.d0*dx(idir))

     sigT = hypertc_kappa * sqrt(Te(3)**5)
     tau = max(4.d0*dt, sigT*Te(3)*courantpar**2*gamma_1/&
        (wCTprim(iw_e,3)*cmax_global**2))

     htc_qrsc = sigT * mag * gradT

     wnew(iw_q) = wnew(iw_q) - qdt * (htc_qrsc + wCTprim(iw_q,3)*third) / tau
#:endif

  else
     !---------------------------------
     ! split sources     
     !---------------------------------

     ! Not yet implemented

  end if
  
end subroutine addsource_nonlocal
#:enddef

#:if GEOM == 'spherical' or GEOM == 'cylindrical'
!> No curvature source terms are needed for ffhd. Its conserved quantities
!> (rho, m_par, e_hd_par, q_par) are all genuine scalars advected along the
!> user-supplied field direction b-hat, fluxed as (quantity) * b-hat; by
!> Gauss's theorem the divergence of a scalar's flux needs no extra terms in
!> any coordinate system (unlike a vector quantity's components, which pick
!> up curvature terms from the position-dependent (r, theta, phi) or
!> (r, z, phi) basis - see HD's/SRHD's/MHD's addsource_geometry). Matches
!> upstream MPI-AMRVAC's ffhd_add_source_geom, which is likewise empty for
!> every coordinate system.
!>
!> Caveat: this does NOT cover the optional PDIVB source term (p * div(b-hat),
!> compile-time ffhd_pdivb=T), whose div(b-hat) in addsource_nonlocal is a
!> plain finite difference over dx(idir) - correct on a slab-uniform mesh, but
!> missing the metric factors a true curvilinear divergence needs (r^2/sin(theta)
!> for spherical, r for cylindrical). PDIVB defaults off and is untouched here.
#:def addsource_geometry()
subroutine addsource_geometry(qdt, wprim, wnew, x, dAdV)
  ${GPU_ROUTINE_SEQ()}$

  real(dp), intent(in)     :: qdt
  real(dp), intent(in)     :: wprim(nw_phys)
  real(dp), intent(in)     :: x(1:ndim)
  real(dp), intent(in)     :: dAdV(1:ndim)
  real(dp), intent(inout)  :: wnew(nw_phys)

end subroutine addsource_geometry
#:enddef
#:endif

#:def to_primitive()
pure subroutine to_primitive(u)
  ${GPU_ROUTINE_SEQ()}$
  real(dp), intent(inout) :: u(nw_phys)

  u(iw_mom(1)) = u(iw_mom(1))/u(iw_rho)

  u(iw_e) = gamma_1 * (u(iw_e) - 0.5_dp * u(iw_rho) * &
     u(iw_mom(1))**2 )

end subroutine to_primitive
#:enddef

#:def to_conservative()  
pure subroutine to_conservative(u)
  ${GPU_ROUTINE_SEQ()}$
  real(dp), intent(inout) :: u(nw_phys)

  ! Compute energy from pressure and kinetic energy
  u(iw_e) = u(iw_e) * inv_gamma_1 + 0.5_dp * u(iw_rho) * &
     u(iw_mom(1))**2

  ! Compute momentum from density and velocity components
  u(iw_mom(1)) = u(iw_rho) * u(iw_mom(1))
  
end subroutine to_conservative
#:enddef

#:def fix_prim_state()
  !> Make a reconstructed primitive face state admissible before it is handed
  !> to the Riemann solver. The MUSCL reconstruction limits every slot of the
  !> state independently, so a face state can come out with a negative
  !> density or pressure even where both neighbouring cells are admissible.
  !> Floor both, at ffhd_small_density and ffhd_small_pressure in &ffhd_list.
  !> Both default to zero, in which case the max() calls only clip values
  !> that were already negative and leave every admissible state bit-for-bit
  !> unchanged. Only the face states are floored, not the cell averages.
  pure subroutine fix_prim_state(u)
    ${GPU_ROUTINE_SEQ()}$
    real(dp), intent(inout) :: u(nw_phys)

    u(iw_rho) = max(u(iw_rho), ffhd_small_density)
    u(iw_e)   = max(u(iw_e),   ffhd_small_pressure)

  end subroutine fix_prim_state
#:enddef

#:def get_flux()
subroutine get_flux(u, xC, flux_dim, flux)
  ${GPU_ROUTINE_SEQ()}$

  real(dp), intent(in)  :: u(nw_phys)
  real(dp), intent(in)  :: xC(1:ndim)
  integer, intent(in)   :: flux_dim
  real(dp), intent(out) :: flux(nw_flux)
  real(dp)              :: mag

  mag = u(iw_b1-1+flux_dim)

  ! Density flux
  flux(iw_rho) = u(iw_rho) * u(iw_mom(1)) * mag

  ! Momentum flux with pressure term
  flux(iw_mom(1)) = (u(iw_rho)*u(iw_mom(1))**2 + u(iw_e)) * mag
  
  ! Energy flux with hyperbolic conduction included
  flux(iw_e) = u(iw_mom(1))*(u(iw_e)*inv_gamma_1 + &
               0.5_dp*u(iw_rho)*u(iw_mom(1))**2 + u(iw_e)) * mag + &
#:if defined('HYPERTC')
  u(iw_q)*mag

  flux(iw_q) = 0.0_dp
#:else
  0.0_dp
#:endif

end subroutine get_flux
#:enddef

#:def get_cmax()  
pure real(dp) function get_cmax(u, x, flux_dim) result(wC)
  ${GPU_ROUTINE_SEQ()}$
  real(dp), intent(in)  :: u(nw_phys)
  real(dp), intent(in)  :: x(1:ndim)
  integer, intent(in)   :: flux_dim
  real(dp)              :: mag

  mag = u(iw_b1-1+flux_dim)

  wC = dsqrt(phys_gamma*u(iw_e)/u(iw_rho)) + abs(u(iw_mom(1))*mag)

end function get_cmax
#:enddef  

#:def get_rho()
pure real(dp) function get_rho(w, x) result(rho)
  ${GPU_ROUTINE_SEQ()}$
  real(dp), intent(in)  :: w(nw_phys)
  real(dp), intent(in)  :: x(1:ndim)

  rho = w(iw_rho)
end function get_rho
#:enddef

#:def get_pthermal()
pure real(dp) function get_pthermal(w, x) result(pth)
  ${GPU_ROUTINE_SEQ()}$
  real(dp), intent(in)  :: w(nw_phys)
  real(dp), intent(in)  :: x(1:ndim)

  pth = gamma_1*(w(iw_e)-0.5_dp*w(iw_mom(1))**2/w(iw_rho))
end function get_pthermal
#:enddef

#:def get_Rfactor()
pure real(dp) function get_Rfactor() result(Rfactor)
  ${GPU_ROUTINE_SEQ()}$
  Rfactor = 1.0d0
end function get_Rfactor
#:enddef

#:endif
