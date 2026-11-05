#:mute
#:include "../mod_gpu_directives.fpp"
#:endmute
#:if PHYS == 'srhd'

#:if defined('N_TRACER')
#:set N_TRACER_ = N_TRACER
#:else
#:set N_TRACER_ = 0
#:endif

#:def phys_vars()

  integer, parameter :: dp = kind(0.0d0)
  integer, parameter, public              :: nw_phys=2+ndim+2+${N_TRACER_}$
  integer, parameter, public              :: nw_flux=2+ndim+${N_TRACER_}$
  
  !> Number of leading variables the MUSCL reconstruction limits (1:nw_recon).
  !> Slots nw_recon+1:nw_phys are left to fix_prim_state, which must rebuild
  !> them before the Riemann solver reads the face state.
  integer, parameter, public              :: nw_recon=nw_flux

  !> Whether synge eos is used
  logical, public                         :: srhd_eos = .false.
  ${GPU_DECLARE_COPYIN('srhd_eos')}$

  !> Index of the density (in the w array) as primitive or conserved
  integer, public                         :: rho_
  integer, public                         :: d_
  ${GPU_DECLARE_CREATE('rho_,d_')}$

  !> Indices of the momentum density
  integer, allocatable, public            :: mom(:)
  ${GPU_DECLARE_CREATE('mom')}$

#:if defined('N_TRACER')
  !> Indices of the tracers
  integer, public                         :: tracer(${N_TRACER_}$)
  ${GPU_DECLARE_CREATE('tracer')}$
#:endif

  !> Index of the energy density
  integer, public                         :: e_
  ${GPU_DECLARE_CREATE('e_')}$

  !> Index of the gas pressure should equal e_
  integer, public                         :: p_
  ${GPU_DECLARE_CREATE('p_')}$

  !> Index of the Lorentz factor
  integer, public     :: lfac_
  ${GPU_DECLARE_CREATE('lfac_')}$

  !> Index of the inertia
  integer, public     :: xi_
  ${GPU_DECLARE_CREATE('xi_')}$

  !> Number of tracer species
  integer, public                         :: srhd_n_tracer = 0
  ${GPU_DECLARE_COPYIN('srhd_n_tracer')}$

  !> The adiabatic index
  double precision, public                :: srhd_gamma = 5.d0/3.0d0
  ${GPU_DECLARE_COPYIN('srhd_gamma')}$

  !> derived values from adiabatic index 
  double precision, public                :: gamma_1,inv_gamma_1,gamma_to_gamma_1
  ${GPU_DECLARE_COPYIN('gamma_1,inv_gamma_1,gamma_to_gamma_1')}$

  !> Helium abundance over Hydrogen
  double precision, public  :: He_abundance=0.1d0
  ${GPU_DECLARE_COPYIN('He_abundance')}$

  !> Smallest gas pressure allowed. Both the lower bracket of the pressure
  !> root-find in con2prim and the floor fix_prim_state applies to a
  !> reconstructed face state. Zero, the default, leaves both unfloored.
  double precision, public                :: srhd_small_pressure = 0.0d0
  !> Smallest rest-mass density allowed
  double precision, public                :: srhd_small_density  = 0.0d0
  ${GPU_DECLARE_COPYIN('srhd_small_pressure,srhd_small_density')}$

  !> Whether particles module is added
  logical, public                         :: srhd_particles = .false.
  ${GPU_DECLARE_COPYIN('srhd_particles')}$

  !> switch for source user
  logical, public                         :: srhd_source_usr = .false.
  ${GPU_DECLARE_COPYIN('srhd_source_usr')}$

#:enddef

#:def read_params()
    !> Read this module's parameters from a file
  subroutine read_params(files)
    use mod_global_parameters
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /srhd_list/ srhd_eos,srhd_gamma,srhd_n_tracer, &
      He_abundance, srhd_source_usr, &
      srhd_small_pressure, srhd_small_density

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, srhd_list, end=111)
111    close(unitpar)
    end do

    if (srhd_small_pressure < 0.0d0) call mpistop(&
       "srhd_small_pressure should be positive.")
    if (srhd_small_density < 0.0d0) call mpistop(&
       "srhd_small_density should be positive.")

    ${GPU_UPDATE_DEVICE('srhd_eos, srhd_gamma, srhd_n_tracer, He_abundance, srhd_source_usr')}$
    ${GPU_UPDATE_DEVICE('srhd_small_pressure, srhd_small_density')}$

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
    double precision :: mp,kB

    !> here no SI_UNIT used by default, to be implemented
    mp = mp_cgs
    kB = kB_cgs
    unit_velocity=const_c

    ! we assume user sets: unit_numberdensity, unit_length, He_abundance
    ! then together with light speed c, all units fixed
    unit_density=(1.0d0+4.0d0*He_abundance)*mp*unit_numberdensity
    unit_pressure=unit_density*unit_velocity**2
    unit_temperature=unit_pressure/((2.0d0+3.0d0*He_abundance)*unit_numberdensity*kB)
    unit_time=unit_length/unit_velocity
    unit_mass=unit_density*unit_length**3

    ${GPU_UPDATE_DEVICE('unit_density, unit_numberdensity, unit_temperature, unit_pressure, unit_velocity, unit_length, unit_time, unit_mass')}$
  end subroutine phys_units
#:enddef
  
#:def phys_init()
    !> Initialize the module
  subroutine phys_init()
    use mod_global_parameters
    integer      :: idir

    call read_params(par_files)
    call phys_units()

    phys_energy  = .true.
    phys_total_energy  = .true.
    phys_gamma = srhd_gamma

    gamma_1=srhd_gamma-1.0d0
    inv_gamma_1=1.0d0/gamma_1
    gamma_to_gamma_1=srhd_gamma/gamma_1
    ${GPU_UPDATE_DEVICE('gamma_1,inv_gamma_1,gamma_to_gamma_1')}$

    phys_internal_e=.false.
    phys_partial_ionization=.false.
    need_global_cmax=.false.

    ! Whether diagonal ghost cells are required for the physics
    phys_req_diagonal = .false.

    ${GPU_UPDATE_DEVICE('physics_type, phys_energy, phys_total_energy, phys_internal_e, phys_gamma, phys_partial_ionization,need_global_cmax,phys_req_diagonal')}$

    use_particles = srhd_particles

    ! Determine flux variables
    rho_ = var_set_rho()
    d_=rho_
    ${GPU_UPDATE_DEVICE('rho_,d_')}$

    allocate(mom(ndir))
    mom(:) = var_set_momentum(ndir)
    ${GPU_UPDATE_DEVICE('mom')}$

    ! Set index of energy variable
    e_ = var_set_energy()
    p_ = e_
    ${GPU_UPDATE_DEVICE('e_,p_')}$

    ! Register tracer fields
#:if defined('N_TRACER')
    #:for i in range(1, N_TRACER_+1)
        tracer(${i}$) = var_set_fluxvar("trc", "trp", ${i}$, need_bc=.false.)
    #:endfor
    ${GPU_UPDATE_DEVICE('tracer')}$
#:endif

    ! Set index for auxiliary variables
    ! MUST be after the possible tracers (which have fluxes)
    xi_  = var_set_auxvar('xi','xi')
    lfac_= var_set_auxvar('lfac','lfac')

    ! muscl_reconstruct_prim limits only slots 1:nw_recon and leaves xi and
    ! lfac to fix_prim_state, so every flux variable has to come before them
    if (nwflux > nw_recon .or. xi_ <= nw_recon .or. lfac_ <= nw_recon) &
         call mpistop("srhd phys_init: flux variables must occupy 1:nw_recon "//&
         "and the auxiliaries xi, lfac come after them")
    ${GPU_UPDATE_DEVICE('xi_,lfac_')}$

    ! set number of variables which need update ghostcells
    nwgc=nwflux+nwaux
    ${GPU_UPDATE_DEVICE('nwgc')}$

    ! Define custom flux types:
    if (.not. allocated(flux_type)) then
       allocate(flux_type(ndir, nw_flux))
       flux_type = flux_default
    else if (any(shape(flux_type) /= [ndir, nw_flux])) then
       call mpistop("phys_check error: flux_type has wrong shape")
    end if
    ${GPU_UPDATE_DEVICE('flux_type')}$

    nvector      = 1 ! No. vector vars
    allocate(iw_vector(nvector))
    iw_vector(1) = mom(1) - 1
    ${GPU_UPDATE_DEVICE('nvector, iw_vector')}$

! use cycle, needs to be dealt with:
!    ! Initialize particles module
!    if (srhd_particles) then
!       call particles_init()
!       phys_req_diagonal = .true.
!    end if

  end subroutine phys_init
#:enddef

#:def phys_get_dt()
  subroutine phys_get_dt(w, x, dx, dtnew)
  ${GPU_ROUTINE_SEQ()}$
    real(dp), intent(in)   :: w(nw_phys), x(1:ndim), dx(1:ndim)
    real(dp), intent(out)  :: dtnew

    dtnew = huge(1.0d0)
    
#:if defined('SOURCE_USR')
    ! TODO: user-set time step limit, also in other physics modules!!!
#:endif    

  end subroutine phys_get_dt
#:enddef  

#:def addsource_local()
subroutine addsource_local(qdt, dtfactor, qtC, wCT, wCTprim, qt, wnew, x, dr, &
    qsourcesplit)
#:if defined('SOURCE_USR')
  use mod_usr, only: addsource_usr
#:endif
  ${GPU_ROUTINE_SEQ()}$

  real(dp), intent(in)     :: qdt, dtfactor, qtC, qt
  real(dp), intent(in)     :: wCT(nw_phys), wCTprim(nw_phys)
  real(dp), intent(in)     :: x(1:ndim), dr(ndim)
  real(dp), intent(inout)  :: wnew(nw_phys)
  logical, intent(in)      :: qsourcesplit

  if (.not. qsourcesplit) then 
     !---------------------------------
     ! unsplit sources
     !---------------------------------

#:if defined('SOURCE_USR')
     call addsource_usr(qdt, qt, wCT, wCTprim, wnew, x, .false.)
#:endif

  !!!else
     !---------------------------------
     ! split sources     
     !---------------------------------
        
  end if

end subroutine addsource_local
#:enddef

#:if GEOM == 'spherical'
!> Curvature source terms for spherical coordinates, ported from upstream
!> MPI-AMRVAC's srhd_add_source_geom. Structurally identical to this fork's
!> HD addsource_geometry (src/hd/mod_hd_templates.fpp): the relativistic
!> spatial stress tensor is T^ij = S^i v^j + p delta^ij, i.e. HD's rho v^i v^j
!> with the conserved momentum density S^i = xi * v^i in place of rho * v^i.
!> wprim(iw_mom(:)) holds the spatial four-velocity (celerity) u^i = lfac*v^i
!> rather than v^i itself (see to_primitive/to_conservative and
!> mod_con2prim.fpp's xi = tau + D + p, v^2 = S^2/xi^2), so
!> S^i v^j = xi * v^i * v^j = xi * u^i * u^j / lfac^2. Every geometric
!> prefactor is the discrete dAdV = (A_upper - A_lower)/dV rather than the
!> continuous 2/r and cot(theta)/r upstream uses directly, matching HD's
!> convention: for the metrics fill_geometry_device builds,
!> dAdV(1) = 2*<1/r> and dAdV(2) = <cot(theta)/r> exactly, so these terms are
!> well balanced and independent of where in the cell bgeo%x sits - which
!> matters, since bgeo%x is the volume barycentre, not the face midpoint.
#:def addsource_geometry()
subroutine addsource_geometry(qdt, wprim, wnew, x, dAdV)
  ${GPU_ROUTINE_SEQ()}$

  real(dp), intent(in)     :: qdt
  !> primitive variables (rho, spatial four-velocity, pressure, xi, lfac)
  !> at the current stage
  real(dp), intent(in)     :: wprim(nw_phys)
  !> cell-centre coordinates (r, theta, phi); unused - every geometric factor
  !> comes from dAdV, see the note above
  real(dp), intent(in)     :: x(1:ndim)
  !> (upper minus lower face area) / cell volume, per direction
  real(dp), intent(in)     :: dAdV(1:ndim)
  real(dp), intent(inout)  :: wnew(nw_phys)
  ! .. local ..
  real(dp)                 :: pth, xi_inv_lfac2, inv_r, cot_r, source

  pth          = wprim(iw_e)
  xi_inv_lfac2 = wprim(xi_) / wprim(lfac_)**2
  inv_r        = 0.5_dp * dAdV(1)   ! <1/r>
  cot_r        = dAdV(2)            ! <cot(theta)/r>

  ! s[m_r] = (2 p + (xi/lfac^2) (u_theta^2 + u_phi^2)) / r
  source = (2.0_dp * pth + xi_inv_lfac2 * (wprim(iw_mom(2))**2 + &
     wprim(iw_mom(3))**2)) * inv_r
  wnew(iw_mom(1)) = wnew(iw_mom(1)) + qdt * source

  ! s[m_theta] = (p cot(theta) + (xi/lfac^2) (u_phi^2 cot(theta) - u_r u_theta)) / r
  source = (pth + xi_inv_lfac2 * wprim(iw_mom(3))**2) * cot_r - &
     xi_inv_lfac2 * wprim(iw_mom(1)) * wprim(iw_mom(2)) * inv_r
  wnew(iw_mom(2)) = wnew(iw_mom(2)) + qdt * source

  ! s[m_phi] = -(xi/lfac^2) u_phi (u_r + u_theta cot(theta)) / r
  source = -xi_inv_lfac2 * wprim(iw_mom(3)) * (wprim(iw_mom(1)) * inv_r + &
     wprim(iw_mom(2)) * cot_r)
  wnew(iw_mom(3)) = wnew(iw_mom(3)) + qdt * source

end subroutine addsource_geometry
#:enddef
#:elif GEOM == 'cylindrical'
!> Curvature source terms for cylindrical (r, z, phi) coordinates. Structurally
!> identical to this fork's HD addsource_geometry (src/hd/mod_hd_templates.fpp)
!> with the conserved momentum density S^i = xi * v^i in place of rho * v^i,
!> and using wprim(iw_mom(:)) = u^i = lfac*v^i as in the spherical branch above
!> (so S^i v^j = xi * u^i * u^j / lfac^2). Only m_r and m_phi pick up curvature
!> terms; m_z does not. Every geometric prefactor is the discrete dAdV,
!> matching HD's convention; for cylindrical dAdV(1) is exactly <1/r> - and,
!> since a cylindrical radial face area is linear in r, exactly 1/r_midpoint
!> as well. It is *not* 1/x(1): bgeo%x holds the volume barycentre.
#:def addsource_geometry()
subroutine addsource_geometry(qdt, wprim, wnew, x, dAdV)
  ${GPU_ROUTINE_SEQ()}$

  real(dp), intent(in)     :: qdt
  !> primitive variables (rho, spatial four-velocity, pressure, xi, lfac)
  !> at the current stage
  real(dp), intent(in)     :: wprim(nw_phys)
  !> cell-centre coordinates (r, z, phi); unused - every geometric factor
  !> comes from dAdV, see the note above
  real(dp), intent(in)     :: x(1:ndim)
  !> (upper minus lower face area) / cell volume, per direction
  real(dp), intent(in)     :: dAdV(1:ndim)
  real(dp), intent(inout)  :: wnew(nw_phys)
  ! .. local ..
  real(dp)                 :: pth, xi_inv_lfac2, inv_r, source

  pth          = wprim(iw_e)
  xi_inv_lfac2 = wprim(xi_) / wprim(lfac_)**2
  inv_r        = dAdV(1)   ! <1/r>

  ! s[m_r] = (p + (xi/lfac^2) u_phi^2) / r
  source = (pth + xi_inv_lfac2 * wprim(iw_mom(3))**2) * inv_r
  wnew(iw_mom(1)) = wnew(iw_mom(1)) + qdt * source

  ! s[m_phi] = -(xi/lfac^2) u_phi u_r / r
  source = -xi_inv_lfac2 * wprim(iw_mom(3)) * wprim(iw_mom(1)) * inv_r
  wnew(iw_mom(3)) = wnew(iw_mom(3)) + qdt * source

end subroutine addsource_geometry
#:enddef
#:endif

#:def get_aux()
  !> The auxiliary variables lfac and xi of a conserved state.
  !>
  !> This is the single definition of lfac and xi: to_primitive and set_aux
  !> both call it, so the two cannot drift apart. The stored values in u are
  !> the initial guess of the root solve.
  pure subroutine get_aux(u, lfac, xi)
    use mod_con2prim
    ${GPU_ROUTINE_SEQ()}$
    real(dp), intent(in)  :: u(nw_phys)
    real(dp), intent(out) :: lfac, xi

    real(dp) :: rho, tau, ssqr

    ssqr=(u(iw_mom(1))**2+u(iw_mom(2))**2+u(iw_mom(3))**2)
    lfac = u(lfac_)
    xi   = u(xi_)
    rho  = u(iw_rho)
    tau  = u(iw_e)
    if(srhd_eos)then
       call con2prim_eos(lfac,xi,rho,ssqr,tau)
    else
       call con2prim(lfac,xi,rho,ssqr,tau)
    endif

  end subroutine get_aux
#:enddef

#:def set_aux()
  !> Overwrite the auxiliary slots (lfac, xi) of a conserved state with the
  !> values consistent with its other components. Nothing else in u changes.
  pure subroutine set_aux(u)
    ${GPU_ROUTINE_SEQ()}$
    real(dp), intent(inout) :: u(nw_phys)
    real(dp) :: lfac, xi

    call get_aux(u, lfac, xi)

    u(lfac_) = lfac
    u(xi_)   = xi

  end subroutine set_aux
#:enddef

#:def to_primitive()
  pure subroutine to_primitive(u)
    ${GPU_ROUTINE_SEQ()}$
    real(dp), intent(inout) :: u(nw_phys)

    real(dp) :: rho,rhoh,pth,E
    real(dp) :: lfac, xi

    call get_aux(u, lfac, xi)
    u(lfac_) = lfac
    u(xi_)   = xi

    rho=u(iw_rho)/u(lfac_)
    rhoh=u(xi_)/u(lfac_)**2
    !! begin: call srhd_get_pressure_eos(rho,rhoh,pth,E)
    if(srhd_eos)then
      E = (rhoh+dsqrt(rhoh**2+(srhd_gamma**2-1.0_dp)*rho**2)) &
         /(srhd_gamma+1.0_dp)
      pth = 0.5d0*gamma_1* (E-rho*(rho/E))
    else
      pth =(rhoh-rho)/gamma_to_gamma_1
    endif 
    !! end: call srhd_get_pressure_eos(rho,rhoh,pth,E)

    u(iw_rho)=rho
    u(iw_mom(1))=u(lfac_)*u(iw_mom(1))/u(xi_)
    u(iw_mom(2))=u(lfac_)*u(iw_mom(2))/u(xi_)
    u(iw_mom(3))=u(lfac_)*u(iw_mom(3))/u(xi_)
    u(iw_e)=pth

  end subroutine to_primitive
#:enddef

#:def to_conservative()  
  pure subroutine to_conservative(u)
    ${GPU_ROUTINE_SEQ()}$
    real(dp), intent(inout) :: u(nw_phys)

    real(dp) :: rho,rhoh,pth
    real(dp) :: E,E_th

    rhoh=(u(iw_mom(1))**2+u(iw_mom(2))**2+u(iw_mom(3))**2)
    u(lfac_)=dsqrt(1.0d0+rhoh)

    rho=u(iw_rho)
    pth=u(iw_e)
    !! begin: call srhd_get_enthalpy_eos(rho,pth,rhoh)
    if(srhd_eos) then
     E_th = pth*inv_gamma_1
     E    = E_th+dsqrt(E_th**2+rho**2)
     ! writing rho/E on purpose, for numerics 
     rhoh = 0.5_dp*((srhd_gamma+1.0_dp)*E &
                   - gamma_1*rho*(rho/E))
    else
     rhoh = rho+gamma_to_gamma_1*pth
    end if
    !! end: call srhd_get_enthalpy_eos(rho,pth,rhoh)

    rhoh=rhoh*u(lfac_)
    u(xi_)=u(lfac_)*rhoh

    u(iw_rho)=u(lfac_)*rho
    u(iw_mom(1))=rhoh*u(iw_mom(1))
    u(iw_mom(2))=rhoh*u(iw_mom(2))
    u(iw_mom(3))=rhoh*u(iw_mom(3))
    u(iw_e)=u(xi_)-pth-u(iw_rho)

  end subroutine to_conservative
#:enddef

#:def fix_prim_state()
  !> Make a reconstructed primitive face state admissible before it is handed
  !> to the Riemann solver.
  !>
  !> Two things happen here, and only the second is optional. The MUSCL
  !> reconstruction limits every slot of the state independently, including the
  !> auxiliaries xi and lfac -- but those are not independent variables, they
  !> are functions of (rho, u^i, p). A reconstructed face state therefore does
  !> not in general satisfy lfac = sqrt(1+|u|^2) or xi = lfac^2*rho*h, while
  !> get_flux, get_cmax and estimate_speeds_minmax all read them. Left alone,
  !> v2 = 1-1/lfac^2 and v_n = u_n/lfac come from different sources, v_n can
  !> exceed v (even exceed 1), the radicand of the relativistic characteristic
  !> speeds goes negative and the HLL wave speeds collapse onto each other --
  !> which removes the solver's dissipation exactly where it is needed. So the
  !> auxiliaries are always recomputed from the reconstructed primitives, by
  !> the same expressions to_conservative uses.
  !>
  !> The floors are the optional part. Both default to zero, in which case the
  !> max() calls only clamp states that were already unphysical; set
  !> srhd_small_pressure (and srhd_small_density) positive in &srhd_list to get
  !> an actual floor, which is the useful configuration -- a state floored to
  !> exactly zero still gives csound2 = gamma*p/(rho*h) a zero denominator.
  pure subroutine fix_prim_state(u)
    ${GPU_ROUTINE_SEQ()}$
    real(dp), intent(inout) :: u(nw_phys)

    real(dp) :: rho, pth, rhoh, E, E_th

    rho = max(u(iw_rho), srhd_small_density)
    pth = max(u(iw_e),   srhd_small_pressure)
    u(iw_rho) = rho
    u(iw_e)   = pth

    u(lfac_) = dsqrt(1.0d0 + u(iw_mom(1))**2 + u(iw_mom(2))**2 &
         + u(iw_mom(3))**2)

    !! begin: call srhd_get_enthalpy_eos(rho,pth,rhoh)
    if (srhd_eos) then
       E_th = pth*inv_gamma_1
       E    = E_th + dsqrt(E_th**2 + rho**2)
       ! writing rho/E on purpose, for numerics
       rhoh = 0.5_dp*((srhd_gamma + 1.0_dp)*E - gamma_1*rho*(rho/E))
    else
       rhoh = rho + gamma_to_gamma_1*pth
    end if
    !! end: call srhd_get_enthalpy_eos(rho,pth,rhoh)

    u(xi_) = u(lfac_)**2 * rhoh

  end subroutine fix_prim_state
#:enddef

#:def get_flux()
  subroutine get_flux(u, xC, flux_dim, flux)
    use mod_global_parameters, only:cmax_global
    ${GPU_ROUTINE_SEQ()}$
    real(dp), intent(in)  :: u(nw_phys)
    real(dp), intent(in)  :: xC(1:ndim)
    integer, intent(in)   :: flux_dim
    real(dp), intent(out) :: flux(nw_flux)

    real(dp) :: pth,vel(1:ndim)

    pth    = u(iw_e)
    vel(1) = u(iw_mom(1)) / u(lfac_)
    vel(2) = u(iw_mom(2)) / u(lfac_)
    vel(3) = u(iw_mom(3)) / u(lfac_)

    ! Density flux
    flux(iw_rho) = u(iw_rho) * u(lfac_) * vel(flux_dim)

    ! Momentum flux with pressure term
    flux(iw_mom(1))        = u(xi_) * vel(1) * vel(flux_dim)
    flux(iw_mom(2))        = u(xi_) * vel(2) * vel(flux_dim)
    flux(iw_mom(3))        = u(xi_) * vel(3) * vel(flux_dim)
    flux(iw_mom(flux_dim)) = flux(iw_mom(flux_dim)) + pth

    ! Energy flux
    flux(iw_e) = vel(flux_dim) * ( u(xi_) - u(iw_rho) * u(lfac_) )

    ! Tracer flux. Note that tracers stay conservative.
#:if defined('N_TRACER')
  #:for i in range(1, N_TRACER_+1)
      flux(tracer(${i}$)) = u(tracer(${i}$)) * vel(flux_dim)
  #:endfor
#:endif

  end subroutine get_flux
#:enddef

#:def get_cmax()
!> Returns maximum local signal speed from primitive state u in direction flux_dim;
!> used in LLF/TVDLF flux estimation.
pure real(dp) function get_cmax(u, x, flux_dim) result(wC)
  ${GPU_ROUTINE_SEQ()}$
  real(dp), intent(in)  :: u(nw_phys)
  real(dp), intent(in)  :: x(1:ndim)
  integer, intent(in)   :: flux_dim

  real(dp) :: tmp1,tmp2,v2,csound2,vidim,cmax,cmin
  real(dp) :: rho,rhoh,pth,E

    ! input u is in primitive form
    ! auxiliaries are filled here
    rho=u(iw_rho)
    rhoh=u(xi_)/u(lfac_)**2
    pth=u(iw_e)
    !!tmp1=u(iw_rho)
    !!tmp2=u(xi_)/u(lfac_)**2.0d0
    !!v2=u(iw_e)
    !! begin: call srhd_get_csound2_prim_eos(tmp1,tmp2,v2,csound2)
    if(srhd_eos) then
       E = (rhoh+dsqrt(rhoh**2+(srhd_gamma**2-1.0_dp)&
            *rho**2))/(srhd_gamma+1.0_dp)
       csound2=(pth*((srhd_gamma+1.0_dp)&
                          +gamma_1*(rho/E)**2))&
                      /(2.0d0*rhoh)
    else
       csound2=srhd_gamma*pth/rhoh
    end if
    !! end call srhd_get_csound2_prim_eos(tmp1,tmp2,v2,csound2)

    v2    = 1.0d0-1.0d0/u(lfac_)**2
    vidim = u(iw_mom(flux_dim))/u(lfac_)
    tmp2  = vidim**2
    tmp1  = 1.0d0-v2*csound2-tmp2*(1.0d0-csound2)
    ! the radicand is positive analytically, as in estimate_speeds_minmax; clip
    ! it so an inadmissible state (c_s^2 > 1 from a collapsed rhoh) yields a
    ! degenerate speed rather than a NaN. get_cmax feeds setdt, where a single
    ! NaN survives the global MPI_ALLREDUCE and takes dt down for every rank.
    tmp2  = dsqrt(max(csound2*(1.0_dp-v2)*tmp1, 0.0_dp))
    tmp1  = vidim*(1.0_dp-csound2)
    cmax  = (tmp1+tmp2)/(1.0_dp-v2*csound2)
    cmin  = (tmp1-tmp2)/(1.0_dp-v2*csound2)
    ! Limit by speed of light
    cmin = max(cmin, - 1.0d0)
    cmin = min(cmin,   1.0d0)
    cmax = max(cmax, - 1.0d0)
    cmax = min(cmax,   1.0d0)
    ! now take extremal value only for dt limit
    wC = max(dabs(cmax),dabs(cmin))

end function get_cmax
#:enddef  

#:def get_rho()
  pure real(dp) function get_rho(w, x) result(rho)
    ${GPU_ROUTINE_SEQ()}$
    real(dp), intent(in)  :: w(nw_phys)
    real(dp), intent(in)  :: x(1:ndim)

    ! TODO: only for local rad loss
  end function get_rho
#:enddef

#:def get_pthermal()
pure real(dp) function get_pthermal(w, x) result(pth)
  ${GPU_ROUTINE_SEQ()}$
  real(dp), intent(in)  :: w(nw_phys)
  real(dp), intent(in)  :: x(1:ndim)

    ! TODO: only for local rad loss
end function get_pthermal
#:enddef

#:def get_Rfactor()
pure real(dp) function get_Rfactor() result(Rfactor)
  ${GPU_ROUTINE_SEQ()}$

    ! TODO: only for local rad loss
end function get_Rfactor
#:enddef

#:def estimate_speeds_minmax()
!> Davis (1988) min/max wave speed estimates for HLL, in their relativistic
!> form: wL = min(lambda^-), wR = max(lambda^+) over the left and right states.
!>
!> The relativistic characteristic speeds in direction n are
!>
!>   lambda^+- = ( v_n (1 - c_s^2)
!>                 +- sqrt( c_s^2 (1 - v^2) (1 - v^2 c_s^2 - v_n^2 (1 - c_s^2)) ) )
!>               / (1 - v^2 c_s^2)
!>
!> (Anile 1989; Mignone & Bodo 2005), which reduce to v_n +- c_s in the
!> Newtonian limit. Unlike the HD and MHD versions of this routine, the
!> velocity is already contained in lambda^+-, so the min/max is taken over the
!> characteristic speeds themselves rather than over v_n +- c.
!>
!> This is the same algebra get_cmax uses; get_cmax collapses the two roots
!> into max(|lambda^+|, |lambda^-|) for the time-step limit, whereas HLL needs
!> them separately and signed. wprim(iw_mom(:)) holds the spatial four-velocity
!> u^i = lfac*v^i, so v^i is recovered as u^i/lfac.
subroutine estimate_speeds_minmax(uL, uR, xC, flux_dim, wL, wR)
  ${GPU_ROUTINE_SEQ()}$
  real(dp), intent(in)  :: uL(nw_phys), uR(nw_phys)
  real(dp), intent(in)  :: xC(ndim)
  integer, intent(in)   :: flux_dim
  real(dp), intent(out) :: wL, wR

  real(dp) :: rho, rhoh, pth, E, csound2
  real(dp) :: v2, vidim, root, denom
  real(dp) :: cminL, cmaxL, cminR, cmaxR

  ! Left state
  rho  = uL(iw_rho)
  rhoh = uL(xi_) / uL(lfac_)**2
  pth  = uL(iw_e)
  if (srhd_eos) then
     E = (rhoh + dsqrt(rhoh**2 + (srhd_gamma**2 - 1.0_dp)*rho**2)) &
          /(srhd_gamma + 1.0_dp)
     csound2 = (pth*((srhd_gamma + 1.0_dp) + gamma_1*(rho/E)**2))/(2.0_dp*rhoh)
  else
     csound2 = srhd_gamma*pth/rhoh
  end if

  v2    = 1.0_dp - 1.0_dp/uL(lfac_)**2
  vidim = uL(iw_mom(flux_dim))/uL(lfac_)
  denom = 1.0_dp - v2*csound2
  ! the radicand is positive analytically; clip the roundoff-level negatives
  ! that would otherwise produce a NaN wave speed
  root  = dsqrt(max(csound2*(1.0_dp - v2) &
       *(1.0_dp - v2*csound2 - vidim**2*(1.0_dp - csound2)), 0.0_dp))
  cmaxL = (vidim*(1.0_dp - csound2) + root)/denom
  cminL = (vidim*(1.0_dp - csound2) - root)/denom

  ! Right state
  rho  = uR(iw_rho)
  rhoh = uR(xi_) / uR(lfac_)**2
  pth  = uR(iw_e)
  if (srhd_eos) then
     E = (rhoh + dsqrt(rhoh**2 + (srhd_gamma**2 - 1.0_dp)*rho**2)) &
          /(srhd_gamma + 1.0_dp)
     csound2 = (pth*((srhd_gamma + 1.0_dp) + gamma_1*(rho/E)**2))/(2.0_dp*rhoh)
  else
     csound2 = srhd_gamma*pth/rhoh
  end if

  v2    = 1.0_dp - 1.0_dp/uR(lfac_)**2
  vidim = uR(iw_mom(flux_dim))/uR(lfac_)
  denom = 1.0_dp - v2*csound2
  root  = dsqrt(max(csound2*(1.0_dp - v2) &
       *(1.0_dp - v2*csound2 - vidim**2*(1.0_dp - csound2)), 0.0_dp))
  cmaxR = (vidim*(1.0_dp - csound2) + root)/denom
  cminR = (vidim*(1.0_dp - csound2) - root)/denom

  wL = min(cminL, cminR)
  wR = max(cmaxL, cmaxR)

  ! no signal can outrun light
  wL = min(max(wL, -1.0_dp), 1.0_dp)
  wR = min(max(wR, -1.0_dp), 1.0_dp)

end subroutine estimate_speeds_minmax
#:enddef


#:endif
