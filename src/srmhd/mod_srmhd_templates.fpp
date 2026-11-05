#:mute
#:include "../mod_gpu_directives.fpp"
#:endmute

#:if PHYS == 'srmhd'

#:if defined('N_TRACER')
#:set N_TRACER_ = N_TRACER
#:else
#:set N_TRACER_ = 0
#:endif


#:def phys_vars()

  integer, parameter :: dp = kind(0.0d0)
  integer, parameter, public              :: nw_phys=2+2*ndim+3 +${N_TRACER_}$
  integer, parameter, public              :: nw_flux=2+2*ndim+1 +${N_TRACER_}$
  
  !> Whether synge eos is used
  logical, public                         :: srmhd_eos = .false.
  ${GPU_DECLARE_COPYIN('srmhd_eos')}$

  !> Number of leading variables the MUSCL reconstruction limits (1:nw_recon).
  !> Slots nw_recon+1:nw_phys are left to fix_prim_state, which must rebuild
  !> them before the Riemann solver reads the face state.
  integer, parameter, public              :: nw_recon=nw_flux

  !> Index of the density (in the w array) as primitive or conserved
  integer, public                         :: rho_
  integer, public                         :: d_
  ${GPU_DECLARE_CREATE('rho_,d_')}$

  !> Indices of the momentum density
  integer, allocatable, public            :: mom(:)
  ${GPU_DECLARE_CREATE('mom')}$

  #:if defined('N_TRACER')
  !> Indices of the tracers
  integer, public                :: tracer(${N_TRACER_}$)
  ${GPU_DECLARE_CREATE('tracer')}$
  #:endif 

  !> Indices of the magnetic field
  integer, allocatable, public            :: mag(:)
  ${GPU_DECLARE_CREATE('mag')}$


  !> Index of the energy density
  integer, public                         :: tau_
  ${GPU_DECLARE_CREATE('tau_')}$

  !> Index of the gas pressure should equal tau_
  integer, public                         :: p_
  ${GPU_DECLARE_CREATE('p_')}$


  !> Index of the auxiliary variable mu (1/(lfac*h))
  integer, public            :: mu_
  ${GPU_DECLARE_CREATE('mu_')}$




  !> Index of the Lorentz factor
  integer, public     :: lfac_
  ${GPU_DECLARE_CREATE('lfac_')}$

  !> Index of GLM psi
  integer, public :: psi_
  ${GPU_DECLARE_CREATE('psi_')}$

  !> GLM-MHD parameter: ratio of the diffusive and advective time scales for div b
  !> taking values within [0, 1]
  double precision, public                :: srmhd_glm_alpha = 0.5d0
  ${GPU_DECLARE_COPYIN('srmhd_glm_alpha')}$


  !> Number of tracer species
  integer, public                         :: srmhd_n_tracer = 0
  ${GPU_DECLARE_COPYIN('srmhd_n_tracer')}$

  !> The adiabatic index
  double precision, public                :: srmhd_gamma = 5.d0/3.0d0
  ${GPU_DECLARE_COPYIN('srmhd_gamma')}$

  !> Helium abundance over Hydrogen
  double precision, public  :: He_abundance=0.1d0
  ${GPU_DECLARE_COPYIN('He_abundance')}$

  !> Whether particles module is added
  logical, public                         :: srmhd_particles = .false.
  ${GPU_DECLARE_COPYIN('srmhd_particles')}$

  !> switch for source user
  logical, public                         :: srmhd_source_usr = .false.
  ${GPU_DECLARE_COPYIN('srmhd_source_usr')}$

  !> Smallest gas pressure allowed: the floor fix_prim_state applies to a
  !> reconstructed face state. Zero, the default, leaves it unfloored.
  double precision, public                :: srmhd_small_pressure = 0.0d0
  !> Smallest rest-mass density allowed
  double precision, public                :: srmhd_small_density  = 1e-14_dp
  ${GPU_DECLARE_COPYIN('srmhd_small_pressure,srmhd_small_density')}$



#:enddef


#:def read_params()
    !> Read this module's parameters from a file
  subroutine read_params(files)
    use mod_global_parameters
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /srmhd_list/ srmhd_eos,srmhd_gamma, &
      He_abundance, srmhd_source_usr, &
      srmhd_small_pressure, srmhd_small_density, srmhd_glm_alpha, srmhd_n_tracer

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, srmhd_list, end=111)
111    close(unitpar)
    end do

    if (srmhd_small_pressure < 0.0d0) call mpistop(&
       "srmhd_small_pressure should be positive.")
    if (srmhd_small_density < 0.0d0) call mpistop(&
       "srmhd_small_density should be positive.")



    ${GPU_UPDATE_DEVICE('srmhd_eos, srmhd_gamma, He_abundance, srmhd_source_usr, srmhd_glm_alpha, srmhd_n_tracer')}$
    ${GPU_UPDATE_DEVICE('srmhd_small_pressure, srmhd_small_density')}$

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
    phys_gamma = srmhd_gamma


    phys_internal_e=.false.
    phys_partial_ionization=.false.
    need_global_cmax=.true.

    ! Whether diagonal ghost cells are required for the physics
    phys_req_diagonal = .false.

 ${GPU_UPDATE_DEVICE('physics_type, phys_energy, phys_total_energy, phys_internal_e, phys_gamma, phys_partial_ionization,need_global_cmax,phys_req_diagonal')}$

    use_particles = srmhd_particles

    ! Determine flux variables
    rho_ = var_set_rho()
    d_=rho_
    ${GPU_UPDATE_DEVICE('rho_,d_')}$

    allocate(mom(ndir))
    mom(:) = var_set_momentum(ndir)
    ${GPU_UPDATE_DEVICE('mom')}$

    ! Set index of energy variable
    tau_ = var_set_energy()
    p_ = tau_
    ${GPU_UPDATE_DEVICE('tau_,p_')}$

    ! set b field indices
    allocate(mag(ndir))
    mag(:) = var_set_bfield(ndir)
    ${GPU_UPDATE_DEVICE('mag')}$

    psi_ = var_set_fluxvar('psi', 'psi', need_bc=.false.)
    ${GPU_UPDATE_DEVICE('psi_')}$

    !> GLM MHD uses split source addition in psi:
    any_source_split = .true.
    ${GPU_UPDATE_DEVICE('any_source_split')}$

     ! Register tracer fields
    #:if defined('N_TRACER')
    #:for i in range(1, N_TRACER_+1)
        tracer(${i}$) = var_set_fluxvar("trc", "trp", ${i}$, need_bc=.false.)
    #:endfor
    ${GPU_UPDATE_DEVICE('tracer')}$
    #:endif


    ! Set index for auxiliary variables
    ! MUST be after the possible tracers (which have fluxes)
    mu_  = var_set_auxvar('mu','mu')
    lfac_= var_set_auxvar('lfac','lfac')

    if (nwflux > nw_recon .or. mu_ <= nw_recon .or. lfac_ <= nw_recon) &
      call mpistop("srhd phys_init: flux variables must occupy 1:nw_recon "//&
      "and the auxiliaries mu, lfac come after them")

    ${GPU_UPDATE_DEVICE('mu_,lfac_')}$




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

    ! BnormLF fix:
    flux_type(:,psi_)=flux_tvdlf
    do idir=1,ndir
       flux_type(idir,mag(idir))=flux_tvdlf
    end do
    ${GPU_UPDATE_DEVICE('flux_type')}$

    nvector      = 2 ! No. vector vars
    allocate(iw_vector(nvector))
    iw_vector(1) = mom(1) - 1
    iw_vector(2) = mag(1) - 1
    ${GPU_UPDATE_DEVICE('nvector, iw_vector')}$

! use cycle, needs to be dealt with:
!    ! Initialize particles module
!    if (srmhd_particles) then
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

#:def get_aux()
  !> The auxiliary variables of a conserved state, together with the pieces
  !> of the con2prim solution that to_primitive needs on top of them.
  !>
  !> This is the single definition of mu and lfac: to_primitive and set_aux
  !> both call it, so the two cannot drift apart. mu comes out of the root
  !> solve; lfac out of the velocity it implies, after that velocity has been
  !> limited to the bound v0 set by the conserved state.
  subroutine get_aux(u, mu, lfac, x, v_sqr, vel)
    use srmhd_con2prim

    ${GPU_ROUTINE_SEQ()}$
    real(dp), intent(in)  :: u(nw_phys)
    real(dp), intent(out) :: mu, lfac, x, v_sqr
    real(dp), intent(out) :: vel(3)
    real(dp) :: d, tau
    real(dp) :: s_sqr, b_sqr, s_dot_b, v0_sqr, v_rescale

    s_sqr = u(iw_mom(1))**2 + u(iw_mom(2))**2 + u(iw_mom(3))**2

    b_sqr = u(iw_mag(1))**2 + u(iw_mag(2))**2 + u(iw_mag(3))**2

    s_dot_b = u(iw_mom(1)) * u(iw_mag(1)) + u(iw_mom(2)) * u(iw_mag(2))& 
             + u(iw_mom(3)) * u(iw_mag(3))

    d=u(iw_rho)

    tau=u(iw_e)

    mu=u(mu_)

    call con2prim(mu, d, tau, s_sqr, b_sqr, s_dot_b)

    x=1/(1 + mu*b_sqr/d) !(26)

    vel(1)= mu*x *(u(iw_mom(1))/d + mu*s_dot_b*u(iw_mag(1))/d**2)

    vel(2)= mu*x *(u(iw_mom(2))/d + mu*s_dot_b*u(iw_mag(2))/d**2)

    vel(3)= mu*x *(u(iw_mom(3))/d + mu*s_dot_b*u(iw_mag(3))/d**2)

    v_sqr = vel(1)**2 + vel(2)**2 + vel(3)**2

    v0_sqr= s_sqr/(h0**2*d**2 + s_sqr)

    if (v_sqr>v0_sqr) then
      v_rescale= sqrt(v0_sqr/v_sqr)

      vel(1) = vel(1)*v_rescale

      vel(2) = vel(2)*v_rescale

      vel(3) = vel(3)*v_rescale

      v_sqr = v0_sqr
    end if

    lfac = 1/sqrt(1-v_sqr)

  end subroutine get_aux
#:enddef

#:def set_aux()
  !> Overwrite the auxiliary slots (mu, lfac) of a conserved state with the
  !> values consistent with its other components. Nothing else in u changes.
  subroutine set_aux(u)
    ${GPU_ROUTINE_SEQ()}$
    real(dp), intent(inout) :: u(nw_phys)
    real(dp) :: mu, lfac, x, v_sqr, vel(3)

    call get_aux(u, mu, lfac, x, v_sqr, vel)

    u(lfac_) = lfac

    u(mu_) = mu

  end subroutine set_aux
#:enddef

#:def to_primitive()
  subroutine to_primitive(u)
    use srmhd_con2prim, only: ideal_eos, synge_eos

    ${GPU_ROUTINE_SEQ()}$
    real(dp), intent(inout) :: u(nw_phys)
    real(dp) :: d, x, tau, lfac
    real(dp) :: mu
    real(dp) :: vel(3)
    real(dp) :: s_sqr, r_bar_sqr, q_bar, b_sqr, s_dot_b, eps, v_sqr

    call get_aux(u, mu, lfac, x, v_sqr, vel)

    s_sqr = u(iw_mom(1))**2 + u(iw_mom(2))**2 + u(iw_mom(3))**2

    b_sqr = u(iw_mag(1))**2 + u(iw_mag(2))**2 + u(iw_mag(3))**2

    s_dot_b = u(iw_mom(1)) * u(iw_mag(1)) + u(iw_mom(2)) * u(iw_mag(2))& 
             + u(iw_mom(3)) * u(iw_mag(3))

    d=u(iw_rho)

    tau=u(iw_e)

    r_bar_sqr = s_sqr*x**2/d**2 + mu * x * (1+x) * s_dot_b**2/d**3 !(38)

    q_bar =  tau/d - 0.5_dp*b_sqr/d - 0.5_dp* mu**2 * x**2 * (s_sqr*b_sqr/d**3-s_dot_b**2/d**3) !(39)


    eps = lfac*(q_bar - mu*r_bar_sqr) + v_sqr*lfac**2/(1 + lfac) !(42)

    if (eps<0) then
      eps=0
    end if

    u(iw_rho) = d/lfac  

    if (srmhd_eos) then
      u(iw_e) = synge_eos(u(iw_rho), eps)
    else
      u(iw_e) = ideal_eos(u(iw_rho), eps)
    end if 

    u(iw_mom(1)) = vel(1)*lfac

    u(iw_mom(2)) = vel(2)*lfac

    u(iw_mom(3)) = vel(3)*lfac

    u(lfac_) = lfac

    u(mu_) = mu
    
  end subroutine to_primitive
#:enddef


#:def to_conservative()
  subroutine to_conservative(u)
    ${GPU_ROUTINE_SEQ()}$
    real(dp), intent(inout) :: u(nw_phys)

    real(dp)    ::  v_sqr, b_sqr, v_dot_b

    real(dp)    :: E_th, E, h




    

    u(lfac_) = sqrt(1+u(iw_mom(1))**2 + u(iw_mom(2))**2 + u(iw_mom(3))**2)

    if (srmhd_eos) then 
      E_th = u(iw_e) * 1/(srmhd_gamma-1)
      E    = E_th + dsqrt(E_th**2+u(iw_rho)**2)
      ! writing rho/E on purpose, for numerics 
      h = 0.5_dp * ((srmhd_gamma+1.0_dp) * E &
          - (srmhd_gamma-1) * (u(iw_rho)/E))

      u(mu_) = 1/u(lfac_)/h
    else
      u(mu_) = 1/u(lfac_)/(1 + srmhd_gamma / (srmhd_gamma-1) / u(iw_rho) * u(iw_e))
    end if


    u(iw_rho)= u(iw_rho)*u(lfac_)

    v_sqr = (u(iw_mom(1))**2 + u(iw_mom(2))**2 + u(iw_mom(3))**2)/u(lfac_)**2

    b_sqr= u(iw_mag(1))**2 + u(iw_mag(2))**2 + u(iw_mag(3))**2

    v_dot_b = (u(iw_mag(1))*u(iw_mom(1)) + u(iw_mag(2))*u(iw_mom(2)) + u(iw_mag(3))*u(iw_mom(3)))/u(lfac_)

    u(iw_mom(1)) = (u(iw_rho)/u(mu_) + b_sqr)*u(iw_mom(1))/u(lfac_) - v_dot_b*u(iw_mag(1))

    u(iw_mom(2)) = (u(iw_rho)/u(mu_) + b_sqr)*u(iw_mom(2))/u(lfac_) - v_dot_b*u(iw_mag(2))

    u(iw_mom(3)) = (u(iw_rho)/u(mu_) + b_sqr)*u(iw_mom(3))/u(lfac_) - v_dot_b*u(iw_mag(3))

    u(iw_e) = u(iw_rho)*(1 / u(mu_) - 1) - u(iw_e) + 0.5_dp*b_sqr * (1 + v_sqr)&
              - 0.5_dp*v_dot_b**2

  
  end subroutine to_conservative
#:enddef

!input primitive
#:def get_flux()
  subroutine get_flux(u, xC, flux_dim, flux)
    use mod_global_parameters, only: cmax_global
    ${GPU_ROUTINE_SEQ()}$
    real(dp), intent(in)  :: u(nw_phys)
    real(dp), intent(in)  :: xC(1:ndim)
    integer, intent(in)   :: flux_dim
    real(dp), intent(out) :: flux(nw_flux)

    
    real(dp) :: vel(3), si
    real(dp) :: v_sqr, b_sqr, v_dot_b


    vel(1) = u(iw_mom(1)) / u(lfac_)
    vel(2) = u(iw_mom(2)) / u(lfac_)
    vel(3) = u(iw_mom(3)) / u(lfac_) 

    v_sqr = vel(1)**2 + vel(2)**2 + vel(3)**2 

    b_sqr= u(iw_mag(1))**2 + u(iw_mag(2))**2 + u(iw_mag(3))**2

    v_dot_b = u(iw_mag(1))*vel(1) + u(iw_mag(2))*vel(2) + u(iw_mag(3))*vel(3)


    si = (u(iw_rho) / u(mu_) * u(lfac_) + b_sqr)*vel(flux_dim) - v_dot_b*u(iw_mag(flux_dim))

    ! density flux
    flux(iw_rho) = u(iw_rho) * u(lfac_) * vel(flux_dim)

    ! momentum flux
    flux(iw_mom(1)) = si*vel(1) - u(iw_mag(1)) * u(iw_mag(flux_dim)) / u(lfac_)**2 - v_dot_b* vel(flux_dim)* u(iw_mag(1))

    flux(iw_mom(2)) = si*vel(2) - u(iw_mag(2)) * u(iw_mag(flux_dim)) / u(lfac_)**2 - v_dot_b* vel(flux_dim)* u(iw_mag(2))

    flux(iw_mom(3)) = si*vel(3) - u(iw_mag(3)) * u(iw_mag(flux_dim)) / u(lfac_)**2 - v_dot_b* vel(flux_dim)* u(iw_mag(3))

    flux(iw_mom(flux_dim)) = flux(iw_mom(flux_dim)) + u(iw_e) + 0.5_dp*(b_sqr/u(lfac_)**2 + v_dot_b**2)

    ! energy flux
    flux(iw_e) = si - vel(flux_dim)*u(iw_rho)*u(lfac_)

    ! magnetic flux
    flux(iw_mag(1)) = vel(flux_dim) * u(iw_mag(1)) - vel(1) * u(iw_mag(flux_dim))

    flux(iw_mag(2)) = vel(flux_dim) * u(iw_mag(2)) - vel(2) * u(iw_mag(flux_dim))

    flux(iw_mag(3)) = vel(flux_dim) * u(iw_mag(3)) - vel(3) * u(iw_mag(flux_dim))

    !GLM psi flux
    flux(iw_mag(flux_dim))=u(psi_)
      !f_i[psi]=Ch^2*b_{i} Eq. 24e and Eq. 38c Dedner et al 2002 JCP, 175, 645
    flux(psi_)=cmax_global**2*u(iw_mag(flux_dim))

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
pure function get_cmax(u, x, flux_dim) result(wC)
  ${GPU_ROUTINE_SEQ()}$
  real(dp), intent(in)  :: u(nw_phys)
  real(dp), intent(in)  :: x(1:ndim)
  integer, intent(in)   :: flux_dim

  real(dp) :: rhoh, b_sqr, c_s_sqr, small_b_sqr, v_dot_b, v_a_sqr
  real(dp) :: cms_sqr
  real(dp) :: A, B, C, root
  real(dp) :: cmax, cmin

  real(dp) :: wC

  rhoh = u(iw_rho) /u(mu_) /u(lfac_)
  
  c_s_sqr = srmhd_gamma*u(iw_e)/rhoh

  b_sqr = u(iw_mag(1))**2 + u(iw_mag(2))**2 + u(iw_mag(3))**2

  v_dot_b = (u(iw_mag(1))*u(iw_mom(1)) + u(iw_mag(2))*u(iw_mom(2)) + u(iw_mag(3))*u(iw_mom(3)))/u(lfac_)

  small_b_sqr = b_sqr/u(lfac_)**2 + v_dot_b**2

  v_a_sqr = small_b_sqr/(small_b_sqr + rhoh)

  cms_sqr = v_a_sqr + c_s_sqr - v_a_sqr*c_s_sqr

  cms_sqr = min(cms_sqr, 1.0d0)

  A = u(lfac_)**2 - (u(lfac_)**2 - 1)*cms_sqr

  B = 2*u(lfac_) * u(iw_mom(flux_dim)) * (1-cms_sqr)

  C = u(iw_mom(flux_dim))**2 - (1 + u(iw_mom(flux_dim))**2)*cms_sqr

  root = sqrt(max(B**2 - 4*A*C, 0.0d0))

  cmax = (B + root)/(2*A)

  cmin = (B - root)/(2*A)

  ! Limit by speed of light
  cmin = max(cmin, - 1.0d0)
  cmin = min(cmin,   1.0d0)
  cmax = max(cmax, - 1.0d0)
  cmax = min(cmax,   1.0d0)

  wC = max(abs(cmin), abs(cmax))

  ! wC=1

end function get_cmax
#:enddef  



#:def addsource_local()
subroutine addsource_local(qdt, dtfactor, qtC, wCT, wCTprim, qt, wnew, x, dr, &
    qsourcesplit)
#:if defined('SOURCE_USR')
  use mod_usr, only: addsource_usr
#:endif
  use mod_global_parameters, only: cmax_global
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

  else
     !---------------------------------
     ! split sources     
     !---------------------------------
        
    wnew(psi_)=wnew(psi_)*dexp(-qdt*cmax_global*srmhd_glm_alpha/minval(dr))
    
  end if

end subroutine addsource_local
#:enddef


#:def fix_prim_state()
  !> Make a reconstructed primitive face state admissible before it is handed
  !> to the Riemann solver.
  !>
  !> Two things happen here, and only the second is optional. The MUSCL
  !> reconstruction limits every slot of the state independently, including the
  !> auxiliaries mu and lfac -- but those are not independent variables, they
  !> are functions of (rho, u^i, p). A reconstructed face state therefore does
  !> not in general satisfy lfac = sqrt(1+|u|^2) or mu = 1/hW, while
  !> get_flux, get_cmax and estimate_speeds_minmax all read them. Left alone,
  !> v2 = 1-1/lfac^2 and v_n = u_n/lfac come from different sources, v_n can
  !> exceed v (even exceed 1), the radicand of the relativistic characteristic
  !> speeds goes negative and the HLL wave speeds collapse onto each other --
  !> which removes the solver's dissipation exactly where it is needed. So the
  !> auxiliaries are always recomputed from the reconstructed primitives, by
  !> the same expressions to_conservative uses.
  !>
  !> The floors are the optional part. srmhd_small_pressure defaults to zero
  !> and srmhd_small_density to 1e-14; set them in &srmhd_list. A pressure
  !> floored to exactly zero still gives the sound speed in get_cmax a zero
  !> numerator, so a positive srmhd_small_pressure is the useful
  !> configuration. Unlike srhd's, these floors are not read by con2prim.
  pure subroutine fix_prim_state(u)
    ${GPU_ROUTINE_SEQ()}$
    real(dp), intent(inout) :: u(nw_phys)

    u(iw_rho) = max(u(iw_rho), srmhd_small_density)
    u(iw_e)   = max(u(iw_e),   srmhd_small_pressure)

    u(lfac_) = sqrt(1+u(iw_mom(1))**2 + u(iw_mom(2))**2 + u(iw_mom(3))**2)

    u(mu_) = 1/u(lfac_)/(1 + srmhd_gamma / (srmhd_gamma-1) / u(iw_rho) * u(iw_e))

  end subroutine fix_prim_state
#:enddef


#:def estimate_speeds_minmax()
subroutine estimate_speeds_minmax(uL, uR, xC, flux_dim, wL, wR)
  ${GPU_ROUTINE_SEQ()}$
  real(dp), intent(in)  :: uL(nw_phys), uR(nw_phys)
  real(dp), intent(in)  :: xC(ndim)
  integer, intent(in)   :: flux_dim
  real(dp), intent(out) :: wL, wR

  real(dp) :: rhoh, b_sqr, c_s_sqr, small_b_sqr, v_dot_b, v_a_sqr
  real(dp) :: cms_sqr
  real(dp) :: A, B, C, root
  real(dp) :: cmaxL, cminL, cmaxR, cminR
  

  !left state 
  rhoh = uL(iw_rho) /uL(mu_) /uL(lfac_)
  
  c_s_sqr = srmhd_gamma*uL(iw_e)/rhoh

  b_sqr = uL(iw_mag(1))**2 + uL(iw_mag(2))**2 + uL(iw_mag(3))**2

  v_dot_b = (uL(iw_mag(1))*uL(iw_mom(1)) + uL(iw_mag(2))*uL(iw_mom(2)) + uL(iw_mag(3))*uL(iw_mom(3)))/uL(lfac_)

  small_b_sqr = b_sqr/uL(lfac_)**2 + v_dot_b**2

  v_a_sqr = small_b_sqr/(small_b_sqr + rhoh)

  cms_sqr = v_a_sqr + c_s_sqr - v_a_sqr*c_s_sqr

  cms_sqr = min(cms_sqr, 1.0d0)

  A = uL(lfac_)**2 - (uL(lfac_)**2 - 1)*cms_sqr

  B = 2*uL(lfac_) * uL(iw_mom(flux_dim)) * (1-cms_sqr)

  C = uL(iw_mom(flux_dim))**2 - (1 + uL(iw_mom(flux_dim))**2)*cms_sqr

  root = sqrt(max(B**2 - 4*A*C, 0.0d0))

  cmaxL = (B + root)/(2*A)

  cminL = (B - root)/(2*A)

  !right state
  rhoh = uR(iw_rho) /uR(mu_) /uR(lfac_)
  
  c_s_sqr = srmhd_gamma*uR(iw_e)/rhoh

  b_sqr = uR(iw_mag(1))**2 + uR(iw_mag(2))**2 + uR(iw_mag(3))**2

  v_dot_b = (uR(iw_mag(1))*uR(iw_mom(1)) + uR(iw_mag(2))*uR(iw_mom(2)) + uR(iw_mag(3))*uR(iw_mom(3)))/uR(lfac_)

  small_b_sqr = b_sqr/uR(lfac_)**2 + v_dot_b**2

  v_a_sqr = small_b_sqr/(small_b_sqr + rhoh)

  cms_sqr = v_a_sqr + c_s_sqr - v_a_sqr*c_s_sqr

  A = uR(lfac_)**2 - (uR(lfac_)**2 - 1)*cms_sqr

  B = 2*uR(lfac_) * uR(iw_mom(flux_dim)) * (1-cms_sqr)

  C = uR(iw_mom(flux_dim))**2 - (1 + uR(iw_mom(flux_dim))**2)*cms_sqr

  root = sqrt(max(B**2 - 4*A*C, 0.0d0))

  cmaxR = (B + root)/(2*A)

  cminR = (B - root)/(2*A)


  wL = min(cminL, cminR)

  wR = max(cmaxL, cmaxR)

  wL = min(max(wL, -1.0_dp), 1.0_dp)
  wR = min(max(wR, -1.0_dp), 1.0_dp)

end subroutine estimate_speeds_minmax
#:enddef

#:if GEOM == 'spherical'
!> Curvature source terms for spherical coordinates, with the pressure
!> gradient folded into the flux divergence, as in this fork's srhd and mhd
!> addsource_geometry. The non-pressure spatial stress is
!>
!>   T^ij = W v^i v^j - (v.B) (B^i v^j + v^i B^j) - B^i B^j / lfac^2,
!>   W    = rho h lfac^2 + B^2,
!>
!> the same tensor get_flux builds, and the total pressure is
!> ptot = p + (B^2/lfac^2 + (v.B)^2)/2. B here is the lab-frame field: this
!> T^ij is the covariant (rho h + b^2) u^i u^j - b^i b^j expanded with
!> b^i = B^i/lfac + lfac (v.B) v^i and b^2 = B^2/lfac^2 + (v.B)^2, where the
!> lfac^2 (v.B)^2 v^i v^j pieces cancel - so W carries the lab-frame B^2, not
!> the comoving b^2, and only ptot takes b^2. With B = 0 these reduce to srhd's
!> terms (W v^i v^j = xi u^i u^j / lfac^2); with lfac = 1 and the (v.B) terms
!> dropped, to mhd's. The induction and GLM psi terms are mhd's unchanged,
!> since the induction equation is the same in both. wprim(iw_mom(:)) holds
!> the spatial four-velocity u^i = lfac*v^i, as in srhd. Every geometric
!> prefactor is the discrete dAdV = (A_upper - A_lower)/dV rather than the
!> continuous 2/r and cot(theta)/r: for the metrics fill_geometry_device
!> builds, dAdV(1) = 2*<1/r> and dAdV(2) = <cot(theta)/r> exactly, so these
!> terms are well balanced and independent of where in the cell bgeo%x sits.
#:def addsource_geometry()
subroutine addsource_geometry(qdt, wprim, wnew, x, dAdV)
  ${GPU_ROUTINE_SEQ()}$

  real(dp), intent(in)     :: qdt
  !> primitive variables (rho, spatial four-velocity, pressure, B, psi, mu,
  !> lfac) at the current stage
  real(dp), intent(in)     :: wprim(nw_phys)
  !> cell-centre coordinates (r, theta, phi); unused - every geometric factor
  !> comes from dAdV, see the note above
  real(dp), intent(in)     :: x(1:ndim)
  !> (upper minus lower face area) / cell volume, per direction
  real(dp), intent(in)     :: dAdV(1:ndim)
  real(dp), intent(inout)  :: wnew(nw_phys)
  ! .. local ..
  real(dp)                 :: v1, v2, v3, b1, b2, b3, inv_lfac2, v_dot_b
  real(dp)                 :: w_tot, ptot, inv_r, cot_r, source
  real(dp)                 :: t22, t33, t12, t13, t23

  inv_lfac2 = 1.0_dp / wprim(lfac_)**2
  v1 = wprim(iw_mom(1)) / wprim(lfac_)
  v2 = wprim(iw_mom(2)) / wprim(lfac_)
  v3 = wprim(iw_mom(3)) / wprim(lfac_)
  b1 = wprim(iw_mag(1))
  b2 = wprim(iw_mag(2))
  b3 = wprim(iw_mag(3))
  v_dot_b = v1 * b1 + v2 * b2 + v3 * b3

  w_tot = wprim(iw_rho) * wprim(lfac_) / wprim(mu_) + b1**2 + b2**2 + b3**2
  ptot  = wprim(iw_e) + 0.5_dp * ((b1**2 + b2**2 + b3**2) * inv_lfac2 + &
     v_dot_b**2)
  inv_r = 0.5_dp * dAdV(1)   ! <1/r>
  cot_r = dAdV(2)            ! <cot(theta)/r>

  t22 = w_tot * v2 * v2 - 2.0_dp * v_dot_b * b2 * v2 - b2 * b2 * inv_lfac2
  t33 = w_tot * v3 * v3 - 2.0_dp * v_dot_b * b3 * v3 - b3 * b3 * inv_lfac2
  t12 = w_tot * v1 * v2 - v_dot_b * (b1 * v2 + v1 * b2) - b1 * b2 * inv_lfac2
  t13 = w_tot * v1 * v3 - v_dot_b * (b1 * v3 + v1 * b3) - b1 * b3 * inv_lfac2
  t23 = w_tot * v2 * v3 - v_dot_b * (b2 * v3 + v2 * b3) - b2 * b3 * inv_lfac2

  ! s[m_r] = (2 ptot + T_theta,theta + T_phi,phi) / r
  source = (2.0_dp * ptot + t22 + t33) * inv_r
  wnew(iw_mom(1)) = wnew(iw_mom(1)) + qdt * source

  ! s[m_theta] = ((ptot + T_phi,phi) cot(theta) - T_r,theta) / r
  source = (ptot + t33) * cot_r - t12 * inv_r
  wnew(iw_mom(2)) = wnew(iw_mom(2)) + qdt * source

  ! s[m_phi] = -(T_r,phi + T_theta,phi cot(theta)) / r
  source = -t13 * inv_r - t23 * cot_r
  wnew(iw_mom(3)) = wnew(iw_mom(3)) + qdt * source

  ! s[B_r] = psi * 2/r (GLM divergence-cleaning coupling only)
  source = 2.0_dp * wprim(psi_) * inv_r
  wnew(iw_mag(1)) = wnew(iw_mag(1)) + qdt * source

  ! s[B_theta] = (v_r B_theta - v_theta B_r + psi cot(theta)) / r
  source = (v1 * b2 - v2 * b1) * inv_r + wprim(psi_) * cot_r
  wnew(iw_mag(2)) = wnew(iw_mag(2)) + qdt * source

  ! s[B_phi] = (v_r B_phi - v_phi B_r
  !             + (v_theta B_phi - v_phi B_theta) cot(theta)) / r
  source = (v1 * b3 - v3 * b1) * inv_r + (v2 * b3 - v3 * b2) * cot_r
  wnew(iw_mag(3)) = wnew(iw_mag(3)) + qdt * source

  ! No geometric source term for D, tau or psi: they are scalars.

end subroutine addsource_geometry
#:enddef
#:elif GEOM == 'cylindrical'
!> Curvature source terms for cylindrical (r, z, phi) coordinates, with
!> the same stress T^ij, total pressure and conventions as the spherical
!> branch above. Only m_r, m_phi, B_phi and B_r (via the GLM psi coupling)
!> pick up curvature terms; m_z and B_z do not. dAdV(1) is exactly <1/r> -
!> and, since a cylindrical radial face area is linear in r, exactly
!> 1/r_midpoint as well. It is *not* 1/x(1): bgeo%x holds the volume
!> barycentre.
#:def addsource_geometry()
subroutine addsource_geometry(qdt, wprim, wnew, x, dAdV)
  ${GPU_ROUTINE_SEQ()}$

  real(dp), intent(in)     :: qdt
  !> primitive variables (rho, spatial four-velocity, pressure, B, psi, mu,
  !> lfac) at the current stage
  real(dp), intent(in)     :: wprim(nw_phys)
  !> cell-centre coordinates (r, z, phi); unused - every geometric factor
  !> comes from dAdV, see the note above
  real(dp), intent(in)     :: x(1:ndim)
  !> (upper minus lower face area) / cell volume, per direction
  real(dp), intent(in)     :: dAdV(1:ndim)
  real(dp), intent(inout)  :: wnew(nw_phys)
  ! .. local ..
  real(dp)                 :: v1, v2, v3, b1, b2, b3, inv_lfac2, v_dot_b
  real(dp)                 :: w_tot, ptot, inv_r, source
  real(dp)                 :: t33, t13

  inv_lfac2 = 1.0_dp / wprim(lfac_)**2
  v1 = wprim(iw_mom(1)) / wprim(lfac_)
  v2 = wprim(iw_mom(2)) / wprim(lfac_)
  v3 = wprim(iw_mom(3)) / wprim(lfac_)
  b1 = wprim(iw_mag(1))
  b2 = wprim(iw_mag(2))
  b3 = wprim(iw_mag(3))
  v_dot_b = v1 * b1 + v2 * b2 + v3 * b3

  w_tot = wprim(iw_rho) * wprim(lfac_) / wprim(mu_) + b1**2 + b2**2 + b3**2
  ptot  = wprim(iw_e) + 0.5_dp * ((b1**2 + b2**2 + b3**2) * inv_lfac2 + &
     v_dot_b**2)
  inv_r = dAdV(1)   ! <1/r>

  t33 = w_tot * v3 * v3 - 2.0_dp * v_dot_b * b3 * v3 - b3 * b3 * inv_lfac2
  t13 = w_tot * v1 * v3 - v_dot_b * (b1 * v3 + v1 * b3) - b1 * b3 * inv_lfac2

  ! s[m_r] = (ptot + T_phi,phi) / r
  source = (ptot + t33) * inv_r
  wnew(iw_mom(1)) = wnew(iw_mom(1)) + qdt * source

  ! s[m_phi] = -T_r,phi / r
  source = -t13 * inv_r
  wnew(iw_mom(3)) = wnew(iw_mom(3)) + qdt * source

  ! s[B_r] = psi / r (GLM divergence-cleaning coupling only)
  source = wprim(psi_) * inv_r
  wnew(iw_mag(1)) = wnew(iw_mag(1)) + qdt * source

  ! s[B_phi] = (B_phi v_r - B_r v_phi) / r
  source = (b3 * v1 - b1 * v3) * inv_r
  wnew(iw_mag(3)) = wnew(iw_mag(3)) + qdt * source

end subroutine addsource_geometry
#:enddef
#:endif





#:endif
