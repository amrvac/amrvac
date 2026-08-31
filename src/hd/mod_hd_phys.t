!> Hydrodynamics physics module
module mod_hd_phys
  use mod_thermal_conduction, only: tc_fluid
  use mod_radiative_cooling, only: rc_fluid
  use mod_thermal_emission, only: te_fluid
  use mod_fld, only: fld_fluid
  use mod_physics
  use mod_eos
  use mod_comm_lib, only: mpistop
  implicit none
  private

  !> Whether an energy equation is used
  logical, public, protected              :: hd_energy = .true.

  !> Whether thermal conduction is added
  logical, public, protected              :: hd_thermal_conduction = .false.
  !> Whether hyperbolic thermal conduction (Cattaneo relaxation) is used.
  !> 1D only : the q-variable is treated as a scalar carrying the
  !> heat flux along the only spatial direction.
  logical, public, protected              :: hd_hyperbolic_thermal_conduction = .false.
  !> Whether saturation is considered for hyperbolic TC
  logical, public, protected              :: hd_htc_sat = .false.
  type(tc_fluid), allocatable, public :: tc_fl
  type(te_fluid), allocatable, public :: te_fl_hd

  !> Whether radiative cooling is added
  logical, public, protected              :: hd_radiative_cooling = .false.
  type(rc_fluid), allocatable, public :: rc_fl

  !> Whether dust is added
  logical, public, protected              :: hd_dust = .false.

  !> Whether dust is added using and implicit update in IMEX
  logical, public, protected              :: hd_dust_implicit = .false.

  !> Whether radiation-gas interaction is handled using flux limited diffusion
  logical, public, protected              :: hd_radiation_fld = .false.
  logical, public, protected              :: hd_fld_pradtensor= .true.
  !> Radiation fluid object (gas-EoS callbacks for FLD), wired in hd_link_eos
  type(fld_fluid), allocatable, public    :: fld_fl

  !> Whether viscosity is added
  logical, public, protected              :: hd_viscosity = .false.

  !> Whether gravity is added
  logical, public, protected              :: hd_gravity = .false.

  !> Whether particles module is added
  logical, public, protected              :: hd_particles = .false.

  !> Whether rotating frame is activated
  logical, public, protected              :: hd_rotating_frame = .false.

  !> Whether CAK radiation line force is activated
  logical, public, protected              :: hd_cak_force = .false.

  !> Number of tracer species
  integer, public, protected              :: hd_n_tracer = 0

  !> Whether plasma is partially ionized

  !> Index of the density (in the w array)
  integer, public, protected              :: rho_

  !> Indices of the momentum density
  integer, allocatable, public, protected :: mom(:)

  !> Indices of the momentum density for the form of better vectorization
  integer, public, protected              :: ^C&m^C_

  !> Indices of the tracers
  integer, allocatable, public, protected :: tracer(:)

  !> Index of the energy density (-1 if not present)
  integer, public, protected              :: e_

  !> Index of the gas pressure (-1 if not present) should equal e_
  integer, public, protected              :: p_

  !> Index of the electron number density for LTE module
  integer, public, protected :: Ne_

  !> Index of the radiation energy (when fld active)
  integer, public, protected              :: r_e

  !> Indices of temperature
  integer, public, protected :: Te_

  !> Index of the FIP passive scalar rho*fip in conserved form, fip in primitive form
  integer, public, protected :: fip_ = -1

  !> Whether FIP passive scalar is enabled
  logical, public, protected :: hd_fip = .false.

  !> Index of the cutoff temperature for the TRAC method
  integer, public, protected              :: Tcoff_

  !> Index of the hyperbolic-TC heat-flux variable (-1 if not present)
  integer, public, protected              :: q_

  !> Thermal-conductivity prefactor in hyperbolic TC, set in hd_physical_units.
  !> Spitzer form: κ(T) = hypertc_kappa · T^{5/2}.
  double precision, public                :: hypertc_kappa

  !> Optional parfile override for hypertc_kappa (e.g. to match a constant-κ
  !> parabolic TC run for benchmarking). Default -1.0 leaves the Spitzer
  !> value computed from physical units.
  double precision, public, protected     :: hd_htc_kappa_override = -1.0d0

  !> Hyperdiffusion coefficient applied to the cell-refreshed q at the
  !> end of each face-recipe substep. 4th-order undivided difference
  !> smoother:  q -= alpha * (qdt/dt) * (q_{i+2} + q_{i-2}
  !>                                     - 4(q_{i+1} + q_{i-1})
  !>                                     + 6 q_i).
  !> MURaM (Rempel 2017) uses 0.02. Damps cell-scale q oscillations
  !> arising from kappa(T) amplifying small persistent T bumps in the
  !> EoS-derived cell-centred T profile; without it, the face-recipe
  !> q field has visible 5-cell-scale wiggles even where T is smooth.
  !> This is the only one of the OLD code's three q-smoothers (HLL
  !> diffusion, Koren reconstruction, hyperdiff) that the face-recipe
  !> needs -- the architectural problems with the other two are not
  !> reintroduced.
  double precision, public, protected     :: hd_htc_hyp_diff = 0.02d0

  !> Face-recipe heat-wave speed scaling: c_HTC,f = hd_htc_beta * c_max,f.
  !> Higher value -> closer to diffusion limit (q tracks Spitzer noise
  !> aggressively, AMR-triggering corona noise). Lower value -> more
  !> hyperbolic (q lags Spitzer target by Delta_t/tau ~ 1/beta^2 per
  !> step, dampening high-frequency noise). At our resolution beta=2-3
  !> is the practical sweet spot: q evolves slowly enough that T-table
  !> round-off noise doesn't propagate, but fast enough that real TR
  !> conduction equilibrates within O(100) hydro timesteps.
  double precision, public, protected     :: hd_htc_beta = 2.0d0

  !> Cowie-McKee saturation coefficient: q_sat = hd_htc_sat_alpha *
  !> rho * c_s^3. Standard convention is alpha ~ 1 (absorbs the
  !> sqrt(m_p/m_e) factor that would appear if c_s were replaced by
  !> the electron thermal speed). Default 1.0.
  double precision, public, protected     :: hd_htc_sat_alpha = 1.0d0

  !> Per-face energy-positivity safety fraction: |q_f^{n+1/2} dt A_f|
  !> <= hd_htc_pos_eta * min(e_int_L V_L, e_int_R V_R). Default 0.5
  !> leaves headroom against simultaneous PdV and cooling decrements.
  double precision, public, protected     :: hd_htc_pos_eta = 0.5d0

  !> Validity-monitor threshold for l_r,f / Delta_x_f. Warn if any face
  !> exceeds this in a given block (printed once per dtsave_log step).
  !> Default 0.1: above this, HTC is modifying the physics beyond pure
  !> Spitzer; above 1.0 the local Spitzer approximation breaks down.
  double precision, public, protected     :: hd_htc_validity_warn = 0.1d0

  !> Gradient deadband: zero out the Spitzer face flux when
  !>   abs(T_R - T_L) / max(T_L, T_R) < hd_htc_gradT_floor
  !> This suppresses sign-flipping q noise driven by EoS-table round-off
  !> (~1e-4 in T) being amplified by huge coronal kappa. Default 1.0e-3
  !> sits 10x above the table noise floor but 5x below typical coronal
  !> gradients (Delta_x / L_T ~ 5e-3 at our resolution), so it kills the
  !> noise without suppressing real conduction. Set to 0.0 to disable.
  double precision, public, protected     :: hd_htc_gradT_floor = 1.0d-3

  !> Running max of l_r,f / dx_f across all face-recipe calls since
  !> simulation start. Inspect post-hoc via debugger or dump alongside
  !> dat files. Not reset per timestep -- monotonic non-decreasing.
  double precision, public                :: hd_htc_validity_max_runtime = 0.0d0

  !> Index into wextra for escape probability column mass
  integer, public, protected              :: iw_colmass = -1

  !> gamma is set in &eos_list and accessed via eos%gamma

  !> The adiabatic constant
  double precision, public                :: hd_adiab = 1.0d0

  !> Whether TRAC method is used
  logical, public, protected              :: hd_trac = .false.
  integer, public, protected              :: hd_trac_type = 1
  integer, public, protected              :: hd_trac_nzones = 1
  double precision, public, protected     :: hd_trac_zone_splits(10) = -1.d0
  !> Johnston 2021 resolution parameter delta (default 0.5)
  double precision, public, protected     :: hd_trac_delta = 0.5d0
  !> Whether well-balanced reconstruction is used (Kaeppeli & Mishra style)
  logical, public, protected              :: hd_well_balanced = .false.


  !> Helium abundance over Hydrogen
  !> He_abundance is set in &eos_list and accessed via eos%He_abundance
  !> Ionization fraction of H
  !> H_ion_fr = H+/(H+ + H)
  double precision, public, protected  :: H_ion_fr=1d0
  !> Ionization fraction of He
  !> He_ion_fr = (He2+ + He+)/(He2+ + He+ + He)
  double precision, public, protected  :: He_ion_fr=1d0
  !> Ratio of number He2+ / number He+ + He2+
  !> He_ion_fr2 = He2+/(He2+ + He+)
  double precision, public, protected  :: He_ion_fr2=1d0
  ! used for eq of state when it is not defined by units,
  ! the units do not contain terms related to ionization fraction
  ! and it is p = RR * rho * T
  double precision, public, protected  :: RR=1d0
  ! remove the below flag  and assume default value = .false.
  ! procedure(sub_get_pthermal), pointer :: hd_get_Rfactor   => null()
  ! Public methods
  public :: hd_phys_init
  public :: hd_kin_en
  public :: hd_get_csound2
  public :: hd_check_params
  public :: hd_check_w
  public :: hd_handle_small_values
  public :: hd_e_to_ei
  public :: hd_ei_to_e
  ! hd_get_Rfactor was the FI pointer; FLD uses phys_get_Rfactor which is
  ! bound by mod_hd_eos:bind_eos_to_source to eos%get_Rfactor.
  ! Begin: following relevant for radiative hydro using FLD
  ! first four are local and only of interest for mod_usr applications
  ! where they can be used in diagnostics
  ! NOTE those with _prim expect primitives on entry
  public :: hd_get_pthermal_plus_pradiation
  public :: hd_get_csrad2
  public :: hd_get_trad
  public :: hd_get_pradiation_from_prim
  ! the following used in FLD modules
  !    as pointer phys_get_csrad2
  public :: hd_get_csrad2_prim
  ! End: following relevant for radiative hydro using FLD
  public :: hd_get_temperature_from_etot

contains

  !> Read this module's parameters from a file
  subroutine hd_read_params(files)
    use mod_global_parameters
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /hd_list/ hd_energy, hd_n_tracer, hd_adiab, &
    hd_dust, hd_dust_implicit, hd_thermal_conduction, &
    hd_hyperbolic_thermal_conduction, hd_htc_sat, &
    hd_htc_kappa_override, hd_htc_hyp_diff, &
    hd_htc_beta, hd_htc_sat_alpha, hd_htc_pos_eta, hd_htc_validity_warn, &
    hd_htc_gradT_floor, &
    hd_radiative_cooling, hd_viscosity, &
    hd_gravity, H_ion_fr, He_ion_fr, He_ion_fr2, &
    SI_unit, hd_particles, hd_rotating_frame, hd_trac, &
    hd_trac_type, hd_trac_nzones, hd_trac_zone_splits, hd_trac_delta, &
    hd_cak_force, hd_well_balanced, &
    hd_radiation_fld, hd_fld_pradtensor,hd_fip

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, hd_list, end=111)
111    close(unitpar)
    end do

  end subroutine hd_read_params

  !> Write this module's parameters to a snapsoht
  subroutine hd_write_info(fh)
    use mod_global_parameters
    integer, intent(in)                 :: fh
    integer, parameter                  :: n_par = 1
    double precision                    :: values(n_par)
    character(len=name_len)             :: names(n_par)
    integer, dimension(MPI_STATUS_SIZE) :: st
    integer                             :: er

    call MPI_FILE_WRITE(fh, n_par, 1, MPI_INTEGER, st, er)

    names(1) = "gamma"
    values(1) = eos%gamma
    call MPI_FILE_WRITE(fh, values, n_par, MPI_DOUBLE_PRECISION, st, er)
    call MPI_FILE_WRITE(fh, names, n_par * name_len, MPI_CHARACTER, st, er)
  end subroutine hd_write_info

  !> Initialize the module
  subroutine hd_phys_init()
    use mod_global_parameters
    use mod_thermal_conduction
    use mod_radiative_cooling
    use mod_dust, only: dust_init
    use mod_viscosity, only: viscosity_init
    use mod_gravity, only: gravity_init
    use mod_rotating_frame, only:rotating_frame_init
    use mod_cak_force, only: cak_init
    use mod_supertimestepping, only: sts_init, add_sts_method,&
            set_conversion_methods_to_head, set_error_handling_to_head
    use mod_eos_PI_tables
    use mod_usr_methods, only: usr_Rfactor, usr_get_heating
    use mod_escape_probability, only: escape_prob_init
    use mod_fld

    integer :: itr, idir

    call hd_read_params(par_files)

    physics_type = "hd"
    phys_energy  = hd_energy
    phys_total_energy  = hd_energy
    phys_internal_e = .false.
    phys_gamma = eos%gamma

    phys_trac=hd_trac
    if(phys_trac) then
      if(ndim .eq. 1) then
        if(hd_trac_type .gt. 2 .and. hd_trac_type .ne. 7) then
          hd_trac_type=1
          if(mype==0) write(*,*) 'WARNING: set hd_trac_type=1'
        end if
        if(hd_trac_type == 7) then
          if(.not. associated(usr_get_heating)) then
            call mpistop("hd_trac_type=7 requires usr_get_heating to be set in mod_usr.t")
          end if
        end if
        phys_trac_type=hd_trac_type
        phys_trac_nzones=hd_trac_nzones
        phys_trac_zone_splits=hd_trac_zone_splits
      else
        ! multi-D: only the local (per-cell) Johnston-2021 TRAC (type 7) is supported here;
        ! the column-sweep variants (1,2) need field-line/vertical tracing not implemented for HD.
        if(hd_trac_type == 7) then
          if(.not. associated(usr_get_heating)) then
            call mpistop("hd_trac_type=7 requires usr_get_heating to be set in mod_usr.t")
          end if
          phys_trac_type=7
        else
          phys_trac=.false.
          if(mype==0) write(*,*) 'WARNING: hd_trac disabled for ndim>=2 (only trac_type=7 supported)'
        end if
      end if
    end if

    ! set default gamma for polytropic/isothermal process
    if(.not.hd_energy) then
      if(hd_thermal_conduction) then
        hd_thermal_conduction=.false.
        if(mype==0) write(*,*) 'WARNING: set hd_thermal_conduction=F when hd_energy=F'
      end if
      if(hd_hyperbolic_thermal_conduction) then
        hd_hyperbolic_thermal_conduction=.false.
        if(mype==0) write(*,*) 'WARNING: set hd_hyperbolic_thermal_conduction=F when hd_energy=F'
      end if
      if(hd_radiative_cooling) then
        hd_radiative_cooling=.false.
        if(mype==0) write(*,*) 'WARNING: set hd_radiative_cooling=F when hd_energy=F'
      end if
    end if
    use_particles = hd_particles

    allocate(start_indices(number_species),stop_indices(number_species))

    ! set the index of the first flux variable for species 1
    start_indices(1)=1
    ! Determine flux variables
    rho_ = var_set_rho()

    allocate(mom(ndir))
    mom(:) = var_set_momentum(ndir)
    m^C_=mom(^C);

    ! Set index of energy variable
    if (hd_energy) then
       e_ = var_set_energy()
       p_ = e_
    else
       e_ = -1
       p_ = -1
    end if

    ! HTC q-variable must be allocated BEFORE any wextra slots (Ne, Te, ...)
    ! so its iw_q index doesn't collide with the wextra indices, which are
    ! assigned from nw (not nwflux). ffHD does the same ordering.
    if(hd_hyperbolic_thermal_conduction) then
      q_ = var_set_q()
      need_global_cmax=.true.
      ! q is a noisy diagnostic-like quantity (face-recipe construction
      ! produces cell-to-cell sign flips in nearly-uniform-T regions
      ! where -kappa*grad(T) is round-off-dominated). Excluding it from
      ! AMR refinement avoids triggering full-domain refinement everywhere.
      ! Lohner already refines on rho and E.
      if (allocated(w_refine_weight)) w_refine_weight(q_) = 0.0d0
    else
      q_=-1
    end if

    if (eos%eos_type == 'LTE') then
      Ne_ = var_set_ne()
      Te_ = var_set_te()
    else if (eos%eos_type == 'PI') then !  PI stores Te via var_set_te (sets iw_te) so the generic mod_eos_PI getters address it like LTE
      Ne_ = -1
      Te_ = var_set_te()
    else
      Ne_ = -1
      Te_ = -1
    end if

    if(hd_radiation_fld)then
       if(hd_cak_force)then
          if(mype==0) then
            write(*,*)'Warning: CAK force addition together with FLD radiation'
          endif
       endif
       if(hd_radiative_cooling)then
          if(mype==0) then
            write(*,*)'Warning: Optically thin cooling together with FLD radiation'
          endif
       endif
       if(hd_dust.and.hd_dust_implicit)then
          call mpistop('implicit dust addition not compatible with FLD radiation')
       endif
       if(.not.hd_energy)then
          call mpistop('using FLD implies the use of an energy equation, set hd_energy=T')
       else
          !> set added variable and equation for radiation energy
          r_e = var_set_radiation_energy()
          phys_get_csrad2          => hd_get_csrad2_prim
          !> Radiation fluid object: its EoS callbacks are wired in hd_link_eos
          allocate(fld_fl)
          !> Initiate radiation-closure module
          call fld_init()
          !> The implicit (MG diffusion) hooks need the fld_fl object, so they
          !> are wired here to physics-module wrappers that inject it.
          if(use_multigrid)then
             phys_implicit_update   => hd_fld_implicit_update
             phys_evaluate_implicit => hd_fld_evaluate_implicit
          endif
       endif
    else
      r_e=-1
    endif

    phys_get_dt              => hd_get_dt
    phys_get_cmax            => hd_get_cmax
    phys_get_tcutoff         => hd_get_tcutoff
    phys_get_cbounds         => hd_get_cbounds
    phys_get_flux            => hd_get_flux
    phys_add_source_geom     => hd_add_source_geom
    phys_add_source          => hd_add_source
    phys_modify_wLR          => hd_modify_wLR
    phys_check_params        => hd_check_params
    phys_check_w             => hd_check_w
    ! phys_get_pthermal is set by hd_link_eos
    phys_get_v               => hd_get_v
    ! phys_get_rho             => hd_get_rho
    phys_write_info          => hd_write_info
    phys_handle_small_values => hd_handle_small_values
    phys_e_to_ei            => hd_e_to_ei
    phys_ei_to_e            => hd_ei_to_e
    phys_get_ei             => hd_get_ei

    ! derive units from basic units
    call hd_physical_units()

    ! Spitzer prefactor for hyperbolic TC in code units. κ(T) = κ_0 · T^{5/2}
    ! with κ_0 = 8e-12 (SI) or 8e-7 (CGS), unless the user provides an override
    ! (e.g. to match a constant-κ parabolic TC run for benchmarking).
    if(hd_hyperbolic_thermal_conduction) then
      if (hd_htc_kappa_override > 0.0d0) then
        hypertc_kappa = hd_htc_kappa_override
        if (mype == 0) write(*,*) ' HTC: using hd_htc_kappa_override =', hypertc_kappa
      else if(SI_unit) then
        hypertc_kappa=8.d-12*unit_temperature**3.5d0/unit_length/unit_density/unit_velocity**3
      else
        hypertc_kappa=8.d-7*unit_temperature**3.5d0/unit_length/unit_density/unit_velocity**3
      end if
    end if

    if (hd_dust) then
        call dust_init(rho_, mom(:), e_)
    endif

    if (hd_fip) then
      fip_ = var_set_fluxvar('rho_fip', 'fip', need_bc=.false.)
    else
      fip_ = -1
    end if

    allocate(tracer(hd_n_tracer))

    ! Set starting index of tracers
    do itr = 1, hd_n_tracer
       tracer(itr) = var_set_fluxvar("trc", "trp", itr, need_bc=.false.)
    end do

    ! set number of variables which need update ghostcells.
    ! The EoS-derived slots Ne/Te are derived state that must stay consistent with (rho,e)
    ! wherever the conserved state is valid, so they have to travel with it. var_set_ne /
    ! var_set_te bump nw but neither nwflux nor nwaux, so the historic nwflux+nwaux silently
    ! excluded them: they were never communicated, only derived, and every ghost cell held
    ! zero until something derived it. Extend the window to whichever of them exist (they are
    ! registered contiguously just above); FI leaves both at -1 and the window is unchanged.
    ! bc_phys is unaffected -- it iterates nwflux+nwaux independently.
    nwgc=nwflux+nwaux
    if (iw_ne > 0) nwgc = max(nwgc, iw_ne)
    if (iw_te > 0) nwgc = max(nwgc, iw_te)

    ! set the index of the last flux variable for species 1
    stop_indices(1)=nwflux

    if(hd_trac) then
      Tcoff_ = var_set_wextra()
      iw_tcoff=Tcoff_
    else
      Tcoff_ = -1
    end if

    !> Cache log10(nH) in wextra for LTE+IonE TC (density invariant during STS)
    if (eos%eos_type == 'LTE' .and. eos%ionE .and. hd_thermal_conduction) then
        iw_log_nH = var_set_wextra()
    end if

    ! phys_get_Rfactor is bound by mod_hd_eos:bind_eos_to_source to eos%get_Rfactor.

    ! initialize thermal conduction module
    if (hd_thermal_conduction) then
      if (.not. hd_energy) &
           call mpistop("thermal conduction needs hd_energy=T")

      call sts_init()
      call tc_init_params(eos%gamma)

      allocate(tc_fl)
      call tc_get_hd_params(tc_fl,tc_params_read_hd)
      call add_sts_method(hd_get_tc_dt_hd,hd_sts_set_source_tc_hd,e_,1,e_,1,.false.)
      if (iw_log_nH > 0) then
        call set_conversion_methods_to_head(hd_e_to_ei_and_cache_log_nH, hd_ei_to_e)
      else
        call set_conversion_methods_to_head(hd_e_to_ei, hd_ei_to_e)
      end if
      call set_error_handling_to_head(hd_tc_handle_small_e)
      ! tc_fl%get_temperature_from_conserved => hd_get_temperature_from_etot
      ! tc_fl%get_temperature_from_eint => hd_get_temperature_from_eint
      ! tc_fl%get_temperature_from_conserved => eos%get_temperature_from_etot
      ! tc_fl%get_temperature_from_eint => eos%get_temperature_from_eint
      ! ! tc_fl%get_rho => hd_get_rho
      ! tc_fl%get_rho => eos%get_rho
      tc_fl%e_ = e_
      tc_fl%Tcoff_ = Tcoff_
    else if (hd_hyperbolic_thermal_conduction .and. hd_trac) then
      ! HTC + TRAC: allocate a stub tc_fl so TRAC type 7 can read tc_k_para.
      ! Use hypertc_kappa (Spitzer prefactor, same form as tc_k_para) so the
      ! TRAC kappa-effective formula reduces to the same expression as in
      ! the PTC branch. No STS init; no flux routine; this is read-only.
      call tc_init_params(eos%gamma)
      allocate(tc_fl)
      tc_fl%tc_k_para = hypertc_kappa
      tc_fl%tc_saturate = hd_htc_sat
      tc_fl%e_ = e_
      tc_fl%Tcoff_ = Tcoff_
    end if

    ! Initialize radiative cooling module
    if (hd_radiative_cooling) then
      if (.not. hd_energy) &
           call mpistop("radiative cooling needs hd_energy=T")
      call radiative_cooling_init_params(eos%gamma,eos%He_abundance)
      allocate(rc_fl)
      rc_fl%fip_ = fip_
      call radiative_cooling_init(rc_fl,rc_params_read)
      rc_fl%e_ = e_
      rc_fl%Tcoff_ = Tcoff_
      ! Initialize escape probability if requested
      if (rc_fl%rad_escape_prob) then
        iw_colmass = var_set_wextra()
        rc_fl%iw_colmass_ = iw_colmass
        phys_escape_prob = .true.
        call escape_prob_init(iw_colmass, rc_fl%rad_modify_sym, rc_fl%rad_escape_height)
      end if
    end if
    allocate(te_fl_hd)
    ! te_fl_hd%get_rho=> hd_get_rho
    te_fl_hd%get_rho=> eos%get_rho
    ! te_fl_hd%get_pthermal=> hd_get_pthermal
    te_fl_hd%get_pthermal=> eos%get_thermal_pressure
    te_fl_hd%get_var_Rfactor => eos%get_Rfactor
    te_fl_hd%get_ne_nH => eos%get_ne_nH


{^IFTHREED
    phys_te_images => hd_te_images
}
    ! Initialize viscosity module
    if (hd_viscosity) call viscosity_init(phys_wider_stencil)

    ! Initialize gravity module
    if (hd_gravity) call gravity_init()

    ! Well-balanced reconstruction: only meaningful with gravity
    if (hd_well_balanced) then
      if (.not. hd_gravity) then
        hd_well_balanced = .false.
        if(mype==0) write(*,*) 'WARNING: set hd_well_balanced=F (requires hd_gravity=T)'
      else
        phys_wb_transform => hd_wb_transform
        phys_wb_inverse   => hd_wb_inverse
        phys_wb_prolong   => hd_wb_prolong
        if (eos%ionE .and. eos%p2eint_method /= 'bisect' &
            .and. eos%method /= 'entropy') then
          eos%p2eint_method = 'bisect'
          if(mype==0) write(*,*) 'WB + ionE: forcing p2eint_method = bisect'
        end if
        if (eos%method == 'entropy' .and. mype == 0) then
          write(*,*) 'WB + ionE + entropy: p2eint_method stays "table"'
          write(*,*) 'eint_from_p_bisect uses legacy log_p table'
          write(*,*) 'not built for entropy method'
        end if
        if(mype==0) write(*,*) 'Well-balanced reconstruction enabled'
      end if
    end if

    ! Initialize rotating_frame module
    if (hd_rotating_frame) call rotating_frame_init()

    ! Initialize CAK radiation force module
    if (hd_cak_force) call cak_init(eos%gamma)


    ! Check whether custom flux types have been defined
    if (.not. allocated(flux_type)) then
       allocate(flux_type(ndir, nw))
       flux_type = flux_default
    else if (any(shape(flux_type) /= [ndir, nw])) then
       call mpistop("phys_check error: flux_type has wrong shape")
    end if

    nvector      = 1 ! No. vector vars
    allocate(iw_vector(nvector))
    iw_vector(1) = mom(1) - 1
    ! ionization-degree table init now lives in eos_finalise (eos% owns
    ! thermodynamic-backend init); see mod_eos_PI.

  end subroutine hd_phys_init

  !> allow the user to control the left/right states (and hence the flux) at
  !> physical boundary interfaces, as in the MHD module
  subroutine hd_modify_wLR(ixI^L,ixO^L,qt,wLC,wRC,wLp,wRp,s,idir)
    use mod_global_parameters
    use mod_usr_methods
    integer, intent(in)             :: ixI^L, ixO^L, idir
    double precision, intent(in)    :: qt
    double precision, intent(inout) :: wLC(ixI^S,1:nw), wRC(ixI^S,1:nw)
    double precision, intent(inout) :: wLp(ixI^S,1:nw), wRp(ixI^S,1:nw)
    type(state)                     :: s

    if(associated(usr_set_wLR)) call usr_set_wLR(ixI^L,ixO^L,qt,wLC,wRC,wLp,wRp,s,idir)

  end subroutine hd_modify_wLR

{^IFTHREED
  subroutine hd_te_images
    use mod_global_parameters
    use mod_thermal_emission
    select case(convert_type)
      case('EIvtiCCmpi','EIvtuCCmpi')
        call get_EUV_image(unitconvert,te_fl_hd)
      case('ESvtiCCmpi','ESvtuCCmpi')
        call get_EUV_spectrum(unitconvert,te_fl_hd)
      case('SIvtiCCmpi','SIvtuCCmpi')
        call get_SXR_image(unitconvert,te_fl_hd)
      case('WIvtiCCmpi','WIvtuCCmpi')
        call get_whitelight_image(unitconvert,te_fl_hd)
      case default
        call mpistop("Error in synthesize emission: Unknown convert_type")
      end select
  end subroutine hd_te_images
}
!!start th cond
  ! wrappers for STS functions in thermal_conductivity module
  ! which take as argument the tc_fluid (defined in the physics module)
  subroutine  hd_sts_set_source_tc_hd(ixI^L,ixO^L,w,x,wres,fix_conserve_at_step,my_dt,igrid,nflux)
    use mod_global_parameters
    use mod_fix_conserve
    use mod_thermal_conduction, only: sts_set_source_tc_hd
    integer, intent(in) :: ixI^L, ixO^L, igrid, nflux
    double precision, intent(in) ::  x(ixI^S,1:ndim)
    double precision, intent(inout) ::  wres(ixI^S,1:nw), w(ixI^S,1:nw)
    double precision, intent(in) :: my_dt
    logical, intent(in) :: fix_conserve_at_step
    call sts_set_source_tc_hd(ixI^L,ixO^L,w,x,wres,fix_conserve_at_step,my_dt,igrid,nflux,tc_fl)
  end subroutine hd_sts_set_source_tc_hd

  function hd_get_tc_dt_hd(w,ixI^L,ixO^L,dx^D,x) result(dtnew)
    !Check diffusion time limit dt < dx_i**2/((gamma-1)*tc_k_para/rho)
    !and   tc_k_para can depend on T=p/rho
    use mod_global_parameters
    use mod_thermal_conduction, only: get_tc_dt_hd
 
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: dx^D, x(ixI^S,1:ndim)
    double precision, intent(in) :: w(ixI^S,1:nw)
    double precision :: dtnew

    dtnew=get_tc_dt_hd(w,ixI^L,ixO^L,dx^D,x,tc_fl)
  end function hd_get_tc_dt_hd

  subroutine hd_tc_handle_small_e(w, x, ixI^L, ixO^L, step)
    ! move this in a different  routine as in mhd if needed in more places
    use mod_global_parameters
    use mod_small_values

    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    integer, intent(in)    :: step

    integer :: idir
    logical :: flag(ixI^S,1:nw)
    character(len=140) :: error_msg

    flag=.false.
    where(w(ixO^S,e_)<small_e) flag(ixO^S,e_)=.true.
    if(any(flag(ixO^S,e_))) then
      select case (small_values_method)
      case ("replace")
        where(flag(ixO^S,e_)) w(ixO^S,e_)=small_e
      case ("average")
        call small_values_average(ixI^L, ixO^L, w, x, flag, e_)
      case default
        ! small values error shows primitive variables
        w(ixO^S,e_)=w(ixO^S,e_)*(eos%gamma - 1.0d0)
        do idir = 1, ndir
           w(ixO^S, iw_mom(idir)) = w(ixO^S, iw_mom(idir))/w(ixO^S,rho_)
        end do
        write(error_msg,"(a,i3)") "Thermal conduction step ", step
        call small_values_error(w, x, ixI^L, ixO^L, flag, error_msg)
      end select
    end if
  end subroutine hd_tc_handle_small_e

    ! fill in tc_fluid fields from namelist
    subroutine tc_params_read_hd(fl)
      use mod_global_parameters, only: unitpar,par_files,unit_temperature
      type(tc_fluid), intent(inout) :: fl
      integer                      :: n
      logical :: tc_saturate=.false.
      logical :: tc_patch_eint=.false.
      double precision :: tc_k_para=0d0
      double precision :: trac_T_floor=1.d4

      namelist /tc_list/ tc_saturate, tc_k_para, trac_T_floor, tc_patch_eint

      do n = 1, size(par_files)
         open(unitpar, file=trim(par_files(n)), status="old")
         read(unitpar, tc_list, end=111)
111      close(unitpar)
      end do
      fl%tc_saturate = tc_saturate
      fl%tc_patch_eint = tc_patch_eint
      fl%tc_k_para = tc_k_para
      fl%trac_T_floor = trac_T_floor / unit_temperature

    end subroutine tc_params_read_hd

  ! subroutine hd_get_rho(w,x,ixI^L,ixO^L,rho)
  !   use mod_global_parameters
  !   integer, intent(in)           :: ixI^L, ixO^L
  !   double precision, intent(in)  :: w(ixI^S,1:nw),x(ixI^S,1:ndim)
  !   double precision, intent(out) :: rho(ixI^S)

  !   rho(ixO^S) = w(ixO^S,rho_) * eos%nH2rhoFactor

  ! end subroutine hd_get_rho

!!end th cond
!!rad cool
    subroutine rc_params_read(fl)
      use mod_global_parameters, only: unitpar,par_files,unit_temperature,unit_length
      use mod_constants, only: bigdouble
      use mod_basic_types, only: std_len
      type(rc_fluid), intent(inout) :: fl
      integer                      :: n
      ! list parameters
      integer :: ncool = 4000
    
      !> Name of cooling curve
      character(len=std_len)  :: coolcurve='JCcorona'
    
      !> Fixed temperature not lower than tlow
      logical    :: Tfix=.false.
    
      !> Lower limit of temperature
      double precision   :: tlow=bigdouble
    
      !> Add cooling source in a split way (.true.) or un-split way (.false.)
      logical    :: rc_split=.false.

      !> Cooling fraction (HEAD addition; used for explicit-mode dt scaling, kept for compat)
      double precision :: cfrac=0.1d0

      !> Master switch for radiative loss modification (spatial + density taper)
      logical :: rad_modify=.false.
      !> Apply spatial taper at both boundaries (default: lower only)
      logical :: rad_modify_sym=.false.
      !> Spatial taper: height from boundary below which taper applies
      double precision :: rad_cut_hgt=0.0d0
      !> Spatial taper: Gaussian decay width
      double precision :: rad_cut_dey=0.15d0
      !> Density taper: threshold above which taper applies
      double precision :: rad_taper_rho=bigdouble
      !> Density taper: Gaussian decay width
      double precision :: rad_taper_dey=0.0d0
      !> Suppress cooling below this temperature (Kelvin) within rad_cut_hgt.
      !> Cells inside rad_cut_hgt with T < rad_suppress_temp get factor=0.
      double precision :: rad_suppress_temp=0.0d0
      !> Enable escape probability cooling modification
      logical :: rad_escape_prob=.false.
      !> Effective opacity for escape probability (code units)
      double precision :: rad_kappa_eff=0.0d0
      !> Temperature above which kappa goes to 0 (Kelvin); 0 = constant kappa
      double precision :: rad_kappa_Tcutoff=0.0d0
      !> Sigmoid sharpness exponent for kappa(T) cutoff
      double precision :: rad_kappa_alpha=4.0d0
      !> Escape probability type: 'slab' or 'voigt'
      character(len=10) :: rad_escape_type='slab'
      !> Exponential cutoff scale: E *= exp(-tau/tau_cutoff); 0 = disabled
      double precision :: rad_escape_tau_cutoff=0.0d0
      !> Max height from footpoint for escape probability column mass (cm); 0 = no limit
      double precision :: rad_escape_height=0.0d0
      !> Variable-c_V Townsend extension (Y_mod): quadrature and sub-intervals
      character(len=8) :: rc_Y_mod_quadrature='boole'
      integer :: rc_Y_mod_N_sub=16
      !> Upstream: cutoff radiative cooling below rad_damp_height (Gaussian damp)
      logical          :: rad_damp=.false.
      double precision :: rad_damp_height=0.5d0
      double precision :: rad_damp_scale=0.15d0
      !> Upstream: Newton-radiative damping (surface treatment)
      logical          :: rad_newton=.false.
      double precision :: rad_newton_trad=0.006d0
      double precision :: rad_newton_rhosurf=1.d4
      double precision :: rad_newton_pthick=25.d0

      namelist /rc_list/ coolcurve, ncool, cfrac, tlow, Tfix, rc_split, &
          rad_modify, rad_modify_sym, rad_suppress_temp, &
          rad_cut_hgt, rad_cut_dey, rad_taper_rho, rad_taper_dey, &
          rad_escape_prob, rad_kappa_eff, rad_kappa_Tcutoff, rad_kappa_alpha, &
          rad_escape_type, rad_escape_tau_cutoff, rad_escape_height, &
          rc_Y_mod_quadrature, rc_Y_mod_N_sub, &
          rad_damp, rad_damp_height, rad_damp_scale, &
          rad_newton, rad_newton_trad, rad_newton_rhosurf, rad_newton_pthick

      do n = 1, size(par_files)
        open(unitpar, file=trim(par_files(n)), status="old")
        read(unitpar, rc_list, end=111)
111     close(unitpar)
      end do

      fl%ncool=ncool
      fl%coolcurve=coolcurve
      fl%tlow=tlow
      fl%Tfix=Tfix
      fl%rc_split=rc_split
      fl%cfrac=cfrac
      fl%rad_modify=rad_modify
      fl%rad_modify_sym=rad_modify_sym
      fl%rad_suppress_temp=rad_suppress_temp
      fl%rad_cut_hgt=rad_cut_hgt
      fl%rad_cut_dey=rad_cut_dey
      fl%rad_taper_rho=rad_taper_rho
      fl%rad_taper_dey=rad_taper_dey
      fl%rad_escape_prob=rad_escape_prob
      fl%rad_kappa_eff=rad_kappa_eff
      fl%rad_kappa_Tcutoff=rad_kappa_Tcutoff/unit_temperature
      fl%rad_kappa_alpha=rad_kappa_alpha
      fl%rad_escape_type=rad_escape_type
      fl%rad_escape_tau_cutoff=rad_escape_tau_cutoff
      fl%rad_escape_height=rad_escape_height/unit_length
      fl%Y_mod_quadrature=rc_Y_mod_quadrature
      fl%Y_mod_N_sub=rc_Y_mod_N_sub
      fl%rad_damp = rad_damp
      fl%rad_damp_height = rad_damp_height
      fl%rad_damp_scale = rad_damp_scale
      fl%rad_newton = rad_newton
      fl%rad_newton_trad = rad_newton_trad
      fl%rad_newton_rhosurf = rad_newton_rhosurf
      fl%rad_newton_pthick = rad_newton_pthick
    end subroutine rc_params_read
!! end rad cool

  subroutine hd_check_params
    use mod_global_parameters
    use mod_geometry, only: coordinate
    use mod_dust, only: dust_check_params, dust_implicit_update, dust_evaluate_implicit
    use mod_particles, only: particles_init
    use mod_particles, only: npayload,nusrpayload,ngridvars,num_particles,physics_type_particles
    use mod_fld

    double precision :: a,b,Xfrac,Yfrac

    ! Initialize particles module, put here so additional gridvars and user payloads are known
    if (hd_particles) then
       call particles_init()
    end if

    if (.not. hd_energy) then
       if (eos%gamma <= 0.0d0) call mpistop ("Error: eos%gamma <= 0")
       if (hd_adiab < 0.0d0) call mpistop  ("Error: hd_adiab < 0")
       small_pressure= hd_adiab*small_density**eos%gamma
    else
       if (eos%gamma <= 0.0d0 .or. eos%gamma == 1.0d0) &
            call mpistop ("Error: eos%gamma <= 0 or eos%gamma == 1.0")
       ! For LTE+ionE, this floor excludes ionisation energy. At tlow (~1000 K)
       ! ionisation is negligible, so the thermal-only floor is adequate.
       small_e = small_pressure * eos%inv_gamma_minus_1
       small_r_e = small_pressure * eos%inv_gamma_minus_1
       ! gamma_minus_1 and inv_gamma_minus_1 are set by eos_init
    end if

    if (hd_dust) call dust_check_params()

    if(hd_dust_implicit) then
        if(.not.use_imex_scheme)then
           call mpistop('select IMEX scheme for implicit dust update')
        endif
        ! implicit dust update
        phys_implicit_update => dust_implicit_update
        phys_evaluate_implicit => dust_evaluate_implicit
    endif

    ! Hyperbolic TC in HD is only implemented for 1D: q is a scalar carrying
    ! the heat flux along the single spatial direction. For ndim>1 the user
    ! should use the ffHD module (B-aligned conduction) or the parabolic
    ! tc_init path. Restrict here rather than silently giving wrong fluxes.
    if (hd_hyperbolic_thermal_conduction .and. ndim /= 1) then
      call mpistop("hd_hyperbolic_thermal_conduction is implemented for ndim=1 only;" // &
                   " for ndim>1 use mod_ffhd or parabolic mod_thermal_conduction.")
    end if
    if (hd_hyperbolic_thermal_conduction .and. hd_thermal_conduction) then
      call mpistop("hd_hyperbolic_thermal_conduction and hd_thermal_conduction are mutually exclusive;" // &
                   " choose one TC implementation.")
    end if

    if(hd_radiation_fld) then
        if(.not.use_imex_scheme)then
           call mpistop('select IMEX scheme for FLD radiation use')
        endif
        if(use_multigrid)then
           call phys_set_mg_bounds()
        else
           if(.not.fld_no_mg)call mpistop('multigrid must have BCs for IMEX and FLD radiation use')
        endif
        if(mype==0)then
           write(*,*)'==FLD SETUP======================'
           write(*,*)'Using FLD with settings:'
           write(*,*)'Using FLD with settings: hd_radiation_fld=',hd_radiation_fld
           write(*,*)'Using FLD with settings: hd_fld_pradtensor=',hd_fld_pradtensor
           write(*,*)'Using FLD with settings: fld_fluxlimiter=',fld_fluxlimiter
           write(*,*)'Using FLD with settings: fld_bound_diff=',fld_bound_diff
           write(*,*)'Using FLD with settings: fld_interaction_method=',fld_interaction_method
           write(*,*)'Using FLD with settings: fld_opacity_law=',fld_opacity_law
           write(*,*)'Using FLD with settings: fld_kappa0=',fld_kappa0
           write(*,*)'Using FLD with settings: fld_opal_table=',fld_opal_table
           write(*,*)'Using FLD with settings: fld_Radforce_split=',fld_Radforce_split
           write(*,*)'Using FLD with settings: fld_bisect_tol=',fld_bisect_tol
           write(*,*)'Using FLD with settings: fld_diff_tol=',fld_diff_tol
           write(*,*)'Using FLD with settings: nth_for_diff_mg=',nth_for_diff_mg
           write(*,*)'      FLD has use_imex_scheme and use_multigrid=',use_imex_scheme,use_multigrid
           write(*,*)'      FLD has fld_no_mg=',fld_no_mg
           if(fld_no_mg)then
              print *,'WARNING: cheating with FLD diffusion ***********************'
              print *,'WARNING: No MG-diffusion for radiative energy at all!!!!!!!'
              print *,'WARNING: cheating with FLD diffusion ***********************'
           endif
           print *,'const_rad_a   =',const_rad_a
           print *,'NORMALIZED arad_norm=',arad_norm
           print *,'NORMALIZED c_norm=',c_norm
           if(fld_cnorm>0.0d0)then
              print *,'WARNING: cheating with c_norm ***********************'
              print *,'WARNING: c_norm reset to=',fld_cnorm
              c_norm=fld_cnorm
              print *,'WARNING: cheating with c_norm ***********************'
           endif
           print *,'const_kappae  =',const_kappae
           if(trim(fld_opacity_law).eq.'const_norm')then
               print *,'NORMALIZED fld_kappa0          =',fld_kappa0
               print *,'physical value (in cgs or SI)  =',fld_kappa0*unit_opacity
           endif
           if(trim(fld_opacity_law).eq.'const')then
               print *,'physical fld_kappa (in cgs or SI) =',fld_kappa0
               print *,'NORMALIZED value                  =',fld_kappa0/unit_opacity
           endif
           write(*,*)'===FLD SETUP====================='
        endif
    endif
    if(mype==0)then
           write(*,*)'====HD run with settings===================='
           write(*,*)'Using mod_hd_phys with settings:'
           write(*,*)'SI_unit=',SI_unit
           write(*,*)'Dimensionality   :',ndim
           write(*,*)'vector components:',ndir
           write(*,*)'coordinate set to type,slab:',coordinate,slab
           write(*,*)'number of variables          nw=',nw
           write(*,*)'    start index         iwstart=',iwstart
           write(*,*)'number of      vector variables=',nvector
           write(*,*)'number of stagger variables nws=',nws
           write(*,*)'number of    variables with BCs=',nwgc
           write(*,*)'number of      vars with fluxes=',nwflux
           write(*,*)'number of   vars with flux + BC=',nwfluxbc
           write(*,*)'number of   auxiliary variables=',nwaux
           write(*,*)'number of extra vars without flux=',nwextra
           write(*,*)'number of extra vars   for wextra=',nw_extra
           write(*,*)'number of auxiliary I/O variables=',nwauxio
           write(*,*)'number of             hd_n_tracer=',hd_n_tracer
           write(*,*)'    hd_energy=',hd_energy
           write(*,*)'    hd_gravity=',hd_gravity
           write(*,*)'    hd_viscosity=',hd_viscosity
           write(*,*)'    hd_radiative_cooling=',hd_radiative_cooling
           write(*,*)'    hd_cak_force=',hd_cak_force
           write(*,*)'    hd_radiation_fld=',hd_radiation_fld
           write(*,*)'    hd_thermal_conduction=',hd_thermal_conduction
           write(*,*)'    hd_hyperbolic_thermal_conduction=',hd_hyperbolic_thermal_conduction
           write(*,*)'    hd_trac=',hd_trac
           write(*,*)'    hd_dust=',hd_dust
           write(*,*)'    hd_rotating_frame=',hd_rotating_frame
           write(*,*)'    hd_particles=',hd_particles
           if(hd_particles) then
              write(*,*) '*****Using particles: npayload,ngridvars :', npayload,ngridvars
              write(*,*) '*****Using particles:        nusrpayload :', nusrpayload
              write(*,*) '*****Using particles:      num_particles :', num_particles
              write(*,*) '*****Using particles: physics_type_particles=',physics_type_particles
           end if
           write(*,*)'number of             ghostcells=',nghostcells
           write(*,*)'number due to phys_wider_stencil=',phys_wider_stencil
           write(*,*)'==========================================='
           print *,'========EOS and UNITS==========='
           print *,'SI_unit       =',SI_unit
           print *,'gamma=',eos%gamma
           print *,'He_abundance  =',eos%He_abundance
           print *,'RR            =',RR
           print *,'========EOS and UNITS==========='
           print *,'unit_time          =',unit_time
           print *,'unit_length        =',unit_length
           print *,'unit_velocity      =',unit_velocity
           print *,'unit_pressure      =',unit_pressure
           print *,'unit_numberdensity =',unit_numberdensity
           print *,'unit_density       =',unit_density
           print *,'unit_temperature   =',unit_temperature
           print *,'unit_mass          =',unit_mass
           print *,'unit_Erad          =',unit_Erad
           print *,'unit_radflux       =',unit_radflux
           print *, 'CHECK that p_u ',unit_pressure,' equals ',unit_density*unit_velocity**2
           print *, 'CHECK that L_u ',unit_length,' equals ',unit_velocity*unit_time
           print *, 'CHECK that M_u',unit_mass,' equals ',unit_density*unit_length**3
           print *, 'density to numberdensity has factor   ',unit_density/unit_numberdensity
           if(SI_unit)then
                print *, '                     compare  this to ',mp_SI*(1.d0+4.d0*eos%He_abundance)
           else
                print *, '                     compare  this to ',mp_cgs*(1.d0+4.d0*eos%He_abundance)
           endif
           print *, 'pressure to n T has factor            ',unit_pressure/(unit_numberdensity*unit_temperature)
           if(SI_unit)then
                print *, '                     compare  this to ',kB_SI*(2.d0+3.d0*eos%He_abundance)
                a=unit_density/unit_numberdensity/mp_SI
                b=unit_pressure/(unit_numberdensity*unit_temperature*kB_SI)
           else
                print *, '                     compare  this to ',kB_cgs*(2.d0+3.d0*eos%He_abundance)
                a=unit_density/unit_numberdensity/mp_cgs
                b=unit_pressure/(unit_numberdensity*unit_temperature*kB_cgs)
           endif
           if(eos%eos_type /= 'LTE')then
              print *, 'mean molecular weight mu is =',a/b,' = ', (1.d0+4.d0*eos%He_abundance)/(2.d0+3.d0*eos%He_abundance)
              Xfrac=1.d0/a
              Yfrac=4.d0*eos%He_abundance/(1.d0+4.d0*eos%He_abundance)
              print *, 'mass fraction hydrogen X is =',1/a,' and this equals ', 1.d0/(1.d0+4.d0*eos%He_abundance)
              print *, 'mass fraction helium   Y is =',Yfrac
              print *, ' check that 1/mu', b/a,' is equal to 2X+3Y/4=',2.d0*Xfrac+3.d0*Yfrac/4.d0
              print *, ' ratio n_e/n_p=',1.d0+2.0d0*eos%He_abundance
           endif
           print *,'========UNITS==========='
    endif

  end subroutine hd_check_params

  subroutine hd_physical_units
    use mod_global_parameters
    double precision :: mp,kB,c_lightspeed,Xfrac,sigma_Telectron
    double precision :: a,b
    ! Derive scaling units
    if(SI_unit) then
      mp=mp_SI
      kB=kB_SI
      const_sigmaSB=sigma_SB_SI
      c_lightspeed=c_SI
      sigma_Telectron=sigma_Te_SI
    else
      mp=mp_cgs
      kB=kB_cgs
      const_sigmaSB=sigma_SB_cgs
      c_lightspeed=const_c
      sigma_Telectron=sigma_Te_cgs
    end if
    ! Normalisation dispatch keyed solely on eos%eos_type (FI is the default, so
    ! legacy parfiles land in the FI/PI absorbed-(a,b), RR=1 branch -- the former
    ! eq_state_units=.true. result).
    if (eos%eos_type == 'LTE') then
      !> Remove the assumed FI normalisation from the units and handle in EoS
      a=1d0
      b=1d0
      eos%nH2rhoFactor = 1d0+4d0*eos%He_abundance
      RR=(2d0+3d0*eos%He_abundance) / (1d0+4d0*eos%He_abundance)
      Xfrac=1.d0/(1.d0+4.d0*eos%He_abundance)
    else
      !> FI / PI: absorbed-(a,b), RR=1 (a=b=1 with RR=1 would be wrong physics
      !> for He>0). PI shares FI's normalisation exactly.
      a=1d0+4d0*eos%He_abundance
      if(eos%eos_type=='PI') then
        b=1d0+H_ion_fr+eos%He_abundance*(He_ion_fr*(He_ion_fr2+1d0)+1d0)
      else
        b=2d0+3d0*eos%He_abundance
      end if
      RR=1d0
      Xfrac=1.d0/a
    end if
    if(unit_density/=1.d0 .or. unit_numberdensity/=1.d0) then
      if(unit_density/=1.d0) then
        unit_numberdensity=unit_density/(a*mp)
      else if(unit_numberdensity/=1.d0) then
        unit_density=a*mp*unit_numberdensity
      end if
      if(unit_temperature/=1.d0) then
        unit_pressure=b*unit_numberdensity*kB*unit_temperature
        unit_velocity=sqrt(unit_pressure/unit_density)
        if(unit_length/=1.d0) then
          unit_time=unit_length/unit_velocity
        else if(unit_time/=1.d0) then
          unit_length=unit_velocity*unit_time
        end if
      else if(unit_pressure/=1.d0) then
        unit_temperature=unit_pressure/(b*unit_numberdensity*kB)
        unit_velocity=sqrt(unit_pressure/unit_density)
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
        unit_velocity=sqrt(unit_pressure/unit_density)
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
    unit_mass = unit_density * unit_length**3

    !> Units needed for radiative flux and opacity as used in FLD
    ! normalized light speed
    c_norm=c_lightspeed/unit_velocity
    ! this is the radiation constant in either cgs or SI units
    const_rad_a=4.d0*const_sigmaSB/c_lightspeed
    ! this is the dimensionless conversion factor for Erad to Trad
    arad_norm=const_rad_a*unit_temperature**4/unit_pressure
    ! This is the Thomson scattering opacity in the correct units
    ! hydrogen mass fraction X=1/a in the absorbed-(a,b) normalisation
    const_kappae=sigma_Telectron*(1.d0+Xfrac)/(2.0d0*mp)
    ! these are the units
    unit_Erad = unit_pressure
    unit_radflux = unit_velocity*unit_Erad
    unit_opacity = one/(unit_density*unit_length)

  end subroutine hd_physical_units

  !> Returns logical argument flag where values are ok
  subroutine hd_check_w(primitive, ixI^L, ixO^L, w, flag)
    use mod_global_parameters
    use mod_dust, only: dust_check_w

    logical, intent(in)          :: primitive
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, nw)
    logical, intent(inout)       :: flag(ixI^S,1:nw)
    
    double precision             :: tmp(ixI^S)
    double precision :: x(ixI^S, 1:ndim)

    flag=.false.
    x(ixI^S,1:ndim)=block%x(ixI^S,1:ndim) !> Rather than redefining the hd_check_w and limiter procedure interfaces
    if (hd_energy) then
       if (primitive) then
          where(w(ixO^S, e_) < small_pressure) flag(ixO^S,e_) = .true.
       else
          ! Inline (gamma-1)*(e - KE) to avoid eos%get_thermal_pressure side
          ! effects (fix_small_values clipping / crash=.true.) that would
          ! suppress the flag or abort before remediation could run.
          tmp(ixO^S)=(eos%gamma-1.0d0)*(w(ixO^S,e_)-&
           half*(^C&w(ixO^S,m^C_)**2+)/w(ixO^S,rho_))
          where(tmp(ixO^S) < small_pressure) flag(ixO^S,e_) = .true.
       endif
       if(hd_radiation_fld)then
          where(w(ixO^S, r_e) < small_r_e) flag(ixO^S,r_e) = .true.
       endif
    end if

    where(w(ixO^S, rho_) < small_density) flag(ixO^S,rho_) = .true.

    if(hd_dust) call dust_check_w(ixI^L,ixO^L,w,x,flag)

  end subroutine hd_check_w

  subroutine hd_bound_fip(primitive, ixI^L, ixO^L, w)
    use mod_global_parameters
    logical, intent(in)             :: primitive
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: rho_safe(ixI^S), fip_prim(ixI^S)

    if (.not. hd_fip) return

    if (primitive) then
      w(ixO^S,fip_) = min(maxfip, max(minfip, w(ixO^S,fip_)))
    else
      rho_safe(ixO^S) = max(w(ixO^S,rho_), small_density)
      fip_prim(ixO^S) = w(ixO^S,fip_) / rho_safe(ixO^S)
      fip_prim(ixO^S) = min(maxfip, max(minfip, fip_prim(ixO^S)))
      w(ixO^S,fip_) = rho_safe(ixO^S) * fip_prim(ixO^S)
    end if
  end subroutine hd_bound_fip

  ! !> Transform primitive variables into conservative ones
  ! subroutine hd_to_conserved(ixI^L, ixO^L, w, x)
  !   use mod_global_parameters
  !   use mod_dust, only: dust_to_conserved
  !   integer, intent(in)             :: ixI^L, ixO^L
  !   double precision, intent(inout) :: w(ixI^S, nw)
  !   double precision, intent(in)    :: x(ixI^S, 1:ndim)

  !   integer :: ix^D

  !   {do ix^DB=ixOmin^DB,ixOmax^DB\}
  !     if (hd_energy) then
  !        ! Calculate total energy from pressure and kinetic energy
  !        w(ix^D,e_)=w(ix^D, e_)*inv_gamma_1+&
  !         half*(^C&w(ix^D,m^C_)**2+)*w(ix^D,rho_)
  !     end if
  !     ! Convert velocity to momentum
  !     ^C&w(ix^D,m^C_)=w(ix^D,rho_)*w(ix^D,m^C_)\
  !   {end do\}

  !   if (hd_dust) then
  !     call dust_to_conserved(ixI^L, ixO^L, w, x)
  !   end if

  ! end subroutine hd_to_conserved

  ! !> Transform conservative variables into primitive ones
  ! subroutine hd_to_primitive(ixI^L, ixO^L, w, x)
  !   use mod_global_parameters
  !   use mod_dust, only: dust_to_primitive
  !   integer, intent(in)             :: ixI^L, ixO^L
  !   double precision, intent(inout) :: w(ixI^S, nw)
  !   double precision, intent(in)    :: x(ixI^S, 1:ndim)

  !   double precision                :: inv_rho
  !   integer :: ix^D

  !   if (fix_small_values) then
  !     call hd_handle_small_values(.false., w, x, ixI^L, ixO^L, 'hd_to_primitive')
  !   end if

  !  {do ix^DB=ixOmin^DB,ixOmax^DB\}
  !     inv_rho = 1.d0/w(ix^D,rho_)
  !     ! Convert momentum to velocity
  !     ^C&w(ix^D,m^C_)=w(ix^D,m^C_)*inv_rho\
  !     ! Calculate pressure = (gamma-1) * (e-ek)
  !     if(hd_energy) then
  !        ! Compute pressure
  !       w(ix^D,p_)=(hd_gamma-1.d0)*(w(ix^D,e_)&
  !                 -half*w(ix^D,rho_)*(^C&w(ix^D,m^C_)**2+))
  !     end if
  !  {end do\}

  !   ! Convert dust momentum to dust velocity
  !   if (hd_dust) then
  !     call dust_to_primitive(ixI^L, ixO^L, w, x)
  !   end if

  ! end subroutine hd_to_primitive

  !> Transform internal energy to total energy
  subroutine hd_ei_to_e(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)

    ! Calculate total energy from internal and kinetic energy
    w(ixO^S,e_)=w(ixO^S,e_)+half*(^C&w(ixO^S,m^C_)**2+)/w(ixO^S,rho_)

  end subroutine hd_ei_to_e

  !> Transform total energy to internal energy
  subroutine hd_e_to_ei(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)

    ! Calculate ei = e - ek
    w(ixO^S,e_)=w(ixO^S,e_)-half*(^C&w(ixO^S,m^C_)**2+)/w(ixO^S,rho_)

  end subroutine hd_e_to_ei

  !> Wrapper: e_to_ei + cache log10(nH) in wextra for LTE TC fast path.
  !> During STS substeps density is invariant, so log10(nH) is computed once
  !> per STS cycle (in sts_before_first_cycle hook) and reused across all substeps.
  subroutine hd_e_to_ei_and_cache_log_nH(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)

    call hd_e_to_ei(ixI^L,ixO^L,w,x)
    block%wextra(ixO^S, iw_log_nH) = dlog10(w(ixO^S, rho_) / eos%nH2rhoFactor)
  end subroutine hd_e_to_ei_and_cache_log_nH

  !> Calculate internal energy from total energy (non-modifying version)
  function hd_get_ei(w, ixI^L, ixO^L) result(ei)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: w(ixI^S, nw)
    double precision                :: ei(ixO^S)

    ! ei = e_total - e_kinetic
    ei(ixO^S) = w(ixO^S,e_) - half*(^C&w(ixO^S,m^C_)**2+)/w(ixO^S,rho_)
  end function hd_get_ei

  !> Calculate v_i = m_i / rho within ixO^L
  subroutine hd_get_v_idim(w, x, ixI^L, ixO^L, idim, v)
    use mod_global_parameters
    integer, intent(in)           :: ixI^L, ixO^L, idim
    double precision, intent(in)  :: w(ixI^S, nw), x(ixI^S, 1:ndim)
    double precision, intent(out) :: v(ixI^S)

    v(ixO^S) = w(ixO^S, mom(idim)) / w(ixO^S, rho_)
  end subroutine hd_get_v_idim

  !> Calculate velocity vector v_i = m_i / rho within ixO^L
  subroutine hd_get_v(w,x,ixI^L,ixO^L,v)
    use mod_global_parameters

    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S,nw), x(ixI^S,1:^ND)
    double precision, intent(out) :: v(ixI^S,1:ndir)

    integer :: idir

    do idir=1,ndir
      v(ixO^S,idir) = w(ixO^S, mom(idir)) / w(ixO^S, rho_)
    end do

  end subroutine hd_get_v

  !> Calculate cmax_idim = csound + abs(v_idim) within ixO^L
  subroutine hd_get_cmax(w, x, ixI^L, ixO^L, idim, cmax)
    use mod_global_parameters
    use mod_dust, only: dust_get_cmax_prim
    use mod_usr_methods, only: usr_set_pthermal

    integer, intent(in)                       :: ixI^L, ixO^L, idim
    ! w in primitive form
    double precision, intent(in)              :: w(ixI^S, nw), x(ixI^S, 1:ndim)
    double precision, intent(inout)           :: cmax(ixI^S)
    double precision                          :: csound2(ixI^S)

    if(hd_energy) then
      call eos%get_csound2(w, x, ixI^L, ixO^L, csound2)
      cmax(ixO^S)=dabs(w(ixO^S,mom(idim)))+dsqrt(csound2(ixO^S))
    else
      if (.not. associated(usr_set_pthermal)) then
        cmax(ixO^S) = hd_adiab * w(ixO^S, rho_)**eos%gamma
      else
        call usr_set_pthermal(w,x,ixI^L,ixO^L,cmax)
      end if
      cmax(ixO^S)=dabs(w(ixO^S,mom(idim)))+dsqrt(eos%gamma*cmax(ixO^S)/w(ixO^S,rho_))
    end if

    if (hd_dust) then
      call dust_get_cmax_prim(w, x, ixI^L, ixO^L, idim, cmax)
    end if
  end subroutine hd_get_cmax

  !> get adaptive cutoff temperature for TRAC (Johnston 2019 ApJL, 873, L22)
  subroutine hd_get_tcutoff(ixI^L,ixO^L,w,x,tco_local,Tmax_local)
    use mod_global_parameters
    use mod_usr_methods, only: usr_get_heating
    use mod_radiative_cooling, only: findL, calc_l_extended
    use mod_eos, only: eos
    use mod_geometry, only: gradient
    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    ! in primitive form
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(out) :: tco_local, Tmax_local

    double precision, parameter :: trac_delta=0.25d0
    double precision :: tmp1(ixI^S),Te(ixI^S),lts(ixI^S), R(ixI^S)
    double precision :: ltrc,ltrp
    integer :: jxO^L,hxO^L
    integer :: jxP^L,hxP^L,ixP^L
    logical :: lrlt(ixI^S)
    ! Johnston 2021 type 7 variables
    double precision :: dTdx, L_T, a_coeff, L1, cooling, net_cool
    double precision :: kappa_par, disc, kappa_TRAC, kappa_eff, Tcoff_eff
    double precision :: dx_over_delta, v_abs, v_n, gnorm, dl_eff
    double precision :: Q_heat(ixI^S), ne(ixI^S), nH_arr(ixI^S)
    double precision :: gradTd(ixI^S,1:ndim), nhat(1:ndim)
    integer :: ix^D, idims

    call eos%get_Rfactor(w,x,ixI^L,ixI^L,R)
    Te(ixI^S)=w(ixI^S,p_)/(R(ixI^S)*w(ixI^S,rho_))

    if (eos%eos_type == 'LTE') then
      Te(ixI^S) = w(ixI^S, Te_)
    endif

    Tco_local=zero
    Tmax_local=maxval(Te(ixO^S))
    select case(hd_trac_type)
    {^IFONED
    case(0)
      block%wextra(ixI^S,Tcoff_)=3.d5/unit_temperature
    case(1)
      hxO^L=ixO^L-1;
      jxO^L=ixO^L+1;
      lts(ixO^S)=0.5d0*dabs(Te(jxO^S)-Te(hxO^S))/Te(ixO^S)
      lrlt=.false.
      where(lts(ixO^S) > trac_delta)
        lrlt(ixO^S)=.true.
      end where
      if(any(lrlt(ixO^S))) then
        Tco_local=maxval(Te(ixO^S), mask=lrlt(ixO^S))
      end if
    case(2)
      !> iijima et al. 2021, LTRAC method
      ltrc=1.5d0
      ltrp=2.5d0
      ixP^L=ixO^L^LADD1;
      hxO^L=ixO^L-1;
      jxO^L=ixO^L+1;
      hxP^L=ixP^L-1;
      jxP^L=ixP^L+1;
      lts(ixP^S)=0.5d0*abs(Te(jxP^S)-Te(hxP^S))/Te(ixP^S)
      lts(ixP^S)=max(one, (exp(lts(ixP^S))/ltrc)**ltrp)
      ! Smoothed Tcoff for interior cells
      lts(ixO^S)=0.25d0*(lts(jxO^S)+two*lts(ixO^S)+lts(hxO^S))
      block%wextra(ixO^S,Tcoff_)=Te(ixO^S)*lts(ixO^S)**0.4d0
      ! Fill one ghost cell on each side with unsmoothed Tcoff.
      ! The thermal conduction routine reads Tcoff at ixO +/- 1;
      ! wextra ghost cells are NOT halo-communicated, so without this
      ! the conduction sees stale values and breaks symmetry.
      block%wextra(ixOmin1-1,Tcoff_)=Te(ixOmin1-1)*lts(ixOmin1-1)**0.4d0
      block%wextra(ixOmax1+1,Tcoff_)=Te(ixOmax1+1)*lts(ixOmax1+1)**0.4d0
    }
    case(7)
      {^IFONED
      !> Johnston et al. 2021 local TRAC (A&A 654, A2) -- 1D
      !> Per-cell kappa_TRAC from steady-state energy balance

      ! Get background heating Q (once per block)
      call usr_get_heating(Q_heat, ixI^L, ixO^L, w, x)

      ! Get n_e and n_H for cooling rate: Q = n_e * n_H * Lambda(T)
      ! For FI: n_e = n_H * neOnH_FI.  For LTE: n_e from Saha EoS.
      call eos%get_ne_nH(ixI^L, ixO^L, w, x, ne, nH_arr)

      hxO^L=ixO^L-1;
      jxO^L=ixO^L+1;
      dx_over_delta = dxlevel(1) / hd_trac_delta

      do ix1=ixOmin1,ixOmax1
        ! Temperature gradient: L_T = T / abs(dT/dx)
        dTdx = abs(Te(ix1+1) - Te(ix1-1)) / (2.d0 * dxlevel(1))
        if(dTdx < smalldouble) then
          ! Uniform temperature -- no broadening needed
          block%wextra(ix1,Tcoff_) = Te(ix1)
          cycle
        end if
        L_T = Te(ix1) / dTdx

        ! Eq.11 mass flux coefficient a = (5/2) p v / T, which is (5/2) kB J for
        ! an ideal gas since p = n kB T.
        v_abs = abs(w(ix1,m1_))
        a_coeff = 2.5d0 * w(ix1,p_) * v_abs / Te(ix1)

        ! Radiative cooling: n_e * n_H * Lambda(T), guarded to the table range (positive in-range
        ! test so a non-finite Te falls to zero rather than indexing findL with int(NaN)).
        if(Te(ix1) > rc_fl%tcoolmin .and. Te(ix1) < rc_fl%tcoolmax) then
          call findL(Te(ix1), L1, rc_fl); cooling = ne(ix1) * nH_arr(ix1) * L1
        else if(Te(ix1) >= rc_fl%tcoolmax) then
          call calc_l_extended(Te(ix1), L1, rc_fl); cooling = ne(ix1) * nH_arr(ix1) * L1
        else
          cooling = 0.d0
        end if

        ! Net cooling - heating
        net_cool = abs(cooling - Q_heat(ix1))

        ! Spitzer conductivity at this T
        kappa_par = tc_fl%tc_k_para * Te(ix1)**2.5d0

        ! Johnston Eq. 11 discriminant
        disc = a_coeff**2 + 4.d0 * tc_fl%tc_k_para * Te(ix1)**1.5d0 * net_cool

        if(L_T <= 2.d0 * dx_over_delta) then
          ! Under-resolved: full TRAC formula (Eq. 11)
          kappa_TRAC = (a_coeff + dsqrt(disc)) / (2.d0 / dx_over_delta)
        else
          ! Over-resolved: limiter only (Eq. 12, drops mass flux)
          kappa_TRAC = dsqrt(4.d0 * tc_fl%tc_k_para * Te(ix1)**1.5d0 * net_cool) &
                       / (2.d0 / dx_over_delta)
        end if

        ! kappa' = max(kappa_TRAC, kappa_par)
        kappa_eff = max(kappa_TRAC, kappa_par)

        ! Convert to effective Tcoff: Tcoff = (kappa'/kappa_0)**(2/5)
        Tcoff_eff = (kappa_eff / tc_fl%tc_k_para)**0.4d0

        ! Store max(Te, Tcoff_eff) -- Tcoff must be >= Te
        block%wextra(ix1,Tcoff_) = max(Te(ix1), Tcoff_eff)
      end do
      ! Fill one ghost cell on each side with nearest interior Tcoff.
      ! The thermal conduction routine reads Tcoff at ixO +/- 1;
      ! wextra ghost cells are NOT halo-communicated, so without this
      ! the conduction sees stale values and breaks symmetry.
      block%wextra(ixOmin1-1,Tcoff_) = block%wextra(ixOmin1,Tcoff_)
      block%wextra(ixOmax1+1,Tcoff_) = block%wextra(ixOmax1,Tcoff_)
      }
      {^NOONED
      ! Johnston et al. 2021 local TRAC, A and A 654 A2 -- multi-D, per-cell.
      ! Isotropic Spitzer conduction, so the TR normal is the unit temperature
      ! gradient direction nhat, and the broadening is applied along it.
      call usr_get_heating(Q_heat, ixI^L, ixO^L, w, x)
      call eos%get_ne_nH(ixI^L, ixO^L, w, x, ne, nH_arr)
      ! default (incl. non-communicated ghost layer): Tcoff=Te -> no broadening; interior overwritten
      block%wextra(ixI^S,Tcoff_) = Te(ixI^S)
      do idims=1,ndim
        call gradient(Te,ixI^L,ixO^L,idims,gradTd(ixI^S,idims))
      end do
      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        gnorm = dsqrt(^D&gradTd({ix^D},^D)**2+ )
        if(gnorm < smalldouble) then
          block%wextra(ix^D,Tcoff_) = Te(ix^D)            ! locally uniform T: no broadening
        else
          ^D&nhat(^D)=gradTd({ix^D},^D)/gnorm\
          ! grid spacing measured along the gradient direction
          dl_eff = 1.d0/dsqrt(^D&(nhat(^D)/block%ds({ix^D},^D))**2+ )
          L_T    = Te(ix^D)/gnorm
          ! Eq.16 enthalpy (mass) flux along the gradient
          v_n      = ^D&w({ix^D},mom(^D))*nhat(^D)+
          a_coeff  = 2.5d0*w(ix^D,p_)*dabs(v_n)/Te(ix^D)
          ! optically-thin cooling n_e n_H Lambda(T), guarded to the table range with a positive
          ! in-range test so a non-finite Te (NaN from an under-resolved collapse) falls to zero
          ! rather than indexing findL with int(NaN) -> out of bounds.
          if(Te(ix^D) > rc_fl%tcoolmin .and. Te(ix^D) < rc_fl%tcoolmax) then
            call findL(Te(ix^D),L1,rc_fl); cooling = L1*ne(ix^D)*nH_arr(ix^D)
          else if(Te(ix^D) >= rc_fl%tcoolmax) then
            call calc_l_extended(Te(ix^D),L1,rc_fl); cooling = L1*ne(ix^D)*nH_arr(ix^D)
          else
            cooling = 0.d0
          end if
          net_cool = dabs(cooling-Q_heat(ix^D))
          kappa_par     = tc_fl%tc_k_para*Te(ix^D)**2.5d0
          disc          = a_coeff**2 + 4.d0*tc_fl%tc_k_para*Te(ix^D)**1.5d0*net_cool
          dx_over_delta = dl_eff/hd_trac_delta
          if(L_T <= 2.d0*dx_over_delta) then
            kappa_TRAC = (a_coeff+dsqrt(disc))/(2.d0/dx_over_delta)
          else
            kappa_TRAC = dsqrt(4.d0*tc_fl%tc_k_para*Te(ix^D)**1.5d0*net_cool)/(2.d0/dx_over_delta)
          end if
          kappa_eff = max(kappa_TRAC,kappa_par)
          Tcoff_eff = (kappa_eff/tc_fl%tc_k_para)**0.4d0
          block%wextra(ix^D,Tcoff_) = max(Te(ix^D),Tcoff_eff)
        end if
      {end do\}
      }
    case default
      call mpistop("hd_trac_type not allowed")
    end select
  end subroutine hd_get_tcutoff

  !> Calculate cmax_idim = csound + abs(v_idim) within ixO^L
  subroutine hd_get_cbounds(wLC, wRC, wLp, wRp, x, ixI^L, ixO^L, idim,Hspeed,cmax, cmin)
    use mod_global_parameters
    use mod_dust, only: dust_get_cmax
    use mod_fld, only: fld_bound_diff
    use mod_variables

    integer, intent(in)             :: ixI^L, ixO^L, idim
    ! conservative left and right status
    double precision, intent(in)    :: wLC(ixI^S, nw), wRC(ixI^S, nw)
    ! primitive left and right status
    double precision, intent(in)    :: wLp(ixI^S, nw), wRp(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    double precision, intent(inout) :: cmax(ixI^S,1:number_species)
    double precision, intent(inout), optional :: cmin(ixI^S,1:number_species)
    double precision, intent(in)    :: Hspeed(ixI^S,1:number_species)

    double precision :: wmean(ixI^S,nw)
    double precision, dimension(ixI^S) :: umean, dmean, csoundL, csoundR, tmp1,tmp2,tmp3
    integer :: ix^D

    select case(boundspeed)
    case (1)
      ! This implements formula (10.52) from "Riemann Solvers and Numerical
      ! Methods for Fluid Dynamics" by Toro.

      tmp1(ixO^S)=dsqrt(wLp(ixO^S,rho_))
      tmp2(ixO^S)=dsqrt(wRp(ixO^S,rho_))
      tmp3(ixO^S)=1.d0/(tmp1(ixO^S)+tmp2(ixO^S))
      umean(ixO^S)=(wLp(ixO^S,mom(idim))*tmp1(ixO^S)+wRp(ixO^S,mom(idim))*tmp2(ixO^S))*tmp3(ixO^S)

      if(hd_energy) then
        call eos%get_csound2(wLp, x, ixI^L, ixO^L, csoundL)
        call eos%get_csound2(wRp, x, ixI^L, ixO^L, csoundR)
        if(hd_radiation_fld.and.fld_bound_diff)then
          csoundL(ixO^S)=csoundL(ixO^S)+(4.0d0/9.0d0)*wLp(ixO^S,r_e)/wLp(ixO^S,rho_)
          csoundR(ixO^S)=csoundR(ixO^S)+(4.0d0/9.0d0)*wRp(ixO^S,r_e)/wRp(ixO^S,rho_)
        endif
      else
         ! note usage of conservatives here
         call hd_get_csound2(wLC,x,ixI^L,ixO^L,csoundL)
         call hd_get_csound2(wRC,x,ixI^L,ixO^L,csoundR)
      end if

      dmean(ixO^S) = (tmp1(ixO^S)*csoundL(ixO^S)+tmp2(ixO^S)*csoundR(ixO^S)) * &
           tmp3(ixO^S) + 0.5d0*tmp1(ixO^S)*tmp2(ixO^S)*tmp3(ixO^S)**2 * &
           (wRp(ixO^S,mom(idim))-wLp(ixO^S,mom(idim)))**2

      dmean(ixO^S)=dsqrt(dmean(ixO^S))
      if(present(cmin)) then
        cmin(ixO^S,1)=umean(ixO^S)-dmean(ixO^S)
        cmax(ixO^S,1)=umean(ixO^S)+dmean(ixO^S)
        if(H_correction) then
          {do ix^DB=ixOmin^DB,ixOmax^DB\}
            cmin(ix^D,1)=sign(one,cmin(ix^D,1))*max(abs(cmin(ix^D,1)),Hspeed(ix^D,1))
            cmax(ix^D,1)=sign(one,cmax(ix^D,1))*max(abs(cmax(ix^D,1)),Hspeed(ix^D,1))
          {end do\}
        end if
      else
        cmax(ixO^S,1)=dabs(umean(ixO^S))+dmean(ixO^S)
      end if

      if (hd_dust) then
        wmean(ixO^S,1:nwflux)=0.5d0*(wLC(ixO^S,1:nwflux)+wRC(ixO^S,1:nwflux))
        call dust_get_cmax(wmean, x, ixI^L, ixO^L, idim, cmax, cmin)
      end if

    case (2)
      !if(hd_energy) then
      !   ! note usage of primitives here
      !   wmean(ixO^S,1:nwflux)=0.5d0*(wLp(ixO^S,1:nwflux)+wRp(ixO^S,1:nwflux))
      !   tmp1(ixO^S)=wmean(ixO^S,mom(idim))
      !   csoundR(ixO^S)=hd_gamma*wmean(ixO^S,p_)/wmean(ixO^S,rho_)
      !else
         ! note usage of conservatives here
         wmean(ixO^S,1:nwflux)=0.5d0*(wLC(ixO^S,1:nwflux)+wRC(ixO^S,1:nwflux))
         tmp1(ixO^S)=wmean(ixO^S,mom(idim))/wmean(ixO^S,rho_)
         call hd_get_csound2(wmean,x,ixI^L,ixO^L,csoundR)
         if(hd_radiation_fld.and.fld_bound_diff)csoundR(ixO^S)=csoundR(ixO^S)+(4.0d0/9.0d0)*wmean(ixO^S,r_e)/wmean(ixO^S,rho_)
      !endif
      csoundR(ixO^S) = dsqrt(csoundR(ixO^S))

      if(present(cmin)) then
        cmax(ixO^S,1)=max(tmp1(ixO^S)+csoundR(ixO^S),zero)
        cmin(ixO^S,1)=min(tmp1(ixO^S)-csoundR(ixO^S),zero)
        if(H_correction) then
          {do ix^DB=ixOmin^DB,ixOmax^DB\}
            cmin(ix^D,1)=sign(one,cmin(ix^D,1))*max(abs(cmin(ix^D,1)),Hspeed(ix^D,1))
            cmax(ix^D,1)=sign(one,cmax(ix^D,1))*max(abs(cmax(ix^D,1)),Hspeed(ix^D,1))
          {end do\}
        end if
      else
        cmax(ixO^S,1)=dabs(tmp1(ixO^S))+csoundR(ixO^S)
      end if

      if (hd_dust) then
        call dust_get_cmax(wmean, x, ixI^L, ixO^L, idim, cmax, cmin)
      end if
    case (3)
      ! Miyoshi 2005 JCP 208, 315 equation (67)
      if(hd_energy) then
        call eos%get_csound2(wLp, x, ixI^L, ixO^L, csoundL)
        call eos%get_csound2(wRp, x, ixI^L, ixO^L, csoundR)
        if(hd_radiation_fld.and.fld_bound_diff)then
          csoundL(ixO^S)=csoundL(ixO^S)+(4.0d0/9.0d0)*wLp(ixO^S,r_e)/wLp(ixO^S,rho_)
          csoundR(ixO^S)=csoundR(ixO^S)+(4.0d0/9.0d0)*wRp(ixO^S,r_e)/wRp(ixO^S,rho_)
        endif
      else
         ! note usage of conservatives here
         call hd_get_csound2(wLC,x,ixI^L,ixO^L,csoundL)
         call hd_get_csound2(wRC,x,ixI^L,ixO^L,csoundR)
      end if
      csoundL(ixO^S)=max(dsqrt(csoundL(ixO^S)),dsqrt(csoundR(ixO^S)))
      if(present(cmin)) then
        cmin(ixO^S,1)=min(wLp(ixO^S,mom(idim)),wRp(ixO^S,mom(idim)))-csoundL(ixO^S)
        cmax(ixO^S,1)=max(wLp(ixO^S,mom(idim)),wRp(ixO^S,mom(idim)))+csoundL(ixO^S)
        if(H_correction) then
          {do ix^DB=ixOmin^DB,ixOmax^DB\}
            cmin(ix^D,1)=sign(one,cmin(ix^D,1))*max(abs(cmin(ix^D,1)),Hspeed(ix^D,1))
            cmax(ix^D,1)=sign(one,cmax(ix^D,1))*max(abs(cmax(ix^D,1)),Hspeed(ix^D,1))
          {end do\}
        end if
      else
        cmax(ixO^S,1)=max(wLp(ixO^S,mom(idim)),wRp(ixO^S,mom(idim)))+csoundL(ixO^S)
      end if
      if (hd_dust) then
        wmean(ixO^S,1:nwflux)=0.5d0*(wLC(ixO^S,1:nwflux)+wRC(ixO^S,1:nwflux))
        call dust_get_cmax(wmean, x, ixI^L, ixO^L, idim, cmax, cmin)
      end if
    case (4)
      !> PVRS pressure-based wave speed estimate (Toro 1999, Section 10.5.2)
      !> Recommended by Coleman 2020 for general EoS given limitations of constant gamma approximation.
      !> Estimates star pressure from linearised Riemann problem, then uses
      !> RH to detect shock vs rarefaction on each side.

      if(hd_energy) then
        call eos%get_csound2(wLp, x, ixI^L, ixO^L, csoundL)
        call eos%get_csound2(wRp, x, ixI^L, ixO^L, csoundR)
        if(hd_radiation_fld.and.fld_bound_diff)then
          csoundL(ixO^S)=csoundL(ixO^S)+(4.0d0/9.0d0)*wLp(ixO^S,r_e)/wLp(ixO^S,rho_)
          csoundR(ixO^S)=csoundR(ixO^S)+(4.0d0/9.0d0)*wRp(ixO^S,r_e)/wRp(ixO^S,rho_)
        endif
      else
        call hd_get_csound2(wLC,x,ixI^L,ixO^L,csoundL)
        call hd_get_csound2(wRC,x,ixI^L,ixO^L,csoundR)
      end if
      csoundL(ixO^S) = dsqrt(csoundL(ixO^S))
      csoundR(ixO^S) = dsqrt(csoundR(ixO^S))
      if(present(cmin)) then
        {do ix^DB=ixOmin^DB,ixOmax^DB\}
          !> PVRS star pressure estimate (Toro Eq. 9.28)
          !> cup = rho_bar * a_bar
          tmp1(ix^D) = 0.25d0*(wLp(ix^D,rho_)+wRp(ix^D,rho_)) &
                      *(csoundL(ix^D)+csoundR(ix^D))
          !> p* = max(0, p_avg + 0.5*(u_L - u_R)*cup)
          tmp2(ix^D) = max(zero, 0.5d0*(wLp(ix^D,e_)+wRp(ix^D,e_)) &
                     + 0.5d0*(wLp(ix^D,mom(idim))-wRp(ix^D,mom(idim))) &
                     *tmp1(ix^D))
          !> Left wave speed: S_L = u_L - a_L * q_L
          if(tmp2(ix^D) > wLp(ix^D,e_) .and. wLp(ix^D,e_) > zero) then
            !> Left shock: q_L from R-H with local Gamma1 = a^2*rho/p
            tmp3(ix^D) = csoundL(ix^D)**2*wLp(ix^D,rho_)/wLp(ix^D,e_)
            dmean(ix^D) = dsqrt(1.0d0 + (tmp3(ix^D)+1.0d0) &
                        /(2.0d0*tmp3(ix^D)) &
                        *(tmp2(ix^D)/wLp(ix^D,e_) - 1.0d0))
          else
            !> Left rarefaction
            dmean(ix^D) = 1.0d0
          end if
          cmin(ix^D,1) = wLp(ix^D,mom(idim)) - csoundL(ix^D)*dmean(ix^D)
          !> Right wave speed: S_R = u_R + a_R * q_R
          if(tmp2(ix^D) > wRp(ix^D,e_) .and. wRp(ix^D,e_) > zero) then
            !> Right shock: q_R from R-H with local Gamma1
            tmp3(ix^D) = csoundR(ix^D)**2*wRp(ix^D,rho_)/wRp(ix^D,e_)
            dmean(ix^D) = dsqrt(1.0d0 + (tmp3(ix^D)+1.0d0) &
                        /(2.0d0*tmp3(ix^D)) &
                        *(tmp2(ix^D)/wRp(ix^D,e_) - 1.0d0))
          else
            !> Right rarefaction
            dmean(ix^D) = 1.0d0
          end if
          cmax(ix^D,1) = wRp(ix^D,mom(idim)) + csoundR(ix^D)*dmean(ix^D)
        {end do\}
        if(H_correction) then
          {do ix^DB=ixOmin^DB,ixOmax^DB\}
            cmin(ix^D,1)=sign(one,cmin(ix^D,1))*max(abs(cmin(ix^D,1)),Hspeed(ix^D,1))
            cmax(ix^D,1)=sign(one,cmax(ix^D,1))*max(abs(cmax(ix^D,1)),Hspeed(ix^D,1))
          {end do\}
        end if
      else
        {do ix^DB=ixOmin^DB,ixOmax^DB\}
          tmp1(ix^D) = 0.25d0*(wLp(ix^D,rho_)+wRp(ix^D,rho_)) &
                      *(csoundL(ix^D)+csoundR(ix^D))
          tmp2(ix^D) = max(zero, 0.5d0*(wLp(ix^D,e_)+wRp(ix^D,e_)) &
                     + 0.5d0*(wLp(ix^D,mom(idim))-wRp(ix^D,mom(idim))) &
                     *tmp1(ix^D))
          if(tmp2(ix^D) > wLp(ix^D,e_) .and. wLp(ix^D,e_) > zero) then
            tmp3(ix^D) = csoundL(ix^D)**2*wLp(ix^D,rho_)/wLp(ix^D,e_)
            dmean(ix^D) = dsqrt(1.0d0 + (tmp3(ix^D)+1.0d0) &
                        /(2.0d0*tmp3(ix^D)) &
                        *(tmp2(ix^D)/wLp(ix^D,e_) - 1.0d0))
          else
            dmean(ix^D) = 1.0d0
          end if
          umean(ix^D) = dabs(wLp(ix^D,mom(idim)) &
                      - csoundL(ix^D)*dmean(ix^D))
          if(tmp2(ix^D) > wRp(ix^D,e_) .and. wRp(ix^D,e_) > zero) then
            tmp3(ix^D) = csoundR(ix^D)**2*wRp(ix^D,rho_)/wRp(ix^D,e_)
            dmean(ix^D) = dsqrt(1.0d0 + (tmp3(ix^D)+1.0d0) &
                        /(2.0d0*tmp3(ix^D)) &
                        *(tmp2(ix^D)/wRp(ix^D,e_) - 1.0d0))
          else
            dmean(ix^D) = 1.0d0
          end if
          cmax(ix^D,1) = max(umean(ix^D), &
                        wRp(ix^D,mom(idim))+csoundR(ix^D)*dmean(ix^D))
        {end do\}
      end if
      if (hd_dust) then
        wmean(ixO^S,1:nwflux)=0.5d0*(wLC(ixO^S,1:nwflux)+wRC(ixO^S,1:nwflux))
        call dust_get_cmax(wmean, x, ixI^L, ixO^L, idim, cmax, cmin)
      end if
    end select

  end subroutine hd_get_cbounds

  !> Calculate the square of the thermal sound speed csound2 within ixO^L.
  !> For conserved w: extracts pthermal first, then applies Gamma_1.
  !> For LTE+IonE: look up Gamma_1 from pressure-indexed table, then cs2 = Gamma_1 * p/rho.
  !> Uses pressure-indexed table (gamma1_from_nH_p) because pressure is continuous
  !> at contact discontinuities, avoiding spurious gamma1 spikes from the eint-indexed table.
  subroutine hd_get_csound2(w,x,ixI^L,ixO^L,csound2)
    use mod_global_parameters
    use mod_timing
    use mod_eos_LTE, only: gamma1_from_nH_p
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: w(ixI^S,nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out)   :: csound2(ixI^S)
    double precision :: pthermal(ixI^S)
    double precision :: nH_val, g1
    double precision :: local_t0
    integer :: ix^D

    !> get_thermal_pressure has its own timing; only time the gamma1 loop here
    call eos%get_thermal_pressure(w, x, ixI^L, ixO^L, pthermal)

    if (eos%ionE) then
      local_t0 = MPI_WTIME()
      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        nH_val = w(ix^D,rho_) / eos%nH2rhoFactor
        g1 = gamma1_from_nH_p(dlog10(nH_val), dlog10(pthermal(ix^D)/nH_val))
        csound2(ix^D) = g1 * pthermal(ix^D) / w(ix^D,rho_)
      {end do\}
      timeeos_csound=timeeos_csound+(MPI_WTIME()-local_t0)
    else
      csound2(ixO^S) = eos%gamma * pthermal(ixO^S) / w(ixO^S,rho_)
    end if

  end subroutine hd_get_csound2


  !> Calculate modified squared sound speed for FLD
  !> NOTE: only for diagnostic purposes, unused subroutine
  subroutine hd_get_csrad2(w,x,ixI^L,ixO^L,csound)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, nw), x(ixI^S,1:ndim)
    double precision, intent(out):: csound(ixI^S)

    double precision :: wprim(ixI^S, nw)

    wprim(ixI^S,1:nw)=w(ixI^S,1:nw)
    call eos%to_primitive(ixI^L,ixO^L,wprim,x)
    call hd_get_csrad2_prim(wprim,x,ixI^L,ixO^L,csound)

  end subroutine hd_get_csrad2

  !> Calculate modified squared sound speed for FLD
  !> NOTE: w is primitive on entry here!
  !> NOTE: used in FLD module as phys_get_csrad2
  subroutine hd_get_csrad2_prim(w,x,ixI^L,ixO^L,csound)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, nw), x(ixI^S,1:ndim)
    double precision, intent(out):: csound(ixI^S)

    integer :: ix^D
    double precision :: inv_rho
    double precision :: prad_tensor(ixI^S, 1:ndim, 1:ndim)
    double precision :: prad_max(ixI^S)
    integer :: idim

   if(hd_fld_pradtensor) then
       call hd_get_pradiation_from_prim(w, x, ixI^L, ixO^L, prad_tensor)
    else
       prad_tensor=zero 
       do idim=1,ndim
         prad_tensor(ixO^S,idim,idim)=w(ixO^S,r_e)/3.0d0
       enddo
    endif

   {do ix^DB=ixOmin^DB,ixOmax^DB \}
      inv_rho=1.d0/w(ix^D,rho_)
      prad_max(ix^D) = (4.0d0/3.0d0)*maxval(prad_tensor(ix^D,:,:))
      csound(ix^D)=(eos%gamma*w(ix^D,p_)+prad_max(ix^D))*inv_rho
   {end do\}

   if(minval(csound(ixO^S))<smalldouble)then
     print *,'issue with squared speed and rad pressure'
     print *,minval(csound(ixO^S))
     print *,minval(prad_max(ixO^S))
     call mpistop("negative squared speed in get_csrad2 for dt")
   endif

  end subroutine hd_get_csrad2_prim

  !> Calculate radiation pressure within ixO^L
  !> NOTE: w is primitive on entry here!
  !> NOTE: used in FLD module as it is called from phys_get_csrad2
  subroutine hd_get_pradiation_from_prim(w, x, ixI^L, ixO^L, prad)
    use mod_global_parameters
    use mod_fld
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: prad(ixI^S, 1:ndim, 1:ndim)

    call fld_get_radpress(w, x, ixI^L, ixO^L, prad, fld_fl)

  end subroutine hd_get_pradiation_from_prim

  !> calculates the sum of the gas pressure and max Prad tensor element
  !> NOTE: only for diagnostic purposes, unused subroutine
  subroutine hd_get_pthermal_plus_pradiation(w, x, ixI^L, ixO^L, pth_plus_prad)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: pth_plus_prad(ixI^S)

    double precision             :: wprim(ixI^S, 1:nw)
    double precision             :: prad_tensor(ixI^S, 1:ndim, 1:ndim)
    double precision             :: prad_max(ixI^S)
    integer :: ix^D

    wprim(ixI^S,1:nw)=w(ixI^S,1:nw)
    call eos%to_primitive(ixI^L,ixO^L,wprim,x)
    call hd_get_pradiation_from_prim(wprim, x, ixI^L, ixO^L, prad_tensor)
    {do ix^D = ixOmin^D,ixOmax^D\}
      prad_max(ix^D) = maxval(prad_tensor(ix^D,:,:))
    {enddo\}
    pth_plus_prad(ixO^S) = wprim(ixO^S,p_) + prad_max(ixO^S)
  end subroutine hd_get_pthermal_plus_pradiation

  !> Calculates radiation temperature
  subroutine hd_get_trad(w, x, ixI^L, ixO^L, trad)
    use mod_global_parameters
    use mod_constants

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: trad(ixI^S)

    trad(ixI^S) = (w(ixI^S,r_e)/arad_norm)**(1.d0/4.d0)

  end subroutine hd_get_trad

  !> Calculate temperature=p/rho when in e_ the  total energy is stored
  subroutine hd_get_temperature_from_etot(w, x, ixI^L, ixO^L, res)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: res(ixI^S)

    double precision :: R(ixI^S)

    call eos%get_Rfactor(w,x,ixI^L,ixO^L,R)
    call eos%get_thermal_pressure(w, x, ixI^L, ixO^L, res)
    res(ixO^S)=res(ixO^S)/(R(ixO^S)*w(ixO^S,rho_))
  end subroutine hd_get_temperature_from_etot

  !> Calculate temperature=p/rho when in e_ the  internal energy is stored
  subroutine hd_get_temperature_from_eint(w, x, ixI^L, ixO^L, res)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: res(ixI^S)

    double precision :: R(ixI^S)

    call eos%get_Rfactor(w,x,ixI^L,ixO^L,R)
    res(ixO^S) = (eos%gamma - 1.0d0) * w(ixO^S, e_)/(w(ixO^S,rho_)*R(ixO^S))
  end subroutine hd_get_temperature_from_eint

  ! Calculate flux f_idim[iw]
  subroutine hd_get_flux(wC, w, x, ixI^L, ixO^L, idim, f)
    use mod_global_parameters
    use mod_dust, only: dust_get_flux_prim

    integer, intent(in)             :: ixI^L, ixO^L, idim
    ! conservative w
    double precision, intent(in)    :: wC(ixI^S, 1:nw)
    ! primitive w
    double precision, intent(in)    :: w(ixI^S, 1:nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    double precision, intent(out)   :: f(ixI^S, nwflux)

    double precision                :: pth(ixI^S)
    integer                         :: ix^DB

    if (hd_energy) then
     {do ix^DB=ixOmin^DB,ixOmax^DB\}
        f(ix^D,rho_)=w(ix^D,mom(idim))*w(ix^D,rho_)
        ! Momentum flux is v_i*m_i, +p in direction idim
        ^C&f(ix^D,m^C_)=w(ix^D,mom(idim))*wC(ix^D,m^C_)\
        f(ix^D,mom(idim))=f(ix^D,mom(idim))+w(ix^D,p_)
        ! Energy flux is v_i*(e + p)
        f(ix^D,e_)=w(ix^D,mom(idim))*(wC(ix^D,e_)+w(ix^D,p_))
     {end do\}

     ! FACE-RECIPE: q is NOT contributing to the Riemann energy flux.
     ! The conductive heat-flux contribution is added by the face-
     ! recipe in add_hypertc_source as a post-Riemann sweep, computed
     ! from cell-centred Te (read from the cached Te_ field, which
     ! update_eos sets correctly each substep using the proper
     ! internal-energy/EoS path). q's own advective flux is zero.
     if(hd_hyperbolic_thermal_conduction) then
      {do ix^DB=ixOmin^DB,ixOmax^DB\}
         f(ix^D,q_)=zero
      {end do\}
     end if
    else
      call eos%get_thermal_pressure(wC, x, ixI^L, ixO^L, pth)
     {do ix^DB=ixOmin^DB,ixOmax^DB\}
        f(ix^D,rho_)=w(ix^D,mom(idim))*w(ix^D,rho_)
        ! Momentum flux is v_i*m_i, +p in direction idim
        ^C&f(ix^D,m^C_)=w(ix^D,mom(idim))*wC(ix^D,m^C_)\
        f(ix^D,mom(idim))=f(ix^D,mom(idim))+pth(ix^D)
     {end do\}
    end if

    if(hd_radiation_fld)then
     {do ix^DB=ixOmin^DB,ixOmax^DB\}
        ! advection of radiation enery v_i*r_e
        f(ix^D,r_e)=w(ix^D,mom(idim))*wC(ix^D,r_e)
     {end do\}
    endif

    do ix1 = 1, hd_n_tracer
       f(ixO^S, tracer(ix1)) = w(ixO^S,mom(idim)) * w(ixO^S, tracer(ix1))
    end do

    if (hd_fip) then
      f(ixO^S,fip_) = w(ixO^S,mom(idim)) * wC(ixO^S,fip_)
    end if

    ! Dust fluxes
    if (hd_dust) then
      call dust_get_flux_prim(w, x, ixI^L, ixO^L, idim, f)
    end if

  end subroutine hd_get_flux

  !> Add geometrical source terms to w
  !>
  !> Notice that the expressions of the geometrical terms depend only on ndir,
  !> not ndim. Eg, they are the same in 2.5D and in 3D, for any geometry.
  !>
  subroutine hd_add_source_geom(qdt, dtfactor, ixI^L, ixO^L, wCT, wprim, w, x)
    use mod_global_parameters
    use mod_usr_methods, only: usr_set_surface, usr_set_pthermal
    use mod_rotating_frame, only: rotating_frame_add_source
    use mod_dust, only: dust_n_species, dust_mom, dust_rho
    use mod_geometry
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, dtfactor, x(ixI^S, 1:ndim)
    double precision, intent(inout) :: wCT(ixI^S, 1:nw), wprim(ixI^S,1:nw),w(ixI^S, 1:nw)
    double precision :: pth(ixI^S), source(ixI^S), minrho
    integer                         :: iw,idir, h1x^L{^NOONED, h2x^L}
    integer :: mr_,mphi_ ! Polar var. names
    integer :: irho, ifluid, n_fluids
    double precision :: exp_factor(ixI^S), del_exp_factor(ixI^S), exp_factor_primitive(ixI^S)

    if (hd_dust) then
       n_fluids = 1 + dust_n_species
    else
       n_fluids = 1
    end if

    select case (coordinate)

    case(Cartesian_expansion)
      !the user provides the functions of exp_factor and del_exp_factor
      if(associated(usr_set_surface)) call usr_set_surface(ixI^L,x,block%dx,exp_factor,del_exp_factor,exp_factor_primitive)
      if(hd_energy) then
        source(ixO^S)=wprim(ixO^S, p_)
      else
        if(.not. associated(usr_set_pthermal)) then
          source(ixO^S)=hd_adiab * wprim(ixO^S, rho_)**eos%gamma
        else
          call usr_set_pthermal(wCT,x,ixI^L,ixO^L,source)
        end if
      end if
      source(ixO^S) = source(ixO^S)*del_exp_factor(ixO^S)/exp_factor(ixO^S)
      w(ixO^S,mom(1)) = w(ixO^S,mom(1)) + qdt*source(ixO^S)

    case (cylindrical)
       do ifluid = 0, n_fluids-1
          ! s[mr]=(pthermal+mphi**2/rho)/radius
          if (ifluid == 0) then
            ! gas
            irho  = rho_
            mr_   = mom(r_)
            if(phi_>0) mphi_ = mom(phi_)
            if(hd_energy) then
              source(ixO^S)=wprim(ixO^S, p_)
            else
              if(.not. associated(usr_set_pthermal)) then
                source(ixO^S)=hd_adiab * wprim(ixO^S, rho_)**eos%gamma
              else
                call usr_set_pthermal(wCT,x,ixI^L,ixO^L,source)
              end if
            end if
            minrho = 0.0d0
          else
            ! dust : no pressure
            irho  = dust_rho(ifluid)
            mr_   = dust_mom(r_, ifluid)
            if(phi_>0) mphi_ = dust_mom(phi_, ifluid)
            source(ixI^S) = zero
            minrho = 0.0d0
          end if
          if(phi_ > 0) then
            where (wCT(ixO^S, irho) > minrho)
              source(ixO^S) = source(ixO^S) + wCT(ixO^S,mphi_)*wprim(ixO^S,mphi_)
              w(ixO^S, mr_) = w(ixO^S, mr_) + qdt*source(ixO^S)/x(ixO^S,r_)
            end where
            ! s[mphi]=(-mphi*vr)/radius
            where (wCT(ixO^S, irho) > minrho)
              source(ixO^S) = -wCT(ixO^S, mphi_) * wprim(ixO^S, mr_)
              w(ixO^S, mphi_) = w(ixO^S, mphi_) + qdt * source(ixO^S) / x(ixO^S, r_)
            end where
          else
            ! s[mr]=2pthermal/radius
            w(ixO^S, mr_) = w(ixO^S, mr_) + qdt * source(ixO^S) / x(ixO^S, r_)
          end if
       end do
    case (spherical)
       if (hd_dust) then
          call mpistop("Dust geom source terms not implemented yet with spherical geometries")
       end if
       mr_   = mom(r_)
       if(phi_>0) mphi_ = mom(phi_)
       h1x^L=ixO^L-kr(1,^D); {^NOONED h2x^L=ixO^L-kr(2,^D);}
       if(hd_energy) then
         pth(ixO^S)=wprim(ixO^S, p_)
       else
         if(.not. associated(usr_set_pthermal)) then
           pth(ixO^S)=hd_adiab * wprim(ixO^S, rho_)**eos%gamma
         else
           call usr_set_pthermal(wCT,x,ixI^L,ixO^L,pth)
         end if
       end if
       ! s[mr]=((vtheta**2+vphi**2)*rho+2*p)/r
       source(ixO^S) = pth(ixO^S) * x(ixO^S, 1) &
            *(block%surfaceC(ixO^S, 1) - block%surfaceC(h1x^S, 1)) &
            /block%dvolume(ixO^S)
       do idir = 2, ndir
         source(ixO^S) = source(ixO^S) + wprim(ixO^S, mom(idir))**2 * wprim(ixO^S, rho_)
       end do
       w(ixO^S, mr_) = w(ixO^S, mr_) + qdt * source(ixO^S) / x(ixO^S, 1)

       {^NOONED
       ! s[mtheta]=-(vr*vtheta*rho)/r+cot(theta)*(vphi**2*rho+p)/r
       source(ixO^S) = pth(ixO^S) * x(ixO^S, 1) &
            * (block%surfaceC(ixO^S, 2) - block%surfaceC(h2x^S, 2)) &
            / block%dvolume(ixO^S)
       if (ndir == 3) then
         source(ixO^S) = source(ixO^S) + (wprim(ixO^S, mom(3))**2 * wprim(ixO^S, rho_)) / tan(x(ixO^S, 2))
       end if
       source(ixO^S) = source(ixO^S) - (wprim(ixO^S, mom(2)) * wprim(ixO^S, mr_)) * wprim(ixO^S, rho_)
       w(ixO^S, mom(2)) = w(ixO^S, mom(2)) + qdt * source(ixO^S) / x(ixO^S, 1)

       if (ndir == 3) then
         ! s[mphi]=-(vphi*vr/rho)/r-cot(theta)*(vtheta*vphi/rho)/r
         source(ixO^S) = -(wprim(ixO^S, mom(3)) * wprim(ixO^S, mr_)) * wprim(ixO^S, rho_)&
                        - (wprim(ixO^S, mom(2)) * wprim(ixO^S, mom(3))) * wprim(ixO^S, rho_) / tan(x(ixO^S, 2))
         w(ixO^S, mom(3)) = w(ixO^S, mom(3)) + qdt * source(ixO^S) / x(ixO^S, 1)
       end if
       }
    end select

    if (hd_rotating_frame) then
       if (hd_dust) then
          call mpistop("Rotating frame not implemented yet with dust")
       else
          call rotating_frame_add_source(qdt,dtfactor,ixI^L,ixO^L,wprim,w,x)
       end if
    end if

  end subroutine hd_add_source_geom

  ! w[iw]= w[iw]+qdt*S[wCT, qtC, x] where S is the source based on wCT within ixO
  subroutine hd_add_source(qdt,dtfactor, ixI^L,ixO^L,wCT,wCTprim,w,x,qsourcesplit,active)
    use mod_global_parameters
    use mod_radiative_cooling, only: radiative_cooling_add_source
    use mod_dust, only: dust_add_source, dust_mom, dust_rho, dust_n_species
    use mod_viscosity, only: viscosity_add_source
    use mod_usr_methods, only: usr_gravity
    use mod_gravity, only: gravity_add_source, grav_split
    use mod_cak_force, only: cak_add_source

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, dtfactor
    double precision, intent(in)    :: wCT(ixI^S, 1:nw),wCTprim(ixI^S,1:nw), x(ixI^S, 1:ndim)
    double precision, intent(inout) :: w(ixI^S, 1:nw)
    logical, intent(in)             :: qsourcesplit
    logical, intent(inout)          :: active

    double precision :: gravity_field(ixI^S, 1:ndim)
    integer :: idust, idim, ix^D

    if(hd_dust .and. .not. hd_dust_implicit) then
      call dust_add_source(qdt,ixI^L,ixO^L,wCT,w,x,qsourcesplit,active)
    end if

    
    ! if (mype == 0) then
    !   {do ix^DB = ixI^LIM^DB\}
    !     if (abs(x(ix^D,1) - xprobmax1/2.0d0) > (xprobmax1/2.0d0 - 1.0d8/unit_length)) then
    !       write(*,*) x(ix^D,1), ' ' , wCT(ix^D,e_)
    !     endif
    !   {end do\}
    ! endif
    if(hd_radiative_cooling) then
      call radiative_cooling_add_source(qdt,ixI^L,ixO^L,wCT,wCTprim,w,x,&
           qsourcesplit,active, rc_fl)
    end if

    if(hd_viscosity) then
      call viscosity_add_source(qdt,ixI^L,ixO^L,wCT,wCTprim,w,x,&
           hd_energy,qsourcesplit,active)
    end if

    if (hd_gravity) then
      call gravity_add_source(qdt,ixI^L,ixO^L,wCT,wCTprim,w,x,&
           hd_energy,qsourcesplit,active)

      if (hd_dust .and. qsourcesplit .eqv. grav_split) then
         active = .true.

         call usr_gravity(ixI^L, ixO^L, wCT, x, gravity_field)
         do idust = 1, dust_n_species
            do idim = 1, ndim
               w(ixO^S, dust_mom(idim, idust)) = w(ixO^S, dust_mom(idim, idust)) &
                    + qdt * gravity_field(ixO^S, idim) * wCT(ixO^S, dust_rho(idust))
            end do
         end do
      end if
    end if

    if (hd_cak_force) then
      call cak_add_source(qdt,ixI^L,ixO^L,wCT,w,x,hd_energy,qsourcesplit,active)
    end if

    ! This is where the radiation force and heating/cooling are added
    if (hd_radiation_fld) then
       call hd_add_radiation_source(qdt,ixI^L,ixO^L,wCT,wCTprim,w,x,qsourcesplit,active)
    endif

    if(eos%eos_type == 'PI') then
      if(.not.qsourcesplit) then
        active = .true.
        call eos%update_eos(ixI^L,ixO^L,w,x)
      end if
    end if

    ! Hyperbolic TC: cell-centred Cattaneo relaxation source for q.
    if(hd_hyperbolic_thermal_conduction .and. .not. qsourcesplit) then
      active = .true.
      call add_hypertc_source(qdt,ixI^L,ixO^L,wCT,w,x,wCTprim)
    end if

  end subroutine hd_add_source

  !> FACE-RECIPE for hyperbolic thermal conduction.
  !>
  !> Replaces the OLD cell-centred Cattaneo source. Architecture:
  !>  1. q is removed from the Riemann pipeline (f(e_) does NOT get +q
  !>     in hd_get_flux); q has no advective flux of its own.
  !>  2. Each face i+1/2 computes a Cattaneo-relaxed face heat flux
  !>     q_f from the local centred T gradient. Te is read from the
  !>     cached Te_ field (set by update_eos using the correct
  !>     internal-energy/EoS path; calling get_temperature_from_eint
  !>     directly on the conservative wCT inflates Te by KE/(R*rho)
  !>     at fast-flow cells, hence the cache read).
  !>  3. Energy update: E_i -= qdt * (q_{i+1/2} - q_{i-1/2}) / dx_i.
  !>  4. Cell-centred q refreshed from the face values for next step.
  !>
  !> Why face-recipe vs HLL-on-q: avoids reconstruction (Koren limiter
  !> applied to q at faces) and HLL wave-speed bias (q weighted by
  !> hydro signal speeds during HLL combination, which has no physical
  !> justification for a scalar conductive flux).
  !>
  !> 1D only.
  subroutine add_hypertc_source(qdt,ixI^L,ixO^L,wCT,w,x,wCTprim)
    use mod_global_parameters

    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: qdt
    double precision, dimension(ixI^S,1:ndim), intent(in) :: x
    double precision, dimension(ixI^S,1:nw), intent(in) :: wCT,wCTprim
    double precision, dimension(ixI^S,1:nw), intent(inout) :: w

    double precision :: Te(ixI^S), Rfactor(ixI^S)
    double precision :: qf_half(ixGlo1:ixGhi1)
    double precision :: qf_full(ixGlo1:ixGhi1)
    double precision :: T_L, T_R, T_cool_L, T_cool_R, T_kap_L, T_kap_R
    double precision :: kappa_L, kappa_R, kappa_f
    double precision :: sigma_T7_f, dT_face, dx_f
    double precision :: rho_f, p_L, p_R, p_f, cs_f, q_sat, q_Sp, q_Sp_star
    double precision :: f_sat, eint_loc_L, eint_loc_R, eint_loc_f
    double precision :: qn_face, tau_f, ratio, decay, ave_factor
    double precision :: q_face_max
    double precision, parameter :: ratio_cap = 50.0d0
    integer :: ix1, nface_lo, nface_hi

    ! Te via cached field (correct EoS path) or primitives (FI fallback).
    if (Te_ > 0) then
      Te(ixI^S) = wCT(ixI^S, Te_)
    else
      call eos%get_Rfactor(wCTprim, x, ixI^L, ixI^L, Rfactor)
      Te(ixI^S) = wCTprim(ixI^S, p_) / (Rfactor(ixI^S) * wCTprim(ixI^S, rho_))
    end if

    {^IFONED
    nface_lo = ixOmin1 - 1
    nface_hi = ixOmax1
    qf_half(:) = 0.0d0
    qf_full(:) = 0.0d0

    do ix1 = nface_lo, nface_hi
      ! Cell-centred T on either side of face, with TRAC clip
      T_L = Te(ix1)
      T_R = Te(ix1 + 1)
      if (hd_trac) then
        T_cool_L = block%wextra(ix1,     Tcoff_)
        T_cool_R = block%wextra(ix1 + 1, Tcoff_)
      else
        T_cool_L = 0.0d0
        T_cool_R = 0.0d0
      end if
      T_kap_L = max(T_L, T_cool_L)
      T_kap_R = max(T_R, T_cool_R)
      kappa_L = hypertc_kappa * dsqrt(T_kap_L**5)
      kappa_R = hypertc_kappa * dsqrt(T_kap_R**5)

      ! Harmonic-mean face kappa: limits to cold side at sharp jumps,
      ! preserves flux conservation across kappa discontinuities.
      if (kappa_L + kappa_R > smalldouble) then
        kappa_f = 2.0d0 * kappa_L * kappa_R / (kappa_L + kappa_R)
      else
        kappa_f = 0.0d0
      end if

      dx_f    = 0.5d0 * (block%ds(ix1, 1) + block%ds(ix1 + 1, 1))
      dT_face = T_R - T_L
      q_Sp    = -kappa_f * dT_face / dx_f

      ! Face-averaged primitives for saturation
      rho_f = 0.5d0 * (wCTprim(ix1, rho_) + wCTprim(ix1 + 1, rho_))
      p_L   = wCTprim(ix1,     p_)
      p_R   = wCTprim(ix1 + 1, p_)
      p_f   = 0.5d0 * (p_L + p_R)
      cs_f  = dsqrt(max(smalldouble, eos%gamma * p_f / max(rho_f, smalldouble)))

      ! Cowie-McKee saturation cap on the Spitzer target
      if (hd_htc_sat) then
        q_sat = 1.5d0 * rho_f * (p_f / max(rho_f, smalldouble))**1.5d0
        if (q_sat > smalldouble) then
          f_sat = 1.0d0 / (1.0d0 + dabs(q_Sp) / q_sat)
        else
          f_sat = 1.0d0
        end if
        q_Sp_star = f_sat * q_Sp
      else
        q_Sp_star = q_Sp
      end if

      ! Internal energy from conservatives (matches eint fix)
      eint_loc_L = wCT(ix1,     e_) &
                 - 0.5d0 * wCT(ix1,     m1_)**2 / max(wCT(ix1,     rho_), smalldouble)
      eint_loc_R = wCT(ix1 + 1, e_) &
                 - 0.5d0 * wCT(ix1 + 1, m1_)**2 / max(wCT(ix1 + 1, rho_), smalldouble)
      eint_loc_L = max(eint_loc_L, smalldouble)
      eint_loc_R = max(eint_loc_R, smalldouble)
      eint_loc_f = 0.5d0 * (eint_loc_L + eint_loc_R)

      ! MHD-style tau formula, using face quantities + global cmax.
      ! sigma_T7_f = kappa_f * T_f acts as the conductivity*temperature
      ! product entering tau. Floor tau at 4*dt for stability.
      sigma_T7_f = kappa_f * 0.5d0 * (T_L + T_R)
      tau_f = max(4.0d0 * dt, &
              sigma_T7_f * courantpar**2 / (eint_loc_f * cmax_global**2))

      ! Cattaneo relaxation: closed-form exponential
      ratio = min(qdt / tau_f, ratio_cap)
      decay = dexp(-ratio)
      if (ratio > 1.0d-6) then
        ave_factor = (1.0d0 - decay) / ratio
      else
        ave_factor = 1.0d0 - 0.5d0 * ratio + ratio * ratio / 6.0d0
      end if

      qn_face = 0.5d0 * (wCT(ix1, q_) + wCT(ix1 + 1, q_))

      qf_full(ix1) = q_Sp_star + (qn_face - q_Sp_star) * decay
      qf_half(ix1) = q_Sp_star + (qn_face - q_Sp_star) * ave_factor

      ! Energy-positivity clip on the n+1/2 flux
      q_face_max = hd_htc_pos_eta * min(eint_loc_L, eint_loc_R) * &
                   block%ds(ix1, 1) / max(qdt, smalldouble)
      qf_half(ix1) = sign(min(dabs(qf_half(ix1)), q_face_max), qf_half(ix1))
    end do

    ! Conservative energy update + cell-centred q refresh
    do ix1 = ixOmin1, ixOmax1
      w(ix1, e_) = w(ix1, e_) &
                 - qdt * (qf_half(ix1) - qf_half(ix1 - 1)) / block%ds(ix1, 1)
      w(ix1, q_) = 0.5d0 * (qf_full(ix1 - 1) + qf_full(ix1))
    end do
    }

  end subroutine add_hypertc_source

  subroutine hd_add_radiation_source(qdt,ixI^L,ixO^L,wCT,wCTprim,w,x,qsourcesplit,active)
    use mod_constants
    use mod_global_parameters
    use mod_usr_methods
    use mod_fld

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, x(ixI^S,1:ndim)
    double precision, intent(in)    :: wCT(ixI^S,1:nw),wCTprim(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    logical, intent(in) :: qsourcesplit
    logical, intent(inout) :: active

    ! add radiation force and work done by it, changes momentum and gas energy
    ! handle photon tiring, heating and cooling exchange between gas and radiation field
    call add_fld_rad_force(qdt,ixI^L,ixO^L,wCT,wCTprim,w,x,qsourcesplit,active,fld_fl)

  end subroutine hd_add_radiation_source

  subroutine hd_get_dt(wprim, ixI^L, ixO^L, dtnew, dx^D, x)
    use mod_global_parameters
    use mod_dust, only: dust_get_dt
    use mod_viscosity, only: viscosity_get_dt
    use mod_gravity, only: gravity_get_dt
    use mod_cak_force, only: cak_get_dt
    use mod_fld, only: fld_radforce_get_dt

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: dx^D, x(ixI^S, 1:^ND)
    double precision, intent(in)    :: wprim(ixI^S, 1:nw)
    double precision, intent(inout) :: dtnew

    dtnew = bigdouble

    if(hd_dust) then
      call dust_get_dt(wprim, ixI^L, ixO^L, dtnew, dx^D, x)
    end if

    if(hd_viscosity) then
      call viscosity_get_dt(wprim,ixI^L,ixO^L,dtnew,dx^D,x)
    end if

    if(hd_gravity) then
      call gravity_get_dt(wprim,ixI^L,ixO^L,dtnew,dx^D,x)
   end if

   if (hd_cak_force) then
     call cak_get_dt(wprim,ixI^L,ixO^L,dtnew,dx^D,x)
   end if

   if(hd_radiation_fld) then
     call fld_radforce_get_dt(wprim,ixI^L,ixO^L,dtnew,dx^D,x,fld_fl)
   endif

  end subroutine hd_get_dt

  !> Wrappers for the FLD implicit (MG diffusion) hooks: phys_implicit_update /
  !> phys_evaluate_implicit have fixed interfaces with no fluid argument, so
  !> these inject the module's fld_fl object into the threaded fld routines.
  subroutine hd_fld_implicit_update(dtfactor,qdt,qtC,psa,psb)
    use mod_global_parameters
    use mod_fld, only: fld_implicit_update
    type(state), target          :: psa(max_blocks)
    type(state), target          :: psb(max_blocks)
    double precision, intent(in) :: qdt
    double precision, intent(in) :: qtC
    double precision, intent(in) :: dtfactor

    call fld_implicit_update(dtfactor,qdt,qtC,psa,psb,fld_fl)
  end subroutine hd_fld_implicit_update

  subroutine hd_fld_evaluate_implicit(qtC,psa)
    use mod_global_parameters
    use mod_fld, only: fld_evaluate_implicit
    type(state), target          :: psa(max_blocks)
    double precision, intent(in) :: qtC

    call fld_evaluate_implicit(qtC,psa,fld_fl)
  end subroutine hd_fld_evaluate_implicit

  function hd_kin_en(w, ixI^L, ixO^L, inv_rho) result(ke)
    use mod_global_parameters, only: nw, ndim
    integer, intent(in)                    :: ixI^L, ixO^L
    double precision, intent(in)           :: w(ixI^S, nw)
    double precision                       :: ke(ixO^S)
    double precision, intent(in), optional :: inv_rho(ixO^S)

    if (present(inv_rho)) then
       ke = 0.5d0 * sum(w(ixO^S, mom(:))**2, dim=ndim+1) * inv_rho
    else
       ke = 0.5d0 * sum(w(ixO^S, mom(:))**2, dim=ndim+1) / w(ixO^S, rho_)
    end if
  end function hd_kin_en

  function hd_inv_rho(w, ixI^L, ixO^L) result(inv_rho)
    use mod_global_parameters, only: nw, ndim
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S, nw)
    double precision              :: inv_rho(ixO^S)

    ! Can make this more robust
    inv_rho = 1.0d0 / w(ixO^S, rho_)
  end function hd_inv_rho

  subroutine hd_handle_small_values(primitive, w, x, ixI^L, ixO^L, subname)
    ! handles hydro (density,pressure,velocity) bootstrapping
    ! any negative dust density is flagged as well (and throws an error)
    ! small_values_method=replace also for dust
    use mod_global_parameters
    use mod_small_values
    use mod_dust, only: dust_n_species, dust_mom, dust_rho
    logical, intent(in)             :: primitive
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    character(len=*), intent(in)    :: subname

    integer :: n,idir
    logical :: flag(ixI^S,1:nw)

    call hd_check_w(primitive, ixI^L, ixO^L, w, flag)

    if (any(flag)) then
      select case (small_values_method)
      case ("replace")
        where(flag(ixO^S,rho_)) w(ixO^S,rho_) = small_density
        do idir = 1, ndir
          if(small_values_fix_iw(mom(idir))) then
            where(flag(ixO^S,rho_)) w(ixO^S, mom(idir)) = 0.0d0
          end if
        end do
        if(hd_radiation_fld)then
          if (small_values_fix_iw(r_e)) then
            where(flag(ixO^S,r_e)) w(ixO^S,r_e) = small_r_e
          end if
        end if
        if(hd_energy)then
          if(small_values_fix_iw(e_)) then
            if(primitive) then
              where(flag(ixO^S,rho_)) w(ixO^S, p_) = small_pressure
            else
              where(flag(ixO^S,rho_)) w(ixO^S, e_) = small_e + hd_kin_en(w,ixI^L,ixO^L)
            endif
          end if
        endif

        if(hd_energy) then
          if(primitive) then
            where(flag(ixO^S,e_)) w(ixO^S,p_) = small_pressure
          else
            where(flag(ixO^S,e_))
              ! Add kinetic energy
              w(ixO^S,e_) = small_e + hd_kin_en(w,ixI^L,ixO^L)
            end where
          end if
        end if

        if(hd_dust)then
           do n=1,dust_n_species
              where(flag(ixO^S,dust_rho(n))) w(ixO^S,dust_rho(n)) = 0.0d0
              do idir = 1, ndir
                  where(flag(ixO^S,dust_rho(n))) w(ixO^S,dust_mom(idir,n)) = 0.0d0
              enddo
           enddo
        endif
      case ("average")
        if(primitive)then
           ! averaging for all primitive fields, including dust
           call small_values_average(ixI^L, ixO^L, w, x, flag)
        else
           ! do averaging of density
           call small_values_average(ixI^L, ixO^L, w, x, flag, rho_)
           if(hd_energy) then
             ! do averaging of pressure
            !  w(ixI^S,p_)=(hd_gamma-1.d0)*(w(ixI^S,e_) &
            !   -0.5d0*sum(w(ixI^S, mom(:))**2, dim=ndim+1)/w(ixI^S,rho_))
             call eos%get_thermal_pressure(w, x, ixI^L, ixO^L, w(ixO^S,p_))
             call small_values_average(ixI^L, ixO^L, w, x, flag, p_)
             do idir = 1, ndir
               w(ixO^S,mom(idir)) = w(ixO^S,mom(idir))/w(ixO^S,rho_) !> Convert to velocity to be compliant with p_to_e
             end do
             call eos%p_to_e(ixI^L, ixO^L, w, x)
             do idir = 1, ndir
               w(ixO^S,mom(idir)) = w(ixO^S,mom(idir))*w(ixO^S,rho_) !> Restore momentum (conserved-form invariant)
             end do
              ! w(ixI^S,e_)=w(ixI^S,p_)/(hd_gamma-1.d0) &
              !  +0.5d0*sum(w(ixI^S, mom(:))**2, dim=ndim+1)/w(ixI^S,rho_)
           end if
           if(hd_radiation_fld) then
              ! do averaging of radiative energy density
              call small_values_average(ixI^L, ixO^L, w, x, flag, r_e)
           endif
           if(hd_dust)then
              do n=1,dust_n_species
                 where(flag(ixO^S,dust_rho(n))) w(ixO^S,dust_rho(n)) = 0.0d0
                 do idir = 1, ndir
                    where(flag(ixO^S,dust_rho(n))) w(ixO^S,dust_mom(idir,n)) = 0.0d0
                 enddo
              enddo
          endif
        endif
      case default
        if(.not.primitive) then
          !convert w to primitive
          ! Calculate pressure = (gamma-1) * (e-ek)
          if(hd_energy) then
            call eos%get_thermal_pressure(w, x, ixI^L, ixO^L, w(ixO^S,p_))
            ! w(ixO^S,p_)=(hd_gamma-1.d0)*(w(ixO^S,e_)-hd_kin_en(w,ixI^L,ixO^L))
          end if
          ! Convert gas momentum to velocity
          do idir = 1, ndir
            w(ixO^S, mom(idir)) = w(ixO^S, mom(idir))/w(ixO^S,rho_)
          end do
        end if
        ! NOTE: dust entries may still have conserved values here
        call small_values_error(w, x, ixI^L, ixO^L, flag, subname)
      end select
    end if
    if (hd_fip) call hd_bound_fip(primitive, ixI^L, ixO^L, w)
  end subroutine hd_handle_small_values

  !> Well-balanced transform: (T, v, q) variable change.
  !>
  !> Replaces (ρ, v, p) with (T, v, q) where:
  !>   T(i) = p(i) / ρ(i)         [temperature — smooth, locally computed]
  !>   q(i) = p(i) / p_eq(i)      [pressure ratio — ≈1 in HSE]
  !> Well-balanced post-prolongation correction for AMR.
  !>
  !> After prolongation interpolates to a fine grid, the pressure does not
  !> satisfy the discrete HSE recurrence at the fine resolution (linear
  !> interpolation of an exponential profile introduces O(dx^2/H^2) error).
  !>
  !> This routine rebuilds p from the multiplicative HSE recurrence using the
  !> interpolated T = p/rho (which IS smooth and well-interpolated). The
  !> density is updated as rho = p_new / T for consistency.
  !>
  !> On entry:  w contains primitive (rho, v, p).
  !> On exit:   w contains HSE-corrected primitive (rho_new, v, p_new).
  !>
  !> The recurrence is anchored at the block midpoint, where the parent
  !> cell-centre value is exact (no interpolation error).
  subroutine hd_wb_prolong(ixI^L, ixO^L, w, x)
    use mod_global_parameters
    use mod_usr_methods, only: usr_gravity

    integer, intent(in)              :: ixI^L, ixO^L
    double precision, intent(inout)  :: w(ixI^S, 1:nw)
    double precision, intent(in)     :: x(ixI^S, 1:ndim)

    double precision :: gravity_field(ixI^S, 1:ndim)
    double precision :: wb_T(ixI^S), p_eq(ixI^S)
    double precision :: dx_idims
    double precision :: alpha(ixI^S), beta(ixI^S)
    {^IFONED
    integer :: ix1, ix_mid
    }

    ! T = p/rho at all fine cells (smooth from interpolation)
    wb_T(ixO^S) = w(ixO^S, p_) / w(ixO^S, rho_)

    ! Get gravity at fine cell centres
    call usr_gravity(ixI^L, ixO^L, w, x, gravity_field)

    dx_idims = dxlevel(1)

    {^IFONED
    ! α_i = (dx/2) · g_i / T_i  (vectorised over ixO)
    alpha(ixOmin1:ixOmax1) = 0.5d0 * dx_idims &
      * gravity_field(ixOmin1:ixOmax1, 1) / wb_T(ixOmin1:ixOmax1)

    ! Anchor at block midpoint (least interpolation error)
    ix_mid = (ixOmin1 + ixOmax1) / 2
    p_eq(ix_mid) = w(ix_mid, p_)

    ! Forward: β_i = (1 + α_i) / (1 - α_{i+1})  (vectorised)
    beta(ix_mid:ixOmax1-1) = (1.0d0 + alpha(ix_mid:ixOmax1-1)) &
      / (1.0d0 - alpha(ix_mid+1:ixOmax1))
    do ix1 = ix_mid + 1, ixOmax1
      p_eq(ix1) = p_eq(ix1 - 1) * beta(ix1 - 1)
    end do

    ! Backward: β_i = (1 - α_{i+1}) / (1 + α_i)  (vectorised)
    beta(ixOmin1:ix_mid-1) = (1.0d0 - alpha(ixOmin1+1:ix_mid)) &
      / (1.0d0 + alpha(ixOmin1:ix_mid-1))
    do ix1 = ix_mid - 1, ixOmin1, -1
      p_eq(ix1) = p_eq(ix1 + 1) * beta(ix1)
    end do

    ! Replace interpolated p with recurrence-derived p, update rho = p/T
    w(ixOmin1:ixOmax1, p_)   = p_eq(ixOmin1:ixOmax1)
    w(ixOmin1:ixOmax1, rho_) = p_eq(ixOmin1:ixOmax1) / wb_T(ixOmin1:ixOmax1)
    }

  end subroutine hd_wb_prolong

  !> Well-balanced transform for reconstruction.
  !>
  !> On entry:  w contains primitive (rho, v, p).
  !> On exit:   w(rho_) = T = p/rho, w(p_) = q = p/p_eq.
  !>            wb_T returns the saved cell-centre temperature.
  !>            wb_phi returns the cell-centre equilibrium pressure.
  !>            wb_phi_face returns the face-centre equilibrium pressure.
  !>
  !> p_eq from multiplicative trapezoidal recurrence (always positive).
  !> In HSE: q = 1, T varies smoothly -> limiter sees flat q -> exact balance.
  subroutine hd_wb_transform(ixI^L, ixO^L, idims, w, x, wb_phi, &
       wb_phi_face, wb_T)
    use mod_global_parameters
    use mod_usr_methods, only: usr_gravity

    integer, intent(in)              :: ixI^L, ixO^L, idims
    double precision, intent(inout)  :: w(ixI^S, 1:nw)
    double precision, intent(in)     :: x(ixI^S, 1:ndim)
    double precision, intent(out)    :: wb_phi(ixI^S)
    double precision, intent(out)    :: wb_phi_face(ixI^S)
    double precision, intent(out)    :: wb_T(ixI^S)

    double precision :: gravity_field(ixI^S, 1:ndim)
    double precision :: dx_idims
    double precision :: alpha(ixI^S), beta(ixI^S)
    {^IFONED
    integer :: ix1
    }

    ! Save T = p/ρ at all cell centers
    wb_T(ixI^S) = w(ixI^S, p_) / w(ixI^S, rho_)

    ! Get gravity acceleration at all cell centers
    call usr_gravity(ixI^L, ixI^L, w, x, gravity_field)

    dx_idims = dxlevel(idims)

    ! Multiplicative trapezoidal recurrence for p_eq.
    ! p_eq(i+1) = p_eq(i) · (1 + α(i)) / (1 - α(i+1))
    ! Always positive for dx < 2H.
    !
    ! Vectorised: precompute α and β factors, then sequential cumulative product.
    {^IFONED
    ! α_i = (dx/2) · g_i / T_i  (vectorised)
    alpha(ixImin1:ixImax1) = 0.5d0 * dx_idims &
      * gravity_field(ixImin1:ixImax1, idims) / wb_T(ixImin1:ixImax1)

    ! β_i = (1 + α_i) / (1 - α_{i+1})  (vectorised)
    beta(ixImin1:ixImax1-1) = (1.0d0 + alpha(ixImin1:ixImax1-1)) &
      / (1.0d0 - alpha(ixImin1+1:ixImax1))

    ! Cumulative product (sequential, unavoidable data dependency)
    wb_phi(ixImin1) = w(ixImin1, p_)
    do ix1 = ixImin1 + 1, ixImax1
      wb_phi(ix1) = wb_phi(ix1 - 1) * beta(ix1 - 1)
    end do
    }

    ! Face equilibrium pressure: isothermal half-step from cell center
    wb_phi_face(ixI^S) = wb_phi(ixI^S) * (1.0d0 + 0.5d0 * dx_idims * &
      gravity_field(ixI^S, idims) / wb_T(ixI^S))

    ! Transform: w(rho_) = T, w(p_) = q = p/p_eq
    w(ixI^S, rho_) = wb_T(ixI^S)
    w(ixI^S, p_) = w(ixI^S, p_) / wb_phi(ixI^S)

  end subroutine hd_wb_transform

  !> Well-balanced inverse: restore physical (rho, v, p) at interfaces.
  !>
  !> On entry:
  !>   wLp/wRp(rho_) = T (reconstructed), wLp/wRp(p_) = q (reconstructed)
  !>   w(rho_) = T (cell-centre), w(p_) = q (cell-centre)
  !>   wb_phi/wb_phi_face = cell/face equilibrium pressures from transform
  !>   wb_T = saved cell-centre temperature from transform
  !>
  !> On exit:
  !>   wLp/wRp(rho_) = rho_face, wLp/wRp(p_) = p_face (physical primitives)
  !>   w(rho_) = rho_cell, w(p_) = p_cell (restored cell-centre primitives)
  !>
  !> Pressure: p_face = q_face * p_eq_face   (multiplicative, well-balanced)
  !> Density:  rho_face = p_face / T_blended
  !>
  !> Blends between shared T (well-balanced, zero HLL dissipation) and
  !> individual limiter-reconstructed T (non-WB, full dissipation) using
  !> a contact detector sigma = |T_faceL - T_faceR| / T_avg.
  !> In HSE: T is smooth → sigma ≈ 0 → pure WB (ρ_L = ρ_R).
  !> At contacts: T jump → sigma > 0 → ρ_L ≠ ρ_R → dissipation restored.
  subroutine hd_wb_inverse(ixI^L, ixL^L, ixR^L, idims, wLp, wRp, w, &
       wb_phi, wb_phi_face, wb_T)
    use mod_global_parameters

    integer, intent(in)              :: ixI^L, ixL^L, ixR^L, idims
    double precision, intent(inout)  :: wLp(ixI^S, 1:nw), wRp(ixI^S, 1:nw)
    double precision, intent(inout)  :: w(ixI^S, 1:nw)
    double precision, intent(in)     :: wb_phi(ixI^S), wb_phi_face(ixI^S)
    double precision, intent(in)     :: wb_T(ixI^S)

    double precision :: T_shared(ixI^S)
    double precision :: T_face_L(ixI^S), T_face_R(ixI^S)
    double precision :: sigma(ixI^S), T_for_rhoL(ixI^S), T_for_rhoR(ixI^S)

    ! Shared face temperature: T_face(i) = ½(T(i) + T(i+1))
    {^IFONED
    T_shared(ixImin1:ixImax1-1) = 0.5d0 * (wb_T(ixImin1:ixImax1-1) &
      + wb_T(ixImin1+1:ixImax1))
    T_shared(ixImax1) = wb_T(ixImax1)
    }

    ! Save limiter-reconstructed T before rho_ slot is overwritten
    ! (rho_ still holds T from the WB transform at this point)
    T_face_L(ixL^S) = wLp(ixL^S, rho_)
    T_face_R(ixR^S) = wRp(ixR^S, rho_)

    ! Contact detector: relative T jump across interface
    ! sigma = 0 in HSE (smooth T), sigma > 0 at contacts (T discontinuity)
    sigma(ixL^S) = dabs(T_face_L(ixL^S) - T_face_R(ixR^S)) &
                 / (0.5d0 * (T_face_L(ixL^S) + T_face_R(ixR^S)) + smalldouble)
    sigma(ixL^S) = min(sigma(ixL^S), 1.0d0)

    ! Blended T for density recovery
    T_for_rhoL(ixL^S) = (1.d0 - sigma(ixL^S)) * T_shared(ixL^S) &
                       + sigma(ixL^S) * T_face_L(ixL^S)
    T_for_rhoR(ixR^S) = (1.d0 - sigma(ixL^S)) * T_shared(ixR^S) &
                       + sigma(ixL^S) * T_face_R(ixR^S)

    ! Pressure inverse: p = q · p_eq_face (multiplicative)
    wLp(ixL^S, p_) = wLp(ixL^S, p_) * wb_phi_face(ixL^S)
    wRp(ixR^S, p_) = wRp(ixR^S, p_) * wb_phi_face(ixR^S)

    ! Density inverse: ρ = p / T_blended
    wLp(ixL^S, rho_) = wLp(ixL^S, p_) / T_for_rhoL(ixL^S)
    wRp(ixR^S, rho_) = wRp(ixR^S, p_) / T_for_rhoR(ixR^S)

    ! Restore cell-center values
    w(ixI^S, p_) = w(ixI^S, p_) * wb_phi(ixI^S)
    w(ixI^S, rho_) = w(ixI^S, p_) / wb_T(ixI^S)

  end subroutine hd_wb_inverse

end module mod_hd_phys
