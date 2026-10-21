!> module radiative cooling -- add optically thin radiative cooling 
!> 
!> only uses the (Townsend) exact integration method, can be used in HD, ffhd, MHD, twofl
!>
!> Assumptions: full ionized plasma dominated by H and He, ionization equilibrium 
!> Formula: Q=-n_H*n_e*f(T), positive f(T) function is pre-computed and tabulated or a piecewise power law
!> Uses the various cooling tables stored in mod_radloss_tables.t
!>
module mod_radiative_cooling

  use mod_global_parameters, only: std_len
  use mod_physics
  use mod_radloss_tables
  use mod_comm_lib, only: mpistop
  implicit none

  !> Per-rank cooling-only compute accumulator for lb_diagnose. Sums the
  !> wall time spent inside radiative_cooling_add_source across all blocks
  !> and substages within one full advance() call. Reset by
  !> mod_advance::advance at the start of each step. Note: this is a subset
  !> of lb_compute_accum (which already encloses the whole advect1 block
  !> loop); it isolates the cooling Newton/exact-solver kernel cost from
  !> the surrounding finite-volume work.
  double precision, public :: lb_cool_accum = 0.0d0

  !> Helium abundance over Hydrogen
  double precision, private    :: He_abundance

  !> The adiabatic index
  double precision, private :: rc_gamma

  !> The adiabatic index minus 1
  double precision, private :: rc_gamma_1

  !> inverse of the adiabatic index minus 1
  double precision, private :: invgam

  !> Voigt escape probability lookup table.
  !> E(tau) is the frequency-integrated single-flight escape probability
  !> from a uniform slab with a Voigt line profile truncated at x_max
  !> Doppler widths (CL12 Sec 2.1, mimics PRD):
  !>   E(tau_0) = [phi(0)/tau_0] * integral_{-xmax}^{xmax} [1 - exp(-tau_0 phi(x)/phi(0))] dx
  !> Precomputed at initialisation via Gauss-Legendre quadrature, then
  !> interpolated at runtime from a table in log10(tau).
  !> Reference Voigt parameter: a_ref = 4.7e-4 at T = 10 kK.
  !> Sensitivity to a is < 20% for tau < 1e6 (see voigt_escape_derivation.py).
  integer, parameter, private :: n_voigt_table = 500
  double precision, parameter, private :: voigt_logtau_min = -2.0d0   ! tau = 0.01
  double precision, parameter, private :: voigt_logtau_max =  8.0d0   ! tau = 1e8
  double precision, parameter, private :: voigt_a_ref = 4.7d-4
  double precision, parameter, private :: voigt_xmax  = 6.0d0
  double precision, private :: voigt_E_table(n_voigt_table)
  double precision, private :: voigt_logtau_step
  logical, private :: voigt_table_ready = .false.

  abstract interface
    subroutine get_subr1(w,x,ixI^L,ixO^L,res)
      use mod_global_parameters
      integer, intent(in)          :: ixI^L, ixO^L
      double precision, intent(in) :: w(ixI^S,nw)
      double precision, intent(in) :: x(ixI^S,1:ndim)
      double precision, intent(out):: res(ixI^S)
    end subroutine get_subr1

    subroutine get_2var_subr(ixI^L, ixO^L, w, x, ne, nH)
      use mod_global_parameters
      integer, intent(in)          :: ixI^L, ixO^L
      double precision, intent(in) :: w(ixI^S, nw)
      double precision, intent(in) :: x(ixI^S, 1:ndim)
      double precision, intent(out):: ne(ixI^S), nH(ixI^S)
    end subroutine get_2var_subr

    !> Optional local multiplier for density-squared radiative losses.
    subroutine local_rho2_factor_subr(ixI^L,ixO^L,w,x,factor)
      use mod_global_parameters
      integer, intent(in) :: ixI^L, ixO^L
      double precision, intent(in) :: w(ixI^S,1:nw), x(ixI^S,1:ndim)
      double precision, intent(out) :: factor(ixI^S)
    end subroutine local_rho2_factor_subr

    !> Scalar EoS inverse, e.g. fl%eint_from_T(log_nH, log_T)
    double precision function eos_scalar2_func(a, b)
      double precision, intent(in) :: a, b
    end function eos_scalar2_func
  end interface

  type rc_fluid

    double precision :: rad_damp_height
    double precision :: rad_damp_scale

    ! these are set in init method
    double precision, allocatable :: tcool(:), Lcool(:), dLdtcool(:)
    double precision, allocatable :: Yc(:)
    double precision  :: tref, lref, tcoolmin,tcoolmax
    double precision  :: lgtcoolmin, lgtcoolmax, lgstep

    ! The piecewise powerlaw (PPL) tabels and variabels
    ! x_* en t_* are given as log_10
    double precision, allocatable :: y_PPL(:), t_PPL(:), l_PPL(:), a_PPL(:)

    !> Lower limit of temperature
    double precision   :: tlow

    !> Index of the energy density
    integer              :: e_
    !> Index of cut off temperature for TRAC
    integer              :: Tcoff_

    ! these are set as parameters
    !> Resolution of temperature in interpolated tables
    integer :: ncool

    integer :: n_PPL

    !> Fixed temperature not lower than tlow
    logical   :: Tfix

    !> Add cooling source in a split way (.true.) or un-split way (.false.)
    logical    :: rc_split

    logical :: isPPL = .false.

    !> Suppress cooling for T below this threshold (Kelvin) within
    !> rad_cut_hgt of footpoints. When > 0, cells inside the spatial
    !> taper region with T < rad_suppress_temp get factor=0; coronal
    !> cooling above that T proceeds unmodified. Default 0 = disabled.
    double precision :: rad_suppress_temp = 0.0d0
    !> Internal: suppress threshold in code units (set from rad_suppress_temp
    !> during init). Not a namelist parameter.
    double precision :: suppress_temp_code = 0.0d0
    !> Master switch for radiative loss modification (spatial + density taper)
    logical :: rad_modify
    !> Apply spatial taper at both boundaries (default: lower only)
    logical :: rad_modify_sym
    !> Spatial taper height (HEAD addition): suppress cooling within rad_cut_hgt of boundary
    double precision :: rad_cut_hgt = 0.0d0
    !> Spatial taper width (HEAD addition): Gaussian dey for the spatial taper
    double precision :: rad_cut_dey = 0.15d0
    !> Cooling-curve dt-fraction (HEAD addition; used by legacy explicit-mode dt scaling)
    double precision :: cfrac = 0.1d0
    !> cutoff radiative cooling below rad_damp_height
    logical :: rad_damp
    ! these are to be set directly
    !> whether background equilibrium contribution is split off
    logical :: has_equi = .false.
    !> whether background equilibrium is compensated in thermal balance
    logical :: subtract_equi = .false.

    double precision, allocatable :: frac_lowFIP(:)
    !> Index of primitive FIP abundance variable, -1 if disabled
    integer :: fip_ = -1
    !> Enable local Newton cooling/heating approximation for optically thick losses
    logical :: rad_newton = .false.
    double precision :: rad_newton_pthick = 25.d0
    double precision :: rad_newton_trad = 0.006d0
    double precision :: rad_newton_rhosurf = 1.d4

    !> Density threshold for Gaussian taper (code units)
    double precision :: rad_taper_rho
    !> Gaussian decay width for density taper
    double precision :: rad_taper_dey

    !> Enable escape probability cooling modification
    logical :: rad_escape_prob = .false.
    !> Effective opacity for escape probability (code units: 1/(density*length))
    double precision :: rad_kappa_eff = 0.0d0
    !> Temperature above which kappa→0 (stored in code units, input in Kelvin); 0 = constant kappa
    double precision :: rad_kappa_Tcutoff = 0.0d0
    !> Sigmoid sharpness exponent for kappa(T) cutoff
    double precision :: rad_kappa_alpha = 4.0d0
    !> Escape probability type: 'slab' = (1-exp(-tau))/tau, 'voigt' = Voigt CRD
    character(len=10) :: rad_escape_type = 'slab'
    !> Exponential cutoff scale: E *= exp(-tau/tau_cutoff); 0 = disabled
    double precision :: rad_escape_tau_cutoff = 0.0d0
    !> Max height from footpoint for column mass integration (code units); 0 = no limit
    double precision :: rad_escape_height = 0.0d0
    !> Index into wextra for column mass (set during init)
    integer :: iw_colmass_ = -1

    !> Name of cooling curve
    character(len=std_len)  :: coolcurve

    procedure (get_subr1), pointer, nopass :: get_rho => null()
    procedure (get_subr1), pointer, nopass :: get_Te => null()
    procedure (get_subr1), pointer, nopass :: get_rho_equi => null()
    procedure (get_subr1), pointer, nopass :: get_pthermal => null()
    procedure (get_subr1), pointer, nopass :: get_pthermal_equi => null()
    procedure (get_subr1), pointer, nopass :: get_var_Rfactor => null()
    procedure (get_2var_subr), pointer, nopass :: get_ne_nH => null()
    procedure (get_2var_subr), pointer, nopass :: get_ne_nH_equi => null()
    procedure (local_rho2_factor_subr), pointer, nopass :: get_rho2_factor => null()
    procedure (get_subr1), pointer, nopass :: get_temperature_equi => null()
    !> EoS snapshots + scalar inverse accessors (set in bind_eos_to_source); let
    !> cooling reach thermodynamics only through this object, never mod_eos.
    logical            :: ionE = .false.
    character(len=20)  :: method = 'tables'
    double precision   :: inv_gamma_minus_1
    double precision   :: nH2rhoFactor
    double precision   :: eion_per_nH
    procedure(eos_scalar2_func), pointer, nopass :: eint_from_T => null()
    procedure(eos_scalar2_func), pointer, nopass :: p2eint      => null()
    procedure(eos_scalar2_func), pointer, nopass :: T_from_eint => null()
    procedure(eos_scalar2_func), pointer, nopass :: y_from_eint => null()

    !> Variable-c_V Townsend extension (Y_mod). Built only when fl%ionE.
    !>
    !>   Y_mod(j, i) is the modified TEF (units of code time) at the
    !>   (log10 nH index j, T index i) grid point. Indexing follows the
    !>   existing AMRVAC convention: Y(ncool) = 0 at the top T = tcoolmax,
    !>   and Y monotonically increases as T decreases. The j axis matches
    !>   the eos%eint_from_T table's nH grid (var1_min..var1_max, dim1).
    !>
    !> Construction: change-of-variables u = e_int/n_H. The integrand
    !>   1/(n_e Lambda) is sampled at composite Simpson or Boole nodes
    !>   between [u(T_i), u(T_{i+1})]. See build_Y_mod_table.
    !>
    !> Lookups: findY_mod(T, nH, fl) returns Y by bilinear interpolation
    !>   in (log_nH, log_T). findT_mod(Y, nH, fl) returns T by bisection
    !>   on the row at the interpolated nH. (A precomputed inverse table
    !>   variant was prototyped during development but proved unusable:
    !>   at extreme nH the per-row Y_max varies so widely that bilinear
    !>   blending between rows mixes physically saturated and unsaturated
    !>   entries. The bisect path is O(log ncool) ≈ 12 iterations anyway.)
    double precision, allocatable :: Y_mod(:,:)
    double precision, allocatable :: Y_mod_max_per_row(:)
    integer :: Y_mod_n_nH = 0
    double precision :: Y_mod_lg_nH_min = 0.0d0
    double precision :: Y_mod_lg_nH_max = 0.0d0
    double precision :: Y_mod_lg_nH_step_inv = 0.0d0
    !> Build flag — set to .true. only after the table has been populated
    !> by build_Y_mod_table (called from bind_eos_to_source after eos_finalise).
    logical :: Y_mod_built = .false.
    !> Quadrature method: 'simpson' (3-point, O(h^4)) or 'boole' (5-point, O(h^6))
    character(len=8) :: Y_mod_quadrature = 'boole'
    !> Number of sub-intervals per [u_i, u_{i+1}] segment for the quadrature
    integer :: Y_mod_N_sub = 16

    !> SPEX-style two-table cooling support.
    !>
    !> The SPEX/SPEX_DM cooling tables follow Schure et al. (2009) which
    !> publishes the cooling function in two parts:
    !>   Lambda_SPEX(T)  -- the cooling rate per n_H^2 (NOT per n_e n_H)
    !>   nenh_SPEX(T)    -- the CIE equilibrium n_e/n_H ratio at temperature T
    !> Reconstructing the volumetric cooling rate is then:
    !>   Q = n_H^2 * nenh_eq(T) * Lambda_SPEX(T)
    !>
    !> All other cooling tables in AMRVAC follow the standard Dere/CHIANTI
    !> convention where Q = n_e * n_H * Lambda(T) and the equilibrium
    !> ionisation balance is baked into Lambda(T) itself.
    !>
    !> Historically, AMRVAC handled the SPEX two-table convention by
    !> absorbing log10(nenh_SPEX) into Lambda_table at construction time, so
    !> that the standard formula Q = n_e n_H * Lambda_table happened to give
    !> the right answer when n_e ~ n_H (the FI assumption with neOnH ~ 1.2).
    !> This trick silently *breaks* in LTE+ionE mode where the simulation
    !> n_e is the actual Saha value: n_e/n_H << 1 at low T, so the formula
    !> double-counts the equilibrium factor and badly under-counts cooling.
    !>
    !> The fix below: do NOT absorb nenh into Lambda_table. Instead, store
    !> the equilibrium array nenh_eq_table on the same tcool grid, and at
    !> runtime use Q = n_H^2 * nenh_eq(T) * Lambda_table whenever
    !> lambda_needs_nenh_table is .true. This honours the published SPEX
    !> convention regardless of the EoS choice.
    logical :: lambda_needs_nenh_table = .false.
    double precision, allocatable :: nenh_eq_table(:)

  end type rc_fluid

  contains

    subroutine radiative_cooling_rho2_factor(ixI^L,ixO^L,w,x,fl,factor)
      use mod_global_parameters
      integer, intent(in) :: ixI^L, ixO^L
      double precision, intent(in) :: w(ixI^S,1:nw), x(ixI^S,1:ndim)
      type(rc_fluid), intent(in) :: fl
      double precision, intent(out) :: factor(ixI^S)
      factor(ixO^S)=one
      if(associated(fl%get_rho2_factor)) &
        call fl%get_rho2_factor(ixI^L,ixO^L,w,x,factor)
    end subroutine radiative_cooling_rho2_factor

    !> Radiative cooling initialization
    subroutine radiative_cooling_init_params(phys_gamma,He_abund)
      use mod_global_parameters
      double precision, intent(in) :: phys_gamma,He_abund

      rc_gamma=phys_gamma
      He_abundance=He_abund
    end subroutine radiative_cooling_init_params

    !> Build the Voigt escape probability lookup table.
    !> Called once (guarded by voigt_table_ready flag).
    !> Uses 64-point Gauss-Legendre quadrature on [0, x_max] to evaluate
    !>   E(tau_0) = (2 phi(0)/tau_0) * integral_0^{x_max} [1-exp(-tau_0 g(x))] dx
    !> where g(x) = phi(x)/phi(0) for the Voigt profile H(a,x) truncated at x_max.
    subroutine voigt_escape_init_table()
      use mod_global_parameters, only: dpi
      implicit none
      integer, parameter :: nquad = 64
      double precision :: xq(nquad), wq(nquad)
      double precision :: logtau, tau0, phi0, gx, integrand, E_val
      double precision :: a_rep
      integer :: i, k

      if(voigt_table_ready) return

      ! Step size in log10(tau)
      voigt_logtau_step = (voigt_logtau_max - voigt_logtau_min) / dble(n_voigt_table - 1)

      ! Get Gauss-Legendre nodes and weights on [0, x_max]
      call voigt_gauss_legendre(0.0d0, voigt_xmax, nquad, xq, wq)

      ! Representative Voigt parameter (a_ref at T = 10 kK)
      a_rep = voigt_a_ref

      ! phi(0) for the Voigt profile: H(a,0) = exp(a^2)*erfc(a)/sqrt(pi) ≈ 1/sqrt(pi)
      phi0 = 1.0d0 / sqrt(dpi)

      ! Precompute g(x) = phi(x)/phi(0) at quadrature nodes
      ! For the Voigt profile: phi(x) = H(a,x), the real part of the
      ! Faddeeva function divided by sqrt(pi).
      ! For small a, H(a,x) ≈ exp(-x^2)/sqrt(pi) + a/(pi*x^2) for |x|>few.
      ! Use the exact Humlicek (1982) rational approximation.

      do i = 1, n_voigt_table
        logtau = voigt_logtau_min + dble(i-1) * voigt_logtau_step
        tau0 = 10.0d0**logtau

        if(tau0 < 1.0d-6) then
          voigt_E_table(i) = 1.0d0
          cycle
        end if

        integrand = 0.0d0
        do k = 1, nquad
          gx = voigt_profile_ratio(a_rep, xq(k))  ! phi(x)/phi(0)
          ! [1 - exp(-tau0 * g(x))]
          if(tau0 * gx > 500.0d0) then
            integrand = integrand + wq(k) * 1.0d0
          else if(tau0 * gx < 1.0d-10) then
            integrand = integrand + wq(k) * tau0 * gx
          else
            integrand = integrand + wq(k) * (1.0d0 - exp(-tau0 * gx))
          end if
        end do

        E_val = 2.0d0 * phi0 / tau0 * integrand
        ! Clamp to [0, 1]
        voigt_E_table(i) = max(0.0d0, min(1.0d0, E_val))
      end do

      voigt_table_ready = .true.

    end subroutine voigt_escape_init_table

    !> Voigt profile ratio phi(x)/phi(0) using Humlicek (1982) Region I/II approx.
    !> For the small-a regime (a < 0.01), this simplifies to:
    !>   H(a,x)/H(a,0) ≈ exp(-x^2) + a*sqrt(pi)/x^2  for |x| > ~2
    !> We use the exact Gaussian core + Lorentzian wing decomposition.
    double precision function voigt_profile_ratio(a, x)
      use mod_global_parameters, only: dpi
      implicit none
      double precision, intent(in) :: a, x
      double precision :: gauss_part, lorentz_part, phi_x, phi_0

      ! phi(0) = H(a,0) ≈ 1/sqrt(pi) * (1 + ...) for small a
      phi_0 = 1.0d0 / sqrt(dpi)

      ! For small a: H(a,x) ≈ exp(-x²)/sqrt(pi) for |x| < ~3
      ! and H(a,x) ≈ a/(pi*x²) for |x| >> 1 where Lorentzian dominates
      ! Use additive approximation: H(a,x) ≈ exp(-x²)/sqrt(pi) + a/(pi*(x²+a²))
      gauss_part = exp(-x*x) / sqrt(dpi)
      if(x*x + a*a > 1.0d-30) then
        lorentz_part = a / (dpi * (x*x + a*a))
      else
        lorentz_part = 0.0d0
      end if
      phi_x = gauss_part + lorentz_part

      voigt_profile_ratio = phi_x / phi_0

    end function voigt_profile_ratio

    !> Look up the Voigt escape probability for a given tau.
    !> Uses linear interpolation in log10(tau) space.
    double precision function voigt_escape_lookup(tau)
      implicit none
      double precision, intent(in) :: tau
      double precision :: logtau, frac
      integer :: idx

      if(tau < 1.0d-6) then
        voigt_escape_lookup = 1.0d0
        return
      end if

      logtau = log10(tau)

      if(logtau <= voigt_logtau_min) then
        voigt_escape_lookup = voigt_E_table(1)
        return
      end if

      if(logtau >= voigt_logtau_max) then
        ! Extrapolate with 1/tau scaling from last table entry
        voigt_escape_lookup = voigt_E_table(n_voigt_table) &
          * (10.0d0**voigt_logtau_max) / tau
        return
      end if

      ! Linear interpolation
      frac = (logtau - voigt_logtau_min) / voigt_logtau_step
      idx = int(frac) + 1
      idx = max(1, min(idx, n_voigt_table - 1))
      frac = frac - dble(idx - 1)

      voigt_escape_lookup = voigt_E_table(idx) * (1.0d0 - frac) &
                          + voigt_E_table(idx + 1) * frac

    end function voigt_escape_lookup

    !> Gauss-Legendre quadrature nodes and weights on [a,b].
    !> Uses the Golub-Welsch algorithm for n points.
    subroutine voigt_gauss_legendre(a, b, n, x, w)
      use mod_global_parameters, only: dpi
      implicit none
      double precision, intent(in) :: a, b
      integer, intent(in) :: n
      double precision, intent(out) :: x(n), w(n)
      double precision :: xi, wi, p0, p1, p2, pp, z, z1
      integer :: i, j, k, m

      m = (n + 1) / 2

      do i = 1, m
        ! Initial guess for i-th root
        z = cos(dpi * (dble(i) - 0.25d0) / (dble(n) + 0.5d0))

        ! Newton iteration
        do j = 1, 100
          p0 = 1.0d0
          p1 = 0.0d0
          do k = 1, n
            p2 = p1
            p1 = p0
            p0 = ((2.0d0*dble(k) - 1.0d0) * z * p1 - (dble(k) - 1.0d0) * p2) / dble(k)
          end do
          ! p0 = P_n(z), derivative:
          pp = dble(n) * (z * p0 - p1) / (z*z - 1.0d0)
          z1 = z
          z = z - p0 / pp
          if(abs(z - z1) < 1.0d-15) exit
        end do

        ! Map from [-1,1] to [a,b]
        xi = 0.5d0 * ((b - a) * z + (b + a))
        wi = (b - a) / ((1.0d0 - z*z) * pp*pp)

        x(i) = xi
        w(i) = wi
        x(n + 1 - i) = a + b - xi
        w(n + 1 - i) = wi
      end do

    end subroutine voigt_gauss_legendre

    subroutine radiative_cooling_init(fl,read_params)
      use mod_global_parameters
      interface 
        subroutine read_params(fl)
          use mod_global_parameters, only: unitpar,par_files
          import rc_fluid
          type(rc_fluid), intent(inout) :: fl
    
        end subroutine read_params
      end interface  

      type(rc_fluid), intent(inout) :: fl
  
      double precision, dimension(:), allocatable :: t_table
      double precision, dimension(:), allocatable :: L_table
      double precision, dimension(:), allocatable :: f_table
      double precision :: ratt, fact1, fact2, fact3, dL1, dL2
      double precision :: tstep, Lstep
      integer :: ntable, i, j
      logical :: jump
      Character(len=65) :: PPL_curves(1:6)

      fl%ncool=4000
      fl%coolcurve='JCcorona'
      fl%tlow=bigdouble
      fl%Tfix=.false.
      fl%rc_split=.false.
      fl%rad_suppress_temp=0.0d0
      fl%rad_cut_hgt=0.0d0
      fl%rad_cut_dey=0.15d0
      fl%rad_modify=.false.
      fl%rad_modify_sym=.false.
      fl%rad_taper_rho=bigdouble
      fl%rad_taper_dey=0.0d0
      fl%rad_damp=.false.
      fl%rad_damp_height=0.5d0
      fl%rad_damp_scale=0.15d0
      call read_params(fl)

      ! Build Voigt escape lookup table if needed (once, shared across fluids)
      if(fl%rad_escape_prob .and. fl%rad_escape_type == 'voigt') then
        call voigt_escape_init_table()
        if(mype == 0) then
          write(*,'(A,I0,A,ES9.2,A,F4.1,A)') &
            ' Voigt escape table: ', n_voigt_table, ' points, a_ref=', &
            voigt_a_ref, ', x_max=', voigt_xmax, ' Doppler widths'
        end if
      end if

      if (fl%fip_ > 0) then
        select case (trim(fl%coolcurve))
        case ('Dere_photo', 'Dere_photo_DM')
        case default
          call mpistop("FIP cooling requires coolcurve='Dere_photo' or 'Dere_photo_DM'")
        end select
      end if

      if(fl%rc_split) any_source_split=.true.

      ! Checks if coolcurve is a piecewise power law (PPL)
      PPL_curves = [Character(len=65) :: 'Hildner','FM', 'Rosner', 'Klimchuk','SPEX_DM_rough','SPEX_DM_fine']
      do i=1,size(PPL_curves)
         if (PPL_curves(i)==fl%coolcurve) then
            fl%isPPL = .true.
         end if
      end do

      ! Init for PPL
      if (fl%isPPL) then
         ! Read in tables and create t_PPL, l_PPL, a_PPL
         select case(fl%coolcurve)
         
         case('Hildner')
            if(mype ==0) &
            print *,'Use Hildner (1974) piecewise power law'
            fl%n_PPL = n_Hildner
            allocate(fl%t_PPL(1:fl%n_PPL+1), fl%l_PPL(1:fl%n_PPL+1))
            allocate(fl%a_PPL(1:fl%n_PPL))
            fl%t_PPL(1:fl%n_PPL+1) = t_Hildner(1:n_Hildner+1)
            fl%a_PPL(1:fl%n_PPL) = a_Hildner(1:n_Hildner)
            fl%l_PPL(1:fl%n_PPL) = 10.d0**x_Hildner(1:n_Hildner) * (10.d0**fl%t_PPL(1:fl%n_PPL))**fl%a_PPL(1:fl%n_PPL)        

         case('FM')
            if(mype==0) &
            print *,'Use Forbes and Malherbe (1991)-like piecewise power law'
            fl%n_PPL = n_FM
            allocate(fl%t_PPL(1:fl%n_PPL+1), fl%l_PPL(1:fl%n_PPL+1))
            allocate(fl%a_PPL(1:fl%n_PPL))
            fl%t_PPL(1:fl%n_PPL+1) = t_FM(1:n_FM+1)
            fl%a_PPL(1:fl%n_PPL) = a_FM(1:n_FM)
            fl%l_PPL(1:fl%n_PPL) = 10.d0**x_FM(1:n_FM) * (10.d0**fl%t_PPL(1:fl%n_PPL))**fl%a_PPL(1:fl%n_PPL) 

         case('Rosner')
            if(mype==0) &
            print *,'Use piecewise power law according to Rosner (1978)'
            if(mype ==0) &
            print *,'and extended by Priest (1982) from Van Der Linden (1991)'
            fl%n_PPL = n_Rosner
            allocate(fl%t_PPL(1:fl%n_PPL+1), fl%l_PPL(1:fl%n_PPL+1))
            allocate(fl%a_PPL(1:fl%n_PPL))
            fl%t_PPL(1:fl%n_PPL+1) = t_Rosner(1:n_Rosner+1)
            fl%a_PPL(1:fl%n_PPL) = a_Rosner(1:n_Rosner)
            fl%l_PPL(1:fl%n_PPL) = 10.d0**x_Rosner(1:n_Rosner) * (10.d0**fl%t_PPL(1:fl%n_PPL))**fl%a_PPL(1:fl%n_PPL)  

         case('Klimchuk')
            if(mype==0) &
            print *,'Use Klimchuk (2008) piecewise power law'
            fl%n_PPL = n_Klimchuk
            allocate(fl%t_PPL(1:fl%n_PPL+1), fl%l_PPL(1:fl%n_PPL+1))
            allocate(fl%a_PPL(1:fl%n_PPL))
            fl%t_PPL(1:fl%n_PPL+1) = t_Klimchuk(1:n_Klimchuk+1)
            fl%a_PPL(1:fl%n_PPL) = a_Klimchuk(1:n_Klimchuk)
            fl%l_PPL(1:fl%n_PPL) = 10.d0**x_Klimchuk(1:n_Klimchuk) * (10.d0**fl%t_PPL(1:fl%n_PPL))**fl%a_PPL(1:fl%n_PPL)  

         case('SPEX_DM_rough')
            if(mype==0) &
            print *,'Use the rough piece wise power law fit to the SPEX_DM curve (2009)'
            fl%n_PPL = n_SPEX_DM_rough
            allocate(fl%t_PPL(1:fl%n_PPL+1), fl%l_PPL(1:fl%n_PPL+1))
            allocate(fl%a_PPL(1:fl%n_PPL))
            fl%t_PPL(1:fl%n_PPL+1) = t_SPEX_DM_rough(1:n_SPEX_DM_rough+1)
            fl%a_PPL(1:fl%n_PPL) = a_SPEX_DM_rough(1:n_SPEX_DM_rough)
            fl%l_PPL(1:fl%n_PPL) = 10.d0**x_SPEX_DM_rough(1:n_SPEX_DM_rough) * (10.d0**fl%t_PPL(1:fl%n_PPL))**fl%a_PPL(1:fl%n_PPL)  

         case('SPEX_DM_fine')
            if(mype==0) &
            print *,'Use the fine, detailed piece wise power law fit to the SPEX_DM curve (2009)'
            fl%n_PPL = n_SPEX_DM_fine
            allocate(fl%t_PPL(1:fl%n_PPL+1), fl%l_PPL(1:fl%n_PPL+1))
            allocate(fl%a_PPL(1:fl%n_PPL))
            fl%t_PPL(1:fl%n_PPL+1) = t_SPEX_DM_fine(1:n_SPEX_DM_fine+1)
            fl%a_PPL(1:fl%n_PPL) = a_SPEX_DM_fine(1:n_SPEX_DM_fine)
            fl%l_PPL(1:fl%n_PPL) = 10.d0**x_SPEX_DM_fine(1:n_SPEX_DM_fine) * (10.d0**fl%t_PPL(1:fl%n_PPL))**fl%a_PPL(1:fl%n_PPL) 

         case default
            call mpistop("This piecewise power law is unknown")
         end select

         ! Go from logarithmic to actual values.
         fl%t_PPL(1:fl%n_PPL+1) = 10.d0**fl%t_PPL(1:fl%n_PPL+1)
         ! Change unit of table if SI is used instead of cgs
         if (SI_unit) fl%l_PPL(1:fl%n_PPL) = fl%l_PPL(1:fl%n_PPL) * 10.0d0**(-13)

         ! Make dimensionless
         fl%t_PPL(1:fl%n_PPL+1) = fl%t_PPL(1:fl%n_PPL+1) / unit_temperature
         fl%l_PPL(1:fl%n_PPL) = fl%l_PPL(1:fl%n_PPL) * unit_numberdensity**2 * unit_time / unit_pressure

         ! Set tref en lref
         fl%l_PPL(fl%n_PPL+1) = fl%l_PPL(fl%n_PPL) * ( fl%t_PPL(fl%n_PPL+1) / fl%t_PPL(fl%n_PPL) )**fl%a_PPL(fl%n_PPL)
         fl%lref = fl%l_PPL(fl%n_PPL+1)
         fl%tref = fl%t_PPL(fl%n_PPL+1)

         ! Set tcoolmin and tcoolmax
         fl%tcoolmin = fl%t_PPL(1)
         fl%tcoolmax = fl%t_PPL(fl%n_PPL+1)
         ! smaller value for lowest temperatures from cooling table and user's choice
         if (fl%tlow==bigdouble) fl%tlow=fl%tcoolmin
         !create y_PPL
         call create_y_PPL(fl)

      else

         ! Init for interpolatable tables
         allocate(fl%tcool(1:fl%ncool), fl%Lcool(1:fl%ncool), fl%dLdtcool(1:fl%ncool))
         allocate(fl%Yc(1:fl%ncool))
         if(fl%fip_ > 0) allocate(fl%frac_lowFIP(1:fl%ncool))

         fl%tcool(1:fl%ncool)    = zero
         fl%Lcool(1:fl%ncool)    = zero
         fl%dLdtcool(1:fl%ncool) = zero

         ! Read in the selected cooling curve
         select case(fl%coolcurve)

         case('JCcorona')
            if(mype ==0) &
            print *,'Use Colgan & Feldman (2008) cooling curve'
            if(mype ==0) &
            print *,'This version only till 10000 K, beware for floor T treatment'
            ntable = n_JCcorona
            allocate(t_table(1:ntable))
            allocate(L_table(1:ntable))
            t_table(1:ntable) = t_JCcorona(1:n_JCcorona)
            L_table(1:ntable) = l_JCcorona(1:n_JCcorona)

         case('DM')
            if(mype ==0) &
            print *,'Use Dalgarno & McCray (1972) cooling curve'
            ntable = n_DM
            allocate(t_table(1:ntable))
            allocate(L_table(1:ntable))
            t_table(1:ntable) = t_DM(1:n_DM)
            L_table(1:ntable) = l_DM(1:n_DM)

         case('MB')
            if(mype ==0) &
            write(*,'(3a)') 'Use MacDonald & Bailey (1981) cooling curve '&
                 ,'as implemented in ZEUS-3D, with the values '&
                 ,'from Dalgarno & McCRay (1972) for low temperatures.'
            ntable = n_MB + 20
            allocate(t_table(1:ntable))
            allocate(L_table(1:ntable))
            t_table(1:ntable) = t_DM(1:21)
            L_table(1:ntable) = l_DM(1:21)
            t_table(22:ntable) = t_MB(2:n_MB)
            L_table(22:ntable) = l_MB(2:n_MB)

         case('MLcosmol')
            if(mype ==0) &
            print *,'Use Mellema & Lundqvist (2002) cooling curve '&
                 ,'for zero metallicity '
            ntable = n_MLcosmol
            allocate(t_table(1:ntable))
            allocate(L_table(1:ntable))
            t_table(1:ntable) = t_MLcosmol(1:n_MLcosmol)
            L_table(1:ntable) = l_MLcosmol(1:n_MLcosmol)

         case('MLwc')
            if(mype ==0) &
            print *,'Use Mellema & Lundqvist (2002) cooling curve '&
                 ,'for WC-star metallicity '
            ntable = n_MLwc
            allocate(t_table(1:ntable))
            allocate(L_table(1:ntable))
            t_table(1:ntable) = t_MLwc(1:n_MLwc)
            L_table(1:ntable) = l_MLwc(1:n_MLwc)

         case('MLsolar1')
            if(mype ==0) &
            print *,'Use Mellema & Lundqvist (2002) cooling curve '&
                 ,'for solar metallicity '
            ntable = n_MLsolar1
            allocate(t_table(1:ntable))
            allocate(L_table(1:ntable))
            t_table(1:ntable) = t_MLsolar1(1:n_MLsolar1)
            L_table(1:ntable) = l_MLsolar1(1:n_MLsolar1)

         case('cloudy_ism')
            if(mype ==0) &
            print *,'Use Cloudy based cooling curve '&
                 ,'for ism metallicity '
            ntable = n_cl_ism
            allocate(t_table(1:ntable))
            allocate(L_table(1:ntable))
            t_table(1:ntable) = t_cl_ism(1:n_cl_ism)
            L_table(1:ntable) = l_cl_ism(1:n_cl_ism)

         case('cloudy_solar')
            if(mype ==0) &
            print *,'Use Cloudy based cooling curve '&
                 ,'for solar metallicity '
            ntable = n_cl_solar
            allocate(t_table(1:ntable))
            allocate(L_table(1:ntable))
            t_table(1:ntable) = t_cl_solar(1:n_cl_solar)
            L_table(1:ntable) = l_cl_solar(1:n_cl_solar)

         case('composite_solar')
            if(mype ==0) then
            print *, 'Use composite cooling curve for solar metallicity:'
            print *, '  T > 12 kK: Dere/Colgan/SPEX weighted average'
            print *, '  7-12 kK:   SPEX CIE through the Lya transition'
            print *, '  T < 7 kK:  cloudy_solar fine-structure/molecular'
            end if
            ntable = n_composite
            allocate(t_table(1:ntable))
            allocate(L_table(1:ntable))
            t_table(1:ntable) = t_composite(1:n_composite)
            L_table(1:ntable) = l_composite(1:n_composite)

         case('SPEX')
            if(mype ==0) &
            print *,'Use SPEX cooling curve (Schure et al. 2009) '&
                 ,'for solar metallicity '
            ntable = n_SPEX
            allocate(t_table(1:ntable))
            allocate(L_table(1:ntable))
            t_table(1:ntable) = t_SPEX(1:n_SPEX)
            L_table(1:ntable) = l_SPEX(1:n_SPEX) + log10(nenh_SPEX(1:n_SPEX))
            ! SPEX two-table convention: the absorbed nenh_SPEX factor in
            ! L_table makes Q = n_e * n_H * Lambda give the published rate
            ! ONLY when n_e ~ n_H (the FI assumption). In LTE+ionE mode the
            ! actual Saha n_e is much smaller than n_H at low T and the
            ! formula double-counts the equilibrium factor. The flag below
            ! tells the LTE+ionE code path to substitute n_H for n_e in the
            ! cooling rate, which recovers the correct published rate
            ! Q = n_H^2 * (nenh_eq * Lambda_SPEX) = n_H^2 * Lambda_table.
            fl%lambda_needs_nenh_table = .true.

         case('SPEX_DM')
            if(mype ==0) then
               print *, 'Use SPEX cooling curve for solar metallicity above 10^4 K. '
               print *, 'At lower temperatures,use Dalgarno & McCray (1972), '
               print *, 'with a pre-set ionization fraction of 10^-3. '
               print *, 'as described by Schure et al. (2009). '
            endif
            ntable = n_SPEX + n_DM_2 - 6
            allocate(t_table(1:ntable))
            allocate(L_table(1:ntable))
            t_table(1:n_DM_2-1) = t_DM_2(1:n_DM_2-1)
            L_table(1:n_DM_2-1) = L_DM_2(1:n_DM_2-1)
            t_table(n_DM_2:ntable) = t_SPEX(6:n_SPEX)
            L_table(n_DM_2:ntable) = l_SPEX(6:n_SPEX) + log10(nenh_SPEX(6:n_SPEX))
            ! Same SPEX two-table convention as the pure SPEX case above.
            ! The DM_2 segment was tabulated assuming y_DM = 10^-3 already
            ! built into the published L_DM_2 values, so it follows the
            ! same "Lambda_table includes the nenh factor" convention as
            ! the SPEX segment above (modulo a constant 10^-3 instead of
            ! the SPEX equilibrium ratio). Treating both halves with the
            ! same n_H^2 * Lambda_table formula in LTE+ionE mode gives
            ! consistent rates and removes the double-counting.
            fl%lambda_needs_nenh_table = .true.

         case('Dere_corona')
            if(mype ==0) &
            print *,'Use Dere (2009) cooling curve for solar corona'
            ntable = n_Dere
            allocate(t_table(1:ntable))
            allocate(L_table(1:ntable))
            t_table(1:ntable) = t_Dere(1:n_Dere)
            L_table(1:ntable) = l_Dere_corona(1:n_Dere)

         case('Dere_corona_DM')
            if(mype==0)&
            print *, 'Combination of Dere_corona (2009) for high temperatures and'
            if(mype==0)&
            print *, 'Dalgarno & McCray (1972), DM2, for low temperatures'
            ntable = n_Dere + n_DM_2 - 1 
            allocate(t_table(1:ntable))
            allocate(L_table(1:ntable))
            t_table(1:n_DM_2-1) = t_DM_2(1:n_DM_2-1)
            L_table(1:n_DM_2-1) = L_DM_2(1:n_DM_2-1)
            t_table(n_DM_2:ntable) = t_Dere(1:n_Dere)
            L_table(n_DM_2:ntable) = l_Dere_corona(1:n_Dere)

         case('Dere_photo')
            if(mype ==0) &
            print *,'Use Dere (2009) cooling curve for solar photophere'
            ntable = n_Dere
            allocate(t_table(1:ntable))
            allocate(L_table(1:ntable))
            if (fl%fip_ > 0) allocate(f_table(1:ntable))
            t_table(1:ntable) = t_Dere(1:n_Dere)
            L_table(1:ntable) = l_Dere_photo(1:n_Dere)
            if (fl%fip_ > 0) f_table(1:ntable) = lowFIP_frac(1:n_Dere)

         case('Dere_photo_DM')
            if(mype==0)&
            print *, 'Combination of Dere_photo (2009) for high temperatures and'
            if(mype==0)&
            print *, 'Dalgarno & McCray (1972), DM2, for low temperatures'
            ntable = n_Dere + n_DM_2 - 1 
            allocate(t_table(1:ntable))
            allocate(L_table(1:ntable))
            if (fl%fip_ > 0) allocate(f_table(1:ntable))
            t_table(1:n_DM_2-1) = t_DM_2(1:n_DM_2-1)
            L_table(1:n_DM_2-1) = L_DM_2(1:n_DM_2-1)
            t_table(n_DM_2:ntable) = t_Dere(1:n_Dere)
            L_table(n_DM_2:ntable) = l_Dere_photo(1:n_Dere)
            if (fl%fip_ > 0) then
              f_table(1:n_DM_2-1) = zero
              f_table(n_DM_2:ntable) = lowFIP_frac(1:n_Dere)
            end if

         case('Colgan')
            if(mype==0) &
            print *, 'Use Colgan (2008) cooling curve'
            ntable = n_Colgan
            allocate(t_table(1:ntable))
            allocate(L_table(1:ntable))
            t_table(1:ntable) = t_Colgan(1:n_Colgan)
            L_table(1:ntable) = l_Colgan(1:n_Colgan)

         case('Colgan_DM')
            if(mype==0)&
            print *, 'Combination of Colgan (2008) for high temperatures and'
            if(mype==0)&
            print *, 'Dalgarno & McCray (1972), DM2, for low temperatures'
            ntable = n_Colgan + n_DM_2
            allocate(t_table(1:ntable))
            allocate(L_table(1:ntable))
            t_table(1:n_DM_2) = t_DM_2(1:n_DM_2)
            L_table(1:n_DM_2) = L_DM_2(1:n_DM_2)
            t_table(n_DM_2+1:ntable) = t_Colgan(1:n_Colgan)
            L_table(n_DM_2+1:ntable) = l_Colgan(1:n_Colgan)

         case default
            call mpistop("This coolingcurve is unknown")
         end select


         ! create cooling table(s) for use in amrvac
         fl%tcoolmax = t_table(ntable)
         fl%tcoolmin = t_table(1)
         ratt = (fl%tcoolmax-fl%tcoolmin)/( dble(fl%ncool-1) + smalldouble)

         fl%tcool(1) = fl%tcoolmin
         fl%Lcool(1) = L_table(1)

         fl%tcool(fl%ncool) = fl%tcoolmax
         fl%Lcool(fl%ncool) = L_table(ntable)

         if (fl%fip_ > 0) then
           fl%frac_lowFIP(1) = f_table(1)
           fl%frac_lowFIP(fl%ncool) = f_table(ntable)
         end if

         do i=2,fl%ncool        ! loop to create one table
           fl%tcool(i) = fl%tcool(i-1)+ratt
           do j=1,ntable-1   ! loop to create one spot on a table
           ! Second order polynomial interpolation, except at the outer edge,
           ! or in case of a large jump.
             if(fl%tcool(i) < t_table(j+1)) then
                if(j.eq. ntable-1 )then
                  fact1 = (fl%tcool(i)-t_table(j+1))     &
                        /(t_table(j)-t_table(j+1)) 
                  fact2 = (fl%tcool(i)-t_table(j))       &
                        /(t_table(j+1)-t_table(j)) 
                  fl%Lcool(i) = L_table(j)*fact1 + L_table(j+1)*fact2 
                  if (fl%fip_ > 0) then
                    fl%frac_lowFIP(i) = f_table(j)*fact1 + f_table(j+1)*fact2
                  end if
                  exit
                else 
                  dL1 = L_table(j+1)-L_table(j)
                  dL2 = L_table(j+2)-L_table(j+1)
                  jump =(max(dabs(dL1),dabs(dL2)) > 2*min(dabs(dL1),dabs(dL2)))
                end if
                if( jump ) then
                  fact1 = (fl%tcool(i)-t_table(j+1))     &
                        /(t_table(j)-t_table(j+1)) 
                  fact2 = (fl%tcool(i)-t_table(j))       &
                        /(t_table(j+1)-t_table(j)) 
                  fl%Lcool(i) = L_table(j)*fact1 + L_table(j+1)*fact2
                  if (fl%fip_ > 0) then
                    fl%frac_lowFIP(i) = f_table(j)*fact1 + f_table(j+1)*fact2
                  end if
                  exit          
                else
                  fact1 = ((fl%tcool(i)-t_table(j+1))     &
                        * (fl%tcool(i)-t_table(j+2)))   &
                        / ((t_table(j)-t_table(j+1)) &
                        * (t_table(j)-t_table(j+2)))
                  fact2 = ((fl%tcool(i)-t_table(j))       &
                        * (fl%tcool(i)-t_table(j+2)))   &
                        / ((t_table(j+1)-t_table(j)) &
                        * (t_table(j+1)-t_table(j+2)))
                  fact3 = ((fl%tcool(i)-t_table(j))       &
                        * (fl%tcool(i)-t_table(j+1)))   &
                        / ((t_table(j+2)-t_table(j)) &
                        * (t_table(j+2)-t_table(j+1)))
                  fl%Lcool(i) = L_table(j)*fact1 + L_table(j+1)*fact2 &
                           + L_table(j+2)*fact3
                  if (fl%fip_ > 0) then
                    fl%frac_lowFIP(i) = f_table(j)*fact1 + f_table(j+1)*fact2 &
                                      + f_table(j+2)*fact3
                  end if
                  exit
                end if
             end if
           end do  ! end loop to find create one spot on a table
         end do    ! end loop to create one table

         ! Go from logarithmic to actual values.
         fl%tcool(1:fl%ncool) = 10.0D0**fl%tcool(1:fl%ncool)
         fl%Lcool(1:fl%ncool) = 10.0D0**fl%Lcool(1:fl%ncool)

         ! Change unit of table if SI is used instead of cgs
         if (SI_unit) fl%Lcool(1:fl%ncool) = fl%Lcool(1:fl%ncool) * 10.0d0**(-13)

         ! Scale both T and Lambda
         fl%tcool(1:fl%ncool) = fl%tcool(1:fl%ncool) / unit_temperature
         fl%Lcool(1:fl%ncool) = fl%Lcool(1:fl%ncool) * unit_numberdensity**2 * unit_time / unit_pressure

         fl%tcoolmin       = fl%tcool(1)+smalldouble  ! avoid pointless interpolation
         ! Convert rad_suppress_temp from Kelvin to code units
         if(fl%rad_suppress_temp > 0.0d0) then
           fl%suppress_temp_code = fl%rad_suppress_temp / unit_temperature
           if(mype == 0) then
             write(*,'(A,ES10.3,A)') ' Cooling suppression active: disabled below T = ', &
               fl%rad_suppress_temp, ' K within rad_cut_hgt'
           end if
         end if
         ! smaller value for lowest temperatures from cooling table and user's choice
         if (fl%tlow==bigdouble) fl%tlow=fl%tcoolmin
         fl%tcoolmax       = fl%tcool(fl%ncool)
         fl%lgtcoolmin = dlog10(fl%tcoolmin)
         fl%lgtcoolmax = dlog10(fl%tcoolmax)
         fl%lgstep = (fl%lgtcoolmax-fl%lgtcoolmin) * 1.d0 / (fl%ncool-1)
         fl%dLdtcool(1)     = (fl%Lcool(2)-fl%Lcool(1))/(fl%tcool(2)-fl%tcool(1))
         fl%dLdtcool(fl%ncool) = (fl%Lcool(fl%ncool)-fl%Lcool(fl%ncool-1))/(fl%tcool(fl%ncool)-fl%tcool(fl%ncool-1))

         do i=2,fl%ncool-1
           fl%dLdtcool(i) = (fl%Lcool(i+1)-fl%Lcool(i-1))/(fl%tcool(i+1)-fl%tcool(i-1))
         end do

         deallocate(t_table)
         deallocate(L_table)
         if (allocated(f_table)) deallocate(f_table)

         fl%tref = fl%tcoolmax
         fl%lref = fl%Lcool(fl%ncool)
         fl%Yc(fl%ncool) = zero
         do i=fl%ncool-1, 1, -1
            fl%Yc(i) = fl%Yc(i+1)
            do j=1,100
               tstep = 1.0d-2*(fl%tcool(i+1)-fl%tcool(i))
               call findL(fl%tcool(i+1)-j*tstep, Lstep, fl)
               fl%Yc(i) = fl%Yc(i) + fl%lref/fl%tref*tstep/Lstep
            end do
         end do
      end if

      rc_gamma_1=rc_gamma-1.d0
      invgam = 1.d0/rc_gamma_1

    end subroutine radiative_cooling_init

    subroutine create_y_PPL(fl)
    !  creates the constants of integration needed for solving
    !  the cooling law exact for a piecewise power law 
    !  In correspondence with eq. A6 of Townsend (2009)
      use mod_global_parameters
      type(rc_fluid) :: fl
      double precision :: y_extra, factor
      integer :: i

      allocate(fl%y_PPL(1:fl%n_PPL+1))

      fl%y_PPL(1:fl%n_PPL+1) = zero

      do i=fl%n_PPL, 1, -1 
          factor = fl%l_PPL(fl%n_PPL+1) * fl%t_PPL(i) / (fl%l_PPL(i) * fl%t_PPL(fl%n_PPL+1))
          if (fl%a_PPL(i) == 1.d0) then 
             y_extra =  log( fl%t_PPL(i) / fl%t_PPL(i+1) )
          else
             y_extra = 1 / (1 - fl%a_PPL(i)) * (1 - ( fl%t_PPL(i) / fl%t_PPL(i+1) )**(fl%a_PPL(i)-1) )
          end if
          fl%y_PPL(i) = fl%y_PPL(i+1) - factor*y_extra
      end do
    end subroutine create_y_PPL

    subroutine getvar_cooling(ixI^L,ixO^L,w,x,coolrate,fl)
    ! Create extra variable to show cooling rate in the output
    ! Uses a simple explicit scheme.
    ! N.B. Since there is no knowledge of the timestep size,
    ! there is no upper limit for the cooling rate.
      use mod_global_parameters
      
      integer, intent(in)          :: ixI^L,ixO^L
      double precision, intent(in) :: x(ixI^S,1:ndim)
      double precision             :: w(ixI^S,1:nw)
      double precision, intent(out):: coolrate(ixI^S)
      type(rc_fluid), intent(in) :: fl

      double precision :: pth(ixI^S),rho(ixI^S)
      double precision :: L1,Te(ixI^S),Rfactor(ixI^S)
      double precision :: ne(ixI^S), nH_arr(ixI^S), rho2_factor(ixI^S)
      double precision :: taper
      integer :: ix^D

      ! call fl%get_pthermal(w,x,ixI^L,ixO^L,pth)
      call fl%get_rho(w,x,ixI^L,ixO^L,rho)
      ! call fl%get_var_Rfactor(w,x,ixI^L,ixO^L,Rfactor)
      ! Te(ixO^S) = pth(ixO^S) / (rho(ixO^S)*Rfactor(ixO^S))
      call fl%get_Te(w,x,ixI^L,ixO^L,Te)
      call fl%get_ne_nH(ixI^L, ixO^L, w, x, ne, nH_arr)
      call radiative_cooling_rho2_factor(ixI^L,ixO^L,w,x,fl,rho2_factor)

      {do ix^DB = ixO^LIM^DB\}
         ! Determine explicit cooling
         if(Te(ix^D) <= fl%tcoolmin) then
           L1 = zero
         else if(Te(ix^D) >= fl%tcoolmax)then
           call calc_l_extended(Te(ix^D),L1,fl)
           L1 = L1*ne(ix^D)*nH_arr(ix^D)
         else
           call findL(Te(ix^D),L1,fl)
           L1 = L1*ne(ix^D)*nH_arr(ix^D)
         end if
         if(slab_uniform .and. fl%rad_damp .and. x(ix^D,ndim) .le. fl%rad_damp_height) then
           L1 = L1*exp(-(x(ix^D,ndim)-fl%rad_damp_height)**2/fl%rad_damp_scale**2)
         end if
         call radiative_cooling_taper(ix^D, x(ix^D,ndim), rho(ix^D), Te(ix^D), fl, taper)
         L1 = L1 * taper
         coolrate(ix^D) = L1*rho2_factor(ix^D)
      {end do\}
    end subroutine getvar_cooling

    subroutine getvar_cooling_exact(qdt, ixI^L, ixO^L, wCT, w, x, coolrate, fl)
    ! Calculates cooling rate using the exact cooling method,
      use mod_global_parameters

      integer, intent(in)           :: ixI^L, ixO^L
      double precision, intent(in)  :: qdt, x(ixI^S, 1:ndim), wCT(ixI^S, 1:nw)
      double precision              :: w(ixI^S, 1:nw)
      double precision, intent(out) :: coolrate(ixI^S)
      type(rc_fluid), intent(in)   :: fl
      double precision              :: y1, y2, l1, tlocal2
      double precision              :: Te(ixI^S), pnew(ixI^S), rho(ixI^S), rhonew(ixI^S)
      double precision              :: emin, Lmax, fact, Rfactor(ixI^S), pth(ixI^S)
      double precision              :: ne(ixI^S), nH_arr(ixI^S), rho2_factor(ixI^S)
      double precision              :: taper
      ! LTE+IonE variables
      double precision              :: nH_val, log_nH, log_p_nH
      double precision              :: eint_current
      double precision              :: y_l, T_l
      integer                       :: ix^D

      call fl%get_pthermal(wCT, x, ixI^L, ixO^L, pth)
      call fl%get_rho(wCT, x, ixI^L, ixO^L, rho)
      call fl%get_var_Rfactor(wCT,x,ixI^L,ixO^L,Rfactor)
      call fl%get_Te(wCT, x, ixI^L, ixO^L, Te)
      call fl%get_ne_nH(ixI^L, ixO^L, wCT, x, ne, nH_arr)
      call radiative_cooling_rho2_factor(ixI^L,ixO^L,wCT,x,fl,rho2_factor)
      ! Te(ixO^S)=pth(ixO^S)/(rho(ixO^S)*Rfactor(ixO^S))

      call fl%get_pthermal(w, x, ixI^L, ixO^L, pnew)
      call fl%get_rho(w, x, ixI^L, ixO^L, rhonew)

      fact=fl%lref*qdt/fl%tref

      {do ix^DB = ixO^LIM^DB\}
         emin = rhonew(ix^D) * fl%tlow * Rfactor(ix^D) * invgam
         if (fl%ionE) then
           nH_val = rhonew(ix^D) / fl%nH2rhoFactor
           log_nH = dlog10(nH_val)
           if (fl%method == 'analytic') then
             T_l = Te(ix^D)
             y_l = wCT(ix^D, iw_ne) / nH_val
             eint_current = fl%inv_gamma_minus_1 * (1.0d0 + y_l) * nH_val * T_l &
                 + y_l * fl%eion_per_nH * nH_val
           else
             log_p_nH = dlog10(pnew(ix^D) / nH_val)
             eint_current = pnew(ix^D) * fl%p2eint(log_nH, log_p_nH)
           end if
           lmax = max(zero, (eint_current - emin) / qdt)
         else
           lmax = max(zero, ( pnew(ix^D)*invgam - emin ) / qdt)
         end if

         ! No cooling if temperature is below floor level.
         ! Assuming Bremsstrahlung if temperature is higher than maximum.
         if( Te(ix^D)<= fl%tcoolmin) then
           l1 = zero
         else if( Te(ix^D)>= fl%tcoolmax ) then
           call calc_l_extended(Te(ix^D), l1, fl)
           if (fl%lambda_needs_nenh_table) then
             l1 = l1 * nH_arr(ix^D) * nH_arr(ix^D)
           else
             l1 = l1 * ne(ix^D) * nH_arr(ix^D)
           end if
           l1 = min(l1*rho2_factor(ix^D), lmax)
         else
           !> Always classical Townsend first. Upgrade to Y_mod only where
           !> ionisation buffering matters (large ΔT, recombination zone).
           call findY(Te(ix^D), y1, fl)
           if (fl%lambda_needs_nenh_table) then
             y2 = y1 + rho2_factor(ix^D)*fact * nH_arr(ix^D) * nH_arr(ix^D) * rc_gamma_1 &
                              / (rho(ix^D) * Rfactor(ix^D))
           else
             y2 = y1 + rho2_factor(ix^D)*fact * ne(ix^D) * nH_arr(ix^D) * rc_gamma_1 &
                              / (rho(ix^D) * Rfactor(ix^D))
           end if
           call findT(tlocal2, y2, fl)

           if (fl%ionE .and. fl%Y_mod_built .and. &
               dabs(Te(ix^D) - tlocal2) > 1.0d-4 * Te(ix^D)) then
             y1 = findY_mod(Te(ix^D), nH_arr(ix^D), fl)
             if (y1 == y1 .and. abs(y1) < huge(1.0d0)) then
               y2 = y1 + rho2_factor(ix^D)*qdt
               tlocal2 = findT_mod(y2, nH_arr(ix^D), fl)
             end if
           end if

           if( tlocal2 <= fl%tcoolmin ) then
             l1 = lmax
           else if (fl%ionE .and. &
                    dabs(Te(ix^D) - tlocal2) > 1.0d-4 * Te(ix^D)) then
             !> Recombination zone: table-based de for variable c_V.
             l1 = ((fl%eint_from_T(log_nH, dlog10(Te(ix^D))) &
                  - fl%eint_from_T(log_nH, dlog10(tlocal2))) * nH_val) / qdt
             l1 = max(l1, zero)
           else
             !> Saturated y or non-ionE: identical Townsend kinetic form.
             l1 = (Te(ix^D)- tlocal2)*rho(ix^D)*Rfactor(ix^D)*invgam/qdt
           end if
           l1 = min(l1, lmax)
         end if
         call radiative_cooling_taper(ix^D, x(ix^D,ndim), rho(ix^D), Te(ix^D), fl, taper)
         l1 = l1 * taper
         if(slab_uniform .and. fl%rad_damp .and. x(ix^D,ndim) .le. fl%rad_damp_height) then
           l1 = l1*exp(-(x(ix^D,ndim)-fl%rad_damp_height)**2/fl%rad_damp_scale**2)
         end if
        coolrate(ix^D) = l1
      {end do\}
    end subroutine getvar_cooling_exact

    subroutine radiative_cooling_add_source(qdt,ixI^L,ixO^L,wCT,wCTprim,w,x,&
         qsourcesplit,active,fl)
    ! w[iw]=w[iw]+qdt*S[wCT,x] where S is the source based on wCT within ixO
      use mod_global_parameters
      integer, intent(in) :: ixI^L, ixO^L
      double precision, intent(in) :: qdt, x(ixI^S,1:ndim), wCT(ixI^S,1:nw), wCTprim(ixI^S,1:nw)
      double precision, intent(inout) :: w(ixI^S,1:nw)
      logical, intent(in) :: qsourcesplit
      logical, intent(inout) :: active
      type(rc_fluid), intent(in) :: fl
      double precision, allocatable, dimension(:^D&) :: Lequi
      double precision :: lb_t0_cool

      if (lb_diagnose) lb_t0_cool = MPI_WTIME()
      if(qsourcesplit .eqv.fl%rc_split) then
       active = .true.
       call cool_exact(qdt,ixI^L,ixO^L,wCT,wCTprim,w,x,fl)
       if(fl%subtract_equi) then
          allocate(Lequi(ixI^S))
          call get_cool_equi(qdt,ixI^L,ixO^L,wCT,w,x,fl,Lequi)
          w(ixO^S,fl%e_) = w(ixO^S,fl%e_)+Lequi(ixO^S)
          deallocate(Lequi)
       endif
       if( fl%Tfix ) call floortemperature(qdt,ixI^L,ixO^L,wCT,w,x,fl)
      end if
      if (lb_diagnose) lb_cool_accum = lb_cool_accum + (MPI_WTIME() - lb_t0_cool)
    end subroutine radiative_cooling_add_source

    subroutine floortemperature(qdt,ixI^L,ixO^L,wCT,w,x,fl) !> This will need revisiting in LTE
    !  Force minimum temperature to a fixed temperature
      use mod_global_parameters
      integer, intent(in)             :: ixI^L, ixO^L
      double precision, intent(in)    :: qdt, x(ixI^S,1:ndim), wCT(ixI^S,1:nw)
      double precision, intent(inout) :: w(ixI^S,1:nw)
      type(rc_fluid), intent(in) :: fl
      double precision :: etherm(ixI^S), rho(ixI^S), Rfactor(ixI^S),emin
      integer :: ix^D

      call fl%get_pthermal(w,x,ixI^L,ixO^L,etherm)  
      call fl%get_rho(w,x,ixI^L,ixO^L,rho)  
      call fl%get_var_Rfactor(wCT,x,ixI^L,ixO^L,Rfactor)
      {do ix^DB = ixO^LIM^DB\}
         emin = rho(ix^D)*fl%tlow*Rfactor(ix^D)
         if(etherm(ix^D) < emin) then
           w(ix^D,fl%e_)=w(ix^D,fl%e_)+(emin-etherm(ix^D))*invgam
         end if
      {end do\}
    end subroutine floortemperature

    subroutine radiative_cooling_taper(ix^D, x_ndim, rho_val, Te_val, fl, factor)
      !> Compute multiplicative taper factor for radiative cooling.
      !> Returns 1.0 when no tapering applies; < 1.0 near boundaries,
      !> at high density, or in optically thick regions (escape probability).
      use mod_global_parameters
      integer, intent(in) :: ix^D
      double precision, intent(in) :: x_ndim, rho_val, Te_val
      type(rc_fluid), intent(in) :: fl
      double precision, intent(out) :: factor
      double precision :: d_boundary, tau, kappa_local

      factor = 1.0d0

      ! Spatial + density taper
      if(slab_uniform .and. fl%rad_modify) then
        ! Spatial taper: distance from nearest relevant boundary
        if(fl%rad_modify_sym) then
          d_boundary = min(x_ndim - xprobmin^ND, xprobmax^ND - x_ndim)
        else
          d_boundary = x_ndim - xprobmin^ND
        end if
        if(d_boundary .le. fl%rad_cut_hgt) then
          if(fl%suppress_temp_code > 0.0d0) then
            ! Temperature suppression: kill cooling for T below threshold.
            ! Coronal cooling above threshold proceeds unmodified.
            if(Te_val .lt. fl%suppress_temp_code) then
              factor = 0.0d0
              return
            end if
          else
            ! Standard: Gaussian taper on ALL cooling near boundary
            factor = factor * exp(-((d_boundary - fl%rad_cut_hgt) / fl%rad_cut_dey)**2)
          end if
        end if

        ! Density taper
        if(rho_val .gt. fl%rad_taper_rho) then
          factor = factor * exp(-((rho_val - fl%rad_taper_rho) / fl%rad_taper_dey)**2)
        end if
      end if

      ! Escape probability: cooling suppression by optical depth
      ! kappa(T) = kappa_0 / (1 + (T/T_cutoff)^alpha) — sigmoid cutoff
      if(fl%rad_escape_prob .and. fl%iw_colmass_ > 0) then
        kappa_local = fl%rad_kappa_eff
        if(fl%rad_kappa_Tcutoff > 0.0d0) then
          kappa_local = kappa_local &
            / (1.0d0 + (Te_val / fl%rad_kappa_Tcutoff)**fl%rad_kappa_alpha)
        end if
        tau = kappa_local * block%wextra(ix^D, fl%iw_colmass_)
        if(tau > 1.0d-6) then
          select case(fl%rad_escape_type)
          case('slab')
            ! Plane-parallel slab: beta(tau) = (1 - exp(-tau))/tau
            factor = factor * (1.0d0 - exp(-tau)) / tau
          case('voigt')
            ! Frequency-integrated escape from a truncated Voigt profile
            ! (CL12 Sec 2.1).  Precomputed lookup table; see voigt_escape_init_table.
            factor = factor * voigt_escape_lookup(tau)
          case default
            call mpistop("Unknown rad_escape_type: use 'slab' or 'voigt'")
          end select
          ! Exponential cutoff at large tau: kills residual cooling
          ! where rho^2 outpaces the escape function decay
          if(fl%rad_escape_tau_cutoff > 0.0d0) then
            factor = factor * exp(-tau / fl%rad_escape_tau_cutoff)
          end if
        end if
      end if
    end subroutine radiative_cooling_taper

    subroutine get_cool_equi(qdt,ixI^L,ixO^L,wCT,w,x,fl,res)
      use mod_global_parameters

      integer, intent(in)             :: ixI^L, ixO^L
      double precision, intent(in)    :: qdt, x(ixI^S,1:ndim), wCT(ixI^S,1:nw)
      double precision, intent(inout) :: w(ixI^S,1:nw)
      type(rc_fluid), intent(in) :: fl
      double precision, intent(out) :: res(ixI^S)

      double precision :: pth(ixI^S),rho(ixI^S),Rfactor(ixI^S),L1,Tlocal2
      double precision :: Te(ixI^S)
      double precision :: emin, Lmax
      double precision :: Y1, Y2
      double precision :: de, emax,fact
      double precision :: ne(ixI^S), nH_arr(ixI^S), rho2_factor(ixI^S)
      double precision :: taper
      ! LTE+IonE variables
      double precision :: nH_val, log_nH, log_p_nH
      double precision :: eint_current
      double precision :: y_l, T_l
      integer :: ix^D

      call fl%get_pthermal_equi(wCT,x,ixI^L,ixO^L,pth)
      call fl%get_rho_equi(wCT,x,ixI^L,ixO^L,rho)
      call fl%get_var_Rfactor(wCT,x,ixI^L,ixO^L,Rfactor)
      ! Te(ixO^S)=pth(ixO^S)/(rho(ixO^S)*Rfactor(ixO^S))
      ! temperature and densities of the background alone, to match pth and rho
      call fl%get_temperature_equi(wCT,x,ixI^L,ixO^L,Te)
      call fl%get_ne_nH_equi(ixI^L, ixO^L, wCT, x, ne, nH_arr)
      call radiative_cooling_rho2_factor(ixI^L,ixO^L,wCT,x,fl,rho2_factor)

      res=0d0

      fact = fl%lref*qdt/fl%tref
      {do ix^DB = ixO^LIM^DB\}
           emin = rho(ix^D)*fl%tlow*Rfactor(ix^D)*invgam
           if (fl%ionE) then
             nH_val = rho(ix^D) / fl%nH2rhoFactor
             log_nH = dlog10(nH_val)
             if (fl%method == 'analytic') then
               T_l = Te(ix^D)
               y_l = wCT(ix^D, iw_ne) / nH_val
               eint_current = 1.5d0 * (1.0d0 + y_l) * nH_val * T_l &
                   + y_l * fl%eion_per_nH * nH_val
             else
               log_p_nH = dlog10(pth(ix^D) / nH_val)
               eint_current = pth(ix^D) * fl%p2eint(log_nH, log_p_nH)
             end if
             Lmax = max(zero, (eint_current - emin) / qdt)
             emax = max(zero, eint_current - emin)
           else
             Lmax = max(zero,(pth(ix^D)*invgam-emin)/qdt)
             emax = max(zero, pth(ix^D)*invgam-emin)
           end if
           !  Determine explicit cooling
           !  If temperature is below floor level, no cooling.
           !  Stop wasting time and go to next gridpoint.
           !  If the temperature is higher than the maximum,
           !  assume Bremsstrahlung
           if( Te(ix^D)<=fl%tcoolmin ) then
             ! res already initialised to 0d0 above; no cooling
           else if( Te(ix^D)>=fl%tcoolmax )then
             call calc_l_extended(Te(ix^D), L1,fl)
             if (fl%lambda_needs_nenh_table) then
               L1 = L1 * nH_arr(ix^D) * nH_arr(ix^D)
             else
               L1 = L1 * ne(ix^D) * nH_arr(ix^D)
             end if
             if(phys_trac) then
               if(Te(ix^D)<block%wextra(ix^D,fl%Tcoff_)) then
                 L1=L1*sqrt((Te(ix^D)/block%wextra(ix^D,fl%Tcoff_))**5)
               end if
             end if
             L1 = min(rho2_factor(ix^D)*L1,Lmax)
             res(ix^D) = L1*qdt
           else
             !> Always classical Townsend first. Upgrade to Y_mod only in the
             !> recombination zone (large ΔT) — saturated y uses the identical
             !> formula as ionE=false, removing per-substep asymmetry.
             call findY(Te(ix^D),Y1,fl)
             if (fl%lambda_needs_nenh_table) then
               Y2 = Y1 + rho2_factor(ix^D)*fact * nH_arr(ix^D) * nH_arr(ix^D) * rc_gamma_1 &
                                / (rho(ix^D) * Rfactor(ix^D))
             else
               Y2 = Y1 + rho2_factor(ix^D)*fact * ne(ix^D) * nH_arr(ix^D) * rc_gamma_1 &
                                / (rho(ix^D) * Rfactor(ix^D))
             end if
             call findT(Tlocal2,Y2,fl)

             if (fl%ionE .and. fl%Y_mod_built .and. &
                 dabs(Te(ix^D) - Tlocal2) > 1.0d-4 * Te(ix^D)) then
               Y1 = findY_mod(Te(ix^D), nH_arr(ix^D), fl)
               if (Y1 == Y1 .and. abs(Y1) < huge(1.0d0)) then
                 Y2 = Y1 + rho2_factor(ix^D)*qdt
                 Tlocal2 = findT_mod(Y2, nH_arr(ix^D), fl)
               end if
             end if

             if(Tlocal2<=fl%tcoolmin) then
               de = emax
             else if (fl%ionE .and. &
                      dabs(Te(ix^D) - Tlocal2) > 1.0d-4 * Te(ix^D)) then
               de = (fl%eint_from_T(log_nH, dlog10(Te(ix^D))) &
                    - fl%eint_from_T(log_nH, dlog10(Tlocal2))) * nH_val
               de = max(de, zero)
             else
               de = (Te(ix^D)-Tlocal2)*rho(ix^D)*Rfactor(ix^D)*invgam
             end if
             if(phys_trac) then
               if(Te(ix^D)<block%wextra(ix^D,fl%Tcoff_)) then
                 de=de*sqrt((Te(ix^D)/block%wextra(ix^D,fl%Tcoff_))**5)
               end if
             end if
             de = min(de,emax)
             res(ix^D) = de
           end if
           call radiative_cooling_taper(ix^D, x(ix^D,ndim), rho(ix^D), Te(ix^D), fl, taper)
           res(ix^D) = res(ix^D) * taper
        {end do\}
    end subroutine get_cool_equi

    subroutine cool_exact(qdt,ixI^L,ixO^L,wCT,wCTprim,w,x,fl)
    !  Cooling routine using exact integration method from Townsend 2009
      use mod_global_parameters
      use mod_physics, only: phys_get_ei
      integer, intent(in)             :: ixI^L, ixO^L
      double precision, intent(in)    :: qdt, x(ixI^S,1:ndim), wCT(ixI^S,1:nw), wCTprim(ixI^S,1:nw)
      double precision, intent(inout) :: w(ixI^S,1:nw)
      type(rc_fluid), intent(in) :: fl
      double precision :: Y1, Y2
      double precision :: L1, pth(ixI^S), Tlocal2, pnew(ixI^S)
      double precision :: rho(ixI^S), Te(ixI^S), rhonew(ixI^S), Rfactor(ixI^S)
      double precision :: emin, Lmax, fact
      double precision :: de, emax
      double precision :: ne(ixI^S), nH_arr(ixI^S), rho2_factor(ixI^S)
      double precision :: taper
      ! LTE+IonE variables
      double precision :: nH_val, log_nH, log_p_nH
      double precision :: eint_current
      double precision :: eint_w(ixI^S)  ! actual internal energy from conserved state
      double precision :: de_thin, de_thick, emax_rem
      double precision :: T1, T2, p1(ixI^S), tau, xi
      double precision :: xi_arr(ixI^S), emax_rem_arr(ixI^S)
      double precision :: cool_fac, fip_prim, frac_lowFIP, fip_factor
      double precision :: y_loc, T_loc
      integer :: ix^D

      call fl%get_rho(wCT,x,ixI^L,ixO^L,rho)
      call fl%get_var_Rfactor(wCT,x,ixI^L,ixO^L,Rfactor)
      call fl%get_Te(wCT,x,ixI^L,ixO^L,Te)
      call fl%get_ne_nH(ixI^L, ixO^L, wCT, x, ne, nH_arr)
      call radiative_cooling_rho2_factor(ixI^L,ixO^L,wCT,x,fl,rho2_factor)
      call fl%get_pthermal(w,x,ixI^L,ixO^L,pnew)
      call fl%get_rho(w,x,ixI^L,ixO^L,rhonew)
      if (fl%ionE) eint_w(ixO^S) = phys_get_ei(w, ixI^L, ixO^L)

      fact = fl%lref*qdt/fl%tref

      xi_arr = one
      emax_rem_arr = zero
      {do ix^DB = ixO^LIM^DB\}
         ! Energy floor: always use FI formula (generous safety margin at low T
         ! where ionE floor would be ~2x lower due to neutral vs FI mean mol. weight)
         emin = rhonew(ix^D)*fl%tlow*Rfactor(ix^D)*invgam
         if (fl%ionE) then
           ! LTE+IonE: EoS quantities for Y-advance; cap from conserved state
           nH_val = rhonew(ix^D) / fl%nH2rhoFactor
           log_nH = dlog10(nH_val)
           if (fl%method == 'analytic') then
             ! Use cached Te_ and Ne_ to compute eint directly
             T_loc = Te(ix^D)
             y_loc = wCT(ix^D, iw_ne) / nH_val
             eint_current = fl%inv_gamma_minus_1 * (1.0d0 + y_loc) * nH_val * T_loc &
                 + y_loc * fl%eion_per_nH * nH_val
           else
             log_p_nH = dlog10(pnew(ix^D) / nH_val)
             eint_current = pnew(ix^D) * fl%p2eint(log_nH, log_p_nH)
           end if
           Lmax = max(zero, eint_w(ix^D) - emin) / qdt
           emax = max(zero, eint_w(ix^D) - emin)
         else
           Lmax = max(zero,pnew(ix^D)*invgam-emin)/qdt
           emax = max(zero,pnew(ix^D)*invgam-emin)
         end if

         ! Skip cells below cooling floor: no de_thin contribution, no rad_newton thick step queued.
         if (Te(ix^D) <= fl%tcoolmin) cycle

         ! Multiplicative cool_fac (upstream): xi (optically thick Newton) * FIP * geometric damping.
         ! Multiplied INTO the integration (Y2 increment / extended-Bremsstrahlung L1) so the
         ! Townsend EI mapping sees the reduced Lambda_eff = cool_fac * Lambda.
         if (fl%rad_newton) then
           xi = exp(-pnew(ix^D) / fl%rad_newton_pthick)
           xi = min(max(xi, zero), one)
         else
           xi = one
         end if
         cool_fac = xi*rho2_factor(ix^D)

         if (fl%fip_ > 0) then
           fip_prim = min(maxfip, max(minfip, wCTprim(ix^D,fl%fip_)))
           ! frac_lowFIP(T) in [0,1]: low-FIP contribution to total Lambda at T
           frac_lowFIP = lowFIP_fraction(Te(ix^D), fl)
           fip_factor = one - frac_lowFIP + fip_prim * frac_lowFIP
           cool_fac = cool_fac * fip_factor
         end if

         ! Geometric (Gaussian) damping near lower boundary (chromosphere/photosphere).
         if (slab_uniform .and. fl%rad_damp .and. x(ix^D,ndim) <= xprobmin1 + fl%rad_damp_height) then
           cool_fac = cool_fac * exp(-(x(ix^D,ndim)-xprobmin1-fl%rad_damp_height)**2/fl%rad_damp_scale**2)
         end if
         {^IFONED
         if (slab_uniform .and. fl%rad_damp .and. x(ix^D,ndim) >= xprobmax1 - fl%rad_damp_height) then
           cool_fac = cool_fac * exp(-(x(ix^D,ndim)-xprobmax1+fl%rad_damp_height)**2/fl%rad_damp_scale**2)
         end if
         }

         if( Te(ix^D)>=fl%tcoolmax )then
           call calc_l_extended(Te(ix^D), L1,fl)
           if (fl%lambda_needs_nenh_table) then
             L1 = L1 * nH_arr(ix^D) * nH_arr(ix^D)
           else
             L1 = L1 * ne(ix^D) * nH_arr(ix^D)
           end if
           L1 = cool_fac * L1
           if(phys_trac) then
             if(Te(ix^D)<block%wextra(ix^D,fl%Tcoff_)) then
               L1=L1*sqrt((Te(ix^D)/block%wextra(ix^D,fl%Tcoff_))**5)
             end if
           end if
           L1 = min(L1,Lmax)
           call radiative_cooling_taper(ix^D, x(ix^D,ndim), rho(ix^D), Te(ix^D), fl, taper)
           L1 = L1 * taper
           de_thin = L1 * qdt
           w(ix^D,fl%e_) = w(ix^D,fl%e_) - de_thin
         else
           !> Always compute the CLASSICAL Townsend advance first. This gives
           !> Tlocal2 using the same path as ionE=false. If the resulting dT
           !> is small (saturated y region), the classical kinetic form is
           !> exact and we're done. If dT is large (recombination zone), we
           !> redo the advance with the Y_mod path to capture variable-c_V
           !> dynamics from ionisation buffering.
           call findY(Te(ix^D),Y1,fl)
           if (fl%lambda_needs_nenh_table) then
             Y2 = Y1 + cool_fac * fact * nH_arr(ix^D) * nH_arr(ix^D) * rc_gamma_1 &
                              / (rho(ix^D) * Rfactor(ix^D))
           else
             Y2 = Y1 + cool_fac * fact * ne(ix^D) * nH_arr(ix^D) * rc_gamma_1 &
                              / (rho(ix^D) * Rfactor(ix^D))
           end if
           call findT(Tlocal2,Y2,fl)

           !> Upgrade to Y_mod only in the recombination zone (large dT).
           !> Guard: if findY_mod overflows (cooling curve near-zero at this T),
           !> keep the classical result -- variable c_V is irrelevant where Lambda ~ 0.
           if (fl%ionE .and. fl%Y_mod_built .and. &
               dabs(Te(ix^D) - Tlocal2) > 1.0d-4 * Te(ix^D)) then
             Y1 = findY_mod(Te(ix^D), nH_arr(ix^D), fl)
             if (Y1 == Y1 .and. abs(Y1) < huge(1.0d0)) then
               Y2 = Y1 + cool_fac * qdt
               Tlocal2 = findT_mod(Y2, nH_arr(ix^D), fl)
             end if
           end if

           if(Tlocal2<=fl%tcoolmin) then
             de = emax
           else if (fl%ionE .and. &
                    dabs(Te(ix^D) - Tlocal2) > 1.0d-4 * Te(ix^D)) then
             !> Recombination zone: use table-based de for variable c_V.
             de = (fl%eint_from_T(log_nH, dlog10(Te(ix^D))) &
                   - fl%eint_from_T(log_nH, dlog10(Tlocal2))) * nH_val
             de = max(de, zero)
           else
             !> Saturated y or non-ionE: classical Townsend kinetic form.
             de = (Te(ix^D)-Tlocal2)*rho(ix^D)*Rfactor(ix^D)*invgam
           end if
           if(phys_trac) then
             if(Te(ix^D)<block%wextra(ix^D,fl%Tcoff_)) then
               de=de*sqrt((Te(ix^D)/block%wextra(ix^D,fl%Tcoff_))**5)
             end if
           end if
           de = min(de,emax)
           call radiative_cooling_taper(ix^D, x(ix^D,ndim), rho(ix^D), Te(ix^D), fl, taper)
           de = de * taper
           de_thin = de
           w(ix^D,fl%e_) = w(ix^D,fl%e_) - de_thin
         end if

         ! Queue remaining energy budget for the optically thick step (upstream rad_newton).
         if (fl%rad_newton) then
           xi_arr(ix^D) = xi
           emax_rem_arr(ix^D) = max(zero, emax - de_thin)
         end if
      {end do\}

      ! ------- (B) OPTICALLY-THICK (NEWTON) PART --------
      if (fl%rad_newton) then
      call fl%get_pthermal(w, x, ixI^L, ixO^L, p1)
      {do ix^DB = ixO^LIM^DB\}
          T1 = p1(ix^D) / (rho(ix^D) * Rfactor(ix^D))
          tau = max(0.1d0 * sqrt( fl%rad_newton_rhosurf / rho(ix^D)), 4.d0 * qdt)
          T2 = fl%rad_newton_trad + (T1 - fl%rad_newton_trad) * exp(-qdt / tau)
          de_thick = min((one - xi_arr(ix^D)) * (T1 - T2) * rho(ix^D) * Rfactor(ix^D) * invgam, emax_rem_arr(ix^D))
          w(ix^D,fl%e_) = w(ix^D,fl%e_) - de_thick
      {end do\}
      end if
    end subroutine cool_exact

    subroutine calc_l_extended (tpoint, lpoint,fl)
    !  Calculate l for t beyond tcoolmax
    !  Assumes Bremsstrahlung for the interpolated tables
    !  Uses the power law for piecewise power laws
      double precision, intent(IN)  :: tpoint
      double precision, intent(OUT) :: lpoint
      type(rc_fluid), intent(in) :: fl

      if(fl%isPPL) then
        lpoint =fl%l_PPL(fl%n_PPL) * ( tpoint / fl%t_PPL(fl%n_PPL) )**fl%a_PPL(fl%n_PPL)
      else
        lpoint = fl%Lcool(fl%ncool) * sqrt( tpoint / fl%tcoolmax)
      end if
    end subroutine calc_l_extended

    double precision function lowFIP_fraction(tpoint, fl)
      use mod_global_parameters
  
      double precision, intent(in) :: tpoint
      type(rc_fluid), intent(in)   :: fl
  
      double precision :: lgtp
      integer :: jl
  
      if (tpoint <= fl%tcool(1)) then
        lowFIP_fraction = fl%frac_lowFIP(1)
        return
      else if (tpoint >= fl%tcool(fl%ncool)) then
        lowFIP_fraction = fl%frac_lowFIP(fl%ncool)
        return
      end if
  
      lgtp = dlog10(tpoint)
      jl   = int((lgtp - fl%lgtcoolmin) / fl%lgstep) + 1
      jl   = max(1, min(fl%ncool-1, jl))
  
      lowFIP_fraction = fl%frac_lowFIP(jl) &
           + (tpoint - fl%tcool(jl)) &
           * (fl%frac_lowFIP(jl+1) - fl%frac_lowFIP(jl)) &
           / (fl%tcool(jl+1) - fl%tcool(jl))
    end function lowFIP_fraction

    subroutine findL (tpoint,Lpoint,fl)
    !  Fast search option to find correct point 
    !  in cooling curve
      use mod_global_parameters

      double precision,intent(IN)   :: tpoint
      double precision, intent(OUT) :: Lpoint
      type(rc_fluid), intent(in) :: fl

      double precision :: lgtp
      integer :: jl,i

      if(fl%isPPL) then
        i = maxloc(fl%t_PPL, dim=1, mask=fl%t_PPL<tpoint)
        Lpoint = fl%l_PPL(i) * (tpoint / fl%t_PPL(i))**fl%a_PPL(i)
      else
        lgtp = dlog10(tpoint)
        jl = int((lgtp - fl%lgtcoolmin) /fl%lgstep) + 1
        Lpoint = fl%Lcool(jl)+ (tpoint-fl%tcool(jl)) &
                  * (fl%Lcool(jl+1)-fl%Lcool(jl)) &
                  / (fl%tcool(jl+1)-fl%tcool(jl))
      end if

    end subroutine findL

    subroutine findY (tpoint,Ypoint,fl)
    !  Fast search option to find correct point in cooling time 
      use mod_global_parameters

      double precision,intent(IN)   :: tpoint
      double precision, intent(OUT) :: Ypoint
      type(rc_fluid), intent(in) :: fl

      double precision :: lgtp
      double precision :: y_extra,factor
      integer :: jl,i

      if(fl%isPPL) then
        i = maxloc(fl%t_PPL, dim=1, mask=fl%t_PPL<tpoint)
        factor = fl%l_PPL(fl%n_PPL+1) * fl%t_PPL(i) / (fl%l_PPL(i) * fl%t_PPL(fl%n_PPL+1))
        if(fl%a_PPL(i)==1.d0) then
          y_extra = log( fl%t_PPL(i) / tpoint )
        else
          y_extra = 1 / (1 - fl%a_PPL(i)) * (1 - ( fl%t_PPL(i) / tpoint )**(fl%a_PPL(i)-1) )
        end if
        Ypoint = fl%y_PPL(i) + factor*y_extra
      else
        lgtp = dlog10(tpoint)
        jl = int((lgtp - fl%lgtcoolmin) / fl%lgstep) + 1
        ! Bounds check: jl must satisfy 1 <= jl <= ncool-1 so jl+1 <= ncool
        if(jl < 1 .or. jl >= fl%ncool) then
          write(*,'(a,es14.6,a,i0,a,2es14.6)') &
            'findY: tpoint=',tpoint,' jl=',jl,' out of bounds [1,ncool-1]; tcoolmin/max=', &
            fl%tcoolmin,fl%tcoolmax
          call mpistop('findY: temperature index out of bounds')
        end if
        Ypoint = fl%Yc(jl)+ (tpoint-fl%tcool(jl)) &
                  * (fl%Yc(jl+1)-fl%Yc(jl)) &
                  / (fl%tcool(jl+1)-fl%tcool(jl))
      end if

    end subroutine findY

    subroutine findT (tpoint,Ypoint,fl)
    !  Fast search option to find correct temperature 
    !  from temporal evolution function. Only possible this way because T is a monotonously
    !  decreasing function for the interpolated tables
    !  Uses eq. A7 from Townsend 2009 for piecewise power laws
      use mod_global_parameters

      double precision,intent(OUT)   :: tpoint
      double precision, intent(IN) :: Ypoint
      type(rc_fluid), intent(in) :: fl

      double precision :: factor
      integer :: jl,jc,jh,i
      
      if(fl%isPPL) then
        i = minloc(fl%y_PPL, dim=1, mask=fl%y_PPL>Ypoint)
        factor =  fl%l_PPL(i) * fl%t_PPL(fl%n_PPL+1) / (fl%l_PPL(fl%n_PPL+1) * fl%t_PPL(i))
        if(fl%a_PPL(i)==1.d0) then
          tpoint = fl%t_PPL(i) * exp( -1.d0 * factor * ( Ypoint - fl%y_PPL(i)))
        else
          tpoint = fl%t_PPL(i) * (1 - (1 - fl%a_PPL(i)) * factor * (Ypoint - fl%y_PPL(i)))**(1 / (1 - fl%a_PPL(i)))
        end if
      else
        if(Ypoint >= fl%Yc(1)) then
          tpoint = fl%tcoolmin
        else if (Ypoint == fl%Yc(fl%ncool)) then
          tpoint = fl%tcoolmax
        else
          jl=0
          jh=fl%ncool+1
          do
            if(jh-jl <= 1) exit
            jc=(jh+jl)/2
            if(Ypoint <= fl%Yc(jc)) then
              jl=jc
            else
              jh=jc
            end if
          end do
          ! Linear interpolation to obtain correct temperature
          tpoint = fl%tcool(jl)+ (Ypoint-fl%Yc(jl)) &
                 * (fl%tcool(jl+1)-fl%tcool(jl)) &
                 / (fl%Yc(jl+1)-fl%Yc(jl))
        end if
      end if
    end subroutine findT

    subroutine finddLdt (tpoint,dLpoint,fl)
    !  Fast search option to find correct point 
    !  in derivative of cooling curve
    !  Does not work for the piecewise power laws
      use mod_global_parameters

      double precision,intent(IN)   :: tpoint
      double precision, intent(OUT) :: dLpoint
      type(rc_fluid), intent(in) :: fl

      double precision :: lgtp
      integer :: jl,jc,jh

      lgtp = dlog10(tpoint)
      jl = int((lgtp -fl%lgtcoolmin) / fl%lgstep) + 1
      dLpoint = fl%dLdtcool(jl)+ (tpoint-fl%tcool(jl)) &
                * (fl%dLdtcool(jl+1)-fl%dLdtcool(jl)) &
                / (fl%tcool(jl+1)-fl%tcool(jl))

!      if (tpoint == tcoolmin) then
!        dLpoint = dLdtcool(1)
!      else if (tpoint == tcoolmax) then
!        dLpoint = dLdtcool(ncool)
!      else
!        jl=0
!        jh=ncool+1  
!        do
!          if (jh-jl <= 1) exit
!          jc=(jh+jl)/2
!          if (tpoint >= tcool(jc)) then
!              jl=jc
!          else
!              jh=jc
!          end if
!        end do
!        ! Linear interpolation to obtain correct cooling derivative
!        dLpoint = dLdtcool(jl)+ (tpoint-tcool(jl)) &
!                  * (dLdtcool(jl+1)-dLdtcool(jl)) &
!                  / (tcool(jl+1)-tcool(jl))
!      end if
    end subroutine finddLdt

    !> ===================================================================
    !> Variable-c_V Townsend extension (Y_mod)
    !> ===================================================================
    !>
    !> The original Townsend (2009) exact integration scheme assumes a
    !> constant heat capacity c_V = ρR/(γ-1). For LTE plasmas with H/He
    !> ionisation, c_V is a strong function of temperature in the
    !> recombination zone (Ibañez 1985, 1992 thermal-instability buffering).
    !>
    !> The derivation defines a modified TEF
    !>     Ỹ(T; ρ) ≡ ∫_T^{T_ref} c_V(T'; ρ) / (n_H n_e(T'; ρ) Λ(T')) dT'
    !> with c_V the LTE volumetric heat capacity *including* the
    !> ionisation-energy reservoir. The integral is recast via the change
    !> of variables u = e_int/n_H so that the discrete construction never
    !> finite-differences c_V — both integrand and energy update are
    !> driven by the same eint_from_T table, guaranteeing bit-consistent
    !> energy conservation.
    !>
    !> Must be called *after* eos_finalise() so that eos%eint_from_T,
    !> eos%T (forward), and eos%neOnH are all in code units. The hook is
    !> bind_eos_to_source() in mod_hd_eos.t / mod_mhd_eos.t.
    !>
    !> Algorithm: change of variables from T to u = e_int/n_H. For each
    !> log10 nH grid point j (taken from the eos%eint_from_T table's nH axis):
    !>   1. Cache u_i = e_int/n_H at every cooling-curve T_i.
    !>   2. Walk i = ncool-1 down to 1 and accumulate
    !>        Y_mod(j, i) = Y_mod(j, i+1) + ∫_{u_i}^{u_{i+1}} du / (n_e Λ)
    !>      using a composite Simpson (3-point, O(h^4)) or Boole (5-point,
    !>      O(h^6)) rule with N_sub sub-intervals.
    !>   3. The integrand evaluates n_e via y_from_nH_eint and Λ via findL
    !>      (with T from T_from_nH_eint).
    !>
    !> Per-row inverse table T_mod_inv(j, k) is also built for the
    !> 'table' inverse method (alternative to bisection).
    subroutine build_Y_mod_table(fl)
      use mod_global_parameters, only: mype
      type(rc_fluid), intent(inout) :: fl

      integer :: i, j, k, n_nH, ncool, N_sub
      double precision :: log_nH_j, nH_j_code, u_lo, u_hi, du_total, du_step
      double precision :: u_s, log_u_s, T_s, y_s, ne_s, Lambda_s, integ
      double precision, allocatable :: u_at_T(:), f_node(:)
      double precision :: Y_max_global, Y_min_global

      ! Preconditions for the variable-c_V Townsend extension:
      !   1. LTE with ionisation energy (otherwise c_V is constant and classical
      !      Townsend already correct)
      !   2. Townsend exact cooling method (Y-advance is its signature)
      !   3. Tabulated cooling curve, not piecewise power law (PPL has its own
      !      y_PPL path and does not allocate fl%tcool / fl%Lcool)
      !   4. eos%eint_from_T table must be built (excludes analytic Saha mode)
      !   5. Build only once per run
      if (.not. fl%ionE) return
      if (fl%isPPL) return
      if (fl%Y_mod_built) return

      ! The (log_nH, log_T) inverse-table grid (extents + n_nH) was snapshotted
      ! into fl by eos_get_eintT_grid in bind_eos_to_source; n_nH=0 means no such
      ! table (analytic/FI). The build below queries fl%eint_from_T which already
      ! dispatches on the EoS method, so we only need the grid extents here.
      n_nH = fl%Y_mod_n_nH
      if (n_nH <= 0) then
        if (mype == 0) write(*,*) ' build_Y_mod_table: no (rho,T) inverse table allocated; skipping'
        return
      end if

      ! Validate quadrature option
      select case (trim(fl%Y_mod_quadrature))
      case ('simpson', 'boole')
        ! ok
      case default
        call mpistop('build_Y_mod_table: rc_Y_mod_quadrature must be simpson or boole')
      end select

      ncool = fl%ncool
      N_sub = max(2, fl%Y_mod_N_sub)
      ! For Boole's rule we need N_sub to be a multiple of 4; round up if not.
      if (trim(fl%Y_mod_quadrature) == 'boole') then
        if (mod(N_sub, 4) /= 0) N_sub = N_sub + (4 - mod(N_sub, 4))
      else
        ! Simpson needs N_sub even
        if (mod(N_sub, 2) /= 0) N_sub = N_sub + 1
      end if

      fl%Y_mod_n_nH = n_nH
      if (n_nH > 1) then
        fl%Y_mod_lg_nH_step_inv = dble(n_nH - 1) &
            / (fl%Y_mod_lg_nH_max - fl%Y_mod_lg_nH_min)
      else
        fl%Y_mod_lg_nH_step_inv = 0.0d0
      end if

      allocate(fl%Y_mod(n_nH, ncool))
      allocate(fl%Y_mod_max_per_row(n_nH))
      allocate(u_at_T(ncool))
      allocate(f_node(0:N_sub))

      do j = 1, n_nH
        log_nH_j = fl%Y_mod_lg_nH_min &
            + dble(j - 1) * (fl%Y_mod_lg_nH_max - fl%Y_mod_lg_nH_min) / dble(max(1, n_nH - 1))
        nH_j_code = 10.0d0**log_nH_j

        ! Cache u_i = e_int/n_H at each cooling-curve temperature
        do i = 1, ncool
          u_at_T(i) = fl%eint_from_T(log_nH_j, dlog10(fl%tcool(i)))
        end do

        fl%Y_mod(j, ncool) = 0.0d0

        ! Step downward in T, accumulating ∫du / (n_e Λ)
        do i = ncool - 1, 1, -1
          u_lo = u_at_T(i)
          u_hi = u_at_T(i + 1)
          du_total = u_hi - u_lo
          if (du_total <= 0.0d0) then
            ! Degenerate segment (shouldn't happen physically); skip
            fl%Y_mod(j, i) = fl%Y_mod(j, i + 1)
            cycle
          end if
          du_step = du_total / dble(N_sub)

          ! Sample integrand at composite quadrature nodes
          do k = 0, N_sub
            u_s = u_lo + dble(k) * du_step
            if (u_s <= 0.0d0) then
              f_node(k) = 0.0d0
              cycle
            end if
            log_u_s = dlog10(u_s)
            T_s = fl%T_from_eint(log_nH_j, log_u_s)
            if (fl%lambda_needs_nenh_table) then
              ! SPEX-style two-table convention: the equilibrium n_e/n_H is
              ! already absorbed into Lambda_table at construction time, so
              ! the cooling rate is Q = n_H^2 * Lambda_table. The integrand
              ! 1/(n_e * Lambda) becomes 1/(n_H * Lambda); equivalently,
              ! substitute n_e -> n_H by setting y_s = 1.
              y_s = 1.0d0
              ne_s = nH_j_code
            else
              y_s = fl%y_from_eint(log_nH_j, log_u_s)
              ne_s = y_s * nH_j_code
            end if
            if (T_s <= fl%tcoolmin) then
              ! Below cooling table: no cooling, so integrand = 0
              f_node(k) = 0.0d0
              cycle
            else if (T_s >= fl%tcoolmax) then
              call calc_l_extended(T_s, Lambda_s, fl)
            else
              call findL(T_s, Lambda_s, fl)
            end if
            if (ne_s * Lambda_s > 0.0d0) then
              f_node(k) = 1.0d0 / (ne_s * Lambda_s)
            else
              f_node(k) = 0.0d0
            end if
          end do

          select case (trim(fl%Y_mod_quadrature))
          case ('boole')
            integ = boole_composite(f_node, N_sub, du_step)
          case default
            integ = simpson_composite(f_node, N_sub, du_step)
          end select

          fl%Y_mod(j, i) = fl%Y_mod(j, i + 1) + integ
        end do

        fl%Y_mod_max_per_row(j) = fl%Y_mod(j, 1)
      end do

      deallocate(u_at_T)
      deallocate(f_node)

      fl%Y_mod_built = .true.

      if (mype == 0) then
        Y_max_global = maxval(fl%Y_mod_max_per_row)
        Y_min_global = minval(fl%Y_mod_max_per_row)
        write(*,'(A,I0,A,I0,A,A,A,I0)') &
          ' Y_mod table built: ', n_nH, ' nH x ', ncool, ' T   quadrature=', &
          trim(fl%Y_mod_quadrature), '   N_sub=', N_sub
        write(*,'(A,F8.4,A,F8.4)') &
          '   log10 nH range = ', fl%Y_mod_lg_nH_min, ' to ', fl%Y_mod_lg_nH_max
        write(*,'(A,ES12.4,A,ES12.4,A)') &
          '   Y_max per row range = [', Y_min_global, ', ', Y_max_global, '] code time'
        write(*,'(A)') '   inverse=bisect (row-interpolated, O(log ncool))'
      end if
    end subroutine build_Y_mod_table

    !> Composite Simpson's rule on (N+1) equally spaced samples (N even).
    !> N must be a positive even integer; h is the step size.
    function simpson_composite(f, N, h) result(s)
      integer, intent(in) :: N
      double precision, intent(in) :: f(0:N), h
      double precision :: s
      integer :: k
      s = f(0) + f(N)
      do k = 1, N - 1, 2
        s = s + 4.0d0 * f(k)
      end do
      do k = 2, N - 2, 2
        s = s + 2.0d0 * f(k)
      end do
      s = s * h / 3.0d0
    end function simpson_composite

    !> Composite Boole's rule on (N+1) equally spaced samples (N a multiple of 4).
    !> Each 4-step block contributes (2h/45)*(7 f0 + 32 f1 + 12 f2 + 32 f3 + 7 f4).
    function boole_composite(f, N, h) result(s)
      integer, intent(in) :: N
      double precision, intent(in) :: f(0:N), h
      double precision :: s
      integer :: k
      s = 0.0d0
      do k = 0, N - 4, 4
        s = s + 7.0d0  * f(k)     &
              + 32.0d0 * f(k + 1) &
              + 12.0d0 * f(k + 2) &
              + 32.0d0 * f(k + 3) &
              + 7.0d0  * f(k + 4)
      end do
      s = s * 2.0d0 * h / 45.0d0
    end function boole_composite

    !> Bisection on a single Y_mod row to find log10(T) such that
    !> Y_mod_row(i) = y_target. The row is monotonically increasing in
    !> *decreasing* i (since Y(ncool)=0 grows toward Y(1)=Y_max).
    !> Returns log10 T (code units).
    function invert_row_bisect(Y_row, t_grid, ncool, y_target) result(log_T_out)
      integer, intent(in) :: ncool
      double precision, intent(in) :: Y_row(ncool), t_grid(ncool), y_target
      double precision :: log_T_out
      integer :: jl, jh, jc
      double precision :: f_lo, f_hi

      if (y_target <= Y_row(ncool)) then
        log_T_out = dlog10(t_grid(ncool))
        return
      end if
      if (y_target >= Y_row(1)) then
        log_T_out = dlog10(t_grid(1))
        return
      end if

      ! Bracket: find jl, jh = jl+1 such that Y_row(jh) <= y_target <= Y_row(jl)
      jl = 1
      jh = ncool
      do
        if (jh - jl <= 1) exit
        jc = (jl + jh) / 2
        if (Y_row(jc) >= y_target) then
          jl = jc
        else
          jh = jc
        end if
      end do
      f_lo = Y_row(jl)
      f_hi = Y_row(jh)
      if (f_lo == f_hi) then
        log_T_out = 0.5d0 * (dlog10(t_grid(jl)) + dlog10(t_grid(jh)))
      else
        ! Linear interpolation in log T (the cooling-curve T grid is uniform in log T)
        log_T_out = dlog10(t_grid(jl)) &
            + (y_target - f_lo) / (f_hi - f_lo) &
              * (dlog10(t_grid(jh)) - dlog10(t_grid(jl)))
      end if
    end function invert_row_bisect

    !> Forward Y_mod lookup: bilinear interpolation in (log10 nH, log10 T)
    !> on the precomputed Y_mod table. Both axes are uniform in log space.
    function findY_mod(Te_loc, nH_loc, fl) result(Y_out)
      double precision, intent(in) :: Te_loc, nH_loc
      type(rc_fluid), intent(in) :: fl
      double precision :: Y_out
      double precision :: log_nH, log_T, ry, rx
      integer :: jy, jy1, jx, jx1
      double precision :: fy, fx

      log_nH = dlog10(nH_loc)
      log_T  = dlog10(Te_loc)

      ! Clamp into table range
      ry = (log_nH - fl%Y_mod_lg_nH_min) * fl%Y_mod_lg_nH_step_inv
      ry = max(0.0d0, min(ry, dble(fl%Y_mod_n_nH - 1)))
      jy  = int(ry)
      jy1 = min(jy + 1, fl%Y_mod_n_nH - 1)
      fy  = ry - dble(jy)

      rx = (log_T - fl%lgtcoolmin) / fl%lgstep
      rx = max(0.0d0, min(rx, dble(fl%ncool - 1)))
      jx  = int(rx)
      jx1 = min(jx + 1, fl%ncool - 1)
      fx  = rx - dble(jx)

      Y_out = (1.0d0 - fy) * ((1.0d0 - fx) * fl%Y_mod(jy + 1,  jx + 1) &
                             +         fx  * fl%Y_mod(jy + 1,  jx1 + 1)) &
            +         fy  * ((1.0d0 - fx) * fl%Y_mod(jy1 + 1, jx + 1) &
                             +         fx  * fl%Y_mod(jy1 + 1, jx1 + 1))
    end function findY_mod

    !> Inverse Y_mod lookup: given Y_target and nH, return T such that
    !> Y_mod(log10 nH, log10 T) = Y_target. Bisection on the cooling-table
    !> i index using values interpolated linearly between the two adjacent
    !> nH rows (same convention as the classical findT). Saturates at the
    !> table extremes when Y_target falls outside [Y(tcoolmax), Y(tcoolmin)].
    function findT_mod(Y_target, nH_loc, fl) result(T_out)
      double precision, intent(in) :: Y_target, nH_loc
      type(rc_fluid), intent(in) :: fl
      double precision :: T_out
      double precision :: log_nH, ry, fy
      integer :: jy, jy1
      double precision :: log_T_lo, log_T_hi
      integer :: jl, jh, jc, ncool
      double precision :: Yc_lo, Yc_hi

      log_nH = dlog10(nH_loc)
      ry = (log_nH - fl%Y_mod_lg_nH_min) * fl%Y_mod_lg_nH_step_inv
      ry = max(0.0d0, min(ry, dble(fl%Y_mod_n_nH - 1)))
      jy  = int(ry)
      jy1 = min(jy + 1, fl%Y_mod_n_nH - 1)
      fy  = ry - dble(jy)
      ncool = fl%ncool

      Yc_lo = (1.0d0 - fy) * fl%Y_mod(jy + 1, 1)      + fy * fl%Y_mod(jy1 + 1, 1)
      Yc_hi = (1.0d0 - fy) * fl%Y_mod(jy + 1, ncool)  + fy * fl%Y_mod(jy1 + 1, ncool)
      if (Y_target >= Yc_lo) then
        T_out = fl%tcoolmin
        return
      end if
      if (Y_target <= Yc_hi) then
        T_out = fl%tcoolmax
        return
      end if
      jl = 1
      jh = ncool
      do
        if (jh - jl <= 1) exit
        jc = (jl + jh) / 2
        if (((1.0d0 - fy) * fl%Y_mod(jy + 1, jc) + fy * fl%Y_mod(jy1 + 1, jc)) &
            >= Y_target) then
          jl = jc
        else
          jh = jc
        end if
      end do
      Yc_lo = (1.0d0 - fy) * fl%Y_mod(jy + 1, jl) + fy * fl%Y_mod(jy1 + 1, jl)
      Yc_hi = (1.0d0 - fy) * fl%Y_mod(jy + 1, jh) + fy * fl%Y_mod(jy1 + 1, jh)
      log_T_lo = dlog10(fl%tcool(jl))
      log_T_hi = dlog10(fl%tcool(jh))
      if (Yc_lo == Yc_hi) then
        T_out = fl%tcool(jl)
      else
        T_out = 10.0d0**(log_T_lo &
              + (Y_target - Yc_lo) / (Yc_hi - Yc_lo) * (log_T_hi - log_T_lo))
      end if
    end function findT_mod

end module mod_radiative_cooling
