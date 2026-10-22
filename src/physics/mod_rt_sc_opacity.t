!=============================================================================
!> Opacity and Planck source of the short-characteristics (SC) radiative transfer
!>
!> The solvers never evaluate the gas themselves. They call the rt_sc_chiB callback, which
!> returns chi = kappa rho and B = a T^4 on a block, with the temperature and gas pressure taken
!> from the EoS, so that the opacity follows whichever EoS the run uses. rt_sc_setup points the
!> callback at one of three models, chosen by rt_sc_opacity:
!>   Rosseland  the tabulated Rosseland mean, shifted to code units when it is read, so that a
!>              cell needs one logarithm each of T and rho and one power of ten;
!>   Kramers    an analytic H-minus and Kramers fit with electron scattering added;
!>   usr        the user routine usr_rt_sc_chiB.
!> The built-in models also apply the transition switch rt_sc_switch. It scales chi and B
!> together, so that the hot, tenuous corona becomes both transparent and dark to the SC
!> transfer, and the thin cooling takes over there. The module depends only on
!> mod_global_parameters and mod_eos
!>
!> Jack Jenkins (10.2026); More or less finished in v1 state for amrvac4.0
!=============================================================================
module mod_rt_sc_opacity
  use mod_global_parameters, only: std_len
  implicit none

  !> The opacity model, either 'Rosseland' (the table), 'Kramers' (the analytic fit) or 'usr',
  !> which uses the usr_rt_sc_chiB callback
  character(len=10), public :: rt_sc_opacity = 'Rosseland'
  !> The Rosseland table file, read in rt_sc_setup when the table model is chosen
  character(len=std_len), public :: rt_sc_kappa_table = 'rosseland_combined.tab'
  !> The built-in opacities keep the SC transfer on where the gas is cool (T < Tcap) or dense
  !> (p > p_TR), and switch it off smoothly in the hot, tenuous corona. The transition bands Tcap
  !> and dTcap are in K, p_TR is in dyn/cm^2 and its width dp_TR in dex
  double precision, public :: rt_sc_Tcap = 2.0d4, rt_sc_dTcap = 5.0d3
  double precision, public :: rt_sc_p_TR = 1.0d0, rt_sc_dp_TR = 0.5d0
  !> log10 of p_TR in code units, set by rt_sc_setup
  double precision, public :: rt_sc_lp_TR = 0.0d0
  !> The Rosseland table, shifted to code units when it is read. It holds log10 kappa on a grid
  !> of log10 T and log10 R, where R = rho/T^3, all three in code units
  double precision, allocatable :: op_logT(:), op_logR(:), op_logK(:,:)
  integer :: nT_op = 0, nR_op = 0

  !> The opacity and source the solvers call. rt_sc_setup points it at one of the built-in
  !> routines, or at usr_rt_sc_chiB for 'usr'
  abstract interface
    subroutine rt_sc_chiB_iface(ixI^L, ixO^L, w, x, chi, bpl)
      use mod_global_parameters
      integer, intent(in)          :: ixI^L, ixO^L
      double precision, intent(in) :: w(ixI^S,1:nw), x(ixI^S,1:ndim)
      double precision, intent(out):: chi(ixI^S), bpl(ixI^S)
    end subroutine rt_sc_chiB_iface
  end interface
  procedure(rt_sc_chiB_iface), pointer :: rt_sc_chiB => null()

  public :: rt_sc_chiB
  public :: rt_sc_chiB_rosseland
  public :: rt_sc_chiB_kramers
  public :: read_opacity_table
  public :: rt_sc_switch

contains

  !> Table Rosseland opacity, chi = kappa_R(rho, T) rho, and source B = a T^4, both multiplied by
  !> the transition switch
  subroutine rt_sc_chiB_rosseland(ixI^L, ixO^L, w, x, chi, bpl)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(out):: chi(ixI^S), bpl(ixI^S)
    call chiB_builtin(ixI^L, ixO^L, w, x, chi, bpl, .true.)
  end subroutine rt_sc_chiB_rosseland

  !> Analytic H-minus and Kramers opacity, with the source and the switch as in
  !> rt_sc_chiB_rosseland
  subroutine rt_sc_chiB_kramers(ixI^L, ixO^L, w, x, chi, bpl)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(out):: chi(ixI^S), bpl(ixI^S)
    call chiB_builtin(ixI^L, ixO^L, w, x, chi, bpl, .false.)
  end subroutine rt_sc_chiB_kramers

  !> chi and B of the built-in opacities, the table when tabulated is set and the analytic fit
  !> otherwise. Temperature and gas pressure come from the EoS, so the opacity follows whichever
  !> EoS the run uses, and the switch scales chi and B together, so that the switched-off gas is
  !> both transparent and dark
  subroutine chiB_builtin(ixI^L, ixO^L, w, x, chi, bpl, tabulated)
    use mod_global_parameters
    use mod_eos, only: eos
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(out):: chi(ixI^S), bpl(ixI^S)
    logical, intent(in)          :: tabulated
    double precision :: Tloc(ixI^S), pgas(ixI^S), Tc, rc, kcgs, scap
    integer :: ix^D
    call eos%get_Te(w, x, ixI^L, ixO^L, Tloc)
    call eos%get_thermal_pressure(w, x, ixI^L, ixO^L, pgas)
    {do ix^DB = ixOmin^DB, ixOmax^DB\}
      Tc = Tloc(ix^D)*unit_temperature
      if(tabulated) then
        ! The table is in code units, so the opacity is ready to be multiplied by the code density
        chi(ix^D) = kappa_code_cell(Tloc(ix^D), w(ix^D,iw_rho))*w(ix^D,iw_rho)
      else
        ! The fit is written in cgs. The H-minus term dominates in cool gas and the Kramers term in
        ! hot gas, and electron scattering at 0.34 cm^2/g sets a floor under both. The density is
        ! kept positive for the fractional powers
        ! kappa = 1e-25 rho^(1/2) T^7.7 + 4e25 rho T^(-3.5) + 0.34, in cgs
        rc   = w(ix^D,iw_rho)*unit_density
        kcgs = 1.0d-25*dsqrt(max(rc,1.0d-30))*Tc**7.7d0 &
             + 4.0d25*max(rc,1.0d-30)*Tc**(-3.5d0)  + 0.34d0
        chi(ix^D) = (kcgs/unit_opacity)*w(ix^D,iw_rho)
      end if
      ! chi = s kappa rho,  B = s a T^4
      bpl(ix^D) = arad_norm*Tloc(ix^D)**4
      scap = rt_sc_switch(Tc, pgas(ix^D))
      chi(ix^D) = chi(ix^D)*scap
      bpl(ix^D) = bpl(ix^D)*scap
    {end do\}
  end subroutine chiB_builtin

  !> Transition switch of the built-in opacities, from the temperature in K and the gas pressure in
  !> code units. The temperature gate is one in cool gas and falls to zero across dTcap above
  !> Tcap, while the pressure gate rises from zero to one across dp_TR dex around p_TR. A cell
  !> keeps its transfer when either holds, so the larger gate wins
  pure double precision function rt_sc_switch(T_K, p_code) result(scap)
    double precision, intent(in) :: T_K, p_code
    double precision :: tsw, psw
    ! s_T = (1 - tanh((T - Tcap)/dTcap))/2
    ! s_p = (1 + tanh((log10 p - log10 p_TR)/dp_TR))/2
    tsw  = 0.5d0*(1.0d0 - tanh((T_K - rt_sc_Tcap)/rt_sc_dTcap))
    psw  = 0.5d0*(1.0d0 + tanh((dlog10(max(p_code,1.0d-99)) - rt_sc_lp_TR)/rt_sc_dp_TR))
    ! s = max(s_T, s_p)
    scap = max(tsw, psw)
  end function rt_sc_switch

  !> Rosseland opacity of one cell from the table, with temperature, density and opacity all in
  !> code units
  double precision function kappa_code_cell(T_code, rho_code) result(kap)
    double precision, intent(in) :: T_code, rho_code
    double precision :: lT, lR
    ! The table axes are log10 T and log10 R, and log10 R = log10 rho - 3 log10 T
    lT  = dlog10(T_code)
    lR  = dlog10(rho_code) - 3.0d0*lT
    kap = 10.0d0**interp_logkappa(lT, lR)
  end function kappa_code_cell

  !> Bilinear interpolation of log10 kappa in (log10 T, log10 R)
  double precision function interp_logkappa(lT, lR) result(logK)
    double precision, intent(in) :: lT, lR
    integer :: iT, iR
    double precision :: tT, tR, ka, kb
    ! op_locate gives the lower corner of the table cell that holds the point
    iT = op_locate(op_logT, nT_op, lT)
    iR = op_locate(op_logR, nR_op, lR)
    ! tT = (lT - lT_i)/(lT_(i+1) - lT_i) and tR likewise, then linear in log R and in log T
    tT = (lT - op_logT(iT)) / (op_logT(iT+1) - op_logT(iT))
    tR = (lR - op_logR(iR)) / (op_logR(iR+1) - op_logR(iR))
    tT = min(max(tT, 0.d0), 1.d0);  tR = min(max(tR, 0.d0), 1.d0)
    ! We interpolate along R on the two bracketing rows of T, and then along T between them
    ka = op_logK(iT,  iR) + tR*(op_logK(iT,  iR+1) - op_logK(iT,  iR))
    kb = op_logK(iT+1,iR) + tR*(op_logK(iT+1,iR+1) - op_logK(iT+1,iR))
    logK = ka + tT*(kb - ka)
  end function interp_logkappa

  !> Index jl of the table interval [xx(jl), xx(jl+1)] that holds xq, found by bisection. A point
  !> beyond either end returns the first or the last interval
  integer function op_locate(xx, n, xq) result(jl)
    integer, intent(in) :: n
    double precision, intent(in) :: xx(n), xq
    integer :: ju, jm
    if (xq <= xx(1)) then; jl = 1; return; end if
    if (xq >= xx(n)) then; jl = n-1; return; end if
    jl = 1; ju = n
    do while (ju - jl > 1)
      jm = (ju + jl)/2
      if (xq >= xx(jm)) then; jl = jm; else; ju = jm; end if
    end do
  end function op_locate

  !> Read the Rosseland table and shift it to code units. The file is in cgs, with T in K and
  !> R = rho/T6^3, where T6 is the temperature in 1e6 K
  subroutine read_opacity_table()
    use mod_global_parameters
    use mod_comm_lib, only: mpistop
    integer :: i, funit
    logical :: alive
    character(len=256) :: dum
    
    inquire(file=trim(rt_sc_kappa_table), exist=alive)
    if (.not. alive) call mpistop('read_opacity_table: Rosseland table not found: '//trim(rt_sc_kappa_table))
    open(newunit=funit, file=trim(rt_sc_kappa_table), status='old')
    read(funit,'(A)') dum
    read(funit,*) nT_op, nR_op
    if(allocated(op_logT)) deallocate(op_logT, op_logR, op_logK)
    allocate(op_logT(nT_op), op_logR(nR_op), op_logK(nT_op,nR_op))
    read(funit,*) op_logR(1:nR_op)
    do i = 1, nT_op
      read(funit,*) op_logT(i), op_logK(i,1:nR_op)
    end do
    close(funit)
    ! Shift the units of the table so it's compliant with the AMRVAC unit system
    op_logT = op_logT - dlog10(unit_temperature)
    op_logR = op_logR - dlog10(unit_density) + 3.0d0*dlog10(1.0d-6*unit_temperature)
    op_logK = op_logK - dlog10(unit_opacity)
    if (mype==0) print *,'read_opacity_table: Rosseland table ',nT_op,' logT x ',nR_op,' logR'
  end subroutine read_opacity_table

end module mod_rt_sc_opacity
