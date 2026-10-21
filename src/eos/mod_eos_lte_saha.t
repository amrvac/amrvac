!=============================================================================
!> Analytic H-only Saha EoS (eos_method == 'analytic').
!>
!> Ground-state Saha for pure hydrogen, solved per call (no tables):
!>   X(T)        = saha_pf * T^{3/2} * exp(-chi/kB/T)
!>   y(=ne/nH)   = 2X / (X + sqrt(X^2 + 4 nH X))
!> Public entry points take and return CODE units; CGS conversion is internal.
!> T(nH, eint/nH) is found by bisection (default) or Newton; p(nH, T) closes
!> the state. build_gamma1_analytic_table tabulates Gamma_1(nH, T) once for the
!> 'exact' gamma1 path; saha_gamma1_from_nH_T reads it back.
!=============================================================================
module mod_eos_lte_saha
    use mod_global_parameters
    use mod_eos_container, only: eos, EOS_ANALYTIC
    use mod_eos_interp,    only: bicubic_lookup, precompute_step_inv
    use mod_comm_lib,      only: mpistop
    use mod_eos_lte_tables, only: precompute_FI_bypass_constants
    implicit none
    private

    public :: saha_eint_from_nH_T, saha_state_from_nH_p
    public :: saha_T_from_nH_eint, saha_gamma1_from_nH_T
    public :: load_analytic_LTE, finalise_analytic_LTE

    ! Analytical H-only Saha constants (module-level parameters)
    ! Prefactor: (2*pi*m_e*k_B/h^2)^{3/2} in CGS (cm^-3)
    ! SI value is 2.4146830395719654e21 m^-3; divide by 1e6 for CGS
    double precision, parameter :: saha_pf = 2.4146830395719654d15
    double precision, parameter :: saha_chi_kB = 157763.42386247337d0
    double precision, parameter :: saha_kB_cgs = 1.380649d-16
    double precision, parameter :: saha_chi_H_cgs = 2.178710282685096d-11

    !> One-shot guard: a bisection that exhausts max_iter without meeting the
    !> tolerance warns once (rank 0) rather than silently returning an
    !> under-converged T every offending cell.
    logical :: saha_nonconv_warned = .false.

contains

    !> Analytic-method load: no binary tables (on-the-fly Saha solve); just
    !> announce the method.
    subroutine load_analytic_LTE()
        if (mype == 0) then
            write(*,*) 'EoS method: analytic (H-only Saha)'
            write(*,*) '  gamma1_method: ', trim(eos%gamma1_method)
        end if
    end subroutine load_analytic_LTE

    !> Analytic-method finalise (H-only Saha): no table unit conversion needed;
    !> build the Gamma1 table only for 'exact' mode, then the shared FI-bypass
    !> constants.
    subroutine finalise_analytic_LTE()
        if (eos%gamma1_method == 'exact') then
            call build_gamma1_analytic_table()
        end if
        call precompute_FI_bypass_constants()
    end subroutine finalise_analytic_LTE

    !> Ionization fraction from Saha in CGS (no unit conversions).
    !> Use this inside bisection loops where nH_cgs is constant.
    double precision function saha_y_cgs(nH_cgs, T_cgs) result(y)
        double precision, intent(in) :: nH_cgs, T_cgs
        double precision :: X_saha, disc

        X_saha = saha_pf * T_cgs * dsqrt(T_cgs) * dexp(-saha_chi_kB / T_cgs)
        disc = dsqrt(X_saha * X_saha + 4.0d0 * nH_cgs * X_saha)
        ! Use numerically stable form to avoid catastrophic cancellation:
        ! When X >> nH (high ionization): standard form loses precision in -X + sqrt(...)
        ! Stable form: y = 2*X / (X + sqrt(X^2 + 4*nH*X)) -> 1 as X -> inf
        y = 2.0d0 * X_saha / (X_saha + disc)

    end function saha_y_cgs

    !> Ionization fraction y = ne/nH from analytical Saha at given (nH, T) in CODE UNITS
    double precision function saha_y_from_nH_T(nH_code, T_code) result(y)
        double precision, intent(in) :: nH_code, T_code

        y = saha_y_cgs(nH_code * unit_numberdensity, T_code * unit_temperature)

    end function saha_y_from_nH_T

    !> Internal energy per nH from analytical Saha, in CODE UNITS.
    !> eint_nH = 1.5*(1+y)*kB*T [+ y*chi_H] in CGS, converted to code units
    double precision function saha_eint_from_nH_T(nH_code, T_code) result(eint_nH_code)
        double precision, intent(in) :: nH_code, T_code
        double precision :: y, T_cgs, eint_nH_cgs

        if (eos%method_id /= EOS_ANALYTIC) call mpistop("saha_eint_from_nH_T requires eos_method=analytic")
        y = saha_y_from_nH_T(nH_code, T_code)
        T_cgs = T_code * unit_temperature

        eint_nH_cgs = 1.5d0 * (1.0d0 + y) * saha_kB_cgs * T_cgs
        if (eos%ionE) then
            eint_nH_cgs = eint_nH_cgs + y * saha_chi_H_cgs
        end if

        eint_nH_code = eint_nH_cgs * unit_numberdensity / unit_pressure

    end function saha_eint_from_nH_T

    !> Given (nH, p) in CODE UNITS, find T and y by solving p = nH*(1+y(T))*T_code.
    !> Uses bisection on T. In code units, kB is absorbed: p = nH*(1+y)*T.
    !> Given (nH, p) in CODE UNITS, find T, y, and eint/nH.
    !> Bisects on p = nH*(1+y(T))*T_code, then computes eint from final T,y.
    !> Returns eint_nH in code units (avoids redundant Saha re-evaluation).
    subroutine saha_state_from_nH_p(nH_code, p_code, T_out, y_out, eint_nH_out)
        double precision, intent(in)  :: nH_code, p_code
        double precision, intent(out) :: T_out, y_out
        double precision, intent(out), optional :: eint_nH_out

        double precision :: nH_cgs, uT, eint_nH_cgs
        double precision :: T_lo, T_hi, T_mid, y_mid, p_eval
        integer :: iter
        integer, parameter :: max_iter = 30
        double precision, parameter :: tol = 1.0d-8

        ! Hoist CGS conversion out of loop
        if (eos%method_id /= EOS_ANALYTIC) call mpistop("saha_state_from_nH_p requires eos_method=analytic")
        nH_cgs = nH_code * unit_numberdensity
        uT = unit_temperature

        ! Brackets in code units: y=0 => T = p/nH, y=1 => T = p/(2*nH)
        T_hi = p_code / nH_code
        T_lo = p_code / (2.0d0 * nH_code)
        T_lo = max(T_lo, 100.0d0 / uT)

        do iter = 1, max_iter
            T_mid = 0.5d0 * (T_lo + T_hi)
            y_mid = saha_y_cgs(nH_cgs, T_mid * uT)
            p_eval = nH_code * (1.0d0 + y_mid) * T_mid

            if (p_eval < p_code) then
                T_lo = T_mid
            else
                T_hi = T_mid
            end if

            if (dabs(T_hi - T_lo) < tol * T_mid) exit
        end do

        if (iter > max_iter .and. .not. saha_nonconv_warned .and. mype == 0) then
            saha_nonconv_warned = .true.
            write(*,*) "WARNING: saha_state_from_nH_p bisection reached max_iter ", &
                 "without converging; using last iterate. Not reported again."
        end if

        T_out = T_mid
        y_out = y_mid

        ! Optionally return eint/nH directly (avoids redundant Saha re-evaluation)
        if (present(eint_nH_out)) then
            eint_nH_cgs = 1.5d0 * (1.0d0 + y_mid) * saha_kB_cgs * (T_mid * uT)
            if (eos%ionE) eint_nH_cgs = eint_nH_cgs + y_mid * saha_chi_H_cgs
            eint_nH_out = eint_nH_cgs * unit_numberdensity / unit_pressure
        end if

    end subroutine saha_state_from_nH_p

    !> Temperature inversion: given (nH, eint/nH) in CODE UNITS, find T in CODE UNITS,
    !> by bisection (guaranteed convergence).
    subroutine saha_T_from_nH_eint(nH_code, eint_nH_code, T_out, y_out)
        double precision, intent(in)  :: nH_code, eint_nH_code
        double precision, intent(out) :: T_out, y_out

        if (eos%method_id /= EOS_ANALYTIC) call mpistop("saha_T_from_nH_eint requires eos_method=analytic")
        call saha_T_bisection(nH_code, eint_nH_code, T_out, y_out)

    end subroutine saha_T_from_nH_eint

    !> Bisection solver for T(nH, eint/nH). Guaranteed convergence.
    !> Brackets: T_lo from fully ionized (y=1), T_hi from neutral (y=0).
    !> Handles the non-monotonic eint(T) when ionE is included by
    !> checking the FI limit first: if y~1 at T_lo, use T_lo directly.
    subroutine saha_T_bisection(nH_code, eint_nH_code, T_out, y_out)
        double precision, intent(in)  :: nH_code, eint_nH_code
        double precision, intent(out) :: T_out, y_out

        double precision :: nH_cgs, eint_nH_cgs, T_lo, T_hi, T_mid
        double precision :: eint_mid, y_mid, y_lo
        integer :: iter
        integer, parameter :: max_iter = 30
        double precision, parameter :: tol = 1.0d-8

        ! Hoist CGS conversions out of loop
        nH_cgs = nH_code * unit_numberdensity
        eint_nH_cgs = eint_nH_code * unit_pressure / unit_numberdensity

        ! Temperature bounds in CGS from energy equation limits
        if (eos%ionE) then
            ! T_lo: assume y=1 (fully ionized) -> eint = 3*kB*T + chi_H
            T_lo = max((eint_nH_cgs - saha_chi_H_cgs) / (3.0d0 * saha_kB_cgs), 100.0d0)
            ! T_hi: assume y=0 (neutral) -> eint = 1.5*kB*T
            T_hi = eint_nH_cgs / (1.5d0 * saha_kB_cgs)
        else
            T_lo = eint_nH_cgs / (3.0d0 * saha_kB_cgs)
            T_hi = eint_nH_cgs / (1.5d0 * saha_kB_cgs)
        end if
        T_lo = max(T_lo, 100.0d0)
        T_hi = max(T_hi, T_lo + 1.0d0)

        ! Guard against non-monotonic eint(T) when ionE is active:
        ! At low nH, hydrogen is fully ionized (y=1) across a wide T range,
        ! making the Saha equation degenerate. The eint function has two roots:
        ! one at low T (ionized, eint ~ chi_H + small thermal) and one at
        ! high T (neutral, eint ~ thermal only). The bisection must stay on
        ! the correct branch.
        ! Fix: check y at T_lo. If y~1, use the FI analytical formula directly.
        if (eos%ionE) then
            y_lo = saha_y_cgs(nH_cgs, T_lo)
            if (y_lo > 0.999d0) then
                ! Fully ionized: eint = 3*kB*T + chi_H => T = (eint - chi_H)/(3*kB)
                T_out = max((eint_nH_cgs - saha_chi_H_cgs) &
                    / (3.0d0 * saha_kB_cgs), 100.0d0) / unit_temperature
                y_out = 1.0d0
                return
            end if
        end if

        ! Bisect in CGS (no unit conversions inside loop)
        do iter = 1, max_iter
            T_mid = 0.5d0 * (T_lo + T_hi)
            y_mid = saha_y_cgs(nH_cgs, T_mid)

            eint_mid = 1.5d0 * (1.0d0 + y_mid) * saha_kB_cgs * T_mid
            if (eos%ionE) eint_mid = eint_mid + y_mid * saha_chi_H_cgs

            if (eint_mid < eint_nH_cgs) then
                T_lo = T_mid
            else
                T_hi = T_mid
            end if

            if (dabs(T_hi - T_lo) < tol * T_mid) exit
        end do

        if (iter > max_iter .and. .not. saha_nonconv_warned .and. mype == 0) then
            saha_nonconv_warned = .true.
            write(*,*) "WARNING: saha_T_bisection reached max_iter without ", &
                 "converging; using last iterate. Not reported again."
        end if

        T_out = T_mid / unit_temperature
        y_out = y_mid

    end subroutine saha_T_bisection

    !> Build 2D Gamma1(nH, T) table from analytical Saha for the 'analytic' EoS method.
    !> Grid: uniform in (log10 nH, log10 T) with 256x256 points.
    !> Gamma1 = 1 + (1+y)*kB / Cv, where Cv = d(eint)/d(T) at constant rho.
    !> Stored in eos%gamma1_p for reuse by the existing csound2 routines
    !> (re-indexed: axis 2 = log10(T) instead of log10(p/nH)).
    subroutine build_gamma1_analytic_table()
        integer, parameter :: ng = 256
        integer :: i, j
        double precision :: log_nH_min, log_nH_max, log_T_min, log_T_max
        double precision :: dx1, dx2, log_nH, log_T
        double precision :: nH_code, T_code, y, yp, ym
        double precision :: Cv, g1_val, g1_min, g1_max
        double precision :: T_cgs, dT_cgs
        double precision :: eint_p, eint_m

        ! Range: log10(nH [CGS]) from 5 to 19, log10(T [K]) from 2.8 to 6.3
        ! Convert to code-unit ranges
        log_nH_min = 5.0d0  - dlog10(unit_numberdensity)
        log_nH_max = 19.0d0 - dlog10(unit_numberdensity)
        log_T_min  = 2.8d0  - dlog10(unit_temperature)
        log_T_max  = 6.3d0  - dlog10(unit_temperature)

        eos%gamma1_p%dim1 = ng
        eos%gamma1_p%dim2 = ng
        eos%gamma1_p%var1_min = log_nH_min
        eos%gamma1_p%var1_max = log_nH_max
        eos%gamma1_p%var2_min = log_T_min
        eos%gamma1_p%var2_max = log_T_max
        eos%gamma1_p%filename = 'computed_gamma1_analytic'

        allocate(eos%gamma1_p%table(ng, ng))

        dx1 = (log_nH_max - log_nH_min) / dble(ng - 1)
        dx2 = (log_T_max  - log_T_min)  / dble(ng - 1)

        g1_min = 1.0d30
        g1_max = -1.0d30

        do j = 1, ng
            log_T = log_T_min + (j - 1) * dx2
            T_code = 10.0d0**log_T
            T_cgs = T_code * unit_temperature
            ! Central difference step: 0.1% of T
            dT_cgs = 1.0d-3 * T_cgs

            do i = 1, ng
                log_nH = log_nH_min + (i - 1) * dx1
                nH_code = 10.0d0**log_nH

                y = saha_y_from_nH_T(nH_code, T_code)

                ! Compute Cv = d(eint/nH)/d(T) via central difference
                yp = saha_y_from_nH_T(nH_code, (T_cgs + dT_cgs) / unit_temperature)
                ym = saha_y_from_nH_T(nH_code, (T_cgs - dT_cgs) / unit_temperature)

                eint_p = 1.5d0 * (1.0d0 + yp) * saha_kB_cgs * (T_cgs + dT_cgs)
                eint_m = 1.5d0 * (1.0d0 + ym) * saha_kB_cgs * (T_cgs - dT_cgs)
                if (eos%ionE) then
                    eint_p = eint_p + yp * saha_chi_H_cgs
                    eint_m = eint_m + ym * saha_chi_H_cgs
                end if

                Cv = (eint_p - eint_m) / (2.0d0 * dT_cgs)

                ! Gamma1 = 1 + (1+y)*kB / Cv
                if (Cv > 0.0d0) then
                    g1_val = 1.0d0 + (1.0d0 + y) * saha_kB_cgs / Cv
                else
                    g1_val = eos%gamma
                end if

                eos%gamma1_p%table(i, j) = g1_val
                g1_min = min(g1_min, g1_val)
                g1_max = max(g1_max, g1_val)
            end do
        end do

        call precompute_step_inv(eos%gamma1_p)

        if (mype == 0) then
            write(*, '(A,F8.4,A,F8.4)') &
                ' Gamma1 analytic table built (nH x T): min = ', g1_min, ', max = ', g1_max
        end if

    end subroutine build_gamma1_analytic_table

    !> Look up Gamma1 from the analytical 2D table (nH, T axes in code units).
    !> For use when eos%method == 'analytic' and gamma1_method == 'exact'.
    double precision function saha_gamma1_from_nH_T(nH_code, T_code) result(g1)
        double precision, intent(in) :: nH_code, T_code
        double precision :: log_nH, log_T

        if (eos%method_id /= EOS_ANALYTIC) call mpistop("saha_gamma1_from_nH_T requires eos_method=analytic")
        log_nH = dlog10(nH_code)
        log_T  = dlog10(T_code)

        g1 = bicubic_lookup(log_nH, log_T, eos%gamma1_p)

    end function saha_gamma1_from_nH_T

end module mod_eos_lte_saha
!> Needs a line after to pass the preprocessor
