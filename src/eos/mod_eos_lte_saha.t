!=============================================================================
!> Analytic hydrogen Saha EoS with neutral helium (eos_method == 'analytic')
!>
!> Ground-state Saha for hydrogen, solved per call without tables:
!>   X(T)        = saha_pf * T^{3/2} * exp(-chi/kB/T)
!>   y(=ne/nH)   = 2X / (X + sqrt(X^2 + 4 nH X))
!> Public entry points take and return code units; the CGS conversion is internal.
!> T(nH, eint/nH) is found by bisection and p(nH, T) closes the state.
!> build_gamma1_analytic_table tabulates Gamma_1(nH, p/nH) once for the 'exact' gamma1 path,
!> on the axes that gamma1_from_nH_p reads for every method
!=============================================================================
module mod_eos_lte_saha
    use mod_global_parameters
    use mod_eos_container, only: eos, EOS_ANALYTIC
    use mod_eos_interp,    only: precompute_step_inv
    use mod_comm_lib,      only: mpistop
    use mod_eos_lte_tables, only: precompute_FI_bypass_constants
    implicit none
    private

    public :: saha_eint_from_nH_T, saha_state_from_nH_p
    public :: saha_T_from_nH_eint
    public :: load_analytic_LTE, finalise_analytic_LTE

    ! Prefactor (2 pi m_e k_B/h^2)^{3/2} in CGS (cm^-3)
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
            write(*,*) 'EoS method: analytic (hydrogen Saha, helium neutral)'
            write(*,*) '  gamma1_method: ', trim(eos%gamma1_method)
        end if
    end subroutine load_analytic_LTE

    !> Analytic-method finalise: no table unit conversion needed;
    !> build the Gamma1 table only for 'exact' mode, then the shared FI-bypass
    !> constants.
    subroutine finalise_analytic_LTE()
        if (eos%gamma1_method == 'exact') then
            call build_gamma1_analytic_table()
        end if
        call precompute_FI_bypass_constants(saha_chi_H_cgs)
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

    !> Ionisation fraction y = ne/nH from Saha at (nH, T) in code units
    double precision function saha_y_from_nH_T(nH_code, T_code) result(y)
        double precision, intent(in) :: nH_code, T_code

        y = saha_y_cgs(nH_code * unit_numberdensity, T_code * unit_temperature)

    end function saha_y_from_nH_T

    !> eint/nH in CGS with helium neutral, 1.5 (1+He+y) kB T, plus y chi_H when ionE
    double precision function saha_eint_cgs(y, T_cgs) result(e)
        double precision, intent(in) :: y, T_cgs
        e = 1.5d0 * (1.0d0 + eos%He_abundance + y) * saha_kB_cgs * T_cgs
        if (eos%ionE) e = e + y * saha_chi_H_cgs
    end function saha_eint_cgs

    !> eint/nH from Saha at (nH, T), in code units
    double precision function saha_eint_from_nH_T(nH_code, T_code) result(eint_nH_code)
        double precision, intent(in) :: nH_code, T_code

        if (eos%method_id /= EOS_ANALYTIC) call mpistop("saha_eint_from_nH_T requires eos_method=analytic")
        eint_nH_code = saha_eint_cgs(saha_y_from_nH_T(nH_code, T_code), T_code * unit_temperature) &
            * unit_numberdensity / unit_pressure

    end function saha_eint_from_nH_T

    !> Given (nH, p) in code units, T, y and optionally eint/nH: bisection on
    !> p = nH*(1+He+y(T))*T, with kB absorbed in the code units
    subroutine saha_state_from_nH_p(nH_code, p_code, T_out, y_out, eint_nH_out)
        double precision, intent(in)  :: nH_code, p_code
        double precision, intent(out) :: T_out, y_out
        double precision, intent(out), optional :: eint_nH_out

        double precision :: nH_cgs, uT
        double precision :: T_lo, T_hi, T_mid, y_mid, p_eval
        integer :: iter
        integer, parameter :: max_iter = 30
        double precision, parameter :: tol = 1.0d-8

        if (eos%method_id /= EOS_ANALYTIC) call mpistop("saha_state_from_nH_p requires eos_method=analytic")
        nH_cgs = nH_code * unit_numberdensity
        uT = unit_temperature

        ! brackets: y = 0 gives T = p/(nH (1+He)), y = 1 gives T = p/(nH (2+He))
        T_hi = p_code / (nH_code * (1.0d0 + eos%He_abundance))
        T_lo = p_code / (nH_code * (2.0d0 + eos%He_abundance))
        T_lo = max(T_lo, 100.0d0 / uT)

        do iter = 1, max_iter
            T_mid = 0.5d0 * (T_lo + T_hi)
            y_mid = saha_y_cgs(nH_cgs, T_mid * uT)
            p_eval = nH_code * (1.0d0 + eos%He_abundance + y_mid) * T_mid

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

        if (present(eint_nH_out)) &
            eint_nH_out = saha_eint_cgs(y_mid, T_mid * uT) * unit_numberdensity / unit_pressure

    end subroutine saha_state_from_nH_p

    !> Temperature and ionisation fraction from (nH, eint/nH) in code units, by bisection
    subroutine saha_T_from_nH_eint(nH_code, eint_nH_code, T_out, y_out)
        double precision, intent(in)  :: nH_code, eint_nH_code
        double precision, intent(out) :: T_out, y_out

        if (eos%method_id /= EOS_ANALYTIC) call mpistop("saha_T_from_nH_eint requires eos_method=analytic")
        call saha_T_bisection(nH_code, eint_nH_code, T_out, y_out)

    end subroutine saha_T_from_nH_eint

    !> Bisection for T(nH, eint/nH) between the fully ionised (y = 1) and neutral (y = 0) limits
    subroutine saha_T_bisection(nH_code, eint_nH_code, T_out, y_out)
        double precision, intent(in)  :: nH_code, eint_nH_code
        double precision, intent(out) :: T_out, y_out

        double precision :: nH_cgs, eint_nH_cgs, T_lo, T_hi, T_mid
        double precision :: y_mid, cv_ion, cv_neu
        integer :: iter
        integer, parameter :: max_iter = 30
        double precision, parameter :: tol = 1.0d-8

        nH_cgs = nH_code * unit_numberdensity
        eint_nH_cgs = eint_nH_code * unit_pressure / unit_numberdensity

        ! temperature bounds from the fully ionised (y = 1) and neutral (y = 0) limits
        cv_ion = 1.5d0 * (2.0d0 + eos%He_abundance) * saha_kB_cgs
        cv_neu = 1.5d0 * (1.0d0 + eos%He_abundance) * saha_kB_cgs
        if (eos%ionE) then
            T_lo = (eint_nH_cgs - saha_chi_H_cgs) / cv_ion
        else
            T_lo = eint_nH_cgs / cv_ion
        end if
        T_lo = max(T_lo, 100.0d0)
        T_hi = max(eint_nH_cgs / cv_neu, T_lo + 1.0d0)

        ! eint(T) is monotonic, so the bracket always holds; where hydrogen is already fully
        ! ionised at T_lo (to 1e-10), T_lo = (eint - chi_H)/cv_ion is the solution
        if (eos%ionE) then
            if (saha_y_cgs(nH_cgs, T_lo) > 1.0d0 - 1.0d-10) then
                T_out = T_lo / unit_temperature
                y_out = 1.0d0
                return
            end if
        end if

        ! Bisect in CGS
        do iter = 1, max_iter
            T_mid = 0.5d0 * (T_lo + T_hi)
            y_mid = saha_y_cgs(nH_cgs, T_mid)
            if (saha_eint_cgs(y_mid, T_mid) < eint_nH_cgs) then
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

    !> Build the analytic Gamma_1 table in eos%gamma1_p on (log10 nH, log10 p/nH), the axes of the
    !> LTE tables, so that gamma1_from_nH_p serves every method. Each node takes T and y from the
    !> Saha solve at (nH, p), and Gamma_1 = chi_rho + chi_T^2 (1+He+y)/c_v, with c_v per nH in kB
    subroutine build_gamma1_analytic_table()
        integer, parameter :: ng = 256
        integer :: i, j
        double precision :: log_nH_min, log_nH_max, log_pn_min, log_pn_max, T_min, T_max
        double precision :: dx1, dx2, nH_code, p_code, T_code, y, chi, dlny, npart, cv
        double precision :: chi_T, chi_rho, g1_val, g1_min, g1_max

        ! log10(nH [CGS]) from 5 to 19, and p/nH = (1+He+y) T over T from 10^2.8 K neutral to
        ! 10^6.3 K ionised, all in code units
        log_nH_min = 5.0d0  - dlog10(unit_numberdensity)
        log_nH_max = 19.0d0 - dlog10(unit_numberdensity)
        T_min = 10.0d0**2.8d0 / unit_temperature
        T_max = 10.0d0**6.3d0 / unit_temperature
        log_pn_min = dlog10((1.0d0 + eos%He_abundance) * T_min)
        log_pn_max = dlog10((2.0d0 + eos%He_abundance) * T_max)

        eos%gamma1_p%dim1 = ng
        eos%gamma1_p%dim2 = ng
        eos%gamma1_p%var1_min = log_nH_min
        eos%gamma1_p%var1_max = log_nH_max
        eos%gamma1_p%var2_min = log_pn_min
        eos%gamma1_p%var2_max = log_pn_max
        eos%gamma1_p%filename = 'computed_gamma1_analytic'
        allocate(eos%gamma1_p%table(ng, ng))

        dx1 = (log_nH_max - log_nH_min) / dble(ng - 1)
        dx2 = (log_pn_max - log_pn_min) / dble(ng - 1)
        g1_min = 1.0d30
        g1_max = -1.0d30
        do j = 1, ng
            do i = 1, ng
                nH_code = 10.0d0**(log_nH_min + (i - 1) * dx1)
                p_code  = nH_code * 10.0d0**(log_pn_min + (j - 1) * dx2)
                call saha_state_from_nH_p(nH_code, p_code, T_code, y)
                chi = saha_chi_kB / (T_code * unit_temperature)
                ! dln y/dln T at fixed nH from y^2/(1-y) ~ T^1.5 exp(-chi_H/kT)/nH
                dlny = (1.0d0 - y) / (2.0d0 - y) * (1.5d0 + chi)
                npart = 1.0d0 + eos%He_abundance + y
                cv = 1.5d0 * npart + 1.5d0 * y * dlny
                if (eos%ionE) cv = cv + y * dlny * chi
                chi_T = 1.0d0 + y * dlny / npart
                chi_rho = 1.0d0 - y * (1.0d0 - y) / ((2.0d0 - y) * npart)
                g1_val = chi_rho + chi_T**2 * npart / cv
                eos%gamma1_p%table(i, j) = g1_val
                g1_min = min(g1_min, g1_val)
                g1_max = max(g1_max, g1_val)
            end do
        end do
        call precompute_step_inv(eos%gamma1_p)

        if (mype == 0) then
            write(*, '(A,F8.4,A,F8.4)') &
                ' Gamma1 analytic table built (nH x p/nH): min = ', g1_min, ', max = ', g1_max
        end if

    end subroutine build_gamma1_analytic_table

end module mod_eos_lte_saha
!> Needs a line after to pass the preprocessor
