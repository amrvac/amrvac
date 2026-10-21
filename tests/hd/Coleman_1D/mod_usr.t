!> Coleman (2020) Riemann problem test suite -- FI and LTE variants.
!>
!> Six Riemann problems through hydrogen and helium ionisation zones,
!> exercising the EoS dispatch on (rho, eint) lookups for LTE and on the
!> ideal-gas closure for FI. Initial conditions are set via the &usr_list
!> namelist (rho_L, rho_R, v_L, v_R, T_L, T_R). Cell states are converted
!> to conserved variables through eos% to keep both FI and LTE paths in
!> sync with mod_eos's dispatch matrix.
!>
!> Output column g1 carries the runtime Gamma_1 (LTE: g1p table; FI: gamma).
!> The Python comparison harness validates against the exact tabular
!> Riemann solution and the Snow (2026) two-gamma jump conditions.

module mod_usr
  use mod_hd
  use mod_eos
  implicit none

  ! Left / right primitive state (code units), read from &usr_list
  double precision :: rho_L, rho_R
  double precision :: v_L,   v_R
  double precision :: T_L,   T_R

  ! Extra output variable index for Gamma_1
  integer :: g1_

contains

  subroutine usr_init()
    use mod_global_parameters

    ! Match Coleman (2020) hydrogen ionisation units (Table 1)
    unit_length        = 1.d9
    unit_temperature   = 1.d6
    unit_numberdensity = 1.d15

    usr_set_parameters => initglobaldata_usr
    usr_init_one_grid  => initonegrid_usr
    usr_modify_output  => set_output_vars

    call set_coordinate_system("Cartesian_1D")
    call hd_activate()

    g1_ = var_set_extravar("g1", "g1")

  end subroutine usr_init

  subroutine usr_params_read(files)
    character(len=*), intent(in) :: files(:)
    integer :: n

    namelist /usr_list/ rho_L, rho_R, v_L, v_R, T_L, T_R

    rho_L = 1.0d0;  rho_R = 0.125d0
    v_L   = 0.0d0;  v_R   = 0.0d0
    T_L   = 1.0d0;  T_R   = 0.8d0

    do n = 1, size(files)
      open(unitpar, file=trim(files(n)), status="old")
      read(unitpar, usr_list, end=109)
109   close(unitpar)
    end do

  end subroutine usr_params_read

  subroutine initglobaldata_usr
    use mod_global_parameters

    call usr_params_read(par_files)

    if (mype == 0) then
      print *, '=== Coleman Riemann Setup ==='
      print *, '  rho_L =', rho_L, '  rho_R =', rho_R
      print *, '  v_L   =', v_L,   '  v_R   =', v_R
      print *, '  T_L   =', T_L,   '  T_R   =', T_R
      print *, '  eos_type =', trim(eos%eos_type)
      print *, '  eos%He_abundance =', eos%He_abundance
    end if

  end subroutine initglobaldata_usr

  subroutine initonegrid_usr(ixI^L, ixO^L, w, x)
    use mod_global_parameters
    use mod_eos_LTE, only: eint_nH_from_T
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, 1:nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)

    double precision :: nH_val, log_nH, log_T, eint_nH_val
    double precision :: p_val, eint_val, T_local, rho_local, v_local
    double precision :: x_mid
    integer :: ix1

    x_mid = 0.5d0 * (xprobmin1 + xprobmax1)

    {do ix^DB = ixO^LIM^DB\}

      if (x(ix1, 1) < x_mid) then
        T_local = T_L; rho_local = rho_L; v_local = v_L
      else
        T_local = T_R; rho_local = rho_R; v_local = v_R
      end if

      w(ix1, iw_rho) = rho_local
      w(ix1, iw_mom(1)) = rho_local * v_local

      select case(eos%eos_type)
      case('LTE')
        nH_val  = rho_local / eos%nH2rhoFactor
        log_nH  = dlog10(nH_val)
        log_T   = dlog10(T_local)
        eint_nH_val = eint_nH_from_T(log_nH, log_T)
        eint_val = nH_val * eint_nH_val
        w(ix1, iw_e) = eint_val + 0.5d0 * rho_local * v_local**2
      case default
        ! FI: ideal gas with eq_state_units=false -> p = rho*T
        p_val = rho_local * T_local
        w(ix1, iw_e) = p_val * eos%inv_gamma_minus_1 &
                     + 0.5d0 * rho_local * v_local**2
      end select

    {end do\}

    call eos%update_eos(ixI^L, ixO^L, w, x)

  end subroutine initonegrid_usr

  subroutine set_output_vars(ixI^L, ixO^L, qt, w, x)
    use mod_global_parameters
    use mod_eos_LTE, only: gamma1_from_nH_p
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qt, x(ixI^S, 1:ndim)
    double precision, intent(inout) :: w(ixI^S, nw)

    double precision :: nH_val, log_nH, p_nH
    integer :: ix1

    call eos%update_eos(ixI^L, ixO^L, w, x)

    select case(eos%eos_type)
    case('LTE')
      {do ix^DB = ixO^LIM^DB\}
        nH_val = w(ix1, iw_rho) / eos%nH2rhoFactor
        log_nH = dlog10(nH_val)
        p_nH = (1.0d0 + eos%He_abundance + w(ix1, iw_ne)/nH_val) &
             * w(ix1, iw_te)
        w(ix1, g1_) = gamma1_from_nH_p(log_nH, dlog10(p_nH))
      {end do\}
    case default
      w(ixO^S, g1_) = eos%gamma
    end select

  end subroutine set_output_vars

end module mod_usr
