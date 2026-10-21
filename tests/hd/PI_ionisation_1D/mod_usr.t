!> Partial-ionisation (eos_type='PI') 1D smoke test.
!>
!> A periodic box with a smooth sinusoidal temperature / density / velocity
!> profile whose temperature spans the hydrogen ionisation zone (~6000-15000 K
!> for the default parameters at unit_temperature=1e6). Evolving it for a few
!> steps exercises the full PI EoS path: to_primitive / to_conserved / p_to_e,
!> get_thermal_pressure, get_csound2, get_Rfactor and the per-substep Te_ update,
!> for both the energy (ionE=T) and no-energy (ionE=F) closures and for each
!> pi_table backend. Purely a regression lock for the PI rewrite -- there is no
!> analytic reference; the test compares its own .log against correct_output/.
module mod_usr
  use mod_hd
  use mod_eos
  implicit none

  double precision :: rho_base, T_base, ampR, ampT, v_amp
  integer :: te_index

contains

  subroutine usr_init()
    use mod_global_parameters

    ! Chromospheric normalisation (matches Coleman_1D units).
    unit_length        = 1.d9
    unit_temperature   = 1.d6
    unit_numberdensity = 1.d15

    usr_set_parameters => initglobaldata_usr
    usr_init_one_grid  => initonegrid_usr

    call set_coordinate_system("Cartesian_1D")
    call hd_activate()

  end subroutine usr_init

  subroutine usr_params_read(files)
    character(len=*), intent(in) :: files(:)
    integer :: n

    namelist /usr_list/ rho_base, T_base, ampR, ampT, v_amp

    rho_base = 1.0d0;  T_base = 0.011d0
    ampR     = 0.2d0;  ampT   = 0.4d0
    v_amp    = 0.2d0

    do n = 1, size(files)
      open(unitpar, file=trim(files(n)), status="old")
      read(unitpar, usr_list, end=111)
111   close(unitpar)
    end do

  end subroutine usr_params_read

  subroutine initglobaldata_usr
    use mod_global_parameters
    use mod_variables, only: cons_wnames
    integer :: iw

    call usr_params_read(par_files)

    ! Locate the PI temperature aux-var by name so the IC can seed it
    ! (the no-energy Te_ update lags via wCT(Te_); it must start consistent).
    ! Done by name to stay portable across the PI implementation.
    te_index = -1
    do iw = 1, nw
      if (trim(cons_wnames(iw)) == 'Te') te_index = iw
    end do

    if (mype == 0) then
      print *, '=== PI ionisation smoke ==='
      print *, '  eos_type =', trim(eos%eos_type), '  pi_table =', trim(eos%pi_table), &
               '  ionE =', eos%ionE
      print *, '  He_abundance =', eos%He_abundance, '  Te index =', te_index
    end if

  end subroutine initglobaldata_usr

  subroutine initonegrid_usr(ixI^L, ixO^L, w, x)
    use mod_global_parameters
    use mod_eos_PI, only: p_eint_from_rho_T_PI
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, 1:nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)

    double precision :: rho_l, T_l, v_l, p_l, eint_l, Rf, twopi
    integer :: ix1

    twopi = 2.0d0 * dpi

    {do ix^DB = ixO^LIM^DB\}
      rho_l = rho_base * (1.0d0 + ampR * dcos(twopi * x(ix1, 1)))
      T_l   = T_base   * (1.0d0 + ampT * dsin(twopi * x(ix1, 1)))
      v_l   = v_amp    * dsin(twopi * x(ix1, 1))

      w(ix1, iw_rho)    = rho_l
      w(ix1, iw_mom(1)) = rho_l * v_l

      ! PI forward map: (rho, T) -> (p, gas eint, R), handling ionE internally.
      call p_eint_from_rho_T_PI(rho_l, T_l, p_l, eint_l, Rf)
      w(ix1, iw_e) = eint_l + 0.5d0 * rho_l * v_l**2

      if (te_index > 0) w(ix1, te_index) = T_l
    {end do\}

  end subroutine initonegrid_usr

end module mod_usr
