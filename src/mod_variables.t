module mod_variables
  use mod_basic_types

  implicit none
  public

  !> Maximum number of variables
  integer, parameter :: max_nvar = 500

  integer, parameter :: hydro_var  = 1
  integer, parameter :: rad_var    = 2
  integer, parameter :: metric_var = 3

  ! ------------------
  ! Metric variables
  ! ------------------
  !> Number of metric variables
  integer            :: nmetric = 0
  !> index of min metric variables
  integer            :: nmetric_lo = 0
  !> index of max metric variables
  integer            :: nmetric_hi = 0

  ! ------------------
  ! hydro variables
  ! ------------------
  !> Number of flux (conservative) variables
  integer :: nwflux = 0
  !> Number of extra variables in w
  integer :: nwextra = 0
  !> Total number of flux variables
  integer :: nw = 0

  !> number of species: each species has different characterictic speeds and should
  !> be used accordingly in mod_finite_volume and mod_finite_difference
  integer :: number_species = 1

  !> the indices in 1:nwflux array are assumed consecutive for each species
  !> this array should be of size number_species and contain the first index in the array of
  !> the number_species
  integer, allocatable :: start_indices(:)
  !> the indices in 1:nwflux array are assumed consecutive for each species
  !> this array should be of size number_species and contain the last index in the array of
  !> the first number_species, the last index for the last one is nwflux
  integer, allocatable :: stop_indices(:)

  ! fixme: the lo hi seems useless now, remove them
  !> Number of conserved hydro variables
  integer            :: nc_hydro = 0
  !> index of min conserved hydro variables
  integer            :: nc_hydro_lo = 0
  !> index of max conserved hydro variables
  integer            :: nc_hydro_hi = 0

  !> Number of conserved rad variables
  integer            :: nc_rad = 0
  !> index of min conserved rad variables
  integer            :: nc_rad_lo = 0
  !> index of max conserved rad variables
  integer            :: nc_rad_hi = 0

  !> number of photons or neutrino species
  integer            :: n_spec = 0
  !> number of energy bins
  integer            :: n_ebin = 0

  !> Number of stagger conservative variables
  integer            :: nws = 0

  !> Number of vector variables (used for writing output)
  integer            :: nvector = 0

  !> Indices of vector variables
  integer            :: iw_vector(max_nvar)

  !> Number of auxiliary variables
  integer :: nwaux

  !> Number of auxiliary variables that are only included in output
  integer :: nwauxio

  !> Types of variables
  integer            :: type_vars(max_nvar)

  logical            :: var_reconstruct(max_nvar) ! fixme: remove me

  !> Number of extra variables in wextra seperated from w
  integer :: nw_extra = 0

  !> Primitive variable names
  character(len=name_len) :: metric_names(16)

  !> Primitive variable names
  character(len=name_len) :: prim_names(max_nvar)

  !> Conservative variable names
  character(len=name_len) :: cons_names(max_nvar)

  !> extra variable names
  character(len=name_len) :: wextra_names(max_nvar)

contains

  !> Set primitive variable
  function var_set_metricvar(name_metric, ix) result(iw)
    character(len=*), intent(in)  :: name_metric !< Primitive name
    integer, intent(in), optional :: ix        !< Optional index (to make var1, var2, ...)
    integer                       :: iw

    ! total number of metric variables
    nmetric = nmetric + 1
    iw      = nmetric

    if (.not. present(ix)) then
      metric_names(nmetric) = name_metric
    else
      write(metric_names(nmetric),"(A,I0)") name_metric, ix
    end if
  end function var_set_metricvar

  !> Set generic flux variable
  function var_set_fluxvar(name_cons, name_prim, ix, var_type) result(iw)
    character(len=*), intent(in)  :: name_cons !< Conservative name
    character(len=*), intent(in)  :: name_prim !< Primitive name
    integer, intent(in), optional :: ix        !< Optional index (to make var1, var2, ...)
    integer, intent(in), optional :: var_type  !< variable type
    integer                       :: iw, var_t

    if ( nwextra > 0 ) &
       error stop "wextra has to be done after defining flux variables and auxiliary variables."

    ! total number of w variables
    nwflux = nwflux + 1
    nw     = nw + 1
    iw     = nw

    if (.not. present(ix)) then
      prim_names(nw) = name_prim
      cons_names(nw) = name_cons
    else
      write(prim_names(nw),"(A,I0)") name_prim, ix
      write(cons_names(nw),"(A,I0)") name_cons, ix
    end if

    var_t = -1
    if (present(var_type)) var_t = var_type
    type_vars(iw) = var_t

    select case (var_t)
    case (hydro_var)
       if ( nc_hydro_lo == 0 ) nc_hydro_lo = iw
       nc_hydro = nc_hydro + 1
       nc_hydro_hi = iw
    case (rad_var)
       if ( nc_rad_lo == 0 ) nc_rad_lo = iw
       nc_rad = nc_rad + 1
       nc_rad_hi = iw
    case default
       ! nothing to do here, this is an aux variable
    end select
  end function var_set_fluxvar

  !> Set extra variable in w, which is not advected and has no boundary conditions.
  !> This has to be done after defining flux variables and auxiliary variables.
  function var_set_extravar(name_vars, ix) result(iw)
    character(len=*), intent(in)  :: name_vars
    integer, intent(in), optional :: ix
    integer                       :: iw

    nwextra = nwextra + 1
    nw      = nw + 1
    iw      = nw

    if (.not. present(ix)) then
      prim_names(iw) = name_vars
      cons_names(iw) = name_vars
    else
      write(cons_names(iw),"(A,I0)") name_vars, ix
      write(prim_names(iw),"(A,I0)") name_vars, ix
    end if
  end function var_set_extravar

  !> Set extra variable in wextra, which is not advected and has no boundary conditions and not output in dat.
  !> This has to be done after defining flux variables and auxiliary variables.
  function var_set_wextra(name_vars, ix) result(iw)
    character(len=*), intent(in)  :: name_vars
    integer, intent(in), optional :: ix
    integer                       :: iw

    nw_extra = nw_extra + 1
    iw       = nw_extra

    if (.not. present(ix)) then
      wextra_names(iw) = name_vars
    else
      write(wextra_names(iw),"(A,I0)") name_vars, ix
    end if
  end function var_set_wextra

end module mod_variables
