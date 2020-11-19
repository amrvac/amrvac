module mod_variables
  use mod_basic_types

  implicit none
  public

  !> Number of flux variables
  integer           :: nwflux = 0

  !> Number of flux variables which need user to specify boundary type
  integer           :: nwfluxbc = 0

  !> Number of auxiliary variables
  integer           :: nwaux = 0

  !> Number of extra variables
  integer           :: nwextra = 0

  !> Total number of variables
  integer           :: nw = 0

  !> Total number of stagger variables
  integer           :: nws = 0

  !> Number of variables which need to be updated in ghost cells
  integer           :: nwgc = 0

  !> Number of vector variables (used for writing output)
  integer           :: nvector = 0

  !> Indices of vector variables
  integer, dimension(:), allocatable :: iw_vector

  ! the number of the first w variable to exchange ghost cells
  integer            :: iwstart=1

  !> Maximum number of variables
  integer, parameter :: max_nw = 50

  !> Primitive variable names
  character(len=name_len) :: prim_wnames(max_nw)

  !> Conservative variable names
  character(len=name_len) :: cons_wnames(max_nw)

  ! Global indices of variables that are often used

  !> Index of the (gas) density
  integer, protected :: iw_rho = -1

  !> Indices of the momentum density
  integer, allocatable, protected :: iw_mom(:)

  !> Index of the energy density
  integer, protected :: iw_e = -1

  !> Index of the internal energy density
  integer, protected :: iw_eaux = -1

  !> Indices of the magnetic field components
  integer, allocatable, protected :: iw_mag(:)

  !> Index of the cutoff temperature for the TRAC method
  integer, protected :: iw_tcoff = -1

contains

  !> Set generic flux variable
  function var_set_fluxvar(name_cons, name_prim, ix, need_bc) result(iw)
    character(len=*), intent(in)  :: name_cons !< Conservative name
    character(len=*), intent(in)  :: name_prim !< Primitive name
    integer, intent(in), optional :: ix        !< Optional index (to make var1, var2, ...)
    logical, intent(in), optional :: need_bc   !< Require boundary condition (default: true)
    integer                       :: iw
    logical                       :: add_bc

    nwflux = nwflux + 1
    nw     = nw + 1
    iw     = nwflux

    add_bc = .true.
    if (present(need_bc)) add_bc = need_bc
    if (add_bc) nwfluxbc = nwfluxbc + 1

    if (.not. present(ix)) then
      prim_wnames(nwflux) = name_cons
      cons_wnames(nwflux) = name_prim
    else
      write(cons_wnames(nwflux),"(A,I0)") name_cons, ix
      write(prim_wnames(nwflux),"(A,I0)") name_prim, ix
    end if
  end function var_set_fluxvar

  !> Set extra variable, which is not advected and has no boundary conditions.
  !> This has to be done after defining flux variables and auxiliary variables.
  function var_set_extravar(name_cons, name_prim, ix) result(iw)
    character(len=*), intent(in)  :: name_cons, name_prim
    integer, intent(in), optional :: ix
    integer                       :: iw

    nwextra = nwextra + 1
    nw      = nw + 1
    iw      = nw

    if (.not. present(ix)) then
      prim_wnames(iw) = name_cons
      cons_wnames(iw) = name_prim
    else
      write(cons_wnames(iw),"(A,I0)") name_cons, ix
      write(prim_wnames(iw),"(A,I0)") name_prim, ix
    end if
  end function var_set_extravar

  !> Set auxiliary variable, which is not advected but has boundary conditions.
  !> This has to be done after defining flux variables.
  function var_set_auxvar(name_cons, name_prim, ix) result(iw)
    character(len=*), intent(in)  :: name_cons, name_prim
    integer, intent(in), optional :: ix
    integer                       :: iw

    nwaux   = nwaux + 1
    nw      = nw + 1
    iw      = nw

    if (.not. present(ix)) then
      prim_wnames(iw) = name_cons
      cons_wnames(iw) = name_prim
    else
      write(cons_wnames(iw),"(A,I0)") name_cons, ix
      write(prim_wnames(iw),"(A,I0)") name_prim, ix
    end if
  end function var_set_auxvar

  !> Set density variable
  function var_set_rho() result(iw)
    integer :: iw

    nwflux              = nwflux + 1
    nwfluxbc            = nwfluxbc + 1
    nw                  = nw + 1
    iw_rho              = nwflux
    iw                  = nwflux
    prim_wnames(nwflux) = 'rho'
    cons_wnames(nwflux) = 'rho'
  end function var_set_rho

  !> Set momentum variables
  function var_set_momentum(ndir) result(iw)
    integer, intent(in) :: ndir
    integer             :: iw(ndir), idir

    if (allocated(iw_mom)) &
         call mpistop("Error: set_mom was already called")
    allocate(iw_mom(ndir))

    do idir = 1, ndir
      nwflux       = nwflux + 1
      nwfluxbc     = nwfluxbc + 1
      nw           = nw + 1
      iw_mom(idir) = nwflux
      iw(idir)     = nwflux
      write(cons_wnames(nwflux),"(A1,I1)") "m", idir
      write(prim_wnames(nwflux),"(A1,I1)") "v", idir
    end do
  end function var_set_momentum

  !> Set energy variable
  function var_set_energy() result(iw)
    integer :: iw

    nwflux              = nwflux + 1
    nwfluxbc            = nwfluxbc + 1
    nw                  = nw + 1
    iw_e                = nwflux
    iw                  = nwflux
    cons_wnames(nwflux) = 'e'
    prim_wnames(nwflux) = 'p'
  end function var_set_energy

  !> Set magnetic field variables
  function var_set_bfield(ndir) result(iw)
    integer, intent(in) :: ndir
    integer             :: iw(ndir), idir

    if (allocated(iw_mag)) &
         call mpistop("Error: set_mag was already called")
    allocate(iw_mag(ndir))

    do idir = 1, ndir
      nwflux       = nwflux + 1
      nwfluxbc     = nwfluxbc + 1
      nw           = nw + 1
      iw_mag(idir) = nwflux
      iw(idir)     = nwflux
      write(cons_wnames(nwflux),"(A1,I1)") "b", idir
      write(prim_wnames(nwflux),"(A1,I1)") "b", idir
    end do
  end function var_set_bfield

  !> Set internal energy variable
  function var_set_internal_energy() result(iw)
    integer :: iw

    nwflux              = nwflux + 1
    nw                  = nw + 1
    iw_eaux             = nwflux
    iw                  = nwflux
    cons_wnames(nwflux) = 'eaux'
    prim_wnames(nwflux) = 'paux'
  end function var_set_internal_energy

  !> Set Tcoff variable for TRAC method
  function var_set_tcoff() result(iw)
    integer :: iw

    nwextra = nwextra + 1
    nw      = nw + 1
    iw      = nw
    iw_tcoff = nw
    cons_wnames(iw) = 'Tcoff'
    prim_wnames(iw) = 'Tcoff'
  end function var_set_tcoff

end module mod_variables
