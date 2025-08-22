module mod_variables
  use mod_basic_types

  implicit none
  public

  !> Number of flux variables
  integer           :: nwflux = 0
  !$acc declare copyin(nwflux)

  !> Number of flux variables which need user to specify boundary type
  integer           :: nwfluxbc = 0
  !$acc declare copyin(nwfluxbc)

  !> Number of auxiliary variables in w
  integer           :: nwaux = 0
  !$acc declare copyin(nwaux)

  !> Number of extra variables in w
  integer           :: nwextra = 0
  !$acc declare copyin(nwextra)

  !> Number of extra variables in wextra seperated from w
  integer           :: nw_extra = 0
  !$acc declare copyin(nw_extra)

  !> Total number of variables
  integer           :: nw = 0
  !$acc declare copyin(nw)

  !> Total number of stagger variables
  integer           :: nws = 0
  !$acc declare copyin(nws)

  !> Number of variables which need to be updated in ghost cells
  integer           :: nwgc = 0
  !$acc declare copyin(nwgc)

  !> Number of vector variables (used for writing output)
  integer           :: nvector = 0
  !$acc declare copyin(nvector)

  !> Indices of vector variables
  integer, dimension(:), allocatable :: iw_vector
  !$acc declare create(iw_vector)

  ! the number of the first w variable to exchange ghost cells
  integer            :: iwstart=1
  !$acc declare copyin(iwstart)
  
  !> Maximum number of variables
  integer, parameter :: max_nw = 50
  !$acc declare copyin(max_nw)

  !> Primitive variable names
  character(len=name_len) :: prim_wnames(max_nw)

  !> Conservative variable names
  character(len=name_len) :: cons_wnames(max_nw)

  ! Global indices of variables that are often used

  !> Index of the (gas) density
  integer :: iw_rho = -1
  !$acc declare copyin(iw_rho)

  !> Indices of the momentum density
  integer, allocatable :: iw_mom(:)
  !$acc declare create(iw_mom)


  !> Index of the energy density
  integer :: iw_e = -1
  !$acc declare copyin(iw_e)

  !> Index of the heat flux
  integer :: iw_q = -1
  !$acc declare copyin(iw_q)

  !> Index of the radiation energy density
  integer :: iw_r_e = -1
  !$acc declare copyin(iw_r_e)


  !> Indices of the magnetic field components
  integer, allocatable, protected :: iw_mag(:)
  !$acc declare create(iw_mag)

  !> Index of the cutoff temperature for the TRAC method
  integer :: iw_tcoff = -1
  !$acc declare copyin(iw_tcoff)

  !> number of species: each species has different characterictic speeds and should
  !> be used accordingly in mod_finite_volume and mod_finite_difference
  integer :: number_species = 1
  !$acc declare copyin(number_species)

  !> index of the var
  !> whose velocity appears in the induction eq.
  integer :: index_v_mag = 1
  !$acc declare copyin(index_v_mag)

  !> the indices in 1:nwflux array are assumed consecutive for each species
  !> this array should be of size number_species and contain the first index in the array of
  !> the number_species
  integer, allocatable :: start_indices(:)
  !$acc declare create(start_indices)

  !> the indices in 1:nwflux array are assumed consecutive for each species
  !> this array should be of size number_species and contain the last index in the array of
  !> the first number_species, the last index for the last one is nwflux
  integer, allocatable :: stop_indices(:)
  !$acc declare create(stop_indices)


  ! indices of equi for the species index_v_mag
  ! these are needed for hlld solver, TODO: consider moving in a separate file
  integer :: iw_equi_rho = -1
  integer :: iw_equi_p = -1
  !$acc declare copyin(iw_equi_rho, iw_equi_p)

contains

  !> Set generic flux variable
  function var_set_fluxvar(name_cons, name_prim, ix, need_bc) result(iw)
    character(len=*), intent(in)  :: name_cons !< Conservative name
    character(len=*), intent(in)  :: name_prim !< Primitive name
    integer, intent(in), optional :: ix !< Optional index (to make var1, var2, ...)
    logical, intent(in), optional :: need_bc !< Require boundary condition (default: true)
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

   !$acc update device(nwflux, nw, nwfluxbc)
  end function var_set_fluxvar

  !> Set extra variable in w, which is not advected and has no boundary conditions.
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
   !$acc update device(nwextra,nw)
  end function var_set_extravar

  !> Set extra variable in wextra, which is not advected and has no boundary conditions and not output in dat.
  !> This has to be done after defining flux variables and auxiliary variables.
  function var_set_wextra() result(iw)
    integer :: iw

    nw_extra = nw_extra + 1
    iw      = nw_extra
   !$acc update device(nw_extra)
  end function var_set_wextra

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
   !$acc update device(nwaux,nw)
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
    !$acc update device(nwflux,nw,nwfluxbc,iw_rho)
  end function var_set_rho

  ! THE INCLUDE files cannot use other modules
  ! mpistop replaced by errormsg, should it exit?
  !> Exit MPI-AMRVAC with an error message
  subroutine errormsg(message)
  
    character(len=*), intent(in) :: message !< The error message
  
    write(*, *) "ERROR for processor"
    write(*, *) trim(message)
  
  
  end subroutine errormsg
  !> Set momentum variables
  function var_set_momentum(ndir) result(iw)
    integer, intent(in) :: ndir
    integer             :: iw(ndir), idir

    if (allocated(iw_mom)) call errormsg("Error: set_mom was already called")
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
    !$acc update device(nwflux,nw,nwfluxbc,iw_mom)
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
    !$acc update device(nwflux,nw,nwfluxbc,iw_e)
  end function var_set_energy

  !> Set heat flux variable (hyperbolic TC treatment)
  function var_set_q() result(iw)
    integer :: iw

    nwflux              = nwflux + 1
    nwfluxbc            = nwfluxbc + 1
    nw                  = nw + 1
    iw_q                = nwflux
    iw                  = nwflux
    cons_wnames(nwflux) = 'q'
    prim_wnames(nwflux) = 'q'
    !$acc update device(nwflux,nw,nwfluxbc,iw_q)
  end function var_set_q

  function var_set_radiation_energy() result(iw)
    integer :: iw

    nwflux              = nwflux + 1
    nwfluxbc            = nwfluxbc + 1
    nw                  = nw + 1
    iw_r_e              = nwflux
    iw                  = nwflux
    cons_wnames(nwflux) = 'r_e'
    prim_wnames(nwflux) = 'r_e'
    !$acc update device(nwflux,nw,nwfluxbc,iw_r_e)
  end function var_set_radiation_energy

  !> Set magnetic field variables
  function var_set_bfield(ndir) result(iw)
    integer, intent(in) :: ndir
    integer             :: iw(ndir), idir

    if (allocated(iw_mag)) call errormsg("Error: set_mag was already called")
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
   !$acc update device(nwflux,nw,nwfluxbc,iw_mag)
  end function var_set_bfield

end module mod_variables
