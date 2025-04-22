!> Hydrodynamics physics module
module mod_hd_phys
  use mod_physics
  use mod_comm_lib, only: mpistop
  implicit none
  private

  !> Whether an energy equation is used
  logical, public, protected              :: hd_energy = .true.
  !$acc declare copyin(hd_energy)

  !> Index of the density (in the w array)
  integer, public, protected              :: rho_
  !$acc declare create(rho_)

  !> Indices of the momentum density
  integer, allocatable, public, protected :: mom(:)
  !$acc declare create(mom)

  !> Index of the energy density (-1 if not present)
  integer, public, protected              :: e_
  !$acc declare create(e_)

  !> Index of the gas pressure (-1 if not present) should equal e_
  integer, public, protected              :: p_
  !$acc declare create(p_)

  !> The adiabatic index
  double precision, public                :: hd_gamma = 5.d0/3.0d0
  !$acc declare copyin(hd_gamma)

  !> The adiabatic constant
  double precision, public                :: hd_adiab = 1.0d0
  !$acc declare copyin(hd_adiab)

  !> Whether plasma is partially ionized
  logical, public, protected              :: hd_partial_ionization = .false.
  !$acc declare copyin(hd_partial_ionization)

  !> Allows overruling default corner filling (for debug mode, since otherwise corner primitives fail)
  logical, public, protected              :: hd_force_diagonal = .false.
  !$acc declare copyin(hd_force_diagonal)
  
  !> Whether particles module is added
  logical, public, protected              :: hd_particles = .false.
  !$acc declare copyin(hd_particles)



  ! Public methods
  public :: hd_read_params
  public :: hd_add_source
  public :: hd_get_cmax_scalar
  
contains

  !> Read this module's parameters from a file
  subroutine hd_read_params(files)
    use mod_global_parameters
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /hd_list/ hd_energy, hd_gamma, hd_adiab, hd_partial_ionization, &
    hd_force_diagonal, hd_particles

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, hd_list, end=111)
111    close(unitpar)
    end do

#ifdef _OPENACC
    !$acc update device(hd_energy, hd_gamma, hd_adiab, hd_partial_ionization, &
    hd_force_diagonal, hd_particles)
#endif
    
  end subroutine hd_read_params

  ! w[iw]= w[iw]+qdt*S[wCT, qtC, x] where S is the source based on wCT within ixO
  subroutine hd_add_source(qdt,dtfactor, ixI^L,ixO^L,wCT,wCTprim,w,x,qsourcesplit,active)
    !$acc routine seq
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, dtfactor
    double precision, intent(in)    :: wCT(ixI^S, 1:nw),wCTprim(ixI^S,1:nw), x(ixI^S, 1:ndim)
    double precision, intent(inout) :: w(ixI^S, 1:nw)
    logical, intent(in)             :: qsourcesplit
    logical, intent(inout)          :: active

  end subroutine hd_add_source

  pure subroutine hd_get_cmax_scalar(w, idim, cmax)
    !$acc routine seq

    integer, intent(in)                       :: idim
    double precision, intent(in)              :: w(:)
    double precision, intent(out)             :: cmax
    double precision                          :: csound
    double precision                          :: v

    v = w(mom(idim)) / w(rho_)

    csound = (hd_gamma - 1.0d0) * (w(e_) - &     
         0.5d0 * ({^D& w(mom(^D))**2|+}) / w(rho_) ) ! p
    
    csound = hd_gamma * csound / w(rho_) ! cs**2
    csound = sqrt(csound)

    cmax = abs(v) + csound
    
  end subroutine hd_get_cmax_scalar

  !> Initialize the module
  subroutine hd_phys_init()
    use mod_global_parameters
    use mod_particles, only: particles_init


    call hd_read_params(par_files)

    physics_type = "hd"
    phys_energy  = hd_energy
    phys_total_energy  = hd_energy
    phys_internal_e = .false.
    phys_gamma = hd_gamma
    phys_partial_ionization=hd_partial_ionization
    !$acc update device(physics_type, phys_energy, phys_total_energy, phys_internal_e, phys_gamma, phys_partial_ionization)

    use_particles = hd_particles
    
    ! Determine flux variables
    rho_ = var_set_rho()
    !$acc update device(rho_)
    
    allocate(mom(ndir))
    mom(:) = var_set_momentum(ndir)
    !$acc update device(mom)

    ! Set index of energy variable
    if (hd_energy) then
       e_ = var_set_energy()
       p_ = e_
    else
       e_ = -1
       p_ = -1
    end if
    !$acc update device(e_,p_)
    
    ! Whether diagonal ghost cells are required for the physics
    phys_req_diagonal = .false.

    if (hd_force_diagonal) then
       ! ensure corners are filled, otherwise divide by zero when getting primitives
       !  --> only for debug purposes
       phys_req_diagonal = .true.
    endif

    ! set number of variables which need update ghostcells
    nwgc=nwflux
    !$acc update device(nwgc)
    
    ! Initialize particles module
    if (hd_particles) then
       call particles_init()
       phys_req_diagonal = .true.
    end if

    nvector      = 1 ! No. vector vars
    allocate(iw_vector(nvector))
    iw_vector(1) = mom(1) - 1
    !$acc update device(nvector, iw_vector, flux_type)
    !$acc update device(phys_req_diagonal)
    
  end subroutine hd_phys_init
  
end module mod_hd_phys
