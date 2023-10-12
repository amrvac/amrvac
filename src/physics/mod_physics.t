!> This module defines the procedures of a physics module. It contains function
!> pointers for the various supported routines. An actual physics module has to
!> set these pointers to its implementation of these routines.
module mod_physics
  use mod_global_parameters, only: name_len, max_nvar
  use mod_physics_hllc
  use mod_physics_roe
  use mod_physics_ppm

  implicit none
  public

  !> String describing the physics type of the simulation
  character(len=name_len) :: physics_type = ""

  !> To use wider stencils in flux calculations. A value of 1 will extend it by
  !> one cell in both directions, in any dimension
  integer :: phys_wider_stencil    = 0
  integer :: phys_extra_ghostcells = 0

  !> Whether the physics routines require diagonal ghost cells, for example for
  !> computing a curl.
  logical :: phys_req_diagonal = .true.

  !> Array per direction per variable, which can be used to specify that certain
  !> fluxes have to be treated differently
  integer, allocatable :: flux_type(:, :)

  !-------------------------------------------------------------------!
  !          Primitive variables (-1 if not present)
  !-------------------------------------------------------------------!

  ! --------- Hydro variables --------- !
  !> Index of the density (in the prim array)
  integer :: rho_ = -1
  !> Indices of the W * velocities ( W = 1 in the newtonian case )
  integer, allocatable :: W_vel(:)
  !> Index of the internal energy (e = rho*eps) or the natural log temperature.
  !> i.e. when tabulated equation of state is use, this index will be used as natural log temperature, otherwise it is internal energy
  integer :: e_ = -1
  !> Index of the electronic frac
  integer :: ye_ = -1

  ! extra variables
  !> Index of the natural log temperature
  integer :: logtemp_ = -1
  !> Index of the lorentz factor W
  integer :: W_ = -1

  ! --------- MHD variables --------- !
  !> Indices of the E field
  integer, allocatable :: Evec(:)
  !> Indices of the B field
  integer, allocatable :: Bvec(:)
  !> Index of the scalar field when using GLM
  integer :: Bphi_ = -1
  integer :: Ephi_ = -1


  !-------------------------------------------------------------------!
  !          Conserved variables (-1 if not present)
  !-------------------------------------------------------------------!
  !> Indicate where is the Density
  integer :: D_ = -1
  !> Indicate where is the electron frac
  integer :: Dye_ = -1
  !> Indices of the momentum density
  integer, allocatable :: mom(:)
  !> Indicate where is the Energy
  integer :: tau_ = -1
  !> Indices of the Econs field
  integer, allocatable :: Econs(:)
  !> Indices of the Bcons field
  integer, allocatable :: Bcons(:)
  !> Index of the scalar field when using GLM
  integer :: Bphi_cons_ = -1
  integer :: Ephi_cons_ = -1

  !-------------------------------------------------------------------!

  !> omitted
  integer, parameter   :: flux_nul              = -1
  !> Indicates a normal flux
  integer, parameter   :: flux_default          = 0
  !> Indicates the flux should be treated with tvdlf
  integer, parameter   :: flux_tvdlf            = 1
  !> Indicates dissipation should be omitted
  integer, parameter   :: flux_no_dissipation   = 2
  !> Indicates the flux should be treated as asymptotic to diffusion limit for energy density
  integer, parameter   :: flux_asym_diffusion_e = 3
  !> Indicates the flux should be treated as asymptotic to diffusion limit for momentum density
  integer, parameter   :: flux_asym_diffusion_f = 4

  !> Type for special methods defined per variable
  type iw_methods
    integer :: test
    !> If this is set, use the routine as a capacity function when adding fluxes
    procedure(sub_get_var), pointer, nopass :: inv_capacity => null()
  end type iw_methods

  !> Special methods defined per variable
  type(iw_methods) :: phys_iw_methods(max_nvar)

  ! hydrodynamics
  procedure(sub_check_params), pointer        :: phys_check_params        => null()
  procedure(sub_convert), pointer             :: phys_to_conserved        => null()
  procedure(sub_convert), pointer             :: phys_to_primitive        => null()
  procedure(sub_modify_wLR), pointer          :: phys_modify_wLR          => null()
  procedure(sub_modify_tvdlfeps), pointer     :: phys_modify_tvdlfeps     => null()
  procedure(sub_get_var_from_state), pointer  :: phys_get_lfac2           => null()
  procedure(sub_get_var_from_state), pointer  :: phys_get_tilde_U         => null()
  procedure(sub_get_vec_from_state), pointer  :: phys_get_tilde_S_i       => null()
  procedure(sub_get_var_from_state), pointer  :: phys_get_tilde_S         => null()
  procedure(sub_get_cmax), pointer            :: phys_get_cmax            => null()
  procedure(sub_get_cbounds), pointer         :: phys_get_cbounds         => null()
  procedure(sub_get_cbounds_and_vct), pointer :: phys_get_cbounds_and_vct => null()
  procedure(sub_get_flux), pointer            :: phys_get_flux            => null()
  procedure(sub_get_dt), pointer              :: phys_get_dt              => null()
  procedure(sub_add_source_geom), pointer     :: phys_add_source_geom     => null()
  procedure(sub_add_source), pointer          :: phys_add_source          => null()
  procedure(sub_global_source), pointer       :: phys_global_source_after => null()
  procedure(sub_special_advance), pointer     :: phys_special_advance     => null()
  procedure(sub_check_prim), pointer          :: phys_check_prim          => null()
  procedure(sub_boundary_adjust), pointer     :: phys_boundary_adjust     => null()
  procedure(sub_write_info), pointer          :: phys_write_info          => null()
  procedure(sub_small_values), pointer        :: phys_handle_small_values => null()
  procedure(sub_get_var_from_state), pointer  :: phys_get_divb            => null()
  procedure(sub_get_var), pointer             :: phys_get_dive            => null()
  procedure(sub_update_faces), pointer        :: phys_update_faces        => null()
  procedure(sub_face_to_center), pointer      :: phys_face_to_center      => null()
  procedure(sub_implicit_update), pointer     :: phys_implicit_update     => null()
  procedure(sub_evaluate_implicit),pointer    :: phys_evaluate_implicit   => null()

  procedure(process_global), pointer          :: phys_clean_divb          => null()
  procedure(process_global), pointer          :: phys_initial_clean_divb  => null()

  ! Called at the beginning of every time step (after determining dt)
  procedure(process_grid), pointer            :: phys_process_grid        => null()
  procedure(process_global), pointer          :: phys_process_global      => null()
  ! Called every time step just after advance (with w^(n+1), it^n, global_time^n)
  procedure(process_grid), pointer            :: phys_process_adv_grid    => null()
  procedure(process_global), pointer          :: phys_process_adv_global  => null()

  abstract interface

     subroutine sub_check_params
     end subroutine sub_check_params

     subroutine sub_boundary_adjust(igrid,psb)
       use mod_global_parameters
       integer, intent(in) :: igrid
       type(state), target :: psb(max_blocks)
     end subroutine sub_boundary_adjust

     subroutine sub_convert(ixI^L, ixO^L, s)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L, ixO^L
       type(state), intent(inout)      :: s
     end subroutine sub_convert

     subroutine sub_modify_wLR(ixI^L, ixO^L, qt, idir, s, sL, sR, sLp, sRp)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L, ixO^L, idir
       double precision, intent(in)    :: qt
       type(state), intent(inout)      :: s, sL, sR, sLp, sRp
     end subroutine sub_modify_wLR

     subroutine sub_modify_tvdlfeps(ixI^L, ixO^L, idir, s, tvdlf_eps)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L, ixO^L, idir
       type(state)                     :: s
       double precision, intent(inout) :: tvdlf_eps(ixI^S, 1:nwflux)
     end subroutine sub_modify_tvdlfeps

     subroutine sub_get_cmax(s, ixI^L, ixO^L, idim, cmax)
       use mod_global_parameters
       type(state), intent(in)         :: s
       integer, intent(in)             :: ixI^L, ixO^L, idim
       double precision, intent(inout) :: cmax(ixI^S)
     end subroutine sub_get_cmax

     subroutine sub_get_cbounds(sL, sR, sLp, sRp, ixI^L, ixO^L, idim, cmax, cmin)
       use mod_global_parameters
       type(state), intent(in)         :: sL, sR, sLp, sRp
       integer, intent(in)             :: ixI^L, ixO^L, idim
       double precision, intent(inout)           :: cmax(ixI^S, number_species)
       double precision, intent(inout), optional :: cmin(ixI^S, number_species)
     end subroutine sub_get_cbounds

     subroutine sub_get_cbounds_and_vct(sL, sR, sLp, sRp, ixI^L, ixO^L, idim, vcts, cmax, cmin)
       use mod_global_parameters
       type(state), intent(in)         :: sL, sR, sLp, sRp
       integer, intent(in)             :: ixI^L, ixO^L, idim
       type(ct_velocity), intent(inout):: vcts
       double precision, intent(inout)           :: cmax(ixI^S, number_species)
       double precision, intent(inout), optional :: cmin(ixI^S, number_species)
     end subroutine sub_get_cbounds_and_vct

     subroutine sub_get_flux(s, sp, ixI^L, ixO^L, idim, f)
       use mod_global_parameters
       type(state), intent(in)         :: s, sp
       integer, intent(in)             :: ixI^L, ixO^L, idim
       double precision, intent(out)   :: f(ixI^S, 1:nwflux)
     end subroutine sub_get_flux

     subroutine sub_add_source_geom(qdt, s, sp, ixI^L, ixO^L, w)
       use mod_global_parameters
       double precision, intent(in)    :: qdt
       type(state), intent(in)         :: s, sp
       integer, intent(in)             :: ixI^L, ixO^L
       double precision, intent(inout) :: w(ixI^S, 1:nwflux)
     end subroutine sub_add_source_geom

     subroutine sub_add_source(qdt, s, sold, ixI^L, ixO^L, w, qsourcesplit, active)
       use mod_global_parameters
       double precision, intent(in)    :: qdt
       type(state), intent(in)         :: s, sold
       integer, intent(in)             :: ixI^L, ixO^L
       double precision, intent(inout) :: w(ixI^S, 1:nwflux)
       logical, intent(in)             :: qsourcesplit
       logical, intent(inout)          :: active !< Needs to be set to true when active
     end subroutine sub_add_source

     !> Add global source terms on complete domain (potentially implicit)
     subroutine sub_global_source(qdt, qt, active)
       double precision, intent(in) :: qdt    !< Current time step
       double precision, intent(in) :: qt     !< Current time
       logical, intent(inout)       :: active !< Output if the source is active
     end subroutine sub_global_source

     !> Add special advance in each advect step
     subroutine sub_special_advance(qt, psa)
       use mod_global_parameters
       double precision, intent(in) :: qt     !< Current time
       type(state), target :: psa(max_blocks) !< Compute based on this state
     end subroutine sub_special_advance

     subroutine sub_get_dt(s, ixI^L, ixO^L, dtnew)
       use mod_global_parameters
       type(state), intent(in)         :: s
       integer, intent(in)             :: ixI^L, ixO^L
       double precision, intent(inout) :: dtnew
     end subroutine sub_get_dt

     subroutine sub_check_prim(s, ixI^L, ixO^L, w_flag)
       use mod_global_parameters
       type(state), intent(in)      :: s
       integer, intent(in)          :: ixI^L, ixO^L
       integer, intent(inout)       :: w_flag(ixG^T)
     end subroutine sub_check_prim

     subroutine sub_write_info(file_handle)
       integer, intent(in) :: file_handle
     end subroutine sub_write_info

     subroutine sub_small_values(s, ixI^L, ixO^L, subname)
       use mod_global_parameters
       type(state), intent(inout)      :: s
       integer, intent(in)             :: ixI^L,ixO^L
       character(len=*), intent(in)    :: subname
     end subroutine sub_small_values

     subroutine sub_get_var_from_state(ixI^L, ixO^L, ps_in, out)
       use mod_global_parameters
       integer, intent(in)           :: ixI^L, ixO^L
       type(state), intent(in)       :: ps_in
       double precision, intent(out) :: out(ixI^S)
     end subroutine sub_get_var_from_state

     subroutine sub_get_vec_from_state(ixI^L, ixO^L, ps_in, out)
       use mod_global_parameters
       integer, intent(in)           :: ixI^L, ixO^L
       type(state), intent(in)       :: ps_in
       double precision, intent(out) :: out(ixI^S, 1:3)
     end subroutine sub_get_vec_from_state

     subroutine sub_get_var(ps_in, ixI^L, ixO^L, out)
       use mod_global_parameters
       type(state), intent(in)       :: ps_in
       integer, intent(in)           :: ixI^L, ixO^L
       double precision, intent(out) :: out(ixI^S)
     end subroutine sub_get_var

     subroutine sub_update_faces(ixI^L,ixO^L,qt,qdt,fC,fE,sCT,s,vcts)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L, ixO^L
       double precision, intent(in)    :: qt,qdt
       ! velocity structure
       type(state)                     :: sCT, s
       type(ct_velocity)               :: vcts
       double precision, intent(in)    :: fC(ixI^S,1:nwflux,1:ndim)
       double precision, intent(inout) :: fE(ixI^S,7-2*ndim:3)
     end subroutine sub_update_faces

     subroutine sub_face_to_center(ixO^L,s)
       use mod_global_parameters
       integer, intent(in)                :: ixO^L
       type(state)                        :: s
     end subroutine sub_face_to_center

     !> for processing after the advance (PIC-MHD, e.g.)
     subroutine process_grid(igrid,level,ixI^L,ixO^L,qt,w,x)
       use mod_global_parameters
       integer, intent(in)             :: igrid,level,ixI^L,ixO^L
       double precision, intent(in)    :: qt,x(ixI^S,1:ndim)
       double precision, intent(inout) :: w(ixI^S,1:nwflux)
     end subroutine process_grid
 
     !> for processing after the advance (PIC-MHD, e.g.)
     subroutine process_global(iit,qt)
       use mod_global_parameters
       integer, intent(in)          :: iit
       double precision, intent(in) :: qt
     end subroutine process_global

     subroutine sub_evaluate_implicit(qtC,psa)
       use mod_global_parameters
       type(state), target :: psa(max_blocks)   
       double precision, intent(in) :: qtC      
     end subroutine sub_evaluate_implicit

     subroutine sub_implicit_update(dtfactor,qdt,qtC,psa,psb)
       use mod_global_parameters
       type(state), target :: psa(max_blocks)   
       type(state), target :: psb(max_blocks)   
       double precision, intent(in) :: qdt
       double precision, intent(in) :: qtC      
       double precision, intent(in) :: dtfactor 
     end subroutine sub_implicit_update

  end interface

contains

  subroutine phys_check()
    use mod_global_parameters, only: nw, ndir
    use mod_comm, only: mpistop

    use mod_physics_hllc, only: phys_hllc_check
    use mod_physics_roe, only: phys_roe_check
    use mod_physics_ppm, only: phys_ppm_check

    if (physics_type == "") call mpistop("Error: no physics module loaded")

    call phys_hllc_check()
    call phys_roe_check()
    call phys_ppm_check()

    ! Checks whether the required physics methods have been defined
    if (.not. associated(phys_check_params)) &
         phys_check_params => dummy_check_params

    if (.not. associated(phys_to_conserved)) &
         call mpistop("Error: phys_to_conserved not defined")

    if (.not. associated(phys_to_primitive)) &
         call mpistop("Error: phys_to_primitive not defined")

    if (.not. associated(phys_modify_wLR)) &
         phys_modify_wLR => dummy_modify_wLR

    if (.not. associated(phys_modify_tvdlfeps)) &
         phys_modify_tvdlfeps => dummy_modify_tvdlfeps

    if (.not. associated(phys_get_cmax)) &
         call mpistop("Error: no phys_get_cmax not defined")

    if (.not. associated(phys_get_cbounds)) &
         call mpistop("Error: no phys_get_cbounds not defined")

    if (.not. associated(phys_get_flux)) &
         call mpistop("Error: no phys_get_flux not defined")

    if (.not. associated(phys_get_dt)) &
         call mpistop("Error: no phys_get_dt not defined")

    if (.not. associated(phys_add_source_geom)) &
         phys_add_source_geom => dummy_add_source_geom

    if (.not. associated(phys_add_source)) &
         phys_add_source => dummy_add_source

    if (.not. associated(phys_check_prim)) &
         phys_check_prim => dummy_check_prim

    if (.not. associated(phys_boundary_adjust)) &
         phys_boundary_adjust => dummy_boundary_adjust

    if (.not. associated(phys_write_info)) &
         phys_write_info => dummy_write_info

    if (.not. associated(phys_handle_small_values)) &
         phys_handle_small_values => dummy_small_values

    if (.not. associated(phys_update_faces)) &
         phys_update_faces => dummy_update_faces

    if (.not. associated(phys_face_to_center)) &
         phys_face_to_center => dummy_face_to_center

  end subroutine phys_check

  subroutine dummy_init_params
  end subroutine dummy_init_params

  subroutine dummy_check_params
  end subroutine dummy_check_params

  subroutine dummy_convert(ixI^L, ixO^L, s)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    type(state), intent(inout)      :: s
  end subroutine dummy_convert

  subroutine dummy_modify_wLR(ixI^L, ixO^L, qt, idir, s, sL, sR, sLp, sRp)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L, idir
    double precision, intent(in)    :: qt
    type(state), intent(inout)      :: s, sL, sR, sLp, sRp
  end subroutine dummy_modify_wLR

  subroutine dummy_modify_tvdlfeps(ixI^L, ixO^L, idir, s, tvdlf_eps)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L, idir
    type(state)                     :: s
    double precision, intent(inout) :: tvdlf_eps(ixI^S, nwflux)
  end subroutine dummy_modify_tvdlfeps

  subroutine dummy_add_source_geom(qdt, s, sp, ixI^L, ixO^L, w)
    use mod_global_parameters
    double precision, intent(in)    :: qdt
    type(state), intent(in)         :: s, sp
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, 1:nwflux)
  end subroutine dummy_add_source_geom

  subroutine dummy_add_source(qdt, s, sold, ixI^L, ixO^L, w, qsourcesplit, active)
    use mod_global_parameters
    double precision, intent(in)    :: qdt
    type(state), intent(in)         :: s, sold
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, 1:nwflux)
    logical, intent(in)             :: qsourcesplit
    logical, intent(inout)          :: active
    ! Don't have to set active, since it starts as .false.
  end subroutine dummy_add_source

  subroutine dummy_check_prim(s, ixI^L, ixO^L, w_flag)
    use mod_global_parameters
    type(state), intent(in)      :: s
    integer, intent(in)          :: ixI^L, ixO^L
    integer, intent(inout)       :: w_flag(ixG^T)
    w_flag(ixO^S) = 0             ! All okay
  end subroutine dummy_check_prim

  subroutine dummy_boundary_adjust(igrid,psb)
    use mod_global_parameters
    integer, intent(in) :: igrid
    type(state), target :: psb(max_blocks)
  end subroutine dummy_boundary_adjust

  subroutine dummy_write_info(fh)
    use mod_global_parameters
    integer, intent(in)                 :: fh !< File handle
    integer, dimension(MPI_STATUS_SIZE) :: st
    integer                             :: er

    ! Number of physics parameters
    integer, parameter                  :: n_par = 0

    call MPI_FILE_WRITE(fh, n_par, 1, MPI_INTEGER, st, er)
  end subroutine dummy_write_info

  subroutine dummy_small_values(s, ixI^L, ixO^L, subname)
    use mod_global_parameters
    type(state), intent(inout)      :: s
    integer, intent(in)             :: ixI^L,ixO^L
    character(len=*), intent(in)    :: subname
  end subroutine dummy_small_values

  subroutine dummy_update_faces(ixI^L,ixO^L,qt,qdt,fC,fE,sCT,s,vcts)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qt,qdt
    type(state)                     :: sCT, s
    type(ct_velocity)               :: vcts
    double precision, intent(in)    :: fC(ixI^S,1:nwflux,1:ndim)
    double precision, intent(inout) :: fE(ixI^S,7-2*ndim:3)
  end subroutine dummy_update_faces

  subroutine dummy_face_to_center(ixO^L,s)
    use mod_global_parameters
    integer, intent(in)                :: ixO^L
    type(state)                        :: s
  end subroutine dummy_face_to_center

  subroutine dummy_evaluate_implicit(qtC,psa)
    use mod_global_parameters
    type(state), target          :: psa(max_blocks)
    double precision, intent(in) :: qtC
    integer :: iigrid, igrid

    ! Just copy in nul state
    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       psa(igrid)%w = 0.0d0*psa(igrid)%w
       if(stagger_grid) psa(igrid)%ws = 0.0d0*psa(igrid)%ws
    end do
    !$OMP END PARALLEL DO

  end subroutine dummy_evaluate_implicit

  subroutine dummy_implicit_update(dtfactor,qdt,qtC,psa,psb)
    use mod_global_parameters
    type(state), target :: psa(max_blocks)   
    type(state), target :: psb(max_blocks)   
    double precision, intent(in) :: qdt      
    double precision, intent(in) :: qtC      
    double precision, intent(in) :: dtfactor 
    integer :: iigrid, igrid

    ! Just copy in psb state when using the scheme without implicit part
    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       psa(igrid)%w = psb(igrid)%w
       if(stagger_grid) psa(igrid)%ws = psb(igrid)%ws
    end do
    !$OMP END PARALLEL DO

  end subroutine dummy_implicit_update


end module mod_physics
