!> This module defines the procedures of a physics module. It contains function
!> pointers for the various supported routines. An actual physics module has to
!> set these pointers to its implementation of these routines.
module mod_physics
  use mod_global_parameters, only: name_len

  implicit none
  public


  double precision :: phys_gamma=5.d0/3.d0
  !$acc declare copyin(phys_gamma)

  !> String describing the physics type of the simulation
  character(len=name_len) :: physics_type = ""
  !$acc declare create(physics_type)

  !> To use wider stencils in flux calculations. A value of 1 will extend it by
  !> one cell in both directions, in any dimension
  integer :: phys_wider_stencil = 0
  !$acc declare copyin(phys_wider_stencil)

  !> Whether the physics routines require diagonal ghost cells, for example for
  !> computing a curl.
  logical :: phys_req_diagonal = .true.
  !$acc declare copyin(phys_req_diagonal)

  !> Solve energy equation or not
  logical :: phys_energy=.false.
  !$acc declare copyin(phys_energy)
  
  !> Solve total energy equation or not
  logical :: phys_total_energy=.false.
  !$acc declare copyin(phys_total_energy)

  !> Solve internal energy instead of total energy
  logical :: phys_internal_e=.false.
  !$acc declare copyin(phys_internal_e)

  !> Solve partially ionized one-fluid plasma
  logical :: phys_partial_ionization=.false.
  !$acc declare copyin(phys_partial_ionization)

  !> if equilibrium pressure is splitted
  logical :: phys_equi_pe=.false.
  !$acc declare copyin(phys_equi_pe)

  
  procedure(sub_check_params), pointer    :: phys_check_params           => null()
  procedure(sub_convert), pointer         :: phys_to_conserved           => null()
  procedure(sub_convert), pointer         :: phys_to_primitive           => null()
  procedure(sub_get_cmax), pointer        :: phys_get_cmax               => null()
  procedure(sub_boundary_adjust), pointer :: phys_boundary_adjust        => null()
  procedure(sub_global_source), pointer   :: phys_global_source_after    => null()
  procedure(sub_implicit_update), pointer :: phys_implicit_update        => null()
  procedure(sub_evaluate_implicit),pointer:: phys_evaluate_implicit      => null()
  procedure(sub_clean_divb), pointer      :: phys_clean_divb             => null()
  procedure(sub_special_advance), pointer :: phys_special_advance        => null()
  ! set the equilibrium variables
  procedure(sub_set_equi_vars), pointer   :: phys_set_equi_vars          => null()
  ! subroutine with no parameters which creates EUV images
  procedure(sub_check_params), pointer    :: phys_te_images              => null()
  procedure(sub_small_values), pointer    :: phys_handle_small_values    => null()
  procedure(sub_face_to_center), pointer  :: phys_face_to_center         => null()
  procedure(sub_write_info), pointer      :: phys_write_info             => null()

  interface
     
     subroutine sub_check_params
     end subroutine sub_check_params

     subroutine sub_convert(ixI^L, ixO^L, w, x)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L, ixO^L
       double precision, intent(inout) :: w(ixI^S, nw)
       double precision, intent(in)    :: x(ixI^S, 1:^ND)
     end subroutine sub_convert

     subroutine sub_get_cmax(w, x, ixI^L, ixO^L, idim, cmax)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L, ixO^L, idim
       double precision, intent(in)    :: w(ixI^S, nw), x(ixI^S, 1:^ND)
       double precision, intent(inout) :: cmax(ixI^S)
     end subroutine sub_get_cmax

     subroutine sub_boundary_adjust(igrid,psb)
       use mod_global_parameters
       integer, intent(in) :: igrid
       type(state), target :: psb(max_blocks)
     end subroutine sub_boundary_adjust
     
     subroutine sub_face_to_center(ixO^L,s)
       use mod_global_parameters
       integer, intent(in)                :: ixO^L
       type(state)                        :: s
     end subroutine sub_face_to_center

     subroutine sub_write_info(file_handle)
       integer, intent(in) :: file_handle
     end subroutine sub_write_info

     !> Add special advance in each advect step
     subroutine sub_special_advance(qt, psa)
       use mod_global_parameters
       double precision, intent(in) :: qt     !< Current time
       type(state), target :: psa(max_blocks) !< Compute based on this state
     end subroutine sub_special_advance

     !> Add global source terms on complete domain (potentially implicit)
     subroutine sub_global_source(qdt, qt, active)
       use mod_global_parameters
       double precision, intent(in) :: qdt    !< Current time step
       double precision, intent(in) :: qt     !< Current time
       logical, intent(inout)       :: active !< Output if the source is active
     end subroutine sub_global_source

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
     
     !> set equilibrium variables other than b0 (e.g. p0 and rho0)
     subroutine sub_set_equi_vars(igrid)
       integer, intent(in) :: igrid
     end subroutine sub_set_equi_vars

     !> clean initial divb
     subroutine sub_clean_divb(qdt, qt, active)
       use mod_global_parameters
       double precision, intent(in) :: qdt    !< Current time step
       double precision, intent(in) :: qt     !< Current time
       logical, intent(inout)       :: active !< Output if the source is active
     end subroutine sub_clean_divb

     subroutine sub_small_values(primitive, w, x, ixI^L, ixO^L, subname)
       use mod_global_parameters
       logical, intent(in)             :: primitive
       integer, intent(in)             :: ixI^L,ixO^L
       double precision, intent(inout) :: w(ixI^S,1:nw)
       double precision, intent(in)    :: x(ixI^S,1:ndim)
       character(len=*), intent(in)    :: subname
     end subroutine sub_small_values
     
  end interface

  contains


  subroutine phys_check()

    use mod_comm_lib, only: mpistop

    if (physics_type == "") call mpistop("Error: no physics module loaded")

    ! Checks whether the required physics methods have been defined
    if (.not. associated(phys_check_params)) &
         phys_check_params => dummy_check_params
    
    ! Checks whether the required physics methods have been defined
    if (.not. associated(phys_to_conserved)) &
         call mpistop("Error: phys_to_conserved not defined")

    if (.not. associated(phys_to_primitive)) &
         call mpistop("Error: phys_to_primitive not defined")

    if (.not. associated(phys_get_cmax)) &
         call mpistop("Error: no phys_get_cmax not defined")

    if (.not. associated(phys_boundary_adjust)) &
         phys_boundary_adjust => dummy_boundary_adjust

    if (.not. associated(phys_evaluate_implicit)) &
         phys_evaluate_implicit => dummy_evaluate_implicit

    if (.not. associated(phys_implicit_update)) &
         phys_implicit_update => dummy_implicit_update

    if (.not. associated(phys_handle_small_values)) &
         phys_handle_small_values => dummy_small_values

    if (.not. associated(phys_face_to_center)) &
         phys_face_to_center => dummy_face_to_center

    if (.not. associated(phys_write_info)) &
         phys_write_info => dummy_write_info

  end subroutine phys_check
    
  subroutine dummy_boundary_adjust(igrid,psb)
    use mod_global_parameters
    integer, intent(in) :: igrid
    type(state), target :: psb(max_blocks)
  end subroutine dummy_boundary_adjust
  
  subroutine dummy_evaluate_implicit(qtC,psa)
    use mod_global_parameters
    type(state), target :: psa(max_blocks)
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

  subroutine dummy_small_values(primitive, w, x, ixI^L, ixO^L, subname)
    use mod_global_parameters
    logical, intent(in)             :: primitive
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    character(len=*), intent(in)    :: subname
  end subroutine dummy_small_values

  subroutine dummy_check_params
  end subroutine dummy_check_params

  subroutine dummy_face_to_center(ixO^L,s)
    use mod_global_parameters
    integer, intent(in)                :: ixO^L
    type(state)                        :: s
  end subroutine dummy_face_to_center

  subroutine dummy_write_info(fh)
    use mod_global_parameters
    integer, intent(in)                 :: fh !< File handle
    integer, dimension(MPI_STATUS_SIZE) :: st
    integer                             :: er

    ! Number of physics parameters
    integer, parameter                  :: n_par = 0

    call MPI_FILE_WRITE(fh, n_par, 1, MPI_INTEGER, st, er)
  end subroutine dummy_write_info
  
end module mod_physics
