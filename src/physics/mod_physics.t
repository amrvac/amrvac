!> This module defines the procedures of a physics module. It contains function
!> pointers for the various supported routines. An actual physics module has to
!> set these pointers to its implementation of these routines.
module mod_physics
  use mod_global_parameters, only: name_len, max_nw
  use mod_physics_hllc
  use mod_physics_roe

  implicit none
  public


  double precision :: phys_gamma=5.d0/3.d0

  !> To use wider stencils in flux calculations. A value of 1 will extend it by
  !> one cell in both directions, in any dimension
  integer :: phys_wider_stencil = 0

  !> Array per direction per variable, which can be used to specify that certain
  !> fluxes have to be treated differently
  integer, allocatable :: flux_type(:, :)

  !> Indicates a normal flux
  integer, parameter   :: flux_default        = 0
  !> Indicates the flux should be treated with tvdlf
  integer, parameter   :: flux_tvdlf          = 1
  !> Indicates dissipation should be omitted
  integer, parameter   :: flux_no_dissipation = 2
  !> Indicates the flux should be specially treated
  integer, parameter   :: flux_special        = 3
  !> Indicates the flux should be treated with hll
  integer, parameter   :: flux_hll        = 4

  !> Solve energy equation or not
  logical :: phys_energy=.false.

  !> Solve total energy equation or not
  logical :: phys_total_energy=.false.

  !> Solve internal energy instead of total energy
  logical :: phys_internal_e=.false.

  !> Solve partially ionized one-fluid plasma
  logical :: phys_partial_ionization=.false.

  !> if equilibrium pressure is splitted
  logical :: phys_equi_pe=.false.

  !> String describing the physics type of the simulation
  character(len=name_len) :: physics_type = ""

  procedure(sub_check_params), pointer    :: phys_check_params           => null()
  procedure(sub_set_mg_bounds), pointer   :: phys_set_mg_bounds          => null()
  procedure(sub_convert), pointer         :: phys_to_conserved           => null()
  procedure(sub_convert), pointer         :: phys_to_primitive           => null()
  procedure(sub_modify_wLR), pointer      :: phys_modify_wLR             => null()
  procedure(sub_get_cmax), pointer        :: phys_get_cmax               => null()
  procedure(sub_get_a2max), pointer       :: phys_get_a2max              => null()
  procedure(sub_get_cs2max), pointer      :: phys_get_cs2max             => null()
  procedure(sub_get_tcutoff), pointer     :: phys_get_tcutoff            => null()
  procedure(sub_trac_after_setdt), pointer:: phys_trac_after_setdt       => null()
  procedure(sub_get_H_speed), pointer     :: phys_get_H_speed            => null()
  procedure(sub_get_cbounds), pointer     :: phys_get_cbounds            => null()
  procedure(sub_get_flux), pointer        :: phys_get_flux               => null()
  procedure(sub_get_v), pointer           :: phys_get_v                  => null()
  procedure(sub_get_rho), pointer         :: phys_get_rho                => null()
  procedure(sub_get_dt), pointer          :: phys_get_dt                 => null()
  procedure(sub_add_source_geom), pointer :: phys_add_source_geom        => null()
  procedure(sub_add_source), pointer      :: phys_add_source             => null()
  procedure(sub_global_source), pointer   :: phys_global_source_after    => null()
  procedure(sub_special_advance), pointer :: phys_special_advance        => null()
  procedure(sub_check_w), pointer         :: phys_check_w                => null()
  procedure(sub_get_pthermal), pointer    :: phys_get_pthermal           => null()
  procedure(sub_get_tgas), pointer        :: phys_get_tgas               => null()
  procedure(sub_get_trad), pointer        :: phys_get_trad               => null()
  procedure(sub_boundary_adjust), pointer :: phys_boundary_adjust        => null()
  procedure(sub_write_info), pointer      :: phys_write_info             => null()
  procedure(sub_small_values), pointer    :: phys_handle_small_values    => null()
  procedure(sub_get_ct_velocity), pointer :: phys_get_ct_velocity        => null()
  procedure(sub_update_faces), pointer    :: phys_update_faces           => null()
  procedure(sub_face_to_center), pointer  :: phys_face_to_center         => null()
  procedure(sub_implicit_update), pointer :: phys_implicit_update        => null()
  procedure(sub_evaluate_implicit),pointer:: phys_evaluate_implicit      => null()
  procedure(sub_clean_divb), pointer      :: phys_clean_divb             => null()
  ! set the equilibrium variables
  procedure(sub_set_equi_vars), pointer   :: phys_set_equi_vars          => null()
  ! subroutine with no parameters which creates EUV images
  procedure(sub_check_params), pointer    :: phys_te_images              => null()
  ! to update temperature variable with partial ionization
  procedure(sub_update_temperature), pointer :: phys_update_temperature  => null()
  procedure(sub_get_auxiliary), pointer         :: phys_get_auxiliary         => null()
  procedure(sub_get_auxiliary_prim), pointer    :: phys_get_auxiliary_prim    => null()

  abstract interface

     subroutine sub_check_params
     end subroutine sub_check_params

     subroutine sub_set_mg_bounds
       use mod_global_parameters
       use mod_usr_methods
     end subroutine sub_set_mg_bounds

     subroutine sub_boundary_adjust(igrid,psb)
       use mod_global_parameters
       integer, intent(in) :: igrid
       type(state), target :: psb(max_blocks)
     end subroutine sub_boundary_adjust

     subroutine sub_convert(ixI^L, ixO^L, w, x)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L, ixO^L
       double precision, intent(inout) :: w(ixI^S, nw)
       double precision, intent(in)    :: x(ixI^S, 1:^ND)
     end subroutine sub_convert

     subroutine sub_modify_wLR(ixI^L, ixO^L, qt, wLC, wRC, wLp, wRp, s, idir)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L, ixO^L, idir
       double precision, intent(in)    :: qt
       double precision, intent(inout) :: wLC(ixI^S,1:nw), wRC(ixI^S,1:nw)
       double precision, intent(inout) :: wLp(ixI^S,1:nw), wRp(ixI^S,1:nw)
       type(state)                     :: s
     end subroutine sub_modify_wLR

     subroutine sub_get_cmax(w, x, ixI^L, ixO^L, idim, cmax)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L, ixO^L, idim
       double precision, intent(in)    :: w(ixI^S, nw), x(ixI^S, 1:^ND)
       double precision, intent(inout) :: cmax(ixI^S)
     end subroutine sub_get_cmax

     subroutine sub_get_a2max(w, x, ixI^L, ixO^L, a2max)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L, ixO^L
       double precision, intent(in)    :: w(ixI^S, nw), x(ixI^S, 1:^ND)
       double precision, intent(inout) :: a2max(ndim)
     end subroutine sub_get_a2max

     subroutine sub_get_cs2max(w, x, ixI^L, ixO^L, cs2max)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L, ixO^L
       double precision, intent(in)    :: w(ixI^S, nw), x(ixI^S, 1:^ND)
       double precision, intent(inout) :: cs2max
     end subroutine sub_get_cs2max

     subroutine sub_get_tcutoff(ixI^L,ixO^L,w,x,tco_local,Tmax_local)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L, ixO^L
       double precision, intent(inout)    :: w(ixI^S, nw)
       double precision, intent(in)    :: x(ixI^S, 1:^ND)
       double precision, intent(out) :: tco_local, Tmax_local
     end subroutine sub_get_tcutoff

     subroutine sub_trac_after_setdt(trac_alfa,tco,T_peak,T_bott)
       double precision, intent(in)    :: trac_alfa,tco,T_peak,T_bott
     end subroutine sub_trac_after_setdt

     subroutine sub_get_v(w,x,ixI^L,ixO^L,v)
       use mod_global_parameters
       integer, intent(in)           :: ixI^L, ixO^L
       double precision, intent(in)  :: w(ixI^S,nw), x(ixI^S,1:^ND)
       double precision, intent(out) :: v(ixI^S,1:ndir)
     end subroutine sub_get_v

     subroutine sub_get_rho(w,x,ixI^L,ixO^L,rho)
       use mod_global_parameters
       integer, intent(in)           :: ixI^L, ixO^L
       double precision, intent(in)  :: w(ixI^S,nw), x(ixI^S,1:^ND)
       double precision, intent(out) :: rho(ixI^S)
     end subroutine sub_get_rho

     subroutine sub_get_H_speed(wprim,x,ixI^L,ixO^L,idim,Hspeed)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L, ixO^L, idim
       double precision, intent(in)    :: wprim(ixI^S, nw)
       double precision, intent(in)    :: x(ixI^S,1:ndim)
       double precision, intent(out)   :: Hspeed(ixI^S,1:number_species)
     end subroutine sub_get_H_speed

     subroutine sub_get_cbounds(wLC, wRC, wLp, wRp, x, ixI^L, ixO^L, idim, Hspeed, cmax, cmin)
       use mod_global_parameters
       use mod_variables
       integer, intent(in)             :: ixI^L, ixO^L, idim
       double precision, intent(in)    :: wLC(ixI^S, nw), wRC(ixI^S, nw)
       double precision, intent(in)    :: wLp(ixI^S, nw), wRp(ixI^S, nw)
       double precision, intent(in)    :: x(ixI^S, 1:^ND)
       double precision, intent(inout) :: cmax(ixI^S,1:number_species)
       double precision, intent(inout), optional :: cmin(ixI^S,1:number_species)
       double precision, intent(in)    :: Hspeed(ixI^S,1:number_species)
     end subroutine sub_get_cbounds

     subroutine sub_get_flux(wC, w, x, ixI^L, ixO^L, idim, f)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L, ixO^L, idim
       double precision, intent(in)    :: wC(ixI^S, 1:nw)
       double precision, intent(in)    :: w(ixI^S, 1:nw)
       double precision, intent(in)    :: x(ixI^S, 1:^ND)
       double precision, intent(out)   :: f(ixI^S, nwflux)
     end subroutine sub_get_flux

     subroutine sub_add_source_geom(qdt, dtfactor, ixI^L, ixO^L, wCT, wprim, w, x)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L, ixO^L
       double precision, intent(in)    :: qdt, dtfactor, x(ixI^S, 1:^ND)
       double precision, intent(inout) :: wCT(ixI^S, 1:nw), wprim(ixI^S, 1:nw), w(ixI^S, 1:nw)
     end subroutine sub_add_source_geom

     subroutine sub_add_source(qdt, dtfactor, ixI^L, ixO^L, wCT, wCTprim, w, x, &
          qsourcesplit, active)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L, ixO^L
       double precision, intent(in)    :: qdt, dtfactor
       double precision, intent(in)    :: wCT(ixI^S, 1:nw), wCTprim(ixI^S,1:nw), x(ixI^S, 1:ndim)
       double precision, intent(inout) :: w(ixI^S, 1:nw)
       logical, intent(in)             :: qsourcesplit
       logical, intent(inout)          :: active !< Needs to be set to true when active
     end subroutine sub_add_source

     !> Add global source terms on complete domain (potentially implicit)
     subroutine sub_global_source(qdt, qt, active)
       use mod_global_parameters
       double precision, intent(in) :: qdt    !< Current time step
       double precision, intent(in) :: qt     !< Current time
       logical, intent(inout)       :: active !< Output if the source is active
     end subroutine sub_global_source

     !> clean initial divb
     subroutine sub_clean_divb(qdt, qt, active)
       use mod_global_parameters
       double precision, intent(in) :: qdt    !< Current time step
       double precision, intent(in) :: qt     !< Current time
       logical, intent(inout)       :: active !< Output if the source is active
     end subroutine sub_clean_divb

     !> set equilibrium variables other than b0 (e.g. p0 and rho0)
     subroutine sub_set_equi_vars(igrid)
       integer, intent(in) :: igrid
     end subroutine sub_set_equi_vars


     !> Add special advance in each advect step
     subroutine sub_special_advance(qt, psa)
       use mod_global_parameters
       double precision, intent(in) :: qt     !< Current time
       type(state), target :: psa(max_blocks) !< Compute based on this state
     end subroutine sub_special_advance

     subroutine sub_get_dt(w, ixI^L, ixO^L, dtnew, dx^D, x)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L, ixO^L
       double precision, intent(in)    :: dx^D, x(ixI^S, 1:^ND)
       double precision, intent(in)    :: w(ixI^S, 1:nw)
       double precision, intent(inout) :: dtnew
     end subroutine sub_get_dt

     subroutine sub_check_w(primitive, ixI^L, ixO^L, w, w_flag)
       use mod_global_parameters
       logical, intent(in)          :: primitive
       integer, intent(in)          :: ixI^L, ixO^L
       double precision, intent(in) :: w(ixI^S,1:nw)
       logical, intent(inout)       :: w_flag(ixI^S,1:nw)
     end subroutine sub_check_w

     subroutine sub_get_pthermal(w,x,ixI^L,ixO^L,pth)
       use mod_global_parameters
       integer, intent(in)          :: ixI^L, ixO^L
       double precision, intent(in) :: w(ixI^S,nw)
       double precision, intent(in) :: x(ixI^S,1:ndim)
       double precision, intent(out):: pth(ixI^S)
     end subroutine sub_get_pthermal

     subroutine sub_get_auxiliary(ixI^L,ixO^L,w,x)
       use mod_global_parameters
       integer, intent(in)          :: ixI^L, ixO^L
       double precision, intent(inout) :: w(ixI^S,nw)
       double precision, intent(in) :: x(ixI^S,1:ndim)
     end subroutine sub_get_auxiliary

     subroutine sub_get_auxiliary_prim(ixI^L,ixO^L,w)
       use mod_global_parameters
       integer, intent(in)          :: ixI^L, ixO^L
       double precision, intent(inout) :: w(ixI^S,nw)
     end subroutine sub_get_auxiliary_prim

     subroutine sub_get_tgas(w,x,ixI^L,ixO^L,tgas)
       use mod_global_parameters
       integer, intent(in)          :: ixI^L, ixO^L
       double precision, intent(in) :: w(ixI^S,nw)
       double precision, intent(in) :: x(ixI^S,1:ndim)
       double precision, intent(out):: tgas(ixI^S)
     end subroutine sub_get_tgas

     subroutine sub_get_trad(w,x,ixI^L,ixO^L,trad)
       use mod_global_parameters
       integer, intent(in)          :: ixI^L, ixO^L
       double precision, intent(in) :: w(ixI^S,nw)
       double precision, intent(in) :: x(ixI^S,1:ndim)
       double precision, intent(out):: trad(ixI^S)
     end subroutine sub_get_trad

     subroutine sub_write_info(file_handle)
       integer, intent(in) :: file_handle
     end subroutine sub_write_info

     subroutine sub_small_values(primitive, w, x, ixI^L, ixO^L, subname)
       use mod_global_parameters
       logical, intent(in)             :: primitive
       integer, intent(in)             :: ixI^L,ixO^L
       double precision, intent(inout) :: w(ixI^S,1:nw)
       double precision, intent(in)    :: x(ixI^S,1:ndim)
       character(len=*), intent(in)    :: subname
     end subroutine sub_small_values

     subroutine sub_get_var(ixI^L, ixO^L, w, out)
       use mod_global_parameters
       integer, intent(in)           :: ixI^L, ixO^L
       double precision, intent(in)  :: w(ixI^S, nw)
       double precision, intent(out) :: out(ixO^S)
     end subroutine sub_get_var

     subroutine sub_get_ct_velocity(vcts,wLp,wRp,ixI^L,ixO^L,idim,cmax,cmin)
       use mod_global_parameters

       integer, intent(in)             :: ixI^L, ixO^L, idim
       double precision, intent(in)    :: wLp(ixI^S, nw), wRp(ixI^S, nw)
       double precision, intent(in)    :: cmax(ixI^S)
       double precision, intent(in), optional :: cmin(ixI^S)
       type(ct_velocity), intent(inout):: vcts
     end subroutine sub_get_ct_velocity

     subroutine sub_update_faces(ixI^L,ixO^L,qt,qdt,wprim,fC,fE,sCT,s,vcts)
       use mod_global_parameters
       integer, intent(in)                :: ixI^L, ixO^L
       double precision, intent(in)       :: qt, qdt
       ! cell-center primitive variables
       double precision, intent(in)       :: wprim(ixI^S,1:nw)
       ! velocity structure
       type(state)                        :: sCT, s
       type(ct_velocity)                  :: vcts
       double precision, intent(in)       :: fC(ixI^S,1:nwflux,1:ndim)
       double precision, intent(inout)    :: fE(ixI^S,sdim:3)
     end subroutine sub_update_faces

     subroutine sub_face_to_center(ixO^L,s)
       use mod_global_parameters
       integer, intent(in)                :: ixO^L
       type(state)                        :: s
     end subroutine sub_face_to_center

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

     subroutine sub_update_temperature(ixI^L,ixO^L,wCT,w,x)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L, ixO^L
       double precision, intent(in)    :: wCT(ixI^S,nw),x(ixI^S,1:ndim)
       double precision, intent(inout) :: w(ixI^S,nw)
     end subroutine sub_update_temperature

   end interface

contains

  subroutine phys_check()
    use mod_global_parameters, only: nw, ndir

    use mod_physics_hllc, only: phys_hllc_check
    use mod_physics_roe, only: phys_roe_check
    use mod_comm_lib, only: mpistop

    if (physics_type == "") call mpistop("Error: no physics module loaded")

    call phys_hllc_check()
    call phys_roe_check()

    ! Checks whether the required physics methods have been defined
    if (.not. associated(phys_check_params)) &
         phys_check_params => dummy_check_params

    if (.not. associated(phys_to_conserved)) &
         call mpistop("Error: phys_to_conserved not defined")

    if (.not. associated(phys_to_primitive)) &
         call mpistop("Error: phys_to_primitive not defined")

    if (.not. associated(phys_modify_wLR)) &
         phys_modify_wLR => dummy_modify_wLR

    if (.not. associated(phys_get_cmax)) &
         call mpistop("Error: no phys_get_cmax not defined")

    if (.not. associated(phys_get_a2max)) &
         phys_get_a2max => dummy_get_a2max

    if (.not. associated(phys_get_cs2max)) &
         phys_get_cs2max => dummy_get_cs2max

    if (.not. associated(phys_get_H_speed)) &
         phys_get_H_speed => dummy_get_H_speed

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

    if (.not. associated(phys_check_w)) &
         phys_check_w => dummy_check_w

    if (.not. associated(phys_get_pthermal)) &
         phys_get_pthermal => dummy_get_pthermal
    if (.not. associated(phys_get_auxiliary)) &
         phys_get_auxiliary => dummy_get_auxiliary
    if (.not. associated(phys_get_auxiliary_prim)) &
         phys_get_auxiliary_prim => dummy_get_auxiliary_prim

    if (.not. associated(phys_boundary_adjust)) &
         phys_boundary_adjust => dummy_boundary_adjust

    if (.not. associated(phys_write_info)) &
         phys_write_info => dummy_write_info

    if (.not. associated(phys_handle_small_values)) &
         phys_handle_small_values => dummy_small_values

    if (.not. associated(phys_get_ct_velocity)) &
         phys_get_ct_velocity => dummy_get_ct_velocity

    if (.not. associated(phys_update_faces)) &
         phys_update_faces => dummy_update_faces

    if (.not. associated(phys_face_to_center)) &
         phys_face_to_center => dummy_face_to_center

    if (.not. associated(phys_evaluate_implicit)) &
         phys_evaluate_implicit => dummy_evaluate_implicit

    if (.not. associated(phys_implicit_update)) &
         phys_implicit_update => dummy_implicit_update

  end subroutine phys_check

  subroutine dummy_init_params
  end subroutine dummy_init_params

  subroutine dummy_check_params
  end subroutine dummy_check_params

  subroutine dummy_modify_wLR(ixI^L, ixO^L, qt, wLC, wRC, wLp, wRp, s, idir)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L, idir
    double precision, intent(in)    :: qt
    double precision, intent(inout) :: wLC(ixI^S,1:nw), wRC(ixI^S,1:nw)
    double precision, intent(inout) :: wLp(ixI^S,1:nw), wRp(ixI^S,1:nw)
    type(state)                     :: s
  end subroutine dummy_modify_wLR

  subroutine dummy_get_H_speed(wprim,x,ixI^L,ixO^L,idim,Hspeed)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L, idim
    double precision, intent(in)    :: wprim(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out)   :: Hspeed(ixI^S,1:number_species)
  end subroutine dummy_get_H_speed

  subroutine dummy_get_a2max(w, x, ixI^L, ixO^L, a2max)
       use mod_global_parameters
       use mod_comm_lib, only: mpistop
       integer, intent(in)             :: ixI^L, ixO^L
       double precision, intent(in)    :: w(ixI^S, nw), x(ixI^S, 1:^ND)
       double precision, intent(inout) :: a2max(ndim)
       call mpistop("Error: entered dummy_get_a2max")
  end subroutine dummy_get_a2max

  subroutine dummy_get_cs2max(w, x, ixI^L, ixO^L, cs2max)
       use mod_global_parameters
       use mod_comm_lib, only: mpistop
       integer, intent(in)             :: ixI^L, ixO^L
       double precision, intent(in)    :: w(ixI^S, nw), x(ixI^S, 1:^ND)
       double precision, intent(inout) :: cs2max
       call mpistop("Error: entered dummy_get_cs2max")
  end subroutine dummy_get_cs2max

  subroutine dummy_add_source_geom(qdt, dtfactor, ixI^L, ixO^L, wCT, wprim, w, x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, dtfactor, x(ixI^S, 1:^ND)
    double precision, intent(inout) :: wCT(ixI^S, 1:nw), wprim(ixI^S,1:nw),w(ixI^S, 1:nw)
  end subroutine dummy_add_source_geom

  subroutine dummy_add_source(qdt, dtfactor, ixI^L, ixO^L, wCT, wCTprim, w, x, &
       qsourcesplit, active)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt,dtfactor
    double precision, intent(in)    :: wCT(ixI^S, 1:nw), wCTprim(ixI^S,1:nw), x(ixI^S, 1:ndim)
    double precision, intent(inout) :: w(ixI^S, 1:nw)
    logical, intent(in)             :: qsourcesplit
    logical, intent(inout)          :: active
    ! Don't have to set active, since it starts as .false.
  end subroutine dummy_add_source

  subroutine dummy_check_w(primitive, ixI^L, ixO^L, w, w_flag)
    use mod_global_parameters
    logical, intent(in)          :: primitive
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S,1:nw)
    logical, intent(inout)       :: w_flag(ixI^S,1:nw)

    w_flag=.false.             ! All okay
  end subroutine dummy_check_w

  subroutine dummy_get_pthermal(w, x, ixI^L, ixO^L, pth)
    use mod_global_parameters
    use mod_comm_lib, only: mpistop

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: pth(ixI^S)

    call mpistop("No get_pthermal method specified")
  end subroutine dummy_get_pthermal

  subroutine dummy_get_auxiliary(ixI^L, ixO^L, w, x)
    use mod_global_parameters
    use mod_comm_lib, only: mpistop

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)

    !call mpistop("No get_auxiliary method specified")
  end subroutine dummy_get_auxiliary

  subroutine dummy_get_auxiliary_prim(ixI^L, ixO^L, w)
    use mod_global_parameters
    use mod_comm_lib, only: mpistop

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)

    call mpistop("No get_auxiliary_prim method specified")
  end subroutine dummy_get_auxiliary_prim

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

  subroutine dummy_small_values(primitive, w, x, ixI^L, ixO^L, subname)
    use mod_global_parameters
    logical, intent(in)             :: primitive
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    character(len=*), intent(in)    :: subname
  end subroutine dummy_small_values

  subroutine dummy_get_ct_velocity(vcts,wLp,wRp,ixI^L,ixO^L,idim,cmax,cmin)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L, idim
    double precision, intent(in)    :: wLp(ixI^S, nw), wRp(ixI^S, nw)
    double precision, intent(in)    :: cmax(ixI^S)
    double precision, intent(in), optional :: cmin(ixI^S)
    type(ct_velocity), intent(inout):: vcts
  end subroutine dummy_get_ct_velocity

  subroutine dummy_update_faces(ixI^L,ixO^L,qt,qdt,wprim,fC,fE,sCT,s,vcts)
    use mod_global_parameters
    integer, intent(in)                :: ixI^L, ixO^L
    double precision, intent(in)       :: qt, qdt
    ! cell-center primitive variables
    double precision, intent(in)       :: wprim(ixI^S,1:nw)
    type(state)                        :: sCT, s
    type(ct_velocity)                  :: vcts
    double precision, intent(in)       :: fC(ixI^S,1:nwflux,1:ndim)
    double precision, intent(inout)    :: fE(ixI^S,sdim:3)
  end subroutine dummy_update_faces

  subroutine dummy_face_to_center(ixO^L,s)
    use mod_global_parameters
    integer, intent(in)                :: ixO^L
    type(state)                        :: s
  end subroutine dummy_face_to_center

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

end module mod_physics
