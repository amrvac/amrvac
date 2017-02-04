!> This module defines the procedures of a physics module. It contains function
!> pointers for the various supported routines. An actual physics module has to
!> set these pointers to its implementation of these routines.
module mod_physics

  use mod_physics_hllc
  use mod_physics_roe
  use mod_physics_ppm

  implicit none
  public

  !> String describing the physics type of the simulation
  character(len=40)    :: physics_type

  !> Array per direction per variable, which can be used to specify that certain
  !> fluxes have to be treated differently
  integer, allocatable :: flux_type(:, :)

  !> Indicates a normal flux
  integer, parameter   :: flux_default        = 0
  !> Indicates the flux should be treated with tvdlf
  integer, parameter   :: flux_tvdlf          = 1
  !> Indicates dissipation should be omitted
  integer, parameter   :: flux_no_dissipation = 2

  !> To use wider stencils in flux calculations. A value of 1 will extend it by
  !> one cell in both directions, in any dimension
  integer :: phys_wider_stencil = 0

  procedure(sub_check_params), pointer    :: phys_check_params           => null()
  procedure(sub_convert), pointer         :: phys_to_conserved           => null()
  procedure(sub_convert), pointer         :: phys_to_primitive           => null()
  procedure(sub_convert), pointer         :: phys_convert_before_prolong => null()
  procedure(sub_convert), pointer         :: phys_convert_after_prolong  => null()
  procedure(sub_convert), pointer         :: phys_convert_before_coarsen => null()
  procedure(sub_convert), pointer         :: phys_convert_after_coarsen  => null()
  procedure(sub_modify_wLR), pointer      :: phys_modify_wLR             => null()
  procedure(sub_get_cmax), pointer        :: phys_get_cmax               => null()
  procedure(sub_get_flux), pointer        :: phys_get_flux               => null()
  procedure(sub_get_dt), pointer          :: phys_get_dt                 => null()
  procedure(sub_add_source_geom), pointer :: phys_add_source_geom        => null()
  procedure(sub_add_source), pointer      :: phys_add_source             => null()
  procedure(sub_get_aux), pointer         :: phys_get_aux                => null()
  procedure(sub_check_w), pointer         :: phys_check_w                => null()
  procedure(sub_get_pthermal), pointer    :: phys_get_pthermal           => null()

  abstract interface

     subroutine sub_check_params()
     end subroutine sub_check_params

     subroutine sub_convert(ixI^L, ixO^L, w, x)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L, ixO^L
       double precision, intent(inout) :: w(ixI^S, nw)
       double precision, intent(in)    :: x(ixI^S, 1:^ND)
     end subroutine sub_convert

     subroutine sub_modify_wLR(wLC, wRC, ixI^L, ixO^L, idir)
       use mod_global_parameters, only: nw
       double precision, intent(inout)    :: wLC(ixI^S,1:nw), wRC(ixI^S,1:nw)
       integer, intent(in)                :: ixI^L, ixO^L, idir
     end subroutine sub_modify_wLR

     subroutine sub_get_cmax(w, x, ixI^L, ixO^L, idim, cmax, cmin)
       use mod_global_parameters
       integer, intent(in)                       :: ixI^L, ixO^L, idim
       double precision, intent(in)              :: w(ixI^S, nw), x(ixI^S, 1:^ND)
       double precision, intent(inout)           :: cmax(ixG^T)
       double precision, intent(inout), optional :: cmin(ixG^T)
     end subroutine sub_get_cmax

     subroutine sub_get_flux(w, x, ixI^L, ixO^L, idim, f)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L, ixO^L, idim
       double precision, intent(in)    :: w(ixI^S, 1:nw), x(ixI^S, 1:^ND)
       double precision, intent(out)   :: f(ixG^T, nwflux)
     end subroutine sub_get_flux

     subroutine sub_add_source_geom(qdt, ixI^L, ixO^L, wCT, w, x)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L, ixO^L
       double precision, intent(in)    :: qdt, x(ixI^S, 1:^ND)
       double precision, intent(inout) :: wCT(ixI^S, 1:nw), w(ixI^S, 1:nw)
     end subroutine sub_add_source_geom

     subroutine sub_add_source(qdt, ixI^L, ixO^L, wCT, w, x, qsourcesplit)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L, ixO^L
       double precision, intent(in)    :: qdt
       double precision, intent(in)    :: wCT(ixI^S, 1:nw), x(ixI^S, 1:ndim)
       double precision, intent(inout) :: w(ixI^S, 1:nw)
       logical, intent(in)             :: qsourcesplit
     end subroutine sub_add_source

     subroutine sub_get_dt(w, ixI^L, ixO^L, dtnew, dx^D, x)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L, ixO^L
       double precision, intent(in)    :: dx^D, x(ixI^S, 1:^ND)
       double precision, intent(in)    :: w(ixI^S, 1:nw)
       double precision, intent(inout) :: dtnew
     end subroutine sub_get_dt

     subroutine sub_get_aux(clipping,w,x,ixI^L,ixO^L,subname)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L, ixO^L
       double precision, intent(in)    :: x(ixI^S,1:ndim)
       double precision, intent(inout) :: w(ixI^S,nw)
       logical, intent(in)             :: clipping
       character(len=*)                :: subname
     end subroutine sub_get_aux

     subroutine sub_check_w(primitive, ixI^L, ixO^L, w, w_flag)
       use mod_global_parameters
       logical, intent(in)          :: primitive
       integer, intent(in)          :: ixI^L, ixO^L
       double precision, intent(in) :: w(ixI^S,nw)
       integer, intent(inout)       :: w_flag(ixG^T)
     end subroutine sub_check_w

     subroutine sub_get_pthermal(w,x,ixI^L,ixO^L,pth)
       use mod_global_parameters
       integer, intent(in)          :: ixI^L, ixO^L
       double precision, intent(in) :: w(ixI^S,nw)
       double precision, intent(in) :: x(ixI^S,1:ndim)
       double precision, intent(out):: pth(ixI^S)
     end subroutine sub_get_pthermal

  end interface

contains

  subroutine phys_check()
    use mod_global_parameters, only: nw, ndir

    use mod_physics_hllc, only: phys_hllc_check
    use mod_physics_roe, only: phys_roe_check
    use mod_physics_ppm, only: phys_ppm_check

    ! Check whether custom flux types have been defined
    if (.not. allocated(flux_type)) then
       allocate(flux_type(ndir, nw))
       flux_type = flux_default
    else if (any(shape(flux_type) /= [ndir, nw])) then
       call mpistop("phys_check error: flux_type has wrong shape")
    end if

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

    if (.not. associated(phys_convert_before_prolong)) &
         phys_convert_before_prolong => dummy_convert

    if (.not. associated(phys_convert_after_prolong)) &
         phys_convert_after_prolong => dummy_convert

    if (.not. associated(phys_convert_before_coarsen)) &
         phys_convert_before_coarsen => dummy_convert

    if (.not. associated(phys_convert_after_coarsen)) &
         phys_convert_after_coarsen => dummy_convert

    if (.not. associated(phys_modify_wLR)) &
         phys_modify_wLR => dummy_modify_wLR

    if (.not. associated(phys_get_cmax)) &
         call mpistop("Error: no phys_get_cmax not defined")

    if (.not. associated(phys_get_flux)) &
         call mpistop("Error: no phys_get_flux not defined")

    if (.not. associated(phys_get_dt)) &
         call mpistop("Error: no phys_get_dt not defined")

    if (.not. associated(phys_add_source_geom)) &
         phys_add_source_geom => dummy_add_source_geom

    if (.not. associated(phys_add_source)) &
         phys_add_source => dummy_add_source

    if (.not. associated(phys_get_aux)) &
         phys_get_aux => dummy_get_aux

    if (.not. associated(phys_check_w)) &
         phys_check_w => dummy_check_w

    if (.not. associated(phys_get_pthermal)) &
         phys_get_pthermal => dummy_get_pthermal

  end subroutine phys_check

  subroutine dummy_init_params()
  end subroutine dummy_init_params

  subroutine dummy_check_params()
  end subroutine dummy_check_params

  subroutine dummy_convert(ixI^L, ixO^L, w, x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:^ND)
  end subroutine dummy_convert

  subroutine dummy_modify_wLR(wLC, wRC, ixI^L, ixO^L, idir)
    use mod_global_parameters, only: nw
    double precision, intent(inout)    :: wLC(ixI^S,1:nw), wRC(ixI^S,1:nw)
    integer, intent(in)                :: ixI^L, ixO^L, idir
  end subroutine dummy_modify_wLR

  subroutine dummy_add_source_geom(qdt, ixI^L, ixO^L, wCT, w, x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, x(ixI^S, 1:^ND)
    double precision, intent(inout) :: wCT(ixI^S, 1:nw), w(ixI^S, 1:nw)
  end subroutine dummy_add_source_geom

  subroutine dummy_add_source(qdt, ixI^L, ixO^L, wCT, w, x, qsourcesplit)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixI^S, 1:nw), x(ixI^S, 1:ndim)
    double precision, intent(inout) :: w(ixI^S, 1:nw)
    logical, intent(in)             :: qsourcesplit
  end subroutine dummy_add_source

  subroutine dummy_get_aux(clipping,w,x,ixI^L,ixO^L,subname)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,nw)
    logical, intent(in)             :: clipping
    character(len=*)                :: subname
  end subroutine dummy_get_aux

  subroutine dummy_check_w(primitive, ixI^L, ixO^L, w, w_flag)
    use mod_global_parameters
    logical, intent(in)          :: primitive
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S,nw)
    integer, intent(inout)       :: w_flag(ixG^T)

    w_flag(ixO^S) = 0             ! All okay
  end subroutine dummy_check_w

  subroutine dummy_get_pthermal(w, x, ixI^L, ixO^L, pth)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: pth(ixI^S)

    call mpistop("No get_pthermal method specified")
  end subroutine dummy_get_pthermal

end module mod_physics
