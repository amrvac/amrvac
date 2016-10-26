!> This module defines the procedures of a physics module. It contains function
!> pointers for the various supported routines. An actual physics module has to
!> set these pointers to its implementation of these routines.
module mod_physics

  implicit none

  character(len=40) :: physics_type

  procedure(sub_init_params), pointer     :: phys_init_params => null()
  procedure(sub_check_params), pointer    :: phys_check_params => null()
  procedure(sub_read_params), pointer     :: phys_read_params => null()
  procedure(sub_convert), pointer         :: phys_to_conserved => null()
  procedure(sub_convert), pointer         :: phys_to_primitive => null()
  procedure(sub_convert), pointer         :: phys_convert_before_prolong => null()
  procedure(sub_convert), pointer         :: phys_convert_after_prolong => null()
  procedure(sub_convert), pointer         :: phys_convert_before_coarsen => null()
  procedure(sub_convert), pointer         :: phys_convert_after_coarsen => null()
  procedure(sub_average), pointer         :: phys_average => null()
  procedure(sub_get_eigenjump), pointer   :: phys_get_eigenjump => null()
  procedure(sub_rtimes), pointer          :: phys_rtimes => null()
  procedure(sub_get_v), pointer           :: phys_get_v => null()
  procedure(sub_get_cmax), pointer        :: phys_get_cmax => null()
  procedure(sub_get_flux), pointer        :: phys_get_flux => null()
  procedure(sub_get_dt), pointer          :: phys_get_dt => null()
  procedure(sub_add_source_geom), pointer :: phys_add_source_geom => null()
  procedure(sub_add_source), pointer      :: phys_add_source => null()
  procedure(sub_get_aux), pointer         :: phys_get_aux => null()
  procedure(sub_check_w), pointer         :: phys_check_w => null()

  abstract interface

     subroutine sub_init_params()
     end subroutine sub_init_params

     subroutine sub_check_params()
     end subroutine sub_check_params

     subroutine sub_read_params(filename, io_state)
       use mod_global_parameters
       character(len=*), intent(in) :: filename
       integer, intent(out)         :: io_state
     end subroutine sub_read_params

     subroutine sub_convert(ixI^L, ixO^L, w, x)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L, ixO^L
       double precision, intent(inout) :: w(ixI^S, nw)
       double precision, intent(in)    :: x(ixI^S, 1:^ND)
     end subroutine sub_convert

     subroutine sub_average(wL, wR, x, ix^L, idim, wroe, workroe)
       use mod_global_parameters
       integer, intent(in)                          :: ix^L, idim
       double precision, dimension(ixG^T, nw)       :: wL, wR, wroe
       double precision, dimension(ixG^T, nworkroe) :: workroe
       double precision, dimension(ixG^T, 1:^ND)   :: x
     end subroutine sub_average

     subroutine sub_get_eigenjump(wL, wR, wC, x, ix^L, il, &
          idim, smalla, a, jump, workroe)
       use mod_global_parameters
       integer                                      :: ix^L, jx^L, ixC^L, il, idim
       double precision, dimension(ixG^T, nw)       :: wL, wR, wC
       double precision, dimension(ixG^T)           :: smalla, a, jump, v
       double precision, dimension(ixG^T, nworkroe) :: workroe
       double precision, dimension(ixG^T, 1:^ND)   :: x
     end subroutine sub_get_eigenjump

     subroutine sub_rtimes(q, w, ix^L, iw, il, idim, rq, workroe)
       use mod_global_parameters
       integer, intent(in)             :: ix^L, iw, il, idim
       double precision, intent(in)    :: w(ixG^T, nw), q(ixG^T)
       double precision, intent(inout) :: rq(ixG^T)
       double precision, intent(inout) :: workroe(ixG^T, nworkroe)
     end subroutine sub_rtimes

     subroutine sub_get_v(w, x, ixI^L, ixO^L, idim, v)
       use mod_global_parameters
       integer, intent(in)           :: ixI^L, ixO^L, idim
       double precision, intent(in)  :: w(ixI^S, nw), x(ixI^S, 1:^ND)
       double precision, intent(out) :: v(ixG^T)
     end subroutine sub_get_v

     subroutine sub_get_cmax(w, x, ixI^L, ixO^L, idim, cmax, cmin)
       use mod_global_parameters
       integer, intent(in)                       :: ixI^L, ixO^L, idim
       double precision, intent(in)              :: w(ixI^S, nw), x(ixI^S, 1:^ND)
       double precision, intent(inout)           :: cmax(ixG^T)
       double precision, optional, intent(inout) :: cmin(ixG^T)
     end subroutine sub_get_cmax

     subroutine sub_get_flux(w, x, ixI^L, ixO^L, iw, idim, f, transport)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L, ixO^L, iw, idim
       double precision, intent(in)    :: w(ixI^S, 1:nw), x(ixI^S, 1:^ND)
       double precision, intent(inout) :: f(ixG^T)
       logical, intent(out)            :: transport
     end subroutine sub_get_flux

     subroutine sub_add_source_geom(qdt, ixI^L, ixO^L, wCT, w, x)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L, ixO^L
       double precision, intent(in)    :: qdt, x(ixI^S, 1:^ND)
       double precision, intent(inout) :: wCT(ixI^S, 1:nw), w(ixI^S, 1:nw)
     end subroutine sub_add_source_geom

     subroutine sub_add_source(qdt, ixI^L, ixO^L, iw^LIM, qtC, wCT, qt, w, x, qsourcesplit)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L, ixO^L, iw^LIM
       double precision, intent(in)    :: qdt, qtC, qt, x(ixI^S, 1:^ND)
       double precision, intent(inout) :: wCT(ixI^S, 1:nw), w(ixI^S, 1:nw)
       logical, intent(in)             :: qsourcesplit
     end subroutine sub_add_source

     subroutine sub_get_dt(w, ixI^L, ixO^L, dtnew, dx^D, x)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L, ixO^L
       double precision, intent(in)    :: dx^D, x(ixI^S, 1:^ND)
       double precision, intent(inout) :: w(ixI^S, 1:nw), dtnew
     end subroutine sub_get_dt

     subroutine sub_get_aux(clipping,w,x,ixI^L,ixO^L,subname)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L, ixO^L
       double precision, intent(in)    :: x(ixI^S,1:ndim)
       double precision, intent(inout) :: w(ixI^S,nw)
       logical, intent(in)             :: clipping
       character(len=*)                :: subname
     end subroutine sub_get_aux

     subroutine sub_check_w(checkprimitive,ixI^L,ixO^L,w,flag)
       use mod_global_parameters
       logical             :: checkprimitive
       integer, intent(in) :: ixI^L, ixO^L
       double precision    :: w(ixI^S,nw)
       logical             :: flag(ixG^T)
     end subroutine sub_check_w

  end interface

contains

  subroutine phys_check_methods()
    ! Checks whether the required physics methods have been defined

    if (.not. associated(phys_init_params)) &
         phys_init_params => dummy_init_params

    if (.not. associated(phys_check_params)) &
         phys_check_params => dummy_check_params

    if (.not. associated(phys_read_params)) &
         phys_read_params => dummy_read_params

    if (.not. associated(phys_to_conserved)) &
         phys_to_conserved => dummy_convert

    if (.not. associated(phys_to_primitive)) &
         phys_to_primitive => dummy_convert

    if (.not. associated(phys_convert_before_prolong)) &
         phys_to_conserved => dummy_convert

    if (.not. associated(phys_convert_after_prolong)) &
         phys_to_primitive => dummy_convert

    if (.not. associated(phys_convert_before_prolong)) &
         phys_to_conserved => dummy_convert

    if (.not. associated(phys_convert_after_prolong)) &
         phys_to_primitive => dummy_convert

    if (.not. associated(phys_average)) &
         call mpistop("Error: no average method has been specified")

    if (.not. associated(phys_get_eigenjump)) &
         call mpistop("Error: no eigenjump method has been specified")

    if (.not. associated(phys_rtimes)) &
         call mpistop("Error: no rtimes method has been specified")

    if (.not. associated(phys_get_v)) &
         call mpistop("Error: no get_v method has been specified")

    if (.not. associated(phys_get_cmax)) &
         call mpistop("Error: no get_cmax method has been specified")

    if (.not. associated(phys_get_flux)) &
         call mpistop("Error: no get_flux method has been specified")

    if (.not. associated(phys_get_dt)) &
         phys_get_dt => dummy_get_dt

    if (.not. associated(phys_add_source_geom)) &
         phys_add_source_geom => dummy_add_source_geom

    if (.not. associated(phys_add_source)) &
         phys_add_source => dummy_add_source

    if (.not. associated(phys_get_aux)) &
         phys_get_aux => dummy_get_aux

    if (.not. associated(phys_check_w)) &
         phys_check_w => dummy_check_w

  end subroutine phys_check_methods

  subroutine dummy_init_params()
  end subroutine dummy_init_params

  subroutine dummy_check_params()
  end subroutine dummy_check_params

  subroutine dummy_read_params(filename, io_state)
    use mod_global_parameters
    character(len=*), intent(in) :: filename
    integer, intent(out)         :: io_state
  end subroutine dummy_read_params

  subroutine dummy_convert(ixI^L, ixO^L, w, x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:^ND)
  end subroutine dummy_convert

  subroutine dummy_add_source_geom(qdt, ixI^L, ixO^L, wCT, w, x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, x(ixI^S, 1:^ND)
    double precision, intent(inout) :: wCT(ixI^S, 1:nw), w(ixI^S, 1:nw)
  end subroutine dummy_add_source_geom

  subroutine dummy_add_source(qdt, ixI^L, ixO^L, iw^LIM, qtC, wCT, qt, w, x, qsourcesplit)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L, iw^LIM
    double precision, intent(in)    :: qdt, qtC, qt, x(ixI^S, 1:^ND)
    double precision, intent(inout) :: wCT(ixI^S, 1:nw), w(ixI^S, 1:nw)
    logical, intent(in)             :: qsourcesplit
  end subroutine dummy_add_source

  subroutine dummy_get_dt(w, ixI^L, ixO^L, dtnew, dx^D, x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: dx^D, x(ixI^S, 1:^ND)
    double precision, intent(inout) :: w(ixI^S, 1:nw), dtnew

    dtnew = bigdouble
  end subroutine dummy_get_dt

  subroutine dummy_get_aux(clipping,w,x,ixI^L,ixO^L,subname)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,nw)
    logical, intent(in)             :: clipping
    character(len=*)                :: subname
  end subroutine dummy_get_aux

  subroutine dummy_check_w(checkprimitive,ixI^L,ixO^L,w,flag)
    use mod_global_parameters
    logical             :: checkprimitive
    integer, intent(in) :: ixI^L, ixO^L
    double precision    :: w(ixI^S,nw)
    logical             :: flag(ixG^T)
    flag(ixO^S)=.true.
  end subroutine dummy_check_w

end module mod_physics
