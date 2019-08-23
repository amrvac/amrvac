module mod_physics_roe

  implicit none
  public

  procedure(sub_average), pointer         :: phys_average => null()
  procedure(sub_get_eigenjump), pointer   :: phys_get_eigenjump => null()
  procedure(sub_rtimes), pointer          :: phys_rtimes => null()

  integer :: nworkroe = -1

  abstract interface
     subroutine sub_average(wL, wR, x, ix^L, idim, wroe, workroe)
       use mod_global_parameters
       import
       integer, intent(in)             :: ix^L, idim
       double precision, intent(in)    :: wL(ixG^T, nw), wR(ixG^T, nw)
       double precision, intent(inout) :: wroe(ixG^T, nw)
       double precision, intent(inout) :: workroe(ixG^T, nworkroe)
       double precision, intent(in)    :: x(ixG^T, 1:^ND)
     end subroutine sub_average

     subroutine sub_get_eigenjump(wL, wR, wC, x, ix^L, il, &
          idim, smalla, a, jump, workroe)
       use mod_global_parameters
       import
       integer, intent(in)                          :: ix^L, il, idim
       double precision, dimension(ixG^T, nw)       :: wL, wR, wC
       double precision, dimension(ixG^T)           :: smalla, a, jump
       double precision, dimension(ixG^T, nworkroe) :: workroe
       double precision, intent(in)                 :: x(ixG^T, 1:^ND)
     end subroutine sub_get_eigenjump

     subroutine sub_rtimes(q, w, ix^L, iw, il, idim, rq, workroe)
       use mod_global_parameters
       import
       integer, intent(in)             :: ix^L, iw, il, idim
       double precision, intent(in)    :: w(ixG^T, nw), q(ixG^T)
       double precision, intent(inout) :: rq(ixG^T)
       double precision, intent(inout) :: workroe(ixG^T, nworkroe)
     end subroutine sub_rtimes
  end interface

contains

  subroutine phys_roe_check()
    if (.not. associated(phys_average)) &
         phys_average => dummy_roe_average

    if (.not. associated(phys_get_eigenjump)) &
         phys_get_eigenjump => dummy_roe_get_eigenjump

    if (.not. associated(phys_rtimes)) &
         phys_rtimes => dummy_roe_rtimes

  end subroutine phys_roe_check

  subroutine dummy_roe_average(wL, wR, x, ix^L, idim, wroe, workroe)
    use mod_global_parameters
    integer, intent(in)             :: ix^L, idim
    double precision, intent(in)    :: wL(ixG^T, nw), wR(ixG^T, nw)
    double precision, intent(inout) :: wroe(ixG^T, nw)
    double precision, intent(inout) :: workroe(ixG^T, nworkroe)
    double precision, intent(in)    :: x(ixG^T, 1:^ND)

    call mpistop("dummy_roe_average should not be called")
  end subroutine dummy_roe_average

  subroutine dummy_roe_get_eigenjump(wL, wR, wC, x, ix^L, il, &
       idim, smalla, a, jump, workroe)
    use mod_global_parameters
    integer, intent(in)                          :: ix^L, il, idim
    double precision, dimension(ixG^T, nw)       :: wL, wR, wC
    double precision, dimension(ixG^T)           :: smalla, a, jump
    double precision, dimension(ixG^T, nworkroe) :: workroe
    double precision, intent(in)                 :: x(ixG^T, 1:^ND)

    call mpistop("dummy_roe_get_eigenjump should not be called")
  end subroutine dummy_roe_get_eigenjump

  subroutine dummy_roe_rtimes(q, w, ix^L, iw, il, idim, rq, workroe)
    use mod_global_parameters
    integer, intent(in)             :: ix^L, iw, il, idim
    double precision, intent(in)    :: w(ixG^T, nw), q(ixG^T)
    double precision, intent(inout) :: rq(ixG^T)
    double precision, intent(inout) :: workroe(ixG^T, nworkroe)

    call mpistop("dummy_roe_rtimes should not be called")
  end subroutine dummy_roe_rtimes

end module mod_physics_roe
