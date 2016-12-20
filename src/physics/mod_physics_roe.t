module mod_physics_roe

  implicit none
  public

  procedure(sub_average), pointer         :: phys_average => null()
  procedure(sub_get_eigenjump), pointer   :: phys_get_eigenjump => null()
  procedure(sub_rtimes), pointer          :: phys_rtimes => null()

  integer :: soundRW_ = -1
  integer :: soundLW_ = -1
  integer :: entropW_ = -1
  integer :: shearW0_ = -1
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
       double precision, dimension(ixG^T)           :: smalla, a, jump, v
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

  subroutine phys_roe_check_methods()
    if (.not. associated(phys_average)) &
         call mpistop("Error: no average method has been specified")

    if (.not. associated(phys_get_eigenjump)) &
         call mpistop("Error: no eigenjump method has been specified")

    if (.not. associated(phys_rtimes)) &
         call mpistop("Error: no rtimes method has been specified")
  end subroutine phys_roe_check_methods

end module mod_physics_roe
