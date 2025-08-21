#:if defined('BFIELD')
module mod_bfield

  use mod_global_parameters
  implicit none
  private
  integer, parameter :: dp=kind(0.0d0)


  public :: magnetic_field, magnetic_field_divergence

contains


 ! analytic formulae for the unit vectors along B
 pure real(dp) function magnetic_field(x, idim) result(field)
    !$acc routine seq
    real(dp), intent(in)    :: x(1:ndim)
    integer, value, intent(in)     :: idim
    ! real(dp)                :: field

    if (idim == 1) field =  0.0_dp
    if (idim == 2) field =  0.0_dp
    if (idim == 3) field = -1.0_dp

  end function magnetic_field

 ! analytic formula for the divergence of the unit vectors along B
 pure real(dp) function magnetic_field_divergence(x) result(field)
    !$acc routine seq
    real(dp), intent(in)    :: x(1:ndim)
    ! real(dp)                :: field

    field =  0.0_dp

  end function magnetic_field_divergence


end module mod_bfield
#:endif
