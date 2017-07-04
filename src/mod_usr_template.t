!> This is a template for a new user problem
module mod_usr

  ! Include a physics module
  use mod_hd

  implicit none

  ! Custom variables can be defined here
  ! ...

contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()
    use mod_usr_methods

    ! A routine for initial conditions is always required
    usr_init_one_grid => initial_conditions

    ! Specify other user routines, for a list see mod_usr_methods.t
    ! ...

    ! Active the physics module
    call hd_activate()
  end subroutine usr_init

  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixG^L, ix^L, w, x)
    use mod_global_parameters

    integer, intent(in)             :: ixG^L, ix^L
    double precision, intent(in)    :: x(ixG^S, ndim)
    double precision, intent(inout) :: w(ixG^S, nw)

    ! Set initial values for w
    ! w(ix^S, :) = ...

  end subroutine initial_conditions

  ! Extra routines can be placed here
  ! ...

end module mod_usr
