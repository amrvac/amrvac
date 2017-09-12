module mod_usr
  use mod_rho

  implicit none

  integer, parameter :: dp = kind(0.0d0)
  integer :: i_err
  integer :: i_sol

contains

  subroutine usr_init()
    use mod_usr_methods
    use mod_variables

    usr_init_one_grid => initonegrid_usr
    usr_process_grid => set_error
    usr_print_log => print_error

    call rho_activate()

    i_err = var_set_extravar("error", "error")
    i_sol = var_set_extravar("solution", "solution")
  end subroutine usr_init

  subroutine initonegrid_usr(ixG^L,ix^L,w,x)
    use mod_global_parameters

    integer, intent(in) :: ixG^L, ix^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    w(ix^S,rho_) = solution(x(ix^S, 1), 0.0_dp)

  end subroutine initonegrid_usr

  elemental function solution(x, t) result(val)
    real(dp), intent(in) :: x, t
    real(dp)             :: val
    real(dp), parameter  :: pi = acos(-1.0_dp)

    val = sin(2 * pi * (x - rho_v(1) * t))
  end function solution

  subroutine print_error()
    use mod_global_parameters
    use mod_input_output, only: get_volume_average
    double precision   :: modes(nw, 2), volume

    call get_volume_average(1, modes(:, 1), volume)
    call get_volume_average(2, modes(:, 2), volume)

    if (mype == 0) then
       write(*, "(A,2E16.8)") " time -- RMSE:", global_time, sqrt(modes(i_err, 2))
    end if
  end subroutine print_error

  subroutine set_error(igrid,level,ixI^L,ixO^L,qt,w,x)
    use mod_global_parameters
    integer, intent(in)             :: igrid,level,ixI^L,ixO^L
    double precision, intent(in)    :: qt,x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    w(ixO^S,i_sol) = solution(x(ixO^S, 1), qt)
    w(ixO^S,i_err) = w(ixO^S,rho_) - w(ixO^S,i_sol)
  end subroutine set_error

end module mod_usr

