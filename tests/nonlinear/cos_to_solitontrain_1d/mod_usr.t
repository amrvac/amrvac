! Initial condition for Zabusky and Kruskal PRL 1965 paper on solitons
module mod_usr
  use mod_nonlinear

  implicit none

  integer :: drho_

contains

  subroutine usr_init()
    use mod_variables

    usr_init_one_grid => initonegrid_usr
    usr_process_grid => set_gradrho
    usr_print_log => print_error

    call nonlinear_activate()

    drho_ = var_set_extravar("drho", "drho")

  end subroutine usr_init

  subroutine initonegrid_usr(ixG^L,ix^L,w,x)
    integer, intent(in) :: ixG^L, ix^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    w(ix^S,rho_) =dcos(dpi*x(ix^S,1))

  end subroutine initonegrid_usr

  subroutine print_error()
    use mod_input_output, only: get_volume_average
    double precision   :: modes(nw,3), volume

    call get_volume_average(1, modes(:, 1), volume)
    call get_volume_average(2, modes(:, 2), volume)
    call get_volume_average(3, modes(:, 3), volume)

    if (mype == 0) then
       ! check on conservation of mass \int \rho dx 
       !                          momentum \int \rho^2/2 dx 
       !                      and energy \int [\rho^3/3+(d\rho/dx)^2 ]dx
       write(*, "(A,4E14.6)") " t mass momentum energy:", global_time, &
            modes(rho_, 1), modes(rho_, 2)/2.0d0, modes(rho_, 3)/3.0d0+modes(drho_, 2)
    end if
  end subroutine print_error

  subroutine set_gradrho(igrid,level,ixI^L,ixO^L,qt,w,x)
    integer, intent(in)             :: igrid,level,ixI^L,ixO^L
    double precision, intent(in)    :: qt,x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: tmp(ixI^S)

    ! here we compute \partial rho\partial x and store it in drho_
    call gradient(w(ixI^S,rho_),ixI^L,ixO^L,1,tmp)
    w(ixO^S,drho_) = tmp(ixO^S)

  end subroutine set_gradrho

end module mod_usr

