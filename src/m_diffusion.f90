#include "cpp_macros.h"
!> Module for implicitly solving diffusion equations
module m_diffusion
  use m_data_structures
  use m_multigrid

  implicit none
  private

  public :: diffusion_solve
  public :: diffusion_solve_vcoeff

contains

  !> Solve a diffusion equation implicitly, assuming a constant diffusion
  !> coefficient. The solution at time t should be stored in mg_iphi, which is
  !> on output replaced by the solution at time t+dt.
  subroutine diffusion_solve(mg, dt, diffusion_coeff, order, max_res)
    use m_helmholtz
    type(mg_t), intent(inout) :: mg
    real(dp), intent(in)      :: dt
    real(dp), intent(in)      :: diffusion_coeff
    integer, intent(in)       :: order
    real(dp), intent(in)      :: max_res
    integer, parameter        :: max_its = 10
    integer                   :: n
    real(dp)                  :: res

    mg%operator_type = mg_helmholtz
    call mg_set_methods(mg)

    select case (order)
    case (1)
       call helmholtz_set_lambda(1/(dt * diffusion_coeff))
       call set_rhs(mg, -1/(dt * diffusion_coeff), 0.0_dp)
    case (2)
       call helmholtz_set_lambda(0.0d0)
       call mg_apply_op(mg, mg_irhs)
       call helmholtz_set_lambda(2/(dt * diffusion_coeff))
       call set_rhs(mg, -2/(dt * diffusion_coeff), -1.0_dp)
    case default
       error stop "diffusion_solve: order should be 1 or 2"
    end select

    ! Start with an FMG cycle
    call mg_fas_fmg(mg, .true., max_res=res)

    ! Add V-cycles if necessary
    do n = 1, max_its
       if (res <= max_res) exit
       call mg_fas_vcycle(mg, max_res=res)
    end do

    if (n == max_its + 1) then
       if (mg%my_rank == 0) &
            print *, "Did you specify boundary conditions correctly?"
       error stop "diffusion_solve: no convergence"
    end if
  end subroutine diffusion_solve

  !> Solve a diffusion equation implicitly, assuming a variable diffusion
  !> coefficient which has been stored in mg_iveps (also on coarse grids). The
  !> solution at time t should be stored in mg_iphi, which is on output replaced
  !> by the solution at time t+dt.
  subroutine diffusion_solve_vcoeff(mg, dt, order, max_res)
    use m_vhelmholtz
    type(mg_t), intent(inout) :: mg
    real(dp), intent(in)      :: dt
    integer, intent(in)       :: order
    real(dp), intent(in)      :: max_res
    integer, parameter        :: max_its = 10
    integer                   :: n
    real(dp)                  :: res

    mg%operator_type = mg_vhelmholtz
    call mg_set_methods(mg)

    select case (order)
    case (1)
       call vhelmholtz_set_lambda(1/dt)
       call set_rhs(mg, -1/dt, 0.0_dp)
    case (2)
       call vhelmholtz_set_lambda(0.0d0)
       call mg_apply_op(mg, mg_irhs)
       call vhelmholtz_set_lambda(2/dt)
       call set_rhs(mg, -2/dt, -1.0_dp)
    case default
       error stop "diffusion_solve: order should be 1 or 2"
    end select

    ! Start with an FMG cycle
    call mg_fas_fmg(mg, .true., max_res=res)

    ! Add V-cycles if necessary
    do n = 1, max_its
       if (res <= max_res) exit
       call mg_fas_vcycle(mg, max_res=res)
    end do

    if (n == max_its + 1) then
       if (mg%my_rank == 0) then
          print *, "Did you specify boundary conditions correctly?"
          print *, "Or is the variation in diffusion too large?"
       end if
       error stop "diffusion_solve: no convergence"
    end if
  end subroutine diffusion_solve_vcoeff

  subroutine set_rhs(mg, f1, f2)
    type(mg_t), intent(inout) :: mg
    real(dp), intent(in)      :: f1, f2
    integer                   :: lvl, i, id, nc

    do lvl = 1, mg%highest_lvl
       nc = mg%box_size_lvl(lvl)
       do i = 1, size(mg%lvls(lvl)%my_leaves)
          id = mg%lvls(lvl)%my_leaves(i)
          mg%boxes(id)%cc(DTIMES(1:nc), mg_irhs) = &
               f1 * mg%boxes(id)%cc(DTIMES(1:nc), mg_iphi) + &
               f2 * mg%boxes(id)%cc(DTIMES(1:nc), mg_irhs)
       end do
    end do
  end subroutine set_rhs

end module m_diffusion
