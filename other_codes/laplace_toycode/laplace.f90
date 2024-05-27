#include "macros.h"
! Toy code for solving a 2D Laplace equation on a uniform grid with OpenACC. The
! grid is divided into blocks with a layer of ghost cells around them.
!
! Author: Jannis Teunissen (jannis.teunissen@cwi.nl)
program ghostcell
  use m_blockgrid
  use m_config

  implicit none
  integer, parameter :: dp = kind(0.0d0)

  type(block_grid_t) :: bg
  type(CFG_t)        :: cfg

  integer, parameter :: n_variables   = 2
  integer            :: n_iterations  = 100
  logical            :: periodic(2)   = .false.
  integer            :: grid_size(2)  = 512
  integer            :: block_size(2) = 16
  integer            :: n_gc          = 1
  logical            :: write_output  = .false.

  integer  :: n
  integer  :: i_new, i_old
  integer  :: t_start, t_end, count_rate
  real(dp) :: t_total, unknowns_per_ns

  call CFG_update_from_arguments(cfg)
  call CFG_add_get(cfg, 'n_iterations', n_iterations, 'Number of iterations')
  call CFG_add_get(cfg, 'periodic', periodic, 'Whether the domain is periodic')
  call CFG_add_get(cfg, 'grid_size', grid_size, 'Size of uniform grid')
  call CFG_add_get(cfg, 'block_size', block_size, 'Size of grid blocks')
  call CFG_add_get(cfg, 'n_gc', n_gc, 'Number of ghost cells')
  call CFG_add_get(cfg, 'write_output', write_output, 'Whether to write output')

  call initialize_grid(grid_size, block_size, n_gc, n_variables, periodic, bg)
  call set_initial_conditions()

  if (write_output) &
       call write_brick_of_values(bg, "test_", "phi", 0, 0.0_dp, 1)

  call system_clock(t_start, count_rate)
  do n = 1, n_iterations
     i_new = mod(n, 2) + 1
     i_old = mod(n+1, 2) + 1
     call update_ghostcells(bg, i_old)
     call update_solution(bg%bx, bg%ilo, bg%ihi, bg%n_vars, bg%n_blocks, &
          i_new, i_old, bg%uu)
  end do
  call system_clock(t_end, count_rate)

  !$acc exit data copyout(bg%uu)

  if (write_output) &
       call write_brick_of_values(bg, "test_", "phi", 1, real(n, dp), i_new)

  t_total = (t_end - t_start)/real(count_rate, dp)
  unknowns_per_ns = 1e-9_dp * n_iterations/t_total * product(bg%nx)

  write(*, "(5(A6,' '),2(A8,' '),2A12)") "nx", "ny", "box_nx", "box_ny", &
       "n_gc", "n_iter", "n_boxes", "t_total", "unknowns/ns"
  write(*, "(5(I6,' '),2(I8,' '),2F12.6)") bg%nx, bg%bx, bg%n_gc, &
       n_iterations, bg%n_blocks, t_total, unknowns_per_ns

contains

  subroutine set_initial_conditions()
    ! Zero everywhere
    !$acc kernels
    bg%uu(IX(:, :, :, :)) = 0.0_dp

    ! Fix boundary value
    bg%bc(:) = 1.0_dp
    !$acc end kernels
  end subroutine set_initial_conditions

  subroutine update_solution(bx, lo, hi, n_vars, n_blocks, i_new, i_old, uu)
    integer, intent(in)     :: n_blocks, bx(2), lo(2), hi(2), i_new, i_old, n_vars
    real(dp), intent(inout) :: uu(IX(lo(1):hi(1), lo(2):hi(2), n_blocks, n_vars))
    integer                 :: n, i, j

    !$acc parallel loop
    do n = 1, n_blocks
       !$acc loop
       do j = 1, bx(2)
          !$acc loop
          do i = 1, bx(1)
             uu(IX(i, j, n, i_new)) = 0.25_dp * ( &
                  uu(IX(i+1, j, n, i_old)) + uu(IX(i-1, j, n, i_old)) + &
                  uu(IX(i, j-1, n, i_old)) + uu(IX(i, j+1, n, i_old)))
          end do
       end do
    end do
  end subroutine update_solution

end program ghostcell
