#include "macros.h"
! Toy code for solving a 2D Laplace equation on a uniform grid with OpenACC. The
! grid is divided into blocks with a layer of ghost cells around them.
!
! Author: Jannis Teunissen (jannis.teunissen@cwi.nl)
program ghostcell
  use m_blockgrid

  implicit none
  integer, parameter :: dp = kind(0.0d0)
  integer, parameter :: NDIM = 2

  integer :: grid_size(NDIM)
  integer :: block_size(NDIM)
  integer :: n_variables
  integer :: n_iterations
  integer :: n_gc
  logical :: write_output

  integer :: n_args
  character(len=100) :: tmp

  type(block_grid_t) :: bg

  integer :: n
  integer :: i_new, i_old
  integer :: t_start, t_end, count_rate
  real(dp) :: t_total, unknowns_per_ns

  ! Set default values
  grid_size = 512
  block_size = 16
  n_gc = 1                      ! number of ghost cells
  n_iterations = 100
  write_output = .false.
  n_variables = 2

  ! Parse command line arguments
  n_args = command_argument_count()

  if (n_args == 0) then
     print *, "Usage: ./laplace_openacc_ghostcell Nx Ny block_Nx block_Ny"
     print *, "...    n_ghostcell n_iterations write_output"
  end if

  if (n_args >= 2) then
     call get_command_argument(1, tmp)
     read(tmp, *) grid_size(1)
     call get_command_argument(2, tmp)
     read(tmp, *) grid_size(2)
  end if

  if (n_args >= 4) then
     call get_command_argument(3, tmp)
     read(tmp, *) block_size(1)
     call get_command_argument(4, tmp)
     read(tmp, *) block_size(2)
  end if

  if (n_args >= 5) then
     call get_command_argument(5, tmp)
     read(tmp, *) n_gc
  end if

  if (n_args >= 6) then
     call get_command_argument(6, tmp)
     read(tmp, *) n_iterations
  end if

  if (n_args >= 7) then
     call get_command_argument(7, tmp)
     read(tmp, *) write_output
  end if

  if (n_args >= 8) then
     error stop "Too many arguments"
  end if

  call initialize_grid(grid_size, block_size, n_gc, 2, bg)
  call set_initial_conditions()

  if (write_output) &
       call write_brick_of_values(bg, "test_", "phi", 0, 0.0_dp, 1)

  ! Copy data structure and allocatable components to device
  !$acc enter data copyin(bg, bg%uu, bg%bnd_x, bg%bnd_y, bg%bc_xlo, &
  !$acc bg%bc_xhi, bg%bc_ylo, bg%bc_yhi)

  call system_clock(t_start, count_rate)
  do n = 1, n_iterations
     i_new = mod(n, 2) + 1
     i_old = mod(n+1, 2) + 1
     call update_ghostcells(bg, i_old)
     call update_solution(bg%bx, bg%n_gc, bg%n_vars, bg%n_blocks, i_new, i_old, bg%uu)
  end do
  call system_clock(t_end, count_rate)

  !$acc exit data copyout(bg%uu)

  if (write_output) &
       call write_brick_of_values(bg, "test_", "phi", 1, real(n, dp), i_new)

  t_total = (t_end - t_start)/real(count_rate, dp)
  unknowns_per_ns = 1e-9_dp * n_iterations/t_total * product(grid_size)

  write(*, "(5(A6,' '),2(A8,' '),2A12)") "nx", "ny", "box_nx", "box_ny", &
       "n_gc", "n_iter", "n_boxes", "t_total", "unknowns/ns"
  write(*, "(5(I6,' '),2(I8,' '),2F12.6)") bg%nx, bg%bx, bg%n_gc, &
       n_iterations, bg%n_blocks, t_total, unknowns_per_ns

contains

  subroutine set_initial_conditions()
    ! Zero everywhere
    bg%uu(IX(:, :, :, :)) = 0.0_dp

    ! Fix boundary value
    bg%bc(:) = 1.0_dp
  end subroutine set_initial_conditions

  subroutine update_solution(bx, n_gc, n_vars, n_blocks, i_new, i_old, uu)
    integer, intent(in) :: n_blocks, bx(2), n_gc, i_new, i_old, n_vars
    real(dp), intent(inout) :: uu(IX(1-n_gc:bx(1)+n_gc, 1-n_gc:bx(2)+n_gc, n_blocks, n_vars))
    integer             :: n, i, j

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
