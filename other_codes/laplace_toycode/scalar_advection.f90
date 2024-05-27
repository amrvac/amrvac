#include "macros.h"
! Author: Jannis Teunissen (jannis.teunissen@cwi.nl)
program scalar_advection
  use m_blockgrid
  use m_config

  implicit none
  integer, parameter :: dp = kind(0.0d0)

  type(block_grid_t) :: bg
  type(CFG_t)        :: cfg

  integer, parameter :: n_variables    = 3
  integer, parameter :: i_phi          = 1
  integer, parameter :: i_fx           = 2
  integer, parameter :: i_fy           = 3
  integer            :: n_iter   = 100
  logical            :: periodic(2)    = .true.
  integer            :: nx(2)   = 512
  integer            :: bx(2)  = 16
  integer            :: n_gc           = 1
  logical            :: write_output   = .false.
  real(dp)           :: v(2)           = [1.0_dp, 1.0_dp]
  real(dp)           :: grid_length(2) = [1.0_dp, 1.0_dp]
  real(dp)           :: dt
  real(dp)           :: cfl_number     = 0.5_dp

  integer            :: n
  integer            :: t_start, t_end, count_rate
  real(dp)           :: t_total, unknowns_per_ns

  call CFG_update_from_arguments(cfg)
  call CFG_add_get(cfg, 'n_iter', n_iter, 'Number of iterations')
  call CFG_add_get(cfg, 'periodic', periodic, 'Whether the domain is periodic')
  call CFG_add_get(cfg, 'nx', nx, 'Size of uniform grid')
  call CFG_add_get(cfg, 'bx', bx, 'Size of grid blocks')
  call CFG_add_get(cfg, 'n_gc', n_gc, 'Number of ghost cells')
  call CFG_add_get(cfg, 'write_output', write_output, 'Whether to write output')
  call CFG_add_get(cfg, 'cfl_number', cfl_number, 'CFL number to use')

  call initialize_grid(nx, grid_length, bx, n_gc, n_variables, &
       periodic, bg)
  call set_initial_conditions()

  dt = cfl_number / sum(v/bg%dr)
  write(*, "(A,E12.4)") "Using dt = ", dt

  if (write_output) then
     !$acc update self(bg%uu)
     call write_brick_of_values(bg, "test_", "phi", 0, 0.0_dp, i_phi)
  end if

  call system_clock(t_start, count_rate)
  do n = 1, n_iter
     call update_ghostcells(bg, i_phi)
     call compute_fluxes(bg%bx, bg%ilo, bg%ihi, bg%n_vars, bg%n_blocks, &
          bg%uu)
     call update_solution(bg%bx, bg%ilo, bg%ihi, bg%n_vars, bg%n_blocks, &
          dt, bg%uu)
  end do
  call system_clock(t_end, count_rate)

  if (write_output) then
     !$acc update self(bg%uu)
     call write_brick_of_values(bg, "test_", "phi", 1, n_iter*dt, i_phi)
  end if

  t_total = (t_end - t_start)/real(count_rate, dp)
  unknowns_per_ns = 1e-9_dp * n_iter/t_total * product(bg%nx)

  write(*, "(5(A6,' '),2(A8,' '),2A12)") "nx", "ny", "box_nx", "box_ny", &
       "n_gc", "n_iter", "n_boxes", "t_total", "unknowns/ns"
  write(*, "(5(I6,' '),2(I8,' '),2F12.6)") bg%nx, bg%bx, bg%n_gc, &
       n_iter, bg%n_blocks, t_total, unknowns_per_ns

contains

  subroutine set_initial_conditions()
    integer  :: n, i, j, ix(2)
    real(dp) :: r(2), radius

    radius = 0.25_dp

    !$acc parallel loop private(r, ix)
    do n = 1, bg%n_blocks
       !$acc loop
       do j = 1, bg%bx(2)
          !$acc loop
          do i = 1, bg%bx(1)
             ix = [i, j] + (bg%id_to_index(:, n) - 1) * bg%bx
             r = bg%dr * (ix - 0.5_dp)

             if (sum((r - 0.5_dp)**2) < radius**2) then
                bg%uu(IX(i, j, n, i_phi)) = 1.0_dp
             else
                bg%uu(IX(i, j, n, i_phi)) = 0.0_dp
             end if
          end do
       end do
    end do
  end subroutine set_initial_conditions

  subroutine compute_fluxes(bx, lo, hi, n_vars, n_blocks, uu)
    integer, intent(in)     :: n_blocks, bx(2), lo(2), hi(2), n_vars
    real(dp), intent(inout) :: uu(IX(lo(1):hi(1), lo(2):hi(2), n_blocks, n_vars))
    integer                 :: n, i, j

    !$acc parallel loop
    do n = 1, n_blocks
       !$acc loop
       do j = 1, bx(2)
          !$acc loop
          do i = 1, bx(1) + 1
             if (v(1) > 0) then
                uu(IX(i, j, n, i_fx)) = v(1) * uu(IX(i-1, j, n, i_phi))
             else
                uu(IX(i, j, n, i_fx)) = v(1) * uu(IX(i, j, n, i_phi))
             end if
          end do
       end do

       !$acc loop
       do j = 1, bx(2) + 1
          !$acc loop
          do i = 1, bx(1)
             if (v(2) > 0) then
                uu(IX(i, j, n, i_fy)) = v(2) * uu(IX(i, j-1, n, i_phi))
             else
                uu(IX(i, j, n, i_fy)) = v(2) * uu(IX(i, j, n, i_phi))
             end if
          end do
       end do
    end do
  end subroutine compute_fluxes

  subroutine update_solution(bx, lo, hi, n_vars, n_blocks, dt, uu)
    integer, intent(in)     :: n_blocks, bx(2), lo(2), hi(2), n_vars
    real(dp), intent(in)    :: dt
    real(dp), intent(inout) :: uu(IX(lo(1):hi(1), lo(2):hi(2), n_blocks, n_vars))
    integer                 :: n, i, j
    real(dp)                :: inv_dr(2)

    inv_dr = 1/bg%dr

    !$acc parallel loop
    do n = 1, n_blocks
       !$acc loop
       do j = 1, bx(2)
          !$acc loop
          do i = 1, bx(1)
             uu(IX(i, j, n, i_phi)) = uu(IX(i, j, n, i_phi)) + dt * (&
                  (uu(IX(i, j, n, i_fx)) - uu(IX(i+1, j, n, i_fx))) * inv_dr(1) + &
                  (uu(IX(i, j, n, i_fy)) - uu(IX(i, j+1, n, i_fy))) * inv_dr(2))
          end do
       end do
    end do
  end subroutine update_solution

end program scalar_advection
