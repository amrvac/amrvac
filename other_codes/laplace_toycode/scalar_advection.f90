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
  integer, parameter :: i_dphi         = 3
  integer            :: n_iter   = 100
  logical            :: periodic(2)    = .true.
  integer            :: nx(2)   = 512
  integer            :: bx(2)  = 16
  integer            :: n_gc           = 2
  logical            :: write_output   = .false.
  real(dp)           :: v(2)           = [1.0_dp, 1.0_dp]
  real(dp)           :: grid_length(2) = [1.0_dp, 1.0_dp]
  real(dp)           :: dt
  real(dp)           :: cfl_number     = 0.5_dp

  integer            :: n, ii
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
  call CFG_add_get(cfg, 'v', v, 'Velocity')

  if (n_gc < 2) error stop "n_gc < 2"

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
     call advance_heuns_method(bg, dt)
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

  subroutine advance_heuns_method(bg, dt)
    type(block_grid_t), intent(inout) :: bg
    real(dp), intent(in)              :: dt

    call forward_euler(bg, bg%bx, bg%ilo, bg%ihi, bg%n_vars, bg%n_blocks, dt, bg%uu, 0, &
            1, [0], [1.0_dp], 1)
    call forward_euler(bg, bg%bx, bg%ilo, bg%ihi, bg%n_vars, bg%n_blocks, dt, bg%uu, 1, &
         2, [0, 1], [0.5_dp, 0.5_dp], 0)
  end subroutine advance_heuns_method

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

  subroutine forward_euler(bg, bx, lo, hi, n_vars, n_blocks, dt, uu, &
       s_deriv, n_prev, s_prev, w_prev, s_out)
    type(block_grid_t), intent(inout) :: bg
    integer, intent(in)     :: n_blocks, bx(2), lo(2), hi(2), n_vars
    real(dp), intent(in)    :: dt
    real(dp), intent(inout) :: uu(IX(lo(1):hi(1), lo(2):hi(2), n_blocks, n_vars))
    integer, intent(in)     :: s_deriv        !< State to compute derivatives from
    integer, intent(in)     :: n_prev         !< Number of previous states
    integer, intent(in)     :: s_prev(n_prev) !< Previous states
    real(dp), intent(in)    :: w_prev(n_prev) !< Weights of previous states
    integer, intent(in)     :: s_out          !< Output state
    integer                 :: n, i, j, m
    real(dp)                :: inv_dr(2)
    real(dp)                :: fx(2), fy(2), tmp(5), newval

    call update_ghostcells(bg, i_phi+s_deriv)

    inv_dr = 1/bg%dr

    !$acc parallel loop private(fx, fy, tmp)
    do n = 1, n_blocks
       !$acc loop
       do j = 1, bx(2)
          !$acc loop
          do i = 1, bx(1)
             tmp = uu(IX(i-2:i+2, j, n, i_phi+s_deriv))
             call muscl_flux(v(1), tmp, fx)

             tmp = uu(IX(i, j-2:j+2, n, i_phi+s_deriv))
             call muscl_flux(v(2), tmp, fy)

             uu(IX(i, j, n, i_dphi)) = dt * ((fx(1) - fx(2)) * inv_dr(1) + &
                  (fy(1) - fy(2)) * inv_dr(2))
          end do
       end do
    end do

    !$acc parallel loop private(newval, m)
    do n = 1, n_blocks
       do m = 1, n_prev
          !$acc loop
          do j = 1, bx(2)
             !$acc loop
             do i = 1, bx(1)
                ! Add weighted previous states
                uu(IX(i, j, n, i_dphi)) = uu(IX(i, j, n, i_dphi)) + &
                     uu(IX(i, j, n, i_phi+s_prev(m))) * w_prev(m)

             end do
          end do
       end do

       !$acc loop
       do j = 1, bx(2)
          !$acc loop
          do i = 1, bx(1)
             uu(IX(i, j, n, i_phi + s_out)) = uu(IX(i, j, n, i_dphi))
          end do
       end do
    end do
  end subroutine forward_euler

  subroutine upwind_flux(v, u, flux)
    !$acc routine seq
    real(dp), intent(in)  :: v
    real(dp), intent(in)  :: u(5)
    real(dp), intent(out) :: flux(2)

    if (v > 0) then
       flux(1) = v * u(2)
       flux(2) = v * u(3)
    else
       flux(1) = v * u(3)
       flux(2) = v * u(4)
    end if
  end subroutine upwind_flux

  subroutine muscl_flux(v, u, flux)
    !$acc routine seq
    real(dp), intent(in)  :: v
    real(dp), intent(in)  :: u(5)
    real(dp), intent(out) :: flux(2)
    real(dp)              :: u_diff(4), uL, uR

    u_diff = u(2:5) - u(1:4)

    uL = u(2) + 0.5_dp * vanleer(u_diff(1), u_diff(2))
    uR = u(3) - 0.5_dp * vanleer(u_diff(2), u_diff(3))
    flux(1) = 0.5 * (v*uL + v*uR - abs(v) * (uR-uL))

    uL = u(3) + 0.5_dp * vanleer(u_diff(2), u_diff(3))
    uR = u(4) - 0.5_dp * vanleer(u_diff(3), u_diff(4))
    flux(2) = 0.5 * (v*uL + v*uR - abs(v) * (uR-uL))

  end subroutine muscl_flux

  pure real(dp) function minmod(a, b)
    real(dp), intent(in) :: a, b

    if (a * b <= 0) then
       minmod = 0.0_dp
    else if (abs(a) < abs(b)) then
       minmod = a
    else
       minmod = b
    end if
  end function minmod

  pure real(dp) function vanleer(a, b) result(phi)
    real(dp), intent(in) :: a, b
    real(dp)             :: ab

    ab = a * b
    if (ab > 0) then
       phi = 2 * ab / (a + b)
    else
       phi = 0
    end if
  end function vanleer

end program scalar_advection
