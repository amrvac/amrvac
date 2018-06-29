#include "cpp_macros.h"
!> Module which contains multigrid procedures for a Helmholtz operator of the
!> form: div(D grad(phi)) - lambda*phi = f, where D has a smooth spatial
!> variation
module m_vhelmholtz
  use m_data_structures

  implicit none
  private

  !> The lambda used for the Helmholtz equation (should be >= 0)
  real(dp), public, protected :: vhelmholtz_lambda = 0.0_dp

  public :: vhelmholtz_set_methods
  public :: vhelmholtz_set_lambda

contains

  subroutine vhelmholtz_set_methods(mg)
    type(mg_t), intent(inout) :: mg

    if (mg%n_extra_vars == 0 .and. mg%is_allocated) then
       error stop "vhelmholtz_set_methods: mg%n_extra_vars == 0"
    else
       mg%n_extra_vars = max(1, mg%n_extra_vars)
    end if

    ! Use Neumann zero boundary conditions for the variable coefficient, since
    ! it is needed in ghost cells.
    mg%bc(:, mg_iveps)%bc_type = mg_bc_neumann
    mg%bc(:, mg_iveps)%bc_value = 0.0_dp

    select case (mg%geometry_type)
    case (mg_cartesian)
       mg%box_op => box_vhelmh

       select case (mg%smoother_type)
       case (mg_smoother_gs, mg_smoother_gsrb)
          mg%box_smoother => box_gs_vhelmh
       case default
          error stop "vhelmholtz_set_methods: unsupported smoother type"
       end select
    case default
       error stop "vhelmholtz_set_methods: unsupported geometry"
    end select

  end subroutine vhelmholtz_set_methods

  subroutine vhelmholtz_set_lambda(lambda)
    real(dp), intent(in) :: lambda

    if (lambda < 0) &
         error stop "vhelmholtz_set_lambda: lambda < 0 not allowed"

    vhelmholtz_lambda = lambda
  end subroutine vhelmholtz_set_lambda

  !> Perform Gauss-Seidel relaxation on box for a Helmholtz operator
  subroutine box_gs_vhelmh(mg, id, nc, redblack_cntr)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: id
    integer, intent(in)       :: nc
    integer, intent(in)       :: redblack_cntr !< Iteration counter
    integer                   :: IJK, i0, di
    logical                   :: redblack
    real(dp)                  :: idr2(2*NDIM), u(2*NDIM)
    real(dp)                  :: a0, a(2*NDIM), c(2*NDIM)

    ! Duplicate 1/dr^2 array to multiply neighbor values
    idr2(1:2*NDIM:2) = 1/mg%dr(:, mg%boxes(id)%lvl)**2
    idr2(2:2*NDIM:2) = idr2(1:2*NDIM:2)
    i0  = 1

    redblack = (mg%smoother_type == mg_smoother_gsrb)
    if (redblack) then
       di = 2
    else
       di = 1
    end if

    ! The parity of redblack_cntr determines which cells we use. If
    ! redblack_cntr is even, we use the even cells and vice versa.
    associate (cc => mg%boxes(id)%cc, n => mg_iphi, &
         i_eps => mg_iveps)
#if NDIM == 2
      do j = 1, nc
         if (redblack) &
              i0 = 2 - iand(ieor(redblack_cntr, j), 1)

         do i = i0, nc, di
            a0     = cc(i, j, i_eps)
            u(1:2) = cc(i-1:i+1:2, j, n)
            a(1:2) = cc(i-1:i+1:2, j, i_eps)
            u(3:4) = cc(i, j-1:j+1:2, n)
            a(3:4) = cc(i, j-1:j+1:2, i_eps)
            c(:)   = 2 * a0 * a(:) / (a0 + a(:)) * idr2

            cc(i, j, n) = &
                 (sum(c(:) * u(:)) - cc(i, j, mg_irhs)) / &
                 (sum(c(:)) + vhelmholtz_lambda)
         end do
      end do
#elif NDIM == 3
      do k = 1, nc
         do j = 1, nc
            if (redblack) &
                 i0 = 2 - iand(ieor(redblack_cntr, k+j), 1)
            do i = i0, nc, di
               a0     = cc(i, j, k, i_eps)
               u(1:2) = cc(i-1:i+1:2, j, k, n)
               a(1:2) = cc(i-1:i+1:2, j, k, i_eps)
               u(3:4) = cc(i, j-1:j+1:2, k, n)
               a(3:4) = cc(i, j-1:j+1:2, k, i_eps)
               u(5:6) = cc(i, j, k-1:k+1:2, n)
               a(5:6) = cc(i, j, k-1:k+1:2, i_eps)
               c(:)   = 2 * a0 * a(:) / (a0 + a(:)) * idr2

               cc(i, j, k, n) = (sum(c(:) * u(:)) - &
                    cc(i, j, k, mg_irhs)) / &
                    (sum(c(:)) + vhelmholtz_lambda)
            end do
         end do
      end do
#endif
    end associate
  end subroutine box_gs_vhelmh

  !> Perform Helmholtz operator on a box
  subroutine box_vhelmh(mg, id, nc, i_out)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: id
    integer, intent(in)       :: nc
    integer, intent(in)       :: i_out !< Index of variable to store Helmholtz in
    integer                   :: IJK
    real(dp)                  :: idr2(2*NDIM), a0, u0, u(2*NDIM), a(2*NDIM)

    ! Duplicate 1/dr^2 array to multiply neighbor values
    idr2(1:2*NDIM:2) = 1/mg%dr(:, mg%boxes(id)%lvl)**2
    idr2(2:2*NDIM:2) = idr2(1:2*NDIM:2)

    associate (cc => mg%boxes(id)%cc, n => mg_iphi, &
         i_eps => mg_iveps)
#if NDIM == 2
      do j = 1, nc
         do i = 1, nc
            a0     = cc(i, j, i_eps)
            a(1:2) = cc(i-1:i+1:2, j, i_eps)
            a(3:4) = cc(i, j-1:j+1:2, i_eps)
            u0     = cc(i, j, n)
            u(1:2) = cc(i-1:i+1:2, j, n)
            u(3:4) = cc(i, j-1:j+1:2, n)

            cc(i, j, i_out) = sum(2 * idr2 * &
                 a0*a(:)/(a0 + a(:)) * (u(:) - u0)) - &
                 vhelmholtz_lambda * u0
         end do
      end do
#elif NDIM == 3
      do k = 1, nc
         do j = 1, nc
            do i = 1, nc
               u0 = cc(i, j, k, n)
               a0 = cc(i, j, k, i_eps)
               u(1:2) = cc(i-1:i+1:2, j, k, n)
               u(3:4) = cc(i, j-1:j+1:2, k, n)
               u(5:6) = cc(i, j, k-1:k+1:2, n)
               a(1:2) = cc(i-1:i+1:2, j, k, i_eps)
               a(3:4) = cc(i, j-1:j+1:2, k, i_eps)
               a(5:6) = cc(i, j, k-1:k+1:2, i_eps)

               cc(i, j, k, i_out) = sum(2 * idr2 * &
                    a0*a(:)/(a0 + a(:)) * (u(:) - u0)) - &
                    vhelmholtz_lambda * u0
            end do
         end do
      end do
#endif
    end associate
  end subroutine box_vhelmh

end module m_vhelmholtz
