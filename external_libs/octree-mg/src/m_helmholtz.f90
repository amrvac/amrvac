#include "cpp_macros.h"
!> Module which contains multigrid procedures for a Helmholtz operator of the
!> form: laplacian(phi) - lambda*phi = f
module m_helmholtz
  use m_data_structures

  implicit none
  private

  !> The lambda used for the Helmholtz equation (should be >= 0)
  real(dp), public, protected :: helmholtz_lambda = 0.0_dp

  public :: helmholtz_set_methods
  public :: helmholtz_set_lambda

contains

  subroutine helmholtz_set_methods(mg)
    type(mg_t), intent(inout) :: mg

    select case (mg%geometry_type)
    case (mg_cartesian)
       mg%box_op => box_helmh

       select case (mg%smoother_type)
       case (mg_smoother_gs, mg_smoother_gsrb)
          mg%box_smoother => box_gs_helmh
       case default
          error stop "helmholtz_set_methods: unsupported smoother type"
       end select
    case default
       error stop "helmholtz_set_methods: unsupported geometry"
    end select

  end subroutine helmholtz_set_methods

  subroutine helmholtz_set_lambda(lambda)
    real(dp), intent(in) :: lambda

    if (lambda < 0) &
         error stop "helmholtz_set_lambda: lambda < 0 not allowed"

    helmholtz_lambda = lambda
  end subroutine helmholtz_set_lambda

  !> Perform Gauss-Seidel relaxation on box for a Helmholtz operator
  subroutine box_gs_helmh(mg, id, nc, redblack_cntr)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: id
    integer, intent(in)       :: nc
    integer, intent(in)       :: redblack_cntr !< Iteration counter
    integer                   :: IJK, i0, di
    real(dp)                  :: idr2(NDIM), fac
    logical                   :: redblack
#if NDIM == 3
    real(dp), parameter       :: sixth = 1/6.0_dp
#endif

    idr2 = 1/mg%dr(:, mg%boxes(id)%lvl)**2
    fac = 1.0_dp / (2 * sum(idr2) + helmholtz_lambda)
    i0  = 1

    redblack = (mg%smoother_type == mg_smoother_gsrb)
    if (redblack) then
       di = 2
    else
       di = 1
    end if

    ! The parity of redblack_cntr determines which cells we use. If
    ! redblack_cntr is even, we use the even cells and vice versa.
    associate (cc => mg%boxes(id)%cc, n => mg_iphi)
#if NDIM == 2
      do j = 1, nc
         if (redblack) &
              i0 = 2 - iand(ieor(redblack_cntr, j), 1)

         do i = i0, nc, di
            cc(i, j, n) = fac * ( &
                 idr2(1) * (cc(i+1, j, n) + cc(i-1, j, n)) + &
                 idr2(2) * (cc(i, j+1, n) + cc(i, j-1, n)) - &
                 cc(i, j, mg_irhs))
         end do
      end do
#elif NDIM == 3
      do k = 1, nc
         do j = 1, nc
            if (redblack) &
                 i0 = 2 - iand(ieor(redblack_cntr, k+j), 1)
            do i = i0, nc, di
               cc(i, j, k, n) = fac * ( &
                    idr2(1) * (cc(i+1, j, k, n) + cc(i-1, j, k, n)) + &
                    idr2(2) * (cc(i, j+1, k, n) + cc(i, j-1, k, n)) + &
                    idr2(3) * (cc(i, j, k+1, n) + cc(i, j, k-1, n)) - &
                    cc(i, j, k, mg_irhs))
            end do
         end do
      end do
#endif
    end associate
  end subroutine box_gs_helmh

  !> Perform Helmholtz operator on a box
  subroutine box_helmh(mg, id, nc, i_out)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: id
    integer, intent(in)       :: nc
    integer, intent(in)       :: i_out !< Index of variable to store Helmholtz in
    integer                   :: IJK
    real(dp)                  :: idr2(NDIM)

    idr2 = 1 / mg%dr(:, mg%boxes(id)%lvl)**2

    associate (cc => mg%boxes(id)%cc, n => mg_iphi)
#if NDIM == 2
      do j = 1, nc
         do i = 1, nc
            cc(i, j, i_out) = &
                 idr2(1) * (cc(i-1, j, n) + cc(i+1, j, n) - 2 * cc(i, j, n)) + &
                 idr2(2) * (cc(i, j-1, n) + cc(i, j+1, n) - 2 * cc(i, j, n)) - &
                 helmholtz_lambda * cc(i, j, n)
         end do
      end do
#elif NDIM == 3
      do k = 1, nc
         do j = 1, nc
            do i = 1, nc
               cc(i, j, k, i_out) = &
                    idr2(1) * (cc(i-1, j, k, n) + cc(i+1, j, k, n) &
                    - 2 * cc(i, j, k, n)) &
                    + idr2(2) * (cc(i, j-1, k, n) + cc(i, j+1, k, n) &
                    - 2 * cc(i, j, k, n)) &
                    + idr2(3) * (cc(i, j, k-1, n) + cc(i, j, k+1, n) &
                    - 2 * cc(i, j, k, n)) - &
                    helmholtz_lambda * cc(i, j, k, n)
            end do
         end do
      end do
#endif
    end associate
  end subroutine box_helmh

end module m_helmholtz
