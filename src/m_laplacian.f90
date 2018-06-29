#include "cpp_macros.h"
!> Module which contains multigrid procedures for a Laplacian operator
module m_laplacian
  use m_data_structures

  implicit none
  private

  public :: laplacian_set_methods

contains

  subroutine laplacian_set_methods(mg)
    type(mg_t), intent(inout) :: mg

    select case (mg%geometry_type)
    case (mg_cartesian)
       mg%box_op => box_lpl

       select case (mg%smoother_type)
       ! case (smoother_jacobi)
       !    mg%box_smoother => box_jacobi_lpl
       case (mg_smoother_gs, mg_smoother_gsrb)
          mg%box_smoother => box_gs_lpl
       case default
          error stop "laplacian_set_methods: unsupported smoother type"
       end select
    case (mg_cylindrical)
       if (NDIM == 3) error stop "Cylindrical 3D not supported yet"

       mg%box_op => box_clpl

       select case (mg%smoother_type)
       case (mg_smoother_gs, mg_smoother_gsrb)
          mg%box_smoother => box_gs_clpl
       case default
          error stop "laplacian_set_methods: unsupported smoother type"
       end select
    case default
       error stop "laplacian_set_methods: unsupported geometry"
    end select

  end subroutine laplacian_set_methods

  !> Perform Gauss-Seidel relaxation on box for a Laplacian operator
  subroutine box_gs_lpl(mg, id, nc, redblack_cntr)
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
    fac = 0.5_dp / sum(idr2)
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
  end subroutine box_gs_lpl

!   !> Perform Jacobi relaxation on box for a Laplacian operator
!   subroutine box_jacobi_lpl(mg, id, nc, cntr)
!     type(mg_t), intent(inout) :: mg
!     integer, intent(in)       :: id
!     integer, intent(in)       :: nc
!     integer, intent(in)       :: cntr !< Not used
!     integer                   :: IJK
!     real(dp), parameter       :: w     = 2.0_dp / 3
!     real(dp)                  :: tmp(DTIMES(0:nc+1))
!     real(dp)                  :: dr2
! #if NDIM == 3
!     real(dp), parameter       :: sixth = 1/6.0_dp
! #endif

!     dr2   = mg%dr(mg%boxes(id)%lvl)**2

!     associate (box => mg%boxes(id))
!       tmp = box%cc(DTIMES(:), mg_iphi)
!       do KJI_DO(1, nc)
! #if NDIM == 2
!          box%cc(i, j, mg_iphi) = (1-w) * box%cc(i, j, mg_iphi) + &
!               0.25_dp * w * ( &
!               tmp(i+1, j) + tmp(i-1, j) + &
!               tmp(i, j+1) + tmp(i, j-1) - &
!               dr2 * box%cc(i, j, mg_irhs))
! #elif NDIM == 3
!          box%cc(i, j, k, mg_iphi) = (1-w) * &
!               box%cc(i, j, k, mg_iphi) + &
!               sixth * w * ( &
!               tmp(i+1, j, k) + tmp(i-1, j, k) + &
!               tmp(i, j+1, k) + tmp(i, j-1, k) + &
!               tmp(i, j, k+1) + tmp(i, j, k-1) - &
!               dr2 * box%cc(i, j, k, mg_irhs))
! #endif
!       end do; CLOSE_DO
!     end associate
!   end subroutine box_jacobi_lpl

  !> Perform Laplacian operator on a box
  subroutine box_lpl(mg, id, nc, i_out)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: id
    integer, intent(in)       :: nc
    integer, intent(in)       :: i_out !< Index of variable to store Laplacian in
    integer                   :: IJK
    real(dp)                  :: idr2(NDIM)

    idr2 = 1 / mg%dr(:, mg%boxes(id)%lvl)**2

    associate (cc => mg%boxes(id)%cc, n => mg_iphi)
#if NDIM == 2
      do j = 1, nc
         do i = 1, nc
            cc(i, j, i_out) = &
                 idr2(1) * (cc(i-1, j, n) + cc(i+1, j, n) - 2 * cc(i, j, n)) + &
                 idr2(2) * (cc(i, j-1, n) + cc(i, j+1, n) - 2 * cc(i, j, n))
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
                    - 2 * cc(i, j, k, n))
            end do
         end do
      end do
#endif
    end associate
  end subroutine box_lpl

  !> Perform Laplacian operator on a box in cylindrical geometry, using (r,z)
  !> and (r,phi,z) coordinates in 2D/3D.
  subroutine box_clpl(mg, id, nc, i_out)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: id
    integer, intent(in)       :: nc
    integer, intent(in)       :: i_out !< Index of variable to store Laplacian in
    integer                   :: IJK
    real(dp)                  :: dr(NDIM), idr2(NDIM)
    real(dp)                  :: r_face(nc+1), r_inv(nc)

    dr     = mg%dr(:, mg%boxes(id)%lvl)
    idr2   = 1 / dr**2
    r_face = mg%boxes(id)%r_min(1) + dr(1) * [(i, i=0,nc)]
    r_inv  = 1/(mg%boxes(id)%r_min(1) + dr(1) * [(i-0.5_dp, i=1,nc)])

    associate (cc => mg%boxes(id)%cc, n => mg_iphi)
#if NDIM == 2
      do j = 1, nc
         do i = 1, nc
            cc(i, j, i_out) = idr2(1) * (&
                 r_face(i) * r_inv(i) * cc(i-1, j, n) + &
                 r_face(i+1) * r_inv(i) * cc(i+1, j, n) &
                 - 2 * cc(i, j, n)) &
                 + idr2(2) * (cc(i, j-1, n) +  cc(i, j+1, n) - 2 * cc(i, j, n))
         end do
      end do
#elif NDIM == 3
      do k = 1, nc
         do j = 1, nc
            do i = 1, nc
               cc(i, j, k, i_out) = idr2(1) * (&
                    r_face(i) * r_inv(i) * cc(i-1, j, k, n) + &
                    r_face(i+1) * r_inv(i) * cc(i+1, j, k, n) &
                    - 2 * cc(i, j, k, n)) + &
                    idr2(2) * r_inv(i)**2 * &
                    (cc(i, j+1, k, n) + cc(i, j-1, k, n) - 2 * cc(i, j, k, n)) + &
                    idr2(3) * &
                    (cc(i, j, k-1, n) + cc(i, j, k+1, n) - 2 * cc(i, j, k, n))
            end do
         end do
      end do
#endif
    end associate
  end subroutine box_clpl

  !> Perform Gauss-Seidel relaxation on box for a Laplacian operator in
  !> cylindrical geometry. TODO: in 3D this does not converge well, maybe it
  !> will for a stretched grid.
  subroutine box_gs_clpl(mg, id, nc, redblack_cntr)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: id
    integer, intent(in)       :: nc
    integer, intent(in)       :: redblack_cntr !< Iteration counter
    integer                   :: IJK, i0, di
    logical                   :: redblack
    real(dp)                  :: idr2(NDIM), dr(NDIM), dr2(NDIM), fac
#if NDIM == 3
    real(dp), parameter       :: sixth = 1/6.0_dp
#endif
    real(dp)                  :: r_face(nc+1), r_inv(nc)

    dr     = mg%dr(:, mg%boxes(id)%lvl)
    dr2    = dr**2
    idr2   = 1/dr**2
    fac    = 0.5_dp / sum(idr2)
    r_face = mg%boxes(id)%r_min(1) + dr(1) * [(i, i=0,nc)]
    r_inv  = 1/(mg%boxes(id)%r_min(1) + dr(1) * [(i-0.5_dp, i=1,nc)])

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
            cc(i, j, n) = fac * (idr2(1) * &
                 (r_face(i+1) * r_inv(i) * cc(i+1, j, n) + &
                 r_face(i) * r_inv(i) * cc(i-1, j, n)) + &
                 idr2(2) * (cc(i, j+1, n) + cc(i, j-1, n)) &
                 - cc(i, j, mg_irhs))
         end do
      end do
#elif NDIM == 3
      do k = 1, nc
         do j = 1, nc
            if (redblack) &
                 i0 = 2 - iand(ieor(redblack_cntr, k+j), 1)
            do i = i0, nc, di
               cc(i, j, k, n) = (idr2(1) * ( &
                    r_face(i+1) * r_inv(i) * cc(i+1, j, k, n) + &
                    r_face(i) * r_inv(i) * cc(i-1, j, k, n)) + &
                    idr2(2) * r_inv(i)**2 * ( &
                    cc(i, j+1, k, n) + cc(i, j-1, k, n)) + &
                    idr2(3) * (cc(i, j, k+1, n) + &
                    cc(i, j, k-1, n)) - cc(i, j, k, mg_irhs)) * &
                    0.5_dp / (idr2(1) + idr2(2) * r_inv(i)**2 + idr2(3))
            end do
         end do
      end do
#endif
    end associate
  end subroutine box_gs_clpl

end module m_laplacian
