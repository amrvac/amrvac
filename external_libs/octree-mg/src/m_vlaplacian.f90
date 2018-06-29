#include "cpp_macros.h"
!> Module which contains multigrid procedures for a variable-coefficient
!> Laplacian operator, assuming the variation is smooth
module m_vlaplacian
  use m_data_structures

  implicit none
  private

  public :: vlaplacian_set_methods

contains

  subroutine vlaplacian_set_methods(mg)
    type(mg_t), intent(inout) :: mg

    if (mg%n_extra_vars == 0 .and. mg%is_allocated) then
       error stop "vlaplacian_set_methods: mg%n_extra_vars == 0"
    else
       mg%n_extra_vars = max(1, mg%n_extra_vars)
    end if

    ! Use Neumann zero boundary conditions for the variable coefficient, since
    ! it is needed in ghost cells.
    mg%bc(:, mg_iveps)%bc_type = mg_bc_neumann
    mg%bc(:, mg_iveps)%bc_value = 0.0_dp

    select case (mg%geometry_type)
    case (mg_cartesian)
       mg%box_op => box_vlpl

       ! TODO: test whether a custom prolongation works better
       ! mg%box_prolong => vlpl_prolong

       select case (mg%smoother_type)
       case (mg_smoother_gs, mg_smoother_gsrb)
          mg%box_smoother => box_gs_vlpl
       case default
          error stop "vlaplacian_set_methods: unsupported smoother type"
       end select

    case default
       error stop "vlaplacian_set_methods: unsupported geometry"
    end select

  end subroutine vlaplacian_set_methods

  !> Perform Gauss-Seidel relaxation on box for a Laplacian operator
  subroutine box_gs_vlpl(mg, id, nc, redblack_cntr)
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
                 (sum(c(:) * u(:)) - cc(i, j, mg_irhs)) / sum(c(:))
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
                    cc(i, j, k, mg_irhs)) / sum(c(:))
            end do
         end do
      end do
#endif
    end associate
  end subroutine box_gs_vlpl

  !> Perform Laplacian operator on a box
  subroutine box_vlpl(mg, id, nc, i_out)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: id
    integer, intent(in)       :: nc
    integer, intent(in)       :: i_out !< Index of variable to store Laplacian in
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
                 a0*a(:)/(a0 + a(:)) * (u(:) - u0))
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
                    a0*a(:)/(a0 + a(:)) * (u(:) - u0))
            end do
         end do
      end do
#endif
    end associate
  end subroutine box_vlpl

  !> Prolong from a parent to a child with index offset dix. This method could
  !> sometimes work better than the default prolongation, which does not take
  !> the variation in epsilon into account
  subroutine vlpl_prolong(mg, p_id, dix, nc, iv, fine)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: p_id             !< Id of parent
    integer, intent(in)       :: dix(NDIM)        !< Offset of child in parent grid
    integer, intent(in)       :: nc               !< Child grid size
    integer, intent(in)       :: iv               !< Prolong from this variable
    real(dp), intent(out)     :: fine(DTIMES(nc)) !< Prolonged values

    integer  :: IJK, hnc
#if NDIM == 2
    integer  :: ic, jc
    real(dp) :: f0, flx, fhx, fly, fhy
    real(dp) :: c0, clx, chx, cly, chy
#elif NDIM == 3
    integer  :: ic, jc, kc
    real(dp) :: f0, flx, fhx, fly, fhy, flz, fhz
    real(dp) :: c0, clx, chx, cly, chy, clz, chz
#endif

    hnc = nc/2

    associate (crs => mg%boxes(p_id)%cc, &
         i_eps => mg_iveps)
#if NDIM == 2
      do j = 1, hnc
         jc = j + dix(2)
         do i = 1, hnc
            ic = i + dix(1)

            f0  = crs(ic, jc, iv)
            flx = crs(ic-1, jc, iv)
            fhx = crs(ic+1, jc, iv)
            fly = crs(ic, jc-1, iv)
            fhy = crs(ic, jc+1, iv)

            c0  = 2 * crs(ic, jc, i_eps)
            clx = crs(ic-1, jc, i_eps)
            chx = crs(ic+1, jc, i_eps)
            cly = crs(ic, jc-1, i_eps)
            chy = crs(ic, jc+1, i_eps)

            fine(2*i-1, 2*j-1) = (f0*c0 + flx*clx + fly*cly) / &
                 (c0 + clx + cly)
            fine(2*i  , 2*j-1) = (f0*c0 + fhx*chx + fly*cly) / &
                 (c0 + chx + cly)
            fine(2*i-1, 2*j)   = (f0*c0 + flx*clx + fhy*chy) / &
                 (c0 + clx + chy)
            fine(2*i  , 2*j)   = (f0*c0 + fhx*chx + fhy*chy) / &
                 (c0 + chx + chy)
         end do
      end do
#elif NDIM == 3
      do k = 1, hnc
         kc = k + dix(3)
         do j = 1, hnc
            jc = j + dix(2)
            do i = 1, hnc
               ic = i + dix(1)

               f0  = crs(ic, jc, kc, iv)
               flx = crs(ic-1, jc, kc, iv)
               fhx = crs(ic+1, jc, kc, iv)
               fly = crs(ic, jc-1, kc, iv)
               fhy = crs(ic, jc+1, kc, iv)
               flz = crs(ic, jc, kc-1, iv)
               fhz = crs(ic, jc, kc+1, iv)

               c0  = crs(ic, jc, kc, i_eps)
               clx = crs(ic-1, jc, kc, i_eps)
               chx = crs(ic+1, jc, kc, i_eps)
               cly = crs(ic, jc-1, kc, i_eps)
               chy = crs(ic, jc+1, kc, i_eps)
               clz = crs(ic, jc, kc-1, i_eps)
               chz = crs(ic, jc, kc+1, i_eps)

               fine(2*i-1, 2*j-1, 2*k-1) = (f0*c0 + flx*clx + fly*cly + flz*clz) / &
                    (c0 + clx + cly + clz)
               fine(2*i, 2*j-1, 2*k-1)   = (f0*c0 + fhx*chx + fly*cly + flz*clz) / &
                    (c0 + chx + cly + clz)
               fine(2*i-1, 2*j, 2*k-1)   = (f0*c0 + flx*clx + fhy*chy + flz*clz) / &
                    (c0 + clx + chy + clz)
               fine(2*i, 2*j, 2*k-1)     = (f0*c0 + fhx*chx + fhy*chy + flz*clz) / &
                    (c0 + chx + chy + clz)
               fine(2*i-1, 2*j-1, 2*k)   = (f0*c0 + flx*clx + fly*cly + fhz*chz) / &
                    (c0 + clx + cly + chz)
               fine(2*i, 2*j-1, 2*k)     = (f0*c0 + fhx*chx + fly*cly + fhz*chz) / &
                    (c0 + chx + cly + chz)
               fine(2*i-1, 2*j, 2*k)     = (f0*c0 + flx*clx + fhy*chy + fhz*chz) / &
                    (c0 + clx + chy + chz)
               fine(2*i, 2*j, 2*k)       = (f0*c0 + fhx*chx + fhy*chy + fhz*chz) / &
                    (c0 + chx + chy + chz)
            end do
         end do
      end do
#endif
    end associate
  end subroutine vlpl_prolong

end module m_vlaplacian
