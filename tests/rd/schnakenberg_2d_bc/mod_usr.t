module mod_usr
  use mod_rd

  implicit none

  integer, parameter :: dp = kind(0.0d0)

contains

  subroutine usr_init()

    usr_init_one_grid => schnakenberg_init
    rd_mg_bc(1, 1)%boundary_cond => my_custom_multigrid_bc
    usr_special_bc => specialbc_usr

    call rd_activate()

  end subroutine usr_init

  ! initialize one grid
  subroutine schnakenberg_init(ixG^L,ix^L,w,x)
    integer, intent(in)             :: ixG^L, ix^L
    double precision, intent(in)    :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)
    double precision                :: x1, x2
    double precision                :: l1, l2

    x1 = xprobmin1 + (xprobmax1 - xprobmin1)/3
    x2 = xprobmin2 + 0.5d0 * (xprobmax2 - xprobmin2)
    l1 = xprobmax1 - xprobmin1
    l2 = xprobmax2 - xprobmin2

    w(ix^S,u_) = sb_alpha + sb_beta + 1d-3 * &
         exp(-100d0 * ((x(ix^S, 1) - x1)**2 + (x(ix^S, 2) - x2)**2))
    w(ix^S,v_) = sb_beta / (sb_alpha + sb_beta)**2

  end subroutine schnakenberg_init

  !> To fill ghost cells near physical boundaries
  subroutine my_custom_multigrid_bc(box, nc, iv, nb, bc_type, bc)
    use mod_multigrid_coupling
    type(mg_box_t), intent(in)    :: box
    integer, intent(in)           :: nc
    integer, intent(in)           :: iv      !< Index of variable
    integer, intent(in)           :: nb      !< Direction
    integer, intent(out)          :: bc_type !< Type of b.c.
    double precision, intent(out) :: bc(nc^DE&) ! ^DE^ means repeat ndim-1 times
    double precision              :: rr(nc^DE&, ^ND)
    integer                       :: i, ixs(ndim-1), nb_dim, n_bc

    ! Type of boundary condition
    bc_type = mg_bc_dirichlet

    ! Get the coordinates at the cell faces at the boundary
    call mg_get_face_coords(box, nb, nc, rr)

    bc = solution(rr(:, 1), rr(:, 2), global_time)

  end subroutine my_custom_multigrid_bc

  subroutine specialbc_usr(qt,ixG^L,ixO^L,iB,w,x)
    integer, intent(in)             :: ixG^L, ixO^L, iB
    double precision, intent(in)    :: qt, x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    w(ixO^S, u_) = solution(x(ixO^S, 1), x(ixO^S, 2), qt)

  end subroutine specialbc_usr

  elemental function solution(x, y, t) result(val)
    real(dp), intent(in) :: x, y, t
    real(dp)             :: val

    val = 0.5_dp
  end function solution

end module mod_usr
