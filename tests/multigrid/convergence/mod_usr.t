module mod_usr
  use mod_rho
  use mod_multigrid_coupling

  implicit none

  integer             :: i_rhs
  integer             :: i_phi
  integer             :: i_err
  real(dp), parameter :: pi = acos(-1.0_dp)

  double precision, parameter :: gauss_ampl      = 1.0d1
  double precision, parameter :: gauss_sigma     = 1.0d-1
  integer, parameter          :: cos_modes(ndim) = [1, 2, 3]

contains

  subroutine usr_init()
    use mod_variables
    integer :: n

    usr_init_one_grid => initial_conditions

    use_multigrid = .true.
    usr_process_global => compute_phi
    usr_process_grid => set_error
    usr_refine_grid => my_refine

    call set_coordinate_system("Cartesian_2D")
    call rho_activate()

    i_phi = var_set_extravar("phi", "phi")
    i_rhs = var_set_extravar("rhs", "rhs")
    i_err = var_set_extravar("err", "err")

    do n = 1, mg_num_neighbors
       mg%bc(n, mg_iphi)%boundary_cond => my_bc
    end do
  end subroutine usr_init

  subroutine initial_conditions(ixG^L, ix^L, w,x)
    integer, intent(in)             :: ixG^L, ix^L
    double precision, intent(in)    :: x(ixG^S, 1:ndim)
    double precision, intent(inout) :: w(ixG^S, 1:nw)

    ! Analytic solution
    w(ix^S, rho_) = solution(x(ix^S, 1), x(ix^S, 2), x(ix^S, 3))

    ! Analytic right-hand side
    w(ix^S, i_rhs) = rhs(x(ix^S, 1), x(ix^S, 2), x(ix^S, 3))
  end subroutine initial_conditions

  elemental function solution(x^D) result(val)
    double precision, intent(in) :: x^D
    double precision             :: val, r(ndim)

    r = [ x^D ] / gauss_sigma
    val = gauss_ampl * exp(-sum(r**2))

    r = [ x^D ] * 2 * pi * cos_modes
    val = val + cos(sum(r))
  end function solution

  elemental function rhs(x^D) result(val)
    double precision, intent(in) :: x^D
    double precision             :: val, r(ndim)

    r = [ x^D ] / gauss_sigma
    val = gauss_ampl * 4/gauss_sigma**2 * &
         (sum(r**2) - 0.5_dp * ndim) * exp(-sum(r**2))

    r = [ x^D ] * 2 * pi * cos_modes
    val = val - 4 * pi**2 * sum(cos_modes**2) * cos(sum(r))
  end function rhs

  !> To fill ghost cells near physical boundaries
  subroutine my_bc(box, nc, iv, nb, bc_type, bc)
    type(mg_box_t), intent(in)    :: box
    integer, intent(in)           :: nc
    integer, intent(in)           :: iv      !< Index of variable
    integer, intent(in)           :: nb      !< Direction
    integer, intent(out)          :: bc_type !< Type of b.c.
    double precision, intent(out) :: bc(nc, nc)
    double precision              :: x(nc, nc, 3)

    call mg_get_face_coords(box, nb, nc, x)
    bc_type = mg_bc_dirichlet
    bc      = solution(x(:, :, 1), x(:, :, 2), x(:, :, 3))
  end subroutine my_bc

  subroutine compute_phi(iit,qt)
    integer, intent(in)          :: iit
    double precision, intent(in) :: qt
    integer                      :: id
    double precision             :: max_res

    call mg_copy_to_tree(i_rhs, mg_irhs, .false., .false.)
    call mg_fas_fmg(mg, .true., max_res)
    ! call mg_fas_vcycle(mg, max_res=max_res)
    if (mype == 0) print *, it, "max residu:", max_res
    call mg_copy_from_tree(mg_iphi, i_phi)
  end subroutine compute_phi

  subroutine set_error(igrid,level,ixI^L,ixO^L,qt,w,x)
    integer, intent(in)             :: igrid,level,ixI^L,ixO^L
    double precision, intent(in)    :: qt,x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    w(ixO^S, i_err) = w(ixO^S, i_phi) - w(ixO^S, rho_)
  end subroutine set_error

  subroutine my_refine(igrid,level,ixI^L,ixO^L,qt,w,x,refine,coarsen)
    use mod_global_parameters
    integer, intent(in)          :: igrid, level, ixI^L, ixO^L
    double precision, intent(in) :: qt, w(ixI^S,1:nw), x(ixI^S,1:ndim)
    integer, intent(inout)       :: refine, coarsen

    refine = 0
    coarsen = 0
    if (all(x(ixO^S, 1) <= 0.1d0) .and. &
         all(x(ixO^S, 2) <= 0.1d0) .and. &
         all(x(ixO^S, 3) <= 0.1d0)) then
       refine = 1
       coarsen = -1
    else
       refine = -1
       coarsen = -1
    end if
  end subroutine my_refine

end module mod_usr

