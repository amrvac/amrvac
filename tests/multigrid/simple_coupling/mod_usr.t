module mod_usr
  use mod_rho
  use mod_multigrid_coupling

  implicit none

  integer             :: i_rhs
  integer             :: i_phi
  integer             :: i_err
  integer             :: i_gradx
  integer             :: i_grady
  integer             :: i_err_gradx
  integer             :: i_err_grady
  real(dp), parameter :: pi = acos(-1.0_dp)

contains

  subroutine usr_init()
    use mod_variables

    usr_init_one_grid => initial_conditions

    use_multigrid = .true.
    usr_process_global => compute_phi
    usr_process_grid => set_error

    call set_coordinate_system("Cartesian_2D")
    call rho_activate()

    i_phi = var_set_extravar("phi", "phi")
    i_rhs = var_set_extravar("rhs", "rhs")
    i_err = var_set_extravar("err", "err")
    i_gradx = var_set_extravar("gradx", "gradx")
    i_grady = var_set_extravar("grady", "grady")
    i_err_gradx = var_set_extravar("err_gradx", "err_gradx")
    i_err_grady = var_set_extravar("err_grady", "err_grady")
  end subroutine usr_init

  subroutine initial_conditions(ixG^L, ix^L, w,x)
    integer, intent(in)             :: ixG^L, ix^L
    double precision, intent(in)    :: x(ixG^S, 1:ndim)
    double precision, intent(inout) :: w(ixG^S, 1:nw)

    w(ix^S, rho_) = sin(2 * pi * x(ix^S, 1)) * &
         sin(2 * pi * x(ix^S, 2))

    ! Analytic right-hand side
    w(ix^S, i_rhs) = -8 * pi**2 * &
         sin(2 * pi * x(ix^S, 1)) * &
         sin(2 * pi * x(ix^S, 2))
  end subroutine initial_conditions

  subroutine compute_phi(iit,qt)
    integer, intent(in)          :: iit
    double precision, intent(in) :: qt
    integer                      :: id
    double precision             :: max_res

    call mg_copy_to_tree(i_rhs, mg_irhs, .false., .false.)
    call mg_fas_fmg(mg, .true., max_res)
    if (mype == 0) print *, it, "max residu:", max_res
    call mg_copy_from_tree_gc(mg_iphi, i_phi)
  end subroutine compute_phi

  subroutine set_error(igrid,level,ixI^L,ixO^L,qt,w,x)
    integer, intent(in)             :: igrid,level,ixI^L,ixO^L
    double precision, intent(in)    :: qt,x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    w(ixO^S, i_err) = w(ixO^S, i_phi) - w(ixO^S, rho_)
    call gradient(w(ixI^S, i_phi), ixI^L, ixO^L, 1, w(ixI^S, i_gradx))
    call gradient(w(ixI^S, i_phi), ixI^L, ixO^L, 2, w(ixI^S, i_grady))

    w(ixO^S, i_err_gradx) = w(ixO^S, i_gradx) - &
         2 * pi * cos(2 * pi * x(ixO^S, 1)) * &
         sin(2 * pi * x(ixO^S, 2))

    w(ixO^S, i_err_grady) = w(ixO^S, i_grady) - &
         2 * pi * sin(2 * pi * x(ixO^S, 1)) * &
         cos(2 * pi * x(ixO^S, 2))
  end subroutine set_error

end module mod_usr

