module mod_usr
  use mod_rho
  use mod_multigrid_coupling

  implicit none

  integer             :: my_rhs
  integer             :: my_eps
  integer             :: my_phi
  real(dp), parameter :: pi = acos(-1.0_dp)

contains

  subroutine usr_init()
    use mod_usr_methods
    use mod_variables

    usr_init_one_grid => initial_conditions

    use_multigrid = .true.
    usr_process_global => compute_phi

    mg%operator_type = mg_vlaplacian
    mg_after_new_tree => set_epsilon

    call set_coordinate_system("Cartesian_2D")
    call rho_activate()

    my_phi = var_set_extravar("phi", "phi")
    my_rhs = var_set_extravar("rhs", "rhs")
    my_eps = var_set_extravar("eps", "eps")

  end subroutine usr_init

  subroutine initial_conditions(ixG^L, ix^L, w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixG^L, ix^L
    double precision, intent(in)    :: x(ixG^S, 1:ndim)
    double precision, intent(inout) :: w(ixG^S, 1:nw)

    w(ix^S, rho_) = sin(2 * pi * x(ix^S, 1)) * &
         sin(2 * pi * x(ix^S, 2))

    ! Analytic right-hand side
    w(ix^S, my_rhs) = -8 * pi**2 * &
         sin(2 * pi * x(ix^S, 1)) * &
         sin(2 * pi * x(ix^S, 2))

    ! Set variable coefficient (should be smooth)
    ! w(ix^S, my_eps) = 2 + &
    !      cos(2 * pi * x(ix^S, 1)) * &
    !      cos(2 * pi * x(ix^S, 2))
    w(ix^S, my_eps) = exp(10*x(ix^S, 1))
  end subroutine initial_conditions

  subroutine set_epsilon()
    call mg_copy_to_tree(my_eps, mg_iveps, .true., .true.)
  end subroutine set_epsilon

  subroutine compute_phi(iit,qt)
    integer, intent(in)          :: iit
    double precision, intent(in) :: qt
    integer                      :: id
    double precision             :: max_res

    call mg_copy_to_tree(my_rhs, mg_irhs)
    call mg_fas_fmg(mg, .true., max_res)
    if (mype == 0) print *, it, "max residu:", max_res
    call mg_copy_from_tree(mg_iphi, my_phi)
  end subroutine compute_phi

end module mod_usr

