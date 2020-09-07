module mod_usr
  use mod_rho
  use mod_multigrid_coupling

  implicit none

  integer             :: i_sol
  integer             :: i_err
  real(dp), parameter :: pi = acos(-1.0_dp)
  real(dp), parameter :: diffusion_coeff = 0.2_dp
  real(dp), parameter :: solution_modes(2) = [1, 1]

contains

  subroutine usr_init()
    use mod_usr_methods
    use mod_variables
    use mod_physics

    usr_init_one_grid => initial_conditions

    use_multigrid = .true.
    phys_global_source_after => diffuse_density
    usr_process_grid => set_error

    mg%bc(:, mg_iphi)%bc_type = mg_bc_neumann
    mg%bc(:, mg_iphi)%bc_value = 0.0d0

    call set_coordinate_system("Cartesian_2D")
    call rho_activate()

    i_sol = var_set_extravar("sol", "sol")
    i_err = var_set_extravar("err", "err")

  end subroutine usr_init

  subroutine initial_conditions(ixG^L, ix^L, w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixG^L, ix^L
    double precision, intent(in)    :: x(ixG^S, 1:ndim)
    double precision, intent(inout) :: w(ixG^S, 1:nw)

    w(ix^S, rho_) = solution(x(ix^S, 1), x(ix^S, 2), 0.0d0)

  end subroutine initial_conditions

  elemental function solution(x, y, t) result(sol)
    real(dp), intent(in) :: x, y, t
    real(dp)             :: sol, tmp(ndim)

    tmp = 2 * pi * solution_modes * [x, y]
    sol = 1 + product(cos(tmp)) * &
         exp(-sum((2 * pi * solution_modes)**2) * diffusion_coeff * t)
  end function solution

  subroutine diffuse_density(qdt, qt, active)
    double precision, intent(in) :: qdt
    double precision, intent(in) :: qt
    logical, intent(inout)       :: active
    double precision             :: max_res

    call mg_copy_to_tree(rho_, mg_iphi)
    call diffusion_solve(mg, qdt, diffusion_coeff, 1, 1d-4)
    call mg_copy_from_tree(mg_iphi, rho_)
    active = .true.
  end subroutine diffuse_density

  subroutine set_error(igrid,level,ixI^L,ixO^L,qt,w,x)
    integer, intent(in)             :: igrid,level,ixI^L,ixO^L
    double precision, intent(in)    :: qt,x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    w(ixO^S,i_sol) = solution(x(ixO^S, 1), x(ixO^S, 2), qt)
    w(ixO^S,i_err) = abs(w(ixO^S,rho_) - w(ixO^S,i_sol))
  end subroutine set_error

end module mod_usr

