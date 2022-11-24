module mod_usr
  use mod_rho
  use mod_multigrid_coupling

  implicit none

  integer             :: i_sol
  integer             :: i_eps1, i_eps2
  integer             :: i_err,i_rel_err
  real(dp), parameter :: sig1 = 1.0d-2
  real(dp), parameter :: sig2 = 1.0d-2
  real(dp), parameter :: diffusion_coeff1 = 1.d0
  real(dp), parameter :: diffusion_coeff2 = 16.d0

contains

  subroutine usr_init()
    use mod_usr_methods
    use mod_variables
    use mod_physics

    integer :: i

    usr_init_one_grid => initial_conditions

    use_multigrid = .true.
    phys_global_source_after => diffuse_density
    usr_process_grid => set_error
    mg_after_new_tree => set_epsilon

    mg%operator_type = mg_ahelmholtz
    mg%bc(:, mg_iphi)%bc_type = mg_bc_continuous

    do i = 1, 2 * ndim
      mg%bc(i, mg_iphi)%boundary_cond => multigrid_bc
    end do

    call set_coordinate_system("Cartesian_2D")
    call rho_activate()

    i_sol = var_set_extravar("sol", "sol")
    i_err = var_set_extravar("err", "err")
    i_rel_err = var_set_extravar("rel_err", "rel_err")
    i_eps1 = var_set_extravar("eps1", "eps1")
    i_eps2 = var_set_extravar("eps2", "eps2")

  end subroutine usr_init

  subroutine initial_conditions(ixG^L, ix^L, w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixG^L, ix^L
    double precision, intent(in)    :: x(ixG^S, 1:ndim)
    double precision, intent(inout) :: w(ixG^S, 1:nw)

    w(ix^S, rho_) = solution(x(ix^S, 1), x(ix^S, 2), 0.0d0)
    w(ix^S, i_eps1) = diffusion_coeff1
    w(ix^S, i_eps2) = diffusion_coeff2

  end subroutine initial_conditions

  !> To fill ghost cells near physical boundaries
  subroutine multigrid_bc(box, nc, iv, nb, bc_type, bc)
    use mod_bc_data
    type(mg_box_t), intent(in)    :: box
    integer, intent(in)           :: nc
    integer, intent(in)           :: iv      !< Index of variable
    integer, intent(in)           :: nb      !< Direction
    integer, intent(out)          :: bc_type !< Type of b.c.

    double precision, intent(out) :: bc(nc)
    double precision              :: rr(nc, 2)
    integer                       :: i, ixs(ndim-1), nb_dim, n_bc

    ! Type of boundary condition
    bc_type      = mg_bc_dirichlet
    ! Index of boundary data (stored per variable, per direction)
    n_bc         = bc_data_ix(rho_, nb)
    ! In which dimension the boundary is
    nb_dim       = mg_neighb_dim(nb)
    ! A list of dimensions *other* than nb_dim
    ixs          = [(i, i=1,ndim-1)]
    ixs(nb_dim:) = ixs(nb_dim:) + 1

    ! Get the coordinates at the cell faces at the boundary
    call mg_get_face_coords(box, nb, nc, rr)

    bc(:) = solution(rr(:,1),rr(:,2),global_time)

  end subroutine multigrid_bc


  elemental function solution(x, y, t) result(sol)
    real(dp), intent(in) :: x, y, t
    real(dp)             :: sol,k1,k2

    k1 = sig1+diffusion_coeff1*t
    k2 = sig2+diffusion_coeff2*t

    sol = dsqrt(1.d0/(4*dpi*k1))*dsqrt(1.d0/(4*dpi*k2))
    sol = sol*dexp(-x**2/(4*k1))*dexp(-y**2/(4*k2))

    sol = 1.d0 + sol

  end function solution

  subroutine diffuse_density(qdt, qt, active)
    use m_octree_mg_2d
    double precision, intent(in) :: qdt
    double precision, intent(in) :: qt
    logical, intent(inout)       :: active
    double precision             :: max_res

    call mg_copy_to_tree(rho_, mg_iphi, .false., .false.)
    call diffusion_solve_acoeff(mg, qdt, 1, 1d-4)
    call mg_copy_from_tree(mg_iphi, rho_)

    active = .true.
  end subroutine diffuse_density

  subroutine set_error(igrid,level,ixI^L,ixO^L,qt,w,x)
    integer, intent(in)             :: igrid,level,ixI^L,ixO^L
    double precision, intent(in)    :: qt,x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    w(ixO^S,i_sol) = solution(x(ixO^S, 1), x(ixO^S, 2), qt)
    w(ixO^S,i_err) = abs(w(ixO^S,rho_) - w(ixO^S,i_sol))
    w(ixO^S,i_rel_err) = abs(w(ixO^S,rho_) - w(ixO^S,i_sol))/w(ixO^S,i_sol)
  end subroutine set_error

  subroutine set_epsilon()
    call mg_copy_to_tree(i_eps1, mg_iveps1, .true., .true.)
    call mg_copy_to_tree(i_eps2, mg_iveps2, .true., .true.)
  end subroutine set_epsilon

end module mod_usr
