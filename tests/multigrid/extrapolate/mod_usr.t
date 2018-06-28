module mod_usr
  use mod_rho
  use mod_multigrid_coupling

  implicit none

  integer :: i_phi
  integer :: mg_boundary_

contains

  subroutine usr_init()

    usr_init_one_grid => initial_conditions
    usr_process_global => compute_phi
    usr_set_parameters => set_boundary_conds

    use_multigrid = .true.

    call set_coordinate_system("Cartesian")
    call rho_activate()

    i_phi = var_set_extravar("phi", "phi")
  end subroutine usr_init

  subroutine set_boundary_conds()
    use mod_bc_data
    integer :: i

    mg%bc(:, mg_iphi)%bc_type = mg_bc_dirichlet

    do i = 1, 2 * ndim
       if (bc_data_ix(rho_, i) /= -1) then
          mg%bc(i, mg_iphi)%boundary_cond => multigrid_bc
       end if
    end do
  end subroutine set_boundary_conds

  subroutine initial_conditions(ixG^L, ix^L, w,x)
    integer, intent(in)             :: ixG^L, ix^L
    double precision, intent(in)    :: x(ixG^S, 1:ndim)
    double precision, intent(inout) :: w(ixG^S, 1:nw)

    w(ix^S, rho_) = 0.0d0
  end subroutine initial_conditions

  !> To fill ghost cells near physical boundaries
  subroutine multigrid_bc(box, nc, iv, nb, bc_type, bc)
    use mod_bc_data
    type(mg_box_t), intent(in)    :: box
    integer, intent(in)           :: nc
    integer, intent(in)           :: iv      !< Index of variable
    integer, intent(in)           :: nb      !< Direction
    integer, intent(out)          :: bc_type !< Type of b.c.
    {^IFTHREED
    double precision, intent(out) :: bc(nc, nc)
    double precision              :: rr(nc, nc, 3)
    }
    {^IFTWOD
    double precision, intent(out) :: bc(nc)
    double precision              :: rr(nc, 2)
    }
    integer                       :: i, ixs(ndim-1), nb_dim, n_bc

    bc_type      = mg_bc_dirichlet
    n_bc         = bc_data_ix(rho_, nb)
    nb_dim       = mg_neighb_dim(nb)
    ixs          = [(i, i=1,ndim-1)]
    ixs(nb_dim:) = ixs(nb_dim:) + 1

    call mg_get_face_coords(box, nb, nc, rr)

    {^IFTWOD
    bc = bc_data_get_2d(n_bc, rr(:, ixs(1)), global_time)
    }
    {^IFTHREED
    bc = bc_data_get_3d(n_bc, rr(:, :, ixs(1)), rr(:, :, ixs(2)), global_time)
    }
  end subroutine multigrid_bc

  subroutine compute_phi(iit,qt)
    integer, intent(in)          :: iit
    double precision, intent(in) :: qt
    integer                      :: id
    double precision             :: max_res

    if (iit == 0) then
       call mg_phi_bc_store(mg)
    end if

    ! call mg_copy_to_tree(i_rhs, mg_irhs, .false., .false.)
    call mg_fas_fmg(mg, .true., max_res)
    if (mype == 0) print *, it, "max residu:", max_res
    call mg_copy_from_tree_gc(mg_iphi, i_phi)
  end subroutine compute_phi

end module mod_usr

