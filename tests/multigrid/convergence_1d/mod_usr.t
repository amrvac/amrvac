module mod_usr
  use mod_rho
  use mod_multigrid_coupling

  implicit none

  integer             :: i_rhs
  integer             :: i_phi
  integer             :: i_err
  integer             :: i_res
  real(dp), parameter :: pi = acos(-1.0_dp)

  double precision, parameter :: gauss_ampl           = 10.0d0
  double precision, parameter :: gauss_sigma          = 0.1d0
  double precision, parameter :: cos_ampl             = 1.0d0
  double precision, parameter :: cos_modes(ndim)      = [0.5d0]
  double precision, parameter :: gauss_r0(ndim)       = [0.0d0]
  character(len=10)           :: multigrid_cycle_name = 'fmg'
  character(len=10)           :: refine_type = 'center'

contains

  subroutine usr_init()
    use mod_variables
    integer :: n

    usr_init_one_grid => initial_conditions

    use_multigrid = .true.
    usr_process_global => compute_phi
    usr_refine_grid => my_refine

    call set_coordinate_system("Cartesian_1D")
    call rho_activate()

    i_phi = var_set_extravar("phi", "phi")
    i_rhs = var_set_extravar("rhs", "rhs")
    i_err = var_set_extravar("err", "err")
    i_res = var_set_extravar("res", "res")

    do n = 1, mg_num_neighbors
       mg%bc(n, mg_iphi)%boundary_cond => my_bc
    end do

    mg%smoother_type = mg_smoother_gsrb
    mg%n_cycle_down = 2
    mg%n_cycle_up = 2
  end subroutine usr_init

  subroutine initial_conditions(ixG^L, ix^L, w,x)
    integer, intent(in)             :: ixG^L, ix^L
    double precision, intent(in)    :: x(ixG^S, 1:ndim)
    double precision, intent(inout) :: w(ixG^S, 1:nw)

    ! Analytic solution
    w(ix^S, rho_) = solution(x(ix^S, 1))

    ! Analytic right-hand side
    w(ix^S, i_rhs) = rhs(x(ix^S, 1))
  end subroutine initial_conditions

  elemental function solution(x^D) result(val)
    double precision, intent(in) :: x^D
    double precision             :: val, r(ndim)

    r = ([ x^D ] - gauss_r0) / gauss_sigma
    val = gauss_ampl * exp(-sum(r**2))

    r = [ x^D ] * 2 * pi * cos_modes
    val = val + cos_ampl * cos(sum(r))
  end function solution

  elemental function rhs(x^D) result(val)
    double precision, intent(in) :: x^D
    double precision             :: val, r(ndim)

    r = ([ x^D ] - gauss_r0) / gauss_sigma
    val = gauss_ampl * 4/gauss_sigma**2 * &
         (sum(r**2) - 0.5_dp * ndim) * exp(-sum(r**2))

    r = [ x^D ] * 2 * pi * cos_modes
    val = val - 4 * cos_ampl * pi**2 * sum(cos_modes**2) * cos(sum(r))
  end function rhs

  !> To fill ghost cells near physical boundaries
  subroutine my_bc(box, nc, iv, nb, bc_type, bc)
    type(mg_box_t), intent(in)    :: box
    integer, intent(in)           :: nc
    integer, intent(in)           :: iv      !< Index of variable
    integer, intent(in)           :: nb      !< Direction
    integer, intent(out)          :: bc_type !< Type of b.c.
    double precision, intent(out) :: bc(1)
    double precision              :: x(nc, 1)

    call mg_get_face_coords(box, nb, nc, x)
    bc_type = mg_bc_dirichlet
    bc      = solution(x(:, 1))
  end subroutine my_bc

  subroutine compute_phi(iit,qt)
    use mod_input_output
    integer, intent(in)          :: iit
    double precision, intent(in) :: qt
    integer                      :: id, iigrid, igrid
    double precision             :: max_res, wmax(nw), modes(nw), vol

    if (iit == 0) then
       call mg_copy_to_tree(i_rhs, mg_irhs, .false., .false.)
       call mg_phi_bc_store(mg)
    end if

    if (multigrid_cycle_name == 'fmg') then
       call mg_fas_fmg(mg, .true., max_res)
    else
       call mg_fas_vcycle(mg, max_res=max_res)
    end if

    call mg_copy_from_tree(mg_iphi, i_phi)
    call mg_copy_from_tree(mg_ires, i_res)

    do iigrid=1,igridstail; igrid=igrids(iigrid);
       call set_error(ixG^LL,ixM^LL,ps(igrid)%w,ps(igrid)%x)
    end do

    call get_global_maxima(wmax)
    call get_volume_average(2, modes, vol)
    if (mype == 0) print *, "CONV", it, max_res, wmax(i_err), sqrt(modes(i_err)), sqrt(modes(i_res))
  end subroutine compute_phi

  subroutine set_error(ixI^L,ixO^L,w,x)
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)

    w(ixO^S, i_err) = abs(w(ixO^S, i_phi) - w(ixO^S, rho_))
  end subroutine set_error

  subroutine my_refine(igrid,level,ixI^L,ixO^L,qt,w,x,refine,coarsen)
    use mod_global_parameters
    integer, intent(in)          :: igrid, level, ixI^L, ixO^L
    double precision, intent(in) :: qt, w(ixI^S,1:nw), x(ixI^S,1:ndim)
    integer, intent(inout)       :: refine, coarsen

    refine = 0
    coarsen = 0

    select case (refine_type)
    case ('corner')
       if (any(abs(x(ixO^S, 1) + 0.25d0) <= 0.1d0)) then
          refine = 1
          coarsen = -1
       else
          refine = -1
          coarsen = -1
       end if
    case ('center')
       if (any(abs(x(ixO^S, 1)) <= 0.1d0)) then
          refine = 1
          coarsen = -1
       else
          refine = -1
          coarsen = -1
       end if
    end select
  end subroutine my_refine

end module mod_usr

