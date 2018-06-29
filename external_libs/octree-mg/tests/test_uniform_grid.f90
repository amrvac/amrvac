#include "../src/cpp_macros.h"
program test_one_level
  use mpi
  use m_octree_mg

  implicit none

  integer             :: n_its          = 10
  integer             :: n_modes(NDIM)  = 5
  integer             :: box_size
  integer             :: domain_size(NDIM)
  real(dp)            :: dr(NDIM)
  real(dp)            :: r_min(NDIM)    = 0.0_dp
  logical             :: periodic(NDIM) = .false.
  integer             :: n_finer        = 0
  real(dp), parameter :: pi             = acos(-1.0_dp)
  character(len=40)   :: arg_string
  integer             :: n, ierr, n_args
  real(dp)            :: t0, t1
  integer             :: i_sol, i_eps
  type(mg_t)          :: mg

  n_args = command_argument_count()
  if (n_args /= NDIM+1) then
     error stop "Usage: ./test_uniform_grid box_size nx ny [nz]"
  end if

  call get_command_argument(1, arg_string)
  read(arg_string, *) box_size

  do n = 1, NDIM
     call get_command_argument(1+n, arg_string)
     read(arg_string, *) domain_size(n)
  end do

  dr =  1.0_dp / domain_size

  mg%n_extra_vars = 2
  i_eps = mg_num_vars + 1
  i_sol = mg_num_vars + 2

  mg%geometry_type = mg_cartesian
  mg%operator_type = mg_laplacian
  mg%smoother_type = mg_smoother_gs

  do n = 1, mg_num_neighbors
     mg%bc(n, mg_iphi)%bc_type = mg_bc_dirichlet
     mg%bc(n, mg_iphi)%bc_value = 0.0_dp
  end do

  call mg_set_methods(mg)
  call mg_comm_init(mg)
  call mg_build_rectangle(mg, domain_size, box_size, dr, r_min, &
       periodic, n_finer)
  call mg_load_balance(mg)

  call mg_allocate_storage(mg)
  call set_solution(mg)
  call compute_rhs_and_reset(mg)

  call print_error(mg, 0)

  t0 = mpi_wtime()
  do n = 1, 10
     call mg_fas_fmg(mg, n > 1)
     call print_error(mg, n)
  end do
  t1 = mpi_wtime()

  if (mg%my_rank == 0) then
     print *, "n_cpu            ", mg%n_cpu
     print *, "problem_size     ", domain_size
     print *, "box_size         ", box_size
     print *, "n_iterations     ", n_its
     print *, "time/iteration   ", (t1-t0) / n_its
     print *, "total_time(s)    ", (t1-t0)
     print *, "unknowns/microsec", 1e-6_dp * n_its * &
          product(real(domain_size, dp)) / (t1-t0)
     print *, ""
  end if
  call mg_timers_show(mg)

  call mpi_barrier(mpi_comm_world, ierr)
  call mpi_finalize(ierr)

contains

  subroutine set_solution(mg)
    type(mg_t), intent(inout) :: mg
    integer                   :: n, id, lvl, nc, IJK
    real(dp)                  :: r(NDIM), sol

    do lvl = mg%lowest_lvl, mg%highest_lvl
       nc = mg%box_size_lvl(lvl)
       do n = 1, size(mg%lvls(lvl)%my_ids)
          id = mg%lvls(lvl)%my_ids(n)
          do KJI_DO(0, nc+1)
             r = mg%boxes(id)%r_min + ([IJK] - 0.5_dp) * mg%dr(:, lvl)
             sol = product(sin(2 * pi * n_modes * r))
             mg%boxes(id)%cc(IJK, i_sol) = sol
             mg%boxes(id)%cc(IJK, i_eps) = max(1.0_dp, 1.0_dp + r(1))
          end do; CLOSE_DO
       end do
    end do
  end subroutine set_solution

  subroutine compute_rhs_and_reset(mg)
    type(mg_t), intent(inout) :: mg
    integer                   :: n, id, lvl, nc

    do lvl = mg%lowest_lvl, mg%highest_lvl
       nc = mg%box_size_lvl(lvl)
       do n = 1, size(mg%lvls(lvl)%my_ids)
          id = mg%lvls(lvl)%my_ids(n)

          mg%boxes(id)%cc(DTIMES(:), mg_iphi) = &
               mg%boxes(id)%cc(DTIMES(:), i_sol)
          call mg%box_op(mg, id, nc, mg_irhs)
          mg%boxes(id)%cc(DTIMES(:), mg_iphi) = 0.0_dp
       end do
    end do
  end subroutine compute_rhs_and_reset

  subroutine print_error(mg, it)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: it
    integer                   :: n, nc, id, lvl, IJK, ierr
    real(dp)                  :: sol, val, err, max_err

    err = 0.0_dp

    do lvl = mg%highest_lvl, mg%highest_lvl
       nc = mg%box_size_lvl(lvl)
       do n = 1, size(mg%lvls(lvl)%my_ids)
          id = mg%lvls(lvl)%my_ids(n)
          do KJI_DO(1, nc)
             sol = mg%boxes(id)%cc(IJK, i_sol)
             val = mg%boxes(id)%cc(IJK, mg_iphi)
             err = max(err, abs(val-sol))
             ! err = max(err, abs(mg%boxes(id)%cc(IJK, mg_ires)))
             ! print *, lvl, id, i, j, sol, val, abs(sol-val)
          end do; CLOSE_DO
       end do
    end do

    call mpi_reduce(err, max_err, 1, MPI_DOUBLE, MPI_MAX, 0, &
         mpi_comm_world, ierr)
    if (mg%my_rank == 0) print *, it, "max err", max_err
  end subroutine print_error

end program test_one_level
