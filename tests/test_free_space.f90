#include "../src/cpp_macros.h"
program test_performance
  use mpi
  use m_octree_mg
  use m_free_space

  implicit none

  real(dp), parameter :: gauss_ampl    = 1.0d0
  real(dp), parameter :: gauss_r0(3)   = [0.5d0, 0.5d0, 0.5d0]
  real(dp), parameter :: gauss_sigma   = 0.1d0
  real(dp), parameter :: pi            = acos(-1.0_dp)
  real(dp), parameter :: domain_len(3) = [1.0_dp, 1.0_dp, 1.0_dp]
  integer, parameter  :: n_its         = 5

  ! Fraction of unknowns the fft solver can use (compared to total number of
  ! multigrid unknowns)
  real(dp)            :: fft_frac = 0.15_dp

  integer             :: box_size
  integer             :: domain_size(NDIM)
  real(dp)            :: dr(NDIM), r_min(NDIM) = 0.0_dp
  logical             :: periodic(NDIM)        = .false.
  integer             :: n_finer               = 0
  character(len=40)   :: arg_string
  integer             :: n, ierr, n_args
  real(dp)            :: t0, t1, max_res
  integer             :: i_sol
  type(mg_t)          :: mg

  n_args = command_argument_count()
  if (n_args < NDIM+1) then
     error stop "Usage: ./test_uniform_grid box_size nx ny nz [fft_frac]"
  end if

  call get_command_argument(1, arg_string)
  read(arg_string, *) box_size

  do n = 1, NDIM
     call get_command_argument(1+n, arg_string)
     read(arg_string, *) domain_size(n)
  end do

  if (n_args > NDIM+1) then
     call get_command_argument(NDIM+2, arg_string)
     read(arg_string, *) fft_frac
  end if

  dr = domain_len / domain_size

  mg%n_extra_vars = 1
  i_sol = mg_num_vars + 1

  mg%geometry_type = mg_cartesian
  mg%operator_type = mg_laplacian
  mg%smoother_type = mg_smoother_gsrb

  call mg_set_methods(mg)
  call mg_comm_init(mg)
  call mg_build_rectangle(mg, domain_size, box_size, dr, r_min, &
       periodic, n_finer)
  call mg_load_balance(mg)
  call mg_allocate_storage(mg)
  call set_rhs_and_solution(mg)

  t0 = mpi_wtime()
  do n = 1, n_its
     call mg_poisson_free_3d(mg, n == 1, fft_frac, .true., max_res)
     call print_error(mg, n, max_res)
  end do
  t1 = mpi_wtime()

  if (mg%my_rank == 0) then
     print *, "fft_frac         ", fft_frac
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

  elemental function solution(x, y, z) result(val)
    real(dp), intent(in) :: x, y, z
    real(dp)             :: val, rnorm
    real(dp), parameter  :: fac = 1/(4 * pi)

    rnorm = norm2([ x, y, z ] - gauss_r0)
    if (rnorm < sqrt(epsilon(1.0d0))) then
       val = 2 * fac * gauss_ampl / (sqrt(pi) * gauss_sigma)
    else
       val = fac * gauss_ampl * erf(rnorm / gauss_sigma) / rnorm
    end if
  end function solution

  elemental function rhs(x, y, z) result(val)
    real(dp), intent(in) :: x, y, z
    real(dp)             :: val, r(NDIM)

    r = ([ x, y, z ] - gauss_r0) / gauss_sigma
    val = -gauss_ampl / (gauss_sigma**3 * pi * sqrt(pi)) * &
         exp(-sum(r**2))
  end function rhs

  subroutine set_rhs_and_solution(mg)
    type(mg_t), intent(inout) :: mg
    real(dp)                  :: r(3)
    integer                   :: n, id, lvl, nc, IJK

    do lvl = mg%lowest_lvl, mg%highest_lvl
       nc = mg%box_size_lvl(lvl)
       do n = 1, size(mg%lvls(lvl)%my_ids)
          id = mg%lvls(lvl)%my_ids(n)
          do KJI_DO(1, nc)
             r = mg%boxes(id)%r_min + ([IJK] - 0.5_dp) * mg%dr(:, lvl)
             mg%boxes(id)%cc(IJK, mg_irhs) = rhs(r(1), r(2), r(3))
             mg%boxes(id)%cc(IJK, i_sol) = solution(r(1), r(2), r(3))
          end do; CLOSE_DO
       end do
    end do
  end subroutine set_rhs_and_solution

  subroutine print_error(mg, it, max_res)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: it
    real(dp), intent(in)      :: max_res
    integer                   :: n, nc, id, lvl, IJK, ierr
    real(dp)                  :: sol, val, err, max_err, err2_sum

    max_err  = 0.0_dp
    err2_sum = 0.0_dp

    do lvl = mg%highest_lvl, mg%highest_lvl
       nc = mg%box_size_lvl(lvl)
       do n = 1, size(mg%lvls(lvl)%my_ids)
          id = mg%lvls(lvl)%my_ids(n)
          do KJI_DO(1, nc)
             sol      = mg%boxes(id)%cc(IJK, i_sol)
             val      = mg%boxes(id)%cc(IJK, mg_iphi)
             err      = abs(val-sol)
             max_err  = max(max_err, err)
             err2_sum = err2_sum + err**2
          end do; CLOSE_DO
       end do
    end do

    call mpi_allreduce(MPI_IN_PLACE, max_err, 1, MPI_DOUBLE, MPI_MAX, &
         mpi_comm_world, ierr)
    call mpi_allreduce(MPI_IN_PLACE, err2_sum, 1, MPI_DOUBLE, MPI_SUM, &
         mpi_comm_world, ierr)
    if (mg%my_rank == 0) print *, it, "max err/err2/res", &
         max_err, sqrt(err2_sum/mg_number_of_unknowns(mg)), max_res
  end subroutine print_error

end program test_performance
