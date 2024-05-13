! Toy code for solving a 2D Laplace equation on a uniform grid with OpenACC. The
! grid is divided into blocks with a layer of ghost cells around them.
!
! Author: Jannis Teunissen (jannis.teunissen@cwi.nl)
program ghostcell
  implicit none
  integer, parameter :: dp = kind(0.0d0)
  integer, parameter :: NDIM = 2

  integer :: grid_size(NDIM)
  integer :: block_size(NDIM)
  integer :: blocks_per_dim(NDIM)
  integer :: n_blocks
  integer :: n_variables
  integer :: n_iterations
  integer :: n_gc
  logical :: write_output

  integer :: n_args
  character(len=100) :: tmp

  integer, allocatable :: grid_id(:, :)
  integer, allocatable :: id_to_index(:, :)
  integer, allocatable :: internal_boundary_x(:, :)
  integer, allocatable :: internal_boundary_y(:, :)

  integer, allocatable :: bc_xlo(:)
  integer, allocatable :: bc_xhi(:)
  integer, allocatable :: bc_ylo(:)
  integer, allocatable :: bc_yhi(:)

  integer :: n_bc(2*NDIM)
  integer :: i, j, n, m
  integer :: i_new, i_old
  integer :: n_boundaries_x, n_boundaries_y

  integer :: t_start, t_end, count_rate

  real(dp) :: t_total
  real(dp), allocatable :: uu(:, :, :, :)

  ! Set default values
  grid_size = 512
  block_size = 16
  n_gc = 1                      ! number of ghost cells
  n_iterations = 100
  write_output = .false.
  n_variables = 2

  ! Parse command line arguments
  n_args = command_argument_count()

  if (n_args == 0) then
     print *, "Usage: ./laplace_openacc_ghostcell Nx Ny block_Nx block_Ny"
     print *, "...    n_ghostcell n_iterations write_output"
  end if

  if (n_args >= 2) then
     call get_command_argument(1, tmp)
     read(tmp, *) grid_size(1)
     call get_command_argument(2, tmp)
     read(tmp, *) grid_size(2)
  end if

  if (n_args >= 4) then
     call get_command_argument(3, tmp)
     read(tmp, *) block_size(1)
     call get_command_argument(4, tmp)
     read(tmp, *) block_size(2)
  end if

  if (n_args >= 5) then
     call get_command_argument(5, tmp)
     read(tmp, *) n_gc
  end if

  if (n_args >= 6) then
     call get_command_argument(6, tmp)
     read(tmp, *) n_iterations
  end if

  if (n_args >= 7) then
     call get_command_argument(7, tmp)
     read(tmp, *) write_output
  end if

  if (n_args >= 8) then
     error stop "Too many arguments"
  end if

  blocks_per_dim = grid_size / block_size

  if (any(blocks_per_dim * block_size /= grid_size)) &
       error stop "Block size should divide grid size"

  n_blocks = product(blocks_per_dim)

  print *, "grid size", grid_size
  print *, "block size", block_size
  print *, "n_blocks", n_blocks
  print *, "n_iterations", n_iterations

  ! Identify each grid block with a number (id) for neighbor detection
  allocate(grid_id(0:blocks_per_dim(1)+1, 0:blocks_per_dim(2)+1))
  allocate(id_to_index(2, n_blocks))

  ! Negative values are used to identify boundaries
  grid_id = -1

  n = 0
  do j = 1, blocks_per_dim(2)
     do i = 1, blocks_per_dim(1)
        n = n + 1
        grid_id(i, j) = n
        id_to_index(:, n) = [i, j]
     end do
  end do

  ! Array of grid blocks
  allocate(uu(1-n_gc:block_size(1)+n_gc, 1-n_gc:block_size(2)+n_gc, &
       n_blocks, n_variables))

  ! Construct list of internal boundaries
  n_boundaries_x = count(grid_id(0:blocks_per_dim(1), :) > 0 .and. &
       grid_id(1:blocks_per_dim(1)+1, :) > 0)
  n_boundaries_y = count(grid_id(:, 0:blocks_per_dim(1)) > 0 .and. &
       grid_id(:, 1:blocks_per_dim(1)+1) > 0)
  allocate(internal_boundary_x(2, n_boundaries_x))
  allocate(internal_boundary_y(2, n_boundaries_y))

  ! Construct list of physical boundaries
  n_bc(1) = count(grid_id(0, 1:blocks_per_dim(2)) < 0)
  n_bc(2) = count(grid_id(blocks_per_dim(1)+1, 1:blocks_per_dim(2)) < 0)
  n_bc(3) = count(grid_id(1:blocks_per_dim(1), 0) < 0)
  n_bc(4) = count(grid_id(1:blocks_per_dim(1), blocks_per_dim(2)+1) < 0)
  allocate(bc_xlo(n_bc(1)), bc_xhi(n_bc(2)))
  allocate(bc_ylo(n_bc(3)), bc_yhi(n_bc(4)))

  n = 0
  m = 0
  do j = 1, blocks_per_dim(2)
     do i = 1, blocks_per_dim(1)
        if (grid_id(i, j) > 0 .and. grid_id(i+1, j) > 0) then
           n = n + 1
           internal_boundary_x(:, n) = [grid_id(i, j), grid_id(i+1, j)]
        end if

        if (grid_id(i, j) > 0 .and. grid_id(i, j+1) > 0) then
           m = m + 1
           internal_boundary_y(:, m) = [grid_id(i, j), grid_id(i, j+1)]
        end if
     end do
  end do

  call set_initial_conditions()

  if (write_output) call write_brick_of_values("test_", 0, 0.0_dp, 1)

  !$acc enter data copyin(uu)

  call system_clock(t_start, count_rate)
  do n = 1, n_iterations
     i_new = mod(n, 2) + 1
     i_old = mod(n+1, 2) + 1
     call update_ghostcells(uu(:, :, :, i_old))
     call update_solution(i_new, i_old)
  end do
  call system_clock(t_end, count_rate)

  !$acc exit data copyout(uu)

  if (write_output) call write_brick_of_values("test_", 1, real(n, dp), i_new)

  t_total = (t_end - t_start)/real(count_rate, dp)
  print *, "Total time: ", t_total
  write(*, "(A,E12.4)") "Unknowns/second: ", &
       n_iterations/t_total * product(grid_size)

contains

  subroutine set_initial_conditions()
    integer  :: i, j, n, id_a, id_b, i_cell
    real(dp) :: bc_left = 1.0_dp, bc_right = 0.0_dp
    real(dp) :: bc_value

    ! Zero everywhere
    uu(:, :, :, :) = 0

    ! Boundary on the left
    do j = 1, blocks_per_dim(2)
       id_a = grid_id(1, j)
       uu(1-n_gc:0, :, id_a, :) = bc_left
    end do

    ! Boundary on the right
    do j = 1, blocks_per_dim(2)
       id_a = grid_id(blocks_per_dim(1), j)
       uu(block_size(1)+1:block_size(1)+n_gc, :, id_a, :) = bc_right
    end do

    ! Linear profile on top and bottom boundary
    do i = 1, blocks_per_dim(1)
       id_a = grid_id(i, 1)
       id_b = grid_id(i, blocks_per_dim(2))

       ! No need for diagonal ghost cells
       do n = 1, block_size(1)
          i_cell = (i-1) * block_size(1) + n
          ! We should have bc_left at i_cell=0, and bc_right at
          ! i_cell=grid_size(1)+1
          bc_value = bc_left + i_cell/(grid_size(1) + 1.0_dp) * (bc_right - bc_left)

          uu(n, 1-n_gc:0, id_a, :) = bc_value
          uu(n, block_size(2)+1:block_size(2)+n_gc, id_b, :) = bc_value
       end do
    end do
  end subroutine set_initial_conditions

  subroutine update_solution(i_new, i_old)
    integer, intent(in) :: i_new, i_old
    integer             :: n, i, j

    !$acc parallel loop
    do n = 1, n_blocks
       !$acc loop
       do j = 1, block_size(2)
          !$acc loop
          do i = 1, block_size(1)
             uu(i, j, n, i_new) = 0.25_dp * ( &
                  uu(i+1, j, n, i_old) + uu(i-1, j, n, i_old) + &
                  uu(i, j-1, n, i_old) + uu(i, j+1, n, i_old))
          end do
       end do
    end do
  end subroutine update_solution

  subroutine update_ghostcells(uus)
    real(dp), intent(inout) :: uus(1-n_gc:block_size(1)+n_gc, &
         1-n_gc:block_size(2)+n_gc, n_blocks)
    integer             :: n, i, j

    !$acc parallel loop async(1)
    do n = 1, n_boundaries_x
       !$acc loop collapse(2)
       do j = 1, block_size(2)
          do i = 1, n_gc
             uus(block_size(1)+i, j, internal_boundary_x(1, n)) = &
                  uus(i, j, internal_boundary_x(2, n))
             uus(-n_gc+i, j, internal_boundary_x(2, n)) = &
                  uus(block_size(1)-n_gc+i, j, internal_boundary_x(1, n))
          end do
       end do
    end do

    !$acc parallel loop async(2)
    do n = 1, n_boundaries_y
       !$acc loop collapse(2)
       do j = 1, n_gc
          do i = 1, block_size(1)
             uus(i, block_size(2)+j, internal_boundary_y(1, n)) = &
                  uus(i, j, internal_boundary_y(2, n))
             uus(i, -n_gc+j, internal_boundary_y(2, n)) = &
                  uus(i, block_size(2)-n_gc+j, internal_boundary_y(1, n))
          end do
       end do
    end do

    !$acc wait

  end subroutine update_ghostcells

  ! Write data in brick-of-values (BOV) format, which is supported by Visit
  subroutine write_brick_of_values(basename, n, time, state)
    character(len=*), intent(in) :: basename
    integer, intent(in)          :: n
    real(dp), intent(in)         :: time
    integer, intent(in)          :: state
    character(len=200)           :: header_name, datfile_name
    integer                      :: my_unit, i, j, id, lo(2), hi(2)
    real(dp), allocatable        :: u_full(:, :)

    write(header_name, "(A,I6.6,A)") basename, n, ".bov"
    write(datfile_name, "(A,I6.6,A)") basename, n, ".dat"

    open(newunit=my_unit, file=trim(header_name), action="write")
    write(my_unit, *) "TIME: ", time
    write(my_unit, *) "DATA_FILE: ", trim(datfile_name)
    write(my_unit, *) "DATA_SIZE: ", grid_size, 1
    if (dp == kind(0.0d0)) then
       write(my_unit, *) "DATA_FORMAT: DOUBLE"
    else if (dp == kind(0.0)) then
       write(my_unit, *) "DATA_FORMAT: FLOAT"
    else
       error stop "Unknown data format"
    end if
    write(my_unit, *) "VARIABLE: phi"
    write(my_unit, *) "DATA_ENDIAN: LITTLE"
    write(my_unit, *) "CENTERING: zonal"
    write(my_unit, *) "BRICK_ORIGIN: 0. 0. 0."
    write(my_unit, *) "BRICK_SIZE: 1. 1. 1."
    write(my_unit, *) "DATA_COMPONENTS: 1"
    close(my_unit)

    allocate(u_full(grid_size(1), grid_size(2)))

    do j = 1, blocks_per_dim(2)
       do i = 1, blocks_per_dim(1)
          id = grid_id(i, j)
          lo = [i-1, j-1] * block_size + 1
          hi = lo + block_size - 1
          u_full(lo(1):hi(1), lo(2):hi(2)) = &
               uu(1:block_size(1), 1:block_size(2), id, state)
       end do
    end do

    open(newunit=my_unit, file=trim(datfile_name), form='unformatted', &
         access='stream', status='replace')
    write(my_unit) u_full
    close(my_unit)
    print *, "Written ", trim(header_name)

  end subroutine write_brick_of_values

end program ghostcell
