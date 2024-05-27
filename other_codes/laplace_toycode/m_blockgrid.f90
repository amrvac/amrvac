#include "macros.h"
! Module for subdividing a grid into blocks
module m_blockgrid

  implicit none
  public

  integer, parameter, private :: dp = kind(0.0d0)

  type block_grid_t
     integer :: nx(2)
     integer :: bx(2)
     integer :: blocks_per_dim(2)
     integer :: n_gc
     integer :: n_blocks
     integer :: n_vars
     integer :: ilo(2), ihi(2)

     integer, allocatable :: id_to_index(:, :)

     integer, allocatable :: bnd_x(:, :)
     integer, allocatable :: bnd_y(:, :)
     integer, allocatable :: bc_xlo(:)
     integer, allocatable :: bc_xhi(:)
     integer, allocatable :: bc_ylo(:)
     integer, allocatable :: bc_yhi(:)

     real(dp), allocatable :: uu(:, :, :, :)
     real(dp)              :: bc(4)
  end type block_grid_t

contains

  subroutine initialize_grid(nx, bx, n_gc, n_vars, periodic, bg)
    integer, intent(in)             :: nx(2)
    integer, intent(in)             :: bx(2)
    integer, intent(in)             :: n_gc
    integer, intent(in)             :: n_vars
    logical, intent(in)             :: periodic(2)
    type(block_grid_t), intent(out) :: bg
    integer                         :: i, j, n, id
    integer                         :: n_bnd(2)
    integer                         :: n_bc(4)
    integer, allocatable            :: grid_id(:, :)

    bg%nx             = nx
    bg%bx             = bx
    bg%n_gc           = n_gc
    bg%blocks_per_dim = nx / bx

    if (any(bg%blocks_per_dim * bx /= nx)) &
         error stop "Block size should divide grid size"

    bg%n_blocks = product(bg%blocks_per_dim)
    bg%n_vars = n_vars

    ! Identify each grid block with a number (id) for neighbor detection
    allocate(grid_id(0:bg%blocks_per_dim(1)+1, 0:bg%blocks_per_dim(2)+1))
    allocate(bg%id_to_index(2, bg%n_blocks))

    ! Negative values are used to identify boundaries
    grid_id = -1

    n = 0
    do j = 1, bg%blocks_per_dim(2)
       do i = 1, bg%blocks_per_dim(1)
          n = n + 1
          grid_id(i, j) = n
          bg%id_to_index(:, n) = [i, j]
       end do
    end do

    if (periodic(1)) then
       grid_id(0, :) = grid_id(bg%blocks_per_dim(1), :)
       grid_id(bg%blocks_per_dim(1)+1, :) = grid_id(1, :)
    end if

    if (periodic(2)) then
       grid_id(:, 0) = grid_id(:, bg%blocks_per_dim(2))
       grid_id(:, bg%blocks_per_dim(2)+1) = grid_id(:, 1)
    end if

    ! Array of grid blocks
    allocate(bg%uu(IX(1-n_gc:bx(1)+n_gc, 1-n_gc:bx(2)+n_gc, bg%n_blocks, bg%n_vars)))
    bg%ilo = 1 - n_gc
    bg%ihi = bg%bx + n_gc

    ! Construct list of internal boundaries
    n_bnd(1) = count(grid_id(0:bg%blocks_per_dim(1), :) > 0 .and. &
         grid_id(1:bg%blocks_per_dim(1)+1, :) > 0)
    n_bnd(2) = count(grid_id(:, 0:bg%blocks_per_dim(1)) > 0 .and. &
         grid_id(:, 1:bg%blocks_per_dim(1)+1) > 0)
    allocate(bg%bnd_x(2, n_bnd(1)))
    allocate(bg%bnd_y(2, n_bnd(2)))

    ! Construct list of physical boundaries
    n_bc(1) = count(grid_id(0, 1:bg%blocks_per_dim(2)) < 0)
    n_bc(2) = count(grid_id(bg%blocks_per_dim(1)+1, 1:bg%blocks_per_dim(2)) < 0)
    n_bc(3) = count(grid_id(1:bg%blocks_per_dim(1), 0) < 0)
    n_bc(4) = count(grid_id(1:bg%blocks_per_dim(1), bg%blocks_per_dim(2)+1) < 0)
    allocate(bg%bc_xlo(n_bc(1)))
    allocate(bg%bc_xhi(n_bc(2)))
    allocate(bg%bc_ylo(n_bc(3)))
    allocate(bg%bc_yhi(n_bc(4)))

    n_bnd(:) = 0
    n_bc(:) = 0

    ! Store internal and physical boundaries
    do j = 1, bg%blocks_per_dim(2)
       do i = 1, bg%blocks_per_dim(1)
          id = grid_id(i, j)
          if (id > 0 .and. grid_id(i+1, j) > 0) then
             n_bnd(1) = n_bnd(1) + 1
             bg%bnd_x(:, n_bnd(1)) = [id, grid_id(i+1, j)]
          end if

          if (id > 0 .and. grid_id(i, j+1) > 0) then
             n_bnd(2) = n_bnd(2) + 1
             bg%bnd_y(:, n_bnd(2)) = [id, grid_id(i, j+1)]
          end if

          if (id > 0 .and. grid_id(i-1, j) <= 0) then
             n_bc(1) = n_bc(1) + 1
             bg%bc_xlo(n_bc(1)) = id
          end if

          if (id > 0 .and. grid_id(i+1, j) <= 0) then
             n_bc(2) = n_bc(2) + 1
             bg%bc_xhi(n_bc(2)) = id
          end if

          if (id > 0 .and. grid_id(i, j-1) <= 0) then
             n_bc(3) = n_bc(3) + 1
             bg%bc_ylo(n_bc(3)) = id
          end if

          if (id > 0 .and. grid_id(i, j+1) <= 0) then
             n_bc(4) = n_bc(4) + 1
             bg%bc_yhi(n_bc(4)) = id
          end if
       end do
    end do

    ! Copy data structure and allocatable components to device
    !$acc enter data copyin(bg)
    !$acc enter data copyin(bg%uu, bg%bnd_x, bg%bnd_y, bg%bc_xlo, bg%bc_xhi, &
    !$acc &bg%bc_ylo, bg%bc_yhi)
  end subroutine initialize_grid

  subroutine update_ghostcells(bg, ivar)
    type(block_grid_t), intent(inout) :: bg
    integer, intent(in)               :: ivar

    call update_ghostcells_args(bg%bx, bg%n_gc, bg%n_vars, bg%n_blocks, ivar, &
         size(bg%bnd_x, 2), bg%bnd_x, size(bg%bnd_y, 2), bg%bnd_y, bg%uu, bg)
  end subroutine update_ghostcells

  subroutine update_ghostcells_args(bx, n_gc, n_vars, n_blocks, ivar, n_bnd_x, bnd_x, &
       n_bnd_y, bnd_y, uu, bg)
    integer, intent(in) :: n_blocks, bx(2), n_gc, n_vars, ivar
    integer, intent(in) :: n_bnd_x, n_bnd_y
    integer, intent(in) :: bnd_x(2, n_bnd_x), bnd_y(2, n_bnd_y)
    real(dp), intent(inout) :: uu(IX(1-n_gc:bx(1)+n_gc, 1-n_gc:bx(2)+n_gc, n_blocks, n_vars))
    type(block_grid_t), intent(inout) :: bg
    integer                           :: n, i, j

    !$acc parallel loop async(1)
    do n = 1, n_bnd_x
       !$acc loop collapse(2)
       do j = 1, bx(2)
          do i = 1, n_gc
             uu(IX(bx(1)+i, j, bnd_x(1, n), ivar)) = &
                  uu(IX(i, j, bnd_x(2, n), ivar))
             uu(IX(-n_gc+i, j, bnd_x(2, n), ivar)) = &
                  uu(IX(bx(1)-n_gc+i, j, bnd_x(1, n), ivar))
          end do
       end do
    end do

    !$acc parallel loop async(2)
    do n = 1, n_bnd_y
       !$acc loop collapse(2)
       do j = 1, n_gc
          do i = 1, bx(1)
             uu(IX(i, bx(2)+j, bnd_y(1, n), ivar)) = &
                  uu(IX(i, j, bnd_y(2, n), ivar))
             uu(IX(i, -n_gc+j, bnd_y(2, n), ivar)) = &
                  uu(IX(i, bx(2)-n_gc+j, bnd_y(1, n), ivar))
          end do
       end do
    end do

    !$acc parallel loop async(3)
    do n = 1, size(bg%bc_xlo)
       !$acc loop collapse(2)
       do j = 1, bx(2)
          do i = 1, n_gc
             uu(IX(-n_gc+i, j, bg%bc_xlo(n), ivar)) = bg%bc(1)
          end do
       end do
    end do

    !$acc parallel loop async(4)
    do n = 1, size(bg%bc_xhi)
       !$acc loop collapse(2)
       do j = 1, bx(2)
          do i = 1, n_gc
             uu(IX(bx(1)+i, j, bg%bc_xhi(n), ivar)) = bg%bc(2)
          end do
       end do
    end do

    !$acc parallel loop async(5)
    do n = 1, size(bg%bc_ylo)
       !$acc loop collapse(2)
       do j = 1, n_gc
          do i = 1, bx(1)
             uu(IX(i, -n_gc+j, bg%bc_ylo(n), ivar)) = bg%bc(3)
          end do
       end do
    end do

    !$acc parallel loop async(6)
    do n = 1, size(bg%bc_yhi)
       !$acc loop collapse(2)
       do j = 1, n_gc
          do i = 1, bx(1)
             uu(IX(i, bx(2)+j, bg%bc_yhi(n), ivar)) = bg%bc(4)
          end do
       end do
    end do

    !$acc wait

  end subroutine update_ghostcells_args

  ! Write data in brick-of-values (BOV) format, which is supported by Visit
  subroutine write_brick_of_values(bg, basename, varname, n, time, ivar)
    type(block_grid_t), intent(in) :: bg
    character(len=*), intent(in)   :: basename
    character(len=*), intent(in)   :: varname
    integer, intent(in)            :: n
    real(dp), intent(in)           :: time
    integer, intent(in)            :: ivar
    character(len=200)             :: header_name, datfile_name
    integer                        :: my_unit, id, lo(2), hi(2)
    real(dp), allocatable          :: u_full(:, :)

    write(header_name, "(A,I6.6,A)") basename, n, ".bov"
    write(datfile_name, "(A,I6.6,A)") basename, n, ".dat"

    open(newunit=my_unit, file=trim(header_name), action="write")
    write(my_unit, *) "TIME: ", time
    write(my_unit, *) "DATA_FILE: ", trim(datfile_name)
    write(my_unit, *) "DATA_SIZE: ", bg%nx, 1
    if (dp == kind(0.0d0)) then
       write(my_unit, *) "DATA_FORMAT: DOUBLE"
    else if (dp == kind(0.0)) then
       write(my_unit, *) "DATA_FORMAT: FLOAT"
    else
       error stop "Unknown data format"
    end if
    write(my_unit, *) "VARIABLE: ", varname
    write(my_unit, *) "DATA_ENDIAN: LITTLE"
    write(my_unit, *) "CENTERING: zonal"
    write(my_unit, *) "BRICK_ORIGIN: 0. 0. 0."
    write(my_unit, *) "BRICK_SIZE: 1. 1. 1."
    write(my_unit, *) "DATA_COMPONENTS: 1"
    close(my_unit)

    allocate(u_full(bg%nx(1), bg%nx(2)))

    do id = 1, bg%n_blocks
       lo = (bg%id_to_index(:, id) - 1) * bg%bx + 1
       hi = lo + bg%bx - 1
       u_full(lo(1):hi(1), lo(2):hi(2)) = &
            bg%uu(IX(1:bg%bx(1), 1:bg%bx(2), id, ivar))
    end do

    open(newunit=my_unit, file=trim(datfile_name), form='unformatted', &
         access='stream', status='replace')
    write(my_unit) u_full
    close(my_unit)
    print *, "Written ", trim(header_name)

  end subroutine write_brick_of_values

end module m_blockgrid
