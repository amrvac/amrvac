!> A Fortran 90 module for creating 1D and 2D lookup tables. These tables can be
!> used to efficiently interpolate one or more values.
!>
!> Author: Jannis Teunissen
module mod_lookup_table
  implicit none
  private

  ! The precision of the real numbers used in the tables
  integer, parameter :: dp = kind(1.0d0)

  ! ** Routines for finding indices in sorted lists **
  public :: find_index_linear
  public :: find_index_bsearch
  public :: find_index_adaptive

  !> The lookup table type. There can be one or more columns, for which values
  !> can be looked up for a given 'x-coordinate'.
  type LT_t
     integer  :: n_points !< The number of points
     integer  :: n_cols   !< The number of columns
     real(dp) :: x_min    !< The minimum lookup coordinate
     real(dp) :: dx       !< The x-spacing in the lookup coordinate
     real(dp) :: inv_dx   !< The inverse x-spacing

     ! The table is stored in two ways, to speed up different types of lookups.
     real(dp), allocatable :: cols_rows(:, :) !< The table in column-major order
     real(dp), allocatable :: rows_cols(:, :) !< The table in row-major order
  end type LT_t

  !> The 2D lookup table type
  type LT2_t
     integer               :: n_points(2) !< The size of the table
     integer               :: n_cols      !< The number of columns/variables
     real(dp)              :: x_min(2)    !< The minimum lookup coordinate
     real(dp)              :: dx(2)       !< The x-spacing in the lookup coordinate
     real(dp)              :: inv_dx(2)   !< The inverse x-spacing
     real(dp), allocatable :: rows_cols(:, :, :)
  end type LT2_t

  !> The 3D lookup table type
  type LT3_t
     integer               :: n_points(3) !< The size of the table
     integer               :: n_cols      !< The number of columns/variables
     real(dp)              :: x_min(3)    !< The minimum lookup coordinate
     real(dp)              :: dx(3)       !< The x-spacing in the lookup coordinate
     real(dp)              :: inv_dx(3)   !< The inverse x-spacing
     real(dp), allocatable :: rows_cols(:, :, :, :)
  end type LT3_t

  !> Type to indicate a location in the lookup table, which can be used to speed
  !> up multiple lookups of different columns.
  type LT_loc_t
     private
     integer  :: low_ix   !< The x-value lies between low_ix and low_ix+1
     real(dp) :: low_frac !< The distance from low_ix (up to low_ix+1), given
                          !< as a real number between 0 and 1.
  end type LT_loc_t

  !> Type to indicate a location in a 2D lookup table
  type LT2_loc_t
     private
     !> The x-value lies between low_ix and low_ix+1
     integer  :: low_ix(2)
     !> The distance from low_ix (up to low_ix+1), given as a real number
     !> between 0 and 1.
     real(dp) :: low_frac(2)
  end type LT2_loc_t

  !> Type to indicate a location in a 3D lookup table
  type LT3_loc_t
     private
     !> The x-value lies between low_ix and low_ix+1
     integer  :: low_ix(3)
     !> The distance from low_ix (up to low_ix+1), given as a real number
     !> between 0 and 1.
     real(dp) :: low_frac(3)
  end type LT3_loc_t

  ! Public types
  public :: LT_t
  public :: LT_loc_t
  public :: LT2_t
  public :: LT2_loc_t
  public :: LT3_t
  public :: LT3_loc_t

  ! Public methods
  public :: LT_create           ! Create a new lookup table
  public :: LT_create_from_data ! Create a new lookup table from existing data
  public :: LT_get_xdata        ! Get the x-values of a table
  public :: LT_get_spaced_data  ! Convert values to regularly spaced
  public :: LT_set_col          ! Set one table column
  public :: LT_add_col          ! Add a column
  public :: LT_get_loc          ! Get the index (row) of a value
  public :: LT_get_col          ! Interpolate one column
  public :: LT_get_mcol         ! Interpolate multiple columns
  public :: LT_get_col_at_loc   ! Get one column at location
  public :: LT_get_mcol_at_loc  ! Get multiple columns at location
  public :: LT_get_data         ! Get all the data of the table
  public :: LT_lin_interp_list  ! Linearly interpolate a list
  public :: LT_to_file          ! Store lookup table in file
  public :: LT_from_file        ! Restore lookup table from file

  ! Public methods
  public :: LT2_create           ! Create a new lookup table
  public :: LT2_create_from_data ! Create a new lookup table from existing data
  public :: LT2_set_col          ! Set one table column
  public :: LT2_get_loc          ! Get the index (row) of a value
  public :: LT2_get_col          ! Interpolate one column
  public :: LT2_get_col_at_loc   ! Get one column at location

  public :: LT3_create           ! Create a new lookup table
  public :: LT3_create_from_data ! Create a new lookup table from existing data
  public :: LT3_get_loc          ! Get the index (row) of a value
  public :: LT3_get_col          ! Interpolate one column
  public :: LT3_get_col_at_loc   ! Get one column at location

contains

  ! ** Routines for finding indices **

  !> Linear search of sorted list for the smallest ix such that list(ix) >= val.
  !> On failure, returns size(list)+1
  pure function find_index_linear(list, val) result(ix)
    real(dp), intent(in) :: list(:) !< Sorted list
    real(dp), intent(in) :: val     !< Value to search for
    integer              :: ix

    do ix = 1, size(list)
       if (list(ix) >= val) exit
    end do
  end function find_index_linear

  !> Binary search of sorted list for the smallest ix such that list(ix) >= val.
  !> On failure, returns size(list)+1
  pure function find_index_bsearch(list, val) result(ix)
    real(dp), intent(in) :: list(:) !< Sorted list
    real(dp), intent(in) :: val     !< Value to search for
    integer              :: ix, i_min, i_max, i_middle

    i_min = 1
    i_max = size(list)

    do while (i_min < i_max)
       ! This safely performs: i_middle = (i_max + i_min) / 2
       i_middle = i_min + ishft(i_max - i_min, -1)

       if (list(i_middle) >= val) then
          i_max = i_middle
       else
          i_min = i_middle + 1
       end if
    end do

    ix = i_min
    if (val > list(ix)) ix = size(list) + 1
  end function find_index_bsearch

  !> Adaptive search (combination of linear and binary search) of sorted list
  !> for the smallest ix such that list(ix) >= val. On failure, returns
  !> size(list)+1
  pure function find_index_adaptive(list, val) result(ix)
    real(dp), intent(in) :: list(:) !< Sorted list
    real(dp), intent(in) :: val     !< Value to search for
    integer              :: ix
    integer, parameter   :: binary_search_limit = 40

    if (size(list) < binary_search_limit) then
       ix = find_index_linear(list, val)
    else
       ix = find_index_bsearch(list, val)
    end if
  end function find_index_adaptive

  ! ** 1D lookup table routines **

  !> This function returns a new lookup table
  function LT_create(x_min, x_max, n_points, n_cols) result(my_lt)
    real(dp), intent(in) :: x_min  !< Minimum x-coordinate
    real(dp), intent(in) :: x_max  !< Maximum x-coordinate
    integer, intent(in)  :: n_points !< How many x-values to store
    integer, intent(in)  :: n_cols !< Number of variables that will be looked up
    type(LT_t)           :: my_lt

    if (x_max <= x_min) stop "set_xdata: x_max should be > x_min"
    if (n_points <= 1)    stop "set_xdata: n_points should be bigger than 1"

    my_lt%n_points = n_points
    my_lt%x_min    = x_min
    my_lt%dx       = (x_max - x_min) / (n_points - 1)
    my_lt%inv_dx   = 1 / my_lt%dx

    allocate(my_lt%cols_rows(n_cols, n_points))
    allocate(my_lt%rows_cols(n_points, n_cols))
    my_lt%cols_rows = 0
    my_lt%rows_cols = 0
    my_lt%n_cols    = n_cols
  end function LT_create

  !> This function returns a new lookup table from regularly spaced data
  function LT_create_from_data(x_min, x_max, spaced_data) result(my_lt)
    real(dp), intent(in) :: x_min             !< Minimum coordinate
    real(dp), intent(in) :: x_max             !< Maximum coordinate
    real(dp), intent(in) :: spaced_data(:, :) !< Input data of shape (n_points, n_cols)
    integer              :: n_points, n_cols
    type(LT_t)           :: my_lt

    n_points = size(spaced_data, 1)
    n_cols   = size(spaced_data, 2)

    if (x_max <= x_min) stop "LT_create error: x_max <= x_min"
    if (n_points <= 1) stop "LT_create error: n_points <= 1"

    my_lt%n_cols   = n_cols
    my_lt%n_points = n_points
    my_lt%x_min    = x_min
    my_lt%dx       = (x_max - x_min) / (n_points - 1)
    my_lt%inv_dx   = 1 / my_lt%dx

    my_lt%rows_cols = spaced_data
    my_lt%cols_rows = transpose(spaced_data)
  end function LT_create_from_data

  !> Returns the x-coordinates of the lookup table
  pure function LT_get_xdata(x_min, dx, n_points) result(xdata)
    real(dp), intent(in) :: x_min, dx
    integer, intent(in)  :: n_points
    real(dp)             :: xdata(n_points)
    integer              :: ix

    do ix = 1, n_points
       xdata(ix) = x_min + (ix-1) * dx
    end do
  end function LT_get_xdata

  !> Linearly interpolate the (x, y) input data to the new_x coordinates
  pure function LT_get_spaced_data(in_x, in_y, new_x) result(out_yy)
    real(dp), intent(in) :: in_x(:), in_y(:), new_x(:)
    real(dp)             :: out_yy(size(new_x))
    integer              :: ix

    do ix = 1, size(new_x)
       call LT_lin_interp_list(in_x, in_y, new_x(ix), out_yy(ix))
    end do
  end function LT_get_spaced_data

  !> Fill the column with index col_ix using the linearly interpolated (x, y)
  !> data
  pure subroutine LT_set_col(my_lt, col_ix, x, y)
    type(LT_t), intent(inout) :: my_lt
    integer, intent(in)       :: col_ix
    real(dp), intent(in)      :: x(:), y(:)
    my_lt%cols_rows(col_ix, :) = LT_get_spaced_data(x, y, &
         LT_get_xdata(my_lt%x_min, my_lt%dx, my_lt%n_points))
    my_lt%rows_cols(:, col_ix) = my_lt%cols_rows(col_ix, :)
  end subroutine LT_set_col

  !> Add a new column by linearly interpolating the (x, y) data
  pure subroutine LT_add_col(my_lt, x, y)
    type(LT_t), intent(inout) :: my_lt
    real(dp), intent(in)      :: x(:), y(:)
    type(LT_t)                :: temp_lt

    temp_lt = my_lt
    deallocate(my_lt%cols_rows)
    deallocate(my_lt%rows_cols)
    allocate(my_lt%cols_rows(my_lt%n_cols+1, my_lt%n_points))
    allocate(my_lt%rows_cols(my_lt%n_points, my_lt%n_cols+1))

    my_lt%cols_rows(1:my_lt%n_cols, :) = temp_lt%cols_rows
    my_lt%rows_cols(:, 1:my_lt%n_cols) = temp_lt%rows_cols
    my_lt%n_cols                       = my_lt%n_cols + 1
    my_lt%cols_rows(my_lt%n_cols, :)   = LT_get_spaced_data(x, y, &
         LT_get_xdata(my_lt%x_min, my_lt%dx, my_lt%n_points))
    my_lt%rows_cols(:, my_lt%n_cols)   = my_lt%cols_rows(my_lt%n_cols, :)
  end subroutine LT_add_col

  !> Get a location in the lookup table
  elemental function LT_get_loc(my_lt, x) result(my_loc)
    type(LT_t), intent(in) :: my_lt
    real(dp), intent(in)   :: x
    type(LT_loc_t)         :: my_loc
    real(dp)               :: frac

    frac            = (x - my_lt%x_min) * my_lt%inv_dx
    my_loc%low_ix   = ceiling(frac)
    my_loc%low_frac = my_loc%low_ix - frac

    ! Check bounds
    if (my_loc%low_ix < 1) then
       my_loc%low_ix   = 1
       my_loc%low_frac = 1
    else if (my_loc%low_ix >= my_lt%n_points) then
       my_loc%low_ix   = my_lt%n_points - 1
       my_loc%low_frac = 0
    end if
  end function LT_get_loc

  !> Get the values of all columns at x
  pure function LT_get_mcol(my_lt, x) result(col_values)
    type(LT_t), intent(in) :: my_lt
    real(dp), intent(in)   :: x
    real(dp)               :: col_values(my_lt%n_cols)
    type(LT_loc_t)         :: loc

    loc        = LT_get_loc(my_lt, x)
    col_values = LT_get_mcol_at_loc(my_lt, loc)
  end function LT_get_mcol

  !> Get the value of a single column at x
  elemental function LT_get_col(my_lt, col_ix, x) result(col_value)
    type(LT_t), intent(in) :: my_lt
    integer, intent(in)    :: col_ix
    real(dp), intent(in)   :: x
    real(dp)               :: col_value
    type(LT_loc_t)         :: loc

    loc       = LT_get_loc(my_lt, x)
    col_value = LT_get_col_at_loc(my_lt, col_ix, loc)
  end function LT_get_col

  !> Get the values of all columns at a location
  pure function LT_get_mcol_at_loc(my_lt, loc) result(col_values)
    type(LT_t), intent(in)     :: my_lt
    type(LT_loc_t), intent(in) :: loc
    real(dp)                   :: col_values(my_lt%n_cols)

    col_values = loc%low_frac * my_lt%cols_rows(:, loc%low_ix) + &
         (1-loc%low_frac) * my_lt%cols_rows(:, loc%low_ix+1)
  end function LT_get_mcol_at_loc

  !> Get the value of a single column at a location
  elemental function LT_get_col_at_loc(my_lt, col_ix, loc) result(col_value)
    type(LT_t), intent(in)     :: my_lt
    integer, intent(in)        :: col_ix
    type(LT_loc_t), intent(in) :: loc
    real(dp)                   :: col_value

    col_value = loc%low_frac * my_lt%rows_cols(loc%low_ix, col_ix) + &
         (1-loc%low_frac) * my_lt%rows_cols(loc%low_ix+1, col_ix)
  end function LT_get_col_at_loc

  !> Get the x-coordinates and the columns of the lookup table
  pure subroutine LT_get_data(my_lt, x_data, cols_rows)
    type(LT_t), intent(in) :: my_lt
    real(dp), intent(out)  :: x_data(:), cols_rows(:, :)

    x_data    = LT_get_xdata(my_lt%x_min, my_lt%dx, my_lt%n_points)
    cols_rows = my_lt%cols_rows
  end subroutine LT_get_data

  !> Compute by use of linear interpolation the value in the middle
  ! of a domain D = [x_list(1) , x_list(size(x_list))].
  ! If x_value is left of domain  D,
  ! then the value becomes the value at the left side of D,
  ! if x_value is right of domain D,
  ! then the value becomes the value at the rigth side of D
  pure subroutine LT_lin_interp_list(x_list, y_list, x_value, y_value)
    real(dp), intent(in)  :: x_list(:), y_list(:)
    real(dp), intent(in)  :: x_value
    real(dp), intent(out) :: y_value

    integer               :: ix, iMin, iMax
    real(dp)              :: temp

    iMin = 1
    iMax = size(x_list)

    if (x_value <= x_list(iMin)) then
       y_value = y_list(iMin)
    else if (x_value >= x_list(iMax)) then
       y_value = y_list(iMax)
    else
       ix = find_index_adaptive(x_list, x_value)
       temp = (x_value - x_list(ix-1)) / (x_list(ix) - x_list(ix-1))
       y_value = (1 - temp) * y_list(ix-1) + temp * y_list(ix)
    end if
  end subroutine LT_lin_interp_list

  !> Write the lookup table to file (in binary, potentially unportable)
  subroutine LT_to_file(my_lt, filename)
    type(LT_t), intent(in)       :: my_lt
    character(len=*), intent(in) :: filename
    integer                      :: my_unit

    open(newunit=my_unit, file=trim(filename), form='UNFORMATTED', &
         access='STREAM', status='REPLACE')
    write(my_unit) my_lt%n_points, my_lt%n_cols
    write(my_unit) my_lt%x_min, my_lt%dx, my_lt%inv_dx
    write(my_unit) my_lt%cols_rows
    close(my_unit)
  end subroutine LT_to_file

  !> Read the lookup table from file (in binary, potentially unportable)
  subroutine LT_from_file(my_lt, filename)
    type(LT_t), intent(inout)    :: my_lt
    character(len=*), intent(in) :: filename
    integer                      :: my_unit

    open(newunit=my_unit, file=trim(filename), form='UNFORMATTED', &
         access='STREAM', status='OLD')
    read(my_unit) my_lt%n_points, my_lt%n_cols
    read(my_unit) my_lt%x_min, my_lt%dx, my_lt%inv_dx

    allocate(my_lt%cols_rows(my_lt%n_cols, my_lt%n_points))
    allocate(my_lt%rows_cols(my_lt%n_points, my_lt%n_cols))

    read(my_unit) my_lt%cols_rows
    my_lt%rows_cols = transpose(my_lt%cols_rows)

    close(my_unit)
  end subroutine LT_from_file

  ! ** 2D lookup table routines **

  !> This function returns a new lookup table
  function LT2_create(x_min, x_max, n_points, n_cols) result(my_lt)
    real(dp), intent(in) :: x_min(2)    !< Minimum coordinate
    real(dp), intent(in) :: x_max(2)    !< Maximum coordinate
    integer, intent(in)  :: n_points(2) !< How many values to store
    integer, intent(in)  :: n_cols      !< Number of variables that will be looked up
    type(LT2_t)          :: my_lt

    if (any(x_max <= x_min)) stop "LT2_create error: x_max <= x_min"
    if (any(n_points <= 1)) stop "LT2_create error: n_points <= 1"

    my_lt%n_points = n_points
    my_lt%x_min    = x_min
    my_lt%dx       = (x_max - x_min) / (n_points - 1)
    my_lt%inv_dx   = 1 / my_lt%dx

    allocate(my_lt%rows_cols(n_points(1), n_points(2), n_cols))
    my_lt%rows_cols = 0
    my_lt%n_cols    = n_cols
  end function LT2_create

  !> This function returns a new lookup table from regularly spaced data
  function LT2_create_from_data(x_min, x_max, spaced_data) result(my_lt)
    real(dp), intent(in) :: x_min(2)             !< Minimum coordinate
    real(dp), intent(in) :: x_max(2)             !< Maximum coordinate
    real(dp), intent(in) :: spaced_data(:, :, :) !< Input data of shape (nx, ny, n_cols)
    integer              :: n_points(2), n_cols
    type(LT2_t)          :: my_lt

    n_points(1) = size(spaced_data, 1)
    n_points(2) = size(spaced_data, 2)
    n_cols      = size(spaced_data, 3)

    if (any(x_max <= x_min)) stop "LT2_create error: x_max <= x_min"
    if (any(n_points <= 1)) stop "LT2_create error: n_points <= 1"

    my_lt%n_cols   = n_cols
    my_lt%n_points = n_points
    my_lt%x_min    = x_min
    my_lt%dx       = (x_max - x_min) / (n_points - 1)
    my_lt%inv_dx   = 1 / my_lt%dx

    my_lt%rows_cols = spaced_data
  end function LT2_create_from_data

  !> Fill the column with index col_ix using linearly interpolated data
  pure subroutine LT2_set_col(my_lt, col_ix, x1, x2, y)
    type(LT2_t), intent(inout) :: my_lt
    integer, intent(in)        :: col_ix
    real(dp), intent(in)       :: x1(:), x2(:), y(:, :)
    real(dp), allocatable      :: tmp(:, :), c1(:), c2(:)
    integer                    :: ix

    allocate(c1(my_lt%n_points(1)),c2(my_lt%n_points(2)))
    c1 = LT_get_xdata(my_lt%x_min(1), my_lt%dx(1), my_lt%n_points(1))
    c2 = LT_get_xdata(my_lt%x_min(2), my_lt%dx(2), my_lt%n_points(2))
    allocate(tmp(my_lt%n_points(1), size(x2)))

    ! Interpolate along first coordinate
    do ix = 1, size(x2)
       tmp(:, ix) = LT_get_spaced_data(x1, y(:, ix), c1)
    end do

    ! Interpolate along second coordinate
    do ix = 1, my_lt%n_points(1)
       my_lt%rows_cols(ix, :, col_ix) = &
            LT_get_spaced_data(x2, tmp(ix, :), c2)
    end do
    deallocate(c1,c2)
  end subroutine LT2_set_col

  !> Get a location in the lookup table
  elemental function LT2_get_loc(my_lt, x1, x2) result(my_loc)
    type(LT2_t), intent(in) :: my_lt
    real(dp), intent(in)    :: x1, x2
    type(LT2_loc_t)         :: my_loc
    real(dp)                :: frac(2)

    frac            = ([x1, x2] - my_lt%x_min) * my_lt%inv_dx
    my_loc%low_ix   = ceiling(frac)
    my_loc%low_frac = my_loc%low_ix - frac

    ! Check bounds
    where (my_loc%low_ix < 1)
       my_loc%low_ix   = 1
       my_loc%low_frac = 1
    end where

    where (my_loc%low_ix >= my_lt%n_points)
       my_loc%low_ix   = my_lt%n_points - 1
       my_loc%low_frac = 0
    end where
  end function LT2_get_loc

  !> Get the value of a single column at x
  elemental function LT2_get_col(my_lt, col_ix, x1, x2) result(col_value)
    type(LT2_t), intent(in) :: my_lt
    integer, intent(in)     :: col_ix
    real(dp), intent(in)    :: x1, x2
    real(dp)                :: col_value
    type(LT2_loc_t)         :: loc

    loc       = LT2_get_loc(my_lt, x1, x2)
    col_value = LT2_get_col_at_loc(my_lt, col_ix, loc)
  end function LT2_get_col

  !> Get the value of a single column at a location
  elemental function LT2_get_col_at_loc(my_lt, col_ix, loc) result(col_value)
    type(LT2_t), intent(in)     :: my_lt
    integer, intent(in)         :: col_ix
    type(LT2_loc_t), intent(in) :: loc
    integer                     :: ix(2)
    real(dp)                    :: w(2, 2)
    real(dp)                    :: col_value

    ! Bilinear interpolation
    w(1, 1) = loc%low_frac(1) * loc%low_frac(2)
    w(2, 1) = (1 - loc%low_frac(1)) * loc%low_frac(2)
    w(1, 2) = loc%low_frac(1) * (1 - loc%low_frac(2))
    w(2, 2) = (1 - loc%low_frac(1)) * (1 - loc%low_frac(2))
    ix = loc%low_ix

    col_value = sum(w * my_lt%rows_cols(ix(1):ix(1)+1, &
         ix(2):ix(2)+1, col_ix))
  end function LT2_get_col_at_loc

  ! ** 3D lookup table routines **

  !> This function returns a new lookup table
  function LT3_create(x_min, x_max, n_points, n_cols) result(my_lt)
    real(dp), intent(in) :: x_min(3)    !< Minimum coordinate
    real(dp), intent(in) :: x_max(3)    !< Maximum coordinate
    integer, intent(in)  :: n_points(3) !< How many values to store
    integer, intent(in)  :: n_cols      !< Number of variables that will be looked up
    type(LT3_t)          :: my_lt

    if (any(x_max <= x_min)) stop "LT3_create error: x_max <= x_min"
    if (any(n_points <= 1)) stop "LT3_create error: n_points <= 1"

    my_lt%n_points = n_points
    my_lt%x_min    = x_min
    my_lt%dx       = (x_max - x_min) / (n_points - 1)
    my_lt%inv_dx   = 1 / my_lt%dx

    allocate(my_lt%rows_cols(n_points(1), n_points(2), n_points(3), n_cols))
    my_lt%rows_cols = 0
    my_lt%n_cols    = n_cols
  end function LT3_create

  !> This function returns a new lookup table from regularly spaced data
  function LT3_create_from_data(x_min, x_max, spaced_data) result(my_lt)
    real(dp), intent(in) :: x_min(3)                !< Minimum coordinate
    real(dp), intent(in) :: x_max(3)                !< Maximum coordinate
    real(dp), intent(in) :: spaced_data(:, :, :, :) !< Input data of shape (nx, ny, nz, n_cols)
    integer              :: n_points(3), n_cols
    type(LT3_t)          :: my_lt

    n_points(1) = size(spaced_data, 1)
    n_points(2) = size(spaced_data, 2)
    n_points(3) = size(spaced_data, 3)
    n_cols      = size(spaced_data, 4)

    if (any(x_max <= x_min)) stop "LT3_create error: x_max <= x_min"
    if (any(n_points <= 1)) stop "LT3_create error: n_points <= 1"

    my_lt%n_cols   = n_cols
    my_lt%n_points = n_points
    my_lt%x_min    = x_min
    my_lt%dx       = (x_max - x_min) / (n_points - 1)
    my_lt%inv_dx   = 1 / my_lt%dx

    my_lt%rows_cols = spaced_data
  end function LT3_create_from_data

  !> Get a location in the lookup table
  elemental function LT3_get_loc(my_lt, x1, x2, x3) result(my_loc)
    type(LT3_t), intent(in) :: my_lt
    real(dp), intent(in)    :: x1, x2, x3
    type(LT3_loc_t)         :: my_loc
    real(dp)                :: frac(3)

    frac            = ([x1, x2, x3] - my_lt%x_min) * my_lt%inv_dx
    my_loc%low_ix   = ceiling(frac)
    my_loc%low_frac = my_loc%low_ix - frac

    ! Check bounds
    where (my_loc%low_ix < 1)
       my_loc%low_ix   = 1
       my_loc%low_frac = 1
    end where

    where (my_loc%low_ix >= my_lt%n_points)
       my_loc%low_ix   = my_lt%n_points - 1
       my_loc%low_frac = 0
    end where
  end function LT3_get_loc

  !> Get the value of a single column at x
  elemental function LT3_get_col(my_lt, col_ix, x1, x2, x3) result(col_value)
    type(LT3_t), intent(in) :: my_lt
    integer, intent(in)     :: col_ix
    real(dp), intent(in)    :: x1, x2, x3
    real(dp)                :: col_value
    type(LT3_loc_t)         :: loc

    loc       = LT3_get_loc(my_lt, x1, x2, x3)
    col_value = LT3_get_col_at_loc(my_lt, col_ix, loc)
  end function LT3_get_col

  !> Get the value of a single column at a location
  elemental function LT3_get_col_at_loc(my_lt, col_ix, loc) result(col_value)
    type(LT3_t), intent(in)     :: my_lt
    integer, intent(in)         :: col_ix
    type(LT3_loc_t), intent(in) :: loc
    integer                     :: ix(3)
    real(dp)                    :: w(2, 2, 2)
    real(dp)                    :: col_value

    ! Bilinear interpolation
    w(1, 1, 1) = loc%low_frac(1) * loc%low_frac(2)
    w(2, 1, 1) = (1 - loc%low_frac(1)) * loc%low_frac(2)
    w(1, 2, 1) = loc%low_frac(1) * (1 - loc%low_frac(2))
    w(2, 2, 1) = (1 - loc%low_frac(1)) * (1 - loc%low_frac(2))
    w(:, :, 2) = w(:, :, 1) * (1 - loc%low_frac(3))
    w(:, :, 1) = w(:, :, 1) * loc%low_frac(3)
    ix = loc%low_ix

    col_value = sum(w * my_lt%rows_cols(ix(1):ix(1)+1, &
         ix(2):ix(2)+1, ix(3):ix(3)+1, col_ix))
  end function LT3_get_col_at_loc

end module mod_lookup_table
