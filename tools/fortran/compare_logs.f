! Purpose: Compare two text files containing a 2D table of number. Return 0 if
! they are the same within a numerical threshold, or 1 if they differ.
!
! Author: Jannis Teunissen
! Usage: ./compare_logs rtol atol a.log b.log
program compare_logs

  integer, parameter    :: dp = kind(0.0d0)
  real(dp), allocatable :: data1(:, :), data2(:, :)
  character(len=200)    :: arg_string
  real(dp)              :: rtol ! Relative tolerance
  real(dp)              :: atol ! Absolute tolerance

  ! Maximum line length
  integer, parameter :: line_len = 400

  ! Maximum number of columns
  integer, parameter :: max_cols = 50

  ! Maximum number of rows
  integer, parameter :: max_rows = 10000

  if (command_argument_count() /= 4) then
     print *, "Usage: ./compare_logs rtol atol a.log b.log"
     print *, "rtol: relative tolerance"
     print *, "atol: absolute tolerance"
     error stop "Incorrect arguments"
  end if

  ! Read relative tolerance
  call get_command_argument(1, arg_string)
  read(arg_string, *) rtol

  ! Read absolute tolerance
  call get_command_argument(2, arg_string)
  read(arg_string, *) atol

  ! Load files
  call get_command_argument(3, arg_string)
  call read_log(trim(arg_string), data1)

  call get_command_argument(4, arg_string)
  call read_log(trim(arg_string), data2)

  if (any(shape(data1) /= shape(data2))) then
     error stop "Data has different shape"
  else if (any(abs(data1 - data2) > atol + &
       0.5_dp * rtol * (abs(data1) + abs(data2)))) then
     !error stop "Difference detected: |d1-d2| > atol + 0.5*rtol*(|d1|+|d2|)"
     call exit(1)               ! Non-standard but does not lead to backtrace
  end if

contains

  subroutine read_log(fname, dd)
    character(len=*), intent(in)         :: fname
    real(dp), allocatable, intent(inout) :: dd(:, :)
    integer                              :: n, n_cols
    integer, parameter                   :: myunit     = 100
    character(len=*), parameter          :: separators = " ,'"""//char(9)
    character(len=line_len)              :: line
    integer                              :: ix_start(max_cols)
    integer                              :: ix_end(max_cols)
    real(dp), allocatable                :: dbuf(:, :)

    open(myunit, file=fname, status="old", action="read")

    ! Skip header lines with # or !
    do
       read(myunit, "(A)") line
       line = adjustl(line)
       if (line(1:1) == '#' .or. line(1:1) == '!') then
          cycle
       else
          exit
       end if
    end do

    call get_fields_string(line, separators, &
         max_cols, n_cols, ix_start, ix_end)

    allocate(dbuf(n_cols, max_rows))

    do n = 1, max_rows
       read(line, *) dbuf(:, n)
       read(myunit, "(A)", end=888) line
    end do

    error stop "Log file has too many rows"

888 continue

    allocate(dd(n_cols, n))
    dd = dbuf(:, 1:n)

  end subroutine read_log

  !> Routine to find the indices of entries in a string
  subroutine get_fields_string(line, delims, n_max, n_found, ixs_start, ixs_end)
    !> The line from which we want to read
    character(len=*), intent(in)  :: line
    !> A string with delimiters. For example delims = " ,'"""//char(9)
    character(len=*), intent(in)  :: delims
    !> Maximum number of entries to read in
    integer, intent(in)           :: n_max
    !> Number of entries found
    integer, intent(inout)        :: n_found
    !> On return, ix_start(i) holds the starting point of entry i
    integer, intent(inout)        :: ixs_start(n_max)
    !> On return, ix_end(i) holds the end point of entry i
    integer, intent(inout)        :: ixs_end(n_max)

    integer                       :: ix, ix_prev

    ix_prev = 0
    n_found = 0

    do while (n_found < n_max)

       ! Find the starting point of the next entry (a non-delimiter value)
       ix = verify(line(ix_prev+1:), delims)
       if (ix == 0) exit

       n_found            = n_found + 1
       ixs_start(n_found) = ix_prev + ix ! This is the absolute position in 'line'

       ! Get the end point of the current entry (next delimiter index minus one)
       ix = scan(line(ixs_start(n_found)+1:), delims) - 1

       if (ix == -1) then              ! If there is no last delimiter,
          ixs_end(n_found) = len(line) ! the end of the line is the endpoint
       else
          ixs_end(n_found) = ixs_start(n_found) + ix
       end if

       ix_prev = ixs_end(n_found) ! We continue to search from here
    end do

  end subroutine get_fields_string

end program compare_logs
