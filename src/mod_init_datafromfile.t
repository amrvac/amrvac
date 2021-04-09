!> Module to set (or derive) initial conditions from user data
!>    We read in a vtk file that provides values on grid
module mod_init_datafromfile
  use mod_lookup_table

  implicit none
  private

  integer, parameter :: max_valuesfromfile = 10

  type read_data_values
     integer                        :: n_variables
     double precision               :: origin(3)
     double precision               :: dx(3)
     integer                        :: n_points(3)
     character(len=40), allocatable :: names(:)
     double precision, allocatable  :: values(:, :, :, :)
  end type read_data_values

  type(LT_t)  :: lt_1d(max_valuesfromfile)
  type(LT2_t) :: lt_2d(max_valuesfromfile)
  type(LT3_t) :: lt_3d(max_valuesfromfile)

  public :: read_data_init
  public :: read_data_set
  public :: read_data_get_1d
  public :: read_data_get_2d
  public :: read_data_get_3d

contains

  subroutine read_data_init()
    use mod_global_parameters

    integer                :: ivar
    double precision       :: xmax(3)
    type(read_data_values) :: bc

    if(usr_filename/='')then
       if(mype==0)then
          print *,'Reading data from input usr_filename=',usr_filename
       endif
       call read_vtk_structured_points(trim(usr_filename), bc)
       xmax=bc%origin + (bc%n_points-1) * bc%dx
       if(mype==0)then
          print *,'Obtained data in ',bc%n_points,' grid with spacings ',bc%dx,' for number of vars=',bc%n_variables
       endif

       do ivar = 1,bc%n_variables
       {^IFONED
          lt_1d(ivar) = LT_create_from_data(bc%origin(1), &
                     xmax(1), bc%values(:, 1, 1:1, ivar))
       }
       {^IFTWOD
          lt_2d(ivar) = LT2_create_from_data(bc%origin(1:ndim), &
                     xmax(1:ndim), bc%values(:, :, 1:1, ivar))
       }
       {^IFTHREED
          print *,'todo'
          !!lt_3d(ivar) = LT3_create_from_data(bc%origin(1:ndim), &
          !!           xmax(1:ndim), bc%values(:, :, :, ivar))
       }
       end do
    endif

  end subroutine read_data_init

  elemental function read_data_get_3d(ivar, x1, x2, x3) result(val)
    integer, intent(in)          :: ivar
    double precision, intent(in) :: x1, x2, x3
    double precision             :: val

    val = LT3_get_col(lt_3d(ivar), 1, x1, x2, x3)
  end function read_data_get_3d

  elemental function read_data_get_2d(ivar, x1, x2) result(val)
    integer, intent(in)          :: ivar
    double precision, intent(in) :: x1, x2
    double precision             :: val

    val = LT2_get_col(lt_2d(ivar), 1, x1, x2)
  end function read_data_get_2d

  elemental function read_data_get_1d(ivar, x1) result(val)
    integer, intent(in)          :: ivar
    double precision, intent(in) :: x1
    double precision             :: val

    val = LT_get_col(lt_1d(ivar), 1, x1)
  end function read_data_get_1d

  !> Set values according to user data
  subroutine read_data_set(ixI^L,ixO^L,x,val,nval)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L, nval
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out)   :: val(ixI^S,1:nval)
    integer                         :: ival

    {^IFONED
    do ival=1,nval
       val(ixOmin1:ixOmax1,ival) = read_data_get_1d(ival,x(ixOmin1:ixOmax1,1))
    enddo
    }
    {^IFTWOD
    do ival=1,nval
       val(ixOmin1:ixOmax1, ixOmin2:ixOmax2, ival) = read_data_get_2d(ival, &
                  x(ixOmin1:ixOmax1, ixOmin2:ixOmax2, 1), x(ixOmin1:ixOmax1, ixOmin2:ixOmax2, 2))
    enddo
    }
    {^IFTHREED
    do ival=1,nval
       val(ixOmin1:ixOmax1, ixOmin2:ixOmax2, ixOmin3:ixOmax3, ival) = read_data_get_3d(ival, &
                  x(ixOmin1:ixOmax1, ixOmin2:ixOmax2, ixOmin3:ixOmax3, 1), &
                  x(ixOmin1:ixOmax1, ixOmin2:ixOmax2, ixOmin3:ixOmax3, 2), &
                  x(ixOmin1:ixOmax1, ixOmin2:ixOmax2, ixOmin3:ixOmax3, 3))
    enddo
    }

  end subroutine read_data_set

  subroutine read_vtk_structured_points(fname, bc)
    character(len=*),        intent(in)  :: fname
    type(read_data_values), intent(out)  :: bc
    double precision, allocatable :: tmp_data(:, :, :, :)
    integer, parameter            :: max_variables = 10
    character(len=40)             :: tmp_names(max_variables)
    character(len=256)            :: line
    character(len=40)             :: word, typename
    integer, parameter            :: my_unit = 123
    integer                       :: n, n_points_total
    integer                       :: n_components

    open(my_unit, file=trim(fname), status="old", action="read")

    ! Header, e.g. # vtk DataFile Version 2.0
    read(my_unit, "(A)") line

    ! Dataset name
    read(my_unit, "(A)") line

    ! ASCII / BINARY
    read(my_unit, "(A)") line

    if (line /= "ASCII") then
       print *, "line: ", trim(line)
       error stop "read_vtk: not an ASCII file"
    end if

    ! DATASET STRUCTURED_POINTS
    read(my_unit, "(A)") line

    if (line /= "DATASET STRUCTURED_POINTS") then
       print *, "line: ", trim(line)
       error stop "read_vtk must have: DATASET STRUCTURED_POINTS"
    end if

    ! DIMENSIONS NX NY NZ
    read(my_unit, "(A)") line
    read(line, *) word, bc%n_points

    if (word /= "DIMENSIONS") then
       print *, "line: ", trim(line)
       error stop "read_vtk expects: DIMENSIONS"
    end if

    ! SPACING DX DY DZ
    read(my_unit, *) word, bc%dx

    if (word /= "SPACING") then
       print *, "line: ", trim(line)
       error stop "read_vtk expects: SPACING"
    end if

    ! ORIGIN 0 0 0
    read(my_unit, *) word, bc%origin
    if (word /= "ORIGIN") then
       print *, "line: ", trim(line)
       error stop "read_vtk expects: ORIGIN"
    end if

    ! POINT_DATA NPOINTS
    read(my_unit, *) word, n_points_total

    if (word /= "POINT_DATA") then
       print *, "line: ", trim(line)
       error stop "read_vtk expects: POINT_DATA n_points"
    end if

    if (n_points_total /= product(bc%n_points)) &
         error stop "read_vtk: n_points not equal to NX*NY*NZ"

    allocate(tmp_data(bc%n_points(1), bc%n_points(2), bc%n_points(3), &
         max_variables))

    ! Read all scalar variables
    do n = 1, max_variables

       ! SCALARS name type ncomponents
       read(my_unit, *, end=900) word, tmp_names(n), typename, n_components

       if (word /= "SCALARS") then
          print *, "line: ", trim(line)
          error stop "read_vtk expects: SCALARS name type ncomponents"
       end if

       if (n_components /= 1) error stop "read_vtk: ncomponents should be 1"

       ! LOOKUP_TABLE default
       read(my_unit, *) word, typename

       if (word /= "LOOKUP_TABLE") then
          print *, "line: ", trim(line)
          error stop "read_vtk expects: LOOKUP_TABLE name"
       end if

       ! Read list of data values
       read(my_unit, *) tmp_data(:, :, :, n)
    end do

    ! Done reading variables
900 continue

    close(my_unit)

    if (n == max_variables + 1) &
         error stop "read_vtk: increase max_variables"

    ! Loop index is one higher than number of variables
    bc%n_variables = n-1
    bc%values      = tmp_data(:, :, :, 1:n-1)
    bc%names       = tmp_names(1:n-1)
  end subroutine read_vtk_structured_points

end module mod_init_datafromfile
