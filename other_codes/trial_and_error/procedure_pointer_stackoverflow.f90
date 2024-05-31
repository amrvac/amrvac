  !Example was reported as working with gfortran-11 at
  !https://stackoverflow.com/questions/70818524/fortran-openacc-invoking-a-function-on-device-using-a-function-pointer


  module modTest2

  use openacc
  implicit none

  
  type :: Container
    sequence
    integer :: n
    integer, allocatable :: arr(:)
  end type Container


  interface Container
    procedure :: new_Container
  end interface 


  abstract interface
    integer function function_template (i)
      integer, intent (in) :: i
    end function function_template
  end interface


  contains
  
    type(Container) function new_Container(n)
      integer, intent(in) :: n

      allocate(new_Container%arr(n))
    end function new_Container    
end module modTest2

program test2

  use modTest2
  implicit none


  integer :: n, x, i
  type(Container) :: c
  procedure(function_template), pointer :: init


  print *, "Enter array size: "
  read *, n
  print *, "Allocating..."
  c = Container(n)
  print *, "Allocation complete!"
  
  
  print *, "Enter initialization type (x): "
  read *, x
  print *, "Initializing..."
  select case (x)
    case (0)
      init => init0
    case default
      init => init1
  end select
  !$acc enter data copyin(c)
  !$acc enter data create(c%arr)
  !$acc parallel loop present(c)
  do i = 1, n
    c%arr(i) = init(i)
  end do
  !$acc exit data copyout(c%arr)
  !$acc exit data delete(c)
  print *, "Initialization complete..."


  do i = 1, n
    print *, i, c%arr(i)
  end do

  
  contains
  
    integer function init0(i)
      !$acc routine
      integer, intent(in) :: i
      init0 = 10*i
    end function init0

    
    integer function init1(i)
      !$acc routine
      integer, intent(in) :: i
      init1 = 20*i
    end function init1
end program test2
