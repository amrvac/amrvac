  module physics
    implicit none

    procedure(sub_update), pointer  :: update   => null()

    abstract interface

       subroutine sub_update(a,b)
         double precision, intent(in)  :: a
         double precision, intent(out) :: b         
       end subroutine sub_update
       
    end interface
    
  end module physics

  module gut
    use physics
    implicit none
    
  contains

    subroutine gut_init()
      update => gut_update
    end subroutine gut_init
    
    subroutine gut_update(a,b)
      !$acc routine
      double precision, intent(in)  :: a
      double precision, intent(out) :: b
      b = 2.0d0 * a
    end subroutine gut_update
    
  end module gut

program procedurepointer
  use physics
  use gut
  implicit none

  integer, parameter               :: n=100
  double precision, dimension(n)   :: x,y
  integer                          :: i
  procedure(sub_update), pointer        :: local_update => Null()
  
  call gut_init()

  x = 666.0d0
  y = 0.0d0

  local_update => contained_update
  !$acc data copyin(x) copyout(y)
  !$acc parallel loop
  do i=1,n
!     y(i) = 2.0d0 * x(i)               ! works
!     call contained_update(x(i),y(i))  ! works
!     call gut_update(x(i),y(i))        ! works 
!     call update(x(i),y(i))            ! breaks: NVFORTRAN-W-0155-Data clause needed for exposed use of pointer update$sd (procedure_pointer.f90: 53)  -- followed by internal compler error if such measures are taken (e.g. attach)

    call local_update(x(i),y(i))        ! breaks: NVFORTRAN-W-0155-Data clause needed for exposed use of pointer local_update$sd (procedure_pointer.f90: 53) -- followed by internal compler error if such measures are taken (e.g. attach)

  end do
  !$acc end data

  print *, x,y

contains

    subroutine contained_update(a,b)
      !$acc routine
      double precision, intent(in)  :: a
      double precision, intent(out) :: b
      b = 2.0d0 * a
    end subroutine contained_update
  
  
end program procedurepointer
