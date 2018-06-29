!!****h* BigDFT/lazy_14
!! NAME
!!   lazy_14
!!
!! FUNCTION
!!   Filters for interpolating scaling functions (order 14)
!!
!! SOURCE
!!
integer, parameter :: m=16
real(kind=8), dimension(-m:m) :: ch = (/ &
     0.d0,0.d0,0.d0,0.0000275373458862304687D0,0.D0,-0.000423073768615722656D0,0.D0,&
     0.00310254096984863281D0,0.D0,-0.0146262645721435547D0,0.D0,&
     0.0511919260025024414D0,0.D0,-0.153575778007507324D0,0.D0,0.614303112030029297D0,&
     1.D0,0.614303112030029297D0,0.D0,-0.153575778007507324D0,0.D0,&
     0.0511919260025024414D0,0.D0,-0.0146262645721435547D0,0.D0,&
     0.00310254096984863281D0,0.D0,-0.000423073768615722656D0,0.D0,&
     0.0000275373458862304687D0,0.d0,0.d0,0.d0&
     /)
real(kind=8), dimension(-m:m) ::  cg,cht,cgt

!******** coefficients for wavelet transform *********************
do i=-m,m
   cht(i)=0.d0
   cg(i)=0.d0
   cgt(i)=0.d0
enddo

! the normalization is chosen such that a constant function remains the same constant 
! on each level of the transform

cht( 0)=1.D0

! g coefficients from h coefficients
do i=-m,m-1
   cg(i+1)=cht(-i)*(-1.d0)**(i+1)
   cgt(i+1)=ch(-i)*(-1.d0)**(i+1)
enddo
!!***

