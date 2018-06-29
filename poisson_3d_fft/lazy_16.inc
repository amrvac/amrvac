!!****h* BigDFT/lazy_16
!! NAME
!!   lazy_16
!!
!! FUNCTION
!!   Filters for interpolating scaling functions (order 16)
!!
!! SOURCE
!!
integer, parameter :: m=18
real(kind=8), dimension(-m:m) :: ch = (/&
     0.d0,0.d0,0.d0,-6.39259815216064453D-6,0.D0,0.000110641121864318848D0,0.D0,&
     -0.000915303826332092285D0,0.D0,0.00484772026538848877D0,0.D0,&
     -0.0186983495950698853D0,0.D0,0.0575909167528152466D0,0.D0,&
     -0.159974768757820129D0,0.D0,0.617045536637306213D0,1.D0,0.617045536637306213D0,&
     0D0,-0.159974768757820129D0,0.D0,0.0575909167528152466D0,0.D0,&
     -0.0186983495950698853D0,0.D0,0.00484772026538848877D0,0.D0,&
     -0.000915303826332092285D0,0.D0,0.000110641121864318848D0,0.D0,&
     -6.39259815216064453D-6,0.d0,0.d0,0.d0&
     /)
real(kind=8), dimension(-m:m) :: cg,cht,cgt

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

