!!****h* BigDFT/lazy_20
!! NAME
!!   lazy_20
!!
!! FUNCTION
!!   Filters for interpolating scaling functions (order 30)
!!
!! SOURCE
!!
integer, parameter :: m=22
real(kind=8), dimension(-m:m) :: ch = (/ &
     0.d0, 0.d0, 0.d0, &
     -3.53757059201598167d-7,0.d0,7.51213519833981991d-6,0.d0, &
     -0.0000766237790230661631d0,0.d0,0.000501001632073894143d0,0.d0, &
     -0.00236837135162204504d0,0.d0,0.00868402828928083181d0,0.d0,&
     -0.0260520848678424954d0,0.d0,0.0677354206563904881d0,0.d0,&
     -0.16933855164097622d0,0.d0,0.620908022683579475d0,1.d0,0.620908022683579475d0,&
     0.d0,-0.16933855164097622d0,0.d0,0.0677354206563904881d0,0.d0,&
     -0.0260520848678424954d0,0.d0,0.00868402828928083181d0,0.d0,&
     -0.00236837135162204504d0,0.d0,0.000501001632073894143d0,0.d0,&
     -0.0000766237790230661631d0,0.d0,7.51213519833981991d-6,0.d0,&
     -3.53757059201598167d-7, &
     0.d0, 0.d0, 0.d0   /)
real(kind=8), dimension(-m:m) ::  cg,cht,cgt

!******** coefficients for wavelet transform *********************
do i=-m,m
   cht(i)=0.d0
   cg(i)=0.d0
   cgt(i)=0.d0
enddo

! the normalization is chosen such that a constant function remains the same constant 
! on each level of the transform

cht( 0)=1.d0

! g coefficients from h coefficients
do i=-m,m-1
   cg(i+1)=cht(-i)*(-1.d0)**(i+1)
   cgt(i+1)=ch(-i)*(-1.d0)**(i+1)
enddo
!!***

