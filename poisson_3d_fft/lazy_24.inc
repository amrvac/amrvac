!!****h* BigDFT/lazy_24
!! NAME
!!   lazy_24
!!
!! FUNCTION
!!   Filters for interpolating scaling functions (order 30)
!!
!! SOURCE
!!
integer, parameter :: m=26
real(kind=8), dimension(-m:m) :: ch = (/ &
     0.d0, 0.d0, 0.d0, &
     -2.00495833269087598d-8,0.d0,5.0505855142546352d-7,0.d0,&
     -6.14044870417274069d-6,0.d0,0.0000480399810385279125d0,0.d0,&
     -0.000272226559218324837d0,0.d0,0.00119360875964957813d0,0.d0,&
     -0.00423188560239395883d0,0.d0,0.0125613112325027032d0,0.d0,&
     -0.032300514597864094d0,0.d0,0.0753678673950162192d0,0.d0,&
     -0.175858357255037845d0,0.d0,0.623497812086043268d0,1.d0,0.623497812086043268d0,&
     0.d0,-0.175858357255037845d0,0.d0,0.0753678673950162192d0,0.d0,&
     -0.032300514597864094d0,0.d0,0.0125613112325027032d0,0.d0,&
     -0.00423188560239395883d0,0.d0,0.00119360875964957813d0,0.d0,&
     -0.000272226559218324837d0,0.d0,0.0000480399810385279125d0,0.d0,&
     -6.14044870417274069d-6,0.d0,5.0505855142546352d-7,0.d0,&
     -2.00495833269087598d-8,&
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

