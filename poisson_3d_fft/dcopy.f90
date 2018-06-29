! copies a vector, x, to a vector, y.
! uses unrolled loops for increments equal to one.
! jack dongarra, linpack, 3/11/78.
! modified 12/3/93, array(1) declarations changed to array(*)
! modified 2018 (Jannis Teunissen): translate to Fortran 90
subroutine dcopy(n,dx,incx,dy,incy)
  double precision dx(*),dy(*)
  integer i,incx,incy,ix,iy,m,mp1,n

  if (n < 0) return
  if (incx == 1 .and. incy == 1) go to 20

  ! code for unequal increments or equal increments
  ! not equal to 1
  ix = 1
  iy = 1
  if (incx < 0) ix = (-n+1)*incx + 1
  if (incy < 0) iy = (-n+1)*incy + 1
  do i = 1, n
     dy(iy) = dx(ix)
     ix     = ix + incx
     iy     = iy + incy
  end do
  return

  ! code for both increments equal to 1
  !
  ! clean-up loop

20 m = mod(n,7)
  if ( m == 0 ) go to 40
  do i = 1,m
     dy(i) = dx(i)
  end do

  if ( n  <  7 ) return

40 mp1 = m + 1
  do i = mp1,n,7
     dy(i)     = dx(i)
     dy(i + 1) = dx(i + 1)
     dy(i + 2) = dx(i + 2)
     dy(i + 3) = dx(i + 3)
     dy(i + 4) = dx(i + 4)
     dy(i + 5) = dx(i + 5)
     dy(i + 6) = dx(i + 6)
  end do

end subroutine dcopy
