!=============================================================================
program xodeint

use mod_odeint
implicit none

external derivs

integer                          :: nok, nbad
double precision                 :: x1, x2, h1
double precision,parameter       :: eps=1.0d-10, hmin=1.0d-12
integer, parameter               :: nvar = 2
double precision                 :: ystart(nvar)
!-----------------------------------------------------------------------------

ystart = 1.0d0
x1     = 1.0d0
x2     = 2.0d0
h1     = 1.0d-6


call odeint(ystart,nvar,x1,x2,eps,h1,hmin,nok,nbad,derivs,rkqs)

print*, ystart, nok, nbad

end program xodeint
!=============================================================================
subroutine derivs(x,y,dydx)

double precision, intent(in)                :: x, y(*)
double precision, intent(out)               :: dydx(*)
!-----------------------------------------------------------------------------

dydx(1) = y(1)
dydx(2) = y(2)

end subroutine derivs
!=============================================================================
