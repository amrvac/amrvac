!{\src2tex{textfont=tt}}
!!****f* ABINIT/invcb
!! NAME
!! invcb
!!
!! FUNCTION
!! Compute a set of inverse cubic roots as fast as possible :
!! rspts(:)=rhoarr(:)$^\frac{1}{3}$
!!
!! COPYRIGHT
!! Copyright (C) 1998-2006 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors .
!!
!! INPUTS
!!  npts=number of real space points on which density is provided
!!  rhoarr(npts)=input data
!!
!! OUTPUT
!!  rspts(npts)=inverse cubic root of rhoarr
!!
!! PARENTS
!!      calc_lifetime,calc_xc_ep,drivexc,xchcth,xcpbe
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 subroutine invcb(rhoarr,rspts,npts)

 use defs_basis

!This section has been created automatically by the script Abilint (TD). Do not modify these by hand.
#ifdef HAVE_FORTRAN_INTERFACES
 use interfaces_01managempi
#endif
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npts
!arrays
 real(dp),intent(in) :: rhoarr(npts)
 real(dp),intent(out) :: rspts(npts)

!Local variables-------------------------------
!scalars
 integer :: ii,ipts
 real(dp),parameter :: c2_27=2.0d0/27.0d0,c5_9=5.0d0/9.0d0,c8_9=8.0d0/9.0d0
 real(dp),parameter :: m1thrd=-third
 real(dp) :: del,prod,rho,rhom1,rhomtrd
 logical :: test
 character(len=500) :: message

! *************************************************************************

!Loop over points : here, brute force algorithm
!do ipts=1,npts
! rspts(ipts)=rhoarr(ipts)**m1thrd
!end do
!

!write(6,*)' invcb : rhoarr, rspts'

 rhomtrd=rhoarr(1)**m1thrd
 rhom1=1.0d0/rhoarr(1)
 rspts(1)=rhomtrd
 do ipts=2,npts
! write(6,*)
! write(6,*)rhoarr(ipts),rspts(ipts)
  rho=rhoarr(ipts)
  prod=rho*rhom1
! If the previous point is too far ...
  if(prod < 0.01d0 .or. prod > 10.0d0 )then
   rhomtrd=rho**m1thrd
   rhom1=1.0d0/rho
  else
   del=prod-1.0d0
   do ii=1,5
!   Choose one of the two next lines, the last one is more accurate
!   rhomtrd=((1.0d0+third*del)/(1.0d0+two_thirds*del))*rhomtrd
    rhomtrd=((1.0d0+c5_9*del)/(1.0d0+del*(c8_9+c2_27*del)))*rhomtrd
    rhom1=rhomtrd*rhomtrd*rhomtrd
    del=rho*rhom1-1.0d0
!   write(6,*)rhomtrd,del
    test = del*del < 1.0d-24
    if(test) exit
   end do
   if( .not. test) then
    write(message,'(a,a,a,a)' ) ch10,&
&    ' invcb : BUG -',ch10,&
&    '  Fast computation of inverse cubic root failed. '
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
   end if
  end if
  rspts(ipts)=rhomtrd
 end do

 end subroutine invcb
!!***
