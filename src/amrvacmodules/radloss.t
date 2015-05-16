!=============================================================================
subroutine addsource_rloss(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

! w[iw]=w[iw]+qdt*S[wCT,qtC,x] where S is the source based on wCT within ixO

!
! PURPOSE: ADDS THE RADIATIVE LOSSES -dt* rho^2 * Q(T) 
!          TO total ENERGY (ok for classical HD and MHD)
!
!  where
!
! T=T0*p/rho
!
! de/dt=Q0*rho**2*khi(T)*T**alpha(T)
!
!  The constants Q0 and T0 depend on the average particle mass and the density,
!  energy density, and temperature units.
!  The functions khi and alpha depend on the temperature, and different options are 
!  pre-implemented and selected using eqpar(qhow_).
!
! To include these subroutines into the AMRVACUSR module write
!
!INCLUDE: amrvacmodules/radloss.t
!
! then call addsource_rloss from the specialsource subroutine like
!
!if(eqpar(qhow_)>smalldouble)&
!   call addsource_rloss(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)
!
! If there are no other special parameters amrvacusrpar.t might look like this:
!
!INTEGER,PARAMETER :: qcoef_ =neqpar+1   ! radiative loss coefficient Q0
!INTEGER,PARAMETER :: tunit_ =neqpar+2   ! temperature unit T0
!INTEGER,PARAMETER :: qhow_  =neqpar+3   ! integer parameter for the loss function 
!INTEGER,PARAMETER :: Tmin_  =neqpar+4   ! T_min parameter for the loss function 
!INTEGER,PARAMETER :: nspecialpar = 4    
!
! note the precise usage of the qcoef/qunit combination, with slight difference in qcoef
! between different cases of qhow_

! To ensure numerical stability call the getdt_rloss subroutine from getdt_special

! CHANGES
! 01.09.2012: moved to amrvacmodules folder and added x in calling interfaces.
! Oliver Porth

include 'amrvacdef.f'

integer, intent(in) :: ixI^L, ixO^L, iw^LIM
double precision, intent(in) :: qdt, qtC, qt, x(ixI^S,1:ndim)
double precision, intent(inout) :: wCT(ixI^S,1:nw), w(ixI^S,1:nw)
double precision :: tmp(ixG^T), tmp1(ixG^T), tmp2(ixG^T), tmp3(ixG^T)
!-----------------------------------------------------------------------------
oktest=index(teststr,'radloss')>=1
if(oktest)write(*,*)'AddSource_RadLoss, energy:',w(ixtest^D,ee_)

if(e_<1) call mpistop("Radiative losses require energy variable indexed by e_>0!")
if(smallp<zero) call mpistop("Radiative losses require setting of smallp to positive value")

! compute old thermal pressure
call getpthermal(w,x,ixI^L,ixO^L,tmp)
! compute old temperature
tmp1(ixO^S)=tmp(ixO^S)/w(ixO^S,rho_)
! store old kinetic(+magnetic) energy contribution
tmp3(ixO^S)=w(ixO^S,e_)-tmp(ixO^S)/(eqpar(gamma_)-one)

call getQ(wCT,x,ixI^L,ixO^L,tmp2)
! subtract rad. losses from old thermal energy alone
! de/dt = -n_e**2 * Q(T)
! only allow cooling down for temperatures above Tmin
! and store in tmp the new value of the internal energy
where(tmp1(ixO^S)>eqpar(Tmin_))
  tmp(ixO^S)=tmp1(ixO^S)*w(ixO^S,rho_)/(eqpar(gamma_)-one) - qdt*wCT(ixO^S,rho_)**2 *tmp2(ixO^S)
elsewhere
  tmp(ixO^S)=tmp1(ixO^S)*w(ixO^S,rho_)/(eqpar(gamma_)-one)
endwhere

! ensure you never trigger negative energy, or energy subtraction from kinetic
! and magnetic part
where (tmp(ixO^S)>zero)
 w(ixO^S,e_) =tmp3(ixO^S)+tmp(ixO^S)
elsewhere
 w(ixO^S,e_) =tmp3(ixO^S)+smallp/(eqpar(gamma_)-one)
endwhere

if(oktest)write(*,*)'Q, energy:',tmp2(ixtest^D),w(ixtest^D,ee_)

return
end subroutine addsource_rloss
!=============================================================================
subroutine getQ(w,x,ixG^L,ix^L,Q)

! 
! Calculates the optically thin loss function Q according to eqpar(Qhow_).
! Qhow_=1 or Qhow_=2
! (1) See Priest, "Solar MHD", (Reidel, Dordrecht, 1982), p. 88
! (2) See Klimchuk & Gary (Ap. J., 448:925-937,1995)
! (3) cooling curve McDonald & Bailey
! (4) cooling curve McDonald & Bailey added Dalgarno&McCray 1972ARA&A..10.375D for low T
! (5) cooling curve McDonald & Bailey added Dalgarno&McCray 1972ARA&A..10.375D for low T
!     implemented using interpolated values in table
! (6) cooling function calculated by Daria Kosenko, 24-9-2008. Implemented
!    by Klara Schure. using SPEX, solar abundances, starting at log10(T)=3.8 to log10(T)=8.1
!    with steps of delta log10(T)=0.1
!

include 'amrvacdef.f'

integer:: ix^L,ixG^L,i,ix^D
double precision :: w(ixG^T,nw),Q(ixG^T)
double precision, intent(in) :: x(ixI^S,1:ndim)
double precision :: tmp(ixG^T)
double precision :: ion_fr = 1.d-4
double precision :: cpar(110),res1(ixG^T),res2(ixG^T),logt(ixG^T),tind1(ixG^T)
!----------------------------------------------------------------------------
oktest=index(teststr,'getQ')>=1
if(oktest)write(*,*)'GetQ: Q0, T0, Qhow:',&
   eqpar(qcoef_),eqpar(tunit_),eqpar(qhow_)

! Calculate temperature in units given by eqpar(tunit_)
! --> back to dimensionfull temperature in Kelvin
call getpthermal(w,x,ixG^L,ix^L,tmp)
tmp(ix^S) = eqpar(tunit_) * tmp(ix^S) / w(ix^S,rho_)

!c from here uncomment if testing/plotting cooling curve
!c tmp(1:2,ixmin2:ixmax2)=10**0.d0
!c do i=3,1002
!c tmp(i,ixmin2:ixmax2)=10**(log10(tmp(i-1,ixmin2:ixmax2))+0.04)
! enddo
!c to here

! Calculate radiative loss based on the method defined by eqpar(Qhow_)
! --> note: Q-values with factor 10**-40 taken out for compilation
!     purposes
Q(ix^S)=zero
select case(nint(eqpar(Qhow_)))
case(1)
     where (tmp(ix^S)<= 10**3.89063)
        Q(ix^S) = 1.D-40*10**(-2.9) * tmp(ix^S)**(11.7)
     endwhere
     where (tmp(ix^S)>= 10**3.89063.and.tmp(ix^S)<= 10**4.30195)
        Q(ix^S) = 10**(-21.307) * tmp(ix^S)**(6.15)
     endwhere
     where (tmp(ix^S)>= 10**4.30195.and.tmp(ix^S)<= 10**4.575)
        Q(ix^S) = 10**(5.15)
     endwhere
     where (tmp(ix^S)>= 10**4.575.and.tmp(ix^S)<= 10**4.9)
        Q(ix^S) = 10**(-4.0) * tmp(ix^S)**(2.0)
     endwhere
     where (tmp(ix^S)>= 10**4.9.and.tmp(ix^S)<= 10**5.4)
        Q(ix^S) = 10**(5.8)
     endwhere
     where (tmp(ix^S)>= 10**5.4.and.tmp(ix^S)<= 10**5.77)
        Q(ix^S) = 10**(16.6) * tmp(ix^S)**(-2.0)
     endwhere
     where (tmp(ix^S)>= 10**5.77.and.tmp(ix^S)<= 10**6.315)
        Q(ix^S) = 10**(5.06)
     endwhere
     where (tmp(ix^S)>= 10**6.315.and.tmp(ix^S)<= 10**7.60457)
        Q(ix^S) = 10**(9.27) * tmp(ix^S)**(-0.66667)
     endwhere
     where (tmp(ix^S)>= 10**7.60457)
        Q(ix^S) = 10**(0.4) * tmp(ix^S)**(0.5)
     endwhere
case(2)
     where (tmp(ix^S)<= 10**6.18)
        Q(ix^S) = 10**(5.28)
     endwhere
     where (tmp(ix^S)>= 10**6.18.and.tmp(ix^S)<= 10**6.55)
        Q(ix^S) = 10**(14.55) * tmp(ix^S)**(-1.5)
   endwhere
     where (tmp(ix^S)>= 10**6.55)
        Q(ix^S) = 10**(2.54) * tmp(ix^S)**(0.33333)
     endwhere
case(3)
     !from zeus: cooling curve McDonald & Bailey
     where (tmp(ix^S)<= 10**2)
        Q(ix^S) = zero
     endwhere
     where (tmp(ix^S)>= 1.d2 .and. tmp(ix^S)<= 1.d4)
        Q(ix^S) = 2.8347d-10 * (tmp(ix^S)-1.d2)**(2.3562)
     endwhere
     where (tmp(ix^S)>= 1.d4 .and.tmp(ix^S)<= 10**4.05)
        Q(ix^S) = tmp(ix^S)**(-0.133)
     endwhere
     where (tmp(ix^S)> 10**4.05 .and. tmp(ix^S)<= 10**4.15)
        Q(ix^S) = tmp(ix^S)**(0.105)
     endwhere
     where (tmp(ix^S)> 10**4.15.and.tmp(ix^S)<= 10**4.25)
        Q(ix^S) = tmp(ix^S)**(0.452)
     endwhere
     where (tmp(ix^S)> 10**4.25.and.tmp(ix^S)<= 10**4.35)
        Q(ix^S) = tmp(ix^S)**(0.715)
     endwhere
     where (tmp(ix^S)> 10**4.35.and.tmp(ix^S)<= 10**4.45)
        Q(ix^S) = tmp(ix^S)**(0.901)
     endwhere
     where (tmp(ix^S)> 10**4.45.and.tmp(ix^S)<= 10**4.55)
        Q(ix^S) = tmp(ix^S)**(1.030)
     endwhere
     where (tmp(ix^S)> 10**4.55.and.tmp(ix^S)<= 10**4.65)
        Q(ix^S) = tmp(ix^S)**(1.082)
     endwhere
     where (tmp(ix^S)> 10**4.65.and.tmp(ix^S)<= 10**4.75)
        Q(ix^S) = tmp(ix^S)**(1.174)
     endwhere
     where (tmp(ix^S)> 10**4.75.and.tmp(ix^S)<= 10**4.85)
        Q(ix^S) = tmp(ix^S)**(1.257)
     endwhere
     where (tmp(ix^S)> 10**4.85.and.tmp(ix^S)<= 10**4.95)
        Q(ix^S) = tmp(ix^S)**(1.362)
     endwhere
     where (tmp(ix^S)> 10**4.95.and.tmp(ix^S)<= 10**5.05)
        Q(ix^S) = tmp(ix^S)**(1.448)
     endwhere
     where (tmp(ix^S)> 10**5.05.and.tmp(ix^S)<= 10**5.15)
        Q(ix^S) = tmp(ix^S)**(1.523)
     endwhere
     where (tmp(ix^S)> 10**5.15.and.tmp(ix^S)<= 10**5.25)
        Q(ix^S) = tmp(ix^S)**(1.569)
     endwhere
     where (tmp(ix^S)> 10**5.25.and.tmp(ix^S)<= 10**5.35)
        Q(ix^S) = tmp(ix^S)**(1.582)
     endwhere
     where (tmp(ix^S)> 10**5.35.and.tmp(ix^S)<= 10**5.45)
        Q(ix^S) = tmp(ix^S)**(1.539)
     endwhere
     where (tmp(ix^S)> 10**5.45.and.tmp(ix^S)<= 10**5.55)
        Q(ix^S) = tmp(ix^S)**(1.430)
     endwhere
     where (tmp(ix^S)> 10**5.55.and.tmp(ix^S)<= 10**5.65)
        Q(ix^S) = tmp(ix^S)**(1.275)
     endwhere
     where (tmp(ix^S)> 10**5.65.and.tmp(ix^S)<= 10**5.75)
        Q(ix^S) = tmp(ix^S)**(1.168)
     endwhere
     where (tmp(ix^S)> 10**5.75.and.tmp(ix^S)<= 10**5.85)
        Q(ix^S) = tmp(ix^S)**(1.092)
     endwhere
     where (tmp(ix^S)> 10**5.85.and.tmp(ix^S)<= 10**5.95)
        Q(ix^S) = tmp(ix^S)**(1.019)
     endwhere
     where (tmp(ix^S)> 10**5.95.and.tmp(ix^S)<= 10**6.25)
        Q(ix^S) = tmp(ix^S)**(1.000)
     endwhere
     where (tmp(ix^S)> 10**6.25.and.tmp(ix^S)<= 10**6.35)
        Q(ix^S) = tmp(ix^S)**(0.987)
     endwhere
     where (tmp(ix^S)> 10**6.35.and.tmp(ix^S)<= 10**6.45)
        Q(ix^S) = tmp(ix^S)**(0.905)
     endwhere
     where (tmp(ix^S)> 10**6.45.and.tmp(ix^S)<= 10**6.55)
        Q(ix^S) = tmp(ix^S)**(0.738)
     endwhere
     where (tmp(ix^S)> 10**6.55.and.tmp(ix^S)<= 10**6.65)
        Q(ix^S) = tmp(ix^S)**(0.603)
     endwhere
     where (tmp(ix^S)> 10**6.65.and.tmp(ix^S)<= 10**7.05)
        Q(ix^S) = tmp(ix^S)**(0.555)
     endwhere
     where (tmp(ix^S)> 10**7.05.and.tmp(ix^S)<= 10**7.15)
        Q(ix^S) = tmp(ix^S)**(0.535)
     endwhere
     where (tmp(ix^S)> 10**7.15.and.tmp(ix^S)<= 10**7.25)
        Q(ix^S) = tmp(ix^S)**(0.425)
     endwhere
     where (tmp(ix^S)> 10**7.25.and.tmp(ix^S)<= 10**7.35)
        Q(ix^S) = tmp(ix^S)**(0.275)
     endwhere
     where (tmp(ix^S)> 10**7.35.and.tmp(ix^S)<= 10**7.45)
        Q(ix^S) = tmp(ix^S)**(0.251)
     endwhere
     where (tmp(ix^S)> 10**7.45.and.tmp(ix^S)<= 10**7.55)
        Q(ix^S) = tmp(ix^S)**(0.232)
     endwhere
     where (tmp(ix^S)> 10**7.55.and.tmp(ix^S)<= 10**7.65)
        Q(ix^S) = tmp(ix^S)**(0.247)
     endwhere
     where (tmp(ix^S)> 10**7.65.and.tmp(ix^S)<= 10**7.75)
        Q(ix^S) = tmp(ix^S)**(0.283)
     endwhere
     where (tmp(ix^S)> 10**7.75.and.tmp(ix^S)<= 10**7.85)
        Q(ix^S) = tmp(ix^S)**(0.322)
     endwhere
     where (tmp(ix^S)> 10**7.85.and.tmp(ix^S)<= 10**7.9)
        Q(ix^S) = tmp(ix^S)**(0.363)
     endwhere
case(4)
     !from zeus: cooling curve McDonald & Bailey
     !added: from Dalgarno&McCray 1972ARA&A..10.375D for low T
     where (tmp(ix^S)<= 10**2)
        Q(ix^S) = zero
     endwhere
     where (tmp(ix^S)>= 1.d2 .and. tmp(ix^S)<= 3.d3)
        Q(ix^S) = 2.2650d-12
     endwhere
     where (tmp(ix^S)>= 3.d3 .and. tmp(ix^S)<= 4.d3)
        Q(ix^S) = 4.2735d-8
     endwhere
     where (tmp(ix^S)>= 4.d3 .and. tmp(ix^S)<= 5.d3)
        Q(ix^S) = 1.5812d-5
     endwhere
     where (tmp(ix^S)>= 5.d3 .and. tmp(ix^S)<= 5.5d3)
        Q(ix^S) = 1.3675d-4
     endwhere
     where (tmp(ix^S)>= 5.5d3 .and. tmp(ix^S)<= 6.d3)
        Q(ix^S) = 8.1197d-4
     endwhere
     where (tmp(ix^S)>= 6.d3 .and. tmp(ix^S)<= 6.5d3)
        Q(ix^S) = 3.7607d-3
     endwhere
     where (tmp(ix^S)>= 6.5d3 .and. tmp(ix^S)<= 7.d3)
        Q(ix^S) = 1.3675d-3
     endwhere
     where (tmp(ix^S)>= 7.d3 .and. tmp(ix^S)<= 7.5d3)
        Q(ix^S) = 4.2735d-2
     endwhere
     where (tmp(ix^S)>= 7.5d3 .and. tmp(ix^S)<= 8.d3)
        Q(ix^S) = 1.1538d-1
     endwhere
     where (tmp(ix^S)>= 8.d3 .and. tmp(ix^S)<= 8.5d3)
        Q(ix^S) = 2.7778d-1
     endwhere
     where (tmp(ix^S)>= 8.5d3 .and. tmp(ix^S)<= 9.d3)
        Q(ix^S) = 5.9829d-1
     endwhere
     where (tmp(ix^S)>= 9.d3 .and. tmp(ix^S)<= 9.5d3)
        Q(ix^S) = 1.2393d0
     endwhere
     where (tmp(ix^S)>= 9.5d3 .and. tmp(ix^S)<= 10.d3)
        Q(ix^S) = 2.3077d0
     endwhere
     where (tmp(ix^S)>= 10.d3 .and. tmp(ix^S)<= 10.5d3)
        Q(ix^S) = 4.0598d0
     endwhere
     where (tmp(ix^S)>= 10.5d3 .and. tmp(ix^S)<= 11.d3)
        Q(ix^S) = 6.7949d0
     endwhere
     where (tmp(ix^S)>= 11.d3 .and. tmp(ix^S)<= 11.5d3)
        Q(ix^S) = 1.0940d0
     endwhere
     where (tmp(ix^S)> 1.15d4 .and. tmp(ix^S)<= 10**4.15)
        Q(ix^S) = tmp(ix^S)**(0.105)
     endwhere
     where (tmp(ix^S)> 10**4.15.and.tmp(ix^S)<= 10**4.25)
        Q(ix^S) = tmp(ix^S)**(0.452)
     endwhere
     where (tmp(ix^S)> 10**4.25.and.tmp(ix^S)<= 10**4.35)
        Q(ix^S) = tmp(ix^S)**(0.715)
     endwhere
     where (tmp(ix^S)> 10**4.35.and.tmp(ix^S)<= 10**4.45)
        Q(ix^S) = tmp(ix^S)**(0.901)
     endwhere
     where (tmp(ix^S)> 10**4.45.and.tmp(ix^S)<= 10**4.55)
        Q(ix^S) = tmp(ix^S)**(1.030)
     endwhere
     where (tmp(ix^S)> 10**4.55.and.tmp(ix^S)<= 10**4.65)
        Q(ix^S) = tmp(ix^S)**(1.082)
     endwhere
     where (tmp(ix^S)> 10**4.65.and.tmp(ix^S)<= 10**4.75)
        Q(ix^S) = tmp(ix^S)**(1.174)
     endwhere
     where (tmp(ix^S)> 10**4.75.and.tmp(ix^S)<= 10**4.85)
        Q(ix^S) = tmp(ix^S)**(1.257)
     endwhere
     where (tmp(ix^S)> 10**4.85.and.tmp(ix^S)<= 10**4.95)
        Q(ix^S) = tmp(ix^S)**(1.362)
     endwhere
     where (tmp(ix^S)> 10**4.95.and.tmp(ix^S)<= 10**5.05)
        Q(ix^S) = tmp(ix^S)**(1.448)
     endwhere
     where (tmp(ix^S)> 10**5.05.and.tmp(ix^S)<= 10**5.15)
        Q(ix^S) = tmp(ix^S)**(1.523)
     endwhere
     where (tmp(ix^S)> 10**5.15.and.tmp(ix^S)<= 10**5.25)
        Q(ix^S) = tmp(ix^S)**(1.569)
     endwhere
     where (tmp(ix^S)> 10**5.25.and.tmp(ix^S)<= 10**5.35)
        Q(ix^S) = tmp(ix^S)**(1.582)
     endwhere
     where (tmp(ix^S)> 10**5.35.and.tmp(ix^S)<= 10**5.45)
        Q(ix^S) = tmp(ix^S)**(1.539)
     endwhere
     where (tmp(ix^S)> 10**5.45.and.tmp(ix^S)<= 10**5.55)
        Q(ix^S) = tmp(ix^S)**(1.430)
     endwhere
     where (tmp(ix^S)> 10**5.55.and.tmp(ix^S)<= 10**5.65)
        Q(ix^S) = tmp(ix^S)**(1.275)
     endwhere
     where (tmp(ix^S)> 10**5.65.and.tmp(ix^S)<= 10**5.75)
        Q(ix^S) = tmp(ix^S)**(1.168)
     endwhere
     where (tmp(ix^S)> 10**5.75.and.tmp(ix^S)<= 10**5.85)
        Q(ix^S) = tmp(ix^S)**(1.092)
     endwhere
     where (tmp(ix^S)> 10**5.85.and.tmp(ix^S)<= 10**5.95)
        Q(ix^S) = tmp(ix^S)**(1.019)
     endwhere
     where (tmp(ix^S)> 10**5.95.and.tmp(ix^S)<= 10**6.25)
        Q(ix^S) = tmp(ix^S)**(1.000)
     endwhere
     where (tmp(ix^S)> 10**6.25.and.tmp(ix^S)<= 10**6.35)
        Q(ix^S) = tmp(ix^S)**(0.987)
     endwhere
     where (tmp(ix^S)> 10**6.35.and.tmp(ix^S)<= 10**6.45)
        Q(ix^S) = tmp(ix^S)**(0.905)
     endwhere
     where (tmp(ix^S)> 10**6.45.and.tmp(ix^S)<= 10**6.55)
        Q(ix^S) = tmp(ix^S)**(0.738)
     endwhere
     where (tmp(ix^S)> 10**6.55.and.tmp(ix^S)<= 10**6.65)
        Q(ix^S) = tmp(ix^S)**(0.603)
     endwhere
     where (tmp(ix^S)> 10**6.65.and.tmp(ix^S)<= 10**7.05)
        Q(ix^S) = tmp(ix^S)**(0.555)
     endwhere
     where (tmp(ix^S)> 10**7.05.and.tmp(ix^S)<= 10**7.15)
        Q(ix^S) = tmp(ix^S)**(0.535)
     endwhere
     where (tmp(ix^S)> 10**7.15.and.tmp(ix^S)<= 10**7.25)
        Q(ix^S) = tmp(ix^S)**(0.425)
     endwhere
     where (tmp(ix^S)> 10**7.25.and.tmp(ix^S)<= 10**7.35)
        Q(ix^S) = tmp(ix^S)**(0.275)
     endwhere
     where (tmp(ix^S)> 10**7.35.and.tmp(ix^S)<= 10**7.45)
        Q(ix^S) = tmp(ix^S)**(0.251)
     endwhere
     where (tmp(ix^S)> 10**7.45.and.tmp(ix^S)<= 10**7.55)
        Q(ix^S) = tmp(ix^S)**(0.232)
     endwhere
     where (tmp(ix^S)> 10**7.55.and.tmp(ix^S)<= 10**7.65)
        Q(ix^S) = tmp(ix^S)**(0.247)
     endwhere
     where (tmp(ix^S)> 10**7.65.and.tmp(ix^S)<= 10**7.75)
        Q(ix^S) = tmp(ix^S)**(0.283)
     endwhere
     where (tmp(ix^S)> 10**7.75.and.tmp(ix^S)<= 10**7.85)
        Q(ix^S) = tmp(ix^S)**(0.322)
     endwhere
     where (tmp(ix^S)> 10**7.85.and.tmp(ix^S)<= 10**7.9)
        Q(ix^S) = tmp(ix^S)**(0.363)
     endwhere
case(5)
  cpar(1:41)=(/ -0.133, 0.105, 0.452, 0.715, 0.901 &
               , 1.030, 1.082, 1.174, 1.257, 1.362 &
               , 1.448, 1.523, 1.569, 1.582, 1.539 &
               , 1.430, 1.275, 1.168, 1.092, 1.019 &
               , 1.000, 1.004, 1.008, 0.987, 0.905 &
               , 0.738, 0.603, 0.555, 0.552, 0.554 &
               , 0.552, 0.535, 0.425, 0.275, 0.251 &
               , 0.232, 0.247, 0.283, 0.322, 0.363, 0.397 /)

{^D&do ix^D=ixmin^D,ixmax^D\}
if (tmp(ix^D)<=1.d4) then
  !abundances used:
  !Oxygen: 4.4d-4
  !Carbon: 3.75d-4
  !Nitrogen: 8.7d-5
  !Silicon: 3.2d-5
  !Iron: 3.2d-5
  !Neon: 2.6d-5
  !Sulphur: 1.4d-5

  !Carbon:
  Q(ix^D) = Q(ix^D)+ion_fr*3.75d-4*(7.9d-20*tmp(ix^D)**(-0.5d0)*exp(-92.d0/tmp(ix^D)))
  !Silicon
  Q(ix^D) = Q(ix^D)+ion_fr*3.2d-5*(1.9d-18*tmp(ix^D)**(-0.5d0)*exp(-4.13d2/tmp(ix^D)))
  Q(ix^D) = Q(ix^D)+3.2d-5*(7.4d-23*exp(-4.13d2/tmp(ix^D)))
  !Iron
  Q(ix^D) = Q(ix^D)+ion_fr*3.2d-5*(1.1d-18*tmp(ix^D)**(-0.5d0)*&
             (exp(-5.54d2/tmp(ix^D))+1.3*exp(-9.61d2/tmp(ix^D))))
  Q(ix^D) = Q(ix^D)+3.2d-5*1.1d-22*(exp(-5.54d2/tmp(ix^D))+1.4*exp(-9.41d2/tmp(ix^D)))
  !Oxygen
  Q(ix^D) = Q(ix^D)+ion_fr*4.4d-4*(1.74d-24*tmp(ix^D)**(0.5d0)*&
             ((1.d0-7.6*tmp(ix^D)**(-0.5d0))*exp(-2.28d2/tmp(ix^D))+&
             0.38*(1.d0-7.7*tmp(ix^D)**(-0.5d0))*exp(-3.26d2/tmp(ix^D))))
  !Oxygen
  Q(ix^D) = Q(ix^D)+ion_fr*3.75d-4*(9.4d-23*tmp(ix^D)**(0.5d0)*exp(-2.27d4/tmp(ix^D)))
  !Nitrogen
  Q(ix^D) = Q(ix^D)+ion_fr*8.7d-5*(8.2d-22*tmp(ix^D)**(0.5d0)*&
             (1.d0-(tmp(ix^D)**2)*2.7d-9)*exp(-2.77d4/tmp(ix^D)))
elseif (tmp(ix^D)>=10**(7.9d0)) then
    Q(ix^D) = 2.6787d-27*tmp(ix^D)**(0.5d0)
else
     logt(ix^D) = log10(tmp(ix^D))
     res1(ix^D) = (1.-(10.*logt(ix^D)-int(10.*logt(ix^D))))
     res2(ix^D) = (1.-(int(10.*logt(ix^D))+1.-10.*logt(ix^D)))
     Q(ix^D) = 1.1167d-23*10**(res1(ix^D)*cpar(int(10.*logt(ix^D))-38)+&
                    res2(ix^D)*cpar(int(10.*logt(ix^D))-37))
endif
{^D&end do\}

case(6)
  ! cooling function calculated by Daria Kosenko, 24-9-2008. Implemented
  ! by Klara Schure.
  ! using SPEX, solar abundances, starting at log10(T)=3.8 to log10(T)=8.1
  ! with steps of delta log10(T)=0.04
  cpar(1:110)=(/ -30.6104, -29.4107, -28.4601, -27.5743&
               , -26.3765, -25.2890, -24.2684, -23.3834&
               , -22.5977, -21.9689, -21.5972, -21.4615&
               , -21.4788, -21.5485, -21.6143, -21.6306&
               , -21.5614, -21.4282, -21.2901, -21.1706&
               , -21.0649, -20.9686, -20.8813, -20.8043&
               , -20.7373, -20.6824, -20.6460, -20.6285&
               , -20.6265, -20.6360, -20.6588, -20.6913&
               , -20.7145, -20.7158, -20.7028, -20.6878&
               , -20.6786, -20.6742, -20.6704, -20.6744&
               , -20.7086, -20.8029, -20.9644, -21.1480&
               , -21.2929, -21.3764, -21.4126, -21.4286&
               , -21.4529, -21.5037, -21.5706, -21.6248&
               , -21.6549, -21.6696, -21.6822, -21.7020&
               , -21.7265, -21.7467, -21.7593, -21.7694&
               , -21.7874, -21.8242, -21.8874, -21.9736&
               , -22.0666, -22.1526, -22.2248, -22.2800&
               , -22.3190, -22.3438, -22.3566, -22.3603&
               , -22.3572, -22.3496, -22.3406, -22.3330&
               , -22.3300, -22.3335, -22.3436, -22.3586&
               , -22.3772, -22.3999, -22.4282, -22.4619&
               , -22.4990, -22.5349, -22.5656, -22.5892&
               , -22.6057, -22.6160, -22.6207, -22.6213&
               , -22.6183, -22.6126, -22.6045, -22.5945&
               , -22.5831, -22.5707, -22.5573, -22.5434&
               , -22.5287, -22.5140, -22.4992, -22.4844&
               , -22.4695, -22.4543, -22.4392, -22.4237&
               , -22.4087, -22.3927 /)

{^D&do ix^D=ixmin^D,ixmax^D\}
Q(ix^D)=zero
if (tmp(ix^D)<=1.d4) then
  !abundances used:
  !Oxygen: 8.5114d-4  !4.4d-4
  !Carbon: 3.6321d-4  !3.75d-4
  !Nitrogen: 1.1220d-4 !8.7d-5
  !Silicon: 3.5481d-5 !3.2d-5
  !Iron: 4.6774d-5  !3.2d-5
  !Neon: 1.2303d-4  !2.6d-5
  !Sulphur: 1.6218d-5  !1.4d-5
  !Carbon:
  Q(ix^D) = Q(ix^D)+ion_fr*3.6321d-4*(7.9d-20*tmp(ix^D)**(-0.5d0)*exp(-92.d0/tmp(ix^D)))
  !Silicon
  Q(ix^D) = Q(ix^D)+ion_fr*3.5481d-5*(1.9d-18*tmp(ix^D)**(-0.5d0)*exp(-4.13d2/tmp(ix^D)))
  Q(ix^D) = Q(ix^D)+3.5481d-5*(7.4d-23*exp(-4.13d2/tmp(ix^D)))
  !Iron
  Q(ix^D) = Q(ix^D)+ion_fr*4.6774d-5*(1.1d-18*tmp(ix^D)**(-0.5d0)*&
             (exp(-5.54d2/tmp(ix^D))+1.3*exp(-9.61d2/tmp(ix^D))))
  Q(ix^D) = Q(ix^D)+4.6774d-5*1.1d-22*(exp(-5.54d2/tmp(ix^D))+1.4*exp(-9.41d2/tmp(ix^D)))
  !Oxygen
  Q(ix^D) = Q(ix^D)+ion_fr*8.5114d-4*(1.74d-24*tmp(ix^D)**(0.5d0)*&
             ((1.d0-7.6*tmp(ix^D)**(-0.5d0))*exp(-2.28d2/tmp(ix^D))+&
             0.38*(1.d0-7.7*tmp(ix^D)**(-0.5d0))*exp(-3.26d2/tmp(ix^D))))
  Q(ix^D) = Q(ix^D)+ion_fr*8.5114d-4*(9.4d-23*tmp(ix^D)**(0.5d0)*exp(-2.27d4/tmp(ix^D)))
  !Nitrogen
  Q(ix^D) = Q(ix^D)+ion_fr*1.1220d-4*(8.2d-22*tmp(ix^D)**(0.5d0)*&
             (1.d0-(tmp(ix^D)**2)*2.7d-9)*exp(-2.77d4/tmp(ix^D)))
elseif (tmp(ix^D)>=10**8.1) then
    Q(ix^D) = 2.6787d-27*tmp(ix^D)**(0.5d0)
else  !from SPEX cooling curve
     logt(ix^D) = log10(tmp(ix^D))
     tind1(ix^D) = (10.0d0*logt(ix^D)-37.6d0)*2.5d0
     res1(ix^D) = (1.0d0-(tind1(ix^D)-int(tind1(ix^D))))
     res2(ix^D) = 1.0d0-res1(ix^D)
     Q(ix^D) =  10.0d0**(res1(ix^D)*cpar(int(tind1(ix^D)))+&
                         res2(ix^D)*cpar(int(tind1(ix^D))+1))
endif
{^D&end do\}

case default
   write(*,*)'Error in GetQ: Unknown value for eqpar(qhow_):',eqpar(qhow_)
end select

!renormalize Q back to used m/l/t^3 (erg/cm3/s) unit
if ((nint(eqpar(Qhow_)) == 1).or.(nint(eqpar(qhow_)) == 2)) then
   Q(ix^S) = Q(ix^S) * eqpar(qcoef_) * 1.d-40
else
   Q(ix^S) = Q(ix^S) * eqpar(qcoef_)
endif

return
end subroutine getQ
!==============================================================================
subroutine getdt_rloss(w,ixG^L,ix^L,dtnew,dx^D,x)

!
! Check diffusion time limit for dt
!

include 'amrvacdef.f'

integer, intent(in)             :: ixG^L, ix^L
double precision, intent(in)    :: dx^D, x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw), dtnew

double precision :: tmp2(ixG^T)
!----------------------------------------------------------------------------
oktest=index(teststr,'getdt')>=1

! Calculate Q
call getQ(w,x,ixG^L,ix^L,tmp2)

! dt< e/(Q*rho**2)
dtnew = dtdiffpar/maxval(w(ix^S,rho_)**2*tmp2(ix^S)/w(ix^S,e_))

if (oktest)write(*,*)'Thin radiative losses dt:',dtnew

return
end subroutine getdt_rloss
!=============================================================================
