!=============================================================================
! Source terms for synchrotron cooling according to Camus et al. 2009
! http://adsabs.harvard.edu/abs/2009MNRAS.400.1241C
! 
! This is for inclusion in the userfile, activate with 
! #define EPSINF in definitions.h
! Currently available for smrhd.  
! Can be adopted for srhdeos where we also have the #define EPSINF 
! You would then need to express the co-moving field in terms of other flow quantities (assuming equipartition or the like).
! Written by Oliver Porth
! Version September 2012.
!=============================================================================
subroutine specialsource(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

! Calculate w(iw)=w(iw)+qdt*SOURCE[wCT,qtC,x] within ixO for all indices
! iw=iwmin...iwmax.  wCT is at time qCT

use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L, iw^LIM
double precision, intent(in) :: qdt, qtC, qt, x(ixI^S,1:ndim)
double precision, intent(inout) :: wCT(ixI^S,1:nw), w(ixI^S,1:nw)
! .. local ..
logical          :: patchw(ixG^T)
double precision :: c2tilde, epstocode
double precision :: bdash2(ixG^T), e2(ixG^T), b2(ixG^T), epsinf(ixG^T)
!-----------------------------------------------------------------------------
patchw = .true.
patchw(ixO^S) = .false.

if(.not. useprimitiveRel)then
   call mpistop('only implemented for useprimitiveRel=T')
endif

! Primitive call because we need velocities, Lorentz factor and primitive epsinf:
call primitive(ixI^L,ixO^L,w,x)
! Check first if we need to cool at all:
if (.not.any(w(ixO^S,epsinf_)>eqpar(epsfloor_))) then
   call conserven(ixI^L,ixO^L,w,patchw)
   return
end if

! Synchrotron constant c2tilde = 2/3*c2
c2tilde = 1.578d-3     ! g^-2 cm^-1 s^3 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Bring to code units:
c2tilde   = c2tilde   / (UNIT_DENSITY**(-2.0D0)*UNIT_LENGTH**(-4.0D0)*UNIT_VELOCITY**(-3.0D0)) 
epstocode = w_convert_factor(epsinf_)  / (UNIT_DENSITY*UNIT_LENGTH**3.0D0*UNIT_VELOCITY**2.0D0) * 4*dpi


! Obtain co-moving magnetic field:
! using the invariant Bdash^2 = B^2-E^2
b2(ixO^S)     = {^C&w(ixO^S,b^C_)**2.0d0|+}
e2(ixO^S) = (w(ixO^S,b2_)*w(ixO^S,u3_)/w(ixO^S,lfac_) - w(ixO^S,b3_)*w(ixO^S,u2_)/w(ixO^S,lfac_))**2.0d0 &
     + (w(ixO^S,b3_)*w(ixO^S,u1_)/w(ixO^S,lfac_) - w(ixO^S,b1_)*w(ixO^S,u3_)/w(ixO^S,lfac_))**2.0d0 &
     + (w(ixO^S,b1_)*w(ixO^S,u2_)/w(ixO^S,lfac_) - w(ixO^S,b2_)*w(ixO^S,u1_)/w(ixO^S,lfac_))**2.0d0
bdash2(ixO^S) = b2(ixO^S) - e2(ixO^S)

! Store the primitive cutoff energy:
epsinf(ixO^S) = w(ixO^S,epsinf_)

call conserven(ixI^L,ixO^L,w,patchw)

where (epsinf(ixO^S)>eqpar(epsfloor_))
   w(ixO^S,Depsinf_) = w(ixO^S,Depsinf_) - qdt*c2tilde &
        * abs(epsinf(ixO^S))**2.0d0*bdash2(ixO^S)*w(ixO^S,rho1_)**(2.0d0/3.0d0) &
        * epstocode
end where

! reset to floor: 
where (w(ixO^S,Depsinf_)<w(ixO^S,lfac_)**(1.0d0/3.0d0)*eqpar(epsfloor_)*w(ixO^S,rho1_)**(2.0d0/3.0d0))
   w(ixO^S,Depsinf_) = w(ixO^S,lfac_)**(1.0d0/3.0d0)*eqpar(epsfloor_)*w(ixO^S,rho1_)**(2.0d0/3.0d0)
end where

end subroutine specialsource
!=============================================================================
subroutine getdt_special(w,ixG^L,ix^L,dtnew,dx^D,x)

! Limit "dt" further if necessary, e.g. due to the special source terms.
! The getdt_courant (CFL condition) and the getdt subroutine in the AMRVACPHYS
! module have already been called.

use mod_global_parameters

integer, intent(in) :: ixG^L, ix^L
double precision, intent(in) :: dx^D, x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw), dtnew
!.. local ..
logical          :: patchw(ixG^T)
double precision :: epstocode, epsinf(ixG^T)
double precision :: bdash2(ixG^T), e2(ixG^T), b2(ixG^T)
double precision :: c2tilde
!
double precision, parameter :: dtfac=0.3d0
!-----------------------------------------------------------------------------

if(.not. useprimitiveRel)then
   call mpistop('only implemented for useprimitiveRel=T')
endif

! Check first if we need to cool further:
if (.not.any(w(ix^S,Depsinf_)>w(ix^S,lfac_)**(1.0d0/3.0d0)*eqpar(epsfloor_)*w(ix^S,Drho1_)**(2.0d0/3.0d0))) then 
   dtnew=bigdouble
   return
end if

patchw = .true.
patchw(ix^S) = .false.

! Synchrotron constant c2tilde = 2/3*c2
c2tilde = 1.578d-3     ! g^-2 cm^-1 s^3 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Bring to code units:
c2tilde   = c2tilde   / (UNIT_DENSITY**(-2.0D0)*UNIT_LENGTH**(-4.0D0)*UNIT_VELOCITY**(-3.0D0)) 
epstocode = w_convert_factor(epsinf_)  / (UNIT_DENSITY*UNIT_LENGTH**3.0D0*UNIT_VELOCITY**2.0D0) * 4*dpi

call primitiven(ixG^L,ix^L,w,patchw)

! Obtain co-moving magnetic field:
! using the invariant Bdash^2 = B^2-E^2
b2(ix^S)     = {^C&w(ix^S,b^C_)**2.0d0|+}
e2(ix^S) = (w(ix^S,b2_)*w(ix^S,u3_)/w(ix^S,lfac_) - w(ix^S,b3_)*w(ix^S,u2_)/w(ix^S,lfac_))**2.0d0 &
     + (w(ix^S,b3_)*w(ix^S,u1_)/w(ix^S,lfac_) - w(ix^S,b1_)*w(ix^S,u3_)/w(ix^S,lfac_))**2.0d0 &
     + (w(ix^S,b1_)*w(ix^S,u2_)/w(ix^S,lfac_) - w(ix^S,b2_)*w(ix^S,u1_)/w(ix^S,lfac_))**2.0d0
bdash2(ix^S) = b2(ix^S) - e2(ix^S)

! Store the primitive cutoff energy:
epsinf(ix^S) = w(ix^S,epsinf_)
call conserven(ixG^L,ix^L,w,patchw)


dtnew=dtfac/c2tilde/epstocode * min(abs(minval(w(ix^S,lfac_)/epsinf(ix^S)/bdash2(ix^S))),bigdouble)

end subroutine getdt_special
!=============================================================================
