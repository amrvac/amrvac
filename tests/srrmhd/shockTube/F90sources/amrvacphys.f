!#############################################################################
! Module amrvacphys/- srrmhd
! Special relativistic resistive MHD as in Komissarov 2007
! 2015-30-07 by Oliver Porth

!=============================================================================
subroutine getdt(w,ixGmin1,ixGmax1,ixmin1,ixmax1,dtnew,dx1,x)

include 'amrvacdef.f'

integer, intent(in) :: ixGmin1,ixGmax1, ixmin1,ixmax1
double precision, intent(in) :: dx1, x(ixGmin1:ixGmax1,1:ndim)
double precision, intent(inout) :: w(ixGmin1:ixGmax1,1:nw), dtnew
!-----------------------------------------------------------------------------

dtnew=bigdouble

end subroutine getdt
!=============================================================================
!##############################################################################
! module amrvacnul.roe.t - nul-routines  
!=============================================================================
subroutine average(wL,wR,x,ixmin1,ixmax1,idims,wroe,workroe)

include 'amrvacdef.f'

integer:: ixmin1,ixmax1,idims,idir
double precision, dimension(ixGlo1:ixGhi1,nw):: wL,wR,wroe
double precision, intent(in)      :: x(ixGlo1:ixGhi1,1:ndim)
double precision, dimension(ixGlo1:ixGhi1,nworkroe):: workroe
!-----------------------------------------------------------------------------

call mpistop("Error:average: no roe solver implemented!")

end subroutine average
!=============================================================================
subroutine geteigenjump(wL,wR,wroe,x,ixmin1,ixmax1,il,idims,smalla,a,jump,&
   workroe)

include 'amrvacdef.f'

integer:: ixmin1,ixmax1,il,idims
double precision, dimension(ixGlo1:ixGhi1,nw):: wL,wR,wroe
double precision, intent(in)         :: x(ixGlo1:ixGhi1,1:ndim)
double precision, dimension(ixGlo1:ixGhi1)   :: smalla,a,jump
double precision, dimension(ixGlo1:ixGhi1,nworkroe) :: workroe
!-----------------------------------------------------------------------------

call mpistop("Error:geteigenjump: no roe solver implemented!")

end subroutine geteigenjump
!=============================================================================
subroutine rtimes(q,wroe,ixmin1,ixmax1,iw,il,idims,rq,workroe)

include 'amrvacdef.f'

integer::          ixmin1,ixmax1,iw,il,idims
double precision:: wroe(ixGlo1:ixGhi1,nw)
double precision, dimension(ixGlo1:ixGhi1):: q,rq
double precision, dimension(ixGlo1:ixGhi1,nworkroe):: workroe
!-----------------------------------------------------------------------------

call mpistop("Error:rtimes: no roe solver implemented!")

end subroutine rtimes
!=============================================================================
! end module amrvacnul.roe.t
!##############################################################################
!##############################################################################
! module amrvacnul.hllc.t - nul-routines  
!=============================================================================
subroutine diffuse_hllcd(ixImin1,ixImax1,ixOmin1,ixOmax1,idims,wLC,wRC,fLC,&
   fRC,patchf)

! when method is hllcd or hllcd1 then: 

! this subroutine is to impose enforce regions where we AVOID HLLC
! and use TVDLF instead: this is achieved by setting patchf to 4 in
! certain regions. An additional input parameter is nxdiffusehllc
! which sets the size of the fallback region.

! This nul version enforces TVDLF everywhere!!!

include 'amrvacdef.f'

integer, intent(in)                                      :: ixImin1,ixImax1,&
   ixOmin1,ixOmax1,idims
double precision, dimension(ixGlo1:ixGhi1,1:nw), intent(in)      :: wRC,wLC
double precision, dimension(ixGlo1:ixGhi1,1:nwflux),intent(in) :: fLC, fRC

integer         , dimension(ixGlo1:ixGhi1), intent(inout)        :: patchf
!-----------------------------------

! enforce TVDLF everywhere
patchf(ixOmin1:ixOmax1) = 4

end subroutine diffuse_hllcd
!=============================================================================
subroutine getlCD(wLC,wRC,fLC,fRC,cmin,cmax,idims,ixImin1,ixImax1,ixOmin1,&
   ixOmax1, whll,Fhll,lambdaCD,patchf)

! Calculate lambda at CD and set the patchf to know the orientation
! of the riemann fan and decide on the flux choice
! We also compute here the HLL flux and w value, for fallback strategy

! In this nul version, we simply compute nothing and ensure TVDLF fallback

include 'amrvacdef.f'

integer, intent(in)                                      :: ixImin1,ixImax1,&
   ixOmin1,ixOmax1,idims
double precision, dimension(ixGlo1:ixGhi1,1:nw), intent(in)      :: wLC,wRC
double precision, dimension(ixGlo1:ixGhi1,1:nwflux), intent(in):: fLC,fRC
double precision, dimension(ixGlo1:ixGhi1), intent(in)           :: cmax,cmin

integer         , dimension(ixGlo1:ixGhi1), intent(inout)        :: patchf

double precision, dimension(ixGlo1:ixGhi1,1:nwflux), intent(out) :: Fhll,whll
double precision, dimension(ixGlo1:ixGhi1), intent(out)            :: lambdaCD
!--------------------------------------------

! Next must normally be computed
Fhll(ixOmin1:ixOmax1,1:nwflux)=zero
whll(ixOmin1:ixOmax1,1:nwflux)=zero
lambdaCD(ixOmin1:ixOmax1)=zero

! this actually ensures fallback to TVDLF
patchf(ixOmin1:ixOmax1)=4

return
end subroutine getlCD
!=============================================================================
subroutine getwCD(wLC,wRC,whll,vLC,vRC,fRC,fLC,Fhll,patchf,lambdaCD,cmin,cmax,&
   ixImin1,ixImax1,ixOmin1,ixOmax1,idims,f)

! compute the intermediate state U*
! only needed where patchf=-1/1

! This nul version simply nullifies all values

include 'amrvacdef.f'

integer, intent(in)                                      :: ixImin1,ixImax1,&
   ixOmin1,ixOmax1,idims
double precision, dimension(ixGlo1:ixGhi1,1:nw), intent(in)      :: wRC,wLC
double precision, dimension(ixGlo1:ixGhi1,1:nwflux), intent(in):: whll, Fhll
double precision, dimension(ixGlo1:ixGhi1), intent(in)           :: vRC, vLC,&
   lambdaCD
double precision, dimension(ixGlo1:ixGhi1), intent(in)           :: cmax,cmin
double precision, dimension(ixGlo1:ixGhi1,1:nwflux), intent(in):: fRC,fLC

integer         , dimension(ixGlo1:ixGhi1), intent(in)           :: patchf

double precision, dimension(ixGlo1:ixGhi1,1:nwflux),intent(out) :: f
!--------------------------------------------

! Next must normally be computed
f(ixOmin1:ixOmax1,1:nwflux)  = zero

return
end subroutine getwCD
!=============================================================================
! end module amrvacnul.hllc.t
!##############################################################################
!=============================================================================
subroutine checkglobaldata
include 'amrvacdef.f'
!-----------------------------------------------------------------------------
minrho = max(zero,smallrho)


govergminone =eqpar(gamma_)/(eqpar(gamma_)-one)
minp  = max(zero,smallp)
smallxi=minrho+minp*govergminone
smalltau = minp/(eqpar(gamma_)-one)



! Check if ssplitdivb = .true.
if (ssplitdivb .ne. .true.) call mpistop&
   ('Please run with ssplitdivb = .true.')

! We require three vector components
if (3 .ne. 3) call mpistop('SRRMHD needs three vector components.')

end subroutine checkglobaldata
!=============================================================================
subroutine initglobaldata

! place to set entropy fixes etc, absent for now

include 'amrvacdef.f'
!-----------------------------------------------------------------------------
eqpar(adiab_)=1.0d0


eqpar(gamma_)=4.0d0/3.0d0



eqpar(kappa_)= 0.5d0

if(strictzero)limitvalue=zero
if(.not.strictzero)limitvalue=smalldouble**2.0d0

end subroutine initglobaldata
!=============================================================================
subroutine checkw(checkprimitive,ixImin1,ixImax1,ixOmin1,ixOmax1,w,flag)

include 'amrvacdef.f'
  
logical, intent(in)          :: checkprimitive
integer, intent(in)          :: ixImin1,ixImax1,ixOmin1,ixOmax1
double precision, intent(in) :: w(ixImin1:ixImax1,1:nw)
logical, intent(out)         :: flag(ixGlo1:ixGhi1)
!-----------------------------------------------------------------------------

flag(ixGlo1:ixGhi1)=.true.


if (checkprimitive) then
  if(useprimitiveRel)then
     ! check   rho>=0, p>=smallp
     flag(ixOmin1:ixOmax1) = (w(ixOmin1:ixOmax1,rho_) > minrho).and. &
                   (w(ixOmin1:ixOmax1,pp_)  >=minp &
)
  else
    ! check  v2 < 1, rho>=0, p>=minp
     ! v2 < 1, rho>0, p>0
     flag(ixOmin1:ixOmax1) = ((w(ixOmin1:ixOmax1,v1_)**2.0d0&
        +w(ixOmin1:ixOmax1,v2_)**2.0d0+w(ixOmin1:ixOmax1,v3_)&
        **2.0d0)< one).and. &
                   (w(ixOmin1:ixOmax1,rho_) > minrho).and. &
                   (w(ixOmin1:ixOmax1,pp_) >=minp&
)
  endif
else
  ! checks on the conservative variables
  flag(ixOmin1:ixOmax1)= (w(ixOmin1:ixOmax1,d_)   > minrho).and. &
               (w(ixOmin1:ixOmax1,tau_) > smalltau&
)
end if



end subroutine checkw
!=============================================================================
subroutine conserve(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x,patchw)

! Transform primitive variables into conservative ones
! (rho,v,p,B,E) ---> (D,S,tau,B,E,lfac,xi)
! call to smallvalues
! --> latter only used for correcting procedure in correctaux
! --> input array patchw for spatially selective transformation

include 'amrvacdef.f'

integer, intent(in)               :: ixImin1,ixImax1,ixOmin1,ixOmax1
double precision, intent(inout)   :: w(ixImin1:ixImax1,nw)
logical, intent(in)               :: patchw(ixGlo1:ixGhi1)
double precision, intent(in)      :: x(ixImin1:ixImax1,1:ndim)

!-----------------------------------------------------------------------------

call conserven(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x,patchw)

if(fixsmall) call smallvalues(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,&
   patchw(ixOmin1:ixOmax1),'conserve')

end subroutine conserve
!=============================================================================
subroutine conserven(ixImin1,ixImax1,ixOmin1,ixOmax1,w,patchw)

! Transform primitive variables into conservative ones
! (rho,v,p,B) ---> (D,S,tau,B,lfac,xi)
! no call to smallvalues
! --> latter only used for correcting procedure in correctaux
! --> input array patchw for spatially selective transformation

include 'amrvacdef.f'

integer, intent(in)               :: ixImin1,ixImax1,ixOmin1,ixOmax1
double precision, intent(inout)   :: w(ixImin1:ixImax1,nw)
logical, intent(in)               :: patchw(ixGlo1:ixGhi1)
! .. local ..
double precision, dimension(ixGlo1:ixGhi1):: ecrossbi, sqrU, sqrE, sqrB, rhoh
!-----------------------------------------------------------------------------

! fill the auxilary variable lfac (lorentz factor)
where(.not.patchw(ixOmin1:ixOmax1))
   sqrU(ixOmin1:ixOmax1)    = w(ixOmin1:ixOmax1,u1_)**2+w(ixOmin1:ixOmax1,&
      u2_)**2+w(ixOmin1:ixOmax1,u3_)**2
   w(ixOmin1:ixOmax1,lfac_) = dsqrt(one+sqrU(ixOmin1:ixOmax1)) 
endwhere

! fill the auxilary variable xi and density D
! with enthalpy w: xi= lfac^2 rhoh
! density: d = lfac * rho
call Enthalpy(w,ixImin1,ixImax1,ixOmin1,ixOmax1,patchw,rhoh)
where(.not.patchw(ixOmin1:ixOmax1))
   w(ixOmin1:ixOmax1,xi_)  = w(ixOmin1:ixOmax1,lfac_)*w(ixOmin1:ixOmax1,&
      lfac_)* rhoh(ixOmin1:ixOmax1)
   w(ixOmin1:ixOmax1,d_)   = w(ixOmin1:ixOmax1,rho_)*w(ixOmin1:ixOmax1,lfac_)
endwhere



! fill the vector S
! s = E x B + xi v
 
call ecrossb(ixImin1,ixImax1,ixOmin1,ixOmax1,1,w,patchw,ecrossbi)
where(.not.patchw(ixOmin1:ixOmax1))
   w(ixOmin1:ixOmax1,s1_) = ecrossbi(ixOmin1:ixOmax1) + w(ixOmin1:ixOmax1,&
      xi_) * w(ixOmin1:ixOmax1,u1_)/w(ixOmin1:ixOmax1,lfac_)
end where

  
call ecrossb(ixImin1,ixImax1,ixOmin1,ixOmax1,2,w,patchw,ecrossbi)
where(.not.patchw(ixOmin1:ixOmax1))
   w(ixOmin1:ixOmax1,s2_) = ecrossbi(ixOmin1:ixOmax1) + w(ixOmin1:ixOmax1,&
      xi_) * w(ixOmin1:ixOmax1,u2_)/w(ixOmin1:ixOmax1,lfac_)
end where

  
call ecrossb(ixImin1,ixImax1,ixOmin1,ixOmax1,3,w,patchw,ecrossbi)
where(.not.patchw(ixOmin1:ixOmax1))
   w(ixOmin1:ixOmax1,s3_) = ecrossbi(ixOmin1:ixOmax1) + w(ixOmin1:ixOmax1,&
      xi_) * w(ixOmin1:ixOmax1,u3_)/w(ixOmin1:ixOmax1,lfac_)
end where


! tau = 1/2 (E^2+B^2) + xi - p - d
! instead of E use tau= E - D
where(.not.patchw(ixOmin1:ixOmax1))
   sqrB(ixOmin1:ixOmax1)    = w(ixOmin1:ixOmax1,B1_)**2+w(ixOmin1:ixOmax1,&
      B2_)**2+w(ixOmin1:ixOmax1,B3_)**2
   sqrE(ixOmin1:ixOmax1)    = w(ixOmin1:ixOmax1,e1_)**2+w(ixOmin1:ixOmax1,&
      e2_)**2+w(ixOmin1:ixOmax1,e3_)**2
   w(ixOmin1:ixOmax1,tau_)  = half * (sqrE(ixOmin1:ixOmax1)&
      +sqrB(ixOmin1:ixOmax1)) + w(ixOmin1:ixOmax1,xi_) - w(ixOmin1:ixOmax1,&
      pp_) - w(ixOmin1:ixOmax1,d_)
endwhere

end subroutine conserven
!=============================================================================
subroutine ecrossb(ixImin1,ixImax1,ixOmin1,ixOmax1,idir,w,patchw,res)

include 'amrvacdef.f'

integer, intent(in)                :: ixImin1,ixImax1,ixOmin1,ixOmax1
double precision, intent(in)       :: w(ixImin1:ixImax1,1:nw)
integer         , intent(in)       :: idir
logical, intent(in)                :: patchw(ixGlo1:ixGhi1)
double precision, intent(out)      :: res(ixGlo1:ixGhi1)
! .. local ..
integer                            :: j,k
!-----------------------------------------------------------------------------

res(ixOmin1:ixOmax1) = zero
do j=1,3
   do k=1,3
      if (j .eq. k .or. j .eq. idir .or. k .eq. idir) cycle
      where(.not.patchw(ixOmin1:ixOmax1))
         res(ixOmin1:ixOmax1) = res(ixOmin1:ixOmax1) + lvc(idir,j,k)&
            *w(ixOmin1:ixOmax1,e0_+j)*w(ixOmin1:ixOmax1,b0_+k)
      endwhere
   end do
end do

end subroutine ecrossb
!=============================================================================
subroutine vcrossb(ixImin1,ixImax1,ixOmin1,ixOmax1,idir,w,x,patchw,res)

! Used with conserved variables, so call to getv needed.
! assuming auxilaries are set.  
include 'amrvacdef.f'

integer, intent(in)                :: ixImin1,ixImax1,ixOmin1,ixOmax1
double precision, intent(in)       :: w(ixImin1:ixImax1,1:nw)
double precision, intent(in)       :: x(ixImin1:ixImax1,1:ndim)
integer         , intent(in)       :: idir
logical, intent(in)                :: patchw(ixGlo1:ixGhi1)
double precision, intent(out)      :: res(ixGlo1:ixGhi1)
! .. local ..
integer                            :: j,k
double precision, dimension(ixGlo1:ixGhi1) :: vj
!-----------------------------------------------------------------------------

res(ixOmin1:ixOmax1) = zero
do j=1,3
   do k=1,3
      if (j .eq. k .or. j .eq. idir .or. k .eq. idir) cycle
      call getv(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,j,vj)
      where(.not.patchw(ixOmin1:ixOmax1))
         res(ixOmin1:ixOmax1) = res(ixOmin1:ixOmax1) + lvc(idir,j,k)&
            *vj(ixOmin1:ixOmax1)*w(ixOmin1:ixOmax1,b0_+k)
      endwhere
   end do
end do

end subroutine vcrossb
!=============================================================================
subroutine primitive(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x)

! Transform conservative variables into primitive ones
! (D,S,tau,B)-->(rho,v,p,B,lfac,xi)

include 'amrvacdef.f'

integer, intent(in)                :: ixImin1,ixImax1,ixOmin1,ixOmax1
double precision, intent(inout)    :: w(ixImin1:ixImax1,1:nw)
double precision, intent(in)       :: x(ixImin1:ixImax1,1:ndim)
! .. local ..
logical, dimension(ixGlo1:ixGhi1)          :: patchw
!-----------------------------------------------------------------------------

patchw(ixOmin1:ixOmax1) = .false.
! calculate lorentz factor and xi from conservatives only
call getaux(.true.,w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,'primitive')
call primitiven(ixImin1,ixImax1,ixOmin1,ixOmax1,w,patchw)



if (tlow>zero) call fixp_usr(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x)


end subroutine primitive
!==============================================================================
subroutine primitiven(ixImin1,ixImax1,ixOmin1,ixOmax1,w,patchw)

! same as primitive, but with padding by patchw
!  --> needed in correctaux, avoiding recursive calling by smallvalues
! assumes updated xi and lfac, and does not use smallxi cutoff
! Transform conservative variables into primitive ones
! (D,S,tau,B)-->(rho,v,p,B,lfac,xi)

include 'amrvacdef.f'

integer, intent(in)                    :: ixImin1,ixImax1,ixOmin1,ixOmax1
double precision, intent(inout)        :: w(ixImin1:ixImax1,1:nw)
logical, intent(in),dimension(ixGlo1:ixGhi1)   :: patchw
! .. local ..
double precision, dimension(ixGlo1:ixGhi1)     :: sqrB, tmpP, ecrossbi
!-----------------------------------------------------------------------------

call Pressuren(w,ixImin1,ixImax1,ixOmin1,ixOmax1,.true.,tmpP,patchw)
where(.not.patchw(ixOmin1:ixOmax1))
   w(ixOmin1:ixOmax1,pp_)=tmpP(ixOmin1:ixOmax1)
end where

where(.not.patchw(ixOmin1:ixOmax1))
   w(ixOmin1:ixOmax1,rho_)=w(ixOmin1:ixOmax1,d_)/w(ixOmin1:ixOmax1,lfac_)
!   sqrB(ixO^S) = {^C&w(ixO^S,b^C_)**2.0d0+}
end where


! u = lfac/xi * (S - E x B)

call ecrossb(ixImin1,ixImax1,ixOmin1,ixOmax1,1,w,patchw,ecrossbi)
where(.not.patchw(ixOmin1:ixOmax1))
   w(ixOmin1:ixOmax1,u1_) = w(ixOmin1:ixOmax1,lfac_)/w(ixOmin1:ixOmax1,xi_) &
      * (w(ixOmin1:ixOmax1,s1_)-ecrossbi(ixOmin1:ixOmax1))
end where


call ecrossb(ixImin1,ixImax1,ixOmin1,ixOmax1,2,w,patchw,ecrossbi)
where(.not.patchw(ixOmin1:ixOmax1))
   w(ixOmin1:ixOmax1,u2_) = w(ixOmin1:ixOmax1,lfac_)/w(ixOmin1:ixOmax1,xi_) &
      * (w(ixOmin1:ixOmax1,s2_)-ecrossbi(ixOmin1:ixOmax1))
end where


call ecrossb(ixImin1,ixImax1,ixOmin1,ixOmax1,3,w,patchw,ecrossbi)
where(.not.patchw(ixOmin1:ixOmax1))
   w(ixOmin1:ixOmax1,u3_) = w(ixOmin1:ixOmax1,lfac_)/w(ixOmin1:ixOmax1,xi_) &
      * (w(ixOmin1:ixOmax1,s3_)-ecrossbi(ixOmin1:ixOmax1))
end where




end subroutine primitiven
!==============================================================================
subroutine e_to_rhos(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x)

include 'amrvacdef.f'
integer:: ixImin1,ixImax1,ixOmin1,ixOmax1
double precision:: w(ixImin1:ixImax1,nw)
double precision, intent(in)      :: x(ixImin1:ixImax1,1:ndim)
!-------------------------------------------------------------------------

call mpistop("e to rhos for SRMHDEOS unavailable")

end subroutine e_to_rhos
!=============================================================================
subroutine rhos_to_e(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x)

include 'amrvacdef.f'

integer:: ixImin1,ixImax1,ixOmin1,ixOmax1
double precision:: w(ixImin1:ixImax1,nw)
double precision, intent(in)      :: x(ixImin1:ixImax1,1:ndim)
!-----------------------------------------------------------------------------

call mpistop("e to rhos for SRMHDEOS unavailable")

end subroutine rhos_to_e
!=============================================================================
subroutine ppmflatcd(ixImin1,ixImax1,ixOmin1,ixOmax1,ixLmin1,ixLmax1,ixRmin1,&
   ixRmax1,w,d2w,drho,dp)

include 'amrvacdef.f'

integer, intent(in)           :: ixImin1,ixImax1,ixOmin1,ixOmax1,ixLmin1,&
   ixLmax1,ixRmin1,ixRmax1
double precision, intent(in)  :: w(ixImin1:ixImax1,nw),d2w(ixGlo1:ixGhi1,&
   1:nwflux)

double precision, intent(inout) :: drho(ixGlo1:ixGhi1),dp(ixGlo1:ixGhi1)

!-----------------------------------------------------------------------------

call mpistop("PPM with flatsh=.true. can not be used with srrmhd !")


end subroutine ppmflatcd
!=============================================================================
subroutine ppmflatsh(ixImin1,ixImax1,ixOmin1,ixOmax1,ixLLmin1,ixLLmax1,&
   ixLmin1,ixLmax1,ixRmin1,ixRmax1,ixRRmin1,ixRRmax1,idims,w,drho,dp,dv)

include 'amrvacdef.f'

integer, intent(in)           :: ixImin1,ixImax1,ixOmin1,ixOmax1,ixLLmin1,&
   ixLLmax1,ixLmin1,ixLmax1,ixRmin1,ixRmax1,ixRRmin1,ixRRmax1
integer, intent(in)           :: idims
double precision, intent(in)  :: w(ixImin1:ixImax1,nw)

double precision, intent(inout) :: drho(ixGlo1:ixGhi1),dp(ixGlo1:ixGhi1),&
   dv(ixGlo1:ixGhi1)

!-----------------------------------------------------------------------------

call mpistop("PPM with flatsh=.true. can not be used with srrmhd !")

end subroutine ppmflatsh
!=============================================================================
subroutine getv(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,idims,v)

! v = (s - E x B)/xi

include 'amrvacdef.f'
  
integer, intent(in)                              :: ixImin1,ixImax1,ixOmin1,&
   ixOmax1,idims
double precision, intent(in)                     :: w(ixImin1:ixImax1,1:nw)
double precision, intent(in)                     :: x(ixImin1:ixImax1,1:ndim)
double precision, dimension(ixGlo1:ixGhi1), intent(out)  :: v
! .. local ..
double precision, dimension(ixGlo1:ixGhi1)               :: ecrossbi
logical, dimension(ixGlo1:ixGhi1)                        :: patchw
!-----------------------------------------------------------------------------

patchw(ixOmin1:ixOmax1) = .false.
call ecrossb(ixImin1,ixImax1,ixOmin1,ixOmax1,idims,w,patchw,ecrossbi)
where(w(ixOmin1:ixOmax1,xi_) <= smallxi)
   v(ixOmin1:ixOmax1) = zero
elsewhere
   v(ixOmin1:ixOmax1) = ( w(ixOmin1:ixOmax1,s0_+idims) - ecrossbi&
      (ixOmin1:ixOmax1) ) / w(ixOmin1:ixOmax1,xi_)
endwhere

end subroutine getv
!=============================================================================
subroutine getcmax(new_cmax,w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,idims,cmax,&
   cmin,needcmin)
  
! Calculate cmax_idim within ixO^L
! new_cmax is not used
  
include 'amrvacdef.f'
  
logical, intent(in)           :: new_cmax,needcmin
integer, intent(in)           :: ixImin1,ixImax1,ixOmin1,ixOmax1,idims
double precision, intent(in)  :: w(ixImin1:ixImax1,1:nw)
double precision, intent(in)  :: x(ixImin1:ixImax1,1:ndim)
double precision, intent(out) :: cmax(ixGlo1:ixGhi1),cmin(ixGlo1:ixGhi1)
!-----------------------------------------------------------------------------

cmax(ixOmin1:ixOmax1) =   one
cmin(ixOmin1:ixOmax1) = - one

end subroutine getcmax
!=============================================================================
subroutine getflux(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,iw,idims,f,transport)

! Calculate non-transport flux f_idim[iw] within ixO^L.
include 'amrvacdef.f'

integer, intent(in)                :: ixImin1,ixImax1,ixOmin1,ixOmax1,iw,idims
double precision, intent(in)       :: w(ixImin1:ixImax1,nw)
double precision, intent(in)       :: x(ixImin1:ixImax1,1:ndim)
double precision, intent(out)      :: f(ixGlo1:ixGhi1)
logical, intent(out)               :: transport
! .. local ..
integer                            :: k
double precision, dimension(ixGlo1:ixGhi1) :: vidims, vc, p, sqrE, sqrB,&
    ecrossbi
!-----------------------------------------------------------------------------
transport=.true.

if (iw==d_) then
   f(ixOmin1:ixOmax1)=zero


! f_i[phib_]=b_{i}
else if(iw==phib_) then
   transport=.false. 

       f(ixOmin1:ixOmax1)=w(ixOmin1:ixOmax1,b0_+idims)

! f_i[psi_]=e_{i}
else if(iw==psi_) then
   transport=.false. 
   f(ixOmin1:ixOmax1)=w(ixOmin1:ixOmax1,e0_+idims)

! F^i[q] = q vi (transport term) + lfac*sigma*(Ei+ v x B - (E dot v) vi)
else if (iw==q_) then
   transport = .true.
   call vcrossb(ixImin1,ixImax1,ixOmin1,ixOmax1,idims,w,x,patchfalse,ecrossbi) !storing vxB in tmpvar for ecrossb
   sqrE(ixOmin1:ixOmax1) =  zero
   do k = 1, 3
      call getv(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,k,vidims)
      sqrE(ixOmin1:ixOmax1) = sqrE(ixOmin1:ixOmax1) + w(ixOmin1:ixOmax1,e0_&
         +k)*vidims(ixOmin1:ixOmax1) !storing EdotV in tmpvar for sqrE
   end do
   call getv(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,idims,vidims)
   f(ixOmin1:ixOmax1) = w(ixOmin1:ixOmax1,lfac_)/eqpar(eta_) &
      * (w(ixOmin1:ixOmax1,e0_+idims) + ecrossbi(ixOmin1:ixOmax1) &
      - sqrE(ixOmin1:ixOmax1) * vidims(ixOmin1:ixOmax1))
else

!f^i[B_c] = lvc(c,i,k) Ek + phi delta(i,c)
if (iw==b1_) then
   transport=.false. 
   if (idims .ne. 1) then 
      do k=1,3
         if (idims .eq. k .or. 1 .eq. k) cycle
         f(ixOmin1:ixOmax1) = lvc(1,idims,k) * w(ixOmin1:ixOmax1,e0_+k)
      end do
   else
      f(ixOmin1:ixOmax1)=w(ixOmin1:ixOmax1,phib_)
   end if
end if
if (iw==b2_) then
   transport=.false. 
   if (idims .ne. 2) then 
      do k=1,3
         if (idims .eq. k .or. 2 .eq. k) cycle
         f(ixOmin1:ixOmax1) = lvc(2,idims,k) * w(ixOmin1:ixOmax1,e0_+k)
      end do
   else
      f(ixOmin1:ixOmax1)=w(ixOmin1:ixOmax1,phib_)
   end if
end if
if (iw==b3_) then
   transport=.false. 
   if (idims .ne. 3) then 
      do k=1,3
         if (idims .eq. k .or. 3 .eq. k) cycle
         f(ixOmin1:ixOmax1) = lvc(3,idims,k) * w(ixOmin1:ixOmax1,e0_+k)
      end do
   else
      f(ixOmin1:ixOmax1)=w(ixOmin1:ixOmax1,phib_)
   end if
end if

!f^i[E_c] = - lvc(c,i,k) Bk + psi delta(i,c)
if (iw==e1_) then
   transport=.false. 
   if (idims .ne. 1) then 
      do k=1,3
         if (idims .eq. k .or. 1 .eq. k) cycle
         f(ixOmin1:ixOmax1) = - lvc(1,idims,k) * w(ixOmin1:ixOmax1,b0_+k)
      end do
   else
      f(ixOmin1:ixOmax1) = w(ixOmin1:ixOmax1,psi_)
   end if
end if
if (iw==e2_) then
   transport=.false. 
   if (idims .ne. 2) then 
      do k=1,3
         if (idims .eq. k .or. 2 .eq. k) cycle
         f(ixOmin1:ixOmax1) = - lvc(2,idims,k) * w(ixOmin1:ixOmax1,b0_+k)
      end do
   else
      f(ixOmin1:ixOmax1) = w(ixOmin1:ixOmax1,psi_)
   end if
end if
if (iw==e3_) then
   transport=.false. 
   if (idims .ne. 3) then 
      do k=1,3
         if (idims .eq. k .or. 3 .eq. k) cycle
         f(ixOmin1:ixOmax1) = - lvc(3,idims,k) * w(ixOmin1:ixOmax1,b0_+k)
      end do
   else
      f(ixOmin1:ixOmax1) = w(ixOmin1:ixOmax1,psi_)
   end if
end if


sqrE(ixOmin1:ixOmax1) =  w(ixOmin1:ixOmax1,e0_+1)**2 + w(ixOmin1:ixOmax1,e0_&
   +2)**2 + w(ixOmin1:ixOmax1,e0_+3)**2 
sqrB(ixOmin1:ixOmax1) =  w(ixOmin1:ixOmax1,b0_+1)**2 + w(ixOmin1:ixOmax1,b0_&
   +2)**2 + w(ixOmin1:ixOmax1,b0_+3)**2 

! f^i[S_c] = - E_c E_i - B_c B_i + xi v_c v_i + (1/2(E**2+B**2)+p) delta(c,i)
if (iw==s1_) then
   transport=.false. 
   call getv(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,idims,vidims)
   call getv(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,1,vc)
   f(ixOmin1:ixOmax1) = - w(ixOmin1:ixOmax1,e1_)*w(ixOmin1:ixOmax1,e0_&
      +idims) - w(ixOmin1:ixOmax1,b1_)*w(ixOmin1:ixOmax1,b0_+idims) &
        + w(ixOmin1:ixOmax1,xi_)*vidims(ixOmin1:ixOmax1)*vc(ixOmin1:ixOmax1)
   if (idims==1) then !!! add (1/2(E**2+B**2)+p)
      call Pressure(w,ixImin1,ixImax1,ixOmin1,ixOmax1,.true.,p)
      f(ixOmin1:ixOmax1)    = f(ixOmin1:ixOmax1) + half*(sqrE&
         (ixOmin1:ixOmax1) + sqrB(ixOmin1:ixOmax1)) + p(ixOmin1:ixOmax1)
   end if
end if
if (iw==s2_) then
   transport=.false. 
   call getv(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,idims,vidims)
   call getv(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,2,vc)
   f(ixOmin1:ixOmax1) = - w(ixOmin1:ixOmax1,e2_)*w(ixOmin1:ixOmax1,e0_&
      +idims) - w(ixOmin1:ixOmax1,b2_)*w(ixOmin1:ixOmax1,b0_+idims) &
        + w(ixOmin1:ixOmax1,xi_)*vidims(ixOmin1:ixOmax1)*vc(ixOmin1:ixOmax1)
   if (idims==2) then !!! add (1/2(E**2+B**2)+p)
      call Pressure(w,ixImin1,ixImax1,ixOmin1,ixOmax1,.true.,p)
      f(ixOmin1:ixOmax1)    = f(ixOmin1:ixOmax1) + half*(sqrE&
         (ixOmin1:ixOmax1) + sqrB(ixOmin1:ixOmax1)) + p(ixOmin1:ixOmax1)
   end if
end if
if (iw==s3_) then
   transport=.false. 
   call getv(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,idims,vidims)
   call getv(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,3,vc)
   f(ixOmin1:ixOmax1) = - w(ixOmin1:ixOmax1,e3_)*w(ixOmin1:ixOmax1,e0_&
      +idims) - w(ixOmin1:ixOmax1,b3_)*w(ixOmin1:ixOmax1,b0_+idims) &
        + w(ixOmin1:ixOmax1,xi_)*vidims(ixOmin1:ixOmax1)*vc(ixOmin1:ixOmax1)
   if (idims==3) then !!! add (1/2(E**2+B**2)+p)
      call Pressure(w,ixImin1,ixImax1,ixOmin1,ixOmax1,.true.,p)
      f(ixOmin1:ixOmax1)    = f(ixOmin1:ixOmax1) + half*(sqrE&
         (ixOmin1:ixOmax1) + sqrB(ixOmin1:ixOmax1)) + p(ixOmin1:ixOmax1)
   end if
end if

! f^i[tau] = tau vi (transport term) + vi*(p -1/2 E^2 -1/2 B^2) + E x B
if (iw==e_) then
   transport = .true.
   call Pressure(w,ixImin1,ixImax1,ixOmin1,ixOmax1,.true.,p)
   call getv(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,idims,vidims)
   call ecrossb(ixImin1,ixImax1,ixOmin1,ixOmax1,idims,w,patchfalse,ecrossbi)
   f(ixOmin1:ixOmax1) = vidims(ixOmin1:ixOmax1) * (p(ixOmin1:ixOmax1) &
      - half*(sqrE(ixOmin1:ixOmax1)+sqrB(ixOmin1:ixOmax1)) ) &
      + ecrossbi(ixOmin1:ixOmax1)
end if

end if

end subroutine getflux
!=============================================================================
subroutine getfluxforhllc(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,iw,idims,f,&
   transport)

! Calculate non-transport flux f_idim[iw] within ixO^L.

include 'amrvacdef.f'

integer, intent(in)::          ixImin1,ixImax1,ixOmin1,ixOmax1,iw,idims
double precision, intent(in)  :: w(ixImin1:ixImax1,nw)
double precision, intent(in)      :: x(ixImin1:ixImax1,1:ndim)
double precision, intent(out) :: f(ixGlo1:ixGhi1,1:nwflux)
logical, intent(out)          :: transport

!-----------------------------------------------------------------------------

call getflux(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,iw,idims,f(ixGlo1:ixGhi1,iw),&
   transport)

end subroutine getfluxforhllc
!=============================================================================
subroutine addgeometry(qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,wCT,w,x)

! Add geometrical source terms to w

include 'amrvacdef.f'

integer, intent(in)                :: ixImin1,ixImax1, ixOmin1,ixOmax1
double precision, intent(in)       :: qdt
double precision, intent(in)       :: x(ixImin1:ixImax1,1:ndim)
double precision, intent(inout)    :: wCT(ixImin1:ixImax1,1:nw),&
    w(ixImin1:ixImax1,1:nw)

integer                            :: iw,ix,idir, h1xmin1,h1xmax1
double precision, dimension(ixGlo1:ixGhi1) :: tmp, sqrVdotB, sqrB, P
logical                            :: angmomfix=.false.
!-----------------------------------------------------------------------------

if(typeaxial /= 'slab')then
  call mpistop('Only slab geometry implemented so far')
endif

end subroutine addgeometry
!=============================================================================
subroutine addsource(qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,iwmin,iwmax,qtC,wCT,&
   qt,w,x,qsourcesplit)

! w[iws]=w[iws]+qdt*S[iws,wCT] where S is the source based on wCT within ixO

include 'amrvacdef.f'

integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1, iwmin,&
   iwmax
double precision, intent(in)    :: qdt, qtC, qt
double precision, intent(in)    :: x(ixImin1:ixImax1,1:ndim)
double precision, intent(in)    :: wCT(ixImin1:ixImax1,1:nw)
double precision, intent(inout) :: w(ixImin1:ixImax1,1:nw)
logical, intent(in)             :: qsourcesplit

double precision :: dx1
!-----------------------------------------------------------------------------

dx1=dxlevel(1);

! Two sources: Sb is added via Strang-splitting, Sa is unsplit.
! Using ssplitdivb to distinguish (needs to be set true!!!)
if(qsourcesplit .eqv. ssplitdivb) then
   call addsource_b(qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,iwmin,iwmax,qtC,wCT,&
      qt,w,x,dx1)
else
   call addsource_a(qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,iwmin,iwmax,qtC,wCT,&
      qt,w,x,dx1)
end if

! now update the auxiliaries to the new state
call getaux(.true.,w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,'addsource')

end subroutine addsource
!=============================================================================
subroutine addsource_a(qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,iwmin,iwmax,qtC,&
   wCT,qt,w,x,dx1)
include 'amrvacdef.f'

integer, intent(in) :: ixImin1,ixImax1, ixOmin1,ixOmax1, iwmin,iwmax
double precision, intent(in) :: qdt, qtC, qt, wCT(ixImin1:ixImax1,1:nw),&
    x(ixImin1:ixImax1,1:ndim)
double precision, intent(in) :: dx1
double precision, intent(inout) :: w(ixImin1:ixImax1,1:nw)
! .. local ..
double precision, dimension(ixGlo1:ixGhi1)   :: vidir
!-----------------------------------------------------------------------------

! S[psi_] = q
w(ixOmin1:ixOmax1,psi_) = w(ixOmin1:ixOmax1,psi_) + qdt * wCT(ixOmin1:ixOmax1,&
   q_)

! S[Ei_] = - q vi

call getv(wCT,x,ixImin1,ixImax1,ixOmin1,ixOmax1,1,vidir)
w(ixOmin1:ixOmax1,e1_) = w(ixOmin1:ixOmax1,e1_) - qdt * wCT(ixOmin1:ixOmax1,&
   q_) * vidir(ixOmin1:ixOmax1)


call getv(wCT,x,ixImin1,ixImax1,ixOmin1,ixOmax1,2,vidir)
w(ixOmin1:ixOmax1,e2_) = w(ixOmin1:ixOmax1,e2_) - qdt * wCT(ixOmin1:ixOmax1,&
   q_) * vidir(ixOmin1:ixOmax1)


call getv(wCT,x,ixImin1,ixImax1,ixOmin1,ixOmax1,3,vidir)
w(ixOmin1:ixOmax1,e3_) = w(ixOmin1:ixOmax1,e3_) - qdt * wCT(ixOmin1:ixOmax1,&
   q_) * vidir(ixOmin1:ixOmax1)


end subroutine addsource_a
!=============================================================================
subroutine addsource_b(qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,iwmin,iwmax,qtC,&
   wCT,qt,w,x,dx1)
include 'amrvacdef.f'

integer, intent(in) :: ixImin1,ixImax1, ixOmin1,ixOmax1, iwmin,iwmax
double precision, intent(in) :: qdt, qtC, qt, wCT(ixImin1:ixImax1,1:nw),&
    x(ixImin1:ixImax1,1:ndim)
double precision, intent(in) :: dx1
double precision, intent(inout) :: w(ixImin1:ixImax1,1:nw)
! .. local ..
double precision, dimension(ixGlo1:ixGhi1,1:ndir)  :: Epar, Eperp, Eperpstar,&
    v
double precision, dimension(ixGlo1:ixGhi1)         :: vabs
!-----------------------------------------------------------------------------

! S[phib_] = -kappa phib
w(ixOmin1:ixOmax1,phib_) = w(ixOmin1:ixOmax1,phib_)*exp(-qdt*eqpar(kappa_))

! S[psi_] = - kappa psi
w(ixOmin1:ixOmax1,psi_) = w(ixOmin1:ixOmax1,psi_)*exp(-qdt*eqpar(kappa_))

! S[E] = - sigma lfac (E + v x B - (Edotv) v)
! Split the source in parallel and perpendicular components
! Then solve the linear evolution equation (assuming constant background state)
! Then add again both components.

 
call getv(wCT,x,ixImin1,ixImax1,ixOmin1,ixOmax1,1,v(ixGlo1:ixGhi1,1))

  
call getv(wCT,x,ixImin1,ixImax1,ixOmin1,ixOmax1,2,v(ixGlo1:ixGhi1,2))

  
call getv(wCT,x,ixImin1,ixImax1,ixOmin1,ixOmax1,3,v(ixGlo1:ixGhi1,3))

vabs(ixOmin1:ixOmax1) = sqrt( v(ixOmin1:ixOmax1,1)**2+ v(ixOmin1:ixOmax1,2)&
   **2+ v(ixOmin1:ixOmax1,3)**2)


Epar(ixOmin1:ixOmax1,1)  = w(ixOmin1:ixOmax1,e1_)*v(ixOmin1:ixOmax1,1)&
   /vabs(ixOmin1:ixOmax1)
Eperp(ixOmin1:ixOmax1,1) = w(ixOmin1:ixOmax1,e1_) - Epar(ixOmin1:ixOmax1,1)
call vcrossb(ixImin1,ixImax1,ixOmin1,ixOmax1,1,wCT,x,patchfalse,&
   Eperpstar(ixGlo1:ixGhi1,1))
Eperpstar(ixOmin1:ixOmax1,1) = - Eperpstar(ixOmin1:ixOmax1,1)


Epar(ixOmin1:ixOmax1,2)  = w(ixOmin1:ixOmax1,e2_)*v(ixOmin1:ixOmax1,2)&
   /vabs(ixOmin1:ixOmax1)
Eperp(ixOmin1:ixOmax1,2) = w(ixOmin1:ixOmax1,e2_) - Epar(ixOmin1:ixOmax1,2)
call vcrossb(ixImin1,ixImax1,ixOmin1,ixOmax1,2,wCT,x,patchfalse,&
   Eperpstar(ixGlo1:ixGhi1,2))
Eperpstar(ixOmin1:ixOmax1,2) = - Eperpstar(ixOmin1:ixOmax1,2)


Epar(ixOmin1:ixOmax1,3)  = w(ixOmin1:ixOmax1,e3_)*v(ixOmin1:ixOmax1,3)&
   /vabs(ixOmin1:ixOmax1)
Eperp(ixOmin1:ixOmax1,3) = w(ixOmin1:ixOmax1,e3_) - Epar(ixOmin1:ixOmax1,3)
call vcrossb(ixImin1,ixImax1,ixOmin1,ixOmax1,3,wCT,x,patchfalse,&
   Eperpstar(ixGlo1:ixGhi1,3))
Eperpstar(ixOmin1:ixOmax1,3) = - Eperpstar(ixOmin1:ixOmax1,3)



! Handle zero velocity separate:
where (vabs(ixOmin1:ixOmax1) .lt. smalldouble)

     w(ixOmin1:ixOmax1,e1_) = w(ixOmin1:ixOmax1,e1_) * exp(-qdt/eqpar(eta_))


     w(ixOmin1:ixOmax1,e2_) = w(ixOmin1:ixOmax1,e2_) * exp(-qdt/eqpar(eta_))


     w(ixOmin1:ixOmax1,e3_) = w(ixOmin1:ixOmax1,e3_) * exp(-qdt/eqpar(eta_))

elsewhere

     Epar(ixOmin1:ixOmax1,1)  = Epar(ixOmin1:ixOmax1,1) * exp(&
        -qdt/(eqpar(eta_)*wCT(ixOmin1:ixOmax1,lfac_)))
     Eperp(ixOmin1:ixOmax1,1) = Eperpstar(ixOmin1:ixOmax1,1) &
        + (Eperp(ixOmin1:ixOmax1,1) - Eperpstar(ixOmin1:ixOmax1,1))&
          * exp(-qdt*wCT(ixOmin1:ixOmax1,lfac_)/eqpar(eta_))
     w(ixOmin1:ixOmax1,e1_)   = Epar(ixOmin1:ixOmax1,1) + Eperp&
        (ixOmin1:ixOmax1,1)


     Epar(ixOmin1:ixOmax1,2)  = Epar(ixOmin1:ixOmax1,2) * exp(&
        -qdt/(eqpar(eta_)*wCT(ixOmin1:ixOmax1,lfac_)))
     Eperp(ixOmin1:ixOmax1,2) = Eperpstar(ixOmin1:ixOmax1,2) &
        + (Eperp(ixOmin1:ixOmax1,2) - Eperpstar(ixOmin1:ixOmax1,2))&
          * exp(-qdt*wCT(ixOmin1:ixOmax1,lfac_)/eqpar(eta_))
     w(ixOmin1:ixOmax1,e2_)   = Epar(ixOmin1:ixOmax1,2) + Eperp&
        (ixOmin1:ixOmax1,2)


     Epar(ixOmin1:ixOmax1,3)  = Epar(ixOmin1:ixOmax1,3) * exp(&
        -qdt/(eqpar(eta_)*wCT(ixOmin1:ixOmax1,lfac_)))
     Eperp(ixOmin1:ixOmax1,3) = Eperpstar(ixOmin1:ixOmax1,3) &
        + (Eperp(ixOmin1:ixOmax1,3) - Eperpstar(ixOmin1:ixOmax1,3))&
          * exp(-qdt*wCT(ixOmin1:ixOmax1,lfac_)/eqpar(eta_))
     w(ixOmin1:ixOmax1,e3_)   = Epar(ixOmin1:ixOmax1,3) + Eperp&
        (ixOmin1:ixOmax1,3)

endwhere


end subroutine addsource_b
!=============================================================================
subroutine getcurrent(w,ixImin1,ixImax1,ixmin1,ixmax1,idirmin,current)

! Calculate idirmin and the idirmin:3 components of the common current array

include 'amrvacdef.f'

integer, parameter:: idirmin0=7-2*ndir
integer :: ixmin1,ixmax1, idirmin, ixImin1,ixImax1
double precision :: w(ixImin1:ixImax1,1:nw)

! For ndir=2 only 3rd component of J can exist, ndir=1 is impossible for MHD
double precision :: current(ixGlo1:ixGhi1,7-2*ndir:3),bvec(ixGlo1:ixGhi1,&
   1:ndir)
!-----------------------------------------------------------------------------

bvec(ixImin1:ixImax1,1:ndir)=w(ixImin1:ixImax1,b0_+1:b0_+ndir)
call curlvector(bvec,ixImin1,ixImax1,ixmin1,ixmax1,current,idirmin,idirmin0,&
   ndir)

end subroutine getcurrent
!=============================================================================
!=============================================================================
! end module amrvacphys/- srrmhd
!#############################################################################
