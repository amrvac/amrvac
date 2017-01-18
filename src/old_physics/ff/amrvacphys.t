!#############################################################################
! Module amrvacphys/- ff
! Maxwells equations plus Ohm's law according to Komissarov, Barkov and Lyutikov 2007
! Also evolving charge density and GLM for divE
! 2015-08-12 by Oliver Porth

! TODO: 
! Check treatment for B<smalldouble
! Adopt perpendicular conductivity depending on E2/B2 to drive E<B.

INCLUDE:amrvacnul/roe.t
INCLUDE:amrvacnul/hllc.t
INCLUDE:amrvacnul/getaux.t
!=============================================================================
subroutine checkglobaldata
use mod_global_parameters
!-----------------------------------------------------------------------------

if (ssplitresis .neqv. .true.) then
     call mpistop('Please run with ssplitresis = .true.')
end if 

! We require three vector components
if (^NC .ne. 3) &
     call mpistop('FF needs three vector components.')

end subroutine checkglobaldata
!=============================================================================
subroutine initglobaldata

! place to set entropy fixes etc, absent for now

use mod_global_parameters
!-----------------------------------------------------------------------------
eqpar(kpar_)   = 100.0d0
eqpar(kperp_)  = zero
eqpar(kappa_)  = eqpar(kpar_)/2.0d0

if (ssplitresis .neqv. .true.) then
   ssplitresis = .true.
   if (mype .eq. 0) write(*,*) 'Overwriting ssplitresis = .true. as required'
end if 


end subroutine initglobaldata
!=============================================================================
subroutine conserve(ixI^L,ixO^L,w,x,patchw)

use mod_global_parameters

integer, intent(in)               :: ixI^L,ixO^L
double precision, intent(inout)   :: w(ixI^S,nw)
logical, intent(in)               :: patchw(ixG^T)
double precision, intent(in)      :: x(ixI^S,1:ndim)
!-----------------------------------------------------------------------------


end subroutine conserve
!=============================================================================
subroutine ecrossb(ixI^L,ixO^L,idir,w,patchw,res)

use mod_global_parameters

integer, intent(in)                :: ixI^L,ixO^L
double precision, intent(in)       :: w(ixI^S,1:nw)
integer         , intent(in)       :: idir
logical, intent(in)                :: patchw(ixG^T)
double precision, intent(out)      :: res(ixG^T)
! .. local ..
integer                            :: j,k
!-----------------------------------------------------------------------------

res(ixO^S) = zero
do j=1,^NC
   do k=1,^NC
      if (j .eq. k .or. j .eq. idir .or. k .eq. idir) cycle
      where(.not.patchw(ixO^S))
         res(ixO^S) = res(ixO^S) + lvc(idir,j,k)*w(ixO^S,e0_+j)*w(ixO^S,b0_+k)
      endwhere
   end do
end do

end subroutine ecrossb
!=============================================================================
subroutine primitive(ixI^L,ixO^L,w,x)

use mod_global_parameters

integer, intent(in)                :: ixI^L,ixO^L
double precision, intent(inout)    :: w(ixI^S,1:nw)
double precision, intent(in)       :: x(ixI^S,1:ndim)
!-----------------------------------------------------------------------------

end subroutine primitive
!==============================================================================
subroutine e_to_rhos(ixI^L,ixO^L,w,x)

use mod_global_parameters
integer:: ixI^L,ixO^L
double precision:: w(ixI^S,nw)
double precision, intent(in)      :: x(ixI^S,1:ndim)
!-------------------------------------------------------------------------

call mpistop("e to rhos for FF unavailable")

end subroutine e_to_rhos
!=============================================================================
subroutine rhos_to_e(ixI^L,ixO^L,w,x)

use mod_global_parameters

integer:: ixI^L,ixO^L
double precision:: w(ixI^S,nw)
double precision, intent(in)      :: x(ixI^S,1:ndim)
!-----------------------------------------------------------------------------

call mpistop("e to rhos for FF unavailable")

end subroutine rhos_to_e
!=============================================================================
subroutine getv(w,x,ixI^L,ixO^L,idims,v)

! Drift velocity
! v = E x B /B^2 

use mod_global_parameters
  
integer, intent(in)                              :: ixI^L,ixO^L,idims
double precision, intent(in)                     :: w(ixI^S,1:nw)
double precision, intent(in)                     :: x(ixI^S,1:ndim)
double precision, dimension(ixG^T), intent(out)  :: v
! .. local ..
double precision, dimension(ixG^T)               :: ecrossbi, b2
logical, dimension(ixG^T)                        :: patchw
!-----------------------------------------------------------------------------

patchw(ixO^S) = .false.
call ecrossb(ixI^L,ixO^L,idims,w,patchw,ecrossbi)
b2(ixO^S) = ({^C& w(ixO^S,b^C_)**2|+})
where( b2(ixO^S) <= smalldouble)
   v(ixO^S) = zero
elsewhere
   v(ixO^S) = ecrossbi(ixO^S) / b2(ixO^S)
endwhere

end subroutine getv
!=============================================================================
subroutine getcmax(new_cmax,w,x,ixI^L,ixO^L,idims,cmax,cmin,needcmin)
  
! Calculate cmax_idim within ixO^L
! new_cmax is not used
  
use mod_global_parameters
  
logical, intent(in)           :: new_cmax,needcmin
integer, intent(in)           :: ixI^L,ixO^L,idims
double precision, intent(in)  :: w(ixI^S,1:nw)
double precision, intent(in)  :: x(ixI^S,1:ndim)
double precision, intent(out) :: cmax(ixG^T),cmin(ixG^T)
!-----------------------------------------------------------------------------

cmax(ixO^S) =   one
cmin(ixO^S) = - one

end subroutine getcmax
!=============================================================================
subroutine getflux(w,x,ixI^L,ixO^L,iw,idims,f,transport)

! Calculate non-transport flux f_idim[iw] within ixO^L.
use mod_global_parameters

integer, intent(in)                :: ixI^L,ixO^L,iw,idims
double precision, intent(in)       :: w(ixI^S,nw)
double precision, intent(in)       :: x(ixI^S,1:ndim)
double precision, intent(out)      :: f(ixG^T)
logical, intent(out)               :: transport
! .. local ..
integer                            :: k
double precision, dimension(ixG^T) :: b2, e2, Epari, Eperpi, EdotB
!-----------------------------------------------------------------------------
transport=.false.

if(iw==phib_) then
{#IFNDEF FCT
       f(ixO^S)=w(ixO^S,b0_+idims)
}{#IFDEF FCT
       f(ixO^S)=zero
}
else if(iw==psi_) then
   f(ixO^S)=w(ixO^S,e0_+idims)       

!f^i[q] = q E x B /B^2 (transport term) + kappa_par Epar_i + kappa_perp Eperp_i
else if(iw==q_) then
   transport=.true.
   b2(ixO^S)   = ({^C& w(ixO^S,b^C_)**2|+})
   e2(ixO^S) = ({^C& w(ixO^S,e^C_)**2|+})
   where (b2(ixO^S) .lt. smalldouble .or. e2(ixO^S) .gt. b2(ixO^S))
      f(ixO^S) = eqpar(kpar_) * w(ixO^S,e0_+idims)
   elsewhere
      EdotB(ixO^S) = {^C& w(ixO^S,e^C_)*w(ixO^S,b^C_) |+}
      Epari(ixO^S)  = EdotB(ixO^S)*w(ixO^S,b0_+idims)/b2(ixO^S)
      Eperpi(ixO^S) = w(ixO^S,e0_+idims) - Epari(ixO^S)
      f(ixO^S) = eqpar(kpar_)*Epari(ixO^S) + eqpar(kperp_)*Eperpi(ixO^S)
   end where
!f^i[B_c] = lvc(c,i,k) Ek + phi delta(i,c)
{else if (iw==b^C_) then
   if (idims .ne. ^C) then
      do k=1,^NC
         if (idims .eq. k .or. ^C .eq. k) cycle
         f(ixO^S) = lvc(^C,idims,k) * w(ixO^S,e0_+k)
      end do
   else
      f(ixO^S)=w(ixO^S,phib_)
   end if
\}
!f^i[E_c] = - lvc(c,i,k) Bk + psi delta(i,c)
{else if (iw==e^C_) then
   if (idims .ne. ^C) then 
      do k=1,^NC
         if (idims .eq. k .or. ^C .eq. k) cycle
         f(ixO^S) = - lvc(^C,idims,k) * w(ixO^S,b0_+k)
      end do
   else
      f(ixO^S)=w(ixO^S,psi_)
   end if
\}
end if

end subroutine getflux
!=============================================================================
subroutine getfluxforhllc(w,x,ixI^L,ixO^L,iw,idims,f,transport)

! Calculate non-transport flux f_idim[iw] within ixO^L.

use mod_global_parameters

integer, intent(in)::          ixI^L,ixO^L,iw,idims
double precision, intent(in)  :: w(ixI^S,nw)
double precision, intent(in)      :: x(ixI^S,1:ndim)
double precision, intent(out) :: f(ixG^T,1:nwflux)
logical, intent(out)          :: transport

!-----------------------------------------------------------------------------

call getflux(w,x,ixI^L,ixO^L,iw,idims,f(ixG^T,iw),transport)

end subroutine getfluxforhllc
!=============================================================================
subroutine addgeometry(qdt,ixI^L,ixO^L,wCT,w,x)

! Add geometrical source terms to w

use mod_global_parameters

integer, intent(in)                :: ixI^L, ixO^L
double precision, intent(in)       :: qdt
double precision, intent(in)       :: x(ixI^S,1:ndim)
double precision, intent(inout)    :: wCT(ixI^S,1:nw), w(ixI^S,1:nw)

integer                            :: iw,ix,idir, h1x^L{^NOONED, h2x^L}
double precision, dimension(ixG^T) :: tmp, sqrVdotB, sqrB, P
logical                            :: angmomfix=.false.
!-----------------------------------------------------------------------------

if(typeaxial /= 'slab')then
  call mpistop('Only slab geometry implemented so far')
endif

end subroutine addgeometry
!=============================================================================
subroutine addsource(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,qsourcesplit)

! w[iws]=w[iws]+qdt*S[iws,wCT] where S is the source based on wCT within ixO

use mod_global_parameters

integer, intent(in)             :: ixI^L, ixO^L, iw^LIM
double precision, intent(in)    :: qdt, qtC, qt
double precision, intent(in)    :: x(ixI^S,1:ndim)
double precision, intent(in)    :: wCT(ixI^S,1:nw)
double precision, intent(inout) :: w(ixI^S,1:nw)
logical, intent(in)             :: qsourcesplit
! .. local ..
double precision                :: dx^D
!-----------------------------------------------------------------------------
dx^D=dxlevel(^D);

! Two sources: Sb is added via Strang-splitting, Sa is unsplit.
if(qsourcesplit) then
   call addsource_b(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,dx^D)
else
   call addsource_a(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,dx^D)
end if

end subroutine addsource
!=============================================================================
subroutine addsource_a(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,dx^D)
use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L, iw^LIM
double precision, intent(in) :: qdt, qtC, qt, wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
double precision, intent(in) :: dx^D
double precision, intent(inout) :: w(ixI^S,1:nw)
! .. local ..
double precision, dimension(ixG^T)   :: vidir
!-----------------------------------------------------------------------------
! S[psi_] = q
w(ixO^S,psi_) = w(ixO^S,psi_) + qdt * wCT(ixO^S,q_)

! S[Ei_] = - q (E x B)/B^2
{^C&
call getv(wCT,x,ixI^L,ixO^L,^C,vidir)
w(ixO^S,e^C_) = w(ixO^S,e^C_) - qdt * wCT(ixO^S,q_) * vidir(ixO^S)
\}

end subroutine addsource_a
!=============================================================================
subroutine addsource_b(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,dx^D)
use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L, iw^LIM
double precision, intent(in) :: qdt, qtC, qt, wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
double precision, intent(in) :: dx^D
double precision, intent(inout) :: w(ixI^S,1:nw)
! .. local ..
double precision, dimension(ixG^T,1:ndir)  :: Epar, Eperp
double precision, dimension(ixG^T)         :: b2, e2, EdotB
double precision, parameter                :: dampfac = 100.0d0
!-----------------------------------------------------------------------------

! S[phib_] = -kappa phib
w(ixO^S,phib_) = w(ixO^S,phib_)*exp(-qdt*eqpar(kappa_))

! S[psi_] = -kappa psi
w(ixO^S,psi_) = w(ixO^S,psi_)*exp(-qdt*eqpar(kappa_))

! S[E] = - kappa_par Epar - kappa_perp Eperp
! Split the source in parallel and perpendicular components
! Then solve the linear evolution equation (assuming constant background state)
! Then add again both components.

b2(ixO^S) = ({^C& w(ixO^S,b^C_)**2|+})
e2(ixO^S) = ({^C& w(ixO^S,e^C_)**2|+})

EdotB(ixO^S) = {^C& w(ixO^S,e^C_)*w(ixO^S,b^C_) |+}
{^C&
where (b2(ixO^S) .ge. smalldouble) 
    Epar(ixO^S,^C)  = EdotB(ixO^S)*w(ixO^S,b^C_)/b2(ixO^S)
 elsewhere
    Epar(ixO^S,^C)  = zero
endwhere
Eperp(ixO^S,^C) = w(ixO^S,e^C_) - Epar(ixO^S,^C)
\}



! Handle zero magnetic field and E>B separate:
! Use increased parallel conductivity for all directions
! (kpar assumed larger than kperp)
where (b2(ixO^S) .lt. smalldouble)
{^C&
     w(ixO^S,e^C_) = w(ixO^S,e^C_) * exp(-dampfac*eqpar(kpar_)*qdt)
\}
elsewhere (e2(ixO^S) .gt. b2(ixO^S))
{^C&
     Epar(ixO^S,^C)  = Epar(ixO^S,^C) * exp(-dampfac*max(eqpar(kpar_),eqpar(kperp_))*qdt)
     Eperp(ixO^S,^C) = Eperp(ixO^S,^C) * exp(-dampfac*max(eqpar(kpar_),eqpar(kperp_))*qdt)
     w(ixO^S,e^C_)   = Epar(ixO^S,^C) + Eperp(ixO^S,^C)
\}
elsewhere
{^C&
     Epar(ixO^S,^C)  = Epar(ixO^S,^C) * exp(-eqpar(kpar_)*qdt)
     Eperp(ixO^S,^C) = Eperp(ixO^S,^C) * exp(-eqpar(kperp_)*qdt)
     w(ixO^S,e^C_)   = Epar(ixO^S,^C) + Eperp(ixO^S,^C)
\}
endwhere


end subroutine addsource_b
!=============================================================================
subroutine getdt(w,ixG^L,ix^L,dtnew,dx^D,x)

! Limit the timestep due to the resistive source addition
  
use mod_global_parameters

integer, intent(in)              :: ixG^L, ix^L
double precision, intent(in)     :: dx^D, x(ixG^S,1:ndim)
double precision, intent(inout)  :: w(ixG^S,1:nw), dtnew
!-----------------------------------------------------------------------------

dtnew = dtdiffpar / max(eqpar(kpar_),eqpar(kperp_))

end subroutine getdt
!=============================================================================
subroutine getcurrent(ixI^L,ixO^L,w,x,primvar,current)

! j = q (E x B)/B^2 + kappa_par Epar + kappa_perp Eperp
use mod_global_parameters

integer, intent(in)                       :: ixO^L, ixI^L
double precision, intent(in)              :: w(ixI^S,1:nw), x(ixI^S,1:ndim)
logical,intent(in)                        :: primvar
double precision, intent(out)             :: current(ixG^T,1:ndir)
! .. local ..
double precision, dimension(ixG^T)        :: vi, b2, e2, EdotB
double precision, dimension(ixG^T,1:ndir) :: Epar, Eperp
!-----------------------------------------------------------------------------

b2(ixO^S) = ({^C& w(ixO^S,b^C_)**2|+})
e2(ixO^S) = ({^C& w(ixO^S,e^C_)**2|+})
EdotB(ixO^S) = {^C& w(ixO^S,e^C_)*w(ixO^S,b^C_) |+}

{^C&
call getv(w,x,ixI^L,ixO^L,^C,vi) ! drift velocity
where (b2(ixO^S) .lt. smalldouble .or. e2(ixO^S) .gt. b2(ixO^S))
   current(ixO^S,^C) = eqpar(kpar_) * w(ixO^S,e^C_)
elsewhere
   Epar(ixO^S,^C)  = EdotB(ixO^S)*w(ixO^S,b^C_)/b2(ixO^S)
   Eperp(ixO^S,^C) = w(ixO^S,e^C_) - Epar(ixO^S,^C)
   current(ixO^S,^C)  = w(ixO^S,q_) * vi(ixO^S) &
        + eqpar(kpar_) * Epar(ixO^S,^C) &
        + eqpar(kperp_) * Eperp(ixO^S,^C)
end where
\}

end subroutine getcurrent
!=============================================================================
!=============================================================================
! end module amrvacphys/- ff
!#############################################################################
