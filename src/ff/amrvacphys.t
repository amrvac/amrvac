!#############################################################################
! Module amrvacphys/- ff
! Maxwells equations plus Ohm's law according to Komissarov, Barkov and Lyutikov 2007
! 2015-08-12 by Oliver Porth

INCLUDE:amrvacnul/roe.t
INCLUDE:amrvacnul/hllc.t
!=============================================================================
subroutine checkglobaldata
include 'amrvacdef.f'
!-----------------------------------------------------------------------------

! Check if ssplitresis = .true.
if (ssplitresis .neqv. .true.) &
     call mpistop('Please run with ssplitresis = .true.')

! We require three vector components
if (^NC .ne. 3) &
     call mpistop('FF needs three vector components.')

end subroutine checkglobaldata
!=============================================================================
subroutine initglobaldata

! place to set entropy fixes etc, absent for now

include 'amrvacdef.f'
!-----------------------------------------------------------------------------
eqpar(kpar_)   = 1.0d-2
eqpar(kperp_)  = 1.0d-6
eqpar(kappa_)  = 0.5d0

end subroutine initglobaldata
!=============================================================================
subroutine checkw(checkprimitive,ixI^L,ixO^L,w,flag)

include 'amrvacdef.f'
  
logical, intent(in)          :: checkprimitive
integer, intent(in)          :: ixI^L,ixO^L
double precision, intent(in) :: w(ixI^S,1:nw)
logical, intent(out)         :: flag(ixG^T)
!-----------------------------------------------------------------------------

flag(ixG^T)=.true.

end subroutine checkw
!=============================================================================
subroutine conserve(ixI^L,ixO^L,w,x,patchw)

include 'amrvacdef.f'

integer, intent(in)               :: ixI^L,ixO^L
double precision, intent(inout)   :: w(ixI^S,nw)
logical, intent(in)               :: patchw(ixG^T)
double precision, intent(in)      :: x(ixI^S,1:ndim)
!-----------------------------------------------------------------------------


end subroutine conserve
!=============================================================================
subroutine ecrossb(ixI^L,ixO^L,idir,w,patchw,res)

include 'amrvacdef.f'

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

include 'amrvacdef.f'

integer, intent(in)                :: ixI^L,ixO^L
double precision, intent(inout)    :: w(ixI^S,1:nw)
double precision, intent(in)       :: x(ixI^S,1:ndim)
!-----------------------------------------------------------------------------

end subroutine primitive
!==============================================================================
subroutine e_to_rhos(ixI^L,ixO^L,w,x)

include 'amrvacdef.f'
integer:: ixI^L,ixO^L
double precision:: w(ixI^S,nw)
double precision, intent(in)      :: x(ixI^S,1:ndim)
!-------------------------------------------------------------------------

call mpistop("e to rhos for FF unavailable")

end subroutine e_to_rhos
!=============================================================================
subroutine rhos_to_e(ixI^L,ixO^L,w,x)

include 'amrvacdef.f'

integer:: ixI^L,ixO^L
double precision:: w(ixI^S,nw)
double precision, intent(in)      :: x(ixI^S,1:ndim)
!-----------------------------------------------------------------------------

call mpistop("e to rhos for FF unavailable")

end subroutine rhos_to_e
!=============================================================================
subroutine ppmflatcd(ixI^L,ixO^L,ixL^L,ixR^L,w,d2w,drho,dp)

include 'amrvacdef.f'

integer, intent(in)           :: ixI^L,ixO^L,ixL^L,ixR^L
double precision, intent(in)  :: w(ixI^S,nw),d2w(ixG^T,1:nwflux)

double precision, intent(inout) :: drho(ixG^T),dp(ixG^T)

!-----------------------------------------------------------------------------

call mpistop("PPM with flatsh=.true. can not be used with FF !")


end subroutine ppmflatcd
!=============================================================================
subroutine ppmflatsh(ixI^L,ixO^L,ixLL^L,ixL^L,ixR^L,ixRR^L,idims,w,drho,dp,dv)

include 'amrvacdef.f'

integer, intent(in)           :: ixI^L,ixO^L,ixLL^L,ixL^L,ixR^L,ixRR^L
integer, intent(in)           :: idims
double precision, intent(in)  :: w(ixI^S,nw)

double precision, intent(inout) :: drho(ixG^T),dp(ixG^T),dv(ixG^T)

!-----------------------------------------------------------------------------

call mpistop("PPM with flatsh=.true. can not be used with FF !")

end subroutine ppmflatsh
!=============================================================================
subroutine getv(w,x,ixI^L,ixO^L,idims,v)

! Drift velocity
! v = E x B /B^2 

include 'amrvacdef.f'
  
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
where( b2(ixO^S) <= smalldouble**2)
   v(ixO^S) = zero
elsewhere
   v(ixO^S) = ecrossbi(ixO^S) / b2(ixO^S)
endwhere

end subroutine getv
!=============================================================================
subroutine getcmax(new_cmax,w,x,ixI^L,ixO^L,idims,cmax,cmin,needcmin)
  
! Calculate cmax_idim within ixO^L
! new_cmax is not used
  
include 'amrvacdef.f'
  
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
include 'amrvacdef.f'

integer, intent(in)                :: ixI^L,ixO^L,iw,idims
double precision, intent(in)       :: w(ixI^S,nw)
double precision, intent(in)       :: x(ixI^S,1:ndim)
double precision, intent(out)      :: f(ixG^T)
logical, intent(out)               :: transport
! .. local ..
integer                            :: k
double precision, dimension(ixG^T) :: vidims, vc, p, sqrE, sqrB, ecrossbi
!-----------------------------------------------------------------------------
transport=.false.

if(iw==phib_) then
{#IFNDEF FCT
       f(ixO^S)=w(ixO^S,b0_+idims)
}{#IFDEF FCT
       f(ixO^S)=zero
}
else

!f^i[B_c] = lvc(c,i,k) Ek + phi delta(i,c)
{if (iw==b^C_) then
   if (idims .ne. ^C) then
      do k=1,^NC
         if (idims .eq. k .or. ^C .eq. k) cycle
         f(ixO^S) = lvc(^C,idims,k) * w(ixO^S,e0_+k)
      end do
   else
      f(ixO^S)=w(ixO^S,phib_)
   end if
end if\}

!f^i[E_c] = - lvc(c,i,k) Bk
{if (iw==e^C_) then
   if (idims .ne. ^C) then 
      do k=1,^NC
         if (idims .eq. k .or. ^C .eq. k) cycle
         f(ixO^S) = - lvc(^C,idims,k) * w(ixO^S,b0_+k)
      end do
   end if
end if\}

end if

end subroutine getflux
!=============================================================================
subroutine getfluxforhllc(w,x,ixI^L,ixO^L,iw,idims,f,transport)

! Calculate non-transport flux f_idim[iw] within ixO^L.

include 'amrvacdef.f'

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

include 'amrvacdef.f'

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

include 'amrvacdef.f'

integer, intent(in)             :: ixI^L, ixO^L, iw^LIM
double precision, intent(in)    :: qdt, qtC, qt
double precision, intent(in)    :: x(ixI^S,1:ndim)
double precision, intent(in)    :: wCT(ixI^S,1:nw)
double precision, intent(inout) :: w(ixI^S,1:nw)
logical, intent(in)             :: qsourcesplit
! .. local ..
!double precision                :: wtmp(ixI^S,1:nw)
double precision                :: dx^D
!-----------------------------------------------------------------------------
dx^D=dxlevel(^D);

! Two sources: Sb is added via Strang-splitting, Sa is unsplit.
! Using ssplitresis to distinguish (needs to be set true!!!)
if(qsourcesplit .eqv. ssplitresis) then
   call addsource_b(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,dx^D)
else
   call addsource_a(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,dx^D)
end if


! now update the auxiliaries to the new state
call getaux(.true.,w,x,ixI^L,ixO^L,'addsource')

end subroutine addsource
!=============================================================================
subroutine addsource_a(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,dx^D)
include 'amrvacdef.f'

integer, intent(in) :: ixI^L, ixO^L, iw^LIM
double precision, intent(in) :: qdt, qtC, qt, wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
double precision, intent(in) :: dx^D
double precision, intent(inout) :: w(ixI^S,1:nw)
! .. local ..
double precision, dimension(ixG^T)   :: vidir, divE
double precision                     :: evec(ixG^T,1:ndir)
!-----------------------------------------------------------------------------

! S[Ei_] = - divE (E x B)/B^2
evec(ixI^S,1:ndir)=w(ixI^S,e0_+1:e0_+ndir)
select case(typediv)
case("central")
   call divvector(evec,ixI^L,ixO^L,divE)
case("limited")
   call divvectorS(evec,ixI^L,ixO^L,divE)
end select


{^C&
call getv(wCT,x,ixI^L,ixO^L,^C,vidir)
w(ixO^S,e^C_) = w(ixO^S,e^C_) - qdt * divE(ixO^S) * vidir(ixO^S)
\}

end subroutine addsource_a
!=============================================================================
subroutine addsource_b(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,dx^D)
include 'amrvacdef.f'

integer, intent(in) :: ixI^L, ixO^L, iw^LIM
double precision, intent(in) :: qdt, qtC, qt, wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
double precision, intent(in) :: dx^D
double precision, intent(inout) :: w(ixI^S,1:nw)
! .. local ..
double precision, dimension(ixG^T,1:ndir)  :: Epar, Eperp
double precision, dimension(ixG^T)         :: babs
!-----------------------------------------------------------------------------

! S[phib_] = -kappa phib
w(ixO^S,phib_) = w(ixO^S,phib_)*exp(-qdt*eqpar(kappa_))

! S[E] = - kappa_par Epar - kappa_perp Eperp
! Split the source in parallel and perpendicular components
! Then solve the linear evolution equation (assuming constant background state)
! Then add again both components.


babs(ixO^S) = sqrt({^C& w(ixO^S,b^C_)**2|+})

{^C&
Epar(ixO^S,^C)  = w(ixO^S,e^C_)*w(ixO^S,b^C_)/babs(ixO^S)
Eperp(ixO^S,^C) = w(ixO^S,e^C_) - Epar(ixO^S,^C)
\}

! Handle zero magnetic field separate:
! Use perpendicular conductivity for all directions
where (babs(ixO^S) .lt. smalldouble)
{^C&
     w(ixO^S,e^C_) = w(ixO^S,e^C_) * exp(-eqpar(kperp_)*qdt)
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

! Limit the timestp due to the resistive source addition
  
include 'amrvacdef.f'

integer, intent(in)              :: ixG^L, ix^L
double precision, intent(in)     :: dx^D, x(ixG^S,1:ndim)
double precision, intent(inout)  :: w(ixG^S,1:nw), dtnew
!-----------------------------------------------------------------------------

dtnew = dtdiffpar / max(eqpar(kpar_),eqpar(kperp_))

end subroutine getdt
!=============================================================================
subroutine getcurrent(ixI^L,ixO^L,w,x,primvar,current)

! j = divE (E x B)/B^2 + kappa_par Epar + kappa_perp Eperp
include 'amrvacdef.f'

integer, intent(in)                       :: ixO^L, ixI^L
double precision, intent(in)              :: w(ixI^S,1:nw), x(ixI^S,1:ndim)
logical,intent(in)                        :: primvar
double precision, intent(out)             :: current(ixG^T,1:ndir)
! .. local ..
double precision, dimension(ixG^T)        :: vi, babs
double precision, dimension(ixG^T)        :: divE
double precision                          :: evec(ixG^T,1:ndir)
double precision, dimension(ixG^T,1:ndir) :: Epar, Eperp
!-----------------------------------------------------------------------------

babs(ixO^S) = sqrt({^C& w(ixO^S,b^C_)**2|+})

evec(ixI^S,1:ndir)=w(ixI^S,e0_+1:e0_+ndir)
select case(typediv)
case("central")
   call divvector(evec,ixI^L,ixO^L,divE)
case("limited")
   call divvectorS(evec,ixI^L,ixO^L,divE)
end select

{^C&
call getv(w,x,ixI^L,ixO^L,^C,vi) ! drift velocity
where (babs(ixO^S) .lt. smalldouble)
     current(ixO^S,^C) = eqpar(kperp_) * w(ixO^S,e^C_)
elsewhere
     Epar(ixO^S,^C)  = w(ixO^S,e^C_)*w(ixO^S,b^C_)/babs(ixO^S)
     Eperp(ixO^S,^C) = w(ixO^S,e^C_) - Epar(ixO^S,^C)
     current(ixO^S,^C)  = divE(ixO^S) * vi(ixO^S) &
          + eqpar(kpar_) * Epar(ixO^S,^C) &
          + eqpar(kperp_) * Eperp(ixO^S,^C)
end where
\}

end subroutine getcurrent
!=============================================================================
!=============================================================================
! end module amrvacphys/- ff
!#############################################################################
