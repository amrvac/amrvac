{#IFDEF GLM
!=============================================================================
subroutine addsource_glm1(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,dx^D)

! Add divB related sources to w within ixO
! corresponding to Dedner JCP 2002, 175, 645 _equation 24_
! giving the EGLM scheme
use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L, iw^LIM
double precision, intent(in) :: qdt, qtC, qt, wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
double precision, intent(in) :: dx^D
double precision, intent(inout) :: w(ixI^S,1:nw)
double precision:: divb(ixG^T),v(ixG^T)
!-----------------------------------------------------------------------------
! We calculate now div B
call getdivb(wCT,ixI^L,ixO^L,divb)

! Psi = Psi - qdt Ch^2/Cp^2 Psi
if (eqpar(Cr_) < zero) then
  w(ixO^S,psi_) = abs(eqpar(Cr_))*w(ixO^S,psi_)
else 
  ! implicit update of psi variable
  w(ixO^S,psi_) = dexp(-qdt*(storeeqparch/eqpar(Cr_)))*w(ixO^S,psi_)
end if

! b = b - qdt v * div b
!{^C&
!call getv(wCT,x,ixI^L,ixO^L,^C,v)
!w(ixO^S,b^C_)=w(ixO^S,b^C_)-qdt*v(ixO^S)*divb(ixO^S)\}

{^C&
   w(ixO^S,s^C_)=w(ixO^S,s^C_)-qdt*wCT(ixO^S,b^C_)*divb(ixO^S)\}

end subroutine addsource_glm1
!=============================================================================
subroutine addsource_glm2(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,dx^D)

! w[iws]=w[iws]+qdt*S[iws,wCT] where S is the source based on wCT within ixO

use mod_global_parameters

integer, intent(in)             :: ixI^L, ixO^L, iw^LIM
double precision, intent(in)    :: dx^D
double precision, intent(in)    :: qdt, qtC, qt
double precision, intent(in)    :: x(ixI^S,1:ndim)
double precision, intent(in)    :: wCT(ixI^S,1:nw)
double precision, intent(inout) :: w(ixI^S,1:nw)


integer, parameter                :: idirmin0=7-2*ndir
integer                           :: idims,idir,jdir,kdir
double precision,dimension(ixG^T) :: divb,Psi,GradPsi,tmp
double precision                  :: Efield(ixG^T,idirmin0:3),vCT(ixG^T,1:^NC)
!-----------------------------------------------------------------------------

! Eq. 45 Dedner JCP 2002, 175, 645
if (eqpar(Cr_) < zero) then
  w(ixO^S,psi_) = abs(eqpar(Cr_))*w(ixO^S,psi_)
else 
  ! implicit update of psi variable
  w(ixO^S,psi_) = dexp(-qdt*(storeeqparch/eqpar(Cr_)))*w(ixO^S,psi_)
end if

! We calculate now auxiliary for wCT, for getv call following
call getaux(.true.,wCT,x,ixI^L,ixO^L,'adddsource_glm2')

! Eq. 38 Dedner JCP 2002, 175, 645
 ! speed vCT=
 Psi(ixO^S)= ({^C&wCT(ixO^S,s^C_)*wCT(ixO^S,b^C_)+})/wCT(ixO^S,xi_) ! VdotB
 tmp(ixO^S) = {^C&wCT(ixO^S,b^C_)**2.0d0+} ! is sqrB
 do idims=1,ndir
  where(wCT(ixO^S,xi_)+tmp(ixO^S)<=smallxi)
    vCT(ixO^S,idims)=zero
  elsewhere
    vCT(ixO^S,idims)=(wCT(ixO^S,s0_+idims)+Psi(ixO^S)*wCT(ixO^S,b0_+idims))/ &
             (wCT(ixO^S,xi_)+tmp(ixO^S))
  endwhere
 end do
 ! Electric field E=B x v
 Efield = zero
 do idir=idirmin0,3; do jdir=1,ndir; do kdir=1,ndir
  if(lvc(idir,jdir,kdir)/=0)then
     tmp(ixO^S)=wCT(ixO^S,b0_+jdir)*vCT(ixO^S,kdir)
     if(lvc(idir,jdir,kdir)==1)then
         Efield(ixO^S,idir)=Efield(ixO^S,idir)+tmp(ixO^S)
     else
         Efield(ixO^S,idir)=Efield(ixO^S,idir)-tmp(ixO^S)
     endif
  endif
 enddo; enddo; enddo
 ! gradient of Psi
 Psi(ixI^S)=wCT(ixI^S,psi_)
 do idims=1,ndim
   select case(typegrad)
    case("central")
     call gradient(Psi,ixO^L,idims,gradPsi)
    case("limited")
     call gradientS(Psi,ixO^L,idims,gradPsi)
   end select
   ! S_new=S_old+qdt*(grad(Psi) x E)
   do idir=1,ndir; do kdir=idirmin0,3
    if (lvc(idir,idims,kdir)/=0) then
      tmp(ixO^S)=qdt*gradPsi(ixO^S)*Efield(ixO^S,kdir)
     if(lvc(idir,idims,kdir)==1)then
       w(ixO^S,s0_+idir)=w(ixO^S,s0_+idir)+tmp(ixO^S)
     else
       w(ixO^S,s0_+idir)=w(ixO^S,s0_+idir)-tmp(ixO^S)
     end if
    end if
   end do; end do
   
   ! Psi=Psi-qdt (v . grad(Psi))
   w(ixO^S,psi_) = w(ixO^S,psi_)&
                   -qdt*vCT(ixO^S,idims)*gradPsi(ixO^S)
 end do

end subroutine addsource_glm2
}
!=============================================================================
subroutine addsource_powel(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,dx^D)

! Add Powell's divB related sources to w within ixO

use mod_global_parameters

integer, intent(in)             :: ixI^L, ixO^L, iw^LIM
double precision, intent(in)    :: qdt, qtC, qt
double precision, intent(in)    :: wCT(ixI^S,1:nw)
double precision, intent(in)    :: x(ixI^S,1:ndim)
double precision, intent(in)    :: dx^D
double precision, intent(inout) :: w(ixI^S,1:nw)

integer :: iw
double precision:: divb(ixG^T),v(ixG^T)
!-----------------------------------------------------------------------------

! We calculate div B
call getdivb(wCT,ixI^L,ixO^L,divb)
divb(ixO^S)=qdt*divb(ixO^S)

! We calculate now auxiliary for wCT, for getv call following
call getaux(.true.,wCT,x,ixI^L,ixO^L,'addsource_divb')

if (typedivbfix=='powel') then
   do iw= iw^LIM
      select case(iw)
      {case(s^C_)
         w(ixO^S,iw)=w(ixO^S,iw)-wCT(ixO^S,b^C_)*divb(ixO^S)\}
      {case(b^C_)
         call getv(wCT,x,ixI^L,ixO^L,^C,v)
         w(ixO^S,iw)=w(ixO^S,iw)-v(ixO^S)*divb(ixO^S)\}
      end select
   end do
end if

end subroutine addsource_powel
!=============================================================================
subroutine addsource_janhunen(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,dx^D)

! Add divB related sources to w within ixO
! corresponding to Janhunen, just the term in the induction equation.

use mod_global_parameters

integer, intent(in)             :: ixI^L, ixO^L, iw^LIM
double precision, intent(in)    :: qdt, qtC, qt
double precision, intent(in)    :: wCT(ixI^S,1:nw)
double precision, intent(in)    :: x(ixI^S,1:ndim)
double precision, intent(in)    :: dx^D
double precision, intent(inout) :: w(ixI^S,1:nw)

integer :: iw
double precision:: divb(ixG^T),v(ixG^T)
!-----------------------------------------------------------------------------

! We calculate div B
call getdivb(wCT,ixI^L,ixO^L,divb)
divb(ixO^S)=qdt*divb(ixO^S)

! We calculate now auxiliary for wCT, for getv call following
call getaux(.true.,wCT,x,ixI^L,ixO^L,'addsource_divb')

if(typedivbfix=='janhunen')then
   do iw= iw^LIM
      select case(iw)
      {case(b^C_)
         call getv(wCT,x,ixI^L,ixO^L,^C,v)
         w(ixO^S,iw)=w(ixO^S,iw)-v(ixO^S)*divb(ixO^S)\}
      end select
   end do
end if

end subroutine addsource_janhunen
!=============================================================================
subroutine addsource_linde(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,dx^D)

! Add Linde's divB related sources to wnew within ixO

use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L, iw^LIM
double precision, intent(in) :: qdt, qtC, qt, x(ixI^S,1:ndim)
double precision, intent(in) :: wCT(ixI^S,1:nw)
double precision, intent(in) :: dx^D
double precision, intent(inout) :: w(ixI^S,1:nw)

integer :: iw, idims, ix^L
double precision :: divb(ixG^T),graddivb(ixG^T)
!-----------------------------------------------------------------------------

! Calculate div B
ix^L=ixO^L^LADD1;
call getdivb(wCT,ixI^L,ix^L,divb)

! Add Linde's diffusive terms, but only to the induction equation
!
do idims=1,ndim
   ! Calculate grad_idim(divb)
   select case(typegrad)
   case("central")
     call gradient(divb,ixO^L,idims,graddivb)
   case("limited")
     call gradientS(divb,ixO^L,idims,graddivb)
   end select

   ! Multiply by Linde's eta*dt = divbdiff*(c_max*dx)*dt = divbdiff*dx**2
   if (slab) then
      graddivb(ixO^S)=graddivb(ixO^S)*divbdiff/(^D&1.0d0/dxlevel(^D)**2+)
   else
      graddivb(ixO^S)=graddivb(ixO^S)*divbdiff &
                      /(^D&1.0d0/mygeo%dx(ixO^S,^D)**2+)
   end if
   do iw= iw^LIM
      if (iw==b0_+idims) then
         ! B_idim += eta*grad_idim(divb)
         w(ixO^S,iw)=w(ixO^S,iw)+graddivb(ixO^S)
      end if
   end do
end do

end subroutine addsource_linde
!=============================================================================
subroutine getdivb(w,ixI^L,ixO^L,divb)

! Calculate div B within ixO

use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L
double precision, intent(in) :: w(ixI^S,1:nw)
double precision :: divb(ixG^T)

double precision :: bvec(ixG^T,1:ndir)
!-----------------------------------------------------------------------------

bvec(ixI^S,1:ndir)=w(ixI^S,b0_+1:b0_+ndir)
select case(typediv)
case("central")
  call divvector(bvec,ixI^L,ixO^L,divb)
case("limited")
  call divvectorS(bvec,ixI^L,ixO^L,divb)
end select

end subroutine getdivb
!=============================================================================
