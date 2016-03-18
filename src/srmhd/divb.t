{#IFDEF GLM
!=============================================================================
subroutine addsource_glm1(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,dx^D)

! Add divB related sources to w within ixO
! corresponding to Dedner JCP 2002, 175, 645 _equation 24_
! giving the EGLM scheme
include 'amrvacdef.f'

integer, intent(in) :: ixI^L, ixO^L, iw^LIM
double precision, intent(in) :: qdt, qtC, qt, wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
double precision, intent(in) :: dx^D
double precision, intent(inout) :: w(ixI^S,1:nw)
double precision:: divb(ixG^T)
!-----------------------------------------------------------------------------
{#IFDEF FCT
! When using FCT, we treat this as a source term
! We calculate now div B consistent with FCT
call getdivb(wCT,ixI^L,ixO^L,divb)
w(ixO^S,psi_) = wCT(ixO^S,psi_) -  cmax_global**2 * qdt * divb(ixO^S) 
}


! Psi = Psi - qdt Ch^2/Cp^2 Psi
if (eqpar(Cr_) < zero) then
  w(ixO^S,psi_) = abs(eqpar(Cr_))*w(ixO^S,psi_)
else 
  ! implicit update of psi variable
  w(ixO^S,psi_) = dexp(-qdt*(cmax_global/eqpar(Cr_)))*w(ixO^S,psi_)
end if

end subroutine addsource_glm1
!=============================================================================
subroutine addsource_glm2(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,dx^D)

! w[iws]=w[iws]+qdt*S[iws,wCT] where S is the source based on wCT within ixO

include 'amrvacdef.f'

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
{#IFDEF FCT
! When using FCT, we treat this as a source term
! We calculate now div B consistent with FCT
call getdivb(wCT,ixI^L,ixO^L,divb)
w(ixO^S,psi_) = wCT(ixO^S,psi_) -  cmax_global**2 * qdt * divb(ixO^S) 
}


! Eq. 45 Dedner JCP 2002, 175, 645
if (eqpar(Cr_) < zero) then
  w(ixO^S,psi_) = abs(eqpar(Cr_))*w(ixO^S,psi_)
else 
  ! implicit update of psi variable
  w(ixO^S,psi_) = dexp(-qdt*(cmax_global/eqpar(Cr_)))*w(ixO^S,psi_)
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
     call gradient(Psi,ixI^L,ixO^L,idims,gradPsi)
    case("limited")
     call gradientS(Psi,ixI^L,ixO^L,idims,gradPsi)
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
   ! e  =e  -qdt (b . grad(Psi))
   ! VERIFIED: same for split strategy B0+B1 (B0field=T has B1.grad(Psi))
   w(ixO^S,e_) = w(ixO^S,e_)&
                      -qdt*wCT(ixO^S,b0_+idims)*gradPsi(ixO^S)
 end do

end subroutine addsource_glm2
}
!=============================================================================
subroutine addsource_powel(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,dx^D)

! Add Powell's divB related sources to w within ixO

include 'amrvacdef.f'

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
      case(e_)
         w(ixO^S,iw)=w(ixO^S,iw)-&
             ({^C&wCT(ixO^S,s^C_)*wCT(ixO^S,b^C_)+})/wCT(ixO^S,xi_)*divb(ixO^S)
      end select
   end do
end if

end subroutine addsource_powel
!=============================================================================
subroutine addsource_janhunen(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,dx^D)

! Add divB related sources to w within ixO
! corresponding to Janhunen, just the term in the induction equation.

include 'amrvacdef.f'

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
   do iw= b1_,b^ND_  !iw^LIM
      select case(iw)
      {case(b^D_)
         call getv(wCT,x,ixI^L,ixO^L,^D,v)
         w(ixO^S,iw)=w(ixO^S,iw)-v(ixO^S)*divb(ixO^S)\}
      end select
   end do
end if

end subroutine addsource_janhunen
!=============================================================================
subroutine addsource_linde(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,dx^D)

! Add Linde's divB related sources to wnew within ixO
include 'amrvacdef.f'

integer, intent(in)             :: ixI^L, ixO^L, iw^LIM
double precision, intent(in)    :: qdt, qtC, qt, x(ixI^S,1:ndim)
double precision, intent(in)    :: dx^D
double precision, intent(inout) :: wCT(ixI^S,1:nw), w(ixI^S,1:nw)
integer :: iw, idims, ix^L, ixp^L, i^D, iside
double precision :: divb(ixI^S),graddivb(ixI^S)
!-----------------------------------------------------------------------------

! Calculate div B
ix^L=ixO^L^LADD1;
call getdivb(wCT,ixI^L,ix^L,divb)
! for AMR stability, retreat one cell layer from the boarders of level jump
ixp^L=ixO^L;
do idims=1,ndim
  select case(idims)
   {case(^D)
      do iside=1,2
        i^DD=kr(^DD,^D)*(2*iside-3);
        if(leveljump(i^DD)) then
          if(iside==1) then
            ixpmin^D=ixOmin^D-i^D
          else
            ixpmax^D=ixOmax^D-i^D
          end if 
        end if
      end do
   \}
  end select
end do

! Add Linde's diffusive terms
do idims=1,ndim
   ! Calculate grad_idim(divb)
   select case(typegrad)
   case("central")
     call gradient(divb,ixI^L,ixp^L,idims,graddivb)
   case("limited")
     call gradientS(divb,ixI^L,ixp^L,idims,graddivb)
   end select
   !ixmin^D=ixpmin^D;
   !ixmax^D=merge(ixmin^D,ixmax^D,kr(idims,^D)==1);
   !call gradient(divb,ixI^L,ix^L,idims,graddivb)
   !ixmin^D=merge(ixmax^D,ixmin^D,kr(idims,^D)==1);
   !ixmax^D=ixpmax^D;
   !call gradient(divb,ixI^L,ix^L,idims,graddivb)
   !ixmin^D=ixpmin^D+kr(idims,^D);
   !ixmax^D=ixpmax^D-kr(idims,^D);
   !call gradientS(divb,ixI^L,ix^L,idims,graddivb)

   ! Multiply by Linde's eta*dt = divbdiff*(c_max*dx)*dt = divbdiff*dx**2
   if (slab) then
      graddivb(ixp^S)=graddivb(ixp^S)*divbdiff/(^D&1.0d0/dxlevel(^D)**2+)
   else
      graddivb(ixp^S)=graddivb(ixp^S)*divbdiff &
                      /(^D&1.0d0/mygeo%dx(ixp^S,^D)**2+)
   end if
   do iw= iw^LIM
      if (iw==b0_+idims) then
         ! B_idim += eta*grad_idim(divb)
         w(ixp^S,iw)=w(ixp^S,iw)+graddivb(ixp^S)
{#IFDEF ENERGY
      else if (iw==e_ .and. typedivbdiff=='all') then
         ! e += B_idim*eta*grad_idim(divb)
         w(ixp^S,iw)=w(ixp^S,iw)+wCT(ixp^S,b0_+idims)*graddivb(ixp^S)
}
      end if
   end do
end do

end subroutine addsource_linde
!=============================================================================
subroutine getdivb(w,ixI^L,ixO^L,divb)

! Calculate div B within ixO

include 'amrvacdef.f'

integer, intent(in)                :: ixI^L, ixO^L
double precision, intent(in)       :: w(ixI^S,1:nw)
double precision                   :: divb(ixG^T)

double precision                   :: bvec(ixG^T,1:ndir)

{#IFDEF FCT
integer                            :: ixC^L, idir, idim, ixJp^L, ic^D, ix^L
integer                            :: ixKp^L, ixJpKp^L, ixJm^L, ixJmKm^L
double precision                   :: divb_corner(ixG^T), sign
}
!-----------------------------------------------------------------------------

bvec(ixI^S,1:ndir)=w(ixI^S,b0_+1:b0_+ndir)

{#IFNDEF FCT
select case(typediv)
case("central")
  call divvector(bvec,ixI^L,ixO^L,divb)
case("limited")
  call divvectorS(bvec,ixI^L,ixO^L,divb)
end select
}{#IFDEF FCT

! For fct, we calculate the divB on the corners according to Toth (2000), 
! eq. (27) and average to the cell centers for output.

{ixCmax^D=ixOmax^D;}
{ixCmin^D=ixOmin^D-1;} ! Extend range by one


! Get the corner-centered divb:

divb_corner(ixC^S) = zero
do idir = 1, ndim ! idir is the component of the field to consider (j)
{do ic^DB=0,1\}
   {ix^L=ixC^L+ic^D;}

   select case(idir)
      {^D&   
           case(^D)
      sign = dble(ic^D*2 - 1)
      \}
   end select

   divb_corner(ixC^S) = divb_corner(ixC^S) &
        + sign * bvec(ix^S,idir)/dxlevel(idir)
   
   {end do\}
end do
divb_corner(ixC^S) = divb_corner(ixC^S) / 2.0d0**(ndim-1)



! Now average back to the cell centers:

divb(ixO^S) = zero
{do ic^DB=-1,0\}
   {ixC^L=ixO^L+ic^D;}

   divb(ixO^S) = divb(ixO^S) &
        + divb_corner(ixC^S)

{end do\}
divb(ixO^S) = divb(ixO^S) / 2.0d0**(ndim)

}
end subroutine getdivb
!=============================================================================
