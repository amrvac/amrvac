!============================================================================
subroutine addsource_grav(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

! w[iw]=w[iw]+qdt*S[wCT,qtC,x] where S is the source based on wCT within ixO

use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L, iw^LIM
double precision, intent(in) :: qdt, qtC, qt, x(ixI^S,1:ndim)
double precision, intent(inout) :: wCT(ixI^S,1:nw), w(ixI^S,1:nw)

integer :: iw, idims
!-----------------------------------------------------------------------------
! add sources from gravity

do idims = 1, ndim
   if (mom(idims) >= iwmin .and. mom(idims) <= iwmax) then
      w(ixO^S,mom(idims))=w(ixO^S,mom(idims)) &
           +qdt*grav_(idims)*wCT(ixO^S,rho_)
   end if

   if (hd_energy .and. e_ >= iwmin .and. e_ <= iwmax) then
      w(ixO^S,e_)=w(ixO^S,e_) &
           +qdt*grav_(idims)*wCT(ixO^S,mom(idims))
   end if
end do

end subroutine addsource_grav
!=============================================================================
subroutine getdt_grav(w,ixG^L,ix^L,dtnew,dx^D,x)

use mod_global_parameters

integer, intent(in) :: ixG^L, ix^L
double precision, intent(in) :: dx^D, x(ixG^S,1:ndim), w(ixG^S,1:nw)
double precision, intent(inout) :: dtnew

double precision:: dxinv(1:ndim), dtgrav
integer:: idims
!----------------------------------------------------------------------------

^D&dxinv(^D)=one/dx^D;
dtgrav=bigdouble
do idims=1,ndim
   if(abs(grav_(idims))>zero)&
   dtgrav=min(dtgrav,one/sqrt(abs(grav_(idims))*dxinv(idims)))
enddo

dtnew=dtgrav

end subroutine getdt_grav
!=============================================================================
