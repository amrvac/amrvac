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
do iw= iw^LIM
   select case (iw)
   case (m^D_)
      ! dm_i/dt= +rho*g_i
      idims=iw-m0_
      if (abs(eqpar(grav0_+idims))>smalldouble) &
          w(ixO^S,m0_+idims)=w(ixO^S,m0_+idims) &
              +qdt*eqpar(grav0_+idims)*wCT(ixO^S,rho_)
{#IFDEF ENERGY
   case (e_)
      ! de/dt= +g_i*m_i
      do idims=1,ndim
         if (abs(eqpar(grav0_+idims))>smalldouble) &
            w(ixO^S,ee_)=w(ixO^S,ee_) &
              +qdt*eqpar(grav0_+idims)*wCT(ixO^S,m0_+idims)
      end do
}
   end select
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
   if(abs(eqpar(grav0_+idims))>zero)&
   dtgrav=min(dtgrav,one/sqrt(abs(eqpar(grav0_+idims))*dxinv(idims)))
enddo

dtnew=dtgrav

end subroutine getdt_grav
!=============================================================================
