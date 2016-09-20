!=============================================================================
! amrvacusr.t.arfff
!=============================================================================
INCLUDE:amrvacnul/usrflags.t
INCLUDE:amrvacmodules/fff.t
!=============================================================================
subroutine initglobaldata_usr
use mod_global_parameters

double precision :: miu0,k_B
integer :: ix1,ix2,ix3,nxbc^D,ixpemin1,ixpemax1
logical, save :: firstusrglobaldata=.true.
!-----------------------------------------------------------------------------
! CGS Unit
k_B=1.3806d-16     ! erg*K^-1
miu0=4.d0*dpi     ! Gauss^2 cm^2 dyne^-1
Lunit=1.d9        ! cm
Teunit=1.d6       ! K
nHunit=1.d9       ! cm^-3
mHunit=1.67262d-24 ! g
runit=1.4d0*mHunit*nHunit ! 2.341668000000000E-015 g*cm^-3
punit=2.3d0*nHunit*k_B*Teunit ! 0.317538000000000 erg*cm^-3
Bunit=dsqrt(miu0*punit)      ! 1.99757357615242 Gauss
vunit=Bunit/dsqrt(miu0*runit) ! 1.16448846777562E007 cm/s = 116.45 km/s
tunit=Lunit/vunit ! 85.8746159942810 s

eqpar(eta_)=0.d0
eqpar(grav1_)=0.d0
eqpar(grav2_)=0.d0
eqpar(grav3_)=-2.74d4*Lunit/vunit**2

! Solar radius
SRadius=6.961d10/Lunit
! base density and temperature
rhob=1.d0
Tiso=1.d6/Teunit
eqpar(adiab_)=Tiso
if(firstusrglobaldata) then
  call init_b_fff_data('hmi20140328_140000.dat',Lunit,Bunit)
  firstusrglobaldata=.false.
endif

end subroutine initglobaldata_usr
!=============================================================================
subroutine initonegrid_usr(ixG^L,ix^L,w,x)
! initialize one grid

use mod_global_parameters

integer, intent(in) :: ixG^L,ix^L
double precision, intent(in) :: x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)

double precision :: Bf(ixG^S,1:ndir),qalpha,llift
integer :: ix^D,idim,idir,ixCo^D, jxCo^D, hxCo^D,na
logical :: patchw(ixG^T)
logical, save :: first=.true.
!----------------------------------------------------------------------------

if(first .and. mype==0 .and. iprob==1) then
  write(*,*)'extrapolating 3D coronal field'
  first=.false.
endif
! use green function method to extrapolate B
qalpha=0.d0
llift=0.d0
call calc_lin_fff(ixG^L,ix^L,Bf,x,qalpha,llift)
w(ix^S,b1_)=Bf(ix^S,1)
w(ix^S,b2_)=Bf(ix^S,2)
w(ix^S,b3_)=Bf(ix^S,3)
w(ix^S,m1_)=zero
w(ix^S,m2_)=zero
w(ix^S,m3_)=zero
w(ix^S,rho_)=rhob*dexp(eqpar(grav3_)*SRadius**2/Tiso*&
               (1.d0/SRadius-1.d0/(x(ix^S,3)+SRadius)))

end subroutine initonegrid_usr
!==============================================================================
subroutine specialbound_usr(qt,ixG^L,ixO^L,iw,iB,w,x)
! special boundary types, user defined

use mod_global_parameters

integer, intent(in) :: ixO^L, iw, iB, ixG^L
double precision, intent(in) :: qt, x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)

double precision :: tmp(ixG^T),coeffrho,invT(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
double precision :: delxdely,delxdelz,delydelx,delydelz,delzdelx,delzdely
integer :: ix^D,ixIM^L
!----------------------------------------------------------------------------
select case(iB)
 case(1)
   w(ixO^S,rho_)=w(ixOmax1+dixB:ixOmax1+1:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)
   w(ixO^S,v1_)=-w(ixOmax1+dixB:ixOmax1+1:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,m1_)/&
                 w(ixOmax1+dixB:ixOmax1+1:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)
   w(ixO^S,v2_)=-w(ixOmax1+dixB:ixOmax1+1:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,m2_)/&
                 w(ixOmax1+dixB:ixOmax1+1:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)
   w(ixO^S,v3_)=-w(ixOmax1+dixB:ixOmax1+1:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,m3_)/&
                 w(ixOmax1+dixB:ixOmax1+1:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)
   do ix1=ixOmax1,ixOmin1,-1
     w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b1_:b3_)=(1.0d0/3.0d0)* &
                (-w(ix1+2,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b1_:b3_) &
           +4.0d0*w(ix1+1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b1_:b3_))
   enddo
   ! 2nd order CD for divB=0 to set normal B component better
   delxdely=dxlevel(1)/dxlevel(2)
   delxdelz=dxlevel(1)/dxlevel(3)
   do ix1=ixOmax1,ixOmin1,-1
     do ix2=ixOmin2+1,ixOmax2-1
       do ix3=ixOmin3+1,ixOmax3-1
         w(ix1,ix2,ix3,b1_)=w(ix1+2,ix2,ix3,b1_) &
          +delxdely*(w(ix1+1,ix2+1,ix3,b2_)-w(ix1+1,ix2-1,ix3,b2_))&
          +delxdelz*(w(ix1+1,ix2,ix3+1,b3_)-w(ix1+1,ix2,ix3-1,b3_))
       enddo
     enddo
   enddo
   call conserve(ixG^L,ixO^L,w,x,patchfalse)
 case(2)
   w(ixO^S,rho_)=w(ixOmin1-1:ixOmin1-dixB:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)
   w(ixO^S,v1_)=-w(ixOmin1-1:ixOmin1-dixB:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,m1_)/&
                 w(ixOmin1-1:ixOmin1-dixB:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)
   w(ixO^S,v2_)=-w(ixOmin1-1:ixOmin1-dixB:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,m2_)/&
                 w(ixOmin1-1:ixOmin1-dixB:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)
   w(ixO^S,v3_)=-w(ixOmin1-1:ixOmin1-dixB:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,m3_)/&
                 w(ixOmin1-1:ixOmin1-dixB:-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)
   do ix1=ixOmin1,ixOmax1
     w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b1_:b3_)=(1.0d0/3.0d0)* &
                (-w(ix1-2,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b1_:b3_) &
           +4.0d0*w(ix1-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b1_:b3_))
   enddo
   ! 2nd order CD for divB=0 to set normal B component better
   delxdely=dxlevel(1)/dxlevel(2)
   delxdelz=dxlevel(1)/dxlevel(3)
   do ix1=ixOmin1,ixOmax1
     do ix2=ixOmin2+1,ixOmax2-1
       do ix3=ixOmin3+1,ixOmax3-1
         w(ix1,ix2,ix3,b1_)=w(ix1-2,ix2,ix3,b1_) &
          -delxdely*(w(ix1-1,ix2+1,ix3,b2_)-w(ix1-1,ix2-1,ix3,b2_))&
          -delxdelz*(w(ix1-1,ix2,ix3+1,b3_)-w(ix1-1,ix2,ix3-1,b3_))
       enddo
     enddo
   enddo
   call conserve(ixG^L,ixO^L,w,x,patchfalse)
 case(3)
   w(ixO^S,rho_)=w(ixOmin1:ixOmax1,ixOmax2+dixB:ixOmax2+1:-1,ixOmin3:ixOmax3,rho_)
   w(ixO^S,v1_)=-w(ixOmin1:ixOmax1,ixOmax2+dixB:ixOmax2+1:-1,ixOmin3:ixOmax3,m1_)/&
                 w(ixOmin1:ixOmax1,ixOmax2+dixB:ixOmax2+1:-1,ixOmin3:ixOmax3,rho_)
   w(ixO^S,v2_)=-w(ixOmin1:ixOmax1,ixOmax2+dixB:ixOmax2+1:-1,ixOmin3:ixOmax3,m2_)/&
                 w(ixOmin1:ixOmax1,ixOmax2+dixB:ixOmax2+1:-1,ixOmin3:ixOmax3,rho_)
   w(ixO^S,v3_)=-w(ixOmin1:ixOmax1,ixOmax2+dixB:ixOmax2+1:-1,ixOmin3:ixOmax3,m3_)/&
                 w(ixOmin1:ixOmax1,ixOmax2+dixB:ixOmax2+1:-1,ixOmin3:ixOmax3,rho_)
   do ix2=ixOmax2,ixOmin2,-1
     w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,b1_:b3_)=(1.0d0/3.0d0)* &
                (-w(ixOmin1:ixOmax1,ix2+2,ixOmin3:ixOmax3,b1_:b3_) &
           +4.0d0*w(ixOmin1:ixOmax1,ix2+1,ixOmin3:ixOmax3,b1_:b3_))
   enddo
   !! 2nd order CD for divB=0 to set normal B component better
   delydelx=dxlevel(2)/dxlevel(1)
   delydelz=dxlevel(2)/dxlevel(3)
   do ix2=ixOmax2,ixOmin2,-1
     do ix1=ixOmin1+1,ixOmax1-1
       do ix3=ixOmin3+1,ixOmax3-1
         w(ix1,ix2,ix3,b2_)=w(ix1,ix2+2,ix3,b2_) &
          +delydelx*(w(ix1+1,ix2+1,ix3,b1_)-w(ix1-1,ix2+1,ix3,b1_))&
          +delydelz*(w(ix1,ix2+1,ix3+1,b3_)-w(ix1,ix2+1,ix3-1,b3_))
       enddo
     enddo
   enddo
   call conserve(ixG^L,ixO^L,w,x,patchfalse)
 case(4)
   w(ixO^S,rho_)=w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-dixB:-1,ixOmin3:ixOmax3,rho_)
   w(ixO^S,v1_)=-w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-dixB:-1,ixOmin3:ixOmax3,m1_)/&
                 w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-dixB:-1,ixOmin3:ixOmax3,rho_)
   w(ixO^S,v2_)=-w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-dixB:-1,ixOmin3:ixOmax3,m2_)/&
                 w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-dixB:-1,ixOmin3:ixOmax3,rho_)
   w(ixO^S,v3_)=-w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-dixB:-1,ixOmin3:ixOmax3,m3_)/&
                 w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-dixB:-1,ixOmin3:ixOmax3,rho_)
   do ix2=ixOmin2,ixOmax2
     w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,b1_:b3_)=(1.0d0/3.0d0)* &
                 (-w(ixOmin1:ixOmax1,ix2-2,ixOmin3:ixOmax3,b1_:b3_)&
            +4.0d0*w(ixOmin1:ixOmax1,ix2-1,ixOmin3:ixOmax3,b1_:b3_))
   enddo
   !! 2nd order CD for divB=0 to set normal B component better
   delydelx=dxlevel(2)/dxlevel(1)
   delydelz=dxlevel(2)/dxlevel(3)
   do ix2=ixOmin2,ixOmax2
     do ix1=ixOmin1+1,ixOmax1-1
       do ix3=ixOmin3+1,ixOmax3-1
         w(ix1,ix2,ix3,b2_)=w(ix1,ix2-2,ix3,b2_) &
          -delydelx*(w(ix1+1,ix2-1,ix3,b1_)-w(ix1-1,ix2-1,ix3,b1_))&
          -delydelz*(w(ix1,ix2-1,ix3+1,b3_)-w(ix1,ix2-1,ix3-1,b3_))
       enddo
     enddo
   enddo
   !endif
   call conserve(ixG^L,ixO^L,w,x,patchfalse)
 case(5)
   do ix3=ixOmax3,ixOmin3,-1
     w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3,b1_:b3_)=(1.0d0/11.0d0)* &
          ( +2.0d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3+3,b1_:b3_) &
            -9.0d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3+2,b1_:b3_) &
           +18.0d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3+1,b1_:b3_))
   enddo
   ! 2nd order CD for divB=0 to set normal B component better
   delzdelx=dxlevel(3)/dxlevel(1)
   delzdely=dxlevel(3)/dxlevel(2)
   do ix3=ixOmax3,ixOmin3,-1
     do ix1=ixOmin1+1,ixOmax1-1
       do ix2=ixOmin2+1,ixOmax2-1
         w(ix1,ix2,ix3,b3_)=w(ix1,ix2,ix3+2,b3_) &
         +delzdelx*(w(ix1+1,ix2,ix3+1,b1_)-w(ix1-1,ix2,ix3+1,b1_))&
         +delzdely*(w(ix1,ix2+1,ix3+1,b2_)-w(ix1,ix2-1,ix3+1,b2_))
       enddo
     enddo
   enddo
   coeffrho=eqpar(grav3_)*SRadius**2/Tiso
   w(ixO^S,rho_)=rhob*dexp(coeffrho*(1.d0/SRadius-1.d0/(x(ixO^S,3)+SRadius)))
   w(ixO^S,v1_)=-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmax3+dixB:ixOmax3+1:-1,m1_)/&
                 w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmax3+dixB:ixOmax3+1:-1,rho_)
   w(ixO^S,v2_)=-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmax3+dixB:ixOmax3+1:-1,m2_)/&
                 w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmax3+dixB:ixOmax3+1:-1,rho_)
   w(ixO^S,v3_)=-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmax3+dixB:ixOmax3+1:-1,m3_)/&
                 w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmax3+dixB:ixOmax3+1:-1,rho_)
   call conserve(ixG^L,ixO^L,w,x,patchfalse)
 case(6)
   ixIM^L=ixO^L;
   ixIMmin3=ixOmin3-1;ixIMmax3=ixOmax3;
   call getggrav(tmp,ixG^L,ixIM^L,x)
   invT(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3-1)/Tiso
   do ix3=ixOmin3,ixOmax3
     invT(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=invT(ixOmin1:ixOmax1,ixOmin2:ixOmax2)+&
     tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3)/Tiso
     w(ix3^%3ixO^S,rho_)=w(ixOmin3-1^%3ixO^S,rho_)*dexp(invT(ixOmin1:ixOmax1,ixOmin2:ixOmax2)*dxlevel(3)*0.5d0)
   enddo
   do ix3=ixOmin3,ixOmax3
     w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3,b1_:b3_)=(1.0d0/3.0d0)* &
           (     -w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3-2,b1_:b3_) &
            +4.d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3-1,b1_:b3_))
   enddo
   delzdelx=dxlevel(3)/dxlevel(1)
   delzdely=dxlevel(3)/dxlevel(2)
   do ix3=ixOmin3,ixOmax3
     do ix1=ixOmin1+1,ixOmax1-1
       do ix2=ixOmin2+1,ixOmax2-1
         w(ix1,ix2,ix3,b3_)=w(ix1,ix2,ix3-2,b3_) &
         -delzdelx*(w(ix1+1,ix2,ix3-1,b1_)-w(ix1-1,ix2,ix3-1,b1_))&
         -delzdely*(w(ix1,ix2+1,ix3-1,b2_)-w(ix1,ix2-1,ix3-1,b2_))
       enddo
     enddo
   enddo
   w(ixO^S,v1_)=-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3-1:ixOmin3-dixB:-1,m1_)/&
                 w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3-1:ixOmin3-dixB:-1,rho_)
   w(ixO^S,v2_)=-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3-1:ixOmin3-dixB:-1,m2_)/&
                 w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3-1:ixOmin3-dixB:-1,rho_)
   w(ixO^S,v3_)=-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3-1:ixOmin3-dixB:-1,m3_)/&
                 w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3-1:ixOmin3-dixB:-1,rho_)
   call conserve(ixG^L,ixO^L,w,x,patchfalse)
 case default
   call mpistop("Special boundary is not defined for this region")
end select

end subroutine specialbound_usr
!==============================================================================
subroutine getggrav(ggrid,ixI^L,ixO^L,x)

use mod_global_parameters

integer, intent(in) :: ixI^L,ixO^L
double precision, intent(in) :: x(ixI^S,1:ndim)
double precision, intent(out) :: ggrid(ixG^T)
!---------------------------------------------------------------------------
! calculate gravity
ggrid(ixO^S)=eqpar(grav3_)*(SRadius/(SRadius+x(ixO^S,3)))**2

end subroutine getggrav
!=============================================================================
subroutine addsource_gravSA(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

! w[iw]=w[iw]+qdt*S[wCT,qtC,x] where S is the source based on wCT within ixO
! gravity distribution along a magnetic loop (a circular arc)

use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L, iw^LIM
double precision, intent(in) :: qdt, qtC, qt, x(ixI^S,1:ndim)
double precision, intent(inout) :: wCT(ixI^S,1:nw), w(ixI^S,1:nw)

double precision :: ggrid(ixG^T)
integer :: iw, idims
!---------------------------------------------------------------------------

call getggrav(ggrid,ixI^L,ixO^L,x)

! add sources from gravity
do iw= iw^LIM
   select case (iw)
   case (m^D_)
     ! dm_i/dt= +rho*g_i
      idims=iw-m0_
      if (abs(eqpar(grav0_+idims))>smalldouble) &
          w(ixO^S,m0_+idims)=w(ixO^S,m0_+idims) &
              +qdt*ggrid(ixO^S)*wCT(ixO^S,rho_)
   end select
end do

end subroutine addsource_gravSA
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
subroutine specialsource(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L, iw^LIM
double precision, intent(in) :: qdt, qtC, qt, x(ixI^S,1:ndim)
double precision, intent(inout) :: w(ixI^S,1:nw), wCT(ixI^S,1:nw)
!-----------------------------------------------------------------------------
call addsource_gravSA(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

end subroutine specialsource
!=============================================================================
subroutine getdt_special(w,ixG^L,ix^L,dtnew,dx^D,x)

use mod_global_parameters

integer, intent(in) :: ixG^L,ix^L
double precision, intent(in) :: dx^D, x(ixG^S,1:ndim)
! note that depending on strictsmall etc, w values may change 
double precision, intent(inout) :: w(ixG^S,1:nw), dtnew
!-----------------------------------------------------------------------------
dtnew=bigdouble

call getdt_grav(w,ixG^L,ix^L,dtnew,dx^D,x)

end subroutine getdt_special
!=============================================================================
subroutine specialsource_impl(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L, iw^LIM
double precision, intent(in) :: qdt, qtC, qt, x(ixI^S,1:ndim)
double precision, intent(inout) :: w(ixI^S,1:nw), wCT(ixI^S,1:nw)
!-----------------------------------------------------------------------------

end subroutine specialsource_impl
!=============================================================================
subroutine getdt_impl(w,ixG^L,ix^L,dtnew,dx^D,x)

use mod_global_parameters

integer, intent(in) :: ixG^L,ix^L
double precision, intent(in) :: dx^D, x(ixG^S,1:ndim)
! note that depending on strictsmall etc, w values may change 
double precision, intent(inout) :: w(ixG^S,1:nw), dtnew
!-----------------------------------------------------------------------------
dtnew=bigdouble
end subroutine getdt_impl
!=============================================================================
subroutine specialeta(w,ixI^L,ix^L,idirmin,x,current,eta)

! Set the common "eta" array for resistive MHD based on w or the
! "current" variable which has components between idirmin and 3.

use mod_global_parameters

integer, intent(in) :: ixI^L, ix^L, idirmin
double precision, intent(in) :: w(ixI^S,nw), x(ixI^S,1:ndim)

double precision :: current(ixG^T,7-2*ndir:3), eta(ixG^T)
!!! common/resist/current,eta
!-----------------------------------------------------------------------------

end subroutine specialeta
!=============================================================================
subroutine specialrefine_grid(igrid,level,ixG^L,ix^L,qt,w,x,refine,coarsen)

! Enforce additional refinement or coarsening
! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.

use mod_global_parameters

integer, intent(in) :: igrid, level, ix^L, ixG^L
double precision, intent(in) :: qt, w(ixG^S,1:nw), x(ixG^S,1:ndim)
integer, intent(inout) :: refine, coarsen
!-----------------------------------------------------------------------------
!if (any(x(ix^S,3)<=xprobmin3+0.05d0)) then
!  refine=1
!  coarsen=-1
!endif

end subroutine specialrefine_grid
!=============================================================================
subroutine specialvarforerrest(ixI^L,ixO^L,iflag,w,var)

! this is the place to compute a local auxiliary variable to be used
! as refinement criterion for the Lohner error estimator only
!  -->it is then requiring and iflag>nw
! note that ixO=ixI=ixG, hence the term local (gradients need special attention!)

use mod_global_parameters

integer, intent(in)          :: ixI^L,ixO^L,iflag
double precision, intent(in) :: w(ixI^S,1:nw)
double precision, intent(out):: var(ixG^T)
!-----------------------------------------------------------------------------

if (iflag >nw)call mpistop(' iflag> nw, make change in parfile or in user file')

var(ixI^S) = zero

end subroutine specialvarforerrest

!=============================================================================
subroutine printlog_special

use mod_global_parameters
!-----------------------------------------------------------------------------

call printlog_default
call userspecialconvert(1)

end subroutine printlog_special
!=============================================================================
subroutine process_grid_usr(igrid,level,ixI^L,ixO^L,qt,w,x)

! this subroutine is ONLY to be used for computing auxiliary variables
! which happen to be non-local (like div v), and are in no way used for
! flux computations. As auxiliaries, they are also not advanced

use mod_global_parameters

integer, intent(in):: igrid,level,ixI^L,ixO^L
double precision, intent(in):: qt,x(ixI^S,1:ndim)
double precision, intent(inout):: w(ixI^S,1:nw)
!-----------------------------------------------------------------------------

end subroutine process_grid_usr
!=============================================================================
subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)

! this subroutine can be used in convert, to add auxiliary variables to the
! converted output file, for further analysis using tecplot, paraview, ....
! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
!
! the array normconv can be filled in the (nw+1:nw+nwauxio) range with
! corresponding normalization values (default value 1)

use mod_global_parameters

integer, intent(in)                :: ixI^L,ixO^L
double precision, intent(in)       :: x(ixI^S,1:ndim)
double precision                   :: w(ixI^S,nw+nwauxio)
double precision                   :: normconv(0:nw+nwauxio)

double precision                   :: qtC,qdt,dxinv(1:ndim)
double precision                   :: tmp(ixG^T),cmax(ixG^T),courantmax(ixG^T)
double precision                   :: qvec(ixG^T,1:ndir),current(ixG^T,1:ndir)
integer                            :: ix^D,idirmin,idims,idir,jdir,kdir
!-----------------------------------------------------------------------------
!! calculate current
if(B0field) then
 ^C&qvec(ixI^S,^C)=w(ixI^S,b^C_)+myB0_cell%w(ixI^S,^C);
else
 ^C&qvec(ixI^S,^C)=w(ixI^S,b^C_);
endif
w(ixO^S,nw+1)=dsqrt(^C&qvec(ixO^S,^C)**2+)
call curlvector(qvec,ixI^L,ixO^L,current,idirmin,1,ndir)
^C&w(ixO^S,nw+1+^C)=current(ixO^S,^C);
! output normalized divB
call getdivb(w,ixI^L,ixO^L,tmp)
w(ixO^S,nw+5)=tmp(ixO^S)
w(ixO^S,nw+6)=0.5d0*dabs(tmp(ixO^S))/w(ixO^S,nw+1)/(^D&1.0d0/dxlevel(^D)+)
! calculate Lorentz force
qvec(ixO^S,1:ndir)=zero
do idir=1,ndir; do jdir=1,ndir; do kdir=idirmin,3
   if(lvc(idir,jdir,kdir)/=0)then
      tmp(ixO^S)=current(ixO^S,jdir)*w(ixO^S,b0_+kdir)
      if(lvc(idir,jdir,kdir)==1)then
         qvec(ixO^S,idir)=qvec(ixO^S,idir)+tmp(ixO^S)
      else
         qvec(ixO^S,idir)=qvec(ixO^S,idir)-tmp(ixO^S)
      endif
   endif
enddo; enddo; enddo
^C&w(ixO^S,nw+6+^C)=qvec(ixO^S,^C);
w(ixO^S,nw+10)=dasin(dsqrt((^C&qvec(ixO^S,^C)**2+)/((^C&current(ixO^S,^C)**2+)+ smalldouble))/w(ixO^S,nw+1))*180.d0/dpi

end subroutine specialvar_output
!=============================================================================
subroutine specialvarnames_output

! newly added variables need to be concatenated with the wnames/primnames string

use mod_global_parameters
!-----------------------------------------------------------------------------
oktest = index(teststr,'printlog')>=1

primnames= TRIM(primnames)//' '//'B j1 j2 j3 divb f_i L1 L2 L3 theta'
wnames=TRIM(wnames)//' '//'B j1 j2 j3 divb f_i L1 L2 L3 theta'

end subroutine specialvarnames_output
!=============================================================================
subroutine specialset_B0(ixI^L,ixO^L,x,wB0)

! Here one can add a steady (time-independent) potential background field

use mod_global_parameters

integer, intent(in)           :: ixI^L,ixO^L
double precision, intent(in)  :: x(ixG^T,1:ndim)
double precision, intent(inout) :: wB0(ixI^S,1:ndir)

double precision :: slopeL, slopeR, slopeC, signC, signR
double precision :: slope(ndir,ndim),xlen^D,res^D
integer :: ix^D,idim,ixCo^D, jxCo^D, hxCo^D,iw
!-----------------------------------------------------------------------------
end subroutine specialset_B0
!=============================================================================
subroutine bc_int(qt,ixG^L,ixO^L,w,x)

! internal boundary, user defined
!
! This subroutine can be used to artificially overwrite ALL conservative 
! variables in a user-selected region of the mesh, and thereby act as
! an internal boundary region. It is called just before external (ghost cell)
! boundary regions will be set by the BC selection. Here, you could e.g. 
! want to introduce an extra variable (nwextra, to be distinguished from nwaux)
! which can be used to identify the internal boundary region location.
! Its effect should always be local as it acts on the mesh.
!

use mod_global_parameters

integer, intent(in) :: ixG^L,ixO^L
double precision, intent(in) :: qt
double precision, intent(inout) :: w(ixG^S,1:nw)
double precision, intent(in) :: x(ixG^S,1:ndim)

! .. local ..
!logical :: patchw(ixG^T)
!----------------------------------------------------------------------------

call mpistop("bc_int not defined")
end subroutine bc_int
!=============================================================================
subroutine userspecialconvert(qunitconvert)

use mod_global_parameters
integer, intent(in) :: qunitconvert
character(len=20):: userconvert_type

integer  :: iigrid,igrid
logical  :: patchw(ixG^T)
!-----------------------------------------------------------------------------

call spatial_integral_w
end subroutine userspecialconvert
!=============================================================================
subroutine spatial_integral_w

use mod_global_parameters

double precision :: dvolume(ixG^T), dsurface(ixG^T),timephy,dvone
double precision, allocatable :: integral_ipe(:), integral_w(:)
double precision, external :: integral_grid

integer           :: nregions,ireg,ncellpe,ncell,idims,hxM^LL,nx^D
integer           :: iigrid,igrid,status(MPI_STATUS_SIZE),ni
character(len=100):: filename,region
character(len=1024) :: line, datastr
logical           :: patchwi(ixG^T),alive
!-----------------------------------------------------------------------------

nregions=1
! number of integrals to perform
ni=3
allocate(integral_ipe(ni),integral_w(ni))
integral_ipe=0.d0
integral_w=0.d0
nx^D=ixMhi^D-ixMlo^D+1;
do ireg=1,nregions
 select case(ireg)
 case(1)
   region='fulldomain'
 case(2)
   region='cropped'
 end select
 ncellpe=0 
 do iigrid=1,igridstail; igrid=igrids(iigrid);
   if(slab) then
     dvone={rnode(rpdx^D_,igrid)|*}
     dvolume(ixM^T)=dvone
     dsurface(ixM^T)=two*(^D&dvone/rnode(rpdx^D_,igrid)+)
   else
     dvolume(ixM^T)=pgeo(igrid)%dvolume(ixM^T)
     dsurface(ixM^T)= ^D&pgeo(igrid)%surfaceC^D(ixM^T)+
     do idims=1,ndim
       hxM^LL=ixM^LL-kr(idims,^D);
       select case(idims)
       {case(^D)
          dsurface(ixM^T)=dsurface(ixM^T)+pgeo(igrid)%surfaceC^D(hxM^T) \}
       end select
     end do
   end if
   ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
   if (.not.slab) mygeo => pgeo(igrid)
   if (B0field) then
      myB0_cell => pB0_cell(igrid)
     {^D&myB0_face^D => pB0_face^D(igrid)\}
   end if
   typelimiter=typelimiter1(node(plevel_,igrid))
   typegradlimiter=typegradlimiter1(node(plevel_,igrid))
   patchwi(ixG^T)=.false.
   select case(region)
   case('fulldomain')
      patchwi(ixM^T)=.true.
      ncellpe=ncellpe+{nx^D*}
   case('cropped')
      call mask_grid(ixG^LL,ixM^LL,pw(igrid)%w,px(igrid)%x,patchwi,ncellpe)
   case default
      call mpistop("region not defined")
   end select
   integral_ipe(1)=integral_ipe(1)+ &
             integral_grid(ixG^LL,ixM^LL,pw(igrid)%w,px(igrid)%x,dvolume,dsurface,1,patchwi)
   integral_ipe(2)=integral_ipe(2)+ &
             integral_grid(ixG^LL,ixM^LL,pw(igrid)%w,px(igrid)%x,dvolume,dsurface,2,patchwi)
   integral_ipe(3)=integral_ipe(3)+ &
             integral_grid(ixG^LL,ixM^LL,pw(igrid)%w,px(igrid)%x,dvolume,dsurface,3,patchwi)
 end do
 call MPI_ALLREDUCE(integral_ipe,integral_w,ni,MPI_DOUBLE_PRECISION,&
                      MPI_SUM,icomm,ierrmpi)
 call MPI_ALLREDUCE(ncellpe,ncell,1,MPI_INTEGER,MPI_SUM,icomm,ierrmpi)
 integral_w(3)=integral_w(3)/dble(ncell)
 timephy=t
 if(mype==0) then
   write(filename,"(a,a,a)") TRIM(filenamelog),TRIM(region),".csv"
   inquire(file=filename,exist=alive)
   if(alive) then
     open(unit=21,file=filename,form='formatted',status='old',access='append')
   else
     open(unit=21,file=filename,form='formatted',status='new')
     write(21,'(a)') 'time, f_i, CWsin, theta'
   endif
   write(datastr,'(es11.4, a)') timephy,','
   line=datastr
   write(datastr,"(es11.4, a)") integral_w(3),','
   line = trim(line)//trim(datastr)
   write(datastr,"(es11.4, a)") integral_w(1)/(integral_w(2)+smalldouble),','
   line = trim(line)//trim(datastr)
   write(datastr,"(es11.4)") dasin(integral_w(1)/(integral_w(2)+smalldouble))*180.d0/dpi
   line = trim(line)//trim(datastr)
   write(21,'(a)') trim(line)
   close(21)
 endif
enddo
deallocate(integral_ipe,integral_w)

end subroutine spatial_integral_w
!=============================================================================
subroutine mask_grid(ixI^L,ixO^L,w,x,patchwi,cellcount)

use mod_global_parameters

integer, intent(in)                :: ixI^L,ixO^L
double precision, intent(in)       :: x(ixI^S,1:ndim)
double precision                   :: w(ixI^S,nw+nwauxio)
logical, intent(inout)             :: patchwi(ixG^T)

double precision  ::  buff
integer                            :: ix^D,cellcount
!-----------------------------------------------------------------------------
buff=0.05d0*(xprobmax1-xprobmin1)
{do ix^DB=ixOmin^DB,ixOmax^DB\}
   if(x(ix^D,1)>xprobmin1+buff .and. x(ix^D,1)<xprobmax1-buff .and. &
      x(ix^D,2)>xprobmin2+buff .and. x(ix^D,2)<xprobmax2-buff) then
     patchwi(ix^D)=.true.
     cellcount=cellcount+1
   else
     patchwi(ix^D)=.false.
   endif
{end do\}
return
end subroutine mask_grid
!=============================================================================
function integral_grid(ixI^L,ixO^L,w,x,dvolume,dsurface,intval,patchwi)

use mod_global_parameters

integer, intent(in)                :: ixI^L,ixO^L,intval
double precision, intent(in)       :: x(ixI^S,1:ndim),dvolume(ixG^T),dsurface(ixG^T)
double precision, intent(in)       :: w(ixI^S,nw)
logical, intent(in) :: patchwi(ixG^T)

double precision, dimension(ixG^T,1:ndir) :: bvec,qvec
double precision :: current(ixG^T,7-2*ndir:3),tmp(ixG^T)
double precision :: integral_grid,mcurrent
integer :: ix^D,idirmin,idir,jdir,kdir
!-----------------------------------------------------------------------------

integral_grid=0.d0
select case(intval)
 case(1)
  ! current times sin theta (between J and B)
  if(B0field) then
    ^C&bvec(ixI^S,^C)=w(ixI^S,b^C_)+myB0_cell%w(ixI^S,^C);
  else
    ^C&bvec(ixI^S,^C)=w(ixI^S,b^C_);
  endif
  call getcurrent(w,ixI^L,ixO^L,idirmin,current)
  ! calculate Lorentz force
  qvec(ixO^S,1:ndir)=zero
  do idir=1,ndir; do jdir=idirmin,3; do kdir=1,ndir
     if(lvc(idir,jdir,kdir)/=0)then
        tmp(ixO^S)=current(ixO^S,jdir)*bvec(ixO^S,kdir)
        if(lvc(idir,jdir,kdir)==1)then
           qvec(ixO^S,idir)=qvec(ixO^S,idir)+tmp(ixO^S)
        else
           qvec(ixO^S,idir)=qvec(ixO^S,idir)-tmp(ixO^S)
        endif
     endif
  enddo; enddo; enddo
  {do ix^DB=ixOmin^DB,ixOmax^DB\}
     if(patchwi(ix^D)) integral_grid=integral_grid+dsqrt((^C&qvec(ix^D,^C)**2+)/&
      (^C&bvec(ix^D,^C)**2+))*dvolume(ix^D)
  {end do\}
 case(2)
  call getcurrent(w,ixI^L,ixO^L,idirmin,current)
  ! integral of current
  {do ix^DB=ixOmin^DB,ixOmax^DB\}
     mcurrent=0.d0
     do idir=7-2*ndir,3
        mcurrent=current(ix^D,idir)**2+mcurrent
     end do
     if(patchwi(ix^D))  integral_grid=integral_grid+dsqrt(mcurrent)*dvolume(ix^D)
  {end do\}
 case(3)
  ! f_i solenoidal property of B: (dvolume |div B|)/(dsurface |B|)
  if(B0field) then
    ^C&bvec(ixI^S,^C)=w(ixI^S,b^C_)+myB0_cell%w(ixI^S,^C);
  else
    ^C&bvec(ixI^S,^C)=w(ixI^S,b^C_);
  endif
  call divvector(bvec,ixI^L,ixO^L,tmp)
  {do ix^DB=ixOmin^DB,ixOmax^DB\}
     if(patchwi(ix^D)) integral_grid=integral_grid+dabs(tmp(ix^D))*&
           dvolume(ix^D)/dsqrt(^C&bvec(ix^D,^C)**2+)/dsurface(ix^D)
  {end do\}
 case default
     call mpistop("intval not defined")
end select

return
end function integral_grid
!=============================================================================
! amrvacusr.t.arfff
!=============================================================================
