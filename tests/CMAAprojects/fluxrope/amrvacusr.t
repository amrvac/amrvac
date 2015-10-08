!=============================================================================
! amrvacusr.t.fluxrope23
!=============================================================================
INCLUDE:amrvacmodules/gravity.t
INCLUDE:amrvacnul/usrflags.t
!=============================================================================
module usr_bc
implicit none
double precision, allocatable, save :: pbc(:),rbc(:)
double precision, allocatable, save :: Bbc(:,:,:),Bbt(:,:,:)
double precision, allocatable, save :: vbc(:,:,:),vbt(:,:,:)
integer, save :: nxbc^D,ngl
end module usr_bc
!=============================================================================
subroutine initglobaldata_usr

include 'amrvacdef.f'

logical, save :: firstusrglobaldata=.true.
!-----------------------------------------------------------------------------
! CGS Unit
k_B=1.3806d-16     ! erg*K^-1
miu0=4.d0*dpi      ! Gauss^2 cm^2 dyne^-1
Lunit=1.d9         ! cm
Teunit=1.d6        ! K
nHunit=1.d9        ! cm^-3
mHunit=1.67262d-24 ! g
runit=1.4d0*mHunit*nHunit     ! 2.341668000000000E-015 g*cm^-3
punit=2.3d0*nHunit*k_B*Teunit ! 0.317538000000000 erg*cm^-3
Bunit=dsqrt(miu0*punit)       ! 1.99757357615242 Gauss
vunit=Bunit/dsqrt(miu0*runit) ! 1.16448846777562E007 cm/s = 116.45 km/s
tunit=Lunit/vunit             ! 85.8746159942810 s
heatunit=punit/tunit          ! 3.697693390805347E-003 erg*cm^-3/s
!! The third component of magnetic field Bz is constant 
eqpar(LB0_)=xprobmax1-xprobmin1

eqpar(grav1_)=0.d0
eqpar(grav2_)=-2.74d4*Lunit/vunit**2
!eqpar(eta_)=-1.d0
!! base value of T, rho 
Tpho=1.d4/Teunit
Ttop=2.0d0
!rpho=2.d13/nHunit 
!rpho=2.552d14/nHunit 
!rpho=1.151d15/nHunit
!rpho=1.d0 !case relaxb
rpho=0.7d0
!rpho=8.95d14/nHunit 
!! spacial resolution of hdstatic solution
gzone=0.2d0
dya=(2.d0*gzone+xprobmax2-xprobmin2)/dble(jmax)
!!parameters for force free field
B0=3.d0/Bunit
! theta is the angle to the plane xy, 90-theta is the angle to the channel
theta=60.d0*dpi/180.d0
! time stop driving
tdstop=1500.d0/tunit

if(firstusrglobaldata) then
!! setup global gravity stratification in one vertical line
  call inithdstatic
  firstusrglobaldata=.false.
endif

end subroutine initglobaldata_usr
!=============================================================================
subroutine inithdstatic
!! initialize the table in a vertical line through the global domain
use usr_bc
include 'amrvacdef.f'

double precision :: Ttr,Fc,kx,ly,v0,post
double precision, allocatable :: xbc(:,:,:)
integer :: j,na,nb,ibc,ix^D
!----------------------------------------------------------------------------

ngl=dixB*2**(mxnest-1)
nxbc1=nint((xprobmax1-xprobmin1)/dx(1,mxnest))+2*ngl
nxbc2=ngl
allocate(Bbc(nxbc1,nxbc2+1,3))
allocate(Bbt(nxbc1,nxbc2+1,3))
allocate(vbc(nxbc1,nxbc2+1,3))
allocate(vbt(nxbc1,nxbc2+1,3))
allocate(xbc(nxbc1,nxbc2+1,2))
kx=dpi/eqpar(LB0_)
ly=kx*dcos(theta)
do ix2=1,nxbc2+1
  do ix1=1,nxbc1
    xbc(ix^D,1)=xprobmin1+(dble(ix1-ngl)-0.5d0)*dx(1,mxnest)
    xbc(ix^D,2)=xprobmin2+(dble(ix2-ngl)-0.5d0)*dx(2,mxnest)
    Bbc(ix^D,1)=-ly/kx*B0*dcos(kx*xbc(ix^D,1))*dexp(-ly*xbc(ix^D,2))
    Bbc(ix^D,2)=B0*dsin(kx*xbc(ix^D,1))*dexp(-ly*xbc(ix^D,2))
    Bbc(ix^D,3)=-dsin(theta)*B0*dcos(kx*xbc(ix^D,1))*dexp(-ly*xbc(ix^D,2))
  end do
end do    
Bbt=Bbc
vbc=0.d0
vbt=0.d0
!v0=-1.2d6/vunit
v0=-7.d5/vunit
post=eqpar(LB0_)/8.d0
where(dabs(xbc(1:nxbc1,1:nxbc2+1,1))<post)
  vbc(1:nxbc1,1:nxbc2+1,1)=v0*xbc(1:nxbc1,1:nxbc2+1,1)/post
elsewhere(xbc(1:nxbc1,1:nxbc2+1,1)>=post)
  vbc(1:nxbc1,1:nxbc2+1,1)=v0*(eqpar(LB0_)*0.5d0-xbc(1:nxbc1,1:nxbc2+1,1))/post/3.d0
elsewhere
  vbc(1:nxbc1,1:nxbc2+1,1)=-v0*(eqpar(LB0_)*0.5d0+xbc(1:nxbc1,1:nxbc2+1,1))/post/3.d0
end where

deallocate(xbc)
end subroutine inithdstatic
!=============================================================================
subroutine initonegrid_usr(ixG^L,ix^L,w,x)

! initialize one grid

include 'amrvacdef.f'

integer, intent(in) :: ixG^L,ix^L
double precision, intent(in) :: x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)

double precision :: kx,ly,Tiso

logical :: patchw(ixG^T)
logical, save :: first=.true.
!----------------------------------------------------------------------------

if (ndir/=3) call mpistop("this is a mhd23 problem!")
if(first)then
  if(mype==0) then
    write(*,*)'Simulating 2.5D flux rope formation in corona'
  endif
  first=.false.
endif
Tiso=1.d0
!w(ixG^S,rho_)=dexp(eqpar(grav2_)*SRadius**2/Tiso*&
!               (1.d0/SRadius-1.d0/(x(ixG^S,2)+SRadius)))
w(ixG^S,rho_)=rpho*dexp(eqpar(grav2_)*x(ixG^S,2)/Tiso)
w(ixG^S,p_)=w(ixG^S,rho_)*1.d0
w(ixG^S,v1_)=zero
w(ixG^S,v2_)=zero
w(ixG^S,v3_)=zero

kx=dpi/eqpar(LB0_)
ly=kx*dcos(theta)
w(ixG^S,b1_)=-ly/kx*B0*dcos(kx*x(ixG^S,1))*dexp(-ly*x(ixG^S,2))
w(ixG^S,b2_)=B0*dsin(kx*x(ixG^S,1))*dexp(-ly*x(ixG^S,2))
w(ixG^S,b3_)=-dsin(theta)*B0*dcos(kx*x(ixG^S,1))*dexp(-ly*x(ixG^S,2))
patchw(ixG^S)=.false.
call conserve(ixG^L,ix^L,w,x,patchw)

{#IFDEF GLM
w(ixG^S,psi_)=0.d0
}

end subroutine initonegrid_usr
!=============================================================================
subroutine specialbound_usr(qt,ixG^L,ixO^L,iw,iB,w,x)

! special boundary types, user defined

use usr_bc
include 'amrvacdef.f'

integer, intent(in) :: ixO^L, iw, iB, ixG^L
double precision, intent(in) :: qt, x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)

double precision :: pth(ixG^T),ggrid(ixG^T)
double precision :: dxa^D,xlen^D,startpos^D,coeffrho,delydelx,cg,kx,ly
integer :: ix^D,idims,ixInt^L,af,ixbc^D
!----------------------------------------------------------------------------
if(energyonly) return
select case(iB)
case(3)
  dxa^D=dx(^D,mxnest);
  af=nint(dxlevel(2)/dxa2)
  startpos^D=xprobmin^D-dble(ngl)*dxa^D;\
  {do ix^DB=ixOmin^DB,ixOmax^DB\}
     xlen^D=x({ix^D,},^D)-startpos^D;\
     if(af>=2) then
       af=af/2
       ixbc^D=nint(xlen^D/dxa^D);\
       w(ix^D,v1_)=sum(vbt(ixbc^D-af+1:ixbc^D+af,1))/dble(size(vbt(ixbc^D-af+&
         1:ixbc^D+af,1)))
       w(ix^D,v2_)=sum(vbt(ixbc^D-af+1:ixbc^D+af,2))/dble(size(vbt(ixbc^D-af+&
         1:ixbc^D+af,2)))
       w(ix^D,v3_)=sum(vbt(ixbc^D-af+1:ixbc^D+af,3))/dble(size(vbt(ixbc^D-af+&
         1:ixbc^D+af,3)))
       w(ix^D,b1_)=sum(Bbt(ixbc^D-af+1:ixbc^D+af,1))/dble(size(Bbt(ixbc^D-af+&
         1:ixbc^D+af,1)))
       w(ix^D,b2_)=sum(Bbt(ixbc^D-af+1:ixbc^D+af,2))/dble(size(Bbt(ixbc^D-af+&
         1:ixbc^D+af,2)))
       w(ix^D,b3_)=sum(Bbt(ixbc^D-af+1:ixbc^D+af,3))/dble(size(Bbt(ixbc^D-af+&
         1:ixbc^D+af,3)))
     else
       ixbc^D=ceiling(xlen^D/dxa^D);\
       w(ix^D,v1_:v3_)=vbt(ixbc^D,1:3)
       w(ix^D,b1_:b3_)=Bbt(ixbc^D,1:3)
     endif
  {end do\}
  !w(ixO^S,v2_) =-w(ixOmin1:ixOmax1,ixOmax2+dixB:ixOmax2+1:-1,m2_)&
  !             /w(ixOmin1:ixOmax1,ixOmax2+dixB:ixOmax2+1:-1,rho_)
  !! 2nd order CD for divB=0 to set normal B component better
  !if(all(dabs(x(ixO^S,1))<0.975d0)) then
  !  delydelx=dxlevel(2)/dxlevel(1)
  !  do ix2=ixOmax2,ixOmin2,-1
  !    do ix1=ixOmin1+1,ixOmax1-1
  !      if(dabs(x(ix1,ix2,1))>0.9975d0) cycle
  !      w(ix1,ix2,b2_)=w(ix1,ix2+2,b2_)+delydelx*(w(ix1+1,ix2+1,b1_)-&
  !        w(ix1-1,ix2+1,b1_))
  !    enddo
  !  enddo
  !endif
  if(qt>=tdstop) then
  ! fixed zero velocity
    w(ixO^S,v1_) =-w(ixOmin1:ixOmax1,ixOmax2+dixB:ixOmax2+1:-1,m1_)&
                  /w(ixOmin1:ixOmax1,ixOmax2+dixB:ixOmax2+1:-1,rho_)
    w(ixO^S,v2_) =-w(ixOmin1:ixOmax1,ixOmax2+dixB:ixOmax2+1:-1,m2_)&
                  /w(ixOmin1:ixOmax1,ixOmax2+dixB:ixOmax2+1:-1,rho_)
    w(ixO^S,v3_) =-w(ixOmin1:ixOmax1,ixOmax2+dixB:ixOmax2+1:-1,m3_)&
                  /w(ixOmin1:ixOmax1,ixOmax2+dixB:ixOmax2+1:-1,rho_)
  endif
  !w(ixO^S,rho_)=dexp(eqpar(grav2_)*SRadius**2*&
  !               (1.d0/SRadius-1.d0/(x(ixO^S,2)+SRadius)))
  w(ixO^S,rho_)=rpho*dexp(eqpar(grav2_)*x(ixO^S,2))
  w(ixO^S,p_)=w(ixO^S,rho_)
  !! fill pth, rho ghost layers according to gravity stratification
  !ixInt^L=ixO^L;
  !ixIntmin2=ixOmax2+1;ixIntmax2=ixOmax2+2;
  !call getpthermal(w,x,ixG^L,ixInt^L,pth)
  !dxa2=2.d0*dxlevel(2)*eqpar(grav2_)
  !do ix2=ixOmax2,ixOmin2,-1
  !  do ix1=ixOmin1,ixOmax1
  !    pth(ix1,ix2)=-dxa2*w(ix1,ix2+1,rho_)+pth(ix1,ix2+2)
  !    w(ix1,ix2,rho_)=pth(ix1,ix2)/pth(ix1,ix2+1)*w(ix1,ix2+1,rho_)
  !  end do
  !end do
  !w(ixO^S,p_)=pth(ixO^S)
  call conserve(ixG^L,ixO^L,w,x,patchfalse)
case(4)
  !! obtain the thermal pressure in top layers and store in pth
  !ixInt^L=ixO^L;
  !ixIntmin2=ixOmin2-2;ixIntmax2=ixOmin2-1;
  !call getpthermal(w,x,ixG^L,ixInt^L,pth)
  !!! fill pth, rho ghost layers according to gravity stratification
  !dxa2=2.d0*dxlevel(2)*eqpar(grav2_)
  !do ix2=ixOmin2,ixOmax2
  !  do ix1=ixOmin1,ixOmax1
  !    pth(ix1,ix2)=dxa2*w(ix1,ix2-1,rho_)+pth(ix1,ix2-2)
  !    w(ix1,ix2,rho_)=pth(ix1,ix2)/pth(ix1,ix2-1)*w(ix1,ix2-1,rho_)
  !  end do
  !end do
  !w(ixO^S,p_)=pth(ixO^S)
  w(ixO^S,rho_)=rpho*dexp(eqpar(grav2_)*x(ixO^S,2))
  w(ixO^S,p_)=w(ixO^S,rho_)
  !do ix2=ixOmin2,ixOmax2
    !w(ixOmin1:ixOmax1,ix2,v1_:v3_)=(1.0d0/25.0d0)* &
    !     ( -3.0d0*w(ixOmin1:ixOmax1,ix2-4,v1_:v3_) &
    !      +16.0d0*w(ixOmin1:ixOmax1,ix2-3,v1_:v3_) &
    !      -36.0d0*w(ixOmin1:ixOmax1,ix2-2,v1_:v3_) &
    !      +48.0d0*w(ixOmin1:ixOmax1,ix2-1,v1_:v3_))
  !  w(ixOmin1:ixOmax1,ix2,v1_:v3_)=w(ixOmin1:ixOmax1,ixOmin2-1,v1_:v3_)
  !enddo
  ! fixed zero velocity
  w(ixO^S,v1_) =-w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-dixB:-1,m1_)&
               /w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-dixB:-1,rho_)
  w(ixO^S,v2_) =-w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-dixB:-1,m2_)&
               /w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-dixB:-1,rho_)
  w(ixO^S,v3_) =-w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-dixB:-1,m3_)&
               /w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-dixB:-1,rho_)
  ! fixed b1 b2 b3
  kx=dpi/eqpar(LB0_)
  ly=kx*dcos(theta)
  w(ixO^S,b1_)=-ly/kx*B0*dcos(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))
  w(ixO^S,b2_)=B0*dsin(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))
  w(ixO^S,b3_)=-dsin(theta)*B0*dcos(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))
  !do ix2=ixOmin2,ixOmax2
  !  w(ixOmin1:ixOmax1,ix2,b1_:b3_)=(1.0d0/3.0d0)* &
  !              (-w(ixOmin1:ixOmax1,ix2-2,b1_:b3_)&
  !         +4.0d0*w(ixOmin1:ixOmax1,ix2-1,b1_:b3_))
  !enddo
  !! 2nd order CD for divB=0 to set normal B component better
  !delydelx=dxlevel(2)/dxlevel(1)
  !do ix2=ixOmin2,ixOmax2
  !  do ix1=ixOmin1+1,ixOmax1-1
  !    w(ix1,ix2,b2_)=w(ix1,ix2-2,b2_)-delydelx*(w(ix1+1,ix2-1,b1_)-&
  !      w(ix1-1,ix2-1,b1_))
  !  enddo
  !enddo
  call conserve(ixG^L,ixO^L,w,x,patchfalse)
case default
   call mpistop("Special boundary is not defined for this region")
end select

end subroutine specialbound_usr
!=============================================================================
subroutine specialsource(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

include 'amrvacdef.f'

integer, intent(in) :: ixI^L, ixO^L, iw^LIM
double precision, intent(in) :: qdt, qtC, qt, wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
double precision, intent(inout) :: w(ixI^S,1:nw)

! coordinate info can be used in region ixO

integer :: iw
!-----------------------------------------------------------------------------
call addsource_grav(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

end subroutine specialsource
!=============================================================================
subroutine getdt_special(w,ixG^L,ix^L,dtnew,dx^D,x)

! Limit "dt" further if necessary, e.g. due to the special source terms.
! The getdt_courant (CFL condition) and the getdt subroutine in the AMRVACPHYS
! module have already been called.

include 'amrvacdef.f'

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

include 'amrvacdef.f'

integer, intent(in) :: ixI^L, ixO^L, iw^LIM
double precision, intent(in) :: qdt, qtC, qt, x(ixI^S,1:ndim)
double precision, intent(inout) :: w(ixI^S,1:nw), wCT(ixI^S,1:nw)
!-----------------------------------------------------------------------------

end subroutine specialsource_impl
!=============================================================================
subroutine getdt_impl(w,ixG^L,ix^L,dtnew,dx^D,x)

include 'amrvacdef.f'

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

include 'amrvacdef.f'

integer, intent(in) :: ixI^L, ix^L, idirmin
double precision, intent(in) :: w(ixI^S,nw), x(ixI^S,1:ndim)

double precision :: current(ixG^T,7-2*ndir:3), eta(ixG^T)
double precision :: Jc,mu0,tmp(ixG^T),etahat
!!! common/resist/current,eta
!-----------------------------------------------------------------------------
Jc=1.2d-8/Bunit*Lunit ! Ampere cm^-2
mu0=3.6d13*tunit/Lunit**2
tmp(ix^S)=dsqrt(current(ix^S,1)**2+current(ix^S,2)**2+current(ix^S,3)**2)
where(tmp(ix^S)<Jc)
  eta(ix^S)=0.d0
elsewhere
  eta(ix^S)=mu0*(tmp(ix^S)/Jc-1.d0)**2
endwhere
etahat=1.8d14*tunit/Lunit**2
where(eta(ix^S)>etahat)
  eta(ix^S)=etahat
endwhere

end subroutine specialeta
!=============================================================================
subroutine specialrefine_grid(igrid,level,ixG^L,ix^L,qt,w,x,refine,coarsen)

! Enforce additional refinement or coarsening
! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.

include 'amrvacdef.f'

integer, intent(in) :: igrid, level, ix^L, ixG^L
double precision, intent(in) :: qt, w(ixG^S,1:nw), x(ixG^S,1:ndim)
integer, intent(inout) :: refine, coarsen
!-----------------------------------------------------------------------------

if (any(x(ix^S,2)<=xprobmin2+0.05d0)) then
  refine=1
  coarsen=-1
endif

end subroutine specialrefine_grid
!=============================================================================
subroutine specialvarforerrest(ixI^L,ixO^L,iflag,w,var)

! this is the place to compute a local auxiliary variable to be used
! as refinement criterion for the Lohner error estimator only
!  -->it is then requiring and iflag>nw
! note that ixO=ixI=ixG, hence the term local (gradients need special attention!)

include 'amrvacdef.f'

integer, intent(in)          :: ixI^L,ixO^L,iflag
double precision, intent(in) :: w(ixI^S,1:nw)
double precision, intent(out):: var(ixG^T)
!-----------------------------------------------------------------------------

if (iflag >nw)call mpistop(' iflag> nw, make change in parfile or in user file')

var(ixI^S) = zero

end subroutine specialvarforerrest
!=============================================================================
subroutine printlog_special

include 'amrvacdef.f'
!-----------------------------------------------------------------------------
oktest = index(teststr,'printlog')>=1

call mpistop("special log file undefined")

end subroutine printlog_special
!=============================================================================
subroutine process_grid_usr(igrid,level,ixI^L,ixO^L,qt,w,x)

! this subroutine is ONLY to be used for computing auxiliary variables
! which happen to be non-local (like div v), and are in no way used for
! flux computations. As auxiliaries, they are also not advanced

include 'amrvacdef.f'

integer, intent(in):: igrid,level,ixI^L,ixO^L
double precision, intent(in):: qt,x(ixI^S,1:ndim)
double precision, intent(inout):: w(ixI^S,1:nw)

double precision :: tmp(ixG^T),tmp1(ixG^T),tmp2(ixG^T)
double precision :: qvec(ixG^T,1:ndir),curlvec(ixG^T,1:ndir)
integer :: idirmin
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

include 'amrvacdef.f'

integer, intent(in)                :: ixI^L,ixO^L
double precision, intent(in)       :: x(ixI^S,1:ndim)
double precision                   :: w(ixI^S,nw+nwauxio)
double precision                   :: normconv(0:nw+nwauxio)

double precision :: tmp(ixG^T),divb(ixG^T),gradp(ixG^T)
double precision :: qvec(ixG^T,1:ndir),curlvec(ixG^T,1:ndir),kwave
integer                            :: idirmin,idims,ix^D,ixB^L,jxC^L
logical          :: patchw(ixG^T)
!-----------------------------------------------------------------------------
! output Te
if(saveprim)then
   tmp(ixO^S)=w(ixO^S,p_)
 else
   call getpthermal(w,x,ixI^L,ixO^L,tmp)
endif
w(ixO^S,nw+1)=tmp(ixO^S)/w(ixO^S,rho_)
!! store current
if(B0field) then
 ^C&qvec(ixI^S,^C)=w(ixI^S,b^C_)+myB0_cell%w(ixI^S,^C);
else
 ^C&qvec(ixI^S,^C)=w(ixI^S,b^C_);
endif
call curlvector(qvec,ixI^L,ixO^L,curlvec,idirmin,1,ndir)
^C&w(ixO^S,nw+1+^C)=curlvec(ixO^S,^C);
! output divB1
call getdivb(w,ixI^L,ixO^L,divb)
w(ixO^S,nw+5)=divb(ixO^S)
! L1
w(ixO^S,nw+6)=curlvec(ixO^S,2)*qvec(ixO^S,3)-curlvec(ixO^S,3)*qvec(ixO^S,2)
w(ixO^S,nw+7)=curlvec(ixO^S,3)*qvec(ixO^S,1)-curlvec(ixO^S,1)*qvec(ixO^S,3)
w(ixO^S,nw+8)=curlvec(ixO^S,1)*qvec(ixO^S,2)-curlvec(ixO^S,2)*qvec(ixO^S,1)
call gradient(tmp,ixI^L,ixO^L,2,gradp)
w(ixO^S,nw+9)=gradp(ixO^S)

end subroutine specialvar_output
!=============================================================================
subroutine specialvarnames_output

! newly added variables need to be concatenated with the wnames/primnames string

include 'amrvacdef.f'
!-----------------------------------------------------------------------------
oktest = index(teststr,'printlog')>=1

primnames=TRIM(primnames)//' '//'Te j1 j2 j3 divb L1 L2 L3 gradp'
   wnames=   TRIM(wnames)//' '//'Te j1 j2 j3 divb L1 L2 L3 gradp'

end subroutine specialvarnames_output
!=============================================================================
subroutine specialset_B0(ixI^L,ixO^L,x,wB0)

! Here one can add a steady (time-independent) potential background field

include 'amrvacdef.f'

integer, intent(in)           :: ixI^L,ixO^L
double precision, intent(in)  :: x(ixG^T,1:ndim)
double precision, intent(inout) :: wB0(ixI^S,1:ndir)
!-----------------------------------------------------------------------------

end subroutine specialset_B0
!=============================================================================
subroutine bc_int(level,qt,ixG^L,ixO^L,w,x)

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

include 'amrvacdef.f'

integer, intent(in) :: ixG^L,ixO^L,level
double precision, intent(in) :: qt
double precision, intent(inout) :: w(ixG^S,1:nw)
double precision, intent(in) :: x(ixG^S,1:ndim)

! .. local ..
!logical :: patchw(ixG^T)
!----------------------------------------------------------------------------

call mpistop("bc_int not defined")

end subroutine bc_int
!=============================================================================
subroutine mask_gridfluxrope(ixI^L,ixO^L,w,x,patchwi)

include 'amrvacdef.f'

integer, intent(in)                :: ixI^L,ixO^L
double precision, intent(in)       :: x(ixI^S,1:ndim)
double precision                   :: w(ixI^S,nw+nwauxio)
logical, intent(inout)             :: patchwi(ixG^T)

integer                            :: ix^D
!-----------------------------------------------------------------------------
{do ix^DB=ixOmin^DB,ixOmax^DB\}
   if(dabs(w(ix^D,b3_))>=1.38d0 .and. x(ix^D,2)>0.1d0) then
     patchwi(ix^D)=.true.
   else
     patchwi(ix^D)=.false.
   endif
{end do\}
return
end subroutine mask_gridfluxrope
!=============================================================================
function integral_grid(ixI^L,ixO^L,w,x,dvolume,iw,patchwi)

include 'amrvacdef.f'

integer, intent(in)                :: ixI^L,ixO^L,iw
double precision, intent(in)       :: x(ixI^S,1:ndim),dvolume(ixG^T)
double precision                   :: w(ixI^S,nw+nwauxio)
logical, intent(in) :: patchwi(ixG^T)

double precision, dimension(ixG^T,1:ndir) :: bvec,qvec,current
double precision :: integral_grid,tmp(ixG^T)
integer :: ix^D
!-----------------------------------------------------------------------------
integral_grid=0.d0
select case(iw)
 case(1)
  ! kinetic energy integration
  {do ix^DB=ixOmin^DB,ixOmax^DB\}
     integral_grid=integral_grid+half*(^C&w(ix^D,m0_+^C)**2+)/w(ix^D,rho_)*&
                   dvolume(ix^D)
  {end do\}
 case(2)
  ! total energy integration
  {do ix^DB=ixOmin^DB,ixOmax^DB\}
     integral_grid=integral_grid+w(ix^D,e_)*dvolume(ix^D)
  {end do\}
 case(9)
  ! magnetic energy integration
  {do ix^DB=ixOmin^DB,ixOmax^DB\}
     integral_grid=integral_grid+half*(^C&w(ix^D,b0_+^C)**2+)*dvolume(ix^D)
  {end do\}
 case(3)
  ! magnetic flux integration
  {do ix^DB=ixOmin^DB,ixOmax^DB\}
     if(patchwi(ix^D)) then
       integral_grid=integral_grid+w(ix^D,b3_)*dvolume(ix^D)
     endif
  {end do\}
end select
return
end function integral_grid
!=============================================================================
subroutine spacial_integral_w

include 'amrvacdef.f'

double precision :: dvolume(ixG^T), timephy
double precision, allocatable :: integral_ipe(:), integral_w(:)
double precision, external :: integral_grid
integer          :: iigrid,igrid,iw,amode,int_fh,status(MPI_STATUS_SIZE),ni
character(len=100):: filename,purpose
character(len=1024) :: line
logical          :: patchwi(ixG^T),alive
!-----------------------------------------------------------------------------
purpose='bflux'
!purpose='Energy'
ni=3
allocate(integral_ipe(ni),integral_w(ni))
integral_ipe=0.d0
integral_w=0.d0

do iigrid=1,igridstail; igrid=igrids(iigrid);
  if(slab) then
    dvolume(ixM^T)={rnode(rpdx^D_,igrid)|*}
  else
    dvolume(ixM^T)=pgeo(igrid)%dvolume(ixM^T)
  end if
  ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
  patchwi=.false.
  select case(purpose)
  case('bflux')
  call mask_gridfluxrope(ixG^LL,ixM^LL,pw(igrid)%w,px(igrid)%x,patchwi)
  ! total magnetic flux in the flux rope
  integral_ipe(1)=integral_ipe(1)+integral_grid(ixG^LL,ixM^LL,pw(igrid)%w,px(igrid)%x,&
    dvolume,3,patchwi)
  case('Energy')
  integral_ipe(1)=integral_ipe(1)+integral_grid(ixG^LL,ixM^LL,pw(igrid)%w,px(igrid)%x,&
    dvolume,2,patchwi)
  integral_ipe(2)=integral_ipe(2)+integral_grid(ixG^LL,ixM^LL,pw(igrid)%w,px(igrid)%x,&
    dvolume,9,patchwi)
  end select
end do
call MPI_ALLREDUCE(integral_ipe,integral_w,ni,MPI_DOUBLE_PRECISION,&
                     MPI_SUM,icomm,ierrmpi)
timephy=t*tunit
if(mype==0) then
  write(filename,"(a,a,a)") TRIM(filenamelog),TRIM(purpose),".int"
  inquire(file=filename,exist=alive)
  if(alive) then
    open(unit=11,file=filename,form='formatted',status='old',access='append')
  else
    open(unit=11,file=filename,form='formatted',status='new')
  endif
  select case(purpose)
   case('bflux')
    write(11,'(3(es12.4))') t,integral_w(1)
   case('mass')
    write(11,'(3(es12.4))') timephy,integral_w
   case('forcefreeness')
    write(11,'(2(es12.4))') timephy,integral_w(1)/integral_w(2)
   case('Energy')
    write(11,'(4(es12.4))') timephy,integral_w(1),integral_w(2),integral_w(3)
  endselect
  close(11)
endif
if(mype==0) print*,'t=',t,'integral_w=',integral_w
deallocate(integral_ipe,integral_w)
end subroutine spacial_integral_w
!=============================================================================
subroutine userspecialconvert(qunitconvert)

include 'amrvacdef.f'
integer, intent(in) :: qunitconvert
character(len=20):: userconvert_type
!-----------------------------------------------------------------------------
call spacial_integral_w
end subroutine userspecialconvert
!=============================================================================
subroutine boundarydriver(qt,qdt,pwa,pwb)
! prepare bottom boundary condition at highest resolution
use usr_bc
include 'amrvacdef.f'

double precision, intent(in)       :: qt,qdt
type(walloc), dimension(ngridshi)  :: pwa,pwb
double precision, allocatable :: wba(:^D&,:),wbb(:^D&,:),En(:^D&,:),curlEn(:^D&,:)
double precision :: dxa^D, xlen^D,startpos^D,ft,tramp1,tramp2,tfstop,rxy,invdx(1:ndim)
integer :: i^D,ix^D,ixbc^D, iigrid, igrid, idims, iside, af, idirmin
!----------------------------------------------------------------------------
if(qt>=tdstop) then
  vbt=0.d0
  return
endif
if(istep==1) Bbc=Bbt
allocate(wba(nxbc1,nxbc2+1,6))
allocate(wbb(nxbc1,nxbc2+1,3))
wba=0.d0
wbb=0.d0
dxa^D=dx(^D,mxnest);
startpos^D=xprobmin^D-dble(ngl)*dxa^D;\
do iigrid=1,igridstail; igrid=igrids(iigrid);
   ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
   do idims=1,ndim
      do iside=1,2
         i^D=kr(^D,idims)*(2*iside-3);
         if(neighbor_type(i^D,igrid)/=1) cycle
         if(idims==2 .and. iside==1) then
           af=nint(dxlevel(2)/dxa2)
           ix2=ixMlo2
           do ix1=ixMlo1,ixMhi1
             xlen^D=px(igrid)%x({ix^D,},^D)-startpos^D;\
             if(af>=2) then
               af=af/2
               ixbc^D=nint(xlen^D/dxa^D);\
               wba(ixbc^D-af+1:ixbc^D+af,1)=pwa(igrid)%w(ix^D,m1_)/pwa(igrid)%w(ix^D,rho_)
               wba(ixbc^D-af+1:ixbc^D+af,2)=pwa(igrid)%w(ix^D,m2_)/pwa(igrid)%w(ix^D,rho_)
               wba(ixbc^D-af+1:ixbc^D+af,3)=pwa(igrid)%w(ix^D,m3_)/pwa(igrid)%w(ix^D,rho_)
               wba(ixbc^D-af+1:ixbc^D+af,4)=pwa(igrid)%w(ix^D,b1_)
               wba(ixbc^D-af+1:ixbc^D+af,5)=pwa(igrid)%w(ix^D,b2_)
               wba(ixbc^D-af+1:ixbc^D+af,6)=pwa(igrid)%w(ix^D,b3_)

               wbb(ixbc^D-af+1:ixbc^D+af,1)=pwb(igrid)%w(ix^D,b1_)
               wbb(ixbc^D-af+1:ixbc^D+af,2)=pwb(igrid)%w(ix^D,b2_)
               wbb(ixbc^D-af+1:ixbc^D+af,3)=pwb(igrid)%w(ix^D,b3_)
             else
               ixbc^D=ceiling(xlen^D/dxa^D);\
               wba(ixbc^D,1)=pwa(igrid)%w(ix^D,m1_)/pwa(igrid)%w(ix^D,rho_)
               wba(ixbc^D,2)=pwa(igrid)%w(ix^D,m2_)/pwa(igrid)%w(ix^D,rho_)
               wba(ixbc^D,3)=pwa(igrid)%w(ix^D,m3_)/pwa(igrid)%w(ix^D,rho_)
               wba(ixbc^D,4)=pwa(igrid)%w(ix^D,b1_)
               wba(ixbc^D,5)=pwa(igrid)%w(ix^D,b2_)
               wba(ixbc^D,6)=pwa(igrid)%w(ix^D,b3_)

               wbb(ixbc^D,1)=pwb(igrid)%w(ix^D,b1_)
               wbb(ixbc^D,2)=pwb(igrid)%w(ix^D,b2_)
               wbb(ixbc^D,3)=pwb(igrid)%w(ix^D,b3_)
             endif
           end do
         end if
      end do
   end do
end do

call MPI_ALLREDUCE(MPI_IN_PLACE,wba,nxbc1*(nxbc2+1)*6,MPI_DOUBLE_PRECISION,&
                     MPI_SUM,icomm,ierrmpi)
call MPI_ALLREDUCE(MPI_IN_PLACE,wbb,nxbc1*(nxbc2+1)*3,MPI_DOUBLE_PRECISION,&
                     MPI_SUM,icomm,ierrmpi)

vbt(:,nxbc2+1,1)=wba(:,nxbc2+1,1)
vbt(:,nxbc2+1,2)=wba(:,nxbc2+1,2)
vbt(:,nxbc2+1,3)=wba(:,nxbc2+1,3)
Bbc(:,nxbc2+1,1)=wba(:,nxbc2+1,4)
Bbc(:,nxbc2+1,2)=wba(:,nxbc2+1,5)
Bbc(:,nxbc2+1,3)=wba(:,nxbc2+1,6)
Bbt(:,nxbc2+1,1)=wbb(:,nxbc2+1,1)
Bbt(:,nxbc2+1,2)=wbb(:,nxbc2+1,2)
Bbt(:,nxbc2+1,3)=wbb(:,nxbc2+1,3)
! 2nd order CD for divB=0 to set normal B component better
rxy=dxa2/dxa1
Bbc(2:nxbc1-1,nxbc2+1,2)=Bbc(2:nxbc1-1,nxbc2-1,2)-rxy*&
   (Bbc(3:nxbc1,nxbc2,1)-Bbc(1:nxbc1-2,nxbc2,1))
!do ix1=2,nxbc1-1
!  Bbc(ix1,nxbc2+1,2)=Bbc(ix1,nxbc2-1,2)-rxy*&
!    (Bbc(ix1+1,nxbc2,1)-Bbc(ix1-1,nxbc2,1))
!end do
deallocate(wba)
allocate(    En(nxbc1,nxbc2+2,3))
allocate(curlEn(nxbc1,nxbc2+2,3))
En=0.d0
curlEn=0.d0
ft=0.d0
tramp1=5.d2/tunit
tramp2=5.d2/tunit
tfstop=tdstop-tramp2
if(qt<tramp1) then
  ft=qt/tramp1
else if(qt<tfstop) then
  ft=1.d0
else if(qt<tfstop+tramp2) then
  ft=(tfstop+tramp2-qt)/tramp2
else
  ft=0.d0
endif
vbt(:,1:nxbc2,1)=vbc(:,1:nxbc2,1)*ft
vbt(:,1:nxbc2,3)=vbc(:,1:nxbc2,3)*ft
!do ix2=1,nxbc2
!  vbt(:,ix2,2)=vbt(:,nxbc2+1,2)
!end do
! test v=0
if(iprob==1) vbt=0.d0
! -E = v x B
En(1:nxbc1,2:nxbc2+2,1)=vbt(1:nxbc1,1:nxbc2+1,2)*Bbc(1:nxbc1,1:nxbc2+1,3)&
                       -vbt(1:nxbc1,1:nxbc2+1,3)*Bbc(1:nxbc1,1:nxbc2+1,2)
En(1:nxbc1,2:nxbc2+2,2)=vbt(1:nxbc1,1:nxbc2+1,3)*Bbc(1:nxbc1,1:nxbc2+1,1)&
                       -vbt(1:nxbc1,1:nxbc2+1,1)*Bbc(1:nxbc1,1:nxbc2+1,3)
En(1:nxbc1,2:nxbc2+2,3)=vbt(1:nxbc1,1:nxbc2+1,1)*Bbc(1:nxbc1,1:nxbc2+1,2)&
                       -vbt(1:nxbc1,1:nxbc2+1,2)*Bbc(1:nxbc1,1:nxbc2+1,1)
En(1:nxbc1,1,:)=En(1:nxbc1,2,:)
! curl (-E)
call curl(En,1,1,nxbc1,nxbc2+2,2,2,nxbc1-1,nxbc2+1,curlEn,idirmin,1,3,mxnest)
! B_N+1 = B_N + dt * curl (-E)
Bbt(ngl+1:nxbc1-ngl,1:nxbc2,:)=Bbc(ngl+1:nxbc1-ngl,1:nxbc2,:)+qdt*&
   curlEn(ngl+1:nxbc1-ngl,2:nxbc2+1,:)
Bbt(2:nxbc1-1,1:nxbc2,:)=Bbc(2:nxbc1-1,1:nxbc2,:)+qdt*&
   curlEn(2:nxbc1-1,2:nxbc2+1,:)
! 2nd order CD for divB=0 to set normal B component better
Bbt(2:nxbc1-1,nxbc2+1,2)=Bbt(2:nxbc1-1,nxbc2-1,2)-rxy*&
   (Bbt(3:nxbc1,nxbc2,1)-Bbt(1:nxbc1-2,nxbc2,1))
! average -E between t_n and t_n+1
En(1:nxbc1,2:nxbc2+2,1)=vbt(1:nxbc1,1:nxbc2+1,2)*Bbt(1:nxbc1,1:nxbc2+1,3)&
                       -vbt(1:nxbc1,1:nxbc2+1,3)*Bbt(1:nxbc1,1:nxbc2+1,2)&
                       +En(1:nxbc1,2:nxbc2+2,1)
En(1:nxbc1,2:nxbc2+2,2)=vbt(1:nxbc1,1:nxbc2+1,3)*Bbt(1:nxbc1,1:nxbc2+1,1)&
                       -vbt(1:nxbc1,1:nxbc2+1,1)*Bbt(1:nxbc1,1:nxbc2+1,3)&
                       +En(1:nxbc1,2:nxbc2+2,2)
En(1:nxbc1,2:nxbc2+2,3)=vbt(1:nxbc1,1:nxbc2+1,1)*Bbt(1:nxbc1,1:nxbc2+1,2)&
                       -vbt(1:nxbc1,1:nxbc2+1,2)*Bbt(1:nxbc1,1:nxbc2+1,1)&
                       +En(1:nxbc1,2:nxbc2+2,3)
En=0.5d0*En
En(:,1,:)=En(:,2,:)
! curl (-E)
call curl(En,1,1,nxbc1,nxbc2+2,2,2,nxbc1-1,nxbc2+1,curlEn,idirmin,1,3,mxnest)
! B_N+1 = B_N + dt * curl (-E)
Bbt(ngl+1:nxbc1-ngl,1:nxbc2,:)=Bbc(ngl+1:nxbc1-ngl,1:nxbc2,:)+qdt*&
   curlEn(ngl+1:nxbc1-ngl,2:nxbc2+1,:)
! 2nd order CD for divB=0 to set normal B component better
!rxy=dxa2/dxa1
!do i2=nxbc2,1,-1
!  Bbt(ngl+1:nxbc1-ngl,i2,2)=Bbt(ngl+1:nxbc1-ngl,i2+2,2)+rxy*&
!   (Bbt(ngl+2:nxbc1-ngl+1,i2+1,1)-Bbt(ngl:nxbc1-ngl-1,i2+1,1))
!end do

do ix1=1,ngl
  Bbt(ix1,:,1)=-Bbt(2*ngl-ix1+1,:,1)
  Bbt(ix1,:,2)= Bbt(2*ngl-ix1+1,:,2)
  Bbt(ix1,:,3)=-Bbt(2*ngl-ix1+1,:,3)
  Bbt(nxbc1-ngl+ix1,:,1)=-Bbt(nxbc1-ngl-ix1+1,:,1)
  Bbt(nxbc1-ngl+ix1,:,2)= Bbt(nxbc1-ngl-ix1+1,:,2)
  Bbt(nxbc1-ngl+ix1,:,3)=-Bbt(nxbc1-ngl-ix1+1,:,3)
end do
deallocate(wbb)
deallocate(En)
deallocate(curlEn)
end subroutine boundarydriver
!=============================================================================
subroutine clean_divb(bfield,ixI^L,ixO^L,level)

! Add Linde's divB related sources to wnew within ixO
include 'amrvacdef.f'

integer, intent(in)             :: ixI^L, ixO^L,level
double precision, intent(inout) :: bfield(ixI^S,1:ndir)
double precision :: divb(ixI^S),graddivb(ixI^S),bdivb(ixI^S,1:ndir),dxp2
integer :: iw, idims
!-----------------------------------------------------------------------------

! Calculate div B
divb=0.d0
call div(bfield,ixI^L,ixO^L,divb,level)

dxp2=(^D&1.0d0/dx(^D,level)**2+)
! Add Linde's diffusive terms
do idims=1,ndim
  ! Calculate grad_idim(divb)
  call grad(divb,ixI^L,ixO^L,idims,graddivb,level)
  ! Multiply by Linde's eta*dt = divbdiff*(c_max*dx)*dt = divbdiff*dx**2
  if (slab) then
     graddivb(ixO^S)=graddivb(ixO^S)*divbdiff/dxp2
  else
     graddivb(ixO^S)=graddivb(ixO^S)*divbdiff &
                     /(^D&1.0d0/mygeo%dx(ixO^S,^D)**2+)
  end if
  do iw= 1,ndir
    ! B_idim += eta*grad_idim(divb)
    bfield(ixO^S,iw)=bfield(ixO^S,iw)+graddivb(ixO^S)
  end do
end do

end subroutine clean_divb
!=============================================================================
subroutine grad(q,ixI^L,ix^L,idir,gradq,level)

! Calculate gradient of a scalar q within ixL in direction idir

include 'amrvacdef.f'

integer :: ixI^L,ix^L, idir,level
double precision :: q(ixI^S), gradq(ixI^S)

double precision :: qC(ixI^S),invdx
integer :: jx^L, hx^L, ixC^L, jxC^L 

!-----------------------------------------------------------------------------

invdx=1.d0/dx(idir,level)
if (slab) then
   jx^L=ix^L+kr(idir,^D);
   hx^L=ix^L-kr(idir,^D);
   gradq(ix^S) = half*(q(jx^S)-q(hx^S))*invdx
else
   hx^L=ix^L-kr(idir,^D);
   ixCmin^D=hxmin^D;ixCmax^D=ixmax^D;
   jxC^L=ixC^L+kr(idir,^D);
   select case(idir)
   {case(^D)
      qC(ixC^S)=mygeo%surfaceC^D(ixC^S)*half*(q(ixC^S)+q(jxC^S))
      gradq(ix^S)=(qC(ix^S)-qC(hx^S))/mygeo%dvolume(ix^S)
      ! Substract difference divergence and gradient
      gradq(ix^S)=gradq(ix^S)-q(ix^S) &
                     *(mygeo%surfaceC^D(ix^S)-mygeo%surfaceC^D(hx^S)) &
                    /mygeo%dvolume(ix^S) \}
   end select
end if

end subroutine grad
!=============================================================================
subroutine div(qvec,ixI^L,ixO^L,divq,level)

! Calculate divergence of a vector qvec within ixL

include 'amrvacdef.f'

integer :: ixI^L,ixO^L,level
double precision :: qvec(ixI^S,1:ndir), divq(ixI^S)

double precision :: qC(ixI^S), invdx(1:ndim)
integer :: jxO^L, hxO^L, ixC^L, jxC^L, idims, ix^L
!-----------------------------------------------------------------------------
ix^L=ixO^L^LADD1;
if (ixImin^D>ixmin^D.or.ixImax^D<ixmax^D|.or.) &
   call mpistop("Error in div: Non-conforming input limits")
invdx(1:ndim)=1.d0/dx(1:ndim,level)
divq(ixO^S)=zero
do idims=1,ndim
   if (slab) then
     jxO^L=ixO^L+kr(idims,^D);
     hxO^L=ixO^L-kr(idims,^D);
     divq(ixO^S)=divq(ixO^S)+half*(qvec(jxO^S,idims)-qvec(hxO^S,idims))*invdx(idims)
   else
     hxO^L=ixO^L-kr(idims,^D);
     ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
     jxC^L=ixC^L+kr(idims,^D);
     select case(idims)
     {case(^D)
        qC(ixC^S)=mygeo%surfaceC^D(ixC^S)*half*(qvec(ixC^S,idims)+qvec(jxC^S,idims))
        divq(ixO^S)=divq(ixO^S)+(qC(ixO^S)-qC(hxO^S))/mygeo%dvolume(ixO^S) \}
      end select
   end if
end do


end subroutine div 
!=============================================================================
subroutine curl(qvec,ixI^L,ixO^L,curlvec,idirmin,idirmin0,ndir0,level)

! Calculate curl of a vector qvec within ixL

include 'amrvacdef.f'

integer :: ixI^L,ixO^L,idirmin,ix^L,idir,jdir,kdir,hxO^L,jxO^L,ndir0,idirmin0,level
double precision :: qvec(ixI^S,1:ndir0),curlvec(ixI^S,idirmin0:3), invdx(1:ndim)
double precision :: tmp(ixI^S),tmp2(ixI^S),surface(ixI^S),mydx(ixI^S)
!-----------------------------------------------------------------------------
ix^L=ixO^L^LADD1;
if (ixImin^D>ixmin^D.or.ixImax^D<ixmax^D|.or.) &
   call mpistop("Error in curl: Non-conforming input limits")

! Calculate curl within ixL: CurlV_i=eps_ijk*d_j V_k
! Curl can have components (idirmin0:3)
! Determine exact value of idirmin while doing the loop.

invdx(:)=1.d0/dx(:,level)
idirmin=4
curlvec(ixO^S,idirmin0:3)=zero

do idir=idirmin0,3; do jdir=1,ndim; do kdir=1,ndir0
   if(lvc(idir,jdir,kdir)/=0)then
      tmp(ix^S)=qvec(ix^S,kdir)
      hxO^L=ixO^L-kr(jdir,^D);
      jxO^L=ixO^L+kr(jdir,^D);
      if(slab)then
         tmp2(ixO^S)=half*(tmp(jxO^S)-tmp(hxO^S))*invdx(jdir)
      else
         ! approximate formula, reduces to slab case
         ! and avoids staggering

         if (kdir .le. ndim) then 
            mydx(ix^S)=mygeo%dx(ix^S,kdir)
         else 
            mydx(ix^S)=one
         end if

         select case(idir)
           {case(^D)
             surface(ixO^S)=mygeo%surface^D(ixO^S)
             tmp2(ixO^S)=half*(mydx(jxO^S)*tmp(jxO^S) &
                              -mydx(hxO^S)*tmp(hxO^S)) &
                     /surface(ixO^S) \}
          end select
      endif
      if(lvc(idir,jdir,kdir)==1)then
         curlvec(ixO^S,idir)=curlvec(ixO^S,idir)+tmp2(ixO^S)
      else
         curlvec(ixO^S,idir)=curlvec(ixO^S,idir)-tmp2(ixO^S)
      endif
      if(idir<idirmin)idirmin=idir
   endif
enddo; enddo; enddo;

end subroutine curl 
!=============================================================================
subroutine write_boundary
! output boundary data for a restart
use usr_bc
include 'amrvacdef.f'

integer :: file_handle,amode
integer, dimension(MPI_STATUS_SIZE) :: statuss
character(len=190) :: boundaryfile
!----------------------------------------------------------------------------
if(mype==0) then
write(boundaryfile,"(a,i4.4,a)") TRIM(filenameout),snapshot,"boundary.dat"
  open(unit=37,file=boundaryfile,status='replace')
  close(37)
  amode=ior(MPI_MODE_CREATE,MPI_MODE_WRONLY)
  call MPI_FILE_OPEN(MPI_COMM_SELF,boundaryfile,amode, &
                       MPI_INFO_NULL,file_handle,ierrmpi)
  call MPI_FILE_WRITE(file_handle,Bbt,nxbc1*(nxbc2+1)*3,&
                       MPI_DOUBLE_PRECISION,statuss,ierrmpi)
  call MPI_FILE_CLOSE(file_handle,ierrmpi)
endif

end subroutine write_boundary
!=============================================================================
subroutine read_boundary
! read boundary data for a restart
use usr_bc
include 'amrvacdef.f'

integer :: file_handle,amode
integer, dimension(MPI_STATUS_SIZE) :: statuss
character(len=190) :: boundaryfile
logical :: bfexist
!----------------------------------------------------------------------------
write(boundaryfile,"(a,i4.4,a)") TRIM(filenameini),snapshotini,"boundary.dat"
inquire(file=boundaryfile, exist=bfexist)
if(.not. bfexist) then
  if (mype==0) write(*,*) "no runtime boundary data found!"
  return
end if

if(mype==0) then
  call MPI_FILE_OPEN(MPI_COMM_SELF,boundaryfile,MPI_MODE_RDONLY,MPI_INFO_NULL,&
                     file_handle,ierrmpi)
  call MPI_FILE_READ(file_handle,Bbt,nxbc1*(nxbc2+1)*3,&
                         MPI_DOUBLE_PRECISION,statuss,ierrmpi)
  call MPI_FILE_CLOSE(file_handle,ierrmpi)
end if
if(npe>1)then
  call MPI_BCAST(Bbt,nxbc1*(nxbc2+1)*3,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
end if
end subroutine read_boundary
{#IFDEF FCT
!=============================================================================
subroutine initvecpot_usr(ixI^L, ixC^L, xC, A)

! initialize the vectorpotential on the corners
! used by b_from_vectorpotential()


include 'amrvacdef.f'

integer, intent(in)                :: ixI^L, ixC^L
double precision, intent(in)       :: xC(ixI^S,1:ndim)
double precision, intent(out)      :: A(ixI^S,1:ndir)

!-----------------------------------------------------------------------------

A(ixC^S,1:ndir) = zero

end subroutine initvecpot_usr
}
!=============================================================================
! amrvacusr.t.fluxrope23
!=============================================================================
