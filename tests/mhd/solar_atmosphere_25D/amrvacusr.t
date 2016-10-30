!=============================================================================
! amrvacusr.t.solaratmosphere23
!=============================================================================
INCLUDE:amrvacmodules/cooling.t
INCLUDE:amrvacmodules/thermalconduction.t
INCLUDE:amrvacnul/usrflags.t
!=============================================================================
module usr_bc
implicit none
double precision, allocatable, save :: pbc(:),rbc(:)
end module usr_bc
!=============================================================================
subroutine initglobaldata_usr
use constants
use mod_global_parameters

double precision :: miu0
logical, save :: firstusrglobaldata=.true.
!-----------------------------------------------------------------------------
! CGS Unit
miu0=4.d0*dpi      ! permeability of free space in Gauss^2 cm^2 dyne^-1

unit_length=1.d9         ! length in cm
unit_temperature=1.d6        ! temperature in K
unit_numberdensity=1.d9        ! number density in cm^-3

unit_density=1.4d0*const_mp*unit_numberdensity     ! 2.341668000000000E-015 g*cm^-3
unit_pressure=2.3d0*unit_numberdensity*const_kB*unit_temperature ! 0.317538000000000 erg*cm^-3
unit_magneticfield=dsqrt(miu0*unit_pressure)       ! 1.99757357615242 Gauss
unit_velocity=unit_magneticfield/dsqrt(miu0*unit_density) ! 1.16448846777562E007 cm/s = 116.45 km/s
unit_time=unit_length/unit_velocity             ! 85.8746159942810 s
heatunit=unit_pressure/unit_time          ! 3.697693390805347E-003 erg*cm^-3/s

eqpar(grav1_)=0.d0
eqpar(grav2_)=-2.74d4*unit_length/unit_velocity**2 ! solar gravity
eqpar(eta_)=0.d0 ! resistivity
! cooling table parameters
eqpar(Tscale_)=1.d0/unit_temperature
eqpar(Lscale_)=unit_density*unit_length/unit_velocity**3/const_mp**2*1.2d0/1.4d0**2
eqpar(Mue_)=one
bQ0=1.d-4/heatunit ! background heating power density
gzone=0.2d0 ! thickness of a ghostzone below the bottom boundary
dya=(2.d0*gzone+xprobmax2-xprobmin2)/dble(jmax) ! cells size of high-resolution 1D solar atmosphere
B0=20.d0/unit_magneticfield ! magnetic field strength at the bottom
theta=60.d0*dpi/180.d0 ! the angle to the plane xy, 90-theta is the angle to the polarity inversion line of the arcade 
kx=dpi/(xprobmax1-xprobmin1)
ly=kx*dcos(theta)
SRadius=69.61d0 ! Solar radius

if(firstusrglobaldata) then
  ! initialize thermal conduction
  call init_thermalconduction
  ! cooling module initalization
  call coolinit
  ! setup global gravity stratification in one vertical line
  call inithdstatic
  firstusrglobaldata=.false.
endif

end subroutine initglobaldata_usr
!=============================================================================
subroutine inithdstatic
!! initialize the table in a vertical line through the global domain
use usr_bc
use mod_thermalconduction, only:kappa
use mod_global_parameters

integer :: j,na,nb,ibc
double precision:: Ta(jmax),gg(jmax)
double precision:: rpho,Ttop,Tpho,wtra,res,rhob,pb,htra,Ttr,Fc,invT
!----------------------------------------------------------------------------
rpho=1.151d15/unit_numberdensity ! number density at the bottom relaxla
Tpho=8.d3/unit_temperature !< temperature of chromosphere
Ttop=1.5d6/unit_temperature !< estimated temperature in the top
htra=0.2d0 !< height of initial transition region
wtra=0.02d0 !< width of initial transition region 
Ttr=1.6d5/unit_temperature !< lowest temperature of upper profile
Fc=2.d5/heatunit/unit_length !< constant thermal conduction flux
do j=1,jmax
   ya(j)=(dble(j)-0.5d0)*dya-gzone
   if(ya(j)>htra) then
     Ta(j)=(3.5d0*Fc/kappa*(ya(j)-htra)+Ttr**3.5d0)**(2.d0/7.d0)
   else
     Ta(j)=Tpho+0.5d0*(Ttop-Tpho)*(tanh((ya(j)-htra-0.027d0)/wtra)+1.d0)
   endif
   gg(j)=eqpar(grav2_)*(SRadius/(SRadius+ya(j)))**2
enddo
!! solution of hydrostatic equation 
nb=int(gzone/dya)
ra(1)=rpho
pa(1)=rpho*Tpho
invT=gg(1)/Ta(1)
invT=0.d0
do j=2,jmax
   invT=invT+(gg(j)/Ta(j)+gg(j-1)/Ta(j-1))*0.5d0
   pa(j)=pa(1)*dexp(invT*dya)
   ra(j)=pa(j)/Ta(j)
end do
!do j=2,jmax
!   pa(j)=(pa(j-1)+dya*(gg(j)+gg(j-1))*ra(j-1)/4.d0)/(one-dya*(gg(j)+gg(j-1))/Ta(j)/4.d0)
!   ra(j)=pa(j)/Ta(j)
!end do
!! initialized rho and p in the fixed bottom boundary
na=floor(gzone/dya+0.5d0)
res=gzone-(dble(na)-0.5d0)*dya
rhob=ra(na)+res/dya*(ra(na+1)-ra(na))
pb=pa(na)+res/dya*(pa(na+1)-pa(na))
allocate(rbc(dixB))
allocate(pbc(dixB))
do ibc=dixB,1,-1
  na=floor((gzone-dx(2,mxnest)*(dble(dixB-ibc+1)-0.5d0))/dya+0.5d0)
  res=gzone-dx(2,mxnest)*(dble(dixB-ibc+1)-0.5d0)-(dble(na)-0.5d0)*dya
  rbc(ibc)=ra(na)+res/dya*(ra(na+1)-ra(na))
  pbc(ibc)=pa(na)+res/dya*(pa(na+1)-pa(na))
end do

if (mype==0) then
 print*,'minra',minval(ra)
 print*,'rhob',rhob
 print*,'pb',pb
endif
    
end subroutine inithdstatic
!=============================================================================
subroutine initonegrid_usr(ixG^L,ix^L,w,x)

! initialize one grid

use mod_global_parameters

integer, intent(in) :: ixG^L,ix^L
double precision, intent(in) :: x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)

double precision :: res
integer :: ix^D,na

logical, save :: first=.true.
!----------------------------------------------------------------------------

{^IFONED   call mpistop("prob prominence23 is 2.5D")}
{^IFTHREED call mpistop("prob prominence23 is 2.5D")}
{^IFGLM
w(ixG^S,psi_)=zero
}
{^IFTWOD
if (ndir/=3) call mpistop("this is a mhd23 problem!")
if(first)then
  if(mype==0) then
    write(*,*)'Simulating 2.5D solar atmosphere'
  endif
  first=.false.
endif
{do ix^DB=ixGmin^DB,ixGmax^DB\}
    na=floor((x(ix^D,2)-xprobmin2+gzone)/dya+0.5d0)
    res=x(ix^D,2)-xprobmin2+gzone-(dble(na)-0.5d0)*dya
    w(ix^D,rho_)=ra(na)+(one-cos(dpi*res/dya))/two*(ra(na+1)-ra(na))
    w(ix^D,p_)  =pa(na)+(one-cos(dpi*res/dya))/two*(pa(na+1)-pa(na))
{end do\}
w(ixG^S,v1_)=zero
w(ixG^S,v2_)=zero
w(ixG^S,v3_)=zero
if(B0field) then
  w(ixG^S,b1_)=zero
  w(ixG^S,b2_)=zero
  w(ixG^S,b3_)=zero
else
  w(ixG^S,b1_)=-B0*dcos(kx*x(ixG^S,1))*dexp(-ly*x(ixG^S,2))*dcos(theta)
  w(ixG^S,b2_)= B0*dsin(kx*x(ixG^S,1))*dexp(-ly*x(ixG^S,2))
  w(ixG^S,b3_)=-B0*dcos(kx*x(ixG^S,1))*dexp(-ly*x(ixG^S,2))*dsin(theta)
endif
call conserve(ixG^L,ixG^L,w,x,patchfalse)
\}

end subroutine initonegrid_usr
!=============================================================================
subroutine specialbound_usr(qt,ixG^L,ixO^L,iw,iB,w,x)

! special boundary types, user defined

use usr_bc
use mod_global_parameters

integer, intent(in) :: ixO^L, iw, iB, ixG^L
double precision, intent(in) :: qt, x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)

double precision :: pth(ixG^S),tmp(ixG^S),ggrid(ixG^S),invT(ixG^S)
double precision :: delydelx
integer :: ix^D,idims,ixInt^L
!----------------------------------------------------------------------------

select case(iB)
case(3)
  !! fixed zero velocity
  w(ixO^S,v1_) =-w(ixOmin1:ixOmax1,ixOmax2+dixB:ixOmax2+1:-1,m1_)&
               /w(ixOmin1:ixOmax1,ixOmax2+dixB:ixOmax2+1:-1,rho_)
  w(ixO^S,v2_) =-w(ixOmin1:ixOmax1,ixOmax2+dixB:ixOmax2+1:-1,m2_)&
               /w(ixOmin1:ixOmax1,ixOmax2+dixB:ixOmax2+1:-1,rho_)
  w(ixO^S,v3_) =-w(ixOmin1:ixOmax1,ixOmax2+dixB:ixOmax2+1:-1,m3_)&
               /w(ixOmin1:ixOmax1,ixOmax2+dixB:ixOmax2+1:-1,rho_)
  !! fixed b1 b2 b3
  if(iprob==0 .or. B0field) then
    w(ixO^S,b1_)=0.d0
    w(ixO^S,b2_)=0.d0
    w(ixO^S,b3_)=0.d0
  else
    w(ixO^S,b1_)=-B0*dcos(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))*dcos(theta)
    w(ixO^S,b2_)= B0*dsin(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))
    w(ixO^S,b3_)=-B0*dcos(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))*dsin(theta)
  endif
  !! fixed gravity stratification of density and pressure pre-determined in initial condition
  do ix2=ixOmin2,ixOmax2
    w(ixOmin1:ixOmax1,ix2,rho_)=rbc(ix2)
    w(ixOmin1:ixOmax1,ix2,p_)=pbc(ix2)
  enddo
  call conserve(ixG^L,ixO^L,w,x,patchfalse)
case(4)
  ixInt^L=ixO^L;
  ixIntmin2=ixOmin2-1;ixIntmax2=ixOmin2-1;
  call getpthermal(w,x,ixG^L,ixInt^L,pth)
  ixIntmin2=ixOmin2-1;ixIntmax2=ixOmax2;
  call getggrav(ggrid,ixG^L,ixInt^L,x)
  !> fill pth, rho ghost layers according to gravity stratification
  invT(ixOmin2-1^%2ixO^S)=w(ixOmin2-1^%2ixO^S,rho_)/pth(ixOmin2-1^%2ixO^S)
  tmp=0.d0
  do ix2=ixOmin2,ixOmax2
    tmp(ixOmin2-1^%2ixO^S)=tmp(ixOmin2-1^%2ixO^S)+0.5d0*&
        (ggrid(ix2^%2ixO^S)+ggrid(ix2-1^%2ixO^S))*invT(ixOmin2-1^%2ixO^S)
    w(ix2^%2ixO^S,p_)=pth(ixOmin2-1^%2ixO^S)*dexp(tmp(ixOmin2-1^%2ixO^S)*dxlevel(2))
    w(ix2^%2ixO^S,rho_)=w(ix2^%2ixO^S,p_)*invT(ixOmin2-1^%2ixO^S)
  enddo
  !> fixed zero velocity
  w(ixO^S,v1_) =-w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-dixB:-1,m1_)&
               /w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-dixB:-1,rho_)
  w(ixO^S,v2_) =-w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-dixB:-1,m2_)&
               /w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-dixB:-1,rho_)
  w(ixO^S,v3_) =-w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-dixB:-1,m3_)&
               /w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-dixB:-1,rho_)
  !> zero normal gradient extrapolation
  do ix2=ixOmin2,ixOmax2
    w(ixOmin1:ixOmax1,ix2,b1_:b3_)=(1.0d0/3.0d0)* &
                (-w(ixOmin1:ixOmax1,ix2-2,b1_:b3_)&
           +4.0d0*w(ixOmin1:ixOmax1,ix2-1,b1_:b3_))
  enddo
  call conserve(ixG^L,ixO^L,w,x,patchfalse)
case default
   call mpistop("Special boundary is not defined for this region")
end select

end subroutine specialbound_usr
!=============================================================================
subroutine getggrav(ggrid,ixI^L,ixO^L,x)
!> calculate gravity

use mod_global_parameters

integer, intent(in) :: ixI^L,ixO^L
double precision, intent(in) :: x(ixI^S,1:ndim)
double precision, intent(out) :: ggrid(ixI^S)
!---------------------------------------------------------------------------
ggrid(ixO^S)=eqpar(grav2_)*(SRadius/(SRadius+x(ixO^S,2)))**2
end subroutine getggrav
!=============================================================================
subroutine addsource_gravSA(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

! w[iw]=w[iw]+qdt*S[wCT,qtC,x] where S is the source based on wCT within ixO
! gravity distribution along a magnetic loop (a circular arc)

use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L, iw^LIM
double precision, intent(in) :: qdt, qtC, qt, x(ixI^S,1:ndim), wCT(ixI^S,1:nw)
double precision, intent(inout) :: w(ixI^S,1:nw)

double precision :: ggrid(ixI^S)
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
   case (e_)
     ! de/dt= +g_i*m_i
      do idims=1,ndim
         if (abs(eqpar(grav0_+idims))>smalldouble) &
            w(ixO^S,ee_)=w(ixO^S,ee_) &
              +qdt*ggrid(ixO^S)*wCT(ixO^S,m0_+idims)
      end do
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
double precision, intent(in) :: qdt, qtC, qt, x(ixI^S,1:ndim), wCT(ixI^S,1:nw)
double precision, intent(inout) :: w(ixI^S,1:nw)

double precision :: lQgrid(ixI^S),bQgrid(ixI^S)
!-----------------------------------------------------------------------------
! add gravity source
call addsource_gravSA(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)
! add radiative cooling source
call addsource_cooling(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)
if(B0field) then
! add source terms when splitting linear force-free field
  call addsource_lfff(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)
endif 
! add global background heating bQ
call getbQ(bQgrid,ixI^L,ixO^L,qtC,wCT,x)
w(ixO^S,e_)=w(ixO^S,e_)+qdt*bQgrid(ixO^S)

end subroutine specialsource
!=============================================================================
subroutine addsource_lfff(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L, iw^LIM
double precision, intent(in) :: qdt, qtC, qt, wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
double precision, intent(inout) :: w(ixI^S,1:nw)

double precision :: current0(ixI^S,1:ndir),v(ixI^S,1:ndir)
!-----------------------------------------------------------------------------
! electric current density of the background linear force-free field
current0(ixO^S,1)= ly*B0*dcos(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))*dsin(theta)
current0(ixO^S,2)=-kx*B0*dsin(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))*dsin(theta)
current0(ixO^S,3)= kx*B0*dcos(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))&
                  -ly*B0*dcos(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))*dcos(theta)

v(ixO^S,1)=wCT(ixO^S,m1_)/wCT(ixO^S,rho_)
v(ixO^S,2)=wCT(ixO^S,m2_)/wCT(ixO^S,rho_)
v(ixO^S,3)=wCT(ixO^S,m3_)/wCT(ixO^S,rho_)

w(ixO^S,e_)=w(ixO^S,e_)-&
qdt*((v(ixO^S,2)*wCT(ixO^S,b3_)-v(ixO^S,3)*wCT(ixO^S,b2_))*current0(ixO^S,1)&
    +(v(ixO^S,3)*wCT(ixO^S,b1_)-v(ixO^S,1)*wCT(ixO^S,b3_))*current0(ixO^S,2)&
    +(v(ixO^S,1)*wCT(ixO^S,b2_)-v(ixO^S,2)*wCT(ixO^S,b1_))*current0(ixO^S,3))

end subroutine addsource_lfff
!=============================================================================
subroutine getdt_special(w,ixG^L,ix^L,dtnew,dx^D,x)

! Limit "dt" further if necessary, e.g. due to the special source terms.
! The getdt_courant (CFL condition) and the getdt subroutine in the AMRVACPHYS
! module have already been called.

use mod_global_parameters

integer, intent(in) :: ixG^L,ix^L
double precision, intent(in) :: dx^D, x(ixG^S,1:ndim)
! note that depending on strictsmall etc, w values may change
double precision, intent(inout) :: w(ixG^S,1:nw), dtnew
!-----------------------------------------------------------------------------
dtnew=bigdouble

call getdt_grav(w,ixG^L,ix^L,dtnew,dx^D,x)

call getdt_cooling(w,ixG^L,ix^L,dtnew,dx^D,x)
end subroutine getdt_special
!=============================================================================
subroutine getbQ(bQgrid,ixI^L,ixO^L,qt,w,x)
! calculate background heating bQ
use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L
double precision, intent(in) :: qt, x(ixI^S,1:ndim), w(ixI^S,1:nw)

double precision :: bQgrid(ixI^S)
!-----------------------------------------------------------------------------

bQgrid(ixO^S)=bQ0*dexp(-x(ixO^S,2)/5.d0)

end subroutine getbQ
!=============================================================================
subroutine specialeta(w,ixI^L,ix^L,idirmin,x,current,eta)

! Set the common "eta" array for resistive MHD based on w or the
! "current" variable which has components between idirmin and 3.

use mod_global_parameters

integer, intent(in) :: ixI^L, ix^L, idirmin
double precision, intent(in) :: w(ixI^S,nw), x(ixI^S,1:ndim)

double precision :: current(ixI^S,7-2*ndir:3), eta(ixI^S)
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
! fix the bottom layer to the highest level
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

use mod_global_parameters

integer, intent(in)          :: ixI^L,ixO^L,iflag
double precision, intent(in) :: w(ixI^S,1:nw)
double precision, intent(out):: var(ixI^S)
!-----------------------------------------------------------------------------

if (iflag >nw)call mpistop(' iflag> nw, make change in parfile or in user file')

var(ixI^S) = zero

end subroutine specialvarforerrest
!=============================================================================
subroutine printlog_special

use mod_global_parameters
!-----------------------------------------------------------------------------

call mpistop("special log file undefined")

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

double precision :: tmp(ixI^S),tmp1(ixI^S),tmp2(ixI^S)
double precision :: qvec(ixI^S,1:ndir),curlvec(ixI^S,1:ndir)
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

use mod_global_parameters

integer, intent(in)                :: ixI^L,ixO^L
double precision, intent(in)       :: x(ixI^S,1:ndim)
double precision                   :: w(ixI^S,nw+nwauxio)
double precision                   :: normconv(0:nw+nwauxio)
double precision :: tmp(ixI^S),tmp1(ixI^S),tmp2(ixI^S),dRdT(ixI^S)
double precision :: lQgrid(ixI^S),bQgrid(ixI^S),divb(ixI^S)

double precision :: qvec(ixI^S,1:ndir),curlvec(ixI^S,1:ndir)
integer :: idirmin,idims,ix^D
!-----------------------------------------------------------------------------
! output temperature
call getpthermal(w,x,ixI^L,ixO^L,tmp)
w(ixO^S,nw+1)=tmp(ixO^S)/w(ixO^S,rho_)
! output Alfven wave speed
if(B0field)then
  w(ixO^S,nw+2)=dsqrt((^C&(w(ixO^S,b^C_)+myB0_cell%w(ixO^S,^C))**2+)/w(ixO^S,rho_))
else
  w(ixO^S,nw+2)=dsqrt((^C&w(ixO^S,b^C_)**2+)/w(ixO^S,rho_))
endif

if(B0field) then
 ^C&qvec(ixI^S,^C)=w(ixI^S,b^C_)+myB0_cell%w(ixI^S,^C);
else
 ^C&qvec(ixI^S,^C)=w(ixI^S,b^C_);
endif

! output divB1
call getdivb(w,ixI^L,ixO^L,divb)
w(ixO^S,nw+3)=0.5d0*divb(ixO^S)/dsqrt(^C&qvec(ixO^S,^C)**2+)/(^D&1.0d0/dxlevel(^D)+)
! output the plasma beta p*2/B**2
if(B0field)then
  w(ixO^S,nw+4)=tmp(ixO^S)*two/(^C&(w(ixO^S,b^C_)+myB0_cell%w(ixO^S,^C))**2+)
else
  w(ixO^S,nw+4)=tmp(ixO^S)*two/(^C&w(ixO^S,b^C_)**2+)
endif
! output heating rate
call getbQ(bQgrid,ixI^L,ixO^L,t,w,x)
w(ixO^S,nw+5)=bQgrid(ixO^S)
! store the cooling rate 
call getvar_cooling(ixI^L,ixO^L,w,x,lQgrid,normconv)
w(ixO^S,nw+6)=lQgrid(ixO^S)

! store current
call curlvector(qvec,ixI^L,ixO^L,curlvec,idirmin,1,ndir)
^C&w(ixO^S,nw+6+^C)=curlvec(ixO^S,^C);

end subroutine specialvar_output
!=============================================================================
subroutine specialvarnames_output

! newly added variables need to be concatenated with the wnames/primnames string

use mod_global_parameters
!-----------------------------------------------------------------------------

primnames=TRIM(primnames)//' '//'Te Valfv divB beta bQ rad j1 j2 j3'
   wnames=   TRIM(wnames)//' '//'Te Valfv divB beta bQ rad j1 j2 j3'

end subroutine specialvarnames_output
!=============================================================================
subroutine specialset_B0(ixI^L,ixO^L,x,wB0)

! Here add a steady (time-independent) potential or linear force-free background field

use mod_global_parameters

integer, intent(in)           :: ixI^L,ixO^L
double precision, intent(in)  :: x(ixI^S,1:ndim)
double precision, intent(inout) :: wB0(ixI^S,1:ndir)
!-----------------------------------------------------------------------------

wB0(ixO^S,1)=wB0(ixO^S,1)-B0*dcos(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))*dcos(theta)
wB0(ixO^S,2)=wB0(ixO^S,2)+B0*dsin(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))
wB0(ixO^S,3)=wB0(ixO^S,3)-B0*dcos(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))*dsin(theta)

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

use mod_global_parameters

integer, intent(in) :: ixG^L,ixO^L,level
double precision, intent(in) :: qt
double precision, intent(inout) :: w(ixG^S,1:nw)
double precision, intent(in) :: x(ixG^S,1:ndim)

! .. local ..
!logical :: patchw(ixG^S)
!----------------------------------------------------------------------------

call mpistop("bc_int not defined")

end subroutine bc_int
!=============================================================================
! amrvacusr.t.solaratmosphere23
!=============================================================================
