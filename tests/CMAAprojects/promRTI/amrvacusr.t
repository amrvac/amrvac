!#############################################################################
! module amrvacusr - promRTideal

INCLUDE:amrvacmodules/gravity.t
INCLUDE:amrvacnul/specialimpl.t
INCLUDE:amrvacnul/usrflags.t
INCLUDE:amrvacmodules/handle_particles.t
INCLUDE:amrvacmodules/integrate_particles.t
!=============================================================================
module usr_bc
implicit none
double precision, allocatable, save :: pbc(:),rbc(:)
end module usr_bc
!=============================================================================
subroutine initglobaldata_usr

use mod_global_parameters
logical, save :: firstinitglobal=.true.
!-----------------------------------------------------------------------------
! CGS unit
k_B=1.3806d-16    ! erg*K^-1
miu0=4.d0*dpi     ! Gauss^2 cm^2 dyne^-1
Lunit=1.d9        ! cm = 10 Mm
UNIT_LENGTH=Lunit
Teunit=1.d6       ! K
nHunit=1.d9       ! cm^-3
mHunit=1.67262d-24 ! g
runit=1.4d0*mHunit*nHunit ! 2.341668000000000E-015 g*cm^-3
UNIT_DENSITY=runit
punit=2.3d0*nHunit*k_B*Teunit ! 0.317538000000000 erg*cm^-3
Bunit=dsqrt(miu0*punit)      ! 1.99757357615242 Gauss
vunit=Bunit/dsqrt(miu0*runit) ! 1.16448846777562E007 cm/s = 116.45 km/s
UNIT_VELOCITY=vunit
tunit=Lunit/vunit ! 85.8746159942810 s
heatunit=punit/tunit ! 3.697693390805347E-003 erg*cm^-3/s

! units for convert
normvar(0) = UNIT_LENGTH
normvar(rho_) = UNIT_DENSITY
{^C&normvar(v^C_)   = UNIT_VELOCITY \}
normvar(pp_)     = UNIT_VELOCITY**2 * UNIT_DENSITY
{^C&normvar(b^C_)   = dsqrt(4.0d0*dpi*normvar(pp_)) \}
normt = UNIT_LENGTH/UNIT_VELOCITY

eqpar(grav1_)=0.d0
eqpar(grav2_)=-2.74d4*Lunit/vunit**2
gzone=dixB*0.1875d0
dr=(2.d0*gzone+xprobmax2-xprobmin2)/dble(jmax)
! Solar radius
SRadius=6.961d10/Lunit
! base density and temperature
rho0=1.d18/nHunit
Tch=8.d3/Teunit
Tco=1.3d6/Teunit
B0=2.d0
bsca=0.3d0
if(firstinitglobal) then
  call inithdstatic
  firstinitglobal=.false.
endif
{#IFDEF GLM
eqpar(Cr_)=-0.2d0
}
        
end subroutine initglobaldata_usr
!=============================================================================
subroutine inithdstatic
!! initialize the table in a vertical line through the global domain
use usr_bc
use mod_global_parameters

integer :: j,na,ibc
double precision:: Ta(jmax),gg(jmax)
double precision:: res,rhob,pb,htra1,wtra1,wtra2,b3
!----------------------------------------------------------------------------
htra1=0.25d0
wtra1=0.01d0
htra2=0.55d0
wtra2=0.01d0
do j=1,jmax
   ya(j)=(dble(j)-0.5d0)*dr-gzone+xprobmin2
   Ta(j)=Tch+0.5d0*(Tco-Tch)*(tanh((ya(j)-htra1)/wtra1)+1.d0)
   if(ya(j)>0.4d0) then
     Ta(j)=Tco-0.5d0*(Tco-Tch)*(tanh((ya(j)-htra2)/wtra2)+1.d0)
   endif
   gg(j)=eqpar(grav2_)
enddo
!! solution of hydrostatic equation 
ra(1)=rho0
pa(1)=rho0*Tch
do j=2,jmax
  if(ya(j)<htra2) then
  pa(j)=(pa(j-1)+dr*(gg(j)+gg(j-1))*ra(j-1)/4.d0)/(one-dr*(gg(j)+gg(j-1))/&
          Ta(j)/4.d0)
!do j=jmax/2+1,jmax
!  b3=B0*dexp(-(ya(j-1)-ya(jmax/2))/bsca)+
!  pa(j)=2.d0*dr*(ra(j-1)*gg(j-1)+b3**2/bsca)+pa(j-2)
!  ra(j)=pa(j)/Ta(j)
!end do
  else
  b3=(B0*dexp(-(ya(j-1)-htra2)/bsca))**2+(B0*dexp(-(ya(j)-htra2)/bsca))**2
  pa(j)=(pa(j-1)+dr*(gg(j)+gg(j-1))*ra(j-1)/4.d0+dr/bsca*b3/2.d0)/(one-dr*(gg(j)+gg(j-1))/&
          Ta(j)/4.d0)
  endif
  ra(j)=pa(j)/Ta(j)
end do
!! initialized rho and p in the fixed bottom boundary
na=floor(gzone/dr+0.5d0)
res=gzone-(dble(na)-0.5d0)*dr
rhob=ra(na)+res/dr*(ra(na+1)-ra(na))
pb=pa(na)+res/dr*(pa(na+1)-pa(na))
allocate(rbc(dixB))
allocate(pbc(dixB))
do ibc=dixB,1,-1
  na=floor((gzone-dx(2,mxnest)*(dble(dixB-ibc+1)-0.5d0))/dr+0.5d0)
  res=gzone-dx(2,mxnest)*(dble(dixB-ibc+1)-0.5d0)-(dble(na)-0.5d0)*dr
  rbc(ibc)=ra(na)+res/dr*(ra(na+1)-ra(na))
  pbc(ibc)=pa(na)+res/dr*(pa(na+1)-pa(na))
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

integer, intent(in) :: ixG^L, ix^L
double precision, intent(in) :: x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)

double precision:: res,y0,pbm,ran,kx,pint,dely,sigma,Av
integer :: ix^D,ndw,na,l,randarr(2)
logical, save:: first=.true.
logical :: patchw(ixG^S)
!----------------------------------------------------------------------------

if (first) then
   if (mype==0) then
      print *,'2.5D MHD Rayleigh-Taylor instability in prominence'
      print *,'  --interface y0:',y0
!      print *,'  --density ratio:',rhodens/rholight
!      print *,'  --disturb wave number:',ndw
   end if
   first=.false.
end if
!w(ixG^S,v2_)=(dsin(kx*(x(ixG^S,1)-xprobmin1))+dsin(kx/2.5d0*(x(ixG^S,1)-xprobmin1))&
!   +dsin(kx/3.5d0*(x(ixG^S,1)-xprobmin1)))/3.d0
!w(ixG^S,v2_)=pbm*w(ixG^S,v2_)*dexp(-(dabs(x(ixG^S,2)-y0)/sigma)**2)
!w(ixG^S,v1_)=0.d0
y0=0.55d0+xprobmin2
sigma=0.02d0
pbm=0.012d0
pint=4.d0
ndw=12
kx=dble(ndw)*2.d0*dpi/(xprobmax1-xprobmin1)


Av=1.d0
w(ixG^S,v1_)=0.d0
{do ix^DB=ixmin^DB,ixmax^DB\}
   na=floor((x(ix^D,2)-xprobmin2+gzone)/dr+0.5d0)
   res=x(ix^D,2)-xprobmin2+gzone-(dble(na)-0.5d0)*dr
   w(ix^D,rho_)=ra(na)+(one-cos(dpi*res/dr))/two*(ra(na+1)-ra(na))
   w(ix^D,p_)  =pa(na)+(one-cos(dpi*res/dr))/two*(pa(na+1)-pa(na))
   do l=1,20
     randarr(1)=l
     randarr(2)=int(x(ix^D,1))
     call random_seed(put=randarr)
     call random_number(ran)
     w(ix^D,v1_)=w(ix^D,v1_)+Av*(-1.d0)**l*dsin(2.d0*dpi*x(ix^D,1)/(dble(l)+ran/2.d0))
   end do
   w(ix^D,v1_)=pbm*w(ix^D,v1_)*dexp(-(dabs(x(ix^D,2)-y0)/sigma)**2)
{end do\}
w(ixG^S,v2_)=0.d0
w(ixG^S,v3_)=0.d0

!where(x(ixG^S,2)>y0)
!   w(ixG^S,rho_)=rhodens
!elsewhere
!   w(ixG^S,rho_)=rholight
!endwhere
!w(ixG^S,p_)=pint+eqpar(grav2_)*w(ixG^S,rho_)*(x(ixG^S,2)-y0)
!where(w(ixG^S,p_)<0.3d0) 
!  w(ixG^S,p_)=0.3d0
!endwhere
if(B0field)then
 w(ix^S,b1_)  =zero
 w(ix^S,b2_)  =zero
 w(ix^S,b3_)  =zero
else
 w(ix^S,b1_)  =0.1d0
! w(ix^S,b1_)  =zero
 w(ix^S,b2_)  =zero
 where(x(ix^S,2)>htra2)
   w(ix^S,b3_)  =B0*dexp(-(x(ix^S,2)-htra2)/bsca)
 elsewhere
   w(ix^S,b3_)  =B0
 endwhere
endif
patchw(ixG^S)=.false.
call conserve(ixG^L,ix^L,w,x,patchw)
{#IFDEF GLM
w(ixG^S,psi_)=0.d0
}

end subroutine initonegrid_usr
!=============================================================================
subroutine specialbound_usr(qt,ixG^L,ixO^L,iw,iB,w,x)
use usr_bc
! special boundary types, user defined

use mod_global_parameters

integer, intent(in) :: ixO^L, iw, iB, ixG^L
double precision, intent(in) :: qt, x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)

double precision :: dx^D,delydelx
double precision :: Teb(ixG^T),pth(ixG^T),cg
logical :: patchw(ixG^T)
integer :: ix^D,idims,ixInt^L
!----------------------------------------------------------------------------
oktest = index(teststr,'specialbound')>=1
if (oktest) write(unitterm,*) ' === specialbound  (in ) : ', &
                'ixO^L : ',ixO^L

select case(iB)
case(3)
   do ix2=ixOmax2,ixOmin2,-1
     w(ixOmin1:ixOmax1,ix2,b1_:b3_)=(1.0d0/11.0d0)* &
          ( +2.0d0*w(ixOmin1:ixOmax1,ix2+3,b1_:b3_) &
            -9.0d0*w(ixOmin1:ixOmax1,ix2+2,b1_:b3_) &
           +18.0d0*w(ixOmin1:ixOmax1,ix2+1,b1_:b3_))
   enddo
   ! 2nd order CD for divB=0 to set normal B component better
   delydelx=dxlevel(2)/dxlevel(1)
   do ix2=ixOmax2,ixOmin2,-1
     do ix1=ixOmin1+1,ixOmax1-1
       w(ix1,ix2,b2_)=w(ix1,ix2+2,b2_) &
       +delydelx*(w(ix1+1,ix2+1,b1_)-w(ix1-1,ix2+1,b1_))
     enddo
   enddo
   !do ix2=ixOmax2,ixOmin2,-1
   !  w(ix2^%2ixO^S,b1_)=w(ixOmax2+1^%2ixO^S,b1_)
   !  w(ix2^%2ixO^S,b2_)=w(ixOmax2+1^%2ixO^S,b2_)
   !  w(ix2^%2ixO^S,b3_)=w(ixOmax2+1^%2ixO^S,b3_)
   !enddo
   w(ixO^S,v1_)=-w(ixOmin1:ixOmax1,ixOmax2+dixB:ixOmax2+1:-1,m1_)&
                 /w(ixOmin1:ixOmax1,ixOmax2+dixB:ixOmax2+1:-1,rho_)
   w(ixO^S,v2_)=-w(ixOmin1:ixOmax1,ixOmax2+dixB:ixOmax2+1:-1,m2_)&
                 /w(ixOmin1:ixOmax1,ixOmax2+dixB:ixOmax2+1:-1,rho_)
   w(ixO^S,v3_)=-w(ixOmin1:ixOmax1,ixOmax2+dixB:ixOmax2+1:-1,m3_)&
                 /w(ixOmin1:ixOmax1,ixOmax2+dixB:ixOmax2+1:-1,rho_)
   do ix2=ixOmin2,ixOmax2
     w(ixOmin1:ixOmax1,ix2,rho_)=rbc(ix2)
     w(ixOmin1:ixOmax1,ix2,p_)=pbc(ix2)
   enddo

   !!! obtain the thermal pressure in top layers and store in pth
   !ixInt^L=ixO^L;
   !ixIntmin2=ixOmax2+1;ixIntmax2=ixOmax2+dixB;
   !call getpthermal(w,x,ixG^L,ixInt^L,pth)
   !cg=dxlevel(2)*eqpar(grav2_)/two
   !!! fill pth, rho ghost layers according to gravity stratification 
   !do ix2=ixOmax2,ixOmin2,-1
   !  Teb(ixOmin1:ixOmax1,ix2+1)=pth(ixOmin1:ixOmax1,ix2+1)/w(ixOmin1:ixOmax1,ix2+1,rho_)
   !  pth(ixOmin1:ixOmax1,ix2)=(pth(ixOmin1:ixOmax1,ix2+1)- &
   !    cg*w(ixOmin1:ixOmax1,ix2+1,rho_))/(one+cg/Teb(ixOmin1:ixOmax1,ix2+1))
   !  w(ixOmin1:ixOmax1,ix2,rho_)=pth(ixOmin1:ixOmax1,ix2)/Teb(ixOmin1:ixOmax1,ix2+1)
   !end do
   !w(ixO^S,p_)=pth(ixO^S)
   patchw(ixO^S)=.false.
   call conserve(ixG^L,ixO^L,w,x,patchw)
case(4)
!! implementation of hydrostatic extrapolation at top boundary
   do ix2=ixOmin2,ixOmax2
     w(ixOmin1:ixOmax1,ix2,b1_:b3_)=(1.0d0/11.0d0)* &
          ( +2.0d0*w(ixOmin1:ixOmax1,ix2-3,b1_:b3_) &
            -9.0d0*w(ixOmin1:ixOmax1,ix2-2,b1_:b3_) &
           +18.0d0*w(ixOmin1:ixOmax1,ix2-1,b1_:b3_))
   enddo
   ! 2nd order CD for divB=0 to set normal B component better
   delydelx=dxlevel(2)/dxlevel(1)
   do ix2=ixOmin2,ixOmax2
     do ix1=ixOmin1+1,ixOmax1-1
       w(ix1,ix2,b2_)=w(ix1,ix2-2,b2_) &
       -delydelx*(w(ix1+1,ix2-1,b1_)-w(ix1-1,ix2-1,b1_))
     enddo
   enddo
   do ix2=ixOmin2,ixOmax2
    ! w(ix2^%2ixO^S,b1_)=w(ixOmin2-1^%2ixO^S,b1_)
    ! w(ix2^%2ixO^S,b2_)=w(ixOmin2-1^%2ixO^S,b2_)
    ! w(ix2^%2ixO^S,b3_)=w(ixOmin2-1^%2ixO^S,b3_)
     w(ix2^%2ixO^S,v1_)=w(ixOmin2-1^%2ixO^S,m1_)/w(ixOmin2-1^%2ixO^S,rho_)
     w(ix2^%2ixO^S,v2_)=w(ixOmin2-1^%2ixO^S,m2_)/w(ixOmin2-1^%2ixO^S,rho_)
     w(ix2^%2ixO^S,v3_)=w(ixOmin2-1^%2ixO^S,m3_)/w(ixOmin2-1^%2ixO^S,rho_)
   enddo
   !! obtain the thermal pressure in top layers and store in pth
   ixInt^L=ixO^L;
   ixIntmin2=ixOmin2-dixB;ixIntmax2=ixOmin2-1;
   call getpthermal(w,x,ixG^L,ixInt^L,pth)
   cg=dxlevel(2)*eqpar(grav2_)/two
   !! fill pth, rho ghost layers according to gravity stratification 
   do ix2=ixOmin2,ixOmax2
     do ix1=ixOmin1,ixOmax1
       Teb(ix1,ix2-1)=pth(ix1,ix2-1)/w(ix1,ix2-1,rho_)
       Teb(ix1,ix2)=Teb(ix1,ix2-1)
       pth(ix1,ix2)=(pth(ix1,ix2-1)+cg*w(ix1,ix2-1,rho_))/(one-cg/Teb(ix1,ix2))
       if(pth(ix1,ix2)<minp) pth(ix1,ix2)=pth(ix1,ix2-1)
       w(ix1,ix2,rho_)=pth(ix1,ix2)/Teb(ix1,ix2)
     end do
   end do
   w(ixO^S,p_)=pth(ixO^S)
   patchw(ixO^S)=.false.
   call conserve(ixG^L,ixO^L,w,x,patchw)
case default
   call mpistop("Special boundary is not defined for this region")
end select
end subroutine specialbound_usr
!=============================================================================
subroutine specialsource(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

! Calculate w(iw)=w(iw)+qdt*SOURCE[wCT,qtC,x] within ixO for all indices
! iw=iwmin...iwmax.  wCT is at time qCT

use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L, iw^LIM
double precision, intent(in) :: qdt, qtC, qt, x(ixI^S,1:ndim)
double precision, intent(inout) :: wCT(ixI^S,1:nw), w(ixI^S,1:nw)
!-----------------------------------------------------------------------------

call addsource_grav(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

end subroutine specialsource
!=============================================================================
subroutine getdt_special(w,ixG^L,ix^L,dtnew,dx^D,x)

! Limit "dt" further if necessary, e.g. due to the special source terms.
! The getdt_courant (CFL condition) and the getdt subroutine in the AMRVACPHYS
! module have already been called.

use mod_global_parameters

integer, intent(in) :: ixG^L, ix^L
double precision, intent(in) :: dx^D, x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw), dtnew
!-----------------------------------------------------------------------------

call getdt_grav(w,ixG^L,ix^L,dtnew,dx^D,x)

end subroutine getdt_special
!=============================================================================
subroutine specialeta(w,ixI^L,ix^L,idirmin,x,current,eta)

! Set the "eta" array for resistive MHD based on w or the
! "current" variable which has components between idirmin and 3.

use mod_global_parameters

integer, intent(in) :: ixI^L, ix^L, idirmin
double precision, intent(in) :: w(ixI^S,nw), x(ixI^S,1:ndim)

double precision :: current(ixG^T,7-2*ndir:3), eta(ixG^T)
!-----------------------------------------------------------------------------

!  eta(ix^S)=...

call mpistop("specialeta is not defined")

end subroutine specialeta
!=============================================================================
subroutine specialrefine_grid(igrid,level,ixG^L,ix^L,qt,w,x,refine,coarsen)

! Enforce additional refinement or coarsening
! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.

use mod_global_parameters

integer, intent(in) :: igrid, level, ixG^L, ix^L
double precision, intent(in) :: qt, w(ixG^S,1:nw), x(ixG^S,1:ndim)
integer, intent(inout) :: refine, coarsen
!-----------------------------------------------------------------------------


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
oktest = index(teststr,'printlog')>=1

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
double precision                   :: w(ixI^S,nw+nwauxio),tmp(ixG^T)
double precision                   :: normconv(0:nw+nwauxio)

double precision                   :: qtC,qdt
integer                            :: ix1
!-----------------------------------------------------------------------------
if(saveprim)then
   tmp(ixO^S)=w(ixO^S,p_)
 else
   call getpthermal(w,x,ixI^L,ixO^L,tmp)
endif
! output the temperature p/rho
w(ixO^S,nw+1)=tmp(ixO^S)/w(ixO^S,rho_)
!! output the plasma beta p*2/B**2
if(B0field)then
  w(ixO^S,nw+2)=tmp(ixO^S)*two/(^C&(w(ixO^S,b^C_)+myB0_cell%w(ixO^S,^C))**2+)
else
  w(ixO^S,nw+2)=tmp(ixO^S)*two/(^C&w(ixO^S,b^C_)**2+)
endif
! output divB1
call getdivb(w,ixI^L,ixO^L,tmp)
w(ixO^S,nw+3)=tmp(ixO^S)

end subroutine specialvar_output
!=============================================================================
subroutine specialvarnames_output

! newly added variables need to be concatenated with the wnames/primnames string

use mod_global_parameters
!-----------------------------------------------------------------------------
oktest = index(teststr,'printlog')>=1

primnames= TRIM(primnames)//' '//'Te beta divb'
wnames=TRIM(wnames)//' '//'Te beta divb'

end subroutine specialvarnames_output
!=============================================================================
subroutine specialset_B0(ixI^L,ixO^L,x,wB0)

! Here one can add a steady (time-independent) potential background field

use mod_global_parameters

integer, intent(in)           :: ixI^L,ixO^L
double precision, intent(in)  :: x(ixI^S,1:ndim)
double precision, intent(inout) :: wB0(ixI^S,1:ndir)
!-----------------------------------------------------------------------------

wB0(ixO^S,1)=wB0(ixO^S,1)+0.2d0*Busr
wB0(ixO^S,2)=wB0(ixO^S,2)+zero
wB0(ixO^S,3)=wB0(ixO^S,3)+Busr

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

end subroutine bc_int
!=============================================================================
function integral_grid(ixI^L,ixO^L,w,x,dvolume,iw,patchwi)

use mod_global_parameters

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
 case(m2_)
  ! kinetic energy integration
  {do ix^DB=ixOmin^DB,ixOmax^DB\}
     integral_grid=integral_grid+half*(^C&w(ix^D,m0_+^C)**2+)/w(ix^D,rho_)*&
                   dvolume(ix^D)
  {end do\}
 case(e_)
  ! total energy integration
  {do ix^DB=ixOmin^DB,ixOmax^DB\}
     integral_grid=integral_grid+w(ix^D,e_)*dvolume(ix^D)
  {end do\}
 case(9)
  ! magnetic energy integration
  {do ix^DB=ixOmin^DB,ixOmax^DB\}
     integral_grid=integral_grid+half*(^C&w(ix^D,b0_+^C)**2+)*dvolume(ix^D)
  {end do\}
end select
return
end function integral_grid
!=============================================================================
subroutine spacial_integral_w
! integration w in a given volume, output *.int ASCII file
use mod_global_parameters

double precision :: dvolume(ixG^T), timep4
double precision, allocatable :: integral_ipe(:), integral_w(:)
double precision, external :: integral_grid
integer          :: iigrid,igrid,iw,amode,int_fh,status(MPI_STATUS_SIZE),ni
character(len=100):: filename,purpose
character(len=1024) :: line
logical          :: patchwi(ixG^T),alive
!-----------------------------------------------------------------------------
purpose='Energy'
ni=8
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
  case('Energy')
  ! total kinetic energy
  integral_ipe(1)=integral_ipe(1)+integral_grid(ixG^LL,ixM^LL,pw(igrid)%w,px(igrid)%x,&
    dvolume,m2_,patchwi)
  ! total magnetic energy
  integral_ipe(2)=integral_ipe(2)+integral_grid(ixG^LL,ixM^LL,pw(igrid)%w,px(igrid)%x,&
    dvolume,9,patchwi)
  ! total energy
  integral_ipe(3)=integral_ipe(3)+integral_grid(ixG^LL,ixM^LL,pw(igrid)%w,px(igrid)%x,&
    dvolume,e_,patchwi)
  end select
end do
call MPI_ALLREDUCE(integral_ipe,integral_w,ni,MPI_DOUBLE_PRECISION,&
                     MPI_SUM,icomm,ierrmpi)
timep4=t**4
if(mype==0) then
  write(filename,"(a,a,a)") TRIM(filenamelog),TRIM(purpose),".int"
  inquire(file=filename,exist=alive)
  if(alive) then
    open(unit=11,file=filename,form='formatted',status='old',access='append')
  else
    open(unit=11,file=filename,form='formatted',status='new')
  endif
  select case(purpose)
   case('Energy')
    ! time, time**4, kinetic, magnetic, total energy
    write(11,'(5(es12.4))') t,timep4,integral_w(1),integral_w(2),integral_w(3)
  endselect
  close(11)
endif
if(mype==0) print*,'t=',t,'integral_w=',integral_w
deallocate(integral_ipe,integral_w)
end subroutine spacial_integral_w
!=============================================================================
subroutine userspecialconvert(qunitconvert)

use mod_global_parameters
integer, intent(in) :: qunitconvert
character(len=20):: userconvert_type
!-----------------------------------------------------------------------------
call spacial_integral_w
end subroutine userspecialconvert
!=============================================================================
subroutine init_particle_integrator()

use mod_particles
use mod_global_parameters
!-----------------------------------------------------------------------------

itmax_particles = 10000000
tmax_particles  = 100.0d0 
ditsave_particles = 1
! save all particles in one file every dtsave_ensemble
dtsave_ensemble   = dtsave(2) * UNIT_LENGTH/UNIT_VELOCITY

losses = .false.

end subroutine init_particle_integrator
!=============================================================================
subroutine init_particles()
! initialise the particles

use constants
use mod_particles
use mod_gridvars
use mod_global_parameters

double precision, dimension(ndir)    :: x
double precision                     :: v(1:ndir)
double precision                     :: lfac
integer                              :: igrid_particle, ipe_particle
integer                              :: Npart,npart_x,npart_y
integer                              :: ipart, idir
logical, allocatable                 :: follow(:)
!-----------------------------------------------------------------------------

call init_gridvars

x(:)=0.0d0

npart_x=20
npart_y=20
Npart=npart_x*npart_y
allocate(follow(Npart))
follow=.false.
do ipart=1,Npart
  if((ipart-1)/npart_x==11) then
    follow(ipart)=.true.
  endif
end do
do while (nparticles .lt. Npart)

   nparticles=nparticles+1
   x(1)=xprobmin1+mod(dble(nparticles)-0.5d0,dble(npart_x))*(xprobmax1-xprobmin1)/dble(npart_x)
   x(2)=xprobmin2+(dble((nparticles-1)/npart_x)+0.5d0)*(xprobmax2-xprobmin2)/dble(npart_y)
   if((nparticles-1)/npart_x==11) follow(nparticles)=.true.
   call find_particle_ipe(x,igrid_particle,ipe_particle)

   particle(nparticles)%igrid  = igrid_particle
   particle(nparticles)%ipe    = ipe_particle

   if (ipe_particle == mype) then

      call push_particle_into_particles_on_mype(nparticles)
      allocate(particle(nparticles)%self)
      particle(nparticles)%self%follow = follow(nparticles)
      particle(nparticles)%self%index  = nparticles
      particle(nparticles)%self%q      = - CONST_e
      particle(nparticles)%self%m      =   CONST_me

      particle(nparticles)%self%t      = 0.0d0
      particle(nparticles)%self%dt     = 0.0d0

      particle(nparticles)%self%x(:) = x

{#IFDEF PARTICLES_ADVECT
      {^C&call interpolate_var(igrid_particle,ixG^LL,ixM^LL,pw(igrid_particle)%w(ixG^T,v^C_),px(igrid_particle)%x(ixG^T,1:ndim),x,v(^C))\}
      {^C&particle(nparticles)%self%u(^C) = v(^C) * UNIT_VELOCITY\}
}
      particle(nparticles)%self%u(:) = 0.d0
      particle(nparticles)%self%payload(1:npayload) = 0.0d0

   end if

end do


end subroutine init_particles
!=============================================================================
logical function user_destroy(myparticle)

use mod_particles, only: particle_node

type(particle_node), intent(in)          :: myparticle
!-----------------------------------------------------------------------------
user_destroy = .false.

end function user_destroy
!=============================================================================
{#IFDEF FCT
subroutine initvecpot_usr(ixI^L, ixC^L, xC, A)

! initialize the vectorpotential on the corners
! used by b_from_vectorpotential()


use mod_global_parameters

integer, intent(in)                :: ixI^L, ixC^L
double precision, intent(in)       :: xC(ixI^S,1:ndim)
double precision, intent(out)      :: A(ixI^S,1:ndir)

!double precision                   :: r(ixG^T)

!-----------------------------------------------------------------------------

!r(ixC^S)=sqrt(xC(ixC^S,1)**2 + xC(ixC^S,2)**2 )

A(ixC^S,1:ndir) = zero
!where (r(ixC^S) .lt. eqpar(rm_))
!   A(ixC^S,3) = - half * eqpar(bm_) / eqpar(rm_) * r(ixC^S)**2
!elsewhere (r(ixC^S) .lt. eqpar(rj_))
!   A(ixC^S,3) = - half * eqpar(bm_) * eqpar(rm_) &
!        - eqpar(bm_) * eqpar(rm_) * log(r(ixC^S)/eqpar(rm_))
!elsewhere
!   A(ixC^S,3) = - half * eqpar(bm_) * eqpar(rm_) &
!        - eqpar(bm_) * eqpar(rm_) * log(eqpar(rj_)/eqpar(rm_))
!end where


end subroutine initvecpot_usr
!=============================================================================
}
!=============================================================================
! amrvacusr.t.promRTideal
!=============================================================================
