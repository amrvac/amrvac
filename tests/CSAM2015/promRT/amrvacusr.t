!#############################################################################
! module amrvacusr - promRTideal

! this setup to simulate Rayleigh-Taylor dynamics in a solar prominence
! It can be used in 2.5D or 3D, and assumes ideal MHD (with external gravity)
! One can use a tracer, and exploit GLM (or not)
! Related literature is in: 
!   `Solar prominences: "double, double ... boil and bubble"', 
!     R. Keppens, X. Cia, & O. Porth, 2015, ApJ Letters 806, L13 (7pp)

INCLUDE:amrvacnul/specialimpl.t
INCLUDE:amrvacnul/usrflags.t
!=============================================================================
module usr_bc
implicit none
double precision, allocatable, save :: pbc(:),rbc(:)
end module usr_bc
!=============================================================================
subroutine initglobaldata_usr

use mod_global_parameters
logical, save :: firstinitglobal=.true.

integer,dimension(:),allocatable:: seed
integer::  seed_size,ix
real:: randphaseP1(1:1000)
!-----------------------------------------------------------------------------
! CGS unit
k_B=1.3806d-16    ! erg*K^-1
miu0=4.d0*dpi     ! Gauss^2 cm^2 dyne^-1
Lunit=1.d9        ! cm = 10 Mm
UNIT_LENGTH=Lunit
Teunit=1.d6       ! K
nHunit=1.d9       ! cm^-3
mHunit=1.67262d-24 ! g
runit=1.4d0*mHunit*nHunit     ! 2.341668000000000E-015 g*cm^-3
UNIT_DENSITY=runit
punit=2.3d0*nHunit*k_B*Teunit ! 0.317538000000000 erg*cm^-3
Bunit=dsqrt(miu0*punit)       ! 1.99757357615242 Gauss
vunit=Bunit/dsqrt(miu0*runit) ! 1.16448846777562E007 cm/s = 116.45 km/s
UNIT_VELOCITY=vunit
tunit=Lunit/vunit     ! 85.8746159942810 s
heatunit=punit/tunit  ! 3.697693390805347E-003 erg*cm^-3/s

! units for convert
if(iprob==-1) then
  normvar(0) = one
else
  normvar(0) = UNIT_LENGTH
endif
normvar(rho_) = UNIT_DENSITY
{^C&normvar(v^C_)   = UNIT_VELOCITY \}
normvar(pp_)     = UNIT_VELOCITY**2 * UNIT_DENSITY
{^C&normvar(b^C_)   = dsqrt(4.0d0*dpi*normvar(pp_)) \}
normt = UNIT_LENGTH/UNIT_VELOCITY

eqpar(grav1_)=0.d0
eqpar(grav2_)=-2.74d4*Lunit/vunit**2 ! solar surface gravity 2.74d2 m*s^-2
{^IFTHREED
eqpar(grav3_)=0.d0
}
eqpar(gamma_)=5.0d0/3.0d0
eqpar(eta_)=0.0d0

! gzone gives the distance from the lower boundary xprobmin2
! to the so-called photosphere where density is fixed to rho0 further on
! this has to be wider than the ghost layer zone for the coarsest mesh....
gzone=0.3d0 ! our grid starts above photosphere
dr=(2.d0*gzone+xprobmax2-xprobmin2)/dble(jmax)
! Solar radius
SRadius=6.961d10/Lunit
! base density and temperature
rho0=1.d17/nHunit
Tch=8.d3/Teunit
Tco=1.8d6/Teunit
Tpromin=6.0d3/Teunit
Tpromax=1.4d4/Teunit
htra1=0.2d0
htra2=1.125d0
htra3=2.0d0
ybot=1.75d0
ytop=2.0d0
bsca=1.5d0
pwidth=0.5d0 ! width in Lunit

   eqpar(nxmodes_)=50
   eqpar(BB1_)=0.1d0
   eqpar(BB3_)=4.d0
   eqpar(eps_)=0.05d0

{#IFDEF GLM
eqpar(Cr_)=-0.2d0
}

if(firstinitglobal) then
  call inithdstatic
  firstinitglobal=.false.
endif
        
randphase(1:1000)=zero
if(eqpar(nxmodes_)>1000) call mpistop('too many modes, edit amrvacusrpar')

if(mype==0)then
    call random_seed(SIZE=seed_size)
    allocate(seed(seed_size))
    call random_seed(GET=seed(1:seed_size))
    call random_number(randphaseP1(1:nint(eqpar(nxmodes_))))
    randphase(1:nint(eqpar(nxmodes_)))=-dpi+two*dpi*dble(randphaseP1(1:nint(eqpar(nxmodes_))))
endif
call MPI_BARRIER(icomm,ierrmpi)
if(npe>1)then
     call MPI_BCAST(randphase,1000,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
endif

if(mype==0)then
    print *,'number of modes=',eqpar(nxmodes_)
    open(123,file='phaseinfo',form='formatted')
    write(123,*) nint(eqpar(nxmodes_))
    do ix=1,nint(eqpar(nxmodes_))
        write(123,"(i4,1es12.4)") ix,randphase(ix)
    enddo
    close(123)
endif

end subroutine initglobaldata_usr
!=============================================================================
subroutine inithdstatic
!! initialize the table in a vertical line through the global domain
use usr_bc
use mod_global_parameters

integer :: j,na,ibc,ix
double precision:: Ta(jmax),gg(jmax),Taext(jmax)
double precision:: res,rhob,pb,wtra1,wtra2,wtra3,b3
!----------------------------------------------------------------------------
wtra1=0.01d0
wtra2=0.01d0
wtra3=0.01d0
do j=1,jmax
   ya(j)=(dble(j)-0.5d0)*dr-gzone+xprobmin2
   if(ya(j)<=0.4d0) then
     Ta(j)=Tch+0.5d0*(Tco-Tch)*(tanh((ya(j)-htra1)/wtra1)+1.d0)
     Taext(j)=Tch+0.5d0*(Tco-Tch)*(tanh((ya(j)-htra1)/wtra1)+1.d0)
   endif
   if(ya(j)>0.4d0.and.ya(j)<ybot) then
     Ta(j)=Tco-0.5d0*(Tco-Tpromin)*(tanh((ya(j)-htra2)/wtra2)+1.d0)
     Taext(j)=Tco
   endif
   if(ya(j)>=ybot.and.ya(j)<ytop) then
     Ta(j)=Tpromin+(Tpromax-Tpromin)*(ya(j)-ybot)/(ytop-ybot)
     Taext(j)=Tco
   endif
   if(ya(j)>ytop) then
     Ta(j)=Tpromax+0.5d0*(Tco-Tpromax)*(tanh((ya(j)-htra3)/wtra3)+1.d0)
     Taext(j)=Tco
   endif
   gg(j)=eqpar(grav2_)*SRadius**2/(SRadius+ya(j))**2
   !!gg(j)=eqpar(grav2_)
enddo
!! solution of hydrostatic equation 
ra(1)=rho0
raext(1)=rho0
pa(1)=rho0*Tch
paext(1)=rho0*Tch
do j=2,jmax
  if(ya(j)<htra2) then
    pa(j)=(pa(j-1)+dr*(gg(j)+gg(j-1))*ra(j-1)/4.d0)/(one-dr*(gg(j)+gg(j-1))/&
          Ta(j)/4.d0)
    paext(j)=(paext(j-1)+dr*(gg(j)+gg(j-1))*raext(j-1)/4.d0)/(one-dr*(gg(j)+gg(j-1))/&
          Taext(j)/4.d0)
  else
    if(ya(j)<ybot) then
    b3=(eqpar(BB3_)*dexp(-(ya(j-1)-htra2)/bsca))**2+(eqpar(BB3_)*dexp(-(ya(j)-htra2)/bsca))**2
    else
    b3=zero
    endif
    pa(j)=(pa(j-1)+dr*(gg(j)+gg(j-1))*ra(j-1)/4.d0+dr/bsca*b3/2.d0)/(one-dr*(gg(j)+gg(j-1))/&
          Ta(j)/4.d0)
    paext(j)=(paext(j-1)+dr*(gg(j)+gg(j-1))*raext(j-1)/4.d0+dr/bsca*b3/2.d0)/(one-dr*(gg(j)+gg(j-1))/&
          Taext(j)/4.d0)
  endif
  ra(j)=pa(j)/Ta(j)
  raext(j)=paext(j)/Taext(j)
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
 open(123,file='pruns',form='formatted')
 write(123,*) jmax
 do ix=1,jmax
    write(123,"(i7,7es12.4)") ix,ya(ix),pa(ix),paext(ix),ra(ix),raext(ix),Ta(ix),Taext(ix)
 enddo
 close(123)
endif

end subroutine inithdstatic
!=============================================================================
subroutine initonegrid_usr(ixG^L,ix^L,w,x)

! initialize one grid

use mod_global_parameters

integer, intent(in) :: ixG^L, ix^L
double precision, intent(in) :: x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)

double precision:: psi(ixG^T),tmp(ixG^T)
double precision:: res,sigma,lxsize,sigma3
integer :: ix^D,na,idims,imode
logical, save:: first=.true.
logical :: patchw(ixG^S)
!----------------------------------------------------------------------------

if (first) then
   if (mype==0) then
      print *,'2.5D or 3D MHD Rayleigh-Taylor instability in prominence'
   end if
   first=.false.
end if
sigma=0.02d0
sigma3=0.1d0

w(ix^S,v1_)=0.d0
w(ix^S,v2_)=0.d0
w(ix^S,v3_)=0.d0
! now add the incompressible perturbations
psi(ixG^T)=zero
lxsize=(xprobmax1-xprobmin1)
do imode=1,nint(eqpar(nxmodes_))
   psi(ixG^T)=psi(ixG^T) &
       +eqpar(eps_)*(dcos(two*dpi*dble(imode)*(x(ixG^T,1)/lxsize)+randphase(imode)) &
                    /dble(imode)) &
                   *dexp(-((x(ixG^T,2)-htra2)/sigma)**2) 
enddo
! compute dpsi/dy
idims=2
select case(typegrad)
    case("central")
     call gradient(psi,ixG^L,ix^L,idims,tmp)
    case("limited")
     call gradientS(psi,ixG^L,ix^L,idims,tmp)
end select
w(ix^S,v1_)=w(ix^S,v1_)-tmp(ix^S){^IFTHREED *dexp(-(x(ixG^T,3)/sigma3)**2)}
! compute dpsi/dx
idims=1
select case(typegrad)
     case("central")
      call gradient(psi,ixG^L,ix^L,idims,tmp)
     case("limited")
      call gradientS(psi,ixG^L,ix^L,idims,tmp)
end select
w(ix^S,v2_)=w(ix^S,v2_)+tmp(ix^S){^IFTHREED *dexp(-(x(ixG^T,3)/sigma3)**2)}


{do ix^DB=ixmin^DB,ixmax^DB\}
   na=floor((x(ix^D,2)-xprobmin2+gzone)/dr+0.5d0)
   res=x(ix^D,2)-xprobmin2+gzone-(dble(na)-0.5d0)*dr
   {^IFTHREED
   if (dabs(x(ix^D,3))<pwidth*half)then
     w(ix^D,rho_)=ra(na)+(one-cos(dpi*res/dr))/two*(ra(na+1)-ra(na))
     w(ix^D,p_)  =pa(na)+(one-cos(dpi*res/dr))/two*(pa(na+1)-pa(na))
   else
     w(ix^D,rho_)=raext(na)+(one-cos(dpi*res/dr))/two*(raext(na+1)-raext(na))
     w(ix^D,p_)  =paext(na)+(one-cos(dpi*res/dr))/two*(paext(na+1)-paext(na))
   endif
   }
   {^IFTWOD
   w(ix^D,rho_)=ra(na)+(one-cos(dpi*res/dr))/two*(ra(na+1)-ra(na))
   w(ix^D,p_)  =pa(na)+(one-cos(dpi*res/dr))/two*(pa(na+1)-pa(na))
   }
{end do\}

w(ix^S,b1_)  =eqpar(BB1_)
w(ix^S,b2_)  =zero
w(ix^S,b3_)  =eqpar(BB3_)
where(x(ix^S,2)>htra2.and.x(ix^S,2)<ybot)
   w(ix^S,b3_)  =eqpar(BB3_)*dexp(-(x(ix^S,2)-htra2)/bsca)
endwhere
where(x(ix^S,2)>=ybot)
   w(ix^S,b3_)  =eqpar(BB3_)*dexp(-(ybot-htra2)/bsca)
endwhere

{#IFDEF TRACER
w(ix^S,tr1_)=zero
where(x(ix^S,2)>htra2.and.x(ix^S,2)<htra3{^IFTHREED .and.(dabs(x(ix^S,3))<pwidth*half)})
 w(ix^S,tr1_)=one
endwhere
where(x(ix^S,2)<htra1)
 w(ix^S,tr1_)=-one
endwhere
}

{#IFDEF GLM
w(ixG^S,psi_)=0.d0
}

patchw(ixG^S)=.false.
call conserve(ixG^L,ix^L,w,x,patchw)

end subroutine initonegrid_usr
!=============================================================================
subroutine specialbound_usr(qt,ixG^L,ixO^L,iw,iB,w,x)
use usr_bc
! special boundary types, user defined

use mod_global_parameters

integer, intent(in) :: ixO^L, iw, iB, ixG^L
double precision, intent(in) :: qt, x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)

double precision :: dx^D,delydelx,gjjm1,delydelz
double precision :: Teb(ixG^T),pth(ixG^T),cg
logical :: patchw(ixG^T)
integer :: ix^D,idims,ixInt^L
!----------------------------------------------------------------------------

select case(iB)
case(3)
   w(ixO^S,b1_)=eqpar(BB1_)
   w(ixO^S,b3_)=eqpar(BB3_)
   w(ixO^S,b2_)=zero
{#IFDEF TRACER
   {^IFTWOD
   do ix2=ixOmax2,ixOmin2,-1
     w(ixOmin1:ixOmax1,ix2,tr1_)=w(ixOmin1:ixOmax1,ixOmax2+1,Dtr1_)/w(ixOmin1:ixOmax1,ixOmax2+1,rho_)
   enddo
   }
   {^IFTHREED
   do ix2=ixOmax2,ixOmin2,-1
     w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,tr1_)= &
     w(ixOmin1:ixOmax1,ixOmax2+1,ixOmin3:ixOmax3,Dtr1_)/w(ixOmin1:ixOmax1,ixOmax2+1,ixOmin3:ixOmax3,rho_)
   enddo
   }
}
{#IFDEF GLM
   {^IFTWOD
   do ix2=ixOmax2,ixOmin2,-1
     w(ixOmin1:ixOmax1,ix2,psi_)=w(ixOmin1:ixOmax1,ixOmax2+1,psi_)
   enddo
   }
   {^IFTHREED
   do ix2=ixOmax2,ixOmin2,-1
     w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,psi_)=w(ixOmin1:ixOmax1,ixOmax2+1,ixOmin3:ixOmax3,psi_)
   enddo
   }
}
   {^IFTWOD
   ! 2nd order CD for divB=0 to set normal B component better
   delydelx=dxlevel(2)/dxlevel(1)
   do ix2=ixOmax2,ixOmin2,-1
     do ix1=ixOmin1+1,ixOmax1-1
       w(ix1,ix2,b2_)=w(ix1,ix2+2,b2_) &
       +delydelx*(w(ix1+1,ix2+1,b1_)-w(ix1-1,ix2+1,b1_))
     enddo
   enddo
   }
   {^IFTHREED
   delydelx=dxlevel(2)/dxlevel(1)
   delydelz=dxlevel(2)/dxlevel(3)
   do ix2=ixOmax2,ixOmin2,-1
      do ix1=ixOmin1+1,ixOmax1-1
      do ix3=ixOmin3+1,ixOmax3-1
       w(ix1,ix2,ix3,b2_)=w(ix1,ix2+2,ix3,b2_)+delydelx*(w(ix1+1,ix2+1,ix3,b1_)-w(ix1-1,ix2+1,ix3,b1_)) &
                                              +delydelz*(w(ix1,ix2+1,ix3+1,b3_)-w(ix1,ix2+1,ix3-1,b3_))
      enddo
      enddo
   enddo
   }

   {^IFTWOD
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
   }
   {^IFTHREED
   w(ixO^S,v1_)=-w(ixOmin1:ixOmax1,ixOmax2+dixB:ixOmax2+1:-1,ixOmin3:ixOmax3,m1_)&
                 /w(ixOmin1:ixOmax1,ixOmax2+dixB:ixOmax2+1:-1,ixOmin3:ixOmax3,rho_)
   w(ixO^S,v2_)=-w(ixOmin1:ixOmax1,ixOmax2+dixB:ixOmax2+1:-1,ixOmin3:ixOmax3,m2_)&
                 /w(ixOmin1:ixOmax1,ixOmax2+dixB:ixOmax2+1:-1,ixOmin3:ixOmax3,rho_)
   w(ixO^S,v3_)=-w(ixOmin1:ixOmax1,ixOmax2+dixB:ixOmax2+1:-1,ixOmin3:ixOmax3,m3_)&
                 /w(ixOmin1:ixOmax1,ixOmax2+dixB:ixOmax2+1:-1,ixOmin3:ixOmax3,rho_)
   do ix2=ixOmin2,ixOmax2
     w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,rho_)=rbc(ix2)
     w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,p_)=pbc(ix2)
   enddo
   }

   patchw(ixO^S)=.false.
   call conserve(ixG^L,ixO^L,w,x,patchw)
case(4)
!! implementation of hydrostatic extrapolation at top boundary
   {^IFTWOD
   do ix2=ixOmin2,ixOmax2
     w(ixOmin1:ixOmax1,ix2,b1_:b3_)=(1.0d0/11.0d0)* &
          ( +2.0d0*w(ixOmin1:ixOmax1,ix2-3,b1_:b3_) &
            -9.0d0*w(ixOmin1:ixOmax1,ix2-2,b1_:b3_) &
           +18.0d0*w(ixOmin1:ixOmax1,ix2-1,b1_:b3_))
   enddo
   }
   {^IFTHREED
   do ix2=ixOmin2,ixOmax2
     w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,b1_:b3_)=(1.0d0/11.0d0)* &
          ( +2.0d0*w(ixOmin1:ixOmax1,ix2-3,ixOmin3:ixOmax3,b1_:b3_) &
            -9.0d0*w(ixOmin1:ixOmax1,ix2-2,ixOmin3:ixOmax3,b1_:b3_) &
           +18.0d0*w(ixOmin1:ixOmax1,ix2-1,ixOmin3:ixOmax3,b1_:b3_))
   enddo
   }
   ! 2nd order CD for divB=0 to set normal B component better
   {^IFTWOD
   delydelx=dxlevel(2)/dxlevel(1)
   do ix2=ixOmin2,ixOmax2
     do ix1=ixOmin1+1,ixOmax1-1
       w(ix1,ix2,b2_)=w(ix1,ix2-2,b2_) &
       -delydelx*(w(ix1+1,ix2-1,b1_)-w(ix1-1,ix2-1,b1_))
     enddo
   enddo
   }
   {^IFTHREED
   delydelx=dxlevel(2)/dxlevel(1)
   delydelz=dxlevel(2)/dxlevel(3)
   do ix2=ixOmin2,ixOmax2,+1
      do ix1=ixOmin1+1,ixOmax1-1
      do ix3=ixOmin3+1,ixOmax3-1
        w(ix1,ix2,ix3,b2_)=w(ix1,ix2-2,ix3,b2_)-delydelx*(w(ix1+1,ix2-1,ix3,b1_)-w(ix1-1,ix2-1,ix3,b1_)) &
                                               -delydelz*(w(ix1,ix2-1,ix3+1,b3_)-w(ix1,ix2-1,ix3-1,b3_))
      enddo
      enddo
   enddo
   }

   do ix2=ixOmin2,ixOmax2
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
   {^IFTWOD
   do ix2=ixOmin2,ixOmax2
     do ix1=ixOmin1,ixOmax1
       Teb(ix1,ix2-1)=pth(ix1,ix2-1)/w(ix1,ix2-1,rho_)
       Teb(ix1,ix2)=Teb(ix1,ix2-1)
       gjjm1=half*(SRadius**2/(SRadius+x(ix1,ix2,2))**2+SRadius**2/(SRadius+x(ix1,ix2-1,2))**2)
       pth(ix1,ix2)=(pth(ix1,ix2-1)+cg*gjjm1*w(ix1,ix2-1,rho_))/(one-cg*gjjm1/Teb(ix1,ix2))
       if(pth(ix1,ix2)<minp) pth(ix1,ix2)=pth(ix1,ix2-1)
       w(ix1,ix2,rho_)=pth(ix1,ix2)/Teb(ix1,ix2)
     end do
   end do
   }
   {^IFTHREED
   do ix2=ixOmin2,ixOmax2
     do ix1=ixOmin1,ixOmax1
     do ix3=ixOmin3,ixOmax3
       Teb(ix1,ix2-1,ix3)=pth(ix1,ix2-1,ix3)/w(ix1,ix2-1,ix3,rho_)
       Teb(ix1,ix2,ix3)=Teb(ix1,ix2-1,ix3)
       gjjm1=half*(SRadius**2/(SRadius+x(ix1,ix2,ix3,2))**2+SRadius**2/(SRadius+x(ix1,ix2-1,ix3,2))**2)
       pth(ix1,ix2,ix3)=(pth(ix1,ix2-1,ix3)+cg*gjjm1*w(ix1,ix2-1,ix3,rho_))/(one-cg*gjjm1/Teb(ix1,ix2,ix3))
       if(pth(ix1,ix2,ix3)<minp) pth(ix1,ix2,ix3)=pth(ix1,ix2-1,ix3)
       w(ix1,ix2,ix3,rho_)=pth(ix1,ix2,ix3)/Teb(ix1,ix2,ix3)
     end do
     end do
   end do
   }
   w(ixO^S,p_)=pth(ixO^S)
{#IFDEF TRACER
   {^IFTWOD
   do ix2=ixOmin2,ixOmax2
     w(ixOmin1:ixOmax1,ix2,tr1_)=w(ixOmin1:ixOmax1,ixOmin2-1,Dtr1_)/w(ixOmin1:ixOmax1,ixOmin2-1,rho_)
   enddo
   }
   {^IFTHREED
   do ix2=ixOmin2,ixOmax2
     w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,tr1_)= &
     w(ixOmin1:ixOmax1,ixOmin2-1,ixOmin3:ixOmax3,Dtr1_)/w(ixOmin1:ixOmax1,ixOmin2-1,ixOmin3:ixOmax3,rho_)
   enddo
   }
}
{#IFDEF GLM
   {^IFTWOD
   do ix2=ixOmin2,ixOmax2
     w(ixOmin1:ixOmax1,ix2,psi_)=w(ixOmin1:ixOmax1,ixOmin2-1,psi_)
   enddo
   }
   {^IFTHREED
   do ix2=ixOmin2,ixOmax2
     w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,psi_)=w(ixOmin1:ixOmax1,ixOmin2-1,ixOmin3:ixOmax3,psi_)
   enddo
   }
}
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

double precision :: bQgrid(ixG^T)
integer :: iw
!-----------------------------------------------------------------------------

call addsource_gravSA(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

end subroutine specialsource
!=============================================================================
subroutine getggrav(ggrid,ixI^L,ixO^L,x)

use mod_global_parameters

integer, intent(in) :: ixI^L,ixO^L
double precision, intent(in) :: x(ixI^S,1:ndim)
double precision, intent(out) :: ggrid(ixG^T)
!---------------------------------------------------------------------------
! calculate gravity
ggrid(ixO^S)=eqpar(grav2_)*(SRadius/(SRadius+x(ixO^S,2)))**2

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
subroutine getdt_special(w,ixG^L,ix^L,dtnew,dx^D,x)

! Limit "dt" further if necessary, e.g. due to the special source terms.
! The getdt_courant (CFL condition) and the getdt subroutine in the AMRVACPHYS
! module have already been called.

use mod_global_parameters

integer, intent(in) :: ixG^L, ix^L
double precision, intent(in) :: dx^D, x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw), dtnew
!-----------------------------------------------------------------------------
dtnew=bigdouble

call getdt_grav(w,ixG^L,ix^L,dtnew,dx^D,x)

end subroutine getdt_special
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
subroutine specialeta(w,ixI^L,ix^L,idirmin,x,current,eta)

! Set the "eta" array for resistive MHD based on w or the
! "current" variable which has components between idirmin and 3.

use mod_global_parameters

integer, intent(in) :: ixI^L, ix^L, idirmin
double precision, intent(in) :: w(ixI^S,nw), x(ixI^S,1:ndim)

double precision :: current(ixG^T,7-2*ndir:3), eta(ixG^T)
!-----------------------------------------------------------------------------

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

double precision :: wloc(ixG^T,1:nw)
!-----------------------------------------------------------------------------

wloc(ixI^S,1:nw)=w(ixI^S,1:nw)
if(saveprim)then
   tmp(ixO^S)=wloc(ixO^S,p_)
 else
   call getpthermal(wloc,x,ixI^L,ixO^L,tmp)
endif

if(iprob==-1)then
  w(ixO^S,nw+1)=wloc(ixO^S,rho_)
  w(ixO^S,nw+2)=tmp(ixO^S)/wloc(ixO^S,rho_)
else
  ! output the temperature p/rho
  w(ixO^S,nw+1)=tmp(ixO^S)/wloc(ixO^S,rho_)*Teunit
  ! output the plasma beta p*2/B**2
  w(ixO^S,nw+2)=tmp(ixO^S)*two/(^C&wloc(ixO^S,b^C_)**2+)
endif

end subroutine specialvar_output
!=============================================================================
subroutine specialvarnames_output

! newly added variables need to be concatenated with the wnames/primnames string

use mod_global_parameters
!-----------------------------------------------------------------------------

if(iprob==-1)then
  primnames= TRIM(primnames)//' '//'rho T'
  wnames=TRIM(wnames)//' '//'rho T'
else
  primnames= TRIM(primnames)//' '//'T beta'
  wnames=TRIM(wnames)//' '//'T beta'
endif

end subroutine specialvarnames_output
!=============================================================================
subroutine specialset_B0(ixI^L,ixO^L,x,wB0)

! Here one can add a steady (time-independent) potential background field

use mod_global_parameters

integer, intent(in)           :: ixI^L,ixO^L
double precision, intent(in)  :: x(ixI^S,1:ndim)
double precision, intent(inout) :: wB0(ixI^S,1:ndir)
!-----------------------------------------------------------------------------

call mpistop("special B0 undefined")

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
subroutine printlog_special

! printlog: calculates volume averaged mean values

use mod_global_parameters

logical :: fileopen
integer :: iigrid, igrid, level, nleafs_level(1:nlevelshi), iw, i
double precision :: wmean(1:nw), volume(1:nlevelshi), volprob, voltotal
double precision :: dvolume(ixG^T), volumeflat(1:nlevelshi)
double precision :: tmpw(ixG^T)
integer :: numlevels, imon, idirmin
integer, dimension(1:nlevelshi) :: isum_send, isum_recv
double precision, dimension(1:nw+2+nlevelshi) :: dsum_send, dsum_recv
double precision :: wmeanmore
double precision :: wmaxte,wminte,wmaxvel,wmaxte_mype,wminte_mype,wmaxvel_mype
double precision :: invgminone
character(len=80) :: filename
character(len=1024) :: line
logical, save :: opened=.false.
integer :: amode, status(MPI_STATUS_SIZE)
!-----------------------------------------------------------------------------
volume(1:mxnest)=zero
volumeflat(1:mxnest)=zero
wmean(1:nw)= zero
nleafs_level(1:mxnest)=0

wmaxte=zero
wmaxte_mype=zero
wmaxvel=zero
wmaxvel_mype=zero
wminte=1.d20
wminte_mype=1.d20
invgminone=one/(eqpar(gamma_)-one)
wmeanmore=zero

do iigrid=1,igridstail; igrid=igrids(iigrid);
   level=node(plevel_,igrid)
   nleafs_level(level)=nleafs_level(level)+1
   volumeflat(level)=volumeflat(level)+ &
          {(rnode(rpxmax^D_,igrid)-rnode(rpxmin^D_,igrid))|*}
   if (slab) then
      dvolume(ixM^T)={rnode(rpdx^D_,igrid)|*}
   else
      dvolume(ixM^T)=pgeo(igrid)%dvolume(ixM^T)
      volume(level)=volume(level)+sum(dvolume(ixM^T))
   end if
   ! set dxlevel for use in gradient evaluation
   ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
   ! just use array for velocity
   tmpw(ixM^T)=dsqrt((pw(igrid)%w(ixM^T,m1_)/pw(igrid)%w(ixM^T,rho_))**2 &
                    +(pw(igrid)%w(ixM^T,m2_)/pw(igrid)%w(ixM^T,rho_))**2 &
                    +(pw(igrid)%w(ixM^T,m3_)/pw(igrid)%w(ixM^T,rho_))**2)
   wmaxvel_mype=max(wmaxvel_mype,maxval(tmpw(ixM^T)))
{#IFDEF TRACER
   wmean(Dtr1_)=wmean(Dtr1_)+sum(dvolume(ixM^T)*pw(igrid)%w(ixM^T,Dtr1_)/pw(igrid)%w(ixM^T,rho_))
}
{#IFDEF GLM
   wmean(psi_)=wmean(psi_)+sum(dvolume(ixM^T)*pw(igrid)%w(ixM^T,psi_))
}
   wmean(rho_)=wmean(rho_)+sum(dvolume(ixM^T)*pw(igrid)%w(ixM^T,rho_))
   ! kinetic energy in x
   wmean(m1_)=wmean(m1_)+half*sum(dvolume(ixM^T)*(pw(igrid)%w(ixM^T,m1_)**2)/pw(igrid)%w(ixM^T,rho_))
   ! kinetic energy in y
   wmean(m2_)=wmean(m2_)+half*sum(dvolume(ixM^T)*(pw(igrid)%w(ixM^T,m2_)**2)/pw(igrid)%w(ixM^T,rho_))
   ! kinetic energy in z
   wmean(m3_)=wmean(m3_)+half*sum(dvolume(ixM^T)*(pw(igrid)%w(ixM^T,m3_)**2)/pw(igrid)%w(ixM^T,rho_))
   ! total energy
   wmean(e_)=wmean(e_)+sum(dvolume(ixM^T)*pw(igrid)%w(ixM^T,e_))
   ! magnetic energy
   wmean(b1_)=wmean(b1_)+half*sum(dvolume(ixM^T)* &
       (pw(igrid)%w(ixM^T,b1_)**2+pw(igrid)%w(ixM^T,b2_)**2+pw(igrid)%w(ixM^T,b3_)**2))
   ! magnetic energy in y
   wmean(b2_)=wmean(b2_)+half*sum(dvolume(ixM^T)*(pw(igrid)%w(ixM^T,b2_)**2))
   ! magnetic energy in z
   wmean(b3_)=wmean(b3_)+half*sum(dvolume(ixM^T)*(pw(igrid)%w(ixM^T,b3_)**2))
   call getpthermal(pw(igrid)%w,px(igrid)%x,ixG^LL,ixM^LL,tmpw)
   wmeanmore=wmeanmore+invgminone*sum(dvolume(ixM^T)*tmpw(ixM^T))
   ! maximal temperature
   wmaxte_mype=max(wmaxte_mype,maxval(tmpw(ixM^T)/pw(igrid)%w(ixM^T,rho_)))
   ! minimal temperature
   wminte_mype=min(wminte_mype,minval(tmpw(ixM^T)/pw(igrid)%w(ixM^T,rho_)))
end do
if (slab) volume(levmin:levmax)=volumeflat(levmin:levmax)

voltotal=sum(volume(levmin:levmax))

call MPI_REDUCE(wmaxte_mype,wmaxte,1,MPI_DOUBLE_PRECISION, &
                MPI_MAX,0,icomm,ierrmpi)
call MPI_REDUCE(wmaxvel_mype,wmaxvel,1,MPI_DOUBLE_PRECISION, &
                MPI_MAX,0,icomm,ierrmpi)

call MPI_REDUCE(wminte_mype,wminte,1,MPI_DOUBLE_PRECISION, &
                MPI_MIN,0,icomm,ierrmpi)

numlevels=levmax-levmin+1
dsum_send(1:nw)=wmean(1:nw)
dsum_send(nw+1)=voltotal
dsum_send(nw+2:nw+1+numlevels)=volumeflat(levmin:levmax)
dsum_send(nw+2+numlevels:nw+2+numlevels)=wmeanmore
call MPI_REDUCE(dsum_send,dsum_recv,nw+2+numlevels,MPI_DOUBLE_PRECISION, &
                MPI_SUM,0,icomm,ierrmpi)
isum_send(1:numlevels)=nleafs_level(levmin:levmax)
call MPI_REDUCE(isum_send,isum_recv,numlevels,MPI_INTEGER, &
                MPI_SUM,0,icomm,ierrmpi)

if (mype==0) then

   wmean(1:nw)=dsum_recv(1:nw)
   wmeanmore=dsum_recv(nw+2+numlevels)
   voltotal=dsum_recv(nw+1)
   volumeflat(levmin:levmax)=dsum_recv(nw+2:nw+1+numlevels)
   nleafs_level(levmin:levmax)=isum_recv(1:numlevels)

   wmean=wmean/voltotal
   wmeanmore=wmeanmore/voltotal

   ! determine coverage in coordinate space
   volprob={(xprobmax^D-xprobmin^D)|*}
   volumeflat(levmin:levmax)=volumeflat(levmin:levmax)/volprob

   if (.not.opened) then
      ! generate filename
      write(filename,"(a,a)") TRIM(filenamelog),".log"

      amode=ior(MPI_MODE_CREATE,MPI_MODE_WRONLY)
      amode=ior(amode,MPI_MODE_APPEND)
      call MPI_FILE_OPEN(MPI_COMM_SELF,filename,amode, &
                         MPI_INFO_NULL,log_fh,ierrmpi)
      opened=.true.
      call MPI_FILE_WRITE(log_fh,fileheadout,len_trim(fileheadout), &
                          MPI_CHARACTER,status,ierrmpi)
      !!call MPI_FILE_WRITE(log_fh,new_line('a'),1,MPI_CHARACTER,status,ierrmpi)
      call MPI_FILE_WRITE(log_fh,achar(10),1,MPI_CHARACTER,status,ierrmpi)

      i=len_trim(wnameslog)-1
      write(wnameslog(i+3:i+14),"(a11)") "d1 d2 d3 d4"
      i=i+17
      do level=1,mxnest
          i=i+3
          write(wnameslog(i:i+1),"(a,i1)") "c",level
      end do
      do level=1,mxnest
          i=i+3
          write(wnameslog(i:i+1),"(a,i1)") "n",level
      end do

      if(residmin>smalldouble) then
        write(line,'(a15,a79)')"it   t  dt res ",wnameslog
      else
        write(line,'(a15,a79)')"it   t   dt    ",wnameslog
      endif

      call MPI_FILE_WRITE(log_fh,line,len_trim(line),MPI_CHARACTER, &
                          status,ierrmpi)
   end if
   !!call MPI_FILE_WRITE(log_fh,new_line('a'),1,MPI_CHARACTER,status,ierrmpi)
   call MPI_FILE_WRITE(log_fh,achar(10),1,MPI_CHARACTER,status,ierrmpi)

   if(residmin>smalldouble) then
      write(line,'(i7,3(e13.5))')it,t,dt,residual
   else
      write(line,'(i7,2(e13.5))')it,t,dt
   endif

   call MPI_FILE_WRITE(log_fh,line,len_trim(line), &
                       MPI_CHARACTER,status,ierrmpi)
   do iw=1,nw
      write(line,'(e13.5)')wmean(iw)
      call MPI_FILE_WRITE(log_fh,line,len_trim(line), &
                          MPI_CHARACTER,status,ierrmpi)
   end do
   do imon=1,1
      write(line,'(e13.5)')wmeanmore
      call MPI_FILE_WRITE(log_fh,line,len_trim(line), &
                          MPI_CHARACTER,status,ierrmpi)
   end do
   write(line,'(e13.5)')wmaxte
   call MPI_FILE_WRITE(log_fh,line,len_trim(line), &
                          MPI_CHARACTER,status,ierrmpi)

   write(line,'(e13.5)')wminte
   call MPI_FILE_WRITE(log_fh,line,len_trim(line), &
                          MPI_CHARACTER,status,ierrmpi)
   write(line,'(e13.5)')wmaxvel
   call MPI_FILE_WRITE(log_fh,line,len_trim(line), &
                          MPI_CHARACTER,status,ierrmpi)
   do level=1,mxnest
      write(line,'(e13.5)')volumeflat(level)
      call MPI_FILE_WRITE(log_fh,line,len_trim(line), &
                          MPI_CHARACTER,status,ierrmpi)
   end do
   do level=1,mxnest
      write(line,'(i6)') nleafs_level(level)
      call MPI_FILE_WRITE(log_fh,line,len_trim(line), &
                          MPI_CHARACTER,status,ierrmpi)
   end do

end if

end subroutine printlog_special
!=============================================================================
subroutine userspecialconvert(qunitconvert)

use mod_global_parameters
integer, intent(in) :: qunitconvert
character(len=20):: userconvert_type

integer  :: iigrid,igrid
logical  :: patchw(ixG^T)
!-----------------------------------------------------------------------------

if(mype==0)then
   print *,'converting to primitives, no specialvarout for integrals'
endif

do iigrid=1,igridstail; igrid=igrids(iigrid)
  !!call primitive(ixG^LL,ixG^LL^LSUB1,pw(igrid)%w,px(igrid)%x)
  call primitive(ixG^LL,ixG^LL,pw(igrid)%w,px(igrid)%x)
end do

call spatial_integral_w


patchw(ixG^T)=.false.
do iigrid=1,igridstail; igrid=igrids(iigrid)
  !!call conserve(ixG^LL,ixG^LL^LSUB1,pw(igrid)%w,px(igrid)%x,patchw)
  call conserve(ixG^LL,ixG^LL,pw(igrid)%w,px(igrid)%x,patchw)
end do

end subroutine userspecialconvert
{#IFDEF TRACER
!=============================================================================
subroutine mask_gridtrpos(ixI^L,ixO^L,w,x,patchwi)

use mod_global_parameters

integer, intent(in)                :: ixI^L,ixO^L
double precision, intent(in)       :: x(ixI^S,1:ndim)
double precision, intent(in)       :: w(ixI^S,1:nw)
logical, intent(inout)             :: patchwi(ixG^T)

double precision :: trtreshold
!-----------------------------------------------------------------------------
! note we have primitives in w-array!
trtreshold=0.95d0
patchwi(ixO^S)=(w(ixO^S,tr1_)>trtreshold)

return
end subroutine mask_gridtrpos
!=============================================================================
subroutine mask_gridtrzer(ixI^L,ixO^L,w,x,patchwi)

use mod_global_parameters

integer, intent(in)                :: ixI^L,ixO^L
double precision, intent(in)       :: x(ixI^S,1:ndim)
double precision, intent(in)       :: w(ixI^S,1:nw)
logical, intent(inout)             :: patchwi(ixG^T)

double precision :: trtreshold,trtreshold2
!-----------------------------------------------------------------------------
! note we have primitives in w-array!
trtreshold=0.05d0
trtreshold2=-0.05d0
patchwi(ixO^S)=(w(ixO^S,tr1_)<trtreshold.and.w(ixO^S,tr1_)>trtreshold2)

return
end subroutine mask_gridtrzer
!=============================================================================
}
subroutine spatial_integral_w

use mod_global_parameters

double precision :: dvolume(ixG^T), timephy,xmom,ymom
double precision, allocatable :: integral_ipe(:), integral_w(:)
double precision, external :: integral_grid

integer           :: nregions,ireg
integer           :: iigrid,igrid,status(MPI_STATUS_SIZE),ni
character(len=100):: filename,region
logical           :: patchwi(ixG^T),alive
!-----------------------------------------------------------------------------
nregions=1

{#IFDEF TRACER
nregions=3
}
do ireg=1,nregions
 select case(ireg)
 case(1)
   region='fulldomain'
{#IFDEF TRACER
 case(2)
   region='trpos'
 case(3)
   region='trzer'
}
 end select
! number of integrals to perform
ni=12
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
  if (.not.slab) mygeo => pgeo(igrid)
  if (B0field) then
     myB0_cell => pB0_cell(igrid)
    {^D&myB0_face^D => pB0_face^D(igrid)\}
  end if
  typelimiter=typelimiter1(node(plevel_,igrid))
  typegradlimiter=typegradlimiter1(node(plevel_,igrid))
  patchwi(ixG^T)=.false.
  select case(region)
{#IFDEF TRACER
  case('trpos')
     call mask_gridtrpos(ixG^LL,ixM^LL,pw(igrid)%w,px(igrid)%x,patchwi)
  case('trzer')
     call mask_gridtrzer(ixG^LL,ixM^LL,pw(igrid)%w,px(igrid)%x,patchwi)
}
  case('fulldomain')
     patchwi(ixM^T)=.true.
  case default
     call mpistop("region not defined")
  end select
  integral_ipe(1)=integral_ipe(1)+ &
            integral_grid(ixG^LL,ixM^LL,pw(igrid)%w,px(igrid)%x,dvolume,1,patchwi)
  integral_ipe(2)=integral_ipe(2)+ &
            integral_grid(ixG^LL,ixM^LL,pw(igrid)%w,px(igrid)%x,dvolume,2,patchwi)
  integral_ipe(3)=integral_ipe(3)+ &
            integral_grid(ixG^LL,ixM^LL,pw(igrid)%w,px(igrid)%x,dvolume,3,patchwi)
  integral_ipe(4)=integral_ipe(4)+ &
            integral_grid(ixG^LL,ixM^LL,pw(igrid)%w,px(igrid)%x,dvolume,4,patchwi)
  integral_ipe(5)=integral_ipe(5)+ &
            integral_grid(ixG^LL,ixM^LL,pw(igrid)%w,px(igrid)%x,dvolume,5,patchwi)
  integral_ipe(6)=integral_ipe(6)+ &
            integral_grid(ixG^LL,ixM^LL,pw(igrid)%w,px(igrid)%x,dvolume,6,patchwi)
  integral_ipe(7)=integral_ipe(7)+ &
            integral_grid(ixG^LL,ixM^LL,pw(igrid)%w,px(igrid)%x,dvolume,7,patchwi)
  integral_ipe(8)=integral_ipe(8)+ &
            integral_grid(ixG^LL,ixM^LL,pw(igrid)%w,px(igrid)%x,dvolume,8,patchwi)
  integral_ipe(9)=integral_ipe(9)+ &
            integral_grid(ixG^LL,ixM^LL,pw(igrid)%w,px(igrid)%x,dvolume,9,patchwi)
  integral_ipe(10)=integral_ipe(10)+ &
            integral_grid(ixG^LL,ixM^LL,pw(igrid)%w,px(igrid)%x,dvolume,10,patchwi)
  integral_ipe(11)=integral_ipe(11)+ &
            integral_grid(ixG^LL,ixM^LL,pw(igrid)%w,px(igrid)%x,dvolume,11,patchwi)
  integral_ipe(12)=integral_ipe(12)+ &
            integral_grid(ixG^LL,ixM^LL,pw(igrid)%w,px(igrid)%x,dvolume,12,patchwi)
end do
call MPI_ALLREDUCE(integral_ipe,integral_w,ni,MPI_DOUBLE_PRECISION,&
                     MPI_SUM,icomm,ierrmpi)
!!timephy=t*tunit
timephy=t
if(mype==0) then
  !!print *,'t=',timephy,'integral_w=',integral_w
  write(filename,"(a,a,a)") TRIM(filenamelog),TRIM(region),".int"
  inquire(file=filename,exist=alive)
  if(alive) then
    open(unit=21,file=filename,form='formatted',status='old',access='append')
  else
    open(unit=21,file=filename,form='formatted',status='new')
  endif
  write(21,'(13(es12.4))') timephy,integral_w(1:12)
  close(21)
endif

deallocate(integral_ipe,integral_w)

enddo
end subroutine spatial_integral_w
!=============================================================================
function integral_grid(ixI^L,ixO^L,w,x,dvolume,intval,patchwi)

use mod_global_parameters

integer, intent(in)                :: ixI^L,ixO^L,intval
double precision, intent(in)       :: x(ixI^S,1:ndim),dvolume(ixG^T)
double precision, intent(in)       :: w(ixI^S,nw)
logical, intent(in) :: patchwi(ixG^T)

double precision, dimension(ixG^T,1:ndir) :: bvec,current
double precision :: integral_grid,gammm1
integer :: ix^D,idirmin,idir,jdir,kdir
!-----------------------------------------------------------------------------
gammm1=eqpar(gamma_)-one

integral_grid=0.d0
select case(intval)
 case(1)
  ! volume integration
  {do ix^DB=ixOmin^DB,ixOmax^DB\}
     if(patchwi(ix^D)) integral_grid=integral_grid+dvolume(ix^D)
  {end do\}
 case(2)
  ! mass integration
  {do ix^DB=ixOmin^DB,ixOmax^DB\}
     if(patchwi(ix^D)) integral_grid=integral_grid+w(ix^D,rho_)*dvolume(ix^D)
  {end do\}
 case(3)
  ! vertical velocities
  {do ix^DB=ixOmin^DB,ixOmax^DB\}
     if(patchwi(ix^D)) integral_grid=integral_grid+dabs(w(ix^D,v2_))*dvolume(ix^D)
  {end do\}
 case(4)
  ! horizontal velocities
  {do ix^DB=ixOmin^DB,ixOmax^DB\}
     if(patchwi(ix^D)) integral_grid=integral_grid+dsqrt(w(ix^D,v1_)**2+w(ix^D,v3_)**2)*dvolume(ix^D)
  {end do\}
 case(5)
  ! T integration
  {do ix^DB=ixOmin^DB,ixOmax^DB\}
     if(patchwi(ix^D)) integral_grid=integral_grid+(w(ix^D,p_)/w(ix^D,rho_))*dvolume(ix^D)
  {end do\}
 case(6)
  ! magnetic energy integration
  if(B0field) then
    ^C&bvec(ixI^S,^C)=w(ixI^S,b^C_)+myB0_cell%w(ixI^S,^C);
  else
    ^C&bvec(ixI^S,^C)=w(ixI^S,b^C_);
  endif
  {do ix^DB=ixOmin^DB,ixOmax^DB\}
     if(patchwi(ix^D)) integral_grid=integral_grid+ &
                half*(^C&bvec(ix^D,^C)**2+)*dvolume(ix^D)
  {end do\}
 case(7)
  ! kinetic energy integration
  {do ix^DB=ixOmin^DB,ixOmax^DB\}
     if(patchwi(ix^D)) integral_grid=integral_grid+ &
                half*w(ix^D,rho_)*(^C&w(ix^D,v0_+^C)**2+)*dvolume(ix^D)
  {end do\}
 case(8)
  ! plasma beta integration
  if(B0field) then
    ^C&bvec(ixI^S,^C)=w(ixI^S,b^C_)+myB0_cell%w(ixI^S,^C);
  else
    ^C&bvec(ixI^S,^C)=w(ixI^S,b^C_);
  endif
  {do ix^DB=ixOmin^DB,ixOmax^DB\}
     if(patchwi(ix^D)) integral_grid=integral_grid+ &
                (two*w(ix^D,p_)/(^C&bvec(ix^D,^C)**2+))*dvolume(ix^D)
  {end do\}
 case(9)
  ! J^2  integration (split case as well, current from B1 only))
  ^C&bvec(ixI^S,^C)=w(ixI^S,b^C_);
  call curlvector(bvec,ixI^L,ixO^L,current,idirmin,1,ndir)
  {do ix^DB=ixOmin^DB,ixOmax^DB\}
     if(patchwi(ix^D)) integral_grid=integral_grid+ &
                         (^C&current(ix^D,^C)**2+)*dvolume(ix^D)
  {end do\}
 case(10)
  ! internal energy integration
  {do ix^DB=ixOmin^DB,ixOmax^DB\}
     if(patchwi(ix^D)) integral_grid=integral_grid+(w(ix^D,p_)/gammm1)*dvolume(ix^D)
  {end do\}
 case(11)
  ! vertical flow with sign
  {do ix^DB=ixOmin^DB,ixOmax^DB\}
     if(patchwi(ix^D)) integral_grid=integral_grid+w(ix^D,v2_)*dvolume(ix^D)
  {end do\}
 case(12)
  ! vertical B component
  if(B0field) then
    ^C&bvec(ixI^S,^C)=w(ixI^S,b^C_)+myB0_cell%w(ixI^S,^C);
  else
    ^C&bvec(ixI^S,^C)=w(ixI^S,b^C_);
  endif
  {do ix^DB=ixOmin^DB,ixOmax^DB\}
     if(patchwi(ix^D)) integral_grid=integral_grid+dabs(bvec(ix^D,2))*dvolume(ix^D)
  {end do\}
 case default
     call mpistop("intval not defined")
end select

return
end function integral_grid
!=============================================================================
! amrvacusr.t.promRTideal
!=============================================================================
