subroutine set_coordinate_system(geom)
  use mod_global_parameters

  character(len=*), intent(in) :: geom

  select case (geom)
  case ("Cartesian","Cartesian_1D","Cartesian_2D","Cartesian_3D")
    ndir = ndim
    typeaxial='slab'
  case ("Cartesian_1.5D")
    if (ndim /= 1) call mpistop("Geometry Cartesian_1.5D but ndim /= 1")
    typeaxial='slab'
    ndir = 2
  case ("Cartesian_1.75D")
    if (ndim /= 1) call mpistop("Geometry Cartesian_1.75D but ndim /= 1")
    typeaxial='slab'
    ndir = 3
  case ("Cartesian_2.5D")
    if (ndim /= 2) call mpistop("Geometry Cartesian_2.5D but ndim /= 2")
    typeaxial='slab'
    ndir = 3
  case ("cylindrical","cylindrical_2D","cylindrical_3D")
    ndir = ndim
    r_   = 1
    z_   = 2
    if(ndir==3) phi_ = 3
    typeaxial='cylindrical'
  case ("cylindrical_2.5D")
    if (ndim /= 2) call mpistop("Geometry cylindrical_2.5D but ndim /= 2")
    ndir = 3
    r_   = 1
    z_   = 2
    phi_ = 3
    typeaxial='cylindrical'
  case ("polar","polar_2D","polar_3D")
    ndir = ndim
    r_   = 1
    phi_ = 2
    if(ndir==3) z_ = 3
    typeaxial='cylindrical'
  case ("polar_2.5D")
    if (ndim /= 2) call mpistop("Geometry polar_2.5D but ndim /= 2")
    ndir = 3
    r_   = 1
    phi_ = 2
    z_   = 3
    typeaxial='cylindrical'
  case ("spherical","spherical_2D","spherical_3D")
    ndir = ndim
    r_   = 1
    if(ndir==3) phi_ = 3
    z_   = -1
    typeaxial='spherical'
  case ("spherical_2.5D")
    if (ndim /= 2) &
         call mpistop("Geometry spherical_2.5D requires ndim == 2")
    ndir = 3
    r_   = 1
    phi_ = 3
    z_   = -1
    typeaxial='spherical'
  case default
    call mpistop("Unknown geometry specified")
  end select
end subroutine set_coordinate_system

subroutine set_pole

  use mod_global_parameters
  
  select case (typeaxial)
  case ("spherical") {^IFTHREED
    ! For spherical grid, check whether phi-direction is periodic
    if(periodB(ndim)) then
      if(phi_/=3) call mpistop("phi_ should be 3 in 3D spherical coord!")
      if(mod(ng3(1),2)/=0) &
        call mpistop("Number of meshes in phi-direction should be even!")
      if(abs(xprobmin2)<smalldouble) then
        if(mype==0) write(unitterm,*) &
          "Will apply pi-periodic conditions at northpole!"
        poleB(1,2)=.true.
      else
        if(mype==0) write(unitterm,*) "There is no northpole!"
      end if
      if(abs(xprobmax2-dpi)<smalldouble) then
        if(mype==0) write(unitterm,*) &
          "Will apply pi-periodic conditions at southpole!"
        poleB(2,2)=.true.
      else
        if(mype==0) write(unitterm,*) "There is no southpole!"
      end if
    end if}
  case ("cylindrical")
    {if (^D == phi_ .and. periodB(^D)) then
      if(mod(ng^D(1),2)/=0) then
        call mpistop("Number of meshes in phi-direction should be even!")
      end if

      if(abs(xprobmin1)<smalldouble) then
        if (mype==0) then
          write(unitterm,*) "Will apply pi-periodic conditions at r=0"
        end if
        poleB(1,1)=.true.
      else
        if (mype==0) then
          write(unitterm,*) "There is no cylindrical axis!"
        end if
      end if
    end if\}
  end select

end subroutine set_pole
!=============================================================================
subroutine putgridgeo(igrid)

use mod_global_parameters

integer, intent(in) :: igrid
!-----------------------------------------------------------------------------
deallocate(pw(igrid)%surfaceC,pw(igrid)%surface,&
     pw(igrid)%dvolume,pw(igrid)%dx,pw(igrid)%dxcoarse,pw(igrid)%dvolumecoarse)

end subroutine putgridgeo
!=============================================================================
subroutine fillgeo(igrid,ixG^L)

use mod_global_parameters

integer, intent(in) :: igrid, ixG^L

integer :: ix^L, ixC^L
double precision :: x(ixG^S,ndim), drs(ixG^S), dx2(ixG^S), dx3(ixG^S)
!-----------------------------------------------------------------------------
ix^L=ixG^L^LSUB1;

select case (typeaxial)
case ("slabstretch")
   drs(ixG^S)=pw(igrid)%dx(ixG^S,1)
   {^NOONED
   dx2(ixG^S)=pw(igrid)%dx(ixG^S,2)}
   {^IFTHREED
   dx3(ixG^S)=pw(igrid)%dx(ixG^S,3)}

   {^IFONED
   ixCmin^D=ixmin^D-kr(^D,1); ixCmax^D=ixmax^D;
   pw(igrid)%surfaceC(ixC^S,1)=1.d0
   pw(igrid)%surface(ixC^S,1) =1.d0
   }
   {^IFTWOD
   ixCmin^D=ixmin^D-kr(^D,1); ixCmax^D=ixmax^D;
   pw(igrid)%surfaceC(ixC^S,1)=dx2(ixC^S)
   pw(igrid)%surface(ixC^S,1) =dx2(ixC^S)
   ixCmin^D=ixmin^D-kr(^D,2); ixCmax^D=ixmax^D;
   pw(igrid)%surfaceC(ixC^S,2)=drs(ixC^S)
   pw(igrid)%surface(ixC^S,2)=drs(ixC^S)
   }
   {^IFTHREED
   ixCmin^D=ixmin^D-kr(^D,1); ixCmax^D=ixmax^D;
   pw(igrid)%surfaceC(ixC^S,1)= dx2(ixC^S)*dx3(ixC^S)
   pw(igrid)%surface(ixC^S,1)=pw(igrid)%surfaceC(ixC^S,1)
   ixCmin^D=ixmin^D-kr(^D,2); ixCmax^D=ixmax^D;
   pw(igrid)%surfaceC(ixC^S,2)= drs(ixC^S)*dx3(ixC^S)
   pw(igrid)%surface(ixC^S,2)=pw(igrid)%surfaceC(ixC^S,2)
   ixCmin^D=ixmin^D-kr(^D,3); ixCmax^D=ixmax^D;
   pw(igrid)%surfaceC(ixC^S,3)= drs(ixC^S)*dx2(ixC^S)
   pw(igrid)%surface(ixC^S,3)=pw(igrid)%surfaceC(ixC^S,3)
   }

case ("spherical")
   x(ixG^S,1)=pw(igrid)%x(ixG^S,1)
   {^NOONED
   x(ixG^S,2)=pw(igrid)%x(ixG^S,2)}

   drs(ixG^S)=pw(igrid)%dx(ixG^S,1)
   {^NOONED
   dx2(ixG^S)=pw(igrid)%dx(ixG^S,2)}
   {^IFTHREED
   dx3(ixG^S)=pw(igrid)%dx(ixG^S,3)}

   ixCmin^D=ixmin^D-kr(^D,1); ixCmax^D=ixmax^D;

   pw(igrid)%surfaceC(ixC^S,1)=(x(ixC^S,1)+half*drs(ixC^S))**2 {^NOONED &
                *two*dsin(x(ixC^S,2))*dsin(half*dx2(ixC^S))}{^IFTHREED*dx3(ixC^S)}

   {^NOONED
   ixCmin^D=ixmin^D-kr(^D,2); ixCmax^D=ixmax^D;
   pw(igrid)%surfaceC(ixC^S,2)=x(ixC^S,1)*drs(ixC^S)&
                *dsin(x(ixC^S,2)+half*dx2(ixC^S))}{^IFTHREED*dx3(ixC^S)}

   {^IFTHREED
   ixCmin^D=ixmin^D-kr(^D,3); ixCmax^D=ixmax^D;
   pw(igrid)%surfaceC(ixC^S,3)=x(ixC^S,1)*drs(ixC^S)*dx2(ixC^S)}

   ixCmin^D=ixmin^D-kr(^D,1); ixCmax^D=ixmax^D;
   pw(igrid)%surface(ixC^S,1)=x(ixC^S,1)**2 {^NOONED &
              *two*dsin(x(ixC^S,2))*dsin(half*dx2(ixC^S))}{^IFTHREED*dx3(ixC^S)}
   {^NOONED
   ixCmin^D=ixmin^D-kr(^D,2); ixCmax^D=ixmax^D;
   pw(igrid)%surface(ixC^S,2)=x(ixC^S,1)*drs(ixC^S)&
                *dsin(x(ixC^S,2))}{^IFTHREED*dx3(ixC^S)}

   {^IFTHREED
   ixCmin^D=ixmin^D-kr(^D,3); ixCmax^D=ixmax^D;
   pw(igrid)%surface(ixC^S,3)=x(ixC^S,1)*drs(ixC^S)*dx2(ixC^S)}

case ("cylindrical")
   x(ixG^S,1)=pw(igrid)%x(ixG^S,1)
   drs(ixG^S)=pw(igrid)%dx(ixG^S,1)
   {^NOONED
   dx2(ixG^S)=pw(igrid)%dx(ixG^S,2)}
   {^IFTHREED
   dx3(ixG^S)=pw(igrid)%dx(ixG^S,3)}

   ixCmin^D=ixmin^D-kr(^D,1); ixCmax^D=ixmax^D;
   pw(igrid)%surfaceC(ixC^S,1)=dabs(x(ixC^S,1)+half*drs(ixC^S)){^DE&*dx^DE(ixC^S) }
   {^NOONED
   ixCmin^D=ixmin^D-kr(^D,2); ixCmax^D=ixmax^D;
   if (z_==2) pw(igrid)%surfaceC(ixC^S,2)=x(ixC^S,1)*drs(ixC^S){^IFTHREED*dx3(ixC^S)}
   if (phi_==2) pw(igrid)%surfaceC(ixC^S,2)=drs(ixC^S){^IFTHREED*dx3(ixC^S)}}
   {^IFTHREED
   ixCmin^D=ixmin^D-kr(^D,3); ixCmax^D=ixmax^D;
   if (z_==3) pw(igrid)%surfaceC(ixC^S,3)=x(ixC^S,1)*drs(ixC^S)*dx2(ixC^S)
   if (phi_==3) pw(igrid)%surfaceC(ixC^S,3)=drs(ixC^S)*dx2(ixC^S)}

   ixCmin^D=ixmin^D-kr(^D,1); ixCmax^D=ixmax^D;
   pw(igrid)%surface(ixC^S,1)=dabs(x(ixC^S,1)){^DE&*dx^DE(ixC^S) }
   {^NOONED
   ixCmin^D=ixmin^D-kr(^D,2); ixCmax^D=ixmax^D;
   if (z_==2) pw(igrid)%surface(ixC^S,2)=x(ixC^S,1)*drs(ixC^S){^IFTHREED*dx3(ixC^S)}
   if (phi_==2) pw(igrid)%surface(ixC^S,2)=drs(ixC^S){^IFTHREED*dx3(ixC^S)}}
   {^IFTHREED
   ixCmin^D=ixmin^D-kr(^D,3); ixCmax^D=ixmax^D;
   if (z_==3) pw(igrid)%surface(ixC^S,3)=x(ixC^S,1)*drs(ixC^S)*dx2(ixC^S)
   if (phi_==3) pw(igrid)%surface(ixC^S,3)=drs(ixC^S)*dx2(ixC^S)}

case default
   call mpistop("Sorry, typeaxial unknown")
end select

end subroutine fillgeo
!=============================================================================
subroutine gradient(q,ixI^L,ix^L,idir,gradq)

! Calculate gradient of a scalar q within ixL in direction idir

use mod_global_parameters

integer :: ixI^L, ix^L, idir
double precision :: q(ixI^S), gradq(ixI^S)

double precision :: qC(ixI^S),invdx
integer :: jx^L, hx^L, ixC^L, jxC^L {#IFDEF FOURTHORDER , lx^L, kx^L}

!-----------------------------------------------------------------------------

invdx=1.d0/dxlevel(idir)
if (slab) then
{#IFNDEF FOURTHORDER
   jx^L=ix^L+kr(idir,^D);
   hx^L=ix^L-kr(idir,^D);
   gradq(ix^S) = half*(q(jx^S)-q(hx^S))*invdx
}{#IFDEF FOURTHORDER
   lx^L=ix^L+2*kr(idir,^D);
   jx^L=ix^L+kr(idir,^D);
   hx^L=ix^L-kr(idir,^D);
   kx^L=ix^L-2*kr(idir,^D);
   gradq(ix^S) = (-q(lx^S) + 8.0d0 * q(jx^S) - 8.0d0 * q(hx^S) + q(kx^S)) &
        /(12.0d0 * dxlevel(idir))
}
else
   hx^L=ix^L-kr(idir,^D);
   ixCmin^D=hxmin^D;ixCmax^D=ixmax^D;
   jxC^L=ixC^L+kr(idir,^D);
   qC(ixC^S)=block%surfaceC(ixC^S,idir)*half*(q(ixC^S)+q(jxC^S))
   gradq(ix^S)=(qC(ix^S)-qC(hx^S))/block%dvolume(ix^S)
end if

end subroutine gradient
!=============================================================================
subroutine gradient_s(q,ixI^L,ixO^L,idims,gradq)

! Calculate gradient of a scalar q within ixL in direction idir

use mod_global_parameters

integer :: ixI^L, ixO^L, idims
double precision :: q(ixI^S), gradq(ixI^S)

double precision :: qC(ixI^S),qpoint(ixI^S),qface(ixI^S),invdx
integer :: ixA^L, ixB^L, ixC^L, jxC^L, ix^D

!-----------------------------------------------------------------------------

invdx=1.d0/dxlevel(idims)
! ixC is cell-corner index
ixCmax^D=ixOmax^D; ixCmin^D=ixOmin^D-1;
if (slab) then
  ! b unit vector at cell corner
  qpoint=0.d0
  {do ix^DB=0,1\}
    ixAmin^D=ixCmin^D+ix^D;
    ixAmax^D=ixCmax^D+ix^D;
    qpoint(ixC^S)=qpoint(ixC^S)+q(ixA^S)
  {end do\}
  qpoint(ixC^S)=qpoint(ixC^S)*0.5d0**ndim
  ! values at cell face
  qface=0.d0
  ixB^L=ixO^L-kr(idims,^D);
  ixAmax^D=ixOmax^D; ixAmin^D=ixBmin^D;
  ixB^L=ixA^L;
  {do ix^DB=0,1 \}
     if({ ix^D==0 .and. ^D==idims | .or.}) then
       ixBmin^D=ixAmin^D-ix^D; 
       ixBmax^D=ixAmax^D-ix^D; 
       qface(ixA^S)=qface(ixA^S)+qpoint(ixB^S)
     end if
  {end do\}
  qface(ixA^S)=qface(ixA^S)*0.5d0**(ndim-1)
  ixB^L=ixO^L-kr(idims,^D);
  gradq(ixO^S)=invdx*(qface(ixO^S)-qface(ixB^S))
end if

end subroutine gradient_s
!=============================================================================
subroutine gradientS(q,ixI^L,ix^L,idir,gradq)

! Calculate gradient of a scalar q within ixL in direction idir
! first use limiter to go from cell center to edge

use mod_global_parameters
use mod_limiter
use mod_ppm

integer :: ixI^L, ix^L, idir
double precision :: q(ixI^S), gradq(ixI^S)

double precision,dimension(ixI^S):: qC,qL,qR,dqC,ldq,rdq
double precision :: invdx
integer :: hx^L,ixC^L,jxC^L,gxC^L,hxC^L

!-----------------------------------------------------------------------------

invdx=1.d0/dxlevel(idir)
hx^L=ix^L-kr(idir,^D);
ixCmin^D=hxmin^D;ixCmax^D=ixmax^D;
jxC^L=ixC^L+kr(idir,^D);
gxCmin^D=ixCmin^D-kr(idir,^D);gxCmax^D=jxCmax^D;
hxC^L=gxC^L+kr(idir,^D);

! set the gradient limiter here
qR(gxC^S) = q(hxC^S)
qL(gxC^S) = q(gxC^S)
if (typegradlimiter/=limiter_ppm) then
   dqC(gxC^S)= qR(gxC^S)-qL(gxC^S)
   call dwlimiter2(dqC,ixI^L,gxC^L,idir,typegradlimiter,ldw=ldq,rdw=rdq)
   qL(ixC^S) = qL(ixC^S) + half*ldq(ixC^S)
   qR(ixC^S) = qR(ixC^S) - half*rdq(jxC^S)
else
   call PPMlimitervar(ixG^LL,ixM^LL,idir,q,q,qL,qR)
endif

if (slab) then
   gradq(ix^S)=half*(qR(ix^S)-qL(hx^S))*invdx
else
   select case(idir)
   {case(^D)
    gradq(ix^S)=(qR(ix^S)-qL(hx^S))/block%dx(ix^S,idir) \}
   end select
end if

end subroutine gradientS
!=============================================================================
subroutine divvector(qvec,ixI^L,ixO^L,divq)

! Calculate divergence of a vector qvec within ixL

use mod_global_parameters

integer :: ixI^L,ixO^L
double precision :: qvec(ixI^S,1:ndir), divq(ixI^S)

double precision :: qC(ixI^S), invdx(1:ndim)
integer :: jxO^L, hxO^L, ixC^L, jxC^L, idims, ix^L {#IFDEF FOURTHORDER , gxO^L, kxO^L}
!-----------------------------------------------------------------------------
{#IFNDEF FOURTHORDER
ix^L=ixO^L^LADD1;
}{#IFDEF FOURTHORDER
ix^L=ixO^L^LADD2;
}
if (ixImin^D>ixmin^D.or.ixImax^D<ixmax^D|.or.) &
   call mpistop("Error in divvector: Non-conforming input limits")
invdx=1.d0/dxlevel
divq(ixO^S)=zero
if (slab) then
  do idims=1,ndim
{#IFNDEF FOURTHORDER
     jxO^L=ixO^L+kr(idims,^D);
     hxO^L=ixO^L-kr(idims,^D);
     divq(ixO^S)=divq(ixO^S)+half*(qvec(jxO^S,idims)-qvec(hxO^S,idims))*invdx(idims)
}{#IFDEF FOURTHORDER
     kxO^L=ixO^L+2*kr(idims,^D);
     jxO^L=ixO^L+kr(idims,^D);
     hxO^L=ixO^L-kr(idims,^D);
     gxO^L=ixO^L-2*kr(idims,^D);
     divq(ixO^S)=divq(ixO^S)+&
          (-qvec(kxO^S,idims) + 8.0d0 * qvec(jxO^S,idims) - 8.0d0 * qvec(hxO^S,idims) + qvec(gxO^S,idims))/(12.0d0 * dxlevel(idims))
}
  end do
else
  do idims=1,ndim
     hxO^L=ixO^L-kr(idims,^D);
     ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
     jxC^L=ixC^L+kr(idims,^D);
     if(stretched_dim(idims)) then
       ! linear interpolation at cell interface along stretched dimension 
       qC(ixC^S)=block%surfaceC(ixC^S,idims)*(qvec(ixC^S,idims)+0.5d0*block%dx(ixC^S,idims)*&
        (qvec(jxC^S,idims)-qvec(ixC^S,idims))/(block%x(jxC^S,idims)-block%x(ixC^S,idims)))
     else
       qC(ixC^S)=block%surfaceC(ixC^S,idims)*half*(qvec(ixC^S,idims)+qvec(jxC^S,idims))
     end if
     divq(ixO^S)=divq(ixO^S)+qC(ixO^S)-qC(hxO^S)
  end do
  divq(ixO^S)=divq(ixO^S)/block%dvolume(ixO^S)
end if


end subroutine divvector 
!=============================================================================
subroutine curlvector(qvec,ixI^L,ixO^L,curlvec,idirmin,idirmin0,ndir0)

! Calculate curl of a vector qvec within ixL
! Using Stokes' theorem for non-Cartesian grids

use mod_global_parameters

integer :: ixI^L,ixO^L,idirmin,ixC^L,jxC^L,idir,jdir,kdir,hxO^L,jxO^L,ndir0,idirmin0
double precision :: qvec(ixI^S,1:ndir0),curlvec(ixI^S,idirmin0:3), invdx(1:ndim)
double precision :: tmp(ixI^S),xC(ixI^S)
!-----------------------------------------------------------------------------

! Calculate curl within ixL: CurlV_i=eps_ijk*d_j V_k
! Curl can have components (idirmin0:3)
! Determine exact value of idirmin while doing the loop.

idirmin=4
curlvec(ixO^S,idirmin0:3)=zero

if(slab) then
  invdx=1.d0/dxlevel
  do idir=idirmin0,3; do jdir=1,ndim; do kdir=1,ndir0
     if(lvc(idir,jdir,kdir)/=0)then
       tmp(ixI^S)=qvec(ixI^S,kdir)
       hxO^L=ixO^L-kr(jdir,^D);
       jxO^L=ixO^L+kr(jdir,^D);
       tmp(ixO^S)=half*(tmp(jxO^S)-tmp(hxO^S))*invdx(jdir)
       if(lvc(idir,jdir,kdir)==1)then
         curlvec(ixO^S,idir)=curlvec(ixO^S,idir)+tmp(ixO^S)
       else
         curlvec(ixO^S,idir)=curlvec(ixO^S,idir)-tmp(ixO^S)
       endif
       if(idir<idirmin)idirmin=idir
     endif
  enddo; enddo; enddo;
else
  do idir=idirmin0,3; do jdir=1,ndim; do kdir=1,ndir0
    if(lvc(idir,jdir,kdir)/=0)then
      select case(typeaxial)
        case('slabstretch')
          select case(idir)
            case(1)
              if(jdir<kdir) then
                ! idir=1,jdir=2,kdir=3
                !! integral along 3rd dimension
                hxO^L=ixO^L-kr(jdir,^D);
                ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                jxC^L=ixC^L+kr(jdir,^D);
                ! qvec(3) at cell interface along 2nd dimension
                tmp(ixC^S)=qvec(ixC^S,kdir)+0.5d0*block%dx(ixC^S,jdir)*&
                 (qvec(jxC^S,kdir)-qvec(ixC^S,kdir))/(block%x(jxC^S,jdir)-block%x(ixC^S,jdir))
                curlvec(ixO^S,idir)=(tmp(ixO^S)-tmp(hxO^S))*block%dx(ixO^S,kdir)
                !! integral along 2nd dimension
                hxO^L=ixO^L-kr(kdir,^D);
                ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                jxC^L=ixC^L+kr(kdir,^D);                
                ! qvec(2) at cell interface along 3rd dimension
                tmp(ixC^S)=qvec(ixC^S,jdir)+0.5d0*block%dx(ixC^S,kdir)*&
                 (qvec(jxC^S,jdir)-qvec(ixC^S,jdir))/(block%x(jxC^S,kdir)-block%x(ixC^S,kdir))
                curlvec(ixO^S,idir)=(curlvec(ixO^S,idir)+(tmp(hxO^S)-tmp(ixO^S))*block%dx(ixO^S,jdir))&
                  /block%surface(ixO^S,idir)
              end if
            case(2)
              if(jdir<kdir) then
                ! idir=2,jdir=1,kdir=3
                !! integral along 1st dimension
                hxO^L=ixO^L-kr(kdir,^D);
                ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                jxC^L=ixC^L+kr(kdir,^D);                
                ! qvec(1) at cell interface along 3rd dimension
                tmp(ixC^S)=qvec(ixC^S,jdir)+0.5d0*block%dx(ixC^S,kdir)*&
                 (qvec(jxC^S,jdir)-qvec(ixC^S,jdir))/(block%x(jxC^S,kdir)-block%x(ixC^S,kdir))
                curlvec(ixO^S,idir)=(tmp(ixO^S)-tmp(hxO^S))*block%dx(ixO^S,jdir)
                !! integral along 3rd dimension
                hxO^L=ixO^L-kr(jdir,^D);
                ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                jxC^L=ixC^L+kr(jdir,^D);
                ! qvec(3) at cell interface along 1st dimension
                tmp(ixC^S)=qvec(ixC^S,kdir)+0.5d0*block%dx(ixC^S,jdir)*&
                 (qvec(jxC^S,kdir)-qvec(ixC^S,kdir))/(block%x(jxC^S,jdir)-block%x(ixC^S,jdir))
                curlvec(ixO^S,idir)=(curlvec(ixO^S,idir)+(tmp(hxO^S)-tmp(ixO^S))*block%dx(ixO^S,kdir))&
                  /block%surface(ixO^S,idir)
              end if
            case(3)
              if(jdir<kdir) then
                ! idir=3,jdir=1,kdir=2
                !! integral along 1st dimension
                hxO^L=ixO^L-kr(kdir,^D);
                ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                jxC^L=ixC^L+kr(kdir,^D);                
                ! qvec(1) at cell interface along 2nd dimension
                tmp(ixC^S)=qvec(ixC^S,jdir)+0.5d0*block%dx(ixC^S,kdir)*&
                 (qvec(jxC^S,jdir)-qvec(ixC^S,jdir))/(block%x(jxC^S,kdir)-block%x(ixC^S,kdir))
                curlvec(ixO^S,idir)=(tmp(hxO^S)-tmp(ixO^S))*block%dx(ixO^S,jdir)
                !! integral along 2nd dimension
                hxO^L=ixO^L-kr(jdir,^D);
                ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                jxC^L=ixC^L+kr(jdir,^D);
                ! qvec(2) at cell interface along 1st dimension
                tmp(ixC^S)=qvec(ixC^S,kdir)+0.5d0*block%dx(ixC^S,jdir)*&
                 (qvec(jxC^S,kdir)-qvec(ixC^S,kdir))/(block%x(jxC^S,jdir)-block%x(ixC^S,jdir))
                curlvec(ixO^S,idir)=(curlvec(ixO^S,idir)+(tmp(ixO^S)-tmp(hxO^S))*block%dx(ixO^S,kdir))&
                   /block%surface(ixO^S,idir)
              end if
         end select
        case('spherical')
          select case(idir)
            case(1)
              if(jdir<kdir) then
                ! idir=1,jdir=2,kdir=3
                !! integral along 3rd dimension
                hxO^L=ixO^L-kr(jdir,^D);
                ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                jxC^L=ixC^L+kr(jdir,^D);
                ! qvec(3) at cell interface along 2nd dimension
                if(stretched_dim(jdir)) then
                  tmp(ixC^S)=qvec(ixC^S,kdir)+0.5d0*block%dx(ixC^S,jdir)*&
                   (qvec(jxC^S,kdir)-qvec(ixC^S,kdir))/(block%x(jxC^S,jdir)-block%x(ixC^S,jdir))
                else
                  tmp(ixC^S)=0.5d0*(qvec(ixC^S,kdir)+qvec(jxC^S,kdir))
                end if
                ! 2nd coordinate at cell interface along 2nd dimension
                xC(ixC^S)=block%x(ixC^S,jdir)+0.5d0*block%dx(ixC^S,jdir)
                curlvec(ixO^S,idir)=(dsin(xC(ixO^S))*tmp(ixO^S)-&
                  dsin(xC(hxO^S))*tmp(hxO^S))*block%dx(ixO^S,kdir)
                !! integral along 2nd dimension
                hxO^L=ixO^L-kr(kdir,^D);
                ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                jxC^L=ixC^L+kr(kdir,^D);                
                ! qvec(2) at cell interface along 3rd dimension
                if(stretched_dim(kdir)) then
                  tmp(ixC^S)=qvec(ixC^S,jdir)+0.5d0*block%dx(ixC^S,kdir)*&
                   (qvec(jxC^S,jdir)-qvec(ixC^S,jdir))/(block%x(jxC^S,kdir)-block%x(ixC^S,kdir))
                else
                  tmp(ixC^S)=0.5d0*(qvec(ixC^S,jdir)+qvec(jxC^S,jdir))
                end if
                curlvec(ixO^S,idir)=(curlvec(ixO^S,idir)+(tmp(hxO^S)-tmp(ixO^S))*block%dx(ixO^S,jdir))&
                  /block%surface(ixO^S,idir)*block%x(ixO^S,idir)
              end if
            case(2)
              if(jdir<kdir) then
                ! idir=2,jdir=1,kdir=3
                !! integral along 1st dimension
                hxO^L=ixO^L-kr(kdir,^D);
                ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                jxC^L=ixC^L+kr(kdir,^D);                
                ! qvec(1) at cell interface along 3rd dimension
                if(stretched_dim(kdir)) then
                  tmp(ixC^S)=qvec(ixC^S,jdir)+0.5d0*block%dx(ixC^S,kdir)*&
                   (qvec(jxC^S,jdir)-qvec(ixC^S,jdir))/(block%x(jxC^S,kdir)-block%x(ixC^S,kdir))
                else
                  tmp(ixC^S)=0.5d0*(qvec(ixC^S,jdir)+qvec(jxC^S,jdir))
                end if
                curlvec(ixO^S,idir)=(tmp(ixO^S)-tmp(hxO^S))*block%dx(ixO^S,1)/block%surface(ixO^S,idir)
                !! integral along 3rd dimension
                hxO^L=ixO^L-kr(jdir,^D);
                ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                jxC^L=ixC^L+kr(jdir,^D);
                ! qvec(3) at cell interface along 1st dimension
                if(stretched_dim(jdir)) then
                  tmp(ixC^S)=qvec(ixC^S,kdir)+0.5d0*block%dx(ixC^S,jdir)*&
                   (qvec(jxC^S,kdir)-qvec(ixC^S,kdir))/(block%x(jxC^S,jdir)-block%x(ixC^S,jdir))
                else
                  tmp(ixC^S)=0.5d0*(qvec(ixC^S,kdir)+qvec(jxC^S,kdir))
                end if
                ! 1st coordinate at cell interface along 1st dimension
                xC(ixC^S)=block%x(ixC^S,jdir)+0.5d0*block%dx(ixC^S,jdir)
                curlvec(ixO^S,idir)=curlvec(ixO^S,idir)+(xC(hxO^S)*tmp(hxO^S)-xC(ixO^S)*tmp(ixO^S))/(block%x(ixO^S,jdir)*block%dx(ixO^S,jdir))
              end if
            case(3)
              if(jdir<kdir) then
                ! idir=3,jdir=1,kdir=2
                !! integral along 1st dimension
                hxO^L=ixO^L-kr(kdir,^D);
                ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                jxC^L=ixC^L+kr(kdir,^D);                
                ! qvec(1) at cell interface along 2nd dimension
                if(stretched_dim(kdir)) then
                  tmp(ixC^S)=qvec(ixC^S,jdir)+0.5d0*block%dx(ixC^S,kdir)*&
                   (qvec(jxC^S,jdir)-qvec(ixC^S,jdir))/(block%x(jxC^S,kdir)-block%x(ixC^S,kdir))
                else
                  tmp(ixC^S)=0.5d0*(qvec(ixC^S,jdir)+qvec(jxC^S,jdir))
                end if
                curlvec(ixO^S,idir)=(tmp(hxO^S)-tmp(ixO^S))*block%dx(ixO^S,jdir)
                !! integral along 2nd dimension
                hxO^L=ixO^L-kr(jdir,^D);
                ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                jxC^L=ixC^L+kr(jdir,^D);
                ! qvec(2) at cell interface along 1st dimension
                if(stretched_dim(jdir)) then
                  tmp(ixC^S)=qvec(ixC^S,kdir)+0.5d0*block%dx(ixC^S,jdir)*&
                   (qvec(jxC^S,kdir)-qvec(ixC^S,kdir))/(block%x(jxC^S,jdir)-block%x(ixC^S,jdir))
                else
                  tmp(ixC^S)=0.5d0*(qvec(ixC^S,kdir)+qvec(jxC^S,kdir))
                end if
                ! 1st coordinate at cell interface along 1st dimension
                xC(ixC^S)=block%x(ixC^S,jdir)+0.5d0*block%dx(ixC^S,jdir)
                curlvec(ixO^S,idir)=(curlvec(ixO^S,idir)+(xC(ixO^S)*tmp(ixO^S)-xC(hxO^S)*tmp(hxO^S))*block%dx(ixO^S,kdir))&
                   /block%surface(ixO^S,idir)
              end if
         end select
        case('cylindrical')
          if(idir==r_) then
            if(jdir==phi_) then
              ! idir=r,jdir=phi,kdir=z
              !! integral along z dimension
              hxO^L=ixO^L-kr(jdir,^D);
              ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
              jxC^L=ixC^L+kr(jdir,^D);
              ! qvec(z) at cell interface along phi dimension
              if(stretched_dim(jdir)) then
                tmp(ixC^S)=qvec(ixC^S,kdir)+0.5d0*block%dx(ixC^S,jdir)*&
                 (qvec(jxC^S,kdir)-qvec(ixC^S,kdir))/(block%x(jxC^S,jdir)-block%x(ixC^S,jdir))
              else
                tmp(ixC^S)=0.5d0*(qvec(ixC^S,kdir)+qvec(jxC^S,kdir))
              end if
              curlvec(ixO^S,idir)=(tmp(ixO^S)-tmp(hxO^S))*block%dx(ixO^S,kdir)
              !! integral along phi dimension
              hxO^L=ixO^L-kr(kdir,^D);
              ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
              jxC^L=ixC^L+kr(kdir,^D);
              ! qvec(phi) at cell interface along z dimension
              if(stretched_dim(kdir)) then
                tmp(ixC^S)=qvec(ixC^S,jdir)+0.5d0*block%dx(ixC^S,kdir)*&
                 (qvec(jxC^S,jdir)-qvec(ixC^S,jdir))/(block%x(jxC^S,kdir)-block%x(ixC^S,kdir))
              else
                tmp(ixC^S)=0.5d0*(qvec(ixC^S,jdir)+qvec(jxC^S,jdir))
              end if
              curlvec(ixO^S,idir)=(curlvec(ixO^S,idir)+(tmp(hxO^S)-tmp(ixO^S))*block%x(ixO^S,idir)*block%dx(ixO^S,jdir))&
                /block%surface(ixO^S,idir)
            end if
          else if(idir==phi_) then
            if(jdir<kdir) then
              ! idir=phi,jdir=r,kdir=z
              !! integral along r dimension
              hxO^L=ixO^L-kr(kdir,^D);
              ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
              jxC^L=ixC^L+kr(kdir,^D);
              ! qvec(r) at cell interface along z dimension
              if(stretched_dim(kdir)) then
                tmp(ixC^S)=qvec(ixC^S,jdir)+0.5d0*block%dx(ixC^S,kdir)*&
                 (qvec(jxC^S,jdir)-qvec(ixC^S,jdir))/(block%x(jxC^S,kdir)-block%x(ixC^S,kdir))
              else
                tmp(ixC^S)=0.5d0*(qvec(ixC^S,jdir)+qvec(jxC^S,jdir))
              end if
              curlvec(ixO^S,idir)=(tmp(ixO^S)-tmp(hxO^S))*block%dx(ixO^S,jdir)
              !! integral along z dimension
              hxO^L=ixO^L-kr(jdir,^D);
              ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
              jxC^L=ixC^L+kr(jdir,^D);
              ! qvec(z) at cell interface along r dimension
              if(stretched_dim(jdir)) then
                tmp(ixC^S)=qvec(ixC^S,kdir)+0.5d0*block%dx(ixC^S,jdir)*&
                 (qvec(jxC^S,kdir)-qvec(ixC^S,kdir))/(block%x(jxC^S,jdir)-block%x(ixC^S,jdir))
              else
                tmp(ixC^S)=0.5d0*(qvec(ixC^S,kdir)+qvec(jxC^S,kdir))
              end if
              curlvec(ixO^S,idir)=(curlvec(ixO^S,idir)+(tmp(hxO^S)-tmp(ixO^S))*block%dx(ixO^S,kdir))&
                /block%surface(ixO^S,idir)
            end if
          else ! idir==z_
            if(jdir<kdir) then
              ! idir=z,jdir=r,kdir=phi
              !! integral along r dimension
              hxO^L=ixO^L-kr(kdir,^D);
              ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
              jxC^L=ixC^L+kr(kdir,^D);
              ! qvec(r) at cell interface along phi dimension
              if(stretched_dim(kdir)) then
                tmp(ixC^S)=qvec(ixC^S,jdir)+0.5d0*block%dx(ixC^S,kdir)*&
                 (qvec(jxC^S,jdir)-qvec(ixC^S,jdir))/(block%x(jxC^S,kdir)-block%x(ixC^S,kdir))
              else
                tmp(ixC^S)=0.5d0*(qvec(ixC^S,jdir)+qvec(jxC^S,jdir))
              end if
              curlvec(ixO^S,idir)=(tmp(hxO^S)-tmp(ixO^S))*block%dx(ixO^S,jdir)
              !! integral along phi dimension
              hxO^L=ixO^L-kr(jdir,^D);
              ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
              jxC^L=ixC^L+kr(jdir,^D);
              ! qvec(phi) at cell interface along r dimension
              if(stretched_dim(jdir)) then
                tmp(ixC^S)=qvec(ixC^S,kdir)+0.5d0*block%dx(ixC^S,jdir)*&
                 (qvec(jxC^S,kdir)-qvec(ixC^S,kdir))/(block%x(jxC^S,jdir)-block%x(ixC^S,jdir))
              else
                tmp(ixC^S)=0.5d0*(qvec(ixC^S,kdir)+qvec(jxC^S,kdir))
              end if
              ! r coordinate at cell interface along r dimension
              xC(ixC^S)=block%x(ixC^S,jdir)+0.5d0*block%dx(ixC^S,jdir)
              curlvec(ixO^S,idir)=(curlvec(ixO^S,idir)+(xC(ixO^S)*tmp(ixO^S)-xC(hxO^S)*tmp(hxO^S))*block%dx(ixO^S,kdir))&
                 /block%surface(ixO^S,idir)
            end if
          end if
        case default
          call mpistop("Sorry, typeaxial unknown in curlvector")
      end select
    endif
  enddo; enddo; enddo;
end if

end subroutine curlvector 
!=============================================================================
subroutine divvectorS(qvec,ixI^L,ixO^L,divq)

! Calculate divergence of a vector qvec within ixL
! using limited extrapolation to cell edges

use mod_global_parameters
use mod_limiter
use mod_ppm

integer :: ixI^L,ixO^L
double precision :: qvec(ixG^T,1:ndir), divq(ixG^T)

double precision,dimension(ixG^T):: qL,qR,dqC,ldq,rdq
double precision :: invdx(1:ndim)

integer :: hxO^L,ixC^L,jxC^L,idims,ix^L,gxC^L,hxC^L
!-----------------------------------------------------------------------------
ix^L=ixO^L^LADD2;

if (ixImin^D>ixmin^D.or.ixImax^D<ixmax^D|.or.) &
   call mpistop("Error in divvectorS: Non-conforming input limits")

invdx=1.d0/dxlevel
divq(ixO^S)=zero
do idims=1,ndim
   hxO^L=ixO^L-kr(idims,^D);
   ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
   jxC^L=ixC^L+kr(idims,^D);
   gxCmin^D=ixCmin^D-kr(idims,^D);gxCmax^D=jxCmax^D;
   hxC^L=gxC^L+kr(idims,^D);
   qR(gxC^S) = qvec(hxC^S,idims)
   qL(gxC^S) = qvec(gxC^S,idims)
   if(typegradlimiter/=limiter_ppm) then
      dqC(gxC^S)= qR(gxC^S)-qL(gxC^S)
      call dwlimiter2(dqC,ixI^L,gxC^L,idims,typegradlimiter,ldw=ldq,rdw=rdq)
      qL(ixC^S) = qL(ixC^S) + half*ldq(ixC^S)
      qR(ixC^S) = qR(ixC^S) - half*rdq(jxC^S)
   else
      dqC(ixI^S)=qvec(ixI^S,idims)
      call PPMlimitervar(ixG^LL,ixM^LL,idims,dqC,dqC,qL,qR)
   endif

   if (slab) then
     divq(ixO^S)=divq(ixO^S)+half*(qR(ixO^S)-qL(hxO^S))*invdx(idims)
   else
     qR(ixC^S)=block%surfaceC(ixC^S,idims)*qR(ixC^S)
     qL(ixC^S)=block%surfaceC(ixC^S,idims)*qL(ixC^S)
     divq(ixO^S)=divq(ixO^S)+qR(ixO^S)-qL(hxO^S)
   end if
end do
if(.not.slab) divq(ixO^S)=divq(ixO^S)/block%dvolume(ixO^S)

end subroutine divvectorS
!=============================================================================
!> cross product of two vectors
subroutine cross_product(ixI^L,ixO^L,a,b,axb)
  use mod_global_parameters

  integer, intent(in) :: ixI^L, ixO^L
  double precision, intent(in) :: a(ixI^S,3), b(ixI^S,3)
  double precision, intent(out) :: axb(ixI^S,3)
!-------------------------------------------------------------------------

  axb(ixO^S,1)=a(ixO^S,2)*b(ixO^S,3)-a(ixO^S,3)*b(ixO^S,2)
  axb(ixO^S,2)=a(ixO^S,3)*b(ixO^S,1)-a(ixO^S,1)*b(ixO^S,3)
  axb(ixO^S,3)=a(ixO^S,1)*b(ixO^S,2)-a(ixO^S,2)*b(ixO^S,1)

end subroutine cross_product
!=============================================================================
subroutine extremaq(ixI^L,ixO^L,q,nshift,qMax,qMin)

use mod_global_parameters

integer,intent(in)           :: ixI^L,ixO^L
double precision, intent(in) :: q(ixI^S)
integer,intent(in)           :: nshift

double precision, intent(out) :: qMax(ixI^S),qMin(ixI^S)

integer           :: ixs^L,ixsR^L,ixsL^L,idims,jdims,kdims,ishift,i,j 
!-------------------------------------------------------------------------
do ishift=1,nshift
 idims=1
 ixsR^L=ixO^L+ishift*kr(idims,^D);
 ixsL^L=ixO^L-ishift*kr(idims,^D);
 if (ishift==1) then
   qMax(ixO^S)=max(q(ixO^S),q(ixsR^S),q(ixsL^S))
   qMin(ixO^S)=min(q(ixO^S),q(ixsR^S),q(ixsL^S))
 else
   qMax(ixO^S)=max(qMax(ixO^S),q(ixsR^S),q(ixsL^S))
   qMin(ixO^S)=min(qMin(ixO^S),q(ixsR^S),q(ixsL^S))
 end if
 {^NOONED
 idims=1
 jdims=idims+1
 do i=-1,1
   ixs^L=ixO^L+i*ishift*kr(idims,^D);
   ixsR^L=ixs^L+ishift*kr(jdims,^D);
   ixsL^L=ixs^L-ishift*kr(jdims,^D);
   qMax(ixO^S)=max(qMax(ixO^S),q(ixsR^S),q(ixsL^S))
   qMin(ixO^S)=min(qMin(ixO^S),q(ixsR^S),q(ixsL^S))
 end do
 }
 {^IFTHREED
 idims=1
 jdims=idims+1
 kdims=jdims+1
 do i=-1,1
   ixs^L=ixO^L+i*ishift*kr(idims,^D);
   do j=-1,1
      ixs^L=ixO^L+j*ishift*kr(jdims,^D);
      ixsR^L=ixs^L+ishift*kr(kdims,^D);
      ixsL^L=ixs^L-ishift*kr(kdims,^D);
      qMax(ixO^S)=max(qMax(ixO^S),q(ixsR^S),q(ixsL^S))
      qMin(ixO^S)=min(qMin(ixO^S),q(ixsR^S),q(ixsL^S))
   end do
 end do
 }
enddo

end subroutine  extremaq
!=============================================================================
subroutine extremaw(ixI^L,ixO^L,w,nshift,wMax,wMin)

use mod_global_parameters

integer,intent(in)            :: ixI^L,ixO^L
double precision, intent(in)  :: w(ixI^S,1:nw)
integer,intent(in)            :: nshift

double precision, intent(out) :: wMax(ixI^S,1:nwflux),wMin(ixI^S,1:nwflux)

integer          :: ixs^L,ixsR^L,ixsL^L,idims,jdims,kdims,ishift,i,j
!-------------------------------------------------------------------------
do ishift=1,nshift
 idims=1
 ixsR^L=ixO^L+ishift*kr(idims,^D);
 ixsL^L=ixO^L-ishift*kr(idims,^D);
 if (ishift==1) then
    wMax(ixO^S,1:nwflux)= &
     max(w(ixO^S,1:nwflux),w(ixsR^S,1:nwflux),w(ixsL^S,1:nwflux))
    wMin(ixO^S,1:nwflux)= &
     min(w(ixO^S,1:nwflux),w(ixsR^S,1:nwflux),w(ixsL^S,1:nwflux))
 else
    wMax(ixO^S,1:nwflux)= &
     max(wMax(ixO^S,1:nwflux),w(ixsR^S,1:nwflux),w(ixsL^S,1:nwflux))
    wMin(ixO^S,1:nwflux)= &
     min(wMin(ixO^S,1:nwflux),w(ixsR^S,1:nwflux),w(ixsL^S,1:nwflux))
 end if
 {^NOONED
 idims=1
 jdims=idims+1
 do i=-1,1
   ixs^L=ixO^L+i*ishift*kr(idims,^D);
   ixsR^L=ixs^L+ishift*kr(jdims,^D);
   ixsL^L=ixs^L-ishift*kr(jdims,^D);
   wMax(ixO^S,1:nwflux)= &
     max(wMax(ixO^S,1:nwflux),w(ixsR^S,1:nwflux),w(ixsL^S,1:nwflux))
   wMin(ixO^S,1:nwflux)= &
     min(wMin(ixO^S,1:nwflux),w(ixsR^S,1:nwflux),w(ixsL^S,1:nwflux))
 end do
 }
 {^IFTHREED
 idims=1
 jdims=idims+1
 kdims=jdims+1
 do i=-1,1
   ixs^L=ixO^L+i*ishift*kr(idims,^D);
   do j=-1,1
      ixs^L=ixO^L+j*ishift*kr(jdims,^D);
      ixsR^L=ixs^L+ishift*kr(kdims,^D);
      ixsL^L=ixs^L-ishift*kr(kdims,^D);
      wMax(ixO^S,1:nwflux)= &
       max(wMax(ixO^S,1:nwflux),w(ixsR^S,1:nwflux),w(ixsL^S,1:nwflux))
      wMin(ixO^S,1:nwflux)= &
       min(wMin(ixO^S,1:nwflux),w(ixsR^S,1:nwflux),w(ixsL^S,1:nwflux))
   end do
 end do
 }
enddo

end subroutine  extremaw
!=============================================================================
subroutine extremaa(ixI^L,ixO^L,a,nshift,aMin)

use mod_global_parameters

integer,intent(in)           :: ixI^L,ixO^L
double precision, intent(in) :: a(ixI^S)
integer,intent(in)           :: nshift

double precision, intent(out) :: aMin(ixI^S)

integer          :: ixs^L,ixsR^L,ixsL^L,idims,jdims,kdims,ishift,i,j
!-------------------------------------------------------------------------
do ishift=1,nshift
  idims=1
  ixsR^L=ixO^L+ishift*kr(idims,^D);
  ixsL^L=ixO^L-ishift*kr(idims,^D);
  aMin(ixO^S)=min(a(ixsR^S),a(ixO^S),a(ixsL^S))
  {^NOONED
  idims=1
  jdims=idims+1
  do i=-1,1
    ixs^L=ixO^L+i*ishift*kr(idims,^D);
    ixsR^L=ixs^L+ishift*kr(jdims,^D);
    ixsL^L=ixs^L-ishift*kr(jdims,^D);
    aMin(ixO^S)=min(aMin(ixO^S),a(ixsR^S),a(ixsL^S))
  end do
  }
  {^IFTHREED
  idims=1
  jdims=idims+1
  kdims=jdims+1
  do i=-1,1
   ixs^L=ixO^L+i*ishift*kr(idims,^D);
   do j=-1,1
      ixs^L=ixO^L+j*ishift*kr(jdims,^D);
      ixsR^L=ixs^L+ishift*kr(kdims,^D);
      ixsL^L=ixs^L-ishift*kr(kdims,^D);
      aMin(ixO^S)=min(aMin(ixO^S),a(ixsR^S),a(ixsL^S))
   end do
  end do
  }
end do

end subroutine extremaa
!=============================================================================
subroutine locate_in_table(xpoint,table,nxtp,ixtp,res)
!  Fast search table to find xpoint's location in the table
!  record the distance between xpoint and the its cloest element table(ixtp)
! INPUT:
!  xpoint : the point you want to locate in the table
!  table : 1D table you want to find location in
!  nxtp : number of elements in the table
! OUTPUT:
!  ixtp : closest element's index in table
!  res : offset(distance) from the closest element in table
use mod_global_parameters

integer, intent(in) :: nxtp
double precision,intent(in)   :: xpoint,table(nxtp)
double precision, intent(out) :: res
integer, intent(out) :: ixtp

integer :: jl,jc,jh
!-----------------------------------------------------------------------------
if(xpoint < table(1)) then
  ixtp=1
  res=(xpoint-table(1))/(table(2)-table(1))
else if (xpoint > table(nxtp)) then
  ixtp=nxtp
  res=(xpoint-table(nxtp))/(table(nxtp)-table(nxtp-1))
else
  jl=0
  jh=nxtp+1
  do
    if (jh-jl <= 1) exit
    jc=(jh+jl)/2
    if (xpoint >= table(jc)) then
        jl=jc
    else
        jh=jc
    end if
  end do
  res=(xpoint-table(jh-1))/(table(jh)-table(jh-1))
  if(res<=0.5d0) then
    ixtp=jh-1
  else
    ixtp=jh
    res=res-1.d0
  endif
end if
end subroutine locate_in_table
!=============================================================================
subroutine divvector_s(qvec,ixI^L,ixO^L,divq)

! Calculate divergence of a vector qvec within ixL

use mod_global_parameters

integer :: ixI^L,ixO^L
double precision :: qvec(ixI^S,1:ndir), divq(ixI^S)

double precision :: invdx(1:ndim)
double precision :: qpoint(ixI^S,1:ndim),qface(ixI^S,1:ndim)
integer :: ixA^L, ixB^L, ixC^L, idims, ix^L,ix^D

ix^L=ixO^L^LADD1;
if (ixImin^D>ixmin^D.or.ixImax^D<ixmax^D|.or.) &
   call mpistop("Error in divvector: Non-conforming input limits")
invdx=1.d0/dxlevel
divq(ixO^S)=zero
! ixC is cell-corner index
ixCmax^D=ixOmax^D; ixCmin^D=ixOmin^D-1;
if (slab) then
  ! vector at cell corner
  qpoint=0.d0
  {do ix^DB=0,1\}
    ixAmin^D=ixCmin^D+ix^D;
    ixAmax^D=ixCmax^D+ix^D;
    qpoint(ixC^S,1:ndim)=qpoint(ixC^S,1:ndim)+qvec(ixA^S,1:ndim)
  {end do\}
  qpoint(ixC^S,1:ndim)=qpoint(ixC^S,1:ndim)*0.5d0**ndim
  ! values at cell face
  qface=0.d0
  do idims=1,ndim
    ixB^L=ixO^L-kr(idims,^D);
    ixAmax^D=ixOmax^D; ixAmin^D=ixBmin^D;
    ixB^L=ixA^L;
    {do ix^DB=0,1 \}
       if({ ix^D==0 .and. ^D==idims | .or.}) then
         ixBmin^D=ixAmin^D-ix^D; 
         ixBmax^D=ixAmax^D-ix^D; 
         qface(ixA^S,idims)=qface(ixA^S,idims)+qpoint(ixB^S,idims)
       end if
    {end do\}
    qface(ixA^S,idims)=qface(ixA^S,idims)*0.5d0**(ndim-1)
  end do
  divq=0.d0
  do idims=1,ndim
    ixB^L=ixO^L-kr(idims,^D);
    divq(ixO^S)=divq(ixO^S)+invdx(idims)*(qface(ixO^S,idims)-qface(ixB^S,idims))
  end do
else
  divq=0.d0
  do idims=1,ndim
    ixB^L=ixO^L-kr(idims,^D);
    ixCmin^D=ixBmin^D;ixCmax^D=ixOmax^D;
    ixA^L=ixC^L+kr(idims,^D);
    qface(ixC^S,idims)=block%surfaceC(ixC^S,idims)*half*(qvec(ixC^S,idims)+qvec(ixA^S,idims))
    divq(ixO^S)=divq(ixO^S)+qface(ixO^S,idims)-qface(ixB^S,idims)
  end do
  divq(ixO^S)=divq(ixO^S)/block%dvolume(ixO^S)
end if

end subroutine divvector_s
