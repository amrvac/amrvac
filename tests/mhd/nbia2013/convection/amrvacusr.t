!=============================================================================
! amrvacusr.t.mhdconvection
!=============================================================================
INCLUDE:amrvacmodules/heatconduct.t
INCLUDE:amrvacmodules/gravity.t
INCLUDE:amrvacmodules/viscosity.t

INCLUDE:amrvacnul/usrflags.t
!=============================================================================
subroutine initglobaldata_usr

use mod_global_parameters

integer:: mpoly
double precision:: zeta0,rhat,sigma,zz0,qchand
double precision:: gamma, qchi, qmpoly, eta2
!-----------------------------------------------------------------------------

! This setup relates to Hurlburt and Toomre
! and calculates the eqpar array from problem-specific input parameters
! ---------------------------------------------------------------------
!
! problem specific input parameters for
! Hurlburt and Toomre magnetoconvection problem
!
! we need (apart from the aspect ratio which is set in par-file):
! i)equilibrium parameters: (1) ratio of specific heats (gamma)
! ------------------------- (2) the polytropic index mpoly (m)
!                            --> sets gravity g
!                           (3) the density contrast chi
!                            --> sets top temperature zz0 (z0)
!                           (4) the Chandrasekhar number qchand (Q)
!                            --> initial field strength B (need nu and eta)
! ii)parameters setting the dissipative coefficients:
! ---------------------------------------------------
!    (1) Prandtl number sigma (viscous/thermal conduction)
!    (2) the Rayleigh number at half-depth rhat (degree of instability)
!    (3) the magnetic Prandtl number zeta0
!      (magnetic diffusion/thermal conduction)
!
! hence fixing   i) gamma, mpoly, chi, qchand
!               ii) sigma, rhat, zeta0
!
! allows calculation of all equation parameters, namely
! the general purpose ones:
! ------------------------
!  eqpar(gamma_) eqpar(eta_) eqpar(kappa_) eqpar(grav1.2_) eqpar(mu_)
!
! and the problem specific one:
! ----------------------------
!  eqpar(temptop_) for setting the top temperature
!  eqpar(bstr_) for setting the B field
!

! equilibrium parameters
gamma=1.66666667d0
mpoly=1
qchi=11.0d0
qchand=72.0d0

qmpoly=dble(mpoly)
zz0=one/(qchi**(one/qmpoly)-one)

if (mype==0) then
  write(*,*)'gamma, mpoly, qchi,qchand'
  write(*,*)gamma, mpoly, qchi, qchand
  write(*,*)'Deduced top dimensionless temperature z0=',zz0
endif

eqpar(temptop_)=zz0
eqpar(gamma_)=gamma
eqpar(grav1_)= zero
eqpar(grav2_)=-(qmpoly+one)
{^IFTHREED
eqpar(grav3_)=zero
}

! dissipative parameters
sigma=1.0d0
rhat=1.0d5
zeta0=0.25d0

if (mype==0) then
  write(*,*) 'sigma,rhat,zeta0:'
  write(*,*) sigma,rhat,zeta0
endif

eta2=(qmpoly+one)*(gamma/(gamma-1)-(qmpoly+one))*((gamma-one)/gamma)*&
     ((zz0+half)**(two*qmpoly-one)/zz0**(two*qmpoly))*&
     (zeta0**2/(sigma*rhat))

if(eta2<=smalldouble)then
   if(mype==0) write(*,*) 'eta2=',eta2
   call mpistop("Negative or too small value for eta**2")
endif

eqpar(eta_)=dsqrt(eta2)
eqpar(mu_)=eqpar(eta_)*sigma/zeta0
eqpar(kappa_)=(gamma/(gamma-one))*eqpar(eta_)/zeta0
eqpar(bstr_)=dsqrt(qchand*eqpar(mu_)*eqpar(eta_))

if (mype==0) then
  write(*,*)'dimensionless values for dissipative coefficients:'
  write(*,*)'resistivity          eta=',eqpar(eta_)
  write(*,*)'viscosity             mu=',eqpar(mu_)
  write(*,*)'thermal conduction kappa=',eqpar(kappa_)
  write(*,*)'dimensionless magnetic field strength:',eqpar(bstr_)
endif

end subroutine initglobaldata_usr
!=============================================================================
subroutine initonegrid_usr(ixG^L,ix^L,w,x)

! initialize one grid

use mod_global_parameters

integer, intent(in) :: ixG^L,ix^L
double precision, intent(in) :: x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)

double precision :: dvx,dvy,dvz,nkx,nky,nkz,zz0,qmpoly

logical :: patchw(ixG^T)
logical, save :: first=.true.
!----------------------------------------------------------------------------

! velocity perturbation parameters
dvx=1.0d-4
dvy=1.0d-4
nkx=two*dpi/3.0d0
nky=two*dpi
{^IFTHREED
dvz=1.0d-4
nkz=two*two*dpi
}

if(first)then
  if(mype==0) then
     write(*,*)'Simulating mhd convection: Hurlburt-Toomre'
     write(*,*) 'velocity perturbation: dvx,dvy,nkx,nky'
     write(*,*) dvx,dvy,nkx,nky
     {^IFTHREED
     write(*,*) '3D setup!'
     write(*,*) dvz,nkz
     }
  endif
  first=.false.
endif

zz0=eqpar(temptop_)
qmpoly=-one-eqpar(grav2_)

w(ix^S,rho_)=((zz0+one-x(ix^S,2))/zz0)**qmpoly
w(ix^S,b1_)=zero
w(ix^S,b2_)=eqpar(bstr_)
{^IFTHREED
w(ix^S,b3_)=zero
}

w(ix^S,p_)= zz0*(((zz0+one-x(ix^S,2))/zz0)**(qmpoly+one))

w(ix^S,v1_)=dvx*dsin(x(ix^S,1)*nkx)*dsin(x(ix^S,2)*nky){^IFTHREED *dsin(x(ix^S,3)*nkz)}
w(ix^S,v2_)=dvy*dsin(x(ix^S,1)*nkx)*dsin(x(ix^S,2)*nky){^IFTHREED *dsin(x(ix^S,3)*nkz)}
{^IFTHREED
w(ix^S,v3_)=dvz*dsin(x(ix^S,1)*nkx)*dsin(x(ix^S,2)*nky)*dsin(x(ix^S,3)*nkz)
}

patchw(ixG^S)=.false.
call conserve(ixG^L,ix^L,w,x,patchw)

end subroutine initonegrid_usr
!=============================================================================
subroutine specialbound_usr(qt,ixG^L,ixO^L,iw,iB,w,x)

! special boundary types, user defined

use mod_global_parameters

integer, intent(in) :: ixO^L, iw, iB, ixG^L
double precision, intent(in) :: qt, x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)

double precision :: tempb(ixG^T)
double precision :: delydelx,delydelz
integer :: ix^D,ixIM^L

logical :: patchw(ixG^T)
!----------------------------------------------------------------------------

select case(iB)
case(3)
  ! special bottom boundary
  ! ensure fixed temperature gradient in ghost layers
  ! use asymm/symm for B-components (ensuring vertical field)
  ! use divB correction from central difference for vertical field
  ! use symm/asymm for v-components (ensuring no flow-through)
  ! density profile: use symmetry

  {^IFTWOD
  ! in dixB rows above the bottom boundary: switch to primitive
  ixIMmin2=ixOmax2+1;ixIMmax2=ixOmax2+dixB;
  ixIMmin1=ixOmin1;ixIMmax1=ixOmax1;
  patchw(ixIM^S)=.false.
  ! primitiven (instead of primitive) to avoid changing internal values
  call primitiven(ixG^L,ixIM^L,w,patchw)
  do ix2=ixOmax2,ixOmin2,-1
     w(ixOmin1:ixOmax1,ix2,rho_)= w(ixOmin1:ixOmax1,2*ixOmax2+1-ix2,rho_)
     w(ixOmin1:ixOmax1,ix2,v1_) = w(ixOmin1:ixOmax1,2*ixOmax2+1-ix2,v1_)
     w(ixOmin1:ixOmax1,ix2,v2_) =-w(ixOmin1:ixOmax1,2*ixOmax2+1-ix2,v2_)
     w(ixOmin1:ixOmax1,ix2,b1_) =-w(ixOmin1:ixOmax1,2*ixOmax2+1-ix2,b1_)
     w(ixOmin1:ixOmax1,ix2,b2_) = w(ixOmin1:ixOmax1,2*ixOmax2+1-ix2,b2_)
  enddo
  ! fill temperature array: extrapolate linearly with fixed dT/dy=-1, from bottom row
  do ix2=ixOmin2,ixOmax2
    tempb(ixOmin1:ixOmax1,ix2)=w(ixOmin1:ixOmax1,ixOmax2+1,p_)/w(ixOmin1:ixOmax1,ixOmax2+1,rho_) &
                             -(x(ixOmin1:ixOmax1,ix2,2)-x(ixOmin1:ixOmax1,ixOmax2+1,2))
  enddo
  ! determine pressure
  w(ixO^S,p_)=w(ixO^S,rho_)*tempb(ixO^S)
  delydelx=(x(ixOmin1,ixOmin2+1,2)-x(ixOmin1,ixOmin2,2))/(x(ixOmin1+1,ixOmin2,1)-x(ixOmin1,ixOmin2,1))
  do ix2=ixOmax2,ixOmin2,-1
     do ix1=ixOmin1+1,ixOmax1-1
         w(ix1,ix2,b2_)=w(ix1,ix2+2,b2_)+delydelx*(w(ix1+1,ix2+1,b1_)-w(ix1-1,ix2+1,b1_))
     enddo
  enddo
  }
  {^IFTHREED
  ! in dixB rows above the bottom boundary: switch to primitive
  ixIMmin2=ixOmax2+1;ixIMmax2=ixOmax2+dixB;
  ixIMmin1=ixOmin1;ixIMmax1=ixOmax1;
  ixIMmin3=ixOmin3;ixIMmax3=ixOmax3;
  patchw(ixIM^S)=.false.
  ! primitiven (instead of primitive) to avoid changing internal values
  call primitiven(ixG^L,ixIM^L,w,patchw)
  do ix2=ixOmax2,ixOmin2,-1
     w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,rho_)= w(ixOmin1:ixOmax1,2*ixOmax2+1-ix2,ixOmin3:ixOmax3,rho_)
     w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,v1_) = w(ixOmin1:ixOmax1,2*ixOmax2+1-ix2,ixOmin3:ixOmax3,v1_)
     w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,v2_) =-w(ixOmin1:ixOmax1,2*ixOmax2+1-ix2,ixOmin3:ixOmax3,v2_)
     w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,v3_) = w(ixOmin1:ixOmax1,2*ixOmax2+1-ix2,ixOmin3:ixOmax3,v3_)
     w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,b1_) =-w(ixOmin1:ixOmax1,2*ixOmax2+1-ix2,ixOmin3:ixOmax3,b1_)
     w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,b2_) = w(ixOmin1:ixOmax1,2*ixOmax2+1-ix2,ixOmin3:ixOmax3,b2_)
     w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,b3_) =-w(ixOmin1:ixOmax1,2*ixOmax2+1-ix2,ixOmin3:ixOmax3,b3_)
  enddo
  ! fill temperature array: extrapolate linearly with fixed dT/dy=-1, from bottom row
  do ix2=ixOmin2,ixOmax2
    tempb(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3)=w(ixOmin1:ixOmax1,ixOmax2+1,ixOmin3:ixOmax3,p_) &
                                              /w(ixOmin1:ixOmax1,ixOmax2+1,ixOmin3:ixOmax3,rho_) &
           -(x(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,2)-x(ixOmin1:ixOmax1,ixOmax2+1,ixOmin3:ixOmax3,2))
  enddo
  ! determine pressure
  w(ixO^S,p_)=w(ixO^S,rho_)*tempb(ixO^S)
  delydelx=(x(ixOmin1,ixOmin2+1,ixOmin3,2)-x(ixOmin1,ixOmin2,ixOmin3,2)) &
          /(x(ixOmin1+1,ixOmin2,ixOmin3,1)-x(ixOmin1,ixOmin2,ixOmin3,1))
  delydelz=(x(ixOmin1,ixOmin2+1,ixOmin3,2)-x(ixOmin1,ixOmin2,ixOmin3,2)) &
          /(x(ixOmin1,ixOmin2,ixOmin3+1,3)-x(ixOmin1,ixOmin2,ixOmin3,3))
  do ix2=ixOmax2,ixOmin2,-1
     do ix1=ixOmin1+1,ixOmax1-1
     do ix3=ixOmin3+1,ixOmax3-1
         w(ix1,ix2,ix3,b2_)=w(ix1,ix2+2,ix3,b2_)+delydelx*(w(ix1+1,ix2+1,ix3,b1_)-w(ix1-1,ix2+1,ix3,b1_)) &
                                                +delydelz*(w(ix1,ix2+1,ix3+1,b3_)-w(ix1,ix2+1,ix3-1,b3_))
     enddo
     enddo
  enddo
  }
  ! now reset the inner mesh values to conservative
  ! use conserven to avoid changing values
  call conserven(ixG^L,ixIM^L,w,patchw)
  ! now switch to conservative in full bottom ghost layer
  patchw(ixO^S)=.false.
  call conserve(ixG^L,ixO^L,w,x,patchw)

case(4)
  ! special top boundary
  ! ensure the fixed temperature in ghost layers
  ! use asymm/symm for B-components (ensuring vertical field)
  ! use divB correction from central difference for vertical field
  ! use symm/asymm for v-components (ensuring no flow-through)
  ! density profile: use symmetry

  {^IFTWOD
  ! in dixB rows below top boundary: switch to primitive
  ixIMmin2=ixOmin2-dixB;ixIMmax2=ixOmin2-1;
  ixIMmin1=ixOmin1;ixIMmax1=ixOmax1;
  patchw(ixIM^S)=.false.
  ! primitiven (instead of primitive) to avoid changing internal values
  call primitiven(ixG^L,ixIM^L,w,patchw)
  do ix2=ixOmin2,ixOmax2,+1
      w(ixOmin1:ixOmax1,ix2,rho_)= w(ixOmin1:ixOmax1,2*ixOmin2-ix2-1,rho_)
      w(ixOmin1:ixOmax1,ix2,v1_) = w(ixOmin1:ixOmax1,2*ixOmin2-ix2-1,v1_)
      w(ixOmin1:ixOmax1,ix2,v2_) =-w(ixOmin1:ixOmax1,2*ixOmin2-ix2-1,v2_)
      w(ixOmin1:ixOmax1,ix2,b1_) =-w(ixOmin1:ixOmax1,2*ixOmin2-ix2-1,b1_)
      w(ixOmin1:ixOmax1,ix2,b2_) = w(ixOmin1:ixOmax1,2*ixOmin2-ix2-1,b2_)
  enddo
  w(ixO^S,p_)=w(ixO^S,rho_)*eqpar(temptop_)
  delydelx=(x(ixOmin1,ixOmin2+1,2)-x(ixOmin1,ixOmin2,2))/(x(ixOmin1+1,ixOmin2,1)-x(ixOmin1,ixOmin2,1))
  do ix2=ixOmin2,ixOmax2,+1
     do ix1=ixOmin1+1,ixOmax1-1
         w(ix1,ix2,b2_)=w(ix1,ix2-2,b2_)-delydelx*(w(ix1+1,ix2-1,b1_)-w(ix1-1,ix2-1,b1_))
     enddo
  enddo
  }
  {^IFTHREED
  ! in dixB rows below top boundary: switch to primitive
  ixIMmin2=ixOmin2-dixB;ixIMmax2=ixOmin2-1;
  ixIMmin1=ixOmin1;ixIMmax1=ixOmax1;
  ixIMmin3=ixOmin3;ixIMmax3=ixOmax3;
  patchw(ixIM^S)=.false.
  ! primitiven (instead of primitive) to avoid changing internal values
  call primitiven(ixG^L,ixIM^L,w,patchw)
  do ix2=ixOmin2,ixOmax2,+1
      w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,rho_)= w(ixOmin1:ixOmax1,2*ixOmin2-ix2-1,ixOmin3:ixOmax3,rho_)
      w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,v1_) = w(ixOmin1:ixOmax1,2*ixOmin2-ix2-1,ixOmin3:ixOmax3,v1_)
      w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,v2_) =-w(ixOmin1:ixOmax1,2*ixOmin2-ix2-1,ixOmin3:ixOmax3,v2_)
      w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,v3_) = w(ixOmin1:ixOmax1,2*ixOmin2-ix2-1,ixOmin3:ixOmax3,v3_)
      w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,b1_) =-w(ixOmin1:ixOmax1,2*ixOmin2-ix2-1,ixOmin3:ixOmax3,b1_)
      w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,b2_) = w(ixOmin1:ixOmax1,2*ixOmin2-ix2-1,ixOmin3:ixOmax3,b2_)
      w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,b3_) =-w(ixOmin1:ixOmax1,2*ixOmin2-ix2-1,ixOmin3:ixOmax3,b3_)
  enddo
  w(ixO^S,p_)=w(ixO^S,rho_)*eqpar(temptop_)
  delydelx=(x(ixOmin1,ixOmin2+1,ixOmin3,2)-x(ixOmin1,ixOmin2,ixOmin3,2)) &
           /(x(ixOmin1+1,ixOmin2,ixOmin3,1)-x(ixOmin1,ixOmin2,ixOmin3,1))
  delydelz=(x(ixOmin1,ixOmin2+1,ixOmin3,2)-x(ixOmin1,ixOmin2,ixOmin3,2)) &
           /(x(ixOmin1,ixOmin2,ixOmin3+1,3)-x(ixOmin1,ixOmin2,ixOmin3,3))
  do ix2=ixOmin2,ixOmax2,+1
     do ix1=ixOmin1+1,ixOmax1-1
     do ix3=ixOmin3+1,ixOmax3-1
         w(ix1,ix2,ix3,b2_)=w(ix1,ix2-2,ix3,b2_)-delydelx*(w(ix1+1,ix2-1,ix3,b1_)-w(ix1-1,ix2-1,ix3,b1_)) &
                                                -delydelz*(w(ix1,ix2-1,ix3+1,b3_)-w(ix1,ix2-1,ix3-1,b3_))
     enddo
     enddo
  enddo
  }
  ! now reset the inner mesh values to conservative
  ! use conserven to avoid changing values
  call conserven(ixG^L,ixIM^L,w,patchw)
  ! now switch to conservative in full bottom layer
  patchw(ixO^S)=.false.
  call conserve(ixG^L,ixO^L,w,x,patchw)

case default
   call mpistop("Special boundary is not defined for this region")
end select

end subroutine specialbound_usr
!=============================================================================
subroutine specialsource(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L, iw^LIM
double precision, intent(in) :: qdt, qtC, qt, wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
double precision, intent(inout) :: w(ixI^S,1:nw)

double precision :: lQgrid(ixG^T),bQgrid(ixG^T)

! coordinate info can be used in region ixO

integer :: iw
!-----------------------------------------------------------------------------
call addsource_grav(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

if(abs(eqpar(mu_))>smalldouble) call addsource_visc(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

if(abs(eqpar(kappa_))>smalldouble) then
  if(.not.sourceimpl)then
    call addsource_heatconduct_mhd(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)
  endif
endif

end subroutine specialsource
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

if(abs(eqpar(mu_))>smalldouble) call getdt_visc(w,ixG^L,ix^L,dtnew,dx^D,x)

if(abs(eqpar(kappa_))>smalldouble) then
   if(.not.sourceimpl)then
      call getdt_heatconduct_mhd(w,ixG^L,ix^L,dtnew,dx^D,x)
   endif
endif

end subroutine getdt_special
!=============================================================================
subroutine specialsource_impl(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L, iw^LIM
double precision, intent(in) :: qdt, qtC, qt, x(ixI^S,1:ndim)
double precision, intent(inout) :: w(ixI^S,1:nw), wCT(ixI^S,1:nw)
!-----------------------------------------------------------------------------

call addsource_heatconduct_mhd(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

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

call getdt_heatconduct_mhd(w,ixG^L,ix^L,dtnew,dx^D,x)

end subroutine getdt_impl
!=============================================================================
subroutine specialeta(w,ixI^L,ix^L,idirmin,x,current,eta)

! Set the common "eta" array for resistive MHD based on w or the
! "current" variable which has components between idirmin and 3.

use mod_global_parameters

integer, intent(in) :: ixI^L, ix^L, idirmin
double precision, intent(in) :: w(ixI^S,nw), x(ixI^S,1:ndim)

double precision :: current(ixG^T,7-2*ndir:3), eta(ixG^T)
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

!double precision :: tmp(ixG^T),tmp1(ixG^T),tmp2(ixG^T)
!double precision :: qvec(ixG^T,1:ndir),curlvec(ixG^T,1:ndir)
!integer :: idirmin
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

double precision :: tmp(ixG^T)
double precision :: divb(ixG^T)
double precision :: qvec(ixG^T,1:ndir),curlvec(ixG^T,7-2*ndir:3)
integer                            :: idirmin
logical          :: patchw(ixG^T)
!-----------------------------------------------------------------------------
! output Te
if(saveprim)then
   tmp(ixO^S)=w(ixO^S,p_)
 else
   call getpthermal(w,x,ixI^L,ixO^L,tmp)
endif
w(ixO^S,nw+1)=tmp(ixO^S)/w(ixO^S,rho_)
! output B 
w(ixO^S,nw+2)=dsqrt(^C&w(ixO^S,b^C_)**2+)
! output divB1
call getdivb(w,ixI^L,ixO^L,divb)
w(ixO^S,nw+3)=divb(ixO^S)
! output the plasma beta p*2/B**2
w(ixO^S,nw+4)=tmp(ixO^S)*two/(^C&w(ixO^S,b^C_)**2+)
! store current
^C&qvec(ixI^S,^C)=w(ixI^S,b^C_);
call curlvector(qvec,ixI^L,ixO^L,curlvec,idirmin,3,ndir)
w(ixO^S,nw+5)=curlvec(ixO^S,3);

end subroutine specialvar_output
!=============================================================================
subroutine specialvarnames_output

! newly added variables need to be concatenated with the wnames/primnames string

use mod_global_parameters
!-----------------------------------------------------------------------------
oktest = index(teststr,'printlog')>=1

primnames=TRIM(primnames)//' '//'Te B divB beta j3'
wnames=   TRIM(wnames)//' '//'Te B divB beta j3'

end subroutine specialvarnames_output
!=============================================================================
subroutine specialset_B0(ixI^L,ixO^L,x,wB0)

! Here one can add a steady (time-independent) potential background field

use mod_global_parameters

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

use mod_global_parameters

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
subroutine userspecialconvert(qunitconvert)

use mod_global_parameters
!double precision, intent(inout) :: w(ixI^S,nw+nwauxio)
!integer, intent(in)           :: ixI^L,ixO^L

integer, intent(in) :: qunitconvert
character(len=20):: userconvert_type
!-----------------------------------------------------------------------------

end subroutine userspecialconvert
!=============================================================================
! amrvacusr.t.mhdconvection
!=============================================================================
