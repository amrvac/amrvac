!=============================================================================
! amrvacusr.t.jetCloudInteraction
!=============================================================================
! Deflection of a jet by a cloud
! Variation/Simplification of De Gouveia del Pino, E., ApJ 526, 862-873 (1999)
!   original paper: 3D SPH simulations of adiabatic versus cooling jets
!                          into gravitationally stratified isothermal cloud
!   our approximation: adiabatic HD, 2D to 3D, non-isothermal (pressure-matched)
!                          cloud, single tracer added to jet
!
! For setup in 2D use:
!$AMRVAC_DIR/setup.pl -d=22 -phi=0 -z=0 -g=16,16    -p=hd -eos=default -nf=1 -ndust=0 -u=nul -arch=default
! For setup in 3D use:
!$AMRVAC_DIR/setup.pl -d=33 -phi=0 -z=0 -g=16,16,16 -p=hd -eos=default -nf=1 -ndust=0 -u=nul -arch=default

INCLUDE:amrvacnul/specialimpl.t
INCLUDE:amrvacnul/speciallog.t
INCLUDE:amrvacnul/specialsource.t
INCLUDE:amrvacnul/usrflags.t
!=============================================================================
subroutine initglobaldata_usr

use mod_global_parameters
!-----------------------------------------------------------------------------
eqpar(gamma_)=1.66666667d0

! jet to cloud density ratio parameter
eqpar(beta_) = 0.04d0
! jet to ambient density contrast
eqpar(eta_)  = 3.0d0
! ambient sound speed (normalized)
eqpar(ca_)   = one
! jet Mach number
eqpar(Ma_)   = 12.0d0
! cloud to jet radii ratio
eqpar(rc_)   = 1.5d0

end subroutine initglobaldata_usr
!=============================================================================
subroutine initonegrid_usr(ixG^L,ix^L,w,x)

! initialize one grid

use mod_global_parameters

integer, intent(in) :: ixG^L,ix^L
double precision, intent(in) :: x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)

double precision:: rinlet(ixG^T)
double precision:: rcloud(ixG^T)

logical :: patchw(ixG^T)
DOUBLE PRECISION :: xc,yc,zc,sigma
!----------------------------------------------------------------------------
{^IFTWOD
    rinlet(ix^S)=abs(x(ix^S,2))
}
{^IFTHREED
    rinlet(ix^S)=dsqrt(x(ix^S,2)**2+x(ix^S,3)**2)
}

where(rinlet(ix^S)<=1.0d0 .and. abs(x(ix^S,1)-xprobmin1)<=2.5d0)
    !=== Set Jet ===!
    w(ix^S,rho_) = 1.0d0
    w(ix^S,v1_)  = eqpar(Ma_)*eqpar(ca_)
    w(ix^S,v2_)  = zero
{^IFTHREED
    w(ix^S,v3_)=zero
}
    w(ix^S,p_)   = eqpar(ca_)**2/(eqpar(gamma_)*eqpar(eta_))
{#IFDEF TRACER
    w(ix^S,tr1_) = 100.0d0
}
elsewhere
    !=== Set Ambient ===!
    w(ix^S,rho_) = 1.0d0/eqpar(eta_)
    w(ix^S,v1_)=zero
    w(ix^S,v2_)=zero
{^IFTHREED
    w(ix^S,v3_)=zero
}
    w(ix^S,p_)   = eqpar(ca_)**2/(eqpar(gamma_)*eqpar(eta_))
{#IFDEF TRACER
    w(ix^S,tr1_) = zero
}
endwhere


!=== Set Cloud ===!
! cloud coordinates xc,yc,zc
xc = zero;yc = 1.2d0;zc = 0.0d0;sigma=0.75d0*eqpar(rc_)

rcloud(ix^S)=(x(ix^S,1)-xc)**2+(x(ix^S,2)-yc)**2
{^IFTHREED
    rcloud(ix^S)= rcloud(ix^S) + (x(ix^S,3)-zc)**2
}

where(dsqrt(rcloud(ix^S))<=eqpar(rc_))
    w(ix^S,rho_) = 1.0d0/eqpar(eta_) + (1.0d0/(eqpar(beta_)**2))*dexp(-rcloud(ix^S)/(sigma*sigma))
endwhere


patchw(ixG^S)=.false.
call conserve(ixG^L,ix^L,w,x,patchw)

end subroutine initonegrid_usr
!=============================================================================
subroutine specialbound_usr(qt,ixG^L,ixO^L,iw,iB,w,x)

! special boundary types, user defined
! user must assign conservative variables in bounderies

use mod_global_parameters

integer, intent(in) :: ixG^L, ixO^L, iw, iB
double precision, intent(in) :: qt, x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)

logical :: patchw(ixG^T)
double precision:: rinlet(ixG^T)
!----------------------------------------------------------------------------


select case(iB)
 case(1)
    ! === Left boundary ===!

  {^IFTWOD
    rinlet(ixO^S)=abs(x(ixO^S,2))
  }
  {^IFTHREED
    rinlet(ixO^S)=dsqrt(x(ixO^S,2)**2+x(ixO^S,3)**2)
  }

  where(rinlet(ixO^S)<1.0d0)
      w(ixO^S,rho_) = one
      w(ixO^S,v1_)  = eqpar(Ma_)*eqpar(ca_)
      w(ixO^S,v2_)  = zero
      w(ixO^S,p_)   = eqpar(ca_)**2/(eqpar(gamma_)*eqpar(eta_))
{#IFDEF TRACER
      w(ixO^S,tr1_) = 100.0d0
}
  elsewhere
      w(ixO^S,rho_) = 1.0d0/eqpar(eta_)
      w(ixO^S,v1_)  = zero
      w(ixO^S,v2_)  = zero
      w(ixO^S,p_)   = eqpar(ca_)**2/(eqpar(gamma_)*eqpar(eta_))
{#IFDEF TRACER
      w(ixO^S,tr1_) = zero
}
  endwhere
  {^IFTHREED
      w(ixO^S,v3_)=zero
  }

  patchw(ixO^S)=.false.
  call conserve(ixG^L,ixO^L,w,x,patchw)
 case default
   call mpistop('boundary not defined')
end select

end subroutine specialbound_usr
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
! amrvacusr.t.jetCloudInteraction
!=============================================================================
