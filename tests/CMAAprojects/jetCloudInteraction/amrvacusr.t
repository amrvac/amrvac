!=============================================================================
! amrvacusr.t.jetDeflection
!=============================================================================
! Deflection of a jet by a cloud
! For setup in 2D use: 
!$AMRVAC_DIR/setup.pl -d=22 -phi=0 -z=0 -g=26,26 -p=hd -eos=default -nf=1 -ndust=0 -u=nul -arch=default

!INCLUDE:amrvacnul/specialini.t
INCLUDE:amrvacnul/specialimpl.t
INCLUDE:amrvacnul/speciallog.t
!INCLUDE:amrvacnul/specialbound.t
INCLUDE:amrvacnul/specialsource.t
INCLUDE:amrvacnul/usrflags.t
!=============================================================================
subroutine initglobaldata_usr

include 'amrvacdef.f'
double precision:: r(0:^NDS)
integer:: i
logical, save :: first=.true.
!-----------------------------------------------------------------------------
eqpar(gamma_)=1.66666667d0

eqpar(beta_) = 0.04d0
eqpar(eta_)  = 3.0d0
eqpar(ca_)   = one
eqpar(Ma_)   = 12.0d0

end subroutine initglobaldata_usr
!=============================================================================
subroutine initonegrid_usr(ixG^L,ix^L,w,x)

! initialize one grid

include 'amrvacdef.f'

integer, intent(in) :: ixG^L,ix^L
double precision, intent(in) :: x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)

double precision:: rinlet(ixG^T)
double precision:: rcloud(ixG^T),pleft,rholeft,vleft,Prat,alfa,cleft
double precision:: pright,rhoright,vright,deltax

logical :: patchw(ixG^T)
logical, save :: first=.true.
DOUBLE PRECISION :: mhydro = 1.67e-24
DOUBLE PRECISION :: kboltz = 1.38e-16
DOUBLE PRECISION :: Tism   = 200.0e+0
DOUBLE PRECISION :: xc,yc,zc,rc,sigma
!----------------------------------------------------------------------------
oktest = index(teststr,'initonegrid_usr')>=1
if (oktest) write(unitterm,*) ' === initonegrid_usr  (in ) : ', &
                'ixG^L : ',ixG^L

{^IFTWOD
    rinlet(ix^S)=abs(x(ix^S,2))
}
{^IFTHREED
    rinlet(ix^S)=dsqrt(x(ix^S,2)**2+x(ix^S,3)**2)
}

where(rinlet(ixG^S)<=1.0d0 .and. abs(x(ixG^S,1)-xprobmin1)<=2.5d0)
    !=== Set Jet ===!
    w(ixG^S,rho_) = 1.0d0
    w(ixG^S,v1_)  = eqpar(Ma_)*eqpar(ca_)
    w(ixG^S,v2_)  = zero
{^IFTHREED
    w(ixG^S,v3_)=zero
}
    w(ixG^S,tr1_) = 100.0d0
    w(ixG^S,p_)   = eqpar(ca_)**2/(eqpar(gamma_)*eqpar(eta_))
elsewhere
    !=== Set Ambient ===!
    w(ixG^S,rho_) = 1.0d0/eqpar(eta_)
    w(ixG^S,v1_)=zero
    w(ixG^S,v2_)=zero
{^IFTHREED
    w(ixG^S,v3_)=zero
}
    w(ixG^S,p_)   = eqpar(ca_)**2/(eqpar(gamma_)*eqpar(eta_))
    w(ixG^S,tr1_) = zero
endwhere



!=== Set Cloud ===!
! cloud coordinates xc,yc,zc
xc = zero;yc = 1.2d0;zc = 0.0d0;rc=1.0d0;sigma=0.75d0

rcloud(ix^S)=(x(ix^S,1)-xc)**2+(x(ix^S,2)-yc)**2
{^IFTHREED
    rcloud(ix^S)= rcloud(ix^S) + (x(ix^S,3)-zc)**2
}
where(dsqrt(rcloud(ix^S))<=3.0*rc)
    w(ix^S,rho_) = 1.0d0/eqpar(eta_) + (1.0d0/(eqpar(beta_)**2))*dexp(-rcloud(ix^S)/(sigma*sigma))
endwhere




patchw(ixG^S)=.false.
call conserve(ixG^L,ix^L,w,x,patchw)

end subroutine initonegrid_usr
!=============================================================================
subroutine specialbound_usr(qt,ixG^L,ixO^L,iw,iB,w,x)

! special boundary types, user defined
! user must assign conservative variables in bounderies

include 'amrvacdef.f'

integer, intent(in) :: ixG^L, ixO^L, iw, iB
double precision, intent(in) :: qt, x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)
integer :: ixI^L, ix2

logical :: patchw(ixG^T)
double precision:: rinlet(ixG^S)
!----------------------------------------------------------------------------


select case(iB)
 case(1)
    ! === Left boundary ===!
  
  {^IFTWOD
    rinlet(ixG^S)=abs(x(ixG^S,2))
  }
  {^IFTHREED
    rinlet(ixG^S)=dsqrt(x(ixG^S,2)**2+x(ixG^S,3)**2)
  }
		
  where(rinlet(ixO^S)<1.0d0)
      w(ixO^S,rho_) = one
      w(ixO^S,v1_)  = eqpar(Ma_)*eqpar(ca_)
      w(ixO^S,v2_)  = zero
      w(ixO^S,tr1_) = 100.0d0
      w(ixO^S,p_)   = eqpar(ca_)**2/(eqpar(gamma_)*eqpar(eta_))
  elsewhere
      w(ixO^S,rho_) = 1.0d0/eqpar(eta_)
      w(ixO^S,v1_)  = zero
      w(ixO^S,v2_)  = zero
      w(ixO^S,tr1_) = zero
      w(ixO^S,p_)   = eqpar(ca_)**2/(eqpar(gamma_)*eqpar(eta_))
  endwhere
  {^IFTHREED
      w(ixO^S,v3_)=zero
  }
  
  call conserve(ixG^L,ixO^L,w,x,patchfalse)
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
! amrvacusr.t.shockcloud
!=============================================================================
