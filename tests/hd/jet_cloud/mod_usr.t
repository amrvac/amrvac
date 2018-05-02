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

module mod_usr
  use mod_hd

  implicit none

  ! jet to cloud density ratio parameter
  double precision :: beta_ = 0.04d0
  ! jet to ambient density contrast
  double precision :: eta_  = 3.0d0
  ! ambient sound speed (normalized)
  double precision :: ca_   = 1.0d0
  ! jet Mach number
  double precision :: Ma_   = 12.0d0
  ! cloud to jet radii ratio
  double precision :: rc_   = 1.5d0

contains

  subroutine usr_init()

    usr_init_one_grid => initonegrid_usr
    usr_special_bc => specialbound_usr

    call hd_activate()

  end subroutine usr_init

  ! initialize one grid
  subroutine initonegrid_usr(ixG^L,ix^L,w,x)

    integer, intent(in) :: ixG^L,ix^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    double precision:: rinlet(ixG^T)
    double precision:: rcloud(ixG^T)

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
       w(ix^S,mom(1))  = Ma_*ca_
       w(ix^S,mom(2))  = 0.0d0
       {^IFTHREED
       w(ix^S,mom(3))=0.0d0
       }
       w(ix^S,e_)   = ca_**2/(hd_gamma*eta_)
       w(ix^S, tracer(1)) = 100.0d0
    elsewhere
       !=== Set Ambient ===!
       w(ix^S,rho_) = 1.0d0/eta_
       w(ix^S,mom(1))=0.0d0
       w(ix^S,mom(2))=0.0d0
       {^IFTHREED
       w(ix^S,mom(3))=0.0d0
       }
       w(ix^S,e_)   = ca_**2/(hd_gamma*eta_)
       w(ix^S,tracer(1)) = 0.0d0
    endwhere


    !=== Set Cloud ===!
    ! cloud coordinates xc,yc,zc
    xc = 0.0d0;yc = 1.2d0;zc = 0.0d0;sigma=0.75d0*rc_

    rcloud(ix^S)=(x(ix^S,1)-xc)**2+(x(ix^S,2)-yc)**2
    {^IFTHREED
    rcloud(ix^S)= rcloud(ix^S) + (x(ix^S,3)-zc)**2
    }

    where(dsqrt(rcloud(ix^S))<=rc_)
       w(ix^S,rho_) = 1.0d0/eta_ + (1.0d0/(beta_**2))*dexp(-rcloud(ix^S)/(sigma*sigma))
    endwhere

    call hd_to_conserved(ixG^L,ix^L,w,x)

  end subroutine initonegrid_usr

  ! special boundary types, user defined
  ! user must assign conservative variables in bounderies
  subroutine specialbound_usr(qt,ixG^L,ixO^L,iB,w,x)
    integer, intent(in) :: ixG^L, ixO^L, iB
    double precision, intent(in) :: qt, x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    logical :: patchw(ixG^T)
    double precision:: rinlet(ixG^T)

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
          w(ixO^S,rho_) = 1.0d0
          w(ixO^S,mom(1))  = Ma_*ca_
          w(ixO^S,mom(2))  = 0.0d0
          w(ixO^S,e_)   = ca_**2/(hd_gamma*eta_)
          w(ixO^S,tracer(1)) = 100.0d0
       elsewhere
          w(ixO^S,rho_) = 1.0d0/eta_
          w(ixO^S,mom(1))  = 0.0d0
          w(ixO^S,mom(2))  = 0.0d0
          w(ixO^S,e_)   = ca_**2/(hd_gamma*eta_)
          w(ixO^S,tracer(1)) = 0.0d0
       endwhere
       {^IFTHREED
       w(ixO^S,mom(3))=0.0d0
       }

       call hd_to_conserved(ixG^L,ixO^L,w,x)
    case default
       call mpistop('boundary not defined')
    end select

  end subroutine specialbound_usr

end module mod_usr
