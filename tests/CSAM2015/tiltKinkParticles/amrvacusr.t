!=============================================================================
! MHD tilt 3D

INCLUDE:amrvacnul/specialimpl.t
INCLUDE:amrvacnul/usrflags.t
INCLUDE:amrvacnul/correctaux_usr.t
INCLUDE:amrvacmodules/integrate_particles.t
INCLUDE:amrvacmodules/handle_particles.t
!=============================================================================
subroutine specialbound_usr(qt,ixG^L,ixO^L,iw,iB,w,x)

! special boundary types, user defined

use mod_global_parameters

integer, intent(in) :: ixG^L, ixO^L, iw, iB
double precision, intent(in) :: qt, x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)

double precision:: delydelx,delxdely,delydelz,delxdelz
integer:: ixIM^L,ix^D
logical :: patchw(ixG^T)
!----------------------------------------------------------------------------

select case(iB)
! implementation of special boundaries
case(1)
      {^IFTHREED ixIMmin3=ixOmin3;ixIMmax3=ixOmax3;}
      ixIMmin2=ixOmin2;ixIMmax2=ixOmax2;
      ixIMmin1=ixOmax1+1;ixIMmax1=ixOmax1+1;
      patchw(ixIM^S)=.false.
      ! primitiven (instead of primitive) to avoid changing internal values
      call primitiven(ixG^L,ixIM^L,w,patchw)
      do ix1=ixOmax1,ixOmin1,-1
        {^IFTWOD
{#IFDEF TRACER
        w(ix1,ixOmin2:ixOmax2,tr1_)= w(ixOmax1+1,ixOmin2:ixOmax2,tr1_)
}
        w(ix1,ixOmin2:ixOmax2,rho_)= w(ixOmax1+1,ixOmin2:ixOmax2,rho_)
        w(ix1,ixOmin2:ixOmax2,p_)  = w(ixOmax1+1,ixOmin2:ixOmax2,p_)
        w(ix1,ixOmin2:ixOmax2,v1_) = w(ixOmax1+1,ixOmin2:ixOmax2,v1_)
        w(ix1,ixOmin2:ixOmax2,v2_) = w(ixOmax1+1,ixOmin2:ixOmax2,v2_)
        w(ix1,ixOmin2:ixOmax2,v3_) = w(ixOmax1+1,ixOmin2:ixOmax2,v3_)
        !w(ix1,ixOmin2:ixOmax2,b1_) = w(ixOmax1+1,ixOmin2:ixOmax2,b1_)
        !w(ix1,ixOmin2:ixOmax2,b2_) = w(ixOmax1+1,ixOmin2:ixOmax2,b2_)
        w(ix1,ixOmin2:ixOmax2,b1_) = two*x(ix1,ixOmin2:ixOmax2,1)*x(ix1,ixOmin2:ixOmax2,2) &
                                    /(x(ix1,ixOmin2:ixOmax2,1)**2+x(ix1,ixOmin2:ixOmax2,2)**2)**2
        w(ix1,ixOmin2:ixOmax2,b2_) = -one+ (x(ix1,ixOmin2:ixOmax2,2)**2-x(ix1,ixOmin2:ixOmax2,1)**2) &
                                          /(x(ix1,ixOmin2:ixOmax2,1)**2+x(ix1,ixOmin2:ixOmax2,2)**2)**2
        w(ix1,ixOmin2:ixOmax2,b3_) = eqpar(b3val_)
        }
        {^IFTHREED
{#IFDEF TRACER
        w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,tr1_)= w(ixOmax1+1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,tr1_)
}
        w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)= w(ixOmax1+1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)
        w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,p_)  = w(ixOmax1+1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,p_)
        w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,v1_) = w(ixOmax1+1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,v1_)
        w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,v2_) = w(ixOmax1+1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,v2_)
        w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,v3_) = w(ixOmax1+1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,v3_)
        !w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b1_) = w(ixOmax1+1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b1_)
        !w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b2_) = w(ixOmax1+1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b2_)
        w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b1_) = two*x(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)*x(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2) &
                                    /(x(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)**2+x(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)**2)**2
        w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b2_) = -one+ (x(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)**2-x(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)**2) &
                                          /(x(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)**2+x(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)**2)**2
        w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b3_) = eqpar(b3val_)
        }
      enddo
      {^IFTWOD
      delxdely=(x(ixOmin1+1,ixOmin2,1)-x(ixOmin1,ixOmin2,1))/(x(ixOmin1,ixOmin2+1,2)-x(ixOmin1,ixOmin2,2))
      do ix1=ixOmax1,ixOmin1,-1
        do ix2=ixOmin2+1,ixOmax2-1
         w(ix1,ix2,b1_)=w(ix1+2,ix2,b1_)+delxdely*(w(ix1+1,ix2+1,b2_)-w(ix1+1,ix2-1,b2_)) 
        enddo
      enddo
      }
      {^IFTHREED
      delxdely=(x(ixOmin1+1,ixOmin2,ixOmin3,1)-x(ixOmin1,ixOmin2,ixOmin3,1))/(x(ixOmin1,ixOmin2+1,ixOmin3,2)-x(ixOmin1,ixOmin2,ixOmin3,2))
      delxdelz=(x(ixOmin1+1,ixOmin2,ixOmin3,1)-x(ixOmin1,ixOmin2,ixOmin3,1))/(x(ixOmin1,ixOmin2,ixOmin3+1,3)-x(ixOmin1,ixOmin2,ixOmin3,3))
      do ix1=ixOmax1,ixOmin1,-1
        do ix2=ixOmin2+1,ixOmax2-1
        do ix3=ixOmin3+1,ixOmax3-1
         w(ix1,ix2,ix3,b1_)=w(ix1+2,ix2,ix3,b1_)+delxdely*(w(ix1+1,ix2+1,ix3,b2_)-w(ix1+1,ix2-1,ix3,b2_)) &
                                                +delxdelz*(w(ix1+1,ix2,ix3+1,b3_)-w(ix1+1,ix2,ix3-1,b3_))
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
case(2)
      {^IFTHREED ixIMmin3=ixOmin3;ixIMmax3=ixOmax3;}
      ixIMmin2=ixOmin2;ixIMmax2=ixOmax2;
      ixIMmin1=ixOmin1-1;ixIMmax1=ixOmin1-1;
      patchw(ixIM^S)=.false.
      ! primitiven (instead of primitive) to avoid changing internal values
      call primitiven(ixG^L,ixIM^L,w,patchw)
      do ix1=ixOmin1,ixOmax1,+1
        {^IFTWOD
{#IFDEF TRACER
        w(ix1,ixOmin2:ixOmax2,tr1_)= w(ixOmin1-1,ixOmin2:ixOmax2,tr1_)
}
        w(ix1,ixOmin2:ixOmax2,rho_)= w(ixOmin1-1,ixOmin2:ixOmax2,rho_)
        w(ix1,ixOmin2:ixOmax2,p_)  = w(ixOmin1-1,ixOmin2:ixOmax2,p_)
        w(ix1,ixOmin2:ixOmax2,v1_) = w(ixOmin1-1,ixOmin2:ixOmax2,v1_)
        w(ix1,ixOmin2:ixOmax2,v2_) = w(ixOmin1-1,ixOmin2:ixOmax2,v2_)
        w(ix1,ixOmin2:ixOmax2,v3_) = w(ixOmin1-1,ixOmin2:ixOmax2,v3_)
        !w(ix1,ixOmin2:ixOmax2,b1_) = w(ixOmin1-1,ixOmin2:ixOmax2,b1_)
        !w(ix1,ixOmin2:ixOmax2,b2_) = w(ixOmin1-1,ixOmin2:ixOmax2,b2_)
        w(ix1,ixOmin2:ixOmax2,b1_) = two*x(ix1,ixOmin2:ixOmax2,1)*x(ix1,ixOmin2:ixOmax2,2) &
                                    /(x(ix1,ixOmin2:ixOmax2,1)**2+x(ix1,ixOmin2:ixOmax2,2)**2)**2
        w(ix1,ixOmin2:ixOmax2,b2_) = -one+ (x(ix1,ixOmin2:ixOmax2,2)**2-x(ix1,ixOmin2:ixOmax2,1)**2) &
                                          /(x(ix1,ixOmin2:ixOmax2,1)**2+x(ix1,ixOmin2:ixOmax2,2)**2)**2
        w(ix1,ixOmin2:ixOmax2,b3_) = eqpar(b3val_)
        }
        {^IFTHREED
{#IFDEF TRACER
        w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,tr1_)= w(ixOmin1-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,tr1_)
}
        w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)= w(ixOmin1-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)
        w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,p_)  = w(ixOmin1-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,p_)
        w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,v1_) = w(ixOmin1-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,v1_)
        w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,v2_) = w(ixOmin1-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,v2_)
        w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,v3_) = w(ixOmin1-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,v3_)
        !w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b1_) = w(ixOmin1-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b1_)
        !w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b2_) = w(ixOmin1-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b2_)
        w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b1_) = two*x(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)*x(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2) &
                                    /(x(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)**2+x(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)**2)**2
        w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b2_) = -one+ (x(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)**2-x(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)**2) &
                                          /(x(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)**2+x(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)**2)**2
        w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b3_) = eqpar(b3val_)
        }
      enddo
      {^IFTWOD
      delxdely=(x(ixOmin1+1,ixOmin2,1)-x(ixOmin1,ixOmin2,1))/(x(ixOmin1,ixOmin2+1,2)-x(ixOmin1,ixOmin2,2))
      do ix1=ixOmin1,ixOmax1,+1
        do ix2=ixOmin2+1,ixOmax2-1
         w(ix1,ix2,b1_)=w(ix1-2,ix2,b1_)-delxdely*(w(ix1-1,ix2+1,b2_)-w(ix1-1,ix2-1,b2_)) 
        enddo
      enddo
      }
      {^IFTHREED
      delxdely=(x(ixOmin1+1,ixOmin2,ixOmin3,1)-x(ixOmin1,ixOmin2,ixOmin3,1))/(x(ixOmin1,ixOmin2+1,ixOmin3,2)-x(ixOmin1,ixOmin2,ixOmin3,2))
      delxdelz=(x(ixOmin1+1,ixOmin2,ixOmin3,1)-x(ixOmin1,ixOmin2,ixOmin3,1))/(x(ixOmin1,ixOmin2,ixOmin3+1,3)-x(ixOmin1,ixOmin2,ixOmin3,3))
      do ix1=ixOmin1,ixOmax1,+1
        do ix2=ixOmin2+1,ixOmax2-1
        do ix3=ixOmin3+1,ixOmax3-1
         w(ix1,ix2,ix3,b1_)=w(ix1-2,ix2,ix3,b1_)-delxdely*(w(ix1-1,ix2+1,ix3,b2_)-w(ix1-1,ix2-1,ix3,b2_)) &
                                                -delxdelz*(w(ix1-1,ix2,ix3+1,b3_)-w(ix1-1,ix2,ix3-1,b3_))
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
case(3)
      ixIMmin2=ixOmax2+1;ixIMmax2=ixOmax2+1;
      ixIMmin1=ixOmin1;ixIMmax1=ixOmax1;
      {^IFTHREED ixIMmin3=ixOmin3;ixIMmax3=ixOmax3;}
      patchw(ixIM^S)=.false.
      ! primitiven (instead of primitive) to avoid changing internal values
      call primitiven(ixG^L,ixIM^L,w,patchw)
      do ix2=ixOmax2,ixOmin2,-1
        {^IFTWOD
{#IFDEF TRACER
        w(ixOmin1:ixOmax1,ix2,tr1_)= w(ixOmin1:ixOmax1,ixOmax2+1,tr1_)
}
        w(ixOmin1:ixOmax1,ix2,rho_)= w(ixOmin1:ixOmax1,ixOmax2+1,rho_)
        w(ixOmin1:ixOmax1,ix2,p_)  = w(ixOmin1:ixOmax1,ixOmax2+1,p_)
        w(ixOmin1:ixOmax1,ix2,v1_) = w(ixOmin1:ixOmax1,ixOmax2+1,v1_)
        w(ixOmin1:ixOmax1,ix2,v2_) = w(ixOmin1:ixOmax1,ixOmax2+1,v2_)
        w(ixOmin1:ixOmax1,ix2,v3_) = w(ixOmin1:ixOmax1,ixOmax2+1,v3_)
        !w(ixOmin1:ixOmax1,ix2,b1_) = w(ixOmin1:ixOmax1,ixOmax2+1,b1_)
        !w(ixOmin1:ixOmax1,ix2,b2_) = w(ixOmin1:ixOmax1,ixOmax2+1,b2_)
        w(ixOmin1:ixOmax1,ix2,b1_) = two*x(ixOmin1:ixOmax1,ix2,1)*x(ixOmin1:ixOmax1,ix2,2) &
                                    /(x(ixOmin1:ixOmax1,ix2,1)**2+x(ixOmin1:ixOmax1,ix2,2)**2)**2
        w(ixOmin1:ixOmax1,ix2,b2_) = -one+ (x(ixOmin1:ixOmax1,ix2,2)**2-x(ixOmin1:ixOmax1,ix2,1)**2) &
                                          /(x(ixOmin1:ixOmax1,ix2,1)**2+x(ixOmin1:ixOmax1,ix2,2)**2)**2
        w(ixOmin1:ixOmax1,ix2,b3_) = eqpar(b3val_)
        }
        {^IFTHREED
{#IFDEF TRACER
        w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,tr1_)= w(ixOmin1:ixOmax1,ixOmax2+1,ixOmin3:ixOmax3,tr1_)
}
        w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,rho_)= w(ixOmin1:ixOmax1,ixOmax2+1,ixOmin3:ixOmax3,rho_)
        w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,p_)  = w(ixOmin1:ixOmax1,ixOmax2+1,ixOmin3:ixOmax3,p_)
        w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,v1_) = w(ixOmin1:ixOmax1,ixOmax2+1,ixOmin3:ixOmax3,v1_)
        w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,v2_) = w(ixOmin1:ixOmax1,ixOmax2+1,ixOmin3:ixOmax3,v2_)
        w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,v3_) = w(ixOmin1:ixOmax1,ixOmax2+1,ixOmin3:ixOmax3,v3_)
        !w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,b1_) = w(ixOmin1:ixOmax1,ixOmax2+1,ixOmin3:ixOmax3,b1_)
        !w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,b2_) = w(ixOmin1:ixOmax1,ixOmax2+1,ixOmin3:ixOmax3,b2_)
        w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,b1_) = two*x(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,1)*x(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,2) &
                                    /(x(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,1)**2+x(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,2)**2)**2
        w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,b2_) = -one+ (x(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,2)**2-x(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,1)**2) &
                                          /(x(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,1)**2+x(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,2)**2)**2
        w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,b3_) = eqpar(b3val_)
        }
      enddo
      {^IFTWOD
      delydelx=(x(ixOmin1,ixOmin2+1,2)-x(ixOmin1,ixOmin2,2))/(x(ixOmin1+1,ixOmin2,1)-x(ixOmin1,ixOmin2,1))
      do ix2=ixOmax2,ixOmin2,-1
        do ix1=ixOmin1+1,ixOmax1-1
         w(ix1,ix2,b2_)=w(ix1,ix2+2,b2_)+delydelx*(w(ix1+1,ix2+1,b1_)-w(ix1-1,ix2+1,b1_)) 
        enddo
      enddo
      }
      {^IFTHREED
      delydelx=(x(ixOmin1,ixOmin2+1,ixOmin3,2)-x(ixOmin1,ixOmin2,ixOmin3,2))/(x(ixOmin1+1,ixOmin2,ixOmin3,1)-x(ixOmin1,ixOmin2,ixOmin3,1))
      delydelz=(x(ixOmin1,ixOmin2+1,ixOmin3,2)-x(ixOmin1,ixOmin2,ixOmin3,2))/(x(ixOmin1,ixOmin2,ixOmin3+1,3)-x(ixOmin1,ixOmin2,ixOmin3,3))
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
      ! now switch to conservative in full bottom layer
      patchw(ixO^S)=.false.
      call conserve(ixG^L,ixO^L,w,x,patchw)
case(4)
      ixIMmin2=ixOmin2-1;ixIMmax2=ixOmin2-1;
      ixIMmin1=ixOmin1;ixIMmax1=ixOmax1;
      {^IFTHREED ixIMmin3=ixOmin3;ixIMmax3=ixOmax3;}
      patchw(ixIM^S)=.false.
      ! primitiven (instead of primitive) to avoid changing internal values
      call primitiven(ixG^L,ixIM^L,w,patchw)
      do ix2=ixOmin2,ixOmax2,+1
        {^IFTWOD 
{#IFDEF TRACER
        w(ixOmin1:ixOmax1,ix2,tr1_)= w(ixOmin1:ixOmax1,ixOmin2-1,tr1_)
}
        w(ixOmin1:ixOmax1,ix2,rho_)= w(ixOmin1:ixOmax1,ixOmin2-1,rho_)
        w(ixOmin1:ixOmax1,ix2,p_)  = w(ixOmin1:ixOmax1,ixOmin2-1,p_)
        w(ixOmin1:ixOmax1,ix2,v1_) = w(ixOmin1:ixOmax1,ixOmin2-1,v1_)
        w(ixOmin1:ixOmax1,ix2,v2_) = w(ixOmin1:ixOmax1,ixOmin2-1,v2_)
        w(ixOmin1:ixOmax1,ix2,v3_) = w(ixOmin1:ixOmax1,ixOmin2-1,v3_)
        !w(ixOmin1:ixOmax1,ix2,b1_) = w(ixOmin1:ixOmax1,ixOmin2-1,b1_)
        !w(ixOmin1:ixOmax1,ix2,b2_) = w(ixOmin1:ixOmax1,ixOmin2-1,b2_)
        w(ixOmin1:ixOmax1,ix2,b1_) = two*x(ixOmin1:ixOmax1,ix2,1)*x(ixOmin1:ixOmax1,ix2,2) &
                                    /(x(ixOmin1:ixOmax1,ix2,1)**2+x(ixOmin1:ixOmax1,ix2,2)**2)**2
        w(ixOmin1:ixOmax1,ix2,b2_) = -one+ (x(ixOmin1:ixOmax1,ix2,2)**2-x(ixOmin1:ixOmax1,ix2,1)**2) &
                                          /(x(ixOmin1:ixOmax1,ix2,1)**2+x(ixOmin1:ixOmax1,ix2,2)**2)**2
        w(ixOmin1:ixOmax1,ix2,b3_) = eqpar(b3val_)
        }
        {^IFTHREED 
{#IFDEF TRACER
        w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,tr1_)= w(ixOmin1:ixOmax1,ixOmin2-1,ixOmin3:ixOmax3,tr1_)
}
        w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,rho_)= w(ixOmin1:ixOmax1,ixOmin2-1,ixOmin3:ixOmax3,rho_)
        w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,p_)  = w(ixOmin1:ixOmax1,ixOmin2-1,ixOmin3:ixOmax3,p_)
        w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,v1_) = w(ixOmin1:ixOmax1,ixOmin2-1,ixOmin3:ixOmax3,v1_)
        w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,v2_) = w(ixOmin1:ixOmax1,ixOmin2-1,ixOmin3:ixOmax3,v2_)
        w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,v3_) = w(ixOmin1:ixOmax1,ixOmin2-1,ixOmin3:ixOmax3,v3_)
        !w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,b1_) = w(ixOmin1:ixOmax1,ixOmin2-1,ixOmin3:ixOmax3,b1_)
        !w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,b2_) = w(ixOmin1:ixOmax1,ixOmin2-1,ixOmin3:ixOmax3,b2_)
        w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,b1_) = two*x(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,1)*x(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,2) &
                                    /(x(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,1)**2+x(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,2)**2)**2
        w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,b2_) = -one+ (x(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,2)**2-x(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,1)**2) &
                                          /(x(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,1)**2+x(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,2)**2)**2
        w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,b3_) = eqpar(b3val_)
        }
      enddo
      {^IFTWOD
      delydelx=(x(ixOmin1,ixOmin2+1,2)-x(ixOmin1,ixOmin2,2))/(x(ixOmin1+1,ixOmin2,1)-x(ixOmin1,ixOmin2,1))
      do ix2=ixOmin2,ixOmax2,+1
        do ix1=ixOmin1+1,ixOmax1-1
         w(ix1,ix2,b2_)=w(ix1,ix2-2,b2_)-delydelx*(w(ix1+1,ix2-1,b1_)-w(ix1-1,ix2-1,b1_)) 
        enddo
      enddo
      }
      {^IFTHREED
      delydelx=(x(ixOmin1,ixOmin2+1,ixOmin3,2)-x(ixOmin1,ixOmin2,ixOmin3,2))/(x(ixOmin1+1,ixOmin2,ixOmin3,1)-x(ixOmin1,ixOmin2,ixOmin3,1))
      delydelz=(x(ixOmin1,ixOmin2+1,ixOmin3,2)-x(ixOmin1,ixOmin2,ixOmin3,2))/(x(ixOmin1,ixOmin2,ixOmin3+1,3)-x(ixOmin1,ixOmin2,ixOmin3,3))
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
subroutine specialsource(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

! Calculate w(iw)=w(iw)+qdt*SOURCE[wCT,qtC,x] within ixO for all indices
! iw=iwmin...iwmax.  wCT is at time qCT

use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L, iw^LIM
double precision, intent(in) :: qdt, qtC, qt, x(ixI^S,1:ndim)
double precision, intent(inout) :: wCT(ixI^S,1:nw), w(ixI^S,1:nw)

! integer :: iw
! double precision :: s(ixG^T)
!-----------------------------------------------------------------------------

! do iw= iw^LIM
!    select case(iw)
!    case(m1_)
!       ! The source is based on the time centered wCT
!       call getmyforce(wCT,ixO^L,s)
!       w(ixO^S,m1_)=w(ixO^S,m1_) + qdt*s(ixO^S)
!    case(e_)
!       call getmyheating(wCT,ixO^L,s)
!       w(ixO^S,e_) =w(ixO^S,e_)  + qdt*s(ixO^S)
!    end select
! end do

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

dtnew=bigdouble

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

! you must set consistent values for integers refine/coarsen:

! refine = -1 enforce to not refine
! refine =  0 doesn't enforce anything
! refine =  1 enforce refinement

! coarsen = -1 enforce to not coarsen
! coarsen =  0 doesn't enforce anything
! coarsen =  1 enforce coarsen

use mod_global_parameters

integer, intent(in) :: igrid, level, ixG^L, ix^L
double precision, intent(in) :: qt, w(ixG^S,1:nw), x(ixG^S,1:ndim)
integer, intent(inout) :: refine, coarsen
double precision:: rval(ixG^T)
!-----------------------------------------------------------------------------

rval(ix^S)=dsqrt(x(ix^S,1)**2+x(ix^S,2)**2)
if (all(rval(ix^S) < one)) refine=1

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
subroutine specialset_B0(ixI^L,ixO^L,x,wB0)

! Here one can add a steady (time-independent) potential background field

use mod_global_parameters

integer, intent(in)           :: ixI^L,ixO^L
double precision, intent(in)  :: x(ixG^T,1:ndim)
double precision, intent(inout) :: wB0(ixI^S,1:ndir)
!-----------------------------------------------------------------------------
call mpistop(' abs(Busr)> 0, make change in parfile or in user file')

wB0(ixO^S,1:ndir)=wB0(ixO^S,1:ndir)
!!wB0(ixO^S,1:ndir)=wB0(ixO^S,1:ndir)+user defined steady potential field

end subroutine specialset_B0

!=============================================================================
subroutine initglobaldata_usr

use constants
use mod_global_parameters
double precision           :: rho0, p0
!-----------------------------------------------------------------------------
eqpar(gamma_)=5.0d0/3.0d0
eqpar(eta_)=0.0001d0
select case (iprob)
  case(1)
eqpar(b3val_)=0.1d0
  case(2)
eqpar(b3val_)=0.0d0
  case(3)
eqpar(b3val_)=0.5d0
  case(4)
eqpar(b3val_)=1.0d0
  case(5)
eqpar(b3val_)=5.0d0
end select



rho0=one
p0=one/eqpar(gamma_)

UNIT_LENGTH   = 1.0d9 !10 megameter
UNIT_VELOCITY = sqrt((CONST_kB*1.0d6)/(CONST_mp)*(rho0/p0)) ! 10^6 K
UNIT_DENSITY  = CONST_mp * 1.0d9 ! 10^9 protons per cubic centimeter

normvar(0)    = UNIT_LENGTH
normvar(rho_) = UNIT_DENSITY
normvar(p_)   = UNIT_DENSITY*UNIT_VELOCITY**2
{^C&normvar(b^C_)   = sqrt(4.0d0*dpi*UNIT_DENSITY)*UNIT_VELOCITY \}
{^C&normvar(v^C_)   = UNIT_VELOCITY \}
normt = UNIT_LENGTH/UNIT_VELOCITY

if (mype .eq. 0) then 
   write(*,*) 'Scaling density, velocity and length:', UNIT_DENSITY, UNIT_VELOCITY, UNIT_LENGTH
   write(*,*) 'Unit time:', normt, ' seconds'
end if
end subroutine initglobaldata_usr
!=============================================================================
subroutine initonegrid_usr(ixG^L,ix^L,w,x)

! initialize one grid

use mod_global_parameters

integer, intent(in) :: ixG^L, ix^L
double precision, intent(in) :: x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)

double precision:: rho0,p0,epsilon,xJ1root,J0valatxJ1root,rlocal,cval
double precision:: rval(ixG^T),xtheta(ixG^T),J1vals(ixG^T),DJ1vals(ixG^T)
double precision:: psi0(ixG^T),phi0(ixG^T),tmp(ixG^T)

double precision:: BJ0,DJ0,BJ1,DJ1,BY0,DY0,BY1,DY1,xkz
integer:: idims,ix^D

logical, save :: first=.true.
logical:: patchw(ixG^T)
!----------------------------------------------------------------------------

epsilon=1.0d-4
rho0=one
p0=one/eqpar(gamma_)

w(ix^S,rho_)=rho0

rval(ixG^T)=dsqrt(x(ixG^T,1)**2+x(ixG^T,2)**2)
xtheta(ixG^T)=datan2(x(ixG^T,2),x(ixG^T,1))

! this is the first root of Bessel J_1
xJ1root=3.831705970207512d0
! evaluate the zero bessel function at the first root of J_1
call JY01A(xJ1root,BJ0,DJ0,BJ1,DJ1,BY0,DY0,BY1,DY1)
J0valatxJ1root=BJ0
! compute the local J1 evaluation
{do ix^DB = ixG^LLIM^DB\}
  rlocal=rval(ix^D)*xJ1root
  call JY01A(rlocal,BJ0,DJ0,BJ1,DJ1,BY0,DY0,BY1,DY1)
  J1vals(ix^D)=BJ1
  DJ1vals(ix^S)=DJ1
{enddo^D&\}

cval=2.0d0/(xJ1root*J0valatxJ1root)
where(rval(ixG^T)<one) 
 psi0(ixG^T)=cval*J1vals(ixG^T)*dcos(xtheta(ixG^T))
elsewhere
 psi0(ixG^T)=(rval(ixG^T)-one/rval(ixG^T))*dcos(xtheta(ixG^T))
endwhere

{^IFTHREED
xkz=two*dpi/(xprobmax3-xprobmin3)
}
phi0(ixG^T)=epsilon*dexp(-rval(ixG^T)**2)

! compute dphi0/dy
idims=2
select case(typegrad)
    case("central")
     call gradient(phi0,ixG^L,ix^L,idims,tmp)
    case("limited")
     call gradientS(phi0,ixG^L,ix^L,idims,tmp)
end select
w(ix^S,v1_)=tmp(ix^S){^IFTHREED *dsin(xkz*x(ix^S,3)) }
! compute dphi0/dx
idims=1
select case(typegrad)
     case("central")
      call gradient(phi0,ixG^L,ix^L,idims,tmp)
     case("limited")
      call gradientS(phi0,ixG^L,ix^L,idims,tmp)
end select
w(ix^S,v2_)=-tmp(ix^S){^IFTHREED *dsin(xkz*x(ix^S,3)) }
w(ix^S,v3_)=zero

!pressure gradient equilibrium
where(rval(ix^S)<one) 
 w(ix^S,p_)=p0+half*(xJ1root**2)*(psi0(ix^S)**2)
elsewhere
 w(ix^S,p_)=p0
endwhere

!w(ix^S,p_)=p0
  
! compute dpsi0/dy
idims=2
select case(typegrad)
    case("central")
     call gradient(psi0,ixG^L,ix^L,idims,tmp)
    case("limited")
     call gradientS(psi0,ixG^L,ix^L,idims,tmp)
end select
w(ix^S,b1_)=tmp(ix^S)
! compute dpsi0/dx
idims=1
select case(typegrad)
     case("central")
      call gradient(psi0,ixG^L,ix^L,idims,tmp)
     case("limited")
      call gradientS(psi0,ixG^L,ix^L,idims,tmp)
end select
w(ix^S,b2_)=-tmp(ix^S)


!force free equilibrium
w(ix^S,b3_)=eqpar(b3val_)


!where(rval(ix^S)<one)
! w(ix^S,b3_)=eqpar(b3val_)+(xJ1root)*(psi0(ix^S))
!elsewhere
! w(ix^S,b3_)=eqpar(b3val_)
!endwhere

{#IFDEF TRACER
w(ix^S,tr1_)=zero
where(rval(ix^S)<one.and.x(ix^S,1)>zero)
 w(ix^S,tr1_)=one
endwhere
where(rval(ix^S)<one.and.x(ix^S,1)<zero)
 w(ix^S,tr1_)=-one
endwhere
}

patchw(ix^S)=.false.
call conserve(ixG^L,ix^L,w,x,patchw)

if(mype==0.and.first)then
{#IFDEF D3
      write(*,*)'Doing 3D MHD, tilt problem'
}
{#IFDEF D2
      write(*,*)'Doing 2D MHD, tilt problem'
}
      first=.false.
endif


return
end subroutine initonegrid_usr
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
! .. local ..
!!double precision:: divb(ixG^T)
!!double precision :: current(ixG^T,7-2*ndir:3)
integer          :: idirmin,idir
double precision :: wloc(ixG^T,1:nw),tmp(ixG^T),bvec(ixG^T,1:ndir)
double precision :: pth(ixG^T), rho(ixG^T)
!-----------------------------------------------------------------------------

wloc(ixI^S,1:nw)=w(ixI^S,1:nw)
! current (split case as well, current from B1 only))
!!^C&bvec(ixI^S,^C)=wloc(ixI^S,b^C_);
!!call curlvector(bvec,ixI^L,ixO^L,current,idirmin,1,ndir)
!w(ixO^S,nw+3)= current(ixO^S,3)
!call getcurrent(w,ixI^L,ixO^L,idirmin,current)
!w(ixO^S,nw+1)=current(ixO^S,3)
!! calculate R.B
!!tmp(ixO^S)=zero
!!do idir=1,ndir; 
!!   tmp(ixO^S)=tmp(ixO^S)+eqpar(eta_)*current(ixO^S,idir)*bvec(ixO^S,idir)
!!enddo;
!!w(ixO^S,nw+1)=tmp(ixO^S)
!!w(ixO^S,nw+2)=dsqrt(current(ixO^S,1)**2+current(ixO^S,2)**2+current(ixO^S,3)**2) 
if(saveprim)then
  {#IFDEF ENERGY
  pth(ixO^S)=wloc(ixO^S,p_)
  }
  {#IFDEF ISO
  pth(ixO^S)=eqpar(adiab_)*wloc(ixO^S,rho_)**eqpar(gamma_)
  }
else
  call getpthermal(wloc,x,ixI^L,ixO^L,pth)
endif
rho(ixO^S)=wloc(ixO^S,rho_)
w(ixO^S,nw+1)=pth(ixO^S)/rho(ixO^S)
normconv(nw+1) = UNIT_VELOCITY**2

end subroutine specialvar_output
!=============================================================================
subroutine specialvarnames_output

! newly added variables need to be concatenated with the varnames/primnames string

use mod_global_parameters
!-----------------------------------------------------------------------------

!!primnames= TRIM(primnames)//' '//'E.B jtot'
!!wnames=TRIM(wnames)//' '//'E.B jtot'
primnames= TRIM(primnames)//' '//'Te'
wnames=TRIM(wnames)//' '//'Te'


end subroutine specialvarnames_output
!=============================================================================
subroutine printlog_special

! printlog: calculates volume averaged mean values

use mod_global_parameters

logical :: fileopen
integer :: iigrid, igrid, level, nleafs_level(1:nlevelshi), iw, i
double precision :: wmean(1:nw), volume(1:nlevelshi), volprob, voltotal
double precision :: dvolume(ixG^T), volumeflat(1:nlevelshi)
double precision :: tmpw(ixG^T)
double precision :: current(ixG^T,7-2*ndir:3)
integer :: numlevels, imon, idirmin
integer, dimension(1:nlevelshi) :: isum_send, isum_recv
double precision, dimension(1:nw+3+nlevelshi) :: dsum_send, dsum_recv
double precision, dimension(1:2) :: wmeanmore
double precision :: wmaxcur,wmincur,wmaxvel,wmaxcur_mype,wmincur_mype,wmaxvel_mype
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

wmaxcur=zero
wmaxcur_mype=zero
wmaxvel=zero
wmaxvel_mype=zero
wmincur=zero
wmincur_mype=zero
invgminone=one/(eqpar(gamma_)-one)
wmeanmore(1:2)=zero

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
   call getcurrent(pw(igrid)%w,ixG^LL,ixM^LL,idirmin,current)
   ! maximal current
   wmaxcur_mype=max(wmaxcur_mype,maxval(current(ixM^T,3)))
   ! minimal current
   wmincur_mype=min(wmincur_mype,minval(current(ixM^T,3)))
   ! just use array for velocity
   tmpw(ixM^T)=dsqrt((pw(igrid)%w(ixM^T,m1_)/pw(igrid)%w(ixM^T,rho_))**2 &
                        +(pw(igrid)%w(ixM^T,m2_)/pw(igrid)%w(ixM^T,rho_))**2 &
                        +(pw(igrid)%w(ixM^T,m3_)/pw(igrid)%w(ixM^T,rho_))**2)
   wmaxvel_mype=max(wmaxvel_mype,maxval(tmpw(ixM^T)))
{#IFDEF TRACER
   wmean(Dtr1_)=wmean(Dtr1_)+sum(dvolume(ixM^T)*pw(igrid)%w(ixM^T,Dtr1_)/pw(igrid)%w(ixM^T,rho_))
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
   wmeanmore(1)=wmeanmore(1)+eqpar(eta_)*sum(dvolume(ixM^T)*current(ixM^T,3)**2)
   call getpthermal(pw(igrid)%w,px(igrid)%x,ixG^LL,ixM^LL,tmpw)
   wmeanmore(2)=wmeanmore(2)+invgminone*sum(dvolume(ixM^T)*tmpw(ixM^T))
end do
if (slab) volume(levmin:levmax)=volumeflat(levmin:levmax)

voltotal=sum(volume(levmin:levmax))

call MPI_REDUCE(wmaxcur_mype,wmaxcur,1,MPI_DOUBLE_PRECISION, &
                MPI_MAX,0,icomm,ierrmpi)
call MPI_REDUCE(wmaxvel_mype,wmaxvel,1,MPI_DOUBLE_PRECISION, &
                MPI_MAX,0,icomm,ierrmpi)

call MPI_REDUCE(wmincur_mype,wmincur,1,MPI_DOUBLE_PRECISION, &
                MPI_MIN,0,icomm,ierrmpi)

numlevels=levmax-levmin+1
dsum_send(1:nw)=wmean(1:nw)
dsum_send(nw+1)=voltotal
dsum_send(nw+2:nw+1+numlevels)=volumeflat(levmin:levmax)
dsum_send(nw+2+numlevels:nw+3+numlevels)=wmeanmore(1:2)
call MPI_REDUCE(dsum_send,dsum_recv,nw+3+numlevels,MPI_DOUBLE_PRECISION, &
                MPI_SUM,0,icomm,ierrmpi)
isum_send(1:numlevels)=nleafs_level(levmin:levmax)
call MPI_REDUCE(isum_send,isum_recv,numlevels,MPI_INTEGER, &
                MPI_SUM,0,icomm,ierrmpi)

if (mype==0) then

   wmean(1:nw)=dsum_recv(1:nw)
   wmeanmore(1:2)=dsum_recv(nw+2+numlevels:nw+3+numlevels)
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
      write(wnameslog(i+3:i+17),"(a14)") "d1 d2 d3 d4 d5"
      i=i+17
      do level=1,mxnest
          i=i+3
          write(wnameslog(i:i+1),"(a,i1)") "c",level
      end do
      do level=1,mxnest
          i=i+3
          write(wnameslog(i:i+1),"(a,i1)") "n",level
      end do
      if (time_accurate) then
         if(residmin>smalldouble) then
           write(line,'(a15,a79)')"it   t  dt res ",wnameslog
         else
           write(line,'(a15,a79)')"it   t   dt    ",wnameslog
         endif
      else
         if(residmin>smalldouble) then
           write(line,'(a7,a79)')"it res ",wnameslog
         else
           write(line,'(a7,a79)')"it     ",wnameslog
         endif
      end if

      call MPI_FILE_WRITE(log_fh,line,len_trim(line),MPI_CHARACTER, &
                          status,ierrmpi)
   end if
   !!call MPI_FILE_WRITE(log_fh,new_line('a'),1,MPI_CHARACTER,status,ierrmpi)
   call MPI_FILE_WRITE(log_fh,achar(10),1,MPI_CHARACTER,status,ierrmpi)

   if (time_accurate) then
      if(residmin>smalldouble) then
         write(line,'(i7,3(e13.5))')it,t,dt,residual
      else
         write(line,'(i7,2(e13.5))')it,t,dt
      endif
   else
      if(residmin>smalldouble) then
         write(line,'(i7,1(e13.5))')it,residual
      else
         write(line,'(i7)')it
      endif
   end if
   call MPI_FILE_WRITE(log_fh,line,len_trim(line), &
                       MPI_CHARACTER,status,ierrmpi)
   do iw=1,nw
      write(line,'(e13.5)')wmean(iw)
      call MPI_FILE_WRITE(log_fh,line,len_trim(line), &
                          MPI_CHARACTER,status,ierrmpi)
   end do
   do imon=1,2
      write(line,'(e13.5)')wmeanmore(imon)
      call MPI_FILE_WRITE(log_fh,line,len_trim(line), &
                          MPI_CHARACTER,status,ierrmpi)
   end do
   write(line,'(e13.5)')wmaxcur
   call MPI_FILE_WRITE(log_fh,line,len_trim(line), &
                          MPI_CHARACTER,status,ierrmpi)

   write(line,'(e13.5)')wmincur
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
! amrvacusr.t.tiltmhd
!=============================================================================
SUBROUTINE JY01A(X,BJ0,DJ0,BJ1,DJ1,BY0,DY0,BY1,DY1)

!       =======================================================
!       Purpose: Compute Bessel functions J0(x), J1(x), Y0(x),
!                Y1(x), and their derivatives
!       Input :  x   --- Argument of Jn(x) & Yn(x) ( x ò 0 )
!       Output:  BJ0 --- J0(x)
!                DJ0 --- J0'(x)
!                BJ1 --- J1(x)
!                DJ1 --- J1'(x)
!                BY0 --- Y0(x)
!                DY0 --- Y0'(x)
!                BY1 --- Y1(x)
!                DY1 --- Y1'(x)
!       =======================================================
!
IMPLICIT DOUBLE PRECISION (A-H,O-Z)

double precision:: X,BJ0,DJ0,BJ1,DJ1,BY0,DY0,BY1,DY1
double precision:: PI,RP2,X2,R,EC,CS0,W0,R0,CS1,W1,R1
double precision:: T1,P0,Q0,T2,P1,Q1,CU
double precision:: A(12),B(12),A1(12),B1(12)

integer:: K,K0

        PI=3.141592653589793D0
        RP2=0.63661977236758D0
        X2=X*X
        IF (X==0.0D0) THEN
           BJ0=1.0D0
           BJ1=0.0D0
           DJ0=0.0D0
           DJ1=0.5D0
           BY0=-1.0D+300
           BY1=-1.0D+300
           DY0=1.0D+300
           DY1=1.0D+300
           RETURN
        ENDIF
        IF (X<=12.0D0) THEN
           BJ0=1.0D0
           R=1.0D0
           DO K=1,30
              R=-0.25D0*R*X2/(K*K)
              BJ0=BJ0+R
              IF (DABS(R)<DABS(BJ0)*1.0D-15) EXIT
           ENDDO
           BJ1=1.0D0
           R=1.0D0
           DO K=1,30
              R=-0.25D0*R*X2/(K*(K+1.0D0))
              BJ1=BJ1+R
              IF (DABS(R)<DABS(BJ1)*1.0D-15) EXIT
           ENDDO
           BJ1=0.5D0*X*BJ1
           EC=DLOG(X/2.0D0)+0.5772156649015329D0
           CS0=0.0D0
           W0=0.0D0
           R0=1.0D0
           DO K=1,30
              W0=W0+1.0D0/K
              R0=-0.25D0*R0/(K*K)*X2
              R=R0*W0
              CS0=CS0+R
              IF (DABS(R)<DABS(CS0)*1.0D-15) EXIT
           ENDDO
           BY0=RP2*(EC*BJ0-CS0)
           CS1=1.0D0
           W1=0.0D0
           R1=1.0D0
           DO K=1,30
              W1=W1+1.0D0/K
              R1=-0.25D0*R1/(K*(K+1))*X2
              R=R1*(2.0D0*W1+1.0D0/(K+1.0D0))
              CS1=CS1+R
              IF (DABS(R)<DABS(CS1)*1.0D-15) EXIT
           ENDDO
           BY1=RP2*(EC*BJ1-1.0D0/X-0.25D0*X*CS1)
        ELSE
         A(1)=-.7031250000000000D-01
         A(2)=.1121520996093750D+00
         A(3)=-.5725014209747314D+00
         A(4)=.6074042001273483D+01
         A(5)=-.1100171402692467D+03
         A(6)=.3038090510922384D+04
         A(7)=-.1188384262567832D+06
         A(8)=.6252951493434797D+07
         A(9)=-.4259392165047669D+09
         A(10)=.3646840080706556D+11
         A(11)=-.3833534661393944D+13
         A(12)=.4854014686852901D+15
!DATA A/-.7031250000000000D-01,.1121520996093750D+00,-.5725014209747314D+00,.6074042001273483D+01, &
!       -.1100171402692467D+03,.3038090510922384D+04,-.1188384262567832D+06,.6252951493434797D+07, &
!       -.4259392165047669D+09,.3646840080706556D+11,-.3833534661393944D+13,.4854014686852901D+15/
         B(1)=.7324218750000000D-01
         B(2)=-.2271080017089844D+00
         B(3)=.1727727502584457D+01
         B(4)=-.2438052969955606D+02
         B(5)=.5513358961220206D+03
         B(6)=-.1825775547429318D+05
         B(7)=.8328593040162893D+06
         B(8)=-.5006958953198893D+08
         B(9)=.3836255180230433D+10
         B(10)=-.3649010818849833D+12
         B(11)=.4218971570284096D+14
         B(12)=-.5827244631566907D+16
!DATA B/ .7324218750000000D-01,-.2271080017089844D+00,.1727727502584457D+01,-.2438052969955606D+02, &
!        .5513358961220206D+03,-.1825775547429318D+05,.8328593040162893D+06,-.5006958953198893D+08, &
!        .3836255180230433D+10,-.3649010818849833D+12,.4218971570284096D+14,-.5827244631566907D+16/
         A1(1)=.1171875000000000D+00
         A1(2)=-.1441955566406250D+00
         A1(3)=.6765925884246826D+00
         A1(4)=-.6883914268109947D+01
         A1(5)=.1215978918765359D+03
         A1(6)=-.3302272294480852D+04
         A1(7)=.1276412726461746D+06
         A1(8)=-.6656367718817688D+07
         A1(9)=.4502786003050393D+09
         A1(10)=-.3833857520742790D+11
         A1(11)=.4011838599133198D+13
         A1(12)=-.5060568503314727D+15
!DATA A1/.1171875000000000D+00,-.1441955566406250D+00,.6765925884246826D+00,-.6883914268109947D+01, &
!        .1215978918765359D+03,-.3302272294480852D+04,.1276412726461746D+06,-.6656367718817688D+07, &
!        .4502786003050393D+09,-.3833857520742790D+11,.4011838599133198D+13,-.5060568503314727D+15/
         B1(1)=-.1025390625000000D+00
         B1(2)=.2775764465332031D+00
         B1(3)=-.1993531733751297D+01
         B1(4)=.2724882731126854D+02
         B1(5)=-.6038440767050702D+03
         B1(6)=.1971837591223663D+05
         B1(7)=-.8902978767070678D+06
         B1(8)=.5310411010968522D+08
         B1(9)=-.4043620325107754D+10
         B1(10)=.3827011346598605D+12
         B1(11)=-.4406481417852278D+14
         B1(12)=.6065091351222699D+16
!DATA B1/-.1025390625000000D+00,.2775764465332031D+00,-.1993531733751297D+01,.2724882731126854D+02, &
!        -.6038440767050702D+03,.1971837591223663D+05,-.8902978767070678D+06,.5310411010968522D+08, &
!        -.4043620325107754D+10,.3827011346598605D+12,-.4406481417852278D+14,.6065091351222699D+16/
           K0=12
           IF (X>=35.0) K0=10
           IF (X>=50.0) K0=8
           T1=X-0.25D0*PI
           P0=1.0D0
           Q0=-0.125D0/X
           DO K=1,K0
              P0=P0+A(K)*X**(-2*K)
              Q0=Q0+B(K)*X**(-2*K-1)
           ENDDO
           CU=DSQRT(RP2/X)
           BJ0=CU*(P0*DCOS(T1)-Q0*DSIN(T1))
           BY0=CU*(P0*DSIN(T1)+Q0*DCOS(T1))
           T2=X-0.75D0*PI
           P1=1.0D0
           Q1=0.375D0/X
           DO K=1,K0
              P1=P1+A1(K)*X**(-2*K)
              Q1=Q1+B1(K)*X**(-2*K-1)
           ENDDO
           CU=DSQRT(RP2/X)
           BJ1=CU*(P1*DCOS(T2)-Q1*DSIN(T2))
           BY1=CU*(P1*DSIN(T2)+Q1*DCOS(T2))
        ENDIF
        DJ0=-BJ1
        DJ1=BJ0-BJ1/X
        DY0=-BY1
        DY1=BY0-BY1/X
        RETURN
        END
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

call spacial_integral_w

!if(mype==0)then
!   print *,'converting back to conservative variables'
!endif

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
trtreshold=0.001d0
patchwi(ixO^S)=(w(ixO^S,tr1_)>trtreshold)

return
end subroutine mask_gridtrpos
!=============================================================================
subroutine mask_gridtrneg(ixI^L,ixO^L,w,x,patchwi)

use mod_global_parameters

integer, intent(in)                :: ixI^L,ixO^L
double precision, intent(in)       :: x(ixI^S,1:ndim)
double precision, intent(in)       :: w(ixI^S,1:nw)
logical, intent(inout)             :: patchwi(ixG^T)

double precision :: trtreshold
!-----------------------------------------------------------------------------
! note we have primitives in w-array!
trtreshold=0.001d0
patchwi(ixO^S)=(w(ixO^S,tr1_)<-trtreshold)

return
end subroutine mask_gridtrneg
}
!=============================================================================
subroutine mask_gridjzpos(ixI^L,ixO^L,w,x,patchwi)

use mod_global_parameters

integer, intent(in)                :: ixI^L,ixO^L
double precision, intent(in)       :: x(ixI^S,1:ndim)
double precision, intent(in)       :: w(ixI^S,1:nw)
logical, intent(inout)             :: patchwi(ixG^T)

double precision :: bvec(ixG^T,1:ndir),current(ixG^T,7-2*ndir:3)
double precision :: currenttreshold
integer                            :: idirmin
!-----------------------------------------------------------------------------
! note we have primitives in w-array!
^C&bvec(ixI^S,^C)=w(ixI^S,b^C_);
call curlvector(bvec,ixI^L,ixO^L,current,idirmin,1,ndir)
currenttreshold=0.1d0
patchwi(ixO^S)=(current(ixO^S,3)>currenttreshold.and.&
               (dabs(x(ixO^S,1))<0.975*xprobmax1).and.&
               (dabs(x(ixO^S,2))<0.975*xprobmax2))

return
end subroutine mask_gridjzpos
!=============================================================================
subroutine mask_gridjzneg(ixI^L,ixO^L,w,x,patchwi)

use mod_global_parameters

integer, intent(in)                :: ixI^L,ixO^L
double precision, intent(in)       :: x(ixI^S,1:ndim)
double precision, intent(in)       :: w(ixI^S,1:nw)
logical, intent(inout)             :: patchwi(ixG^T)

double precision :: bvec(ixG^T,1:ndir),current(ixG^T,7-2*ndir:3)
double precision :: currenttreshold
integer                            :: idirmin
!-----------------------------------------------------------------------------
! note we have primitives in w-array!
^C&bvec(ixI^S,^C)=w(ixI^S,b^C_);
call curlvector(bvec,ixI^L,ixO^L,current,idirmin,1,ndir)
currenttreshold=0.1d0
patchwi(ixO^S)=(current(ixO^S,3)<-currenttreshold.and.&
               (dabs(x(ixO^S,1))<0.975*xprobmax1).and.&
               (dabs(x(ixO^S,2))<0.975*xprobmax2))

return
end subroutine mask_gridjzneg
!=============================================================================
subroutine spacial_integral_w

use mod_global_parameters

double precision :: dvolume(ixG^T), timephy,xmom,ymom
double precision, allocatable :: integral_ipe(:), integral_w(:)
double precision, external :: integral_grid

integer           :: nregions,ireg
integer           :: iigrid,igrid,status(MPI_STATUS_SIZE),ni
character(len=100):: filename,region
logical           :: patchwi(ixG^T),alive
!-----------------------------------------------------------------------------
nregions=3

{#IFDEF TRACER
nregions=5
}
do ireg=1,nregions
 select case(ireg)
 case(1)
   region='jzpos'
 case(2)
   region='jzneg'
 case(3)
   region='fulldomain'
{#IFDEF TRACER
 case(4)
   region='trpos'
 case(5)
   region='trneg'
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
  case('trneg')
     call mask_gridtrneg(ixG^LL,ixM^LL,pw(igrid)%w,px(igrid)%x,patchwi)
}
  case('jzpos')
     call mask_gridjzpos(ixG^LL,ixM^LL,pw(igrid)%w,px(igrid)%x,patchwi)
  case('jzneg')
     call mask_gridjzneg(ixG^LL,ixM^LL,pw(igrid)%w,px(igrid)%x,patchwi)
  case('fulldomain')
     patchwi(ixM^T)=.true.
  case default
     call mpistop("region not defined")
  end select
  ! volume of region
  integral_ipe(1)=integral_ipe(1)+ &
            integral_grid(ixG^LL,ixM^LL,pw(igrid)%w,px(igrid)%x,dvolume,1,patchwi)
  ! total mass of region
  integral_ipe(2)=integral_ipe(2)+ &
            integral_grid(ixG^LL,ixM^LL,pw(igrid)%w,px(igrid)%x,dvolume,2,patchwi)
  ! x-jzmoment of region
  integral_ipe(3)=integral_ipe(3)+ &
            integral_grid(ixG^LL,ixM^LL,pw(igrid)%w,px(igrid)%x,dvolume,3,patchwi)
  ! y-jzmoment of region
  integral_ipe(4)=integral_ipe(4)+ &
            integral_grid(ixG^LL,ixM^LL,pw(igrid)%w,px(igrid)%x,dvolume,4,patchwi)
  ! total T of region
  integral_ipe(5)=integral_ipe(5)+ &
            integral_grid(ixG^LL,ixM^LL,pw(igrid)%w,px(igrid)%x,dvolume,5,patchwi)
  ! total magnetic energy of region
  integral_ipe(6)=integral_ipe(6)+ &
            integral_grid(ixG^LL,ixM^LL,pw(igrid)%w,px(igrid)%x,dvolume,6,patchwi)
  ! total kinetic energy of region
  integral_ipe(7)=integral_ipe(7)+ &
            integral_grid(ixG^LL,ixM^LL,pw(igrid)%w,px(igrid)%x,dvolume,7,patchwi)
  ! total plasma beta of region
  integral_ipe(8)=integral_ipe(8)+ &
            integral_grid(ixG^LL,ixM^LL,pw(igrid)%w,px(igrid)%x,dvolume,8,patchwi)
  ! total current of region
  integral_ipe(9)=integral_ipe(9)+ &
            integral_grid(ixG^LL,ixM^LL,pw(igrid)%w,px(igrid)%x,dvolume,9,patchwi)
  ! total internal energy of region
  integral_ipe(10)=integral_ipe(10)+ &
            integral_grid(ixG^LL,ixM^LL,pw(igrid)%w,px(igrid)%x,dvolume,10,patchwi)
  ! total z-current of region
  integral_ipe(11)=integral_ipe(11)+ &
            integral_grid(ixG^LL,ixM^LL,pw(igrid)%w,px(igrid)%x,dvolume,11,patchwi)
  ! total Bz of region
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
  if(integral_w(1)<smalldouble) then
    xmom=0.0d0
    ymom=-1.0d0
  else
    xmom=integral_w(3)/integral_w(1)
    ymom=integral_w(4)/integral_w(1)
  endif
  write(21,'(15(es12.4))') timephy,integral_w(1:12),xmom,ymom
  close(21)
endif

deallocate(integral_ipe,integral_w)

enddo
end subroutine spacial_integral_w
!=============================================================================
function integral_grid(ixI^L,ixO^L,w,x,dvolume,intval,patchwi)

use mod_global_parameters

integer, intent(in)                :: ixI^L,ixO^L,intval
double precision, intent(in)       :: x(ixI^S,1:ndim),dvolume(ixG^T)
double precision, intent(in)       :: w(ixI^S,nw)
logical, intent(in) :: patchwi(ixG^T)

double precision, dimension(ixG^T,1:ndir) :: bvec,qvec,current
double precision :: integral_grid,tmp(ixG^T),tmp2(ixG^T),gammm1
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
  ! x-moment integration
  {do ix^DB=ixOmin^DB,ixOmax^DB\}
     if(patchwi(ix^D)) integral_grid=integral_grid+x(ix^D,1)*dvolume(ix^D)
  {end do\}
 case(4)
  ! y-moment integration
  {do ix^DB=ixOmin^DB,ixOmax^DB\}
     if(patchwi(ix^D)) integral_grid=integral_grid+x(ix^D,2)*dvolume(ix^D)
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
  ! eta J^2  integration (split case as well, current from B1 only))
  ^C&bvec(ixI^S,^C)=w(ixI^S,b^C_);
  call curlvector(bvec,ixI^L,ixO^L,current,idirmin,1,ndir)
  {do ix^DB=ixOmin^DB,ixOmax^DB\}
     if(patchwi(ix^D)) integral_grid=integral_grid+ &
                         eqpar(eta_)*(^C&current(ix^D,^C)**2+)*dvolume(ix^D)
  {end do\}
 case(10)
  ! internal energy integration
  {do ix^DB=ixOmin^DB,ixOmax^DB\}
     if(patchwi(ix^D)) integral_grid=integral_grid+(w(ix^D,p_)/gammm1)*dvolume(ix^D)
  {end do\}
 case(11)
  ! current z component
  ! current (split case as well, current from B1 only))
  ^C&bvec(ixI^S,^C)=w(ixI^S,b^C_);
  call curlvector(bvec,ixI^L,ixO^L,current,idirmin,1,ndir)
  {do ix^DB=ixOmin^DB,ixOmax^DB\}
     if(patchwi(ix^D)) integral_grid=integral_grid+current(ix^D,3)*dvolume(ix^D)
  {end do\}
 case(12)
  ! B z component
  if(B0field) then
    ^C&bvec(ixI^S,^C)=w(ixI^S,b^C_)+myB0_cell%w(ixI^S,^C);
  else
    ^C&bvec(ixI^S,^C)=w(ixI^S,b^C_);
  endif
  {do ix^DB=ixOmin^DB,ixOmax^DB\}
     if(patchwi(ix^D)) integral_grid=integral_grid+bvec(ix^D,3)*dvolume(ix^D)
  {end do\}
 case default
     call mpistop("intval not defined")
end select

return
end function integral_grid
!=============================================================================
{#IFDEF SPECIALTOLERANCE
subroutine special_tolerance(xlocal,tolerance)

!PURPOSE: use different tolerance in special regions for AMR to
!reduce/increase resolution there where nothing/something interesting happens.

use mod_global_parameters

double precision, intent(in) :: xlocal(1:ndim)
double precision, intent(inout) :: tolerance

double precision :: bczone^D,addtol,tol_add
!-----------------------------------------------------------------------------
!amplitude of additional tolerance
addtol=0.1d0
! thickness of near-boundary region
bczone1=0.2d0*(xprobmax1-xprobmin1)
! linear changing of additional tolerance
if(xlocal(1)-xprobmin1 < bczone1 .or. xprobmax1-xlocal(1) < bczone1) then
tol_add=(1.d0-min(xlocal(1)-xprobmin1,xprobmax1-xlocal(1))/bczone1)*addtol
endif
!! for tilt problem: ultimately reaches y boundaries
bczone2=0.1d0*(xprobmax2-xprobmin2)
if(xlocal(2)-xprobmin2 < bczone2 .or. xprobmax2-xlocal(2) < bczone2) then
tol_add=(1.d0-min(xlocal(2)-xprobmin2,xprobmax2-xlocal(2))/bczone2)*addtol
endif
!! for tilt problem: following not used, third direction is periodic
!!bczone3=0.2d0*(xprobmax3-xprobmin3)
!!if(xprobmax3-xlocal(3) < bczone3) then
!!  tol_add=(1.d0-(xprobmax3-xlocal(3))/bczone3)*addtol
!!endif
tolerance=tolerance+tol_add

end subroutine special_tolerance
!=============================================================================
}
