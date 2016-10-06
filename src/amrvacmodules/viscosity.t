!=============================================================================
!
!    THE FOLLOWING SUBROUTINES ADD VISCOUS SOURCE TERMS AND CHECK DT
!
!------------------------------------------------------------------------------
!    Viscous forces in the momentum equations:
!
!    d m_i/dt +=  - div (eqpar(mu_) * PI)
!
!    !! Viscous work in the energy equation:
!
!    !! de/dt    += - div (v . eqpar(mu_) * PI)
!
!    where the PI stress tensor is
!
!    PI_i,j = - (dv_j/dx_i + dv_i/dx_j) + (2/3)*Sum_k dv_k/dx_k
!
!    where eqpar(mu_) is the dynamic viscosity coefficient (g cm^-1 s^-1). 
!    Positive value for mu defines a constant viscosity, and use 0 for no viscosity. 
!
!------------------------------------------------------------------------------
!    This can be added via amrvacusr.t as follows:
!
!    INCLUDE:amrvacmodules/viscosity.t
!    
!    In specialsource:
!    if(abs(eqpar(mu_))>smalldouble) call addsource_visc(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)
!    
!    In getdt_special:
!    if(abs(eqpar(mu_))>smalldouble) call getdt_visc(w,ixG^L,ix^L,dtnew,dx^D,x)
!
!    In initglobaldata_usr:
!    eqpar(mu_)    = x
!
!    And add mu_ as user-definable parameter in amrvacusrpar.t
!
!==============================================================================
subroutine addsource_visc(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

! Add viscosity source to w within ixO 

use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L, iw^LIM
double precision, intent(in) :: qdt, qtC, qt, x(ixI^S,1:ndim)
double precision, intent(inout) :: wCT(ixI^S,1:nw), w(ixI^S,1:nw)

integer:: ix^L,idims,idirs,jdirs,iw
double precision:: lambda(ixI^S,ndir,ndir),tmp(ixI^S),tmp2(ixI^S)
!-----------------------------------------------------------------------------
! standard case, textbook viscosity
if(eqpar(mu_)>zero)then

! Calculating viscosity sources 
{#IFNDEF FOURTHORDER
! involves second derivatives, two extra layers
ix^L=ixO^L^LADD2;
if({ ixImin^D>ixmin^D .or. ixImax^D<ixmax^D|.or.})&
  call mpistop("error for viscous source addition, 2 layers needed")
ix^L=ixO^L^LADD1;
}
{#IFDEF FOURTHORDER
! involves second derivatives, four extra layers
ix^L=ixO^L^LADD4;
if({ ixImin^D>ixmin^D .or. ixImax^D<ixmax^D|.or.})&
  call mpistop("error for viscous source addition, &
  requested fourth order gradients: 4 layers needed")
ix^L=ixO^L^LADD2;
}


! construct lambda tensor: lambda_ij = gradv_ij + gradv_ji
! initialize
lambda(ix^S,1:ndir,1:ndir)=zero

!next construct 
do idims=1,ndim; do idirs=1,ndir
! Calculate velocity gradient tensor within ixL: gradv= grad v, 
! thus gradv_ij=d_j v_i
  tmp(ixI^S)=wCT(ixI^S,m0_+idirs)/wCT(ixI^S,rho_)
  select case(typegrad)
  case("central")
    call gradient(tmp,ixI^L,ix^L,idims,tmp2)
  case("limited")
    call gradientS(tmp,ixI^L,ix^L,idims,tmp2)
  end select
  lambda(ix^S,idims,idirs)= lambda(ix^S,idims,idirs)+ tmp2(ix^S)
  lambda(ix^S,idirs,idims)= lambda(ix^S,idirs,idims)+ tmp2(ix^S)
enddo; enddo;

! Multiply lambda with viscosity coefficient and dt

do idirs=1,ndir; do jdirs=1,ndir
  lambda(ix^S,idirs,jdirs)=lambda(ix^S,idirs,jdirs)*eqpar(mu_)*qdt
enddo; enddo;

!calculate div v term through trace action separately
tmp(ix^S)=(^C&lambda(ix^S,^C,^C)+)/3

!substract trace from diagonal elements
do idirs=1,ndir
   lambda(ix^S,idirs,idirs)=lambda(ix^S,idirs,idirs)-tmp(ix^S)
enddo

! Calculate source terms from viscosity 
do iw= iw^LIM
   select case(iw)
   case(m^C_)
      ! dm/dt= +div(mu*[d_j v_i+d_i v_j]-(2*mu/3)* div v * kr) 
      ! hence m_j=m_j+d_i tensor_ji
      idirs=iw-m0_
      do idims=1,ndim 
            tmp(ix^S)=lambda(ix^S,idirs,idims)
            select case(typegrad)
            case("central")
              call gradient(tmp,ixI^L,ixO^L,idims,tmp2)
            case("limited")
              call gradientS(tmp,ixI^L,ixO^L,idims,tmp2)
            end select
            w(ixO^S,iw)=w(ixO^S,iw)+tmp2(ixO^S)
      enddo
{#IFDEF ENERGY
   case(e_)
      ! de/dt= +div(v.dot.[mu*[d_j v_i+d_i v_j]-(2*mu/3)* div v *kr])
      ! thus e=e+d_i v_j tensor_ji
      do idims=1,ndim
         tmp(ix^S)=(^C&wCT(ix^S,m^C_)*lambda(ix^S,^C,idims)+)/wCT(ix^S,rho_)
         call gradient(tmp,ixI^L,ixO^L,idims,tmp2)
         w(ixO^S,ee_)=w(ixO^S,ee_)+tmp2(ixO^S)
     enddo
}
   end select ! iw
end do        ! iiw

end if

end subroutine addsource_visc
!=============================================================================
subroutine getdt_visc(w,ixG^L,ix^L,dtnew,dx^D,x)

! Check diffusion time limit for dt < dtdiffpar * dx**2 / (mu/rho)

use mod_global_parameters

integer, intent(in) :: ixG^L, ix^L
double precision, intent(in) :: dx^D, x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw), dtnew

double precision :: tmp(ixG^S)
double precision:: dtdiff_visc, dxinv2(1:ndim)
integer:: idims

!-----------------------------------------------------------------------------
dtnew=bigdouble

if(abs(eqpar(mu_))<smalldouble) return

! Calculate the kinematic viscosity tmp=mu/rho

tmp(ix^S)=eqpar(mu_)/w(ix^S,rho_)

^D&dxinv2(^D)=one/dx^D**2;
do idims=1,ndim
   dtdiff_visc=dtdiffpar/maxval(tmp(ix^S)*dxinv2(idims))
   ! limit the time step
   dtnew=min(dtnew,dtdiff_visc)
enddo

end subroutine getdt_visc
!=============================================================================
