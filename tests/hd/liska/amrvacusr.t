!=============================================================================
! amrvacusr.t.liska

INCLUDE:amrvacnul/speciallog.t
INCLUDE:amrvacnul/specialbound.t
INCLUDE:amrvacnul/specialsource.t
INCLUDE:amrvacnul/specialimpl.t
INCLUDE:amrvacnul/usrflags.t
INCLUDE:amrvacnul/correctaux_usr.t
!=============================================================================
subroutine initglobaldata_usr

use mod_global_parameters
{#IFDEF DUST
double precision:: qrhodust, qmin_ar, qmax_ar, r(0:^NDS)
integer:: i
}
!-----------------------------------------------------------------------------
eqpar(gamma_)=1.4d0
{#IFDEF ISO
eqpar(adiab_)=1.0d0
}
{#IFDEF DUST
qrhodust  = 2.9d0   ! density in dust in g/cm^3, typical for Carbon
                    ! silicon would be 3.3d0
qmax_ar=1.00d-4  ! maximal dust size in cm
qmin_ar=0.01d-4  ! minimal dust size in cm

rhodust(1:^NDS)= qrhodust
! here the dust sizes are defined.
! first dust sizes in Ndust bins, with all bins having equal total mass.
! To do this, assume the particle distribution goes as r^-3.5
r(0) = qmin_ar
do i=1,^NDS
    r(i) = (dsqrt(r(i-1)) +(dsqrt(qmax_ar)- &
        dsqrt(qmin_ar))/^NDS)**2.0d0
    dsdust(i) = r(i)-r(i-1)
end do
! now calculate the weigthed mean size of each bin, again assuming n goes as r^-3.5
do i=1,^NDS
    sdust(i) = (5.0d0/3.0d0)*(r(i)**(-1.5d0) - r(i-1)**(-1.5d0)) &
        /(r(i)**(-2.5d0) - r(i-1)**(-2.5d0))
end do
}

end subroutine initglobaldata_usr
!=============================================================================
subroutine initonegrid_usr(ixG^L,ix^L,w,x)

! initialize one grid 

use mod_global_parameters

integer, intent(in) :: ixG^L, ix^L
double precision, intent(in) :: x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)

double precision :: front, rho0, rho1, pr0, pr1
integer          :: iw, ix^D, kxO^L
logical          :: patchw(ixG^S)
!----------------------------------------------------------------------------

{^IFONED call mpistop("This is a multi-D HD problem") }
patchw=.false.

front=0.15d0
rho0=1.0d0
rho1=0.125d0
pr0=1.0d0
pr1=0.14d0
where(( ^D&x(ixG^S,^D)+ )>=front+smalldouble)
   w(ixG^S,rho_)=rho0
   w(ixG^S,m1_)=zero
   w(ixG^S,m2_)=zero
{#IFDEF ENERGY
   w(ixG^S,e_)=pr0/(eqpar(gamma_)-one)
}
{#IFDEF TRACER
{^FL& w(ixG^S,Dtr^FL_)=zero\}
}
elsewhere
   w(ixG^S,rho_)=rho1
   w(ixG^S,m1_)=zero
   w(ixG^S,m2_)=zero
{#IFDEF ENERGY
   w(ixG^S,e_)=pr1/(eqpar(gamma_)-one)
}
{#IFDEF TRACER
{^FL& w(ixG^S,Dtr^FL_)=w(ixG^S,rho_)\}
}
endwhere

{#IFDEF DUST
{^DS& w(ixG^S,rhod^DS_)=0.001d0*w(ixG^S,rho_)\}
{^DS& w(ixG^S,m1d^DS_)=zero \}
{^DS& w(ixG^S,m2d^DS_)=zero \}
}


! Smoothing of the initial condition:
{do ix^DB= ix^LIM^D\}
           {kxOmin^D= ix^D-1;
            kxOmax^D= ix^D+1;\}
            call primitiven(ixG^L,kxO^L,w,patchw)
            do iw = 1,nw
                   w(ix^D,iw)=sum(w(kxO^S,iw))&
                           /9.
            end do
            call conserven(ixG^L,kxO^L,w,patchw)
{enddo^D&\}

end subroutine initonegrid_usr
!=============================================================================
! amrvacusr.t.liska
!=============================================================================
