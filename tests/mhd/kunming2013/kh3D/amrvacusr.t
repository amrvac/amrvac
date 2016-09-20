!=============================================================================
! amrvacusr.t.khmhd

! INCLUDE:amrvacnul/specialini.t
! INCLUDE:amrvacnul/speciallog.t
INCLUDE:amrvacnul/specialbound.t
INCLUDE:amrvacnul/specialsource.t
INCLUDE:amrvacnul/specialimpl.t
INCLUDE:amrvacnul/usrflags.t
INCLUDE:amrvacnul/correctaux_usr.t
!=============================================================================
subroutine initglobaldata_usr

use mod_global_parameters
!-----------------------------------------------------------------------------
eqpar(gamma_)=5.0d0/3.0d0
eqpar(eta_)=zero

end subroutine initglobaldata_usr
!=============================================================================
subroutine initonegrid_usr(ixG^L,ix^L,w,x)

! initialize one grid

use mod_global_parameters

integer, intent(in) :: ixG^L, ix^L
double precision, intent(in) :: x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)

double precision:: mpol,k1,Rjet,width,qv,B0,dv
double precision, dimension(ixG^T) :: r,phi

logical, save :: first=.true.
logical:: patchw(ixG^T)
!----------------------------------------------------------------------------
oktest = index(teststr,'initonegrid_usr')>=1
if (oktest) write(unitterm,*) ' === initonegrid_usr  (in ) : ', &
      'ixG^L : ',ixG^L

! KH in 3D, Keppens & Toth
r(ixG^S)=dsqrt(x(ixG^S,2)**2+x(ixG^S,3)**2)
phi(ixG^S)=datan2(x(ixG^S,3),x(ixG^S,2))

mpol=one
k1=two*dpi/(xprobmax1-xprobmin1)
Rjet=0.5d0
width=0.1d0*Rjet
qv=0.645d0
B0=0.129d0
dv=0.01d0

w(ix^S,rho_)=one
w(ix^S,p_)=one
w(ix^S,b1_)=B0
w(ix^S,b2_)=zero
w(ix^S,b3_)=zero
w(ix^S,v1_)=qv*dtanh((r(ix^S)-Rjet)/width)
w(ix^S,v2_)=dv*dexp(-((r(ix^S)-Rjet)/(4.0d0*width))**2) &
               *dcos(mpol*phi(ix^S))*dsin(k1*x(ix^S,1))*dcos(phi(ix^S))
w(ix^S,v3_)=dv*dexp(-((r(ix^S)-Rjet)/(4.0d0*width))**2) &
               *dcos(mpol*phi(ix^S))*dsin(k1*x(ix^S,1))*dsin(phi(ix^S))

patchw(ix^S)=.false.
call conserve(ixG^L,ix^L,w,x,patchw)

if(first)then
      if(mype==0)then
       write(*,*)'Doing 3D ideal MHD, KH jet problem'
      endif
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
!-----------------------------------------------------------------------------

end subroutine specialvar_output
!=============================================================================
subroutine specialvarnames_output

! newly added variables need to be concatenated with the varnames/primnames string

use mod_global_parameters
!-----------------------------------------------------------------------------

end subroutine specialvarnames_output
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
! amrvacusr.t.khmhd
!=============================================================================
