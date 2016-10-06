!=============================================================================
! amrvacusr.t.nul
!=============================================================================
!INCLUDE:amrvacnul/speciallog.t
INCLUDE:amrvacnul/specialbound.t
INCLUDE:amrvacnul/specialsource.t
INCLUDE:amrvacnul/specialimpl.t
!INCLUDE:amrvacnul/usrflags.t
INCLUDE:amrvacnul/correctaux_usr.t
!=============================================================================
subroutine initglobaldata_usr

use mod_global_parameters
!-----------------------------------------------------------------------------

eqpar(kappa_)  = 100.0d0
eqpar(kpar_)   = 100.0d0
eqpar(kperp_)  = 0.0d0

eqpar(C_)     = 0.01d0
eqpar(vcoll_) = 0.1d0

end subroutine initglobaldata_usr
!=============================================================================
subroutine initonegrid_usr(ixG^L,ix^L,w,x)

! initialize one grid within ix^L

use mod_global_parameters

integer, intent(in) :: ixG^L, ix^L
double precision, intent(in) :: x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)
! .. local ..
double precision, dimension(ixG^T)  :: bphi1, bphi2, bz1, bz2, r1, r2, sinphi1, sinphi2, cosphi1, cosphi2
integer                             :: ix^D
logical patchw(ixG^T)
double precision, parameter     :: alphat = 3.8317d0
!
patchw(ixG^S) = .false.
!-----------------------------------------------------------------------------
! Do some geometry:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
r1(ixG^S) = sqrt((x(ixG^S,1)+one)**2+x(ixG^S,2)**2)
r2(ixG^S) = sqrt((x(ixG^S,1)-one)**2+x(ixG^S,2)**2)

sinphi1(ixG^S) = x(ixG^S,2)/r1(ixG^S)
cosphi1(ixG^S) = (x(ixG^S,1)+one)/r1(ixG^S)

sinphi2(ixG^S) = x(ixG^S,2)/r2(ixG^S)
cosphi2(ixG^S) = (x(ixG^S,1)-one)/r2(ixG^S)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

w(ixG^S,1:nw) = zero

w(ixG^S,phib_) = 0.0d0
w(ixG^S,psi_)  = 0.0d0

! Toroidal magnetic field:
{^D& do ix^D=ixGmin^D,ixGmax^D\}
    bphi1(ix^D) = bessel_j1(alphat*r1(ix^D))
    bphi2(ix^D) = bessel_j1(alphat*r2(ix^D))
{^D& end do\}

! Poloidal (z) magnetic field:
{^D& do ix^D=ixGmin^D,ixGmax^D\}
    bz1(ix^D) = sqrt(bessel_j0(alphat*r1(ix^D))**2 + eqpar(C_))
    bz2(ix^D) = sqrt(bessel_j0(alphat*r2(ix^D))**2 + eqpar(C_))
{^D& end do\}

! Fill and reset outside of r=1:
where(r1(ixG^S) .lt. one)
   w(ixG^S,b1_) = - sinphi1(ixG^S) * bphi1(ixG^S)
   w(ixG^S,b2_) = + cosphi1(ixG^S) * bphi1(ixG^S)
   w(ixG^S,b3_) = bz1(ixG^S)
   w(ixG^S,e1_)  =   zero
   w(ixG^S,e2_)  =   w(ixG^S,b3_)*eqpar(vcoll_)
   w(ixG^S,e3_)  = - w(ixG^S,b2_)*eqpar(vcoll_)
elsewhere(r2(ixG^S) .lt. one)
   w(ixG^S,b1_) = - sinphi2(ixG^S) * bphi2(ixG^S)
   w(ixG^S,b2_) = + cosphi2(ixG^S) * bphi2(ixG^S)
   w(ixG^S,b3_) = bz2(ixG^S)
   w(ixG^S,e1_)  =   zero
   w(ixG^S,e2_)  = -  w(ixG^S,b3_)*eqpar(vcoll_)
   w(ixG^S,e3_)  = + w(ixG^S,b2_)*eqpar(vcoll_)
elsewhere
   w(ixG^S,b1_) = zero
   w(ixG^S,b2_) = zero
   w(ixG^S,b3_) = sqrt(bessel_j0(alphat)**2 + eqpar(C_))
   w(ixG^S,e1_)  =   zero
   w(ixG^S,e2_)  =   zero
   w(ixG^S,e3_)  =   zero
end where

call divvector(w(ixG^T,e1_:e3_),ixG^LL,ix^L,w(ixG^T,q_))

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
double precision,dimension(ixG^T,1:ndir) :: current
double precision,dimension(ixG^T,1:ndir) :: curlb
double precision,dimension(ixG^T)        :: vidir, divE
double precision                         :: evec(ixG^T,1:ndir)
integer                                  :: idirmin
integer, parameter                       :: idirmin0=1
!-----------------------------------------------------------------------------

call getcurrent(ixI^L,ixO^L,w,x,saveprim,current)
w(ixO^S,nw+1) = current(ixO^S,3)

! Also get curl of b:
call curlvector(w(ixG^T,b1_:b3_),ixI^L,ixO^L,curlb,idirmin,idirmin0,ndir)
w(ixO^S,nw+2) = curlb(ixO^S,3)

! Get drift velocity:
{^C&
call getv(w,x,ixI^L,ixO^L,^C,vidir)
w(ixO^S,nw+2+^C) = vidir(ixO^S)
\}

! S[Ei_] = - divE (E x B)/B^2
evec(ixI^S,1:ndir)=w(ixI^S,e0_+1:e0_+ndir)
call divvector(evec,ixI^L,ixO^L,divE)
w(ixO^S,nw+2+ndir+1) = divE(ixO^S)

end subroutine specialvar_output
!=============================================================================
subroutine specialvarnames_output

! newly added variables to be concatenated with the primnames/wnames string

use mod_global_parameters
!-----------------------------------------------------------------------------

primnames= TRIM(primnames)//' '//'jz'//' '//'curlbz'//' '//'v1'//' '//'v2'//' '//'v3'//' '//'divE'
wnames=TRIM(wnames)//' '//'jz'//' '//'curlbz'//' '//'v1'//' '//'v2'//' '//'v3'//' '//'divE'

end subroutine specialvarnames_output
!=============================================================================

















!=============================================================================
!================= JUST DUMMIES ==============================================
!=============================================================================
subroutine flag_grid_usr(qt,ixG^L,ixO^L,w,x,flag)

use mod_global_parameters

integer, intent(in)             :: ixG^L, ixO^L
integer, intent(inout)          :: flag
double precision, intent(in)    :: qt
double precision, intent(inout) :: w(ixG^S,1:nw)
double precision, intent(in)    :: x(ixG^S,1:ndim)

! flag=-1 : Treat all cells active, omit deactivation (onentry, default)
! flag=0  : Treat as normal domain
! flag=1  : Treat as passive, but reduce by safety belt
! flag=2  : Always treat as passive

!-----------------------------------------------------------------------------
      
end subroutine flag_grid_usr
!=============================================================================
{#IFDEF FCT
subroutine initvecpot_usr(ixI^L, ixC^L, xC, A)

! initialize the vectorpotential on the corners
! used by b_from_vectorpotential()


use mod_global_parameters

integer, intent(in)                :: ixI^L, ixC^L
double precision, intent(in)       :: xC(ixI^S,1:ndim)
double precision, intent(out)      :: A(ixI^S,1:ndir)

!double precision                   :: r(ixG^T)

!-----------------------------------------------------------------------------

!r(ixC^S)=sqrt(xC(ixC^S,1)**2 + xC(ixC^S,2)**2 )

A(ixC^S,1:ndir) = zero
!where (r(ixC^S) .lt. eqpar(rm_))
!   A(ixC^S,3) = - half * eqpar(bm_) / eqpar(rm_) * r(ixC^S)**2
!elsewhere (r(ixC^S) .lt. eqpar(rj_))
!   A(ixC^S,3) = - half * eqpar(bm_) * eqpar(rm_) &
!        - eqpar(bm_) * eqpar(rm_) * log(r(ixC^S)/eqpar(rm_))
!elsewhere
!   A(ixC^S,3) = - half * eqpar(bm_) * eqpar(rm_) &
!        - eqpar(bm_) * eqpar(rm_) * log(eqpar(rj_)/eqpar(rm_))
!end where


end subroutine initvecpot_usr
!=============================================================================
}



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
{#IFDEF PROCESSGLOBAL
!==============================================================================

subroutine process_global_usr(iit,qt)
!
! This subroutine is called at the beginning of each time step 
! by each processor. No communication is specified, so the user
! has to implement MPI routines if information has to be shared
!

use mod_global_parameters

integer, intent(in)          :: iit
double precision, intent(in) :: qt

!-----------------------------------------------


return
end subroutine process_global_usr
}
!=============================================================================
{#IFDEF UCONVERT
subroutine userspecialconvert(qunitconvert)
! Allow user to use their own data-postprocessing procedures

use mod_global_parameters
integer, intent(in) :: qunitconvert
character(len=20):: userconvert_type
!-----------------------------------------------------------------------------

end subroutine userspecialconvert
!=============================================================================
}
{#IFDEF TRANSFORMW
subroutine transformw_usr(w,wtf,eqpar_tf,ixI^L,ixO^L)
! regenerate w and eqpar arrays to output into *tf.dat, e.g., add/remove e_
! variable
use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L
double precision, intent(in) :: w(ixI^S,1:nw)
double precision, intent(out):: wtf(ixI^S,1:nwtf)
double precision, intent(out):: eqpar_tf(neqpartf)

!-----------------------------------------------------------------------------

end subroutine transformw_usr
!=============================================================================
}
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
addtol=0.3d0
! thickness of near-boundary region
bczone1=0.2d0*(xprobmax1-xprobmin1)
! linear changing of additional tolerance
if(xlocal(1)-xprobmin1 < bczone1 .or. xprobmax1-xlocal(1) < bczone1) then
  tol_add=(1.d0-min(xlocal(1)-xprobmin1,xprobmax1-xlocal(1))/bczone1)*addtol
endif
bczone2=0.2d0*(xprobmax2-xprobmin2)
if(xlocal(2)-xprobmin2 < bczone2 .or. xprobmax2-xlocal(2) < bczone2) then
  tol_add=(1.d0-min(xlocal(2)-xprobmin2,xprobmax2-xlocal(2))/bczone2)*addtol
endif
bczone3=0.2d0*(xprobmax3-xprobmin3)
if(xprobmax3-xlocal(3) < bczone3) then
  tol_add=(1.d0-(xprobmax3-xlocal(3))/bczone3)*addtol
endif
tolerance=tolerance+tol_add

end subroutine special_tolerance
!=============================================================================
}

!=============================================================================
! amrvacusr.t.nul
!=============================================================================
