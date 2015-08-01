!=============================================================================
! amrvacusr.t.nul
!=============================================================================
!INCLUDE:amrvacnul/speciallog.t
INCLUDE:amrvacnul/specialbound.t
INCLUDE:amrvacnul/specialsource.t
INCLUDE:amrvacnul/specialimpl.t
INCLUDE:amrvacnul/usrflags.t
INCLUDE:amrvacnul/correctaux_usr.t
!=============================================================================
subroutine initglobaldata_usr

include 'amrvacdef.f'
!-----------------------------------------------------------------------------

eqpar(gamma_) = 4.0d0/3.0d0
eqpar(eta_)   = 0.01d0
eqpar(kappa_) = one/0.18d0

eqpar(beta0_) = 0.005d0
eqpar(rho0_)  = 0.005d0
eqpar(vcoll_) = 0.01d0

end subroutine initglobaldata_usr
!=============================================================================
subroutine initonegrid_usr(ixG^L,ix^L,w,x)

! initialize one grid within ix^L

include 'amrvacdef.f'

integer, intent(in) :: ixG^L, ix^L
double precision, intent(in) :: x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)
! .. local ..
double precision, dimension(ixG^T)  :: bphi1, bphi2, bz1, bz2, p1, p2, r1, r2, sinphi1, sinphi2, cosphi1, cosphi2
double precision                    :: pt0
integer                             :: ix^D
logical patchw(ixG^T)
double precision, parameter     :: ct = 0.186d0, alphat = 3.83d0, Ft = -1.10d0, p0 = -0.287087d0
!
patchw(ix^S) = .false.
!-----------------------------------------------------------------------------
! Do some geometry:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
r1(ix^S) = sqrt((x(ix^S,1)+one)**2+x(ix^S,2)**2)
r2(ix^S) = sqrt((x(ix^S,1)-one)**2+x(ix^S,2)**2)

sinphi1(ix^S) = x(ix^S,2)/r1(ix^S)
cosphi1(ix^S) = (x(ix^S,1)+one)/r1(ix^S)

sinphi2(ix^S) = x(ix^S,2)/r2(ix^S)
cosphi2(ix^S) = (x(ix^S,1)-one)/r2(ix^S)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Calculate pt0 from central plasma beta:
pt0 = eqpar(beta0_)/2.0d0 - p0

w(ix^S,1:nw) = zero
w(ix^S,rho_) = eqpar(rho0_)

w(ix^S,psi_)  = 0.0d0
w(ix^S,phib_) = 0.0d0
w(ix^S,q_)    = 0.0d0


! Toroidal magnetic field:
{^D& do ix^D=ixmin^D,ixmax^D\}
    bphi1(ix^D) = ct*alphat*bessel_j1(alphat*r1(ix^D))
    bphi2(ix^D) = ct*alphat*bessel_j1(alphat*r2(ix^D))
{^D& end do\}

! Poloidal (z) magnetic field:
{^D& do ix^D=ixmin^D,ixmax^D\}
    bz1(ix^D) = alphat * ( ct*bessel_j0(alphat*r1(ix^D)) - Ft/alphat**2 )
    bz2(ix^D) = alphat * ( ct*bessel_j0(alphat*r2(ix^D)) - Ft/alphat**2 )
{^D& end do\}

! Thermal pressure:
{^D& do ix^D=ixmin^D,ixmax^D\}
    p1(ix^D) = Ft * (ct * bessel_j0(alphat*r1(ix^D))-Ft/alphat**2) + pt0
    p2(ix^D) = Ft * (ct * bessel_j0(alphat*r2(ix^D))-Ft/alphat**2) + pt0
{^D& end do\}

! Fill and reset outside of r=1:
where(r1(ix^S) .lt. one)
   w(ix^S,pp_) = p1(ix^S)
   w(ix^S,b1_) = - sinphi1(ix^S) * bphi1(ix^S)
   w(ix^S,b2_) = + cosphi1(ix^S) * bphi1(ix^S)
   w(ix^S,b3_) = bz1(ix^S)
   w(ix^S,u1_) = + eqpar(vcoll_)
elsewhere(r2(ix^S) .lt. one)
   w(ix^S,pp_) = p2(ix^S)
   w(ix^S,b1_) = - sinphi2(ix^S) * bphi2(ix^S)
   w(ix^S,b2_) = + cosphi2(ix^S) * bphi2(ix^S)
   w(ix^S,b3_) = bz2(ix^S)
   w(ix^S,u1_) = - eqpar(vcoll_)
elsewhere
   w(ix^S,pp_) = pt0
   w(ix^S,b1_) = zero
   w(ix^S,b2_) = zero
   w(ix^S,b3_) = zero
end where

w(ix^S,lfac_)= sqrt(one + {^C& w(ix^S,u^C_)**2|+})
w(ix^S,e1_)  = w(ix^S,b2_)*w(ix^S,u3_) - w(ix^S,b3_)*w(ix^S,u2_)
w(ix^S,e2_)  = w(ix^S,b3_)*w(ix^S,u1_) - w(ix^S,b1_)*w(ix^S,u3_)
w(ix^S,e3_)  = w(ix^S,b1_)*w(ix^S,u2_) - w(ix^S,b2_)*w(ix^S,u1_)
{^C& w(ix^S,e^C_)  = w(ix^S,e^C_)/w(ix^S,lfac_)\}

call conserve(ixG^L,ix^L,w,x,patchw)
end subroutine initonegrid_usr
!=============================================================================
subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)

! this subroutine can be used in convert, to add auxiliary variables to the
! converted output file, for further analysis using tecplot, paraview, ....
! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
!
! the array normconv can be filled in the (nw+1:nw+nwauxio) range with 
! corresponding normalization values (default value 1)

include 'amrvacdef.f'

integer, intent(in)                :: ixI^L,ixO^L
double precision, intent(in)       :: x(ixI^S,1:ndim)
double precision                   :: w(ixI^S,nw+nwauxio)
double precision                   :: normconv(0:nw+nwauxio)

integer:: idirmin

double precision :: current(ixG^T,7-2*ndir:3)
!-----------------------------------------------------------------------------

call getcurrent(w,ixI^L,ixO^L,idirmin,current)

w(ixO^S,nw+1)=current(ixO^S,3)

end subroutine specialvar_output
!=============================================================================
subroutine specialvarnames_output

! newly added variables to be concatenated with the primnames/wnames string

include 'amrvacdef.f'
!-----------------------------------------------------------------------------
oktest = index(teststr,'printlog')>=1

primnames= TRIM(primnames)//' '//'jz'
wnames=TRIM(wnames)//' '//'jz'

end subroutine specialvarnames_output
!=============================================================================


















!=============================================================================
!================= JUST DUMMIES ==============================================
!=============================================================================
{#IFDEF FCT
subroutine initvecpot_usr(ixI^L, ixC^L, xC, A)

! initialize the vectorpotential on the corners
! used by b_from_vectorpotential()


include 'amrvacdef.f'

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

include 'amrvacdef.f'
!-----------------------------------------------------------------------------
oktest = index(teststr,'printlog')>=1

call mpistop("special log file undefined")

end subroutine printlog_special
!=============================================================================
subroutine process_grid_usr(igrid,level,ixI^L,ixO^L,qt,w,x)

! this subroutine is ONLY to be used for computing auxiliary variables
! which happen to be non-local (like div v), and are in no way used for
! flux computations. As auxiliaries, they are also not advanced

include 'amrvacdef.f'

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

include 'amrvacdef.f'

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

include 'amrvacdef.f'
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
include 'amrvacdef.f'

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
include 'amrvacdef.f'

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
