!=============================================================================
subroutine initglobaldata_usr

use mod_global_parameters
!-----------------------------------------------------------------------------

end subroutine initglobaldata_usr
!=============================================================================
subroutine initonegrid_usr(ixI^L,ixO^L,w,x)

! initialize one grid within ixO^L

use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L
double precision, intent(in) :: x(ixI^S,1:ndim)
double precision, intent(inout) :: w(ixI^S,1:nw)
logical patchw(ixG^T)
patchw(ixO^S) = .false.
!-----------------------------------------------------------------------------

w(ixO^S,1:nw)=zero

w(ixO^S,rho_)=1.0d0

{#IFDEF ENERGY
   w(ixO^S,pp_) = 1.0d0
}

call conserve(ixI^L,ixO^L,w,x,patchw)
end subroutine initonegrid_usr
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
