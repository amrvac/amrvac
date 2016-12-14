!=============================================================================
! amrvacusr.t.RM2D
INCLUDE:amrvacnul/speciallog.t
INCLUDE:amrvacnul/specialbound.t
INCLUDE:amrvacnul/specialsource.t
INCLUDE:amrvacnul/specialimpl.t
INCLUDE:amrvacnul/usrflags.t
!=============================================================================
subroutine initglobaldata_usr
use mod_hd_phys
use mod_global_parameters
!-----------------------------------------------------------------------------

end subroutine initglobaldata_usr
!=============================================================================
subroutine initonegrid_usr(ixG^L,ix^L,w,x)

  ! initialize one grid
  use mod_physics
  use mod_hd_phys
  use mod_global_parameters

integer, intent(in) :: ixG^L, ix^L
double precision, intent(in) :: x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)
logical,save :: first=.true.
!----------------------------------------------------------------------------
if(first)then
  if(mype==0) then
     write(*,*)'Simulate 2D RM problem'
  endif
  first=.false.
endif

{^IFONED call mpistop("This is a multi-D HD problem") }

! the left top  quadrant
where(dabs(x(ix^S,1))<0.5d0.and.x(ix^S,2)>0.5d0.and.x(ix^S,2)<1.d0)
   w(ix^S,e_)=0.3d0
   w(ix^S,rho_)=0.5323d0
   w(ix^S,mom(1))=1.206d0
   w(ix^S,mom(2))=zero
end where

! the right top  quadrant
where(dabs(x(ix^S,1))>0.5d0.and.x(ix^S,1)<1.d0.and.x(ix^S,2)>0.5d0.and.x(ix^S,2)<1.d0)
   w(ix^S,e_)=1.5d0
   w(ix^S,rho_)=1.5d0
   w(ix^S,mom(1))=zero
   w(ix^S,mom(2))=zero
end where

! the left bottom quadrant
where(dabs(x(ix^S,1))<0.5d0.and.x(ix^S,2)<0.5d0)
   w(ix^S,e_)=1.5d0
   w(ix^S,rho_)=1.5d0
   w(ix^S,mom(1))=zero
   w(ix^S,mom(2))=zero
end where

! the right bottom quadrant
where(dabs(x(ix^S,1))>0.5d0.and.x(ix^S,1)<1.d0.and.x(ix^S,2)<0.5d0)
   w(ix^S,e_)=0.3d0
   w(ix^S,rho_)=0.5323d0
   w(ix^S,mom(1))=zero
   w(ix^S,mom(2))=1.206d0
end where

call phys_to_conserved(ixG^L,ix^L,w,x)

end subroutine initonegrid_usr
!=============================================================================
! amrvacusr.t.RM2D
!=============================================================================
