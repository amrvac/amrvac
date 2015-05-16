!=============================================================================
subroutine specialbound_usr(qt,ixI^L,ixO^L,iw,iB,w,x)

! special boundary types, user defined
! user must assign conservative variables in bounderies

include 'amrvacdef.f'

integer, intent(in) :: ixI^L, ixO^L, iw, iB
double precision, intent(in) :: qt, x(ixI^S,1:ndim)
double precision, intent(inout) :: w(ixI^S,1:nw)
!----------------------------------------------------------------------------

call mpistop("specialbound not defined")
! just to give an example for 3D MHD :
!select case(iB)
! case(1)
!   ! min boundary in the 1st dimension 
!   w(ixO^S,rho_)=1.d0
!   w(ixO^S,v1_)=0.d0
!   w(ixO^S,v2_)=0.d0
!   w(ixO^S,v3_)=0.d0
!   w(ixO^S,p_)=2.d0
!   w(ixO^S,b1_)=1.d0
!   w(ixO^S,b2_)=0.d0
!   w(ixO^S,b3_)=0.d0
!   call conserve(ixI^L,ixO^L,w,x,patchfalse)
! case(2)
!   ! max boundary in the 1st dimension
!   w(ixO^S,rho_)=1.d0
!   w(ixO^S,v1_)=0.d0
!   w(ixO^S,v2_)=0.d0
!   w(ixO^S,v3_)=0.d0
!   w(ixO^S,p_)=2.d0
!   w(ixO^S,b1_)=1.d0
!   w(ixO^S,b2_)=0.d0
!   w(ixO^S,b3_)=0.d0
!   call conserve(ixI^L,ixO^L,w,x,patchfalse)
! case(3)
!   ! min boundary in the 2nd dimension 
!   w(ixO^S,rho_)=1.d0
!   w(ixO^S,v1_)=0.d0
!   w(ixO^S,v2_)=0.d0
!   w(ixO^S,v3_)=0.d0
!   w(ixO^S,p_)=2.d0
!   w(ixO^S,b1_)=1.d0
!   w(ixO^S,b2_)=0.d0
!   w(ixO^S,b3_)=0.d0
!   call conserve(ixI^L,ixO^L,w,x,patchfalse)
! case(4)
!   ! max boundary in the 2nd dimension
!   w(ixO^S,rho_)=1.d0
!   w(ixO^S,v1_)=0.d0
!   w(ixO^S,v2_)=0.d0
!   w(ixO^S,v3_)=0.d0
!   w(ixO^S,p_)=2.d0
!   w(ixO^S,b1_)=1.d0
!   w(ixO^S,b2_)=0.d0
!   w(ixO^S,b3_)=0.d0
!   call conserve(ixI^L,ixO^L,w,x,patchfalse)
! case(5)
!   ! min boundary in the 3rd dimension 
!   w(ixO^S,rho_)=1.d0
!   w(ixO^S,v1_)=0.d0
!   w(ixO^S,v2_)=0.d0
!   w(ixO^S,v3_)=0.d0
!   w(ixO^S,p_)=2.d0
!   w(ixO^S,b1_)=1.d0
!   w(ixO^S,b2_)=0.d0
!   w(ixO^S,b3_)=0.d0
!   call conserve(ixI^L,ixO^L,w,x,patchfalse)
! case(6)
!   ! max boundary in the 3rd dimension
!   w(ixO^S,rho_)=1.d0
!   w(ixO^S,v1_)=0.d0
!   w(ixO^S,v2_)=0.d0
!   w(ixO^S,v3_)=0.d0
!   w(ixO^S,p_)=2.d0
!   w(ixO^S,b1_)=1.d0
!   w(ixO^S,b2_)=0.d0
!   w(ixO^S,b3_)=0.d0
!   call conserve(ixI^L,ixO^L,w,x,patchfalse)
!end select

end subroutine specialbound_usr
!=============================================================================
subroutine bc_int(level,qt,ixI^L,ixO^L,w,x)

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

integer, intent(in) :: ixI^L,ixO^L,level
double precision, intent(in) :: qt
double precision, intent(inout) :: w(ixI^S,1:nw)
double precision, intent(in) :: x(ixI^S,1:ndim)

! .. local ..
!logical :: patchw(ixG^T)
!----------------------------------------------------------------------------

call mpistop("bc_int not defined")

! just to give an example for relativistic MHD
!  -----------------------------------------
!patchw(ixO^S)=.true.
!where (({^D&x(ixO^S,^D)**2+})<half**2.0d0) 
!    patchw(ixO^S) = .false.
!  ^C&w(ixO^S,v^C_)=zero;
!  ^C&w(ixO^S,b^C_)=zero;
!    w(ixO^S,b3_) = one
!    w(ixO^S,v1_) = 0.99
!    w(ixO^S,rho_) = 1.d0
!    w(ixO^S,pp_)  = 2.0d0
!    w(ixO^S,lfac_)=one/dsqrt(one-({^C&w(ixO^S,v^C_)**2.0d0+}))
!end where
!!if (useprimitiveRel) then
!!  where (({^D&x(ixO^S,^D)**2+})<half**2.0d0) 
!!  {^C&w(ixO^S,u^C_)=w(ixO^S,lfac_)*w(ixO^S,v^C_);\}
!!  end where
!!endif
!call conserve(ixI^L,ixO^L,w,x,patchw)

end subroutine bc_int
!=============================================================================
