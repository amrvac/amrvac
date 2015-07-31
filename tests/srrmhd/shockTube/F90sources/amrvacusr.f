!=============================================================================
! amrvacusr.t.nul
!=============================================================================
!=============================================================================
subroutine printlog_special

include 'amrvacdef.f'
!-----------------------------------------------------------------------------
oktest = index(teststr,'printlog')>=1

call mpistop("special log file undefined")

end subroutine printlog_special
!=============================================================================
subroutine process_grid_usr(igrid,level,ixImin1,ixImax1,ixOmin1,ixOmax1,qt,w,&
   x)

! this subroutine is ONLY to be used for computing auxiliary variables
! which happen to be non-local (like div v), and are in no way used for
! flux computations. As auxiliaries, they are also not advanced

include 'amrvacdef.f'

integer, intent(in):: igrid,level,ixImin1,ixImax1,ixOmin1,ixOmax1
double precision, intent(in):: qt,x(ixImin1:ixImax1,1:ndim)
double precision, intent(inout):: w(ixImin1:ixImax1,1:nw)
!-----------------------------------------------------------------------------

end subroutine process_grid_usr

!=============================================================================
subroutine specialvar_output(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x,normconv)

! this subroutine can be used in convert, to add auxiliary variables to the
! converted output file, for further analysis using tecplot, paraview, ....
! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
!
! the array normconv can be filled in the (nw+1:nw+nwauxio) range with 
! corresponding normalization values (default value 1)

include 'amrvacdef.f'

integer, intent(in)                :: ixImin1,ixImax1,ixOmin1,ixOmax1
double precision, intent(in)       :: x(ixImin1:ixImax1,1:ndim)
double precision                   :: w(ixImin1:ixImax1,nw+nwauxio)
double precision                   :: normconv(0:nw+nwauxio)
!-----------------------------------------------------------------------------
oktest = index(teststr,'printlog')>=1

call mpistop("special output file undefined")

! Example: assuming nwauxio=1 at convert stage and desire to see -w(1)
! w(ixO^S,nw+1)=-w(ixO^S,1)

end subroutine specialvar_output
!=============================================================================
subroutine specialvarnames_output

! newly added variables to be concatenated with the primnames/wnames string

include 'amrvacdef.f'
!-----------------------------------------------------------------------------
oktest = index(teststr,'printlog')>=1

call mpistop("special wnames and primnames undefined")

! Example : as above in specialvar_output, assuming relativistic HD here...
! primnames= TRIM(primnames)//' '//'-rho'
! wnames=TRIM(wnames)//' '//'-d'

end subroutine specialvarnames_output
!=============================================================================



!=============================================================================
subroutine specialbound_usr(qt,ixImin1,ixImax1,ixOmin1,ixOmax1,iw,iB,w,x)

! special boundary types, user defined
! user must assign conservative variables in bounderies

include 'amrvacdef.f'

integer, intent(in) :: ixImin1,ixImax1, ixOmin1,ixOmax1, iw, iB
double precision, intent(in) :: qt, x(ixImin1:ixImax1,1:ndim)
double precision, intent(inout) :: w(ixImin1:ixImax1,1:nw)
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
subroutine bc_int(level,qt,ixImin1,ixImax1,ixOmin1,ixOmax1,w,x)

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

integer, intent(in) :: ixImin1,ixImax1,ixOmin1,ixOmax1,level
double precision, intent(in) :: qt
double precision, intent(inout) :: w(ixImin1:ixImax1,1:nw)
double precision, intent(in) :: x(ixImin1:ixImax1,1:ndim)

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
!=============================================================================
subroutine specialsource(qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,iwmin,iwmax,qtC,&
   wCT,qt,w,x)

! Calculate w(iw)=w(iw)+qdt*SOURCE[wCT,qtC,x] within ixO for all indices
! iw=iwmin...iwmax.  wCT is at time qCT

include 'amrvacdef.f'

integer, intent(in) :: ixImin1,ixImax1, ixOmin1,ixOmax1, iwmin,iwmax
double precision, intent(in) :: qdt, qtC, qt, x(ixImin1:ixImax1,1:ndim)
double precision, intent(inout) :: wCT(ixImin1:ixImax1,1:nw),&
    w(ixImin1:ixImax1,1:nw)

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
subroutine getdt_special(w,ixImin1,ixImax1,ixOmin1,ixOmax1,dtnew,dx1,x)

! Limit "dt" further if necessary, e.g. due to the special source terms.
! The getdt_courant (CFL condition) and the getdt subroutine in the AMRVACPHYS
! module have already been called.

include 'amrvacdef.f'

integer, intent(in) :: ixImin1,ixImax1, ixOmin1,ixOmax1
double precision, intent(in) :: dx1, x(ixImin1:ixImax1,1:ndim)
double precision, intent(inout) :: w(ixImin1:ixImax1,1:nw), dtnew
!-----------------------------------------------------------------------------

dtnew=bigdouble

end subroutine getdt_special
!=============================================================================
subroutine specialeta(w,ixImin1,ixImax1,ixOmin1,ixOmax1,idirmin,x,current,eta)

! Set the "eta" array for resistive MHD based on w or the
! "current" variable which has components between idirmin and 3.

include 'amrvacdef.f'

integer, intent(in) :: ixImin1,ixImax1, ixOmin1,ixOmax1, idirmin
double precision, intent(in) :: w(ixImin1:ixImax1,nw), x(ixImin1:ixImax1,&
   1:ndim)

double precision :: current(ixGlo1:ixGhi1,7-2*ndir:3), eta(ixGlo1:ixGhi1)
!-----------------------------------------------------------------------------

!  eta(ix^S)=...

call mpistop("specialeta is not defined")

end subroutine specialeta
!=============================================================================
subroutine specialrefine_grid(igrid,level,ixImin1,ixImax1,ixOmin1,ixOmax1,qt,&
   w,x,refine,coarsen)

! Enforce additional refinement or coarsening
! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.

! you must set consistent values for integers refine/coarsen:

! refine = -1 enforce to not refine
! refine =  0 doesn't enforce anything
! refine =  1 enforce refinement

! coarsen = -1 enforce to not coarsen
! coarsen =  0 doesn't enforce anything
! coarsen =  1 enforce coarsen

include 'amrvacdef.f'

integer, intent(in) :: igrid, level, ixImin1,ixImax1, ixOmin1,ixOmax1
double precision, intent(in) :: qt, w(ixImin1:ixImax1,1:nw),&
    x(ixImin1:ixImax1,1:ndim)
integer, intent(inout) :: refine, coarsen
!-----------------------------------------------------------------------------

! e.g. refine for negative first coordinate x < 0 as
!
! if (any(x(ix^S,1) < zero)) refine=1

end subroutine specialrefine_grid
!=============================================================================
subroutine specialvarforerrest(ixImin1,ixImax1,ixOmin1,ixOmax1,iflag,w,var)

! this is the place to compute a local auxiliary variable to be used
! as refinement criterion for the Lohner error estimator only
!  -->it is then requiring and iflag>nw
! note that ixO=ixI=ixG, hence the term local (gradients need special attention!)

include 'amrvacdef.f'

integer, intent(in)          :: ixImin1,ixImax1,ixOmin1,ixOmax1,iflag
double precision, intent(in) :: w(ixImin1:ixImax1,1:nw)
double precision, intent(out):: var(ixGlo1:ixGhi1)
!-----------------------------------------------------------------------------

if (iflag >nw)call mpistop&
   (' iflag> nw, make change in parfile or in user file')

var(ixImin1:ixImax1) = zero 

end subroutine specialvarforerrest
!=============================================================================
subroutine specialset_B0(ixImin1,ixImax1,ixOmin1,ixOmax1,x,wB0)

! Here one can add a steady (time-independent) potential background field

include 'amrvacdef.f'

integer, intent(in)           :: ixImin1,ixImax1,ixOmin1,ixOmax1
double precision, intent(in)  :: x(ixGlo1:ixGhi1,1:ndim)
double precision, intent(inout) :: wB0(ixImin1:ixImax1,1:ndir)
!-----------------------------------------------------------------------------
call mpistop(' abs(Busr)> 0, make change in parfile or in user file')

wB0(ixOmin1:ixOmax1,1:ndir)=wB0(ixOmin1:ixOmax1,1:ndir)
!!wB0(ixO^S,1:ndir)=wB0(ixO^S,1:ndir)+user defined steady potential field

end subroutine specialset_B0
!=============================================================================
!=============================================================================
subroutine specialsource_impl(qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,iwmin,iwmax,&
   qtC,wCT,qt,w,x)

include 'amrvacdef.f'

integer, intent(in) :: ixImin1,ixImax1, ixOmin1,ixOmax1, iwmin,iwmax
double precision, intent(in) :: qdt, qtC, qt, x(ixImin1:ixImax1,1:ndim)
double precision, intent(inout) :: w(ixImin1:ixImax1,1:nw),&
    wCT(ixImin1:ixImax1,1:nw)
!-----------------------------------------------------------------------------

end subroutine specialsource_impl
!=============================================================================
subroutine getdt_impl(w,ixGmin1,ixGmax1,ixmin1,ixmax1,dtnew,dx1,x)

include 'amrvacdef.f'

integer, intent(in) :: ixGmin1,ixGmax1,ixmin1,ixmax1
double precision, intent(in) :: dx1, x(ixGmin1:ixGmax1,1:ndim)
! note that depending on strictsmall etc, w values may change
double precision, intent(inout) :: w(ixGmin1:ixGmax1,1:nw), dtnew
!-----------------------------------------------------------------------------
dtnew=bigdouble

end subroutine getdt_impl
!=============================================================================
subroutine fixp_usr(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x)
include 'amrvacdef.f'

integer, intent(in)                :: ixImin1,ixImax1,ixOmin1,ixOmax1
double precision, intent(inout)    :: w(ixImin1:ixImax1,1:nw)
double precision, intent(in)       :: x(ixImin1:ixImax1,1:ndim)
!----------------------------------------------------------------------------


end subroutine fixp_usr
!=============================================================================
subroutine flag_grid_usr(qt,ixGmin1,ixGmax1,ixOmin1,ixOmax1,w,x,flag)

include 'amrvacdef.f'

integer, intent(in)             :: ixGmin1,ixGmax1, ixOmin1,ixOmax1
integer, intent(inout)          :: flag
double precision, intent(in)    :: qt
double precision, intent(inout) :: w(ixGmin1:ixGmax1,1:nw)
double precision, intent(in)    :: x(ixGmin1:ixGmax1,1:ndim)

! flag=-1 : Treat all cells active, omit deactivation (onentry, default)
! flag=0  : Treat as normal domain
! flag=1  : Treat as passive, but reduce by safety belt
! flag=2  : Always treat as passive

!-----------------------------------------------------------------------------
      
end subroutine flag_grid_usr
!=============================================================================
subroutine correctaux_usr(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x,patchierror,&
   subname)

include 'amrvacdef.f'

integer, intent(in)            :: ixImin1,ixImax1, ixOmin1,ixOmax1
integer, intent(inout)         :: patchierror(ixGlo1:ixGhi1)
character(len=*), intent(in)   :: subname
double precision, intent(inout):: w(ixImin1:ixImax1,1:nw)
double precision, intent(in)   :: x(ixImin1:ixImax1,1:ndim)

! correct solution from analytic case

end subroutine correctaux_usr
!==========================================================================================
!=============================================================================
subroutine initglobaldata_usr

include 'amrvacdef.f'
!-----------------------------------------------------------------------------

eqpar(gamma_) = 2.0d0
eqpar(eta_)   = 1.0d-3

end subroutine initglobaldata_usr
!=============================================================================
subroutine initonegrid_usr(ixGmin1,ixGmax1,ixmin1,ixmax1,w,x)

! initialize one grid within ix^L

include 'amrvacdef.f'

integer, intent(in) :: ixGmin1,ixGmax1, ixmin1,ixmax1
double precision, intent(in) :: x(ixGmin1:ixGmax1,1:ndim)
double precision, intent(inout) :: w(ixGmin1:ixGmax1,1:nw)
! .. loacl ..
!double precision,dimension(ixG^T)       :: ec
!-----------------------------------------------------------------------------

w(ixmin1:ixmax1,psi_)  = 0.0d0
w(ixmin1:ixmax1,phib_) = 0.0d0
w(ixmin1:ixmax1,q_)    = 0.0d0


where(x(ixmin1:ixmax1,1) .lt. 0.0d0) 
   w(ixmin1:ixmax1,rho_) = 1.0d0
   w(ixmin1:ixmax1,pp_)  = 10.0d0

   w(ixmin1:ixmax1,u1_)  = 4.925d0
   w(ixmin1:ixmax1,u2_)  = 0.0d0
   w(ixmin1:ixmax1,u3_)  = 0.0d0 

   w(ixmin1:ixmax1,b1_)  = 5.0d0
   w(ixmin1:ixmax1,b2_)  = 15.08d0
   w(ixmin1:ixmax1,b3_)  = 0.0d0
elsewhere
   w(ixmin1:ixmax1,rho_) = 7.930d0
   w(ixmin1:ixmax1,pp_)  = 274.1d0

   w(ixmin1:ixmax1,u1_)  = 0.6209d0
   w(ixmin1:ixmax1,u2_)  = 0.1009d0
   w(ixmin1:ixmax1,u3_)  = 0.0d0 

   w(ixmin1:ixmax1,b1_)  = 5.0d0
   w(ixmin1:ixmax1,b2_)  = 28.92d0
   w(ixmin1:ixmax1,b3_)  = 0.0d0
end where

call conserve(ixGmin1,ixGmax1,ixmin1,ixmax1,w,x,patchfalse)
 
call vcrossb(ixGmin1,ixGmax1,ixmin1,ixmax1,1,w,x,patchfalse,w(ixGlo1:ixGhi1,&
   e1_))
w(ixmin1:ixmax1,e1_) = - w(ixmin1:ixmax1,e1_)

  
call vcrossb(ixGmin1,ixGmax1,ixmin1,ixmax1,2,w,x,patchfalse,w(ixGlo1:ixGhi1,&
   e2_))
w(ixmin1:ixmax1,e2_) = - w(ixmin1:ixmax1,e2_)

  
call vcrossb(ixGmin1,ixGmax1,ixmin1,ixmax1,3,w,x,patchfalse,w(ixGlo1:ixGhi1,&
   e3_))
w(ixmin1:ixmax1,e3_) = - w(ixmin1:ixmax1,e3_)

end subroutine initonegrid_usr
!=============================================================================

!=============================================================================
! amrvacusr.t.nul
!=============================================================================
