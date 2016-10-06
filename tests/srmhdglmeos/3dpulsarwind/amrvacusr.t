!=============================================================================
! amrvacusr.t.nul
!=============================================================================
!INCLUDE:amrvacnul.specialini.t
!INCLUDE:amrvacnul.speciallog.t
!INCLUDE:amrvacnul.specialbound.t
!INCLUDE:amrvacnul.specialsource.t
!=============================================================================
! amrvacusr.t.nul
!=============================================================================

subroutine initglobaldata_usr

use mod_global_parameters

{^IFGLM
eqpar(Cr_)=0.18d0
eqpar(Ch_)=0.5d0
}

eqpar(f0_)   = 1.0d0
eqpar(sig0_) = 1.0d1
eqpar(lor_)  = 5.0d0
eqpar(xi0_)  = 0.1d0
eqpar(ve_)   = 0.017d0
eqpar(rhoe_) = 1.0d4
eqpar(r0_)   = 0.5d0
eqpar(rc_)   = 2.0d0
eqpar(gamma_)= 4.0d0/3.0d0


end subroutine initglobaldata_usr

subroutine initonegrid_usr(ixG^L,ix^L,w,x)

use mod_global_parameters

integer, intent(in) :: ixG^L, ix^L
double precision, intent(in) :: x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)

! .. local ..
double precision  :: ftot(ixG^S), bphi(ixG^S,1:3), sinphi(ixG^S), cosphi(ixG^S), rcyl(ixG^S), &
     costheta(ixG^S), vR(ixG^S), rsphere(ixG^S), theta(ixG^S),AR(ixG^S,1:3)
logical :: patchw(ixG^T)
integer :: idirmin
!----------------------------------------------------------------------------
patchw(ixG^S)=.false.

{^IFGLM
  w(ixG^S,psi_) = zero
}

! Do some geometry: 

   rcyl(ixG^S)     = dsqrt(dabs(x(ixG^S,1))**2.0d0 + dabs(x(ixG^S,2))**2.0d0)
   rsphere(ixG^S)  = dsqrt({^C&dabs(x(ixG^S,^C))**2.0d0+})
   costheta(ixG^S) = x(ixG^S,3)/rsphere(ixG^S)
! Unique in [0,\pi): 
   theta(ixG^S)    = acos(costheta(ixG^S))  
   sinphi(ixG^S)   = x(ixG^S,2)/rcyl(ixG^S)
   cosphi(ixG^S)   = x(ixG^S,1)/rcyl(ixG^S)

! Calculate the vectorpotential:
   AR(ixG^S,1) = dsqrt(eqpar(f0_))*eqpar(xi0_)*x(ixG^S,1)/rsphere(ixG^S)**2.0d0 &
	* (x(ixG^S,3) + 2.0d0/dpi * (rcyl(ixG^S) - theta(ixG^S) * x(ixG^S,3)))
   AR(ixG^S,2) = dsqrt(eqpar(f0_))*eqpar(xi0_)*x(ixG^S,2)/rsphere(ixG^S)**2.0d0 &
	* (x(ixG^S,3) + 2.0d0/dpi * (rcyl(ixG^S) - theta(ixG^S) * x(ixG^S,3)))
   AR(ixG^S,3) = dsqrt(eqpar(f0_))*eqpar(xi0_)*x(ixG^S,3)/rsphere(ixG^S)**2.0d0 &
	* (x(ixG^S,3) + 2.0d0/dpi * (rcyl(ixG^S) - theta(ixG^S) * x(ixG^S,3)))


      ^C&bphi(ixG^S,^C)=zero;

   call curlvector(AR(ixG^S,1:3),ixG^L,ix^L,bphi(ixG^S,1:3),idirmin,1,3)


! Fill the external medium:  

  w(ix^S,rho_)=eqpar(rhoe_)
  w(ix^S,pp_) =1.0d0

   w(ix^S,v1_) =   eqpar(ve_) * x(ix^S,1)/rsphere(ix^S)
   w(ix^S,v2_) =   eqpar(ve_) * x(ix^S,2)/rsphere(ix^S)
   w(ix^S,v3_) =   eqpar(ve_) * x(ix^S,3)/rsphere(ix^S)

where (({^C&dabs(x(ix^S,^C))**2.0d0+}) < eqpar(rc_)**2.0d0)

! Omit the singularity at the origin:

   w(ix^S,lfac_) = (eqpar(lor_)-one)*(rsphere(ix^S)**4.0d0/(rsphere(ix^S)**4.0d0+(eqpar(r0_)/2.0d0)**4d0))+one
   vR(ix^S)   = dsqrt(one-one/w(ix^S,lfac_)**2.0d0)

!   bphi(ix^S) = dsqrt(eqpar(f0_))*eqpar(xi0_)/rsphere(ix^S)*&
!           sintheta(ix^S)*(one-2.0d0*theta(ix^S)/dpi)*x(ix^S,3)/dabs(x(ix^S,3))

   ftot(ix^S)  = eqpar(f0_)/rsphere(ix^S) &
      * (sin(theta(ix^S))**2.0d0+one/eqpar(sig0_))

   w(ix^S,rho_)= one
!(ftot(ix^S)-{^C& bphi(ix^S,^C)**2.0d0+}) &
!      / w(ix^S,lfac_)**2.0d0

   w(ix^S,pp_) = w(ix^S,rho_)/100.d0

! Fill vectors by rotating coordinate systems: 

   w(ix^S,b1_) =   bphi(ix^S,1)
   w(ix^S,b2_) =   bphi(ix^S,2)
   w(ix^S,b3_) =   bphi(ix^S,3)

   w(ix^S,v1_) =   vR(ix^S) * x(ix^S,1)/rsphere(ix^S)
   w(ix^S,v2_) =   vR(ix^S) * x(ix^S,2)/rsphere(ix^S)
   w(ix^S,v3_) =   vR(ix^S) * x(ix^S,3)/rsphere(ix^S)

end where


w(ix^S,lfac_)=one/dsqrt(one-({^C&w(ix^S,v^C_)**2.0d0+}))
if (useprimitiveRel) then
   {^C&w(ix^S,u^C_)=w(ix^S,lfac_)*w(ix^S,v^C_);\}
endif


  call conserve(ixG^L,ix^L,w,patchw)
end subroutine initonegrid_usr

!=============================================================================
subroutine bc_int(time,w,x,ixG^L,ix^L)

use mod_global_parameters

integer, intent(in) :: ixG^L, ix^L
double precision, intent(in) :: time
double precision, intent(inout) :: w(ixG^S,1:nw)
double precision, intent(in) :: x(ixG^S,1:ndim)

! .. local ..
double precision  :: ftot(ixG^S), bphi(ixG^S,1:3), sinphi(ixG^S), cosphi(ixG^S), rcyl(ixG^S), &
     costheta(ixG^S), vR(ixG^S), rsphere(ixG^S), theta(ixG^S), AR(ixG^S,1:3)
logical :: patchw(ixG^S)
integer :: idirmin
!----------------------------------------------------------------------------
patchw(ixG^S)=.true.

!unfortunately can't use where because of call to curlvector.
!where (({^C&dabs(x(ix^S,^C))**2.0d0+}) < (1.1d0*eqpar(r0_))**2.0d0)

if (minval({^C&dabs(x(ix^S,^C))**2.0d0+}) < (eqpar(r0_))**2.0d0) then

! Do some geometry: 

   rcyl(ixG^S)     = dsqrt(dabs(x(ixG^S,1))**2.0d0 + dabs(x(ixG^S,2))**2.0d0)
   rsphere(ixG^S)  = dsqrt({^C&dabs(x(ixG^S,^C))**2.0d0+})
   costheta(ixG^S) = x(ixG^S,3)/rsphere(ixG^S)
! Unique in [0,\pi): 
   theta(ixG^S)    = acos(costheta(ixG^S))  
   sinphi(ixG^S)   = x(ixG^S,2)/rcyl(ixG^S)
   cosphi(ixG^S)   = x(ixG^S,1)/rcyl(ixG^S)

! Calculate the vectorpotential:
   AR(ixG^S,1) = dsqrt(eqpar(f0_))*eqpar(xi0_)*x(ixG^S,1)/rsphere(ixG^S)**2.0d0 &
	* (x(ixG^S,3) + 2.0d0/dpi * (rcyl(ixG^S) - theta(ixG^S) * x(ixG^S,3)))
   AR(ixG^S,2) = dsqrt(eqpar(f0_))*eqpar(xi0_)*x(ixG^S,2)/rsphere(ixG^S)**2.0d0 &
	* (x(ixG^S,3) + 2.0d0/dpi * (rcyl(ixG^S) - theta(ixG^S) * x(ixG^S,3)))
   AR(ixG^S,3) = dsqrt(eqpar(f0_))*eqpar(xi0_)*x(ixG^S,3)/rsphere(ixG^S)**2.0d0 &
	* (x(ixG^S,3) + 2.0d0/dpi * (rcyl(ixG^S) - theta(ixG^S) * x(ixG^S,3)))
      ^C&bphi(ixG^S,^C)=zero;
   call curlvector(AR(ixG^S,1:3),ixG^L,ix^L,bphi(ixG^S,1:3),idirmin,1,3)


!end where

where (({^C&dabs(x(ix^S,^C))**2.0d0+}) < (eqpar(r0_))**2.0d0)
patchw(ix^S)=.false.  ! To conserve at the end of the routine


! Omit the singularity at the origin:

   w(ix^S,lfac_) = (eqpar(lor_)-one)*(rsphere(ix^S)**4.0d0/(rsphere(ix^S)**4.0d0+(eqpar(r0_)/2.0d0)**4d0))+one
   vR(ix^S)   = dsqrt(one-one/w(ix^S,lfac_)**2.0d0)

!   bphi(ix^S) = dsqrt(eqpar(f0_))*eqpar(xi0_)/rsphere(ix^S)*&
!           sintheta(ix^S)*(one-2.0d0*theta(ix^S)/dpi)*x(ix^S,3)/dabs(x(ix^S,3))
   ftot(ix^S)  = eqpar(f0_)/rsphere(ix^S) &
      * (sin(theta(ix^S))**2.0d0+one/eqpar(sig0_))

! Fill scalars: 

   w(ix^S,rho_)= one
!(ftot(ix^S)-{^C& bphi(ix^S,^C)**2.0d0+}) &
!      / w(ix^S,lfac_)**2.0d0

   w(ix^S,pp_) = w(ix^S,rho_)/100.d0

! Fill vectors by rotating coordinate systems: 

   w(ix^S,b1_) =   bphi(ix^S,1)
   w(ix^S,b2_) =   bphi(ix^S,2)
   w(ix^S,b3_) =   bphi(ix^S,3)

   w(ix^S,v1_) =   vR(ix^S) * x(ix^S,1)/rsphere(ix^S)
   w(ix^S,v2_) =   vR(ix^S) * x(ix^S,2)/rsphere(ix^S)
   w(ix^S,v3_) =   vR(ix^S) * x(ix^S,3)/rsphere(ix^S)

   w(ix^S,lfac_)=one/dsqrt(one-({^C&w(ix^S,v^C_)**2.0d0+}))
   {^C&w(ix^S,u^C_)=w(ix^S,lfac_)*w(ix^S,v^C_);\}

{^IFGLM
   w(ix^S,psi_) = zero
}

end where

	call conserve(ixG^L,ix^L,w,patchw)
end if


end subroutine bc_int
!=============================================================================
subroutine specialbound_usr(qt,ixG^L,ix^L,iw,iB,w,x)

! special boundary types, user defined

use mod_global_parameters

integer, intent(in) :: ixG^L, ix^L, iw, iB
double precision, intent(in) :: qt, x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)

! .. local ..
 logical          :: patchw(ixG^T)
 double precision :: ftot(ixG^T)
!----------------------------------------------------------------------------
patchw(ixG^S)=.false.

^C&w(ix^S,v^C_)=zero;
^C&w(ix^S,b^C_)=zero;


call conserve(ixG^L,ix^L,w,patchw)
end subroutine specialbound_usr
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
!--------------------------------------------------------------


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

if (iflag >nw)call mpistop(' iflag> nw, make change in &
      parfile or in user file')

var(ixI^S) = zero 

end subroutine specialvarforerrest
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
subroutine printlog_special

use mod_global_parameters
!-----------------------------------------------------------------------------
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
double precision:: divb(ixI^S)
!-----------------------------------------------------------------------------

call getdivb(w,ixI^L,ixO^L,divb)
 w(ixO^S,nw+1)=divb(ixO^S)

end subroutine specialvar_output
!=============================================================================
subroutine specialvarnames_output

! newly added variables need to be concatenated with the varnames/primnames string

use mod_global_parameters
!-----------------------------------------------------------------------------

! Example : as above in specialvar_output, assuming relativistic HD here...
 primnames= TRIM(primnames)//' '//'divb'
 wnames=TRIM(wnames)//' '//'divb'

end subroutine specialvarnames_output
!=============================================================================
subroutine Global_useroutput
use mod_global_parameters
end subroutine Global_useroutput
!=============================================================================
subroutine process_allgrid_usr(iit,qt)
use mod_global_parameters

! .. scalars ..
integer,intent(in)          :: iit
double precision, intent(in):: qt

!-----------------------------------------------------------------------------
end subroutine process_allgrid_usr
!=============================================================================
