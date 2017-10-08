!=============================================================================
! amrvacusr.t.movstar
! 
! update : 31/08/2015
! made by Allard Jan van Marle, small modifications by R. Keppens
! 2D axially symmetric Hydro setup, with optically thin radiative cooling
!    simulates the wind-blown bubble of a mass-loosing star, 
!     moving with respect to ISM
!
! For setup in 2D use:
!$AMRVAC_DIR/setup.pl -d=22 -phi=0 -z=2 -g=14,14 -p=hd -eos=default -nf=0 -ndust=0 -u=nul -arch=default



INCLUDE:amrvacnul/correctaux_usr.t

!=============================================================================
subroutine initglobaldata_usr

use mod_global_parameters

integer :: fileid

character(79) :: filename
!--------------------------------------------------------------------------

eqpar(gamma_)=5.0d0/3.0d0

select case( iprob )
 case(1) ! O-star in cold medium
   eqpar(Mdot1_)  = 1.0d-6*msol/year 
   eqpar(vwind1_) = 2.0d8
   eqpar(Twind1_) = 1.0d4
   eqpar(rhoISM_) = (10.0d0)**(-22.5)
   eqpar(vISM_)   = 5.0d6  
   eqpar(TISM_)   = 1.0d2
   eqpar(Rstar_)  = 5.0d13
   eqpar(Rwind_)  = 2.0d-1
 case(2)  ! RSG in cold medium  
   eqpar(Mdot1_)  = 1.0d-4*msol/year 
   eqpar(vwind1_) = 2.0d6
   eqpar(Twind1_) = 1.0d3
   eqpar(rhoISM_) = (10.0d0)**(-23.5)
   eqpar(vISM_)   = 5.0d6  
   eqpar(TISM_)   = 1.0d2
   eqpar(Rstar_)  = 5.0d13
   eqpar(Rwind_)  = 2.0d-1
 case(3)  ! WR in cold medium  
   eqpar(Mdot1_)  = 1.0d-5*msol/year 
   eqpar(vwind1_) = 2.5d8
   eqpar(Twind1_) = 2.0d4
   eqpar(rhoISM_) = (10.0d0)**(-22.5)
   eqpar(vISM_)   = 5.0d6  
   eqpar(TISM_)   = 1.0d2
   eqpar(Rstar_)  = 5.0d13
   eqpar(Rwind_)  = 2.0d-1
 case default
     call mpistop("This problem has not been defined")
end select 

!==========================================================================
!= Above this line nothing is normalized
!==========================================================================

length_convert_factor    = 3.0857D18
w_convert_factor(rho_) = 10.0**(-25)
w_convert_factor(mom(1))  = 1.0d7
{^NOONEC w_convert_factor(mom(2))    = w_convert_factor(mom(1))}
w_convert_factor(p_)   = w_convert_factor(rho_)*w_convert_factor(mom(1))*w_convert_factor(mom(1))
time_convert_factor         = length_convert_factor/w_convert_factor(mom(1))

eqpar(Rstar_) = eqpar(Rstar_) / length_convert_factor
eqpar(Tscale_) = (1.0D0/(w_convert_factor(mom(1))**2.0d0)) &
               * kboltz/mhydro
eqpar(Lscale_) =  w_convert_factor(rho_)*time_convert_factor/((mhydro*w_convert_factor(mom(1)))**2.0)
eqpar(Mue_) = 1.0D0

if(mype == 0) then
   fileid = 11
   write(filename,"(a,a)") TRIM(filenamelog),".scale"
   open(fileid,file=filename, status='unknown')
   write(fileid,1004) 'time_convert_factor:        ', time_convert_factor
   write(fileid,1004) 'length_convert_factor:   ', length_convert_factor
   write(fileid,1004) 'w_convert_factor(mom(1)): ', w_convert_factor(mom(1))
   write(fileid,1004) 'w_convert_factor(rho_):', w_convert_factor(rho_)    
   write(fileid,1004) 'w_convert_factor(p_):  ', w_convert_factor(p_)    
   write(fileid,*) 
   write(fileid,1004) 'accel         ', w_convert_factor(mom(1))*w_convert_factor(mom(1))/length_convert_factor
   write(fileid,*)
   write(fileid,1002) 1.0d0/eqpar(Tscale_)
   write(fileid,1003) eqpar(Lscale_)
   write(fileid,*) 
   close(fileid)
endif

call coolinit

tlow=eqpar(TISM_)*eqpar(Tscale_)

1001 format(1x,i4,1x,1pe12.5)
1002 format('Temperature unit: ', 1x1pe12.5)
1003 format('Luminosity scale: ', 1x1pe12.5)
1004 format(a14,1x,1pe12.5)
 
end subroutine initglobaldata_usr
!===============================================================================
subroutine initonegrid_usr(ixG^L,ix^L,w,x)

! initialize one grid 

use mod_global_parameters
integer, intent(in) :: ixG^L, ix^L

double precision, intent(in) :: x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)

logical :: interpw(1:nw)
double precision :: rad(ixG^T), cosTh(ixG^T), sinTh(ixG^T)

logical :: patchw(ixG^T)
!--------------------------------------------------------------------------

rad(ix^S) = dsqrt(x(ix^S,1)**2 + x(ix^S,2)**2)
cosTh(ix^S) = x(ix^S,2)/rad(ix^S)
sinTh(ix^S) = x(ix^S,1)/rad(ix^S)

where ( rad(ix^S)>= eqpar(Rwind_) ) 
      w(ix^S,rho_) = eqpar(rhoISM_)/w_convert_factor(rho_)
      w(ix^S,mom(1))  = zero
{^NOONEC                      
       w(ix^S,mom(2))  = -eqpar(vISM_)  /  w_convert_factor(mom(2))
}
      w(ix^S,p_)   = w(ix^S,rho_)*eqpar(Tism_)*eqpar(Tscale_)
   elsewhere
      w(ix^S,rho_) = eqpar(Mdot1_)/(4.0D0*dpi*eqpar(vwind1_) * (rad(ix^S)*length_convert_factor)**2 ) &
                   / w_convert_factor(rho_)
      w(ix^S,mom(1))  = (eqpar(vwind1_) /w_convert_factor(mom(1))) * sinTh(ix^S)
{^NOONEC
      w(ix^S,mom(2))  = (eqpar(vwind1_) /w_convert_factor(mom(2))) * cosTh(ix^S)
}
      w(ix^S,p_)     =  w(ix^S,rho_)*eqpar(Twind1_)* eqpar(Tscale_)
end where

patchw(ix^S)=.false.
call conserve(ixG^L,ix^L,w,x,patchw)

return
end subroutine initonegrid_usr
!===============================================================================
subroutine specialbound_usr(qt,ixG^L,ixO^L,iw,iB,w,x)

! special boundary types, user defined

use mod_global_parameters

integer, intent(in) :: ixG^L, ixO^L, iw, iB 
double precision, intent(in) :: qt,x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)

logical :: patchw(ixG^T)
!-----------------------------------------------------------------------------
 
select case(iB)
case(4)
      w(ixO^S,rho_)   = eqpar(rhoISM_)/w_convert_factor(rho_)
      w(ixO^S,mom(1))    = zero
      w(ixO^S,mom(2))    = -eqpar(vISM_)/w_convert_factor(mom(2))
      w(ixO^S,p_)   = w(ixO^S,rho_)*eqpar(TISM_)*eqpar(Tscale_)
      patchw(ixO^S)=.false.
      call conserve(ixG^L,ixO^L,w,x,patchw)
case default    
  call mpistop("This boundary is not supposed to be special")
end select

return
end subroutine specialbound_usr
!=============================================================================
subroutine specialsource(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

! Calculate w(iw)=w(iw)+qdt*SOURCE[wCT,qtC,x] within ixO for all indices
! iw=iwmin...iwmax.  wCT is at time qCT

use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L, iw^LIM
double precision, intent(in) :: qdt, qtC, qt, x(ixI^S,1:ndim)
double precision, intent(inout) :: wCT(ixI^S,1:nw), w(ixI^S,1:nw)

double precision :: rad(ixG^T), cosTh(ixG^T), sinTh(ixG^T), Rw
!-----------------------------------------------------------------------------

Rw = eqpar(Rwind_)
rad(ixO^S) = dsqrt(x(ixO^S,1)**2 + x(ixO^S,2)**2)
cosTh(ixO^S) = x(ixO^S,2)/rad(ixO^S)
sinTh(ixO^S) = x(ixO^S,1)/rad(ixO^S)

where ( rad(ixO^S)< Rw ) 
      w(ixO^S,rho_)  = eqpar(Mdot1_)/(4.0D0*dpi*eqpar(vwind1_)* (rad(ixO^S)*length_convert_factor)**2 ) &
                      / w_convert_factor(rho_)
      w(ixO^S,mom(1))   = (eqpar(vwind1_) /w_convert_factor(mom(1))) *sinTh(ixO^S)*w(ixO^S,rho_)
{^NOONEC
      w(ixO^S,mom(2))   = (eqpar(vwind1_) /w_convert_factor(mom(1))) *cosTh(ixO^S)*w(ixO^S,rho_)
}
      w(ixO^S,e_)    = w(ixO^S,rho_)*eqpar(Twind1_)*eqpar(Tscale_)/(eqpar(gamma_)-one)+ &
           half*(^C&w(ixO^S,m^C_)**2.0d0+)/w(ixO^S,rho_)
end where

call addsource_cooling(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

return
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

call getdt_cooling(w,ixG^L,ix^L,dtnew,dx^D,x)

return
end subroutine getdt_special
!=============================================================================
subroutine specialeta(w,ixI^L,ix^L,idirmin,x,current,eta)

! Set the "eta" array for resistive MHD based on w or the
! "current" variable which has components between idirmin and 3.

use mod_global_parameters

integer, intent(in) :: ixI^L, ix^L, idirmin
double precision, intent(in) :: w(ixI^S,nw), x(ixI^S,1:ndim)

double precision :: current(ixG^T,7-2*ndir:3), eta(ixG^T)
!-----------------------------------------------------------------------------

call mpistop("specialeta is not defined")

end subroutine specialeta

!=============================================================================
subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)

! this subroutine can be used in convert, to add auxiliary variables to the
! converted output file, for further analysis using tecplot, paraview, ....

use mod_global_parameters

integer, intent(in)                :: ixI^L,ixO^L
double precision, intent(in)       :: x(ixI^S,1:ndim)
double precision                   :: w(ixI^S,nw+nwauxio)
double precision                   :: normconv(0:nw+nwauxio)

double precision :: wloc(ixG^T,1:nw),tmp(ixG^T)
!-----------------------------------------------------------------------------

wloc(ixI^S,1:nw)=w(ixI^S,1:nw)
if(saveprim)then
   tmp(ixO^S)=wloc(ixO^S,p_)
 else
   call getpthermal(wloc,x,ixI^L,ixO^L,tmp)
endif

if(iprob==-1)then
  w(ixO^S,nw+1)=dlog10(wloc(ixO^S,rho_))
  w(ixO^S,nw+2)=tmp(ixO^S)/wloc(ixO^S,rho_)
else
  w(ixO^S,nw+1)=dlog10(wloc(ixO^S,rho_)*w_convert_factor(rho_))
  ! output the temperature p/rho
  w(ixO^S,nw+2)=(tmp(ixO^S)/wloc(ixO^S,rho_))/eqpar(Tscale_)
endif


end subroutine specialvar_output
!=============================================================================
subroutine specialvarnames_output

! newly added variables need to be concatenated with the varnames/primnames string

use mod_global_parameters
!-----------------------------------------------------------------------------

if(iprob==-1)then
  primnames= TRIM(primnames)//' '//'logrhoN TN'
  w_names=TRIM(w_names)//' '//'logrhoN TN'
else
  primnames= TRIM(primnames)//' '//'logrho T'
  w_names=TRIM(w_names)//' '//'logrho T'
endif

end subroutine specialvarnames_output
!=============================================================================
subroutine specialrefine_grid(igrid,level,ixG^L,ix^L,qt,w,x,refine,coarsen)

! Enforce additional refinement or coarsening
! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.

use mod_global_parameters

integer, intent(in) :: igrid, level, ixG^L, ix^L
double precision, intent(in) :: qt, w(ixG^S,1:nw), x(ixG^S,1:ndim)
integer, intent(inout) :: refine, coarsen

double precision:: Rw, rad(ixG^T)
!-----------------------------------------------------------------------------

Rw = eqpar(Rwind_)
rad(ix^S) = dsqrt(x(ix^S,1)**2 + x(ix^S,2)**2)
if( any(rad(ix^S) <= 1.5*Rw) ) then
      coarsen = -1
      refine  = 1
endif

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

if (iflag >nw+1)call mpistop(' iflag error')
var(ixI^S) = dsqrt(w(ixI^S,mom(1))**2 + w(ixI^S,mom(2))**2 )/w(ixI^S,rho_)

end subroutine specialvarforerrest
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
subroutine specialset_B0(ixI^L,ixO^L,x,wB0)

! Here one can add a steady (time-independent) potential background field

use mod_global_parameters

integer, intent(in)           :: ixI^L,ixO^L
double precision, intent(in)  :: x(ixG^T,1:ndim)
double precision, intent(inout) :: wB0(ixI^S,1:ndir)
!-----------------------------------------------------------------------------
call mpistop(' abs(Busr)> 0, make change in parfile or in user file')

wB0(ixO^S,1:ndir)=wB0(ixO^S,1:ndir)
!!wB0(ixO^S,1:ndir)=wB0(ixO^S,1:ndir)+user defined steady potential field1

end subroutine specialset_B0
!=============================================================================
subroutine bc_int(qt,ixG^L,ixO^L,w,x)

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

use mod_global_parameters

integer, intent(in) :: ixG^L,ixO^L
double precision, intent(in) :: qt
double precision, intent(inout) :: w(ixG^S,1:nw)
double precision, intent(in) :: x(ixG^S,1:ndim)
!----------------------------------------------------------------------------

call mpistop("bc_int not defined")

end subroutine bc_int
!=============================================================================
! amrvacusr.t.movstar
!=============================================================================
