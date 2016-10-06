!=============================================================================
! amrvacusr.t.raytest
! Project : testing raytracing module
! 
! update : 27/0/2012
! made by Allard 
! configuration :
! setamrvac -d=22, -phi=0 -z=0 -g=12,12 -p=hd -u=ray_pne


! INCLUDE:amrvacnul/specialini.t
!INCLUDE:amrvacnul/speciallog.t
INCLUDE:amrvacnul/specialimpl.t
!INCLUDE:amrvacnul/specialsource.t
INCLUDE:amrvacmodules/cooling.t
INCLUDE:amrvacmodules/raytracing.t

!=============================================================================
subroutine initglobaldata_usr



use mod_raytracing

use mod_global_parameters

integer :: i, fileid
double precision, allocatable :: evol(:,:)

double precision :: Rstar, Mstar

character(79) :: filename,fname


!-----------------------------------------------------------------

eqpar(gamma_) = 5.0D0/3.0D0


eqpar(rhoISM_) = 1.0d5*mhydro
eqpar(TISM_)   = 1.0d2

eqpar(Fstar_) =  1.0d46 ! 0.0d0      


eqpar(alphah_)  = 1.0d-13  ! 2.0d-13! 1.0d-13 

 eqpar(Mdot1_) = 1.0d-5 * msol/year
 eqpar(Mdot2_) = 1.0d-7 * msol/year
! eqpar(Mdot2_) = 1.0d-7 * msol/year
 
 eqpar(vel1_) = 1.0d6
 eqpar(vel2_) = 1.0d8


!------------------------------------------------------------------

{^IFTWOD
 nrays          =  1152 ! 768 !384 ! 192  ! 96   ! total number of rays
! nrays2         = nrays
! nrays1         = 1
! nrays3         = 1
 npoints        =  769 !385 ! 193 ! 97    ! number of points per ray
 nu             =  1              ! number of ray-specific variables
 raysource      = 'x1dir'
 raytype        = 'fixed'
 regular_rays   = .true.          ! true if raysource is x1dir, x2dir or x3dir
                                  ! overrides the raytogrid setting
 gridtoray      = 'nearest'       ! how to interpolate grid to ray
 raytogrid      = 'nearestray'    ! how to interpolate ray to grid
 raycomm        = 'distribute_block'
 writerays      = .false.          ! true if you want to output the rays
 quicksearch    = 'basic'         ! check if rays go anywhere near the grid
 checkintersect = .true.          ! see quicksearch
 rayphysics     = 'stromgren'     !   'columndens'
 ion_           = 1               ! cdens_         = 1
}

{^IFTHREED 
 nrays          = 2304 ! total number of rays
! nrays2         = 48
! nrays1         = 1
! nrays3         = nrays/nrays2
 npoints        =  193 ! 769 !385 ! 193 ! 97    ! number of points per ray
 nu             =  1              ! number of ray-specific variables
 raysource      = 'x1dir'
 raytype        = 'fixed'
 regular_rays   = .false.          ! true if raysource is x1dir, x2dir or x3dir
                                  ! overrides the raytogrid setting
 gridtoray      = 'nearest'       ! how to interpolate grid to ray
 raytogrid      = 'nearest'    ! how to interpolate ray to grid
 raycomm        = 'distribute_block'
 writerays      = .true.          ! true if you want to output the rays
 quicksearch    = 'basic'         ! check if rays go anywhere near the grid
 checkintersect = .true.          ! see quicksearch
 rayphysics     = 'stromgren'     !   'columndens'
 ion_           = 1               ! cdens_         = 1 
} 
 
 
if (.not. (fixprocess) ) call MPISTOP( 'getting the grid data to the rays would be a good idea!' )

!-----------------------------------------------------------------------------

!==========================================================================
!= Above this line nothing is normalized
!==========================================================================


select case( iprob)
 case( 1 )
   normvar(0)    = 3.0857D17
   normvar(rho_) = 10.0**(-18.5)
   normvar(v1_)  = 1.0d6
   normvar(p_)   = normvar(rho_)*normvar(v1_)*normvar(v1_)

{^NOONEC normvar(v2_)    = normvar(v1_)}
{^IFTHREEC normvar(v3_)  = normvar(v1_)}
   normt           = normvar(0)/normvar(v1_)
 case( 2 )
   normvar(0)    = 3.0857D17
   normvar(rho_) = 10.0**(-20.0)
   normvar(v1_)  = 1.0d6
   normvar(p_)   = normvar(rho_)*normvar(v1_)*normvar(v1_)

{^NOONEC normvar(v2_)    = normvar(v1_)}
{^IFTHREEC normvar(v3_)  = normvar(v1_)}
   normt           = normvar(0)/normvar(v1_)

end select

eqpar(Tscale_) = (1.0D0/(normvar(v1_)**2.0d0)) &
               * kboltz/mhydro
eqpar(Lscale_) =  normvar(rho_)*normt/((mhydro*normvar(v1_))**2.0)

eqpar(Mue_) = 1.0d0

!eqpar(Fstar_) = (eqpar(Fstar_)/(4.0d0*dpi*eqpar(alphah_)) ) / (normvar(0)**3) ! Not sure about this part yet
eqpar(Fstar_) = (eqpar(Fstar_)/(4.0d0*dpi) ) / (normvar(0)**3)


if(mype == 0) then
   fileid = 11
   write(filename,"(a,a)") TRIM(filenamelog),".scale"
   open(fileid,file=filename, status='unknown')
!   do i=1,nspecialpar
!      write(fileid,1001) neqpar+i, eqpar(neqpar+i)
!   end do
   write(fileid,1004) 'normt:        ', normt
   write(fileid,1004) 'normvar(0):   ', normvar(0)
   write(fileid,1004) 'normvar(v1_): ', normvar(v1_)
   write(fileid,1004) 'normvar(rho_):', normvar(rho_)    
   write(fileid,1004) 'normvar(p_):  ', normvar(p_)    
   write(fileid,*)
   write(fileid,1002) eqpar(Tscale_)
   write(fileid,1003) eqpar(Lscale_)

   close(fileid)
endif


 call random_reseed
  

 call coolinit
 Tlow = 1.0d2*eqpar(Tscale_)


1001 format(1x,i4,1x,1pe12.5)
1002 format('Temperature scale: ', 1x,1pe12.5)
1003 format('Luminosity scale: ', 1x,1pe12.5)
1004 format(a14,1x,1pe12.5)
return
end subroutine initglobaldata_usr
!===============================================================================
subroutine initonegrid_usr(ixG^L,ix^L,w,x)
use mod_raytracing

! initialize one grid 

use mod_global_parameters
integer, intent(in) :: ixG^L, ix^L

double precision, intent(in) :: x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)

logical :: interpw(1:nw)

logical :: patchw(ixG^T)

integer :: ix^D
real    :: rnum

!--------------------------------------------------------------------------


   where( x(ix^S,1) > xprobmin1*2.0d0 )
 
      w(ix^S,rho_) = (eqpar(Mdot1_)/(4.0d0*dpi*eqpar(vel1_)*((x(ix^S,1)*normvar(0))**2)) )/ normvar(rho_)
      w(ix^S,v1_)  = eqpar(vel1_)/normvar(v1_)

{^NOONEC
      w(ix^S,v2_)  = zero   
}
{^IFTHREEC
      w(ix^S,v3_)  = zero   
}
      w(ix^S,e_)   = w(ix^S,rho_) * 1.0d2 *eqpar(Tscale_)

   elsewhere
      w(ix^S,rho_) = (eqpar(Mdot2_)/(4.0d0*dpi*eqpar(vel2_)*((x(ix^S,1)*normvar(0))**2)) )/ normvar(rho_)
      w(ix^S,v1_)  = eqpar(vel2_)/normvar(v1_)

{^NOONEC
      w(ix^S,v2_)  = zero   
}
{^IFTHREEC
      w(ix^S,v3_)  = zero   
}
      w(ix^S,e_)   =  w(ix^S,rho_) * 1.0d4 *eqpar(Tscale_)
     
   end where


!{^NOONED
!{do ix^D = ix^LIM^D\}
!  call random_number(rnum)
!  if(rnum > half ) then
!      w(ix^D,1:nw) = w(ix^D,1:nw)*1.01d0
!  else if(rnum< half) then
!    w(ix^D,1:nw) = w(ix^D,1:nw)*0.99d0
!  endif
!{enddo^D&\}
!}



 patchw(ix^S)=.false.
 call conserve(ixG^L,ix^L,w,x,patchw)




return
end subroutine initonegrid_usr
!===============================================================================
subroutine specialbound_usr(qt,ixG^L,ixO^L,iw,iB,w,x)

! special boundary types, user defined

use mod_global_parameters

double precision, parameter :: gam1=(5.0d0/3.0d0)-1.0d0


integer, intent(in) :: ixG^L, ixO^L, iw, iB 
double precision, intent(in) :: qt,x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)

integer :: rightnei(1),leftnei(1),righti,lefti
double precision :: dmdt_t,vwind_t, Twind_t, Omega_t, dmdisk_t
double precision :: dtime,deltatime


logical :: patchw(ixG^T)

integer :: ix^D

double precision :: mdot, vwind, twind
real    :: rnum(ixG^S)

!-----------------------------------------------------------------------------

select case(iB)
case(1)
    w(ixO^S,rho_) = eqpar(Mdot2_)/(4.0D0*dpi*eqpar(vel2_)*(x(ixO^S,1)*normvar(0))**2)& 
    / normvar(rho_) 
    w(ixO^S,v1_)  = eqpar(vel2_)/normvar(v1_)
{^NOONEC
    w(ixO^S,v2_)  = zero
}
{^IFTHREEC
    w(ixO^S,v3_)  = zero   
}
    !      w(ixO^S,p_)   =  w(ixO^S,rho_)*1.0d4*eqpar(Tscale_)
    w(ixO^S,p_)   =  w(ixO^S,rho_)*2.0d4*eqpar(Tscale_)

    patchw(ixO^S)=.false.
    call conserve(ixG^L,ixO^L,w,x,patchw)                
end select




return

end subroutine specialbound_usr
!=============================================================================
subroutine specialsource(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

! Calculate w(iw)=w(iw)+qdt*SOURCE[wCT,qtC,x] within ixO for all indices
! iw=iwmin...iwmax.  wCT is at time qCT
use mod_raytracing

use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L, iw^LIM
double precision, intent(in) :: qdt, qtC, qt, x(ixI^S,1:ndim)
double precision, intent(inout) :: wCT(ixI^S,1:nw), w(ixI^S,1:nw)

integer :: ix^D
double precision :: tau(ixG^T),ion(ixG^T)

logical :: patchw(ixG^T)

integer :: iray, iraynear
logical :: patchray(1:nrays)
double precision  :: xpoint^D
double precision  :: rtoray(1:nrays),rmin          



!-----------------------------------------------------------------------------

   call addsource_cooling(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

   ion(ixO^S) = zero
   call primitive(ixI^L,ixO^L,w,x)
      
   call getraydata( ixI^L,ixO^L,x,ion ) 
   
{do ix^D = ixO^LIM^D\}
   if( ion(ix^D ) > 0.5d0 ) then
      w(ix^D,p_) = max( w(ix^D,p_),w(ix^D,rho_)*1.0d4*eqpar(Tscale_) )
   endif
{enddo^D&\}

    patchw(ixO^S)=.false.
    call conserve(ixI^L,ixO^L,w,x,patchw)


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

double precision :: dt_tdtime,deltatime
integer :: rightnei(1),leftnei(1),righti,lefti

!-----------------------------------------------------------------------------

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

!  eta(ix^S)=...

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
double precision                   :: tmp(ixG^T)

double precision                   :: divb(ixG^T),cmax(ixG^T)
integer                            :: jxO^L,ixOO^L,idims,ix^D
double precision                   :: wtmp(ixG^T,1:nw)
!-----------------------------------------------------------------------------



end subroutine specialvar_output

!=============================================================================
subroutine specialvarnames_output

! newly added variables need to be concatenated with the varnames/primnames string

use mod_global_parameters
!-----------------------------------------------------------------------------
oktest = index(teststr,'printlog')>=1

!call mpistop("special varnames and primnames undefined")

! Example : as above in specialvar_output, assuming relativistic HD here...



end subroutine specialvarnames_output
!=============================================================================
subroutine specialrefine_grid(igrid,level,ixG^L,ix^L,qt,w,x,refine,coarsen)

! Enforce additional refinement or coarsening
! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.

use mod_global_parameters

integer, intent(in) :: igrid, level, ixG^L, ix^L
double precision, intent(in) :: qt, w(ixG^S,1:nw), x(ixG^S,1:ndim)
integer, intent(inout) :: refine, coarsen
!-----------------------------------------------------------------------------

!if( any(x(ixG^S,1) <= xprobmin1*1.05) ) then
if( any(x(ixG^S,1) <= (xprobmin1*1.01)) ) then
if( level <= mxnest ) then
   coarsen = -1
   refine  = 1

endif 
{^IFTHREED
else
 coarsen = 1
 refine = -1
}
endif




return
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

if (iflag >nw+1)call mpistop(' iflag> nw, make change in parfile or in user file')
{^IFTWOC
var(ixI^S) = dsqrt(w(ixI^S,m1_)**2 + w(ixI^S,m2_)**2 )/w(ixI^S,rho_)
}
{^IFTHREEC
var(ixI^S) = dsqrt(w(ixI^S,m1_)**2 + w(ixI^S,m2_)**2 +  w(ixI^S,m3_)**2)/w(ixI^S,rho_)
}

end subroutine specialvarforerrest
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

 call  getgriddata(ixI^L,ixO^L,x,w)

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
!call conserve(ixG^L,ixO^L,w,x,patchw)

end subroutine bc_int
!=============================================================================
subroutine random_reseed
!
! Gives new seeding every timestep to avoid 
! resusing random number sequence.
!

use mod_global_parameters

integer, dimension(:), allocatable :: seed
integer :: kseed, clock, i
!-----------------------------------------------------------------------------

call random_seed(size=kseed)

allocate(seed(1:kseed))

call system_clock(COUNT=clock)

seed=clock + (mype+1)*(/ (i-1,i=1,kseed) /)

call random_seed(put=seed(1:kseed))



deallocate(seed)

return
end subroutine random_reseed
!=============================================================================
!=============================================================================


subroutine fixp_usr(ixI^L,ixO^L,w)
use mod_global_parameters

integer, intent(in)                :: ixI^L,ixO^L
double precision, intent(inout)    :: w(ixI^S,1:nw)
logical :: patchp(ixI^S)
!----------------------------------------------------------------------------
!patchp(ixO^S) = .true.


end subroutine fixp_usr
!=============================================================================
      subroutine flag_grid_usr(qt,ixG^L,ixO^L,w,x,flag)

      use mod_global_parameters

      integer, intent(in)             :: ixG^L, ixO^L
      integer, intent(inout)          :: flag
      double precision, intent(in)    :: qt
      double precision, intent(inout) :: w(ixG^S,1:nw)
      double precision, intent(in)    :: x(ixG^S,1:ndim)
double precision                :: ro, ri

      ! flag=-1 : Treat all cells active, omit deactivation (onentry, default)
      ! flag=0  : Treat as normal domain
      ! flag=1  : Treat as passive, but reduce by safety belt
      ! flag=2  : Always treat as passive

!-----------------------------------------------------------------------------

!flag=0
! march out with r = v * t
!ro = (eqpar(vel2_)/normvar(v1_) * qt + xprobmin1*10.0d0)*zero +xprobmax1*2.0d0
!ri = eqpar(vel1_)/normvar(v1_) * qt + xprobmin1-1e-3
!if ( .not. any(x(ixO^S,1)<ro) .or. any(x(ixO^S,1) < ri)) then
!   flag = 1
!end if

      end subroutine flag_grid_usr
!=============================================================================

!=============================================================================
! amrvacusr.t.ramram2
!=============================================================================


