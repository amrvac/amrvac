!=============================================================================
! amrvacusr.t.movstar
! Demonstration for moving star model
!
! update : 11/2014
! made by Allard 
! configuration :
! setamrvac -d=22, -z=2, -phi =0, -g=24,24 -p=hd 

! INCLUDE:amrvacnul/specialini.t
!INCLUDE:amrvacnul/speciallog.t
INCLUDE:amrvacnul/specialimpl.t
!INCLUDE:amrvacnul/specialsource.t
INCLUDE:amrvacmodules/cooling.t

!=============================================================================
subroutine initglobaldata_usr

use mod_global_parameters

integer :: i, fileid
double precision, allocatable :: evol(:,:)
character(79) :: filename,fname


eqpar(gamma_) = 5.0D0/3.0D0

select case ( iprob )

case( 0 ) ! No movement

eqpar(Mdot1_)  = 1.0d-6*msol/year 
eqpar(vwind1_) = 1.0d8
eqpar(Twind1_) = 1.0d4
eqpar(rhoISM_) = (10.0d0)**(-23.5)
eqpar(vISM_)   = 0.0d0
eqpar(TISM_)   = 1.0d2
!
eqpar(Rstar_)  = 5.0d13
eqpar(Rwind_)  = 2.0d-2*xprobmax1

case( 1) ! O-star

eqpar(Mdot1_)  = 1.0d-6*msol/year 
eqpar(vwind1_) = 2.0d8
eqpar(Twind1_) = 1.0d4
eqpar(rhoISM_) = (10.0d0)**(-23.5)
eqpar(vISM_)   = 5.0d6  
eqpar(TISM_)   = 1.0d2
!
eqpar(Rstar_)  = 5.0d13
eqpar(Rwind_)  = 2.0d-2* xprobmax1
 
case(2) ! RGB star

eqpar(Mdot1_)  = 1.0d-6*msol/year 
eqpar(vwind1_) = 1.5d6
eqpar(Twind1_) = 1.0d2
eqpar(rhoISM_) = (10.0d0)**(-23.5)
eqpar(vISM_)   = 2.5d6  
eqpar(TISM_)   = 1.0d2
!
eqpar(Rstar_)  = 5.0d13
eqpar(Rwind_)  = 2.0d-2* xprobmax1

case default
 call mpistop( "This problem is not defined" )
end select


!==========================================================================
!= Above this line nothing is normalized
!==========================================================================



normvar(0)      = 3.0857D18
normvar(rho_)   = 10.0**(-25)
normvar(v1_)    = 1.0d7


normvar(p_)   = normvar(rho_)*normvar(v1_)*normvar(v1_)
{^NOONEC normvar(v2_)    = normvar(v1_)}
{^IFTHREEC normvar(v3_)  = normvar(v1_)}
normt           = normvar(0)/normvar(v1_)

eqpar(Rstar_) = eqpar(Rstar_) / normvar(0)

eqpar(Tscale_) = (1.0D0/(normvar(v1_)**2.0d0)) &
               * kboltz/mhydro
eqpar(Lscale_) =  normvar(rho_)*normt/((mhydro*normvar(v1_))**2.0)
               
eqpar(Mue_) = 1.0D0




if(mype == 0) then
   fileid = 11
   write(filename,"(a,a)") TRIM(filenamelog),".scale"
   open(fileid,file=filename, status='unknown')
!   do i=1,nspecialpar
!      write(fileid,1001) neqpar+i, eqpar(neqpar+i)
!   end do
   write(fileid,1004) 'normt: zero   !       ', normt
   write(fileid,1004) 'normvar(0):   ', normvar(0)
   write(fileid,1004) 'normvar(v1_): ', normvar(v1_)
   write(fileid,1004) 'normvar(rho_):', normvar(rho_)    
   write(fileid,1004) 'normvar(p_):  ', normvar(p_)    
   write(fileid,*) 
   write(fileid,1004) 'accel         ', normvar(v1_)*normvar(v1_)/normvar(0)
   write(fileid,*)
 
   write(fileid,1002) eqpar(Tscale_)
   write(fileid,1003) eqpar(Lscale_)
   write(fileid,*) 
   close(fileid)
endif

call coolinit

tlow=eqpar(TISM_)*eqpar(Tscale_)

1001 format(1x,i4,1x,1pe12.5)
1002 format('Temperature scale: ', 1x1pe12.5)
1003 format('Luminosity scale: ', 1x1pe12.5)
1004 format(a14,1x,1pe12.5)
return
end subroutine initglobaldata_usr
!===============================================================================
subroutine initonegrid_usr(ixG^L,ix^L,w,x)

! initialize one grid 

use mod_global_parameters
integer, intent(in) :: ixG^L, ix^L

double precision, intent(in) :: x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)

double precision :: rad(ixG^S), cosTh(ix^S), sinTh(ix^S)

integer :: idust

logical :: patchw(ixG^T)
!-----------------------------------------------------------------------------


rad(ix^S) = dsqrt(x(ix^S,1)**2 + x(ix^S,2)**2)
cosTh(ix^S) = x(ix^S,2)/rad(ix^S)
sinTh(ix^S) = x(ix^S,1)/rad(ix^S)

   where ( rad(ix^S)>= eqpar(Rwind_) ) 
      w(ix^S,rho_) = eqpar(rhoISM_)/normvar(rho_)
      w(ix^S,v1_)  = zero
{^NOONEC                      
       w(ix^S,v2_)  = -eqpar(vISM_)  /  normvar(v2_)
}
{^IFTHREEC
      w(ix^S,v3_)  = zero
}
      w(ix^S,p_)   = w(ix^S,rho_)*eqpar(Tism_)*eqpar(Tscale_)

   elsewhere
      w(ix^S,rho_) = eqpar(Mdot1_)/(4.0D0*dpi*eqpar(vwind1_) * (rad(ix^S)*normvar(0))**2 ) &
                   / normvar(rho_)
      w(ix^S,v1_)  = (eqpar(vwind1_) /normvar(v1_)) * sinTh(ix^S)
{^NOONEC
      w(ix^S,v2_)  = (eqpar(vwind1_) /normvar(v2_)) * cosTh(ix^S)
}

{^IFTHREEC
      w(ix^S,v3_)  = zero   
}
!     Assume that at this point Tdust=Tgas; Probably not too far off.

      w(ix^S,p_)     =  w(ix^S,rho_)*eqpar(Twind1_)* eqpar(Tscale_)
   end where




patchw(ix^S)=.false.
call conserve(ixG^L,ix^L,w,patchw)




if(loglimit(m1_)) call mpistop("loglimit of vector quantity is a bad idea!")
{^NOONEC
if(loglimit(m2_)) call mpistop("loglimit of vector quantity is a bad idea!")
}
{^IFTHREEC
if(loglimit(m3_)) call mpistop("loglimit of vector quantity is a bad idea!")
}





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


logical :: patchw(ixG^T)

integer :: idust

!-----------------------------------------------------------------------------


select case(iB)
case(4)
      w(ixO^S,rho_)   = eqpar(rhoISM_)/normvar(rho_)
      w(ixO^S,v1_)    = zero
 {^NOONEC     
      w(ixO^S,v2_)    = -eqpar(vISM_)/normvar(v2_)
  }    
{^IFTHREEC
      w(ixO^S,v3_)    = zero
}
      w(ixO^S,p_)   = w(ixO^S,rho_)*eqpar(TISM_)*eqpar(Tscale_)
      patchw(ixO^S)=.false.
      call conserve(ixG^L,ixO^L,w,x,patchw)
end select





return

2001 format(5(1x,1pe12.5))
end subroutine specialbound_usr
!=============================================================================
subroutine specialsource(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

! Calculate w(iw)=w(iw)+qdt*SOURCE[wCT,qtC,x] within ixO for all indices
! iw=iwmin...iwmax.  wCT is at time qCT

use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L, iw^LIM
double precision, intent(in) :: qdt, qtC, qt, x(ixI^S,1:ndim)
double precision, intent(inout) :: wCT(ixI^S,1:nw), w(ixI^S,1:nw)

double precision :: rad(ixI^S), cosTh(ixI^S), sinTh(ixI^S), Rw

!-----------------------------------------------------------------------------

Rw = eqpar(Rwind_)

select case( iprob )
 case( 1 ) 
 if (it==0) then 
    Rw = eqpar(Rwind_)
 else
    Rw = eqpar(Rwind_)* min( 3.0d0, (one + qtC*normt/(1.0d4*year)) )
 endif 
end select




rad(ixO^S) = dsqrt(x(ixO^S,1)**2 + x(ixO^S,2)**2)
cosTh(ixO^S) = x(ixO^S,2)/rad(ixO^S)
sinTh(ixO^S) = x(ixO^S,1)/rad(ixO^S)


where ( rad(ixO^S)< Rw ) 
      w(ixO^S,rho_)   = eqpar(Mdot1_)/(4.0D0*dpi*eqpar(vwind1_)* (rad(ixO^S)*normvar(0))**2 ) &
                      / normvar(rho_)
      w(ixO^S,m1_)    = (eqpar(vwind1_) /normvar(v1_)) *sinTh(ixO^S)*w(ixO^S,rho_)
{^NOONEC
      w(ixO^S,m2_)    = (eqpar(vwind1_) /normvar(v1_)) *cosTh(ixO^S)*w(ixO^S,rho_)
}

{^IFTHREEC
      w(ixO^S,m3_)    = zero   
}

      w(ixO^S,p_)     = w(ixO^S,rho_)*eqpar(Twind1_)*eqpar(Tscale_)
       w(ixO^S,e_)    = w(ixO^S,p_)/(eqpar(gamma_)-one)+ &
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

double precision :: dt_tdtime,deltatime
integer :: rightnei(1),leftnei(1),righi,lefti

!-----------------------------------------------------------------------------

 call getdt_cooling(w,ixG^L,ix^L,dtnew,dx^D,x)

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
subroutine specialrefine_grid(igrid,level,ixG^L,ix^L,qt,w,x,my_refine,my_coarsen)

! Enforce additional refinement or coarsening
! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.

use mod_global_parameters

integer, intent(in) :: igrid, level, ixG^L, ix^L
double precision, intent(in) :: qt, w(ixG^S,1:nw), x(ixG^S,1:ndim)
integer, intent(inout) :: my_refine, my_coarsen

double precision :: rad(ixG^S)
!-----------------------------------------------------------------------------

rad(ix^S) = dsqrt(x(ix^S,1)**2 + x(ix^S,2)**2)


     if( any(rad(ix^S) <= (eqpar(Rwind_)*1.5d0)) ) then ! Max resolution in wind initialization zone
	my_refine  = 1
        my_coarsen = -1
     else if( any(dabs(x(ix^S,2)) > dabs(xprobmax1)) ) then ! Lower resolution in the tail
        if( level > 4 ) then
	 my_refine  = -1
         my_coarsen = 1
        else if( level==4 ) then
 	 my_refine  = -1
         my_coarsen = 0
        endif 
     endif

     if(it==0)then                                   ! Limit initial resolution 
     if( .not. (any(rad(ix^S) <= (eqpar(Rwind_)*2.5d0))) ) then
        if( level > 1 ) then
	 my_refine  = -1
         my_coarsen = 1
        else if( level==1 ) then
 	 my_refine  = -1
         my_coarsen = 0
        endif
     endif
     endif

return
end subroutine specialrefine_grid
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
double precision :: H(ixI^S), alpha_T(ixI^S), Tg(ixI^S), Td(ixI^S)
double precision :: ndens(ixI^S)

double precision :: Pg(ixI^S)

return

end subroutine specialvar_output

!=============================================================================
subroutine specialvarnames_output

! newly added variables need to be concatenated with the varnames/primnames string

use mod_global_parameters
!-----------------------------------------------------------------------------
primnames= TRIM(primnames)//' '//'Td'
wnames=TRIM(wnames)//' '//'Td'

return
end subroutine specialvarnames_output
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

if (iflag ==nw+1) then
 {^IFTWOC   var(ixI^S) = dsqrt(w(ixI^S,m1_)**2 + w(ixI^S,m2_)**2)/w(ixI^S,rho_)}
 {^IFTHREEC var(ixI^S) = dsqrt(w(ixI^S,m1_)**2 + w(ixI^S,m2_)**2+ w(ixI^S,m3_)**2)/w(ixI^S,rho_)}
endif

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
!!wB0(ixO^S,1:ndir)=wB0(ixO^S,1:ndir)+user defined steady potential field

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
!call conserve(ixG^L,ixO^L,w,patchw)

end subroutine bc_int
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
