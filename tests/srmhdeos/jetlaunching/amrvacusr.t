!=============================================================================
! amrvacusr.t.srmhdeos23JetLaun
! configuration
! setamrvac -d=23 -phi=3 -z=2 -g=16,16 -p=srmhdeos -u=srmhdeos23JetLaun
! last modified July 17, 2009
! Hendrik van Eerten / Zakaria Meliani

! if not commented the following include statements provide dummy routines
!INCLUDE:amrvacnul.specialini.t
!INCLUDE:amrvacnul.specialsource.t
!INCLUDE:amrvacnul.speciallog.t
!INCLUDE:amrvacnul.specialbound.t

!===============================================================================

subroutine initglobaldata_usr

  use mod_global_parameters
!----------------------------------------------------------------------------
{^IFGLM
eqpar(Ch_)= .5d0
eqpar(Cr_)= 0.18d0
}

! polytropic index
eqpar(gamma_)=5.0d0/3.0d0

! inlet parameters
eqpar(rhoin_)     = 1.0d0
eqpar(beta_)      = 1.0d-2
eqpar(vphim_)     = 5.0d-1
eqpar(Bmax_)      = 1.0d0

eqpar(rhoe_)      = 1.0d-1

! radii
eqpar(rs_)        = 10.0d0
eqpar(r1_)        = 1.0d0
eqpar(r2_)        = 1.0d6

end subroutine initglobaldata_usr
!===============================================================================

subroutine initonegrid_usr(ixG^L,ixO^L,w,x)
  ! initialize one grid 

  use mod_global_parameters ! standard definitions, always needed

  ! .. scalars ..
  integer         :: ixG^L,ixO^L
  ! .. arrays ..
  double precision, intent(in) :: x(ixG^S,1:ndim)
  double precision, intent(inout) :: w(ixG^S,1:nw)

  ! .. Local ..
  double precision :: d0,p0,b0
  ! .. Local Logical ..
  logical, save :: first=.true. ! used to write settings to file only once
  logical :: patchw(ixG^T)
  !----------------------------------------------------------------------------
  
  patchw(ixG^S)=.false.
  ! initialize the speed
  ^C&w(ixO^S,v^C_)=zero;
  ! density and presssure
      d0    = eqpar(rhoe_)
      b0    = (eqpar(r1_)**2.0d0+eqpar(rs_)**2.0d0)*eqpar(Bmax_)
      p0    = half*eqpar(beta_)*eqpar(Bmax_)**2.0d0
      w(ixO^S,rho_) = d0
      w(ixO^S,pp_)  = p0

  w(ixO^S,br_) = b0*x(ixO^S,r_)*(x(ixO^S,r_)**2.0d0+(x(ixO^S,z_)+eqpar(rs_))**2.0d0)**(-1.5d0)
  w(ixO^S,bz_) = b0*(x(ixO^S,z_)+eqpar(rs_))*(x(ixO^S,r_)**2.0d0+(x(ixO^S,z_)+eqpar(rs_))**2.0d0)**(-1.5d0)

  w(ixO^S,bphi_) = zero 
  w(ixO^S,lfac_)=one/dsqrt(one-({^C&w(ixO^S,v^C_)**2.0d0+}))
  {^IFGLM w(ixO^S,psi_)=zero}

  if(first.and.mype==0) call write_setting ! output settings

  ! translate the primitive variables to conserved variables
  ! only those entries that are 'false'. Since we use primitive
  ! variables everywhere, we set all elements of this boolean array to .false.
  ! amrvac offers to possibility to use the spatial components of the four
  ! velocity instead of the proper velocity. This is the case if
  ! useprimitiveRel is .true. and in that case those entries that are 
  ! initialised using the w entries need to be modified accordingly.
  if (useprimitiveRel) then
    where(.not.patchw(ixO^S))
      {^C&w(ixO^S,u^C_)=w(ixO^S,lfac_)*w(ixO^S,v^C_);\}
    end where
  endif
  call conserve(ixG^L,ixO^L,w,patchw)

  return

end subroutine initonegrid_usr
!=============================================================================
!=============================================================================
subroutine specialbound_usr(qt,ixI^L,ixO^L,iw,iB,w,x)

! special boundary types, user defined

use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L, iw, iB
double precision, intent(in) :: qt, x(ixI^S,1:ndim)
double precision, intent(inout) :: w(ixI^S,1:nw)

! .. local ..
 double precision :: pin, b0
 double precision, dimension(ixG^T)     :: Omega, rd, v1, v2, v3, VdotB, sqrU, lor,rhoh,&
      E_Th,E,Pa,vp,csound2,pprof
 integer          :: ix^D,ixL^D,ixR^D,ixIn^L,ixC^L,ixl^L,ixR^L
 logical          :: patchw(ixG^T)
!----------------------------------------------------------------------------

   patchw(ixG^T)=.false.


      b0    = (eqpar(r1_)**2.0d0+eqpar(rs_)**2.0d0)*eqpar(Bmax_)
  !  presssure
      pin       = half*eqpar(beta_)*eqpar(Bmax_)**2.0d0
      pprof(ixO^S)=pin

      rd(ixO^S) = x(ixO^S,r_)

! All for getting the injection speed:
!============================================================
      E_Th(ixO^S) = pprof(ixO^S)/(eqpar(gamma_)-one) 
      E(ixO^S) = E_Th(ixO^S) + dsqrt(E_Th(ixO^S)**2.0d0+eqpar(rhoin_)**2.0d0)
      rhoh(ixO^S) = half*((eqpar(gamma_)+one) * E(ixO^S)-&
               (eqpar(gamma_)-one)* eqpar(rhoin_)*(eqpar(rhoin_)/E(ixO^S)))


      E(ixO^S)=(rhoh(ixO^S)&
           +dsqrt(rhoh(ixO^S)**2.0d0+(eqpar(gamma_)**2.0d0&
           -one)*eqpar(rhoin_)**2.0d0))/(eqpar(gamma_)&
           +one)
      Pa(ixO^S) = (eqpar(gamma_)-one)&
           /2.0d0* (E(ixO^S)-eqpar(rhoin_)**2.0d0/E(ixO^S))
      csound2(ixO^S)=(Pa(ixO^S)*((eqpar(gamma_)+one)+(eqpar(gamma_)-one)&
           *eqpar(rhoin_)/E(ixO^S))**2.0d0)/(2.0*rhoh(ixO^S))
!============================================================

      where(dabs(rd(ixO^S))<eqpar(r2_) .and. dabs(rd(ixO^S))>=eqpar(r1_))
         Omega(ixO^S) = eqpar(vphim_)/eqpar(r1_)*(dabs(rd(ixO^S))/eqpar(r1_))**(-1.5)
         vp(ixO^S) = dsqrt(csound2(ixO^S))
      end where
      where(dabs(rd(ixO^S))<eqpar(r1_))
         Omega(ixO^S) = eqpar(vphim_)/eqpar(r1_)
!         Omega(ixO^S) = 0.0d0
!         vp(ixO^S)    = 0.0d0
!         pprof(ixO^S) = 1.0d2*pin*(one-dabs(rd(ixO^S))/eqpar(r1_))+pin
      end where


!      where(dabs(rd(ixO^S))<eqpar(r2_))

!         w(ixO^S,b1_) =  b0*x(ixO^S,r_)*(x(ixO^S,r_)**2.0d0+(x(ixO^S,z_)+eqpar(rs_))**2.0d0)**(-1.5d0)
         w(ixO^S,b2_) = b0*(x(ixO^S,z_)+eqpar(rs_))*(x(ixO^S,r_)**2.0d0+(x(ixO^S,z_)+eqpar(rs_))**2.0d0)**(-1.5d0)


! vp||BP:
         v2(ixO^S) = vp(ixO^S)*((one+(w(ixO^S,b1_)/w(ixO^S,b2_))**2.0d0))**(-0.5d0)
         v1(ixO^S) = (vp(ixO^S)**2.0d0-v2(ixO^S)**2.0d0)**0.5d0
         v1(ixO^S) = v1(ixO^S)*sign(one,w(ixO^S,b1_))

         v3(ixO^S)  = rd(ixO^S)*Omega(ixO^S)
         v3(ixO^S)  =v3(ixO^S)+w(ixO^S,b3_)*vp(ixO^S)/(w(ixO^S,b1_)**2.0d0+w(ixO^S,b2_)**2.0d0)**0.5d0

         lor(ixO^S)=(one-({v^C(ixO^S)**2.0d0 +}))**(-0.5d0)


         w(ixO^S,rho_)=eqpar(rhoin_)
         w(ixO^S,pp_)=pprof(ixO^S)
         w(ixO^S,v1_)=v1(ixO^S)
         w(ixO^S,v2_)=v2(ixO^S)
         w(ixO^S,v3_)=v3(ixO^S)
         w(ixO^S,lfac_)=lor(ixO^S)


!end where

  if (useprimitiveRel) then
      {^C&w(ixO^S,u^C_)=w(ixO^S,lfac_)*w(ixO^S,v^C_);\}
  endif

!print*,'r:'
!print*,x(ixO^S,r_)
!print*,'rho:'
!print*,w(ixO^S,rho_)
!print*,'u1:'
!print*,w(ixO^S,u1_)
!print*,'u2:'
!print*,w(ixO^S,u2_)
!print*,'u3:'
!print*,w(ixO^S,u3_)
!print*,'pp:'
!print*,w(ixO^S,pp_)
!print*,'b1:'
!print*,w(ixO^S,b1_)
!print*,'b2:'
!print*,w(ixO^S,b2_)
!print*,'b3:'
!print*,w(ixO^S,b3_)
!print*,'lfac:'
!print*,w(ixO^S,lfac_)


{^IFGLM w(ixO^S,psi_)=zero}
call conserve(ixI^L,ixO^L,w,patchw)

!print*,'d:'
!print*,w(ixO^S,d_)

end subroutine specialbound_usr
!=============================================================================

!=============================================================================

subroutine specialsource(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

! w[iw]=w[iw]+qdt*S[wCT,qtC,x] where S is the source based on wCT within ixO
! use ieee_arithmetic ! DEBUG FEATURE

use mod_global_parameters

integer, intent(in)                :: ixI^L, ixO^L, iw^LIM
double precision, intent(in)       :: qdt, qtC, qt, x(ixI^S,1:ndim)
double precision, intent(inout)    :: wCT(ixI^S,1:nw), w(ixI^S,1:nw)

! .. Local..
double precision                   :: xg^D

double precision,dimension(ixG^T)        :: Radius{^IFREL ,Rho}
double precision,dimension(ixG^T,1:ndim) :: dRdx,Acc

!-----------------------------------------------------------------------------
  ! set the auxiliary variables

end subroutine specialsource

!=============================================================================

subroutine getdt_special(w,ixI^L,ixO^L,dtnew,dx^D,x)

  ! Limit "dt" further if necessary, e.g. due to the special source terms.
  ! The getdt_courant (CFL condition) and the getdt subroutine in the AMRVACPHYS
  ! module have already been called.

  use mod_global_parameters

  integer, intent(in) :: ixI^L, ixO^L
  double precision, intent(in) :: dx^D, x(ixI^S,1:ndim)
  double precision, intent(inout) :: w(ixI^S,1:nw), dtnew

  ! .. local variables ..
  double precision,dimension(ixG^T)        :: Radius
  double precision,dimension(ixG^T,1:ndim) :: dRdx
  double precision                         :: xg^D,dxinv(1:ndim), dtgrav
  integer                                  :: idims
  !-----------------------------------------------------------------------------


end subroutine getdt_special

!=============================================================================

subroutine specialeta(w,ixI^L,ix^L,idirmin,x,current,eta)
  ! DUMMY SUBROUTINE

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

  ! special refinement is necessary because the entire BM shock width tends
  ! to be initially smaller than a single grid cell. If we don't force extra
  ! refinement then the shock simply won't be seen
  use mod_global_parameters

  integer, intent(in) :: igrid, level, ixG^L, ix^L
  double precision, intent(in) :: qt, w(ixG^S,1:nw), x(ixG^S,1:ndim)
  integer, intent(inout) :: my_refine, my_coarsen

  !--------------------------------------------------------------

  if (minval(x(ix^S,2)) <= one) my_refine=1


end subroutine specialrefine_grid
!=============================================================================
!=============================================================================
subroutine specialvarforerrest(ixI^L,ixO^L,iflag,w,var)
! DUMMY SUBROUTINE
  use mod_global_parameters

  integer, intent(in)          :: ixI^L,ixO^L,iflag
  double precision, intent(in) :: w(ixI^S,1:nw)
  double precision, intent(out):: var(ixG^T)
  
  !-----------------------------------------------------------------------------
  
  if (iflag >nw) then
    call mpistop(' iflag> nw, made change in parfile or in user file')
  end if
  var(ixI^S) = zero 

end subroutine specialvarforerrest
!=============================================================================
subroutine process_allgrid_usr(iit,qt)
use mod_global_parameters

! .. scalars ..
integer,intent(in)          :: iit
double precision, intent(in):: qt

! .. local variables ..
integer          :: igrid,iigrid,indice_vec(^ND,2),save_mype
double precision :: save_maxlfac_recv
!-----------------------------------------------------------------------------

end subroutine process_allgrid_usr
!=============================================================================

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
!-----------------------------------------------------------------------------

call mpistop("special output file undefined")

! Example: assuming nwauxio=1 at convert stage and desire to see -w(1)
! w(ixO^S,nw+1)=-w(ixO^S,1)

end subroutine specialvar_output
!=============================================================================
subroutine specialvarnames_output(wnameio)

! newly added variables need to be concatenated with the varnames/primnames string

use mod_global_parameters
character(len=10) ::  wnameio(1:nw+nwauxio)
!-----------------------------------------------------------------------------

call mpistop("special varnames and primnames undefined")

! Example : as above in specialvar_output, assuming relativistic HD here...
! wnameio(nw+1)='rho'
end subroutine specialvarnames_output
!=============================================================================
subroutine Global_useroutput
use mod_global_parameters
  ! .. local scalars ..
  integer                    :: i,unituser
  ! .. local character
  character(len=80) :: filenameusr
  ! .. local array ..
  logical :: fileopen
  logical    first
  data  first/.true./

!-------------------------------------------------------------------
end subroutine Global_useroutput
!=============================================================================

Subroutine write_setting
  use mod_global_parameters

  ! .. local scalars ..
  integer                    :: i,unituser
  ! .. local character
  character(len=80) :: filenameusr
  ! .. local array ..
  logical :: fileopen
  logical    first
  data  first/.true./
  
  !-------------------------------------------------------------------------
  
  if  (.not.convert .and. it==0 .and. first .and. mype==0) then

    unituser=99
    write(filenameusr,"(a,a,a)") TRIM(filenameout),'setting',".init"
    inquire(unituser,opened=fileopen)
    if(.not.fileopen) open(unituser,file=filenameusr,status='unknown')
    write(unituser,*) '****** Description of usr file'
    write(unituser,*) ' relativistic Jet Launching'
    write(unituser,*) 'polytropic index', eqpar(gamma_)
    write(unituser,*) 'transition inner-outer jet', eqpar(r1_)
    write(unituser,*) 'jetradius inner-outer jet', eqpar(r2_)
    
    write(unituser,*)'******************************************'
    close(unituser)
  end if

  first=.false. ! make sure we don't print this every time we generate a grid
  ! ('grids' are what is called 'blocks' in BLAST, whereas in BLAST we call
  ! the entire collection of 'blocks' the single 'grid'. Hence filenames like
  ! grid.cpp etc. Every block in blast / grid in amrvac contains a number of 
  ! cells )

end subroutine write_setting
 
!=============================================================================
! amrvacusr.t.srmhdeos23JetLaun
!=============================================================================
