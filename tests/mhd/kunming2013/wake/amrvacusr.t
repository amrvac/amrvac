!=============================================================================
! amrvacusr.t.wakemhd

! INCLUDE:amrvacnul/specialini.t
! INCLUDE:amrvacnul/speciallog.t
INCLUDE:amrvacnul/specialbound.t
INCLUDE:amrvacnul/specialsource.t
INCLUDE:amrvacnul/specialimpl.t
INCLUDE:amrvacnul/usrflags.t
INCLUDE:amrvacnul/correctaux_usr.t
!=============================================================================
subroutine initglobaldata_usr

include 'amrvacdef.f'
!-----------------------------------------------------------------------------
eqpar(gamma_)=5.0d0/3.0d0
eqpar(eta_)=0.0001d0

end subroutine initglobaldata_usr
!=============================================================================
subroutine initonegrid_usr(ixG^L,ix^L,w,x)

! initialize one grid

include 'amrvacdef.f'

integer, intent(in) :: ixG^L, ix^L
double precision, intent(in) :: x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)

double precision:: dA,dM,width,dv,sigma,alfa,Lx,rho0,p0

logical, save :: first=.true.
logical:: patchw(ixG^T)
!----------------------------------------------------------------------------
oktest = index(teststr,'initonegrid_usr')>=1
if (oktest) write(unitterm,*) ' === initonegrid_usr  (in ) : ', &
      'ixG^L : ',ixG^L

   dA=0.2d0
   dM=3.0d0
   width=one
   dv=0.01d0
   sigma=one
   alfa=0.35d0
   Lx=two*dpi/alfa
   rho0=one
   w(ix^S,rho_)=rho0
   if(first .and. mype==0)then
      write(*,*)'Doing 2.5D resistive MHD, Wake problem'
      write(*,*)'A  - M - gamma:',dA,dM,eqpar(gamma_)
      write(*,*)'dv sigma:',dv,sigma
      write(*,*)'alfa Lx:',alfa,Lx
      write(*,*)'Assuming eta set, using value:',eqpar(eta_)
      first=.false.
   endif
   w(ix^S,m1_)=rho0*(one-one/dcosh(x(ix^S,2)))
   w(ix^S,m2_)=rho0*(dv*(dsin(two*dpi*x(ix^S,1)/Lx)*dtanh(x(ix^S,2))&
               +dcos(two*dpi*x(ix^S,1)/Lx)/dcosh(x(ix^S,2)))&
               *dexp(-(x(ix^S,2)/sigma)**2))
   w(ix^S,b1_)=dA*dtanh(x(ix^S,2)/width)
   w(ix^S,b2_)=zero
   if (ndir/=3) call mpistop("this is a mhd23 problem!")
   w(ix^S,m^NC_)=zero
   w(ix^S,b^NC_)=dA/dcosh(x(ix^S,2)/width)
   p0=one/(dM**2)/eqpar(gamma_)
   w(ix^S,e_)=p0/(eqpar(gamma_)-1)+&
     half*((^C&w(ix^S,m^C_)**2+)/rho0+(^C&w(ix^S,b^C_)**2+))

return
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
! .. local ..
double precision :: current(ixG^T,7-2*ndir:3)
integer          :: idirmin
double precision :: wloc(ixI^S,1:nw)
double precision :: pth(ixG^T), rho(ixG^T)
!-----------------------------------------------------------------------------

wloc(ixI^S,1:nw)=w(ixI^S,1:nw)
if(saveprim)then
  pth(ixO^S)=wloc(ixO^S,p_)
else
  call getpthermal(wloc,x,ixI^L,ixO^L,pth)
endif
rho(ixO^S)=wloc(ixO^S,rho_)
w(ixO^S,nw+1)=pth(ixO^S)/rho(ixO^S)

call getcurrent(wloc,ixI^L,ixO^L,idirmin,current)
w(ixO^S,nw+2)=current(ixO^S,3)

end subroutine specialvar_output
!=============================================================================
subroutine specialvarnames_output

! newly added variables need to be concatenated with the varnames/primnames string

include 'amrvacdef.f'
!-----------------------------------------------------------------------------

primnames= TRIM(primnames)//' '//'T'
primnames= TRIM(primnames)//' '//'jz'
wnames=TRIM(wnames)//' '//'T'
wnames=TRIM(wnames)//' '//'jz'


end subroutine specialvarnames_output
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
!=============================================================================
! amrvacusr.t.wakemhd
!=============================================================================
