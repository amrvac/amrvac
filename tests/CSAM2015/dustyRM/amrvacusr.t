!=============================================================================
! amrvacusr.t.dustyrimhd
! Richtmyer-Meshkov hydro with dust (2D)
!  RM through planar shock impinging on inclined density discontinuity
! setting for 2D hydro with 4 dust species:
! -d=22 -phi=0 -z=0 -g=14,14 -p=hd -eos=default -nf=0 -ndust=4 -u=nul -arch=default

INCLUDE:amrvacnul/specialbound.t
INCLUDE:amrvacnul/specialsource.t
INCLUDE:amrvacnul/specialimpl.t
INCLUDE:amrvacnul/usrflags.t
!=============================================================================
subroutine initglobaldata_usr

include 'amrvacdef.f'

double precision:: r(0:^NDS)
integer:: i
logical, save :: first=.true.
!-----------------------------------------------------------------------------
eqpar(gamma_)=1.4d0
 
normvar(0)     = 0.001d0*3.08567758d18  ! 0.001 pc distance (in cm)
normvar(rho_)  = 1.0d-20                ! normalization for rho g/cc
normvar(v1_)   = 1.0d7                  ! normalization for speed cm/s
normvar(v2_)   = normvar(v1_)
normt          = normvar(0)/normvar(v1_)
normvar(p_)    = normvar(rho_)*(normvar(v1_)**2)


{#IFDEF DUST
eqpar(mu_)=2.3d0 ! molecular hydrogen

rhodust(1:^NDS) = 3.3d0    ! dust grain density
eqpar(min_ar_)  = 1.0d-7   ! minimum dust grain size (cm)
eqpar(max_ar_)  = 500.0d-7 ! maximum dust grain size (cm)
{normvar(rhod^DS_)   = normvar(rho_)\}
{^DS&{^C&normvar(v^Cd^DS_) = normvar(v^C_);}\}

! === rescale dust quantities to dimensionless scale === !
rhodust(1:^NDS) = rhodust(1:^NDS)/normvar(rhod1_)
eqpar(min_ar_)  = eqpar(min_ar_)/normvar(0)
eqpar(max_ar_)  = eqpar(max_ar_)/normvar(0) 
!-------------------------------

! here the dust sizes are defined. Ndust bins, with all bins having equal total mass.
! To do this, assume the particle distribution goes as r^-3.5

r(0) = eqpar(min_ar_)
do i=1,^NDS
    r(i) = (dsqrt(r(i-1)) +(dsqrt(eqpar(max_ar_))- &
        dsqrt(eqpar(min_ar_)))/^NDS)**2.0d0
    dsdust(i) = r(i)-r(i-1)
end do
! now calculate the weighted mean size of each bin, again assuming n goes as r^-3.5
do i=1,^NDS
    sdust(i) = (5.0d0/3.0d0)*(r(i)**(-1.5d0) - r(i-1)**(-1.5d0)) &
        /(r(i)**(-2.5d0) - r(i-1)**(-2.5d0))
end do

if(first)then
  if(mype==0)then
    do i=1,^NDS
        write(*,*) 'Dust type ',i,': grain radius r=',sdust(i)*normvar(0)
    end do
  endif
  first=.false.
endif
}

end subroutine initglobaldata_usr
!=============================================================================
subroutine initonegrid_usr(ixG^L,ix^L,w,x)

! initialize one grid 

include 'amrvacdef.f'

integer, intent(in) :: ixG^L, ix^L
double precision, intent(in) :: x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)

double precision :: xshock, xbound, rhopost, vpost, ppost, alpha
double precision :: Prat,alfa,c_pre
double precision :: M,eta
logical, save :: first=.true.
!----------------------------------------------------------------------------

{^IFONED   call mpistop("This is a 2D HD problem: Richtmyer Meshkov")}
{^IFTHREED call mpistop("This is a 2D HD problem: Richtmyer Meshkov")}
{^IFTWOD

! === Shock and CD position === !
xshock=0.4d0
xbound=0.5d0 ! xbound must be bigger than xshock!
! inclination angle for CD
alpha = (3.14159265d0/4.0d0)

! Mach number of planar shock
M = 2.0d0
! inclination angle for CD
alpha = (3.14159265d0/4.0d0)
! density ratio across CD
eta=3.0d0

! compute the RH related states across the planar shock
Prat=one/(one+(M**2-one)*two*eqpar(gamma_)/(eqpar(gamma_)+one))
alfa=(eqpar(gamma_)+one)/(eqpar(gamma_)-one)
c_pre = one ! pre shock sound speed
rhopost=eqpar(gamma_)*(alfa+Prat)/(alfa*Prat+one)
ppost=one/Prat
vpost=c_pre*M*(one-(alfa*Prat+one)/(alfa+Prat))

if (mype==0.and.first) then
  print *,'============================================================='
  print *,' HD Richtmyer Meshkov simulation '
  print *,'============================================================='
  print *,' Mach number of shock: ',M
  print *,' post-shock density: ',rhopost
  print *,' post-shock velocity:',vpost
  print *,' post-shock pressure:',ppost
  print *,' Density ratio at CD: ',eta
  print *,'============================================================='
  first=.false.
endif


where(x(ix^S,1)>xshock.and.(x(ix^S,1)>x(ix^S,2)/dtan(alpha)+xbound))
   ! pre shock region
   w(ix^S,rho_)=eqpar(gamma_)*eta
   w(ix^S,m1_)=zero
   w(ix^S,m2_)=zero
   w(ix^S,e_)=one/(eqpar(gamma_)-one)
{#IFDEF DUST
   {^DS& w(ix^S,rhod^DS_)=(0.01d0*eqpar(gamma_)*eta)/^NDS\}
}
endwhere
where(x(ix^S,1)>xshock.and.(x(ix^S,1)<=x(ix^S,2)/dtan(alpha)+xbound))
   ! pre shock region
   w(ix^S,rho_)=eqpar(gamma_)
   w(ix^S,m1_)=zero
   w(ix^S,m2_)=zero
   w(ix^S,e_)=one/(eqpar(gamma_)-one)
{#IFDEF DUST
   {^DS& w(ix^S,rhod^DS_)=zero\}
}
endwhere
where(x(ix^S,1)<=xshock)
   ! post shock region
   w(ix^S,rho_)= rhopost
   w(ix^S,m1_) = rhopost*vpost
   w(ix^S,m2_) = zero
   w(ix^S,e_)  = ppost/(eqpar(gamma_)-one)+0.5d0*rhopost*vpost**2
{#IFDEF DUST
   {^DS& w(ix^S,rhod^DS_)=zero\}
}
endwhere
}

{#IFDEF DUST
{^DS& w(ix^S,v1d^DS_)=zero \}
{^DS& w(ix^S,v2d^DS_)=zero \}
}

end subroutine initonegrid_usr
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

double precision                   :: gradrho(ixG^T),rho(ixG^T),drho(ixG^T)
double precision                   :: kk,kk0,grhomax,kk1
integer                            :: idims
!-----------------------------------------------------------------------------

! Example: assuming nwauxio=1 at convert stage 

rho(ixI^S)=w(ixI^S,rho_)
gradrho(ixO^S)=zero
do idims=1,ndim
   select case(typegrad)
   case("central")
     call gradient(rho,ixI^L,ixO^L,idims,drho)
   case("limited")
     call gradientS(rho,ixI^L,ixO^L,idims,drho)
   end select
   gradrho(ixO^S)=gradrho(ixO^S)+drho(ixO^S)**2.0d0
enddo
gradrho(ixO^S)=dsqrt(gradrho(ixO^S))
kk=5.0d0
kk0=0.001d0
kk1=1.0d0
grhomax=2000.0d0

! putting the schlierplot of density in nwauxio=1
w(ixO^S,nw+1)=dexp(-kk*(gradrho(ixO^S)-kk0*grhomax)/(kk1*grhomax-kk0*grhomax))

end subroutine specialvar_output
!=============================================================================
subroutine specialvarnames_output

! newly added variables to be concatenated with the primnames/wnames string

include 'amrvacdef.f'
!-----------------------------------------------------------------------------
oktest = index(teststr,'printlog')>=1

!!call mpistop("special wnames and primnames undefined")

! Example : as above in specialvar_output
primnames= TRIM(primnames)//' '//'schlierrho'
wnames=TRIM(wnames)//' '//'schlierrho'

end subroutine specialvarnames_output
!=============================================================================
! amrvacusr.t.dustyrimhd
!=============================================================================
