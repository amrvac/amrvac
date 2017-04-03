!=============================================================================
! amrvacusr.t.rimhd22

! 
!





! -d=22 -phi=0 -z=0 -g=10,10 -p=hd -eos=default -nf=0 -ndust=4 -u=nul -arch=default

!=============================================================================
subroutine initglobaldata_usr

use mod_global_parameters
double precision:: r(0:^NDS)
integer:: i
logical, save :: first=.true.
!-----------------------------------------------------------------------------
eqpar(gamma_)=1.4d0
 
length_convert_factor     = 0.001d0*3.08567758d18  ! 0.1 pc          ! normalization for distance
w_convert_factor(rho_)  = 1.0d-20              ! normalization for rho
w_convert_factor(mom(1))   = 1.0d7                ! normalization for speed

w_convert_factor(mom(2))   = w_convert_factor(mom(1))
{^IFTHREED
w_convert_factor(mom(3))   = w_convert_factor(mom(1))
}
time_convert_factor          = length_convert_factor/w_convert_factor(mom(1))
w_convert_factor(p_)    = w_convert_factor(rho_)*(w_convert_factor(mom(1))**2)




{#IFDEF DUST
eqpar(mu_)=2.3d0 ! moleculair waterstof

rhodust(1:^NDS) = 3.3d0    ! dust grain density
eqpar(min_ar_)  = 1.0d-7   ! minimum dust grain size (cm)
eqpar(max_ar_)  = 500.0d-7 ! maximum dust grain size (cm)
{w_convert_factor(rhod^DS_)   = w_convert_factor(rho_)\}
{^DS&{^C&w_convert_factor(v^Cd^DS_) = w_convert_factor(v^C_);}\}

! === rescale dust quantities to dimensionless scale === !
rhodust(1:^NDS)  = rhodust(1:^NDS)/w_convert_factor(rhod1_)
eqpar(min_ar_)= eqpar(min_ar_)/length_convert_factor
eqpar(max_ar_)= eqpar(max_ar_)/length_convert_factor 


!-------------------------------

! here the dust sizes are defined. Note that this is
! now done differently for the method used in (van Marle et al. (2011)).

! first dust sizes in Ndust bins, with all bins having equal total mass.
! To do this, assume the particle distribution goes as r^-3.5

r(0) = eqpar(min_ar_)
do i=1,^NDS
    r(i) = (dsqrt(r(i-1)) +(dsqrt(eqpar(max_ar_))- &
        dsqrt(eqpar(min_ar_)))/^NDS)**2.0d0
    dsdust(i) = r(i)-r(i-1)
end do
! now calculate the weigthed mean size of each bin, again assuming n goes as r^-3.5
do i=1,^NDS
    sdust(i) = (5.0d0/3.0d0)*(r(i)**(-1.5d0) - r(i-1)**(-1.5d0)) &
        /(r(i)**(-2.5d0) - r(i-1)**(-2.5d0))
end do

if(first)then
  if(mype==0)then
    do i=1,^NDS
        write(*,*) 'Dust type ',i,': grain radius r=',sdust(i)*length_convert_factor
    end do
  endif
  first=.false.
endif
}



end subroutine initglobaldata_usr
!=============================================================================
subroutine initonegrid_usr(ixG^L,ix^L,w,x)

! initialize one grid 

use mod_global_parameters

integer, intent(in) :: ixG^L, ix^L
double precision, intent(in) :: x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)

double precision :: xshock, xbound, rhopost, vpost, ppost, alpha
double precision :: Prat,alfa,c_pre
double precision :: M,eta
logical, save :: first=.true.
!----------------------------------------------------------------------------

{^IFONED call mpistop("This is a 2D HD problem: Richtmyer Meshkov")}
{^IFTHREED call mpistop("This is a 2D HD problem: Richtmyer Meshkov")}
{^IFTWOD

! === Shock and CD position === !
xshock=1.4d0
xbound=1.5d0 ! xbound must be bigger than xshock!

M = 2.0d0
eta = 3.0d0
alpha = (3.14159265d0/4.0d0)

! compute the RH related states
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
  print *,' post-shock density: ',rhopost
  print *,' post-shock velocity:',vpost
  print *,' post-shock pressure:',ppost
  print *,'============================================================='
  first=.false.
endif


where(x(ix^S,1)>xshock.and.(x(ix^S,1)>x(ix^S,2)/dtan(alpha)+xbound))
   ! pre shock region
   w(ix^S,rho_)=eqpar(gamma_)*eta
   w(ix^S,mom(1))=zero
   w(ix^S,mom(2))=zero
   w(ix^S,e_)=one/(eqpar(gamma_)-one)
{#IFDEF DUST
   {^DS& w(ix^S,rhod^DS_)=(0.01d0*eqpar(gamma_)*eta)/^NDS\}
}
endwhere
where(x(ix^S,1)>xshock.and.(x(ix^S,1)<=x(ix^S,2)/dtan(alpha)+xbound))
   ! pre shock region
   w(ix^S,rho_)=eqpar(gamma_)
   w(ix^S,mom(1))=zero
   w(ix^S,mom(2))=zero
   w(ix^S,e_)=one/(eqpar(gamma_)-one)
{#IFDEF DUST
   {^DS& w(ix^S,rhod^DS_)=zero\}
}
endwhere
where(x(ix^S,1)<=xshock)
   ! post shock region
   w(ix^S,rho_)= rhopost
   w(ix^S,mom(1)) = rhopost*vpost
   w(ix^S,mom(2)) = zero
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
!w(ixO^S,nw+2)=dlog10(w(ixO^S,rho_))

!!call mpistop("special output file undefined")


end subroutine specialvar_output
!=============================================================================
subroutine specialvarnames_output

! newly added variables to be concatenated with the primnames/w_names string

use mod_global_parameters
!-----------------------------------------------------------------------------
! Example : as above in specialvar_output
primnames= TRIM(primnames)//' '//'schlierrho'
w_names=TRIM(w_names)//' '//'schlierrho'

end subroutine specialvarnames_output
!=============================================================================
! amrvacusr.t.rimhd22
!=============================================================================
