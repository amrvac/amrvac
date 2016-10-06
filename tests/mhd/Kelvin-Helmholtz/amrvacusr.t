!=============================================================================
! amrvacusr.t.testmhd

! INCLUDE:amrvacnul/specialini.t
!INCLUDE:amrvacnul/speciallog.t
INCLUDE:amrvacnul/specialbound.t
INCLUDE:amrvacnul/specialsource.t
INCLUDE:amrvacnul/specialimpl.t
INCLUDE:amrvacnul/usrflags.t
!=============================================================================
subroutine initglobaldata_usr

use mod_global_parameters
!-----------------------------------------------------------------------------
eqpar(gamma_)=5.0d0/3.0d0
eqpar(eta_)=zero
{#IFDEF GLM
eqpar(Cr_)=-0.18d0
}
select case(iprob)
 case(3)
   eqpar(eta_)=0.00001d0
 case(18)
   eqpar(eta_)=0.0001d0
endselect 

end subroutine initglobaldata_usr
!=============================================================================
subroutine initonegrid_usr(ixG^L,ix^L,w,x)

! initialize one grid 

use mod_global_parameters

integer, intent(in) :: ixG^L,ix^L
double precision, intent(in) :: x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)

logical::          first
integer::          n1
double precision:: rholeft,rhoright,eleft,eright,m^Cleft,m^Cright,slocx^D
double precision:: v^Cleft,pleft,v^Cright,pright,gamm
double precision:: bx,byleft,bzleft,byright,bzright
double precision:: qv,width,dv,k1,sigma,rho0,p0,dA,dM,alfa,Lx,yjet
double precision:: rhodisk,rr0,rr1,x1disk,x2disk,v0
double precision:: cs,machs,macha,ff,Rjet
double precision, dimension(ixG^T) :: r,fslope,rinv 

integer,dimension(:),allocatable:: seed
integer::  ix2,seed_size
real::             ranx1d(ixGlo1:ixGhi1)
double precision:: ranx(ixG^T)

logical, dimension(ixG^T)           :: patchw(ixG^T)

data first/.true./
!----------------------------------------------------------------------------

if (typephys/='mhd') then
   call mpistop("test problems are all MHD problems: set typephys!")
end if

select case(iprob)
 case(1)
   ! setamrvac -d=22 -phi=0 -z=0 -g=16,16 -p=mhd -u=testmhd
   {^IFONED   call mpistop("prob 1 is 2D") }
   {^IFTHREED call mpistop("prob 1 is 2D") }
   qv=0.645d0
   width=0.05d0
   dv=0.01d0
   k1=6.d0*dpi
   sigma=0.20d0
   {^IFTWOD
   w(ixG^S,rho_)=one
   w(ixG^S,e_)=one
   w(ixG^S,b1_)=0.129d0
   w(ixG^S,b2_)=zero
   w(ixG^S,m1_)=qv*tanh((x(ixG^S,2)-(xprobmax2+xprobmin2)/two)/width)
   w(ixG^S,m2_)=dv*sin(k1*(x(ixG^S,1)-xprobmin1)/(xprobmax1-xprobmin1))*&
                exp(-((x(ixG^S,2)-(xprobmax2+xprobmin2)/two)/sigma)**2)
   w(ixG^S,e_)=w(ixG^S,e_)/(eqpar(gamma_)-1)+&
     half*(w(ixG^S,rho_)*(^C&w(ixG^S,m^C_)**2+)+(^C&w(ixG^S,b^C_)**2+))
   if(first .and. mype==0)then
      write(*,*)'Doing 2D MHD, Kelvin-Helmholtz problem, uniform density'
      write(*,*)'qv, width, dv, k1, sigma:'
      write(*,*)qv,width,dv,k1,sigma
      first=.false.
   endif
   \}
 case(2)
   ! setamrvac -d=22 -phi=0 -z=0 -g=16,16 -p=mhd -u=testmhd
   {^IFONED   call mpistop("prob 2 is 2D") }
   {^IFTHREED call mpistop("prob 2 is 2D") }
   qv=0.645d0
   width=0.05d0
   dv=0.01d0
   k1=6.d0*dpi
   sigma=0.20d0
   {^IFTWOD
   where(x(ixG^S,2)>=(xprobmax2+xprobmin2)/two)
      w(ixG^S,rho_)=0.5d0
   elsewhere
      w(ixG^S,rho_)=1.d0
   endwhere
   w(ixG^S,p_)=one
   w(ixG^S,b1_)=0.129d0
   w(ixG^S,b2_)=zero
   w(ixG^S,v1_)=qv*tanh((x(ixG^S,2)-(xprobmax2+xprobmin2)/two)/width)
   w(ixG^S,v2_)=dv*sin(k1*(x(ixG^S,1)-xprobmin1)/(xprobmax1-xprobmin1))*&
                exp(-((x(ixG^S,2)-(xprobmax2+xprobmin2)/two)/sigma)**2)
{#IFDEF TRACER
   where((x(ixG^S,1)-1.05)**2+(x(ixG^S,2)-2.05)**2<0.1**2)
     w(ixG^S,tr1_)=one
   elsewhere((x(ixG^S,1)-2.05)**2+(x(ixG^S,2)-2.05)**2<0.1**2)
     w(ixG^S,tr1_)=2.d0
   endwhere
}
   patchw(ixG^S)=.false.
   call conserve(ixG^L,ix^L,w,x,patchw)
   if(first .and. mype==0)then
      write(*,*)'Doing 2D MHD, Kelvin-Helmholtz problem, two density layers'
      write(*,*)'qv, width, dv, k1, sigma:'
      write(*,*)qv,width,dv,k1,sigma
      first=.false.
   endif
   \}
 case(3)
   ! setamrvac -d=22 -phi=0 -z=0 -g=16,16 -p=mhd -u=testmhd
   {^IFONED   call mpistop("prob 3 is 2D") }
   {^IFTHREED call mpistop("prob 3 is 2D") }
   {^IFTWOD
   w(ixG^S,rho_)=one
   w(ixG^S,e_)=one
   where(x(ixG^S,2)>=(xprobmax2-xprobmin2)/two)
      w(ixG^S,b1_)=0.129d0
   elsewhere
      w(ixG^S,b1_)=-0.129d0
   endwhere
   w(ixG^S,b2_)=zero
   qv=0.645d0
   width=0.05d0
   dv=0.01d0
   k1=two*dpi
   sigma=0.20d0
   if(first .and. mype==0)then
      write(*,*)'Doing 2D MHD, Reversed Kelvin-Helmholtz problem'
      write(*,*)'qv, width, dv, k1, sigma:'
      write(*,*)qv,width,dv,k1,sigma
      write(*,*)'Assuming eta set, using value:',eqpar(eta_)
      first=.false.
   endif
   w(ixG^S,m1_)=qv*tanh((x(ixG^S,2)-(xprobmax2+xprobmin2)/two)/width)
   w(ixG^S,m2_)=dv*sin(k1*(x(ixG^S,1)-xprobmin1)/(xprobmax1-xprobmin1))*&
                exp(-((x(ixG^S,2)-(xprobmax2+xprobmin2)/two)/sigma)**2)
   w(ixG^S,e_)=w(ixG^S,e_)/(eqpar(gamma_)-1)+&
     half*(w(ixG^S,rho_)*(^C&w(ixG^S,m^C_)**2+)+(^C&w(ixG^S,b^C_)**2+))
{#IFDEF TRACER
   where((x(ixG^S,1)-1.05)**2+(x(ixG^S,2)-2.05)**2<0.1**2)
     w(ixG^S,Dtr1_)=w(ixG^S,rho_)
   endwhere
}
   \}
case default
            write(unitterm,*)'Undefined Iprob in Userfile ',iprob
   Call mpistop(' --- initonegrid_usr ---')
end  select 

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

double precision :: divb(ixG^T)
!-----------------------------------------------------------------------------

! output divB1
call getdivb(w,ixI^L,ixO^L,divb)
w(ixO^S,nw+1)=divb(ixO^S)
! Example: assuming nwauxio=1 at convert stage and desire to see -w(1)
! w(ixO^S,nw+1)=-w(ixO^S,1)

end subroutine specialvar_output
!=============================================================================
subroutine specialvarnames_output

! newly added variables to be concatenated with the primnames/wnames string

use mod_global_parameters
!-----------------------------------------------------------------------------

! Example : as above in specialvar_output, assuming relativistic HD here...
primnames= TRIM(primnames)//' '//'divB'
wnames=TRIM(wnames)//' '//'divB'

end subroutine specialvarnames_output
!=============================================================================
subroutine userspecialconvert(qunitconvert)

! Allow user to use their own data-converting procedures

use mod_global_parameters
integer, intent(in) :: qunitconvert
character(len=20):: userconvert_type
!-----------------------------------------------------------------------------

end subroutine userspecialconvert
!=============================================================================
! amrvacusr.t.testmhd
!=============================================================================
