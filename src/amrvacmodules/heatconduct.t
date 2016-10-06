!##############################################################################
!> module heatconduct.t -- heat conduction for HD and MHD
!! 10.07.2011 developed by Chun Xia and Rony Keppens
!! 01.09.2012 moved to modules folder by Oliver Porth
!! 13.10.2013 optimized further by Chun Xia
!! 12.03.2014 implemented RKL2 super timestepping scheme to reduce iterations 
!! and improve stability and accuracy up to second order in time by Chun Xia.
!! 23.08.2014 implemented saturation and perpendicular TC by Chun Xia
!! 
!! PURPOSE: 
!! IN MHD ADD THE HEAT CONDUCTION SOURCE TO THE ENERGY EQUATION
!! S=DIV(KAPPA_i,j . GRAD_j T)
!! where KAPPA_i,j = kappa b_i b_j + kappe (I - b_i b_j)
!! b_i b_j = B_i B_j / B**2, I is the unit matrix, and i, j= 1, 2, 3 for 3D
!! IN HD ADD THE HEAT CONDUCTION SOURCE TO THE ENERGY EQUATION
!! S=DIV(kappa . GRAD T)
!! USAGE:
!! 1. in amrvacusr.t add: 
!!   a. INCLUDE:amrvacmodules/heatconduct.t
!!   b. in subroutine initglobaldata_usr, add conductivity:
!!      eqpar(kappa_)=kappa0*Teunit**3.5d0/Lunit/runit/vunit**3 
!!     kappa0=1.d-6 erg/cm/s/K**3.5 or 1.d-11 J/m/s/K**3.5  is Spitzer 
!!     conductivity for solar corona. Teunit, Lunit, runit, and vunit are the 
!!     unit of temperature, length, density, and velocity. Teunit and vunit are 
!!     dependent, e.g., vunit=sqrt(R Teunit).
!! 2. in amrvacusrpar.t add kappa_ 
!! 3. in definitions.h :
!!    #define TCRKL2
!! 4. in the methodlist of amrvac.par add:
!!    conduction=.true.
!! Saturation: (default off)
!!    in the methodlist of amrvac.par add:
!!    TCsaturate=.true.
!!    phi coefficient of saturated flux
!!    TCphi=1.d0
!! Addition thermal conduction perpendicular to magnetic field (Orlando 2008 ApJ
!! 678 247)
!!    in definition.h :
!!    #define TCPERPENDICULAR
!!    in subroutine initglobaldata_usr, add perpendicular conductivity:
!!eqpar(kappe_)=kappa1*nHunit**2/dsqrt(Teunit)/Bunit**2*Teunit/Lunit/runit/vunit**3
!!    kappa1=3.3d-16 erg Gauss**2 cm**5 /s/K**0.5 or 
!!           3.3d-41 J T**2 m**5 /s/K**0.5
!!    (nHunit and Bunit are the unit of number density and magnetic field.)
!!    in amrvacusrpar.t add kappe_ 
!=============================================================================
subroutine thermalconduct_RKL2(s,qdt,qt)
! Meyer 2012 MNRAS 422,2102
use mod_global_parameters

integer,intent(in) :: s
double precision, intent(in) :: qdt,qt

double precision :: omega1,cmu,cmut,cnu,cnut
double precision, allocatable :: bj(:)
integer:: iigrid, igrid,j
logical :: evenstep
!-----------------------------------------------------------------------------
if(e_<1) call mpistop("Thermal conduction requires by e_>0!")
bcphys=.false.

do iigrid=1,igridstail; igrid=igrids(iigrid);
   allocate(pw1(igrid)%w(ixG^T,1:nw))
   allocate(pw2(igrid)%w(ixG^T,1:nw))
   allocate(pw3(igrid)%w(ixG^T,1:nw))
   pw1(igrid)%w=pw(igrid)%w
   pw2(igrid)%w=pw(igrid)%w
   pw3(igrid)%w=pwold(igrid)%w
end do

allocate(bj(0:s))
bj(0)=1.d0/3.d0
bj(1)=bj(0)
if(s>1) then
  omega1=4.d0/dble(s**2+s-2)
  cmut=omega1/3.d0
else
  cmut=1.d0
endif

!$OMP PARALLEL DO PRIVATE(igrid)
do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
  if (.not.slab) mygeo => pgeo(igrid)
  if (B0field) then
     myB0_cell => pB0_cell(igrid)
     {^D&myB0_face^D => pB0_face^D(igrid)\}
  end if
  typelimiter=typelimiter1(node(plevel_,igrid))
  typegradlimiter=typegradlimiter1(node(plevel_,igrid))
  ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
  call evolve_step1(cmut,qdt,ixG^LL,ixM^LL,pw1(igrid)%w,pw(igrid)%w,&
                    px(igrid)%x,pw3(igrid)%w)
end do
!$OMP END PARALLEL DO
call getbc(qt,0d0,ixG^LL,pw1,pwCoarse,pgeo,pgeoCoarse,.false.,e_-1,1)
if(s==1) then
  do iigrid=1,igridstail; igrid=igrids(iigrid);
    pw(igrid)%w(ixG^T,e_)=pw1(igrid)%w(ixG^T,e_)
    deallocate(pw1(igrid)%w)
    deallocate(pw2(igrid)%w)
    deallocate(pw3(igrid)%w)
  end do
  bcphys=.true.
  return
endif
evenstep=.true.
do j=2,s
  bj(j)=dble(j**2+j-2)/dble(2*j*(j+1))
  cmu=dble(2*j-1)/dble(j)*bj(j)/bj(j-1)
  cmut=omega1*cmu
  cnu=dble(1-j)/dble(j)*bj(j)/bj(j-2)
  cnut=(bj(j-1)-1.d0)*cmut
  if(evenstep) then
!$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
      if(.not.slab) mygeo => pgeo(igrid)
      if(B0field) then
        myB0_cell => pB0_cell(igrid)
        {^D&myB0_face^D => pB0_face^D(igrid)\}
      end if
      typelimiter=typelimiter1(node(plevel_,igrid))
      typegradlimiter=typegradlimiter1(node(plevel_,igrid))
      ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
      call evolve_stepj(cmu,cmut,cnu,cnut,qdt,ixG^LL,ixM^LL,pw1(igrid)%w,&
                        pw2(igrid)%w,pw(igrid)%w,px(igrid)%x,pw3(igrid)%w)
    end do
!$OMP END PARALLEL DO
    call getbc(qt,0d0,ixG^LL,pw2,pwCoarse,pgeo,pgeoCoarse,.false.,e_-1,1)
    evenstep=.false.
  else
!$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
      if(.not.slab) mygeo => pgeo(igrid)
      if(B0field) then
         myB0_cell => pB0_cell(igrid)
        {^D&myB0_face^D => pB0_face^D(igrid)\}
      end if
      typelimiter=typelimiter1(node(plevel_,igrid))
      typegradlimiter=typegradlimiter1(node(plevel_,igrid))
      ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
      call evolve_stepj(cmu,cmut,cnu,cnut,qdt,ixG^LL,ixM^LL,pw2(igrid)%w,&
                        pw1(igrid)%w,pw(igrid)%w,px(igrid)%x,pw3(igrid)%w)
    end do
!$OMP END PARALLEL DO
    call getbc(qt,0d0,ixG^LL,pw1,pwCoarse,pgeo,pgeoCoarse,.false.,e_-1,1)
    evenstep=.true.
  end if 
end do
if(evenstep) then
  do iigrid=1,igridstail; igrid=igrids(iigrid);
    pw(igrid)%w(ixG^T,e_)=pw1(igrid)%w(ixG^T,e_)
    deallocate(pw1(igrid)%w)
    deallocate(pw2(igrid)%w)
    deallocate(pw3(igrid)%w)
  end do 
else
  do iigrid=1,igridstail; igrid=igrids(iigrid);
    pw(igrid)%w(ixG^T,e_)=pw2(igrid)%w(ixG^T,e_)
    deallocate(pw1(igrid)%w)
    deallocate(pw2(igrid)%w)
    deallocate(pw3(igrid)%w)
  end do 
end if
deallocate(bj)

bcphys=.true.
end subroutine thermalconduct_RKL2
!=============================================================================
subroutine evolve_stepj(qcmu,qcmut,qcnu,qcnut,qdt,ixI^L,ixO^L,w1,w2,w,x,wold)

use mod_global_parameters

integer, intent(in) :: ixI^L,ixO^L
double precision, intent(in) :: qcmu,qcmut,qcnu,qcnut,qdt
double precision, intent(in) :: w1(ixI^S,1:nw),w(ixI^S,1:nw),wold(ixI^S,1:nw)
double precision, intent(in) :: x(ixI^S,1:ndim)
double precision, intent(inout) :: w2(ixI^S,1:nw)

double precision :: tmp(ixI^S),tmp1(ixI^S),tmp2(ixI^S)
!-----------------------------------------------------------------------------
{^IFMHDPHYS
call heatconduct_mhd(tmp,tmp1,tmp2,ixI^L,ixO^L,w1,x)
}
{^IFHDPHYS
call heatconduct_hd(tmp,tmp1,tmp2,ixI^L,ixO^L,w1,x)
}
w2(ixO^S,e_)=qcmu*w1(ixO^S,e_)+qcnu*w2(ixO^S,e_)+(1.d0-qcmu-qcnu)*w(ixO^S,e_)&
            +qcmut*qdt*tmp(ixO^S)+qcnut*wold(ixO^S,e_)

if(fixsmall) call smallvalues(w2,x,ixI^L,ixO^L,'evolve_stepj')
end subroutine evolve_stepj
!=============================================================================
subroutine evolve_step1(qcmut,qdt,ixI^L,ixO^L,w1,w,x,wold)

use mod_global_parameters

integer, intent(in) :: ixI^L,ixO^L
double precision, intent(in) :: qcmut, qdt, w(ixI^S,1:nw), x(ixI^S,1:ndim)
double precision, intent(out) ::w1(ixI^S,1:nw),wold(ixI^S,1:nw)

double precision :: tmp(ixI^S),tmp1(ixI^S),tmp2(ixI^S),Te(ixI^S)
integer :: lowindex(ndim), ix^D
!-----------------------------------------------------------------------------
{^IFMHDPHYS
call heatconduct_mhd(tmp,tmp1,tmp2,ixI^L,ixO^L,w,x)
}
{^IFHDPHYS
call heatconduct_hd(tmp,tmp1,tmp2,ixI^L,ixO^L,w,x)
}

wold(ixO^S,e_)=qdt*tmp(ixO^S)
! update internal energy
tmp1(ixO^S) = tmp1(ixO^S) + qcmut*wold(ixO^S,e_)

! ensure you never trigger negative pressure 
! hence code up energy change with respect to kinetic and magnetic
! part(nonthermal)
if(smallT>0.d0) then
  Te(ixO^S)=tmp1(ixO^S)*(eqpar(gamma_)-1.d0)/w(ixO^S,rho_)
endif
if(strictsmall) then
  if(smallT>0.d0 .and. any(Te(ixO^S)<smallT)) then
    lowindex=minloc(Te(ixO^S))
    ^D&lowindex(^D)=lowindex(^D)+ixOmin^D-1;
    write(*,*)'too small temperature = ',minval(Te(ixO^S)),'at x=',&
   x(^D&lowindex(^D),1:ndim),lowindex,' with limit=',smallT,&
   ' on time=',t,' step=',it,' where density=',w(^D&lowindex(^D),rho_),&
   ' velocity=',&
    dsqrt((^C&w(^D&lowindex(^D),m^C_)**2+)/w(^D&lowindex(^D),rho_)**2)
{^IFMHDPHYS
    write(*,*)'magnetic field=',dsqrt((^C&w(^D&lowindex(^D),b^C_)**2+))
}
    call mpistop("==evolve_step1: too small temperature==")
  end if
  if(any(tmp1(ixO^S)<smalle)) then
    lowindex=minloc(tmp1(ixO^S))
    ^D&lowindex(^D)=lowindex(^D)+ixOmin^D-1;
    write(*,*)'too small internal energy = ',minval(tmp1(ixO^S)),'at x=',&
   x(^D&lowindex(^D),1:ndim),lowindex,' with limit=',smalle,&
   ' on time=',t,' step=',it,' where density=',w(^D&lowindex(^D),rho_),&
   ' velocity=',&
    dsqrt((^C&w(^D&lowindex(^D),m^C_)**2+)/w(^D&lowindex(^D),rho_)**2)
{^IFMHDPHYS
    write(*,*)'magnetic field=',dsqrt((^C&w(^D&lowindex(^D),b^C_)**2+))
}
    call mpistop("==evolve_step1: too small internal energy==")
  end if
  w1(ixO^S,e_) = tmp2(ixO^S)+tmp1(ixO^S)
else
 {do ix^DB=ixOmin^DB,ixOmax^DB\}
    if(smallT>0.d0) then
      if(Te(ix^D)<smallT) then
        w1(ix^D,e_) = tmp2(ix^D)+smallT*w(ix^D,rho_)/(eqpar(gamma_)-1.d0)
      else
        w1(ix^D,e_) = tmp2(ix^D)+tmp1(ix^D)
      end if
    else
      if(tmp1(ix^D)<smalle) then
        w1(ix^D,e_) = tmp2(ix^D)+smalle
      else
        w1(ix^D,e_) = tmp2(ix^D)+tmp1(ix^D)
      end if
    end if
 {end do\}
end if
!if(fixsmall) call smallvalues(w1,x,ixI^L,ixO^L,'evolve_step1')

end subroutine evolve_step1
!=============================================================================
subroutine addsource_heatconduct_mhd(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

use mod_global_parameters

integer, intent(in) :: ixI^L,ixO^L,iw^LIM
double precision, intent(in) :: qdt,qtC,qt, x(ixI^S,1:ndim),wCT(ixI^S,1:nw)
double precision, intent(inout) ::w(ixI^S,1:nw)
double precision :: tmp(ixI^S),tmp1(ixI^S),tmp2(ixI^S)
integer :: ix^D
integer, dimension(ndim)       :: lowindex
!------------------------------------------------------------------------------

if(e_<1) call mpistop("Heat conduction requires by e_>0!")
!heat conduct rate tmp, old internal energy tmp1, old nonthermal energy tmp2
call heatconduct_mhd(tmp,tmp1,tmp2,ixI^L,ixO^L,w,x)
! update internal energy
tmp1(ixO^S) = tmp1(ixO^S) + qdt*tmp(ixO^S)

! code up energy change with respect to kinetic and magnetic part(nonthermal)
if(strictsmall) then
  if(any(tmp1(ixO^S)<smalle)) then
    lowindex=minloc(tmp1(ixO^S))
    ^D&lowindex(^D)=lowindex(^D)+ixOmin^D-1;
    write(*,*)'too small internal energy = ',minval(tmp1(ixO^S)),'at x=',&
  x(^D&lowindex(^D),1:ndim),lowindex,' with limit=',smalle,&
 ' on time=',t,' step=',it,' where density=',w(^D&lowindex(^D),rho_),&
 ' velocity=',&
 dsqrt((^C&w(^D&lowindex(^D),m^C_)**2+)/w(^D&lowindex(^D),rho_)**2),&
    ' magnetic field=',dsqrt((^C&w(^D&lowindex(^D),b^C_)**2+))
    call mpistop("=strictsmall in sourceheatcond_mhd: too small energy==")
  else
    w(ixO^S,e_) = tmp2(ixO^S)+tmp1(ixO^S)
  endif
else
 {do ix^DB=ixOmin^DB,ixOmax^DB\}
    if(tmp1(ix^D)>smalle) then
      w(ix^D,e_) = tmp2(ix^D)+tmp1(ix^D)
    else
      w(ix^D,e_) = tmp2(ix^D)+smalle
    end if
 {end do\}
endif

end subroutine addsource_heatconduct_mhd
!=============================================================================
subroutine heatconduct_mhd(tmp,tmp1,tmp2,ixI^L,ixO^L,w,x)

use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L
double precision, intent(in) ::  x(ixI^S,1:ndim), w(ixI^S,1:nw)
!! tmp store the heat conduction energy changing rate
double precision, intent(out) :: tmp(ixI^S),tmp1(ixI^S),tmp2(ixI^S)

double precision, dimension(ixI^S,1:ndir) :: mf,qvec
double precision, dimension(ixI^S,1:ndim) :: gradT,qvecsat
{#IFDEF TCPERPENDICULAR
double precision, dimension(ixI^S,1:ndim) :: qvec_per, qvec_max
}
double precision, dimension(ixI^S) :: B2inv, BgradT,cs3,qflux,qsatflux
integer, dimension(ndim) :: lowindex
integer :: ix^L,idims,ix^D
logical :: Bnull(ixI^S)
!-----------------------------------------------------------------------------
ix^L=ixO^L^LADD2;
if (ixI^L^LTix^L|.or.|.or.) &
   call mpistop("Need two extra layers for thermal conduction")

! store kinetic+magnetic energy before addition of heat conduction source
tmp2(ixI^S)=half*( (^C&w(ixI^S,m^C_)**2+)/w(ixI^S,rho_)+ &
  (^C&w(ixI^S,b^C_)**2+) )

! Calculate einternal=e-0.5*(2ek+2eb)
tmp1(ixI^S)=w(ixI^S,e_)-tmp2(ixI^S)
! Clip off negative pressure if smallp is set
if(strictsmall) then
   if(any(tmp1(ixI^S)<smalle)) then
     lowindex=minloc(tmp1(ixI^S))
     ^D&lowindex(^D)=lowindex(^D)+ixImin^D-1;
     write(*,*)'too low internal energy = ',minval(tmp1(ixI^S)),' at x=',&
     x(^D&lowindex(^D),1:ndim),lowindex,' with limit=',smalle,' on time=',t,&
     ' step=',it,' where density=',w(^D&lowindex(^D),rho_),' velocity=',&
     dsqrt((^C&w(^D&lowindex(^D),m^C_)**2+)/w(^D&lowindex(^D),rho_)**2),&
     ' magnetic field=',dsqrt((^C&w(^D&lowindex(^D),b^C_)**2+))
     call mpistop("=== strictsmall in heatcond_mhd: low internal energy ===")
   end if
else
{do ix^DB=ixImin^DB,ixImax^DB\}
   if(tmp1(ix^D)<smalle) then
      tmp1(ix^D)=smalle
   end if
{end do\}
end if
! compute the temperature
tmp(ixI^S)=tmp1(ixI^S)*(eqpar(gamma_)-one)/w(ixI^S,rho_)
if(smallT>0.d0) then
  if(strictsmall) then
     if(any(tmp(ixI^S)<smallT)) then
       lowindex=minloc(tmp(ixI^S))
       ^D&lowindex(^D)=lowindex(^D)+ixImin^D-1;
       write(*,*)'too low temperature = ',minval(tmp(ixI^S)),' at x=',&
       x(^D&lowindex(^D),1:ndim),lowindex,' with limit=',smallT,' on time=',t,&
       ' step=',it,' where density=',w(^D&lowindex(^D),rho_),' velocity=',&
       dsqrt((^C&w(^D&lowindex(^D),m^C_)**2+)/w(^D&lowindex(^D),rho_)**2),&
       ' magnetic field=',dsqrt((^C&w(^D&lowindex(^D),b^C_)**2+))
       call mpistop("=== strictsmall in heatcond_mhd: low temperature ===")
     end if
  else
  {do ix^DB=ixImin^DB,ixImax^DB\}
     if(tmp(ix^D)<smallT) then
        tmp(ix^D)=smallT
     end if
  {end do\}
  end if
end if
do idims=1,ndim
  ! idirth component of gradient of temperature at cell center
  select case(typegrad)
  case("central")
    ix^L=ixI^L^LSUB1;
    call gradient(tmp,ixI^L,ix^L,idims,B2inv)
  case("limited")
    ix^L=ixI^L^LSUB2;
    call gradientS(tmp,ixI^L,ix^L,idims,B2inv)
  end select
  ! get grad T in all directions
  gradT(ix^S,idims)=B2inv(ix^S)
end do
! B
if(B0field) then
 ^C&mf(ix^S,^C)=w(ix^S,b^C_)+myB0_cell%w(ix^S,^C);
else
 ^C&mf(ix^S,^C)=w(ix^S,b^C_);
end if
! B^-2
B2inv(ix^S)= ^C&mf(ix^S,^C)**2+
Bnull(ixI^S)=.false.
where(B2inv(ix^S)/=0.d0)
  B2inv(ix^S)=1.d0/B2inv(ix^S)
elsewhere
  Bnull(ix^S)=.true.
end where

BgradT(ix^S)=(^D&mf(ix^S,^D)*gradT(ix^S,^D)+)*B2inv(ix^S)
cs3(ix^S)=eqpar(kappa_)*dsqrt(tmp(ix^S))*tmp(ix^S)**2*BgradT(ix^S)

do idims=1,ndim
  qvec(ix^S,idims)=mf(ix^S,idims)*cs3(ix^S)
end do

if(TCsaturate) then
  ! Cowie and Mckee 1977 ApJ, 211, 135
  cs3(ix^S)=dsqrt(tmp(ix^S)**3)
  ! unsigned saturated TC flux = 5 phi rho c**3
  qsatflux(ix^S)=5.d0*TCphi*w(ix^S,rho_)*cs3(ix^S)
  ! strength of classic TC flux
  qflux(ix^S)=dsqrt(^D&qvec(ix^S,^D)**2+)
  ! sign(b * Grad Te) 5 phi rho c**3 Bi/B 
  {do ix^DB=ixmin^DB,ixmax^DB\}
    if(.false. .and. qflux(ix^D)>qsatflux(ix^D) .and. idims==1) write(*,*)& 
      'it',it,' ratio=',qflux(ix^D)/qsatflux(ix^D),' TC saturated at ',&
      x(ix^D,:),' rho',w(ix^D,rho_),' Te',tmp(ix^D)
    if(qflux(ix^D)>qsatflux(ix^D)) then
    ! saturated TC flux = sign(b * Grad Te) 5 phi rho c**3
      qsatflux(ix^D)=sign(1.d0,BgradT(ix^D))*qsatflux(ix^D)*dsqrt(B2inv(ix^D))
      do idims=1,ndim
        qvec(ix^D,idims)=qsatflux(ix^D)*mf(ix^D,idims)
      end do
    end if
  {end do\}
end if
{#IFDEF TCPERPENDICULAR
! van der Linden and Goossens 1991 SoPh 131, 79; Orlando et al 2008 ApJ 678, 274
do idims=1,ndim
  ! q_per = kappe n^2 B^-2 Te^-0.5 (Grad Te - (e_b . Grad Te )e_b) 
  qvec_per(ix^S,idims)=eqpar(kappe_)*w(ix^S,rho_)**2*B2inv(ix^S)/dsqrt(tmp(ix^S))&
                       *(gradT(ix^S,idims)-BgradT(ix^S)*mf(ix^S,idims))
end do
! maximal thermal conducting (saturated) flux
do idims=1,ndim
  qvec_max(ix^S,idims)=eqpar(kappa_)*dsqrt(tmp(ix^S)**5)*gradT(ix^S,idims)
end do
if(TCsaturate) then
  qsatflux(ix^S)=5.d0*TCphi*w(ix^S,rho_)*cs3(ix^S)
  qflux(ix^S)=dsqrt(^D&qvec_max(ix^S,^D)**2+)
  {do ix^DB=ixmin^DB,ixmax^DB\}
    if(.false. .and. qflux(ix^D)>qsatflux(ix^D) .and. idims==1) write(*,*) & 
      'it',it,' ratio=',qflux(ix^D)/qsatflux(ix^D),' TC_PER saturated at ',&
      x(ix^D,:),' rho',w(ix^D,rho_),' Te',tmp(ix^D)
    if(qflux(ix^D)>qsatflux(ix^D)) then
      qsatflux(ix^D)=qsatflux(ix^D)/qflux(ix^D)
      do idims=1,ndim
        qvec_max(ix^D,idims)=qsatflux(ix^D)*qvec_max(ix^D,idims)
      end do
    end if
  {end do\}
end if
! maximal thermal conducting flux perpendicular to magnetic field
qvec_max(ix^S,1:ndim)=qvec_max(ix^S,1:ndim)-qvec(ix^S,1:ndim)
qsatflux(ix^S)= ^D&qvec_max(ix^S,^D)**2+
qflux(ix^S)= ^D&qvec_per(ix^S,^D)**2+
{do ix^DB=ixmin^DB,ixmax^DB\}
   if(qflux(ix^D)>qsatflux(ix^D)) then
     qvec(ix^D,1:ndim)=qvec(ix^D,1:ndim)+qvec_max(ix^D,1:ndim)
   else
     qvec(ix^D,1:ndim)=qvec(ix^D,1:ndim)+qvec_per(ix^D,1:ndim)
   end if   
{end do\}
}
{do ix^DB=ixmin^DB,ixmax^DB\}
  if(Bnull(ix^D)) then
    qvec(ix^D,1:ndim)=eqpar(kappa_)*dsqrt(tmp(ix^D)**5)*gradT(ix^D,1:ndim)
    if(TCsaturate) then
      ! unsigned saturated TC flux = 5 phi rho c**3
      qsatflux(ix^D)=5.d0*TCphi*w(ix^D,rho_)*cs3(ix^D)
      qflux(ix^D)=dsqrt(sum(qvec(ix^D,1:ndim)**2))
      if(qflux(ix^D)>qsatflux(ix^D)) then
        qsatflux(ix^D)=qsatflux(ix^D)/qflux(ix^D)
        do idims=1,ndim
          qvec(ix^D,idims)=qsatflux(ix^D)*qvec(ix^D,idims)
        end do
      end if
    end if
  end if  
{end do\}

select case(typediv)
   case("central")
      call divvector(qvec,ixI^L,ixO^L,tmp)
   case("limited")
      call divvectorS(qvec,ixI^L,ixO^L,tmp)
end select

end subroutine heatconduct_mhd
!=============================================================================
subroutine getdt_heatconduct_mhd(w,ixG^L,ix^L,dtnew,dx^D,x)

!Check diffusion time limit dt < dtTCpar*dx_i**2/((gamma-1)*kappa_i/rho)
!where                      kappa_i=kappa*B_i**2/B**2
!and                        T=p/rho

use mod_global_parameters

integer, intent(in) :: ixG^L, ix^L
double precision, intent(in) :: dx^D, x(ixG^S,1:ndim)
! note that depending on strictsmall etc, w values may change 
! through call to getpthermal
double precision, intent(inout) :: w(ixG^S,1:nw), dtnew

double precision :: dxinv(1:ndim),mf(ixG^S,1:ndir)
double precision :: tmp2(ixG^S),tmp(ixG^S),Te(ixG^S),B2inv(ixG^S)
double precision :: dtdiff_tcond, dtdiff_tsat
integer          :: idim,ix^D
integer, dimension(ndim)       :: lowindex
!-----------------------------------------------------------------------------
^D&dxinv(^D)=one/dx^D;

! Calculate pressure=(gamma-1)*(e-0.5*(2ek+2eb))
tmp(ix^S)=(eqpar(gamma_)-one)*(w(ix^S,e_)- &
       half*(({^C&w(ix^S,m^C_)**2+})/w(ix^S,rho_)&
       +{ ^C&w(ix^S,b^C_)**2+}))
! Clip off negative pressure if smallp is set
if(strictsmall) then
  if(any(tmp(ix^S)<minp)) then
    lowindex=minloc(tmp(ix^S))
    ^D&lowindex(^D)=lowindex(^D)+ixmin^D-1;
    write(*,*)'low pressure = ',minval(tmp(ix^S)),' at x=',&
    x(^D&lowindex(^D),1:ndim),lowindex,' with limit=',minp,' on time=',t,&
    ' step=',it
   call mpistop("=== strictsmall in getdt_heatconduct_mhd: low pressure ===")
  end if
else
{do ix^DB=ixmin^DB,ixmax^DB\}
   if(tmp(ix^D)<minp) then
      tmp(ix^D)=minp
   end if
{end do\}
end if
!temperature
Te(ix^S)=tmp(ix^S)/w(ix^S,rho_)
!kappa_i
tmp(ix^S)=eqpar(kappa_)*dsqrt(Te(ix^S)**5)
!(gamma-1)*kappa_i/rho
tmp(ix^S)=(eqpar(gamma_)-one)*tmp(ix^S)/w(ix^S,rho_)

! B
if(B0field) then
 ^C&mf(ix^S,^C)=w(ix^S,b^C_)+myB0_cell%w(ix^S,^C);
else
 ^C&mf(ix^S,^C)=w(ix^S,b^C_);
end if
! B^-2
B2inv(ix^S)= ^C&mf(ix^S,^C)**2+
where(B2inv(ix^S)/=0.d0)
  B2inv(ix^S)=1.d0/B2inv(ix^S)
end where
do idim=1,ndim
   ! B_i**2/B**2
   where(B2inv(ix^S)/=0.d0)
     tmp2(ix^S)=mf(ix^S,idim)**2*B2inv(ix^S)
   elsewhere
     tmp2(ix^S)=1.d0
   end where
   ! dt< dtTCpar * dx_idim**2/((gamma-1)*kappa_i/rho*B_i**2/B**2)
   dtdiff_tcond=dtTcpar/maxval(tmp(ix^S)*tmp2(ix^S)*dxinv(idim)**2)
   if(TCsaturate) then
     ! dt< dtTCpar* dx_idim**2/((gamma-1)*sqrt(Te)*5*phi)
     ! with an empirical coefficient dx_idim
     dtdiff_tsat=dtTCpar/maxval((eqpar(gamma_)-1.d0)*dsqrt(Te(ix^S))*&
                 5.d0*TCphi*dxinv(idim)**2)
     ! choose the slower flux (bigger time step) between classic and saturated
     dtdiff_tcond=max(dtdiff_tcond,dtdiff_tsat)
   end if
   ! limit the time step
   dtnew=min(dtnew,dtdiff_tcond)
end do
dtnew=dtnew/dble(ndim)

end subroutine getdt_heatconduct_mhd
!=============================================================================
subroutine addsource_heatconduct_hd(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

use mod_global_parameters

integer, intent(in) :: ixI^L,ixO^L,iw^LIM
double precision, intent(in) :: qdt,qtC,qt, x(ixI^S,1:ndim),wCT(ixI^S,1:nw)
double precision, intent(inout) ::w(ixI^S,1:nw)

double precision :: tmp(ixI^S),tmp1(ixI^S),tmp2(ixI^S)
integer:: ix^L,idim,ix^D
integer, dimension(ndim)       :: lowindex
!------------------------------------------------------------------------------
if(e_<1) call mpistop("Thermal conduction requires e_>0!")
!get heat conduct rate in tmp, old thermal energy in tmp1, old nonthermal energy in tmp2
call heatconduct_hd(tmp,tmp1,tmp2,ixI^L,ixO^L,w,x)
! update internal energy
tmp1(ixO^S) = tmp1(ixO^S) + qdt*tmp(ixO^S)

if(strictsmall) then
  if(any(tmp1(ixO^S)<smalle)) then
    lowindex=minloc(tmp1(ixO^S))
    ^D&lowindex(^D)=lowindex(^D)+ixOmin^D-1;
    write(*,*)'too small internal energy = ',minval(tmp1(ixO^S)),' at x=',&
    x(^D&lowindex(^D),1:ndim),lowindex,' with limit=',smalle,&
    ' on time=',t
    call mpistop("=== strictsmall in heatconduct: too small internal energy ===")
  end if
endif

! ensure you never trigger negative pressure 
! hence code up energy change with respect to kinetic and magnetic part
{do ix^DB=ixOmin^DB,ixOmax^DB\}
 if(tmp1(ix^D)>smalle) then
   w(ix^D,e_) = tmp2(ix^D)+tmp1(ix^D)
 else
   w(ix^D,e_) = tmp2(ix^D)+smalle
 endif
{end do\}

return
end subroutine addsource_heatconduct_hd
!=============================================================================
subroutine heatconduct_hd(tmp,tmp1,tmp2,ixI^L,ixO^L,w,x)

use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L
double precision, intent(in) ::  x(ixI^S,1:ndim), w(ixI^S,1:nw)
!! tmp store the heat conduction energy changing rate
double precision, intent(out) :: tmp(ixI^S),tmp1(ixI^S),tmp2(ixI^S)
double precision :: qvec(ixI^S,1:ndir),Te(ixI^S),qflux(ixI^S),qsatflux(ixI^S)
integer:: ix^L,idims,ix^D
integer, dimension(ndim)       :: lowindex
!-----------------------------------------------------------------------------

ix^L=ixO^L^LADD2;
if (ixI^L^LTix^L|.or.|.or.) &
   call mpistop("Need two extra layers for thermal conduction")

! store old kinetic energy
tmp2(ixI^S)=half*(^C&w(ixI^S,m^C_)**2+ )/w(ixI^S,rho_)
! Calculate pressure=(gamma-1)*(e-0.5*(2ek+2eb))
tmp(ixI^S)=(eqpar(gamma_)-one)*(w(ixI^S,e_)-tmp2(ixI^S))
! Clip off negative pressure if smallp is set
if(strictsmall) then
  if(any(tmp(ixI^S)<minp)) then
    lowindex=minloc(tmp(ixI^S))
    write(*,*)'low pressure = ',minval(tmp(ixI^S)),' at x=',&
    x(^D&lowindex(^D),1:ndim),lowindex,' with limit=',minp,' on time=',t
    call mpistop("=== strictsmall in heatconduct: low pressure ===")
  end if
else
   {do ix^DB=ixImin^DB,ixImax^DB\}
     if(tmp(ix^D)<minp) then
      tmp(ix^D)=minp
     end if
   {end do\}
end if
! store old internal energy
tmp1(ixI^S)=tmp(ixI^S)/(eqpar(gamma_)-one)

! compute temperature before source addition
Te(ixI^S)=tmp(ixI^S)/w(ixI^S,rho_)
if(smallT>0.d0) then
  if(strictsmall) then
     if(any(Te(ixI^S)<smallT)) then
       lowindex=minloc(Te(ixI^S))
       ^D&lowindex(^D)=lowindex(^D)+ixImin^D-1;
       write(*,*)'too low temperature = ',minval(Te(ixI^S)),' at x=',&
       x(^D&lowindex(^D),1:ndim),lowindex,' with limit=',smallT,' on time=',t,&
       ' step=',it,' where density=',w(^D&lowindex(^D),rho_),' velocity=',&
       dsqrt((^C&w(^D&lowindex(^D),m^C_)**2+)/w(^D&lowindex(^D),rho_)**2)
       call mpistop("=== strictsmall in heatcond_hd: low temperature ===")
     end if
  else
  {do ix^DB=ixImin^DB,ixImax^DB\}
     if(Te(ix^D)<smallT) then
        Te(ix^D)=smallT
     end if
  {end do\}
  end if
end if
! compute grad T and store grad T vector
do idims=1,ndim
  ! idirth component of gradient of temperature at cell center
   select case(typegrad)
   case("central")
     ix^L=ixI^L^LSUB1;
     call gradient(Te,ixI^L,ix^L,idims,tmp)
   case("limited")
     ix^L=ixI^L^LSUB2;
     call gradientS(Te,ixI^L,ix^L,idims,tmp)
   end select
   qvec(ix^S,idims)=tmp(ix^S)*eqpar(kappa_)*dsqrt(Te(ix^S)**5)
end do

if(TCsaturate) then
  ! unsigned saturated TC flux = 5 phi rho c**3
  qsatflux(ix^S)=5.d0*TCphi*w(ix^S,rho_)*dsqrt(Te(ix^S)**3)
  qflux(ix^S)=dsqrt(^D&qvec(ix^S,^D)**2+)
  {do ix^DB=ixmin^DB,ixmax^DB\}
    if(qflux(ix^D)>qsatflux(ix^D)) then
      qsatflux(ix^D)=qsatflux(ix^D)/qflux(ix^D)
      do idims=1,ndim
        qvec(ix^D,idims)=qsatflux(ix^D)*qvec(ix^D,idims)
      end do
    end if
  {end do\}
end if
select case(typediv)
   case("central")
      call divvector(qvec,ixI^L,ixO^L,tmp)
   case("limited")
      call divvectorS(qvec,ixI^L,ixO^L,tmp)
end select
end subroutine heatconduct_hd
!=============================================================================
subroutine getdt_heatconduct_hd(w,ixG^L,ix^L,dtnew,dx^D,x)

! Check diffusion time limit dt < dtTCpar * dx_i**2 / ((gamma-1)*kappa_i/rho)

use mod_global_parameters

integer, intent(in) :: ixG^L, ix^L
double precision, intent(in) :: dx^D, x(ixG^S,1:ndim)
! note that depending on strictsmall etc, w values may change 
! through call to getpthermal
double precision, intent(inout) :: w(ixG^S,1:nw), dtnew

double precision :: dxinv(1:ndim), tmp(ixG^S), Te(ixG^S)
double precision :: dtdiff_tcond,dtdiff_tsat
integer          :: idim,ix^D
integer, dimension(ndim)       :: lowindex
!-----------------------------------------------------------------------------
^D&dxinv(^D)=one/dx^D;

! Determine pressure and then the direction independent part of
tmp(ix^S)=( ^C&w(ix^S,m^C_)**2+ )/w(ix^S,rho_)
! Calculate pressure=(gamma-1)*(e-0.5*(2ek))
tmp(ix^S)=(eqpar(gamma_)-one)*(w(ix^S,e_)-half*tmp(ix^S))
! Clip off negative pressure if smallp is set
if(strictsmall) then
  if(any(tmp(ix^S)<minp)) then
    lowindex=minloc(tmp(ix^S))
    ^D&lowindex(^D)=lowindex(^D)+ixmin^D-1;
    write(*,*)'low pressure = ',minval(tmp(ix^S)),' at x=',&
    x(^D&lowindex(^D),1:ndim),lowindex,' with limit=',minp,' on time=',t
    call mpistop("=== strictsmall in getdt_heatconduct_hd: low pressure ===")
  end if
else
  {do ix^DB=ixmin^DB,ixmax^DB\}
    if(tmp(ix^D)<minp) then
     tmp(ix^D)=minp
    end if
  {end do\}
endif
Te(ix^S)=tmp(ix^S)/w(ix^S,rho_)
tmp(ix^S)=(eqpar(gamma_)-one)*eqpar(kappa_)*dsqrt((Te(ix^S))**5)/w(ix^S,rho_)

do idim=1,ndim
   ! dt< dtTCpar * dx_idim**2/((gamma-1)*kappa_idim/rho)
   dtdiff_tcond=dtTCpar/maxval(tmp(ix^S)*dxinv(idim)**2)
   if(TCsaturate) then
     ! dt< dtTCpar* dx_idim**2/((gamma-1)*sqrt(Te)*5*phi)
     dtdiff_tsat=dtTCpar/maxval((eqpar(gamma_)-1.d0)*dsqrt(Te(ix^S))*&
                 5.d0*TCphi*dxinv(idim)**2)
     ! choose the slower flux (bigger time scale) between classic and saturated
     dtdiff_tcond=max(dtdiff_tcond,dtdiff_tsat)
   end if
   ! limit the time step
   dtnew=min(dtnew,dtdiff_tcond)
enddo
dtnew=dtnew/dble(ndim)

end subroutine getdt_heatconduct_hd
!=============================================================================
