!{\src2tex{textfont=tt}}
!!****f* ABINIT/xcpbe
!! NAME
!! xcpbe
!!
!! FUNCTION
!! Treat XC functionals closely linked with the Perdew-Wang 92 LSD
!! and the PBE GGA.
!!
!! For a series of values of the density and, if GGA, the square of the
!! gradient of the density, return the associated Exc energy,
!! potential, and, in case of response-function, functions needed
!! to build the XC kernel.
!!
!! If option==2, Exchange-correlation functional from Perdew-Burke-Ernzerhof,
!! Phys.Rev.Lett. 77, 3866 (1996).
!! If option==1, Reduces to Perdew-Wang LSD , PRB45,13244 (1992).
!! If option==-1 or -2, take only exchange part of PW (-1) or PBE (-2) functionals.
!! If option==3 or 4, take exchange plus RPA correlation
!!   part of LSD PW (3) or GGA PBE (4) functionals.
!! If option==5, revPBE functional of Zhang and Yang, PRL 80, 890 (1998)
!! If option==6, RPBE functional of Hammer, Hansen and Norskov, PRB 59, 7413 (1999)
!! If option==7, WC functional of Wu and Cohen, PRB 73, 235116 (2006)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2007 ABINIT group (XG,MF,LG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  exexch= choice of local exact exchange. Active if exexch=3
!!  npts= number of points to be computed
!!  nspden=1 for unpolarized, 2 for spin-polarized
!!  grho2_updn(npts,ngr2)=square of the gradient of the spin-up,
!!     and, if nspden==2, spin-down, and total density (Hartree/Bohr**2),
!!     only used if gradient corrected functional (option=2,-2 and 4 or beyond)
!!  option= see above
!!  order=its absolute value gives the maximal derivative of Exc to be computed.
!!  rho_updn(npts,nspden)=spin-up and spin-down density (Hartree/bohr**3)
!!  ndvxci= size of dvxci(npts,ndvxci)
!!  ngr2= size of grho2_updn(npts,ngr2)
!!
!! OUTPUT
!!
!!  d2vxci=third derivative of the xc energy with respect to the density, only
!!    only if local-density approximation, and non-spin-polarized case
!!    (nsppol=1)
!!  dvxcdgr(npts,3)=partial derivative of the exchange-correlation
!!    energy (exci*$\rho$) with respect to the spin-up (dvxcdgr(:,1)),
!!    spin-down (dvxcdgr(:,2)), or total spin (dvxcdgr(:,3)) gradients of the density
!!    divided by the norm of the gradient (the definition changed in v3.3)
!!
!!  dvxci=partial second derivatives of the xc energy, only if abs(order)>1
!!   In case of local energy functional (option=1,-1 or 3):
!!    dvxci(npts,1+nspden)=              (Hartree*bohr^3)
!!     if(nspden=1 .and. order==2): dvxci(:,1)=dvxc/d$\rho$ , dvxc(:,2) empty
!!     if(nspden=1 .and. order==-2): also compute dvxci(:,2)=dvxc($\uparrow$)/d$\rho(\downarrow)$
!!     if(nspden=2): dvxci(:,1)=dvxc($\uparrow$)/d$\rho(\downarrow)$,
!!                   dvxci(:,2)=dvxc($\uparrow$)/d$\rho(\downarrow)$,
!!                   dvxci(:,3)=dvxc($\downarrow$)/d$\rho(\downarrow)$
!!   In case of gradient corrected functional (option=2,-2, 4, 5, 6):
!!    dvxci(npts,15)=
!!     dvxci(:,1)= d2Ex/drho_up drho_up
!!     dvxci(:,2)= d2Ex/drho_dn drho_dn
!!     dvxci(:,3)= dEx/d(abs(grad(rho_up))) / abs(grad(rho_up))
!!     dvxci(:,4)= dEx/d(abs(grad(rho_dn))) / abs(grad(rho_dn))
!!     dvxci(:,5)= d2Ex/d(abs(grad(rho_up))) drho_up / abs(grad(rho_up))
!!     dvxci(:,6)= d2Ex/d(abs(grad(rho_dn))) drho_dn / abs(grad(rho_dn))
!!     dvxci(:,7)= 1/abs(grad(rho_up)) * d/drho_up (dEx/d(abs(grad(rho_up))) /abs(grad(rho_up)))
!!     dvxci(:,8)= 1/abs(grad(rho_dn)) * d/drho_dn (dEx/d(abs(grad(rho_dn))) /abs(grad(rho_dn)))
!!     dvxci(:,9)= d2Ec/drho_up drho_up
!!     dvxci(:,10)=d2Ec/drho_up drho_dn
!!     dvxci(:,11)=d2Ec/drho_dn drho_dn
!!     dvxci(:,12)=dEc/d(abs(grad(rho))) / abs(grad(rho))
!!     dvxci(:,13)=d2Ec/d(abs(grad(rho))) drho_up / abs(grad(rho))
!!     dvxci(:,14)=d2Ec/d(abs(grad(rho))) drho_dn / abs(grad(rho))
!!     dvxci(:,15)=1/abs(grad(rho)) * d/drho (dEc/d(abs(grad(rho))) /abs(grad(rho)))
!!
!!  exci(npts)=exchange-correlation energy density (hartree)
!!  vxci(npts,nspden)=partial derivative of the exchange-correlation energy (exci*$\rho$)
!!    with respect to the spin-down (vxci(:,1)) and spin-up (vxci(:,2) densities
!! Normalization: Exc=$\int (exc(r)*\rho (r) d^3 r)$ for $\rho$(r)=electron density.
!!
!! TODO
!!  WARNING: option=4 not yet implemented.
!!
!! NOTES
!!
!! PARENTS
!!      drivexc
!!
!! CHILDREN
!!      invcb
!!
!! SOURCE
!!$#if defined HAVE_CONFIG_H
!!$#include "config.h"
!!$#endif

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine xcpbe(exci,npts,nspden,option,order,rho_updn,vxci,ndvxci,ngr2,& !Mandatory Arguments
&                d2vxci,dvxcdgr,dvxci,exexch,grho2_updn)                          !Optional Arguments

 use defs_basis

!This section has been created automatically by the script Abilint (TD). Do not modify these by hand.
#ifdef HAVE_FORTRAN_INTERFACES
 use interfaces_01manage_mpi
 use interfaces_13xc, except_this_one => xcpbe
#endif
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ndvxci,ngr2,npts,nspden,option,order
 integer,intent(in),optional :: exexch
!arrays
 real(dp),intent(in) :: rho_updn(npts,nspden)
 real(dp),intent(in),optional :: grho2_updn(npts,ngr2)
 real(dp),intent(out) :: exci(npts),vxci(npts,nspden)
 real(dp),intent(out),optional :: d2vxci(npts),dvxcdgr(npts,3)
 real(dp),intent(out),optional :: dvxci(npts,ndvxci)

!Local variables-------------------------------
! The "accurate" value of mu is taken from the PBE code
!scalars
 integer,save :: initialized=0
 integer :: ipts,ispden
 real(dp),parameter :: alpha_zeta2=1.0_dp-1.0e-6_dp,alpha_zeta=1.0_dp-1.0e-6_dp
 real(dp),parameter :: b_wc=0.123456790123_dp,beta=0.066725_dp
 real(dp),parameter :: beta_inv=1.0_dp/beta,c_wc=0.00793746933516_dp
 real(dp),parameter :: fsec_inv=1.0_dp/1.709921_dp
 real(dp),parameter :: kappa_pbe=0.804_dp,kappa_revpbe=1.245_dp
 real(dp),parameter :: mu=0.2195149727645171_dp,mu_divkappa_pbe=mu/kappa_pbe
 real(dp),parameter :: mu_divkappa_revpbe=mu/kappa_revpbe
 real(dp),parameter :: rsfac=0.6203504908994000_dp,tolgrad=tol10
 real(dp),save :: beta_gamma,coeff_tt,factf_zeta,factfp_zeta,gamma,gamma_inv
 real(dp),save :: sixpi2_1_3,sixpi2m1_3,sq_rsfac,sq_rsfac_inv,threefourth_divpi
 real(dp),save :: twom1_3
 real(dp) :: aa,arg_rr,bb,cc,coeff_aa,coeff_qq,coeffss,d2aa_drs2,d2aa_drsdzeta
 real(dp) :: d2aa_dzeta2,d2bb_drs2,d2bb_drsdzeta,d2bb_dzeta2,d2cc_dbb2,d_wc
 real(dp) :: d2cc_drs2,d2cc_drsdzeta,d2cc_dzeta2,d2ecrs0_drs2,d2ecrs1_drs2
 real(dp) :: d2ecrs_drdn2,d2ecrs_drdndrup,d2ecrs_drho2,d2ecrs_drs2
 real(dp) :: d2ecrs_drsdzeta,d2ecrs_drup2,d2ecrs_dzeta2,d2fxdg2,d2fxdn2
 real(dp) :: d2fxdndg,d2fxdss2,d2fzeta4_dzeta2,d2gcrs_drs2,d2hh_drs2
 real(dp) :: d2hh_drsdtt,d2hh_drsdzeta,d2hh_dtt2,d2hh_dttdzeta,d2hh_dzeta2
 real(dp) :: d2macrs_drs2,d2pade_drs2,d2pade_drsdtt,d2pade_drsdzeta,d2pade_dtt2
 real(dp) :: d2pade_dttdzeta,d2pade_dxx2,d2pade_dzeta2,d2qq_drs2,d2qq_drsdtt
 real(dp) :: d2qq_drsdzeta,d2qq_dtt2,d2qq_dttdzeta,d2qq_dzeta2,d2rhohh_drho2
 real(dp) :: d2rhohh_drhodg,d2rr_dqq2,d2rr_drs2,d2rr_drsdtt,d2rr_drsdzeta
 real(dp) :: d2rr_dtt2,d2rr_dttdzeta,d2rr_dzeta2,d2rs_dn2,d2ssdn2,d2ssdndg
 real(dp) :: d2vcrs_drs2,d2xx_drs2,d2xx_drsdtt,d2xx_drsdzeta,d2xx_dttdzeta
 real(dp) :: d2xx_dzeta2,d3ecrs0_drs3,daa_drs,daa_dzeta,dbb_drs,dbb_dzeta
 real(dp) :: dcc_dbb,dcc_drs,dcc_dzeta,decrs0_drs,decrs1_drs,decrs_drs
 real(dp) :: decrs_dzeta,dfxdg,dfxdn,dfxdss,dfzeta4_dzeta,dgcrs_drs
 real(dp) :: dhh_drs,dhh_dtt,dhh_dzeta,div_rr,divss,dmacrs_drs,dpade_drs
 real(dp) :: dpade_dtt,dpade_dxx,dpade_dzeta,dqq_drs,dqq_dtt,dqq_dzeta
 real(dp) :: drhohh_drho,drr_dqq,drr_drs,drr_dtt,drr_dzeta,drs_dn,dssdg,dssdn
 real(dp) :: dtt_dg,dvcrs_drs,dxx_drs,dxx_dtt,dxx_dzeta,ec0_a1,ec0_aa,ec0_b1
 real(dp) :: ec0_b2,ec0_b3,ec0_b4,ec0_den,ec0_f1,ec0_f2,ec0_log,ec0_q0,ec0_q1
 real(dp) :: ec0_q1p,ec0_q1pp,ec0_q1ppp,ec1_a1,ec1_aa,ec1_b1,ec1_b2,ec1_b3
 real(dp) :: ec1_b4,ec1_den,ec1_log,ec1_q0,ec1_q1,ec1_q1p,ec1_q1pp,ecrs,ecrs0
 real(dp) :: ecrs1,ex_gga,ex_lsd,exc,exp_pbe,expss,f_zeta,factfpp_zeta
 real(dp) :: fp_zeta,fpp_zeta,fx,gamphi3inv,gcrs,grrho2,hh,kappa,mac_a1
 real(dp) :: mac_aa,mac_b1,mac_b2,mac_b3,mac_b4,mac_den,mac_log,mac_q0,mac_q1
 real(dp) :: mac_q1p,mac_q1pp,macrs,mu_divkappa,p1_wc,p2_wc,pade,pade_den
 real(dp) :: phi3_zeta,phi_logder,phi_zeta,phi_zeta_inv,phip_zeta,phipp_zeta,qq
 real(dp) :: rho,rho_inv,rhomot
 real(dp) :: rhotmo6,rhotmot,rhoto6,rhotot,rhotot_inv,rr,rs,rsm1_2,sqr_rs
 real(dp) :: sqr_sqr_rs,ss,tt,vxcadd,xx,zeta,zeta4,zetm_1_3,zetp_1_3
 character(len=500) :: message
!arrays
 real(dp),allocatable :: rho_updnm1_3(:,:),rhoarr(:),rhom1_3(:),zetm(:)
 real(dp),allocatable :: zetmm1_3(:),zetp(:),zetpm1_3(:)

! *************************************************************************

!DEBUG
!write(6,*)' xcpbe : enter'
!write(6,*)' nspden=',nspden
!ENDDEBUG

 d_wc=mu-b_wc

 if (option<=-3 .or. option==0 .or. option==4 .or. option>=8 ) then
  write(message, '(a,a,a,a,i12,a)' ) ch10,&
&   ' xcpbe : BUG -',ch10,&
&   '  Option must be 1, 2, 3, 5, 6, 7, -1 or -2 ; argument was ',option,'.'
!  call wrtout(06,message,'COLL')
!  call leave_new('COLL')
 end if
!Checks the compatibility between the presence of dvxci and ndvxci
 if(ndvxci /=0 .neqv. present(dvxci))then
  write(message, '(4a)' )ch10,&
&  ' xcpbe : BUG -',ch10,&
&  '  If ndvxci/=0 there must the optional arguments dvxci'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if
!Checks the compatibility between the inputs and the presence of the optional arguments
 if(ndvxci /= 0 .and. abs(order) <= 1)then
  write(message, '(6a,i8,a)' )ch10,&
&  ' xcpbe : BUG -',ch10,&
&  '  The order does not require the presence of dvxci',ch10,&
&  '  that is allowed when |order|>1, while we have',&
&       order,'.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if
 if(ndvxci /= 0 .and. (&
&      ((option == 1 .or. option == -1 .or. option == 3) .and. ndvxci /= nspden + 1)&
&      .or. (option == -2 .and. ndvxci /= 8)&
&      .or. ((option == 2 .or. option == 5 .or. option == 6 .or. option == 7) .and. ndvxci /= 15)&
&   ))then
  write(message, '(16a,i8,a,i8,a,i8)' )ch10,&
&  ' xcpbe : BUG -',ch10,&
&  '  The option is not consistent with the value of ndvxci',ch10,&
&  '  Allowed values are:',ch10,&
&  '  ndvxci     option',ch10,&
&  ' nspden+1    1,-1,3',ch10,&
&  '    8          -2',ch10,&
&  '    15       2, 5,6,7',ch10,&
&  '  While we have: order=',order,'option=',option,'nspden=',nspden
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if
 if (present(grho2_updn)) then
  if (ngr2/=2*nspden-1 ) then
   write(message, '(4a)' ) ch10,&
&    ' xcpbe : BUG -',ch10,&
&    '  ngr2 must be 2*nspden-1 !'
   call wrtout(06,message,'COLL')
   call leave_new('COLL')
  end if
 end if

 if ((option == 1 .or. option == -1 .or. option ==3) .and.  (present(grho2_updn) .or. present(dvxcdgr))) then
  write(message, '(a,a,a,a,a,a,i6,a)' )ch10,&
&  ' xcpbe : BUG -',ch10,&
&  '  The option chosen does not need the presence',ch10,&
&  '  of the gradient, or of the array dvxcdgr in the input, needed if option/=1,-1,3 , while we have',&
&       option,'.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if
 if (order /= 3 .and. present(d2vxci)) then
  write(message, '(a,a,a,a,a,a,i6,a)' )ch10,&
&  ' xcpbe : BUG -',ch10,&
&  '  The order chosen does not need the presence',ch10,&
&  '  of the array d2vxci, needed if order=3 , while we have',&
&       order,'.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if


 if(initialized==0)then
    twom1_3=two**(-third)
    sixpi2_1_3=(six*pi**2)**third
    sixpi2m1_3=one/sixpi2_1_3
    threefourth_divpi=three_quarters*piinv
    gamma=(one-log(two))*piinv**2
    gamma_inv=one/gamma
    beta_gamma=beta*gamma_inv
    factf_zeta= one / ( two**(four/three)-two )
    factfp_zeta= four_thirds * factf_zeta * alpha_zeta2
    coeff_tt= one/ (four*four*piinv*(three*pi**2)**third)
    ! coeff_tt= two * sqrt(four*piinv*(three*pi**2)**third)
    sq_rsfac=sqrt(rsfac)
    sq_rsfac_inv=one/sq_rsfac
    initialized=1
 end if

 !Parameters for the Perdew-Wang 92 LSD as well as LSD-RPA,
 !see Table I of Phys.Rev.B 45,13244 (1992)
 ec0_aa=0.031091_dp  ; ec1_aa=0.015545_dp ; mac_aa=0.016887_dp
 if(option/=3 .and. option/=4)then
    ec0_a1=0.21370_dp  ; ec1_a1=0.20548_dp  ; mac_a1=0.11125_dp
    ec0_b1=7.5957_dp   ; ec1_b1=14.1189_dp  ; mac_b1=10.357_dp
    ec0_b2=3.5876_dp   ; ec1_b2=6.1977_dp   ; mac_b2=3.6231_dp
    ec0_b3=1.6382_dp   ; ec1_b3=3.3662_dp   ; mac_b3=0.88026_dp
    ec0_b4=0.49294_dp  ; ec1_b4=0.62517_dp  ; mac_b4=0.49671_dp
 else  ! RPA values
    ec0_a1=0.082477_dp ; ec1_a1=0.035374_dp ; mac_a1=0.028829_dp
    ec0_b1=5.1486_dp   ; ec1_b1=6.4869_dp   ; mac_b1=10.357_dp
    ec0_b2=1.6483_dp   ; ec1_b2=1.3083_dp   ; mac_b2=3.6231_dp
    ec0_b3=0.23647_dp  ; ec1_b3=0.11518_dp  ; mac_b3=0.479_dp
    ec0_b4=0.20614_dp  ; ec1_b4=0.082349_dp ; mac_b4=0.112279_dp
 end if

 if(option/=5)then
    kappa=kappa_pbe
    mu_divkappa=mu_divkappa_pbe
 else
    kappa=kappa_revpbe
    mu_divkappa=mu_divkappa_revpbe
 end if


 !DEBUG
 !Finite-difference debugging, do not take away
 !Note : here work with collinear gradients. Might be generalized ...
 !debug=2  ! Choose 1 (rho grads) or 2 (grho grads)
 !factor=1.0_dp
 !zeta_mean=0.98_dp
 !zeta_mean=zero
 !delta=0.000025*factor
 !delta=0.0000125*factor
 !if(debug/=0)then
 ! do ipts=1,npts,5
 !  rho=ipts*0.01_dp*factor
 !  rho_up=rho*(1.0_dp+zeta_mean)*0.5_dp
 !  rho_dn=rho*(1.0_dp-zeta_mean)*0.5_dp
 !  rho_upp=rho_up+delta
 !  rho_upm=rho_up-delta
 !  rho_dnp=rho_dn+delta
 !  rho_dnm=rho_dn-delta
 !! Here, vary rho
 !  if(debug==1)then
 !   rho_updn(ipts  ,1)=rho_up ; rho_updn(ipts  ,2)=rho_dn
 !   rho_updn(ipts+1,1)=rho_upp; rho_updn(ipts+1,2)=rho_dn
 !   rho_updn(ipts+2,1)=rho_upm; rho_updn(ipts+2,2)=rho_dn
 !   rho_updn(ipts+3,1)=rho_up ; rho_updn(ipts+3,2)=rho_dnp
 !   rho_updn(ipts+4,1)=rho_up ; rho_updn(ipts+4,2)=rho_dnm
 !   grho2_updn(ipts:ipts+4,1)=(0.2_dp*factor)**2     ! grad2 of spin up density
 !   grho2_updn(ipts:ipts+4,2)=(0.2_dp*factor)**2     ! grad2 of spin down density
 !   grho2_updn(ipts:ipts+4,3)=(0.3_dp*factor)**2     ! grad2 of total density
 !  else
 !!  Here, vary grho (interchange rho and grho)
 !   grho2_updn(ipts  ,1)=rho_up**2 ; grho2_updn(ipts  ,2)=rho_dn**2
 !   grho2_updn(ipts+1,1)=rho_upp**2; grho2_updn(ipts+1,2)=rho_dn**2
 !   grho2_updn(ipts+2,1)=rho_upm**2; grho2_updn(ipts+2,2)=rho_dn**2
 !   grho2_updn(ipts+3,1)=rho_up**2 ; grho2_updn(ipts+3,2)=rho_dnp**2
 !   grho2_updn(ipts+4,1)=rho_up**2 ; grho2_updn(ipts+4,2)=rho_dnm**2
 !   grho2_updn(ipts  ,3)=(ipts*0.01_dp*factor)**2
 !   grho2_updn(ipts+1,3)=(ipts*0.01_dp*factor+delta)**2
 !   grho2_updn(ipts+2,3)=(ipts*0.01_dp*factor-delta)**2
 !   grho2_updn(ipts+3,3)=(ipts*0.01_dp*factor+delta)**2   ! identical to ipts+1
 !   grho2_updn(ipts+4,3)=(ipts*0.01_dp*factor-delta)**2   ! identical to ipts+2
 !   rho_updn(ipts:ipts+4,1)=0.2_dp*factor*(1.0_dp+zeta_mean)*0.5_dp    ! spin up density
 !   rho_updn(ipts:ipts+4,2)=0.2_dp*factor*(1.0_dp-zeta_mean)*0.5_dp    ! spin down density
 !  end if
 ! end do
 !end if
 !Usual option :
 !nspden=2 ; order=2
 !GGA
 !nspden=2 ; order=1
 !Might take also, although finite difference later is meaningless
 !nspden=1 ; order=-2
 !ENDDEBUG

 if(order**2 >1)then
    factfpp_zeta= third * factfp_zeta * alpha_zeta2
 end if


 allocate(rhoarr(npts),rhom1_3(npts),rho_updnm1_3(npts,2))
 allocate(zetm(npts),zetmm1_3(npts),zetp(npts),zetpm1_3(npts))

 do ispden=1,nspden
    call invcb(rho_updn(:,ispden),rho_updnm1_3(:,ispden),npts)
 end do


 if(nspden==1)then
    rhoarr(:)=two*rho_updn(:,1)
    rhom1_3(:)=twom1_3*rho_updnm1_3(:,1)
    rho_updnm1_3(:,2)=rho_updnm1_3(:,1)
 else
    rhoarr(:)=rho_updn(:,1)+rho_updn(:,2)
    call invcb(rhoarr,rhom1_3,npts)
    do ipts=1,npts
       rhotmot=rhom1_3(ipts)
       rhotot_inv=rhotmot*rhotmot*rhotmot
       zeta=(rho_updn(ipts,1)-rho_updn(ipts,2))*rhotot_inv
       zetp(ipts)=1.0_dp+zeta*alpha_zeta
       zetm(ipts)=1.0_dp-zeta*alpha_zeta
    end do
    call invcb(zetp,zetpm1_3,npts)
    call invcb(zetm,zetmm1_3,npts)
 end if


 if (order==3 .and. nspden == 1) d2vxci(:)=0._dp

!!!Loop unrolling summary
! Completely unrolled for spin non-polarized case
! To be optimized for spin-polarized cases
! The loops are unrolled as follows:
! nspden=1
!       order^2<=1
!             option=2,5
!             option=6
!             option=-1
!             option=-2
!             option=1
!             option=3
!       order=3
!             option=2,5
!             option=6
!             option=-1
!             option=-2
!             option=1
!             option=3
!       order=-2
!             option=2,5
!             option=6
!             option=-1
!             option=-2
!             option=1
!             option=3
!       order^2>1
!             option=2,5
!             option=6
!             option=-1
!             option=-2
!             option=1
!             option=3
! nspden=2
!       order^2<=1
!       order^2>1 (with if statements inside distinguishing between order=3 or -2)
!!!End loop unrolling summary

!we separate different cases, depending on nspden
 if (nspden==1) then
    !we separate different cases, depending on order
    if (order**2<=1) then
       !we separate different cases, depending on option
       if(option==2 .or. option==5)then

          do ipts=1,npts

             rhotot=rhoarr(ipts)
             rhotmot=rhom1_3(ipts)
             rhotot_inv=rhotmot*rhotmot*rhotmot
             rhotmo6=sqrt(rhotmot)
             rhoto6=rhotot*rhotmot*rhotmot*rhotmo6
             ! -----------------------------------------------------------------------
             ! First take care of the exchange part of the functional

             exc=zero
             dvxcdgr(ipts,3)=zero
             !loop over the spin
             ispden=1
             rho   =rho_updn(ipts,ispden)
             rhomot=rho_updnm1_3(ipts,ispden)
             ex_lsd= - threefourth_divpi * sixpi2_1_3*rhomot*rhomot*rho
             !   Perdew-Burke-Ernzerhof GGA, exchange part
             rho_inv=rhomot*rhomot*rhomot
             coeffss=quarter*sixpi2m1_3*sixpi2m1_3*rho_inv*rho_inv*rhomot*rhomot
             ss=grho2_updn(ipts,ispden)*coeffss
             divss=one/(one+mu_divkappa*ss)
             dfxdss= mu*divss*divss
             d2fxdss2=-mu*two*mu_divkappa*divss*divss*divss
             fx    = one+kappa*(one-divss)
             ex_gga= ex_lsd*fx
             dssdn=-eight*third*ss*rho_inv
             dfxdn  = dfxdss*dssdn
             vxci(ipts,ispden)=ex_lsd*(four_thirds*fx+rho*dfxdn)
             !   The new definition (v3.3) includes the division by the norm of the gradient
             dssdg =two*coeffss
             dfxdg=dfxdss*dssdg
             dvxcdgr(ipts,ispden)=ex_lsd*rho*dfxdg
             exc=exc+ex_gga*rho

             !end of loop over the spin
             !  If non spin-polarized, treat spin down contribution now, similar to spin up
             exc=exc*2
             exci(ipts)=exc*rhotot_inv
             if(present(exexch).and.exexch==1) cycle
             ! -----------------------------------------------------------------------------
             ! Then takes care of the LSD correlation part of the functional


             rs=rsfac*rhotmot
             sqr_rs=sq_rsfac*rhotmo6
             rsm1_2=sq_rsfac_inv*rhoto6

             !  Formulas A6-A8 of PW92LSD
             ec0_q0=-2.0_dp*ec0_aa*(1.0_dp+ec0_a1*rs)
             ec0_q1=2.0_dp*ec0_aa*(ec0_b1*sqr_rs+ec0_b2*rs+ec0_b3*rs*sqr_rs+ec0_b4*rs*rs)
             ec0_q1p=ec0_aa*(ec0_b1*rsm1_2+2._dp*ec0_b2+3._dp*ec0_b3*sqr_rs+4._dp*ec0_b4*rs)
             ec0_den=1.0_dp/(ec0_q1*ec0_q1+ec0_q1)
             !  ec0_log=log( 1.0_dp + 1.0_dp / ec0_q1 )
             ec0_log=-log( ec0_q1*ec0_q1*ec0_den )
             ecrs0=ec0_q0*ec0_log
             decrs0_drs= -2.0_dp*ec0_aa*ec0_a1*ec0_log - ec0_q0*ec0_q1p *ec0_den

             ecrs=ecrs0
             decrs_drs=decrs0_drs
             decrs_dzeta=0.0_dp
             zeta=0.0_dp

             !  Add LSD correlation functional to GGA exchange functional
             exci(ipts)=exci(ipts)+ecrs
             vxci(ipts,1)=vxci(ipts,1)+ecrs-rs*third*decrs_drs


             !  -----------------------------------------------------------------------------
             !  Eventually add the GGA correlation part of the PBE functional
             !  Note : the computation of the potential in the spin-unpolarized
             !  case could be optimized much further. Other optimizations are left to do.

             phi_zeta=1.0_dp
             phip_zeta=0.0_dp
             phi_zeta_inv=1.0_dp
             phi_logder=0.0_dp
             phi3_zeta=1.0_dp
             gamphi3inv=gamma_inv
             phipp_zeta=-two*ninth*alpha_zeta*alpha_zeta

             !   From ec to bb
             bb=ecrs*gamphi3inv
             dbb_drs=decrs_drs*gamphi3inv
             dbb_dzeta=gamphi3inv*(decrs_dzeta-three*ecrs*phi_logder)

             !   From bb to cc
             exp_pbe=exp(-bb)
             cc=one/(exp_pbe-one)
             dcc_dbb=cc*cc*exp_pbe
             dcc_drs=dcc_dbb*dbb_drs
             dcc_dzeta=dcc_dbb*dbb_dzeta

             !   From cc to aa
             coeff_aa=beta*gamma_inv*phi_zeta_inv*phi_zeta_inv
             aa=coeff_aa*cc
             daa_drs=coeff_aa*dcc_drs
             daa_dzeta=-two*aa*phi_logder+coeff_aa*dcc_dzeta

             !   Introduce tt : do not assume that the spin-dependent gradients are collinear
             grrho2=four*grho2_updn(ipts,1)
             dtt_dg=two*rhotot_inv*rhotot_inv*rhotmot*coeff_tt
             !   Note that tt is (the t variable of PBE divided by phi) squared
             tt=half*grrho2*dtt_dg

             !   Get xx from aa and tt
             xx=aa*tt
             dxx_drs=daa_drs*tt
             dxx_dzeta=daa_dzeta*tt
             dxx_dtt=aa

             !   From xx to pade
             pade_den=one/(one+xx*(one+xx))
             pade=(one+xx)*pade_den
             dpade_dxx=-xx*(two+xx)*pade_den**2
             dpade_drs=dpade_dxx*dxx_drs
             dpade_dtt=dpade_dxx*dxx_dtt
             dpade_dzeta=dpade_dxx*dxx_dzeta

             !   From pade to qq
             coeff_qq=tt*phi_zeta_inv*phi_zeta_inv
             qq=coeff_qq*pade
             dqq_drs=coeff_qq*dpade_drs
             dqq_dtt=pade*phi_zeta_inv*phi_zeta_inv+coeff_qq*dpade_dtt
             dqq_dzeta=coeff_qq*(dpade_dzeta-two*pade*phi_logder)

             !   From qq to rr
             arg_rr=one+beta*gamma_inv*qq
             div_rr=one/arg_rr
             rr=gamma*log(arg_rr)
             drr_dqq=beta*div_rr
             drr_drs=drr_dqq*dqq_drs
             drr_dtt=drr_dqq*dqq_dtt
             drr_dzeta=drr_dqq*dqq_dzeta

             !   From rr to hh
             hh=phi3_zeta*rr
             dhh_drs=phi3_zeta*drr_drs
             dhh_dtt=phi3_zeta*drr_dtt
             dhh_dzeta=phi3_zeta*(drr_dzeta+three*rr*phi_logder)

             !   The GGA correlation energy is added
             exci(ipts)=exci(ipts)+hh

             !   Change of variables : from (rs,zeta,tt) to (rhoup,rhodn,grrho)

             !   From hh to the derivative of the energy wrt the density
             drhohh_drho=hh-third*rs*dhh_drs-zeta*dhh_dzeta-seven*third*tt*dhh_dtt
             vxci(ipts,1)=vxci(ipts,1)+drhohh_drho

             !   From hh to the derivative of the energy wrt to the gradient of the
             !   density, divided by the gradient of the density
             !   (The v3.3 definition includes the division by the norm of the gradient)
             dvxcdgr(ipts,3)=rhotot*dtt_dg*dhh_dtt

             !  End condition of GGA

             ! Correlation has been added
             ! -----------------------------------------------------------------------------

             !  vxci(ipts,2)=vxci(ipts,1)
             dvxcdgr(ipts,2)=dvxcdgr(ipts,1)

          end do
       else if((option==6) .or. (option==7)) then

          do ipts=1,npts

             rhotot=rhoarr(ipts)
             rhotmot=rhom1_3(ipts)
             rhotot_inv=rhotmot*rhotmot*rhotmot
             rhotmo6=sqrt(rhotmot)
             rhoto6=rhotot*rhotmot*rhotmot*rhotmo6
             ! -----------------------------------------------------------------------
             ! First take care of the exchange part of the functional

             exc=zero
             dvxcdgr(ipts,3)=zero
             !loop over the spin
             ispden=1
             rho   =rho_updn(ipts,ispden)
             rhomot=rho_updnm1_3(ipts,ispden)
             ex_lsd= - threefourth_divpi * sixpi2_1_3*rhomot*rhomot*rho
             !   Perdew-Burke-Ernzerhof GGA, exchange part
             rho_inv=rhomot*rhomot*rhomot
             coeffss=quarter*sixpi2m1_3*sixpi2m1_3*rho_inv*rho_inv*rhomot*rhomot
             ss=grho2_updn(ipts,ispden)*coeffss

             ! This is RPBE modification
             if (option==6) then
                divss=exp(-mu_divkappa*ss)
                dfxdss= mu*divss
                d2fxdss2=-mu*mu_divkappa*divss

                fx    = one+kappa*(one-divss)
                ex_gga= ex_lsd*fx
                dssdn=-eight*third*ss*rho_inv
                dfxdn  = dfxdss*dssdn
                vxci(ipts,ispden)=ex_lsd*(four_thirds*fx+rho*dfxdn)
                !   The new definition (v3.3) includes the division by the norm of the gradient
                dssdg =two*coeffss
                dfxdg=dfxdss*dssdg
                dvxcdgr(ipts,ispden)=ex_lsd*rho*dfxdg
                exc=exc+ex_gga*rho
             ! This is the Wu and Cohen modification
             else
                expss=exp(-ss)
                p1_wc=b_wc+(mu-b_wc)*(one-ss)*expss+two*c_wc*ss/(one+c_wc*ss*ss)
                p2_wc=d_wc*(ss-two)*expss+two*c_wc/(one+c_wc*ss*ss)-&
&                     four*c_wc*c_wc*ss*ss/((one+c_wc*ss*ss)*(one+c_wc*ss*ss))
                divss=one/(one+(b_wc*ss+d_wc*ss*expss+log(one+c_wc*ss*ss))/kappa)
                dfxdss=p1_wc*divss*divss
                d2fxdss2=p2_wc*divss*divss-two*divss*divss*divss*p1_wc*p1_wc/kappa

                fx    = one+kappa*(one-divss)
                ex_gga= ex_lsd*fx
                dssdn=-eight*third*ss*rho_inv
                dfxdn  = dfxdss*dssdn
                vxci(ipts,ispden)=ex_lsd*(four_thirds*fx+rho*dfxdn)
                !   The new definition (v3.3) includes the division by the norm of the gradient
                dssdg =two*coeffss
                dfxdg=dfxdss*dssdg
                dvxcdgr(ipts,ispden)=ex_lsd*rho*dfxdg
                exc=exc+ex_gga*rho
             end if

             !end of loop over the spin
             !  If non spin-polarized, treat spin down contribution now, similar to spin up
             exc=exc*2
             exci(ipts)=exc*rhotot_inv
             ! -----------------------------------------------------------------------------
             ! Then takes care of the LSD correlation part of the functional


             rs=rsfac*rhotmot
             sqr_rs=sq_rsfac*rhotmo6
             rsm1_2=sq_rsfac_inv*rhoto6

             !  Formulas A6-A8 of PW92LSD
             ec0_q0=-2.0_dp*ec0_aa*(1.0_dp+ec0_a1*rs)
             ec0_q1=2.0_dp*ec0_aa*(ec0_b1*sqr_rs+ec0_b2*rs+ec0_b3*rs*sqr_rs+ec0_b4*rs*rs)
             ec0_q1p=ec0_aa*(ec0_b1*rsm1_2+2._dp*ec0_b2+3._dp*ec0_b3*sqr_rs+4._dp*ec0_b4*rs)
             ec0_den=1.0_dp/(ec0_q1*ec0_q1+ec0_q1)
             !  ec0_log=log( 1.0_dp + 1.0_dp / ec0_q1 )
             ec0_log=-log( ec0_q1*ec0_q1*ec0_den )
             ecrs0=ec0_q0*ec0_log
             decrs0_drs= -2.0_dp*ec0_aa*ec0_a1*ec0_log - ec0_q0*ec0_q1p *ec0_den

             ecrs=ecrs0
             decrs_drs=decrs0_drs
             decrs_dzeta=0.0_dp
             zeta=0.0_dp

             !  Add LSD correlation functional to GGA exchange functional
             exci(ipts)=exci(ipts)+ecrs
             vxci(ipts,1)=vxci(ipts,1)+ecrs-rs*third*decrs_drs


             !  -----------------------------------------------------------------------------
             !  Eventually add the GGA correlation part of the PBE functional
             !  Note : the computation of the potential in the spin-unpolarized
             !  case could be optimized much further. Other optimizations are left to do.

             phi_zeta=1.0_dp
             phip_zeta=0.0_dp
             phi_zeta_inv=1.0_dp
             phi_logder=0.0_dp
             phi3_zeta=1.0_dp
             gamphi3inv=gamma_inv
             phipp_zeta=-two*ninth*alpha_zeta*alpha_zeta

             !   From ec to bb
             bb=ecrs*gamphi3inv
             dbb_drs=decrs_drs*gamphi3inv
             dbb_dzeta=gamphi3inv*(decrs_dzeta-three*ecrs*phi_logder)

             !   From bb to cc
             exp_pbe=exp(-bb)
             cc=one/(exp_pbe-one)
             dcc_dbb=cc*cc*exp_pbe
             dcc_drs=dcc_dbb*dbb_drs
             dcc_dzeta=dcc_dbb*dbb_dzeta

             !   From cc to aa
             coeff_aa=beta*gamma_inv*phi_zeta_inv*phi_zeta_inv
             aa=coeff_aa*cc
             daa_drs=coeff_aa*dcc_drs
             daa_dzeta=-two*aa*phi_logder+coeff_aa*dcc_dzeta

             !   Introduce tt : do not assume that the spin-dependent gradients are collinear
             grrho2=four*grho2_updn(ipts,1)
             dtt_dg=two*rhotot_inv*rhotot_inv*rhotmot*coeff_tt
             !   Note that tt is (the t variable of PBE divided by phi) squared
             tt=half*grrho2*dtt_dg

             !   Get xx from aa and tt
             xx=aa*tt
             dxx_drs=daa_drs*tt
             dxx_dzeta=daa_dzeta*tt
             dxx_dtt=aa

             !   From xx to pade
             pade_den=one/(one+xx*(one+xx))
             pade=(one+xx)*pade_den
             dpade_dxx=-xx*(two+xx)*pade_den**2
             dpade_drs=dpade_dxx*dxx_drs
             dpade_dtt=dpade_dxx*dxx_dtt
             dpade_dzeta=dpade_dxx*dxx_dzeta

             !   From pade to qq
             coeff_qq=tt*phi_zeta_inv*phi_zeta_inv
             qq=coeff_qq*pade
             dqq_drs=coeff_qq*dpade_drs
             dqq_dtt=pade*phi_zeta_inv*phi_zeta_inv+coeff_qq*dpade_dtt
             dqq_dzeta=coeff_qq*(dpade_dzeta-two*pade*phi_logder)

             !   From qq to rr
             arg_rr=one+beta*gamma_inv*qq
             div_rr=one/arg_rr
             rr=gamma*log(arg_rr)
             drr_dqq=beta*div_rr
             drr_drs=drr_dqq*dqq_drs
             drr_dtt=drr_dqq*dqq_dtt
             drr_dzeta=drr_dqq*dqq_dzeta

             !   From rr to hh
             hh=phi3_zeta*rr
             dhh_drs=phi3_zeta*drr_drs
             dhh_dtt=phi3_zeta*drr_dtt
             dhh_dzeta=phi3_zeta*(drr_dzeta+three*rr*phi_logder)

             !   The GGA correlation energy is added
             exci(ipts)=exci(ipts)+hh

             !   Change of variables : from (rs,zeta,tt) to (rhoup,rhodn,grrho)

             !   From hh to the derivative of the energy wrt the density
             drhohh_drho=hh-third*rs*dhh_drs-zeta*dhh_dzeta-seven*third*tt*dhh_dtt
             vxci(ipts,1)=vxci(ipts,1)+drhohh_drho

             !   From hh to the derivative of the energy wrt to the gradient of the
             !   density, divided by the gradient of the density
             !   (The v3.3 definition includes the division by the norm of the gradient)
             dvxcdgr(ipts,3)=rhotot*dtt_dg*dhh_dtt

             !  End condition of GGA

             ! Correlation has been added
             ! -----------------------------------------------------------------------------

             !  vxci(ipts,2)=vxci(ipts,1)
             dvxcdgr(ipts,2)=dvxcdgr(ipts,1)

          end do


       else if (option==-1) then

          do ipts=1,npts

             rhotot=rhoarr(ipts)
             rhotmot=rhom1_3(ipts)
             rhotot_inv=rhotmot*rhotmot*rhotmot
             rhotmo6=sqrt(rhotmot)
             rhoto6=rhotot*rhotmot*rhotmot*rhotmo6
             ! -----------------------------------------------------------------------
             ! First take care of the exchange part of the functional

             exc=zero
             !loop over the spin
             ispden=1
             rho   =rho_updn(ipts,ispden)
             rhomot=rho_updnm1_3(ipts,ispden)
             ex_lsd= - threefourth_divpi * sixpi2_1_3*rhomot*rhomot*rho
             !   Perdew_Wang 91 LSD
             vxci(ipts,ispden)=four_thirds*ex_lsd
             exc=exc+ex_lsd*rho

             !end of loop over the spin
             !  If non spin-polarized, treat spin down contribution now, similar to spin up
             exc=exc*2
             exci(ipts)=exc*rhotot_inv
          end do

       else if(option==-2) then


          do ipts=1,npts

             rhotot=rhoarr(ipts)
             rhotmot=rhom1_3(ipts)
             rhotot_inv=rhotmot*rhotmot*rhotmot
             rhotmo6=sqrt(rhotmot)
             rhoto6=rhotot*rhotmot*rhotmot*rhotmo6
             ! -----------------------------------------------------------------------
             ! First take care of the exchange part of the functional

             exc=zero
             dvxcdgr(ipts,3)=zero
             !loop over the spin
             ispden=1
             rho   =rho_updn(ipts,ispden)
             rhomot=rho_updnm1_3(ipts,ispden)
             ex_lsd= - threefourth_divpi * sixpi2_1_3*rhomot*rhomot*rho
             !   Perdew-Burke-Ernzerhof GGA, exchange part
             rho_inv=rhomot*rhomot*rhomot
             coeffss=quarter*sixpi2m1_3*sixpi2m1_3*rho_inv*rho_inv*rhomot*rhomot
             ss=grho2_updn(ipts,ispden)*coeffss
             divss=one/(one+mu_divkappa*ss)
             dfxdss= mu*divss*divss
             d2fxdss2=-mu*two*mu_divkappa*divss*divss*divss
             fx    = one+kappa*(one-divss)
             ex_gga= ex_lsd*fx
             dssdn=-eight*third*ss*rho_inv
             dfxdn  = dfxdss*dssdn
             vxci(ipts,ispden)=ex_lsd*(four_thirds*fx+rho*dfxdn)
             !   The new definition (v3.3) includes the division by the norm of the gradient
             dssdg =two*coeffss
             dfxdg=dfxdss*dssdg
             dvxcdgr(ipts,ispden)=ex_lsd*rho*dfxdg
             exc=exc+ex_gga*rho

             !end of loop over the spin
             !  If non spin-polarized, treat spin down contribution now, similar to spin up
             exc=exc*2
             exci(ipts)=exc*rhotot_inv

             ! Correlation has been added
             ! -----------------------------------------------------------------------------

             !  vxci(ipts,2)=vxci(ipts,1)
             dvxcdgr(ipts,2)=dvxcdgr(ipts,1)

          end do

       else if(option==1)then


          do ipts=1,npts

             rhotot=rhoarr(ipts)
             rhotmot=rhom1_3(ipts)
             rhotot_inv=rhotmot*rhotmot*rhotmot
             rhotmo6=sqrt(rhotmot)
             rhoto6=rhotot*rhotmot*rhotmot*rhotmo6
             ! -----------------------------------------------------------------------
             ! First take care of the exchange part of the functional

             exc=zero
             !loop over the spin
             ispden=1
             rho   =rho_updn(ipts,ispden)
             rhomot=rho_updnm1_3(ipts,ispden)
             ex_lsd= - threefourth_divpi * sixpi2_1_3*rhomot*rhomot*rho
             !   Perdew_Wang 91 LSD
             vxci(ipts,ispden)=four_thirds*ex_lsd
             exc=exc+ex_lsd*rho

             !end of loop over the spin
             !  If non spin-polarized, treat spin down contribution now, similar to spin up
             exc=exc*2
             exci(ipts)=exc*rhotot_inv
             ! -----------------------------------------------------------------------------
             ! Then takes care of the LSD correlation part of the functional


             rs=rsfac*rhotmot
             sqr_rs=sq_rsfac*rhotmo6
             rsm1_2=sq_rsfac_inv*rhoto6

             !  Formulas A6-A8 of PW92LSD
             ec0_q0=-2.0_dp*ec0_aa*(1.0_dp+ec0_a1*rs)
             ec0_q1=2.0_dp*ec0_aa*(ec0_b1*sqr_rs+ec0_b2*rs+ec0_b3*rs*sqr_rs+ec0_b4*rs*rs)
             ec0_q1p=ec0_aa*(ec0_b1*rsm1_2+2._dp*ec0_b2+3._dp*ec0_b3*sqr_rs+4._dp*ec0_b4*rs)
             ec0_den=1.0_dp/(ec0_q1*ec0_q1+ec0_q1)
             !  ec0_log=log( 1.0_dp + 1.0_dp / ec0_q1 )
             ec0_log=-log( ec0_q1*ec0_q1*ec0_den )
             ecrs0=ec0_q0*ec0_log
             decrs0_drs= -2.0_dp*ec0_aa*ec0_a1*ec0_log - ec0_q0*ec0_q1p *ec0_den

             ecrs=ecrs0
             decrs_drs=decrs0_drs
             decrs_dzeta=0.0_dp
             zeta=0.0_dp

             !  Add LSD correlation functional to GGA exchange functional
             exci(ipts)=exci(ipts)+ecrs
             vxci(ipts,1)=vxci(ipts,1)+ecrs-rs*third*decrs_drs

             ! Correlation has been added
             ! -----------------------------------------------------------------------------

          end do

       else if (option==3) then
          do ipts=1,npts

             rhotot=rhoarr(ipts)
             rhotmot=rhom1_3(ipts)
             rhotot_inv=rhotmot*rhotmot*rhotmot
             rhotmo6=sqrt(rhotmot)
             rhoto6=rhotot*rhotmot*rhotmot*rhotmo6
             ! -----------------------------------------------------------------------
             ! First take care of the exchange part of the functional

             exc=zero
             !loop over the spin
             ispden=1
             rho   =rho_updn(ipts,ispden)
             rhomot=rho_updnm1_3(ipts,ispden)
             ex_lsd= - threefourth_divpi * sixpi2_1_3*rhomot*rhomot*rho
             !   Perdew_Wang 91 LSD
             vxci(ipts,ispden)=four_thirds*ex_lsd
             exc=exc+ex_lsd*rho

             !end of loop over the spin
             !  If non spin-polarized, treat spin down contribution now, similar to spin up
             exc=exc*2
             exci(ipts)=exc*rhotot_inv
             ! -----------------------------------------------------------------------------
             ! Then takes care of the LSD correlation part of the functional


             rs=rsfac*rhotmot
             sqr_rs=sq_rsfac*rhotmo6
             rsm1_2=sq_rsfac_inv*rhoto6

             !  Formulas A6-A8 of PW92LSD
             ec0_q0=-2.0_dp*ec0_aa*(1.0_dp+ec0_a1*rs)
             sqr_sqr_rs=max(1.e-15_dp,sqrt(sqr_rs))
             ec0_q1=2.0_dp*ec0_aa*(ec0_b1*sqr_rs+ec0_b2*rs+ec0_b3*rs*sqr_rs+ec0_b4*rs*rs/sqr_sqr_rs)
             ec0_q1p=ec0_aa*(ec0_b1*rsm1_2+2._dp*ec0_b2+3._dp*ec0_b3*sqr_rs+3.5_dp*ec0_b4*rs/sqr_sqr_rs)
             ec0_den=1.0_dp/(ec0_q1*ec0_q1+ec0_q1)
             !  ec0_log=log( 1.0_dp + 1.0_dp / ec0_q1 )
             ec0_log=-log( ec0_q1*ec0_q1*ec0_den )
             ecrs0=ec0_q0*ec0_log
             decrs0_drs= -2.0_dp*ec0_aa*ec0_a1*ec0_log - ec0_q0*ec0_q1p *ec0_den

             ecrs=ecrs0
             decrs_drs=decrs0_drs
             decrs_dzeta=0.0_dp
             zeta=0.0_dp

             !  Add LSD correlation functional to GGA exchange functional
             exci(ipts)=exci(ipts)+ecrs
             vxci(ipts,1)=vxci(ipts,1)+ecrs-rs*third*decrs_drs

             ! Correlation has been added
             ! -----------------------------------------------------------------------------

          end do


       end if

    else if (order==3) then
       !separate cases with respect to option
       if(option==2 .or. option==5) then

          do ipts=1,npts

             rhotot=rhoarr(ipts)
             rhotmot=rhom1_3(ipts)
             rhotot_inv=rhotmot*rhotmot*rhotmot
             rhotmo6=sqrt(rhotmot)
             rhoto6=rhotot*rhotmot*rhotmot*rhotmo6
             ! -----------------------------------------------------------------------
             ! First take care of the exchange part of the functional

             exc=zero
             dvxcdgr(ipts,3)=zero
             !loop over the spin
             ispden=1
             rho   =rho_updn(ipts,ispden)
             rhomot=rho_updnm1_3(ipts,ispden)
             ex_lsd= - threefourth_divpi * sixpi2_1_3*rhomot*rhomot*rho
             !   Perdew-Burke-Ernzerhof GGA, exchange part
             rho_inv=rhomot*rhomot*rhomot
             coeffss=quarter*sixpi2m1_3*sixpi2m1_3*rho_inv*rho_inv*rhomot*rhomot
             ss=grho2_updn(ipts,ispden)*coeffss
             divss=one/(one+mu_divkappa*ss)
             dfxdss= mu*divss*divss
             d2fxdss2=-mu*two*mu_divkappa*divss*divss*divss
             fx    = one+kappa*(one-divss)
             ex_gga= ex_lsd*fx
             dssdn=-eight*third*ss*rho_inv
             dfxdn  = dfxdss*dssdn
             vxci(ipts,ispden)=ex_lsd*(four_thirds*fx+rho*dfxdn)
             !   The new definition (v3.3) includes the division by the norm of the gradient
             dssdg =two*coeffss
             dfxdg=dfxdss*dssdg
             dvxcdgr(ipts,ispden)=ex_lsd*rho*dfxdg
             exc=exc+ex_gga*rho

             !    Perdew-Burke-Ernzerhof GGA, exchange part
             !    Components 3 or 4
             dvxci(ipts,2+ispden)=dvxcdgr(ipts,ispden)
             !    Components 1 or 2
             d2ssdn2=-11.0_dp*third*dssdn*rho_inv
             d2fxdn2=d2fxdss2*dssdn**2+dfxdss*d2ssdn2
             dvxci(ipts,ispden)=third*rho_inv*vxci(ipts,ispden)+&
                  &     ex_lsd*(seven*third*dfxdn+rho*d2fxdn2)
             !    Components 5 or 6
             d2ssdndg=-eight*third*dssdg*rho_inv
             d2fxdndg=d2fxdss2*dssdn*dssdg+dfxdss*d2ssdndg
             dvxci(ipts,4+ispden)=ex_lsd*(four_thirds*dfxdg+rho*d2fxdndg)
             !    Components 7 or 8
             d2fxdg2=d2fxdss2*dssdg**2
             dvxci(ipts,6+ispden)=ex_lsd*rho*d2fxdg2
             !    For the time being, treat non-spin-polarized like spin-polarized
             dvxci(ipts,2)=dvxci(ipts,1)
             dvxci(ipts,4)=dvxci(ipts,3)
             dvxci(ipts,6)=dvxci(ipts,5)
             dvxci(ipts,8)=dvxci(ipts,7)

             !end of loop over the spin
             !  If non spin-polarized, treat spin down contribution now, similar to spin up
             exc=exc*2
             exci(ipts)=exc*rhotot_inv
             ! -----------------------------------------------------------------------------
             ! Then takes care of the LSD correlation part of the functional


             rs=rsfac*rhotmot
             sqr_rs=sq_rsfac*rhotmo6
             rsm1_2=sq_rsfac_inv*rhoto6

             !  Formulas A6-A8 of PW92LSD
             ec0_q0=-2.0_dp*ec0_aa*(1.0_dp+ec0_a1*rs)
             ec0_q1=2.0_dp*ec0_aa*(ec0_b1*sqr_rs+ec0_b2*rs+ec0_b3*rs*sqr_rs+ec0_b4*rs*rs)
             ec0_q1p=ec0_aa*(ec0_b1*rsm1_2+2._dp*ec0_b2+3._dp*ec0_b3*sqr_rs+4._dp*ec0_b4*rs)
             ec0_den=1.0_dp/(ec0_q1*ec0_q1+ec0_q1)
             !  ec0_log=log( 1.0_dp + 1.0_dp / ec0_q1 )
             ec0_log=-log( ec0_q1*ec0_q1*ec0_den )
             ecrs0=ec0_q0*ec0_log
             decrs0_drs= -2.0_dp*ec0_aa*ec0_a1*ec0_log - ec0_q0*ec0_q1p *ec0_den
             ec0_q1pp=0.5_dp*ec0_aa*(-ec0_b1*rsm1_2**3+3._dp*ec0_b3*rsm1_2+8._dp*ec0_b4)
             d2ecrs0_drs2= 4.0_dp*ec0_aa*ec0_a1*ec0_q1p*ec0_den            &
                  &                -ec0_q0*ec0_q1pp*ec0_den                        &
                  &                +ec0_q0*ec0_q1p**2*ec0_den**2*(2._dp*ec0_q1+1.0_dp)
             ec0_q1ppp = 0.75_dp*ec0_aa*(rsm1_2**5)*(ec0_b1-ec0_b3*rs)
             ec0_f1 = 1._dp/(ec0_q1*ec0_q1*(1._dp + ec0_q1))
             ec0_f2 = 1._dp/(ec0_q1*(1+ec0_q1))
             d3ecrs0_drs3 = 6._dp*ec0_q1p*ec0_f1*(-2._dp*ec0_aa*ec0_a1*ec0_q1p + &
                  &        ec0_q0*ec0_q1pp) - &
                  &        ec0_f2*(-6._dp*ec0_aa*ec0_a1*ec0_q1pp + ec0_q0*ec0_q1ppp + &
                  &        ec0_f2*(3._dp*ec0_q1p*(-2._dp*ec0_aa*ec0_a1*ec0_q1p + ec0_q0*ec0_q1pp) + &
                  &        ec0_f2*2._dp*ec0_q0*(ec0_q1p**3)*(1._dp + 3._dp*ec0_q1*(1._dp + ec0_q1))))

             mac_q0=-2.0_dp*mac_aa*(1.0_dp+mac_a1*rs)
             mac_q1=2.0_dp*mac_aa*(mac_b1*sqr_rs+mac_b2*rs+mac_b3*rs*sqr_rs+mac_b4*rs*rs)
             mac_q1p=mac_aa*(mac_b1*rsm1_2+2._dp*mac_b2+3._dp*mac_b3*sqr_rs+4._dp*mac_b4*rs)
             mac_den=1.0_dp/(mac_q1*mac_q1+mac_q1)
             mac_log=-log( mac_q1*mac_q1*mac_den )
             macrs=mac_q0*mac_log
             dmacrs_drs= -2.0_dp*mac_aa*mac_a1*mac_log - mac_q0*mac_q1p*mac_den

             ecrs=ecrs0
             decrs_drs=decrs0_drs
             decrs_dzeta=0.0_dp
             d2ecrs_drs2=d2ecrs0_drs2
             d2ecrs_dzeta2=alpha_zeta**2*(-macrs)
             d2ecrs_drsdzeta=zero
             zeta=0.0_dp


             !  Add LSD correlation functional to GGA exchange functional
             exci(ipts)=exci(ipts)+ecrs
             vxci(ipts,1)=vxci(ipts,1)+ecrs-rs*third*decrs_drs

             dvcrs_drs=third*(2._dp*decrs_drs-rs*d2ecrs_drs2)
             !   And d(vxc)/d(rho)=(-rs/(3*rho))*d(vxc)/d(rs)
             d2ecrs_drho2= -rs**4*(four_pi*third)*third*dvcrs_drs
             dvxci(ipts,9)=d2ecrs_drho2
             dvxci(ipts,10)=d2ecrs_drho2
             dvxci(ipts,11)=d2ecrs_drho2

             !  -----------------------------------------------------------------------------
             !  Eventually add the GGA correlation part of the PBE functional
             !  Note : the computation of the potential in the spin-unpolarized
             !  case could be optimized much further. Other optimizations are left to do.

             phi_zeta=1.0_dp
             phip_zeta=0.0_dp
             phi_zeta_inv=1.0_dp
             phi_logder=0.0_dp
             phi3_zeta=1.0_dp
             gamphi3inv=gamma_inv
             phipp_zeta=-two*ninth*alpha_zeta*alpha_zeta

             !   From ec to bb
             bb=ecrs*gamphi3inv
             dbb_drs=decrs_drs*gamphi3inv
             dbb_dzeta=gamphi3inv*(decrs_dzeta-three*ecrs*phi_logder)
             d2bb_drs2=d2ecrs_drs2*gamphi3inv
             d2bb_drsdzeta=gamphi3inv*(d2ecrs_drsdzeta-three*decrs_drs*phi_logder)
             d2bb_dzeta2=gamphi3inv*(d2ecrs_dzeta2-six*decrs_dzeta*phi_logder+&
                  &     12.0_dp*ecrs*phi_logder*phi_logder-three*ecrs*phi_zeta_inv*phipp_zeta)

             !   From bb to cc
             exp_pbe=exp(-bb)
             cc=one/(exp_pbe-one)
             dcc_dbb=cc*cc*exp_pbe
             dcc_drs=dcc_dbb*dbb_drs
             dcc_dzeta=dcc_dbb*dbb_dzeta
             d2cc_dbb2=cc*cc*exp_pbe*(two*cc*exp_pbe-one)
             d2cc_drs2=d2cc_dbb2*dbb_drs*dbb_drs+dcc_dbb*d2bb_drs2
             d2cc_drsdzeta=d2cc_dbb2*dbb_drs*dbb_dzeta+dcc_dbb*d2bb_drsdzeta
             d2cc_dzeta2=d2cc_dbb2*dbb_dzeta*dbb_dzeta+dcc_dbb*d2bb_dzeta2

             !   From cc to aa
             coeff_aa=beta*gamma_inv*phi_zeta_inv*phi_zeta_inv
             aa=coeff_aa*cc
             daa_drs=coeff_aa*dcc_drs
             daa_dzeta=-two*aa*phi_logder+coeff_aa*dcc_dzeta
             d2aa_drs2=coeff_aa*d2cc_drs2
             d2aa_drsdzeta=-two*daa_drs*phi_logder+coeff_aa*d2cc_drsdzeta
             d2aa_dzeta2=aa*(-two*phi_zeta_inv*phipp_zeta+six*phi_logder*phi_logder)+&
                  &     coeff_aa*(-four*dcc_dzeta*phi_logder+d2cc_dzeta2)

             !   Introduce tt : do not assume that the spin-dependent gradients are collinear
             grrho2=four*grho2_updn(ipts,1)
             dtt_dg=two*rhotot_inv*rhotot_inv*rhotmot*coeff_tt
             !   Note that tt is (the t variable of PBE divided by phi) squared
             tt=half*grrho2*dtt_dg

             !   Get xx from aa and tt
             xx=aa*tt
             dxx_drs=daa_drs*tt
             dxx_dzeta=daa_dzeta*tt
             dxx_dtt=aa
             d2xx_drs2=d2aa_drs2*tt
             d2xx_drsdzeta=d2aa_drsdzeta*tt
             d2xx_drsdtt=daa_drs
             d2xx_dttdzeta=daa_dzeta
             d2xx_dzeta2=d2aa_dzeta2*tt

             !   From xx to pade
             pade_den=one/(one+xx*(one+xx))
             pade=(one+xx)*pade_den
             dpade_dxx=-xx*(two+xx)*pade_den**2
             dpade_drs=dpade_dxx*dxx_drs
             dpade_dtt=dpade_dxx*dxx_dtt
             dpade_dzeta=dpade_dxx*dxx_dzeta
             d2pade_dxx2=two*(-one+xx*xx*(three+xx))*pade_den*pade_den*pade_den
             d2pade_drs2=d2pade_dxx2*dxx_drs*dxx_drs+dpade_dxx*d2xx_drs2
             d2pade_drsdtt=d2pade_dxx2*dxx_drs*dxx_dtt+dpade_dxx*d2xx_drsdtt
             d2pade_drsdzeta=d2pade_dxx2*dxx_drs*dxx_dzeta+dpade_dxx*d2xx_drsdzeta
             d2pade_dtt2=d2pade_dxx2*dxx_dtt*dxx_dtt
             d2pade_dttdzeta=d2pade_dxx2*dxx_dtt*dxx_dzeta+dpade_dxx*d2xx_dttdzeta
             d2pade_dzeta2=d2pade_dxx2*dxx_dzeta*dxx_dzeta+dpade_dxx*d2xx_dzeta2

             !   From pade to qq
             coeff_qq=tt*phi_zeta_inv*phi_zeta_inv
             qq=coeff_qq*pade
             dqq_drs=coeff_qq*dpade_drs
             dqq_dtt=pade*phi_zeta_inv*phi_zeta_inv+coeff_qq*dpade_dtt
             dqq_dzeta=coeff_qq*(dpade_dzeta-two*pade*phi_logder)
             d2qq_drs2=coeff_qq*d2pade_drs2
             d2qq_drsdtt=phi_zeta_inv*phi_zeta_inv*(dpade_drs+tt*d2pade_drsdtt)
             d2qq_drsdzeta=coeff_qq*(d2pade_drsdzeta-two*dpade_drs*phi_logder)
             d2qq_dtt2=phi_zeta_inv*phi_zeta_inv*(two*dpade_dtt+tt*d2pade_dtt2)
             d2qq_dttdzeta=phi_zeta_inv*phi_zeta_inv*(dpade_dzeta-two*pade*phi_logder)+&
                  &     coeff_qq*(d2pade_dttdzeta-two*dpade_dtt*phi_logder)
             d2qq_dzeta2=coeff_qq*( d2pade_dzeta2-four*dpade_dzeta*phi_logder &
                  &                          +six*pade*phi_logder*phi_logder            &
                  &                          -two*pade*phi_zeta_inv*phipp_zeta)

             !   From qq to rr
             arg_rr=one+beta*gamma_inv*qq
             div_rr=one/arg_rr
             rr=gamma*log(arg_rr)
             drr_dqq=beta*div_rr
             drr_drs=drr_dqq*dqq_drs
             drr_dtt=drr_dqq*dqq_dtt
             drr_dzeta=drr_dqq*dqq_dzeta
             d2rr_dqq2=-div_rr**2*beta*beta*gamma_inv
             d2rr_drs2=d2rr_dqq2*dqq_drs*dqq_drs+drr_dqq*d2qq_drs2
             d2rr_drsdtt=d2rr_dqq2*dqq_drs*dqq_dtt+drr_dqq*d2qq_drsdtt
             d2rr_drsdzeta=d2rr_dqq2*dqq_drs*dqq_dzeta+drr_dqq*d2qq_drsdzeta
             d2rr_dtt2=d2rr_dqq2*dqq_dtt*dqq_dtt+drr_dqq*d2qq_dtt2
             d2rr_dttdzeta=d2rr_dqq2*dqq_dtt*dqq_dzeta+drr_dqq*d2qq_dttdzeta
             d2rr_dzeta2=d2rr_dqq2*dqq_dzeta*dqq_dzeta+drr_dqq*d2qq_dzeta2

             !   From rr to hh
             hh=phi3_zeta*rr
             dhh_drs=phi3_zeta*drr_drs
             dhh_dtt=phi3_zeta*drr_dtt
             dhh_dzeta=phi3_zeta*(drr_dzeta+three*rr*phi_logder)
             d2hh_drs2=phi3_zeta*d2rr_drs2
             d2hh_drsdtt=phi3_zeta*d2rr_drsdtt
             d2hh_drsdzeta=phi3_zeta*(d2rr_drsdzeta+three*drr_drs*phi_logder)
             d2hh_dtt2=phi3_zeta*d2rr_dtt2
             d2hh_dttdzeta=phi3_zeta*(d2rr_dttdzeta+three*drr_dtt*phi_logder)
             d2hh_dzeta2=phi3_zeta*(six*rr*phi_logder*phi_logder+&
                  &                           six*phi_logder*drr_dzeta+d2rr_dzeta2)  &
                  &               +three*phi_zeta*phi_zeta*rr*phipp_zeta


             !   The GGA correlation energy is added
             exci(ipts)=exci(ipts)+hh

             !   Change of variables : from (rs,zeta,tt) to (rhoup,rhodn,grrho)



             !   From hh to the derivative of the energy wrt the density
             drhohh_drho=hh-third*rs*dhh_drs-zeta*dhh_dzeta-seven*third*tt*dhh_dtt
             vxci(ipts,1)=vxci(ipts,1)+drhohh_drho

             !   From hh to the derivative of the energy wrt to the gradient of the
             !   density, divided by the gradient of the density
             !   (The v3.3 definition includes the division by the norm of the gradient)
             dvxcdgr(ipts,3)=rhotot*dtt_dg*dhh_dtt

             d2rhohh_drho2=rhotot_inv*&
                  &     (-two*ninth*rs*dhh_drs +seven*four*ninth*tt*dhh_dtt &
                  &      +ninth*rs*rs*d2hh_drs2+zeta*zeta*d2hh_dzeta2+(seven*third*tt)**2*d2hh_dtt2 &
                  &      +two*third*rs*zeta*d2hh_drsdzeta+two*seven*ninth*rs*tt*d2hh_drsdtt &
                  &      +two*seven*third*tt*zeta*d2hh_dttdzeta)
             d2rhohh_drhodg=dtt_dg*(-four*third*dhh_dtt-third*rs*d2hh_drsdtt &
                  &                           -zeta*d2hh_dttdzeta-seven*third*tt*d2hh_dtt2)

             !    Component 12 : first derivative with respect to the gradient
             !    of the density, div by the grad of the density
             dvxci(ipts,12)=dvxcdgr(ipts,3)
             !    Components 9, 10 and 11 : second derivatives with respect to the spin-density
             !    Note that there is already a contribution from LSDA
             dvxci(ipts,9)=dvxci(ipts,9)+d2rhohh_drho2+rhotot_inv*           &
                  &     ( d2hh_dzeta2*(one-two*zeta) &
                  &      -two*third*rs*d2hh_drsdzeta-14.0_dp*third*tt*d2hh_dttdzeta)
             dvxci(ipts,10)=dvxci(ipts,10)+d2rhohh_drho2-rhotot_inv*d2hh_dzeta2
             dvxci(ipts,11)=dvxci(ipts,11)+d2rhohh_drho2+rhotot_inv*           &
                  &     ( d2hh_dzeta2*(one+two*zeta) &
                  &      +two*third*rs*d2hh_drsdzeta+14.0_dp*third*tt*d2hh_dttdzeta)
             !    Components 13 and 14 : second derivatives with respect to spin density
             !    and gradient, divided by the gradient
             dvxci(ipts,13)=d2rhohh_drhodg+dtt_dg*d2hh_dttdzeta
             dvxci(ipts,14)=d2rhohh_drhodg-dtt_dg*d2hh_dttdzeta
             !    Component 15 : derivative of the (derivative wrt the gradient div by the grad),
             !    divided by the grad
             dvxci(ipts,15)=rhotot*d2hh_dtt2*dtt_dg*dtt_dg


             !  End condition of GGA

             ! Correlation has been added
             ! -----------------------------------------------------------------------------

             !  vxci(ipts,2)=vxci(ipts,1)
             dvxcdgr(ipts,2)=dvxcdgr(ipts,1)

          end do
       else if ((option==6) .or. (option==7)) then

          do ipts=1,npts

             rhotot=rhoarr(ipts)
             rhotmot=rhom1_3(ipts)
             rhotot_inv=rhotmot*rhotmot*rhotmot
             rhotmo6=sqrt(rhotmot)
             rhoto6=rhotot*rhotmot*rhotmot*rhotmo6
             ! -----------------------------------------------------------------------
             ! First take care of the exchange part of the functional

             exc=zero
             dvxcdgr(ipts,3)=zero
             !loop over the spin
             ispden=1
             rho   =rho_updn(ipts,ispden)
             rhomot=rho_updnm1_3(ipts,ispden)
             ex_lsd= - threefourth_divpi * sixpi2_1_3*rhomot*rhomot*rho
             !   Perdew-Burke-Ernzerhof GGA, exchange part
             rho_inv=rhomot*rhomot*rhomot
             coeffss=quarter*sixpi2m1_3*sixpi2m1_3*rho_inv*rho_inv*rhomot*rhomot
             ss=grho2_updn(ipts,ispden)*coeffss

             if (option==6) then
                divss=exp(-mu_divkappa*ss)
                dfxdss= mu*divss
                d2fxdss2=-mu*mu_divkappa*divss

                fx    = one+kappa*(one-divss)
                ex_gga= ex_lsd*fx
                dssdn=-eight*third*ss*rho_inv
                dfxdn  = dfxdss*dssdn
                vxci(ipts,ispden)=ex_lsd*(four_thirds*fx+rho*dfxdn)
                !   The new definition (v3.3) includes the division by the norm of the gradient
                dssdg =two*coeffss
                dfxdg=dfxdss*dssdg
                dvxcdgr(ipts,ispden)=ex_lsd*rho*dfxdg
                exc=exc+ex_gga*rho
             ! This is the Wu and Cohen modification
             else
                expss=exp(-ss)
                p1_wc=b_wc+(mu-b_wc)*(one-ss)*expss+two*c_wc*ss/(one+c_wc*ss*ss)
                p2_wc=d_wc*(ss-two)*expss+two*c_wc/(one+c_wc*ss*ss)-&
&                     four*c_wc*c_wc*ss*ss/((one+c_wc*ss*ss)*(one+c_wc*ss*ss))
                divss=one/(one+(b_wc*ss+d_wc*ss*expss+log(one+c_wc*ss*ss))/kappa)
                dfxdss=p1_wc*divss*divss
                d2fxdss2=p2_wc*divss*divss-two*divss*divss*divss*p1_wc*p1_wc/kappa

                fx    = one+kappa*(one-divss)
                ex_gga= ex_lsd*fx
                dssdn=-eight*third*ss*rho_inv
                dfxdn  = dfxdss*dssdn
                vxci(ipts,ispden)=ex_lsd*(four_thirds*fx+rho*dfxdn)
                !   The new definition (v3.3) includes the division by the norm of the gradient
                dssdg =two*coeffss
                dfxdg=dfxdss*dssdg
                dvxcdgr(ipts,ispden)=ex_lsd*rho*dfxdg
                exc=exc+ex_gga*rho
             end if

             !    Perdew-Burke-Ernzerhof GGA, exchange part
             !    Components 3 or 4
             dvxci(ipts,2+ispden)=dvxcdgr(ipts,ispden)
             !    Components 1 or 2
             d2ssdn2=-11.0_dp*third*dssdn*rho_inv
             d2fxdn2=d2fxdss2*dssdn**2+dfxdss*d2ssdn2
             dvxci(ipts,ispden)=third*rho_inv*vxci(ipts,ispden)+&
                  &     ex_lsd*(seven*third*dfxdn+rho*d2fxdn2)
             !    Components 5 or 6
             d2ssdndg=-eight*third*dssdg*rho_inv
             d2fxdndg=d2fxdss2*dssdn*dssdg+dfxdss*d2ssdndg
             dvxci(ipts,4+ispden)=ex_lsd*(four_thirds*dfxdg+rho*d2fxdndg)
             !    Components 7 or 8
             d2fxdg2=d2fxdss2*dssdg**2
             dvxci(ipts,6+ispden)=ex_lsd*rho*d2fxdg2
             !    For the time being, treat non-spin-polarized like spin-polarized
             dvxci(ipts,2)=dvxci(ipts,1)
             dvxci(ipts,4)=dvxci(ipts,3)
             dvxci(ipts,6)=dvxci(ipts,5)
             dvxci(ipts,8)=dvxci(ipts,7)

             !end of loop over the spin
             !  If non spin-polarized, treat spin down contribution now, similar to spin up
             exc=exc*2
             exci(ipts)=exc*rhotot_inv
             ! -----------------------------------------------------------------------------
             ! Then takes care of the LSD correlation part of the functional


             rs=rsfac*rhotmot
             sqr_rs=sq_rsfac*rhotmo6
             rsm1_2=sq_rsfac_inv*rhoto6

             !  Formulas A6-A8 of PW92LSD
             ec0_q0=-2.0_dp*ec0_aa*(1.0_dp+ec0_a1*rs)
             ec0_q1=2.0_dp*ec0_aa*(ec0_b1*sqr_rs+ec0_b2*rs+ec0_b3*rs*sqr_rs+ec0_b4*rs*rs)
             ec0_q1p=ec0_aa*(ec0_b1*rsm1_2+2._dp*ec0_b2+3._dp*ec0_b3*sqr_rs+4._dp*ec0_b4*rs)
             ec0_den=1.0_dp/(ec0_q1*ec0_q1+ec0_q1)
             !  ec0_log=log( 1.0_dp + 1.0_dp / ec0_q1 )
             ec0_log=-log( ec0_q1*ec0_q1*ec0_den )
             ecrs0=ec0_q0*ec0_log
             decrs0_drs= -2.0_dp*ec0_aa*ec0_a1*ec0_log - ec0_q0*ec0_q1p *ec0_den
             ec0_q1pp=0.5_dp*ec0_aa*(-ec0_b1*rsm1_2**3+3._dp*ec0_b3*rsm1_2+8._dp*ec0_b4)
             d2ecrs0_drs2= 4.0_dp*ec0_aa*ec0_a1*ec0_q1p*ec0_den            &
                  &                -ec0_q0*ec0_q1pp*ec0_den                        &
                  &                +ec0_q0*ec0_q1p**2*ec0_den**2*(2._dp*ec0_q1+1.0_dp)
             ec0_q1ppp = 0.75_dp*ec0_aa*(rsm1_2**5)*(ec0_b1-ec0_b3*rs)
             ec0_f1 = 1._dp/(ec0_q1*ec0_q1*(1._dp + ec0_q1))
             ec0_f2 = 1._dp/(ec0_q1*(1+ec0_q1))
             d3ecrs0_drs3 = 6._dp*ec0_q1p*ec0_f1*(-2._dp*ec0_aa*ec0_a1*ec0_q1p + &
                  &        ec0_q0*ec0_q1pp) - &
                  &        ec0_f2*(-6._dp*ec0_aa*ec0_a1*ec0_q1pp + ec0_q0*ec0_q1ppp + &
                  &        ec0_f2*(3._dp*ec0_q1p*(-2._dp*ec0_aa*ec0_a1*ec0_q1p + ec0_q0*ec0_q1pp) + &
                  &        ec0_f2*2._dp*ec0_q0*(ec0_q1p**3)*(1._dp + 3._dp*ec0_q1*(1._dp + ec0_q1))))

             mac_q0=-2.0_dp*mac_aa*(1.0_dp+mac_a1*rs)
             mac_q1=2.0_dp*mac_aa*(mac_b1*sqr_rs+mac_b2*rs+mac_b3*rs*sqr_rs+mac_b4*rs*rs)
             mac_q1p=mac_aa*(mac_b1*rsm1_2+2._dp*mac_b2+3._dp*mac_b3*sqr_rs+4._dp*mac_b4*rs)
             mac_den=1.0_dp/(mac_q1*mac_q1+mac_q1)
             mac_log=-log( mac_q1*mac_q1*mac_den )
             macrs=mac_q0*mac_log
             dmacrs_drs= -2.0_dp*mac_aa*mac_a1*mac_log - mac_q0*mac_q1p*mac_den

             ecrs=ecrs0
             decrs_drs=decrs0_drs
             decrs_dzeta=0.0_dp
             d2ecrs_drs2=d2ecrs0_drs2
             d2ecrs_dzeta2=alpha_zeta**2*(-macrs)
             d2ecrs_drsdzeta=zero
             zeta=0.0_dp


             !  Add LSD correlation functional to GGA exchange functional
             exci(ipts)=exci(ipts)+ecrs
             vxci(ipts,1)=vxci(ipts,1)+ecrs-rs*third*decrs_drs

             dvcrs_drs=third*(2._dp*decrs_drs-rs*d2ecrs_drs2)
             !   And d(vxc)/d(rho)=(-rs/(3*rho))*d(vxc)/d(rs)
             d2ecrs_drho2= -rs**4*(four_pi*third)*third*dvcrs_drs
             dvxci(ipts,9)=d2ecrs_drho2
             dvxci(ipts,10)=d2ecrs_drho2
             dvxci(ipts,11)=d2ecrs_drho2

             !  -----------------------------------------------------------------------------
             !  Eventually add the GGA correlation part of the PBE functional
             !  Note : the computation of the potential in the spin-unpolarized
             !  case could be optimized much further. Other optimizations are left to do.

             phi_zeta=1.0_dp
             phip_zeta=0.0_dp
             phi_zeta_inv=1.0_dp
             phi_logder=0.0_dp
             phi3_zeta=1.0_dp
             gamphi3inv=gamma_inv
             phipp_zeta=-two*ninth*alpha_zeta*alpha_zeta

             !   From ec to bb
             bb=ecrs*gamphi3inv
             dbb_drs=decrs_drs*gamphi3inv
             dbb_dzeta=gamphi3inv*(decrs_dzeta-three*ecrs*phi_logder)
             d2bb_drs2=d2ecrs_drs2*gamphi3inv
             d2bb_drsdzeta=gamphi3inv*(d2ecrs_drsdzeta-three*decrs_drs*phi_logder)
             d2bb_dzeta2=gamphi3inv*(d2ecrs_dzeta2-six*decrs_dzeta*phi_logder+&
                  &     12.0_dp*ecrs*phi_logder*phi_logder-three*ecrs*phi_zeta_inv*phipp_zeta)

             !   From bb to cc
             exp_pbe=exp(-bb)
             cc=one/(exp_pbe-one)
             dcc_dbb=cc*cc*exp_pbe
             dcc_drs=dcc_dbb*dbb_drs
             dcc_dzeta=dcc_dbb*dbb_dzeta
             d2cc_dbb2=cc*cc*exp_pbe*(two*cc*exp_pbe-one)
             d2cc_drs2=d2cc_dbb2*dbb_drs*dbb_drs+dcc_dbb*d2bb_drs2
             d2cc_drsdzeta=d2cc_dbb2*dbb_drs*dbb_dzeta+dcc_dbb*d2bb_drsdzeta
             d2cc_dzeta2=d2cc_dbb2*dbb_dzeta*dbb_dzeta+dcc_dbb*d2bb_dzeta2

             !   From cc to aa
             coeff_aa=beta*gamma_inv*phi_zeta_inv*phi_zeta_inv
             aa=coeff_aa*cc
             daa_drs=coeff_aa*dcc_drs
             daa_dzeta=-two*aa*phi_logder+coeff_aa*dcc_dzeta
             d2aa_drs2=coeff_aa*d2cc_drs2
             d2aa_drsdzeta=-two*daa_drs*phi_logder+coeff_aa*d2cc_drsdzeta
             d2aa_dzeta2=aa*(-two*phi_zeta_inv*phipp_zeta+six*phi_logder*phi_logder)+&
                  &     coeff_aa*(-four*dcc_dzeta*phi_logder+d2cc_dzeta2)

             !   Introduce tt : do not assume that the spin-dependent gradients are collinear
             grrho2=four*grho2_updn(ipts,1)
             dtt_dg=two*rhotot_inv*rhotot_inv*rhotmot*coeff_tt
             !   Note that tt is (the t variable of PBE divided by phi) squared
             tt=half*grrho2*dtt_dg

             !   Get xx from aa and tt
             xx=aa*tt
             dxx_drs=daa_drs*tt
             dxx_dzeta=daa_dzeta*tt
             dxx_dtt=aa
             d2xx_drs2=d2aa_drs2*tt
             d2xx_drsdzeta=d2aa_drsdzeta*tt
             d2xx_drsdtt=daa_drs
             d2xx_dttdzeta=daa_dzeta
             d2xx_dzeta2=d2aa_dzeta2*tt

             !   From xx to pade
             pade_den=one/(one+xx*(one+xx))
             pade=(one+xx)*pade_den
             dpade_dxx=-xx*(two+xx)*pade_den**2
             dpade_drs=dpade_dxx*dxx_drs
             dpade_dtt=dpade_dxx*dxx_dtt
             dpade_dzeta=dpade_dxx*dxx_dzeta
             d2pade_dxx2=two*(-one+xx*xx*(three+xx))*pade_den*pade_den*pade_den
             d2pade_drs2=d2pade_dxx2*dxx_drs*dxx_drs+dpade_dxx*d2xx_drs2
             d2pade_drsdtt=d2pade_dxx2*dxx_drs*dxx_dtt+dpade_dxx*d2xx_drsdtt
             d2pade_drsdzeta=d2pade_dxx2*dxx_drs*dxx_dzeta+dpade_dxx*d2xx_drsdzeta
             d2pade_dtt2=d2pade_dxx2*dxx_dtt*dxx_dtt
             d2pade_dttdzeta=d2pade_dxx2*dxx_dtt*dxx_dzeta+dpade_dxx*d2xx_dttdzeta
             d2pade_dzeta2=d2pade_dxx2*dxx_dzeta*dxx_dzeta+dpade_dxx*d2xx_dzeta2

             !   From pade to qq
             coeff_qq=tt*phi_zeta_inv*phi_zeta_inv
             qq=coeff_qq*pade
             dqq_drs=coeff_qq*dpade_drs
             dqq_dtt=pade*phi_zeta_inv*phi_zeta_inv+coeff_qq*dpade_dtt
             dqq_dzeta=coeff_qq*(dpade_dzeta-two*pade*phi_logder)
             d2qq_drs2=coeff_qq*d2pade_drs2
             d2qq_drsdtt=phi_zeta_inv*phi_zeta_inv*(dpade_drs+tt*d2pade_drsdtt)
             d2qq_drsdzeta=coeff_qq*(d2pade_drsdzeta-two*dpade_drs*phi_logder)
             d2qq_dtt2=phi_zeta_inv*phi_zeta_inv*(two*dpade_dtt+tt*d2pade_dtt2)
             d2qq_dttdzeta=phi_zeta_inv*phi_zeta_inv*(dpade_dzeta-two*pade*phi_logder)+&
                  &     coeff_qq*(d2pade_dttdzeta-two*dpade_dtt*phi_logder)
             d2qq_dzeta2=coeff_qq*( d2pade_dzeta2-four*dpade_dzeta*phi_logder &
                  &                          +six*pade*phi_logder*phi_logder            &
                  &                          -two*pade*phi_zeta_inv*phipp_zeta)

             !   From qq to rr
             arg_rr=one+beta*gamma_inv*qq
             div_rr=one/arg_rr
             rr=gamma*log(arg_rr)
             drr_dqq=beta*div_rr
             drr_drs=drr_dqq*dqq_drs
             drr_dtt=drr_dqq*dqq_dtt
             drr_dzeta=drr_dqq*dqq_dzeta
             d2rr_dqq2=-div_rr**2*beta*beta*gamma_inv
             d2rr_drs2=d2rr_dqq2*dqq_drs*dqq_drs+drr_dqq*d2qq_drs2
             d2rr_drsdtt=d2rr_dqq2*dqq_drs*dqq_dtt+drr_dqq*d2qq_drsdtt
             d2rr_drsdzeta=d2rr_dqq2*dqq_drs*dqq_dzeta+drr_dqq*d2qq_drsdzeta
             d2rr_dtt2=d2rr_dqq2*dqq_dtt*dqq_dtt+drr_dqq*d2qq_dtt2
             d2rr_dttdzeta=d2rr_dqq2*dqq_dtt*dqq_dzeta+drr_dqq*d2qq_dttdzeta
             d2rr_dzeta2=d2rr_dqq2*dqq_dzeta*dqq_dzeta+drr_dqq*d2qq_dzeta2

             !   From rr to hh
             hh=phi3_zeta*rr
             dhh_drs=phi3_zeta*drr_drs
             dhh_dtt=phi3_zeta*drr_dtt
             dhh_dzeta=phi3_zeta*(drr_dzeta+three*rr*phi_logder)
             d2hh_drs2=phi3_zeta*d2rr_drs2
             d2hh_drsdtt=phi3_zeta*d2rr_drsdtt
             d2hh_drsdzeta=phi3_zeta*(d2rr_drsdzeta+three*drr_drs*phi_logder)
             d2hh_dtt2=phi3_zeta*d2rr_dtt2
             d2hh_dttdzeta=phi3_zeta*(d2rr_dttdzeta+three*drr_dtt*phi_logder)
             d2hh_dzeta2=phi3_zeta*(six*rr*phi_logder*phi_logder+&
                  &                           six*phi_logder*drr_dzeta+d2rr_dzeta2)  &
                  &               +three*phi_zeta*phi_zeta*rr*phipp_zeta


             !   The GGA correlation energy is added
             exci(ipts)=exci(ipts)+hh

             !   Change of variables : from (rs,zeta,tt) to (rhoup,rhodn,grrho)

             !   From hh to the derivative of the energy wrt the density
             drhohh_drho=hh-third*rs*dhh_drs-zeta*dhh_dzeta-seven*third*tt*dhh_dtt
             vxci(ipts,1)=vxci(ipts,1)+drhohh_drho

             !   From hh to the derivative of the energy wrt to the gradient of the
             !   density, divided by the gradient of the density
             !   (The v3.3 definition includes the division by the norm of the gradient)
             dvxcdgr(ipts,3)=rhotot*dtt_dg*dhh_dtt

             d2rhohh_drho2=rhotot_inv*&
                  &     (-two*ninth*rs*dhh_drs +seven*four*ninth*tt*dhh_dtt &
                  &      +ninth*rs*rs*d2hh_drs2+zeta*zeta*d2hh_dzeta2+(seven*third*tt)**2*d2hh_dtt2 &
                  &      +two*third*rs*zeta*d2hh_drsdzeta+two*seven*ninth*rs*tt*d2hh_drsdtt &
                  &      +two*seven*third*tt*zeta*d2hh_dttdzeta)
             d2rhohh_drhodg=dtt_dg*(-four*third*dhh_dtt-third*rs*d2hh_drsdtt &
                  &                           -zeta*d2hh_dttdzeta-seven*third*tt*d2hh_dtt2)

             !    Component 12 : first derivative with respect to the gradient
             !    of the density, div by the grad of the density
             dvxci(ipts,12)=dvxcdgr(ipts,3)
             !    Components 9, 10 and 11 : second derivatives with respect to the spin-density
             !    Note that there is already a contribution from LSDA
             dvxci(ipts,9)=dvxci(ipts,9)+d2rhohh_drho2+rhotot_inv*           &
                  &     ( d2hh_dzeta2*(one-two*zeta) &
                  &      -two*third*rs*d2hh_drsdzeta-14.0_dp*third*tt*d2hh_dttdzeta)
             dvxci(ipts,10)=dvxci(ipts,10)+d2rhohh_drho2-rhotot_inv*d2hh_dzeta2
             dvxci(ipts,11)=dvxci(ipts,11)+d2rhohh_drho2+rhotot_inv*           &
                  &     ( d2hh_dzeta2*(one+two*zeta) &
                  &      +two*third*rs*d2hh_drsdzeta+14.0_dp*third*tt*d2hh_dttdzeta)
             !    Components 13 and 14 : second derivatives with respect to spin density
             !    and gradient, divided by the gradient
             dvxci(ipts,13)=d2rhohh_drhodg+dtt_dg*d2hh_dttdzeta
             dvxci(ipts,14)=d2rhohh_drhodg-dtt_dg*d2hh_dttdzeta
             !    Component 15 : derivative of the (derivative wrt the gradient div by the grad),
             !    divided by the grad
             dvxci(ipts,15)=rhotot*d2hh_dtt2*dtt_dg*dtt_dg


             !  End condition of GGA


             ! Correlation has been added
             ! -----------------------------------------------------------------------------

             !  vxci(ipts,2)=vxci(ipts,1)
             dvxcdgr(ipts,2)=dvxcdgr(ipts,1)
          end do

       else if (option==-1) then

          do ipts=1,npts

             rhotot=rhoarr(ipts)
             rhotmot=rhom1_3(ipts)
             rhotot_inv=rhotmot*rhotmot*rhotmot
             rhotmo6=sqrt(rhotmot)
             rhoto6=rhotot*rhotmot*rhotmot*rhotmo6
             ! -----------------------------------------------------------------------
             ! First take care of the exchange part of the functional

             exc=zero
             !loop over the spin
             ispden=1
             rho   =rho_updn(ipts,ispden)
             rhomot=rho_updnm1_3(ipts,ispden)
             ex_lsd= - threefourth_divpi * sixpi2_1_3*rhomot*rhomot*rho
             !   Perdew_Wang 91 LSD
             vxci(ipts,ispden)=four_thirds*ex_lsd
             exc=exc+ex_lsd*rho
             !    Perdew_Wang 91 LSD
             dvxci(ipts,2*ispden-1)=-four_thirds*third*&
                  &                   threefourth_divpi*sixpi2_1_3*rhomot*rhomot
             dvxci(ipts,2)=zero
             !    If non-spin-polarized, first component of dvxci is second
             !    derivative with respect to TOTAL density.
             dvxci(ipts,1)=dvxci(ipts,1)*half
             ! Compute the second derivative of vx
             ! vx^(2) = -2*vx^(1)/(3*rhotot)
             d2vxci(ipts) = -2._dp*dvxci(ipts,1)/(3._dp*rhotot)
             !end of loop over the spin
             !  If non spin-polarized, treat spin down contribution now, similar to spin up
             exc=exc*2
             exci(ipts)=exc*rhotot_inv
             ! -----------------------------------------------------------------------------
             ! Then takes care of the LSD correlation part of the functional

             ! Correlation has been added
             ! -----------------------------------------------------------------------------

          end do
       else if (option==-2) then

          do ipts=1,npts

             rhotot=rhoarr(ipts)
             rhotmot=rhom1_3(ipts)
             rhotot_inv=rhotmot*rhotmot*rhotmot
             rhotmo6=sqrt(rhotmot)
             rhoto6=rhotot*rhotmot*rhotmot*rhotmo6
             ! -----------------------------------------------------------------------
             ! First take care of the exchange part of the functional

             exc=zero
             dvxcdgr(ipts,3)=zero
             !loop over the spin
             ispden=1
             rho   =rho_updn(ipts,ispden)
             rhomot=rho_updnm1_3(ipts,ispden)
             ex_lsd= - threefourth_divpi * sixpi2_1_3*rhomot*rhomot*rho
             !   Perdew-Burke-Ernzerhof GGA, exchange part
             rho_inv=rhomot*rhomot*rhomot
             coeffss=quarter*sixpi2m1_3*sixpi2m1_3*rho_inv*rho_inv*rhomot*rhomot
             ss=grho2_updn(ipts,ispden)*coeffss
             divss=one/(one+mu_divkappa*ss)
             dfxdss= mu*divss*divss
             d2fxdss2=-mu*two*mu_divkappa*divss*divss*divss
             fx    = one+kappa*(one-divss)
             ex_gga= ex_lsd*fx
             dssdn=-eight*third*ss*rho_inv
             dfxdn  = dfxdss*dssdn
             vxci(ipts,ispden)=ex_lsd*(four_thirds*fx+rho*dfxdn)
             !   The new definition (v3.3) includes the division by the norm of the gradient
             dssdg =two*coeffss
             dfxdg=dfxdss*dssdg
             dvxcdgr(ipts,ispden)=ex_lsd*rho*dfxdg
             exc=exc+ex_gga*rho

             !    Perdew-Burke-Ernzerhof GGA, exchange part
             !    Components 3 or 4
             dvxci(ipts,2+ispden)=dvxcdgr(ipts,ispden)
             !    Components 1 or 2
             d2ssdn2=-11.0_dp*third*dssdn*rho_inv
             d2fxdn2=d2fxdss2*dssdn**2+dfxdss*d2ssdn2
             dvxci(ipts,ispden)=third*rho_inv*vxci(ipts,ispden)+&
                  &     ex_lsd*(seven*third*dfxdn+rho*d2fxdn2)
             !    Components 5 or 6
             d2ssdndg=-eight*third*dssdg*rho_inv
             d2fxdndg=d2fxdss2*dssdn*dssdg+dfxdss*d2ssdndg
             dvxci(ipts,4+ispden)=ex_lsd*(four_thirds*dfxdg+rho*d2fxdndg)
             !    Components 7 or 8
             d2fxdg2=d2fxdss2*dssdg**2
             dvxci(ipts,6+ispden)=ex_lsd*rho*d2fxdg2
             !    For the time being, treat non-spin-polarized like spin-polarized
             dvxci(ipts,2)=dvxci(ipts,1)
             dvxci(ipts,4)=dvxci(ipts,3)
             dvxci(ipts,6)=dvxci(ipts,5)
             dvxci(ipts,8)=dvxci(ipts,7)

             !end of loop over the spin
             !  If non spin-polarized, treat spin down contribution now, similar to spin up
             exc=exc*2
             exci(ipts)=exc*rhotot_inv
             ! -----------------------------------------------------------------------------
             ! Then takes care of the LSD correlation part of the functional

             ! Correlation has been added
             ! -----------------------------------------------------------------------------

             !  vxci(ipts,2)=vxci(ipts,1)
             dvxcdgr(ipts,2)=dvxcdgr(ipts,1)

          end do
       else if(option==1) then

          do ipts=1,npts


             rhotot=rhoarr(ipts)
             rhotmot=rhom1_3(ipts)
             rhotot_inv=rhotmot*rhotmot*rhotmot
             rhotmo6=sqrt(rhotmot)
             rhoto6=rhotot*rhotmot*rhotmot*rhotmo6
             ! -----------------------------------------------------------------------
             ! First take care of the exchange part of the functional

             exc=zero
             !loop over the spin
             ispden=1
             rho   =rho_updn(ipts,ispden)
             rhomot=rho_updnm1_3(ipts,ispden)
             ex_lsd= - threefourth_divpi * sixpi2_1_3*rhomot*rhomot*rho
             !   Perdew_Wang 91 LSD
             vxci(ipts,ispden)=four_thirds*ex_lsd
             exc=exc+ex_lsd*rho

             !    Perdew_Wang 91 LSD
             dvxci(ipts,2*ispden-1)=-four_thirds*third*&
                  &                   threefourth_divpi*sixpi2_1_3*rhomot*rhomot
             dvxci(ipts,2)=zero
             !    If non-spin-polarized, first component of dvxci is second
             !    derivative with respect to TOTAL density.
             dvxci(ipts,1)=dvxci(ipts,1)*half
             ! Compute the second derivative of vx
             ! vx^(2) = -2*vx^(1)/(3*rhotot)
             d2vxci(ipts) = -2._dp*dvxci(ipts,1)/(3._dp*rhotot)
             !end of loop over the spin
             !  If non spin-polarized, treat spin down contribution now, similar to spin up
             exc=exc*2
             exci(ipts)=exc*rhotot_inv
             ! -----------------------------------------------------------------------------
             ! Then takes care of the LSD correlation part of the functional

             rs=rsfac*rhotmot
             sqr_rs=sq_rsfac*rhotmo6
             rsm1_2=sq_rsfac_inv*rhoto6

             !  Formulas A6-A8 of PW92LSD
             ec0_q0=-2.0_dp*ec0_aa*(1.0_dp+ec0_a1*rs)
             ec0_q1=2.0_dp*ec0_aa*(ec0_b1*sqr_rs+ec0_b2*rs+ec0_b3*rs*sqr_rs+ec0_b4*rs*rs)
             ec0_q1p=ec0_aa*(ec0_b1*rsm1_2+2._dp*ec0_b2+3._dp*ec0_b3*sqr_rs+4._dp*ec0_b4*rs)
             ec0_den=1.0_dp/(ec0_q1*ec0_q1+ec0_q1)
             !  ec0_log=log( 1.0_dp + 1.0_dp / ec0_q1 )
             ec0_log=-log( ec0_q1*ec0_q1*ec0_den )
             ecrs0=ec0_q0*ec0_log
             decrs0_drs= -2.0_dp*ec0_aa*ec0_a1*ec0_log - ec0_q0*ec0_q1p *ec0_den
             ec0_q1pp=0.5_dp*ec0_aa*(-ec0_b1*rsm1_2**3+3._dp*ec0_b3*rsm1_2+8._dp*ec0_b4)
             d2ecrs0_drs2= 4.0_dp*ec0_aa*ec0_a1*ec0_q1p*ec0_den            &
                  &                -ec0_q0*ec0_q1pp*ec0_den                        &
                  &                +ec0_q0*ec0_q1p**2*ec0_den**2*(2._dp*ec0_q1+1.0_dp)
             ec0_q1ppp = 0.75_dp*ec0_aa*(rsm1_2**5)*(ec0_b1-ec0_b3*rs)
             ec0_f1 = 1._dp/(ec0_q1*ec0_q1*(1._dp + ec0_q1))
             ec0_f2 = 1._dp/(ec0_q1*(1+ec0_q1))
             d3ecrs0_drs3 = 6._dp*ec0_q1p*ec0_f1*(-2._dp*ec0_aa*ec0_a1*ec0_q1p + &
                  &        ec0_q0*ec0_q1pp) - &
                  &        ec0_f2*(-6._dp*ec0_aa*ec0_a1*ec0_q1pp + ec0_q0*ec0_q1ppp + &
                  &        ec0_f2*(3._dp*ec0_q1p*(-2._dp*ec0_aa*ec0_a1*ec0_q1p + ec0_q0*ec0_q1pp) + &
                  &        ec0_f2*2._dp*ec0_q0*(ec0_q1p**3)*(1._dp + 3._dp*ec0_q1*(1._dp + ec0_q1))))

             mac_q0=-2.0_dp*mac_aa*(1.0_dp+mac_a1*rs)
             mac_q1=2.0_dp*mac_aa*(mac_b1*sqr_rs+mac_b2*rs+mac_b3*rs*sqr_rs+mac_b4*rs*rs)
             mac_q1p=mac_aa*(mac_b1*rsm1_2+2._dp*mac_b2+3._dp*mac_b3*sqr_rs+4._dp*mac_b4*rs)
             mac_den=1.0_dp/(mac_q1*mac_q1+mac_q1)
             mac_log=-log( mac_q1*mac_q1*mac_den )
             macrs=mac_q0*mac_log
             dmacrs_drs= -2.0_dp*mac_aa*mac_a1*mac_log - mac_q0*mac_q1p*mac_den

             ecrs=ecrs0
             decrs_drs=decrs0_drs
             decrs_dzeta=0.0_dp
             d2ecrs_drs2=d2ecrs0_drs2
             d2ecrs_dzeta2=alpha_zeta**2*(-macrs)
             d2ecrs_drsdzeta=zero
             zeta=0.0_dp


             !  Add LSD correlation functional to GGA exchange functional
             exci(ipts)=exci(ipts)+ecrs
             vxci(ipts,1)=vxci(ipts,1)+ecrs-rs*third*decrs_drs

             dvcrs_drs=third*(2._dp*decrs_drs-rs*d2ecrs_drs2)
             !   And d(vxc)/d(rho)=(-rs/(3*rho))*d(vxc)/d(rs)
             d2ecrs_drho2= -rs**4*(four_pi*third)*third*dvcrs_drs
             dvxci(ipts,1)=dvxci(ipts,1)+d2ecrs_drho2
             d2vcrs_drs2 = third*(d2ecrs_drs2 - rs*d3ecrs0_drs3)
             drs_dn = -1._dp*four_pi*ninth*rs**4
             d2rs_dn2 = 64._dp*pi*pi*(rs**7)/81._dp
             d2vxci(ipts) = d2vxci(ipts) + d2vcrs_drs2*drs_dn*drs_dn + &
                  &        dvcrs_drs*d2rs_dn2

             !  -----------------------------------------------------------------------------
             !  Eventually add the GGA correlation part of the PBE functional
             !  Note : the computation of the potential in the spin-unpolarized
             !  case could be optimized much further. Other optimizations are left to do.


             ! Correlation has been added
             ! -----------------------------------------------------------------------------
          end do
       else if (option==3) then

          do ipts=1,npts

             rhotot=rhoarr(ipts)
             rhotmot=rhom1_3(ipts)
             rhotot_inv=rhotmot*rhotmot*rhotmot
             rhotmo6=sqrt(rhotmot)
             rhoto6=rhotot*rhotmot*rhotmot*rhotmo6
             ! -----------------------------------------------------------------------
             ! First take care of the exchange part of the functional

             exc=zero
             !loop over the spin
             ispden=1
             rho   =rho_updn(ipts,ispden)
             rhomot=rho_updnm1_3(ipts,ispden)
             ex_lsd= - threefourth_divpi * sixpi2_1_3*rhomot*rhomot*rho
             !   Perdew_Wang 91 LSD
             vxci(ipts,ispden)=four_thirds*ex_lsd
             exc=exc+ex_lsd*rho


             !    Perdew_Wang 91 LSD
             dvxci(ipts,2*ispden-1)=-four_thirds*third*&
                  &                   threefourth_divpi*sixpi2_1_3*rhomot*rhomot
             dvxci(ipts,2)=zero
             !    If non-spin-polarized, first component of dvxci is second
             !    derivative with respect to TOTAL density.
             dvxci(ipts,1)=dvxci(ipts,1)*half
             ! Compute the second derivative of vx
             ! vx^(2) = -2*vx^(1)/(3*rhotot)
             d2vxci(ipts) = -2._dp*dvxci(ipts,1)/(3._dp*rhotot)
             !end of loop over the spin
             !  If non spin-polarized, treat spin down contribution now, similar to spin up
             exc=exc*2
             exci(ipts)=exc*rhotot_inv
             ! -----------------------------------------------------------------------------
             ! Then takes care of the LSD correlation part of the functional

             rs=rsfac*rhotmot
             sqr_rs=sq_rsfac*rhotmo6
             rsm1_2=sq_rsfac_inv*rhoto6

             !  Formulas A6-A8 of PW92LSD
             ec0_q0=-2.0_dp*ec0_aa*(1.0_dp+ec0_a1*rs)
             sqr_sqr_rs=max(1.e-15_dp,sqrt(sqr_rs))
             ec0_q1=2.0_dp*ec0_aa*(ec0_b1*sqr_rs+ec0_b2*rs+ec0_b3*rs*sqr_rs+ec0_b4*rs*rs/sqr_sqr_rs)
             ec0_q1p=ec0_aa*(ec0_b1*rsm1_2+2._dp*ec0_b2+3._dp*ec0_b3*sqr_rs+3.5_dp*ec0_b4*rs/sqr_sqr_rs)
             ec0_den=1.0_dp/(ec0_q1*ec0_q1+ec0_q1)
             !  ec0_log=log( 1.0_dp + 1.0_dp / ec0_q1 )
             ec0_log=-log( ec0_q1*ec0_q1*ec0_den )
             ecrs0=ec0_q0*ec0_log
             decrs0_drs= -2.0_dp*ec0_aa*ec0_a1*ec0_log - ec0_q0*ec0_q1p *ec0_den
             ec0_q1pp=0.5_dp*ec0_aa*(-ec0_b1*rsm1_2**3+3._dp*ec0_b3*rsm1_2+8._dp*ec0_b4)
             d2ecrs0_drs2= 4.0_dp*ec0_aa*ec0_a1*ec0_q1p*ec0_den            &
                  &                -ec0_q0*ec0_q1pp*ec0_den                        &
                  &                +ec0_q0*ec0_q1p**2*ec0_den**2*(2._dp*ec0_q1+1.0_dp)
             ec0_q1ppp = 0.75_dp*ec0_aa*(rsm1_2**5)*(ec0_b1-ec0_b3*rs)
             ec0_f1 = 1._dp/(ec0_q1*ec0_q1*(1._dp + ec0_q1))
             ec0_f2 = 1._dp/(ec0_q1*(1+ec0_q1))
             d3ecrs0_drs3 = 6._dp*ec0_q1p*ec0_f1*(-2._dp*ec0_aa*ec0_a1*ec0_q1p + &
                  &        ec0_q0*ec0_q1pp) - &
                  &        ec0_f2*(-6._dp*ec0_aa*ec0_a1*ec0_q1pp + ec0_q0*ec0_q1ppp + &
                  &        ec0_f2*(3._dp*ec0_q1p*(-2._dp*ec0_aa*ec0_a1*ec0_q1p + ec0_q0*ec0_q1pp) + &
                  &        ec0_f2*2._dp*ec0_q0*(ec0_q1p**3)*(1._dp + 3._dp*ec0_q1*(1._dp + ec0_q1))))

             mac_q0=-2.0_dp*mac_aa*(1.0_dp+mac_a1*rs)
             mac_q1=2.0_dp*mac_aa*(mac_b1*sqr_rs+mac_b2*rs+mac_b3*rs*sqr_rs+mac_b4*rs*rs)
             mac_q1p=mac_aa*(mac_b1*rsm1_2+2._dp*mac_b2+3._dp*mac_b3*sqr_rs+4._dp*mac_b4*rs)
             mac_den=1.0_dp/(mac_q1*mac_q1+mac_q1)
             mac_log=-log( mac_q1*mac_q1*mac_den )
             macrs=mac_q0*mac_log
             dmacrs_drs= -2.0_dp*mac_aa*mac_a1*mac_log - mac_q0*mac_q1p*mac_den

             ecrs=ecrs0
             decrs_drs=decrs0_drs
             decrs_dzeta=0.0_dp
             d2ecrs_drs2=d2ecrs0_drs2
             d2ecrs_dzeta2=alpha_zeta**2*(-macrs)
             d2ecrs_drsdzeta=zero
             zeta=0.0_dp


             !  Add LSD correlation functional to GGA exchange functional
             exci(ipts)=exci(ipts)+ecrs
             vxci(ipts,1)=vxci(ipts,1)+ecrs-rs*third*decrs_drs

             dvcrs_drs=third*(2._dp*decrs_drs-rs*d2ecrs_drs2)
             !   And d(vxc)/d(rho)=(-rs/(3*rho))*d(vxc)/d(rs)
             d2ecrs_drho2= -rs**4*(four_pi*third)*third*dvcrs_drs

             dvxci(ipts,1)=dvxci(ipts,1)+d2ecrs_drho2
             d2vcrs_drs2 = third*(d2ecrs_drs2 - rs*d3ecrs0_drs3)
             drs_dn = -1._dp*four_pi*ninth*rs**4
             d2rs_dn2 = 64._dp*pi*pi*(rs**7)/81._dp
             d2vxci(ipts) = d2vxci(ipts) + d2vcrs_drs2*drs_dn*drs_dn + &
                  &        dvcrs_drs*d2rs_dn2


          end do

       end if

    else if(order==-2) then
       if(option==2 .or. option==5) then
          do ipts=1,npts

             rhotot=rhoarr(ipts)
             rhotmot=rhom1_3(ipts)
             rhotot_inv=rhotmot*rhotmot*rhotmot
             rhotmo6=sqrt(rhotmot)
             rhoto6=rhotot*rhotmot*rhotmot*rhotmo6
             ! -----------------------------------------------------------------------
             ! First take care of the exchange part of the functional

             exc=zero
             dvxcdgr(ipts,3)=zero
             !loop over the spin
             ispden=1
             rho   =rho_updn(ipts,ispden)
             rhomot=rho_updnm1_3(ipts,ispden)
             ex_lsd= - threefourth_divpi * sixpi2_1_3*rhomot*rhomot*rho
             !   Perdew-Burke-Ernzerhof GGA, exchange part
             rho_inv=rhomot*rhomot*rhomot
             coeffss=quarter*sixpi2m1_3*sixpi2m1_3*rho_inv*rho_inv*rhomot*rhomot
             ss=grho2_updn(ipts,ispden)*coeffss
             divss=one/(one+mu_divkappa*ss)
             dfxdss= mu*divss*divss
             d2fxdss2=-mu*two*mu_divkappa*divss*divss*divss
             fx    = one+kappa*(one-divss)
             ex_gga= ex_lsd*fx
             dssdn=-eight*third*ss*rho_inv
             dfxdn  = dfxdss*dssdn
             vxci(ipts,ispden)=ex_lsd*(four_thirds*fx+rho*dfxdn)
             !   The new definition (v3.3) includes the division by the norm of the gradient
             dssdg =two*coeffss
             dfxdg=dfxdss*dssdg
             dvxcdgr(ipts,ispden)=ex_lsd*rho*dfxdg
             exc=exc+ex_gga*rho

             !    Perdew-Burke-Ernzerhof GGA, exchange part
             !    Components 3 or 4
             dvxci(ipts,2+ispden)=dvxcdgr(ipts,ispden)
             !    Components 1 or 2
             d2ssdn2=-11.0_dp*third*dssdn*rho_inv
             d2fxdn2=d2fxdss2*dssdn**2+dfxdss*d2ssdn2
             dvxci(ipts,ispden)=third*rho_inv*vxci(ipts,ispden)+&
                  &     ex_lsd*(seven*third*dfxdn+rho*d2fxdn2)
             !    Components 5 or 6
             d2ssdndg=-eight*third*dssdg*rho_inv
             d2fxdndg=d2fxdss2*dssdn*dssdg+dfxdss*d2ssdndg
             dvxci(ipts,4+ispden)=ex_lsd*(four_thirds*dfxdg+rho*d2fxdndg)
             !    Components 7 or 8
             d2fxdg2=d2fxdss2*dssdg**2
             dvxci(ipts,6+ispden)=ex_lsd*rho*d2fxdg2
             !    For the time being, treat non-spin-polarized like spin-polarized
             dvxci(ipts,2)=dvxci(ipts,1)
             dvxci(ipts,4)=dvxci(ipts,3)
             dvxci(ipts,6)=dvxci(ipts,5)
             dvxci(ipts,8)=dvxci(ipts,7)
             !end of loop over the spin
             !  If non spin-polarized, treat spin down contribution now, similar to spin up
             exc=exc*2
             exci(ipts)=exc*rhotot_inv
             ! -----------------------------------------------------------------------------
             ! Then takes care of the LSD correlation part of the functional


             rs=rsfac*rhotmot
             sqr_rs=sq_rsfac*rhotmo6
             rsm1_2=sq_rsfac_inv*rhoto6

             !  Formulas A6-A8 of PW92LSD
             ec0_q0=-2.0_dp*ec0_aa*(1.0_dp+ec0_a1*rs)
             ec0_q1=2.0_dp*ec0_aa*(ec0_b1*sqr_rs+ec0_b2*rs+ec0_b3*rs*sqr_rs+ec0_b4*rs*rs)
             ec0_q1p=ec0_aa*(ec0_b1*rsm1_2+2._dp*ec0_b2+3._dp*ec0_b3*sqr_rs+4._dp*ec0_b4*rs)
             ec0_den=1.0_dp/(ec0_q1*ec0_q1+ec0_q1)
             !  ec0_log=log( 1.0_dp + 1.0_dp / ec0_q1 )
             ec0_log=-log( ec0_q1*ec0_q1*ec0_den )
             ecrs0=ec0_q0*ec0_log
             decrs0_drs= -2.0_dp*ec0_aa*ec0_a1*ec0_log - ec0_q0*ec0_q1p *ec0_den
             ec0_q1pp=0.5_dp*ec0_aa*(-ec0_b1*rsm1_2**3+3._dp*ec0_b3*rsm1_2+8._dp*ec0_b4)
             d2ecrs0_drs2= 4.0_dp*ec0_aa*ec0_a1*ec0_q1p*ec0_den            &
                  &                -ec0_q0*ec0_q1pp*ec0_den                        &
                  &                +ec0_q0*ec0_q1p**2*ec0_den**2*(2._dp*ec0_q1+1.0_dp)

             mac_q0=-2.0_dp*mac_aa*(1.0_dp+mac_a1*rs)
             mac_q1=2.0_dp*mac_aa*(mac_b1*sqr_rs+mac_b2*rs+mac_b3*rs*sqr_rs+mac_b4*rs*rs)
             mac_q1p=mac_aa*(mac_b1*rsm1_2+2._dp*mac_b2+3._dp*mac_b3*sqr_rs+4._dp*mac_b4*rs)
             mac_den=1.0_dp/(mac_q1*mac_q1+mac_q1)
             mac_log=-log( mac_q1*mac_q1*mac_den )
             macrs=mac_q0*mac_log
             dmacrs_drs= -2.0_dp*mac_aa*mac_a1*mac_log - mac_q0*mac_q1p*mac_den

             ecrs=ecrs0
             decrs_drs=decrs0_drs
             decrs_dzeta=0.0_dp
             d2ecrs_drs2=d2ecrs0_drs2
             d2ecrs_dzeta2=alpha_zeta**2*(-macrs)
             d2ecrs_drsdzeta=zero
             zeta=0.0_dp


             !  Add LSD correlation functional to GGA exchange functional
             exci(ipts)=exci(ipts)+ecrs
             vxci(ipts,1)=vxci(ipts,1)+ecrs-rs*third*decrs_drs

             dvcrs_drs=third*(2._dp*decrs_drs-rs*d2ecrs_drs2)
             !   And d(vxc)/d(rho)=(-rs/(3*rho))*d(vxc)/d(rs)
             d2ecrs_drho2= -rs**4*(four_pi*third)*third*dvcrs_drs
             dvxci(ipts,9)=d2ecrs_drho2
             dvxci(ipts,10)=d2ecrs_drho2
             dvxci(ipts,11)=d2ecrs_drho2

             !  -----------------------------------------------------------------------------
             !  Eventually add the GGA correlation part of the PBE functional
             !  Note : the computation of the potential in the spin-unpolarized
             !  case could be optimized much further. Other optimizations are left to do.

             phi_zeta=1.0_dp
             phip_zeta=0.0_dp
             phi_zeta_inv=1.0_dp
             phi_logder=0.0_dp
             phi3_zeta=1.0_dp
             gamphi3inv=gamma_inv
             phipp_zeta=-two*ninth*alpha_zeta*alpha_zeta

             !   From ec to bb
             bb=ecrs*gamphi3inv
             dbb_drs=decrs_drs*gamphi3inv
             dbb_dzeta=gamphi3inv*(decrs_dzeta-three*ecrs*phi_logder)
             d2bb_drs2=d2ecrs_drs2*gamphi3inv
             d2bb_drsdzeta=gamphi3inv*(d2ecrs_drsdzeta-three*decrs_drs*phi_logder)
             d2bb_dzeta2=gamphi3inv*(d2ecrs_dzeta2-six*decrs_dzeta*phi_logder+&
                  &     12.0_dp*ecrs*phi_logder*phi_logder-three*ecrs*phi_zeta_inv*phipp_zeta)

             !   From bb to cc
             exp_pbe=exp(-bb)
             cc=one/(exp_pbe-one)
             dcc_dbb=cc*cc*exp_pbe
             dcc_drs=dcc_dbb*dbb_drs
             dcc_dzeta=dcc_dbb*dbb_dzeta
             d2cc_dbb2=cc*cc*exp_pbe*(two*cc*exp_pbe-one)
             d2cc_drs2=d2cc_dbb2*dbb_drs*dbb_drs+dcc_dbb*d2bb_drs2
             d2cc_drsdzeta=d2cc_dbb2*dbb_drs*dbb_dzeta+dcc_dbb*d2bb_drsdzeta
             d2cc_dzeta2=d2cc_dbb2*dbb_dzeta*dbb_dzeta+dcc_dbb*d2bb_dzeta2

             !   From cc to aa
             coeff_aa=beta*gamma_inv*phi_zeta_inv*phi_zeta_inv
             aa=coeff_aa*cc
             daa_drs=coeff_aa*dcc_drs
             daa_dzeta=-two*aa*phi_logder+coeff_aa*dcc_dzeta
             d2aa_drs2=coeff_aa*d2cc_drs2
             d2aa_drsdzeta=-two*daa_drs*phi_logder+coeff_aa*d2cc_drsdzeta
             d2aa_dzeta2=aa*(-two*phi_zeta_inv*phipp_zeta+six*phi_logder*phi_logder)+&
                  &     coeff_aa*(-four*dcc_dzeta*phi_logder+d2cc_dzeta2)

             !   Introduce tt : do not assume that the spin-dependent gradients are collinear
             grrho2=four*grho2_updn(ipts,1)
             dtt_dg=two*rhotot_inv*rhotot_inv*rhotmot*coeff_tt
             !   Note that tt is (the t variable of PBE divided by phi) squared
             tt=half*grrho2*dtt_dg

             !   Get xx from aa and tt
             xx=aa*tt
             dxx_drs=daa_drs*tt
             dxx_dzeta=daa_dzeta*tt
             dxx_dtt=aa
             d2xx_drs2=d2aa_drs2*tt
             d2xx_drsdzeta=d2aa_drsdzeta*tt
             d2xx_drsdtt=daa_drs
             d2xx_dttdzeta=daa_dzeta
             d2xx_dzeta2=d2aa_dzeta2*tt

             !   From xx to pade
             pade_den=one/(one+xx*(one+xx))
             pade=(one+xx)*pade_den
             dpade_dxx=-xx*(two+xx)*pade_den**2
             dpade_drs=dpade_dxx*dxx_drs
             dpade_dtt=dpade_dxx*dxx_dtt
             dpade_dzeta=dpade_dxx*dxx_dzeta
             d2pade_dxx2=two*(-one+xx*xx*(three+xx))*pade_den*pade_den*pade_den
             d2pade_drs2=d2pade_dxx2*dxx_drs*dxx_drs+dpade_dxx*d2xx_drs2
             d2pade_drsdtt=d2pade_dxx2*dxx_drs*dxx_dtt+dpade_dxx*d2xx_drsdtt
             d2pade_drsdzeta=d2pade_dxx2*dxx_drs*dxx_dzeta+dpade_dxx*d2xx_drsdzeta
             d2pade_dtt2=d2pade_dxx2*dxx_dtt*dxx_dtt
             d2pade_dttdzeta=d2pade_dxx2*dxx_dtt*dxx_dzeta+dpade_dxx*d2xx_dttdzeta
             d2pade_dzeta2=d2pade_dxx2*dxx_dzeta*dxx_dzeta+dpade_dxx*d2xx_dzeta2

             !   From pade to qq
             coeff_qq=tt*phi_zeta_inv*phi_zeta_inv
             qq=coeff_qq*pade
             dqq_drs=coeff_qq*dpade_drs
             dqq_dtt=pade*phi_zeta_inv*phi_zeta_inv+coeff_qq*dpade_dtt
             dqq_dzeta=coeff_qq*(dpade_dzeta-two*pade*phi_logder)
             d2qq_drs2=coeff_qq*d2pade_drs2
             d2qq_drsdtt=phi_zeta_inv*phi_zeta_inv*(dpade_drs+tt*d2pade_drsdtt)
             d2qq_drsdzeta=coeff_qq*(d2pade_drsdzeta-two*dpade_drs*phi_logder)
             d2qq_dtt2=phi_zeta_inv*phi_zeta_inv*(two*dpade_dtt+tt*d2pade_dtt2)
             d2qq_dttdzeta=phi_zeta_inv*phi_zeta_inv*(dpade_dzeta-two*pade*phi_logder)+&
                  &     coeff_qq*(d2pade_dttdzeta-two*dpade_dtt*phi_logder)
             d2qq_dzeta2=coeff_qq*( d2pade_dzeta2-four*dpade_dzeta*phi_logder &
                  &                          +six*pade*phi_logder*phi_logder            &
                  &                          -two*pade*phi_zeta_inv*phipp_zeta)

             !   From qq to rr
             arg_rr=one+beta*gamma_inv*qq
             div_rr=one/arg_rr
             rr=gamma*log(arg_rr)
             drr_dqq=beta*div_rr
             drr_drs=drr_dqq*dqq_drs
             drr_dtt=drr_dqq*dqq_dtt
             drr_dzeta=drr_dqq*dqq_dzeta
             d2rr_dqq2=-div_rr**2*beta*beta*gamma_inv
             d2rr_drs2=d2rr_dqq2*dqq_drs*dqq_drs+drr_dqq*d2qq_drs2
             d2rr_drsdtt=d2rr_dqq2*dqq_drs*dqq_dtt+drr_dqq*d2qq_drsdtt
             d2rr_drsdzeta=d2rr_dqq2*dqq_drs*dqq_dzeta+drr_dqq*d2qq_drsdzeta
             d2rr_dtt2=d2rr_dqq2*dqq_dtt*dqq_dtt+drr_dqq*d2qq_dtt2
             d2rr_dttdzeta=d2rr_dqq2*dqq_dtt*dqq_dzeta+drr_dqq*d2qq_dttdzeta
             d2rr_dzeta2=d2rr_dqq2*dqq_dzeta*dqq_dzeta+drr_dqq*d2qq_dzeta2

             !   From rr to hh
             hh=phi3_zeta*rr
             dhh_drs=phi3_zeta*drr_drs
             dhh_dtt=phi3_zeta*drr_dtt
             dhh_dzeta=phi3_zeta*(drr_dzeta+three*rr*phi_logder)
             d2hh_drs2=phi3_zeta*d2rr_drs2
             d2hh_drsdtt=phi3_zeta*d2rr_drsdtt
             d2hh_drsdzeta=phi3_zeta*(d2rr_drsdzeta+three*drr_drs*phi_logder)
             d2hh_dtt2=phi3_zeta*d2rr_dtt2
             d2hh_dttdzeta=phi3_zeta*(d2rr_dttdzeta+three*drr_dtt*phi_logder)
             d2hh_dzeta2=phi3_zeta*(six*rr*phi_logder*phi_logder+&
                  &                           six*phi_logder*drr_dzeta+d2rr_dzeta2)  &
                  &               +three*phi_zeta*phi_zeta*rr*phipp_zeta

             !   The GGA correlation energy is added
             exci(ipts)=exci(ipts)+hh

             !   Change of variables : from (rs,zeta,tt) to (rhoup,rhodn,grrho)

             !   From hh to the derivative of the energy wrt the density
             drhohh_drho=hh-third*rs*dhh_drs-zeta*dhh_dzeta-seven*third*tt*dhh_dtt
             vxci(ipts,1)=vxci(ipts,1)+drhohh_drho

             !   From hh to the derivative of the energy wrt to the gradient of the
             !   density, divided by the gradient of the density
             !   (The v3.3 definition includes the division by the norm of the gradient)
             dvxcdgr(ipts,3)=rhotot*dtt_dg*dhh_dtt

             d2rhohh_drho2=rhotot_inv*&
                  &     (-two*ninth*rs*dhh_drs +seven*four*ninth*tt*dhh_dtt &
                  &      +ninth*rs*rs*d2hh_drs2+zeta*zeta*d2hh_dzeta2+(seven*third*tt)**2*d2hh_dtt2 &
                  &      +two*third*rs*zeta*d2hh_drsdzeta+two*seven*ninth*rs*tt*d2hh_drsdtt &
                  &      +two*seven*third*tt*zeta*d2hh_dttdzeta)
             d2rhohh_drhodg=dtt_dg*(-four*third*dhh_dtt-third*rs*d2hh_drsdtt &
                  &                           -zeta*d2hh_dttdzeta-seven*third*tt*d2hh_dtt2)

             !    Component 12 : first derivative with respect to the gradient
             !    of the density, div by the grad of the density
             dvxci(ipts,12)=dvxcdgr(ipts,3)
             !    Components 9, 10 and 11 : second derivatives with respect to the spin-density
             !    Note that there is already a contribution from LSDA
             dvxci(ipts,9)=dvxci(ipts,9)+d2rhohh_drho2+rhotot_inv*           &
                  &     ( d2hh_dzeta2*(one-two*zeta) &
                  &      -two*third*rs*d2hh_drsdzeta-14.0_dp*third*tt*d2hh_dttdzeta)
             dvxci(ipts,10)=dvxci(ipts,10)+d2rhohh_drho2-rhotot_inv*d2hh_dzeta2
             dvxci(ipts,11)=dvxci(ipts,11)+d2rhohh_drho2+rhotot_inv*           &
                  &     ( d2hh_dzeta2*(one+two*zeta) &
                  &      +two*third*rs*d2hh_drsdzeta+14.0_dp*third*tt*d2hh_dttdzeta)
             !    Components 13 and 14 : second derivatives with respect to spin density
             !    and gradient, divided by the gradient
             dvxci(ipts,13)=d2rhohh_drhodg+dtt_dg*d2hh_dttdzeta
             dvxci(ipts,14)=d2rhohh_drhodg-dtt_dg*d2hh_dttdzeta
             !    Component 15 : derivative of the (derivative wrt the gradient div by the grad),
             !    divided by the grad
             dvxci(ipts,15)=rhotot*d2hh_dtt2*dtt_dg*dtt_dg


             !  End condition of GGA


             ! Correlation has been added
             ! -----------------------------------------------------------------------------

             !  vxci(ipts,2)=vxci(ipts,1)
             dvxcdgr(ipts,2)=dvxcdgr(ipts,1)

          end do

       else if ((option==6) .or. (option==7)) then
          do ipts=1,npts

             rhotot=rhoarr(ipts)
             rhotmot=rhom1_3(ipts)
             rhotot_inv=rhotmot*rhotmot*rhotmot
             rhotmo6=sqrt(rhotmot)
             rhoto6=rhotot*rhotmot*rhotmot*rhotmo6
             ! -----------------------------------------------------------------------
             ! First take care of the exchange part of the functional

             exc=zero
             dvxcdgr(ipts,3)=zero
             !loop over the spin
             ispden=1
             rho   =rho_updn(ipts,ispden)
             rhomot=rho_updnm1_3(ipts,ispden)
             ex_lsd= - threefourth_divpi * sixpi2_1_3*rhomot*rhomot*rho
             !   Perdew-Burke-Ernzerhof GGA, exchange part
             rho_inv=rhomot*rhomot*rhomot
             coeffss=quarter*sixpi2m1_3*sixpi2m1_3*rho_inv*rho_inv*rhomot*rhomot
             ss=grho2_updn(ipts,ispden)*coeffss

             if (option==6) then
                divss=exp(-mu_divkappa*ss)
                dfxdss= mu*divss
                d2fxdss2=-mu*mu_divkappa*divss

                fx    = one+kappa*(one-divss)
                ex_gga= ex_lsd*fx
                dssdn=-eight*third*ss*rho_inv
                dfxdn  = dfxdss*dssdn
                vxci(ipts,ispden)=ex_lsd*(four_thirds*fx+rho*dfxdn)
                !   The new definition (v3.3) includes the division by the norm of the gradient
                dssdg =two*coeffss
                dfxdg=dfxdss*dssdg
                dvxcdgr(ipts,ispden)=ex_lsd*rho*dfxdg
                exc=exc+ex_gga*rho
             ! This is the Wu and Cohen modification
             else
                expss=exp(-ss)
                p1_wc=b_wc+(mu-b_wc)*(one-ss)*expss+two*c_wc*ss/(one+c_wc*ss*ss)
                p2_wc=d_wc*(ss-two)*expss+two*c_wc/(one+c_wc*ss*ss)-&
&                     four*c_wc*c_wc*ss*ss/((one+c_wc*ss*ss)*(one+c_wc*ss*ss))
                divss=one/(one+(b_wc*ss+d_wc*ss*expss+log(one+c_wc*ss*ss))/kappa)
                dfxdss=p1_wc*divss*divss
                d2fxdss2=p2_wc*divss*divss-two*divss*divss*divss*p1_wc*p1_wc/kappa

                fx    = one+kappa*(one-divss)
                ex_gga= ex_lsd*fx
                dssdn=-eight*third*ss*rho_inv
                dfxdn  = dfxdss*dssdn
                vxci(ipts,ispden)=ex_lsd*(four_thirds*fx+rho*dfxdn)
                !   The new definition (v3.3) includes the division by the norm of the gradient
                dssdg =two*coeffss
                dfxdg=dfxdss*dssdg
                dvxcdgr(ipts,ispden)=ex_lsd*rho*dfxdg
                exc=exc+ex_gga*rho
             end if

             !    Perdew-Burke-Ernzerhof GGA, exchange part
             !    Components 3 or 4
             dvxci(ipts,2+ispden)=dvxcdgr(ipts,ispden)
             !    Components 1 or 2
             d2ssdn2=-11.0_dp*third*dssdn*rho_inv
             d2fxdn2=d2fxdss2*dssdn**2+dfxdss*d2ssdn2
             dvxci(ipts,ispden)=third*rho_inv*vxci(ipts,ispden)+&
                  &     ex_lsd*(seven*third*dfxdn+rho*d2fxdn2)
             !    Components 5 or 6
             d2ssdndg=-eight*third*dssdg*rho_inv
             d2fxdndg=d2fxdss2*dssdn*dssdg+dfxdss*d2ssdndg
             dvxci(ipts,4+ispden)=ex_lsd*(four_thirds*dfxdg+rho*d2fxdndg)
             !    Components 7 or 8
             d2fxdg2=d2fxdss2*dssdg**2
             dvxci(ipts,6+ispden)=ex_lsd*rho*d2fxdg2
             !    For the time being, treat non-spin-polarized like spin-polarized
             dvxci(ipts,2)=dvxci(ipts,1)
             dvxci(ipts,4)=dvxci(ipts,3)
             dvxci(ipts,6)=dvxci(ipts,5)
             dvxci(ipts,8)=dvxci(ipts,7)
             !end of loop over the spin
             !  If non spin-polarized, treat spin down contribution now, similar to spin up
             exc=exc*2
             exci(ipts)=exc*rhotot_inv
             ! -----------------------------------------------------------------------------
             ! Then takes care of the LSD correlation part of the functional


             rs=rsfac*rhotmot
             sqr_rs=sq_rsfac*rhotmo6
             rsm1_2=sq_rsfac_inv*rhoto6

             !  Formulas A6-A8 of PW92LSD
             ec0_q0=-2.0_dp*ec0_aa*(1.0_dp+ec0_a1*rs)
             ec0_q1=2.0_dp*ec0_aa*(ec0_b1*sqr_rs+ec0_b2*rs+ec0_b3*rs*sqr_rs+ec0_b4*rs*rs)
             ec0_q1p=ec0_aa*(ec0_b1*rsm1_2+2._dp*ec0_b2+3._dp*ec0_b3*sqr_rs+4._dp*ec0_b4*rs)
             ec0_den=1.0_dp/(ec0_q1*ec0_q1+ec0_q1)
             !  ec0_log=log( 1.0_dp + 1.0_dp / ec0_q1 )
             ec0_log=-log( ec0_q1*ec0_q1*ec0_den )
             ecrs0=ec0_q0*ec0_log
             decrs0_drs= -2.0_dp*ec0_aa*ec0_a1*ec0_log - ec0_q0*ec0_q1p *ec0_den
             ec0_q1pp=0.5_dp*ec0_aa*(-ec0_b1*rsm1_2**3+3._dp*ec0_b3*rsm1_2+8._dp*ec0_b4)
             d2ecrs0_drs2= 4.0_dp*ec0_aa*ec0_a1*ec0_q1p*ec0_den            &
                  &                -ec0_q0*ec0_q1pp*ec0_den                        &
                  &                +ec0_q0*ec0_q1p**2*ec0_den**2*(2._dp*ec0_q1+1.0_dp)

             mac_q0=-2.0_dp*mac_aa*(1.0_dp+mac_a1*rs)
             mac_q1=2.0_dp*mac_aa*(mac_b1*sqr_rs+mac_b2*rs+mac_b3*rs*sqr_rs+mac_b4*rs*rs)
             mac_q1p=mac_aa*(mac_b1*rsm1_2+2._dp*mac_b2+3._dp*mac_b3*sqr_rs+4._dp*mac_b4*rs)
             mac_den=1.0_dp/(mac_q1*mac_q1+mac_q1)
             mac_log=-log( mac_q1*mac_q1*mac_den )
             macrs=mac_q0*mac_log
             dmacrs_drs= -2.0_dp*mac_aa*mac_a1*mac_log - mac_q0*mac_q1p*mac_den

             ecrs=ecrs0
             decrs_drs=decrs0_drs
             decrs_dzeta=0.0_dp
             d2ecrs_drs2=d2ecrs0_drs2
             d2ecrs_dzeta2=alpha_zeta**2*(-macrs)
             d2ecrs_drsdzeta=zero
             zeta=0.0_dp


             !  Add LSD correlation functional to GGA exchange functional
             exci(ipts)=exci(ipts)+ecrs
             vxci(ipts,1)=vxci(ipts,1)+ecrs-rs*third*decrs_drs

             dvcrs_drs=third*(2._dp*decrs_drs-rs*d2ecrs_drs2)
             !   And d(vxc)/d(rho)=(-rs/(3*rho))*d(vxc)/d(rs)
             d2ecrs_drho2= -rs**4*(four_pi*third)*third*dvcrs_drs
             dvxci(ipts,9)=d2ecrs_drho2
             dvxci(ipts,10)=d2ecrs_drho2
             dvxci(ipts,11)=d2ecrs_drho2

             !  -----------------------------------------------------------------------------
             !  Eventually add the GGA correlation part of the PBE functional
             !  Note : the computation of the potential in the spin-unpolarized
             !  case could be optimized much further. Other optimizations are left to do.

             phi_zeta=1.0_dp
             phip_zeta=0.0_dp
             phi_zeta_inv=1.0_dp
             phi_logder=0.0_dp
             phi3_zeta=1.0_dp
             gamphi3inv=gamma_inv
             phipp_zeta=-two*ninth*alpha_zeta*alpha_zeta

             !   From ec to bb
             bb=ecrs*gamphi3inv
             dbb_drs=decrs_drs*gamphi3inv
             dbb_dzeta=gamphi3inv*(decrs_dzeta-three*ecrs*phi_logder)
             d2bb_drs2=d2ecrs_drs2*gamphi3inv
             d2bb_drsdzeta=gamphi3inv*(d2ecrs_drsdzeta-three*decrs_drs*phi_logder)
             d2bb_dzeta2=gamphi3inv*(d2ecrs_dzeta2-six*decrs_dzeta*phi_logder+&
                  &     12.0_dp*ecrs*phi_logder*phi_logder-three*ecrs*phi_zeta_inv*phipp_zeta)

             !   From bb to cc
             exp_pbe=exp(-bb)
             cc=one/(exp_pbe-one)
             dcc_dbb=cc*cc*exp_pbe
             dcc_drs=dcc_dbb*dbb_drs
             dcc_dzeta=dcc_dbb*dbb_dzeta
             d2cc_dbb2=cc*cc*exp_pbe*(two*cc*exp_pbe-one)
             d2cc_drs2=d2cc_dbb2*dbb_drs*dbb_drs+dcc_dbb*d2bb_drs2
             d2cc_drsdzeta=d2cc_dbb2*dbb_drs*dbb_dzeta+dcc_dbb*d2bb_drsdzeta
             d2cc_dzeta2=d2cc_dbb2*dbb_dzeta*dbb_dzeta+dcc_dbb*d2bb_dzeta2

             !   From cc to aa
             coeff_aa=beta*gamma_inv*phi_zeta_inv*phi_zeta_inv
             aa=coeff_aa*cc
             daa_drs=coeff_aa*dcc_drs
             daa_dzeta=-two*aa*phi_logder+coeff_aa*dcc_dzeta
             d2aa_drs2=coeff_aa*d2cc_drs2
             d2aa_drsdzeta=-two*daa_drs*phi_logder+coeff_aa*d2cc_drsdzeta
             d2aa_dzeta2=aa*(-two*phi_zeta_inv*phipp_zeta+six*phi_logder*phi_logder)+&
                  &     coeff_aa*(-four*dcc_dzeta*phi_logder+d2cc_dzeta2)

             !   Introduce tt : do not assume that the spin-dependent gradients are collinear
             grrho2=four*grho2_updn(ipts,1)
             dtt_dg=two*rhotot_inv*rhotot_inv*rhotmot*coeff_tt
             !   Note that tt is (the t variable of PBE divided by phi) squared
             tt=half*grrho2*dtt_dg

             !   Get xx from aa and tt
             xx=aa*tt
             dxx_drs=daa_drs*tt
             dxx_dzeta=daa_dzeta*tt
             dxx_dtt=aa
             d2xx_drs2=d2aa_drs2*tt
             d2xx_drsdzeta=d2aa_drsdzeta*tt
             d2xx_drsdtt=daa_drs
             d2xx_dttdzeta=daa_dzeta
             d2xx_dzeta2=d2aa_dzeta2*tt

             !   From xx to pade
             pade_den=one/(one+xx*(one+xx))
             pade=(one+xx)*pade_den
             dpade_dxx=-xx*(two+xx)*pade_den**2
             dpade_drs=dpade_dxx*dxx_drs
             dpade_dtt=dpade_dxx*dxx_dtt
             dpade_dzeta=dpade_dxx*dxx_dzeta
             d2pade_dxx2=two*(-one+xx*xx*(three+xx))*pade_den*pade_den*pade_den
             d2pade_drs2=d2pade_dxx2*dxx_drs*dxx_drs+dpade_dxx*d2xx_drs2
             d2pade_drsdtt=d2pade_dxx2*dxx_drs*dxx_dtt+dpade_dxx*d2xx_drsdtt
             d2pade_drsdzeta=d2pade_dxx2*dxx_drs*dxx_dzeta+dpade_dxx*d2xx_drsdzeta
             d2pade_dtt2=d2pade_dxx2*dxx_dtt*dxx_dtt
             d2pade_dttdzeta=d2pade_dxx2*dxx_dtt*dxx_dzeta+dpade_dxx*d2xx_dttdzeta
             d2pade_dzeta2=d2pade_dxx2*dxx_dzeta*dxx_dzeta+dpade_dxx*d2xx_dzeta2

             !   From pade to qq
             coeff_qq=tt*phi_zeta_inv*phi_zeta_inv
             qq=coeff_qq*pade
             dqq_drs=coeff_qq*dpade_drs
             dqq_dtt=pade*phi_zeta_inv*phi_zeta_inv+coeff_qq*dpade_dtt
             dqq_dzeta=coeff_qq*(dpade_dzeta-two*pade*phi_logder)
             d2qq_drs2=coeff_qq*d2pade_drs2
             d2qq_drsdtt=phi_zeta_inv*phi_zeta_inv*(dpade_drs+tt*d2pade_drsdtt)
             d2qq_drsdzeta=coeff_qq*(d2pade_drsdzeta-two*dpade_drs*phi_logder)
             d2qq_dtt2=phi_zeta_inv*phi_zeta_inv*(two*dpade_dtt+tt*d2pade_dtt2)
             d2qq_dttdzeta=phi_zeta_inv*phi_zeta_inv*(dpade_dzeta-two*pade*phi_logder)+&
                  &     coeff_qq*(d2pade_dttdzeta-two*dpade_dtt*phi_logder)
             d2qq_dzeta2=coeff_qq*( d2pade_dzeta2-four*dpade_dzeta*phi_logder &
                  &                          +six*pade*phi_logder*phi_logder            &
                  &                          -two*pade*phi_zeta_inv*phipp_zeta)

             !   From qq to rr
             arg_rr=one+beta*gamma_inv*qq
             div_rr=one/arg_rr
             rr=gamma*log(arg_rr)
             drr_dqq=beta*div_rr
             drr_drs=drr_dqq*dqq_drs
             drr_dtt=drr_dqq*dqq_dtt
             drr_dzeta=drr_dqq*dqq_dzeta
             d2rr_dqq2=-div_rr**2*beta*beta*gamma_inv
             d2rr_drs2=d2rr_dqq2*dqq_drs*dqq_drs+drr_dqq*d2qq_drs2
             d2rr_drsdtt=d2rr_dqq2*dqq_drs*dqq_dtt+drr_dqq*d2qq_drsdtt
             d2rr_drsdzeta=d2rr_dqq2*dqq_drs*dqq_dzeta+drr_dqq*d2qq_drsdzeta
             d2rr_dtt2=d2rr_dqq2*dqq_dtt*dqq_dtt+drr_dqq*d2qq_dtt2
             d2rr_dttdzeta=d2rr_dqq2*dqq_dtt*dqq_dzeta+drr_dqq*d2qq_dttdzeta
             d2rr_dzeta2=d2rr_dqq2*dqq_dzeta*dqq_dzeta+drr_dqq*d2qq_dzeta2

             !   From rr to hh
             hh=phi3_zeta*rr
             dhh_drs=phi3_zeta*drr_drs
             dhh_dtt=phi3_zeta*drr_dtt
             dhh_dzeta=phi3_zeta*(drr_dzeta+three*rr*phi_logder)
             d2hh_drs2=phi3_zeta*d2rr_drs2
             d2hh_drsdtt=phi3_zeta*d2rr_drsdtt
             d2hh_drsdzeta=phi3_zeta*(d2rr_drsdzeta+three*drr_drs*phi_logder)
             d2hh_dtt2=phi3_zeta*d2rr_dtt2
             d2hh_dttdzeta=phi3_zeta*(d2rr_dttdzeta+three*drr_dtt*phi_logder)
             d2hh_dzeta2=phi3_zeta*(six*rr*phi_logder*phi_logder+&
                  &                           six*phi_logder*drr_dzeta+d2rr_dzeta2)  &
                  &               +three*phi_zeta*phi_zeta*rr*phipp_zeta

             !   The GGA correlation energy is added
             exci(ipts)=exci(ipts)+hh

             !   Change of variables : from (rs,zeta,tt) to (rhoup,rhodn,grrho)

             !   From hh to the derivative of the energy wrt the density
             drhohh_drho=hh-third*rs*dhh_drs-zeta*dhh_dzeta-seven*third*tt*dhh_dtt
             vxci(ipts,1)=vxci(ipts,1)+drhohh_drho

             !   From hh to the derivative of the energy wrt to the gradient of the
             !   density, divided by the gradient of the density
             !   (The v3.3 definition includes the division by the norm of the gradient)
             dvxcdgr(ipts,3)=rhotot*dtt_dg*dhh_dtt

             d2rhohh_drho2=rhotot_inv*&
                  &     (-two*ninth*rs*dhh_drs +seven*four*ninth*tt*dhh_dtt &
                  &      +ninth*rs*rs*d2hh_drs2+zeta*zeta*d2hh_dzeta2+(seven*third*tt)**2*d2hh_dtt2 &
                  &      +two*third*rs*zeta*d2hh_drsdzeta+two*seven*ninth*rs*tt*d2hh_drsdtt &
                  &      +two*seven*third*tt*zeta*d2hh_dttdzeta)
             d2rhohh_drhodg=dtt_dg*(-four*third*dhh_dtt-third*rs*d2hh_drsdtt &
                  &                           -zeta*d2hh_dttdzeta-seven*third*tt*d2hh_dtt2)

             !    Component 12 : first derivative with respect to the gradient
             !    of the density, div by the grad of the density
             dvxci(ipts,12)=dvxcdgr(ipts,3)
             !    Components 9, 10 and 11 : second derivatives with respect to the spin-density
             !    Note that there is already a contribution from LSDA
             dvxci(ipts,9)=dvxci(ipts,9)+d2rhohh_drho2+rhotot_inv*           &
                  &     ( d2hh_dzeta2*(one-two*zeta) &
                  &      -two*third*rs*d2hh_drsdzeta-14.0_dp*third*tt*d2hh_dttdzeta)
             dvxci(ipts,10)=dvxci(ipts,10)+d2rhohh_drho2-rhotot_inv*d2hh_dzeta2
             dvxci(ipts,11)=dvxci(ipts,11)+d2rhohh_drho2+rhotot_inv*           &
                  &     ( d2hh_dzeta2*(one+two*zeta) &
                  &      +two*third*rs*d2hh_drsdzeta+14.0_dp*third*tt*d2hh_dttdzeta)
             !    Components 13 and 14 : second derivatives with respect to spin density
             !    and gradient, divided by the gradient
             dvxci(ipts,13)=d2rhohh_drhodg+dtt_dg*d2hh_dttdzeta
             dvxci(ipts,14)=d2rhohh_drhodg-dtt_dg*d2hh_dttdzeta
             !    Component 15 : derivative of the (derivative wrt the gradient div by the grad),
             !    divided by the grad
             dvxci(ipts,15)=rhotot*d2hh_dtt2*dtt_dg*dtt_dg


             !  End condition of GGA


             ! Correlation has been added
             ! -----------------------------------------------------------------------------

             !  vxci(ipts,2)=vxci(ipts,1)
             dvxcdgr(ipts,2)=dvxcdgr(ipts,1)

          end do

       else if (option==-1) then
          do ipts=1,npts

             rhotot=rhoarr(ipts)
             rhotmot=rhom1_3(ipts)
             rhotot_inv=rhotmot*rhotmot*rhotmot
             rhotmo6=sqrt(rhotmot)
             rhoto6=rhotot*rhotmot*rhotmot*rhotmo6
             ! -----------------------------------------------------------------------
             ! First take care of the exchange part of the functional

             exc=zero
             !loop over the spin
             ispden=1
             rho   =rho_updn(ipts,ispden)
             rhomot=rho_updnm1_3(ipts,ispden)
             ex_lsd= - threefourth_divpi * sixpi2_1_3*rhomot*rhomot*rho
             !   Perdew_Wang 91 LSD
             vxci(ipts,ispden)=four_thirds*ex_lsd
             exc=exc+ex_lsd*rho

             !    Perdew_Wang 91 LSD
             dvxci(ipts,2*ispden-1)=-four_thirds*third*&
                  &                   threefourth_divpi*sixpi2_1_3*rhomot*rhomot
             dvxci(ipts,2)=zero
             !    If non-spin-polarized, first component of dvxci is second
             !    derivative with respect to TOTAL density.
             dvxci(ipts,1)=dvxci(ipts,1)*half
             ! Compute the second derivative of vx
             ! vx^(2) = -2*vx^(1)/(3*rhotot)
             !end of loop over the spin
             !  If non spin-polarized, treat spin down contribution now, similar to spin up
             exc=exc*2
             exci(ipts)=exc*rhotot_inv


             ! Correlation has been added
             ! -----------------------------------------------------------------------------


          end do

       else if (option==-2) then
          do ipts=1,npts

             rhotot=rhoarr(ipts)
             rhotmot=rhom1_3(ipts)
             rhotot_inv=rhotmot*rhotmot*rhotmot
             rhotmo6=sqrt(rhotmot)
             rhoto6=rhotot*rhotmot*rhotmot*rhotmo6
             ! -----------------------------------------------------------------------
             ! First take care of the exchange part of the functional

             exc=zero
             dvxcdgr(ipts,3)=zero
             !loop over the spin
             ispden=1
             rho   =rho_updn(ipts,ispden)
             rhomot=rho_updnm1_3(ipts,ispden)
             ex_lsd= - threefourth_divpi * sixpi2_1_3*rhomot*rhomot*rho
             !   Perdew-Burke-Ernzerhof GGA, exchange part
             rho_inv=rhomot*rhomot*rhomot
             coeffss=quarter*sixpi2m1_3*sixpi2m1_3*rho_inv*rho_inv*rhomot*rhomot
             ss=grho2_updn(ipts,ispden)*coeffss
             divss=one/(one+mu_divkappa*ss)
             dfxdss= mu*divss*divss
             d2fxdss2=-mu*two*mu_divkappa*divss*divss*divss
             fx    = one+kappa*(one-divss)
             ex_gga= ex_lsd*fx
             dssdn=-eight*third*ss*rho_inv
             dfxdn  = dfxdss*dssdn
             vxci(ipts,ispden)=ex_lsd*(four_thirds*fx+rho*dfxdn)
             !   The new definition (v3.3) includes the division by the norm of the gradient
             dssdg =two*coeffss
             dfxdg=dfxdss*dssdg
             dvxcdgr(ipts,ispden)=ex_lsd*rho*dfxdg
             exc=exc+ex_gga*rho

             !    Perdew-Burke-Ernzerhof GGA, exchange part
             !    Components 3 or 4
             dvxci(ipts,2+ispden)=dvxcdgr(ipts,ispden)
             !    Components 1 or 2
             d2ssdn2=-11.0_dp*third*dssdn*rho_inv
             d2fxdn2=d2fxdss2*dssdn**2+dfxdss*d2ssdn2
             dvxci(ipts,ispden)=third*rho_inv*vxci(ipts,ispden)+&
                  &     ex_lsd*(seven*third*dfxdn+rho*d2fxdn2)
             !    Components 5 or 6
             d2ssdndg=-eight*third*dssdg*rho_inv
             d2fxdndg=d2fxdss2*dssdn*dssdg+dfxdss*d2ssdndg
             dvxci(ipts,4+ispden)=ex_lsd*(four_thirds*dfxdg+rho*d2fxdndg)
             !    Components 7 or 8
             d2fxdg2=d2fxdss2*dssdg**2
             dvxci(ipts,6+ispden)=ex_lsd*rho*d2fxdg2
             !    For the time being, treat non-spin-polarized like spin-polarized
             dvxci(ipts,2)=dvxci(ipts,1)
             dvxci(ipts,4)=dvxci(ipts,3)
             dvxci(ipts,6)=dvxci(ipts,5)
             dvxci(ipts,8)=dvxci(ipts,7)
             !end of loop over the spin
             !  If non spin-polarized, treat spin down contribution now, similar to spin up
             exc=exc*2
             exci(ipts)=exc*rhotot_inv
             ! -----------------------------------------------------------------------------
             ! Then takes care of the LSD correlation part of the functional


             ! Correlation has been added
             ! -----------------------------------------------------------------------------

             !  vxci(ipts,2)=vxci(ipts,1)
             dvxcdgr(ipts,2)=dvxcdgr(ipts,1)

          end do

       else if (option==1) then
          do ipts=1,npts

             rhotot=rhoarr(ipts)
             rhotmot=rhom1_3(ipts)
             rhotot_inv=rhotmot*rhotmot*rhotmot
             rhotmo6=sqrt(rhotmot)
             rhoto6=rhotot*rhotmot*rhotmot*rhotmo6
             ! -----------------------------------------------------------------------
             ! First take care of the exchange part of the functional

             exc=zero
             !loop over the spin
             ispden=1
             rho   =rho_updn(ipts,ispden)
             rhomot=rho_updnm1_3(ipts,ispden)
             ex_lsd= - threefourth_divpi * sixpi2_1_3*rhomot*rhomot*rho
             !   Perdew_Wang 91 LSD
             vxci(ipts,ispden)=four_thirds*ex_lsd
             exc=exc+ex_lsd*rho


             !    Perdew_Wang 91 LSD
             dvxci(ipts,2*ispden-1)=-four_thirds*third*&
                  &                   threefourth_divpi*sixpi2_1_3*rhomot*rhomot
             dvxci(ipts,2)=zero
             !    If non-spin-polarized, first component of dvxci is second
             !    derivative with respect to TOTAL density.
             dvxci(ipts,1)=dvxci(ipts,1)*half
             ! Compute the second derivative of vx
             ! vx^(2) = -2*vx^(1)/(3*rhotot)
             !end of loop over the spin
             !  If non spin-polarized, treat spin down contribution now, similar to spin up
             exc=exc*2
             exci(ipts)=exc*rhotot_inv
             ! -----------------------------------------------------------------------------
             ! Then takes care of the LSD correlation part of the functional

             rs=rsfac*rhotmot
             sqr_rs=sq_rsfac*rhotmo6
             rsm1_2=sq_rsfac_inv*rhoto6

             !  Formulas A6-A8 of PW92LSD
             ec0_q0=-2.0_dp*ec0_aa*(1.0_dp+ec0_a1*rs)

             ec0_q1=2.0_dp*ec0_aa*(ec0_b1*sqr_rs+ec0_b2*rs+ec0_b3*rs*sqr_rs+ec0_b4*rs*rs)
             ec0_q1p=ec0_aa*(ec0_b1*rsm1_2+2._dp*ec0_b2+3._dp*ec0_b3*sqr_rs+4._dp*ec0_b4*rs)
             ec0_den=1.0_dp/(ec0_q1*ec0_q1+ec0_q1)
             !  ec0_log=log( 1.0_dp + 1.0_dp / ec0_q1 )
             ec0_log=-log( ec0_q1*ec0_q1*ec0_den )
             ecrs0=ec0_q0*ec0_log
             decrs0_drs= -2.0_dp*ec0_aa*ec0_a1*ec0_log - ec0_q0*ec0_q1p *ec0_den
             ec0_q1pp=0.5_dp*ec0_aa*(-ec0_b1*rsm1_2**3+3._dp*ec0_b3*rsm1_2+8._dp*ec0_b4)
             d2ecrs0_drs2= 4.0_dp*ec0_aa*ec0_a1*ec0_q1p*ec0_den            &
                  &                -ec0_q0*ec0_q1pp*ec0_den                        &
                  &                +ec0_q0*ec0_q1p**2*ec0_den**2*(2._dp*ec0_q1+1.0_dp)

             mac_q0=-2.0_dp*mac_aa*(1.0_dp+mac_a1*rs)
             mac_q1=2.0_dp*mac_aa*(mac_b1*sqr_rs+mac_b2*rs+mac_b3*rs*sqr_rs+mac_b4*rs*rs)
             mac_q1p=mac_aa*(mac_b1*rsm1_2+2._dp*mac_b2+3._dp*mac_b3*sqr_rs+4._dp*mac_b4*rs)
             mac_den=1.0_dp/(mac_q1*mac_q1+mac_q1)
             mac_log=-log( mac_q1*mac_q1*mac_den )
             macrs=mac_q0*mac_log
             dmacrs_drs= -2.0_dp*mac_aa*mac_a1*mac_log - mac_q0*mac_q1p*mac_den

             ecrs=ecrs0
             decrs_drs=decrs0_drs
             decrs_dzeta=0.0_dp
             d2ecrs_drs2=d2ecrs0_drs2
             d2ecrs_dzeta2=alpha_zeta**2*(-macrs)
             d2ecrs_drsdzeta=zero
             zeta=0.0_dp


             !  Add LSD correlation functional to GGA exchange functional
             exci(ipts)=exci(ipts)+ecrs
             vxci(ipts,1)=vxci(ipts,1)+ecrs-rs*third*decrs_drs

             dvcrs_drs=third*(2._dp*decrs_drs-rs*d2ecrs_drs2)
             !   And d(vxc)/d(rho)=(-rs/(3*rho))*d(vxc)/d(rs)
             d2ecrs_drho2= -rs**4*(four_pi*third)*third*dvcrs_drs

             dvxci(ipts,1)=dvxci(ipts,1)+d2ecrs_drho2
             dvxci(ipts,2)=dvxci(ipts,2)+d2ecrs_drho2-d2ecrs_dzeta2*rhotot_inv


             ! Correlation has been added
             ! -----------------------------------------------------------------------------

          end do

       else if (option==3) then
          do ipts=1,npts

             rhotot=rhoarr(ipts)
             rhotmot=rhom1_3(ipts)
             rhotot_inv=rhotmot*rhotmot*rhotmot
             rhotmo6=sqrt(rhotmot)
             rhoto6=rhotot*rhotmot*rhotmot*rhotmo6
             ! -----------------------------------------------------------------------
             ! First take care of the exchange part of the functional

             exc=zero
             !loop over the spin
             ispden=1
             rho   =rho_updn(ipts,ispden)
             rhomot=rho_updnm1_3(ipts,ispden)
             ex_lsd= - threefourth_divpi * sixpi2_1_3*rhomot*rhomot*rho
             !   Perdew_Wang 91 LSD
             vxci(ipts,ispden)=four_thirds*ex_lsd
             exc=exc+ex_lsd*rho

             !    Perdew_Wang 91 LSD
             dvxci(ipts,2*ispden-1)=-four_thirds*third*&
                  &                   threefourth_divpi*sixpi2_1_3*rhomot*rhomot
             dvxci(ipts,2)=zero
             !    If non-spin-polarized, first component of dvxci is second
             !    derivative with respect to TOTAL density.
             dvxci(ipts,1)=dvxci(ipts,1)*half
             ! Compute the second derivative of vx
             ! vx^(2) = -2*vx^(1)/(3*rhotot)
             !end of loop over the spin
             !  If non spin-polarized, treat spin down contribution now, similar to spin up
             exc=exc*2
             exci(ipts)=exc*rhotot_inv
             ! -----------------------------------------------------------------------------
             ! Then takes care of the LSD correlation part of the functional


             rs=rsfac*rhotmot
             sqr_rs=sq_rsfac*rhotmo6
             rsm1_2=sq_rsfac_inv*rhoto6

             !  Formulas A6-A8 of PW92LSD
             ec0_q0=-2.0_dp*ec0_aa*(1.0_dp+ec0_a1*rs)
             sqr_sqr_rs=max(1.e-15_dp,sqrt(sqr_rs))
             ec0_q1=2.0_dp*ec0_aa*(ec0_b1*sqr_rs+ec0_b2*rs+ec0_b3*rs*sqr_rs+ec0_b4*rs*rs/sqr_sqr_rs)
             ec0_q1p=ec0_aa*(ec0_b1*rsm1_2+2._dp*ec0_b2+3._dp*ec0_b3*sqr_rs+3.5_dp*ec0_b4*rs/sqr_sqr_rs)
             ec0_den=1.0_dp/(ec0_q1*ec0_q1+ec0_q1)
             !  ec0_log=log( 1.0_dp + 1.0_dp / ec0_q1 )
             ec0_log=-log( ec0_q1*ec0_q1*ec0_den )
             ecrs0=ec0_q0*ec0_log
             decrs0_drs= -2.0_dp*ec0_aa*ec0_a1*ec0_log - ec0_q0*ec0_q1p *ec0_den
             ec0_q1pp=0.5_dp*ec0_aa*(-ec0_b1*rsm1_2**3+3._dp*ec0_b3*rsm1_2+8._dp*ec0_b4)
             d2ecrs0_drs2= 4.0_dp*ec0_aa*ec0_a1*ec0_q1p*ec0_den            &
                  &                -ec0_q0*ec0_q1pp*ec0_den                        &
                  &                +ec0_q0*ec0_q1p**2*ec0_den**2*(2._dp*ec0_q1+1.0_dp)

             mac_q0=-2.0_dp*mac_aa*(1.0_dp+mac_a1*rs)
             mac_q1=2.0_dp*mac_aa*(mac_b1*sqr_rs+mac_b2*rs+mac_b3*rs*sqr_rs+mac_b4*rs*rs)
             mac_q1p=mac_aa*(mac_b1*rsm1_2+2._dp*mac_b2+3._dp*mac_b3*sqr_rs+4._dp*mac_b4*rs)
             mac_den=1.0_dp/(mac_q1*mac_q1+mac_q1)
             mac_log=-log( mac_q1*mac_q1*mac_den )
             macrs=mac_q0*mac_log
             dmacrs_drs= -2.0_dp*mac_aa*mac_a1*mac_log - mac_q0*mac_q1p*mac_den

             ecrs=ecrs0
             decrs_drs=decrs0_drs
             decrs_dzeta=0.0_dp
             d2ecrs_drs2=d2ecrs0_drs2
             d2ecrs_dzeta2=alpha_zeta**2*(-macrs)
             d2ecrs_drsdzeta=zero
             zeta=0.0_dp


             !  Add LSD correlation functional to GGA exchange functional
             exci(ipts)=exci(ipts)+ecrs
             vxci(ipts,1)=vxci(ipts,1)+ecrs-rs*third*decrs_drs

             dvcrs_drs=third*(2._dp*decrs_drs-rs*d2ecrs_drs2)
             !   And d(vxc)/d(rho)=(-rs/(3*rho))*d(vxc)/d(rs)
             d2ecrs_drho2= -rs**4*(four_pi*third)*third*dvcrs_drs
             dvxci(ipts,1)=dvxci(ipts,1)+d2ecrs_drho2
             dvxci(ipts,2)=dvxci(ipts,2)+d2ecrs_drho2-d2ecrs_dzeta2*rhotot_inv

             ! Correlation has been added
             ! -----------------------------------------------------------------------------

          end do

       end if


    else if (order**2>1) then
       !separate cases depending to option
       if(option==2 .or. option==5) then
          do ipts=1,npts

             rhotot=rhoarr(ipts)
             rhotmot=rhom1_3(ipts)
             rhotot_inv=rhotmot*rhotmot*rhotmot
             rhotmo6=sqrt(rhotmot)
             rhoto6=rhotot*rhotmot*rhotmot*rhotmo6
             ! -----------------------------------------------------------------------
             ! First take care of the exchange part of the functional

             exc=zero
             dvxcdgr(ipts,3)=zero
             !loop over the spin
             ispden=1
             rho   =rho_updn(ipts,ispden)
             rhomot=rho_updnm1_3(ipts,ispden)
             ex_lsd= - threefourth_divpi * sixpi2_1_3*rhomot*rhomot*rho
             !   Perdew-Burke-Ernzerhof GGA, exchange part
             rho_inv=rhomot*rhomot*rhomot
             coeffss=quarter*sixpi2m1_3*sixpi2m1_3*rho_inv*rho_inv*rhomot*rhomot
             ss=grho2_updn(ipts,ispden)*coeffss

             divss=one/(one+mu_divkappa*ss)
             dfxdss= mu*divss*divss
             d2fxdss2=-mu*two*mu_divkappa*divss*divss*divss

             fx    = one+kappa*(one-divss)
             ex_gga= ex_lsd*fx
             dssdn=-eight*third*ss*rho_inv
             dfxdn  = dfxdss*dssdn
             vxci(ipts,ispden)=ex_lsd*(four_thirds*fx+rho*dfxdn)
             !   The new definition (v3.3) includes the division by the norm of the gradient
             dssdg =two*coeffss
             dfxdg=dfxdss*dssdg
             dvxcdgr(ipts,ispden)=ex_lsd*rho*dfxdg
             exc=exc+ex_gga*rho

             !    Perdew-Burke-Ernzerhof GGA, exchange part
             !    Components 3 or 4
             dvxci(ipts,2+ispden)=dvxcdgr(ipts,ispden)
             !    Components 1 or 2
             d2ssdn2=-11.0_dp*third*dssdn*rho_inv
             d2fxdn2=d2fxdss2*dssdn**2+dfxdss*d2ssdn2
             dvxci(ipts,ispden)=third*rho_inv*vxci(ipts,ispden)+&
                  &     ex_lsd*(seven*third*dfxdn+rho*d2fxdn2)
             !    Components 5 or 6
             d2ssdndg=-eight*third*dssdg*rho_inv
             d2fxdndg=d2fxdss2*dssdn*dssdg+dfxdss*d2ssdndg
             dvxci(ipts,4+ispden)=ex_lsd*(four_thirds*dfxdg+rho*d2fxdndg)
             !    Components 7 or 8
             d2fxdg2=d2fxdss2*dssdg**2
             dvxci(ipts,6+ispden)=ex_lsd*rho*d2fxdg2
             !    For the time being, treat non-spin-polarized like spin-polarized
             dvxci(ipts,2)=dvxci(ipts,1)
             dvxci(ipts,4)=dvxci(ipts,3)
             dvxci(ipts,6)=dvxci(ipts,5)
             dvxci(ipts,8)=dvxci(ipts,7)
             !end of loop over the spin
             !  If non spin-polarized, treat spin down contribution now, similar to spin up
             exc=exc*2
             exci(ipts)=exc*rhotot_inv
             ! -----------------------------------------------------------------------------
             ! Then takes care of the LSD correlation part of the functional

             rs=rsfac*rhotmot
             sqr_rs=sq_rsfac*rhotmo6
             rsm1_2=sq_rsfac_inv*rhoto6

             !  Formulas A6-A8 of PW92LSD
             ec0_q0=-2.0_dp*ec0_aa*(1.0_dp+ec0_a1*rs)
             ec0_q1=2.0_dp*ec0_aa*(ec0_b1*sqr_rs+ec0_b2*rs+ec0_b3*rs*sqr_rs+ec0_b4*rs*rs)
             ec0_q1p=ec0_aa*(ec0_b1*rsm1_2+2._dp*ec0_b2+3._dp*ec0_b3*sqr_rs+4._dp*ec0_b4*rs)
             ec0_den=1.0_dp/(ec0_q1*ec0_q1+ec0_q1)
             !  ec0_log=log( 1.0_dp + 1.0_dp / ec0_q1 )
             ec0_log=-log( ec0_q1*ec0_q1*ec0_den )
             ecrs0=ec0_q0*ec0_log
             decrs0_drs= -2.0_dp*ec0_aa*ec0_a1*ec0_log - ec0_q0*ec0_q1p *ec0_den
             ec0_q1pp=0.5_dp*ec0_aa*(-ec0_b1*rsm1_2**3+3._dp*ec0_b3*rsm1_2+8._dp*ec0_b4)
             d2ecrs0_drs2= 4.0_dp*ec0_aa*ec0_a1*ec0_q1p*ec0_den            &
                  &                -ec0_q0*ec0_q1pp*ec0_den                        &
                  &                +ec0_q0*ec0_q1p**2*ec0_den**2*(2._dp*ec0_q1+1.0_dp)
             mac_q0=-2.0_dp*mac_aa*(1.0_dp+mac_a1*rs)
             mac_q1=2.0_dp*mac_aa*(mac_b1*sqr_rs+mac_b2*rs+mac_b3*rs*sqr_rs+mac_b4*rs*rs)
             mac_q1p=mac_aa*(mac_b1*rsm1_2+2._dp*mac_b2+3._dp*mac_b3*sqr_rs+4._dp*mac_b4*rs)
             mac_den=1.0_dp/(mac_q1*mac_q1+mac_q1)
             mac_log=-log( mac_q1*mac_q1*mac_den )
             macrs=mac_q0*mac_log
             dmacrs_drs= -2.0_dp*mac_aa*mac_a1*mac_log - mac_q0*mac_q1p*mac_den
             ecrs=ecrs0
             decrs_drs=decrs0_drs
             decrs_dzeta=0.0_dp
             d2ecrs_drs2=d2ecrs0_drs2
             d2ecrs_dzeta2=alpha_zeta**2*(-macrs)
             d2ecrs_drsdzeta=zero
             zeta=0.0_dp

             !  Add LSD correlation functional to GGA exchange functional
             exci(ipts)=exci(ipts)+ecrs
             vxci(ipts,1)=vxci(ipts,1)+ecrs-rs*third*decrs_drs

             dvcrs_drs=third*(2._dp*decrs_drs-rs*d2ecrs_drs2)
             !   And d(vxc)/d(rho)=(-rs/(3*rho))*d(vxc)/d(rs)
             d2ecrs_drho2= -rs**4*(four_pi*third)*third*dvcrs_drs
             dvxci(ipts,9)=d2ecrs_drho2
             dvxci(ipts,10)=d2ecrs_drho2
             dvxci(ipts,11)=d2ecrs_drho2


             !  -----------------------------------------------------------------------------
             !  Eventually add the GGA correlation part of the PBE functional
             !  Note : the computation of the potential in the spin-unpolarized
             !  case could be optimized much further. Other optimizations are left to do.

             phi_zeta=1.0_dp
             phip_zeta=0.0_dp
             phi_zeta_inv=1.0_dp
             phi_logder=0.0_dp
             phi3_zeta=1.0_dp
             gamphi3inv=gamma_inv
             phipp_zeta=-two*ninth*alpha_zeta*alpha_zeta

             !   From ec to bb
             bb=ecrs*gamphi3inv
             dbb_drs=decrs_drs*gamphi3inv
             dbb_dzeta=gamphi3inv*(decrs_dzeta-three*ecrs*phi_logder)
             d2bb_drs2=d2ecrs_drs2*gamphi3inv
             d2bb_drsdzeta=gamphi3inv*(d2ecrs_drsdzeta-three*decrs_drs*phi_logder)
             d2bb_dzeta2=gamphi3inv*(d2ecrs_dzeta2-six*decrs_dzeta*phi_logder+&
                  &     12.0_dp*ecrs*phi_logder*phi_logder-three*ecrs*phi_zeta_inv*phipp_zeta)

             !   From bb to cc
             exp_pbe=exp(-bb)
             cc=one/(exp_pbe-one)
             dcc_dbb=cc*cc*exp_pbe
             dcc_drs=dcc_dbb*dbb_drs
             dcc_dzeta=dcc_dbb*dbb_dzeta
             d2cc_dbb2=cc*cc*exp_pbe*(two*cc*exp_pbe-one)
             d2cc_drs2=d2cc_dbb2*dbb_drs*dbb_drs+dcc_dbb*d2bb_drs2
             d2cc_drsdzeta=d2cc_dbb2*dbb_drs*dbb_dzeta+dcc_dbb*d2bb_drsdzeta
             d2cc_dzeta2=d2cc_dbb2*dbb_dzeta*dbb_dzeta+dcc_dbb*d2bb_dzeta2

             !   From cc to aa
             coeff_aa=beta*gamma_inv*phi_zeta_inv*phi_zeta_inv
             aa=coeff_aa*cc
             daa_drs=coeff_aa*dcc_drs
             daa_dzeta=-two*aa*phi_logder+coeff_aa*dcc_dzeta
             d2aa_drs2=coeff_aa*d2cc_drs2
             d2aa_drsdzeta=-two*daa_drs*phi_logder+coeff_aa*d2cc_drsdzeta
             d2aa_dzeta2=aa*(-two*phi_zeta_inv*phipp_zeta+six*phi_logder*phi_logder)+&
                  &     coeff_aa*(-four*dcc_dzeta*phi_logder+d2cc_dzeta2)

             !   Introduce tt : do not assume that the spin-dependent gradients are collinear
             grrho2=four*grho2_updn(ipts,1)
             dtt_dg=two*rhotot_inv*rhotot_inv*rhotmot*coeff_tt
             !   Note that tt is (the t variable of PBE divided by phi) squared
             tt=half*grrho2*dtt_dg

             !   Get xx from aa and tt
             xx=aa*tt
             dxx_drs=daa_drs*tt
             dxx_dzeta=daa_dzeta*tt
             dxx_dtt=aa
             d2xx_drs2=d2aa_drs2*tt
             d2xx_drsdzeta=d2aa_drsdzeta*tt
             d2xx_drsdtt=daa_drs
             d2xx_dttdzeta=daa_dzeta
             d2xx_dzeta2=d2aa_dzeta2*tt

             !   From xx to pade
             pade_den=one/(one+xx*(one+xx))
             pade=(one+xx)*pade_den
             dpade_dxx=-xx*(two+xx)*pade_den**2
             dpade_drs=dpade_dxx*dxx_drs
             dpade_dtt=dpade_dxx*dxx_dtt
             dpade_dzeta=dpade_dxx*dxx_dzeta
             d2pade_dxx2=two*(-one+xx*xx*(three+xx))*pade_den*pade_den*pade_den
             d2pade_drs2=d2pade_dxx2*dxx_drs*dxx_drs+dpade_dxx*d2xx_drs2
             d2pade_drsdtt=d2pade_dxx2*dxx_drs*dxx_dtt+dpade_dxx*d2xx_drsdtt
             d2pade_drsdzeta=d2pade_dxx2*dxx_drs*dxx_dzeta+dpade_dxx*d2xx_drsdzeta
             d2pade_dtt2=d2pade_dxx2*dxx_dtt*dxx_dtt
             d2pade_dttdzeta=d2pade_dxx2*dxx_dtt*dxx_dzeta+dpade_dxx*d2xx_dttdzeta
             d2pade_dzeta2=d2pade_dxx2*dxx_dzeta*dxx_dzeta+dpade_dxx*d2xx_dzeta2


             !   From pade to qq
             coeff_qq=tt*phi_zeta_inv*phi_zeta_inv
             qq=coeff_qq*pade
             dqq_drs=coeff_qq*dpade_drs
             dqq_dtt=pade*phi_zeta_inv*phi_zeta_inv+coeff_qq*dpade_dtt
             dqq_dzeta=coeff_qq*(dpade_dzeta-two*pade*phi_logder)
             d2qq_drs2=coeff_qq*d2pade_drs2
             d2qq_drsdtt=phi_zeta_inv*phi_zeta_inv*(dpade_drs+tt*d2pade_drsdtt)
             d2qq_drsdzeta=coeff_qq*(d2pade_drsdzeta-two*dpade_drs*phi_logder)
             d2qq_dtt2=phi_zeta_inv*phi_zeta_inv*(two*dpade_dtt+tt*d2pade_dtt2)
             d2qq_dttdzeta=phi_zeta_inv*phi_zeta_inv*(dpade_dzeta-two*pade*phi_logder)+&
                  &     coeff_qq*(d2pade_dttdzeta-two*dpade_dtt*phi_logder)
             d2qq_dzeta2=coeff_qq*( d2pade_dzeta2-four*dpade_dzeta*phi_logder &
                  &                          +six*pade*phi_logder*phi_logder            &
                  &                          -two*pade*phi_zeta_inv*phipp_zeta)

             !   From qq to rr
             arg_rr=one+beta*gamma_inv*qq
             div_rr=one/arg_rr
             rr=gamma*log(arg_rr)
             drr_dqq=beta*div_rr
             drr_drs=drr_dqq*dqq_drs
             drr_dtt=drr_dqq*dqq_dtt
             drr_dzeta=drr_dqq*dqq_dzeta
             d2rr_dqq2=-div_rr**2*beta*beta*gamma_inv
             d2rr_drs2=d2rr_dqq2*dqq_drs*dqq_drs+drr_dqq*d2qq_drs2
             d2rr_drsdtt=d2rr_dqq2*dqq_drs*dqq_dtt+drr_dqq*d2qq_drsdtt
             d2rr_drsdzeta=d2rr_dqq2*dqq_drs*dqq_dzeta+drr_dqq*d2qq_drsdzeta
             d2rr_dtt2=d2rr_dqq2*dqq_dtt*dqq_dtt+drr_dqq*d2qq_dtt2
             d2rr_dttdzeta=d2rr_dqq2*dqq_dtt*dqq_dzeta+drr_dqq*d2qq_dttdzeta
             d2rr_dzeta2=d2rr_dqq2*dqq_dzeta*dqq_dzeta+drr_dqq*d2qq_dzeta2

             !   From rr to hh
             hh=phi3_zeta*rr
             dhh_drs=phi3_zeta*drr_drs
             dhh_dtt=phi3_zeta*drr_dtt
             dhh_dzeta=phi3_zeta*(drr_dzeta+three*rr*phi_logder)
             d2hh_drs2=phi3_zeta*d2rr_drs2
             d2hh_drsdtt=phi3_zeta*d2rr_drsdtt
             d2hh_drsdzeta=phi3_zeta*(d2rr_drsdzeta+three*drr_drs*phi_logder)
             d2hh_dtt2=phi3_zeta*d2rr_dtt2
             d2hh_dttdzeta=phi3_zeta*(d2rr_dttdzeta+three*drr_dtt*phi_logder)
             d2hh_dzeta2=phi3_zeta*(six*rr*phi_logder*phi_logder+&
                  &                           six*phi_logder*drr_dzeta+d2rr_dzeta2)  &
                  &               +three*phi_zeta*phi_zeta*rr*phipp_zeta


             !   The GGA correlation energy is added
             exci(ipts)=exci(ipts)+hh

             !   Change of variables : from (rs,zeta,tt) to (rhoup,rhodn,grrho)


             !   From hh to the derivative of the energy wrt the density
             drhohh_drho=hh-third*rs*dhh_drs-zeta*dhh_dzeta-seven*third*tt*dhh_dtt
             vxci(ipts,1)=vxci(ipts,1)+drhohh_drho

             !   From hh to the derivative of the energy wrt to the gradient of the
             !   density, divided by the gradient of the density
             !   (The v3.3 definition includes the division by the norm of the gradient)
             dvxcdgr(ipts,3)=rhotot*dtt_dg*dhh_dtt

             d2rhohh_drho2=rhotot_inv*&
                  &     (-two*ninth*rs*dhh_drs +seven*four*ninth*tt*dhh_dtt &
                  &      +ninth*rs*rs*d2hh_drs2+zeta*zeta*d2hh_dzeta2+(seven*third*tt)**2*d2hh_dtt2 &
                  &      +two*third*rs*zeta*d2hh_drsdzeta+two*seven*ninth*rs*tt*d2hh_drsdtt &
                  &      +two*seven*third*tt*zeta*d2hh_dttdzeta)
             d2rhohh_drhodg=dtt_dg*(-four*third*dhh_dtt-third*rs*d2hh_drsdtt &
                  &                           -zeta*d2hh_dttdzeta-seven*third*tt*d2hh_dtt2)

             !    Component 12 : first derivative with respect to the gradient
             !    of the density, div by the grad of the density
             dvxci(ipts,12)=dvxcdgr(ipts,3)
             !    Components 9, 10 and 11 : second derivatives with respect to the spin-density
             !    Note that there is already a contribution from LSDA
             dvxci(ipts,9)=dvxci(ipts,9)+d2rhohh_drho2+rhotot_inv*           &
                  &     ( d2hh_dzeta2*(one-two*zeta) &
                  &      -two*third*rs*d2hh_drsdzeta-14.0_dp*third*tt*d2hh_dttdzeta)
             dvxci(ipts,10)=dvxci(ipts,10)+d2rhohh_drho2-rhotot_inv*d2hh_dzeta2
             dvxci(ipts,11)=dvxci(ipts,11)+d2rhohh_drho2+rhotot_inv*           &
                  &     ( d2hh_dzeta2*(one+two*zeta) &
                  &      +two*third*rs*d2hh_drsdzeta+14.0_dp*third*tt*d2hh_dttdzeta)
             !    Components 13 and 14 : second derivatives with respect to spin density
             !    and gradient, divided by the gradient
             dvxci(ipts,13)=d2rhohh_drhodg+dtt_dg*d2hh_dttdzeta
             dvxci(ipts,14)=d2rhohh_drhodg-dtt_dg*d2hh_dttdzeta
             !    Component 15 : derivative of the (derivative wrt the gradient div by the grad),
             !    divided by the grad
             dvxci(ipts,15)=rhotot*d2hh_dtt2*dtt_dg*dtt_dg


             !  End condition of GGA

             ! Correlation has been added
             ! -----------------------------------------------------------------------------

             !  vxci(ipts,2)=vxci(ipts,1)
             dvxcdgr(ipts,2)=dvxcdgr(ipts,1)

          end do

       else if ((option==6) .or. (option==7)) then
          do ipts=1,npts

             rhotot=rhoarr(ipts)
             rhotmot=rhom1_3(ipts)
             rhotot_inv=rhotmot*rhotmot*rhotmot
             rhotmo6=sqrt(rhotmot)
             rhoto6=rhotot*rhotmot*rhotmot*rhotmo6
             ! -----------------------------------------------------------------------
             ! First take care of the exchange part of the functional

             exc=zero
             dvxcdgr(ipts,3)=zero
             !loop over the spin
             ispden=1
             rho   =rho_updn(ipts,ispden)
             rhomot=rho_updnm1_3(ipts,ispden)
             ex_lsd= - threefourth_divpi * sixpi2_1_3*rhomot*rhomot*rho
             !   Perdew-Burke-Ernzerhof GGA, exchange part
             rho_inv=rhomot*rhomot*rhomot
             coeffss=quarter*sixpi2m1_3*sixpi2m1_3*rho_inv*rho_inv*rhomot*rhomot
             ss=grho2_updn(ipts,ispden)*coeffss

             if (option==6) then
                divss=exp(-mu_divkappa*ss)
                dfxdss= mu*divss
                d2fxdss2=-mu*mu_divkappa*divss

                fx    = one+kappa*(one-divss)
                ex_gga= ex_lsd*fx
                dssdn=-eight*third*ss*rho_inv
                dfxdn  = dfxdss*dssdn
                vxci(ipts,ispden)=ex_lsd*(four_thirds*fx+rho*dfxdn)
                !   The new definition (v3.3) includes the division by the norm of the gradient
                dssdg =two*coeffss
                dfxdg=dfxdss*dssdg
                dvxcdgr(ipts,ispden)=ex_lsd*rho*dfxdg
                exc=exc+ex_gga*rho
             ! This is the Wu and Cohen modification
             else
                expss=exp(-ss)
                p1_wc=b_wc+(mu-b_wc)*(one-ss)*expss+two*c_wc*ss/(one+c_wc*ss*ss)
                p2_wc=d_wc*(ss-two)*expss+two*c_wc/(one+c_wc*ss*ss)-&
&                     four*c_wc*c_wc*ss*ss/((one+c_wc*ss*ss)*(one+c_wc*ss*ss))
                divss=one/(one+(b_wc*ss+d_wc*ss*expss+log(one+c_wc*ss*ss))/kappa)
                dfxdss=p1_wc*divss*divss
                d2fxdss2=p2_wc*divss*divss-two*divss*divss*divss*p1_wc*p1_wc/kappa

                fx    = one+kappa*(one-divss)
                ex_gga= ex_lsd*fx
                dssdn=-eight*third*ss*rho_inv
                dfxdn  = dfxdss*dssdn
                vxci(ipts,ispden)=ex_lsd*(four_thirds*fx+rho*dfxdn)
                !   The new definition (v3.3) includes the division by the norm of the gradient
                dssdg =two*coeffss
                dfxdg=dfxdss*dssdg
                dvxcdgr(ipts,ispden)=ex_lsd*rho*dfxdg
                exc=exc+ex_gga*rho
             end if

             !    Perdew-Burke-Ernzerhof GGA, exchange part
             !    Components 3 or 4
             dvxci(ipts,2+ispden)=dvxcdgr(ipts,ispden)
             !    Components 1 or 2
             d2ssdn2=-11.0_dp*third*dssdn*rho_inv
             d2fxdn2=d2fxdss2*dssdn**2+dfxdss*d2ssdn2
             dvxci(ipts,ispden)=third*rho_inv*vxci(ipts,ispden)+&
                  &     ex_lsd*(seven*third*dfxdn+rho*d2fxdn2)
             !    Components 5 or 6
             d2ssdndg=-eight*third*dssdg*rho_inv
             d2fxdndg=d2fxdss2*dssdn*dssdg+dfxdss*d2ssdndg
             dvxci(ipts,4+ispden)=ex_lsd*(four_thirds*dfxdg+rho*d2fxdndg)
             !    Components 7 or 8
             d2fxdg2=d2fxdss2*dssdg**2
             dvxci(ipts,6+ispden)=ex_lsd*rho*d2fxdg2
             !    For the time being, treat non-spin-polarized like spin-polarized
             dvxci(ipts,2)=dvxci(ipts,1)
             dvxci(ipts,4)=dvxci(ipts,3)
             dvxci(ipts,6)=dvxci(ipts,5)
             dvxci(ipts,8)=dvxci(ipts,7)
             !end of loop over the spin
             !  If non spin-polarized, treat spin down contribution now, similar to spin up
             exc=exc*2
             exci(ipts)=exc*rhotot_inv
             ! -----------------------------------------------------------------------------
             ! Then takes care of the LSD correlation part of the functional


             rs=rsfac*rhotmot
             sqr_rs=sq_rsfac*rhotmo6
             rsm1_2=sq_rsfac_inv*rhoto6

             !  Formulas A6-A8 of PW92LSD
             ec0_q0=-2.0_dp*ec0_aa*(1.0_dp+ec0_a1*rs)
             ec0_q1=2.0_dp*ec0_aa*(ec0_b1*sqr_rs+ec0_b2*rs+ec0_b3*rs*sqr_rs+ec0_b4*rs*rs)
             ec0_q1p=ec0_aa*(ec0_b1*rsm1_2+2._dp*ec0_b2+3._dp*ec0_b3*sqr_rs+4._dp*ec0_b4*rs)
             ec0_den=1.0_dp/(ec0_q1*ec0_q1+ec0_q1)
             !  ec0_log=log( 1.0_dp + 1.0_dp / ec0_q1 )
             ec0_log=-log( ec0_q1*ec0_q1*ec0_den )
             ecrs0=ec0_q0*ec0_log
             decrs0_drs= -2.0_dp*ec0_aa*ec0_a1*ec0_log - ec0_q0*ec0_q1p *ec0_den
             ec0_q1pp=0.5_dp*ec0_aa*(-ec0_b1*rsm1_2**3+3._dp*ec0_b3*rsm1_2+8._dp*ec0_b4)
             d2ecrs0_drs2= 4.0_dp*ec0_aa*ec0_a1*ec0_q1p*ec0_den            &
                  &                -ec0_q0*ec0_q1pp*ec0_den                        &
                  &                +ec0_q0*ec0_q1p**2*ec0_den**2*(2._dp*ec0_q1+1.0_dp)
             mac_q0=-2.0_dp*mac_aa*(1.0_dp+mac_a1*rs)
             mac_q1=2.0_dp*mac_aa*(mac_b1*sqr_rs+mac_b2*rs+mac_b3*rs*sqr_rs+mac_b4*rs*rs)
             mac_q1p=mac_aa*(mac_b1*rsm1_2+2._dp*mac_b2+3._dp*mac_b3*sqr_rs+4._dp*mac_b4*rs)
             mac_den=1.0_dp/(mac_q1*mac_q1+mac_q1)
             mac_log=-log( mac_q1*mac_q1*mac_den )
             macrs=mac_q0*mac_log
             dmacrs_drs= -2.0_dp*mac_aa*mac_a1*mac_log - mac_q0*mac_q1p*mac_den
             ecrs=ecrs0
             decrs_drs=decrs0_drs
             decrs_dzeta=0.0_dp
             d2ecrs_drs2=d2ecrs0_drs2
             d2ecrs_dzeta2=alpha_zeta**2*(-macrs)
             d2ecrs_drsdzeta=zero
             zeta=0.0_dp

             !  Add LSD correlation functional to GGA exchange functional
             exci(ipts)=exci(ipts)+ecrs
             vxci(ipts,1)=vxci(ipts,1)+ecrs-rs*third*decrs_drs

             dvcrs_drs=third*(2._dp*decrs_drs-rs*d2ecrs_drs2)
             !   And d(vxc)/d(rho)=(-rs/(3*rho))*d(vxc)/d(rs)
             d2ecrs_drho2= -rs**4*(four_pi*third)*third*dvcrs_drs
             dvxci(ipts,9)=d2ecrs_drho2
             dvxci(ipts,10)=d2ecrs_drho2
             dvxci(ipts,11)=d2ecrs_drho2


             !  -----------------------------------------------------------------------------
             !  Eventually add the GGA correlation part of the PBE functional
             !  Note : the computation of the potential in the spin-unpolarized
             !  case could be optimized much further. Other optimizations are left to do.

             phi_zeta=1.0_dp
             phip_zeta=0.0_dp
             phi_zeta_inv=1.0_dp
             phi_logder=0.0_dp
             phi3_zeta=1.0_dp
             gamphi3inv=gamma_inv
             phipp_zeta=-two*ninth*alpha_zeta*alpha_zeta

             !   From ec to bb
             bb=ecrs*gamphi3inv
             dbb_drs=decrs_drs*gamphi3inv
             dbb_dzeta=gamphi3inv*(decrs_dzeta-three*ecrs*phi_logder)
             d2bb_drs2=d2ecrs_drs2*gamphi3inv
             d2bb_drsdzeta=gamphi3inv*(d2ecrs_drsdzeta-three*decrs_drs*phi_logder)
             d2bb_dzeta2=gamphi3inv*(d2ecrs_dzeta2-six*decrs_dzeta*phi_logder+&
                  &     12.0_dp*ecrs*phi_logder*phi_logder-three*ecrs*phi_zeta_inv*phipp_zeta)

             !   From bb to cc
             exp_pbe=exp(-bb)
             cc=one/(exp_pbe-one)
             dcc_dbb=cc*cc*exp_pbe
             dcc_drs=dcc_dbb*dbb_drs
             dcc_dzeta=dcc_dbb*dbb_dzeta
             d2cc_dbb2=cc*cc*exp_pbe*(two*cc*exp_pbe-one)
             d2cc_drs2=d2cc_dbb2*dbb_drs*dbb_drs+dcc_dbb*d2bb_drs2
             d2cc_drsdzeta=d2cc_dbb2*dbb_drs*dbb_dzeta+dcc_dbb*d2bb_drsdzeta
             d2cc_dzeta2=d2cc_dbb2*dbb_dzeta*dbb_dzeta+dcc_dbb*d2bb_dzeta2

             !   From cc to aa
             coeff_aa=beta*gamma_inv*phi_zeta_inv*phi_zeta_inv
             aa=coeff_aa*cc
             daa_drs=coeff_aa*dcc_drs
             daa_dzeta=-two*aa*phi_logder+coeff_aa*dcc_dzeta
             d2aa_drs2=coeff_aa*d2cc_drs2
             d2aa_drsdzeta=-two*daa_drs*phi_logder+coeff_aa*d2cc_drsdzeta
             d2aa_dzeta2=aa*(-two*phi_zeta_inv*phipp_zeta+six*phi_logder*phi_logder)+&
                  &     coeff_aa*(-four*dcc_dzeta*phi_logder+d2cc_dzeta2)

             !   Introduce tt : do not assume that the spin-dependent gradients are collinear
             grrho2=four*grho2_updn(ipts,1)
             dtt_dg=two*rhotot_inv*rhotot_inv*rhotmot*coeff_tt
             !   Note that tt is (the t variable of PBE divided by phi) squared
             tt=half*grrho2*dtt_dg

             !   Get xx from aa and tt
             xx=aa*tt
             dxx_drs=daa_drs*tt
             dxx_dzeta=daa_dzeta*tt
             dxx_dtt=aa
             d2xx_drs2=d2aa_drs2*tt
             d2xx_drsdzeta=d2aa_drsdzeta*tt
             d2xx_drsdtt=daa_drs
             d2xx_dttdzeta=daa_dzeta
             d2xx_dzeta2=d2aa_dzeta2*tt

             !   From xx to pade
             pade_den=one/(one+xx*(one+xx))
             pade=(one+xx)*pade_den
             dpade_dxx=-xx*(two+xx)*pade_den**2
             dpade_drs=dpade_dxx*dxx_drs
             dpade_dtt=dpade_dxx*dxx_dtt
             dpade_dzeta=dpade_dxx*dxx_dzeta
             d2pade_dxx2=two*(-one+xx*xx*(three+xx))*pade_den*pade_den*pade_den
             d2pade_drs2=d2pade_dxx2*dxx_drs*dxx_drs+dpade_dxx*d2xx_drs2
             d2pade_drsdtt=d2pade_dxx2*dxx_drs*dxx_dtt+dpade_dxx*d2xx_drsdtt
             d2pade_drsdzeta=d2pade_dxx2*dxx_drs*dxx_dzeta+dpade_dxx*d2xx_drsdzeta
             d2pade_dtt2=d2pade_dxx2*dxx_dtt*dxx_dtt
             d2pade_dttdzeta=d2pade_dxx2*dxx_dtt*dxx_dzeta+dpade_dxx*d2xx_dttdzeta
             d2pade_dzeta2=d2pade_dxx2*dxx_dzeta*dxx_dzeta+dpade_dxx*d2xx_dzeta2


             !   From pade to qq
             coeff_qq=tt*phi_zeta_inv*phi_zeta_inv
             qq=coeff_qq*pade
             dqq_drs=coeff_qq*dpade_drs
             dqq_dtt=pade*phi_zeta_inv*phi_zeta_inv+coeff_qq*dpade_dtt
             dqq_dzeta=coeff_qq*(dpade_dzeta-two*pade*phi_logder)
             d2qq_drs2=coeff_qq*d2pade_drs2
             d2qq_drsdtt=phi_zeta_inv*phi_zeta_inv*(dpade_drs+tt*d2pade_drsdtt)
             d2qq_drsdzeta=coeff_qq*(d2pade_drsdzeta-two*dpade_drs*phi_logder)
             d2qq_dtt2=phi_zeta_inv*phi_zeta_inv*(two*dpade_dtt+tt*d2pade_dtt2)
             d2qq_dttdzeta=phi_zeta_inv*phi_zeta_inv*(dpade_dzeta-two*pade*phi_logder)+&
                  &     coeff_qq*(d2pade_dttdzeta-two*dpade_dtt*phi_logder)
             d2qq_dzeta2=coeff_qq*( d2pade_dzeta2-four*dpade_dzeta*phi_logder &
                  &                          +six*pade*phi_logder*phi_logder            &
                  &                          -two*pade*phi_zeta_inv*phipp_zeta)

             !   From qq to rr
             arg_rr=one+beta*gamma_inv*qq
             div_rr=one/arg_rr
             rr=gamma*log(arg_rr)
             drr_dqq=beta*div_rr
             drr_drs=drr_dqq*dqq_drs
             drr_dtt=drr_dqq*dqq_dtt
             drr_dzeta=drr_dqq*dqq_dzeta
             d2rr_dqq2=-div_rr**2*beta*beta*gamma_inv
             d2rr_drs2=d2rr_dqq2*dqq_drs*dqq_drs+drr_dqq*d2qq_drs2
             d2rr_drsdtt=d2rr_dqq2*dqq_drs*dqq_dtt+drr_dqq*d2qq_drsdtt
             d2rr_drsdzeta=d2rr_dqq2*dqq_drs*dqq_dzeta+drr_dqq*d2qq_drsdzeta
             d2rr_dtt2=d2rr_dqq2*dqq_dtt*dqq_dtt+drr_dqq*d2qq_dtt2
             d2rr_dttdzeta=d2rr_dqq2*dqq_dtt*dqq_dzeta+drr_dqq*d2qq_dttdzeta
             d2rr_dzeta2=d2rr_dqq2*dqq_dzeta*dqq_dzeta+drr_dqq*d2qq_dzeta2

             !   From rr to hh
             hh=phi3_zeta*rr
             dhh_drs=phi3_zeta*drr_drs
             dhh_dtt=phi3_zeta*drr_dtt
             dhh_dzeta=phi3_zeta*(drr_dzeta+three*rr*phi_logder)
             d2hh_drs2=phi3_zeta*d2rr_drs2
             d2hh_drsdtt=phi3_zeta*d2rr_drsdtt
             d2hh_drsdzeta=phi3_zeta*(d2rr_drsdzeta+three*drr_drs*phi_logder)
             d2hh_dtt2=phi3_zeta*d2rr_dtt2
             d2hh_dttdzeta=phi3_zeta*(d2rr_dttdzeta+three*drr_dtt*phi_logder)
             d2hh_dzeta2=phi3_zeta*(six*rr*phi_logder*phi_logder+&
                  &                           six*phi_logder*drr_dzeta+d2rr_dzeta2)  &
                  &               +three*phi_zeta*phi_zeta*rr*phipp_zeta


             !   The GGA correlation energy is added
             exci(ipts)=exci(ipts)+hh

             !   Change of variables : from (rs,zeta,tt) to (rhoup,rhodn,grrho)


             !   From hh to the derivative of the energy wrt the density
             drhohh_drho=hh-third*rs*dhh_drs-zeta*dhh_dzeta-seven*third*tt*dhh_dtt
             vxci(ipts,1)=vxci(ipts,1)+drhohh_drho

             !   From hh to the derivative of the energy wrt to the gradient of the
             !   density, divided by the gradient of the density
             !   (The v3.3 definition includes the division by the norm of the gradient)
             dvxcdgr(ipts,3)=rhotot*dtt_dg*dhh_dtt

             d2rhohh_drho2=rhotot_inv*&
                  &     (-two*ninth*rs*dhh_drs +seven*four*ninth*tt*dhh_dtt &
                  &      +ninth*rs*rs*d2hh_drs2+zeta*zeta*d2hh_dzeta2+(seven*third*tt)**2*d2hh_dtt2 &
                  &      +two*third*rs*zeta*d2hh_drsdzeta+two*seven*ninth*rs*tt*d2hh_drsdtt &
                  &      +two*seven*third*tt*zeta*d2hh_dttdzeta)
             d2rhohh_drhodg=dtt_dg*(-four*third*dhh_dtt-third*rs*d2hh_drsdtt &
                  &                           -zeta*d2hh_dttdzeta-seven*third*tt*d2hh_dtt2)

             !    Component 12 : first derivative with respect to the gradient
             !    of the density, div by the grad of the density
             dvxci(ipts,12)=dvxcdgr(ipts,3)
             !    Components 9, 10 and 11 : second derivatives with respect to the spin-density
             !    Note that there is already a contribution from LSDA
             dvxci(ipts,9)=dvxci(ipts,9)+d2rhohh_drho2+rhotot_inv*           &
                  &     ( d2hh_dzeta2*(one-two*zeta) &
                  &      -two*third*rs*d2hh_drsdzeta-14.0_dp*third*tt*d2hh_dttdzeta)
             dvxci(ipts,10)=dvxci(ipts,10)+d2rhohh_drho2-rhotot_inv*d2hh_dzeta2
             dvxci(ipts,11)=dvxci(ipts,11)+d2rhohh_drho2+rhotot_inv*           &
                  &     ( d2hh_dzeta2*(one+two*zeta) &
                  &      +two*third*rs*d2hh_drsdzeta+14.0_dp*third*tt*d2hh_dttdzeta)
             !    Components 13 and 14 : second derivatives with respect to spin density
             !    and gradient, divided by the gradient
             dvxci(ipts,13)=d2rhohh_drhodg+dtt_dg*d2hh_dttdzeta
             dvxci(ipts,14)=d2rhohh_drhodg-dtt_dg*d2hh_dttdzeta
             !    Component 15 : derivative of the (derivative wrt the gradient div by the grad),
             !    divided by the grad
             dvxci(ipts,15)=rhotot*d2hh_dtt2*dtt_dg*dtt_dg


             !  End condition of GGA

             ! Correlation has been added
             ! -----------------------------------------------------------------------------

             !  vxci(ipts,2)=vxci(ipts,1)
             dvxcdgr(ipts,2)=dvxcdgr(ipts,1)

          end do

       else if (option==-1) then
          do ipts=1,npts

             rhotot=rhoarr(ipts)
             rhotmot=rhom1_3(ipts)
             rhotot_inv=rhotmot*rhotmot*rhotmot
             rhotmo6=sqrt(rhotmot)
             rhoto6=rhotot*rhotmot*rhotmot*rhotmo6
             ! -----------------------------------------------------------------------
             ! First take care of the exchange part of the functional

             exc=zero
             !loop over the spin
             ispden=1
             rho   =rho_updn(ipts,ispden)
             rhomot=rho_updnm1_3(ipts,ispden)
             ex_lsd= - threefourth_divpi * sixpi2_1_3*rhomot*rhomot*rho
             !   Perdew_Wang 91 LSD
             vxci(ipts,ispden)=four_thirds*ex_lsd
             exc=exc+ex_lsd*rho

             !    Perdew_Wang 91 LSD
             dvxci(ipts,2*ispden-1)=-four_thirds*third*&
                  &                   threefourth_divpi*sixpi2_1_3*rhomot*rhomot
             dvxci(ipts,2)=zero
             !    If non-spin-polarized, first component of dvxci is second
             !    derivative with respect to TOTAL density.
             dvxci(ipts,1)=dvxci(ipts,1)*half
             ! Compute the second derivative of vx
             ! vx^(2) = -2*vx^(1)/(3*rhotot)
             !end of loop over the spin
             !  If non spin-polarized, treat spin down contribution now, similar to spin up
             exc=exc*2
             exci(ipts)=exc*rhotot_inv
             ! -----------------------------------------------------------------------------
             ! Then takes care of the LSD correlation part of the functional

          end do

       else if (option==-2) then
          do ipts=1,npts

             rhotot=rhoarr(ipts)
             rhotmot=rhom1_3(ipts)
             rhotot_inv=rhotmot*rhotmot*rhotmot
             rhotmo6=sqrt(rhotmot)
             rhoto6=rhotot*rhotmot*rhotmot*rhotmo6
             ! -----------------------------------------------------------------------
             ! First take care of the exchange part of the functional

             exc=zero
             dvxcdgr(ipts,3)=zero
             !loop over the spin
             ispden=1
             rho   =rho_updn(ipts,ispden)
             rhomot=rho_updnm1_3(ipts,ispden)
             ex_lsd= - threefourth_divpi * sixpi2_1_3*rhomot*rhomot*rho
             !   Perdew-Burke-Ernzerhof GGA, exchange part
             rho_inv=rhomot*rhomot*rhomot
             coeffss=quarter*sixpi2m1_3*sixpi2m1_3*rho_inv*rho_inv*rhomot*rhomot
             ss=grho2_updn(ipts,ispden)*coeffss
             divss=one/(one+mu_divkappa*ss)
             dfxdss= mu*divss*divss
             d2fxdss2=-mu*two*mu_divkappa*divss*divss*divss
             fx    = one+kappa*(one-divss)
             ex_gga= ex_lsd*fx
             dssdn=-eight*third*ss*rho_inv
             dfxdn  = dfxdss*dssdn
             vxci(ipts,ispden)=ex_lsd*(four_thirds*fx+rho*dfxdn)
             !   The new definition (v3.3) includes the division by the norm of the gradient
             dssdg =two*coeffss
             dfxdg=dfxdss*dssdg
             dvxcdgr(ipts,ispden)=ex_lsd*rho*dfxdg
             exc=exc+ex_gga*rho

             !    Perdew-Burke-Ernzerhof GGA, exchange part
             !    Components 3 or 4
             dvxci(ipts,2+ispden)=dvxcdgr(ipts,ispden)
             !    Components 1 or 2
             d2ssdn2=-11.0_dp*third*dssdn*rho_inv
             d2fxdn2=d2fxdss2*dssdn**2+dfxdss*d2ssdn2
             dvxci(ipts,ispden)=third*rho_inv*vxci(ipts,ispden)+&
                  &     ex_lsd*(seven*third*dfxdn+rho*d2fxdn2)
             !    Components 5 or 6
             d2ssdndg=-eight*third*dssdg*rho_inv
             d2fxdndg=d2fxdss2*dssdn*dssdg+dfxdss*d2ssdndg
             dvxci(ipts,4+ispden)=ex_lsd*(four_thirds*dfxdg+rho*d2fxdndg)
             !    Components 7 or 8
             d2fxdg2=d2fxdss2*dssdg**2
             dvxci(ipts,6+ispden)=ex_lsd*rho*d2fxdg2
             !    For the time being, treat non-spin-polarized like spin-polarized
             dvxci(ipts,2)=dvxci(ipts,1)
             dvxci(ipts,4)=dvxci(ipts,3)
             dvxci(ipts,6)=dvxci(ipts,5)
             dvxci(ipts,8)=dvxci(ipts,7)
             !end of loop over the spin
             !  If non spin-polarized, treat spin down contribution now, similar to spin up
             exc=exc*2
             exci(ipts)=exc*rhotot_inv
             ! -----------------------------------------------------------------------------
             ! Then takes care of the LSD correlation part of the functional

             ! Correlation has been added
             ! -----------------------------------------------------------------------------

             !  vxci(ipts,2)=vxci(ipts,1)
             dvxcdgr(ipts,2)=dvxcdgr(ipts,1)

          end do
       else if (option==1) then
          do ipts=1,npts

             rhotot=rhoarr(ipts)
             rhotmot=rhom1_3(ipts)
             rhotot_inv=rhotmot*rhotmot*rhotmot
             rhotmo6=sqrt(rhotmot)
             rhoto6=rhotot*rhotmot*rhotmot*rhotmo6
             ! -----------------------------------------------------------------------
             ! First take care of the exchange part of the functional

             exc=zero
             !loop over the spin
             ispden=1
             rho   =rho_updn(ipts,ispden)
             rhomot=rho_updnm1_3(ipts,ispden)
             ex_lsd= - threefourth_divpi * sixpi2_1_3*rhomot*rhomot*rho
             !   Perdew_Wang 91 LSD
             vxci(ipts,ispden)=four_thirds*ex_lsd
             exc=exc+ex_lsd*rho

             !    Perdew_Wang 91 LSD
             dvxci(ipts,2*ispden-1)=-four_thirds*third*&
                  &                   threefourth_divpi*sixpi2_1_3*rhomot*rhomot
             dvxci(ipts,2)=zero
             !    If non-spin-polarized, first component of dvxci is second
             !    derivative with respect to TOTAL density.
             dvxci(ipts,1)=dvxci(ipts,1)*half
             ! Compute the second derivative of vx
             ! vx^(2) = -2*vx^(1)/(3*rhotot)
             !end of loop over the spin
             !  If non spin-polarized, treat spin down contribution now, similar to spin up
             exc=exc*2
             exci(ipts)=exc*rhotot_inv
             ! -----------------------------------------------------------------------------
             ! Then takes care of the LSD correlation part of the functional


             rs=rsfac*rhotmot
             sqr_rs=sq_rsfac*rhotmo6
             rsm1_2=sq_rsfac_inv*rhoto6

             !  Formulas A6-A8 of PW92LSD
             ec0_q0=-2.0_dp*ec0_aa*(1.0_dp+ec0_a1*rs)
             ec0_q1=2.0_dp*ec0_aa*(ec0_b1*sqr_rs+ec0_b2*rs+ec0_b3*rs*sqr_rs+ec0_b4*rs*rs)
             ec0_q1p=ec0_aa*(ec0_b1*rsm1_2+2._dp*ec0_b2+3._dp*ec0_b3*sqr_rs+4._dp*ec0_b4*rs)
             ec0_den=1.0_dp/(ec0_q1*ec0_q1+ec0_q1)
             !  ec0_log=log( 1.0_dp + 1.0_dp / ec0_q1 )
             ec0_log=-log( ec0_q1*ec0_q1*ec0_den )
             ecrs0=ec0_q0*ec0_log
             decrs0_drs= -2.0_dp*ec0_aa*ec0_a1*ec0_log - ec0_q0*ec0_q1p *ec0_den
             ec0_q1pp=0.5_dp*ec0_aa*(-ec0_b1*rsm1_2**3+3._dp*ec0_b3*rsm1_2+8._dp*ec0_b4)
             d2ecrs0_drs2= 4.0_dp*ec0_aa*ec0_a1*ec0_q1p*ec0_den            &
                  &                -ec0_q0*ec0_q1pp*ec0_den                        &
                  &                +ec0_q0*ec0_q1p**2*ec0_den**2*(2._dp*ec0_q1+1.0_dp)
             mac_q0=-2.0_dp*mac_aa*(1.0_dp+mac_a1*rs)
             mac_q1=2.0_dp*mac_aa*(mac_b1*sqr_rs+mac_b2*rs+mac_b3*rs*sqr_rs+mac_b4*rs*rs)
             mac_q1p=mac_aa*(mac_b1*rsm1_2+2._dp*mac_b2+3._dp*mac_b3*sqr_rs+4._dp*mac_b4*rs)
             mac_den=1.0_dp/(mac_q1*mac_q1+mac_q1)
             mac_log=-log( mac_q1*mac_q1*mac_den )
             macrs=mac_q0*mac_log
             dmacrs_drs= -2.0_dp*mac_aa*mac_a1*mac_log - mac_q0*mac_q1p*mac_den
             ecrs=ecrs0
             decrs_drs=decrs0_drs
             decrs_dzeta=0.0_dp
             d2ecrs_drs2=d2ecrs0_drs2
             d2ecrs_dzeta2=alpha_zeta**2*(-macrs)
             d2ecrs_drsdzeta=zero
             zeta=0.0_dp

             !  Add LSD correlation functional to GGA exchange functional
             exci(ipts)=exci(ipts)+ecrs
             vxci(ipts,1)=vxci(ipts,1)+ecrs-rs*third*decrs_drs

             dvcrs_drs=third*(2._dp*decrs_drs-rs*d2ecrs_drs2)
             !   And d(vxc)/d(rho)=(-rs/(3*rho))*d(vxc)/d(rs)
             d2ecrs_drho2= -rs**4*(four_pi*third)*third*dvcrs_drs

             dvxci(ipts,1)=dvxci(ipts,1)+d2ecrs_drho2

          end do
       else if (option==3) then
          do ipts=1,npts

             rhotot=rhoarr(ipts)
             rhotmot=rhom1_3(ipts)
             rhotot_inv=rhotmot*rhotmot*rhotmot
             rhotmo6=sqrt(rhotmot)
             rhoto6=rhotot*rhotmot*rhotmot*rhotmo6
             ! -----------------------------------------------------------------------
             ! First take care of the exchange part of the functional

             exc=zero
             !loop over the spin
             ispden=1
             rho   =rho_updn(ipts,ispden)
             rhomot=rho_updnm1_3(ipts,ispden)
             ex_lsd= - threefourth_divpi * sixpi2_1_3*rhomot*rhomot*rho
             !   Perdew_Wang 91 LSD
             vxci(ipts,ispden)=four_thirds*ex_lsd
             exc=exc+ex_lsd*rho

             !    Perdew_Wang 91 LSD
             dvxci(ipts,2*ispden-1)=-four_thirds*third*&
                  &                   threefourth_divpi*sixpi2_1_3*rhomot*rhomot
             dvxci(ipts,2)=zero
             !    If non-spin-polarized, first component of dvxci is second
             !    derivative with respect to TOTAL density.
             dvxci(ipts,1)=dvxci(ipts,1)*half
             ! Compute the second derivative of vx
             ! vx^(2) = -2*vx^(1)/(3*rhotot)
             !end of loop over the spin
             !  If non spin-polarized, treat spin down contribution now, similar to spin up
             exc=exc*2
             exci(ipts)=exc*rhotot_inv
             ! -----------------------------------------------------------------------------
             ! Then takes care of the LSD correlation part of the functional


             rs=rsfac*rhotmot
             sqr_rs=sq_rsfac*rhotmo6
             rsm1_2=sq_rsfac_inv*rhoto6

             !  Formulas A6-A8 of PW92LSD
             ec0_q0=-2.0_dp*ec0_aa*(1.0_dp+ec0_a1*rs)
             sqr_sqr_rs=max(1.e-15_dp,sqrt(sqr_rs))
             ec0_q1=2.0_dp*ec0_aa*(ec0_b1*sqr_rs+ec0_b2*rs+ec0_b3*rs*sqr_rs+ec0_b4*rs*rs/sqr_sqr_rs)
             ec0_q1p=ec0_aa*(ec0_b1*rsm1_2+2._dp*ec0_b2+3._dp*ec0_b3*sqr_rs+3.5_dp*ec0_b4*rs/sqr_sqr_rs)
             ec0_den=1.0_dp/(ec0_q1*ec0_q1+ec0_q1)
             !  ec0_log=log( 1.0_dp + 1.0_dp / ec0_q1 )
             ec0_log=-log( ec0_q1*ec0_q1*ec0_den )
             ecrs0=ec0_q0*ec0_log
             decrs0_drs= -2.0_dp*ec0_aa*ec0_a1*ec0_log - ec0_q0*ec0_q1p *ec0_den
             ec0_q1pp=0.5_dp*ec0_aa*(-ec0_b1*rsm1_2**3+3._dp*ec0_b3*rsm1_2+8._dp*ec0_b4)
             d2ecrs0_drs2= 4.0_dp*ec0_aa*ec0_a1*ec0_q1p*ec0_den            &
                  &                -ec0_q0*ec0_q1pp*ec0_den                        &
                  &                +ec0_q0*ec0_q1p**2*ec0_den**2*(2._dp*ec0_q1+1.0_dp)
             mac_q0=-2.0_dp*mac_aa*(1.0_dp+mac_a1*rs)
             mac_q1=2.0_dp*mac_aa*(mac_b1*sqr_rs+mac_b2*rs+mac_b3*rs*sqr_rs+mac_b4*rs*rs)
             mac_q1p=mac_aa*(mac_b1*rsm1_2+2._dp*mac_b2+3._dp*mac_b3*sqr_rs+4._dp*mac_b4*rs)
             mac_den=1.0_dp/(mac_q1*mac_q1+mac_q1)
             mac_log=-log( mac_q1*mac_q1*mac_den )
             macrs=mac_q0*mac_log
             dmacrs_drs= -2.0_dp*mac_aa*mac_a1*mac_log - mac_q0*mac_q1p*mac_den
             ecrs=ecrs0
             decrs_drs=decrs0_drs
             decrs_dzeta=0.0_dp
             d2ecrs_drs2=d2ecrs0_drs2
             d2ecrs_dzeta2=alpha_zeta**2*(-macrs)
             d2ecrs_drsdzeta=zero
             zeta=0.0_dp

             !  Add LSD correlation functional to GGA exchange functional
             exci(ipts)=exci(ipts)+ecrs
             vxci(ipts,1)=vxci(ipts,1)+ecrs-rs*third*decrs_drs

             dvcrs_drs=third*(2._dp*decrs_drs-rs*d2ecrs_drs2)
             !   And d(vxc)/d(rho)=(-rs/(3*rho))*d(vxc)/d(rs)
             d2ecrs_drho2= -rs**4*(four_pi*third)*third*dvcrs_drs
             dvxci(ipts,1)=dvxci(ipts,1)+d2ecrs_drho2

          end do
       end if

    end if
 else if(nspden==2) then

    !we separate different cases depending on order

    if (order**2<=1) then

       do ipts=1,npts

          rhotot=rhoarr(ipts)
          rhotmot=rhom1_3(ipts)
          rhotot_inv=rhotmot*rhotmot*rhotmot
          rhotmo6=sqrt(rhotmot)
          rhoto6=rhotot*rhotmot*rhotmot*rhotmo6
          ! -----------------------------------------------------------------------
          ! First take care of the exchange part of the functional

          exc=zero
          if (present(dvxcdgr)) dvxcdgr(ipts,3)=zero
          do ispden=1,nspden
             rho   =rho_updn(ipts,ispden)
             rhomot=rho_updnm1_3(ipts,ispden)
             ex_lsd= - threefourth_divpi * sixpi2_1_3*rhomot*rhomot*rho
             if(option==1 .or. option==-1 .or. option==3)then
                !   Perdew_Wang 91 LSD
                vxci(ipts,ispden)=four_thirds*ex_lsd
                if(present(dvxcdgr)) dvxcdgr(ipts,ispden)=0.0_dp
                exc=exc+ex_lsd*rho
             else
                !   Perdew-Burke-Ernzerhof GGA, exchange part
                rho_inv=rhomot*rhomot*rhomot
                coeffss=quarter*sixpi2m1_3*sixpi2m1_3*rho_inv*rho_inv*rhomot*rhomot
                ss=grho2_updn(ipts,ispden)*coeffss
                if(option/=6)then
                   divss=one/(one+mu_divkappa*ss)
                   dfxdss= mu*divss*divss
                   d2fxdss2=-mu*two*mu_divkappa*divss*divss*divss
                else  ! This is RPBE modification
                   divss=exp(-mu_divkappa*ss)
                   dfxdss= mu*divss
                   d2fxdss2=-mu*mu_divkappa*divss
                end if
                fx    = one+kappa*(one-divss)
                ex_gga= ex_lsd*fx
                dssdn=-eight*third*ss*rho_inv
                dfxdn  = dfxdss*dssdn
                vxci(ipts,ispden)=ex_lsd*(four_thirds*fx+rho*dfxdn)
                !   The new definition (v3.3) includes the division by the norm of the gradient
                dssdg =two*coeffss
                dfxdg=dfxdss*dssdg
                dvxcdgr(ipts,ispden)=ex_lsd*rho*dfxdg
                exc=exc+ex_gga*rho
             end if

          end do
          exci(ipts)=exc*rhotot_inv
          if(present(exexch).and.exexch==1)cycle

          ! -----------------------------------------------------------------------------
          ! Then takes care of the LSD correlation part of the functional

          if(option>0)then

             rs=rsfac*rhotmot
             sqr_rs=sq_rsfac*rhotmo6
             rsm1_2=sq_rsfac_inv*rhoto6

             !  Formulas A6-A8 of PW92LSD
             ec0_q0=-2.0_dp*ec0_aa*(1.0_dp+ec0_a1*rs)
             if(option/=3 .and. option/=4)then
                ec0_q1=2.0_dp*ec0_aa*(ec0_b1*sqr_rs+ec0_b2*rs+ec0_b3*rs*sqr_rs+ec0_b4*rs*rs)
                ec0_q1p=ec0_aa*(ec0_b1*rsm1_2+2._dp*ec0_b2+3._dp*ec0_b3*sqr_rs+4._dp*ec0_b4*rs)
             else
                sqr_sqr_rs=max(1.e-15_dp,sqrt(sqr_rs))
                ec0_q1=2.0_dp*ec0_aa*(ec0_b1*sqr_rs+ec0_b2*rs+ec0_b3*rs*sqr_rs+ec0_b4*rs*rs/sqr_sqr_rs)
                ec0_q1p=ec0_aa*(ec0_b1*rsm1_2+2._dp*ec0_b2+3._dp*ec0_b3*sqr_rs+3.5_dp*ec0_b4*rs/sqr_sqr_rs)
             end if
             ec0_den=1.0_dp/(ec0_q1*ec0_q1+ec0_q1)
             !  ec0_log=log( 1.0_dp + 1.0_dp / ec0_q1 )
             ec0_log=-log( ec0_q1*ec0_q1*ec0_den )
             ecrs0=ec0_q0*ec0_log
             decrs0_drs= -2.0_dp*ec0_aa*ec0_a1*ec0_log - ec0_q0*ec0_q1p *ec0_den

             mac_q0=-2.0_dp*mac_aa*(1.0_dp+mac_a1*rs)
             mac_q1=2.0_dp*mac_aa*(mac_b1*sqr_rs+mac_b2*rs+mac_b3*rs*sqr_rs+mac_b4*rs*rs)
             mac_q1p=mac_aa*(mac_b1*rsm1_2+2._dp*mac_b2+3._dp*mac_b3*sqr_rs+4._dp*mac_b4*rs)
             mac_den=1.0_dp/(mac_q1*mac_q1+mac_q1)
             mac_log=-log( mac_q1*mac_q1*mac_den )
             macrs=mac_q0*mac_log
             dmacrs_drs= -2.0_dp*mac_aa*mac_a1*mac_log - mac_q0*mac_q1p*mac_den

             zeta=(rho_updn(ipts,1)-rho_updn(ipts,2))*rhotot_inv
             ec1_q0=-2.0_dp*ec1_aa*(1.0_dp+ec1_a1*rs)
             if(option/=3 .and. option/=4)then
                ec1_q1=2.0_dp*ec1_aa*(ec1_b1*sqr_rs+ec1_b2*rs+ec1_b3*rs*sqr_rs+ec1_b4*rs*rs)
                ec1_q1p=ec1_aa*(ec1_b1*rsm1_2+2._dp*ec1_b2+3._dp*ec1_b3*sqr_rs+4._dp*ec1_b4*rs)
             else
                ec1_q1=2.0_dp*ec1_aa*(ec1_b1*sqr_rs+ec1_b2*rs+ec1_b3*rs*sqr_rs+ec1_b4*rs*rs/sqr_sqr_rs)
                ec1_q1p=ec1_aa*(ec1_b1*rsm1_2+2._dp*ec1_b2+3._dp*ec1_b3*sqr_rs+3.5_dp*ec1_b4*rs/sqr_sqr_rs)
             end if
             ec1_den=1.0_dp/(ec1_q1*ec1_q1+ec1_q1)
             !   ec1_log=log( 1.0_dp + 1.0_dp / ec1_q1 )
             ec1_log=-log( ec1_q1*ec1_q1*ec1_den )
             ecrs1=ec1_q0*ec1_log
             decrs1_drs= -2.0_dp*ec1_aa*ec1_a1*ec1_log - ec1_q0*ec1_q1p *ec1_den

             !   alpha_zeta is introduced in order to remove singularities for fully
             !   polarized systems.
             zetp_1_3=(1.0_dp+zeta*alpha_zeta)*zetpm1_3(ipts)**2
             zetm_1_3=(1.0_dp-zeta*alpha_zeta)*zetmm1_3(ipts)**2

             f_zeta=( (1.0_dp+zeta*alpha_zeta2)*zetp_1_3 +                      &
                  &            (1.0_dp-zeta*alpha_zeta2)*zetm_1_3 - 2.0_dp ) * factf_zeta
             fp_zeta=( zetp_1_3 - zetm_1_3 ) * factfp_zeta
             zeta4=zeta**4

             gcrs=ecrs1-ecrs0+macrs*fsec_inv
             !   ecrs=ecrs0+f_zeta*(-macrs*(1.0_dp-zeta4)*fsec_inv+(ecrs1-ecrs0)*zeta4)
             ecrs=ecrs0+f_zeta*(zeta4*gcrs-macrs*fsec_inv)

             dgcrs_drs=decrs1_drs-decrs0_drs+dmacrs_drs*fsec_inv
             !    decrs_drs=decrs0_drs+f_zeta*&
             !&        (-dmacrs_drs*(1.0_dp-zeta4)*fsec_inv+(decrs1_drs-decrs0_drs)*zeta4)
             decrs_drs=decrs0_drs+f_zeta*(zeta4*dgcrs_drs-dmacrs_drs*fsec_inv)
             dfzeta4_dzeta=4.0_dp*zeta**3*f_zeta+fp_zeta*zeta4
             decrs_dzeta=dfzeta4_dzeta*gcrs-fp_zeta*macrs*fsec_inv

             !  Add LSD correlation functional to GGA exchange functional
             exci(ipts)=exci(ipts)+ecrs
             vxcadd=ecrs-rs*third*decrs_drs-zeta*decrs_dzeta
             vxci(ipts,1)=vxci(ipts,1)+vxcadd+decrs_dzeta
             vxci(ipts,2)=vxci(ipts,2)+vxcadd-decrs_dzeta

             !  -----------------------------------------------------------------------------
             !  Eventually add the GGA correlation part of the PBE functional
             !  Note : the computation of the potential in the spin-unpolarized
             !  case could be optimized much further. Other optimizations are left to do.

             if(option==2 .or. option==5 .or. option==6 .or. option==7)then
                !    The definition of phi has been slightly changed, because
                !    the original PBE one gives divergent behaviour for fully
                !    polarized points
                !    zetpm1_3=(1.0_dp+zeta*alpha_zeta)**(-third)
                !    zetmm1_3=(1.0_dp-zeta*alpha_zeta)**(-third)
                phi_zeta=( zetpm1_3(ipts)*(1.0_dp+zeta*alpha_zeta)+ &
                     &               zetmm1_3(ipts)*(1.0_dp-zeta*alpha_zeta)   )*0.5_dp
                phip_zeta=(zetpm1_3(ipts)-zetmm1_3(ipts))*third*alpha_zeta
                phi_zeta_inv=1.0_dp/phi_zeta
                phi_logder=phip_zeta*phi_zeta_inv
                phi3_zeta=phi_zeta*phi_zeta*phi_zeta
                gamphi3inv=gamma_inv*phi_zeta_inv*phi_zeta_inv*phi_zeta_inv

                !   From ec to bb
                bb=ecrs*gamphi3inv
                dbb_drs=decrs_drs*gamphi3inv
                dbb_dzeta=gamphi3inv*(decrs_dzeta-three*ecrs*phi_logder)
                !   From bb to cc
                exp_pbe=exp(-bb)
                cc=one/(exp_pbe-one)
                dcc_dbb=cc*cc*exp_pbe
                dcc_drs=dcc_dbb*dbb_drs
                dcc_dzeta=dcc_dbb*dbb_dzeta

                !   From cc to aa
                coeff_aa=beta*gamma_inv*phi_zeta_inv*phi_zeta_inv
                aa=coeff_aa*cc
                daa_drs=coeff_aa*dcc_drs
                daa_dzeta=-two*aa*phi_logder+coeff_aa*dcc_dzeta
                !   Introduce tt : do not assume that the spin-dependent gradients are collinear
                grrho2=grho2_updn(ipts,3)
                dtt_dg=two*rhotot_inv*rhotot_inv*rhotmot*coeff_tt
                !   Note that tt is (the t variable of PBE divided by phi) squared
                tt=half*grrho2*dtt_dg

                !   Get xx from aa and tt
                xx=aa*tt
                dxx_drs=daa_drs*tt
                dxx_dzeta=daa_dzeta*tt
                dxx_dtt=aa
                !   From xx to pade
                pade_den=one/(one+xx*(one+xx))
                pade=(one+xx)*pade_den
                dpade_dxx=-xx*(two+xx)*pade_den**2
                dpade_drs=dpade_dxx*dxx_drs
                dpade_dtt=dpade_dxx*dxx_dtt
                dpade_dzeta=dpade_dxx*dxx_dzeta

                !   From pade to qq
                coeff_qq=tt*phi_zeta_inv*phi_zeta_inv
                qq=coeff_qq*pade
                dqq_drs=coeff_qq*dpade_drs
                dqq_dtt=pade*phi_zeta_inv*phi_zeta_inv+coeff_qq*dpade_dtt
                dqq_dzeta=coeff_qq*(dpade_dzeta-two*pade*phi_logder)

                !   From qq to rr
                arg_rr=one+beta*gamma_inv*qq
                div_rr=one/arg_rr
                rr=gamma*log(arg_rr)
                drr_dqq=beta*div_rr
                drr_drs=drr_dqq*dqq_drs
                drr_dtt=drr_dqq*dqq_dtt
                drr_dzeta=drr_dqq*dqq_dzeta

                !   From rr to hh
                hh=phi3_zeta*rr
                dhh_drs=phi3_zeta*drr_drs
                dhh_dtt=phi3_zeta*drr_dtt
                dhh_dzeta=phi3_zeta*(drr_dzeta+three*rr*phi_logder)

                !   The GGA correlation energy is added
                exci(ipts)=exci(ipts)+hh

                !   Change of variables : from (rs,zeta,tt) to (rhoup,rhodn,grrho)

                !   From hh to the derivative of the energy wrt the density
                drhohh_drho=hh-third*rs*dhh_drs-zeta*dhh_dzeta-seven*third*tt*dhh_dtt
                vxci(ipts,1)=vxci(ipts,1)+drhohh_drho+dhh_dzeta
                vxci(ipts,2)=vxci(ipts,2)+drhohh_drho-dhh_dzeta


                !   From hh to the derivative of the energy wrt to the gradient of the
                !   density, divided by the gradient of the density
                !   (The v3.3 definition includes the division by the norm of the gradient)
                dvxcdgr(ipts,3)=rhotot*dtt_dg*dhh_dtt

                !  End condition of GGA
             end if

          else  ! no correlation

             ! End condition of including correlation, and not only exchange
          end if

          ! Correlation has been added
          ! -----------------------------------------------------------------------------

       end do
    else

       do ipts=1,npts

          rhotot=rhoarr(ipts)
          rhotmot=rhom1_3(ipts)
          rhotot_inv=rhotmot*rhotmot*rhotmot
          rhotmo6=sqrt(rhotmot)
          rhoto6=rhotot*rhotmot*rhotmot*rhotmo6
          ! -----------------------------------------------------------------------
          ! First take care of the exchange part of the functional

          exc=zero
          if (present(dvxcdgr)) dvxcdgr(ipts,3)=zero
          do ispden=1,nspden
             rho   =rho_updn(ipts,ispden)
             rhomot=rho_updnm1_3(ipts,ispden)
             ex_lsd= - threefourth_divpi * sixpi2_1_3*rhomot*rhomot*rho
             if(option==1 .or. option==-1 .or. option==3)then
                !   Perdew_Wang 91 LSD
                vxci(ipts,ispden)=four_thirds*ex_lsd
                if(present(dvxcdgr)) dvxcdgr(ipts,ispden)=0.0_dp
                exc=exc+ex_lsd*rho
             else
                !   Perdew-Burke-Ernzerhof GGA, exchange part
                rho_inv=rhomot*rhomot*rhomot
                coeffss=quarter*sixpi2m1_3*sixpi2m1_3*rho_inv*rho_inv*rhomot*rhomot
                ss=grho2_updn(ipts,ispden)*coeffss
                if(option/=6)then
                   divss=one/(one+mu_divkappa*ss)
                   dfxdss= mu*divss*divss
                   d2fxdss2=-mu*two*mu_divkappa*divss*divss*divss
                else  ! This is RPBE modification
                   divss=exp(-mu_divkappa*ss)
                   dfxdss= mu*divss
                   d2fxdss2=-mu*mu_divkappa*divss
                end if
                fx    = one+kappa*(one-divss)
                ex_gga= ex_lsd*fx
                dssdn=-eight*third*ss*rho_inv
                dfxdn  = dfxdss*dssdn
                vxci(ipts,ispden)=ex_lsd*(four_thirds*fx+rho*dfxdn)
                !   The new definition (v3.3) includes the division by the norm of the gradient
                dssdg =two*coeffss
                dfxdg=dfxdss*dssdg
                dvxcdgr(ipts,ispden)=ex_lsd*rho*dfxdg
                exc=exc+ex_gga*rho
             end if

             if(option==1 .or. option==-1 .or. option==3)then
                !    Perdew_Wang 91 LSD
                dvxci(ipts,2*ispden-1)=-four_thirds*third*&
                     &                   threefourth_divpi*sixpi2_1_3*rhomot*rhomot
                dvxci(ipts,2)=zero
                !    If non-spin-polarized, first component of dvxci is second
                !    derivative with respect to TOTAL density.
                ! Compute the second derivative of vx
                ! vx^(2) = -2*vx^(1)/(3*rhotot)
             else
                !    Perdew-Burke-Ernzerhof GGA, exchange part
                !    Components 3 or 4
                dvxci(ipts,2+ispden)=dvxcdgr(ipts,ispden)
                !    Components 1 or 2
                d2ssdn2=-11.0_dp*third*dssdn*rho_inv
                d2fxdn2=d2fxdss2*dssdn**2+dfxdss*d2ssdn2
                dvxci(ipts,ispden)=third*rho_inv*vxci(ipts,ispden)+&
                     &     ex_lsd*(seven*third*dfxdn+rho*d2fxdn2)
                !    Components 5 or 6
                d2ssdndg=-eight*third*dssdg*rho_inv
                d2fxdndg=d2fxdss2*dssdn*dssdg+dfxdss*d2ssdndg
                dvxci(ipts,4+ispden)=ex_lsd*(four_thirds*dfxdg+rho*d2fxdndg)
                !    Components 7 or 8
                d2fxdg2=d2fxdss2*dssdg**2
                dvxci(ipts,6+ispden)=ex_lsd*rho*d2fxdg2
                !    For the time being, treat non-spin-polarized like spin-polarized
             end if
          end do
          exci(ipts)=exc*rhotot_inv
          ! -----------------------------------------------------------------------------
          ! Then takes care of the LSD correlation part of the functional

          if(option>0)then

             rs=rsfac*rhotmot
             sqr_rs=sq_rsfac*rhotmo6
             rsm1_2=sq_rsfac_inv*rhoto6

             !  Formulas A6-A8 of PW92LSD
             ec0_q0=-2.0_dp*ec0_aa*(1.0_dp+ec0_a1*rs)
             if(option/=3 .and. option/=4)then
                ec0_q1=2.0_dp*ec0_aa*(ec0_b1*sqr_rs+ec0_b2*rs+ec0_b3*rs*sqr_rs+ec0_b4*rs*rs)
                ec0_q1p=ec0_aa*(ec0_b1*rsm1_2+2._dp*ec0_b2+3._dp*ec0_b3*sqr_rs+4._dp*ec0_b4*rs)
             else
                sqr_sqr_rs=max(1.e-15_dp,sqrt(sqr_rs))
                ec0_q1=2.0_dp*ec0_aa*(ec0_b1*sqr_rs+ec0_b2*rs+ec0_b3*rs*sqr_rs+ec0_b4*rs*rs/sqr_sqr_rs)
                ec0_q1p=ec0_aa*(ec0_b1*rsm1_2+2._dp*ec0_b2+3._dp*ec0_b3*sqr_rs+3.5_dp*ec0_b4*rs/sqr_sqr_rs)
             end if
             ec0_den=1.0_dp/(ec0_q1*ec0_q1+ec0_q1)
             !  ec0_log=log( 1.0_dp + 1.0_dp / ec0_q1 )
             ec0_log=-log( ec0_q1*ec0_q1*ec0_den )
             ecrs0=ec0_q0*ec0_log
             decrs0_drs= -2.0_dp*ec0_aa*ec0_a1*ec0_log - ec0_q0*ec0_q1p *ec0_den
             ec0_q1pp=0.5_dp*ec0_aa*(-ec0_b1*rsm1_2**3+3._dp*ec0_b3*rsm1_2+8._dp*ec0_b4)
             d2ecrs0_drs2= 4.0_dp*ec0_aa*ec0_a1*ec0_q1p*ec0_den            &
                  &                -ec0_q0*ec0_q1pp*ec0_den                        &
                  &                +ec0_q0*ec0_q1p**2*ec0_den**2*(2._dp*ec0_q1+1.0_dp)
             if (order==3) then
                ec0_q1ppp = 0.75_dp*ec0_aa*(rsm1_2**5)*(ec0_b1-ec0_b3*rs)
                ec0_f1 = 1._dp/(ec0_q1*ec0_q1*(1._dp + ec0_q1))
                ec0_f2 = 1._dp/(ec0_q1*(1+ec0_q1))
                d3ecrs0_drs3 = 6._dp*ec0_q1p*ec0_f1*(-2._dp*ec0_aa*ec0_a1*ec0_q1p + &
                 &        ec0_q0*ec0_q1pp) - &
                 &        ec0_f2*(-6._dp*ec0_aa*ec0_a1*ec0_q1pp + ec0_q0*ec0_q1ppp + &
                 &        ec0_f2*(3._dp*ec0_q1p*(-2._dp*ec0_aa*ec0_a1*ec0_q1p + ec0_q0*ec0_q1pp) + &
                 &        ec0_f2*2._dp*ec0_q0*(ec0_q1p**3)*(1._dp + 3._dp*ec0_q1*(1._dp + ec0_q1))))
             end if

             mac_q0=-2.0_dp*mac_aa*(1.0_dp+mac_a1*rs)
             mac_q1=2.0_dp*mac_aa*(mac_b1*sqr_rs+mac_b2*rs+mac_b3*rs*sqr_rs+mac_b4*rs*rs)
             mac_q1p=mac_aa*(mac_b1*rsm1_2+2._dp*mac_b2+3._dp*mac_b3*sqr_rs+4._dp*mac_b4*rs)
             mac_den=1.0_dp/(mac_q1*mac_q1+mac_q1)
             mac_log=-log( mac_q1*mac_q1*mac_den )
             macrs=mac_q0*mac_log
             dmacrs_drs= -2.0_dp*mac_aa*mac_a1*mac_log - mac_q0*mac_q1p*mac_den

             zeta=(rho_updn(ipts,1)-rho_updn(ipts,2))*rhotot_inv
             ec1_q0=-2.0_dp*ec1_aa*(1.0_dp+ec1_a1*rs)
             if(option/=3 .and. option/=4)then
                ec1_q1=2.0_dp*ec1_aa*(ec1_b1*sqr_rs+ec1_b2*rs+ec1_b3*rs*sqr_rs+ec1_b4*rs*rs)
                ec1_q1p=ec1_aa*(ec1_b1*rsm1_2+2._dp*ec1_b2+3._dp*ec1_b3*sqr_rs+4._dp*ec1_b4*rs)
             else
                ec1_q1=2.0_dp*ec1_aa*(ec1_b1*sqr_rs+ec1_b2*rs+ec1_b3*rs*sqr_rs+ec1_b4*rs*rs/sqr_sqr_rs)
                ec1_q1p=ec1_aa*(ec1_b1*rsm1_2+2._dp*ec1_b2+3._dp*ec1_b3*sqr_rs+3.5_dp*ec1_b4*rs/sqr_sqr_rs)
             end if
             ec1_den=1.0_dp/(ec1_q1*ec1_q1+ec1_q1)
             !   ec1_log=log( 1.0_dp + 1.0_dp / ec1_q1 )
             ec1_log=-log( ec1_q1*ec1_q1*ec1_den )
             ecrs1=ec1_q0*ec1_log
             decrs1_drs= -2.0_dp*ec1_aa*ec1_a1*ec1_log - ec1_q0*ec1_q1p *ec1_den

             !   alpha_zeta is introduced in order to remove singularities for fully
             !   polarized systems.
             zetp_1_3=(1.0_dp+zeta*alpha_zeta)*zetpm1_3(ipts)**2
             zetm_1_3=(1.0_dp-zeta*alpha_zeta)*zetmm1_3(ipts)**2

             f_zeta=( (1.0_dp+zeta*alpha_zeta2)*zetp_1_3 +                      &
                  &            (1.0_dp-zeta*alpha_zeta2)*zetm_1_3 - 2.0_dp ) * factf_zeta
             fp_zeta=( zetp_1_3 - zetm_1_3 ) * factfp_zeta
             zeta4=zeta**4

             gcrs=ecrs1-ecrs0+macrs*fsec_inv
             !   ecrs=ecrs0+f_zeta*(-macrs*(1.0_dp-zeta4)*fsec_inv+(ecrs1-ecrs0)*zeta4)
             ecrs=ecrs0+f_zeta*(zeta4*gcrs-macrs*fsec_inv)

             dgcrs_drs=decrs1_drs-decrs0_drs+dmacrs_drs*fsec_inv
             !    decrs_drs=decrs0_drs+f_zeta*&
             !&        (-dmacrs_drs*(1.0_dp-zeta4)*fsec_inv+(decrs1_drs-decrs0_drs)*zeta4)
             decrs_drs=decrs0_drs+f_zeta*(zeta4*dgcrs_drs-dmacrs_drs*fsec_inv)
             dfzeta4_dzeta=4.0_dp*zeta**3*f_zeta+fp_zeta*zeta4
             decrs_dzeta=dfzeta4_dzeta*gcrs-fp_zeta*macrs*fsec_inv

             ec1_q1pp=0.5_dp*ec1_aa*(-ec1_b1*rsm1_2**3+3._dp*ec1_b3*rsm1_2+8._dp*ec1_b4)
             d2ecrs1_drs2= 4.0_dp*ec1_aa*ec1_a1*ec1_q1p*ec1_den            &
                  &                 -ec1_q0*ec1_q1pp*ec1_den                        &
                  &                 +ec1_q0*ec1_q1p**2*ec1_den**2*(2._dp*ec1_q1+1.0_dp)

             mac_q1pp=0.5_dp*mac_aa*(-mac_b1*rsm1_2**3+3._dp*mac_b3*rsm1_2+8._dp*mac_b4)
             d2macrs_drs2= 4.0_dp*mac_aa*mac_a1*mac_q1p*mac_den            &
                  &                 -mac_q0*mac_q1pp*mac_den                        &
                  &                 +mac_q0*mac_q1p**2*mac_den**2*(2._dp*mac_q1+1.0_dp)

             d2gcrs_drs2=d2ecrs1_drs2-d2ecrs0_drs2+d2macrs_drs2*fsec_inv
             fpp_zeta=(zetpm1_3(ipts)**2+zetmm1_3(ipts)**2) * factfpp_zeta
             d2fzeta4_dzeta2=12.0_dp*zeta**2*f_zeta  &
                  &                   + 8.0_dp*zeta**3*fp_zeta &
                  &                   +       zeta4  *fpp_zeta

             d2ecrs_drs2=d2ecrs0_drs2+&
                  &        f_zeta*(zeta4*d2gcrs_drs2-d2macrs_drs2*fsec_inv)
             d2ecrs_drsdzeta=dfzeta4_dzeta*dgcrs_drs-fp_zeta*dmacrs_drs*fsec_inv
             d2ecrs_dzeta2=d2fzeta4_dzeta2*gcrs-fpp_zeta*macrs*fsec_inv

             !   End condition of abs(order)>1
             !  Add LSD correlation functional to GGA exchange functional
             exci(ipts)=exci(ipts)+ecrs
             vxcadd=ecrs-rs*third*decrs_drs-zeta*decrs_dzeta
             vxci(ipts,1)=vxci(ipts,1)+vxcadd+decrs_dzeta
             vxci(ipts,2)=vxci(ipts,2)+vxcadd-decrs_dzeta

             dvcrs_drs=third*(2._dp*decrs_drs-rs*d2ecrs_drs2)
             !   And d(vxc)/d(rho)=(-rs/(3*rho))*d(vxc)/d(rs)
             d2ecrs_drho2= -rs**4*(four_pi*third)*third*dvcrs_drs
             d2ecrs_drup2=d2ecrs_drho2+&
                  &       two*(-third*rs*d2ecrs_drsdzeta)*(1._dp-zeta)*rhotot_inv+ &
                  &       d2ecrs_dzeta2*(1._dp-zeta)**2*rhotot_inv
             d2ecrs_drdndrup=d2ecrs_drho2+&
                  &       2.0_dp*(-third*rs*d2ecrs_drsdzeta)*(-zeta)*rhotot_inv+ &
                  &       d2ecrs_dzeta2*(1._dp-zeta)*(-1._dp-zeta)*rhotot_inv
             d2ecrs_drdn2=d2ecrs_drho2+&
                  &       2.0_dp*(-third*rs*d2ecrs_drsdzeta)*(-1._dp-zeta)*rhotot_inv+ &
                  &       d2ecrs_dzeta2*(-1._dp-zeta)**2*rhotot_inv
             if(option==1 .or. option==-1 .or. option==3)then
                dvxci(ipts,1)=dvxci(ipts,1)+d2ecrs_drup2
                dvxci(ipts,2)=dvxci(ipts,2)+d2ecrs_drdndrup
                dvxci(ipts,3)=dvxci(ipts,3)+d2ecrs_drdn2
             else
                dvxci(ipts,9)=d2ecrs_drup2
                dvxci(ipts,10)=d2ecrs_drdndrup
                dvxci(ipts,11)=d2ecrs_drdn2
             end if

             !  -----------------------------------------------------------------------------
             !  Eventually add the GGA correlation part of the PBE functional
             !  Note : the computation of the potential in the spin-unpolarized
             !  case could be optimized much further. Other optimizations are left to do.

             if(option==2 .or. option==5 .or. option==6 .or. option==7)then
                !    The definition of phi has been slightly changed, because
                !    the original PBE one gives divergent behaviour for fully
                !    polarized points
                !    zetpm1_3=(1.0_dp+zeta*alpha_zeta)**(-third)
                !    zetmm1_3=(1.0_dp-zeta*alpha_zeta)**(-third)
                phi_zeta=( zetpm1_3(ipts)*(1.0_dp+zeta*alpha_zeta)+ &
                     &               zetmm1_3(ipts)*(1.0_dp-zeta*alpha_zeta)   )*0.5_dp
                phip_zeta=(zetpm1_3(ipts)-zetmm1_3(ipts))*third*alpha_zeta
                phi_zeta_inv=1.0_dp/phi_zeta
                phi_logder=phip_zeta*phi_zeta_inv
                phi3_zeta=phi_zeta*phi_zeta*phi_zeta
                gamphi3inv=gamma_inv*phi_zeta_inv*phi_zeta_inv*phi_zeta_inv
                phipp_zeta=-alpha_zeta*alpha_zeta*ninth*&
                     &      (zetpm1_3(ipts)*zetpm1_3(ipts)*zetpm1_3(ipts)*zetpm1_3(ipts) + &
                     &       zetmm1_3(ipts)*zetmm1_3(ipts)*zetmm1_3(ipts)*zetmm1_3(ipts)  )

                !   From ec to bb
                bb=ecrs*gamphi3inv
                dbb_drs=decrs_drs*gamphi3inv
                dbb_dzeta=gamphi3inv*(decrs_dzeta-three*ecrs*phi_logder)
                d2bb_drs2=d2ecrs_drs2*gamphi3inv
                d2bb_drsdzeta=gamphi3inv*(d2ecrs_drsdzeta-three*decrs_drs*phi_logder)
                d2bb_dzeta2=gamphi3inv*(d2ecrs_dzeta2-six*decrs_dzeta*phi_logder+&
                     &     12.0_dp*ecrs*phi_logder*phi_logder-three*ecrs*phi_zeta_inv*phipp_zeta)

                !   From bb to cc
                exp_pbe=exp(-bb)
                cc=one/(exp_pbe-one)
                dcc_dbb=cc*cc*exp_pbe
                dcc_drs=dcc_dbb*dbb_drs
                dcc_dzeta=dcc_dbb*dbb_dzeta
                d2cc_dbb2=cc*cc*exp_pbe*(two*cc*exp_pbe-one)
                d2cc_drs2=d2cc_dbb2*dbb_drs*dbb_drs+dcc_dbb*d2bb_drs2
                d2cc_drsdzeta=d2cc_dbb2*dbb_drs*dbb_dzeta+dcc_dbb*d2bb_drsdzeta
                d2cc_dzeta2=d2cc_dbb2*dbb_dzeta*dbb_dzeta+dcc_dbb*d2bb_dzeta2

                !   From cc to aa
                coeff_aa=beta*gamma_inv*phi_zeta_inv*phi_zeta_inv
                aa=coeff_aa*cc
                daa_drs=coeff_aa*dcc_drs
                daa_dzeta=-two*aa*phi_logder+coeff_aa*dcc_dzeta
                d2aa_drs2=coeff_aa*d2cc_drs2
                d2aa_drsdzeta=-two*daa_drs*phi_logder+coeff_aa*d2cc_drsdzeta
                d2aa_dzeta2=aa*(-two*phi_zeta_inv*phipp_zeta+six*phi_logder*phi_logder)+&
                     &     coeff_aa*(-four*dcc_dzeta*phi_logder+d2cc_dzeta2)

                !   Introduce tt : do not assume that the spin-dependent gradients are collinear
                grrho2=grho2_updn(ipts,3)
                dtt_dg=two*rhotot_inv*rhotot_inv*rhotmot*coeff_tt
                !   Note that tt is (the t variable of PBE divided by phi) squared
                tt=half*grrho2*dtt_dg

                !   Get xx from aa and tt
                xx=aa*tt
                dxx_drs=daa_drs*tt
                dxx_dzeta=daa_dzeta*tt
                dxx_dtt=aa
                d2xx_drs2=d2aa_drs2*tt
                d2xx_drsdzeta=d2aa_drsdzeta*tt
                d2xx_drsdtt=daa_drs
                d2xx_dttdzeta=daa_dzeta
                d2xx_dzeta2=d2aa_dzeta2*tt

                !   From xx to pade
                pade_den=one/(one+xx*(one+xx))
                pade=(one+xx)*pade_den
                dpade_dxx=-xx*(two+xx)*pade_den**2
                dpade_drs=dpade_dxx*dxx_drs
                dpade_dtt=dpade_dxx*dxx_dtt
                dpade_dzeta=dpade_dxx*dxx_dzeta
                d2pade_dxx2=two*(-one+xx*xx*(three+xx))*pade_den*pade_den*pade_den
                d2pade_drs2=d2pade_dxx2*dxx_drs*dxx_drs+dpade_dxx*d2xx_drs2
                d2pade_drsdtt=d2pade_dxx2*dxx_drs*dxx_dtt+dpade_dxx*d2xx_drsdtt
                d2pade_drsdzeta=d2pade_dxx2*dxx_drs*dxx_dzeta+dpade_dxx*d2xx_drsdzeta
                d2pade_dtt2=d2pade_dxx2*dxx_dtt*dxx_dtt
                d2pade_dttdzeta=d2pade_dxx2*dxx_dtt*dxx_dzeta+dpade_dxx*d2xx_dttdzeta
                d2pade_dzeta2=d2pade_dxx2*dxx_dzeta*dxx_dzeta+dpade_dxx*d2xx_dzeta2

                !   From pade to qq
                coeff_qq=tt*phi_zeta_inv*phi_zeta_inv
                qq=coeff_qq*pade
                dqq_drs=coeff_qq*dpade_drs
                dqq_dtt=pade*phi_zeta_inv*phi_zeta_inv+coeff_qq*dpade_dtt
                dqq_dzeta=coeff_qq*(dpade_dzeta-two*pade*phi_logder)
                d2qq_drs2=coeff_qq*d2pade_drs2
                d2qq_drsdtt=phi_zeta_inv*phi_zeta_inv*(dpade_drs+tt*d2pade_drsdtt)
                d2qq_drsdzeta=coeff_qq*(d2pade_drsdzeta-two*dpade_drs*phi_logder)
                d2qq_dtt2=phi_zeta_inv*phi_zeta_inv*(two*dpade_dtt+tt*d2pade_dtt2)
                d2qq_dttdzeta=phi_zeta_inv*phi_zeta_inv*(dpade_dzeta-two*pade*phi_logder)+&
                     &     coeff_qq*(d2pade_dttdzeta-two*dpade_dtt*phi_logder)
                d2qq_dzeta2=coeff_qq*( d2pade_dzeta2-four*dpade_dzeta*phi_logder &
                     &                          +six*pade*phi_logder*phi_logder            &
                     &                          -two*pade*phi_zeta_inv*phipp_zeta)

                !   From qq to rr
                arg_rr=one+beta*gamma_inv*qq
                div_rr=one/arg_rr
                rr=gamma*log(arg_rr)
                drr_dqq=beta*div_rr
                drr_drs=drr_dqq*dqq_drs
                drr_dtt=drr_dqq*dqq_dtt
                drr_dzeta=drr_dqq*dqq_dzeta
                d2rr_dqq2=-div_rr**2*beta*beta*gamma_inv
                d2rr_drs2=d2rr_dqq2*dqq_drs*dqq_drs+drr_dqq*d2qq_drs2
                d2rr_drsdtt=d2rr_dqq2*dqq_drs*dqq_dtt+drr_dqq*d2qq_drsdtt
                d2rr_drsdzeta=d2rr_dqq2*dqq_drs*dqq_dzeta+drr_dqq*d2qq_drsdzeta
                d2rr_dtt2=d2rr_dqq2*dqq_dtt*dqq_dtt+drr_dqq*d2qq_dtt2
                d2rr_dttdzeta=d2rr_dqq2*dqq_dtt*dqq_dzeta+drr_dqq*d2qq_dttdzeta
                d2rr_dzeta2=d2rr_dqq2*dqq_dzeta*dqq_dzeta+drr_dqq*d2qq_dzeta2

                !   From rr to hh
                hh=phi3_zeta*rr
                dhh_drs=phi3_zeta*drr_drs
                dhh_dtt=phi3_zeta*drr_dtt
                dhh_dzeta=phi3_zeta*(drr_dzeta+three*rr*phi_logder)
                d2hh_drs2=phi3_zeta*d2rr_drs2
                d2hh_drsdtt=phi3_zeta*d2rr_drsdtt
                d2hh_drsdzeta=phi3_zeta*(d2rr_drsdzeta+three*drr_drs*phi_logder)
                d2hh_dtt2=phi3_zeta*d2rr_dtt2
                d2hh_dttdzeta=phi3_zeta*(d2rr_dttdzeta+three*drr_dtt*phi_logder)
                d2hh_dzeta2=phi3_zeta*(six*rr*phi_logder*phi_logder+&
                     &                           six*phi_logder*drr_dzeta+d2rr_dzeta2)  &
                     &               +three*phi_zeta*phi_zeta*rr*phipp_zeta

                !   The GGA correlation energy is added
                exci(ipts)=exci(ipts)+hh

                !   Change of variables : from (rs,zeta,tt) to (rhoup,rhodn,grrho)


                !   From hh to the derivative of the energy wrt the density
                drhohh_drho=hh-third*rs*dhh_drs-zeta*dhh_dzeta-seven*third*tt*dhh_dtt
                vxci(ipts,1)=vxci(ipts,1)+drhohh_drho+dhh_dzeta
                vxci(ipts,2)=vxci(ipts,2)+drhohh_drho-dhh_dzeta


                !   From hh to the derivative of the energy wrt to the gradient of the
                !   density, divided by the gradient of the density
                !   (The v3.3 definition includes the division by the norm of the gradient)
                dvxcdgr(ipts,3)=rhotot*dtt_dg*dhh_dtt

                d2rhohh_drho2=rhotot_inv*&
                     &     (-two*ninth*rs*dhh_drs +seven*four*ninth*tt*dhh_dtt &
                     &      +ninth*rs*rs*d2hh_drs2+zeta*zeta*d2hh_dzeta2+(seven*third*tt)**2*d2hh_dtt2 &
                     &      +two*third*rs*zeta*d2hh_drsdzeta+two*seven*ninth*rs*tt*d2hh_drsdtt &
                     &      +two*seven*third*tt*zeta*d2hh_dttdzeta)
                d2rhohh_drhodg=dtt_dg*(-four*third*dhh_dtt-third*rs*d2hh_drsdtt &
                     &                           -zeta*d2hh_dttdzeta-seven*third*tt*d2hh_dtt2)

                !    Component 12 : first derivative with respect to the gradient
                !    of the density, div by the grad of the density
                dvxci(ipts,12)=dvxcdgr(ipts,3)
                !    Components 9, 10 and 11 : second derivatives with respect to the spin-density
                !    Note that there is already a contribution from LSDA
                dvxci(ipts,9)=dvxci(ipts,9)+d2rhohh_drho2+rhotot_inv*           &
                     &     ( d2hh_dzeta2*(one-two*zeta) &
                     &      -two*third*rs*d2hh_drsdzeta-14.0_dp*third*tt*d2hh_dttdzeta)
                dvxci(ipts,10)=dvxci(ipts,10)+d2rhohh_drho2-rhotot_inv*d2hh_dzeta2
                dvxci(ipts,11)=dvxci(ipts,11)+d2rhohh_drho2+rhotot_inv*           &
                     &     ( d2hh_dzeta2*(one+two*zeta) &
                     &      +two*third*rs*d2hh_drsdzeta+14.0_dp*third*tt*d2hh_dttdzeta)
                !    Components 13 and 14 : second derivatives with respect to spin density
                !    and gradient, divided by the gradient
                dvxci(ipts,13)=d2rhohh_drhodg+dtt_dg*d2hh_dttdzeta
                dvxci(ipts,14)=d2rhohh_drhodg-dtt_dg*d2hh_dttdzeta
                !    Component 15 : derivative of the (derivative wrt the gradient div by the grad),
                !    divided by the grad
                dvxci(ipts,15)=rhotot*d2hh_dtt2*dtt_dg*dtt_dg

                !  End condition of GGA
             end if

          else  ! no correlation

             if(ndvxci > 8)then
                !   Must zero the correlation part of the xc kernel
                dvxci(:,9:15)=zero
             end if

             ! End condition of including correlation, and not only exchange
          end if

          ! Correlation has been added
          ! -----------------------------------------------------------------------------

       end do
    end if
 else
    ! Disallowed value for nspden
    write(message, '(a,a,a,a,a,a,i12,a)' ) ch10,&
         &   ' xcpbe: BUG -',ch10,&
         &   '  Argument nspden must be 1 or 2; ',ch10,&
         &   '  Value provided as argument was ',nspden,'.'
    call wrtout(06,message,'COLL')
    call leave_new('COLL')
 end if

!DEBUG
!Finite-difference debugging, do not take away
!if(debug/=0)then
! do ipts=1,npts,5

!  rho=rho_updn(ipts,1)+rho_updn(ipts,2)
!  write(6, '(a,i5,a,es16.8)' ) ' Point number',ipts,' with rho=',rho

!! For rho
!  if(debug==1)then
!   write(6, '(3es16.8)' )exci(ipts)*rho,vxci(ipts,1),vxci(ipts,2)
!  else
!!  For grho2
!   write(6, '(4es16.8)' )exci(ipts)*rho,dvxcdgr(ipts,1),&
! &  dvxcdgr(ipts,2),dvxcdgr(ipts,3)
!  end if

!  write(6, '(4es16.8)' )dvxci(ipts,1:4)  ! For exchange
!  write(6, '(4es16.8)' )dvxci(ipts,5:8)  ! For exchange
!  write(6, '(4es16.8)' )dvxci(ipts,9:12) ! For correlation
!  write(6, '(4es16.8)' )dvxci(ipts,13:15)! For correlation

!  if(debug==1)then
!!  For rho
!   write(6, '(3es16.8)' )exci(ipts)*rho,&
!&      ( exci(ipts+1)*(rho+delta) - exci(ipts+2)*(rho-delta) )/2._dp/delta,&
!&      ( exci(ipts+3)*(rho+delta) - exci(ipts+4)*(rho-delta) )/2._dp/delta
!   write(6, '(3es16.8)' )&
!&    ( vxci(ipts+1,1) - vxci(ipts+2,1) )/2._dp/delta,&
!&    ( vxci(ipts+3,1) - vxci(ipts+4,1) )/2._dp/delta,&
!&    ( vxci(ipts+3,2) - vxci(ipts+4,2) )/2._dp/delta
!   write(6, '(4es16.8)' )&
!&    ( dvxcdgr(ipts+1,1) - dvxcdgr(ipts+2,1) )/2._dp/delta,&
!&    ( dvxcdgr(ipts+3,2) - dvxcdgr(ipts+4,2) )/2._dp/delta,&
!&    ( dvxcdgr(ipts+1,3) - dvxcdgr(ipts+2,3) )/2._dp/delta,&
!&    ( dvxcdgr(ipts+3,3) - dvxcdgr(ipts+4,3) )/2._dp/delta
!  else
!!  For grho2  (should distinguish exchange and correlation ...)
!   grr=sqrt(grho2_updn(ipts,1)) ! Analysis of exchange
!   grr=sqrt(grho2_updn(ipts,3)) ! Analysis of correlation
!   write(6, '(3es16.8)' )exci(ipts)*rho,&
!&      ( exci(ipts+1)*rho - exci(ipts+2)*rho )/2._dp/delta/grr,&
!&      ( exci(ipts+3)*rho - exci(ipts+4)*rho )/2._dp/delta/grr
!   write(6, '(3es16.8)' )&
! &    ( vxci(ipts+1,1) - vxci(ipts+2,1) )/2._dp/delta/grr,&
! &    ( vxci(ipts+3,1) - vxci(ipts+4,1) )/2._dp/delta/grr,&
! &    ( vxci(ipts+3,2) - vxci(ipts+4,2) )/2._dp/delta/grr
!   write(6, '(4es16.8)' )&
! &    ( dvxcdgr(ipts+1,1) - dvxcdgr(ipts+2,1) )/2._dp/delta/grr,&
! &    ( dvxcdgr(ipts+3,2) - dvxcdgr(ipts+4,2) )/2._dp/delta/grr,&
! &    ( dvxcdgr(ipts+1,3) - dvxcdgr(ipts+2,3) )/2._dp/delta/grr,&
! &    ( dvxcdgr(ipts+3,3) - dvxcdgr(ipts+4,3) )/2._dp/delta/grr
!  end if
! end do
! stop
!end if
!ENDDEBUG

 deallocate(rhoarr,rhom1_3,rho_updnm1_3)
 deallocate(zetm,zetmm1_3,zetp,zetpm1_3)

!DEBUG
!write(6,*)' xcpbe : exit'
!write(6,*)' nspden=',nspden
!if(order==2)stop
!ENDDEBUG

end subroutine xcpbe
!!***
