!=============================================================================
! amrvacusr.t.binaryAccretion  simulating binary wind accretion,  11/2014
!=============================================================================
!INCLUDE:amrvacnul/specialini.t
!INCLUDE:amrvacnul/speciallog.t
!INCLUDE:amrvacnul/specialbound.t
!INCLUDE:amrvacnul/specialsource.t
INCLUDE:amrvacmodules/cooling.t
INCLUDE:amrvacnul/specialimpl.t
!INCLUDE:amrvacnul/usrflags.t

! update : 20/11/2014
! made by Tom 
! configuration :
! 2D:
!$AMRVAC_DIR/setup.pl -d=23 -phi=0 -z=0 -g=10,10 -p=hd -eos=default -nf=1 -ndust=0 -u=nul
! 3D:
!$AMRVAC_DIR/setup.pl -d=33 -phi=0 -z=0 -g=10,10,10 -p=hd -eos=default -nf=1 -ndust=0 -u=nul
! typically use -arch=default on your machine or -arch=thinking on the hpc cluster
!=============================================================================
subroutine initglobaldata_usr

include 'amrvacdef.f'

integer :: fileid, i
character*79 :: filename

double precision:: qrhodust {#IFDEF DUST ,r(0:^NDS)}
double precision:: frat,rbalminu,rbalplus
!-----------------------------------------------------------------------------
eqpar(gamma_)=4.0d0/3.0d0   ! ratio of specific heats, monoatomic ideal gas
eqpar(Mue_)=one 

! total mass accreted by the binary accretor
acc_mass_mype = zero
acc_mass_global = zero
acc_mass_tot = zero


! ism parameter, density in cgs
eqpar(rhoism_) = 1.0d-16

! star parameters, all cgs units here
! A is the donor star
eqpar(vinfA_)  = 1.5d6       ! terminal velocity 
eqpar(rinfA_)  = 200.0*rsolar 
eqpar(TwindA_) = 3.0d3
eqpar(rhoinfA_)= &
      1.0d-6*msolar/(sinyear*4.0d0*dpi*eqpar(rinfA_)**2*eqpar(vinfA_))

! B is the accretor, these parameters are not relevant in this setup 
eqpar(vinfB_)  = 1.5d6       
eqpar(rinfB_)  = 5.0d0*rsolar 
eqpar(TwindB_) = 3.0d3
eqpar(rhoinfB_)= &
      1.0d-6*msolar/(sinyear*4.0d0*dpi*eqpar(rinfB_)**2*eqpar(vinfB_))


! binary parameters
eqpar(MassA_)=3.0d0*msolar   ! mass for star A
eqpar(MassB_)=1.5d0*msolar   ! mass for star B
eqpar(Porb_) =(895.3d0/365.0d0)*sinyear       ! orbital period for binary
eqpar(ecc_)  =0.00d0         ! eccentricity for orbit of binary    

eqpar(semimajor_)=(ggrav/(4.0D0*dpi*dpi)* &
                   (eqpar(MassA_)+eqpar(MassB_))*eqpar(Porb_)**2)**(1.0d0/3.0d0)

eqpar(tfix_)=0.0*eqpar(Porb_) ! time during which the binary is not moved

! normalization
normt         = eqpar(Porb_)
normvar(0)    = eqpar(semimajor_)
normvar(rho_) = 1.0D-16    
{^C&normvar(v^C_)  = normvar(0)/normt\}
normvar(e_)   = normvar(rho_)*normvar(v1_)*normvar(v1_)
normvar(tr1_)  = normvar(rho_)
! following is actually 1/T_0, with T_0 reference temperature in cgs
eqpar(Tscale_) = (1.0D0/(normvar(v1_)**2.0d0)) &
               * kboltz/mhydro
! luminosity scale
eqpar(Lscale_) = normvar(rho_)*normt/((mhydro*normvar(v1_))**2.0)
   

if(mype == 0) then
   fileid = 10
   write(filename,"(a,a)") TRIM(filenamelog),".scale"
   open(fileid,file=filename, status='unknown')
   write(fileid,*) 'normt:        ', normt
   write(fileid,*) 'normvar(0):   ', normvar(0)
   write(fileid,*) 'normvar(v1_): ', normvar(v1_)
   write(fileid,*) 'normvar(rho_):', normvar(rho_)    
   write(fileid,*) 'normvar(e_):  ', normvar(p_)    
   write(fileid,*) 
   write(fileid,*) 'accel        ', normvar(v1_)*normvar(v1_)/normvar(0)
   write(fileid,*)
endif

! switch to normalized parameters
eqpar(rhoism_)= eqpar(rhoism_)  / normvar(rho_)

eqpar(rhoinfA_)=eqpar(rhoinfA_) / normvar(rho_)
eqpar(vinfA_)  =eqpar(vinfA_)   / normvar(v1_)
eqpar(rinfA_)  =eqpar(rinfA_)   / normvar(0) 
                            
eqpar(rhoinfB_)=eqpar(rhoinfB_) / normvar(rho_)
eqpar(vinfB_)  =eqpar(vinfB_)   / normvar(v1_) 
eqpar(rinfB_)  =eqpar(rinfB_)   / normvar(0)
                            
eqpar(MassA_)= eqpar(MassA_)    / (normvar(rho_)*(normvar(0)**3))  
eqpar(MassB_)= eqpar(MassB_)    / (normvar(rho_)*(normvar(0)**3))
eqpar(Porb_) = eqpar(Porb_)     /  normt 
eqpar(tfix_) = eqpar(tfix_)     /  normt 
eqpar(semimajor_)= eqpar(semimajor_)/normvar(0)

eqpar(TwindA_) = eqpar(TwindA_)*eqpar(Tscale_)
eqpar(TwindB_) = eqpar(TwindB_)*eqpar(Tscale_)



! dimensionless ISM pressure: temperature of ISM set to 10000 K
eqpar(pism_)=eqpar(rhoism_)*1.0d4*eqpar(Tscale_)



if(mype == 0) then
   write(fileid,*)
   write(fileid,*) trim(eqparname)
   write(fileid,*) trim(specialparname)
   do i=1,neqpar+nspecialpar
      write(fileid,*)'eqpar(',i,'): ', eqpar(i)
   enddo
   close(fileid)
endif

call coolinit
! the following says that 10000 K, measured in our T units, 
! is the minimal temperature, no cooling below it
tcoolmin=1.0d4*eqpar(Tscale_)
! enforcing a floor temperature in process
tlow=1.0d4*eqpar(Tscale_)


call MPI_BARRIER(icomm,ierrmpi)

if(npe>1)then
    call MPI_BCAST(nxmapA,1,MPI_INTEGER,0,icomm,ierrmpi)
    call MPI_BCAST(rmapA,jmax,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
    call MPI_BCAST(logrhomapA,jmax,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
    call MPI_BCAST(vmapA,jmax,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
    call MPI_BCAST(nxmapB,1,MPI_INTEGER,0,icomm,ierrmpi)
    call MPI_BCAST(rmapB,jmax,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
    call MPI_BCAST(logrhomapB,jmax,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
    call MPI_BCAST(vmapB,jmax,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
endif

return
end subroutine initglobaldata_usr
!=============================================================================
subroutine initonegrid_usr(ixG^L,ix^L,w,x)

! initialize one grid within ix^L

include 'amrvacdef.f'

integer, intent(in) :: ixG^L, ix^L
double precision, intent(in) :: x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)

double precision :: tstart,xstarA^D,xstarB^D,mfaca,mfacb
double precision :: RtostarA(ixG^T),RtostarB(ixG^T)

logical :: patchw(ixG^T)

integer:: iorbstep, norbstep
logical, save :: first=.true.
!-----------------------------------------------------------------------------

if(mype==0.and.first)then
  print *,'============================================================='
  print *,' HD simulation for stellar binary wind interaction'
  print *,'============================================================='
  print *,'equation and ISM parameters are:'
  print *,' gamma set to:', eqpar(gamma_),' Pism=',eqpar(pism_) 
  print *,' density for ISM:',eqpar(rhoism_)
  print *,' hence dimensionless T for ISM is=',eqpar(pism_)/eqpar(rhoism_)
  print *,' while dimensionless T for star A and B=',eqpar(TwindA_),eqpar(TwindB_)
  print *,' and no cooling below dimensionless temperature =',tlow
  print *,'============================================================='
  print *,'stellar parameters for star A: mass, vinf, rinf, rhoinf:' 
  print *,eqpar(MassA_),eqpar(vinfA_),eqpar(rinfA_),eqpar(rhoinfA_)
  print *,'============================================================='
  print *,'stellar parameters for star B: mass, vinf, rinf, rhoinf:' 
  print *,eqpar(MassB_),eqpar(vinfB_),eqpar(rinfB_),eqpar(rhoinfB_)
  print *,'============================================================='
  print *,'further orbital parameters: Porb, ecc, semimajor:', eqpar(Porb_),eqpar(ecc_),eqpar(semimajor_)
  print *,'============================================================='
  print *,' periastron separation is:',eqpar(semimajor_)*(one-eqpar(ecc_))
  mfaca=eqpar(MassB_)/(eqpar(MassA_)+eqpar(MassB_))
  mfacb=eqpar(MassA_)/(eqpar(MassA_)+eqpar(MassB_))
  print *,' periastron A x-location should be:', &
       -mfaca*eqpar(semimajor_)*&
       (one-eqpar(ecc_))
  print *,' periastron B x-location should be:',&
       mfacb*eqpar(semimajor_)*&
       (one-eqpar(ecc_))
  if(.false.)then
    print *,'============================='
    print *,'Testing orbit calculation'
    print *,' half orbit A-location should be:', &
        mfaca*eqpar(semimajor_)*&
        (one+eqpar(ecc_))
    print *,' half orbit B-location should be:',&
        -mfacb*eqpar(semimajor_)*&
        (one+eqpar(ecc_))
    print *,'============================='
    norbstep=8
    do iorbstep=0,norbstep-1
       tstart=eqpar(Porb_)*dble(iorbstep)/dble(norbstep-1)
       call getbinarylocations(tstart,xstarA^D,xstarB^D)
       print *,'fract. orbit t=',tstart/eqpar(Porb_),' pos. star A ',xstarA^D
       print *,'fract. orbit t=',tstart/eqpar(Porb_),' pos. star B ',xstarB^D
    enddo
    print *,'============================='
  endif
  print *,'Note we run with temperature fix set to: Tfix=',Tfix,' Tfixprim=',Tfixprim
  first=.false.
endif



! initialize everywhere with static ISM
w(ix^S,rho_)=eqpar(rhoism_)
{^C& w(ix^S,v^C_)=zero \}
w(ix^S,p_)=eqpar(pism_)
w(ix^S,tr1_)=0.0d0

! switch to conservative variables
patchw(ixG^S)=.false.
call conserve(ixG^L,ix^L,w,x,patchw)



tstart=zero
call getRtostar(ixG^L,tstart,w,x,RtostarA,RtostarB,xstarA^D,xstarB^D)
call putwindzone(ixG^L,tstart,w,x,RtostarA,RtostarB,xstarA^D,xstarB^D)


end subroutine initonegrid_usr
!=============================================================================
double precision function myinterpolation(yy1,yy2,mu)
include 'amrvacdef.f'
double precision :: yy1,yy2,mu
!----------------------------------------------------------------------------

!!myinterpolation=yy1*(one-half*(one-dcos(mu*dpi)))+yy2*half*(one-dcos(mu*dpi))
myinterpolation=yy1*(one-mu)+yy2*mu

end function myinterpolation
!=============================================================================
subroutine getRtostar(ixI^L,qt,w,x,RtostarA,RtostarB,xstarA^D,xstarB^D)

include 'amrvacdef.f'

integer, intent(in) :: ixI^L
double precision, intent(in) :: qt, x(ixI^S,1:ndim)
double precision, intent(inout) :: w(ixI^S,1:nw)
double precision, intent(out) :: RtostarA(ixG^T), RtostarB(ixG^T)
double precision, intent(out) :: xstarA^D,xstarB^D

double precision :: xstarAnew^D,xstarBnew^D
logical :: ingridA, ingridB, relocate
integer :: inearest^D(1:ndim)
!-----------------------------------------------------------------------------

call getbinarylocations(qt,xstarA^D,xstarB^D)

! Fix star location in center of nearest gridpoint 
! if the star is inside this grid

relocate=.false.

if(relocate)then
call staringrid(ixI^L,x,xstarA^D,ingridA)
call staringrid(ixI^L,x,xstarB^D,ingridB)

if(ingridA)then
  {^D&inearest^D = minloc(DABS(x(ixI^S,^D)-xstarA^D)) \}
  xstarAnew1 = x(inearest1(1),inearest2(1){^IFTHREED,inearest3(1)},1)
  xstarAnew2 = x(inearest1(1),inearest2(2){^IFTHREED,inearest3(1)},2)
  {^IFTHREED
  xstarAnew3 = x(inearest1(1),inearest2(1),inearest3(3),3)
  if(dabs(xstarA3-xstarAnew3)>1.0d-6)then 
    print *,'relocating A in direction 3 not ok'
    print *,'qt=',qt
    print *,x(ixImin1,ixImin2,ixImin3:ixImax3,3)
    print *,'versus xstarA3=',xstarA3
    if(qt>zero)call mpistop("fix orbital plane!")
  endif
  }
  if(dabs(xstarA1-xstarAnew1)>0.01d0.or.dabs(xstarA2-xstarAnew2)>0.01d0)then 
    print *,'relocating A in direction 1 to ',xstarAnew1,' from ',xstarA1
    print *,'             in direction 2 to ',xstarAnew2,' from ',xstarA2
    if(qt>zero)call mpistop("fix orbit versus planar resolution!")
  endif 
  xstarA^D=xstarAnew^D;
endif

if(ingridB)then
  {^D&inearest^D = minloc(DABS(x(ixI^S,^D)-xstarB^D)) \}
  xstarBnew1 = x(inearest1(1),inearest2(1){^IFTHREED,inearest3(1)},1)
  xstarBnew2 = x(inearest1(1),inearest2(2){^IFTHREED,inearest3(1)},2)
  {^IFTHREED
  xstarBnew3 = x(inearest1(1),inearest2(1),inearest3(3),3)
  if(dabs(xstarB3-xstarBnew3)>1.0d-6)then 
    print *,'relocating B in direction 3 not ok'
    print *,'qt=',qt
    print *,x(ixImin1,ixImin2,ixImin3:ixImax3,3)
    print *,'versus xstarB3=',xstarB3
    if(qt>zero)call mpistop("fix orbital plane!")
  endif
  }
  if(dabs(xstarB1-xstarBnew1)>0.01d0.or.dabs(xstarB2-xstarBnew2)>0.01d0)then 
    print *,'relocating B in direction 1 to ',xstarBnew1,' from ',xstarB1
    print *,'             in direction 2 to ',xstarBnew2,' from ',xstarB2
    if(qt>zero)call mpistop("fix orbit versus planar resolution!")
  endif 
  xstarB^D=xstarBnew^D;
endif
endif


! avoid divisions by zero
RtostarA(ixI^S)=max( dsqrt( {^D&(x(ixI^S,^D)-xstarA^D)**2+} ), 0.01d0*eqpar(rinfA_))
RtostarB(ixI^S)=max( dsqrt( {^D&(x(ixI^S,^D)-xstarB^D)**2+} ), 0.01d0*eqpar(rinfB_))


end subroutine getRtostar
!=============================================================================
subroutine putwindzone(ixI^L,qt,w,x,RtostarA,RtostarB,xstarA^D,xstarB^D)

include 'amrvacdef.f'

integer, intent(in) :: ixI^L
double precision, intent(in) :: qt, x(ixI^S,1:ndim)
double precision, intent(inout) :: w(ixI^S,1:nw)
double precision, intent(in) :: RtostarA(ixG^T), RtostarB(ixG^T)
double precision, intent(in) :: xstarA^D,xstarB^D

integer :: iw
double precision :: xA^CinvR(ixG^T), xB^CinvR(ixG^T)
double precision :: outrangeA,outrangeB
!-----------------------------------------------------------------------------

{^IFTWOD
xA3invR(ixI^S)=zero
xB3invR(ixI^S)=zero
}
{^D&
xA^DinvR(ixI^S)=(x(ixI^S,^D)-xstarA^D)/RtostarA(ixI^S)
xB^DinvR(ixI^S)=(x(ixI^S,^D)-xstarB^D)/RtostarB(ixI^S)
\}

if(qt==0.0) then
   outrangeA = 10.0d0
   outrangeB = 1.0d0
else
   outrangeA = 1.0d0
   outrangeB = 1.0d0
end if


do iw= 1,nw
    select case(iw)
    case(tr1_)
       where(RtostarA(ixI^S)<=outrangeA*eqpar(rinfA_))
            w(ixI^S,tr1_)= 1000.0d0
       end where
       where(RtostarB(ixI^S)<=outrangeB*eqpar(rinfB_))
            w(ixI^S,tr1_)=-1000.0d0
       end where      
    case(rho_)
       ! ensure a density variation consistent with a fixed mass loss rate as function of radius
       ! throughout the sphere of radius rinf
       where(RtostarA(ixI^S)<=outrangeA*eqpar(rinfA_))
            w(ixI^S,rho_)= (eqpar(rinfA_)/RtostarA(ixI^S))**2*eqpar(rhoinfA_)
       end where
       where(RtostarB(ixI^S)<=outrangeB*eqpar(rinfB_))
            !w(ixI^S,rho_)= (eqpar(rinfB_)/RtostarB(ixI^S))**2*eqpar(rhoinfB_)
       end where      
    {case(m^C_)
       where(RtostarA(ixI^S)<=outrangeA*eqpar(rinfA_))
            w(ixI^S,m^C_)= (eqpar(rinfA_)/RtostarA(ixI^S))**2*eqpar(rhoinfA_)*eqpar(vinfA_)*xA^CinvR(ixI^S)
       end where
       where(RtostarB(ixI^S)<=outrangeB*eqpar(rinfB_))
            !w(ixI^S,m^C_)= (eqpar(rinfB_)/RtostarB(ixI^S))**2*eqpar(rhoinfB_)*eqpar(vinfB_)*xB^CinvR(ixI^S)
       end where  \}
    case(e_)
       where(RtostarA(ixI^S)<=outrangeA*eqpar(rinfA_))
          w(ixI^S,e_)= (eqpar(rinfA_)/RtostarA(ixI^S))**2*eqpar(rhoinfA_) &
                      *(half*eqpar(vinfA_)**2 + eqpar(TwindA_)/(eqpar(gamma_)-one))
       end where
       where(RtostarB(ixI^S)<=outrangeB*eqpar(rinfB_))
          !w(ixI^S,e_)= (eqpar(rinfB_)/RtostarB(ixI^S))**2*eqpar(rhoinfB_) &
          !            *(half*eqpar(vinfB_)**2 + eqpar(TwindB_)/(eqpar(gamma_)-one))
       end where
   end select
end do

end subroutine putwindzone
!=============================================================================
subroutine specialsource(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

! Calculate w(iw)=w(iw)+qdt*SOURCE[wCT,qtC,x] within ixO for all indices
! iw=iwmin...iwmax.  wCT is at time qCT

include 'amrvacdef.f'

integer, intent(in) :: ixI^L, ixO^L, iw^LIM
double precision, intent(in) :: qdt, qtC, qt, x(ixI^S,1:ndim)
double precision, intent(inout) :: wCT(ixI^S,1:nw), w(ixI^S,1:nw)

double precision:: R2tostarA(ixG^T),R2tostarB(ixG^T),xstarA^D,xstarB^D
double precision:: accMass(ixG^T),newvdust^DS(ixG^T,1:ndir)
logical :: patchstars(ixG^T),patchshock(ixG^T)
double precision :: pmassA,pmassB,rsmooth,ggravN,R_bondi
!-----------------------------------------------------------------------------

ggravN = ggrav / (1.0d0 / (normvar(rho_)*normt*normt))
call getbinarylocations(qt,xstarA^D,xstarB^D)
pmassA = eqpar(MassB_)*ggravN
pmassB = eqpar(MassB_)*ggravN
rsmooth = eqpar(rinfB_)
R2tostarA(ixO^S) = {^D&(x(ixO^S,^D)-xstarA^D)**2+}
R2tostarB(ixO^S) = {^D&(x(ixO^S,^D)-xstarB^D)**2+} 

if(qt<eqpar(Porb_)) then
    pmassA = pmassA*dsqrt(qt/eqpar(Porb_))
    pmassB = pmassB*dsqrt(qt/eqpar(Porb_))
endif

! add sources from gravity
!where(R2tostarA(ixO^S)>eqpar(rinfA_)**2)
    w(ixI^S,tr1_) = zero
    ! de/dt= +g_i*m_i
    {^D&w(ixO^S,ee_)=w(ixO^S,ee_) &
            -qdt*pmassA*wCT(ixO^S,m0_+^D)*(x(ixO^S,^D)-xstarA^D)&
            /(rsmooth**2 + R2tostarA(ixO^S))**(1.5d0);}
    ! dm_i/dt= +rho*g_i            
    {^D&w(ixO^S,m0_+^D)=w(ixO^S,m0_+^D)-qdt*pmassA*wCT(ixO^S,rho_)* &
                (x(ixO^S,^D)-xstarA^D)&
                /(rsmooth**2 + R2tostarA(ixO^S))**(1.5d0);}        
!end where
!where(R2tostarA(ixO^S)>eqpar(rinfB_)**2)
    w(ixI^S,tr1_) = zero       
    ! de/dt= +g_i*m_i
    {^D&w(ixO^S,ee_)=w(ixO^S,ee_) &
            -qdt*pmassB*wCT(ixO^S,m0_+^D)*(x(ixO^S,^D)-xstarB^D)&
            /(rsmooth**2 + R2tostarB(ixO^S))**(1.5d0);}
    ! dm_i/dt= +rho*g_i
    {^D&w(ixO^S,m0_+^D)=w(ixO^S,m0_+^D)-qdt*pmassB*wCT(ixO^S,rho_)* &
                (x(ixO^S,^D)-xstarB^D)&
                /(rsmooth**2 + R2tostarB(ixO^S))**(1.5d0);}
!end where


! Bondi accretion
if (.true.) then
    R_bondi = 2.0d0*pmassB/(cmax_global**2 + 0.1)
    R_bondi = max(min(R_bondi,10.0*eqpar(rinfB_)),eqpar(rinfB_))
    where (R2tostarB(ixO^S)<R_bondi**2) 
        accMass(ixO^S) = qdt*3.0d0*wCT(ixO^S,rho_) &
		    *pmassB**2/(R_bondi*cmax_global)**3
		w(ixO^S,rho_) = max(w(ixO^S,rho_),3.0*minrho)
		where((w(ixO^S,rho_)-accMass(ixO^S))>3.0*minrho)
	        w(ixO^S,rho_) = w(ixO^S,rho_) - accMass(ixO^S)
	    elsewhere
	        accMass(ixO^S) = w(ixO^S,rho_) - 3.0*minrho
	        w(ixO^S,rho_) = 3.0*minrho
	    end where
	    accMass(ixO^S) = {^D&dxlevel(^D)*}*accMass(ixO^S)	        
    elsewhere
        accMass(ixO^S) = zero
    end where
    acc_mass_mype=acc_mass_mype+SUM(accMass(ixO^S))
end if

call addsource_cooling(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

return
end subroutine specialsource
!=============================================================================
subroutine getbinarylocations(qt,xstarA^D,xstarB^D)

! this subroutine determines the position of the binary component stars
! at time qt, we compute the x and y locations of both stars

! we assume initial position such that A and B are at periastron on time zero
! both along the x-axis (time of periastron passage and longitude of periastron set to zero)

include 'amrvacdef.f'

double precision, intent(in) :: qt
double precision, intent(out) :: xstarA^D, xstarB^D

integer :: maxiterb, iter
double precision:: epsilon, ksiold, ksinew, dAB, phiAB, dABx^D, time
!-----------------------------------------------------------------------------

! accuracy asked for fixed point iteration on Kepler's equation
epsilon=1.0d-8
! maximum number of iterations performed
maxiterb=1000

! ensure time between 0 en Porb, fix locations up to tfix
if (qt<=eqpar(tfix_)) then
  time=zero
else
  time=qt-eqpar(tfix_)
  do while((time-eqpar(Porb_))>epsilon)
     time=time-eqpar(Porb_)
  enddo
endif

! fixed point iteration on Kepler's equation
iter=0
! start guess is zero value
ksiold=0.0d0
do while(iter<=maxiterb)
     ksinew=eqpar(ecc_)*dsin(ksiold)+2.0d0*dpi*time/eqpar(Porb_)
     if(dabs(ksinew-ksiold)<epsilon) exit
     ksiold=ksinew
     iter=iter+1
enddo

! stop when fixed point iteration failed
if(iter==maxiterb)then
       print *,'qt,iter,ksiold,ksinew,epsilon'
       print *,qt,iter,ksiold,ksinew,epsilon
     call mpistop("iteration for Kepler failed within maxiterb iterations")
endif

if(time>eqpar(Porb_)/2.0d0) then
	!ksinew=ksinew+dpi
endif

! distance from A to B
dAB=eqpar(semimajor_)*(one-eqpar(ecc_)*dcos(ksinew))
! polar angle from A to B, as viewed from star A
phiAB=dacos((dcos(ksinew)-eqpar(ecc_))/(one-eqpar(ecc_)*dcos(ksinew)))


!if(time>eqpar(Porb_)/2.0d0)dAB=-dAB
if(time>eqpar(Porb_)/2.0d0)phiAB=2.0*dpi-phiAB
if(time>1.0) write(*,*) time,phiAB,qt

! x and y components of vector from star A to star B
dABx1=dAB*dcos(phiAB)
dABx2=dAB*dsin(phiAB)

{^IFTHREED dABx3=zero }

! positions of stars A and B
xstarA^D=-eqpar(MassB_)*dABx^D/(eqpar(MassA_)+eqpar(MassB_));
xstarB^D=+eqpar(MassA_)*dABx^D/(eqpar(MassA_)+eqpar(MassB_));

end subroutine getbinarylocations
!============================================================================
subroutine getdt_special(w,ixG^L,ix^L,dtnew,dx^D,x)

! Limit "dt" further if necessary, e.g. due to the special source terms.
! The getdt_courant (CFL condition) and the getdt subroutine in the AMRVACPHYS
! module have already been called.

include 'amrvacdef.f'

integer, intent(in) :: ixG^L, ix^L
double precision, intent(in) :: dx^D, x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw), dtnew

double precision :: dtcool
double precision:: dxinv(1:ndim), dtgrav,dtacc
double precision:: R2tostarA(ixG^T),R2tostarB(ixG^T),xstarA^D,xstarB^D
double precision :: pmassA,pmassB,rsmooth,ggravN,R_bondi,maxDeltaRho
!-----------------------------------------------------------------------------

dtnew=bigdouble
! use a CFL like condition based on (spherical) orbital revolution rate
dtnew=min(dtnew,dtdiffpar*dx1/(2.0d0*dpi*eqpar(semimajor_)/eqpar(Porb_)))
dtnew=min(dtnew,dtdiffpar*dx2/(2.0d0*dpi*eqpar(semimajor_)/eqpar(Porb_)))

call getdt_cooling(w,ixG^L,ix^L,dtcool,dx^D,x)
dtnew = min(dtnew,dtcool)

if(dtnew<dtmin)then
   write(unitterm,*)"-------------------------------------"
   write(unitterm,*)"Warning: found getdt_special related time step too small! dtnew=",dtnew
   write(unitterm,*)"on grid with index:", saveigrid," grid level=",node(plevel_,saveigrid)
   write(unitterm,*)"grid corners are=",{^D&rnode(rpxmin^D_,saveigrid), rnode(rpxmax^D_,saveigrid)}
   write(unitterm,*)"orbital dt1=",dtdiffpar*dx1/(2.0d0*dpi*eqpar(semimajor_)/eqpar(Porb_))
   write(unitterm,*)"orbital dt2=",dtdiffpar*dx2/(2.0d0*dpi*eqpar(semimajor_)/eqpar(Porb_))
   write(unitterm,*)"dtcool=",dtcool
   write(unitterm,*)"on processor:", mype
   write(unitterm,*)"-------------------------------------"
endif


! dt condition imposed by gravity
ggravN = ggrav / (1.0d0 / (normvar(rho_)*normt*normt))
pmassB = eqpar(MassB_)*ggravN
rsmooth = eqpar(rinfB_)
call getbinarylocations(t,xstarA^D,xstarB^D)
R2tostarB(ixG^S) = {^D&(x(ixG^S,^D)-xstarB^D)**2+} 
if(t<eqpar(Porb_)) then
    pmassB = pmassB*dsqrt(t/eqpar(Porb_))
endif
dtgrav=bigdouble
{^D&dtgrav=min(dtgrav,dsqrt(minval(dx^D*(rsmooth**2 + R2tostarB(ixG^S))**(1.5d0) &
                /(pmassB*dabs(x(ixG^S,^D)-xstarB^D)))));}
dtnew = min(dtnew,dtdiffpar*dtgrav)



! Bondi accretion
if (.true.) then
    R_bondi = 2.0d0*pmassB/(cmax_global**2 + 0.1)
    R_bondi = max(min(R_bondi,10.0*eqpar(rinfB_)),eqpar(rinfB_))
    dtacc=bigdouble
    dtacc=0.1*(R_bondi*cmax_global)**3 / &
            (3.0d0*pmassB**2)
    dtnew = min(dtnew,dtacc)
end if


end subroutine getdt_special
!=============================================================================
subroutine specialeta(w,ixI^L,ix^L,idirmin,x,current,eta)

! Set the "eta" array for resistive MHD based on w or the
! "current" variable which has components between idirmin and 3.

include 'amrvacdef.f'

integer, intent(in) :: ixI^L, ix^L, idirmin
double precision, intent(in) :: w(ixI^S,nw), x(ixI^S,1:ndim)

double precision :: current(ixG^T,7-2*ndir:3), eta(ixG^T)
!-----------------------------------------------------------------------------

call mpistop("specialeta is not defined")

end subroutine specialeta
!=============================================================================
subroutine specialrefine_grid(igrid,level,ixG^L,ix^L,qt,w,x,refine,coarsen)

! Enforce additional refinement or coarsening
! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.

! you must set consistent values for integers refine/coarsen:

! refine = -1 enforce to not refine
! refine =  0 doesn't enforce anything
! refine =  1 enforce refinement

! coarsen = -1 enforce to not coarsen
! coarsen =  0 doesn't enforce anything
! coarsen =  1 enforce coarsen

include 'amrvacdef.f'

integer, intent(in) :: igrid, level, ixG^L, ix^L
double precision, intent(in) :: qt, x(ixG^S,1:ndim)
integer, intent(inout) :: refine, coarsen
double precision, intent(inout) :: w(ixG^S,1:nw)
double precision :: wtemp(ixG^S,1:nw)

double precision :: xstarA^D,xstarB^D
double precision :: RtostarA(ixG^T), RtostarB(ixG^T)
integer :: donorlevel,outlevels,plus
logical :: patchstars(ixG^T),patchshock(ixG^T)
!-----------------------------------------------------------------------------

! AMR level of the donor star
donorlevel = 4


call getRtostar(ixG^L,qt,w,x,RtostarA,RtostarB,xstarA^D,xstarB^D)

! enforce coarsest level at large distances from both stars
if (all(RtostarA(ix^S)>0.8d0*eqpar(semimajor_)).and.all(RtostarB(ix^S)>0.8d0*eqpar(semimajor_))) then
  refine=-1
  coarsen=1
else
  ! enforce maximal refinement within fraction of radius about each star
  if (any(RtostarB(ix^S)<=(0.05d0*eqpar(semimajor_)))) then
     refine=1
     coarsen=-1
  else
     if( level > donorlevel ) then
       refine = -1
       coarsen = 1
     else if( level == donorlevel ) then
       refine = -1
       coarsen = 0
     else if( level < donorlevel ) then
       refine=1
       coarsen=-1
     endif
  endif
endif


end subroutine specialrefine_grid
!=============================================================================
subroutine specialvarforerrest(ixI^L,ixO^L,iflag,w,var)

! this is the place to compute a local auxiliary variable to be used
! as refinement criterion for the Lohner error estimator only
!  -->it is then requiring an iflag>nw
! note that ixO=ixI=ixG, hence the term local (gradients need special attention!)

include 'amrvacdef.f'

integer, intent(in)          :: ixI^L,ixO^L,iflag
double precision, intent(in) :: w(ixI^S,1:nw)
double precision, intent(out):: var(ixG^T)
!-----------------------------------------------------------------------------

! we will refine based on velocity magnitude V only
select case(iflag)
  case(nw+1)
    var(ixI^S) = dsqrt( ^C&w(ixI^S,m^C_)**2+ ) / w(ixI^S,rho_)
  {^DS&case(nw+1+^DS)
    var(ixI^S) = dsqrt( ^C&w(ixI^S,m^Cd^DS_)**2+ ) \}
  case default
    call mpistop("no such specialvarforerrest")
end select

end subroutine specialvarforerrest
!=============================================================================
subroutine specialset_B0(ixI^L,ixO^L,x,wB0)

! Here one can add a steady (time-independent) potential background field

include 'amrvacdef.f'

integer, intent(in)           :: ixI^L,ixO^L
double precision, intent(in)  :: x(ixG^T,1:ndim)
double precision, intent(inout) :: wB0(ixI^S,1:ndir)
!-----------------------------------------------------------------------------
call mpistop(' abs(Busr)> 0, make change in parfile or in user file')

wB0(ixO^S,1:ndir)=wB0(ixO^S,1:ndir)
!!wB0(ixO^S,1:ndir)=wB0(ixO^S,1:ndir)+user defined steady potential field

end subroutine specialset_B0
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

double precision:: RtostarA(ixG^T),RtostarB(ixG^T),xstarA^D,xstarB^D
double precision:: tmprho(ixG^T),newvdust^DS(ixG^T,1:ndir)
logical :: patchstars(ixG^T),patchshock(ixG^T)
logical, save :: first=.true.
!-----------------------------------------------------------------------------

print *,'Problem!'


end subroutine process_grid_usr
!=================================================================
subroutine myfloortemperature(ixI^L,ixO^L,qt,w,x)
!
!  Force minimum temperature to a fixed temperature
!
include 'amrvacdef.f'

integer, intent(in)             :: ixI^L, ixO^L
double precision, intent(in)    :: qt, x(ixI^S,1:ndim)
double precision, intent(inout) :: w(ixI^S,1:nw)

double precision :: etherm(ixG^T), emin, tfloor
integer :: ix^D
!-----------------------------------------------------------------------------
tfloor = tlow

if(tfloor>zero)then

  call getpthermal(w,ixI^L,ixO^L,etherm)

  {do ix^DB = ixO^LIM^DB\}
    emin         = w(ix^D,rho_)*tfloor/(eqpar(gamma_)-one)
    etherm(ix^D) = etherm(ix^D)/(eqpar(gamma_)-one)
    if( etherm(ix^D) < emin ) w(ix^D,e_)     = w(ix^D,e_) - etherm(ix^D) + emin
  {enddo^D&\}

endif

return
end subroutine myfloortemperature
!=============================================================================
!!subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv,level)
subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)

! this subroutine can be used in convert, to add auxiliary variables to the
! converted output file, for further analysis using tecplot, paraview, ....
! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
!
! the array normconv can be filled in the (nw+1:nw+nwauxio) range with 
! corresponding normalization values (default value 1)

include 'amrvacdef.f'

integer, intent(in)                :: ixI^L,ixO^L
!!integer, intent(in)                :: level
double precision, intent(in)       :: x(ixI^S,1:ndim)
double precision                   :: w(ixI^S,nw+nwauxio)
double precision                   :: normconv(0:nw+nwauxio)

double precision:: p(ixG^T),wloc(ixI^S,nw)
double precision:: RtostarA(ixG^T),RtostarB(ixG^T),xstarA^D,xstarB^D
double precision:: tmprho(ixG^T)
logical :: patchstars(ixG^T),patchshock(ixG^T)
double precision :: rcool(ixG^T)
double precision :: R_bondi,pmassB,ggravN
!-----------------------------------------------------------------------------
call getRtostar(ixI^L,t,wloc,x,RtostarA,RtostarB,xstarA^D,xstarB^D)
w(ixO^S,tr1_)=zero
where(RtostarA(ixO^S)<eqpar(rinfA_))
    w(ixO^S,tr1_)= 1000.0d0
end where
where(RtostarB(ixO^S)<eqpar(rinfB_))
    w(ixO^S,tr1_)= -1000.0d0
end where

ggravN = ggrav / (1.0d0 / (normvar(rho_)*normt*normt))
pmassB = eqpar(MassB_)*ggravN
if(t<eqpar(Porb_)) then
    pmassB = pmassB*dsqrt(t/eqpar(Porb_))
endif

R_bondi = 2.0d0*pmassB/(cmax_global**2 + 0.1)
R_bondi = max(min(R_bondi,10.0*eqpar(rinfB_)),eqpar(rinfB_))
where (RtostarB(ixO^S)<R_bondi)
	w(ixO^S,tr1_)=w(ixO^S,tr1_)-100.0d0
end where

if(nwauxio>0)then
  if(saveprim)then

    if(Tfixprim)then
      w(ixO^S,nw+1)=dlog10(w(ixO^S,p_)/eqpar(Tscale_))
    else
      w(ixO^S,nw+1)=dlog10(w(ixO^S,p_)/w(ixO^S,rho_)/eqpar(Tscale_))
    endif
    w(ixO^S,nw+2)= dsqrt( ^C&w(ixO^S,v^C_)**2+ )*normvar(v1_) 
  
    call getvar_cooling(ixI^L,ixO^L,w,x,rcool,normconv)
    w(ixO^S,nw+3)=rcool(ixO^S)

  else
   wloc(ixI^S,1:nw)=w(ixI^S,1:nw)
   call getpthermal(wloc,ixI^L,ixO^L,p)
   w(ixO^S,nw+1)=dlog10(p(ixO^S)/w(ixO^S,rho_)/eqpar(Tscale_))
  endif
endif

end subroutine specialvar_output
!=============================================================================
subroutine specialvarnames_output

! newly added variables need to be concatenated with the varnames/primnames string

include 'amrvacdef.f'

logical, save :: firstloc=.true.
!-----------------------------------------------------------------------------

if(nwauxio>0.and.firstloc) then
! note to self: commented out 2D and 3D stuff

 primnames= TRIM(primnames)//' '//'logT vmag'
 wnames=TRIM(wnames)//' '//'logT vmag'


 primnames= TRIM(primnames)//' '//'cool'
 wnames=TRIM(wnames)//' '//'cool'

 firstloc=.false.
endif

end subroutine specialvarnames_output
!=============================================================================
subroutine specialbound_usr(qt,ixG^L,ixO^L,iw,iB,w,x)

! special boundary types, user defined

include 'amrvacdef.f'

integer, intent(in) :: ixG^L, ixO^L, iw, iB
double precision, intent(in) :: qt, x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)
!----------------------------------------------------------------------------

call mpistop("specialbound not defined")
!!!w(ixO^S,iw)=zero

end subroutine specialbound_usr
!=============================================================================
subroutine bc_int(level,qt,ixG^L,ixO^L,w,x)

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

include 'amrvacdef.f'

integer, intent(in) :: ixG^L,ixO^L,level
double precision, intent(in) :: qt
double precision, intent(inout) :: w(ixG^S,1:nw)
double precision, intent(in) :: x(ixG^S,1:ndim)

double precision:: RtostarA(ixG^T),RtostarB(ixG^T),xstarA^D,xstarB^D
!-----------------------------------------------------------------------------

call getRtostar(ixG^L,qt,w,x,RtostarA,RtostarB,xstarA^D,xstarB^D)

if(any(RtostarA(ixG^S)<=1.5*eqpar(rinfA_)).or.any(RtostarB(ixG^S)<=1.5*eqpar(rinfB_))) then
   call putwindzone(ixG^L,qt,w,x,RtostarA,RtostarB,xstarA^D,xstarB^D)
endif

end subroutine bc_int
!============================================================================
subroutine staringrid(ixI^L,x,xstar^D,ingrid)
!
! Returns a true logical if the star is inside this grid
!
include 'amrvacdef.f'

integer, intent(in) :: ixI^L
double precision, intent(in) :: x(ixI^S,1:ndim)
double precision, intent(in) :: xstar^D

logical, intent(out) :: ingrid
!-----------------------------------------------------------------------------
ingrid = .false.

if( xstar1>=minval(x(ixI^S,1)) .and. xstar1<=maxval(x(ixI^S,1)) )then
  if( xstar2>=minval(x(ixI^S,2)) .and. xstar2<=maxval(x(ixI^S,2)) ) then
   {^IFTHREED
    if( xstar3>=minval(x(ixI^S,3)) .and. xstar3<=maxval(x(ixI^S,3)) ) then
   }
    ingrid = .true.
   {^IFTHREED
    endif
   }
  endif
endif

end subroutine staringrid
!=======================================================================================
subroutine flag_grid_usr(qt,ixG^L,ixO^L,w,x,flag)

include 'amrvacdef.f'

integer, intent(in)             :: ixG^L, ixO^L
integer, intent(inout)          :: flag
double precision, intent(in)    :: qt
double precision, intent(inout) :: w(ixG^S,1:nw)
double precision, intent(in)    :: x(ixG^S,1:ndim)

! flag=-1 : Treat all cells active, omit deactivation (onentry, default)
! flag=0  : Treat as normal domain
! flag=1  : Treat as passive, but reduce by safety belt
! flag=2  : Always treat as passive

double precision:: RtostarA(ixG^T),RtostarB(ixG^T),xstarA^D,xstarB^D
!-----------------------------------------------------------------------------

flag=-1

!flag = 0
!call getRtostar(ixG^L,qt,w,x,RtostarA,RtostarB,xstarA^D,xstarB^D)
!if(all(RtostarA(ixG^S)<=0.5*eqpar(rinfA_)).or.all(RtostarB(ixG^S)<=0.5*&
!  eqpar(rinfB_))) then
!   flag = 2
!end if


end subroutine flag_grid_usr
!=============================================================================
! amrvacusr.t.binaryhdWR98a
!=============================================================================
