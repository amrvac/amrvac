module mod_srmhd_con2prim


 use mod_global_parameters
 use mod_srmhd_parameters
 use mod_srmhd_eos
 implicit none

contains
subroutine srmhd_con2prim(d,s,tau,b,lfac,xi,ierror)
  ! made by Z. MELIANI 14/02/2018
! This subroutine srmhd_con2prim is re-writed by Z. Meliani

! (D,S,tau,B) --> compute auxiliaries lfac and xi
implicit none

double precision:: lfac,xi
double precision:: d,s(1:ndir),tau,b(1:ndir)
integer         :: ierror
      
double precision:: sdotb,ssqr,bsqr,sb2
double precision:: xi1,xi2,xii,xih,xil,dxi,f,df
double precision:: temp,vsqr,fl,fh,lfacl,lfach,lfaci
double precision:: er,er1,xiprev,dplus,lfacmax
double precision:: fv,dfv,v(1:ndir),pressure
double precision:: lastxiok,lastlfacok,lastf
double precision:: dlfac,p,dpdxi
integer:: ni,nit,niiter

logical:: finished,tests,testsl,testsh
!-----------------------------------------------------------------------------


ierror=0

! ierror=1 : error on entry: must have D>=0, tau>=smallp/(gamma-1)
! ierror=2 : srmhd_maxiterationNR reached without convergence
! ierror=3 : no range for solution was find
! ierror=4 : final xi <small_xi
! ierror=5 : lower bound non-negative f(xi) value
! ierror=6 : v^2>1
! ierror=7 :  solution out from the range?
! ierror=8 : nonmonotonic function f?

if(d<small_density .or. tau<small_e)then
  ierror=1
  return
endif

bsqr=sum(b**2.0d0)
ssqr=sum(s**2.0d0)

! handle the case of no flow first
if(ssqr<=small_vec2)then
  call srmhd_get_h_noflux(d,tau+d-0.5*bsqr,xi) 
  lfac=1.0d0
  return
endif

! Hydro case: handle exactly as before (using p to iterate on)
!{#IFDEF ENERGY
! opedit: I don't know how this would be done if we don't have energy so removing for now.
if(bsqr<=small_vec2)then
  call srmhd_srhdcon2prim(d,s,tau,lfac,xi,ierror)
  return 
end if ! Hydro
!}

temp=sum(s*b)
if(temp<zero)then
  sdotb=-min(dsqrt(ssqr*bsqr),-temp)
else
  sdotb=min(dsqrt(ssqr*bsqr),temp)
endif
sb2=sdotb**2.0

! the maximal allowed Lorentz factor
lfacmax=1.0d0/dsqrt(1.0d0-srmhd_maxspeed**2.0d0)
! starting point for NR on xi
! xi1 is meant to be the lower bound for the bracket to be used in NR
dplus=d+srmhd_gamma*small_pressure/(srmhd_gamma-1.0d0)
xi1=dplus

! compute v^2(xi1)
vsqr = (ssqr + (2.0*xi1+ bsqr)*sb2/(xi1**2.0))&
      /((xi1+bsqr)**2.0)

!=== find new xi1, in the case v2(xi1) > maxvel^2 = (1-srmhd_maxdspeed)^2 
! locate xi corresponding to maximal velocity allowed, namely 1-srmhd_maxdspeed
niiter=-1
if(vsqr > srmhd_maxspeed**2.0d0 )then 
   er1=1.0d0
   xiprev=xi1
LoopVmax:  do ni = 1,srmhd_maxiterationNR
      if(ni>srmhd_maxiterationNR/2)then
        xi1=half*(xi1+xiprev)
        er1=10.0d0*er1
      endif

      ! v^2(xi1) - maxvel^2
      fv= ((ssqr + (2.0*xi1 + bsqr)*sb2/(xi1**2.0))&
        /((xi1+bsqr)**2.0))-srmhd_maxspeed**2.0d0
      ! d(v^2(xi)-maxvel^2)/dxi
      dfv= -two * (sb2*(3.0d0*xi1*(xi1+bsqr)+bsqr*bsqr)+ssqr*xi1**3)/ &
                     ((xi1*(xi1+bsqr))**3.0d0)
      if(fv<0)exit LoopVmax
      if(fv*dfv==zero) then
         if(fv==zero)then
            exit LoopVmax
         else
            !print *,'stop: dfv becomes zero, non-monotonic function of xi!!'
            ierror=8
            return
         endif
      else
        xiprev=xi1
        xi1   =xi1 -fv/dfv
        if(fv*dfv>zero)then
          ! xi-iterate decreased
          ! restrict to left
          xi1=max(xi1,dplus)
        else ! fv*dfv <0 
          ! xi-iterate increased
          ! restrict to right
          xi1=min(xi1,tau+d)
        endif
      endif
      er=dabs(fv/dfv)/xi1
      if((er<srmhd_tolerNR*er1).or.(dabs(fv/dfv)<srmhd_absaccnr))exit LoopVmax
      niiter=ni
   enddo LoopVmax
endif 
 !=============================================!

if(niiter==srmhd_maxiterationNR)then
 ! could not find consistent value of lower bound for xi compliant with maxvel
 print *,' could not find value of lower bound for xi compliant with maxvel'
 print *,'xi1=',xi1,'dplus=',dplus,' tau+d=',tau+d
 ierror=2
 return
endif

vsqr = (ssqr + (2.0*xi1+ bsqr)*sb2/(xi1**2.0))&
      /((xi1+bsqr)**2.0)
if(vsqr>=1.0)then
 if(fv*dfv>0.0)then
  xi1=xi1/2
 else
  xi1=2.0*xi1
 end if
end if
! we now compute f(xi1) and lfac(xi1)

call srmhd_funcd(xi1,fl,df,lfacl,d,ssqr,tau,bsqr,sdotb,ierror)
if(ierror /=0)return
testsl=(xi1>=dplus*lfacl.and.lfacl<=lfacmax)
 
if(fl>zero)then
  print *,'warning: lower bound non-negative f(xi) value!'
  print*,'iteration  = ',it
  print *,'xi1=',xi1,' vs dplus=',dplus,'dplus*lfacl=',dplus*lfacl&
         ,'lfacl=',lfacl,'lfacmax=',lfacmax
  print *,'fl=',fl,' niiter=',niiter,' vsqr(dplus)=',vsqr
  print*,' ierror = 5',it,mype,saveigrid,global_time
  ierror=5
  return
endif
 
!--------------------------------------------------------!
! xi2 is meant to be the maximal bound for the bracket to be used in NR
! for this we take the value of E=tau+d, increased by smallp

xi2= max(tau+d+small_pressure - half *bsqr,10.0d0*xi1)
niiter=-1

LoopxiMax : do ni=1,srmhd_maxiterationNR
    ! we now compute f(xi2) and lfac(xi2)
    call srmhd_funcd(xi2,fh,df,lfach,d,ssqr,tau,bsqr,sdotb,ierror)
    if(ierror /=0)return
    testsh=(xi2>=dplus*lfach.and.lfach<=lfacmax)

    ! maximal bound found when f(xi1) opposite sign of f(xi2)
    !, enforce consistency tests on xi2
    if (testsh.and.fh *fl <=zero) exit LoopxiMax
    !==== Zak 17/05 fast convergence ====!
    xi1=xi2
    fl=fh
    lfacl=lfach
    testsl=testsh
    !====================================!
    xi2=two*xi2
    niiter=ni
end do   LoopxiMax
!--------------------------------------------------------!

if(niiter == srmhd_maxiterationNR .or. (fh*fl>zero))then
  ! could not find upper bound for NR on xi
  print *,'could not find upper bound for NR on xi'
  print *,'niiter=',niiter,' versus srmhd_maxiterationNR=',&
          srmhd_maxiterationNR
  print *,'xi1=',xi1,' fl=',fl,' lfacl=',lfacl,' vs        dplus=',dplus
  print *,'xi2=',xi2,' fh=',fh,' lfach=',lfach,' vs tau+d+smallp=',&
          tau+d+small_pressure
  ierror = 3
  return
end if

finished=.false.
if(fl==zero)then
  xii=xi1
  lfaci=lfacl
  finished=(xi1>=dplus*lfacl.and.lfacl<=lfacmax)
endif

if(fh==zero)then
  xii=xi2
  lfaci=lfach
  finished=(xi2>=dplus*lfach.and.lfach<=lfacmax)
endif

if(finished)then
  xi=xii
  lfac=lfaci
  return
end if


xil=xi1
xih=xi2
! opedit: Try to take previous values:
if (xil .le. xi .and. xi .le. xih) then
   xii = xi
   lfaci = lfac
!   print*, 'taking previous values for initial guess, lfac, xi:', lfaci, xii
else
   xii=half*(xih+xil)    !Initialize the guess for rootfinder
   call srmhd_funcd(xii,f,df,lfaci,d,ssqr,tau,bsqr,sdotb,ierror)
!   print*, 'previous values not good for initial guess, lfac, xi:', lfaci, xii
end if

if(ierror /=0)return

er1 = 1.0d0
nit = 0
niiter=-1
xiprev=xii
lastxiok=-1.0d0
lastlfacok=-1.0d0
lastf=-1.0d0

!--- Start iteration ---!
LoopNRRMHD :  do ni=1,srmhd_maxiterationNR 
      nit = nit + 1
      if(nit>srmhd_maxiterationNR/2)then
        ! mix the last  value for convergence
        xii=half*(xii+xiprev)
        ! relax accuracy requirement
        er1=10.0d0*er1
        ! following avoids decrease of accuracy requirement 
        ! *every* iteration step beyond srmhd_maxiterationNR/2
        nit = nit - srmhd_maxiterationNR/10
      endif
      call srmhd_funcd(xii,f,df,lfaci,d,ssqr,tau,bsqr,sdotb,ierror) 
      if(ierror /=0)return
      tests=(xii>=dplus*lfaci.and.lfaci<=lfacmax)
      if(tests)then
        lastxiok=xii
        lastlfacok=lfaci
        lastf=f
      endif
      if(f*df==zero) then
         if(f==zero) then
            tests=(xii>=dplus*lfaci.and.lfaci<=lfacmax)
              !Print*,'ERROR IN ONE',' XI  ',xii,dplus*lfaci,' lfac ',lfaci,lfacmax
            if(tests)then
               exit ! zero found and consistency checks fullfilled
            else

              ierror=7
              !if(lfaci<=lfacmax .and. dabs(xii-dplus*lfaci)/xii<1.0d-5)then
                ! change in the density
              !  xii=(d+2.0d0*(1.0d0+small_pressure)*&
              !         srmhd_gamma*small_pressure/(srmhd_gamma-1.0d0))*lfaci
              !  vsqr  = 1.0d0-1.0d0/lfaci**2.d0
              !  dlfac = zero 
              !  call FuncPressure(xii,lfaci,d,ssqr,tau,dlfac,p,dpdxi)
              !  tau=xii-d-p+half*bsqr*(1.0d0+vsqr)-half*sdotb**2.0d0/xii
              !  end if
              return
            endif
         else
            print *,'stop: df becomes zero, non-monotonic function of xi!!'
            ierror=8
            return
         endif
      else
        xiprev=xii
        xii   =xii -f/df
        if(f*df>zero)then
          ! xi-iterate decreased
          ! restrict to left
          xii=max(xii,xil)
        else
         ! xi-iterate increased
         ! restrict to right
          xii=min(xii,xih)
        endif
      endif
      er=dabs(f/df)/xii
      if((er<srmhd_tolerNR*er1).or.(dabs(f/df)<srmhd_absaccnr))then
            call srmhd_funcd(xii,f,df,lfaci,d,ssqr,tau,bsqr,sdotb,ierror) 
            tests=(xii>=dplus*lfaci.and.lfaci<=lfacmax)
            if(tests)then
               exit LoopNRRMHD ! converged solution with tests ensured
            else
              ierror=7
              !if(lfaci<=lfacmax .and. dabs(xii-dplus*lfaci)/xii<1.0d-5)then
              !  ! change in the density
              !  xii=(d+lfaci*2.0d0*(1.0d0+small_pressure)*&
              !         srmhd_gamma*small_pressure/(srmhd_gamma-1.0d0))*lfaci
              !  vsqr  = 1.0d0-1.0d0/lfaci**2.0d0
              !  dlfac = zero 
              !  call FuncPressure(xii,lfaci,d,ssqr,tau,dlfac,p,dpdxi)
              !  tau=xii-d-p+half*bsqr*(1.0d0+vsqr)-half*sdotb**2.0d0/xii
              !  print*,'ierror=7:: ',xii,lfaci,d,ssqr,tau,p
              !else
              ! ierror=7
              !  end if
              return
             
            endif
      endif
      niiter=ni

enddo LoopNRRMHD

if(niiter==srmhd_maxiterationNR) then
     ! no convergence of NR for xi, although zero bracketed
     print *,'no convergence of NR for xi, although zero bracketed'
     print *,'er=',er,'srmhd_tolerNR=',srmhd_tolerNR,'er1=',er1,'df=',df,'srmhd_absaccnr=',srmhd_absaccnr
     print *,'xii=',xii,' f=',f
     print *,'brackets xil=',xil,' and xih=',xih,' with fl fh=',fl,fh
     print *,'lastxiok=',lastxiok,'lastlfacok=',lastlfacok,'lastf=',lastf
     ierror=2
     return
endif ! niiter==srmhd_maxiterationNR

 !===============================!
 ! final values for auxiliary variables are now passed to w-array
  xi=xii
  lfac=lfaci
 !===============================!
 ! we now perform some additional consistency checks
  
  if(xi<small_xi)then
    print *,'xi smaller than small_xi!!! '
    ierror=4
    return
  endif

 ! compute v*(xi+B^2)
if(ssqr .ne. zero)then
  temp=dsqrt(sum((s(:) + b(:)*sdotb/xi)*(s(:) + b(:)*sdotb/xi)))
  v=dsqrt(1.0d0-1.0d0/(lfac**2))*(s + b*sdotb/xi)/temp
  if((sum(v)**2.0d0)>1.0d0) then
    print *,'bsqr=',bsqr,'v=',v(:),'v^2=',sum(v(:)**2)
    print *,'v(xi+B^2)=',temp
    print *,'xi=',xi,'lfac=',lfac
    print *,' v^2>1!!! '
    ierror=6
    return
  end if
end if

end subroutine srmhd_con2prim
!=============================================================================
subroutine srmhd_srhdcon2prim(d,s,tau,lfac,xi,ierror)
! This subroutine srmhd_srhdcon2prim is re-writed by Z. Meliani 14/02/2018

! this is a copy of the HD iteration, where we solve for p via NR, and modified
! to give xi on output

double precision :: xi,lfac
double precision :: d,s(1:ndir),tau


integer          :: ierror

integer          :: ni,niiter
double precision :: pcurrent,pnew,pL,pR
double precision :: er,er1,ff,df,dp,v(1:ndir)
double precision :: psmall_pressure,lfac2inv,pLabs,pRabs,pprev
double precision :: s2overcubeG2rh,sqrs
double precision :: xicurrent
double precision :: oldff1,oldff2
double precision :: Nff
double precision :: pleft,pright,pnewi
integer          ::nit,n2it,ni2,ni3
double precision :: h,dhdp,rho,drhodp
!-----------------------------------------------------------------------------

ierror=0
! ierror=0 : ok
!
! ierror<>0
!
! ierror=1 : error on entry: must have D>=small_density, tau>=smallt_e
! ierror=2 : srmhd_maxiterationNR reached without convergence
! ierror=3 : final pressure value < smallp or xi<small_xi during iteration
! ierror=4 : final v^2=1 hence problem as lfac=1/0
! ierror=5 : nonmonotonic function f?
! ierror=7 : stop due to srmhd_checkNR violation

if(d<small_density .or. tau<small_e) then
  ierror=1
  return
endif

! incase input pressure is not available or random value: replace by smallp
sqrs=sum(s(:)**2.0d0)


! left and right brackets for p-range
psmall_pressure=dsqrt(sqrs)/srmhd_maxspeed-tau-d
pLabs=max(small_pressure,psmall_pressure)
pRabs=1.0d99
! start value from input
pcurrent=pLabs

er1=1.0d0
pprev=pcurrent

! Fudge Parameters
oldff1=1.0d7  ! High number
oldff2=1.0d9  ! High number bigger then oldff1
n2it = 0
nit  = 0



LoopNR:  do ni=1,srmhd_maxiterationNR
     nit = nit + 1
     !============= Controle ~1~=============!
     if(nit>srmhd_maxiterationNR/4)then
        !print *,'ni,er,p',ni,er,pcurrent
        ! mix pressure value for convergence
        pcurrent=half*(pcurrent+pprev)
        ! relax accuracy requirement
        er1=10.*er1
        nit = nit - srmhd_maxiterationNR/10
     endif
     !=======================================!

     niiter=ni  
     xicurrent=tau+d+pcurrent

     if(xicurrent<small_xi) then
       print*,'!--- srmhd/ in -- con2prim --- at SRHD case!'
       print *,'stop: too small xi iterate:',xicurrent
       print *,'for pressure iterate p',pcurrent
       print *,'pressure bracket pLabs pRabs',pLabs,pRabs
       print *,'iteration number:',ni
       print *,'values for d,s,tau,s2:',d,s,tau,sqrs
       ierror=3 
       return
     endif

     
     lfac2inv=1.0d0 - sqrs/xicurrent**2.0d0
     if(lfac2inv>zero) then
       lfac=1.0d0/dsqrt(lfac2inv)
     else
       print*,'!--- smrhd in :-- con2prim --- at SRHD case!'
       print *,'stop: negative or zero factor 1-v2:',lfac2inv
       print *,'for pressure iterate p',pcurrent
       print *,'pressure bracket pL pR',pL,pR
       print *,'absolute pressure bracket pLabs pRabs',pLabs,pRabs
       print *,'iteration number:',ni
       print *,'values for d,s,tau,s2:',d,s,tau,sqrs
       print *,'values for v,xi:',v,xicurrent
       call mpistop('is thye end ')
       ierror=4
       return
     endif
       
     s2overcubeG2rh=sqrs/(xicurrent**3.0d0)
     rho    = d/lfac 
     drhodp = -s2overcubeG2rh * d
     !== ZM calculation done using the EOS ==!
     call srmhd_get_val_h_dhdp(rho,pcurrent,drhodp,h,dhdp=dhdp)
!     call FuncEnthalpy(rho,pcurrent,tau,sqrs,xicurrent,&
!        s2overcubeG2rh,h,dhdp,ierror)
     !=======================================!   
     ff=-xicurrent*lfac2inv + h 
     df=- two*sqrs/xicurrent**2.0d0  + dhdp - lfac2inv

     if (ff*df==zero) then
        if (ff==zero) then
            exit ! zero found
        else
            print *,'stop: df becomes zero, non-monotonic f(p)!!'
            ierror=5
            return
        endif
     else 
        pnew=pcurrent-ff/df
        if (ff*df>zero) then
            ! pressure iterate has decreased
            ! restrict to left 
            pnew=max(pnew,pLabs)
        else  ! ff*df<0
            ! pressure iterate has increased
            ! restrict to right 
            pnew=min(pnew,pRabs)
        endif
     endif
        

     ! handle special case where NR incorrectly believes in convergence
     if(pnew == pLabs .and. pcurrent==pnew .and. &
        abs(ff)> srmhd_absaccnr .and. sqrs > zero)then
        pnewi=pnew
        ! try 2 higher pressure values to locate a sign change for f(p)
LoopCor:  do ni2=1,2
           pcurrent=pnewi*500.0d0
           xicurrent=tau+d+pcurrent
           v=s/xicurrent
           lfac2inv=1.0d0 - sum(v**2.0d0)
        
           if(lfac2inv>zero)then
             lfac=1.0d0/dsqrt(lfac2inv)
           else
             ierror=4
             return
           endif

           s2overcubeG2rh=-sqrs/(xicurrent**3.0d0)
           !==== Calculate enthalpy and derivative ====!
           rho=d/lfac
           call srmhd_get_val_h_dhdp(rho,pcurrent,drhodp,h)
!   call Bisection_Enthalpy(pcurrent,lfac2inv,d,s,&
!                                  tau,sqrs,xicurrent,h,ierror)
           Nff=-xicurrent*lfac2inv + h

       !== Save old value of pressure ==!
       pnewi=pcurrent
       !================================!

       !== find the interval where is the root ==!
       if(Nff * ff <=zero)then
          pnew=pcurrent
          exit LoopCor
       endif
       !=========================================!
        enddo LoopCor

        !== No possible solution, correct all including the conservatives ==!
        if( Nff*ff>zero)then
     
           ! following is in accord with trick done in smallvalues
           d   = 2.0d0*(1.0d0 + 10.0d0 * small_density) * small_density
           tau = 2.0d0*(1.0d0 + 10.0d0 * small_e) * small_e
           s   =zero
           pcurrent     = (srmhd_gamma-1.0d0)*tau
           xi=tau+d+pcurrent
           lfac = 1.0d0

           if(srmhd_checkNR)ierror=7
           ! leave the do loop here
           return
        endif
     endif
     !===============================================!
     dp=pcurrent-pnew
     er=two*dabs(dp)/(pnew+pcurrent)
     if(((er<srmhd_tolerNR*er1).or.(dabs(dp)<srmhd_absaccnr))) exit LoopNR
     !===============================================!

     ! For very small values of pressure, NR algorithm is not efficient to
     ! find root, use Euler algorithm to find precise value of pressure
     if((dabs(oldff2-ff) < 1.0d-8 .or. &
        niiter >= srmhd_maxiterationNR-srmhd_maxiterationNR/20).and.&
           ff * oldff1 < zero    .and.  dabs(ff)>srmhd_absaccnr)then

       n2it=n2it+1
       if(n2it<=3) pcurrent=half*(pnew+pcurrent)
       if(n2it>3)then
         pright =pcurrent
         pleft=pprev
         pcurrent=half*(pleft+pright)
 Dicho:  do ni3=1,srmhd_maxiterationNR
           !===================!
           xicurrent=tau+d+pcurrent
           v=s/xicurrent
           lfac2inv=1.0d0 - sum(v**2.0d0)
           if(lfac2inv>zero)then
             lfac=1.0d0/dsqrt(lfac2inv)
           else
             ierror=4
             return
           endif
           !===================!


       !== ZM calculation done using the EOS ==!
           rho=d/lfac
           call srmhd_get_val_h_dhdp(rho,pnew,drhodp,h)
!       call Bisection_Enthalpy(pnew,lfac2inv,d,s,&
!                                   tau,sqrs,xicurrent,h,ierror)
           Nff=-xicurrent*lfac2inv + h 
        !=======================================!
       !==== Iterate ====!
       if(ff * Nff < zero)then
             pleft=pcurrent
         else
        pright=pcurrent
           endif

       pcurrent=half*(pleft+pright)
       !==================!

       !=== The iteration converge ===!
           if(2.0d0*dabs(pleft-pright)/(pleft+pright)< srmhd_absaccnr &
          .or. dabs(ff)<srmhd_absaccnr)then
              pnew=pcurrent
          exit LoopNR
            endif
        !==============================!
        !======================================!
         enddo    Dicho
       endif

     else
       !====== There is no problems, continue the NR iteration ======!
       pprev=pcurrent
       pcurrent=pnew
       !=============================================================!
     endif 

 
     !=== keep the values of the 2 last ff ===!
     oldff2=oldff1
     oldff1=ff
     !========================================!
enddo LoopNR
if(niiter==srmhd_maxiterationNR)then
   ierror=2
   return
endif

if(pcurrent<small_pressure) then
   ierror=3
   return
endif

!------------------------------!
xi=tau+d+pcurrent
v = s/xi
lfac2inv=1.0d0 - sum(v**2.0d0)
if(lfac2inv>zero) then
   lfac=1.0d0/dsqrt(lfac2inv)
else
   ierror=4
   return
endif

end subroutine srmhd_srhdcon2prim
!=============================================================================
subroutine srmhd_funcd(xi,F,dF,lfac,d,ssqr,tau,bsqr,sdotb,ierror)


double precision, intent(in)  :: xi,d,ssqr,tau,bsqr,sdotb
double precision, intent(out) :: F,dF,lfac
integer, intent(inout)        :: ierror
  
double precision  :: dlfacdxi,sb2,h,dhdxi,rho,drhodxi
double precision  :: vsqr,p,dpdxi,invlfac
!-----------------------------------------------------------------------------

sb2 = sdotb*sdotb
vsqr = (ssqr + (2.0*(xi + bsqr)*sb2)/(xi**2.0))/((xi+bsqr)**2)
if (vsqr<1.0d0) then

   lfac = 1.0d0/dsqrt(1.0d0-vsqr)
   invlfac = 1.0/lfac
   rho=d/lfac
   h=xi/lfac

   dlfacdxi = -lfac**3.0d0*(sb2*(3.0d0*xi*(xi+bsqr)+bsqr*bsqr)&
               +ssqr*xi**3.0d0)&
              /((xi*(xi+bsqr))**3.0d0)

   drhodxi  = -rho*dlfacdxi*invlfac
   dhdxi    = invlfac**2.0*(1.0-2.0*xi*invlfac*dlfacdxi)

   call srmhd_get_val_p_dpdxi(rho,h,drhodxi,dhdxi,p,dpdxi)   


   F  = xi-tau-d+half*bsqr*(1.0d0+vsqr)-half*sb2/xi**2.0d0-p
   dF = 1.0d0 + sb2/xi**3.0d0  -dpdxi+dlfacdxi*bsqr/lfac**3.0d0
else 
  ! print *,'Warning: erroneous input to funcd since vsrq=',vsqr,' >=1'
   print *,'input values d, ssqr, tau, bsqr, sdotb:',d,ssqr,tau,bsqr,sdotb
   print*,'ierror ==6 ',it,mype,global_time,saveigrid 
   ierror =6
   call mpistop('at srmhd_funcd')
   return
end if

end subroutine srmhd_funcd
!=============================================================================

end module mod_srmhd_con2prim
