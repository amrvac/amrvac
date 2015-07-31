subroutine con2prim(lfac,xi,d,s^C,tau,b^C,e^C,ierror)

! (D,S,tau,B,E) --> compute auxiliaries lfac and xi

include 'amrvacdef.f'

double precision:: lfac,xi
double precision:: d,s^C,tau,b^C,e^C
integer         :: ierror
      
double precision:: ssqr,bsqr,esqr,edotb,sdotecrossb,ecrossb2,ecrossb(1:^NC),e(1:^NC),b(1:^NC)
double precision:: xi1,xi2,xii,xih,xil,dxi,f,df, tmp
double precision:: temp,vsqr,fl,fh,lfacl,lfach,lfaci
double precision:: er,er1,xiprev,dplus,lfacmax
double precision:: fv,dfv,v^C,pressure
double precision:: lastxiok,lastlfacok,lastf
double precision:: dlfac,p,dpdxi
integer:: ni,nit,niiter,i,j,k

logical:: finished,tests,testsl,testsh
!-----------------------------------------------------------------------------

ierror=0

! ierror=1 : error on entry: must have D>=0, tau>=smallp/(gamma-1)
! ierror=2 : maxitnr reached without convergence
! ierror=3 : no range for solution was find
! ierror=4 : final xi <smallxi
! ierror=5 : lower bound non-negative f(xi) value
! ierror=6 : v^2>1
! ierror=7 :  solution out from the range?
! ierror=8 : nonmonotonic function f?

if(d<minrho .or. tau<smalltau)then
  ierror=1
  return
endif

bsqr={b^C**2+}
esqr={e^C**2+}
ssqr={s^C**2+}
edotb={e^C*b^C+}
ecrossb2=esqr*bsqr - edotb**2

{^C& e(^C) = e^C;}
{^C& b(^C) = b^C;}
ecrossb(:) = zero
do i=1,^NC
   do j=1,^NC
      do k=1,^NC
         if (i .eq. j .or. i .eq. k .or. j .eq. k) cycle
         ecrossb(i) = ecrossb(i) + lvc(i,j,k)*e(j)*b(k)
      end do
   end do
end do
sdotecrossb = {s^C*ecrossb(^C)+}

! Hydro case: handle exactly as before (using p to iterate on)
{#IFDEF ENERGY
! opedit: I don't know how this would be done if we don't have energy so removing for now.
if(bsqr<=limitvalue .and. esqr<=limitvalue)then
  call con2primHydro(lfac,xi,d,s^C,tau,ierror)
  return 
end if ! Hydro
}

! the maximal allowed Lorentz factor
lfacmax=one/dsqrt(one-(one-dmaxvel)**2.0d0)
! starting point for NR on xi
! xi1 is meant to be the lower bound for the bracket to be used in NR
dplus=d+eqpar(gamma_)*minp/(eqpar(gamma_)-one)
xi1=dplus

! compute v^2(xi1)
tmp = (ssqr+ecrossb2-two*sdotecrossb)
vsqr = tmp/xi1**2

!=== find new xi1, in the case v2(xi1) > maxvel^2 = (1-dmaxvel)^2 
! locate xi corresponding to maximal velocity allowed, namely 1-dmaxvel
if(vsqr > (one-dmaxvel)**2.0d0 )then 
   xi1 = sqrt(tmp)/(1-dmaxvel)
endif 
 !=============================================!

if(niiter==maxitnr)then
 ! could not find consistent value of lower bound for xi compliant with maxvel
 print *,' could not find value of lower bound for xi compliant with maxvel'
 print *,'xi1=',xi1,'dplus=',dplus,' tau+d=',tau+d
 ierror=2
 return
endif

! we now compute f(xi1) and lfac(xi1)

call funcd(xi1,fl,df,lfacl,d,tau,ssqr,bsqr,esqr,ecrossb2,sdotecrossb,ierror)
if(ierror /=0)return
testsl=(xi1>=dplus*lfacl.and.lfacl<=lfacmax)
 
vsqr = (ssqr+ecrossb2-two*sdotecrossb)/xi1**2
if(fl>zero)then
  print *,'warning: lower bound non-negative f(xi) value!'
  print*,'iteration  = ',it
  print *,'xi1=',xi1,' vs dplus=',dplus,'dplus*lfacl=',dplus*lfacl&
         ,'lfacl=',lfacl,'lfacmax=',lfacmax
  print *,'fl=',fl,' niiter=',niiter,' vsqr(dplus)=',vsqr
  print*,' ierror = 5',it,mype,saveigrid,t,fl,df,xi1,lfacl,vsqr,ierror
  ierror=5
  return
endif
 
!--------------------------------------------------------!
! xi2 is meant to be the maximal bound for the bracket to be used in NR
! for this we take the value of E=tau+d, increased by smallp

xi2= max(tau+d+minp - half *bsqr,10.0d0*xi1)
niiter=-1

LoopxiMax : do ni=1,maxitnr
    ! we now compute f(xi2) and lfac(xi2)
    call funcd(xi2,fh,df,lfach,d,tau,ssqr,bsqr,esqr,ecrossb2,sdotecrossb,ierror)
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

if(niiter == maxitnr .or. (fh*fl>zero))then
  ! could not find upper bound for NR on xi
  print *,'could not find upper bound for NR on xi'
  print *,'niiter=',niiter,' versus maxitnr=',maxitnr
  print *,'xi1=',xi1,' fl=',fl,' lfacl=',lfacl,' vs        dplus=',dplus
  print *,'xi2=',xi2,' fh=',fh,' lfach=',lfach,' vs tau+d+smallp=',tau+d+smallp
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
   call funcd(xii,f,df,lfaci,d,tau,ssqr,bsqr,esqr,ecrossb2,sdotecrossb,ierror)
!   print*, 'previous values not good for initial guess, lfac, xi:', lfaci, xii
end if

if(ierror /=0)return

er1 = one
nit = 0
niiter=-1
xiprev=xii
lastxiok=-one
lastlfacok=-one
lastf=-one

!--- Start iteration ---!
LoopNRRMHD :  do ni=1,maxitnr 
      nit = nit + 1
      if(nit>maxitnr/2)then
        ! mix the last  value for convergence
        xii=half*(xii+xiprev)
        ! relax accuracy requirement
        er1=10.0d0*er1
        ! following avoids decrease of accuracy requirement 
        ! *every* iteration step beyond maxitnr/2
        nit = nit - maxitnr/10
      endif
      call funcd(xii,f,df,lfaci,d,tau,ssqr,bsqr,esqr,ecrossb2,sdotecrossb,ierror) 
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
            if(tests)then
               exit ! zero found and consistency checks fullfilled
            else
              !Print*,'ERROR IN ONE',' XI  ',xii,dplus*lfaci,' lfac ',lfaci,lfacmax
              ierror=7
              !if(lfaci<=lfacmax .and. dabs(xii-dplus*lfaci)/xii<1.0d-5)then
                ! change in the density
              !  xii=(d+2.0d0*(1.0d0+minp)*&
              !         eqpar(gamma_)*minp/(eqpar(gamma_)-one))*lfaci
              !  vsqr  = one-one/lfaci**2.d0
              !  dlfac = zero 
              !  call FuncPressure(xii,lfaci,d,ssqr,tau,dlfac,p,dpdxi)
              !  tau=xii-d-p+half*bsqr*(one+vsqr)-half*sdotb**2.0d0/xii
              !end if
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
      if((er<tolernr*er1).or.(dabs(f/df)<absaccnr))then
            call funcd(xii,f,df,lfaci,d,tau,ssqr,bsqr,esqr,ecrossb2,sdotecrossb,ierror) 
            tests=(xii>=dplus*lfaci.and.lfaci<=lfacmax)
            if(tests)then
               exit LoopNRRMHD ! converged solution with tests ensured
            else
              ierror=7
              !if(lfaci<=lfacmax .and. dabs(xii-dplus*lfaci)/xii<1.0d-5)then
              !  ! change in the density
              !  xii=(d+lfaci*2.0d0*(1.0d0+minp)*&
              !         eqpar(gamma_)*minp/(eqpar(gamma_)-one))*lfaci
              !  vsqr  = one-one/lfaci**2.0d0
              !  dlfac = zero 
              !  call FuncPressure(xii,lfaci,d,ssqr,tau,dlfac,p,dpdxi)
              !  tau=xii-d-p+half*bsqr*(one+vsqr)-half*sdotb**2.0d0/xii
              !  print*,'ierror=7:: ',xii,lfaci,d,ssqr,tau,p
              !else
              ! ierror=7
              !end if
              return
             
            endif
      endif
      niiter=ni

enddo LoopNRRMHD

if(niiter==maxitnr) then
     ! no convergence of NR for xi, although zero bracketed
     print *,'no convergence of NR for xi, although zero bracketed'
     print *,'er=',er,'tolernr=',tolernr,'er1=',er1,'df=',df,'absaccnr=',absaccnr
     print *,'xii=',xii,' f=',f
     print *,'brackets xil=',xil,' and xih=',xih,' with fl fh=',fl,fh
     print *,'lastxiok=',lastxiok,'lastlfacok=',lastlfacok,'lastf=',lastf
     ierror=2
     return
endif ! niiter==maxitnr

 !===============================!
 ! final values for auxiliary variables are now passed to w-array
  xi=xii
  lfac=lfaci
 !===============================!
 ! we now perform some additional consistency checks
  
  if(xi<smallxi)then
    print *,'xi smaller than smallxi!!! '
    ierror=4
    return
  endif

 ! compute vsqr for checking
vsqr = (ssqr+ecrossb2-two*sdotecrossb)/xi**2
if(vsqr>one) then
    print *,'xi=',xi,'lfac=',lfac
    print *,' v^2>1!!! '
    ierror=6
    return
end if

end subroutine con2prim
{#IFDEF ENERGY
!=============================================================================
subroutine con2primHydro(lfac,xi,d,s^C,tau,ierror)
!use ieee_arithmetic
include 'amrvacdef.f'

! this is a copy of the HD iteration, where we solve for p via NR, and modified
! to give xi on output

double precision :: xi,lfac
double precision :: d,s^C,tau
integer          :: ierror

integer          :: ni,niiter
double precision :: pcurrent,pnew,pL,pR
double precision :: er,er1,ff,df,dp,v^C
double precision :: pmin,lfac2inv,pLabs,pRabs,pprev
double precision :: s2overcubeG2rh,sqrs
double precision :: xicurrent
double precision :: oldff1,oldff2
double precision :: Nff
double precision :: pleft,pright,pnewi
integer          ::nit,n2it,ni2,ni3
double precision :: h,dhdp
!-----------------------------------------------------------------------------

ierror=0
! ierror=0 : ok
!
! ierror<>0
!
! ierror=1 : error on entry: must have D>=minrho, tau>=smalltau
! ierror=2 : maxitnr reached without convergence
! ierror=3 : final pressure value < smallp or xi<smallxi during iteration
! ierror=4 : final v^2=1 hence problem as lfac=1/0
! ierror=5 : nonmonotonic function f?
! ierror=7 : stop due to strictnr violation

if(d<minrho .or. tau<smalltau) then
  ierror=1
  return
endif

! incase input pressure is not available or random value: replace by smallp

sqrs={s^C**2.0d0+}


! left and right brackets for p-range
pmin=dsqrt(sqrs)/(one-dmaxvel)-tau-d
pLabs=max(minp,pmin)
pRabs=1.0d99
! start value from input
pcurrent=pLabs

er1=one
pprev=pcurrent

! Fudge Parameters
oldff1=1.0d7  ! High number
oldff2=1.0d9  ! High number bigger then oldff1
n2it = 0
nit  = 0



LoopNR:  do ni=1,maxitnr
     nit = nit + 1
     !============= Controle ~1~=============!
     if(nit>maxitnr/4)then
        !print *,'ni,er,p',ni,er,pcurrent
        ! mix pressure value for convergence
        pcurrent=half*(pcurrent+pprev)
        ! relax accuracy requirement
        er1=10.*er1
        nit = nit - maxitnr/10
     endif
     !=======================================!

     niiter=ni  
     xicurrent=tau+d+pcurrent

     if(xicurrent<smallxi) then
       print*,'!--- amrvacphys/t.srhd-- con2prim ---!'
       print *,'stop: too small xi iterate:',xicurrent
       print *,'for pressure iterate p',pcurrent
       print *,'pressure bracket pLabs pRabs',pLabs,pRabs
       print *,'iteration number:',ni
       print *,'values for d,s,tau,s2:',d,s^C,tau,sqrs
       ierror=3 
       return
     endif

     {v^C=s^C/xicurrent\}
     lfac2inv=one - ({v^C**2.0d0+})
     if(lfac2inv>zero) then
       lfac=one/dsqrt(lfac2inv)
     else
       print*,'!--- amrvacphys/t.srhd-- con2prim ---!'
       print *,'stop: negative or zero factor 1-v2:',lfac2inv
       print *,'for pressure iterate p',pcurrent
       print *,'pressure bracket pL pR',pL,pR
       print *,'absolute pressure bracket pLabs pRabs',pLabs,pRabs
       print *,'iteration number:',ni
       print *,'values for d,s,tau,s2:',d,s^C,tau,sqrs
       print *,'values for v,xi:',v^C,xicurrent
       ierror=4
       return
     endif
       
     s2overcubeG2rh=sqrs/(xicurrent**3.0d0)
     !== ZM calculation done using the EOS ==!
     call FuncEnthalpy(pcurrent,lfac2inv,d,s^C,tau,sqrs,xicurrent,&
		s2overcubeG2rh,h,dhdp,ierror)
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
        abs(ff)> absaccnr .and. sqrs > zero)then
        pnewi=pnew
        ! try 2 higher pressure values to locate a sign change for f(p)
LoopCor:  do ni2=1,2
	   !=====================!
	   pcurrent=pnewi*500.0d0
  	   xicurrent=tau+d+pcurrent
     	   {v^C=s^C/xicurrent\}
     	   lfac2inv=one - ({v^C**2.0d0+})
 	   !=====================!
        
	   !=====================!
    	   if(lfac2inv>zero)then
	      lfac=one/dsqrt(lfac2inv)
	   else
              ierror=4
              return
	   endif
 	   !=====================!

	   s2overcubeG2rh=-sqrs/(xicurrent**3.0d0)
	   !==== Calculate enthalpy and derivative ====!
	   call Bisection_Enthalpy(pcurrent,lfac2inv,d,s^C,&
                                   tau,sqrs,xicurrent,h,ierror)
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
           d   = 2.0d0*(one + 10.0d0 * minrho) * minrho
           tau = 2.0d0*(one + 10.0d0 * smalltau) * smalltau
           {^C&s^C =zero;}
           pcurrent     = (eqpar(gamma_)-one)*tau
           xi=tau+d+pcurrent
           lfac = one

           if(strictnr)ierror=7
           ! leave the do loop here
           return
        endif
     endif
     !===============================================!
     dp=pcurrent-pnew
     er=two*dabs(dp)/(pnew+pcurrent)
     if(((er<tolernr*er1).or.(dabs(dp)<absaccnr))) exit LoopNR
     !===============================================!

     ! For very small values of pressure, NR algorithm is not efficient to
     ! find root, use Euler algorithm to find precise value of pressure
     if((dabs(oldff2-ff) < 1.0d-8 .or. niiter >= maxitnr-maxitnr/20).and.&
           ff * oldff1 < zero    .and.  dabs(ff)>absaccnr)then

       n2it=n2it+1
       if(n2it<=3) pcurrent=half*(pnew+pcurrent)
       if(n2it>3)then
         pright =pcurrent
         pleft=pprev
         pcurrent=half*(pleft+pright)
 Dicho:  do ni3=1,maxitnr
           !===================!
           xicurrent=tau+d+pcurrent
           {v^C=s^C/xicurrent\}
           lfac2inv=one - ({v^C**2.0d0+})
           if(lfac2inv>zero)then
             lfac=one/dsqrt(lfac2inv)
           else
             ierror=4
             return
           endif
           !===================!


	   !== ZM calculation done using the EOS ==!
	   call Bisection_Enthalpy(pnew,lfac2inv,d,s^C,&
                                   tau,sqrs,xicurrent,h,ierror)
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
           if(2.0d0*dabs(pleft-pright)/(pleft+pright)< absaccnr &
	      .or. dabs(ff)<absaccnr)then
              pnew=pcurrent
	      exit LoopNR
            endif
	    !==============================!

	    !=== conserve the last value of Nff ===!
	    ff=Nff
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
  
if(niiter==maxitnr)then
   ierror=2
   return
endif

if(pcurrent<minp) then
   ierror=3
   return
endif

!------------------------------!
xi=tau+d+pcurrent
{v^C = s^C/xi\}
lfac2inv=one - ({v^C**2.0d0+})
if(lfac2inv>zero) then
   lfac=one/dsqrt(lfac2inv)
else
   ierror=4
   return
endif

end subroutine con2primHydro
}
!=============================================================================
subroutine funcd(xi,F,dF,lfac,d,tau,ssqr,bsqr,esqr,ecrossb2,sdotecrossb,ierror)

include 'amrvacdef.f'

double precision, intent(in)  :: xi,d,tau,bsqr,esqr,ssqr,ecrossb2,sdotecrossb
double precision, intent(out) :: F,dF,lfac
integer, intent(inout)        :: ierror
  
double precision  :: dlfac,sb2
double precision  :: vsqr,p,dpdxi
!-----------------------------------------------------------------------------

vsqr = (ssqr+ecrossb2-two*sdotecrossb)/xi**2

if (vsqr<one) then
   lfac  = one/dsqrt(one-vsqr)
   dlfac = - (vsqr/xi*lfac**3)

   !===== Pressure, calculate using EOS =====!
   call FuncPressure(xi,lfac,d,ssqr,tau,dlfac,p,dpdxi)
   !=========================================!
   print*,p,dpdxi

   F  = xi-tau-d + half*(esqr+bsqr)-p
   dF = one - dpdxi

else 
  ! print *,'Warning: erroneous input to funcd since vsrq=',vsqr,' >=1'
  ! print *,'input values d, ssqr, tau, bsqr, sdotb:',d,ssqr,tau,bsqr,sdotb
   print*,'ierror ==6 ',it,mype,t,saveigrid 
   ierror =6
   return
end if

end subroutine funcd
!=============================================================================
