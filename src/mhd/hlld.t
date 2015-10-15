subroutine hlld(method,qdt,ixI^L,ixO^L,idim^LIM, &
                     qtC,wCT,qt,wnew,wold,fC,dx^D,x)

! method=='hlld'  --> 2nd order HLLD scheme.
! method=='hlld1' --> 1st order HLLD scheme.
! method=='hlldd' --> 2nd order HLLD+tvdlf scheme.
! method=='hlldd1'--> 1st order HLLD+tvdlf scheme.

include 'amrvacdef.f'

character(len=*), intent(in)                         :: method
double precision, intent(in)                         :: qdt, qtC, qt, dx^D
integer, intent(in)                                  :: ixI^L, ixO^L, idim^LIM
double precision, dimension(ixI^S,1:ndim), intent(in) ::  x
double precision, dimension(ixI^S,1:ndim)             ::  xi
double precision, dimension(ixI^S,1:nw)               :: wCT, wnew, wold
double precision, dimension(ixI^S,1:nwflux,1:ndim)  :: fC

double precision, dimension(ixI^S,1:nw)            :: wLC, wRC
double precision, dimension(ixI^S)                 :: vLC, vRC
double precision, dimension(ixI^S)                 :: cmaxC,cminC

double precision, dimension(1:ndim)                :: dxinv
double precision                                   :: dxdim

integer, dimension(ixI^S)                          :: patchf
integer :: idims, iw, ix^L, hxO^L, ixC^L, jxC^L, kxC^L, kxR^L
logical :: transport, new_cmax, CmaxMeanState, logiB
logical, dimension(ixI^S) :: patchw

!=== specific to HLLD and HLLDD ===!
double precision, dimension(ixI^S,1:nwflux)     :: fLC, fRC
double precision, dimension(ixI^S,1:nwflux)     :: whll, Fhll, fD
double precision, dimension(ixI^S,-1:1)         :: lambdaD 
!-----------------------------------------------------------------------------

call mpistop('hlld is still just a dummy, not implemented...')

CmaxMeanState = (typetvdlf=='cmaxmean')
logiB=(BnormLF.and.b0_>0)

if (idimmax>idimmin .and. typelimited=='original' .and. &
   method/='hlld1' .and. method/='hlldd1')&
   call mpistop("Error in hlld: Unsplit dim. and original is limited")

! The flux calculation contracts by one in the idim direction it is applied.
! The limiter contracts the same directions by one more, so expand ixO by 2.
ix^L=ixO^L;
do idims= idim^LIM
   ix^L=ix^L^LADD2*kr(idims,^D);
end do
if (ixI^L^LTix^L|.or.|.or.) &
   call mpistop("Error in hlld : Nonconforming input limits")

if ((method/='hlld'.and.method/='hlldd1').and.useprimitive) then  
   call primitive(ixI^L,ixI^L,wCT,x)
endif 

^D&dxinv(^D)=-qdt/dx^D;
do idims= idim^LIM
   if (B0field) then
      select case (idims)
      {case (^D)
         myB0 => myB0_face^D\}
      end select
   end if

   hxO^L=ixO^L-kr(idims,^D);
   ! ixC is centered index in the idim direction from ixOmin-1/2 to ixOmax+1/2
   ixCmax^D=ixOmax^D; ixCmin^D=hxOmin^D;

   ! Calculate wRC=uR_{j+1/2} and wLC=uL_j+1/2 (eq.4.38a,b)
   jxC^L=ixC^L+kr(idims,^D);

   ! enlarged for ppm purposes
   kxCmin^D=ixImin^D; kxCmax^D=ixImax^D-kr(idims,^D); ![1,15]
   kxR^L=kxC^L+kr(idims,^D);                          ![2,16]

   wRC(kxC^S,1:nwflux)=wCT(kxR^S,1:nwflux)
   wLC(kxC^S,1:nwflux)=wCT(kxC^S,1:nwflux)

   ! Get interface positions:
   xi(kxC^S,1:ndim) = x(kxC^S,1:ndim)
   xi(kxC^S,idims) = half* ( x(kxR^S,idims)+x(kxC^S,idims) )

   ! for hlld and hlldd (second order schemes): apply limiting
   if (method=='hlld'.or.method=='hlldd') then
      dxdim=-qdt/dxinv(idims)
      select case (typelimited)
      case ('previous')
         call upwindLR(ixI^L,ixC^L,ixC^L,idims,wold,wCT,wLC,wRC,x,.true.,dxdim)
      case ('predictor')
         call upwindLR(ixI^L,ixC^L,ixC^L,idims,wCT,wCT,wLC,wRC,x,.false.,dxdim)
      case ('original')
         call upwindLR(ixI^L,ixC^L,ixC^L,idims,wnew,wCT,wLC,wRC,x,.true.,dxdim)
      case default
         call mpistop("Error in hlld: no such base for limiter")
      end select
   end if

   ! For the high order hlld scheme the limiter is based on
   ! the maximum eigenvalue, it is calculated in advance.
   if (CmaxMeanState) then
         ! determine mean state and store in wLC
         wLC(ixC^S,1:nwflux)= &
               half*(wLC(ixC^S,1:nwflux)+wRC(ixC^S,1:nwflux))
         ! get auxilaries for mean state
         if (nwaux>0) then
            call getaux(.true.,wLC,xi,ixG^LL,ixC^L,'hlld_cmaxmeanstate')
         end if
         new_cmax=.true.
         call getcmax(new_cmax,wLC,xi,ixG^LL,ixC^L,idims,cmaxC,cminC,.true.)

         ! We regain wLC for further use
         wLC(ixC^S,1:nwflux)=two*wLC(ixC^S,1:nwflux)-wRC(ixC^S,1:nwflux)
   else
         ! get auxilaries for L and R states
         if (nwaux>0.and.(.not.(useprimitive))) then
            call getaux(.true.,wLC,xi,ixG^LL,ixC^L,'hlld_wLC')
            call getaux(.true.,wRC,xi,ixG^LL,ixC^L,'hlld_wRC')
         end if
         new_cmax=.true.
         ! to save memory, use cmaxC and lambdaD for cmacRC and cmaxLC respectively
         ! to save memory, use vLC   and vRC      for cminRC and cminLC respectively
         call getcmax(new_cmax,wLC,xi,ixG^LL,ixC^L,idims,cmaxC,vLC,.true.)
         call getcmax(new_cmax,wRC,xi,ixG^LL,ixC^L,idims,cminC,vRC,.true.)
         ! now take the maximum of left and right states
         cmaxC(ixC^S)=max(cminC(ixC^S),cmaxC(ixC^S))
         cminC(ixC^S)=min(vRC(ixC^S),vLC(ixC^S))
   end if

   patchf(ixC^S) =  1
   where(cminC(ixC^S) >= zero)
        patchf(ixC^S) = -2
   elsewhere(cmaxC(ixC^S) <= zero)
        patchf(ixC^S) =  2
   endwhere

   if ((nwaux>0) .and. CmaxMeanState .and.(.not.(useprimitive))) then
      call getaux(.true.,wLC,xi,ixG^LL,ixC^L,'hlld_wLC_B')
      call getaux(.true.,wRC,xi,ixG^LL,ixC^L,'hlld_wRC_B')
   end if

   ! Calculate velocities for transport fluxes
   if(any(patchf(ixC^S)/= 2).or.(logiB)) call getv(wLC,xi,ixG^LL,ixC^L,idims,vLC)
   if(any(patchf(ixC^S)/=-2).or.(logiB)) call getv(wRC,xi,ixG^LL,ixC^L,idims,vRC)


   ! Calculate fLC=f(uL_j+1/2) and fRC=f(uR_j+1/2) for each iw
   do iw=1,nwflux
     if(any(patchf(ixC^S)/= 2).or.(logiB.and.(iw==b0_+idims{#IFDEF GLM .or.iw==psi_}))) &
       call getfluxforhllc(wLC,xi,ixG^LL,ixC^L,iw,idims,fLC,transport)
     if(any(patchf(ixC^S)/=-2).or.(logiB.and.(iw==b0_+idims{#IFDEF GLM .or.iw==psi_}))) &
       call getfluxforhllc(wRC,xi,ixG^LL,ixC^L,iw,idims,fRC,transport)
     if (transport) then
       if(any(patchf(ixC^S)/= 2).or.(logiB.and.(iw==b0_+idims{#IFDEF GLM .or.iw==psi_}))) &
         fLC(ixC^S,iw)=fLC(ixC^S,iw)+vLC(ixC^S)*wLC(ixC^S,iw)
       if(any(patchf(ixC^S)/=-2).or.(logiB.and.(iw==b0_+idims{#IFDEF GLM .or. iw==psi_}))) &
         fRC(ixC^S,iw)=fRC(ixC^S,iw)+vRC(ixC^S)*wRC(ixC^S,iw)
     end if
   end do
   ! Use more diffusive scheme, is actually TVDLF and selected by patchf=4 
   if(method=='hlldd' .or. method=='hlldd1') &
     call diffuse_hlldd(ixG^LL,ixC^L,idims,wLC,wRC,fLC,fRC,patchf)

   !---- calculate speed lambda at CD ----!
   if(any(patchf(ixC^S)==1)) &
     call getlD(wLC,wRC,fLC,fRC,cminC,cmaxC,idims,ixG^LL,ixC^L, &
            whll,Fhll,lambdaD,patchf)

   ! now patchf may be -1 or 1 due to getlD 
   if(any(abs(patchf(ixC^S))== 1))then
      !======== flux at intermediate state ========!
      call getwD(wLC,wRC,whll,x,vLC,vRC,fRC,fLC,Fhll,patchf,lambdaD,&
                  cminC,cmaxC,ixG^LL,ixC^L,idims,fD)
   endif ! Calculate the CD flux


   do iw=1,nwflux
     if (logiB.and.(iw==b0_+idims{#IFDEF GLM .or.iw==psi_})) then
       fLC(ixC^S,iw) = half*((fLC(ixC^S,iw)+fRC(ixC^S,iw)) &
                       -tvdlfeps*max(cmaxC(ixC^S)&
                       ,dabs(cminC(ixC^S)))*(wRC(ixC^S,iw)-wLC(ixC^S,iw)))
     else
       where(patchf(ixC^S)==-3)
        fLC(ixC^S,iw)=fLC(ixC^S,iw)
       elsewhere(abs(patchf(ixC^S))<=2)
        fLC(ixC^S,iw)=fD(ixC^S,iw)
       elsewhere(patchf(ixC^S)==3)
        fLC(ixC^S,iw)=fRC(ixC^S,iw)
       elsewhere(patchf(ixC^S)==3)
        ! fallback option, reducing to HLL flux
        fLC(ixC^S,iw)=Fhll(ixC^S,iw)
       elsewhere(patchf(ixC^S)==4)
        ! fallback option, reducing to TVDLF flux
        fLC(ixC^S,iw) = half*((fLC(ixC^S,iw)+fRC(ixC^S,iw)) &
          -tvdlfeps*max(cmaxC(ixC^S),dabs(cminC(ixC^S)))*(wRC(ixC^S,iw)-wLC(ixC^S,iw)))
       endwhere
     endif

     if (slab) then
         fC(ixC^S,iw,idims)=dxinv(idims)*fLC(ixC^S,iw)
         wnew(ixO^S,iw)=wnew(ixO^S,iw)+ &
              (fC(ixO^S,iw,idims)-fC(hxO^S,iw,idims))
     else
         select case (idims)
         {case (^D)
            fC(ixC^S,iw,^D)=-qdt*mygeo%surfaceC^D(ixC^S)*fLC(ixC^S,iw)
            wnew(ixO^S,iw)=wnew(ixO^S,iw)+ &
              (fC(ixO^S,iw,^D)-fC(hxO^S,iw,^D))/mygeo%dvolume(ixO^S)\}
         end select
      end if

   end do ! Next iw

end do ! Next idims

if ((method/='hlld1'.and.method/='hlldd1').and.useprimitive) then  
   patchw(ixI^S)=.false.
   call conserve(ixI^L,ixI^L,wCT,x,patchw)
endif

if (.not.slab.and.idimmin==1) call addgeometry(qdt,ixI^L,ixO^L,wCT,wnew,x)
!if (sourceunsplit) call addsource2(qdt*dble(idimmax-idimmin+1)/dble(ndim), &
!                                   ixI^L,ixO^L,1,nw,qtC,wCT,qt,wnew,x)
call addsource2(qdt*dble(idimmax-idimmin+1)/dble(ndim), &
                                   ixI^L,ixO^L,1,nw,qtC,wCT,qt,wnew,x,.false.)

end subroutine hlld
