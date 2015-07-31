!============================================================================
subroutine PPMlimitervar(ixImin1,ixImax1,ixmin1,ixmax1,idims,q,qCT,qLC,qRC)

! references:
! Mignone et al 2005, ApJS 160, 199,
! Miller and Colella 2002, JCP 183, 26
! Fryxell et al. 2000 ApJ, 131, 273 (Flash)
! baciotti Phd (http://www.aei.mpg.de/~baiotti/Baiotti_PhD.pdf)
! version : april 2009
! author: zakaria.meliani@wis.kuleuven.be

include 'amrvacdef.f'

integer, intent(in)             :: ixImin1,ixImax1, ixmin1,ixmax1, idims
double precision, intent(in)    :: q(ixImin1:ixImax1),qCT(ixImin1:ixImax1)

double precision, intent(inout) :: qRC(ixGlo1:ixGhi1),qLC(ixGlo1:ixGhi1)

double precision,dimension(ixGlo1:ixGhi1)  :: dqC,d2qC,ldq
double precision,dimension(ixGlo1:ixGhi1)  :: qMin,qMax,tmp

integer   :: lxCmin1,lxCmax1,lxRmin1,lxRmax1
integer   :: ixLLmin1,ixLLmax1,ixLmin1,ixLmax1,ixOmin1,ixOmax1,ixRmin1,&
   ixRmax1,ixRRmin1,ixRRmax1
integer   :: hxLmin1,hxLmax1,hxCmin1,hxCmax1,hxRmin1,hxRmax1
integer   :: kxLLmin1,kxLLmax1,kxLmin1,kxLmax1,kxCmin1,kxCmax1,kxRmin1,&
   kxRmax1,kxRRmin1,kxRRmax1
!--------------------------------------------------------------------------
ixOmin1=ixmin1-kr(idims,1);ixOmax1=ixmax1+kr(idims,1); !ixO[ixMmin1-1,ixMmax1+1]
ixLmin1=ixOmin1-kr(idims,1);ixLmax1=ixOmax1-kr(idims,1); !ixL[ixMmin1-2,ixMmax1]
ixLLmin1=ixLmin1-kr(idims,1);ixLLmax1=ixLmax1-kr(idims,1); !ixLL[ixMmin1-3,ixMmax1-1]
ixRmin1=ixOmin1+kr(idims,1);ixRmax1=ixOmax1+kr(idims,1); !ixR=[iMmin1,ixMmax+2]
ixRRmin1=ixRmin1+kr(idims,1);ixRRmax1=ixRmax1+kr(idims,1); !ixRR=[iMmin1+1,ixMmax+3]

hxCmin1=ixOmin1;hxCmax1=ixmax1; ! hxC = [ixMmin-1,ixMmax]
hxLmin1=hxCmin1-kr(idims,1);hxLmax1=hxCmax1-kr(idims,1); !hxL = [ixMmin-2,ixMmax-1]
hxRmin1=hxCmin1+kr(idims,1);hxRmax1=hxCmax1+kr(idims,1); !hxR = [ixMmin,ixMmax+1]

kxCmin1=ixLLmin1; kxCmax1=ixRmax1; ! kxC=[iMmin1-3,ixMmax1+2]
kxLmin1=kxCmin1-kr(idims,1);kxLmax1=kxCmax1-kr(idims,1); !kxL=[iMmin1-4,ixMmax1+1]
kxRmin1=kxCmin1+kr(idims,1);kxRmax1=kxCmax1+kr(idims,1); !kxR=[iMmin1-2,ixMmax1+3]

lxCmin1=ixLLmin1-kr(idims,1);lxCmax1=ixRRmax1;! ixC=[iMmin1-4,ixMmax1+3]
lxRmin1=lxCmin1+kr(idims,1);lxRmax1=lxCmax1+kr(idims,1); !lxR=[iMmin1-3,ixMmax1+4]


dqC(lxCmin1:lxCmax1)=q(lxRmin1:lxRmax1)-q(lxCmin1:lxCmax1)
! Eq. 64,  Miller and Colella 2002, JCP 183, 26
d2qC(kxCmin1:kxCmax1)=half*(q(kxRmin1:kxRmax1)-q(kxLmin1:kxLmax1))
where(dqC(kxCmin1:kxCmax1)*dqC(kxLmin1:kxLmax1)>zero)
   ! Store the sign of d2qC in qMin
   qMin(kxCmin1:kxCmax1)= sign(one,d2qC(kxCmin1:kxCmax1))
   ! Eq. 65,  Miller and Colella 2002, JCP 183, 26
   ldq(kxCmin1:kxCmax1)= qMin(kxCmin1:kxCmax1)*min(dabs(d2qC&
      (kxCmin1:kxCmax1)),2.0d0*dabs(dqC(kxLmin1:kxLmax1)),2.0d0&
      *dabs(dqC(kxCmin1:kxCmax1)))
elsewhere
   ldq(kxCmin1:kxCmax1)=zero
end where

! Eq. 66, Miller and Colella 2002, JCP 183, 26
qLC(ixOmin1:ixOmax1)=qLC(ixOmin1:ixOmax1)+half*dqC(ixOmin1:ixOmax1)&
   +(ldq(ixOmin1:ixOmax1)-ldq(ixRmin1:ixRmax1))/6.0d0
if (flatppm)then
   qRC(ixLmin1:ixLmax1)=qRC(ixLLmin1:ixLLmax1)+(half*dqC(ixLmin1:ixLmax1)&
      +(ldq(ixLmin1:ixLmax1)-ldq(ixOmin1:ixOmax1))/6.0d0)
else
   qRC(ixLmin1:ixLmax1)=qRC(ixLmin1:ixLmax1) -(half*dqC(ixLmin1:ixLmax1)&
      +(ldq(ixLmin1:ixLmax1)-ldq(ixOmin1:ixOmax1))/6.0d0)
endif

! make sure that min wCT(i)<wLC(i)<wCT(i+1) same for wRC(i)
call extremaq(ixImin1,ixImax1,kxCmin1,kxCmax1,qCT,1,qMax,qMin)

! Eq. B8, page 217, Mignone et al 2005, ApJS
qRC(ixLmin1:ixLmax1)=max(qMin(ixOmin1:ixOmax1),min(qMax(ixOmin1:ixOmax1),&
   qRC(ixLmin1:ixLmax1)))
qLC(ixOmin1:ixOmax1)=max(qMin(ixOmin1:ixOmax1),min(qMax(ixOmin1:ixOmax1),&
   qLC(ixOmin1:ixOmax1)))

! Eq. B9, page 217, Mignone et al 2005, ApJS
where((qRC(ixLmin1:ixLmax1)-qCT(ixOmin1:ixOmax1))*(qCT(ixOmin1:ixOmax1)&
   -qLC(ixOmin1:ixOmax1))<=zero)
   qRC(ixLmin1:ixLmax1)=qCT(ixOmin1:ixOmax1)
   qLC(ixOmin1:ixOmax1)=qCT(ixOmin1:ixOmax1)
end where

qMax(ixOmin1:ixOmax1)=(qLC(ixOmin1:ixOmax1)-qRC(ixLmin1:ixLmax1))&
   *(qCT(ixOmin1:ixOmax1)-(qLC(ixOmin1:ixOmax1)+qRC(ixLmin1:ixLmax1))/2.0d0)
qMin(ixOmin1:ixOmax1)=(qLC(ixOmin1:ixOmax1)-qRC(ixLmin1:ixLmax1))**2.0d0/6.0d0
tmp(ixLmin1:ixLmax1)=qRC(ixLmin1:ixLmax1)

! Eq. B10, page 218, Mignone et al 2005, ApJS
where(qMax(hxRmin1:hxRmax1)>qMin(hxRmin1:hxRmax1))
   qRC(hxCmin1:hxCmax1)= 3.0d0*qCT(hxRmin1:hxRmax1)-2.0d0*qLC(hxRmin1:hxRmax1)
end where
! Eq. B11, page 218, Mignone et al 2005, ApJS
where(qMax(hxCmin1:hxCmax1)<-qMin(hxCmin1:hxCmax1))
   qLC(hxCmin1:hxCmax1)= 3.0d0*qCT(hxCmin1:hxCmax1)-2.0d0*tmp(hxLmin1:hxLmax1)
end where

end subroutine PPMlimitervar
!============================================================================
subroutine PPMlimiter(ixImin1,ixImax1,ixmin1,ixmax1,idims,w,wCT,wLC,wRC)

! references:
! Mignone et al 2005, ApJS 160, 199, 
! Miller and Colella 2002, JCP 183, 26 
! Fryxell et al. 2000 ApJ, 131, 273 (Flash)
! baciotti Phd (http://www.aei.mpg.de/~baiotti/Baiotti_PhD.pdf)
! version : april 2009
! author: zakaria.meliani@wis.kuleuven.be

include 'amrvacdef.f'

integer, intent(in)             :: ixImin1,ixImax1, ixmin1,ixmax1, idims
double precision, intent(in)    :: w(ixImin1:ixImax1,1:nw),&
   wCT(ixImin1:ixImax1,1:nw)

double precision, intent(inout) :: wRC(ixGlo1:ixGhi1,1:nw),wLC(ixGlo1:ixGhi1,&
   1:nw) 

double precision,dimension(ixGlo1:ixGhi1,1:nwflux)  :: dwC,d2wC,ldw
double precision,dimension(ixGlo1:ixGhi1,1:nwflux)  :: wMin,wMax,tmp
double precision,dimension(ixGlo1:ixGhi1) :: aa, ab, ac, dv
double precision,dimension(ixGlo1:ixGhi1,1:ndim) ::  exi

integer   :: lxCmin1,lxCmax1,lxRmin1,lxRmax1
integer   :: ixLLmin1,ixLLmax1,ixLmin1,ixLmax1,ixOmin1,ixOmax1,ixRmin1,&
   ixRmax1,ixRRmin1,ixRRmax1
integer   :: hxLmin1,hxLmax1,hxCmin1,hxCmax1,hxRmin1,hxRmax1
integer   :: kxLLmin1,kxLLmax1,kxLmin1,kxLmax1,kxCmin1,kxCmax1,kxRmin1,&
   kxRmax1,kxRRmin1,kxRRmax1
integer   :: iw, idimss

double precision, parameter :: betamin=0.75d0, betamax=0.85d0,Zmin&
   =0.25d0, Zmax=0.75d0,eta1=20.0d0,eta2=0.05d0,eps=0.01d0,kappa=0.1d0
!--------------------------------------------------------------------------
ixOmin1=ixmin1-kr(idims,1);ixOmax1=ixmax1+kr(idims,1); !ixO[ixMmin1-1,ixMmax1+1]
ixLmin1=ixOmin1-kr(idims,1);ixLmax1=ixOmax1-kr(idims,1); !ixL[ixMmin1-2,ixMmax1]
ixLLmin1=ixLmin1-kr(idims,1);ixLLmax1=ixLmax1-kr(idims,1); !ixLL[ixMmin1-3,ixMmax1-1]
ixRmin1=ixOmin1+kr(idims,1);ixRmax1=ixOmax1+kr(idims,1); !ixR=[iMmin1,ixMmax+2]
ixRRmin1=ixRmin1+kr(idims,1);ixRRmax1=ixRmax1+kr(idims,1); !ixRR=[iMmin1+1,ixMmax+3]

hxCmin1=ixOmin1;hxCmax1=ixmax1; ! hxC = [ixMmin-1,ixMmax]
hxLmin1=hxCmin1-kr(idims,1);hxLmax1=hxCmax1-kr(idims,1); !hxL = [ixMmin-2,ixMmax-1]
hxRmin1=hxCmin1+kr(idims,1);hxRmax1=hxCmax1+kr(idims,1); !hxR = [ixMmin,ixMmax+1]

kxCmin1=ixLLmin1; kxCmax1=ixRmax1; ! kxC=[iMmin1-3,ixMmax1+2]
kxLmin1=kxCmin1-kr(idims,1);kxLmax1=kxCmax1-kr(idims,1); !kxL=[iMmin1-4,ixMmax1+1]
kxRmin1=kxCmin1+kr(idims,1);kxRmax1=kxCmax1+kr(idims,1); !kxR=[iMmin1-2,ixMmax1+3]

lxCmin1=ixLLmin1-kr(idims,1);lxCmax1=ixRRmax1;! ixC=[iMmin1-4,ixMmax1+3]
lxRmin1=lxCmin1+kr(idims,1);lxRmax1=lxCmax1+kr(idims,1); !lxR=[iMmin1-3,ixMmax1+4]

dwC(lxCmin1:lxCmax1,1:nwflux)=w(lxRmin1:lxRmax1,1:nwflux)-w(lxCmin1:lxCmax1,&
   1:nwflux)
! Eq. 64,  Miller and Colella 2002, JCP 183, 26 
d2wC(kxCmin1:kxCmax1,1:nwflux)=half*(w(kxRmin1:kxRmax1,1:nwflux)&
   -w(kxLmin1:kxLmax1,1:nwflux))
where(dwC(kxCmin1:kxCmax1,1:nwflux)*dwC(kxLmin1:kxLmax1,1:nwflux)>zero)
   ! Store the sign of dwC in wMin
   wMin(kxCmin1:kxCmax1,1:nwflux)= sign(one,d2wC(kxCmin1:kxCmax1,1:nwflux))
   ! Eq. 65,  Miller and Colella 2002, JCP 183, 26 
   ldw(kxCmin1:kxCmax1,1:nwflux)= wMin(kxCmin1:kxCmax1,1:nwflux)&
      *min(dabs(d2wC(kxCmin1:kxCmax1,1:nwflux)),2.0d0*dabs(dwC&
      (kxLmin1:kxLmax1,1:nwflux)),2.0d0*dabs(dwC(kxCmin1:kxCmax1,1:nwflux)))
elsewhere
   ldw(kxCmin1:kxCmax1,1:nwflux)=zero
endwhere

! Eq. 66,  Miller and Colella 2002, JCP 183, 26 
wLC(ixOmin1:ixOmax1,1:nwflux)=wLC(ixOmin1:ixOmax1,1:nwflux)&
   +half*dwC(ixOmin1:ixOmax1,1:nwflux)+(ldw(ixOmin1:ixOmax1,1:nwflux)&
   -ldw(ixRmin1:ixRmax1,1:nwflux))/6.0d0
if(flatppm)then 
   wRC(ixLmin1:ixLmax1,1:nwflux)=wRC(ixLLmin1:ixLLmax1,1:nwflux)&
      +half*dwC(ixLmin1:ixLmax1,1:nwflux)+(ldw(ixLmin1:ixLmax1,1:nwflux)&
      -ldw(ixOmin1:ixOmax1,1:nwflux))/6.0d0
else
   wRC(ixLmin1:ixLmax1,1:nwflux)=wRC(ixLmin1:ixLmax1,1:nwflux)&
      -(half*dwC(ixLmin1:ixLmax1,1:nwflux)+(ldw(ixLmin1:ixLmax1,1:nwflux)&
      -ldw(ixOmin1:ixOmax1,1:nwflux))/6.0d0)
endif
 
! make sure that min wCT(i)<wLC(i)<wCT(i+1) same for wRC(i)
call extremaw(ixImin1,ixImax1,kxCmin1,kxCmax1,wCT,1,wMax,wMin)

! Eq. B8, page 217, Mignone et al 2005, ApJS
wRC(ixLmin1:ixLmax1,1:nwflux)=max(wMin(ixOmin1:ixOmax1,1:nwflux),&
   min(wMax(ixOmin1:ixOmax1,1:nwflux),wRC(ixLmin1:ixLmax1,1:nwflux))) 
wLC(ixOmin1:ixOmax1,1:nwflux)=max(wMin(ixOmin1:ixOmax1,1:nwflux),&
   min(wMax(ixOmin1:ixOmax1,1:nwflux),wLC(ixOmin1:ixOmax1,1:nwflux)))
  

! Eq. B9, page 217, Mignone et al 2005, ApJS
where((wRC(ixLmin1:ixLmax1,1:nwflux)-wCT(ixOmin1:ixOmax1,1:nwflux))&
   *(wCT(ixOmin1:ixOmax1,1:nwflux)-wLC(ixOmin1:ixOmax1,1:nwflux))<=zero)
   wRC(ixLmin1:ixLmax1,1:nwflux)=wCT(ixOmin1:ixOmax1,1:nwflux)
   wLC(ixOmin1:ixOmax1,1:nwflux)=wCT(ixOmin1:ixOmax1,1:nwflux)
end where

wMax(ixOmin1:ixOmax1,1:nwflux)=(wLC(ixOmin1:ixOmax1,1:nwflux)&
   -wRC(ixLmin1:ixLmax1,1:nwflux))*(wCT(ixOmin1:ixOmax1,1:nwflux)&
   -(wLC(ixOmin1:ixOmax1,1:nwflux)+wRC(ixLmin1:ixLmax1,1:nwflux))/2.0d0)
wMin(ixOmin1:ixOmax1,1:nwflux)=(wLC(ixOmin1:ixOmax1,1:nwflux)&
   -wRC(ixLmin1:ixLmax1,1:nwflux))**2.0d0/6.0d0
tmp(ixLmin1:ixLmax1,1:nwflux)=wRC(ixLmin1:ixLmax1,1:nwflux)
! Eq. B10, page 218, Mignone et al 2005, ApJS
where(wMax(hxRmin1:hxRmax1,1:nwflux)>wMin(hxRmin1:hxRmax1,1:nwflux))
   wRC(hxCmin1:hxCmax1,1:nwflux)= 3.0d0*wCT(hxRmin1:hxRmax1,1:nwflux)&
      -2.0d0*wLC(hxRmin1:hxRmax1,1:nwflux)
endwhere
! Eq. B11, page 218, Mignone et al 2005, ApJS
where(wMax(hxCmin1:hxCmax1,1:nwflux)<-wMin(hxCmin1:hxCmax1,1:nwflux))
   wLC(hxCmin1:hxCmax1,1:nwflux)= 3.0d0*wCT(hxCmin1:hxCmax1,1:nwflux)&
      -2.0d0*tmp(hxLmin1:hxLmax1,1:nwflux)
endwhere

! flattening at the contact discontinuities
if(flatcd)then
 call ppmflatcd(ixImin1,ixImax1,kxCmin1,kxCmax1,kxLmin1,kxLmax1,kxRmin1,&
    kxRmax1,wCT,d2wC,aa,ab)
 if(any(kappa*aa(kxCmin1:kxCmax1)>=ab(kxCmin1:kxCmax1)))then
  do iw=1,nwflux
   where(kappa*aa(kxCmin1:kxCmax1)>=ab(kxCmin1:kxCmax1).and. &
      dabs(dwC(kxCmin1:kxCmax1,iw))>smalldouble)
      wMax(kxCmin1:kxCmax1,iw) = wCT(kxRmin1:kxRmax1,iw)-2.0d0&
         *wCT(kxCmin1:kxCmax1,iw)+wCT(kxLmin1:kxLmax1,iw)
   end where

   where(wMax(ixRmin1:ixRmax1,iw)*wMax(ixLmin1:ixLmax1,iw)<zero &
      .and.dabs(wCT(ixRmin1:ixRmax1,iw)-wCT(ixLmin1:ixLmax1,iw))&
      -eps*min(dabs(wCT(ixRmin1:ixRmax1,iw)),dabs(wCT(ixLmin1:ixLmax1,iw)))&
      >zero .and. kappa*aa(ixOmin1:ixOmax1)>=ab(ixOmin1:ixOmax1)&
      .and. dabs(dwC(ixOmin1:ixOmax1,iw))>smalldouble)

     ac(ixOmin1:ixOmax1)=(wCT(ixLLmin1:ixLLmax1,iw)-wCT(ixRRmin1:ixRRmax1,iw)&
        +4.0d0*dwC(ixOmin1:ixOmax1,iw))/(12.0d0*dwC(ixOmin1:ixOmax1,iw))
     wMin(ixOmin1:ixOmax1,iw)=max(zero,min(eta1*(ac(ixOmin1:ixOmax1)&
        -eta2),one))
   elsewhere
     wMin(ixOmin1:ixOmax1,iw)=zero
   end where

   where(wMin(hxCmin1:hxCmax1,iw)>zero)
     wLC(hxCmin1:hxCmax1,iw) = wLC(hxCmin1:hxCmax1,iw)*(one&
        -wMin(hxCmin1:hxCmax1,iw))+(wCT(hxCmin1:hxCmax1,iw)&
        +half*ldw(hxCmin1:hxCmax1,iw))*wMin(hxCmin1:hxCmax1,iw)
   end where
   where(wMin(hxRmin1:hxRmax1,iw)>zero)
     wRC(hxCmin1:hxCmax1,iw) = wRC(hxCmin1:hxCmax1,iw)*(one&
        -wMin(hxRmin1:hxRmax1,iw))+(wCT(hxRmin1:hxRmax1,iw)&
        -half*ldw(hxRmin1:hxRmax1,iw))*wMin(hxRmin1:hxRmax1,iw)
   end where
  end do
 endif
endif

! flattening at the shocks
if(flatsh)then
  ! following MILLER and COLELLA 2002 JCP 183, 26
  kxCmin1=ixmin1-2; kxCmax1=ixmax1+2; ! kxC=[ixMmin1-2,ixMmax1+2]
  do idimss=1,ndim
   kxLmin1=kxCmin1-kr(idimss,1);kxLmax1=kxCmax1-kr(idimss,1); !kxL=[ixMmin1-3,ixMmax1+1]
   kxRmin1=kxCmin1+kr(idimss,1);kxRmax1=kxCmax1+kr(idimss,1); !kxR=[ixMmin1-1,ixMmax1+3]
   kxLLmin1=kxLmin1-kr(idimss,1);kxLLmax1=kxLmax1-kr(idimss,1); !kxLL=[ixMmin-4,ixMmax]
   kxRRmin1=kxRmin1+kr(idimss,1);kxRRmax1=kxRmax1+kr(idimss,1); !kxRR=[ixMmin,ixMmax+4]

   call ppmflatsh(ixImin1,ixImax1,kxCmin1,kxCmax1,kxLLmin1,kxLLmax1,kxLmin1,&
      kxLmax1,kxRmin1,kxRmax1,kxRRmin1,kxRRmax1,idimss,wCT,aa,ab,dv)

   ! eq. B17, page 218, Mignone et al 2005, ApJS (had(Xi1))
   ac(kxCmin1:kxCmax1) = max(zero,min(one,(betamax-aa(kxCmin1:kxCmax1))&
      /(betamax-betamin)))
   ! eq. B18, page 218, Mignone et al 2005, ApJS (had(Xi1))
   ! recycling aa(ixL^S)
   where (dabs(dv(kxCmin1:kxCmax1))<smalldouble)
      aa(kxCmin1:kxCmax1) = max(ac(kxCmin1:kxCmax1), min(one,(Zmax&
         -ab(kxCmin1:kxCmax1))/(Zmax-Zmin)))
   elsewhere
      aa(kxCmin1:kxCmax1) = one
   endwhere

     call extremaa(ixImin1,ixImax1,ixOmin1,ixOmax1,aa,1,ab)
    
  enddo
  
  ! recycling wMax
  do iw=1,nwflux
     where(dabs(ab(ixOmin1:ixOmax1)-one)>smalldouble)
       wMax(ixOmin1:ixOmax1,iw) = (one-ab(ixOmin1:ixOmax1))&
          *wCT(ixOmin1:ixOmax1,iw)
     endwhere

     where(dabs(ab(hxCmin1:hxCmax1)-one)>smalldouble)
       wLC(hxCmin1:hxCmax1,iw) = ab(hxCmin1:hxCmax1)*wLC(hxCmin1:hxCmax1,iw)&
          +wMax(hxCmin1:hxCmax1,iw)
     endwhere

     where(dabs(ab(hxRmin1:hxRmax1)-one)>smalldouble)
       wRC(hxCmin1:hxCmax1,iw) = ab(hxRmin1:hxRmax1)*wRC(hxCmin1:hxCmax1,iw)&
          +wMax(hxRmin1:hxRmax1,iw)
     endwhere
  enddo
endif

end subroutine PPMlimiter
!============================================================================
subroutine MP5limiter(ixImin1,ixImax1,iLmin1,iLmax1,idims,w,wCT,wLC,wRC)
! MP5 limiter from Suresh & Huynh 1997
! Following the convention of Mignone et al. 2010.
! Needs at least three ghost cells.  Set dixB=3.

include 'amrvacdef.f'

integer, intent(in)             :: ixImin1,ixImax1, iLmin1,iLmax1, idims
double precision, intent(in)    :: w(ixImin1:ixImax1,1:nw),&
   wCT(ixImin1:ixImax1,1:nw)

double precision, intent(inout) :: wRC(ixGlo1:ixGhi1,1:nw),wLC(ixGlo1:ixGhi1,&
   1:nw) 
! .. local ..
integer                         :: iLmmin1,iLmmax1, iLmmmin1,iLmmmax1,&
    iLpmin1,iLpmax1, iLppmin1,iLppmax1, iLpppmin1,iLpppmax1
integer                         :: idmin1,idmax1, idpmin1,idpmax1, idppmin1,&
   idppmax1, idmmin1,idmmax1, iemin1,iemax1, iemmin1,iemmax1, iepmin1,iepmax1,&
    ieppmin1,ieppmax1
integer                         :: iw
double precision, dimension(ixGlo1:ixGhi1,1:nw)  :: f, fmp, fmin, fmax, ful,&
    dm4, d, fmd, flc, flim
double precision, dimension(ixGlo1:ixGhi1,1:nw)  :: wRCtmp, wLCtmp
double precision, dimension(ixGlo1:ixGhi1) :: tmp, tmp2, tmp3, a, b, c
logical, dimension(ixGlo1:ixGhi1)       :: flagL, flagR
double precision, parameter     :: eps=1.0d-12, alpha=4.0d0
!double precision                :: alpha
!----------------------------------------------------------------------------

! Variable alpha:
!alpha = float(nstep)/courantpar - one

! Left side:
! range to process:
!iLmin^D=ixmin^D-kr(idims,^D);iLmax^D=ixmax^D;

!{#IFDEF HALL
   ! For Hall, we need one more reconstructed layer since currents are computed in getflux:
   ! also add one ghost zone!
!   {iL^L=iL^L^LADD1;}
!}

! iL^L holds the indices of interfaces to reconstruct to.  Convention is that a center index holds the _right-side_ interface.  

iLmmin1=iLmin1-kr(idims,1);iLmmax1=iLmax1-kr(idims,1);
iLmmmin1=iLmmin1-kr(idims,1);iLmmmax1=iLmmax1-kr(idims,1);
iLpmin1=iLmin1+kr(idims,1);iLpmax1=iLmax1+kr(idims,1);
iLppmin1=iLpmin1+kr(idims,1);iLppmax1=iLpmax1+kr(idims,1);

f(iLmin1:iLmax1,1:nwflux) = 1.0d0/60.0d0 * (2.0d0* w(iLmmmin1:iLmmmax1,&
   1:nwflux) - 13.0d0* w(iLmmin1:iLmmax1,1:nwflux) + 47.0d0&
   * w(iLmin1:iLmax1,1:nwflux) + 27.0d0* w(iLpmin1:iLpmax1,1:nwflux) &
   - 3.0d0*  w(iLppmin1:iLppmax1,1:nwflux))

! get fmp and ful:
do iw=1,nwflux
   a(iLmin1:iLmax1) = w(iLpmin1:iLpmax1,iw)-w(iLmin1:iLmax1,iw)
   b(iLmin1:iLmax1) = alpha*(w(iLmin1:iLmax1,iw)-w(iLmmin1:iLmmax1,iw))
   call minmod(ixImin1,ixImax1,iLmin1,iLmax1,a,b,tmp)
   fmp(iLmin1:iLmax1,iw) = w(iLmin1:iLmax1,iw) + tmp(iLmin1:iLmax1)
   ful(iLmin1:iLmax1,iw) = w(iLmin1:iLmax1,iw) + b(iLmin1:iLmax1)
end do ! iw loop

! get dm4:
idmax1=iLmax1; idmin1=iLmin1-kr(idims,1);
idmmin1=idmin1-kr(idims,1);idmmax1=idmax1-kr(idims,1);
idpmin1=idmin1+kr(idims,1);idpmax1=idmax1+kr(idims,1);

iemax1=idmax1+kr(idims,1); iemin1=idmin1;
iemmin1=iemin1-kr(idims,1);iemmax1=iemax1-kr(idims,1);
iepmin1=iemin1+kr(idims,1);iepmax1=iemax1+kr(idims,1);

d(iemin1:iemax1,1:nwflux) = w(iepmin1:iepmax1,1:nwflux)-2.0d0&
   *w(iemin1:iemax1,1:nwflux)+w(iemmin1:iemmax1,1:nwflux)

do iw=1,nwflux
   a(idmin1:idmax1) = 4.0d0*d(idmin1:idmax1,iw)-d(idpmin1:idpmax1,iw)
   b(idmin1:idmax1) = 4.0d0*d(idpmin1:idpmax1,iw)-d(idmin1:idmax1,iw)
   call minmod(ixImin1,ixImax1,idmin1,idmax1,a,b,tmp)
   a(idmin1:idmax1) = d(idmin1:idmax1,iw)
   b(idmin1:idmax1) = d(idpmin1:idpmax1,iw)
   call minmod(ixImin1,ixImax1,idmin1,idmax1,a,b,tmp2)
   call minmod(ixImin1,ixImax1,idmin1,idmax1,tmp,tmp2,tmp3)
   dm4(idmin1:idmax1,iw) = tmp3(idmin1:idmax1)
end do

! get fmd:
fmd(iLmin1:iLmax1,1:nwflux) = (w(iLmin1:iLmax1,1:nwflux)+w(iLpmin1:iLpmax1,&
   1:nwflux))/2.0d0-dm4(iLmin1:iLmax1,1:nwflux)/2.0d0

!get flc: 
flc(iLmin1:iLmax1,1:nwflux) = half*(3.0d0*w(iLmin1:iLmax1,1:nwflux) &
   - w(iLmmin1:iLmmax1,1:nwflux)) + 4.0d0/3.0d0*dm4(iLmmin1:iLmmax1,1:nwflux)

fmin(iLmin1:iLmax1,1:nwflux) = max(min(w(iLmin1:iLmax1,1:nwflux),&
   w(iLpmin1:iLpmax1,1:nwflux),fmd(iLmin1:iLmax1,1:nwflux)),&
   min(w(iLmin1:iLmax1,1:nwflux),ful(iLmin1:iLmax1,1:nwflux),&
   flc(iLmin1:iLmax1,1:nwflux)))

fmax(iLmin1:iLmax1,1:nwflux) = min(max(w(iLmin1:iLmax1,1:nwflux),&
   w(iLpmin1:iLpmax1,1:nwflux),fmd(iLmin1:iLmax1,1:nwflux)),&
   max(w(iLmin1:iLmax1,1:nwflux),ful(iLmin1:iLmax1,1:nwflux),&
   flc(iLmin1:iLmax1,1:nwflux)))

do iw=1,nwflux
   a(iLmin1:iLmax1) = fmin(iLmin1:iLmax1,iw)
   b(iLmin1:iLmax1) = f(iLmin1:iLmax1,iw)
   c(iLmin1:iLmax1) = fmax(iLmin1:iLmax1,iw)
   call median(ixImin1,ixImax1,iLmin1,iLmax1,a,b,c,tmp)
   flim(iLmin1:iLmax1,iw) = tmp(iLmin1:iLmax1)
end do



! check case
where ((f(iLmin1:iLmax1,1:nwflux)-w(iLmin1:iLmax1,1:nwflux))&
   *(f(iLmin1:iLmax1,1:nwflux)-fmp(iLmin1:iLmax1,1:nwflux)) .le. eps)
   wLCtmp(iLmin1:iLmax1,1:nwflux) = f(iLmin1:iLmax1,1:nwflux)
elsewhere
   wLCtmp(iLmin1:iLmax1,1:nwflux) = flim(iLmin1:iLmax1,1:nwflux)
end where

! Right side:
! the interpolation from the right is obtained when the left-hand formula is applied to
! data mirrored about the interface.  
! thus substitute: 
! i-2 -> i+3
! i-1 -> i+2
! i   -> i+1
! i+1 -> i
! i+2 -> i-1

iLpppmin1=iLppmin1+kr(idims,1);iLpppmax1=iLppmax1+kr(idims,1);

f(iLmin1:iLmax1,1:nwflux) = 1.0d0/60.0d0 * (2.0d0* w(iLpppmin1:iLpppmax1,&
   1:nwflux) - 13.0d0* w(iLppmin1:iLppmax1,1:nwflux) + 47.0d0&
   * w(iLpmin1:iLpmax1,1:nwflux) + 27.0d0* w(iLmin1:iLmax1,1:nwflux) &
   - 3.0d0*  w(iLmmin1:iLmmax1,1:nwflux))

! get fmp and ful:
do iw=1,nwflux
   a(iLmin1:iLmax1) = w(iLmin1:iLmax1,iw)-w(iLpmin1:iLpmax1,iw)
   b(iLmin1:iLmax1) = alpha*(w(iLpmin1:iLpmax1,iw)-w(iLppmin1:iLppmax1,iw))
   call minmod(ixImin1,ixImax1,iLmin1,iLmax1,a,b,tmp)
   fmp(iLmin1:iLmax1,iw) = w(iLpmin1:iLpmax1,iw) + tmp(iLmin1:iLmax1)
   ful(iLmin1:iLmax1,iw) = w(iLpmin1:iLpmax1,iw) + b(iLmin1:iLmax1)
end do ! iw loop

! get dm4:
idmax1=iLmax1+kr(idims,1); idmin1=iLmin1;
idmmin1=idmin1-kr(idims,1);idmmax1=idmax1-kr(idims,1);
idpmin1=idmin1+kr(idims,1);idpmax1=idmax1+kr(idims,1);

iemax1=idmax1; iemin1=idmin1-kr(idims,1);
iemmin1=iemin1-kr(idims,1);iemmax1=iemax1-kr(idims,1);
iepmin1=iemin1+kr(idims,1);iepmax1=iemax1+kr(idims,1);
ieppmin1=iepmin1+kr(idims,1);ieppmax1=iepmax1+kr(idims,1);

d(iemin1:iemax1,1:nwflux) = w(iemin1:iemax1,1:nwflux)-2.0d0&
   *w(iepmin1:iepmax1,1:nwflux)+w(ieppmin1:ieppmax1,1:nwflux)

do iw=1,nwflux
   a(idmin1:idmax1) = 4.0d0*d(idmin1:idmax1,iw)-d(idmmin1:idmmax1,iw)
   b(idmin1:idmax1) = 4.0d0*d(idmmin1:idmmax1,iw)-d(idmin1:idmax1,iw)
   call minmod(ixImin1,ixImax1,idmin1,idmax1,a,b,tmp)
   a(idmin1:idmax1) = d(idmin1:idmax1,iw)
   b(idmin1:idmax1) = d(idmmin1:idmmax1,iw)
   call minmod(ixImin1,ixImax1,idmin1,idmax1,a,b,tmp2)
   call minmod(ixImin1,ixImax1,idmin1,idmax1,tmp,tmp2,tmp3)
   dm4(idmin1:idmax1,iw) = tmp3(idmin1:idmax1)
end do

! get fmd:
fmd(iLmin1:iLmax1,1:nwflux) = (w(iLmin1:iLmax1,1:nwflux)+w(iLpmin1:iLpmax1,&
   1:nwflux))/2.0d0-dm4(iLmin1:iLmax1,1:nwflux)/2.0d0

!get flc: 
flc(iLmin1:iLmax1,1:nwflux) = half*(3.0d0*w(iLpmin1:iLpmax1,1:nwflux) &
   - w(iLppmin1:iLppmax1,1:nwflux)) + 4.0d0/3.0d0*dm4(iLpmin1:iLpmax1,&
   1:nwflux)

fmin(iLmin1:iLmax1,1:nwflux) = max(min(w(iLpmin1:iLpmax1,1:nwflux),&
   w(iLmin1:iLmax1,1:nwflux),fmd(iLmin1:iLmax1,1:nwflux)),min(w&
   (iLpmin1:iLpmax1,1:nwflux),ful(iLmin1:iLmax1,1:nwflux),flc(iLmin1:iLmax1,&
   1:nwflux)))

fmax(iLmin1:iLmax1,1:nwflux) = min(max(w(iLpmin1:iLpmax1,1:nwflux),&
   w(iLmin1:iLmax1,1:nwflux),fmd(iLmin1:iLmax1,1:nwflux)),max(w&
   (iLpmin1:iLpmax1,1:nwflux),ful(iLmin1:iLmax1,1:nwflux),flc(iLmin1:iLmax1,&
   1:nwflux)))

do iw=1,nwflux
   a(iLmin1:iLmax1) = fmin(iLmin1:iLmax1,iw)
   b(iLmin1:iLmax1) = f(iLmin1:iLmax1,iw)
   c(iLmin1:iLmax1) = fmax(iLmin1:iLmax1,iw)
   call median(ixImin1,ixImax1,iLmin1,iLmax1,a,b,c,tmp)
   flim(iLmin1:iLmax1,iw) = tmp(iLmin1:iLmax1)
end do

! check case
where ((f(iLmin1:iLmax1,1:nwflux)-w(iLpmin1:iLpmax1,1:nwflux))&
   *(f(iLmin1:iLmax1,1:nwflux)-fmp(iLmin1:iLmax1,1:nwflux))  .le. eps)
   wRCtmp(iLmin1:iLmax1,1:nwflux) = f(iLmin1:iLmax1,1:nwflux)
elsewhere
   wRCtmp(iLmin1:iLmax1,1:nwflux) = flim(iLmin1:iLmax1,1:nwflux)
end where

! Since limiter not TVD, negative pressures or densities could result.  
! Fall back to flat interpolation (minmod would also work). 
call checkw(useprimitive,ixGlo1,ixGhi1,iLmin1,iLmax1,wLCtmp,flagL)
call checkw(useprimitive,ixGlo1,ixGhi1,iLmin1,iLmax1,wRCtmp,flagR)

do iw=1,nwflux
   where (flagL(iLmin1:iLmax1).and.flagR(iLmin1:iLmax1))
      wLC(iLmin1:iLmax1,iw)=wLCtmp(iLmin1:iLmax1,iw)
      wRC(iLmin1:iLmax1,iw)=wRCtmp(iLmin1:iLmax1,iw)
   end where
end do


end subroutine MP5limiter
!============================================================================
subroutine MP5limiterL(ixImin1,ixImax1,iLmin1,iLmax1,idims,w,wLC)
! MP5 limiter from Suresh & Huynh 1997
! Following the convention of Mignone et al. 2010.
! Needs at least three ghost cells.  Set dixB=3.

include 'amrvacdef.f'

integer, intent(in)             :: ixImin1,ixImax1, iLmin1,iLmax1, idims
double precision, intent(in)    :: w(ixImin1:ixImax1,1:nw)

double precision, intent(inout) :: wLC(ixGlo1:ixGhi1,1:nw) 
! .. local ..
integer                         :: iLmmin1,iLmmax1, iLmmmin1,iLmmmax1,&
    iLpmin1,iLpmax1, iLppmin1,iLppmax1
integer                         :: idmin1,idmax1, idpmin1,idpmax1, idppmin1,&
   idppmax1, idmmin1,idmmax1, iemin1,iemax1, iemmin1,iemmax1, iepmin1,iepmax1,&
    ieppmin1,ieppmax1
integer                         :: iw
double precision, dimension(ixGlo1:ixGhi1,1:nw)  :: f, fmp, fmin, fmax, ful,&
    dm4, d, fmd, flc, flim
double precision, dimension(ixGlo1:ixGhi1) :: tmp, tmp2, tmp3, a, b, c
double precision, parameter     :: eps=1.0d-12, alpha=4.0d0
!double precision                :: alpha
!----------------------------------------------------------------------------

! Variable alpha:
!alpha = float(nstep)/courantpar - one

! Left side:


iLmmin1=iLmin1-kr(idims,1);iLmmax1=iLmax1-kr(idims,1);
iLmmmin1=iLmmin1-kr(idims,1);iLmmmax1=iLmmax1-kr(idims,1);
iLpmin1=iLmin1+kr(idims,1);iLpmax1=iLmax1+kr(idims,1);
iLppmin1=iLpmin1+kr(idims,1);iLppmax1=iLpmax1+kr(idims,1);

f(iLmin1:iLmax1,1:nwflux) = 1.0d0/60.0d0 * (2.0d0* w(iLmmmin1:iLmmmax1,&
   1:nwflux) - 13.0d0* w(iLmmin1:iLmmax1,1:nwflux) + 47.0d0&
   * w(iLmin1:iLmax1,1:nwflux) + 27.0d0* w(iLpmin1:iLpmax1,1:nwflux) &
   - 3.0d0*  w(iLppmin1:iLppmax1,1:nwflux))

! get fmp and ful:
do iw=1,nwflux
   a(iLmin1:iLmax1) = w(iLpmin1:iLpmax1,iw)-w(iLmin1:iLmax1,iw)
   b(iLmin1:iLmax1) = alpha*(w(iLmin1:iLmax1,iw)-w(iLmmin1:iLmmax1,iw))
   call minmod(ixImin1,ixImax1,iLmin1,iLmax1,a,b,tmp)
   fmp(iLmin1:iLmax1,iw) = w(iLmin1:iLmax1,iw) + tmp(iLmin1:iLmax1)
   ful(iLmin1:iLmax1,iw) = w(iLmin1:iLmax1,iw) + b(iLmin1:iLmax1)
end do ! iw loop

! get dm4:
idmax1=iLmax1; idmin1=iLmin1-kr(idims,1);
idmmin1=idmin1-kr(idims,1);idmmax1=idmax1-kr(idims,1);
idpmin1=idmin1+kr(idims,1);idpmax1=idmax1+kr(idims,1);

iemax1=idmax1+kr(idims,1); iemin1=idmin1;
iemmin1=iemin1-kr(idims,1);iemmax1=iemax1-kr(idims,1);
iepmin1=iemin1+kr(idims,1);iepmax1=iemax1+kr(idims,1);

d(iemin1:iemax1,1:nwflux) = w(iepmin1:iepmax1,1:nwflux)-2.0d0&
   *w(iemin1:iemax1,1:nwflux)+w(iemmin1:iemmax1,1:nwflux)

do iw=1,nwflux
   a(idmin1:idmax1) = 4.0d0*d(idmin1:idmax1,iw)-d(idpmin1:idpmax1,iw)
   b(idmin1:idmax1) = 4.0d0*d(idpmin1:idpmax1,iw)-d(idmin1:idmax1,iw)
   call minmod(ixImin1,ixImax1,idmin1,idmax1,a,b,tmp)
   a(idmin1:idmax1) = d(idmin1:idmax1,iw)
   b(idmin1:idmax1) = d(idpmin1:idpmax1,iw)
   call minmod(ixImin1,ixImax1,idmin1,idmax1,a,b,tmp2)
   call minmod(ixImin1,ixImax1,idmin1,idmax1,tmp,tmp2,tmp3)
   dm4(idmin1:idmax1,iw) = tmp3(idmin1:idmax1)
end do

! get fmd:
fmd(iLmin1:iLmax1,1:nwflux) = (w(iLmin1:iLmax1,1:nwflux)+w(iLpmin1:iLpmax1,&
   1:nwflux))/2.0d0-dm4(iLmin1:iLmax1,1:nwflux)/2.0d0

!get flc: 
flc(iLmin1:iLmax1,1:nwflux) = half*(3.0d0*w(iLmin1:iLmax1,1:nwflux) &
   - w(iLmmin1:iLmmax1,1:nwflux)) + 4.0d0/3.0d0*dm4(iLmmin1:iLmmax1,1:nwflux)

fmin(iLmin1:iLmax1,1:nwflux) = max(min(w(iLmin1:iLmax1,1:nwflux),&
   w(iLpmin1:iLpmax1,1:nwflux),fmd(iLmin1:iLmax1,1:nwflux)),&
   min(w(iLmin1:iLmax1,1:nwflux),ful(iLmin1:iLmax1,1:nwflux),&
   flc(iLmin1:iLmax1,1:nwflux)))

fmax(iLmin1:iLmax1,1:nwflux) = min(max(w(iLmin1:iLmax1,1:nwflux),&
   w(iLpmin1:iLpmax1,1:nwflux),fmd(iLmin1:iLmax1,1:nwflux)),&
   max(w(iLmin1:iLmax1,1:nwflux),ful(iLmin1:iLmax1,1:nwflux),&
   flc(iLmin1:iLmax1,1:nwflux)))

do iw=1,nwflux
   a(iLmin1:iLmax1) = fmin(iLmin1:iLmax1,iw)
   b(iLmin1:iLmax1) = f(iLmin1:iLmax1,iw)
   c(iLmin1:iLmax1) = fmax(iLmin1:iLmax1,iw)
   call median(ixImin1,ixImax1,iLmin1,iLmax1,a,b,c,tmp)
   flim(iLmin1:iLmax1,iw) = tmp(iLmin1:iLmax1)
end do


! check case
where ((f(iLmin1:iLmax1,1:nwflux)-w(iLmin1:iLmax1,1:nwflux))&
   *(f(iLmin1:iLmax1,1:nwflux)-fmp(iLmin1:iLmax1,1:nwflux)) .le. eps)
   wLC(iLmin1:iLmax1,1:nwflux) = f(iLmin1:iLmax1,1:nwflux)
elsewhere
   wLC(iLmin1:iLmax1,1:nwflux) = flim(iLmin1:iLmax1,1:nwflux)
end where


end subroutine MP5limiterL
!============================================================================
subroutine MP5limiterR(ixImin1,ixImax1,iLmin1,iLmax1,idims,w,wRC)
! MP5 limiter from Suresh & Huynh 1997
! Following the convention of Mignone et al. 2010.
! Needs at least three ghost cells.  Set dixB=3.

include 'amrvacdef.f'

integer, intent(in)             :: ixImin1,ixImax1, iLmin1,iLmax1, idims
double precision, intent(in)    :: w(ixImin1:ixImax1,1:nw)

double precision, intent(inout) :: wRC(ixGlo1:ixGhi1,1:nw)
! .. local ..
integer                         :: iLmmin1,iLmmax1, iLpmin1,iLpmax1, iLppmin1,&
   iLppmax1, iLpppmin1,iLpppmax1
integer                         :: idmin1,idmax1, idpmin1,idpmax1, idppmin1,&
   idppmax1, idmmin1,idmmax1, iemin1,iemax1, iemmin1,iemmax1, iepmin1,iepmax1,&
    ieppmin1,ieppmax1
integer                         :: iw
double precision, dimension(ixGlo1:ixGhi1,1:nw)  :: f, fmp, fmin, fmax, ful,&
    dm4, d, fmd, flc, flim
double precision, dimension(ixGlo1:ixGhi1) :: tmp, tmp2, tmp3, a, b, c
double precision, parameter     :: eps=1.0d-12, alpha=4.0d0
!double precision                :: alpha
!----------------------------------------------------------------------------
! Right side:
! the interpolation from the right is obtained when the left-hand formula is applied to
! data mirrored about the interface.  
! thus substitute: 
! i-2 -> i+3
! i-1 -> i+2
! i   -> i+1
! i+1 -> i
! i+2 -> i-1

iLmmin1=iLmin1-kr(idims,1);iLmmax1=iLmax1-kr(idims,1);
iLpmin1=iLmin1+kr(idims,1);iLpmax1=iLmax1+kr(idims,1);
iLppmin1=iLpmin1+kr(idims,1);iLppmax1=iLpmax1+kr(idims,1);
iLpppmin1=iLppmin1+kr(idims,1);iLpppmax1=iLppmax1+kr(idims,1);

f(iLmin1:iLmax1,1:nwflux) = 1.0d0/60.0d0 * (2.0d0* w(iLpppmin1:iLpppmax1,&
   1:nwflux) - 13.0d0* w(iLppmin1:iLppmax1,1:nwflux) + 47.0d0&
   * w(iLpmin1:iLpmax1,1:nwflux) + 27.0d0* w(iLmin1:iLmax1,1:nwflux) &
   - 3.0d0*  w(iLmmin1:iLmmax1,1:nwflux))

! get fmp and ful:
do iw=1,nwflux
   a(iLmin1:iLmax1) = w(iLmin1:iLmax1,iw)-w(iLpmin1:iLpmax1,iw)
   b(iLmin1:iLmax1) = alpha*(w(iLpmin1:iLpmax1,iw)-w(iLppmin1:iLppmax1,iw))
   call minmod(ixImin1,ixImax1,iLmin1,iLmax1,a,b,tmp)
   fmp(iLmin1:iLmax1,iw) = w(iLpmin1:iLpmax1,iw) + tmp(iLmin1:iLmax1)
   ful(iLmin1:iLmax1,iw) = w(iLpmin1:iLpmax1,iw) + b(iLmin1:iLmax1)
end do ! iw loop

! get dm4:
idmax1=iLmax1+kr(idims,1); idmin1=iLmin1;
idmmin1=idmin1-kr(idims,1);idmmax1=idmax1-kr(idims,1);
idpmin1=idmin1+kr(idims,1);idpmax1=idmax1+kr(idims,1);

iemax1=idmax1; iemin1=idmin1-kr(idims,1);
iemmin1=iemin1-kr(idims,1);iemmax1=iemax1-kr(idims,1);
iepmin1=iemin1+kr(idims,1);iepmax1=iemax1+kr(idims,1);
ieppmin1=iepmin1+kr(idims,1);ieppmax1=iepmax1+kr(idims,1);

d(iemin1:iemax1,1:nwflux) = w(iemin1:iemax1,1:nwflux)-2.0d0&
   *w(iepmin1:iepmax1,1:nwflux)+w(ieppmin1:ieppmax1,1:nwflux)

do iw=1,nwflux
   a(idmin1:idmax1) = 4.0d0*d(idmin1:idmax1,iw)-d(idmmin1:idmmax1,iw)
   b(idmin1:idmax1) = 4.0d0*d(idmmin1:idmmax1,iw)-d(idmin1:idmax1,iw)
   call minmod(ixImin1,ixImax1,idmin1,idmax1,a,b,tmp)
   a(idmin1:idmax1) = d(idmin1:idmax1,iw)
   b(idmin1:idmax1) = d(idmmin1:idmmax1,iw)
   call minmod(ixImin1,ixImax1,idmin1,idmax1,a,b,tmp2)
   call minmod(ixImin1,ixImax1,idmin1,idmax1,tmp,tmp2,tmp3)
   dm4(idmin1:idmax1,iw) = tmp3(idmin1:idmax1)
end do

! get fmd:
fmd(iLmin1:iLmax1,1:nwflux) = (w(iLmin1:iLmax1,1:nwflux)+w(iLpmin1:iLpmax1,&
   1:nwflux))/2.0d0-dm4(iLmin1:iLmax1,1:nwflux)/2.0d0

!get flc: 
flc(iLmin1:iLmax1,1:nwflux) = half*(3.0d0*w(iLpmin1:iLpmax1,1:nwflux) &
   - w(iLppmin1:iLppmax1,1:nwflux)) + 4.0d0/3.0d0*dm4(iLpmin1:iLpmax1,&
   1:nwflux)

fmin(iLmin1:iLmax1,1:nwflux) = max(min(w(iLpmin1:iLpmax1,1:nwflux),&
   w(iLmin1:iLmax1,1:nwflux),fmd(iLmin1:iLmax1,1:nwflux)),min(w&
   (iLpmin1:iLpmax1,1:nwflux),ful(iLmin1:iLmax1,1:nwflux),flc(iLmin1:iLmax1,&
   1:nwflux)))

fmax(iLmin1:iLmax1,1:nwflux) = min(max(w(iLpmin1:iLpmax1,1:nwflux),&
   w(iLmin1:iLmax1,1:nwflux),fmd(iLmin1:iLmax1,1:nwflux)),max(w&
   (iLpmin1:iLpmax1,1:nwflux),ful(iLmin1:iLmax1,1:nwflux),flc(iLmin1:iLmax1,&
   1:nwflux)))

do iw=1,nwflux
   a(iLmin1:iLmax1) = fmin(iLmin1:iLmax1,iw)
   b(iLmin1:iLmax1) = f(iLmin1:iLmax1,iw)
   c(iLmin1:iLmax1) = fmax(iLmin1:iLmax1,iw)
   call median(ixImin1,ixImax1,iLmin1,iLmax1,a,b,c,tmp)
   flim(iLmin1:iLmax1,iw) = tmp(iLmin1:iLmax1)
end do

! check case
where ((f(iLmin1:iLmax1,1:nwflux)-w(iLpmin1:iLpmax1,1:nwflux))&
   *(f(iLmin1:iLmax1,1:nwflux)-fmp(iLmin1:iLmax1,1:nwflux))  .le. eps)
   wRC(iLmin1:iLmax1,1:nwflux) = f(iLmin1:iLmax1,1:nwflux)
elsewhere
   wRC(iLmin1:iLmax1,1:nwflux) = flim(iLmin1:iLmax1,1:nwflux)
end where

end subroutine Mp5limiterR
!============================================================================
subroutine minmod(ixImin1,ixImax1,ixOmin1,ixOmax1,a,b,minm)

include 'amrvacdef.f'

integer, intent(in)          :: ixImin1,ixImax1, ixOmin1,ixOmax1
double precision, intent(in) :: a(ixImin1:ixImax1), b(ixImin1:ixImax1)
double precision, intent(out):: minm(ixGlo1:ixGhi1)

minm(ixOmin1:ixOmax1) = (sign(one,a(ixOmin1:ixOmax1))+sign(one,&
   b(ixOmin1:ixOmax1)))/2.0d0 * min(abs(a(ixOmin1:ixOmax1)),&
   abs(b(ixOmin1:ixOmax1)))

end subroutine minmod
!============================================================================
subroutine median(ixImin1,ixImax1,ixOmin1,ixOmax1,a,b,c,med)

include 'amrvacdef.f'

integer, intent(in)          :: ixImin1,ixImax1, ixOmin1,ixOmax1
double precision, intent(in) :: a(ixImin1:ixImax1), b(ixImin1:ixImax1),&
    c(ixImin1:ixImax1)
double precision, intent(out):: med(ixGlo1:ixGhi1)
! .. local ..
double precision             :: tmp1(ixImin1:ixImax1),tmp2(ixImin1:ixImax1)

tmp1(ixOmin1:ixOmax1) = b(ixOmin1:ixOmax1) - a(ixOmin1:ixOmax1)
tmp2(ixOmin1:ixOmax1) = c(ixOmin1:ixOmax1) - a(ixOmin1:ixOmax1)

med(ixOmin1:ixOmax1) = a(ixOmin1:ixOmax1) + (sign(one,tmp1(ixOmin1:ixOmax1))&
   +sign(one,tmp2(ixOmin1:ixOmax1)))/2.0d0 * min(abs(tmp1(ixOmin1:ixOmax1)),&
   abs(tmp2(ixOmin1:ixOmax1)))

end subroutine median
!=============================================================================
subroutine hancock(qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,idimmin,idimmax,qtC,&
   wCT,qt,wnew,dx1,x)

! The non-conservative Hancock predictor for TVDLF

! on entry:
! input available on ixI^L=ixG^L asks for output on ixO^L=ixG^L^LSUBdixB

! one entry: (predictor): wCT -- w_n        wnew -- w_n   qdt=dt/2

! on exit :  (predictor): wCT -- w_n        wnew -- w_n+1/2


! FCT not implemented here

include 'amrvacdef.f'

integer, intent(in) :: ixImin1,ixImax1, ixOmin1,ixOmax1, idimmin,idimmax
double precision, intent(in) :: qdt, qtC, qt, dx1, x(ixImin1:ixImax1,1:ndim)
double precision, intent(inout) :: wCT(ixImin1:ixImax1,1:nw),&
    wnew(ixImin1:ixImax1,1:nw)

double precision, dimension(ixGlo1:ixGhi1,1:nw) :: wLC, wRC
double precision, dimension(ixGlo1:ixGhi1) :: fLC, fRC
double precision, dimension(ixGlo1:ixGhi1) :: vLC, vRC
double precision :: dxinv(1:ndim),dxdim(1:ndim)
integer :: idims, iw, ixmin1,ixmax1, hxOmin1,hxOmax1, ixtestmin1,ixtestmax1
logical :: transport
logical, dimension(ixGlo1:ixGhi1) :: patchw
!-----------------------------------------------------------------------------

oktest=.false.
if(oktest.and.mype==0) then
   print *,'======Hancock predictor: qdt, qtC, qt:',qdt,qtC,qt
   ixtestmin1=ixImin1;ixtestmax1=ixImax1;
!   ixtestmin1=400;ixtestmax1=400+2*dixB;
   print *,'reporting in ranges:',ixtestmin1,ixtestmax1
endif

! Expand limits in each idims direction in which fluxes are added
ixmin1=ixOmin1;ixmax1=ixOmax1;
do idims= idimmin,idimmax
   ixmin1=ixmin1-kr(idims,1);ixmax1=ixmax1+kr(idims,1);
end do
if (ixImin1>ixmin1.or.ixImax1<ixmax1) call mpistop&
   ("Error in Hancock: Nonconforming input limits")

if (useprimitive) then  
   if(oktest.and.mype==0) then
      do iw=1,nw 
         print *,'iw,wCT before useprimitive:',iw,wCT(ixtestmin1:ixtestmax1,iw)
      enddo
   endif
   call primitive(ixImin1,ixImax1,ixImin1,ixImax1,wCT,x)
endif 

if(oktest.and.mype==0) then
    do iw=1,nw 
       print *,'iw,wCT: before dimensional loop',iw,wCT(ixtestmin1:ixtestmax1,iw)
    enddo
endif

dxinv(1)=-qdt/dx1;
dxdim(1)=dx1;
do idims= idimmin,idimmax
   if (B0field) then
      select case (idims)
      case (1)
         myB0 => myB0_face1
      end select
   end if

   ! Calculate w_j+g_j/2 and w_j-g_j/2
   ! First copy all variables, then upwind wLC and wRC.
   ! wLC is to the left of ixO, wRC is to the right of wCT.
   hxOmin1=ixOmin1-kr(idims,1);hxOmax1=ixOmax1-kr(idims,1);

   wRC(hxOmin1:hxOmax1,1:nwflux)=wCT(ixOmin1:ixOmax1,1:nwflux)
   wLC(ixOmin1:ixOmax1,1:nwflux)=wCT(ixOmin1:ixOmax1,1:nwflux)

   if(oktest.and.mype==0) then
    do iw=1,nw 
       print *,'iw,wRC: before upwindLR',iw,wRC(ixtestmin1:ixtestmax1,iw)
    enddo
    do iw=1,nw 
       print *,'iw,wLC: before upwindLR',iw,wLC(ixtestmin1:ixtestmax1,iw)
    enddo
   endif

   call upwindLR(ixImin1,ixImax1,ixOmin1,ixOmax1,hxOmin1,hxOmax1,idims,wCT,&
      wCT,wLC,wRC,x,.false.,dxdim(idims))

   if(oktest.and.mype==0) then
    do iw=1,nw 
       print *,'iw,wRC: after upwindLR',iw,wRC(ixtestmin1:ixtestmax1,iw)
    enddo
    do iw=1,nw 
       print *,'iw,wLC: after upwindLR',iw,wLC(ixtestmin1:ixtestmax1,iw)
    enddo
   endif

   if (nwaux>0.and.(.not.(useprimitive))) then
   !!if (nwaux>0) then
      call getaux(.true.,wLC,x,ixGlo1,ixGhi1,ixOmin1,ixOmax1,'hancock_wLC')
      call getaux(.true.,wRC,x,ixGlo1,ixGhi1,hxOmin1,hxOmax1,'hancock_wRC')
   end if

   ! Calculate vLC and vRC velocities
   call getv(wRC,x,ixGlo1,ixGhi1,hxOmin1,hxOmax1,idims,vRC)
   call getv(wLC,x,ixGlo1,ixGhi1,ixOmin1,ixOmax1,idims,vLC)

   ! Advect w(iw)
   do iw=1,nwflux
      ! Calculate the fLC and fRC fluxes
      call getflux(wRC,x,ixGlo1,ixGhi1,hxOmin1,hxOmax1,iw,idims,fRC,transport)
      call getflux(wLC,x,ixGlo1,ixGhi1,ixOmin1,ixOmax1,iw,idims,fLC,transport)
      if (transport) then
         fRC(hxOmin1:hxOmax1)=fRC(hxOmin1:hxOmax1)+vRC(hxOmin1:hxOmax1)&
            *wRC(hxOmin1:hxOmax1,iw)
         fLC(ixOmin1:ixOmax1)=fLC(ixOmin1:ixOmax1)+vLC(ixOmin1:ixOmax1)&
            *wLC(ixOmin1:ixOmax1,iw)
      end if

      if(oktest.and.mype==0) then
           print *,'iw: updating wnew',iw,wnew(ixtestmin1:ixtestmax1,iw)
           print *,'iw: updating wnew with fluxes fLC',fLC(ixtestmin1:ixtestmax1)
           print *,'iw: updating wnew with fluxes fRC',fRC(ixtestmin1:ixtestmax1)
           if(.not.slab) then
             print *,'volumes:',mygeo%dvolume(ixtestmin1:ixtestmax1)
             print *,'surfaceC1:',mygeo%surfaceC1(ixtestmin1:ixtestmax1)
           endif
      endif

      ! Advect w(iw)
      if (slab) then
         wnew(ixOmin1:ixOmax1,iw)=wnew(ixOmin1:ixOmax1,iw)+dxinv(idims)&
            * (fLC(ixOmin1:ixOmax1)-fRC(hxOmin1:hxOmax1))
      else
         select case (idims)
         case (1)
            wnew(ixOmin1:ixOmax1,iw)=wnew(ixOmin1:ixOmax1,iw)&
               -qdt/mygeo%dvolume(ixOmin1:ixOmax1) &
                  *(mygeo%surfaceC1(ixOmin1:ixOmax1)*fLC(ixOmin1:ixOmax1) &
                   -mygeo%surfaceC1(hxOmin1:hxOmax1)*fRC(hxOmin1:hxOmax1))
         end select
      end if
   end do
end do ! next idims

if (useprimitive) then  
    patchw(ixImin1:ixImax1)=.false.
    call conserve(ixImin1,ixImax1,ixImin1,ixImax1,wCT,x,patchw)
    if(oktest.and.mype==0) then
      do iw=1,nw 
         print *,'iw,wCT after conserve:',iw,wCT(ixtestmin1:ixtestmax1,iw)
      enddo
    endif
else
   if(nwaux>0) call getaux(.true.,wCT,x,ixImin1,ixImax1,ixImin1,ixImax1,&
      'hancock_wCT')
endif

if (.not.slab.and.idimmin==1) call addgeometry(qdt,ixImin1,ixImax1,ixOmin1,&
   ixOmax1,wCT,wnew,x)
call addsource2(qdt*dble(idimmax-idimmin+1)/dble(ndim), ixImin1,ixImax1,&
   ixOmin1,ixOmax1,1,nw,qtC,wCT,qt,wnew,x,.false.)

end subroutine hancock
!============================================================================
subroutine upwindLR(ixImin1,ixImax1,ixLmin1,ixLmax1,ixRmin1,ixRmax1,idims,w,&
   wCT,wLC,wRC,x,needprim,dxdim)

! Determine the upwinded wLC(ixL) and wRC(ixR) from w. 
! the wCT is only used when PPM is exploited.

include 'amrvacdef.f'

integer, intent(in) :: ixImin1,ixImax1, ixLmin1,ixLmax1, ixRmin1,ixRmax1,&
    idims
logical, intent(in) :: needprim
double precision, intent(in) :: dxdim
double precision, dimension(ixImin1:ixImax1,1:nw) :: w, wCT
double precision, dimension(ixGlo1:ixGhi1,1:nw) :: wLC, wRC
double precision, dimension(ixGlo1:ixGhi1,1:ndim) :: x

integer :: jxRmin1,jxRmax1, ixCmin1,ixCmax1, jxCmin1,jxCmax1, iw, ixtestmin1,&
   ixtestmax1
double precision :: wLtmp(ixGlo1:ixGhi1,1:nw), wRtmp(ixGlo1:ixGhi1,1:nw)
double precision :: ldw(ixGlo1:ixGhi1), dwC(ixGlo1:ixGhi1)
logical, dimension(ixGlo1:ixGhi1) :: flagL, flagR
logical, dimension(ixGlo1:ixGhi1) :: patchw, patchwLC, patchwRC

character*79 :: savetypelimiter
!-----------------------------------------------------------------------------

oktest=.false.
if(oktest.and.mype==0) then
   print *,'======in upwindLR'
   ixtestmin1=ixImin1;ixtestmax1=ixImax1;
!   ixtestmin1=400;ixtestmax1=400+2*dixB;
   print *,'reporting in ranges:',ixtestmin1,ixtestmax1
endif

! Transform w,wL,wR to primitive variables
if (needprim.and.useprimitive) then
   call primitive(ixImin1,ixImax1,ixImin1,ixImax1,w,x)
end if

if(typelimiter/='ppm' .and. typelimiter /= 'mp5')then
 jxRmin1=ixRmin1+kr(idims,1);jxRmax1=ixRmax1+kr(idims,1);
 ixCmax1=jxRmax1; ixCmin1=ixLmin1-kr(idims,1);
 jxCmin1=ixCmin1+kr(idims,1);jxCmax1=ixCmax1+kr(idims,1);

 do iw=1,nwflux
   if (loglimit(iw)) then
      w(ixCmin1:jxCmax1,iw)=dlog10(w(ixCmin1:jxCmax1,iw))
      wLC(ixLmin1:ixLmax1,iw)=dlog10(wLC(ixLmin1:ixLmax1,iw))
      wRC(ixRmin1:ixRmax1,iw)=dlog10(wRC(ixRmin1:ixRmax1,iw))
   end if

   dwC(ixCmin1:ixCmax1)=w(jxCmin1:jxCmax1,iw)-w(ixCmin1:ixCmax1,iw)

   savetypelimiter=typelimiter
   if(savetypelimiter=='koren') typelimiter='korenL'
   if(savetypelimiter=='cada')  typelimiter='cadaL'
   if(savetypelimiter=='cada3') typelimiter='cada3L'
   call dwlimiter2(dwC,ixCmin1,ixCmax1,iw,idims,ldw,dxdim)

   wLtmp(ixLmin1:ixLmax1,iw)=wLC(ixLmin1:ixLmax1,iw)+half*ldw(ixLmin1:ixLmax1)
   if(savetypelimiter=='koren')then
     typelimiter='korenR'
     call dwlimiter2(dwC,ixCmin1,ixCmax1,iw,idims,ldw,dxdim)
   endif
   if(savetypelimiter=='cada')then
     typelimiter='cadaR'
     call dwlimiter2(dwC,ixCmin1,ixCmax1,iw,idims,ldw,dxdim)
   endif
   if(savetypelimiter=='cada3')then
     typelimiter='cada3R'
     call dwlimiter2(dwC,ixCmin1,ixCmax1,iw,idims,ldw,dxdim)
   endif
   wRtmp(ixRmin1:ixRmax1,iw)=wRC(ixRmin1:ixRmax1,iw)-half*ldw(jxRmin1:jxRmax1)
   typelimiter=savetypelimiter

   if (loglimit(iw)) then
      w(ixCmin1:jxCmax1,iw)=10.0d0**w(ixCmin1:jxCmax1,iw)
      wLtmp(ixLmin1:ixLmax1,iw)=10.0d0**wLtmp(ixLmin1:ixLmax1,iw)
      wRtmp(ixRmin1:ixRmax1,iw)=10.0d0**wRtmp(ixRmin1:ixRmax1,iw)
   end if
 end do

 call checkw(useprimitive,ixGlo1,ixGhi1,ixLmin1,ixLmax1,wLtmp,flagL)
 call checkw(useprimitive,ixGlo1,ixGhi1,ixRmin1,ixRmax1,wRtmp,flagR)

 do iw=1,nwflux
   where (flagL(ixLmin1:ixLmax1).and.flagR(ixRmin1:ixRmax1))
      wLC(ixLmin1:ixLmax1,iw)=wLtmp(ixLmin1:ixLmax1,iw)
      wRC(ixRmin1:ixRmax1,iw)=wRtmp(ixRmin1:ixRmax1,iw)
   end where

   if (loglimit(iw)) then
      where (.not.(flagL(ixLmin1:ixLmax1).and.flagR(ixRmin1:ixRmax1)))
         wLC(ixLmin1:ixLmax1,iw)=10.0d0**wLC(ixLmin1:ixLmax1,iw)
         wRC(ixRmin1:ixRmax1,iw)=10.0d0**wRC(ixRmin1:ixRmax1,iw)
      end where
   end if
 enddo
else if (typelimiter .eq. 'ppm') then
 call PPMlimiter(ixImin1,ixImax1,ixMlo1,ixMhi1,idims,w,wCT,wLC,wRC)
else
 call MP5limiter(ixGlo1,ixGhi1,ixLmin1,ixLmax1,idims,w,wCT,wLC,wRC)
endif

! Transform w,wL,wR back to conservative variables
if (useprimitive) then
   if(needprim)then 
      patchw(ixImin1:ixImax1)=.false.
      call conserve(ixImin1,ixImax1,ixImin1,ixImax1,w,x,patchw)
   endif
   patchwLC(ixGlo1:ixGhi1)=.false.
   patchwRC(ixGlo1:ixGhi1)=.false.
   call conserve(ixGlo1,ixGhi1,ixLmin1,ixLmax1,wLC,x,patchwLC)
   call conserve(ixGlo1,ixGhi1,ixRmin1,ixRmax1,wRC,x,patchwRC)
end if

end subroutine upwindLR
!============================================================================
subroutine dwlimiter2(dwC,ixCmin1,ixCmax1,iw,idims,ldw,dxdim)

! Limit the centered dwC differences within ixC for iw in direction idim.
! The limiter is chosen according to typelimiter.

! Note that this subroutine is called from upwindLR (hence from methods 
! like tvdlf, hancock, hll(c) etc) or directly from tvd.t,
! but also from the gradientS and divvectorS subroutines in geometry.t
! Accordingly, the typelimiter here corresponds to one of typelimiter1
! or one of typegradlimiter1.

! note: there is no iw dependence here...

include 'amrvacdef.f'

integer, intent(in) :: ixCmin1,ixCmax1, iw, idims
double precision, intent(in) :: dxdim
double precision, intent(in) :: dwC(ixGlo1:ixGhi1)
double precision, intent(out) :: ldw(ixGlo1:ixGhi1)

double precision :: tmp(ixGlo1:ixGhi1)
integer :: ixOmin1,ixOmax1, hxOmin1,hxOmax1
double precision, parameter :: qsmall=1.d-12, qsmall2=2.d-12

! cada limiter parameter values
double precision, parameter :: cadalfa=0.5d0, cadbeta=2.0d0, cadgamma=1.6d0
! full third order cada limiter
double precision :: rdelinv
double precision :: ldwA(ixGlo1:ixGhi1),ldwB(ixGlo1:ixGhi1),&
   tmpeta(ixGlo1:ixGhi1)
double precision, parameter :: cadepsilon=1.d-14, invcadepsilon&
   =1.d14,cadradius=0.1d0
!-----------------------------------------------------------------------------

! Contract indices in idim for output.
ixOmin1=ixCmin1+kr(idims,1); ixOmax1=ixCmax1;
hxOmin1=ixOmin1-kr(idims,1);hxOmax1=ixOmax1-kr(idims,1);

! Store the sign of dwC in tmp
tmp(ixOmin1:ixOmax1)=sign(one,dwC(ixOmin1:ixOmax1))
rdelinv=one/(cadradius*dxdim)**2

select case (typelimiter)
case ('minmod')
   ! Minmod limiter eq(3.51e) and (eq.3.38e) with omega=1
   ldw(ixOmin1:ixOmax1)=tmp(ixOmin1:ixOmax1)* max(zero,min(dabs(dwC&
      (ixOmin1:ixOmax1)),tmp(ixOmin1:ixOmax1)*dwC(hxOmin1:hxOmax1)))
case ('woodward')
   ! Woodward and Collela limiter (eq.3.51h), a factor of 2 is pulled out
   ldw(ixOmin1:ixOmax1)=two*tmp(ixOmin1:ixOmax1)* max(zero,&
      min(dabs(dwC(ixOmin1:ixOmax1)),tmp(ixOmin1:ixOmax1)*dwC&
      (hxOmin1:hxOmax1),tmp(ixOmin1:ixOmax1)*quarter*(dwC(hxOmin1:hxOmax1)&
      +dwC(ixOmin1:ixOmax1))))
case ('mcbeta')
   ! Woodward and Collela limiter, with factor beta
   ldw(ixOmin1:ixOmax1)=tmp(ixOmin1:ixOmax1)* max(zero,min(mcbeta&
      *dabs(dwC(ixOmin1:ixOmax1)),mcbeta*tmp(ixOmin1:ixOmax1)&
      *dwC(hxOmin1:hxOmax1),tmp(ixOmin1:ixOmax1)*half*(dwC(hxOmin1:hxOmax1)&
      +dwC(ixOmin1:ixOmax1))))
case ('superbee')
   ! Roes superbee limiter (eq.3.51i)
   ldw(ixOmin1:ixOmax1)=tmp(ixOmin1:ixOmax1)* max(zero,min(two&
      *dabs(dwC(ixOmin1:ixOmax1)),tmp(ixOmin1:ixOmax1)*dwC(hxOmin1:hxOmax1)),&
      min(dabs(dwC(ixOmin1:ixOmax1)),two*tmp(ixOmin1:ixOmax1)&
      *dwC(hxOmin1:hxOmax1)))
case ('vanleer')
  ! van Leer limiter (eq 3.51f), but a missing delta2=1.D-12 is added
  ldw(ixOmin1:ixOmax1)=two*max(dwC(hxOmin1:hxOmax1)*dwC(ixOmin1:ixOmax1),&
     zero) /(dwC(ixOmin1:ixOmax1)+dwC(hxOmin1:hxOmax1)+qsmall)
case ('albada')
  ! Albada limiter (eq.3.51g) with delta2=1D.-12
  ldw(ixOmin1:ixOmax1)=(dwC(hxOmin1:hxOmax1)*(dwC(ixOmin1:ixOmax1)**2+qsmall)&
     +dwC(ixOmin1:ixOmax1)*(dwC(hxOmin1:hxOmax1)**2+qsmall))&
     /(dwC(ixOmin1:ixOmax1)**2+dwC(hxOmin1:hxOmax1)**2+qsmall2)
case ('korenR')
   ! Barry Koren Right variant
   ldw(ixOmin1:ixOmax1)=tmp(ixOmin1:ixOmax1)* max(zero,min(two&
      *dabs(dwC(ixOmin1:ixOmax1)),two*tmp(ixOmin1:ixOmax1)*dwC&
      (hxOmin1:hxOmax1),(two*dwC(hxOmin1:hxOmax1)*tmp(ixOmin1:ixOmax1)&
      +dabs(dwC(ixOmin1:ixOmax1)))*third))
case ('korenL')
   ! Barry Koren Left variant
   ldw(ixOmin1:ixOmax1)=tmp(ixOmin1:ixOmax1)* max(zero,min(two&
      *dabs(dwC(ixOmin1:ixOmax1)),two*tmp(ixOmin1:ixOmax1)*dwC&
      (hxOmin1:hxOmax1),(dwC(hxOmin1:hxOmax1)*tmp(ixOmin1:ixOmax1)&
      +two*dabs(dwC(ixOmin1:ixOmax1)))*third))
case ('cadaR')
   ! Cada Right variant
   ldw(ixOmin1:ixOmax1)=tmp(ixOmin1:ixOmax1)* max(zero,min((two&
      *dwC(hxOmin1:hxOmax1)*tmp(ixOmin1:ixOmax1)+dabs(dwC(ixOmin1:ixOmax1)))&
      *third, max(-cadalfa*dabs(dwC(ixOmin1:ixOmax1)),                     &
      min(cadbeta*dabs(dwC(ixOmin1:ixOmax1)),                  (two&
      *dwC(hxOmin1:hxOmax1)*tmp(ixOmin1:ixOmax1)+dabs(dwC(ixOmin1:ixOmax1)))&
      *third, cadgamma*tmp(ixOmin1:ixOmax1)*dwC(hxOmin1:hxOmax1)))))
case ('cadaL')
   ! Cada Left variant
   ldw(ixOmin1:ixOmax1)=tmp(ixOmin1:ixOmax1)* max(zero,min((two&
      *dabs(dwC(ixOmin1:ixOmax1))+tmp(ixOmin1:ixOmax1)*dwC(hxOmin1:hxOmax1))&
      *third, max(-cadalfa*tmp(ixOmin1:ixOmax1)*dwC(hxOmin1:hxOmax1),&
                           min(cadbeta*tmp(ixOmin1:ixOmax1)&
      *dwC(hxOmin1:hxOmax1),                  (two*dabs(dwC(ixOmin1:ixOmax1))&
      +tmp(ixOmin1:ixOmax1)*dwC(hxOmin1:hxOmax1))*third, cadgamma&
      *dabs(dwC(ixOmin1:ixOmax1))))))
case ('cada3R')
   tmpeta(ixOmin1:ixOmax1)=(dwC(ixOmin1:ixOmax1)**2+dwC(hxOmin1:hxOmax1)&
      **2)*rdelinv
   ldwA(ixOmin1:ixOmax1)=(two*dwC(hxOmin1:hxOmax1)+dwC(ixOmin1:ixOmax1))*third
   ldwB(ixOmin1:ixOmax1)=tmp(ixOmin1:ixOmax1)* max(zero,min((two&
      *dwC(hxOmin1:hxOmax1)*tmp(ixOmin1:ixOmax1)+dabs(dwC(ixOmin1:ixOmax1)))&
      *third, max(-cadalfa*dabs(dwC(ixOmin1:ixOmax1)),                     &
      min(cadbeta*dabs(dwC(ixOmin1:ixOmax1)),                  (two&
      *dwC(hxOmin1:hxOmax1)*tmp(ixOmin1:ixOmax1)+dabs(dwC(ixOmin1:ixOmax1)))&
      *third, cadgamma*tmp(ixOmin1:ixOmax1)*dwC(hxOmin1:hxOmax1)))))
   where(tmpeta(ixOmin1:ixOmax1)<=one-cadepsilon)
     ldw(ixOmin1:ixOmax1)=ldwA(ixOmin1:ixOmax1)
   elsewhere(tmpeta(ixOmin1:ixOmax1)>=one+cadepsilon)
     ldw(ixOmin1:ixOmax1)=ldwB(ixOmin1:ixOmax1)
   elsewhere
     tmp(ixOmin1:ixOmax1)=(tmpeta(ixOmin1:ixOmax1)-one)*invcadepsilon
     ldw(ixOmin1:ixOmax1)=half*( (one-tmp(ixOmin1:ixOmax1))&
        *ldwA(ixOmin1:ixOmax1) +(one+tmp(ixOmin1:ixOmax1))*ldwB&
        (ixOmin1:ixOmax1))
   endwhere
case ('cada3L')
   tmpeta(ixOmin1:ixOmax1)=(dwC(ixOmin1:ixOmax1)**2+dwC(hxOmin1:hxOmax1)&
      **2)*rdelinv
   ldwA(ixOmin1:ixOmax1)=(two*dwC(ixOmin1:ixOmax1)+dwC(hxOmin1:hxOmax1))*third
   ldwB(ixOmin1:ixOmax1)=tmp(ixOmin1:ixOmax1)* max(zero,min((two&
      *dabs(dwC(ixOmin1:ixOmax1))+tmp(ixOmin1:ixOmax1)*dwC(hxOmin1:hxOmax1))&
      *third, max(-cadalfa*tmp(ixOmin1:ixOmax1)*dwC(hxOmin1:hxOmax1),&
                           min(cadbeta*tmp(ixOmin1:ixOmax1)&
      *dwC(hxOmin1:hxOmax1),                  (two*dabs(dwC(ixOmin1:ixOmax1))&
      +tmp(ixOmin1:ixOmax1)*dwC(hxOmin1:hxOmax1))*third, cadgamma&
      *dabs(dwC(ixOmin1:ixOmax1))))))
   where(tmpeta(ixOmin1:ixOmax1)<=one-cadepsilon)
     ldw(ixOmin1:ixOmax1)=ldwA(ixOmin1:ixOmax1)
   elsewhere(tmpeta(ixOmin1:ixOmax1)>=one+cadepsilon)
     ldw(ixOmin1:ixOmax1)=ldwB(ixOmin1:ixOmax1)
   elsewhere
     tmp(ixOmin1:ixOmax1)=(tmpeta(ixOmin1:ixOmax1)-one)*invcadepsilon
     ldw(ixOmin1:ixOmax1)=half*( (one-tmp(ixOmin1:ixOmax1))&
        *ldwA(ixOmin1:ixOmax1) +(one+tmp(ixOmin1:ixOmax1))*ldwB&
        (ixOmin1:ixOmax1))
   endwhere
case default
   write(*,*)'Unknown limiter:',typelimiter
   call mpistop("Error in dwLimiter: No such TVD limiter")
end select

end subroutine dwlimiter2
!=============================================================================
subroutine tvdlf(method,qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,idimmin,idimmax,&
    qtC,wCT,qt,wnew,wold,fC,dx1,x)

! method=='tvdmu'  --> 2nd order (may be 3rd order in 1D) TVD-MUSCL scheme.
! method=='tvdmu1' --> 1st order TVD-MUSCL scheme (upwind per charact. var.)
! method=='tvdlf'  --> 2nd order TVD-Lax-Friedrich scheme.
! method=='tvdlf1' --> 1st order TVD-Lax-Friedrich scheme.

include 'amrvacdef.f'

character(len=*), intent(in)                         :: method
double precision, intent(in)                         :: qdt, qtC, qt, dx1
integer, intent(in)                                  :: ixImin1,ixImax1,&
    ixOmin1,ixOmax1, idimmin,idimmax
double precision, dimension(ixImin1:ixImax1,1:ndim), intent(in) ::  x
double precision, dimension(ixImin1:ixImax1,1:ndim)             ::  xi
double precision, dimension(ixImin1:ixImax1,1:nw)               :: wCT, wnew,&
    wold
double precision, dimension(ixImin1:ixImax1,1:nwflux,1:ndim)        :: fC

double precision, dimension(ixGlo1:ixGhi1,1:nw) :: wLC, wRC
double precision, dimension(ixGlo1:ixGhi1)      :: fLC, fRC, vLC, vRC
double precision, dimension(ixGlo1:ixGhi1)      :: cmaxC, cmaxRC, cmaxLC,ie,&
   ieL,ieR,ieCT
double precision :: dxinv(1:ndim),dxdim(1:ndim)
integer :: idims, iw, ixmin1,ixmax1, hxOmin1,hxOmax1, ixCmin1,ixCmax1,&
    ixCRmin1,ixCRmax1, kxCmin1,kxCmax1, kxRmin1,kxRmax1, ixtestmin1,ixtestmax1
logical :: transport, new_cmax, CmaxMeanState
logical, dimension(ixGlo1:ixGhi1) :: patchw
logical, save :: evolvepressure=.false.
!-----------------------------------------------------------------------------

CmaxMeanState = (typetvdlf=='cmaxmean')

oktest=.false.
if(oktest.and.mype==0) then
   print *,'======tvdlf: qdt, qtC, qt:',qdt,qtC,qt
   ixtestmin1=ixImin1;ixtestmax1=ixImax1;
   ixtestmin1=1;ixtestmax1=2*dixB;
   print *,'reporting in ranges:',ixtestmin1,ixtestmax1
endif

!!call mpistop("tijdelijke stop")

if (idimmax>idimmin .and. typelimited=='original' .and. method&
   /='tvdlf1')call mpistop&
   ("Error in TVDMUSCLF: Unsplit dim. and original is limited")

! The flux calculation contracts by one in the idim direction it is applied.
! The limiter contracts the same directions by one more, so expand ixO by 2.
ixmin1=ixOmin1;ixmax1=ixOmax1;
do idims= idimmin,idimmax
   ixmin1=ixmin1-2*kr(idims,1);ixmax1=ixmax1+2*kr(idims,1);
end do
if (ixImin1>ixmin1.or.ixImax1<ixmax1) call mpistop&
   ("Error in TVDLF: Nonconforming input limits")





if ((method=='tvdlf').and.useprimitive) then  
   ! second order methods with primitive limiting: 
   ! this call ensures wCT is primitive with updated auxiliaries
   call primitive(ixImin1,ixImax1,ixImin1,ixImax1,wCT,x)
endif 


dxinv(1)=-qdt/dx1;
dxdim(1)=dx1;
do idims= idimmin,idimmax
   if (B0field) then
      select case (idims)
      case (1)
         myB0 => myB0_face1
      end select
   end if

   hxOmin1=ixOmin1-kr(idims,1);hxOmax1=ixOmax1-kr(idims,1);
   ! ixC is centered index in the idim direction from ixOmin-1/2 to ixOmax+1/2
!   ixCmax^D=ixOmax^D; ixCmin^D=hxOmin^D;


   ixCmax1=ixOmax1; ixCmin1=hxOmin1;


   kxCmin1=ixImin1; kxCmax1=ixImax1-kr(idims,1);
   kxRmin1=kxCmin1+kr(idims,1);kxRmax1=kxCmax1+kr(idims,1);

   ixCRmin1=ixCmin1;ixCRmax1=ixCmax1;

   wRC(kxCmin1:kxCmax1,1:nwflux)=wCT(kxRmin1:kxRmax1,1:nwflux)
   wLC(kxCmin1:kxCmax1,1:nwflux)=wCT(kxCmin1:kxCmax1,1:nwflux)

   ! Get interface positions:
   xi(kxCmin1:kxCmax1,idims) = half* ( x(kxRmin1:kxRmax1,idims)&
      +x(kxCmin1:kxCmax1,idims) )

   ! for tvdlf (second order scheme): apply limiting
   if (method=='tvdlf') then
      select case (typelimited)
      case ('previous')
         call upwindLR(ixImin1,ixImax1,ixCRmin1,ixCRmax1,ixCRmin1,ixCRmax1,&
            idims,wold,wCT,wLC,wRC,x,.true.,dxdim(idims))
      case ('predictor')
         call upwindLR(ixImin1,ixImax1,ixCRmin1,ixCRmax1,ixCRmin1,ixCRmax1,&
            idims,wCT,wCT,wLC,wRC,x,.false.,dxdim(idims))
      case ('original')
         call upwindLR(ixImin1,ixImax1,ixCRmin1,ixCRmax1,ixCRmin1,ixCRmax1,&
            idims,wnew,wCT,wLC,wRC,x,.true.,dxdim(idims))
      case default
         call mpistop("Error in TVDMUSCLF: no such base for limiter")
      end select
   end if

   ! For the high order Lax-Friedrich TVDLF scheme the limiter is based on
   ! the maximum eigenvalue, it is calculated in advance.
   if (CmaxMeanState) then
      ! determine mean state and store in wLC
      wLC(ixCmin1:ixCmax1,1:nwflux)= half*(wLC(ixCmin1:ixCmax1,1:nwflux)&
         +wRC(ixCmin1:ixCmax1,1:nwflux))
      ! get auxilaries for mean state
      if (nwaux>0) then
         call getaux(.true.,wLC,xi,ixGlo1,ixGhi1,ixCmin1,ixCmax1,&
            'tvdlf_cmaxmeanstate')
      end if
      new_cmax=.true.
      call getcmax(new_cmax,wLC,xi,ixGlo1,ixGhi1,ixCmin1,ixCmax1,idims,cmaxC,&
         cmaxLC,.false.)
      
      ! We regain wLC for further use
      wLC(ixCmin1:ixCmax1,1:nwflux)=two*wLC(ixCmin1:ixCmax1,1:nwflux)&
         -wRC(ixCmin1:ixCmax1,1:nwflux)
      if (nwaux>0) then
         call getaux(.true.,wLC,xi,ixGlo1,ixGhi1,ixCmin1,ixCmax1,&
            'tvdlf_wLC_A')
      endif
      if (nwaux>0.and.(.not.(useprimitive).or.method=='tvdlf1')) then
         call getaux(.true.,wRC,xi,ixGlo1,ixGhi1,ixCmin1,ixCmax1,&
            'tvdlf_wRC_A')
      end if
   else
      ! get auxilaries for L and R states
      if (nwaux>0.and.(.not.(useprimitive).or.method=='tvdlf1')) then
         !!if (nwaux>0) then
         call getaux(.true.,wLC,xi,ixGlo1,ixGhi1,ixCmin1,ixCmax1,'tvdlf_wLC')
         call getaux(.true.,wRC,xi,ixGlo1,ixGhi1,ixCmin1,ixCmax1,'tvdlf_wRC')
      end if
      new_cmax=.true.
      call getcmax(new_cmax,wLC,xi,ixGlo1,ixGhi1,ixCmin1,ixCmax1,idims,cmaxLC,&
         cmaxC,.false.)
      call getcmax(new_cmax,wRC,xi,ixGlo1,ixGhi1,ixCmin1,ixCmax1,idims,cmaxRC,&
         cmaxC,.false.)
      ! now take the maximum of left and right states
      cmaxC(ixCmin1:ixCmax1)=max(cmaxRC(ixCmin1:ixCmax1),cmaxLC&
         (ixCmin1:ixCmax1))
   end if
   
   ! Calculate velocities for transport fluxes
   call getv(wLC,xi,ixGlo1,ixGhi1,ixCmin1,ixCmax1,idims,vLC)
   call getv(wRC,xi,ixGlo1,ixGhi1,ixCmin1,ixCmax1,idims,vRC)





   ! Calculate fLC=f(uL_j+1/2) and fRC=f(uR_j+1/2) for each iw
   do iw=1,nwflux
      call getflux(wLC,xi,ixGlo1,ixGhi1,ixCmin1,ixCmax1,iw,idims,fLC,&
         transport)
      call getflux(wRC,xi,ixGlo1,ixGhi1,ixCmin1,ixCmax1,iw,idims,fRC,&
         transport)
      if (transport) then
         fLC(ixCmin1:ixCmax1)=fLC(ixCmin1:ixCmax1)+vLC(ixCmin1:ixCmax1)&
            *wLC(ixCmin1:ixCmax1,iw)
         fRC(ixCmin1:ixCmax1)=fRC(ixCmin1:ixCmax1)+vRC(ixCmin1:ixCmax1)&
            *wRC(ixCmin1:ixCmax1,iw)
      end if
      ! To save memory we use fLC to store (F_L+F_R)/2=half*(fLC+fRC)
      fLC(ixCmin1:ixCmax1)=half*(fLC(ixCmin1:ixCmax1)+fRC(ixCmin1:ixCmax1))

      ! Add TVDLF dissipation to the flux
      ! To save memory we use fRC to store -cmax*half*(w_R-w_L)
      fRC(ixCmin1:ixCmax1)=-tvdlfeps*cmaxC(ixCmin1:ixCmax1)*half&
         *(wRC(ixCmin1:ixCmax1,iw)-wLC(ixCmin1:ixCmax1,iw))
      ! fLC contains physical+dissipative fluxes
      fLC(ixCmin1:ixCmax1)=fLC(ixCmin1:ixCmax1)+fRC(ixCmin1:ixCmax1)

      if (slab) then
         fC(ixCmin1:ixCmax1,iw,idims)=fLC(ixCmin1:ixCmax1)
      else
         select case (idims)
         case (1)
            fC(ixCmin1:ixCmax1,iw,1)=mygeo%surfaceC1(ixCmin1:ixCmax1)&
               *fLC(ixCmin1:ixCmax1)
         end select
      end if

   end do ! Next iw

end do ! Next idims



!Now update the state:
do idims= idimmin,idimmax
   hxOmin1=ixOmin1-kr(idims,1);hxOmax1=ixOmax1-kr(idims,1);
   do iw=1,nwflux

      ! Multiply the fluxes by -dt/dx since Flux fixing expects this
      if (slab) then
         fC(ixImin1:ixImax1,iw,idims)=dxinv(idims)*fC(ixImin1:ixImax1,iw,&
            idims)
         wnew(ixOmin1:ixOmax1,iw)=wnew(ixOmin1:ixOmax1,iw) &
            + (fC(ixOmin1:ixOmax1,iw,idims)-fC(hxOmin1:hxOmax1,iw,idims))
      else
         select case (idims)
         case (1)
            fC(ixImin1:ixImax1,iw,1)=-qdt*fC(ixImin1:ixImax1,iw,idims)
            wnew(ixOmin1:ixOmax1,iw)=wnew(ixOmin1:ixOmax1,iw) &
              + (fC(ixOmin1:ixOmax1,iw,1)-fC(hxOmin1:hxOmax1,iw,1))&
                 /mygeo%dvolume(ixOmin1:ixOmax1)
         end select
      end if

   end do ! Next iw

end do ! Next idims




if ((method=='tvdlf').and.useprimitive) then  
    patchw(ixImin1:ixImax1)=.false.
    call conserve(ixImin1,ixImax1,ixImin1,ixImax1,wCT,x,patchw)
else
   if(nwaux>0) call getaux(.true.,wCT,x,ixImin1,ixImax1,ixImin1,ixImax1,&
      'tvdlf_wCT')
endif

if (.not.slab.and.idimmin==1) call addgeometry(qdt,ixImin1,ixImax1,ixOmin1,&
   ixOmax1,wCT,wnew,x)
call addsource2(qdt*dble(idimmax-idimmin+1)/dble(ndim),ixImin1,ixImax1,&
   ixOmin1,ixOmax1,1,nw,qtC,wCT,qt,wnew,x,.false.)

end subroutine tvdlf
!=============================================================================
subroutine hll(method,qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,idimmin,idimmax,&
    qtC,wCT,qt,wnew,wold,fC,dx1,x)

! method=='hll'  --> 2nd order HLL scheme.
! method=='hll1' --> 1st order HLL scheme.

include 'amrvacdef.f'

character(len=*), intent(in)                         :: method
double precision, intent(in)                         :: qdt, qtC, qt, dx1
integer, intent(in)                                  :: ixImin1,ixImax1,&
    ixOmin1,ixOmax1, idimmin,idimmax
double precision, dimension(ixImin1:ixImax1,1:ndim), intent(in) ::  x
double precision, dimension(ixImin1:ixImax1,1:ndim)             :: xi
double precision, dimension(ixImin1:ixImax1,1:nw)               :: wCT, wnew,&
    wold
double precision, dimension(ixImin1:ixImax1,1:nwflux,1:ndim)  :: fC

double precision, dimension(ixGlo1:ixGhi1,1:nw) :: wLC, wRC
double precision, dimension(ixGlo1:ixGhi1)      :: fLC, fRC, vLC, vRC
double precision, dimension(ixGlo1:ixGhi1)      :: cmaxC, cmaxRC, cmaxLC
double precision, dimension(ixGlo1:ixGhi1)      :: cminC, cminRC, cminLC
double precision, dimension(1:ndim)     :: dxinv, dxdim
integer, dimension(ixGlo1:ixGhi1)               :: patchf
integer :: idims, iw, ixmin1,ixmax1, hxOmin1,hxOmax1, ixCmin1,ixCmax1,&
    ixCRmin1,ixCRmax1, jxCmin1,jxCmax1, kxCmin1,kxCmax1, kxRmin1,kxRmax1
logical :: transport, new_cmax, CmaxMeanState, logiB
logical, dimension(ixGlo1:ixGhi1) :: patchw
!-----------------------------------------------------------------------------

CmaxMeanState = (typetvdlf=='cmaxmean')
logiB=(BnormLF.and.b0_>0)

if (idimmax>idimmin .and. typelimited=='original' .and. method&
   /='hll1')call mpistop("Error in hll: Unsplit dim. and original is limited")


! The flux calculation contracts by one in the idim direction it is applied.
! The limiter contracts the same directions by one more, so expand ixO by 2.
ixmin1=ixOmin1;ixmax1=ixOmax1;
do idims= idimmin,idimmax
   ixmin1=ixmin1-2*kr(idims,1);ixmax1=ixmax1+2*kr(idims,1);
end do
if (ixImin1>ixmin1.or.ixImax1<ixmax1) call mpistop&
   ("Error in hll : Nonconforming input limits")

if (method=='hll'.and.useprimitive) then  
   ! second order methods with primitive limiting: 
   ! this call ensures wCT is primitive with updated auxiliaries
   call primitive(ixImin1,ixImax1,ixImin1,ixImax1,wCT,x)
endif 

dxinv(1)=-qdt/dx1;
dxdim(1)=dx1;
do idims= idimmin,idimmax
   if (B0field) then
      select case (idims)
      case (1)
         myB0 => myB0_face1
      end select
   end if

   hxOmin1=ixOmin1-kr(idims,1);hxOmax1=ixOmax1-kr(idims,1);
   ! ixC is centered index in the idim direction from ixOmin-1/2 to ixOmax+1/2

   ixCmax1=ixOmax1; ixCmin1=hxOmin1;


   ! Calculate wRC=uR_{j+1/2} and wLC=uL_j+1/2 
   jxCmin1=ixCmin1+kr(idims,1);jxCmax1=ixCmax1+kr(idims,1);

   kxCmin1=ixImin1; kxCmax1=ixImax1-kr(idims,1); 
   kxRmin1=kxCmin1+kr(idims,1);kxRmax1=kxCmax1+kr(idims,1)
   

   ixCRmin1=ixCmin1;ixCRmax1=ixCmax1;

   wRC(kxCmin1:kxCmax1,1:nwflux)=wCT(kxRmin1:kxRmax1,1:nwflux)
   wLC(kxCmin1:kxCmax1,1:nwflux)=wCT(kxCmin1:kxCmax1,1:nwflux)

   ! Get interface positions:
   xi(kxCmin1:kxCmax1,1:ndim) = x(kxCmin1:kxCmax1,1:ndim)
   xi(kxCmin1:kxCmax1,idims) = half* ( x(kxRmin1:kxRmax1,idims)&
      +x(kxCmin1:kxCmax1,idims) )

   ! for hll (second order scheme): apply limiting
   if (method=='hll') then
      select case (typelimited)
      case ('previous')
         call upwindLR(ixImin1,ixImax1,ixCRmin1,ixCRmax1,ixCRmin1,ixCRmax1,&
            idims,wold,wCT,wLC,wRC,x,.true.,dxdim(idims))
      case ('predictor')
         call upwindLR(ixImin1,ixImax1,ixCRmin1,ixCRmax1,ixCRmin1,ixCRmax1,&
            idims,wCT,wCT,wLC,wRC,x,.false.,dxdim(idims))
      case ('original')
         call upwindLR(ixImin1,ixImax1,ixCRmin1,ixCRmax1,ixCRmin1,ixCRmax1,&
            idims,wnew,wCT,wLC,wRC,x,.true.,dxdim(idims))
      case default
         call mpistop("Error in hll: no such base for limiter")
      end select
   end if

   ! For the high order hll scheme the limiter is based on
   ! the maximum eigenvalue, it is calculated in advance.
   if (CmaxMeanState) then
         ! determine mean state and store in wLC
         wLC(ixCmin1:ixCmax1,1:nwflux)= half*(wLC(ixCmin1:ixCmax1,1:nwflux)&
            +wRC(ixCmin1:ixCmax1,1:nwflux))
         ! get auxilaries for mean state
         if (nwaux>0) then
            call getaux(.true.,wLC,xi,ixGlo1,ixGhi1,ixCmin1,ixCmax1,&
               'hll_cmaxmeanstate')
         end if
         new_cmax=.true.
         call getcmax(new_cmax,wLC,xi,ixGlo1,ixGhi1,ixCmin1,ixCmax1,idims,&
            cmaxC,cminC,.true.)

         ! We regain wLC for further use
         wLC(ixCmin1:ixCmax1,1:nwflux)=two*wLC(ixCmin1:ixCmax1,1:nwflux)&
            -wRC(ixCmin1:ixCmax1,1:nwflux)
         if (nwaux>0) then
            call getaux(.true.,wLC,xi,ixGlo1,ixGhi1,ixCmin1,ixCmax1,&
               'hll_wLC_B')
         endif
         if (nwaux>0.and.(.not.(useprimitive).or.method=='hll1')) then
            call getaux(.true.,wRC,xi,ixGlo1,ixGhi1,ixCmin1,ixCmax1,&
               'hll_wRC_B')
         end if
    else
         ! get auxilaries for L and R states
         if (nwaux>0.and.(.not.(useprimitive).or.method=='hll1')) then
         !!if (nwaux>0) then
            call getaux(.true.,wLC,xi,ixGlo1,ixGhi1,ixCmin1,ixCmax1,'hll_wLC')
            call getaux(.true.,wRC,xi,ixGlo1,ixGhi1,ixCmin1,ixCmax1,'hll_wRC')
         end if
         new_cmax=.true.
         call getcmax(new_cmax,wLC,xi,ixGlo1,ixGhi1,ixCmin1,ixCmax1,idims,&
            cmaxLC,cminLC,.true.)
         call getcmax(new_cmax,wRC,xi,ixGlo1,ixGhi1,ixCmin1,ixCmax1,idims,&
            cmaxRC,cminRC,.true.)
         ! now take the maximum of left and right states
         cmaxC(ixCmin1:ixCmax1)=max(cmaxRC(ixCmin1:ixCmax1),&
            cmaxLC(ixCmin1:ixCmax1))
         cminC(ixCmin1:ixCmax1)=min(cminRC(ixCmin1:ixCmax1),&
            cminLC(ixCmin1:ixCmax1))
   end if 

   patchf(ixCmin1:ixCmax1) =  1
   where(cminC(ixCmin1:ixCmax1) >= zero)
        patchf(ixCmin1:ixCmax1) = -2
   elsewhere(cmaxC(ixCmin1:ixCmax1) <= zero)
        patchf(ixCmin1:ixCmax1) =  2
   endwhere

   ! Calculate velocities for transport fluxes
   if(any(patchf(ixCmin1:ixCmax1)/= 2).or.(logiB)) call getv(wLC,xi,ixGlo1,&
      ixGhi1,ixCmin1,ixCmax1,idims,vLC)
   if(any(patchf(ixCmin1:ixCmax1)/=-2).or.(logiB)) call getv(wRC,xi,ixGlo1,&
      ixGhi1,ixCmin1,ixCmax1,idims,vRC)




   ! Calculate fLC=f(uL_j+1/2) and fRC=f(uR_j+1/2) for each iw
   do iw=1,nwflux
     if(any(patchf(ixCmin1:ixCmax1)/= 2).or.logiB.and.(iw==b0_+idims)) then 
        call getflux(wLC,xi,ixGlo1,ixGhi1,ixCmin1,ixCmax1,iw,idims,fLC,&
           transport)
        if (transport) fLC(ixCmin1:ixCmax1)=fLC(ixCmin1:ixCmax1)&
           +vLC(ixCmin1:ixCmax1)*wLC(ixCmin1:ixCmax1,iw)
     end if
     if(any(patchf(ixCmin1:ixCmax1)/=-2).or.logiB.and.(iw==b0_+idims)) then 
        call getflux(wRC,xi,ixGlo1,ixGhi1,ixCmin1,ixCmax1,iw,idims,fRC,&
           transport)
        if (transport) fRC(ixCmin1:ixCmax1)=fRC(ixCmin1:ixCmax1)&
           +vRC(ixCmin1:ixCmax1)*wRC(ixCmin1:ixCmax1,iw)
     end if


     if (logiB.and.(iw==b0_+idims)) then
         ! flat B norm using tvdlf
         fLC(ixCmin1:ixCmax1)= half*((fLC(ixCmin1:ixCmax1)+fRC&
            (ixCmin1:ixCmax1)) -tvdlfeps*max(cmaxC(ixCmin1:ixCmax1),&
            dabs(cminC(ixCmin1:ixCmax1)))*(wRC(ixCmin1:ixCmax1,iw)&
            -wLC(ixCmin1:ixCmax1,iw)))
     else
       where(patchf(ixCmin1:ixCmax1)==1)
         ! Add hll dissipation to the flux
         fLC(ixCmin1:ixCmax1)=(cmaxC(ixCmin1:ixCmax1)*fLC(ixCmin1:ixCmax1)&
            -cminC(ixCmin1:ixCmax1)*fRC(ixCmin1:ixCmax1) +tvdlfeps&
            *cminC(ixCmin1:ixCmax1)*cmaxC(ixCmin1:ixCmax1)*(wRC&
            (ixCmin1:ixCmax1,iw)-wLC(ixCmin1:ixCmax1,iw)))/(cmaxC&
            (ixCmin1:ixCmax1)-cminC(ixCmin1:ixCmax1))
       elsewhere(patchf(ixCmin1:ixCmax1)== 2)
         fLC(ixCmin1:ixCmax1)=fRC(ixCmin1:ixCmax1)
       elsewhere(patchf(ixCmin1:ixCmax1)==-2)
         fLC(ixCmin1:ixCmax1)=fLC(ixCmin1:ixCmax1)
       endwhere
     endif

     if (slab) then
         fC(ixCmin1:ixCmax1,iw,idims)=fLC(ixCmin1:ixCmax1)
     else
         select case (idims)
         case (1)
            fC(ixCmin1:ixCmax1,iw,1)=mygeo%surfaceC1(ixCmin1:ixCmax1)&
               *fLC(ixCmin1:ixCmax1)
         end select
     end if

   end do ! Next iw
end do ! Next idims



do idims= idimmin,idimmax
   hxOmin1=ixOmin1-kr(idims,1);hxOmax1=ixOmax1-kr(idims,1);
   do iw=1,nwflux

      ! Multiply the fluxes by -dt/dx since Flux fixing expects this
      if (slab) then
         fC(ixImin1:ixImax1,iw,idims)=dxinv(idims)*fC(ixImin1:ixImax1,iw,&
            idims)
         wnew(ixOmin1:ixOmax1,iw)=wnew(ixOmin1:ixOmax1,iw) &
            + (fC(ixOmin1:ixOmax1,iw,idims)-fC(hxOmin1:hxOmax1,iw,idims))
      else
         select case (idims)
         case (1)
            fC(ixImin1:ixImax1,iw,1)=-qdt*fC(ixImin1:ixImax1,iw,idims)
            wnew(ixOmin1:ixOmax1,iw)=wnew(ixOmin1:ixOmax1,iw) &
              + (fC(ixOmin1:ixOmax1,iw,1)-fC(hxOmin1:hxOmax1,iw,1))&
                 /mygeo%dvolume(ixOmin1:ixOmax1)
         end select
      end if

   end do ! Next iw

end do ! Next idims

if (method=='hll'.and.useprimitive) then  
   patchw(ixImin1:ixImax1)=.false.
   call conserve(ixImin1,ixImax1,ixImin1,ixImax1,wCT,x,patchw)
else
   if(nwaux>0) call getaux(.true.,wCT,x,ixImin1,ixImax1,ixImin1,ixImax1,&
      'hll_wCT')
endif

if (.not.slab.and.idimmin==1) call addgeometry(qdt,ixImin1,ixImax1,ixOmin1,&
   ixOmax1,wCT,wnew,x)
call addsource2(qdt*dble(idimmax-idimmin+1)/dble(ndim), ixImin1,ixImax1,&
   ixOmin1,ixOmax1,1,nw,qtC,wCT,qt,wnew,x,.false.)

end subroutine hll
!=============================================================================
subroutine hllc(method,qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,idimmin,idimmax,&
    qtC,wCT,qt,wnew,wold,fC,dx1,x)

! method=='hllc'  --> 2nd order HLLC scheme.
! method=='hllc1' --> 1st order HLLC scheme.
! method=='hllcd' --> 2nd order HLLC+tvdlf scheme.
! method=='hllcd1'--> 1st order HLLC+tvdlf scheme.

include 'amrvacdef.f'

character(len=*), intent(in)                         :: method
double precision, intent(in)                         :: qdt, qtC, qt, dx1
integer, intent(in)                                  :: ixImin1,ixImax1,&
    ixOmin1,ixOmax1, idimmin,idimmax
double precision, dimension(ixImin1:ixImax1,1:ndim), intent(in) ::  x
double precision, dimension(ixImin1:ixImax1,1:ndim)             ::  xi
double precision, dimension(ixImin1:ixImax1,1:nw)               :: wCT, wnew,&
    wold
double precision, dimension(ixImin1:ixImax1,1:nwflux,1:ndim)  :: fC

double precision, dimension(ixGlo1:ixGhi1,1:nw)            :: wLC, wRC
double precision, dimension(ixGlo1:ixGhi1)                 :: vLC, vRC
double precision, dimension(ixGlo1:ixGhi1)                 :: cmaxC,cminC

double precision, dimension(1:ndim)                :: dxinv, dxdim

integer, dimension(ixGlo1:ixGhi1)                          :: patchf
integer :: idims, iw, ixmin1,ixmax1, hxOmin1,hxOmax1, ixCmin1,ixCmax1,&
    ixCRmin1,ixCRmax1, jxCmin1,jxCmax1, kxCmin1,kxCmax1, kxRmin1,kxRmax1
logical :: transport, new_cmax, CmaxMeanState, logiB, firstordermethod
logical, dimension(ixGlo1:ixGhi1) :: patchw

!=== specific to HLLC and HLLCD ===!
double precision, dimension(ixGlo1:ixGhi1,1:nwflux)     :: fLC, fRC
double precision, dimension(ixGlo1:ixGhi1,1:nwflux)     :: whll, Fhll, fCD
double precision, dimension(ixGlo1:ixGhi1)              :: lambdaCD 
!-----------------------------------------------------------------------------

CmaxMeanState = (typetvdlf=='cmaxmean')
logiB=(BnormLF.and.b0_>0)

if (idimmax>idimmin .and. typelimited=='original' .and. method&
   /='hllc1' .and. method/='hllcd1')call mpistop&
   ("Error in hllc: Unsplit dim. and original is limited")



! The flux calculation contracts by one in the idim direction it is applied.
! The limiter contracts the same directions by one more, so expand ixO by 2.
ixmin1=ixOmin1;ixmax1=ixOmax1;
do idims= idimmin,idimmax
   ixmin1=ixmin1-2*kr(idims,1);ixmax1=ixmax1+2*kr(idims,1);
end do
if (ixImin1>ixmin1.or.ixImax1<ixmax1) call mpistop&
   ("Error in hllc : Nonconforming input limits")

if ((method=='hllc'.or.method=='hllcd').and.useprimitive) then
   call primitive(ixImin1,ixImax1,ixImin1,ixImax1,wCT,x)
endif 
firstordermethod=(method=='hllc1'.or.method=='hllcd1')

dxinv(1)=-qdt/dx1;
dxdim(1)=dx1;
do idims= idimmin,idimmax
   if (B0field) then
      select case (idims)
      case (1)
         myB0 => myB0_face1
      end select
   end if

   hxOmin1=ixOmin1-kr(idims,1);hxOmax1=ixOmax1-kr(idims,1);
   ! ixC is centered index in the idim direction from ixOmin-1/2 to ixOmax+1/2

   ixCmax1=ixOmax1; ixCmin1=hxOmin1;


   ! Calculate wRC=uR_{j+1/2} and wLC=uL_j+1/2 
   jxCmin1=ixCmin1+kr(idims,1);jxCmax1=ixCmax1+kr(idims,1);

   ! enlarged for ppm purposes
   kxCmin1=ixImin1; kxCmax1=ixImax1-kr(idims,1); 
   kxRmin1=kxCmin1+kr(idims,1);kxRmax1=kxCmax1+kr(idims,1)
   

   ixCRmin1=ixCmin1;ixCRmax1=ixCmax1;

   wRC(kxCmin1:kxCmax1,1:nwflux)=wCT(kxRmin1:kxRmax1,1:nwflux)
   wLC(kxCmin1:kxCmax1,1:nwflux)=wCT(kxCmin1:kxCmax1,1:nwflux)

   ! Get interface positions:
   xi(kxCmin1:kxCmax1,1:ndim) = x(kxCmin1:kxCmax1,1:ndim)
   xi(kxCmin1:kxCmax1,idims) = half* ( x(kxRmin1:kxRmax1,idims)&
      +x(kxCmin1:kxCmax1,idims) )

   ! for hllc and hllcd (second order schemes): apply limiting
   if (method=='hllc'.or.method=='hllcd') then
      select case (typelimited)
      case ('previous')
         call upwindLR(ixImin1,ixImax1,ixCRmin1,ixCRmax1,ixCRmin1,ixCRmax1,&
            idims,wold,wCT,wLC,wRC,x,.true.,dxdim(idims))
      case ('predictor')
         call upwindLR(ixImin1,ixImax1,ixCRmin1,ixCRmax1,ixCRmin1,ixCRmax1,&
            idims,wCT,wCT,wLC,wRC,x,.false.,dxdim(idims))
      case ('original')
         call upwindLR(ixImin1,ixImax1,ixCRmin1,ixCRmax1,ixCRmin1,ixCRmax1,&
            idims,wnew,wCT,wLC,wRC,x,.true.,dxdim(idims))
      case default
         call mpistop("Error in hllc: no such base for limiter")
      end select
   end if

   ! For the high order hllc scheme the limiter is based on
   ! the maximum eigenvalue, it is calculated in advance.
   if (CmaxMeanState) then
         ! determine mean state and store in wLC
         wLC(ixCmin1:ixCmax1,1:nwflux)= half*(wLC(ixCmin1:ixCmax1,1:nwflux)&
            +wRC(ixCmin1:ixCmax1,1:nwflux))
         ! get auxilaries for mean state
         if (nwaux>0) then
            call getaux(.true.,wLC,xi,ixGlo1,ixGhi1,ixCmin1,ixCmax1,&
               'hllc_cmaxmeanstate')
         end if
         new_cmax=.true.
         call getcmax(new_cmax,wLC,xi,ixGlo1,ixGhi1,ixCmin1,ixCmax1,idims,&
            cmaxC,cminC,.true.)

         ! We regain wLC for further use
         wLC(ixCmin1:ixCmax1,1:nwflux)=two*wLC(ixCmin1:ixCmax1,1:nwflux)&
            -wRC(ixCmin1:ixCmax1,1:nwflux)
         if (nwaux>0) then
            call getaux(.true.,wLC,xi,ixGlo1,ixGhi1,ixCmin1,ixCmax1,&
               'hllc_wLC_B')
         end if
         if (nwaux>0.and.(.not.(useprimitive).or.firstordermethod)) then
            call getaux(.true.,wRC,xi,ixGlo1,ixGhi1,ixCmin1,ixCmax1,&
               'hllc_wRC_B')
         end if
   else
         ! get auxilaries for L and R states
         if (nwaux>0.and.(.not.(useprimitive).or.firstordermethod)) then
         !!if (nwaux>0) then
            call getaux(.true.,wLC,xi,ixGlo1,ixGhi1,ixCmin1,ixCmax1,&
               'hllc_wLC')
            call getaux(.true.,wRC,xi,ixGlo1,ixGhi1,ixCmin1,ixCmax1,&
               'hllc_wRC')
         end if
         new_cmax=.true.
         ! to save memory, use cmaxC and lambdaCD for cmacRC and cmaxLC respectively
         ! to save memory, use vLC   and vRC      for cminRC and cminLC respectively
         call getcmax(new_cmax,wLC,xi,ixGlo1,ixGhi1,ixCmin1,ixCmax1,idims,&
            cmaxC,vLC,.true.)
         call getcmax(new_cmax,wRC,xi,ixGlo1,ixGhi1,ixCmin1,ixCmax1,idims,&
            lambdaCD,vRC,.true.)
         ! now take the maximum of left and right states
         cmaxC(ixCmin1:ixCmax1)=max(lambdaCD(ixCmin1:ixCmax1),&
            cmaxC(ixCmin1:ixCmax1))
         cminC(ixCmin1:ixCmax1)=min(vRC(ixCmin1:ixCmax1),vLC(ixCmin1:ixCmax1))
   end if

   patchf(ixCmin1:ixCmax1) =  1
   where(cminC(ixCmin1:ixCmax1) >= zero)
        patchf(ixCmin1:ixCmax1) = -2
   elsewhere(cmaxC(ixCmin1:ixCmax1) <= zero)
        patchf(ixCmin1:ixCmax1) =  2
   endwhere

   ! Calculate velocities for transport fluxes
   if(any(patchf(ixCmin1:ixCmax1)/= 2).or.(logiB)) call getv(wLC,xi,ixGlo1,&
      ixGhi1,ixCmin1,ixCmax1,idims,vLC)
   if(any(patchf(ixCmin1:ixCmax1)/=-2).or.(logiB)) call getv(wRC,xi,ixGlo1,&
      ixGhi1,ixCmin1,ixCmax1,idims,vRC)



   ! Calculate fLC=f(uL_j+1/2) and fRC=f(uR_j+1/2) for each iw
   do iw=1,nwflux
     if(any(patchf(ixCmin1:ixCmax1)/= 2).or.(logiB.and.(iw==b0_+idims))) then
        call getfluxforhllc(wLC,xi,ixGlo1,ixGhi1,ixCmin1,ixCmax1,iw,idims,fLC,&
           transport)
        if (transport)  fLC(ixCmin1:ixCmax1,iw)=fLC(ixCmin1:ixCmax1,iw)&
           +vLC(ixCmin1:ixCmax1)*wLC(ixCmin1:ixCmax1,iw)
     end if
     if(any(patchf(ixCmin1:ixCmax1)/=-2).or.(logiB.and.(iw==b0_+idims))) then
        call getfluxforhllc(wRC,xi,ixGlo1,ixGhi1,ixCmin1,ixCmax1,iw,idims,fRC,&
           transport)
        if (transport)   fRC(ixCmin1:ixCmax1,iw)=fRC(ixCmin1:ixCmax1,iw)&
           +vRC(ixCmin1:ixCmax1)*wRC(ixCmin1:ixCmax1,iw)
     end if
   end do

   ! Use more diffusive scheme, is actually TVDLF and selected by patchf=4 
   if(method=='hllcd' .or. method=='hllcd1') call diffuse_hllcd(ixGlo1,ixGhi1,&
      ixCmin1,ixCmax1,idims,wLC,wRC,fLC,fRC,patchf)

   !---- calculate speed lambda at CD ----!
   if(any(patchf(ixCmin1:ixCmax1)==1)) call getlCD(wLC,wRC,fLC,fRC,cminC,&
      cmaxC,idims,ixGlo1,ixGhi1,ixCmin1,ixCmax1, whll,Fhll,lambdaCD,patchf)

   ! now patchf may be -1 or 1 due to getlCD 
   if(any(abs(patchf(ixCmin1:ixCmax1))== 1))then
      !======== flux at intermediate state ========!
      call getwCD(wLC,wRC,whll,vLC,vRC,fRC,fLC,Fhll,patchf,lambdaCD,cminC,&
         cmaxC,ixGlo1,ixGhi1,ixCmin1,ixCmax1,idims,fCD)
   endif ! Calculate the CD flux

   do iw=1,nwflux
     if (logiB.and.(iw==b0_+idims)) then
       fLC(ixCmin1:ixCmax1,iw) = half*((fLC(ixCmin1:ixCmax1,iw)&
          +fRC(ixCmin1:ixCmax1,iw)) -tvdlfeps*max(cmaxC(ixCmin1:ixCmax1),&
          dabs(cminC(ixCmin1:ixCmax1)))*(wRC(ixCmin1:ixCmax1,iw)&
          -wLC(ixCmin1:ixCmax1,iw)))
     else
       where(patchf(ixCmin1:ixCmax1)==-2)
        fLC(ixCmin1:ixCmax1,iw)=fLC(ixCmin1:ixCmax1,iw)
       elsewhere(abs(patchf(ixCmin1:ixCmax1))==1)
        fLC(ixCmin1:ixCmax1,iw)=fCD(ixCmin1:ixCmax1,iw)
       elsewhere(patchf(ixCmin1:ixCmax1)==2)
        fLC(ixCmin1:ixCmax1,iw)=fRC(ixCmin1:ixCmax1,iw)
       elsewhere(patchf(ixCmin1:ixCmax1)==3)
        ! fallback option, reducing to HLL flux
        fLC(ixCmin1:ixCmax1,iw)=Fhll(ixCmin1:ixCmax1,iw)
       elsewhere(patchf(ixCmin1:ixCmax1)==4)
        ! fallback option, reducing to TVDLF flux
        fLC(ixCmin1:ixCmax1,iw) = half*((fLC(ixCmin1:ixCmax1,iw)&
           +fRC(ixCmin1:ixCmax1,iw)) -tvdlfeps*max(cmaxC(ixCmin1:ixCmax1),&
           dabs(cminC(ixCmin1:ixCmax1)))*(wRC(ixCmin1:ixCmax1,iw)&
           -wLC(ixCmin1:ixCmax1,iw)))
       endwhere
     endif

     if (slab) then
         fC(ixCmin1:ixCmax1,iw,idims)=fLC(ixCmin1:ixCmax1,iw)
     else
         select case (idims)
         case (1)
            fC(ixCmin1:ixCmax1,iw,1)=mygeo%surfaceC1(ixCmin1:ixCmax1)&
               *fLC(ixCmin1:ixCmax1,iw)
         end select
      end if

   end do ! Next iw

end do ! Next idims





!Now update the state
do idims= idimmin,idimmax
   hxOmin1=ixOmin1-kr(idims,1);hxOmax1=ixOmax1-kr(idims,1);
   do iw=1,nwflux

      ! Multiply the fluxes by -dt/dx since Flux fixing expects this
      if (slab) then
         fC(ixImin1:ixImax1,iw,idims)=dxinv(idims)*fC(ixImin1:ixImax1,iw,&
            idims)
         wnew(ixOmin1:ixOmax1,iw)=wnew(ixOmin1:ixOmax1,iw) &
            + (fC(ixOmin1:ixOmax1,iw,idims)-fC(hxOmin1:hxOmax1,iw,idims))
      else
         select case (idims)
         case (1)
            fC(ixImin1:ixImax1,iw,1)=-qdt*fC(ixImin1:ixImax1,iw,idims)
            wnew(ixOmin1:ixOmax1,iw)=wnew(ixOmin1:ixOmax1,iw) &
              + (fC(ixOmin1:ixOmax1,iw,1)-fC(hxOmin1:hxOmax1,iw,1))&
                 /mygeo%dvolume(ixOmin1:ixOmax1)
         end select
      end if

   end do ! Next iw

end do ! Next idims

if ((method=='hllc'.or.method=='hllcd').and.useprimitive) then  
   patchw(ixImin1:ixImax1)=.false.
   call conserve(ixImin1,ixImax1,ixImin1,ixImax1,wCT,x,patchw)
else
   if(nwaux>0) call getaux(.true.,wCT,x,ixImin1,ixImax1,ixImin1,ixImax1,&
      'hllc_wCT')
endif

if (.not.slab.and.idimmin==1) call addgeometry(qdt,ixImin1,ixImax1,ixOmin1,&
   ixOmax1,wCT,wnew,x)
call addsource2(qdt*dble(idimmax-idimmin+1)/dble(ndim), ixImin1,ixImax1,&
   ixOmin1,ixOmax1,1,nw,qtC,wCT,qt,wnew,x,.false.)

end subroutine hllc
!=============================================================================
subroutine tvdmusclf(method,qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,idimmin,&
   idimmax, qtC,wCT,qt,wnew,wold,fC,dx1,x)

! method=='tvdmu'  --> 2nd order (may be 3rd order in 1D) TVD-MUSCL scheme.
! method=='tvdmu1' --> 1st order TVD-MUSCL scheme (upwind per charact. var.)
! FCT not implemented here.
include 'amrvacdef.f'

character(len=*), intent(in)                         :: method
double precision, intent(in)                         :: qdt, qtC, qt, dx1
integer, intent(in)                                  :: ixImin1,ixImax1,&
    ixOmin1,ixOmax1, idimmin,idimmax
double precision, dimension(ixImin1:ixImax1,1:ndim), intent(in) ::  x
double precision, dimension(ixImin1:ixImax1,1:ndim)             ::  xi
double precision, dimension(ixImin1:ixImax1,1:nw)               :: wCT, wnew,&
    wold
double precision, dimension(ixImin1:ixImax1,1:nwflux,1:ndim)        :: fC

double precision, dimension(ixGlo1:ixGhi1,1:nw) :: wLC, wRC
double precision, dimension(ixGlo1:ixGhi1)      :: fLC, fRC, vLC, vRC
double precision, dimension(ixGlo1:ixGhi1)      :: cmaxC, cmaxRC, cmaxLC,ie,&
   ieL,ieR,ieCT
double precision :: dxinv(1:ndim),dxdim(1:ndim)
integer :: idims, iw, ixmin1,ixmax1, hxOmin1,hxOmax1, ixCmin1,ixCmax1,&
    ixCRmin1,ixCRmax1, kxCmin1,kxCmax1, kxRmin1,kxRmax1, ixtestmin1,ixtestmax1
logical :: transport, new_cmax, CmaxMeanState
logical, dimension(ixGlo1:ixGhi1) :: patchw
logical, save :: evolvepressure=.false.
!-----------------------------------------------------------------------------

CmaxMeanState = (typetvdlf=='cmaxmean')

oktest=.false.
if(oktest.and.mype==0) then
   print *,'======tvdmusl: qdt, qtC, qt:',qdt,qtC,qt
   ixtestmin1=ixImin1;ixtestmax1=ixImax1;
   ixtestmin1=1;ixtestmax1=2*dixB;
   print *,'reporting in ranges:',ixtestmin1,ixtestmax1
endif

!!call mpistop("tijdelijke stop")

if (idimmax>idimmin .and. typelimited=='original' .and. method&
   /='tvdmu1')call mpistop&
   ("Error in TVDMUSCLF: Unsplit dim. and original is limited")

! The flux calculation contracts by one in the idim direction it is applied.
! The limiter contracts the same directions by one more, so expand ixO by 2.
ixmin1=ixOmin1;ixmax1=ixOmax1;
do idims= idimmin,idimmax
   ixmin1=ixmin1-2*kr(idims,1);ixmax1=ixmax1+2*kr(idims,1);
end do
if (ixImin1>ixmin1.or.ixImax1<ixmax1) call mpistop&
   ("Error in TVDMUSCLF: Nonconforming input limits")





if ((method=='tvdmu').and.useprimitive) then  
   ! second order methods with primitive limiting: 
   ! this call ensures wCT is primitive with updated auxiliaries
   call primitive(ixImin1,ixImax1,ixImin1,ixImax1,wCT,x)
endif 


dxinv(1)=-qdt/dx1;
dxdim(1)=dx1;
do idims= idimmin,idimmax
   if (B0field) then
      select case (idims)
      case (1)
         myB0 => myB0_face1
      end select
   end if

   hxOmin1=ixOmin1-kr(idims,1);hxOmax1=ixOmax1-kr(idims,1);
   ! ixC is centered index in the idim direction from ixOmin-1/2 to ixOmax+1/2
   ixCmax1=ixOmax1; ixCmin1=hxOmin1;

   kxCmin1=ixImin1; kxCmax1=ixImax1-kr(idims,1);
   kxRmin1=kxCmin1+kr(idims,1);kxRmax1=kxCmax1+kr(idims,1);

   ixCRmin1=ixCmin1;ixCRmax1=ixCmax1;

   wRC(kxCmin1:kxCmax1,1:nwflux)=wCT(kxRmin1:kxRmax1,1:nwflux)
   wLC(kxCmin1:kxCmax1,1:nwflux)=wCT(kxCmin1:kxCmax1,1:nwflux)

   ! Get interface positions:
   xi(kxCmin1:kxCmax1,1:ndim) = x(kxCmin1:kxCmax1,1:ndim)
   xi(kxCmin1:kxCmax1,idims) = half* ( x(kxRmin1:kxRmax1,idims)&
      +x(kxCmin1:kxCmax1,idims) )

   ! for tvdlf and tvdmu (second order schemes): apply limiting
   if (method=='tvdmu') then
      select case (typelimited)
      case ('previous')
         call upwindLR(ixImin1,ixImax1,ixCRmin1,ixCRmax1,ixCRmin1,ixCRmax1,&
            idims,wold,wCT,wLC,wRC,x,.true.,dxdim(idims))
      case ('predictor')
         call upwindLR(ixImin1,ixImax1,ixCRmin1,ixCRmax1,ixCRmin1,ixCRmax1,&
            idims,wCT,wCT,wLC,wRC,x,.false.,dxdim(idims))
      case ('original')
         call upwindLR(ixImin1,ixImax1,ixCRmin1,ixCRmax1,ixCRmin1,ixCRmax1,&
            idims,wnew,wCT,wLC,wRC,x,.true.,dxdim(idims))
      case default
         call mpistop("Error in TVDMUSCLF: no such base for limiter")
      end select
   end if

   ! handle all other methods than tvdlf, namely tvdmu and tvdmu1 here
   if (nwaux>0.and.(.not.(useprimitive).or.method=='tvdmu1')) then
      !!if (nwaux>0) then
      call getaux(.true.,wLC,xi,ixGlo1,ixGhi1,ixCmin1,ixCmax1,'tvdlf_wLC_B')
      call getaux(.true.,wRC,xi,ixGlo1,ixGhi1,ixCmin1,ixCmax1,'tvdlf_wRC_B')
   end if
   

   ! Calculate velocities for transport fluxes
   call getv(wLC,xi,ixGlo1,ixGhi1,ixCmin1,ixCmax1,idims,vLC)
   call getv(wRC,xi,ixGlo1,ixGhi1,ixCmin1,ixCmax1,idims,vRC)




   ! Calculate fLC=f(uL_j+1/2) and fRC=f(uR_j+1/2) for each iw
   do iw=1,nwflux
      call getflux(wLC,xi,ixGlo1,ixGhi1,ixCmin1,ixCmax1,iw,idims,fLC,&
         transport)
      call getflux(wRC,xi,ixGlo1,ixGhi1,ixCmin1,ixCmax1,iw,idims,fRC,&
         transport)
      if (transport) then
         fLC(ixCmin1:ixCmax1)=fLC(ixCmin1:ixCmax1)+vLC(ixCmin1:ixCmax1)&
            *wLC(ixCmin1:ixCmax1,iw)
         fRC(ixCmin1:ixCmax1)=fRC(ixCmin1:ixCmax1)+vRC(ixCmin1:ixCmax1)&
            *wRC(ixCmin1:ixCmax1,iw)
      end if
      ! To save memory we use fLC to store (F_L+F_R)/2=half*(fLC+fRC)
      fLC(ixCmin1:ixCmax1)=half*(fLC(ixCmin1:ixCmax1)+fRC(ixCmin1:ixCmax1))

      if (slab) then
         fC(ixCmin1:ixCmax1,iw,idims)=dxinv(idims)*fLC(ixCmin1:ixCmax1)
         wnew(ixOmin1:ixOmax1,iw)=wnew(ixOmin1:ixOmax1,iw)+ &
            (fC(ixOmin1:ixOmax1,iw,idims)-fC(hxOmin1:hxOmax1,iw,idims))
      else
         select case (idims)
         case (1)
            fC(ixCmin1:ixCmax1,iw,1)=-qdt*mygeo%surfaceC1(ixCmin1:ixCmax1)&
               *fLC(ixCmin1:ixCmax1)
            wnew(ixOmin1:ixOmax1,iw)=wnew(ixOmin1:ixOmax1,iw)+ &
              (fC(ixOmin1:ixOmax1,iw,1)-fC(hxOmin1:hxOmax1,iw,1))&
                 /mygeo%dvolume(ixOmin1:ixOmax1)
         end select
      end if

   end do ! Next iw

   ! For the MUSCL scheme apply the characteristic based limiter
   if (method=='tvdmu'.or.method=='tvdmu1') call tvdlimit2(method,qdt,ixImin1,&
      ixImax1,ixCmin1,ixCmax1,ixOmin1,ixOmax1,idims,wLC,wRC,wnew,x,fC,dx1)

end do ! Next idims



if ((method=='tvdmu').and.useprimitive) then  
    patchw(ixImin1:ixImax1)=.false.
    call conserve(ixImin1,ixImax1,ixImin1,ixImax1,wCT,x,patchw)
else
   if(nwaux>0) call getaux(.true.,wCT,x,ixImin1,ixImax1,ixImin1,ixImax1,&
      'tvdlf_wCT')
endif

if (.not.slab.and.idimmin==1) call addgeometry(qdt,ixImin1,ixImax1,ixOmin1,&
   ixOmax1,wCT,wnew,x)
call addsource2(qdt*dble(idimmax-idimmin+1)/dble(ndim),ixImin1,ixImax1,&
   ixOmin1,ixOmax1,1,nw,qtC,wCT,qt,wnew,x,.false.)

end subroutine tvdmusclf
!=============================================================================

