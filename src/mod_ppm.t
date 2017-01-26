module mod_ppm

  implicit none
  private

  public :: PPMlimiter
  public :: PPMlimitervar

contains

  subroutine PPMlimitervar(ixI^L,ix^L,idims,q,qCT,qLC,qRC)

    ! references:
    ! Mignone et al 2005, ApJS 160, 199,
    ! Miller and Colella 2002, JCP 183, 26
    ! Fryxell et al. 2000 ApJ, 131, 273 (Flash)
    ! baciotti Phd (http://www.aei.mpg.de/~baiotti/Baiotti_PhD.pdf)
    ! version : april 2009
    ! author: zakaria.meliani@wis.kuleuven.be

    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ix^L, idims
    double precision, intent(in)    :: q(ixI^S),qCT(ixI^S)

    double precision, intent(inout) :: qRC(ixG^T),qLC(ixG^T)

    double precision,dimension(ixG^T)  :: dqC,d2qC,ldq
    double precision,dimension(ixG^T)  :: qMin,qMax,tmp

    integer   :: lxC^L,lxR^L
    integer   :: ixLL^L,ixL^L,ixO^L,ixR^L,ixRR^L
    integer   :: hxL^L,hxC^L,hxR^L
    integer   :: kxLL^L,kxL^L,kxC^L,kxR^L,kxRR^L
    !--------------------------------------------------------------------------
    ixOmin^D=ixmin^D-kr(idims,^D);ixOmax^D=ixmax^D+kr(idims,^D);!ixO[ixMmin1-1,ixMmax1+1]
    ixL^L=ixO^L-kr(idims,^D);                                   !ixL[ixMmin1-2,ixMmax1]
    ixLL^L=ixL^L-kr(idims,^D);                                  !ixLL[ixMmin1-3,ixMmax1-1]
    ixR^L=ixO^L+kr(idims,^D);                                   !ixR=[iMmin1,ixMmax+2]
    ixRR^L=ixR^L+kr(idims,^D);                                  !ixRR=[iMmin1+1,ixMmax+3]

    hxCmin^D=ixOmin^D;hxCmax^D=ixmax^D; ! hxC = [ixMmin-1,ixMmax]
    hxL^L=hxC^L-kr(idims,^D);           ! hxL = [ixMmin-2,ixMmax-1]
    hxR^L=hxC^L+kr(idims,^D);           ! hxR = [ixMmin,ixMmax+1]

    kxCmin^D=ixLLmin^D; kxCmax^D=ixRmax^D; ! kxC=[iMmin1-3,ixMmax1+2]
    kxL^L=kxC^L-kr(idims,^D);              ! kxL=[iMmin1-4,ixMmax1+1]
    kxR^L=kxC^L+kr(idims,^D);              ! kxR=[iMmin1-2,ixMmax1+3]

    lxCmin^D=ixLLmin^D-kr(idims,^D);lxCmax^D=ixRRmax^D;! ixC=[iMmin1-4,ixMmax1+3]
    lxR^L=lxC^L+kr(idims,^D);                          ! lxR=[iMmin1-3,ixMmax1+4]


    dqC(lxC^S)=q(lxR^S)-q(lxC^S)
    ! Eq. 64,  Miller and Colella 2002, JCP 183, 26
    d2qC(kxC^S)=half*(q(kxR^S)-q(kxL^S))
    where(dqC(kxC^S)*dqC(kxL^S)>zero)
       ! Store the sign of d2qC in qMin
       qMin(kxC^S)= sign(one,d2qC(kxC^S))
       ! Eq. 65,  Miller and Colella 2002, JCP 183, 26
       ldq(kxC^S)= qMin(kxC^S)*min(dabs(d2qC(kxC^S)),&
            2.0d0*dabs(dqC(kxL^S)),&
            2.0d0*dabs(dqC(kxC^S)))
    elsewhere
       ldq(kxC^S)=zero
    end where

    ! Eq. 66, Miller and Colella 2002, JCP 183, 26
    qLC(ixO^S)=qLC(ixO^S)+half*dqC(ixO^S)&
         +(ldq(ixO^S)-ldq(ixR^S))/6.0d0
    if (flatppm)then
       qRC(ixL^S)=qRC(ixLL^S)+(half*dqC(ixL^S)&
            +(ldq(ixL^S)-ldq(ixO^S))/6.0d0)
    else
       qRC(ixL^S)=qRC(ixL^S) -(half*dqC(ixL^S)&
            +(ldq(ixL^S)-ldq(ixO^S))/6.0d0)
    endif

    ! make sure that min wCT(i)<wLC(i)<wCT(i+1) same for wRC(i)
    call extremaq(ixI^L,kxC^L,qCT,1,qMax,qMin)

    ! Eq. B8, page 217, Mignone et al 2005, ApJS
    qRC(ixL^S)=max(qMin(ixO^S),min(qMax(ixO^S),qRC(ixL^S)))
    qLC(ixO^S)=max(qMin(ixO^S),min(qMax(ixO^S),qLC(ixO^S)))

    ! Eq. B9, page 217, Mignone et al 2005, ApJS
    where((qRC(ixL^S)-qCT(ixO^S))*(qCT(ixO^S)-qLC(ixO^S))<=zero)
       qRC(ixL^S)=qCT(ixO^S)
       qLC(ixO^S)=qCT(ixO^S)
    end where

    qMax(ixO^S)=(qLC(ixO^S)-qRC(ixL^S))&
         *(qCT(ixO^S)-(qLC(ixO^S)+qRC(ixL^S))/2.0d0)
    qMin(ixO^S)=(qLC(ixO^S)-qRC(ixL^S))**2.0d0/6.0d0
    tmp(ixL^S)=qRC(ixL^S)

    ! Eq. B10, page 218, Mignone et al 2005, ApJS
    where(qMax(hxR^S)>qMin(hxR^S))
       qRC(hxC^S)= 3.0d0*qCT(hxR^S)-2.0d0*qLC(hxR^S)
    end where
    ! Eq. B11, page 218, Mignone et al 2005, ApJS
    where(qMax(hxC^S)<-qMin(hxC^S))
       qLC(hxC^S)= 3.0d0*qCT(hxC^S)-2.0d0*tmp(hxL^S)
    end where

  end subroutine PPMlimitervar

  subroutine PPMlimiter(ixI^L,ix^L,idims,w,wCT,wLC,wRC)

    ! references:
    ! Mignone et al 2005, ApJS 160, 199, 
    ! Miller and Colella 2002, JCP 183, 26 
    ! Fryxell et al. 2000 ApJ, 131, 273 (Flash)
    ! baciotti Phd (http://www.aei.mpg.de/~baiotti/Baiotti_PhD.pdf)
    ! version : april 2009
    ! author: zakaria.meliani@wis.kuleuven.be

    use mod_global_parameters
    use mod_physics, only: phys_ppm_flatcd, phys_ppm_flatsh

    integer, intent(in)             :: ixI^L, ix^L, idims
    double precision, intent(in)    :: w(ixI^S,1:nw),wCT(ixI^S,1:nw)

    double precision, intent(inout) :: wRC(ixG^T,1:nw),wLC(ixG^T,1:nw) 

    double precision,dimension(ixG^T,1:nwflux)  :: dwC,d2wC,ldw
    double precision,dimension(ixG^T,1:nwflux)  :: wMin,wMax,tmp
    double precision,dimension(ixG^T) :: aa, ab, ac, dv
    double precision,dimension(ixG^T,1:ndim) ::  exi

    integer   :: lxC^L,lxR^L
    integer   :: ixLL^L,ixL^L,ixO^L,ixR^L,ixRR^L
    integer   :: hxL^L,hxC^L,hxR^L
    integer   :: kxLL^L,kxL^L,kxC^L,kxR^L,kxRR^L
    integer   :: iw, idimss

    double precision, parameter :: betamin=0.75d0, betamax=0.85d0,&
         Zmin=0.25d0, Zmax=0.75d0,&
         eta1=20.0d0,eta2=0.05d0,eps=0.01d0,kappa=0.1d0
    !--------------------------------------------------------------------------
    ixOmin^D=ixmin^D-kr(idims,^D);ixOmax^D=ixmax^D+kr(idims,^D);!ixO[ixMmin1-1,ixMmax1+1]
    ixL^L=ixO^L-kr(idims,^D);                                   !ixL[ixMmin1-2,ixMmax1]
    ixLL^L=ixL^L-kr(idims,^D);                                  !ixLL[ixMmin1-3,ixMmax1-1]
    ixR^L=ixO^L+kr(idims,^D);                                   !ixR=[iMmin1,ixMmax+2]
    ixRR^L=ixR^L+kr(idims,^D);                                  !ixRR=[iMmin1+1,ixMmax+3]

    hxCmin^D=ixOmin^D;hxCmax^D=ixmax^D; ! hxC = [ixMmin-1,ixMmax]
    hxL^L=hxC^L-kr(idims,^D);           ! hxL = [ixMmin-2,ixMmax-1]
    hxR^L=hxC^L+kr(idims,^D);           ! hxR = [ixMmin,ixMmax+1]

    kxCmin^D=ixLLmin^D; kxCmax^D=ixRmax^D; ! kxC=[iMmin1-3,ixMmax1+2]
    kxL^L=kxC^L-kr(idims,^D);              ! kxL=[iMmin1-4,ixMmax1+1]
    kxR^L=kxC^L+kr(idims,^D);              ! kxR=[iMmin1-2,ixMmax1+3]

    lxCmin^D=ixLLmin^D-kr(idims,^D);lxCmax^D=ixRRmax^D;! ixC=[iMmin1-4,ixMmax1+3]
    lxR^L=lxC^L+kr(idims,^D);                          ! lxR=[iMmin1-3,ixMmax1+4]

    dwC(lxC^S,1:nwflux)=w(lxR^S,1:nwflux)-w(lxC^S,1:nwflux)
    ! Eq. 64,  Miller and Colella 2002, JCP 183, 26 
    d2wC(kxC^S,1:nwflux)=half*(w(kxR^S,1:nwflux)-w(kxL^S,1:nwflux))
    where(dwC(kxC^S,1:nwflux)*dwC(kxL^S,1:nwflux)>zero)
       ! Store the sign of dwC in wMin
       wMin(kxC^S,1:nwflux)= sign(one,d2wC(kxC^S,1:nwflux))
       ! Eq. 65,  Miller and Colella 2002, JCP 183, 26 
       ldw(kxC^S,1:nwflux)= wMin(kxC^S,1:nwflux)*min(dabs(d2wC(kxC^S,1:nwflux)),&
            2.0d0*dabs(dwC(kxL^S,1:nwflux)),&
            2.0d0*dabs(dwC(kxC^S,1:nwflux)))
    elsewhere
       ldw(kxC^S,1:nwflux)=zero
    endwhere

    ! Eq. 66,  Miller and Colella 2002, JCP 183, 26 
    wLC(ixO^S,1:nwflux)=wLC(ixO^S,1:nwflux)+half*dwC(ixO^S,1:nwflux)&
         +(ldw(ixO^S,1:nwflux)-ldw(ixR^S,1:nwflux))/6.0d0
    if(flatppm)then 
       wRC(ixL^S,1:nwflux)=wRC(ixLL^S,1:nwflux)+half*dwC(ixL^S,1:nwflux)&
            +(ldw(ixL^S,1:nwflux)-ldw(ixO^S,1:nwflux))/6.0d0
    else
       wRC(ixL^S,1:nwflux)=wRC(ixL^S,1:nwflux)-(half*dwC(ixL^S,1:nwflux)&
            +(ldw(ixL^S,1:nwflux)-ldw(ixO^S,1:nwflux))/6.0d0)
    endif

    ! make sure that min wCT(i)<wLC(i)<wCT(i+1) same for wRC(i)
    call extremaw(ixI^L,kxC^L,wCT,1,wMax,wMin)

    ! Eq. B8, page 217, Mignone et al 2005, ApJS
    wRC(ixL^S,1:nwflux)=max(wMin(ixO^S,1:nwflux)&
         ,min(wMax(ixO^S,1:nwflux),wRC(ixL^S,1:nwflux))) 
    wLC(ixO^S,1:nwflux)=max(wMin(ixO^S,1:nwflux)&
         ,min(wMax(ixO^S,1:nwflux),wLC(ixO^S,1:nwflux)))


    ! Eq. B9, page 217, Mignone et al 2005, ApJS
    where((wRC(ixL^S,1:nwflux)-wCT(ixO^S,1:nwflux))&
         *(wCT(ixO^S,1:nwflux)-wLC(ixO^S,1:nwflux))<=zero)
       wRC(ixL^S,1:nwflux)=wCT(ixO^S,1:nwflux)
       wLC(ixO^S,1:nwflux)=wCT(ixO^S,1:nwflux)
    end where

    wMax(ixO^S,1:nwflux)=(wLC(ixO^S,1:nwflux)-wRC(ixL^S,1:nwflux))*&
         (wCT(ixO^S,1:nwflux)-&
         (wLC(ixO^S,1:nwflux)+wRC(ixL^S,1:nwflux))/2.0d0)
    wMin(ixO^S,1:nwflux)=(wLC(ixO^S,1:nwflux)-wRC(ixL^S,1:nwflux))**2.0d0/6.0d0
    tmp(ixL^S,1:nwflux)=wRC(ixL^S,1:nwflux)
    ! Eq. B10, page 218, Mignone et al 2005, ApJS
    where(wMax(hxR^S,1:nwflux)>wMin(hxR^S,1:nwflux))
       wRC(hxC^S,1:nwflux)= 3.0d0*wCT(hxR^S,1:nwflux)&
            -2.0d0*wLC(hxR^S,1:nwflux)
    endwhere
    ! Eq. B11, page 218, Mignone et al 2005, ApJS
    where(wMax(hxC^S,1:nwflux)<-wMin(hxC^S,1:nwflux))
       wLC(hxC^S,1:nwflux)= 3.0d0*wCT(hxC^S,1:nwflux)&
            -2.0d0*tmp(hxL^S,1:nwflux)
    endwhere

    ! flattening at the contact discontinuities
    if(flatcd)then
       call phys_ppm_flatcd(ixI^L,kxC^L,kxL^L,kxR^L,wCT,d2wC,aa,ab)
       if(any(kappa*aa(kxC^S)>=ab(kxC^S)))then
          do iw=1,nwflux
             where(kappa*aa(kxC^S)>=ab(kxC^S).and. dabs(dwC(kxC^S,iw))>smalldouble)
                wMax(kxC^S,iw) = wCT(kxR^S,iw)-2.0d0*wCT(kxC^S,iw)+wCT(kxL^S,iw)
             end where

             where(wMax(ixR^S,iw)*wMax(ixL^S,iw)<zero .and.&
                  dabs(wCT(ixR^S,iw)-wCT(ixL^S,iw))&
                  -eps*min(dabs(wCT(ixR^S,iw)),dabs(wCT(ixL^S,iw)))>zero &
                  .and. kappa*aa(ixO^S)>=ab(ixO^S)&
                  .and. dabs(dwC(ixO^S,iw))>smalldouble)

                ac(ixO^S)=(wCT(ixLL^S,iw)-wCT(ixRR^S,iw)+4.0d0*dwC(ixO^S,iw))&
                     /(12.0d0*dwC(ixO^S,iw))
                wMin(ixO^S,iw)=max(zero,min(eta1*(ac(ixO^S)-eta2),one))
             elsewhere
                wMin(ixO^S,iw)=zero
             end where

             where(wMin(hxC^S,iw)>zero)
                wLC(hxC^S,iw) = wLC(hxC^S,iw)*(one-wMin(hxC^S,iw))&
                     +(wCT(hxC^S,iw)+half*ldw(hxC^S,iw))*wMin(hxC^S,iw)
             end where
             where(wMin(hxR^S,iw)>zero)
                wRC(hxC^S,iw) = wRC(hxC^S,iw)*(one-wMin(hxR^S,iw))&
                     +(wCT(hxR^S,iw)-half*ldw(hxR^S,iw))*wMin(hxR^S,iw)
             end where
          end do
       endif
    endif

    ! flattening at the shocks
    if(flatsh)then
       ! following MILLER and COLELLA 2002 JCP 183, 26
       kxCmin^D=ixmin^D-2; kxCmax^D=ixmax^D+2; ! kxC=[ixMmin1-2,ixMmax1+2]
       do idimss=1,ndim
          kxL^L=kxC^L-kr(idimss,^D); ! kxL=[ixMmin1-3,ixMmax1+1]
          kxR^L=kxC^L+kr(idimss,^D); ! kxR=[ixMmin1-1,ixMmax1+3]
          kxLL^L=kxL^L-kr(idimss,^D);! kxLL=[ixMmin-4,ixMmax]
          kxRR^L=kxR^L+kr(idimss,^D);! kxRR=[ixMmin,ixMmax+4]

          call phys_ppm_flatsh(ixI^L,kxC^L,kxLL^L,kxL^L,kxR^L,kxRR^L,idimss,wCT,aa,ab,dv)

          ! eq. B17, page 218, Mignone et al 2005, ApJS (had(Xi1))
          ac(kxC^S) = max(zero,min(one,(betamax-aa(kxC^S))/(betamax-betamin)))
          ! eq. B18, page 218, Mignone et al 2005, ApJS (had(Xi1))
          ! recycling aa(ixL^S)
          where (dabs(dv(kxC^S))<smalldouble)
             aa(kxC^S) = max(ac(kxC^S), min(one,(Zmax-ab(kxC^S))/(Zmax-Zmin)))
          elsewhere
             aa(kxC^S) = one
          endwhere

          {^IFONED call extremaa(ixI^L,ixO^L,aa,1,ab)}
          {^NOONED call extremaa(ixI^L,ixO^L,aa,1,exi(ixG^T,idimss))}
       enddo
       {^NOONED ab(ixO^S)=min(^D&exi(ixO^S,^D))}
       ! recycling wMax
       do iw=1,nwflux
          where(dabs(ab(ixO^S)-one)>smalldouble)
             wMax(ixO^S,iw) = (one-ab(ixO^S))*wCT(ixO^S,iw)
          endwhere

          where(dabs(ab(hxC^S)-one)>smalldouble)
             wLC(hxC^S,iw) = ab(hxC^S)*wLC(hxC^S,iw)+wMax(hxC^S,iw)
          endwhere

          where(dabs(ab(hxR^S)-one)>smalldouble)
             wRC(hxC^S,iw) = ab(hxR^S)*wRC(hxC^S,iw)+wMax(hxR^S,iw)
          endwhere
       enddo
    endif

  end subroutine PPMlimiter

end module mod_ppm
