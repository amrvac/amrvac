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
    ! old version : april 2009
    ! author: zakaria.meliani@wis.kuleuven.be
    ! current version : March 2023 by Chun Xia

    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ix^L, idims
    double precision, intent(in)    :: q(ixI^S),qCT(ixI^S)

    double precision, intent(inout) :: qRC(ixI^S),qLC(ixI^S)

    double precision,dimension(ixI^S)  :: dqC,d2qC,ldq
    double precision,dimension(ixI^S)  :: qMin,qMax,tmp

    integer   :: lxC^L,lxR^L,lxL^L
    integer   :: ixL^L,ixO^L,ixR^L,ixOL^L,ixOR^L
    integer   :: hxL^L,hxC^L,hxR^L
    integer   :: kxL^L,kxC^L,kxR^L

    ixOmin^D=ixmin^D-kr(idims,^D);ixOmax^D=ixmax^D+kr(idims,^D);!ixO[ixMmin-1,ixMmax+1]
    ixOLmin^D=ixOmin^D-kr(idims,^D);ixOLmax^D=ixOmax^D;         !ixOL[ixMmin-2,ixMmax+1]
    ixOR^L=ixOL^L+kr(idims,^D);                                 !ixOR=[iMmin-1,ixMmax+2]
    ixL^L=ixO^L-kr(idims,^D);                                   !ixL[ixMmin-2,ixMmax]
    ixR^L=ixO^L+kr(idims,^D);                                   !ixR=[iMmin,ixMmax+2]

    hxCmin^D=ixOmin^D;hxCmax^D=ixmax^D; ! hxC = [ixMmin-1,ixMmax]
    hxL^L=hxC^L-kr(idims,^D);           ! hxL = [ixMmin-2,ixMmax-1]
    hxR^L=hxC^L+kr(idims,^D);           ! hxR = [ixMmin,ixMmax+1]

    kxCmin^D=ixLmin^D-1; kxCmax^D=ixRmax^D; ! kxC=[iMmin-3,ixMmax+2]
    kxR^L=kxC^L+kr(idims,^D);              ! kxR=[iMmin-2,ixMmax+3]

    lxCmin^D=ixOmin^D-kr(idims,^D);lxCmax^D=ixOmax^D+kr(idims,^D);! lxC=[iMmin-2,ixMmax+2]
    lxL^L=lxC^L-kr(idims,^D);                          ! lxL=[iMmin-3,ixMmax+1]
    lxR^L=lxC^L+kr(idims,^D);                          ! lxR=[iMmin-1,ixMmax+3]

    dqC(kxC^S)=q(kxR^S)-q(kxC^S)
    ! Eq. 64,  Miller and Colella 2002, JCP 183, 26 
    d2qC(kxC^S)=half*(q(lxR^S)-q(lxL^S))
    where(dqC(lxC^S)*dqC(lxL^S)>zero)
       ! Store the sign of dwC in wMin
       qMin(kxC^S)= sign(one,d2qC(lxC^S))
       ! Eq. 65,  Miller and Colella 2002, JCP 183, 26 
       ldq(lxC^S)= qMin(lxC^S)*min(dabs(d2qC(lxC^S)),&
            2.0d0*dabs(dqC(lxL^S)),&
            2.0d0*dabs(dqC(lxC^S)))
    else where
       ldq(lxC^S)=zero
    end where

    ! Eq. 66, Miller and Colella 2002, JCP 183, 26
    qLC(ixOL^S)=qLC(ixOL^S)+half*dqC(ixOL^S)&
         +(ldq(ixOL^S)-ldq(ixOR^S))/6.0d0
    qRC(ixL^S)=qLC(ixL^S)

    ! make sure that min wCT(i)<wLC(i)<wCT(i+1) same for wRC(i)
    call extremaq(ixI^L,ixO^L,qCT,1,qMax,qMin)

    ! Eq. B8, page 217, Mignone et al 2005, ApJS
    qRC(ixL^S)=max(qMin(ixO^S),min(qMax(ixO^S),qRC(ixL^S)))
    qLC(ixO^S)=max(qMin(ixO^S),min(qMax(ixO^S),qLC(ixO^S)))

    ! Eq. B9, page 217, Mignone et al 2005, ApJS
    where((qRC(ixL^S)-qCT(ixO^S))*(qCT(ixO^S)-qLC(ixO^S))<=zero)
       qRC(ixL^S)=qCT(ixO^S)
       qLC(ixO^S)=qCT(ixO^S)
    end where

    qMax(ixO^S)=(qLC(ixO^S)-qRC(ixL^S))&
         *(qCT(ixO^S)-half*(qLC(ixO^S)+qRC(ixL^S)))
    qMin(ixO^S)=(qLC(ixO^S)-qRC(ixL^S))**2/6.0d0

    tmp(hxL^S)=qRC(hxL^S)
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
    ! old version : april 2009
    ! author: zakaria.meliani@wis.kuleuven.be
    ! current version : March 2023 by Chun Xia
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ix^L, idims
    double precision, intent(in)    :: w(ixI^S,1:nw),wCT(ixI^S,1:nw)

    double precision, intent(inout) :: wRC(ixI^S,1:nw),wLC(ixI^S,1:nw) 

    double precision,dimension(ixI^S,1:nw_recon)  :: dwC,d2wC,ldw
    double precision,dimension(ixI^S,1:nw_recon)  :: wMin,wMax,tmp
    double precision,dimension(ixI^S) :: aa, ab, ac, ki

    double precision, parameter :: betamin=0.75d0, betamax=0.85d0,&
         Zmin=0.25d0, Zmax=0.75d0,&
         eta1=20.0d0,eta2=0.05d0,eps=0.01d0,kappa=0.1d0
    integer   :: lxC^L,lxL^L,lxR^L
    integer   :: ixLL^L,ixL^L,ixO^L,ixR^L,ixRR^L,ixOL^L,ixOR^L
    integer   :: hxL^L,hxC^L,hxR^L
    integer   :: kxLL^L,kxL^L,kxC^L,kxR^L,kxRR^L
    integer   :: iw, idimss

    ixOmin^D=ixmin^D-kr(idims,^D);ixOmax^D=ixmax^D+kr(idims,^D);!ixO[ixMmin-1,ixMmax+1]
    ixOLmin^D=ixOmin^D-kr(idims,^D);ixOLmax^D=ixOmax^D;         !ixOL[ixMmin-2,ixMmax+1]
    ixOR^L=ixOL^L+kr(idims,^D);                                 !ixOR=[iMmin-1,ixMmax+2]
    ixL^L=ixO^L-kr(idims,^D);                                   !ixL[ixMmin-2,ixMmax]
    ixR^L=ixO^L+kr(idims,^D);                                   !ixR=[iMmin,ixMmax+2]

    hxCmin^D=ixOmin^D;hxCmax^D=ixmax^D; ! hxC = [ixMmin-1,ixMmax]
    hxL^L=hxC^L-kr(idims,^D);           ! hxL = [ixMmin-2,ixMmax-1]
    hxR^L=hxC^L+kr(idims,^D);           ! hxR = [ixMmin,ixMmax+1]

    kxCmin^D=ixLmin^D-1; kxCmax^D=ixRmax^D; ! kxC=[iMmin-3,ixMmax+2]
    kxR^L=kxC^L+kr(idims,^D);              ! kxR=[iMmin-2,ixMmax+3]

    lxCmin^D=ixOmin^D-kr(idims,^D);lxCmax^D=ixOmax^D+kr(idims,^D);! lxC=[iMmin-2,ixMmax+2]
    lxL^L=lxC^L-kr(idims,^D);                          ! lxL=[iMmin-3,ixMmax+1]
    lxR^L=lxC^L+kr(idims,^D);                          ! lxR=[iMmin-1,ixMmax+3]

    dwC(kxC^S,1:nw_recon)=w(kxR^S,1:nw_recon)-w(kxC^S,1:nw_recon)
    ! Eq. 64,  Miller and Colella 2002, JCP 183, 26 
    d2wC(lxC^S,1:nw_recon)=half*(w(lxR^S,1:nw_recon)-w(lxL^S,1:nw_recon))
    where(dwC(lxC^S,1:nw_recon)*dwC(lxL^S,1:nw_recon)>zero)
       ! Store the sign of dwC in wMin
       wMin(lxC^S,1:nw_recon)= sign(one,d2wC(lxC^S,1:nw_recon))
       ! Eq. 65,  Miller and Colella 2002, JCP 183, 26 
       ldw(lxC^S,1:nw_recon)= wMin(lxC^S,1:nw_recon)*min(dabs(d2wC(lxC^S,1:nw_recon)),&
            2.0d0*dabs(dwC(lxL^S,1:nw_recon)),&
            2.0d0*dabs(dwC(lxC^S,1:nw_recon)))
    else where
       ldw(lxC^S,1:nw_recon)=zero
    end where

    ! Eq. 66,  Miller and Colella 2002, JCP 183, 26 
    wLC(ixOL^S,1:nw_recon)=wLC(ixOL^S,1:nw_recon)+half*dwC(ixOL^S,1:nw_recon)&
         +(ldw(ixOL^S,1:nw_recon)-ldw(ixOR^S,1:nw_recon))/6.0d0
    wRC(ixL^S,1:nw_recon)=wLC(ixL^S,1:nw_recon)

    ! make sure that min w(i)<wLC(i)<w(i+1) same for wRC(i)
    call extremaw(ixI^L,ixO^L,w,1,wMax,wMin)

    ! Eq. B8, page 217, Mignone et al 2005, ApJS
    wRC(ixL^S,1:nw_recon)=max(wMin(ixO^S,1:nw_recon)&
         ,min(wMax(ixO^S,1:nw_recon),wRC(ixL^S,1:nw_recon))) 
    wLC(ixO^S,1:nw_recon)=max(wMin(ixO^S,1:nw_recon)&
         ,min(wMax(ixO^S,1:nw_recon),wLC(ixO^S,1:nw_recon)))

    ! Eq. B9, page 217, Mignone et al 2005, ApJS
    where((wRC(ixL^S,1:nw_recon)-wCT(ixO^S,1:nw_recon))&
         *(wCT(ixO^S,1:nw_recon)-wLC(ixO^S,1:nw_recon))<=zero)
       wRC(ixL^S,1:nw_recon)=wCT(ixO^S,1:nw_recon)
       wLC(ixO^S,1:nw_recon)=wCT(ixO^S,1:nw_recon)
    end where

    wMax(ixO^S,1:nw_recon)=(wLC(ixO^S,1:nw_recon)-wRC(ixL^S,1:nw_recon))*&
         (wCT(ixO^S,1:nw_recon)-half*(wLC(ixO^S,1:nw_recon)+wRC(ixL^S,1:nw_recon)))
    wMin(ixO^S,1:nw_recon)=(wLC(ixO^S,1:nw_recon)-wRC(ixL^S,1:nw_recon))**2/6.0d0
    tmp(hxL^S,1:nw_recon)=wRC(hxL^S,1:nw_recon)
    ! Eq. B10, page 218, Mignone et al 2005, ApJS
    where(wMax(hxR^S,1:nw_recon)>wMin(hxR^S,1:nw_recon))
       wRC(hxC^S,1:nw_recon)= 3.0d0*wCT(hxR^S,1:nw_recon)&
            -2.0d0*wLC(hxR^S,1:nw_recon)
    end where
    ! Eq. B11, page 218, Mignone et al 2005, ApJS
    where(wMax(hxC^S,1:nw_recon)<-wMin(hxC^S,1:nw_recon))
       wLC(hxC^S,1:nw_recon)= 3.0d0*wCT(hxC^S,1:nw_recon)&
            -2.0d0*tmp(hxL^S,1:nw_recon)
    end where

    ! flattening at the contact discontinuities
    if(flatcd)then
      ixRR^L=ixR^L+kr(idims,^D);               !ixRR=[iMmin+1,ixMmax+3]
      kxL^L=kxC^L-kr(idims,^D);                ! kxL=[iMmin-4,ixMmax+1]
      call ppm_flatcd(ixI^L,kxC^L,kxL^L,kxR^L,wCT,d2wC,aa,ab)
      if(any(kappa*aa(kxC^S)>=ab(kxC^S)))then
        do iw=1,nw_recon
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
          else where
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
      end if
    end if

    ! flattening at the shocks
    if(flatsh)then
       ! following MILLER and COLELLA 2002 JCP 183, 26
       kxCmin^D=ixmin^D-2; kxCmax^D=ixmax^D+2; ! kxC=[ixMmin-2,ixMmax+2]
       ki=bigdouble
       do idimss=1,ndim
          kxL^L=kxC^L-kr(idimss,^D); ! kxL=[ixMmin-3,ixMmax+1]
          kxR^L=kxC^L+kr(idimss,^D); ! kxR=[ixMmin-1,ixMmax+3]
          kxLL^L=kxL^L-kr(idimss,^D);! kxLL=[ixMmin-4,ixMmax]
          kxRR^L=kxR^L+kr(idimss,^D);! kxRR=[ixMmin,ixMmax+4]

          call ppm_flatsh(ixI^L,kxC^L,kxLL^L,kxL^L,kxR^L,kxRR^L,idimss,wCT,aa,ab)

          ! eq. B17, page 218, Mignone et al 2005, ApJS (Xi1 min)
          ac(kxC^S) = max(zero,min(one,(betamax-aa(kxC^S))/(betamax-betamin)))
          ! eq. B18, page 218, Mignone et al 2005, ApJS (Xi1)
          ! recycling aa(ixL^S)
          where(wCT(kxR^S,iw_mom(idimss))<wCT(kxL^S,iw_mom(idimss)))
             aa(kxC^S) = max(ac(kxC^S), min(one,(Zmax-ab(kxC^S))/(Zmax-Zmin)))
          else where
             aa(kxC^S) = one
          end where
          ixL^L=ixO^L-kr(idimss,^D); ! ixL=[ixMmin-2,ixMmax]
          ixR^L=ixO^L+kr(idimss,^D); ! ixR=[ixMmin,ixMmax+2]
          ki(ixO^S)=min(ki(ixO^S),aa(ixL^S),aa(ixO^S),aa(ixR^S))
       end do
       ! recycling wMax
       do iw=1,nw_recon
          where(dabs(ab(ixO^S)-one)>smalldouble)
             wMax(ixO^S,iw) = (one-ab(ixO^S))*wCT(ixO^S,iw)
          end where

          where(dabs(ab(hxC^S)-one)>smalldouble)
             wLC(hxC^S,iw) = ab(hxC^S)*wLC(hxC^S,iw)+wMax(hxC^S,iw)
          end where

          where(dabs(ab(hxR^S)-one)>smalldouble)
             wRC(hxC^S,iw) = ab(hxR^S)*wRC(hxC^S,iw)+wMax(hxR^S,iw)
          end where
       end do
    end if

  end subroutine PPMlimiter

  subroutine ppm_flatcd(ixI^L,ixO^L,ixL^L,ixR^L,w,d2w,drho,dp)
    use mod_global_parameters
    use mod_physics, only: phys_gamma

    integer, intent(in)             :: ixI^L, ixO^L, ixL^L, ixR^L
    double precision, intent(in)    :: w(ixI^S, nw), d2w(ixG^T, 1:nw_recon)
    double precision, intent(inout) :: drho(ixG^T), dp(ixG^T)

    drho(ixO^S) = phys_gamma*abs(d2w(ixO^S, iw_rho))&
         /min(w(ixL^S, iw_rho), w(ixR^S, iw_rho))
    dp(ixO^S) = abs(d2w(ixO^S, iw_e))/min(w(ixL^S, iw_e), w(ixR^S, iw_e))

  end subroutine ppm_flatcd

  subroutine ppm_flatsh(ixI^L,ixO^L,ixLL^L,ixL^L,ixR^L,ixRR^L,idims,w,drho,dp)
    use mod_global_parameters
    use mod_physics, only: phys_gamma

    integer, intent(in)             :: ixI^L, ixO^L, ixLL^L, ixL^L, ixR^L, ixRR^L
    integer, intent(in)             :: idims
    double precision, intent(in)    :: w(ixI^S, nw)
    double precision, intent(inout) :: drho(ixI^S), dp(ixI^S)

    ! eq. B15, page 218, Mignone and Bodo 2005, ApJS (beta1)
    where (abs(w(ixRR^S, iw_e)-w(ixLL^S, iw_e))>smalldouble)
       drho(ixO^S) = abs((w(ixR^S, iw_e)-w(ixL^S, iw_e))&
            /(w(ixRR^S, iw_e)-w(ixLL^S, iw_e)))
    else where
       drho(ixO^S) = zero
    end where

    ! eq. 76, page 48, Miller and Colella 2002, JCoPh, adjusted
    dp(ixO^S) = abs(w(ixR^S, iw_e)-w(ixL^S, iw_e))&
         /(phys_gamma*(w(ixR^S, iw_e)+w(ixL^S, iw_e)))

  end subroutine ppm_flatsh

  subroutine extremaq(ixI^L,ixO^L,q,nshift,qMax,qMin)

    use mod_global_parameters

    integer,intent(in)           :: ixI^L,ixO^L
    double precision, intent(in) :: q(ixI^S)
    integer,intent(in)           :: nshift

    double precision, intent(out) :: qMax(ixI^S),qMin(ixI^S)

    integer           :: ixs^L,ixsR^L,ixsL^L,idims,jdims,kdims,ishift,i,j

    do ishift=1,nshift
      idims=1
      ixsR^L=ixO^L+ishift*kr(idims,^D);
      ixsL^L=ixO^L-ishift*kr(idims,^D);
      if (ishift==1) then
        qMax(ixO^S)=max(q(ixO^S),q(ixsR^S),q(ixsL^S))
        qMin(ixO^S)=min(q(ixO^S),q(ixsR^S),q(ixsL^S))
      else
        qMax(ixO^S)=max(qMax(ixO^S),q(ixsR^S),q(ixsL^S))
        qMin(ixO^S)=min(qMin(ixO^S),q(ixsR^S),q(ixsL^S))
      end if
      {^NOONED
      idims=1
      jdims=idims+1
      do i=-1,1
        ixs^L=ixO^L+i*ishift*kr(idims,^D);
        ixsR^L=ixs^L+ishift*kr(jdims,^D);
        ixsL^L=ixs^L-ishift*kr(jdims,^D);
        qMax(ixO^S)=max(qMax(ixO^S),q(ixsR^S),q(ixsL^S))
        qMin(ixO^S)=min(qMin(ixO^S),q(ixsR^S),q(ixsL^S))
      end do
      }
      {^IFTHREED
      idims=1
      jdims=idims+1
      kdims=jdims+1
      do i=-1,1
        ixs^L=ixO^L+i*ishift*kr(idims,^D);
        do j=-1,1
          ixs^L=ixO^L+j*ishift*kr(jdims,^D);
          ixsR^L=ixs^L+ishift*kr(kdims,^D);
          ixsL^L=ixs^L-ishift*kr(kdims,^D);
          qMax(ixO^S)=max(qMax(ixO^S),q(ixsR^S),q(ixsL^S))
          qMin(ixO^S)=min(qMin(ixO^S),q(ixsR^S),q(ixsL^S))
        end do
      end do
      }
    enddo

  end subroutine  extremaq

  subroutine extremaw(ixI^L,ixO^L,w,nshift,wMax,wMin)
    use mod_global_parameters

    integer,intent(in)            :: ixI^L,ixO^L
    double precision, intent(in)  :: w(ixI^S,1:nw)
    integer,intent(in)            :: nshift

    double precision, intent(out) :: wMax(ixI^S,1:nw_recon),wMin(ixI^S,1:nw_recon)

    integer          :: ixs^L,ixsR^L,ixsL^L,idims,jdims,kdims,ishift,i,j

    do ishift=1,nshift
      idims=1
      ixsR^L=ixO^L+ishift*kr(idims,^D);
      ixsL^L=ixO^L-ishift*kr(idims,^D);
      if (ishift==1) then
        wMax(ixO^S,1:nw_recon)= &
             max(w(ixO^S,1:nw_recon),w(ixsR^S,1:nw_recon),w(ixsL^S,1:nw_recon))
        wMin(ixO^S,1:nw_recon)= &
             min(w(ixO^S,1:nw_recon),w(ixsR^S,1:nw_recon),w(ixsL^S,1:nw_recon))
      else
        wMax(ixO^S,1:nw_recon)= &
             max(wMax(ixO^S,1:nw_recon),w(ixsR^S,1:nw_recon),w(ixsL^S,1:nw_recon))
        wMin(ixO^S,1:nw_recon)= &
             min(wMin(ixO^S,1:nw_recon),w(ixsR^S,1:nw_recon),w(ixsL^S,1:nw_recon))
      end if
      {^NOONED
      idims=1
      jdims=idims+1
      do i=-1,1
        ixs^L=ixO^L+i*ishift*kr(idims,^D);
        ixsR^L=ixs^L+ishift*kr(jdims,^D);
        ixsL^L=ixs^L-ishift*kr(jdims,^D);
        wMax(ixO^S,1:nw_recon)= &
             max(wMax(ixO^S,1:nw_recon),w(ixsR^S,1:nw_recon),w(ixsL^S,1:nw_recon))
        wMin(ixO^S,1:nw_recon)= &
             min(wMin(ixO^S,1:nw_recon),w(ixsR^S,1:nw_recon),w(ixsL^S,1:nw_recon))
      end do
      }
      {^IFTHREED
      idims=1
      jdims=idims+1
      kdims=jdims+1
      do i=-1,1
        ixs^L=ixO^L+i*ishift*kr(idims,^D);
        do j=-1,1
          ixs^L=ixO^L+j*ishift*kr(jdims,^D);
          ixsR^L=ixs^L+ishift*kr(kdims,^D);
          ixsL^L=ixs^L-ishift*kr(kdims,^D);
          wMax(ixO^S,1:nw_recon)= &
               max(wMax(ixO^S,1:nw_recon),w(ixsR^S,1:nw_recon),w(ixsL^S,1:nw_recon))
          wMin(ixO^S,1:nw_recon)= &
               min(wMin(ixO^S,1:nw_recon),w(ixsR^S,1:nw_recon),w(ixsL^S,1:nw_recon))
        end do
      end do
      }
    enddo

  end subroutine  extremaw

end module mod_ppm
