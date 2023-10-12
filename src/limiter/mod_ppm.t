module mod_ppm

  implicit none
  private

  public :: PPMlimiter
  public :: PPMflattening

contains

  subroutine PPMlimiter(ixI^L,hxC^L,idims,w,wL,wR,extrema)

    ! references:
    ! Mignone et al 2005, ApJS 160, 199, 
    ! Miller and Colella 2002, JCP 183, 26 
    ! Fryxell et al. 2000 ApJ, 131, 273 (Flash)
    ! baciotti Phd (http://www.aei.mpg.de/~baiotti/Baiotti_PhD.pdf)
    ! version : Nov 2022
    ! author: zakaria.meliani@wis.kuleuven.be, and 
    !         patrick.cheong@berkeley.edu

    use mod_global_parameters

    integer, intent(in)             :: ixI^L, hxC^L, idims
    double precision, intent(in)    :: w(ixI^S)
    double precision, intent(inout) :: wR(ixI^S),wL(ixI^S) 
    logical, intent(in)             :: extrema

    double precision, dimension(ixI^S)  :: dwC,d2wC,ldw
    double precision, dimension(ixI^S)  :: wMin,wMax,tmp
    integer   :: lxC^L,lxR^L
    integer   :: ixLL^L,ixL^L,ixO^L,ixR^L,ixRR^L
    integer   :: hxL^L,hxR^L
    integer   :: kxL^L,kxC^L,kxR^L

    ixOmin^D=hxCmin^D; ixOmax^D=hxCmax^D+kr(idims,^D); ! ixO[ixMmin1-1,ixMmax1+1]
    ixL^L=ixO^L-kr(idims,^D);                          ! ixL[ixMmin1-2,ixMmax1]
    ixLL^L=ixL^L-kr(idims,^D);                         ! ixLL[ixMmin1-3,ixMmax1-1]
    ixR^L=ixO^L+kr(idims,^D);                          ! ixR=[iMmin1,ixMmax+2]
    ixRR^L=ixR^L+kr(idims,^D);                         ! ixRR=[iMmin1+1,ixMmax+3]

    hxL^L=hxC^L-kr(idims,^D);                          ! hxL = [ixMmin-2,ixMmax-1]
    hxR^L=hxC^L+kr(idims,^D);                          ! hxR = [ixMmin,ixMmax+1]

    kxCmin^D=ixLLmin^D; kxCmax^D=ixRmax^D;             ! kxC=[iMmin1-3,ixMmax1+2]
    kxL^L=kxC^L-kr(idims,^D);                          ! kxL=[iMmin1-4,ixMmax1+1]
    kxR^L=kxC^L+kr(idims,^D);                          ! kxR=[iMmin1-2,ixMmax1+3]

    lxCmin^D=ixLLmin^D-kr(idims,^D);lxCmax^D=ixRRmax^D;! ixC=[iMmin1-4,ixMmax1+3]
    lxR^L=lxC^L+kr(idims,^D);                          ! lxR=[iMmin1-3,ixMmax1+4]

    dwC(lxC^S)=w(lxR^S)-w(lxC^S)
    ! Eq. 64,  Miller and Colella 2002, JCP 183, 26 
    d2wC(kxC^S)=half*(w(kxR^S)-w(kxL^S))
    where(dwC(kxC^S)*dwC(kxL^S)>zero)
       ! Store the sign of dwC in wMin
       wMin(kxC^S)= sign(one,d2wC(kxC^S))
       ! Eq. 65,  Miller and Colella 2002, JCP 183, 26 
       ldw(kxC^S)= wMin(kxC^S)*min(dabs(d2wC(kxC^S)),&
            2.0d0*dabs(dwC(kxL^S)),&
            2.0d0*dabs(dwC(kxC^S)))
    elsewhere
       ldw(kxC^S)=zero
    endwhere

    ! Eq. 66,  Miller and Colella 2002, JCP 183, 26 
    wL(ixO^S)=wL(ixO^S)+(half*dwC(ixO^S)+(ldw(ixO^S)-ldw(ixR^S))/6.0d0)
    wR(ixL^S)=wR(ixL^S)-(half*dwC(ixL^S)-(ldw(ixL^S)-ldw(ixO^S))/6.0d0)
   
    if (extrema) then
       ! Note: the following part needs surrounding 9 or 27 cells in 2D or 3D
       ! make sure that min w(i)<wL(i)<w(i+1) same for wR(i)
       call extremaw(ixI^L,kxC^L,w,1,wMax,wMin)
       ! Eq. B8, page 217, Mignone et al 2005, ApJS
       wR(ixL^S)=max(wMin(ixO^S),min(wMax(ixO^S),wR(ixL^S))) 
       wL(ixO^S)=max(wMin(ixO^S),min(wMax(ixO^S),wL(ixO^S)))
    end if

    ! Eq. B9, page 217, Mignone et al 2005, ApJS
    where((wR(ixL^S)-w(ixO^S))&
         *(w(ixO^S)-wL(ixO^S))<=zero)
       wR(ixL^S)=w(ixO^S)
       wL(ixO^S)=w(ixO^S)
    end where

    wMax(ixO^S)=(wL(ixO^S)-wR(ixL^S))*&
         (w(ixO^S)-&
         (wL(ixO^S)+wR(ixL^S))/2.0d0)
    wMin(ixO^S)=(wL(ixO^S)-wR(ixL^S))**2.0d0/6.0d0
    tmp(ixL^S)=wR(ixL^S)
    ! Eq. B10, page 218, Mignone et al 2005, ApJS
    where(wMax(hxR^S)>wMin(hxR^S))
       wR(hxC^S)= 3.0d0*w(hxR^S)-2.0d0*wL(hxR^S)
    endwhere
    ! Eq. B11, page 218, Mignone et al 2005, ApJS
    where(wMax(hxC^S)<-wMin(hxC^S))
       wL(hxC^S)= 3.0d0*w(hxC^S)-2.0d0*tmp(hxL^S)
    endwhere

  end subroutine PPMlimiter

  subroutine PPMflattening(ixI^L,hxC^L,idims,w,wL,wR,flat_cd,flat_sh)

    use mod_global_parameters
    use mod_physics, only: phys_ppm_flatsh, phys_ppm_flatcd

    integer, intent(in)             :: ixI^L, hxC^L, idims
    double precision, intent(in)    :: w(ixI^S,1:nw)
    double precision, intent(inout) :: wR(ixI^S,1:nw),wL(ixI^S,1:nw) 
    logical, intent(in)             :: flat_cd, flat_sh

    double precision, dimension(ixI^S,1:nw)  :: dwC,d2wC,ldw
    double precision, dimension(ixI^S,1:nw)  :: wMin,wMax
    double precision, dimension(ixI^S)          :: aa, ab, ac, dv
    double precision, dimension(ixI^S,1:ndim)   ::  exi
    double precision, parameter :: betamin=0.75d0, betamax=0.85d0,&
                                   Zmin=0.25d0, Zmax=0.75d0,&
                                   eta1=20.0d0,eta2=0.05d0,eps1=0.01d0,kappa=0.1d0
    integer   :: lxC^L,lxR^L
    integer   :: ix^L
    integer   :: ixLL^L,ixL^L,ixO^L,ixR^L,ixRR^L
    integer   :: hxL^L,hxR^L
    integer   :: kxLL^L,kxL^L,kxC^L,kxR^L,kxRR^L
    integer   :: iw, idimss

    if ( (.not.flat_cd) .and. (.not. flat_sh) ) return

    ! Note: flatteneing need to be done after all variables are reconstructed,
    ! we apply it after all PPM are done

    ixOmin^D=hxCmin^D; ixOmax^D=hxCmax^D+kr(idims,^D); ! ixO[ixMmin1-1,ixMmax1+1]
    ixL^L=ixO^L-kr(idims,^D);                          ! ixL[ixMmin1-2,ixMmax1]
    ixLL^L=ixL^L-kr(idims,^D);                         ! ixLL[ixMmin1-3,ixMmax1-1]
    ixR^L=ixO^L+kr(idims,^D);                          ! ixR=[iMmin1,ixMmax+2]
    ixRR^L=ixR^L+kr(idims,^D);                         ! ixRR=[iMmin1+1,ixMmax+3]

    hxL^L=hxC^L-kr(idims,^D);                          ! hxL = [ixMmin-2,ixMmax-1]
    hxR^L=hxC^L+kr(idims,^D);                          ! hxR = [ixMmin,ixMmax+1]

    kxCmin^D=ixLLmin^D; kxCmax^D=ixRmax^D;             ! kxC=[iMmin1-3,ixMmax1+2]
    kxL^L=kxC^L-kr(idims,^D);                          ! kxL=[iMmin1-4,ixMmax1+1]
    kxR^L=kxC^L+kr(idims,^D);                          ! kxR=[iMmin1-2,ixMmax1+3]

    ixmin^D=ixRmin^D; ixmax^D=ixLmax^D;

    lxCmin^D=ixLLmin^D-kr(idims,^D);lxCmax^D=ixRRmax^D;! ixC=[iMmin1-4,ixMmax1+3]
    lxR^L=lxC^L+kr(idims,^D);                          ! lxR=[iMmin1-3,ixMmax1+4]

    dwC(lxC^S,1:nw)=w(lxR^S,1:nw)-w(lxC^S,1:nw)
    ! Eq. 64,  Miller and Colella 2002, JCP 183, 26 
    d2wC(kxC^S,1:nw)=half*(w(kxR^S,1:nw)-w(kxL^S,1:nw))
    where(dwC(kxC^S,1:nw)*dwC(kxL^S,1:nw)>zero)
       ! Store the sign of dwC in wMin
       wMin(kxC^S,1:nw)= sign(one,d2wC(kxC^S,1:nw))
       ! Eq. 65,  Miller and Colella 2002, JCP 183, 26 
       ldw(kxC^S,1:nw)= wMin(kxC^S,1:nw)*min(dabs(d2wC(kxC^S,1:nw)),&
            2.0d0*dabs(dwC(kxL^S,1:nw)),&
            2.0d0*dabs(dwC(kxC^S,1:nw)))
    elsewhere
       ldw(kxC^S,1:nw)=zero
    endwhere

    ! flattening at the contact discontinuities
    if(flat_cd)then
       ! Note that aa here is drho, ab is dpress
       call phys_ppm_flatcd(ixI^L,kxC^L,kxL^L,kxR^L,w,d2wC,aa,ab)
       if(any(kappa*aa(kxC^S)>=ab(kxC^S)))then
          do iw=1,nw
             ! metric cannot be modified here
             if ( iw >= nmetric_lo .and. iw <= nmetric_hi ) cycle

             !where(kappa*aa(kxC^S)>=ab(kxC^S).and. dabs(dwC(kxC^S,iw))>smalldouble)
                wMax(kxC^S,iw) = w(kxR^S,iw)-2.0d0*w(kxC^S,iw)+w(kxL^S,iw)
             !end where
   
             where(wMax(ixR^S,iw)*wMax(ixL^S,iw)<zero &
                  .and. dabs(w(ixR^S,iw)-w(ixL^S,iw)) &
                        -eps1*min(dabs(w(ixR^S,iw)),dabs(w(ixL^S,iw)))>zero &
                  .and. kappa*aa(ixO^S)>=ab(ixO^S) &
                  .and. dabs(dwC(ixO^S,iw))>smalldouble)
   
                ac(ixO^S)=(w(ixLL^S,iw)-w(ixRR^S,iw)+4.0d0*dwC(ixO^S,iw))&
                     /(12.0d0*dwC(ixO^S,iw))
                wMin(ixO^S,iw)=max(zero,min(eta1*(ac(ixO^S)-eta2),one))
             elsewhere
                wMin(ixO^S,iw)=zero
             end where
   
             where(wMin(hxC^S,iw)>zero)
                wL(hxC^S,iw) = wL(hxC^S,iw)*(one-wMin(hxC^S,iw))&
                     +(w(hxC^S,iw)+half*ldw(hxC^S,iw))*wMin(hxC^S,iw)
             end where
             where(wMin(hxR^S,iw)>zero)
                wR(hxC^S,iw) = wR(hxC^S,iw)*(one-wMin(hxR^S,iw))&
                     +(w(hxR^S,iw)-half*ldw(hxR^S,iw))*wMin(hxR^S,iw)
             end where
          end do
       endif
    endif

    ! flattening at the shocks
    if(flat_sh)then
       ! following MILLER and COLELLA 2002 JCP 183, 26
       kxCmin^D=ixmin^D-2; kxCmax^D=ixmax^D+2; ! kxC=[ixMmin1-2,ixMmax1+2]
       do idimss=1,ndim
          kxL^L=kxC^L-kr(idimss,^D); ! kxL=[ixMmin1-3,ixMmax1+1]
          kxR^L=kxC^L+kr(idimss,^D); ! kxR=[ixMmin1-1,ixMmax1+3]
          kxLL^L=kxL^L-kr(idimss,^D);! kxLL=[ixMmin-4,ixMmax]
          kxRR^L=kxR^L+kr(idimss,^D);! kxRR=[ixMmin,ixMmax+4]

          ! Note that aa here is betai, ab is zi
          call phys_ppm_flatsh(ixI^L,kxC^L,kxLL^L,kxL^L,kxR^L,kxRR^L,idimss,w,aa,ab,dv)

          ! eq. B17, page 218, Mignone et al 2005, ApJS (had(Xi1))
          ac(kxC^S) = max(zero,min(one,(betamax-aa(kxC^S))/(betamax-betamin)))
          ! eq. B18, page 218, Mignone et al 2005, ApJS (had(Xi1))
          ! recycling aa(ixL^S)
          where (dv(kxC^S)<smalldouble)
             aa(kxC^S) = max(ac(kxC^S), min(one,(Zmax-ab(kxC^S))/(Zmax-Zmin)))
          elsewhere
             aa(kxC^S) = one
          endwhere

          {^IFONED call extremaa(ixI^L,ixO^L,aa,1,ab)}
          {^NOONED call extremaa(ixI^L,ixO^L,aa,1,exi(ixI^S,idimss))}
       end do
       {^NOONED ab(ixO^S)=min(^D&exi(ixO^S,^D))}
       ! recycling wMax
       do iw=1,nw
          ! metric cannot be modified here
          if ( iw >= nmetric_lo .and. iw <= nmetric_hi ) cycle

          where(dabs(ab(ixO^S)-one)>smalldouble)
             wMax(ixO^S,iw) = (one-ab(ixO^S))*w(ixO^S,iw)
          endwhere
   
          where(dabs(ab(hxC^S)-one)>smalldouble)
             wL(hxC^S,iw) = ab(hxC^S)*wL(hxC^S,iw)+wMax(hxC^S,iw)
          endwhere
   
          where(dabs(ab(hxR^S)-one)>smalldouble)
             wR(hxC^S,iw) = ab(hxR^S)*wR(hxC^S,iw)+wMax(hxR^S,iw)
          endwhere
       end do
    endif
  end subroutine PPMflattening

  subroutine extremaa(ixI^L,ixO^L,a,nshift,aMin)
    use mod_global_parameters

    integer,intent(in)            :: ixI^L,ixO^L
    double precision, intent(in)  :: a(ixI^S)
    integer,intent(in)            :: nshift
    double precision, intent(out) :: aMin(ixI^S)

    integer          :: ixs^L,ixsR^L,ixsL^L,idims,jdims,kdims,ishift,i,j

    do ishift=1,nshift
      idims=1
      ixsR^L=ixO^L+ishift*kr(idims,^D);
      ixsL^L=ixO^L-ishift*kr(idims,^D);
      aMin(ixO^S)=min(a(ixsR^S),a(ixO^S),a(ixsL^S))
      {^NOONED
      idims=1
      jdims=idims+1
      do i=-1,1
        ixs^L=ixO^L+i*ishift*kr(idims,^D);
        ixsR^L=ixs^L+ishift*kr(jdims,^D);
        ixsL^L=ixs^L-ishift*kr(jdims,^D);
        aMin(ixO^S)=min(aMin(ixO^S),a(ixsR^S),a(ixsL^S))
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
          aMin(ixO^S)=min(aMin(ixO^S),a(ixsR^S),a(ixsL^S))
        end do
      end do
      }
    end do

  end subroutine extremaa

  subroutine extremaw(ixI^L,ixO^L,w,nshift,wMax,wMin)
    use mod_global_parameters

    integer,intent(in)            :: ixI^L,ixO^L
    double precision, intent(in)  :: w(ixI^S)
    integer,intent(in)            :: nshift

    double precision, intent(out) :: wMax(ixI^S),wMin(ixI^S)

    integer          :: ixs^L,ixsR^L,ixsL^L,idims,jdims,kdims,ishift,i,j

    do ishift=1,nshift
      idims=1
      ixsR^L=ixO^L+ishift*kr(idims,^D);
      ixsL^L=ixO^L-ishift*kr(idims,^D);
      if (ishift==1) then
        wMax(ixO^S)= max(w(ixO^S),w(ixsR^S),w(ixsL^S))
        wMin(ixO^S)= min(w(ixO^S),w(ixsR^S),w(ixsL^S))
      else
        wMax(ixO^S)= max(wMax(ixO^S),w(ixsR^S),w(ixsL^S))
        wMin(ixO^S)= min(wMin(ixO^S),w(ixsR^S),w(ixsL^S))
      end if
      {^NOONED
      idims=1
      jdims=idims+1
      do i=-1,1
        ixs^L=ixO^L+i*ishift*kr(idims,^D);
        ixsR^L=ixs^L+ishift*kr(jdims,^D);
        ixsL^L=ixs^L-ishift*kr(jdims,^D);
        wMax(ixO^S)= max(wMax(ixO^S),w(ixsR^S),w(ixsL^S))
        wMin(ixO^S)= min(wMin(ixO^S),w(ixsR^S),w(ixsL^S))
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
          wMax(ixO^S)= max(wMax(ixO^S),w(ixsR^S),w(ixsL^S))
          wMin(ixO^S)= min(wMin(ixO^S),w(ixsR^S),w(ixsL^S))
        end do
      end do
      }
    end do

  end subroutine  extremaw

end module mod_ppm
