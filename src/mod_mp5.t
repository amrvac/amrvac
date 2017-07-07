module mod_mp5

  implicit none
  private

  public :: MP5limiter
  public :: MP5limiterL
  public :: MP5limiterR

contains

  subroutine MP5limiter(ixI^L,iL^L,idims,w,wLC,wRC)
    ! MP5 limiter from Suresh & Huynh 1997
    ! Following the convention of Mignone et al. 2010.
    ! Needs at least three ghost cells.  Set nghostcells=3.

    use mod_global_parameters
    use mod_physics, only: phys_check_w

    integer, intent(in)             :: ixI^L, iL^L, idims
    double precision, intent(in)    :: w(ixI^S,1:nw)

    double precision, intent(inout) :: wRC(ixI^S,1:nw),wLC(ixI^S,1:nw) 
    ! .. local ..
    integer                         :: iLm^L, iLmm^L, iLp^L, iLpp^L, iLppp^L
    integer                         :: id^L, idp^L, idpp^L, idm^L, ie^L, iem^L, iep^L, iepp^L
    integer                         :: iw
    double precision, dimension(ixI^S,1:nw)  :: f, fmp, fmin, fmax, ful, dm4, d, fmd, flc, flim
    double precision, dimension(ixI^S,1:nw)  :: wRCtmp, wLCtmp
    double precision, dimension(ixI^S) :: tmp, tmp2, tmp3, a, b, c
    integer                         :: flagL(ixI^S), flagR(ixI^S)
    double precision, parameter     :: eps=0.d0, alpha=4.0d0
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

    iLm^L=iL^L-kr(idims,^D);
    iLmm^L=iLm^L-kr(idims,^D);
    iLp^L=iL^L+kr(idims,^D);
    iLpp^L=iLp^L+kr(idims,^D);

    f(iL^S,1:nwflux) = 1.0d0/60.0d0 * (&
         2.0d0* w(iLmm^S,1:nwflux) &
         - 13.0d0* w(iLm^S,1:nwflux) &
         + 47.0d0* w(iL^S,1:nwflux) &
         + 27.0d0* w(iLp^S,1:nwflux) &
         - 3.0d0*  w(iLpp^S,1:nwflux))

    ! get fmp and ful:
    do iw=1,nwflux
       a(iL^S) = w(iLp^S,iw)-w(iL^S,iw)
       b(iL^S) = alpha*(w(iL^S,iw)-w(iLm^S,iw))
       call minmod(ixI^L,iL^L,a,b,tmp)
       fmp(iL^S,iw) = w(iL^S,iw) + tmp(iL^S)
       ful(iL^S,iw) = w(iL^S,iw) + b(iL^S)
    end do ! iw loop

    ! get dm4:
    idmax^D=iLmax^D; idmin^D=iLmin^D-kr(idims,^D);
    idm^L=id^L-kr(idims,^D);
    idp^L=id^L+kr(idims,^D);

    iemax^D=idmax^D+kr(idims,^D); iemin^D=idmin^D;
    iem^L=ie^L-kr(idims,^D);
    iep^L=ie^L+kr(idims,^D);

    d(ie^S,1:nwflux) = w(iep^S,1:nwflux)-2.0d0*w(ie^S,1:nwflux)+w(iem^S,1:nwflux)

    do iw=1,nwflux
       a(id^S) = 4.0d0*d(id^S,iw)-d(idp^S,iw)
       b(id^S) = 4.0d0*d(idp^S,iw)-d(id^S,iw)
       call minmod(ixI^L,id^L,a,b,tmp)
       a(id^S) = d(id^S,iw)
       b(id^S) = d(idp^S,iw)
       call minmod(ixI^L,id^L,a,b,tmp2)
       call minmod(ixI^L,id^L,tmp,tmp2,tmp3)
       dm4(id^S,iw) = tmp3(id^S)
    end do

    ! get fmd:
    fmd(iL^S,1:nwflux) = (w(iL^S,1:nwflux)+w(iLp^S,1:nwflux))/2.0d0-dm4(iL^S,1:nwflux)/2.0d0

    !get flc: 
    flc(iL^S,1:nwflux) = half*(3.0d0*w(iL^S,1:nwflux) &
         - w(iLm^S,1:nwflux)) + 4.0d0/3.0d0*dm4(iLm^S,1:nwflux)

    fmin(iL^S,1:nwflux) = max(min(w(iL^S,1:nwflux),w(iLp^S,1:nwflux),fmd(iL^S,1:nwflux)),&
         min(w(iL^S,1:nwflux),ful(iL^S,1:nwflux),flc(iL^S,1:nwflux)))

    fmax(iL^S,1:nwflux) = min(max(w(iL^S,1:nwflux),w(iLp^S,1:nwflux),fmd(iL^S,1:nwflux)),&
         max(w(iL^S,1:nwflux),ful(iL^S,1:nwflux),flc(iL^S,1:nwflux)))

    do iw=1,nwflux
       a(iL^S) = fmin(iL^S,iw)
       b(iL^S) = f(iL^S,iw)
       c(iL^S) = fmax(iL^S,iw)
       call median(ixI^L,iL^L,a,b,c,tmp)
       flim(iL^S,iw) = tmp(iL^S)
    end do



    ! check case
    where ((f(iL^S,1:nwflux)-w(iL^S,1:nwflux))*(f(iL^S,1:nwflux)-fmp(iL^S,1:nwflux)) .le. eps)
       wLCtmp(iL^S,1:nwflux) = f(iL^S,1:nwflux)
    elsewhere
       wLCtmp(iL^S,1:nwflux) = flim(iL^S,1:nwflux)
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

    iLppp^L=iLpp^L+kr(idims,^D);

    f(iL^S,1:nwflux) = 1.0d0/60.0d0 * (&
         2.0d0* w(iLppp^S,1:nwflux) &
         - 13.0d0* w(iLpp^S,1:nwflux) &
         + 47.0d0* w(iLp^S,1:nwflux) &
         + 27.0d0* w(iL^S,1:nwflux) &
         - 3.0d0*  w(iLm^S,1:nwflux))

    ! get fmp and ful:
    do iw=1,nwflux
       a(iL^S) = w(iL^S,iw)-w(iLp^S,iw)
       b(iL^S) = alpha*(w(iLp^S,iw)-w(iLpp^S,iw))
       call minmod(ixI^L,iL^L,a,b,tmp)
       fmp(iL^S,iw) = w(iLp^S,iw) + tmp(iL^S)
       ful(iL^S,iw) = w(iLp^S,iw) + b(iL^S)
    end do ! iw loop

    ! get dm4:
    idmax^D=iLmax^D+kr(idims,^D); idmin^D=iLmin^D;
    idm^L=id^L-kr(idims,^D);
    idp^L=id^L+kr(idims,^D);

    iemax^D=idmax^D; iemin^D=idmin^D-kr(idims,^D);
    iem^L=ie^L-kr(idims,^D);
    iep^L=ie^L+kr(idims,^D);
    iepp^L=iep^L+kr(idims,^D);

    d(ie^S,1:nwflux) = w(ie^S,1:nwflux)-2.0d0*w(iep^S,1:nwflux)+w(iepp^S,1:nwflux)

    do iw=1,nwflux
       a(id^S) = 4.0d0*d(id^S,iw)-d(idm^S,iw)
       b(id^S) = 4.0d0*d(idm^S,iw)-d(id^S,iw)
       call minmod(ixI^L,id^L,a,b,tmp)
       a(id^S) = d(id^S,iw)
       b(id^S) = d(idm^S,iw)
       call minmod(ixI^L,id^L,a,b,tmp2)
       call minmod(ixI^L,id^L,tmp,tmp2,tmp3)
       dm4(id^S,iw) = tmp3(id^S)
    end do

    ! get fmd:
    fmd(iL^S,1:nwflux) = (w(iL^S,1:nwflux)+w(iLp^S,1:nwflux))/2.0d0-dm4(iL^S,1:nwflux)/2.0d0

    !get flc: 
    flc(iL^S,1:nwflux) = half*(3.0d0*w(iLp^S,1:nwflux) &
         - w(iLpp^S,1:nwflux)) + 4.0d0/3.0d0*dm4(iLp^S,1:nwflux)

    fmin(iL^S,1:nwflux) = max(min(w(iLp^S,1:nwflux),w(iL^S,1:nwflux),fmd(iL^S,1:nwflux)),&
         min(w(iLp^S,1:nwflux),ful(iL^S,1:nwflux),flc(iL^S,1:nwflux)))

    fmax(iL^S,1:nwflux) = min(max(w(iLp^S,1:nwflux),w(iL^S,1:nwflux),fmd(iL^S,1:nwflux)),&
         max(w(iLp^S,1:nwflux),ful(iL^S,1:nwflux),flc(iL^S,1:nwflux)))

    do iw=1,nwflux
       a(iL^S) = fmin(iL^S,iw)
       b(iL^S) = f(iL^S,iw)
       c(iL^S) = fmax(iL^S,iw)
       call median(ixI^L,iL^L,a,b,c,tmp)
       flim(iL^S,iw) = tmp(iL^S)
    end do

    ! check case
    where ((f(iL^S,1:nwflux)-w(iLp^S,1:nwflux))*(f(iL^S,1:nwflux)-fmp(iL^S,1:nwflux))  .le. eps)
       wRCtmp(iL^S,1:nwflux) = f(iL^S,1:nwflux)
    elsewhere
       wRCtmp(iL^S,1:nwflux) = flim(iL^S,1:nwflux)
    end where

    ! Since limiter not TVD, negative pressures or densities could result.  
    ! Fall back to flat interpolation (minmod would also work). 
    call phys_check_w(.true.,ixG^LL,iL^L,wLCtmp,flagL)
    call phys_check_w(.true.,ixG^LL,iL^L,wRCtmp,flagR)

    do iw=1,nwflux
       where (flagL(iL^S) == 0 .and. flagR(iL^S) == 0)
          wLC(iL^S,iw)=wLCtmp(iL^S,iw)
          wRC(iL^S,iw)=wRCtmp(iL^S,iw)
       end where
    end do


  end subroutine MP5limiter
  !============================================================================
  subroutine MP5limiterL(ixI^L,iL^L,idims,w,wLC)
    ! MP5 limiter from Suresh & Huynh 1997
    ! Following the convention of Mignone et al. 2010.
    ! Needs at least three ghost cells.  Set nghostcells=3.

    use mod_global_parameters

    integer, intent(in)             :: ixI^L, iL^L, idims
    double precision, intent(in)    :: w(ixI^S,1:nw)

    double precision, intent(inout) :: wLC(ixI^S,1:nw) 
    ! .. local ..
    integer                         :: iLm^L, iLmm^L, iLp^L, iLpp^L
    integer                         :: id^L, idp^L, idpp^L, idm^L, ie^L, iem^L, iep^L, iepp^L
    integer                         :: iw
    double precision, dimension(ixI^S,1:nw)  :: f, fmp, fmin, fmax, ful, dm4, d, fmd, flc, flim
    double precision, dimension(ixI^S) :: tmp, tmp2, tmp3, a, b, c
    double precision, parameter     :: eps=0.d0, alpha=4.0d0
    !double precision                :: alpha
    !----------------------------------------------------------------------------

    ! Variable alpha:
    !alpha = float(nstep)/courantpar - one

    ! Left side:


    iLm^L=iL^L-kr(idims,^D);
    iLmm^L=iLm^L-kr(idims,^D);
    iLp^L=iL^L+kr(idims,^D);
    iLpp^L=iLp^L+kr(idims,^D);

    f(iL^S,1:nwflux) = 1.0d0/60.0d0 * (&
         2.0d0* w(iLmm^S,1:nwflux) &
         - 13.0d0* w(iLm^S,1:nwflux) &
         + 47.0d0* w(iL^S,1:nwflux) &
         + 27.0d0* w(iLp^S,1:nwflux) &
         - 3.0d0*  w(iLpp^S,1:nwflux))

    ! get fmp and ful:
    do iw=1,nwflux
       a(iL^S) = w(iLp^S,iw)-w(iL^S,iw)
       b(iL^S) = alpha*(w(iL^S,iw)-w(iLm^S,iw))
       call minmod(ixI^L,iL^L,a,b,tmp)
       fmp(iL^S,iw) = w(iL^S,iw) + tmp(iL^S)
       ful(iL^S,iw) = w(iL^S,iw) + b(iL^S)
    end do ! iw loop

    ! get dm4:
    idmax^D=iLmax^D; idmin^D=iLmin^D-kr(idims,^D);
    idm^L=id^L-kr(idims,^D);
    idp^L=id^L+kr(idims,^D);

    iemax^D=idmax^D+kr(idims,^D); iemin^D=idmin^D;
    iem^L=ie^L-kr(idims,^D);
    iep^L=ie^L+kr(idims,^D);

    d(ie^S,1:nwflux) = w(iep^S,1:nwflux)-2.0d0*w(ie^S,1:nwflux)+w(iem^S,1:nwflux)

    do iw=1,nwflux
       a(id^S) = 4.0d0*d(id^S,iw)-d(idp^S,iw)
       b(id^S) = 4.0d0*d(idp^S,iw)-d(id^S,iw)
       call minmod(ixI^L,id^L,a,b,tmp)
       a(id^S) = d(id^S,iw)
       b(id^S) = d(idp^S,iw)
       call minmod(ixI^L,id^L,a,b,tmp2)
       call minmod(ixI^L,id^L,tmp,tmp2,tmp3)
       dm4(id^S,iw) = tmp3(id^S)
    end do

    ! get fmd:
    fmd(iL^S,1:nwflux) = (w(iL^S,1:nwflux)+w(iLp^S,1:nwflux))/2.0d0-dm4(iL^S,1:nwflux)/2.0d0

    !get flc: 
    flc(iL^S,1:nwflux) = half*(3.0d0*w(iL^S,1:nwflux) &
         - w(iLm^S,1:nwflux)) + 4.0d0/3.0d0*dm4(iLm^S,1:nwflux)

    fmin(iL^S,1:nwflux) = max(min(w(iL^S,1:nwflux),w(iLp^S,1:nwflux),fmd(iL^S,1:nwflux)),&
         min(w(iL^S,1:nwflux),ful(iL^S,1:nwflux),flc(iL^S,1:nwflux)))

    fmax(iL^S,1:nwflux) = min(max(w(iL^S,1:nwflux),w(iLp^S,1:nwflux),fmd(iL^S,1:nwflux)),&
         max(w(iL^S,1:nwflux),ful(iL^S,1:nwflux),flc(iL^S,1:nwflux)))

    do iw=1,nwflux
       a(iL^S) = fmin(iL^S,iw)
       b(iL^S) = f(iL^S,iw)
       c(iL^S) = fmax(iL^S,iw)
       call median(ixI^L,iL^L,a,b,c,tmp)
       flim(iL^S,iw) = tmp(iL^S)
    end do


    ! check case
    where ((f(iL^S,1:nwflux)-w(iL^S,1:nwflux))*(f(iL^S,1:nwflux)-fmp(iL^S,1:nwflux)) .le. eps)
       wLC(iL^S,1:nwflux) = f(iL^S,1:nwflux)
    elsewhere
       wLC(iL^S,1:nwflux) = flim(iL^S,1:nwflux)
    end where


  end subroutine MP5limiterL
  !============================================================================
  subroutine MP5limiterR(ixI^L,iL^L,idims,w,wRC)
    ! MP5 limiter from Suresh & Huynh 1997
    ! Following the convention of Mignone et al. 2010.
    ! Needs at least three ghost cells.  Set nghostcells=3.

    use mod_global_parameters

    integer, intent(in)             :: ixI^L, iL^L, idims
    double precision, intent(in)    :: w(ixI^S,1:nw)

    double precision, intent(inout) :: wRC(ixI^S,1:nw)
    ! .. local ..
    integer                         :: iLm^L, iLp^L, iLpp^L, iLppp^L
    integer                         :: id^L, idp^L, idpp^L, idm^L, ie^L, iem^L, iep^L, iepp^L
    integer                         :: iw
    double precision, dimension(ixI^S,1:nw)  :: f, fmp, fmin, fmax, ful, dm4, d, fmd, flc, flim
    double precision, dimension(ixI^S) :: tmp, tmp2, tmp3, a, b, c
    double precision, parameter     :: eps=0.d0, alpha=4.0d0
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

    iLm^L=iL^L-kr(idims,^D);
    iLp^L=iL^L+kr(idims,^D);
    iLpp^L=iLp^L+kr(idims,^D);
    iLppp^L=iLpp^L+kr(idims,^D);

    f(iL^S,1:nwflux) = 1.0d0/60.0d0 * (&
         2.0d0* w(iLppp^S,1:nwflux) &
         - 13.0d0* w(iLpp^S,1:nwflux) &
         + 47.0d0* w(iLp^S,1:nwflux) &
         + 27.0d0* w(iL^S,1:nwflux) &
         - 3.0d0*  w(iLm^S,1:nwflux))

    ! get fmp and ful:
    do iw=1,nwflux
       a(iL^S) = w(iL^S,iw)-w(iLp^S,iw)
       b(iL^S) = alpha*(w(iLp^S,iw)-w(iLpp^S,iw))
       call minmod(ixI^L,iL^L,a,b,tmp)
       fmp(iL^S,iw) = w(iLp^S,iw) + tmp(iL^S)
       ful(iL^S,iw) = w(iLp^S,iw) + b(iL^S)
    end do ! iw loop

    ! get dm4:
    idmax^D=iLmax^D+kr(idims,^D); idmin^D=iLmin^D;
    idm^L=id^L-kr(idims,^D);
    idp^L=id^L+kr(idims,^D);

    iemax^D=idmax^D; iemin^D=idmin^D-kr(idims,^D);
    iem^L=ie^L-kr(idims,^D);
    iep^L=ie^L+kr(idims,^D);
    iepp^L=iep^L+kr(idims,^D);

    d(ie^S,1:nwflux) = w(ie^S,1:nwflux)-2.0d0*w(iep^S,1:nwflux)+w(iepp^S,1:nwflux)

    do iw=1,nwflux
       a(id^S) = 4.0d0*d(id^S,iw)-d(idm^S,iw)
       b(id^S) = 4.0d0*d(idm^S,iw)-d(id^S,iw)
       call minmod(ixI^L,id^L,a,b,tmp)
       a(id^S) = d(id^S,iw)
       b(id^S) = d(idm^S,iw)
       call minmod(ixI^L,id^L,a,b,tmp2)
       call minmod(ixI^L,id^L,tmp,tmp2,tmp3)
       dm4(id^S,iw) = tmp3(id^S)
    end do

    ! get fmd:
    fmd(iL^S,1:nwflux) = (w(iL^S,1:nwflux)+w(iLp^S,1:nwflux))/2.0d0-dm4(iL^S,1:nwflux)/2.0d0

    !get flc: 
    flc(iL^S,1:nwflux) = half*(3.0d0*w(iLp^S,1:nwflux) &
         - w(iLpp^S,1:nwflux)) + 4.0d0/3.0d0*dm4(iLp^S,1:nwflux)

    fmin(iL^S,1:nwflux) = max(min(w(iLp^S,1:nwflux),w(iL^S,1:nwflux),fmd(iL^S,1:nwflux)),&
         min(w(iLp^S,1:nwflux),ful(iL^S,1:nwflux),flc(iL^S,1:nwflux)))

    fmax(iL^S,1:nwflux) = min(max(w(iLp^S,1:nwflux),w(iL^S,1:nwflux),fmd(iL^S,1:nwflux)),&
         max(w(iLp^S,1:nwflux),ful(iL^S,1:nwflux),flc(iL^S,1:nwflux)))

    do iw=1,nwflux
       a(iL^S) = fmin(iL^S,iw)
       b(iL^S) = f(iL^S,iw)
       c(iL^S) = fmax(iL^S,iw)
       call median(ixI^L,iL^L,a,b,c,tmp)
       flim(iL^S,iw) = tmp(iL^S)
    end do

    ! check case
    where ((f(iL^S,1:nwflux)-w(iLp^S,1:nwflux))*(f(iL^S,1:nwflux)-fmp(iL^S,1:nwflux))  .le. eps)
       wRC(iL^S,1:nwflux) = f(iL^S,1:nwflux)
    elsewhere
       wRC(iL^S,1:nwflux) = flim(iL^S,1:nwflux)
    end where

  end subroutine Mp5limiterR
  !============================================================================
  subroutine minmod(ixI^L,ixO^L,a,b,minm)

    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: a(ixI^S), b(ixI^S)
    double precision, intent(out):: minm(ixI^S)

    minm(ixO^S) = (sign(one,a(ixO^S))+sign(one,b(ixO^S)))/2.0d0 * &
         min(abs(a(ixO^S)),abs(b(ixO^S)))

  end subroutine minmod
  !============================================================================
  subroutine median(ixI^L,ixO^L,a,b,c,med)

    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: a(ixI^S), b(ixI^S), c(ixI^S)
    double precision, intent(out):: med(ixI^S)
    ! .. local ..
    double precision             :: tmp1(ixI^S),tmp2(ixI^S)

    tmp1(ixO^S) = b(ixO^S) - a(ixO^S); tmp2(ixO^S) = c(ixO^S) - a(ixO^S)

    med(ixO^S) = a(ixO^S) + (sign(one,tmp1(ixO^S))+sign(one,tmp2(ixO^S)))/2.0d0 * &
         min(abs(tmp1(ixO^S)),abs(tmp2(ixO^S)))

  end subroutine median

end module mod_mp5
