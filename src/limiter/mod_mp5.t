!> Module containing the MP5 (fifth order) flux scheme
module mod_mp5

  implicit none
  private

  double precision, parameter :: mp5_eps_machine = 1.0d-18
  double precision, parameter :: mp5_alpha       = 4.0d0

  public :: MP5limiter
  public :: MP5limiterL
  public :: MP5limiterR

contains

  !> MP5 limiter from Suresh & Huynh 1997 
  !> Following the convention of Mignone et al. 2010. 
  !> Needs at least three ghost cells
  subroutine MP5limiter(ixI^L,iL^L,idims,w,wLC,wRC,rec)
    use mod_global_parameters

    integer, intent(in)                :: ixI^L, iL^L, idims
    double precision, intent(in)       :: w(ixI^S)
    double precision, intent(inout)    :: wRC(ixI^S),wLC(ixI^S) 
    logical, intent(in)                :: rec
    ! .. local ..
    integer                            :: iLm^L, iLmm^L, iLp^L, iLpp^L, iLppp^L
    integer                            :: id^L, idp^L, idpp^L, idm^L, ie^L, iem^L, iep^L, iepp^L
    integer                            :: flagL(ixI^S), flagR(ixI^S)
    double precision, dimension(ixI^S) :: f, fmp, fmin, fmax, ful, dm4, d, fmd, flc, flim
    double precision, dimension(ixI^S) :: tmp, tmp2, tmp3, a, b, c
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

    if (rec) then
      f(iL^S) = 1.0d0/60.0d0 * (  &
             2.0d0* w(iLmm^S)     &
           - 13.0d0* w(iLm^S)     &
           + 47.0d0* w(iL^S)      &
           + 27.0d0* w(iLp^S)     &
           - 3.0d0*  w(iLpp^S))
    else
      f(iL^S) = 1.0d0/128.0d0 * ( &
             3.0d0* w(iLmm^S)     &
           - 20.0d0* w(iLm^S)     &
           + 90.0d0* w(iL^S)      &
           + 60.0d0* w(iLp^S)     &
           - 5.0d0*  w(iLpp^S))
    end if

    ! get fmp and ful:
    a(iL^S) = w(iLp^S)-w(iL^S)
    b(iL^S) = mp5_alpha*(w(iL^S)-w(iLm^S))
    call minmod(ixI^L,iL^L,a,b,tmp)
    fmp(iL^S) = w(iL^S) + tmp(iL^S)
    ful(iL^S) = w(iL^S) + b(iL^S)

    ! get dm4:
    idmax^D=iLmax^D; idmin^D=iLmin^D-kr(idims,^D);
    idm^L=id^L-kr(idims,^D);
    idp^L=id^L+kr(idims,^D);

    iemax^D=idmax^D+kr(idims,^D); iemin^D=idmin^D;
    iem^L=ie^L-kr(idims,^D);
    iep^L=ie^L+kr(idims,^D);

    d(ie^S) = w(iep^S)-2.0d0*w(ie^S)+w(iem^S)

    a(id^S) = 4.0d0*d(id^S)-d(idp^S)
    b(id^S) = 4.0d0*d(idp^S)-d(id^S)
    call minmod(ixI^L,id^L,a,b,tmp)
    a(id^S) = d(id^S)
    b(id^S) = d(idp^S)
    call minmod(ixI^L,id^L,a,b,tmp2)
    call minmod(ixI^L,id^L,tmp,tmp2,tmp3)
    dm4(id^S) = tmp3(id^S)

    ! get fmd:
    fmd(iL^S) = (w(iL^S)+w(iLp^S))/2.0d0-dm4(iL^S)/2.0d0

    !get flc: 
    flc(iL^S) = half*(3.0d0*w(iL^S) &
         - w(iLm^S)) + 4.0d0/3.0d0*dm4(iLm^S)

    fmin(iL^S) = max(min(w(iL^S),w(iLp^S),fmd(iL^S)),&
         min(w(iL^S),ful(iL^S),flc(iL^S)))

    fmax(iL^S) = min(max(w(iL^S),w(iLp^S),fmd(iL^S)),&
         max(w(iL^S),ful(iL^S),flc(iL^S)))

    call median(ixI^L,iL^L,fmin,f,fmax,flim)

    ! check case
    where ((f(iL^S)-w(iL^S))*(f(iL^S)-fmp(iL^S)) <= mp5_eps_machine)
       wLC(iL^S) = f(iL^S)
    elsewhere
       wLC(iL^S) = flim(iL^S)
    end where

    ! Right side:
    ! the interpolation from the right is obtained when the left-hand fromula is applied to
    ! data mirrored about the interface.  
    ! thus substitute: 
    ! i-2 -> i+3
    ! i-1 -> i+2
    ! i   -> i+1
    ! i+1 -> i
    ! i+2 -> i-1

    iLppp^L=iLpp^L+kr(idims,^D);

    if (rec) then
      f(iL^S) = 1.0d0/60.0d0 * (  &
              2.0d0* w(iLppp^S)   &
           - 13.0d0* w(iLpp^S)    &
           + 47.0d0* w(iLp^S)     &
           + 27.0d0* w(iL^S)      &
           - 3.0d0*  w(iLm^S))
    else
      f(iL^S) = 1.0d0/128.0d0 * ( &
              3.0d0* w(iLppp^S)   &
           - 20.0d0* w(iLpp^S)    &
           + 90.0d0* w(iLp^S)     &
           + 60.0d0* w(iL^S)      &
           - 5.0d0*  w(iLm^S))
    end if

    ! get fmp and ful:
    a(iL^S) = w(iL^S)-w(iLp^S)
    b(iL^S) = mp5_alpha*(w(iLp^S)-w(iLpp^S))
    call minmod(ixI^L,iL^L,a,b,tmp)
    fmp(iL^S) = w(iLp^S) + tmp(iL^S)
    ful(iL^S) = w(iLp^S) + b(iL^S)

    ! get dm4:
    idmax^D=iLmax^D+kr(idims,^D); idmin^D=iLmin^D;
    idm^L=id^L-kr(idims,^D);
    idp^L=id^L+kr(idims,^D);

    iemax^D=idmax^D; iemin^D=idmin^D-kr(idims,^D);
    iem^L=ie^L-kr(idims,^D);
    iep^L=ie^L+kr(idims,^D);
    iepp^L=iep^L+kr(idims,^D);

    d(ie^S) = w(ie^S)-2.0d0*w(iep^S)+w(iepp^S)

    a(id^S) = 4.0d0*d(id^S)-d(idm^S)
    b(id^S) = 4.0d0*d(idm^S)-d(id^S)
    call minmod(ixI^L,id^L,a,b,tmp)
    a(id^S) = d(id^S)
    b(id^S) = d(idm^S)
    call minmod(ixI^L,id^L,a,b,tmp2)
    call minmod(ixI^L,id^L,tmp,tmp2,tmp3)
    dm4(id^S) = tmp3(id^S)

    ! get fmd:
    fmd(iL^S) = (w(iL^S)+w(iLp^S))/2.0d0-dm4(iL^S)/2.0d0

    !get flc: 
    flc(iL^S) = half*(3.0d0*w(iLp^S) &
         - w(iLpp^S)) + 4.0d0/3.0d0*dm4(iLp^S)

    fmin(iL^S) = max(min(w(iLp^S),w(iL^S),fmd(iL^S)),&
         min(w(iLp^S),ful(iL^S),flc(iL^S)))

    fmax(iL^S) = min(max(w(iLp^S),w(iL^S),fmd(iL^S)),&
         max(w(iLp^S),ful(iL^S),flc(iL^S)))

    call median(ixI^L,iL^L,fmin,f,fmax,flim)

    ! check case
    where ((f(iL^S)-w(iLp^S))*(f(iL^S)-fmp(iL^S)) <= mp5_eps_machine)
       wRC(iL^S) = f(iL^S)
    elsewhere
       wRC(iL^S) = flim(iL^S)
    end where

    ! fixme: refine here
    ! Since limiter not TVD, negative pressures or densities could result.  
    ! Fall back to flat interpolation (minmod would also work). 
    !call phys_check_prim(ixG^LL,iL^L,wLCtmp,flagL)
    !call phys_check_prim(ixG^LL,iL^L,wRCtmp,flagR)

    !where (flagL(iL^S) == 0 .and. flagR(iL^S) == 0)
    !   wLC(iL^S)=wLCtmp(iL^S)
    !   wRC(iL^S)=wRCtmp(iL^S)
    !end where
  end subroutine MP5limiter

  subroutine MP5limiterL(rec_from,rec_to,n_w,ixI^L,iL^L,idims,w,wLC,rec)

    use mod_global_parameters
    use mod_physics, only: phys_check_prim

    integer, intent(in)             :: rec_from, rec_to, n_w
    integer, intent(in)             :: ixI^L, iL^L, idims
    double precision, intent(in)    :: w(ixI^S,1:n_w)
    double precision, intent(inout) :: wLC(ixI^S,1:n_w) 
    logical, intent(in)             :: rec

    ! .. local ..
    integer                         :: iLm^L, iLmm^L, iLp^L, iLpp^L
    integer                         :: id^L, idp^L, idpp^L, idm^L, ie^L, iem^L, iep^L, iepp^L
    integer                         :: iw
    double precision, dimension(ixI^S,1:n_w)  :: f, fmp, fmin, fmax, ful, dm4, d, fmd, flc, flim
    double precision, dimension(ixI^S) :: tmp, tmp2, tmp3, a, b, c

    ! Variable alpha:
    !alpha = float(nstep)/courantpar - one

    ! Left side:

    iLm^L=iL^L-kr(idims,^D);
    iLmm^L=iLm^L-kr(idims,^D);
    iLp^L=iL^L+kr(idims,^D);
    iLpp^L=iLp^L+kr(idims,^D);

    if (rec) then
      f(iL^S,rec_from:rec_to) = 1.0d0/60.0d0 * ( &
              2.0d0* w(iLmm^S,rec_from:rec_to)   &
           - 13.0d0* w(iLm^S,rec_from:rec_to)    &
           + 47.0d0* w(iL^S,rec_from:rec_to)     &
           + 27.0d0* w(iLp^S,rec_from:rec_to)    &
           - 3.0d0*  w(iLpp^S,rec_from:rec_to))
    else
      f(iL^S,rec_from:rec_to) = 1.0d0/128.0d0 * ( &
              3.0d0* w(iLmm^S,rec_from:rec_to)    &
           - 20.0d0* w(iLm^S,rec_from:rec_to)     &
           + 90.0d0* w(iL^S,rec_from:rec_to)      &
           + 60.0d0* w(iLp^S,rec_from:rec_to)     &
           - 5.0d0*  w(iLpp^S,rec_from:rec_to))
    end if

    ! get fmp and ful:
    do iw=rec_from,rec_to
       a(iL^S) = w(iLp^S,iw)-w(iL^S,iw)
       b(iL^S) = mp5_alpha*(w(iL^S,iw)-w(iLm^S,iw))
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

    d(ie^S,rec_from:rec_to) = w(iep^S,rec_from:rec_to)-2.0d0*w(ie^S,rec_from:rec_to)+w(iem^S,rec_from:rec_to)

    do iw=rec_from,rec_to
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
    fmd(iL^S,rec_from:rec_to) = (w(iL^S,rec_from:rec_to)+w(iLp^S,rec_from:rec_to))/2.0d0-dm4(iL^S,rec_from:rec_to)/2.0d0

    !get flc: 
    flc(iL^S,rec_from:rec_to) = half*(3.0d0*w(iL^S,rec_from:rec_to) &
         - w(iLm^S,rec_from:rec_to)) + 4.0d0/3.0d0*dm4(iLm^S,rec_from:rec_to)

    fmin(iL^S,rec_from:rec_to) = max(min(w(iL^S,rec_from:rec_to),w(iLp^S,rec_from:rec_to),fmd(iL^S,rec_from:rec_to)),&
         min(w(iL^S,rec_from:rec_to),ful(iL^S,rec_from:rec_to),flc(iL^S,rec_from:rec_to)))

    fmax(iL^S,rec_from:rec_to) = min(max(w(iL^S,rec_from:rec_to),w(iLp^S,rec_from:rec_to),fmd(iL^S,rec_from:rec_to)),&
         max(w(iL^S,rec_from:rec_to),ful(iL^S,rec_from:rec_to),flc(iL^S,rec_from:rec_to)))

    do iw=rec_from,rec_to
       a(iL^S) = fmin(iL^S,iw)
       b(iL^S) = f(iL^S,iw)
       c(iL^S) = fmax(iL^S,iw)
       call median(ixI^L,iL^L,a,b,c,tmp)
       flim(iL^S,iw) = tmp(iL^S)
    end do

    ! check case
    where ((f(iL^S,rec_from:rec_to)-w(iL^S,rec_from:rec_to))*(f(iL^S,rec_from:rec_to)-fmp(iL^S,rec_from:rec_to)) <= mp5_eps_machine)
       wLC(iL^S,rec_from:rec_to) = f(iL^S,rec_from:rec_to)
    elsewhere
       wLC(iL^S,rec_from:rec_to) = flim(iL^S,rec_from:rec_to)
    end where

  end subroutine MP5limiterL

  subroutine MP5limiterR(rec_from,rec_to,n_w,ixI^L,iL^L,idims,w,wRC,rec)
    use mod_global_parameters

    integer, intent(in)             :: rec_from, rec_to, n_w
    integer, intent(in)             :: ixI^L, iL^L, idims
    double precision, intent(in)    :: w(ixI^S,1:n_w)
    double precision, intent(inout) :: wRC(ixI^S,1:n_w)
    logical, intent(in)             :: rec
    ! .. local ..
    integer                         :: iLm^L, iLp^L, iLpp^L, iLppp^L
    integer                         :: id^L, idp^L, idpp^L, idm^L, ie^L, iem^L, iep^L, iepp^L
    integer                         :: iw
    double precision, dimension(ixI^S,1:n_w)  :: f, fmp, fmin, fmax, ful, dm4, d, fmd, flc, flim
    double precision, dimension(ixI^S) :: tmp, tmp2, tmp3, a, b, c
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

    if (rec) then
      f(iL^S,rec_from:rec_to) = 1.0d0/60.0d0 * ( &
              2.0d0* w(iLppp^S,rec_from:rec_to)  &
           - 13.0d0* w(iLpp^S,rec_from:rec_to)   &
           + 47.0d0* w(iLp^S,rec_from:rec_to)    &
           + 27.0d0* w(iL^S,rec_from:rec_to)     &
           - 3.0d0*  w(iLm^S,rec_from:rec_to))
    else
      f(iL^S,rec_from:rec_to) = 1.0d0/128.0d0 * ( &
              3.0d0* w(iLppp^S,rec_from:rec_to)   &
           - 20.0d0* w(iLpp^S,rec_from:rec_to)    &
           + 90.0d0* w(iLp^S,rec_from:rec_to)     &
           + 60.0d0* w(iL^S,rec_from:rec_to)      &
           - 5.0d0*  w(iLm^S,rec_from:rec_to))
    end if

    ! get fmp and ful:
    do iw=rec_from,rec_to
       a(iL^S) = w(iL^S,iw)-w(iLp^S,iw)
       b(iL^S) = mp5_alpha*(w(iLp^S,iw)-w(iLpp^S,iw))
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

    d(ie^S,rec_from:rec_to) = w(ie^S,rec_from:rec_to)-2.0d0*w(iep^S,rec_from:rec_to)+w(iepp^S,rec_from:rec_to)

    do iw=rec_from,rec_to
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
    fmd(iL^S,rec_from:rec_to) = (w(iL^S,rec_from:rec_to)+w(iLp^S,rec_from:rec_to))/2.0d0-dm4(iL^S,rec_from:rec_to)/2.0d0

    !get flc: 
    flc(iL^S,rec_from:rec_to) = half*(3.0d0*w(iLp^S,rec_from:rec_to) &
         - w(iLpp^S,rec_from:rec_to)) + 4.0d0/3.0d0*dm4(iLp^S,rec_from:rec_to)

    fmin(iL^S,rec_from:rec_to) = max(min(w(iLp^S,rec_from:rec_to),w(iL^S,rec_from:rec_to),fmd(iL^S,rec_from:rec_to)),&
         min(w(iLp^S,rec_from:rec_to),ful(iL^S,rec_from:rec_to),flc(iL^S,rec_from:rec_to)))

    fmax(iL^S,rec_from:rec_to) = min(max(w(iLp^S,rec_from:rec_to),w(iL^S,rec_from:rec_to),fmd(iL^S,rec_from:rec_to)),&
         max(w(iLp^S,rec_from:rec_to),ful(iL^S,rec_from:rec_to),flc(iL^S,rec_from:rec_to)))

    do iw=rec_from,rec_to
       a(iL^S) = fmin(iL^S,iw)
       b(iL^S) = f(iL^S,iw)
       c(iL^S) = fmax(iL^S,iw)
       call median(ixI^L,iL^L,a,b,c,tmp)
       flim(iL^S,iw) = tmp(iL^S)
    end do

    ! check case
    where ((f(iL^S,rec_from:rec_to)-w(iLp^S,rec_from:rec_to))*(f(iL^S,rec_from:rec_to)-fmp(iL^S,rec_from:rec_to)) <= mp5_eps_machine)
       wRC(iL^S,rec_from:rec_to) = f(iL^S,rec_from:rec_to)
    elsewhere
       wRC(iL^S,rec_from:rec_to) = flim(iL^S,rec_from:rec_to)
    end where

  end subroutine Mp5limiterR

  subroutine minmod(ixI^L,ixO^L,a,b,minm)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: a(ixI^S), b(ixI^S)
    double precision, intent(out):: minm(ixI^S)
    minm(ixO^S) = (sign(one,a(ixO^S))+sign(one,b(ixO^S)))/2.0d0 * &
         min(abs(a(ixO^S)),abs(b(ixO^S)))
  end subroutine minmod

  subroutine median(ixI^L,ixO^L,a,b,c,med)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: a(ixI^S), b(ixI^S), c(ixI^S)
    double precision, intent(out):: med(ixI^S)
    double precision             :: tmp1(ixI^S),tmp2(ixI^S)
    tmp1(ixO^S) = b(ixO^S) - a(ixO^S); tmp2(ixO^S) = c(ixO^S) - a(ixO^S)
    med(ixO^S) = a(ixO^S) + (sign(one,tmp1(ixO^S))+sign(one,tmp2(ixO^S)))/2.0d0 * &
         min(abs(tmp1(ixO^S)),abs(tmp2(ixO^S)))
  end subroutine median

end module mod_mp5
