module mod_venk
  ! Venkatakrishnan limiter
  !
  ! 2019.10.11 coded up by nanami;
  !
  ! see Venkatakrishnan 1993 for this limiter;


  implicit none
  private

  public :: venklimiter

contains

  subroutine venklimiter(ixI^L,iL^L,idims,dxdim,w,wLp,wRp)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, iL^L, idims
    double precision, intent(in)    :: dxdim
    double precision, intent(in)    :: w(ixI^S)
    double precision, intent(inout) :: wRp(ixI^S),wLp(ixI^S) 
    !> local
    integer                         :: iM^L, iMm^L, iMp^L
    integer                         :: iLm^L, iLp^L, iLpp^L
    double precision                :: wmax(ixI^S),wmin(ixI^S)
    double precision                :: westp(ixI^S),westm(ixI^S)
    double precision                :: phi1(ixI^S),phi2(ixI^S)
    double precision                :: phi3(ixI^S),phi4(ixI^S)
    double precision                :: phim(ixI^S),phip(ixI^S),phi(ixI^S)
    double precision                :: deltap(ixI^S),deltam(ixI^S)
    double precision                :: eps2
    double precision, parameter     :: venk_omega = 1.d-12
    double precision, parameter     :: venk_k = 0.3d0

    iM^L=iL^L;
    iMmax^D=iLmax^D+kr(idims,^D);
    iLm^L=iL^L-kr(idims,^D);
    iMm^L=iM^L-kr(idims,^D);
    iLp^L=iL^L+kr(idims,^D);
    iMp^L=iM^L+kr(idims,^D);
    iLpp^L=iLp^L+kr(idims,^D);

    eps2 = (venk_k * dxdim) ** 3

    wmax(iM^S) = max(w(iMm^S), w(iM^S), w(iMp^S))
    wmin(iM^S) = min(w(iMm^S), w(iM^S), w(iMp^S))
    !> use central difference approximation as (eq.5) and take phi = 1
    westp(iM^S) = w(iM^S) + (w(iMp^S) - w(iMm^S)) * 0.25d0
    westm(iM^S) = w(iM^S) - (w(iMp^S) - w(iMm^S)) * 0.25d0

    !> (eq.30) & (eq.31)
    deltap = 0
    deltam = 0
    phi1 = 0
    phi2 = 0
    phip = 1
    where(westp(iM^S) .gt. w(iM^S))
      deltap(iM^S) = wmax(iM^S) - w(iM^S)
      deltam(iM^S) = westp(iM^S) - w(iM^S)
      deltam(iM^S) = sign(dabs(deltam(iM^S)) + venk_omega, deltam(iM^S))
      phi1(iM^S) = (deltap(iM^S)**2 + eps2) + 2.d0 * deltap(iM^S) * deltam(iM^S)
      phi2(iM^S) = deltap(iM^S)**2 + 2.d0 * deltam(iM^S)**2 + deltap(iM^S) * deltam(iM^S) + eps2
      phip(iM^S) = phi1(iM^S) / phi2(iM^S)
    elsewhere(westp(iM^S) .lt. w(iM^S))
      deltap(iM^S) = wmin(iM^S) - w(iM^S)
      deltam(iM^S) = westp(iM^S) - w(iM^S)
      deltam(iM^S) = sign(dabs(deltam(iM^S)) + venk_omega, deltam(iM^S))
      phi1(iM^S) = (deltap(iM^S)**2 + eps2) + 2.d0 * deltap(iM^S) * deltam(iM^S)
      phi2(iM^S) = deltap(iM^S)**2 + 2.d0 * deltam(iM^S)**2 + deltap(iM^S) * deltam(iM^S) + eps2
      phip(iM^S) = phi1(iM^S) / phi2(iM^S)
    elsewhere
      phip(iM^S) = one
   endwhere

    deltap = 0
    deltam = 0
    phi3 = 0
    phi4 = 0
    phim = 0
    where(westm(iM^S) .lt. w(iM^S))
      deltap(iM^S) = - (wmax(iM^S) - w(iM^S))
      deltam(iM^S) = westm(iM^S) - w(iM^S)
      deltam(iM^S) = sign(dabs(deltam(iM^S)) + venk_omega, deltam(iM^S))
      phi3(iM^S) = (deltap(iM^S)**2 + eps2) + 2.d0 * deltap(iM^S) * deltam(iM^S)
      phi4(iM^S) = deltap(iM^S)**2 + 2.d0 * deltam(iM^S)**2 + deltap(iM^S) * deltam(iM^S) + eps2
      phim(iM^S) = phi3(iM^S) / phi4(iM^S)
    elsewhere(westm(iM^S) .gt. w(iM^S))
      deltap(iM^S) = - (wmin(iM^S) - w(iM^S))
      deltam(iM^S) = westm(iM^S) - w(iM^S)
      deltam(iM^S) = sign(dabs(deltam(iM^S)) + venk_omega, deltam(iM^S))
      phi3(iM^S) = (deltap(iM^S)**2 + eps2) + 2.d0 * deltap(iM^S) * deltam(iM^S)
      phi4(iM^S) = deltap(iM^S)**2 + 2.d0 * deltam(iM^S)**2 + deltap(iM^S) * deltam(iM^S) + eps2
      phim(iM^S) = phi3(iM^S) / phi4(iM^S)
    elsewhere
      phim(iM^S) = one
    endwhere
    !> (eq.3)
    phi(iM^S) = min(phim(iM^S),phip(iM^S))
    !> (eq.5)
    wLp(iL^S) = w(iL^S) + 0.25d0 * (w(iLp^S)-w(iLm^S)) * phi(iL^S)
    wRp(iL^S) = w(iLp^S) - 0.25d0 * (w(iLpp^S)-w(iL^S)) * phi(iLp^S)

  end subroutine venklimiter

end module mod_venk
