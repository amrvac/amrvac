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

  subroutine venklimiter(ixI^L,iL^L,idims,dxdim,w,wLC,wRC)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, iL^L, idims
    double precision, intent(in)    :: dxdim
    double precision, intent(in)    :: w(ixI^S,1:nw)
    double precision, intent(inout) :: wRC(ixI^S,1:nw),wLC(ixI^S,1:nw) 
    !> local
    integer                         :: iM^L, iMm^L, iMp^L
    integer                         :: iLm^L, iLp^L, iLpp^L
    double precision                :: wmax(ixI^S,1:nw),wmin(ixI^S,1:nw)
    double precision                :: westp(ixI^S,1:nw),westm(ixI^S,1:nw)
    double precision                :: phi1(ixI^S,1:nw),phi2(ixI^S,1:nw)
    double precision                :: phi3(ixI^S,1:nw),phi4(ixI^S,1:nw)
    double precision                :: phim(ixI^S,1:nw),phip(ixI^S,1:nw),phi(ixI^S,1:nw)
    double precision                :: deltap(ixI^S,1:nw),deltam(ixI^S,1:nw)
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

    wmax(iM^S,1:nwflux) = max(w(iMm^S,1:nwflux), w(iM^S,1:nwflux), w(iMp^S,1:nwflux))
    wmin(iM^S,1:nwflux) = min(w(iMm^S,1:nwflux), w(iM^S,1:nwflux), w(iMp^S,1:nwflux))
    !> use central difference approximation as (eq.5) and take phi = 1
    westp(iM^S,1:nwflux) = w(iM^S,1:nwflux) + (w(iMp^S,1:nwflux) - w(iMm^S,1:nwflux)) * 0.25d0
    westm(iM^S,1:nwflux) = w(iM^S,1:nwflux) - (w(iMp^S,1:nwflux) - w(iMm^S,1:nwflux)) * 0.25d0

    !> (eq.30) & (eq.31)
    deltap = 0
    deltam = 0
    phi1 = 0
    phi2 = 0
    phip = 1
    where(westp(iM^S,1:nwflux) .gt. w(iM^S,1:nwflux))
      deltap(iM^S,1:nwflux) = wmax(iM^S,1:nwflux) - w(iM^S,1:nwflux)
      deltam(iM^S,1:nwflux) = westp(iM^S,1:nwflux) - w(iM^S,1:nwflux)
      deltam(iM^S,1:nwflux) = sign(dabs(deltam(iM^S,1:nwflux)) + venk_omega, deltam(iM^S,1:nwflux))
      phi1(iM^S,1:nwflux) = (deltap(iM^S,1:nwflux)**2 + eps2) + 2.d0 * deltap(iM^S,1:nwflux) * deltam(iM^S,1:nwflux)
      phi2(iM^S,1:nwflux) = deltap(iM^S,1:nwflux)**2 + 2.d0 * deltam(iM^S,1:nwflux)**2 + deltap(iM^S,1:nwflux) * deltam(iM^S,1:nwflux) + eps2
      phip(iM^S,1:nwflux) = phi1(iM^S,1:nwflux) / phi2(iM^S, 1:nwflux)
    elsewhere(westp(iM^S,1:nwflux) .lt. w(iM^S,1:nwflux))
      deltap(iM^S,1:nwflux) = wmin(iM^S,1:nwflux) - w(iM^S,1:nwflux)
      deltam(iM^S,1:nwflux) = westp(iM^S,1:nwflux) - w(iM^S,1:nwflux)
      deltam(iM^S,1:nwflux) = sign(dabs(deltam(iM^S,1:nwflux)) + venk_omega, deltam(iM^S,1:nwflux))
      phi1(iM^S,1:nwflux) = (deltap(iM^S,1:nwflux)**2 + eps2) + 2.d0 * deltap(iM^S,1:nwflux) * deltam(iM^S,1:nwflux)
      phi2(iM^S,1:nwflux) = deltap(iM^S,1:nwflux)**2 + 2.d0 * deltam(iM^S,1:nwflux)**2 + deltap(iM^S,1:nwflux) * deltam(iM^S,1:nwflux) + eps2
      phip(iM^S,1:nwflux) = phi1(iM^S,1:nwflux) / phi2(iM^S, 1:nwflux)
    elsewhere
      phip(iM^S,1:nwflux) = one
   endwhere

    deltap = 0
    deltam = 0
    phi3 = 0
    phi4 = 0
    phim = 0
    where(westm(iM^S,1:nwflux) .lt. w(iM^S,1:nwflux))
      deltap(iM^S,1:nwflux) = - (wmax(iM^S,1:nwflux) - w(iM^S,1:nwflux))
      deltam(iM^S,1:nwflux) = westm(iM^S,1:nwflux) - w(iM^S,1:nwflux)
      deltam(iM^S,1:nwflux) = sign(dabs(deltam(iM^S,1:nwflux)) + venk_omega, deltam(iM^S,1:nwflux))
      phi3(iM^S,1:nwflux) = (deltap(iM^S,1:nwflux)**2 + eps2) + 2.d0 * deltap(iM^S,1:nwflux) * deltam(iM^S,1:nwflux)
      phi4(iM^S,1:nwflux) = deltap(iM^S,1:nwflux)**2 + 2.d0 * deltam(iM^S,1:nwflux)**2 + deltap(iM^S,1:nwflux) * deltam(iM^S,1:nwflux) + eps2
      phim(iM^S,1:nwflux) = phi3(iM^S,1:nwflux) / phi4(iM^S, 1:nwflux)
    elsewhere(westm(iM^S,1:nwflux) .gt. w(iM^S,1:nwflux))
      deltap(iM^S,1:nwflux) = - (wmin(iM^S,1:nwflux) - w(iM^S,1:nwflux))
      deltam(iM^S,1:nwflux) = westm(iM^S,1:nwflux) - w(iM^S,1:nwflux)
      deltam(iM^S,1:nwflux) = sign(dabs(deltam(iM^S,1:nwflux)) + venk_omega, deltam(iM^S,1:nwflux))
      phi3(iM^S,1:nwflux) = (deltap(iM^S,1:nwflux)**2 + eps2) + 2.d0 * deltap(iM^S,1:nwflux) * deltam(iM^S,1:nwflux)
      phi4(iM^S,1:nwflux) = deltap(iM^S,1:nwflux)**2 + 2.d0 * deltam(iM^S,1:nwflux)**2 + deltap(iM^S,1:nwflux) * deltam(iM^S,1:nwflux) + eps2
      phim(iM^S,1:nwflux) = phi3(iM^S,1:nwflux) / phi4(iM^S, 1:nwflux)
    elsewhere
      phim(iM^S,1:nwflux) = one
    endwhere
    !> (eq.3)
    phi(iM^S,1:nwflux) = min(phim(iM^S,1:nwflux),phip(iM^S,1:nwflux))
    !> (eq.5)
    wLC(iL^S,1:nwflux) = w(iL^S,1:nwflux) + 0.25d0 * (w(iLp^S,1:nwflux)-w(iLm^S,1:nwflux)) * phi(iL^S,1:nwflux)
    wRC(iL^S,1:nwflux) = w(iLp^S,1:nwflux) - 0.25d0 * (w(iLpp^S,1:nwflux)-w(iL^S,1:nwflux)) * phi(iLp^S,1:nwflux)

  end subroutine venklimiter

end module mod_venk
