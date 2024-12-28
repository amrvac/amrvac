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
    double precision                :: wmax(ixI^S,1:nw),wmin(ixI^S,1:nw)
    double precision                :: westp(ixI^S,1:nw),westm(ixI^S,1:nw)
    double precision                :: phi1(ixI^S,1:nw),phi2(ixI^S,1:nw)
    double precision                :: phi3(ixI^S,1:nw),phi4(ixI^S,1:nw)
    double precision                :: phim(ixI^S,1:nw),phip(ixI^S,1:nw),phi(ixI^S,1:nw)
    double precision                :: deltap(ixI^S,1:nw),deltam(ixI^S,1:nw)
    double precision                :: eps2
    double precision, parameter     :: venk_omega = 1.d-12
    double precision, parameter     :: venk_k = 0.3d0
    integer                         :: iM^L, iMm^L, iMp^L
    integer                         :: iLm^L, iLp^L, iLpp^L

    iM^L=iL^L;
    iMmax^D=iLmax^D+kr(idims,^D);
    iLm^L=iL^L-kr(idims,^D);
    iMm^L=iM^L-kr(idims,^D);
    iLp^L=iL^L+kr(idims,^D);
    iMp^L=iM^L+kr(idims,^D);
    iLpp^L=iLp^L+kr(idims,^D);

    eps2 = (venk_k * dxdim) ** 3

    wmax(iM^S,1:nw_recon) = max(w(iMm^S,1:nw_recon), w(iM^S,1:nw_recon), w(iMp^S,1:nw_recon))
    wmin(iM^S,1:nw_recon) = min(w(iMm^S,1:nw_recon), w(iM^S,1:nw_recon), w(iMp^S,1:nw_recon))
    !> use central difference approximation as (eq.5) and take phi = 1
    westp(iM^S,1:nw_recon) = w(iM^S,1:nw_recon) + (w(iMp^S,1:nw_recon) - w(iMm^S,1:nw_recon)) * 0.25d0
    westm(iM^S,1:nw_recon) = w(iM^S,1:nw_recon) - (w(iMp^S,1:nw_recon) - w(iMm^S,1:nw_recon)) * 0.25d0

    !> (eq.30) & (eq.31)
    deltap = 0
    deltam = 0
    phi1 = 0
    phi2 = 0
    phip = 1
    where(westp(iM^S,1:nw_recon) .gt. w(iM^S,1:nw_recon))
      deltap(iM^S,1:nw_recon) = wmax(iM^S,1:nw_recon) - w(iM^S,1:nw_recon)
      deltam(iM^S,1:nw_recon) = westp(iM^S,1:nw_recon) - w(iM^S,1:nw_recon)
      deltam(iM^S,1:nw_recon) = sign(dabs(deltam(iM^S,1:nw_recon)) + venk_omega, deltam(iM^S,1:nw_recon))
      phi1(iM^S,1:nw_recon) = (deltap(iM^S,1:nw_recon)**2 + eps2) + 2.d0 * deltap(iM^S,1:nw_recon) * deltam(iM^S,1:nw_recon)
      phi2(iM^S,1:nw_recon) = deltap(iM^S,1:nw_recon)**2 + 2.d0 * deltam(iM^S,1:nw_recon)**2 + deltap(iM^S,1:nw_recon) * deltam(iM^S,1:nw_recon) + eps2
      phip(iM^S,1:nw_recon) = phi1(iM^S,1:nw_recon) / phi2(iM^S, 1:nw_recon)
    elsewhere(westp(iM^S,1:nw_recon) .lt. w(iM^S,1:nw_recon))
      deltap(iM^S,1:nw_recon) = wmin(iM^S,1:nw_recon) - w(iM^S,1:nw_recon)
      deltam(iM^S,1:nw_recon) = westp(iM^S,1:nw_recon) - w(iM^S,1:nw_recon)
      deltam(iM^S,1:nw_recon) = sign(dabs(deltam(iM^S,1:nw_recon)) + venk_omega, deltam(iM^S,1:nw_recon))
      phi1(iM^S,1:nw_recon) = (deltap(iM^S,1:nw_recon)**2 + eps2) + 2.d0 * deltap(iM^S,1:nw_recon) * deltam(iM^S,1:nw_recon)
      phi2(iM^S,1:nw_recon) = deltap(iM^S,1:nw_recon)**2 + 2.d0 * deltam(iM^S,1:nw_recon)**2 + deltap(iM^S,1:nw_recon) * deltam(iM^S,1:nw_recon) + eps2
      phip(iM^S,1:nw_recon) = phi1(iM^S,1:nw_recon) / phi2(iM^S, 1:nw_recon)
    elsewhere
      phip(iM^S,1:nw_recon) = one
   endwhere

    deltap = 0
    deltam = 0
    phi3 = 0
    phi4 = 0
    phim = 0
    where(westm(iM^S,1:nw_recon) .lt. w(iM^S,1:nw_recon))
      deltap(iM^S,1:nw_recon) = - (wmax(iM^S,1:nw_recon) - w(iM^S,1:nw_recon))
      deltam(iM^S,1:nw_recon) = westm(iM^S,1:nw_recon) - w(iM^S,1:nw_recon)
      deltam(iM^S,1:nw_recon) = sign(dabs(deltam(iM^S,1:nw_recon)) + venk_omega, deltam(iM^S,1:nw_recon))
      phi3(iM^S,1:nw_recon) = (deltap(iM^S,1:nw_recon)**2 + eps2) + 2.d0 * deltap(iM^S,1:nw_recon) * deltam(iM^S,1:nw_recon)
      phi4(iM^S,1:nw_recon) = deltap(iM^S,1:nw_recon)**2 + 2.d0 * deltam(iM^S,1:nw_recon)**2 + deltap(iM^S,1:nw_recon) * deltam(iM^S,1:nw_recon) + eps2
      phim(iM^S,1:nw_recon) = phi3(iM^S,1:nw_recon) / phi4(iM^S, 1:nw_recon)
    elsewhere(westm(iM^S,1:nw_recon) .gt. w(iM^S,1:nw_recon))
      deltap(iM^S,1:nw_recon) = - (wmin(iM^S,1:nw_recon) - w(iM^S,1:nw_recon))
      deltam(iM^S,1:nw_recon) = westm(iM^S,1:nw_recon) - w(iM^S,1:nw_recon)
      deltam(iM^S,1:nw_recon) = sign(dabs(deltam(iM^S,1:nw_recon)) + venk_omega, deltam(iM^S,1:nw_recon))
      phi3(iM^S,1:nw_recon) = (deltap(iM^S,1:nw_recon)**2 + eps2) + 2.d0 * deltap(iM^S,1:nw_recon) * deltam(iM^S,1:nw_recon)
      phi4(iM^S,1:nw_recon) = deltap(iM^S,1:nw_recon)**2 + 2.d0 * deltam(iM^S,1:nw_recon)**2 + deltap(iM^S,1:nw_recon) * deltam(iM^S,1:nw_recon) + eps2
      phim(iM^S,1:nw_recon) = phi3(iM^S,1:nw_recon) / phi4(iM^S, 1:nw_recon)
    elsewhere
      phim(iM^S,1:nw_recon) = one
    endwhere
    !> (eq.3)
    phi(iM^S,1:nw_recon) = min(phim(iM^S,1:nw_recon),phip(iM^S,1:nw_recon))
    !> (eq.5)
    wLC(iL^S,1:nw_recon) = w(iL^S,1:nw_recon) + 0.25d0 * (w(iLp^S,1:nw_recon)-w(iLm^S,1:nw_recon)) * phi(iL^S,1:nw_recon)
    wRC(iL^S,1:nw_recon) = w(iLp^S,1:nw_recon) - 0.25d0 * (w(iLpp^S,1:nw_recon)-w(iL^S,1:nw_recon)) * phi(iLp^S,1:nw_recon)

  end subroutine venklimiter

end module mod_venk
