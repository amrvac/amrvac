module mod_limiter
  use mod_ppm
  use mod_mp5

  implicit none
  public

contains

  !> Limit the centered dwC differences within ixC for iw in direction idim.
  !> The limiter is chosen according to typelimiter.
  !>
  !> Note that this subroutine is called from upwindLR (hence from methods
  !> like tvdlf, hancock, hll(c) etc) or directly from tvd.t,
  !> but also from the gradientS and divvectorS subroutines in geometry.t
  !> Accordingly, the typelimiter here corresponds to one of typelimiter1
  !> or one of typegradlimiter1.
  !>
  !> note: there is no iw dependence here...
  subroutine dwlimiter2(dwC,ixI^L,ixC^L,iw,idims,ldw,dxdim)

    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixC^L, iw, idims
    double precision, intent(in) :: dxdim
    double precision, intent(in) :: dwC(ixI^S)
    double precision, intent(out) :: ldw(ixI^S)

    double precision :: tmp(ixI^S)
    integer :: ixO^L, hxO^L
    double precision, parameter :: qsmall=1.d-12, qsmall2=2.d-12

    ! cada limiter parameter values
    double precision, parameter :: cadalfa=0.5d0, cadbeta=2.0d0, cadgamma=1.6d0
    ! full third order cada limiter
    double precision :: rdelinv
    double precision :: ldwA(ixI^S),ldwB(ixI^S),tmpeta(ixI^S)
    double precision, parameter :: cadepsilon=1.d-14, invcadepsilon=1.d14,cadradius=0.1d0
    !-----------------------------------------------------------------------------

    ! Contract indices in idim for output.
    ixOmin^D=ixCmin^D+kr(idims,^D); ixOmax^D=ixCmax^D;
    hxO^L=ixO^L-kr(idims,^D);

    ! Store the sign of dwC in tmp
    tmp(ixO^S)=sign(one,dwC(ixO^S))
    rdelinv=one/(cadradius*dxdim)**2

    select case (typelimiter)
    case ('minmod')
       ! Minmod limiter eq(3.51e) and (eq.3.38e) with omega=1
       ldw(ixO^S)=tmp(ixO^S)* &
            max(zero,min(dabs(dwC(ixO^S)),tmp(ixO^S)*dwC(hxO^S)))
    case ('woodward')
       ! Woodward and Collela limiter (eq.3.51h), a factor of 2 is pulled out
       ldw(ixO^S)=two*tmp(ixO^S)* &
            max(zero,min(dabs(dwC(ixO^S)),tmp(ixO^S)*dwC(hxO^S),&
            tmp(ixO^S)*quarter*(dwC(hxO^S)+dwC(ixO^S))))
    case ('mcbeta')
       ! Woodward and Collela limiter, with factor beta
       ldw(ixO^S)=tmp(ixO^S)* &
            max(zero,min(mcbeta*dabs(dwC(ixO^S)),mcbeta*tmp(ixO^S)*dwC(hxO^S),&
            tmp(ixO^S)*half*(dwC(hxO^S)+dwC(ixO^S))))
    case ('superbee')
       ! Roes superbee limiter (eq.3.51i)
       ldw(ixO^S)=tmp(ixO^S)* &
            max(zero,min(two*dabs(dwC(ixO^S)),tmp(ixO^S)*dwC(hxO^S)),&
            min(dabs(dwC(ixO^S)),two*tmp(ixO^S)*dwC(hxO^S)))
    case ('vanleer')
       ! van Leer limiter (eq 3.51f), but a missing delta2=1.D-12 is added
       ldw(ixO^S)=two*max(dwC(hxO^S)*dwC(ixO^S),zero) &
            /(dwC(ixO^S)+dwC(hxO^S)+qsmall)
    case ('albada')
       ! Albada limiter (eq.3.51g) with delta2=1D.-12
       ldw(ixO^S)=(dwC(hxO^S)*(dwC(ixO^S)**2+qsmall)&
            +dwC(ixO^S)*(dwC(hxO^S)**2+qsmall))&
            /(dwC(ixO^S)**2+dwC(hxO^S)**2+qsmall2)
    case ('korenR')
       ! Barry Koren Right variant
       ldw(ixO^S)=tmp(ixO^S)* &
            max(zero,min(two*dabs(dwC(ixO^S)),two*tmp(ixO^S)*dwC(hxO^S),&
            (two*dwC(hxO^S)*tmp(ixO^S)+dabs(dwC(ixO^S)))*third))
    case ('korenL')
       ! Barry Koren Left variant
       ldw(ixO^S)=tmp(ixO^S)* &
            max(zero,min(two*dabs(dwC(ixO^S)),two*tmp(ixO^S)*dwC(hxO^S),&
            (dwC(hxO^S)*tmp(ixO^S)+two*dabs(dwC(ixO^S)))*third))
    case ('cadaR')
       ! Cada Right variant
       ldw(ixO^S)=tmp(ixO^S)* &
            max(zero,min((two*dwC(hxO^S)*tmp(ixO^S)+dabs(dwC(ixO^S)))*third, &
            max(-cadalfa*dabs(dwC(ixO^S)),                     &
            min(cadbeta*dabs(dwC(ixO^S)),                  &
            (two*dwC(hxO^S)*tmp(ixO^S)+dabs(dwC(ixO^S)))*third, &
            cadgamma*tmp(ixO^S)*dwC(hxO^S)))))
    case ('cadaL')
       ! Cada Left variant
       ldw(ixO^S)=tmp(ixO^S)* &
            max(zero,min((two*dabs(dwC(ixO^S))+tmp(ixO^S)*dwC(hxO^S))*third, &
            max(-cadalfa*tmp(ixO^S)*dwC(hxO^S),                     &
            min(cadbeta*tmp(ixO^S)*dwC(hxO^S),                  &
            (two*dabs(dwC(ixO^S))+tmp(ixO^S)*dwC(hxO^S))*third, &
            cadgamma*dabs(dwC(ixO^S))))))
    case ('cada3R')
       tmpeta(ixO^S)=(dwC(ixO^S)**2+dwC(hxO^S)**2)*rdelinv
       ldwA(ixO^S)=(two*dwC(hxO^S)+dwC(ixO^S))*third
       ldwB(ixO^S)=tmp(ixO^S)* &
            max(zero,min((two*dwC(hxO^S)*tmp(ixO^S)+dabs(dwC(ixO^S)))*third, &
            max(-cadalfa*dabs(dwC(ixO^S)),                     &
            min(cadbeta*dabs(dwC(ixO^S)),                  &
            (two*dwC(hxO^S)*tmp(ixO^S)+dabs(dwC(ixO^S)))*third, &
            cadgamma*tmp(ixO^S)*dwC(hxO^S)))))
       where(tmpeta(ixO^S)<=one-cadepsilon)
          ldw(ixO^S)=ldwA(ixO^S)
       elsewhere(tmpeta(ixO^S)>=one+cadepsilon)
          ldw(ixO^S)=ldwB(ixO^S)
       elsewhere
          tmp(ixO^S)=(tmpeta(ixO^S)-one)*invcadepsilon
          ldw(ixO^S)=half*( (one-tmp(ixO^S))*ldwA(ixO^S) &
               +(one+tmp(ixO^S))*ldwB(ixO^S))
       endwhere
    case ('cada3L')
       tmpeta(ixO^S)=(dwC(ixO^S)**2+dwC(hxO^S)**2)*rdelinv
       ldwA(ixO^S)=(two*dwC(ixO^S)+dwC(hxO^S))*third
       ldwB(ixO^S)=tmp(ixO^S)* &
            max(zero,min((two*dabs(dwC(ixO^S))+tmp(ixO^S)*dwC(hxO^S))*third, &
            max(-cadalfa*tmp(ixO^S)*dwC(hxO^S),                     &
            min(cadbeta*tmp(ixO^S)*dwC(hxO^S),                  &
            (two*dabs(dwC(ixO^S))+tmp(ixO^S)*dwC(hxO^S))*third, &
            cadgamma*dabs(dwC(ixO^S))))))
       where(tmpeta(ixO^S)<=one-cadepsilon)
          ldw(ixO^S)=ldwA(ixO^S)
       elsewhere(tmpeta(ixO^S)>=one+cadepsilon)
          ldw(ixO^S)=ldwB(ixO^S)
       elsewhere
          tmp(ixO^S)=(tmpeta(ixO^S)-one)*invcadepsilon
          ldw(ixO^S)=half*( (one-tmp(ixO^S))*ldwA(ixO^S) &
               +(one+tmp(ixO^S))*ldwB(ixO^S))
       endwhere
    case default
       write(*,*)'Unknown limiter:',typelimiter
       call mpistop("Error in dwLimiter: No such TVD limiter")
    end select

  end subroutine dwlimiter2

end module mod_limiter
