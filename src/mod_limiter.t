!> Module with slope/flux limiters
module mod_limiter
  use mod_ppm
  use mod_mp5

  implicit none
  public

contains

  pure logical function limiter_symmetric(typelim)
    character(len=*), intent(in) :: typelim

    select case (typelim)
    case ("koren", "cada", "cada3")
       limiter_symmetric = .false.
    case default
       limiter_symmetric = .true.
    end select
  end function limiter_symmetric

  !> Limit the centered dwC differences within ixC for iw in direction idim.
  !> The limiter is chosen according to typelimiter.
  !>
  !> Note that this subroutine is called from upwindLR (hence from methods
  !> like tvdlf, hancock, hll(c) etc) or directly from tvd.t,
  !> but also from the gradientS and divvectorS subroutines in geometry.t
  !> Accordingly, the typelimiter here corresponds to one of limiter
  !> or one of gradient_limiter.
  subroutine dwlimiter2(dwC,ixI^L,ixC^L,idims,dxdim,qtypelimiter,ldw,rdw)

    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixC^L, idims
    double precision, intent(in) :: dxdim
    double precision, intent(in) :: dwC(ixI^S)
    character(len=std_len), intent(in) :: qtypelimiter
    !> Result using left-limiter (same as right for symmetric)
    double precision, intent(out), optional :: ldw(ixI^S)
    !> Result using right-limiter (same as right for symmetric)
    double precision, intent(out), optional :: rdw(ixI^S)

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

    select case (qtypelimiter)
    case ('minmod')
       ! Minmod limiter eq(3.51e) and (eq.3.38e) with omega=1
       tmp(ixO^S)=tmp(ixO^S)* &
            max(zero,min(dabs(dwC(ixO^S)),tmp(ixO^S)*dwC(hxO^S)))
       if (present(ldw)) ldw = tmp
       if (present(rdw)) rdw = tmp
    case ('woodward')
       ! Woodward and Collela limiter (eq.3.51h), a factor of 2 is pulled out
       tmp(ixO^S)=two*tmp(ixO^S)* &
            max(zero,min(dabs(dwC(ixO^S)),tmp(ixO^S)*dwC(hxO^S),&
            tmp(ixO^S)*quarter*(dwC(hxO^S)+dwC(ixO^S))))
       if (present(ldw)) ldw = tmp
       if (present(rdw)) rdw = tmp
    case ('mcbeta')
       ! Woodward and Collela limiter, with factor beta
       tmp(ixO^S)=tmp(ixO^S)* &
            max(zero,min(mcbeta*dabs(dwC(ixO^S)),mcbeta*tmp(ixO^S)*dwC(hxO^S),&
            tmp(ixO^S)*half*(dwC(hxO^S)+dwC(ixO^S))))
       if (present(ldw)) ldw = tmp
       if (present(rdw)) rdw = tmp
    case ('superbee')
       ! Roes superbee limiter (eq.3.51i)
       tmp(ixO^S)=tmp(ixO^S)* &
            max(zero,min(two*dabs(dwC(ixO^S)),tmp(ixO^S)*dwC(hxO^S)),&
            min(dabs(dwC(ixO^S)),two*tmp(ixO^S)*dwC(hxO^S)))
       if (present(ldw)) ldw = tmp
       if (present(rdw)) rdw = tmp
    case ('vanleer')
       ! van Leer limiter (eq 3.51f), but a missing delta2=1.D-12 is added
       tmp(ixO^S)=two*max(dwC(hxO^S)*dwC(ixO^S),zero) &
            /(dwC(ixO^S)+dwC(hxO^S)+qsmall)
       if (present(ldw)) ldw = tmp
       if (present(rdw)) rdw = tmp
    case ('albada')
       ! Albada limiter (eq.3.51g) with delta2=1D.-12
       tmp(ixO^S)=(dwC(hxO^S)*(dwC(ixO^S)**2+qsmall)&
            +dwC(ixO^S)*(dwC(hxO^S)**2+qsmall))&
            /(dwC(ixO^S)**2+dwC(hxO^S)**2+qsmall2)
       if (present(ldw)) ldw = tmp
       if (present(rdw)) rdw = tmp
    case ('koren')
       if (present(ldw)) then
          ldw(ixO^S)=tmp(ixO^S)* &
               max(zero,min(two*dabs(dwC(ixO^S)),two*tmp(ixO^S)*dwC(hxO^S),&
               (dwC(hxO^S)*tmp(ixO^S)+two*dabs(dwC(ixO^S)))*third))
       end if
       if (present(rdw)) then
          rdw(ixO^S)=tmp(ixO^S)* &
               max(zero,min(two*dabs(dwC(ixO^S)),two*tmp(ixO^S)*dwC(hxO^S),&
               (two*dwC(hxO^S)*tmp(ixO^S)+dabs(dwC(ixO^S)))*third))
       end if
    case ('cada')
       if (present(rdw)) then
          ! Cada Right variant
          rdw(ixO^S)=tmp(ixO^S)* &
               max(zero,min((two*dwC(hxO^S)*tmp(ixO^S)+dabs(dwC(ixO^S)))*third, &
               max(-cadalfa*dabs(dwC(ixO^S)),                     &
               min(cadbeta*dabs(dwC(ixO^S)),                  &
               (two*dwC(hxO^S)*tmp(ixO^S)+dabs(dwC(ixO^S)))*third, &
               cadgamma*tmp(ixO^S)*dwC(hxO^S)))))
       end if
       if (present(ldw)) then
          ! Cada Left variant
          ldw(ixO^S)=tmp(ixO^S)* &
               max(zero,min((two*dabs(dwC(ixO^S))+tmp(ixO^S)*dwC(hxO^S))*third, &
               max(-cadalfa*tmp(ixO^S)*dwC(hxO^S),                     &
               min(cadbeta*tmp(ixO^S)*dwC(hxO^S),                  &
               (two*dabs(dwC(ixO^S))+tmp(ixO^S)*dwC(hxO^S))*third, &
               cadgamma*dabs(dwC(ixO^S))))))
       end if
    case ('cada3')
       rdelinv=one/(cadradius*dxdim)**2
       tmpeta(ixO^S)=(dwC(ixO^S)**2+dwC(hxO^S)**2)*rdelinv

       if (present(ldw)) then
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
       end if

       if (present(rdw)) then
          ldwA(ixO^S)=(two*dwC(hxO^S)+dwC(ixO^S))*third
          ldwB(ixO^S)=tmp(ixO^S)* &
               max(zero,min((two*dwC(hxO^S)*tmp(ixO^S)+dabs(dwC(ixO^S)))*third, &
               max(-cadalfa*dabs(dwC(ixO^S)),                     &
               min(cadbeta*dabs(dwC(ixO^S)),                  &
               (two*dwC(hxO^S)*tmp(ixO^S)+dabs(dwC(ixO^S)))*third, &
               cadgamma*tmp(ixO^S)*dwC(hxO^S)))))
          where(tmpeta(ixO^S)<=one-cadepsilon)
             rdw(ixO^S)=ldwA(ixO^S)
          elsewhere(tmpeta(ixO^S)>=one+cadepsilon)
             rdw(ixO^S)=ldwB(ixO^S)
          elsewhere
             tmp(ixO^S)=(tmpeta(ixO^S)-one)*invcadepsilon
             rdw(ixO^S)=half*( (one-tmp(ixO^S))*ldwA(ixO^S) &
                  +(one+tmp(ixO^S))*ldwB(ixO^S))
          endwhere
       end if

    case default
       write(*,*)'Unknown limiter:',qtypelimiter
       call mpistop("Error in dwLimiter: No such TVD limiter")
    end select

  end subroutine dwlimiter2

end module mod_limiter
