!> Module with slope/flux limiters
module mod_limiter
  use mod_ppm
  use mod_mp5
  use mod_weno
  use mod_venk

  implicit none
  public

  !> radius of the asymptotic region [0.001, 10], larger means more accurate in smooth
  !> region but more overshooting at discontinuities
  double precision :: cada3_radius
  double precision :: schmid_rad^D
  integer, parameter :: limiter_minmod = 1
  integer, parameter :: limiter_woodward = 2
  integer, parameter :: limiter_mcbeta = 3
  integer, parameter :: limiter_superbee = 4
  integer, parameter :: limiter_vanleer = 5
  integer, parameter :: limiter_albada = 6
  integer, parameter :: limiter_koren = 7
  integer, parameter :: limiter_cada = 8
  integer, parameter :: limiter_cada3 = 9
  integer, parameter :: limiter_schmid = 10
  integer, parameter :: limiter_venk = 11
  ! Special cases
  integer, parameter :: limiter_ppm = 12
  integer, parameter :: limiter_mp5 = 13
  integer, parameter :: limiter_weno3  = 14
  integer, parameter :: limiter_wenoyc3  = 15
  integer, parameter :: limiter_weno5 = 16
  integer, parameter :: limiter_weno5nm = 17
  integer, parameter :: limiter_wenoz5  = 18
  integer, parameter :: limiter_wenoz5nm = 19
  integer, parameter :: limiter_wenozp5  = 20
  integer, parameter :: limiter_wenozp5nm = 21
  integer, parameter :: limiter_weno5cu6 = 22
  integer, parameter :: limiter_teno5ad = 23
  integer, parameter :: limiter_weno7 = 24
  integer, parameter :: limiter_mpweno7 = 25
  integer, parameter :: limiter_exeno7 = 26

contains

  integer function limiter_type(namelim)
    character(len=*), intent(in) :: namelim

    select case (namelim)
    case ('minmod')
       limiter_type = limiter_minmod
    case ('woodward')
       limiter_type = limiter_woodward
    case ('mcbeta')
       limiter_type = limiter_mcbeta
    case ('superbee')
       limiter_type = limiter_superbee
    case ('vanleer')
       limiter_type = limiter_vanleer
    case ('albada')
       limiter_type = limiter_albada
    case ('koren')
       limiter_type = limiter_koren
    case ('cada')
       limiter_type = limiter_cada
    case ('cada3')
       limiter_type = limiter_cada3
    case ('schmid1')
       limiter_type = limiter_schmid
    case ('schmid2')
       limiter_type = limiter_schmid
    case('venk')
       limiter_type = limiter_venk
    case ('ppm')
       limiter_type = limiter_ppm
    case ('mp5')
       limiter_type = limiter_mp5
    case ('weno3')
       limiter_type = limiter_weno3
    case ('wenoyc3')
       limiter_type = limiter_wenoyc3
    case ('weno5')
       limiter_type = limiter_weno5
    case ('weno5nm')
       limiter_type = limiter_weno5nm
    case ('wenoz5')
       limiter_type = limiter_wenoz5
    case ('wenoz5nm')
       limiter_type = limiter_wenoz5nm
    case ('wenozp5')
       limiter_type = limiter_wenozp5
    case ('wenozp5nm')
       limiter_type = limiter_wenozp5nm
    case ('weno5cu6')
       limiter_type = limiter_weno5cu6
    case ('teno5ad')
       limiter_type = limiter_teno5ad
    case ('weno7')
       limiter_type = limiter_weno7
    case ('mpweno7')
       limiter_type = limiter_mpweno7
    case ('exeno7')
       limiter_type = limiter_exeno7

    case default
       limiter_type = -1
       write(*,*) 'Unknown limiter: ', namelim
       call mpistop("No such limiter")
    end select
  end function limiter_type

  pure logical function limiter_symmetric(typelim)
    integer, intent(in) :: typelim

    select case (typelim)
    case (limiter_koren, limiter_cada, limiter_cada3)
       limiter_symmetric = .false.
    case default
       limiter_symmetric = .true.
    end select
  end function limiter_symmetric

  !> Limit the centered dwC differences within ixC for iw in direction idim.
  !> The limiter is chosen according to typelim.
  !>
  !> Note that this subroutine is called from upwindLR (hence from methods
  !> like tvdlf, hancock, hll(c) etc) or directly from tvd.t,
  !> but also from the gradientS and divvectorS subroutines in geometry.t
  !> Accordingly, the typelim here corresponds to one of limiter
  !> or one of gradient_limiter.
  subroutine dwlimiter2(dwC,ixI^L,ixC^L,idims,typelim,ldw,rdw,a2max)

    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixC^L, idims
    double precision, intent(in) :: dwC(ixI^S)
    integer, intent(in) :: typelim
    !> Result using left-limiter (same as right for symmetric)
    double precision, intent(out), optional :: ldw(ixI^S)
    !> Result using right-limiter (same as left for symmetric)
    double precision, intent(out), optional :: rdw(ixI^S)
    double precision, intent(in), optional :: a2max

    double precision :: tmp(ixI^S), tmp2(ixI^S)
    integer :: ixO^L, hxO^L
    double precision, parameter :: qsmall=1.d-12, qsmall2=2.d-12
    double precision, parameter :: eps = sqrt(epsilon(1.0d0))

    ! mcbeta limiter parameter value
    double precision, parameter :: c_mcbeta=1.4d0
    ! cada limiter parameter values
    double precision, parameter :: cadalfa=0.5d0, cadbeta=2.0d0, cadgamma=1.6d0
    ! full third order cada limiter
    double precision :: rdelinv
    double precision :: ldwA(ixI^S),ldwB(ixI^S),tmpeta(ixI^S)
    double precision, parameter :: cadepsilon=1.d-14, invcadepsilon=1.d14,cada3_radius=0.1d0
    integer :: ix^D
    !-----------------------------------------------------------------------------

    ! Contract indices in idim for output.
    ixOmin^D=ixCmin^D+kr(idims,^D); ixOmax^D=ixCmax^D;
    hxO^L=ixO^L-kr(idims,^D);

    ! About the notation: the conventional argument theta (the ratio of slopes)
    ! would be given by dwC(ixO^S)/dwC(hxO^S). However, in the end one
    ! multiplies phi(theta) by dwC(hxO^S), which is incorporated in the
    ! equations below. The minmod limiter can for example be written as:
    ! A:
    ! max(0.0d0, min(1.0d0, dwC(ixO^S)/dwC(hxO^S))) * dwC(hxO^S)
    ! B:
    ! tmp(ixO^S)*max(0.0d0,min(abs(dwC(ixO^S)),tmp(ixO^S)*dwC(hxO^S)))
    ! where tmp(ixO^S)=sign(1.0d0,dwC(ixO^S))

    select case (typelim)
    case (limiter_minmod)
       ! Minmod limiter eq(3.51e) and (eq.3.38e) with omega=1
       tmp(ixO^S)=sign(one,dwC(ixO^S))
       tmp(ixO^S)=tmp(ixO^S)* &
            max(zero,min(abs(dwC(ixO^S)),tmp(ixO^S)*dwC(hxO^S)))
       if (present(ldw)) ldw(ixO^S) = tmp(ixO^S)
       if (present(rdw)) rdw(ixO^S) = tmp(ixO^S)
    case (limiter_woodward)
       ! Woodward and Collela limiter (eq.3.51h), a factor of 2 is pulled out
       tmp(ixO^S)=sign(one,dwC(ixO^S))
       tmp(ixO^S)=2*tmp(ixO^S)* &
            max(zero,min(abs(dwC(ixO^S)),tmp(ixO^S)*dwC(hxO^S),&
            tmp(ixO^S)*quarter*(dwC(hxO^S)+dwC(ixO^S))))
       if (present(ldw)) ldw(ixO^S) = tmp(ixO^S)
       if (present(rdw)) rdw(ixO^S) = tmp(ixO^S)
    case (limiter_mcbeta)
       ! Woodward and Collela limiter, with factor beta
       tmp(ixO^S)=sign(one,dwC(ixO^S))
       tmp(ixO^S)=tmp(ixO^S)* &
            max(zero,min(c_mcbeta*abs(dwC(ixO^S)),c_mcbeta*tmp(ixO^S)*dwC(hxO^S),&
            tmp(ixO^S)*half*(dwC(hxO^S)+dwC(ixO^S))))
       if (present(ldw)) ldw(ixO^S) = tmp(ixO^S)
       if (present(rdw)) rdw(ixO^S) = tmp(ixO^S)
    case (limiter_superbee)
       ! Roes superbee limiter (eq.3.51i)
       tmp(ixO^S)=sign(one,dwC(ixO^S))
       tmp(ixO^S)=tmp(ixO^S)* &
            max(zero,min(2*abs(dwC(ixO^S)),tmp(ixO^S)*dwC(hxO^S)),&
            min(abs(dwC(ixO^S)),2*tmp(ixO^S)*dwC(hxO^S)))
       if (present(ldw)) ldw(ixO^S) = tmp(ixO^S)
       if (present(rdw)) rdw(ixO^S) = tmp(ixO^S)
    case (limiter_vanleer)
       ! van Leer limiter (eq 3.51f), but a missing delta2=1.D-12 is added
       tmp(ixO^S)=2*max(dwC(hxO^S)*dwC(ixO^S),zero) &
            /(dwC(ixO^S)+dwC(hxO^S)+qsmall)
       if (present(ldw)) ldw(ixO^S) = tmp(ixO^S)
       if (present(rdw)) rdw(ixO^S) = tmp(ixO^S)
    case (limiter_albada)
       ! Albada limiter (eq.3.51g) with delta2=1D.-12
       tmp(ixO^S)=(dwC(hxO^S)*(dwC(ixO^S)**2+qsmall)&
            +dwC(ixO^S)*(dwC(hxO^S)**2+qsmall))&
            /(dwC(ixO^S)**2+dwC(hxO^S)**2+qsmall2)
       if (present(ldw)) ldw(ixO^S) = tmp(ixO^S)
       if (present(rdw)) rdw(ixO^S) = tmp(ixO^S)
    case (limiter_koren)
       tmp(ixO^S)=sign(one,dwC(ixO^S))
       tmp2(ixO^S)=min(2*abs(dwC(ixO^S)),2*tmp(ixO^S)*dwC(hxO^S))
       if (present(ldw)) then
          ldw(ixO^S)=tmp(ixO^S)* &
               max(zero,min(tmp2(ixO^S),&
               (dwC(hxO^S)*tmp(ixO^S)+2*abs(dwC(ixO^S)))*third))
       end if
       if (present(rdw)) then
          rdw(ixO^S)=tmp(ixO^S)* &
                max(zero,min(tmp2(ixO^S),&
               (2*dwC(hxO^S)*tmp(ixO^S)+abs(dwC(ixO^S)))*third))
       end if
    case (limiter_cada)
       ! This limiter has been rewritten in the usual form, and uses a division
       ! of the gradients.
       if (present(ldw)) then
          ! Cada Left variant
          ! Compute theta, but avoid division by zero
          tmp(ixO^S)=dwC(hxO^S)/(dwC(ixO^S) + sign(eps, dwC(ixO^S)))
          tmp2(ixO^S)=(2+tmp(ixO^S))*third
          ldw(ixO^S)= max(zero,min(tmp2(ixO^S), &
               max(-cadalfa*tmp(ixO^S), &
               min(cadbeta*tmp(ixO^S), tmp2(ixO^S), &
               cadgamma)))) * dwC(ixO^S)
       end if

       if (present(rdw)) then
          ! Cada Right variant
          tmp(ixO^S)=dwC(ixO^S)/(dwC(hxO^S) + sign(eps, dwC(hxO^S)))
          tmp2(ixO^S)=(2+tmp(ixO^S))*third
          rdw(ixO^S)= max(zero,min(tmp2(ixO^S), &
               max(-cadalfa*tmp(ixO^S), &
               min(cadbeta*tmp(ixO^S), tmp2(ixO^S), &
               cadgamma)))) * dwC(hxO^S)
       end if
    case (limiter_cada3)
       rdelinv=one/(cada3_radius*dxlevel(idims))**2
       tmpeta(ixO^S)=(dwC(ixO^S)**2+dwC(hxO^S)**2)*rdelinv
       if (present(ldw)) then
          tmp(ixO^S)=dwC(hxO^S)/(dwC(ixO^S) + sign(eps, dwC(ixO^S)))
          ldwA(ixO^S)=(two+tmp(ixO^S))*third
          where(tmpeta(ixO^S)<=one-cadepsilon)
             ldw(ixO^S)=ldwA(ixO^S)
          elsewhere(tmpeta(ixO^S)>=one+cadepsilon)
             ldwB(ixO^S)= max(zero,min(ldwA(ixO^S), max(-cadalfa*tmp(ixO^S), &
               min(cadbeta*tmp(ixO^S), ldwA(ixO^S), cadgamma))))
             ldw(ixO^S)=ldwB(ixO^S)
          elsewhere
             ldwB(ixO^S)= max(zero,min(ldwA(ixO^S), max(-cadalfa*tmp(ixO^S), &
               min(cadbeta*tmp(ixO^S), ldwA(ixO^S), cadgamma))))
             tmp2(ixO^S)=(tmpeta(ixO^S)-one)*invcadepsilon
             ldw(ixO^S)=half*( (one-tmp2(ixO^S))*ldwA(ixO^S) &
                  +(one+tmp2(ixO^S))*ldwB(ixO^S))
          endwhere
          ldw(ixO^S)=ldw(ixO^S) * dwC(ixO^S)
       end if

       if (present(rdw)) then
          tmp(ixO^S)=dwC(ixO^S)/(dwC(hxO^S) + sign(eps, dwC(hxO^S)))
          ldwA(ixO^S)=(two+tmp(ixO^S))*third
          where(tmpeta(ixO^S)<=one-cadepsilon)
             rdw(ixO^S)=ldwA(ixO^S)
          elsewhere(tmpeta(ixO^S)>=one+cadepsilon)
             ldwB(ixO^S)= max(zero,min(ldwA(ixO^S), max(-cadalfa*tmp(ixO^S), &
               min(cadbeta*tmp(ixO^S), ldwA(ixO^S), cadgamma))))
             rdw(ixO^S)=ldwB(ixO^S)
          elsewhere
             ldwB(ixO^S)= max(zero,min(ldwA(ixO^S), max(-cadalfa*tmp(ixO^S), &
               min(cadbeta*tmp(ixO^S), ldwA(ixO^S), cadgamma))))
             tmp2(ixO^S)=(tmpeta(ixO^S)-one)*invcadepsilon
             rdw(ixO^S)=half*( (one-tmp2(ixO^S))*ldwA(ixO^S) &
                  +(one+tmp2(ixO^S))*ldwB(ixO^S))
          endwhere
          rdw(ixO^S)=rdw(ixO^S) * dwC(hxO^S)
       end if
    case(limiter_schmid)
      tmpeta(ixO^S)=(sqrt(0.4d0*(dwC(ixO^S)**2+dwC(hxO^S)**2)))&
        /((a2max+cadepsilon)*dxlevel(idims)**2)
      if(present(ldw)) then
        tmp(ixO^S)=dwC(hxO^S)/(dwC(ixO^S)+sign(eps,dwC(ixO^S)))
        ldwA(ixO^S)=(two+tmp(ixO^S))*third
        where(tmpeta(ixO^S)<=one-cadepsilon)
          ldw(ixO^S)=ldwA(ixO^S)
        else where(tmpeta(ixO^S)>=one+cadepsilon)
          ldwB(ixO^S)=max(zero,min(ldwA(ixO^S),max(-tmp(ixO^S),&
             min(cadbeta*tmp(ixO^S),ldwA(ixO^S),1.5d0))))
          ldw(ixO^S)=ldwB(ixO^S)
        else where
          ldwB(ixO^S)=max(zero,min(ldwA(ixO^S),max(-tmp(ixO^S),&
             min(cadbeta*tmp(ixO^S),ldwA(ixO^S),1.5d0))))
          tmp2(ixO^S)=(tmpeta(ixO^S)-one)*invcadepsilon
          ldw(ixO^S)=half*((one-tmp2(ixO^S))*ldwA(ixO^S)&
            +(one+tmp2(ixO^S))*ldwB(ixO^S))
        end where
        ldw(ixO^S)=ldw(ixO^S)*dwC(ixO^S)
      end if
      if(present(rdw)) then
        tmp(ixO^S)=dwC(ixO^S)/(dwC(hxO^S)+sign(eps,dwC(hxO^S)))
        ldwA(ixO^S)=(two+tmp(ixO^S))*third
        where(tmpeta(ixO^S)<=one-cadepsilon)
          rdw(ixO^S)=ldwA(ixO^S)
        else where(tmpeta(ixO^S)>=one+cadepsilon)
          ldwB(ixO^S)=max(zero,min(ldwA(ixO^S),max(-tmp(ixO^S),&
             min(cadbeta*tmp(ixO^S),ldwA(ixO^S),1.5d0))))
          rdw(ixO^S)=ldwB(ixO^S)
        else where
          ldwB(ixO^S)=max(zero,min(ldwA(ixO^S),max(-tmp(ixO^S),&
            min(cadbeta*tmp(ixO^S), ldwA(ixO^S), 1.5d0))))
          tmp2(ixO^S)=(tmpeta(ixO^S)-one)*invcadepsilon
          rdw(ixO^S)=half*((one-tmp2(ixO^S))*ldwA(ixO^S)&
            +(one+tmp2(ixO^S))*ldwB(ixO^S))
        end where
        rdw(ixO^S)=rdw(ixO^S)*dwC(hxO^S)
       end if
    case default
       call mpistop("Error in dwLimiter: unknown limiter")
    end select

  end subroutine dwlimiter2

end module mod_limiter
