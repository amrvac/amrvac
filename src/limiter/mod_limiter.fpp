!> Module with slope/flux limiters
module mod_limiter
  implicit none
  public

  !> radius of the asymptotic region [0.001, 10], larger means more accurate in smooth
  !> region but more overshooting at discontinuities
  double precision :: cada3_radius
  double precision :: schmid_rad1,schmid_rad2,schmid_rad3
  !$acc declare create(cada3_radius, schmid_rad1,schmid_rad2,schmid_rad3)
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

contains

  integer function limiter_type(namelim)
    use mod_comm_lib, only: mpistop
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

    case default
       limiter_type = -1
       write(*,*) 'Unknown limiter: ', namelim
       call mpistop("No such limiter")
    end select
  end function limiter_type


  !> Limit the centered dwC differences within ixC for iw in direction idim.
  !> The limiter is chosen according to typelim.
  !>
  !> Note that this subroutine is called from upwindLR (hence from methods
  !> like tvdlf, hancock, hll(c) etc) or directly from tvd.t,
  !> but also from the gradientS and divvectorS subroutines in geometry.t
  !> Accordingly, the typelim here corresponds to one of limiter
  !> or one of gradient_limiter.
  subroutine dwlimiter2(dwC,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idims,typelim,ldw,rdw,&
     a2max)

    use mod_global_parameters
    use mod_comm_lib, only: mpistop

    integer, value, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3, idims
    double precision, intent(in) :: dwC(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    integer, intent(in) :: typelim
    !> Result using left-limiter (same as right for symmetric)
    double precision, intent(out), optional :: ldw(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3)
    !> Result using right-limiter (same as left for symmetric)
    double precision, intent(out), optional :: rdw(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3)
    double precision, value, intent(in), optional :: a2max

    double precision :: tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
        tmp2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    integer :: ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, hxOmin1,&
       hxOmin2,hxOmin3,hxOmax1,hxOmax2,hxOmax3
    double precision, parameter :: qsmall=1.d-12, qsmall2=2.d-12
    double precision, parameter :: eps = sqrt(epsilon(1.0d0))

    ! mcbeta limiter parameter value
    double precision, parameter :: c_mcbeta=1.4d0
    ! cada limiter parameter values
    double precision, parameter :: cadalfa=0.5d0, cadbeta=2.0d0,&
        cadgamma=1.6d0
    ! full third order cada limiter
    double precision :: rdelinv
    double precision :: ldwA(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       ldwB(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       tmpeta(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision, parameter :: cadepsilon=1.d-14, invcadepsilon=1.d14,&
       cada3_radius=0.1d0
    integer :: ix1,ix2,ix3
    !-----------------------------------------------------------------------------

    ! Contract indices in idim for output.
    ixOmin1=ixCmin1+kr(idims,1);ixOmin2=ixCmin2+kr(idims,2)
    ixOmin3=ixCmin3+kr(idims,3); ixOmax1=ixCmax1;ixOmax2=ixCmax2
    ixOmax3=ixCmax3;
    hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
    hxOmin3=ixOmin3-kr(idims,3);hxOmax1=ixOmax1-kr(idims,1)
    hxOmax2=ixOmax2-kr(idims,2);hxOmax3=ixOmax3-kr(idims,3);

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
       tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=sign(one,&
          dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
       tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)=tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)* max(zero,min(abs(dwC(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,ixOmin3:ixOmax3)),tmp(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,ixOmin3:ixOmax3)*dwC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
          hxOmin3:hxOmax3)))
       if (present(ldw)) ldw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)
       if (present(rdw)) rdw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)
    case (limiter_woodward)
       ! Woodward and Collela limiter (eq.3.51h), a factor of 2 is pulled out
       tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=sign(one,&
          dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
       tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)=2*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)* max(zero,min(abs(dwC(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,ixOmin3:ixOmax3)),tmp(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,ixOmin3:ixOmax3)*dwC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
          hxOmin3:hxOmax3),tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)*quarter*(dwC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
          hxOmin3:hxOmax3)+dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3))))
       if (present(ldw)) ldw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)
       if (present(rdw)) rdw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)
    case (limiter_mcbeta)
       ! Woodward and Collela limiter, with factor beta
       tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=sign(one,&
          dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
       tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)=tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)* max(zero,min(c_mcbeta*abs(dwC(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,ixOmin3:ixOmax3)),c_mcbeta*tmp(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,ixOmin3:ixOmax3)*dwC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
          hxOmin3:hxOmax3),tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)*half*(dwC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
          hxOmin3:hxOmax3)+dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3))))
       if (present(ldw)) ldw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)
       if (present(rdw)) rdw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)
    case (limiter_superbee)
       ! Roes superbee limiter (eq.3.51i)
       tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=sign(one,&
          dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
       tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)=tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)* max(zero,min(2*abs(dwC(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,ixOmin3:ixOmax3)),tmp(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,ixOmin3:ixOmax3)*dwC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
          hxOmin3:hxOmax3)),min(abs(dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)),2*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)*dwC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
          hxOmin3:hxOmax3)))
       if (present(ldw)) ldw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)
       if (present(rdw)) rdw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)
    case (limiter_vanleer)
       ! van Leer limiter (eq 3.51f), but a missing delta2=1.D-12 is added
       tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)=2*max(dwC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
          hxOmin3:hxOmax3)*dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3),zero) /(dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)+dwC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
          hxOmin3:hxOmax3)+qsmall)
       if (present(ldw)) ldw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)
       if (present(rdw)) rdw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)
    case (limiter_albada)
       ! Albada limiter (eq.3.51g) with delta2=1D.-12
       tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)=(dwC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
          hxOmin3:hxOmax3)*(dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)**2+qsmall)+dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)*(dwC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
          hxOmin3:hxOmax3)**2+qsmall))/(dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)**2+dwC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
          hxOmin3:hxOmax3)**2+qsmall2)
       if (present(ldw)) ldw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)
       if (present(rdw)) rdw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)
    case (limiter_koren)
       tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=sign(one,&
          dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
       tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)=min(2*abs(dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)),2*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)*dwC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
          hxOmin3:hxOmax3))
       if (present(ldw)) then
          ldw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)=tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)* max(zero,min(tmp2(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3),(dwC(hxOmin1:hxOmax1,&
             hxOmin2:hxOmax2,hxOmin3:hxOmax3)*tmp(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3)+2*abs(dwC(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3)))*third))
       end if
       if (present(rdw)) then
          rdw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)=tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)* max(zero,min(tmp2(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3),(2*dwC(hxOmin1:hxOmax1,&
             hxOmin2:hxOmax2,hxOmin3:hxOmax3)*tmp(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3)+abs(dwC(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3)))*third))
       end if
    case (limiter_cada)
       ! This limiter has been rewritten in the usual form, and uses a division
       ! of the gradients.
       if (present(ldw)) then
          ! Cada Left variant
          ! Compute theta, but avoid division by zero
          tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)=dwC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
             hxOmin3:hxOmax3)/(dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) + sign(eps, dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)))
          tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)=(2+tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3))*third
          ldw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)= max(zero,&
             min(tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
              max(-cadalfa*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3), min(cadbeta*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3), tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3), cadgamma)))) * dwC(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3)
       end if

       if (present(rdw)) then
          ! Cada Right variant
          tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)=dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)/(dwC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
             hxOmin3:hxOmax3) + sign(eps, dwC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
             hxOmin3:hxOmax3)))
          tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)=(2+tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3))*third
          rdw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)= max(zero,&
             min(tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
              max(-cadalfa*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3), min(cadbeta*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3), tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3), cadgamma)))) * dwC(hxOmin1:hxOmax1,&
             hxOmin2:hxOmax2,hxOmin3:hxOmax3)
       end if
    case (limiter_cada3)
       rdelinv=one/(cada3_radius*dxlevel(idims))**2
       tmpeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)=(dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)**2+dwC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
          hxOmin3:hxOmax3)**2)*rdelinv
       if (present(ldw)) then
          tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)=dwC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
             hxOmin3:hxOmax3)/(dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) + sign(eps, dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)))
          ldwA(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)=(two+tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3))*third
          where(tmpeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)<=one-cadepsilon)
             ldw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3)=ldwA(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3)
          elsewhere(tmpeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)>=one+cadepsilon)
             ldwB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)= max(zero,&
                min(ldwA(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
                 max(-cadalfa*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3), min(cadbeta*tmp(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2,ixOmin3:ixOmax3), ldwA(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2,ixOmin3:ixOmax3), cadgamma))))
             ldw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3)=ldwB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3)
          elsewhere
             ldwB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)= max(zero,&
                min(ldwA(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
                 max(-cadalfa*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3), min(cadbeta*tmp(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2,ixOmin3:ixOmax3), ldwA(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2,ixOmin3:ixOmax3), cadgamma))))
             tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3)=(tmpeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3)-one)*invcadepsilon
             ldw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3)=half*( (one-tmp2(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2,ixOmin3:ixOmax3))*ldwA(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2,ixOmin3:ixOmax3) +(one+tmp2(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2,ixOmin3:ixOmax3))*ldwB(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2,ixOmin3:ixOmax3))
          endwhere
          ldw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)=ldw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)
       end if

       if (present(rdw)) then
          tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)=dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)/(dwC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
             hxOmin3:hxOmax3) + sign(eps, dwC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
             hxOmin3:hxOmax3)))
          ldwA(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)=(two+tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3))*third
          where(tmpeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)<=one-cadepsilon)
             rdw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3)=ldwA(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3)
          elsewhere(tmpeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)>=one+cadepsilon)
             ldwB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)= max(zero,&
                min(ldwA(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
                 max(-cadalfa*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3), min(cadbeta*tmp(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2,ixOmin3:ixOmax3), ldwA(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2,ixOmin3:ixOmax3), cadgamma))))
             rdw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3)=ldwB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3)
          elsewhere
             ldwB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)= max(zero,&
                min(ldwA(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
                 max(-cadalfa*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3), min(cadbeta*tmp(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2,ixOmin3:ixOmax3), ldwA(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2,ixOmin3:ixOmax3), cadgamma))))
             tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3)=(tmpeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3)-one)*invcadepsilon
             rdw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3)=half*( (one-tmp2(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2,ixOmin3:ixOmax3))*ldwA(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2,ixOmin3:ixOmax3) +(one+tmp2(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2,ixOmin3:ixOmax3))*ldwB(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2,ixOmin3:ixOmax3))
          endwhere
          rdw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)=rdw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * dwC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
             hxOmin3:hxOmax3)
       end if
    case(limiter_schmid)
      tmpeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=(sqrt(0.4d0*(dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)**2+dwC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
         hxOmin3:hxOmax3)**2)))/((a2max+cadepsilon)*dxlevel(idims)**2)
      if(present(ldw)) then
        tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)=dwC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
           hxOmin3:hxOmax3)/(dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)+sign(eps,dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)))
        ldwA(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)=(two+tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3))*third
        where(tmpeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)<=one-cadepsilon)
          ldw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)=ldwA(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)
        else where(tmpeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)>=one+cadepsilon)
          ldwB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=max(zero,&
             min(ldwA(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
             max(-tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
             min(cadbeta*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
             ldwA(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),1.5d0))))
          ldw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)=ldwB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)
        else where
          ldwB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=max(zero,&
             min(ldwA(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
             max(-tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
             min(cadbeta*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
             ldwA(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),1.5d0))))
          tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)=(tmpeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)-one)*invcadepsilon
          ldw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)=half*((one-tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3))*ldwA(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)+(one+tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3))*ldwB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3))
        end where
        ldw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)=ldw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)*dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
      end if
      if(present(rdw)) then
        tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)=dwC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)/(dwC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
           hxOmin3:hxOmax3)+sign(eps,dwC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
           hxOmin3:hxOmax3)))
        ldwA(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)=(two+tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3))*third
        where(tmpeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)<=one-cadepsilon)
          rdw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)=ldwA(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)
        else where(tmpeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)>=one+cadepsilon)
          ldwB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=max(zero,&
             min(ldwA(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
             max(-tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
             min(cadbeta*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
             ldwA(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),1.5d0))))
          rdw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)=ldwB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)
        else where
          ldwB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=max(zero,&
             min(ldwA(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
             max(-tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
             min(cadbeta*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
              ldwA(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3), 1.5d0))))
          tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)=(tmpeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)-one)*invcadepsilon
          rdw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)=half*((one-tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3))*ldwA(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)+(one+tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3))*ldwB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3))
        end where
        rdw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)=rdw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)*dwC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
           hxOmin3:hxOmax3)
       end if
    case default
       call mpistop("Error in dwLimiter: unknown limiter")
    end select

  end subroutine dwlimiter2


  
end module mod_limiter
