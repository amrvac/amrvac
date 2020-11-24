!> Module with Roe-type Riemann solver for hydrodynamics
module mod_hd_roe
  use mod_hd_phys
  use mod_physics_roe

  implicit none
  private

  integer :: soundRW_ = -1
  integer :: soundLW_ = -1
  integer :: entropW_ = -1
  integer :: shearW0_ = -1

  public :: hd_roe_init

contains

  subroutine hd_roe_init()
    use mod_global_parameters, only: entropycoef, nw

    integer :: iw

    if (hd_energy) then
       ! Characteristic waves
       soundRW_ = 1
       soundLW_ = 2
       entropW_ = 3
       shearW0_ = 3
       nworkroe = 3

       phys_average => hd_average
       phys_get_eigenjump => hd_get_eigenjump
       phys_rtimes => hd_rtimes
    else
       ! Characteristic waves
       soundRW_ = 1
       soundLW_ = 2
       shearW0_ = 2
       nworkroe = 1

       phys_average => hd_average_iso
       phys_get_eigenjump => hd_get_eigenjump_iso
       phys_rtimes => hd_rtimes_iso
    end if

    allocate(entropycoef(nw))

    do iw = 1, nw
       if (iw == soundRW_ .or. iw == soundLW_) then
          ! TODO: Jannis: what's this?
          entropycoef(iw) = 0.2d0
       else
          entropycoef(iw) = -1.0d0
       end if
    end do

  end subroutine hd_roe_init

  !> Calculate the Roe average of w, assignment of variables:
  !> rho -> rho, m -> v, e -> h
  subroutine hd_average(wL,wR,x,ix^L,idim,wroe,workroe)
    use mod_global_parameters
    integer, intent(in)             :: ix^L, idim
    double precision, intent(in)    :: wL(ixG^T, nw), wR(ixG^T, nw)
    double precision, intent(inout) :: wroe(ixG^T, nw)
    double precision, intent(inout) :: workroe(ixG^T, nworkroe)
    double precision, intent(in)    :: x(ixG^T, 1:^ND)
    integer                         :: idir

    ! call average2(wL,wR,x,ix^L,idim,wroe,workroe(ixG^T,1),workroe(ixG^T,2))
    workroe(ix^S, 1) = sqrt(wL(ix^S,rho_))
    workroe(ix^S, 2) = sqrt(wR(ix^S,rho_))

    ! The averaged density is sqrt(rhoL*rhoR)
    wroe(ix^S,rho_)  = workroe(ix^S, 1)*workroe(ix^S, 2)

    ! Now the ratio sqrt(rhoL/rhoR) is put into workroe(ix^S, 1)
    workroe(ix^S, 1) = workroe(ix^S, 1)/workroe(ix^S, 2)

    ! Roe-average velocities
    do idir = 1, ndir
       wroe(ix^S,mom(idir)) = (wL(ix^S,mom(idir))/wL(ix^S,rho_) * workroe(ix^S, 1)+&
            wR(ix^S,mom(idir))/wR(ix^S,rho_))/(one+workroe(ix^S, 1))
    end do

    ! Calculate enthalpyL, then enthalpyR, then Roe-average. Use tmp2 for pressure.
    call hd_get_pthermal(wL,x,ixG^LL,ix^L, workroe(ixG^T, 2))

    wroe(ix^S,e_)    = (workroe(ix^S, 2)+wL(ix^S,e_))/wL(ix^S,rho_)

    call hd_get_pthermal(wR,x,ixG^LL,ix^L, workroe(ixG^T, 2))

    workroe(ix^S, 2) = (workroe(ix^S, 2)+wR(ix^S,e_))/wR(ix^S,rho_)
    wroe(ix^S,e_)    = (wroe(ix^S,e_)*workroe(ix^S, 1) + workroe(ix^S, 2))/(one+workroe(ix^S, 1))
  end subroutine hd_average

  subroutine average2(wL,wR,x,ix^L,idim,wroe,tmp,tmp2)

    ! Calculate the Roe average of w, assignment of variables:
    ! rho -> rho, m -> v, e -> h

    use mod_global_parameters

    integer                                             :: ix^L,idim,idir
    double precision, dimension(ixG^T,nw)               :: wL,wR,wroe
    double precision, dimension(ixG^T,ndim), intent(in) :: x
    double precision, dimension(ixG^T)                  :: tmp,tmp2
    !-----------------------------------------------------------------------------


  end subroutine average2

  !> Calculate the il-th characteristic speed and the jump in the il-th 
  !> characteristic variable in the idim direction within ixL. 
  !> The eigenvalues and the L=R**(-1) matrix is calculated from wroe. 
  !> jump(il)=Sum_il L(il,iw)*(wR(iw)-wL(iw))
  subroutine hd_get_eigenjump(wL,wR,wroe,x,ix^L,il,idim,smalla,a,jump,workroe)
    use mod_global_parameters

    integer, intent(in) :: ix^L,il,idim
    double precision, dimension(ixG^T,nw):: wL,wR,wroe
    double precision, dimension(ixG^T,ndim), intent(in) :: x
    double precision, dimension(ixG^T)   :: smalla,a,jump
    double precision, dimension(ixG^T,nworkroe) :: workroe

    call geteigenjump2(wL,wR,wroe,x,ix^L,il,idim,smalla,a,jump, &
         workroe(ixG^T,1),workroe(ixG^T,2),workroe(ixG^T,3))

  end subroutine hd_get_eigenjump

  subroutine geteigenjump2(wL,wR,wroe,x,ix^L,il,idim,smalla,a,jump,&
       csound,dpperc2,dvperc)

    ! Calculate the il-th characteristic speed and the jump in the il-th 
    ! characteristic variable in the idim direction within ixL. 
    ! The eigenvalues and the L=R**(-1) matrix is calculated from wroe. 
    ! jump(il)=Sum_il L(il,iw)*(wR(iw)-wL(iw))

    use mod_global_parameters
    use mod_tvd

    integer                                             :: ix^L,il,idim,idir
    double precision, dimension(ixG^T,nw)               :: wL,wR,wroe
    double precision, dimension(ixG^T,ndim), intent(in) :: x
    double precision, dimension(ixG^T)                  :: smalla,a,jump,tmp,tmp2
    double precision, dimension(ixG^T)                  :: csound,dpperc2,dvperc
    double precision                                    :: kin_en(ixG^T)

    if(il==1)then
       !First calculate the square of the sound speed: c**2=(gamma-1)*(h-0.5*v**2)
       kin_en(ix^S) = 0.5d0 * sum(wroe(ix^S, mom(:))**2, dim=^ND+1)
       csound(ix^S)=(hd_gamma-one)*(wroe(ix^S,e_) - kin_en(ix^S))
       ! Make sure that csound**2 is positive
       csound(ix^S)=max(hd_gamma*smalldouble/wroe(ix^S,rho_),csound(ix^S))

       ! Calculate (pR-pL)/c**2
       ! To save memory we use tmp amnd tmp2 for pL and pR (hd_get_pthermal is OK)
       call hd_get_pthermal(wL,x,ixG^LL,ix^L,tmp)
       call hd_get_pthermal(wR,x,ixG^LL,ix^L,tmp2)
       dpperc2(ix^S)=(tmp2(ix^S)-tmp(ix^S))/csound(ix^S)

       !Now get the correct sound speed
       csound(ix^S)=sqrt(csound(ix^S))

       ! Calculate (vR_idim-vL_idim)/c
       dvperc(ix^S)=(wR(ix^S,mom(idim))/wR(ix^S,rho_)-&
            wL(ix^S,mom(idim))/wL(ix^S,rho_))/csound(ix^S)

    endif

    if (il == soundRW_) then
       a(ix^S)=wroe(ix^S,mom(idim))+csound(ix^S)
       jump(ix^S)=half*(dpperc2(ix^S)+wroe(ix^S,rho_)*dvperc(ix^S))
    else if (il == soundLW_) then
       a(ix^S)=wroe(ix^S,mom(idim))-csound(ix^S)
       jump(ix^S)=half*(dpperc2(ix^S)-wroe(ix^S,rho_)*dvperc(ix^S))
    else if (il == entropW_) then
       a(ix^S)=wroe(ix^S,mom(idim))
       jump(ix^S)=-dpperc2(ix^S)+wR(ix^S,rho_)-wL(ix^S,rho_)
    else
       !Determine the direction of the shear wave
       idir=il-shearW0_; if(idir>=idim)idir=idir+1
       a(ix^S)=wroe(ix^S,mom(idim))
       jump(ix^S)=wroe(ix^S,rho_)*&
            (wR(ix^S,mom(idir))/wR(ix^S,rho_)-wL(ix^S,mom(idir))/wL(ix^S,rho_))
    end if

    ! Calculate "smalla" or modify "a" based on the "typeentropy" switch
    ! Put left and right eigenvalues, if needed, into tmp and tmp2
    ! OK, since subroutines hd_get_pthermal and entropyfix do not use tmp and tmp2

    select case(typeentropy(il))
    case('yee')
       ! Based on Yee JCP 68,151 eq 3.23
       smalla(ix^S)=entropycoef(il)
    case('harten','powell')
       ! Based on Harten & Hyman JCP 50, 235 and Zeeuw & Powell JCP 104,56
       if (il == soundRW_) then
          call hd_get_pthermal(wL,x,ixG^LL,ix^L,tmp)
          tmp(ix^S)=wL(ix^S,mom(idim))/wL(ix^S,rho_)&
               + sqrt(hd_gamma*tmp(ix^S)/wL(ix^S,rho_))
          call hd_get_pthermal(wR,x,ixG^LL,ix^L,tmp2)
          tmp2(ix^S)=wR(ix^S,mom(idim))/wR(ix^S,rho_)&
               + sqrt(hd_gamma*tmp2(ix^S)/wR(ix^S,rho_))
       else if (il == soundLW_) then
          call hd_get_pthermal(wL,x,ixG^LL,ix^L,tmp)
          tmp(ix^S)=wL(ix^S,mom(idim))/wL(ix^S,rho_)&
               - sqrt(hd_gamma*tmp(ix^S)/wL(ix^S,rho_))
          call hd_get_pthermal(wR,x,ixG^LL,ix^L,tmp2)
          tmp2(ix^S)=wR(ix^S,mom(idim))/wR(ix^S,rho_)&
               - sqrt(hd_gamma*tmp2(ix^S)/wR(ix^S,rho_))
       else
          tmp(ix^S) =wL(ix^S,mom(idim))/wL(ix^S,rho_)
          tmp2(ix^S)=wR(ix^S,mom(idim))/wR(ix^S,rho_)
       end if
    end select

    call entropyfix(ix^L,il,tmp,tmp2,a,smalla)

  end subroutine geteigenjump2

  !> Multiply q by R(il,iw), where R is the right eigenvalue matrix at wroe
  subroutine hd_rtimes(q,wroe,ix^L,iw,il,idim,rq,workroe)
    use mod_global_parameters

    integer, intent(in)             :: ix^L,iw,il,idim
    double precision, intent(in)    :: wroe(ixG^T,nw)
    double precision, intent(in)    :: q(ixG^T)
    double precision, intent(inout) :: rq(ixG^T)
    double precision, intent(inout) :: workroe(ixG^T,nworkroe)
    !-----------------------------------------------------------------------------
    call rtimes2(q,wroe,ix^L,iw,il,idim,rq,workroe(ixG^T,1))
  end subroutine hd_rtimes

  subroutine rtimes2(q,wroe,ix^L,iw,il,idim,rq,csound)

    ! Multiply q by R(il,iw), where R is the right eigenvalue matrix at wroe

    use mod_global_parameters

    integer                            :: ix^L,iw,il,idim,idir
    double precision                   :: wroe(ixG^T,nw)
    double precision, dimension(ixG^T) :: q,rq,csound
    logical                            :: shearwave

    shearwave=il>shearW0_
    idir=idim
    if(shearwave)then
       ! Direction of shearwave increases with il plus idir==idim is jumped over
       idir=il-shearW0_; if(idir>=idim)idir=idir+1
    endif

    if (iw == rho_) then
       if(shearwave)then 
          rq(ix^S)=zero
       else
          rq(ix^S)=q(ix^S)
       endif
    else if (iw == e_) then
       if (il == soundRW_) then
          rq(ix^S)=q(ix^S)*(wroe(ix^S,e_)+wroe(ix^S,mom(idim))*csound(ix^S))
       else if (il == soundLW_) then
          rq(ix^S)=q(ix^S)*(wroe(ix^S,e_)-wroe(ix^S,mom(idim))*csound(ix^S))
       else if (il == entropW_) then
          rq(ix^S)=q(ix^S) * 0.5d0 * sum(wroe(ix^S, mom(:))**2, dim=^ND+1)
       else
          rq(ix^S)=q(ix^S)*wroe(ix^S,mom(idir))
       end if
    else
       if(iw==mom(idim))then
          if (il == soundRW_) then
             rq(ix^S)=q(ix^S)*(wroe(ix^S,mom(idim))+csound(ix^S))
          else if (il == soundLW_) then
             rq(ix^S)=q(ix^S)*(wroe(ix^S,mom(idim))-csound(ix^S))
          else if (il == entropW_) then
             rq(ix^S)=q(ix^S)*wroe(ix^S,mom(idim))
          else
             rq(ix^S)=zero
          end if
       else
          if(shearwave)then
             if(iw==mom(idir))then
                rq(ix^S)=q(ix^S)
             else
                rq(ix^S)=zero
             endif
          else
             rq(ix^S)=q(ix^S)*wroe(ix^S,iw)
          endif
       endif
    end if

  end subroutine rtimes2

  subroutine hd_average_iso(wL,wR,x,ix^L,idim,wroe,workroe)

    ! Calculate the Roe average of w, assignment of variables:
    ! rho -> rho, m -> v
    use mod_global_parameters

    integer, intent(in)             :: ix^L, idim
    double precision, intent(in)    :: wL(ixG^T, nw), wR(ixG^T, nw)
    double precision, intent(inout) :: wroe(ixG^T, nw)
    double precision, intent(inout) :: workroe(ixG^T, nworkroe)
    double precision, intent(in)    :: x(ixG^T, 1:^ND)

    call average2_iso(wL,wR,x,ix^L,idim,wroe,workroe(ixG^T,1))

  end subroutine hd_average_iso

  subroutine average2_iso(wL,wR,x,ix^L,idim,wroe,tmp)

    ! Calculate the Roe average of w, assignment of variables:
    ! rho -> rho, m -> v

    use mod_global_parameters

    integer:: ix^L,idim,idir
    double precision, dimension(ixG^T,nw):: wL,wR,wroe
    double precision, intent(in)    :: x(ixG^T,1:ndim)
    double precision, dimension(ixG^T):: tmp
    !-----------------------------------------------------------------------------

    select case (typeaverage)
    case ('arithmetic')
       ! This is the simple arithmetic average
       wroe(ix^S,rho_)=half*(wL(ix^S,rho_)+wR(ix^S,rho_))
       do idir = 1, ndir
          wroe(ix^S, mom(idir)) =  half * (wL(ix^S,mom(idir))/wL(ix^S,rho_) + &
               wR(ix^S,mom(idir))/wR(ix^S,rho_))
       end do
    case ('roe','default')
       ! Calculate the Roe-average
       wroe(ix^S,rho_)=sqrt(wL(ix^S,rho_)*wR(ix^S,rho_))
       ! Roe-average velocities
       tmp(ix^S)=sqrt(wL(ix^S,rho_)/wR(ix^S,rho_))
       do idir=1,ndir
          wroe(ix^S,mom(idir))=(wL(ix^S,mom(idir))/wL(ix^S,rho_)*tmp(ix^S)+&
               wR(ix^S,mom(idir))/wR(ix^S,rho_))/(one+tmp(ix^S))
       end do
    end select

  end subroutine average2_iso

  subroutine hd_get_eigenjump_iso(wL,wR,wroe,x,ix^L,il,idim,smalla,a,jump,workroe)

    ! Calculate the il-th characteristic speed and the jump in the il-th 
    ! characteristic variable in the idim direction within ixL. 
    ! The eigenvalues and the L=R**(-1) matrix is calculated from wroe. 
    ! jump(il)=Sum_il L(il,iw)*(wR(iw)-wL(iw))

    use mod_global_parameters

    integer, intent(in) :: ix^L,il,idim
    double precision, dimension(ixG^T,nw):: wL,wR,wroe
    double precision, intent(in)    :: x(ixG^T,1:ndim)
    double precision, dimension(ixG^T)   :: smalla,a,jump
    double precision, dimension(ixG^T,nworkroe) :: workroe
    !-----------------------------------------------------------------------------
    call geteigenjump2_iso(wL,wR,wroe,x,ix^L,il,idim,smalla,a,jump,workroe(ixG^T,1))

  end subroutine hd_get_eigenjump_iso

  subroutine geteigenjump2_iso(wL,wR,wroe,x,ix^L,il,idim,smalla,a,jump,csound)

    ! Calculate the il-th characteristic speed and the jump in the il-th 
    ! characteristic variable in the idim direction within ixL. 
    ! The eigenvalues and the L=R**(-1) matrix is calculated from wroe. 
    ! jump(il)=Sum_il L(il,iw)*(wR(iw)-wL(iw))

    use mod_global_parameters
    use mod_tvd

    integer:: ix^L,il,idim,idir
    double precision, dimension(ixG^T,nw):: wL,wR,wroe
    double precision, intent(in)    :: x(ixG^T,1:ndim)
    double precision, dimension(ixG^T)   :: smalla,a,jump,tmp,tmp2
    double precision, dimension(ixG^T)   :: csound
    DOUBLE PRECISION,PARAMETER:: qsmall=1.D-6
    !-----------------------------------------------------------------------------

    select case (typeaverage)
    case ('arithmetic')
       call hd_get_pthermal(wroe,x,ixG^LL,ix^L,csound)
       csound(ix^S) = sqrt(hd_gamma*csound(ix^S)/wroe(ix^S,rho_))
       ! This is the original simple Roe-solver
       if (il == soundRW_) then
          a(ix^S)=wroe(ix^S,mom(idim))+csound(ix^S)
          jump(ix^S)=half*((one-wroe(ix^S,mom(idim))/csound(ix^S))*&
               (wR(ix^S,rho_)-wL(ix^S,rho_))&
               +(wR(ix^S,mom(idim))-wL(ix^S,mom(idim)))/csound(ix^S))
       else if (il == soundLW_) then
          a(ix^S)=wroe(ix^S,mom(idim))-csound(ix^S)
          jump(ix^S)=half*((one+wroe(ix^S,mom(idim))/csound(ix^S))*&
               (wR(ix^S,rho_)-wL(ix^S,rho_))&
               -(wR(ix^S,mom(idim))-wL(ix^S,mom(idim)))/csound(ix^S))
       else
          ! Determine direction of shear wave
          idir=il-shearW0_; if(idir>=idim)idir=idir+1
          a(ix^S)=wroe(ix^S,mom(idim))
          jump(ix^S)=-wroe(ix^S,mom(idir))*(wR(ix^S,rho_)-wL(ix^S,rho_))&
               +(wR(ix^S,mom(idir))-wL(ix^S,mom(idir)))
       end if
    case ('roe','default')
       call hd_get_pthermal(wroe,x,ixG^LL,ix^L,csound)
       call hd_get_pthermal(wL,x,ixG^LL,ix^L,tmp)
       call hd_get_pthermal(wR,x,ixG^LL,ix^L,tmp2)
       where(abs(wL(ix^S,rho_)-wR(ix^S,rho_))<=qsmall*(wL(ix^S,rho_)+wR(ix^S,rho_)))
          csound(ix^S) = sqrt(hd_gamma*csound(ix^S)/wroe(ix^S,rho_))
       elsewhere
          csound(ix^S) =  sqrt(hd_gamma*(tmp2(ix^S)-tmp(ix^S))/&
               (wR(ix^S,rho_)-wL(ix^S,rho_)))
       end where
       ! This is the Roe solver by Glaister
       ! based on P. Glaister JCP 93, 477-480 (1991)
       if (il == soundRW_) then
          a(ix^S)=wroe(ix^S,mom(idim))+csound(ix^S)
          jump(ix^S)=half*((wR(ix^S,rho_)-wL(ix^S,rho_))+&
               wroe(ix^S,rho_)/csound(ix^S)*(wR(ix^S,mom(idim))/wR(ix^S,rho_)-&
               wL(ix^S,mom(idim))/wL(ix^S,rho_)))
       else if (il == soundLW_) then
          a(ix^S)=wroe(ix^S,mom(idim))-csound(ix^S)
          jump(ix^S)=half*((wR(ix^S,rho_)-wL(ix^S,rho_))-&
               wroe(ix^S,rho_)/csound(ix^S)*(wR(ix^S,mom(idim))/wR(ix^S,rho_)-&
               wL(ix^S,mom(idim))/wL(ix^S,rho_)))
       else
          ! Determine direction of shear wave
          idir=il-shearW0_; if(idir>=idim)idir=idir+1
          a(ix^S)=wroe(ix^S,mom(idim))
          jump(ix^S)=wroe(ix^S,rho_)*(wR(ix^S,mom(idir))/wR(ix^S,rho_)-&
               wL(ix^S,mom(idir))/wL(ix^S,rho_))
       end if
    end select

    ! Calculate "smalla" or modify "a" based on the "typeentropy" switch
    ! Use tmp and tmp2 for the left and right eigenvalues if needed
    select case(typeentropy(il))
    case('yee')
       ! Based on Yee JCP 68,151 eq 3.23
       smalla(ix^S)=entropycoef(il)
    case('harten','powell')
       ! Based on Harten & Hyman JCP 50, 235 and Zeeuw & Powell JCP 104,56
       if (il == soundRW_) then
          call hd_get_pthermal(wL,x,ixG^LL,ix^L,tmp)
          tmp(ix^S) =wL(ix^S,mom(idim))/wL(ix^S,rho_)&
               + sqrt(hd_gamma*tmp(ix^S)/wL(ix^S,rho_))
          call hd_get_pthermal(wR,x,ixG^LL,ix^L,tmp2)
          tmp2(ix^S)=wR(ix^S,mom(idim))/wR(ix^S,rho_)&
               + sqrt(hd_gamma*tmp2(ix^S)/wR(ix^S,rho_))
       else if (il == soundLW_) then
          call hd_get_pthermal(wL,x,ixG^LL,ix^L,tmp)
          tmp(ix^S) =wL(ix^S,mom(idim))/wL(ix^S,rho_)&
               - sqrt(hd_gamma*tmp(ix^S)/wL(ix^S,rho_))
          call hd_get_pthermal(wR,x,ixG^LL,ix^L,tmp2)
          tmp2(ix^S)=wR(ix^S,mom(idim))/wR(ix^S,rho_)&
               - sqrt(hd_gamma*tmp2(ix^S)/wR(ix^S,rho_))
       else
          tmp(ix^S) =wL(ix^S,mom(idim))/wL(ix^S,rho_)
          tmp2(ix^S)=wR(ix^S,mom(idim))/wR(ix^S,rho_)
       end if
    end select

    call entropyfix(ix^L,il,tmp,tmp2,a,smalla)

  end subroutine geteigenjump2_iso

  subroutine hd_rtimes_iso(q,wroe,ix^L,iw,il,idim,rq,workroe)

    ! Multiply q by R(il,iw), where R is the right eigenvalue matrix at wroe

    use mod_global_parameters

    integer, intent(in)             :: ix^L,iw,il,idim
    double precision, intent(in)    :: wroe(ixG^T,nw)
    double precision, intent(in)    :: q(ixG^T)
    double precision, intent(inout) :: rq(ixG^T)
    double precision, intent(inout) :: workroe(ixG^T,nworkroe)
    !-----------------------------------------------------------------------------

    call rtimes2_iso(q,wroe,ix^L,iw,il,idim,rq,workroe(ixG^T,1))

  end subroutine hd_rtimes_iso

  subroutine rtimes2_iso(q,wroe,ix^L,iw,il,idim,rq,csound)

    ! Multiply q by R(il,iw), where R is the right eigenvalue matrix at wroe

    use mod_global_parameters

    integer::          ix^L,iw,il,idim,idir
    double precision:: wroe(ixG^T,nw)
    double precision, dimension(ixG^T):: q,rq,csound

    if(iw==rho_)then
       if (il == soundRW_ .or. il == soundLW_) then
          rq(ix^S)=q(ix^S)
       else
          rq(ix^S)=zero
       end if
    else if(iw==mom(idim))then
       if (il == soundRW_) then
          rq(ix^S)=q(ix^S)*(wroe(ix^S,mom(idim))+csound(ix^S))
       else if (il == soundLW_) then
          rq(ix^S)=q(ix^S)*(wroe(ix^S,mom(idim))-csound(ix^S))
       else
          rq(ix^S)=zero
       end if
    else
       if (il == soundRW_ .or. il == soundLW_) then
          rq(ix^S)=q(ix^S)*wroe(ix^S,iw)
       else
          !Determine direction of shear wave
          idir=il-shearW0_; if(idir>=idim)idir=idir+1
          if(iw==mom(idir)) then
             rq(ix^S)=q(ix^S)
          else
             rq(ix^S)=zero
          endif
       end if
    endif

  end subroutine rtimes2_iso

end module mod_hd_roe
