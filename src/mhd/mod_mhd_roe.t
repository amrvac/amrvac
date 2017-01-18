!> Subroutines for Roe-type Riemann solver for HD
module mod_mhd_roe
  use mod_mhd_phys
  use mod_physics_roe

  implicit none
  private

  integer, parameter :: fastRW_ = 3,fastLW_=4,slowRW_=5,slowLW_=6 ! Characteristic
  integer, parameter :: entroW_ = 8,diverW_=7,alfvRW_=1,alfvLW_=2 ! waves

  public :: mhd_roe_init

contains

  subroutine mhd_roe_init()
    use mod_global_parameters, only: entropycoef, nw
    integer :: iw

    nworkroe=15
    allocate(entropycoef(nw))

    do il = 1, nw
       select case(il)
       case(fastRW_,fastLW_,slowRW_,slowLW_)
          entropycoef(il) = 0.2d0
       case(alfvRW_,alfvLW_)
          entropycoef(il) = 0.4d0
       case default
          entropycoef(il) = -1.0d0
       end select
    end do

  end subroutine mhd_roe_init

  ! Eight-wave MHD Riemann solver. See Powell, Notes on the eigensystem, Gombosi
  ! Calculate the wroe average of primitive variables in wL and wR, assignment:
  ! rho -> sqrho, m -> v, e -> p, B_idim -> B_idim, B_idir -> beta_idir
  ! Calculate also alpha_f,alpha_s,c_f,c_s,csound2,dp,rhodv
  !
  ! wL,wR,wroe are all interface centered quantities
  subroutine average(wL,wR,x,ix^L,idim,wroe,workroe)
    use mod_global_parameters

    integer                                     :: ix^L,idim,iw
    double precision, dimension(ixG^T,nw)       :: wL,wR,wroe
    double precision, intent(in)                :: x(ixG^T,1:ndim)
    double precision, dimension(ixG^T,nworkroe) :: workroe

    call average2(wL,wR,x,ix^L,idim,wroe,workroe(ixG^T,1),workroe(ixG^T,2), &
         workroe(ixG^T,3),workroe(ixG^T,4),workroe(ixG^T,5),workroe(ixG^T,6), &
         workroe(ixG^T,7),workroe(ixG^T,8))

  end subroutine average

  ! Eight-wave MHD Riemann solver. See Powell, Notes on the eigensystem, Gombosi
  ! Calculate the wroe average of primitive variables in wL and wR, assignment:
  ! rho -> sqrho, m -> v, e -> p, B_idim -> B_idim, B_idir -> beta_idir
  ! Calculate also alpha_f,alpha_s,c_f,c_s,csound2,dp,rhodv
  !
  ! wL,wR,wroe are all interface centered quantities
  subroutine average2(wL,wR,x,ix^L,idim,wroe,cfast,cslow,afast,aslow,csound2,dp, &
       rhodv,tmp)
    use mod_global_parameters

    integer                               :: ix^L,idim,idir,jdir,iw
    double precision, dimension(ixG^T,nw) :: wL,wR,wroe
    double precision, intent(in)          :: x(ixG^T,1:ndim)
    double precision, dimension(ixG^T)    :: cfast,cslow,afast,aslow,csound2,dp, &
         rhodv,tmp

    if (ndir==1) call mpistop("MHD with d=11 is the same as HD")

    !Averaging primitive variables
    wroe(ix^S,rho_)=half*(wL(ix^S,rho_)+wR(ix^S,rho_))
    {^C&
         wroe(ix^S,v^C_)=half*(wL(ix^S,m^C_)/wL(ix^S,rho_)+wR(ix^S,m^C_)/wR(ix^S,rho_))
    wroe(ix^S,b^C_)=half*(wL(ix^S,b^C_)+wR(ix^S,b^C_))
    \}

    ! Use afast and aslow for pressures pL and pR
    call getpthermal(wL,x,ixG^LL,ix^L,afast)
    call getpthermal(wR,x,ixG^LL,ix^L,aslow)
    TODO Why energy iso gamma?
    {#IFDEF ENERGY
    wroe(ix^S,pp_)=half*(afast(ix^S)+aslow(ix^S))
    }

    {#IFDEF ISO
    dp(ix^S)=aslow(ix^S)-afast(ix^S)
    }
    {#IFDEF GAMMA
    if(useprimitive.or.eqpar(gamma_)<=zero)then
       ! dp=pR-pL
       dp(ix^S)=aslow(ix^S)-afast(ix^S)
    else
       !CONSERVATIVE dp=(g-1)*(de-v*dm+0.5*v**2*drho-0.5*d(B**2))
       dp(ix^S)=(eqpar(gamma_)-one)*(wR(ix^S,e_)-wL(ix^S,e_)&
            -(^C&wroe(ix^S,m^C_)*(wR(ix^S,m^C_)-wL(ix^S,m^C_))+)&
            +half*(^C&wroe(ix^S,m^C_)**2+)*(wR(ix^S,rho_)-wL(ix^S,rho_))&
            -half*(^C&wR(ix^S,b^C_)**2-wL(ix^S,b^C_)**2+))
    endif
    }
    !CONSERVATIVE rho*dv_idim=dm_idim-v_idim*drho
    rhodv(ix^S)=wR(ix^S,m0_+idim)-wL(ix^S,m0_+idim)-&
         wroe(ix^S,m0_+idim)*(wR(ix^S,rho_)-wL(ix^S,rho_))

    !Calculate csound2,cfast,cslow,alphafast and alphaslow

    ! get csound**2
    call getcsound2prim(wroe,x,ixG^LL,ix^L,csound2)

    ! aa=B**2/rho+a**2
    cfast(ix^S)=( ^C&wroe(ix^S,b^C_)**2+ )/wroe(ix^S,rho_)+csound2(ix^S)

    ! cs**2=0.5*(aa+dsqrt(aa**2-4*a**2*(b_i**2/rho)))
    cslow(ix^S)=half*(cfast(ix^S)-dsqrt(cfast(ix^S)**2-&
         4d0*csound2(ix^S)*wroe(ix^S,b0_+idim)**2/wroe(ix^S,rho_)))

    ! cf**2=aa-cs**2
    cfast(ix^S)=cfast(ix^S)-cslow(ix^S)

    ! alpha_f**2=(a**2-cs**2)/(cf**2-cs**2)
    afast(ix^S)=(csound2(ix^S)-cslow(ix^S))/(cfast(ix^S)-cslow(ix^S))
    afast(ix^S)=min(one,max(afast(ix^S),zero))

    ! alpha_s=dsqrt(1-alpha_f**2)
    aslow(ix^S)=dsqrt(one-afast(ix^S))

    ! alpha_f=dsqrt(alpha_f**2)
    afast(ix^S)=dsqrt(afast(ix^S))

    ! cf=dsqrt(cf**2)
    cfast(ix^S)=dsqrt(cfast(ix^S))

    ! cs=dsqrt(cs**2)
    cslow(ix^S)=dsqrt(cslow(ix^S))

    !Replace the primitive variables with more useful quantities:
    ! rho -> dsqrt(rho)
    wroe(ix^S,rho_)=dsqrt(wroe(ix^S,rho_))

    ! Avoid sgn(b_idim)==0
    where(dabs(wroe(ix^S,b0_+idim))<smalldouble)&
         wroe(ix^S,b0_+idim)=smalldouble
    ! B_idir,jdir -> beta_idir,jdir
    idir=idim+1-ndir*(idim/ndir)
    if(ndir==2)then
       where(wroe(ix^S,b0_+idir)>=zero)
          wroe(ix^S,b0_+idir)=one
       elsewhere
          wroe(ix^S,b0_+idir)=-one
       end where
    else
       !beta_j=B_j/dsqrt(B_i**2+B_j**2); beta_i=B_i/dsqrt(B_i**2+B_j**2)
       jdir=idir+1-ndir*(idir/ndir)
       tmp(ix^S)=dsqrt(wroe(ix^S,b0_+idir)**2+wroe(ix^S,b0_+jdir)**2)
       where(tmp(ix^S)>smalldouble)
          wroe(ix^S,b0_+idir)=wroe(ix^S,b0_+idir)/tmp(ix^S)
          wroe(ix^S,b0_+jdir)=wroe(ix^S,b0_+jdir)/tmp(ix^S)
       elsewhere
          wroe(ix^S,b0_+idir)=dsqrt(half)
          wroe(ix^S,b0_+jdir)=dsqrt(half)
          !!wroe(ix^S,b0_+idir)=zero
          !!wroe(ix^S,b0_+jdir)=zero
       end where
    endif

  end subroutine average2

  ! Calculate the il-th characteristic speed and the jump in the il-th
  ! characteristic variable in the idim direction within ixL.
  ! The eigenvalues and the l=r**(-1) matrix is calculated from wroe.
  ! jump(il)=Sum_il l(il,iw)*(wR(iw)-wL(iw)), where w are the conservative
  ! variables. However part of the summation is done in advance and saved into
  ! bdv,bdb,dp and dv variables. "smalla" contains a lower limit for "a" to be
  ! used in the entropy fix.
  !
  ! All the variables are centered on the cell interface, thus the 
  ! "*C" notation is omitted for sake of brevity.
  subroutine geteigenjump(wL,wR,wroe,x,ix^L,il,idim,smalla,a,jump,workroe)
    use mod_global_parameters

    integer                                     :: ix^L,il,idim
    double precision, dimension(ixG^T,nw)       :: wL,wR,wroe
    double precision, intent(in)                :: x(ixG^T,1:ndim)
    double precision, dimension(ixG^T)          :: smalla,a,jump
    double precision, dimension(ixG^T,nworkroe) :: workroe

    call geteigenjump2(wL,wR,wroe,x,ix^L,il,idim,smalla,a,jump, &
         workroe(ixG^T,1),workroe(ixG^T,2), &
         workroe(ixG^T,3),workroe(ixG^T,4),workroe(ixG^T,5),workroe(ixG^T,6), &
         workroe(ixG^T,7),workroe(ixG^T,8),workroe(ixG^T,9),workroe(ixG^T,10), &
         workroe(ixG^T,11),workroe(ixG^T,12),workroe(ixG^T,13))

  end subroutine geteigenjump

  ! Calculate the il-th characteristic speed and the jump in the il-th
  ! characteristic variable in the idim direction within ixL.
  ! The eigenvalues and the l=r**(-1) matrix is calculated from wroe.
  ! jump(il)=Sum_il l(il,iw)*(wR(iw)-wL(iw)), where w are the conservative
  ! variables. However part of the summation is done in advance and saved into
  ! bdv,bdb,dp and dv variables. "smalla" contains a lower limit for "a" to be
  ! used in the entropy fix.
  !
  ! All the variables are centered on the cell interface, thus the 
  ! "*C" notation is omitted for sake of brevity.
  subroutine geteigenjump2(wL,wR,wroe,x,ix^L,il,idim,smalla,a,jump, &
       cfast,cslow,afast,aslow,csound2,dp,rhodv,bdv,bdb,cs2L,cs2R,cs2ca2L,cs2ca2R)
    use mod_global_parameters

    integer                               :: ix^L,il,idim,idir,jdir
    double precision, dimension(ixG^T,nw) :: wL,wR,wroe
    double precision, intent(in)          :: x(ixG^T,1:ndim)
    double precision, dimension(ixG^T)    :: smalla,a,jump
    double precision, dimension(ixG^T)    :: cfast,cslow,afast,aslow,csound2,dp,rhodv
    double precision, dimension(ixG^T)    :: bdv,bdb
    double precision, dimension(ixG^T)    :: aL,aR,cs2L,cs2R,cs2ca2L,cs2ca2R

    idir=idim+1-ndir*(idim/ndir)
    jdir=idir+1-ndir*(idir/ndir)

    if(il==fastRW_)then
       !Fast and slow waves use bdv=sqrho**2*sign(bx)*(betay*dvy+betaz*dvz)
       !                        bdb=sqrho*a*          (betay*dBy+betaz*dBz)
       bdv(ix^S)=wroe(ix^S,b0_+idir)* &
            (wR(ix^S,m0_+idir)/wR(ix^S,rho_)-wL(ix^S,m0_+idir)/wL(ix^S,rho_))
       if(ndir==3)bdv(ix^S)=bdv(ix^S)+wroe(ix^S,b0_+jdir)* &
            (wR(ix^S,m0_+jdir)/wR(ix^S,rho_)-wL(ix^S,m0_+jdir)/wL(ix^S,rho_))
       bdv(ix^S)=bdv(ix^S)*sign(wroe(ix^S,rho_)**2,wroe(ix^S,b0_+idim))

       bdb(ix^S)=wroe(ix^S,b0_+idir)*(wR(ix^S,b0_+idir)-wL(ix^S,b0_+idir))
       if(ndir==3)bdb(ix^S)=bdb(ix^S)+&
            wroe(ix^S,b0_+jdir)*(wR(ix^S,b0_+jdir)-wL(ix^S,b0_+jdir))
       bdb(ix^S)=bdb(ix^S)*dsqrt(csound2(ix^S))*wroe(ix^S,rho_)
    endif

    if(il==alfvRW_)then
       !Alfven waves use      bdv=0.5*sqrho**2*      (betaz*dvy-betay*dvz)
       !                      bdb=0.5*sqrho*sign(bx)*(betaz*dBy-betay*dBz)
       bdv(ix^S)=wroe(ix^S,b0_+jdir)* &
            (wR(ix^S,m0_+idir)/wR(ix^S,rho_)-wL(ix^S,m0_+idir)/wL(ix^S,rho_)) &
            -wroe(ix^S,b0_+idir)* &
            (wR(ix^S,m0_+jdir)/wR(ix^S,rho_)-wL(ix^S,m0_+jdir)/wL(ix^S,rho_))
       bdb(ix^S)=wroe(ix^S,b0_+jdir)*(wR(ix^S,b0_+idir)-wL(ix^S,b0_+idir)) &
            -wroe(ix^S,b0_+idir)*(wR(ix^S,b0_+jdir)-wL(ix^S,b0_+jdir))
       bdv(ix^S)=bdv(ix^S)*half*wroe(ix^S,rho_)**2
       bdb(ix^S)=bdb(ix^S)*half*sign(wroe(ix^S,rho_),wroe(ix^S,b0_+idim))
    endif

    select case(il)
    case(fastRW_)
       a(ix^S)=wroe(ix^S,m0_+idim)+cfast(ix^S)
       jump(ix^S)=half/csound2(ix^S)*(&
            afast(ix^S)*(+cfast(ix^S)*rhodv(ix^S)+dp(ix^S))&
            +aslow(ix^S)*(-cslow(ix^S)*bdv(ix^S)+bdb(ix^S)))
    case(fastLW_)
       a(ix^S)=wroe(ix^S,m0_+idim)-cfast(ix^S)
       jump(ix^S)=half/csound2(ix^S)*(&
            afast(ix^S)*(-cfast(ix^S)*rhodv(ix^S)+dp(ix^S))&
            +aslow(ix^S)*(+cslow(ix^S)*bdv(ix^S)+bdb(ix^S)))
    case(slowRW_)
       a(ix^S)=wroe(ix^S,m0_+idim)+cslow(ix^S)
       jump(ix^S)=half/csound2(ix^S)*(&
            aslow(ix^S)*(+cslow(ix^S)*rhodv(ix^S)+dp(ix^S))&
            +afast(ix^S)*(+cfast(ix^S)*bdv(ix^S)-bdb(ix^S)))
    case(slowLW_)
       a(ix^S)=wroe(ix^S,m0_+idim)-cslow(ix^S)
       jump(ix^S)=half/csound2(ix^S)*(&
            aslow(ix^S)*(-cslow(ix^S)*rhodv(ix^S)+dp(ix^S))&
            +afast(ix^S)*(-cfast(ix^S)*bdv(ix^S)-bdb(ix^S)))
    case(entroW_)
       a(ix^S)=wroe(ix^S,m0_+idim)
       jump(ix^S)=wR(ix^S,rho_)-wL(ix^S,rho_)-dp(ix^S)/csound2(ix^S)
    case(diverW_)
       if(divbwave)then
          a(ix^S)=wroe(ix^S,m0_+idim)
          jump(ix^S)=wR(ix^S,b0_+idim)-wL(ix^S,b0_+idim)
       else
          a(ix^S)=zero
          jump(ix^S)=zero
       endif
    case(alfvRW_)
       a(ix^S)=wroe(ix^S,m0_+idim)+dabs(wroe(ix^S,b0_+idim))/wroe(ix^S,rho_)
       jump(ix^S)=+bdv(ix^S)-bdb(ix^S)
    case(alfvLW_)
       a(ix^S)=wroe(ix^S,m0_+idim)-dabs(wroe(ix^S,b0_+idim))/wroe(ix^S,rho_)
       jump(ix^S)=-bdv(ix^S)-bdb(ix^S)
    end select

    ! Calculate "smalla" or modify "a" based on the "typeentropy" switch

    select case(typeentropy(il))
    case('yee')
       ! Based on Yee JCP 68,151 eq 3.23
       smalla(ix^S)=entropycoef(il)
    case('harten','powell', 'ratio')
       ! Based on Harten & Hyman JCP 50, 235 and Zeeuw & Powell JCP 104,56
       ! Initialize left and right eigenvalues by velocities
       aL(ix^S)= wL(ix^S,m0_+idim)/wL(ix^S,rho_)
       aR(ix^S)= wR(ix^S,m0_+idim)/wR(ix^S,rho_)
       ! Calculate the final "aL" and "aR"
       select case(il)
       case(fastRW_)
          ! These quantities will be used for all the fast and slow waves
          ! Calculate soundspeed**2 and cs**2+ca**2.
          call getcsound2(wL,x,ixG^LL,ix^L,cs2L)
          call getcsound2(wR,x,ixG^LL,ix^L,cs2R)
          cs2ca2L(ix^S)=cs2L(ix^S)+(^C&wL(ix^S,b^C_)**2+)/wL(ix^S,rho_)
          cs2ca2R(ix^S)=cs2R(ix^S)+(^C&wR(ix^S,b^C_)**2+)/wR(ix^S,rho_)
          ! Save the discriminants into cs2L and cs2R
          cs2L(ix^S)=&
               dsqrt(cs2ca2L(ix^S)**2-4d0*cs2L(ix^S)*wL(ix^S,b0_+idim)**2/wL(ix^S,rho_))
          cs2R(ix^S)=&
               dsqrt(cs2ca2R(ix^S)**2-4d0*cs2R(ix^S)*wR(ix^S,b0_+idim)**2/wR(ix^S,rho_))

          ! The left and right eigenvalues for the fast wave going to right
          aL(ix^S)=aL(ix^S) + dsqrt(half*(cs2ca2L(ix^S) + cs2L(ix^S)))
          aR(ix^S)=aR(ix^S) + dsqrt(half*(cs2ca2R(ix^S) + cs2R(ix^S)))
       case(fastLW_)
          aL(ix^S)=aL(ix^S) - dsqrt(half*(cs2ca2L(ix^S) + cs2L(ix^S)))
          aR(ix^S)=aR(ix^S) - dsqrt(half*(cs2ca2R(ix^S) + cs2R(ix^S)))
       case(slowRW_)
          aL(ix^S)=aL(ix^S) + dsqrt(half*(cs2ca2L(ix^S) - cs2L(ix^S)))
          aR(ix^S)=aR(ix^S) + dsqrt(half*(cs2ca2R(ix^S) - cs2R(ix^S)))
       case(slowLW_)
          aL(ix^S)=aL(ix^S) - dsqrt(half*(cs2ca2L(ix^S) - cs2L(ix^S)))
          aR(ix^S)=aR(ix^S) - dsqrt(half*(cs2ca2R(ix^S) - cs2R(ix^S)))
       case(entroW_,diverW_)
          ! These propagate by the velocity
       case(alfvRW_)
          ! Store the Alfven speeds into cs2ca2L and cs2ca2R
          cs2ca2L(ix^S)=dabs(wL(ix^S,b0_+idim))/dsqrt(wL(ix^S,rho_))
          cs2ca2R(ix^S)=dabs(wR(ix^S,b0_+idim))/dsqrt(wR(ix^S,rho_))

          aL(ix^S)=aL(ix^S) + cs2ca2L(ix^S)
          aR(ix^S)=aR(ix^S) + cs2ca2R(ix^S)
       case(alfvLW_)
          aL(ix^S)=aL(ix^S) - cs2ca2L(ix^S)
          aR(ix^S)=aR(ix^S) - cs2ca2R(ix^S)
       end select
    end select

    call entropyfix(ix^L,il,aL,aR,a,smalla)

  end subroutine geteigenjump2

  ! Multiply q by R(il,iw), where R is the right eigenvalue matrix at wroe
  subroutine rtimes(q,wroe,ix^L,iw,il,idim,rq,workroe)
    use mod_global_parameters

    integer                                     :: ix^L,iw,il,idim
    double precision                            :: wroe(ixG^T,nw)
    double precision, dimension(ixG^T)          :: q,rq
    double precision, dimension(ixG^T,nworkroe) :: workroe

    call rtimes2(q,wroe,ix^L,iw,il,idim,rq,&
         workroe(ixG^T,1),workroe(ixG^T,2), &
         workroe(ixG^T,3),workroe(ixG^T,4),workroe(ixG^T,5),workroe(ixG^T,6), &
         workroe(ixG^T,7),workroe(ixG^T,14),workroe(ixG^T,15))

  end subroutine rtimes

  ! Multiply q by R(il,iw), where R is the right eigenvalue matrix at wroe
  subroutine rtimes2(q,wroe,ix^L,iw,il,idim,rq, &
       cfast,cslow,afast,aslow,csound2,dp,rhodv,bv,v2a2)
    use mod_global_parameters

    integer                            :: ix^L,iw,il,idim,idir,jdir
    double precision                   :: wroe(ixG^T,nw)
    double precision, dimension(ixG^T) :: q,rq
    double precision, dimension(ixG^T) :: cfast,cslow,afast,aslow,csound2,dp,rhodv
    double precision, dimension(ixG^T) :: bv,v2a2

    idir=idim+1-ndir*(idim/ndir)
    jdir=idir+1-ndir*(idir/ndir)

    select case(iw)
    case(rho_)
       select case(il)
       case(fastRW_,fastLW_)
          rq(ix^S)=q(ix^S)*afast(ix^S)
       case(slowRW_,slowLW_)
          rq(ix^S)=q(ix^S)*aslow(ix^S)
       case(entroW_)
          rq(ix^S)=q(ix^S)
       case(diverW_,alfvRW_,alfvLW_)
          rq(ix^S)=zero
       end select
       {#IFDEF GAMMA
    case(e_)
       if(il==fastRW_)then
          if(eqpar(gamma_)>zero)then
!!!IDEAL GAS
             ! Store 0.5*v**2+(2-gamma)/(gamma-1)*a**2
             v2a2(ix^S)=half*(^C&wroe(ix^S,m^C_)**2+)+ &
                  (two-eqpar(gamma_))/(eqpar(gamma_)-one)*csound2(ix^S)
          else
             call mpistop("Correct rTimes for NONIDEAL gas in src/mhd/roe.t")
             ! Express the partial derivative de/dp using wroe
             !v2a2(ix^S)=half*(^C&wroe(ix^S,m^C_)**2+)+ &
             !  (??dedp(ix^S)??-1)*csound2(ix^S)
          endif
          ! Store sgn(bx)*(betay*vy+betaz*vz) in bv
          bv(ix^S)=wroe(ix^S,b0_+idir)*wroe(ix^S,m0_+idir)
          if(ndir==3)bv(ix^S)=bv(ix^S)+wroe(ix^S,b0_+jdir)*wroe(ix^S,m0_+jdir)
          bv(ix^S)=bv(ix^S)*sign(one,wroe(ix^S,b0_+idim))
       else if(il==alfvRW_)then
          !Store betaz*vy-betay*vz in bv
          bv(ix^S)=(wroe(ix^S,b0_+jdir)*wroe(ix^S,m0_+idir)-&
               wroe(ix^S,b0_+idir)*wroe(ix^S,m0_+jdir))
       endif
       select case(il)
       case(fastRW_)
          rq(ix^S)=q(ix^S)*(-aslow(ix^S)*cslow(ix^S)*bv(ix^S)+afast(ix^S)*&
               (v2a2(ix^S)+cfast(ix^S)*(cfast(ix^S)+wroe(ix^S,m0_+idim))))
       case(fastLW_)
          rq(ix^S)=q(ix^S)*(+aslow(ix^S)*cslow(ix^S)*bv(ix^S)+afast(ix^S)*&
               (v2a2(ix^S)+cfast(ix^S)*(cfast(ix^S)-wroe(ix^S,m0_+idim))))
       case(slowRW_)
          rq(ix^S)=q(ix^S)*(+afast(ix^S)*cfast(ix^S)*bv(ix^S)+aslow(ix^S)*&
               (v2a2(ix^S)+cslow(ix^S)*(cslow(ix^S)+wroe(ix^S,m0_+idim))))
       case(slowLW_)
          rq(ix^S)=q(ix^S)*(-afast(ix^S)*cfast(ix^S)*bv(ix^S)+aslow(ix^S)*&
               (v2a2(ix^S)+cslow(ix^S)*(cslow(ix^S)-wroe(ix^S,m0_+idim))))
       case(entroW_)
          rq(ix^S)= q(ix^S)*half*(^C&wroe(ix^S,m^C_)**2+)
       case(diverW_)
          if(divbwave)then
             rq(ix^S)= q(ix^S)*wroe(ix^S,b0_+idim)
          else
             rq(ix^S)= zero
          endif
       case(alfvRW_)
          rq(ix^S)=+q(ix^S)*bv(ix^S)
       case(alfvLW_)
          rq(ix^S)=-q(ix^S)*bv(ix^S)
       end select
       }
    case(m^C_)
       if(iw==m0_+idim)then
          select case(il)
          case(fastRW_)
             rq(ix^S)=q(ix^S)*afast(ix^S)*(wroe(ix^S,iw)+cfast(ix^S))
          case(fastLW_)
             rq(ix^S)=q(ix^S)*afast(ix^S)*(wroe(ix^S,iw)-cfast(ix^S))
          case(slowRW_)
             rq(ix^S)=q(ix^S)*aslow(ix^S)*(wroe(ix^S,iw)+cslow(ix^S))
          case(slowLW_)
             rq(ix^S)=q(ix^S)*aslow(ix^S)*(wroe(ix^S,iw)-cslow(ix^S))
          case(entroW_)
             rq(ix^S)=q(ix^S)*wroe(ix^S,iw)
          case(diverW_,alfvLW_,alfvRW_)
             rq(ix^S)=zero
          end select
       else
          select case(il)
          case(fastRW_)
             rq(ix^S)=q(ix^S)*(afast(ix^S)*wroe(ix^S,iw)-aslow(ix^S)*&
                  cslow(ix^S)*wroe(ix^S,b0_-m0_+iw)*sign(one,wroe(ix^S,b0_+idim)))
          case(fastLW_)
             rq(ix^S)=q(ix^S)*(afast(ix^S)*wroe(ix^S,iw)+aslow(ix^S)*&
                  cslow(ix^S)*wroe(ix^S,b0_-m0_+iw)*sign(one,wroe(ix^S,b0_+idim)))
          case(slowRW_)
             rq(ix^S)=q(ix^S)*(aslow(ix^S)*wroe(ix^S,iw)+afast(ix^S)*&
                  cfast(ix^S)*wroe(ix^S,b0_-m0_+iw)*sign(one,wroe(ix^S,b0_+idim)))
          case(slowLW_)
             rq(ix^S)=q(ix^S)*(aslow(ix^S)*wroe(ix^S,iw)-afast(ix^S)*&
                  cfast(ix^S)*wroe(ix^S,b0_-m0_+iw)*sign(one,wroe(ix^S,b0_+idim)))
          case(entroW_)
             rq(ix^S)=q(ix^S)*wroe(ix^S,iw)
          case(diverW_)
             rq(ix^S)=zero
          case(alfvRW_)
             if(iw==m0_+idir)then
                rq(ix^S)=+q(ix^S)*wroe(ix^S,b0_+jdir)
             else
                rq(ix^S)=-q(ix^S)*wroe(ix^S,b0_+idir)
             endif
          case(alfvLW_)
             if(iw==m0_+idir)then
                rq(ix^S)=-q(ix^S)*wroe(ix^S,b0_+jdir)
             else
                rq(ix^S)=+q(ix^S)*wroe(ix^S,b0_+idir)
             endif
          end select
       end if ! iw=m_idir,m_jdir
    case(b^C_)
       if(iw==b0_+idim)then
          if(il==diverW_ .and. divbwave)then
             rq(ix^S)=q(ix^S)
          else
             rq(ix^S)=zero
          endif
       else
          select case(il)
          case(fastRW_,fastLW_)
             rq(ix^S)=+q(ix^S)*aslow(ix^S)*dsqrt(csound2(ix^S))*wroe(ix^S,iw)&
                  /wroe(ix^S,rho_)
          case(slowRW_,slowLW_)
             rq(ix^S)=-q(ix^S)*afast(ix^S)*dsqrt(csound2(ix^S))*wroe(ix^S,iw)&
                  /wroe(ix^S,rho_)
          case(entroW_,diverW_)
             rq(ix^S)=zero
          case(alfvRW_,alfvLW_)
             if(iw==b0_+idir)then
                rq(ix^S)=-q(ix^S)*wroe(ix^S,b0_+jdir)&
                     /sign(wroe(ix^S,rho_),wroe(ix^S,b0_+idim))
             else
                rq(ix^S)=+q(ix^S)*wroe(ix^S,b0_+idir)&
                     /sign(wroe(ix^S,rho_),wroe(ix^S,b0_+idim))
             end if
          end select
       end if ! iw=b_idir,b_jdir
    end select

  end subroutine rtimes2

end module mod_mhd_roe
