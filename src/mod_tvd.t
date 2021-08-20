!> Subroutines for TVD-MUSCL schemes
module mod_tvd

  implicit none
  private

  public :: tvdlimit
  public :: tvdlimit2
  public :: entropyfix

contains

  subroutine tvdlimit(method,qdt,ixI^L,ixO^L,idim^LIM,s,qt,snew,fC,dx^D,x)
    use mod_global_parameters

    integer, intent(in) :: method
    double precision, intent(in) :: qdt, qt, dx^D
    integer, intent(in) :: ixI^L, ixO^L, idim^LIM
    double precision, dimension(ixI^S,nw) :: w, wnew
    type(state) :: s, snew
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision :: fC(ixI^S,1:nwflux,1:ndim)

    integer :: idims, ixIC^L, jxIC^L
    double precision, dimension(ixI^S,nw) :: wR, wL

    associate(w=>s%w,wnew=>snew%w)
    do idims= idim^LIM
       ixICmax^D=ixOmax^D+kr(idims,^D); ixICmin^D=ixOmin^D-2*kr(idims,^D);
       wL(ixIC^S,1:nw)=w(ixIC^S,1:nw)
       jxIC^L=ixIC^L+kr(idims,^D);
       wR(ixIC^S,1:nw)=w(jxIC^S,1:nw)
       call tvdlimit2(method,qdt,ixI^L,ixIC^L,ixO^L,idims,wL,wR,wnew,x,fC,dx^D)
    end do
    end associate
  end subroutine tvdlimit

  subroutine tvdlimit2(method,qdt,ixI^L,ixIC^L,ixO^L,idims,wL,wR,wnew,x,fC,dx^D)

    ! Limit the flow variables in wnew according to typetvd. 
    ! wroeC is based on wL and wR.
    ! If method=fs_tvd an extra adtdx**2*jumpC is added to phiC for 2nd order
    ! accuracy in time.

    use mod_global_parameters
    use mod_physics_roe

    integer, intent(in) :: method
    double precision, intent(in) :: qdt, dx^D
    integer, intent(in) :: ixI^L, ixIC^L, ixO^L, idims
    double precision, dimension(ixG^T,nw) :: wL, wR
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision :: wnew(ixI^S,1:nw)
    double precision :: fC(ixI^S,1:nwflux,1:ndim)

    double precision:: workroe(ixG^T,1:nworkroe)
    double precision, dimension(ixG^T,nw) :: wroeC
    double precision, dimension(ixG^T) :: phiC, rphiC, jumpC, adtdxC, smallaC
    double precision :: dxinv(1:ndim)
    integer :: hxO^L, ixC^L, jxC^L, jxIC^L, iw, il
    !-----------------------------------------------------------------------------

    hxO^L=ixO^L-kr(idims,^D);
    ixCmax^D=ixOmax^D; ixCmin^D=hxOmin^D; 

    jxC^L=ixC^L+kr(idims,^D);
    jxIC^L=ixIC^L+kr(idims,^D);

    call phys_average(wL,wR,x,ixIC^L,idims,wroeC,workroe)

    ^D&dxinv(^D)=qdt/dx^D;

    ! A loop on characteristic variables to calculate the dissipative flux phiC.
    do il=1,nwflux
       !Calculate the jump in the il-th characteristic variable: L(wroe)*dw
       call phys_get_eigenjump(wL,wR,wroeC,x,ixIC^L,il,idims,smallaC,adtdxC,jumpC,workroe)

       ! Normalize the eigenvalue "a" (and its limit "smalla" if needed):
       if (slab_uniform) then
          adtdxC(ixIC^S)=adtdxC(ixIC^S)*dxinv(idims)
          if (typeentropy(il)=='harten' .or. typeentropy(il)=='powell')&
            smallaC(ixIC^S)=smallaC(ixIC^S)*dxinv(idims)
       else
          adtdxC(ixIC^S)=adtdxC(ixIC^S)*qdt*block%surfaceC(ixIC^S,idims)*&
             2.0d0/(block%dvolume(ixIC^S)+block%dvolume(jxIC^S))
          if (typeentropy(il)=='harten' .or. typeentropy(il)=='powell')&
            smallaC(ixIC^S)=smallaC(ixIC^S)*qdt*block%surfaceC(ixIC^S,idims)*&
             2.0d0/(block%dvolume(ixIC^S)+block%dvolume(jxIC^S))
       endif

       ! Calculate the flux limiter function phi
       call getphi(method,jumpC,adtdxC,smallaC,ixI^L,ixIC^L,ixC^L,il,idims,phiC)

       !Add R(iw,il)*phiC(il) to each variable iw in wnew
       do iw=1,nwflux
          call phys_rtimes(phiC,wroeC,ixC^L,iw,il,idims,rphiC,workroe)

          if (slab_uniform) then
             rphiC(ixC^S)=rphiC(ixC^S)*half
             fC(ixC^S,iw,idims)=fC(ixC^S,iw,idims)+rphiC(ixC^S)
             wnew(ixO^S,iw)=wnew(ixO^S,iw)+rphiC(ixO^S)-rphiC(hxO^S)
          else
             rphiC(ixC^S)=rphiC(ixC^S)*quarter* &
                   (block%dvolume(ixC^S)+block%dvolume(jxC^S))
             fC(ixC^S,iw,idims)=fC(ixC^S,iw,idims)+rphiC(ixC^S)
             wnew(ixO^S,iw)=wnew(ixO^S,iw)+(rphiC(ixO^S)-rphiC(hxO^S)) &
                                            /block%dvolume(ixO^S)
          endif
       end do  !iw
    end do     !il

  end subroutine tvdlimit2

  subroutine getphi(method,jumpC,adtdxC,smallaC,ixI^L,ixIC^L,ixC^L,il,idims,phiC)

    ! Calculate the dissipative flux from jumpC=L*dw and adtdx=eigenvalue*dt/dx.
    ! Add Lax-Wendroff type correction if method=fs_tvd.
    ! Limit according to method and typetvd.
    use mod_limiter
    use mod_global_parameters

    integer, intent(in) :: method
    integer, intent(in) :: ixI^L, ixIC^L, ixC^L, il, idims
    double precision, dimension(ixG^T) :: jumpC, adtdxC, smallaC, phiC

    double precision, dimension(ixG^T) :: ljumpC, tmp
    integer :: jxC^L, ix^L, hx^L, typelimiter
    !-----------------------------------------------------------------------------

    typelimiter=type_limiter(block%level)
    if(method==fs_tvdmu)then
       ! In the MUSCL scheme phi=|a|*jump, apply entropy fix to it
       if(typeentropy(il)=='nul'.or.typeentropy(il)=='ratio')then
          phiC(ixC^S)=abs(adtdxC(ixC^S))*jumpC(ixC^S)
       else
          where(abs(adtdxC(ixC^S))>=smallaC(ixC^S))
             phiC(ixC^S)=abs(adtdxC(ixC^S))*jumpC(ixC^S)
          elsewhere
             phiC(ixC^S)=half*(smallaC(ixC^S)+adtdxC(ixC^S)**2/smallaC(ixC^S))&
                  *jumpC(ixC^S)
          endwhere
       endif
       ! That's all for the MUSCL scheme
       return
    endif

    if(method==fs_tvd)then
       !Entropy fix to |a|-a**2
       if(typeentropy(il)=='nul'.or.typeentropy(il)=='ratio')then
          phiC(ixIC^S)=abs(adtdxC(ixIC^S))-adtdxC(ixIC^S)**2
       else
          where(abs(adtdxC(ixIC^S))>=smallaC(ixIC^S))
             phiC(ixIC^S)=abs(adtdxC(ixIC^S))-adtdxC(ixIC^S)**2
          elsewhere
             phiC(ixIC^S)=half*smallaC(ixIC^S)+&
                  (half/smallaC(ixIC^S)-one)*adtdxC(ixIC^S)**2
          endwhere
       endif
    else
       !Entropy fix to |a|
       if(typeentropy(il)=='nul'.or.typeentropy(il)=='ratio')then
          phiC(ixIC^S)=abs(adtdxC(ixIC^S))
       else
          where(abs(adtdxC(ixIC^S))>=smallaC(ixIC^S))
             phiC(ixIC^S)=abs(adtdxC(ixIC^S))
          elsewhere
             phiC(ixIC^S)=half*smallaC(ixIC^S)+&
                  half/smallaC(ixIC^S)*adtdxC(ixIC^S)**2
          endwhere
       endif
    endif

    jxC^L=ixC^L+kr(idims,^D);
    hxmin^D=ixICmin^D; hxmax^D=ixICmax^D-kr(idims,^D);
    ix^L=hx^L+kr(idims,^D);

    if (.not. limiter_symmetric(typelimiter)) then
       call mpistop("TVD only supports symmetric limiters")
    end if

    select case(typetvd)
    case('roe')
       call dwlimiter2(jumpC,ixI^L,ixIC^L,idims,typelimiter,ldw=ljumpC)
       where(adtdxC(ixC^S)<=0)
          phiC(ixC^S)=phiC(ixC^S)*(jumpC(ixC^S)-ljumpC(jxC^S))
       elsewhere
          phiC(ixC^S)=phiC(ixC^S)*(jumpC(ixC^S)-ljumpC(ixC^S))
       end where
       !extra (a*lambda)**2*delta
       if(method==fs_tvd)phiC(ixC^S)=phiC(ixC^S)+adtdxC(ixC^S)**2*jumpC(ixC^S)
    case('sweby')
       !Sweby eqs.4.11-4.15, but no 0.5 ?!
       phiC(ixIC^S)=phiC(ixIC^S)*jumpC(ixIC^S)
       call dwlimiter2(phiC,ixI^L,ixIC^L,idims,typelimiter,ldw=ljumpC)
       where(adtdxC(ixC^S)<=0)
          phiC(ixC^S)=phiC(ixC^S)-ljumpC(jxC^S)
       elsewhere
          phiC(ixC^S)=phiC(ixC^S)-ljumpC(ixC^S)
       end where
       !extra (a*lambda)**2*delta
       if(method==fs_tvd)phiC(ixC^S)=phiC(ixC^S)+adtdxC(ixC^S)**2*jumpC(ixC^S)
    case('yee')
       !eq.3.51 with correction
       call dwlimiter2(jumpC,ixI^L,ixIC^L,idims,typelimiter,ldw=ljumpC)

       !Use phiC as 0.5*(|nu|-nu**2) eq.3.45e for tvd otherwise 0.5*|nu|
       phiC(ixC^S)=half*phiC(ixC^S)
       !gamma*lambda eq.3.51d, use tmp to store agdtdxC
       where(abs(jumpC(ixC^S))>smalldouble)
          tmp(ixC^S)=adtdxC(ixC^S)+phiC(ixC^S)*&
               (ljumpC(jxC^S)-ljumpC(ixC^S))/jumpC(ixC^S)
       elsewhere
          tmp(ixC^S)=adtdxC(ixC^S)
       end where

       !eq.3.51a
       if(typeentropy(il)=='nul'.or.typeentropy(il)=='ratio')then
          phiC(ixC^S)=-phiC(ixC^S)*(ljumpC(jxC^S)+ljumpC(ixC^S))+&
               abs(tmp(ixC^S))*jumpC(ixC^S)
       else
          where(abs(tmp(ixC^S))>=smallaC(ixC^S))
             phiC(ixC^S)=-phiC(ixC^S)*(ljumpC(jxC^S)+ljumpC(ixC^S))+&
                  abs(tmp(ixC^S))*jumpC(ixC^S)
          elsewhere
             phiC(ixC^S)=-phiC(ixC^S)*(ljumpC(jxC^S)+ljumpC(ixC^S))+&
                  (half*smallaC(ixC^S)+half/smallaC(ixC^S)*&
                  tmp(ixC^S)**2)*jumpC(ixC^S)
          endwhere
       endif
    case('harten')
       !See Ryu, section 2.3
       !Use phiC as 0.5*(|nu|-nu**2)*jumpC eq.3.45b,e
       phiC(ixIC^S)=half*phiC(ixIC^S)*jumpC(ixIC^S)
       call dwlimiter2(phiC,ixI^L,ixIC^L,idims,typelimiter,ldw=ljumpC)

       !gamma*lambda eq.3.45d, use tmp as agdtdxC
       where(abs(jumpC(ixC^S))>smalldouble)
          tmp(ixC^S)=adtdxC(ixC^S)+&
               (ljumpC(jxC^S)-ljumpC(ixC^S))/jumpC(ixC^S)
       elsewhere
          tmp(ixC^S)=adtdxC(ixC^S)
       end where
       !eq.3.45a with correction
       if(typeentropy(il)=='nul'.or.typeentropy(il)=='ratio')then
          phiC(ixC^S)=-ljumpC(jxC^S)-ljumpC(ixC^S)+jumpC(ixC^S)*&
               abs(tmp(ixC^S))
       else
          where(abs(tmp(ixC^S))>=smallaC(ixC^S))
             phiC(ixC^S)=-ljumpC(jxC^S)-ljumpC(ixC^S)+jumpC(ixC^S)*&
                  abs(tmp(ixC^S))
          elsewhere
             phiC(ixC^S)=-ljumpC(jxC^S)-ljumpC(ixC^S)+jumpC(ixC^S)*&
                  (half*smallaC(ixC^S)+half/smallaC(ixC^S)*tmp(ixC^S)**2)
          endwhere
       endif
       !extra -(a*lambda)**2*delta
    case default
       call mpistop("Error in TVDLimit: Unknown TVD type")
    end select

  end subroutine getphi

  subroutine entropyfix(ix^L,il,aL,aR,a,smalla)

    ! Apply entropyfix based on typeentropy(il),aL,aR, and a
    ! Calculate "smalla" (Harten,Powell) or modify "a" (ratio)

    use mod_global_parameters

    integer, intent(in) :: ix^L, il
    double precision, dimension(ixG^T) :: aL, aR, a, smalla
    !-----------------------------------------------------------------------------

    select case(typeentropy(il))
    case('harten')
       smalla(ix^S)=max(zero,a(ix^S)-aL(ix^S),aR(ix^S)-a(ix^S))
    case('powell')
       smalla(ix^S)=max(zero,two*(aR(ix^S)-aL(ix^S)))
       !!case('ratio')
       !!   where(aL(ix^S)<zero .and. aR(ix^S)>zero)&
       !!      a(ix^S)=a(ix^S)-2*aR(ix^S)*aL(ix^S)/(aR(ix^S)-aL(ix^S))
    case('yee')
       ! This has been done in geteigenjump already
    case('nul')
       ! No entropyfix is applied
    case default
       call mpistop("No such type of entropy fix")
    end select

  end subroutine entropyfix

end module mod_tvd
