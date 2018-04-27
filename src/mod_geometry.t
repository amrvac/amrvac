!> Module with geometry-related routines (e.g., divergence, curl)
module mod_geometry
  implicit none
  public

contains

  subroutine set_pole
    use mod_global_parameters

    select case (typeaxial)
    case ("spherical") {^IFTHREED
      ! For spherical grid, check whether phi-direction is periodic
      if(periodB(ndim)) then
        if(phi_/=3) call mpistop("phi_ should be 3 in 3D spherical coord!")
        if(mod(ng3(1),2)/=0) &
             call mpistop("Number of meshes in phi-direction should be even!")
        if(abs(xprobmin2)<smalldouble) then
          if(mype==0) write(unitterm,*) &
               "Will apply pi-periodic conditions at northpole!"
          poleB(1,2)=.true.
        else
          if(mype==0) write(unitterm,*) "There is no northpole!"
        end if
        if(abs(xprobmax2-dpi)<smalldouble) then
          if(mype==0) write(unitterm,*) &
               "Will apply pi-periodic conditions at southpole!"
          poleB(2,2)=.true.
        else
          if(mype==0) write(unitterm,*) "There is no southpole!"
        end if
      end if}
    case ("cylindrical")
      {
      if (^D == phi_ .and. periodB(^D)) then
        if(mod(ng^D(1),2)/=0) then
          call mpistop("Number of meshes in phi-direction should be even!")
        end if

        if(abs(xprobmin1)<smalldouble) then
          if (mype==0) then
            write(unitterm,*) "Will apply pi-periodic conditions at r=0"
          end if
          poleB(1,1)=.true.
        else
          if (mype==0) then
            write(unitterm,*) "There is no cylindrical axis!"
          end if
        end if
      end if\}
    end select

  end subroutine set_pole

  !> Deallocate geometry-related variables
  subroutine putgridgeo(igrid)
    use mod_global_parameters
    integer, intent(in) :: igrid

    deallocate(pw(igrid)%surfaceC,pw(igrid)%surface,pw(igrid)%dvolume,pw(igrid)%dx,&
         pw(igrid)%dxcoarse,pw(igrid)%ds,pw(igrid)%dvolumecoarse)

  end subroutine putgridgeo

  subroutine fillgeo(igrid,ixG^L)
    use mod_global_parameters

    integer, intent(in) :: igrid, ixG^L

    integer :: ix^L, ixC^L
    double precision :: x(ixG^S,ndim), drs(ixG^S), dx2(ixG^S), dx3(ixG^S)
    !-----------------------------------------------------------------------------
    ix^L=ixG^L^LSUB1;

    select case (typeaxial)
    case ("slabstretch")
      drs(ixG^S)=pw(igrid)%dx(ixG^S,1)
      {^NOONED
      dx2(ixG^S)=pw(igrid)%dx(ixG^S,2)}
      {^IFTHREED
      dx3(ixG^S)=pw(igrid)%dx(ixG^S,3)}

      {^IFONED
      ixCmin^D=ixmin^D-kr(^D,1); ixCmax^D=ixmax^D;
      pw(igrid)%surfaceC(ixC^S,1)=1.d0
      pw(igrid)%surface(ixC^S,1) =1.d0
      }
      {^IFTWOD
      ixCmin^D=ixmin^D-kr(^D,1); ixCmax^D=ixmax^D;
      pw(igrid)%surfaceC(ixC^S,1)=dx2(ixC^S)
      pw(igrid)%surface(ixC^S,1) =dx2(ixC^S)
      ixCmin^D=ixmin^D-kr(^D,2); ixCmax^D=ixmax^D;
      pw(igrid)%surfaceC(ixC^S,2)=drs(ixC^S)
      pw(igrid)%surface(ixC^S,2)=drs(ixC^S)
      }
      {^IFTHREED
      ixCmin^D=ixmin^D-kr(^D,1); ixCmax^D=ixmax^D;
      pw(igrid)%surfaceC(ixC^S,1)= dx2(ixC^S)*dx3(ixC^S)
      pw(igrid)%surface(ixC^S,1)=pw(igrid)%surfaceC(ixC^S,1)
      ixCmin^D=ixmin^D-kr(^D,2); ixCmax^D=ixmax^D;
      pw(igrid)%surfaceC(ixC^S,2)= drs(ixC^S)*dx3(ixC^S)
      pw(igrid)%surface(ixC^S,2)=pw(igrid)%surfaceC(ixC^S,2)
      ixCmin^D=ixmin^D-kr(^D,3); ixCmax^D=ixmax^D;
      pw(igrid)%surfaceC(ixC^S,3)= drs(ixC^S)*dx2(ixC^S)
      pw(igrid)%surface(ixC^S,3)=pw(igrid)%surfaceC(ixC^S,3)
      }

    case ("spherical")
      x(ixG^S,1)=pw(igrid)%x(ixG^S,1)
      {^NOONED
      x(ixG^S,2)=pw(igrid)%x(ixG^S,2)}

      drs(ixG^S)=pw(igrid)%dx(ixG^S,1)
      {^NOONED
      dx2(ixG^S)=pw(igrid)%dx(ixG^S,2)}
      {^IFTHREED
      dx3(ixG^S)=pw(igrid)%dx(ixG^S,3)}

      ixCmin^D=ixmin^D-kr(^D,1); ixCmax^D=ixmax^D;

      pw(igrid)%surfaceC(ixC^S,1)=(x(ixC^S,1)+half*drs(ixC^S))**2 {^NOONED &
           *two*dsin(x(ixC^S,2))*dsin(half*dx2(ixC^S))}{^IFTHREED*dx3(ixC^S)}

      {^NOONED
      ixCmin^D=ixmin^D-kr(^D,2); ixCmax^D=ixmax^D;
      pw(igrid)%surfaceC(ixC^S,2)=x(ixC^S,1)*drs(ixC^S)&
           *dsin(x(ixC^S,2)+half*dx2(ixC^S))}{^IFTHREED*dx3(ixC^S)}

      {^IFTHREED
      ixCmin^D=ixmin^D-kr(^D,3); ixCmax^D=ixmax^D;
      pw(igrid)%surfaceC(ixC^S,3)=x(ixC^S,1)*drs(ixC^S)*dx2(ixC^S)}

      ixCmin^D=ixmin^D-kr(^D,1); ixCmax^D=ixmax^D;
      pw(igrid)%surface(ixC^S,1)=x(ixC^S,1)**2 {^NOONED &
           *two*dsin(x(ixC^S,2))*dsin(half*dx2(ixC^S))}{^IFTHREED*dx3(ixC^S)}
      {^NOONED
      ixCmin^D=ixmin^D-kr(^D,2); ixCmax^D=ixmax^D;
      pw(igrid)%surface(ixC^S,2)=x(ixC^S,1)*drs(ixC^S)&
           *dsin(x(ixC^S,2))}{^IFTHREED*dx3(ixC^S)}

      {^IFTHREED
      ixCmin^D=ixmin^D-kr(^D,3); ixCmax^D=ixmax^D;
      pw(igrid)%surface(ixC^S,3)=x(ixC^S,1)*drs(ixC^S)*dx2(ixC^S)}

    case ("cylindrical")
      x(ixG^S,1)=pw(igrid)%x(ixG^S,1)
      drs(ixG^S)=pw(igrid)%dx(ixG^S,1)
      {^NOONED
      dx2(ixG^S)=pw(igrid)%dx(ixG^S,2)}
      {^IFTHREED
      dx3(ixG^S)=pw(igrid)%dx(ixG^S,3)}

      ixCmin^D=ixmin^D-kr(^D,1); ixCmax^D=ixmax^D;
      pw(igrid)%surfaceC(ixC^S,1)=dabs(x(ixC^S,1)+half*drs(ixC^S)){^DE&*dx^DE(ixC^S) }
      {^NOONED
      ixCmin^D=ixmin^D-kr(^D,2); ixCmax^D=ixmax^D;
      if (z_==2) pw(igrid)%surfaceC(ixC^S,2)=x(ixC^S,1)*drs(ixC^S){^IFTHREED*dx3(ixC^S)}
      if (phi_==2) pw(igrid)%surfaceC(ixC^S,2)=drs(ixC^S){^IFTHREED*dx3(ixC^S)}}
      {^IFTHREED
      ixCmin^D=ixmin^D-kr(^D,3); ixCmax^D=ixmax^D;
      if (z_==3) pw(igrid)%surfaceC(ixC^S,3)=x(ixC^S,1)*drs(ixC^S)*dx2(ixC^S)
      if (phi_==3) pw(igrid)%surfaceC(ixC^S,3)=drs(ixC^S)*dx2(ixC^S)}

      ixCmin^D=ixmin^D-kr(^D,1); ixCmax^D=ixmax^D;
      pw(igrid)%surface(ixC^S,1)=dabs(x(ixC^S,1)){^DE&*dx^DE(ixC^S) }
      {^NOONED
      ixCmin^D=ixmin^D-kr(^D,2); ixCmax^D=ixmax^D;
      if (z_==2) pw(igrid)%surface(ixC^S,2)=x(ixC^S,1)*drs(ixC^S){^IFTHREED*dx3(ixC^S)}
      if (phi_==2) pw(igrid)%surface(ixC^S,2)=drs(ixC^S){^IFTHREED*dx3(ixC^S)}}
      {^IFTHREED
      ixCmin^D=ixmin^D-kr(^D,3); ixCmax^D=ixmax^D;
      if (z_==3) pw(igrid)%surface(ixC^S,3)=x(ixC^S,1)*drs(ixC^S)*dx2(ixC^S)
      if (phi_==3) pw(igrid)%surface(ixC^S,3)=drs(ixC^S)*dx2(ixC^S)}

    case default
      call mpistop("Sorry, typeaxial unknown")
    end select

  end subroutine fillgeo

  !> Calculate gradient of a scalar q within ixL in direction idir
  subroutine gradient(q,ixI^L,ixO^L,idir,gradq)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L, idir
    double precision, intent(in)    :: q(ixI^S)
    double precision, intent(inout) :: gradq(ixI^S)
    double precision                :: x(ixI^S,1:ndim)
    integer                         :: jxO^L, hxO^L

    x(ixI^S,1:ndim)=block%x(ixI^S,1:ndim)

    hxO^L=ixO^L-kr(idir,^D);
    jxO^L=ixO^L+kr(idir,^D);
    if(slab) then
      gradq(ixO^S)=half*(q(jxO^S)-q(hxO^S))/dxlevel(idir)
    else
      select case(typeaxial)
      case('slabstretch')
        gradq(ixO^S)=(q(jxO^S)-q(hxO^S))/(x(jxO^S,idir)-x(hxO^S,idir))
      case('spherical')
        select case(idir)
        case(1)
          gradq(ixO^S)=(q(jxO^S)-q(hxO^S))/((x(jxO^S,1)-x(hxO^S,1)))
          {^NOONED
        case(2)
          gradq(ixO^S)=(q(jxO^S)-q(hxO^S))/((x(jxO^S,2)-x(hxO^S,2))*x(ixO^S,1))
          }
          {^IFTHREED
        case(3)
          gradq(ixO^S)=(q(jxO^S)-q(hxO^S))/((x(jxO^S,3)-x(hxO^S,3))*x(ixO^S,1)*dsin(x(ixO^S,2)))
          }
        end select
      case('cylindrical')
        if(idir==phi_) then
          gradq(ixO^S)=(q(jxO^S)-q(hxO^S))/((x(jxO^S,phi_)-x(hxO^S,phi_))*x(ixO^S,r_))
        else
          gradq(ixO^S)=(q(jxO^S)-q(hxO^S))/(x(jxO^S,idir)-x(hxO^S,idir))
        end if
      case default
        call mpistop('Unknown geometry')
      end select
    end if

  end subroutine gradient

  !> Calculate gradient of a scalar q within ixL in direction idir
  !> first use limiter to go from cell center to edge
  subroutine gradientS(q,ixI^L,ixO^L,idir,gradq)
    use mod_global_parameters
    use mod_limiter
    use mod_ppm

    integer, intent(in)                :: ixI^L, ixO^L, idir
    double precision, intent(in)       :: q(ixI^S)
    double precision, intent(inout)    :: gradq(ixI^S)
    double precision ,dimension(ixI^S) :: qC,qL,qR,dqC,ldq,rdq

    double precision :: x(ixI^S,1:ndim)
    double precision :: invdx
    integer          :: hxO^L,ixC^L,jxC^L,gxC^L,hxC^L

    x(ixI^S,1:ndim)=block%x(ixI^S,1:ndim)

    invdx=1.d0/dxlevel(idir)
    hxO^L=ixO^L-kr(idir,^D);
    ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
    jxC^L=ixC^L+kr(idir,^D);
    gxCmin^D=ixCmin^D-kr(idir,^D);gxCmax^D=jxCmax^D;
    hxC^L=gxC^L+kr(idir,^D);

    ! set the gradient limiter here
    qR(gxC^S) = q(hxC^S)
    qL(gxC^S) = q(gxC^S)
    if (typegradlimiter/=limiter_ppm) then
      dqC(gxC^S)= qR(gxC^S)-qL(gxC^S)
      call dwlimiter2(dqC,ixI^L,gxC^L,idir,typegradlimiter,ldw=ldq,rdw=rdq)
      qL(ixC^S) = qL(ixC^S) + half*ldq(ixC^S)
      qR(ixC^S) = qR(ixC^S) - half*rdq(jxC^S)
    else
      call PPMlimitervar(ixI^L,ixO^L,idir,q,q,qL,qR)
    endif

    if(slab) then
      gradq(ixO^S)=half*(qR(ixO^S)-qL(hxO^S))*invdx
    else
      gradq(ixO^S)=(qR(ixO^S)-qL(hxO^S))/block%dx(ixO^S,idir)
      select case(typeaxial)
      case('spherical')
        select case(idir)
        case(2)
          gradq(ixO^S)=gradq(ixO^S)/x(ixO^S,1)
          {^IFTHREED
        case(3)
          gradq(ixO^S)=gradq(ixO^S)/(x(ixO^S,1)*dsin(x(ixO^S,2)))
          }
        end select
      case('cylindrical')
        if(idir==phi_) gradq(ixO^S)=gradq(ixO^S)/x(ixO^S,1)
      end select
    end if

  end subroutine gradientS

  !> Calculate divergence of a vector qvec within ixL
  subroutine divvector(qvec,ixI^L,ixO^L,divq, fourthorder)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(in)    :: qvec(ixI^S,1:ndir)
    double precision, intent(inout) :: divq(ixI^S)
    logical, intent(in), optional   :: fourthorder !< Default: false
    logical                         :: use_4th_order
    double precision                :: qC(ixI^S), invdx(1:ndim)
    integer                         :: jxO^L, hxO^L, ixC^L, jxC^L
    integer                         :: idims, ix^L, gxO^L, kxO^L

    use_4th_order = .false.
    if (present(fourthorder)) use_4th_order = fourthorder

    if (use_4th_order) then
      if (.not. slab) &
           call mpistop("divvector: 4th order only supported for slab geometry")
      ! Fourth order, stencil width is two
      ix^L=ixO^L^LADD2;
    else
      ! Second order, stencil width is one
      ix^L=ixO^L^LADD1;
    end if

    if (ixImin^D>ixmin^D.or.ixImax^D<ixmax^D|.or.) &
         call mpistop("Error in divvector: Non-conforming input limits")

    invdx=1.d0/dxlevel
    divq(ixO^S)=0.0d0

    if (slab) then
      do idims=1,ndim
        if (.not. use_4th_order) then
          ! Use second order scheme
          jxO^L=ixO^L+kr(idims,^D);
          hxO^L=ixO^L-kr(idims,^D);
          divq(ixO^S)=divq(ixO^S)+half*(qvec(jxO^S,idims) &
               - qvec(hxO^S,idims))*invdx(idims)
        else
          ! Use fourth order scheme
          kxO^L=ixO^L+2*kr(idims,^D);
          jxO^L=ixO^L+kr(idims,^D);
          hxO^L=ixO^L-kr(idims,^D);
          gxO^L=ixO^L-2*kr(idims,^D);
          divq(ixO^S)=divq(ixO^S)+&
               (-qvec(kxO^S,idims) + 8.0d0 * qvec(jxO^S,idims) - 8.0d0 * &
               qvec(hxO^S,idims) + qvec(gxO^S,idims))/(12.0d0 * dxlevel(idims))
        end if
      end do
    else
      do idims=1,ndim
        hxO^L=ixO^L-kr(idims,^D);
        ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
        jxC^L=ixC^L+kr(idims,^D);
        if(stretched_dim(idims)) then
          ! linear interpolation at cell interface along stretched dimension
          qC(ixC^S)=block%surfaceC(ixC^S,idims)*(qvec(ixC^S,idims)+0.5d0*block%dx(ixC^S,idims)*&
               (qvec(jxC^S,idims)-qvec(ixC^S,idims))/(block%x(jxC^S,idims)-block%x(ixC^S,idims)))
        else
          qC(ixC^S)=block%surfaceC(ixC^S,idims)*half*(qvec(ixC^S,idims)+qvec(jxC^S,idims))
        end if
        divq(ixO^S)=divq(ixO^S)+qC(ixO^S)-qC(hxO^S)
      end do
      divq(ixO^S)=divq(ixO^S)/block%dvolume(ixO^S)
    end if


  end subroutine divvector

  !> Calculate divergence of a vector qvec within ixL
  !> using limited extrapolation to cell edges
  subroutine divvectorS(qvec,ixI^L,ixO^L,divq)
    use mod_global_parameters
    use mod_limiter
    use mod_ppm

    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: qvec(ixI^S,1:ndir)
    double precision, intent(inout)    :: divq(ixI^S)
    double precision, dimension(ixI^S) :: qL,qR,dqC,ldq,rdq

    double precision :: invdx(1:ndim)
    integer          :: hxO^L,ixC^L,jxC^L,idims,ix^L,gxC^L,hxC^L

    ix^L=ixO^L^LADD2;

    if (ixImin^D>ixmin^D.or.ixImax^D<ixmax^D|.or.) &
         call mpistop("Error in divvectorS: Non-conforming input limits")

    invdx=1.d0/dxlevel
    divq(ixO^S)=zero
    do idims=1,ndim
      hxO^L=ixO^L-kr(idims,^D);
      ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
      jxC^L=ixC^L+kr(idims,^D);
      gxCmin^D=ixCmin^D-kr(idims,^D);gxCmax^D=jxCmax^D;
      hxC^L=gxC^L+kr(idims,^D);
      qR(gxC^S) = qvec(hxC^S,idims)
      qL(gxC^S) = qvec(gxC^S,idims)
      if(typegradlimiter/=limiter_ppm) then
        dqC(gxC^S)= qR(gxC^S)-qL(gxC^S)
        call dwlimiter2(dqC,ixI^L,gxC^L,idims,typegradlimiter,ldw=ldq,rdw=rdq)
        qL(ixC^S) = qL(ixC^S) + half*ldq(ixC^S)
        qR(ixC^S) = qR(ixC^S) - half*rdq(jxC^S)
      else
        dqC(ixI^S)=qvec(ixI^S,idims)
        call PPMlimitervar(ixI^L,ixO^L,idims,dqC,dqC,qL,qR)
      endif

      if (slab) then
        divq(ixO^S)=divq(ixO^S)+half*(qR(ixO^S)-qL(hxO^S))*invdx(idims)
      else
        qR(ixC^S)=block%surfaceC(ixC^S,idims)*qR(ixC^S)
        qL(ixC^S)=block%surfaceC(ixC^S,idims)*qL(ixC^S)
        divq(ixO^S)=divq(ixO^S)+qR(ixO^S)-qL(hxO^S)
      end if
    end do
    if(.not.slab) divq(ixO^S)=divq(ixO^S)/block%dvolume(ixO^S)

  end subroutine divvectorS

  !> Calculate curl of a vector qvec within ixL
  !> Using Stokes' theorem for non-Cartesian grids
  subroutine curlvector(qvec,ixI^L,ixO^L,curlvec,idirmin,idirmin0,ndir0)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L,ixO^L
    integer, intent(in)             :: ndir0, idirmin0
    integer, intent(inout)          :: idirmin
    double precision, intent(in)    :: qvec(ixI^S,1:ndir0)
    double precision, intent(inout) :: curlvec(ixI^S,idirmin0:3)

    integer          :: ixC^L,jxC^L,idir,jdir,kdir,hxO^L,jxO^L
    double precision :: invdx(1:ndim)
    double precision :: tmp(ixI^S),tmp2(ixI^S),xC(ixI^S)
    !-----------------------------------------------------------------------------

    ! Calculate curl within ixL: CurlV_i=eps_ijk*d_j V_k
    ! Curl can have components (idirmin0:3)
    ! Determine exact value of idirmin while doing the loop.

    idirmin=4
    curlvec(ixO^S,idirmin0:3)=zero

    if(slab) then
      invdx=1.d0/dxlevel
      do idir=idirmin0,3; do jdir=1,ndim; do kdir=1,ndir0
        if(lvc(idir,jdir,kdir)/=0)then
          tmp(ixI^S)=qvec(ixI^S,kdir)
          hxO^L=ixO^L-kr(jdir,^D);
          jxO^L=ixO^L+kr(jdir,^D);
          tmp(ixO^S)=half*(tmp(jxO^S)-tmp(hxO^S))*invdx(jdir)
          if(lvc(idir,jdir,kdir)==1)then
            curlvec(ixO^S,idir)=curlvec(ixO^S,idir)+tmp(ixO^S)
          else
            curlvec(ixO^S,idir)=curlvec(ixO^S,idir)-tmp(ixO^S)
          endif
          if(idir<idirmin)idirmin=idir
        endif
      enddo; enddo; enddo;
    else
      if(.false.) then
        do idir=idirmin0,3; do jdir=1,ndim; do kdir=1,ndir0
          if(lvc(idir,jdir,kdir)/=0)then
            tmp(ixI^S)=qvec(ixI^S,kdir)
            hxO^L=ixO^L-kr(jdir,^D);
            jxO^L=ixO^L+kr(jdir,^D);
            select case(typeaxial)
            case('slabstretch')
              tmp2(ixO^S)=(tmp(jxO^S)-tmp(hxO^S))/(block%x(jxO^S,jdir)-block%x(hxO^S,jdir))
            case('spherical')
              select case(jdir)
              case(1)
                tmp(ixI^S)=tmp(ixI^S)*block%x(ixI^S,1)
                tmp2(ixO^S)=(tmp(jxO^S)-tmp(hxO^S))/((block%x(jxO^S,1)-block%x(hxO^S,1))*block%x(ixO^S,1))
                {^NOONED    case(2)
                if(idir==1) tmp(ixI^S)=tmp(ixI^S)*dsin(block%x(ixI^S,2))
                tmp2(ixO^S)=(tmp(jxO^S)-tmp(hxO^S))/((block%x(jxO^S,2)-block%x(hxO^S,2))*block%x(ixO^S,1))
                if(idir==1) tmp2(ixO^S)=tmp2(ixO^S)/dsin(block%x(ixO^S,2))
                }
                {^IFTHREED  case(3)
                tmp2(ixO^S)=(tmp(jxO^S)-tmp(hxO^S))/((block%x(jxO^S,3)-block%x(hxO^S,3))*block%x(ixO^S,1)*dsin(block%x(ixO^S,2)))
                }
              end select
            case('cylindrical')
              if(z_==3) then
                select case(jdir)
                case(1)
                  if(idir==z_) tmp(ixI^S)=tmp(ixI^S)*block%x(ixI^S,1)
                  tmp2(ixO^S)=(tmp(jxO^S)-tmp(hxO^S))/(block%x(jxO^S,1)-block%x(hxO^S,1))
                  if(idir==z_) tmp2(ixO^S)=tmp2(ixO^S)/block%x(ixO^S,1)
                  {^NOONED      case(2)
                  tmp2(ixO^S)=(tmp(jxO^S)-tmp(hxO^S))/((block%x(jxO^S,2)-block%x(hxO^S,2))*block%x(ixO^S,1))
                  }
                  {^IFTHREED    case(3)
                  tmp2(ixO^S)=(tmp(jxO^S)-tmp(hxO^S))/(block%x(jxO^S,3)-block%x(hxO^S,3))
                  }
                end select
              end if
              if(phi_==3) then
                select case(jdir)
                case(1)
                  if(idir==z_) tmp(ixI^S)=tmp(ixI^S)*block%x(ixI^S,1)
                  tmp2(ixO^S)=(tmp(jxO^S)-tmp(hxO^S))/(block%x(jxO^S,1)-block%x(hxO^S,1))
                  if(idir==z_) tmp2(ixO^S)=tmp2(ixO^S)/block%x(ixO^S,1)
                  {^NOONED      case(2)
                  tmp2(ixO^S)=(tmp(jxO^S)-tmp(hxO^S))/(block%x(jxO^S,2)-block%x(hxO^S,2))
                  }
                  {^IFTHREED    case(3)
                  tmp2(ixO^S)=(tmp(jxO^S)-tmp(hxO^S))/((block%x(jxO^S,3)-block%x(hxO^S,3))*block%x(ixO^S,1))
                  }
                end select
              end if
            end select
            if(lvc(idir,jdir,kdir)==1)then
              curlvec(ixO^S,idir)=curlvec(ixO^S,idir)+tmp2(ixO^S)
            else
              curlvec(ixO^S,idir)=curlvec(ixO^S,idir)-tmp2(ixO^S)
            endif
            if(idir<idirmin)idirmin=idir
          endif
        enddo; enddo; enddo;
      else
        do idir=idirmin0,3; do jdir=1,ndim; do kdir=1,ndir0
          if(lvc(idir,jdir,kdir)/=0)then
            select case(typeaxial)
            case('slabstretch')
              select case(idir)
              case(1)
                if(jdir<kdir) then
                  ! idir=1,jdir=2,kdir=3
                  !! integral along 3rd dimension
                  hxO^L=ixO^L-kr(jdir,^D);
                  ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                  jxC^L=ixC^L+kr(jdir,^D);
                  ! qvec(3) at cell interface along 2nd dimension
                  tmp(ixC^S)=qvec(ixC^S,kdir)+0.5d0*block%dx(ixC^S,jdir)*&
                       (qvec(jxC^S,kdir)-qvec(ixC^S,kdir))/(block%x(jxC^S,jdir)-block%x(ixC^S,jdir))
                  curlvec(ixO^S,idir)=(tmp(ixO^S)-tmp(hxO^S))*block%dx(ixO^S,kdir)
                  !! integral along 2nd dimension
                  hxO^L=ixO^L-kr(kdir,^D);
                  ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                  jxC^L=ixC^L+kr(kdir,^D);
                  ! qvec(2) at cell interface along 3rd dimension
                  tmp(ixC^S)=qvec(ixC^S,jdir)+0.5d0*block%dx(ixC^S,kdir)*&
                       (qvec(jxC^S,jdir)-qvec(ixC^S,jdir))/(block%x(jxC^S,kdir)-block%x(ixC^S,kdir))
                  curlvec(ixO^S,idir)=(curlvec(ixO^S,idir)+(tmp(hxO^S)-tmp(ixO^S))*block%dx(ixO^S,jdir))&
                       /block%surface(ixO^S,idir)
                end if
              case(2)
                if(jdir<kdir) then
                  ! idir=2,jdir=1,kdir=3
                  !! integral along 1st dimension
                  hxO^L=ixO^L-kr(kdir,^D);
                  ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                  jxC^L=ixC^L+kr(kdir,^D);
                  ! qvec(1) at cell interface along 3rd dimension
                  tmp(ixC^S)=qvec(ixC^S,jdir)+0.5d0*block%dx(ixC^S,kdir)*&
                       (qvec(jxC^S,jdir)-qvec(ixC^S,jdir))/(block%x(jxC^S,kdir)-block%x(ixC^S,kdir))
                  curlvec(ixO^S,idir)=(tmp(ixO^S)-tmp(hxO^S))*block%dx(ixO^S,jdir)
                  !! integral along 3rd dimension
                  hxO^L=ixO^L-kr(jdir,^D);
                  ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                  jxC^L=ixC^L+kr(jdir,^D);
                  ! qvec(3) at cell interface along 1st dimension
                  tmp(ixC^S)=qvec(ixC^S,kdir)+0.5d0*block%dx(ixC^S,jdir)*&
                       (qvec(jxC^S,kdir)-qvec(ixC^S,kdir))/(block%x(jxC^S,jdir)-block%x(ixC^S,jdir))
                  curlvec(ixO^S,idir)=(curlvec(ixO^S,idir)+(tmp(hxO^S)-tmp(ixO^S))*block%dx(ixO^S,kdir))&
                       /block%surface(ixO^S,idir)
                end if
              case(3)
                if(jdir<kdir) then
                  ! idir=3,jdir=1,kdir=2
                  !! integral along 1st dimension
                  hxO^L=ixO^L-kr(kdir,^D);
                  ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                  jxC^L=ixC^L+kr(kdir,^D);
                  ! qvec(1) at cell interface along 2nd dimension
                  tmp(ixC^S)=qvec(ixC^S,jdir)+0.5d0*block%dx(ixC^S,kdir)*&
                       (qvec(jxC^S,jdir)-qvec(ixC^S,jdir))/(block%x(jxC^S,kdir)-block%x(ixC^S,kdir))
                  curlvec(ixO^S,idir)=(tmp(hxO^S)-tmp(ixO^S))*block%dx(ixO^S,jdir)
                  !! integral along 2nd dimension
                  hxO^L=ixO^L-kr(jdir,^D);
                  ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                  jxC^L=ixC^L+kr(jdir,^D);
                  ! qvec(2) at cell interface along 1st dimension
                  tmp(ixC^S)=qvec(ixC^S,kdir)+0.5d0*block%dx(ixC^S,jdir)*&
                       (qvec(jxC^S,kdir)-qvec(ixC^S,kdir))/(block%x(jxC^S,jdir)-block%x(ixC^S,jdir))
                  curlvec(ixO^S,idir)=(curlvec(ixO^S,idir)+(tmp(ixO^S)-tmp(hxO^S))*block%dx(ixO^S,kdir))&
                       /block%surface(ixO^S,idir)
                end if
              end select
            case('spherical')
              select case(idir)
              case(1)
                if(jdir<kdir) then
                  ! idir=1,jdir=2,kdir=3
                  !! integral along 3rd dimension
                  hxO^L=ixO^L-kr(jdir,^D);
                  ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                  jxC^L=ixC^L+kr(jdir,^D);
                  ! qvec(3) at cell interface along 2nd dimension
                  if(stretched_dim(jdir)) then
                    tmp(ixC^S)=qvec(ixC^S,kdir)+0.5d0*block%dx(ixC^S,jdir)*&
                         (qvec(jxC^S,kdir)-qvec(ixC^S,kdir))/(block%x(jxC^S,jdir)-block%x(ixC^S,jdir))
                  else
                    tmp(ixC^S)=0.5d0*(qvec(ixC^S,kdir)+qvec(jxC^S,kdir))
                  end if
                  ! 2nd coordinate at cell interface along 2nd dimension
                  xC(ixC^S)=block%x(ixC^S,jdir)+0.5d0*block%dx(ixC^S,jdir)
                  curlvec(ixO^S,idir)=(dsin(xC(ixO^S))*tmp(ixO^S)-&
                       dsin(xC(hxO^S))*tmp(hxO^S))*block%dx(ixO^S,kdir)
                  !! integral along 2nd dimension
                  hxO^L=ixO^L-kr(kdir,^D);
                  ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                  jxC^L=ixC^L+kr(kdir,^D);
                  ! qvec(2) at cell interface along 3rd dimension
                  if(stretched_dim(kdir)) then
                    tmp(ixC^S)=qvec(ixC^S,jdir)+0.5d0*block%dx(ixC^S,kdir)*&
                         (qvec(jxC^S,jdir)-qvec(ixC^S,jdir))/(block%x(jxC^S,kdir)-block%x(ixC^S,kdir))
                  else
                    tmp(ixC^S)=0.5d0*(qvec(ixC^S,jdir)+qvec(jxC^S,jdir))
                  end if
                  curlvec(ixO^S,idir)=(curlvec(ixO^S,idir)+(tmp(hxO^S)-tmp(ixO^S))*block%dx(ixO^S,jdir))&
                       /block%surface(ixO^S,idir)*block%x(ixO^S,idir)
                end if
              case(2)
                if(jdir<kdir) then
                  ! idir=2,jdir=1,kdir=3
                  !! integral along 1st dimension
                  hxO^L=ixO^L-kr(kdir,^D);
                  ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                  jxC^L=ixC^L+kr(kdir,^D);
                  ! qvec(1) at cell interface along 3rd dimension
                  if(stretched_dim(kdir)) then
                    tmp(ixC^S)=qvec(ixC^S,jdir)+0.5d0*block%dx(ixC^S,kdir)*&
                         (qvec(jxC^S,jdir)-qvec(ixC^S,jdir))/(block%x(jxC^S,kdir)-block%x(ixC^S,kdir))
                  else
                    tmp(ixC^S)=0.5d0*(qvec(ixC^S,jdir)+qvec(jxC^S,jdir))
                  end if
                  curlvec(ixO^S,idir)=(tmp(ixO^S)-tmp(hxO^S))*block%dx(ixO^S,1)
                  !! integral along 3rd dimension
                  hxO^L=ixO^L-kr(jdir,^D);
                  ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                  jxC^L=ixC^L+kr(jdir,^D);
                  ! qvec(3) at cell interface along 1st dimension
                  if(stretched_dim(jdir)) then
                    tmp(ixC^S)=qvec(ixC^S,kdir)+0.5d0*block%dx(ixC^S,jdir)*&
                         (qvec(jxC^S,kdir)-qvec(ixC^S,kdir))/(block%x(jxC^S,jdir)-block%x(ixC^S,jdir))
                  else
                    tmp(ixC^S)=0.5d0*(qvec(ixC^S,kdir)+qvec(jxC^S,kdir))
                  end if
                  ! 1st coordinate at cell interface along 1st dimension
                  xC(ixC^S)=block%x(ixC^S,jdir)+0.5d0*block%dx(ixC^S,jdir)
                  curlvec(ixO^S,idir)=(curlvec(ixO^S,idir)+(xC(hxO^S)*tmp(hxO^S)-xC(ixO^S)*tmp(ixO^S))*&
                       dsin(block%x(ixO^S,idir))*block%dx(ixO^S,kdir))/block%surface(ixO^S,idir)
                end if
              case(3)
                if(jdir<kdir) then
                  ! idir=3,jdir=1,kdir=2
                  !! integral along 1st dimension
                  hxO^L=ixO^L-kr(kdir,^D);
                  ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                  jxC^L=ixC^L+kr(kdir,^D);
                  ! qvec(1) at cell interface along 2nd dimension
                  if(stretched_dim(kdir)) then
                    tmp(ixC^S)=qvec(ixC^S,jdir)+0.5d0*block%dx(ixC^S,kdir)*&
                         (qvec(jxC^S,jdir)-qvec(ixC^S,jdir))/(block%x(jxC^S,kdir)-block%x(ixC^S,kdir))
                  else
                    tmp(ixC^S)=0.5d0*(qvec(ixC^S,jdir)+qvec(jxC^S,jdir))
                  end if
                  curlvec(ixO^S,idir)=(tmp(hxO^S)-tmp(ixO^S))*block%dx(ixO^S,jdir)
                  !! integral along 2nd dimension
                  hxO^L=ixO^L-kr(jdir,^D);
                  ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                  jxC^L=ixC^L+kr(jdir,^D);
                  ! qvec(2) at cell interface along 1st dimension
                  if(stretched_dim(jdir)) then
                    tmp(ixC^S)=qvec(ixC^S,kdir)+0.5d0*block%dx(ixC^S,jdir)*&
                         (qvec(jxC^S,kdir)-qvec(ixC^S,kdir))/(block%x(jxC^S,jdir)-block%x(ixC^S,jdir))
                  else
                    tmp(ixC^S)=0.5d0*(qvec(ixC^S,kdir)+qvec(jxC^S,kdir))
                  end if
                  ! 1st coordinate at cell interface along 1st dimension
                  xC(ixC^S)=block%x(ixC^S,jdir)+0.5d0*block%dx(ixC^S,jdir)
                  curlvec(ixO^S,idir)=(curlvec(ixO^S,idir)+(xC(ixO^S)*tmp(ixO^S)-xC(hxO^S)*tmp(hxO^S))*block%dx(ixO^S,kdir))&
                       /block%surface(ixO^S,idir)
                end if
              end select
            case('cylindrical')
              if(idir==r_) then
                if(jdir==phi_) then
                  ! idir=r,jdir=phi,kdir=z
                  !! integral along z dimension
                  hxO^L=ixO^L-kr(jdir,^D);
                  ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                  jxC^L=ixC^L+kr(jdir,^D);
                  ! qvec(z) at cell interface along phi dimension
                  if(stretched_dim(jdir)) then
                    tmp(ixC^S)=qvec(ixC^S,kdir)+0.5d0*block%dx(ixC^S,jdir)*&
                         (qvec(jxC^S,kdir)-qvec(ixC^S,kdir))/(block%x(jxC^S,jdir)-block%x(ixC^S,jdir))
                  else
                    tmp(ixC^S)=0.5d0*(qvec(ixC^S,kdir)+qvec(jxC^S,kdir))
                  end if
                  curlvec(ixO^S,idir)=(tmp(ixO^S)-tmp(hxO^S))*block%dx(ixO^S,kdir)
                  !! integral along phi dimension
                  hxO^L=ixO^L-kr(kdir,^D);
                  ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                  jxC^L=ixC^L+kr(kdir,^D);
                  ! qvec(phi) at cell interface along z dimension
                  if(stretched_dim(kdir)) then
                    tmp(ixC^S)=qvec(ixC^S,jdir)+0.5d0*block%dx(ixC^S,kdir)*&
                         (qvec(jxC^S,jdir)-qvec(ixC^S,jdir))/(block%x(jxC^S,kdir)-block%x(ixC^S,kdir))
                  else
                    tmp(ixC^S)=0.5d0*(qvec(ixC^S,jdir)+qvec(jxC^S,jdir))
                  end if
                  curlvec(ixO^S,idir)=(curlvec(ixO^S,idir)+(tmp(hxO^S)-tmp(ixO^S))*block%x(ixO^S,idir)*block%dx(ixO^S,jdir))&
                       /block%surface(ixO^S,idir)
                end if
              else if(idir==phi_) then
                if(jdir<kdir) then
                  ! idir=phi,jdir=r,kdir=z
                  !! integral along r dimension
                  hxO^L=ixO^L-kr(kdir,^D);
                  ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                  jxC^L=ixC^L+kr(kdir,^D);
                  ! qvec(r) at cell interface along z dimension
                  if(stretched_dim(kdir)) then
                    tmp(ixC^S)=qvec(ixC^S,jdir)+0.5d0*block%dx(ixC^S,kdir)*&
                         (qvec(jxC^S,jdir)-qvec(ixC^S,jdir))/(block%x(jxC^S,kdir)-block%x(ixC^S,kdir))
                  else
                    tmp(ixC^S)=0.5d0*(qvec(ixC^S,jdir)+qvec(jxC^S,jdir))
                  end if
                  curlvec(ixO^S,idir)=(tmp(ixO^S)-tmp(hxO^S))*block%dx(ixO^S,jdir)
                  !! integral along z dimension
                  hxO^L=ixO^L-kr(jdir,^D);
                  ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                  jxC^L=ixC^L+kr(jdir,^D);
                  ! qvec(z) at cell interface along r dimension
                  if(stretched_dim(jdir)) then
                    tmp(ixC^S)=qvec(ixC^S,kdir)+0.5d0*block%dx(ixC^S,jdir)*&
                         (qvec(jxC^S,kdir)-qvec(ixC^S,kdir))/(block%x(jxC^S,jdir)-block%x(ixC^S,jdir))
                  else
                    tmp(ixC^S)=0.5d0*(qvec(ixC^S,kdir)+qvec(jxC^S,kdir))
                  end if
                  curlvec(ixO^S,idir)=(curlvec(ixO^S,idir)+(tmp(hxO^S)-tmp(ixO^S))*block%dx(ixO^S,kdir))&
                       /block%surface(ixO^S,idir)
                end if
              else ! idir==z_
                if(jdir<kdir) then
                  ! idir=z,jdir=r,kdir=phi
                  !! integral along r dimension
                  hxO^L=ixO^L-kr(kdir,^D);
                  ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                  jxC^L=ixC^L+kr(kdir,^D);
                  ! qvec(r) at cell interface along phi dimension
                  if(stretched_dim(kdir)) then
                    tmp(ixC^S)=qvec(ixC^S,jdir)+0.5d0*block%dx(ixC^S,kdir)*&
                         (qvec(jxC^S,jdir)-qvec(ixC^S,jdir))/(block%x(jxC^S,kdir)-block%x(ixC^S,kdir))
                  else
                    tmp(ixC^S)=0.5d0*(qvec(ixC^S,jdir)+qvec(jxC^S,jdir))
                  end if
                  curlvec(ixO^S,idir)=(tmp(hxO^S)-tmp(ixO^S))*block%dx(ixO^S,jdir)
                  !! integral along phi dimension
                  hxO^L=ixO^L-kr(jdir,^D);
                  ixCmin^D=hxOmin^D;ixCmax^D=ixOmax^D;
                  jxC^L=ixC^L+kr(jdir,^D);
                  ! qvec(phi) at cell interface along r dimension
                  if(stretched_dim(jdir)) then
                    tmp(ixC^S)=qvec(ixC^S,kdir)+0.5d0*block%dx(ixC^S,jdir)*&
                         (qvec(jxC^S,kdir)-qvec(ixC^S,kdir))/(block%x(jxC^S,jdir)-block%x(ixC^S,jdir))
                  else
                    tmp(ixC^S)=0.5d0*(qvec(ixC^S,kdir)+qvec(jxC^S,kdir))
                  end if
                  ! r coordinate at cell interface along r dimension
                  xC(ixC^S)=block%x(ixC^S,jdir)+0.5d0*block%dx(ixC^S,jdir)
                  curlvec(ixO^S,idir)=(curlvec(ixO^S,idir)+(xC(ixO^S)*tmp(ixO^S)-xC(hxO^S)*tmp(hxO^S))*block%dx(ixO^S,kdir))&
                       /block%surface(ixO^S,idir)
                end if
              end if
            case default
              call mpistop("Sorry, typeaxial unknown in curlvector")
            end select
            if(idir<idirmin)idirmin=idir
          endif
        enddo; enddo; enddo;
      end if
    end if

  end subroutine curlvector

end module mod_geometry
