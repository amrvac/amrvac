!> Module with finite difference methods for fluxes
module mod_finite_difference

  implicit none
  private

  public :: fd
  public :: centdiff
  public :: centdiff4

contains

  subroutine fd(method,qdt,ixI^L,ixO^L,idim^LIM, &
       qtC,sCT,qt,snew,sold,fC,dx^D,x)
    use mod_physics
    use mod_source, only: addsource2
    use mod_global_parameters

    character(len=*), intent(in)                                     :: method
    double precision, intent(in)                                     :: qdt, qtC, qt, dx^D
    integer, intent(in)                                              :: ixI^L, ixO^L, idim^LIM
    double precision, dimension(ixI^S,1:ndim), intent(in)            :: x

    type(state)                                                      :: sCT, snew, sold
    double precision, dimension(ixI^S,1:nwflux,1:ndim), intent(out)  :: fC

    double precision, dimension(ixI^S,1:nwflux)                      :: fCT
    double precision, dimension(ixI^S,1:nw)                          :: fm, fp, fmR, fpL, wprim
    double precision, dimension(ixI^S)                               :: v
    double precision                                                 :: dxinv(1:ndim)
    logical                                                          :: transport
    integer                                                          :: idims, iw, ixC^L, ix^L, hxO^L, ixCR^L

    associate(wCT=>sCT%w,wnew=>snew%w,wold=>sold%w)

    ^D&dxinv(^D)=-qdt/dx^D;
    do idims= idim^LIM

       block%iw0=idims

       ! Get fluxes for the whole grid (mesh+nghostcells)
       {^D& ixCmin^D = ixOmin^D - nghostcells * kr(idims,^D)\}
       {^D& ixCmax^D = ixOmax^D + nghostcells * kr(idims,^D)\}

       hxO^L=ixO^L-kr(idims,^D);
       ! ix is centered index in the idim direction from ixOmin-1/2 to ixOmax+1/2
       ixmax^D=ixOmax^D; ixmin^D=hxOmin^D;

       {^D& ixCRmin^D = ixCmin^D + kr(idims,^D)*phys_wider_stencil\}
       {^D& ixCRmax^D = ixCmax^D - kr(idims,^D)*phys_wider_stencil\}

       wprim=wCT
       call phys_to_primitive(ixG^LL,ixCR^L,wprim,x)
       call phys_get_flux(wCT,wprim,x,ixG^LL,ixCR^L,idims,fCT)

       do iw=1,nwflux
          ! Lax-Friedrich splitting:
          fp(ixCR^S,iw) = half * (fCT(ixCR^S,iw) + tvdlfeps * cmax_global * wCT(ixCR^S,iw))
          fm(ixCR^S,iw) = half * (fCT(ixCR^S,iw) - tvdlfeps * cmax_global * wCT(ixCR^S,iw))
       end do ! iw loop

       ! now do the reconstruction of fp and fm:
       call reconstructL(ixI^L,ix^L,idims,fp,fpL)
       call reconstructR(ixI^L,ix^L,idims,fm,fmR)

       do iw=1,nwflux
          if (slab_uniform) then
             fC(ix^S,iw,idims) = dxinv(idims) * (fpL(ix^S,iw) + fmR(ix^S,iw))
             wnew(ixO^S,iw)=wnew(ixO^S,iw)+ &
                  (fC(ixO^S,iw,idims)-fC(hxO^S,iw,idims))
          else
             fC(ix^S,iw,idims)=-qdt*block%surfaceC(ix^S,idims) * (fpL(ix^S,iw) + fmR(ix^S,iw))
             wnew(ixO^S,iw)=wnew(ixO^S,iw)+ &
                  (fC(ixO^S,iw,idims)-fC(hxO^S,iw,idims))/block%dvolume(ixO^S)
          end if
       end do ! iw loop

    end do !idims loop

    if (.not.slab.and.idimmin==1) call phys_add_source_geom(qdt,ixI^L,ixO^L,wCT,wnew,x)
    call addsource2(qdt*dble(idimmax-idimmin+1)/dble(ndim), &
         ixI^L,ixO^L,1,nw,qtC,wCT,qt,wnew,x,.false.)

    ! check and optionally correct unphysical values
    call phys_handle_small_values(.false.,wnew,x,ixI^L,ixO^L,'fd')
    end associate
  end subroutine fd

  subroutine reconstructL(ixI^L,iL^L,idims,w,wLC)
    use mod_global_parameters
    use mod_mp5
    use mod_limiter

    integer, intent(in)             :: ixI^L, iL^L, idims
    double precision, intent(in)    :: w(ixI^S,1:nw)

    double precision, intent(out)   :: wLC(ixI^S,1:nw) 

    double precision                :: ldw(ixI^S), dwC(ixI^S)
    integer                         :: jxR^L, ixC^L, jxC^L, kxC^L, iw

    select case (typelimiter)
    case (limiter_mp5)
       call MP5limiterL(ixI^L,iL^L,idims,w,wLC)
    case default 

       kxCmin^D=ixImin^D; kxCmax^D=ixImax^D-kr(idims,^D);

       wLC(kxC^S,1:nwflux) = w(kxC^S,1:nwflux)

       jxR^L=iL^L+kr(idims,^D);

       ixCmax^D=jxRmax^D; ixCmin^D=iLmin^D-kr(idims,^D);
       jxC^L=ixC^L+kr(idims,^D);

       do iw=1,nwflux
          dwC(ixC^S)=w(jxC^S,iw)-w(ixC^S,iw)

          call dwlimiter2(dwC,ixI^L,ixC^L,idims,typelimiter,ldw=ldw)

          wLC(iL^S,iw)=wLC(iL^S,iw)+half*ldw(iL^S)
       end do
    end select

  end subroutine reconstructL

  subroutine reconstructR(ixI^L,iL^L,idims,w,wRC)
    use mod_global_parameters
    use mod_mp5
    use mod_limiter

    integer, intent(in)             :: ixI^L, iL^L, idims
    double precision, intent(in)    :: w(ixI^S,1:nw)

    double precision, intent(out)   :: wRC(ixI^S,1:nw) 

    double precision                :: rdw(ixI^S), dwC(ixI^S)
    integer                         :: jxR^L, ixC^L, jxC^L, kxC^L, kxR^L, iw

    select case (typelimiter)
    case (limiter_mp5)
       call MP5limiterR(ixI^L,iL^L,idims,w,wRC)
    case default 

       kxCmin^D=ixImin^D; kxCmax^D=ixImax^D-kr(idims,^D);
       kxR^L=kxC^L+kr(idims,^D);

       wRC(kxC^S,1:nwflux)=w(kxR^S,1:nwflux)

       jxR^L=iL^L+kr(idims,^D);
       ixCmax^D=jxRmax^D; ixCmin^D=iLmin^D-kr(idims,^D);
       jxC^L=ixC^L+kr(idims,^D);

       do iw=1,nwflux
          dwC(ixC^S)=w(jxC^S,iw)-w(ixC^S,iw)
          call dwlimiter2(dwC,ixI^L,ixC^L,idims,typelimiter,rdw=rdw)

          wRC(iL^S,iw)=wRC(iL^S,iw)-half*rdw(jxR^S)
       end do
    end select

  end subroutine reconstructR

  subroutine centdiff(qdt,ixI^L,ixO^L,idim^LIM,qtC,sCT,qt,s,fC,dx^D,x)

    ! Advance the iws flow variables from global_time to global_time+qdt within ixO^L by centered 
    ! differencing in space the dw/dt+dF_i(w)/dx_i=S type equation. 
    ! wCT contains the time centered variables at time qtC for flux and source.
    ! w is the old value at qt on input and the new value at qt+qdt on output.
    use mod_physics
    use mod_limiter
    use mod_global_parameters
    use mod_source, only: addsource2

    integer, intent(in) :: ixI^L, ixO^L, idim^LIM
    double precision, intent(in) :: qdt, qtC, qt, dx^D
    double precision :: wprim(ixI^S,1:nw)
    type(state)      :: sCT, s
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision :: fC(ixI^S,1:nwflux,1:ndim)

    double precision :: f(ixI^S, nwflux)
    double precision :: dxinv(1:ndim)
    integer :: idims, iw, ix^L, hxO^L, ixC^L, jxC^L
    logical :: transport

    associate(wCT=>sCT%w,w=>s%w)
    ! An extra layer is needed in each direction for which fluxes are added.
    ix^L=ixO^L;
    do idims= idim^LIM
       ix^L=ix^L^LADDkr(idims,^D);
    end do
    if (ixI^L^LTix^L|.or.|.or.) then
       call mpistop("Error in CentDiff: Non-conforming input limits")
    end if

    wprim=wCT
    call phys_to_primitive(ixI^L,ixI^L,wprim,x)

    ^D&dxinv(^D)=-qdt/dx^D;

    ! Add fluxes to w
    do idims= idim^LIM
       ix^L=ixO^L^LADDkr(idims,^D); ixCmin^D=ixmin^D; ixCmax^D=ixOmax^D;
       jxC^L=ixC^L+kr(idims,^D); 
       hxO^L=ixO^L-kr(idims,^D);

       call phys_get_flux(wCT,wprim,x,ixI^L,ix^L,idims,f)

       do iw=1,nwflux
          ! Center flux to interface
          fC(ixC^S,iw,idims)=half*(f(ixC^S, iw)+f(jxC^S, iw))
          if (slab_uniform) then
             fC(ixC^S,iw,idims)=dxinv(idims)*fC(ixC^S,iw,idims)
             w(ixO^S,iw)=w(ixO^S,iw)+(fC(ixO^S,iw,idims)-fC(hxO^S,iw,idims))
          else
             fC(ixC^S,iw,idims)=-qdt*block%surfaceC(ixC^S,idims)*fC(ixC^S,iw,idims)
             w(ixO^S,iw)=w(ixO^S,iw)+ &
                  (fC(ixO^S,iw,idims)-fC(hxO^S,iw,idims))/block%dvolume(ixO^S)
          end if
       end do    !next iw
    end do       !next idims

    if (.not.slab.and.idimmin==1) call phys_add_source_geom(qdt,ixI^L,ixO^L,wCT,w,x)
    call addsource2(qdt*dble(idimmax-idimmin+1)/dble(ndim), &
         ixI^L,ixO^L,1,nw,qtC,wCT,qt,w,x,.false.)

    ! check and optionally correct unphysical values
    call phys_handle_small_values(.false.,w,x,ixI^L,ixO^L,'centdiff')
    
    end associate
  end subroutine centdiff

  subroutine centdiff4(qdt,ixI^L,ixO^L,idim^LIM,qtC,sCT,qt,s,fC,dx^D,x)

    ! Advance the flow variables from global_time to global_time+qdt within ixO^L by
    ! fourth order centered differencing in space 
    ! for the dw/dt+dF_i(w)/dx_i=S type equation.
    ! wCT contains the time centered variables at time qtC for flux and source.
    ! w is the old value at qt on input and the new value at qt+qdt on output.
    use mod_physics
    use mod_finite_volume, only: reconstruct_LR
    use mod_global_parameters
    use mod_source, only: addsource2

    integer, intent(in) :: ixI^L, ixO^L, idim^LIM
    double precision, intent(in) :: qdt, qtC, qt, dx^D
    type(state)      :: sCT, s
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision :: fC(ixI^S,1:nwflux,1:ndim)

    double precision :: v(ixI^S,ndim), f(ixI^S, nwflux)

    double precision, dimension(ixI^S,1:nw) :: wprim, wLC, wRC
    ! left and right constructed status in primitive form, needed for better performance
    double precision, dimension(ixI^S,1:nw) :: wLp, wRp
    double precision, dimension(ixI^S)      :: vLC, phi, cmaxLC, cmaxRC

    double precision :: dxinv(1:ndim)
    integer :: idims, iw, ix^L, hxO^L, ixC^L, jxC^L, hxC^L, kxC^L, kkxC^L, kkxR^L
    logical :: transport, new_cmax, patchw(ixI^S)

    associate(wCT=>sCT%w,w=>s%w)
    ! two extra layers are needed in each direction for which fluxes are added.
    ix^L=ixO^L;
    do idims= idim^LIM
       ix^L=ix^L^LADD2*kr(idims,^D);
    end do

    if (ixI^L^LTix^L|.or.|.or.) then
       call mpistop("Error in CentDiff4: Non-conforming input limits")
    end if

    ^D&dxinv(^D)=-qdt/dx^D;

    wprim=wCT
    call phys_to_primitive(ixI^L,ixI^L,wprim,x)

    ! Add fluxes to w
    do idims= idim^LIM
       block%iw0=idims

       ix^L=ixO^L^LADD2*kr(idims,^D); 
       hxO^L=ixO^L-kr(idims,^D);

       ixCmin^D=hxOmin^D; ixCmax^D=ixOmax^D;
       hxC^L=ixC^L-kr(idims,^D); 
       jxC^L=ixC^L+kr(idims,^D); 
       kxC^L=ixC^L+2*kr(idims,^D); 

       kkxCmin^D=ixImin^D; kkxCmax^D=ixImax^D-kr(idims,^D);
       kkxR^L=kkxC^L+kr(idims,^D);
       wRp(kkxC^S,1:nwflux)=wprim(kkxR^S,1:nwflux)
       wLp(kkxC^S,1:nwflux)=wprim(kkxC^S,1:nwflux)

       ! apply limited reconstruction for left and right status at cell interfaces
       call reconstruct_LR(ixI^L,ixC^L,ixC^L,idims,wprim,wLC,wRC,wLp,wRp,x)

       ! Calculate velocities from upwinded values
       call phys_get_cmax(wLC,x,ixG^LL,ixC^L,idims,cmaxLC)
       call phys_get_cmax(wRC,x,ixG^LL,ixC^L,idims,cmaxRC)
       ! now take the maximum of left and right states
       vLC(ixC^S)=max(cmaxRC(ixC^S),cmaxLC(ixC^S))

       call phys_get_flux(wCT,wprim,x,ixI^L,ix^L,idims,f)

       do iw=1,nwflux
          ! Center flux to interface
          ! f_i+1/2= (-f_(i+2) +7 f_(i+1) + 7 f_i - f_(i-1))/12
          fC(ixC^S,iw,idims)=(-f(kxC^S, iw)+7.0d0*(f(jxC^S, iw) + &
               f(ixC^S, iw))-f(hxC^S, iw))/12.0d0
          ! add rempel dissipative flux, only second order version for now
          fC(ixC^S,iw,idims)=fC(ixC^S,iw,idims)-half*vLC(ixC^S) &
               *(wRC(ixC^S,iw)-wLC(ixC^S,iw))

          if (slab_uniform) then
             fC(ixC^S,iw,idims)=dxinv(idims)*fC(ixC^S,iw,idims)
             ! result: f_(i+1/2)-f_(i-1/2) = [-f_(i+2)+8(f_(i+1)-f_(i-1))+f_(i-2)]/12
             w(ixO^S,iw)=w(ixO^S,iw)+(fC(ixO^S,iw,idims)-fC(hxO^S,iw,idims))
          else
             fC(ixC^S,iw,idims)=-qdt*block%surfaceC(ixC^S,idims)*fC(ixC^S,iw,idims)
             w(ixO^S,iw)=w(ixO^S,iw)+ &
                  (fC(ixO^S,iw,idims)-fC(hxO^S,iw,idims))/block%dvolume(ixO^S)
          end if
       end do    !next iw
    end do       !next idims

    if (.not.slab.and.idimmin==1) &
         call phys_add_source_geom(qdt,ixI^L,ixO^L,wCT,w,x)
    call addsource2(qdt*dble(idimmax-idimmin+1)/dble(ndim), &
         ixI^L,ixO^L,1,nw,qtC,wCT,qt,w,x,.false.)

    ! check and optionally correct unphysical values
    call phys_handle_small_values(.false.,w,x,ixI^L,ixO^L,'centdiff4')
    end associate
  end subroutine centdiff4

end module mod_finite_difference
