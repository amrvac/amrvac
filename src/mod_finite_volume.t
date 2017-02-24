module mod_finite_volume
  implicit none
  private

  public :: finite_volume 
  public :: hancock
  public :: upwindLR

contains

  !> The non-conservative Hancock predictor for TVDLF
  !>
  !> on entry:
  !> input available on ixI^L=ixG^L asks for output on ixO^L=ixG^L^LSUBnghostcells
  !> one entry: (predictor): wCT -- w_n        wnew -- w_n   qdt=dt/2
  !> on exit :  (predictor): wCT -- w_n        wnew -- w_n+1/2
  subroutine hancock(qdt,ixI^L,ixO^L,idim^LIM,qtC,wCT,qt,wnew,dx^D,x)
    use mod_physics
    use mod_global_parameters
    use mod_source, only: addsource2

    integer, intent(in) :: ixI^L, ixO^L, idim^LIM
    double precision, intent(in) :: qdt, qtC, qt, dx^D, x(ixI^S,1:ndim)
    double precision, intent(inout) :: wCT(ixI^S,1:nw), wnew(ixI^S,1:nw)

    double precision, dimension(ixI^S,1:nw) :: wprim, wLC, wRC
    double precision :: fLC(ixI^S, nwflux), fRC(ixI^S, nwflux)
    double precision :: dxinv(1:ndim),dxdim(1:ndim)
    integer :: idims, iw, ix^L, hxO^L

    ! Expand limits in each idims direction in which fluxes are added
    ix^L=ixO^L;
    do idims= idim^LIM
       ix^L=ix^L^LADDkr(idims,^D);
    end do
    if (ixI^L^LTix^L|.or.|.or.) &
         call mpistop("Error in Hancock: Nonconforming input limits")

    wprim=wCT
    call phys_to_primitive(ixI^L,ixI^L,wprim,x)

    ^D&dxinv(^D)=-qdt/dx^D;
    ^D&dxdim(^D)=dx^D;
    do idims= idim^LIM
       block%iw0=idims
       ! Calculate w_j+g_j/2 and w_j-g_j/2
       ! First copy all variables, then upwind wLC and wRC.
       ! wLC is to the left of ixO, wRC is to the right of wCT.
       hxO^L=ixO^L-kr(idims,^D);

       wRC(hxO^S,1:nwflux)=wprim(ixO^S,1:nwflux)
       wLC(ixO^S,1:nwflux)=wprim(ixO^S,1:nwflux)

       call upwindLR(ixI^L,ixO^L,hxO^L,idims,wprim,wprim,wLC,wRC,x,.false.,dxdim(idims))

       ! Calculate the fLC and fRC fluxes
       call phys_get_flux(wRC,x,ixI^L,hxO^L,idims,fRC)
       call phys_get_flux(wLC,x,ixI^L,ixO^L,idims,fLC)

       ! Advect w(iw)
       do iw=1,nwflux
          if (slab) then
             wnew(ixO^S,iw)=wnew(ixO^S,iw)+dxinv(idims)* &
                  (fLC(ixO^S, iw)-fRC(hxO^S, iw))
          else
             select case (idims)
                {case (^D)
                wnew(ixO^S,iw)=wnew(ixO^S,iw)-qdt/block%dvolume(ixO^S) &
                     *(block%surfaceC^D(ixO^S)*fLC(ixO^S, iw) &
                     -block%surfaceC^D(hxO^S)*fRC(hxO^S, iw))\}
             end select
          end if
       end do
       block%iw0=0
    end do ! next idims

    if (.not.slab.and.idimmin==1) call phys_add_source_geom(qdt,ixI^L,ixO^L,wCT,wnew,x)
    call addsource2(qdt*dble(idimmax-idimmin+1)/dble(ndim), &
         ixI^L,ixO^L,1,nw,qtC,wCT,qt,wnew,x,.false.)

  end subroutine hancock

  !> finite volume method
  subroutine finite_volume(method,qdt,ixI^L,ixO^L,idim^LIM, &
       qtC,wCT,qt,wnew,wold,fC,dx^D,x)

    use mod_physics
    use mod_global_parameters
    use mod_tvd, only:tvdlimit2
    use mod_source, only: addsource2

    character(len=*), intent(in)                         :: method
    double precision, intent(in)                         :: qdt, qtC, qt, dx^D
    integer, intent(in)                                  :: ixI^L, ixO^L, idim^LIM
    double precision, dimension(ixI^S,1:ndim), intent(in) ::  x
    double precision, dimension(ixI^S,1:ndim)             :: xi
    double precision, dimension(ixI^S,1:nw)               :: wCT, wnew, wold
    double precision, dimension(ixI^S,1:nwflux,1:ndim)  :: fC

    double precision, dimension(ixI^S,1:nw) :: wprim, wLC, wRC, wmean
    double precision, dimension(ixI^S, nwflux) :: fLC, fRC
    double precision, dimension(ixI^S)      :: cmaxC, cmaxRC, cmaxLC
    double precision, dimension(ixI^S)      :: cminC, cminRC, cminLC
    double precision, dimension(1:ndim)     :: dxinv, dxdim
    integer, dimension(ixI^S)               :: patchf
    integer :: idims, iw, ix^L, hxO^L, ixC^L, ixCR^L, jxC^L, kxC^L, kxR^L
    logical :: CmaxMeanState

    CmaxMeanState = (typetvdlf=='cmaxmean')

    if (idimmax>idimmin .and. typelimited=='original')&
         call mpistop("Error in fv: Unsplit dim. and original is limited")


    ! The flux calculation contracts by one in the idim direction it is applied.
    ! The limiter contracts the same directions by one more, so expand ixO by 2.
    ix^L=ixO^L;
    do idims= idim^LIM
       ix^L=ix^L^LADD2*kr(idims,^D);
    end do
    if (ixI^L^LTix^L|.or.|.or.) &
         call mpistop("Error in fv : Nonconforming input limits")

    wprim=wCT
    call phys_to_primitive(ixI^L,ixI^L,wprim,x)

    ^D&dxinv(^D)=-qdt/dx^D;
    ^D&dxdim(^D)=dx^D;
    do idims= idim^LIM
       block%iw0=idims

       hxO^L=ixO^L-kr(idims,^D);
       ! ixC is centered index in the idim direction from ixOmin-1/2 to ixOmax+1/2
       ixCmax^D=ixOmax^D; ixCmin^D=hxOmin^D;

       kxCmin^D=ixImin^D; kxCmax^D=ixImax^D-kr(idims,^D);
       kxR^L=kxC^L+kr(idims,^D);

       wRC(kxC^S,1:nwflux)=wprim(kxR^S,1:nwflux)
       wLC(kxC^S,1:nwflux)=wprim(kxC^S,1:nwflux)

       ! Determine stencil size
       {ixCRmin^D = ixCmin^D - phys_wider_stencil\}
       {ixCRmax^D = ixCmax^D + phys_wider_stencil\}

       ! for second order scheme: apply limiting
       select case (typelimited)
       case ('previous')
          call upwindLR(ixI^L,ixCR^L,ixCR^L,idims,wold,wprim,wLC,wRC,x,.true.,dxdim(idims))
       case ('predictor')
          call upwindLR(ixI^L,ixCR^L,ixCR^L,idims,wprim,wprim,wLC,wRC,x,.false.,dxdim(idims))
       case ('original')
          call upwindLR(ixI^L,ixCR^L,ixCR^L,idims,wnew,wprim,wLC,wRC,x,.true.,dxdim(idims))
       case default
          call mpistop("Error in reconstruction: no such base for limiter")
       end select

       ! For the high order scheme the limiter is based on
       ! the maximum eigenvalue, it is calculated in advance.
       if (CmaxMeanState) then
          ! determine mean state and store in wprim
          wmean(ixC^S,1:nwflux)=half*(wLC(ixC^S,1:nwflux)+wRC(ixC^S,1:nwflux))
          if(method=='tvdlf'.or.method=='tvdmu') then
            call phys_get_cmax(wmean,x,ixI^L,ixC^L,idims,cmaxC)
          else
            call phys_get_cmax(wmean,x,ixI^L,ixC^L,idims,cmaxC,cminC)
          end if
       else
          ! now take the maximum of left and right states
          ! S.F. Davis, SIAM J. Sci. Statist. Comput. 1988, 9, 445
          if(method=='tvdlf'.or.method=='tvdmu') then
            call phys_get_cmax(wLC,x,ixI^L,ixC^L,idims,cmaxLC)
            call phys_get_cmax(wRC,x,ixI^L,ixC^L,idims,cmaxRC)
            cmaxC(ixC^S)=max(cmaxRC(ixC^S),cmaxLC(ixC^S))
          else
            call phys_get_cmax(wLC,x,ixI^L,ixC^L,idims,cmaxLC,cminLC)
            call phys_get_cmax(wRC,x,ixI^L,ixC^L,idims,cmaxRC,cminRC)
            cmaxC(ixC^S)=max(cmaxRC(ixC^S),cmaxLC(ixC^S))
            cminC(ixC^S)=min(cminRC(ixC^S),cminLC(ixC^S))
          end if
       end if

       call phys_modify_wLR(wLC, wRC, ixI^L, ixC^L, idims)

       call phys_get_flux(wLC,x,ixI^L,ixC^L,idims,fLC)
       call phys_get_flux(wRC,x,ixI^L,ixC^L,idims,fRC)

       ! use approximate Riemann solver to get flux at interfaces
       call get_Riemann_flux()

       block%iw0=0

    end do ! Next idims


    do idims= idim^LIM
       hxO^L=ixO^L-kr(idims,^D);
       do iw=1,nwflux

          ! Multiply the fluxes by -dt/dx since Flux fixing expects this
          if (slab) then
             fC(ixI^S,iw,idims)=dxinv(idims)*fC(ixI^S,iw,idims)
             wnew(ixO^S,iw)=wnew(ixO^S,iw) &
                  + (fC(ixO^S,iw,idims)-fC(hxO^S,iw,idims))
          else
             select case (idims)
                {case (^D)
                fC(ixI^S,iw,^D)=-qdt*fC(ixI^S,iw,idims)
                wnew(ixO^S,iw)=wnew(ixO^S,iw) &
                     + (fC(ixO^S,iw,^D)-fC(hxO^S,iw,^D))/block%dvolume(ixO^S)\}
             end select
          end if

       end do ! Next iw
       ! For the MUSCL scheme apply the characteristic based limiter
       if (method=='tvdmu') &
            call tvdlimit2(method,qdt,ixI^L,ixC^L,ixO^L,idims,wLC,wRC,wnew,x,fC,dx^D)

    end do ! Next idims

    if (.not.slab.and.idimmin==1) &
         call phys_add_source_geom(qdt,ixI^L,ixO^L,wCT,wnew,x)

    call addsource2(qdt*dble(idimmax-idimmin+1)/dble(ndim), &
         ixI^L,ixO^L,1,nw,qtC,wCT,qt,wnew,x,.false.)

  contains

    subroutine get_Riemann_flux()
      implicit none
      !=== specific to HLLC and HLLCD ===!
      double precision, dimension(ixI^S,1:nwflux)     :: whll, Fhll, fCD
      double precision, dimension(ixI^S)              :: lambdaCD
      select case(method)
      case('tvdmu')
        do iw=1,nwflux
           ! To save memory we use fLC to store (F_L+F_R)/2=half*(fLC+fRC)
           fLC(ixC^S, iw)=half*(fLC(ixC^S, iw)+fRC(ixC^S, iw))
           if (slab) then
              fC(ixC^S,iw,idims)=fLC(ixC^S, iw)
           else
              select case (idims)
                 {case (^D)
                 fC(ixC^S,iw,^D)=block%surfaceC^D(ixC^S)*fLC(ixC^S, iw)\}
              end select
           end if
        end do
      case('tvdlf')
        ! Calculate fLC=f(uL_j+1/2) and fRC=f(uR_j+1/2) for each iw
        do iw=1,nwflux

           ! To save memory we use fLC to store (F_L+F_R)/2=half*(fLC+fRC)
           fLC(ixC^S, iw)=half*(fLC(ixC^S, iw)+fRC(ixC^S, iw))

           ! Add TVDLF dissipation to the flux
           if (flux_type(idims, iw) == flux_no_dissipation) then
              fRC(ixC^S, iw)=0.d0
           else
              ! To save memory we use fRC to store -cmax*half*(w_R-w_L)
              fRC(ixC^S, iw)=-tvdlfeps*cmaxC(ixC^S)*half*(wRC(ixC^S,iw)-wLC(ixC^S,iw))
           end if

           ! fLC contains physical+dissipative fluxes
           fLC(ixC^S, iw)=fLC(ixC^S, iw)+fRC(ixC^S, iw)

           if (slab) then
              fC(ixC^S,iw,idims)=fLC(ixC^S, iw)
           else
              select case (idims)
                 {case (^D)
                 fC(ixC^S,iw,^D)=block%surfaceC^D(ixC^S)*fLC(ixC^S, iw)\}
              end select
           end if

        end do ! Next iw
      case('hll')
         patchf(ixC^S) =  1
         where(cminC(ixC^S) >= zero)
            patchf(ixC^S) = -2
         elsewhere(cmaxC(ixC^S) <= zero)
            patchf(ixC^S) =  2
         endwhere

         ! Calculate fLC=f(uL_j+1/2) and fRC=f(uR_j+1/2) for each iw
         do iw=1,nwflux
            if (flux_type(idims, iw) == flux_tvdlf) then
               fLC(ixC^S, iw) = half*((fLC(ixC^S, iw) + fRC(ixC^S, iw)) &
                    -tvdlfeps*max(cmaxC(ixC^S), dabs(cminC(ixC^S))) * &
                    (wRC(ixC^S,iw)-wLC(ixC^S,iw)))
            else
               where(patchf(ixC^S)==1)
                  ! Add hll dissipation to the flux
                  fLC(ixC^S, iw) = (cmaxC(ixC^S)*fLC(ixC^S, iw)-cminC(ixC^S) * fRC(ixC^S, iw) &
                       +tvdlfeps*cminC(ixC^S)*cmaxC(ixC^S)*(wRC(ixC^S,iw)-wLC(ixC^S,iw)))&
                       /(cmaxC(ixC^S)-cminC(ixC^S))
               elsewhere(patchf(ixC^S)== 2)
                  fLC(ixC^S, iw)=fRC(ixC^S, iw)
               elsewhere(patchf(ixC^S)==-2)
                  fLC(ixC^S, iw)=fLC(ixC^S, iw)
               endwhere
            endif

            if (slab) then
               fC(ixC^S,iw,idims)=fLC(ixC^S, iw)
            else
               select case (idims)
                  {case (^D)
                  fC(ixC^S,iw,^D)=block%surfaceC^D(ixC^S)*fLC(ixC^S, iw)\}
               end select
            end if

         end do ! Next iw
      case('hllc','hllcd')
        patchf(ixC^S) =  1
        where(cminC(ixC^S) >= zero)
           patchf(ixC^S) = -2
        elsewhere(cmaxC(ixC^S) <= zero)
           patchf(ixC^S) =  2
        endwhere
        ! Use more diffusive scheme, is actually TVDLF and selected by patchf=4
        if(method=='hllcd') &
             call phys_diffuse_hllcd(ixI^L,ixC^L,idims,wLC,wRC,fLC,fRC,patchf)

        !---- calculate speed lambda at CD ----!
        if(any(patchf(ixC^S)==1)) &
             call phys_get_lCD(wLC,wRC,fLC,fRC,cminC,cmaxC,idims,ixI^L,ixC^L, &
             whll,Fhll,lambdaCD,patchf)

        ! now patchf may be -1 or 1 due to phys_get_lCD
        if(any(abs(patchf(ixC^S))== 1))then
           !======== flux at intermediate state ========!
           call phys_get_wCD(wLC,wRC,whll,fRC,fLC,Fhll,patchf,lambdaCD,&
                cminC,cmaxC,ixI^L,ixC^L,idims,fCD)
        endif ! Calculate the CD flux

        do iw=1,nwflux
           if (flux_type(idims, iw) == flux_tvdlf) then
              fLC(ixC^S,iw) = 0.5d0 * (fLC(ixC^S,iw) + fRC(ixC^S,iw) - tvdlfeps * &
                   max(cmaxC(ixC^S), abs(cminC(ixC^S))) * &
                   (wRC(ixC^S,iw) - wLC(ixC^S,iw)))
           else
              where(patchf(ixC^S)==-2)
                 fLC(ixC^S,iw)=fLC(ixC^S,iw)
              elsewhere(abs(patchf(ixC^S))==1)
                 fLC(ixC^S,iw)=fCD(ixC^S,iw)
              elsewhere(patchf(ixC^S)==2)
                 fLC(ixC^S,iw)=fRC(ixC^S,iw)
              elsewhere(patchf(ixC^S)==3)
                 ! fallback option, reducing to HLL flux
                 fLC(ixC^S,iw)=Fhll(ixC^S,iw)
              elsewhere(patchf(ixC^S)==4)
                 ! fallback option, reducing to TVDLF flux
                 fLC(ixC^S,iw) = half*((fLC(ixC^S,iw)+fRC(ixC^S,iw)) &
                      -tvdlfeps * max(cmaxC(ixC^S), dabs(cminC(ixC^S))) * &
                      (wRC(ixC^S,iw)-wLC(ixC^S,iw)))
              endwhere
           end if

           if (slab) then
              fC(ixC^S,iw,idims)=fLC(ixC^S,iw)
           else
              select case (idims)
                 {case (^D)
                 fC(ixC^S,iw,^D)=block%surfaceC^D(ixC^S)*fLC(ixC^S,iw)\}
              end select
           end if

        end do ! Next iw
      end select
    end subroutine get_Riemann_flux

  end subroutine finite_volume
  !> Determine the upwinded wLC(ixL) and wRC(ixR) from w.
  !> the wCT is only used when PPM is exploited.
  subroutine upwindLR(ixI^L,ixL^L,ixR^L,idims,w,wCT,wLC,wRC,x,needprim,dxdim)
    use mod_physics
    use mod_global_parameters
    use mod_limiter

    integer, intent(in) :: ixI^L, ixL^L, ixR^L, idims
    logical, intent(in) :: needprim
    double precision, intent(in) :: dxdim
    double precision, dimension(ixI^S,1:nw) :: w, wCT
    double precision, dimension(ixI^S,1:nw) :: wLC, wRC
    double precision, dimension(ixI^S,1:ndim) :: x

    integer :: jxR^L, ixC^L, jxC^L, iw
    double precision :: wLtmp(ixI^S,1:nw), wRtmp(ixI^S,1:nw)
    double precision :: ldw(ixI^S), dwC(ixI^S)
    integer :: flagL(ixI^S), flagR(ixI^S)
    character(std_len) :: savetypelimiter

    ! Transform w,wL,wR to primitive variables
    if (needprim) then
       call phys_to_primitive(ixI^L,ixI^L,w,x)
    end if

    if(typelimiter/='ppm' .and. typelimiter /= 'mp5')then
       jxR^L=ixR^L+kr(idims,^D);
       ixCmax^D=jxRmax^D; ixCmin^D=ixLmin^D-kr(idims,^D);
       jxC^L=ixC^L+kr(idims,^D);

       do iw=1,nwflux
          if (loglimit(iw)) then
             w(ixCmin^D:jxCmax^D,iw)=dlog10(w(ixCmin^D:jxCmax^D,iw))
             wLC(ixL^S,iw)=dlog10(wLC(ixL^S,iw))
             wRC(ixR^S,iw)=dlog10(wRC(ixR^S,iw))
          end if

          dwC(ixC^S)=w(jxC^S,iw)-w(ixC^S,iw)

          savetypelimiter=typelimiter
          if(savetypelimiter=='koren') typelimiter='korenL'
          if(savetypelimiter=='cada')  typelimiter='cadaL'
          if(savetypelimiter=='cada3') typelimiter='cada3L'
          call dwlimiter2(dwC,ixI^L,ixC^L,iw,idims,ldw,dxdim)

          wLtmp(ixL^S,iw)=wLC(ixL^S,iw)+half*ldw(ixL^S)
          if(savetypelimiter=='koren')then
             typelimiter='korenR'
             call dwlimiter2(dwC,ixI^L,ixC^L,iw,idims,ldw,dxdim)
          endif
          if(savetypelimiter=='cada')then
             typelimiter='cadaR'
             call dwlimiter2(dwC,ixI^L,ixC^L,iw,idims,ldw,dxdim)
          endif
          if(savetypelimiter=='cada3')then
             typelimiter='cada3R'
             call dwlimiter2(dwC,ixI^L,ixC^L,iw,idims,ldw,dxdim)
          endif
          wRtmp(ixR^S,iw)=wRC(ixR^S,iw)-half*ldw(jxR^S)
          typelimiter=savetypelimiter

          if (loglimit(iw)) then
             w(ixCmin^D:jxCmax^D,iw)=10.0d0**w(ixCmin^D:jxCmax^D,iw)
             wLtmp(ixL^S,iw)=10.0d0**wLtmp(ixL^S,iw)
             wRtmp(ixR^S,iw)=10.0d0**wRtmp(ixR^S,iw)
          end if
       end do

       call phys_check_w(.true., ixI^L, ixL^L, wLtmp, flagL)
       call phys_check_w(.true., ixI^L, ixR^L, wRtmp, flagR)

       do iw=1,nwflux
          where (flagL(ixL^S) == 0 .and. flagR(ixR^S) == 0)
             wLC(ixL^S,iw)=wLtmp(ixL^S,iw)
             wRC(ixR^S,iw)=wRtmp(ixR^S,iw)
          end where

          ! Elsewhere, we still need to convert back when using loglimit
          if (loglimit(iw)) then
             where (flagL(ixL^S) /= 0 .or. flagR(ixR^S) /= 0)
                wLC(ixL^S,iw)=10.0d0**wLC(ixL^S,iw)
                wRC(ixR^S,iw)=10.0d0**wRC(ixR^S,iw)
             end where
          end if
       enddo
    else if (typelimiter .eq. 'ppm') then
       call PPMlimiter(ixI^L,ixM^LL,idims,w,wCT,wLC,wRC)
    else
       call MP5limiter(ixI^L,ixL^L,idims,w,wCT,wLC,wRC)
    endif

    ! Transform w,wL,wR back to conservative variables
    if(needprim)then
       call phys_to_conserved(ixI^L,ixI^L,w,x)
    endif
    call phys_to_conserved(ixI^L,ixL^L,wLC,x)
    call phys_to_conserved(ixI^L,ixR^L,wRC,x)

  end subroutine upwindLR

end module mod_finite_volume
