!> Module with finite difference methods for fluxes
module mod_finite_difference

  implicit none
  private

  public :: fd
  public :: centdiff

contains

  subroutine fd(qdt,dtfactor,ixI^L,ixO^L,idims^LIM,qtC,sCT,qt,snew,fC,fE,dxs,x)
    use mod_physics
    use mod_source, only: addsource2
    use mod_finite_volume, only: reconstruct_LR
    use mod_global_parameters
    use mod_usr_methods

    double precision, intent(in)                                     :: qdt, dtfactor, qtC, qt, dxs(ndim)
    integer, intent(in)                                              :: ixI^L, ixO^L, idims^LIM
    double precision, dimension(ixI^S,1:ndim), intent(in)            :: x

    type(state)                                                      :: sCT, snew
    double precision, dimension(ixI^S,1:nwflux,1:ndim), intent(out)  :: fC
    double precision, dimension(ixI^S,sdim:3)                        :: fE

    double precision, dimension(ixI^S,1:nwflux)                      :: fCT
    double precision, dimension(ixI^S,1:nw)                          :: fm, fp, fmR, fpL, wprim
    ! left and right constructed status in conservative form
    double precision, dimension(ixI^S,1:nw) :: wLC, wRC
    ! left and right constructed status in primitive form, needed for better performance
    double precision, dimension(ixI^S,1:nw) :: wLp, wRp
    double precision, dimension(ixI^S)      :: cmaxC, cminC
    double precision, dimension(ixI^S)      :: Hspeed
    double precision, dimension(ixO^S)      :: inv_volume
    double precision, dimension(1:ndim)     :: dxinv
    logical :: transport, active
    integer :: idims, iw, ixC^L, ix^L, hxO^L, kxC^L, kxR^L
    type(ct_velocity) :: vcts

    associate(wCT=>sCT%w,wnew=>snew%w)

    fC=0.d0
    wprim=wCT
    call phys_to_primitive(ixI^L,ixI^L,wprim,x)

    b0i = 0 ! fd uses centered values in phys_get_flux  
    do idims= idims^LIM

       !b0i=idims

       ! Get fluxes for the whole grid (mesh+nghostcells)
       {^D& ixmin^D = ixOmin^D - nghostcells * kr(idims,^D)\}
       {^D& ixmax^D = ixOmax^D + nghostcells * kr(idims,^D)\}

       hxO^L=ixO^L-kr(idims,^D);

       if(stagger_grid) then
         ! ct needs all transverse cells
         ixCmax^D=ixOmax^D+nghostcells-nghostcells*kr(idims,^D); ixCmin^D=hxOmin^D-nghostcells+nghostcells*kr(idims,^D);
         ixmax^D=ixmax^D+nghostcells-nghostcells*kr(idims,^D); ixmin^D=ixmin^D-nghostcells+nghostcells*kr(idims,^D);
         kxCmin^D=ixImin^D; kxCmax^D=ixImax^D-kr(idims,^D);
         kxR^L=kxC^L+kr(idims,^D);
         ! wRp and wLp are defined at the same locations, and will correspond to
         ! the left and right reconstructed values at a cell face. Their indexing
         ! is similar to cell-centered values, but in direction idims they are
         ! shifted half a cell towards the 'lower' direction.
         wRp(kxC^S,1:nw)=wprim(kxR^S,1:nw)
         wLp(kxC^S,1:nw)=wprim(kxC^S,1:nw)
       else
         ! ixC is centered index in the idims direction from ixOmin-1/2 to ixOmax+1/2
         ixCmax^D=ixOmax^D; ixCmin^D=hxOmin^D;
       end if

       call phys_get_flux(wCT,wprim,x,ixG^LL,ix^L,idims,fCT)

       do iw=iwstart,nwflux
          ! Lax-Friedrich splitting:
          fp(ix^S,iw) = half * (fCT(ix^S,iw) + tvdlfeps * cmax_global * wCT(ix^S,iw))
          fm(ix^S,iw) = half * (fCT(ix^S,iw) - tvdlfeps * cmax_global * wCT(ix^S,iw))
       end do ! iw loop

       ! now do the reconstruction of fp and fm:
       call reconstructL(ixI^L,ixC^L,idims,fp,fpL)
       call reconstructR(ixI^L,ixC^L,idims,fm,fmR)

       fC(ixC^S,iwstart:nwflux,idims) = fpL(ixC^S,iwstart:nwflux) + fmR(ixC^S,iwstart:nwflux)

       if(stagger_grid) then
         ! apply limited reconstruction for left and right status at cell interfaces
         call reconstruct_LR(ixI^L,ixC^L,ixC^L,idims,wprim,wLC,wRC,wLp,wRp,x,dxs(idims))
         if(H_correction) then
           call phys_get_H_speed(wprim,x,ixI^L,ixO^L,idims,Hspeed)
         end if
         call phys_get_cbounds(wLC,wRC,wLp,wRp,x,ixI^L,ixC^L,idims,Hspeed,cmaxC,cminC)
         call phys_get_ct_velocity(vcts,wLp,wRp,ixI^L,ixC^L,idims,cmaxC,cminC)
       end if

    end do !idims loop
    !b0i=0

    if(stagger_grid) call phys_update_faces(ixI^L,ixO^L,qt,qdt,wprim,fC,fE,sCT,snew,vcts)

    if(slab_uniform) then
      dxinv=-qdt/dxs
      do idims= idims^LIM
        hxO^L=ixO^L-kr(idims,^D);
        if(local_timestep) then
          do iw=iwstart,nwflux
            fC(ixI^S,iw,idims) = -block%dt(ixI^S)*dtfactor/dxs(idims) * fC(ixI^S,iw,idims)
            wnew(ixO^S,iw)=wnew(ixO^S,iw)+(fC(ixO^S,iw,idims)-fC(hxO^S,iw,idims))
          end do ! iw loop
        else
          do iw=iwstart,nwflux
            fC(ixI^S,iw,idims) = dxinv(idims) * fC(ixI^S,iw,idims)
            wnew(ixO^S,iw)=wnew(ixO^S,iw)+(fC(ixO^S,iw,idims)-fC(hxO^S,iw,idims))
          end do ! iw loop
        endif
      end do ! Next idims
    else
      inv_volume=1.d0/block%dvolume(ixO^S)
      do idims= idims^LIM
        hxO^L=ixO^L-kr(idims,^D);
        if(local_timestep) then
          do iw=iwstart,nwflux
            fC(ixI^S,iw,idims)=-block%dt(ixI^S)*dtfactor*fC(ixI^S,iw,idims)*block%surfaceC(ixI^S,idims)
            wnew(ixO^S,iw)=wnew(ixO^S,iw)+ &
                 (fC(ixO^S,iw,idims)-fC(hxO^S,iw,idims))*inv_volume(ixO^S)
          end do ! iw loop
        else
          do iw=iwstart,nwflux
            fC(ixI^S,iw,idims)=-qdt*fC(ixI^S,iw,idims)*block%surfaceC(ixI^S,idims)
            wnew(ixO^S,iw)=wnew(ixO^S,iw)+ &
                 (fC(ixO^S,iw,idims)-fC(hxO^S,iw,idims))*inv_volume(ixO^S)
          end do ! iw loop
        endif ! local_timestep
      end do ! Next idims
    end if

    if (.not.slab.and.idimsmin==1) call phys_add_source_geom(qdt,dtfactor,ixI^L,ixO^L,wCT,wprim,wnew,x)

    if(stagger_grid) call phys_face_to_center(ixO^L,snew)

    ! check and optionally correct unphysical values
    if(fix_small_values) then
       call phys_handle_small_values(.false.,wnew,x,ixI^L,ixO^L,'fd')
       if(crash) then
         ! replace erroneous values with values at previous step
         wnew=wCT
         if(stagger_grid) snew%ws=sCT%ws
       end if
    endif

    call addsource2(qdt*dble(idimsmax-idimsmin+1)/dble(ndim), &
         dtfactor*dble(idimsmax-idimsmin+1)/dble(ndim), &
         ixI^L,ixO^L,1,nw,qtC,wCT,wprim,qt,wnew,x,.false.,active)

    end associate

  end subroutine fd

  subroutine reconstructL(ixI^L,iL^L,idims,w,wLC)
    use mod_global_parameters
    use mod_mp5
    use mod_limiter
    use mod_comm_lib, only: mpistop

    integer, intent(in)             :: ixI^L, iL^L, idims
    double precision, intent(in)    :: w(ixI^S,1:nw)

    double precision, intent(out)   :: wLC(ixI^S,1:nw) 

    double precision                :: ldw(ixI^S), dwC(ixI^S)
    integer                         :: jxR^L, ixC^L, jxC^L, kxC^L, iw
    double precision                :: a2max

    select case (type_limiter(block%level))
    case (limiter_mp5)
       call MP5limiterL(ixI^L,iL^L,idims,w,wLC)
    case (limiter_weno5)
       call WENO5limiterL(ixI^L,iL^L,idims,w,wLC,1)
    case (limiter_weno5nm)
       call WENO5NMlimiterL(ixI^L,iL^L,idims,w,wLC,1)
    case (limiter_wenoz5)
       call WENO5limiterL(ixI^L,iL^L,idims,w,wLC,2)
    case (limiter_wenoz5nm)
       call WENO5NMlimiterL(ixI^L,iL^L,idims,w,wLC,2)
    case (limiter_wenozp5)
       call WENO5limiterL(ixI^L,iL^L,idims,w,wLC,3)
    case (limiter_wenozp5nm)
       call WENO5NMlimiterL(ixI^L,iL^L,idims,w,wLC,3)
    case default 

       kxCmin^D=ixImin^D; kxCmax^D=ixImax^D-kr(idims,^D);

       wLC(kxC^S,iwstart:nwflux) = w(kxC^S,iwstart:nwflux)

       jxR^L=iL^L+kr(idims,^D);

       ixCmax^D=jxRmax^D; ixCmin^D=iLmin^D-kr(idims,^D);
       jxC^L=ixC^L+kr(idims,^D);

       do iw=iwstart,nwflux
          dwC(ixC^S)=w(jxC^S,iw)-w(ixC^S,iw)

           if(need_global_a2max) then 
             a2max=a2max_global(idims)
           else
             select case(idims)
             case(1)
               a2max=schmid_rad1
             {^IFTWOD
             case(2)
               a2max=schmid_rad2}
             {^IFTHREED
             case(2)
               a2max=schmid_rad2
             case(3)
               a2max=schmid_rad3}
             case default
               call mpistop("idims is wrong in mod_limiter")
             end select
           end if

          call dwlimiter2(dwC,ixI^L,ixC^L,idims,type_limiter(block%level),ldw=ldw,a2max=a2max)

          wLC(iL^S,iw)=wLC(iL^S,iw)+half*ldw(iL^S)
       end do
    end select

  end subroutine reconstructL

  subroutine reconstructR(ixI^L,iL^L,idims,w,wRC)
    use mod_global_parameters
    use mod_mp5
    use mod_limiter
    use mod_comm_lib, only: mpistop

    integer, intent(in)             :: ixI^L, iL^L, idims
    double precision, intent(in)    :: w(ixI^S,1:nw)

    double precision, intent(out)   :: wRC(ixI^S,1:nw) 

    double precision                :: rdw(ixI^S), dwC(ixI^S)
    integer                         :: jxR^L, ixC^L, jxC^L, kxC^L, kxR^L, iw
    double precision                :: a2max

    select case (type_limiter(block%level))
    case (limiter_mp5)
       call MP5limiterR(ixI^L,iL^L,idims,w,wRC)
    case (limiter_weno5)
       call WENO5limiterR(ixI^L,iL^L,idims,w,wRC,1)
    case (limiter_weno5nm)
       call WENO5NMlimiterR(ixI^L,iL^L,idims,w,wRC,1)
    case (limiter_wenoz5)
       call WENO5limiterR(ixI^L,iL^L,idims,w,wRC,2)
    case (limiter_wenoz5nm)
       call WENO5NMlimiterR(ixI^L,iL^L,idims,w,wRC,2)
    case (limiter_wenozp5)
       call WENO5limiterR(ixI^L,iL^L,idims,w,wRC,3)
    case (limiter_wenozp5nm)
       call WENO5NMlimiterR(ixI^L,iL^L,idims,w,wRC,3)
    case default 

       kxCmin^D=ixImin^D; kxCmax^D=ixImax^D-kr(idims,^D);
       kxR^L=kxC^L+kr(idims,^D);

       wRC(kxC^S,iwstart:nwflux)=w(kxR^S,iwstart:nwflux)

       jxR^L=iL^L+kr(idims,^D);
       ixCmax^D=jxRmax^D; ixCmin^D=iLmin^D-kr(idims,^D);
       jxC^L=ixC^L+kr(idims,^D);

       do iw=iwstart,nwflux
          dwC(ixC^S)=w(jxC^S,iw)-w(ixC^S,iw)

           if(need_global_a2max) then
             a2max=a2max_global(idims)
           else
             select case(idims)
             case(1)
               a2max=schmid_rad1
             {^IFTWOD
             case(2)
               a2max=schmid_rad2}
             {^IFTHREED
             case(2)
               a2max=schmid_rad2
             case(3)
               a2max=schmid_rad3}
             case default
               call mpistop("idims is wrong in mod_limiter")
             end select
           end if

          call dwlimiter2(dwC,ixI^L,ixC^L,idims,type_limiter(block%level),rdw=rdw,a2max=a2max)

          wRC(iL^S,iw)=wRC(iL^S,iw)-half*rdw(jxR^S)
       end do
    end select

  end subroutine reconstructR

  subroutine centdiff(method,qdt,dtfactor,ixI^L,ixO^L,idims^LIM,qtC,sCT,qt,s,fC,fE,dxs,x)

    ! Advance the flow variables from global_time to global_time+qdt within ixO^L by
    ! fourth order centered differencing in space 
    ! for the dw/dt+dF_i(w)/dx_i=S type equation.
    ! wCT contains the time centered variables at time qtC for flux and source.
    ! w is the old value at qt on input and the new value at qt+qdt on output.
    use mod_physics
    use mod_finite_volume, only: reconstruct_LR
    use mod_global_parameters
    use mod_source, only: addsource2
    use mod_usr_methods
    use mod_variables
    use mod_comm_lib, only: mpistop

    integer, intent(in) :: method
    integer, intent(in) :: ixI^L, ixO^L, idims^LIM
    double precision, intent(in) :: qdt, dtfactor, qtC, qt, dxs(ndim)
    type(state)      :: sCT, s
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision :: fC(ixI^S,1:nwflux,1:ndim)
    double precision :: fE(ixI^S,sdim:3)

    double precision :: v(ixI^S,ndim), f(ixI^S, nwflux)

    double precision, dimension(ixI^S,1:nw) :: wprim, wLC, wRC
    ! left and right constructed status in primitive form, needed for better performance
    double precision, dimension(ixI^S,1:nw) :: wLp, wRp
    double precision, dimension(ixI^S)      :: vLC, cmaxLC, cmaxRC
    double precision, dimension(ixI^S,1:number_species)      :: cmaxC
    double precision, dimension(ixI^S,1:number_species)      :: cminC
    double precision, dimension(ixI^S)      :: Hspeed
    double precision, dimension(ixO^S)      :: inv_volume

    double precision :: dxinv(1:ndim)
    integer :: idims, iw, ix^L, hxO^L, ixC^L, jxC^L, hxC^L, kxC^L, kkxC^L, kkxR^L
    type(ct_velocity) :: vcts
    logical :: transport, new_cmax, patchw(ixI^S), active

    associate(wCT=>sCT%w,w=>s%w)
    ! two extra layers are needed in each direction for which fluxes are added.
    ix^L=ixO^L;
    do idims= idims^LIM
       ix^L=ix^L^LADD2*kr(idims,^D);
    end do

    if (ixI^L^LTix^L|.or.|.or.) then
       call mpistop("Error in centdiff: Non-conforming input limits")
    end if

    fC=0.d0
    wprim=wCT
    call phys_to_primitive(ixI^L,ixI^L,wprim,x)

    ! get fluxes
    do idims= idims^LIM
       b0i=idims

       ix^L=ixO^L^LADD2*kr(idims,^D); 
       hxO^L=ixO^L-kr(idims,^D);

       if(stagger_grid) then
         ! ct needs all transverse cells
         ixCmax^D=ixOmax^D+nghostcells-nghostcells*kr(idims,^D);
         ixCmin^D=hxOmin^D-nghostcells+nghostcells*kr(idims,^D);
         ixmax^D=ixmax^D+nghostcells-nghostcells*kr(idims,^D);
         ixmin^D=ixmin^D-nghostcells+nghostcells*kr(idims,^D);
       else
         ! ixC is centered index in the idims direction from ixOmin-1/2 to ixOmax+1/2
         ixCmax^D=ixOmax^D; ixCmin^D=hxOmin^D;
       end if
       hxC^L=ixC^L-kr(idims,^D); 
       jxC^L=ixC^L+kr(idims,^D); 
       kxC^L=ixC^L+2*kr(idims,^D); 

       kkxCmin^D=ixImin^D; kkxCmax^D=ixImax^D-kr(idims,^D);
       kkxR^L=kkxC^L+kr(idims,^D);
       wRp(kkxC^S,1:nwflux)=wprim(kkxR^S,1:nwflux)
       wLp(kkxC^S,1:nwflux)=wprim(kkxC^S,1:nwflux)

       ! apply limited reconstruction for left and right status at cell interfaces
       call reconstruct_LR(ixI^L,ixC^L,ixC^L,idims,wprim,wLC,wRC,wLp,wRp,x,dxs(idims))

       if(stagger_grid) then
         if(H_correction) then
           call phys_get_H_speed(wprim,x,ixI^L,ixO^L,idims,Hspeed)
         end if
         call phys_get_cbounds(wLC,wRC,wLp,wRp,x,ixI^L,ixC^L,idims,Hspeed,cmaxC,cminC)
         call phys_get_ct_velocity(vcts,wLp,wRp,ixI^L,ixC^L,idims,cmaxC,cminC)
       end if

       ! Calculate velocities from upwinded values
       call phys_get_cmax(wLp,x,ixG^LL,ixC^L,idims,cmaxLC)
       call phys_get_cmax(wRp,x,ixG^LL,ixC^L,idims,cmaxRC)
       ! now take the maximum of left and right states
       vLC(ixC^S)=max(cmaxRC(ixC^S),cmaxLC(ixC^S))

       call phys_get_flux(wCT,wprim,x,ixI^L,ix^L,idims,f)

       ! Center flux to interface
       if(method==fs_cd) then
          fC(ixC^S,iwstart:nwflux,idims)=half*(f(ixC^S,iwstart:nwflux)+f(jxC^S,iwstart:nwflux))
       else
          ! f_i+1/2= (-f_(i+2) +7 f_(i+1) + 7 f_i - f_(i-1))/12
          fC(ixC^S,iwstart:nwflux,idims)=(-f(kxC^S,iwstart:nwflux)+7.0d0*(f(jxC^S,iwstart:nwflux) + &
               f(ixC^S,iwstart:nwflux))-f(hxC^S,iwstart:nwflux))/12.0d0

          do iw=iwstart,nwflux
             ! add rempel dissipative flux, only second order version for now
             fC(ixC^S,iw,idims)=fC(ixC^S,iw,idims)-tvdlfeps*half*vLC(ixC^S) &
                  *(wRC(ixC^S,iw)-wLC(ixC^S,iw))
          end do
       end if

    end do       !next idims
    b0i=0

    if(stagger_grid) call phys_update_faces(ixI^L,ixO^L,qt,qdt,wprim,fC,fE,sCT,s,vcts)

    if(slab_uniform) then
      dxinv=-qdt/dxs
      do idims= idims^LIM
        hxO^L=ixO^L-kr(idims,^D);
        if(local_timestep) then
          do iw=iwstart,nwflux
            fC(ixI^S,iw,idims)=-block%dt(ixI^S)*dtfactor/dxs(idims)*fC(ixI^S,iw,idims)
            ! result: f_(i+1/2)-f_(i-1/2) = [-f_(i+2)+8(f_(i+1)-f_(i-1))+f_(i-2)]/12
            w(ixO^S,iw)=w(ixO^S,iw)+(fC(ixO^S,iw,idims)-fC(hxO^S,iw,idims))
          end do    !next iw
        else
          do iw=iwstart,nwflux
            fC(ixI^S,iw,idims)=dxinv(idims)*fC(ixI^S,iw,idims)
            ! result: f_(i+1/2)-f_(i-1/2) = [-f_(i+2)+8(f_(i+1)-f_(i-1))+f_(i-2)]/12
            w(ixO^S,iw)=w(ixO^S,iw)+(fC(ixO^S,iw,idims)-fC(hxO^S,iw,idims))
          end do    !next iw
        endif
      end do ! Next idims
    else
      inv_volume=1.d0/block%dvolume
      do idims= idims^LIM
        hxO^L=ixO^L-kr(idims,^D);
        if(local_timestep) then
          do iw=iwstart,nwflux
            fC(ixI^S,iw,idims)=-block%dt(ixI^S)*dtfactor*block%surfaceC(ixI^S,idims)*fC(ixI^S,iw,idims)
            w(ixO^S,iw)=w(ixO^S,iw)+ &
                 (fC(ixO^S,iw,idims)-fC(hxO^S,iw,idims))*inv_volume(ixO^S)
          end do    !next iw
        else
          do iw=iwstart,nwflux
            fC(ixI^S,iw,idims)=-qdt*block%surfaceC(ixI^S,idims)*fC(ixI^S,iw,idims)
            w(ixO^S,iw)=w(ixO^S,iw)+ &
                 (fC(ixO^S,iw,idims)-fC(hxO^S,iw,idims))*inv_volume(ixO^S)
          end do    !next iw
        endif 
      end do ! Next idims
    end if

    if (.not.slab.and.idimsmin==1) call phys_add_source_geom(qdt,dtfactor,ixI^L,ixO^L,wCT,wprim,w,x)

    if(stagger_grid) call phys_face_to_center(ixO^L,s)

    ! check and optionally correct unphysical values
    if(fix_small_values) then
       call phys_handle_small_values(.false.,w,x,ixI^L,ixO^L,'centdiff')
    endif

    call addsource2(qdt*dble(idimsmax-idimsmin+1)/dble(ndim), &
         dtfactor*dble(idimsmax-idimsmin+1)/dble(ndim), &
         ixI^L,ixO^L,1,nw,qtC,wCT,wprim,qt,w,x,.false.,active)

    end associate
  end subroutine centdiff

end module mod_finite_difference
