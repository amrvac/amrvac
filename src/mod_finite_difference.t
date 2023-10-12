!> Module with finite difference methods for fluxes
module mod_finite_difference

  implicit none
  private

  public :: fd
  public :: finite_difference

contains

  subroutine fd(qdt,ixI^L,ixO^L,idims^LIM,qtC,sCT,qt,snew,sold,fC,fE,dx^D)
    use mod_physics
    use mod_geometry
    use mod_source, only: addsource2
    use mod_finite_volume, only: reconstruct_LR, metric_interpolation, get_theta
    use mod_limiter
    use mod_global_parameters
    use mod_usr_methods

    double precision, intent(in)                        :: qdt, qtC, qt, dx^D
    integer, intent(in)                                 :: ixI^L, ixO^L, idims^LIM
    type(state)                                         :: sCT, snew, sold
    double precision, dimension(ixI^S,1:nwflux,1:ndim)  :: fC
    double precision, dimension(ixI^S,7-2*ndim:3)       :: fE
    ! cell-face location coordinates
    type(mesh_t), target                                :: mesh_i
    ! cell-face metric
    type(metric_t), target                              :: metric_i
    ! primitive state at cell center
    type(state)                                         :: sCTp
    ! left and right constructed status in conservative form
    type(state)                                         :: sL, sR
    ! left and right constructed status in primitive form, needed for better performance
    type(state)                                         :: sLp, sRp

    integer                                             :: idims, iw, ixC^L, ix^L, hxO^L, kxC^L, kxR^L
    double precision, dimension(ixI^S,1:nwflux)         :: fCT
    double precision, dimension(ixI^S,1:nwflux)         :: fm, fp, fmR, fpL
    ! cell-face location coordinates
    double precision, dimension(ixI^S,1:ndim)           :: xi
    double precision, dimension(ixI^S,1:number_species) :: cmaxC, cminC
    double precision, dimension(1:ndim)                 :: dxinv, dxdim
    ! when CT is being used
    type(ct_velocity)                                   :: vcts
    ! variables that related to positivity preserving limiter
    ! low order flux
    double precision, dimension(:^D&,:), allocatable  :: fC_low

    sCTp = sCT
    call phys_to_primitive(ixI^L,ixI^L,sCTp)

    sL  = sCT  ; sR  = sCT
    sLp = sCTp ; sRp = sCTp

    allocate(mesh_i%x(ixI^S, 1:ndim))
    allocate(metric_i%vars(ixI^S, 1:nmetric))
    sL%mesh   => mesh_i    ; sR%mesh   => mesh_i
    sL%metric => metric_i  ; sR%metric => metric_i
    sLp%mesh   => mesh_i   ; sRp%mesh   => mesh_i
    sLp%metric => metric_i ; sRp%metric => metric_i

    associate( w_new => snew%w,     w_old => sold%w,   &
               x     => sCT%mesh%x, xi    => mesh_i%x, &
               wCT   => sCT%w,      wprim => sCTp%w,   &
               wLC   => sL%w,       wRC   => sR%w,     &
               wLp   => sLp%w,      wRp   => sRp%w     &
             )

    if (positivity_preserving) then
       allocate(fC_low(ixI^S,1:nwflux))
       fC_low=0.0d0
    end if

    fC=0.d0
    !if (stagger_grid) then
    !   primL=0.0d0
    !   primR=0.0d0
    !end if

    ^D&dxdim(^D)=dx^D;
    do idims= idims^LIM

       ! get cell-face coordinates from cell-center
       mesh_i%x(ixI^S,1:ndim) = sCT%mesh%x(ixI^S,1:ndim)
       mesh_i%x(ixI^S,idims)  = mesh_i%x(ixI^S,idims) + 0.5d0 * sCT%mesh%dx(ixI^S,idims)

       ! Get fluxes for the whole grid (mesh+nghostcells)
       {^D& ixmin^D = ixOmin^D - nghostcells * kr(idims,^D)\}
       {^D& ixmax^D = ixOmax^D + nghostcells * kr(idims,^D)\}

       hxO^L=ixO^L-kr(idims,^D);

       if(stagger_grid) then
         ! ct needs all transverse cells
         ixCmax^D=ixOmax^D+nghostcells-nghostcells*kr(idims,^D); 
         ixCmin^D=hxOmin^D-nghostcells+nghostcells*kr(idims,^D);
         ixmax^D=ixmax^D+nghostcells-nghostcells*kr(idims,^D); 
         ixmin^D=ixmin^D-nghostcells+nghostcells*kr(idims,^D);
         kxCmin^D=ixImin^D; kxCmax^D=ixImax^D-kr(idims,^D);
         kxR^L=kxC^L+kr(idims,^D);
         ! wRp and wLp are defined at the same locations, and will correspond to
         ! the left and right reconstructed values at a cell face. Their indexing
         ! is similar to cell-centered values, but in direction idims they are
         ! shifted half a cell towards the 'lower' direction.
         wRp(kxC^S,1:nw) = wprim(kxR^S,1:nw)
         wLp(kxC^S,1:nw) = wprim(kxC^S,1:nw)
       else
         ! ixC is centered index in the idims direction from ixOmin-1/2 to ixOmax+1/2
         ixCmax^D=ixOmax^D; ixCmin^D=hxOmin^D;
       end if

       call phys_get_flux(sCT,sCTp,ixG^LL,ix^L,idims,fCT)

       do iw=1,nwflux
          ! Lax-Friedrich splitting:
          fp(ix^S,iw) = half * (fCT(ix^S,iw) + tvdlfeps * cmax_global * wCT(ix^S,iw))
          fm(ix^S,iw) = half * (fCT(ix^S,iw) - tvdlfeps * cmax_global * wCT(ix^S,iw))
       end do ! iw loop

       ! now do the reconstruction of fp and fm:
       call reconstructL(type_limiter(sCT%mesh%level),ixI^L,ixC^L,idims,fp,fpL)
       call reconstructR(type_limiter(sCT%mesh%level),ixI^L,ixC^L,idims,fm,fmR)

       fC(ixC^S,1:nwflux,idims) = fpL(ixC^S,1:nwflux) + fmR(ixC^S,1:nwflux)

       if (positivity_preserving) then
          call reconstructL(limiter_woodward,ixI^L,ixC^L,idims,fp,fpL)
          call reconstructR(limiter_woodward,ixI^L,ixC^L,idims,fm,fmR)
          fC_low(ixC^S,1:nwflux) = fpL(ixC^S,1:nwflux) + fmR(ixC^S,1:nwflux)
          call positivity_preserving_limiter()
       end if

       !if(associated(usr_set_flux)) call usr_set_flux(ixI^L,ixC^L,qt,wLC,wRC,wLp,wRp,sCT,idims,fC)

       if(stagger_grid) then
         ! apply limited reconstruction for left and right status at cell interfaces
         ! get cell-face metric
         ! fixme: maybe this part can be faster
         call metric_interpolation(ixI^L,ixC^L,idims,nmetric,sCT%metric%vars,x,metric_i%vars,xi)
         call reconstruct_LR(type_limiter(sCT%mesh%level), &
           ixI^L,ixC^L,ixC^L,idims,sCTp,sL,sR,sLp,sRp,dxdim(idims))
         ! special modification of left and right status before flux evaluation
         call phys_modify_wLR(ixI^L,ixC^L,qt,idims,sCT,sL,sR,sLp,sRp)
         call phys_get_cbounds_and_vct(sL,sR,sLp,sRp,ixI^L,ixC^L,idims,vcts,cmaxC,cminC)
       end if

    end do !idims loop

    if (stagger_grid) call phys_update_faces(ixI^L,ixO^L,qt,qdt,fC,fE,sCTp,snew,vcts)

    if(slab_uniform) then
      ^D&dxinv(^D)=-qdt/dx^D;
      do idims= idims^LIM
        hxO^L=ixO^L-kr(idims,^D);
        do iw=1,nwflux
          fC(ixI^S,iw,idims) = dxinv(idims) * fC(ixI^S,iw,idims)
          w_new(ixO^S,iw)=w_new(ixO^S,iw)+(fC(ixO^S,iw,idims)-fC(hxO^S,iw,idims))
        end do ! iw loop
      end do ! Next idims
    else
      do idims= idims^LIM
        hxO^L=ixO^L-kr(idims,^D);
        do iw=1,nwflux
          fC(ixI^S,iw,idims)=-qdt*fC(ixI^S,iw,idims)*block%mesh%surfaceC(ixI^S,idims)
          w_new(ixO^S,iw)=w_new(ixO^S,iw)+ &
               (fC(ixO^S,iw,idims)-fC(hxO^S,iw,idims))/block%mesh%dvolume(ixO^S)
        end do ! iw loop
      end do ! Next idims
    end if

    if (.not.slab.and.idimsmin==1) &
         call phys_add_source_geom(qdt,sCT,sCTp,ixI^L,ixO^L,w_new)

    if(stagger_grid) call phys_face_to_center(ixO^L,snew)

    ! check and optionally correct unphysical values
    !if(fix_small_values) then
    !   call phys_handle_small_values(.false.,cons_new,x,ixI^L,ixO^L,'fd')
    !endif

    call addsource2(qdt*dble(idimsmax-idimsmin+1)/dble(ndim), &
         ixI^L,ixO^L,1,nwflux,qtC,sCTp,sold,qt,w_new,.false.)

    if (positivity_preserving) then
       deallocate(fC_low)
    end if

    deallocate(mesh_i%x)
    deallocate(metric_i%vars)

    end associate

    contains

    subroutine positivity_preserving_limiter()
      use mod_hd_parameters, only: small_rho, small_e
      use mod_geometry, only: get_radii_pt
      implicit none
      double precision, dimension(ixI^S)   :: eps_D, eps_tau
      double precision, dimension(ixI^S)   :: inv_lambda
      double precision, dimension(ixI^S)   :: theta, theta_tau
      double precision, dimension(ixI^S)   :: flow_dA, fhigh_dA
      double precision, parameter          :: eps_fac = 1.0d-6
      double precision                     :: radii_pt, x_tmp(1:^ND)
      integer                              :: ix^D

      {do ix^D = ixI^LIM^D \}
         x_tmp         = mesh_i%x(ix^D,:)
         radii_pt      = get_radii_pt(x_tmp)
         eps_D(ix^D)   = eps_fac * small_rho
         eps_tau(ix^D) = eps_fac * small_e 
      {end do^D&\}

      ! get theta for D
      if (slab_uniform) then
         inv_lambda(ixI^S) = dxdim(idims)/(qdt)
         flow_dA(ixI^S)  = fC_low(ixI^S,D_)
         fhigh_dA(ixI^S) = fC(ixI^S,D_,idims)
      else
         inv_lambda(ixI^S) = sCT%mesh%dvolume(ixI^S)/(qdt)
         flow_dA(ixI^S)  = fC_low(ixI^S,D_)   * sCT%mesh%surfaceC(ixI^S,idims)
         fhigh_dA(ixI^S) = fC(ixI^S,D_,idims) * sCT%mesh%surfaceC(ixI^S,idims)
      end if
      call get_theta(ixI^L,ixC^L,idims,eps_D,inv_lambda,sCT%w(ixI^S,D_),flow_dA,fhigh_dA,theta)

      ! get theta for tau
      if (slab_uniform) then
         flow_dA(ixI^S)  = fC_low(ixI^S,tau_)
         fhigh_dA(ixI^S) = fC(ixI^S,tau_,idims)
      else
         flow_dA(ixI^S)  = fC_low(ixI^S,tau_)   * sCT%mesh%surfaceC(ixI^S,idims)
         fhigh_dA(ixI^S) = fC(ixI^S,tau_,idims) * sCT%mesh%surfaceC(ixI^S,idims)
      end if
      ! Note that cons(tau_) may be -ve for tabulate EOS
      call get_theta(ixI^L,ixC^L,idims,eps_tau,inv_lambda,sCT%w(ixI^S,tau_),flow_dA,fhigh_dA,theta_tau)

      theta(ixC^S) = min(theta(ixC^S),theta_tau(ixC^S))

      do iw = 1,nwflux
         ! note that pp limiter cannot act on stagger grid
         if ( stagger_grid ) then
            if ( (iw>=Bcons(1)) .and. (iw<=Bcons(ndim)) ) cycle
         end if
         fC(ixC^S,iw,idims) = theta(ixC^S)*fC(ixC^S,iw,idims) &
              + (1.0d0-theta(ixC^S))*fC_low(ixC^S,iw)
      end do ! Next iw
    end subroutine positivity_preserving_limiter

  end subroutine fd

  subroutine finite_difference(method,qdt,ixI^L,ixO^L,idims^LIM,qtC,sCT,qt,snew,sold,fC,fE,dx^D)
    use mod_physics
    use mod_geometry
    use mod_source, only: addsource2
    use mod_finite_volume, only: reconstruct_LR, metric_interpolation, get_theta
    use mod_limiter
    use mod_global_parameters
    use mod_usr_methods

    integer, intent(in)                                 :: method
    double precision, intent(in)                        :: qdt, qtC, qt, dx^D
    integer, intent(in)                                 :: ixI^L, ixO^L, idims^LIM
    type(state)                                         :: sCT, snew, sold
    double precision, dimension(ixI^S,1:nwflux,1:ndim)  :: fC
    double precision, dimension(ixI^S,7-2*ndim:3)       :: fE

    ! cell-face location coordinates
    type(mesh_t), target                                :: mesh_i
    ! cell-face metric
    type(metric_t), target                              :: metric_i
    ! primitive state at cell center
    type(state)                                         :: sCTp
    ! left and right constructed status in conservative form
    type(state)                                         :: sL, sR
    ! left and right constructed status in primitive form, needed for better performance
    type(state)                                         :: sLp, sRp

    integer, dimension(ixI^S)                           :: patchf
    integer                                             :: idims, iw, ii
    integer                                             :: ix^L, hxO^L, ixC^L, ixCR^L, kxC^L, kxR^L
    double precision, dimension(ixI^S)                  :: inv_volume
    double precision, dimension(1:ndim)                 :: dxinv, dxdim
    double precision, dimension(ixI^S,1:nwflux)         :: fCT, fLC, fRC
    double precision, dimension(ixI^S,1:number_species) :: tvdlf_eps
    double precision, dimension(ixI^S,1:number_species) :: cmaxC, cminC
    ! when CT is being used
    type(ct_velocity)                                   :: vcts
    ! variables that related to positivity preserving limiter
    ! low order flux
    double precision, dimension(:^D&,:), allocatable  :: fC_low

    sCTp = sCT
    call phys_to_primitive(ixI^L,ixI^L,sCTp)

    sL  = sCT  ; sR  = sCT
    sLp = sCTp ; sRp = sCTp

    allocate(mesh_i%x(ixI^S, 1:ndim))
    allocate(metric_i%vars(ixI^S, 1:nmetric))
    sL%mesh   => mesh_i    ; sR%mesh   => mesh_i
    sL%metric => metric_i  ; sR%metric => metric_i
    sLp%mesh   => mesh_i   ; sRp%mesh   => mesh_i
    sLp%metric => metric_i ; sRp%metric => metric_i

    associate( w_new => snew%w,     w_old => sold%w,   &
               x     => sCT%mesh%x, xi    => mesh_i%x, &
               wCT   => sCT%w,      wprim => sCTp%w,   &
               wLC   => sL%w,       wRC   => sR%w,     &
               wLp   => sLp%w,      wRp   => sRp%w     &
             )

    if (positivity_preserving) then
       allocate(fC_low(ixI^S,1:nwflux))
       fC_low=0.0d0
    end if

    fC  = 0.d0
    fLC = 0.d0; fRC = 0.d0
    tvdlf_eps = 1.0d0
    !sLp%w=0.0d0; sRp%w=0.0d0

    ! The flux calculation contracts by one in the idims direction it is applied.
    ! The limiter contracts the same directions by one more, so expand ixO by 2.
    ix^L=ixO^L;
    do idims= idims^LIM
       ix^L=ix^L^LADD2*kr(idims,^D);
    end do
    if (ixI^L^LTix^L|.or.|.or.) &
         call mpistop("Error in fv : Nonconforming input limits")

    ^D&dxdim(^D)=dx^D;

    do idims= idims^LIM

       ! Assemble indices
       hxO^L=ixO^L-kr(idims,^D);
       kxCmin^D=ixImin^D; kxCmax^D=ixImax^D-kr(idims,^D);
       kxR^L=kxC^L+kr(idims,^D);

       ixCmin^D=hxOmin^D-(nghostcells-phys_extra_ghostcells)*(1-kr(idims,^D));
       ixCmax^D=ixOmax^D+(nghostcells-phys_extra_ghostcells)*(1-kr(idims,^D)); 

       ! Determine stencil size, indices to reconstruct to
       {ixCRmin^D = max(ixCmin^D - phys_wider_stencil,ixGlo^D)\}
       {ixCRmax^D = min(ixCmax^D + phys_wider_stencil,ixGhi^D)\}

       ! get cell-face coordinates from cell-center
       mesh_i%x(ixI^S,1:ndim) = sCT%mesh%x(ixI^S,1:ndim)
       mesh_i%x(ixI^S,idims)  = mesh_i%x(ixI^S,idims) + 0.5d0 * sCT%mesh%dx(ixI^S,idims)

       ! modify tvdlf_eps, the prefactor of the dissipation term
       call phys_modify_tvdlfeps(ixI^L,ixC^L,idims,sCTp,tvdlf_eps)

       ! wR and wL are defined at the same locations, and will correspond to
       ! the left and right reconstructed values at a cell face. Their indexing
       ! is similar to cell-centered values, but in direction idims they are
       ! shifted half a cell towards the 'higher' direction.
       wRp(kxC^S,1:nw) = wprim(kxR^S,1:nw)
       wLp(kxC^S,1:nw) = wprim(kxC^S,1:nw)

       ! get cell-face metric
       !metric_i%vars(kxC^S,1:nmetric) = sCT%metric%vars(kxR^S,1:nmetric)
       call metric_interpolation(ixI^L,ixCR^L,idims,nmetric,sCT%metric%vars,x,metric_i%vars,xi)

       ! low order TVD-Lax-Friedrich flux first.
       if (positivity_preserving) then
         ! fixme: make the limiter as an option here
         call reconstruct_LR(limiter_minmod, &
           ixI^L,ixCR^L,ixCR^L,idims,sCTp,sL,sR,sLp,sRp,dxdim(idims))
         ! special modification of left and right status before flux evaluation
         call phys_modify_wLR(ixI^L,ixCR^L,qt,idims,sCT,sL,sR,sLp,sRp)
         call phys_get_flux(sL,sLp,ixI^L,ixC^L,idims,fLC)
         call phys_get_flux(sR,sRp,ixI^L,ixC^L,idims,fRC)
         call phys_get_cbounds(sL,sR,sLp,sRp,ixI^L,ixC^L,idims,cmaxC)
         do ii=1,number_species
           call get_Riemann_flux_tvdlf(start_indices(ii),stop_indices(ii))
         end do
         do iw=1,nwflux
           fC_low(ixC^S,iw)=fC(ixC^S,iw,idims)
         end do
       end if ! end pp limiter

       ! apply limited reconstruction for left and right status at cell interfaces
       call reconstruct_LR(type_limiter(sCT%mesh%level), &
         ixI^L,ixCR^L,ixCR^L,idims,sCTp,sL,sR,sLp,sRp,dxdim(idims))

       ! special modification of left and right status before flux evaluation
       call phys_modify_wLR(ixI^L,ixCR^L,qt,idims,sCT,sL,sR,sLp,sRp)

       ! evaluate physical fluxes according to reconstructed status
       call phys_get_flux(sL,sLp,ixI^L,ixC^L,idims,fLC)
       call phys_get_flux(sR,sRp,ixI^L,ixC^L,idims,fRC)

       ! estimating bounds for the minimum and maximum signal velocities
       if(method==fs_fd_tvdlf) then
         if (.not.stagger_grid) then
           call phys_get_cbounds(sL,sR,sLp,sRp,ixI^L,ixC^L,idims,cmaxC)
         else
           call phys_get_cbounds_and_vct(sL,sR,sLp,sRp,ixI^L,ixC^L,idims,vcts,cmaxC)
         end if
       else
         if (.not.stagger_grid) then
           call phys_get_cbounds(sL,sR,sLp,sRp,ixI^L,ixC^L,idims,cmaxC,cminC)
         else
           call phys_get_cbounds_and_vct(sL,sR,sLp,sRp,ixI^L,ixC^L,idims,vcts,cmaxC,cminC)
         end if
       end if
       
       ! use approximate Riemann solver to get flux at interfaces
       select case(method)
       case(fs_fd_tvdlf)
         do ii=1,number_species
           call get_Riemann_flux_tvdlf(start_indices(ii),stop_indices(ii))
         end do
       case(fs_fd_hll)
         do ii=1,number_species
           call get_Riemann_flux_hll(start_indices(ii),stop_indices(ii))
         end do
       case default
         call mpistop('unkown Riemann flux in finite volume')
       end select

       call phys_get_flux(sCT,sCTp,ixI^L,ixI^L,idims,fCT)
       call high_order_flux(ixI^L,ixC^L,idims,fC(ixI^S,1:nwflux,idims),fCT)

       ! If use positivity preserving limiter, with fC and fC_low, work out the modify the flux
       if (positivity_preserving) call positivity_preserving_limiter()

    end do ! Next idims

    if (stagger_grid) call phys_update_faces(ixI^L,ixO^L,qt,qdt,fC,fE,sCTp,snew,vcts)

    do idims= idims^LIM
       do iw=1,nwflux
         if ( flux_type(idims, iw) == flux_nul ) fC(ixI^S,iw,idims)=0.0d0
       end do
    end do ! Next idims

    if(slab_uniform) then
      ^D&dxinv(^D)=-qdt/dx^D;
      do idims= idims^LIM
        hxO^L=ixO^L-kr(idims,^D);

        ! Multiply the fluxes by -dt/dx since Flux fixing expects this
        fC(ixI^S,1:nwflux,idims)=dxinv(idims)*fC(ixI^S,1:nwflux,idims)

        w_new(ixO^S,1:nwflux)=w_new(ixO^S,1:nwflux) &
                   + (fC(ixO^S,1:nwflux,idims)-fC(hxO^S,1:nwflux,idims))
      end do ! Next idims
    else
      inv_volume(ixO^S) = 1.d0/sCT%mesh%dvolume(ixO^S)
      do idims= idims^LIM
        hxO^L=ixO^L-kr(idims,^D);

        do iw=1,nwflux
          fC(ixI^S,iw,idims)=-qdt*fC(ixI^S,iw,idims)*sCT%mesh%surfaceC(ixI^S,idims)
          w_new(ixO^S,iw)=w_new(ixO^S,iw) &
                + inv_volume(ixO^S) * (fC(ixO^S,iw,idims)-fC(hxO^S,iw,idims)) 
        end do
      end do ! Next idims
    end if

    if (.not.slab.and.idimsmin==1) &
         call phys_add_source_geom(qdt,sCT,sCTp,ixI^L,ixO^L,w_new)

    if (stagger_grid) call phys_face_to_center(ixO^L,snew)

    call addsource2(qdt*dble(idimsmax-idimsmin+1)/dble(ndim), &
         ixI^L,ixO^L,1,nwflux,qtC,sCTp,sold,qt,w_new,.false.)

    if (positivity_preserving) then
       deallocate(fC_low)
    end if

    deallocate(mesh_i%x)
    deallocate(metric_i%vars)

  end associate

  contains

    subroutine get_Riemann_flux_tvdlf(iws,iwe)
      integer, intent(in) :: iws,iwe
      double precision    :: fac(ixC^S)

      fac(ixC^S) = -0.5d0*tvdlfeps*cmaxC(ixC^S,ii)
      ! Calculate fLC=f(uL_j+1/2) and fRC=f(uR_j+1/2) for each iw
      do iw=iws,iwe
         ! To save memory we use fLC to store (F_L+F_R)/2=0.5d0*(fLC+fRC)
         fLC(ixC^S, iw)=0.5d0*(fLC(ixC^S, iw)+fRC(ixC^S, iw))

         select case (flux_type(idims, iw))
         case (flux_no_dissipation)
            ! do nothing
         case (flux_asym_diffusion_e,flux_asym_diffusion_f)
            ! note: tvdlf seems cannot be used in radiation hydrodynamics
            fLC(ixC^S, iw)=fLC(ixC^S, iw) + tvdlf_eps(ixC^S, ii)*fac(ixC^S)*(sR%w(ixC^S,iw)-sL%w(ixC^S,iw))
         case default
            ! Add TVDLF dissipation to the flux
            fLC(ixC^S, iw)=fLC(ixC^S, iw) + fac(ixC^S)*(sR%w(ixC^S,iw)-sL%w(ixC^S,iw))
         end select

         fC(ixC^S,iw,idims)=fLC(ixC^S, iw)
      end do ! Next iw
    end subroutine get_Riemann_flux_tvdlf

    subroutine get_Riemann_flux_hll(iws,iwe)
      integer, intent(in) :: iws,iwe
      integer             :: ix^D

      ! Calculate fLC=f(uL_j+1/2) and fRC=f(uR_j+1/2) for each iw
      do iw=iws,iwe

         select case (flux_type(idims, iw))
         case (flux_tvdlf)
            fC(ixC^S,iw,idims) = 0.5d0*(fLC(ixC^S, iw) + fRC(ixC^S, iw) &
                 -tvdlfeps*max(cmaxC(ixC^S, ii), dabs(cminC(ixC^S, ii))) * &
                 (sR%w(ixC^S,iw)-sL%w(ixC^S,iw)))
         case default
            {do ix^DB=ixCmin^DB,ixCmax^DB\} 
              if(cminC(ix^D,ii) >= 0.0d0) then
                fC(ix^D,iw,idims)=fLC(ix^D,iw)
              else if(cmaxC(ix^D,ii) <= 0.0d0) then
                fC(ix^D,iw,idims)=fRC(ix^D,iw)
              else
                ! Add hll dissipation to the flux
                fC(ix^D,iw,idims)=(cmaxC(ix^D,ii)*fLC(ix^D,iw)-cminC(ix^D,ii)*fRC(ix^D,iw)&
                      +tvdlfeps*cminC(ix^D,ii)&
                       *cmaxC(ix^D,ii)*(sR%w(ix^D,iw)-sL%w(ix^D,iw)))&
                      /(cmaxC(ix^D,ii)-cminC(ix^D,ii))
              end if
            {end do\}
         end select

      end do ! Next iw
    end subroutine get_Riemann_flux_hll

    subroutine high_order_flux(ixI^L,ixC^L,idims,fC,fCT)
      integer, intent(in)                :: ixI^L, ixC^L, idims
      double precision, intent(in)       :: fCT(ixI^S,1:nwflux)
      double precision, intent(inout)    :: fC(ixI^S,1:nwflux) 
      double precision, dimension(ixI^S) :: df2, df4
      integer                            :: ixCpp^L, ixCp^L, ixCm^L, ixCmm^L
      ixCp^L=ixC^L+kr(idims,^D);
      ixCpp^L=ixCp^L+kr(idims,^D);
      ixCm^L=ixC^L-kr(idims,^D);

      do iw = 1,nwflux
         df2(ixC^S) = ( fCT(ixC^S,iw) - 2.0d0 * fC(ixC^S,iw) + fCT(ixCp^S,iw) ) / 6.0d0
         df4(ixC^S) = ( fCT(ixCm^S,iw) - 9.0d0 * fCT(ixC^S,iw) &
                      + 1.6d1 * fC(ixC^S,iw) &
                      + fCT(ixCpp^S,iw) - 9.0d0 * fCT(ixCp^S,iw) &
                          ) / 1.8d2
         fC(ixC^S,iw) = fC(ixC^S,iw) - df2(ixC^S) + df4(ixC^S)
      end do ! Next iw

      ! Del Zanna 
      !ixCmm^L=ixCm^L-kr(idims,^D);
      !do iw = 1,nwflux
      !   df2(ixC^S) = ( fC(ixCm^S,iw) - 2.0d0 * fC(ixC^S,iw) + fC(ixCp^S,iw) ) / 2.4d1
      !   !df4(ixC^S) = ( fC(ixCmm^S,iw) - 4.0d0 * fC(ixCm^S,iw) &
      !   !             + 6.0d0 * fC(ixC^S,iw) &
      !   !             + fC(ixCpp^S,iw) - 4.0d0 * fC(ixCp^S,iw) &
      !   !                 ) * 3.0d0 / 6.4d2
      !   !fC(ixC^S,iw) = fC(ixC^S,iw) - df2(ixC^S) !+ df4(ixC^S)
      !end do ! Next iw
    end subroutine high_order_flux

    subroutine positivity_preserving_limiter()
      use mod_hd_parameters, only: small_rho, small_e
      use mod_geometry, only: get_radii_pt
      implicit none
      double precision, dimension(ixI^S)   :: eps_D, eps_tau
      double precision, dimension(ixI^S)   :: inv_lambda
      double precision, dimension(ixI^S)   :: theta, theta_tau
      double precision, dimension(ixI^S)   :: flow_dA, fhigh_dA
      double precision, parameter          :: eps_fac = 1.0d-6
      double precision                     :: radii_pt, x_tmp(1:^ND)
      integer                              :: ix^D

      {do ix^D = ixI^LIM^D \}
         x_tmp         = mesh_i%x(ix^D,:)
         radii_pt      = get_radii_pt(x_tmp)
         eps_D(ix^D)   = eps_fac * small_rho
         eps_tau(ix^D) = eps_fac * small_e
      {end do^D&\}

      ! get theta for D
      if (slab_uniform) then
         inv_lambda(ixI^S) = dxdim(idims)/(qdt)
         flow_dA(ixI^S)  = fC_low(ixI^S,D_)
         fhigh_dA(ixI^S) = fC(ixI^S,D_,idims)
      else
         inv_lambda(ixI^S) = sCT%mesh%dvolume(ixI^S)/(qdt)
         flow_dA(ixI^S)  = fC_low(ixI^S,D_)   * sCT%mesh%surfaceC(ixI^S,idims)
         fhigh_dA(ixI^S) = fC(ixI^S,D_,idims) * sCT%mesh%surfaceC(ixI^S,idims)
      end if
      call get_theta(ixI^L,ixC^L,idims,eps_D,inv_lambda,sCT%w(ixI^S,D_),flow_dA,fhigh_dA,theta)

      ! get theta for tau
      if (slab_uniform) then
         flow_dA(ixI^S)  = fC_low(ixI^S,tau_)
         fhigh_dA(ixI^S) = fC(ixI^S,tau_,idims)
      else
         flow_dA(ixI^S)  = fC_low(ixI^S,tau_)   * sCT%mesh%surfaceC(ixI^S,idims)
         fhigh_dA(ixI^S) = fC(ixI^S,tau_,idims) * sCT%mesh%surfaceC(ixI^S,idims)
      end if
      ! Note that cons(tau_) may be -ve for tabulate EOS
      call get_theta(ixI^L,ixC^L,idims,eps_tau,inv_lambda,sCT%w(ixI^S,tau_),flow_dA,fhigh_dA,theta_tau)

      theta(ixC^S) = min(theta(ixC^S),theta_tau(ixC^S))

      do iw = 1,nwflux
         ! note that pp limiter cannot act on stagger grid
         if ( stagger_grid ) then
            if ( (iw>=Bcons(1)) .and. (iw<=Bcons(ndim)) ) cycle
         end if
         fC(ixC^S,iw,idims) = theta(ixC^S)*fC(ixC^S,iw,idims) &
              + (1.0d0-theta(ixC^S))*fC_low(ixC^S,iw)
      end do ! Next iw
    end subroutine positivity_preserving_limiter

  end subroutine finite_difference

  subroutine reconstructL(typelimiter_in,ixI^L,iL^L,idims,f,fL)
    use mod_global_parameters
    use mod_limiter
    integer, intent(in)             :: typelimiter_in
    integer, intent(in)             :: ixI^L, iL^L, idims
    double precision, intent(in)    :: f(ixI^S,1:nwflux)
    double precision, intent(out)   :: fL(ixI^S,1:nwflux) 
    double precision                :: ldw(ixI^S), dwC(ixI^S)
    integer                         :: jxR^L, ixC^L, jxC^L, kxC^L, iw

    select case (typelimiter_in)
    case (limiter_mp5)
       call MP5limiterL(1,nwflux,nwflux,ixI^L,iL^L,idims,f,fL,.true.)
    case (limiter_weno5)
       call WENO5limiterL(1,nwflux,nwflux,ixI^L,iL^L,idims,f,fL,1,.true.)
    case (limiter_weno5nm)
       call WENO5NMlimiterL(1,nwflux,nwflux,ixI^L,iL^L,idims,f,fL,1,.true.)
    case (limiter_wenoz5)
       call WENO5limiterL(1,nwflux,nwflux,ixI^L,iL^L,idims,f,fL,2,.true.)
    case (limiter_wenoz5nm)
       call WENO5NMlimiterL(1,nwflux,nwflux,ixI^L,iL^L,idims,f,fL,2,.true.)
    case (limiter_wenozp5)
       call WENO5limiterL(1,nwflux,nwflux,ixI^L,iL^L,idims,f,fL,3,.true.)
    case (limiter_wenozp5nm)
       call WENO5NMlimiterL(1,nwflux,nwflux,ixI^L,iL^L,idims,f,fL,3,.true.)
    case default 
       kxCmin^D=ixImin^D; kxCmax^D=ixImax^D-kr(idims,^D);
       fL(kxC^S,1:nwflux) = f(kxC^S,1:nwflux)
       jxR^L=iL^L+kr(idims,^D);
       ixCmax^D=jxRmax^D; ixCmin^D=iLmin^D-kr(idims,^D);
       jxC^L=ixC^L+kr(idims,^D);
       do iw=1,nwflux
          dwC(ixC^S)=f(jxC^S,iw)-f(ixC^S,iw)
          call dwlimiter2(dwC,ixI^L,ixC^L,idims,typelimiter_in,ldw)
          fL(iL^S,iw)=fL(iL^S,iw)+half*ldw(iL^S)
       end do
    end select
  end subroutine reconstructL

  subroutine reconstructR(typelimiter_in,ixI^L,iL^L,idims,f,fR)
    use mod_global_parameters
    use mod_limiter
    integer, intent(in)             :: typelimiter_in
    integer, intent(in)             :: ixI^L, iL^L, idims
    double precision, intent(in)    :: f(ixI^S,1:nwflux)
    double precision, intent(out)   :: fR(ixI^S,1:nwflux) 
    double precision                :: rdw(ixI^S), dwC(ixI^S)
    integer                         :: jxR^L, ixC^L, jxC^L, kxC^L, kxR^L, iw

    select case (typelimiter_in)
    case (limiter_mp5)
       call MP5limiterR(1,nwflux,nwflux,ixI^L,iL^L,idims,f,fR,.true.)
    case (limiter_weno5)
       call WENO5limiterR(1,nwflux,nwflux,ixI^L,iL^L,idims,f,fR,1,.true.)
    case (limiter_weno5nm)
       call WENO5NMlimiterR(1,nwflux,nwflux,ixI^L,iL^L,idims,f,fR,1,.true.)
    case (limiter_wenoz5)
       call WENO5limiterR(1,nwflux,nwflux,ixI^L,iL^L,idims,f,fR,2,.true.)
    case (limiter_wenoz5nm)
       call WENO5NMlimiterR(1,nwflux,nwflux,ixI^L,iL^L,idims,f,fR,2,.true.)
    case (limiter_wenozp5)
       call WENO5limiterR(1,nwflux,nwflux,ixI^L,iL^L,idims,f,fR,3,.true.)
    case (limiter_wenozp5nm)
       call WENO5NMlimiterR(1,nwflux,nwflux,ixI^L,iL^L,idims,f,fR,3,.true.)
    case default 
       kxCmin^D=ixImin^D; kxCmax^D=ixImax^D-kr(idims,^D);
       kxR^L=kxC^L+kr(idims,^D);
       fR(kxC^S,1:nwflux)=f(kxR^S,1:nwflux)
       jxR^L=iL^L+kr(idims,^D);
       ixCmax^D=jxRmax^D; ixCmin^D=iLmin^D-kr(idims,^D);
       jxC^L=ixC^L+kr(idims,^D);
       do iw=1,nwflux
          dwC(ixC^S)=f(jxC^S,iw)-f(ixC^S,iw)
          call dwlimiter2(dwC,ixI^L,ixC^L,idims,typelimiter_in,rdw)
          fR(iL^S,iw)=fR(iL^S,iw)-half*rdw(jxR^S)
       end do
    end select
  end subroutine reconstructR

end module mod_finite_difference
