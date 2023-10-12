!> Module with finite volume methods for fluxes
module mod_finite_volume
  implicit none
  private

  public :: finite_volume
  public :: reconstruct_LR
  public :: metric_interpolation
  public :: get_theta

contains

  !> finite volume method
  ! fixme: replace dx^D to dxs
  subroutine finite_volume(method,qdt,ixI^L,ixO^L,idims^LIM,qtC,sCT,qt,snew,sold,fC,fE,dx^D)
    use mod_physics
    use mod_geometry
    use mod_global_parameters
    use mod_limiter
    use mod_tvd, only:tvdlimit2
    use mod_source, only: addsource2
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
    double precision, dimension(ixI^S,1:nwflux)         :: fLC, fRC
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

       if(stagger_grid) then
         ! ct needs all transverse cells
         ixCmin^D=hxOmin^D-(nghostcells-phys_extra_ghostcells)*(1-kr(idims,^D));
         ixCmax^D=ixOmax^D+(nghostcells-phys_extra_ghostcells)*(1-kr(idims,^D)); 
       else
         ! ixC is centered index in the idims direction from ixOmin-1/2 to ixOmax+1/2
         ixCmin^D=hxOmin^D;
         ixCmax^D=ixOmax^D; 
       end if
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
       if(method==fs_tvdlf.or.method==fs_tvdmu) then
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
       case(fs_hll)
         do ii=1,number_species
           call get_Riemann_flux_hll(start_indices(ii),stop_indices(ii))
         end do
       case(fs_tvdlf)
         do ii=1,number_species
           call get_Riemann_flux_tvdlf(start_indices(ii),stop_indices(ii))
         end do
       case(fs_tvdmu)
         call get_Riemann_flux_tvdmu()
       case default
         call mpistop('unkown Riemann flux in finite volume')
       end select

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

        ! For the MUSCL scheme apply the characteristic based limiter
        if (method==fs_tvdmu) &
            call tvdlimit2(method,qdt,ixI^L,ixC^L,ixO^L,idims,wLC,wRC,w_new,sCT%mesh%x,fC,dx^D)
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
        ! For the MUSCL scheme apply the characteristic based limiter
        if (method==fs_tvdmu) &
            call tvdlimit2(method,qdt,ixI^L,ixC^L,ixO^L,idims,wLC,wRC,w_new,sCT%mesh%x,fC,dx^D)
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

    subroutine get_Riemann_flux_tvdmu()
      do iw=1,nwflux
         ! To save memory we use fLC to store (F_L+F_R)/2=0.5d0*(fLC+fRC)
         fLC(ixC^S, iw)=0.5d0*(fLC(ixC^S, iw)+fRC(ixC^S, iw))
         fC(ixC^S,iw,idims)=fLC(ixC^S, iw)
      end do
    end subroutine get_Riemann_flux_tvdmu

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
         case (flux_asym_diffusion_e)
            {do ix^DB=ixCmin^DB,ixCmax^DB\} 
              if(cminC(ix^D,ii) >= 0.0d0) then
                fC(ix^D,iw,idims)=fLC(ix^D,iw)
              else if(cmaxC(ix^D,ii) <= 0.0d0) then
                fC(ix^D,iw,idims)=fRC(ix^D,iw)
              else
                ! Add hll dissipation to the flux
                fC(ix^D,iw,idims)=(cmaxC(ix^D,ii)*fLC(ix^D,iw)-cminC(ix^D,ii)*fRC(ix^D,iw)&
                      +tvdlf_eps(ix^D,ii)*tvdlfeps*cminC(ix^D,ii)&
                       *cmaxC(ix^D,ii)*(sR%w(ix^D,iw)-sL%w(ix^D,iw)))&
                      /(cmaxC(ix^D,ii)-cminC(ix^D,ii))
              end if
            {end do\}
         case (flux_asym_diffusion_f)
            {do ix^DB=ixCmin^DB,ixCmax^DB\} 
              if(cminC(ix^D,ii) >= 0.0d0) then
                fC(ix^D,iw,idims)=fLC(ix^D,iw)
              else if(cmaxC(ix^D,ii) <= 0.0d0) then
                fC(ix^D,iw,idims)=fRC(ix^D,iw)
              else
                ! Add hll dissipation to the flux
                fC(ix^D,iw,idims)= (1.0d0-tvdlf_eps(ix^D,ii)**2)*0.5d0*(fLC(ix^D,iw)+fRC(ix^D,iw)) &
                  + ( tvdlf_eps(ix^D,ii)**2*(cmaxC(ix^D,ii)*fLC(ix^D,iw)-cminC(ix^D,ii)*fRC(ix^D,iw)) &
                      +tvdlf_eps(ix^D,ii)*tvdlfeps*cminC(ix^D,ii)*cmaxC(ix^D,ii) &
                        *(sR%w(ix^D,iw)-sL%w(ix^D,iw)) ) &
                         /(cmaxC(ix^D,ii)-cminC(ix^D,ii))
              end if
            {end do\}
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

    ! fixme: make this limiter as a standalone subroutine
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

  end subroutine finite_volume

  subroutine get_theta(ixI^L,ixC^L,idims,eps,inv_lambda,u,flow_dA,fhigh_dA,theta)
    use mod_global_parameters, only: kr, smalldouble
    integer, intent(in)                             :: ixI^L, ixC^L, idims
    double precision, dimension(ixI^S), intent(in)  :: eps
    double precision, dimension(ixI^S), intent(in)  :: inv_lambda, u
    double precision, dimension(ixI^S), intent(in)  :: flow_dA, fhigh_dA
    double precision, dimension(ixI^S), intent(out) :: theta

    integer                                         :: ixCp^L, ixOp^L
    double precision, dimension(ixI^S)              :: tmp, thp, thm
    double precision, dimension(ixI^S)              :: diff_fdA

    ! Note: here we assume that u( i=0 ) is given
    ixCp^L=ixC^L+kr(idims,^D);
    ixOpmin^D=ixCmin^D; ixOpmax^D=ixCpmax^D;
    
    thm(ixC^S) = 1.0d0
    thp(ixC^S) = 1.0d0
    
    tmp(ixOp^S) = 0.5d0*inv_lambda(ixOp^S)*(u(ixOp^S)-eps(ixOp^S))/dble(^ND)

    diff_fdA(ixC^S) = -flow_dA(ixC^S) + fhigh_dA(ixC^S)
    where (diff_fdA(ixC^S) == 0.0d0)
       diff_fdA(ixC^S) = smalldouble ! avoid flow = fhight case
    end where
    
    where (tmp(ixC^S) < fhigh_dA(ixC^S))
       thm(ixC^S) = tmp(ixC^S) - flow_dA(ixC^S)
       thm(ixC^S) = thm(ixC^S) / (diff_fdA(ixC^S))
    end where
    
    where (tmp(ixCp^S) < -fhigh_dA(ixC^S))
       thp(ixC^S) = - tmp(ixCp^S) - flow_dA(ixC^S)
       thp(ixC^S) = thp(ixC^S) / (diff_fdA(ixC^S))
    end where

    theta(ixC^S) = min(thm(ixC^S),thp(ixC^S))
    theta(ixC^S) = min(max(theta(ixC^S),0.0d0),1.0d0)
  end subroutine get_theta

  !> Determine the upwinded consL(ixL) and consR(ixR) from w.
  subroutine reconstruct_LR(typelimiter_in,ixI^L,ixL^L,ixR^L,idims,s,sL,sR,sLp,sRp,dxdim)
      ! fixme: remove dxdim?
    use mod_physics
    use mod_geometry
    use mod_global_parameters
    use mod_limiter

    integer, intent(in)          :: typelimiter_in
    integer, intent(in)          :: ixI^L, ixL^L, ixR^L, idims
    ! cell center w in primitive form
    type(state), intent(inout)   :: s
    ! left and right constructed status in conservative form
    type(state), intent(inout)   :: sL, sR
    ! left and right constructed status in primitive form
    type(state), intent(inout)   :: sLp, sRp
    double precision, intent(in) :: dxdim

    integer            :: jxR^L, ixC^L, jxC^L, iw
    double precision   :: ldw(ixI^S), rdw(ixI^S), dwC(ixI^S)

    associate( x   => s%mesh%x, xi => sL%mesh%x, &
               w   => s%w,                       &
               wLC => sL%w,     wRC => sR%w,     &
               wLp => sLp%w,    wRp => sRp%w     &
             )

    do iw=1, nwflux
       select case (typelimiter_in)
       case (limiter_pc)
         ! piecewise constant, nothing to do here
       case (limiter_venk)
          call venklimiter(ixI^L,ixL^L,idims,dxdim,w(ixI^S,iw),wLp(ixI^S,iw),wRp(ixI^S,iw)) 
       case (limiter_mp5)
          call MP5limiter(ixI^L,ixL^L,idims,w(ixI^S,iw),wLp(ixI^S,iw),wRp(ixI^S,iw),.false.)
       case (limiter_weno3)
          call WENO3limiter(ixI^L,ixL^L,idims,dxdim,w(ixI^S,iw),wLp(ixI^S,iw),wRp(ixI^S,iw),1,.false.)
       case (limiter_wenoyc3)
          call WENO3limiter(ixI^L,ixL^L,idims,dxdim,w(ixI^S,iw),wLp(ixI^S,iw),wRp(ixI^S,iw),2,.false.)
       case (limiter_weno5)
          call WENO5limiter(ixI^L,ixL^L,idims,dxdim,w(ixI^S,iw),wLp(ixI^S,iw),wRp(ixI^S,iw),1,.false.)
       case (limiter_weno5nm)
          call WENO5NMlimiter(ixI^L,ixL^L,idims,dxdim,w(ixI^S,iw),wLp(ixI^S,iw),wRp(ixI^S,iw),1,.false.)
       case (limiter_wenoz5)
          call WENO5limiter(ixI^L,ixL^L,idims,dxdim,w(ixI^S,iw),wLp(ixI^S,iw),wRp(ixI^S,iw),2,.false.)
       case (limiter_wenoz5nm)
          call WENO5NMlimiter(ixI^L,ixL^L,idims,dxdim,w(ixI^S,iw),wLp(ixI^S,iw),wRp(ixI^S,iw),2,.false.)
       case (limiter_wenozp5)
          call WENO5limiter(ixI^L,ixL^L,idims,dxdim,w(ixI^S,iw),wLp(ixI^S,iw),wRp(ixI^S,iw),3,.false.)
       case (limiter_wenozp5nm)
          call WENO5NMlimiter(ixI^L,ixL^L,idims,dxdim,w(ixI^S,iw),wLp(ixI^S,iw),wRp(ixI^S,iw),3,.false.)
       case (limiter_weno5cu6)
          call WENO5CU6limiter(ixI^L,ixL^L,idims,w(ixI^S,iw),wLp(ixI^S,iw),wRp(ixI^S,iw),.false.)
       case (limiter_teno5ad)
          call TENO5ADlimiter(ixI^L,ixL^L,idims,dxdim,w(ixI^S,iw),wLp(ixI^S,iw),wRp(ixI^S,iw),.false.)
       case (limiter_weno7)
          call WENO7limiter(ixI^L,ixL^L,idims,w(ixI^S,iw),wLp(ixI^S,iw),wRp(ixI^S,iw),1,.false.)
       case (limiter_mpweno7)
          call WENO7limiter(ixI^L,ixL^L,idims,w(ixI^S,iw),wLp(ixI^S,iw),wRp(ixI^S,iw),2,.false.)
       case (limiter_exeno7)
          call exENO7limiter(ixI^L,ixL^L,idims,w(ixI^S,iw),wLp(ixI^S,iw),wRp(ixI^S,iw),.false.)
       case (limiter_ppm)
          call PPMlimiter(ixI^L,ixL^L,idims,w(ixI^S,iw),wLp(ixI^S,iw),wRp(ixI^S,iw),PPM_extrema)
       case default
          jxR^L=ixR^L+kr(idims,^D);
          ixCmax^D=jxRmax^D; ixCmin^D=ixLmin^D-kr(idims,^D);
          jxC^L=ixC^L+kr(idims,^D);

          if (loglimit(iw)) then
             w(ixCmin^D:jxCmax^D,iw)=dlog10(w(ixCmin^D:jxCmax^D,iw))
             wLp(ixL^S,iw)=dlog10(wLp(ixL^S,iw))
             wRp(ixR^S,iw)=dlog10(wRp(ixR^S,iw))
          end if

          ! ----- original version ------
          dwC(ixC^S)=( w(jxC^S,iw)-w(ixC^S,iw) )
          ! limit flux from left and/or right
          call dwlimiter2(dwC,ixI^L,ixC^L,idims,typelimiter_in,ldw,rdw)
          wLp(ixL^S,iw)=wLp(ixL^S,iw)+0.5d0*ldw(ixL^S)
          wRp(ixR^S,iw)=wRp(ixR^S,iw)-0.5d0*rdw(jxR^S)
          ! ----------------------------

          ! ----- bary center version ------
          !!dxi(ixC^S) = dxdim!xi(jxC^S,idims) - xi(ixC^S,idims)
          !dwC(ixC^S)=( w(jxC^S,iw)-w(ixC^S,iw) ) &
          !            / (xbar(jxC^S,idims) - xbar(ixC^S,idims)) &
          !             * dxdim
          !
          !! limit flux from left and/or right
          !call dwlimiter2(dwC,ixI^L,ixC^L,idims,typelimiter_in,ldw,rdw)
          !wLp(ixL^S,iw)=w(ixL^S,iw)+ldw(ixL^S)*(xi(ixL^S,idims)-xbar(ixL^S,idims))/dxdim
          !wRp(ixR^S,iw)=w(jxR^S,iw)+rdw(jxR^S)*(xi(ixR^S,idims)-xbar(jxR^S,idims))/dxdim
          ! ----------------------------

          if (loglimit(iw)) then
             w(ixCmin^D:jxCmax^D,iw)=10.0d0**w(ixCmin^D:jxCmax^D,iw)
             wLp(ixL^S,iw)=10.0d0**wLp(ixL^S,iw)
             wRp(ixR^S,iw)=10.0d0**wRp(ixR^S,iw)
          end if
       end select

    end do

    if ( fix_small_values ) then
       call phys_handle_small_values(sLp,ixI^L,ixL^L,'reconstruction L')
       call phys_handle_small_values(sRp,ixI^L,ixR^L,'reconstruction R')
    end if

    !if (typelimiter_in == limiter_ppm) &
    !      call PPMflattening(ixI^L,ixL^L,idims,w,wL,wR,PPM_flatcd,PPM_flatsh)

    ! convert to conserved variables
    sL = sLp; sR = sRp
    call phys_to_conserved(ixI^L,ixL^L,sL)
    call phys_to_conserved(ixI^L,ixR^L,sR)
   
    end associate
  end subroutine reconstruct_LR

  subroutine metric_interpolation(ixI^L,ixO^L,idims,nvar,vars,x,varsi,xi)
    use mod_global_parameters
    use mod_interpolation

    integer, intent(in) :: ixI^L, ixO^L, idims, nvar
    double precision, dimension(ixI^S,1:nvar), intent(in)  :: vars
    double precision, dimension(ixI^S,1:nvar), intent(out) :: varsi
    double precision, dimension(ixI^S,1:ndim), intent(in)  :: x, xi

    ! local vars
    integer :: ix^D, iw
    integer :: n_lo, n_hi
    double precision, allocatable :: xx(:), yy(:)

    if ( mod(fv_n_interp,2) == 0 ) then
       n_lo = fv_n_interp / 2 - 1
       n_hi = fv_n_interp / 2
    else
       n_lo = (fv_n_interp+1) / 2 - 1
       n_hi = (fv_n_interp-1) / 2
    end if

    allocate(xx(fv_n_interp))
    allocate(yy(fv_n_interp))

    do iw = 1, nvar

       {^IFONED
       {do ix^D = ixO^LIM^D \}
          xx = x(ix1-n_lo:ix1+n_hi, 1)
          yy = vars(ix1-n_lo:ix1+n_hi, iw)
          call lagrange_interpolation(xx,yy,&
                  xi(ix^D, idims), varsi(ix^D, iw), (xi(ix^D, idims)<0.0d0))
       {end do^D&\}
       }
       {^IFTWOD
       select case (idims)
       case (1)
          {do ix^D = ixO^LIM^D \}
             xx = x(ix1-n_lo:ix1+n_hi, ix2, 1)
             yy = vars(ix1-n_lo:ix1+n_hi, ix2, iw)
             call lagrange_interpolation(xx,yy,&
                     xi(ix^D, idims), varsi(ix^D, iw), (xi(ix^D, idims)<0.0d0))
          {end do^D&\}
       case (2)
          {do ix^D = ixO^LIM^D \}
             xx = x(ix1, ix2-n_lo:ix2+n_hi, 2)
             yy = vars(ix1, ix2-n_lo:ix2+n_hi, iw)
             call lagrange_interpolation(xx,yy,&
                     xi(ix^D, idims), varsi(ix^D, iw), (xi(ix^D, idims)<0.0d0))
          {end do^D&\}
       case default
          call mpistop(" idims can only be 1 or 2 in 2D")
       end select
       }
       {^IFTHREED
       select case (idims)
       case (1)
          {do ix^D = ixO^LIM^D \}
             xx = x(ix1-n_lo:ix1+n_hi, ix2, ix3, 1)
             yy = vars(ix1-n_lo:ix1+n_hi, ix2, ix3, iw)
             call lagrange_interpolation(xx,yy,&
                     xi(ix^D, idims), varsi(ix^D, iw), (xi(ix^D, idims)<0.0d0))
          {end do^D&\}
       case (2)
          {do ix^D = ixO^LIM^D \}
             xx = x(ix1, ix2-n_lo:ix2+n_hi, ix3, 2)
             yy = vars(ix1, ix2-n_lo:ix2+n_hi, ix3, iw)
             call lagrange_interpolation(xx,yy,&
                     xi(ix^D, idims), varsi(ix^D, iw), (xi(ix^D, idims)<0.0d0))
          {end do^D&\}
       case (3)
          {do ix^D = ixO^LIM^D \}
             xx = x(ix1, ix2, ix3-n_lo:ix3+n_hi, 3)
             yy = vars(ix1, ix2, ix3-n_lo:ix3+n_hi, iw)
             call lagrange_interpolation(xx,yy,&
                     xi(ix^D, idims), varsi(ix^D, iw), (xi(ix^D, idims)<0.0d0))
          {end do^D&\}
       case default
          call mpistop(" idims can only be 1, 2 or 3 in 3D")
       end select
       }
    end do

    deallocate(xx)
    deallocate(yy)
  end subroutine metric_interpolation

end module mod_finite_volume
