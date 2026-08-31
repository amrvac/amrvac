!> Module for flux limited diffusion (FLD)-approximation in Radiation-(Magneto)hydrodynamics simulations
!>
!> Full description of RHD-FLD in
!> Moens N., Sundqvist J.O., El Mellah I., Poniatowski L., Teunissen J. & Keppens R. 2022, A&A 657, A81
!> Radiation-hydrodynamics with MPI-AMRVAC . Flux-limited diffusion
!> doi:10.1051/0004-6361/202141023
!>
!> Full description for RMHD-FLD in
!> N. Narechania, R. Keppens, A. ud-Doula, N. Moens & J. Sundqvist 2025, A&A 696, A131
!> doi:10.1051/0004-6361/202452208
!> Radiation-magnetohydrodynamics with MPI-AMRVAC using flux-limited diffusion

module mod_fld
    use mod_comm_lib, only: mpistop
    use mod_geometry
    implicit none
    !> source split for energy interact and radforce:
    logical :: fld_Radforce_split = .false.
    !> switch to handle photon tiring explicit or implicit
    logical :: fld_tiring_explicit = .false.
    !> Opacity value when using constant opacity
    double precision, public :: fld_kappa0 = 0.0d0
    !> Tolerance for bisection method for Energy sourceterms
    !> This is a percentage of the minimum of gas- and radiation energy
    double precision, public :: fld_bisect_tol = 1.d-4
    !> Tolerance for radiative Energy diffusion
    double precision, public :: fld_diff_tol = 1.d-4
    !> handling ramp-up phase with explicit diffusion dt limit, slowly boosted
    logical :: fld_slowsteps = .false.
    double precision, public :: fld_boost_dt = 0.0d0
    !> switches for using changed cmax-cmin bounds
    logical :: fld_bound_diff = .true.
    !> switches for local debug purposes
    logical :: fld_force_MG_converged = .true.
    logical :: fld_debug,fld_no_mg
    double precision, public :: fld_cnorm = 0.0d0
    !> switches for opacity
    character(len=40)  :: fld_opacity_law = 'const'
    character(len=40) :: fld_opal_table = 'Y09800' 
    !> flux limiter choice
    character(len=40) :: fld_fluxlimiter = 'Pomraning'
    !> diffusion coefficient for multigrid method
    integer :: i_diff_mg
    !> diffusion coefficient stencil control
    integer :: nth_for_diff_mg
    !> Which method to find the root for the energy interaction polynomial
    character(len=40) :: fld_interaction_method = 'Halley'
    !> Abstract interface for the gas-EoS getters the radiation fluid needs
    !> (same shape as thermal_conduction's get_var_subr).
    abstract interface
      subroutine fld_get_var(w,x,ixI^L,ixO^L,res)
        use mod_global_parameters
        integer, intent(in)          :: ixI^L, ixO^L
        double precision, intent(in) :: w(ixI^S,nw)
        double precision, intent(in) :: x(ixI^S,1:ndim)
        double precision, intent(out):: res(ixI^S)
      end subroutine fld_get_var
    end interface

    !> Radiation fluid object: gas-EoS callbacks the FLD module needs, wired by
    !> the physics module at link time (mirrors tc_fluid / rc_fluid). The
    !> instance lives in mod_(m)hd_phys and is threaded in here as `fl`.
    type fld_fluid
      !> adiabatic index (constant for FI; a snapshot of eos%gamma)
      double precision :: gamma
      !> gas temperature from primitive pressure, T = p/(R*rho)
      procedure(fld_get_var), pointer, nopass :: get_tgas => null()
      !> R factor (mean-molecular-weight gas constant) for the emission term
      procedure(fld_get_var), pointer, nopass :: get_Rfactor => null()
    end type fld_fluid

    !> public methods
    !> these are called in mod_hd_phys or mod_mhd_phys
    public :: fld_fluid
    public :: fld_init
    public :: fld_get_radpress
    public :: fld_get_diffcoef_central
    public :: add_fld_rad_force
    public :: fld_radforce_get_dt
    public :: fld_get_local_invtauc
    !> wired to physics-module wrappers (which inject fld_fl) in mod_(m)hd_phys
    public :: fld_implicit_update
    public :: fld_evaluate_implicit
    public :: fld_set_mg_bounds
    ! these are made public for mod_usr purposes and diagnostics
    public :: fld_get_radflux
    public :: fld_get_fluxlimiter
    public :: fld_get_fluxlimiter_prim
    public :: fld_get_opacity_prim
  contains

  !> Reading in fld-list parameters from .par file
  subroutine fld_params_read(files)
    use mod_global_parameters, only: unitpar
    use mod_constants
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /fld_list/ fld_kappa0, fld_Radforce_split, &
    fld_bisect_tol, fld_diff_tol, fld_opacity_law, fld_fluxlimiter, &
    fld_interaction_method, fld_opal_table, nth_for_diff_mg, &
    fld_debug,fld_no_mg,fld_force_MG_converged, fld_slowsteps,fld_boost_dt, &
    fld_tiring_explicit,fld_cnorm,fld_bound_diff

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, fld_list, end=111)
       111    close(unitpar)
    end do
  end subroutine fld_params_read

  !> Initialising FLD-module
  !> Read opacities
  !> Initialise Multigrid and adimensionalise kappa
  subroutine fld_init()
    use mod_global_parameters
    use mod_variables
    use mod_physics
    use mod_opal_opacity, only: init_opal_table
    use mod_multigrid_coupling

    nth_for_diff_mg=1
    fld_debug=.false.
    fld_no_mg=.false.
    ! initialize constant opacity with free electron Thomson scattering value
    fld_kappa0=const_kappae
    call fld_params_read(par_files)
    ! sanity checks on input
    if(fld_kappa0<smalldouble)then
       if(mype==0) print *,'fld_kappa0=',fld_kappa0
       call mpistop("please set the constant opacity to a reasonable value")
    endif
    if(fld_bisect_tol<smalldouble)then
       if(mype==0) print *,'fld_bisect_tol=',fld_bisect_tol
       call mpistop("convergence tolerance for root solver too strict")
    endif
    if(fld_tiring_explicit)then
       if(mype==0) print *,'Will do photon tiring explicit!!!'
    endif
    if(.not.fld_no_mg)then
       if(fld_diff_tol<smalldouble)then
         if(mype==0) print *,'fld_diff_tol=',fld_diff_tol
         call mpistop("convergence tolerance for MG solver too strict")
      endif
      select case(nth_for_diff_mg)
        case(1)
         ! no need for stencil extension
        case(2)
         ! need for stencil extension
         phys_wider_stencil=1 
       case default
         call mpistop("nth_for_diff_mg must be 1 or 2")
       end select
       ! phys_implicit_update / phys_evaluate_implicit are wired in the physics
       ! module (mod_(m)hd_phys) to wrappers that inject the fld_fl object;
       ! set_mg_bounds needs no fluid object so it is bound here directly.
       phys_set_mg_bounds     => fld_set_mg_bounds
       ! store the diffusion coefficient as extra variable (needed in mg vhelmholtz)
       i_diff_mg = var_set_extravar("D", "D")
       use_multigrid = .true.
       ! use multigrid to solve a helmholtz equation with variable coefficient
       ! this is stored also as extra variable in the mg solver as mg_iveps
       mg%n_extra_vars = 1
       mg%operator_type = mg_vhelmholtz
       ! choice of smoother: can be mg_smoother_gs or gsrb (latter recommended)
       mg%smoother_type = mg_smoother_gsrb
    endif
    !> Read in opacity table if necesary
    if(trim(fld_opacity_law) .eq. 'opal') then
      if(SI_unit)call mpistop("adjust opal module with SI-cgs conversions for SI - or use cgs!")
      call init_opal_table(fld_opal_table)
    endif
  end subroutine fld_init

  !> Set the boundaries for the diffusion of E
  subroutine fld_set_mg_bounds
    use mod_global_parameters
    use mod_multigrid_coupling
    use mod_usr_methods

    integer :: iB

    ! Set boundary conditions for the multigrid solver
    do iB = 1, 2*ndim
       select case (typeboundary(iw_r_e, iB))
       case (bc_symm)
          ! d/dx u = 0
          mg%bc(iB, mg_iphi)%bc_type = mg_bc_neumann
          mg%bc(iB, mg_iphi)%bc_value = 0.0_dp
       case (bc_asymm)
          ! u = 0
          mg%bc(iB, mg_iphi)%bc_type = mg_bc_dirichlet
          mg%bc(iB, mg_iphi)%bc_value = 0.0_dp
       case (bc_cont)
          ! d/dx u = 0
          mg%bc(iB, mg_iphi)%bc_type = mg_bc_neumann
          mg%bc(iB, mg_iphi)%bc_value = 0.0_dp
       case (bc_periodic)
          ! Nothing to do here, this is picked up through periodB variable
       case (bc_special)
          if(mype==0)then
             print *,'Special boundary for Erad needs specific user-set MG BC treatment'
             print *,'                 and this could be through usr_special_mg_bc call'
          endif
          if (associated(usr_special_mg_bc)) then 
             call usr_special_mg_bc(iB)
          endif
       case default
          call mpistop("divE_multigrid warning: unknown b.c. ")
       end select
       ! Neumann on diffusion coefficient is needed on all lower grid levels
       ! d/dx u = 0
       mg%bc(iB, mg_iveps)%bc_type = mg_bc_neumann
       mg%bc(iB, mg_iveps)%bc_value = 0.0_dp
    end do

  end subroutine fld_set_mg_bounds


  !> w[iw]=w[iw]+qdt*S[wCT,qtC,x] where S is the source based on wCT within ixO
  !> This subroutine handles the radiation force and its work added explicitly
  !> and the energy interaction term combined with photon tiring using an implicit update
  subroutine add_fld_rad_force(qdt,ixI^L,ixO^L,wCT,wCTprim,w,x,qsourcesplit,active,fl)
    use mod_constants
    use mod_global_parameters
    use mod_geometry
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, x(ixI^S,1:ndim)
    double precision, intent(in)    :: wCT(ixI^S,1:nw),wCTprim(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    logical, intent(in) :: qsourcesplit
    logical, intent(inout) :: active
    type(fld_fluid), intent(in)     :: fl

    integer :: idir,jdir,nth_for_fld,ix^D
    double precision, dimension(ixI^S) :: a1,a2,a3,c0,c1,kappa
    double precision, dimension(ixI^S) :: e_gas,E_rad,tmp
    double precision, dimension(ixI^S,1:ndim,1:ndim) :: div_v,edd

    !> Calculate and add sourceterms
    if(qsourcesplit .eqv. fld_Radforce_split) then
      active = .true.
      nth_for_fld=2
      ! store here lambda in a1 and fld_R in a2
      call fld_get_eddington(wCTprim,x,ixI^L,ixO^L,edd,a1,a2,nth_for_fld,fl)

      !> Photon tiring : calculate tensor grad v (named div_v here)
      ! NOTE: This is ok for uniform Cartesian only!!!!!
      ! TODO: introduce gradient of vector in geometry module and call that one
      do idir = 1,ndim
        do jdir = 1,ndim
          call gradient(wCTprim(ixI^S,iw_mom(jdir)),ixI^L,ixO^L,idir,tmp)
          div_v(ixO^S,idir,jdir) = tmp(ixO^S)
        enddo
      enddo
      ! perform contraction fe : grad(v) with fe eddington tensor
      {^IFONED
      a3(ixO^S) = div_v(ixO^S,1,1)*edd(ixO^S,1,1)
      }
      {^IFTWOD
      a3(ixO^S) = div_v(ixO^S,1,1)*edd(ixO^S,1,1) &
                      + div_v(ixO^S,1,2)*edd(ixO^S,1,2) &
                      + div_v(ixO^S,2,1)*edd(ixO^S,2,1) &
                      + div_v(ixO^S,2,2)*edd(ixO^S,2,2)
      }
      {^IFTHREED
      a3(ixO^S) = div_v(ixO^S,1,1)*edd(ixO^S,1,1) &
                      + div_v(ixO^S,1,2)*edd(ixO^S,1,2) &
                      + div_v(ixO^S,1,3)*edd(ixO^S,1,3) &
                      + div_v(ixO^S,2,1)*edd(ixO^S,2,1) &
                      + div_v(ixO^S,2,2)*edd(ixO^S,2,2) &
                      + div_v(ixO^S,2,3)*edd(ixO^S,2,3) &
                      + div_v(ixO^S,3,1)*edd(ixO^S,3,1) &
                      + div_v(ixO^S,3,2)*edd(ixO^S,3,2) &
                      + div_v(ixO^S,3,3)*edd(ixO^S,3,3)
      }

      do idir = 1,ndim
        call gradient(wCTprim(ixI^S,iw_r_e),ixI^L,ixO^L,idir,tmp,nth_for_fld)
        ! Radiation force = kappa*rho/c *Flux = lambda gradE
        ! recycle grad E to store -lambda (grad E)_i
        tmp(ixO^S) = -a1(ixO^S)*tmp(ixO^S)
        !> Momentum equation source term
        w(ixO^S,iw_mom(idir)) = w(ixO^S,iw_mom(idir))+ qdt*tmp(ixO^S)
        !> Energy equation source term 
        w(ixO^S,iw_e) = w(ixO^S,iw_e) + qdt*wCTprim(ixO^S,iw_mom(idir))*tmp(ixO^S)
      enddo
      !> photon tiring when handled explicitly
      if(fld_tiring_explicit)then
          w(ixO^S,iw_r_e) = w(ixO^S,iw_r_e) - qdt*wCTprim(ixO^S,iw_r_e)*a3(ixO^S)
          a3(ixO^S)=0.0d0
      endif

      call get_and_check_egas_Erad_from_conserved(w,ixI^L,ixO^L,e_gas,E_rad)
      ! BEGIN radiative exchange part with or without photon tiring included
         call fld_get_opacity_prim(wCTprim,x,ixI^L,ixO^L,kappa,fl)
         !> Coefficients for the polynomial in Moens et al. 2022, eq 37. but with photon tiring (a3)
         !> FI emission term: a1 carries arad*T^4 with T=(gamma-1)*e_gas/(rho*R) folded
         !> into the e_gas^4 root-find below. gamma and R come from the EoS via fl.
         call fl%get_Rfactor(wCT,x,ixI^L,ixO^L,tmp)
         a1(ixO^S) = qdt*kappa(ixO^S)*c_norm*arad_norm*(fl%gamma-one)**4/(wCT(ixO^S,iw_rho)**3*tmp(ixO^S)**4)
         a2(ixO^S) = c_norm*kappa(ixO^S)*wCT(ixO^S,iw_rho)*qdt
         a3(ixO^S) = a3(ixO^S)*qdt

         c0(ixO^S) = ((one+a2(ixO^S)+a3(ixO^S))*e_gas(ixO^S)+a2(ixO^S)*E_rad(ixO^S))/a1(ixO^S)/(one+a3(ixO^S))
         c1(ixO^S) = (one+a2(ixO^S)+a3(ixO^S))/a1(ixO^S)/(one+a3(ixO^S))

         !> Loop over every cell for rootfinding method
         {do ix^D = ixOmin^D,ixOmax^D\}
           select case(fld_interaction_method)
           case('Bisect')
            call Bisection_method(e_gas(ix^D),c0(ix^D),c1(ix^D))
           case('Newton')
            call Newton_method(e_gas(ix^D),c0(ix^D),c1(ix^D))
           case('Halley')
            call Halley_method(e_gas(ix^D),c0(ix^D),c1(ix^D))
           case default
            call mpistop('root-method not known')
           end select 
         {enddo\}

         E_rad(ixO^S) = (a1(ixO^S)*e_gas(ixO^S)**4.d0+E_rad(ixO^S))/(one+a2(ixO^S)+a3(ixO^S))

         if(check_small_values.and..not.fix_small_values)then
          {do ix^DB= ixOmin^DB,ixOmax^DB\}
           if(e_gas(ix^D)<small_e.or.E_rad(ix^D)<small_r_e) then
              write(*,*) "Error in FLD add_fld_rad_force: small value"
              write(*,*) "of internal or radiation energy density after exchange"
              write(*,*) "Iteration: ", it, " Time: ", global_time
              write(*,*) "Location: ", x(ix^D,:)
              write(*,*) "Cell number: ", ix^D
              write(*,*) "internal energy density  is=",e_gas(ix^D)," versus   small_e=",small_e
              write(*,*) "radiation energy density is=",E_rad(ix^D)," versus small_r_e=",small_r_e
              call mpistop("FLD error:May need to turn on fixes")
           end if
          {end do\}
         endif

        if(fix_small_values)then
          {do ix^D = ixOmin^D,ixOmax^D\ }
            e_gas(ix^D) = max(e_gas(ix^D),small_e)
            E_rad(ix^D) = max(E_rad(ix^D),small_r_e)
          {enddo\}
        endif

      ! END radiative exchange part with or without photon tiring included

      !> Update gas-energy in w, internal + kinetic
      w(ixO^S,iw_e) = e_gas(ixO^S)
      w(ixO^S,iw_e) = w(ixO^S,iw_e)+half*sum(w(ixO^S,iw_mom(:))**2,dim=ndim+1)/w(ixO^S,iw_rho)
      if(allocated(iw_mag)) then
        w(ixO^S,iw_e) = w(ixO^S,iw_e)+half*sum(w(ixO^S,iw_mag(:))**2,dim=ndim+1)
      endif
      !> Update rad-energy in w
      w(ixO^S,iw_r_e) = E_rad(ixO^S)
    end if
  end subroutine add_fld_rad_force

  !> maybe needed for restricting the cmax/cmin in diffusion limit
  !> CURRENTLY UNUSED
  !> NOTE: w is primitive on entry
  subroutine fld_get_local_invtauc(w,ixI^L,ixO^L,dx^D,x,invtauc,fl)
    use mod_global_parameters
    use mod_usr_methods
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: dx^D, x(ixI^S,1:ndim), w(ixI^S,1:nw)
    double precision, intent(out)   :: invtauc(ixI^S)
    type(fld_fluid), intent(in)     :: fl

    double precision :: dxinv(1:ndim)
    double precision :: kappa(ixI^S)

    call fld_get_opacity_prim(w, x, ixI^L, ixO^L, kappa, fl)
    ^D&dxinv(^D)=one/dx^D;
    invtauc(ixO^S)=4.0d0*minval(dxinv(1:ndim))/(3.0d0*kappa(ixO^S)*w(ixO^S,iw_rho))
  end subroutine fld_get_local_invtauc

  !> get dt limit for radiation force and FLD explicit source additions
  !> NOTE: w is primitive on entry
  subroutine fld_radforce_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x,fl)
    use mod_global_parameters
    use mod_usr_methods
    use mod_geometry
    use mod_physics, only: phys_get_csrad2

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: dx^D, x(ixI^S,1:ndim), w(ixI^S,1:nw)
    double precision, intent(inout) :: dtnew
    type(fld_fluid), intent(in)     :: fl

    integer          :: idim,idims,nth_for_fld
    double precision :: dxinv(1:ndim), max_grav
    double precision :: lambda(ixI^S),fld_R(ixI^S)
    double precision :: tmp(ixI^S),boostfactor
    double precision :: max_diff,max_diff_coef,max_diff_dim,dtdifflimit
    double precision :: cmax(ixI^S),cmaxtot(ixI^S),courantmaxtots

    if(fld_debug)print *,'DT limit on entry to radforce_get_dt=',dtnew
    nth_for_fld=2
    call fld_get_fluxlimiter_prim(w,x,ixI^L,ixO^L,lambda,fld_R,nth_for_fld,fl)
    if(slab_uniform) then
      ^D&dxinv(^D)=one/dx^D;
      do idim = 1, ndim
        call gradient(w(ixI^S,iw_r_e),ixI^L,ixO^L,idim,tmp,nth_for_fld)
        max_grav = maxval(dabs(-lambda(ixO^S)*tmp(ixO^S)/w(ixO^S,iw_rho)))
        max_grav = max(max_grav, epsilon(1.0d0))
        dtnew = min(dtnew, 1.0d0 / dsqrt(max_grav * dxinv(idim)))
      end do
    else
      do idim = 1, ndim
        call gradient(w(ixI^S,iw_r_e),ixI^L,ixO^L,idim,tmp,nth_for_fld)
        max_grav = maxval(dabs(-lambda(ixO^S)*tmp(ixO^S)/w(ixO^S,iw_rho))/block%ds(ixO^S,idim))
        max_grav = max(max_grav, epsilon(1.0d0))
        dtnew = min(dtnew, 1.0d0 / dsqrt(max_grav))
      end do
    endif
    if(fld_debug)print *,'DT limit after RADFORCE eff grav=',dtnew

    if(fld_slowsteps) then
      boostfactor=1.0d0
      call fld_get_opacity_prim(w, x, ixI^L, ixO^L, tmp, fl)
      max_diff=0.0d0
      if(slab_uniform) then
         max_diff_coef = maxval(c_norm*lambda(ixO^S)/(w(ixO^S,iw_rho)*tmp(ixO^S)))
         ^D&dxinv(^D)=one/dx^D;
         do idim = 1, ndim
            max_diff_dim = max(2.0d0*ndim*max_diff_coef*dxinv(idim)**2, epsilon(1.0d0))
            max_diff = max(max_diff,max_diff_dim)
         end do
      else
         do idim = 1, ndim
            max_diff_coef = maxval(c_norm*lambda(ixO^S)/(w(ixO^S,iw_rho)*tmp(ixO^S))/block%ds(ixO^S,idim)**2)
            max_diff_dim = max(2.0d0*ndim*max_diff_coef, epsilon(1.0d0))
            max_diff = max(max_diff,max_diff_dim)
         end do
      endif
      dtdifflimit=1.0d0/max_diff
      if(fld_debug) print *,'DT limit from diffusion is actually=',dtdifflimit
      if(slowsteps>it-it_init+1) then
         boostfactor=1.0d0+fld_boost_dt*(dble(it-it_init+1)/dble(slowsteps))
      else
         boostfactor=1.0d0+fld_boost_dt
      endif
      dtdifflimit=dtdifflimit*boostfactor
      if(fld_debug) print *,'DT limit from diffusion is boosted to=',dtdifflimit
      dtnew=min(dtnew,dtdifflimit)
      if(fld_debug) print *,'DT limit after diffusion is =',dtnew
    endif

    ! here we interface back to fld_get_radpress
    call phys_get_csrad2(w,x,ixI^L,ixO^L,tmp)
    if(slab_uniform) then
       ^D&dxinv(^D)=one/dx^D;
       do idims=1,ndim
          cmax(ixO^S)=dabs(w(ixO^S,iw_mom(idims)))+dsqrt(tmp(ixO^S))
          if(idims==1) then
            cmaxtot(ixO^S)=cmax(ixO^S)*dxinv(idims)
          else
            cmaxtot(ixO^S)=cmaxtot(ixO^S)+cmax(ixO^S)*dxinv(idims)
          end if
       end do
    else
       do idims=1,ndim
          cmax(ixO^S)=dabs(w(ixO^S,iw_mom(idims)))+dsqrt(tmp(ixO^S))
          if(idims==1) then
            cmaxtot(ixO^S)=cmax(ixO^S)/block%ds(ixO^S,idims)
          else
            cmaxtot(ixO^S)=cmaxtot(ixO^S)+cmax(ixO^S)/block%ds(ixO^S,idims)
          end if
       end do
    end if
    ! courantmaxtots='max(summed c/dx)'
    courantmaxtots=maxval(cmaxtot(ixO^S))
    if(fld_debug)print *,'DT limit RADFORCE CSRAD=',courantpar/courantmaxtots
    if(courantmaxtots>smalldouble) dtnew=min(dtnew,courantpar/courantmaxtots)
    if(fld_debug)print *,'DT limit FINALLY ENFORCED IS NOW=',dtnew

  end subroutine fld_radforce_get_dt

  !> Sets the opacity in the w-array
  !> by calling mod_opal_opacity
  !> NOTE: assumes primitives in w
  !> NOTE: assuming opacity is local, not ok with cak line force
  subroutine fld_get_opacity_prim(w, x, ixI^L, ixO^L, fld_kappa, fl)
    use mod_global_parameters
    use mod_usr_methods
    use mod_opal_opacity
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out) :: fld_kappa(ixI^S)
    type(fld_fluid), intent(in)  :: fl

    integer :: ix^D
    double precision :: rho0,Temp0,kapp0
    double precision :: Temp(ixI^S)

    select case (trim(fld_opacity_law))
      case('const_norm')
        fld_kappa(ixO^S) = fld_kappa0
      case('const')
        fld_kappa(ixO^S) = fld_kappa0/unit_opacity
      case('opal')
        call fl%get_tgas(w,x,ixI^L,ixO^L,Temp)
        {do ix^D=ixOmin^D,ixOmax^D\ }
          rho0 = w(ix^D,iw_rho)*unit_density
          Temp0 = Temp(ix^D)*unit_temperature
          call set_opal_opacity(rho0,Temp0,kapp0)
          fld_kappa(ix^D) = kapp0/unit_opacity
        {enddo\ }
      case('special')
        if (.not. associated(usr_special_opacity)) then
          call mpistop("special opacity not defined")
        endif
        call usr_special_opacity(ixI^L, ixO^L, w, x, fld_kappa)
      case default
        call mpistop("Doesn't know opacity law")
      end select
  end subroutine fld_get_opacity_prim

  !> Returns Radiation Pressure as tensor
  !> NOTE: w is primitive on entry
  subroutine fld_get_radpress(w, x, ixI^L, ixO^L, rad_pressure, fl)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: rad_pressure(ixI^S,1:ndim,1:ndim)
    type(fld_fluid), intent(in)  :: fl
    integer :: i,j,nth
    double precision             :: eddington_tensor(ixI^S,1:ndim,1:ndim)
    double precision             :: lambda(ixI^S),fld_R(ixI^S)

    ! always use 4th order CD here
    nth=2
    call fld_get_eddington(w, x, ixI^L, ixO^L, eddington_tensor, lambda, fld_R, nth, fl)
    if(fld_debug.and..false.)then
       print *,'In get_radPress with nth=',nth,' on ixO=',ixO^L
       print *,'Max and Min value of fe'
       print *,maxval(eddington_tensor(ixO^S,1:ndim,1:ndim))
       print *,minval(eddington_tensor(ixO^S,1:ndim,1:ndim))
       print *,'Max and Min value of Erad'
       print *,maxval(w(ixO^S,iw_r_e))
       print *,minval(w(ixO^S,iw_r_e))
       print *,'End get_radPress'
    endif
    do i=1,ndim
      do j=1,ndim
        rad_pressure(ixO^S,i,j) = eddington_tensor(ixO^S,i,j)*w(ixO^S,iw_r_e)
      enddo
    enddo
  end subroutine fld_get_radpress

  !> This subroutine calculates flux limiter lambda according to fld_fluxlimiter
  !> It also calculates fld_R which is ratio of radiation scaleheight and mean free path
  !> NOTE: nth and ixI and ixO not free to choose here: TODO
  subroutine fld_get_fluxlimiter(w,x,ixI^L,ixO^L,fld_lambda,fld_R,nth,fl)
    use mod_global_parameters
    use mod_geometry
    use mod_usr_methods
    use mod_physics, only: phys_to_primitive
    integer, intent(in)           :: ixI^L,ixO^L,nth
    double precision, intent(in)  :: w(ixI^S,1:nw)
    double precision, intent(in)  :: x(ixI^S,1:ndim)
    double precision, intent(out) :: fld_R(ixI^S),fld_lambda(ixI^S)
    type(fld_fluid), intent(in)   :: fl

    double precision :: wprim(ixI^S,1:nw)

    wprim(ixI^S,1:nw)=w(ixI^S,1:nw)
    call phys_to_primitive(ixI^L,ixI^L,wprim,x)
    call fld_get_fluxlimiter_prim(wprim,x,ixI^L,ixO^L,fld_lambda,fld_R,nth,fl)

  end subroutine fld_get_fluxlimiter

  !> This subroutine calculates flux limiter lambda according to fld_fluxlimiter
  !> It also calculates fld_R which is ratio of radiation scaleheight and mean free path
  !> NOTE: this one operates on primitives
  !> NOTE: nth and ixI and ixO not free to choose 
  subroutine fld_get_fluxlimiter_prim(w,x,ixI^L,ixO^L,fld_lambda,fld_R,nth,fl)
    use mod_global_parameters
    use mod_geometry
    use mod_usr_methods
    integer, intent(in)           :: ixI^L,ixO^L,nth
    double precision, intent(in)  :: w(ixI^S,1:nw)
    double precision, intent(in)  :: x(ixI^S,1:ndim)
    double precision, intent(out) :: fld_R(ixI^S),fld_lambda(ixI^S)
    type(fld_fluid), intent(in)   :: fl

    integer :: idir, ix^D
    double precision :: kappa(ixI^S),normgrad2(ixI^S)
    double precision :: grad_r_e(ixI^S)

    select case(fld_fluxlimiter)
    case('Diffusion')
      ! optically thick limit
      fld_lambda(ixO^S) = 1.d0/3.d0
      fld_R(ixO^S) = zero
    case('FreeStream')
      ! optically thin limit
      normgrad2(ixO^S) = zero
      do idir=1,ndim
        call gradient(w(ixI^S,iw_r_e),ixI^L,ixO^L,idir,grad_r_e,nth)
        normgrad2(ixO^S) = normgrad2(ixO^S)+grad_r_e(ixO^S)**2
      end do
      call fld_get_opacity_prim(w,x,ixI^L,ixO^L,kappa,fl)
      ! Calculate R everywhere
      ! |grad E|/(rho kappa E)
      fld_R(ixO^S) = dsqrt(normgrad2(ixO^S))/(kappa(ixO^S)*w(ixO^S,iw_rho)*w(ixO^S,iw_r_e))
      where(normgrad2(ixO^S)<smalldouble**2)
        ! treat uniform case as diffusion limit
        fld_R(ixO^S)=zero
        fld_lambda(ixO^S) = 1.0d0/3.0d0
      elsewhere
        fld_lambda(ixO^S) = one/fld_R(ixO^S)
      endwhere
    case('Pomraning')
      ! Calculate R everywhere
      ! |grad E|/(rho kappa E)
      normgrad2(ixO^S) = zero
      do idir = 1,ndim
        call gradient(w(ixI^S,iw_r_e),ixI^L,ixO^L,idir,grad_r_e,nth)
        normgrad2(ixO^S) = normgrad2(ixO^S) + grad_r_e(ixO^S)**2
      end do
      call fld_get_opacity_prim(w,x,ixI^L,ixO^L,kappa,fl)
      fld_R(ixO^S) = dsqrt(normgrad2(ixO^S))/(kappa(ixO^S)*w(ixO^S,iw_rho)*w(ixO^S,iw_r_e))
      ! Calculate the flux limiter, lambda
      ! Levermore and Pomraning: lambda = (2 + R)/(6 + 3R + R^2)
      fld_lambda(ixO^S) = (2.d0+fld_R(ixO^S))/(6.d0+3*fld_R(ixO^S)+fld_R(ixO^S)**2)
    case('Minerbo')
      ! Calculate R everywhere
      ! |grad E|/(rho kappa E)
      normgrad2(ixO^S) = zero
      do idir = 1,ndim
        call gradient(w(ixI^S,iw_r_e),ixI^L,ixO^L,idir,grad_r_e,nth)
        normgrad2(ixO^S) = normgrad2(ixO^S) + grad_r_e(ixO^S)**2
      end do
      call fld_get_opacity_prim(w, x, ixI^L, ixO^L, kappa, fl)
      fld_R(ixO^S) = dsqrt(normgrad2(ixO^S))/(kappa(ixO^S)*w(ixO^S,iw_rho)*w(ixO^S,iw_r_e))
      ! Calculate the flux limiter, lambda
      ! Minerbo:
      {do ix^D = ixOmin^D,ixOmax^D\ }
        if(fld_R(ix^D) .lt. 3.d0/2.d0) then
          fld_lambda(ix^D) = 2.d0/(3.d0+dsqrt(9.d0+12.d0*fld_R(ix^D)**2))
        else
          fld_lambda(ix^D) = 1.d0/(1.d0+fld_R(ix^D)+dsqrt(1.d0+2.d0*fld_R(ix^D)))
        endif
      {enddo\}
    case('special')
      if (.not. associated(usr_special_fluxlimiter)) then
        call mpistop("special fluxlimiter not defined")
      endif
      call usr_special_fluxlimiter(ixI^L,ixO^L,w,x,fld_lambda,fld_R)
    case default
      call mpistop('Fluxlimiter unknown')
    end select

  end subroutine fld_get_fluxlimiter_prim

  !> Calculate Radiation Flux 
  !> NOTE: only for diagnostics purposes (w conservative on entry)
  !> This returns cell centered values for radiation flux 
  subroutine fld_get_radflux(w, x, ixI^L, ixO^L, rad_flux, fl)
    use mod_global_parameters
    use mod_geometry
    use mod_physics, only: phys_to_primitive
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out) :: rad_flux(ixI^S, 1:ndim)
    type(fld_fluid), intent(in)  :: fl

    integer :: idir,nth_for_fld
    double precision :: wprim(ixI^S,1:nw)
    double precision :: grad_r_e(ixI^S)
    double precision :: kappa(ixI^S), lambda(ixI^S), fld_R(ixI^S)

    wprim(ixI^S,1:nw)=w(ixI^S,1:nw)
    call phys_to_primitive(ixI^L,ixI^L,wprim,x)

    call fld_get_opacity_prim(wprim, x, ixI^L, ixO^L, kappa, fl)
    ! always use 4th order CD here
    nth_for_fld=2
    call fld_get_fluxlimiter_prim(wprim, x, ixI^L, ixO^L, lambda, fld_R, nth_for_fld, fl)
    !> Calculate the Flux using the fld closure relation
    !> F = -c*lambda/(kappa*rho) *grad E
    do idir = 1,ndim
      call gradient(wprim(ixI^S,iw_r_e),ixI^L,ixO^L,idir,grad_r_e,nth_for_fld)
      rad_flux(ixO^S,idir)=-(c_norm*lambda(ixO^S)/(kappa(ixO^S)*wprim(ixO^S,iw_rho)))*grad_r_e(ixO^S)
    end do
  end subroutine fld_get_radflux

  !> Calculate Eddington-tensor (where w is primitive)
  !> also feeds back the flux limiter lambda and R
  subroutine fld_get_eddington(w, x, ixI^L, ixO^L, eddington_tensor, lambda, fld_R, nth, fl)
    use mod_global_parameters
    use mod_geometry
    integer, intent(in)          :: ixI^L, ixO^L, nth
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out) :: eddington_tensor(ixI^S,1:ndim,1:ndim)
    double precision, intent(out) :: lambda(ixI^S),fld_R(ixI^S)
    type(fld_fluid), intent(in)  :: fl

    integer :: idir,jdir
    double precision :: normgrad2(ixI^S)
    double precision :: tmp(ixI^S),grad_r_e(ixI^S,1:ndim)
    double precision :: nn_regularized(ixI^S,1:ndim,1:ndim)

    normgrad2(ixO^S) = zero
    do idir = 1,ndim
      call gradient(w(ixI^S, iw_r_e),ixI^L,ixO^L,idir,tmp,nth)
      grad_r_e(ixO^S,idir)=tmp(ixO^S)
      normgrad2(ixO^S)=normgrad2(ixO^S)+tmp(ixO^S)**2
    end do
    do idir = 1,ndim
      do jdir = 1,ndim
       if(idir==jdir)then
          nn_regularized(ixO^S,idir,jdir)=(grad_r_e(ixO^S,idir)*grad_r_e(ixO^S,jdir)+smalldouble**2)/(normgrad2(ixO^S)+smalldouble**2)
       else
          nn_regularized(ixO^S,idir,jdir)=(grad_r_e(ixO^S,idir)*grad_r_e(ixO^S,jdir))/(normgrad2(ixO^S)+smalldouble**2)
       endif
      enddo
    enddo
    ! get lambda and R
    call fld_get_fluxlimiter_prim(w,x,ixI^L,ixO^L,lambda,fld_R,nth,fl)
    ! store f_e= lambda + lambda^2 R^2
    tmp(ixO^S) = lambda(ixO^S)+(lambda(ixO^S)*fld_R(ixO^S))**2
    do idir = 1,ndim
      ! first compute the isotropic (diagonal) part
      eddington_tensor(ixO^S,idir,idir) = half*(one-tmp(ixO^S))
    enddo
    do idir = 1,ndim
      do jdir = 1,ndim
        ! initialize off-diagonal part here
        if(idir .ne. jdir) eddington_tensor(ixO^S,idir,jdir) = zero
        ! add part depending on unit vectors along gradient E
        eddington_tensor(ixO^S,idir,jdir) = eddington_tensor(ixO^S,idir,jdir)+&
            half*(3.d0*tmp(ixO^S)-one)*nn_regularized(ixO^S,idir,jdir)
      enddo
    enddo
  end subroutine fld_get_eddington

  subroutine get_and_check_egas_Erad_from_conserved(w,ixI^L,ixO^L,e_gas,E_rad)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(out) :: e_gas(ixI^S),E_rad(ixI^S)
    integer :: ix^D

      !> e_gas is the INTERNAL ENERGY without KINETIC ENERGY
      e_gas(ixO^S) = w(ixO^S,iw_e)-half*sum(w(ixO^S,iw_mom(:))**2,dim=ndim+1)/w(ixO^S,iw_rho)
      if(allocated(iw_mag)) then
        e_gas(ixO^S) = e_gas(ixO^S)-half*sum(w(ixO^S,iw_mag(:))**2,dim=ndim+1)
      endif
      E_rad(ixO^S) = w(ixO^S,iw_r_e)

      if(check_small_values.and..not.fix_small_values)then
         {do ix^DB= ixOmin^DB,ixOmax^DB\}
           if(e_gas(ix^D)<small_e.or.E_rad(ix^D)<small_r_e) then
              write(*,*) "Error in FLD get_egas_Erad: small value"
              write(*,*) "Iteration: ", it, " Time: ", global_time
              write(*,*) "Cell number: ", ix^D
              write(*,*) "internal energy density  is=",e_gas(ix^D)," versus   small_e=",small_e
              write(*,*) "radiation energy density is=",E_rad(ix^D)," versus small_r_e=",small_r_e
              call mpistop("FLD error:May need to turn on fixes")
           end if
         {end do\}
      endif

      if(fix_small_values)then
        {do ix^D = ixOmin^D,ixOmax^D\ }
          e_gas(ix^D) = max(e_gas(ix^D),small_e)
          E_rad(ix^D) = max(E_rad(ix^D),small_r_e)
        {enddo\}
      endif

  end subroutine get_and_check_egas_Erad_from_conserved

  !> Calling all subroutines to perform the multigrid method
  !> Communicates rad_e and diff_coeff to multigrid library
  !> Advance psa=psb+dtfactor*qdt*F_im(psa)
  subroutine fld_implicit_update(dtfactor,qdt,qtC,psa,psb,fl)
    use mod_global_parameters
    use mod_forest
    use mod_multigrid_coupling
    use mod_input_output, only: get_global_maxima,get_global_minima

    type(state), target :: psa(max_blocks)
    type(state), target :: psb(max_blocks)
    double precision, intent(in) :: qdt
    double precision, intent(in) :: qtC
    double precision, intent(in) :: dtfactor
    type(fld_fluid), intent(in)  :: fl

    integer, parameter           :: max_its = 100
    integer                      :: n,ixO^L,ix^D
    double precision             :: res, max_residual, mg_lambda, fac
    double precision :: wmax(nw),wmin(nw)
    double precision :: wmaxb(nw),wminb(nw)
    integer :: iigrid, igrid

    if(fld_no_mg)return
    
    ! we need first to compute the (variable) diffusion coefficient on entire grid
    !   this must be done in mesh+1 ghostcell layer
    call update_diffcoeff(psa,fl)

    ! now we multiply the diffusion coefficient with dtfactor*dt on entire mesh+1 domain
    ixO^L=ixM^LL^LADD1;
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       {do ix^D = ixOmin^D,ixOmax^D\ }
          psa(igrid)%w(ix^D,i_diff_mg)= psa(igrid)%w(ix^D,i_diff_mg)*dtfactor*qdt
       {enddo\}
    end do

    fac = 1.0d0
    max_residual = fld_diff_tol

    if(fld_debug)then
       call get_global_maxima(wmax,psa)
       call get_global_minima(wmin,psa)
       call get_global_maxima(wmaxb,psb)
       call get_global_minima(wminb,psb)
       if(mype==0)then
          ! the MG needs to be scaled such that everything is order unity
          print *,'Currently at time=',global_time,' time step=',qdt,' dtfactor=',dtfactor
          print *,'at start of MG solver, we have fld_diff_tol       =',fld_diff_tol
          print *,'at start of MG solver, we have LHS E_rad range as :',wmax(iw_r_e),wmin(iw_r_e)
          print *,'at start of MG solver, we have Diff coeff range as:',wmax(i_diff_mg),wmin(i_diff_mg)
          print *,'at start of MG solver, we have density range as   :',wmax(iw_rho),wmin(iw_rho)
          print *,'at start of MG solver, we have RHS E_rad range as :',wmaxb(iw_r_e),wminb(iw_r_e)
          print *,'at start of MG solver, we have qdt as',qdt,' and max_residual=',max_residual
          print *,'at start of MG solver, ratio coeffs on level 1  =',wmax(i_diff_mg)/(dx(1,1)**2)
          print *,'at start of MG solver, ratio coeffs on level max=',wmax(i_diff_mg)/(dx(1,refine_max_level)**2)
       endif
    endif

    call mg_set_methods(mg)
    if(.not. mg%is_allocated) call mpistop("multigrid tree not allocated yet")

    ! Here we handle the global helmholtz problem with variable coefficient
    ! The equation we solve is div([D]^n nabla Erad^(n+1)) -(1/dt)Erad^(n+1)=-(1/dt)Erad^n
    ! we reformulate to div([D x dt]^n nabla Erad^(n+1)) -Erad^(n+1)=-Erad^n
    ! Helmholtz equation is div(eps nabla phi) -mg_lambda phi = f
    !    hence eps is our variable coefficient dt x D=fld_lambda*c/(kappa*rho)
    !    hence phi is Erad and mg_lambda is unity
    call vhelmholtz_set_lambda(fac)
    ! copy in the (variable) diffusion coefficient in mg_iveps
    ! NOTE: this copies also the ghostcell values for this coefficient to mg
    call mg_copy_to_tree(i_diff_mg, mg_iveps, factor=fac, state_from=psa)
    ! copy in the Erad variable in mg_iphi (the one we solve for)
    call mg_copy_to_tree(iw_r_e, mg_iphi, factor=fac, state_from=psa)
    ! copy in RHS f factor as Erad with factor -1
    call mg_copy_to_tree(iw_r_e, mg_irhs, factor=-fac, state_from=psb)
    ! becuase the variable coefficient is needed on lower grid levels in mg
    ! we need to restrict and adopt BCs for this variable: Neumann is set 
    call mg_restrict(mg, mg_iveps)
    call mg_fill_ghost_cells(mg, mg_iveps)
    ! Now try solving with MG
    call mg_fas_fmg(mg, .true., max_res=res)
    if(fld_debug.and.mype==0)print *,'MG residual obtained with FMG is =',res
    do n = 1, max_its
      if(res < max_residual) exit
      call mg_fas_vcycle(mg, max_res=res)
      if(fld_debug.and.mype==0)print *,'MG residual obtained with Vcycle =',res,' at step =',n
    end do
    if(fld_debug.and.mype==0)print *,'FINAL MG residual obtained is =',res
    if(res .ge. max_residual) then
      if (mype == 0) then 
        write(*,*) it, ' residual from MG ', res
        write(*,*) it, ' max_residual in MG ', max_residual
        write(*,*) it, ' dtfactor*qdt in MG ', qdt*dtfactor
        print *,'Currently at time=',global_time,' time step=',qdt,' dtfactor=',dtfactor
      endif
      if(fld_force_MG_converged) then
          call mpistop("no convergence in MG")
      else
          if(mype==0) write(*,*) 'WARNING for it=',it,' NO CONVERGENCE IN MG but still carry on'
      endif
    end if
    ! copy back the Erad variable in iw_r_e
    call mg_copy_from_tree_gc(mg_iphi, iw_r_e, state_to=psa)

   if(check_small_values.and..not.fix_small_values)then
       ixO^L=ixM^LL;
       do iigrid=1,igridstail; igrid=igrids(iigrid);
         {do ix^DB= ixOmin^DB,ixOmax^DB\}
           if(psa(igrid)%w(ix^D,iw_r_e)<small_r_e) then
              write(*,*) "Error in FLD fld_implicit_update: small value"
              write(*,*) "of radiation energy density after MG"
              write(*,*) "Iteration: ", it, " Time: ", global_time
              write(*,*) "Location: ", psa(igrid)%x(ix^D,:)," on grid",igrid
              write(*,*) "Cell number: ", ix^D
              write(*,*) "radiation energy density is=",psa(igrid)%w(ix^D,iw_r_e)," versus small_r_e=",small_r_e
              call mpistop("FLD error:May need to turn on fixes")
           end if
         {end do\}
       end do
    endif

    if(fix_small_values)then
       ixO^L=ixM^LL;
       do iigrid=1,igridstail; igrid=igrids(iigrid);
         {do ix^D = ixOmin^D,ixOmax^D\ }
           psa(igrid)%w(ix^D,iw_r_e)= max(psa(igrid)%w(ix^D,iw_r_e),small_r_e)
         {enddo\}
       end do
    endif

  end subroutine fld_implicit_update

  subroutine update_diffcoeff(psa,fl)
    use mod_global_parameters
    type(state), target :: psa(max_blocks)
    type(fld_fluid), intent(in) :: fl

    integer :: iigrid, igrid, ixO^L

    ! we will need diffusion coefficients in i-1 i i+1
    ixO^L=ixM^LL^LADD1;
    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
        call fld_get_diffcoef_central(psa(igrid)%w, psa(igrid)%x, ixG^LL, ixO^L, fl)
    end do
    !$OMP END PARALLEL DO

  end subroutine update_diffcoeff

  !> inplace update of psa==>F_im(psa)
  subroutine fld_evaluate_implicit(qtC,psa,fl)
    use mod_global_parameters
    type(state), target :: psa(max_blocks)
    double precision, intent(in) :: qtC
    type(fld_fluid), intent(in)  :: fl
    integer :: iigrid, igrid
    integer :: ixO^L

    if(fld_no_mg)return

    call update_diffcoeff(psa,fl)

    ixO^L=ixM^LL;
    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
       call evaluate_diffterm_onegrid(ixG^LL,ixO^L,psa(igrid)%w,psa(igrid)%x)
    end do
    !$OMP END PARALLEL DO

  end subroutine fld_evaluate_implicit

  !> inplace update of psa==>F_im(psa)
  subroutine evaluate_diffterm_onegrid(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)

    double precision :: divF(ixI^S)
    integer :: idir, jxO^L, hxO^L

    if(.not.slab)call mpistop("laplacian coded up for uniform cartesian grid")

    ! Here we use diffusion coefficients in positions i-1 i i+1 and exploit harmonic means i.e. 2ab/(a+b)
    ! since this is how the multigrid library handles the div(eps nabla phi) term in Helmholtz equation
    !   div(eps nabla phi) - mg_lambda phi = f
    divF(ixO^S) = 0.d0
    do idir = 1,ndim
       hxO^L=ixO^L-kr(idir,^D);
       jxO^L=ixO^L+kr(idir,^D);
       divF(ixO^S) = divF(ixO^S) + &
           (w(jxO^S,iw_r_e)*two*w(ixO^S,i_diff_mg)*w(jxO^S,i_diff_mg)/(w(ixO^S,i_diff_mg) + w(jxO^S,i_diff_mg)) &
          -w(ixO^S,iw_r_e)*(two*w(ixO^S,i_diff_mg)*w(jxO^S,i_diff_mg)/(w(ixO^S,i_diff_mg) + w(jxO^S,i_diff_mg)) &
                           +two*w(ixO^S,i_diff_mg)*w(hxO^S,i_diff_mg)/(w(ixO^S,i_diff_mg) + w(hxO^S,i_diff_mg))) &
           +w(hxO^S,iw_r_e)*two*w(ixO^S,i_diff_mg)*w(hxO^S,i_diff_mg)/(w(ixO^S,i_diff_mg) + w(hxO^S,i_diff_mg)))/dxlevel(idir)**2
       ! below uses artihmetic mean, different from mg method
       !divF(ixO^S) = divF(ixO^S) + &
       !    (w(jxO^S,iw_r_e)*half*(w(ixO^S,i_diff_mg) + w(jxO^S,i_diff_mg)) &
       !   -w(ixO^S,iw_r_e)*(half*(w(ixO^S,i_diff_mg) + w(jxO^S,i_diff_mg)) &
       !                    +half*(w(ixO^S,i_diff_mg) + w(hxO^S,i_diff_mg))) &
       !    +w(hxO^S,iw_r_e)*half*(w(ixO^S,i_diff_mg) + w(hxO^S,i_diff_mg)))/dxlevel(idir)**2
    enddo
    ! only the E variable is handled implicitly, all else must be zero here
    w(ixO^S,1:nw)=zero
    w(ixO^S,iw_r_e) = divF(ixO^S)

  end subroutine evaluate_diffterm_onegrid

  !> Calculates cell-centered diffusion coefficient to be used in multigrid
  subroutine fld_get_diffcoef_central(w, x, ixI^L, ixO^L, fl)
    use mod_global_parameters
    use mod_geometry
    use mod_usr_methods
    use mod_physics, only: phys_to_primitive
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    type(fld_fluid), intent(in)  :: fl

    double precision :: wprim(ixI^S,1:nw)
    double precision :: kappa(ixI^S),lambda(ixI^S),fld_R(ixI^S)

    wprim(ixI^S,1:nw)=w(ixI^S,1:nw)
    ! ensure entries in entire ixI range
    call phys_to_primitive(ixI^L,ixI^L,wprim,x)
    call fld_get_opacity_prim(wprim, x, ixI^L, ixO^L, kappa, fl)
    ! note that we use central difference here (last argument is 1 or 2)
    call fld_get_fluxlimiter_prim(wprim, x, ixI^L, ixO^L, lambda, fld_R, nth_for_diff_mg, fl)
    w(ixO^S,i_diff_mg) = c_norm*lambda(ixO^S)/(kappa(ixO^S)*wprim(ixO^S,iw_rho))
    if(associated(usr_special_diffcoef)) call usr_special_diffcoef(w, wprim, x, ixI^L, ixO^L)
    if(minval(w(ixO^S,i_diff_mg))<smalldouble) then
        print *,'min diffcoef=',minval(w(ixO^S,i_diff_mg))
        call mpistop("too small diffusion coefficient")
    endif
    if(maxval(w(ixO^S,i_diff_mg))>bigdouble)  call mpistop("too large diffusion coefficient")

    if(fld_debug.and..false.)then
      print *,'setting diffcoefs with data on',ixI^L
      print *,'min diffcoef=',minval(w(ixO^S,i_diff_mg))
      print *,'min lambda kappa rho fld_R'
      print *,minval(lambda(ixO^S))
      print *,minval(kappa(ixO^S))
      print *,minval(wprim(ixO^S,iw_rho))
      print *,minval(fld_R(ixO^S))
      print *,'max diffcoef=',maxval(w(ixO^S,i_diff_mg))
      print *,'max lambda kappa rho fld_R'
      print *,maxval(lambda(ixO^S))
      print *,maxval(kappa(ixO^S))
      print *,maxval(wprim(ixO^S,iw_rho))
      print *,maxval(fld_R(ixO^S))
      print *,'done setting diffcoefs in slot',i_diff_mg,' on range',ixO^L
    endif

  end subroutine fld_get_diffcoef_central


  !> Find the root of the 4th degree polynomial using the bisection method
  subroutine Bisection_method(e_gas, c0, c1)
    use mod_global_parameters
    double precision, intent(in)    :: c0, c1
    double precision, intent(inout) :: e_gas
    double precision :: bisect_a, bisect_b, bisect_c
    integer :: n, max_its

    n = 0
    max_its = 100
    bisect_a = zero
    bisect_b = min(dabs(c0/c1),dabs(c0)**(1.d0/4.d0))+smalldouble
    do while(dabs(bisect_b-bisect_a) .ge. fld_bisect_tol)
      bisect_c = (bisect_a + bisect_b)/two
      n = n +1
      if(n .gt. max_its) then
        call mpistop('No convergece in bisection scheme')
      endif
      if(Polynomial_Bisection(bisect_a, c0, c1)*&
      Polynomial_Bisection(bisect_b, c0, c1) .lt. zero) then
        if(Polynomial_Bisection(bisect_a, c0, c1)*&
        Polynomial_Bisection(bisect_c, c0, c1) .lt. zero) then
          bisect_b = bisect_c
        elseif(Polynomial_Bisection(bisect_b, c0, c1)*&
        Polynomial_Bisection(bisect_c, c0, c1) .lt. zero) then
          bisect_a = bisect_c
        elseif(Polynomial_Bisection(bisect_a, c0, c1) .eq. zero) then
          bisect_b = bisect_a
          bisect_c = bisect_a
          goto 2435
        elseif(Polynomial_Bisection(bisect_b, c0, c1) .eq. zero) then
          bisect_a = bisect_b
          bisect_c = bisect_b
          goto 2435
        elseif(Polynomial_Bisection(bisect_c, c0, c1) .eq. zero) then
          bisect_a = bisect_c
          bisect_b = bisect_c
          goto 2435
        else
          call mpistop("Problem with fld bisection method")
        endif
      elseif(Polynomial_Bisection(bisect_a, c0, c1) &
        - Polynomial_Bisection(bisect_b, c0, c1) .lt. fld_bisect_tol*Polynomial_Bisection(bisect_a, c0, c1)) then
        goto 2435
      else
        bisect_a = e_gas
        bisect_b = e_gas
        if(fld_debug)print*, "IGNORING GAS-RAD ENERGY EXCHANGE ", c0, c1
        if(fld_debug)print*, Polynomial_Bisection(bisect_a, c0, c1), Polynomial_Bisection(bisect_b, c0, c1)
        call mpistop('issues in bisection scheme')
        if(Polynomial_Bisection(bisect_a, c0, c1) .le. smalldouble) then
          bisect_b = bisect_a
        elseif(Polynomial_Bisection(bisect_a, c0, c1) .le. smalldouble) then
          bisect_a = bisect_b
        endif
        goto 2435
      endif
    enddo
      2435 e_gas = (bisect_a + bisect_b)/two
  end subroutine Bisection_method

  !> Find the root of the 4th degree polynomial using the Newton method
  subroutine Newton_method(e_gas, c0, c1)
    use mod_global_parameters
    double precision, intent(in)    :: c0, c1
    double precision, intent(inout) :: e_gas
    double precision :: xval, yval, der, deltax
    integer :: ii

    yval = bigdouble
    xval = e_gas
    der = one
    deltax = one
    ii = 0
    !> Compare error with dx = dx/dy dy
    do while(dabs(deltax) .gt. fld_bisect_tol)
      yval = Polynomial_Bisection(xval, c0, c1)
      der  = dPolynomial_Bisection(xval, c0, c1)
      deltax = -yval/der
      xval = xval + deltax
      ii = ii + 1
      if(ii .gt. 1d3) then
        if(fld_debug)print*, 'skip to bisection algorithm'
        call Bisection_method(e_gas, c0, c1)
        return
      endif
    enddo
    e_gas = xval
  end subroutine Newton_method

  !> Find the root of the 4th degree polynomial using the Halley method
  subroutine Halley_method(e_gas, c0, c1)
    use mod_global_parameters
    double precision, intent(in)    :: c0, c1
    double precision, intent(inout) :: e_gas
    double precision :: xval, yval, der, dder, deltax
    integer :: ii

    yval = bigdouble
    xval = e_gas
    der = one
    dder = one
    deltax = one
    ii = 0
    !> Compare error with dx = dx/dy dy
    do while (dabs(deltax) .gt. fld_bisect_tol)
      yval = Polynomial_Bisection(xval, c0, c1)
      der  = dPolynomial_Bisection(xval, c0, c1)
      dder = ddPolynomial_Bisection(xval, c0, c1)
      deltax = -two*yval*der/(two*der**2 - yval*dder)
      xval = xval + deltax
      ii = ii + 1
      if(ii .gt. 1d3) then
        if(fld_debug)print*, 'skip to Newton algorithm'
        call Newton_method(e_gas, c0, c1)
        return
      endif
    enddo
    e_gas = xval
  end subroutine Halley_method

  !> Evaluate polynomial at argument e_gas
  function Polynomial_Bisection(e_gas, c0, c1) result(val)
    use mod_global_parameters
    double precision, intent(in) :: e_gas
    double precision, intent(in) :: c0, c1
    double precision :: val

    val = e_gas**4.d0 + c1*e_gas - c0
  end function Polynomial_Bisection

  !> Evaluate first derivative of polynomial at argument e_gas
  function dPolynomial_Bisection(e_gas, c0, c1) result(der)
    use mod_global_parameters
    double precision, intent(in) :: e_gas
    double precision, intent(in) :: c0, c1
    double precision :: der

    der = 4.d0*e_gas**3.d0 + c1
  end function dPolynomial_Bisection

  !> Evaluate second derivative of polynomial at argument e_gas
  function ddPolynomial_Bisection(e_gas, c0, c1) result(dder)
    use mod_global_parameters
    double precision, intent(in) :: e_gas
    double precision, intent(in) :: c0, c1
    double precision :: dder

    dder = 4.d0*3.d0*e_gas**2.d0
  end function ddPolynomial_Bisection
end module mod_fld
