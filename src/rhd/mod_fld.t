!> Nicolas Moens
!> Module for including flux limited diffusion (FLD)-approximation in Radiation-hydrodynamics simulations using mod_rhd
!> Based on Turner and stone 2001. See
!> [1]Moens, N., Sundqvist, J. O., El Mellah, I., Poniatowski, L., Teunissen, J., and Keppens, R.,
!> “Radiation-hydrodynamics with MPI-AMRVAC . Flux-limited diffusion”,
!> <i>Astronomy and Astrophysics</i>, vol. 657, 2022. doi:10.1051/0004-6361/202141023.
!> For more information.

module mod_fld
    implicit none

    !> source split for energy interact and radforce:
    logical :: fld_Eint_split = .false.
    logical :: fld_Radforce_split = .false.

    !> Opacity per unit of unit_density
    double precision, public :: fld_kappa0 = 0.34d0

    !> mean particle mass
    double precision, public :: fld_mu = 0.6d0

    !> Tolerance for bisection method for Energy sourceterms
    !> This is a percentage of the minimum of gas- and radiation energy
    double precision, public :: fld_bisect_tol = 1.d-4

    !> Tolerance for adi method for radiative Energy diffusion
    double precision, public :: fld_diff_tol = 1.d-4

    !> Number for splitting the diffusion module
    double precision, public :: diff_crit

    !> Use constant Opacity?
    character(len=8)  :: fld_opacity_law = 'const'
    character(len=50) :: fld_opal_table = 'Y09800' !>'xxxxxx'

    !> Diffusion limit lambda = 0.33
    character(len=16) :: fld_fluxlimiter = 'Pomraning'

    !> diffusion coefficient for multigrid method
    integer :: i_diff_mg

    !> Which method to solve diffusion part
    character(len=8) :: fld_diff_scheme = 'mg'

    !> Which method to find the root for the energy interaction polynomial
    character(len=8) :: fld_interaction_method = 'Halley'

    !> Take running average for Diffusion coefficient
    logical :: diff_coef_filter = .false.
    integer :: size_D_filter = 1

    !> Take a running average over the fluxlimiter
    logical :: flux_lim_filter = .false.
    integer :: size_L_filter = 1

    !> Use or dont use lineforce opacities
    logical :: Lineforce_opacities = .false.

    !> Resume run when multigrid returns error
    logical :: diffcrash_resume = .true.

    !> Index for Flux weighted opacities
    integer, allocatable, public :: i_opf(:)

    !> A copy of rhd_Gamma
    double precision, private, protected :: fld_gamma

    !> running timestep for diffusion solver, initialised as zero
    double precision :: dt_diff = 0.d0

    !> public methods
    !> these are called in mod_rhd_phys
    public :: get_fld_rad_force
    public :: get_fld_energy_interact
    public :: fld_radforce_get_dt
    public :: fld_init
    public :: fld_get_radflux
    public :: fld_get_radpress
    public :: fld_get_fluxlimiter
    public :: fld_get_opacity

  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!! GENERAL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Reading in fld-list parameters from .par file
  subroutine fld_params_read(files)
    use mod_global_parameters, only: unitpar
    use mod_constants
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /fld_list/ fld_kappa0, fld_Eint_split, fld_Radforce_split, &
    fld_bisect_tol, fld_diff_tol,&
    fld_opacity_law, fld_fluxlimiter, fld_diff_scheme, fld_interaction_method, &
    diff_coef_filter, size_D_filter, flux_lim_filter, size_L_filter, &
    lineforce_opacities, diffcrash_resume, fld_opal_table

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, fld_list, end=111)
       111    close(unitpar)
    end do

  end subroutine fld_params_read

  !> Initialising FLD-module:
  !> Read opacities
  !> Initialise Multigrid
  !> adimensionalise kappa
  !> Add extra variables to w-array, flux, kappa, eddington Tensor
  !> Lambda and R
  !> ...
  subroutine fld_init(He_abundance, rhd_radiation_diffusion, rhd_gamma)
    use mod_global_parameters
    use mod_variables
    use mod_physics
    use mod_opal_opacity, only: init_opal_table
    use mod_multigrid_coupling

    double precision, intent(in) :: He_abundance, rhd_gamma
    logical, intent(in) :: rhd_radiation_diffusion
    double precision :: sigma_thomson
    integer :: idir,jdir

    character(len=1) :: ind_1
    character(len=1) :: ind_2
    character(len=2) :: cmp_f
    character(len=5) :: cmp_e

    !> read par files
    call fld_params_read(par_files)

    !> Set lineforce opacities as variable
    if (lineforce_opacities) then
      allocate(i_opf(ndim))
      do idir = 1,ndim
        write(ind_1,'(I1)') idir
        cmp_f = 'k' // ind_1
        i_opf(idir) = var_set_extravar(cmp_f,cmp_f)
      enddo
    endif

    if (rhd_radiation_diffusion) then
      if (fld_diff_scheme .eq. 'mg') then

        use_multigrid = .true.

        phys_implicit_update => Diffuse_E_rad_mg
        phys_evaluate_implicit => Evaluate_E_rad_mg

        mg%n_extra_vars = 1
        mg%operator_type = mg_vhelmholtz

      endif
    endif

    i_diff_mg = var_set_extravar("D", "D")

    !> Need mean molecular weight
    fld_mu = (1.+4*He_abundance)/(2.+3.*He_abundance)

    !> set rhd_gamma
    fld_gamma = rhd_gamma

    !> Read in opacity table if necesary
    if (fld_opacity_law .eq. 'opal') call init_opal_table(He_abundance,fld_opal_table)
    if ((fld_opacity_law .eq. 'thomson') .or. (fld_opacity_law .eq. 'fastwind'))  then
      sigma_thomson = 6.6524585d-25
      fld_kappa0 = sigma_thomson/const_mp * (1.+2.*He_abundance)/(1.+4.*He_abundance)
    endif
  end subroutine fld_init

  !> w[iw]=w[iw]+qdt*S[wCT,qtC,x] where S is the source based on wCT within ixO
  !> This subroutine handles the radiation force
  subroutine get_fld_rad_force(qdt,ixI^L,ixO^L,wCT,w,x,&
       energy,qsourcesplit,active)
    use mod_constants
    use mod_global_parameters
    use mod_usr_methods
    use mod_geometry

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, x(ixI^S,1:ndim)
    double precision, intent(in)    :: wCT(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    logical, intent(in) :: energy,qsourcesplit
    logical, intent(inout) :: active

    double precision :: radiation_forceCT(ixO^S,1:ndim)
    double precision :: kappaCT(ixO^S)
    double precision :: rad_fluxCT(ixO^S,1:ndim)

    double precision :: div_v(ixI^S,1:ndim,1:ndim)
    double precision :: edd(ixO^S,1:ndim,1:ndim)
    double precision :: nabla_vP(ixO^S)
    double precision :: vel(ixI^S), grad_v(ixI^S), grad0_v(ixO^S)
    double precision :: grad_E(ixO^S)

    integer :: idir, jdir

    double precision :: fld_R(ixO^S), lambda(ixO^S)

    !> Calculate and add sourceterms
    if(qsourcesplit .eqv. fld_Radforce_split) then
      active = .true.

      call fld_get_opacity(wCT, x, ixI^L, ixO^L, kappaCT)
      call fld_get_radflux(wCT, x, ixI^L, ixO^L, rad_fluxCT)
      call fld_get_fluxlimiter(wCT, x, ixI^L, ixO^L, lambda, fld_R)

      do idir = 1,ndim
        !> Radiation force = kappa*rho/c *Flux = lambda gradE
        radiation_forceCT(ixO^S,idir) = kappaCT(ixO^S)*rad_fluxCT(ixO^S,idir)/(const_c/unit_velocity)

        ! call gradientO(wCT(ixI^S,iw_r_e),x,ixI^L,ixO^L,idir,grad_E,nghostcells)
        ! radiation_forceCT(ixO^S,idir) = lambda(ixO^S)*grad_E(ixO^S)

        !> Momentum equation source term
        w(ixO^S,iw_mom(idir)) = w(ixO^S,iw_mom(idir)) &
            + qdt * wCT(ixO^S,iw_rho)*radiation_forceCT(ixO^S,idir)
        ! w(ixO^S,iw_mom(idir)) = w(ixO^S,iw_mom(idir)) &
        !     + qdt *radiation_forceCT(ixO^S,idir)

        ! if (energy .and. .not. block%e_is_internal) then
          !> Energy equation source term (kinetic energy)
          w(ixO^S,iw_e) = w(ixO^S,iw_e) &
              + qdt * wCT(ixO^S,iw_mom(idir))*radiation_forceCT(ixO^S,idir)
          ! w(ixO^S,iw_e) = w(ixO^S,iw_e) &
          !     + qdt * wCT(ixO^S,iw_mom(idir))/wCT(ixO^S,iw_rho)*radiation_forceCT(ixO^S,idir)
        ! endif
      enddo

      !> Photon tiring
      !> calculate tensor div_v
      !> !$OMP PARALLEL DO
      do idir = 1,ndim
        do jdir = 1,ndim
          vel(ixI^S) = wCt(ixI^S,iw_mom(jdir))/wCt(ixI^S,iw_rho)

          call gradient(vel,ixI^L,ixO^L,idir,grad_v)
          div_v(ixO^S,idir,jdir) = grad_v(ixO^S)

          ! call gradientO(vel,x,ixI^L,ixO^L,idir,grad0_v,nghostcells)
          ! div_v(ixO^S,idir,jdir) = grad0_v(ixO^S)
        enddo
      enddo
      !> !$OMP END PARALLEL DO

      call fld_get_eddington(wCt, x, ixI^L, ixO^L, edd)

      !> VARIABLE NAMES DIV ARE ACTUALLY GRADIENTS
      {^IFONED
      nabla_vP(ixO^S) = div_v(ixO^S,1,1)*edd(ixO^S,1,1)  &
      }

      {^IFTWOD
      !>eq 34 Turner and stone (Only 2D)
      nabla_vP(ixO^S) = div_v(ixO^S,1,1)*edd(ixO^S,1,1)  &
                   + div_v(ixO^S,1,2)*edd(ixO^S,1,2)  &
                   + div_v(ixO^S,2,1)*edd(ixO^S,2,1)  &
                   + div_v(ixO^S,2,2)*edd(ixO^S,2,2)
      }

      {^IFTHREED
      nabla_vP(ixO^S) = div_v(ixO^S,1,1)*edd(ixO^S,1,1)  &
                   + div_v(ixO^S,1,2)*edd(ixO^S,1,2)  &
                   + div_v(ixO^S,1,3)*edd(ixO^S,1,3)  &
                   + div_v(ixO^S,2,1)*edd(ixO^S,2,1)  &
                   + div_v(ixO^S,2,2)*edd(ixO^S,2,2)  &
                   + div_v(ixO^S,2,3)*edd(ixO^S,2,3)  &
                   + div_v(ixO^S,3,1)*edd(ixO^S,3,1)  &
                   + div_v(ixO^S,3,2)*edd(ixO^S,3,2)  &
                   + div_v(ixO^S,3,3)*edd(ixO^S,3,3)
      }

      w(ixO^S,iw_r_e) = w(ixO^S,iw_r_e) &
          - qdt * nabla_vP(ixO^S)*wCt(ixO^S,iw_r_e)

    end if

  end subroutine get_fld_rad_force


  subroutine fld_radforce_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
    use mod_global_parameters
    use mod_usr_methods

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: dx^D, x(ixI^S,1:ndim), w(ixI^S,1:nw)
    double precision, intent(inout) :: dtnew

    double precision :: radiation_force(ixO^S,1:ndim)
    double precision :: rad_flux(ixO^S,1:ndim)
    double precision :: kappa(ixO^S)
    double precision                :: dxinv(1:ndim), max_grad
    integer                         :: idim

    ^D&dxinv(^D)=one/dx^D;

    call fld_get_opacity(w, x, ixI^L, ixO^L, kappa)
    call fld_get_radflux(w, x, ixI^L, ixO^L, rad_flux)

    do idim = 1, ndim
      radiation_force(ixO^S,idim) = w(ixO^S,iw_rho)*kappa(ixO^S)*rad_flux(ixO^S,idim)/(const_c/unit_velocity)
      max_grad = maxval(abs(radiation_force(ixO^S,idim)))
      max_grad = max(max_grad, epsilon(1.0d0))
      dtnew = min(dtnew, courantpar / sqrt(max_grad * dxinv(idim)))
    end do

  end subroutine fld_radforce_get_dt

  !> w[iw]=w[iw]+qdt*S[wCT,qtC,x] where S is the source based on wCT within ixO
  !> This subroutine handles the energy exchange between gas and radiation
  subroutine get_fld_energy_interact(qdt,ixI^L,ixO^L,wCT,w,x,&
       energy,qsourcesplit,active)
    use mod_constants
    use mod_global_parameters
    use mod_usr_methods

    use mod_physics, only: phys_get_pthermal  !needed to get temp

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, x(ixI^S,1:ndim)
    double precision, intent(in)    :: wCT(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    logical, intent(in) :: energy,qsourcesplit
    logical, intent(inout) :: active

    !> Calculate and add sourceterms
    if(qsourcesplit .eqv. fld_Eint_split) then
      active = .true.
      !> Add energy sourceterms
      call Energy_interaction(w, w, x, ixI^L, ixO^L)
    end if
  end subroutine get_fld_energy_interact

  !> Sets the opacity in the w-array
  !> by calling mod_opal_opacity
  subroutine fld_get_opacity(w, x, ixI^L, ixO^L, fld_kappa)
    use mod_global_parameters
    use mod_physics, only: phys_get_pthermal
    use mod_physics, only: phys_get_tgas
    use mod_usr_methods
    use mod_opal_opacity

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out) :: fld_kappa(ixO^S)
    double precision :: Temp(ixI^S), pth(ixI^S), a2(ixO^S)
    double precision :: rho0,Temp0,n,sigma_b
    double precision :: akram, bkram
    double precision :: vth(ixO^S), gradv(ixI^S), eta(ixO^S), t(ixO^S)

    integer :: i,j,ix^D, idir

    select case (fld_opacity_law)
      case('const')
        fld_kappa(ixO^S) = fld_kappa0/unit_opacity
      case('thomson')
        fld_kappa(ixO^S) = fld_kappa0/unit_opacity
      case('kramers')
        rho0 = half !> Take lower value of rho in domain
        fld_kappa(ixO^S) = fld_kappa0/unit_opacity*((w(ixO^S,iw_rho)/rho0))
      case('bump')
        !> Opacity bump
        rho0 = 0.2d0 !0.5d-1
        n = 7.d0
        sigma_b = 2.d-2
        !fld_kappa(ixO^S) = fld_kappa0/unit_opacity*(one + n*dexp(-((rho0  - w(ixO^S,iw_rho))**two)/rho0))
        fld_kappa(ixO^S) = fld_kappa0/unit_opacity*(one + n*dexp(-one/sigma_b*(dlog(w(ixO^S,iw_rho)/rho0))**two))
      case('non_iso')
        call phys_get_pthermal(w,x,ixI^L,ixO^L,Temp)
        Temp(ixO^S)=Temp(ixO^S)/w(ixO^S,iw_rho)

        rho0 = 0.5d0 !> Take lower value of rho in domain
        Temp0 = one
        n = -7.d0/two
        fld_kappa(ixO^S) = fld_kappa0/unit_opacity*(w(ixO^S,iw_rho)/rho0)*(Temp(ixO^S)/Temp0)**n
      case('fastwind')
        call phys_get_pthermal(w,x,ixI^L,ixO^L,pth)
        a2(ixO^S) = pth(ixO^S)/w(ixO^S,iw_rho)*unit_velocity**2.d0

        akram = 13.1351597305
        bkram = -4.5182188206

        fld_kappa(ixO^S) = fld_kappa0/unit_opacity &
        * (1.d0+10.d0**akram*w(ixO^S,iw_rho)*unit_density*(a2(ixO^S)/1.d12)**bkram)

        {do ix^D=ixOmin^D,ixOmax^D\ }
          !> Hard limit on kappa
          fld_kappa(ix^D) = min(fld_kappa(ix^D),2.3d0*fld_kappa0/unit_opacity)

          !> Limit kappa through T
          ! fld_kappa(ix^D) = fld_kappa0/unit_opacity &
          ! * (1.d0+10.d0**akram*w(ix^D,iw_rho)*unit_density &
          ! * (max(a2(ix^D),const_kB*5.9d4/(fld_mu*const_mp))/1.d12)**bkram)
        {enddo\ }

      case('opal')
        call phys_get_tgas(w,x,ixI^L,ixO^L,Temp)
        {do ix^D=ixOmin^D,ixOmax^D\ }
            rho0 = w(ix^D,iw_rho)*unit_density
            Temp0 = Temp(ix^D)*unit_temperature
            call set_opal_opacity(rho0,Temp0,n)
            fld_kappa(ix^D) = n/unit_opacity
        {enddo\ }

      case('special')
        if (.not. associated(usr_special_opacity)) then
          call mpistop("special opacity not defined")
        endif
        call usr_special_opacity(ixI^L, ixO^L, w, x, fld_kappa)

      case default
        call mpistop("Doesn't know opacity law")
      end select
  end subroutine fld_get_opacity

  !> Calculate fld flux limiter
  !> This subroutine calculates flux limiter lambda using the prescription
  !> stored in fld_fluxlimiter.
  !> It also calculates the ratio of radiation scaleheight and mean free path
  subroutine fld_get_fluxlimiter(w, x, ixI^L, ixO^L, fld_lambda, fld_R)
    use mod_global_parameters
    use mod_geometry
    use mod_usr_methods

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out) :: fld_R(ixO^S), fld_lambda(ixO^S)
    double precision :: kappa(ixO^S)
    double precision ::  normgrad2(ixO^S)
    double precision :: grad_r_e(ixI^S), grad_r_eO(ixO^S)
    double precision :: rad_e(ixI^S)
    integer :: idir, ix^D

    double precision :: tmp_L(ixI^S), filtered_L(ixI^S)
    integer :: filter, idim

    select case (fld_fluxlimiter)
    case('Diffusion')
      fld_lambda(ixO^S) = 1.d0/3.d0
      fld_R(ixO^S) = zero

    case('FreeStream')
      !> Calculate R everywhere
      !> |grad E|/(rho kappa E)
      normgrad2(ixO^S) = 0.d0 !smalldouble

      rad_e(ixI^S) = w(ixI^S, iw_r_e)
      do idir = 1,ndim
        call gradientO(rad_e,x,ixI^L,ixO^L,idir,grad_r_eO,nghostcells)
        normgrad2(ixO^S) = normgrad2(ixO^S) + grad_r_eO(ixO^S)**2

        ! call gradient(rad_e,ixI^L,ixO^L,idir,grad_r_e)
        ! normgrad2(ixO^S) = normgrad2(ixO^S) + grad_r_e(ixO^S)**2
      end do

      call fld_get_opacity(w, x, ixI^L, ixO^L, kappa)

      fld_R(ixO^S) = dsqrt(normgrad2(ixO^S))/(kappa(ixO^S)*w(ixO^S,iw_rho)*w(ixO^S,iw_r_e))

      !> Calculate the flux limiter, lambda
      fld_lambda(ixO^S) = one/fld_R(ixO^S)

    case('Pomraning')
      !> Calculate R everywhere
      !> |grad E|/(rho kappa E)
      normgrad2(ixO^S) = 0.d0 !smalldouble !*w(ixO^S,iw_r_e)

      rad_e(ixI^S) = w(ixI^S, iw_r_e)
      do idir = 1,ndim
        call gradientO(rad_e,x,ixI^L,ixO^L,idir,grad_r_eO,nghostcells)
        normgrad2(ixO^S) = normgrad2(ixO^S) + grad_r_eO(ixO^S)**2

        ! call gradient(rad_e,ixI^L,ixO^L,idir,grad_r_e)
        ! normgrad2(ixO^S) = normgrad2(ixO^S) + grad_r_e(ixO^S)**2
      end do

      call fld_get_opacity(w, x, ixI^L, ixO^L, kappa)

      fld_R(ixO^S) = dsqrt(normgrad2(ixO^S))/(kappa(ixO^S)*w(ixO^S,iw_rho)*w(ixO^S,iw_r_e))

      !> Calculate the flux limiter, lambda
      !> Levermore and Pomraning: lambda = (2 + R)/(6 + 3R + R^2)
      fld_lambda(ixO^S) = (2.d0+fld_R(ixO^S))/(6.d0+3*fld_R(ixO^S)+fld_R(ixO^S)**2.d0)

    case('Pomraning2')
      call mpistop("Pomraning2 is not quite working, use Pomraning or Minerbo")
      !> Calculate R everywhere
      !> |grad E|/(rho kappa E)
      normgrad2(ixO^S) = 0.d0 !smalldouble

      rad_e(ixI^S) = w(ixI^S, iw_r_e)
      do idir = 1,ndim
        call gradientO(rad_e,x,ixI^L,ixO^L,idir,grad_r_eO,nghostcells)
        normgrad2(ixO^S) = normgrad2(ixO^S) + grad_r_eO(ixO^S)**2

        ! call gradient(rad_e,ixI^L,ixO^L,idir,grad_r_e)
        ! normgrad2(ixO^S) = normgrad2(ixO^S) + grad_r_e(ixO^S)**2
      end do

      call fld_get_opacity(w, x, ixI^L, ixO^L, kappa)

      fld_R(ixO^S) = dsqrt(normgrad2(ixO^S))/(kappa(ixO^S)*w(ixO^S,iw_rho)*w(ixO^S,iw_r_e))

      !> Calculate the flux limiter, lambda
      !> Levermore and Pomraning: lambda = 1/R(coth(R)-1/R)
      fld_lambda(ixO^S) = one/fld_R(ixO^S)*(one/dtanh(fld_R(ixO^S)) - one/fld_R(ixO^S))
      ! fld_lambda(ixO^S) = one/fld_R(ixO^S)*(dcosh(fld_R(ixO^S))/dsinh(fld_R(ixO^S)) - one/fld_R(ixO^S))

      !>WHAT HAPPENS WHEN R=0 (full diffusion) => 1/R = NAN => dtanh(1/R) =????
      where(dtanh(fld_R(ixO^S)) .ne. dtanh(fld_R(ixO^S)))
        fld_lambda(ixO^S) = 1.d0/3.d0
      endwhere

    case('Minerbo')
      !> Calculate R everywhere
      !> |grad E|/(rho kappa E)
      normgrad2(ixO^S) = 0.d0 !smalldouble

      rad_e(ixI^S) = w(ixI^S, iw_r_e)
      do idir = 1,ndim
        call gradientO(rad_e,x,ixI^L,ixO^L,idir,grad_r_eO,nghostcells)
        normgrad2(ixO^S) = normgrad2(ixO^S) + grad_r_eO(ixO^S)**2

        ! call gradient(rad_e,ixI^L,ixO^L,idir,grad_r_e)
        ! normgrad2(ixO^S) = normgrad2(ixO^S) + grad_r_e(ixO^S)**2
      end do

      call fld_get_opacity(w, x, ixI^L, ixO^L, kappa)

      fld_R(ixO^S) = dsqrt(normgrad2(ixO^S))/(kappa(ixO^S)*w(ixO^S,iw_rho)*w(ixO^S,iw_r_e))

      !> Calculate the flux limiter, lambda
      !> Minerbo:
      {do ix^D=ixOmin^D,ixOmax^D\ }
          if (fld_R(ix^D) .lt. 3.d0/2.d0) then
            fld_lambda(ix^D) = 2.d0/(3.d0 + dsqrt(9.d0 + 12.d0*fld_R(ix^D)**2.d0))
          else
            fld_lambda(ix^D) = 1.d0/(1.d0 + fld_R(ix^D) + dsqrt(1.d0 + 2.d0*fld_R(ix^D)))
          endif
      {enddo\}

    case('special')
      if (.not. associated(usr_special_fluxlimiter)) then
        call mpistop("special fluxlimiter not defined")
      endif
      call usr_special_fluxlimiter(ixI^L, ixO^L, w, x, fld_lambda, fld_R)
    case default
      call mpistop('Fluxlimiter unknown')
    end select


    if (flux_lim_filter) then
      if (size_L_filter .lt. 1) call mpistop("D filter of size < 1 makes no sense")
      if (size_L_filter .gt. nghostcells) call mpistop("D filter of size > nghostcells makes no sense")

      tmp_L(ixO^S) = fld_lambda(ixO^S)
      filtered_L(ixO^S) = zero

      do filter = 1,size_L_filter
        {do ix^D = ixOmin^D+size_D_filter,ixOmax^D-size_L_filter\}
          do idim = 1,ndim
            filtered_L(ix^D) = filtered_L(ix^D) &
                             + tmp_L(ix^D+filter*kr(idim,^D)) &
                             + tmp_L(ix^D-filter*kr(idim,^D))
          enddo
        {enddo\}
      enddo

      {do ix^D = ixOmin^D+size_D_filter,ixOmax^D-size_D_filter\}
        tmp_L(ix^D) = (tmp_L(ix^D)+filtered_L(ix^D))/(1+2*size_L_filter*ndim)
      {enddo\}

      fld_lambda(ixO^S) = tmp_L(ixO^S)
    endif

  end subroutine fld_get_fluxlimiter

  !> Calculate Radiation Flux
  !> stores radiation flux in w-array
  subroutine fld_get_radflux(w, x, ixI^L, ixO^L, rad_flux)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out) :: rad_flux(ixO^S, 1:ndim)

    double precision :: grad_r_e(ixI^S), grad_r_eO(ixO^S)
    double precision :: rad_e(ixI^S)
    double precision :: kappa(ixO^S), lambda(ixO^S), fld_R(ixO^S)
    integer :: ix^D, idir

    rad_e(ixI^S) = w(ixI^S, iw_r_e)

    call fld_get_opacity(w, x, ixI^L, ixO^L, kappa)
    call fld_get_fluxlimiter(w, x, ixI^L, ixO^L, lambda, fld_R)

    !> Calculate the Flux using the fld closure relation
    !> F = -c*lambda/(kappa*rho) *grad E
    do idir = 1,ndim
      call gradientO(rad_e,x,ixI^L,ixO^L,idir,grad_r_eO,nghostcells)
      rad_flux(ixO^S, idir) = -(const_c/unit_velocity)*lambda(ixO^S)/(kappa(ixO^S)*w(ixO^S,iw_rho))*grad_r_eO(ixO^S)

      ! call gradient(rad_e,ixI^L,ixO^L,idir,grad_r_e)
      ! rad_flux(ixO^S, idir) = -(const_c/unit_velocity)*lambda(ixO^S)/(kappa(ixO^S)*w(ixO^S,iw_rho))*grad_r_e(ixO^S)
    end do

  end subroutine fld_get_radflux

  !> Calculate Eddington-tensor
  !> Stores Eddington-tensor in w-array
  subroutine fld_get_eddington(w, x, ixI^L, ixO^L, eddington_tensor)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out) :: eddington_tensor(ixO^S,1:ndim,1:ndim)
    double precision :: tnsr2(ixO^S,1:ndim,1:ndim)
    double precision :: normgrad2(ixO^S), f(ixO^S)
    double precision :: grad_r_e(ixI^S, 1:ndim), rad_e(ixI^S)
    double precision :: grad_r_eO(ixO^S, 1:ndim)
    double precision :: lambda(ixO^S), fld_R(ixO^S)
    integer :: i,j, idir,jdir

    !> Calculate R everywhere
    !> |grad E|/(rho kappa E)
    normgrad2(ixO^S) = zero

    rad_e(ixI^S) = w(ixI^S, iw_r_e)
    do idir = 1,ndim
      call gradientO(rad_e,x,ixI^L,ixO^L,idir,grad_r_eO(ixO^S, idir),nghostcells)
      normgrad2(ixO^S) = normgrad2(ixO^S) + grad_r_eO(ixO^S,idir)**2

      ! call gradient(rad_e,ixI^L,ixO^L,idir,grad_r_e(ixI^S,idir))
      ! normgrad2(ixO^S) = normgrad2(ixO^S) + grad_r_e(ixO^S,idir)**two
    end do

    call fld_get_fluxlimiter(w, x, ixI^L, ixO^L, lambda, fld_R)

    !> Calculate radiation pressure
    !> P = (lambda + lambda^2 R^2)*E
    f(ixO^S) = lambda(ixO^S) + lambda(ixO^S)**two * fld_R(ixO^S)**two
    f(ixO^S) = one/two*(one-f(ixO^S)) + one/two*(3.d0*f(ixO^S) - one)

    {^IFONED
    eddington_tensor(ixO^S,1,1) = f(ixO^S)
    }

    {^NOONED
    do idir = 1,ndim
      eddington_tensor(ixO^S,idir,idir) = half*(one-f(ixO^S))
    enddo

    do idir = 1,ndim
      do jdir = 1,ndim
        if (idir .ne. jdir) eddington_tensor(ixO^S,idir,jdir) = zero
        tnsr2(ixO^S,idir,jdir) =  half*(3.d0*f(ixO^S) - 1)&
        *grad_r_eO(ixO^S,idir)*grad_r_eO(ixO^S,jdir)/normgrad2(ixO^S)
      enddo
    enddo

    ! do idir = 1,ndim
    !   do jdir = 1,ndim
    !     if (idir .ne. jdir) eddington_tensor(ixO^S,idir,jdir) = zero
    !     tnsr2(ixO^S,idir,jdir) =  half*(3.d0*f(ixO^S) - 1)&
    !     *grad_r_e(ixO^S,idir)*grad_r_e(ixO^S,jdir)/normgrad2(ixO^S)
    !   enddo
    ! enddo

    do idir = 1,ndim
      do jdir = 1,ndim
        where ((tnsr2(ixO^S,idir,jdir) .eq. tnsr2(ixO^S,idir,jdir)) &
          .and. (normgrad2(ixO^S) .gt. smalldouble))
          eddington_tensor(ixO^S,idir,jdir) = eddington_tensor(ixO^S,idir,jdir) + tnsr2(ixO^S,idir,jdir)
        endwhere
      enddo
    enddo
    }

  end subroutine fld_get_eddington

  !> Calculate Radiation Pressure
  !> Returns Radiation Pressure as tensor
  subroutine fld_get_radpress(w, x, ixI^L, ixO^L, rad_pressure)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision             :: eddington_tensor(ixO^S,1:ndim,1:ndim)
    double precision, intent(out):: rad_pressure(ixO^S,1:ndim,1:ndim)

    integer i,j

    call fld_get_eddington(w, x, ixI^L, ixO^L, eddington_tensor)

    do i=1,ndim
      do j=1,ndim
        rad_pressure(ixO^S,i,j) = eddington_tensor(ixO^S,i,j)* w(ixO^S,iw_r_e)
      enddo
    enddo
  end subroutine fld_get_radpress

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!! Multigrid diffusion
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calling all subroutines to perform the multigrid method
  !> Communicates rad_e and diff_coeff to multigrid library
  subroutine Diffuse_E_rad_mg(dtfactor,qdt,qtC,psa,psb)
    use mod_global_parameters
    use mod_forest
    ! use mod_ghostcells_update
    use mod_multigrid_coupling
    ! use mod_physics, only: phys_set_mg_bounds, phys_req_diagonal

    type(state), target :: psa(max_blocks) !< Advance psa=psb+dtfactor*qdt*F_im(psa)
    type(state), target :: psb(max_blocks)
    double precision, intent(in) :: qdt
    double precision, intent(in) :: qtC
    double precision, intent(in) :: dtfactor

    integer                      :: n
    double precision             :: res, max_residual, lambda
    integer, parameter           :: max_its = 100

    integer                        :: iw_to,iw_from
    integer                        :: iigrid, igrid, id
    integer                        :: nc, lvl, i
    type(tree_node), pointer       :: pnode
    real(dp)                       :: fac, facD

    ! print*, '1'

    !> Set diffusion timestep, add previous timestep if mg did not converge:
    if (it == 0) dt_diff = 0
    dt_diff = dt_diff + qdt

    ! Avoid setting a very restrictive limit to the residual when the time step
    ! is small (as the operator is ~ 1/(D * qdt))
    if (qdt < dtmin) then
        if(mype==0)then
            print *,'skipping implicit solve: dt too small!'
            print *,'Currently at time=',global_time,' time step=',qdt,' dtmin=',dtmin
        endif
        return
    endif
    ! max_residual = 1d-7/qdt
    max_residual = fld_diff_tol !1d-7/qdt

    mg%operator_type = mg_vhelmholtz
    mg%smoother_type = mg_smoother_gs
    call mg_set_methods(mg)

    if (.not. mg%is_allocated) call mpistop("multigrid tree not allocated yet")

!   lambda = 1.d0/(dtfactor * qdt)
    lambda = 1.d0/(dtfactor * dt_diff)
    call vhelmholtz_set_lambda(lambda)

    call update_diffcoeff(psb)

    fac = 1.d0
    facD = 1.d0

    ! print*, '2'

    !This is mg_copy_to_tree from psb state
    call mg_copy_to_tree(i_diff_mg, mg_iveps, factor=facD, state_from=psb)
    !This is mg_copy_to_tree from psb state
    !!!  replaces::  call mg_copy_to_tree(su_, mg_irhs, factor=-lambda)
    ! iw_from=i_diff_mg
    ! iw_to=mg_iveps
    ! fac=facD
    ! do iigrid=1,igridstail; igrid=igrids(iigrid);
    !    pnode => igrid_to_node(igrid, mype)%node
    !    id    =  pnode%id
    !    lvl   =  mg%boxes(id)%lvl
    !    nc    =  mg%box_size_lvl(lvl)
    !    ! Include one layer of ghost cells on grid leaves
    !    {^IFONED
    !    mg%boxes(id)%cc(0:nc+1, iw_to) = fac * &
    !         psb(igrid)%w(ixMlo1-1:ixMhi1+1, iw_from)
    !    }
    !    {^IFTWOD
    !    mg%boxes(id)%cc(0:nc+1, 0:nc+1, iw_to) = fac * &
    !         psb(igrid)%w(ixMlo1-1:ixMhi1+1, ixMlo2-1:ixMhi2+1, iw_from)
    !    }
    !    {^IFTHREED
    !    mg%boxes(id)%cc(0:nc+1, 0:nc+1, 0:nc+1, iw_to) = fac * &
    !         psb(igrid)%w(ixMlo1-1:ixMhi1+1, ixMlo2-1:ixMhi2+1, &
    !         ixMlo3-1:ixMhi3+1, iw_from)
    !    }
    ! end do

    ! print*, '3'

    !This is mg_copy_to_tree from psb state
    call mg_copy_to_tree(iw_r_e, mg_iphi, factor=fac, state_from=psb)
    !This is mg_copy_to_tree from psb state
    !!!  replaces::  call mg_copy_to_tree(su_, mg_irhs, factor=-lambda)
    ! iw_from=iw_r_e
    ! iw_to=mg_iphi
    ! fac=fac
    ! do iigrid=1,igridstail; igrid=igrids(iigrid);
    !    pnode => igrid_to_node(igrid, mype)%node
    !    id    =  pnode%id
    !    lvl   =  mg%boxes(id)%lvl
    !    nc    =  mg%box_size_lvl(lvl)
    !    ! Include one layer of ghost cells on grid leaves
    !    {^IFONED
    !    mg%boxes(id)%cc(0:nc+1, iw_to) = fac * &
    !         psb(igrid)%w(ixMlo1-1:ixMhi1+1, iw_from)
    !    }
    !    {^IFTWOD
    !    mg%boxes(id)%cc(0:nc+1, 0:nc+1, iw_to) = fac * &
    !         psb(igrid)%w(ixMlo1-1:ixMhi1+1, ixMlo2-1:ixMhi2+1, iw_from)
    !    }
    !    {^IFTHREED
    !    mg%boxes(id)%cc(0:nc+1, 0:nc+1, 0:nc+1, iw_to) = fac * &
    !         psb(igrid)%w(ixMlo1-1:ixMhi1+1, ixMlo2-1:ixMhi2+1, &
    !         ixMlo3-1:ixMhi3+1, iw_from)
    !    }
    ! end do

    ! print*, '4'


    !>replace call set_rhs(mg, -1/dt, 0.0_dp)
    call mg_copy_to_tree(iw_r_e, mg_irhs, factor=-1/(dtfactor*dt_diff), state_from=psb)
    !This is mg_copy_to_tree from psb state
    !!!  replaces::  call mg_copy_to_tree(su_, mg_irhs, factor=-lambda)
    ! iw_from=iw_r_e
    ! iw_to=mg_irhs
    ! fac=-1/(dtfactor*dt_diff)
    ! do iigrid=1,igridstail; igrid=igrids(iigrid);
    !    pnode => igrid_to_node(igrid, mype)%node
    !    id    =  pnode%id
    !    lvl   =  mg%boxes(id)%lvl
    !    nc    =  mg%box_size_lvl(lvl)
    !    ! Include one layer of ghost cells on grid leaves
    !    {^IFONED
    !    mg%boxes(id)%cc(0:nc+1, iw_to) = fac * &
    !         psb(igrid)%w(ixMlo1-1:ixMhi1+1, iw_from)
    !    }
    !    {^IFTWOD
    !    mg%boxes(id)%cc(0:nc+1, 0:nc+1, iw_to) = fac * &
    !         psb(igrid)%w(ixMlo1-1:ixMhi1+1, ixMlo2-1:ixMhi2+1, iw_from)
    !    }
    !    {^IFTHREED
    !    mg%boxes(id)%cc(0:nc+1, 0:nc+1, 0:nc+1, iw_to) = fac * &
    !         psb(igrid)%w(ixMlo1-1:ixMhi1+1, ixMlo2-1:ixMhi2+1, &
    !         ixMlo3-1:ixMhi3+1, iw_from)
    !    }
    ! end do

    ! call phys_set_mg_bounds()

    ! print*, '5'


    if (time_advance)then
      call mg_restrict(mg, mg_iveps)
      call mg_fill_ghost_cells(mg, mg_iveps)
    endif

    ! print*, '6'


    call mg_fas_fmg(mg, .true., max_res=res)
    do n = 1, max_its
      !print*, n, res
      if (res < max_residual) exit
      call mg_fas_vcycle(mg, max_res=res)
    end do

    ! print*, '7'


    if (res .le. 0.d0) then
      if (diffcrash_resume) then
        if (mg%my_rank == 0) &
        write(*,*) it, ' resiudal zero ', res
        return
      endif
      if (mg%my_rank == 0) then
        print*, res
        error stop "Diffusion residual to zero"
      endif
    endif

    ! print*, '8'


    if (n == max_its + 1) then
      if (diffcrash_resume) then
        if (mg%my_rank == 0) &
        write(*,*) it, ' resiudal high ', res
        return
      endif
       if (mg%my_rank == 0) then
          print *, "Did you specify boundary conditions correctly?"
          print *, "Or is the variation in diffusion too large?"
          print*, n, res
          print *, mg%bc(1, mg_iphi)%bc_value, mg%bc(2, mg_iphi)%bc_value
       end if
       error stop "diffusion_solve: no convergence"
    end if

    ! print*, '9'


    !> Reset dt_diff when diffusion worked out
!0887 dt_diff = 0.d0
    dt_diff = 0.d0

    ! !This is mg_copy_from_tree_gc for psa state
    call mg_copy_from_tree_gc(mg_iphi, iw_r_e, state_to=psa)
    !This is mg_copy_from_tree_gc for psa state
    !!! replaces:: call mg_copy_from_tree_gc(mg_iphi, su_)
    ! iw_from=mg_iphi
    ! iw_to=iw_r_e
    ! do iigrid=1,igridstail; igrid=igrids(iigrid);
    !    pnode => igrid_to_node(igrid, mype)%node
    !    id    =  pnode%id
    !    lvl   =  mg%boxes(id)%lvl
    !    nc    =  mg%box_size_lvl(lvl)
    !    ! Include one layer of ghost cells on grid leaves
    !    {^IFONED
    !    psa(igrid)%w(ixMlo1-1:ixMhi1+1, iw_to) = &
    !         mg%boxes(id)%cc(0:nc+1, iw_from)
    !    }
    !    {^IFTWOD
    !    psa(igrid)%w(ixMlo1-1:ixMhi1+1, ixMlo2-1:ixMhi2+1, iw_to) = &
    !         mg%boxes(id)%cc(0:nc+1, 0:nc+1, iw_from)
    !    }
    !    {^IFTHREED
    !    psa(igrid)%w(ixMlo1-1:ixMhi1+1, ixMlo2-1:ixMhi2+1, &
    !         ixMlo3-1:ixMhi3+1, iw_to) = &
    !         mg%boxes(id)%cc(0:nc+1, 0:nc+1, 0:nc+1, iw_from)
    !    }
    ! end do

    ! print*, '10'


  end subroutine Diffuse_E_rad_mg


  !> inplace update of psa==>F_im(psa)
  subroutine Evaluate_E_rad_mg(qtC,psa)
    use mod_global_parameters
    use mod_ghostcells_update
    use mod_physics, only: phys_req_diagonal

    type(state), target :: psa(max_blocks)
    double precision, intent(in) :: qtC

    integer :: iigrid, igrid, level
    integer :: ixO^L

    call update_diffcoeff(psa)

    ixO^L=ixG^LL^LSUB1;
    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
       ! call fld_get_diffcoef_central(psa(igrid)%w, psa(igrid)%w, psa(igrid)%x, ixG^LL, ixO^L)
       call put_diffterm_onegrid(ixG^LL,ixO^L,psa(igrid)%w)
    end do
    !$OMP END PARALLEL DO

    ! enforce boundary conditions for psa
    call getbc(qtC,0.d0,psa,1,nwflux+nwaux,phys_req_diagonal)

  end subroutine Evaluate_E_rad_mg

  !> inplace update of psa==>F_im(psa)
  subroutine put_diffterm_onegrid(ixI^L,ixO^L,w)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, 1:nw)

    double precision :: gradE(ixO^S),divF(ixO^S)
    double precision :: divF_h(ixO^S),divF_j(ixO^S)
    double precision :: diff_term(ixO^S)

    integer                       :: idir, jxO^L, hxO^L

    ! call mpistop("phys_evaluate_implicit not implemented for FLD")

    divF(ixO^S) = 0.d0
    do idir = 1,ndim
      hxO^L=ixO^L-kr(idir,^D);
      jxO^L=ixO^L+kr(idir,^D);

      divF_h(ixO^S) = w(ixO^S,i_diff_mg)*w(hxO^S,i_diff_mg)/(w(ixO^S,i_diff_mg) + w(hxO^S,i_diff_mg))*(w(hxO^S,iw_r_e) - w(ixO^S,iw_r_e))
      divF_j(ixO^S) = w(ixO^S,i_diff_mg)*w(jxO^S,i_diff_mg)/(w(ixO^S,i_diff_mg) + w(jxO^S,i_diff_mg))*(w(jxO^S,iw_r_e) - w(ixO^S,iw_r_e))
      divF(ixO^S) = divF(ixO^S) + 2.d0*(divF_h(ixO^S) + divF_j(ixO^S))/dxlevel(idir)**2
    enddo

    w(ixO^S,iw_r_e) = divF(ixO^S)

  end subroutine put_diffterm_onegrid


  !> Calculates cell-centered diffusion coefficient to be used in multigrid
  subroutine fld_get_diffcoef_central(w, wCT, x, ixI^L, ixO^L)
    use mod_global_parameters
    use mod_geometry
    use mod_usr_methods

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: wCT(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)

    double precision :: kappa(ixO^S), lambda(ixO^S), fld_R(ixO^S)

    double precision :: max_D(ixI^S), grad_r_e(ixI^S), rad_e(ixI^S)
    integer :: idir,i,j, ix^D


    call fld_get_opacity(wCT, x, ixI^L, ixO^L, kappa)
    call fld_get_fluxlimiter(wCT, x, ixI^L, ixO^L, lambda, fld_R)

    !> calculate diffusion coefficient
    w(ixO^S,i_diff_mg) = (const_c/unit_velocity)*lambda(ixO^S)/(kappa(ixO^S)*wCT(ixO^S,iw_rho))

    where (w(ixO^S,i_diff_mg) .lt. 0.d0) &
      w(ixO^S,i_diff_mg) = smalldouble

    if (diff_coef_filter) then
      !call mpistop('Hold your bloody horses, not implemented yet ')
      call fld_smooth_diffcoef(w, ixI^L, ixO^L)
    endif

    if (associated(usr_special_diffcoef)) &
      call usr_special_diffcoef(w, wCT, x, ixI^L, ixO^L)

  end subroutine fld_get_diffcoef_central

  subroutine update_diffcoeff(psa)
    use mod_global_parameters

    type(state), target :: psa(max_blocks)

    ! double precision :: wCT(ixG^LL,1:nw)
    integer :: iigrid, igrid, level
    integer :: ixO^L

    ixO^L=ixG^LL^LSUB1;
    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);

       ! wCT = psa(igrid)%w
        call fld_get_diffcoef_central(psa(igrid)%w, psa(igrid)%w, psa(igrid)%x, ixG^LL, ixO^L)
    end do
    !$OMP END PARALLEL DO

  end subroutine update_diffcoeff

  !> Use running average on Diffusion coefficient
  subroutine fld_smooth_diffcoef(w, ixI^L, ixO^L)
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, 1:nw)

    double precision :: tmp_D(ixI^S), filtered_D(ixI^S)
    integer :: ix^D, filter, idim

    if (size_D_filter .lt. 1) call mpistop("D filter of size < 1 makes no sense")
    if (size_D_filter .gt. nghostcells) call mpistop("D filter of size > nghostcells makes no sense")

    tmp_D(ixO^S) = w(ixO^S,i_diff_mg)
    filtered_D(ixO^S) = zero

    do filter = 1,size_D_filter
      {do ix^D = ixOmin^D+size_D_filter,ixOmax^D-size_D_filter\}
        do idim = 1,ndim
          filtered_D(ix^D) = filtered_D(ix^D) &
                           + tmp_D(ix^D+filter*kr(idim,^D)) &
                           + tmp_D(ix^D-filter*kr(idim,^D))
        enddo
      {enddo\}
    enddo

    {do ix^D = ixOmin^D+size_D_filter,ixOmax^D-size_D_filter\}
      tmp_D(ix^D) = (tmp_D(ix^D)+filtered_D(ix^D))/(1+2*size_D_filter*ndim)
    {enddo\}

    w(ixO^S,i_diff_mg) = tmp_D(ixO^S)
  end subroutine fld_smooth_diffcoef


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!! Gas-Rad Energy interaction
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> This subroutine calculates the radiation heating, radiation cooling
  !> and photon tiring using an implicit scheme.
  !> These sourceterms are applied using the root-finding of a 4th order polynomial
  !> This routine loops over every cell in the domain
  !> and computes the coefficients of the polynomials in every cell
  subroutine Energy_interaction(w, wCT, x, ixI^L, ixO^L)
    use mod_global_parameters
    use mod_geometry
    use mod_physics
    use mod_usr_methods

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(in)    :: wCT(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: a1(ixO^S), a2(ixO^S)
    double precision :: c0(ixO^S), c1(ixO^S)
    double precision :: e_gas(ixO^S), E_rad(ixO^S)
    double precision :: kappa(ixO^S)
    double precision :: sigma_b

    integer :: i,j,idir,ix^D

    !> e_gas is the INTERNAL ENERGY without KINETIC ENERGY
    ! if (.not. block%e_is_internal) then
      e_gas(ixO^S) = wCT(ixO^S,iw_e) - half*sum(wCT(ixO^S, iw_mom(:))**2, dim=ndim+1)/wCT(ixO^S, iw_rho)
    ! else
    !   e_gas(ixO^S) = wCT(ixO^S,iw_e)
    ! endif

    {do ix^D=ixOmin^D,ixOmax^D\ }
      e_gas(ix^D) = max(e_gas(ix^D),small_pressure/(fld_gamma-1.d0))
    {enddo\}

    E_rad(ixO^S) = wCT(ixO^S,iw_r_e)

    if (associated(usr_special_opacity_qdot)) then
      call usr_special_opacity_qdot(ixI^L,ixO^L,w,x,kappa)
    else
      call fld_get_opacity(wCT, x, ixI^L, ixO^L, kappa)
    endif

    sigma_b = const_rad_a*const_c/4.d0*(unit_temperature**4.d0)/(unit_velocity*unit_pressure)

    if (fld_interaction_method .eq. 'Instant') then
      a1(ixO^S) = const_rad_a*(fld_mu*const_mp/const_kB*(fld_gamma-1))**4/(wCT(ixO^S,iw_rho)*unit_density)**4 &
                  /unit_pressure**3
      a2(ixO^S) = e_gas(ixO^S) + E_rad(ixO^S)

      c0(ixO^S) = a2(ixO^S)/a1(ixO^S)
      c1(ixO^S) = 1.d0/a1(ixO^S)
    else
      !> Calculate coefficients for polynomial
      a1(ixO^S) = 4*kappa(ixO^S)*sigma_b*(fld_gamma-one)**4/wCT(ixO^S,iw_rho)**3*dt
      a2(ixO^S) = (const_c/unit_velocity)*kappa(ixO^S)*wCT(ixO^S,iw_rho)*dt

      c0(ixO^S) = ((one + a2(ixO^S))*e_gas(ixO^S) + a2(ixO^S)*E_rad(ixO^S))/a1(ixO^S)
      c1(ixO^S) = (one + a2(ixO^S))/a1(ixO^S)
    endif

    !> Loop over every cell for rootfinding method
    {do ix^D=ixOmin^D,ixOmax^D\ }
      select case(fld_interaction_method)
      case('Bisect')
        call Bisection_method(e_gas(ix^D), E_rad(ix^D), c0(ix^D), c1(ix^D))
      case('Newton')
        call Newton_method(e_gas(ix^D), E_rad(ix^D), c0(ix^D), c1(ix^D))
      case('Halley')
        call Halley_method(e_gas(ix^D), E_rad(ix^D), c0(ix^D), c1(ix^D))
      case('Instant')
        call Halley_method(e_gas(ix^D), E_rad(ix^D), c0(ix^D), c1(ix^D))
      case default
        call mpistop('root-method not known')
      end select
    {enddo\}

    if (fld_interaction_method .eq. 'Instant') then
      E_rad(ixO^S) = a2(ixO^S) - e_gas(ixO^S)
    else
      !> advance E_rad
      E_rad(ixO^S) = (a1(ixO^S)*e_gas(ixO^S)**4.d0 + E_rad(ixO^S))/(one + a2(ixO^S))
    endif

    !> new w = w + dt f(wCT)
    !> e_gas,E_rad = wCT + dt f(wCT)
    !> dt f(wCT) = e_gas,E_rad - wCT
    !> new w = w +  e_gas,E_rad - wCT

    !> WAIT A SECOND?! DOES THIS EVEN WORK WITH THESE KIND OF IMPLICIT METHODS?
    !> NOT QUITE SURE TO BE HONEST
    !> IS IT POSSIBLE TO SHUT DOWN SOURCE SPLITTING FOR JUST THIS TERM?
    !> FIX BY PERFORMING Energy_interaction on (w,w,...)

    ! call mpistop('This still has to be fixed somehow')

    !> Update gas-energy in w, internal + kinetic
    ! w(ixO^S,iw_e) = w(ixO^S,iw_e) + e_gas(ixO^S) - wCT(ixO^S,iw_e)
    ! {do ix^D=ixOmin^D,ixOmax^D\ }
    !   e_gas(ix^D) = max(e_gas(ix^D),small_pressure/(fld_gamma-1.d0))
    ! {enddo\}

    !{do ix^D=ixOmin^D,ixOmax^D\ }
    !  w(ix^D,iw_e) = max(e_gas(ix^D),small_pressure/(fld_gamma - 1))
    !{enddo\}
    w(ixO^S,iw_e) = e_gas(ixO^S)

    !> Beginning of module substracted wCT Ekin
    !> So now add wCT Ekin
    ! if (.not. block%e_is_internal) then
    w(ixO^S,iw_e) = w(ixO^S,iw_e) + half*sum(wCT(ixO^S, iw_mom(:))**2, dim=ndim+1)/wCT(ixO^S, iw_rho)
    ! else
    !   w(ixO^S,iw_e) = w(ixO^S,iw_e)
    ! endif

    !> Update rad-energy in w
    ! w(ixO^S,iw_r_e) = w(ixO^S,iw_r_e) + E_rad(ixO^S) - wCT(ixO^S,iw_r_e)
    w(ixO^S,iw_r_e) = E_rad(ixO^S)

  end subroutine Energy_interaction


  !> Find the root of the 4th degree polynomial using the bisection method
  subroutine Bisection_method(e_gas, E_rad, c0, c1)
    use mod_global_parameters

    double precision, intent(in)    :: c0, c1
    double precision, intent(in)    :: E_rad
    double precision, intent(inout) :: e_gas

    double precision :: bisect_a, bisect_b, bisect_c
    integer :: n, max_its

    n = 0
    max_its = 100

    bisect_a = zero
    bisect_b = min(abs(c0/c1),abs(c0)**(1.d0/4.d0))+smalldouble

    ! do while (abs(Polynomial_Bisection(bisect_b, c0, c1)-Polynomial_Bisection(bisect_a, c0, c1))&
    !    .ge. fld_bisect_tol*min(e_gas,E_rad))
    do while (abs(bisect_b-bisect_a) .ge. fld_bisect_tol*min(e_gas,E_rad))
      bisect_c = (bisect_a + bisect_b)/two

      n = n +1
      if (n .gt. max_its) then
        goto 2435
        call mpistop('No convergece in bisection scheme')
      endif

      if (Polynomial_Bisection(bisect_a, c0, c1)*&
      Polynomial_Bisection(bisect_b, c0, c1) .lt. zero) then

        if (Polynomial_Bisection(bisect_a, c0, c1)*&
        Polynomial_Bisection(bisect_c, c0, c1) .lt. zero) then
          bisect_b = bisect_c
        elseif (Polynomial_Bisection(bisect_b, c0, c1)*&
        Polynomial_Bisection(bisect_c, c0, c1) .lt. zero) then
          bisect_a = bisect_c
        elseif (Polynomial_Bisection(bisect_a, c0, c1) .eq. zero) then
          bisect_b = bisect_a
          bisect_c = bisect_a
          goto 2435
        elseif (Polynomial_Bisection(bisect_b, c0, c1) .eq. zero) then
          bisect_a = bisect_b
          bisect_c = bisect_b
          goto 2435
        elseif (Polynomial_Bisection(bisect_c, c0, c1) .eq. zero) then
          bisect_a = bisect_c
          bisect_b = bisect_c
          goto 2435
        else
          call mpistop("Problem with fld bisection method")
        endif
      elseif (Polynomial_Bisection(bisect_a, c0, c1) &
        - Polynomial_Bisection(bisect_b, c0, c1) .lt. fld_bisect_tol*Polynomial_Bisection(bisect_a, c0, c1)) then
        goto 2435
      else
        bisect_a = e_gas
        bisect_b = e_gas
        print*, "IGNORING GAS-RAD ENERGY EXCHANGE ", c0, c1

        print*, Polynomial_Bisection(bisect_a, c0, c1), Polynomial_Bisection(bisect_b, c0, c1)

        if (Polynomial_Bisection(bisect_a, c0, c1) .le. smalldouble) then
          bisect_b = bisect_a
        elseif (Polynomial_Bisection(bisect_a, c0, c1) .le. smalldouble) then
          bisect_a = bisect_b
        endif

        goto 2435

      endif
    enddo

      2435 e_gas = (bisect_a + bisect_b)/two
  end subroutine Bisection_method

  !> Find the root of the 4th degree polynomial using the Newton-Ralphson method
  subroutine Newton_method(e_gas, E_rad, c0, c1)
    use mod_global_parameters

    double precision, intent(in)    :: c0, c1
    double precision, intent(in)    :: E_rad
    double precision, intent(inout) :: e_gas

    double precision :: xval, yval, der, deltax

    integer :: ii

    yval = bigdouble
    xval = e_gas
    der = one
    deltax = one

    ii = 0
    !> Compare error with dx = dx/dy dy
    do while (abs(deltax) .gt. fld_bisect_tol)
      yval = Polynomial_Bisection(xval, c0, c1)
      der = dPolynomial_Bisection(xval, c0, c1)
      deltax = -yval/der
      xval = xval + deltax
      ii = ii + 1
      if (ii .gt. 1d3) then
        print*, 'skip to bisection algorithm'
        call Bisection_method(e_gas, E_rad, c0, c1)
        return
      endif
    enddo

    e_gas = xval
  end subroutine Newton_method

  !> Find the root of the 4th degree polynomial using the Halley method
  subroutine Halley_method(e_gas, E_rad, c0, c1)
    use mod_global_parameters

    double precision, intent(in)    :: c0, c1
    double precision, intent(in)    :: E_rad
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
    do while (abs(deltax) .gt. fld_bisect_tol)
      yval = Polynomial_Bisection(xval, c0, c1)
      der = dPolynomial_Bisection(xval, c0, c1)
      dder = ddPolynomial_Bisection(xval, c0, c1)
      deltax = -two*yval*der/(two*der**2 - yval*dder)
      xval = xval + deltax
      ii = ii + 1
      if (ii .gt. 1d3) then
        ! call mpistop('Halley did not convergggge')
        call Newton_method(e_gas, E_rad, c0, c1)
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

  !> Calculate gradient of a scalar q within ixL in direction idir
  !> difference with gradient is gradq(ixO^S), NOT gradq(ixI^S)
  subroutine gradientO(q,x,ixI^L,ixO^L,idir,gradq,n)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L, idir
    integer, intent(in)             :: n

    double precision, intent(in)    :: q(ixI^S), x(ixI^S,1:ndim)
    double precision, intent(out)   :: gradq(ixO^S)
    integer                         :: jxO^L, hxO^L

    ! hxO^L=ixO^L-n*kr(idir,^D);
    ! jxO^L=ixO^L+n*kr(idir,^D);
    !
    ! if (n .gt. nghostcells) call mpistop("gradientO stencil too wide")
    !
    ! gradq(ixO^S)=(q(jxO^S)-q(hxO^S))/(2*n*dxlevel(idir))
    ! gradq(ixO^S)=(q(jxO^S)-q(hxO^S))/(x(jxO^S,idir)-x(hxO^S,idir))

    !> Using higher order derivatives with wider stencil according to:
    !> https://en.wikipedia.org/wiki/Finite_difference_coefficient

     if (n .gt. nghostcells) then
       call mpistop("gradientO stencil too wide")
     elseif (n .eq. 1) then
       hxO^L=ixO^L-kr(idir,^D);
       jxO^L=ixO^L+kr(idir,^D);
       gradq(ixO^S)=(q(jxO^S)-q(hxO^S))/(x(jxO^S,idir)-x(hxO^S,idir))
     elseif (n .eq. 2) then
       gradq(ixO^S) = 0.d0
       !> coef 2/3
       hxO^L=ixO^L-kr(idir,^D);
       jxO^L=ixO^L+kr(idir,^D);
       gradq(ixO^S) = gradq(ixO^S) + 2.d0/3.d0*(q(jxO^S)-q(hxO^S))
       !> coef -1/12
       hxO^L=ixO^L-2*kr(idir,^D);
       jxO^L=ixO^L+2*kr(idir,^D);
       gradq(ixO^S) = gradq(ixO^S) - 1.d0/12.d0*(q(jxO^S)-q(hxO^S))
       !> divide by dx
       gradq(ixO^S) = gradq(ixO^S)/dxlevel(idir)
     elseif (n .eq. 3) then
       gradq(ixO^S) = 0.d0
       !> coef 3/4
       hxO^L=ixO^L-kr(idir,^D);
       jxO^L=ixO^L+kr(idir,^D);
       gradq(ixO^S) = gradq(ixO^S) + 3.d0/4.d0*(q(jxO^S)-q(hxO^S))
       !> coef -3/20
       hxO^L=ixO^L-2*kr(idir,^D);
       jxO^L=ixO^L+2*kr(idir,^D);
       gradq(ixO^S) = gradq(ixO^S) - 3.d0/20.d0*(q(jxO^S)-q(hxO^S))
       !> coef 1/60
       hxO^L=ixO^L-3*kr(idir,^D);
       jxO^L=ixO^L+3*kr(idir,^D);
       gradq(ixO^S) = gradq(ixO^S) + 1.d0/60.d0*(q(jxO^S)-q(hxO^S))
       !> divide by dx
       gradq(ixO^S) = gradq(ixO^S)/dxlevel(idir)
     else
       call mpistop("gradientO stencil unknown")
     endif

  end subroutine gradientO

end module mod_fld
