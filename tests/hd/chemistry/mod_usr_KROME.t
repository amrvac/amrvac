!> This is a template for a new user problem
module mod_usr

  ! Include a physics module
  use mod_hd

  implicit none

  ! Custom variables can be defined here
  ! ...
  ! Hydrodynamic variables
  real(dp) :: temperature_star = 2.5e3_dp ! in Kelvin
  real(dp) :: radius_star = au_SI ! in AU
  real(dp) :: mass_star = mSun_SI ! in kg
  real(dp) :: density_star = 1.e-6_dp ! in kg
  real(dp) :: pulsation_period = 300._dp * day ! in s
  real(dp) :: pulsation_amplitude = 2.5e3_dp ! in m/s

  ! Chemistry variables
  real(dp) :: cosmic_ray_rate = 1.36e-17_dp ! in s^-1
  real(dp) :: He_init = 3.11e-1_dp ! mass fraction
  real(dp) :: C_init = 2.63e-3_dp ! mass fraction
  real(dp) :: N_init = 1.52e-3_dp ! mass fraction
  real(dp) :: O_init = 9.60e-3_dp ! mass fraction
  real(dp) :: S_init = 3.97e-4_dp ! mass fraction
  real(dp) :: Fe_init = 1.17e-3_dp ! mass fraction
  real(dp) :: Si_init = 6.54e-4_dp ! mass fraction
  real(dp) :: Mg_init = 5.16e-4_dp ! mass fraction
  real(dp) :: Na_init = 3.38e-5_dp ! mass fraction
  real(dp) :: P_init = 8.17e-6_dp ! mass fraction
  real(dp) :: F_init = 4.06e-7_dp ! mass fraction
  real(dp) :: H_init = 1.0_dp - He_init - C_init - N_init - O_init - S_init &
   - Fe_init - Si_init - Mg_init - Na_init - P_init - F_init  ! mass fraction

  ! Local variables indices
  ! TODO: make them also public, protected????
  integer :: gamma
  integer :: temperature

contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()
    use mod_usr_methods
    use mod_global_parameters
    use krome_main
    use krome_user

    character(len=16) :: species(krome_nmols)

    ! Choose coordinate system as 1D Cartesian
    ! NOTE: not necessary when typeaxial is defined in .par file
    call set_coordinate_system("Cartesian_1D")
    ! Initialize krome
    call krome_init()
    ! Set krome variables
    call krome_set_user_crflux(cosmic_ray_rate)

    ! A routine for initial conditions is always required
    usr_init_one_grid => initial_conditions
    hd_usr_gamma      => get_gamma
    usr_process_grid  => calculate_local_variables
    special_bc        => special_boundary
    usr_gravity       => gravitation_acceleration
    usr_get_dt        => pulsation_get_dt
    usr_source        => chemical_evolution
    usr_add_aux_names => additional_output_names
    usr_aux_output    => additional_output
    ! Specify other user routines, for a list see mod_usr_methods.t
    ! ...

    ! Add own variables
    gamma = var_set_extravar("gamma", "gamma")
    temperature = var_set_extravar("temperature", "temperature")

    ! Active the physics module
    call hd_activate()

    ! get names of all species (incl dummy's of krome)
    species = krome_get_names()

    ! Rename tracers
    do idx = 1, hd_n_tracer
       prim_wnames(tracer(idx)) = species(idx)
    end do

  end subroutine usr_init


  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixI^L, ixO^L, w, x)
    use mod_global_parameters
    use krome_user
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    real(dp) :: temperature_init(ixI^S), mu_init
    real(dp) :: alpha, beta
    integer :: ix^D

    ! Set initial values for w
    ! w(ix^S, :) = ...
    alpha = 10.0_dp
    beta = 0.5_dp
    ! initial density profile
    w(ixI^S, rho_) = density_star * ( x(ixI^S,1) )**( -alpha )
    ! initial temperature profile
    temperature_init(ixI^S) = temperature_star * ( x(ixI^S,1) )**( -beta )
    ! initial velocity profile
    w(ixI^S, mom(:) ) = 0.0_dp


    ! set all number densities to zero
    w(ixI^S, tracer(:)) = 0.0_dp
    ! set initial mass fractions
    w(ixI^S, tracer(krome_idx_H)) = H_init
    w(ixI^S, tracer(krome_idx_HE)) = He_init
    w(ixI^S, tracer(krome_idx_C)) = C_init
    w(ixI^S, tracer(krome_idx_N)) = N_init
    w(ixI^S, tracer(krome_idx_O)) = O_init
    w(ixI^S, tracer(krome_idx_S)) = S_init
    w(ixI^S, tracer(krome_idx_FE)) = Fe_init
    w(ixI^S, tracer(krome_idx_SI)) = Si_init
    w(ixI^S, tracer(krome_idx_MG)) = Mg_init
    w(ixI^S, tracer(krome_idx_NA)) = Na_init
    w(ixI^S, tracer(krome_idx_P)) = P_init
    w(ixI^S, tracer(krome_idx_F)) = F_init

    ! loop over all cells
    {do ix^DB = ixI^LIM^DB\}

      ! convert to number densities (in cgs)
      w(ix^D, tracer(:)) = krome_x2n( w(ix^D, tracer(:)), w(ix^D, rho_)*1.e-3_dp )
      ! get gamma
      w(ix^D, gamma) = krome_get_gamma( w(ix^D, tracer(:)), temperature_init(ix^D))

    {enddo^D&\}

    ! get mean molecular weight
    ! is everywhere the same due one initial composition
    ! TODO: make space dependent if initial composition is as well
    mu_init = krome_get_mu( w(ixImin1, tracer(:)) )

    ! initial pressure profile
    w(ixI^S, p_)	= w(ixI^S, rho_) * temperature_init(ixI^S) &
                  * kB_SI /( mu_init * mp_SI )

    call hd_to_conserved(ixI^L, ixO^L, w, x)

  end subroutine initial_conditions

  ! Extra routines can be placed here
  ! ...

  !> returns current local value for Gamma
  subroutine get_gamma(w, ixI^L, ixO^L, var)
    use mod_global_parameters, only: nw, ndim
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S, nw)
    double precision, intent(out) :: var(ixO^S)

    var(ixO^S) = w(ixO^S, gamma)

  end subroutine get_gamma


  !> calculates local variables: temperature and gamma
  !> Note that temperature and gamme depend on each other
  !> therfore the temperature is calculated with the previous value of gamma
  !> this temperature is then used to calculate a new value of gamme
  !> which also depends on the new chemical composition (tracers).
  !> This method is valid when gamma locally doesn't change significantly,
  !> from the previous time step, which is okay because the possible
  !> range of gamma is quite narrow [1.0, 1.67] and an abrupt local change
  !> is highly unlikely.
  subroutine calculate_local_variables(igrid,level,ixI^L,ixO^L,qt,w,x)
    use mod_global_parameters
    use krome_user
    integer, intent(in)             :: igrid,level,ixI^L,ixO^L
    double precision, intent(in)    :: qt,x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    real(dp) :: e_internal(ixI^S), mu(ixI^S)
    integer :: ix^D

    ! get internal energy
    ! TODO: this can be implemented in mod_hd_phys just like hd_kin_en
    e_internal(ixI^S) = w(ixI^S,e_) - hd_kin_en(w, ixI^L, ixO^L)

    ! loop over all cells
    {do ix^DB = ixI^LIM^DB\}
      ! calculate mean molecular weight (tracer in cgs)
    	mu(ix^D) = krome_get_mu( w(ix^D, tracer(:) ) )

      ! calculate temperature
      w(ix^D, temperature) = (e_internal(ix^D) / w(ix^D, rho_) ) &
           * ( w(ix^D, gamma) - 1.0_dp ) * mu(ix^D)) &
           * ( mp_SI / kB_SI )

      ! calculate gamma (tracer in cgs)
      w(ix^D, gamma) = krome_get_gamma( w(ix^D, tracer(:)), w(ix^D, temperature))

    {enddo^D&\}

  end subroutine calculate_local_variables


  subroutine special_boundary(qt,ixI^L,ixO^L,iB,w,x)
    use mod_global_parameters
    use krome_user
    integer, intent(in)             :: ixI^L, ixO^L, iB
    double precision, intent(in)    :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    real(dp) :: mu_star

    ! density at boundary
    w(ixI^S, rho_) = density_star
    ! temperature at boundary
    temperature_bc(ixI^S) = temperature_star
    ! velocity at boundary
    w(ixI^S, mom(:) ) = 0.0_dp
    w(ixI^S, mom(1) ) = pulsation_amplitude * &
                    dsin( 2.0_dpr * dpi * qt / pulsation_period )

    ! set all number densities to zero
    w(ixI^S, tracer(:)) = 0.0_dp
    ! set initial mass fractions
    w(ixI^S, tracer(krome_idx_H)) = H_init
    w(ixI^S, tracer(krome_idx_HE)) = He_init
    w(ixI^S, tracer(krome_idx_C)) = C_init
    w(ixI^S, tracer(krome_idx_N)) = N_init
    w(ixI^S, tracer(krome_idx_O)) = O_init
    w(ixI^S, tracer(krome_idx_S)) = S_init
    w(ixI^S, tracer(krome_idx_FE)) = Fe_init
    w(ixI^S, tracer(krome_idx_SI)) = Si_init
    w(ixI^S, tracer(krome_idx_MG)) = Mg_init
    w(ixI^S, tracer(krome_idx_NA)) = Na_init
    w(ixI^S, tracer(krome_idx_P)) = P_init
    w(ixI^S, tracer(krome_idx_F)) = F_init

    ! loop over all cells
    {do ix^DB = ixI^LIM^DB\}

      ! convert to number densities (in cgs)
      w(ix^D, tracer(:)) = krome_x2n( w(ix^D, tracer(:)), w(ix^D, rho_)*1.e-3_dp )
      ! get gamma
      w(ix^D, gamma) = krome_get_gamma( w(ix^D, tracer(:)), temperature_bc(ix^D))

    {enddo^D&\}

    ! get mean molecular weight
    ! is everywhere the same due one initial composition
    ! TODO: make space dependent if initial composition is as well
    mu_bc = krome_get_mu( w(ixImin1, tracer(:)) )

    ! initial pressure profile
    w(ixI^S, p_)	= w(ixI^S, rho_) * temperature_init(ixI^S) &
                  * kB_SI /( mu_bc * mp_SI )

    call hd_to_conserved(ixI^L, ixO^L, w, x)

  end subroutine special_boundary


  !> Calculate gravitational acceleration in each dimension
  subroutine gravitation_acceleration(ixI^L,ixO^L,wCT,x,gravity_field)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(in)    :: wCT(ixI^S,1:nw)
    double precision, intent(out)   :: gravity_field(ixI^S,ndim)

    ! gravitational acceleration
    ! 1D only in radial direction
    gravity_field(ixO^S, 1) = -grav_SI * mass_star /( (x(ixO^S,1)**2.0_dp))

  end subroutine gravitation_acceleration

  !> Limit "dt" further if necessary, e.g. due to the special source terms.
  !> The getdt_courant (CFL condition) and the getdt subroutine in the AMRVACPHYS
  !> module have already been called.
  subroutine pulsation_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: dx^D, x(ixI^S,1:ndim)
    double precision, intent(in)    :: w(ixI^S,1:nw)
    double precision, intent(inout) :: dtnew

    ! make sure the pulsation at the inner boundary is sampled properly
    dtnew = pulsation_period * 1e-2_dp

  end subroutine pulsation_get_dt

  !> Calculate w(iw)=w(iw)+qdt*SOURCE[wCT,qtC,x] within ixO for all indices
  !> iw=iwmin...iwmax.  wCT is at time qCT
  subroutine chemical_evolution(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)
    use mod_global_parameters
    use krome_main
    use krome_user
    integer, intent(in)             :: ixI^L, ixO^L, iw^LIM
    double precision, intent(in)    :: qdt, qtC, qt
    double precision, intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    integer :: ix^D
    real(dp) :: temperature_local
    real(dp) :: tracer_local(1:hd_n_tracer)
    real(dp) :: mu_old, mu_new
    real(dp) :: delta_energy !< change in energy due to new chemical composition
    real(dp) :: delta_tracer(1:hd_n_tracer) !< change in chemical composition
    ! if we have soure splitting also AFTER flux calculations (e.g. 'sfs')
    ! then we need another call to 'calculate_local_variables'
    ! but this is not a 'usr_process_grid'

    !> loop over all grid cells
    {do ix^DB = ixO^LIM^DB\}

      !TODO: is it better/faster to define an old and new local variable instead
      ! of calling the w-array and wCT-array?
      tracer_local = wCT(ix^D, tracer(:))
      temperature_local = wCT(ix^D, temperature)
      mu_old = krome_get_mu(tracer_local(:))

      !> update tracers and temperature by solving chemistry over qdt
      !> returns updated values for input
      call krome(tracer_local(:), temperature_local, qdt)

      ! update local variables
      w(ix^D, gamma) = krome_get_gamma(tracer_local(:), temperature_local)
      w(ix^D, temperature) = temperature_local
      mu_new = krome_get_mu(tracer_local(:))

      !delta =  new - old
      delta_tracer = tracer_local - wCT(ix^D, tracer(:))
      delta_energy = wCT(ix^D, rho_) * kB_SI/mp_SI * &
                ( temperature_local/(mu_new * (w(ix^D, gamma) - 1.0_dp) ) &
                - wCT(ix^D, temperature)/(mu_old * (wCT(ix^D, gamma)) - 1.0_dp))

      !w(iw) = w(iw) + delta_w
      w(ix^D, tracer(:)) = w(ix^D, tracer(:)) + delta_tracer
      w(ix^D, e_) = w(ix^D, e_) + delta_energy

    {enddo^D&\}

  end subroutine chemical_evolution


  !TODO: doesn't this need 'intent(out)' ??
  subroutine additional_output(ixI^L,ixO^L,w,x,normconv)
    use mod_global_parameters
    use krome_user
    integer, intent(in)          :: ixI^L,ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision             :: w(ixI^S,nw+nwauxio)
    double precision             :: normconv(0:nw+nwauxio)

    real(dp) :: mass_fractions(1:hd_n_tracer)

    !loop over all grid cells
    {do ix^DB = ixO^LIM^DB\}

        mass_fractions(:) = krome_n2x(w(ix^D, tracer(:), w(ix^D,rho_)*1.0e-3_dp)
        cooling_rates(:) = krome_get_cooling_array(w(ix^D, tracer(:)),w(ix^D, temperature))
        heating_rates(:) = krome_get_heating_array(w(ix^D, tracer(:)),w(ix^D, temperature))

        w(ix^D,nw+1)= sum(mfracs(:))
        w(ix^D,nw+2)= krome_get_mu(w(ix^D, tracer(:)))
        w(ix^D,nw+3) = cooling_rates(krome_idx_cool_h2)
        w(ix^D,nw+4) = cooling_rates(krome_idx_cool_atomic)
        w(ix^D,nw+5) = cooling_rates(krome_idx_cool_z)
        w(ix^D,nw+6) = cooling_rates(krome_idx_cool_compton)
        w(ix^D,nw+7) = cooling_rates(krome_idx_cool_cont)
        w(ix^D,nw+8) = cooling_rates(krome_idx_cool_co)
        w(ix^D,nw+9) = heating_rates(krome_idx_heat_chem)

    {enddo^D&\}

  end subroutine additional_output


  !> Add names for the auxiliary variables space-separated
  !TODO: doesn't this need 'intent(inout)' ??
  subroutine additional_output_names(varnames)
    use mod_global_parameters
    character(len=*) :: varnames

    varnames = 'ntot mu coolH2 coolAtom coolZ coolCIE coolCO heatChem heatCR'


  end subroutine additional_output_names

  subroutine krome_modify_wLR(wLC, wLp, wRC, wRp, ixI^L, ixO^L, idir)
    use mod_global_parameters
    integer, intent(in)                :: ixI^L, ixO^L, idir
    double precision, intent(inout)    :: wLp(ixI^S,1:nw), wRp(ixI^S,1:nw)
    double precision, intent(inout)    :: wLC(ixI^S,1:nw), wRC(ixI^S,1:nw)

    print *, "Hello"
  end subroutine krome_modify_wLR

end module mod_usr
