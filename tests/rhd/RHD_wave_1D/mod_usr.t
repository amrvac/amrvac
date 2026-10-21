module mod_usr

  ! Include a physics module
  use mod_hd
  use mod_eos, only: eos
  use mod_fld

  implicit none

  ! input variables with units
  double precision :: rho0,eg0
  ! input variables that are dimensionless
  double precision :: rho_kappa_lambda,A_rho

  ! derived vars
  double precision :: omega,wavenumber,csound_norm,wvlength
  double precision :: A_v, A_p
  double precision :: p0,rho0_norm,eg0_norm,Er0_norm,p0_norm,T0_norm

  ! Storing additional var in the dat file
  integer :: Tgas_,Trad_,pres_,vel_,amr_


contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()

    ! Note how we here must set three values that in turn define M-L-T
    unit_density       = 3.216d-9     ! g cm^-3
    unit_length        = 7.7736318d11 ! cm
    unit_velocity      = 2.32d6       ! cm/s

    call usr_params_read(par_files)

    ! computing parameters and setting Boundary conditions for MG solver
    usr_set_parameters => set_params_and_mg_boundary_conds

    ! A routine for initial conditions is always required
    usr_init_one_grid => initial_conditions

   ! to add selected variables to the .dat file
    usr_modify_output => set_output_vars

    ! Drive the wave using an internal boundary
    usr_internal_bc => Initialize_Wave

    ! Choose coordinate system as n-D Cartesian
    {^IFONED call set_coordinate_system("Cartesian_1D")}
    {^IFTWOD call set_coordinate_system("Cartesian_2D")}
    {^IFTHREED call set_coordinate_system("Cartesian_3D")}

    ! Active the physics module
    call hd_activate()

    ! to add selected variables to the .dat file
    Tgas_ = var_set_extravar("Tgas", "Tgas")
    Trad_ = var_set_extravar("Trad", "Trad")
    pres_ = var_set_extravar("pres", "pres")
    vel_ = var_set_extravar("vel", "vel")
    amr_ = var_set_extravar("level", "level")

  end subroutine usr_init


  !> Read parameters from a file
  subroutine usr_params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /wave_list/ rho0, eg0, rho_kappa_lambda, A_rho

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       rewind(unitpar)
       read(unitpar, wave_list, end=113)
113    close(unitpar)
    end do

    p0 = eg0*(eos%gamma - one)
    if(mype==0)then
    print *,'============================================'
    write(*,*) 'INPUT GIVEN IN cgs units is'
    write(*,*) 'input density, gas energy density=',rho0,eg0
    write(*,*) 'derived gas pressure is =',p0
    write(*,*) 'input dimensionless amplitude for density is   =',A_rho
    write(*,*) 'input dimensionless parameter for wavelength   =',rho_kappa_lambda
    print *,'============================================'
    endif

  end subroutine usr_params_read

  subroutine set_params_and_mg_boundary_conds()
    use mod_global_parameters
    use mod_fld

    ! here we normalize all input values and compute the equilibrium parameters
    rho0_norm  = rho0/unit_density
    p0_norm   = p0/unit_pressure
    eg0_norm = eg0/unit_pressure
    T0_norm   = p0_norm/(rho0_norm*RR)
    Er0_norm = arad_norm*T0_norm**4.d0
    
    csound_norm=dsqrt(eos%gamma*p0/rho0)/unit_velocity

    if(mype==0)then
    write(*,*) 'derived dimensionless rad/gas energy ratio =',Er0_norm/eg0_norm
    write(*,*) 'derived ratio r=E/(4gamma e) =',Er0_norm/(4.d0*eos%gamma*eg0_norm)
    write(*,*) 'derived Boltzmann ratio Bo=4 gamma cs e/(cE) =',4.0d0*eos%gamma*csound_norm*eg0_norm/(c_norm*Er0_norm)
    endif

    ! here we compute the wave-related parameters
    ! all computed things here are dimensionless
    select case(trim(fld_opacity_law))
      case('const_norm')
         wvlength=rho_kappa_lambda/(rho0_norm*fld_kappa0)
      case('const')
          wvlength=rho_kappa_lambda/(rho0*fld_kappa0*unit_length)
      case default
         call mpistop("unknown opacity law")
    end select
    wavenumber=2.0d0*dpi/wvlength
    omega=wavenumber*csound_norm
    A_v   = A_rho*omega/wavenumber
    A_p   = eos%gamma*p0_norm*A_rho

    if(mype==0)then
    print *,'wavelength in physical units is=',wvlength*unit_length
    write(*,*) 'derived dimensionless wavelength for wave =',wvlength
    write(*,*) 'derived dimensionless omega and k for wave =',omega,wavenumber
    write(*,*) 'check on input and derivation based on ',rho_kappa_lambda,rho0,unit_density,fld_kappa0
    write(*,*) 'dispersion relation asks for omega=',omega,' to equal k= ',wavenumber,' times cs0=',csound_norm
    write(*,*) 'and sets amplitude perturbations A_rho-v-p=',A_rho,A_v,A_p
    endif

  end subroutine set_params_and_mg_boundary_conds

  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixI^L, ixO^L, w, x)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    logical, save:: first=.true.

    ! Set initial uniform values 
    w(ixI^S, rho_)  = rho0_norm
    w(ixI^S, mom(:))= 0.d0
    w(ixI^S, e_)    = eg0_norm
    w(ixI^S, r_e)   = Er0_norm

  if(mype==0.and.first)then
    print *,'===IN INITIAL CONDITIONS========================================='
    write(*,*) 'converted to normalized values'
    write(*,*) 'normalized density, pressure, energy density gas-rad =',rho0_norm,p0_norm,eg0_norm,Er0_norm
    print *,'===========================================-====================='
    print *,'========GLOBAL values==========='
    print *,'rho0=',rho0
    print *,'p0=',p0
    print *,'eg0=',eg0
    print *,'========GLOBAL values==========='
    first=.false.
  endif

  end subroutine initial_conditions

  subroutine Initialize_Wave(level,qt,qdt,ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L,ixO^L,level
    double precision, intent(in)    :: qt, qdt
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)

    double precision :: vel(ixI^S),pres(ixI^S)

    where (x(ixI^S,1) .lt. wvlength)
      vel(ixI^S)       = A_v*dsin(wavenumber*x(ixI^S,1)-omega*qt)
      pres(ixI^S)      = p0_norm+A_p*dsin(wavenumber*x(ixI^S,1)-omega*qt)
      w(ixI^S, rho_)   = rho0_norm + A_rho*dsin(wavenumber*x(ixI^S,1)-omega*qt)
      {^IFTWOD
      w(ixI^S, mom(2)) = 0.d0
      }
      {^IFTHREED
      w(ixI^S, mom(3)) = 0.d0
      }
      w(ixI^S, mom(1)) = w(ixI^S,rho_)*vel(ixI^S)
      w(ixI^S, e_)     = pres(ixI^S)/(eos%gamma-1.0d0)+half*w(ixI^S,rho_)*vel(ixI^S)**2
      w(ixI^S, r_e)    = arad_norm*(pres(ixI^S)/(w(ixI^S,rho_)*RR))**4
    endwhere

  end subroutine Initialize_Wave

  subroutine set_output_vars(ixI^L,ixO^L,qt,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(in)    :: qt,x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: Trad(ixI^S),Tgas(ixI^S),pth(ixI^S)

    call eos%get_thermal_pressure(w,x,ixI^L,ixO^L,pth)
    call hd_get_trad(w,x,ixI^L,ixO^L,Trad)
    call eos%get_temperature_from_etot(w,x,ixI^L,ixO^L,Tgas)
    w(ixO^S,Tgas_)=Tgas(ixO^S)
    w(ixO^S,Trad_)=Trad(ixO^S)
    w(ixO^S,pres_)=pth(ixO^S)
    w(ixO^S,vel_)=w(ixO^S,mom(1))/w(ixO^S,rho_)
    ! output the AMR level (assuming uniform grid blocks)
    w(ixO^S,amr_)=dlog(((xprobmax1-xprobmin1)/domain_nx1)/dxlevel(1))/dlog(2.0d0)+1.0d0

  end subroutine set_output_vars

end module mod_usr
