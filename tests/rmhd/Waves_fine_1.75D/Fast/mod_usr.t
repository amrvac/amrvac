module mod_usr
  use mod_mhd 
  use mod_eos, only: eos

  implicit none

  ! input variables with units
  double precision :: rho0,p0,Bx0,By0,Bz0
  ! input variables that are dimensionless
  double precision :: rho_kappa_lambda,ampl

  ! derived vars
  double precision :: omega,wavenumber,csound_norm,wvlength
  double precision :: A_rho,A_vx, A_p, A_vy, A_vz, A_By, A_Bz, A_e, A_Er
  double precision :: eintg0,rho0_norm,eg0_norm,Er0_norm,p0_norm,T0_norm
  double precision :: Bx0_n,By0_n,Bz0_n
  double precision :: p0_by_Er0, p0_by_p0_B, p0_B
  double precision :: va, ca, vax, vax_norm, cslow_n, cfast_n, va_norm

  ! Storing additional var in the dat file
  integer :: Tgas_,Trad_,pres_,velx_,amr_

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

    ! Choose coordinate system 
    call set_coordinate_system("Cartesian_1.75D")

    ! Active the physics module
    call mhd_activate() 

    ! to add selected variables to the .dat file
    Tgas_ = var_set_extravar("Tgas", "Tgas")
    Trad_ = var_set_extravar("Trad", "Trad")
    pres_ = var_set_extravar("pres", "pres")
    velx_ = var_set_extravar("velx", "velx")
    amr_ = var_set_extravar("level", "level")

  end subroutine usr_init

  subroutine set_params_and_mg_boundary_conds()
    use mod_global_parameters
    use mod_fld

    ! here we normalize all input values and compute the equilibrium parameters
    rho0_norm = rho0/unit_density
    p0_norm   = p0/unit_pressure
    T0_norm   = p0_norm/(rho0_norm*RR)
    Er0_norm  = arad_norm*T0_norm**4.d0
    Bx0_n = Bx0/unit_magneticfield
    By0_n = By0/unit_magneticfield
    Bz0_n = Bz0/unit_magneticfield

    eg0_norm =p0_norm/(eos%gamma-1.0d0)+half*(Bx0_n**2+By0_n**2+Bz0_n**2) 

    csound_norm=dsqrt(eos%gamma*p0/rho0)/unit_velocity
    vax_norm=dsqrt(Bx0_n**2/rho0_norm)

    va_norm = dsqrt((Bx0_n**2 + By0_n**2+ Bz0_n**2)/rho0_norm)

    va=va_norm
    ca=csound_norm
    vax=vax_norm
    cslow_n = dsqrt(0.5d0*((va*va + ca*ca) - dsqrt((va*va + ca*ca)*(va*va + ca*ca) - 4.0d0*ca*ca*vax*vax)))
    cfast_n = dsqrt(0.5d0*((va*va + ca*ca) + dsqrt((va*va + ca*ca)*(va*va + ca*ca) - 4.0d0*ca*ca*vax*vax)))

    p0_B = half*(Bx0_n**2+By0_n**2+Bz0_n**2)
    p0_by_Er0 = p0_norm/Er0_norm 
    p0_by_p0_B = p0_norm/p0_B 

    if(mype==0)then
      write(*,*) 'derived dimensionless rad/gas energy ratio =',Er0_norm/eg0_norm
      write(*,*) 'derived ratio r=E/(4gamma e) =',Er0_norm/(4.d0*eos%gamma*eg0_norm)
      write(*,*) 'derived Boltzmann ratio Bo=4 gamma cs e/(cE) =',4.0d0*eos%gamma*csound_norm*eg0_norm/(c_norm*Er0_norm)
      print*, 'derived ratio p0_by_Er0', p0_by_Er0
      print*, 'derived ratio (inverse plasma beta) p0_by_p0_B', p0_by_p0_B
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
    omega=wavenumber*cfast_n

    A_rho = ampl
    A_vx = omega/(wavenumber*rho0_norm)*A_rho
    A_p = (ca**2)*A_rho
    A_By = By0_n*(ampl/rho0_norm)/(1.0d0 - (vax**2)*(wavenumber**2)/(omega**2))
    A_Bz = Bz0_n*(ampl/rho0_norm)/(1.0d0 - (vax**2)*(wavenumber**2)/(omega**2))
    if (Bx0_n .eq. 0) then
      A_vy = 0.d0 
      A_vz = 0.d0
    else
      A_vy = (By0_n/Bx0_n)*(omega/wavenumber)*(ampl/rho0_norm)/(1.0d0 - (omega**2)/((vax**2)*(wavenumber**2)))
      A_vz = (Bz0_n/Bx0_n)*(omega/wavenumber)*(ampl/rho0_norm)/(1.0d0 - (omega**2)/((vax**2)*(wavenumber**2)))
    endif
    A_Er = 0.0d0 
    A_e = A_p/(eos%gamma-one) + By0_n*A_By + Bz0_n*A_Bz
  
    if(mype==0)then
    print *,'wavelength in physical units is=',wvlength*unit_length
    write(*,*) 'derived dimensionless wavelength for wave =',wvlength
    write(*,*) 'derived dimensionless omega and k for wave =',omega,wavenumber
    write(*,*) 'normalized sound speed is   =',csound_norm
    write(*,*) 'normalized alfven x speed is=',vax_norm
    write(*,*) 'normalized alfven   speed is=',va_norm
    write(*,*) 'normalized slow     speed is=',cslow_n
    write(*,*) 'normalized fast     speed is=',cfast_n
    write(*,*) 'check on input and derivation based on ',rho_kappa_lambda,rho0,unit_density,fld_kappa0
    write(*,*) 'dispersion relation asks for omega=',omega,' to equal k= ',wavenumber,' times vax=',vax_norm
    write(*,*) 'and sets amplitude perturbations A_rho-A_p=',A_rho,A_p
    write(*,*) 'and sets amplitude perturbations A_vx,A_vy,A_vz=',A_vx,A_vy,A_vz
    write(*,*) 'and sets amplitude perturbations A_By,A_Bz=',A_By,A_Bz
    write(*,*) 'and sets amplitude perturbations A_e,A_Er=',A_e,A_Er
    endif

  end subroutine set_params_and_mg_boundary_conds
      
  !> Read parameters from a file
  subroutine usr_params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /wave_list/ rho0, p0, Bx0, By0, Bz0, rho_kappa_lambda, ampl

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       rewind(unitpar)
       read(unitpar, wave_list, end=113)
113    close(unitpar)
    end do

    eintg0 = p0/(eos%gamma - one)
    if(mype==0)then
    print *,'============================================'
    write(*,*) 'INPUT GIVEN IN cgs units is'
    write(*,*) 'input density, gas pressure =',rho0,p0
    write(*,*) 'input magnetic field values =',Bx0,By0,Bz0
    write(*,*) 'derived internal gas energy density is =',eintg0
    write(*,*) 'input dimensionless amplitude is   =',ampl
    write(*,*) 'input dimensionless parameter for wavelength   =',rho_kappa_lambda
    print *,'============================================'
    endif

  end subroutine usr_params_read

  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixI^L, ixO^L, w, x)
    use mod_global_parameters
    use mod_fld

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    logical, save:: first=.true.

    ! Set initial uniform values
    w(ixI^S, rho_)  = rho0_norm
    w(ixI^S, mom(:))= 0.d0
    w(ixI^S, e_)    = eg0_norm
    w(ixI^S, r_e)   = Er0_norm
    w(ixI^S,mag(1)) = Bx0_n
    w(ixI^S,mag(2)) = By0_n
    w(ixI^S,mag(3)) = Bz0_n

  if(mype==0.and.first)then
    print *,'===IN INITIAL CONDITIONS========================================='
    write(*,*) 'converted to normalized values'
    write(*,*) 'normalized density, pressure, energy density gas-rad =',rho0_norm,p0_norm,eg0_norm,Er0_norm
    write(*,*) 'normalized temperature gas-rad                       =',T0_norm
    write(*,*) 'normalized magnetic field values                     =',Bx0_n,By0_n,Bz0_n
    print *,'================================================================='
    print *,'========GLOBAL values==========='
    print *,'rho0=',rho0
    print *,'p0=',p0
    print *,'eintg0=',eintg0
    write(*,*) ' magnetic field =',Bx0,By0,Bz0
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

    double precision :: velx(ixI^S),pres(ixI^S)
    double precision :: vely(ixI^S),velz(ixI^S)

    where (x(ixI^S,1) .lt. wvlength)
      velx(ixI^S)       = A_vx*dsin(wavenumber*x(ixI^S,1)-omega*qt)
      vely(ixI^S)       = A_vy*dsin(wavenumber*x(ixI^S,1)-omega*qt)
      velz(ixI^S)       = A_vz*dsin(wavenumber*x(ixI^S,1)-omega*qt)
      pres(ixI^S)      = p0_norm+A_p*dsin(wavenumber*x(ixI^S,1)-omega*qt)
      w(ixI^S, rho_)   = rho0_norm + A_rho*dsin(wavenumber*x(ixI^S,1)-omega*qt)
      w(ixI^S, mom(1)) = w(ixI^S,rho_)*velx(ixI^S)
      w(ixI^S, mom(2)) = w(ixI^S,rho_)*vely(ixI^S)
      w(ixI^S, mom(3)) = w(ixI^S,rho_)*velz(ixI^S)
      w(ixI^S, r_e)    = arad_norm*(pres(ixI^S)/(w(ixI^S,rho_)*RR))**4
      !!w(ixI^S, r_e)    = Er0_norm+A_Er*dsin(wavenumber*x(ixI^S,1)-omega*qt)
      w(ixI^S, mag(1)) = Bx0_n
      w(ixI^S, mag(2)) = By0_n + A_By*dsin(wavenumber*x(ixI^S,1)-omega*qt)
      w(ixI^S, mag(3)) = Bz0_n + A_Bz*dsin(wavenumber*x(ixI^S,1)-omega*qt)
      w(ixI^S, e_)     = pres(ixI^S)/(eos%gamma-1.0d0)+half*w(ixI^S,rho_)* &
                               (velx(ixI^S)**2+vely(ixI^S)**2+velz(ixI^S)**2) &
         +half*(Bx0_n*Bx0_n+w(ixI^S, mag(2))*w(ixI^S, mag(2)) + w(ixI^S, mag(3))*w(ixI^S, mag(3)))
      !!w(ixI^S, e_)    = eg0_norm+A_e*dsin(wavenumber*x(ixI^S,1)-omega*qt)
    endwhere

  end subroutine Initialize_Wave

  subroutine set_output_vars(ixI^L,ixO^L,qt,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(in)    :: qt,x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: Trad(ixI^S),Tgas(ixI^S),pth(ixI^S)

    call eos%get_thermal_pressure(w,x,ixI^L,ixO^L,pth)
    call mhd_get_trad(w,x,ixI^L,ixO^L,Trad)
    call eos%get_temperature_from_etot(w,x,ixI^L,ixO^L,Tgas)
    w(ixO^S,Tgas_)=Tgas(ixO^S)
    w(ixO^S,Trad_)=Trad(ixO^S)
    w(ixO^S,pres_)=pth(ixO^S)
    w(ixO^S,velx_)=w(ixO^S,mom(1))/w(ixO^S,rho_)
    ! output the AMR level (assuming uniform grid blocks)
    w(ixO^S,amr_)=dlog(((xprobmax1-xprobmin1)/domain_nx1)/dxlevel(1))/dlog(2.0d0)+1.0d0

  end subroutine set_output_vars

end module mod_usr
