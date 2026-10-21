module mod_usr

  ! Include a physics module
  use mod_hd
  use mod_fld
  use mod_eos, only: eos

  implicit none

    double precision :: rho0,v0,T0,T1,wdth
    double precision :: rho0_norm,T0_norm,T1_norm,v0_norm

  ! Storing additional var in the dat file
  integer :: Tgas_,Trad_,pres_,vel_,amr_

contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()

    ! Note how we here must set three values that in turn define M-L-T
    unit_density       =1.2d0        ! 1.2 g cm^-3
    unit_temperature   =1.d7         ! 10^7 K
    unit_length        =24.d0        ! 24 cm

    call usr_params_read(par_files)

    ! A routine for initial conditions is always required
    usr_init_one_grid => initial_conditions

   ! to add selected variables to the .dat file
    usr_modify_output => set_output_vars

    call set_coordinate_system("Cartesian_1D")
    ! Active the physics module
    call hd_activate()

    ! to add selected variables to the .dat file
    Tgas_ = var_set_extravar("Tgas", "Tgas")
    Trad_ = var_set_extravar("Trad", "Trad")
    pres_ = var_set_extravar("pres", "pres")
    vel_ = var_set_extravar("vel", "vel")
    amr_ = var_set_extravar("level", "level")

  end subroutine usr_init

  subroutine usr_params_read(files)
  use mod_global_parameters, only: unitpar
  character(len=*), intent(in) :: files(:)
  integer                      :: n

  namelist /test_list/ rho0, v0, T0, T1, wdth

  do n = 1, size(files)
     open(unitpar, file=trim(files(n)), status="old")
     read(unitpar, test_list, end=111)
     111    close(unitpar)
  end do

  if(mype==0)then
    print *,'============================================'
    write(*,*) 'INPUT GIVEN IN cgs units is'
    write(*,*) 'input density, velocity, temperatures0-1 =',rho0,v0,T0,T1
    write(*,*) 'input dimensionless (length) factor is   =',wdth
    print *,'============================================'
  endif

  end subroutine usr_params_read

  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixI^L, ixO^L, w, x)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: temp(ixI^S)
    double precision :: a,b,Xfrac,Yfrac
    logical, save:: first=.true.

    rho0_norm=rho0/unit_density
    T0_norm=T0/unit_temperature
    T1_norm=T1/unit_temperature
    v0_norm = v0/unit_velocity

    temp(ixI^S) = T0_norm + (T1_norm-T0_norm)*dexp(-x(ixI^S,1)**2.d0/(2*(wdth**2)))

    w(ixI^S,rho_) = rho0_norm*T0_norm/temp(ixI^S) &
                   + arad_norm * (T0_norm**4.d0/temp(ixI^S) - temp(ixI^S)**3.d0)
    w(ixI^S,mom(:)) = 0.d0
    w(ixI^S,mom(1)) = v0_norm
    w(ixI^S,p_) = temp(ixI^S)*w(ixI^S,rho_)
    w(ixI^S,r_e) = arad_norm*(temp(ixI^S)**4.d0)

    call eos%to_conserved(ixI^L,ixI^L,w,x)

  if(mype==0.and.first)then
    print *,'===IN INITIAL CONDITIONS========================================='
    write(*,*) 'converted to normalized values'
    write(*,*) 'normalized density, velocity, temperatures0-1 =',rho0_norm,v0_norm,T0_norm,T1_norm
    print *,'=======================---------------------====================='
    print *,'========GLOBAL values==========='
    print *,'rho0=',rho0
    print *,'T0,T1=',T0,T1
    print *,'v0=',v0
    print *,'========GLOBAL values==========='
    first=.false.
  endif

  end subroutine initial_conditions

  subroutine set_output_vars(ixI^L,ixO^L,qt,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(in)    :: qt,x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: Trad(ixI^S),Tgas(ixI^S),pth(ixI^S)

    call hd_get_pthermal(w,x,ixI^L,ixO^L,pth)
    call hd_get_trad(w,x,ixI^L,ixO^L,Trad)
    call hd_get_temperature_from_etot(w,x,ixI^L,ixO^L,Tgas)
    w(ixO^S,Tgas_)=Tgas(ixO^S)
    w(ixO^S,Trad_)=Trad(ixO^S)
    w(ixO^S,pres_)=pth(ixO^S)
    w(ixO^S,vel_)=w(ixO^S,mom(1))/w(ixO^S,rho_)
    ! output the AMR level (assuming uniform grid blocks)
    w(ixO^S,amr_)=dlog(((xprobmax1-xprobmin1)/domain_nx1)/dxlevel(1))/dlog(2.0d0)+1.0d0

  end subroutine set_output_vars


end module mod_usr
