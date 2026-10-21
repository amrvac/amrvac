module mod_usr

  ! Include a physics module
  use mod_hd
  use mod_eos, only: eos
  use mod_fld

  implicit none

  double precision :: e_eq, Tgas_eq
  double precision :: rho0
  double precision :: factor0
  double precision :: E_r0
  ! Storing additional var in the dat file
  integer :: delta_e_, delta_T_, delta_Trel_


contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods
    use mod_constants

    ! Note how we here must set three values that in turn define M-L-T
    ! Because we check on deviation from unity, we trick it here with ALMOST unity
    unit_density       =1.0000001d0  ! 1 g cm^-3
    unit_temperature   =1.0000001d0  ! 1 K
    unit_time          =1.0000001d0  ! 1 sec
    ! better chocie for this problem is
    unit_density       =1.d-7        ! g cm^-3
    unit_temperature   =1.d7         ! 10^7 K
    unit_time          =1.0000001d0  ! 1 sec
 
    call usr_params_read(par_files)

    ! A routine for initial conditions is always required
    usr_init_one_grid => initial_conditions

   ! to add selected variables to the .dat file
    usr_modify_output => set_output_vars

    call set_coordinate_system("Cartesian_2D")

    ! Active the physics module
    call hd_activate()

    ! to add selected variables to the .dat file
    delta_e_ = var_set_extravar("delta_e", "delta_e")
    delta_T_ = var_set_extravar("delta_T", "delta_T")
    delta_Trel_ = var_set_extravar("delta_Trel", "delta_Trel")

  end subroutine usr_init


subroutine usr_params_read(files)
  use mod_global_parameters, only: unitpar
  character(len=*), intent(in) :: files(:)
  integer                      :: n

  namelist /test_list/ rho0, E_r0, factor0

  do n = 1, size(files)
     open(unitpar, file=trim(files(n)), status="old")
     read(unitpar, test_list, end=111)
     111    close(unitpar)
  end do

  if(mype==0)then
    print *,'============================================'
    write(*,*) 'INPUT GIVEN IN cgs units is'
    write(*,*) 'input density, Erad, egas=',rho0,E_r0
    write(*,*) 'input dimensionless factor is=',factor0
    print *,'============================================'
  endif

end subroutine usr_params_read

  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixI^L, ixO^L, w, x)
    use mod_global_parameters
    use mod_constants

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S, ndim)
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision :: a,b,Xfrac,Yfrac,ne_np
    logical, save:: first=.true.

    Tgas_eq = ((E_r0/unit_Erad)/arad_norm)**(1.d0/4.d0)
    e_eq= RR*(rho0/unit_density)*Tgas_eq/(eos%gamma-one)

    ! Set initial values for w
    w(ixI^S, rho_) = rho0/unit_density
    w(ixI^S, mom(:)) = zero
    w(ixI^S,r_e) = E_r0/unit_Erad
    w(ixI^S, e_) = factor0*e_eq


  if(mype==0.and.first)then
    write(*,*)'initial dimensionless density             =',rho0/unit_density
    write(*,*)'initial dimensionless radiation E0        =',E_r0/unit_pressure
    write(*,*)'initial dimensionless gas energy density e=',factor0*e_eq
    write(*,*)'initial ratio to equilibrium density e    =',factor0
    write(*,*)'dimensionless equilibrium e=',e_eq
    print *,'with units e_eq=',e_eq*unit_pressure,' so log10 value is=',dlog10(e_eq*unit_pressure)
    write(*,*)'equilibrium dimensionless Tgas_eq=',Tgas_eq
    print *,'with units Tgas_eq=',Tgas_eq*unit_temperature
    print *,'========GLOBAL values==========='
    print *,'rho0=',rho0
    print *,'E_r0=',E_r0
    print *,'factor0=',factor0
    print *,'fld_kappa0    =',fld_kappa0
    print *,'========GLOBAL values==========='
    first=.false.
  endif

  end subroutine initial_conditions

  subroutine set_output_vars(ixI^L,ixO^L,qt,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(in)    :: qt,x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: Trad(ixI^S),Tgas(ixI^S)

    w(ixO^S,delta_e_)=w(ixO^S,e_)/e_eq
    call hd_get_trad(w,x,ixI^L,ixO^L,Trad)
    call eos%get_temperature_from_etot(w,x,ixI^L,ixO^L,Tgas)
    w(ixO^S,delta_T_)=Tgas(ixO^S)-Trad(ixO^S)
    w(ixO^S,delta_Trel_)=(Tgas(ixO^S)-Trad(ixO^S))/Tgas_eq
    ! output the AMR level (assuming uniform grid blocks)
    !w(ixO^S,amr_)=dlog(((xprobmax1-xprobmin1)/domain_nx1)/dxlevel(1))/dlog(2.0d0)+1.0d0

  end subroutine set_output_vars

end module mod_usr
