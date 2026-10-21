!=============================================================================
! 2D HD setup for thermal instability
! Setup has no gravity, may have conduction, has radiative cooling and heating
!=============================================================================

module mod_usr
  use mod_hd
  use mod_eos, only: eos
  implicit none
  double precision              :: rho0, p0, A0

  ! Storing additional variables in the dat file
  integer :: temp_, velx_, vely_, press_, thcond_, rad_, heat_, netheat_, amr_

contains

  subroutine usr_init()
    unit_length        = 1.d9 ! in cm (1 Mm in cm) 
    unit_temperature   = 1.d6 ! in K (1Mk)
    unit_numberdensity = 1.d9 ! in cm^-3
    
    call usr_params_read(par_files)

    usr_init_one_grid   => initonegrid_usr 
    usr_source          => special_source  

    ! to add selected variables to the .dat file
    usr_modify_output => set_output_vars

    call set_coordinate_system("Cartesian_2D")
    call hd_activate()

    ! to add selected variables to the .dat file
    temp_ = var_set_extravar("Temp", "Temp")
    velx_ = var_set_extravar("velx", "velx")
    vely_ = var_set_extravar("vely", "vely")
    press_ = var_set_extravar("pr", "pr")
    thcond_ = var_set_extravar("thcond", "thcond")
    rad_ = var_set_extravar("rad", "rad")
    heat_ = var_set_extravar("heat", "heat")
    netheat_ = var_set_extravar("netheat", "netheat")
    amr_ = var_set_extravar("level", "level")

  end subroutine usr_init

  !> Read this module's parameters from a file
  subroutine usr_params_read(files)
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /usr_list/ A0, rho0, p0

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, usr_list, end=111)
111    close(unitpar)
    end do

    if(rho0*p0<smalldouble) call mpistop("Wrong input values: rho0-p0 must be positive")

  end subroutine usr_params_read

   subroutine initonegrid_usr(ixG^L,ix^L,w,x)    
    
    integer, intent(in) :: ixG^L, ix^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    double precision :: factorrho(ixG^S)
    double precision :: k1,k2

    logical, save :: first=.true.
   

    k1 = 2.d0 * dpi / (xprobmax1 - xprobmin1) 
    k2 = 2.d0 * dpi / (xprobmax2 - xprobmin2) 
    
    ! 1. Density: Perturbed in a cosine wave
    factorrho(ix^S)=1.0d0-A0*dcos(k1*x(ix^S,1))*dcos(k2*x(ix^S,2))
    if(any(factorrho(ix^S)<smalldouble)) call mpistop("Change A0: negative initial density")
    w(ix^S,rho_)=rho0*factorrho(ix^S)

    ! 2. Velocity: no perturbation (since entropy mode)
    w(ix^S,mom(1))= 0.0d0
    w(ix^S,mom(2))= 0.0d0
       
    ! 3. Pressure : isobaric (entropy mode)
    w(ix^S,p_)=p0

    if(mype==0.and.first)then
        write(*,*)'Uniform medium with pressure p0=',p0
        write(*,*)'Uniform medium with density rho0=',rho0
        write(*,*)'Uniform medium with temperature=',p0/rho0
        write(*,*)'Uniform medium with perturbation amplitude=',A0
        write(*,*)'perturbing density only'
        write(*,*)'k1 and k2  is=',k1,k2
        write(*,*)'domain=',xprobmin1,xprobmax1,xprobmin2,xprobmax2
        first=.false.
    endif
    call eos%to_conserved(ixG^L,ix^L,w,x)
  
  end subroutine initonegrid_usr

  subroutine special_source(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)
     use mod_radiative_cooling, only: getvar_cooling

    integer, intent(in) :: ixI^L, ixO^L, iw^LIM
    double precision, intent(in) :: qdt, qtC, qt
    double precision, intent(in) :: x(ixI^S,1:ndim), wCT(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: bQgrid(ixI^S),winit_nopert(ixI^S,1:nw)
    
    if(hd_radiative_cooling)then
       winit_nopert(ixI^S,rho_)= rho0 
       winit_nopert(ixI^S,mom(1))=0.0d0
       winit_nopert(ixI^S,mom(2))=0.0d0
       winit_nopert(ixI^S,e_)=p0/(eos%gamma - 1.0d0)   
       ! this is adding heating that balances the cooling from the uniform background
       call getvar_cooling(ixI^L,ixO^L,winit_nopert,x,bQgrid,rc_fl)
       w(ixO^S,e_)=w(ixO^S,e_)+qdt*bQgrid(ixO^S)
    endif
    
  end subroutine special_source

  subroutine set_output_vars(ixI^L,ixO^L,qt,w,x)
    use mod_radiative_cooling
    use mod_thermal_conduction
    use mod_global_parameters
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(in)    :: qt,x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    !!double precision :: local_dt
    double precision :: wlocal(ixI^S,1:nw),wcond(ixI^S,1:nw)
    double precision :: bQgrid(ixI^S),bQgridB(ixI^S),winit_nopert(ixI^S,1:nw)
    
    wlocal(ixI^S,1:nw)=w(ixI^S,1:nw)
    
    ! output thermal conduction TC
    wcond(ixI^S,1:nw)=0.0d0
    if(hd_thermal_conduction)then
      ! 5 last arguments unused here: fix_conserve_at_step to false, my_dt=1, igrid/nflux to 1
      !!!convert wlocal e to ei here!!!
      call hd_e_to_ei(ixI^L,ixI^L,wlocal,x)
      call sts_set_source_tc_hd(ixI^L,ixO^L,wlocal,x,wcond,.false.,1.d0,1,1,tc_fl)
      call hd_ei_to_e(ixI^L,ixI^L,wlocal,x)
      w(ixO^S,thcond_)=wcond(ixO^S,tc_fl%e_)
    endif
    ! store the cooling and heating balance
    if(hd_radiative_cooling)then
       winit_nopert(ixI^S,rho_)= rho0 
       winit_nopert(ixI^S,mom(1))=0.0d0
       winit_nopert(ixI^S,mom(2))=0.0d0
       winit_nopert(ixI^S,e_)=p0/(eos%gamma - 1.0d0)   
       call getvar_cooling(ixI^L,ixO^L,winit_nopert,x,bQgrid,rc_fl)
       w(ixO^S,heat_)=bQgrid(ixO^S)
       !!using the global timestep dt here, could be unset at last timesave!!!
       !!local_dt=max(dt,1.0d-6)
       !!call getvar_cooling_exact(local_dt,ixI^L,ixO^L,wlocal,wlocal,x,bQgridB,rc_fl)
       call getvar_cooling(ixI^L,ixO^L,wlocal,x,bQgridB,rc_fl)
       w(ixO^S,rad_)=bQgridB(ixO^S)
       w(ixO^S,netheat_)=(bQgrid(ixO^S)-bQgridB(ixO^S))
    endif  

    ! now also add extra primitive variables in output dat file 
    call eos%to_primitive(ixI^L,ixO^L,wlocal,x)
    w(ixO^S,temp_)=wlocal(ixO^S,p_)/wlocal(ixO^S,rho_)
    w(ixO^S,press_)=wlocal(ixO^S,p_)
    w(ixO^S,velx_)=wlocal(ixO^S,mom(1))
    w(ixO^S,vely_)=wlocal(ixO^S,mom(2))
    ! output the AMR level (assuming uniform grid blocks)
    w(ixO^S,amr_)=dlog(((xprobmax1-xprobmin1)/domain_nx1)/dxlevel(1))/dlog(2.0d0)+1.0d0

  end subroutine set_output_vars
end module mod_usr
