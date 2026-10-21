module mod_usr

  ! Include a physics module
  use mod_mhd
  use mod_eos, only: eos
  use mod_fld
  use mod_multigrid_coupling

  implicit none

  ! input values in cgs units
  double precision :: rho1,rho2,vx1,Bx1,T1
  double precision :: vy1, vz1, By1, Bz1 

  ! derived values from RH jumps
  double precision :: vx2,T2,Bx2
  double precision :: vy2, vz2, By2, Bz2 

  ! normalized values
  double precision :: rho1_norm,rho2_norm,vx1_norm,vx2_norm,T1_norm,T2_norm,Bx1_norm,Bx2_norm
  double precision :: vy1_norm, vz1_norm, By1_norm, Bz1_norm, Btotsq1, vdotB1, vsq1, fac1
  double precision :: vy2_norm, vz2_norm, By2_norm, Bz2_norm, Btotsq2, vdotB2, vsq2, fac2
  double precision :: p1_norm,p2_norm,eg1_norm,eg2_norm,Er1_norm,Er2_norm
  double precision :: momc_n, T2A_norm, T2B_norm

  ! for computing the balanced temperature
  double precision :: c1A_norm,c0A_norm,val1,val2
  double precision :: c1B_norm,c0B_norm

  ! Storing additional var in the dat file
  integer :: Tgas_,Trad_,pres_,velx_,vely_,amr_

contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()

    ! Note how we here must set three values that in turn define M-L-T
    unit_density      =0.000116963d0       ! g cm^-3
    unit_velocity     =9.82668d7         !  cm s^-1
    unit_length       =1.d5         ! 10^5 cm

    call usr_params_read(par_files)

    usr_set_parameters => set_params_and_mg_boundary_conds

    usr_init_one_grid => initial_conditions

    ! Manually refine grid near shock
    usr_refine_grid => refine_shock

   ! to add selected variables to the .dat file
    usr_modify_output => set_output_vars

    ! Boundary conditions
    usr_special_bc => boundary_conditions

    ! custom output
    !usr_aux_output      => specialvar_output
    !usr_add_aux_names   => specialvarnames_output

    call set_coordinate_system("Cartesian_1.75D")

    ! Active the physics module
    call mhd_activate() 

    ! to add selected variables to the .dat file
    Tgas_ = var_set_extravar("Tgas", "Tgas")
    Trad_ = var_set_extravar("Trad", "Trad")
    pres_ = var_set_extravar("pres", "pres")
    velx_ = var_set_extravar("velx", "velx")
    vely_ = var_set_extravar("vely", "vely")
    amr_ = var_set_extravar("level", "level")

  end subroutine usr_init

  subroutine set_params_and_mg_boundary_conds
    use mod_global_parameters

    vx2=rho1*vx1/rho2

    rho1_norm = rho1/unit_density
    rho2_norm = rho2/unit_density
    T1_norm = T1/unit_temperature
    p1_norm=T1_norm*rho1_norm*RR
    Er1_norm = arad_norm*T1_norm**4.d0
    vx1_norm = vx1/unit_velocity
    vx2_norm = vx2/unit_velocity
    vy1_norm = vy1/unit_velocity
    vz1_norm = vz1/unit_velocity
    momc_n=rho1*vx1/(unit_density*unit_velocity)

    ! enforce Bx constant
    Bx2=Bx1
    Bx1_norm = Bx1/unit_magneticfield
    Bx2_norm = Bx1_norm
    By1_norm = By1/unit_magneticfield
    Bz1_norm = Bz1/unit_magneticfield

    fac1=Bx1_norm**2-momc_n**2/rho1_norm
    fac2=Bx2_norm**2-momc_n**2/rho2_norm

    if(momc_n.ne.zero)then
      if(fac2==zero)call mpistop('case not allowed yet')
      By2_norm=By1_norm*fac1/fac2
      Bz2_norm=Bz1_norm*fac1/fac2
      vy2_norm=vy1_norm-(Bx1_norm/momc_n)*(By1_norm-By2_norm)
      vz2_norm=vz1_norm-(Bx1_norm/momc_n)*(Bz1_norm-Bz2_norm)
    else
      if(mype==0)print *,'fixing all tangential velocity and magnetic equal'
      vy2_norm=vy1_norm
      vz2_norm=vz1_norm
      By2_norm=By1_norm
      Bz2_norm=Bz1_norm
    endif

    vy2=vy2_norm*unit_velocity
    vz2=vz2_norm*unit_velocity
    By2=By2_norm*unit_magneticfield
    Bz2=Bz2_norm*unit_magneticfield

    Btotsq1=Bx1_norm**2+By1_norm**2+Bz1_norm**2
    Btotsq2=Bx2_norm**2+By2_norm**2+Bz2_norm**2

    vsq1=vx1_norm**2+vy1_norm**2+vz1_norm**2
    vsq2=vx2_norm**2+vy2_norm**2+vz2_norm**2

    vdotB1=vx1_norm*Bx1_norm+vy1_norm*By1_norm+vz1_norm*Bz1_norm
    vdotB2=vx2_norm*Bx2_norm+vy2_norm*By2_norm+vz2_norm*Bz2_norm

    if(momc_n.eq.zero)then
    ! this computes T2 ensuring that the momentum flux is exactly equal
    ! use the Halley root finder for the 4th order polynomial in T2 (normalized)
    c1A_norm=3.0d0*RR*rho2_norm/arad_norm
    c0A_norm=((rho2_norm-rho1_norm)*momc_n**2/(rho2_norm*rho1_norm) &
             +rho1_norm*T1_norm*RR+Er1_norm/3.0d0+half*(Btotsq1-Btotsq2))*3.0d0/arad_norm
    T2A_norm=T1_norm ! this is a bad guess
    call Halley_method(T2A_norm,c0A_norm,c1A_norm)
    if(mype==0)print *,'Momentum balance needs T2 (normalized)=',T2A_norm
    T2_norm=T2A_norm
    else
    ! this computes T2 ensuring that the energy flux is exactly equal
    ! use the Halley root finder for the 4th order polynomial in T2 (normalized)
    c1B_norm=3.0d0*RR*rho2_norm*eos%gamma/(arad_norm*4.0d0*(eos%gamma-1.d0))
    c0B_norm=(3.0d0/(4.0d0*arad_norm))* &
       ((eos%gamma*RR*rho1_norm*T1_norm/(eos%gamma-1.d0) &
         +4.0d0*Er1_norm/3.0d0+half*rho1_norm*vsq1+Btotsq1)*rho2_norm/rho1_norm &
       -half*rho2_norm*vsq2-Btotsq2-(rho2_norm/momc_n)*Bx1_norm*(vdotB1-vdotB2))
    T2B_norm=T1_norm ! this is a bad guess
    call Halley_method(T2B_norm,c0B_norm,c1B_norm)
    if(mype==0)print *,'Energy balance needs T2 (normalized)=',T2B_norm
    T2_norm=T2B_norm
    endif

    T2=T2_norm*unit_temperature
    p2_norm=T2_norm*rho2_norm*RR
    Er2_norm = arad_norm*T2_norm**4.d0
    eg1_norm=p1_norm/(eos%gamma-1.0d0)+half*rho1_norm*vsq1+half*Btotsq1
    eg2_norm=p2_norm/(eos%gamma-1.0d0)+half*rho2_norm*vsq2+half*Btotsq2

    if(mype==0)then
     print*, 'M_1: ', vx1_norm/dsqrt(eos%gamma*p1_norm/rho1_norm) 
     print*, 'M_2: ', vx2_norm/dsqrt(eos%gamma*p2_norm/rho2_norm) 
     print*, 'RMHD-quantities  : ', 'Left', ' | ', 'Right'
     print*, 'density          :', rho1_norm, ' | ', rho2_norm
     print*, 'velocity         :', vx1_norm,vy1_norm,vz2_norm, ' | ', vx2_norm,vy2_norm,vz2_norm 
     print*, 'momentum x       :', rho1_norm*vx1_norm, ' | ', rho2_norm*vx2_norm
     print*, 'gas pressure     :', p1_norm, ' | ', p2_norm
     print*, 'gas energy       :', eg1_norm, ' | ', eg2_norm
     print*, 'radiation energy :', Er1_norm, ' | ', Er2_norm
     print*, 'magnetic field   :', Bx1_norm,By1_norm,Bz1_norm, ' | ', Bx2_norm,By2_norm,Bz2_norm
     print*, 'RMHD-fluxes: ', 'Left', ' | ', 'Right'
     print*, 'density flux   =', rho1_norm*vx1_norm, ' | ', rho2_norm*vx2_norm
     print*, 'momentum flux  =', (rho1_norm*vx1_norm*vx1_norm+p1_norm+Er1_norm/3+half*Btotsq1), ' | ',&
                                 (rho2_norm*vx2_norm*vx2_norm+p2_norm+Er2_norm/3+half*Btotsq2) 
     print*, 'gas energy flux=', p1_norm*vx1_norm + eg1_norm*vx1_norm, ' | ', &
                                 p2_norm*vx2_norm + eg2_norm*vx2_norm
     print*, 'rad energy flux=', Er1_norm*vx1_norm, ' | ', Er2_norm*vx2_norm
     print*, 'total energy   =', (p1_norm+eg1_norm+Er1_norm)*vx1_norm, ' | ', &
                                 (p2_norm+eg2_norm+Er2_norm)*vx2_norm
     print*, 'energy flux    =', ((p1_norm+eg1_norm+4.0*Er1_norm/3.0+half*Btotsq1)*vx1_norm &
                                 - Bx1_norm*(vx1_norm*Bx1_norm+vy1_norm*By1_norm+vz1_norm*Bz1_norm)), ' | ',&
                                 ((p2_norm+eg2_norm+4.0*Er2_norm/3.0+half*Btotsq2)*vx2_norm &
                                 - Bx2_norm*(vx2_norm*Bx2_norm+vy2_norm*By2_norm+vz2_norm*Bz2_norm)) 
     print*, 'momentum_y', (rho1_norm*vx1_norm*vy1_norm - Bx1_norm*By1_norm), ' | ', &
                           (rho2_norm*vx2_norm*vy2_norm - Bx2_norm*By2_norm) 
     print*, 'momentum_z', (rho1_norm*vx1_norm*vz1_norm - Bx1_norm*Bz1_norm), ' | ', &
                           (rho2_norm*vx2_norm*vz2_norm - Bx2_norm*Bz2_norm) 
     print*, 'magfield_y', (By1_norm*vx1_norm - Bx1_norm*vy1_norm), ' | ', &
                           (By2_norm*vx2_norm - Bx2_norm*vy2_norm) 
     print*, 'magfield_z', (Bz1_norm*vx1_norm - Bx1_norm*vz1_norm), ' | ', &
                           (Bz2_norm*vx2_norm - Bx2_norm*vz2_norm) 
    endif

    mg%bc(1, mg_iphi)%bc_type = mg_bc_dirichlet
    mg%bc(1, mg_iphi)%bc_value = Er1_norm
    !mg%bc(2, mg_iphi)%bc_type = mg_bc_dirichlet
    !mg%bc(2, mg_iphi)%bc_value = Er2_norm
    mg%bc(2, mg_iphi)%bc_type = mg_bc_neumann
    mg%bc(2, mg_iphi)%bc_value = 0.0d0

  end subroutine set_params_and_mg_boundary_conds

  !> Read parameters from a file
  subroutine usr_params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /shock_list/ rho1, rho2, vx1, T1, Bx1, vy1, vz1, By1, Bz1

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       rewind(unitpar)
       read(unitpar, shock_list, end=113)
113    close(unitpar)
    end do

  if(mype==0)then
    print *,'============================================'
    write(*,*) 'INPUT GIVEN IN cgs units is'
    write(*,*) 'input density 1-2     =',rho1,rho2
    write(*,*) 'input Temperature 1 =',T1
    write(*,*) 'input velocity 1      =',vx1,vy1,vz1
    write(*,*) 'input B field 1      =',Bx1,By1,Bz1
    print *,'============================================'
  endif

  end subroutine usr_params_read

  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixI^L, ixO^L, w, x)

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    logical, save:: first=.true.

    w(ixI^S,rho_)   = rho1_norm
    w(ixI^S,mom(1)) = rho1_norm*vx1_norm
    w(ixI^S,mom(2)) = rho1_norm*vy1_norm
    w(ixI^S,mom(3)) = rho1_norm*vz1_norm
    w(ixI^S,e_)     = eg1_norm
    w(ixI^S,r_e)    = Er1_norm
    w(ixI^S,mag(1)) = Bx1_norm
    w(ixI^S,mag(2)) = By1_norm 
    w(ixI^S,mag(3)) = Bz1_norm

    where (x(ixI^S,1) .gt. 0.d0)
      w(ixI^S,rho_)   = rho2_norm
      w(ixI^S,mom(1)) = rho2_norm*vx2_norm
      w(ixI^S,mom(2)) = rho2_norm*vy2_norm
      w(ixI^S,mom(3)) = rho2_norm*vz2_norm
      w(ixI^S,e_)     = eg2_norm
      w(ixI^S,r_e)    = Er2_norm
      w(ixI^S,mag(1)) = Bx2_norm
      w(ixI^S,mag(2)) = By2_norm
      w(ixI^S,mag(3)) = Bz2_norm 
    end where


  if(mype==0.and.first)then
    print *,'===IN INITIAL CONDITIONS========================================='
    write(*,*) 'converted to normalized values'
    write(*,*) 'normalized densities          =',rho1_norm,rho2_norm
    write(*,*) 'normalized velocities x       =',vx1_norm,vx2_norm
    write(*,*) 'normalized velocities y       =',vy1_norm,vy2_norm
    write(*,*) 'normalized velocities z       =',vz1_norm,vz2_norm
    write(*,*) 'normalized temperatures       =',T1_norm,T2_norm
    write(*,*) 'normalized pressures          =',p1_norm,p2_norm
    write(*,*) 'normalized gas energies       =',eg1_norm,eg2_norm
    write(*,*) 'normalized radiation energies =',Er1_norm,Er2_norm
    write(*,*) 'normalized B x                =',Bx1_norm,Bx2_norm
    write(*,*) 'normalized B y                =',By1_norm,By2_norm
    write(*,*) 'normalized B z                =',Bz1_norm,Bz2_norm
    print *,'================================================================='
    print *,'========GLOBAL values==========='
    print *,'rho 1-2=',rho1,rho2
    print *,'T   1-2=',T1,T2
    print *,'v   1-2=',vx1,vx2,vy1,vy2,vz1,vz2
    print *,'B   1-2=',Bx1,Bx2,By1,By2,Bz1,Bz2
    print *,'================================================================='
    print*, 'RMHD-fluxes for stationary shock: ', 'Left', ' | ', 'Right'
    print*, 'density flux must be exactly equal  :', rho1_norm*vx1_norm, ' | ', rho2_norm*vx2_norm
    print*, 'momentum flux must be exact   equal :',rho1_norm*vx1_norm**2+p1_norm+Er1_norm/3+half*Btotsq1, &
                                              ' | ',rho2_norm*vx2_norm**2+p2_norm+Er2_norm/3+half*Btotsq2
   val1= (p1_norm+eg1_norm+4.0d0*Er1_norm/3.0d0+half*Btotsq1)*vx1_norm-Bx1_norm*vdotB1
    val2= (p2_norm+eg2_norm+4.0d0*Er2_norm/3.0d0+half*Btotsq2)*vx2_norm-Bx2_norm*vdotB2
    print*, 'energy flux must be equal           :', val1, ' | ', val2 
    print*, 'with radiation and gas T equal on each side: LEFT  is', T1_norm,(Er1_norm/arad_norm)**0.25d0
    print*, 'with radiation and gas T equal on each side: RIGHT is', T2_norm,(Er2_norm/arad_norm)**0.25d0
    print *,'================================================================='
    first=.false.
  endif

  end subroutine initial_conditions


  subroutine boundary_conditions(qdt,qt,ixI^L,ixB^L,iB,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixB^L, iB
    double precision, intent(in)    :: qdt, qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    select case (iB)
    case(1)
      w(ixB^S,rho_)   = rho1_norm
      w(ixB^S,mom(1)) = rho1_norm*vx1_norm
      w(ixB^S,mom(2)) = rho1_norm*vy1_norm
      w(ixB^S,mom(3)) = rho1_norm*vz1_norm 
      w(ixB^S,e_)     = eg1_norm
      w(ixB^S,r_e)    = Er1_norm
      w(ixB^S,mag(1)) = Bx1_norm
      w(ixB^S,mag(2)) = By1_norm
      w(ixB^S,mag(3)) = Bz1_norm
    case(2)
      w(ixB^S,rho_)   = rho2_norm
      w(ixB^S,mom(1)) = rho2_norm*vx2_norm
      w(ixB^S,mom(2)) = rho2_norm*vy2_norm
      w(ixB^S,mom(3)) = rho2_norm*vz2_norm 
      w(ixB^S,e_)     = eg2_norm
      w(ixB^S,r_e)    = Er2_norm
      w(ixB^S,mag(1)) = Bx2_norm
      w(ixB^S,mag(2)) = By2_norm
      w(ixB^S,mag(3)) = Bz2_norm
    case default
      call mpistop('boundary not known')
    end select
  end subroutine boundary_conditions

  subroutine refine_shock(igrid,level,ixG^L,ix^L,qt,w,x,refine,coarsen)
    ! Enforce additional refinement or coarsening
    ! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.
    ! you must set consistent values for integers refine/coarsen:
    ! refine = -1 enforce to not refine
    ! refine =  0 doesn't enforce anything
    ! refine =  1 enforce refinement
    ! coarsen = -1 enforce to not coarsen
    ! coarsen =  0 doesn't enforce anything
    ! coarsen =  1 enforce coarsen
    use mod_global_parameters

    integer, intent(in) :: igrid, level, ixG^L, ix^L
    double precision, intent(in) :: qt, w(ixG^S,1:nw), x(ixG^S,1:ndim)
    integer, intent(inout) :: refine, coarsen

    if (any(dabs(x(ixG^S,1)) < 0.1d0)) refine=1

  end subroutine refine_shock

  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)

    double precision :: wlocal(ixI^S,nw)
    double precision :: lamb(ixI^S), R(ixI^S)

    wlocal(ixI^S,1:nw)=w(ixI^S,1:nw)
    call fld_get_fluxlimiter(wlocal,x,ixI^L,ixO^L,lamb,R,2,fld_fl)
    w(ixO^S,nw+1)=lamb(ixO^S)
    w(ixO^S,nw+2)=R(ixO^S)
  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
  character(len=*) :: varnames
  varnames='Lambda R'

  end subroutine specialvarnames_output

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
    w(ixO^S,vely_)=w(ixO^S,mom(2))/w(ixO^S,rho_)
    ! output the AMR level (assuming uniform grid blocks)
    w(ixO^S,amr_)=dlog(((xprobmax1-xprobmin1)/domain_nx1)/dxlevel(1))/dlog(2.0d0)+1.0d0

  end subroutine set_output_vars

end module mod_usr
