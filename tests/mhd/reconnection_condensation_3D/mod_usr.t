module mod_usr
  use mod_mhd
  implicit none

  ! initial parameters of system 
  double precision, parameter :: La = 20.d0, T0 = 1.d6, n0 = 1.d9, sol_grav = 274.0d2
  double precision            :: unit_grav, unit_currentdensity, unit_eta, Ls, as, alpha
  double precision            :: beta, Tini, ncor, grav_norm, miu0_cgs, t1, t2, t1_user, t2_user
  double precision            :: v00, v00_usr, solar_radius, B_a, perturbation_str
  double precision            :: T_threshold_refine, T_threshold_refine_norm, special_dt_min
  integer                     :: icase, number_field_lines, number_field_lines_forward_backward
  logical                     :: perturbation
  double precision            :: t1_low_refine, t2_high_refine

  ! Additional variables
  integer                     :: j1_, j2_, j3_, normalDivB_, eta_, Temp_, mtsn2_, pgrad2_, H_, Lambda_

contains

  subroutine usr_init()    
    call usr_params_read(par_files)

    ! Choose independent normalization units if using dimensionless variables.
    unit_length               = 1.d8 ! cm
    unit_temperature          = 1.d6 ! K
    unit_numberdensity        = 1.d9 ! cm-3

    ! Define auxiliary variables that will be used in other subroutines
    usr_set_parameters        => initparam_usr

    ! A routine for initial conditions is always required
    usr_init_one_grid         => initonegrid_usr

    ! We use magnetic field splitting, so define background magnetic field and current density
    usr_set_B0                => specialset_B0
    usr_set_J0                => specialset_J0
    usr_init_vector_potential => specialset_A0

    ! Special boundary conditions 
    usr_special_bc            => specialbound_usr

    ! Gravity
    usr_gravity               => calc_gravity
    
    ! Resistivity to drive magnetic reconnection
    usr_special_resistivity   => special_eta

    ! Background heating in energy equation
    usr_source                => special_source

    ! Refine criteria
    usr_refine_grid           => special_refine_grid

    ! Calculate variables along magnetic field lines
    usr_set_field_w           => special_field_w
    
    ! Add additional variables to the vtu output.
    usr_aux_output           => specialvar_output
    usr_add_aux_names        => specialvarnames_output       

    ! Custom processing
    usr_special_convert       => custom_processing

    ! Add additional variables to the dat output 
    usr_modify_output         => set_output_vars

    call set_coordinate_system("Cartesian_3D")
    call mhd_activate()

    ! Add additional variables to the dat file
    j1_                       = var_set_extravar("j1","j1")
    j2_                       = var_set_extravar("j2","j2") 
    j3_                       = var_set_extravar("j3","j3")
    !normalDivB_               = var_set_extravar("normalDivB","normalDivB")  
    !eta_                      = var_set_extravar("eta","eta")
    Temp_                     = var_set_extravar("Te","Te")
    mtsn2_                    = var_set_extravar("mtsn2","mtsn2")
    pgrad2_                   = var_set_extravar("pgrad2","pgrad2")
    H_                        = var_set_extravar("H", "H")
    Lambda_                   = var_set_extravar("Lambda","Lambda")
  end subroutine usr_init

!======================================================================================
!> Read this module's parameters from a file.par
!======================================================================================
subroutine usr_params_read(files)
  character(len=*), intent(in) :: files(:)
  integer                      :: n

  namelist /usr_list/ v00_usr !Should match your naming and elements of your usr.par namelist
  namelist /usr_list/ perturbation
  namelist /usr_list/ perturbation_str
  namelist /usr_list/ T_threshold_refine
  namelist /usr_list/ special_dt_min
  namelist /usr_list/ t1_user
  namelist /usr_list/ t2_user
  
  do n = 1, size(files)
     open(unitpar, file=trim(files(n)), status="old")
     read(unitpar, usr_list, end=111)
     111 close(unitpar)
  end do

end subroutine usr_params_read

  subroutine initparam_usr()
    ! Units
    !miu0_cgs             = 1.d0!miu0_si / 1.d1
    !unit_density        = 1.4d0*mh_cgs*unit_numberdensity               
    !unit_magneticfield  = dsqrt(miu0_cgs*unit_pressure)
    !unit_velocity       = unit_magneticfield/dsqrt(miu0_cgs*unit_density)  
    !unit_time           = unit_length/unit_velocity     
    !unit_pressure       = 2.3d0*unit_numberdensity*kB_cgs*unit_temperature                 
    unit_mass           = unit_density / unit_numberdensity
    unit_currentdensity = unit_magneticfield/unit_length/4.d0/dpi
    unit_grav           = unit_velocity**2 / unit_length          
    unit_eta            = unit_length * unit_velocity

    ! Solar constants
    solar_radius        = 6.957d10 / unit_length


    ! Auxiliary variables

    ! We multiply with 2.3d0 because this factor is within unit_pressure
    Tini                = T0 / unit_temperature             
    
    ncor                = n0 / unit_numberdensity
    grav_norm           = sol_grav / unit_grav

    Ls                  = Tini / grav_norm
    as                  = Ls

    alpha               = -2.d0 * La / (dpi * as)
    beta                = -dsqrt(1.d0 - alpha**2)

    ! Footpoint motion variables
    t1                  = t1_user / unit_time
    t2                  = t2_user/ unit_time
    v00                 = v00_usr / unit_velocity
    B_a                 = Busr / unit_magneticfield
    
    ! Refinement variables
    T_threshold_refine_norm = T_threshold_refine / unit_temperature
    t1_low_refine           = 42.d0 * 60.d0 / unit_time
    t2_high_refine          = 65.d0 * 60.d0 / unit_time

    ! Read number of field lines from first line in file
    open(unit=25, file="seedPoints.txt", status="old", action="read")
    read(25, *) number_field_lines
    close(unit=25)

    ! One entire streamline consists of one forward and one backward calculated streamline
    number_field_lines_forward_backward = 2 * number_field_lines

    if (mype == 0) then
      print*, 'unit_numberdensity  = ', unit_numberdensity
      print*, 'unit_density        = ', unit_density
      print*, 'unit_pressure       = ', unit_pressure
      print*, 'unit_velocity       = ', unit_velocity
      print*, 'unit_magneticfield  = ', unit_magneticfield
      print*, 'unit_time           = ', unit_time
      print*, 'unit_grav           = ', unit_grav
      print*, 'unit_mass           = ', unit_mass
      print*, 'unit_currentdensity = ', unit_currentdensity
      print*, 'unit_eta            = ', unit_eta
      print*, 'Ls                  = ', Ls
    end if
  end subroutine initparam_usr

  subroutine initonegrid_usr(ixI^L, ixO^L, w, x)
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision                :: argx(ixI^S), argy(ixI^S)

    if(mhd_glm) w(ixO^S,psi_)=0.d0

    argx(ixI^S)         = dpi * x(ixI^S, 1) / (2.0d0 * La)
    argy(ixI^S)         = -x(ixI^S, 2) / as

    ! Density under stratisfied gravity
    !w(ixI^S, rho_)      = ncor * exp( -(solar_radius / Ls) + (1.d0 / Ls) * (solar_radius**2 / (solar_radius + x(ixI^S, 2))))

    w(ixI^S, rho_)      = ncor * dexp( -x(ixI^S, 2) / Ls)

    ! static equilibrium
    w(ixI^S, mom(1))    = 0.d0
    w(ixI^S, mom(2))    = 0.d0
    w(ixI^S, mom(3))    = 0.d0


    ! pressure from ideal gas law
    w(ixI^S, p_)        = w(ixI^S, rho_) * Tini 

    ! magnetic field splitting so zero here
    if (B0field) then
      w(ixI^S, mag(1))  = 0.0d0
      w(ixI^S, mag(2))  = 0.0d0
      w(ixI^S, mag(3))  = 0.0d0

    else if(stagger_grid) then
      call b_from_vector_potential(block%ixGs^L,ixI^L,ixO^L,block%ws,x)
      call mhd_face_to_center(ixO^L,block)

      w(ixI^S,mag(3)) =  beta * B_a * cos(argx(ixI^S)) * exp(argy(ixI^S))

    else
      w(ixI^S, mag(1))  =  alpha * B_a * cos(argx(ixI^S)) * exp(argy(ixI^S))
      w(ixI^S, mag(2))  =          B_a * sin(argx(ixI^S)) * exp(argy(ixI^S))
      w(ixI^S, mag(3))  =   beta * B_a * cos(argx(ixI^S)) * exp(argy(ixI^S))
    end if

    ! convert primitive variables to conservative variables
    call mhd_to_conserved(ixI^L,ixO^L,w,x)
  end subroutine initonegrid_usr

  subroutine specialset_B0(ixI^L,ixO^L,x,wB0)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: wB0(ixI^S,1:ndir)

    double precision                :: argx(ixI^S), argy(ixI^S)
    
    argx(ixI^S)   = dpi * x(ixI^S, 1) / (2.0d0 * La)
    argy(ixI^S)   = -x(ixI^S, 2) / as

    wB0(ixI^S, 1) = alpha * B_a * cos(argx(ixI^S)) * exp(argy(ixI^S))
    wB0(ixI^S, 2) =         B_a * sin(argx(ixI^S)) * exp(argy(ixI^S))
    wB0(ixI^S, 3) =  beta * B_a * cos(argx(ixI^S)) * exp(argy(ixI^S))

  end subroutine specialset_B0

  subroutine specialset_J0(ixI^L,ixO^L,x,wJ0)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: wJ0(ixI^S,1:ndir)

    double precision                :: argx(ixI^S), argy(ixI^S)

    argx(ixI^S)   = dpi * x(ixI^S, 1) / (2.0d0 * La)
    argy(ixI^S)   = -x(ixI^S, 2) / as

    ! miu0_cgs is equal to one since we are now working in a unitless environment
    wJ0(ixI^S, 1) =                  - beta / as * B_a * cos(argx(ixI^S)) * exp(argy(ixI^S))
    wJ0(ixI^S, 2) =        - B_a * beta / (as * alpha) * sin(argx(ixI^S)) * exp(argy(ixI^S))
    wJ0(ixI^S, 3) =   (alpha - alpha**(-1)) * B_a / as * cos(argx(ixI^S)) * exp(argy(ixI^S))
  end subroutine specialset_J0

  ! Necessary for constraint transport and keeping magnetic field divergence-less
  subroutine specialset_A0(ixI^L, ixC^L, xC, A, idir)
    use mod_global_parameters
    integer, intent(in)                :: ixI^L, ixC^L, idir
    double precision, intent(in)       :: xC(ixI^S,1:ndim)
    double precision, intent(out)      :: A(ixI^S)

    double precision                   :: argx(ixI^S), argy(ixI^S)
  
    argx(ixC^S) = dpi * xC(ixC^S, 1) / (2.0d0 * La)
    argy(ixC^S) = -xC(ixC^S, 2) / as

    if (idir == 2) then
      A(ixC^S)  = -beta * as * B_a * alpha * sin(argx(ixC^S)) * exp(argy(ixC^S))
    else if (idir == 3) then
      A(ixC^S)  =        -as * B_a * alpha * cos(argx(ixC^S)) * exp(argy(ixC^S))
    else
      A(ixC^S)  = 0.d0
    end if
  end subroutine specialset_A0

  subroutine specialbound_usr(qt,ixI^L,ixO^L,iB,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L
    integer, intent(in)             :: ixO^L
    integer, intent(in)             :: iB
    double precision, intent(in)    :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision                :: v0, argx(ixI^S), argy(ixI^S), argz(ixI^S)
    double precision                :: Qp(ixI^S), pth(ixI^S)
    integer                         :: ix, idir, ixOs^L, ixA^L

    if(mhd_glm) w(ixO^S,psi_)=0.d0

    call v0func(global_time, v0)

    argx(ixO^S)                     = dpi * x(ixO^S, 1) / La
    argy(ixO^S)                     = -xprobmax2 / as
    argz(ixO^S)                     = -(x(ixO^S, 3) / La)**2

    whichboundary: select case(iB) 

    case(1)
      w(ixO^S, mom(:))              = 0.d0

      do ix = ixOmax1,ixOmin1,-1
        w(ix^%1ixO^S,rho_)      = w(ixOmax1+1^%1ixO^S,rho_)
        w(ix^%1ixO^S,p_)        = w(ix^%1ixO^S,rho_) * Tini

        ! Second order zero gradient extrapolation
        w(ix^%1ixO^S, mag(1))   =   1.d0 / 3.d0 * (4.d0 * w(ix+1^%1ixO^S, mag(1)) &
                                                        - w(ix+2^%1ixO^S, mag(1)))
        w(ix^%1ixO^S, mag(2))   =   1.d0 / 3.d0 * (4.d0 * w(ix+1^%1ixO^S, mag(2)) &
                                                        - w(ix+2^%1ixO^S, mag(2)))
        w(ix^%1ixO^S, mag(3))   =   1.d0 / 3.d0 * (4.d0 * w(ix+1^%1ixO^S, mag(3)) &
                                                        - w(ix+2^%1ixO^S, mag(3)))
      end do
    
    case(2)
      w(ixO^S, mom(:))              = 0.d0

      do ix = ixOmin1,ixOmax1
        w(ix^%1ixO^S,rho_)      = w(ixOmin1-1^%1ixO^S,rho_)
        w(ix^%1ixO^S,p_)        = w(ix^%1ixO^S,rho_) * Tini

        ! Second order zero gradient extrapolation
        w(ix^%1ixO^S, mag(1))   =   1.d0 / 3.d0 * (4.d0 * w(ix-1^%1ixO^S, mag(1)) &
                                                        - w(ix-2^%1ixO^S, mag(1)))
        w(ix^%1ixO^S, mag(2))   =   1.d0 / 3.d0 * (4.d0 * w(ix-1^%1ixO^S, mag(2)) &
                                                        - w(ix-2^%1ixO^S, mag(2)))
        w(ix^%1ixO^S, mag(3))   =   1.d0 / 3.d0 * (4.d0 * w(ix-1^%1ixO^S, mag(3)) &
                                                        - w(ix-2^%1ixO^S, mag(3)))
      end do

    case(3)                         ! Bottom boundary
      ! Density under Stratified Gravity 
      !w(ixO^S, rho_)                = ncor * exp( -(solar_radius / Ls) + (1.d0 / Ls) * (solar_radius**2 / (solar_radius + x(ixO^S, 2))))
      w(ixO^S, rho_)                = ncor * exp( -x(ixO^S, 2) / Ls) 
      w(ixO^S, p_)                  = w(ixO^S, rho_) * Tini

      w(ixO^S, mom(1))              = -v0 * sin(argx(ixO^S)) * dexp(argz(ixO^S))
      w(ixO^S, mom(2))              = 0.d0
      w(ixO^S, mom(3))              = 0.d0!w(ixO^S, mom(1))

      if (stagger_grid) then
        do idir=1,nws
          if(idir==2) cycle
          ixOsmax^D=ixOmax^D;
          ixOsmin^D=ixOmin^D-kr(^D,idir);
          do ix=ixOsmax2,ixOsmin2,-1
            block%ws(ix^%2ixOs^S, idir)       = 1.d0 / 3.d0 * (4.d0 * block%ws(ix+1^%2ixOs^S, idir)&
                                                                    - block%ws(ix+2^%2ixOs^S, idir))

          end do
        end do

        ixOs^L=ixO^L-kr(2,^D);
        block%ws(ixOs^S,2)=zero

        do ix=ixOsmax2,ixOsmin2,-1
          call get_divb(w,ixI^L,ixO^L,Qp)
          block%ws(ix^%2ixOs^S,2)=Qp(ix+1^%2ixO^S)*block%dvolume(ix+1^%2ixO^S)&
            /block%surfaceC(ix^%2ixOs^S,2)
        end do
        call mhd_face_to_center(ixO^L,block)

      else
        do ix=ixOmax2,ixOmin2,-1
          w(ix^%2ixO^S, mag(1))        =  1.d0 / 3.d0 * (4.d0 * w(ix+1^%2ixO^S, mag(1)) &
                                                              - w(ix+2^%2ixO^S, mag(1)))
          w(ix^%2ixO^S, mag(2))        =  1.d0 / 3.d0 * (4.d0 * w(ix+1^%2ixO^S, mag(2)) &
                                                              - w(ix+2^%2ixO^S, mag(2)))
          w(ix^%2ixO^S, mag(3))        =  1.d0 / 3.d0 * (4.d0 * w(ix+1^%2ixO^S, mag(3)) &
                                                              - w(ix+2^%2ixO^S, mag(3)))

          if(mhd_glm) w(ix^%2ixO^S,psi_)=w(ixOmax2+1^%2ixO^S,psi_)
        end do
      end if

    case(4)                         ! Top boundary
      w(ixO^S, mom(1))              = 0.0d0
      w(ixO^S, mom(2))              = 0.0d0
      w(ixO^S, mom(3))              = 0.0d0

      do ix = ixOmin2, ixOmax2
        w(ix^%2ixO^S,rho_)      = w(ixOmin2-1^%2ixO^S,rho_)
        w(ix^%2ixO^S,p_)        = w(ix^%2ixO^S,rho_) * Tini

        ! Second order zero gradient extrapolation
        w(ix^%2ixO^S, mom(1))   = 1.d0 / 3.d0 * (4.d0  &
                                  * w(ix-1^%2ixO^S,mom(1)) / w(ix-1^%2ixO^S,rho_)  &
                                  - w(ix-2^%2ixO^S,mom(1)) / w(ix-2^%2ixO^S,rho_))      
        w(ix^%2ixO^S, mom(2))   = 1.d0 / 3.d0 * (4.d0  &
                                  * w(ix-1^%2ixO^S,mom(2)) / w(ix-1^%2ixO^S,rho_)  &
                                  - w(ix-2^%2ixO^S,mom(2)) / w(ix-2^%2ixO^S,rho_))   
        w(ix^%2ixO^S, mom(3))   = 1.d0 / 3.d0 * (4.d0  &
                                  * w(ix-1^%2ixO^S,mom(3)) / w(ix-1^%2ixO^S,rho_)  &
                                  - w(ix-2^%2ixO^S,mom(3)) / w(ix-2^%2ixO^S,rho_))
      end do

      ! No inflow
      w(ixO^S, mom(2)) = max(w(ixO^S, mom(2)),zero)

      ! Second order zero gradient extrapolation 
      if (stagger_grid) then
        do idir=1,nws
          if(idir==2) cycle
          ixOsmax^D=ixOmax^D;
          ixOsmin^D=ixOmin^D-kr(^D,idir);
          do ix=ixOsmin2,ixOsmax2
            !block%ws(ix^%2ixOs^S, idir)       =   3.d0 * block%ws(ix-1^%2ixOs^S, idir) &
            !                                     -3.d0 * block%ws(ix-2^%2ixOs^S, idir) &
            !                                           + block%ws(ix-3^%2ixOs^S, idir) 

            block%ws(ix^%2ixOs^S, idir)       = 1.d0 / 3.d0 * (4.d0 * block%ws(ix-1^%2ixOs^S, idir)&
                                                                    - block%ws(ix-2^%2ixOs^S, idir))

          end do
        end do
        ixOs^L=ixO^L;
        block%ws(ixOs^S,2)=zero
        do ix=ixOsmin2,ixOsmax2
          call get_divb(w,ixI^L,ixO^L,Qp)
          block%ws(ix^%2ixOs^S,2)=-Qp(ix^%2ixO^S)*block%dvolume(ix^%2ixO^S)&
            /block%surfaceC(ix^%2ixOs^S,2)
        end do
        call mhd_face_to_center(ixO^L,block)


      else
        do ix=ixOmin2,ixOmax2
          w(ix^%2ixO^S, mag(1))        =   1.d0 / 3.d0 * (4.d0 * w(ix-1^%2ixO^S, mag(1)) &
                                                               - w(ix-2^%2ixO^S, mag(1)))
          w(ix^%2ixO^S, mag(2))        =   1.d0 / 3.d0 * (4.d0 * w(ix-1^%2ixO^S, mag(2)) &
                                                               - w(ix-2^%2ixO^S, mag(2)))
          w(ix^%2ixO^S, mag(3))        =   1.d0 / 3.d0 * (4.d0 * w(ix-1^%2ixO^S, mag(3)) &
                                                               - w(ix-2^%2ixO^S, mag(3)))
        end do
      end if
    
    case(5)
      w(ixO^S, mom(:))              = 0.d0

      do ix = ixOmax3,ixOmin3,-1
        w(ix^%3ixO^S,rho_)      = w(ixOmax3+1^%3ixO^S,rho_)
        w(ix^%3ixO^S,p_)        = w(ix^%3ixO^S,rho_) * Tini

        ! Second order zero gradient extrapolation
        w(ix^%3ixO^S, mag(1))   =   1.d0 / 3.d0 * (4.d0 * w(ix+1^%3ixO^S, mag(1)) &
                                                        - w(ix+2^%3ixO^S, mag(1)))
        w(ix^%3ixO^S, mag(2))   =   1.d0 / 3.d0 * (4.d0 * w(ix+1^%3ixO^S, mag(2)) &
                                                        - w(ix+2^%3ixO^S, mag(2)))
        w(ix^%3ixO^S, mag(3))   =   1.d0 / 3.d0 * (4.d0 * w(ix+1^%3ixO^S, mag(3)) &
                                                        - w(ix+2^%3ixO^S, mag(3)))
      end do

    case(6)
      w(ixO^S, mom(:))              = 0.d0

      do ix = ixOmin3,ixOmax3
        w(ix^%3ixO^S,rho_)      = w(ixOmin3-1^%3ixO^S,rho_)
        w(ix^%3ixO^S,p_)        = w(ix^%3ixO^S,rho_) * Tini

        ! Second order zero gradient extrapolation
        w(ix^%3ixO^S, mag(1))   =   1.d0 / 3.d0 * (4.d0 * w(ix-1^%3ixO^S, mag(1)) &
                                                        - w(ix-2^%3ixO^S, mag(1)))
        w(ix^%3ixO^S, mag(2))   =   1.d0 / 3.d0 * (4.d0 * w(ix-1^%3ixO^S, mag(2)) &
                                                        - w(ix-2^%3ixO^S, mag(2)))
        w(ix^%3ixO^S, mag(3))   =   1.d0 / 3.d0 * (4.d0 * w(ix-1^%3ixO^S, mag(3)) &
                                                        - w(ix-2^%3ixO^S, mag(3)))
      end do
    end select whichboundary
  
    call mhd_to_conserved(ixI^L,ixO^L,w,x)

  end subroutine specialbound_usr


  subroutine v0func(global_time, v0)
    double precision, intent(in)  :: global_time
    double precision, intent(out) :: v0

    v0 = 0.d0
    if (global_time < t1) then
      v0  = v00
    
    else if ((t1 <= global_time) .and. (global_time < t2)) then
      v0  = v00 * (t2 - global_time) / (t2 - t1)

    else
      v0  = 0.0d0
    end if
  end subroutine v0func

  subroutine calc_gravity(ixI^L,ixO^L,wCT,x,gravity_field) 
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(in)    :: wCT(ixI^S,1:nw)
    double precision, intent(out)   :: gravity_field(ixI^S,ndim)

    ! gravity_field is a vector with only a vertical y-component
    gravity_field           = 0.d0
    gravity_field(ixI^S, 2) = -grav_norm !* (solar_radius / (solar_radius + x(ixI^S, 2)))**2
  end subroutine calc_gravity

  subroutine special_eta(w,ixI^L,ixO^L,idirmin,x,current,eta)
    use mod_global_parameters
    integer, intent(in)           :: ixI^L, ixO^L, idirmin
    double precision, intent(in)  :: w(ixI^S,nw), x(ixI^S,1:ndim)
    double precision              :: current(ixI^S,7-2*ndir:3), eta(ixI^S)

    double precision              :: jcrit, eta0, etamax, index_point(ixI^S), lightspeed

    lightspeed            = 29979245800.d0 

    ! Divide by lightspeed because KY2017 define index_point = c / 4pi * curl(B)
    jcrit                 = 25.d0 / (lightspeed*unit_currentdensity)
    eta0                  = 3.6d13 / unit_eta
    etamax                = 1.8d14 / unit_eta
    index_point(ixI^S)    = dsqrt(current(ixI^S,1)**2 + current(ixI^S,2)**2 + current(ixI^S,3)**2)

    eta(ixI^S)            = 0.d0
    if (global_time <= t2) then
      where (  index_point(ixO^S) >  jcrit) eta(ixO^S) = eta0 * (index_point(ixO^S)/jcrit - 1.d0)**2
      where (eta(ixO^S) > etamax) eta(ixO^S) = etamax

      where ( (x(ixO^S,1) <  -9.d0) .or. (x(ixO^S,1) >  9.d0)) eta(ixO^S) = 0.d0 * x(ixO^S,1)
      where ( (x(ixO^S,3) < -40.d0) .or. (x(ixO^S,3) > 40.d0)) eta(ixO^S) = 0.d0 * x(ixO^S,3)
    end if
  end subroutine special_eta

  subroutine special_source(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)
    !use mod_radiative_cooling
    integer, intent(in) :: ixI^L, ixO^L, iw^LIM
    double precision, intent(in) :: qdt, qtC, qt
    double precision, intent(in) :: x(ixI^S,1:ndim), wCT(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: bQgrid(ixI^S)

    bQgrid(ixI^S) = 0.d0
    call getbQ(bQgrid,ixI^L,ixO^L,qtC,wCT,x)

    w(ixO^S,e_)=w(ixO^S,e_)+qdt*bQgrid(ixO^S)
  end subroutine special_source

  subroutine getbQ(bQgrid,ixI^L,ixO^L,qt,w,x)
    ! calculate background heating bQ
      use mod_radiative_cooling
      integer, intent(in)             :: ixI^L, ixO^L
      double precision, intent(in)    :: qt, x(ixI^S,1:ndim), w(ixI^S,1:nw)
      double precision, intent (out)  :: bQgrid(ixI^S)
  
      double precision                :: alpha_A, lambda_val1, lambda_val2
  
      call findL(Tini,lambda_val1,rc_fl)

      alpha_A   =  ncor**2 * lambda_val1
      bQgrid(ixO^S) = alpha_A * dexp(-2.d0 * x(ixO^S, 2) / as)
  end subroutine getbQ

    !> Enforce additional refinement or coarsening
   !> One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.
   !> you must set consistent values for integers refine/coarsen:
   !> refine = -1 enforce to not refine
   !> refine =  0 doesn't enforce anything
   !> refine =  1 enforce refinement
   !> coarsen = -1 enforce to not coarsen
   !> coarsen =  0 doesn't enforce anything
   !> coarsen =  1 enforce coarsen
   !> e.g. refine for negative first coordinate x < 0 as
   !> if (any(x(ix^S,1) < zero)) refine=1
  subroutine special_refine_grid(igrid,level,ixI^L,ixO^L,qt,w,x,refine,coarsen)
    use mod_global_parameters
    use mod_mhd_phys
    integer, intent(in)          :: igrid, level, ixI^L, ixO^L
    double precision, intent(in) :: qt, w(ixI^S,1:nw), x(ixI^S,1:ndim)
    integer, intent(inout)       :: refine, coarsen

    double precision :: pth(ixI^S), temp(ixI^S)

    ! Coarsen the entire grid except for prominence regions
    refine    = -1
    coarsen   =  1

    yif: if (any(x(ixO^S,2) <= 25.d0)) then
      ylvlif: if (level > 3) then
        refine  = -1
        coarsen =  1

      else if (level == 3) then
        refine  = 0
        coarsen = 0

      else 
        refine  = 1
        coarsen = -1
      end if ylvlif
    end if yif

    call mhd_get_pthermal(w,x,ixI^L,ixO^L,pth)
    temp(ixI^S)=pth(ixI^S)/w(ixI^S,rho_)
    tempif: if (any(temp(ixO^S)<=T_threshold_refine_norm)) then
        refine    =  1
        coarsen   = -1
    end if tempif 

    xif: if (all(abs(x(ixO^S,1)) >= 7.5d0)) then
      refine  = -1
      coarsen =  1
    end if xif

    yif2: if (all(x(ixO^S,2) >= 25.d0)) then
      refine  = -1
      coarsen =  1
    end if yif2

    zif: if (any(abs(x(ixO^S,3)) >= 40.d0)) then
      refine  = -1
      coarsen =  1
    end if zif

  end subroutine special_refine_grid

  subroutine trace_Bfield_parallel(step_length,max_points_in_field_line,number_variables_along_field_line,number_line_properties)
    use mod_point_searching
    use mod_trace_field

    integer :: max_points_in_field_line,number_variables_along_field_line,number_line_properties
    double precision :: step_length

    double precision :: xs^L, x_seed, y_seed, z_seed
    double precision :: xpp(ndim),wpp(nw)
    double precision :: line_array(number_field_lines_forward_backward,max_points_in_field_line,ndim),variables_along_line_array(number_field_lines_forward_backward,max_points_in_field_line,number_variables_along_field_line),line_properties_array(number_field_lines_forward_backward,number_line_properties+1)
    integer :: number_valid_points_array(number_field_lines_forward_backward)
    integer :: index_line,index_point
    logical :: forwardm(number_field_lines_forward_backward)
    character(len=std_len) :: ftype,tcondi,fname
    
    ! Read seedpoints from file 
    open(unit=13, file="seedPoints.txt",status="old", action="read")

    ! First line in file tells the amount of field lines
    read(13, *) 
    do index_line = 1,number_field_lines
      read(13, 12) x_seed, y_seed, z_seed
      12 format(3F10.5)

      line_array(2*index_line-1,1,:)  = [x_seed, y_seed, z_seed]
      line_array(2*index_line,1,:)    = [x_seed, y_seed, z_seed]

      forwardm(2*index_line-1)        = .false.
      forwardm(2*index_line)          = .true.
    end do
    close(unit=13)

    ftype='Bfield'
    tcondi='user1'

    call trace_field_multi(line_array,variables_along_line_array,line_properties_array,step_length,number_field_lines_forward_backward,max_points_in_field_line,number_variables_along_field_line,number_line_properties,forwardm,ftype,tcondi)
    
    number_valid_points_array(:)=int(line_properties_array(:,1))
    call writeFieldLines(line_array,variables_along_line_array,line_properties_array, number_field_lines_forward_backward, max_points_in_field_line, number_variables_along_field_line, number_line_properties, number_valid_points_array, step_length)
  
  end subroutine trace_Bfield_parallel
  
  subroutine special_field_w(igrid,ip,xf,wP,wL,max_points_in_field_line,number_variables_along_field_line,number_line_properties,step_length,forward,ftype,tcondi)
    use mod_global_parameters
    use mod_point_searching

    integer, intent(in)                 :: igrid,ip,max_points_in_field_line,number_variables_along_field_line,number_line_properties
    double precision, intent(in)        :: xf(max_points_in_field_line,ndim)
    double precision, intent(inout)     :: wP(max_points_in_field_line,number_variables_along_field_line),wL(1+number_line_properties)
    double precision, intent(in)        :: step_length
    logical, intent(in)                 :: forward
    character(len=std_len), intent(in)  :: ftype,tcondi

    double precision :: coordinate_point(1:ndim), variable_value_at_point(1:nw)
    integer          :: ipe, ixc^D

    coordinate_point(1:ndim)=xf(ip,1:ndim)
    variable_value_at_point = zero

    variable_value_at_point = zero
    call get_point_w_ingrid(igrid,coordinate_point,variable_value_at_point,'primitive')
    wP(ip,1)= variable_value_at_point(rho_)
    wP(ip,2)= variable_value_at_point(mom(1))
    wP(ip,3)= variable_value_at_point(mom(2))
    wP(ip,4)= variable_value_at_point(mom(3))
    wP(ip,5)= variable_value_at_point(p_)
    wP(ip,6)= variable_value_at_point(mag(1))
    wP(ip,7)= variable_value_at_point(mag(2))
    wP(ip,8)= variable_value_at_point(mag(3))
    wP(ip,9)= variable_value_at_point(p_) / variable_value_at_point(rho_)                


    ! Current density
    wP(ip,10)=variable_value_at_point(j1_)
    wP(ip,11)=variable_value_at_point(j2_)
    wP(ip,12)=variable_value_at_point(j3_)

    ! Heating and cooling
    wP(ip,13)=variable_value_at_point(H_)
    wP(ip,14)=variable_value_at_point(Lambda_)

    ! ip is the index of the point. 
    if (ip==1) then
      wL(2)=0.d0
    
    ! Here, we calculate the length of the field line for each point
    else
      wL(2)=wL(2)+step_length
    endif

  end subroutine special_field_w

  subroutine writeFieldLines(line_array,variables_along_line_array,line_properties_array, number_field_lines_forward_backward, max_points_in_field_line, number_variables_along_field_line, number_line_properties, number_valid_points_array, step_length)
    integer, intent(in)          :: number_field_lines_forward_backward, max_points_in_field_line
    integer, intent(in)          :: number_variables_along_field_line, number_line_properties
    integer, intent(in)          :: number_valid_points_array(number_field_lines_forward_backward)
    double precision, intent(in) :: step_length
    double precision, intent(in) :: line_array(number_field_lines_forward_backward,max_points_in_field_line,ndim),variables_along_line_array(number_field_lines_forward_backward,max_points_in_field_line,number_variables_along_field_line),line_properties_array(number_field_lines_forward_backward,number_line_properties+1)

    character(len=std_len) :: fname1, fname2
    integer                :: index_line, index_point


    if (mype==0) then

      ! output field line data
      write(fname1, '(a,a,i4.4,a)') trim(base_filename), '_Blines_', snapshotnext,'.txt'
      open(UNIT=1, FILE=fname1, STATUS="REPLACE", ACTION="WRITE") 

      ! Header
      write(1, 99) 'line','point', 'x', 'y', 'z', 'rho', 'v1', 'v2', 'v3', 'p',&
                    'b1', 'b2', 'b3', 'Te', 'j1', 'j2', 'j3', 'H', 'Lambda'
      99 FORMAT(2A10, A12, 16A18)

      do index_line=1,number_field_lines_forward_backward
        do index_point=1,number_valid_points_array(index_line)
          write(1,88) index_line,&
                      index_point,&
                      line_array(index_line,index_point,1),&
                      line_array(index_line,index_point,2),&
                      line_array(index_line,index_point,3),&
                      variables_along_line_array(index_line,index_point,1),&
                      variables_along_line_array(index_line,index_point,2),&
                      variables_along_line_array(index_line,index_point,3),&
                      variables_along_line_array(index_line,index_point,4),&
                      variables_along_line_array(index_line,index_point,5),&
                      variables_along_line_array(index_line,index_point,6),&
                      variables_along_line_array(index_line,index_point,7),&
                      variables_along_line_array(index_line,index_point,8),&
                      variables_along_line_array(index_line,index_point,9),&
                      variables_along_line_array(index_line,index_point,10),&
                      variables_along_line_array(index_line,index_point,11),&
                      variables_along_line_array(index_line,index_point,12),&
                      variables_along_line_array(index_line,index_point,13),&
                      variables_along_line_array(index_line,index_point,14)
          88 FORMAT (2I10, 17E18.8)
        enddo
      enddo
      close(1)

      write(fname2, '(a,a,i4.4,a)') trim(base_filename), '_Blines_summary_', snapshotnext,'.txt'
      open(UNIT=2, FILE=fname2, STATUS="REPLACE", ACTION="WRITE") 

      ! General line data
      write(2, 100) 'line', 'n_point', 'length'
      100 FORMAT(2A10, A14)

      do index_line=1,number_field_lines_forward_backward
        write(2,89) index_line,&
                    number_valid_points_array(index_line),&
                    line_properties_array(index_line, 2)
        89 FORMAT (2I10, 12E18.8)
      enddo
      close(2)
    end if 
  end subroutine writeFieldLines

  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
    use mod_global_parameters
    use mod_mhd_phys
    use mod_radiative_cooling
    integer, intent(in)          :: ixI^L,ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision             :: w(ixI^S,nw+nwauxio)
    double precision             :: normconv(0:nw+nwauxio)
    double precision             :: current_o(ixI^S,3), normalDivB(ixI^S), divB(ixI^S), eta(ixI^S)
    double precision             :: pth(ixI^S), bQgrid(ixI^S), Qe(ixI^S), mtsn(ixI^S,1:ndir)
    double precision             :: pgrad2(ixI^S)
    integer                      :: idirmin
    logical                      :: fourthorder

    call get_current(w,ixI^L,ixO^L,idirmin,current_o)
    w(ixO^S,nw+1)=current_o(ixO^S,1) 
    w(ixO^S,nw+2)=current_o(ixO^S,2)
    w(ixO^S,nw+3)=current_o(ixO^S,3)

    !call get_normalized_divb(w,ixI^L,ixO^L,normalDivB)
    !w(ixO^S,nw+4)=normalDivB(ixO^S)

    !call special_eta(w,ixI^L,ixO^L,idirmin,x,current_o,eta)
    !w(ixO^S,nw+5)=eta(ixO^S) 

    call getmagnetictension(ixI^L,ixO^L,w,x,mtsn)
    w(ixO^S,nw+4)=mtsn(ixO^S,2)

    call gradient(w(ixI^S,p_),ixI^L,ixO^L,2,pgrad2) 
    w(ixO^S,nw+5)=pgrad2(ixO^S)

    ! output temperature
    call mhd_get_pthermal(w,x,ixI^L,ixO^L,pth)
    w(ixO^S,nw+6)=pth(ixO^S)/w(ixO^S,rho_)

    call getbQ(bQgrid,ixI^L,ixO^L,global_time,w,x)
    w(ixO^S,nw+7)=bQgrid(ixO^S)

    if(mhd_radiative_cooling) call getvar_cooling(ixI^L,ixO^L,w,x,Qe, rc_fl)
    w(ixO^S,nw+8)=Qe(ixO^S)

  end subroutine specialvar_output
  
  subroutine specialvarnames_output(varnames)
    use mod_global_parameters
    character(len=*) :: varnames
    varnames = 'j1 j2 j3 mtsn2 pgrad2 Te H Lambda'
  end subroutine specialvarnames_output

  ! Use aiconvert and set 'convert_type' in .par file to 'user'
  subroutine custom_processing(qunitconvert)
    use mod_global_parameters
    integer, intent(in) :: qunitconvert
    character(len=20)   :: userconvert_type

    integer             :: max_points_in_field_line,number_variables_along_field_line
    integer             :: number_line_properties
    double precision    :: max_length_field_line,step_length

    !---------------- Magnetic Streamlines -----------------!    
    ! Step length of field lines
    step_length=(xprobmax1-xprobmin1)/(domain_nx1*2**(refine_max_level-1))
    max_length_field_line=5.d0*sqrt((xprobmax1-xprobmin1)**2 + (xprobmax2-xprobmin2)**2 + &
                                     (xprobmax3-xprobmin3)**2)
    max_points_in_field_line=floor(max_length_field_line/step_length)+1

    ! Number of variables to compute along the field line
    number_variables_along_field_line=14

    ! Number of properties along the entire field line, e.g. maximum density
    number_line_properties=1

    call trace_Bfield_parallel(step_length,max_points_in_field_line,number_variables_along_field_line,number_line_properties)

  end subroutine custom_processing

  subroutine set_output_vars(ixI^L,ixO^L,qt,w,x)
    use mod_global_parameters
    use mod_radiative_cooling
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(in)    :: qt,x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision                :: normconv(0:nw+nwauxio)
    double precision                :: current_o(ixI^S,3), normalDivB(ixI^S), divB(ixI^S), eta(ixI^S)
    double precision                :: pth(ixI^S), bQgrid(ixI^S), Qe(ixI^S), mtsn(ixI^S,1:ndir)
    double precision                :: pgrad2(ixI^S)
    integer                         :: idirmin
    logical                         :: fourthorder

    !---------------- Current -----------------!
    call get_current(w,ixI^L,ixO^L,idirmin,current_o)
    w(ixO^S,j1_)=current_o(ixO^S,1) 
    w(ixO^S,j2_)=current_o(ixO^S,2)
    w(ixO^S,j3_)=current_o(ixO^S,3)


    !---------------- Divergence -----------------!
    !call get_normalized_divb(w,ixI^L,ixO^L,normalDivB)
    !w(ixO^S,normalDivB_)=normalDivB(ixO^S)


    !---------------- Resisitvity -----------------!
    !call special_eta(w,ixI^L,ixO^L,idirmin,x,current_o,eta)
    !w(ixO^S,eta_)=eta(ixO^S) 


    !---------------- Magnetic Tension -----------------!
    call getmagnetictension(ixI^L,ixO^L,w,x,mtsn)
    w(ixO^S,mtsn2_)=mtsn(ixO^S,2)


    !---------------- Pressure Gradient -----------------!
    call gradient(w(ixI^S,p_),ixI^L,ixO^L,2,pgrad2) 
    w(ixO^S,pgrad2_)=pgrad2(ixO^S)


    !---------------- Temperature -----------------!
    call mhd_get_pthermal(w,x,ixI^L,ixO^L,pth)
    w(ixO^S,Temp_)=pth(ixO^S)/w(ixO^S,rho_)


    !---------------- Background Heating -----------------!
    call getbQ(bQgrid,ixI^L,ixO^L,global_time,w,x)
    w(ixO^S,H_)=bQgrid(ixO^S)


    !---------------- Radiative cooling -----------------!
    if(mhd_radiative_cooling) call getvar_cooling(ixI^L,ixO^L,w,x,Qe, rc_fl)
    w(ixO^S,Lambda_)=Qe(ixO^S)    

  end subroutine set_output_vars

  subroutine getmagnetictension(ixI^L,ixO^L,w,x,mtsn)

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S,1:nw)
    double precision, intent(out) :: mtsn(ixI^S,1:ndir)

    double precision :: Bmag(ixI^S),bvec(ixI^S,1:ndir),tmp(ixI^S),x(ixI^S,1:ndim)
    integer :: idims, idir

    if(B0field) then
      !Bmag(ixI^S)=dsqrt(sum((w(ixI^S,mag(:))+block%B0(ixI^S,:,0))**2,dim=ndim))
      do idir=1,ndir
        bvec(ixI^S,idir)=(w(ixI^S,mag(idir))+block%B0(ixI^S,idir,0))
      end do
    else
      !Bmag(ixI^S)=dsqrt(sum(w(ixI^S,mag(:))**2,dim=ndim))
      do idir=1,ndir
        bvec(ixI^S,idir)=w(ixI^S,mag(idir))
      end do
    endif

    mtsn=0.d0
    ! calculate local curvature of magnetic field
    do idims=1,ndim
      call gradient(bvec(ixI^S,1),ixI^L,ixO^L,idims,tmp) 
      mtsn(ixO^S,1)=mtsn(ixO^S,1)+bvec(ixO^S,idims)*tmp(ixO^S) 
      call gradient(bvec(ixI^S,2),ixI^L,ixO^L,idims,tmp) 
      mtsn(ixO^S,2)=mtsn(ixO^S,2)+bvec(ixO^S,idims)*tmp(ixO^S)
      call gradient(bvec(ixI^S,3),ixI^L,ixO^L,idims,tmp) 
      mtsn(ixO^S,3)=mtsn(ixO^S,3)+bvec(ixO^S,idims)*tmp(ixO^S)
    end do 
  end subroutine getmagnetictension

end module mod_usr
