module mod_radiative_cooling

  #:if defined('COOLING')

  #:mute
  #:include "../hd/mod_hd_templates.fpp"
  #:include "../ffhd/mod_ffhd_templates.fpp"
  #:endmute

! For the interpolatable tables: these tables contain log_10 temperature values and corresponding
! log_10 luminosity values. The simulation-dependent temperature and luminosity
! scaling parameters are supposed to be provided in the user file. 
! All tables have been extended to at least T=10^9 K using a pure Bremsstrahlung 
! relationship of Lambda~sqrt(T). This to ensure that a purely explicit calculation 
! without timestep check is only used for extremely high temperatures. 
! (Except for the SPEX curve, which is more complicated and therefore simply stops  
! at the official upper limit of log(T) = 8.16)

  use mod_global_parameters
  use mod_physics_vars
  use mod_comm_lib, only: mpistop

  implicit none

  !> to be pushed to the device
  !> Helium abundance over Hydrogen
  double precision, private    :: He_abundance_rc

  !> The adiabatic index
  double precision, private :: rc_gamma

  !> The adiabatic index minus 1
  double precision, private :: rc_gamma_1

  !> inverse of the adiabatic index minus 1
  double precision, private :: invgam

  !$acc declare create(rc_gamma_1, invgam)

  type rc_fluid

    double precision :: rad_cut_hgt
    double precision :: rad_cut_dey

    ! these are set in init method
    double precision, allocatable :: tcool(:), Lcool(:)
    double precision, allocatable :: Yc(:)
    double precision  :: tref, lref, tcoolmin,tcoolmax
    double precision  :: lgtcoolmin, lgtcoolmax, lgstep

    !> Lower limit of temperature
    double precision   :: tlow

    !> Coefficent of cooling time step
    double precision   :: cfrac

    !> Index of cut off temperature for TRAC
    integer              :: Tcoff_

    ! these are set as parameters
    !> Resolution of temperature in interpolated tables
    integer :: ncool

    integer :: n_PPL

    !> Fixed temperature not lower than tlow
    logical   :: Tfix

    !> cutoff radiative cooling below rad_cut_hgt
    logical :: rad_cut

    !> Name of cooling curve
    character(len=std_len)  :: coolcurve

    !> Name of cooling method
    character(len=std_len)  :: coolmethod

  end type rc_fluid

  type(rc_fluid)   :: rc_fl
  !$acc declare create(rc_fl)

  ! Interpolatable tables

  double precision :: t_JCcorona(1:45), t_DM_2(1:76), t_Colgan(1:55) 

  double precision :: l_JCcorona(1:45), l_DM_2(1:76), l_Colgan(1:55) 

  integer          :: n_JCcorona, n_DM_2 , n_Colgan

  data    n_JCcorona / 45 /

  data    t_JCcorona / 4.00000, 4.14230, 4.21995, 4.29761, 4.37528, &
                       4.45294, 4.53061, 4.60827, 4.68593, 4.76359, &
                       4.79705, 4.83049, 4.86394, 4.89739, 4.93084, &
                       4.96428, 4.99773, 5.03117, 5.06461, 5.17574, &
                       5.28684, 5.39796, 5.50907, 5.62018, 5.73129, &
                       5.84240, 5.95351, 6.06461, 6.17574, 6.28684, &
                       6.39796, 6.50907, 6.62018, 6.73129, 6.84240, &
                       6.95351, 7.06461, 7.17574, 7.28684, 7.39796, &
                       7.50907, 7.62018, 7.73129, 7.84240, 7.95351  /

  data    l_JCcorona / -200.18883, -100.78630, -30.60384, -22.68481, -21.76445, &
                       -21.67936, -21.54218, -21.37958, -21.25172, -21.17584, &
                       -21.15783, -21.14491, -21.13527, -21.12837, -21.12485, &
                       -21.12439, -21.12642, -21.12802, -21.12548, -21.08965, &
                       -21.08812, -21.19542, -21.34582, -21.34839, -21.31701, &
                       -21.29072, -21.28900, -21.34104, -21.43122, -21.62448, &
                       -21.86694, -22.02897, -22.08051, -22.06057, -22.01973, &
                       -22.00000, -22.05161, -22.22175, -22.41452, -22.52581, &
                       -22.56914, -22.57486, -22.56151, -22.53969, -22.51490  /

  !
  ! To be used together with the SPEX table for the SPEX_DM option
  ! Assuming an ionization fraction of 10^-3
  !
  data    n_DM_2 / 76 /

  data    t_DM_2 / 1.00, 1.04, 1.08, 1.12, 1.16, 1.20 &
                 , 1.24, 1.28, 1.32, 1.36, 1.40 &
                 , 1.44, 1.48, 1.52, 1.56, 1.60 &
                 , 1.64, 1.68, 1.72, 1.76, 1.80 &
                 , 1.84, 1.88, 1.92, 1.96, 2.00 &
                 , 2.04, 2.08, 2.12, 2.16, 2.20 &
                 , 2.24, 2.28, 2.32, 2.36, 2.40 &
                 , 2.44, 2.48, 2.52, 2.56, 2.60 & 
                 , 2.64, 2.68, 2.72, 2.76, 2.80 &
                 , 2.84, 2.88, 2.92, 2.96, 3.00 &
                 , 3.04, 3.08, 3.12, 3.16, 3.20 &
                 , 3.24, 3.28, 3.32, 3.36, 3.40 &
                 , 3.44, 3.48, 3.52, 3.56, 3.60 & 
                 , 3.64, 3.68, 3.72, 3.76, 3.80 &
                 , 3.84, 3.88, 3.92, 3.96, 4.00 /

  data    l_DM_2 / -30.0377, -29.7062, -29.4055, -29.1331, -28.8864, -28.6631 &
                 , -28.4614, -28.2791, -28.1146, -27.9662, -27.8330 &
                 , -27.7129, -27.6052, -27.5088, -27.4225, -27.3454 &
                 , -27.2767, -27.2153, -27.1605, -27.1111, -27.0664 &
                 , -27.0251, -26.9863, -26.9488, -26.9119, -26.8742 &
                 , -26.8353, -26.7948, -26.7523, -26.7080, -26.6619 &
                 , -26.6146, -26.5666, -26.5183, -26.4702, -26.4229 &
                 , -26.3765, -26.3317, -26.2886, -26.2473, -26.2078 &
                 , -26.1704, -26.1348, -26.1012, -26.0692, -26.0389 &
                 , -26.0101, -25.9825, -25.9566, -25.9318, -25.9083 &
                 , -25.8857, -25.8645, -25.8447, -25.8259, -25.8085 &
                 , -25.7926, -25.7778, -25.7642, -25.7520, -25.7409 &
                 , -25.7310, -25.7222, -25.7142, -25.7071, -25.7005 &
                 , -25.6942, -25.6878, -25.6811, -25.6733, -25.6641 &
                 , -25.6525, -25.6325, -25.6080, -25.5367, -25.4806  /

  data    n_Colgan / 55 /

  data    t_Colgan / 4.06460772, 4.14229559, 4.21995109, 4.29760733, 4.37527944, 4.45293587, &
                     4.53060946, 4.60826923, 4.68592974, 4.76359269, 4.79704583, 4.83049243, &
                     4.86394114, 4.89738514, 4.93083701, 4.96428321, 4.99773141, 5.03116600, &
                     5.06460772, 5.17574368, 5.28683805, 5.39795738, 5.50906805, 5.62017771, &
                     5.73129054, 5.84240328, 5.95351325, 6.06460772, 6.17574368, 6.28683805, &
                     6.39795738, 6.50906805, 6.62017771, 6.73129054, 6.84240328, 6.95351325, &
                     7.06460772, 7.17574368, 7.28683805, 7.39795738, 7.50906805, 7.62017771, &
                     7.73129054, 7.84240328, 7.95351325, 8.06460772, 8.17574368, 8.28683805, &
                     8.39795738, 8.50906805, 8.62017771, 8.73129054, 8.84240328, 8.95351325, & 
                     9.06460772                                                              /

  data    l_Colgan / -22.18883401, -21.78629635, -21.60383554, -21.68480662, -21.76444630, &
                     -21.67935529, -21.54217864, -21.37958284, -21.25171892, -21.17584161, &
                     -21.15783402, -21.14491111, -21.13526945, -21.12837453, -21.12485189, &
                     -21.12438898, -21.12641785, -21.12802448, -21.12547760, -21.08964778, &
                     -21.08812360, -21.19542445, -21.34582346, -21.34839251, -21.31700703, &
                     -21.29072156, -21.28900309, -21.34104468, -21.43122351, -21.62448270, &
                     -21.86694036, -22.02897478, -22.08050874, -22.06057061, -22.01973295, &
                     -22.00000434, -22.05161149, -22.22175466, -22.41451671, -22.52581288, &
                     -22.56913516, -22.57485721, -22.56150512, -22.53968863, -22.51490350, &
                     -22.48895932, -22.46071057, -22.42908363, -22.39358639, -22.35456791, &
                     -22.31261375, -22.26827428, -22.22203698, -22.17422996, -22.12514145  /

  contains

  !> get physical quantities from simulation result in specific physics module
  @:get_rho()
  @:get_pthermal()
  @:get_Rfactor()

    !> Radiative cooling initialization
    subroutine radiative_cooling_init_params(phys_gamma,He_abund)
      use mod_global_parameters
      double precision, intent(in) :: phys_gamma,He_abund

      rc_gamma=phys_gamma
      He_abundance_rc=He_abund
    end subroutine radiative_cooling_init_params


    !> the params require to be read in phys template file phys_init
    subroutine radiative_cooling_init(fl)
      use mod_global_parameters

      type(rc_fluid), intent(inout) :: fl
  
      double precision, dimension(:), allocatable :: t_table
      double precision, dimension(:), allocatable :: L_table
      double precision :: ratt, fact1, fact2, fact3, dL1, dL2
      double precision :: tstep, Lstep
      integer :: ntable, i, j
      logical :: jump

      call rc_params_read(fl)

      ! Init for interpolatable tables
      allocate(fl%tcool(1:fl%ncool), fl%Lcool(1:fl%ncool))
      allocate(fl%Yc(1:fl%ncool))

      fl%tcool(1:fl%ncool)    = zero
      fl%Lcool(1:fl%ncool)    = zero

      ! Read in the selected cooling curve
      select case(fl%coolcurve)
      case('JCcorona')
        if(mype ==0) &
        print *,'Use Colgan & Feldman (2008) cooling curve'
        if(mype ==0) &
        print *,'This version only till 10000 K, beware for floor T treatment'
        ntable = n_JCcorona
        allocate(t_table(1:ntable))
        allocate(L_table(1:ntable))
        t_table(1:ntable) = t_JCcorona(1:n_JCcorona)
        L_table(1:ntable) = l_JCcorona(1:n_JCcorona)

      case('Colgan_DM')
        if(mype==0)&
        print *, 'Combination of Colgan (2008) for high temperatures and'
        if(mype==0)&
        print *, 'Dalgarno & McCray (1972), DM2, for low temperatures'
        ntable = n_Colgan + n_DM_2
        allocate(t_table(1:ntable))
        allocate(L_table(1:ntable))
        t_table(1:n_DM_2) = t_DM_2(1:n_DM_2)
        L_table(1:n_DM_2) = L_DM_2(1:n_DM_2)
        t_table(n_DM_2+1:ntable) = t_Colgan(1:n_Colgan)
        L_table(n_DM_2+1:ntable) = l_Colgan(1:n_Colgan)

      case default
        call mpistop("This coolingcurve is unknown")
      end select

      ! create cooling table(s) for use in amrvac
      fl%tcoolmax = t_table(ntable)
      fl%tcoolmin = t_table(1)
      ratt = (fl%tcoolmax-fl%tcoolmin)/(dble(fl%ncool-1) + smalldouble)

      fl%tcool(1) = fl%tcoolmin
      fl%Lcool(1) = L_table(1)

      fl%tcool(fl%ncool) = fl%tcoolmax
      fl%Lcool(fl%ncool) = L_table(ntable)

      do i=2,fl%ncool        ! loop to create one table
        fl%tcool(i) = fl%tcool(i-1)+ratt
        do j=1,ntable-1   ! loop to create one spot on a table
        ! Second order polynomial interpolation, except at the outer edge,
        ! or in case of a large jump.
          if(fl%tcool(i) < t_table(j+1)) then
            if(j.eq. ntable-1 )then
              fact1 = (fl%tcool(i)-t_table(j+1))     &
                    /(t_table(j)-t_table(j+1)) 
              fact2 = (fl%tcool(i)-t_table(j))       &
                    /(t_table(j+1)-t_table(j)) 
              fl%Lcool(i) = L_table(j)*fact1 + L_table(j+1)*fact2 
              exit
            else 
              dL1 = L_table(j+1)-L_table(j)
              dL2 = L_table(j+2)-L_table(j+1)
              jump =(max(dabs(dL1),dabs(dL2)) > 2*min(dabs(dL1),dabs(dL2)))
            end if
            if( jump ) then
              fact1 = (fl%tcool(i)-t_table(j+1))     &
                    /(t_table(j)-t_table(j+1)) 
              fact2 = (fl%tcool(i)-t_table(j))       &
                    /(t_table(j+1)-t_table(j)) 
              fl%Lcool(i) = L_table(j)*fact1 + L_table(j+1)*fact2
              exit          
            else
              fact1 = ((fl%tcool(i)-t_table(j+1))     &
                    * (fl%tcool(i)-t_table(j+2)))   &
                    / ((t_table(j)-t_table(j+1)) &
                    * (t_table(j)-t_table(j+2)))
              fact2 = ((fl%tcool(i)-t_table(j))       &
                    * (fl%tcool(i)-t_table(j+2)))   &
                    / ((t_table(j+1)-t_table(j)) &
                    * (t_table(j+1)-t_table(j+2)))
              fact3 = ((fl%tcool(i)-t_table(j))       &
                    * (fl%tcool(i)-t_table(j+1)))   &
                    / ((t_table(j+2)-t_table(j)) &
                    * (t_table(j+2)-t_table(j+1)))
              fl%Lcool(i) = L_table(j)*fact1 + L_table(j+1)*fact2 &
                        + L_table(j+2)*fact3
              exit
            end if
          end if
        end do  ! end loop to find create one spot on a table
      end do    ! end loop to create one table

      ! Go from logarithmic to actual values.
      fl%tcool(1:fl%ncool) = 10.0D0**fl%tcool(1:fl%ncool)
      fl%Lcool(1:fl%ncool) = 10.0D0**fl%Lcool(1:fl%ncool)

      ! Change unit of table if SI is used instead of cgs
      if (si_unit) fl%Lcool(1:fl%ncool) = fl%Lcool(1:fl%ncool) * 10.0d0**(-13)

      ! Scale both T and Lambda
      fl%tcool(1:fl%ncool) = fl%tcool(1:fl%ncool) / unit_temperature
      fl%Lcool(1:fl%ncool) = fl%Lcool(1:fl%ncool) * unit_numberdensity**2 * unit_time / unit_pressure * (1.d0+2.d0*He_abundance_rc) 

      fl%tcoolmin       = fl%tcool(1)+smalldouble  ! avoid pointless interpolation
      ! smaller value for lowest temperatures from cooling table and user's choice
      if (fl%tlow==bigdouble) fl%tlow=fl%tcoolmin
      fl%tcoolmax       = fl%tcool(fl%ncool)
      fl%lgtcoolmin = dlog10(fl%tcoolmin)
      fl%lgtcoolmax = dlog10(fl%tcoolmax)
      fl%lgstep = (fl%lgtcoolmax-fl%lgtcoolmin) * 1.d0 / (fl%ncool-1)

      deallocate(t_table)
      deallocate(L_table)

      if( fl%coolmethod == 'exact' ) then
        fl%tref = fl%tcoolmax
        fl%lref = fl%Lcool(fl%ncool)
        fl%Yc(fl%ncool) = zero
        do i=fl%ncool-1, 1, -1
            fl%Yc(i) = fl%Yc(i+1)
            do j=1,100
              tstep = 1.0d-2*(fl%tcool(i+1)-fl%tcool(i))
              call findL(fl%tcool(i+1)-j*tstep, Lstep, fl)
              fl%Yc(i) = fl%Yc(i) + fl%lref/fl%tref*tstep/Lstep
            end do
        end do
      end if

      rc_gamma_1=rc_gamma-1.d0
      invgam = 1.d0/rc_gamma_1
      !$acc update device(rc_gamma_1, invgam)

    contains

      subroutine rc_params_read(fl)

        use mod_global_parameters, only: unitpar,par_files
        use mod_constants, only: bigdouble
        type(rc_fluid), intent(inout) :: fl
        integer                      :: n
        integer :: ncool = 4000
        double precision :: cfrac=0.1d0
    
        !> Name of cooling curve
        character(len=std_len)  :: coolcurve='JCcorona'
    
        !> Name of cooling method
        character(len=std_len)  :: coolmethod='exact'
    
        !> Fixed temperature not lower than tlow
        logical    :: Tfix=.false.
    
        !> Lower limit of temperature
        double precision   :: tlow=bigdouble
    
        logical    :: rad_cut=.false.
        double precision :: rad_cut_hgt=0.5d0
        double precision :: rad_cut_dey=0.15d0

        namelist /rc_list/ coolcurve, coolmethod, ncool, cfrac, tlow, Tfix, rad_cut, rad_cut_hgt, rad_cut_dey

        do n = 1, size(par_files)
          open(unitpar, file=trim(par_files(n)), status="old")
          read(unitpar, rc_list, end=111)
  111       close(unitpar)
        end do

        fl%ncool=ncool
        fl%coolcurve=coolcurve
        fl%coolmethod=coolmethod
        fl%tlow=tlow
        fl%Tfix=Tfix
        fl%cfrac=cfrac
        fl%rad_cut=rad_cut
        fl%rad_cut_hgt=rad_cut_hgt
        fl%rad_cut_dey=rad_cut_dey

      end subroutine rc_params_read

    end subroutine radiative_cooling_init

    subroutine radiative_cooling_add_source(qdt,wCT,wCTprim,wnew,x,fl)
      !$acc routine seq

    ! w[iw]=w[iw]+qdt*S[wCT,x] where S is the source based on wCT within ixO
      double precision, intent(in) :: qdt, wCT(nw_phys), wCTprim(nw_phys)
      double precision, intent(inout) :: wnew(nw_phys)
      double precision, intent(in) :: x(1:ndim)
      type(rc_fluid), intent(in) :: fl

      select case(fl%coolmethod)
      case ('exact')   
        call cool_exact(qdt,wCT,wCTprim,wnew,x,fl)
      case default
        call mpistop("This cooling method is unknown")
      end select
      if( fl%Tfix ) call floortemperature(qdt,wCT,wnew,x,fl)

    end subroutine radiative_cooling_add_source

    subroutine floortemperature(qdt,wCT,wnew,x,fl)
      !$acc routine seq
    !  Force minimum temperature to a fixed temperature
      double precision, intent(in)    :: qdt, wCT(nw_phys)
      double precision, intent(inout) :: wnew(nw_phys)
      double precision, intent(in)    :: x(1:ndim)
      type(rc_fluid), intent(in) :: fl
      double precision :: etherm, rho, emin, Rfactor

      etherm = get_pthermal(wCT,x)  
      rho = get_rho(wCT,x)  
      Rfactor = get_Rfactor()
      emin = rho*fl%tlow*Rfactor
      if(etherm < emin) then
        wnew(iw_e)=wnew(iw_e)+(emin-etherm)*invgam
      end if
    end subroutine floortemperature

    subroutine cool_exact(qdt,wCT,wCTprim,wnew,x,fl)
      !$acc routine seq
    !  Cooling routine using exact integration method from Townsend 2009
      double precision, intent(in)    :: qdt, wCT(nw_phys), wCTprim(nw_phys), x(1:ndim)
      double precision, intent(inout) :: wnew(nw_phys)
      type(rc_fluid), intent(in) :: fl
      double precision :: Y1, Y2
      double precision :: L1, Tlocal2, pnew
      double precision :: rho, Te, rhonew, Rfactor
      double precision :: emin, Lmax, fact
      double precision :: de, emax

      rho = get_rho(wCT,x)
      Rfactor = get_Rfactor()
      Te=wCTprim(iw_e)/(rho*Rfactor)
      pnew = get_pthermal(wCT,x)
      rhonew = get_rho(wnew,x)

      fact = fl%lref*qdt/fl%tref

      emin = rhonew*fl%tlow*Rfactor*invgam
      Lmax = max(zero,pnew*invgam-emin)/qdt
      emax = max(zero,pnew*invgam-emin)
      !  Determine explicit cooling
      !  If temperature is below floor level, no cooling.
      !  Stop wasting time and go to next gridpoint.
      !  If the temperature is higher than the maximum,
      !  assume Bremsstrahlung
      if( Te<=fl%tcoolmin ) then
        L1 = zero
      else if( Te>=fl%tcoolmax )then
        call calc_l_extended(Te, L1, fl)
        L1 = L1*rho**2

        L1 = min(L1,Lmax)
        if(slab_uniform .and. fl%rad_cut .and. x(ndim) .le. fl%rad_cut_hgt) then
          L1 = L1*exp(-(x(ndim)-fl%rad_cut_hgt)**2/fl%rad_cut_dey**2)
        end if
        wnew(iw_e) = wnew(iw_e)-L1*qdt
      else
        call findY(Te,Y1,fl)
        Y2 = Y1 + fact*rho*rc_gamma_1
        call findT(Tlocal2,Y2,fl)
        if(Tlocal2<=fl%tcoolmin) then
          de = emax
        else
          de = (Te-Tlocal2)*rho*Rfactor*invgam
        end if
        de = min(de,emax)
        if(slab_uniform .and. fl%rad_cut .and. x(ndim) .le. fl%rad_cut_hgt) then
          de = de*exp(-(x(ndim)-fl%rad_cut_hgt)**2/fl%rad_cut_dey**2)
        end if
        wnew(iw_e) = wnew(iw_e)-de
      end if
    end subroutine cool_exact

    subroutine calc_l_extended (tpoint, lpoint,fl)
      !$acc routine seq
    !  Calculate l for t beyond tcoolmax
    !  Assumes Bremsstrahlung for the interpolated tables
    !  Uses the power law for piecewise power laws
      double precision, intent(IN)  :: tpoint
      double precision, intent(OUT) :: lpoint
      type(rc_fluid), intent(in) :: fl

      lpoint = fl%Lcool(fl%ncool) * dsqrt( tpoint / fl%tcoolmax)

    end subroutine calc_l_extended

    subroutine findL (tpoint,Lpoint,fl)
    !  Fast search option to find correct point 
    !  in cooling curve
      !$acc routine seq

      double precision,intent(IN)   :: tpoint
      double precision, intent(OUT) :: Lpoint
      type(rc_fluid), intent(in) :: fl

      double precision :: lgtp
      integer :: jl,jc,jh,i

      lgtp = dlog10(tpoint)
      jl = int((lgtp - fl%lgtcoolmin) /fl%lgstep) + 1
      Lpoint = fl%Lcool(jl)+ (tpoint-fl%tcool(jl)) &
                * (fl%Lcool(jl+1)-fl%Lcool(jl)) &
                / (fl%tcool(jl+1)-fl%tcool(jl))

    end subroutine findL

    subroutine findY (tpoint,Ypoint,fl)
    !  Fast search option to find correct point in cooling time (TEF)
      !$acc routine seq

      double precision,intent(IN)   :: tpoint
      double precision, intent(OUT) :: Ypoint
      type(rc_fluid), intent(in) :: fl

      double precision :: lgtp
      double precision :: y_extra,factor
      integer :: jl,jc,jh,i

      lgtp = dlog10(tpoint)
      jl = int((lgtp - fl%lgtcoolmin) / fl%lgstep) + 1
      Ypoint = fl%Yc(jl)+ (tpoint-fl%tcool(jl)) &
                * (fl%Yc(jl+1)-fl%Yc(jl)) &
                / (fl%tcool(jl+1)-fl%tcool(jl))

    end subroutine findY

    subroutine findT (tpoint,Ypoint,fl)
    !  Fast search option to find correct temperature 
    !  from temporal evolution function. Only possible this way because T is a monotonously
    !  decreasing function for the interpolated tables
    !  Uses eq. A7 from Townsend 2009 for piecewise power laws
      !$acc routine seq

      double precision,intent(OUT)   :: tpoint
      double precision, intent(IN) :: Ypoint
      type(rc_fluid), intent(in) :: fl

      double precision :: factor
      integer :: jl,jc,jh,i

      if(Ypoint >= fl%Yc(1)) then
        tpoint = fl%tcoolmin
      else if (Ypoint == fl%Yc(fl%ncool)) then
        tpoint = fl%tcoolmax
      else
        jl=0
        jh=fl%ncool+1
        do
          if(jh-jl <= 1) exit
          jc=(jh+jl)/2
          if(Ypoint <= fl%Yc(jc)) then
            jl=jc
          else
            jh=jc
          end if
        end do
        ! Linear interpolation to obtain correct temperature
        tpoint = fl%tcool(jl)+ (Ypoint-fl%Yc(jl)) &
                * (fl%tcool(jl+1)-fl%tcool(jl)) &
                / (fl%Yc(jl+1)-fl%Yc(jl))
      end if

    end subroutine findT

    subroutine getvar_cooling_exact(qdt, wCT, w, x, coolrate, fl)
      !$acc routine seq
      ! Calculates cooling rate using the exact cooling method,
      ! for usage in eg. source_terms subroutine.
      ! The TEF must be known, so this routine can only be used
      ! together with the "exact" cooling method.
  
      use mod_global_parameters
      double precision, intent(in)    :: qdt, wCT(nw_phys), w(nw_phys), x(1:ndim)
      double precision, intent(inout) :: coolrate
      type(rc_fluid), intent(in) :: fl
      double precision :: Y1, Y2
      double precision :: L1, Tlocal2, pnew
      double precision :: rho, pth, Te, rhonew, Rfactor
      double precision :: emin, Lmax, fact
      double precision :: de, emax

      rho = get_rho(wCT,x)
      pth = get_pthermal(wCT,x)
      Rfactor = get_Rfactor()
      Te=pth/(rho*Rfactor)
      pnew = get_pthermal(wCT,x)
      rhonew = get_rho(w,x)

      fact = fl%lref*qdt/fl%tref

      emin = rhonew*fl%tlow*Rfactor*invgam
      Lmax = max(zero,pnew*invgam-emin)/qdt
      emax = max(zero,pnew*invgam-emin)
      !  Determine explicit cooling
      !  If temperature is below floor level, no cooling.
      !  Stop wasting time and go to next gridpoint.
      !  If the temperature is higher than the maximum,
      !  assume Bremsstrahlung
      if( Te<=fl%tcoolmin ) then
        de = zero
      else if( Te>=fl%tcoolmax )then
        call calc_l_extended(Te, L1, fl)
        L1 = L1*rho**2

        L1 = min(L1,Lmax)
        if(slab_uniform .and. fl%rad_cut .and. x(ndim) .le. fl%rad_cut_hgt) then
          L1 = L1*exp(-(x(ndim)-fl%rad_cut_hgt)**2/fl%rad_cut_dey**2)
        end if
        de = L1*qdt
      else
        call findY(Te,Y1,fl)
        Y2 = Y1 + fact*rho*rc_gamma_1
        call findT(Tlocal2,Y2,fl)
        if(Tlocal2<=fl%tcoolmin) then
          de = emax
        else
          de = (Te-Tlocal2)*rho*Rfactor*invgam
        end if
        de = min(de,emax)
        if(slab_uniform .and. fl%rad_cut .and. x(ndim) .le. fl%rad_cut_hgt) then
          de = de*exp(-(x(ndim)-fl%rad_cut_hgt)**2/fl%rad_cut_dey**2)
        end if
      end if
      coolrate = de/qdt
      end subroutine getvar_cooling_exact
    
  #:endif

end module mod_radiative_cooling
