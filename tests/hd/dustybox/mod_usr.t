module mod_usr
  use mod_hd

  double precision            :: min_ar_ = 5.0d-8
  double precision            :: max_ar_ = 250.0d-8
  double precision, parameter :: rho1_   = 1.0d-20
  double precision, parameter :: vel1_   = 5.0d4
  double precision, parameter :: T1_     = 1.0d2

contains

  subroutine usr_init()

    unit_length=1.0d18        ! in cm
    unit_numberdensity=1.0d3  ! in g/cc
    unit_velocity=1.0d6       ! in cm/s

    call set_coordinate_system("Cartesian_1D")

    usr_init_one_grid => initonegrid_usr
    usr_set_parameters => initglobaldata_usr

    call hd_activate()

  end subroutine usr_init

  subroutine initglobaldata_usr
    use mod_dust
    integer          :: i
    double precision :: r(0:dust_n_species)
    logical, save    :: first = .true.

    hd_gamma                 = 5.0d0/3.0d0
    dust_density(:)          = 3.3d0   ! specific density of dust in g/cc

    dust_density(:) = dust_density(:)/unit_density
    min_ar_         = min_ar_/unit_length
    max_ar_         = max_ar_/unit_length

    ! first dust sizes in Ndust bins, with all bins having equal total mass.
    ! To do this, assume the particle distribution goes as r^-3.5
    r(0) = min_ar_
    do i=1,dust_n_species
      r(i) = (dsqrt(r(i-1)) +(dsqrt(max_ar_)- &
           dsqrt(min_ar_))/dust_n_species)**2.0d0
    end do

    ! now calculate the weigthed mean size of each bin, again assuming n goes as r^-3.5
    do i=1,dust_n_species
      dust_size(i) = (5.0d0/3.0d0)*(r(i)**(-1.5d0) - r(i-1)**(-1.5d0)) &
           /(r(i)**(-2.5d0) - r(i-1)**(-2.5d0))
    end do


      if(first)then
         if(mype==0)then
            write(*,*) '*****************************************'
            if(SI_unit)then
               write(*,*) 'Units system in SI'
            else
               write(*,*) 'Units system in cgs'
            endif
            write(*,*) 'He_abundance is       =',He_abundance
            write(*,*) 'unit length is        =',unit_length
            write(*,*) 'unit number density is=',unit_numberdensity
            write(*,*) 'unit velocity is      =',unit_velocity
            write(*,*) 'unit time is          =',unit_time
            write(*,*) 'unit density is       =',unit_density
            write(*,*) 'unit pressure is      =',unit_pressure
            write(*,*) 'unit temperature is   =',unit_temperature
            write(*,*) 'specific heat ratio is=',hd_gamma
            write(*,*) '*****************************************'
            write(*,*) 'Dust included using ',dust_n_species,' dust species'
            write(*,*) 'Dust bins all have specific density rhop ',dust_density(1)
            write(*,*) '   in cgs units specific density is rhop ',dust_density(1)*unit_density
            write(*,*) 'Dust bins between min=',min_ar_,' and max=',max_ar_
            do i=1,dust_n_species
              write(*,*) 'Dust type ',i,': grain radius r              =', dust_size(i)*unit_length
              write(*,*) 'Dust type ',i,': dimensionless grain radius r=', dust_size(i)
              write(*,*) 'Dust type ',i,': dimensionless rhop x r       =', dust_size(i)*dust_density(i)
            end do
            write(*,*) '*****************************************'
         endif
         first=.false.
      endif

    do i=1,nw
      if(loglimit(i) .and. (.not. (i==rho_ .or. i==p_))) then
        call mpistop('Bad idea to limit the logarithm of possibly negative numbers')
      endif
    enddo

  end subroutine initglobaldata_usr

  subroutine initonegrid_usr(ixG^L,ix^L,w,x)
    use mod_dust
    integer, intent(in) :: ixG^L, ix^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    integer :: n

    w(ix^S,mom(1:ndir))   = 0.0d0
    w(ix^S,rho_)          = rho1_/unit_density
    w(ix^S,p_)            = (rho1_* T1_ / (unit_density*unit_temperature))

    do n = 1, dust_n_species
      w(ix^S, dust_rho(n))          = 0.01d0*rho1_/(unit_density*dust_n_species)
      w(ix^S, dust_mom(1:ndir, n))  = vel1_/unit_velocity
    end do

    call hd_to_conserved(ixG^L,ix^L,w,x)
  end subroutine initonegrid_usr

end module mod_usr
