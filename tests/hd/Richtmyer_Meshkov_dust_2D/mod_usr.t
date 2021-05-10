!=============================================================================
! Richtmyer-Meshkov hydro with dust (2D)
! RM through planar shock impinging on inclined density discontinuity
! setting for 2D hydro with 4 dust species
module mod_usr
  use mod_hd

  double precision :: min_ar
  double precision :: max_ar

contains

  subroutine usr_init()

    unit_length = 0.001d0*3.08567758d18 ! 0.001 pc distance (in cm)
    unit_numberdensity = 1.0d3          ! 10^3 particles per cc
    unit_velocity = 1.0d6               ! normalization for speed cm/s

    call set_coordinate_system("Cartesian_2D")

    usr_init_one_grid => initonegrid_usr
    usr_set_parameters => initglobaldata_usr

    call hd_activate()
  end subroutine usr_init

  subroutine initglobaldata_usr
    use mod_dust
    double precision :: r(0:dust_n_species)
    integer          :: i
    logical, save    :: first = .true.

    hd_gamma=1.4d0

    if (dust_n_species > 0) then
      dust_density(1:dust_n_species)   = 3.3d0    ! dust grain density (g/cc)
      min_ar              = 1.0d-7   ! minimum dust grain size (cm)
      max_ar              = 500.0d-7 ! maximum dust grain size (cm)

      ! === rescale dust quantities to dimensionless scale === !
      dust_density(1:dust_n_species) = dust_density(1:dust_n_species)/unit_density
      min_ar=min_ar/unit_length
      max_ar=max_ar/unit_length

      ! here the dust sizes are defined. Ndust bins, with all bins having equal total mass.
      ! To do this, assume the particle distribution goes as r^-3.5
      r(0) = min_ar
      do i=1,dust_n_species
        r(i) = (dsqrt(r(i-1)) +(dsqrt(max_ar)- &
                    dsqrt(min_ar))/dust_n_species)**2.0d0
      end do

      ! now calculate the weighted mean size of each bin, again assuming n goes as r^-3.5
      ! here we just take average over r
      do i=1,dust_n_species
         dust_size(i) = (5.0d0/3.0d0)*(r(i)**(-1.5d0) - r(i-1)**(-1.5d0)) &
              /(r(i)**(-2.5d0) - r(i-1)**(-2.5d0))
      end do
      ! here we weigh distribution over r^2
      !do i=1,dust_n_species
      !   dust_size(i) = dsqrt(5.0d0*(r(i)**(-0.5d0) - r(i-1)**(-0.5d0)) &
      !                             /(r(i)**(-2.5d0) - r(i-1)**(-2.5d0)))
      !end do

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
            write(*,*) 'Dust inluded using ',dust_n_species,' dust species'
            write(*,*) 'Dust bins all have specific density rhop ',dust_density(1)
            write(*,*) '   in cgs units specific density is rhop ',dust_density(1)*unit_density
            write(*,*) 'Dust bins between min=',min_ar,' and max=',max_ar
            write(*,*) ' in cgs min=',min_ar*unit_length,' and max=',max_ar*unit_length
            do i=0,dust_n_species
               print *,r(i)
            enddo
            do i=1,dust_n_species
              write(*,*) 'Dust type ',i,': grain radius r              =', dust_size(i)*unit_length
              write(*,*) 'Dust type ',i,': dimensionless grain radius r=', dust_size(i)
              write(*,*) 'Dust type ',i,': dimensionless rhop x r      =', dust_size(i)*dust_density(i)
            end do
            write(*,*) '*****************************************'
         endif
         first=.false.
      endif
    endif

  end subroutine initglobaldata_usr

  ! initialize one grid
  subroutine initonegrid_usr(ixG^L,ix^L,w,x)
    use mod_dust
    integer, intent(in)             :: ixG^L, ix^L
    double precision, intent(in)    :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    integer                         :: n
    double precision                :: xshock, xbound, rhopost, vpost, ppost, alpha
    double precision                :: Prat,alfa,c_pre
    double precision                :: M,eta
    logical, save                   :: first = .true.
    !----------------------------------------------------------------------------

    {^IFONED   call mpistop("This is a 2D HD problem: Richtmyer Meshkov")}
    {^IFTHREED call mpistop("This is a 2D HD problem: Richtmyer Meshkov")}

    ! === Shock and CD position === !
    xshock=0.4d0
    xbound=0.5d0 ! xbound must be bigger than xshock!
    ! inclination angle for CD
    alpha = (3.14159265d0/4.0d0)

    ! Mach number of planar shock
    M     = 2.0d0
    ! inclination angle for CD
    alpha = (3.14159265d0/4.0d0)
    ! density ratio across CD
    eta   = 3.0d0

    ! compute the RH related states across the planar shock
    Prat    = one/(one+(M**2-one)*two*hd_gamma/(hd_gamma+one))
    alfa    = (hd_gamma+one)/(hd_gamma-one)
    c_pre   = one ! pre shock sound speed
    rhopost = hd_gamma*(alfa+Prat)/(alfa*Prat+one)
    ppost   = one/Prat
    vpost   = c_pre*M*(one-(alfa*Prat+one)/(alfa+Prat))

    if (mype==0.and.first) then
      print *,'============================================================='
      print *,' HD Richtmyer Meshkov simulation '
      print *,'============================================================='
      print *,' Mach number of shock: ',M
      print *,' post-shock density: ',rhopost
      print *,' post-shock velocity:',vpost
      print *,' post-shock pressure:',ppost
      print *,' Density ratio at CD: ',eta
      print *,'============================================================='
      first=.false.
    endif

    do n = 1, dust_n_species
      w(ix^S, dust_mom(:, n)) = 0.0d0

      where(x(ix^S,1)>xshock.and.(x(ix^S,1)>x(ix^S,2)/dtan(alpha)+xbound))
        w(ix^S,dust_rho(n)) = (0.01d0*hd_gamma*eta)/dust_n_species
      elsewhere(x(ix^S,1)>xshock.and.(x(ix^S,1)<=x(ix^S,2)/dtan(alpha)+xbound))
        w(ix^S,dust_rho(n)) = 0.0d0
      elsewhere
        w(ix^S,dust_rho(n)) = 0.0d0
      end where
    end do

    where(x(ix^S,1)>xshock.and.(x(ix^S,1)>x(ix^S,2)/dtan(alpha)+xbound))
      ! pre shock region
      w(ix^S,rho_)=hd_gamma*eta
      w(ix^S,mom(1))=zero
      w(ix^S,mom(2))=zero
      w(ix^S,e_)=one/(hd_gamma-one)
    elsewhere(x(ix^S,1)>xshock.and.(x(ix^S,1)<=x(ix^S,2)/dtan(alpha)+xbound))
      ! pre shock region
      w(ix^S,rho_)=hd_gamma
      w(ix^S,mom(1))=zero
      w(ix^S,mom(2))=zero
      w(ix^S,e_)=one/(hd_gamma-one)
    elsewhere
      ! post shock region
      w(ix^S,rho_)= rhopost
      w(ix^S,mom(1)) = rhopost*vpost
      w(ix^S,mom(2)) = zero
      w(ix^S,e_)  = ppost/(hd_gamma-one)+0.5d0*rhopost*vpost**2
    endwhere

  end subroutine initonegrid_usr

end module mod_usr
