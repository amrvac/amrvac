module mod_usr
  use mod_hd

  double precision            :: min_ar_ = 5.0d-8
  double precision            :: max_ar_ = 250.0d-8
  double precision, parameter :: rho1_   = 1.0d-20
  double precision, parameter :: vel1_   = 5.0d4
  double precision, parameter :: T1_     = 1.0d2

contains

  subroutine usr_init()

    call set_coordinate_system("Cartesian_1D")

    usr_init_one_grid => initonegrid_usr
    usr_set_parameters => initglobaldata_usr

    call hd_activate()

  end subroutine usr_init

  subroutine initglobaldata_usr
    use mod_dust
    integer          :: i
    double precision :: r(0:dust_n_species)
    double precision :: dsdust(dust_n_species)
    logical, save    :: first = .true.

    hd_gamma                 = 5.0d0/3.0d0
    dust_density(:)          = 3.3d0   ! density in the dust
    length_convert_factor    = 1.0d18  ! normalization for distance
    w_convert_factor(rho_)   = 1.0d-21 ! normalization for rho
    w_convert_factor(mom(1)) = 1.0d7   ! normalization for speed

    time_convert_factor  = length_convert_factor/w_convert_factor(mom(1))
    w_convert_factor(p_) = w_convert_factor(rho_)*(w_convert_factor(mom(1))**2)

    w_convert_factor(dust_rho(:))    = w_convert_factor(rho_)
    w_convert_factor(dust_mom(1, :)) = w_convert_factor(mom(1))

    dust_density(:) = dust_density(:)/w_convert_factor(dust_rho(1))
    min_ar_         = min_ar_/length_convert_factor
    max_ar_         = max_ar_/length_convert_factor

    ! if not using "dustmethod='linear'", define dust_density(1:dust_n_species),
    ! sdust(1:dust_n_species) and gas_mu
    !-------------------------------

    ! here the dust sizes are defined. Note that this is
    ! now done differently for the method used in (van Marle et al. (2011)).

    ! first dust sizes in Ndust bins, with all bins having equal total mass.
    ! To do this, assume the particle distribution goes as r^-3.5

    r(0) = min_ar_
    do i=1,dust_n_species
      r(i) = (sqrt(r(i-1)) +(sqrt(max_ar_)- &
           sqrt(min_ar_))/dust_n_species)**2.0d0
      dsdust(i) = r(i)-r(i-1)
    end do

    ! now calculate the weigthed mean size of each bin, again assuming n goes as r^-3.5
    do i=1,dust_n_species
      dust_size(i) = (5.0d0/3.0d0)*(r(i)**(-1.5d0) - r(i-1)**(-1.5d0)) &
           /(r(i)**(-2.5d0) - r(i-1)**(-2.5d0))
    end do

    if(first)then
      if(mype==0)then
        do i=1,dust_n_species
          write(*,*) 'Dust type ',i,': grain radius r=',dust_size(i)
        end do
      endif
      first=.false.
    endif

    do i=1,nw
      if(loglimit(i) .and. (.not. (i==rho_ .or. i==p_))) then
        call mpistop('Bad idea to take the logarithm of a negative number')
      endif
    enddo

  end subroutine initglobaldata_usr

  subroutine initonegrid_usr(ixG^L,ix^L,w,x)
    use mod_dust
    integer, intent(in) :: ixG^L, ix^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    integer :: n
    double precision :: delr
    double precision :: vsound
    !----------------------------------------------------------------------------

    {^IFTWOD call mpistop("This is a 1D HDDust Riemann problem!!!")}
    {^IFTHREED call mpistop("This is a 1D HDDust Riemann problem!!!")}

    w(ix^S,mom(1))         = 0.0d0
    w(ix^S,rho_)           = rho1_/w_convert_factor(rho_)
    w(ixG^S,p_)            = rho1_*(kboltzmann_cgs/(hydrogen_mass_cgs*gas_mu)) &
         * T1_ / w_convert_factor(p_)

    do n = 1, dust_n_species
      w(ixG^S, dust_rho(n)) = 0.01d0*rho1_/(w_convert_factor(rho_)*dust_n_species)
      w(ixG^S, dust_mom(1, n))  = vel1_/w_convert_factor(mom(1))
    end do

    call hd_to_conserved(ixG^L,ix^L,w,x)
  end subroutine initonegrid_usr

end module mod_usr
