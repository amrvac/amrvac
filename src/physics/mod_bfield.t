!> Minimal magnetic-field-only state for static field construction.
!>
!> The B-only path deliberately registers only the magnetic components as
!> conservative/primitive variables.  It is intended for static potential
!> (and other field-construction) cases; it does not provide a time-evolving
!> induction equation.  The legacy MHD-embedded magnetofriction path is
!> intentionally unaffected by this module.
module mod_bfield
  implicit none
  public

  !> Optional ionisation fractions used only by the unit normalisation.
  double precision :: bfield_h_ion_fr=1.d0
  double precision :: bfield_he_ion_fr=1.d0
  double precision :: bfield_he_ion_fr2=1.d0

contains

  !> Activate the B-only physics/state path.
  subroutine bfield_activate()
    call bfield_phys_init()
  end subroutine bfield_activate

  subroutine bfield_read_params(files)
    use mod_global_parameters

    character(len=*), intent(in) :: files(:)
    integer :: n

    namelist /bfield_list/ bfield_h_ion_fr,bfield_he_ion_fr,&
       bfield_he_ion_fr2,SI_unit

    do n=1,size(files)
      open(unitpar,file=trim(files(n)),status='old')
      read(unitpar,bfield_list,end=111)
111   close(unitpar)
    end do
  end subroutine bfield_read_params

  subroutine bfield_phys_init()
    use mod_global_parameters
    use mod_physics
    use mod_functions_bfield, only: mag
    use mod_eos_container, only: eos
    use mod_comm_lib, only: mpistop

    call bfield_read_params(par_files)

    physics_type='bfield'
    phys_energy=.false.
    stagger_grid=.false.

    allocate(mag(ndir))
    mag(:)=var_set_bfield(ndir)

    iwstart=mag(1)
    allocate(start_indices(number_species),stop_indices(number_species))
    start_indices(1)=mag(1)
    stop_indices(1)=nwflux
    nwgc=nwflux+nwaux
    nws=0

    nvector=1
    allocate(iw_vector(nvector))
    iw_vector(1)=mag(1)-1

    if(.not.allocated(flux_type)) then
      allocate(flux_type(ndir,nwflux))
      flux_type=flux_default
    else if(any(shape(flux_type)/=[ndir,nwflux])) then
      call mpistop('bfield physics: flux_type has wrong shape')
    end if

    phys_to_conserved=>bfield_to_conserved
    phys_to_primitive=>bfield_to_primitive
    phys_get_cmax=>bfield_get_cmax
    phys_get_cbounds=>bfield_get_cbounds
    phys_get_flux=>bfield_get_flux
    phys_get_dt=>bfield_get_dt
    phys_add_source_geom=>bfield_add_source_geom
    phys_check_params=>bfield_check_params

    ! Keep the same normalisation as the MHD PotentialField case so that a
    ! boundary file in Gauss maps to identical AMRVAC magnetic values.
    call bfield_physical_units()
  end subroutine bfield_phys_init

  subroutine bfield_check_params
    use mod_global_parameters
    use mod_physics, only: physics_type
    use mod_geometry, only: coordinate

    if(mype==0) then
      write(*,*) '====B-only field construction settings==============='
      write(*,*) 'physics_type=',trim(physics_type)
      write(*,*) 'Dimensionality=',ndim,' vector components=',ndir
      write(*,*) 'number of variables nw=',nw
      write(*,*) 'start index iwstart=',iwstart
      write(*,*) 'number of variables with BCs=',nwgc
      write(*,*) 'number of variables with fluxes=',nwflux
      write(*,*) 'coordinate set to type,slab:',coordinate,slab
      write(*,*) 'unit_magneticfield=',unit_magneticfield
      write(*,*) '======================================================='
    end if
  end subroutine bfield_check_params

  subroutine bfield_to_conserved(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(inout) :: w(ixI^S,nw)
    double precision, intent(in) :: x(ixI^S,1:^ND)
  end subroutine bfield_to_conserved

  subroutine bfield_to_primitive(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(inout) :: w(ixI^S,nw)
    double precision, intent(in) :: x(ixI^S,1:^ND)
  end subroutine bfield_to_primitive

  subroutine bfield_get_cmax(w,x,ixI^L,ixO^L,idim,cmax)
    use mod_global_parameters
    integer, intent(in) :: ixI^L,ixO^L,idim
    double precision, intent(in) :: w(ixI^S,nw),x(ixI^S,1:^ND)
    double precision, intent(inout) :: cmax(ixI^S)

    cmax(ixO^S)=zero
  end subroutine bfield_get_cmax

  subroutine bfield_get_cbounds(wLC,wRC,wLp,wRp,x,ixI^L,ixO^L,idim,&
       Hspeed,cmax,cmin)
    use mod_global_parameters
    integer, intent(in) :: ixI^L,ixO^L,idim
    double precision, intent(in) :: wLC(ixI^S,nw),wRC(ixI^S,nw)
    double precision, intent(in) :: wLp(ixI^S,nw),wRp(ixI^S,nw)
    double precision, intent(in) :: x(ixI^S,1:^ND)
    double precision, intent(in) :: Hspeed(ixI^S,1:number_species)
    double precision, intent(inout) :: cmax(ixI^S,1:number_species)
    double precision, intent(inout), optional :: cmin(ixI^S,1:number_species)

    cmax(ixO^S,1)=zero
    if(present(cmin)) cmin(ixO^S,1)=zero
  end subroutine bfield_get_cbounds

  subroutine bfield_get_flux(wC,w,x,ixI^L,ixO^L,idim,f)
    use mod_global_parameters
    integer, intent(in) :: ixI^L,ixO^L,idim
    double precision, intent(in) :: wC(ixI^S,nw),w(ixI^S,nw)
    double precision, intent(in) :: x(ixI^S,1:^ND)
    double precision, intent(out) :: f(ixI^S,nwflux)

    f(ixO^S,1:nwflux)=zero
  end subroutine bfield_get_flux

  subroutine bfield_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
    use mod_global_parameters
    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: dx^D,x(ixI^S,1:^ND)
    double precision, intent(in) :: w(ixI^S,1:nw)
    double precision, intent(inout) :: dtnew

    dtnew=bigdouble
  end subroutine bfield_get_dt

  subroutine bfield_add_source_geom(qdt,dtfactor,ixI^L,ixO^L,wCT,wprim,w,x)
    use mod_global_parameters
    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: qdt,dtfactor,x(ixI^S,1:^ND)
    double precision, intent(inout) :: wCT(ixI^S,1:nw),wprim(ixI^S,1:nw),&
       w(ixI^S,1:nw)
  end subroutine bfield_add_source_geom

  subroutine bfield_physical_units()
    use mod_global_parameters
    use mod_eos_container, only: eos

    double precision :: mp,kB,miu0,c_lightspeed
    double precision :: a,b,he
    character(len=std_len) :: eos_type

    if(SI_unit) then
      mp=mp_SI
      kB=kB_SI
      miu0=miu0_SI
      const_sigmaSB=sigma_SB_SI
      c_lightspeed=c_SI
    else
      mp=mp_cgs
      kB=kB_cgs
      miu0=4.d0*dpi
      const_sigmaSB=sigma_SB_cgs
      c_lightspeed=const_c
    end if

    he=0.1d0
    eos_type='FI'
    if(allocated(eos)) then
      he=eos%He_abundance
      eos_type=eos%eos_type
    end if

    if(eos_type=='LTE') then
      a=1.d0
      b=1.d0
    else
      a=1.d0+4.d0*he
      if(eos_type=='PI') then
        b=1.d0+bfield_h_ion_fr+he*(bfield_he_ion_fr*&
           (bfield_he_ion_fr2+1.d0)+1.d0)
      else
        b=2.d0+3.d0*he
      end if
    end if

    if(unit_density/=1.d0 .or. unit_numberdensity/=1.d0) then
      if(unit_density/=1.d0) then
        unit_numberdensity=unit_density/(a*mp)
      else if(unit_numberdensity/=1.d0) then
        unit_density=a*mp*unit_numberdensity
      end if
      if(unit_temperature/=1.d0) then
        unit_pressure=b*unit_numberdensity*kB*unit_temperature
        unit_velocity=sqrt(unit_pressure/unit_density)
        unit_magneticfield=sqrt(miu0*unit_pressure)
        if(unit_length/=1.d0) then
          unit_time=unit_length/unit_velocity
        else if(unit_time/=1.d0) then
          unit_length=unit_velocity*unit_time
        end if
      else if(unit_magneticfield/=1.d0) then
        unit_pressure=unit_magneticfield**2/miu0
        unit_temperature=unit_pressure/(b*unit_numberdensity*kB)
        unit_velocity=sqrt(unit_pressure/unit_density)
        if(unit_length/=1.d0) then
          unit_time=unit_length/unit_velocity
        else if(unit_time/=1.d0) then
          unit_length=unit_velocity*unit_time
        end if
      else if(unit_pressure/=1.d0) then
        unit_temperature=unit_pressure/(b*unit_numberdensity*kB)
        unit_velocity=sqrt(unit_pressure/unit_density)
        unit_magneticfield=sqrt(miu0*unit_pressure)
        if(unit_length/=1.d0) then
          unit_time=unit_length/unit_velocity
        else if(unit_time/=1.d0) then
          unit_length=unit_velocity*unit_time
        end if
      else if(unit_velocity/=1.d0) then
        unit_pressure=unit_density*unit_velocity**2
        unit_temperature=unit_pressure/(b*unit_numberdensity*kB)
        unit_magneticfield=sqrt(miu0*unit_pressure)
        if(unit_length/=1.d0) then
          unit_time=unit_length/unit_velocity
        else if(unit_time/=1.d0) then
          unit_length=unit_velocity*unit_time
        end if
      else if(unit_time/=1.d0) then
        unit_velocity=unit_length/unit_time
        unit_pressure=unit_density*unit_velocity**2
        unit_magneticfield=sqrt(miu0*unit_pressure)
        unit_temperature=unit_pressure/(b*unit_numberdensity*kB)
      end if
    else if(unit_temperature/=1.d0) then
      if(unit_magneticfield/=1.d0) then
        unit_pressure=unit_magneticfield**2/miu0
        unit_numberdensity=unit_pressure/(b*unit_temperature*kB)
        unit_density=a*mp*unit_numberdensity
        unit_velocity=sqrt(unit_pressure/unit_density)
        if(unit_length/=1.d0) then
          unit_time=unit_length/unit_velocity
        else if(unit_time/=1.d0) then
          unit_length=unit_velocity*unit_time
        end if
      else if(unit_pressure/=1.d0) then
        unit_magneticfield=sqrt(miu0*unit_pressure)
        unit_numberdensity=unit_pressure/(b*unit_temperature*kB)
        unit_density=a*mp*unit_numberdensity
        unit_velocity=sqrt(unit_pressure/unit_density)
        if(unit_length/=1.d0) then
          unit_time=unit_length/unit_velocity
        else if(unit_time/=1.d0) then
          unit_length=unit_velocity*unit_time
        end if
      end if
    else if(unit_magneticfield/=1.d0) then
      if(unit_velocity/=1.d0) then
        unit_pressure=unit_magneticfield**2/miu0
        unit_numberdensity=unit_pressure/(b*unit_temperature*kB)
        unit_density=a*mp*unit_numberdensity
        unit_temperature=unit_pressure/(b*unit_numberdensity*kB)
        if(unit_length/=1.d0) then
          unit_time=unit_length/unit_velocity
        else if(unit_time/=1.d0) then
          unit_length=unit_velocity*unit_time
        end if
      else if(unit_time/=0.d0) then
        unit_pressure=unit_magneticfield**2/miu0
        unit_velocity=unit_length/unit_time
        unit_density=unit_pressure/unit_velocity**2
        unit_numberdensity=unit_density/(a*mp)
        unit_temperature=unit_pressure/(b*unit_numberdensity*kB)
      end if
    else if(unit_pressure/=1.d0) then
      if(unit_velocity/=1.d0) then
        unit_magneticfield=sqrt(miu0*unit_pressure)
        unit_density=unit_pressure/unit_velocity**2
        unit_numberdensity=unit_density/(a*mp)
        unit_temperature=unit_pressure/(b*unit_numberdensity*kB)
        if(unit_length/=1.d0) then
          unit_time=unit_length/unit_velocity
        else if(unit_time/=1.d0) then
          unit_length=unit_velocity*unit_time
        end if
      else if(unit_time/=0.d0) then
        unit_magneticfield=sqrt(miu0*unit_pressure)
        unit_velocity=unit_length/unit_time
        unit_density=unit_pressure/unit_velocity**2
        unit_numberdensity=unit_density/(a*mp)
        unit_temperature=unit_pressure/(b*unit_numberdensity*kB)
      end if
    end if

    c_norm=c_lightspeed/unit_velocity
    unit_charge=unit_magneticfield*unit_length**2/unit_velocity/miu0
    if(.not.SI_unit) unit_charge=unit_charge*const_c
    unit_mass=unit_density*unit_length**3
  end subroutine bfield_physical_units

end module mod_bfield
