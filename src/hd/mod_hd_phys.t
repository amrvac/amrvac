!> Hydrodynamics module
module mod_hd_phys

  implicit none
  private

  !> Whether an energy equation is used
  logical, public, protected              :: hd_energy = .true.

  !> Whether thermal conduction is added
  logical, public, protected              :: hd_thermal_conduction = .false.

  !> Whether radiative cooling is added
  logical, public, protected              :: hd_radiative_cooling = .false.

  !> Whether dust is added
  logical, public, protected              :: hd_dust = .false.

  !> Whether viscosity is added
  logical, public, protected              :: hd_viscosity = .false.

  !> Whether gravity is added
  logical, public, protected              :: hd_gravity = .false.

  !> Whether particles module is added
  logical, public, protected              :: hd_particles = .false.

  !> Number of tracer species
  integer, public, protected              :: hd_n_tracer = 0

  !> Index of the density (in the w array)
  integer, public, protected              :: rho_

  !> Indices of the momentum density
  integer, allocatable, public, protected :: mom(:)

  !> Indices of the tracers
  integer, allocatable, public, protected :: tracer(:)

  !> Index of the energy density (-1 if not present)
  integer, public, protected              :: e_

  !> Index of the gas pressure (-1 if not present) should equal e_
  integer, public, protected              :: p_

  !> The adiabatic index
  double precision, public                :: hd_gamma = 5.d0/3.0d0

  !> The adiabatic constant
  double precision, public                :: hd_adiab = 1.0d0

  !> The smallest allowed energy
  double precision, protected             :: smalle

  !> The smallest allowed density
  double precision, protected             :: minrho

  !> The smallest allowed pressure
  double precision, protected             :: minp

  !> Helium abundance over Hydrogen
  double precision, public, protected  :: He_abundance=0.1d0

  ! Public methods
  public :: hd_phys_init
  public :: hd_kin_en
  public :: hd_get_pthermal
  public :: hd_get_v
  public :: hd_to_conserved
  public :: hd_to_primitive

contains

  !> Read this module's parameters from a file
  subroutine hd_read_params(files)
    use mod_global_parameters
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /hd_list/ hd_energy, hd_n_tracer, hd_gamma, hd_adiab, &
    hd_dust, hd_thermal_conduction, hd_radiative_cooling, hd_viscosity, &
    hd_gravity, He_abundance, SI_unit, hd_particles

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, hd_list, end=111)
111    close(unitpar)
    end do

  end subroutine hd_read_params

  !> Write this module's parameters to a snapsoht
  subroutine hd_write_info(fh)
    use mod_global_parameters
    integer, intent(in)                 :: fh
    integer, parameter                  :: n_par = 1
    double precision                    :: values(n_par)
    character(len=name_len)             :: names(n_par)
    integer, dimension(MPI_STATUS_SIZE) :: st
    integer                             :: er

    call MPI_FILE_WRITE(fh, n_par, 1, MPI_INTEGER, st, er)

    names(1) = "gamma"
    values(1) = hd_gamma
    call MPI_FILE_WRITE(fh, values, n_par, MPI_DOUBLE_PRECISION, st, er)
    call MPI_FILE_WRITE(fh, names, n_par * name_len, MPI_CHARACTER, st, er)
  end subroutine hd_write_info

  !> Initialize the module
  subroutine hd_phys_init()
    use mod_global_parameters
    use mod_thermal_conduction
    use mod_radiative_cooling
    use mod_dust, only: dust_init
    use mod_viscosity, only: viscosity_init
    use mod_gravity, only: gravity_init
    use mod_particles, only: particles_init
    use mod_physics

    integer :: itr, idir

    call hd_read_params(par_files)

    physics_type = "hd"
    phys_energy  = hd_energy
    use_particles = hd_particles

    ! Determine flux variables
    rho_ = var_set_rho()

    allocate(mom(ndir))
    mom(:) = var_set_momentum(ndir)

    ! Set index of energy variable
    if (hd_energy) then
       e_ = var_set_energy()
       p_ = e_
    else
       e_ = -1
       p_ = -1
    end if

    allocate(tracer(hd_n_tracer))

    ! Set starting index of tracers
    do itr = 1, hd_n_tracer
       tracer(itr) = var_set_fluxvar("trc", "trp", itr)
    end do

    phys_get_dt          => hd_get_dt
    phys_get_cmax        => hd_get_cmax
    phys_get_flux        => hd_get_flux
    phys_add_source_geom => hd_add_source_geom
    phys_add_source      => hd_add_source
    phys_to_conserved    => hd_to_conserved
    phys_to_primitive    => hd_to_primitive
    phys_check_params    => hd_check_params
    phys_check_w         => hd_check_w
    phys_get_pthermal    => hd_get_pthermal
    phys_write_info      => hd_write_info

    ! derive units from basic units
    call hd_physical_units()

    if (hd_dust) call dust_init(rho_, mom(:), e_)

    ! initialize thermal conduction module
    if (hd_thermal_conduction) then
      if (.not. hd_energy) &
           call mpistop("thermal conduction needs hd_energy=T")
      call thermal_conduction_init(hd_gamma)
    end if

    ! Initialize radiative cooling module
    if (hd_radiative_cooling) then
      if (.not. hd_energy) &
           call mpistop("radiative cooling needs hd_energy=T")
      call radiative_cooling_init(hd_gamma,He_abundance)
    end if

    ! Initialize viscosity module
    if (hd_viscosity) call viscosity_init()

    ! Initialize gravity module
    if (hd_gravity) call gravity_init()

    ! Initialize particles module
    if (hd_particles) call particles_init()

    ! Check whether custom flux types have been defined
    if (.not. allocated(flux_type)) then
       allocate(flux_type(ndir, nw))
       flux_type = flux_default
    else if (any(shape(flux_type) /= [ndir, nw])) then
       call mpistop("phys_check error: flux_type has wrong shape")
    end if

    nvector      = 1 ! No. vector vars
    allocate(iw_vector(nvector))
    iw_vector(1) = mom(1) - 1

  end subroutine hd_phys_init

  subroutine hd_check_params
    use mod_global_parameters
    use mod_dust, only: dust_check_params

    minrho = max(0.0d0, small_density)

    if (.not. hd_energy) then
       if (hd_gamma <= 0.0d0) call mpistop ("Error: hd_gamma <= 0")
       if (hd_adiab <= 0.0d0) call mpistop ("Error: hd_adiab <= 0")
       minp   = hd_adiab*minrho**hd_gamma
    else
       if (hd_gamma <= 0.0d0 .or. hd_gamma == 1.0d0) &
            call mpistop ("Error: hd_gamma <= 0 or hd_gamma == 1.0")
       minp   = max(0.0d0, small_pressure)
       smalle = minp/(hd_gamma - 1.0d0)
    end if

    if (hd_dust) call dust_check_params()

  end subroutine hd_check_params

  subroutine hd_physical_units
    use mod_global_parameters
    double precision :: mp,kB
    ! Derive scaling units
    if(SI_unit) then
      mp=mp_SI
      kB=kB_SI
    else
      mp=mp_cgs
      kB=kB_cgs
    end if
    if(unit_velocity==0) then
      unit_density=(1.d0+4.d0*He_abundance)*mp*unit_numberdensity
      unit_pressure=(2.d0+3.d0*He_abundance)*unit_numberdensity*kB*unit_temperature
      unit_velocity=dsqrt(unit_pressure/unit_density)
      unit_time=unit_length/unit_velocity
    else
      unit_density=(1.d0+4.d0*He_abundance)*mp*unit_numberdensity
      unit_pressure=unit_density*unit_velocity**2
      unit_temperature=unit_pressure/((2.d0+3.d0*He_abundance)*unit_numberdensity*kB)
      unit_time=unit_length/unit_velocity
    end if

  end subroutine hd_physical_units

  !> Returns 0 in argument flag where values are ok
  subroutine hd_check_w(primitive, ixI^L, ixO^L, w, flag)
    use mod_global_parameters

    logical, intent(in)          :: primitive
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, nw)
    integer, intent(inout)       :: flag(ixI^S)
    double precision             :: tmp(ixI^S)

    flag(ixO^S) = 0
    where(w(ixO^S, rho_) < minrho) flag(ixO^S) = rho_

    if (hd_energy) then
       if (primitive) then
          where(w(ixO^S, e_) < minp) flag(ixO^S) = e_
       else
          tmp(ixO^S) = (hd_gamma - 1.0d0)*(w(ixO^S, e_) - &
               hd_kin_en(w, ixI^L, ixO^L))
          where(tmp(ixO^S) < minp) flag(ixO^S) = e_
       endif
    end if

  end subroutine hd_check_w

  !> Transform primitive variables into conservative ones
  subroutine hd_to_conserved(ixI^L, ixO^L, w, x)
    use mod_global_parameters
    use mod_dust, only: dust_to_conserved
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    double precision                :: invgam
    integer                         :: idir, itr

    ! Convert velocity to momentum
    do idir = 1, ndir
       w(ixO^S, mom(idir)) = w(ixO^S, rho_) * w(ixO^S, mom(idir))
    end do

    if (hd_energy) then
       invgam = 1.d0/(hd_gamma - 1.0d0)
       ! Calculate total energy from pressure and kinetic energy
       w(ixO^S, e_) = w(ixO^S, e_) * invgam + hd_kin_en(w, ixI^L, ixO^L)
    end if

    do itr = 1, hd_n_tracer
       w(ixO^S, tracer(itr)) = w(ixO^S, rho_) * w(ixO^S, tracer(itr))
    end do

    call dust_to_conserved(ixI^L, ixO^L, w, x)

    call handle_small_values(.false., w, x, ixI^L, ixO^L, 'hd_to_conserved')

  end subroutine hd_to_conserved

  !> Transform conservative variables into primitive ones
  subroutine hd_to_primitive(ixI^L, ixO^L, w, x)
    use mod_global_parameters
    use mod_dust, only: dust_to_primitive
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    integer                         :: itr, idir

    if (hd_energy) then
       ! Compute pressure
       w(ixO^S, e_) = (hd_gamma - 1.0d0) * (w(ixO^S, e_) - &
            hd_kin_en(w, ixI^L, ixO^L))
    end if

    ! Convert momentum to velocity
    do idir = 1, ndir
       w(ixO^S, mom(idir)) = w(ixO^S, mom(idir)) * hd_inv_rho(w, ixI^L, ixO^L)
    end do

    do itr = 1, hd_n_tracer
       w(ixO^S, tracer(itr)) = w(ixO^S, tracer(itr)) * hd_inv_rho(w, ixI^L, ixO^L)
    end do

    ! Convert dust momentum to dust velocity
    call dust_to_primitive(ixI^L, ixO^L, w, x)

    call handle_small_values(.true., w, x, ixI^L, ixO^L, 'hd_to_primitive')

  end subroutine hd_to_primitive

  subroutine e_to_rhos(ixI^L, ixO^L, w, x)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L
    double precision             :: w(ixI^S, nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)

    if (hd_energy) then
       w(ixO^S, e_) = (hd_gamma - 1.0d0) * w(ixO^S, rho_)**(1.0d0 - hd_gamma) * &
            (w(ixO^S, e_) - hd_kin_en(w, ixI^L, ixO^L))
    else
       call mpistop("energy from entropy can not be used with -eos = iso !")
    end if
  end subroutine e_to_rhos

  subroutine rhos_to_e(ixI^L, ixO^L, w, x)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L
    double precision             :: w(ixI^S, nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)

    if (hd_energy) then
       w(ixO^S, e_) = w(ixO^S, rho_)**(hd_gamma - 1.0d0) * w(ixO^S, e_) &
            / (hd_gamma - 1.0d0) + hd_kin_en(w, ixI^L, ixO^L)
    else
       call mpistop("entropy from energy can not be used with -eos = iso !")
    end if
  end subroutine rhos_to_e

  !> Calculate v_i = m_i / rho within ixO^L
  subroutine hd_get_v(w, x, ixI^L, ixO^L, idim, v)
    use mod_global_parameters
    integer, intent(in)           :: ixI^L, ixO^L, idim
    double precision, intent(in)  :: w(ixI^S, nw), x(ixI^S, 1:ndim)
    double precision, intent(out) :: v(ixI^S)

    v(ixO^S) = w(ixO^S, mom(idim)) * hd_inv_rho(w, ixI^L, ixO^L)
  end subroutine hd_get_v

  !> Calculate cmax_idim = csound + abs(v_idim) within ixO^L
  subroutine hd_get_cmax(w, x, ixI^L, ixO^L, idim, cmax, cmin)
    use mod_global_parameters
    use mod_dust, only: dust_get_cmax

    integer, intent(in)                       :: ixI^L, ixO^L, idim
    double precision, intent(in)              :: w(ixI^S, nw), x(ixI^S, 1:ndim)
    double precision, intent(inout)           :: cmax(ixI^S)
    double precision, intent(inout), optional :: cmin(ixI^S)
    double precision                          :: csound(ixI^S)
    double precision                          :: v(ixI^S)

    call hd_get_v(w, x, ixI^L, ixO^L, idim, v)
    call hd_get_pthermal(w, x, ixI^L, ixO^L, csound)

    csound(ixO^S) = sqrt(hd_gamma * csound(ixO^S) * &
         hd_inv_rho(w, ixI^L, ixO^L))

    if (present(cmin)) then
       cmax(ixO^S) = max(v(ixO^S)+csound(ixO^S), zero)
       cmin(ixO^S) = min(v(ixO^S)-csound(ixO^S), zero)
    else
       cmax(ixO^S) = abs(v(ixO^S))+csound(ixO^S)
    end if

    if (hd_dust) &
         call dust_get_cmax(w, x, ixI^L, ixO^L, idim, cmax, cmin)
  end subroutine hd_get_cmax

  !> Calculate thermal pressure=(gamma-1)*(e-0.5*m**2/rho) within ixO^L
  subroutine hd_get_pthermal(w, x, ixI^L, ixO^L, pth)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: pth(ixI^S)

    if (hd_energy) then
       pth(ixO^S) = (hd_gamma - 1.0d0) * (w(ixO^S, e_) - &
            hd_kin_en(w, ixI^L, ixO^L))
    else
       pth(ixO^S) = hd_adiab * w(ixO^S, rho_)**hd_gamma
    end if

  end subroutine hd_get_pthermal

  ! Calculate non-transport flux f_idim[iw] within ixO^L.
  subroutine hd_get_flux(w, x, ixI^L, ixO^L, idim, f)
    use mod_global_parameters
    use mod_dust, only: dust_get_flux

    integer, intent(in)             :: ixI^L, ixO^L, idim
    double precision, intent(in)    :: w(ixI^S, 1:nw), x(ixI^S, 1:ndim)
    double precision, intent(out)   :: f(ixI^S, nwflux)
    double precision                :: pth(ixI^S), v(ixI^S)
    integer                         :: idir, itr

    call hd_get_pthermal(w, x, ixI^L, ixO^L, pth)
    call hd_get_v(w, x, ixI^L, ixO^L, idim, v)

    f(ixO^S, rho_) = v(ixO^S) * w(ixO^S, rho_)

    ! Momentum flux is v_i*m_i, +p in direction idim
    do idir = 1, ndir
      f(ixO^S, mom(idir)) = v(ixO^S) * w(ixO^S, mom(idir))
    end do

    f(ixO^S, mom(idim)) = f(ixO^S, mom(idim)) + pth(ixO^S)

    if(hd_energy) then
      ! Energy flux is v_i*e + v*p ! Check? m_i/rho*p
      f(ixO^S, e_) = v(ixO^S) * (w(ixO^S, e_) + pth(ixO^S))
    end if

    do itr = 1, hd_n_tracer
       f(ixO^S, tracer(itr)) = v(ixO^S) * w(ixO^S, tracer(itr))
    end do

    ! Dust fluxes
    call dust_get_flux(w, x, ixI^L, ixO^L, idim, f)

  end subroutine hd_get_flux

  !> Add geometrical source terms to w
  !>
  !> Notice that the expressions of the geometrical terms depend only on ndir,
  !> not ndim. Eg, they are the same in 2.5D and in 3D, for any geometry.
  !>
  !> Ileyk : to do :
  !>     - give the possibility to set angmomfix=.true.
  !>     - address the source term for the dust
  subroutine hd_add_source_geom(qdt, ixI^L, ixO^L, wCT, w, x)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, x(ixI^S, 1:ndim)
    double precision, intent(inout) :: wCT(ixI^S, 1:nw), w(ixI^S, 1:nw)
    ! to change and to set as a parameter in the parfile once the possibility to
    ! solve the equations in an angular momentum conserving form has been
    ! implemented (change tvdlf.t eg)
    double precision :: tmp(ixI^S),tmp1(ixI^S)
    integer                         :: iw,idir, h1x^L{^NOONED, h2x^L}
    integer :: mr_,mphi_ ! Polar var. names
    logical                         :: angmomfix = .false.

    mr_=mom(1); mphi_=mom(1)-1+phi_ ! Polar var. names

    select case (typeaxial)
    case ("cylindrical")
       ! s[mr]=(pthermal+mphi**2/rho)/radius
       call hd_get_pthermal(wCT,x,ixI^L,ixO^L,tmp)
       if(phi_>0) then
         tmp(ixO^S)=tmp(ixO^S)+wCT(ixO^S,mphi_)**2/wCT(ixO^S,rho_)
         w(ixO^S,mr_)=w(ixO^S,mr_)+qdt*tmp(ixO^S)/x(ixO^S,1)
         ! s[mphi]=(-mphi*mr/rho)/radius
         ! Ileyk : beware the index permutation : mphi=2 if -phi=2 (2.5D
         ! (r,theta) grids) BUT mphi=3 if -phi=3 (for 2.5D (r,z) grids)
         if(.not. angmomfix) then
           tmp(ixO^S)=-wCT(ixO^S,mphi_)*wCT(ixO^S,mr_)/wCT(ixO^S,rho_)
         else
           tmp(ixO^S)=0.d0
         end if
         ! no geometrical source term if angular momentum conserving form of
         ! the equations
         w(ixO^S,mphi_)=w(ixO^S,mphi_)+qdt*tmp(ixO^S)/x(ixO^S,1)
       else
         ! s[mr]=2pthermal/radius
         w(ixO^S,mr_)=w(ixO^S,mr_)+qdt*tmp(ixO^S)/x(ixO^S,1)
       end if
    case ("spherical")
       h1x^L=ixO^L-kr(1,^D); {^NOONED h2x^L=ixO^L-kr(2,^D);}
       ! s[mr]=((mtheta**2+mphi**2)/rho+2*p)/r
       call hd_get_pthermal(wCT,x,ixI^L,ixO^L,tmp1)
       tmp(ixO^S)=tmp1(ixO^S)*x(ixO^S,1) &
            *(block%surfaceC1(ixO^S)-block%surfaceC1(h1x^S)) &
            /block%dvolume(ixO^S)
       if(ndir>1) then
         do idir=2,ndir
           tmp(ixO^S)=tmp(ixO^S)+wCT(ixO^S,mom(idir))**2/wCT(ixO^S,rho_)
         end do
       end if
       w(ixO^S,mr_)=w(ixO^S,mr_)+qdt*tmp(ixO^S)/x(ixO^S,1)

       {^NOONED
       ! s[mtheta]=-(mr*mtheta/rho)/r+cot(theta)*(mphi**2/rho+p)/r
       tmp(ixO^S)=tmp1(ixO^S)*x(ixO^S,1) &
            *(block%surfaceC2(ixO^S)-block%surfaceC2(h2x^S)) &
            /block%dvolume(ixO^S)
       if(ndir==3) tmp(ixO^S)=tmp(ixO^S)+(wCT(ixO^S,mom(3))**2/wCT(ixO^S,rho_))/tan(x(ixO^S,2))
       if (.not. angmomfix) tmp(ixO^S)=tmp(ixO^S)-(wCT(ixO^S,mom(2))*wCT(ixO^S,mr_))/wCT(ixO^S,rho_)
       w(ixO^S,mom(2))=w(ixO^S,mom(2))+qdt*tmp(ixO^S)/x(ixO^S,1)

       if(ndir==3) then
         ! s[mphi]=-(mphi*mr/rho)/r-cot(theta)*(mtheta*mphi/rho)/r
         if(.not. angmomfix) then 
           tmp(ixO^S)=-(wCT(ixO^S,mom(3))*wCT(ixO^S,mr_))/wCT(ixO^S,rho_)&
                      -(wCT(ixO^S,mom(2))*wCT(ixO^S,mom(3)))/wCT(ixO^S,rho_)/tan(x(ixO^S,2))
           w(ixO^S,mom(3))=w(ixO^S,mom(3))+qdt*tmp(ixO^S)/x(ixO^S,1)
         end if
       end if
       }
    end select

  end subroutine hd_add_source_geom

  ! w[iw]= w[iw]+qdt*S[wCT, qtC, x] where S is the source based on wCT within ixO
  subroutine hd_add_source(qdt,ixI^L,ixO^L,wCT,w,x,qsourcesplit)
    use mod_global_parameters
    use mod_radiative_cooling, only: radiative_cooling_add_source
    use mod_dust, only: dust_add_source
    use mod_viscosity, only: viscosity_add_source
    use mod_gravity, only: gravity_add_source

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixI^S, 1:nw), x(ixI^S, 1:ndim)
    double precision, intent(inout) :: w(ixI^S, 1:nw)
    logical, intent(in)             :: qsourcesplit

    if(hd_dust) then
      call dust_add_source(qdt,ixI^L,ixO^L,wCT,w,x,qsourcesplit)
    end if

    if(hd_radiative_cooling) then
      call radiative_cooling_add_source(qdt,ixI^L,ixO^L,wCT,w,x,qsourcesplit)
    end if

    if(hd_viscosity) then
      call viscosity_add_source(qdt,ixI^L,ixO^L,wCT,w,x,hd_energy,qsourcesplit)
    end if

    if(hd_gravity) then
      call gravity_add_source(qdt,ixI^L,ixO^L,wCT,w,x,hd_energy,qsourcesplit)
    end if

  end subroutine hd_add_source

  subroutine hd_get_dt(w, ixI^L, ixO^L, dtnew, dx^D, x)
    use mod_global_parameters
    use mod_dust, only: dust_get_dt
    use mod_radiative_cooling, only: cooling_get_dt
    use mod_viscosity, only: viscosity_get_dt 
    use mod_gravity, only: gravity_get_dt

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: dx^D, x(ixI^S, 1:^ND)
    double precision, intent(in)    :: w(ixI^S, 1:nw)
    double precision, intent(inout) :: dtnew

    dtnew = bigdouble

    if(hd_dust) then
      call dust_get_dt(w, ixI^L, ixO^L, dtnew, dx^D, x)
    end if

    if(hd_radiative_cooling) then
      call cooling_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
    end if

    if(hd_viscosity) then
      call viscosity_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
    end if

    if(hd_gravity) then
      call gravity_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
    end if

  end subroutine hd_get_dt

  function hd_kin_en(w, ixI^L, ixO^L) result(ke)
    use mod_global_parameters, only: nw, ndim
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S, nw)
    double precision              :: ke(ixO^S)

    ke = 0.5d0 * sum(w(ixO^S, mom(:))**2, dim=ndim+1) * &
         hd_inv_rho(w, ixI^L, ixO^L)
  end function hd_kin_en

  function hd_inv_rho(w, ixI^L, ixO^L) result(inv_rho)
    use mod_global_parameters, only: nw, ndim
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S, nw)
    double precision              :: inv_rho(ixO^S)

    ! Can make this more robust
    inv_rho = 1.0d0 / w(ixO^S, rho_)
  end function hd_inv_rho

  subroutine handle_small_values(primitive, w, x, ixI^L, ixO^L, subname)
    use mod_global_parameters
    use mod_small_values
    logical, intent(in)             :: primitive
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    character(len=*), intent(in)    :: subname

    double precision :: smallone
    integer :: idir, flag(ixI^S)

    call hd_check_w(primitive, ixI^L, ixO^L, w, flag)

    if (any(flag(ixO^S) /= 0)) then
       select case (small_values_method)
       case ("replace")
          where(flag(ixO^S) /= 0) w(ixO^S,rho_) = minrho

          do idir = 1, ndir
             where(flag(ixO^S) /= 0) w(ixO^S, mom(idir)) = 0.0d0
          end do

          if (hd_energy) then
             if(primitive) then
               smallone = minp
             else
               smallone = smalle
             end if
             where(flag(ixO^S) /= 0) w(ixO^S,e_) = smallone
          end if
       case ("average")
          call small_values_average(ixI^L, ixO^L, w, x, flag)
       case default
          call small_values_error(w, x, ixI^L, ixO^L, flag, subname)
       end select
    end if
  end subroutine handle_small_values

end module mod_hd_phys
