!> Magneto-hydrodynamics module
module mod_mhd_phys
  use mod_global_parameters, only: std_len
  implicit none
  private

  !> Whether an energy equation is used
  logical, public, protected              :: mhd_energy = .true.

  !> Whether thermal conduction is used
  logical, public, protected              :: mhd_thermal_conduction = .false.

  !> Whether radiative cooling is added
  logical, public, protected              :: mhd_radiative_cooling = .false.

  !> Whether viscosity is added
  logical, public, protected              :: mhd_viscosity= .false.

  !> Whether gravity is added
  logical, public, protected              :: mhd_gravity= .false.

  !> Whether Hall-MHD is used
  logical, public, protected              :: mhd_Hall = .false.

  !> Whether particles module is added
  logical, public, protected              :: mhd_particles = .false.

  !> Whether MHD-GLM is used
  logical, public, protected              :: mhd_glm = .false.

  !> TODO: describe and set value
  double precision, public                :: mhd_glm_Cr = -0.2d0

  !> MHD fourth order
  logical, public, protected              :: mhd_4th_order = .false.

  !> Number of tracer species
  integer, public, protected              :: mhd_n_tracer = 0

  !> Index of the density (in the w array)
  integer, public, protected              :: rho_

  !> Indices of the momentum density
  integer, allocatable, public, protected :: mom(:)

  !> Index of the energy density (-1 if not present)
  integer, public, protected              :: e_

  !> Index of the gas pressure (-1 if not present) should equal e_
  integer, public, protected              :: p_

  !> Indices of the magnetic field
  integer, allocatable, public, protected :: mag(:)

  !> Indices of the GLM psi
  integer, public, protected :: psi_

  !> Indices of the tracers
  integer, allocatable, public, protected :: tracer(:)

  !> The adiabatic index
  double precision, public                :: mhd_gamma = 5.d0/3.0d0

  !> The adiabatic constant
  double precision, public                :: mhd_adiab = 1.0d0

  !> The MHD resistivity
  double precision, public                :: mhd_eta = 0.0d0

  !> The MHD hyper-resistivity
  double precision, public                :: mhd_eta_hyper = 0.0d0

  !> TODO: what is this?
  double precision, public                :: mhd_etah = 0.0d0

  !> The smallest allowed energy
  double precision, protected             :: smalle

  !> The smallest allowed density
  double precision, protected             :: minrho

  !> The smallest allowed pressure
  double precision, protected             :: minp

  !> The number of waves
  integer :: nwwave=8

  !> Method type to clean divergence of B
  character(len=std_len) :: typedivbfix  = 'linde'

  !> Coefficient of diffusive divB cleaning
  double precision :: divbdiff     = 0.5d0

  !> Update all equations due to divB cleaning
  character(len=std_len) ::    typedivbdiff = 'all'

  !> Use a compact way to add resistivity
  logical :: compactres   = .false.

  !> Add divB wave in Roe solver
  logical, public :: divbwave     = .true.

  !> Helium abundance over Hydrogen
  double precision, public, protected  :: He_abundance=0.1d0

  !> To control divB=0 fix for boundary
  logical, public, protected :: boundary_divbfix(2*^ND)=.true.

  !> B0 field is force-free
  logical, public, protected :: B0field_forcefree=.true.
  ! Public methods
  public :: mhd_phys_init
  public :: mhd_kin_en
  public :: mhd_get_pthermal
  public :: mhd_get_v
  public :: mhd_to_conserved
  public :: mhd_to_primitive
  public :: mhd_get_csound2
  public :: get_divb
  public :: get_current

contains

  !> Read this module"s parameters from a file
  subroutine mhd_read_params(files)
    use mod_global_parameters
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /mhd_list/ mhd_energy, mhd_n_tracer, mhd_gamma, mhd_adiab,&
      mhd_eta, mhd_eta_hyper, mhd_etah, mhd_glm, mhd_glm_Cr, &
      mhd_thermal_conduction, mhd_radiative_cooling, mhd_Hall, mhd_gravity,&
      mhd_viscosity, mhd_4th_order, typedivbfix, divbdiff, typedivbdiff, compactres, &
      divbwave, He_abundance, SI_unit, B0field, B0field_forcefree, Bdip, Bquad, Boct, &
      Busr,mhd_particles, boundary_divbfix

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, mhd_list, end=111)
111    close(unitpar)
    end do

  end subroutine mhd_read_params

  !> Write this module's parameters to a snapsoht
  subroutine mhd_write_info(fh)
    use mod_global_parameters
    integer, intent(in)                 :: fh
    integer, parameter                  :: n_par = 1
    double precision                    :: values(n_par)
    character(len=name_len)             :: names(n_par)
    integer, dimension(MPI_STATUS_SIZE) :: st
    integer                             :: er

    call MPI_FILE_WRITE(fh, n_par, 1, MPI_INTEGER, st, er)

    names(1) = "gamma"
    values(1) = mhd_gamma
    call MPI_FILE_WRITE(fh, values, n_par, MPI_DOUBLE_PRECISION, st, er)
    call MPI_FILE_WRITE(fh, names, n_par * name_len, MPI_CHARACTER, st, er)
  end subroutine mhd_write_info

  subroutine mhd_phys_init()
    use mod_global_parameters
    use mod_thermal_conduction
    use mod_radiative_cooling
    use mod_viscosity, only: viscosity_init
    use mod_gravity, only: gravity_init
    use mod_particles, only: particles_init
    use mod_physics

    integer :: itr, idir

    call mhd_read_params(par_files)

    physics_type = "mhd"
    phys_energy=mhd_energy
    use_particles=mhd_particles

    ! Determine flux variables
    rho_ = var_set_rho()

    allocate(mom(ndir))
    mom(:) = var_set_momentum(ndir)

    ! Set index of energy variable
    if (mhd_energy) then
      nwwave = 8
      e_     = var_set_energy() ! energy density
      p_     = e_               ! gas pressure
    else
      nwwave = 7
      e_     = -1
      p_     = -1
    end if

    allocate(mag(ndir))
    mag(:) = var_set_bfield(ndir)

    if (mhd_glm) then
      psi_ = var_set_fluxvar('psi', 'psi')
    else
      psi_ = -1
    end if

    allocate(tracer(mhd_n_tracer))

    ! Set starting index of tracers
    do itr = 1, mhd_n_tracer
      tracer(itr) = var_set_fluxvar("trc", "trp", itr)
    end do

    nvector      = 2 ! No. vector vars
    allocate(iw_vector(nvector))
    iw_vector(1) = mom(1) - 1   ! TODO: why like this?
    iw_vector(2) = mag(1) - 1   ! TODO: why like this?

    ! Check whether custom flux types have been defined
    if (.not. allocated(flux_type)) then
       allocate(flux_type(ndir, nw))
       flux_type = flux_default
    else if (any(shape(flux_type) /= [ndir, nw])) then
       call mpistop("phys_check error: flux_type has wrong shape")
    end if
    do idir=1,ndir
      if(ndim>1) flux_type(idir,mag(idir))=flux_tvdlf
    end do
    if(mhd_glm .and. ndim>1) flux_type(:,psi_)=flux_tvdlf

    phys_get_dt          => mhd_get_dt
    phys_get_cmax        => mhd_get_cmax
    phys_get_flux        => mhd_get_flux
    phys_add_source_geom => mhd_add_source_geom
    phys_add_source      => mhd_add_source
    phys_to_conserved    => mhd_to_conserved
    phys_to_primitive    => mhd_to_primitive
    phys_check_params    => mhd_check_params
    phys_check_w         => mhd_check_w
    phys_get_pthermal    => mhd_get_pthermal
    phys_boundary_adjust => mhd_boundary_adjust
    phys_write_info      => mhd_write_info

    ! derive units from basic units
    call mhd_physical_units()

    if(.not. mhd_energy .and. mhd_thermal_conduction) then
      call mpistop("thermal conduction needs mhd_energy=T")
    end if
    if(.not. mhd_energy .and. mhd_radiative_cooling) then
      call mpistop("radiative cooling needs mhd_energy=T")
    end if

    ! initialize thermal conduction module
    if (mhd_thermal_conduction) then
       call thermal_conduction_init(mhd_gamma)
    end if

    ! Initialize radiative cooling module
    if (mhd_radiative_cooling) then
      call radiative_cooling_init(mhd_gamma,He_abundance)
    end if

    ! Initialize viscosity module
    if(mhd_viscosity) then
      call viscosity_init()
    end if

    ! Initialize gravity module
    if(mhd_gravity) then
      call gravity_init()
    end if

    ! Initialize particles module
    if(mhd_particles) then
      call particles_init()
    end if

    if (mhd_glm) then
       ! Solve the Riemann problem for the linear 2x2 system for normal
       ! B-field and GLM_Psi according to Dedner 2002:
       phys_modify_wLR => glmSolve
    end if

    ! For Hall, we need one more reconstructed layer since currents are computed
    ! in getflux: assuming one additional ghost layer (two for FOURTHORDER) was
    ! added in nghostcells.
    if (mhd_hall) then
       if (mhd_4th_order) then
          phys_wider_stencil = 2
       else
          phys_wider_stencil = 1
       end if
    end if

  end subroutine mhd_phys_init

  subroutine mhd_check_params
    use mod_global_parameters

    minrho = max(0.0d0, small_density)

    if (.not. mhd_energy) then
       if (mhd_gamma <= 0.0d0) call mpistop ("Error: mhd_gamma <= 0")
       if (mhd_adiab < 0.0d0) call mpistop ("Error: mhd_adiab < 0")
       minp   = mhd_adiab*minrho**mhd_gamma
    else
       if (mhd_gamma <= 0.0d0 .or. mhd_gamma == 1.0d0) &
            call mpistop ("Error: mhd_gamma <= 0 or mhd_gamma == 1")
       minp   = max(0.0d0, small_pressure)
       smalle = minp/(mhd_gamma - 1.0d0)
    end if

  end subroutine mhd_check_params

  subroutine mhd_physical_units()
    use mod_global_parameters
    double precision :: mp,kB,miu0
    ! Derive scaling units
    if(SI_unit) then
      mp=mp_SI
      kB=kB_SI
      miu0=miu0_SI
    else
      mp=mp_cgs
      kB=kB_cgs
      miu0=4.d0*dpi
    end if
    if(unit_velocity==0) then
      unit_density=(1.d0+4.d0*He_abundance)*mp*unit_numberdensity
      unit_pressure=(2.d0+3.d0*He_abundance)*unit_numberdensity*kB*unit_temperature
      unit_velocity=dsqrt(unit_pressure/unit_density)
      unit_magneticfield=dsqrt(miu0*unit_pressure)
      unit_time=unit_length/unit_velocity
    else
      unit_density=(1.d0+4.d0*He_abundance)*mp*unit_numberdensity
      unit_pressure=unit_density*unit_velocity**2
      unit_temperature=unit_pressure/((2.d0+3.d0*He_abundance)*unit_numberdensity*kB)
      unit_magneticfield=dsqrt(miu0*unit_pressure)
      unit_time=unit_length/unit_velocity
    end if

  end subroutine mhd_physical_units

  subroutine mhd_check_w(primitive,ixI^L,ixO^L,w,flag)
    use mod_global_parameters

    logical, intent(in) :: primitive
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S,nw)
    integer, intent(inout) :: flag(ixI^S)
    double precision :: tmp(ixI^S)

    flag(ixO^S)=0
    where(w(ixO^S, rho_) < minrho) flag(ixO^S) = rho_

    if (mhd_energy) then
      if (primitive)then
        where(w(ixO^S, e_) < minp) flag(ixO^S) = e_
      else
        ! Calculate pressure=(gamma-1)*(e-0.5*(2ek+2eb))
        tmp(ixO^S)=(mhd_gamma-1.d0)*(w(ixO^S,e_) - &
           mhd_kin_en(w,ixI^L,ixO^L)-mhd_mag_en(w,ixI^L,ixO^L))
        where(tmp(ixO^S) < minp) flag(ixO^S) = e_
      end if
    end if
  end subroutine mhd_check_w

  !> Transform primitive variables into conservative ones
  subroutine mhd_to_conserved(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    integer                         :: idir, itr

    ! Convert velocity to momentum
    do idir = 1, ndir
       w(ixO^S, mom(idir)) = w(ixO^S, rho_) * w(ixO^S, mom(idir))
    end do

    if (mhd_energy) then
       ! Calculate total energy from pressure, kinetic and magnetic energy
       w(ixO^S,e_)=w(ixO^S,p_)/(mhd_gamma-1.d0) + &
            mhd_kin_en(w, ixI^L, ixO^L) + mhd_mag_en(w, ixI^L, ixO^L)
    end if

    do itr = 1, mhd_n_tracer
       w(ixO^S, tracer(itr)) = w(ixO^S, rho_) * w(ixO^S, tracer(itr))
    end do

    call handle_small_values(.false., w, x, ixI^L, ixO^L,'mhd_to_conserved')
  end subroutine mhd_to_conserved

  !> Transform conservative variables into primitive ones
  subroutine mhd_to_primitive(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    integer                         :: itr, idir

    if (mhd_energy) then
       ! Calculate pressure = (gamma-1) * (e-0.5*(2ek+2eb))
       w(ixO^S, e_) = (mhd_gamma - 1.0d0) * (w(ixO^S, e_) &
            - mhd_kin_en(w, ixI^L, ixO^L) &
            - mhd_mag_en(w, ixI^L, ixO^L))
    end if

    ! Convert momentum to velocity
    do idir = 1, ndir
       w(ixO^S, mom(idir)) = w(ixO^S, mom(idir)) * mhd_inv_rho(w, ixI^L, ixO^L)
    end do

    do itr = 1, mhd_n_tracer
       w(ixO^S, tracer(itr)) = w(ixO^S, tracer(itr)) * mhd_inv_rho(w, ixI^L, ixO^L)
    end do

    call handle_small_values(.true., w, x, ixI^L, ixO^L,'mhd_to_primitive')
  end subroutine mhd_to_primitive

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

    call mhd_check_w(primitive, ixI^L, ixO^L, w, flag)

    if (any(flag(ixO^S) /= 0)) then
       select case (small_values_method)
       case ("replace")
          where(flag(ixO^S) /= 0) w(ixO^S,rho_) = minrho

          do idir = 1, ndir
             where(flag(ixO^S) /= 0) w(ixO^S, mom(idir)) = 0.0d0
          end do

          if (mhd_energy) then
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

  !> Convert energy to entropy
  subroutine e_to_rhos(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision,intent(inout)  :: w(ixI^S,nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)

    if (mhd_energy) then
      w(ixO^S, e_) = (mhd_gamma - 1.0d0) * w(ixO^S, rho_)**(1.0d0 - mhd_gamma) * &
            (w(ixO^S, e_) - mhd_kin_en(w, ixI^L, ixO^L) &
            - mhd_mag_en(w, ixI^L, ixO^L))
    else
      call mpistop("e_to_rhos can not be used without energy equation!")
    end if
  end subroutine e_to_rhos

  !> Convert entropy to energy
  subroutine rhos_to_e(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in) :: ixI^L, ixO^L
    double precision :: w(ixI^S,nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)

    if (mhd_energy) then
       w(ixO^S, e_) = w(ixO^S, rho_)**(mhd_gamma - 1.0d0) * w(ixO^S, e_) &
            / (mhd_gamma - 1.0d0) + mhd_kin_en(w, ixI^L, ixO^L) + &
            mhd_mag_en(w, ixI^L, ixO^L)
    else
       call mpistop("rhos_to_e can not be used without energy equation!")
    end if
  end subroutine rhos_to_e

  !> Calculate v vector
  subroutine mhd_get_v(w,x,ixI^L,ixO^L,v)
    use mod_global_parameters

    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S,nw), x(ixI^S,1:ndim)
    double precision, intent(out) :: v(ixI^S,ndir)

    integer :: idir

    do idir=1,ndir
      v(ixO^S,idir) = w(ixO^S, mom(idir)) * mhd_inv_rho(w, ixI^L, ixO^L)
    end do

  end subroutine mhd_get_v

  !> Calculate cmax_idim=csound+abs(v_idim) within ixO^L
  subroutine mhd_get_cmax(w,x,ixI^L,ixO^L,idim,cmax,cmin)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L, idim
    double precision, intent(in) :: w(ixI^S, nw), x(ixI^S,1:ndim)
    double precision, intent(inout)           :: cmax(ixI^S)
    double precision, intent(inout), optional :: cmin(ixI^S)

    double precision :: csound(ixI^S), cfast2(ixI^S), AvMinCs2(ixI^S), v(ixI^S), kmax

    call mhd_get_csound2(w,x,ixI^L,ixO^L,csound)
    ! store |B|^2 in v
    v(ixO^S)        = mhd_mag_en_all(w,ixI^L,ixO^L)
    cfast2(ixO^S)   = v(ixO^S) / w(ixO^S,rho_)+csound(ixO^S)
    AvMinCs2(ixO^S) = cfast2(ixO^S)**2-4.0d0*csound(ixO^S) &
         * mhd_mag_i_all(w,ixI^L,ixO^L,idim)**2 &
         / w(ixO^S,rho_)

    where(AvMinCs2(ixO^S)<zero)
       AvMinCs2(ixO^S)=zero
    end where

    AvMinCs2(ixO^S)=sqrt(AvMinCs2(ixO^S))

    if (.not. MHD_Hall) then
       csound(ixO^S) = sqrt(half*(cfast2(ixO^S)+AvMinCs2(ixO^S)))
    else
       ! take the Hall velocity into account:
       ! most simple estimate, high k limit:
       ! largest wavenumber supported by grid: Nyquist (in practise can reduce by some factor)
       kmax = dpi/min({dxlevel(^D)},bigdouble)*half
       csound(ixO^S) = max(sqrt(half*(cfast2(ixO^S)+AvMinCs2(ixO^S))), &
            mhd_etah * sqrt(v(ixO^S))/w(ixO^S,rho_)*kmax)
    end if

    v(ixO^S)=w(ixO^S,mom(idim))/w(ixO^S,rho_)

    if (present(cmin))then
       cmax(ixO^S)=max(v(ixO^S)+csound(ixO^S),zero)
       cmin(ixO^S)=min(v(ixO^S)-csound(ixO^S),zero)
    else
       cmax(ixO^S) = abs(v(ixO^S))+csound(ixO^S)
    end if

  end subroutine mhd_get_cmax

  !> Calculate thermal pressure=(gamma-1)*(e-0.5*m**2/rho-b**2/2) within ixO^L
  subroutine mhd_get_pthermal(w,x,ixI^L,ixO^L,pth)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S,nw)
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(out):: pth(ixI^S)

    if (mhd_energy) then
       pth(ixO^S)=(mhd_gamma-1.0d0)*(w(ixO^S,e_)&
          - mhd_kin_en(w,ixI^L,ixO^L)&
          - mhd_mag_en(w,ixI^L,ixO^L))
    else
       pth(ixO^S)=mhd_adiab*w(ixO^S,rho_)**mhd_gamma
    end if
  end subroutine mhd_get_pthermal

  !> Calculate the square of the thermal sound speed csound2 within ixO^L.
  !> csound2=gamma*p/rho
  subroutine mhd_get_csound2(w,x,ixI^L,ixO^L,csound2)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: w(ixI^S,nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out)   :: csound2(ixI^S)

    if (mhd_energy) then
       call mhd_get_pthermal(w,x,ixI^L,ixO^L,csound2)
       csound2(ixO^S)=mhd_gamma*csound2(ixO^S)/w(ixO^S,rho_)
    else
       csound2(ixO^S)=mhd_gamma*mhd_adiab*w(ixO^S,rho_)**(mhd_gamma-one)
    end if
  end subroutine mhd_get_csound2

  !> Calculate total pressure within ixO^L including magnetic pressure
  subroutine mhd_get_p_total(w,x,ixI^L,ixO^L,p)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: w(ixI^S,nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out)   :: p(ixI^S)

    integer, dimension(ixI^S)       :: patchierror
    integer, dimension(ndim)       :: lowpindex

    call mhd_get_pthermal(w,x,ixI^L,ixO^L,p)

    p(ixO^S) = p(ixO^S) + 0.5d0 * sum(w(ixO^S, mag(:))**2, dim=ndim+1)

  end subroutine mhd_get_p_total

  !> Calculate non-transport fluxes within ixO^L.
  subroutine mhd_get_flux(w,x,ixI^L,ixO^L,idim,f)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L, idim
    double precision, intent(in) :: w(ixI^S,nw)
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision,intent(out) :: f(ixI^S,nwflux)

    double precision             :: ptotal(ixI^S),tmp(ixI^S), v(ixI^S,ndir)
    double precision, allocatable:: vHall(:^D&,:)
    integer                      :: idirmin, iw, idir

    call mhd_get_v(w,x,ixI^L,ixO^L,v)

    if (mhd_Hall) then
      allocate(vHall(ixI^S,1:ndir))
      call mhd_getv_Hall(w,x,ixI^L,ixO^L,vHall)
    end if

    if(B0field) tmp(ixO^S)=sum(block%B0(ixO^S,:,idim)*w(ixO^S,mag(:)),dim=ndim+1)

    call mhd_get_p_total(w,x,ixI^L,ixO^L,ptotal)

    ! Get flux of density
    f(ixO^S,rho_)=v(ixO^S,idim)*w(ixO^S,rho_)

    ! Get flux of tracer
    do iw=1,mhd_n_tracer
      f(ixO^S,tracer(iw))=v(ixO^S,idim)*w(ixO^S,tracer(iw))
    end do

    ! Get flux of momentum
    ! f_i[m_k]=v_i*m_k-b_k*b_i [+ptotal if i==k]
    do idir=1,ndir
      if(idim==idir) then
        f(ixO^S,mom(idir))=ptotal(ixO^S)-w(ixO^S,mag(idim))*w(ixO^S,mag(idir))
        if(B0field) f(ixO^S,mom(idir))=f(ixO^S,mom(idir))+tmp(ixO^S)
      else
        f(ixO^S,mom(idir))= -w(ixO^S,mag(idir))*w(ixO^S,mag(idim))
      end if
      if (B0field) then
        f(ixO^S,mom(idir))=f(ixO^S,mom(idir))&
             -w(ixO^S,mag(idir))*block%B0(ixO^S,idim,idim)&
             -w(ixO^S,mag(idim))*block%B0(ixO^S,idir,idim)
      end if
      f(ixO^S,mom(idir))=f(ixO^S,mom(idir))+v(ixO^S,idim)*w(ixO^S,mom(idir))
    end do

    ! Get flux of energy
    ! f_i[e]=v_i*e+v_i*ptotal-b_i*(b_k*v_k)
    if(mhd_energy) then
      f(ixO^S,e_)=v(ixO^S,idim)*(w(ixO^S,e_)+ptotal(ixO^S))- &
            w(ixO^S,mag(idim))*sum(w(ixO^S,mag(:))*v(ixO^S,:),dim=ndim+1)

      if (B0field) then
        f(ixO^S,e_) = f(ixO^S,e_) &
             + v(ixO^S,idim) * tmp(ixO^S) &
             - sum(v(ixO^S,:)*w(ixO^S,mag(:)),dim=ndim+1) * block%B0(ixO^S,idim,idim)
      end if

      if (mhd_Hall) then
        ! f_i[e]= f_i[e] + vHall_i*(b_k*b_k) - b_i*(vHall_k*b_k)
        if (mhd_etah>zero) then
          f(ixO^S,e_) = f(ixO^S,e_) + vHall(ixO^S,idim) * &
                 sum(w(ixO^S, mag(:))**2,dim=ndim+1) &
               - w(ixO^S,mag(idim)) * sum(vHall(ixO^S,:)*w(ixO^S,mag(:)),dim=ndim+1)
          if (B0field) then
            f(ixO^S,e_) = f(ixO^S,e_) &
                 + vHall(ixO^S,idim) * tmp(ixO^S) &
                 - sum(vHall(ixO^S,:)*w(ixO^S,mag(:)),dim=ndim+1) * block%B0(ixO^S,idim,idim)
          end if
        end if
      end if

    end if

    ! compute flux of magnetic field
    ! f_i[b_k]=v_i*b_k-v_k*b_i
    do idir=1,ndir
      if (idim==idir) then
        ! f_i[b_i] should be exactly 0, so we do not use the transport flux
        if (mhd_glm) then
           f(ixO^S,mag(idir))=w(ixO^S,psi_)
        else
           f(ixO^S,mag(idir))=zero
        end if
      else
        f(ixO^S,mag(idir))=v(ixO^S,idim)*w(ixO^S,mag(idir))-w(ixO^S,mag(idim))*v(ixO^S,idir)

        if (B0field) then
          f(ixO^S,mag(idir))=f(ixO^S,mag(idir))&
                +v(ixO^S,idim)*block%B0(ixO^S,idir,idim)&
                -v(ixO^S,idir)*block%B0(ixO^S,idim,idim)
        end if

        if (mhd_Hall) then
          ! f_i[b_k] = f_i[b_k] + vHall_i*b_k - vHall_k*b_i
          if (mhd_etah>zero) then
            if (B0field) then
              f(ixO^S,mag(idir)) = f(ixO^S,mag(idir)) &
                   - vHall(ixO^S,idir)*(w(ixO^S,mag(idim))+block%B0(ixO^S,idim,idim)) &
                   + vHall(ixO^S,idim)*(w(ixO^S,mag(idir))+block%B0(ixO^S,idir,idim))
            else
              f(ixO^S,mag(idir)) = f(ixO^S,mag(idir)) &
                   - vHall(ixO^S,idir)*w(ixO^S,mag(idim)) &
                   + vHall(ixO^S,idim)*w(ixO^S,mag(idir))
            end if
          end if
        end if

      end if
    end do

    if (mhd_glm) then
      !f_i[psi]=Ch^2*b_{i} Eq. 24e and Eq. 38c Dedner et al 2002 JCP, 175, 645
      f(ixO^S,psi_)  = cmax_global**2*w(ixO^S,mag(idim))
    end if

  end subroutine mhd_get_flux

  !> w[iws]=w[iws]+qdt*S[iws,wCT] where S is the source based on wCT within ixO
  subroutine mhd_add_source(qdt,ixI^L,ixO^L,wCT,w,x,qsourcesplit)
    use mod_global_parameters
    use mod_radiative_cooling, only: radiative_cooling_add_source
    use mod_viscosity, only: viscosity_add_source
    use mod_gravity, only: gravity_add_source

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    logical, intent(in)             :: qsourcesplit

    if (.not. qsourcesplit) then
       ! Sources for resistivity in eqs. for e, B1, B2 and B3
       if (dabs(mhd_eta)>smalldouble)then
          if (.not.slab) call mpistop("no resistivity in non-slab geometry")
          if (compactres)then
             call add_source_res1(qdt,ixI^L,ixO^L,wCT,w,x)
          else
             call add_source_res2(qdt,ixI^L,ixO^L,wCT,w,x)
          end if
       end if

       if (mhd_eta_hyper>0.d0)then
          call add_source_hyperres(qdt,ixI^L,ixO^L,wCT,w,x)
       end if
    end if

    if(mhd_radiative_cooling) then
      call radiative_cooling_add_source(qdt,ixI^L,ixO^L,wCT,w,x,qsourcesplit)
    end if

    if(mhd_viscosity) then
      call viscosity_add_source(qdt,ixI^L,ixO^L,wCT,w,x,mhd_energy,qsourcesplit)
    end if

    if(mhd_gravity) then
      call gravity_add_source(qdt,ixI^L,ixO^L,wCT,w,x,mhd_energy,qsourcesplit)
    end if

    if(B0field .and. .not.qsourcesplit) then
      call add_source_B0split(qdt,ixI^L,ixO^L,wCT,w,x)
    end if

    {^NOONED
    if (qsourcesplit) then
       ! Sources related to div B
       select case (typedivbfix)
       case ('glm1')
          call add_source_glm1(qdt,ixI^L,ixO^L,wCT,w,x)
       case ('glm2')
          call add_source_glm2(qdt,ixI^L,ixO^L,wCT,w,x)
       case ('glm3')
          call add_source_glm3(qdt,ixI^L,ixO^L,wCT,w,x)
       case ('powel', 'powell')
          call add_source_powel(qdt,ixI^L,ixO^L,wCT,w,x)
       case ('janhunen')
          call add_source_janhunen(qdt,ixI^L,ixO^L,wCT,w,x)
       case ('linde')
          call add_source_linde(qdt,ixI^L,ixO^L,wCT,w,x)
       case ('lindejanhunen')
          call add_source_linde(qdt,ixI^L,ixO^L,wCT,w,x)
          call add_source_janhunen(qdt,ixI^L,ixO^L,wCT,w,x)
       case ('lindepowel')
          call add_source_linde(qdt,ixI^L,ixO^L,wCT,w,x)
          call add_source_powel(qdt,ixI^L,ixO^L,wCT,w,x)
       case ('none')
         ! Do nothing
       case default
         call mpistop('Unknown divB fix')
       end select
    end if
    }
  end subroutine mhd_add_source

  !> Source terms after split off time-independent magnetic field
  subroutine add_source_B0split(qdt,ixI^L,ixO^L,wCT,w,x)
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: qdt, wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: a(ixI^S,3), b(ixI^S,3), axb(ixI^S,3)
    integer :: idir

    a=0.d0
    b=0.d0
    ! for force-free field J0xB0 =0
    if(.not.B0field_forcefree) then
      ! store B0 magnetic field in b
      b(ixO^S,:)=block%B0(ixO^S,:,0)

      ! store J0 current in a
      do idir=7-2*ndir,3
        a(ixO^S,idir)=block%J0(ixO^S,idir)
      end do
      call cross_product(ixI^L,ixO^L,a,b,axb)
      axb(ixO^S,:)=axb(ixO^S,:)*qdt
      ! add J0xB0 source term in momentum equations
      do idir=1,ndir
        w(ixO^S,mom(idir))=w(ixO^S,mom(idir))+axb(ixO^S,idir)
      end do
    end if

    if(mhd_energy) then
      a=0.d0
      ! for free-free field -(vxB0) dot J0 =0
      b(ixO^S,:)=wCT(ixO^S,mag(:))
      ! store full magnetic field B0+B1 in b
      if(.not.B0field_forcefree) b(ixO^S,:)=b(ixO^S,:)+block%B0(ixO^S,:,0)
      ! store velocity in a
      do idir=1,ndir
        a(ixO^S,idir)=wCT(ixO^S,mom(idir))/wCT(ixO^S,rho_)
      end do
      call cross_product(ixI^L,ixO^L,a,b,axb)
      axb(ixO^S,:)=axb(ixO^S,:)*qdt
      ! add -(vxB) dot J0 source term in energy equation
      do idir=7-2*ndir,3
        w(ixO^S,e_)=w(ixO^S,e_)-axb(ixO^S,idir)*block%J0(ixO^S,idir)
      end do
    end if

  end subroutine add_source_B0split

  !> Add resistive source to w within ixO Uses 3 point stencil (1 neighbour) in
  !> each direction, non-conservative. If the fourthorder precompiler flag is
  !> set, uses fourth order central difference for the laplacian. Then the
  !> stencil is 5 (2 neighbours).
  subroutine add_source_res1(qdt,ixI^L,ixO^L,wCT,w,x)
    use mod_global_parameters
    use mod_usr_methods

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt
    double precision, intent(in) :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    integer :: ixA^L,idir,jdir,kdir,idirmin,idim,jxO^L,hxO^L,ix
    integer :: lxO^L, kxO^L

    double precision :: tmp(ixI^S),tmp2(ixI^S)

    ! For ndir=2 only 3rd component of J can exist, ndir=1 is impossible for MHD
    double precision :: current(ixI^S,7-2*ndir:3),eta(ixI^S)
    double precision :: gradeta(ixI^S,1:ndim)


    ! Calculating resistive sources involve one extra layer
    if (mhd_4th_order) then
      ixA^L=ixO^L^LADD2;
    else
      ixA^L=ixO^L^LADD1;
    end if

    if (ixImin^D>ixAmin^D.or.ixImax^D<ixAmax^D|.or.) &
         call mpistop("Error in add_source_res1: Non-conforming input limits")

    ! Calculate current density and idirmin
    call get_current(wCT,ixI^L,ixO^L,idirmin,current)

    if (mhd_eta>zero)then
       eta(ixA^S)=mhd_eta
       gradeta(ixO^S,1:ndim)=zero
    else
       call usr_special_resistivity(wCT,ixI^L,ixA^L,idirmin,x,current,eta)
       ! assumes that eta is not function of current?
       do idim=1,ndim
          call gradient(eta,ixI^L,ixO^L,idim,tmp)
          gradeta(ixO^S,idim)=tmp(ixO^S)
       end do
    end if

    do idir=1,ndir
       ! Put B_idir into tmp2 and eta*Laplace B_idir into tmp
       if (mhd_4th_order) then
         tmp(ixO^S)=zero
         tmp2(ixI^S)=wCT(ixI^S,mag(idir))
         do idim=1,ndim
            lxO^L=ixO^L+2*kr(idim,^D);
            jxO^L=ixO^L+kr(idim,^D);
            hxO^L=ixO^L-kr(idim,^D);
            kxO^L=ixO^L-2*kr(idim,^D);
            tmp(ixO^S)=tmp(ixO^S)+&
                 (-tmp2(lxO^S)+16.0d0*tmp2(jxO^S)-30.0d0*tmp2(ixO^S)+16.0d0*tmp2(hxO^S)-tmp2(kxO^S)) &
                 /(12.0d0 * dxlevel(idim)**2)
         end do
       else
         tmp(ixO^S)=zero
         tmp2(ixI^S)=wCT(ixI^S,mag(idir))
         do idim=1,ndim
            jxO^L=ixO^L+kr(idim,^D);
            hxO^L=ixO^L-kr(idim,^D);
            tmp(ixO^S)=tmp(ixO^S)+&
                 (tmp2(jxO^S)-2.0d0*tmp2(ixO^S)+tmp2(hxO^S))/dxlevel(idim)**2
         end do
       end if

       ! Multiply by eta
       tmp(ixO^S)=tmp(ixO^S)*eta(ixO^S)

       ! Subtract grad(eta) x J = eps_ijk d_j eta J_k if eta is non-constant
       if (mhd_eta<zero)then
          do jdir=1,ndim; do kdir=idirmin,3
             if (lvc(idir,jdir,kdir)/=0)then
                if (lvc(idir,jdir,kdir)==1)then
                   tmp(ixO^S)=tmp(ixO^S)-gradeta(ixO^S,jdir)*current(ixO^S,kdir)
                else
                   tmp(ixO^S)=tmp(ixO^S)+gradeta(ixO^S,jdir)*current(ixO^S,kdir)
                end if
             end if
          end do; end do
       end if

       ! Add sources related to eta*laplB-grad(eta) x J to B and e
       w(ixO^S,mag(idir))=w(ixO^S,mag(idir))+qdt*tmp(ixO^S)
       if (mhd_energy) then
          w(ixO^S,e_)=w(ixO^S,e_)+qdt*tmp(ixO^S)*wCT(ixO^S,mag(idir))
       end if

    end do ! idir

    if (mhd_energy) then
       ! de/dt+=eta*J**2
       tmp(ixO^S)=zero
       do idir=idirmin,3
          tmp(ixO^S)=tmp(ixO^S)+current(ixO^S,idir)**2
       end do
       w(ixO^S,e_)=w(ixO^S,e_)+qdt*eta(ixO^S)*tmp(ixO^S)

       call handle_small_values(.false.,w,x,ixI^L,ixO^L,'add_source_res1')
    end if
  end subroutine add_source_res1

  !> Add resistive source to w within ixO
  !> Uses 5 point stencil (2 neighbours) in each direction, conservative
  subroutine add_source_res2(qdt,ixI^L,ixO^L,wCT,w,x)
    use mod_global_parameters
    use mod_usr_methods

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    integer :: ixA^L,idir,jdir,kdir,idirmin,iw,idim,idirmin1

    double precision :: tmp(ixI^S),tmp2(ixI^S)

    ! For ndir=2 only 3rd component of J can exist, ndir=1 is impossible for MHD
    double precision :: current(ixI^S,7-2*ndir:3),eta(ixI^S),curlj(ixI^S,1:3)
    double precision :: tmpvec(ixI^S,1:3),tmpvec2(ixI^S,1:ndir)

    ixA^L=ixO^L^LADD2;

    if (ixImin^D>ixAmin^D.or.ixImax^D<ixAmax^D|.or.) &
         call mpistop("Error in add_source_res2: Non-conforming input limits")

    ixA^L=ixO^L^LADD1;
    ! Calculate current density within ixL: J=curl B, thus J_i=eps_ijk*d_j B_k
    ! Determine exact value of idirmin while doing the loop.
    call get_current(wCT,ixI^L,ixA^L,idirmin,current)

    if (mhd_eta>zero)then
       eta(ixA^S)=mhd_eta
    else
       call usr_special_resistivity(wCT,ixI^L,ixA^L,idirmin,x,current,eta)
    end if

    ! dB/dt= -curl(J*eta), thus B_i=B_i-eps_ijk d_j Jeta_k
    tmpvec(ixA^S,1:ndir)=zero
    do jdir=idirmin,3
       tmpvec(ixA^S,jdir)=current(ixA^S,jdir)*eta(ixA^S)*qdt
    end do
    call curlvector(tmpvec,ixI^L,ixO^L,curlj,idirmin1,1,3)
    do idir=1,ndir
      w(ixO^S,mag(idir)) = w(ixO^S,mag(idir))-curlj(ixO^S,idir)
    end do

    if (mhd_energy) then
       ! de/dt= +div(B x Jeta)
       tmpvec2(ixA^S,1:ndir)=zero
       do idir=1,ndir; do jdir=1,ndir; do kdir=idirmin,3
          if (lvc(idir,jdir,kdir)/=0)then
             tmp(ixA^S)=wCT(ixA^S,mag(jdir))*current(ixA^S,kdir)*eta(ixA^S)*qdt
             if (lvc(idir,jdir,kdir)==1)then
                tmpvec2(ixA^S,idir)=tmpvec2(ixA^S,idir)+tmp(ixA^S)
             else
                tmpvec2(ixA^S,idir)=tmpvec2(ixA^S,idir)-tmp(ixA^S)
             end if
          end if
       end do; end do; end do

       call divvector(tmpvec2,ixI^L,ixO^L,tmp)

       w(ixO^S,e_)=w(ixO^S,e_)+tmp(ixO^S)

       call handle_small_values(.false.,w,x,ixI^L,ixO^L,'add_source_res2')
    end if
  end subroutine add_source_res2

  !> Add Hyper-resistive source to w within ixO
  !> Uses 9 point stencil (4 neighbours) in each direction.
  subroutine add_source_hyperres(qdt,ixI^L,ixO^L,wCT,w,x)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    !.. local ..
    double precision                :: current(ixI^S,7-2*ndir:3)
    double precision                :: tmpvec(ixI^S,1:3),tmpvec2(ixI^S,1:3),tmp(ixI^S),ehyper(ixI^S,1:3)
    integer                         :: ixA^L,idir,jdir,kdir,idirmin,idirmin1

    ixA^L=ixO^L^LADD3;
    if (ixImin^D>ixAmin^D.or.ixImax^D<ixAmax^D|.or.) &
         call mpistop("Error in add_source_hyperres: Non-conforming input limits")

    call get_current(wCT,ixI^L,ixA^L,idirmin,current)
    tmpvec(ixA^S,1:ndir)=zero
    do jdir=idirmin,3
       tmpvec(ixA^S,jdir)=current(ixA^S,jdir)
    end do

    ixA^L=ixO^L^LADD2;
    call curlvector(tmpvec,ixI^L,ixA^L,tmpvec2,idirmin1,1,3)

    ixA^L=ixO^L^LADD1;
    tmpvec(ixA^S,1:ndir)=zero
    call curlvector(tmpvec2,ixI^L,ixA^L,tmpvec,idirmin1,1,3)
    ehyper(ixA^S,1:ndir) = - tmpvec(ixA^S,1:ndir)*mhd_eta_hyper

    ixA^L=ixO^L;
    tmpvec2(ixA^S,1:ndir)=zero
    call curlvector(ehyper,ixI^L,ixA^L,tmpvec2,idirmin1,1,3)

    do idir=1,ndir
      w(ixO^S,mag(idir)) = w(ixO^S,mag(idir))-tmpvec2(ixO^S,idir)*qdt
    end do

    if (mhd_energy) then
       ! de/dt= +div(B x Ehyper)
       ixA^L=ixO^L^LADD1;
       tmpvec2(ixA^S,1:ndir)=zero
       do idir=1,ndir; do jdir=1,ndir; do kdir=idirmin,3
          tmpvec2(ixA^S,idir) = tmpvec(ixA^S,idir)&
               + lvc(idir,jdir,kdir)*wCT(ixA^S,mag(jdir))*ehyper(ixA^S,kdir)
       end do; end do; end do
       tmp(ixO^S)=zero
       call divvector(tmpvec2,ixI^L,ixO^L,tmp)
       w(ixO^S,e_)=w(ixO^S,e_)+tmp(ixO^S)*qdt

       call handle_small_values(.false.,w,x,ixI^L,ixO^L,'add_source_hyperres')
    end if
  end subroutine add_source_hyperres

  subroutine add_source_glm1(qdt,ixI^L,ixO^L,wCT,w,x)
    ! Add divB related sources to w within ixO
    ! corresponding to Dedner JCP 2002, 175, 645 _equation 24_
    ! giving the EGLM-MHD scheme
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: qdt, wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision:: divb(ixI^S)
    integer          :: idim,idir
    double precision :: gradPsi(ixI^S)

    ! We calculate now div B
    call get_divb(wCT,ixI^L,ixO^L,divb)

    ! Psi = Psi - qdt Ch^2/Cp^2 Psi
    if (mhd_glm_Cr < zero) then
      w(ixO^S,psi_) = abs(mhd_glm_Cr)*wCT(ixO^S,psi_)
    else
      ! implicit update of psi variable
      w(ixO^S,psi_) = dexp(-qdt*(cmax_global/mhd_glm_Cr))*wCT(ixO^S,psi_)
    end if

    ! gradient of Psi
    do idim=1,ndim
       select case(typegrad)
       case("central")
          call gradient(wCT(ixI^S,psi_),ixI^L,ixO^L,idim,gradPsi)
       case("limited")
          call gradientS(wCT(ixI^S,psi_),ixI^L,ixO^L,idim,gradPsi)
       end select
       if (mhd_energy) then
       ! e  = e  -qdt (b . grad(Psi))
         w(ixO^S,e_) = w(ixO^S,e_)-qdt*wCT(ixO^S,mag(idim))*gradPsi(ixO^S)
       end if
    end do

    ! m = m - qdt b div b
    do idir=1,ndir
      w(ixO^S,mom(idir))=w(ixO^S,mom(idir))-qdt*mhd_mag_i_all(w,ixI^L,ixO^L,idir)*divb(ixO^S)
    end do

    call handle_small_values(.false.,w,x,ixI^L,ixO^L,'add_source_glm1')

  end subroutine add_source_glm1

  subroutine add_source_glm2(qdt,ixI^L,ixO^L,wCT,w,x)
    ! Add divB related sources to w within ixO
    ! corresponding to Dedner JCP 2002, 175, 645 _equation 38_
    ! giving the non conservative EGLM-MHD scheme.
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: qdt,   wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision:: divb(ixI^S),v(ixI^S,1:ndir)
    integer          :: idim,idir
    double precision :: gradPsi(ixI^S)

    ! calculate now div B
    call get_divb(wCT,ixI^L,ixO^L,divb)

    ! calculate velocity
    call mhd_get_v(wCT,x,ixI^L,ixO^L,v)

    ! Psi = Psi - qdt Ch^2/Cp^2 Psi
    if (mhd_glm_Cr < zero) then
      w(ixO^S,psi_) = abs(mhd_glm_Cr)*wCT(ixO^S,psi_)
    else
      ! implicit update of psi variable
      w(ixO^S,psi_) = dexp(-qdt*(cmax_global/mhd_glm_Cr))*wCT(ixO^S,psi_)
    end if

    ! gradient of Psi
    do idim=1,ndim
       select case(typegrad)
       case("central")
         call gradient(wCT(ixI^S,psi_),ixI^L,ixO^L,idim,gradPsi)
       case("limited")
         call gradientS(wCT(ixI^S,psi_),ixI^L,ixO^L,idim,gradPsi)
       end select

      ! Psi=Psi - qdt (v . grad(Psi))
      w(ixO^S,psi_) = w(ixO^S,psi_)-qdt*v(ixO^S,idim)*gradPsi(ixO^S)

      if (mhd_energy) then
      ! e  = e  - qdt (b . grad(Psi))
        w(ixO^S,e_) = w(ixO^S,e_)&
             -qdt*wCT(ixO^S,mag(idim))*gradPsi(ixO^S)
      end if
    end do

    if (mhd_energy) then
    ! e = e - qdt (v . b) * div b
       w(ixO^S,e_)=w(ixO^S,e_) - qdt * divb(ixO^S) * &
            sum(v(ixO^S,:)*wCT(ixO^S,mag(:)),dim=ndim+1)
    end if
    ! b = b - qdt v * div b
    do idir=1,ndir
      w(ixO^S,mag(idir))=w(ixO^S,mag(idir))-qdt*v(ixO^S,idir)*divb(ixO^S)
    end do

    ! m = m - qdt b div b
    do idir=1,ndir
      w(ixO^S,mom(idir))=w(ixO^S,mom(idir))-qdt*mhd_mag_i_all(w,ixI^L,ixO^L,idir)*divb(ixO^S)
    end do

    call handle_small_values(.false.,w,x,ixI^L,ixO^L,'add_source_glm2')

  end subroutine add_source_glm2

  subroutine add_source_glm3(qdt,ixI^L,ixO^L,wCT,w,x)
    ! Add divB related sources to w within ixO
    ! corresponding to Dedner JCP 2002, 175, 645 _equation (1a), (1b), (4), (1d), 19
    ! conservative hyperbolic mixed GLM-MHD with no additional source terms.
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt,   wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    ! Psi = Psi - qdt Ch^2/Cp^2 Psi
    if (mhd_glm_Cr < zero) then
      w(ixO^S,psi_) = abs(mhd_glm_Cr)*w(ixO^S,psi_)
    else
      ! implicit update of psi variable
      w(ixO^S,psi_) = dexp(-qdt*(cmax_global/mhd_glm_Cr))*w(ixO^S,psi_)
    end if

  end subroutine add_source_glm3

  !> Add divB related sources to w within ixO corresponding to Powel
  subroutine add_source_powel(qdt,ixI^L,ixO^L,wCT,w,x)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt,   wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision                :: divb(ixI^S),v(ixI^S,1:ndir)
    integer                         :: idir

    ! We calculate now div B
    call get_divb(wCT,ixI^L,ixO^L,divb)

    ! calculate velocity
    call mhd_get_v(wCT,x,ixI^L,ixO^L,v)

    if (mhd_energy) then
      ! e = e - qdt (v . b) * div b
      w(ixO^S,e_)=w(ixO^S,e_)-&
           qdt*sum(v(ixO^S,:)*wCT(ixO^S,mag(:)),dim=ndim+1)*divb(ixO^S)
    end if

    ! b = b - qdt v * div b
    do idir=1,ndir
      w(ixO^S,mag(idir))=w(ixO^S,mag(idir))-qdt*v(ixO^S,idir)*divb(ixO^S)
    end do

    ! m = m - qdt b div b
    do idir=1,ndir
      w(ixO^S,mom(idir))=w(ixO^S,mom(idir))-qdt*mhd_mag_i_all(w,ixI^L,ixO^L,idir)*divb(ixO^S)
    end do

    call handle_small_values(.false.,w,x,ixI^L,ixO^L,'add_source_powel')

  end subroutine add_source_powel

  subroutine add_source_janhunen(qdt,ixI^L,ixO^L,wCT,w,x)
    ! Add divB related sources to w within ixO
    ! corresponding to Janhunen, just the term in the induction equation.
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt,   wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision                :: divb(ixI^S)
    integer                         :: idir

    ! We calculate now div B
    call get_divb(wCT,ixI^L,ixO^L,divb)

    ! b = b - qdt v * div b
    do idir=1,ndir
      w(ixO^S,mag(idir))=w(ixO^S,mag(idir))-qdt*wCT(ixO^S,mom(idir))/wCT(ixO^S,rho_)*divb(ixO^S)
    end do

    call handle_small_values(.false.,w,x,ixI^L,ixO^L,'add_source_janhunen')

  end subroutine add_source_janhunen

  subroutine add_source_linde(qdt,ixI^L,ixO^L,wCT,w,x)
    ! Add Linde's divB related sources to wnew within ixO
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    integer :: idim, idir, ixp^L, i^D, iside
    double precision :: divb(ixI^S),graddivb(ixI^S)


    ! Calculate div B
    ixp^L=ixO^L^LADD1;
    call get_divb(wCT,ixI^L,ixp^L,divb)
    ! for AMR stability, retreat one cell layer from the boarders of level jump
    ixp^L=ixO^L;
    do idim=1,ndim
      select case(idim)
       {case(^D)
          do iside=1,2
            i^DD=kr(^DD,^D)*(2*iside-3);
            if (leveljump(i^DD)) then
              if (iside==1) then
                ixpmin^D=ixOmin^D-i^D
              else
                ixpmax^D=ixOmax^D-i^D
              end if
            end if
          end do
       \}
      end select
    end do

    ! Add Linde's diffusive terms
    do idim=1,ndim
       ! Calculate grad_idim(divb)
       select case(typegrad)
       case("central")
         call gradient(divb,ixI^L,ixp^L,idim,graddivb)
       case("limited")
         call gradientS(divb,ixI^L,ixp^L,idim,graddivb)
       end select

       ! Multiply by Linde's eta*dt = divbdiff*(c_max*dx)*dt = divbdiff*dx**2
       if (slab) then
          graddivb(ixp^S)=graddivb(ixp^S)*divbdiff/(^D&1.0d0/dxlevel(^D)**2+)
       else
          graddivb(ixp^S)=graddivb(ixp^S)*divbdiff &
                          /(^D&1.0d0/block%dx(ixp^S,^D)**2+)
       end if

       w(ixp^S,mag(idim))=w(ixp^S,mag(idim))+graddivb(ixp^S)

       if (mhd_energy .and. typedivbdiff=='all') then
         ! e += B_idim*eta*grad_idim(divb)
         w(ixp^S,e_)=w(ixp^S,e_)+wCT(ixp^S,mag(idim))*graddivb(ixp^S)
       end if
    end do

    call handle_small_values(.false.,w,x,ixI^L,ixO^L,'add_source_linde')

  end subroutine add_source_linde

  subroutine get_divb(w,ixI^L,ixO^L,divb)

    ! Calculate div B within ixO

    use mod_global_parameters

    integer, intent(in)                :: ixI^L, ixO^L
    double precision, intent(in)       :: w(ixI^S,1:nw)
    double precision                   :: divb(ixI^S)

    double precision                   :: bvec(ixI^S,1:ndir)

    bvec(ixI^S,:)=w(ixI^S,mag(:))

    select case(typediv)
    case("central")
      call divvector(bvec,ixI^L,ixO^L,divb)
    case("limited")
      call divvectorS(bvec,ixI^L,ixO^L,divb)
    end select
  end subroutine get_divb

  !> Calculate idirmin and the idirmin:3 components of the common current array
  !> make sure that dxlevel(^D) is set correctly.
  subroutine get_current(w,ixI^L,ixO^L,idirmin,current)
    use mod_global_parameters

    integer :: idirmin0
    integer :: ixO^L, idirmin, ixI^L
    double precision :: w(ixI^S,1:nw)
    integer :: idir

    ! For ndir=2 only 3rd component of J can exist, ndir=1 is impossible for MHD
    double precision :: current(ixI^S,7-2*ndir:3),bvec(ixI^S,1:ndir)

    idirmin0 = 7-2*ndir

    if (B0field) then
       do idir = 1, ndir
          bvec(ixI^S,idir)=w(ixI^S,mag(idir))+block%B0(ixI^S,idir,0)
       end do
    else
       do idir = 1, ndir
          bvec(ixI^S,idir)=w(ixI^S,mag(idir))
       end do
    end if

    call curlvector(bvec,ixI^L,ixO^L,current,idirmin,idirmin0,ndir)

  end subroutine get_current

  !> If resistivity is not zero, check diffusion time limit for dt
  subroutine mhd_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
    use mod_global_parameters
    use mod_usr_methods
    use mod_radiative_cooling, only: cooling_get_dt
    use mod_viscosity, only: viscosity_get_dt 
    use mod_gravity, only: gravity_get_dt

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: dtnew
    double precision, intent(in)    :: dx^D
    double precision, intent(in)    :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)

    integer                       :: idirmin,idim
    double precision              :: dxarr(ndim)
    double precision              :: current(ixI^S,7-2*ndir:3),eta(ixI^S)

    dtnew = bigdouble

    ^D&dxarr(^D)=dx^D;
    if (mhd_eta>zero)then
       dtnew=dtdiffpar*minval(dxarr(1:ndim))**2/mhd_eta
    else if (mhd_eta<zero)then
       call get_current(w,ixI^L,ixO^L,idirmin,current)
       call usr_special_resistivity(w,ixI^L,ixO^L,idirmin,x,current,eta)
       dtnew=bigdouble
       do idim=1,ndim
          dtnew=min(dtnew,&
               dtdiffpar/(smalldouble+maxval(eta(ixO^S)/dxarr(idim)**2)))
       end do
    end if

    if (mhd_eta_hyper>zero)then
       dtnew=min(dtdiffpar*minval(dxarr(1:ndim))**4/mhd_eta_hyper,dtnew)
    end if

    if(mhd_radiative_cooling) then
      call cooling_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
    end if

    if(mhd_viscosity) then
      call viscosity_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
    end if

    if(mhd_gravity) then
      call gravity_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
    end if

  end subroutine mhd_get_dt

  ! Add geometrical source terms to w
  subroutine mhd_add_source_geom(qdt,ixI^L,ixO^L,wCT,w,x)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: wCT(ixI^S,1:nw), w(ixI^S,1:nw)

    integer          :: iw,idir, h1x^L{^NOONED, h2x^L}
    double precision :: tmp(ixI^S),tmp1(ixI^S),tmp2(ixI^S)
    logical          :: angmomfix=.false.

    integer :: mr_,mphi_ ! Polar var. names
    integer :: br_,bphi_

    mr_=mom(1); mphi_=mom(1)-1+phi_  ! Polar var. names
    br_=mag(1); bphi_=mag(1)-1+phi_

    select case (typeaxial)
    case ('cylindrical')
       call mhd_get_p_total(wCT,x,ixI^L,ixO^L,tmp)
       if(phi_>0) then
         w(ixO^S,mr_)=w(ixO^S,mr_)+qdt/x(ixO^S,1)*(tmp(ixO^S)-&
                   wCT(ixO^S,bphi_)**2+wCT(ixO^S,mphi_)**2/wCT(ixO^S,rho_))
         w(ixO^S,mphi_)=w(ixO^S,mphi_)+qdt/x(ixO^S,1)*(&
                  -wCT(ixO^S,mphi_)*wCT(ixO^S,mr_)/wCT(ixO^S,rho_) &
                  +wCT(ixO^S,bphi_)*wCT(ixO^S,br_))
         w(ixO^S,bphi_)=w(ixO^S,bphi_)+qdt/x(ixO^S,1)*&
                  (wCT(ixO^S,bphi_)*wCT(ixO^S,mr_) &
                  -wCT(ixO^S,br_)*wCT(ixO^S,mphi_)) &
                  /wCT(ixO^S,rho_)
       else
         w(ixO^S,mr_)=w(ixO^S,mr_)+qdt/x(ixO^S,1)*tmp(ixO^S)
       end if
       if(mhd_glm) w(ixO^S,br_)=w(ixO^S,br_)+qdt*wCT(ixO^S,psi_)/x(ixO^S,1)
    case ('spherical')
       h1x^L=ixO^L-kr(1,^D); {^NOONED h2x^L=ixO^L-kr(2,^D);}
       call mhd_get_p_total(wCT,x,ixI^L,ixO^L,tmp1)
       tmp(ixO^S)=tmp1(ixO^S)
       if(B0field) then
         tmp2(ixO^S)=sum(block%B0(ixO^S,:,0)*wCT(ixO^S,mag(:)),dim=ndim+1)
         tmp(ixO^S)=tmp(ixO^S)+tmp2(ixO^S)
       end if
       ! m1
       tmp(ixO^S)=tmp(ixO^S)*x(ixO^S,1) &
                  *(block%surfaceC1(ixO^S)-block%surfaceC1(h1x^S))/block%dvolume(ixO^S)
       if(ndir>1) then
         do idir=2,ndir
           tmp(ixO^S)=tmp(ixO^S)+wCT(ixO^S,mom(idir))**2/wCT(ixO^S,rho_)-wCT(ixO^S,mag(idir))**2
           if(B0field) tmp(ixO^S)=tmp(ixO^S)-2.0d0*block%B0(ixO^S,idir,0)*wCT(ixO^S,mag(idir))
         end do
       end if
       w(ixO^S,mom(1))=w(ixO^S,mom(1))+qdt*tmp(ixO^S)/x(ixO^S,1)
       ! b1
       if(mhd_glm) then
         w(ixO^S,mag(1))=w(ixO^S,mag(1))+qdt/x(ixO^S,1)*2.0d0*wCT(ixO^S,psi_)
       end if

       {^NOONED
       ! m2
       tmp(ixO^S)=tmp1(ixO^S)
       if(B0field) then
         tmp(ixO^S)=tmp(ixO^S)+tmp2(ixO^S)
       end if
       ! This will make hydrostatic p=const an exact solution
       w(ixO^S,mom(2))=w(ixO^S,mom(2))+qdt*tmp(ixO^S) &
            *(block%surfaceC2(ixO^S)-block%surfaceC2(h2x^S)) &
            /block%dvolume(ixO^S)
       tmp(ixO^S)=-(wCT(ixO^S,mom(1))*wCT(ixO^S,mom(2))/wCT(ixO^S,rho_) &
            -wCT(ixO^S,mag(1))*wCT(ixO^S,mag(2)))
       if (B0field) then
          tmp(ixO^S)=tmp(ixO^S)+block%B0(ixO^S,1,0)*wCT(ixO^S,mag(2)) &
               +wCT(ixO^S,mag(1))*block%B0(ixO^S,2,0)
       end if
       if(ndir==3) then
         tmp(ixO^S)=tmp(ixO^S)+(wCT(ixO^S,mom(3))**2/wCT(ixO^S,rho_) &
              -wCT(ixO^S,mag(3))**2)*dcos(x(ixO^S,2))/dsin(x(ixO^S,2))
         if (B0field) then
            tmp(ixO^S)=tmp(ixO^S)-2.0d0*block%B0(ixO^S,3,0)*wCT(ixO^S,mag(3))&
                 *dcos(x(ixO^S,2))/dsin(x(ixO^S,2))
         end if
       end if
       w(ixO^S,mom(2))=w(ixO^S,mom(2))+qdt*tmp(ixO^S)/x(ixO^S,1)
       ! b2
       tmp(ixO^S)=(wCT(ixO^S,mom(1))*wCT(ixO^S,mag(2)) &
            -wCT(ixO^S,mom(2))*wCT(ixO^S,mag(1)))/wCT(ixO^S,rho_)
       if(B0field) then
         tmp(ixO^S)=tmp(ixO^S)+(wCT(ixO^S,mom(1))*block%B0(ixO^S,2,0) &
              -wCT(ixO^S,mom(2))*block%B0(ixO^S,1,0))/wCT(ixO^S,rho_)
       end if
       if(mhd_glm) then
         tmp(ixO^S)=tmp(ixO^S) &
              + dcos(x(ixO^S,2))/dsin(x(ixO^S,2))*wCT(ixO^S,psi_)
       end if
       w(ixO^S,mag(2))=w(ixO^S,mag(2))+qdt*tmp(ixO^S)/x(ixO^S,1)
       }

       if(ndir==3) then 
         ! m3
         if(.not.angmomfix) then
           tmp(ixO^S)=-(wCT(ixO^S,mom(3))*wCT(ixO^S,mom(1))/wCT(ixO^S,rho_) &
                -wCT(ixO^S,mag(3))*wCT(ixO^S,mag(1))) {^NOONED &
                -(wCT(ixO^S,mom(2))*wCT(ixO^S,mom(3))/wCT(ixO^S,rho_) &
                -wCT(ixO^S,mag(2))*wCT(ixO^S,mag(3))) &
                *dcos(x(ixO^S,2))/dsin(x(ixO^S,2)) }
           if (B0field) then
              tmp(ixO^S)=tmp(ixO^S)+block%B0(ixO^S,1,0)*wCT(ixO^S,mag(3)) &
                   +wCT(ixO^S,mag(1))*block%B0(ixO^S,3,0) {^NOONED &
                   +(block%B0(ixO^S,2,0)*wCT(ixO^S,mag(3)) &
                   +wCT(ixO^S,mag(2))*block%B0(ixO^S,3,0)) &
                   *dcos(x(ixO^S,2))/dsin(x(ixO^S,2)) }
           end if
           w(ixO^S,mom(3))=w(ixO^S,mom(3))+qdt*tmp(ixO^S)/x(ixO^S,1)
         end if
         ! b3
         tmp(ixO^S)=(wCT(ixO^S,mom(1))*wCT(ixO^S,mag(3)) &
              -wCT(ixO^S,mom(3))*wCT(ixO^S,mag(1)))/wCT(ixO^S,rho_) {^NOONED &
              -(wCT(ixO^S,mom(3))*wCT(ixO^S,mag(2)) &
              -wCT(ixO^S,mom(2))*wCT(ixO^S,mag(3)))*dcos(x(ixO^S,2)) &
              /(wCT(ixO^S,rho_)*dsin(x(ixO^S,2))) }
         if (B0field) then
            tmp(ixO^S)=tmp(ixO^S)+(wCT(ixO^S,mom(1))*block%B0(ixO^S,3,0) &
                 -wCT(ixO^S,mom(3))*block%B0(ixO^S,1,0))/wCT(ixO^S,rho_){^NOONED &
                 -(wCT(ixO^S,mom(3))*block%B0(ixO^S,2,0) &
                 -wCT(ixO^S,mom(2))*block%B0(ixO^S,3,0))*dcos(x(ixO^S,2)) &
                 /(wCT(ixO^S,rho_)*dsin(x(ixO^S,2))) }
         end if
         w(ixO^S,mag(3))=w(ixO^S,mag(3))+qdt*tmp(ixO^S)/x(ixO^S,1)
       end if
    end select
  end subroutine mhd_add_source_geom

  !> Compute 2 times total magnetic energy
  function mhd_mag_en_all(w, ixI^L, ixO^L) result(mge)
    use mod_global_parameters
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S, nw)
    double precision              :: mge(ixO^S)

    if (B0field) then
      mge = sum((w(ixO^S, mag(:))+block%B0(ixO^S,:,block%iw0))**2, dim=ndim+1)
    else
      mge = sum(w(ixO^S, mag(:))**2, dim=ndim+1)
    end if
  end function mhd_mag_en_all

  !> Compute full magnetic field by direction
  function mhd_mag_i_all(w, ixI^L, ixO^L,idir) result(mgf)
    use mod_global_parameters
    integer, intent(in)           :: ixI^L, ixO^L, idir
    double precision, intent(in)  :: w(ixI^S, nw)
    double precision              :: mgf(ixO^S)

    if (B0field) then
      mgf = w(ixO^S, mag(idir))+block%B0(ixO^S,idir,block%iw0)
    else
      mgf = w(ixO^S, mag(idir))
    end if
  end function mhd_mag_i_all

  !> Compute evolving magnetic energy
  function mhd_mag_en(w, ixI^L, ixO^L) result(mge)
    use mod_global_parameters, only: nw, ndim
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S, nw)
    double precision              :: mge(ixO^S)

    mge = 0.5d0 * sum(w(ixO^S, mag(:))**2, dim=ndim+1)
  end function mhd_mag_en

  !> compute kinetic energy
  function mhd_kin_en(w, ixI^L, ixO^L) result(ke)
    use mod_global_parameters, only: nw, ndim
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S, nw)
    double precision              :: ke(ixO^S)

    ke = 0.5d0 * sum(w(ixO^S, mom(:))**2, dim=ndim+1) * &
         mhd_inv_rho(w, ixI^L, ixO^L)
  end function mhd_kin_en

  function mhd_inv_rho(w, ixI^L, ixO^L) result(inv_rho)
    use mod_global_parameters, only: nw, ndim
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S, nw)
    double precision              :: inv_rho(ixO^S)

    ! Can make this more robust
    inv_rho = 1.0d0 / w(ixO^S, rho_)
  end function mhd_inv_rho

  subroutine mhd_getv_Hall(w,x,ixI^L,ixO^L,vHall)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: w(ixI^S,nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: vHall(ixI^S,1:3)

    integer          :: idir, idirmin
    double precision :: current(ixI^S,7-2*ndir:3)

    ! Calculate current density and idirmin
    call get_current(w,ixI^L,ixO^L,idirmin,current)
    vHall(ixI^S,1:3) = zero
    vHall(ixO^S,idirmin:3) = - mhd_etah*current(ixO^S,idirmin:3)
    do idir = idirmin, 3
       vHall(ixO^S,idir) = vHall(ixO^S,idir)/w(ixO^S,rho_)
    end do

  end subroutine mhd_getv_Hall

  subroutine mhd_getdt_Hall(w,x,ixI^L,ixO^L,dx^D,dthall)
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in)    :: dx^D
    double precision, intent(in)    :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out)   :: dthall
    !.. local ..
    double precision :: dxarr(ndim)
    double precision :: bmag(ixI^S)

    dthall=bigdouble

    ! because we have that in cmax now:
    return

    ^D&dxarr(^D)=dx^D;

    if (.not. B0field) then
       bmag(ixO^S)=sqrt(sum(w(ixO^S,mag(:))**2, dim=ndim+1))
       bmag(ixO^S)=sqrt(sum((w(ixO^S,mag(:)) + block%B0(ixO^S,1:ndir,block%iw0))**2))
    end if

    dthall=dtdiffpar*minval(dxarr(1:ndim))**2.0d0/(mhd_etah*maxval(bmag(ixO^S)/w(ixO^S,rho_)))
  end subroutine mhd_getdt_Hall

  !> This implements eq. (42) in Dedner et al. 2002 JcP 175
  !> Gives the Riemann solution on the interface
  !> for the normal B component and Psi in the GLM-MHD system.
  !> 23/04/2013 Oliver Porth
  subroutine glmSolve(wLC,wRC,ixI^L,ixO^L,idir)
    use mod_global_parameters
    double precision, intent(inout) :: wLC(ixI^S,1:nw), wRC(ixI^S,1:nw)
    integer, intent(in)             :: ixI^L, ixO^L, idir
    double precision                :: dB(ixI^S), dPsi(ixI^S)

    dB(ixO^S)   = wRC(ixO^S,mag(idir)) - wLC(ixO^S,mag(idir))
    dPsi(ixO^S) = wRC(ixO^S,psi_) - wLC(ixO^S,psi_)

    wLC(ixO^S,mag(idir))   = 0.5d0 * (wRC(ixO^S,mag(idir)) + wLC(ixO^S,mag(idir))) &
         - half/cmax_global * dPsi(ixO^S)
    wLC(ixO^S,psi_)       = 0.5d0 * (wRC(ixO^S,psi_) + wLC(ixO^S,psi_)) &
         - half*cmax_global * dB(ixO^S)

    wRC(ixO^S,mag(idir)) = wLC(ixO^S,mag(idir))
    wRC(ixO^S,psi_) = wLC(ixO^S,psi_)

  end subroutine glmSolve

  subroutine mhd_boundary_adjust()
    use mod_global_parameters
    integer :: iB, idim, iside, iigrid, igrid 
    integer :: ixG^L, ixO^L, i^D

    ixG^L=ixG^LL;
     do iigrid=1,igridstail; igrid=igrids(iigrid);
        if(.not.phyboundblock(igrid)) cycle
        saveigrid=igrid
        block=>pw(igrid)
        ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
        do idim=1,ndim
           ! to avoid using as yet unknown corner info in more than 1D, we
           ! fill only interior mesh ranges of the ghost cell ranges at first,
           ! and progressively enlarge the ranges to include corners later
           do iside=1,2
              i^D=kr(^D,idim)*(2*iside-3);
              if (neighbor_type(i^D,igrid)/=1) cycle
              iB=(idim-1)*2+iside
              if(.not.boundary_divbfix(iB)) cycle
              if(any(typeboundary(:,iB)=="special")) then
                select case (idim)
                {case (^D)
                   if (iside==2) then
                      ! maximal boundary
                      ixOmin^DD=ixGmax^D+1-nghostcells^D%ixOmin^DD=ixGmin^DD;
                      ixOmax^DD=ixGmax^DD;
                   else
                      ! minimal boundary
                      ixOmin^DD=ixGmin^DD;
                      ixOmax^DD=ixGmin^D-1+nghostcells^D%ixOmax^DD=ixGmax^DD;
                   end if \}
                end select
                call fixdivB_boundary(ixG^L,ixO^L,pw(igrid)%wb,pw(igrid)%x,iB)
              end if
           end do
        end do
     end do

  end subroutine mhd_boundary_adjust

  subroutine fixdivB_boundary(ixG^L,ixO^L,w,x,iB)
    use mod_global_parameters

    integer, intent(in) :: ixG^L,ixO^L,iB
    double precision, intent(inout) :: w(ixG^S,1:nw)
    double precision, intent(in) :: x(ixG^S,1:ndim)

    double precision :: dx1x2,dx1x3,dx2x1,dx2x3,dx3x1,dx3x2
    integer :: ix^D,ixF^L

    select case(iB)
     case(1)
       ! 2nd order CD for divB=0 to set normal B component better
       call mhd_to_primitive(ixG^L,ixO^L,w,x)
       {^IFTWOD
       ixFmin1=ixOmin1+1
       ixFmax1=ixOmax1+1
       ixFmin2=ixOmin2+1
       ixFmax2=ixOmax2-1
       if(slab) then
         dx1x2=dxlevel(1)/dxlevel(2)
         do ix1=ixFmax1,ixFmin1,-1
           w(ix1-1,ixFmin2:ixFmax2,mag(1))=w(ix1+1,ixFmin2:ixFmax2,mag(1)) &
            +dx1x2*(w(ix1,ixFmin2+1:ixFmax2+1,mag(2))-&
                    w(ix1,ixFmin2-1:ixFmax2-1,mag(2)))
         enddo
       else
         do ix1=ixFmax1,ixFmin1,-1
           w(ix1-1,ixFmin2:ixFmax2,mag(1))=( (w(ix1+1,ixFmin2:ixFmax2,mag(1))+&
             w(ix1,ixFmin2:ixFmax2,mag(1)))*block%surfaceC1(ix1,ixFmin2:ixFmax2)&
           +(w(ix1,ixFmin2+1:ixFmax2+1,mag(2))+w(ix1,ixFmin2:ixFmax2,mag(2)))*&
             block%surfaceC2(ix1,ixFmin2:ixFmax2)&
           -(w(ix1,ixFmin2:ixFmax2,mag(2))+w(ix1,ixFmin2-1:ixFmax2-1,mag(2)))*&
             block%surfaceC2(ix1,ixFmin2-1:ixFmax2-1) )&
            /block%surfaceC1(ix1-1,ixFmin2:ixFmax2)-w(ix1,ixFmin2:ixFmax2,mag(1))
         end do
       end if
       }
       {^IFTHREED
       ixFmin1=ixOmin1+1
       ixFmax1=ixOmax1+1
       ixFmin2=ixOmin2+1
       ixFmax2=ixOmax2-1
       ixFmin3=ixOmin3+1
       ixFmax3=ixOmax3-1
       if(slab) then
         dx1x2=dxlevel(1)/dxlevel(2)
         dx1x3=dxlevel(1)/dxlevel(3)
         do ix1=ixFmax1,ixFmin1,-1
           w(ix1-1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(1))=&
                     w(ix1+1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(1)) &
             +dx1x2*(w(ix1,ixFmin2+1:ixFmax2+1,ixFmin3:ixFmax3,mag(2))-&
                     w(ix1,ixFmin2-1:ixFmax2-1,ixFmin3:ixFmax3,mag(2))) &
             +dx1x3*(w(ix1,ixFmin2:ixFmax2,ixFmin3+1:ixFmax3+1,mag(3))-&
                     w(ix1,ixFmin2:ixFmax2,ixFmin3-1:ixFmax3-1,mag(3)))
         end do
       else
         do ix1=ixFmax1,ixFmin1,-1
           w(ix1-1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(1))=&
          ( (w(ix1+1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(1))+&
             w(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(1)))*&
             block%surfaceC1(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3)&
           +(w(ix1,ixFmin2+1:ixFmax2+1,ixFmin3:ixFmax3,mag(2))+&
             w(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(2)))*&
             block%surfaceC2(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3)&
           -(w(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(2))+&
             w(ix1,ixFmin2-1:ixFmax2-1,ixFmin3:ixFmax3,mag(2)))*&
             block%surfaceC2(ix1,ixFmin2-1:ixFmax2-1,ixFmin3:ixFmax3)&
           +(w(ix1,ixFmin2:ixFmax2,ixFmin3+1:ixFmax3+1,mag(3))+&
             w(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(3)))*&
             block%surfaceC3(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3)&
           -(w(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(3))+&
             w(ix1,ixFmin2:ixFmax2,ixFmin3-1:ixFmax3-1,mag(3)))*&
             block%surfaceC3(ix1,ixFmin2:ixFmax2,ixFmin3-1:ixFmax3-1) )&
            /block%surfaceC1(ix1-1,ixFmin2:ixFmax2,ixFmin3:ixFmax3)-&
             w(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(1))
         end do
       end if
       }
       call mhd_to_conserved(ixG^L,ixO^L,w,x)
     case(2)
       call mhd_to_primitive(ixG^L,ixO^L,w,x)
       {^IFTWOD
       ixFmin1=ixOmin1-1
       ixFmax1=ixOmax1-1
       ixFmin2=ixOmin2+1
       ixFmax2=ixOmax2-1
       if(slab) then
         dx1x2=dxlevel(1)/dxlevel(2)
         do ix1=ixFmin1,ixFmax1
           w(ix1+1,ixFmin2:ixFmax2,mag(1))=w(ix1-1,ixFmin2:ixFmax2,mag(1)) &
            -dx1x2*(w(ix1,ixFmin2+1:ixFmax2+1,mag(2))-&
                    w(ix1,ixFmin2-1:ixFmax2-1,mag(2)))
         enddo
       else
         do ix1=ixFmin1,ixFmax1
           w(ix1+1,ixFmin2:ixFmax2,mag(1))=( (w(ix1-1,ixFmin2:ixFmax2,mag(1))+&
             w(ix1,ixFmin2:ixFmax2,mag(1)))*block%surfaceC1(ix1-1,ixFmin2:ixFmax2)&
           -(w(ix1,ixFmin2+1:ixFmax2+1,mag(2))+w(ix1,ixFmin2:ixFmax2,mag(2)))*&
             block%surfaceC2(ix1,ixFmin2:ixFmax2)&
           +(w(ix1,ixFmin2:ixFmax2,mag(2))+w(ix1,ixFmin2-1:ixFmax2-1,mag(2)))*&
             block%surfaceC2(ix1,ixFmin2-1:ixFmax2-1) )&
            /block%surfaceC1(ix1,ixFmin2:ixFmax2)-w(ix1,ixFmin2:ixFmax2,mag(1))
         end do
       end if
       }
       {^IFTHREED
       ixFmin1=ixOmin1-1
       ixFmax1=ixOmax1-1
       ixFmin2=ixOmin2+1
       ixFmax2=ixOmax2-1
       ixFmin3=ixOmin3+1
       ixFmax3=ixOmax3-1
       if(slab) then
         dx1x2=dxlevel(1)/dxlevel(2)
         dx1x3=dxlevel(1)/dxlevel(3)
         do ix1=ixFmin1,ixFmax1
           w(ix1+1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(1))=&
                     w(ix1-1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(1)) &
             -dx1x2*(w(ix1,ixFmin2+1:ixFmax2+1,ixFmin3:ixFmax3,mag(2))-&
                     w(ix1,ixFmin2-1:ixFmax2-1,ixFmin3:ixFmax3,mag(2))) &
             -dx1x3*(w(ix1,ixFmin2:ixFmax2,ixFmin3+1:ixFmax3+1,mag(3))-&
                     w(ix1,ixFmin2:ixFmax2,ixFmin3-1:ixFmax3-1,mag(3)))
         end do
       else
         do ix1=ixFmin1,ixFmax1
           w(ix1+1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(1))=&
          ( (w(ix1-1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(1))+&
             w(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(1)))*&
             block%surfaceC1(ix1-1,ixFmin2:ixFmax2,ixFmin3:ixFmax3)&
           -(w(ix1,ixFmin2+1:ixFmax2+1,ixFmin3:ixFmax3,mag(2))+&
             w(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(2)))*&
             block%surfaceC2(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3)&
           +(w(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(2))+&
             w(ix1,ixFmin2-1:ixFmax2-1,ixFmin3:ixFmax3,mag(2)))*&
             block%surfaceC2(ix1,ixFmin2-1:ixFmax2-1,ixFmin3:ixFmax3)&
           -(w(ix1,ixFmin2:ixFmax2,ixFmin3+1:ixFmax3+1,mag(3))+&
             w(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(3)))*&
             block%surfaceC3(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3)&
           +(w(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(3))+&
             w(ix1,ixFmin2:ixFmax2,ixFmin3-1:ixFmax3-1,mag(3)))*&
             block%surfaceC3(ix1,ixFmin2:ixFmax2,ixFmin3-1:ixFmax3-1) )&
            /block%surfaceC1(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3)-&
             w(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(1))
         end do
       end if
       }
       call mhd_to_conserved(ixG^L,ixO^L,w,x)
     case(3)
       call mhd_to_primitive(ixG^L,ixO^L,w,x)
       {^IFTWOD
       ixFmin1=ixOmin1+1
       ixFmax1=ixOmax1-1
       ixFmin2=ixOmin2+1
       ixFmax2=ixOmax2+1
       if(slab) then
         dx2x1=dxlevel(2)/dxlevel(1)
         do ix2=ixFmax2,ixFmin2,-1
           w(ixFmin1:ixFmax1,ix2-1,mag(2))=w(ixFmin1:ixFmax1,ix2+1,mag(2)) &
            +dx2x1*(w(ixFmin1+1:ixFmax1+1,ix2,mag(1))-&
                    w(ixFmin1-1:ixFmax1-1,ix2,mag(1)))
         enddo
       else
         do ix2=ixFmax2,ixFmin2,-1
           w(ixFmin1:ixFmax1,ix2-1,mag(2))=( (w(ixFmin1:ixFmax1,ix2+1,mag(2))+&
             w(ixFmin1:ixFmax1,ix2,mag(2)))*block%surfaceC2(ixFmin1:ixFmax1,ix2)&
           +(w(ixFmin1+1:ixFmax1+1,ix2,mag(1))+w(ixFmin1:ixFmax1,ix2,mag(1)))*&
             block%surfaceC1(ixFmin1:ixFmax1,ix2)&
           -(w(ixFmin1:ixFmax1,ix2,mag(1))+w(ixFmin1-1:ixFmax1-1,ix2,mag(1)))*&
             block%surfaceC1(ixFmin1-1:ixFmax1-1,ix2) )&
            /block%surfaceC2(ixFmin1:ixFmax1,ix2-1)-w(ixFmin1:ixFmax1,ix2,mag(2))
         end do
       end if
       }
       {^IFTHREED
       ixFmin1=ixOmin1+1
       ixFmax1=ixOmax1-1
       ixFmin3=ixOmin3+1
       ixFmax3=ixOmax3-1
       ixFmin2=ixOmin2+1
       ixFmax2=ixOmax2+1
       if(slab) then
         dx2x1=dxlevel(2)/dxlevel(1)
         dx2x3=dxlevel(2)/dxlevel(3)
         do ix2=ixFmax2,ixFmax2,-1
           w(ixFmin1:ixFmax1,ix2-1,ixFmin3:ixFmax3,mag(2))=w(ixFmin1:ixFmax1,&
             ix2+1,ixFmin3:ixFmax3,mag(2)) &
             +dx2x1*(w(ixFmin1+1:ixFmax1+1,ix2,ixFmin3:ixFmax3,mag(1))-&
                     w(ixFmin1-1:ixFmax1-1,ix2,ixFmin3:ixFmax3,mag(1))) &
             +dx2x3*(w(ixFmin1:ixFmax1,ix2,ixFmin3+1:ixFmax3+1,mag(3))-&
                     w(ixFmin1:ixFmax1,ix2,ixFmin3-1:ixFmax3-1,mag(3)))
         end do
       else
         do ix2=ixFmax2,ixFmin2,-1
           w(ixFmin1:ixFmax1,ix2-1,ixFmin3:ixFmax3,mag(2))=&
          ( (w(ixFmin1:ixFmax1,ix2+1,ixFmin3:ixFmax3,mag(2))+&
             w(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,mag(2)))*&
             block%surfaceC2(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3)&
           +(w(ixFmin1+1:ixFmax1+1,ix2,ixFmin3:ixFmax3,mag(1))+&
             w(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,mag(1)))*&
             block%surfaceC1(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3)&
           -(w(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,mag(1))+&
             w(ixFmin1-1:ixFmax1-1,ix2,ixFmin3:ixFmax3,mag(1)))*&
             block%surfaceC1(ixFmin1-1:ixFmax1-1,ix2,ixFmin3:ixFmax3)&
           +(w(ixFmin1:ixFmax1,ix2,ixFmin3+1:ixFmax3+1,mag(3))+&
             w(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,mag(3)))*&
             block%surfaceC3(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3)&
           -(w(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,mag(3))+&
             w(ixFmin1:ixFmax1,ix2,ixFmin3-1:ixFmax3-1,mag(3)))*&
             block%surfaceC3(ixFmin1:ixFmax1,ix2,ixFmin3-1:ixFmax3-1) )&
            /block%surfaceC2(ixFmin1:ixFmax1,ix2-1,ixFmin3:ixFmax3)-&
             w(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,mag(2))
         end do
       end if
       }
       call mhd_to_conserved(ixG^L,ixO^L,w,x)
     case(4)
       call mhd_to_primitive(ixG^L,ixO^L,w,x)
       {^IFTWOD
       ixFmin1=ixOmin1+1
       ixFmax1=ixOmax1-1
       ixFmin2=ixOmin2-1
       ixFmax2=ixOmax2-1
       if(slab) then
         dx2x1=dxlevel(2)/dxlevel(1)
         do ix2=ixFmin2,ixFmax2
           w(ixFmin1:ixFmax1,ix2+1,mag(2))=w(ixFmin1:ixFmax1,ix2-1,mag(2)) &
            -dx2x1*(w(ixFmin1+1:ixFmax1+1,ix2,mag(1))-&
                    w(ixFmin1-1:ixFmax1-1,ix2,mag(1)))
         end do
       else
         do ix2=ixFmin2,ixFmax2
           w(ixFmin1:ixFmax1,ix2+1,mag(2))=( (w(ixFmin1:ixFmax1,ix2-1,mag(2))+&
             w(ixFmin1:ixFmax1,ix2,mag(2)))*block%surfaceC2(ixFmin1:ixFmax1,ix2-1)&
           -(w(ixFmin1+1:ixFmax1+1,ix2,mag(1))+w(ixFmin1:ixFmax1,ix2,mag(1)))*&
             block%surfaceC1(ixFmin1:ixFmax1,ix2)&
           +(w(ixFmin1:ixFmax1,ix2,mag(1))+w(ixFmin1-1:ixFmax1-1,ix2,mag(1)))*&
             block%surfaceC1(ixFmin1-1:ixFmax1-1,ix2) )&
            /block%surfaceC2(ixFmin1:ixFmax1,ix2)-w(ixFmin1:ixFmax1,ix2,mag(2))
         end do
       end if
       }
       {^IFTHREED
       ixFmin1=ixOmin1+1
       ixFmax1=ixOmax1-1
       ixFmin3=ixOmin3+1
       ixFmax3=ixOmax3-1
       ixFmin2=ixOmin2-1
       ixFmax2=ixOmax2-1
       if(slab) then
         dx2x1=dxlevel(2)/dxlevel(1)
         dx2x3=dxlevel(2)/dxlevel(3)
         do ix2=ixFmin2,ixFmax2
           w(ixFmin1:ixFmax1,ix2+1,ixFmin3:ixFmax3,mag(2))=w(ixFmin1:ixFmax1,&
             ix2-1,ixFmin3:ixFmax3,mag(2)) &
             -dx2x1*(w(ixFmin1+1:ixFmax1+1,ix2,ixFmin3:ixFmax3,mag(1))-&
                     w(ixFmin1-1:ixFmax1-1,ix2,ixFmin3:ixFmax3,mag(1))) &
             -dx2x3*(w(ixFmin1:ixFmax1,ix2,ixFmin3+1:ixFmax3+1,mag(3))-&
                     w(ixFmin1:ixFmax1,ix2,ixFmin3-1:ixFmax3-1,mag(3)))
         end do
       else
         do ix2=ixFmin2,ixFmax2
           w(ixFmin1:ixFmax1,ix2+1,ixFmin3:ixFmax3,mag(2))=&
          ( (w(ixFmin1:ixFmax1,ix2-1,ixFmin3:ixFmax3,mag(2))+&
             w(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,mag(2)))*&
             block%surfaceC2(ixFmin1:ixFmax1,ix2-1,ixFmin3:ixFmax3)&
           -(w(ixFmin1+1:ixFmax1+1,ix2,ixFmin3:ixFmax3,mag(1))+&
             w(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,mag(1)))*&
             block%surfaceC1(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3)&
           +(w(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,mag(1))+&
             w(ixFmin1-1:ixFmax1-1,ix2,ixFmin3:ixFmax3,mag(1)))*&
             block%surfaceC1(ixFmin1-1:ixFmax1-1,ix2,ixFmin3:ixFmax3)&
           -(w(ixFmin1:ixFmax1,ix2,ixFmin3+1:ixFmax3+1,mag(3))+&
             w(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,mag(3)))*&
             block%surfaceC3(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3)&
           +(w(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,mag(3))+&
             w(ixFmin1:ixFmax1,ix2,ixFmin3-1:ixFmax3-1,mag(3)))*&
             block%surfaceC3(ixFmin1:ixFmax1,ix2,ixFmin3-1:ixFmax3-1) )&
            /block%surfaceC2(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3)-&
             w(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,mag(2))
         end do
       end if
       }
       call mhd_to_conserved(ixG^L,ixO^L,w,x)
     {^IFTHREED
     case(5)
       call mhd_to_primitive(ixG^L,ixO^L,w,x)
       ixFmin1=ixOmin1+1
       ixFmax1=ixOmax1-1
       ixFmin2=ixOmin2+1
       ixFmax2=ixOmax2-1
       ixFmin3=ixOmin3+1
       ixFmax3=ixOmax3+1
       if(slab) then
         dx3x1=dxlevel(3)/dxlevel(1)
         dx3x2=dxlevel(3)/dxlevel(2)
         do ix3=ixFmax3,ixFmax3,-1
           w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3-1,mag(3))=w(ixFmin1:ixFmax1,&
             ixFmin2:ixFmax2,ix3+1,mag(3)) &
             +dx3x1*(w(ixFmin1+1:ixFmax1+1,ixFmin2:ixFmax2,ix3,mag(1))-&
                     w(ixFmin1-1:ixFmax1-1,ixFmin2:ixFmax2,ix3,mag(1))) &
             +dx3x2*(w(ixFmin1:ixFmax1,ixFmin2+1:ixFmax2+1,ix3,mag(2))-&
                     w(ixFmin1:ixFmax1,ixFmin2-1:ixFmax2-1,ix3,mag(2)))
         end do
       else
         do ix3=ixFmax3,ixFmin3,-1
           w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3-1,mag(3))=&
          ( (w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3+1,mag(3))+&
             w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,mag(3)))*&
             block%surfaceC3(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3)&
           +(w(ixFmin1+1:ixFmax1+1,ixFmin2:ixFmax2,ix3,mag(1))+&
             w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,mag(1)))*&
             block%surfaceC1(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3)&
           -(w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,mag(1))+&
             w(ixFmin1-1:ixFmax1-1,ixFmin2:ixFmax2,ix3,mag(1)))*&
             block%surfaceC1(ixFmin1-1:ixFmax1-1,ixFmin2:ixFmax2,ix3)&
           +(w(ixFmin1:ixFmax1,ixFmin2+1:ixFmax2+1,ix3,mag(2))+&
             w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,mag(2)))*&
             block%surfaceC2(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3)&
           -(w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,mag(2))+&
             w(ixFmin1:ixFmax1,ixFmin2-1:ixFmax2-1,ix3,mag(2)))*&
             block%surfaceC2(ixFmin1:ixFmax1,ixFmin2-1:ixFmax2-1,ix3) )&
            /block%surfaceC3(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3-1)-&
             w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,mag(3))
         end do
       end if
       call mhd_to_conserved(ixG^L,ixO^L,w,x)
     case(6)
       call mhd_to_primitive(ixG^L,ixO^L,w,x)
       ixFmin1=ixOmin1+1
       ixFmax1=ixOmax1-1
       ixFmin2=ixOmin2+1
       ixFmax2=ixOmax2-1
       ixFmin3=ixOmin3-1
       ixFmax3=ixOmax3-1
       if(slab) then
         dx3x1=dxlevel(3)/dxlevel(1)
         dx3x2=dxlevel(3)/dxlevel(2)
         do ix3=ixFmin3,ixFmax3
           w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3+1,mag(3))=w(ixFmin1:ixFmax1,&
             ixFmin2:ixFmax2,ix3-1,mag(3)) &
             -dx3x1*(w(ixFmin1+1:ixFmax1+1,ixFmin2:ixFmax2,ix3,mag(1))-&
                     w(ixFmin1-1:ixFmax1-1,ixFmin2:ixFmax2,ix3,mag(1))) &
             -dx3x2*(w(ixFmin1:ixFmax1,ixFmin2+1:ixFmax2+1,ix3,mag(2))-&
                     w(ixFmin1:ixFmax1,ixFmin2-1:ixFmax2-1,ix3,mag(2)))
         end do
       else
         do ix3=ixFmin3,ixFmax3
           w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3+1,mag(3))=&
          ( (w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3-1,mag(3))+&
             w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,mag(3)))*&
             block%surfaceC3(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3-1)&
           -(w(ixFmin1+1:ixFmax1+1,ixFmin2:ixFmax2,ix3,mag(1))+&
             w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,mag(1)))*&
             block%surfaceC1(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3)&
           +(w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,mag(1))+&
             w(ixFmin1-1:ixFmax1-1,ixFmin2:ixFmax2,ix3,mag(1)))*&
             block%surfaceC1(ixFmin1-1:ixFmax1-1,ixFmin2:ixFmax2,ix3)&
           -(w(ixFmin1:ixFmax1,ixFmin2+1:ixFmax2+1,ix3,mag(2))+&
             w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,mag(2)))*&
             block%surfaceC2(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3)&
           +(w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,mag(2))+&
             w(ixFmin1:ixFmax1,ixFmin2-1:ixFmax2-1,ix3,mag(2)))*&
             block%surfaceC2(ixFmin1:ixFmax1,ixFmin2-1:ixFmax2-1,ix3) )&
            /block%surfaceC3(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3)-&
             w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,mag(3))
         end do
       end if
       call mhd_to_conserved(ixG^L,ixO^L,w,x)
     }
     case default
       call mpistop("Special boundary is not defined for this region")
    end select
  
   end subroutine fixdivB_boundary

end module mod_mhd_phys
