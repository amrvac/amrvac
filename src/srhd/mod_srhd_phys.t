!> Special Relativistic Hydrodynamics (with EOS) physics module
module mod_srhd_phys
  use mod_physics
  use mod_constants
  use mod_comm_lib, only: mpistop
  implicit none
  private

  !> Whether particles module is added
  logical, public, protected              :: srhd_particles = .false.

  !> Number of tracer species
  integer, public,protected              :: srhd_n_tracer = 0

  !> Index of the density (in the w array)
  integer, public,protected              :: rho_
  integer, public,protected              :: d_ 

  !> Indices of the momentum density
  integer, allocatable, public, protected :: mom(:)

  !> Indices of the tracers
  integer, allocatable, public, protected :: tracer(:)

  !> Index of the energy density 
  integer, public,protected              :: e_
  !> Index of the gas pressure should equal e_
  integer, public,protected              :: p_

  !> Index of the Lorentz factor
  integer, public,protected     :: lfac_

  !> Index of the inertia
  integer, public,protected     :: xi_

  !> Whether synge eos is used
  logical, public                         :: srhd_eos = .false.

  !> The adiabatic index and derived values
  double precision, public                :: srhd_gamma = 5.d0/3.0d0
  double precision, public                :: gamma_1,inv_gamma_1,gamma_to_gamma_1

  !> The smallest allowed energy
  double precision, public             :: small_e
  !> The smallest allowed inertia
  double precision, public             :: small_xi

  !> Allows overruling default corner filling (for debug mode, otherwise corner primitives fail)
  logical, public, protected              :: srhd_force_diagonal = .false.

  !> Helium abundance over Hydrogen
  double precision, public, protected  :: He_abundance=0.0d0

  !> parameters for NR in con2prim
  integer, public                  :: maxitnr   = 100
  double precision, public         :: absaccnr  = 1.0d-8
  double precision, public         :: tolernr   = 1.0d-9
  double precision, public         :: dmaxvel   = 1.0d-7
  double precision, public         :: lfacmax
  double precision, public :: minp, minrho, smalltau, smallxi

  ! Public methods
  public :: srhd_phys_init
  public :: srhd_get_pthermal
  public :: srhd_get_auxiliary
  public :: srhd_get_auxiliary_prim
  public :: srhd_get_csound2
  public :: srhd_to_conserved
  public :: srhd_to_primitive
  public :: srhd_check_params
  public :: srhd_check_w 
  public :: srhd_get_Geff_eos
  public :: srhd_get_enthalpy_eos
contains

  !> Read this module's parameters from a file
  subroutine srhd_read_params(files)
    use mod_global_parameters
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /srhd_list/ srhd_n_tracer, srhd_eos, srhd_gamma, &
                    srhd_particles, srhd_force_diagonal, &
                    SI_unit, He_abundance

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, srhd_list, end=111)
111    close(unitpar)
    end do

  end subroutine srhd_read_params

  !> Write this module's parameters to a snapshot
  subroutine srhd_write_info(fh)
    use mod_global_parameters
    integer, intent(in)                 :: fh
    integer, parameter                  :: n_par = 1
    double precision                    :: values(n_par)
    character(len=name_len)             :: names(n_par)
    integer, dimension(MPI_STATUS_SIZE) :: st
    integer                             :: er

    call MPI_FILE_WRITE(fh, n_par, 1, MPI_INTEGER, st, er)

    names(1) = "gamma"
    values(1) = srhd_gamma
    call MPI_FILE_WRITE(fh, values, n_par, MPI_DOUBLE_PRECISION, st, er)
    call MPI_FILE_WRITE(fh, names, n_par * name_len, MPI_CHARACTER, st, er)

  end subroutine srhd_write_info

  !> Initialize the module
  subroutine srhd_phys_init()
    use mod_global_parameters
    use mod_particles, only: particles_init
    integer :: itr,idir

    call srhd_read_params(par_files)

    physics_type = "srhd"
    phys_energy  = .true.
    phys_total_energy  = .true.
    phys_gamma = srhd_gamma

    ! unused physics options
    phys_internal_e = .false.
    phys_partial_ionization=.false.
    phys_trac=.false.

    use_particles = srhd_particles

    ! note: number_species is 1 for srhd
    allocate(start_indices(number_species),stop_indices(number_species))
    ! set the index of the first flux variable for species 1
    start_indices(1)=1

    ! Determine flux variables
    rho_ = var_set_rho()
    d_=rho_

    allocate(mom(ndir))
    mom(:) = var_set_momentum(ndir)

    ! Set index of energy variable
    e_ = var_set_energy()
    p_ = e_

    ! Whether diagonal ghost cells are required for the physics
    phys_req_diagonal = .false.

    ! derive units from basic units
    call srhd_physical_units()

    if (srhd_force_diagonal) then
       ! ensure corners are filled, otherwise divide by zero when getting primitives
       !  --> only for debug purposes
       phys_req_diagonal = .true.
    endif

    allocate(tracer(srhd_n_tracer))

    ! Set starting index of tracers
    do itr = 1, srhd_n_tracer
       tracer(itr) = var_set_fluxvar("trc", "trp", itr, need_bc=.false.)
    end do

    ! Set index for auxiliary variables
    ! MUST be after the possible tracers (which have fluxes)
    xi_  = var_set_auxvar('xi','xi')
    lfac_= var_set_auxvar('lfac','lfac')

    ! set number of variables which need update ghostcells
    nwgc=nwflux+nwaux

    ! set the index of the last flux variable for species 1
    stop_indices(1)=nwflux

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

    ! dummy for now, no extra source terms precoded
    phys_add_source          => srhd_add_source
    phys_get_dt              => srhd_get_dt
    ! copied in from HD/MHD, for certain limiters
    phys_get_a2max           => srhd_get_a2max

    ! actual srhd routines
    phys_check_params        => srhd_check_params
    phys_check_w             => srhd_check_w
    phys_get_cmax            => srhd_get_cmax
    phys_get_cbounds         => srhd_get_cbounds
    phys_get_flux            => srhd_get_flux
    phys_add_source_geom     => srhd_add_source_geom
    phys_to_conserved        => srhd_to_conserved
    phys_to_primitive        => srhd_to_primitive
    phys_get_pthermal        => srhd_get_pthermal
    phys_get_auxiliary       => srhd_get_auxiliary
    phys_get_auxiliary_prim  => srhd_get_auxiliary_prim
    phys_get_pthermal        => srhd_get_pthermal
    phys_get_v               => srhd_get_v
    phys_write_info          => srhd_write_info
    phys_handle_small_values => srhd_handle_small_values

    ! Initialize particles module
    if (srhd_particles) then
       call particles_init()
       phys_req_diagonal = .true.
    end if

  end subroutine srhd_phys_init

  subroutine srhd_check_params
    use mod_global_parameters

    if (srhd_gamma <= 0.0d0 .or. srhd_gamma == 1.0d0) &
            call mpistop ("Error: srhd_gamma <= 0 or srhd_gamma == 1")
    ! additional useful values
    gamma_1=srhd_gamma-1.0d0
    inv_gamma_1=1.0d0/gamma_1
    gamma_to_gamma_1=srhd_gamma/gamma_1
   
    ! the following sets small_e and small_xi from small_density/small_pressure
    ! according to the srhd_eos used
    call srhd_get_smallvalues_eos

    if(mype==0)then
       write(*,*)'------------------------------------------------------------'
       write(*,*)'Using EOS set via srhd_eos=',srhd_eos
       write(*,*)'Maximal lorentz factor (via dmaxvel) is=',lfacmax
       write(*,*)'Use fixes set through check/fix small values:', check_small_values,fix_small_values
       write(*,*)'Controlled with small pressure/density:', small_pressure,small_density
       write(*,*)'Derived small values: xi and e ',small_xi,small_e
       write(*,*)'------------------------------------------------------------'
    endif

  end subroutine srhd_check_params

  subroutine srhd_physical_units
    use mod_global_parameters
    double precision :: mp,kB

    ! Derive scaling units
    if(SI_unit) then
      mp=mp_SI
      kB=kB_SI
      unit_velocity=c_SI
    else
      mp=mp_cgs
      kB=kB_cgs
      unit_velocity=const_c
    end if
    if(unit_numberdensity*unit_length<=0.0d0)then
       call mpistop("Abort: must set positive values for unit length and numberdensity")
    endif
    ! we assume user sets: unit_numberdensity, unit_length, He_abundance
    ! then together with light speed c, all units fixed
    unit_density=(1.0d0+4.0d0*He_abundance)*mp*unit_numberdensity
    unit_pressure=unit_density*unit_velocity**2
    unit_temperature=unit_pressure/((2.0d0+3.0d0*He_abundance)*unit_numberdensity*kB)
    unit_time=unit_length/unit_velocity
    unit_mass = unit_density*unit_length**3

  end subroutine srhd_physical_units

  !> Returns logical argument flag T where values are not ok
  subroutine srhd_check_w(primitive, ixI^L, ixO^L, w, flag)
    use mod_global_parameters
    logical, intent(in)          :: primitive
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, nw)
    logical, intent(inout)       :: flag(ixI^S,1:nw)

    flag=.false.

    ! NOTE: we should not check or use nwaux variables here
    if(primitive) then
       where(w(ixO^S,rho_) < small_density) flag(ixO^S,rho_) = .true.
       where(w(ixO^S,p_) < small_pressure) flag(ixO^S,p_) = .true.
    else
       where(w(ixO^S,d_) < small_density) flag(ixO^S,d_) = .true.
       where(w(ixO^S,e_) < small_e) flag(ixO^S,e_) = .true.
    endif

  end subroutine srhd_check_w

  !> Returns logical argument flag T where auxiliary values are not ok
  subroutine srhd_check_w_aux(ixI^L, ixO^L, w, flag)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, nw)
    logical, intent(inout)       :: flag(ixI^S,1:nw)

    flag=.false.

    where(w(ixO^S,xi_) < small_xi) flag(ixO^S,xi_) = .true.
    where(w(ixO^S,lfac_) < one) flag(ixO^S,lfac_) = .true.

    if(any(flag(ixO^S,xi_)))then
       write(*,*)'auxiliary xi too low: abort'
       call mpistop('auxiliary  check failed')
    endif
    if(any(flag(ixO^S,lfac_)))then
       write(*,*)'auxiliary lfac too low: abort'
       call mpistop('auxiliary  check failed')
    endif

  end subroutine srhd_check_w_aux

  !> Set auxiliary variables lfac and xi from a primitive state
  !> only used when handle_small_values average on primitives
  subroutine srhd_get_auxiliary_prim(ixI^L,ixO^L,w)
    use mod_global_parameters
    integer, intent(in)                :: ixI^L, ixO^L
    double precision, intent(inout)    :: w(ixI^S, nw)
    double precision, dimension(ixO^S) :: rho,rhoh,pth

    ! assume four-velocity in momentum vector (i.e. lfac*v)
    rho(ixO^S)    = sum(w(ixO^S, mom(:))**2, dim=ndim+1)
    w(ixO^S,lfac_) = dsqrt(1.0d0+rho(ixO^S))

    rho(ixO^S)=w(ixO^S,rho_)
    pth(ixO^S)=w(ixO^S,p_)
    ! compute rho*h (enthalpy h) from density-pressure
    call srhd_get_enthalpy_eos(ixO^L,rho,pth,rhoh)

    ! fill auxiliary variable xi= lfac^2 rhoh
    w(ixO^S,xi_) = w(ixO^S,lfac_)**2.0d0*rhoh(ixO^S)

  end subroutine srhd_get_auxiliary_prim

  !> Compute auxiliary variables lfac and xi from a conservative state
  !> using srhd_con2prim to calculate enthalpy and lorentz factor
  subroutine srhd_get_auxiliary(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    implicit none

    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)

    integer                        :: ix^D,ierror,idir
    integer                        :: flag_error(ixO^S)
    double precision               :: ssqr

    if(srhd_eos)then
      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        ierror=0
        ssqr=0.0d0
        do idir=1,ndir
          ssqr= ssqr+w(ix^D,mom(idir))**2
        enddo
        if(w(ix^D,d_)<small_density)then
           print *,'entering con2prim with', w(ix^D,lfac_),w(ix^D,xi_), &
                    w(ix^D,d_),ssqr,w(ix^D,e_)
           print *,'in position',ix^D,x(ix^D,1:ndim)
           print *,small_density,small_e,small_pressure,small_xi
           if(check_small_values) call mpistop('small density on entry con2prim')
           if(fix_small_values) w(ix^D,d_)=small_density
        endif
        if(w(ix^D,e_)<small_e)then
           print *,'entering con2prim with', w(ix^D,lfac_),w(ix^D,xi_), &
                    w(ix^D,d_),ssqr,w(ix^D,e_)
           print *,'in position',ix^D,x(ix^D,1:ndim)
           print *,small_density,small_e,small_pressure,small_xi
           if(check_small_values) call mpistop('small energy on entry con2prim')
           if(fix_small_values) w(ix^D,e_)=small_e
        endif
        call con2prim_eos(w(ix^D,lfac_),w(ix^D,xi_), &
                    w(ix^D,d_),ssqr,w(ix^D,e_),ierror)
        flag_error(ix^D) = ierror
      {enddo^D&\}
    else
      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        ierror=0
        ssqr=0.0d0
        do idir=1,ndir
          ssqr= ssqr+w(ix^D,mom(idir))**2
        enddo
        if(w(ix^D,d_)<small_density)then
           print *,'entering con2prim with', w(ix^D,lfac_),w(ix^D,xi_), &
                    w(ix^D,d_),ssqr,w(ix^D,e_)
           print *,'in position',ix^D,x(ix^D,1:ndim)
           print *,small_density,small_e,small_pressure,small_xi
           if(check_small_values) call mpistop('small density on entry con2prim')
           if(fix_small_values) w(ix^D,d_)=small_density
        endif
        if(w(ix^D,e_)<small_e)then
           print *,'entering con2prim with', w(ix^D,lfac_),w(ix^D,xi_), &
                    w(ix^D,d_),ssqr,w(ix^D,e_)
           print *,'in position',ix^D,x(ix^D,1:ndim)
           print *,small_density,small_e,small_pressure,small_xi
           if(check_small_values) call mpistop('small energy on entry con2prim')
           if(fix_small_values) w(ix^D,e_)=small_e
        endif
        call con2prim(w(ix^D,lfac_),w(ix^D,xi_), &
                    w(ix^D,d_),ssqr,w(ix^D,e_),ierror)
        flag_error(ix^D) = ierror
      {enddo^D&\}
    endif

    if(check_small_values)then
     if(any(flag_error(ixO^S)/=0))then
         print *,flag_error(ixO^S)
         call mpistop('Problem when getting auxiliaries')
    !    call srhd_handle_small_values(.false.,w,x,ixI^L,ixO^L,'srhd_get_auxiliary')
     end if 
    end if 

  end subroutine srhd_get_auxiliary

  !> Transform primitive variables into conservative ones
  subroutine srhd_to_conserved(ixI^L, ixO^L, w, x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    integer                         :: idir,itr
    double precision, dimension(ixO^S) :: rhoh,rho,pth

    ! assume four-velocity in momentum vector (i.e. lfac*v)
    ! use rhoh slot for temporary array
    rhoh(ixO^S)    = sum(w(ixO^S, mom(:))**2, dim=ndim+1)
    w(ixO^S,lfac_) = dsqrt(1.0d0+rhoh(ixO^S))

    rho(ixO^S)=w(ixO^S,rho_)
    pth(ixO^S)=w(ixO^S,p_)
    ! compute rho*h (enthalpy h) from density-pressure
    call srhd_get_enthalpy_eos(ixO^L,rho,pth,rhoh)

    ! compute rhoh*lfac (recycle rhoh)
    rhoh(ixO^S)= rhoh(ixO^S)*w(ixO^S,lfac_)
    ! fill auxiliary variable xi= lfac^2 rhoh
    w(ixO^S,xi_) = w(ixO^S,lfac_)*rhoh(ixO^S)

    ! set conservative density: d = lfac * rho
    w(ixO^S,d_)=w(ixO^S,lfac_)*rho(ixO^S)

    ! Convert four-velocity (lfac*v) to momentum (xi*v=[rho*h*lfac^2]*v)
    do idir = 1, ndir
       w(ixO^S, mom(idir)) = rhoh(ixO^S)*w(ixO^S, mom(idir))
    end do 

    ! set tau = xi-p-d energy variable
    w(ixO^S,e_) = w(ixO^S,xi_)-pth(ixO^S)-w(ixO^S,d_)

    do itr=1,srhd_n_tracer
       w(ixO^S,tracer(itr)) = w(ixO^S,d_)*w(ixO^S,tracer(itr))
    end do

  end subroutine srhd_to_conserved

  !> Transform conservative variables into primitive ones
  subroutine srhd_to_primitive(ixI^L, ixO^L, w, x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    integer                         :: idir,itr
    double precision, dimension(ixO^S) :: rho,rhoh,E
    double precision, dimension(ixI^S) :: pth
    character(len=30)                  :: subname_loc

    ! get auxiliary variables lfac and xi from conserved set
    call srhd_get_auxiliary(ixI^L,ixO^L,w,x)

    ! from d to rho (d=rho*lfac)
    rho(ixO^S) = w(ixO^S,d_)/w(ixO^S,lfac_)

    ! compute pressure
    ! deduce rho*h from xi/lfac^2
    rhoh(ixO^S) = w(ixO^S,xi_)/w(ixO^S,lfac_)**2.0d0
    call srhd_get_pressure_eos(ixI^L,ixO^L,rho,rhoh,pth,E)

    w(ixO^S,rho_)=rho(ixO^S)
    ! from xi*v to U=lfac*v (four-velocity) 
    do idir=1,ndir
      w(ixO^S,mom(idir)) = w(ixO^S,lfac_)*w(ixO^S,mom(idir))&
                           /w(ixO^S,xi_)
    end do
    w(ixO^S,p_)=pth(ixO^S)

    do itr=1,srhd_n_tracer
       w(ixO^S,tracer(itr)) = w(ixO^S,tracer(itr)) &
                               /(rho(ixO^S)*w(ixO^S,lfac_))
    end do

  end subroutine srhd_to_primitive

  !> Calculate v vector from conservatives
  subroutine srhd_get_v(w,x,ixI^L,ixO^L,v)
    use mod_global_parameters
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S, nw), x(ixI^S,1:ndim)
    double precision, intent(out) :: v(ixI^S,1:ndir)
    integer :: idir

    ! get v from xi*v
    do idir=1,ndir
      v(ixO^S,idir) = w(ixO^S, mom(idir))/w(ixO^S,xi_)
    end do 

  end subroutine srhd_get_v

  !> Calculate the square of the thermal sound speed csound2 within ixO^L.
  !> here computed from conservative set WITH ADDED rho*h
  !> local version: does not do con2prim
  subroutine srhd_get_csound2_rhoh(w,x,ixI^L,ixO^L,rhoh,csound2)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: w(ixI^S,nw),rhoh(ixO^S)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out)   :: csound2(ixO^S)

    double precision                :: rho(ixO^S)

    rho=w(ixO^S,d_)/w(ixO^S,lfac_)
    call srhd_get_csound2_eos(ixI^L,ixO^L,rho,rhoh,csound2)

  end subroutine srhd_get_csound2_rhoh

  !> Calculate the square of the thermal sound speed csound2 within ixO^L.
  !> here computed from conservative set and uses con2prim!!!
  !> public version!
  subroutine srhd_get_csound2(w,x,ixI^L,ixO^L,csound2)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S,nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out)   :: csound2(ixO^S)

    double precision                :: rho(ixO^S),rhoh(ixO^S)

    ! get auxiliary variables lfac and xi from conserved set
    call srhd_get_auxiliary(ixI^L,ixO^L,w,x)
    ! quantify rho, rho*h
    rho  = w(ixO^S,d_)/w(ixO^S,lfac_)
    rhoh = w(ixO^S,xi_)/w(ixO^S,lfac_)**2.0d0
    call srhd_get_csound2_eos(ixI^L,ixO^L,rho,rhoh,csound2)

  end subroutine srhd_get_csound2

  !> Calculate thermal pressure p within ixO^L
  !> must follow after update conservative with auxiliaries
  subroutine srhd_get_pthermal(w, x, ixI^L, ixO^L, pth)
    use mod_global_parameters
    use mod_small_values, only: trace_small_values

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    double precision, intent(out)   :: pth(ixI^S)

    integer                      :: iw, ix^D
    double precision             :: rho(ixO^S),rhoh(ixO^S),E(ixO^S)

    ! quantify rho, rho*h, tau and get pthermal
    rho  = w(ixO^S,d_)/w(ixO^S,lfac_)
    rhoh = w(ixO^S,xi_)/w(ixO^S,lfac_)**2.0d0
    call srhd_get_pressure_eos(ixI^L,ixO^L,rho,rhoh,pth,E)

  end subroutine srhd_get_pthermal

  !> dummy addsource subroutine
  ! w[iw]= w[iw]+qdt*S[wCT, qtC, x] where S is the source based on wCT within ixO
  subroutine srhd_add_source(qdt,dtfactor,ixI^L,ixO^L,wCT,wCTprim,w,x,qsourcesplit,active)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt,dtfactor
    double precision, intent(in)    :: wCT(ixI^S, 1:nw),wCTprim(ixI^S,1:nw), x(ixI^S, 1:ndim)
    double precision, intent(inout) :: w(ixI^S, 1:nw)
    logical, intent(in)             :: qsourcesplit
    logical, intent(inout)          :: active

  end subroutine srhd_add_source

  !> dummy get_dt subroutine
  subroutine srhd_get_dt(w, ixI^L, ixO^L, dtnew, dx^D, x)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: dx^D, x(ixI^S, 1:^ND)
    double precision, intent(in)    :: w(ixI^S, 1:nw)
    double precision, intent(inout) :: dtnew

    dtnew = bigdouble

  end subroutine srhd_get_dt

  !> Calculate cmax_idim within ixO^L
  !> used especially for setdt CFL limit
  subroutine srhd_get_cmax(w, x, ixI^L, ixO^L, idim, cmax)
    use mod_global_parameters

    integer, intent(in)                       :: ixI^L, ixO^L, idim
    double precision, intent(in)              :: w(ixI^S, nw), x(ixI^S, 1:ndim)
    double precision, intent(inout)           :: cmax(ixI^S)

    double precision :: wc(ixI^S,nw)
    double precision, dimension(ixO^S)        :: csound2,tmp1,tmp2,v2
    double precision, dimension(ixI^S)        :: vidim, cmin

    logical       :: flag(ixI^S,1:nw)

    !!call srhd_check_w_aux(ixI^L, ixO^L, w, flag)

    ! input w is in primitive form TODO use it
    wc=w
    call srhd_to_conserved(ixI^L, ixO^L, wc, x)
    ! auxiliaries are filled here
    tmp1(ixO^S)=wc(ixO^S,xi_)/wc(ixO^S,lfac_)**2.0d0
    v2(ixO^S)=1.0d0-1.0d0/wc(ixO^S,lfac_)**2
    call srhd_get_csound2_rhoh(wc,x,ixI^L,ixO^L,tmp1,csound2)
    vidim(ixO^S) = wc(ixO^S, mom(idim))/wc(ixO^S, xi_)
    tmp2(ixO^S)=vidim(ixO^S)**2.0d0
    tmp1(ixO^S)=1.0d0-v2(ixO^S)*csound2(ixO^S) &
                        -tmp2(ixO^S)*(1.0d0-csound2(ixO^S))
    tmp2(ixO^S)=dsqrt(csound2(ixO^S)*(one-v2(ixO^S))*tmp1(ixO^S))
    tmp1(ixO^S)=vidim(ixO^S)*(one-csound2(ixO^S))
    cmax(ixO^S)=(tmp1(ixO^S)+tmp2(ixO^S))/(one-v2(ixO^S)*csound2(ixO^S))
    cmin(ixO^S)=(tmp1(ixO^S)-tmp2(ixO^S))/(one-v2(ixO^S)*csound2(ixO^S))
    ! Limit by speed of light
    cmin(ixO^S) = max(cmin(ixO^S), - 1.0d0)
    cmin(ixO^S) = min(cmin(ixO^S),   1.0d0)
    cmax(ixO^S) = max(cmax(ixO^S), - 1.0d0)
    cmax(ixO^S) = min(cmax(ixO^S),   1.0d0)
    ! now take extremal value only for dt limit
    cmax(ixO^S) = max(dabs(cmax(ixO^S)),dabs(cmin(ixO^S)))

  end subroutine srhd_get_cmax

  subroutine srhd_get_a2max(w,x,ixI^L,ixO^L,a2max)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: a2max(ndim)
    double precision :: a2(ixI^S,ndim,nw)
    integer :: gxO^L,hxO^L,jxO^L,kxO^L,i

    a2=zero
    do i = 1,ndim
      !> 4th order
      hxO^L=ixO^L-kr(i,^D);
      gxO^L=hxO^L-kr(i,^D);
      jxO^L=ixO^L+kr(i,^D);
      kxO^L=jxO^L+kr(i,^D);
      a2(ixO^S,i,1:nwflux)=dabs(-w(kxO^S,1:nwflux)+16.d0*w(jxO^S,1:nwflux)&
        -30.d0*w(ixO^S,1:nwflux)+16.d0*w(hxO^S,1:nwflux)-w(gxO^S,1:nwflux))
      a2max(i)=maxval(a2(ixO^S,i,1:nwflux))/12.d0/dxlevel(i)**2
    end do

  end subroutine srhd_get_a2max

  !> local version for recycling code when computing cmax-cmin
  subroutine srhd_get_cmax_loc(ixI^L,ixO^L,vidim,csound2,v2,cmax,cmin)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in)             :: vidim(ixI^S)
    double precision, intent(in), dimension(ixO^S) :: csound2
    double precision, intent(in)             :: v2(ixI^S)
    double precision, intent(out)          :: cmax(ixI^S)
    double precision, intent(out)          :: cmin(ixI^S)

    double precision, dimension(ixI^S):: tmp1,tmp2

    tmp2(ixO^S)=vidim(ixO^S)**2.0d0
    tmp1(ixO^S)=1.0d0-v2(ixO^S)*csound2(ixO^S) &
                        -tmp2(ixO^S)*(1.0d0-csound2(ixO^S))
    tmp2(ixO^S)=dsqrt(csound2(ixO^S)*(one-v2(ixO^S))*tmp1(ixO^S))
    tmp1(ixO^S)=vidim(ixO^S)*(one-csound2(ixO^S))
    cmax(ixO^S)=(tmp1(ixO^S)+tmp2(ixO^S))/(one-v2(ixO^S)*csound2(ixO^S))
    ! Limit by speed of light
    cmax(ixO^S) = max(cmax(ixO^S), - 1.0d0)
    cmax(ixO^S) = min(cmax(ixO^S),   1.0d0)
    cmin(ixO^S)=(tmp1(ixO^S)-tmp2(ixO^S))/(one-v2(ixO^S)*csound2(ixO^S))
    ! Limit by speed of light
    cmin(ixO^S) = max(cmin(ixO^S), - 1.0d0)
    cmin(ixO^S) = min(cmin(ixO^S),   1.0d0)

  end subroutine srhd_get_cmax_loc

  !> Estimating bounds for the minimum and maximum signal velocities
  !> here we will not use Hspeed at all (one species only)
  subroutine srhd_get_cbounds(wLC,wRC,wLp,wRp,x,ixI^L,ixO^L,idim,Hspeed,cmax,cmin)
    use mod_global_parameters
    use mod_variables

    integer, intent(in)             :: ixI^L, ixO^L, idim
    ! conservative left and right status
    double precision, intent(in)    :: wLC(ixI^S, nw), wRC(ixI^S, nw)
    ! primitive left and right status
    double precision, intent(in)    :: wLp(ixI^S, nw), wRp(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    double precision, intent(inout) :: cmax(ixI^S,1:number_species)
    double precision, intent(inout), optional :: cmin(ixI^S,1:number_species)
    double precision, intent(in)    :: Hspeed(ixI^S,1:number_species)

    double precision :: wmean(ixI^S,nw)
    double precision, dimension(ixO^S) :: csound2,tmp1,tmp2,tmp3
    double precision, dimension(ixI^S) :: vidim,cmaxL,cmaxR,cminL,cminR,v2

    logical       :: flag(ixI^S,1:nw)

    select case(boundspeed)
    case(1) ! we do left-right first and take maximals
      !!call srhd_check_w(.true.,ixI^L, ixO^L, wLp, flag)
      !!call srhd_check_w_aux(ixI^L, ixO^L, wLp, flag)
      tmp1=wLp(ixO^S,rho_)
      tmp2=wLp(ixO^S,xi_)/wLp(ixO^S,lfac_)**2.0d0
      tmp3=wLp(ixO^S,p_)
      call srhd_get_csound2_prim_eos(ixO^L,tmp1,tmp2,tmp3,csound2)
      vidim(ixO^S) = wLp(ixO^S, mom(idim))/wLp(ixO^S, lfac_)
      v2(ixO^S) = 1.0d0-1.0d0/wLp(ixO^S,lfac_)**2
      call srhd_get_cmax_loc(ixI^L,ixO^L,vidim,csound2,v2,cmaxL,cminL)

      !!call srhd_check_w(.true.,ixI^L, ixO^L, wRp, flag)
      !!call srhd_check_w_aux(ixI^L, ixO^L, wRp, flag)
      tmp1=wRp(ixO^S,rho_)
      tmp2=wRp(ixO^S,xi_)/wRp(ixO^S,lfac_)**2.0d0
      tmp3=wRp(ixO^S,p_)
      call srhd_get_csound2_prim_eos(ixO^L,tmp1,tmp2,tmp3,csound2)
      vidim(ixO^S) = wRp(ixO^S, mom(idim))/wRp(ixO^S, lfac_)
      v2(ixO^S) = 1.0d0-1.0d0/wRp(ixO^S,lfac_)**2
      call srhd_get_cmax_loc(ixI^L,ixO^L,vidim,csound2,v2,cmaxR,cminR)

      if(present(cmin))then
        ! for HLL
        cmax(ixO^S,1)=max(cmaxL(ixO^S),cmaxR(ixO^S))
        cmin(ixO^S,1)=min(cminL(ixO^S),cminR(ixO^S))
      else
        ! for TVDLF
        cmaxL(ixO^S)=max(cmaxL(ixO^S),dabs(cminL(ixO^S)))
        cmaxR(ixO^S)=max(cmaxR(ixO^S),dabs(cminR(ixO^S)))
        cmax(ixO^S,1)=max(cmaxL(ixO^S),cmaxR(ixO^S))
      endif
    case(2) ! this is cmaxmean from conservatives
      ! here we do arithmetic mean of conservative vars
      wmean(ixO^S,1:nwflux)=0.5d0*(wLC(ixO^S,1:nwflux)+wRC(ixO^S,1:nwflux))
      ! get auxiliary variables
      call srhd_get_auxiliary(ixI^L,ixO^L,wmean,x)
      ! here tmp1 is rhoh
      tmp1=wmean(ixO^S,xi_)/wmean(ixO^S,lfac_)**2.0d0
      call srhd_get_csound2_rhoh(wmean,x,ixI^L,ixO^L,tmp1,csound2)
      vidim(ixO^S) = wmean(ixO^S, mom(idim))/wmean(ixO^S, xi_)
      v2(ixO^S)=1.0d0-1.0d0/wmean(ixO^S,lfac_)**2
      call srhd_get_cmax_loc(ixI^L,ixO^L,vidim,csound2,v2,cmaxL,cminL)
      if(present(cmin)) then
        cmax(ixO^S,1)=cmaxL(ixO^S)
        cmin(ixO^S,1)=cminL(ixO^S)
      else
        cmax(ixO^S,1)=max(cmaxL(ixO^S),dabs(cminL(ixO^S)))
      endif
    case(3) ! this is cmaxmean from primitives
      ! here we do arithmetic mean of primitive vars
      wmean(ixO^S,1:nwflux)=0.5d0*(wLp(ixO^S,1:nwflux)+wRp(ixO^S,1:nwflux))
      ! get auxiliary variables for wmean (primitive array)
      call srhd_get_auxiliary_prim(ixI^L,ixO^L,wmean)
      ! here tmp1 is rhoh
      tmp1=wmean(ixO^S,rho_)
      tmp2=wmean(ixO^S,xi_)/wmean(ixO^S,lfac_)**2.0d0
      tmp3=wmean(ixO^S,p_)
      call srhd_get_csound2_prim_eos(ixO^L,tmp1,tmp2,tmp3,csound2)
      vidim(ixO^S) = wmean(ixO^S, mom(idim))/wmean(ixO^S, lfac_)
      v2(ixO^S) = 1.0d0-1.0d0/wmean(ixO^S,lfac_)**2
      call srhd_get_cmax_loc(ixI^L,ixO^L,vidim,csound2,v2,cmaxL,cminL)
      if(present(cmin)) then
        cmax(ixO^S,1)=cmaxL(ixO^S)
        cmin(ixO^S,1)=cminL(ixO^S)
      else
        cmax(ixO^S,1)=max(cmaxL(ixO^S),dabs(cminL(ixO^S)))
      endif
    end select

  end subroutine srhd_get_cbounds

  !> Calculate fluxes within ixO^L.
  subroutine srhd_get_flux(wC,wP,x,ixI^L,ixO^L,idim,f)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L, idim
    ! conservative w
    double precision, intent(in) :: wC(ixI^S,nw)
    ! primitive w
    double precision, intent(in) :: wP(ixI^S,nw)
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision,intent(out) :: f(ixI^S,nwflux)

    double precision             :: pth(ixI^S)
    double precision             :: v(ixI^S,1:ndir)
    integer                      :: iw,idir

    pth(ixO^S)=wP(ixO^S,p_)
    do idir=1,ndir
      v(ixO^S,idir) = wP(ixO^S, mom(idir))/wP(ixO^S,lfac_)
    end do 

    ! Get flux of density d, namely D*v 
    f(ixO^S,d_)=v(ixO^S,idim)*wC(ixO^S,rho_)

    ! Get flux of tracer
    do iw=1,srhd_n_tracer
      f(ixO^S,tracer(iw))=v(ixO^S,idim)*wC(ixO^S,tracer(iw))
    end do

    ! Get flux of momentum
    ! f_i[m_k]=v_i*m_k [+pth if i==k]
    do idir=1,ndir
      f(ixO^S,mom(idir))= v(ixO^S,idim)*wC(ixO^S,mom(idir))
    end do 
    f(ixO^S,mom(idim))=pth(ixO^S)+f(ixO^S,mom(idim))

    ! Get flux of energy
    ! f_i[e]=v_i*e+v_i*pth
    f(ixO^S,e_)=v(ixO^S,idim)*(wC(ixO^S,e_) + pth(ixO^S))

  end subroutine srhd_get_flux

  !> Add geometrical source terms to w
  subroutine srhd_add_source_geom(qdt, dtfactor, ixI^L, ixO^L, wCT, wprim, w, x)
    use mod_global_parameters
    use mod_usr_methods, only: usr_set_surface
    use mod_geometry
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, dtfactor, x(ixI^S, 1:ndim)
    double precision, intent(inout) :: wCT(ixI^S, 1:nw), wprim(ixI^S, 1:nw), w(ixI^S, 1:nw)

    double precision :: pth(ixI^S), source(ixI^S), v(ixI^S,1:ndir)
    integer                         :: idir, h1x^L{^NOONED, h2x^L}
    integer :: mr_,mphi_ ! Polar var. names
    double precision :: exp_factor(ixI^S), del_exp_factor(ixI^S), exp_factor_primitive(ixI^S)

    select case (coordinate)

    case(Cartesian_expansion)
      !the user provides the functions of exp_factor and del_exp_factor
      if(associated(usr_set_surface)) &
        call usr_set_surface(ixI^L,x,block%dx,exp_factor,del_exp_factor,exp_factor_primitive)
      ! get auxiliary variables lfac and xi from conserved set
      call srhd_get_auxiliary(ixI^L,ixO^L,wCT,x)
      call srhd_get_pthermal(wCT, x, ixI^L, ixO^L, source)
      source(ixO^S) = source(ixO^S)*del_exp_factor(ixO^S)/exp_factor(ixO^S)
      w(ixO^S,mom(1)) = w(ixO^S,mom(1)) + qdt*source(ixO^S)

    case (cylindrical)
          mr_   = mom(r_)
          ! get auxiliary variables lfac and xi from conserved set
          call srhd_get_auxiliary(ixI^L,ixO^L,wCT,x)
          call srhd_get_pthermal(wCT, x, ixI^L, ixO^L, source)
          if (phi_ > 0) then
             mphi_ = mom(phi_)
             source(ixO^S) = source(ixO^S) + wCT(ixO^S, mphi_)*wprim(ixO^S,mom(phi_))
             w(ixO^S, mr_) = w(ixO^S, mr_) + qdt * source(ixO^S) / x(ixO^S, r_)
             source(ixO^S) = -wCT(ixO^S, mphi_) * wprim(ixO^S,mom(r_))
             w(ixO^S, mphi_) = w(ixO^S, mphi_) + qdt * source(ixO^S) / x(ixO^S, r_)
          else
             w(ixO^S, mr_) = w(ixO^S, mr_) + qdt * source(ixO^S) / x(ixO^S, r_)
          end if

    case (spherical)
       mr_   = mom(r_)
       h1x^L=ixO^L-kr(1,^D); {^NOONED h2x^L=ixO^L-kr(2,^D);}
       ! s[mr]=((stheta*vtheta+sphi*vphi)+2*p)/r
       ! get auxiliary variables lfac and xi from conserved set
       call srhd_get_auxiliary(ixI^L,ixO^L,wCT,x)
       call srhd_get_pthermal(wCT, x, ixI^L, ixO^L, pth)
       source(ixO^S) = pth(ixO^S) * x(ixO^S, 1) &
            *(block%surfaceC(ixO^S, 1) - block%surfaceC(h1x^S, 1)) &
            /block%dvolume(ixO^S)
       if (ndir > 1) then
         do idir = 2, ndir
           source(ixO^S) = source(ixO^S) + wCT(ixO^S, mom(idir))*wprim(ixO^S,mom(idir))
         end do
       end if
       w(ixO^S, mr_) = w(ixO^S, mr_) + qdt * source(ixO^S) / x(ixO^S, 1)

       {^NOONED
       ! s[mtheta]=-(stheta*vr)/r+cot(theta)*(sphi*vphi+p)/r
       source(ixO^S) = pth(ixO^S) * x(ixO^S, 1) &
            * (block%surfaceC(ixO^S, 2) - block%surfaceC(h2x^S, 2)) &
            / block%dvolume(ixO^S)
       if (ndir == 3) then
          source(ixO^S) = source(ixO^S) + (wCT(ixO^S, mom(3))*v(ixO^S,ndir)) / dtan(x(ixO^S, 2))
       end if
       source(ixO^S) = source(ixO^S) - wCT(ixO^S, mom(2)) * wprim(ixO^S, mom(1))
       w(ixO^S, mom(2)) = w(ixO^S, mom(2)) + qdt * source(ixO^S) / x(ixO^S, 1)

       if (ndir == 3) then
         ! s[mphi]=-(sphi*vr)/r-cot(theta)*(sphi*vtheta)/r
         source(ixO^S) = -(wCT(ixO^S, mom(3)) * wprim(ixO^S, mom(1))) &
                        - (wCT(ixO^S, mom(3)) * wprim(ixO^S, mom(2))) / dtan(x(ixO^S, 2))
         w(ixO^S, mom(3)) = w(ixO^S, mom(3)) + qdt * source(ixO^S) / x(ixO^S, 1)
       end if
       }
    end select

  end subroutine srhd_add_source_geom

  !> handles bootstrapping
  subroutine srhd_handle_small_values(primitive, w, x, ixI^L, ixO^L, subname)
    use mod_global_parameters
    use mod_small_values
    logical, intent(in)             :: primitive
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    character(len=*), intent(in)    :: subname

    integer :: n,idir
    logical :: flag(ixI^S,1:nw),flagall(ixI^S)

    call srhd_check_w(primitive, ixI^L, ixO^L, w, flag)

    if (any(flag)) then
      select case (small_values_method)
      case ("replace")
        ! any faulty cell is replaced by physical lower limit
        flagall(ixO^S)=(flag(ixO^S,rho_).or.flag(ixO^S,e_)) 

        where(flagall(ixO^S))
           ! D or rho: no difference primitive-conservative
           w(ixO^S,rho_) = small_density
           !w(ixO^S,lfac_)= 1.0d0
           !w(ixO^S,xi_)  = small_xi
        endwhere
        !do idir = 1, ndir
        !   where(flagall(ixO^S)) w(ixO^S, mom(idir)) = 0.0d0
        !end do
        if(primitive) then
            where(flagall(ixO^S)) w(ixO^S, p_) = small_pressure
        else
            where(flagall(ixO^S)) w(ixO^S, e_) = small_e 
        endif

      case ("average")
        ! note: in small_values_average we use 
        ! small_values_fix_iw(1:nw) and small_values_daverage       
        ! when fails, may use small_pressure/small_density
        if(primitive)then
           ! averaging for all primitive fields (p, lfac*v, tau))
           call small_values_average(ixI^L, ixO^L, w, x, flag)
           ! update the auxiliaries from primitives
           call srhd_get_auxiliary_prim(ixI^L,ixO^L,w)
        else
           ! do averaging of density d
           call small_values_average(ixI^L, ixO^L, w, x, flag, d_)
           ! do averaging of energy tau
           call small_values_average(ixI^L, ixO^L, w, x, flag, e_)
           ! and now hope for the best....
        endif
      case default
        if(.not.primitive) then
          ! note that we throw error here, which assumes w is primitive
          write(*,*) "handle_small_values default: note reporting conservatives!"
        end if
        call small_values_error(w, x, ixI^L, ixO^L, flag, subname)
      end select
    end if

  end subroutine srhd_handle_small_values

  !> calculate effective gamma
  subroutine srhd_get_Geff_eos(w,ixI^L,ixO^L,varconserve,Geff)
    use mod_global_parameters, only: nw
  !================== IMPORTANT ==================!
  !This subroutine is used with conserved variables in w when varconserve=T
  !This subroutine is used with primitive variables in w when varconserve=F
  !   both cases assume updated auxiliary variables xi_ en lfac_
  !===============================================!
    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: w(ixI^S, 1:nw)
    logical, intent(in)                :: varconserve
    double precision, intent(out)      :: Geff(ixI^S)

    double precision, dimension(ixO^S) :: pth,rho,E_th,E

    if (srhd_eos) then
      if (varconserve) then
        pth(ixO^S)=w(ixO^S,xi_)-w(ixO^S,e_)-w(ixO^S,d_)
        rho(ixO^S)=w(ixO^S,d_)/w(ixO^S,lfac_)
        E_th = pth*inv_gamma_1
        E    = E_th+dsqrt(E_th**2+rho**2)
        Geff(ixO^S) = srhd_gamma-half*gamma_1 *          &
                             (one-(rho(ixO^S)/E(ixO^S))**2)
      else
        ! primitives available
        E_th = w(ixO^S,p_)*inv_gamma_1
        E    = E_th+dsqrt(E_th**2+w(ixO^S,rho_)**2)
        Geff(ixO^S) = srhd_gamma-half*gamma_1 *          &
                             (one-(w(ixO^S,rho_)/E(ixO^S))**2)
      end if
    else
      Geff(ixO^S) = srhd_gamma
    endif

  end subroutine srhd_get_Geff_eos

  !> Compute the small value limits
  subroutine srhd_get_smallvalues_eos
    use mod_global_parameters, only: small_pressure, small_density
    implicit none
    ! local small values
    double precision :: LsmallE,Lsmallp,Lsmallrho

    ! the maximal allowed Lorentz factor
    lfacmax=one/dsqrt(one-(one-dmaxvel)**2)
    minrho=small_density
    minp=small_pressure
    if(small_density*small_pressure<=0.0d0)then
       call mpistop("must set finite values small-density/pressure for small value treatments")
    endif
    if(srhd_eos)then
       Lsmallp=(one+10.d0*small_pressure)*small_pressure
       Lsmallrho=(one+10.d0*small_density)*small_density
       !!Lsmallp=small_pressure
       !!Lsmallrho=small_density
       LsmallE=Lsmallp*inv_gamma_1+&
                dsqrt((Lsmallp*inv_gamma_1)**2+Lsmallrho**2)
       small_xi=half*((srhd_gamma+one)*LsmallE-&
                      gamma_1*Lsmallrho*(Lsmallrho/LsmallE))
       small_e=small_xi-Lsmallp-Lsmallrho
    else
       small_xi=small_density+gamma_to_gamma_1*small_pressure
       small_e =small_pressure*inv_gamma_1
    endif
    smallxi=small_xi
    smalltau=small_e

  end subroutine srhd_get_smallvalues_eos

  !> Compute the enthalpy rho*h from rho and pressure p
  subroutine srhd_get_enthalpy_eos(ixO^L,rho,p,rhoh)
    use mod_global_parameters
    integer, intent(in)                :: ixO^L
    double precision, intent(in)       :: rho(ixO^S),p(ixO^S)
    double precision, intent(out)      :: rhoh(ixO^S)

    double precision, dimension(ixO^S) :: E_th,E
    integer :: ix^D

    if(srhd_eos) then
     E_th = p*inv_gamma_1
     E    = E_th+dsqrt(E_th**2+rho**2)
     ! writing rho/E on purpose, for numerics 
     rhoh = half*((srhd_gamma+one)*E &
                   - gamma_1*rho*(rho/E))
    else
     rhoh = rho+gamma_to_gamma_1*p
    end if

    if (check_small_values) then
      {do ix^DB= ixO^LIM^DB\}
         if(rhoh(ix^D)<small_xi) then
           write(*,*) "local pressure and density",p(ix^D),rho(ix^D)
           write(*,*) "Error: small value of enthalpy rho*h=",rhoh(ix^D),&
                " encountered when call srhd_get_enthalpy_eos"
           call mpistop('enthalpy below small_xi: stop (may need to turn on fixes)')
         end if
      {enddo^D&\}
    end if

    if (fix_small_values) then
      {do ix^DB= ixO^LIM^DB\}
         if(rhoh(ix^D)<small_xi) then
            rhoh(ix^D)=small_xi
         endif
      {enddo^D&\}
    endif

  end subroutine srhd_get_enthalpy_eos

  !> Calculate thermal pressure p from density rho and enthalpy rho*h 
  !> will provide p (and E if srhd_eos)
  subroutine srhd_get_pressure_eos(ixI^L,ixO^L,rho,rhoh,p,E)
    use mod_global_parameters
    integer, intent(in)            :: ixI^L, ixO^L
    double precision, intent(in)   :: rho(ixO^S),rhoh(ixO^S)
    double precision, intent(out)  :: p(ixI^S)
    double precision, intent(out)  :: E(ixO^S)
    integer :: ix^D

    if(srhd_eos) then
     E = (rhoh+dsqrt(rhoh**2+(srhd_gamma**2-one)*rho**2)) &
         /(srhd_gamma+one)
     p(ixO^S) = half*gamma_1* (E-rho*(rho/E))
    else 
     p(ixO^S) = (rhoh-rho)/gamma_to_gamma_1
    end if 

    if (check_small_values) then
      {do ix^DB= ixO^LIM^DB\}
         if(p(ix^D)<small_pressure) then
           write(*,*) "local enthalpy rho*h and density rho",rhoh(ix^D),rho(ix^D)
           if(srhd_eos) write(*,*) 'E, rho^2/E, difference', &
                       E(ix^D),rho(ix^D)**2/E(ix^D),E(ix^D)-rho(ix^D)**2/E(ix^D)
           write(*,*) "Error: small value of gas pressure",p(ix^D),&
                " encountered when call srhd_get_pressure_eos"
           call mpistop('pressure below small_pressure: stop (may need to turn on fixes)')
         end if
      {enddo^D&\}
    end if

    if (fix_small_values) then
      {do ix^DB= ixO^LIM^DB\}
         if(p(ix^D)<small_pressure) then
            p(ix^D)=small_pressure
            if(srhd_eos)E(ix^D)=max(small_e,E(ix^D))
         endif
      {enddo^D&\}
    endif

  end subroutine srhd_get_pressure_eos

  !> Calculate the square of the thermal sound speed csound2 within ixO^L.
  !> available rho - rho*h 
  subroutine srhd_get_csound2_eos(ixI^L,ixO^L,rho,rhoh,csound2)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: rho(ixO^S),rhoh(ixO^S)
    double precision, intent(out)   :: csound2(ixO^S)

    double precision                :: p(ixI^S)
    double precision                :: E(ixO^S)
    integer :: ix^D

    call srhd_get_pressure_eos(ixI^L,ixO^L,rho,rhoh,p,E)
    if(srhd_eos) then
       csound2(ixO^S)=(p(ixO^S)*((srhd_gamma+one)&
                          +gamma_1*(rho(ixO^S)/E(ixO^S))**2))&
                      /(2.0d0*rhoh(ixO^S))
    else
       csound2(ixO^S)=srhd_gamma*p(ixO^S)/rhoh(ixO^S)
    end if

    if (check_small_values) then
      {do ix^DB= ixO^LIM^DB\}
         if(csound2(ix^D)>=1.0d0.or.csound2(ix^D)<=0.0d0) then
           write(*,*) "sound speed error with p - rho - rhoh",p(ix^D),rhoh(ix^D),rho(ix^D)
           if(srhd_eos) write(*,*) 'and E', E(ix^D)
           write(*,*) "Error: value of csound2",csound2(ix^D),&
                " encountered when call srhd_get_csound2_eos"
           call mpistop('sound speed stop (may need to turn on fixes)')
         end if
      {enddo^D&\}
    end if

    if (fix_small_values) then
      {do ix^DB= ixO^LIM^DB\}
         if(csound2(ix^D)>=1.0d0) then
            csound2(ix^D)=1.0d0-1.0d0/lfacmax**2
         endif
         if(csound2(ix^D)<=0.0d0) then
            csound2(ix^D)=srhd_gamma*small_pressure/small_xi
         endif
      {enddo^D&\}
    endif

  end subroutine srhd_get_csound2_eos

  !> Calculate the square of the thermal sound speed csound2 within ixO^L.
  !> available rho - rho*h - p
  subroutine srhd_get_csound2_prim_eos(ixO^L,rho,rhoh,p,csound2)
    use mod_global_parameters
    integer, intent(in)             :: ixO^L
    double precision, intent(in)    :: rho(ixO^S),rhoh(ixO^S),p(ixO^S)
    double precision, intent(out)   :: csound2(ixO^S)

    double precision                :: E(ixO^S)
    integer :: ix^D

    if(srhd_eos) then
       E = (rhoh+dsqrt(rhoh**2+(srhd_gamma**2-one)&
            *rho**2))/(srhd_gamma+one)
       csound2(ixO^S)=(p*((srhd_gamma+one)&
                          +gamma_1*(rho/E)**2))&
                      /(2.0d0*rhoh)
    else
       csound2(ixO^S)=srhd_gamma*p/rhoh
    end if

    if (check_small_values) then
      {do ix^DB= ixO^LIM^DB\}
         if(csound2(ix^D)>=1.0d0.or.csound2(ix^D)<=0.0d0) then
           write(*,*) "sound speed error with p - rho - rhoh",p(ix^D),rhoh(ix^D),rho(ix^D)
           if(srhd_eos) write(*,*) 'and E', E(ix^D)
           write(*,*) "Error: value of csound2",csound2(ix^D),&
                " encountered when call srhd_get_csound2_prim_eos"
           call mpistop('sound speed stop (may need to turn on fixes)')
         end if
      {enddo^D&\}
    end if

    if (fix_small_values) then
      {do ix^DB= ixO^LIM^DB\}
         if(csound2(ix^D)>=1.0d0) then
            csound2(ix^D)=1.0d0-1.0d0/lfacmax**2
         endif
         if(csound2(ix^D)<=0.0d0) then
            csound2(ix^D)=srhd_gamma*small_pressure/small_xi
         endif
      {enddo^D&\}
    endif

  end subroutine srhd_get_csound2_prim_eos

  !> con2prim: (D,S**2,tau) --> compute auxiliaries lfac and xi
  subroutine con2prim_eos(lfac,xi,myd,myssqr,mytau,ierror)
    use mod_con2prim_vars

    double precision, intent(in)    :: myd, myssqr, mytau
    double precision, intent(inout) :: lfac, xi
    integer, intent(inout)          :: ierror

    ! .. local ..
    double precision:: f,df,lfacl
    !------------------------------------------------------------------

    ! Save the input-state in mod_con2prim_vars 
    d = myd; ssqr = myssqr; tau = mytau; 

    ierror=0

    ! Check if guess is close enough: gives f,df,lfacl
    if(xi>smallxi)then
    call funcd_eos(xi,f,df,lfacl,d,ssqr,tau,ierror)
    if (ierror == 0 .and. dabs(f/df)<absaccnr) then
       xi   = xi - f/df
       lfac = lfacl
       return
    else
       ierror = 0
    end if
    else
      write(*,*)'entering con2prim_eos with xi=',xi
    end if

    ! ierror=1 : must have D>=minrho, tau>=smalltau
    !if(d<minrho .or. tau<smalltau)then
    !   ierror=1
    !   return
    !endif

    call con2primHydro_eos(lfac,xi,d,ssqr,tau,ierror)

  end subroutine con2prim_eos

  subroutine funcd_eos(xi,f,df,mylfac,d,ssqr,tau,ierror)
    double precision, intent(in)  :: xi,d,ssqr,tau
    double precision, intent(out) :: f,df,mylfac
    integer, intent(inout)        :: ierror

    ! .. local ..
    double precision  :: dlfac
    double precision  :: vsqr,p,dpdxi
    !-----------------------------------------------------------------

    vsqr = ssqr/xi**2

    if (vsqr<one) then
       mylfac = one/dsqrt(one-vsqr)
       dlfac = -mylfac**3*ssqr/(xi**3)
       !===== Pressure, calculate using EOS =====!
       call FuncPressure_eos(xi,mylfac,d,dlfac,p,dpdxi)
       !=========================================!
       f  = xi-tau-d-p
       df = one-dpdxi
    else
       ! print *,'Erroneous input to funcd since vsqr=',vsqr,' >=1'
       ! print *,'input values d, ssqr, tau:',d,ssqr,tau
       ierror =6
       return
    end if

  end subroutine funcd_eos

  !> SRHD iteration solves for p via NR, and then gives xi as output
  subroutine con2primHydro_eos(lfac,xi,d,sqrs,tau,ierror)
    double precision, intent(out) :: xi,lfac
    double precision, intent(in)  :: d,sqrs,tau
    integer,intent(inout)         :: ierror

    ! .. local ..
    integer          :: ni,niiter,nit,n2it,ni3
    double precision :: pcurrent,pnew
    double precision :: er,er1,ff,df,dp,v2
    double precision :: pmin,lfac2inv,pLabs,pRabs,pprev
    double precision :: s2overcubeG2rh
    double precision :: xicurrent,h,dhdp
    double precision :: oldff1,oldff2,Nff
    double precision :: pleft,pright
    !---------------------------------------------------------------------

    ierror=0
    ! ierror=0 : ok
    !  we already checked D>=minrho, tau>=smalltau (ierror=1)
    !
    ! ierror<>0
    !
    ! ierror=2 : maxitnr reached without convergence
    ! ierror=3 : final pressure value < smallp or xi<smallxi during iteration
    ! ierror=4 : final v^2=1 hence problem as lfac=1/0
    ! ierror=5 : nonmonotonic function f (as df=0)

    ! left and right brackets for p-range
    pmin=dsqrt(sqrs)/(one-dmaxvel)-tau-d
    pLabs=max(minp,pmin)
    pRabs=1.0d99
    ! start value from input
    pcurrent=pLabs

    er1=one
    pprev=pcurrent

    ! Fudge Parameters
    oldff1=1.0d7  ! High number
    oldff2=1.0d9  ! High number bigger then oldff1
    n2it = 0
    nit  = 0

    LoopNR:  do ni=1,maxitnr
       nit = nit + 1
       !=== Relax NR iteration accuracy=======!
       if(nit>maxitnr/4)then
          ! mix pressure value for convergence
          pcurrent=half*(pcurrent+pprev)
          ! relax accuracy requirement
          er1=10.0d0*er1
          nit = nit - maxitnr/10
       endif
       !=======================================!

       niiter=ni
       xicurrent=tau+d+pcurrent

       if(xicurrent<smallxi) then
          !        print *,'stop: too small xi iterate:',xicurrent
          !        print *,'for pressure iterate p',pcurrent
          !        print *,'pressure bracket pLabs pRabs',pLabs,pRabs
          !        print *,'iteration number:',ni
          !        print *,'values for d,s,tau,s2:',d,sqrs,tau,sqrs
          ierror=3
          return
       endif

       v2=sqrs/xicurrent**2
       lfac2inv=one - v2
       if(lfac2inv>zero) then
          lfac=one/dsqrt(lfac2inv)
       else
          !        print *,'stop: negative or zero factor 1-v2:',lfac2inv
          !        print *,'for pressure iterate p',pcurrent
          !        print *,'absolute pressure bracket pLabs pRabs',pLabs,pRabs
          !        print *,'iteration number:',ni
          !        print *,'values for d,s,tau,s2:',d,sqrs,tau,sqrs
          !        print *,'values for v2,xi:',v2,xicurrent
          ierror=4
          return
       endif

       s2overcubeG2rh=sqrs/(xicurrent**3)
       !== calculation done using the EOS ==!
       call FuncEnthalpy_eos(pcurrent,lfac2inv,d,sqrs,xicurrent,&
                          s2overcubeG2rh,h,dhdp)
       !=======================================!
       ff=-xicurrent*lfac2inv + h
       df=- two*sqrs/xicurrent**2  + dhdp - lfac2inv

       if (ff*df==zero) then
          if (ff==zero) then
             exit ! zero found
          else
             !     print *,'stop: df becomes zero, non-monotonic f(p)!'
             ierror=5
             return
          endif
       else
          pnew=pcurrent-ff/df
          if (ff*df>zero) then
             ! pressure iterate has decreased
             ! restrict to left
             pnew=max(pnew,pLabs)
          else  ! ff*df<0
             ! pressure iterate has increased
             ! restrict to right
             pnew=min(pnew,pRabs)
          endif
       endif

       !===============================================!
       dp=pcurrent-pnew
       er=two*dabs(dp)/(pnew+pcurrent)
       if(((er<tolernr*er1).or.(dabs(dp)<absaccnr))) exit LoopNR
       !===============================================!

       ! For very small values of pressure, NR algorithm is not efficient to
       ! find root, use Euler algorithm to find precise value of pressure
       if((dabs(oldff2-ff) < 1.0d-8 .or. niiter >= maxitnr-maxitnr/20).and.&
            ff * oldff1 < zero    .and.  dabs(ff)>absaccnr)then

          n2it=n2it+1
          if(n2it<=3) pcurrent=half*(pnew+pcurrent)
          if(n2it>3)then
             pright =pcurrent
             pleft=pprev
             pcurrent=half*(pleft+pright)
             Dicho:  do ni3=1,maxitnr
                !===================!
                xicurrent=tau+d+pcurrent
                v2=sqrs/xicurrent**2
                lfac2inv=one - v2

                if(lfac2inv>zero)then
                   lfac=one/dsqrt(lfac2inv)
                else
                   ierror=4
                   return
                endif
                !===================!

                !== calculation done using the EOS ==!
                call Bisection_Enthalpy_eos(pnew,lfac2inv,d,xicurrent,h)
                Nff=-xicurrent*lfac2inv + h
                !=======================================!
                !==== Iterate ====!
                if(ff * Nff < zero)then
                   pleft=pcurrent
                else
                   pright=pcurrent
                endif

                pcurrent=half*(pleft+pright)
                !==================!

                !=== The iteration converged ===!
                if(2.0d0*dabs(pleft-pright)/(pleft+pright)< absaccnr &
                     .or. dabs(ff)<absaccnr)then
                   pnew=pcurrent
                   exit LoopNR
                endif
                !==============================!

                !==============================!

                !=== conserve the last value of Nff ===!
                ff=Nff
                !======================================!
             enddo    Dicho
          endif

       else
          !====== There is no problems, continue the NR iteration ======!
          pprev=pcurrent
          pcurrent=pnew
          !=============================================================!
       endif


       !=== keep the values of the 2 last ff ===!
       oldff2=oldff1
       oldff1=ff
       !========================================!
    enddo LoopNR

    if(niiter==maxitnr)then
       ierror=2
       return
    endif

    if(pcurrent<minp) then
       ierror=3
       return
    endif

    !------------------------------!
    xi=tau+d+pcurrent
    v2=sqrs/xicurrent**2
    lfac2inv=one - v2
    if(lfac2inv>zero) then
       lfac=one/dsqrt(lfac2inv)
    else
       ierror=4
       return
    endif

  end subroutine con2primHydro_eos

  !> pointwise evaluations used in con2prim
  !> compute pointwise value for pressure p and dpdxi
  subroutine FuncPressure_eos(xicurrent,lfac,d,dlfacdxi,p,dpdxi)

    double precision, intent(in)         :: xicurrent,lfac,d,dlfacdxi
    double precision, intent(out)        :: p,dpdxi
    ! .. local ..
    double precision                     :: rho,h,E,dhdxi,rhotoE
    double precision                     :: dpdchi,dEdxi

    ! rhoh here called h
    h=xicurrent/(lfac**2)
    rho=d/lfac
    E = (h+dsqrt(h**2+(srhd_gamma**2-one)*rho**2)) &
              /(srhd_gamma+one)
    ! output pressure
    rhotoE = rho/E
    p = half*gamma_1*(E-rho*rhotoE)

    dhdxi = one/(lfac**2)-2.0d0*xicurrent/(lfac**2)*dlfacdxi/lfac

    dEdxi=(dhdxi+(h*dhdxi-(srhd_gamma**2-one)*rho**2*dlfacdxi/lfac)&
        /dsqrt(h**2+(srhd_gamma**2-one)*rho**2))&
        /(srhd_gamma+one)

    ! output pressure derivative to xi
    dpdxi=half*gamma_1*(2.0d0*rho*rhotoE*dlfacdxi/lfac+&
          (one+rhotoE**2)*dEdxi)

  end subroutine FuncPressure_eos

  !> pointwise evaluations used in con2prim
  !> returns enthalpy rho*h (h) and derivative d(rho*h)/dp (dhdp)
  subroutine FuncEnthalpy_eos(pcurrent,lfac2inv,d,sqrs,xicurrent,dv2d2p,h,dhdp)

    double precision, intent(in) :: pcurrent,lfac2inv,d,sqrs,xicurrent,dv2d2p
    double precision, intent(out):: h,dhdp

    ! local
    double precision:: rho,E_th,E,dE_thdp,dEdp

    rho=d*dsqrt(lfac2inv)
    E_th = pcurrent*inv_gamma_1
    E = (E_th + dsqrt(E_th**2+rho**2))
    !== Enthalpy ==!
    h = half*((srhd_gamma+one)*E-gamma_1*rho*(rho/E))
    !=== Derivative of thermal energy ===!
    dE_thdp = one*inv_gamma_1
    !=== Derivative of internal energy ===!
    dEdp = dE_thdp * (one+E_th/dsqrt(E_th**2+rho**2))&
              +  d**2*dv2d2p/dsqrt(E_th**2+rho**2)
    !====== Derivative of Enthalpy ======!
    dhdp = half*((srhd_gamma+one)*dEdp + &
              gamma_1*(rho*(rho/E))*(-2.0d0*dv2d2p/lfac2inv+dEdp/E))
  end subroutine FuncEnthalpy_eos

  !> pointwise evaluations used in con2prim
  !> returns enthalpy rho*h (h)
  subroutine Bisection_Enthalpy_eos(pcurrent,lfac2inv,d,xicurrent,h)

    double precision, intent(in) :: pcurrent,lfac2inv,d,xicurrent
    double precision, intent(out):: h

    ! local
    double precision:: rho,E_th,E

    rho=d*dsqrt(lfac2inv)
    E_th = pcurrent*inv_gamma_1
    E = (E_th + dsqrt(E_th**2+rho**2))
    !== Enthalpy ==!
    h = half*((srhd_gamma+one)*E-gamma_1*rho*(rho/E))

    return
  end subroutine Bisection_Enthalpy_eos

  !> con2prim: (D,S**2,tau) --> compute auxiliaries lfac and xi
  subroutine con2prim(lfac,xi,myd,myssqr,mytau,ierror)
    use mod_con2prim_vars

    double precision, intent(in)    :: myd, myssqr, mytau
    double precision, intent(inout) :: lfac, xi
    integer, intent(inout)          :: ierror

    ! .. local ..
    double precision:: f,df,lfacl
    !------------------------------------------------------------------

    ! Save the input-state in mod_con2prim_vars 
    d = myd; ssqr = myssqr; tau = mytau; 

    ierror=0

    ! Check if guess is close enough: gives f,df,lfacl
    if(xi>smallxi)then
    call funcd(xi,f,df,lfacl,d,ssqr,tau,ierror)
    if (ierror == 0 .and. dabs(f/df)<absaccnr) then
       xi   = xi - f/df
       lfac = lfacl
       return
    else
       ierror = 0
    end if
    else
       write(*,*) 'entering con2prim with xi=',xi
    end if

    ! ierror=1 : must have D>=minrho, tau>=smalltau
    !if(d<minrho .or. tau<smalltau)then
    !   ierror=1
    !   return
    !endif

    call con2primHydro(lfac,xi,d,ssqr,tau,ierror)

  end subroutine con2prim

  subroutine funcd(xi,f,df,mylfac,d,ssqr,tau,ierror)
    double precision, intent(in)  :: xi,d,ssqr,tau
    double precision, intent(out) :: f,df,mylfac
    integer, intent(inout)        :: ierror

    ! .. local ..
    double precision  :: dlfac
    double precision  :: vsqr,p,dpdxi
    !-----------------------------------------------------------------

    vsqr = ssqr/xi**2

    if (vsqr<one) then
       mylfac = one/dsqrt(one-vsqr)
       dlfac = -mylfac**3*ssqr/(xi**3)
       !===== Pressure, calculate using EOS =====!
       call FuncPressure(xi,mylfac,d,dlfac,p,dpdxi)
       !=========================================!
       f  = xi-tau-d-p
       df = one-dpdxi
    else
       ! print *,'Erroneous input to funcd since vsqr=',vsqr,' >=1'
       ! print *,'input values d, ssqr, tau:',d,ssqr,tau
       ierror =6
       return
    end if

  end subroutine funcd

  !> SRHD iteration solves for p via NR, and then gives xi as output
  subroutine con2primHydro(lfac,xi,d,sqrs,tau,ierror)
    double precision, intent(out) :: xi,lfac
    double precision, intent(in)  :: d,sqrs,tau
    integer,intent(inout)         :: ierror

    ! .. local ..
    integer          :: ni,niiter,nit,n2it,ni3
    double precision :: pcurrent,pnew
    double precision :: er,er1,ff,df,dp,v2
    double precision :: pmin,lfac2inv,pLabs,pRabs,pprev
    double precision :: s2overcubeG2rh
    double precision :: xicurrent,h,dhdp
    double precision :: oldff1,oldff2,Nff
    double precision :: pleft,pright
    !---------------------------------------------------------------------

    ierror=0
    ! ierror=0 : ok
    !  we already checked D>=minrho, tau>=smalltau (ierror=1)
    !
    ! ierror<>0
    !
    ! ierror=2 : maxitnr reached without convergence
    ! ierror=3 : final pressure value < smallp or xi<smallxi during iteration
    ! ierror=4 : final v^2=1 hence problem as lfac=1/0
    ! ierror=5 : nonmonotonic function f (as df=0)

    ! left and right brackets for p-range
    pmin=dsqrt(sqrs)/(one-dmaxvel)-tau-d
    pLabs=max(minp,pmin)
    pRabs=1.0d99
    ! start value from input
    pcurrent=pLabs

    er1=one
    pprev=pcurrent

    ! Fudge Parameters
    oldff1=1.0d7  ! High number
    oldff2=1.0d9  ! High number bigger then oldff1
    n2it = 0
    nit  = 0

    LoopNR:  do ni=1,maxitnr
       nit = nit + 1
       !=== Relax NR iteration accuracy=======!
       if(nit>maxitnr/4)then
          ! mix pressure value for convergence
          pcurrent=half*(pcurrent+pprev)
          ! relax accuracy requirement
          er1=10.0d0*er1
          nit = nit - maxitnr/10
       endif
       !=======================================!

       niiter=ni
       xicurrent=tau+d+pcurrent

       if(xicurrent<smallxi) then
          !        print *,'stop: too small xi iterate:',xicurrent
          !        print *,'for pressure iterate p',pcurrent
          !        print *,'pressure bracket pLabs pRabs',pLabs,pRabs
          !        print *,'iteration number:',ni
          !        print *,'values for d,s,tau,s2:',d,sqrs,tau,sqrs
          ierror=3
          return
       endif

       v2=sqrs/xicurrent**2
       lfac2inv=one - v2
       if(lfac2inv>zero) then
          lfac=one/dsqrt(lfac2inv)
       else
          !        print *,'stop: negative or zero factor 1-v2:',lfac2inv
          !        print *,'for pressure iterate p',pcurrent
          !        print *,'absolute pressure bracket pLabs pRabs',pLabs,pRabs
          !        print *,'iteration number:',ni
          !        print *,'values for d,s,tau,s2:',d,sqrs,tau,sqrs
          !        print *,'values for v2,xi:',v2,xicurrent
          ierror=4
          return
       endif

       s2overcubeG2rh=sqrs/(xicurrent**3)
       !== calculation done using the EOS ==!
       call FuncEnthalpy(pcurrent,lfac2inv,d,sqrs,xicurrent,&
                          s2overcubeG2rh,h,dhdp)
       !=======================================!
       ff=-xicurrent*lfac2inv + h
       df=- two*sqrs/xicurrent**2  + dhdp - lfac2inv

       if (ff*df==zero) then
          if (ff==zero) then
             exit ! zero found
          else
             !     print *,'stop: df becomes zero, non-monotonic f(p)!'
             ierror=5
             return
          endif
       else
          pnew=pcurrent-ff/df
          if (ff*df>zero) then
             ! pressure iterate has decreased
             ! restrict to left
             pnew=max(pnew,pLabs)
          else  ! ff*df<0
             ! pressure iterate has increased
             ! restrict to right
             pnew=min(pnew,pRabs)
          endif
       endif

       !===============================================!
       dp=pcurrent-pnew
       er=two*dabs(dp)/(pnew+pcurrent)
       if(((er<tolernr*er1).or.(dabs(dp)<absaccnr))) exit LoopNR
       !===============================================!

       ! For very small values of pressure, NR algorithm is not efficient to
       ! find root, use Euler algorithm to find precise value of pressure
       if((dabs(oldff2-ff) < 1.0d-8 .or. niiter >= maxitnr-maxitnr/20).and.&
            ff * oldff1 < zero    .and.  dabs(ff)>absaccnr)then

          n2it=n2it+1
          if(n2it<=3) pcurrent=half*(pnew+pcurrent)
          if(n2it>3)then
             pright =pcurrent
             pleft=pprev
             pcurrent=half*(pleft+pright)
             Dicho:  do ni3=1,maxitnr
                !===================!
                xicurrent=tau+d+pcurrent
                v2=sqrs/xicurrent**2
                lfac2inv=one - v2

                if(lfac2inv>zero)then
                   lfac=one/dsqrt(lfac2inv)
                else
                   ierror=4
                   return
                endif
                !===================!

                !== calculation done using the EOS ==!
                call Bisection_Enthalpy(pnew,lfac2inv,d,xicurrent,h)
                Nff=-xicurrent*lfac2inv + h
                !=======================================!
                !==== Iterate ====!
                if(ff * Nff < zero)then
                   pleft=pcurrent
                else
                   pright=pcurrent
                endif

                pcurrent=half*(pleft+pright)
                !==================!

                !=== The iteration converged ===!
                if(2.0d0*dabs(pleft-pright)/(pleft+pright)< absaccnr &
                     .or. dabs(ff)<absaccnr)then
                   pnew=pcurrent
                   exit LoopNR
                endif
                !==============================!

                !==============================!

                !=== conserve the last value of Nff ===!
                ff=Nff
                !======================================!
             enddo    Dicho
          endif

       else
          !====== There is no problems, continue the NR iteration ======!
          pprev=pcurrent
          pcurrent=pnew
          !=============================================================!
       endif


       !=== keep the values of the 2 last ff ===!
       oldff2=oldff1
       oldff1=ff
       !========================================!
    enddo LoopNR

    if(niiter==maxitnr)then
       ierror=2
       return
    endif

    if(pcurrent<minp) then
       ierror=3
       return
    endif

    !------------------------------!
    xi=tau+d+pcurrent
    v2=sqrs/xicurrent**2
    lfac2inv=one - v2
    if(lfac2inv>zero) then
       lfac=one/dsqrt(lfac2inv)
    else
       ierror=4
       return
    endif

  end subroutine con2primHydro

  !> pointwise evaluations used in con2prim
  !> compute pointwise value for pressure p and dpdxi
  subroutine FuncPressure(xicurrent,lfac,d,dlfacdxi,p,dpdxi)

    double precision, intent(in)         :: xicurrent,lfac,d,dlfacdxi
    double precision, intent(out)        :: p,dpdxi
    ! .. local ..
    double precision                     :: rho,h,E,dhdxi,rhotoE
    double precision                     :: dpdchi,dEdxi

    ! rhoh here called h
    h=xicurrent/(lfac**2)
    rho=d/lfac
    ! output pressure
    p = (h - rho)/gamma_to_gamma_1
    dpdchi = one/gamma_to_gamma_1
    dpdxi = dpdchi * one/lfac**2
    ! zero case dlfacdxi implies zero velocity (ssqr=0)
    if (dlfacdxi /= 0.0d0) &
          dpdxi = dpdxi  + dpdchi * ((d*lfac-2.0d0*xicurrent)/lfac**3) * dlfacdxi

  end subroutine FuncPressure

  !> pointwise evaluations used in con2prim
  !> returns enthalpy rho*h (h) and derivative d(rho*h)/dp (dhdp)
  subroutine FuncEnthalpy(pcurrent,lfac2inv,d,sqrs,xicurrent,dv2d2p,h,dhdp)

    double precision, intent(in) :: pcurrent,lfac2inv,d,sqrs,xicurrent,dv2d2p
    double precision, intent(out):: h,dhdp

    ! local
    double precision:: rho,E_th,E,dE_thdp,dEdp

    rho=d*dsqrt(lfac2inv)
    h = rho + gamma_to_gamma_1 * pcurrent
    dhdp = gamma_to_gamma_1 + d/dsqrt(lfac2inv)*sqrs/xicurrent**3
  end subroutine FuncEnthalpy

  !> pointwise evaluations used in con2prim
  !> returns enthalpy rho*h (h)
  subroutine Bisection_Enthalpy(pcurrent,lfac2inv,d,xicurrent,h)

    double precision, intent(in) :: pcurrent,lfac2inv,d,xicurrent
    double precision, intent(out):: h

    ! local
    double precision:: rho,E_th,E

    rho=d*dsqrt(lfac2inv)
    h = rho + gamma_to_gamma_1 * pcurrent

    return
  end subroutine Bisection_Enthalpy

end module mod_srhd_phys
