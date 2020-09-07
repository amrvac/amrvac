!> Module containing the physics routines for scalar nonlinear equation
module mod_nonlinear_phys

  implicit none
  private

  !> index of the single scalar unknown
  integer, protected, public :: rho_       = 1

  !> switch between burgers (i.e. rho**2) 
  !> or nonconvex flux (i.e. rho**3)
  integer, protected, public :: nonlinear_flux_type = 1

  !> whether the KdV source term is added
  logical, protected, public :: kdv_source_term = .false.

  !> Whether particles module is added
  logical, public, protected              :: nonlinear_particles = .false.

  ! Public methods
  public :: nonlinear_phys_init
  public :: nonlinear_get_v

contains

  !> Read this module's parameters from a file
  subroutine nonlinear_params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /nonlinear_list/ nonlinear_flux_type, kdv_source_term, nonlinear_particles

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status='old')
       read(unitpar, nonlinear_list, end=111)
111    close(unitpar)
    end do

  end subroutine nonlinear_params_read

  !> Write this module's parameters to a snapshot
  subroutine nonlinear_write_info(fh)
  ! for nonlinear scalar equation, nothing to write
  ! note: this is info only stored at end of dat files, 
  !       is never read/used for restarts, only expects
  !       an integer (number of parameters) and 
  !       corresponding double values and character names
  !       and is meant for use in the python tools
    use mod_global_parameters
    integer, intent(in)                 :: fh
    integer, dimension(MPI_STATUS_SIZE) :: st
    integer                             :: er

    ! Write zero parameters
    call MPI_FILE_WRITE(fh, 0, 1, MPI_INTEGER, st, er)

  end subroutine nonlinear_write_info

  subroutine nonlinear_phys_init()
    use mod_global_parameters
    use mod_physics
    use mod_kdv, only: kdv_init
    use mod_particles, only: particles_init

    call nonlinear_params_read(par_files)

    physics_type = "nonlinear"
    phys_energy  = .false.
    ! Whether diagonal ghost cells are required for the physics
    phys_req_diagonal = .false.
    use_particles = nonlinear_particles

    rho_ = var_set_rho()

    ! Check whether custom flux types have been defined
    if (.not. allocated(flux_type)) then
       allocate(flux_type(ndir, nw))
       flux_type = flux_default
    else if (any(shape(flux_type) /= [ndir, nw])) then
       call mpistop("phys_check error: flux_type has wrong shape")
    end if

    phys_get_cmax        => nonlinear_get_cmax
    phys_get_cbounds     => nonlinear_get_cbounds
    phys_get_flux        => nonlinear_get_flux
    phys_add_source_geom => nonlinear_add_source_geom
    phys_add_source      => nonlinear_add_source
    phys_to_conserved    => nonlinear_to_conserved
    phys_to_primitive    => nonlinear_to_primitive
    phys_get_dt          => nonlinear_get_dt
    phys_write_info      => nonlinear_write_info

    if (kdv_source_term) call kdv_init()

    ! Initialize particles module
    if (nonlinear_particles) then
       call particles_init()
       phys_req_diagonal = .true.
    end if

  end subroutine nonlinear_phys_init

  subroutine nonlinear_to_conserved(ixI^L, ixO^L, w, x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:^ND)

    ! Do nothing (primitive and conservative are equal for nonlinear module)
  end subroutine nonlinear_to_conserved

  subroutine nonlinear_to_primitive(ixI^L, ixO^L, w, x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:^ND)

    ! Do nothing (primitive and conservative are equal for nonlinear module)
  end subroutine nonlinear_to_primitive

  subroutine nonlinear_get_v(w, x, ixI^L, ixO^L, idim, v)
    use mod_global_parameters
    integer, intent(in)           :: ixI^L, ixO^L, idim
    double precision, intent(in)  :: w(ixI^S, nw), x(ixI^S, 1:^ND)
    double precision, intent(out) :: v(ixI^S)

   select case(nonlinear_flux_type)
    case(1)
       v(ixO^S)=w(ixO^S,rho_)
    case(2)
       v(ixO^S)=3.0d0*w(ixO^S,rho_)**2
    case default
       call mpistop('Undefined fluxtype: set nonlinear_flux_type to 1 or 2')
    end select

  end subroutine nonlinear_get_v

  subroutine nonlinear_get_cmax(w, x, ixI^L, ixO^L, idim, cmax)
    use mod_global_parameters
    integer, intent(in)                       :: ixI^L, ixO^L, idim
    double precision, intent(in)              :: w(ixI^S, nw), x(ixI^S, 1:^ND)
    double precision, intent(inout)           :: cmax(ixI^S)

    call nonlinear_get_v(w, x, ixI^L, ixO^L, idim, cmax)

    cmax(ixO^S) = abs(cmax(ixO^S))

  end subroutine nonlinear_get_cmax

  subroutine nonlinear_get_cbounds(wLC, wRC, wLp, wRp, x, ixI^L, ixO^L, idim, cmax, cmin)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L, idim
    double precision, intent(in)    :: wLC(ixI^S, nw), wRC(ixI^S,nw)
    double precision, intent(in)    :: wLp(ixI^S, nw), wRp(ixI^S,nw)
    double precision, intent(in)    :: x(ixI^S, 1:^ND)
    double precision, intent(inout) :: cmax(ixI^S)
    double precision, intent(inout), optional :: cmin(ixI^S)

    double precision :: wmean(ixI^S,nw)

    ! since get_v depends on w, the first argument should be some average over the
    ! left and right state
    wmean(ixO^S,1:nwflux)=0.5d0*(wLC(ixO^S,1:nwflux)+wRC(ixO^S,1:nwflux))
    call nonlinear_get_v(wmean, x, ixI^L, ixO^L, idim, cmax)

    if (present(cmin)) then
       cmin(ixO^S) = min(cmax(ixO^S), zero)
       cmax(ixO^S) = max(cmax(ixO^S), zero)
    else
       cmax(ixO^S) = maxval(abs(cmax(ixO^S)))
    end if

  end subroutine nonlinear_get_cbounds

  subroutine nonlinear_get_dt(w, ixI^L, ixO^L, dtnew, dx^D, x)
    use mod_global_parameters
    use mod_kdv, only: kdv_get_dt

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: dx^D, x(ixI^S, 1:^ND)
    double precision, intent(in)    :: w(ixI^S, 1:nw)
    double precision, intent(inout) :: dtnew

    dtnew = bigdouble

    if(kdv_source_term) then
      call kdv_get_dt(w, ixI^L, ixO^L, dtnew, dx^D, x)
    endif
  end subroutine nonlinear_get_dt

  ! here we select the flux according to the nonlinear_flux_type parameter
  subroutine nonlinear_get_flux(wC, w, x, ixI^L, ixO^L, idim, f)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L, idim
    double precision, intent(in)    :: wC(ixI^S, 1:nw)
    double precision, intent(in)    :: w(ixI^S, 1:nw)
    double precision, intent(in)    :: x(ixI^S, 1:^ND)
    double precision, intent(out)   :: f(ixI^S, nwflux)

    select case(nonlinear_flux_type)
    case(1)
       f(ixO^S,rho_)=half*w(ixO^S,rho_)**2
    case(2)
       f(ixO^S,rho_)=w(ixO^S,rho_)**3
    case default
       call mpistop('Undefined fluxtype: set nonlinear_flux_type to 1 or 2')
    end select

  end subroutine nonlinear_get_flux

  subroutine nonlinear_add_source_geom(qdt, ixI^L, ixO^L, wCT, w, x)

    ! Add geometrical source terms to w
    ! There are no geometrical source terms 

    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: qdt, x(ixI^S, 1:^ND)
    double precision, intent(inout) :: wCT(ixI^S, 1:nw), w(ixI^S, 1:nw)

  end subroutine nonlinear_add_source_geom

  subroutine nonlinear_add_source(qdt,ixI^L,ixO^L,wCT,w,x,qsourcesplit,active)
    use mod_global_parameters
    use mod_kdv, only: kdv_add_source

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixI^S, 1:nw), x(ixI^S, 1:ndim)
    double precision, intent(inout) :: w(ixI^S, 1:nw)
    logical, intent(in)             :: qsourcesplit
    logical, intent(inout)          :: active

    if(kdv_source_term) then
      call kdv_add_source(qdt,ixI^L,ixO^L,wCT,w,x,qsourcesplit,active)
    end if

  end subroutine nonlinear_add_source

end module mod_nonlinear_phys
