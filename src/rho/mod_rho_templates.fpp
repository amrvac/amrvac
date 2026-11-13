#:mute
#:include "../mod_gpu_directives.fpp"
#:endmute
#:if PHYS == 'rho'

!> Scalar advection, d(rho)/dt + div(rho v) = 0, with a constant Cartesian
!> velocity rho_v from &rho_list. Port of MPI-AMRVAC's src/rho/mod_rho_phys.t.
!> Primitive and conservative variables coincide.
!>
!> In a spherical (r, theta, phi) or cylindrical (r, z, phi) build the flux
!> uses the local component of rho_v normal to the face, evaluated at the
!> position the finite-volume kernel hands get_flux: the face in the normal
!> direction, the cell's volume barycentre in the others. That is second
!> order, but - unlike upstream's exact face averages - not exactly
!> divergence free, so a uniform rho in a uniform flow is kept uniform only
!> to O(h^2), not to round-off. The components depend on the angles alone,
!> so a logarithmic radius (LOG_RADIUS) needs nothing extra: the kernel
!> already hands over physical radii.

!> rho_vn(x, d, vn): statements setting vn to the component of the constant
!> Cartesian velocity rho_v along coordinate direction d at position x
!> (physical coordinates, angles in radians). Expanded inline into every
!> routine that needs it, so that no extra device routine is involved.
#:def rho_vn(x, d, vn)
#:if GEOM == 'spherical'
  select case (${d}$)
  case (1)   ! v_r
     ${vn}$ = sin(${x}$(2))*(rho_v(1)*cos(${x}$(3)) + rho_v(2)*sin(${x}$(3))) &
          + rho_v(3)*cos(${x}$(2))
  case (2)   ! v_theta
     ${vn}$ = cos(${x}$(2))*(rho_v(1)*cos(${x}$(3)) + rho_v(2)*sin(${x}$(3))) &
          - rho_v(3)*sin(${x}$(2))
  case default   ! v_phi
     ${vn}$ = -rho_v(1)*sin(${x}$(3)) + rho_v(2)*cos(${x}$(3))
  end select
#:elif GEOM == 'cylindrical'
  select case (${d}$)
  case (1)   ! v_r
     ${vn}$ = rho_v(1)*cos(${x}$(3)) + rho_v(2)*sin(${x}$(3))
  case (2)   ! v_z
     ${vn}$ = rho_v(3)
  case default   ! v_phi
     ${vn}$ = -rho_v(1)*sin(${x}$(3)) + rho_v(2)*cos(${x}$(3))
  end select
#:else
  ${vn}$ = rho_v(${d}$)
#:endif
#:enddef

#:def phys_vars()

  integer, parameter :: dp = kind(0.0d0)
  integer, parameter, public              :: nw_phys=1
  integer, parameter, public              :: nw_flux=1
  !> Number of leading variables the MUSCL reconstruction limits (1:nw_recon)
  integer, parameter, public              :: nw_recon=nw_flux

  !> Index of the advected density (in the w array)
  integer, public                         :: rho_
  ${GPU_DECLARE_CREATE('rho_')}$

  !> Constant Cartesian advection velocity
  double precision, public                :: rho_v(1:ndim) = 1.0d0
  ${GPU_DECLARE_COPYIN('rho_v')}$

#:enddef

#:def read_params()
  !> Read this module's parameters from a file
  subroutine read_params(files)
    use mod_global_parameters
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /rho_list/ rho_v

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, rho_list, end=111)
111    close(unitpar)
    end do

    if (all(rho_v == 0.0d0)) call mpistop("rho_v = 0: nothing is advected &
       &and the CFL time step would be infinite")

    ${GPU_UPDATE_DEVICE('rho_v')}$

  end subroutine read_params
#:enddef

#:def phys_activate()
  subroutine phys_activate()
    call phys_init()
  end subroutine phys_activate
#:enddef

#:def phys_units()
  subroutine phys_units()
    ! Linear advection is scale free: the units stay at their default of 1.
  end subroutine phys_units
#:enddef

#:def phys_init()
  !> Initialize the module
  subroutine phys_init()
    use mod_global_parameters

    call read_params(par_files)
    call phys_units()

    phys_energy       = .false.
    phys_total_energy = .false.
    phys_internal_e   = .false.
    ${GPU_UPDATE_DEVICE('physics_type, phys_energy, phys_total_energy, phys_internal_e')}$

    use_particles = .false.

    ! the single flux variable
    rho_ = var_set_rho()
    ${GPU_UPDATE_DEVICE('rho_')}$

    ! the unsplit finite-volume update needs no corner ghost cells
    phys_req_diagonal = .false.
    ${GPU_UPDATE_DEVICE('phys_req_diagonal')}$

    ! set number of variables which need update ghostcells
    nwgc = nwflux
    ${GPU_UPDATE_DEVICE('nwgc')}$

    if (.not. allocated(flux_type)) then
       allocate(flux_type(ndir, nw_flux))
       flux_type = flux_default
    else if (any(shape(flux_type) /= [ndir, nw_flux])) then
       call mpistop("phys_check error: flux_type has wrong shape")
    end if
    ${GPU_UPDATE_DEVICE('flux_type')}$

    ! no vector variables
    nvector = 0
    ${GPU_UPDATE_DEVICE('nvector')}$

  end subroutine phys_init
#:enddef

#:def phys_get_dt()
  !> No time-step limit beyond the CFL condition built from get_cmax
  subroutine phys_get_dt(w, x, dx, dtnew)
    ${GPU_ROUTINE_SEQ()}$
    real(dp), intent(in)   :: w(nw_phys), x(1:ndim), dx(1:ndim)
    real(dp), intent(out)  :: dtnew

    dtnew = huge(1.0d0)

  end subroutine phys_get_dt
#:enddef

#:def addsource_local()
!> The transport equation has no source terms
subroutine addsource_local(qdt, dtfactor, qtC, wCT, wCTprim, qt, wnew, x, dr, &
    qsourcesplit)
  ${GPU_ROUTINE_SEQ()}$
  real(dp), intent(in)     :: qdt, dtfactor, qtC, qt
  real(dp), intent(in)     :: wCT(nw_phys), wCTprim(nw_phys)
  real(dp), intent(in)     :: x(1:ndim), dr(ndim)
  real(dp), intent(inout)  :: wnew(nw_phys)
  logical, intent(in)      :: qsourcesplit

end subroutine addsource_local
#:enddef

#:if GEOM == 'spherical' or GEOM == 'cylindrical'
#:def addsource_geometry()
!> No curvature terms: rho is a scalar, and with the exact face areas and
!> volumes of fill_geometry_device the flux divergence of a scalar needs no
!> correction in any coordinate system (as for ffhd). Defined, though, since
!> the default in mod_physics_dummies stops a curvilinear run.
subroutine addsource_geometry(qdt, wprim, wnew, x, dAdV)
  ${GPU_ROUTINE_SEQ()}$
  real(dp), intent(in)     :: qdt
  real(dp), intent(in)     :: wprim(nw_phys)
  real(dp), intent(in)     :: x(1:ndim)
  real(dp), intent(in)     :: dAdV(1:ndim)
  real(dp), intent(inout)  :: wnew(nw_phys)

end subroutine addsource_geometry
#:enddef
#:endif

#:def to_primitive()
!> Primitive and conservative variables are identical
pure subroutine to_primitive(u)
  ${GPU_ROUTINE_SEQ()}$
  real(dp), intent(inout) :: u(nw_phys)

end subroutine to_primitive
#:enddef

#:def to_conservative()
!> Primitive and conservative variables are identical
pure subroutine to_conservative(u)
  ${GPU_ROUTINE_SEQ()}$
  real(dp), intent(inout) :: u(nw_phys)

end subroutine to_conservative
#:enddef

#:def get_flux()
subroutine get_flux(u, xC, flux_dim, flux)
  ${GPU_ROUTINE_SEQ()}$
  real(dp), intent(in)  :: u(nw_phys)
  real(dp), intent(in)  :: xC(ndim)
  integer, intent(in)   :: flux_dim
  real(dp), intent(out) :: flux(nw_flux)
  real(dp)              :: vn

  @:rho_vn(xC, flux_dim, vn)
  flux(iw_rho) = u(iw_rho) * vn

end subroutine get_flux
#:enddef

#:def get_cmax()
pure real(dp) function get_cmax(u, x, flux_dim) result(wC)
  ${GPU_ROUTINE_SEQ()}$
  real(dp), intent(in)  :: u(nw_phys)
  real(dp), intent(in)  :: x(1:ndim)
  integer, intent(in)   :: flux_dim
  real(dp)              :: vn

  @:rho_vn(x, flux_dim, vn)
  wC = abs(vn)

end function get_cmax
#:enddef

#:def estimate_speeds_minmax()
!> The single characteristic speed is v_n. Bounding it by zero makes the HLL
!> flux exact upwinding: F = rho_L v_n for v_n >= 0, rho_R v_n otherwise.
subroutine estimate_speeds_minmax(uL, uR, xC, flux_dim, wL, wR)
  ${GPU_ROUTINE_SEQ()}$
  real(dp), intent(in)  :: uL(nw_phys), uR(nw_phys)
  real(dp), intent(in)  :: xC(ndim)
  integer, intent(in)   :: flux_dim
  real(dp), intent(out) :: wL, wR
  real(dp)              :: vn

  @:rho_vn(xC, flux_dim, vn)
  wL = min(vn, 0.0_dp)
  wR = max(vn, 0.0_dp)

end subroutine estimate_speeds_minmax
#:enddef

#:endif
