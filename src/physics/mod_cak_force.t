!> Module to include CAK radiation line force in (magneto)hydrodynamic models
!> Computes both the force from free electrons and the force from an ensemble of
!> lines (various possibilities for the latter).
!> There is an option to only simulate the pure radial CAK force (with various
!> corrections applied) as well as the full vector CAK force. Depending on the
!> chosen option additional output are the CAK line force component(s) and,
!> when doing a 1-D radial force, the finite disc factor.
!>
!> USAGE:
!>
!>  1. Include a cak_list in the .par file and activate (m)hd_cak_force in the
!>     (m)hd_list
!>  2. Create a mod_usr.t file for the problem with appropriate initial and
!>     boundary conditions
!>  3. In the mod_usr.t header call the mod_cak_force module to have access to
!>     global variables from mod_cak_force, which may be handy for printing or
!>     the computation of other variables inside mod_usr.t
!>  4. In usr_init of mod_usr.t call the set_cak_force_norm routine and pass
!>     along the stellar radius and wind temperature---this is needed to
!>     correctly compute the (initial) force normalisation inside mod_cak_force
!>  5. Ensure that the order of calls in usr_init is similar as for test problem
!>     CAKwind_spherical_1D: first reading usr.par list; then set unit scales;
!>     then call (M)HD_activate; then call set_cak_force_norm. This order avoids
!>     an incorrect force normalisation and code crash
!>
!> Developed by Florian Driessen (2022)
module mod_cak_force
  use mod_physics, only: phys_get_pthermal, physics_type
  implicit none

  !> Line-ensemble parameters in the Gayley (1995) formalism
  real(8), public :: cak_alpha, gayley_qbar, gayley_q0

  !> Switch to choose between the 1-D CAK line force options
  integer :: cak_1d_opt

  ! Avoid magic numbers in code for 1-D CAK line force option
  integer, parameter, private :: radstream=0, fdisc=1, fdisc_cutoff=2

  !> To treat source term in split or unsplit (default) fashion
  logical :: cak_split=.false.

  !> To activate the original CAK 1-D line force computation
  logical :: cak_1d_force=.false.

  !> To activate the vector CAK line force computation
  logical :: cak_vector_force=.false.

  !> To activate the pure radial vector CAK line force computation
  logical :: fix_vector_force_1d=.false.

  !> Amount of rays in radiation polar and radiation azimuthal direction
  integer :: nthetaray, nphiray

  !> Ray positions + weights for impact parameter and azimuthal radiation angle
  real(8), allocatable, private :: ay(:), wy(:), aphi(:), wphi(:)

  !> The adiabatic index
  real(8), private :: cak_gamma

  !> Variables needed to compute force normalisation fnorm in initialisation
  real(8), private :: lum, dlum, drstar, dke, dclight

  !> To enforce a floor temperature when doing adiabatic (M)HD
  real(8), private :: tfloor

  !> Extra slots to store quantities in w-array
  integer :: gcak1_, gcak2_, gcak3_, fdf_

  !> Public method
  public :: set_cak_force_norm
  
contains

  !> Read this module's parameters from a file
  subroutine cak_params_read(files)
    use mod_global_parameters, only: unitpar

    character(len=*), intent(in) :: files(:)

    ! Local variable
    integer :: n

    namelist /cak_list/ cak_alpha, gayley_qbar, gayley_q0, cak_1d_opt, &
                        cak_split, cak_1d_force, cak_vector_force, &
                        nphiray, nthetaray, fix_vector_force_1d

    do n = 1,size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, cak_list, end=111)
       111 close(unitpar)
    enddo

  end subroutine cak_params_read

  !> Initialize the module
  subroutine cak_init(phys_gamma)
    use mod_global_parameters

    real(8), intent(in) :: phys_gamma

    cak_gamma = phys_gamma

    ! Set some defaults when user does not
    cak_alpha   = 0.65d0
    gayley_qbar = 2000.0d0
    gayley_q0   = 2000.0d0
    cak_1d_opt  = 1
    nthetaray   = 6
    nphiray     = 6

    call cak_params_read(par_files)

    if (cak_1d_force) then
      gcak1_ = var_set_extravar("gcak1", "gcak1")
      fdf_   = var_set_extravar("fdfac", "fdfac")
    endif

    if (cak_vector_force) then
      gcak1_ = var_set_extravar("gcak1", "gcak1")
      gcak2_ = var_set_extravar("gcak2", "gcak2")
      gcak3_ = var_set_extravar("gcak3", "gcak3")
      call rays_init(nthetaray,nphiray)
    endif

    ! Some sanity checks
    if ((cak_alpha <= 0.0d0) .or. (cak_alpha > 1.0d0)) then
      call mpistop('CAK error: choose alpha in [0,1[')
    endif

    if ((gayley_qbar < 0.0d0) .or. (gayley_q0 < 0.0d0)) then
      call mpistop('CAK error: chosen Qbar or Q0 is < 0')
    endif

    if (cak_1d_force .and. cak_vector_force) then
      call mpistop('CAK error: choose either 1-D or vector force')
    endif

  end subroutine cak_init

  !> Compute some (unitless) variables for CAK force normalisation
  subroutine set_cak_force_norm(rstar,twind)
    use mod_global_parameters
    use mod_constants

    real(8), intent(in) :: rstar, twind

    lum     = 4.0d0*dpi * rstar**2.0d0 * const_sigma * twind**4.0d0
    dke     = const_kappae * unit_density * unit_length
    dclight = const_c/unit_velocity
    dlum    = lum/(unit_density * unit_length**5.0d0 / unit_time**3.0d0)
    drstar  = rstar/unit_length
    tfloor  = twind/unit_temperature

  end subroutine set_cak_force_norm
  
  !> w[iw]=w[iw]+qdt*S[wCT,qtC,x] where S is the source based on wCT within ixO
  subroutine cak_add_source(qdt,ixI^L,ixO^L,wCT,w,x,energy,qsourcesplit,active)
    use mod_global_parameters

    integer, intent(in)    :: ixI^L, ixO^L
    real(8), intent(in)    :: qdt, x(ixI^S,1:ndim), wCT(ixI^S,1:nw)
    real(8), intent(inout) :: w(ixI^S,1:nw)
    logical, intent(in)    :: energy, qsourcesplit
    logical, intent(inout) :: active

    ! Local variables
    integer :: idir
    real(8) :: gl(ixO^S,1:3), ge(ixO^S), etherm(ixI^S), emin(ixI^S)

    ! By default add source in unsplit fashion together with the fluxes
    if (qsourcesplit .eqv. cak_split) then

      ! Thomson force
      call get_gelectron(ixI^L,ixO^L,wCT,x,ge)

      ! CAK line force
      gl(ixO^S,1:3) = 0.0d0

      if (cak_1d_force) then
        call get_cak_force_radial(ixI^L,ixO^L,wCT,w,x,gl)
      elseif (cak_vector_force) then
        call get_cak_force_vector(ixI^L,ixO^L,wCT,w,x,gl)
      else
        call mpistop("No valid force option")
      endif

      ! Update conservative vars: w = w + qdt*gsource
      do idir = 1,ndir
        if (idir == 1) gl(ixO^S,idir) = gl(ixO^S,idir) + ge(ixO^S)
        
        w(ixO^S,iw_mom(idir)) = w(ixO^S,iw_mom(idir)) &
                                + qdt * gl(ixO^S,idir) * wCT(ixO^S,iw_rho)
                                
        if (energy) then
          w(ixO^S,iw_e) = w(ixO^S,iw_e) + qdt * gl(ixO^S,idir) * wCT(ixO^S,iw_mom(idir))
          
          ! Impose fixed floor temperature to mimic stellar heating
          call phys_get_pthermal(w,x,ixI^L,ixO^L,etherm)
          etherm(ixO^S) = etherm(ixO^S) / (cak_gamma - 1.0d0)
          emin(ixO^S)   = w(ixO^S,iw_rho)*tfloor / (cak_gamma - 1.0d0)
          
          where (etherm < emin)
            w(ixO^S,iw_e) = w(ixO^S,iw_e) - etherm(ixO^S) + emin(ixO^S)
          endwhere
        endif
      enddo
    endif

  end subroutine cak_add_source

  !> 1-D CAK line force in the Gayley line-ensemble distribution parametrisation
  subroutine get_cak_force_radial(ixI^L,ixO^L,wCT,w,x,gcak)
    use mod_global_parameters

    integer, intent(in)    :: ixI^L, ixO^L
    real(8), intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    real(8), intent(inout) :: w(ixI^S,1:nw)
    real(8), intent(inout) :: gcak(ixO^S,1:3)
  
    ! Local variables
    real(8) :: vr(ixI^S), dvrdr(ixO^S)
    real(8) :: beta_fd(ixO^S), fdfac(ixO^S), taus(ixO^S), ge(ixO^S)
  
    vr(ixI^S) = wCT(ixI^S,iw_mom(1)) / wCT(ixI^S,iw_rho)
    call get_velocity_gradient(ixI^L,ixO^L,vr,x,1,dvrdr)

    if (physics_type == 'hd') then
      ! Monotonic flow to avoid multiple resonances and radiative coupling
      dvrdr(ixO^S) = abs(dvrdr(ixO^S))
    else
      ! Allow material to fallback to the star in a magnetosphere model
      dvrdr(ixO^S) = max(dvrdr(ixO^S), smalldouble)
    endif
  
    ! Thomson force
    call get_gelectron(ixI^L,ixO^L,wCT,x,ge)

    ! Sobolev optical depth for line ensemble (tau = Qbar * t_r) and the force
    select case (cak_1d_opt)
    case(radstream, fdisc)
      taus(ixO^S)   = gayley_qbar * dke * dclight * wCT(ixO^S,iw_rho)/dvrdr(ixO^S)
      gcak(ixO^S,1) = gayley_qbar/(1.0d0 - cak_alpha) &
                      * ge(ixO^S)/taus(ixO^S)**cak_alpha

    case(fdisc_cutoff)
      taus(ixO^S)   = gayley_q0 * dke * dclight * wCT(ixO^S,iw_rho)/dvrdr(ixO^S)
      gcak(ixO^S,1) = gayley_qbar * ge(ixO^S)                                  &
                      * ( (1.0d0 + taus(ixO^S))**(1.0d0 - cak_alpha) - 1.0d0 ) &
                      / ( (1.0d0 - cak_alpha) * taus(ixO^S) )
    case default
      call mpistop("Error in force computation.")
    end select

    ! Finite disk factor parameterisation (Owocki & Puls 1996)
    beta_fd(ixO^S) = ( 1.0d0 - vr(ixO^S)/(x(ixO^S,1) * dvrdr(ixO^S)) ) &
                      * (drstar/x(ixO^S,1))**2.0d0

    select case (cak_1d_opt)
    case(radstream)
      fdfac(ixO^S) = 1.0d0
    case(fdisc, fdisc_cutoff)
      where (beta_fd(ixO^S) >= 1.0d0)
        fdfac(ixO^S) = 1.0d0/(1.0d0 + cak_alpha)
      elsewhere (beta_fd(ixO^S) < -1.0d10)
        fdfac(ixO^S) = abs(beta_fd(ixO^S))**cak_alpha / (1.0d0 + cak_alpha)
      elsewhere (abs(beta_fd) > 1.0d-3)
        fdfac(ixO^S) = (1.0d0 - (1.0d0 - beta_fd(ixO^S))**(1.0d0 + cak_alpha)) &
                       / (beta_fd(ixO^S)*(1.0d0 + cak_alpha))
      elsewhere
        fdfac(ixO^S) = 1.0d0 - 0.5d0*cak_alpha*beta_fd(ixO^S) &
                       * (1.0d0 + 1.0d0/3.0d0 * (1.0d0 - cak_alpha)*beta_fd(ixO^S))
      endwhere
    end select

    ! Correct radial line force for finite disc (if applicable)
    gcak(ixO^S,1) = gcak(ixO^S,1) * fdfac(ixO^S)
    gcak(ixO^S,2) = 0.0d0
    gcak(ixO^S,3) = 0.0d0
      
    ! Fill the nwextra slots for output
    w(ixO^S,gcak1_) = gcak(ixO^S,1)
    w(ixO^S,fdf_)   = fdfac(ixO^S)
    
  end subroutine get_cak_force_radial

  !> Vector CAK line force in the Gayley line-ensemble distribution parametrisation
  subroutine get_cak_force_vector(ixI^L,ixO^L,wCT,w,x,gcak)
    use mod_global_parameters
    use mod_usr_methods

    ! Subroutine arguments
    integer, intent(in)    :: ixI^L, ixO^L
    real(8), intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    real(8), intent(inout) :: w(ixI^S,1:nw)
    real(8), intent(inout) :: gcak(ixO^S,1:3)

    ! Local variables
    integer :: ix^D, itray, ipray
    real(8) :: a1, a2, a3, wyray, y, wpray, phiray, wtot, mustar, dvndn
    real(8) :: costp, costp2, sintp, cospp, sinpp, cott0
    real(8) :: vr(ixI^S), vt(ixI^S), vp(ixI^S)
    real(8) :: vrr(ixI^S), vtr(ixI^S), vpr(ixI^S)
    real(8) :: dvrdr(ixO^S), dvtdr(ixO^S), dvpdr(ixO^S)
    real(8) :: dvrdt(ixO^S), dvtdt(ixO^S), dvpdt(ixO^S)
    real(8) :: dvrdp(ixO^S), dvtdp(ixO^S), dvpdp(ixO^S)
    
    ! Initialisation to have full velocity strain tensor expression at all times
    vt(ixO^S) = 0.0d0; vtr(ixO^S) = 0.0d0
    vp(ixO^S) = 0.0d0; vpr(ixO^S) = 0.0d0
    cott0 = 0.0d0
    dvrdr(ixO^S) = 0.0d0; dvtdr(ixO^S) = 0.0d0; dvpdr(ixO^S) = 0.0d0
    dvrdt(ixO^S) = 0.0d0; dvtdt(ixO^S) = 0.0d0; dvpdt(ixO^S) = 0.0d0
    dvrdp(ixO^S) = 0.0d0; dvtdp(ixO^S) = 0.0d0; dvpdp(ixO^S) = 0.0d0

    ! Populate velocity field(s) depending on dimensions and directions
    vr(ixI^S)  = wCT(ixI^S,iw_mom(1)) / wCT(ixI^S,iw_rho)
    vrr(ixI^S) = vr(ixI^S) / x(ixI^S,1)

    {^NOONED
    vt(ixI^S)  = wCT(ixI^S,iw_mom(2)) / wCT(ixI^S,iw_rho)
    vtr(ixI^S) = vt(ixI^S) / x(ixI^S,1)
    
    if (ndir > 2) then
      vp(ixI^S)  = wCT(ixI^S,iw_mom(3)) / wCT(ixI^S,iw_rho)
      vpr(ixI^S) = vp(ixI^S) / x(ixI^S,1)
    endif
    }
    
    ! Derivatives of velocity field in each coordinate direction (r=1,t=2,p=3)
    call get_velocity_gradient(ixI^L,ixO^L,vr,x,1,dvrdr)
    
    {^NOONED
    call get_velocity_gradient(ixI^L,ixO^L,vr,x,2,dvrdt)
    call get_velocity_gradient(ixI^L,ixO^L,vt,x,1,dvtdr)
    call get_velocity_gradient(ixI^L,ixO^L,vt,x,2,dvtdt)

    if (ndir > 2) then
      call get_velocity_gradient(ixI^L,ixO^L,vp,x,1,dvpdr)
      call get_velocity_gradient(ixI^L,ixO^L,vp,x,2,dvpdt)
    endif
    }
    {^IFTHREED
    call get_velocity_gradient(ixI^L,ixO^L,vr,x,3,dvrdp)
    call get_velocity_gradient(ixI^L,ixO^L,vt,x,3,dvtdp)
    call get_velocity_gradient(ixI^L,ixO^L,vp,x,3,dvpdp)
    }

    ! Get total acceleration from all rays at a certain grid point
    {do ix^DB=ixOmin^DB,ixOmax^DB\}
      ! Loop over the rays; first theta then phi radiation angle
      ! Get weights from current ray and their position
      do itray = 1,nthetaray
        wyray  = wy(itray)
        y      = ay(itray)

        do ipray = 1,nphiray
          wpray = wphi(ipray)
          phiray  = aphi(ipray)

          ! Redistribute the phi rays by a small offset
          ! if (mod(itp,3) == 1) then
          !   phip = phip + dphi/3.0d0
          ! elseif (mod(itp,3) == 2) then
          !   phip = phip - dphi/3.0d0
          ! endif
          
          ! === Geometrical factors ===
          ! Make y quadrature linear in mu, not mu**2; better for gtheta,gphi
          ! y -> mu quadrature is preserved; y=0 <=> mu=1; y=1 <=> mu=mustar
          mustar = sqrt(max(1.0d0 - (drstar/x(ix^D,1))**2.0d0, 0.0d0))
          costp  = 1.0d0 - y*(1.0d0 - mustar)
          costp2 = costp*costp
          sintp  = sqrt(max(1.0d0 - costp2, 0.0d0))
          sinpp  = sin(phiray)
          cospp  = cos(phiray)
          {^NOONED cott0  = cos(x(ix^D,2))/sin(x(ix^D,2))}

          ! More weight close to star, less farther away
          wtot  = wyray * wpray * (1.0d0 - mustar)

          ! Convenients a la Cranmer & Owocki (1995)
          a1 = costp
          a2 = sintp * cospp
          a3 = sintp * sinpp

          ! Get total velocity gradient for one ray with given (theta', phi')
          dvndn = a1*a1 * dvrdr(ix^D) + a2*a2 * (dvtdt(ix^D) + vrr(ix^D))  &
                 + a3*a3 * (dvpdp(ix^D) + cott0 * vtr(ix^D) + vrr(ix^D))   &
                 + a1*a2 * (dvtdr(ix^D) + dvrdt(ix^D) - vtr(ix^D))         &
                 + a1*a3 * (dvpdr(ix^D) + dvrdp(ix^D) - vpr(ix^D))         &
                 + a2*a3 * (dvpdt(ix^D) + dvtdp(ix^D) - cott0 * vpr(ix^D))

          ! No multiple resonances in CAK
          dvndn = abs(dvndn)

          ! Convert gradient back from wind coordinates (r',theta',phi') to
          ! stellar coordinates (r,theta,phi)
          gcak(ix^D,1) = gcak(ix^D,1) + (dvndn/wCT(ix^D,iw_rho))**cak_alpha * a1 * wtot
          gcak(ix^D,2) = gcak(ix^D,2) + (dvndn/wCT(ix^D,iw_rho))**cak_alpha * a2 * wtot
          gcak(ix^D,3) = gcak(ix^D,3) + (dvndn/wCT(ix^D,iw_rho))**cak_alpha * a3 * wtot
        enddo
      enddo
    {enddo\}

    ! Normalisation for line force
    ! NOTE: extra 1/pi factor comes from integration in radiation Phi angle
    gcak(ixO^S,:) = (dke*gayley_qbar)**(1.0d0 - cak_alpha)/(1.0d0 - cak_alpha)    &
                    * dlum/(4.0d0*dpi*drstar**2.0d0 * dclight**(1.0d0+cak_alpha)) &
                    * gcak(ixO^S,:)/dpi

    if (fix_vector_force_1d) then
      gcak(ixO^S,2) = 0.0d0
      gcak(ixO^S,3) = 0.0d0
    endif
          
    ! Fill the nwextra slots for output
    w(ixO^S,gcak1_) = gcak(ixO^S,1)
    w(ixO^S,gcak2_) = gcak(ixO^S,2)
    w(ixO^S,gcak3_) = gcak(ixO^S,3)

  end subroutine get_cak_force_vector
  
  !> Compute continuum radiation force from Thomson scattering
  subroutine get_gelectron(ixI^L,ixO^L,w,x,ge)
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    real(8), intent(in) :: w(ixI^S,1:nw), x(ixI^S,1:ndim)
    real(8), intent(out):: ge(ixO^S)

    ge(ixO^S) = dke * dlum/(4.0d0*dpi * dclight * x(ixO^S,1)**2.0d0)

  end subroutine get_gelectron

  !> Check time step for total radiation contribution
  subroutine cak_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
    use mod_global_parameters

    integer, intent(in)    :: ixI^L, ixO^L
    real(8), intent(in)    :: dx^D, x(ixI^S,1:ndim)
    real(8), intent(in)    :: w(ixI^S,1:nw)
    real(8), intent(inout) :: dtnew
    
    ! Local variables
    real(8) :: tdumr(ixO^S), tdumt(ixO^S), ge(ixO^S)
    real(8) :: dt_cakr, dt_cakt, dt_cakp

    call get_gelectron(ixI^L,ixO^L,w,x,ge)

    dt_cakr = bigdouble
    dt_cakt = bigdouble
    dt_cakp = bigdouble

    ! Get dt from line force that is saved in the w-array in nwextra slot
    tdumr(ixO^S) = sqrt( block%dx(ixO^S,1) / abs(ge(ixO^S) + w(ixO^S,gcak1_)) )
    dt_cakr      = courantpar * minval(tdumr(ixO^S))
    
    {^NOONED
    if (cak_vector_force) then
      tdumt(ixO^S) = sqrt( block%dx(ixO^S,1) * block%dx(ixO^S,2) / abs(w(ixO^S,gcak2_)) )
      dt_cakt      = courantpar * minval(tdumt(ixO^S))

      {^IFTHREED
      tdumt(ixO^S) = sqrt( block%dx(ixO^S,1) * sin(block%dx(ixO^S,3)) / abs(w(ixO^S,gcak3_)) )
      dt_cakt      = courantpar * minval(tdumt(ixO^S))
      }
    endif
    }

    dtnew = min(dtnew,dt_cakr,dt_cakt,dt_cakp)

  end subroutine cak_get_dt

  !> Compute velocity gradient in direction 'idir' on a non-uniform grid
  subroutine get_velocity_gradient(ixI^L,ixO^L,vfield,x,idir,grad_vn)
    use mod_global_parameters

    integer, intent(in)  :: ixI^L, ixO^L, idir
    real(8), intent(in)  :: vfield(ixI^S), x(ixI^S,1:ndim)
    real(8), intent(out) :: grad_vn(ixO^S)

    ! Local variables
    real(8) :: forw(ixO^S), backw(ixO^S), cent(ixO^S)
    integer :: jrx^L, hrx^L{^NOONED,jtx^L, htx^L}{^IFTHREED,jpx^L, hpx^L}

    ! Index +1 (j) and index -1 (h) in radial direction; kr(dir,dim)=1, dir=dim
    jrx^L=ixO^L+kr(1,^D);
    hrx^L=ixO^L-kr(1,^D);

    {^NOONED
    ! Index +1 (j) and index -1 (h) in polar direction
    jtx^L=ixO^L+kr(2,^D);
    htx^L=ixO^L-kr(2,^D);
    }

    {^IFTHREED
    ! Index +1 (j) and index -1 (h) in azimuthal direction
    jpx^L=ixO^L+kr(3,^D);
    hpx^L=ixO^L-kr(3,^D);
    }

    ! grad(v.n) on non-uniform grid according to Sundqvist & Veronis (1970)
    select case (idir)
    case(1) ! Radial forward, backward, and central derivatives
      forw(ixO^S)  = (x(ixO^S,1) - x(hrx^S,1)) * vfield(jrx^S) &
                      / ((x(jrx^S,1) - x(ixO^S,1)) * (x(jrx^S,1) - x(hrx^S,1)))

      backw(ixO^S) = -(x(jrx^S,1) - x(ixO^S,1)) * vfield(hrx^S) &
                      / ((x(ixO^S,1) - x(hrx^S,1)) * (x(jrx^S,1) - x(hrx^S,1)))

      cent(ixO^S)  = (x(jrx^S,1) + x(hrx^S,1) - 2.0d0*x(ixO^S,1)) * vfield(ixO^S) &
                      / ((x(ixO^S,1) - x(hrx^S,1)) * (x(jrx^S,1) - x(ixO^S,1)))
    {^NOONED
    case(2) ! Polar forward, backward, and central derivatives
      forw(ixO^S)  = (x(ixO^S,2) - x(htx^S,2)) * vfield(jtx^S) &
                      / (x(ixO^S,1) * (x(jtx^S,2) - x(ixO^S,2)) * (x(jtx^S,2) - x(htx^S,2)))

      backw(ixO^S) = -(x(jtx^S,2) - x(ixO^S,2)) * vfield(htx^S) &
                      / ( x(ixO^S,1) * (x(ixO^S,2) - x(htx^S,2)) * (x(jtx^S,2) - x(htx^S,2)))

      cent(ixO^S)  = (x(jtx^S,2) + x(htx^S,2) - 2.0d0*x(ixO^S,2)) * vfield(ixO^S) &
                      / ( x(ixO^S,1) * (x(ixO^S,2) - x(htx^S,2)) * (x(jtx^S,2) - x(ixO^S,2)))
    }
    {^IFTHREED
    case(3) ! Azimuthal forward, backward, and central derivatives
      forw(ixO^S)  = (x(ixO^S,3) - x(hpx^S,3)) *  vfield(jpx^S) &
                      / ( x(ixO^S,1)*sin(x(ixO^S,2)) * (x(jpx^S,3) - x(ixO^S,3)) * (x(jpx^S,3) - x(hpx^S,3)))

      backw(ixO^S) = -(x(jpx^S,3) - x(ixO^S,3)) *  vfield(hpx^S) &
                      / ( x(ixO^S,1)*sin(x(ixO^S,2)) * (x(ixO^S,3) - x(hpx^S,3)) * (x(jpx^S,3) - x(hpx^S,3)))

      cent(ixO^S)  = (x(jpx^S,3) + x(hpx^S,3) - 2.0d0*x(ixO^S,3)) *  vfield(ixO^S) &
                      / ( x(ixO^S,1)*sin(x(ixO^S,2)) * (x(ixO^S,3) - x(hpx^S,3)) * (x(jpx^S,3) - x(ixO^S,3)))
    }
    end select

    ! Total gradient for given velocity field
    grad_vn(ixO^S) = backw(ixO^S) + cent(ixO^S) + forw(ixO^S)

  end subroutine get_velocity_gradient

  !> Initialise (theta',phi') radiation angles coming from stellar disc
  subroutine rays_init(ntheta_point,nphi_point)
    use mod_global_parameters

    ! Subroutine arguments
    integer, intent(in) :: ntheta_point, nphi_point

    ! Local variables
    real(8) :: ymin, ymax, phipmin, phipmax, adum
    integer :: ii

    ! Minimum and maximum range of theta and phi rays
    ! NOTE: theta points are cast into y-space
    ymin    = 0.0d0
    ymax    = 1.0d0
    phipmin = -dpi !0.0d0
    phipmax = dpi !2.0d0*dpi
    ! dphi    = (phipmax - phipmin) / nphi_point

    if (mype == 0) then
      allocate(ay(ntheta_point))
      allocate(wy(ntheta_point))
      allocate(aphi(nphi_point))
      allocate(wphi(nphi_point))

      ! theta and phi ray positions and weights: Gauss-Legendre
      call gauss_legendre_quadrature(ymin,ymax,ntheta_point,ay,wy)
      call gauss_legendre_quadrature(phipmin,phipmax,nphi_point,aphi,wphi)

      ! theta rays and weights: uniform
      ! dth = 1.0d0 / nthetap
      ! adum = ymin + 0.5d0*dth
      ! do ip = 1,nthetap
      !   ay(ip) = adum
      !   wy(ip) = 1.0d0/nthetap
      !   adum = adum + dth
      !   !print*,'phipoints'
      !   !print*,ip,aphi(ip),wphi(ip),dphi
      ! enddo

      ! phi ray position and weights: uniform
      ! adum = phipmin + 0.5d0*dphi
      ! do ii = 1,nphi_point
      !   aphi(ii) = adum
      !   wphi(ii) = 1.0d0/nphi_point
      !   adum     = adum + dphi
      ! enddo

      print*, '==========================='
      print*, '    Radiation ray setup    '
      print*, '==========================='
      print*, 'Theta ray points + weights '
      do ii = 1,ntheta_point
        print*,ii,ay(ii),wy(ii)
      enddo
      print*
      print*, 'Phi ray points + weights   '
      do ii = 1,nphi_point
        print*,ii,aphi(ii),wphi(ii)
      enddo
      print*
    endif

    call MPI_BARRIER(icomm,ierrmpi)

    !===========================
    ! Broadcast what mype=0 read
    !===========================
    if (npe > 1) then
      call MPI_BCAST(ntheta_point,1,MPI_INTEGER,0,icomm,ierrmpi)
      call MPI_BCAST(nphi_point,1,MPI_INTEGER,0,icomm,ierrmpi)

      if (mype /= 0) then
        allocate(ay(ntheta_point))
        allocate(wy(ntheta_point))
        allocate(aphi(nphi_point))
        allocate(wphi(nphi_point))
      endif

      call MPI_BCAST(ay,ntheta_point,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
      call MPI_BCAST(wy,ntheta_point,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
      call MPI_BCAST(aphi,nphi_point,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
      call MPI_BCAST(wphi,nphi_point,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
    endif

  end subroutine rays_init
  
  !> Fast Gauss-Legendre N-point quadrature algorithm by G. Rybicki
  subroutine gauss_legendre_quadrature(xlow,xhi,n,x,w)
    ! Given the lower and upper limits of integration xlow and xhi, and given n,
    ! this routine returns arrays x and w of length n, containing the abscissas
    ! and weights of the Gauss-Legendre N-point quadrature
    use mod_global_parameters

    ! Subroutine arguments
    real(8), intent(in)  :: xlow, xhi
    integer, intent(in)  :: n
    real(8), intent(out) :: x(n), w(n)

    ! Local variables
    integer :: i, j, m
    real(8) :: p1, p2, p3, pp, xl, xm, z, z1
    real(8), parameter :: error=3.0d-14

    m = (n + 1)/2
    xm = 0.5d0*(xhi + xlow)
    xl = 0.5d0*(xhi - xlow)

    do i = 1,m
      z = cos( dpi * (i - 0.25d0)/(n + 0.5d0) )
      z1 = 2.0d0 * z

      do while (abs(z1 - z) > error)
        p1 = 1.0d0
        p2 = 0.0d0

        do j = 1,n
          p3 = p2
          p2 = p1
          p1 = ( (2.0d0*j - 1.0d0)*z*p2 - (j - 1.0d0)*p3 )/j
        enddo

        pp = n*(z*p1 - p2) / (z*z - 1.0d0)
        z1 = z
        z = z1 - p1/pp
      enddo

      x(i)     = xm - xl*z
      x(n+1-i) = xm + xl*z
      w(i)     = 2.0d0*xl / ((1.0d0 - z*z) * pp*pp)
      w(n+1-i) = w(i)
    enddo

  end subroutine gauss_legendre_quadrature

end module mod_cak_force
