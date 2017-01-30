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

  !> Whether Hall-MHD is used
  logical, public, protected              :: mhd_Hall = .false.

  !> Whether MHD-GLM is used
  logical, public, protected              :: mhd_glm = .false.

  !> TODO: describe and set value
  double precision, public, protected     :: mhd_glm_Cr = -2.0d0

  !> MHD fourth order
  logical, public, protected              :: mhd_4th_order = .false.

  !> Number of tracer species
  integer, public, protected              :: mhd_n_tracer = 0

  !> Index of the density (in the w array)
  integer, public, parameter              :: rho_ = 1

  !> Indices of the momentum density
  integer, allocatable, public, protected :: mom(:)

  !> Index of the energy density (-1 if not present)
  integer, public, protected              :: e_

  !> Index of the gas pressure (-1 if not present) should equal e_
  integer, public, protected              :: p_

  !> Indices of the magnetic field
  integer, allocatable, public, protected :: mag(:)

  !> Indices of the GLM psi
  integer, allocatable, public, protected :: psi_

  !> Indices of the tracers
  integer, allocatable, public, protected :: tracer(:)

  !> The number of flux variables in this module
  integer, public, protected              :: mhd_nwflux

  !> The adiabatic index
  double precision, public, protected     :: mhd_gamma = 5/3.0d0

  !> The adiabatic constant
  double precision, public, protected     :: mhd_adiab = 1.0d0

  !> The MHD resistivity
  double precision, public, protected     :: mhd_eta = 0.0d0

  !> The MHD hyper-resistivity
  double precision, public, protected     :: mhd_eta_hyper = 0.0d0

  !> TODO: what is this?
  double precision, public, protected     :: mhd_etah = 0.0d0

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

  ! Public methods
  public :: mhd_phys_init
  public :: mhd_kin_en
  public :: mhd_get_pthermal
  public :: mhd_get_v
  public :: mhd_to_conserved
  public :: mhd_to_primitive
  public :: mhd_get_csound2
  public :: get_divb

contains

  !> Read this module"s parameters from a file
  subroutine mhd_read_params(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /mhd_list/ mhd_energy, mhd_n_tracer, mhd_gamma, mhd_adiab,&
      mhd_eta, mhd_eta_hyper, mhd_etah, mhd_glm, mhd_glm_Cr, &
      mhd_thermal_conduction, mhd_radiative_cooling, mhd_Hall, &
      mhd_4th_order, typedivbfix, divbdiff, typedivbdiff, compactres, divbwave

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, mhd_list, end=111)
111    close(unitpar)
    end do

  end subroutine mhd_read_params

  subroutine mhd_phys_init()
    use mod_global_parameters
    use mod_thermal_conduction
    use mod_radiative_cooling
    use mod_physics

    integer :: itr, idir

    call mhd_read_params(par_files)

    physics_type = "mhd"

    ! Determine flux variables
    nwflux = 1                  ! rho (density)

    allocate(mom(ndir))
    do idir = 1, ndir
       nwflux    = nwflux + 1
       mom(idir) = nwflux       ! momentum density
    end do

    ! Set index of energy variable
    if (mhd_energy) then
       nwwave=8
       nwflux = nwflux + 1
       e_     = nwflux          ! energy density
       p_     = nwflux          ! gas pressure
    else
       nwwave=7
       e_ = -1
       p_ = -1
    end if

    allocate(mag(ndir))
    do idir = 1, ndir
       nwflux    = nwflux + 1
       mag(idir) = nwflux       ! magnetic field
    end do

    if (mhd_glm) then
       nwflux = nwflux + 1
       psi_   = nwflux
    else
       psi_ = -1
    end if

    allocate(tracer(mhd_n_tracer))

    ! Set starting index of tracers
    do itr = 1, mhd_n_tracer
       nwflux = nwflux + 1
       tracer(itr) = nwflux     ! tracers
    end do

    mhd_nwflux = nwflux

    nwaux   = 0
    nwextra = 0
    nw      = nwflux + nwaux + nwextra
    nflag_  = nw + 1

    nvector      = 2 ! No. vector vars
    allocate(iw_vector(nvector))
    iw_vector(1) = mom(1) - 1   ! TODO: why like this?
    iw_vector(2) = mag(1) - 1   ! TODO: why like this?

    phys_get_dt          => mhd_get_dt
    phys_get_cmax        => mhd_get_cmax
    phys_get_flux        => mhd_get_flux
    phys_add_source_geom => mhd_add_source_geom
    phys_to_conserved    => mhd_to_conserved
    phys_to_primitive    => mhd_to_primitive
    phys_check_params    => mhd_check_params
    phys_check_w         => mhd_check_w
    phys_get_pthermal    => mhd_get_pthermal

    ! initialize thermal conduction module
    if (mhd_thermal_conduction) then
       call thermal_conduction_init(mhd_gamma)
    end if

    ! Initialize radiative cooling module
    if (mhd_radiative_cooling) then
      call radiative_cooling_init(mhd_gamma)
    end if

    if (mhd_glm) then
       ! Solve the Riemann problem for the linear 2x2 system for normal
       ! B-field and GLM_Psi according to Dedner 2002:
       phys_modify_wLR => glmSolve
    end if

    ! For Hall, we need one more reconstructed layer since currents are computed
    ! in getflux: assuming one additional ghost layer (two for FOURTHORDER) was
    ! added in dixB.
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

    minrho = max(0.0d0, smallrho)

    if (.not. mhd_energy) then
       if (mhd_gamma <= 0.0d0) call mpistop ("Error: mhd_gamma <= 0")
       if (mhd_adiab < 0.0d0) call mpistop ("Error: mhd_adiab < 0")
       minp   = mhd_adiab*minrho**mhd_gamma
    else
       if (mhd_gamma <= 0.0d0 .or. mhd_gamma == 1.0d0) &
            call mpistop ("Error: mhd_gamma <= 0 or mhd_gamma == 1")
       minp   = max(0.0d0, smallp)
       smalle = minp/(mhd_gamma - 1.0d0)
    end if

  end subroutine mhd_check_params

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

    ! if (fixsmall) call smallvalues(w,x,ixI^L,ixO^L,"conserve")
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

    ! call handle_small_values(.true., w, x, ixI^L, ixO^L)
  end subroutine mhd_to_primitive

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

    call mhd_get_p_total(w,x,ixI^L,ixO^L,ptotal)

    ! Get flux of density
    f(ixO^S,rho_)=v(ixO^S,idim)*w(ixO^S,rho_)

    ! Get flux of momentum
    ! f_i[m_k]=v_i*m_k-b_k*b_i [+ptotal if i==k]
    do idir=1,ndir
      if(idim==idir) then
        f(ixO^S,mom(idir))=ptotal(ixO^S)-w(ixO^S,mag(idim))*w(ixO^S,mag(idir))
      else
        f(ixO^S,mom(idir))= -w(ixO^S,mag(idir))*w(ixO^S,mag(idim))
      end if
      if (B0field) then
        f(ixO^S,mom(idir))=f(ixO^S,mom(idir))&
             -w(ixO^S,mag(idir))*myB0%w(ixO^S,idim)&
             -w(ixO^S,mag(idim))*myB0%w(ixO^S,idir)
      end if
      f(ixO^S,mom(idir))=v(ixO^S,idim)*w(ixO^S,mom(idir))
    end do

    ! Get flux of energy
    ! f_i[e]=v_i*e+v_i*ptotal-b_i*(b_k*v_k)
    if(mhd_energy) then
      f(ixO^S,e_)=v(ixO^S,idim)*(w(ixO^S,e_)+ptotal(ixO^S))- &
            w(ixO^S,mag(idim))*sum(w(ixO^S,mag(:))*v(ixO^S,:),dim=ndim+1)

      if (B0field) then
        tmp(ixO^S)=sum(myB0%w(ixO^S,:)*w(ixO^S,mag(:)),dim=ndim+1)
        f(ixO^S,e_) = f(ixO^S,e_) &
             + v(ixO^S,idim) * tmp(ixO^S) &
             - sum(v(ixO^S,:)*w(ixO^S,mag(:))**2,dim=ndim+1) * myB0%w(ixO^S,idim)
      end if

      if (mhd_Hall) then
        ! f_i[e]= f_i[e] + vHall_i*(b_k*b_k) - b_i*(vHall_k*b_k)
        if (mhd_etah>zero) then
          f(ixO^S,e_) = f(ixO^S,e_) + vHall(ixO^S,idim) * &
                 sum(w(ixO^S, mag(:))**2,dim=ndim+1) &
               - w(ixO^S,mag(idim)) * sum(vHall(ixO^S,:)*w(ixO^S,mag(:))**2,dim=ndim+1)
          if (B0field) then
            f(ixO^S,e_) = f(ixO^S,e_) &
                 + vHall(ixO^S,idim) * tmp(ixO^S) &
                 - sum(vHall(ixO^S,:)*w(ixO^S,mag(:))**2,dim=ndim+1) * myB0%w(ixO^S,idim)
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
                +v(ixO^S,idim)*myB0%w(ixO^S,idir)&
                -v(ixO^S,idir)*myB0%w(ixO^S,idim)
        end if

        if (mhd_Hall) then
          ! f_i[b_k] = f_i[b_k] + vHall_i*b_k - vHall_k*b_i
          if (mhd_etah>zero) then
            if (B0field) then
              f(ixO^S,mag(idir)) = f(ixO^S,mag(idir)) &
                   - vHall(ixO^S,idir)*(w(ixO^S,mag(idim))+myB0%w(ixO^S,idim)) &
                   + vHall(ixO^S,idim)*(w(ixO^S,mag(idir))+myB0%w(ixO^S,idir))
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
  subroutine mhd_add_source(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,qsourcesplit)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L, iw^LIM
    double precision, intent(in)    :: qdt, qtC, qt
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

    if (mhd_radiative_cooling) then
       call radiative_cooling_add_source(qdt,ixI^L,ixO^L,wCT,w,x,qsourcesplit)
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
       case ('powel')
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
       end select
    end if
    }
  end subroutine mhd_add_source

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

       if (fixsmall) call smallvalues(w,x,ixI^L,ixO^L,"add_source_res1")
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

       if (fixsmall) call smallvalues(w,x,ixI^L,ixO^L,"add_source_res2")
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
    !-----------------------------------------------------------------------------
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

       if (fixsmall) call smallvalues(w,x,ixI^L,ixO^L,"add_source_hyperres")
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
    ! since this option changes energy: smallvalues call
    if (fixsmall) call smallvalues(w,x,ixI^L,ixO^L,"add_source_glm")
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
    if (mhd_energy) then
    ! since this option changes energy: smallvalues call
    if (fixsmall) call smallvalues(w,x,ixI^L,ixO^L,"add_source_glm2")
    end if
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

    ! since this option changes energy: smallvalues call
    if (fixsmall) call smallvalues(w,x,ixI^L,ixO^L,"add_source_powel")

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
                          /(^D&1.0d0/mygeo%dx(ixp^S,^D)**2+)
       end if

       do idir=1,ndir
         w(ixp^S,mag(idir))=w(ixp^S,mag(idir))+graddivb(ixp^S)
       end do

       if (mhd_energy .and. typedivbdiff=='all') then
         ! e += B_idim*eta*grad_idim(divb)
         w(ixp^S,e_)=w(ixp^S,e_)+wCT(ixp^S,mag(idim))*graddivb(ixp^S)
       end if
    end do

    if (fixsmall) call smallvalues(w,x,ixI^L,ixp^L,"add_source_linde")

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
          bvec(ixI^S,idir)=w(ixI^S,mag(idir))+myB0_cell%w(ixI^S,idir)
       end do
    else
       do idir = 1, ndir
          bvec(ixI^S,idir)=w(ixI^S,mag(idir))
       end do
    end if

    !bvec(ixI^S,1:ndir)=w(ixI^S,b0_+1:b0_+ndir)
    call curlvector(bvec,ixI^L,ixO^L,current,idirmin,idirmin0,ndir)

  end subroutine get_current

  !> If resistivity is not zero, check diffusion time limit for dt
  subroutine mhd_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
    use mod_global_parameters
    use mod_usr_methods
    use mod_radiative_cooling, only: cooling_get_dt

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

  end subroutine mhd_get_dt

  ! Add geometrical source terms to w
  subroutine mhd_add_source_geom(qdt,ixI^L,ixO^L,wCT,w,x)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: wCT(ixI^S,1:nw), w(ixI^S,1:nw)

    integer          :: iw,idir, h1x^L{^NOONED, h2x^L}
    double precision :: tmp(ixI^S)
    logical          :: angmomfix=.false.

    ! TODO
    ! INTEGER,PARAMETER:: mr_=m0_+r_,mphi_=m0_+phi_,mz_=m0_+z_  ! Polar var. names
    ! integer,parameter:: br_=b0_+r_,bphi_=b0_+phi_,bz_=b0_+z_


    ! select case (typeaxial)
    ! case ('slab')
    !    ! No source terms in slab symmetry
    ! case ('cylindrical')
    !    do iw=1,nwflux
    !       select case (iw)
    !          ! s[mr]=(ptotal-Bphi**2+mphi**2/rho)/radius
    !       case (mr_)
    !          call mhd_get_p_total(wCT,x,ixI^L,ixO^L,tmp)
    !          w(ixO^S,iw)=w(ixO^S,iw)+qdt*tmp(ixO^S)/x(ixO^S,1)
    !          tmp(ixO^S)=zero
    !          {^IFPHI
    !          tmp(ixO^S)= &
    !               -wCT(ixO^S,bphi_)**2+wCT(ixO^S,mphi_)**2/wCT(ixO^S,rho_)

    !          ! s[mphi]=(-mphi*mr/rho+Bphi*Br)/radius
    !       case (mphi_)
    !          tmp(ixO^S)= &
    !               -wCT(ixO^S,mphi_)*wCT(ixO^S,mr_)/wCT(ixO^S,rho_) &
    !               +wCT(ixO^S,bphi_)*wCT(ixO^S,br_)

    !          ! s[Bphi]=((Bphi*mr-Br*mphi)/rho)/radius
    !       case (bphi_)
    !          tmp(ixO^S)=(wCT(ixO^S,bphi_)*wCT(ixO^S,mr_) &
    !               -wCT(ixO^S,br_)*wCT(ixO^S,mphi_)) &
    !               /wCT(ixO^S,rho_)
    !          }
    !          {#IFDEF GLM      ! s[br]=psi/radius
    !       case (br_)
    !          tmp(ixO^S)=wCT(ixO^S,psi_)
    !          }
    !       end select

    !       ! Divide by radius and add to w
    !       if (iw==mr_{#IFDEF GLM .or.iw==br_}{^IFPHI .or.iw==mphi_.or.iw==bphi_}) then
    !          w(ixO^S,iw)=w(ixO^S,iw)+qdt*tmp(ixO^S)/x(ixO^S,1)
    !       end if

    !    end do
    ! case ('spherical')
    !    h1x^L=ixO^L-kr(1,^D); {^NOONED h2x^L=ixO^L-kr(2,^D);}
    !    do iw=1,nwflux
    !       select case (iw)
    !          ! s[m1]=((mtheta**2+mphi**2)/rho+2*ptotal-(Btheta**2+Bphi**2))/r
    !       case (m1_)
    !          call mhd_get_p_total(wCT,x,ixI^L,ixO^L,tmp)
    !          if (B0field) then
    !             tmp(ixO^S)=tmp(ixO^S)+{^C&myB0_cell%w(ixO^S,^C)*wCT(ixO^S,mag(idir))+}
    !          end if
    !          ! For nonuniform Cartesian grid this provides hydrostatic equil.
    !          tmp(ixO^S)=tmp(ixO^S)*x(ixO^S,1) &
    !               *(mygeo%surfaceC1(ixO^S)-mygeo%surfaceC1(h1x^S)) &
    !               /mygeo%dvolume(ixO^S){&^CE&
    !               +wCT(ixO^S,m^CE_)**2/wCT(ixO^S,rho_)-wCT(ixO^S,b^CE_)**2 }
    !          if (B0field.and.ndir>1) then
    !             tmp(ixO^S)=tmp(ixO^S){^CE&-2.0d0*myB0_cell%w(ixO^S,^CE) &
    !                  *wCT(ixO^S,b^CE_)|}
    !          end if
    !          {^NOONEC
    !          ! s[m2]=-(mr*mtheta/rho-Br*Btheta)/r
    !          !       + cot(theta)*(mphi**2/rho+(p+0.5*B**2)-Bphi**2)/r
    !       case (m2_)
    !          }
    !          {^NOONED
    !          call mhd_get_p_total(wCT,x,ixI^L,ixO^L,tmp)
    !          if (B0field) then
    !             tmp(ixO^S)=tmp(ixO^S)+{^C&myB0_cell%w(ixO^S,^C)*wCT(ixO^S,mag(idir))+}
    !          end if
    !          ! This will make hydrostatic p=const an exact solution
    !          w(ixO^S,iw)=w(ixO^S,iw)+qdt*tmp(ixO^S) &
    !               *(mygeo%surfaceC2(ixO^S)-mygeo%surfaceC2(h2x^S)) &
    !               /mygeo%dvolume(ixO^S)
    !          }
    !          {^NOONEC
    !          tmp(ixO^S)=-(wCT(ixO^S,m1_)*wCT(ixO^S,m2_)/wCT(ixO^S,rho_) &
    !               -wCT(ixO^S,b1_)*wCT(ixO^S,b2_))
    !          if (B0field) then
    !             tmp(ixO^S)=tmp(ixO^S)+myB0_cell%w(ixO^S,1)*wCT(ixO^S,b2_) &
    !                  +wCT(ixO^S,b1_)*myB0_cell%w(ixO^S,2)
    !          end if
    !          }
    !          {^IFTHREEC
    !          {^NOONED
    !          tmp(ixO^S)=tmp(ixO^S)+(wCT(ixO^S,m3_)**2/wCT(ixO^S,rho_) &
    !               -wCT(ixO^S,b3_)**2)*dcos(x(ixO^S,2)) &
    !               /dsin(x(ixO^S,2))
    !          if (B0field) then
    !             tmp(ixO^S)=tmp(ixO^S)-2.0d0*myB0_cell%w(ixO^S,3)*wCT(ixO^S,b3_)&
    !                  *dcos(x(ixO^S,2))/dsin(x(ixO^S,2))
    !          end if
    !          }
    !          ! s[m3]=-(mphi*mr/rho-Bphi*Br)/r
    !          !       -cot(theta)*(mtheta*mphi/rho-Btheta*Bphi)/r
    !       case (m3_)
    !          if (.not.angmomfix) then
    !             tmp(ixO^S)=-(wCT(ixO^S,m3_)*wCT(ixO^S,m1_)/wCT(ixO^S,rho_) &
    !                  -wCT(ixO^S,b3_)*wCT(ixO^S,b1_)) {^NOONED &
    !                  -(wCT(ixO^S,m2_)*wCT(ixO^S,m3_)/wCT(ixO^S,rho_) &
    !                  -wCT(ixO^S,b2_)*wCT(ixO^S,b3_)) &
    !                  *dcos(x(ixO^S,2))/dsin(x(ixO^S,2)) }
    !             if (B0field) then
    !                tmp(ixO^S)=tmp(ixO^S)+myB0_cell%w(ixO^S,1)*wCT(ixO^S,b3_) &
    !                     +wCT(ixO^S,b1_)*myB0_cell%w(ixO^S,3) {^NOONED &
    !                     +(myB0_cell%w(ixO^S,2)*wCT(ixO^S,b3_) &
    !                     +wCT(ixO^S,b2_)*myB0_cell%w(ixO^S,3)) &
    !                     *dcos(x(ixO^S,2))/dsin(x(ixO^S,2)) }
    !             end if
    !          end if
    !          }
    !          {#IFDEF GLM
    !          ! s[b1]=2*psi/r
    !       case (b1_)
    !          tmp(ixO^S)=2.0d0*wCT(ixO^S,psi_)
    !          }
    !          {^NOONEC
    !          ! s[b2]=(mr*Btheta-mtheta*Br)/rho/r
    !          !       + cot(theta)*psi/r
    !       case (b2_)
    !          tmp(ixO^S)=(wCT(ixO^S,m1_)*wCT(ixO^S,b2_) &
    !               -wCT(ixO^S,m2_)*wCT(ixO^S,b1_))/wCT(ixO^S,rho_)
    !          if (B0field) then
    !             tmp(ixO^S)=tmp(ixO^S)+(wCT(ixO^S,m1_)*myB0_cell%w(ixO^S,2) &
    !                  -wCT(ixO^S,m2_)*myB0_cell%w(ixO^S,1))/wCT(ixO^S,rho_)
    !          end if
    !          {#IFDEF GLM
    !          tmp(ixO^S)=tmp(ixO^S) &
    !               + dcos(x(ixO^S,2))/dsin(x(ixO^S,2))*wCT(ixO^S,psi_)
    !          }
    !          }
    !          {^IFTHREEC
    !          ! s[b3]=(mr*Bphi-mphi*Br)/rho/r
    !          !       -cot(theta)*(mphi*Btheta-mtheta*Bphi)/rho/r
    !       case (b3_)
    !          tmp(ixO^S)=(wCT(ixO^S,m1_)*wCT(ixO^S,b3_) &
    !               -wCT(ixO^S,m3_)*wCT(ixO^S,b1_))/wCT(ixO^S,rho_) {^NOONED &
    !               -(wCT(ixO^S,m3_)*wCT(ixO^S,b2_) &
    !               -wCT(ixO^S,m2_)*wCT(ixO^S,b3_))*dcos(x(ixO^S,2)) &
    !               /(wCT(ixO^S,rho_)*dsin(x(ixO^S,2))) }
    !          if (B0field) then
    !             tmp(ixO^S)=tmp(ixO^S)+(wCT(ixO^S,m1_)*myB0_cell%w(ixO^S,3) &
    !                  -wCT(ixO^S,m3_)*myB0_cell%w(ixO^S,1))/wCT(ixO^S,rho_){^NOONED &
    !                  -(wCT(ixO^S,m3_)*myB0_cell%w(ixO^S,2) &
    !                  -wCT(ixO^S,m2_)*myB0_cell%w(ixO^S,3))*dcos(x(ixO^S,2)) &
    !                  /(wCT(ixO^S,rho_)*dsin(x(ixO^S,2))) }
    !          end if
    !          }
    !       end select
    !       ! Divide by radius and add to w
    !       if (iw==m1_{#IFDEF GLM .or.iw==b1_}{^NOONEC.or.iw==m2_.or.iw==b2_}&
    !            {^IFTHREEC .or.iw==b3_ .or.(iw==m3_.and..not.angmomfix)}) &
    !            w(ixO^S,iw)=w(ixO^S,iw)+qdt*tmp(ixO^S)/x(ixO^S,1)
    !    end do
    ! end select
  end subroutine mhd_add_source_geom

  !> Compute 2 times total magnetic energy
  function mhd_mag_en_all(w, ixI^L, ixO^L) result(mge)
    use mod_global_parameters
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S, nw)
    double precision              :: mge(ixO^S)

    if (B0field) then
      mge = sum((w(ixO^S, mag(:))+myB0%w(ixO^S,:))**2, dim=ndim+1)
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
      mgf = w(ixO^S, mag(idir))+myB0%w(ixO^S,idir)
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
       bmag(ixO^S)=sqrt(sum((w(ixO^S,mag(:)) + myB0%w(ixO^S,1:ndir))**2))
    end if

    dthall=dtdiffpar*minval(dxarr(1:ndim))**2.0d0/(mhd_etah*maxval(bmag(ixO^S)/w(ixO^S,rho_)))
  end subroutine mhd_getdt_Hall
  ! TODO: rewrite
  subroutine smallvalues(w,x,ixI^L,ixO^L,subname)

    use mod_global_parameters

    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    character(len=*), intent(in)    ::subname
    !.. local ..
    integer                         :: posvec(ndim)
    integer, dimension(ixI^S)       :: patchierror
    double precision                :: pth(ixI^S), Te(ixI^S)
    !-----------------------------------------------------------------------------

    ! {#IFDEF ENERGY
    ! pth(ixO^S)=(mhd_gamma-one)*(w(ixO^S,e_)- &
    !      half*(({^C&w(ixO^S,m^C_)**2+})/w(ixO^S,rho_)&
    !      +{ ^C&w(ixO^S,mag(idir))**2+}))
    ! if (smallT>0.d0) then
    !    Te(ixO^S)=pth(ixO^S)/w(ixO^S,rho_)
    ! else
    !    Te(ixO^S)=zero
    ! end if

    ! if (strictsmall) then
    !    if (smallT>0.d0 .and. any(Te(ixO^S) <=smallT)) then
    !       print *,'SMALLVALUES of temperature under strictsmall problem From:  ', &
    !            subname,' iteration=', it,' time=',t
    !       posvec(1:ndim)=minloc(Te(ixO^S))
    !       ^D&posvec(^D)=posvec(^D)+ixOmin^D-1;
    !       write(*,*)'minimum temperature= ', minval(Te(ixO^S)),' with limit=',smallT,&
    !            ' at x=', x(^D&posvec(^D),1:ndim),' array index=',posvec,' where rho=',&
    !            w({^D&posvec(^D)},rho_),', velocity v=',&
    !            ^C&w({^D&posvec(^D)},m^C_)/w({^D&posvec(^D)},rho_),&
    !            ', and magnetic field B=',^C&w({^D&posvec(^D)},mag(idir)),&
    !            ' w(1:nwflux)=',w(^D&posvec(^D),1:nwflux)
    !       call mpistop("Smallvalues of temperature with strictsmall=T failed")
    !    end if
    !    if (any(pth(ixO^S) <=minp)) then
    !       print *,'SMALLVALUES of pressure under strictsmall problem From:  ', &
    !            subname,' iteration=', it,' time=',t
    !       posvec(1:ndim)=minloc(pth(ixO^S))
    !       ^D&posvec(^D)=posvec(^D)+ixOmin^D-1;
    !       write(*,*)'minimum pressure = ', minval(pth(ixO^S)),' with limit=',minp,&
    !            ' at x=', x(^D&posvec(^D),1:ndim),' array index=',posvec,' where rho=',&
    !            w({^D&posvec(^D)},rho_),', velocity v=',&
    !            ^C&w({^D&posvec(^D)},m^C_)/w({^D&posvec(^D)},rho_),&
    !            ', and magnetic field B=',^C&w({^D&posvec(^D)},mag(idir)),&
    !            ' w(1:nwflux)=',w(^D&posvec(^D),1:nwflux)
    !       call mpistop("Smallvalues of pressure with strictsmall=T failed")
    !    end if
    !    if (any(w(ixO^S,e_) <=smalle)) then
    !       print *,'SMALLVALUES of energy under strictsmall problem From:  ', &
    !            subname,' iteration=', it,' time=',t
    !       posvec(1:ndim)=minloc(w(ixO^S,e_))
    !       ^D&posvec(^D)=posvec(^D)+ixOmin^D-1;
    !       write(*,*)'minimum e =', minval(w(ixO^S,e_)),' with limit=',smalle,&
    !            ' at x=', x(^D&posvec(^D),1:ndim),' array index=',posvec,' where E_k=',&
    !            half*(^C&w(^D&posvec(^D),m^C_)**2+)/w(^D&posvec(^D),rho_),&
    !            ' E_total=',w(^D&posvec(^D),e_),&
    !            ' w(1:nwflux)=',w(^D&posvec(^D),1:nwflux)
    !       call mpistop("Smallvalues of energy with strictsmall=T failed")
    !    end if
    !    if (any(w(ixO^S,rho_) <=minrho)) then
    !       print *,'SMALLVALUES of density under strictsmall problem From:  ', &
    !            subname,' iteration=', it,' time=',t
    !       posvec(1:ndim)=minloc(w(ixO^S,rho_))
    !       ^D&posvec(^D)=posvec(^D)+ixOmin^D-1;
    !       write(*,*)'minimum rho =', minval(w(ixO^S,rho_)),' with limit=',minrho,&
    !            ' at x=', x(^D&posvec(^D),1:ndim),' array index=',posvec,' where E_k=',&
    !            half*(^C&w(^D&posvec(^D),m^C_)**2+)/w(^D&posvec(^D),rho_),&
    !            ' E_total=',w(^D&posvec(^D),e_),&
    !            ' w(1:nwflux)=',w(^D&posvec(^D),1:nwflux)
    !       call mpistop("Smallvalues of density with strictsmall=T failed")
    !    end if
    ! else
    !    if (strictgetaux)then
    !       where(w(ixO^S,rho_) < minrho)
    !          w(ixO^S,rho_)=minrho
    !          {^C&w(ixO^S,m^C_) =zero;}
    !       end where
    !       where(pth(ixO^S) < minp)
    !          w(ixO^S,e_)=minp/(mhd_gamma-one)+&
    !               (({^C&w(ixO^S,m^C_)**2+})/w(ixO^S,rho_)+{^C&w(ixO^S,mag(idir))**2+})*half
    !       end where
    !       where(Te(ixO^S) < smallT)
    !          w(ixO^S,e_)=smallT*w(ixO^S,rho_)/(mhd_gamma-one)+&
    !               (({^C&w(ixO^S,m^C_)**2+})/w(ixO^S,rho_)+{^C&w(ixO^S,mag(idir))**2+})*half
    !       end where
    !    else
    !       where(w(ixO^S,rho_) < minrho .or. w(ixO^S,e_) < smalle&
    !            .or. pth(ixO^S) < minp .or. Te(ixO^S) < smallT)
    !          patchierror(ixO^S)=-1
    !       elsewhere
    !          patchierror(ixO^S)=0
    !       end where
    !       call correctaux(ixI^L,ixO^L,w,x,patchierror,subname)
    !    end if
    ! end if
    ! }

    ! {#IFDEF ISO
    ! if (any(w(ixO^S,rho_) < minrho)) then
    !    if (strictsmall)then
    !       write(*,*)'SMALLVALUES of density under strictsmall problem From:  ', &
    !            subname,' iteration ', it,' time ',t
    !       posvec(1:ndim)=minloc(w(ixO^S,rho_))
    !       ^D&posvec(^D)=posvec(^D)+ixOmin^D-1;
    !       write(*,*)'minimum rho =', minval(w(ixO^S,rho_)),' with limit=',minrho,&
    !            ' at x=', x(^D&posvec(^D),1:ndim),' array index=',posvec,' where E_k=',&
    !            half*(^C&w(^D&posvec(^D),m^C_)**2+)/w(^D&posvec(^D),rho_),&
    !            ' w(1:nwflux)=',w(^D&posvec(^D),1:nwflux)
    !       call mpistop("Smallvalues of density with strictsmall=T failed")
    !    else
    !       if (strictgetaux)then
    !          where(w(ixO^S,rho_) < minrho)
    !             w(ixO^S,rho_)  = 2.0*(1.0d0 + 10.0d0 * minrho)*minrho
    !          end where
    !       else
    !          where(w(ixO^S,rho_) < minrho)
    !             patchierror(ixO^S)=-1
    !          elsewhere
    !             patchierror(ixO^S)=0
    !          end where
    !          call correctaux(ixI^L,ixO^L,w,x,patchierror,subname)
    !       end if
    !    end if ! strict
    ! end if
    ! }
  end subroutine smallvalues

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

end module mod_mhd_phys
