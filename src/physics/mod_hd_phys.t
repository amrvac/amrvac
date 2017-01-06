! TODO:
! * Keep/remove strictgetaux?
! * Can we make methods robust without fixes?
! * Generic names for momentum and other indices?
! * Remove "subname"
! * Use check_w more generally?
! * Single method to quit and show problematic values

!> Hydrodynamics module
module mod_hd_phys
  use mod_physics

  implicit none
  private

  !> Whether an energy equation is used
  logical, public, protected              :: hd_energy = .true.

  !> Number of tracer species
  integer, public, protected              :: hd_n_tracer = 0

  !> Index of the density (in the w array)
  integer, public, parameter              :: rho_ = 1

  !> Indices of the momentum density
  integer, allocatable, public, protected :: mom(:)

  !> Indices of the tracers
  integer, allocatable, public, protected :: tracer(:)

  !> Index of the energy density (-1 if not present)
  integer, public, protected              :: e_

  !> The number of flux variables in this module
  integer, public, protected              :: hd_nwflux

  !> The adiabatic index
  double precision, public, protected     :: hd_gamma = 5/3.0d0

  !> The adiabatic constant
  double precision, public, protected     :: hd_adiab = 1.0d0

  !> The smallest allowed energy
  double precision, protected             :: smalle

  !> The smallest allowed density
  double precision, protected             :: minrho

  !> The smallest allowed pressure
  double precision, protected             :: minp

  ! Public methods
  public :: hd_phys_init
  public :: hd_kin_en
  public :: hd_get_pthermal
  public :: hd_get_v
  public :: hd_to_conserved
  public :: hd_to_primitive

contains

  !> Read this module"s parameters from a file
  subroutine hd_params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /hd_list/ hd_energy, hd_n_tracer, hd_gamma, hd_adiab

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, hd_list, end=111)
111    close(unitpar)
    end do

  end subroutine hd_params_read

  !> Initialize the module
  subroutine hd_phys_init()
    use mod_global_parameters

    integer :: itr, idir

    call hd_params_read(par_files)

    physics_type = "hd"

    ! Determine flux variables
    nwflux = 1                  ! rho (density)

    allocate(mom(ndir))
    do idir = 1, ndir
       nwflux    = nwflux + 1
       mom(idir) = nwflux       ! momentum density
    end do

    ! Set index of energy variable
    if (hd_energy) then
       nwflux = nwflux + 1
       e_     = nwflux          ! energy density
    else
       e_ = -1
    end if

    allocate(tracer(hd_n_tracer))

    ! Set starting index of tracers
    do itr = 1, hd_n_tracer
       nwflux = nwflux + 1
       tracer(itr) = nwflux     ! tracers
    end do

    hd_nwflux = nwflux

    ! call dust_init()

    nwaux   = 0
    nwextra = 0
    nw      = nwflux + nwaux + nwextra
    nflag_  = nw + 1

    nvector      = 1 ! No. vector vars
    allocate(iw_vector(nvector))
    iw_vector(1) = mom(1) - 1   ! TODO: why like this?

    phys_get_v           => hd_get_v
    phys_get_dt          => hd_get_dt
    phys_get_cmax        => hd_get_cmax
    phys_get_flux        => hd_get_flux
    phys_add_source_geom => hd_add_source_geom
    phys_to_conserved    => hd_to_conserved
    phys_to_primitive    => hd_to_primitive
    phys_check_params    => hd_check_params
    phys_check_w         => hd_check_w

  end subroutine hd_phys_init

  subroutine hd_check_params
    use mod_global_parameters

    if (.not. hd_energy) then
       if (hd_gamma <= 0.0d0) call mpistop ("hd_gamma negative not ok")
       if (hd_adiab < 0.0d0) call mpistop ("adiab strict negative not ok")
       minrho = max(0.0d0, smallrho)
       minp   = hd_adiab*minrho**hd_gamma
    else
       if (hd_gamma <= 0.0d0 .or. hd_gamma == 1.0d0) &
            call mpistop ("hd_gamma negative or 1 not ok")
       minp   = max(0.0d0, smallp)
       minrho = max(0.0d0, smallrho)
       smalle = minp/(hd_gamma - 1.0d0)
    end if

    ! TODO: dust
    ! if (dust_num_species > 0) then
    !    if (eqpar(mu_)<= 0.0d0) call mpistop ("mu (molecular weight) negative not ok")
    !    minrhod = max(0.0d0, smallrhod)
    ! end if

  end subroutine hd_check_params

  !> @todo can perhaps use this function more generally?
  subroutine hd_check_w(checkprimitive, ixI^L, ixO^L, w, flag)
    use mod_global_parameters

    logical             :: checkprimitive
    integer, intent(in) :: ixI^L, ixO^L
    double precision    :: w(ixI^S, nw)
    logical             :: flag(ixG^T)
    double precision    :: tmp(ixG^T)

    flag(ixG^T) =.true.

    if (hd_energy) then
       if (checkprimitive) then
          flag(ixO^S) =(w(ixO^S, e_)>= minp .and. w(ixO^S, rho_)>= minrho)
       else
          tmp(ixO^S) =(hd_gamma - 1.0d0)*(w(ixO^S, e_) - &
               hd_kin_en(w, ixI^L, ixO^L))
          flag(ixO^S) =(tmp(ixO^S)>= minp .and. w(ixO^S, rho_)>= minrho)
       endif
    else
       if (hd_adiab > 0.0d0) then
          flag(ixO^S) = (w(ixO^S, rho_) >= 0.0d0)
       endif
    end if

  end subroutine hd_check_w

  ! Transform primitive variables into conservative ones
  recursive subroutine hd_to_conserved(ixI^L, ixO^L, w, x, fix)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    logical, intent(in), optional   :: fix
    double precision                :: invgam
    integer                         :: idir, itr
    logical                         :: apply_fixes

    apply_fixes = .true.
    if (present(fix)) apply_fixes = fix

    invgam = 1.d0/(hd_gamma - 1.0d0)

    ! Convert velocity to momentum
    do idir = 1, ndir
       w(ixO^S, mom(idir)) = w(ixO^S, rho_) * w(ixO^S, mom(idir))
    end do

    if (hd_energy) then
       ! Calculate total energy from pressure and kinetic energy
       w(ixO^S, e_) = w(ixO^S, e_) * invgam + hd_kin_en(w, ixI^L, ixO^L)
    end if

    do itr = 1, hd_n_tracer
       w(ixO^S, tracer(itr)) = w(ixO^S, rho_) * w(ixO^S, tracer(itr))
    end do

    ! call dust_conserve(...)

    if (apply_fixes .and. fixsmall) &
         call smallvalues(w, x, ixI^L, ixO^L, "hd_to_conserved")

  end subroutine hd_to_conserved

  !> Transform conservative variables into primitive ones
  recursive subroutine hd_to_primitive(ixI^L, ixO^L, w, x, fix)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    logical, intent(in), optional   :: fix
    integer, dimension(ixG^T)       :: patchierror
    integer, dimension(ndim)        :: lowpindex
    integer                         :: itr, idir
    logical                         :: apply_fixes

    apply_fixes = .true.
    if (present(fix)) apply_fixes = fix

    if (apply_fixes .and. fixsmall) &
         call smallvalues(w, x, ixI^L, ixO^L, "hd_to_primitive")

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

    if (apply_fixes) then
       if (strictsmall) then
          if (any(w(ixO^S, e_)<minp)) then
             call mpistop("=== primitive pressure problem===")
          end if
       else if (strictgetaux) then
          ! TODO: check
          where(w(ixO^S, e_)<minp)
             w(ixO^S, e_) = minp
          endwhere
       else
          where(w(ixO^S, e_)<minp)
             patchierror(ixO^S) = 1
          else where
             patchierror(ixO^S) = 0
          end where

          if (any(patchierror(ixO^S)/= 0)) &
               call correctaux(ixI^L, ixO^L, w, x, patchierror, "hd_to_primitive")
       end if
    end if

    ! Convert dust momentum to dust velocity
    ! call dust_primitive(...)
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
    double precision, intent(in)  :: w(ixI^S, nw), x(ixI^S, 1:^ND)
    double precision, intent(out) :: v(ixG^T)

    v(ixO^S) = w(ixO^S, mom(idim)) * hd_inv_rho(w, ixI^L, ixO^L)

    ! Jannis: Remove this case?
    ! if (hd_energy .or. hd_adiab > 0.0d0) then
    ! else
    !    ! case of zero temperature: allow zero density
    !    where(w(ixO^S, rho_)/= zero)
    !       v(ixO^S) = w(ixO^S, mom(idims))/w(ixO^S, rho_)
    !    elsewhere
    !       v(ixO^S) = zero
    !    endwhere
    ! end if
  end subroutine hd_get_v

  ! Calculate cmax_idim = csound + abs(v_idim) within ixO^L
  subroutine hd_get_cmax(w, x, ixI^L, ixO^L, idim, cmax, cmin)
    use mod_global_parameters

    integer, intent(in)                       :: ixI^L, ixO^L, idim
    double precision, intent(in)              :: w(ixI^S, nw), x(ixI^S, 1:^ND)
    double precision, intent(inout)           :: cmax(ixG^T)
    double precision, intent(inout), optional :: cmin(ixG^T)
    double precision                          :: csound(ixG^T)
    double precision                          :: v(ixG^T)

    call hd_get_v(w, x, ixI^L, ixO^L, idim, v)

    ! if (hd_energy .or. hd_adiab > 0.0d0) then
    call hd_get_pthermal(w, x, ixI^L, ixO^L, csound)
    csound(ixO^S) = sqrt(hd_gamma * csound(ixO^S) * &
         hd_inv_rho(w, ixI^L, ixO^L))

    if (present(cmin)) then
       cmax(ixO^S) = max(v(ixO^S)+csound(ixO^S), zero)
       cmin(ixO^S) = min(v(ixO^S)-csound(ixO^S), zero)
    else
       cmax(ixO^S) = abs(v(ixO^S))+csound(ixO^S)
    end if

    ! else
    ! Jannis: removed this case
    ! case of zero temperature: allow zero density
    ! cmax(ixO^S) = max(v(ixO^S), zero)
    ! cmin(ixO^S) = min(v(ixO^S), zero)
    ! end if

    ! TODO
    ! call dust_get_cmax(...)
  end subroutine hd_get_cmax

  ! Calculate thermal pressure=(gamma-1)*(e-0.5*m**2/rho) within ixO^L
  subroutine hd_get_pthermal(w, x, ixI^L, ixO^L, pth)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L
    double precision             :: w(ixI^S, nw), pth(ixG^T)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    integer, dimension(ixG^T)    :: patchierror
    integer                      :: lowpindex(ndim), idir

    if (fixsmall) call smallvalues(w, x, ixI^L, ixO^L,"hd_get_pthermal")

    if (hd_energy) then
       pth(ixO^S) = (hd_gamma - 1.0d0) * (w(ixO^S, e_) - &
            hd_kin_en(w, ixI^L, ixO^L))
    else
       pth(ixO^S) = hd_adiab*w(ixO^S, rho_)**hd_gamma
    end if

    if (any(pth(ixO^S)<minp)) then
       if (strictsmall) then
          lowpindex = minloc(pth(ixO^S))
          ^D&lowpindex(^D)=lowpindex(^D)+ixOmin^D-1;
          write(*,*) "too small pressure = ", minval(pth(ixO^S))," with limit=", minp,&
               " at x=", x(^D&lowpindex(^D), 1:ndim), lowpindex
          if (hd_energy) then
             write(*,*)" E_total=", w(^D&lowpindex(^D), e_)
          end if
          write(*,*)" w(1:nwflux) =", w(^D&lowpindex(^D), 1:nwflux)," when t=", t," it=", it
          call mpistop("=== strictsmall in hd_get_pthermal ===")
       else if (strictgetaux) then
          if (.not. hd_energy) then
             where(pth(ixO^S)<minp)
                w(ixO^S, rho_) = minrho
             end where

             do idir = 1, ndir
                where(pth(ixO^S)<minp)
                   w(ixO^S, mom(idir)) = zero
                end where
             end do
          end if
          where(pth(ixO^S)<minp)
             pth(ixO^S) = minp
          endwhere
       else
          where(pth(ixO^S)<minp)
             patchierror(ixO^S) = 1
          else where
             patchierror(ixO^S) = 0
          end where

          if (any(patchierror(ixO^S)/= 0)) then
             call correctaux(ixI^L, ixO^L, w, x, patchierror,"hd_get_pthermal")
             if (hd_energy) then
                where(patchierror(ixO^S)/= 0)
                   pth(ixO^S) = (hd_gamma - 1.0d0) * (w(ixO^S, e_)- &
                        hd_kin_en(w, ixI^L, ixO^L))
                end where
             else
                where(patchierror(ixO^S)/= 0)
                   pth(ixO^S) = hd_adiab*w(ixO^S, rho_)**hd_gamma
                end where
             end if

          end if
       end if
    end if

  end subroutine hd_get_pthermal

  ! Calculate non-transport flux f_idim[iw] within ixO^L.
  subroutine hd_get_flux(w, x, ixI^L, ixO^L, iw, idim, f, transport)
    use mod_global_parameters
    ! use mod_dust, only: dust_get_flux

    integer, intent(in)             :: ixI^L, ixO^L, iw, idim
    double precision, intent(in)    :: w(ixI^S, 1:nw), x(ixI^S, 1:^ND)
    double precision, intent(inout) :: f(ixG^T)
    logical, intent(out)            :: transport

    ! TODO: reorganize this, maybe compute all fluxes in a direction at once?

    transport = .true.

    if (iw == mom(idim)) then
       ! f_i[m_i]= v_i*m_i + p
       call hd_get_pthermal(w, x, ixI^L, ixO^L, f)
    else if (iw == e_) then
       ! f_i[e]= v_i*e + m_i/rho*p
       call hd_get_pthermal(w, x, ixI^L, ixO^L, f)
       f(ixO^S) = w(ixO^S, mom(idim))/w(ixO^S, rho_)*f(ixO^S)
       ! else if (iw > hd_nwflux) then
       !    ! A dust flux
       !    ! call dust_get_flux(w, x, ixI^L, ixO^L, iw, idim, f, transport)
    else
       f(ixO^S) = zero
    endif

  end subroutine hd_get_flux

  !> Add geometrical source terms to w
  !>
  !> Notice that the expressions of the geometrical terms depend only on ndir,
  !> not ndim. Eg, they are the same in 2.5D and in 3D, for any geometry.
  !>
  !> Ileyk : to do :
  !>     - give the possibility to set angmomfix=.true.
  !>     - address the source term for the dust
  subroutine hd_add_source_geom(qdt, ixI^L, ixO^L, wCT, w, x) ! - - - - - - - - - - - -
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, x(ixI^S, 1:ndim)
    double precision, intent(inout) :: wCT(ixI^S, 1:nw), w(ixI^S, 1:nw)
    double precision                :: tmp(ixI^S)
    integer                         :: iw,idir,h1x^L,h2x^L
    ! to change and to set as a parameter in the parfile once the possibility to
    ! solve the equations in an angular momentum conserving form has been
    ! implemented (change tvdlf.t eg)
    logical                         :: angmomfix = .false.

    call mpistop("Not implemented yet")

    ! if (ndir==3) then

    !    select case (typeaxial)
    !    case ("slab")
    !       ! No source terms in slab symmetry
    !    case ("cylindrical")
    !       ! s[mr]=(pthermal+mphi**2/rho)/radius
    !       call hd_get_pthermal(wCT,x,ixI^L,ixO^L,tmp)
    !       tmp(ixO^S)=tmp(ixO^S)+wCT(ixO^S,mphi_)**2/wCT(ixO^S,rho_)
    !       w(ixO^S,mr_)=w(ixO^S,mr_)+qdt*tmp(ixO^S)/x(ixO^S,1)
    !       tmp(ixO^S)=zero
    !       ! s[mphi]=(-mphi*mr/rho)/radius
    !       !
    !       ! Ileyk : beware the index permutation : mphi=2 if -phi=2 (2.5D
    !       ! (r,theta) grids) BUT mphi=3 if -phi=3 (for 2.5D (r,z) grids)
    !       if (.not. angmomfix) then
    !          tmp(ixO^S)=-wCT(ixO^S,mphi_)*wCT(ixO^S,mr_)/wCT(ixO^S,rho_)
    !       end if
    !       ! no geometrical source term if angular momentum conserving form of
    !       ! the equations
    !       w(ixO^S,mphi_)=w(ixO^S,mphi_)+qdt*tmp(ixO^S)/x(ixO^S,1)
    !       tmp(ixO^S)=zero
    !    case ("spherical")
    !       h1x^L=ixO^L-kr(1,^D); h2x^L=ixO^L-kr(2,^D)
    !       ! s[mr]=((mtheta**2+mphi**2)/rho+2*p)/r
    !       call hd_get_pthermal(wCT,x,ixI^L,ixO^L,tmp)
    !       tmp(ixO^S)=tmp(ixO^S)*x(ixO^S,1) &
    !            *(mygeo%surfaceC1(ixO^S)-mygeo%surfaceC1(h1x^S)) &
    !            /mygeo%dvolume(ixO^S)
    !       tmp(ixO^S)=tmp(ixO^S)+wCT(ixO^S,mphi_  )**2/wCT(ixO^S,rho_)
    !       tmp(ixO^S)=tmp(ixO^S)+wCT(ixO^S,mtheta_)**2/wCT(ixO^S,rho_)
    !       w(ixO^S,mr_)=w(ixO^S,mr_)+qdt*tmp(ixO^S)/x(ixO^S,1)
    !       tmp(ixO^S)=zero
    !       ! s[mtheta]=-(mr*mtheta/rho)/r+cot(theta)*(mphi**2/rho+p)/r
    !       call hd_get_pthermal(wCT,x,ixI^L,ixO^L,tmp)
    !       tmp(ixO^S)=tmp(ixO^S)*x(ixO^S,1) &
    !            *(mygeo%surfaceC2(ixO^S)-mygeo%surfaceC2(h2x^S)) &
    !            /mygeo%dvolume(ixO^S)
    !       tmp(ixO^S)=tmp(ixO^S)+(wCT(ixO^S,mphi_)**2/wCT(ixO^S,rho_))/dtan(x(ixO^S,2))
    !       if (.not. angmomfix) tmp(ixO^S)=tmp(ixO^S)-(wCT(ixO^S,mtheta_)*wCT(ixO^S,mr_))/wCT(ixO^S,rho_)
    !       w(ixO^S,mtheta_)=w(ixO^S,mtheta_)+qdt*tmp(ixO^S)/x(ixO^S,1)
    !       tmp(ixO^S)=zero
    !       ! s[mphi]=-(mphi*mr/rho)/r-cot(theta)*(mtheta*mphi/rho)/r
    !       if (.not. angmomfix) tmp(ixO^S)=          -(wCT(ixO^S,mtheta_)*wCT(ixO^S,mphi_))/wCT(ixO^S,rho_)
    !       tmp(ixO^S)=tmp(ixO^S)/dtan(x(ixO^S,2))
    !       if (.not. angmomfix) tmp(ixO^S)=tmp(ixO^S)-(wCT(ixO^S,mphi_  )*wCT(ixO^S,mr_  ))/wCT(ixO^S,rho_)
    !       w(ixO^S,mphi_)=w(ixO^S,mphi_)+qdt*tmp(ixO^S)/x(ixO^S,1)
    !       tmp(ixO^S)=zero
    !    case default
    !       call mpistop("typeaxial is slab, cylindrical or spherical")
    !    end select

    ! elseif (ndir==2) then

    !    select case (typeaxial)
    !    case ("slab")
    !       ! No source terms in slab symmetry
    !    case ("cylindrical")
    !       ! (r,phi) : same as ndir==3
    !       ! phi true if and only if -d=22 -phi=2 (and typeaxial==cyl)
    !       if (phi) then ! Ileyk : new argument "phi" here for the parfile. Make sense just if typeaxial==cyl.
    !          ! s[mr]=(pthermal+mphi**2/rho)/radius
    !          call hd_get_pthermal(wCT,x,ixI^L,ixO^L,tmp)
    !          tmp(ixO^S)=tmp(ixO^S)+wCT(ixO^S,mphi_)**2/wCT(ixO^S,rho_)
    !          w(ixO^S,mr_)=w(ixO^S,mr_)+qdt*tmp(ixO^S)/x(ixO^S,1)
    !          tmp(ixO^S)=zero
    !          ! s[mphi]=(-mphi*mr/rho)/radius
    !          if (.not. angmomfix) tmp(ixO^S)=-wCT(ixO^S,mphi_)*wCT(ixO^S,mr_)/wCT(ixO^S,rho_)
    !          w(ixO^S,mphi_)=w(ixO^S,mphi_)+qdt*tmp(ixO^S)/x(ixO^S,1)
    !          tmp(ixO^S)=zero
    !          ! (r,z) : no mphi, just the pressure in the geom. source term
    !       else
    !          ! s[mr]=pthermal/radius
    !          call hd_get_pthermal(wCT,x,ixI^L,ixO^L,tmp)
    !          w(ixO^S,mr_)=w(ixO^S,mr_)+qdt*tmp(ixO^S)/x(ixO^S,1)
    !          tmp(ixO^S)=zero
    !       endif
    !    case ("spherical") ! (r,theta), w/ theta the colatitude. No mphi
    !       ! s[mr]=((mtheta**2)/rho+2*p)/r
    !       call hd_get_pthermal(wCT,x,ixI^L,ixO^L,tmp)
    !       tmp(ixO^S)=tmp(ixO^S)*x(ixO^S,1) &
    !            *(mygeo%surfaceC1(ixO^S)-mygeo%surfaceC1(h1x^S)) &
    !            /mygeo%dvolume(ixO^S)
    !       tmp(ixO^S)=tmp(ixO^S)+wCT(ixO^S,mtheta_)**2/wCT(ixO^S,rho_)
    !       w(ixO^S,mr_)=w(ixO^S,mr_)+qdt*tmp(ixO^S)/x(ixO^S,1)
    !       tmp(ixO^S)=zero
    !       ! s[mtheta]=-(mr*mtheta/rho)/r+cot(theta)*p/r
    !       call hd_get_pthermal(wCT,x,ixI^L,ixO^L,tmp)
    !       tmp(ixO^S)=tmp(ixO^S)*x(ixO^S,1) &
    !            *(mygeo%surfaceC2(ixO^S)-mygeo%surfaceC2(h2x^S)) &
    !            /mygeo%dvolume(ixO^S)
    !       if (.not. angmomfix) tmp(ixO^S)=tmp(ixO^S)-(wCT(ixO^S,mtheta_)*wCT(ixO^S,mr_))/wCT(ixO^S,rho_)
    !       w(ixO^S,mtheta_)=w(ixO^S,mtheta_)+qdt*tmp(ixO^S)/x(ixO^S,1)
    !       tmp(ixO^S)=zero
    !    case default
    !       call mpistop("typeaxial is slab, cylindrical or spherical")
    !    end select

    ! elseif (ndir==1) then

    !    select case (typeaxial)
    !    case ("slab")
    !       ! No source terms in slab symmetry
    !    case ("cylindrical")
    !       ! s[mr]=pthermal/radius
    !       call hd_get_pthermal(wCT,x,ixI^L,ixO^L,tmp)
    !       w(ixO^S,mr_)=w(ixO^S,mr_)+qdt*tmp(ixO^S)/x(ixO^S,1)
    !       tmp(ixO^S)=zero
    !    case ("spherical")
    !       ! s[mr]=2pthermal/radius
    !       call hd_get_pthermal(wCT,x,ixI^L,ixO^L,tmp)
    !       tmp(ixO^S)=tmp(ixO^S)*x(ixO^S,1) &
    !            *(mygeo%surfaceC1(ixO^S)-mygeo%surfaceC1(h1x^S)) &
    !            /mygeo%dvolume(ixO^S)
    !       w(ixO^S,mr_)=w(ixO^S,mr_)+qdt*tmp(ixO^S)/x(ixO^S,1)
    !       tmp(ixO^S)=zero
    !    case default
    !       call mpistop("typeaxial is slab, cylindrical or spherical")
    !    end select

    ! endif

    ! if (fixsmall) call smallvalues(w, x, ixI^L, ixO^L,"addgeometry")

  end subroutine hd_add_source_geom

  ! w[iw]= w[iw]+qdt*S[wCT, qtC, x] where S is the source based on wCT within ixO
  subroutine hd_add_source(qdt, ixI^L, ixO^L, iw^LIM, qtC, wCT, qt, w, x, qsourcesplit)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L, iw^LIM
    double precision, intent(in)    :: qdt, qtC, qt, x(ixI^S, 1:ndim)
    double precision, intent(inout) :: wCT(ixI^S, 1:nw), w(ixI^S, 1:nw)
    logical, intent(in)             :: qsourcesplit

    double precision, dimension(ixG^T, 1:ndir, 1:^NDS) :: fdrag
    integer                                           :: idir

    ! TODO
    ! call dust_add_source(qdt, ixI^L, ixO^L, iw^LIM, &
    !      qtC, wCT, qt, w, x, qsourcesplit)
  end subroutine hd_add_source

  subroutine hd_get_dt(w, ixI^L, ixO^L, dtnew, dx^D, x)
    use mod_global_parameters

    integer, intent(in)                        :: ixI^L, ixO^L
    double precision, intent(in)               :: dx^D, x(ixI^S, 1:ndim)
    double precision, intent(inout)            :: w(ixI^S, 1:nw), dtnew
    double precision :: dtdust

    dtnew = bigdouble

    ! TODO
    ! call dust_get_dt(w, ixI^L, ixO^L, dtnew, dx^D, x)
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

  subroutine smallvalues(w,x,ixI^L,ixO^L,subname)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    character(len=*), intent(in)    :: subname
    integer                         :: posvec(ndim)
    integer, dimension(ixG^T)       :: patchierror
    double precision                :: pth(ixG^T)

    patchierror(ixO^S) = 0

    if (hd_energy) then
       pth(ixO^S) = (hd_gamma - 1.0d0) * (w(ixO^S, e_) - &
            hd_kin_en(w, ixI^L, ixO^L))

       where(w(ixO^S,rho_) < minrho .or. w(ixO^S,e_) < smalle &
            .or. pth(ixO^S)<minp)
          patchierror(ixO^S)=-1
       end where
    else
       where(w(ixO^S,rho_) < minrho)
          patchierror(ixO^S)=-1
       end where
    end if

    if (any(patchierror(ixO^S) /= 0)) then
       if (strictsmall) then
          write(*,*)"SMALLVALUES under strictsmall problem From:  ", &
               subname," iteration ", it," time ",t
          call mpistop("Smallvalues with strictsmall=T failed")
       else if(strictgetaux)then
          call mpistop("Smallvalues with strictgetaux=T not implemented")
          ! where(w(ixO^S,rho_) < minrho .or. w(ixO^S,e_) < smalle&
          !      .or. pth(ixO^S)<minp)
          !    w(ixO^S,rho_)  = 2.0*(1.0d0 + 10.0d0 * minrho)*minrho
          !    w(ixO^S,e_)    = 2.0*(1.0d0 + 10.0d0 * minp)*smalle
          !    {^C&w(ixO^S,m^C_) =zero;}
          ! end where
       else
          call correctaux(ixI^L,ixO^L,w,x,patchierror,subname)
       end if
    end if
  end subroutine smallvalues

  subroutine correctaux(ixI^L,ixO^L,w,x,patchierror,subname)
    use mod_global_parameters

    integer, intent(in)         :: ixI^L, ixO^L
    integer, intent(in)         :: patchierror(ixG^T)
    character(len=*), intent(in)   :: subname
    double precision, intent(inout):: w(ixI^S,1:nw)
    double precision, intent(in)      :: x(ixI^S,1:ndim)
    integer        :: iw, kxO^L, ix^D, i

    {do ix^DB= ixO^LIM^DB\}
    ! point with local failure identified by patchierror
    ! patchierror=-1 then from smallvalues call
    ! patchierror=1  then from getpthermal or primitive call
    if (patchierror(ix^D)/=0) then
       ! verify in cube with border width nflatgetaux the presence
       ! of cells where all went ok
       do i=1,nflatgetaux
          {kxOmin^D= max(ix^D-i,ixOmin^D);
          kxOmax^D= min(ix^D+i,ixOmax^D);\}
          ! in case cells are fine within smaller cube than 
          ! the userset nflatgetaux: use that smaller cube
          if (any(patchierror(kxO^S)==0)) exit
       end do
       if (any(patchierror(kxO^S)==0))then
          ! within surrounding cube, cells without problem were found
          if (useprimitive .and. subname/="hd_to_primitive") then
             call hd_to_primitive(ixI^L,kxO^L,w,x, .false.)
          end if
          ! faulty cells are corrected by averaging here
          ! only average those which were ok and replace faulty cells
          do iw = 1,nw
             w(ix^D,iw)=sum(w(kxO^S,iw),patchierror(kxO^S)==0)&
                  /count(patchierror(kxO^S)==0)
          end do
          if (useprimitive .and. subname/="hd_to_primitive") then
             ! in addition to those switched to primitive variables
             ! above, also switch the corrected variables
             call hd_to_conserved(ixI^L,kxO^L,w,x, .false.)
          end if
       else
          ! no cells without error were found in cube of size nflatgetaux
          ! --> point of no recovery
          print*,"Getaux error:",patchierror(ix^D),"ix^D=",ix^D
          !print*,"New ","rho=",w(ix^D,rho_),"m=", &
          !        {^C&w(ix^D,m^C_)},"e=",w(ix^D,e_)
          !print*,"position  ", px(saveigrid)%x(ix^D, 1:ndim)
          print*,"Called from: ",subname
          if (patchierror(ix^D)<0) then
             call mpistop("-correctaux from smallvalues-----")
          else
             call mpistop("-correctaux from primitive or getpthermal--")
          end if
       end if
    end if
    {enddo^D&\}

  end subroutine correctaux

end module mod_hd_phys
