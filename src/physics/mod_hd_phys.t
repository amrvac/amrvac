module mod_hd_phys
  use mod_physics

  implicit none
  private

  logical, public, protected              :: hd_energy   = .true.
  integer, public, protected              :: hd_n_tracer = 0
  integer, public, parameter              :: rho_        = 1
  integer, public, parameter              :: m0_         = rho_
  integer, allocatable, public, protected :: mom(:)
  integer, allocatable, public, protected :: tracer(:)
  integer, public, protected              :: e_
  double precision, public, protected     :: hd_gamma    = 1.4d0
  double precision, public, protected     :: hd_adiab    = 1.0d0
  double precision, protected             :: smalle, minrho, minp

  ! Public methods
  public :: hd_phys_init
  public :: hd_kin_en
  public :: hd_get_pthermal
  public :: hd_get_v

contains

  subroutine hd_params_read(par_files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: par_files(:)
    integer                      :: n

    namelist /hd_list/ hd_energy, hd_n_tracer, hd_gamma, hd_adiab

    do n = 1, size(par_files)
       open(unitpar, file=trim(par_files(n)), status='old')
       read(unitpar, hd_list, end=111)
111    close(unitpar)
    end do

  end subroutine hd_params_read

  subroutine hd_phys_init(par_files)
    use mod_global_parameters

    integer :: itr, idir

    character(len=*), intent(in) :: par_files(:)

    call hd_params_read(par_files)

    physics_type = 'hd'

    ! The number of flux variables
    ! TODO: add dust
    if (hd_energy) then
       nwflux = 2 + hd_n_tracer + ndir
    else
       nwflux = 1 + hd_n_tracer + ndir
    end if

    nwaux   = 0
    nwextra = 0
    nw      = nwflux + nwaux + nwextra

    nvector      = 1 ! No. vector vars
    allocate(iw_vector(nvector))
    iw_vector(1) = m0_

    ! TODO: dust related velocity vectors not handled here

    nflag_ = nw + 1

    allocate(mom(ndir))
    do idir = 1, ndir
       mom(idir) = m0_ + idir
    end do

    ! Set index of energy variable
    if (hd_energy) then
       e_ = mom(ndir) + 1
    else
       e_ = -1
    end if

    allocate(tracer(hd_n_tracer))

    ! Set starting index of tracers
    do itr = 1, hd_n_tracer
       tracer(itr) = max(e_, mom(ndir)) + itr
    end do

    ! TODO
    ! call dust_activate()

    phys_get_v           => hd_get_v
    phys_get_cmax        => hd_get_cmax
    phys_get_flux        => hd_get_flux
    phys_add_source_geom => hd_add_source_geom
    phys_to_conserved => hd_to_conserved
    phys_to_primitive => hd_to_primitive
    phys_check_params => hd_check_params

  end subroutine hd_phys_init

  subroutine hd_read_params(file_unit, success)
    integer, intent(in) :: file_unit
    logical, intent(out) :: success

    namelist /hd_list/ hd_energy, hd_n_tracer, hd_gamma, hd_adiab

    ! Success indicates if read succeeds
    success = .false.
    read(file_unit, hd_list, end=101)
    success = .true.

101 return
  end subroutine hd_read_params

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

  ! set default values for entropy fixes for 'yee' type
  subroutine initglobaldata

    if (hd_energy) then
       hd_gamma = 5.d0/3.d0
    else
       hd_gamma = 1.d0
       hd_adiab = 1.d0
    end if

    ! if (dust_num_species > 0) then
    !    eqpar(mu_) = 1.0d0
    !    mhcgspar = 1.6733D-24
    !    kbcgspar = 1.38065D-16
    ! end if

  end subroutine initglobaldata

  ! Calculate auxilary variables ixO^L from non-auxiliary entries in w
  ! clipping can be set to .true. to e.g. correct unphysical pressures,
  ! densities, v > c,  etc.
  subroutine getaux(clipping, w, x, ixI^L, ixO^L, subname)
    use mod_global_parameters

    logical                         :: clipping
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    character(len=*)                :: subname

  end subroutine getaux

  subroutine checkw(checkprimitive, ixI^L, ixO^L, w, flag)
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

  end subroutine checkw

  ! Transform primitive variables into conservative ones
  subroutine hd_to_conserved(ixI^L, ixO^L, w, x)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    double precision                :: invgam
    integer                         :: idir, itr

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

    if (fixsmall) call smallvalues(w, x, ixI^L, ixO^L, "hd_to_conserved")

  end subroutine hd_to_conserved

  ! Transform conservative variables into primitive ones
  subroutine hd_to_primitive(ixI^L, ixO^L, w, x)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    integer, dimension(ixG^T)       :: patchierror
    integer, dimension(ndim)        :: lowpindex
    integer                         :: itr, idir

    if (fixsmall) call smallvalues(w, x, ixI^L, ixO^L, "hd_to_primitive")

    if (hd_energy) then
       ! Compute pressure
       w(ixO^S, e_) = (hd_gamma - 1.0d0) * (w(ixO^S, e_) - &
            hd_kin_en(w, ixI^L, ixO^L))
    end if

    ! Convert momentum to velocity
    do idir = 1, ndir
       w(ixO^S, mom(idir)) = w(ixO^S, mom(idir))/w(ixO^S, rho_);
    end do

    ! We got rho, Dtr, now we can get the tracers
    do itr = 1, hd_n_tracer
       w(ixO^S, tracer(itr)) = w(ixO^S, tracer(itr)) / w(ixO^S, rho_)
    end do

    ! if (strictsmall) then
    !    if (any(w(ixO^S, p_)<minp)) then
    !       lowpindex = minloc(w(ixO^S, p_))
    !       ^D&lowpindex(^D) = lowpindex(^D)+ixOmin^D-1;
    !       write(*,*)'too small pressure = ', minval(w(ixO^S, p_)),' with limit=', minp,&
    !            ' at x=', x(^D&lowpindex(^D), 1:ndim), lowpindex,' where E_k=',&
    !            half*(^C&w(^D&lowpindex(^D), mom(idir))**2+)/w(^D&lowpindex(^D), rho_),&
    !            ' E_total=', w(^D&lowpindex(^D), e_),&
    !            ' w(1:nwflux) =', w(^D&lowpindex(^D), 1:nwflux),' when t=', t,' it=', it
    !       call mpistop("=== primitive pressure problem===")
    !    end if
    ! else if (strictgetaux) then
    !    ! TODO: check
    !    where(w(ixO^S, p_)<minp)
    !       w(ixO^S, p_) = minp
    !    endwhere
    ! else
    !    where(w(ixO^S, p_)<minp)
    !       patchierror(ixO^S) = 1
    !    else where
    !       patchierror(ixO^S) = 0
    !    end where

    !    if (any(patchierror(ixO^S)/= 0)) &
    !         call correctaux(ixI^L, ixO^L, w, x, patchierror,'primitive')
    ! end if


    ! Jannis: Removed (hd_adiab == 0.0)
    ! case of zero temperature: allow zero density
    ! where(w(ixO^S, rho_)/= zero)
    !    ^C&w(ixO^S, v^C_) = w(ixO^S, mom(idir))/w(ixO^S, rho_);
    ! elsewhere
    !    ^C&w(ixO^S, v^C_) = zero;
    ! endwhere

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

  ! Calculate v_i = m_i/rho within ixO^L
  subroutine hd_get_v(w, x, ixI^L, ixO^L, idim, v)
    use mod_global_parameters

    integer, intent(in)           :: ixI^L, ixO^L, idim
    double precision, intent(in)  :: w(ixI^S, nw), x(ixI^S, 1:^ND)
    double precision, intent(out) :: v(ixG^T)

    ! if (hd_energy .or. hd_adiab > 0.0d0) then
    v(ixO^S) = w(ixO^S, mom(idim)) / w(ixO^S, rho_)
    ! Jannis: Removed this case
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
    csound(ixO^S) = sqrt(hd_gamma*csound(ixO^S)/w(ixO^S, rho_))

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
  subroutine hd_get_pthermal(w, x, ixI^L, ixO^L, p)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L
    double precision             :: w(ixI^S, nw), p(ixG^T)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    integer, dimension(ixG^T)    :: patchierror
    integer                      :: lowpindex(ndim), idir

    if (fixsmall) call smallvalues(w, x, ixI^L, ixO^L,"hd_get_pthermal")

    if (hd_energy) then
       p(ixO^S) = (hd_gamma - 1.0d0) * (w(ixO^S, e_) - &
            hd_kin_en(w, ixI^L, ixO^L))
    else
       p(ixO^S) = hd_adiab*w(ixO^S, rho_)**hd_gamma
    end if

    if (any(p(ixO^S)<minp)) then
       if (strictsmall) then
          lowpindex = minloc(p(ixO^S))
          ^D&lowpindex(^D)=lowpindex(^D)+ixOmin^D-1;
          write(*,*)'too small pressure = ', minval(p(ixO^S)),' with limit=', minp,&
               ' at x=', x(^D&lowpindex(^D), 1:ndim), lowpindex
          if (hd_energy) then
             write(*,*)' E_total=', w(^D&lowpindex(^D), e_)
          end if
          write(*,*)' w(1:nwflux) =', w(^D&lowpindex(^D), 1:nwflux),' when t=', t,' it=', it
          call mpistop("=== strictsmall in hd_get_pthermal ===")
       else
          if (strictgetaux) then
             where(p(ixO^S)<minp)
                p(ixO^S) = minp
             endwhere

             if (.not. hd_energy) then
                where(p(ixO^S)<minp)
                   w(ixO^S, rho_) = minrho
                end where

                do idir = 1, ndir
                   where(p(ixO^S)<minp)
                      w(ixO^S, mom(idir)) = zero
                   end where
                end do
             end if
          else

             where(p(ixO^S)<minp)
                patchierror(ixO^S) = 1
             else where
                patchierror(ixO^S) = 0
             end where

             if (any(patchierror(ixO^S)/= 0)) then
                call correctaux(ixI^L, ixO^L, w, x, patchierror,'hd_get_pthermal')
                if (hd_energy) then
                   where(patchierror(ixO^S)/= 0)
                      p(ixO^S) =(hd_gamma - 1.0d0)*(w(ixO^S, e_)- &
                           hd_kin_en(w, ixI^L, ixO^L))
                   end where
                else
                   where(patchierror(ixO^S)/= 0)
                      p(ixO^S) = hd_adiab*w(ixO^S, rho_)**hd_gamma
                   end where
                end if

             end if
          end if
       end if
    end if

  end subroutine hd_get_pthermal

  ! Calculate non-transport flux f_idim[iw] within ixO^L.
  subroutine hd_get_flux(w, x, ixI^L, ixO^L, iw, idim, f, transport)
    use mod_global_parameters

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
    else
       f(ixO^S) = zero
    endif

    ! TODO: dust
  end subroutine hd_get_flux

  ! Add geometrical source terms to w
  ! Notice that the expressions of the geometrical terms depend only on ndir,
  ! not ndim. Eg, they are the same in 2.5D and in 3D, for any geometry.
  ! Ileyk : to do :
  !     - give the possibility to set angmomfix=.true.
  !     - address the source term for the dust
  subroutine hd_add_source_geom(qdt, ixI^L, ixO^L, wCT, w, x) ! - - - - - - - - - - - -
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, x(ixI^S, 1:ndim)
    double precision, intent(inout) :: wCT(ixI^S, 1:nw), w(ixI^S, 1:nw)
    double precision :: tmp(ixI^S)
    integer :: iw,idir,h1x^L,h2x^L
    ! to change and to set as a parameter in the parfile once the possibility to
    ! solve the equations in an angular momentum conserving form has been
    ! implemented (change tvdlf.t eg)
    logical :: angmomfix=.false.

    call mpistop("Not implemented yet")

    ! if (ndir==3) then

    !    select case (typeaxial)
    !    case ('slab')
    !       ! No source terms in slab symmetry
    !    case ('cylindrical')
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
    !    case ('spherical')
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
    !    case ('slab')
    !       ! No source terms in slab symmetry
    !    case ('cylindrical')
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
    !    case ('spherical') ! (r,theta), w/ theta the colatitude. No mphi
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
    !    case ('slab')
    !       ! No source terms in slab symmetry
    !    case ('cylindrical')
    !       ! s[mr]=pthermal/radius
    !       call hd_get_pthermal(wCT,x,ixI^L,ixO^L,tmp)
    !       w(ixO^S,mr_)=w(ixO^S,mr_)+qdt*tmp(ixO^S)/x(ixO^S,1)
    !       tmp(ixO^S)=zero
    !    case ('spherical')
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

  subroutine getdt(w, ixI^L, ixO^L, dtnew, dx^D, x)
    use mod_global_parameters

    integer, intent(in)                        :: ixI^L, ixO^L
    double precision, intent(in)               :: dx^D, x(ixI^S, 1:ndim)
    double precision, intent(inout)            :: w(ixI^S, 1:nw), dtnew
    double precision :: dtdust

    ! TODO
    ! call dust_get_dt(w, ixI^L, ixO^L, dtnew, dx^D, x)
    dtnew = bigdouble
  end subroutine getdt

  function hd_kin_en(w, ixI^L, ixO^L) result(ke)
    use mod_global_parameters, only: nw, ndim
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S, nw)
    double precision              :: ke(ixO^S)

    ke = 0.5d0 * sum(w(ixO^S, mom(:))**2, dim=ndim+1) / &
         w(ixO^S, rho_)
  end function hd_kin_en

  subroutine smallvalues(w,x,ixI^L,ixO^L,subname)

    use mod_global_parameters

    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    character(len=*), intent(in)    :: subname
    !.. local ..
    integer                         :: posvec(ndim)
    integer, dimension(ixG^T)       :: patchierror
    double precision                :: pth(ixG^T)
    !-----------------------------------------------------------------------------

    call mpistop("smallvalues not implemented yet")
    ! TODO
    ! if (hd_energy) then
    !    pth(ixO^S)=(eqpar(gamma_)-one)*(w(ixO^S,e_)- &
    !         half*(({^C&w(ixO^S,m^C_)**2+})/w(ixO^S,rho_)))
    !    if(any(w(ixO^S,rho_) < minrho).or.any(w(ixO^S,e_) < smalle)&
    !         .or.any(pth(ixO^S) < minp))then
    !       if(strictsmall)then
    !          write(*,*)'SMALLVALUES under strictsmall problem From:  ', &
    !               subname,' iteration ', it,' time ',t
    !          if(any(w(ixO^S,rho_) < minrho)) then
    !             posvec(1:ndim)=minloc(w(ixO^S,rho_))
    !             ^D&posvec(^D)=posvec(^D)+ixOmin^D-1;
    !             write(*,*)'minimum rho =', minval(w(ixO^S,rho_)),' with limit=',minrho,&
    !                  ' at x=', x(^D&posvec(^D),1:ndim),' array index=',posvec,' where E_k=',&
    !                  half*(^C&w(^D&posvec(^D),m^C_)**2+)/w(^D&posvec(^D),rho_),&
    !                  ' E_total=',w(^D&posvec(^D),e_),&
    !                  ' w(1:nwflux)=',w(^D&posvec(^D),1:nwflux)
    !          endif
    !          if(any(w(ixO^S,e_) < smalle)) then
    !             posvec(1:ndim)=minloc(w(ixO^S,e_))
    !             ^D&posvec(^D)=posvec(^D)+ixOmin^D-1;
    !             write(*,*)'minimum e =', minval(w(ixO^S,e_)),' with limit=',smalle,&
    !                  ' at x=', x(^D&posvec(^D),1:ndim),' array index=',posvec,' where E_k=',&
    !                  half*(^C&w(^D&posvec(^D),m^C_)**2+)/w(^D&posvec(^D),rho_),&
    !                  ' E_total=',w(^D&posvec(^D),e_),&
    !                  ' w(1:nwflux)=',w(^D&posvec(^D),1:nwflux)
    !          endif
    !          if(any(pth(ixO^S) < minp)) then
    !             posvec(1:ndim)=minloc(pth(ixO^S))
    !             ^D&posvec(^D)=posvec(^D)+ixOmin^D-1;
    !             write(*,*)'minimum pressure = ', minval(pth(ixO^S)),' with limit=',minp,&
    !                  ' at x=',x(^D&posvec(^D),1:ndim),' array index=',posvec,' where rho=',&
    !                  w({^D&posvec(^D)},rho_),', velocity v=',&
    !                  ^C&w({^D&posvec(^D)},m^C_)/w({^D&posvec(^D)},rho_),&
    !                  ' w(1:nwflux)=',w(^D&posvec(^D),1:nwflux)
    !          endif
    !          call mpistop("Smallvalues with strictsmall=T failed")
    !       else
    !          if(strictgetaux)then
    !             where(w(ixO^S,rho_) < minrho .or. w(ixO^S,e_) < smalle&
    !                  .or. pth(ixO^S)<minp)
    !                w(ixO^S,rho_)  = 2.0*(1.0d0 + 10.0d0 * minrho)*minrho
    !                w(ixO^S,e_)    = 2.0*(1.0d0 + 10.0d0 * minp)*smalle
    !                {^C&w(ixO^S,m^C_) =zero;}
    !             end where
    !          else
    !             where(w(ixO^S,rho_) < minrho .or. w(ixO^S,e_) < smalle&
    !                  .or. pth(ixO^S)<minp)
    !                patchierror(ixO^S)=-1
    !             elsewhere
    !                patchierror(ixO^S)=0
    !             end where
    !             call correctaux(ixI^L,ixO^L,w,x,patchierror,subname)
    !          end if
    !       end if ! strict
    !    end if
    ! else                        ! No energy
    !    if(any(w(ixO^S,rho_) < minrho))then
    !       if(strictsmall)then
    !          write(*,*)'SMALLVALUES under strictsmall problem From:  ', &
    !               subname,' iteration ', it,' time ',t
    !          posvec(1:ndim)=minloc(w(ixO^S,rho_))
    !          ^D&posvec(^D)=posvec(^D)+ixOmin^D-1;
    !          write(*,*)'minimum rho =', minval(w(ixO^S,rho_)),' with limit=',minrho,&
    !               ' at x=', x(^D&posvec(^D),1:ndim),' array index=',posvec,' where E_k=',&
    !               half*(^C&w(^D&posvec(^D),m^C_)**2+)/w(^D&posvec(^D),rho_),&
    !               ' w(1:nwflux)=',w(^D&posvec(^D),1:nwflux)
    !          call mpistop("Smallvalues with strictsmall=T failed")
    !       else
    !          if(strictgetaux)then
    !             where(w(ixO^S,rho_) < minrho )
    !                w(ixO^S,rho_)  = 2.0*(1.0d0 + 10.0d0 * minrho)*minrho
    !                {^C&w(ixO^S,m^C_) =zero;}
    !             end where
    !          else
    !             where(w(ixO^S,rho_) < minrho)
    !                patchierror(ixO^S)=-1
    !             elsewhere
    !                patchierror(ixO^S)=0
    !             end where
    !             call correctaux(ixI^L,ixO^L,w,x,patchierror,subname)
    !          end if
    !       end if ! strict
    !    end if
    ! end if
  end subroutine smallvalues

  subroutine correctaux(ixI^L,ixO^L,w,x,patchierror,subname)

    use mod_global_parameters

    integer, intent(in)         :: ixI^L, ixO^L
    integer, intent(in)         :: patchierror(ixG^T)
    character(len=*), intent(in)   :: subname
    double precision, intent(inout):: w(ixI^S,1:nw)
    double precision, intent(in)      :: x(ixI^S,1:ndim)
    integer        :: iw, kxO^L, ix^D, i
    logical        :: patchw(ixG^T)
    !-----------------------------------------------------------------------------

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
          call mpistop("correctaux not implemented yet")
          ! ! within surrounding cube, cells without problem were found
          ! if (useprimitive.and.subname/='primitive') then
          !    patchw(kxO^S)=(patchierror(kxO^S)/=0)
          !    call hd_to_primitive(ixI^L,kxO^L,w,patchw)
          ! end if
          ! ! faulty cells are corrected by averaging here
          ! ! only average those which were ok and replace faulty cells
          ! do iw = 1,nw
          !    w(ix^D,iw)=sum(w(kxO^S,iw),patchierror(kxO^S)==0)&
          !         /count(patchierror(kxO^S)==0)
          ! end do
          ! if (useprimitive.and.subname/='primitive') then
          !    ! in addition to those switched to primitive variables
          !    ! above, also switch the corrected variables
          !    patchw(ix^D)=.false.
          !    call conserven(ixI^L,kxO^L,w,patchw)
          ! end if
       else
          ! no cells without error were found in cube of size nflatgetaux
          ! --> point of no recovery
          print*,'Getaux error:',patchierror(ix^D),'ix^D=',ix^D
          !print*,'New ','rho=',w(ix^D,rho_),'m=', &
          !        {^C&w(ix^D,m^C_)},'e=',w(ix^D,e_)
          !print*,'position  ', px(saveigrid)%x(ix^D, 1:ndim)
          print*,'Called from: ',subname
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
