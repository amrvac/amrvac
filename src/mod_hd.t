module mod_hd
  use mod_physics

  implicit none
  private

  logical :: hd_energy
  integer :: hd_n_tracer

  integer, parameter   :: rho_ = 1
  integer, parameter   :: m0_  = rho_
  integer, allocatable :: mom(:)
  integer              :: e_

  integer, allocatable :: tracer(:)
  integer, allocatable :: rho_dust(:)
  integer, allocatable :: m_dust(:, :)

  double precision              :: gamma
  double precision              :: hd_adiab
  double precision              :: smalle, minrho, minp
  double precision              :: minrhod
  double precision, allocatable :: sdust(:), dsdust(:), rhodust(:), mhcgspar, kbcgspar
  double precision              :: Lstar, Tdust

  ! Public methods
  public :: hd_activate

contains

  subroutine hd_activate
    use mod_global_parameters
    integer :: ix

    physics_type = 'hd'

    ! The number of flux variables (dust module can later add to this)
    if (hd_energy) then
       nwflux = 2 + hd_n_tracer + ndir
    else
       nwflux = 1 + hd_n_tracer + ndir
    end if

    nwaux   = 0
    nwextra = 0
    nw      = nwflux + nwaux + nwextra

    ! Note: dust related velocity vectors not handled here
    nvector      = 1 ! No. vector vars
    allocate(iw_vector(nvector))
    iw_vector(1) = m0_

    nflag_ = nw + 1

    allocate(mom(ndir))
    mom(:) = [{m0_+^C}]

    ix = m^NC_

    ! Set index of energy variable
    if (hd_energy) then
       e_ = ix + 1
       ix = ix + 1
    else
       e_ = -1
    end if

    ! Set starting index of tracers
    if (hd_n_tracer > 0) then
       tracer_ = ix + 1
       ix      = ix + hd_n_tracer
    else
       tracer_ = -1
    end if

    if (hd_energy) then
       ! Characteristic waves
       soundRW_ = 1
       soundLW_ = 2
       entropW_ = 3
       shearW0_ = 3
       nworkroe = 3
    else
       ! Characteristic waves
       soundRW_ = 1
       soundLW_ = 2
       shearW0_ = 2
       nworkroe = 1
    end if

    TODO
    activate dust module

  end subroutine hd_activate

  subroutine checkglobaldata
    use mod_global_parameters

    if (.not. hd_energy) then
       if (hd_gamma <= 0.0d0) call mpistop ("gamma negative not ok")
       if (hd_adiab < 0.0d0) call mpistop ("adiab strict negative not ok")
       minrho = max(0.0d0, smallrho)
       minp   = hd_adiab*minrho**hd_gamma
    else
       if (hd_gamma <= 0.0d0 .or. hd_gamma == 1.0d0) &
            call mpistop ("gamma negative or 1 not ok")
       minp   = max(0.0d0, smallp)
       minrho = max(0.0d0, smallrho)
       smalle = minp/(hd_gamma - 1.0d0)
    end if

    if (dust_num_species > 0) then
       if (eqpar(mu_)<= 0.0d0) call mpistop ("mu (molecular weight) negative not ok")
       minrhod = max(0.0d0, smallrhod)
    end if

  end subroutine checkglobaldata

  ! set default values for entropy fixes for 'yee' type
  subroutine initglobaldata
    use mod_global_parameters

    integer :: il

    if (hd_energy) then
       hd_gamma = 5.d0/3.d0
    else
       hd_gamma = 1.d0
       hd_adiab = 1.d0
    end if

    if (dust_num_species > 0) then
       eqpar(mu_) = 1.0d0
       mhcgspar = 1.6733D-24
       kbcgspar = 1.38065D-16
    end if

    do il = 1, nw
       select case(il)
       case(soundRW_, soundLW_)
          entropycoef(il) = 0.2d0
       case default
          entropycoef(il) = -1.0d0
       end select
    end do

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
          flag(ixO^S) =(w(ixO^S, p_)>= minp .and. w(ixO^S, rho_)>= minrho)
       else
          tmp(ixO^S) =(hd_gamma - 1.0d0)*(w(ixO^S, e_)- &
               half*( ^C&w(ixO^S, mom(i_dir))**2+ )/w(ixO^S, rho_))
          flag(ixO^S) =(tmp(ixO^S)>= minp .and. w(ixO^S, rho_)>= minrho)
       endif
    else
       if (hd_adiab > 0.0d0) then
          flag(ixO^S) = (w(ixO^S, rho_) >= 0.0d0)
       endif
    end if

  end subroutine checkw

  ! Transform primitive variables into conservative ones
  subroutine conserve(ixI^L, ixO^L, w, x, patchw)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L
    double precision             :: w(ixI^S, nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    logical                      :: patchw(ixG^T)
    double precision             :: invgam

    invgam = 1.d0/(hd_gamma - 1.0d0)

    if (hd_energy) then
       ! Calculate total energy from pressure and kinetic energy
       w(ixO^S, e_) = w(ixO^S, p_)*invgam+ &
            half*w(ixO^S, rho_)*(^C&w(ixO^S, v^C_)**2+)
    end if

    ! Convert velocity to momentum
    TODO: loop
    ^C&w(ixO^S, mom(i_dir)) = w(ixO^S, rho_)*w(ixO^S, v^C_);

    if (hd_n_tracer > 0) then
       TODO: loop
       {^FL&w(ixO^S, tracer(i_tr)) = w(ixO^S, rho_)*w(ixO^S, tracer(i_tr))\}
    end if

    if (dust_num_species > 0) then
       call dust_conserve(...)
       ! Convert dust velocity to dust momentum
       !{^DS&{^C&w(ixO^S, m_dust(i_dim, i_dust)) = w(ixO^S, rho_dust(i_dust))*w(ixO^S, v^Cd^DS_);}\}
    end if

    if (fixsmall) call smallvalues(w, x, ixI^L, ixO^L,"conserve")

  end subroutine conserve

  ! Transform conservative variables into primitive ones
  subroutine primitive(ixI^L, ixO^L, w, x)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L
    double precision             :: w(ixI^S, nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    integer, dimension(ixG^T)    :: patchierror
    integer, dimension(ndim)     :: lowpindex
    double precision             :: kin_en(ix0^S)

    if (fixsmall) call smallvalues(w, x, ixI^L, ixO^L,"primitive")

    if (hd_energy) then
       ! compute pressure
       kin_en = 0.5d0 * sum(w(ixO^S, mom(:))**2, dim=^ND+1) / w(ixO^S, rho_)
       w(ixO^S, p_) = (hd_gamma - 1.0d0) * (w(ixO^S, e_) - kin_en

       ! Convert momentum to velocity
       do i_dir = 1, ndir
          w(ixO^S, mom(idir)) = w(ixO^S, mom(i_dir))/w(ixO^S, rho_);
       end do

       if (strictsmall) then
          if (any(w(ixO^S, p_)<minp)) then
             lowpindex = minloc(w(ixO^S, p_))
             ^D&lowpindex(^D) = lowpindex(^D)+ixOmin^D-1;
             write(*,*)'too small pressure = ', minval(w(ixO^S, p_)),' with limit=', minp,&
                  ' at x=', x(^D&lowpindex(^D), 1:ndim), lowpindex,' where E_k=',&
                  half*(^C&w(^D&lowpindex(^D), mom(i_dir))**2+)/w(^D&lowpindex(^D), rho_),&
                  ' E_total=', w(^D&lowpindex(^D), e_),&
                  ' w(1:nwflux) =', w(^D&lowpindex(^D), 1:nwflux),' when t=', t,' it=', it
             call mpistop("=== primitive pressure problem===")
          end if
       else if (strictgetaux) then
          ! TODO: check
          where(w(ixO^S, p_)<minp)
             w(ixO^S, p_) = minp
          endwhere
       else
          where(w(ixO^S, p_)<minp)
             patchierror(ixO^S) = 1
          else where
             patchierror(ixO^S) = 0
          end where

          if (any(patchierror(ixO^S)/= 0)) &
               call correctaux(ixI^L, ixO^L, w, x, patchierror,'primitive')
       end if
    end if

    if (.not. hd_energy) then

       if (hd_adiab > 0.0d0) then
          ! Convert momentum to velocity
          ^C&w(ixO^S, v^C_) = w(ixO^S, mom(i_dir))/w(ixO^S, rho_);
       else
          ! case of zero temperature: allow zero density
          where(w(ixO^S, rho_)/= zero)
             ^C&w(ixO^S, v^C_) = w(ixO^S, mom(i_dir))/w(ixO^S, rho_);
          elsewhere
             ^C&w(ixO^S, v^C_) = zero;
          endwhere
       endif
    end if

    if (hd_n_tracer > 0) then
       ! We got rho, Dtr, now we can get the tracers
       do i_tr = 1, hd_n_tracer
          w(ixO^S, tracer(i_tr)) = w(ixO^S, tracer(i_tr)) / w(ixO^S, rho_)\}
       end do
    end if

    ! Convert dust momentum to dust velocity
    call dust_primitive(...)
    ! do i_dust = 1, dust_num_species
    !    where(w(ixO^S, rho_dust(i_dust))>minrhod)
    !       do i_dir = 1, ndir
    !          w(ixO^S, m_dust(i_dir, i_dust)) = w(ixO^S, m_dust(i_dir, i_dust)) /&
    !               w(ixO^S, rho_dust(i_dust))
    !       end do
    !    elsewhere
    !       w(ixO^S, m_dust(:, i_dust)) = 0.0d0
    !    end where
    ! end if

  end subroutine primitive

  subroutine e_to_rhos(ixI^L, ixO^L, w, x)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L
    double precision             :: w(ixI^S, nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision             :: kin_en(ix0^S)

    if (hd_energy) then
       kin_en = 0.5d0 * sum(w(ixO^S, mom(:))**2, dim=^ND+1) / w(ixO^S, rho_)
       w(ixO^S, e_) = (hd_gamma - 1.0d0) * w(ixO^S, rho_)**(1.0d0 - hd_gamma) * &
            (w(ixO^S, e_) - kin_en)
    else
       call mpistop("energy from entropy can not be used with -eos = iso !")
    end if
  end subroutine e_to_rhos

  subroutine rhos_to_e(ixI^L, ixO^L, w, x)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L
    double precision             :: w(ixI^S, nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision             :: kin_en(ix0^S)

    if (hd_energy) then
       kin_en = 0.5d0 * sum(w(ixO^S, mom(:))**2, dim=^ND+1) / w(ixO^S, rho_)
       w(ixO^S, e_) = w(ixO^S, rho_)**(hd_gamma - 1.0d0) * w(ixO^S, e_) &
            / (hd_gamma - 1.0d0) + kin_en
    else
       call mpistop("entropy from energy can not be used with -eos = iso !")
    end if
  end subroutine rhos_to_e

  subroutine ppmflatcd(ixI^L, ixO^L, ixL^L, ixR^L, w, d2w, drho, dp)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L, ixL^L, ixR^L
    double precision, intent(in)    :: w(ixI^S, nw), d2w(ixG^T, 1:nwflux)
    double precision, intent(inout) :: drho(ixG^T), dp(ixG^T)

    if (hd_energy) then
       if (useprimitive) then
          drho(ixO^S) = hd_gamma*abs(d2w(ixO^S, rho_))&
               /min(w(ixL^S, rho_), w(ixR^S, rho_))
          dp(ixO^S) = abs(d2w(ixO^S, p_))/min(w(ixL^S, p_), w(ixR^S, p_))
       end if
    else
       call mpistop("PPM with flatcd=.true. can not be used with -eos = iso !")
    end if
  end subroutine ppmflatcd

  subroutine ppmflatsh(ixI^L, ixO^L, ixLL^L, ixL^L, ixR^L, ixRR^L, idims, w, drho, dp, dv)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L, ixLL^L, ixL^L, ixR^L, ixRR^L
    integer, intent(in)             :: idims
    double precision, intent(in)    :: w(ixI^S, nw)
    double precision, intent(inout) :: drho(ixG^T), dp(ixG^T), dv(ixG^T)
    double precision                :: v(ixG^T)

    if (hd_energy) then
       if (useprimitive) then
          ! eq. B15, page 218, Mignone and Bodo 2005, ApJS (beta1)
          where (abs(w(ixRR^S, p_)-w(ixLL^S, p_))>smalldouble)
             drho(ixO^S) = abs((w(ixR^S, p_)-w(ixL^S, p_))&
                  /(w(ixRR^S, p_)-w(ixLL^S, p_)))
          else where
             drho(ixO^S) = zero
          end where

          !  eq. B76, page 48, Miller and Collela 2002, JCP 183, 26
          !  use "dp" to save squared sound speed, assuming primitives
          dp(ixO^S) =(hd_gamma*w(ixO^S, p_)/w(ixO^S, rho_))

          dp(ixO^S) = abs(w(ixR^S, p_)-w(ixL^S, p_))&
               /(w(ixO^S, rho_)*dp(ixO^S))
          v(ixI^S)  = w(ixI^S, v0_+idims)
          call gradient(v, ixI^L, ixO^L, idims, dv)
       end if
    else
       call mpistop("PPM with flatsh=.true. can not be used with -eos = iso !")
    end if
  end subroutine ppmflatsh

  ! Calculate v_i = m_i/rho within ixO^L
  subroutine getv(w, x, ixI^L, ixO^L, idims, v)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L, idims
    double precision             :: w(ixI^S, nw), v(ixG^T)
    double precision, intent(in) :: x(ixI^S, 1:ndim)

    if (hd_energy .or. hd_adiab > 0.0d0) then
       v(ixO^S) = w(ixO^S, mom(idims))/w(ixO^S, rho_)
    else
       ! case of zero temperature: allow zero density
       where(w(ixO^S, rho_)/= zero)
          v(ixO^S) = w(ixO^S, mom(idims))/w(ixO^S, rho_)
       elsewhere
          v(ixO^S) = zero
       endwhere
    end if
  end subroutine getv

  ! Calculate cmax_idim = csound + abs(v_idim) within ixO^L
  subroutine getcmax(w, x, ixI^L, ixO^L, idims, cmax, cmin)
    use mod_global_parameters

    logical                      :: new_cmax, needcmin
    integer, intent(in)          :: ixI^L, ixO^L, idims
    double precision             :: w(ixI^S, nw), cmax(ixG^T), cmin(ixG^T)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision             :: csound(ixG^T){#IFDEF DUST , speeddust(ixG^T, 1:^NDS)}
    double precision             :: v(ixG^T)

    call getv(w, x, ixI^L, ixO^L, idims, v)

    if (hd_energy .or. hd_adiab > 0.0d0) then
       call getpthermal(w, x, ixI^L, ixO^L, csound)
       csound(ixO^S) = sqrt(hd_gamma*csound(ixO^S)/w(ixO^S, rho_))

       cmax(ixO^S)   = max(v(ix0^S)+csound(ixO^S), zero)
       cmin(ixO^S)   = min(v(ix0^S)-csound(ixO^S), zero)
    else
       ! case of zero temperature: allow zero density
       cmax(ixO^S) = max(v(ix0^S), zero)
       cmin(ixO^S) = min(v(ix0^S), zero)
    end if

    if (dust_num_species > 0) then
       call dust_get_cmax(todo)
    end if

  end subroutine getcmax

  ! Calculate thermal pressure=(gamma-1)*(e-0.5*m**2/rho) within ixO^L
  subroutine getpthermal(w, x, ixI^L, ixO^L, p)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L
    double precision             :: w(ixI^S, nw), p(ixG^T)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    integer, dimension(ixG^T)    :: patchierror
    integer, dimension(ndim)     :: lowpindex

    if (fixsmall) call smallvalues(w, x, ixI^L, ixO^L,"getpthermal")

    if (hd_energy) then
       p(ixO^S) = (hd_gamma - 1.0d0) * (w(ixO^S, e_) - &
            half * ({^C&w(ixO^S, mom(i_dir))**2+})/w(ixO^S, rho_))
    else
       p(ixO^S) = hd_adiab*w(ixO^S, rho_)**hd_gamma
    end if

    if (any(p(ixO^S)<minp)) then
       if (strictsmall) then
          lowpindex = minloc(p(ixO^S))
          ^D&lowpindex(^D) = lowpindex(^D)+ixOmin^D-1;
          write(*,*)'too small pressure = ', minval(p(ixO^S)),' with limit=', minp,&
               ' at x=', x(^D&lowpindex(^D), 1:ndim), lowpindex,' where E_k=',&
               half*(^C&w(^D&lowpindex(^D), mom(i_dir))**2+)/w(^D&lowpindex(^D), rho_)
          if (hd_energy) then
             write(*,*)' E_total=', w(^D&lowpindex(^D), e_)
          end if
          write(*,*)' w(1:nwflux) =', w(^D&lowpindex(^D), 1:nwflux),' when t=', t,' it=', it
          call mpistop("=== strictsmall in getpthermal ===")
       else
          if (strictgetaux) then
             where(p(ixO^S)<minp)
                p(ixO^S) = minp
                if (.not. hd_energy) then
                   w(ixO^S, rho_) = minrho
                   {^C&w(ixO^S, mom(i_dir)) = zero;}
                end if
             endwhere
          else
             where(p(ixO^S)<minp)
                patchierror(ixO^S) = 1
             else where
                patchierror(ixO^S) = 0
             end where
             if (any(patchierror(ixO^S)/= 0)) then
                call correctaux(ixI^L, ixO^L, w, x, patchierror,'getpthermal')
                where(patchierror(ixO^S)/= 0)
                   if (hd_energy) then
                      p(ixO^S) =(hd_gamma - 1.0d0)*(w(ixO^S, e_)- &
                           half*({^C&w(ixO^S, mom(i_dir))**2+})/w(ixO^S, rho_))
                   end if
                   if (.not. hd_energy) then
                      p(ixO^S) = hd_adiab*w(ixO^S, rho_)**hd_gamma
                   end if
                end where
             end if
          end if
       end if
    end if

  end subroutine getpthermal

  ! Calculate non-transport flux f_idim[iw] within ixO^L.
  subroutine getflux(w, x, ixI^L, ixO^L, idims, f, transport)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L, iw, idims
    double precision             :: w(ixI^S, nw), f(ixG^T), tmp(ixG^T)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    logical                      :: transport

    TODO: reorganize this, maybe compute all fluxes in a direction at once?

    transport = .true.

    if (iw == mom(idims)) then
       ! f_i[m_i]= v_i*m_i + p
       call getpthermal(w, x, ixI^L, ixO^L, f)
    else if (iw == e_) then
       ! f_i[e]= v_i*e + m_i/rho*p
       call getpthermal(w, x, ixI^L, ixO^L, f)
       f(ixO^S) = w(ixO^S, mom(idims))/w(ixO^S, rho_)*f(ixO^S)

    TODO: move this to dust module, see above comment
    else if (iw == rho_dust(i_dust)) then TODO
       where(w(ixO^S, rho_dust(i_dust))>minrhod)
          f(ixO^S) = w(ixO^S,(rho_dust(i_dust))+idims*^NDS);
       else where
          f(ixO^S) = zero
       end where
       transport = .false.
    else if (iw == m_dust(i_dim, i_dust)) then TODO
       ! use tmp for speeddust here
       where(w(ixO^S, rho_dust(i_dust))>minrhod)
          tmp(ixO^S) = w(ixO^S,(rho_dust(i_dust))+idims*^NDS)/w(ixO^S, rho_dust(i_dust));
       else where
          tmp(ixO^S) = zero;
       end where
       f(ixO^S) = w(ixO^S, iw)*tmp(ixO^S)
       transport = .false.
    else
       f(ixO^S) = zero
    endif
  end subroutine getflux

  ! Add geometrical source terms to w
  ! Notice that the expressions of the geometrical terms depend only on ndir,
  ! not ndim. Eg, they are the same in 2.5D and in 3D, for any geometry.
  ! Ileyk : to do :
  !     - give the possibility to set angmomfix=.true.
  !     - address the source term for the dust
  subroutine addgeometry(qdt, ixI^L, ixO^L, wCT, w, x) ! - - - - - - - - - - - -
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

    if (ndir==3) then

      select case (typeaxial)
      case ('slab')
        ! No source terms in slab symmetry
      case ('cylindrical')
        ! s[mr]=(pthermal+mphi**2/rho)/radius
        call getpthermal(wCT,x,ixI^L,ixO^L,tmp)
        tmp(ixO^S)=tmp(ixO^S)+wCT(ixO^S,mphi_)**2/wCT(ixO^S,rho_)
        w(ixO^S,mr_)=w(ixO^S,mr_)+qdt*tmp(ixO^S)/x(ixO^S,1)
        tmp(ixO^S)=zero
        ! s[mphi]=(-mphi*mr/rho)/radius
        ! Ileyk : beware the index permutation : mphi=2 if -phi=2 (2.5D (r,theta) grids) BUT mphi=3 if -phi=3 (for 2.5D (r,z) grids)
        if (.not. angmomfix) tmp(ixO^S)=-wCT(ixO^S,mphi_)*wCT(ixO^S,mr_)/wCT(ixO^S,rho_) ! no geometrical source term if angular momentum conserving form of the equations
        w(ixO^S,mphi_)=w(ixO^S,mphi_)+qdt*tmp(ixO^S)/x(ixO^S,1)
        tmp(ixO^S)=zero
      case ('spherical')
        h1x^L=ixO^L-kr(1,^D); h2x^L=ixO^L-kr(2,^D)
        ! s[mr]=((mtheta**2+mphi**2)/rho+2*p)/r
        call getpthermal(wCT,x,ixI^L,ixO^L,tmp)
        tmp(ixO^S)=tmp(ixO^S)*x(ixO^S,1) &
                  *(mygeo%surfaceC1(ixO^S)-mygeo%surfaceC1(h1x^S)) &
                  /mygeo%dvolume(ixO^S)
        tmp(ixO^S)=tmp(ixO^S)+wCT(ixO^S,mphi_  )**2/wCT(ixO^S,rho_)
      	tmp(ixO^S)=tmp(ixO^S)+wCT(ixO^S,mtheta_)**2/wCT(ixO^S,rho_)
      	w(ixO^S,mr_)=w(ixO^S,mr_)+qdt*tmp(ixO^S)/x(ixO^S,1)
        tmp(ixO^S)=zero
        ! s[mtheta]=-(mr*mtheta/rho)/r+cot(theta)*(mphi**2/rho+p)/r
        call getpthermal(wCT,x,ixI^L,ixO^L,tmp)
        tmp(ixO^S)=tmp(ixO^S)*x(ixO^S,1) &
                  *(mygeo%surfaceC2(ixO^S)-mygeo%surfaceC2(h2x^S)) &
                  /mygeo%dvolume(ixO^S)
        tmp(ixO^S)=tmp(ixO^S)+(wCT(ixO^S,mphi_)**2/wCT(ixO^S,rho_))/dtan(x(ixO^S,2))
        if (.not. angmomfix) tmp(ixO^S)=tmp(ixO^S)-(wCT(ixO^S,mtheta_)*wCT(ixO^S,mr_))/wCT(ixO^S,rho_)
        w(ixO^S,mtheta_)=w(ixO^S,mtheta_)+qdt*tmp(ixO^S)/x(ixO^S,1)
        tmp(ixO^S)=zero
        ! s[mphi]=-(mphi*mr/rho)/r-cot(theta)*(mtheta*mphi/rho)/r
        if (.not. angmomfix) tmp(ixO^S)=          -(wCT(ixO^S,mtheta_)*wCT(ixO^S,mphi_))/wCT(ixO^S,rho_)
        tmp(ixO^S)=tmp(ixO^S)/dtan(x(ixO^S,2))
        if (.not. angmomfix) tmp(ixO^S)=tmp(ixO^S)-(wCT(ixO^S,mphi_  )*wCT(ixO^S,mr_  ))/wCT(ixO^S,rho_)
        w(ixO^S,mphi_)=w(ixO^S,mphi_)+qdt*tmp(ixO^S)/x(ixO^S,1)
        tmp(ixO^S)=zero
      case default
          call mpistop("typeaxial is slab, cylindrical or spherical")
      end select

    elseif (ndir==2) then

      select case (typeaxial)
      case ('slab')
        ! No source terms in slab symmetry
      case ('cylindrical')
        ! (r,phi) : same as ndir==3
        ! phi true if and only if -d=22 -phi=2 (and typeaxial==cyl)
        if (phi) then ! Ileyk : new argument "phi" here for the parfile. Make sense just if typeaxial==cyl.
          ! s[mr]=(pthermal+mphi**2/rho)/radius
          call getpthermal(wCT,x,ixI^L,ixO^L,tmp)
          tmp(ixO^S)=tmp(ixO^S)+wCT(ixO^S,mphi_)**2/wCT(ixO^S,rho_)
          w(ixO^S,mr_)=w(ixO^S,mr_)+qdt*tmp(ixO^S)/x(ixO^S,1)
          tmp(ixO^S)=zero
          ! s[mphi]=(-mphi*mr/rho)/radius
          if (.not. angmomfix) tmp(ixO^S)=-wCT(ixO^S,mphi_)*wCT(ixO^S,mr_)/wCT(ixO^S,rho_)
          w(ixO^S,mphi_)=w(ixO^S,mphi_)+qdt*tmp(ixO^S)/x(ixO^S,1)
          tmp(ixO^S)=zero
        ! (r,z) : no mphi, just the pressure in the geom. source term
        else
          ! s[mr]=pthermal/radius
          call getpthermal(wCT,x,ixI^L,ixO^L,tmp)
          w(ixO^S,mr_)=w(ixO^S,mr_)+qdt*tmp(ixO^S)/x(ixO^S,1)
          tmp(ixO^S)=zero
        endif
      case ('spherical') ! (r,theta), w/ theta the colatitude. No mphi
        ! s[mr]=((mtheta**2)/rho+2*p)/r
        call getpthermal(wCT,x,ixI^L,ixO^L,tmp)
        tmp(ixO^S)=tmp(ixO^S)*x(ixO^S,1) &
                  *(mygeo%surfaceC1(ixO^S)-mygeo%surfaceC1(h1x^S)) &
                  /mygeo%dvolume(ixO^S)
      	tmp(ixO^S)=tmp(ixO^S)+wCT(ixO^S,mtheta_)**2/wCT(ixO^S,rho_)
      	w(ixO^S,mr_)=w(ixO^S,mr_)+qdt*tmp(ixO^S)/x(ixO^S,1)
        tmp(ixO^S)=zero
        ! s[mtheta]=-(mr*mtheta/rho)/r+cot(theta)*p/r
        call getpthermal(wCT,x,ixI^L,ixO^L,tmp)
        tmp(ixO^S)=tmp(ixO^S)*x(ixO^S,1) &
                  *(mygeo%surfaceC2(ixO^S)-mygeo%surfaceC2(h2x^S)) &
                  /mygeo%dvolume(ixO^S)
        if (.not. angmomfix) tmp(ixO^S)=tmp(ixO^S)-(wCT(ixO^S,mtheta_)*wCT(ixO^S,mr_))/wCT(ixO^S,rho_)
        w(ixO^S,mtheta_)=w(ixO^S,mtheta_)+qdt*tmp(ixO^S)/x(ixO^S,1)
        tmp(ixO^S)=zero
      case default
          call mpistop("typeaxial is slab, cylindrical or spherical")
      end select

    elseif (ndir==1) then

      select case (typeaxial)
      case ('slab')
        ! No source terms in slab symmetry
      case ('cylindrical')
        ! s[mr]=pthermal/radius
        call getpthermal(wCT,x,ixI^L,ixO^L,tmp)
        w(ixO^S,mr_)=w(ixO^S,mr_)+qdt*tmp(ixO^S)/x(ixO^S,1)
        tmp(ixO^S)=zero
      case ('spherical')
        ! s[mr]=2pthermal/radius
        call getpthermal(wCT,x,ixI^L,ixO^L,tmp)
        tmp(ixO^S)=tmp(ixO^S)*x(ixO^S,1) &
                  *(mygeo%surfaceC1(ixO^S)-mygeo%surfaceC1(h1x^S)) &
                  /mygeo%dvolume(ixO^S)
        w(ixO^S,mr_)=w(ixO^S,mr_)+qdt*tmp(ixO^S)/x(ixO^S,1)
        tmp(ixO^S)=zero
      case default
          call mpistop("typeaxial is slab, cylindrical or spherical")
      end select

    endif

    if (fixsmall) call smallvalues(w, x, ixI^L, ixO^L,"addgeometry")

  end subroutine addgeometry ! - - - - - - - - - - - - - - - - - - - - - - - - -

  ! w[iw]= w[iw]+qdt*S[wCT, qtC, x] where S is the source based on wCT within ixO
  subroutine addsource(qdt, ixI^L, ixO^L, iw^LIM, qtC, wCT, qt, w, x, qsourcesplit)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L, iw^LIM
    double precision, intent(in)    :: qdt, qtC, qt, x(ixI^S, 1:ndim)
    double precision, intent(inout) :: wCT(ixI^S, 1:nw), w(ixI^S, 1:nw)
    logical, intent(in)             :: qsourcesplit

    double precision, dimension(ixG^T, 1:ndir, 1:^NDS) :: fdrag
    integer                                           :: idir

    if (dust_num_species > 0) then
       call dust_add_source(qdt, ixI^L, ixO^L, iw^LIM, &
            qtC, wCT, qt, w, x, qsourcesplit)
    end if

  end subroutine addsource

  subroutine getdt(w, ixI^L, ixO^L, dtnew, dx^D, x)
    use mod_global_parameters

    integer, intent(in)                        :: ixI^L, ixO^L
    double precision, intent(in)               :: dx^D, x(ixI^S, 1:ndim)
    double precision, intent(inout)            :: w(ixI^S, 1:nw), dtnew

    integer                                    :: idims, idust, idir
    double precision, dimension(1:^NDS)        :: dtdust
    double precision, dimension(ixG^T)         :: vt2, deltav, tstop, ptherm, vdust, vgas
    double precision, dimension(ixG^T, 1:^NDS) :: alpha_T
    double precision                           :: K

    if (dust_num_species > 0) then
       call dust_get_dt(w, ixI^L, ixO^L, dtnew, dx^D, x)
    end if

  end subroutine getdt

end module mod_hd
