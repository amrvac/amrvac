module mod_hd
  use mod_physics

  implicit none
  private

  logical :: hd_energy
  integer :: hd_n_dust
  integer :: hd_n_tracer

  integer, parameter :: rho_ = 1
  integer, parameter :: m0_  = rho_
  integer, allocatable :: mom(:)
  integer :: e_

  integer, allocatable :: rho_dust(:)
  integer, allocatable :: m_dust(:, :)

  double precision :: gamma

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
    ! Todo: document
    nwflux       = ^NC * (1 + hd_n_dust) + 2 + hd_n_dust + hd_n_tracer
    if (.not. hd_energy) nwflux = nwflux - 1

    nwaux   = 0
    nwextra = 0
    nw      = nwflux + nwaux + nwextra

    ! Note: dust related velocity vectors not handled here
    nvector      = 1 ! No. vector vars
    allocate(iw_vector(nvector))
    iw_vector(1) = m0_

    nflag_ = nw + 1

    allocate(mom(^NC))
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

    ! Set starting index of dust species
    do n = 1, hd_n_dust
       TODO
       dust_var(n, ) = ix + 1
       do j = 1, ^NC
          dust_mdens(n, j) = ix + 1 + j
       end do
       ix = ix + 1 + ^NC
    end do

    if (hd_n_dust > 0) then
       dust_ = ix + 1
       ix    = ix + hd_n_dust * (1 + ^NC)
    else
       dust_ = -1
    end if

    mr_= m0_+r_, mphi_= m0_+phi_, mz_= m0_+z_  ! Polar var. names

    if (hd_n_dust > 0) then
       TODO
       v1d^DS_= rhod^NDS_+^DS
       {^NOONEC
       v2d^DS_= m1d^NDS_+^DS }
       {^IFTHREEC
       v3d^DS_= m2d^NDS_+^DS }
       mrd^DS_= m0d^DS_+(^NDS*r_)              ! Polar var. names
       mphid^DS_= m0d^DS_+(^NDS*phi_)          ! Polar var. names
       mzd^DS_= m0d^DS_+(^NDS*z_)              ! Polar var. names
    end if

    if (hd_energy) then
       ! Characteristic waves
       soundRW_= 1
       soundLW_= 2
       entropW_= 3
       shearW0_= 3
       nworkroe = 3
    else
       ! Characteristic waves
       soundRW_= 1
       soundLW_= 2
       shearW0_= 2
       nworkroe = 1
    end if

  end subroutine hd_activate

  subroutine checkglobaldata
    use mod_global_parameters

    if (.not. hd_energy) then
       if (hd_gamma <= zero) call mpistop ("gamma negative not ok")
       if (hd_adiab < zero) call mpistop ("adiab strict negative not ok")
       minrho= max(zero, smallrho)
       minp = hd_adiab*minrho**hd_gamma
    else
       if (hd_gamma <= zero .or. hd_gamma== one) call mpistop ("gamma negative or 1 not ok")
       minp  = max(zero, smallp)
       minrho= max(zero, smallrho)
       smalle= minp/(hd_gamma-one)
    end if

    if (hd_n_dust > 0) then
       if (eqpar(mu_)<= zero) call mpistop ("mu (molecular weight) negative not ok")
       minrhod= max(zero, smallrhod)
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

    if (hd_n_dust > 0) then
       eqpar(mu_) = one
       mhcgspar = 1.6733D-24
       kbcgspar = 1.38065D-16
    end if

    do il = 1, nw
       select case(il)
       case(soundRW_, soundLW_)
          entropycoef(il) = 0.2d0
       case default
          entropycoef(il) = -one
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
          tmp(ixO^S) =(hd_gamma-one)*(w(ixO^S, e_)- &
               half*( ^C&w(ixO^S, mom(i_dir))**2+ )/w(ixO^S, rho_))
          flag(ixO^S) =(tmp(ixO^S)>= minp .and. w(ixO^S, rho_)>= minrho)
       endif
    else
       if (hd_adiab /= zero) then
          flag(ixO^S) = (w(ixO^S, rho_) >= zero)
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

    invgam = 1.d0/(hd_gamma-one)

    if (hd_energy) then
       ! Calculate total energy from pressure and kinetic energy
       w(ixO^S, e_) = w(ixO^S, p_)*invgam+ &
            half*w(ixO^S, rho_)*(^C&w(ixO^S, v^C_)**2+)
    end if

    ! Convert velocity to momentum
    ^C&w(ixO^S, mom(i_dir)) = w(ixO^S, rho_)*w(ixO^S, v^C_);

    if (hd_n_tracer > 0) then
       {^FL&w(ixO^S, tracer(i_tr)) = w(ixO^S, rho_)*w(ixO^S, tracer(i_tr))\}
    end if

    if (hd_n_dust > 0) then
       ! Convert dust velocity to dust momentum
       {^DS&{^C&w(ixO^S, m_dust(i_dim, i_dust)) = w(ixO^S, rho_dust(i_dust))*w(ixO^S, v^Cd^DS_);}\}
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

    if (fixsmall) call smallvalues(w, x, ixI^L, ixO^L,"primitive")

    if (hd_energy) then
       ! compute pressure
       w(ixO^S, p_) =(hd_gamma-one)*(w(ixO^S, e_)- &
            half*( ^C&w(ixO^S, mom(i_dir))**2+ )/w(ixO^S, rho_))
    end if

    ! Convert dust momentum to dust velocity
    do i_dust = 1, hd_n_dust
       where(w(ixO^S, rho_dust(i_dust))>minrhod)
          do i_dir = 1, ^NC
          w(ixO^S, m_dust(i_dir, i_dust)) = w(ixO^S, m_dust(i_dir, i_dust)) /&
               w(ixO^S, rho_dust(i_dust))
          end do
       elsewhere
          w(ixO^S, m_dust(:, i_dust)) = 0.0d0
       end where
    end if

    if (hd_energy) then
       ! Convert momentum to velocity
       do i_dir = 1, ^NC
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
       else
          if (strictgetaux) then
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
    end if

    if (.not. hd_energy) then

       if (hd_adiab > zero) then
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

  end subroutine primitive

  subroutine e_to_rhos(ixI^L, ixO^L, w, x)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L
    double precision             :: w(ixI^S, nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision             :: sum_m2(ix0^S)

    if (hd_energy) then
       sum_m2 = sum(w(ixO^S, mom(:))**2, dim=^ND+1)
       w(ixO^S, e_) = (hd_gamma-one) * w(ixO^S, rho_)**(one-hd_gamma) * &
            (w(ixO^S, e_) - half * sum_m2 / w(ixO^S, rho_))
    else
       call mpistop("energy from entropy can not be used with -eos = iso !")
    end if
  end subroutine e_to_rhos

  subroutine rhos_to_e(ixI^L, ixO^L, w, x)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L
    double precision             :: w(ixI^S, nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision             :: sum_m2(ix0^S)

    if (hd_energy) then
       sum_m2 = sum(w(ixO^S, mom(:))**2, dim=^ND+1)
       w(ixO^S, e_) = w(ixO^S, rho_)**(hd_gamma-one) * w(ixO^S, rhos_) &
            / (hd_gamma - one) + half * sum_m2 / w(ixO^S, rho_)
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
          drho(ixO^S) = hd_gamma*dabs(d2w(ixO^S, rho_))&
               /min(w(ixL^S, rho_), w(ixR^S, rho_))
          dp(ixO^S) = dabs(d2w(ixO^S, p_))/min(w(ixL^S, p_), w(ixR^S, p_))
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
          where (dabs(w(ixRR^S, p_)-w(ixLL^S, p_))>smalldouble)
             drho(ixO^S) = dabs((w(ixR^S, p_)-w(ixL^S, p_))&
                  /(w(ixRR^S, p_)-w(ixLL^S, p_)))
          else where
             drho(ixO^S) = zero
          end where

          !  eq. B76, page 48, Miller and Collela 2002, JCP 183, 26
          !  use "dp" to save squared sound speed, assuming primitives
          dp(ixO^S) =(hd_gamma*w(ixO^S, p_)/w(ixO^S, rho_))

          dp(ixO^S) = dabs(w(ixR^S, p_)-w(ixL^S, p_))&
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

    if (hd_energy) then
       v(ixO^S) = w(ixO^S, mom(idims))/w(ixO^S, rho_)
    else
       if (hd_adiab > zero) then
          v(ixO^S) = w(ixO^S, mom(idims))/w(ixO^S, rho_)
       else
          ! case of zero temperature: allow zero density
          where(w(ixO^S, rho_)/= zero)
             v(ixO^S) = w(ixO^S, mom(idims))/w(ixO^S, rho_)
          elsewhere
             v(ixO^S) = zero
          endwhere
       endif
    end if
  end subroutine getv

  ! Calculate cmax_idim = csound + abs(v_idim) within ixO^L
  subroutine getcmax(new_cmax, w, x, ixI^L, ixO^L, idims, cmax, cmin, needcmin)
    use mod_global_parameters

    logical                      :: new_cmax, needcmin
    integer, intent(in)          :: ixI^L, ixO^L, idims
    double precision             :: w(ixI^S, nw), cmax(ixG^T), cmin(ixG^T)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision             :: csound(ixG^T){#IFDEF DUST , speeddust(ixG^T, 1:^NDS)}

    if (hd_energy) then
       call getpthermal(w, x, ixI^L, ixO^L, csound)
       csound(ixO^S) = sqrt(hd_gamma*csound(ixO^S)/w(ixO^S, rho_))
       if (needcmin) then
          cmax(ixO^S) = max(w(ixO^S, mom(idims))/w(ixO^S, rho_)+csound(ixO^S), &
               zero)
          cmin(ixO^S) = min(w(ixO^S, mom(idims))/w(ixO^S, rho_)-csound(ixO^S),&
               zero)
       else
          cmax(ixO^S) = max(csound(ixO^S)+dabs(w(ixO^S, mom(idims))/w(ixO^S, rho_)), &
               zero)
       endif
    else
       if (hd_adiab > zero) then
          call getpthermal(w, x, ixI^L, ixO^L, csound)
          csound(ixO^S) = sqrt(hd_gamma*csound(ixO^S)/w(ixO^S, rho_))
          if (needcmin) then
             cmax(ixO^S) = max(w(ixO^S, mom(idims))/w(ixO^S, rho_)+csound(ixO^S),&
                  zero)
             cmin(ixO^S) = min(w(ixO^S, mom(idims))/w(ixO^S, rho_)-csound(ixO^S),&
                  zero)
          else
             cmax(ixO^S) = max(csound(ixO^S)+abs(w(ixO^S, mom(idims))/w(ixO^S, rho_)), &
                  zero)
          endif
       else
          ! case of zero temperature: allow zero density
          if (needcmin) then
             where(w(ixO^S, rho_)/= zero)
                cmax(ixO^S) = max(w(ixO^S, mom(idims))/w(ixO^S, rho_), zero)
                cmin(ixO^S) = min(w(ixO^S, mom(idims))/w(ixO^S, rho_), zero)
             elsewhere
                cmax(ixO^S) = 0.0d0
                cmin(ixO^S) = 0.0d0
             endwhere
          else
             where(w(ixO^S, rho_)/= zero)
                cmax(ixO^S) = abs(w(ixO^S, mom(idims))/w(ixO^S, rho_))
             elsewhere
                cmax(ixO^S) = 0.0d0
             endwhere
          endif
       endif
    end if

    ! TODO: update using dust speed (call subroutine)
    do i_dust = 1, hd_n_dust
       where(w(ixO^S, rho_dust(i_dust)) > minrhod)
          speeddust(ixO^S, i_dust) = w(ixO^S, m_dust(i_dims, i_dust)) / &
               w(ixO^S, rho_dust(i_dust));
       elsewhere
          speeddust(ixO^S, i_dust) = zero;
       end where
    end do

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
       p(ixO^S) =(hd_gamma-one)*(w(ixO^S, e_)- &
            half*({^C&w(ixO^S, mom(i_dir))**2+})/w(ixO^S, rho_))
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
                      p(ixO^S) =(hd_gamma-one)*(w(ixO^S, e_)- &
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
  subroutine getflux(w, x, ixI^L, ixO^L, iw, idims, f, transport)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L, iw, idims
    double precision             :: w(ixI^S, nw), f(ixG^T), tmp(ixG^T)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    logical                      :: transport

    transport = .true.

    if (iw == mom(idims)) then
       ! f_i[m_i]= v_i*m_i + p
       call getpthermal(w, x, ixI^L, ixO^L, f)
    else if (iw == e_) then
       ! f_i[e]= v_i*e + m_i/rho*p
       call getpthermal(w, x, ixI^L, ixO^L, f)
       f(ixO^S) = w(ixO^S, mom(idims))/w(ixO^S, rho_)*f(ixO^S)
    else if (iw == rho_dust(i_dust)) then
       where(w(ixO^S, rho_dust(i_dust))>minrhod)
          f(ixO^S) = w(ixO^S,(rho_dust(i_dust))+idims*^NDS);
       else where
          f(ixO^S) = zero
       end where
       transport = .false.
    else if (iw == m_dust(i_dim, i_dust)) then
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
  subroutine addgeometry(qdt, ixI^L, ixO^L, wCT, w, x)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, x(ixI^S, 1:ndim)
    double precision, intent(inout) :: wCT(ixI^S, 1:nw), w(ixI^S, 1:nw)

    select case (typeaxial)
    case ('slab')
       ! No source terms in slab symmetry
    case default
       call mpistop("addgeometry not yet implemented")
    end select

    if (fixsmall) call smallvalues(w, x, ixI^L, ixO^L,"addgeometry")

  end subroutine addgeometry

  ! w[iw]= w[iw]+qdt*S[wCT, qtC, x] where S is the source based on wCT within ixO
  subroutine addsource(qdt, ixI^L, ixO^L, iw^LIM, qtC, wCT, qt, w, x, qsourcesplit)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L, iw^LIM
    double precision, intent(in)    :: qdt, qtC, qt, x(ixI^S, 1:ndim)
    double precision, intent(inout) :: wCT(ixI^S, 1:nw), w(ixI^S, 1:nw)
    logical, intent(in)             :: qsourcesplit

    double precision, dimension(ixG^T, 1:^NC, 1:^NDS) :: fdrag
    integer                                           :: idir

    if (hd_n_dust > 0) then
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

    dtnew = bigdouble
    dtdust = bigdouble

    !-----------------------------
    !get dt related to dust and gas stopping time (Laibe 2011)
    !-----------------------------

    select case( TRIM(dustmethod) )

    case( 'Kwok' ) ! assume sticking coefficient equals 0.25
       dtdust(1:^NDS) = bigdouble

       call getpthermal(w, x, ixI^L, ixO^L, ptherm)
       vt2(ixO^S) = 3.0d0*ptherm(ixO^S)/w(ixO^S, rho_)

       ! Tgas, mu = mean molecular weight
       ptherm(ixO^S) = ( ptherm(ixO^S)*normvar(p_)*mhcgspar*eqpar(mu_))/(w(ixO^S, rho_)*normvar(rho_)*kbcgspar)

       do idir = 1,^NC
          call getv(w, x, ixI^L, ixO^L, idir, vgas)

          TODO
          where(w(ixO^S, rho_dust(i_dust))>minrhod)
             vdust(ixO^S)  = w(ixO^S,(rho_dust(i_dust))+idir*^NDS)/w(ixO^S, rho_dust(i_dust))
             deltav(ixO^S) = (vgas(ixO^S)-vdust(ixO^S))
             tstop(ixO^S)  = 4.0d0*(rhodust(^DS)*sdust(^DS))/ &
                  (3.0d0*(0.75d0)*dsqrt(vt2(ixO^S) + &
                  deltav(ixO^S)**2)*(w(ixO^S, rho_dust(i_dust)) + &
                  w(ixO^S, rho_)))
          else where
             tstop(ixO^S) = bigdouble
          end where


          dtdust(^DS) = min(minval(tstop(ixO^S)), dtdust(^DS))
       enddo

       dtnew = min(minval(dtdiffpar*dtdust(1:^NDS)), dtnew)

    case( 'sticking' ) ! Calculate sticking coefficient based on the gas temperature
       dtdust(1:^NDS) = bigdouble

       call getpthermal(w, x, ixI^L, ixO^L, ptherm)
       vt2(ixO^S) = 3.0d0*ptherm(ixO^S)/w(ixO^S, rho_)


       ! Sticking coefficient
       call get_sticking(w, x, ixI^L, ixO^L, alpha_T )

       ! Tgas, mu = mean molecular weight
       ptherm(ixO^S) = ( ptherm(ixO^S)*normvar(p_)*mhcgspar*eqpar(mu_))/(w(ixO^S, rho_)*normvar(rho_)*kbcgspar)



       do idir = 1,^NC
          call getv(w, x, ixI^L, ixO^L, idir, vgas)

          TODO
          where(w(ixO^S, rho_dust(i_dust))>minrhod)
             vdust(ixO^S)  = w(ixO^S,(rho_dust(i_dust))+idir*^NDS)/w(ixO^S, rho_dust(i_dust))
             deltav(ixO^S) = (vgas(ixO^S)-vdust(ixO^S))
             tstop(ixO^S)  = 4.0d0*(rhodust(^DS)*sdust(^DS))/ &
                  (3.0d0*(one-alpha_T(ixO^S,^DS))*dsqrt(vt2(ixO^S) + &
                  deltav(ixO^S)**2)*(w(ixO^S, rho_dust(i_dust)) + &
                  w(ixO^S, rho_)))
          else where
             tstop(ixO^S) = bigdouble
          end where


          dtdust(^DS) = min(minval(tstop(ixO^S)), dtdust(^DS))
          \}
       enddo

       dtnew = min(minval(dtdiffpar*dtdust(1:^NDS)), dtnew)



    case('linear') !linear with Deltav, for testing (see Laibe & Price 2011)
       K = 3.4d5/^NDS
       dtdust(1:^NDS) = bigdouble

       TODO
       where(w(ixO^S, rho_dust(i_dust))>minrhod)
          tstop(ixO^S)  = (w(ixO^S, rho_dust(i_dust))*w(ixO^S, rho_))/ &
               (K*(w(ixO^S, rho_dust(i_dust)) + w(ixO^S, rho_)))
       else where
          tstop(ixO^S) = bigdouble
       end where


       dtdust(^DS) = min(minval(tstop(ixO^S)), dtdust(^DS))
       \}

       dtnew = min(minval(dtdiffpar*dtdust(1:^NDS)), dtnew)
    case('none')
       ! no dust timestep
    case default
       call mpistop( "=== This dust method has not been implemented===" )
    end select


    if (dtnew < dtmin) then
       write(unitterm,*)"-------------------------------------"
       write(unitterm,*)"Warning: found DUST related time step too small! dtnew=", dtnew
       write(unitterm,*)"on grid with index:", saveigrid," grid level=", node(plevel_, saveigrid)
       write(unitterm,*)"grid corners are=",{^D&rnode(rpxmin^D_, saveigrid), rnode(rpxmax^D_, saveigrid)}
       write(unitterm,*)" dtdust =", dtdust(1:^NDS)
       write(unitterm,*)"on processor:", mype
       write(unitterm,*)"-------------------------------------"
    endif

  end subroutine getdt
}

end module mod_hd
