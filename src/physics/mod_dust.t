module mod_dust

  implicit none
  private

  integer, public, protected :: n_dust
  integer, allocatable :: rho_dust(:)
  integer, allocatable :: m_dust(:, :)

  double precision, allocatable :: sdust(:), dsdust(:), rhodust(:), mhcgspar, kbcgspar
  double precision              :: Lstar, Tdust

  ! Public methods
  public :: dust_init

contains

  subroutine dust_init()
    use mod_global_parameters

    integer :: n

    allocate(rho_dust(n_dust))
    allocate(m_dust(ndir, n_dust))

    ! Set starting index of dust species
    do n = 1, n_dust
       nwflux = nwflux + 1
       rho_dust(n) = nwflux     ! Dust density

       do idir = 1, ndir
          nwflux = nwflux + 1
          m_dust(idir, n) = nwflux ! Dust momentum
       end do
    end do

    ! set default values for entropy fixes for 'yee' type
    ! subroutine initglobaldata

    !   if (hd_energy) then
    !      hd_gamma = 5.d0/3.d0
    !   else
    !      hd_gamma = 1.d0
    !      hd_adiab = 1.d0
    !   end if

    !   ! if (dust_num_species > 0) then
    !   !    eqpar(mu_) = 1.0d0
    !   !    mhcgspar = 1.6733D-24
    !   !    kbcgspar = 1.38065D-16
    !   ! end if

    ! end subroutine initglobaldata

  end subroutine dust_init

  subroutine dust_to_conserved(ixI^L, ixO^L, w, x)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    integer                         :: n, idir

    do n = 1, n_dust
       ! Convert velocity to momentum
       do idir = 1, ndir
          w(ixO^S, m_dust(i_dim, n)) = w(ixO^S, rho_dust(n)) * &
               w(ixO^S, m_dust(i_dim, n))
       end do
    end do
  end subroutine dust_to_conserved

  subroutine dust_to_primitive(ixI^L, ixO^L, w, x)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    integer                         :: n, idir

    do n = 1, n_dust
       ! Convert momentum to velocity
       do idir = 1, ndir
          where (w(ixO^S, rho_dust(n)) > minrhod)
             w(ixO^S, m_dust(i_dim, n)) = w(ixO^S, rho_dust(n)) / &
                  w(ixO^S, m_dust(i_dim, n))
          elsewhere
             w(ixO^S, m_dust(i_dim, n)) = 0
          end where
       end do
    end do
  end subroutine dust_to_primitive

  subroutine dust_get_flux(w, x, ixI^L, ixO^L, iw, idim, f, transport)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L, iw, idim
    double precision, intent(in)    :: w(ixI^S, 1:nw), x(ixI^S, 1:^ND)
    double precision, intent(inout) :: f(ixG^T)
    double precision                :: vdust(ixG^T)
    logical, intent(out)            :: transport
    integer                         :: n

    transport = .false.
    f(ixO^S) = 0.0d0

    do n = 1, n_dust
       if (iw == rho_dust(n)) then ! A dust density
          where (w(ixO^S, rho_dust(n)) > minrhod)
             f(ixO^S) = w(ixO^S, m_dust(idim, n))
          elsewhere             ! TODO: remove?
             f(ixO^S) = 0.0d0
          end where
          exit
       else if (iw == m_dust(idim, n)) then ! A dust momentum component
          f(ixO^S) = w(ixO^S, iw) * get_vdust(w, ixI^L, ixO^L, idim, n)
          exit
       end if
    end do
  end subroutine dust_get_flux

  function get_vdust(w, ixI^L, ixO^L, idim, n) result(vdust)
    use mod_global_parameters, only: nw
    integer, intent(in)           :: ixI^L, ixO^L, idim, n
    double precision, intent(in)  :: w(ixI^S, nw)
    double precision              :: vdust(ixO^S)

    where (w(ixO^S, rho_dust(n)) > minrhod)
       vdust = w(ixO^S, m_dust(idim, n)) / w(ixO^S, rho_dust(n))
    elsewhere
       vdust = 0.0d0;
    end where
  end function get_vdust

  ! Force dust density to zero if rho_dust <= minrhod
  subroutine set_dusttozero(qdt, ixI^L, ixO^L, iw^LIM, qtC, wCT, qt, w, x)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L, iw^LIM
    double precision, intent(in)    :: qdt, qtC, qt, x(ixI^S, 1:ndim)
    double precision, intent(inout) :: wCT(ixI^S, 1:nw), w(ixI^S, 1:nw)
    integer                         :: n, idir

    do n = 1, n_dust
       where (w(ixO^S, rho_dust(n)) <= minrhod)
          w(ixO^S, rho_dust(n)) = 0.0d0
       end where

       do idir = 1, ndir
          where (w(ixO^S, rho_dust(n)) <= minrhod)
             w(ixO^S, m_dust(idir, n)) = 0.0d0
          end where
       end do
    end do
  end subroutine set_dusttozero

  ! Calculate drag force based on Epstein's or Stokes' law
  ! From Kwok 1975, page 584 (between eqn 8 and 9)
  subroutine get_3d_dragforce(ixI^L, ixO^L, w, x, fdrag)
    use mod_global_parameters

    integer, intent(in)                                            :: ixI^L, ixO^L
    double precision, intent(in)                                   :: x(ixI^S, 1:ndim)
    double precision, intent(inout)                                :: w(ixI^S, 1:nw)
    double precision, intent(out), dimension(ixG^T, 1:^NC, 1:n_dust) :: fdrag

    double precision, dimension(ixG^T)           :: vt2, deltav, fd, ptherm, vdust, vgas, Tgas
    double precision, dimension(ixG^T, 1:n_dust) :: alpha_T
    integer                                      :: n, idir
    double precision                             :: K

    ! TODO: include pthermal in generic physics?
    call hd_get_pthermal(w, x, ixI^L, ixO^L, ptherm)

    vt2(ixO^S) = 3.0d0*ptherm(ixO^S)/w(ixO^S, rho_)

    select case( TRIM(dustmethod) )
    case ('Kwok') ! assume sticking coefficient equals 0.25

       do idir = 1, ndir
          call hd_get_v(w, x, ixI^L, ixO^L, idir, vgas)

          do n = 1, n_dust
             ! TODO: simplify
             where(w(ixO^S, rho_dust(n)) > minrhod)
                vdust(ixO^S)  = w(ixO^S, m_dust(idir, n)) / w(ixO^S, rho_dust(n))
                deltav(ixO^S) = (vgas(ixO^S)-vdust(ixO^S))

                ! 0.75 from sticking coefficient
                fd(ixO^S)     = 0.75d0*w(ixO^S, rho_dust(n))*w(ixO^S, rho_)*deltav(ixO^S) &
                     / (rhodust(n) * sdust(n))

                ! 0.75 from spherical grainvolume
                fd(ixO^S)     = -fd(ixO^S)*0.75d0*dsqrt(vt2(ixO^S) + deltav(ixO^S)**2)
             elsewhere
                fd(ixO^S) = 0.0d0
             end where
             fdrag(ixO^S, idir, n) = fd(ixO^S)
          end do
       end do

    case ('sticking') ! Calculate sticking coefficient based on the gas and dust temperatures
       !  Equation from Decin et al. 2006
       if (.not. hd_energy) call mpistop("dust sticking requires hd_energy")

       Tgas(ixO^S) = ( ptherm(ixO^S)*normvar(e_)*mhcgspar) / &
            (w(ixO^S, rho_)*normvar(rho_)*kbcgspar)
       call get_sticking(w, x, ixI^L, ixO^L, alpha_T)

       do idir = 1, ndir
          call hd_get_v(w, x, ixI^L, ixO^L, idir, vgas)

          do n = 1, n_dust
             where(w(ixO^S, rho_dust(n))>minrhod)
                vdust(ixO^S)  = w(ixO^S,m_dust(idir, n)) / w(ixO^S, rho_dust(n))
                deltav(ixO^S) = (vgas(ixO^S)-vdust(ixO^S))
                fd(ixO^S)     = (one-alpha_T(ixO^S,n)) * w(ixO^S, rho_dust(n))*w(ixO^S, rho_) * &
                     deltav(ixO^S) / (rhodust(n)*sdust(n))
                fd(ixO^S)     = -fd(ixO^S) * 0.75d0 * dsqrt(vt2(ixO^S) + deltav(ixO^S)**2)
             else where
                fd(ixO^S) = 0.0d0
             end where
             fdrag(ixO^S, idir,n) = fd(ixO^S)
          end do
       end do
    case('linear') !linear with Deltav, for testing (see Laibe & Price 2011)
       K = 3.4d5 / n_dust
       do idir = 1, ndir
          call hd_get_v(w, x, ixI^L, ixO^L, idir, vgas)

          do n = 1, n_dust
             where(w(ixO^S, rho_dust(n))>minrhod)
                vdust(ixO^S)  = w(ixO^S,m_dust(idir, n))/w(ixO^S, rho_dust(n))
                deltav(ixO^S) = (vgas(ixO^S)-vdust(ixO^S))

                fd(ixO^S)     = -K*deltav(ixO^S)
             else where
                fd(ixO^S) = 0.0d0
             end where
             fdrag(ixO^S, idir,n) = fd(ixO^S)
          end do
       end do
    case('none')
       fdrag(ixO^S, idir, :) = 0.0d0
    case default
       call mpistop( "=== This dust method has not been implemented===" )
    end select

  end subroutine get_3d_dragforce

  subroutine get_sticking(w, x, ixI^L, ixO^L, alpha_T)
    !
    !  get sticking coefficient
    !
    !  Assume cgs units, and use of normvar(0:nw) array for conversion
    !
    !
    !  Equation from Decin et al. 2006
    !
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    double precision, intent(inout) :: w(ixI^S, 1:nw)


    double precision, intent(out)   :: alpha_T(ixG^T, 1:n_dust)

    double precision                :: Tgas(ixG^T)
    !-----------------------------------------------------------------------------

    call getpthermal(w, x, ixI^L, ixO^L, Tgas)
    call get_tdust(w, x, ixI^L, ixO^L, alpha_T)

    Tgas(ixO^S) = (Tgas(ixO^S)*normvar(e_)*mhcgspar) / &
         (w(ixO^S, rho_) * normvar(rho_) * kbcgspar)

    do n = 1, n_dust
       alpha_T(ixO^S,n) =  max(0.35d0 * exp(-dsqrt((Tgas(ixO^S) + &
            alpha_T(ixO^S,n))/5.0d2))+0.1d0, smalldouble)
    end do
  end subroutine get_sticking

  ! Returns dust temperature (in K), either as constant or based
  ! on equ. 5.41, 5.42 and 5.44 from Tielens (2005)
  !
  !  Note that this calculation assumes cgs!!!! with conversion between
  !  physical and scaled quantities done through the normvar(0:nw) array!!!!
  !
  !  It takes as input the stellar luminosoity in solar units
  !  and/or a fixed dust temperature in Kelvin
  subroutine get_tdust(w, x, ixI^L, ixO^L, Td)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    double precision, intent(inout) :: w(ixI^S, 1:nw)
    double precision, intent(out)   :: Td(ixG^T, 1:n_dust)
    double precision                :: G0(ixO^S)

    select case( TRIM(dusttemp) )
    case( 'constant' )
       Td(ixO^S, :) = Tdust
    case( 'ism' )
       select case( trim(dustspecies) )
       case( 'graphite' )
          Td(ixO^S, :) = 15.8d0*((0.0001d0/(sdust(:)*normvar(0)))**0.06d0)
       case( 'silicate' )
          Td(ixO^S, :) = 13.6d0*((0.0001d0/(sdust(:)*normvar(0)))**0.06d0)
       case default
          call mpistop( "=== Dust species undetermined===" )
       end select
    case( 'stellar' )
       select case( trim(typeaxial) )
       case( 'spherical' )
          G0(ixO^S) = max(x(ixO^S, 1)*normvar(0), smalldouble)
       case( 'cylindrical' )
          G0(ixO^S) = max(dsqrt(x(ixO^S, 1)**2 + x(ixO^S, 2)**2)*normvar(0), smalldouble)
       case( 'slab' )
          {^IFTHREED
          G0(ixO^S) = max(dsqrt((x(ixO^S, 1)-x1ptms)**2 + (x(ixO^S, 2)-x2ptms)**2  &
               + (x(ixO^S, 3)-x3ptms)**2)*normvar(0), smalldouble)
          }
       end select

       G0(ixO^S) = 2.1d4*(Lstar/1.0d8)*((3.0857d17/G0(ixO^S))**2)

       select case( trim(dustspecies) )
       case( 'graphite' )
          Td(ixO^S, :) = 61.0d0*((0.0001d0/(sdust(:)*normvar(0)))**0.06d0) &
               *(G0(ixO^S)**(one/5.8d0))
       case( 'silicate' )
          Td(ixO^S, :) = 50.0d0*((0.0001d0/(sdust(:)*normvar(0)))**0.06d0) &
               *(G0(ixO^S)**(one/6.0d0))
       case default
          call mpistop( "=== Dust species undetermined===" )
       end select
    case default
       call mpistop( "=== Dust temperature undetermined===" )
    end select

  end subroutine get_tdust

  ! w[iw]= w[iw]+qdt*S[wCT, qtC, x] where S is the source based on wCT within ixO
  subroutine addsource(qdt, ixI^L, ixO^L, iw^LIM, qtC, wCT, qt, w, x, qsourcesplit)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L, iw^LIM
    double precision, intent(in)    :: qdt, qtC, qt, x(ixI^S, 1:ndim)
    double precision, intent(inout) :: wCT(ixI^S, 1:nw), w(ixI^S, 1:nw)
    logical, intent(in)             :: qsourcesplit

    double precision, dimension(ixG^T, 1:^NC, 1:n_dust) :: fdrag
    integer                                             :: n, idir

    select case( TRIM(dustmethod) )
    case( 'none' )
       !do nothing here
    case default !all regular dust methods here
       if (qsourcesplit .eqv. ssplitdust) then
          call get_3d_dragforce(ixI^L, ixO^L, wCT, x, fdrag)

          do idir = 1, ndir
             fdrag(ixO^S, idir, :) = fdrag(ixO^S, idir, :) * qdt

             w(ixO^S, mom(idir))  = w(ixO^S, mom(idir)) + &
                  sum(fdrag(ixO^S, idir, :), dim=ndim+2)

             if (hd_energy) then
                w(ixO^S, e_) = w(ixO^S, e_) + (wCT(ixO^S, mom(idir)) / &
                     wCT(ixO^S, rho_)) * sum(fdrag(ixO^S, idir, :), dim=ndim+2)
             end if
             w(ixO^S,m_dust(idir, :)) = w(ixO^S,m_dust(idir, :)) - fdrag(ixO^S, idir, :)
          end do

          if ( dustzero ) call set_dusttozero(qdt, ixI^L, ixO^L, iw^LIM, qtC, wCT, qt, w, x)
       endif
    end select

  end subroutine addsource

  subroutine dust_get_dt(w, ixI^L, ixO^L, dtnew, dx^D, x)
    use mod_global_parameters

    integer, intent(in)                        :: ixI^L, ixO^L
    double precision, intent(in)               :: dx^D, x(ixI^S, 1:ndim)
    double precision, intent(inout)            :: w(ixI^S, 1:nw), dtnew

    integer                                    :: idims, idust, idir
    double precision, dimension(1:n_dust)        :: dtdust
    double precision, dimension(ixG^T)         :: vt2, deltav, tstop, ptherm, vdust, vgas
    double precision, dimension(ixG^T, 1:n_dust) :: alpha_T
    double precision                           :: K

    !-----------------------------
    !get dt related to dust and gas stopping time (Laibe 2011)
    !-----------------------------

    select case( TRIM(dustmethod) )

    case( 'Kwok' ) ! assume sticking coefficient equals 0.25
       dtdust(:) = bigdouble

       call getpthermal(w, x, ixI^L, ixO^L, ptherm)
       vt2(ixO^S) = 3.0d0*ptherm(ixO^S)/w(ixO^S, rho_)

       ! Tgas, mu = mean molecular weight
       ptherm(ixO^S) = ( ptherm(ixO^S)*normvar(e_)*mhcgspar*eqpar(mu_))/(w(ixO^S, rho_)*normvar(rho_)*kbcgspar)

       do idir = 1, ndir
          call hd_get_v(w, x, ixI^L, ixO^L, idir, vgas)

          do n = 1, n_dust
             where(w(ixO^S, rho_dust(n))>minrhod)
                vdust(ixO^S)  = w(ixO^S,m_dust(idir, n))/w(ixO^S, rho_dust(n))
                deltav(ixO^S) = (vgas(ixO^S)-vdust(ixO^S))
                tstop(ixO^S)  = 4.0d0*(rhodust(n)*sdust(n))/ &
                     (3.0d0*(0.75d0)*dsqrt(vt2(ixO^S) + &
                     deltav(ixO^S)**2)*(w(ixO^S, rho_dust(n)) + &
                     w(ixO^S, rho_)))
             else where
                tstop(ixO^S) = bigdouble
             end where

             dtdust(n) = min(minval(tstop(ixO^S)), dtdust(n))
          end do
       end do

       dtnew = min(minval(dtdiffpar*dtdust(:)), dtnew)

    case( 'sticking' ) ! Calculate sticking coefficient based on the gas temperature
       dtdust(:) = bigdouble

       call getpthermal(w, x, ixI^L, ixO^L, ptherm)
       vt2(ixO^S) = 3.0d0*ptherm(ixO^S)/w(ixO^S, rho_)

       ! Sticking coefficient
       call get_sticking(w, x, ixI^L, ixO^L, alpha_T )

       ! Tgas, mu = mean molecular weight
       ptherm(ixO^S) = ( ptherm(ixO^S)*normvar(e_) * mhcgspar*eqpar(mu_)) / &
            (w(ixO^S, rho_)*normvar(rho_)*kbcgspar)

       do idir = 1, ndir
          call hd_get_v(w, x, ixI^L, ixO^L, idir, vgas)

          do n = 1, n_dust
             where(w(ixO^S, rho_dust(n))>minrhod)
                vdust(ixO^S)  = w(ixO^S,m_dust(idir, n))/w(ixO^S, rho_dust(n))
                deltav(ixO^S) = (vgas(ixO^S)-vdust(ixO^S))
                tstop(ixO^S)  = 4.0d0*(rhodust(n)*sdust(n))/ &
                     (3.0d0*(one-alpha_T(ixO^S,n))*dsqrt(vt2(ixO^S) + &
                     deltav(ixO^S)**2)*(w(ixO^S, rho_dust(n)) + &
                     w(ixO^S, rho_)))
             else where
                tstop(ixO^S) = bigdouble
             end where

             dtdust(n) = min(minval(tstop(ixO^S)), dtdust(n))
          end do
       end do

       dtnew = min(minval(dtdiffpar*dtdust(:)), dtnew)

    case('linear') !linear with Deltav, for testing (see Laibe & Price 2011)
       K = 3.4d5/n_dust
       dtdust(:) = bigdouble

       do n = 1, n_dust
          where(w(ixO^S, rho_dust(n))>minrhod)
             tstop(ixO^S)  = (w(ixO^S, rho_dust(n))*w(ixO^S, rho_))/ &
                  (K*(w(ixO^S, rho_dust(n)) + w(ixO^S, rho_)))
          else where
             tstop(ixO^S) = bigdouble
          end where

          dtdust(n) = min(minval(tstop(ixO^S)), dtdust(n))
       end do

       dtnew = min(minval(dtdiffpar*dtdust(:)), dtnew)
    case('none')
       ! no dust timestep
       dtnew = bigdouble
    case default
       call mpistop( "=== This dust method has not been implemented===" )
    end select

    if (dtnew < dtmin) then
       write(unitterm,*)"-------------------------------------"
       write(unitterm,*)"Warning: found DUST related time step too small! dtnew=", dtnew
       write(unitterm,*)"on grid with index:", saveigrid," grid level=", node(plevel_, saveigrid)
       write(unitterm,*)"grid corners are=",{^D&rnode(rpxmin^D_, saveigrid), rnode(rpxmax^D_, saveigrid)}
       write(unitterm,*)" dtdust =", dtdust(:)
       write(unitterm,*)"on processor:", mype
       write(unitterm,*)"-------------------------------------"
    endif

  end subroutine dust_get_dt

  subroutine dust_get_cmax(w, x, ixI^L, ixO^L, idim, cmax, cmin)
    use mod_global_parameters

    integer, intent(in)                       :: ixI^L, ixO^L, idim
    double precision, intent(in)              :: w(ixI^S, nw), x(ixI^S, 1:^ND)
    double precision, intent(inout)           :: cmax(ixG^T)
    double precision, intent(inout), optional :: cmin(ixG^T)
    double precision                          :: csound(ixG^T)
    double precision                          :: vdust(ixG^T, n_dust)
    integer :: n

    do n = 1, n_dust
       vdust(ixO^S, n) = get_vdust(w, x, ixI^L, ixO^L, idim, n)
    end do

    if (present(cmin)) then
       cmin(ixO^S) = min(0.0d0, minval(vdust(ixO^S, n), dim=ndim+1))
       cmax(ixO^S) = max(0.0d0, maxval(vdust(ixO^S, n), dim=ndim+1))
    else
       cmax(ixO^S) = maxval(abs(vdust(ixO^S, n)), dim=ndim+1)
    end if

  end subroutine dust_get_cmax

end module mod_dust
