module mod_dust

  integer :: dust_num_species
  double precision, allocatable :: sdust(:), dsdust(:), rhodust(:), mhcgspar, kbcgspar
  double precision              :: Lstar, Tdust


contains

  subroutine dust_activate()

    TODO

    ! Set starting index of dust species
    do n = 1, dust_num_species
       rho_dust(n) = ix + 1
       m_dust(:, n) = [ix+2,]
       ix = ix + 1 + ndir
    end do

    if (dust_num_species > 0) then
       dust_ = ix + 1
       ix    = ix + dust_num_species * (1 + ndir)
    else
       dust_ = -1
    end if

  end subroutine dust_activate

  subroutine dust_conserve
    !{^DS&{^C&w(ixO^S, m_dust(i_dim, i_dust)) = w(ixO^S, rho_dust(i_dust))*w(ixO^S, v^Cd^DS_);}\}
  end subroutine dust_conserve

  subroutine dust_primitive
    ! do i_dust = 1, dust_num_species
    !    where(w(ixO^S, rho_dust(i_dust))>minrhod)
    !       do idir = 1, ndir
    !          w(ixO^S, m_dust(idir, i_dust)) = w(ixO^S, m_dust(idir, i_dust)) /&
    !               w(ixO^S, rho_dust(i_dust))
    !       end do
    !    elsewhere
    !       w(ixO^S, m_dust(:, i_dust)) = 0.0d0
    !    end where
    ! end if
  end subroutine dust_primitive

  subroutine getflux(w,x,ixI^L,ixO^L,iw,idims,f,transport)
    ! TODO: move this to dust module, see above comment
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
  end subroutine getflux
  
  !  Force dust density to zero if rho_dust <= minrhod
  subroutine set_dusttozero(qdt, ixI^L, ixO^L, iw^LIM, qtC, wCT, qt, w, x)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L, iw^LIM
    double precision, intent(in)    :: qdt, qtC, qt, x(ixI^S, 1:ndim)
    double precision, intent(inout) :: wCT(ixI^S, 1:nw), w(ixI^S, 1:nw)

    TODO
    where( w(ixO^S, rho_dust(i_dust))<= minrhod )
       w(ixO^S, rho_dust(i_dust)) = 0.0d0
       w(ixO^S, m_dust(:, i_dust)) = 0.0d0
    end where

  end subroutine set_dusttozero

  ! Calculate drag force based on Epstein's or Stokes' law
  ! From Kwok 1975, page 584 (between eqn 8 and 9)
  subroutine get_3d_dragforce(ixI^L, ixO^L, w, x, fdrag)
    use mod_global_parameters

    integer, intent(in)                                            :: ixI^L, ixO^L
    double precision, intent(in)                                   :: x(ixI^S, 1:ndim)
    double precision, intent(inout)                                :: w(ixI^S, 1:nw)
    double precision, intent(out), dimension(ixG^T, 1:^NC, 1:^NDS) :: fdrag

    double precision, dimension(ixG^T)         :: vt2, deltav, fd, ptherm, vdust, vgas, Tgas
    double precision, dimension(ixG^T, 1:^NDS) :: alpha_T
    integer                                    :: idir
    double precision                           :: K

    call getpthermal(w, x, ixI^L, ixO^L, ptherm)

    vt2(ixO^S) = 3.0d0*ptherm(ixO^S)/w(ixO^S, rho_)

    select case( TRIM(dustmethod) )
    case( 'Kwok' ) ! assume sticking coefficient equals 0.25

       do idir = 1,^NC
          call getv(w, x, ixI^L, ixO^L, idir, vgas)

          todo
          where(w(ixO^S, rho_dust(i_dust))>minrhod)
             vdust(ixO^S)  = w(ixO^S,(rho_dust(i_dust))+idir*^NDS)/w(ixO^S, rho_dust(i_dust))
             deltav(ixO^S) = (vgas(ixO^S)-vdust(ixO^S))
             ! 0.75 from sticking coefficient
             fd(ixO^S)     = 0.75d0*w(ixO^S, rho_dust(i_dust))*w(ixO^S, rho_)*deltav(ixO^S) &
                  / (rhodust(^DS)*sdust(^DS))
             ! 0.75 from spherical grainvolume
             fd(ixO^S)     = -fd(ixO^S)*0.75d0*dsqrt(vt2(ixO^S) + deltav(ixO^S)**2)
          else where
             fd(ixO^S) = zero
          end where
          fdrag(ixO^S, idir,^DS) = fd(ixO^S)
       enddo

    case( 'sticking' ) ! Calculate sticking coefficient based on the gas and dust temperatures
       !
       !  Equation from Decin et al. 2006
       !
       TODO: only allowed if (hd_energy)
       Tgas(ixO^S) = ( ptherm(ixO^S)*normvar(p_)*mhcgspar)/(w(ixO^S, rho_)*normvar(rho_)*kbcgspar)
       call get_sticking(w, x, ixI^L, ixO^L, alpha_T)

       do idir = 1,^NC
          call getv(w, x, ixI^L, ixO^L, idir, vgas)

          TODO
          where(w(ixO^S, rho_dust(i_dust))>minrhod)
             vdust(ixO^S)  = w(ixO^S,(rho_dust(i_dust))+idir*^NDS)/w(ixO^S, rho_dust(i_dust))
             deltav(ixO^S) = (vgas(ixO^S)-vdust(ixO^S))
             fd(ixO^S)     = (one-alpha_T(ixO^S,^DS))*w(ixO^S, rho_dust(i_dust))*w(ixO^S, rho_)* &
                  deltav(ixO^S) / (rhodust(^DS)*sdust(^DS))
             fd(ixO^S)     = -fd(ixO^S)*0.75d0*dsqrt(vt2(ixO^S) + deltav(ixO^S)**2)
          else where
             fd(ixO^S) = zero
          end where
          fdrag(ixO^S, idir,^DS) = fd(ixO^S)
       enddo
    case('linear') !linear with Deltav, for testing (see Laibe & Price 2011)
       K = 3.4d5/^NDS
       do idir = 1,^NC
          call getv(w, x, ixI^L, ixO^L, idir, vgas)

          TODO
          where(w(ixO^S, rho_dust(i_dust))>minrhod)
             vdust(ixO^S)  = w(ixO^S,(rho_dust(i_dust))+idir*^NDS)/w(ixO^S, rho_dust(i_dust))
             deltav(ixO^S) = (vgas(ixO^S)-vdust(ixO^S))

             fd(ixO^S)     = -K*deltav(ixO^S)
          else where
             fd(ixO^S) = zero
          end where
          fdrag(ixO^S, idir,^DS) = fd(ixO^S)
       enddo
    case('none')
       {^DS&fdrag(ixO^S, idir,^DS) = zero\}
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


    double precision, intent(out)   :: alpha_T(ixG^T, 1:^NDS)

    double precision                :: Tgas(ixG^T)
    !-----------------------------------------------------------------------------

    call getpthermal(w, x, ixI^L, ixO^L, Tgas)
    call get_tdust(w, x, ixI^L, ixO^L, alpha_T)

    Tgas(ixO^S) = (Tgas(ixO^S)*normvar(p_)*mhcgspar)/(w(ixO^S, rho_)*normvar(rho_)*kbcgspar)

    {^DS&alpha_T(ixO^S,^DS) = max(0.35d0*dexp(-dsqrt((Tgas(ixO^S)+alpha_T(ixO^S,^DS))/5.0d2))+0.1d0, smalldouble)\}

    return
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
    double precision, intent(out)   :: Td(ixG^T, 1:^NDS)
    double precision                :: G0(ixO^S)

    select case( TRIM(dusttemp) )
    case( 'constant' )
       Td(ixO^S, 1:^NDS) = Tdust
    case( 'ism' )
       todo
       select case( trim(dustspecies) )
       case( 'graphite' )
          Td(ixO^S,^DS) = 15.8d0*((0.0001d0/(sdust(^DS)*normvar(0)))**0.06d0)
       case( 'silicate' )
          Td(ixO^S,^DS) = 13.6d0*((0.0001d0/(sdust(^DS)*normvar(0)))**0.06d0)
       case default
          call mpistop( "=== Dust species undetermined===" )
       end select
       \}
    case( 'stellar' )
       select case( TRIM(typeaxial) )
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
       todo
       select case( trim(dustspecies) )
       case( 'graphite' )
          Td(ixO^S,^DS) = 61.0d0*((0.0001d0/(sdust(^DS)*normvar(0)))**0.06d0) &
               *(G0(ixO^S)**(one/5.8d0))
       case( 'silicate' )
          Td(ixO^S,^DS) = 50.0d0*((0.0001d0/(sdust(^DS)*normvar(0)))**0.06d0) &
               *(G0(ixO^S)**(one/6.0d0))
       case default
          call mpistop( "=== Dust species undetermined===" )
       end select
       \}
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

    double precision, dimension(ixG^T, 1:^NC, 1:^NDS) :: fdrag
    integer                                           :: idir

    select case( TRIM(dustmethod) )
    case( 'none' )
       !do nothing here
    case default !all regular dust methods here
       if (qsourcesplit .eqv. ssplitdust) then
          call get_3d_dragforce(ixI^L, ixO^L, wCT, x, fdrag)
          do idir = 1, ndir
             fdrag(ixO^S, idir, 1:^NDS) = fdrag(ixO^S, idir, 1:^NDS)*qdt
             TODO
             w(ixO^S, m0_+idir)  = w(ixO^S, m0_+idir)  + fdrag(ixO^S, idir,^DS)
             if (hd_energy) then
                w(ixO^S, e_) = w(ixO^S, e_)        + (wCT(ixO^S, m0_+idir)/ &
                     wCT(ixO^S, rho_))*fdrag(ixO^S, idir,^DS)
             end if
             w(ixO^S,(rho_dust(i_dust))+idir*^NDS) = w(ixO^S,(rho_dust(i_dust))+idir*^NDS) - fdrag(ixO^S, idir,^DS)
          end do

          if ( dustzero ) call set_dusttozero(qdt, ixI^L, ixO^L, iw^LIM, qtC, wCT, qt, w, x)
       endif
    end select

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

  subroutine dust_get_cmax(new_cmax, w, x, ixI^L, ixO^L, idims, cmax, cmin, needcmin)
    use mod_global_parameters

    logical                      :: new_cmax, needcmin
    integer, intent(in)          :: ixI^L, ixO^L, idims
    double precision             :: w(ixI^S, nw), cmax(ixG^T), cmin(ixG^T)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision             :: csound(ixG^T){#IFDEF DUST , speeddust(ixG^T, 1:^NDS)}

    do i_dust = 1, hd_n_dust
       where(w(ixO^S, rho_dust(i_dust)) > minrhod)
          speeddust(ixO^S, i_dust) = w(ixO^S, m_dust(i_dims, i_dust)) / &
               w(ixO^S, rho_dust(i_dust));
       elsewhere
          speeddust(ixO^S, i_dust) = zero;
       end where
    end do

  end subroutine dust_get_cmax


end module mod_dust