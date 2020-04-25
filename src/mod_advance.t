!> Module containing all the time stepping schemes
module mod_advance

  implicit none
  private

  logical :: firstsweep, lastsweep

  !> Whether to conserve fluxes at the current sub-step
  logical :: fix_conserve_at_step = .true.

  public :: advance
  public :: process
  public :: process_advanced

contains

  !> Advance all the grids over one time step, including all sources
  subroutine advance(iit)
    use mod_global_parameters
    use mod_particles, only: handle_particles
    use mod_source, only: add_split_source

    integer, intent(in) :: iit

    integer :: iigrid, igrid, idimsplit

    ! old solution values at t_n-1 no longer needed: make copy of w(t_n)
    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       pso(igrid)%w(ixG^T,1:nwflux+nwaux)=ps(igrid)%w(ixG^T,1:nwflux+nwaux)
       if(stagger_grid) pso(igrid)%ws=ps(igrid)%ws
    end do
    !$OMP END PARALLEL DO

    {#IFDEF RAY
    call update_rays
    }

    ! split source addition
    call add_split_source(prior=.true.)

    firstsweep=.true.
    if (dimsplit) then
       if ((iit/2)*2==iit .or. typedimsplit=='xy') then
          ! do the sweeps in order of increasing idim,
          do idimsplit=1,ndim
             lastsweep= idimsplit==ndim
             call advect(idimsplit,idimsplit)
          end do
       else
          ! If the parity of "iit" is odd and typedimsplit=xyyx,
          ! do sweeps backwards
          do idimsplit=ndim,1,-1
             lastsweep= idimsplit==1
             call advect(idimsplit,idimsplit)
          end do
       end if
    else
       ! Add fluxes from all directions at once
       lastsweep= .true.
       call advect(1,ndim)
    end if

    ! split source addition
    call add_split_source(prior=.false.)

    if(use_particles) call handle_particles

  end subroutine advance

  !> Advance all grids over one time step, but without taking dimensional
  !> splitting or split source terms into account
  subroutine advect(idim^LIM)
    use mod_global_parameters
    use mod_fix_conserve

    integer, intent(in) :: idim^LIM
    integer             :: iigrid, igrid

    call init_comm_fix_conserve(idim^LIM,nwflux)
    fix_conserve_at_step = time_advance .and. levmax>levmin

    ! copy w instead of wold because of potential use of dimsplit or sourcesplit
    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       ps1(igrid)%w=ps(igrid)%w
       if(stagger_grid) ps1(igrid)%ws=ps(igrid)%ws
    end do
    !$OMP END PARALLEL DO

    istep = 0

    select case (time_integrator)
    case ("onestep")
       call advect1(flux_scheme,one,idim^LIM,global_time,ps1,global_time,ps)

    case ("twostep")
      ! predictor step
       fix_conserve_at_step = .false.
       call advect1(typepred1,half, idim^LIM,global_time,ps,global_time,ps1)

       ! corrector step
       fix_conserve_at_step = time_advance .and. levmax>levmin
       call advect1(flux_scheme,one,idim^LIM,global_time+half*dt,ps1,global_time,ps)

    case ("twostep_trapezoidal")
       ! Explicit trapezoidal rule / Euler's method / Heun's method
       ! In the future, this could be implemented using just w and w1
       call advect1(flux_scheme, 1.0d0, idim^LIM, global_time, ps, global_time, ps1)

       do iigrid=1,igridstail; igrid=igrids(iigrid);
          ps2(igrid)%w = ps(igrid)%w
          ps(igrid)%w = ps1(igrid)%w
          if(stagger_grid) then
            ps2(igrid)%ws = ps(igrid)%ws
            ps(igrid)%ws = ps1(igrid)%ws
          end if
       end do

       call advect1(flux_scheme, 1.0d0, idim^LIM, global_time+dt, ps1, global_time+dt, ps)

       do iigrid=1,igridstail; igrid=igrids(iigrid);
          ps(igrid)%w = 0.5d0 * (ps(igrid)%w + ps2(igrid)%w)
       end do
    case ("threestep")
       ! three step Runge-Kutta in accordance with Gottlieb & Shu 1998
       call advect1(flux_scheme,one, idim^LIM,global_time,ps,global_time,ps1)
       do iigrid=1,igridstail; igrid=igrids(iigrid);
          ps2(igrid)%w(ixG^T,1:nwflux)=0.75d0*ps(igrid)%w(ixG^T,1:nwflux)+0.25d0*&
               ps1(igrid)%w(ixG^T,1:nwflux)
          if (nw>nwflux) ps2(igrid)%w(ixG^T,nwflux+1:nw) = &
               ps(igrid)%w(ixG^T,nwflux+1:nw)
          if(stagger_grid) ps2(igrid)%ws=0.75d0*ps(igrid)%ws+0.25d0*ps1(igrid)%ws
       end do

       call advect1(flux_scheme,0.25d0, idim^LIM,global_time+dt,ps1,global_time+dt*0.25d0,ps2)

       !$OMP PARALLEL DO PRIVATE(igrid)
       do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
          ps(igrid)%w(ixG^T,1:nwflux)=1.0d0/3.0d0*ps(igrid)%w(ixG^T,1:nwflux)+&
               2.0d0/3.0d0*ps2(igrid)%w(ixG^T,1:nwflux)
          if(stagger_grid) ps(igrid)%ws=1.0d0/3.0d0*ps(igrid)%ws+2.0d0/3.0d0*ps2(igrid)%ws
       end do
       !$OMP END PARALLEL DO
       call advect1(flux_scheme,2.0d0/3.0d0, idim^LIM,global_time+dt/2.0d0,ps2,global_time+dt/3.0d0,ps)

    case ("ssprk43")
       ! Strong stability preserving 4 stage RK 3rd order method by Ruuth and Spiteri
       !
       ! Ruuth & Spiteri
       ! J. S C, 17 (2002) p. 211 - 220
       !
       ! supposed to be stable up to CFL=2.
       ! don't use time-dependent sources since i did not bother to set those intermediate times.
       ! oliver.

       ! === First step ===
       call advect1(flux_scheme,0.5d0, idim^LIM,global_time,ps,global_time,ps1)

       ! === Second step ===
       do iigrid=1,igridstail; igrid=igrids(iigrid);
          ps2(igrid)%w(ixG^T,1:nwflux)=ps1(igrid)%w(ixG^T,1:nwflux)
          if (nw>nwflux) ps2(igrid)%w(ixG^T,nwflux+1:nw) = &
               ps(igrid)%w(ixG^T,nwflux+1:nw)
          if(stagger_grid) ps2(igrid)%ws=ps1(igrid)%ws
       end do
       call advect1(flux_scheme,0.5d0, idim^LIM,global_time,ps1,global_time+dt,ps2)

       ! === Third step ===
       do iigrid=1,igridstail; igrid=igrids(iigrid);
          ps3(igrid)%w(ixG^T,1:nwflux)=2.0d0/3.0d0 * ps(igrid)%w(ixG^T,1:nwflux) &
               + 1.0d0/3.0d0 * ps2(igrid)%w(ixG^T,1:nwflux)
          if (nw>nwflux) ps3(igrid)%w(ixG^T,nwflux+1:nw) = &
               ps(igrid)%w(ixG^T,nwflux+1:nw)
          if(stagger_grid) ps3(igrid)%ws=2.0d0/3.0d0*ps(igrid)%ws+1.0d0/3.0d0*ps2(igrid)%ws
       end do
       call advect1(flux_scheme,1.0d0/6.0d0, idim^LIM,global_time,ps2,global_time+dt,ps3)

       ! === Fourth step ===
       !$OMP PARALLEL DO PRIVATE(igrid)
       do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
          ps(igrid)%w(ixG^T,1:nwflux)=ps3(igrid)%w(ixG^T,1:nwflux)
          if(stagger_grid) ps3(igrid)%ws=ps(igrid)%ws
       end do
       !$OMP END PARALLEL DO
       call advect1(flux_scheme,0.5d0, idim^LIM, global_time,ps3, global_time+dt,ps)

    case ("ssprk54")
       ! Strong stability preserving 5 stage RK 4th order method by Ruuth and Spiteri
       !
       ! SIAM J. NUMER. ANAL.
       ! Vol. 40, No. 2, pp. 469–491
       ! c 2002 Society for Industrial and Applied Mathematics
       ! A NEW CLASS OF OPTIMAL
       ! HIGH-ORDER STRONG-STABILITY-PRESERVING
       ! TIME DISCRETIZATION METHODS
       !
       ! E.g. Table A.2
       !
       ! I have the actual coefficients however from the overview article by Gottlieb, JoSC 25 (2005)
       ! ON HIGH ORDER STRONG STABILITY PRESERVING RUNGE-KUTTA AND MULTI STEP TIME DISCRETIZATIONS
       !
       ! there are slight differences in the coefficients (~8th digit after the .)
       ! This is SSP till CFL number 1.508 which makes the effective CFL per step ceff=0.377,
       ! in contrast to the ceff=0.3333 of the classical RK3.
       !
       ! coded by oliver on 11/05/2013.  Enjoy!

       ! === First step ===
       call advect1(flux_scheme,0.391752226571890d0, idim^LIM, global_time,ps,global_time,ps1)

       ! === Second step ===
       do iigrid=1,igridstail; igrid=igrids(iigrid);
          ps2(igrid)%w(ixG^T,1:nwflux)=0.444370493651235d0 * ps(igrid)%w(ixG^T,1:nwflux) &
               + 0.555629506348765d0 * ps1(igrid)%w(ixG^T,1:nwflux)
          if (nw>nwflux) ps2(igrid)%w(ixG^T,nwflux+1:nw) = ps(igrid)%w(ixG^T,nwflux+1:nw)
          if(stagger_grid) ps2(igrid)%ws=0.444370493651235d0*ps(igrid)%ws+0.555629506348765d0*ps1(igrid)%ws
       end do
       call advect1(flux_scheme,0.368410593050371d0, idim^LIM, global_time+0.391752226571890d0*dt,ps1, global_time+0.2176690962611688d0*dt,ps2)

       ! === Third step ===
       do iigrid=1,igridstail; igrid=igrids(iigrid);
          ps3(igrid)%w(ixG^T,1:nwflux)=0.620101851488403d0 * ps(igrid)%w(ixG^T,1:nwflux) &
               + 0.379898148511597d0 * ps2(igrid)%w(ixG^T,1:nwflux)
          if (nw>nwflux) ps3(igrid)%w(ixG^T,nwflux+1:nw) = ps(igrid)%w(ixG^T,nwflux+1:nw)
          if(stagger_grid) ps3(igrid)%ws=0.620101851488403d0 * ps(igrid)%ws &
               + 0.379898148511597d0 * ps2(igrid)%ws
       end do
       call advect1(flux_scheme,0.251891774271694d0, idim^LIM, global_time+0.5860796893115398d0*dt,ps2, global_time+0.222650588849706d0*dt,ps3)

       ! === Fourth step ===
       do iigrid=1,igridstail; igrid=igrids(iigrid);
          ps4(igrid)%w(ixG^T,1:nwflux)=0.178079954393132d0 * ps(igrid)%w(ixG^T,1:nwflux) &
               + 0.821920045606868d0 * ps3(igrid)%w(ixG^T,1:nwflux)
          if (nw>nwflux) ps4(igrid)%w(ixG^T,nwflux+1:nw) = ps(igrid)%w(ixG^T,nwflux+1:nw)
          if(stagger_grid) ps4(igrid)%ws=0.178079954393132d0 * ps(igrid)%ws &
               + 0.821920045606868d0 * ps3(igrid)%ws
       end do
       call advect1(flux_scheme,0.544974750228521d0, idim^LIM, global_time+0.4745423631214d0*dt,ps3, global_time+0.390035880739132d0*dt,ps4)
       ! Now recover back the dt*L(u3), store in w1:
       !$OMP PARALLEL DO PRIVATE(igrid)
       do iigrid=1,igridstail; igrid=igrids(iigrid);
          ps1(igrid)%w(ixG^T,1:nwflux) = ( ps4(igrid)%w(ixG^T,1:nwflux) &
               - (0.178079954393132d0 * ps(igrid)%w(ixG^T,1:nwflux) &
               + 0.821920045606868d0 * ps3(igrid)%w(ixG^T,1:nwflux)) ) / 0.544974750228521d0
          if(stagger_grid) ps1(igrid)%ws = ( ps4(igrid)%ws &
               - (0.178079954393132d0 * ps(igrid)%ws &
               + 0.821920045606868d0 * ps3(igrid)%ws) ) / 0.544974750228521d0
       end do
       !$OMP END PARALLEL DO

       ! === Fifth step ===
       !$OMP PARALLEL DO PRIVATE(igrid)
       do iigrid=1,igridstail; igrid=igrids(iigrid);
          ps(igrid)%w(ixG^T,1:nwflux)= 0.517231671970585d0 * ps2(igrid)%w(ixG^T,1:nwflux) &
               + 0.096059710526147d0 * ps3(igrid)%w(ixG^T,1:nwflux) &
               + 0.063692468666290d0 * ps1(igrid)%w(ixG^T,1:nwflux) &
               + 0.386708617503269d0 * ps4(igrid)%w(ixG^T,1:nwflux)
          if(stagger_grid) ps(igrid)%ws= 0.517231671970585d0 * ps2(igrid)%ws &
               + 0.096059710526147d0 * ps3(igrid)%ws &
               + 0.063692468666290d0 * ps1(igrid)%ws &
               + 0.386708617503269d0 * ps4(igrid)%ws
       end do
       !$OMP END PARALLEL DO
       call advect1(flux_scheme,0.226007483236906d0, idim^LIM, global_time+0.935010630967653d0*dt,ps4, global_time+0.710300048096804d0*dt,ps)


    case ("rk4")
       ! classical RK4 four step scheme
       call advect1(flux_scheme,0.5d0, idim^LIM,global_time,ps,global_time,ps1)
       do iigrid=1,igridstail; igrid=igrids(iigrid);
          ps2(igrid)%w(ixG^T,1:nwflux)=ps(igrid)%w(ixG^T,1:nwflux)
          if(stagger_grid) ps2(igrid)%ws=ps(igrid)%ws
       end do
       call advect1(flux_scheme,0.5d0, idim^LIM,global_time+dt/2d0,ps1,global_time,ps2)
       do iigrid=1,igridstail; igrid=igrids(iigrid);
          ps3(igrid)%w(ixG^T,1:nwflux)=ps(igrid)%w(ixG^T,1:nwflux)
          if(stagger_grid) ps3(igrid)%ws=ps(igrid)%ws
       end do
       call advect1(flux_scheme,one,   idim^LIM,global_time+dt/2d0,ps2,global_time,ps3)
       !$OMP PARALLEL DO PRIVATE(igrid)
       do iigrid=1,igridstail; igrid=igrids(iigrid);
          ps1(igrid)%w(ixG^T,1:nwflux)=(ps1(igrid)%w(ixG^T,1:nwflux) &
               +two*ps2(igrid)%w(ixG^T,1:nwflux) &
               +ps3(igrid)%w(ixG^T,1:nwflux) &
               -4.0d0*ps(igrid)%w(ixG^T,1:nwflux))/3.0d0
          if(stagger_grid) ps1(igrid)%ws=(ps1(igrid)%ws &
               +two*ps2(igrid)%ws &
               +ps3(igrid)%ws &
               -4.0d0*ps(igrid)%ws)/3.0d0
       end do
       !$OMP END PARALLEL DO
       call advect1(flux_scheme,1.0d0/6.0d0, idim^LIM,global_time+dt,ps3,global_time,ps)
       !$OMP PARALLEL DO PRIVATE(igrid)
       do iigrid=1,igridstail; igrid=igrids(iigrid);
          ps(igrid)%w(ixG^T,1:nwflux)=ps1(igrid)%w(ixG^T,1:nwflux)+&
               ps(igrid)%w(ixG^T,1:nwflux)
          if(stagger_grid) ps(igrid)%ws=ps1(igrid)%ws+ps(igrid)%ws
       end do
       !$OMP END PARALLEL DO
    case ("fourstep")
       ! four step scheme, variant Hans De Sterck
       call advect1(flux_scheme,0.12d0, idim^LIM,global_time,ps,global_time,ps1)
       do iigrid=1,igridstail; igrid=igrids(iigrid);
          ps2(igrid)%w(ixG^T,1:nwflux)=ps(igrid)%w(ixG^T,1:nwflux)
          if(stagger_grid) ps2(igrid)%ws=ps(igrid)%ws
       end do
       call advect1(flux_scheme,0.25d0, idim^LIM,global_time+dt*0.12d0,ps1,global_time,ps2)
       !$OMP PARALLEL DO PRIVATE(igrid)
       do iigrid=1,igridstail; igrid=igrids(iigrid);
          ps1(igrid)%w(ixG^T,1:nwflux)=ps(igrid)%w(ixG^T,1:nwflux)
          if(stagger_grid) ps1(igrid)%ws=ps(igrid)%ws
       end do
       !$OMP END PARALLEL DO
       call advect1(flux_scheme,0.5d0,  idim^LIM,global_time+dt/4d0   ,ps2,global_time,ps1)
       call advect1(flux_scheme,one,    idim^LIM,global_time+dt/2d0   ,ps1,global_time,ps)
    case ("jameson")
       ! four step scheme, variant jameson
       call advect1(flux_scheme,0.25d0, idim^LIM,global_time,ps,global_time,ps1)
       do iigrid=1,igridstail; igrid=igrids(iigrid);
          ps2(igrid)%w(ixG^T,1:nwflux)=ps(igrid)%w(ixG^T,1:nwflux)
          if(stagger_grid) ps2(igrid)%ws=ps(igrid)%ws
       end do
       call advect1(flux_scheme,(1.0d0/3.0d0), idim^LIM,global_time+dt*0.25d0,ps1,global_time,ps2)
       !$OMP PARALLEL DO PRIVATE(igrid)
       do iigrid=1,igridstail; igrid=igrids(iigrid);
          ps1(igrid)%w(ixG^T,1:nwflux)=ps(igrid)%w(ixG^T,1:nwflux)
          if(stagger_grid) ps1(igrid)%ws=ps(igrid)%ws
       end do
       !$OMP END PARALLEL DO
       call advect1(flux_scheme,0.5d0,  idim^LIM,global_time+dt/3d0,ps2,global_time,ps1)
       call advect1(flux_scheme,one,    idim^LIM,global_time+dt/2d0,ps1,global_time,ps)
    case default
       write(unitterm,*) "time_integrator=",time_integrator
       write(unitterm,*) "Error in advect: Unknown time integration method"
       call mpistop("Correct time_integrator")
    end select

    firstsweep=.false.
  end subroutine advect

  !> Integrate all grids by one partial step
  subroutine advect1(method,dtfactor,idim^LIM,qtC,psa,qt,psb)
    use mod_global_parameters
    use mod_ghostcells_update
    use mod_fix_conserve
    use mod_physics

    integer, intent(in) :: idim^LIM
    type(state), target :: psa(max_blocks) !< Compute fluxes based on this state
    type(state), target :: psb(max_blocks) !< Update solution on this state
    double precision, intent(in) :: dtfactor !< Advance over dtfactor * dt
    double precision, intent(in) :: qtC
    double precision, intent(in) :: qt
    character(len=*), intent(in) :: method(nlevelshi)

    double precision :: qdt
    integer :: iigrid, igrid, level

    logical :: setigrid

    istep = istep+1

    ! opedit: Just advance the active grids:
    !$OMP PARALLEL DO PRIVATE(igrid,level,qdt)
    do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
       level=node(plevel_,igrid)
       qdt=dtfactor*dt_grid(igrid)

       call process1_grid(method(level),igrid,qdt,ixG^LL,idim^LIM,qtC,&
            psa(igrid),qt,psb(igrid),pso(igrid))
    end do
    !$OMP END PARALLEL DO

    ! opedit: Send flux for all grids, expects sends for all
    ! nsend_fc(^D), set in connectivity.t.

    if (fix_conserve_global .and. fix_conserve_at_step) then
      call recvflux(idim^LIM)
      call sendflux(idim^LIM)
      call fix_conserve(psb,idim^LIM,1,nwflux)
      if(stagger_grid) call fix_edges(psb,idim^LIM)
    end if

    if(stagger_grid) then
      ! Now fill the cell-center values for the staggered variables
      !$OMP PARALLEL DO PRIVATE(igrid)
      do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
        call phys_face_to_center(ixG^LL,psb(igrid))
      end do
      !$OMP END PARALLEL DO
    end if

    ! For all grids: fill ghost cells
    qdt = dtfactor*dt
    call getbc(qt+qdt,qdt,psb,1,nwflux+nwaux,phys_req_diagonal)

  end subroutine advect1

  !> Prepare to advance a single grid over one partial time step
  subroutine process1_grid(method,igrid,qdt,ixI^L,idim^LIM,qtC,sCT,qt,s,sold)
    use mod_global_parameters
    use mod_fix_conserve

    character(len=*), intent(in) :: method
    integer, intent(in) :: igrid, ixI^L, idim^LIM
    double precision, intent(in) :: qdt, qtC, qt
    type(state), target          :: sCT, s, sold

    double precision :: dx^D
    ! cell face flux
    double precision :: fC(ixI^S,1:nwflux,1:ndim)
    ! cell edge flux
    double precision :: fE(ixI^S,7-2*ndim:3)

    dx^D=rnode(rpdx^D_,igrid);
    ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
    saveigrid=igrid

    block0=>sCT
    block=>s
    typelimiter=type_limiter(node(plevel_,igrid))
    typegradlimiter=type_gradient_limiter(node(plevel_,igrid))

    call advect1_grid(method,qdt,ixI^L,idim^LIM,qtC,sCT,qt,s,sold,fC,fE,dx^D, &
         ps(igrid)%x)


    ! opedit: Obviously, flux is stored only for active grids.
    ! but we know in fix_conserve wether there is a passive neighbor
    ! but we know in conserve_fix wether there is a passive neighbor
    ! via neighbor_active(i^D,igrid) thus we skip the correction for those.
    ! This violates strict conservation when the active/passive interface
    ! coincides with a coarse/fine interface.
    if (fix_conserve_global .and. fix_conserve_at_step) then
      call store_flux(igrid,fC,idim^LIM,nwflux)
      if(stagger_grid) call store_edge(igrid,ixI^L,fE,idim^LIM)
    end if

  end subroutine process1_grid

  !> Advance a single grid over one partial time step
  subroutine advect1_grid(method,qdt,ixI^L,idim^LIM,qtC,sCT,qt,s,sold,fC,fE,dx^D,x)

    !  integrate one grid by one partial step
    use mod_finite_volume
    use mod_finite_difference
    use mod_tvd
    use mod_source, only: addsource2
    use mod_global_parameters

    character(len=*), intent(in) :: method
    integer, intent(in) :: ixI^L, idim^LIM
    double precision, intent(in) :: qdt, qtC, qt, dx^D, x(ixI^S,1:ndim)
    type(state), target          :: sCT, s, sold
    double precision :: fC(ixI^S,1:nwflux,1:ndim)
    double precision :: fE(ixI^S,7-2*ndim:3)

    integer :: ixO^L

    ixO^L=ixI^L^LSUBnghostcells;
    select case (method)
    case ('tvdmu','tvdlf','hll','hllc','hllcd','hlld')
       call finite_volume(method,qdt,ixI^L,ixO^L,idim^LIM,qtC,sCT,qt,s,sold,fC,fE,dx^D,x)
    case ('cd','cd4')
       call centdiff(method,qdt,ixI^L,ixO^L,idim^LIM,qtC,sCT,qt,s,fC,fE,dx^D,x)
    case ('hancock')
       call hancock(qdt,ixI^L,ixO^L,idim^LIM,qtC,sCT,qt,s,dx^D,x)
    case ('fd')
       call fd(method,qdt,ixI^L,ixO^L,idim^LIM,qtC,sCT,qt,s,fC,fE,dx^D,x)
    case ('tvd')
       call centdiff('cd',qdt,ixI^L,ixO^L,idim^LIM,qtC,sCT,qt,s,fC,fE,dx^D,x)
       call tvdlimit(method,qdt,ixI^L,ixO^L,idim^LIM,sCT,qt+qdt,s,fC,dx^D,x)
    case ('source')
       call addsource2(qdt*dble(idimmax-idimmin+1)/dble(ndim),&
            ixI^L,ixO^L,1,nw,qtC,sCT%w,qt,s%w,x,.false.)
    case ('nul')
       ! There is nothing to do
    case default
       write(unitterm,*)'Error in advect1_grid:',method,' is unknown!'
       call mpistop("")
    end select

  end subroutine advect1_grid

  !> process is a user entry in time loop, before output and advance
  !>         allows to modify solution, add extra variables, etc.
  !> Warning: CFL dt already determined (and is not recomputed)! 
  subroutine process(iit,qt)
    use mod_usr_methods, only: usr_process_grid, usr_process_global
    use mod_global_parameters
    use mod_ghostcells_update
    ! .. scalars ..
    integer,intent(in)          :: iit
    double precision, intent(in):: qt

    integer:: iigrid, igrid,level

    if (associated(usr_process_global)) then
       call usr_process_global(iit,qt)
    end if

    if (associated(usr_process_grid)) then
      !$OMP PARALLEL DO PRIVATE(igrid,level)
      do iigrid=1,igridstail; igrid=igrids(iigrid);
         level=node(plevel_,igrid)
         ! next few lines ensure correct usage of routines like divvector etc
         ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
         block=>ps(igrid)
         typelimiter=type_limiter(node(plevel_,igrid))
         typegradlimiter=type_gradient_limiter(node(plevel_,igrid))
         call usr_process_grid(igrid,level,ixG^LL,ixM^LL,qt,ps(igrid)%w,ps(igrid)%x)
      end do
      !$OMP END PARALLEL DO
      call getbc(qt,dt,ps,1,nwflux+nwaux)
    end if
  end subroutine process

  !> process_advanced is user entry in time loop, just after advance
  !>           allows to modify solution, add extra variables, etc.
  !>           added for handling two-way coupled PIC-MHD
  !> Warning: w is now at global_time^(n+1), global time and iteration at global_time^n, it^n
  subroutine process_advanced(iit,qt)
    use mod_usr_methods, only: usr_process_adv_grid, &
                               usr_process_adv_global
    use mod_global_parameters
    use mod_ghostcells_update
    ! .. scalars ..
    integer,intent(in)          :: iit
    double precision, intent(in):: qt

    integer:: iigrid, igrid,level

    if (associated(usr_process_adv_global)) then
       call usr_process_adv_global(iit,qt)
    end if

    if (associated(usr_process_adv_grid)) then
      !$OMP PARALLEL DO PRIVATE(igrid,level)
      do iigrid=1,igridstail; igrid=igrids(iigrid);
         level=node(plevel_,igrid)
         ! next few lines ensure correct usage of routines like divvector etc
         ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
         block=>ps(igrid)
         typelimiter=type_limiter(node(plevel_,igrid))
         typegradlimiter=type_gradient_limiter(node(plevel_,igrid))

         call usr_process_adv_grid(igrid,level,ixG^LL,ixM^LL, &
              qt,ps(igrid)%w,ps(igrid)%x)
      end do
      !$OMP END PARALLEL DO
    call getbc(qt,dt,ps,1,nwflux+nwaux)
    end if
  end subroutine process_advanced

end module mod_advance
