!> Module containing all the time stepping schemes
module mod_advance

  implicit none
  private

  !> Whether to conserve fluxes at the current sub-step
  logical :: fix_conserve_at_step = .true.

  !> Per-rank compute-only accumulator for lb_diagnose (sums the iigrid
  !> block-loop wall time across all substages within one full advance call,
  !> excluding ghostcell exchanges and flux conservation collectives).
  double precision :: lb_compute_accum = 0.0d0

  public :: advance
  public :: process
  public :: process_advanced

contains

  !> Advance all the grids over one time step, including all sources
  subroutine advance(iit)
    use mod_global_parameters
    use mod_particles, only: handle_particles
    use mod_eos
    use mod_source, only: add_split_source
    use mod_supertimestepping, only: lb_tc_accum
    use mod_radiative_cooling, only: lb_cool_accum

    integer, intent(in) :: iit

    integer :: iigrid, igrid, idimsplit

    ! ---- per-rank load-balance diagnostic (gated by lb_diagnose) ----
    double precision :: t_advance_start, t_advance_local
    double precision :: t_dummy(1)
    double precision, allocatable, save :: t_all_ranks(:)
    double precision, allocatable, save :: c_all_ranks(:)
    double precision, allocatable, save :: tc_all_ranks(:)
    double precision, allocatable, save :: cool_all_ranks(:)
    double precision, allocatable, save :: rt_all_ranks(:)
    double precision :: cmax, cmean, ratioc
    double precision :: tcmax, tcmean, ratiotc
    double precision :: coolmax, coolmean, ratiocool
    double precision :: rtmax, rtmean, ratiort
    integer, save :: lb_log_unit = -1
    logical, save :: lb_first_call = .true.
    integer :: ipe
    double precision :: tmax, tmean, ratio
    character(len=256) :: lb_log_name
    ! ----------------------------------------------------------------

    if (lb_diagnose) then
      if (lb_first_call .and. mype==0) then
        allocate(t_all_ranks(npe))
        allocate(c_all_ranks(npe))
        allocate(tc_all_ranks(npe))
        allocate(cool_all_ranks(npe))
        allocate(rt_all_ranks(npe))
        write(lb_log_name,'(a,a)') trim(base_filename), 'rank_timing.log'
        open(newunit=lb_log_unit, file=trim(lb_log_name), status='replace', action='write')
        write(lb_log_unit,'(a)',advance='no') '# it time '
        do ipe=0,npe-1
          write(lb_log_unit,'(a,i0,a)',advance='no') 't_rank',ipe,' '
        end do
        do ipe=0,npe-1
          write(lb_log_unit,'(a,i0,a)',advance='no') 'c_rank',ipe,' '
        end do
        do ipe=0,npe-1
          write(lb_log_unit,'(a,i0,a)',advance='no') 't_tc_rank',ipe,' '
        end do
        do ipe=0,npe-1
          write(lb_log_unit,'(a,i0,a)',advance='no') 't_cool_rank',ipe,' '
        end do
        do ipe=0,npe-1
          write(lb_log_unit,'(a,i0,a)',advance='no') 't_rt_rank',ipe,' '
        end do
        write(lb_log_unit,'(a)',advance='no') 'tmax tmean R cmax cmean Rc '
        write(lb_log_unit,'(a)',advance='no') 'tcmax tcmean Rtc coolmax coolmean Rcool '
        write(lb_log_unit,'(a)') 'rtmax rtmean Rrt'
        flush(lb_log_unit)
        lb_first_call = .false.
      else if (lb_first_call) then
        lb_first_call = .false.
      end if
      t_advance_start = MPI_WTIME()
      lb_compute_accum = 0.0d0
      lb_tc_accum = 0.0d0
      lb_cool_accum = 0.0d0
    end if

    ! Per-block cost reset for the cost-weighted load balancer.
    ! Cleared every step before the iigrid loops fill it.
    ! Seed with the sweep cost measured by rt_sc_solve() earlier this step, then
    ! clear that accumulator for the next solve. The hydro timers below add on
    ! top, so the partitioner sees transfer and hydro cost for the same step.
    if (lb_automatic) then
      block_cost = block_cost_rt
      block_cost_rt = 0.0d0
    end if

    ! split source addition
    call add_split_source(prior=.true.) !> calculates temperature based on conservative state

    if(dimsplit) then
       if((iit/2)*2==iit .or. typedimsplit=='xy') then
          ! do the sweeps in order of increasing idim,
          do idimsplit=1,ndim
             call advect(idimsplit,idimsplit)
          end do
       else
          ! If the parity of "iit" is odd and typedimsplit=xyyx,
          ! do sweeps backwards
          do idimsplit=ndim,1,-1
             call advect(idimsplit,idimsplit)
          end do
       end if
    else
       ! Add fluxes from all directions at once
       call advect(1,ndim)
    end if

    ! split source addition
    call add_split_source(prior=.false.) !> calculates temperature based on conservative state

    if(use_particles) call handle_particles

    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
      call eos%update_eos(ixG^LL,ixG^LL,ps(igrid)%w,ps(igrid)%x)
    end do
    !$OMP END PARALLEL DO

    ! Cost-weighted load balancer: per-block measurements in block_cost
    ! (per-rank, per-igrid) are folded into the global Morton-indexed
    ! costlist via EWMA blend inside get_Morton_range_costed when
    ! load_balance is next invoked. No end-of-advance work needed here.

    if (lb_diagnose) then
      t_advance_local = MPI_WTIME() - t_advance_start
      if (mype==0) then
        call MPI_GATHER(t_advance_local,1,MPI_DOUBLE_PRECISION, &
                        t_all_ranks,1,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
        call MPI_GATHER(lb_compute_accum,1,MPI_DOUBLE_PRECISION, &
                        c_all_ranks,1,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
        call MPI_GATHER(lb_tc_accum,1,MPI_DOUBLE_PRECISION, &
                        tc_all_ranks,1,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
        call MPI_GATHER(lb_cool_accum,1,MPI_DOUBLE_PRECISION, &
                        cool_all_ranks,1,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
        call MPI_GATHER(lb_rt_accum,1,MPI_DOUBLE_PRECISION, &
                        rt_all_ranks,1,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
        tmax = maxval(t_all_ranks)
        tmean = sum(t_all_ranks)/dble(npe)
        cmax = maxval(c_all_ranks)
        cmean = sum(c_all_ranks)/dble(npe)
        tcmax = maxval(tc_all_ranks)
        tcmean = sum(tc_all_ranks)/dble(npe)
        coolmax = maxval(cool_all_ranks)
        coolmean = sum(cool_all_ranks)/dble(npe)
        rtmax = maxval(rt_all_ranks)
        rtmean = sum(rt_all_ranks)/dble(npe)
        if (tmean > 0.0d0) then
          ratio = tmax/tmean
        else
          ratio = 1.0d0
        end if
        if (cmean > 0.0d0) then
          ratioc = cmax/cmean
        else
          ratioc = 1.0d0
        end if
        if (tcmean > 0.0d0) then
          ratiotc = tcmax/tcmean
        else
          ratiotc = 1.0d0
        end if
        if (coolmean > 0.0d0) then
          ratiocool = coolmax/coolmean
        else
          ratiocool = 1.0d0
        end if
        if (rtmean > 0.0d0) then
          ratiort = rtmax/rtmean
        else
          ratiort = 1.0d0
        end if
        write(lb_log_unit,'(i10,1x,es16.8,1x)',advance='no') it, global_time
        do ipe=1,npe
          write(lb_log_unit,'(es14.6,1x)',advance='no') t_all_ranks(ipe)
        end do
        do ipe=1,npe
          write(lb_log_unit,'(es14.6,1x)',advance='no') c_all_ranks(ipe)
        end do
        do ipe=1,npe
          write(lb_log_unit,'(es14.6,1x)',advance='no') tc_all_ranks(ipe)
        end do
        do ipe=1,npe
          write(lb_log_unit,'(es14.6,1x)',advance='no') cool_all_ranks(ipe)
        end do
        do ipe=1,npe
          write(lb_log_unit,'(es14.6,1x)',advance='no') rt_all_ranks(ipe)
        end do
        write(lb_log_unit,'(15(es14.6,1x))') &
             tmax, tmean, ratio, cmax, cmean, ratioc, &
             tcmax, tcmean, ratiotc, coolmax, coolmean, ratiocool, &
             rtmax, rtmean, ratiort
        flush(lb_log_unit)
      else
        call MPI_GATHER(t_advance_local,1,MPI_DOUBLE_PRECISION, &
                        t_dummy,1,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
        call MPI_GATHER(lb_compute_accum,1,MPI_DOUBLE_PRECISION, &
                        t_dummy,1,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
        call MPI_GATHER(lb_tc_accum,1,MPI_DOUBLE_PRECISION, &
                        t_dummy,1,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
        call MPI_GATHER(lb_cool_accum,1,MPI_DOUBLE_PRECISION, &
                        t_dummy,1,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
        call MPI_GATHER(lb_rt_accum,1,MPI_DOUBLE_PRECISION, &
                        t_dummy,1,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
      end if
      ! The sweep accumulates before advance is entered, so it is cleared here,
      ! after the gather, rather than in the reset block at the top.
      lb_rt_accum = 0.0d0
    end if

  end subroutine advance

  !> Advance all grids over one time step, but without taking dimensional
  !> splitting or split source terms into account
  subroutine advect(idim^LIM)
    use mod_global_parameters
    use mod_fix_conserve
    use mod_ghostcells_update
    use mod_comm_lib, only: mpistop

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

    select case (t_stepper)
    case (onestep)
       select case (t_integrator)
       case (Forward_Euler)
          call advect1(flux_method,one,idim^LIM,global_time,ps1,global_time,ps)

       case (IMEX_Euler)
          call advect1(flux_method,one,idim^LIM,global_time,ps,global_time,ps1)
          call global_implicit_update(one,dt,global_time+dt,ps,ps1)

       case (IMEX_SP)
          call global_implicit_update(one,dt,global_time,ps,ps1)
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail; igrid=igrids(iigrid);
             ps1(igrid)%w=ps(igrid)%w
             if(stagger_grid) ps1(igrid)%ws=ps(igrid)%ws
          end do
          !$OMP END PARALLEL DO
          call advect1(flux_method,one,idim^LIM,global_time,ps1,global_time,ps)

       case default
          call mpistop("unkown onestep time_integrator in advect")
       end select

    case (twostep)
       select case (t_integrator)
       case (Predictor_Corrector)
          ! PC or explicit midpoint
          ! predictor step
          fix_conserve_at_step = .false.
          call advect1(typepred1,half,idim^LIM,global_time,ps,global_time,ps1)
          ! corrector step
          fix_conserve_at_step = time_advance .and. levmax>levmin
          call advect1(flux_method,one,idim^LIM,global_time+half*dt,ps1,global_time,ps)

       case (RK2_alf)
          ! RK2 with alfa parameter, where rk_a21=alfa
          call advect1(flux_method,rk_a21, idim^LIM,global_time,ps,global_time,ps1)
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps(igrid)%w = ps(igrid)%w+rk_b1*(ps1(igrid)%w-ps(igrid)%w)/rk_a21
             if(stagger_grid) ps(igrid)%ws = ps(igrid)%ws+(one-rk_b2)*(ps1(igrid)%ws-ps(igrid)%ws)/rk_a21
          end do
          !$OMP END PARALLEL DO
          call advect1(flux_method,rk_b2,idim^LIM,global_time+rk_a21*dt,ps1,global_time+rk_b1*dt,ps)

       case (ssprk2)
          ! ssprk2 or Heun's method
          call advect1(flux_method,one, idim^LIM,global_time,ps,global_time,ps1)
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps(igrid)%w = half*ps(igrid)%w+half*ps1(igrid)%w
             if(stagger_grid) ps(igrid)%ws = half*ps(igrid)%ws+half*ps1(igrid)%ws
          end do
          !$OMP END PARALLEL DO
          call advect1(flux_method,half,idim^LIM,global_time+dt,ps1,global_time+half*dt,ps)

       case (IMEX_Midpoint)
          call advect1(flux_method,half, idim^LIM,global_time,ps,global_time,ps1)
          call global_implicit_update(half,dt,global_time+half*dt,ps2,ps1)
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps(igrid)%w = ps(igrid)%w+2.0d0*(ps2(igrid)%w-ps1(igrid)%w)
             if(stagger_grid) ps(igrid)%ws = ps(igrid)%ws+2.0d0*(ps2(igrid)%ws-ps1(igrid)%ws)
          end do
          !$OMP END PARALLEL DO
          call advect1(flux_method,one, idim^LIM,global_time+half*dt,ps2,global_time,ps)

       case (IMEX_Trapezoidal)
          call advect1(flux_method,one, idim^LIM,global_time,ps,global_time,ps1)
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps2(igrid)%w = half*(ps(igrid)%w+ps1(igrid)%w)
             if(stagger_grid) ps2(igrid)%ws = half*(ps(igrid)%ws+ps1(igrid)%ws)
          end do
          !$OMP END PARALLEL DO
          call evaluate_implicit(global_time,ps)
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps1(igrid)%w = ps1(igrid)%w+half*dt*ps(igrid)%w
             if(stagger_grid) ps1(igrid)%ws = ps1(igrid)%ws+half*dt*ps(igrid)%ws
          end do
          !$OMP END PARALLEL DO
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps(igrid)%w = ps2(igrid)%w+half*dt*ps(igrid)%w
             if(stagger_grid) ps(igrid)%ws = ps2(igrid)%ws+half*dt*ps(igrid)%ws
          end do
          !$OMP END PARALLEL DO
          call getbc(global_time+dt,dt,ps1,iwstart,nwgc)
          call global_implicit_update(half,dt,global_time+dt,ps2,ps1)
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps(igrid)%w = ps(igrid)%w+ps2(igrid)%w-ps1(igrid)%w
             if(stagger_grid) ps(igrid)%ws = ps(igrid)%ws+ps2(igrid)%ws-ps1(igrid)%ws
          end do
          !$OMP END PARALLEL DO
          call advect1(flux_method,half, idim^LIM,global_time+dt,ps2,global_time+half*dt,ps)

       case (IMEX_222)
          ! One-parameter family of schemes (parameter is imex222_lambda) from
          ! Pareschi&Russo 2005, which is L-stable (for default lambda) and
          ! asymptotically SSP.
          ! See doi.org/10.1007/s10915-004-4636-4 (table II)
          ! See doi.org/10.1016/j.apnum.2016.10.018 for interesting values of lambda

          ! Preallocate ps2 as y^n for the implicit update
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps2(igrid)%w = ps(igrid)%w
             if(stagger_grid) ps2(igrid)%ws = ps(igrid)%ws
          end do
          !$OMP END PARALLEL DO
          ! Solve xi1 = y^n + lambda.dt.F_im(xi1)
          call global_implicit_update(imex222_lambda, dt, global_time, ps2, ps)

          ! Set ps1 = y^n + dt.F_ex(xi1)
          call advect1(flux_method, one, idim^LIM, global_time, ps2, global_time, ps1)
          ! Set ps2 = dt.F_im(xi1)        (is at t^n)
          ! Set ps  = y^n + dt/2 . F(xi1) (is at t^n+dt/2)
          ! Set ps1 = y^n + dt.F_ex(xi1) + (1-2.lambda).dt.F_im(xi1) and enforce BC (at t^n+dt)
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps2(igrid)%w = (ps2(igrid)%w - ps(igrid)%w) / imex222_lambda
             if(stagger_grid) ps2(igrid)%ws = (ps2(igrid)%ws - ps(igrid)%ws) / imex222_lambda
          end do
          !$OMP END PARALLEL DO
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps(igrid)%w = half*(ps(igrid)%w + ps1(igrid)%w + ps2(igrid)%w)
             if(stagger_grid) ps(igrid)%ws = half*(ps(igrid)%ws + ps1(igrid)%ws + ps2(igrid)%ws)
          end do
          !$OMP END PARALLEL DO
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps1(igrid)%w = ps1(igrid)%w + (1.0d0 - 2.0d0*imex222_lambda)*ps2(igrid)%w
             if(stagger_grid) ps1(igrid)%ws = ps1(igrid)%ws + (1.0d0 - 2.0d0*imex222_lambda)*ps2(igrid)%ws
          end do
          !$OMP END PARALLEL DO
          call getbc(global_time+dt,dt,ps1,iwstart,nwgc)

          ! Preallocate ps2 as xi1 for the implicit update (is at t^n)
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps2(igrid)%w = 2.0d0*ps2(igrid)%w - ps1(igrid)%w - imex222_lambda*ps2(igrid)%w
             if(stagger_grid) ps2(igrid)%ws = 2.0d0*ps2(igrid)%ws - ps1(igrid)%ws - imex222_lambda*ps2(igrid)%ws
          end do
          !$OMP END PARALLEL DO
          ! Solve xi2 = (ps1) + lambda.dt.F_im(xi2)
          call global_implicit_update(imex222_lambda, dt, global_time, ps2, ps1)

          ! Add dt/2.F_im(xi2) to ps
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps(igrid)%w = ps(igrid)%w + (ps2(igrid)%w - ps1(igrid)%w) / (2.0d0 * imex222_lambda)
             if(stagger_grid) ps(igrid)%ws = ps(igrid)%ws + (ps2(igrid)%ws - ps1(igrid)%ws) / (2.0d0 * imex222_lambda)
          end do
          !$OMP END PARALLEL DO
          ! Set ps = y^n + dt/2.(F(xi1)+F(xi2)) = y^(n+1)
          call advect1(flux_method, half, idim^LIM, global_time+dt, ps2, global_time+half*dt, ps)

       case default
          call mpistop("unkown twostep time_integrator in advect")
       end select

    case (threestep)
       select case (t_integrator)
       case (ssprk3)
          ! this is SSPRK(3,3) Gottlieb-Shu 1998 or SSP(3,2) depending on ssprk_order (3 vs 2)
          call advect1(flux_method,rk_beta11, idim^LIM,global_time,ps,global_time,ps1)
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps2(igrid)%w=rk_alfa21*ps(igrid)%w+rk_alfa22*ps1(igrid)%w
             if(stagger_grid) ps2(igrid)%ws=rk_alfa21*ps(igrid)%ws+rk_alfa22*ps1(igrid)%ws
          end do
          !$OMP END PARALLEL DO
          call advect1(flux_method,rk_beta22, idim^LIM,global_time+rk_c2*dt,ps1,global_time+rk_alfa22*rk_c2*dt,ps2)
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps(igrid)%w=rk_alfa31*ps(igrid)%w+rk_alfa33*ps2(igrid)%w
             if(stagger_grid) ps(igrid)%ws=rk_alfa31*ps(igrid)%ws+rk_alfa33*ps2(igrid)%ws
          end do
          !$OMP END PARALLEL DO
          call advect1(flux_method,rk_beta33, idim^LIM,global_time+rk_c3*dt,ps2,global_time+(1.0d0-rk_beta33)*dt,ps)

       case (RK3_BT)
          ! this is a general threestep RK according to its Butcher Table
          call advect1(flux_method,rk3_a21, idim^LIM,global_time,ps,global_time,ps1)
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps3(igrid)%w=(ps1(igrid)%w-ps(igrid)%w)/rk3_a21
             if(stagger_grid) ps3(igrid)%ws=(ps1(igrid)%ws-ps(igrid)%ws)/rk3_a21
          end do
          !$OMP END PARALLEL DO
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps2(igrid)%w=ps(igrid)%w+rk3_a31*ps3(igrid)%w
             if(stagger_grid) ps2(igrid)%ws=ps(igrid)%ws+rk3_a31*ps3(igrid)%ws
          end do
          !$OMP END PARALLEL DO
          call advect1(flux_method,rk3_a32, idim^LIM,global_time+rk3_c2*dt,ps1,global_time+rk3_a31*dt,ps2)
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps(igrid)%w=ps(igrid)%w+rk3_b1*ps3(igrid)%w &
                  +rk3_b2*(ps2(igrid)%w-(ps(igrid)%w+rk3_a31*ps3(igrid)%w))/rk3_a32
             if(stagger_grid)then
                 ps(igrid)%ws=ps(igrid)%ws+rk3_b1*ps3(igrid)%ws &
                   +rk3_b2*(ps2(igrid)%ws-(ps(igrid)%ws+rk3_a31*ps3(igrid)%ws))/rk3_a32
             endif
          end do
          !$OMP END PARALLEL DO
          call advect1(flux_method,rk3_b3, idim^LIM,global_time+rk3_c3*dt,ps2,global_time+(1.0d0-rk3_b3)*dt,ps)

       case (IMEX_ARS3)
          ! this is IMEX scheme ARS3
          call advect1(flux_method,ars_gamma, idim^LIM,global_time,ps,global_time,ps1)
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps4(igrid)%w=(ps1(igrid)%w-ps(igrid)%w)/ars_gamma
             if(stagger_grid) ps4(igrid)%ws=(ps1(igrid)%ws-ps(igrid)%ws)/ars_gamma
          end do
          !$OMP END PARALLEL DO
          call global_implicit_update(ars_gamma,dt,global_time+ars_gamma*dt,ps2,ps1)
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps1(igrid)%w=(ps2(igrid)%w-ps1(igrid)%w)/ars_gamma
             if(stagger_grid) ps1(igrid)%ws=(ps2(igrid)%ws-ps1(igrid)%ws)/ars_gamma
          end do
          !$OMP END PARALLEL DO
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps3(igrid)%w=ps(igrid)%w+(ars_gamma-1.0d0)*ps4(igrid)%w+(1.0d0-2.0d0*ars_gamma)*ps1(igrid)%w
             if(stagger_grid) then
                ps3(igrid)%ws=ps(igrid)%ws+(ars_gamma-1.0d0)*ps4(igrid)%ws+(1.0d0-2.0d0*ars_gamma)*ps1(igrid)%ws
             endif
          end do
          !$OMP END PARALLEL DO
          call advect1(flux_method,2.0d0*(1.0d0-ars_gamma), idim^LIM,global_time+ars_gamma*dt,ps2,global_time+(ars_gamma-1.0d0)*dt,ps3)
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps2(igrid)%w=ps1(igrid)%w+(ps3(igrid)%w-(ps(igrid)%w+ &
               (ars_gamma-1.0d0)*ps4(igrid)%w+(1.0d0-2.0d0*ars_gamma)*ps1(igrid)%w))/(2.0d0*(1.0d0-ars_gamma))
             if(stagger_grid) then
             ps2(igrid)%ws=ps1(igrid)%ws+(ps3(igrid)%ws-(ps(igrid)%ws+ &
               (ars_gamma-1.0d0)*ps4(igrid)%ws+(1.0d0-2.0d0*ars_gamma)*ps1(igrid)%ws))/(2.0d0*(1.0d0-ars_gamma))
             endif
          end do
          !$OMP END PARALLEL DO
          call global_implicit_update(ars_gamma,dt,global_time+(1.0d0-ars_gamma)*dt,ps4,ps3)
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps(igrid)%w=ps(igrid)%w+half*ps2(igrid)%w &
                +half*(ps4(igrid)%w-ps3(igrid)%w)/ars_gamma
             if(stagger_grid) then
                ps(igrid)%ws=ps(igrid)%ws+half*ps2(igrid)%ws &
                    +half*(ps4(igrid)%ws-ps3(igrid)%ws)/ars_gamma
             endif
          end do
          !$OMP END PARALLEL DO
          call advect1(flux_method,half, idim^LIM,global_time+(1.0d0-ars_gamma)*dt,ps4,global_time+half*dt,ps)

       case (IMEX_232)
          ! this is IMEX_ARK(2,3,2) or IMEX_SSP(2,3,2)
          call advect1(flux_method,imex_a21, idim^LIM,global_time,ps,global_time,ps1)
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps4(igrid)%w=(ps1(igrid)%w-ps(igrid)%w)/imex_a21
             ps3(igrid)%w=ps(igrid)%w
             if(stagger_grid) then
               ps4(igrid)%ws=(ps1(igrid)%ws-ps(igrid)%ws)/imex_a21
               ps3(igrid)%ws=ps(igrid)%ws
             endif
          end do
          !$OMP END PARALLEL DO
          call evaluate_implicit(global_time,ps3)
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps1(igrid)%w=ps1(igrid)%w+imex_ha21*dt*ps3(igrid)%w
             if(stagger_grid) ps1(igrid)%ws=ps1(igrid)%ws+imex_ha21*dt*ps3(igrid)%ws
          end do
          !$OMP END PARALLEL DO
          call getbc(global_time+imex_a21*dt,dt,ps1,iwstart,nwgc)
          call global_implicit_update(imex_ha22,dt,global_time+imex_c2*dt,ps2,ps1)
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps(igrid)%w=ps(igrid)%w+imex_a31*ps4(igrid)%w &
                +imex_b1*dt*ps3(igrid)%w+imex_b2*(ps2(igrid)%w-ps1(igrid)%w)/imex_ha22
             if(stagger_grid) then
             ps(igrid)%ws=ps(igrid)%ws+imex_a31*ps4(igrid)%ws &
                +imex_b1*dt*ps3(igrid)%ws+imex_b2*(ps2(igrid)%ws-ps1(igrid)%ws)/imex_ha22
             endif
          end do
          !$OMP END PARALLEL DO
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps3(igrid)%w=ps1(igrid)%w-imex_a21*ps4(igrid)%w &
                -imex_ha21*dt*ps3(igrid)%w+imex_b1*dt*ps3(igrid)%w
             if(stagger_grid) then
             ps3(igrid)%ws=ps1(igrid)%ws-imex_a21*ps4(igrid)%ws &
                -imex_ha21*dt*ps3(igrid)%ws+imex_b1*dt*ps3(igrid)%ws
             endif
          end do
          !$OMP END PARALLEL DO
          call advect1(flux_method,imex_a32, idim^LIM,global_time+imex_c2*dt,ps2,global_time+imex_a31*dt,ps)
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps2(igrid)%w=(ps(igrid)%w-ps3(igrid)%w-imex_a31*ps4(igrid)%w)/imex_a32 &
                +(1.0d0-imex_b2/imex_a32)*(ps2(igrid)%w-ps1(igrid)%w)/imex_ha22
             if(stagger_grid) then
             ps2(igrid)%ws=(ps(igrid)%ws-ps3(igrid)%ws-imex_a31*ps4(igrid)%ws)/imex_a32 &
                +(1.0d0-imex_b2/imex_a32)*(ps2(igrid)%ws-ps1(igrid)%ws)/imex_ha22
             endif
          end do
          !$OMP END PARALLEL DO
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps1(igrid)%w=ps3(igrid)%w+imex_b1*ps4(igrid)%w+imex_b2*ps2(igrid)%w
             if(stagger_grid) then
             ps1(igrid)%ws=ps3(igrid)%ws+imex_b1*ps4(igrid)%ws+imex_b2*ps2(igrid)%ws
             endif
          end do
          !$OMP END PARALLEL DO
          call global_implicit_update(imex_b3,dt,global_time+imex_c3*dt,ps2,ps)
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps(igrid)%w=ps1(igrid)%w+ps2(igrid)%w-ps(igrid)%w
             if(stagger_grid) then
             ps(igrid)%ws=ps1(igrid)%ws+ps2(igrid)%ws-ps(igrid)%ws
             endif
          end do
          !$OMP END PARALLEL DO
          call advect1(flux_method,imex_b3, idim^LIM,global_time+imex_c3*dt,ps2,global_time+(1.0d0-imex_b3)*dt,ps)

       case (IMEX_CB3a)
          ! Third order IMEX scheme with low-storage implementation (4 registers).
          ! From Cavaglieri&Bewley 2015, see doi.org/10.1016/j.jcp.2015.01.031
          ! (scheme called "IMEXRKCB3a" there). Uses 3 explicit and 2 implicit stages.
          ! Parameters are in imex_bj, imex_cj (same for implicit/explicit),
          ! imex_aij (implicit tableau) and imex_haij (explicit tableau).
          call advect1(flux_method, imex_ha21, idim^LIM, global_time, ps, global_time, ps1)
          call global_implicit_update(imex_a22, dt, global_time+imex_c2*dt, ps2, ps1)
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps3(igrid)%w = ps(igrid)%w + imex_a32/imex_a22 * (ps2(igrid)%w - ps1(igrid)%w)
             ps(igrid)%w  = ps(igrid)%w + imex_b2 /imex_a22 * (ps2(igrid)%w - ps1(igrid)%w)
             ps1(igrid)%w = ps3(igrid)%w
             if(stagger_grid) ps3(igrid)%ws = ps(igrid)%ws + imex_a32/imex_a22 * (ps2(igrid)%ws - ps1(igrid)%ws)
             if(stagger_grid) ps(igrid)%ws  = ps(igrid)%ws + imex_b2 /imex_a22 * (ps2(igrid)%ws - ps1(igrid)%ws)
             if(stagger_grid) ps1(igrid)%ws = ps3(igrid)%ws
          end do
          !$OMP END PARALLEL DO
          call advect1(flux_method, imex_ha32, idim^LIM, global_time+imex_c2*dt, ps2, global_time, ps3)
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps(igrid)%w = ps(igrid)%w + imex_b2 /imex_ha32 * (ps3(igrid)%w - ps1(igrid)%w)
             if(stagger_grid) ps(igrid)%ws = ps(igrid)%ws + imex_b2 /imex_ha32 * (ps3(igrid)%ws - ps1(igrid)%ws)
          end do
          call global_implicit_update(imex_a33, dt, global_time+imex_c3*dt, ps1, ps3)
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps(igrid)%w = ps(igrid)%w + imex_b3 /imex_a33 * (ps1(igrid)%w - ps3(igrid)%w)
             if(stagger_grid) ps(igrid)%ws = ps(igrid)%ws + imex_b3 /imex_a33 * (ps1(igrid)%ws - ps3(igrid)%ws)
          end do
          !$OMP END PARALLEL DO
          call advect1(flux_method, imex_b3, idim^LIM, global_time+imex_c3*dt, ps1, global_time+imex_b2*dt, ps)

       case default
          call mpistop("unkown threestep time_integrator in advect")
       end select

    case (fourstep)
       select case (t_integrator)
       case (ssprk4)
          ! SSPRK(4,3) or SSP(4,2) depending on ssprk_order (3 vs 2)
          ! ssprk43: Strong stability preserving 4 stage RK 3rd order by Ruuth and Spiteri
          !    Ruuth & Spiteri J. S C, 17 (2002) p. 211 - 220
          !    supposed to be stable up to CFL=2.
          ! ssp42: stable up to CFL=3
          call advect1(flux_method,rk_beta11, idim^LIM,global_time,ps,global_time,ps1)
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps2(igrid)%w=rk_alfa21*ps(igrid)%w+rk_alfa22*ps1(igrid)%w
             if(stagger_grid) ps2(igrid)%ws=rk_alfa21*ps(igrid)%ws+rk_alfa22*ps1(igrid)%ws
          end do
          !$OMP END PARALLEL DO
          call advect1(flux_method,rk_beta22, idim^LIM,global_time+rk_c2*dt,ps1,global_time+rk_alfa22*rk_c2*dt,ps2)
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps1(igrid)%w=rk_alfa31*ps(igrid)%w+rk_alfa33*ps2(igrid)%w
             if(stagger_grid) ps1(igrid)%ws=rk_alfa31*ps(igrid)%ws+rk_alfa33*ps2(igrid)%ws
          end do
          !$OMP END PARALLEL DO
          call advect1(flux_method,rk_beta33, idim^LIM,global_time+rk_c3*dt,ps2,global_time+rk_alfa33*rk_c3*dt,ps1)
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps(igrid)%w=rk_alfa41*ps(igrid)%w+rk_alfa44*ps1(igrid)%w
             if(stagger_grid) ps(igrid)%ws=rk_alfa41*ps(igrid)%ws+rk_alfa44*ps1(igrid)%ws
          end do
          !$OMP END PARALLEL DO
          call advect1(flux_method,rk_beta44, idim^LIM,global_time+rk_c4*dt,ps1,global_time+(1.0d0-rk_beta44)*dt,ps)

       case (rk4)
          ! the standard RK(4,4) method
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps2(igrid)%w=ps(igrid)%w
             ps3(igrid)%w=ps(igrid)%w
             if(stagger_grid) then
                ps2(igrid)%ws=ps(igrid)%ws
                ps3(igrid)%ws=ps(igrid)%ws
             endif
          end do
          !$OMP END PARALLEL DO
          call advect1(flux_method,half, idim^LIM,global_time,ps,global_time,ps1)
          call advect1(flux_method,half, idim^LIM,global_time+half*dt,ps1,global_time,ps2)
          call advect1(flux_method,1.0d0, idim^LIM,global_time+half*dt,ps2,global_time,ps3)
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps(igrid)%w=(1.0d0/3.0d0)*(-ps(igrid)%w+ps1(igrid)%w+2.0d0*ps2(igrid)%w+ps3(igrid)%w)
             if(stagger_grid) ps(igrid)%ws=(1.0d0/3.0d0) &
                 *(-ps(igrid)%ws+ps1(igrid)%ws+2.0d0*ps2(igrid)%ws+ps3(igrid)%ws)
          end do
          !$OMP END PARALLEL DO
          call advect1(flux_method,1.0d0/6.0d0, idim^LIM,global_time+dt,ps3,global_time+dt*5.0d0/6.0d0,ps)

       case default
          call mpistop("unkown fourstep time_integrator in advect")
       end select

    case (fivestep)
       select case (t_integrator)
       case (ssprk5)
          ! SSPRK(5,4) by Ruuth and Spiteri
          !bcexch = .false.
          call advect1(flux_method,rk_beta11, idim^LIM,global_time,ps,global_time,ps1)
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps2(igrid)%w=rk_alfa21*ps(igrid)%w+rk_alfa22*ps1(igrid)%w
             if(stagger_grid) ps2(igrid)%ws=rk_alfa21*ps(igrid)%ws+rk_alfa22*ps1(igrid)%ws
          end do
          !$OMP END PARALLEL DO
          call advect1(flux_method,rk_beta22, idim^LIM,global_time+rk_c2*dt,ps1,global_time+rk_alfa22*rk_c2*dt,ps2)
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps1(igrid)%w=rk_alfa31*ps(igrid)%w+rk_alfa33*ps2(igrid)%w
             if(stagger_grid) ps1(igrid)%ws=rk_alfa31*ps(igrid)%ws+rk_alfa33*ps2(igrid)%ws
          end do
          !$OMP END PARALLEL DO
          call advect1(flux_method,rk_beta33, idim^LIM,global_time+rk_c3*dt,ps2,global_time+rk_alfa33*rk_c3*dt,ps1)
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps3(igrid)%w=rk_alfa53*ps2(igrid)%w+rk_alfa54*ps1(igrid)%w
             if(stagger_grid) ps3(igrid)%ws=rk_alfa53*ps2(igrid)%ws+rk_alfa54*ps1(igrid)%ws
          end do
          !$OMP END PARALLEL DO
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps2(igrid)%w=rk_alfa41*ps(igrid)%w+rk_alfa44*ps1(igrid)%w
             if(stagger_grid) ps2(igrid)%ws=rk_alfa41*ps(igrid)%ws+rk_alfa44*ps1(igrid)%ws
          end do
          !$OMP END PARALLEL DO
          call advect1(flux_method,rk_beta44, idim^LIM,global_time+rk_c4*dt,ps1,global_time+rk_alfa44*rk_c4*dt,ps2)
          !$OMP PARALLEL DO PRIVATE(igrid)
          do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
             ps(igrid)%w=ps3(igrid)%w+rk_alfa55*ps2(igrid)%w &
                +(rk_beta54/rk_beta44)*(ps2(igrid)%w-(rk_alfa41*ps(igrid)%w+rk_alfa44*ps1(igrid)%w))
             if(stagger_grid) then
             ps(igrid)%ws=ps3(igrid)%ws+rk_alfa55*ps2(igrid)%ws &
                +(rk_beta54/rk_beta44)*(ps2(igrid)%ws-(rk_alfa41*ps(igrid)%ws+rk_alfa44*ps1(igrid)%ws))
             endif
          end do
          !$OMP END PARALLEL DO
          !bcexch = .true.
          call advect1(flux_method,rk_beta55, idim^LIM,global_time+rk_c5*dt,ps2,global_time+(1.0d0-rk_beta55)*dt,ps)

       case default
          call mpistop("unkown fivestep time_integrator in advect")
       end select

    case default
       call mpistop("unkown time_stepper in advect")
    end select

  end subroutine advect

  !> Implicit global update step within IMEX schemes, advance psa=psb+dtfactor*qdt*F_im(psa)
  subroutine global_implicit_update(dtfactor,qdt,qtC,psa,psb)
    use mod_global_parameters
    use mod_ghostcells_update
    use mod_physics, only: phys_implicit_update

    type(state), target :: psa(max_blocks)   !< Compute implicit part from this state and update it
    type(state), target :: psb(max_blocks)   !< Will be unchanged, as on entry
    double precision, intent(in) :: qdt      !< overall time step dt
    double precision, intent(in) :: qtC      !< Both states psa and psb at this time level
    double precision, intent(in) :: dtfactor !< Advance psa=psb+dtfactor*qdt*F_im(psa)

    integer                        :: iigrid, igrid

    !> First copy all variables from a to b, this is necessary to account for
    ! quantities in w with no implicit sourceterm
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       psa(igrid)%w = psb(igrid)%w
       if(stagger_grid) psa(igrid)%ws = psb(igrid)%ws
    end do

    if (associated(phys_implicit_update)) then
       call phys_implicit_update(dtfactor,qdt,qtC,psa,psb)
    end if

    ! enforce boundary conditions for psa
    call getbc(qtC,0.d0,psa,iwstart,nwgc)

  end subroutine global_implicit_update

  !> Evaluate Implicit part in place, i.e. psa==>F_im(psa)
  subroutine evaluate_implicit(qtC,psa)
    use mod_global_parameters
    use mod_physics, only: phys_evaluate_implicit
    type(state), target :: psa(max_blocks)   !< Compute implicit part from this state and update it
    double precision, intent(in) :: qtC      !< psa at this time level

    if (associated(phys_evaluate_implicit)) then
       call phys_evaluate_implicit(qtC,psa)
    end if
  end subroutine evaluate_implicit

  !> Integrate all grids by one partial step
  subroutine advect1(method,dtfactor,idim^LIM,qtC,psa,qt,psb)
    use mod_global_parameters
    use mod_ghostcells_update
    use mod_fix_conserve
    use mod_physics
    use mod_eos

    integer, intent(in) :: idim^LIM
    type(state), target :: psa(max_blocks) !< Compute fluxes based on this state
    type(state), target :: psb(max_blocks) !< Update solution on this state
    double precision, intent(in) :: dtfactor !< Advance over dtfactor * dt
    double precision, intent(in) :: qtC
    double precision, intent(in) :: qt
    integer, intent(in) :: method(nlevelshi)

    double precision :: qdt
    double precision :: lb_t0_advect1, lb_t0_block
    integer :: iigrid, igrid

    istep = istep+1

    if(associated(phys_special_advance)) then
      call phys_special_advance(qtC,psa)
    end if

    qdt=dtfactor*dt
    ! opedit: Just advance the active grids:
    if (lb_diagnose) lb_t0_advect1 = MPI_WTIME()
    !$OMP PARALLEL DO PRIVATE(igrid,lb_t0_block)
    do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
      if (lb_automatic) lb_t0_block = MPI_WTIME()
      block=>ps(igrid)
      ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
      call eos%update_eos(ixG^LL,ixG^LL,psa(igrid)%w,psa(igrid)%x)
      call advect1_grid(igrid,method(block%level),qdt,dtfactor,ixG^LL,idim^LIM,&
        qtC,psa(igrid),qt,psb(igrid),rnode(rpdx1_:rnodehi,igrid),ps(igrid)%x)
      if (lb_automatic) block_cost(igrid) = block_cost(igrid) + (MPI_WTIME() - lb_t0_block)
    end do
    !$OMP END PARALLEL DO
    if (lb_diagnose) lb_compute_accum = lb_compute_accum + (MPI_WTIME() - lb_t0_advect1)

    ! opedit: Send flux for all grids, expects sends for all
    ! nsend_fc(^D), set in connectivity.t.

    if (fix_conserve_global .and. fix_conserve_at_step) then
      call recvflux(idim^LIM)
      call sendflux(idim^LIM)
      call fix_conserve(psb,idim^LIM,1,nwflux)
      if(stagger_grid) then
        call fix_edges(psb,idim^LIM)
        ! fill the cell-center values from the updated staggered variables
        !$OMP PARALLEL DO PRIVATE(igrid)
        do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
          call phys_face_to_center(ixM^LL,psb(igrid))
        end do
        !$OMP END PARALLEL DO
      end if
    end if

    ! For all grids: fill ghost cells
    call getbc(qt+qdt,qdt,psb,iwstart,nwgc)

  end subroutine advect1

  !> Advance a single grid over one partial time step
  subroutine advect1_grid(igrid,method,qdt,dtfactor,ixI^L,idim^LIM,qtC,sCT,qt,s,dxs,x)

    !  integrate one grid by one partial step
    use mod_finite_volume
    use mod_finite_difference
    use mod_tvd
    use mod_source, only: addsource2
    use mod_physics, only: phys_to_primitive
    use mod_global_parameters
    use mod_comm_lib, only: mpistop
    use mod_fix_conserve

    integer, intent(in) :: igrid,method
    integer, intent(in) :: ixI^L, idim^LIM
    double precision, intent(in) :: qdt, dtfactor, qtC, qt, dxs(ndim), x(ixI^S,1:ndim)
    type(state), target          :: sCT, s

    ! cell face flux
    double precision :: fC(ixI^S,1:nwflux,1:ndim)
    ! cell edge flux
    double precision :: fE(ixI^S,sdim:3)
    double precision :: wprim(ixI^S,1:nw)
    integer :: ixO^L

    ! for mf module
    if(iwstart>1) fC=0.d0

    ixO^L=ixI^L^LSUBnghostcells;
    select case (method)
    case (fs_hll,fs_hllc,fs_hllcd,fs_hlld,fs_tvdlf,fs_tvdmu)
       call finite_volume(method,qdt,dtfactor,ixI^L,ixO^L,idim^LIM,qtC,sCT,qt,s,fC,fE,dxs,x)
    case (fs_cd,fs_cd4)
       call centdiff(method,qdt,dtfactor,ixI^L,ixO^L,idim^LIM,qtC,sCT,qt,s,fC,fE,dxs,x)
    case (fs_hancock)
       call hancock(qdt,dtfactor,ixI^L,ixO^L,idim^LIM,qtC,sCT,qt,s,dxs,x)
    case (fs_fd)
       call fd(qdt,dtfactor,ixI^L,ixO^L,idim^LIM,qtC,sCT,qt,s,fC,fE,dxs,x)
    case (fs_tvd)
       call centdiff(fs_cd,qdt,dtfactor,ixI^L,ixO^L,idim^LIM,qtC,sCT,qt,s,fC,fE,dxs,x)
       call tvdlimit(method,qdt,ixI^L,ixO^L,idim^LIM,sCT,qt+qdt,s,fC,dxs,x)
    case (fs_source)
       fC=0.d0
       fE=0.d0
       wprim=sCT%w
       call phys_to_primitive(ixI^L,ixI^L,wprim,x)
       call addsource2(qdt*dble(idimmax-idimmin+1)/dble(ndim),&
            dtfactor*dble(idimmax-idimmin+1)/dble(ndim),&
            ixI^L,ixO^L,1,nw,qtC,sCT%w,wprim,qt,s%w,x,.false.)
    case (fs_nul)
       ! There is nothing to do
    case default
       call mpistop("unknown flux scheme in advect1_grid")
    end select

    ! opedit: Obviously, flux is stored only for active grids.
    ! but we know in fix_conserve wether there is a passive neighbor
    ! but we know in conserve_fix wether there is a passive neighbor
    ! via neighbor_active(i^D,igrid) thus we skip the correction for those.
    ! This violates strict conservation when the active/passive interface
    ! coincides with a coarse/fine interface.
    if (fix_conserve_global .and. fix_conserve_at_step) then
      call store_flux(igrid,fC,idim^LIM,nwflux)
      if(stagger_grid) call store_edge(igrid,ixG^LL,fE,idim^LIM)
    end if

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

    integer:: iigrid, igrid

    if (associated(usr_process_global)) then
       call usr_process_global(iit,qt)
    end if

    if (associated(usr_process_grid)) then
      !$OMP PARALLEL DO PRIVATE(igrid)
      do iigrid=1,igridstail; igrid=igrids(iigrid);
         ! next few lines ensure correct usage of routines like divvector etc
         ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
         block=>ps(igrid)
         call usr_process_grid(igrid,node(plevel_,igrid),ixG^LL,ixM^LL, &
              qt,ps(igrid)%w,ps(igrid)%x)
      end do
      !$OMP END PARALLEL DO
      call getbc(qt,dt,ps,iwstart,nwgc)
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

    integer:: iigrid, igrid

    if (associated(usr_process_adv_global)) then
       call usr_process_adv_global(iit,qt)
    end if

    if (associated(usr_process_adv_grid)) then
      !$OMP PARALLEL DO PRIVATE(igrid)
      do iigrid=1,igridstail; igrid=igrids(iigrid);
         ! next few lines ensure correct usage of routines like divvector etc
         ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
         block=>ps(igrid)

         call usr_process_adv_grid(igrid,node(plevel_,igrid),ixG^LL,ixM^LL, &
              qt,ps(igrid)%w,ps(igrid)%x)
      end do
      !$OMP END PARALLEL DO
      call getbc(qt,dt,ps,iwstart,nwgc)
    end if
  end subroutine process_advanced

end module mod_advance
