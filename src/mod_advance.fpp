#:mute
#:include "mod_gpu_directives.fpp"
#:endmute
!> Module containing all the time stepping schemes
!>
!> advect implements every explicit time integrator of MPI-AMRVAC, selected by
!> time_stepper and time_integrator in &methodlist: Forward_Euler;
!> Predictor_Corrector, ssprk2 and RK2_alfa; ssprk3 and RK3_BT; ssprk4 and
!> rk4; ssprk5.  Each is a sequence of advect1 stages, which add dtfactor*dt
!> times the spatial operator on one time level of bg to another, and of
!> bg_assign/bg_accum combinations of the levels.  How many levels a scheme
!> holds, and so how many bg is allocated with, is time_levels in
!> mod_global_parameters.  See the "Time integration" page for the schemes,
!> their stability limits and memory cost.
module mod_advance

  implicit none
  private

  !> Whether to conserve fluxes at the current sub-step
  logical :: fix_conserve_at_step = .true.

  public :: advance
  public :: process
  public :: process_advanced

  !> bgd = a1*bg1 [+ a2*bg2 [+ a3*bg3]]: a pointwise combination of time
  !> levels between stages, generated below for one to three sources
  interface bg_assign
     module procedure bg_assign1, bg_assign2, bg_assign3
  end interface bg_assign

  !> bgd = s*bgd + a1*bg1 [+ a2*bg2 [+ a3*bg3]]: the in-place form of
  !> bg_assign, generated below for one to three sources
  interface bg_accum
     module procedure bg_accum1, bg_accum2, bg_accum3
  end interface bg_accum

contains

  !> Advance all the grids over one time step, including all sources
  subroutine advance(iit)
    use mod_global_parameters
    use mod_particles, only: handle_particles
    use mod_source, only: add_split_source
    use mod_set_aux, only: set_aux_all

    integer, intent(in) :: iit

    integer :: iigrid, igrid, idimsplit
    
    ! split source addition
    call add_split_source(prior=.true.)

    if (dimsplit) then
       if ((iit/2)*2==iit .or. typedimsplit=='xy') then
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
    call add_split_source(prior=.false.)

    ! the update only wrote the flux variables: refresh the auxiliaries of the
    ! final state, once per step however many sweeps and split sources it took
    ! (needs the ghost cells the last stage filled)
    call set_aux_all()

    if(use_particles) call handle_particles

    ! do iigrid=1,igridstail; igrid=igrids(iigrid);
    !    ${GPU_EXIT_DATA_DELETE('ps(igrid)%x, ps1(igrid)%w, ps2(igrid)%w) copyout(ps(igrid)%w')}$
    ! end do
    
  end subroutine advance

  !> Advance all grids over one time step, but without taking dimensional
  !> splitting or split source terms into account
  subroutine advect(idimmin,idimmax)
    use mod_global_parameters
    use mod_fix_conserve
    use mod_ghostcells_update
    use mod_physics, only: phys_req_diagonal
    use mod_comm_lib, only: mpistop
    use mod_physicaldata

    integer, intent(in) :: idimmin,idimmax
    integer             :: iigrid, igrid, ix1,ix2,ix3, iw

    call init_comm_fix_conserve(idimmin,idimmax,nwflux)
    fix_conserve_at_step = time_advance .and. levmax>levmin
    
    ! copy w instead of wold because of potential use of dimsplit or sourcesplit
    ${GPU_PARALLEL_LOOP_GANG("private(igrid)")}$
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       ${GPU_LOOP_VECTOR("collapse(ndim+1)")}$
       do iw = 1, nw
          do ix3 = ixGlo3, ixGhi3 
             do ix2 = ixGlo2, ixGhi2 
                do ix1 = ixGlo1, ixGhi1 
                   bg(2)%w(ix1,ix2,ix3,iw,igrid) = bg(1)%w(ix1,ix2,ix3,iw,igrid)
                end do
             end do
          end do
       end do
    end do
    
    if(stagger_grid) then
       ${GPU_PARALLEL_LOOP_GANG("private(igrid)")}$
       do iigrid=1,igridstail; igrid=igrids(iigrid);
          ${GPU_LOOP_VECTOR("collapse(ndim+1)")}$
          do iw = 1, nws
             do ix3 = ps(igrid)%ixGsmin3, ps(igrid)%ixGsmax3 
                do ix2 = ps(igrid)%ixGsmin2, ps(igrid)%ixGsmax2 
                   do ix1 = ps(igrid)%ixGsmin1, ps(igrid)%ixGsmax1 
                      ps1(igrid)%ws(ix1,ix2,ix3,iw) = ps(igrid)%ws(ix1,ix2,ix3,iw)
                   end do
                end do
             end do
          end do
       end do
    end if

    istep = 0

    ! Each scheme is a sequence of advect1 stages and pointwise combinations
    ! of the time levels, transcribed from upstream MPI-AMRVAC's advect with
    ! ps -> bg(1), ps1 -> bg(2), ps2 -> bg(3).  How many levels a scheme may
    ! use is fixed by time_levels in mod_global_parameters.
    select case (t_stepper)

    case (onestep)
       select case (t_integrator)
       case (Forward_Euler)
          call advect1(1.0d0, idimmin,idimmax, global_time,2, global_time,1)

       case default
          call mpistop("unkown onestep time_integrator in advect")
       end select

    case (twostep)
       select case (t_integrator)
       case (Predictor_Corrector)
          ! PC or explicit midpoint; the predictor is not refluxed
          fix_conserve_at_step = .false.
          call advect1(0.5d0, idimmin,idimmax, global_time,1, global_time,2)
          fix_conserve_at_step = time_advance .and. levmax>levmin
          call advect1(1.0d0, idimmin,idimmax, global_time+0.5d0*dt,2,&
               global_time,1)

       case (RK2_alf)
          ! RK2 with alfa parameter, where rk_a21=alfa
          call advect1(rk_a21, idimmin,idimmax, global_time,1, global_time,2)
          call bg_accum(bg(1), 1.0d0-rk_b1/rk_a21, rk_b1/rk_a21, bg(2))
          call advect1(rk_b2, idimmin,idimmax, global_time+rk_a21*dt,2,&
               global_time+rk_b1*dt,1)

       case (ssprk2)
          ! ssprk2 or Heun's method
          call advect1(1.0d0, idimmin,idimmax, global_time,1, global_time,2)
          call bg_accum(bg(1), 0.5d0, 0.5d0, bg(2))
          call advect1(0.5d0, idimmin,idimmax, global_time+dt,2,&
               global_time+0.5d0*dt,1)

       case default
          call mpistop("unkown twostep time_integrator in advect")
       end select

    case (threestep)
       select case (t_integrator)
       case (ssprk3)
          ! this is SSPRK(3,3) Gottlieb-Shu 1998 or SSP(3,2) depending on ssprk_order (3 vs 2)
          call advect1(rk_beta11, idimmin,idimmax, global_time,1,&
               global_time,2)
          call bg_assign(bg(3), rk_alfa21, bg(1), rk_alfa22, bg(2))
          call advect1(rk_beta22, idimmin,idimmax, global_time+rk_c2*dt,2,&
               global_time+rk_alfa22*rk_c2*dt,3)
          call bg_accum(bg(1), rk_alfa31, rk_alfa33, bg(3))
          call advect1(rk_beta33, idimmin,idimmax, global_time+rk_c3*dt,3,&
               global_time+(1.0d0-rk_beta33)*dt,1)

       case (RK3_BT)
          ! general threestep RK from its Butcher table (rk3_switch), with
          ! level 4 holding dt*k1
          call advect1(rk3_a21, idimmin,idimmax, global_time,1, global_time,2)
          call bg_assign(bg(4), 1.0d0/rk3_a21, bg(2), -1.0d0/rk3_a21, bg(1))
          call bg_assign(bg(3), 1.0d0, bg(1), rk3_a31, bg(4))
          call advect1(rk3_a32, idimmin,idimmax, global_time+rk3_c2*dt,2,&
               global_time+rk3_a31*dt,3)
          ! ps = ps + b1*ps3 + b2*(ps2 - (ps + a31*ps3))/a32, upstream's form
          call bg_accum(bg(1), 1.0d0-rk3_b2/rk3_a32,&
               rk3_b1-rk3_b2*rk3_a31/rk3_a32, bg(4), rk3_b2/rk3_a32, bg(3))
          call advect1(rk3_b3, idimmin,idimmax, global_time+rk3_c3*dt,3,&
               global_time+(1.0d0-rk3_b3)*dt,1)

       case default
          call mpistop("unkown threestep time_integrator in advect")
       end select

    case (fourstep)
       select case (t_integrator)
       case (ssprk4)
          ! SSPRK(4,3) by Spiteri & Ruuth, or SSP(4,2), depending on
          ! ssprk_order (3 vs 2); stable up to CFL 2 and 3 respectively
          call advect1(rk_beta11, idimmin,idimmax, global_time,1,&
               global_time,2)
          call bg_assign(bg(3), rk_alfa21, bg(1), rk_alfa22, bg(2))
          call advect1(rk_beta22, idimmin,idimmax, global_time+rk_c2*dt,2,&
               global_time+rk_alfa22*rk_c2*dt,3)
          call bg_assign(bg(2), rk_alfa31, bg(1), rk_alfa33, bg(3))
          call advect1(rk_beta33, idimmin,idimmax, global_time+rk_c3*dt,3,&
               global_time+rk_alfa33*rk_c3*dt,2)
          call bg_accum(bg(1), rk_alfa41, rk_alfa44, bg(2))
          call advect1(rk_beta44, idimmin,idimmax, global_time+rk_c4*dt,2,&
               global_time+(1.0d0-rk_beta44)*dt,1)

       case (rk4)
          ! the classical RK(4,4); levels 2, 3, 4 become u+dt/2*k1,
          ! u+dt/2*k2, u+dt*k3
          call bg_assign(bg(3), 1.0d0, bg(1))
          call bg_assign(bg(4), 1.0d0, bg(1))
          call advect1(0.5d0, idimmin,idimmax, global_time,1, global_time,2)
          call advect1(0.5d0, idimmin,idimmax, global_time+0.5d0*dt,2,&
               global_time,3)
          call advect1(1.0d0, idimmin,idimmax, global_time+0.5d0*dt,3,&
               global_time,4)
          ! ps = (-ps + ps1 + 2*ps2 + ps3)/3
          call bg_accum(bg(1), -1.0d0/3.0d0, 1.0d0/3.0d0, bg(2),&
               2.0d0/3.0d0, bg(3), 1.0d0/3.0d0, bg(4))
          call advect1(1.0d0/6.0d0, idimmin,idimmax, global_time+dt,4,&
               global_time+dt*5.0d0/6.0d0,1)

       case default
          call mpistop("unkown fourstep time_integrator in advect")
       end select

    case (fivestep)
       select case (t_integrator)
       case (ssprk5)
          ! SSPRK(5,4) by Spiteri & Ruuth
          call advect1(rk_beta11, idimmin,idimmax, global_time,1,&
               global_time,2)
          call bg_assign(bg(3), rk_alfa21, bg(1), rk_alfa22, bg(2))
          call advect1(rk_beta22, idimmin,idimmax, global_time+rk_c2*dt,2,&
               global_time+rk_alfa22*rk_c2*dt,3)
          call bg_assign(bg(2), rk_alfa31, bg(1), rk_alfa33, bg(3))
          call advect1(rk_beta33, idimmin,idimmax, global_time+rk_c3*dt,3,&
               global_time+rk_alfa33*rk_c3*dt,2)
          call bg_assign(bg(4), rk_alfa53, bg(3), rk_alfa54, bg(2))
          call bg_assign(bg(3), rk_alfa41, bg(1), rk_alfa44, bg(2))
          call advect1(rk_beta44, idimmin,idimmax, global_time+rk_c4*dt,2,&
               global_time+rk_alfa44*rk_c4*dt,3)
          ! ps = ps3 + alfa55*ps2 + (beta54/beta44)*(ps2 - (alfa41*ps +
          ! alfa44*ps1)), upstream's form, with the bracket expanded
          call bg_accum(bg(1), -rk_alfa41*rk_beta54/rk_beta44, 1.0d0, bg(4),&
               rk_alfa55+rk_beta54/rk_beta44, bg(3),&
               -rk_alfa44*rk_beta54/rk_beta44, bg(2))
          call advect1(rk_beta55, idimmin,idimmax, global_time+rk_c5*dt,3,&
               global_time+(1.0d0-rk_beta55)*dt,1)

       case default
          call mpistop("unkown fivestep time_integrator in advect")
       end select

    case default
       call mpistop("unkown time_stepper in advect")
    end select

  end subroutine advect

#:for n in [1, 2, 3]
  #:set args = ', '.join('a{0}, bg{0}'.format(i) for i in range(1, n+1))
  #:set coefs = ', '.join('a{0}'.format(i) for i in range(1, n+1))
  #:set srcs = ', '.join('bg{0}'.format(i) for i in range(1, n+1))
  #:set sum = ' &\n                        + '.join('a{0} * bg{0}%w(ix1,ix2,ix3,iw,igrid)'.format(i) for i in range(1, n+1))
  !> bgd = sum of ${n}$ scaled time levels, on every active block, ghost
  !> cells included.  The destination is never read, so it may hold anything
  !> beforehand.
  subroutine bg_assign${n}$(bgd, ${args}$)
    use mod_global_parameters

    type(block_grid_t), intent(inout) :: bgd
    double precision, intent(in)      :: ${coefs}$
    type(block_grid_t), intent(in)    :: ${srcs}$

    integer :: iigrid, igrid, ix1,ix2,ix3, iw

    ${GPU_PARALLEL_LOOP_GANG("private(igrid)")}$
    do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
       ${GPU_LOOP_VECTOR("collapse(ndim+1)")}$
       do iw = 1, nw
          do ix3 = ixGlo3, ixGhi3
             do ix2 = ixGlo2, ixGhi2
                do ix1 = ixGlo1, ixGhi1
                   bgd%w(ix1,ix2,ix3,iw,igrid) = ${sum}$
                end do
             end do
          end do
       end do
    end do

  end subroutine bg_assign${n}$

  !> bgd = s*bgd plus the sum of ${n}$ scaled time levels, on every active
  !> block, ghost cells included.  A separate shape from bg_assign so that
  !> the destination is never also passed as a source.
  subroutine bg_accum${n}$(bgd, s, ${args}$)
    use mod_global_parameters

    type(block_grid_t), intent(inout) :: bgd
    double precision, intent(in)      :: s, ${coefs}$
    type(block_grid_t), intent(in)    :: ${srcs}$

    integer :: iigrid, igrid, ix1,ix2,ix3, iw

    ${GPU_PARALLEL_LOOP_GANG("private(igrid)")}$
    do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
       ${GPU_LOOP_VECTOR("collapse(ndim+1)")}$
       do iw = 1, nw
          do ix3 = ixGlo3, ixGhi3
             do ix2 = ixGlo2, ixGhi2
                do ix1 = ixGlo1, ixGhi1
                   bgd%w(ix1,ix2,ix3,iw,igrid) = s * bgd%w(ix1,ix2,ix3,iw,igrid) &
                        + ${sum}$
                end do
             end do
          end do
       end do
    end do

  end subroutine bg_accum${n}$

#:endfor

  !> Implicit global update step within IMEX schemes, advance psa=psb+dtfactor*qdt*F_im(psa)
  subroutine global_implicit_update(dtfactor,qdt,qtC,psa,psb)
    use mod_global_parameters
    use mod_ghostcells_update
    use mod_physics, only: phys_implicit_update, phys_req_diagonal

    type(state), target :: psa(max_blocks) !< Compute implicit part from this state and update it
    type(state), target :: psb(max_blocks)   !< Will be unchanged, as on entry
    double precision, intent(in) :: qdt      !< overall time step dt
    double precision, intent(in) :: qtC !< Both states psa and psb at this time level
    double precision, intent(in) :: dtfactor !< Advance psa=psb+dtfactor*qdt*F_im(psa)

    integer                        :: iigrid, igrid

    !> First copy all variables from a to b, this is necessary to account for
    ! quantities is w with no implicit sourceterm
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       psa(igrid)%w = psb(igrid)%w
    end do

    if (associated(phys_implicit_update)) then
       call phys_implicit_update(dtfactor,qdt,qtC,psa,psb)
    end if

    ! enforce boundary conditions for psa
    call getbc(qtC,0.d0,psa,iwstart,nwgc,phys_req_diagonal)

  end subroutine global_implicit_update

  !> Evaluate Implicit part in place, i.e. psa==>F_im(psa)
  subroutine evaluate_implicit(qtC,psa)
    use mod_global_parameters
    use mod_physics, only: phys_evaluate_implicit

    type(state), target :: psa(max_blocks) !< Compute implicit part from this state and update it
    double precision, intent(in) :: qtC      !< psa at this time level

    if (associated(phys_evaluate_implicit)) then
       call phys_evaluate_implicit(qtC,psa)
    end if

  end subroutine evaluate_implicit

  !> Integrate all grids by one partial step: time level kb += dtfactor*dt
  !> times the flux divergence and sources evaluated on time level ka
  subroutine advect1(dtfactor,idimmin,idimmax,qtC,ka,qt,kb)
    use mod_global_parameters
    use mod_ghostcells_update
    use mod_fix_conserve
    use mod_physics
    use mod_finite_volume, only: finite_volume_local
    use mod_usr_methods

    integer, intent(in)          :: idimmin,idimmax
    double precision, intent(in) :: dtfactor !< Advance over dtfactor * dt
    double precision, intent(in) :: qtC !< Time of level ka
    integer, intent(in)          :: ka  !< Compute fluxes based on this level
    double precision, intent(in) :: qt  !< Time of level kb before the update
    integer, intent(in)          :: kb  !< Update solution on this level

    integer :: ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision             :: qdt
    integer                      :: iigrid, igrid
    type(state), dimension(:), pointer, contiguous :: psb

    istep = istep+1

    ! AGILE doesn't happen in our test case
    ! if(associated(phys_special_advance)) then
    !  call phys_special_advance(qtC,psa)
    ! end if

    qdt=dtfactor*dt
    ! FIXME: AGILE The following replaces the `call advect1_grid` loop. The
    ! `advect1_grid` variable is a function pointer to the configured solver.
    ! Since NVidia doesn't support function pointers, we hard-code to use
    ! `finite_volume_all` here. In the future this needs to be replaced with
    ! some logic to select the desired method.

    ixOmin1=ixGlo1+nghostcells;ixOmin2=ixGlo2+nghostcells
    ixOmin3=ixGlo3+nghostcells;ixOmax1=ixGhi1-nghostcells
    ixOmax2=ixGhi2-nghostcells;ixOmax3=ixGhi3-nghostcells;
    call finite_volume_local( &
        qdt, dtfactor, & !some scalars related to time stepping
        ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,ixOmin1,ixOmin2,ixOmin3,&
           ixOmax1,ixOmax2,ixOmax3, idimmin,idimmax, & !bounds for some arrays
        qtC, &                          ! scalar related to time stepping
        bg(ka), &                       ! first block grid
        qt,  &                          ! scalar related to time stepping
        bg(kb) &                        ! second block grid
        )

    ! Off for the predictor of Predictor_Corrector, and with a single level
    if (fix_conserve_at_step) then
       call recvflux(idimmin,idimmax)
       call sendflux(idimmin,idimmax)
       call fix_conserve(kb,idimmin,idimmax,1,nwflux)
    end if

    !   if(stagger_grid) then
    !     call fix_edges(psb,idimmin,idimmax)
    !     ! fill the cell-center values from the updated staggered variables
    !     do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
    !       call phys_face_to_center(ixMlo1,ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3,&
    !          psb(igrid))
    !     end do
    !   end if
    ! end if

    
    ! For all grids: fill ghost cells; getbc finds level kb from the views
    psb => ps_level(kb)
    call getbc(qt+qdt,qdt,psb,iwstart,nwgc,phys_req_diagonal)

  end subroutine advect1

  !> process is a user entry in time loop, before output and advance
  !>         allows to modify solution, add extra variables, etc.
  !> Warning: CFL dt already determined (and is not recomputed)!
  subroutine process(iit,qt)
    use mod_usr_methods, only: usr_process_grid, usr_process_global
    use mod_global_parameters
    use mod_geometry, only: sync_positions_host
    use mod_ghostcells_update
    use mod_physics, only: phys_req_diagonal
    ! .. scalars ..
    integer,intent(in)          :: iit
    double precision, intent(in):: qt

    integer:: iigrid, igrid

    if (associated(usr_process_global)) then
       call usr_process_global(iit,qt)
    end if

    if (associated(usr_process_grid)) then
      ! usr_process_grid reads ps(igrid)%x on the host
      call sync_positions_host()
      do iigrid=1,igridstail; igrid=igrids(iigrid);
         ! next few lines ensure correct usage of routines like divvector etc
         dxlevel(1)=rnode(rpdx1_,igrid);dxlevel(2)=rnode(rpdx2_,igrid)
         dxlevel(3)=rnode(rpdx3_,igrid);
         block=>ps(igrid)
         call usr_process_grid(igrid,node(plevel_,igrid),ixGlo1,ixGlo2,ixGlo3,&
            ixGhi1,ixGhi2,ixGhi3,ixMlo1,ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3, qt,&
            ps(igrid)%w,ps(igrid)%x)
      end do
      call getbc(qt,dt,ps,iwstart,nwgc,phys_req_diagonal)
    end if
  end subroutine process

  !> process_advanced is user entry in time loop, just after advance
  !>           allows to modify solution, add extra variables, etc.
  !>           added for handling two-way coupled PIC-MHD
  !> Warning: w is now at global_time^(n+1), global time and iteration at global_time^n, it^n
  subroutine process_advanced(iit,qt)
    use mod_usr_methods, only: usr_process_adv_grid, usr_process_adv_global
    use mod_global_parameters
    use mod_geometry, only: sync_positions_host
    use mod_ghostcells_update
    use mod_physics, only: phys_req_diagonal
    ! .. scalars ..
    integer,intent(in)          :: iit
    double precision, intent(in):: qt

    integer:: iigrid, igrid

    if (associated(usr_process_adv_global)) then
       call usr_process_adv_global(iit,qt)
    end if

    if (associated(usr_process_adv_grid)) then
      ! usr_process_adv_grid reads ps(igrid)%x on the host
      call sync_positions_host()
      do iigrid=1,igridstail; igrid=igrids(iigrid);
         ! next few lines ensure correct usage of routines like divvector etc
         dxlevel(1)=rnode(rpdx1_,igrid);dxlevel(2)=rnode(rpdx2_,igrid)
         dxlevel(3)=rnode(rpdx3_,igrid);
         block=>ps(igrid)

         call usr_process_adv_grid(igrid,node(plevel_,igrid),ixGlo1,ixGlo2,&
            ixGlo3,ixGhi1,ixGhi2,ixGhi3,ixMlo1,ixMlo2,ixMlo3,ixMhi1,ixMhi2,&
            ixMhi3, qt,ps(igrid)%w,ps(igrid)%x)
      end do
      call getbc(qt,dt,ps,iwstart,nwgc,phys_req_diagonal)
    end if
  end subroutine process_advanced

end module mod_advance
