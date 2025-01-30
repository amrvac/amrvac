!> Module containing all the time stepping schemes
module mod_advance

  implicit none
  private

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

    ! !$acc update device(ps(1:max_blocks))
    ! do iigrid=1,igridstail; igrid=igrids(iigrid);
    !    !$acc enter data copyin(ps(igrid)%w, ps(igrid)%x) create(ps1(igrid)%w, ps2(igrid)%w)
    ! end do
    
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

    if(use_particles) call handle_particles

    ! do iigrid=1,igridstail; igrid=igrids(iigrid);
    !    !$acc exit data delete(ps(igrid)%x, ps1(igrid)%w, ps2(igrid)%w) copyout(ps(igrid)%w)
    ! end do
    
  end subroutine advance

  !> Advance all grids over one time step, but without taking dimensional
  !> splitting or split source terms into account
  subroutine advect(idim^LIM)
    use mod_global_parameters
    use mod_fix_conserve
    use mod_ghostcells_update
    use mod_physics, only: phys_req_diagonal
    use mod_comm_lib, only: mpistop
    use mod_physicaldata

    integer, intent(in) :: idim^LIM
    integer             :: iigrid, igrid, ix^D, iw

    call init_comm_fix_conserve(idim^LIM,nwflux)
    fix_conserve_at_step = time_advance .and. levmax>levmin
    
    ! OpenACC data region to manage data movement
!    !$acc update device(bg(1)%w, bg(2)%w)
!    print*,'bg start', bg(1)%w(28,41,3,1), '  *** ', bg(2)%w(28,41,3,1), ' *** ', bg(3)%w(28,41,3,1)

    ! copy w instead of wold because of potential use of dimsplit or sourcesplit
    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       !$acc parallel loop collapse(ndim+1) present(bg)
       do iw = 1, nw
          {^D& do ix^DB = ixGlo^DB, ixGhi^DB \}
          bg(2)%w(ix^D,iw,igrid) = bg(1)%w(ix^D,iw,igrid)
          {^D& end do \}
       end do
          if(stagger_grid) then
             !$acc kernels
             ps1(igrid)%ws=ps(igrid)%ws
             !$acc end kernels
          end if
    end do
    !$OMP END PARALLEL DO

!    !$acc update self(bg(1)%w, bg(2)%w)
    istep = 0

     select case (t_stepper)
  
    case (threestep)
       select case (t_integrator)
       ! AGILE this is our integrator (default threestep)
       case (ssprk3)
          ! this is SSPRK(3,3) Gottlieb-Shu 1998 or SSP(3,2) depending on ssprk_order (3 vs 2)
         
          ! TODO call advect1 with bg(2) instead of ps1 ??? 
          call advect1(flux_method,rk_beta11, idim^LIM,global_time,ps,bg(1),global_time,ps1,bg(2))
!         !$acc update device(bg(1)%w, bg(2)%w, bg(3)%w)
          !$OMP PARALLEL DO PRIVATE(igrid)
           do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
              !$acc parallel loop collapse(ndim+1) present(bg)
              do iw = 1, nw
                 {^D& do ix^DB = ixGlo^DB, ixGhi^DB \}
                 bg(3)%w(ix^D,iw,igrid) = rk_alfa21 * bg(1)%w(ix^D,iw,igrid) + rk_alfa22 * bg(2)%w(ix^D,iw,igrid)
                 {^D& end do \}
              end do
              if(stagger_grid) ps2(igrid)%ws=rk_alfa21*ps(igrid)%ws+rk_alfa22*ps1(igrid)%ws
           end do
          !$OMP END PARALLEL DO
!          !$acc update self(bg(1)%w, bg(2)%w, bg(3)%w)

          ! TODO call advect1 with bg(3) instead of ps2 ??? 
          call advect1(flux_method,rk_beta22, idim^LIM,global_time+rk_c2*dt,ps1,bg(2),global_time+rk_alfa22*rk_c2*dt,ps2,bg(3))

!          !$acc update device(bg(1)%w, bg(2)%w, bg(3)%w)
          !$OMP PARALLEL DO PRIVATE(igrid)
           do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
              !$acc parallel loop collapse(ndim+1) present(bg)
              do iw = 1, nw
                 {^D& do ix^DB = ixGlo^DB, ixGhi^DB \}
                 bg(1)%w(ix^D,iw,igrid) = rk_alfa31 * bg(1)%w(ix^D,iw,igrid) + rk_alfa33 * bg(3)%w(ix^D,iw,igrid)
                 {^D& end do \}
              end do
              if(stagger_grid) ps(igrid)%ws=rk_alfa31*ps(igrid)%ws+rk_alfa33*ps2(igrid)%ws
           end do
          !$OMP END PARALLEL DO
!          !$acc update self(bg(1)%w, bg(2)%w, bg(3)%w)

          ! TODO call advect1 with bg(1) instead of ps ??? 
          call advect1(flux_method,rk_beta33, &
                idim^LIM,global_time+rk_c3*dt,ps2,bg(3),global_time+(1.0d0-rk_beta33)*dt,ps,bg(1))
!          !$acc update host(bg(1)%w, bg(2)%w, bg(3)%w)
!          print*, 'advect, 6. called afvect1 with bg(3) and bg(1)', bg(1)%w(28,41,3,1), bg(2)%w(28,41,3,1), bg(3)%w(28,41,3,1)
  
        case default
           call mpistop("unkown threestep time_integrator in advect")
        end select

    case default
       call mpistop("unkown time_stepper in advect")
    end select

  end subroutine advect

  !> Implicit global update step within IMEX schemes, advance psa=psb+dtfactor*qdt*F_im(psa)
  subroutine global_implicit_update(dtfactor,qdt,qtC,psa,psb)
    use mod_global_parameters
    use mod_ghostcells_update
    use mod_physics, only: phys_implicit_update, phys_req_diagonal

    type(state), target :: psa(max_blocks)   !< Compute implicit part from this state and update it
    type(state), target :: psb(max_blocks)   !< Will be unchanged, as on entry
    double precision, intent(in) :: qdt      !< overall time step dt
    double precision, intent(in) :: qtC      !< Both states psa and psb at this time level
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

    type(state), target :: psa(max_blocks)   !< Compute implicit part from this state and update it
    double precision, intent(in) :: qtC      !< psa at this time level

    if (associated(phys_evaluate_implicit)) then
       call phys_evaluate_implicit(qtC,psa)
    end if

  end subroutine evaluate_implicit

  !> Integrate all grids by one partial step
  subroutine advect1(method,dtfactor,idim^LIM,qtC,psa,bga,qt,psb,bgb)
    use mod_global_parameters
    use mod_ghostcells_update
    use mod_fix_conserve
    use mod_physics
    use mod_finite_volume_all, only: finite_volume_all, finite_volume_local

    integer, intent(in) :: idim^LIM
    integer :: ixO^L
    type(state), target :: psa(max_blocks) !< Compute fluxes based on this state
    type(state), target :: psb(max_blocks) !< Update solution on this state
    type(block_grid_t), target :: bga !< Compute fluxes based on this state
    type(block_grid_t), target :: bgb !< Update solution on this state
    double precision, intent(in) :: dtfactor !< Advance over dtfactor * dt
    double precision, intent(in) :: qtC
    double precision, intent(in) :: qt
    integer, intent(in) :: method(nlevelshi)

    ! cell face flux
    double precision :: fC(ixG^T,1:nwflux,1:ndim)
    ! cell edge flux
    double precision :: fE(ixG^T,sdim:3)
    !$acc declare create(fC,fE)
    double precision :: qdt
    integer :: iigrid, igrid

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

    ixO^L=ixG^LL^LSUBnghostcells;
    call finite_volume_local( &
        fs_hll, &          ! fs_hll
        qdt, dtfactor, &                ! some scalars related to time stepping
        ixG^LL,ixO^L, idim^LIM, &      ! bounds for some arrays
        qtC, &                          ! scalar related to time stepping
        bga, &                          ! first block grid
        qt,  &                          ! scalar related to time stepping
        bgb, &                          ! second block grid
        fC, fE &                        ! fluxes
    )
    !print*,'advect1, call fva', bg(1)%w(28,41,3,1)
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
!    do iigrid=1,igridstail; igrid=igrids(iigrid);
!       !$acc update self(psb(igrid)%w)
!    end do
    call getbc(qt+qdt,qdt,psb,iwstart,nwgc,phys_req_diagonal)
!    do iigrid=1,igridstail; igrid=igrids(iigrid);
!       !$acc update device(psb(igrid)%w)
!    end do

  end subroutine advect1

  !> Advance a single grid over one partial time step
  subroutine advect1_grid(method,qdt,dtfactor,ixI^L,idim^LIM,qtC,sCT,qt,s,fC,fE,dxs,x)

    !  integrate one grid by one partial step
    use mod_finite_volume_all
    use mod_finite_volume
    use mod_finite_difference
    use mod_tvd
    use mod_source, only: addsource2
    use mod_physics, only: phys_to_primitive
    use mod_global_parameters
    use mod_comm_lib, only: mpistop

    integer, intent(in) :: method
    integer, intent(in) :: ixI^L, idim^LIM
    double precision, intent(in) :: qdt, dtfactor, qtC, qt, dxs(ndim), x(ixI^S,1:ndim)
    type(state), target          :: sCT, s
    double precision :: fC(ixI^S,1:nwflux,1:ndim), wprim(ixI^S,1:nw)
    double precision :: fE(ixI^S,sdim:3)

    integer :: ixO^L

    ixO^L=ixI^L^LSUBnghostcells;
    
    select case (method)
    case (fs_hll,fs_hllc,fs_hllcd,fs_hlld,fs_tvdlf,fs_tvdmu)
       call finite_volume(method,qdt,dtfactor,ixI^L,ixO^L,idim^LIM,qtC,sCT,qt,s,fC,fE,dxs,x)
    case (fs_cd,fs_cd4)
       call centdiff(method,qdt,dtfactor,ixI^L,ixO^L,idim^LIM,qtC,sCT,qt,s,fC,fE,dxs,x)
    !case (fs_hancock)
    !   call hancock(qdt,dtfactor,ixI^L,ixO^L,idim^LIM,qtC,sCT,qt,s,dxs,x)
    case (fs_fd)
       call fd(qdt,dtfactor,ixI^L,ixO^L,idim^LIM,qtC,sCT,qt,s,fC,fE,dxs,x)
    case (fs_tvd)
       call centdiff(fs_cd,qdt,dtfactor,ixI^L,ixO^L,idim^LIM,qtC,sCT,qt,s,fC,fE,dxs,x)
       call tvdlimit(method,qdt,ixI^L,ixO^L,idim^LIM,sCT,qt+qdt,s,fC,dxs,x)
    case (fs_source)
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

  end subroutine advect1_grid

  !> process is a user entry in time loop, before output and advance
  !>         allows to modify solution, add extra variables, etc.
  !> Warning: CFL dt already determined (and is not recomputed)!
  subroutine process(iit,qt)
    use mod_usr_methods, only: usr_process_grid, usr_process_global
    use mod_global_parameters
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
      !$OMP PARALLEL DO PRIVATE(igrid)
      do iigrid=1,igridstail; igrid=igrids(iigrid);
         ! next few lines ensure correct usage of routines like divvector etc
         ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
         block=>ps(igrid)
         call usr_process_grid(igrid,node(plevel_,igrid),ixG^LL,ixM^LL, &
              qt,ps(igrid)%w,ps(igrid)%x)
      end do
      !$OMP END PARALLEL DO
      call getbc(qt,dt,ps,iwstart,nwgc,phys_req_diagonal)
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
    use mod_physics, only: phys_req_diagonal
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
      call getbc(qt,dt,ps,iwstart,nwgc,phys_req_diagonal)
    end if
  end subroutine process_advanced

end module mod_advance
