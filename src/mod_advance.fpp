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
    !$OMP PARALLEL DO PRIVATE(igrid)
    !$acc parallel loop present(bg, bg(1), bg(2)) private(igrid)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       !$acc loop collapse(ndim+1)
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
    !$OMP END PARALLEL DO
    
    if(stagger_grid) then
       !$OMP PARALLEL DO PRIVATE(igrid)
       !$acc parallel loop present(ps1, ps) private(igrid)
       do iigrid=1,igridstail; igrid=igrids(iigrid);
          !$acc loop collapse(ndim+1)
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
    !$OMP END PARALLEL DO

 istep = 0

 select case (t_stepper)
  
    case (threestep)
       select case (t_integrator)
       ! AGILE this is our integrator (default threestep)
       case (ssprk3)
          ! this is SSPRK(3,3) Gottlieb-Shu 1998 or SSP(3,2) depending on ssprk_order (3 vs 2)
         
          call advect1(flux_method,rk_beta11, idimmin,idimmax,global_time,ps,&
             bg(1),global_time,ps1,bg(2))

          !$OMP PARALLEL DO PRIVATE(igrid)
          !$acc parallel loop present(bg, ps2, ps1, ps) private(igrid)
           do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
              !$acc loop collapse(ndim+1)
              do iw = 1, nw
                  do ix3 = ixGlo3, ixGhi3 
                   do ix2 = ixGlo2, ixGhi2 
                   do ix1 = ixGlo1, ixGhi1 
                 bg(3)%w(ix1,ix2,ix3,iw,igrid) = rk_alfa21 * bg(1)%w(ix1,ix2,&
                    ix3,iw,igrid) + rk_alfa22 * bg(2)%w(ix1,ix2,ix3,iw,igrid)
                  end do 
                   end do 
                   end do 
              end do
              if(stagger_grid) ps2(igrid)%ws=rk_alfa21*ps(igrid)%ws+&
                 rk_alfa22*ps1(igrid)%ws
           end do
          !$OMP END PARALLEL DO

          call advect1(flux_method,rk_beta22, idimmin,idimmax,&
             global_time+rk_c2*dt,ps1,bg(2),global_time+rk_alfa22*rk_c2*dt,ps2,&
             bg(3))

          !$OMP PARALLEL DO PRIVATE(igrid)
          !$acc parallel loop present(bg, ps2, ps) private(igrid)
           do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
              !$acc loop collapse(ndim+1)
              do iw = 1, nw
                  do ix3 = ixGlo3, ixGhi3 
                   do ix2 = ixGlo2, ixGhi2 
                   do ix1 = ixGlo1, ixGhi1 
                 bg(1)%w(ix1,ix2,ix3,iw,igrid) = rk_alfa31 * bg(1)%w(ix1,ix2,&
                    ix3,iw,igrid) + rk_alfa33 * bg(3)%w(ix1,ix2,ix3,iw,igrid)
                  end do 
                   end do 
                   end do 
              end do
              if(stagger_grid) ps(igrid)%ws=rk_alfa31*ps(igrid)%ws+&
                 rk_alfa33*ps2(igrid)%ws
           end do
          !$OMP END PARALLEL DO

          call advect1(flux_method,rk_beta33, idimmin,idimmax,&
             global_time+rk_c3*dt,ps2,bg(3),global_time+(1.0d0-rk_beta33)*dt,&
             ps,bg(1))
  
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

  !> Integrate all grids by one partial step
  subroutine advect1(method,dtfactor,idimmin,idimmax,qtC,psa,bga,qt,psb,bgb)
    use mod_global_parameters
    use mod_ghostcells_update
    use mod_fix_conserve
    use mod_physics
    use mod_finite_volume, only: finite_volume_local

    integer, intent(in)          :: idimmin,idimmax
    integer :: ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    type(state), target          :: psa(max_blocks) !< Compute fluxes based on this state
    type(state), target          :: psb(max_blocks) !< Update solution on this state
    type(block_grid_t)           :: bga !< Compute fluxes based on this state
    type(block_grid_t)           :: bgb !< Update solution on this state
    double precision, intent(in) :: dtfactor !< Advance over dtfactor * dt
    double precision, intent(in) :: qtC
    double precision, intent(in) :: qt
    integer, intent(in)          :: method(nlevelshi)

    ! cell face flux
    double precision             :: fC(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3,1:nwflux,1:ndim)
    ! cell edge flux
    double precision             :: fE(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3,sdim:3)
    !$acc declare create(fC,fE)
    double precision             :: qdt
    integer                      :: iigrid, igrid

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
    call finite_volume_local( fs_hll, &          ! fs_hll
        qdt, dtfactor, & !some scalars related to time stepping
        ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,ixOmin1,ixOmin2,ixOmin3,&
           ixOmax1,ixOmax2,ixOmax3, idimmin,idimmax, & !bounds for some arrays
        qtC, &                          ! scalar related to time stepping
        bga, &                          ! first block grid
        qt,  &                          ! scalar related to time stepping
        bgb, &                          ! second block grid
        fC, fE &                        ! fluxes
        )
    
    if (fix_conserve_global .and. fix_conserve_at_step) then
      call recvflux(idimmin,idimmax)
      call sendflux(idimmin,idimmax)
      call fix_conserve(psb,idimmin,idimmax,1,nwflux)
      if(stagger_grid) then
        call fix_edges(psb,idimmin,idimmax)
        ! fill the cell-center values from the updated staggered variables
        !$OMP PARALLEL DO PRIVATE(igrid)
        do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
          call phys_face_to_center(ixMlo1,ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3,&
             psb(igrid))
        end do
        !$OMP END PARALLEL DO
      end if
    end if

    ! For all grids: fill ghost cells
    call getbc(qt+qdt,qdt,psb,iwstart,nwgc,phys_req_diagonal)

  end subroutine advect1

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
         dxlevel(1)=rnode(rpdx1_,igrid);dxlevel(2)=rnode(rpdx2_,igrid)
         dxlevel(3)=rnode(rpdx3_,igrid);
         block=>ps(igrid)
         call usr_process_grid(igrid,node(plevel_,igrid),ixGlo1,ixGlo2,ixGlo3,&
            ixGhi1,ixGhi2,ixGhi3,ixMlo1,ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3, qt,&
            ps(igrid)%w,ps(igrid)%x)
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
    use mod_usr_methods, only: usr_process_adv_grid, usr_process_adv_global
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
         dxlevel(1)=rnode(rpdx1_,igrid);dxlevel(2)=rnode(rpdx2_,igrid)
         dxlevel(3)=rnode(rpdx3_,igrid);
         block=>ps(igrid)

         call usr_process_adv_grid(igrid,node(plevel_,igrid),ixGlo1,ixGlo2,&
            ixGlo3,ixGhi1,ixGhi2,ixGhi3,ixMlo1,ixMlo2,ixMlo3,ixMhi1,ixMhi2,&
            ixMhi3, qt,ps(igrid)%w,ps(igrid)%x)
      end do
      !$OMP END PARALLEL DO
      call getbc(qt,dt,ps,iwstart,nwgc,phys_req_diagonal)
    end if
  end subroutine process_advanced

end module mod_advance
