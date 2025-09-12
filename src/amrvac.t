!> AMRVAC solves a set of hyperbolic equations
!> \f$\vec{u}_t + \nabla_x \cdot \vec{f}(\vec{u}) = \vec{s}\f$
!> using adaptive mesh refinement.
program amrvac

  use mod_global_parameters
  use mod_input_output
  use mod_usr_methods
  use mod_ghostcells_update
  use mod_usr
  use mod_initialize
  use mod_initialize_amr, only: initlevelone, modify_IC
  use mod_initialize_amr, only: improve_initial_condition
  use mod_selectgrids, only: selectgrids
  use mod_particles
  use mod_fix_conserve
  use mod_advance, only: process
  use mod_multigrid_coupling
  use mod_convert, only: init_convert
  use mod_physics
  use mod_amr_grid, only: resettree, settree, resettree_convert
  use mod_trac, only: initialize_trac_after_settree
  use mod_convert_files, only: generate_plotfile
  use mod_comm_lib, only: comm_start, comm_finalize,mpistop


  double precision :: time0, time_in
  logical,save     :: part_file_exists=.false.


  call comm_start()

  time0 = MPI_WTIME()
  time_advance = .false.
  time_bc      = zero

  ! read command line arguments first
  call read_arguments()

  ! init_convert is called before usr_init as user might associate a convert method
  call init_convert()
  ! the user_init routine should load a physics module
  call usr_init()

  call initialize_amrvac()

  if (restart_from_file /= undefined) then
     ! restart from previous file or dat file conversion
     ! get input data from previous AMRVAC run

     ! read in dat file
     call read_snapshot()

     ! rewrite it=0 snapshot when restart from it=0 state 
     if(it==0.and.itsave(1,2)==0) snapshotnext=snapshotnext-1

     if (reset_time) then
       ! reset it and global time to original value
       it           = it_init
       global_time  = time_init
       ! reset snapshot number
       snapshotnext=0
     end if

     if (reset_it) then
       ! reset it to original value
       it           = it_init
     end if

     ! allow use to read extra data before filling boundary condition
     if (associated(usr_process_grid) .or. &
          associated(usr_process_global)) then
        call process(it,global_time)
     end if

     ! modify initial condition
     if (firstprocess) then
       ! update ghost cells for all need-boundary variables before modification
       call getbc(global_time,0.d0,ps,1,nwflux+nwaux)
       call modify_IC
     end if

     ! select active grids
     call selectgrids

     ! update ghost cells for all need-boundary variables
     call getbc(global_time,0.d0,ps,1,nwflux+nwaux)

     ! reset AMR grid
     if (reset_grid) then
       call settree
     else
       ! set up boundary flux conservation arrays
       if (levmax>levmin) call allocateBflux
     end if

     ! all blocks refined to the same level for output
     if(convert .and. level_io>0 .or. level_io_min.ne.1 .or. level_io_max.ne.nlevelshi) &
       call resettree_convert

     {^NOONED
     ! improve initial condition after restart and modification
     if(firstprocess) call improve_initial_condition()
     }

     if (use_multigrid) call mg_setup_multigrid()

     if(use_particles) then
       call read_particles_snapshot(part_file_exists)
       call init_gridvars()
       if (.not. part_file_exists) call particles_create()
       if(convert) then
         call handle_particles()
         call finish_gridvars()
         call time_spent_on_particles()
         call comm_finalize
         stop
       end if
     end if

     if(convert) then
       if (npe/=1.and.(.not.(index(convert_type,'mpi')>=1)) &
            .and. convert_type .ne. 'user')  &
            call mpistop("non-mpi conversion only uses 1 cpu")
       if(mype==0.and.level_io>0) write(unitterm,*)'reset tree to fixed level=',level_io

       !here requires -1 snapshot
       if (autoconvert .or. snapshotnext>0) snapshotnext = snapshotnext - 1

       if(associated(phys_special_advance)) then
         ! e.g. calculate MF velocity from magnetic field
         call phys_special_advance(global_time,ps)
       end if

       call generate_plotfile
       call comm_finalize
       stop
     end if

  else

     ! form and initialize all grids at level one
     call initlevelone

     ! set up and initialize finer level grids, if needed
     call settree

     if (use_multigrid) call mg_setup_multigrid()

     ! improve initial condition
     call improve_initial_condition()

     ! select active grids
     call selectgrids

     if (use_particles) then
       call init_gridvars()
       call particles_create()
     end if

  end if

  ! initialize something base on tree information
  call initialize_trac_after_settree

  if (mype==0) then
     print*,'-------------------------------------------------------------------------------'
     write(*,'(a,f17.3,a)')' Startup phase took : ',MPI_WTIME()-time0,' sec'
     print*,'-------------------------------------------------------------------------------'
  end if

  ! an interface to allow user to do special things before the main loop
  if (associated(usr_before_main_loop)) &
       call usr_before_main_loop()

  ! do time integration of all grids on all levels
  call timeintegration()

  if (mype==0) then
     print*,'-------------------------------------------------------------------------------'
     write(*,'(a,f17.3,a)')' Finished AMRVAC in : ',MPI_WTIME()-time0,' sec'
     print*,'-------------------------------------------------------------------------------'
  end if

  call comm_finalize

contains

  subroutine timeintegration()
    use mod_timing
    use mod_advance, only: advance, process, process_advanced
    use mod_forest, only: nleafs_active
    use mod_global_parameters
    use mod_input_output, only: saveamrfile
    use mod_input_output_helper, only: save_now
    use mod_ghostcells_update
    use mod_dt, only: setdt


    double precision :: time_last_print, time_write0, time_write, time_before_advance, dt_loop
    integer(kind=8) ncells_update
    integer :: level, ifile, fixcount, ncells_block, igrid, iigrid
    logical :: crashall

    time_in=MPI_WTIME()
    time_last_print = -bigdouble
    fixcount=1

    n_saves(filelog_:fileout_) = snapshotini

    do ifile=nfile,1,-1
       if(resume_previous_run) then
         tsavelast(ifile)=aint((global_time+smalldouble)/dtsave(ifile))*dtsave(ifile)
         itsavelast(ifile)=it/ditsave(ifile)*ditsave(ifile)
       else
         tsavelast(ifile)=global_time
         itsavelast(ifile)=it
       end if
    end do

    ! the next two are used to keep track of the performance during runtime:
    itTimeLast=it
    timeLast=MPI_WTIME()

    !  ------ start of integration loop. ------------------
    if (mype==0) then
      write(*, '(A,ES9.2,A)') ' Start integrating, print status every ', &
           time_between_print, ' seconds'
      write(*, '(A4,A10,A12,A12,A12)') '  #', 'it', 'time', 'dt', 'wc-time(s)'
    end if

    timeloop0=MPI_WTIME()
    time_bc=0.d0
    time_write=0.d0
    ncells_block={(ixGhi^D-2*nghostcells)*}
    ncells_update=0
    dt_loop=0.d0

    time_advance=.true.

    time_evol : do

       time_before_advance=MPI_WTIME()
       ! set time step
       call setdt()

       ! Optionally call a user method that can modify the grid variables at the
       ! beginning of a time step
       if (associated(usr_process_grid) .or. &
            associated(usr_process_global)) then
          call process(it,global_time)
       end if

       ! Check if output needs to be written
       do ifile=nfile,1,-1
         save_file(ifile) = timetosave(ifile)
       end do

       timeio0=MPI_WTIME()

       if (timeio0 - time_last_print > time_between_print) then
         time_last_print = timeio0
         if (mype == 0) then
           write(*, '(A4,I10,ES12.4,ES12.4,ES12.4)') " #", &
                it, global_time, dt, timeio0 - time_in
         end if
       end if

       ! output data
       if (any(save_file)) then
         if(associated(usr_modify_output)) then
           ! Users can modify or set variables before output is written
           do iigrid=1,igridstail; igrid=igrids(iigrid);
             ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
             block=>ps(igrid)
             call usr_modify_output(ixG^LL,ixM^LL,global_time,ps(igrid)%w,ps(igrid)%x)
           end do
         end if
         time_write=0.d0
         do ifile=nfile,1,-1
           if (save_file(ifile)) then
             time_write0=MPI_WTIME()
             call saveamrfile(ifile)
             time_write=time_write+MPI_WTIME()-time_write0
           end if
         end do
       end if

       ! output a snapshot when user write a file named 'savenow' in the same
       ! folder as the executable amrvac
       if (mype==0) inquire(file='savenow',exist=save_now)
       if (npe>1) call MPI_BCAST(save_now,1,MPI_LOGICAL,0,icomm,ierrmpi)

       if (save_now) then
          if(mype==0) write(*,'(a,i7,a,i7,a,es12.4)') ' save a snapshot No.',&
               snapshotnext,' at it=',it,' global_time=',global_time
          call saveamrfile(1)
          call saveamrfile(2)
          call MPI_FILE_DELETE('savenow',MPI_INFO_NULL,ierrmpi)
       end if
       timeio_tot=timeio_tot+MPI_WTIME()-timeio0

       pass_wall_time=MPI_WTIME()-time0+dt_loop+4.d0*time_write >=wall_time_max

       ! exit time loop if time is up
       if (it>=it_max .or. global_time>=time_max .or. pass_wall_time .or. final_dt_exit) exit time_evol

       ! solving equations
       call advance(it)

       ! if met unphysical values, output the last good status and stop the run
       call MPI_ALLREDUCE(crash,crashall,1,MPI_LOGICAL,MPI_LOR,icomm,ierrmpi)
       if (crashall) then
         call saveamrfile(1)
         call saveamrfile(2)
         if(mype==0) write(*,*) "Error: small value encountered, run crash."
         call MPI_ABORT(icomm, iigrid, ierrmpi)
       end if

       ! Optionally call a user method that can modify the grid variables at the
       ! end of a time step: this is for two-way coupling to PIC, e.g.
       if (associated(usr_process_adv_grid) .or. &
            associated(usr_process_adv_global)) then
          call process_advanced(it,global_time)
       end if

       ! update AMR mesh and tree
       timegr0=MPI_WTIME()
       if(ditregrid>1) then
          if(fixcount<ditregrid) then
             fixcount=fixcount+1
          else
             if (refine_max_level>1 .and. .not.(fixgrid())) call resettree
             fixcount=1
          end if
       else
          if (refine_max_level>1 .and. .not.(fixgrid())) call resettree
       end if
       timegr_tot=timegr_tot+(MPI_WTIME()-timegr0)

       ! update time variables
       it = it + 1
       global_time = global_time + dt

       if(it>9000000)then
          it = slowsteps+it_init
          itsavelast(:)=0
       end if

       ! count updated cells
       ncells_update=ncells_update+ncells_block*nleafs_active

       ! time lapses in one loop
       dt_loop=MPI_WTIME()-time_before_advance
    end do time_evol

    if(use_particles) then
      call finish_gridvars()
    end if

    time_advance=.false.

    timeloop=MPI_WTIME()-timeloop0

    if (mype==0) then
       write(*,'(a,f12.3,a)')' Total timeloop took        : ',timeloop,' sec'
       write(*,'(a,f12.3,a)')' Time spent on AMR          : ',timegr_tot,' sec'
       write(*,'(a,f12.2,a)')'                  Percentage: ',100.0*timegr_tot/timeloop,' %'
       write(*,'(a,f12.3,a)')' Time spent on IO in loop   : ',timeio_tot,' sec'
       write(*,'(a,f12.2,a)')'                  Percentage: ',100.0*timeio_tot/timeloop,' %'
       write(*,'(a,f12.3,a)')' Time spent on ghost cells  : ',time_bc,' sec'
       write(*,'(a,f12.2,a)')'                  Percentage: ',100.0*time_bc/timeloop,' %'
       write(*,'(a,f12.3,a)')' Time spent on computing    : ',timeloop-timeio_tot-timegr_tot-time_bc,' sec'
       write(*,'(a,f12.2,a)')'                  Percentage: ',100.0*(timeloop-timeio_tot-timegr_tot-time_bc)/timeloop,' %'
       write(*,'(a,es12.3 )')' Cells updated / proc / sec : ',dble(ncells_update)*dble(nstep)/dble(npe)/timeloop
    end if

    ! output end state
    timeio0=MPI_WTIME()
    do ifile=nfile,1,-1
       if(itsavelast(ifile)<it) call saveamrfile(ifile)
    end do
    if (mype==0) call MPI_FILE_CLOSE(log_fh,ierrmpi)
    timeio_tot=timeio_tot+(MPI_WTIME()-timeio0)

    if (mype==0) then
       write(*,'(a,f12.3,a)')' Total time spent on IO     : ',timeio_tot,' sec'
       write(*,'(a,f12.3,a)')' Total timeintegration took : ',MPI_WTIME()-time_in,' sec'
       write(*, '(A4,I10,ES12.3,ES12.3,ES12.3)') " #", &
            it, global_time, dt, timeio0 - time_in
    end if

    {#IFDEF RAY
    call time_spent_on_rays
    }

    if(use_particles) call time_spent_on_particles

    if (use_multigrid) call mg_timers_show(mg)
  end subroutine timeintegration

  !> Save times are defined by either tsave(isavet(ifile),ifile) or
  !> itsave(isaveit(ifile),ifile) or dtsave(ifile) or ditsave(ifile)
  !> tsavestart(ifile) determines first start time. This only affects
  !> read out times determined by dtsave(ifiles).
  !> Other conditions may be included.
  logical function timetosave(ifile)
    use mod_global_parameters

    integer:: ifile
    logical:: oksave

    oksave=.false.
    if (it==itsave(isaveit(ifile),ifile)) then
       oksave=.true.
       isaveit(ifile)=isaveit(ifile)+1
    end if
    if (it==itsavelast(ifile)+ditsave(ifile)) oksave=.true.

    if (global_time>=tsave(isavet(ifile),ifile).and.global_time-dt<tsave(isavet(ifile),ifile)) then
       oksave=.true.
       isavet(ifile)=isavet(ifile)+1
    end if

    if(global_time>=tsavestart(ifile)-smalldouble)then
      if (global_time>=tsavelast(ifile)+dtsave(ifile)-smalldouble)then
         oksave=.true.
         n_saves(ifile) = n_saves(ifile) + 1
      endif
    endif

    if (oksave) then
       tsavelast(ifile) =global_time
       itsavelast(ifile)=it
    end if
    timetosave=oksave

    return
  end function timetosave

  !> Return true if the AMR grid should not be adapted any more. This is
  !> controlled by tfixgrid or itfixgrid. Other conditions may be included.
  logical function fixgrid()
    use mod_global_parameters

    fixgrid= (global_time>=tfixgrid .or. it>=itfixgrid)
  end function fixgrid

end program amrvac
