!> AMRVAC solves a set of hyperbolic equations
!> \f$\vec{u}_t + \nabla_x \cdot \vec{f}(\vec{u}) = \vec{s}\f$
!> using adaptive mesh refinement.
program amrvac

  ! following line may avoid mystery problems in the combination
  !  where we use ifort compiler and MPT (SGIs MPI implementation)
  !DEC$ ATTRIBUTES NOINLINE :: read_snapshot

  {#IFDEF PARTICLES
  use mod_gridvars, only: init_gridvars, finish_gridvars
  }
  use mod_global_parameters
  use mod_input_output
  use mod_physics
  use mod_usr_methods, only: usr_before_main_loop
  use mod_ghostcells_update
  use mod_usr

  integer          :: itin
  double precision :: time0, time_in, tin

  call comm_start()

  time_advance = .false.
  time0        = MPI_WTIME()
  time_bc      = zero

  call read_arguments()

  call usr_init()
  call phys_check()

  call read_par_files()

  call initialize_vars()
  call init_comm_types()

  ! Begin of the code
  ! -----------------
  if (snapshotini/=-1) then
     ! restart from previous file or dat file conversion
     ! get input data from previous VAC/AMRVAC run
     itin=it
     tin=global_time

     ! order: first physics dependent part then user for defaults.
     call phys_check_params()

     {#IFDEF RAY
     call init_rays
     }

     if (typeparIO==1)then
        call read_snapshot
     else
        call read_snapshotnopar
     end if

     {#IFDEF PARTICLES
     call init_tracerparticles
     call getbc(global_time,0.d0,pw,0,nwflux+nwaux)
     call init_gridvars
     call read_particles_snapshot
     call finish_gridvars
     }

     if (convert) then
        if (npe/=1.and.(.not.(index(convert_type,'mpi')>=1)) &
             .and. convert_type .ne. 'user')  &
             call mpistop("non-mpi conversion only uses 1 cpu")
        call generate_plotfile
        {#IFDEF PARTICLES
        call finish_tracerparticles
        }
        call comm_finalize
        stop
     end if

     if (restart_reset_time) then
       ! Reset it and global_time to original values, and start writing
       ! snapshots at index 0
       it           = itin
       global_time  = tin
       snapshotnext = 0
     end if

     ! modify initial condition
     if (firstprocess) call modify_IC
     ! reset AMR grid
     if (resetgrid) call settree
     ! Select active grids
     call selectgrids

  else
     ! order: first physics dependent part then user for defaults.
     call phys_check_params()

     {#IFDEF RAY
     call init_rays
     }

     ! form and initialize all grids at level one
     call initlevelone
     ! set up and initialize finer level grids, if needed
     call settree

     {#IFDEF PARTICLES
     call init_tracerparticles
     call getbc(global_time,0.d0,pw,0,nwflux+nwaux)
     call init_gridvars
     call init_particles
     call finish_gridvars
     }

  end if

  ! set up boundary flux conservation arrays
  if (levmax>levmin) call allocateBflux

  if (mype==0) then
     print*,'-------------------------------------------------------------------------------'
     write(*,'(a,f17.3,a)')' Startup phase took : ',MPI_WTIME()-time0,' sec'
     print*,'-------------------------------------------------------------------------------'
  end if

  ! an interface to allow user to do special things before the main loop
  if (associated(usr_before_main_loop)) &
       call usr_before_main_loop()

  time_advance=.true.

  ! do time integration of all grids on all levels
  call timeintegration()

  if (mype==0) then
     print*,'-------------------------------------------------------------------------------'
     write(*,'(a,f17.3,a)')' Finished AMRVAC in : ',MPI_WTIME()-time0,' sec'
     print*,'-------------------------------------------------------------------------------'
  end if

  {#IFDEF PARTICLES
  time_advance=.false.
  call finish_tracerparticles
  }
  call comm_finalize

contains

  subroutine timeintegration()
    use mod_timing
    use mod_advance, only: advance, process
    use mod_forest, only: nleafs_active
    use mod_global_parameters
    use mod_input_output, only: saveamrfile
    use mod_ghostcells_update

    integer :: level, ifile, fixcount, ncells_update, ncells_block
    logical :: alive

    time_in=MPI_WTIME()
    fixcount=1

    n_saves(filelog_:fileout_) = snapshotini

    do ifile=nfile,1,-1
       tsavelast(ifile)=global_time
       itsavelast(ifile)=it
    end do

    itmin=it
    ! the next two are used to keep track of the performance during runtime:
    itTimeLast=it
    timeLast=MPI_WTIME()

    call getbc(global_time, 0.d0, pw, 0, nwflux+nwaux)

    !  ------ start of integration loop. ------------------
    timeloop0=MPI_WTIME()
    timefirstbc=timeloop0-time_in
    time_bc=0.d0
    if (mype==0) then
       write(*,'(a,f12.3,a)')&
            ' BCs before Advance took : ',timefirstbc,' sec'
    end if
    ncells_block={(ixGhi^D-2*nghostcells)*}
    ncells_update=0
    time_evol : do
       call setdt()
       if(fixprocess) call process(it,global_time)

       timeio0=MPI_WTIME()
       do ifile=nfile,1,-1
          if(timetosave(ifile)) call saveamrfile(ifile)
       end do
       ! save a snapshot when a file name 'savenow' is present
       if(mype==0) inquire(file='savenow',exist=alive)
       if(npe>1) call MPI_BCAST(alive,1,MPI_LOGICAL,0,icomm,ierrmpi)
       if(alive) then
          if(mype==0) write(*,'(a,i7,a,i7,a,es12.4)') ' save a snapshot No.',&
               snapshot,' at it=',it,' global_time=',global_time
          call saveamrfile(1)
          call saveamrfile(2)
          call MPI_FILE_DELETE('savenow',MPI_INFO_NULL,ierrmpi)
       endif
       timeio_tot=timeio_tot+(MPI_WTIME()-timeio0)

       ! exit time loop criteria
       if (it>=itmax) exit time_evol
       if (global_time>=time_max) exit time_evol

       call advance(it)

       ! resetting of tree BEFORE IO and setdt
       timegr0=MPI_WTIME()
       if(ditregrid>1) then
          if(fixcount<ditregrid) then
             fixcount=fixcount+1
          else
             if (refine_max_level>1 .and. .not.(fixgrid(0))) call resettree
             fixcount=1
          endif
       else
          if (refine_max_level>1 .and. .not.(fixgrid(0))) call resettree
       endif
       timegr_tot=timegr_tot+(MPI_WTIME()-timegr0)

       it = it + 1
       global_time = global_time + dt
       if(addmpibarrier) call MPI_BARRIER(icomm,ierrmpi)

       if(it>9000000)then
          it = slowsteps+10
          itmin=0
          itsavelast(:)=0
       end if
       ncells_update=ncells_update+ncells_block*nleafs_active
    end do time_evol

    timeloop=MPI_WTIME()-timeloop0

    if (mype==0) then
       write(*,'(a,f12.3,a)')' Total timeloop took        : ',timeloop,' sec'
       write(*,'(a,f12.3,a)')' Time spent on Regrid+Update: ',timegr_tot,' sec'
       write(*,'(a,f12.2,a)')'                  Percentage: ',100.0*timegr_tot/timeloop,' %'
       write(*,'(a,f12.3,a)')' Time spent on IO in loop   : ',timeio_tot,' sec'
       write(*,'(a,f12.2,a)')'                  Percentage: ',100.0*timeio_tot/timeloop,' %'
       write(*,'(a,f12.3,a)')' Time spent on BC           : ',time_bc,' sec'
       write(*,'(a,f12.2,a)')'                  Percentage: ',100.0*time_bc/timeloop,' %'
       write(*,'(a,f12.3,a)')' Time spent on run          : ',timeloop-timeio_tot,' sec'
       write(*,'(a,es16.8,a)')' Cells_updated / cpu / sec  : ',dble(ncells_update*nstep)/timeloop/dble(npe)
    end if

    timeio0=MPI_WTIME()
    do ifile=nfile,1,-1
       if(itsavelast(ifile)<it)call saveamrfile(ifile)
    enddo
    if (mype==0) call MPI_FILE_CLOSE(log_fh,ierrmpi)
    timeio_tot=timeio_tot+(MPI_WTIME()-timeio0)

    {#IFDEF RAY
    call time_spent_on_rays
    }
    {#IFDEF PARTICLES
    call time_spent_on_particles
    }

    if (mype==0) then
       write(*,'(a,f12.3,a)')' Total time spent on IO     : ',timeio_tot,' sec'
       write(*,'(a,f12.3,a)')' Total timeintegration took : ',MPI_WTIME()-time_in,' sec'
    end if

  end subroutine timeintegration

  !> Save times are defined by either tsave(isavet(ifile),ifile) or
  !> itsave(isaveit(ifile),ifile) or dtsave(ifile) or ditsave(ifile)
  !> Other conditions may be included.
  logical function timetosave(ifile)
    use mod_global_parameters

    integer:: ifile
    logical:: oksave
    !-----------------------------------------------------------------------------
    oksave=.false.
    if (it==itsave(isaveit(ifile),ifile)) then
       oksave=.true.
       isaveit(ifile)=isaveit(ifile)+1
    end if
    if (it==itsavelast(ifile)+ditsave(ifile)) oksave=.true.

    if (global_time>=tsave(isavet(ifile),ifile)) then
       oksave=.true.
       isavet(ifile)=isavet(ifile)+1
    end if
    if (global_time>=tsavelast(ifile)+dtsave(ifile)-smalldouble)then
       oksave=.true.
       n_saves(ifile) = n_saves(ifile) + 1
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
  !> @todo Fix dummy argument?
  logical function fixgrid(dummy)
    use mod_global_parameters
    integer :: dummy              !< Unused dummy argument

    fixgrid= (global_time>=tfixgrid .or. it>=itfixgrid)
  end function fixgrid

end program amrvac
