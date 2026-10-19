module mod_timing
! module file for counters of wallclock time.
! 14.05.2012 Oliver Porth
! to use the timers, just
! use mod_timing
  use mpi, only: MPI_WTIME
  !$ use omp_lib, only: omp_in_parallel
  implicit none
  save
  double precision       :: time_in, timeio0, timeio_tot=0.0d0
  double precision       :: timegr0, timegr_tot=0.0d0, timeloop, timeloop0
  double precision       :: timeeos_update=0.0d0
  double precision       :: timeeos_csound=0.0d0
  double precision       :: timeeos_conv=0.0d0
  double precision       :: timeeos_pthermal=0.0d0
  double precision       :: timeeos_Tfromei=0.0d0
  double precision       :: time_wb_transform=0.0d0, time_wb_inverse=0.0d0
  double precision       :: time_wb_recon=0.0d0
  double precision       :: time_rt=0.0d0, time_sts=0.0d0
  double precision       :: tpartc=0.0d0, tpartc_io=0.0d0, tpartc_int=0.0d0, tpartc_com=0.0d0, tpartc_grid=0.0d0
  double precision       :: tpartc0, tpartc_int_0, tpartc_com0, tpartc_io_0, tpartc_grid_0
  double precision       :: timeLast
  double precision       :: tw_advance=0.0d0, tw_setdt=0.0d0
  double precision       :: tw_regrid=0.0d0, tw_save=0.0d0
  double precision       :: tw_process=0.0d0, tw_sts=0.0d0
  double precision       :: tw_sts_prev=0.0d0, tw_tmp=0.0d0
  double precision       :: tw_eos_hydro=0.0d0, tw_eos_hydro_prev=0.0d0
  double precision       :: tw_eos_tc=0.0d0, tw_eos_tc_prev=0.0d0
  integer                :: itTimeLast
  integer, parameter     :: timing_unit=55
  logical                :: timing_log_opened=.false.
  !> Enable writing of detailed timing breakdown log (set via savelist)
  logical                :: write_timing_log=.false.
  !> Write timing log every N iterations (default 10)
  integer                :: timing_log_interval=10
  !> nesting depth of the EoS timers: only the outermost call accumulates, so an EoS call
  !> made inside another is counted once, under the category of the outer one
  integer, private          :: eos_depth=0
  double precision, private :: eos_t0

contains

  !> Open an EoS timing window; inside an OpenMP parallel region the EoS is not timed
  subroutine eos_timer_start()
    !$ if(omp_in_parallel()) return
    if(eos_depth==0) eos_t0=MPI_WTIME()
    eos_depth=eos_depth+1
  end subroutine eos_timer_start

  !> Close an EoS timing window, adding its time to acc when it is the outermost
  subroutine eos_timer_stop(acc)
    double precision, intent(inout) :: acc
    !$ if(omp_in_parallel()) return
    eos_depth=eos_depth-1
    if(eos_depth==0) acc=acc+(MPI_WTIME()-eos_t0)
  end subroutine eos_timer_stop

  !> All EoS time so far; windows that enclose EoS calls subtract its change
  double precision function time_eos()
    time_eos=timeeos_update+timeeos_Tfromei+timeeos_csound+timeeos_conv+timeeos_pthermal
  end function time_eos

  !> Open a timing window, whose length leaves out the EoS and ghost-cell time spent inside it,
  !> since those are counted under their own categories
  subroutine timer_open(t0)
    use mod_global_parameters, only: time_bc
    double precision, intent(out) :: t0(3)
    t0 = [MPI_WTIME(), time_eos(), time_bc]
  end subroutine timer_open

  !> Wall time since timer_open, less the EoS and ghost-cell time spent in between
  double precision function timer_close(t0)
    use mod_global_parameters, only: time_bc
    double precision, intent(in) :: t0(3)
    timer_close = (MPI_WTIME()-t0(1)) - (time_eos()-t0(2)) - (time_bc-t0(3))
  end function timer_close

  !> Zero the category timers of the time loop, so that the summary covers the loop only
  subroutine timers_reset()
    use mod_global_parameters, only: time_bc
    time_bc=0.d0; time_rt=0.d0; time_sts=0.d0
    timeeos_update=0.d0; timeeos_Tfromei=0.d0; timeeos_csound=0.d0
    timeeos_conv=0.d0; timeeos_pthermal=0.d0
    time_wb_transform=0.d0; time_wb_recon=0.d0; time_wb_inverse=0.d0
  end subroutine timers_reset

end module mod_timing
