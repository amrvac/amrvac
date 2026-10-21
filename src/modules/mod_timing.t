module mod_timing
! module file for counters of wallclock time.
! 14.05.2012 Oliver Porth
! to use the timers, just 
! use mod_timing
  implicit none
  save
  double precision       :: time_in, timeio0, timeio_tot=0.0d0
  double precision       :: timegr0, timegr_tot=0.0d0, timeloop, timeloop0
  double precision       :: timeeos0, timeeos_tot=0.0d0
  double precision       :: timeeos_update=0.0d0
  double precision       :: timeeos_csound=0.0d0
  double precision       :: timeeos_conv=0.0d0
  double precision       :: timeeos_pthermal=0.0d0
  double precision       :: timeeos_Tfromei=0.0d0
  double precision       :: time_wb_transform=0.0d0, time_wb_inverse=0.0d0
  double precision       :: time_wb_recon=0.0d0
  double precision       :: tpartc=0.0d0, tpartc_io=0.0d0, tpartc_int=0.0d0, tpartc_com=0.0d0, tpartc_grid=0.0d0
  double precision       :: tpartc0, tpartc_int_0, tpartc_com0, tpartc_io_0, tpartc_grid_0
  double precision       :: timeLast
  double precision       :: tw_advance=0.0d0, tw_setdt=0.0d0
  double precision       :: tw_regrid=0.0d0, tw_save=0.0d0
  double precision       :: tw_process=0.0d0, tw_tc_total=0.0d0
  double precision       :: tw_tc_prev=0.0d0, tw_tmp=0.0d0
  double precision       :: tw_eos_hydro=0.0d0, tw_eos_hydro_prev=0.0d0
  double precision       :: tw_eos_tc=0.0d0, tw_eos_tc_prev=0.0d0
  double precision       :: time_htc_total=0.0d0, time_htc0=0.0d0
  integer                :: itTimeLast
  integer, parameter     :: timing_unit=55
  logical                :: timing_log_opened=.false.
  !> Enable writing of detailed timing breakdown log (set via savelist)
  logical                :: write_timing_log=.false.
  !> Write timing log every N iterations (default 10)
  integer                :: timing_log_interval=10
end module mod_timing
