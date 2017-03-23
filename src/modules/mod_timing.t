!=============================================================================
module mod_timing
! module file for counters of wallclock time.
! 14.05.2012 Oliver Porth
! to use the timers, just 
! use mod_timing
!============================================================================
implicit none
save
double precision       :: time_in, timeio0, timeio_tot=0.0d0
double precision       :: timegr0, timegr_tot=0.0d0, timeloop, timeloop0
double precision       :: tpartc=0.0d0, tpartc_io=0.0d0, tpartc_int=0.0d0, tpartc_com=0.0d0, tpartc_grid=0.0d0
double precision       :: tpartc0, tpartc_int_0, tpartc_com0, tpartc_io_0, tpartc_grid_0

integer                :: itTimeLast
double precision       :: timeLast
end module mod_timing
!=============================================================================
