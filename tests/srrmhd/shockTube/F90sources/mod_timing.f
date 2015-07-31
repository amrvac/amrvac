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
double precision       :: timegr0, timegr_tot=0.0d0, timeloop, timeloop0,&
    timefirstbc


integer                :: itTimeLast
double precision       :: timeLast
end module mod_timing
!=============================================================================
