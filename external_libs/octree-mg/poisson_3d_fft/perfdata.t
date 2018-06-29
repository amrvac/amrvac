!value of ncache for the FFT routines
integer, parameter :: ncache_optimal=%ncache%*1024
!flag that states if we must perform the timings or not
!definitions of the timing variambles just in case
integer, parameter :: timing_flag=%flag%
integer :: count_time1,count_time2,count_rate,count_max,number_time,index_time
real(kind=8) :: time_b0,time_b1,serial_time,parallel_time
