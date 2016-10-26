!> This is a custom analysis routine. When it is called is controlled by
!> dtsave(5), ditsave(5) etc. in the savelist part of your .par file.
subroutine write_analysis()
  include 'amrvacdef.f'

  integer, parameter          :: n_modes   = 2
  logical, save               :: file_open = .false.
  character(len=*), parameter :: fname = "analysis.txt"
  integer                     :: power
  double precision            :: modes(nw, n_modes), volume
  double precision            :: kin_en_y

  ! Compute the volume average of w**1 and w**2
  do power = 1, n_modes
     call get_volume_average(power, modes(:, power), volume)
  end do

  ! Only the root task writes output
  if (mype == 0) then

     ! On the first call, open the analysis output file and write a header
     if (.not. file_open) then
        open(unitanalysis, file=fname)
        write(unitanalysis, *) "# time kin_en_y"
        file_open = .true.
     end if

     ! Compute the kinetic energy in the y-direction
     kin_en_y = modes(m2_, 2) / modes(rho_, 1)

     ! Write them to the analysis file and 'flush' the unit, so that the file is
     ! immediately updated
     write(unitanalysis, *) t, kin_en_y
     flush(unitanalysis)
  end if

end subroutine write_analysis
