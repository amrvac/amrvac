!> This is a custom analysis routine. When it is called is controlled by
!> dtsave(5), ditsave(5) etc. in the savelist part of your .par file.
subroutine write_analysis()
  use mod_global_parameters
  use mod_input_output

  integer, parameter          :: n_modes   = 2
  logical, save               :: file_open = .false.
  character(len=*), parameter :: fname = "analysis.txt"
  integer                     :: power
  double precision            :: volume, kin_en_y

  call get_volume_average_func(get_kin_en_y, kin_en_y, volume)

  ! Only the root task writes output
  if (mype == 0) then

     ! On the first call, open the analysis output file and write a header
     if (.not. file_open) then
        open(unitanalysis, file=fname)
        write(unitanalysis, *) "# time kin_en_y"
        file_open = .true.
     end if

     ! Write them to the analysis file and 'flush' the unit, so that the file is
     ! immediately updated
     write(unitanalysis, *) t, kin_en_y
     flush(unitanalysis)
  end if

contains

  pure function get_kin_en_y(w_vec, w_size) result(val)
    integer, intent(in)          :: w_size
    double precision, intent(in) :: w_vec(w_size)
    double precision             :: val

    val = w_vec(m2_)**2 / w_vec(rho_)
  end function get_kin_en_y

end subroutine write_analysis
