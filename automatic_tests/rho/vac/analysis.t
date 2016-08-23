!> Write the first and second mode to a special log file
subroutine write_analysis
  include 'amrvacdef.f'

  integer, parameter :: n_modes   = 2
  character(len=40)  :: fmt_string
  logical, save      :: file_open = .false.
  integer            :: power
  double precision   :: modes(nw, n_modes), volume

  do power = 1, n_modes
     call volume_average(power, modes(:, power), volume)
  end do

  if (mype == 0) then
     if (.not. file_open) then
        open(unit=unitanalysis, file=trim(inifile)//".log", status='unknown')
        file_open = .true.
        write(unitanalysis, *) "# time mean(w) mean(w**2)"
     end if

     write(fmt_string, "(A,I0,A)") "(", nw * n_modes + 1, "es14.6)"
     write(unitanalysis, fmt_string) t, modes
  end if

end subroutine write_analysis

!> Compute mean(w**power) over the leaves of the grid. The first mode (power=1)
!> corresponds to to the mean, the second to the mean squared values etc.
subroutine volume_average(power, mode, volume)
  use constants
  include 'amrvacdef.f'

  integer, intent(in)           :: power     !< Which mode to compute
  double precision, intent(out) :: mode(nw)  !< The computed mode
  double precision, intent(out) :: volume    !< The total grid volume
  integer                       :: iigrid, igrid, level, iw
  double precision              :: wsum(nw+1)
  double precision              :: dvolume(ixG^T)
  double precision              :: dsum_recv(1:nw+1)

  wsum(:) = 0

  ! Loop over all the grids
  do iigrid = 1, igridstail
     igrid = igrids(iigrid)
     level = node(plevel_,igrid)

     ! Determine the volume of the grid cells
     if (slab) then
        dvolume(ixM^T) = {rnode(rpdx^D_,igrid)|*}
     else
        dvolume(ixM^T) = pgeo(igrid)%dvolume(ixM^T)
     end if

     ! Store total volume in last element
     wsum(nw+1) = wsum(nw+1) + sum(dvolume(ixM^T))

     ! Compute the modes of the cell-centered variables, weighted by volume
     do iw = 1, nw
        wsum(iw) = wsum(iw) + &
             sum(dvolume(ixM^T)*pw(igrid)%w(ixM^T,iw)**power)
     end do
  end do

  ! Make the information available on all tasks
  call MPI_ALLREDUCE(wsum, dsum_recv, nw+1, MPI_DOUBLE_PRECISION, &
       MPI_SUM, 0, icomm, ierrmpi)

  ! Set the volume and the average
  volume = wsum(nw+1)
  mode   = wsum(1:nw) / volume

end subroutine volume_average
