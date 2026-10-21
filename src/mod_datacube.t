module mod_datacube
  use, intrinsic :: ieee_arithmetic
  use mod_constants, only : dpi
  implicit none
  double precision :: globaltime
  integer :: i, found
  double precision :: ix2, ix3
  integer :: pos_index
  integer, parameter :: nclt = 90, nlons = 180

contains

  subroutine get_datacube_arrays(target_time, cme_number, base_file, found, clts, lons, mask, vr, vp, vt, br, bt, bp, md, temp)
    implicit none
    double precision, intent(in) :: target_time
    integer, intent(in) :: cme_number
    character(len=*), intent(in) :: base_file
    integer, intent(out) :: found
    double precision,    allocatable, intent(out) :: clts(:), lons(:)
    real,    allocatable, intent(out) :: vr(:),  vt(:), vp(:)
    real,    allocatable, intent(out) :: br(:),  bt(:), bp(:)
    real,    allocatable, intent(out) :: md(:),  temp(:),  mask(:)
    character(len=256) :: binfile, cltsfile, lonsfile
    character(len=256) :: output_dir, fullpath, cme_str
    integer :: unit, ios, size(9), size_2(9)
    character(len=256) :: cmd, line, filename, bestfile
    character(len=256), allocatable :: tmpfile
    integer :: ierr, nfiles, n, i, j, ndata,  pos
    double precision :: file_time
    logical :: file_exists, in_mask
    integer :: unit_tmp, file_num
    character(len=256) :: base_name, number_str
    character(len=256), allocatable :: files(:)
    integer :: clt_size, lon_size
    double precision :: min_diff, latest_time
character(len=20) :: time_str
integer :: time_int
    ! Initialize outputs
    found = 0
    if (allocated(clts)) deallocate(clts)
    if (allocated(lons)) deallocate(lons)
    if (allocated(mask)) deallocate(mask)
    if (allocated(vr)) deallocate(vr)
    if (allocated(vp)) deallocate(vp)
    if (allocated(vt)) deallocate(vt)
    if (allocated(br)) deallocate(br)
    if (allocated(bt)) deallocate(bt)
    if (allocated(bp)) deallocate(bp)
    if (allocated(md)) deallocate(md)
    if (allocated(temp)) deallocate(temp)
    ! --- Step 1: List all  files in datacube folder ---
    pos = index(base_file, '/', back=.true.)
    if (pos > 0) then
        output_dir = base_file(1:pos)
    else
        output_dir = './'
    end if

    write(cme_str, '(I0)') cme_number
    fullpath = trim(output_dir)//'datacube/cme_'//trim(cme_str)
    tmpfile=trim(fullpath)//'/file_list.txt'

    ! Check if temp file exists and has content
    inquire(file=trim(tmpfile), exist=file_exists)
    if (.not. file_exists) then
        print *, "No file list file found in directory: ", trim(fullpath)
        return
    end if

    ! Count number of files
    open(newunit=unit_tmp, file=tmpfile, status='old', action='read', iostat=ierr)
    if (ierr /= 0) then
        print *, "Error opening temp file"
        return
    end if

    nfiles = 0
    do
        read(unit_tmp, '(A)', iostat=ios) line
        if (ios /= 0) exit
        if (trim(line) /= '') nfiles = nfiles + 1
    end do

    if (nfiles == 0) then
        close(unit_tmp)
        print *, "No valid .bin files found"
        return
    end if

    rewind(unit_tmp)

    ! Read file names
    allocate(files(nfiles), stat=ierr)
    if (ierr /= 0) then
        print *, "Error allocating files array"
        close(unit_tmp)
        return
    end if

    j = 1
    do
        read(unit_tmp, '(A)', iostat=ios) line
        if (ios /= 0) exit
        if (trim(line) /= '') then
            files(j) = trim(line)
            j = j + 1
        end if
    end do
    close(unit_tmp)

    ! --- Step 2: Find latest file time ---
    latest_time = -1.0d0
    do j = 1, nfiles
        pos = index(files(j), '/', back=.true.)
        if (pos == 0) then
            base_name = trim(files(j))
        else
            base_name = trim(files(j)(pos+1:))
        end if

        pos = index(base_name, '.bin')
        if (pos > 1) then
            number_str = base_name(1:pos-1)
            read(number_str, *, iostat=ierr) file_num
            file_time = file_num/60.0d0
            if (ierr == 0) then
                if (file_time > latest_time) latest_time = file_time
            end if
        end if
    end do
    ! --- Step 3: Check if input time is within range ---
    if (target_time > latest_time) then
        ! print *, "Target time exceeds latest available data time"
        deallocate(files)
        return
    end if

    ! --- Step 4: Find closest file by time ---
    min_diff = huge(min_diff)
    bestfile = ""
    do j = 1, nfiles
        pos = index(files(j), '/', back=.true.)
        if (pos == 0) then
            base_name = trim(files(j))
        else
            base_name = trim(files(j)(pos+1:))
        end if

        pos = index(base_name, '.bin')
        if (pos > 1) then
            number_str = base_name(1:pos-1)
            read(number_str, *, iostat=ierr) file_time
            file_time = file_time/60.0d0
            if (ierr == 0) then
                if (abs(file_time - target_time) < min_diff) then
                    min_diff = abs(file_time - target_time)
                    bestfile = files(j)
                end if
            end if
        end if
    end do

    if (len_trim(bestfile) == 0) then
        print *, "Could not determine best matching file"
        deallocate(files)
        return
    end if

  ! --- Step 5: Read datacube files ---

    ! Extract output directory from base_file path
    pos = index(base_file, '/', back=.true.)
    if (pos > 0) then
      output_dir = base_file(1:pos)
    else
      output_dir = './'
    end if

  write(cme_str, '(I0)') cme_number

  binfile = trim(output_dir) // 'datacube/cme_' // trim(cme_str) // '/'  // trim(bestfile)
    inquire(file=binfile, exist=file_exists)
    if (.not. file_exists) then
      print *, "Binary datacube file does not exist: ", trim(binfile)
      return
    end if

    ! Open binary file to read
    open(newunit=unit, file=binfile, status='old', access='stream', form='unformatted', action='read', iostat=ios)
    if (ios /= 0) then
      print *, "Error opening datacube binary file: ", trim(binfile)
      return
    end if
    read(unit) size

    ! Allocate all arrays of length size
    allocate(mask(size(1)), vr(size(2)), vp(size(3)), vt(size(4)), br(size(5)), bt(size(6)), bp(size(7)), md(size(8)), temp(size(9)), stat=ios)
    if (ios /= 0) then
      print *, "Error allocating arrays of size ", size
      close(unit)
      return
    end if

    ! Read arrays sequentially in the fixed order
    read(unit) mask
    read(unit) vr
    read(unit) vt
    read(unit) vp
    read(unit) br
    read(unit) bt
    read(unit) bp
    read(unit) md
    read(unit) temp

    close(unit)


    ! Read clts and lons arrays from separate binary files
    clt_size = 90
    lon_size = 180


    allocate(clts(clt_size), lons(lon_size))
    do i = 1, (clt_size)
        clts(i) = dble(i-1) * dpi / dble(clt_size)
    end do

    do i = 1, (lon_size)
        lons(i) = -dpi + dble(i-1) * 2*dpi / dble(lon_size)
    end do
    found = 1  ! success
  end subroutine get_datacube_arrays



  subroutine get_value(clt, lon_in, ts_magnetogram, clts, lons, mask, vr, vp, vt, br, bt, bp, temp, md, ndata, &
                      found, vr_val, vp_val, vt_val, br_val, bt_val, bp_val, md_val, temp_val)
    implicit none
    double precision, intent(in) :: clt, lon_in, ts_magnetogram
    double precision, intent(in) :: clts(:), lons(:)
    real, intent(in) :: mask(:), vr(:), vp(:), vt(:), br(:), bt(:), bp(:), temp(:), md(:)
    integer, intent(in) :: ndata
    integer, intent(out) :: found
    double precision, intent(out) :: vr_val, vp_val, vt_val, br_val, bt_val, bp_val, md_val, temp_val
    integer :: i, j, k
    real :: val_1(8), val_2(8), val_min(8), val_max(8)
    double precision :: synoptic_correction, lon
    
    ! Correct synoptic rotation
    synoptic_correction = 2*dpi/27.27/24.0d0
    
    lon = modulo(lon_in + synoptic_correction*ts_magnetogram + dpi, 2*dpi) -dpi
    ! Initialize outputs
    found = 0

    vr_val   = ieee_value(vr_val, ieee_quiet_nan)
    vp_val   = ieee_value(vp_val, ieee_quiet_nan)
    vt_val   = ieee_value(vt_val, ieee_quiet_nan)
    br_val   = ieee_value(br_val, ieee_quiet_nan)
    bt_val   = ieee_value(bt_val, ieee_quiet_nan)
    bp_val   = ieee_value(bp_val, ieee_quiet_nan)
    md_val   = ieee_value(md_val, ieee_quiet_nan)
    temp_val = ieee_value(temp_val, ieee_quiet_nan)

   
    ! Get position index
    i = get_position_index(clt, lon, clts, lons)
    if (i == -1) return

    ! Find matching data point
    do j = 1, ndata
        if (nint(mask(j)) == i) then
            found = 1
            vr_val   = vr(min(j, size(vr)))
            vp_val   = vp(min(j, size(vp)))
            vt_val   = vt(min(j, size(vt)))
            br_val   = br(min(j, size(br)))
            bt_val   = bt(min(j, size(bt)))
            bp_val   = bp(min(j, size(bp)))
            md_val   = md(min(j, size(md)))
            temp_val = temp(min(j, size(temp)))
            exit
        end if
    end do


  end subroutine get_value


  function get_position_index(clt_val, lon_val, clts, lons) result(position_index)
    implicit none
    double precision, intent(in) :: clt_val, lon_val
    double precision, intent(in) :: clts(:), lons(:)
    integer :: position_index
    integer :: clt_index, lon_index, n
    double precision :: current_diff, min_diff

    ! Find closest clt index
    clt_index = 1
    min_diff = abs(clts(1) - real(clt_val))
    do n = 2, size(clts)
        current_diff = abs(clts(n) - real(clt_val))
        if (current_diff < min_diff) then
            min_diff = current_diff
            clt_index = n
        end if
    end do

    ! Find closest lon index
    lon_index = 1
    min_diff = abs(lons(1) - real(lon_val))
    do n = 2, size(lons)
        current_diff = abs(lons(n) - real(lon_val))
        if (current_diff < min_diff) then
            min_diff = current_diff
            lon_index = n
        end if
    end do

    ! Convert to 0-based index
    lon_index = lon_index - 1
    clt_index = clt_index - 1

    if (clt_index >= 0 .and. lon_index >= 0) then
        position_index = lon_index + clt_index * size(lons)
    else
        position_index = -1
    end if
  end function get_position_index

  function round5(x) result(r)
  implicit none
  real, intent(in) :: x   ! input variable (double precision)
  real :: r               ! result
  real, parameter :: scale = 1.0d5

  ! Multiply, round to nearest integer, then divide back
  r = nint(x * scale) / scale
end function round5

end module mod_datacube
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   

                                                                                                                            
