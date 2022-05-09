module mod_usr
  use mod_twofl

  implicit none


  type arrayptr
    double precision,dimension(:), pointer :: p
    character(len=30) ::  namevar
  end type arrayptr
  
  double precision, parameter ::  pi = 3.14159

  !user defined vars
  double precision ::  H_ion_fr=1d-2


contains

  function get_number_lines(filename)
    integer, parameter :: MAX_LINE_LENGTH=1024
    character(len=MAX_LINE_LENGTH) line
    character(len=*), intent(in) :: filename
    integer, parameter :: lunit=338
    integer :: get_number_lines

    get_number_lines = 0
    OPEN (lunit, file = filename)
    DO
        READ (lunit,*, END=10) line
        if (line(1:1)=='#') cycle
        get_number_lines = get_number_lines + 1
    END DO
    10 CLOSE (lunit)
  end function get_number_lines

  subroutine read_formatted_file(filename,nrows,ncols,arr)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: nrows, ncols
    double precision, allocatable, dimension(:,:) :: arr

    integer, parameter :: MAX_LINE_LENGTH=1024
    character(len=MAX_LINE_LENGTH) line
    integer :: i,j,io

    integer, parameter :: reading_unit=100

    double precision, dimension(ncols) :: test_array
  

    OPEN(reading_unit, file=filename)
    
    i=1
    DO
      READ(reading_unit,'(A)',iostat=io) line
      IF (io/=0) exit
      if (line(1:1)=='#') cycle
      READ(line,*,iostat=io) test_array(1:ncols)
      if(io==-1) exit
      arr(i,:) = test_array
      if(i==nrows) exit
      i=i+1
    ENDDO
    CLOSE(reading_unit)
  end subroutine read_formatted_file

  function get_number_columns(filename)
    character(len=*), intent(in) :: filename
    integer :: get_number_columns
  
    integer, parameter :: MAX_NUM_OF_COLS=30
    integer, parameter :: MAX_LINE_LENGTH=1024
    character(len=MAX_LINE_LENGTH) line
    double precision, dimension(MAX_NUM_OF_COLS) :: test_array
    
    integer, parameter :: reading_unit=100
    integer i, io

    get_number_columns = -1

    OPEN(reading_unit, file=filename)
    
    DO
      READ(reading_unit,'(A)',iostat=io) line
      IF (io/=0) exit
      if (line(1:1)=='#') cycle
      do i=1,MAX_NUM_OF_COLS
        READ(line,*,iostat=io) test_array(1:i)
        if(io==-1) exit
      enddo
      if (get_number_columns .eq. -1) then 
        get_number_columns =  (i-1)
      else
        if (get_number_columns .ne. i-1) then
           print*, "Number of columns not the same for all the rows"
           exit
         endif 
      endif
    ENDDO
    CLOSE(reading_unit)
    
  end  function get_number_columns

    subroutine usrspecial_convert(qunitconvert)
    integer, intent(in) :: qunitconvert
    character(len=20):: userconvert_type

    end subroutine usrspecial_convert


  !> Read this module's parameters from a file
  subroutine usr_params_read(files)
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /usr_list/H_ion_fr 

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, usr_list, end=111)
111    close(unitpar)
    end do


  end subroutine usr_params_read



  subroutine usr_init()
    !use mod_variables

    call usr_params_read(par_files)

    if(mype .eq. 0) then
      print*, "H_ion_fr ", H_ion_fr
    endif

    usr_init_one_grid => initonegrid_usr
    usr_var_for_errest => special_errest

    usr_special_convert => usrspecial_convert
    call set_coordinate_system("Cartesian_2D")

    call twofl_activate()

  end subroutine usr_init
 
  !refine taking into account only vertical velocity
  subroutine special_errest(ixI^L,ixO^L,iflag,w,x,var)
    integer, intent(in)           :: ixI^L,ixO^L,iflag
    double precision, intent(in)  :: w(ixI^S,1:nw),x(ixI^S,1:ndim)
    double precision, intent(out) :: var(ixI^S)

    var(ixO^S) = sqrt((w(ixO^S,mom_c(1))/w(ixO^S,rho_c_)-w(ixO^S,mom_n(1))/w(ixO^S,rho_n_))**2+&
                     (w(ixO^S,mom_c(2))/w(ixO^S,rho_c_)-w(ixO^S,mom_n(2))/w(ixO^S,rho_n_))**2)

  end subroutine special_errest

 

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    ! initialize one grid
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision,parameter :: rho1 = 1d0
    double precision,parameter :: rho2 = 1.5d0
    double precision,parameter :: VV = 0.2
    double precision,parameter :: v1 = -1.5d0/2.5d0
    double precision,parameter :: v2 = 1d0/2.5d0
    double precision,parameter :: beta = 2d2
    double precision :: pe 
    double precision :: Bx 
    double precision :: ampl 
    double precision :: rho_nc_fr 
    character(len=*), parameter :: rfilename="random.txt"
    double precision, allocatable, dimension(:,:) :: rData
    integer :: nn, ncols, nrows
    double precision ::  pert(ixO^S)

    pe = 1d0/twofl_gamma
    Bx = sqrt(2d0/(twofl_gamma * beta))
    rho_nc_fr = (1d0 - H_ion_fr)/H_ion_fr

    w(ixO^S,mom_c(2)) = 0d0
    w(ixO^S,mom_n(2)) = 0d0
    where(x(ixO^S,2)<0d0)
      w(ixO^S,rho_c_) = rho1/(1d0 + rho_nc_fr)    
      w(ixO^S,rho_n_) = (rho1 * rho_nc_fr)/(1d0 + rho_nc_fr)    
      w(ixO^S,mom_c(1)) = v1 * VV
      w(ixO^S,mom_n(1)) = v1 * VV
    elsewhere
      w(ixO^S,rho_c_) = rho2/(1d0 + rho_nc_fr)    
      w(ixO^S,rho_n_) = (rho2 * rho_nc_fr)/(1d0 + rho_nc_fr)    
      w(ixO^S,mom_c(1)) = v2 * VV
      w(ixO^S,mom_n(1)) = v2 * VV
    endwhere

    w(ixO^S,e_c_) = (pe * Rc *  w(ixO^S,rho_c_))/ (Rc *  w(ixO^S,rho_c_) + Rn *  w(ixO^S,rho_n_))
    w(ixO^S,e_n_) = (pe * Rn *  w(ixO^S,rho_n_))/ (Rc *  w(ixO^S,rho_c_) + Rn *  w(ixO^S,rho_n_))
    !w(ixO^S,mag(1)) = Bx * sqrt(rho_nc_fr)
    w(ixO^S,mag(1)) = Bx 
    w(ixO^S,mag(2)) = 0d0

    ampl = 0.1d-2*sqrt(twofl_gamma * pe/rho1)

    nrows = get_number_lines(rfilename)
    ncols = get_number_columns(rfilename)
    if(ncols .ne. 2) then
      call mpistop("Number columns in random file must be 2")
    endif
    allocate(rData(nrows,ncols))
    call read_formatted_file(rfilename,nrows,ncols,rData)

    pert(ixO^S) = 0d0
    do nn=1,nrows
     pert(ixO^S) =  pert(ixO^S) + rData(nn,1)  * sin(2*pi*x(ixO^S,1)*nn/(xprobmax1-xprobmin1) + rData(nn,2))
    enddo

    where(x(ixO^S,2).ge. -0.5 * dxlevel(2) .and. x(ixO^S,2).lt. 0.5 * dxlevel(2))
      !pert in neutrals
      w(ixO^S,mom_n(2)) = ampl *  pert(ixO^S)
      !same pert in charges
      w(ixO^S,mom_c(2)) = w(ixO^S,mom_n(2)) 
    endwhere

    deallocate(rData)

    call twofl_to_conserved(ixI^L,ixO^L,w,x)


  end subroutine initonegrid_usr








end module mod_usr
