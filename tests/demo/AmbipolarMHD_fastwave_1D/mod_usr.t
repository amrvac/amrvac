module mod_usr
  use mod_mhd


  implicit none


  double precision:: usr_grav
  !!user defined

  double precision:: Period=5
  double precision:: ampl=1d-3
  logical :: usePml = .false.
  double precision :: xL  = 1.7
  double precision :: sigmaPml =259 
  
  logical :: useFilter = .false.
  integer :: niter_filter = 1

  logical :: maskAmbi = .true.
  double precision :: xLambi = 1.65
  double precision :: wLambi = 0.01
  logical :: ambi_mask_smooth  = .true.
  integer, parameter :: MASK_DISC = 1
  integer, parameter :: MASK_TANH = 2
  integer, parameter :: MASK_PROF = 3
  integer :: ambi_mask_method = MASK_DISC
  logical :: random_ambi = .false.



  double precision :: fracEqP = 0.5 
  logical :: useFracEqP = .true.
  double precision :: b0Param = 5.0


  character(len=*), parameter :: VALCfilename="valc.txt"
  character(len=*), parameter :: equi_filename="myEqui.txt"
  integer, parameter :: equi_zz=1, equi_pe=2, equi_rho=3, equi_b=4 !indices in equi generated file


  type arrayptr
    double precision,dimension(:), pointer :: p
    character(len=30) ::  namevar
  end type arrayptr



contains

  !> Read this module's parameters from a file
  subroutine usr_params_read(files)
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /usr_list/ xL,sigmaPml,niter_filter,maskAmbi,xLambi, Period, ampl, &
              wLambi,ambi_mask_method,ambi_mask_smooth,usePml,useFilter,random_ambi,&
              fracEqP, useFracEqP,b0Param

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, usr_list, end=111)
111    close(unitpar)
    end do

  end subroutine usr_params_read


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




  subroutine write_formatted_file(filename,vars)
    character(len=*), intent(in) :: filename
    type(arrayptr), intent(in), dimension(:) :: vars
    integer :: i,j
    integer, parameter :: reading_unit=100
    OPEN(reading_unit, file=filename)
    write(reading_unit,'(A)',advance="no") "#"
    !!this should also work if the format is 'number of vars(A)'instead of just '(A)' which puts one value at one line
    !write(reading_unit,'(A)') (vars(i)%namevar,i=1,size(vars))
    write(reading_unit,*) (vars(i)%namevar,i=1,size(vars))
    do i=1,size(vars(1)%p)
      !write(reading_unit,'(E2.5)') (vars(j)%p(i),j=1,size(vars))
      write(reading_unit,*) (vars(j)%p(i),j=1,size(vars))
    enddo
    CLOSE(reading_unit)
  end subroutine write_formatted_file




!> get the number of columns from file filename
!> It ignores lines starting with #
!> It also tests that the number of columns is the same for all the rows and exits when it encounters
!> the first row with different number of columns 
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





  subroutine usr_set_equi_vars2(xnorm, pe, rho, b0)
    use mod_interpolation
    !use interp1D
    double precision, intent(in) :: xnorm(:)
    double precision, dimension(:) :: pe,rho,b0
 
    double precision, allocatable  :: xx(:)
    integer :: nrows, ncols
    !!!we already know there are 10 columns
    double precision, allocatable, dimension(:,:) :: valcData
    integer, parameter :: zz_ = 1
    integer, parameter :: te_ = 4
    integer, parameter :: nn_ = 6
    integer, parameter :: ne_ = 7
    integer, parameter :: rho_ = 10
    !!constants see mod_constants
    double precision, parameter ::  MH=1.67339e-27
    double precision, parameter :: BK=1.380662e-23
    !double precision, parameter ::  G=2.7398e2
    double precision, parameter ::  G=2.74d2  !!!!use exactly the same value as in usr
    double precision, parameter ::  MU0=4.0e-7*dpi
    !!constants end 

    double precision :: B00, nn00, nc00, pe_c00, pe_n00, pe00, sumTemp, expArg,dz,const
    integer :: nz,k, eqP

    double precision, allocatable, dimension(:) :: zzValc
    double precision, allocatable, dimension(:) :: tempValc
    integer, allocatable, dimension(:) :: ind

    double precision, allocatable :: te(:)

  
    if(size(xnorm)<2) return


    nrows = get_number_lines(VALCfilename)
    ncols = get_number_columns(VALCfilename)

    nz = size(xnorm)
    allocate(xx(nz))
    xx = xnorm * unit_length !convert normalized x to SI units (in m)

    allocate(valcData(nrows,ncols))
    call read_formatted_file(VALCfilename,nrows,ncols,valcData)
    allocate(zzValc(nrows))
    allocate(ind(nrows))
    allocate(tempValc(nrows))

    allocate(te(nz))
    zzValc=valcData(:,zz_)*1e3 !zz in the data in valc is km (transform to SI units: m)
 
    !print*, "ZZvalc ", zzValc
    ind = rargsort(zzValc)
    zzValc=zzValc(ind)

    call interp_cubic_spline(zzValc,valcData(ind,te_),nrows,xx,te,nz)
    
    !!reuse pe array to make the interpolation
    call interp_cubic_spline(zzValc,valcData(ind,nn_),nrows,xx,pe,nz)
    nn00 = pe(1) * 1e6 !!VALC norm
    !!reuse pe array to make the interpolation
    call interp_cubic_spline(zzValc,valcData(ind,ne_),nrows,xx,pe,nz)
    nc00 = 2*pe(1) * 1e6!!VALC norm



    pe_n00 = nn00 * BK * te(1)
    pe_c00 = nc00 *  BK * te(1)
    pe00 = pe_n00 + pe_c00
    pe00 = pe00/unit_pressure 
    te=te/unit_temperature
    dz=xnorm(2)-xnorm(1)
    sumTemp = 0d0
    do k=1,nz
      pe(k)=pe00*exp(usr_grav * dz * sumTemp)
      rho(k) = pe(k)/ te(k)
      sumTemp = sumTemp + 1.0/ te(k)
    enddo
    if (useFracEqP) then
      eqP = 1 + int((nz-1) * fracEqP)
      b0(:) = sqrt(2 * pe(eqP) )
    else
      b0(:) = b0Param
      b0=b0/unit_magneticfield
    endif  


  
    deallocate(valcData,zzValc,ind, tempValc)
    deallocate(te)
  end subroutine usr_set_equi_vars2


    subroutine usrspecial_convert(qunitconvert)
    integer, intent(in) :: qunitconvert
    character(len=20):: userconvert_type

    end subroutine usrspecial_convert


  subroutine set_equi_vars2(xnorm, pe, rho)
    use mod_interpolation
    !use interp1D
    double precision, intent(in) :: xnorm(:)
    double precision, dimension(:) :: pe,rho

    integer :: nrows, ncols, mz
    !!!we already know there are 10 columns
    double precision, allocatable, dimension(:,:) :: valcData
    nrows = get_number_lines(equi_filename)
    ncols = get_number_columns(equi_filename)

    mz = size(xnorm)

    allocate(valcData(nrows,ncols))
    call read_formatted_file(equi_filename,nrows,ncols,valcData)

    call interp_cubic_spline(valcData(:,equi_zz),valcData(:,equi_rho),nrows,xnorm,rho,mz)
    call interp_cubic_spline(valcData(:,equi_zz),valcData(:,equi_pe),nrows,xnorm,pe,mz)
   deallocate(valcData)
  end subroutine set_equi_vars2


  subroutine set_equi_vars2_b0(xnorm,bx0)
    use mod_interpolation
    !use interp1D
    double precision, intent(in) :: xnorm(:)
    double precision, dimension(:) :: bx0

    integer :: nrows, ncols, mz
    double precision, allocatable, dimension(:,:) :: valcData
    nrows = get_number_lines(equi_filename)
    ncols = get_number_columns(equi_filename)

    mz = size(xnorm)

    allocate(valcData(nrows,ncols))
    call read_formatted_file(equi_filename,nrows,ncols,valcData)

    call interp_cubic_spline(valcData(:,equi_zz),valcData(:,equi_b),nrows,xnorm,bx0,mz)
   deallocate(valcData)
  end subroutine set_equi_vars2_b0


  subroutine init_equi_vars(xnorm)
    double precision, allocatable, target, intent(in) :: xnorm(:)
    double precision, allocatable, dimension(:),target :: pe,rho,b0
    type(arrayptr), dimension(4) :: write_vars
    integer :: mz 
 
    mz= size(xnorm)

    allocate(pe(mz),rho(mz),b0(mz))
    call usr_set_equi_vars2(xnorm,pe,rho,b0)

    write_vars(equi_zz)%namevar = "z"
    write_vars(equi_pe)%namevar  = "pe"
    write_vars(equi_rho)%namevar  = "rho"
    write_vars(equi_b)%namevar  = "b0"

    write_vars(equi_zz)%p =>xnorm
    write_vars(equi_pe)%p =>pe
    write_vars(equi_rho)%p =>rho
    write_vars(equi_b)%p =>b0

    call  write_formatted_file(equi_filename,write_vars)
  
    deallocate(pe,rho,b0)
  end subroutine init_equi_vars


  subroutine dump_units()
    character(len=*), parameter :: units_filename="units.dat"
    type(arrayptr), dimension(7) :: write_vars
    integer, parameter :: I_LEN =1, I_TIME=2, I_DENS=3, I_VEL=4, I_PRES=5, I_MAGFIELD=6, I_TEMP=7

    write_vars(I_LEN)%namevar = "unit_length"
    write_vars(I_TIME)%namevar = "unit_time"
    write_vars(I_DENS)%namevar = "unit_density"
    write_vars(I_VEL)%namevar = "unit_velocity"
    write_vars(I_PRES)%namevar = "unit_pressure"
    write_vars(I_MAGFIELD)%namevar = "unit_magneticfield"
    write_vars(I_TEMP)%namevar = "unit_temperature"

    allocate(write_vars(I_LEN)%p(1), write_vars(I_TIME)%p(1), write_vars(I_DENS)%p(1), write_vars(I_VEL)%p(1), write_vars(I_PRES)%p(1), write_vars(I_MAGFIELD)%p(1), write_vars(I_TEMP)%p(1))
    write_vars(I_LEN)%p(1) = unit_length
    write_vars(I_TIME)%p(1) = unit_time
    write_vars(I_DENS)%p(1) = unit_density
    write_vars(I_VEL)%p(1) = unit_velocity
    write_vars(I_PRES)%p(1) = unit_pressure
    write_vars(I_MAGFIELD)%p(1) = unit_magneticfield
    write_vars(I_TEMP)%p(1) = unit_temperature
    call  write_formatted_file(units_filename,write_vars)

  deallocate(write_vars(I_LEN)%p, write_vars(I_TIME)%p, write_vars(I_DENS)%p, write_vars(I_VEL)%p, write_vars(I_PRES)%p, write_vars(I_MAGFIELD)%p, write_vars(I_TEMP)%p)
  end subroutine dump_units


  subroutine usr_init()
    !use mod_variables

    unit_length        = 1.d6                                         ! m
    unit_temperature   = 5d3                                         ! K
    unit_numberdensity = 1.d20                                        !/m^3

    call usr_params_read(par_files)

    if(mype .eq. 0) then
      print*, "Period ", Period
      print*, "Gamma ", mhd_gamma
      print*, "Amplitude ", ampl
    endif

    usr_init_one_grid => initonegrid_usr
    usr_set_equi_vars => special_set_equi_vars
    usr_set_B0              => specialset_B0
    usr_set_J0              => specialset_J0
    usr_special_bc    => specialbound_usr
    usr_gravity         => gravity

    usr_set_parameters  => init_params_usr

    !!FILTER
    if(useFilter) then
      usr_process_adv_grid => special_process_filter
      !usr_process_grid => special_process_filter
    endif
    if (usePml) then
      !!as int bc
      !usr_internal_bc => special_process
      usr_process_grid => special_process
    endif
    if (maskAmbi) then
      usr_mask_ambipolar => special_ambipolar
    endif

    !usr_aux_output    => specialvar_output
    usr_special_convert => usrspecial_convert
    !usr_add_aux_names => specialvarnames_output
    call set_coordinate_system("Cartesian_1.75D")

    call mhd_activate()

  end subroutine usr_init

    !> Here one can add a steady (time-independent) equi vars
    subroutine special_set_equi_vars(ixI^L,ixO^L,x,w0)
      use mod_global_parameters
      integer, intent(in)             :: ixI^L,ixO^L
      double precision, intent(in)    :: x(ixI^S,1:ndim)
      double precision, intent(inout) :: w0(ixI^S,1:number_equi_vars)

    call set_equi_vars2(x(ixO^S,1), w0(ixO^S,equi_pe0_), w0(ixO^S,equi_rho0_))
    end subroutine special_set_equi_vars


  subroutine specialset_B0(ixI^L,ixO^L,x,wB0)
  ! Here add a time-independent background magnetic field
    integer, intent(in)           :: ixI^L,ixO^L
    double precision, intent(in)  :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: wB0(ixI^S,1:ndir)

    call set_equi_vars2_b0(x(ixO^S,1), wB0(ixO^S,3))

    wB0(ixO^S,1:2)=0d0
  end subroutine specialset_B0

  subroutine specialset_J0(ixI^L,ixO^L,x,wJ0)
  ! Here add a time-independent background current density 
    integer, intent(in)           :: ixI^L,ixO^L
    double precision, intent(in)  :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: wJ0(ixI^S,7-2*ndir:ndir)
    wJ0(ixI^S,:)=zero


  end subroutine specialset_J0
  
  subroutine special_ambipolar(ixI^L,ixO^L,w,x,ambiCoef)
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(in) :: w(ixI^S,1:nw)
    double precision, intent(inout) :: ambiCoef(ixI^S)

    !!METHOD 3
    double precision, allocatable, dimension(:^D&) :: tmp
    double precision, DIMENSION(8), PARAMETER :: prof=(/  1.0000000D0, 0.98985000D0, 0.89172687D0, 0.65310198D0, &
             0.34690427D0, 0.10828125D0, 0.010160819D0, 1.1351600D-05  /) 
    integer :: ixG^L
    logical :: found
    !!METHOD 3 end

    !integer, parameter :: bottomS=15  
    !integer, parameter :: topS=15
    !integer :: nn
    integer :: ii

    double precision  :: myRand


    !!this works better than tanh


    if(random_ambi) then
      if(mype .eq. 0) then
        myRand = rand(24)
      endif
      call MPI_BCAST(myRand,1,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
      xLambi = xLambi +  myRand * 0.15
    endif
    if(ambi_mask_method .eq. MASK_DISC .or. ambi_mask_method .eq. MASK_TANH) then
      if(ambi_mask_method .eq. MASK_DISC) then
        !!METHOD 1
        where(x(:,1) .ge. xLambi) 
          ambiCoef=0.0
        endwhere
      else
        !!METHOD 2
        ambiCoef(ixO^S) = ambiCoef(ixO^S) * 0.5 * (1.0 - tanh(  (x(ixO^S,1)-xLambi)/(wLambi) ))
  
      endif
      if(ambi_mask_smooth) then
        forall (ii = ixO^S)
          ambiCoef(ii) = (ambiCoef(ii-1) + ambiCoef(ii) + ambiCoef(ii+1))/3d0  
        endforall
      endif
    elseif (ambi_mask_method .eq. MASK_PROF) then   
 
    !!METHOD 3
      ixG^L=ixO^L^LADD1;
      allocate(tmp(ixG^S))
  
      if(x(ixOmax1,1) < xLambi) then
         tmp(ixG^S)=1d0
      else
        if(x(ixOmin1,1) < xLambi) then
          found = .false.
          ii=ixOmin1+1
          do while (.not. found .and. ii .le. ixOmax1)
            if(x(ii,1) .ge. xLambi) then
              found = .true.
              !print* , " II ", ii, " x ", x(ii,1)
              if(ii+7 .ge. ixGmax1) then
                call mpistop("No points por ambi tmp")
              endif
              tmp(ixGmin1:ii-1)=1d0
              tmp(ii:ii+7) = prof(1:8)
              tmp(ii+8:ixGmax1)=0d0
            endif 
            ii=ii+1
           enddo
        else
          tmp(ixG^S)=0D0
        endif 
      endif
  
      if(ambi_mask_smooth) then 
        tmp(ixG^S) = ambiCoef(ixG^S) * tmp(ixG^S)
        !!SMOOTH no temp array
        do ii = ixOmin1, ixOmax1
          ambiCoef(ii) = (tmp(ii-1) + tmp(ii) + tmp(ii+1))/3d0  
        enddo
      else
        ambiCoef(ixO^S) = ambiCoef(ixO^S) * tmp(ixO^S)
      endif
      deallocate(tmp)
    endif  

  end subroutine special_ambipolar

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    ! initialize one grid
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
  
    w(ixO^S,1:nw) = 0d0

  end subroutine initonegrid_usr

   subroutine gradient1(w,ixI^L, ixO^L,temp3)
    double precision :: w(ixI^L)
    integer, intent(in) :: ixO^L, ixI^L
    double precision, intent(inout):: temp3(ixI^S)
 
    integer :: ixG^L

    ixG^L=ixO^L;

    if (ixOmin1 .eq. ixImin1) then
      ixGmin1=ixOmin1+1
      temp3(ixOmin1)=0d0
    endif
    if (ixOmax1 .eq. ixImax1) then
      ixGmax1=ixOmax1-1
      temp3(ixOmax1)=0d0
    endif

    !!second order only needs one ghost
    call gradient(w ,ixI^L,ixG^L,1,temp3)

   end subroutine gradient1




  subroutine setLowerBoundary(ixI^L, ixO^L,w,x,time)
    integer, intent(in) :: ixO^L, ixI^L
    double precision, intent(in) :: x(ixI^S,1:ndim),time
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: omega
    complex, parameter :: ic = dcmplx(0,1)
    double complex :: RR(ixO^S), VV(ixO^S), PP(ixO^S), BB(ixO^S), k(ixO^S),wave(ixO^S)!,wave2(ixO^S)
    
    double precision,allocatable :: vA02(:), c02(:),a(:), b(:), c(:),  temp1(:), temp3(:)
    double precision,allocatable, dimension(:) :: pe0, rho0, bx0
    integer :: ixG^L


    omega = 2*dpi/(Period /unit_time)!normalized omega
    !print*, "OMEGA ", omega



    !ixG^L=ixO^L^LADD1;
    ixG^L=ixO^L;
    ixGmax1=ixOmax1+1


    allocate(vA02(ixG^S), c02(ixG^S),a(ixG^S), b(ixG^S), c(ixG^S),  temp1(ixG^S), temp3(ixG^S))
    allocate(rho0(ixG^S), pe0(ixG^S),bx0(ixG^S))

    call set_equi_vars2(x(ixG^S,1), pe0(ixG^S), rho0(ixG^S))
    call set_equi_vars2_b0(x(ixG^S,1), bx0(ixG^S))

    c02(ixG^S) = mhd_gamma * pe0(ixG^S)/rho0(ixG^S)
    vA02(ixG^S) = bx0(ixG^S)**2/rho0(ixG^S)
    !print*, "c02 ", c02(ixG^S) * unit_velocity**2
    !print*, "va02 ", vA02(ixG^S) * unit_velocity**2
    !print*, "rho ", w(ixG^S,rho_)* unit_density
    a(ixG^S) = c02(ixG^S)+vA02(ixG^S)
    
    temp1(ixG^S) = rho0(ixG^S) * a(ixG^S)
    !print*, "a rho ", temp1(ixG^S) * unit_velocity**2 * unit_density
    call gradient1(temp1 ,ixG^L,ixO^L,temp3)
    !print* , "UNIT_length " , unit_length
    !print* , "UNIT_length dx " , dxlevel(1) * unit_length
    !print*, "grad a rho ", temp3(ixO^S)* unit_velocity**2 * unit_density/unit_length
!    print*, dxlevel 
    b(ixO^S) = temp3(ixO^S)/rho0(ixO^S) 
    c(ixO^S) = -omega**2

    !set delta in k
    k(ixO^S) = sqrt(dcmplx(-b(ixO^S)**2-4.0*a(ixO^S)*c(ixO^S),0))
    !print*, "sqrt delta ", k(ixO^S)
    !print*, "a ", a(ixO^S) * unit_velocity**2
    !print*, "b ", b(ixO^S)* (unit_velocity**2/unit_length)
    !print*, "c ", c(ixO^S)/unit_time**2
    k(ixO^S) = (-ic * b(ixO^S) + k(ixO^S))/(2*a(ixO^S))
    !print*, "-B ", (-ic * b(ixO^S))

!    print*, "K ", k(ixO^S)/unit_length

!    print*, "A ", A
!    print*, "c ", sqrt(c02(ixO^S))
    VV(ixO^S) = dcmplx(0,ampl * sqrt(c02(ixO^S)))
!    print*, "VV ", VV(ixO^S) * unit_velocity

    temp1(ixG^S)=rho0(ixG^S)
    call gradient1(temp1 ,ixG^L,ixO^L,temp3)
!    print*, "grad rho  ", temp3(ixO^S)
    temp3(ixO^S) =  temp3(ixO^S)/rho0(ixO^S) 
!    print*, "grad rho/rho  ", temp3(ixO^S)
    RR(ixO^S) =  rho0(ixO^S)* VV(ixO^S) * (k(ixO^S) + ic * temp3(ixO^S))/omega

!    print*, "rho0 ", w(ixO^S,rho_)* unit_density
!    print*, "RR ", RR(ixO^S) * unit_density

    temp1(ixG^S)=pe0(ixG^S)
    call gradient1(temp1 ,ixG^L,ixO^L,temp3)
!    print*, "grad p  ", temp3(ixO^S)
    temp3(ixO^S) =  temp3(ixO^S)/pe0(ixO^S)
    PP(ixO^S) =  pe0(ixO^S)* VV(ixO^S) * (k(ixO^S) * mhd_gamma + ic * temp3(ixO^S))/omega
!    print*, "p0 ", w(ixO^S,p_)* unit_pressure
!    print*, "PP ", PP(ixO^S) * unit_pressure

    temp1(ixG^S)=bx0(ixG^S)
    call gradient1(temp1 ,ixI^L,ixG^L,temp3)
!    print*, "SETbc grad bx  ", temp3(ixO^S)
    temp3(ixO^S) =  temp3(ixO^S)/bx0(ixO^S)
    BB(ixO^S) =  bx0(ixO^S)* VV(ixO^S) * (k(ixO^S)  + ic * temp3(ixO^S))/omega

!    print*, "bx ", w(ixO^S,mag(3))* unit_magneticfield
!    print*, "BB ", BB(ixO^S) * unit_magneticfield
!
!    print*, "x ", x(ixO^S,1)
!    print*, "k x ", k(ixO^S)*x(ixO^S,1)
!    print*, "omega * time ", omega * time
!    print*, "argExp ",(ic *  (omega * time - k(ixO^S)*(x(ixO^S,1)-xprobmin1)) )
    wave(ixO^S) = exp(ic *  (omega * time - k(ixO^S)*(x(ixO^S,1) -xprobmin1)) )
    !print*, "WAVE ", wave(ixO^S)
    

!    print*, "SETbc c02 ",  c02(1:3)
    deallocate(vA02, c02,a , b, c,  temp1, temp3)

    !print*, "EQUI IN GHoSTS rho ", block%equi_vars(ixO^S,equi_rho_c0_)
    !print*, "EQUI IN GHoSTS pe ", block%equi_vars(ixO^S,equi_pe_c0_)
    !print*, "EQUI IN GHoSTS B0x ", block%B0(ixO^S,3,0)

    w(ixO^S,rho_)=real(wave(ixO^S)*RR(ixO^S))
    w(ixO^S,e_)=real(wave(ixO^S)*PP(ixO^S))
    w(ixO^S,mag(3))=real(wave(ixO^S)*BB(ixO^S))
    w(ixO^S,mom(1))= real(wave(ixO^S)*VV(ixO^S))
    w(ixO^S,mom(2:3))= 0d0
    w(ixO^S,mag(1:2))= 0d0
!    print*, "SETbc for  ", x(ixG^S,1), " at time ", time
!    print*, "SETbc vals ",  w(1:3,mag(3))
!    print*, "SETbc bx0 ",  bx0(1:3)
!    print*, "SETbc k ",  k(1:3)
!    print*, "SETbc wave ",  wave(1:3)
!    print*, "SETbc VV ",  VV(1:3)
!    print*, "SETbc BB ",  BB(1:3)

    deallocate(pe0, rho0, bx0)
!    print*, "SET"
!    print*, "x ", x(ixO^S,1)
!    print*, w(ixO^S,mom(1))* unit_velocity
!!    print*, w(ixO^S,rho_)* unit_density
!!    print*, w(ixO^S,p_)* unit_pressure
!!    print*, w(ixO^S,mag(3))* unit_magneticfield
!    print*, real(wave(ixO^S)*RR(ixO^S))* unit_density
!    print*, real(wave(ixO^S)*PP(ixO^S))* unit_pressure
!    print*, real(wave(ixO^S)*BB(ixO^S))* unit_magneticfield
!    print*, "END SET"

     ! print*, "INDICES ", ixO^L, " wsize ", size(w,1), size(w,2)

      call mhd_to_conserved(ixI^L,ixO^L,w,x)

  end subroutine setLowerBoundary


  subroutine gravity(ixI^L,ixO^L,wCT,x,gravity_field)
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(in)    :: wCT(ixI^S,1:nw)
    double precision, intent(out)   :: gravity_field(ixI^S,ndim)


    gravity_field=usr_grav

  end subroutine gravity

  
  subroutine special_process_filter(igrid,level,ixI^L,ixO^L, qt,w,x)
    use mod_global_parameters
    integer, intent(in) :: igrid, level, ixI^L, ixO^L
    double precision, intent(in) :: qt,  x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    !integer, save :: lastit = -1


    if(mod(it,niter_filter).eq.0) then ! .and. it .ne. lastit) then

      call filter_x(ixI^L,ixO^L,w)
      !lastit = it
   endif

end  subroutine special_process_filter

  !as process
    subroutine special_process(igrid,level,ixI^L,ixO^L,qt, w,x)
      use mod_global_parameters
      integer, intent(in)             :: ixI^L,ixO^L,igrid, level


  !as int_bc
  !  subroutine special_process(level,qt,ixI^L,ixO^L,w,x)
  !    use mod_global_parameters
  !    integer, intent(in)             :: ixI^L,ixO^L,level

      double precision, intent(in)    :: qt
      double precision, intent(inout) :: w(ixI^S,1:nw)
      double precision, intent(in)    :: x(ixI^S,1:ndim)


    double precision,allocatable, dimension(:) :: pe0, rho0, bx0,c0
    double precision,allocatable, dimension(:) :: mask
    double precision,allocatable, dimension(:) :: tmp

    integer :: i, ixG^L

    !double precision,parameter :: xLMin = 1.65
    !double precision :: xL


    logical :: found

!!!!!METHOD1!!!!!
!    xL = xLMin +  rand(24) * 0.35
!    !if(mod(it,100) .eq. 0) then
!    !  xL = 1.8
!    !else
!    !  xL = xLMin + (10 - mod(it,11)) * 0.01 
!    !endif
!!!!!METHOD1end!!!!!

    !print*, "XL", xL
    ixG^L=ixO^L;

!!!!!METHOD2!!!!!
    !xL=1.7
    !xL=1.75
!!!!!METHOD2!!!!!

    if(x(ixOmax1,1) < xL) then
      found = .false.
    else

      if(x(ixOmin1,1) < xL) then
        found = .false.
        i=ixOmin1+1
        do while (.not. found .and. i .le. ixOmax1)
          if(x(i,1) .ge. xL) then
            ixGmin1=i
            found = .true.
          endif 
          i=i+1
         enddo
      else
        found = .true.
      endif 
    endif
 

   if(found) then

      allocate(rho0(ixG^S), pe0(ixG^S),bx0(ixG^S),c0(ixG^S))
      call set_equi_vars2(x(ixG^S,1), pe0(ixG^S), rho0(ixG^S))
      call set_equi_vars2_b0(x(ixG^S,1), bx0(ixG^S))
      c0(ixG^S) = sqrt(bx0(ixG^S)/rho0(ixG^S)) + sqrt(mhd_gamma * pe0(ixG^S)/rho0(ixG^S))

      !!METHOD1
!      w(ixG^S,rho_)=rho0(ixG^S) 
!      w(ixG^S,p_)=pe0(ixG^S)
!      w(ixG^S,mag(3))=bx0(ixG^S) 
!      w(ixG^S,mom(1))= 0.0
      !!METHOD1end

      !!METHOD2
      allocate(mask(ixG^S))
      !mask(ixG^S) = 0.5*(1d0-tanh((x(ixG^S,1) - 1.9)/0.2))
      do i = ixGmin1, ixGmax1
        mask(i) = (1.0 - (x(i,1)-x(ixGmin1,1))/(xprobmax1-x(ixGmin1,1)))**2
        !mask(i) = (1.0 - c0(ixG^S)/c0(ixGmin1))**2
        !mask(i) = (1.0 - c0(i)/c0(ixGmin1))
        !mask(i) = (1.0 - (xprobmax1-x(i,1))/(xprobmax1-xL))
      enddo
        !print*, "mask ", mask(ixG^S)
        !print*, "old mask ", (1.0 - (xprobmax1-x(ixG^S,1))/(xprobmax1-x(ixGmin1,1)))**2


      !!mnacha
      !mask(ixG^S) = 1.0 - sigmaPml * mask(ixG^S) * c0(ixGmin1) * dt/dx(level,1)  !!!  NOW
      mask(ixG^S) = exp(- sigmaPml * mask(ixG^S) * c0(ixG^S) * dt/dx(level,1) )




      !for xL=0.7
      !mask(ixG^S) = 1.0 - 70 * mask(ixG^S) * c0Ref * dt/dx(level,1)
      !print*, "M ", mask(ixG^S)
      allocate(tmp(ixG^S))
      !mask(ixG^S) = 0.5*(1d0-tanh((x(ixG^S,1) - 1.7)/0.01))
      !print*, mask(ixG^S)
      call mhd_to_primitive(ixI^L,ixG^L,w,x)
      tmp(ixG^S)= w(ixG^S,rho_) - rho0(ixG^S)
      w(ixG^S,rho_)=rho0(ixG^S) + tmp(ixG^S)*mask(ixG^S)
      tmp(ixG^S)= w(ixG^S,e_) - pe0(ixG^S)
      w(ixG^S,e_)=pe0(ixG^S)+ tmp(ixG^S)*mask(ixG^S)
      tmp(ixG^S)= w(ixG^S,mag(3)) - bx0(ixG^S)
      w(ixG^S,mag(3))=bx0(ixG^S) + tmp(ixG^S)*mask(ixG^S)
      tmp(ixG^S)= w(ixG^S,mom(1)) 
      w(ixG^S,mom(1))= tmp(ixG^S)*mask(ixG^S)
      deallocate(mask)
      deallocate(tmp)
      !!METHOD2


      w(ixG^S,mag(1:2))= 0.0
      w(ixG^S,mom(2:3))= 0.0


      call mhd_to_conserved(ixI^L,ixG^L,w,x)
      deallocate(rho0, pe0,bx0)
    endif

 end subroutine special_process



    !3 ghostcells needed   
  subroutine filter_x(ixI^L,ixO^L,w)
    double precision, intent(inout)    :: w(ixI^S,1:nw)
    integer, intent(in) :: ixI^L, ixO^L
 
    double precision, DIMENSION(7), PARAMETER :: d=(/  -0.015625000d0,0.093750000d0, -0.234375000d0, &
         0.312500000d0,-0.234375000d0, 0.093750000d0, -0.015625000d0  /) ! Parameters for the filter
    
    double precision    :: dum(ixI^S,1:nw)
    integer :: ii
    

    !print*, "FILTER ", ixI^L, ixO^L

     dum=d(4)*w
     do ii=1,3
        dum=dum+d(ii+4)*(cshift(w,ii) + cshift(w,-ii))
     enddo
  
  
     !dum(ixImin1:ixOmin1-1,:) = 0.0d0
     !dum(ixOmax1+1:ixImax1,:) = 0.0d0

    w(ixO^S,1:nw) = w(ixO^S,1:nw) - dum(ixO^S,1:nw)

  contains

  !====================================================================
  !> shift: shift a one dimensional array by k elements
  !>
  !> If k is positive, shift to the right k positions
  !> If k is negative, shift to the left
  !> 
  !> Cshift (Fortran intrinsic fuction) was much slower when I (angelv)
  !> tried it (March 2012), so we leave this one instead
  !==================================================================== 
  PURE FUNCTION shift(array,k)
    double precision, DIMENSION(:,:), INTENT(IN) :: array
    INTEGER, INTENT(IN) :: k

    double precision, DIMENSION(SIZE(array,1), size(array,2)) :: shift
    INTEGER :: n

    !__________________________________________

    n=SIZE(array,1)

    !shift=0.0_fp !TODO needed?
    IF (k >  0) THEN
       shift(k+1:n,:)=array(1:n-k,:)
       shift(1:k,:)=array(n-k+1:n,:)
    ElSE
       shift(1:n+k,:)=array(-k+1:n,:)
       shift(n+k+1:n,:)=array(1:-k,:)
    ENDIF
  END FUNCTION shift
  end subroutine filter_x
!



  subroutine init_params_usr()

    double precision :: dx
    double precision, allocatable,target :: xx(:)
    integer :: nx,k

    usr_grav=-2.74d2*unit_length/unit_velocity**2 ! solar gravity

     !!from the 2fluid code
     !! I calculate this eta as mean(cjbb * density**2) and made it adim 
!    >>> eta = 3.176883735537642e-09
!    >>> L=1000000.0000000000
!    >>> t=110.06703398489996
!    >>> d=1.6726217770000000E-007
!    >>> u=t*d
!    >>> u
!    1.8410051797294276e-05
!    >>> eta/u
!    0.00017256245503907534


    call dump_units()


    if(mype .eq. 0) then
      dx = (xprobmax1-xprobmin1)/(domain_nx1-1)
      nx = domain_nx1 + 2*nghostcells
      allocate(xx(nx))
  
      do k=0,nx-1
        xx(k+1) = xprobmin1-(nghostcells - 0.5d0)*dx + k*dx
      enddo
  
      call init_equi_vars(xx)
      deallocate(xx)
    endif

    call MPI_BARRIER(icomm,ierrmpi)

  end subroutine init_params_usr



  subroutine specialbound_usr(qt,ixI^L,ixO^L,iB,w,x)
    ! special boundary types, user defined
    integer, intent(in) :: ixO^L, iB, ixI^L
    double precision, intent(in) :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)


    select case(iB)
    ! special left boundary
    case(1)
      call setLowerBoundary(ixI^L, ixO^L,w,x,qt)
    end select

  end subroutine specialbound_usr


  function rargsort(a) result(b)
    ! Returns the indices that would sort an array.
    !
    ! Arguments
    ! ---------
    !
    double precision, intent(in):: a(:)   ! array of numbers
    integer :: b(size(a))         ! indices into the array 'a' that sort it
    !
    ! Example
    ! -------
    !
    ! rargsort([4.1_dp, 2.1_dp, 2.05_dp, -1.5_dp, 4.2_dp]) ! Returns [4, 3, 2, 1, 5]
    
    integer :: N                           ! number of numbers/vectors
    integer :: i,imin                      ! indices: i, i of smallest
    integer :: temp1                       ! temporary
    double precision:: temp2
    double precision:: a2(size(a))
    a2 = a
    N=size(a)
    do i = 1, N
        b(i) = i
    end do
    do i = 1, N-1
        ! find ith smallest in 'a'
        imin = minloc(a2(i:),1) + i - 1
        ! swap to position i in 'a' and 'b', if not already there
        if (imin /= i) then
            temp2 = a2(i); a2(i) = a2(imin); a2(imin) = temp2
            temp1 = b(i); b(i) = b(imin); b(imin) = temp1
        end if
    end do
  end function




  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
    use mod_mhd_phys
    ! this subroutine can be used in convert, to add auxiliary variables to the
    ! converted output file, for further analysis using tecplot, paraview, ....
    ! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
    !
    ! the array normconv can be filled in the (nw+1:nw+nwauxio) range with 
    ! corresponding normalization values (default value 1)
    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)
    
    w(ixO^S,nw+1)=1d0
    call multiplyAmbiCoef(ixI^L,ixO^L,w(ixI^S,nw+1),w,x)

  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
    character(len=*) :: varnames
    varnames='ambi'

  end subroutine specialvarnames_output


end module mod_usr
