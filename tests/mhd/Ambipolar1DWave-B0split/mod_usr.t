module mod_usr
  use mod_mhd


  implicit none


  !!user defined
  double precision:: Period=5
  double precision:: ampl=1d-3
  
  logical :: maskAmbi = .true.
  double precision :: xLambi = 1.65
  double precision :: wLambi = 0.01
  logical :: ambi_mask_smooth  = .true.
  integer, parameter :: MASK_DISC = 1
  integer, parameter :: MASK_TANH = 2
  integer :: ambi_mask_method = MASK_DISC


  character(len=*), parameter :: VALCfilename="valc.txt"
  character(len=*), parameter :: equi_filename="myEqui.txt"
  integer, parameter :: equi_zz=1, equi_pe=2, equi_rho=3, equi_bx=4 !indices in equi generated file
  double precision:: usr_grav


  type arrayptr
    double precision,dimension(:), pointer :: p
    character(len=30) ::  namevar
  end type arrayptr



contains

  !> Read this module's parameters from a file
  subroutine usr_params_read(files)
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /usr_list/ Period,ampl,maskAmbi,xLambi,&
              wLambi,ambi_mask_method,ambi_mask_smooth

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




  !!this has been obtained from the two-fluid model
  !by adding neutrals and charges
  subroutine set_equi_vars_valc(xnorm, pe, rho, bx0)
    use mod_interpolation
    double precision, intent(in) :: xnorm(:)
    double precision, dimension(:) :: pe,rho,bx0
 
    double precision, allocatable  :: xx(:)
    integer :: nrows, ncols
    double precision, allocatable, dimension(:,:) :: valcData
    !the indices of columns
    integer, parameter :: zz_ = 1
    integer, parameter :: te_ = 4
    integer, parameter :: nn_ = 6
    integer, parameter :: ne_ = 7
    integer, parameter :: rho_ = 10
    !!constants see mod_constants
    double precision, parameter ::  MH=1.67339e-27
    double precision, parameter :: BK=1.380662e-23
    double precision, parameter ::  G=2.74d2  !!!!use exactly the same value as in usr_grav, but unnormalized
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

    eqP = 1+int((nz-1)/2)
    bx0(:) = sqrt(2 * pe(eqP) )
  
    deallocate(valcData,zzValc,ind, tempValc)
    deallocate(te)
  end subroutine set_equi_vars_valc

  subroutine set_equi_vars2(xnorm, pe, rho)
    use mod_interpolation
    !use interp1D
    double precision, intent(in) :: xnorm(:)
    double precision, dimension(:) :: pe,rho
 
    integer :: nrows, ncols, nz
    !!!we already know there are 10 columns
    double precision, allocatable, dimension(:,:) :: valcData
    nrows = get_number_lines(equi_filename)
    ncols = get_number_columns(equi_filename)

    nz = size(xnorm)

    allocate(valcData(nrows,ncols))
    call read_formatted_file(equi_filename,nrows,ncols,valcData)

    call interp_cubic_spline(valcData(:,equi_zz),valcData(:,equi_rho),nrows,xnorm,rho,nz)
    call interp_cubic_spline(valcData(:,equi_zz),valcData(:,equi_pe),nrows,xnorm,pe,nz)
 
   deallocate(valcData)
  end subroutine set_equi_vars2

  subroutine set_equi_vars2_b0(xnorm, bx0)
    use mod_interpolation
    !use interp1D
    double precision, intent(in) :: xnorm(:)
    double precision, dimension(:) :: bx0
 
    integer :: nrows, ncols, nz
    !!!we already know there are 10 columns
    double precision, allocatable, dimension(:,:) :: valcData
    nrows = get_number_lines(equi_filename)
    ncols = get_number_columns(equi_filename)

    nz = size(xnorm)

    allocate(valcData(nrows,ncols))
    call read_formatted_file(equi_filename,nrows,ncols,valcData)

    call interp_cubic_spline(valcData(:,equi_zz),valcData(:,equi_bx),nrows,xnorm,bx0,nz)
 
   deallocate(valcData)
  end subroutine set_equi_vars2_b0

  subroutine init_equi_vars(xnorm)
    double precision, allocatable, target, intent(in) :: xnorm(:)
    double precision, allocatable, dimension(:),target :: pe,rho,bx0
    type(arrayptr), dimension(4) :: write_vars
    integer :: nz 
 
    nz= size(xnorm)

    allocate(pe(nz),rho(nz),bx0(nz))
    call set_equi_vars_valc(xnorm,pe,rho,bx0)

    write_vars(equi_zz)%namevar = "z"
    write_vars(equi_pe)%namevar  = "pe"
    write_vars(equi_rho)%namevar  = "rho"
    write_vars(equi_bx)%namevar  = "bx0"

    write_vars(equi_zz)%p =>xnorm
    write_vars(equi_pe)%p =>pe
    write_vars(equi_rho)%p =>rho
    write_vars(equi_bx)%p =>bx0

    call  write_formatted_file(equi_filename,write_vars)
  
    deallocate(pe,rho,bx0)
  end subroutine init_equi_vars




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
    usr_special_bc    => specialbound_usr
    usr_gravity         => gravity
    usr_set_B0              => specialset_B0
    usr_set_J0              => specialset_J0

    usr_set_parameters  => init_params_usr
    !usr_process_grid => special_process_filter

    if (maskAmbi) then
      usr_mask_ambipolar => special_ambipolar
    endif

    usr_aux_output    => specialvar_output
    usr_add_aux_names => specialvarnames_output
    call set_coordinate_system("Cartesian_1.75D")

    call mhd_activate()

  end subroutine usr_init

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

    integer :: ii

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
    endif
  end subroutine special_ambipolar

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    ! initialize one grid
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    call set_equi_vars2(x(ixO^S,1), w(ixO^S,p_), w(ixO^S,rho_))
    call mhd_to_conserved(ixI^L,ixO^L,w,x)
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


  !>symm boundary condition for the perturbation only. The equilibrium is stratitified, not uniform.
  subroutine setUpperBoundary(ixI^L, ixO^L,w,x,time)
    integer, intent(in) :: ixO^L, ixI^L
    double precision, intent(in) :: x(ixI^S,1:ndim),time
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision,allocatable, dimension(:) :: pe0, rho0
    integer :: ixG^L

    ixGmin1 = ixOmin1-nghostcells
    ixGmax1 = ixOmin1-1
    call mhd_to_primitive(ixI^L,ixG^L,w,x)

    ixGmax1 = ixOmin1+nghostcells-1


    allocate(rho0(ixG^S), pe0(ixG^S))
    
    call set_equi_vars2(x(ixG^S,1), pe0(ixG^S), rho0(ixG^S))

     w(ixO^S,rho_) =   (w(ixOmin^D-1:ixOmin^D-nghostcells:-1^D%ixO^S,rho_) - rho0(ixOmin^D-1:ixOmin^D-nghostcells:-1^D%ixO^S)) + rho0(ixO^S)
     w(ixO^S,p_) = (w(ixOmin^D-1:ixOmin^D-nghostcells:-1^D%ixO^S,p_)- pe0(ixOmin^D-1:ixOmin^D-nghostcells:-1^D%ixO^S)) + pe0(ixO^S)


    deallocate(pe0, rho0)
    call mhd_to_conserved(ixI^L,ixG^L,w,x)
  end subroutine setUpperBoundary


  subroutine setLowerBoundary(ixI^L, ixO^L,w,x,time)
    integer, intent(in) :: ixO^L, ixI^L
    double precision, intent(in) :: x(ixI^S,1:ndim),time
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: omega
    double complex, parameter :: ic = dcmplx(0,1)
    double complex :: RR(ixO^S), VV(ixO^S), PP(ixO^S), BB(ixO^S), k(ixO^S),wave(ixO^S)
    
    double precision,allocatable :: vA02(:), c02(:),a(:), b(:), c(:),  temp1(:), temp3(:)
    double precision,allocatable, dimension(:) :: pe0, rho0, bx0
    integer :: ixG^L


    omega = 2*dpi/(Period /unit_time)!normalized omega
    ixG^L=ixO^L;
    ixGmax1=ixOmax1+1


    allocate(vA02(ixG^S), c02(ixG^S),a(ixG^S), b(ixG^S), c(ixG^S),  temp1(ixG^S), temp3(ixG^S))
    allocate(rho0(ixG^S), pe0(ixG^S),bx0(ixG^S))

    call set_equi_vars2(x(ixG^S,1), pe0(ixG^S), rho0(ixG^S))
    call set_equi_vars2_b0(x(ixG^S,1), bx0(ixG^S))

    c02(ixG^S) = mhd_gamma * pe0(ixG^S)/rho0(ixG^S)
    vA02(ixG^S) = bx0(ixG^S)**2/rho0(ixG^S)
    a(ixG^S) = c02(ixG^S)+vA02(ixG^S)
    
    temp1(ixG^S) = rho0(ixG^S) * a(ixG^S)
    call gradient1(temp1 ,ixG^L,ixO^L,temp3)
    b(ixO^S) = temp3(ixO^S)/rho0(ixO^S) 
    c(ixO^S) = -omega**2

    !set delta in k
    k(ixO^S) = sqrt(dcmplx(-b(ixO^S)**2-4.0*a(ixO^S)*c(ixO^S),0))
    k(ixO^S) = (-ic * b(ixO^S) + k(ixO^S))/(2*a(ixO^S))
    VV(ixO^S) = dcmplx(0,ampl * sqrt(c02(ixO^S)))

    temp1(ixG^S)=rho0(ixG^S)
    call gradient1(temp1 ,ixG^L,ixO^L,temp3)
    temp3(ixO^S) =  temp3(ixO^S)/rho0(ixO^S) 
    RR(ixO^S) =  rho0(ixO^S)* VV(ixO^S) * (k(ixO^S) + ic * temp3(ixO^S))/omega


    temp1(ixG^S)=pe0(ixG^S)
    call gradient1(temp1 ,ixG^L,ixO^L,temp3)
    temp3(ixO^S) =  temp3(ixO^S)/pe0(ixO^S)
    PP(ixO^S) =  pe0(ixO^S)* VV(ixO^S) * (k(ixO^S) * mhd_gamma + ic * temp3(ixO^S))/omega

    temp1(ixG^S)=bx0(ixG^S)
    call gradient1(temp1 ,ixI^L,ixG^L,temp3)
    temp3(ixO^S) =  temp3(ixO^S)/bx0(ixO^S)
    BB(ixO^S) =  bx0(ixO^S)* VV(ixO^S) * (k(ixO^S)  + ic * temp3(ixO^S))/omega

    wave(ixO^S) = exp(ic *  (omega * time - k(ixO^S)*(x(ixO^S,1) -xprobmin1)) )
    

    deallocate(vA02, c02,a , b, c,  temp1, temp3)

    w(ixO^S,rho_)=rho0(ixO^S) + real(wave(ixO^S)*RR(ixO^S))
    w(ixO^S,p_)=pe0(ixO^S) + real(wave(ixO^S)*PP(ixO^S))
    w(ixO^S,mag(3))=real(wave(ixO^S)*BB(ixO^S))
    w(ixO^S,mom(1))= real(wave(ixO^S)*VV(ixO^S))

    deallocate(pe0, rho0, bx0)

      call mhd_to_conserved(ixI^L,ixO^L,w,x)

  end subroutine setLowerBoundary


  subroutine gravity(ixI^L,ixO^L,wCT,x,gravity_field)
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(in)    :: wCT(ixI^S,1:nw)
    double precision, intent(out)   :: gravity_field(ixI^S,ndim)


    gravity_field=usr_grav

  end subroutine gravity

  


  subroutine init_params_usr()

    double precision :: dx
    double precision, allocatable,target :: xx(:)
    integer :: nx,k

    usr_grav=-2.74d2*unit_length/unit_velocity**2 ! solar gravity

    if(mype .eq. 0) then
      dx = (xprobmax1-xprobmin1)/(domain_nx1)
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
    ! special right boundary
    case(2)
      call setUpperBoundary(ixI^L, ixO^L,w,x,qt)

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
