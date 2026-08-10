!> Potential-field initial condition from a Python V1 data-driven boundary frame.
!>
!> The boundary frame is written by tools/python/notebooks/data_driven_pipeline.ipynb
!> with layout: snapshot_time, nx, ny, dx, dy, Bx, By, Bz.
module mod_usr
  use mod_mhd
  use mod_lfff
  implicit none

  double precision :: lalpha,llift
  character(len=256) :: boundary_filename
  character(len=16) :: potential_field_method
  character(len=16) :: fft_top_boundary
  character(len=16) :: lfff_flux_treatment
  double precision :: lfff_max_flux_imbalance
  integer :: fft_padding_factor
  logical, save :: firstusrglobaldata=.true.

contains

  subroutine usr_init()
    use mod_usr_methods

    ! Standard coronal normalization used throughout the AMRVAC data-driven
    ! workflow. MHD derives the magnetic-field unit from these three units.
    unit_length        = 1.d9 ! cm
    unit_temperature   = 1.d6 ! K
    unit_numberdensity = 1.d9 ! cm^-3

    call set_coordinate_system('Cartesian_3D')

    usr_set_parameters => initglobaldata_usr
    usr_init_one_grid  => initonegrid_usr
    usr_improve_initial_condition => improve_initial_condition_usr

    call mhd_activate()
  end subroutine usr_init

  subroutine usr_params_read(files)
    use mod_global_parameters

    character(len=*), intent(in) :: files(:)
    integer :: n

    namelist /usr_list/ boundary_filename,lalpha,llift,potential_field_method,&
       fft_padding_factor,fft_top_boundary,&
       lfff_flux_treatment,lfff_max_flux_imbalance

    boundary_filename = 'boundary_single/B_0001.dat'
    lalpha = 0.d0
    llift = 0.d0
    potential_field_method = 'fft'
    fft_padding_factor = 2
    fft_top_boundary = 'open'
    lfff_flux_treatment = 'strict'
    lfff_max_flux_imbalance = 0.1d0
    do n=1,size(files)
      open(unitpar,file=trim(files(n)),status='old')
      read(unitpar,usr_list,end=111)
111   close(unitpar)
    end do
  end subroutine usr_params_read

  subroutine initglobaldata_usr()
    use mod_global_parameters

    call usr_params_read(par_files)

    if(refine_max_level/=1) &
       call mpistop('PotentialField requires refine_max_level=1')

    potential_field_method = lowercase(trim(adjustl(potential_field_method)))
    select case(trim(potential_field_method))
    case('green')
      continue
    case('fft')
      if(llift/=0.d0) call mpistop('FFT potential field requires llift=0')
      fft_top_boundary = lowercase(trim(adjustl(fft_top_boundary)))
      if(trim(fft_top_boundary)/='open' .and. &
         trim(fft_top_boundary)/='closed') &
         call mpistop("fft_top_boundary must be 'open' or 'closed'")
      lfff_flux_treatment = lowercase(trim(adjustl(lfff_flux_treatment)))
      if(trim(lfff_flux_treatment)/='strict' .and. &
         trim(lfff_flux_treatment)/='subtract_mean') &
         call mpistop("lfff_flux_treatment must be 'strict' or 'subtract_mean'")
      if(lfff_max_flux_imbalance<0.d0 .or. &
         lfff_max_flux_imbalance>1.d0) &
         call mpistop('lfff_max_flux_imbalance must be between zero and one')
    case default
      call mpistop("potential_field_method must be 'green' or 'fft'")
    end select
    if(fft_padding_factor<1) call mpistop('fft_padding_factor must be at least one')

    if(firstusrglobaldata) then
      call init_b_fff_data_driven_boundary(trim(boundary_filename),&
         unit_length,unit_magneticfield)
      firstusrglobaldata = .false.
    end if
  end subroutine initglobaldata_usr

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    use mod_global_parameters

    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: Bf(ixI^S,1:ndir)
    logical, save :: first=.true.

    if(mype==0 .and. first) then
      write(*,*) 'Initializing potential-field grids from data-driven boundary'
      first = .false.
    end if

    select case(trim(potential_field_method))
    case('green')
      call calc_lin_fff(ixI^L,ixO^L,Bf,x,lalpha,llift)
      w(ixO^S,mag(:)) = Bf(ixO^S,1:3)
    case('fft')
      ! The global FFT solver fills the magnetic field after the AMR tree is
      ! complete. Initialize only the local non-magnetic state here.
      w(ixO^S,mag(:)) = zero
    end select
    w(ixO^S,mom(:)) = zero
    w(ixO^S,rho_) = one
  end subroutine initonegrid_usr

  subroutine improve_initial_condition_usr()
    use mod_global_parameters
    use mod_ghostcells_update, only: getbc

    if(trim(potential_field_method)/='fft') return

    if(mype==0) write(*,*) 'Initializing FFT potential/LFFF magnetic field'
    ! The observed magnetogram is located at z=0 on the first lower ghost-cell
    ! center, half a finest-cell spacing below the physical domain face.
    call extrapolate_potential_fft(mag(:),fft_padding_factor,&
       0.5d0*dx(3,refine_max_level),lalpha,fft_top_boundary,&
       lfff_flux_treatment,lfff_max_flux_imbalance)
    call getbc(global_time,0.d0,ps,iwstart,nwgc)
  end subroutine improve_initial_condition_usr

  pure function lowercase(input) result(output)
    character(len=*), intent(in) :: input
    character(len=len(input)) :: output
    integer :: i,code

    output=input
    do i=1,len(input)
      code=iachar(input(i:i))
      if(code>=iachar('A') .and. code<=iachar('Z')) &
         output(i:i)=achar(code+iachar('a')-iachar('A'))
    end do
  end function lowercase

end module mod_usr
