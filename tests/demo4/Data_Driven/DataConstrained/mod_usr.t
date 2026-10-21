!> Fixed-bottom data-constrained MHD initialized from a magnetic snapshot.
!>
!> The restart snapshot supplies only the magnetic field. Density, velocity,
!> pressure, and (when present) total energy are rebuilt from a user-selected
!> atmosphere when firstprocess=.true.  A normal checkpoint restart must leave
!> firstprocess and the reset flags false, so the evolved plasma is preserved.
module mod_usr
  use mod_mhd
  use mod_eos, only: eos
  implicit none

  integer, parameter :: max_line_len=1024
  double precision, allocatable, save :: Bx0(:,:),By0(:,:),Bz0(:,:)
  double precision, allocatable, save :: rho_bottom(:,:,:),p_bottom(:,:,:)
  double precision, allocatable, save :: atmosphere_z(:),atmosphere_rho(:,:)
  ! atmosphere_rho(:,1:3) stores rho, p, and T.

  character(len=256) :: boundary_filename
  character(len=256) :: relaxed_atmosphere_file
  character(len=32)  :: physics_model,atmosphere_model,atmosphere_source
  character(len=32)  :: temperature_curve
  integer, save :: nx_boundary=0,ny_boundary=0,n_atmosphere=0
  integer, save :: pad_boundary1=0,pad_boundary2=0
  double precision :: coronal_temperature_k
  double precision :: rho_reference_height_cm
  double precision :: rho_reference_numberdensity_cm3
  double precision :: heating_amplitude_cgs,heating_scale_height_cm
  double precision, save :: solar_gravity,solar_radius,heating_amplitude
  double precision, save :: heating_scale_height
  logical, save :: is_thermodynamic=.false.

contains

  subroutine usr_init()
    use mod_usr_methods

    unit_length        = 1.d9 ! cm
    unit_temperature   = 1.d6 ! K
    unit_numberdensity = 1.d9 ! cm^-3

    call set_coordinate_system('Cartesian_3D')

    usr_set_parameters => initglobaldata_usr
    usr_init_one_grid  => initonegrid_usr
    usr_transform_w    => transform_w_from_magnetic_snapshot
    usr_special_bc     => specialbound_usr
    usr_refine_grid    => special_refine_grid
    usr_gravity        => gravity_usr
    usr_source         => heating_source_usr

    call mhd_activate()
  end subroutine usr_init

  subroutine usr_params_read(files)
    use mod_comm_lib, only: mpistop
    use mod_global_parameters

    character(len=*), intent(in) :: files(:)
    integer :: n,ios

    namelist /usr_list/ boundary_filename,physics_model,atmosphere_model, &
      atmosphere_source,temperature_curve,coronal_temperature_k, &
      rho_reference_height_cm,rho_reference_numberdensity_cm3, &
      relaxed_atmosphere_file,heating_amplitude_cgs, &
      heating_scale_height_cm

    boundary_filename = 'boundary_single/B_0001.dat'
    physics_model = 'zero_beta'
    atmosphere_model = 'uniform'
    atmosphere_source = 'hydrostatic'
    temperature_curve = 'AL-C7'
    coronal_temperature_k = 1.d6
    rho_reference_height_cm = 0.d0
    rho_reference_numberdensity_cm3 = 1.d9
    relaxed_atmosphere_file = ''
    heating_amplitude_cgs = 1.d-4
    heating_scale_height_cm = 5.d9

    do n=1,size(files)
      open(unitpar,file=trim(files(n)),status='old',iostat=ios)
      if(ios/=0) cycle
      read(unitpar,usr_list,iostat=ios)
      close(unitpar)
      if(ios>0) call mpistop('failed to read usr_list')
    end do
  end subroutine usr_params_read

  subroutine initglobaldata_usr()
    use mod_comm_lib, only: mpistop
    use mod_global_parameters

    call usr_params_read(par_files)
    if(restart_from_file==undefined) &
      call mpistop('DataConstrained must restart from a magnetic snapshot')

    solar_gravity = -2.74d4*unit_length/unit_velocity**2
    solar_radius = 6.957d10/unit_length
    heating_amplitude = heating_amplitude_cgs/unit_pressure*unit_time
    heating_scale_height = heating_scale_height_cm/unit_length

    call validate_configuration()
    call init_b_data_constrained_boundary(trim(boundary_filename),unit_magneticfield)
    call init_atmosphere()
    call init_bottom_thermodynamic_boundary()
    if(firstprocess) call write_initial_atmosphere()
  end subroutine initglobaldata_usr

  subroutine validate_configuration()
    use mod_comm_lib, only: mpistop
    use mod_global_parameters

    if(coronal_temperature_k<=0.d0) &
      call mpistop('coronal_temperature_k must be positive')
    if(rho_reference_numberdensity_cm3<=0.d0) &
      call mpistop('rho_reference_numberdensity_cm3 must be positive')
    if(heating_scale_height_cm<=0.d0) &
      call mpistop('heating_scale_height_cm must be positive')
    if(heating_amplitude_cgs<0.d0) &
      call mpistop('heating_amplitude_cgs must be non-negative')

    select case(trim(physics_model))
    case('zero_beta')
      if(mhd_energy .or. mhd_adiab/=0.d0 .or. mhd_gravity) &
        call mpistop('zero_beta par must disable energy, pressure, and gravity')
      if(trim(atmosphere_source)=='relaxed_table' .and. &
         len_trim(relaxed_atmosphere_file)==0) &
        call mpistop('relaxed_table requires relaxed_atmosphere_file')
      if(mype==0 .and. (trim(atmosphere_model)/='uniform' .or. &
         trim(atmosphere_source)=='relaxed_table')) then
        write(*,*) 'WARNING: zero_beta atmosphere is only a density prescription.'
        write(*,*) 'Pressure and gravity are off; the atmosphere is not in HSE.'
      end if
    case('isothermal')
      if(mhd_energy .or. mhd_adiab<=0.d0 .or. .not.mhd_gravity) &
        call mpistop('invalid energy, mhd_adiab, or gravity for isothermal')
      if(trim(atmosphere_model)/='corona' .or. trim(atmosphere_source)/='hydrostatic') &
        call mpistop('isothermal supports only corona + hydrostatic')
      if(abs(mhd_adiab-coronal_temperature_k/unit_temperature)>1.d-10) &
        call mpistop('mhd_adiab differs from the requested isothermal T')
    case('adiabatic')
      if(.not.mhd_energy .or. .not.mhd_gravity) &
        call mpistop('adiabatic par must enable energy and gravity')
      if(mhd_thermal_conduction .or. mhd_radiative_cooling) &
        call mpistop('adiabatic par must disable conduction and radiative cooling')
    case('thermodynamic')
      if(.not.mhd_energy .or. .not.mhd_gravity .or. &
         .not.mhd_thermal_conduction .or. .not.mhd_radiative_cooling) &
        call mpistop('thermodynamic requires energy/gravity/conduction/cooling')
      is_thermodynamic=.true.
    case default
      call mpistop('unknown physics_model in usr_list')
    end select

    select case(trim(atmosphere_model))
    case('uniform','corona','chromosphere')
    case default
      call mpistop('atmosphere_model must be uniform, corona, or chromosphere')
    end select
    select case(trim(atmosphere_source))
    case('hydrostatic')
      if(trim(atmosphere_model)=='uniform' .and. trim(physics_model)/='zero_beta') &
        call mpistop('uniform atmosphere is supported only for zero_beta')
    case('relaxed_table')
      if(trim(physics_model)=='isothermal') &
        call mpistop('isothermal does not support relaxed_table')
      if(len_trim(relaxed_atmosphere_file)==0) &
        call mpistop('relaxed_table requires relaxed_atmosphere_file')
    case default
      call mpistop('atmosphere_source must be hydrostatic or relaxed_table')
    end select
  end subroutine validate_configuration

  subroutine init_b_data_constrained_boundary(boundaryname,qBunit)
    use mod_global_parameters
    use mod_comm_lib, only: mpistop
    use mod_data_driven_boundary, only: read_data_driven_boundary_frame

    character(len=*), intent(in) :: boundaryname
    double precision, intent(in) :: qBunit

    double precision :: snapshot_time,dx_km,dy_km,Bmax
    double precision, allocatable :: bframe(:,:,:)
    integer :: bnx,bny,amr_factor,nx_physical,ny_physical

    call read_data_driven_boundary_frame(boundaryname,snapshot_time, &
      bnx,bny,dx_km,dy_km,bframe)
    if(bnx<=0 .or. bny<=0) call mpistop('invalid data-driven boundary size')

    amr_factor = 2**(refine_max_level-1)
    nx_physical = domain_nx1*amr_factor
    ny_physical = domain_nx2*amr_factor
    if(bnx<nx_physical .or. mod(bnx-nx_physical,2)/=0) &
      call mpistop('boundary nx is incompatible with the finest data-constrained grid')
    if(bny<ny_physical .or. mod(bny-ny_physical,2)/=0) &
      call mpistop('boundary ny is incompatible with the finest data-constrained grid')
    pad_boundary1 = (bnx-nx_physical)/2
    pad_boundary2 = (bny-ny_physical)/2

    if(allocated(Bx0)) deallocate(Bx0,By0,Bz0)
    allocate(Bx0(bnx,bny),By0(bnx,bny),Bz0(bnx,bny))
    Bx0 = bframe(:,:,1)/qBunit
    By0 = bframe(:,:,2)/qBunit
    Bz0 = bframe(:,:,3)/qBunit
    deallocate(bframe)
    nx_boundary = bnx
    ny_boundary = bny
    Bmax = maxval(abs(Bz0))
    if(Bmax<=0.d0) call mpistop('zero Bz in data-driven boundary frame')

    if(mype==0) then
      write(*,*) 'data-constrained bottom boundary:',trim(boundaryname)
      write(*,*) 'snapshot_time [s]:',snapshot_time
      write(*,*) 'boundary frame:',nx_boundary,'by',ny_boundary,'pixels'
      write(*,*) 'horizontal padding:',pad_boundary1,pad_boundary2
      write(*,*) 'dx, dy [km]:',dx_km,dy_km
      write(*,*) 'max |Bz| in code units:',Bmax
    end if
  end subroutine init_b_data_constrained_boundary

  subroutine init_atmosphere()
    use mod_global_parameters
    use mod_solar_atmosphere, only: get_atm_para

    double precision, allocatable :: grav(:)
    double precision :: dz,zmin,zmax,zref,rhoref,temp0,integral
    integer :: j

    if(trim(atmosphere_source)=='relaxed_table') then
      call read_relaxed_atmosphere()
      return
    end if

    dz=dx(3,refine_max_level)
    zmin=xprobmin3-(dble(nghostcells)+0.5d0)*dz
    zmax=xprobmax3+(dble(nghostcells)+0.5d0)*dz
    zref=rho_reference_height_cm/unit_length
    zmin=min(zmin,zref-dz)
    zmax=max(zmax,zref+dz)
    n_atmosphere=max(2,ceiling((zmax-zmin)/dz)+1)
    dz=(zmax-zmin)/dble(n_atmosphere-1)

    if(allocated(atmosphere_z)) deallocate(atmosphere_z,atmosphere_rho)
    allocate(atmosphere_z(n_atmosphere),atmosphere_rho(n_atmosphere,3),grav(n_atmosphere))
    do j=1,n_atmosphere
      atmosphere_z(j)=zmin+dble(j-1)*dz
      grav(j)=gravity_at_height(atmosphere_z(j))
    end do

    rhoref=rho_reference_numberdensity_cm3/unit_numberdensity
    select case(trim(atmosphere_model))
    case('uniform')
      atmosphere_rho(:,1)=rhoref
      atmosphere_rho(:,2)=0.d0
      atmosphere_rho(:,3)=0.d0
    case('corona')
      temp0=coronal_temperature_k/unit_temperature
      do j=1,n_atmosphere
        integral=solar_gravity*solar_radius**2/temp0* &
          (1.d0/(solar_radius+zref)-1.d0/(solar_radius+atmosphere_z(j)))
        atmosphere_rho(j,1)=rhoref*exp(integral)
        atmosphere_rho(j,2)=temp0*atmosphere_rho(j,1)
        atmosphere_rho(j,3)=temp0
      end do
    case('chromosphere')
      call get_atm_para(atmosphere_z,atmosphere_rho(:,1),atmosphere_rho(:,2), &
        grav,n_atmosphere,trim(temperature_curve),zref,rhoref,atmosphere_rho(:,3), &
        clamp_low_T=.true.)
    end select
    deallocate(grav)
  end subroutine init_atmosphere

  subroutine read_relaxed_atmosphere()
    use mod_comm_lib, only: mpistop
    use mod_global_parameters

    character(len=max_line_len) :: line
    integer :: iu,ios,n,j
    double precision :: ztmp,rhotmp,ptmp

    n=0
    open(newunit=iu,file=trim(relaxed_atmosphere_file),status='old', &
      action='read',iostat=ios)
    if(ios/=0) call mpistop('cannot open relaxed_atmosphere_file')
    do
      read(iu,'(A)',iostat=ios) line
      if(ios/=0) exit
      line=adjustl(line)
      if(len_trim(line)==0 .or. line(1:1)=='#') cycle
      read(line,*,iostat=ios) ztmp,rhotmp,ptmp
      if(ios/=0) call mpistop('bad row in relaxed atmosphere table')
      n=n+1
    end do
    close(iu)
    if(n<2) call mpistop('relaxed atmosphere table needs at least two rows')

    if(allocated(atmosphere_z)) deallocate(atmosphere_z,atmosphere_rho)
    n_atmosphere=n
    allocate(atmosphere_z(n),atmosphere_rho(n,3))
    open(newunit=iu,file=trim(relaxed_atmosphere_file),status='old', &
      action='read',iostat=ios)
    if(ios/=0) call mpistop('cannot reopen relaxed_atmosphere_file')
    j=0
    do
      read(iu,'(A)',iostat=ios) line
      if(ios/=0) exit
      line=adjustl(line)
      if(len_trim(line)==0 .or. line(1:1)=='#') cycle
      read(line,*,iostat=ios) ztmp,rhotmp,ptmp
      if(ios/=0 .or. rhotmp<=0.d0 .or. ptmp<0.d0) &
        call mpistop('invalid rho or p in relaxed atmosphere table')
      if(trim(physics_model)/='zero_beta' .and. ptmp<=0.d0) &
        call mpistop('relaxed atmosphere pressure must be positive')
      j=j+1
      atmosphere_z(j)=ztmp
      atmosphere_rho(j,1)=rhotmp
      atmosphere_rho(j,2)=ptmp
      if(ptmp>0.d0) then
        atmosphere_rho(j,3)=ptmp/rhotmp
      else
        atmosphere_rho(j,3)=0.d0
      end if
    end do
    close(iu)

    do j=2,n
      if(atmosphere_z(j)<=atmosphere_z(j-1)) &
        call mpistop('relaxed atmosphere height must be strictly increasing')
    end do
    if(atmosphere_z(1)>xprobmin3 .or. atmosphere_z(n)<xprobmax3) &
      call mpistop('relaxed atmosphere table does not cover the physical domain')
  end subroutine read_relaxed_atmosphere

  subroutine init_bottom_thermodynamic_boundary()
    use mod_global_parameters

    double precision :: z,rho0,p0,temp0,dz
    integer :: i,j,k

    if(allocated(rho_bottom)) deallocate(rho_bottom,p_bottom)
    allocate(rho_bottom(nx_boundary,ny_boundary,nghostcells))
    allocate(p_bottom(nx_boundary,ny_boundary,nghostcells))
    dz=dx(3,refine_max_level)
    do k=1,nghostcells
      z=xprobmin3-(dble(k)-0.5d0)*dz
      call atmosphere_at_height(z,rho0,p0,temp0)
      do j=1,ny_boundary
        do i=1,nx_boundary
          rho_bottom(i,j,k)=rho0
          p_bottom(i,j,k)=p0
        end do
      end do
    end do
  end subroutine init_bottom_thermodynamic_boundary

  subroutine atmosphere_at_height(z,rho0,p0,temp0)
    double precision, intent(in) :: z
    double precision, intent(out) :: rho0,p0,temp0
    double precision :: f
    integer :: lo,hi,mid

    if(z<=atmosphere_z(1)) then
      rho0=atmosphere_rho(1,1)
      p0=atmosphere_rho(1,2)
      temp0=atmosphere_rho(1,3)
      return
    else if(z>=atmosphere_z(n_atmosphere)) then
      rho0=atmosphere_rho(n_atmosphere,1)
      p0=atmosphere_rho(n_atmosphere,2)
      temp0=atmosphere_rho(n_atmosphere,3)
      return
    end if

    lo=1
    hi=n_atmosphere
    do while(hi-lo>1)
      mid=(lo+hi)/2
      if(atmosphere_z(mid)<=z) then
        lo=mid
      else
        hi=mid
      end if
    end do
    f=(z-atmosphere_z(lo))/(atmosphere_z(hi)-atmosphere_z(lo))
    rho0=exp((1.d0-f)*log(atmosphere_rho(lo,1))+f*log(atmosphere_rho(hi,1)))
    if(atmosphere_rho(lo,2)>0.d0 .and. atmosphere_rho(hi,2)>0.d0) then
      p0=exp((1.d0-f)*log(atmosphere_rho(lo,2))+f*log(atmosphere_rho(hi,2)))
    else
      p0=(1.d0-f)*atmosphere_rho(lo,2)+f*atmosphere_rho(hi,2)
    end if
    temp0=(1.d0-f)*atmosphere_rho(lo,3)+f*atmosphere_rho(hi,3)
  end subroutine atmosphere_at_height

  subroutine transform_w_from_magnetic_snapshot(ixI^L,ixO^L,nw_in,w_in,x,w_out)
    use mod_comm_lib, only: mpistop
    use mod_global_parameters

    integer, intent(in) :: ixI^L,ixO^L,nw_in
    double precision, intent(in) :: w_in(ixI^S,1:nw_in),x(ixI^S,1:ndim)
    double precision, intent(out) :: w_out(ixI^S,1:nw)

    w_out=0.d0
    select case(nw_in)
    case(7)
      ! Energyless MHD/MF ordering: rho, v1:3, B1:3.
      w_out(ixO^S,mag(1:3))=w_in(ixO^S,5:7)
    case(8)
      ! Energy MHD ordering: rho, v1:3, e/p, B1:3.
      w_out(ixO^S,mag(1:3))=w_in(ixO^S,6:8)
    case default
      call mpistop('magnetic snapshot must contain 7 or 8 MHD variables')
    end select
  end subroutine transform_w_from_magnetic_snapshot

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    use mod_global_parameters

    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: Bsave(ixI^S,1:ndir)
    double precision :: rho0,p0,temp0
    integer :: ix1,ix2,ix3

    Bsave(ixO^S,1:ndir)=w(ixO^S,mag(1:ndir))
    do ix3=ixOmin3,ixOmax3
      do ix2=ixOmin2,ixOmax2
        do ix1=ixOmin1,ixOmax1
          call atmosphere_at_height(x(ix1,ix2,ix3,3),rho0,p0,temp0)
          w(ix1,ix2,ix3,rho_)=rho0
          w(ix1,ix2,ix3,mom(1:3))=0.d0
          w(ix1,ix2,ix3,mag(1:3))=Bsave(ix1,ix2,ix3,1:3)
          if(mhd_energy) w(ix1,ix2,ix3,p_)=p0
        end do
      end do
    end do
    call eos%to_conserved(ixI^L,ixO^L,w,x)
  end subroutine initonegrid_usr

  subroutine specialbound_usr(qdt,qt,ixI^L,ixO^L,iB,w,x)
    use mod_comm_lib, only: mpistop
    use mod_global_parameters

    integer, intent(in) :: ixI^L,ixO^L,iB
    double precision, intent(in) :: qdt,qt,x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: dxb1,dxb2,xlen1,xlen2
    integer :: ix1,ix2,ix3,ixbc1,ixbc2,k

    if(iB/=5) call mpistop('Only the lower x3 boundary is special in this demo')

    dxb1=dx(1,refine_max_level)
    dxb2=dx(2,refine_max_level)
    do ix3=ixOmin3,ixOmax3
      k=max(1,min(nghostcells,ixOmax3-ix3+1))
      do ix2=ixOmin2,ixOmax2
        do ix1=ixOmin1,ixOmax1
          xlen1=x(ix1,ix2,ix3,1)-xprobmin1+dble(pad_boundary1)*dxb1
          xlen2=x(ix1,ix2,ix3,2)-xprobmin2+dble(pad_boundary2)*dxb2
          ixbc1=max(1,min(nx_boundary,ceiling(xlen1/dxb1)))
          ixbc2=max(1,min(ny_boundary,ceiling(xlen2/dxb2)))
          w(ix1,ix2,ix3,rho_)=rho_bottom(ixbc1,ixbc2,k)
          w(ix1,ix2,ix3,mom(1:3))=0.d0
          w(ix1,ix2,ix3,mag(1))=Bx0(ixbc1,ixbc2)
          w(ix1,ix2,ix3,mag(2))=By0(ixbc1,ixbc2)
          w(ix1,ix2,ix3,mag(3))=Bz0(ixbc1,ixbc2)
          if(mhd_energy) w(ix1,ix2,ix3,p_)=p_bottom(ixbc1,ixbc2,k)
        end do
      end do
    end do

    ! Keep the observed vector field in the nearest ghost layer.  Extrapolate
    ! farther ghost layers from it and the interior, matching the established
    ! data-driven boundary treatment.
    do ix3=ixOmax3-1,ixOmin3,-1
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3,mag(1):mag(3)) = &
        (-3.d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3+4,mag(1):mag(3)) &
        +16.d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3+3,mag(1):mag(3)) &
        -36.d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3+2,mag(1):mag(3)) &
        +48.d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3+1,mag(1):mag(3)))/25.d0
    end do
    call eos%to_conserved(ixI^L,ixO^L,w,x)
  end subroutine specialbound_usr

  subroutine gravity_usr(ixI^L,ixO^L,wCT,x,gravity_field)
    use mod_global_parameters

    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim),wCT(ixI^S,1:nw)
    double precision, intent(out) :: gravity_field(ixI^S,ndim)

    gravity_field=0.d0
    gravity_field(ixO^S,3)=solar_gravity*(solar_radius/(solar_radius+x(ixO^S,3)))**2
  end subroutine gravity_usr

  double precision function gravity_at_height(z)
    double precision, intent(in) :: z
    gravity_at_height=solar_gravity*(solar_radius/(solar_radius+z))**2
  end function gravity_at_height

  subroutine heating_source_usr(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)
    use mod_global_parameters

    integer, intent(in) :: ixI^L,ixO^L,iw^LIM
    double precision, intent(in) :: qdt,qtC,qt,x(ixI^S,1:ndim),wCT(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    if(.not.is_thermodynamic) return
    w(ixO^S,e_)=w(ixO^S,e_)+qdt*heating_amplitude* &
      exp(-max(x(ixO^S,3),0.d0)/heating_scale_height)
  end subroutine heating_source_usr

  subroutine special_refine_grid(igrid,level,ixI^L,ixO^L,qt,w,x,refine,coarsen)
    use mod_global_parameters

    integer, intent(in) :: igrid,level,ixI^L,ixO^L
    double precision, intent(in) :: qt,w(ixI^S,1:nw),x(ixI^S,1:ndim)
    integer, intent(inout) :: refine,coarsen

    if(block%is_physical_boundary(5)) then
      refine=1
      coarsen=-1
    end if
  end subroutine special_refine_grid

  subroutine write_initial_atmosphere()
    use mod_global_parameters

    integer :: iu,ios,j
    double precision :: heat

    if(mype/=0) return
    open(newunit=iu,file='output/initial_atmosphere.dat',status='replace', &
      action='write',iostat=ios)
    if(ios/=0) then
      write(*,*) 'WARNING: cannot write output/initial_atmosphere.dat'
      return
    end if
    write(iu,'(A)') '# code units: z rho p T gravity heating'
    do j=1,n_atmosphere
      heat=0.d0
      if(is_thermodynamic) heat=heating_amplitude* &
        exp(-max(atmosphere_z(j),0.d0)/heating_scale_height)
      write(iu,'(6(ES24.16,1X))') atmosphere_z(j),atmosphere_rho(j,1), &
        atmosphere_rho(j,2),atmosphere_rho(j,3), &
        gravity_at_height(atmosphere_z(j)),heat
    end do
    close(iu)
  end subroutine write_initial_atmosphere

end module mod_usr
