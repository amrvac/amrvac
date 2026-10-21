!> B-only time-dependent magnetofriction from unified boundary frames.
!>
!> The initial magnetic field can come from either a PotentialField snapshot
!> or a legacy-MHD MagnetofrictionalRelaxation snapshot.  Only B is imported
!> into the standalone mf physics; its three ``mom`` variables are the
!> artificial magnetofrictional velocity, not MHD momentum or plasma velocity.
module mod_usr
  use mod_mf
  use mod_data_driven_boundary, only: data_driven_boundary_series, &
    read_data_driven_boundary_series, interpolate_data_driven_boundary_scaled
  implicit none

  type(data_driven_boundary_series), save :: boundary_series
  character(len=256) :: boundary_series_dir
  double precision :: driving_time_scale
  integer :: boundary_frame_count
  integer, save :: nx_boundary=0,ny_boundary=0
  integer, save :: pad_boundary1=0,pad_boundary2=0

contains

  subroutine usr_init()
    use mod_usr_methods

    ! Keep the restart and the boundary magnetogram on the same standard
    ! coronal normalization as the PotentialField case. Standalone mf derives all
    ! remaining units, including the magnetic-field unit, from these values.
    unit_length        = 1.d9 ! cm
    unit_temperature   = 1.d6 ! K
    unit_numberdensity = 1.d9 ! cm^-3

    call set_coordinate_system('Cartesian_3D')

    usr_set_parameters => initglobaldata_usr
    usr_init_one_grid  => initonegrid_usr
    usr_transform_w    => transform_w_from_magnetic_snapshot
    usr_special_bc     => specialbound_usr
    usr_refine_grid    => special_refine_grid

    call mf_activate()
  end subroutine usr_init

  subroutine usr_params_read(files)
    use mod_global_parameters

    character(len=*), intent(in) :: files(:)
    integer :: n

    namelist /usr_list/ boundary_series_dir,boundary_frame_count,driving_time_scale

    boundary_series_dir = '../MagneticBoundary/Evolution/Sequence'
    boundary_frame_count = 2
    driving_time_scale = 12.d0

    do n=1,size(files)
      open(unitpar,file=trim(files(n)),status='old')
      read(unitpar,usr_list,end=111)
111   close(unitpar)
    end do
  end subroutine usr_params_read

  subroutine initglobaldata_usr()
    use mod_comm_lib, only: mpistop
    use mod_global_parameters

    call usr_params_read(par_files)
    if(restart_from_file==undefined) &
      call mpistop('TimeDependentMagnetofriction must restart from a magnetic snapshot')

    call init_b_mfr_data_driven_series()
  end subroutine initglobaldata_usr

  subroutine init_b_mfr_data_driven_series()
    use mod_global_parameters
    use mod_comm_lib, only: mpistop

    integer :: amr_factor,nx_physical,ny_physical

    if(driving_time_scale<=0.d0) call mpistop('driving_time_scale must be positive')
    call read_data_driven_boundary_series(trim(boundary_series_dir),boundary_series, &
      expected_nframe=boundary_frame_count)

    amr_factor = 2**(refine_max_level-1)
    nx_physical = domain_nx1*amr_factor
    ny_physical = domain_nx2*amr_factor
    if(boundary_series%nx<nx_physical .or. mod(boundary_series%nx-nx_physical,2)/=0) &
       call mpistop('boundary nx is incompatible with the finest MF grid')
    if(boundary_series%ny<ny_physical .or. mod(boundary_series%ny-ny_physical,2)/=0) &
       call mpistop('boundary ny is incompatible with the finest MF grid')
    pad_boundary1 = (boundary_series%nx-nx_physical)/2
    pad_boundary2 = (boundary_series%ny-ny_physical)/2
    nx_boundary = boundary_series%nx
    ny_boundary = boundary_series%ny

    if(mype==0) then
      write(*,*) 'time-dependent MF boundary series:',trim(boundary_series_dir)
      write(*,*) 'frames:',boundary_series%nframe
      write(*,*) 'time range [s]:',boundary_series%times(1),boundary_series%times(boundary_series%nframe)
      write(*,*) 'boundary frame:',nx_boundary,'by',ny_boundary,'pixels'
      write(*,*) 'horizontal padding:',pad_boundary1,pad_boundary2
      write(*,*) 'dx, dy [km]:',boundary_series%dx_km,boundary_series%dy_km
      write(*,*) 'driving time scale:',driving_time_scale
    end if
  end subroutine init_b_mfr_data_driven_series

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    use mod_global_parameters

    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    ! In standalone mf these are artificial frictional velocities.  A
    ! potential/MFR restart contributes only its magnetic field.
    w(ixO^S,mom(1:3))=0.d0
    if(mf_glm) w(ixO^S,psi_)=0.d0
  end subroutine initonegrid_usr

  subroutine transform_w_from_magnetic_snapshot(ixI^L,ixO^L,nw_in,w_in,x,w_out)
    use mod_comm_lib, only: mpistop
    use mod_global_parameters

    integer, intent(in) :: ixI^L,ixO^L,nw_in
    double precision, intent(in) :: w_in(ixI^S,1:nw_in),x(ixI^S,1:ndim)
    double precision, intent(out) :: w_out(ixI^S,1:nw)

    w_out=0.d0
    select case(nw_in)
    case(7)
      ! Energyless MHD ordering: rho, v1:3, B1:3.
      w_out(ixO^S,mag(1:3))=w_in(ixO^S,5:7)
    case(8)
      ! Energy MHD ordering: rho, v1:3, e/p, B1:3.
      w_out(ixO^S,mag(1:3))=w_in(ixO^S,6:8)
    case default
      call mpistop('TMF input snapshot must contain 7 or 8 MHD variables')
    end select
    w_out(ixO^S,mom(1:3))=0.d0
    if(mf_glm) w_out(ixO^S,psi_)=0.d0
  end subroutine transform_w_from_magnetic_snapshot

  subroutine specialbound_usr(qdt,qt,ixI^L,ixO^L,iB,w,x)
    use mod_global_parameters

    integer, intent(in) :: ixI^L,ixO^L,iB
    double precision, intent(in) :: qdt,qt,x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: dxb1,dxb2,xlen1,xlen2
    double precision, allocatable :: bframe(:,:,:)
    integer :: ix1,ix2,ix3,ixbc1,ixbc2

    select case(iB)
    case(5)
      if(mf_glm) w(ixO^S,psi_)=0.d0
      allocate(bframe(nx_boundary,ny_boundary,3))
      call interpolate_data_driven_boundary_scaled(boundary_series,qt,unit_time, &
        driving_time_scale,bframe)
      bframe=bframe/unit_magneticfield
      ! No physical velocity series is supplied.  Continue the standalone-mf
      ! artificial velocity smoothly through the moving lower boundary rather
      ! than imposing the antisymmetric, zero-face condition used by static
      ! magnetofrictional relaxation.
      do ix3=ixOmax3,ixOmin3,-1
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3,mom(1):mom(3)) = &
           0.12d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3+5,mom(1):mom(3)) &
          -0.76d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3+4,mom(1):mom(3)) &
          +2.08d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3+3,mom(1):mom(3)) &
          -3.36d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3+2,mom(1):mom(3)) &
          +2.92d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3+1,mom(1):mom(3))
      end do

      ! Boundary frames include horizontal ghost padding. Map AMRVAC cell
      ! centers using the finest-level boundary coordinate system.
      dxb1 = dx(1,refine_max_level)
      dxb2 = dx(2,refine_max_level)
      ! The observed B is imposed only in the ghost layer adjacent to the
      ! physical domain.  More distant ghost layers are extrapolated below.
      do ix3=ixOmax3,ixOmax3
        do ix2=ixOmin2,ixOmax2
          do ix1=ixOmin1,ixOmax1
            xlen1 = x(ix1,ix2,ix3,1)-xprobmin1+dble(pad_boundary1)*dxb1
            xlen2 = x(ix1,ix2,ix3,2)-xprobmin2+dble(pad_boundary2)*dxb2
            ixbc1 = ceiling(xlen1/dxb1)
            ixbc2 = ceiling(xlen2/dxb2)
            ixbc1 = max(1,min(nx_boundary,ixbc1))
            ixbc2 = max(1,min(ny_boundary,ixbc2))
            w(ix1,ix2,ix3,mag(1)) = bframe(ixbc1,ixbc2,1)
            w(ix1,ix2,ix3,mag(2)) = bframe(ixbc1,ixbc2,2)
            w(ix1,ix2,ix3,mag(3)) = bframe(ixbc1,ixbc2,3)
          end do
        end do
      end do
      deallocate(bframe)

      do ix3=ixOmax3-1,ixOmin3,-1
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3,mag(1):mag(3)) = &
           0.12d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3+5,mag(1):mag(3)) &
          -0.76d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3+4,mag(1):mag(3)) &
          +2.08d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3+3,mag(1):mag(3)) &
          -3.36d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3+2,mag(1):mag(3)) &
          +2.92d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3+1,mag(1):mag(3))
      end do
    case default
      call mpistop('Only the lower x3 boundary is special in this demo')
    end select
  end subroutine specialbound_usr

  subroutine special_refine_grid(igrid,level,ixI^L,ixO^L,qt,w,x,refine,coarsen)
    use mod_global_parameters

    integer, intent(in) :: igrid,level,ixI^L,ixO^L
    double precision, intent(in) :: qt,w(ixI^S,1:nw),x(ixI^S,1:ndim)
    integer, intent(inout) :: refine,coarsen

    ! Every block touching the driven lower boundary must be refined to the
    ! maximum level so its cell centers map one-to-one onto the magnetogram.
    if(minval(x(ixO^S,3))<xprobmin3+dx(3,level)) then
      refine = 1
      coarsen = -1
    end if
  end subroutine special_refine_grid

end module mod_usr
