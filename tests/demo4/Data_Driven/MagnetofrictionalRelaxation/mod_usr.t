!> Magnetofrictional relaxation from a Python V1 data-driven boundary frame.
!>
!> This case is intended to restart from the PotentialField demo output.  The
!> bottom magnetic field is fixed from the same boundary frame written by
!> tools/python/notebooks/data_driven_pipeline.ipynb:
!> snapshot_time, nx, ny, dx, dy, Bx, By, Bz.
module mod_usr
  use mod_mhd
  implicit none

  double precision, allocatable, save :: Bx0(:,:),By0(:,:),Bz0(:,:)
  character(len=256) :: boundary_filename
  integer, save :: nx_boundary=0,ny_boundary=0
  integer, save :: pad_boundary1=0,pad_boundary2=0

contains

  subroutine usr_init()
    use mod_usr_methods

    ! Keep the restart and the boundary magnetogram on the same standard
    ! coronal normalization as the PotentialField case. MHD derives all
    ! remaining units, including the magnetic-field unit, from these values.
    unit_length        = 1.d9 ! cm
    unit_temperature   = 1.d6 ! K
    unit_numberdensity = 1.d9 ! cm^-3

    call set_coordinate_system('Cartesian_3D')

    usr_set_parameters => initglobaldata_usr
    usr_init_one_grid  => initonegrid_usr
    usr_special_bc     => specialbound_usr
    usr_refine_grid    => special_refine_grid

    call mhd_activate()
  end subroutine usr_init

  subroutine usr_params_read(files)
    use mod_global_parameters

    character(len=*), intent(in) :: files(:)
    integer :: n

    namelist /usr_list/ boundary_filename

    boundary_filename = 'boundary_single/B_0001.dat'

    do n=1,size(files)
      open(unitpar,file=trim(files(n)),status='old')
      read(unitpar,usr_list,end=111)
111   close(unitpar)
    end do
  end subroutine usr_params_read

  subroutine initglobaldata_usr()
    use mod_global_parameters

    call usr_params_read(par_files)

    call init_b_mfr_data_driven_boundary(trim(boundary_filename),unit_magneticfield)
  end subroutine initglobaldata_usr

  subroutine init_b_mfr_data_driven_boundary(boundaryname,qBunit)
    use mod_global_parameters
    use mod_comm_lib, only: mpistop
    use mod_data_driven_boundary, only: read_data_driven_boundary_frame

    character(len=*), intent(in) :: boundaryname
    double precision, intent(in) :: qBunit

    double precision :: snapshot_time,dx_km,dy_km,Bmax
    double precision, allocatable :: bframe(:,:,:)
    integer :: bnx,bny,amr_factor,nx_physical,ny_physical

    call read_data_driven_boundary_frame(boundaryname,snapshot_time,bnx,bny,dx_km,dy_km,bframe)
    if(bnx<=0 .or. bny<=0) call mpistop('invalid data-driven boundary size')

    amr_factor = 2**(refine_max_level-1)
    nx_physical = domain_nx1*amr_factor
    ny_physical = domain_nx2*amr_factor
    if(bnx<nx_physical .or. mod(bnx-nx_physical,2)/=0) &
       call mpistop('boundary nx is incompatible with the finest MF grid')
    if(bny<ny_physical .or. mod(bny-ny_physical,2)/=0) &
       call mpistop('boundary ny is incompatible with the finest MF grid')
    pad_boundary1 = (bnx-nx_physical)/2
    pad_boundary2 = (bny-ny_physical)/2

    if(allocated(Bx0)) deallocate(Bx0)
    if(allocated(By0)) deallocate(By0)
    if(allocated(Bz0)) deallocate(Bz0)
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
      write(*,*) 'magnetofrictional bottom boundary:',trim(boundaryname)
      write(*,*) 'snapshot_time [s]:',snapshot_time
      write(*,*) 'boundary frame:',nx_boundary,'by',ny_boundary,'pixels'
      write(*,*) 'horizontal padding:',pad_boundary1,pad_boundary2
      write(*,*) 'dx, dy [km]:',dx_km,dy_km
      write(*,*) 'max |Bz| in code units:',Bmax
    end if
  end subroutine init_b_mfr_data_driven_boundary

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    use mod_comm_lib, only: mpistop
    use mod_global_parameters

    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    call mpistop('MagnetofrictionalRelaxation must restart from the PotentialField snapshot')
  end subroutine initonegrid_usr

  subroutine specialbound_usr(qdt,qt,ixI^L,ixO^L,iB,w,x)
    use mod_global_parameters

    integer, intent(in) :: ixI^L,ixO^L,iB
    double precision, intent(in) :: qdt,qt,x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: dxb1,dxb2,xlen1,xlen2
    integer :: ix1,ix2,ix3,ixbc1,ixbc2,af

    select case(iB)
    case(5)
      ! Density is a fixed positive normalization used only by the legacy
      ! magnetofrictional CFL estimate. The momentum reflection pins the
      ! frictional velocity to zero at the lower physical face.
      w(ixO^S,rho_) = one
      w(ixO^S,mom(1)) = -w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmax3+nghostcells:ixOmax3+1:-1,mom(1))
      w(ixO^S,mom(2)) = -w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmax3+nghostcells:ixOmax3+1:-1,mom(2))
      w(ixO^S,mom(3)) = -w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmax3+nghostcells:ixOmax3+1:-1,mom(3))

      ! Boundary frames include horizontal ghost padding. Map AMRVAC cell
      ! centers using the finest-level boundary coordinate system.
      dxb1 = dx(1,refine_max_level)
      dxb2 = dx(2,refine_max_level)
      do ix3=ixOmin3,ixOmax3
        do ix2=ixOmin2,ixOmax2
          do ix1=ixOmin1,ixOmax1
            xlen1 = x(ix1,ix2,ix3,1)-xprobmin1+dble(pad_boundary1)*dxb1
            xlen2 = x(ix1,ix2,ix3,2)-xprobmin2+dble(pad_boundary2)*dxb2
            ixbc1 = ceiling(xlen1/dxb1)
            ixbc2 = ceiling(xlen2/dxb2)
            ixbc1 = max(1,min(nx_boundary,ixbc1))
            ixbc2 = max(1,min(ny_boundary,ixbc2))
            w(ix1,ix2,ix3,mag(1)) = Bx0(ixbc1,ixbc2)
            w(ix1,ix2,ix3,mag(2)) = By0(ixbc1,ixbc2)
            w(ix1,ix2,ix3,mag(3)) = Bz0(ixbc1,ixbc2)
          end do
        end do
      end do

      af = 1
      do ix3=ixOmax3-af,ixOmin3,-1
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3,mag(1):mag(3)) = &
          ( -3.d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3+4,mag(1):mag(3)) &
           +16.d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3+3,mag(1):mag(3)) &
           -36.d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3+2,mag(1):mag(3)) &
           +48.d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ix3+1,mag(1):mag(3)))/25.d0
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
