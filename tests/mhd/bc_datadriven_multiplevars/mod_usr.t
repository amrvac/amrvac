module mod_usr
  use mod_mhd

  implicit none


contains

  !> Read this module's parameters from a file
  subroutine usr_params_read(files)
    character(len=*), intent(in) :: files(:)
    integer                      :: n

!    namelist /usr_list/ 
!
!    do n = 1, size(files)
!       open(unitpar, file=trim(files(n)), status="old")
!       read(unitpar, usr_list, end=111)
!111    close(unitpar)
!    end do

  end subroutine usr_params_read


  subroutine usr_init()
    call usr_params_read(par_files)

    unit_length=1d9 !Mm  
    unit_temperature=1.d6 ! in K
    unit_numberdensity=1.d9 ! in cm^-3

    usr_set_parameters=> initglobaldata_usr
    usr_init_one_grid => initonegrid_usr
    !usr_process_adv_grid => special_process

    call set_coordinate_system("Cartesian_3D")
    call usr_params_read(par_files)

    call mhd_activate()

  end subroutine usr_init


!  !!COPY VALUES FROM lower bc in dim 1 where it is set from vtk data 
!    subroutine special_process(igrid,level,ixI^L,ixO^L,qt,w,x)
!      use mod_global_parameters
!      integer, intent(in)             :: igrid,level,ixI^L,ixO^L
!      double precision, intent(in)    :: qt,x(ixI^S,1:ndim)
!      double precision, intent(inout) :: w(ixI^S,1:nw)
!      integer :: ix1
!
!      w(ixO^S,1:nw)=0.0d0
!      w(ixO^S,rho_)=1.0
!      w(ixO^S,e_)=1.0
!      do ix1=ixOmin1,ixOmax1
!        w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1):mag(3))=w(1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1):mag(3))
!      enddo
!    end subroutine special_process






  subroutine initglobaldata_usr

    if(mype==0) then
         print*, "UNIT TIME ", unit_time
         print*, "UNIT TEMP ", unit_temperature
         print*, "UNIT PRES ", unit_pressure
         print*, "UNIT DENS ", unit_density
         print*, "UNIT MAG  ", unit_magneticfield
         print*, "UNIT VELOCITY  ", unit_velocity
    endif


 end subroutine initglobaldata_usr


  subroutine initonegrid_usr(ixG^L,ix^L,w,x)
    use mod_global_parameters
    use mod_bc_data
    integer, intent(in) :: ixG^L,ix^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    double precision :: arr2d(ixmin2:ixmax2,ixmin3:ixmax3)
    logical, save :: first=.true.

    integer :: ix1, idir

    w(ix^S,1:nw)=0.0d0
    w(ix^S,rho_)=1.0
    w(ix^S,e_)=1.0

    !! set the three components of the mag field
    do idir=1,3
      arr2d(ixmin2:ixmax2, ixmin3:ixmax3) = bc_data_get_3d(bc_data_ix(mag(idir), 1), &
           x(ixmin1, ixmin2:ixmax2, ixmin3:ixmax3, 2), &
           x(ixmin1, ixmin2:ixmax2, ixmin3:ixmax3, 3), 0d0)
  
     do ix1=ixmin1,ixmax1
        w(ix1,ixmin2:ixmax2,ixmin3:ixmax3,mag(idir))=arr2d(ixmin2:ixmax2, ixmin3:ixmax3) 
     enddo  
   enddo 
  
  end subroutine initonegrid_usr




end module mod_usr
