!> Program to extrapolate linear force-free fields in 3D Cartesian coordinates,
!> based on exact Green function method (Chiu & Hilton 1977 ApJ 212,873).
!>
!> Usage:
!> 1 In the subroutine usr_set_parameters of mod_usr.t:
!>  To extrapolate a linear force free field from a observed magnetogram 
!>  prepared in a data file, e.g., 'hmiM720sxxxx.dat' replace 
!>  call init_bc_fff_data('hmiM720sxxxx.dat',unit_length,unit_magneticfield)
!>  'hmiM720sxxxx.dat' must be a binary file containing nx1,nx2,xc1,xc2,dxm1,
!>  dxm2, Bz0(nx1,nx2). Integers nx1 and nx2 give the resolution of the 
!>  uniform-grid magentogram. Others are double-precision floats. xc1 and xc2
!>  are coordinates of the central point of the magnetogram. dxm1 and dxm2 
!>  are the cell sizes for each direction, Bz0 is the vertical conponent 
!>  of magetic field on the solar surface from observations.
!>2 In the subroutine usr_init_one_grid of mod_usr.t,
!>  add lines like:
!>
!>  double precision :: Bf(ixG^S,1:ndir), alpha, zshift
!>
!>  alpha=0.d0     ! potential field
!>  !alpha=0.08d0  ! non-potential linear force-free field
!>  zshift=0.05d0  ! lift your box zshift heigher to the bottom magnetogram
!>  call calc_lin_fff(ixG^L,ix^L,Bf,x,alpha,zshift) 
!>
!>3 Notice that the resolution of input magnetogram must be better than the best
!>  resolution of your AMR grid to have a good behavior close to the bottom layer
module mod_lfff
  implicit none
  
  integer, save :: nx1,nx2
  double precision, save :: Bzmax,darea
  double precision, allocatable, save :: Bz0(:,:)
  double precision, allocatable, save :: xa1(:),xa2(:)
  
contains

  subroutine init_b_fff_data(magnetogramname,qLunit,qBunit)
    use mod_global_parameters

    double precision, intent(in) :: qLunit,qBunit
    double precision :: xc1,xc2,dxm1,dxm2
    integer, dimension(MPI_STATUS_SIZE) :: statuss
    integer :: file_handle,i
    character(len=*), intent(in) :: magnetogramname
    logical :: aexist
    ! nx1,nx2 are numbers of cells for each direction
    ! xc1,xc2 are coordinates of the central point of the magnetogram
    ! dxm1,dxm2 are cell sizes for each direction
    ! Bz0 is the 2D Bz magnetogram
    inquire(file=magnetogramname,exist=aexist)
    if(.not. aexist) then
      if(mype==0) write(*,'(2a)') "can not find file:",magnetogramname
      call mpistop("no input magnetogram----init_b_fff_data")
    end if
    call MPI_FILE_OPEN(icomm,magnetogramname,MPI_MODE_RDONLY,MPI_INFO_NULL,&
                       file_handle,ierrmpi)
    call MPI_FILE_READ_ALL(file_handle,nx1,1,MPI_INTEGER,statuss,ierrmpi)
    call MPI_FILE_READ_ALL(file_handle,nx2,1,MPI_INTEGER,statuss,ierrmpi)
    allocate(Bz0(nx1,nx2))
    call MPI_FILE_READ_ALL(file_handle,xc1,1,MPI_DOUBLE_PRECISION,statuss,ierrmpi)
    call MPI_FILE_READ_ALL(file_handle,xc2,1,MPI_DOUBLE_PRECISION,statuss,ierrmpi)
    call MPI_FILE_READ_ALL(file_handle,dxm1,1,MPI_DOUBLE_PRECISION,statuss,ierrmpi)
    call MPI_FILE_READ_ALL(file_handle,dxm2,1,MPI_DOUBLE_PRECISION,statuss,ierrmpi)
    call MPI_FILE_READ_ALL(file_handle,Bz0,nx1*nx2,MPI_DOUBLE_PRECISION,&
                           statuss,ierrmpi)
    call MPI_FILE_CLOSE(file_handle,ierrmpi)
    allocate(xa1(nx1))
    allocate(xa2(nx2))
    xa1(nx1/2)=xc1
    xa2(nx2/2)=xc2
    do i=nx1/2+1,nx1
      xa1(i)=xa1(nx1/2)+dble(i-nx1/2)*dxm1
    enddo
    do i=nx1/2-1,1,-1
      xa1(i)=xa1(nx1/2)+dble(i-nx1/2)*dxm1
    enddo
    do i=nx2/2+1,nx2
      xa2(i)=xa2(nx2/2)+dble(i-nx2/2)*dxm2
    enddo
    do i=nx2/2-1,1,-1
      xa2(i)=xa2(nx2/2)+dble(i-nx2/2)*dxm2
    enddo
    ! declare and define global variables Lunit and Bunit to be your length unit in
    ! cm and magnetic strength unit in Gauss first
    dxm1=dxm1/qLunit
    dxm2=dxm2/qLunit
    xa1=xa1/qLunit
    xa2=xa2/qLunit
    darea=dxm1*dxm2
    Bz0=Bz0/qBunit
    Bzmax=maxval(dabs(Bz0(:,:)))
    
    ! normalize b
    Bz0=Bz0/Bzmax
    if(mype==0) then
      print*,'magnetogram xrange:',minval(xa1),maxval(xa1)
      print*,'magnetogram yrange:',minval(xa2),maxval(xa2)
    end if
    
    if(mype==0) then
      print*,'extrapolating 3D force-free field from an observed Bz '
      print*,'magnetogram of',nx1,'by',nx2,'pixels. Bzmax=',Bzmax
    endif
  
  end subroutine init_b_fff_data

{^IFTHREED
  subroutine calc_lin_fff(ixI^L,ixO^L,Bf,x,alpha,zshift,idir)
  ! PURPOSE: 
  ! Calculation to determine linear FFF from the field on 
  ! the lower boundary (Chiu and Hilton 1977 ApJ 212,873). 
  ! NOTE: Only works for Cartesian coordinates 
  ! INPUT: Bf,x
  ! OUTPUT: updated b in w 
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    integer, optional, intent(in) :: idir
    double precision, intent(in) :: x(ixI^S,1:ndim),alpha,zshift
    double precision, intent(inout) :: Bf(ixI^S,1:ndir)

    double precision, dimension(ixO^S) :: cos_az,sin_az,zk,bigr,r,r2,r3,cos_ar,sin_ar,g,dgdz
    double precision :: gx(ixO^S,1:ndim),twopiinv
    integer :: idim,ixp1,ixp2

    Bf=0.d0
    twopiinv = 0.5d0/dpi*Bzmax*darea
    ! get cos and sin arrays
    zk(ixO^S)=x(ixO^S,3)-xprobmin3+zshift
    cos_az(ixO^S)=dcos(alpha*zk(ixO^S))
    sin_az(ixO^S)=dsin(alpha*zk(ixO^S))
    ! looping Bz0 pixels
    do ixp2=1,nx2
      do ixp1=1,nx1
        bigr(ixO^S)=dsqrt((x(ixO^S,1)-xa1(ixp1))**2+&
                          (x(ixO^S,2)-xa2(ixp2))**2)
        r2=bigr**2+zk**2
        r=dsqrt(r2)
        r3=r**3
        cos_ar=dcos(alpha*r)
        sin_ar=dsin(alpha*r)
        where(bigr/=0.d0)
          bigr=1.d0/bigr
        end where
        where(r/=0.d0)
          r=1.d0/r
        end where
        where(r2/=0.d0)
          r2=1.d0/r2
        end where
        where(r3/=0.d0)
          r3=1.d0/r3
        end where
        g=(zk*cos_ar*r-cos_az)*bigr
        dgdz=(cos_ar*(r-zk**2*r3)-alpha*zk**2*sin_ar*r2+alpha*sin_az)*bigr
        do idim=1,ndim
          if(present(idir).and.idim/=idir) cycle
          select case(idim)
          case(1)
            Bf(ixO^S,1)=Bf(ixO^S,1)+Bz0(ixp1,ixp2)*((x(ixO^S,1)-xa1(ixp1))*dgdz(ixO^S)&
                     +alpha*g(ixO^S)*(x(ixO^S,2)-xa2(ixp2)))*bigr(ixO^S)
          case(2)
            Bf(ixO^S,2)=Bf(ixO^S,2)+Bz0(ixp1,ixp2)*((x(ixO^S,2)-xa2(ixp2))*dgdz(ixO^S)&
                     -alpha*g(ixO^S)*(x(ixO^S,1)-xa1(ixp1)))*bigr(ixO^S)
          case(3)
            Bf(ixO^S,3)=Bf(ixO^S,3)+Bz0(ixp1,ixp2)*(zk(ixO^S)*cos_ar(ixO^S)*r3(ixO^S)+alpha*&
                                        zk(ixO^S)*sin_ar(ixO^S)*r2(ixO^S))
          end select
        end do
      end do
    end do
    Bf(ixO^S,:)=Bf(ixO^S,:)*twopiinv

  end subroutine calc_lin_fff

  subroutine get_potential_field_potential(ixI^L,ixO^L,potential,x,zshift)
  ! PURPOSE: 
  ! Calculation to determine linear FFF from the field on 
  ! the lower boundary (Chiu and Hilton 1977 ApJ 212,873). 
  ! NOTE: Only works for Cartesian coordinates 
  ! INPUT: Bf,x
  ! OUTPUT: updated b in w 
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim),zshift
    double precision, intent(inout) :: potential(ixI^S)

    double precision, dimension(ixO^S) :: zk,bigr
    integer :: ixp1,ixp2

    zk(ixO^S)=x(ixO^S,3)-xprobmin3+zshift
    potential=0.d0
    ! looping Bz0 pixels
    do ixp2=1,nx2
      do ixp1=1,nx1
        bigr(ixO^S)=dsqrt((x(ixO^S,1)-xa1(ixp1))**2+&
                          (x(ixO^S,2)-xa2(ixp2))**2+&
                          zk(ixO^S)**2)
        potential(ixO^S)=potential(ixO^S)+0.5d0*Bz0(ixp1,ixp2)/bigr*darea/dpi
      end do
    end do
  end subroutine get_potential_field_potential
}
end module mod_lfff
