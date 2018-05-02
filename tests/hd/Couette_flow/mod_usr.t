! IC = 1 (and 2) : planar Couette flow (driving upper edge while lower edge is set fixed)
! IC = 3 : polar Couette flow (driving outer cylinder while inner cylinder is set fixed)

module mod_usr
  use mod_hd
  use mod_viscosity

  implicit none
  double precision :: Reynolds
  integer :: IC
  double precision :: v0, rho0

contains

  subroutine usr_init()
    character(len=30) :: geometry

    usr_set_parameters => initglobaldata_usr
    usr_init_one_grid  => initonegrid_usr
    usr_special_bc     => my_bounds
    usr_internal_bc    => no_vr
    usr_aux_output     => specialvar_output
    usr_add_aux_names  => specialvarnames_output

    call hd_activate()
    call params_read(par_files)

    if (IC==1 .or. IC==2) geometry='Cartesian_2D'
    if (IC==3) geometry='polar_2D'
    call set_coordinate_system(geometry)

  end subroutine usr_init

  !> Read parameters from a file
  subroutine params_read(files)
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /my_list/ Reynolds, IC

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, my_list, end=111)
111    close(unitpar)
    end do

  end subroutine params_read


  subroutine initglobaldata_usr
    double precision :: vc_tau, Mach=1.d0 ! fiducial value of the Mach # (see below)
    v0     = one ! default value
    rho0   = one ! default value
    vc_mu  = 1.d0 / Reynolds ! ( rho0 * v0 * (xprobmax2-xprobmin2) ) / Reynolds
    if ( (IC==1 .or. IC==2) .and. (xprobmax2-xprobmin2/=one) ) call mpistop("Beware normalization")
    vc_tau = one / vc_mu ! (xprobmax2-xprobmin2)**two/(vc_mu/rho0) ! characteristic viscous time scale
    time_max = 0.5d0 * vc_tau ! relaxation time
    print*, '- - -'
    print*, 'Reynolds # ', Reynolds
    print*, 'The simulation duration has been set to ', time_max
    print*, '- - -'
    tsave(1,2)=time_max/1000.d0
    tsave(2,2)=time_max/100.d0
    tsave(3,2)=time_max/10.d0

    ! This does not matter as long as there is no gradient of pressure,
    ! which should always be the case since the flow is uniform and we work
    ! with a polytropic assumption
    hd_adiab = one / Mach**two

    if (hd_energy) call mpistop("This is not the classic Couette flow")

  end subroutine initglobaldata_usr

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    if (IC==1) then
      w(ixO^S,rho_)   = rho0
      w(ixO^S,mom(1)) = zero
      w(ixO^S,mom(2)) = zero
    elseif (IC==2) then
      w(ixO^S,rho_)   = rho0
      where (x(ixO^S,2)>0.5d0)
        w(ixO^S,mom(1)) = rho0 * v0
      elsewhere
        w(ixO^S,mom(1)) = zero
      end where
      w(ixO^S,mom(2)) = zero
    elseif (IC==3) then
      w(ixO^S,rho_)   = rho0
      ! where (x(ixO^S,1)>1.5d0)
      !   w(ixO^S,mom(2)) = rho0 * v0
      ! elsewhere
        w(ixO^S,mom(2)) = zero
      ! end where
      w(ixO^S,mom(1)) = zero
    endif

  end subroutine initonegrid_usr

  ! Top boundary sheared, bottom fixed
  subroutine my_bounds(qt,ixG^L,ixB^L,iB,w,x)

    integer, intent(in) :: ixG^L, ixB^L, iB
    double precision, intent(in) :: qt, x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)
    integer :: i

    if (IC==1 .or. IC==2) then
      select case (iB)
      case(1) ! left continuous
        w(ixBmax1,:,rho_)   = w(ixBmax1+1,:,rho_)
        w(ixBmax1,:,mom(1)) = w(ixBmax1+1,:,mom(1))
        w(ixBmax1,:,mom(2)) = zero !w(ixBmax1+1,:,mom(2))
        do i=ixBmin1,ixBmax1-1
          w(i,:,:)      = w(ixBmax1,:,:)
        enddo
      case(2) ! right no inflow
        w(ixBmin1,:,rho_)   = w(ixBmin1-1,:,rho_)
        where (w(ixBmin1-1,:,mom(1))>zero)
          w(ixBmin1,:,mom(1)) = w(ixBmin1-1,:,mom(1))
        elsewhere
          w(ixBmin1,:,mom(1)) = zero
        end where
        w(ixBmin1,:,mom(2)) = zero !w(ixBmin1-1,:,mom(2))
        do i=ixBmin1+1,ixBmax1
          w(i,:,:)      = w(ixBmin1,:,:)
        enddo
      case(3) ! bottom
        w(:,ixBmax2,rho_)   = w(:,ixBmax2+1,rho_) ! continuous
        w(:,ixBmax2,mom(1)) = zero
        w(:,ixBmax2,mom(2)) = zero
        do i=ixBmin2,ixBmax2-1
          w(:,i,:)      = w(:,ixBmax2,:)
        enddo
      case(4) ! top
        w(:,ixBmin2,rho_)   = w(:,ixBmin2-1,rho_) ! continuous
        w(:,ixBmin2,mom(1)) = w(:,ixBmin2-1,rho_) * v0
        w(:,ixBmin2,mom(2)) = zero
        do i=ixBmin2+1,ixBmax2
          w(:,i,:)      = w(:,ixBmin2,:)
        enddo
      case default
        call mpistop('BC not implemented')
      end select
    elseif (IC==3) then
      select case (iB)
      case(1) ! inner fixed
        w(ixBmax1,:,rho_)   =  w(ixBmax1+1,:,rho_)
        w(ixBmax1,:,mom(1)) = -w(ixBmax1+1,:,mom(1))
        w(ixBmax1,:,mom(2)) =  zero
        do i=ixBmin1,ixBmax1-1
          w(i,:,:)      =  w(ixBmax1,:,:)
        enddo
      case(2) ! outer rotating
        w(ixBmin1,:,rho_)   =  w(ixBmin1-1,:,rho_)
        w(ixBmin1,:,mom(1)) = -w(ixBmin1-1,:,mom(1))
        w(ixBmin1,:,mom(2)) =  w(ixBmin1-1,:,rho_) * v0
        do i=ixBmin1+1,ixBmax1
          w(i,:,:)      =  w(ixBmin1,:,:)
        enddo
      case default
        call mpistop('BC not implemented')
      end select
    endif

  end subroutine my_bounds

  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
  use mod_physics
  integer, intent(in)                :: ixI^L,ixO^L
  double precision, intent(in)       :: x(ixI^S,1:ndim)
  double precision                   :: w(ixI^S,nw+nwauxio)
  double precision                   :: normconv(0:nw+nwauxio)
  double precision :: v(ixI^S,ndir), divV(ixI^S)
  integer :: i

  do i=1,ndir
    v(ixI^S,i)=w(ixI^S,mom(i))
  enddo
  call divvector(v,ixI^L,ixO^L,divV)
  w(ixO^S,nw+1) = divV(ixO^S)

  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
    character(len=*) :: varnames
    varnames='divV'

  end subroutine specialvarnames_output

  subroutine no_vr(level,qt,ixI^L,ixO^L,w,x)

    integer, intent(in) :: ixI^L,ixO^L,level
    double precision, intent(in) :: qt
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in) :: x(ixI^S,1:ndim)

    if (IC/=3) call mpistop("Makes sense only for polar case where radial speed is neglected")
    w(ixI^S,mom(1))=zero

  end subroutine no_vr

end module mod_usr
