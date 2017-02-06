module mod_gravity
  implicit none

  !> source split or not
  logical :: grav_split= .false.

  !> Index of the density (in the w array)
  integer, private, parameter              :: rho_ = 1

  !> Indices of the momentum density
  integer, allocatable, private, protected :: mom(:)

  !> Index of the energy density (-1 if not present)
  integer, private, protected              :: e_

contains
  !> Read this module's parameters from a file
  subroutine grav_params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /grav_list/ grav_split

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, grav_list, end=111)
111    close(unitpar)
    end do

  end subroutine grav_params_read

  !> Initialize the module
  subroutine gravity_init()
    use mod_global_parameters
    integer :: nwx,idir

    call grav_params_read(par_files)

    ! Determine flux variables
    nwx = 1                  ! rho (density)

    allocate(mom(ndir))
    do idir = 1, ndir
       nwx    = nwx + 1
       mom(idir) = nwx       ! momentum density
    end do

    nwx = nwx + 1
    e_     = nwx          ! energy density

  end subroutine gravity_init

  !> w[iw]=w[iw]+qdt*S[wCT,qtC,x] where S is the source based on wCT within ixO
  subroutine gravity_add_source(qdt,ixI^L,ixO^L,wCT,w,x,energy,qsourcesplit)
    use mod_global_parameters
    use mod_usr_methods

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, x(ixI^S,1:ndim)
    double precision, intent(in)    :: wCT(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    logical, intent(in) :: energy,qsourcesplit

    integer                         :: idim

    double precision :: gravity_field(ixI^S,ndim)

    if(qsourcesplit .eqv. grav_split) then
  
      if (.not. associated(usr_gravity)) then
        call mpistop("gravity_add_source: usr_gravity not defined")
      else
        call usr_gravity(ixI^L,ixO^L,wCT,x,gravity_field)
      end if
  
      do idim = 1, ndim
        w(ixO^S,mom(idim)) = w(ixO^S,mom(idim)) &
              + qdt * gravity_field(ixO^S,idim) * wCT(ixO^S,rho_)
        if(energy) then
          w(ixO^S,e_)=w(ixO^S,e_) &
              + qdt * gravity_field(ixO^S,idim) * wCT(ixO^S,mom(idim))
        end if
      end do
    end if

  end subroutine gravity_add_source

  subroutine gravity_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
    use mod_global_parameters
    use mod_usr_methods

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: dx^D, x(ixI^S,1:ndim), w(ixI^S,1:nw)
    double precision, intent(inout) :: dtnew

    double precision                :: dxinv(1:ndim)
    integer                         :: idim

    double precision :: gravity_field(ixI^S,ndim)

    ^D&dxinv(^D)=one/dx^D;

    if (.not. associated(usr_gravity)) then
      call mpistop("gravity_add_source: usr_gravity not defined")
    else
      call usr_gravity(ixI^L,ixO^L,w,x,gravity_field)
    end if

    do idim = 1, ndim
      dtnew = min(dtnew, minval(1.0d0 / &
           sqrt(abs(gravity_field(ixO^S,idim)) * dxinv(idim))))
    end do

  end subroutine gravity_get_dt

end module mod_gravity
