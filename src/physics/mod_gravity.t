!> Module for including gravity in (magneto)hydrodynamics simulations
module mod_gravity
  implicit none

  !> source split or not
  logical :: grav_split= .false.

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
    use mod_usr_methods
    use mod_comm_lib, only: mpistop
    integer :: nwx,idir

    if (.not. associated(usr_gravity)) then
      write(*,*) "mod_usr.t: please define usr_gravity before (m)hd activate"
      write(*,*) "like the phys_gravity in mod_usr_methods.t"
      call mpistop("usr_gravity not defined or pointed")
    end if
    call grav_params_read(par_files)
    if(grav_split) any_source_split=.true.
  end subroutine gravity_init

  !> w[iw]=w[iw]+qdt*S[wCT,qtC,x] where S is the source based on wCT within ixO
  subroutine gravity_add_source(qdt,ixI^L,ixO^L,wCT,wCTprim,w,x,&
       energy,rhov,qsourcesplit,active)
    use mod_global_parameters
    use mod_usr_methods
    use mod_comm_lib, only: mpistop
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, x(ixI^S,1:ndim)
    double precision, intent(in)    :: wCT(ixI^S,1:nw),wCTprim(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    logical, intent(in)             :: energy,rhov,qsourcesplit
    logical, intent(inout)          :: active

    double precision                :: gravity_field(ixI^S,ndim)
    integer                         :: idim

    if(qsourcesplit .eqv. grav_split) then
      active = .true.
      call usr_gravity(ixI^L,ixO^L,wCT,x,gravity_field)
      do idim = 1, ndim
        if(size(iw_mom)<ndim) then
          w(ixO^S,iw_mom(1)) = w(ixO^S,iw_mom(1)) &
                + qdt*gravity_field(ixO^S,idim)*wCT(ixO^S,iw_rho)*block%B0(ixO^S,idim,0)
        else
          w(ixO^S,iw_mom(idim)) = w(ixO^S,iw_mom(idim)) &
                + qdt*gravity_field(ixO^S,idim)*wCT(ixO^S,iw_rho)
        endif
        if(energy) then
          if(rhov) then
            if(size(iw_mom)<ndim) call mpistop("rhov is not supported for ffHD")
            w(ixO^S,iw_e)=w(ixO^S,iw_e) &
              +qdt*gravity_field(ixO^S,idim)*wCTprim(ixO^S,iw_mom(idim))*wCTprim(ixO^S,iw_rho)
          else
            if(size(iw_mom)<ndim) then
              w(ixO^S,iw_e)=w(ixO^S,iw_e) &
                + qdt*gravity_field(ixO^S,idim)*wCT(ixO^S,iw_mom(1))*block%B0(ixO^S,idim,0)
            else
              w(ixO^S,iw_e)=w(ixO^S,iw_e) &
                + qdt*gravity_field(ixO^S,idim)*wCT(ixO^S,iw_mom(idim))
            endif
          end if
        end if
      end do
    end if
  end subroutine gravity_add_source

  subroutine gravity_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
    use mod_global_parameters
    use mod_usr_methods
    use mod_geometry

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: dx^D, x(ixI^S,1:ndim), w(ixI^S,1:nw)
    double precision, intent(inout) :: dtnew

    double precision                :: dxinv(1:ndim), max_grav
    double precision :: gravity_field(ixI^S,ndim)
    integer                         :: idim

    call usr_gravity(ixI^L,ixO^L,w,x,gravity_field)
    if(slab_uniform) then
      ^D&dxinv(^D)=one/dx^D;
      do idim = 1, ndim
        max_grav = maxval(abs(gravity_field(ixO^S,idim)))
        max_grav = max(max_grav, epsilon(1.0d0))
        dtnew = min(dtnew, 1.0d0 / sqrt(max_grav * dxinv(idim)))
      end do
    else
      do idim = 1, ndim
        max_grav = maxval(abs(gravity_field(ixO^S,idim))/block%ds(ixO^S,idim))
        max_grav = max(max_grav, epsilon(1.0d0))
        dtnew = min(dtnew, 1.0d0 / sqrt(max_grav))
      end do
    endif
  end subroutine gravity_get_dt
end module mod_gravity
