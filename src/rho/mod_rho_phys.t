!> Module containing the physics routines for scalar advection
module mod_rho_phys

  implicit none
  private

  integer, protected, public          :: rho_       = 1

  !> Whether particles module is added
  logical, public, protected              :: rho_particles = .false.

  double precision, protected, public :: rho_v(^ND) = 1.0d0

  ! Public methods
  public :: rho_phys_init
  public :: rho_get_v

contains

  subroutine rho_params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /rho_list/ rho_v, rho_particles

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status='old')
       read(unitpar, rho_list, end=111)
111    close(unitpar)
    end do

  end subroutine rho_params_read

  !> Write this module's parameters to a snapsoht
  subroutine rho_write_info(fh)
    use mod_global_parameters
    integer, intent(in)                 :: fh
    integer, parameter                  :: n_par = ^ND
    double precision                    :: values(n_par)
    character(len=name_len)             :: names(n_par)
    integer, dimension(MPI_STATUS_SIZE) :: st
    integer                             :: er
    integer                             :: idim

    call MPI_FILE_WRITE(fh, n_par, 1, MPI_INTEGER, st, er)

    do idim=1,ndim
      write(names(idim),'(a,i1)') "v",idim
      values(idim) = rho_v(idim)
    end do
    call MPI_FILE_WRITE(fh, values, n_par, MPI_DOUBLE_PRECISION, st, er)
    call MPI_FILE_WRITE(fh, names, n_par * name_len, MPI_CHARACTER, st, er)
  end subroutine rho_write_info

  subroutine rho_phys_init()
    use mod_global_parameters
    use mod_physics
    use mod_particles, only: particles_init

    call rho_params_read(par_files)

    physics_type = "rho"
    phys_energy  = .false.
    phys_req_diagonal = .false.
    use_particles = rho_particles

    allocate(start_indices(number_species),stop_indices(number_species))
    ! set the index of the first flux variable for species 1
    start_indices(1)=1
    rho_ = var_set_rho()

    ! set number of variables which need update ghostcells
    nwgc=nwflux

    ! set the index of the last flux variable for species 1
    stop_indices(1)=nwflux

    ! Check whether custom flux types have been defined
    if (.not. allocated(flux_type)) then
       allocate(flux_type(ndir, nw))
       flux_type = flux_default
    else if (any(shape(flux_type) /= [ndir, nw])) then
       call mpistop("phys_check error: flux_type has wrong shape")
    end if

    phys_get_cmax        => rho_get_cmax
    phys_get_cbounds     => rho_get_cbounds
    phys_get_flux        => rho_get_flux
    phys_get_v_idim      => rho_get_v_idim
    phys_add_source_geom => rho_add_source_geom
    phys_to_conserved    => rho_to_conserved
    phys_to_primitive    => rho_to_primitive
    phys_get_dt          => rho_get_dt
    phys_write_info      => rho_write_info

    ! Initialize particles module
    if (rho_particles) then
       call particles_init()
       phys_req_diagonal = .true.
    end if

  end subroutine rho_phys_init

  subroutine rho_to_conserved(ixI^L, ixO^L, w, x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:^ND)

    ! Do nothing (primitive and conservative are equal for rho module)
  end subroutine rho_to_conserved

  subroutine rho_to_primitive(ixI^L, ixO^L, w, x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:^ND)

    ! Do nothing (primitive and conservative are equal for rho module)
  end subroutine rho_to_primitive

  subroutine rho_get_v(w, x, ixI^L, ixO^L, idim, v, centered)
    use mod_global_parameters
    use mod_geometry
    logical, intent(in)           :: centered
    integer, intent(in)           :: ixI^L, ixO^L, idim
    double precision, intent(in)  :: w(ixI^S, nw), x(ixI^S, 1:^ND)
    double precision, intent(out) :: v(ixI^S)

    double precision :: dtheta, dphi, halfdtheta, halfdphi, invdtheta, invdphi
    {^IFTHREED
    double precision :: appcosphi(ixI^S), appsinphi(ixI^S), &
                        appcosthe(ixI^S), appsinthe(ixI^S)
    }

    select case (coordinate)
    case (cylindrical)
       {^IFONED 
       call mpistop("advection in 1D cylindrical not available")
       }
       {^IFTWOD 
       ! advection in 2D cylindrical: polar grid: to v_R, v_varphi
       if(centered)then
          select case (idim)
             case (1) ! radial velocity
                v(ixO^S) = rho_v(1)*dcos(x(ixO^S,2))+rho_v(2)*dsin(x(ixO^S,2))
             case (2) ! v_varphi
                v(ixO^S) =-rho_v(1)*dsin(x(ixO^S,2))+rho_v(2)*dcos(x(ixO^S,2))
          end select
       else
          ! assumed uniform in varphi=theta
          dtheta=x(ixOmin1,ixOmin2+1,2)-x(ixOmin1,ixOmin2,2)
          halfdtheta=0.5d0*dtheta
          invdtheta=1.0d0/dtheta
          select case (idim)
             case (1) ! radial velocity
                v(ixO^S) =( rho_v(1)*( dsin(x(ixO^S,2)+halfdtheta) &
                                      -dsin(x(ixO^S,2)-halfdtheta)) &
                           +rho_v(2)*(-dcos(x(ixO^S,2)+halfdtheta) &
                                      +dcos(x(ixO^S,2)-halfdtheta)))*invdtheta
             case (2) ! v_varphi
                v(ixO^S) =-rho_v(1)*dsin(x(ixO^S,2)+halfdtheta) &
                          +rho_v(2)*dcos(x(ixO^S,2)+halfdtheta)
          end select
       endif
       }
       {^IFTHREED
       ! advection in 3D cylindrical: convert to v_R, v_Z, v_varphi
       if(centered)then
          select case (idim)
             case (1) ! v_R velocity
                v(ixO^S) = rho_v(1)*dcos(x(ixO^S,3))+rho_v(2)*dsin(x(ixO^S,3))
             case (2) ! v_Z velocity
                v(ixO^S) = rho_v(3)
             case (3) ! v_varphi velocity
                v(ixO^S) =-rho_v(1)*dsin(x(ixO^S,3))+rho_v(2)*dcos(x(ixO^S,3))
          end select
       else
          ! assumed uniform in varphi=theta
          dtheta=x(ixOmin1,ixOmin2,ixOmin3+1,3)-x(ixOmin1,ixOmin2,ixOmin3,3)
          halfdtheta=0.5d0*dtheta
          invdtheta=1.0d0/dtheta
          select case (idim)
             case (1) ! radial velocity
                v(ixO^S) =( rho_v(1)*( dsin(x(ixO^S,3)+halfdtheta) &
                                      -dsin(x(ixO^S,3)-halfdtheta)) &
                           +rho_v(2)*(-dcos(x(ixO^S,3)+halfdtheta) &
                                      +dcos(x(ixO^S,3)-halfdtheta)))*invdtheta
             case (2) ! v_Z velocity
                v(ixO^S) = rho_v(3)
             case (3) ! v_varphi
                v(ixO^S) =-rho_v(1)*dsin(x(ixO^S,3)+halfdtheta) &
                          +rho_v(2)*dcos(x(ixO^S,3)+halfdtheta)
          end select
       endif
       }
    case (spherical)
       {^IFONED 
       call mpistop("advection in 1D spherical not available")
       }
       {^IFTWOD 
       call mpistop("advection in 2D spherical not available")
       }
       {^IFTHREED
       ! advection in 3D spherical: convert to v_r, v_theta, v_phi
       if(centered)then
          select case (idim)
             case (1) ! radial velocity
                v(ixO^S) = rho_v(1)*dsin(x(ixO^S,2))*dcos(x(ixO^S,3)) &
                          +rho_v(2)*dsin(x(ixO^S,2))*dsin(x(ixO^S,3)) &
                          +rho_v(3)*dcos(x(ixO^S,2)) 
             case (2) ! theta velocity
                v(ixO^S) = rho_v(1)*dcos(x(ixO^S,2))*dcos(x(ixO^S,3)) &
                          +rho_v(2)*dcos(x(ixO^S,2))*dsin(x(ixO^S,3)) &
                          -rho_v(3)*dsin(x(ixO^S,2)) 
             case (3) ! phi velocity
                v(ixO^S) =-rho_v(1)*dsin(x(ixO^S,3)) &
                          +rho_v(2)*dcos(x(ixO^S,3)) 
          end select
       else
          ! assumed uniform in theta and phi
          dtheta=x(ixOmin1,ixOmin2+1,ixOmin3,2)-x(ixOmin1,ixOmin2,ixOmin3,2)
          dphi=x(ixOmin1,ixOmin2,ixOmin3+1,3)-x(ixOmin1,ixOmin2,ixOmin3,3)
          halfdtheta=0.5d0*dtheta
          halfdphi=0.5d0*dphi
          invdtheta=1.0d0/dtheta
          invdphi=1.0d0/dphi
          select case (idim)
             case (1) ! radial velocity
                appcosphi(ixO^S)=( dsin(x(ixO^S,3)+halfdphi) &
                                  -dsin(x(ixO^S,3)-halfdphi))*invdphi 
                appsinphi(ixO^S)=(-dcos(x(ixO^S,3)+halfdphi) &
                                  +dcos(x(ixO^S,3)-halfdphi))*invdphi 
                appcosthe(ixO^S)=(dsin(x(ixO^S,2)+halfdtheta)**2 &
                                 -dsin(x(ixO^S,2)-halfdtheta)**2) &
                          /(4.0d0*dabs(dsin(x(ixO^S,2)))*dsin(halfdtheta))
                appsinthe(ixO^S)= &
                           (-dsin(x(ixO^S,2)+halfdtheta)*dcos(x(ixO^S,2)+halfdtheta)  & 
                            +dsin(x(ixO^S,2)-halfdtheta)*dcos(x(ixO^S,2)-halfdtheta)  & 
                            +dtheta)/(4.0d0*dabs(dsin(x(ixO^S,2)))*dsin(halfdtheta))
                v(ixO^S) = rho_v(1)*appsinthe(ixO^S)*appcosphi(ixO^S) &
                          +rho_v(2)*appsinthe(ixO^S)*appsinphi(ixO^S) &
                          +rho_v(3)*appcosthe(ixO^S) 
             case (2) ! theta velocity
                appcosphi(ixO^S)=( dsin(x(ixO^S,3)+halfdphi) &
                                  -dsin(x(ixO^S,3)-halfdphi))*invdphi 
                appsinphi(ixO^S)=(-dcos(x(ixO^S,3)+halfdphi) &
                                  +dcos(x(ixO^S,3)-halfdphi))*invdphi 
                v(ixO^S) = rho_v(1)*dcos(x(ixO^S,2)+halfdtheta)*appcosphi(ixO^S) &
                          +rho_v(2)*dcos(x(ixO^S,2)+halfdtheta)*appsinphi(ixO^S) &
                          -rho_v(3)*dsin(x(ixO^S,2)+halfdtheta) 
             case (3) ! phi velocity
                v(ixO^S) =-rho_v(1)*dsin(x(ixO^S,3)+halfdphi) &
                          +rho_v(2)*dcos(x(ixO^S,3)+halfdphi) 
          end select
       endif
       }
    case default
       v(ixO^S) = rho_v(idim)
    end select
  end subroutine rho_get_v

  !> Calculate simple v component
  subroutine rho_get_v_idim(w,x,ixI^L,ixO^L,idim,v)
    use mod_global_parameters

    integer, intent(in)           :: ixI^L, ixO^L, idim
    double precision, intent(in)  :: w(ixI^S,nw), x(ixI^S,1:ndim)
    double precision, intent(out) :: v(ixI^S)

    v(ixO^S) = rho_v(idim)

  end subroutine rho_get_v_idim

  subroutine rho_get_cmax(w, x, ixI^L, ixO^L, idim, cmax)
    use mod_global_parameters
    integer, intent(in)                       :: ixI^L, ixO^L, idim
    double precision, intent(in)              :: w(ixI^S, nw), x(ixI^S, 1:^ND)
    double precision, intent(inout)           :: cmax(ixI^S)

    call rho_get_v(w, x, ixI^L, ixO^L, idim, cmax, .true.)

    cmax(ixO^S) = abs(cmax(ixO^S))

  end subroutine rho_get_cmax

  subroutine rho_get_cbounds(wLC, wRC, wLp, wRp, x, ixI^L, ixO^L, idim,Hspeed, cmax, cmin)
    use mod_global_parameters
    use mod_variables
    integer, intent(in)             :: ixI^L, ixO^L, idim
    double precision, intent(in)    :: wLC(ixI^S, nw), wRC(ixI^S,nw)
    double precision, intent(in)    :: wLp(ixI^S, nw), wRp(ixI^S,nw)
    double precision, intent(in)    :: x(ixI^S, 1:^ND)
    double precision, intent(inout) :: cmax(ixI^S,1:number_species)
    double precision, intent(inout), optional :: cmin(ixI^S,1:number_species)
    double precision, intent(in)    :: Hspeed(ixI^S)

    ! If get_v depends on w, the first argument should be some average over the
    ! left and right state
    call rho_get_v(wLC, x, ixI^L, ixO^L, idim, cmax(ixI^S,1), .false.)

    if (present(cmin)) then
       cmin(ixO^S,1) = min(cmax(ixO^S,1), zero)
       cmax(ixO^S,1) = max(cmax(ixO^S,1), zero)
    else
       cmax(ixO^S,1) = maxval(abs(cmax(ixO^S,1)))
    end if

  end subroutine rho_get_cbounds

  subroutine rho_get_dt(w, ixI^L, ixO^L, dtnew, dx^D, x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: dx^D, x(ixI^S, 1:^ND)
    double precision, intent(in)    :: w(ixI^S, 1:nw)
    double precision, intent(inout) :: dtnew

    dtnew = bigdouble
  end subroutine rho_get_dt

  ! There is nothing to add to the transport flux in the transport equation
  subroutine rho_get_flux(wC, w, x, ixI^L, ixO^L, idim, f)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L, idim
    double precision, intent(in)    :: wC(ixI^S, 1:nw)
    double precision, intent(in)    :: w(ixI^S, 1:nw)
    double precision, intent(in)    :: x(ixI^S, 1:^ND)
    double precision, intent(out)   :: f(ixI^S, nwflux)
    double precision                :: v(ixI^S)

    call rho_get_v(wC, x, ixI^L, ixO^L, idim, v, .false.)

    f(ixO^S, rho_) = w(ixO^S, rho_) * v(ixO^S)
  end subroutine rho_get_flux

  subroutine rho_add_source_geom(qdt, ixI^L, ixO^L, wCT, w, x)

    ! Add geometrical source terms to w
    ! There are no geometrical source terms in the transport equation

    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: qdt, x(ixI^S, 1:^ND)
    double precision, intent(inout) :: wCT(ixI^S, 1:nw), w(ixI^S, 1:nw)

  end subroutine rho_add_source_geom

end module mod_rho_phys
