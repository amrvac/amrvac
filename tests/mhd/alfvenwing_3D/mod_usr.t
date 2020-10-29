module mod_usr

  ! Include a physics module
  use mod_mhd

  implicit none
  ! Input values
  double precision :: velocity_x, velocity_y, velocity_z, mag_x, mag_y, mag_z, Collfreq, pressure, density, radius
  double precision :: Vnorm, Bnorm, Collfnorm, Pnorm, Dnorm
contains

  !> Read this module's parameters from a file
  subroutine usr_params_read(files)
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /usr_list/  velocity_x,velocity_y,velocity_z,mag_x,mag_y,mag_z, Collfreq, pressure, density, radius
! Values are in SI units
    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, usr_list, end=111)
111    close(unitpar)
    end do

 !Normalization values
  Vnorm     = sqrt(velocity_x**2+velocity_y**2*velocity_z**2) !normalization value for the plasma velocity
  Bnorm     = sqrt(miu0_si*density*Vnorm**2) !normalization value for the magnetic field
  Collfnorm = Vnorm/radius !normalization value for the ion-neutral collision frequency
  Pnorm     = density*Vnorm**2 !normalization value for the thermal plasma pressure
  Dnorm     = density !normalization value for the plasma density

  end subroutine usr_params_read


  subroutine usr_init()
    call usr_params_read(par_files)
    call set_coordinate_system("Cartesian_3D")

    usr_init_one_grid => initonegrid_usr
    usr_special_bc => specialbound_usr
    usr_source => specialsource
    usr_refine_grid => specialrefine_grid
    usr_get_dt        => get_dt_ionosphere

    call mhd_activate()
  end subroutine usr_init

  subroutine initglobaldata_usr

    use mod_global_parameters

    mhd_gamma=5.0d0/3.0d0
    mhd_eta=zero

  end subroutine initglobaldata_usr

  subroutine initonegrid_usr(ixG^L,ix^L,w,x)

    ! initialize one grid

    use mod_global_parameters

    integer, intent(in) :: ixG^L,ix^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)
    double precision:: ss(ixG^S), HH(ixG^S), Req(ixG^S), rhoeq(ixG^S), peq(ixG^S)

    logical :: patchw(ixG^T)
    logical, save :: first=.true.

    {^IFONED call mpistop("This is a 3D MHD problem")}
    {^IFTWOD call mpistop("This is a 3D MHD problem")}

    ! values are dimensionless    
    w(ix^S,rho_)=density/dnorm
    w(ix^S,mom(1))=velocity_x/Vnorm
    w(ix^S,mom(2))=velocity_y/Vnorm
    w(ix^S,mom(3))=velocity_z/Vnorm
    w(ix^S,p_)=pressure/Pnorm
    w(ix^S,mag(1))=mag_x/Bnorm
    w(ix^S,mag(2))=mag_y/Bnorm
    w(ix^S,mag(3))=mag_z/Bnorm

    if(first .and. mype==0)then
        write(*,*)'Doing Alfven wing test, 3D ideal MHD'
        write(*,*)'velocity (vx,vy,vz) in m/s:',velocity_x,velocity_y,velocity_z
        write(*,*)'velocity (vx,vy,vz) normalized:',velocity_x/Vnorm,velocity_y/Vnorm,velocity_z/Vnorm
        write(*,*)'magnetic field (Bx,By,Bz) in nT:',mag_x,mag_y,mag_z
        write(*,*)'magnetic field (Bx,By,Bz) normalized:',mag_x/Bnorm,mag_y/Bnorm,mag_z/Bnorm
        write(*,*)'Collfreq (Hz), pressure (Pa), density (kg/m^3):',Collfreq, pressure, density
        write(*,*)'Collfreq, pressure, density normalized:',Collfreq/Collfnorm, pressure/Pnorm, density/dnorm
        first=.false.
    endif

    ! now let the code compute the conservative variables itself
    call mhd_to_conserved(ixG^L,ix^L, w, x)

  end subroutine initonegrid_usr

  ! special boundary types, user defined
  subroutine specialbound_usr(qt,ixG^L,ixO^L,iB,w,x)
    use mod_global_parameters

    integer, intent(in) :: ixO^L, iB, ixG^L
    integer:: ix^D
    double precision, intent(in) :: qt, x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)
    logical :: patchw(ixG^T)


    select case(iB)

    case(1)

      w(ixO^S,rho_)=density/dnorm
      w(ixO^S,p_)=pressure/Pnorm 
      w(ixO^S,mom(1))=velocity_x/Vnorm
      w(ixO^S,mom(2))=velocity_y/Vnorm
      w(ixO^S,mom(3))=velocity_z/Vnorm
      w(ixO^S,mag(1))=mag_x/Bnorm
      w(ixO^S,mag(2))=mag_y/Bnorm
      w(ixO^S,mag(3))=mag_z/Bnorm
      if(mhd_glm) w(ixO^S,psi_)=0.d0
      call mhd_to_conserved(ixG^L,ixO^L,w,x)
    case default
       call mpistop("Special boundary is not defined for this region")
    end select

end subroutine specialbound_usr


subroutine specialsource(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L, iw^LIM
    double precision, intent(in) :: qdt, qtC, qt
    double precision, intent(in) :: x(ixI^S,1:ndim), wCT(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: r(ixO^S),collfreq_n

       collfreq_n=Collfreq/Collfnorm

      r(ixO^S)=sqrt(x(ixO^S,1)**2+x(ixO^S,2)**2+x(ixO^S,3)**2)

     where(r(ixO^S)<1.0d0)
      w(ixO^S,mom(1))=w(ixO^S,mom(1))-qdt*collfreq_n*wCT(ixO^S,mom(1))
      w(ixO^S,mom(2))=w(ixO^S,mom(2))-qdt*collfreq_n*wCT(ixO^S,mom(2))
      w(ixO^S,mom(3))=w(ixO^S,mom(3))-qdt*collfreq_n*wCT(ixO^S,mom(3))

      w(ixO^S,e_)=w(ixO^S,e_)+qdt*0.5d0*collfreq_n/wCT(ixO^S,rho_)*(-wCT(ixO^S,mom(1))*wCT(ixO^S,mom(1))& 
                                                                  -wCT(ixO^S,mom(2))*wCT(ixO^S,mom(2))&
                                                                  -wCT(ixO^S,mom(3))*wCT(ixO^S,mom(3)) ) 
     end where

  end subroutine specialsource


  subroutine get_dt_ionosphere(w,ixG^L,ix^L,dtnew,dx^D,x)

    use mod_global_parameters

    integer, intent(in) :: ixG^L, ix^L
    double precision, intent(in) :: dx^D, x(ixG^S,1:ndim), w(ixG^S,1:nw)
    double precision, intent(inout) :: dtnew
    double precision:: collfreq_n

    collfreq_n=Collfreq/Collfnorm
    dtnew=bigdouble
    dtnew=min(dtnew,dtdiffpar*one/collfreq_n)

  end subroutine get_dt_ionosphere



  subroutine specialrefine_grid(igrid,level,ixG^L,ix^L,qt,w,x,refine,coarsen)

    ! Enforce additional refinement or coarsening
    ! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.

    use mod_global_parameters

    integer, intent(in) :: igrid, level, ixG^L, ix^L
    double precision, intent(in) :: qt, w(ixG^S,1:nw), x(ixG^S,1:ndim)
    integer, intent(inout) :: refine, coarsen

    if (any(abs(x(ix^S,1)) < 2.5d0) .and. any(abs(x(ix^S,2)) < 2.5d0) .and. any(abs(x(ix^S,3)) < 6d0)) then 
       refine=1 
    else
       coarsen=1
    endif

  end subroutine specialrefine_grid

end module mod_usr
