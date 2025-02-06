module mod_usr
  use mod_hd
  use mod_viscosity
  implicit none

  ! input control values: Reynolds and Mach number, derived pressure
  double precision :: Re,Ma,p0,RR

contains

  subroutine usr_init()
    use mod_variables

    call usr_params_read(par_files)

    usr_set_parameters => initglobaldata_usr
    usr_modify_output  => set_internal_cylinder
    usr_init_one_grid  => initonegrid_usr
    usr_special_bc     => specialbound_usr
    usr_internal_bc    => no_vel
    usr_refine_grid    => specialrefine_grid
    usr_aux_output     => specialvar_output
    usr_add_aux_names  => specialvarnames_output

    call set_coordinate_system("Cartesian")

    call hd_activate()

  end subroutine usr_init

  !> Read this module's parameters from a file
  subroutine usr_params_read(files)
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /usr_list/ Re, Ma

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, usr_list, end=111)
111    close(unitpar)
    end do

  end subroutine usr_params_read

  subroutine initglobaldata_usr
    character(len=20) :: printsettingformat

    printsettingformat='(1x,A50,ES15.7,A7)'
    
    vc_mu=1.0/Re
    p0=1.0d0/(hd_gamma*Ma**2)
    RR=0.5d0

    if(mype==0) then
      write(*,*) "Karman street setup:"
      write(*,printsettingformat) "Mach number ",Ma," input"
      write(*,printsettingformat) "Reynolds number ",Re," input"
      write(*,printsettingformat) "viscosity coefficient ",vc_mu," derived as 1/Re"
      write(*,printsettingformat) "gamma ",hd_gamma," input"
      write(*,printsettingformat) "pressure ",p0," derived as 1/(gamma M2)"
    end if

  end subroutine initglobaldata_usr

  !> Initialize one grid
  subroutine initonegrid_usr(ixG^L,ix^L,w,x)
    integer, intent(in)             :: ixG^L, ix^L
    double precision, intent(in)    :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    double precision :: rad(ixG^S),costhe(ixG^S),cos2theta(ixG^S),rad2(ixG^S),RR2
    integer :: idims
    logical                         :: first = .true.

    if (first) then
       if (mype==0) then
          print *,'A flow around an infinitely long cylinder'
          print *,'Re=', Re
          print *,'Ma=', Ma
          print *,'Initial flow is potential'
          print *,'radius of cylinder is ', RR
       end if
       first=.false.
    end if

    w(ix^S,rho_)  =1.0d0
    RR2=RR**2
    rad2(ix^S)=x(ix^S,1)**2+x(ix^S,2)**2
    rad(ix^S)=dsqrt(x(ix^S,1)**2+x(ix^S,2)**2)
    costhe(ix^S)=x(ix^S,1)/rad(ix^S)
    cos2theta(ix^S)=2.0d0*costhe(ix^S)**2-1.0d0
    w(ix^S,mom(1))=1.0d0+RR2/rad2(ix^S)-2.0d0*x(ix^S,1)**2*RR2/rad2(ix^S)**2
    w(ix^S,mom(2))=-2.0d0*x(ix^S,1)*x(ix^S,2)*RR2/rad2(ix^S)**2
    w(ix^S,p_)=0.5d0*(2.0d0*cos2theta(ix^S)*RR2/rad2(ix^S)-RR2**2/rad2(ix^S)**2)+p0

    where(rad(ix^S)<RR)
      w(ix^S,mom(1))=0.0d0
      w(ix^S,mom(2))=0.0d0
      w(ix^S,rho_)  =1.0d0
      w(ix^S,p_)    =p0
    endwhere

    call hd_to_conserved(ixG^L,ix^L,w,x)

  end subroutine initonegrid_usr

  subroutine specialbound_usr(qt,ixI^L,ixO^L,iB,w,x)
    ! special boundary types, user defined
    integer, intent(in) :: ixO^L, iB, ixI^L
    double precision, intent(in) :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: rad(ixI^S),costhe(ixI^S),cos2theta(ixI^S),rad2(ixI^S),RR2

    select case(iB)
    ! special left boundary
    ! fix inflow properties
    case(1)
       RR2=RR**2
       rad2(ixO^S)=x(ixO^S,1)**2+x(ixO^S,2)**2
       rad(ixO^S)=dsqrt(x(ixO^S,1)**2+x(ixO^S,2)**2)
       costhe(ixO^S)=x(ixO^S,1)/rad(ixO^S)
       cos2theta(ixO^S)=2.0d0*costhe(ixO^S)**2-1.0d0
       w(ixO^S,mom(1))=1.0d0+RR2/rad2(ixO^S)-2.0d0*x(ixO^S,1)**2*RR2/rad2(ixO^S)**2
       w(ixO^S,mom(2))=-2.0d0*x(ixO^S,1)*x(ixO^S,2)*RR2/rad2(ixO^S)**2
       w(ixO^S,p_)=0.5d0*(2.0d0*cos2theta(ixO^S)*RR2/rad2(ixO^S)-RR2**2/rad2(ixO^S)**2)+p0
       w(ixO^S,rho_)   = 1.0d0
       !!w(ixO^S,p_)     = p0
       !!w(ixO^S,mom(1)) = 1.0d0
       !!w(ixO^S,mom(2)) = zero
       call hd_to_conserved(ixI^L,ixO^L,w,x)
    case default
       call mpistop("Special boundary is not defined for this region")
    end select

  end subroutine specialbound_usr

  subroutine specialrefine_grid(igrid,level,ixG^L,ix^L,qt,w,x,refine,coarsen)
    ! Enforce additional refinement or coarsening
    ! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.
    ! you must set consistent values for integers refine/coarsen:
    ! refine = -1 enforce to not refine
    ! refine =  0 doesn't enforce anything
    ! refine =  1 enforce refinement
    ! coarsen = -1 enforce to not coarsen
    ! coarsen =  0 doesn't enforce anything
    ! coarsen =  1 enforce coarsen
    integer, intent(in) :: igrid, level, ixG^L, ix^L
    double precision, intent(in) :: qt, w(ixG^S,1:nw), x(ixG^S,1:ndim)
    integer, intent(inout) :: refine, coarsen
    double precision:: R(ixG^S)

    R(ix^S)=dsqrt(x(ix^S,1)**2+x(ix^S,2)**2)

    if (any(R(ix^S) <= 2.0d0*RR)) refine=1

  end subroutine specialrefine_grid

  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
  ! this subroutine can be used in convert, to add auxiliary variables to the
  ! converted output file, for further analysis using tecplot, paraview, ....
  ! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
  ! the array normconv can be filled in the (nw+1:nw+nwauxio) range with
  ! corresponding normalization values (default value 1)
    use mod_radiative_cooling
    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)

    double precision :: pth(ixI^S),gradrho(ixI^S),drho(ixI^S),vrot(ixI^S),tmp(ixI^S)
    double precision                   :: kk,grhomax,kk1
    double precision :: wlocal(ixI^S,1:nw)
    integer                            :: idims

    wlocal(ixI^S,1:nw)=w(ixI^S,1:nw)
    ! output temperature
    call hd_get_pthermal(wlocal,x,ixI^L,ixO^L,pth)
    w(ixO^S,nw+1)=pth(ixO^S)/w(ixO^S,rho_)

    ! output Mach number V/c_s
    w(ixO^S,nw+2)=dsqrt(wlocal(ixO^S,mom(1))**2+wlocal(ixO^S,mom(2))**2) &
                  /dsqrt(hd_gamma*pth(ixO^S)*w(ixO^S,rho_))

    ! output vorticity
    vrot(ixO^S)=zero
    idims=1
    tmp(ixI^S)=wlocal(ixI^S,mom(2))/wlocal(ixI^S,rho_)
    call gradient(tmp,ixI^L,ixO^L,idims,drho)
    vrot(ixO^S)=vrot(ixO^S)+drho(ixO^S)
    idims=2
    tmp(ixI^S)=wlocal(ixI^S,mom(1))/wlocal(ixI^S,rho_)
    call gradient(tmp,ixI^L,ixO^L,idims,drho)
    vrot(ixO^S)=vrot(ixO^S)-drho(ixO^S)
    w(ixO^S,nw+3)=vrot(ixO^S)

    ! output schlieren plot
    gradrho(ixO^S)=zero
    do idims=1,ndim
       select case(typegrad)
          case("central")
           call gradient(wlocal(ixI^S,rho_),ixI^L,ixO^L,idims,drho)
          case("limited")
           call gradientL(wlocal(ixI^S,rho_),ixI^L,ixO^L,idims,drho)
       end select
       gradrho(ixO^S)=gradrho(ixO^S)+drho(ixO^S)**2.0d0
    enddo
    gradrho(ixO^S)=dsqrt(gradrho(ixO^S))
    kk=5.0d0
    kk1=0.0001d0
    ! need the global maximum here, otherwise see the block structure reflected...
    !grhomax=max(1000.0d0,maxval(gradrho(ixO^S)))
    grhomax=100.0d0
    w(ixO^S,nw+4)=dexp(-kk*(gradrho(ixO^S)/grhomax-kk1))

  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
    ! newly added variables need to be concatenated with the w_names/primnames string
    character(len=*) :: varnames

    varnames='Te Mach omega Schlier'

  end subroutine specialvarnames_output

  subroutine set_internal_cylinder(ixI^L,ixO^L,qt,w,x)
    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: qt
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision :: rad(ixI^S)

    rad(ixO^S)=dsqrt(x(ixO^S,1)**2+x(ixO^S,2)**2)
    where (rad(ixO^S)<RR)
        w(ixO^S,mom(1)) = 0.d0
        w(ixO^S,mom(2)) = 0.d0
        w(ixO^S,p_)     = p0/(hd_gamma-1.0d0)
        w(ixO^S,rho_)   = 1.d0
    end where
  end subroutine set_internal_cylinder

  subroutine no_vel(level,qt,ixI^L,ixO^L,w,x)
    integer, intent(in) :: ixI^L,ixO^L,level
    double precision, intent(in) :: qt
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision :: rad2(ixI^S)
    rad2(ixO^S)=x(ixO^S,1)**2+x(ixO^S,2)**2
    where (rad2(ixO^S)<RR**2)
        w(ixO^S,mom(1)) = 0.d0
        w(ixO^S,mom(2)) = 0.d0
    end where
  end subroutine no_vel

end module mod_usr
