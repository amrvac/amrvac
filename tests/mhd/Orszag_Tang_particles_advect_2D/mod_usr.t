module mod_usr
  use mod_mhd
  implicit none

contains

  subroutine usr_init()
    unit_length=1.d9
    unit_numberdensity=1.d9
    unit_velocity=1.d7
    usr_init_one_grid => initonegrid_usr
    usr_aux_output    => specialvar_output
    usr_add_aux_names => specialvarnames_output 

    call set_coordinate_system("Cartesian_2D")

    call mhd_activate()

  end subroutine usr_init

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    ! initialize one grid
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision:: rho0,p0,b0
    logical, save :: first=.true.

    rho0=25.0d0/36.0d0/dpi
    p0=5.0d0/12.0d0/dpi
    b0=1.0d0/sqrt(4.0d0*dpi)

    w(ixO^S,rho_)=rho0
    w(ixO^S,mom(1))=-sin(2.0d0*dpi*x(ixO^S,2))
    w(ixO^S,mom(2))= sin(2.0d0*dpi*x(ixO^S,1))
    w(ixO^S,p_)=p0
    w(ixO^S,mag(1))=-b0*sin(2.0d0*dpi*x(ixO^S,2))
    w(ixO^S,mag(2))= b0*sin(4.0d0*dpi*x(ixO^S,1))

    call mhd_to_conserved(ixI^L,ixO^L,w,x)

    if(first .and. mype==0 )then
      write(*,*)'Doing 2.5D ideal MHD, Orszag Tang problem with tracing particles'
      write(*,*)'rho - p - b:',rho0,p0,b0
      first=.false.
    endif

  end subroutine initonegrid_usr

  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
    ! this subroutine can be used in convert, to add auxiliary variables to the
    ! converted output file, for further analysis using tecplot, paraview, ....
    ! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
    !
    ! the array normconv can be filled in the (nw+1:nw+nwauxio) range with 
    ! corresponding normalization values (default value 1)
    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)

    double precision:: divb(ixI^S)
    double precision :: current(ixI^S,7-2*ndir:3)
    integer          :: idirmin
    
    call get_divb(w,ixI^L,ixO^L,divb)
    w(ixO^S,nw+1)=divb(ixO^S)
    
    call get_current(w,ixI^L,ixO^L,idirmin,current)
    w(ixO^S,nw+2)=current(ixO^S,3)

  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
    character(len=*) :: varnames
    varnames='divb jz'

  end subroutine specialvarnames_output

end module mod_usr
