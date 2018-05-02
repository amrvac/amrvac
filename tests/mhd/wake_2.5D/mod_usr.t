module mod_usr
  use mod_mhd
  implicit none

contains

  subroutine usr_init()
    usr_init_one_grid => initonegrid_usr
    usr_aux_output    => specialvar_output
    usr_add_aux_names => specialvarnames_output 

    call set_coordinate_system("Cartesian_2.5D")

    call mhd_activate()

  end subroutine usr_init
  subroutine initonegrid_usr(ixG^L,ix^L,w,x)
    ! initialize one grid
    integer, intent(in) :: ixG^L, ix^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    double precision:: dA,dM,width,dv,sigma,alfa,Lx,rho0,p0
    logical, save :: first=.true.

    dA=0.2d0
    dM=3.0d0
    width=one
    dv=0.01d0
    sigma=one
    alfa=0.35d0
    Lx=two*dpi/alfa
    rho0=one
    w(ix^S,rho_)=rho0
    if(first .and. mype==0)then
       write(*,*)'Doing 2.5D resistive MHD, Wake problem'
       write(*,*)'A  - M - gamma:',dA,dM,mhd_gamma
       write(*,*)'dv sigma:',dv,sigma
       write(*,*)'alfa Lx:',alfa,Lx
       write(*,*)'Assuming eta set, using value:',mhd_eta
       first=.false.
    endif
    w(ix^S,mom(1))=rho0*(one-one/dcosh(x(ix^S,2)))
    w(ix^S,mom(2))=rho0*(dv*(dsin(two*dpi*x(ix^S,1)/Lx)*dtanh(x(ix^S,2))&
                +dcos(two*dpi*x(ix^S,1)/Lx)/dcosh(x(ix^S,2)))&
                *dexp(-(x(ix^S,2)/sigma)**2))
    w(ix^S,mag(1))=dA*dtanh(x(ix^S,2)/width)
    w(ix^S,mag(2))=zero
    w(ix^S,mom(3))=zero
    w(ix^S,mag(3))=dA/dcosh(x(ix^S,2)/width)
    p0=one/(dM**2)/mhd_gamma
    w(ix^S,e_)=p0/(mhd_gamma-1)+&
      half*(sum(w(ix^S,mom(:))**2,dim=ndim+1)/rho0+sum(w(ix^S,mag(:))**2,dim=ndim+1))

  end subroutine initonegrid_usr

  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
  integer, intent(in)                :: ixI^L,ixO^L
  double precision, intent(in)       :: x(ixI^S,1:ndim)
  double precision                   :: w(ixI^S,nw+nwauxio)
  double precision                   :: normconv(0:nw+nwauxio)

  double precision :: current(ixI^S,7-2*ndir:3)
  double precision :: wloc(ixI^S,1:nw)
  double precision :: pth(ixI^S)
  integer          :: idirmin

  wloc(ixI^S,1:nw)=w(ixI^S,1:nw)
  call mhd_get_pthermal(wloc,x,ixI^L,ixO^L,pth)
  w(ixO^S,nw+1)=pth(ixO^S)/wloc(ixO^S,rho_)
  
  call get_current(wloc,ixI^L,ixO^L,idirmin,current)
  w(ixO^S,nw+2)=current(ixO^S,3)
  
  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
    character(len=*) :: varnames
    varnames='T jz'

  end subroutine specialvarnames_output

end module mod_usr
