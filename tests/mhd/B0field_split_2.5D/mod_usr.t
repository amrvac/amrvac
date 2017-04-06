module mod_usr
  use mod_mhd
  implicit none
  double precision :: theta, kx, ly, vc

contains

  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods

    call set_coordinate_system("Cartesian_2.5D")

    usr_set_parameters  => initglobaldata_usr
    usr_init_one_grid   => initonegrid_usr
    usr_special_bc      => specialbound_usr
    usr_set_B0          => specialset_B0
    usr_set_J0          => specialset_J0
    usr_aux_output      => specialvar_output
    usr_add_aux_names   => specialvarnames_output 

    call mhd_activate()
  end subroutine usr_init

  subroutine initglobaldata_usr()
    use mod_global_parameters

    theta=60.d0*dpi/180.d0 ! the angle to the plane xy, 90-theta is the angle to the polarity inversion line of the arcade 
    kx=2.d0*dpi/(xprobmax1-xprobmin1)
    ly=kx*dcos(theta)
    if(iprob==1) then
      vc=0.d0
    else
      vc=(xprobmax1-xprobmin1)/2.d0
    end if

  end subroutine initglobaldata_usr

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    ! initialize one grid
    use mod_global_parameters

    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: res
    integer :: ix^D,na
    logical, save :: first=.true.

    if(first)then
      if(mype==0) then
        write(*,*)'Simulating 2.5D static magnetic arcade'
      endif
      first=.false.
    endif
    w(ixO^S,rho_)=1.d0
    w(ixO^S,p_)=1.d0
    w(ixO^S,mom(:))=zero
    w(ixO^S,mom(1))=vc
    if(B0field) then
      w(ixO^S,mag(:))=zero
    else
      w(ixO^S,mag(1))=-Busr*dcos(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))*dcos(theta)
      w(ixO^S,mag(2))= Busr*dsin(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))
      w(ixO^S,mag(3))=-Busr*dcos(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))*dsin(theta)
    endif
    call mhd_to_conserved(ixI^L,ixO^L,w,x)

  end subroutine initonegrid_usr

  subroutine specialbound_usr(qt,ixI^L,ixO^L,iB,w,x)
    ! special boundary types, user defined
    use mod_global_parameters

    integer, intent(in) :: ixO^L, iB, ixI^L
    double precision, intent(in) :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: pth(ixI^S),tmp(ixI^S),ggrid(ixI^S),invT(ixI^S)
    double precision :: delydelx
    integer :: ix^D,idir,ixInt^L

    select case(iB)
    case(3)
      w(ixO^S,rho_)=1.d0
      w(ixO^S,p_)=1.d0
      !! fixed zero velocity
      do idir=2,ndir
        w(ixO^S,mom(idir)) =-w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,mom(idir))&
                   /w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,rho_)
      end do
      w(ixO^S,mom(1))=vc
      w(ixO^S,mom(2))=0.d0
      w(ixO^S,mom(3))=0.d0
      if(B0field) then
        w(ixO^S,mag(1))=-Busr*dcos(kx*(x(ixO^S,1)-vc*qt))*dexp(-ly*x(ixO^S,2))*dcos(theta)
        w(ixO^S,mag(2))= Busr*dsin(kx*(x(ixO^S,1)-vc*qt))*dexp(-ly*x(ixO^S,2))
        w(ixO^S,mag(3))=-Busr*dcos(kx*(x(ixO^S,1)-vc*qt))*dexp(-ly*x(ixO^S,2))*dsin(theta)
        w(ixO^S,mag(:))=w(ixO^S,mag(:))-block%B0(ixO^S,:,0)
      else
        w(ixO^S,mag(1))=-Busr*dcos(kx*(x(ixO^S,1)-vc*qt))*dexp(-ly*x(ixO^S,2))*dcos(theta)
        w(ixO^S,mag(2))= Busr*dsin(kx*(x(ixO^S,1)-vc*qt))*dexp(-ly*x(ixO^S,2))
        w(ixO^S,mag(3))=-Busr*dcos(kx*(x(ixO^S,1)-vc*qt))*dexp(-ly*x(ixO^S,2))*dsin(theta)
      endif
      call mhd_to_conserved(ixI^L,ixO^L,w,x)
    case(4)
      w(ixO^S,rho_)=1.d0
      w(ixO^S,p_)=1.d0
      ! fixed zero velocity
      do idir=2,ndir
        w(ixO^S,mom(idir)) =-w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-nghostcells:-1,mom(idir))&
                     /w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-nghostcells:-1,rho_)
      end do
      w(ixO^S,mom(1))=vc
      w(ixO^S,mom(2))=0.d0
      w(ixO^S,mom(3))=0.d0
      if(B0field) then
        w(ixO^S,mag(1))=-Busr*dcos(kx*(x(ixO^S,1)-vc*qt))*dexp(-ly*x(ixO^S,2))*dcos(theta)
        w(ixO^S,mag(2))= Busr*dsin(kx*(x(ixO^S,1)-vc*qt))*dexp(-ly*x(ixO^S,2))
        w(ixO^S,mag(3))=-Busr*dcos(kx*(x(ixO^S,1)-vc*qt))*dexp(-ly*x(ixO^S,2))*dsin(theta)
        w(ixO^S,mag(:))=w(ixO^S,mag(:))-block%B0(ixO^S,:,0)
      else
        w(ixO^S,mag(1))=-Busr*dcos(kx*(x(ixO^S,1)-vc*qt))*dexp(-ly*x(ixO^S,2))*dcos(theta)
        w(ixO^S,mag(2))= Busr*dsin(kx*(x(ixO^S,1)-vc*qt))*dexp(-ly*x(ixO^S,2))
        w(ixO^S,mag(3))=-Busr*dcos(kx*(x(ixO^S,1)-vc*qt))*dexp(-ly*x(ixO^S,2))*dsin(theta)
      endif
      call mhd_to_conserved(ixI^L,ixO^L,w,x)
    case default
       call mpistop("Special boundary is not defined for this region")
    end select
    
  end subroutine specialbound_usr

  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
  ! this subroutine can be used in convert, to add auxiliary variables to the
  ! converted output file, for further analysis using tecplot, paraview, ....
  ! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
  ! the array normconv can be filled in the (nw+1:nw+nwauxio) range with
  ! corresponding normalization values (default value 1)
    use mod_global_parameters

    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)

    double precision :: pth(ixI^S),B2(ixI^S),divb(ixI^S)
    double precision :: Btotal(ixI^S,1:ndir)
    integer :: idir

    ! output temperature
    call mhd_get_pthermal(w,x,ixI^L,ixO^L,pth)
    w(ixO^S,nw+1)=pth(ixO^S)/w(ixO^S,rho_)

    do idir=1,ndir
      if(B0field) then
        Btotal(ixI^S,idir)=w(ixI^S,mag(idir))+block%B0(ixI^S,idir,0)
      else
        Btotal(ixI^S,idir)=w(ixI^S,mag(idir))
      endif
    end do
    ! B^2
    B2(ixO^S)=sum((Btotal(ixO^S,:))**2,dim=ndim+1)

    ! output Alfven wave speed B/sqrt(rho)
    w(ixO^S,nw+2)=dsqrt(B2(ixO^S)/w(ixO^S,rho_))

    ! output divB1
    call divvector(Btotal,ixI^L,ixO^L,divb)
    w(ixO^S,nw+3)=0.5d0*divb(ixO^S)/dsqrt(B2(ixO^S))/(^D&1.0d0/dxlevel(^D)+)
    ! output the plasma beta p*2/B**2
    w(ixO^S,nw+4)=pth(ixO^S)*two/B2(ixO^S)

  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
  ! newly added variables need to be concatenated with the w_names/primnames string
    use mod_global_parameters
    character(len=*) :: varnames

    varnames='Te Alfv divB beta'

  end subroutine specialvarnames_output

  subroutine specialset_B0(ixI^L,ixO^L,x,wB0)
  ! Here add a time-independent background magnetic field
    use mod_global_parameters

    integer, intent(in)           :: ixI^L,ixO^L
    double precision, intent(in)  :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: wB0(ixI^S,1:ndir)

    wB0(ixO^S,1)=wB0(ixO^S,1)-Busr*dcos(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))*dcos(theta)
    wB0(ixO^S,2)=wB0(ixO^S,2)+Busr*dsin(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))
    wB0(ixO^S,3)=wB0(ixO^S,3)-Busr*dcos(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))*dsin(theta)

  end subroutine specialset_B0

  subroutine specialset_J0(ixI^L,ixO^L,x,wJ0)
  ! Here add a time-independent background current density 
    use mod_global_parameters

    integer, intent(in)           :: ixI^L,ixO^L
    double precision, intent(in)  :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: wJ0(ixI^S,7-2*ndir:ndir)

    wJ0(ixO^S,1)= ly*Busr*dcos(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))*dsin(theta)
    wJ0(ixO^S,2)=-kx*Busr*dsin(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))*dsin(theta)
    wJ0(ixO^S,3)= kx*Busr*dcos(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))&
                 -ly*Busr*dcos(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))*dcos(theta)

  end subroutine specialset_J0

end module mod_usr
