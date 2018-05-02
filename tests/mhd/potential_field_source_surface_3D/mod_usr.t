module mod_usr
  use mod_mhd
  use mod_pfss
  implicit none
  double precision :: usr_grav,rhob,Tiso

contains

  subroutine usr_init()
    unit_length=6.955d10 !cm
    unit_numberdensity=1.d9 !cm^-3
    unit_temperature=1.d6 !K

    usr_set_parameters => initglobaldata_usr
    usr_init_one_grid  => initonegrid_usr
    usr_gravity        => gravity
    usr_aux_output     => specialvar_output
    usr_add_aux_names  => specialvarnames_output

    call set_coordinate_system('spherical_3D')
    call mhd_activate()

  end subroutine usr_init

  subroutine initglobaldata_usr 
    usr_grav=-2.74d4*unit_length/unit_velocity**2

    ! base density and temperature
    rhob=1.d0
    Tiso=1.d6/unit_temperature

    trunc=.true.
    lmax=720
    call harm_coef('hmisyncr2119.dat')

  end subroutine initglobaldata_usr 

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    !initialize one grid within ixI^L
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    
    double precision :: rss,bpf(ixI^S,1:ndir),bpfss(ixI^S,1:ndir),xrss(ixI^S,1:ndim)
    double precision :: fmass,BB(ixI^S)
    integer :: ix1
    logical, save:: first=.true.

    rss=R_s
    if(mype==0 .and. first) then
      print*,'Global Sun PFSS model'
      first=.false.
    endif
    
    w(ixO^S,rho_)=rhob*dexp(usr_grav*R_0**2/Tiso*(1.d0/R_0-1.d0/(x(ixO^S,1))))
    
    if(B0field) then
      w(ixO^S,mag(1))=0.d0
      w(ixO^S,mag(2))=0.d0
      w(ixO^S,mag(3))=0.d0
    else 
      xrss=x
      xrss(ixI^S,r_)=rss
      call pfss(ixI^L,ixO^L,bpfss,xrss) 
      if(any(x(ixO^S,r_)<rss)) then
        call pfss(ixI^L,ixO^L,bpf,x) 
        w(ixO^S,mag(1):mag(3))=bpf(ixO^S,1:3)
        where(x(ixO^S,r_)>=rss)
          w(ixO^S,mag(1))=bpfss(ixO^S,1)*(rss/x(ixO^S,r_))**2
          w(ixO^S,mag(2))=0.d0
          w(ixO^S,mag(3))=0.d0
        endwhere
      else
        w(ixO^S,mag(1))=bpfss(ixO^S,1)*(rss/x(ixO^S,r_))**2
        w(ixO^S,mag(2))=0.d0
        w(ixO^S,mag(3))=0.d0
      endif
    endif
    w(ixO^S,mag(1):mag(3))=w(ixO^S,mag(1):mag(3))/unit_magneticfield
    
    w(ixO^S,mom(:))=0.d0

  end subroutine initonegrid_usr

  subroutine gravity(ixI^L,ixO^L,wCT,x,gravity_field)
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(in)    :: wCT(ixI^S,1:nw)
    double precision, intent(out)   :: gravity_field(ixI^S,ndim)
    double precision                :: ggrid(ixI^S)

    gravity_field=0.d0
    call getggrav(ggrid,ixI^L,ixO^L,x)
    gravity_field(ixO^S,1)=ggrid(ixO^S)

  end subroutine gravity

  subroutine getggrav(ggrid,ixI^L,ixO^L,x)
    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(out) :: ggrid(ixI^S)

    integer :: ix1

    do ix1=ixOmin1,ixOmax1
      ggrid(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=usr_grav*(R_0/&
          x(ix1,ixOmin2,ixOmin3,1))**2
    enddo

  end subroutine getggrav

  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)
    double precision                   :: csound2(ixI^S)

    w(ixO^S,nw+1)=w(ixO^S,mag(1))
    w(ixO^S,nw+2)=w(ixO^S,mag(2))
    w(ixO^S,nw+3)=w(ixO^S,mag(3))

  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
    character(len=*) varnames

    varnames='br bth bph'

  end subroutine specialvarnames_output

end module mod_usr
