! blast wave
module mod_usr
  use mod_mhd
  implicit none

contains

  subroutine usr_init()
    usr_init_one_grid => initonegrid_usr
    usr_aux_output    => specialvar_output
    usr_add_aux_names => specialvarnames_output 
    usr_set_B0        => specialset_B0

    call set_coordinate_system("spherical")
    call mhd_activate()

  end subroutine usr_init

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
  ! initialize one grid
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: rbs,xc^D,xcc^D
    double precision :: xcart(ixI^S,1:ndim),Bloc(ixI^S,ndir)
    integer :: ixC^L,idir
    logical, save:: first=.true.

    if (first) then
       if (mype==0) then
          print *,'3D MHD blast wave in spherical coordinate'
       end if
       first=.false.
    end if
    rbs=0.2d0
    w(ixO^S,rho_)=1.d0
    w(ixO^S,p_)=1.d0
    xc1=(xprobmin1+xprobmax1)*0.5d0
    !xc2=(xprobmin2+xprobmax2)*0.5d0
    xc2=xprobmin2+(xprobmax2-xprobmin2)*0.1d0
    xc3=(xprobmin3+xprobmax3)*0.5d0
    xcc1=xc1*dsin(xc2)*dcos(xc3)
    xcc2=xc1*dsin(xc2)*dsin(xc3)
    xcc3=xc1*dcos(xc2)
    xcart(ixO^S,1)=x(ixO^S,1)*dsin(x(ixO^S,2))*dcos(x(ixO^S,3))
    xcart(ixO^S,2)=x(ixO^S,1)*dsin(x(ixO^S,2))*dsin(x(ixO^S,3))
    xcart(ixO^S,3)=x(ixO^S,1)*dcos(x(ixO^S,2))
    where({^D&(xcart(ixO^S,^D)-xcc^D)**2+}<rbs**2)
      w(ixO^S,p_)=100.d0
    endwhere
    w(ixO^S,mom(:))=0.d0
    if(B0field) then
      w(ixO^S,mag(:))=0.d0
    else if(stagger_grid) then
      do idir=1,ndim
        xcart=x
        ixCmax^D=ixOmax^D;
        ixCmin^D=ixOmin^D-kr(idir,^D);
        xcart(ixC^S,idir)=x(ixC^S,idir)+0.5d0*block%dx(ixC^S,idir)
        call get_B(ixI^L,ixC^L,Bloc,xcart)
        block%ws(ixC^S,idir)=Bloc(ixC^S,idir)
      end do
      call mhd_face_to_center(ixO^L,block)
    else
      call get_B(ixI^L,ixO^L,Bloc,x)
      w(ixO^S,mag(:))=Bloc(ixO^S,:)
    end if

    if(mhd_glm) w(ixO^S,psi_)=0.d0

    call mhd_to_conserved(ixI^L,ixO^L,w,x)

  end subroutine initonegrid_usr

  subroutine get_B(ixI^L,ixO^L,B,x)
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(out) :: B(ixI^S,1:ndir)

    B(ixO^S,1)=Busr*(dsin(x(ixO^S,2))*dsin(x(ixO^S,3))+&
                          dcos(x(ixO^S,2)))
    B(ixO^S,2)=Busr*(dcos(x(ixO^S,2))*dsin(x(ixO^S,3))-&
                          dsin(x(ixO^S,2)))
    B(ixO^S,3)=Busr*(dcos(x(ixO^S,3)))

  end subroutine get_B

  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
  ! this subroutine can be used in convert, to add auxiliary variables to the
  ! converted output file, for further analysis using tecplot, paraview, ....
  ! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
  !
  ! the array normconv can be filled in the (nw+1:nw+nwauxio) range with
  ! corresponding normalization values (default value 1)
    use mod_physics
    use mod_mhd_phys

    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)

    double precision                   :: tmp(ixI^S) 

    call phys_get_pthermal(w,x,ixI^L,ixO^L,tmp)
    ! output the temperature p/rho
    w(ixO^S,nw+1)=tmp(ixO^S)/w(ixO^S,rho_)
    !! output the plasma beta p*2/B**2
    if(B0field)then
      w(ixO^S,nw+2)=tmp(ixO^S)*two/sum((w(ixO^S,mag(:))+&
                    block%B0(ixO^S,:,0))**2,dim=ndim+1)
    else
      w(ixO^S,nw+2)=tmp(ixO^S)*two/sum(w(ixO^S,mag(:))**2,dim=ndim+1)
    endif
    ! output divB1
    call get_divb(w,ixI^L,ixO^L,tmp)
    w(ixO^S,nw+3)=tmp(ixO^S)
    
  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
  ! newly added variables need to be concatenated with the w_names/primnames string
    character(len=*) :: varnames
    varnames='Te beta divb'

  end subroutine specialvarnames_output

  subroutine specialset_B0(ixI^L,ixO^L,x,wB0)
  ! Here one can add a steady (time-independent) potential background field
    integer, intent(in)           :: ixI^L,ixO^L
    double precision, intent(in)  :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: wB0(ixI^S,1:ndir)

    double precision :: Bloc(ixI^S,1:ndir)

    call get_B(ixI^L,ixO^L,Bloc,x)
    wB0(ixO^S,:)=wB0(ixO^S,:)+Bloc(ixO^S,:)

  end subroutine specialset_B0

end module mod_usr
