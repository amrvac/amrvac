! blast wave
module mod_usr
  use mod_mhd
  implicit none

contains

  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods

    usr_init_one_grid => initonegrid_usr
    usr_aux_output    => specialvar_output
    usr_add_aux_names => specialvarnames_output 
    usr_set_B0        => specialset_B0

    call set_coordinate_system("Cartesian")
    call mhd_activate()

  end subroutine usr_init

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
  ! initialize one grid
    use mod_global_parameters
    use mod_physics

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: rbs,xc1,xc2,Bloc(ixI^S,1:ndir)
    logical, save:: first=.true.

    if (first) then
       if (mype==0) then
          print *,'2.5D MHD blast wave in Cartesian coordinate'
       end if
       first=.false.
    end if
    w(ixO^S,rho_)=1.d0
    w(ixO^S,p_)=1.d0
    rbs=0.2d0
    xc1=(xprobmin1+xprobmax1)*0.5d0
    xc2=(xprobmin2+xprobmax2)*0.5d0
    where((x(ixO^S,1)-xc1)**2+(x(ixO^S,2)-xc2)**2<rbs**2)
      w(ixO^S,p_)=100.d0
    endwhere
    if(B0field) then
      w(ixO^S,mag(:))=0.d0
    else
      call get_B(ixI^L,ixO^L,Bloc,x)
      w(ixO^S,mag(:))=Bloc(ixO^S,:)
    end if
    w(ixO^S,mom(:))=0.d0

    if(mhd_glm) w(ixO^S,psi_)=0.d0

    call phys_to_conserved(ixI^L,ixO^L,w,x)

  end subroutine initonegrid_usr

  subroutine get_B(ixI^L,ixO^L,B,x)
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(out) :: B(ixI^S,1:ndir)

    B(ixO^S,:)=Busr

  end subroutine get_B

  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
  ! this subroutine can be used in convert, to add auxiliary variables to the
  ! converted output file, for further analysis using tecplot, paraview, ....
  ! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
  !
  ! the array normconv can be filled in the (nw+1:nw+nwauxio) range with
  ! corresponding normalization values (default value 1)
    use mod_global_parameters
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
                    myB0_cell%w(ixO^S,:))**2,dim=ndim+1)
    else
      w(ixO^S,nw+2)=tmp(ixO^S)*two/sum(w(ixO^S,mag(:))**2,dim=ndim+1)
    endif
    ! output divB1
    call get_divb(w,ixI^L,ixO^L,tmp)
    w(ixO^S,nw+3)=tmp(ixO^S)
    
  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
  ! newly added variables need to be concatenated with the w_names/primnames string
    use mod_global_parameters
    character(len=*) :: varnames

    varnames='Te beta divb'

  end subroutine specialvarnames_output

  subroutine specialset_B0(ixI^L,ixO^L,x,wB0)
  ! Here one can add a steady (time-independent) potential background field
    use mod_global_parameters

    integer, intent(in)           :: ixI^L,ixO^L
    double precision, intent(in)  :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: wB0(ixI^S,1:ndir)

    double precision :: Bloc(ixI^S,1:ndir)

    call get_B(ixI^L,ixO^L,Bloc,x)
    wB0(ixO^S,:)=wB0(ixO^S,:)+Bloc(ixO^S,:)

  end subroutine specialset_B0

end module mod_usr
