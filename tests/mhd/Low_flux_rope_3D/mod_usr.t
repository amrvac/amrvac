module mod_usr
  use mod_mhd
  implicit none
  double precision :: miu,k,zc

contains

  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods

    usr_set_parameters=> initglobaldata_usr
    usr_init_one_grid => initonegrid_usr
    usr_special_bc    => specialbound_usr
    usr_aux_output    => specialvar_output
    usr_add_aux_names => specialvarnames_output 
    usr_set_B0        => specialset_B0
    usr_set_J0        => specialset_J0

    call set_coordinate_system("Cartesian")
    call mhd_activate()

  end subroutine usr_init

  subroutine initglobaldata_usr
    use mod_global_parameters
    miu=1.d0
    k=0.5d0
    zc=0.5d0*(xprobmax3+xprobmin3)

  end subroutine initglobaldata_usr

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
  ! initialize one grid
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: Bloc(ixI^S,1:ndir)
    logical, save:: first=.true.

    if (first) then
       if (mype==0) then
          print *,'3D MHD BC Low force-free flux rope in Cartesian coordinate'
       end if
       first=.false.
    end if
    w(ixO^S,rho_)=1.d0
    w(ixO^S,p_)=1.d0
    w(ixO^S,mom(:))=zero
    if(B0field) then
      w(ixO^S,mag(:))=0.d0
    else
      call specialset_B0(ixI^L,ixO^L,x,Bloc)
      w(ixO^S,mag(:))=Bloc(ixO^S,:)
    end if
    if(mhd_glm) w(ixO^S,psi_)=0.d0
    call mhd_to_conserved(ixI^L,ixO^L,w,x)

  end subroutine initonegrid_usr

  subroutine specialbound_usr(qt,ixI^L,ixO^L,iB,w,x)
    ! special boundary types, user defined
    use mod_global_parameters

    integer, intent(in) :: ixO^L, iB, ixI^L
    double precision, intent(in) :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: pth(ixI^S),Bloc(ixI^S,1:ndir)
    integer :: ix^D, ixA^L

    select case(iB)
    case(1)
      if(B0field) then
        w(ixO^S,mag(:))=0.d0
      else
        call specialset_B0(ixI^L,ixO^L,x,Bloc)
        w(ixO^S,mag(:))=Bloc(ixO^S,:)
      end if
      ixA^L=ixO^L;
      ixAmin1=ixOmax1+1;ixAmax1=ixOmax1+1;
      call mhd_get_pthermal(w,x,ixI^L,ixA^L,pth)
      do ix1=ixOmin1,ixOmax1
        w(ix1^%1ixO^S,rho_)=w(ixOmax1+1^%1ixO^S,rho_)
        w(ix1^%1ixO^S,mom(1))=w(ixOmax1+1^%1ixO^S,mom(1))/w(ixOmax1+1^%1ixO^S,rho_)
        w(ix1^%1ixO^S,mom(2))=w(ixOmax1+1^%1ixO^S,mom(2))/w(ixOmax1+1^%1ixO^S,rho_)
        w(ix1^%1ixO^S,mom(3))=w(ixOmax1+1^%1ixO^S,mom(3))/w(ixOmax1+1^%1ixO^S,rho_)
        w(ix1^%1ixO^S,p_)=pth(ixOmax1+1^%1ixO^S)
      enddo
      call mhd_to_conserved(ixI^L,ixO^L,w,x)
    case(2)
      if(B0field) then
        w(ixO^S,mag(:))=0.d0
      else
        call specialset_B0(ixI^L,ixO^L,x,Bloc)
        w(ixO^S,mag(:))=Bloc(ixO^S,:)
      end if
      ixA^L=ixO^L;
      ixAmin1=ixOmin1-1;ixAmax1=ixOmin1-1;
      call mhd_get_pthermal(w,x,ixI^L,ixA^L,pth)
      do ix1=ixOmin1,ixOmax1
        w(ix1^%1ixO^S,rho_)=w(ixOmin1-1^%1ixO^S,rho_)
        w(ix1^%1ixO^S,mom(1))=w(ixOmin1-1^%1ixO^S,mom(1))/w(ixOmin1-1^%1ixO^S,rho_)
        w(ix1^%1ixO^S,mom(2))=w(ixOmin1-1^%1ixO^S,mom(2))/w(ixOmin1-1^%1ixO^S,rho_)
        w(ix1^%1ixO^S,mom(3))=w(ixOmin1-1^%1ixO^S,mom(3))/w(ixOmin1-1^%1ixO^S,rho_)
        w(ix1^%1ixO^S,p_)=pth(ixOmin1-1^%1ixO^S)
      enddo
      call mhd_to_conserved(ixI^L,ixO^L,w,x)
    case(3)
      if(B0field) then
        w(ixO^S,mag(:))=0.d0
      else
        call specialset_B0(ixI^L,ixO^L,x,Bloc)
        w(ixO^S,mag(:))=Bloc(ixO^S,:)
      end if
      ixA^L=ixO^L;
      ixAmin2=ixOmax2+1;ixAmax2=ixOmax2+1;
      call mhd_get_pthermal(w,x,ixI^L,ixA^L,pth)
      do ix2=ixOmin2,ixOmax2
        w(ix2^%2ixO^S,rho_)=w(ixOmax2+1^%2ixO^S,rho_)
        w(ix2^%2ixO^S,mom(1))=w(ixOmax2+1^%2ixO^S,mom(1))/w(ixOmax2+1^%2ixO^S,rho_)
        w(ix2^%2ixO^S,mom(2))=w(ixOmax2+1^%2ixO^S,mom(2))/w(ixOmax2+1^%2ixO^S,rho_)
        w(ix2^%2ixO^S,mom(3))=w(ixOmax2+1^%2ixO^S,mom(3))/w(ixOmax2+1^%2ixO^S,rho_)
        w(ix2^%2ixO^S,p_)=pth(ixOmax2+1^%2ixO^S)
      enddo
      call mhd_to_conserved(ixI^L,ixO^L,w,x)
    case(4)
      if(B0field) then
        w(ixO^S,mag(:))=0.d0
      else
        call specialset_B0(ixI^L,ixO^L,x,Bloc)
        w(ixO^S,mag(:))=Bloc(ixO^S,:)
      end if
      ixA^L=ixO^L;
      ixAmin2=ixOmin2-1;ixAmax2=ixOmin2-1;
      call mhd_get_pthermal(w,x,ixI^L,ixA^L,pth)
      do ix2=ixOmin2,ixOmax2
        w(ix2^%2ixO^S,rho_)=w(ixOmin2-1^%2ixO^S,rho_)
        w(ix2^%2ixO^S,mom(1))=w(ixOmin2-1^%2ixO^S,mom(1))/w(ixOmin2-1^%2ixO^S,rho_)
        w(ix2^%2ixO^S,mom(2))=w(ixOmin2-1^%2ixO^S,mom(2))/w(ixOmin2-1^%2ixO^S,rho_)
        w(ix2^%2ixO^S,mom(3))=w(ixOmin2-1^%2ixO^S,mom(3))/w(ixOmin2-1^%2ixO^S,rho_)
        w(ix2^%2ixO^S,p_)=pth(ixOmin2-1^%2ixO^S)
      enddo
      call mhd_to_conserved(ixI^L,ixO^L,w,x)
    case(5)
      if(B0field) then
        w(ixO^S,mag(:))=0.d0
      else
        call specialset_B0(ixI^L,ixO^L,x,Bloc)
        w(ixO^S,mag(:))=Bloc(ixO^S,:)
      end if
      ixA^L=ixO^L;
      ixAmin3=ixOmax3+1;ixAmax3=ixOmax3+1;
      call mhd_get_pthermal(w,x,ixI^L,ixA^L,pth)
      do ix3=ixOmin3,ixOmax3
        w(ix3^%3ixO^S,rho_)=w(ixOmax3+1^%3ixO^S,rho_)
        w(ix3^%3ixO^S,mom(1))=w(ixOmax3+1^%3ixO^S,mom(1))/w(ixOmax3+1^%3ixO^S,rho_)
        w(ix3^%3ixO^S,mom(2))=w(ixOmax3+1^%3ixO^S,mom(2))/w(ixOmax3+1^%3ixO^S,rho_)
        w(ix3^%3ixO^S,mom(3))=w(ixOmax3+1^%3ixO^S,mom(3))/w(ixOmax3+1^%3ixO^S,rho_)
        w(ix3^%3ixO^S,p_)=pth(ixOmax3+1^%3ixO^S)
      enddo
      call mhd_to_conserved(ixI^L,ixO^L,w,x)
    case(6)
      if(B0field) then
        w(ixO^S,mag(:))=0.d0
      else
        call specialset_B0(ixI^L,ixO^L,x,Bloc)
        w(ixO^S,mag(:))=Bloc(ixO^S,:)
      end if
      ixA^L=ixO^L;
      ixAmin3=ixOmin3-1;ixAmax3=ixOmin3-1;
      call mhd_get_pthermal(w,x,ixI^L,ixA^L,pth)
      do ix3=ixOmin3,ixOmax3
        w(ix3^%3ixO^S,rho_)=w(ixOmin3-1^%3ixO^S,rho_)
        w(ix3^%3ixO^S,mom(1))=w(ixOmin3-1^%3ixO^S,mom(1))/w(ixOmin3-1^%3ixO^S,rho_)
        w(ix3^%3ixO^S,mom(2))=w(ixOmin3-1^%3ixO^S,mom(2))/w(ixOmin3-1^%3ixO^S,rho_)
        w(ix3^%3ixO^S,mom(3))=w(ixOmin3-1^%3ixO^S,mom(3))/w(ixOmin3-1^%3ixO^S,rho_)
        w(ix3^%3ixO^S,p_)=pth(ixOmin3-1^%3ixO^S)
      enddo
      call mhd_to_conserved(ixI^L,ixO^L,w,x)
    case default
      call mpistop("Special boundary is not defined for this region")
    end select
  end subroutine specialbound_usr

  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
    use mod_global_parameters

    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)

    double precision                   :: tmp(ixI^S),current(ixI^S,7-2*ndir:3) 
    integer :: idirmin

    call mhd_get_pthermal(w,x,ixI^L,ixO^L,tmp)
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
    ! output current
    call get_current(w,ixI^L,ixO^L,idirmin,current)
    w(ixO^S,nw+4)=current(ixO^S,1)
    w(ixO^S,nw+5)=current(ixO^S,2)
    w(ixO^S,nw+6)=current(ixO^S,3)

  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
  ! newly added variables need to be concatenated with the w_names/primnames string
    use mod_global_parameters
    character(len=*) :: varnames

    varnames='Te beta divb j1 j2 j3'

  end subroutine specialvarnames_output

  subroutine specialset_B0(ixI^L,ixO^L,x,wB0)
  ! Here one can add a steady (time-independent) background magnetic field
    use mod_global_parameters

    integer, intent(in)           :: ixI^L,ixO^L
    double precision, intent(in)  :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: wB0(ixI^S,1:ndir)

    double precision :: f(ixI^S),zz(ixI^S)

    zz(ixO^S)=k*(x(ixO^S,3)-zc)+(1.d0-miu**2)/(1.d0+miu**2)
    f(ixO^S)=Busr/(4.d0*miu**2/(1.d0+miu**2)**2+k**2*x(ixO^S,2)**2+zz(ixO^S)**2)
    wB0(ixO^S,1)=4.d0*miu/(1.d0+miu**2)*f(ixO^S)
    wB0(ixO^S,2)=2.d0*zz(ixO^S)*f(ixO^S)
    wB0(ixO^S,3)=-2.d0*k*x(ixO^S,2)*f(ixO^S)

  end subroutine specialset_B0

  subroutine specialset_J0(ixI^L,ixO^L,x,wJ0)
  ! Here one can add a steady (time-independent) background current
    use mod_global_parameters

    integer, intent(in)           :: ixI^L,ixO^L
    double precision, intent(in)  :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: wJ0(ixI^S,1:ndir)

    double precision :: f(ixI^S),zz(ixI^S)

    zz(ixO^S)=k*(x(ixO^S,3)-zc)+(1.d0-miu**2)/(1.d0+miu**2)
    f(ixO^S)=Busr/(4.d0*miu**2/(1.d0+miu**2)**2+k**2*x(ixO^S,2)**2+zz(ixO^S)**2)
    wJ0(ixO^S,1)=-4.d0*k*f(ixO^S)+4.d0*k/Busr*f(ixO^S)**2*(k**2*x(ixO^S,2)**2+zz(ixO^S)**2)
    wJ0(ixO^S,2)=8.d0*miu*k/(1.d0+miu**2)/Busr*zz(ixO^S)*f(ixO^S)**2
    wJ0(ixO^S,3)=-8.d0*miu*k/(1.d0+miu**2)/Busr*x(ixO^S,2)*f(ixO^S)**2

  end subroutine specialset_J0

end module mod_usr
