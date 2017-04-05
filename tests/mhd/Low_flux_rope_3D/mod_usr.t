module mod_usr
  use mod_mhd
  implicit none

contains

  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods

    usr_init_one_grid => initonegrid_usr
    usr_special_bc    => specialbound_usr
    usr_aux_output    => specialvar_output
    usr_add_aux_names => specialvarnames_output 
    usr_set_B0        => specialset_B0

    call set_coordinate_system("Cartesian")
    call mhd_activate()

  end subroutine usr_init

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
    call get_B(ixI^L,ixO^L,Bloc,x)
    w(ixO^S,mag(:))=Bloc(ixO^S,:)
    if(mhd_glm) w(ixO^S,psi_)=0.d0
    call mhd_to_conserved(ixI^L,ixO^L,w,x)

  end subroutine initonegrid_usr

  subroutine get_B(ixI^L,ixO^L,B,x)
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(out) :: B(ixI^S,1:ndir)

    double precision :: miu,park,parf(ixI^S),zz(ixI^S),zmin

    miu=0.95d0
    park=0.5d0
    zmin=3.d0
    zz(ixO^S)=x(ixO^S,3)-zmin
    parf(ixO^S)=4.d0*miu*miu/(1.d0+miu*miu)**2+park*park*x(ixO^S,1)**2+&
       (park*zz(ixO^S)+(1.d0-miu*miu)/(1.d0+miu*miu))**2
    B(ixO^S,1)=2.d0*Busr*(park*zz(ixO^S)+(1.d0-miu*miu)/(1.d0+miu*miu))/parf(ixO^S)
    B(ixO^S,2)=4.d0*miu*Busr/(1.d0+miu*miu)/parf(ixO^S)
    B(ixO^S,3)=-2.d0*Busr*park*x(ixO^S,1)/parf(ixO^S)

  end subroutine get_B

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
      call get_B(ixI^L,ixO^L,Bloc,x)
      w(ixO^S,mag(:))=Bloc(ixO^S,:)
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
      call get_B(ixI^L,ixO^L,Bloc,x)
      w(ixO^S,mag(:))=Bloc(ixO^S,:)
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
      call get_B(ixI^L,ixO^L,Bloc,x)
      w(ixO^S,mag(:))=Bloc(ixO^S,:)
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
      call get_B(ixI^L,ixO^L,Bloc,x)
      w(ixO^S,mag(:))=Bloc(ixO^S,:)
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
      call get_B(ixI^L,ixO^L,Bloc,x)
      w(ixO^S,mag(:))=Bloc(ixO^S,:)
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
      call get_B(ixI^L,ixO^L,Bloc,x)
      w(ixO^S,mag(:))=Bloc(ixO^S,:)
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
    w(ixO^S,nw+4)=current(ixO^S,2)

  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
  ! newly added variables need to be concatenated with the w_names/primnames string
    use mod_global_parameters
    character(len=*) :: varnames

    varnames='Te beta divb j2'

  end subroutine specialvarnames_output

  subroutine specialset_B0(ixI^L,ixO^L,x,wB0)
  ! Here one can add a steady (time-independent) background magnetic field
    use mod_global_parameters

    integer, intent(in)           :: ixI^L,ixO^L
    double precision, intent(in)  :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: wB0(ixI^S,1:ndir)

    double precision :: Bloc(ixI^S,1:ndir)

    call get_B(ixI^L,ixO^L,Bloc,x)
    wB0(ixO^S,:)=wB0(ixO^S,:)+Bloc(ixO^S,:)

  end subroutine specialset_B0

end module mod_usr
