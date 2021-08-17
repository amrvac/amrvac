module mod_usr
  use mod_mhd
  implicit none
  double precision :: miu,k,zc

contains

  subroutine usr_init()
    usr_set_parameters=> initglobaldata_usr
    usr_init_one_grid => initonegrid_usr
    usr_special_bc    => specialbound_usr
    usr_aux_output    => specialvar_output
    usr_add_aux_names => specialvarnames_output 
    usr_set_B0        => specialset_B0
    usr_set_J0        => specialset_J0
    usr_special_convert=> usrspecial_convert

    call set_coordinate_system("Cartesian")
    call mhd_activate()

  end subroutine usr_init

  subroutine initglobaldata_usr
    miu=1.d0
    k=0.5d0
    zc=0.5d0*(xprobmax3+xprobmin3)

  end subroutine initglobaldata_usr

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
  ! initialize one grid
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
    integer, intent(in)          :: ixI^L,ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision             :: w(ixI^S,nw+nwauxio)
    double precision             :: normconv(0:nw+nwauxio)

    double precision             :: tmp(ixI^S),current(ixI^S,7-2*ndir:3) 
    double precision             :: bvec(ixI^S,1:ndir),qvec(ixI^S,1:ndir)
    integer :: idirmin, idir, jdir, kdir

    call mhd_get_pthermal(w,x,ixI^L,ixO^L,tmp)
    ! output the temperature p/rho
    w(ixO^S,nw+1)=tmp(ixO^S)/w(ixO^S,rho_)
    ! output the plasma beta p*2/B**2
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
    ! output Lorenz force
    if(B0field)then
      bvec(ixO^S,:)=w(ixO^S,mag(:))+block%B0(ixO^S,:,0)
    else
      bvec(ixO^S,:)=w(ixO^S,mag(:))
    end if
    qvec=0.d0
    do idir=1,ndir; do jdir=idirmin,3; do kdir=1,ndir
      if(lvc(idir,jdir,kdir)/=0)then
        tmp(ixO^S)=current(ixO^S,jdir)*bvec(ixO^S,kdir)
        if(lvc(idir,jdir,kdir)==1)then
           qvec(ixO^S,idir)=qvec(ixO^S,idir)+tmp(ixO^S)
        else
           qvec(ixO^S,idir)=qvec(ixO^S,idir)-tmp(ixO^S)
        endif
      endif
    enddo; enddo; enddo
    w(ixO^S,nw+7)=qvec(ixO^S,1)
    w(ixO^S,nw+8)=qvec(ixO^S,2)
    w(ixO^S,nw+9)=qvec(ixO^S,3)
    ! output the Alfven speed
    if(B0field)then
      w(ixO^S,nw+10)=sqrt(sum((w(ixO^S,mag(:))+block%B0(ixO^S,:,0))**2,dim=ndim+1)/w(ixO^S,rho_))
    else
      w(ixO^S,nw+10)=sqrt(sum(w(ixO^S,mag(:))**2,dim=ndim+1)/w(ixO^S,rho_))
    endif

  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
  ! newly added variables need to be concatenated with the w_names/primnames string
    character(len=*) :: varnames
    varnames='Te beta divb j1 j2 j3 L1 L2 L3 Alfven'

  end subroutine specialvarnames_output

  subroutine specialset_B0(ixI^L,ixO^L,x,wB0)
  ! Here one can add a steady (time-independent) background magnetic field
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

  subroutine usrspecial_convert(qunitconvert)
    integer, intent(in) :: qunitconvert
    character(len=20):: userconvert_type
  
    call spatial_integral_w
  end subroutine usrspecial_convert

  subroutine spatial_integral_w
    double precision :: dvolume(ixG^T), dsurface(ixG^T),timephy,dvone
    double precision, allocatable :: integral_ipe(:), integral_w(:)

    integer           :: nregions,ireg,ncellpe,ncell,idims,hxM^LL,nx^D
    integer           :: iigrid,igrid,status(MPI_STATUS_SIZE),ni
    character(len=100):: filename,region
    character(len=1024) :: line, datastr
    logical           :: patchwi(ixG^T),alive

    nregions=1
    ! number of integrals to perform
    ni=3
    allocate(integral_ipe(ni),integral_w(ni))
    integral_ipe=0.d0
    integral_w=0.d0
    nx^D=ixMhi^D-ixMlo^D+1;
    do ireg=1,nregions
      select case(ireg)
      case(1)
        region='fulldomain'
      case(2)
        region='cropped'
      end select
      ncellpe=0 
      do iigrid=1,igridstail; igrid=igrids(iigrid);
        block=>ps(igrid)
        if(slab) then
          dvone={rnode(rpdx^D_,igrid)|*}
          dvolume(ixM^T)=dvone
          dsurface(ixM^T)=two*(^D&dvone/rnode(rpdx^D_,igrid)+)
        else
          dvolume(ixM^T)=ps(igrid)%dvolume(ixM^T)
          dsurface(ixM^T)= sum(ps(igrid)%surfaceC(ixM^T,:),dim=ndim+1)
          do idims=1,ndim
            hxM^LL=ixM^LL-kr(idims,^D);
            dsurface(ixM^T)=dsurface(ixM^T)+ps(igrid)%surfaceC(hxM^T,idims)
          end do
        end if
        ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
        patchwi(ixG^T)=.false.
        select case(region)
        case('cropped')
           call mask_grid(ixG^LL,ixM^LL,ps(igrid)%w,ps(igrid)%x,patchwi,ncellpe)
        case('fulldomain')
           patchwi(ixM^T)=.true.
           ncellpe=ncellpe+{nx^D*}
        case default
           call mpistop("region not defined")
        end select
        integral_ipe(1)=integral_ipe(1)+ &
                  integral_grid(ixG^LL,ixM^LL,ps(igrid)%w,ps(igrid)%x,dvolume,dsurface,4,patchwi)
        integral_ipe(2)=integral_ipe(2)+ &
                  integral_grid(ixG^LL,ixM^LL,ps(igrid)%w,ps(igrid)%x,dvolume,dsurface,5,patchwi)
        integral_ipe(3)=integral_ipe(3)+ &
                  integral_grid(ixG^LL,ixM^LL,ps(igrid)%w,ps(igrid)%x,dvolume,dsurface,6,patchwi)
      end do
      call MPI_ALLREDUCE(integral_ipe,integral_w,ni,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,icomm,ierrmpi)
      !call MPI_ALLREDUCE(ncellpe,ncell,1,MPI_INTEGER,MPI_SUM,icomm,ierrmpi)
      timephy=global_time
      if(mype==0) then
        write(filename,"(a,a,a)") TRIM(base_filename),TRIM(region),"aftime.csv"
        inquire(file=filename,exist=alive)
        if(alive) then
          open(unit=21,file=filename,form='formatted',status='old',access='append')
        else
          open(unit=21,file=filename,form='formatted',status='new')
          write(21,'(a)') 'time, emagnetic, einternal, current'
        endif
        write(datastr,'(es13.6, a)') timephy,','
        line=datastr
        write(datastr,"(es13.6, a)") integral_w(1),','
        line = trim(line)//trim(datastr)
        write(datastr,"(es13.6, a)") integral_w(2),','
        line = trim(line)//trim(datastr)
        write(datastr,"(es13.6)") integral_w(3)
        line = trim(line)//trim(datastr)
        write(21,'(a)') trim(line)
        close(21)
      endif
    enddo
    deallocate(integral_ipe,integral_w)
  end subroutine spatial_integral_w

  subroutine mask_grid(ixI^L,ixO^L,w,x,patchwi,cellcount)
    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    logical, intent(inout)             :: patchwi(ixG^T)

    double precision  ::  buff
    integer                            :: ix^D,cellcount

    buff=0.05d0*(xprobmax1-xprobmin1)
    {do ix^DB=ixOmin^DB,ixOmax^DB\}
       if(x(ix^D,1)>xprobmin1+buff .and. x(ix^D,1)<xprobmax1-buff .and. &
          x(ix^D,2)>xprobmin2+buff .and. x(ix^D,2)<xprobmax2-buff) then
         patchwi(ix^D)=.true.
         cellcount=cellcount+1
       else
         patchwi(ix^D)=.false.
       endif
    {end do\}
    return

  end subroutine mask_grid

  function integral_grid(ixI^L,ixO^L,w,x,dvolume,dsurface,intval,patchwi)
    integer, intent(in)                :: ixI^L,ixO^L,intval
    double precision, intent(in)       :: x(ixI^S,1:ndim),dvolume(ixG^T),dsurface(ixG^T)
    double precision, intent(in)       :: w(ixI^S,nw)
    logical, intent(in) :: patchwi(ixG^T)
    
    double precision, dimension(ixG^T,1:ndir) :: bvec,qvec
    double precision :: current(ixG^T,7-2*ndir:3),tmp(ixG^T)
    double precision :: integral_grid,mcurrent
    integer :: ix^D,idirmin,idir,jdir,kdir

    integral_grid=0.d0
    select case(intval)
     case(1)
      ! magnetic energy
      if(B0field)then
        tmp(ixO^S)=0.5d0*sum((w(ixO^S,mag(:))+&
                      block%B0(ixO^S,:,0))**2,dim=ndim+1)
      else
        tmp(ixO^S)=0.5d0*sum(w(ixO^S,mag(:))**2,dim=ndim+1)
      endif
      {do ix^DB=ixOmin^DB,ixOmax^DB\}
         if(patchwi(ix^D)) integral_grid=integral_grid+tmp(ix^D)*dvolume(ix^D)
      {end do\}
     case(2)
      ! internal energy
      call mhd_get_pthermal(w,x,ixI^L,ixO^L,tmp)
      {do ix^DB=ixOmin^DB,ixOmax^DB\}
         if(patchwi(ix^D))  integral_grid=integral_grid+tmp(ix^D)/(mhd_gamma-1.d0)*dvolume(ix^D)
      {end do\}
     case(3)
      ! current strength
      call get_current(w,ixI^L,ixO^L,idirmin,current)
      {do ix^DB=ixOmin^DB,ixOmax^DB\}
         if(patchwi(ix^D)) integral_grid=integral_grid+dsqrt(sum(current(ix^D,:)**2))*dvolume(ix^D)
      {end do\}
     case(4)
      ! Alfven crossing time x
      if(B0field)then
        tmp(ixO^S)=sum((w(ixO^S,mag(:))+&
                      block%B0(ixO^S,:,0))**2,dim=ndim+1)
      else
        tmp(ixO^S)=sum(w(ixO^S,mag(:))**2,dim=ndim+1)
      endif
      do ix3=ixOmin3,ixOmin3
      do ix2=ixOmin2,ixOmin2
      do ix1=ixOmin1,ixOmax1
         if(patchwi(ix^D)) integral_grid=integral_grid+dxlevel(1)/dsqrt(tmp(ix^D)/w(ix^D,rho_))/64.d0
      end do
      end do
      end do
     case(5)
      ! Alfven crossing time y
      if(B0field)then
        tmp(ixO^S)=sum((w(ixO^S,mag(:))+&
                      block%B0(ixO^S,:,0))**2,dim=ndim+1)
      else
        tmp(ixO^S)=sum(w(ixO^S,mag(:))**2,dim=ndim+1)
      endif
      do ix3=ixOmin3,ixOmin3
      do ix2=ixOmin2,ixOmax2
      do ix1=ixOmin1,ixOmin1
         if(patchwi(ix^D)) integral_grid=integral_grid+dxlevel(2)/dsqrt(tmp(ix^D)/w(ix^D,rho_))/64.d0
      end do
      end do
      end do
     case(6)
      ! Alfven crossing time z
      if(B0field)then
        tmp(ixO^S)=sum((w(ixO^S,mag(:))+&
                      block%B0(ixO^S,:,0))**2,dim=ndim+1)
      else
        tmp(ixO^S)=sum(w(ixO^S,mag(:))**2,dim=ndim+1)
      endif
      do ix3=ixOmin3,ixOmax3
      do ix2=ixOmin2,ixOmin2
      do ix1=ixOmin1,ixOmin1
         if(patchwi(ix^D)) integral_grid=integral_grid+dxlevel(3)/dsqrt(tmp(ix^D)/w(ix^D,rho_))/64.d0
      end do
      end do
      end do
     case default
         call mpistop("intval not defined")
    end select
    
    return
  end function integral_grid

end module mod_usr
