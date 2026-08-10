!> module mod_pfss.t -- potential field source surface model
!> PURPOSE : to extrapolate global potential magnetic field of the sun from
!>           synoptic magnetograms
!> 2013.11.04 Developed by S. Moschou and C. Xia
!> 2014.04.01 Allow to change source surface (C. Xia)
!> PRECONDITIONS: 
!>  1. 3D spherical coordinates
!>  2. A synoptic magnetogram in a binary file contains nphi, ntheta, 
!>     theta(ntheta), phi(nphi), B_r(nphi,ntheta) succesively.
!>  3. nphi, ntheta are long integers and other arrays are double precision.
!>     theta contains  decreasing radians with increasing indice (Pi to 0) 
!>     phi contains increasing radians with increasing indice (0 to 2*Pi)
!>  4. By default, theta points are interpreted as cell centers and harmonic
!>     coefficients use cell-area quadrature in mu=cos(theta). To reproduce the
!>     old Gauss-Legendre quadrature, set
!>     pfss_theta_quadrature='gauss_legendre' before calling harm_coef.
!>     If a mapname.coef file already exists, it is reused and this setting does
!>     not take effect until that coefficient file is regenerated.
!> USAGE:
!>   example for a magnetogram with name 'mdicr2020.dat':
!>
!>   subroutine initglobaldata_usr
!>     ...
!>     R_0=1.d0 ! dimensionless Solar radius (default 1.0)
!>     R_s=2.5d0 ! dimensionless radius of source surface (default 2.5)
!>     lmax=120     ! use a fixed value instead of the value determined by the 
!>                    resolution of input magnetogram
!>     trunc=.true. ! use less spherical harmonics at larger distances
!>     call harm_coef('mdicr2020.dat')
!>   end subroutine initglobaldata_usr
!>
!>   subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
!>     ...
!>     double precision :: bpf(ixI^S,1:ndir)
!>     ...
!>     call pfss(ixI^L,ixO^L,bpf,x)
!>     w(ix^S,mag(:))=bpf(ix^S,:)
!>
!>   end subroutine initonegrid_usr
module mod_pfss
  use mod_fft, only: fft_raw
  implicit none
  private

  double complex, allocatable :: flm(:,:),Alm(:,:),Blm(:,:)
  double precision, allocatable :: Rlm(:,:), xrg(:)
  double precision, public :: R_s=2.5d0, R_0=1.d0
  integer, allocatable :: lmaxarray(:)
  integer, public :: lmax=0
  logical, public :: trunc=.false.
  character(len=20), public :: pfss_theta_quadrature='cell_area'
 
  public :: harm_coef
{^IFTHREED
  public :: pfss
}
  
contains

  subroutine harm_coef(mapname)
    use mod_global_parameters
    use mod_comm_lib, only: mpistop

    double precision, allocatable :: b_r0(:,:)
    double precision, allocatable :: theta(:),phi(:),cfwm(:)
    double precision :: rsl,xrl,dxr
    integer :: xm,ym,l,m,amode,file_handle,il,ir,nlarr,nsh
    integer, dimension(MPI_STATUS_SIZE) :: statuss
    logical :: aexist
    character(len=*) :: mapname
    character(len=80) :: fharmcoef

    fharmcoef=mapname//'.coef'
    inquire(file=fharmcoef, exist=aexist)
    if(aexist) then
      if(mype==0) write(*,'(2a)') &
        'Using existing PFSS coefficient file; pfss_theta_quadrature is ignored: ',&
        trim(fharmcoef)
      if(mype==0) then
        call MPI_FILE_OPEN(MPI_COMM_SELF,fharmcoef,MPI_MODE_RDONLY, &
                             MPI_INFO_NULL,file_handle,ierrmpi)
        call MPI_FILE_READ(file_handle,lmax,1,MPI_INTEGER,statuss,ierrmpi)
        allocate(flm(0:lmax,0:lmax))
        call MPI_FILE_READ(file_handle,flm,(lmax+1)*(lmax+1),&
                             MPI_DOUBLE_COMPLEX,statuss,ierrmpi)
        call MPI_FILE_CLOSE(file_handle,ierrmpi)
      end if
      call MPI_BARRIER(icomm,ierrmpi)
      if(npe>0)  call MPI_BCAST(lmax,1,MPI_INTEGER,0,icomm,ierrmpi)
      if(mype/=0) allocate(flm(0:lmax,0:lmax))
      call MPI_BARRIER(icomm,ierrmpi)
      if(npe>0) call MPI_BCAST(flm,(lmax+1)*(lmax+1),MPI_DOUBLE_COMPLEX,0,icomm,&
         ierrmpi)
    else
      if(mype==0) then
        inquire(file=mapname,exist=aexist)
        if(.not. aexist) then
          if(mype==0) write(*,'(2a)') "can not find file:",mapname
          call mpistop("no input magnetogram found")
        end if
        call MPI_FILE_OPEN(MPI_COMM_SELF,mapname,MPI_MODE_RDONLY,MPI_INFO_NULL,&
                           file_handle,ierrmpi)
        call MPI_FILE_READ(file_handle,xm,1,MPI_INTEGER,statuss,ierrmpi)
        call MPI_FILE_READ(file_handle,ym,1,MPI_INTEGER,statuss,ierrmpi)
        if(lmax==0) lmax=min(2*ym/3,xm/3)
    
        allocate(b_r0(xm,ym))
        allocate(theta(ym))
        allocate(phi(xm))
        call MPI_FILE_READ(file_handle,theta,ym,MPI_DOUBLE_PRECISION,&
                           statuss,ierrmpi)
        call MPI_FILE_READ(file_handle,phi,xm,MPI_DOUBLE_PRECISION,&
                           statuss,ierrmpi)
        call MPI_FILE_READ(file_handle,b_r0,xm*ym,MPI_DOUBLE_PRECISION,&
                           statuss,ierrmpi)
        call MPI_FILE_CLOSE(file_handle,ierrmpi)
        print*,'nphi,ntheta',xm,ym
        print*,'theta range:',minval(theta),maxval(theta)
        print*,'phi range:',minval(phi),maxval(phi)
        print*,'Brmax,Brmin',maxval(b_r0),minval(b_r0)
        allocate(cfwm(ym))
        select case(trim(pfss_theta_quadrature))
        case('cell_area')
          call cfweights_cell_area(ym,theta,cfwm)
        case('gauss_legendre')
          call cfweights_gauss_legendre(ym,dcos(theta),cfwm)
        case default
          call mpistop("Unknown pfss_theta_quadrature")
        end select
        allocate(flm(0:lmax,0:lmax))
        call coef(b_r0,xm,ym,dcos(theta),dsin(theta),cfwm)
        deallocate(b_r0)
        deallocate(theta)
        deallocate(phi)
        amode=ior(MPI_MODE_CREATE,MPI_MODE_WRONLY)
        call MPI_FILE_OPEN(MPI_COMM_SELF,fharmcoef,amode, &
                             MPI_INFO_NULL,file_handle,ierrmpi)
        call MPI_FILE_WRITE(file_handle,lmax,1,MPI_INTEGER,statuss,ierrmpi)
        call MPI_FILE_WRITE(file_handle,flm,(lmax+1)*(lmax+1),&
                             MPI_DOUBLE_COMPLEX,statuss,ierrmpi)
        call MPI_FILE_CLOSE(file_handle,ierrmpi)
      endif
      call MPI_BARRIER(icomm,ierrmpi)
      if(npe>1) call MPI_BCAST(lmax,1,MPI_INTEGER,0,icomm,ierrmpi)
      if(mype/=0) allocate(flm(0:lmax,0:lmax))
      call MPI_BARRIER(icomm,ierrmpi)
      if(npe>1) call MPI_BCAST(flm,(lmax+1)*(lmax+1),MPI_DOUBLE_COMPLEX,0,&
                     icomm,ierrmpi)
    end if
    if(mype==0) print*,'lmax=',lmax,'trunc=',trunc
    nlarr=501
    allocate(lmaxarray(nlarr))
    allocate(xrg(nlarr))
    lmaxarray=lmax
    if(trunc) then
      dxr=(R_s-R_0)/dble(nlarr-1)
      do ir=1,nlarr
        xrg(ir)=dxr*dble(ir-1)+R_0
        do il=0,lmax
          xrl=xrg(ir)**il
          if(xrl > 1.d6) then
            lmaxarray(ir)=il
            exit
          end if
        end do
      end do
    endif
    ! calculate global Alm Blm Rlm 
    allocate(Alm(0:lmax,0:lmax))
    allocate(Blm(0:lmax,0:lmax))
    allocate(Rlm(0:lmax,0:lmax))
    Alm=(0.d0,0.d0)
    Blm=(0.d0,0.d0)
    do l=0,lmax
      do m=0,l
        rsl=R_s**(-(2*l+1))
        Rlm(l,m)=dsqrt(dble(l**2-m**2)/dble(4*l**2-1))
        Blm(l,m)=-flm(l,m)/(1.d0+dble(l)+dble(l)*rsl)
        Alm(l,m)=-rsl*Blm(l,m)
      end do
    end do

  end subroutine harm_coef

  subroutine cfweights_cell_area(ym,theta,cfwm)
    use mod_global_parameters

    integer, intent(in) :: ym
    double precision, intent(in) :: theta(ym)
    double precision, intent(out) :: cfwm(ym)

    double precision,dimension(ym) :: miu
    double precision :: edge_l,edge_r
    integer :: i

    miu=dcos(theta)
    if(miu(1)<=miu(ym)) then
      do i=1,ym
        if(i==1) then
          edge_l=-1.d0
        else
          edge_l=0.5d0*(miu(i-1)+miu(i))
        end if
        if(i==ym) then
          edge_r=1.d0
        else
          edge_r=0.5d0*(miu(i)+miu(i+1))
        end if
        cfwm(i)=dabs(edge_r-edge_l)*(2.d0*dpi)
      end do
    else
      do i=1,ym
        if(i==1) then
          edge_l=1.d0
        else
          edge_l=0.5d0*(miu(i-1)+miu(i))
        end if
        if(i==ym) then
          edge_r=-1.d0
        else
          edge_r=0.5d0*(miu(i)+miu(i+1))
        end if
        cfwm(i)=dabs(edge_r-edge_l)*(2.d0*dpi)
      end do
    end if

  end subroutine cfweights_cell_area

  subroutine cfweights_gauss_legendre(ym,miu,cfwm)
    use mod_global_parameters

    integer, intent(in) :: ym
    double precision, intent(in) :: miu(ym)
    double precision, intent(out) :: cfwm(ym)

    double precision,dimension(ym) :: Pl,Pm2,Pm1,Pprime,sintheta
    double precision :: lr
    integer :: l

    sintheta=dsqrt(1.d0-miu**2)

    Pm2=1.d0
    Pm1=miu

    do l=2,ym-1
      lr=1.d0/dble(l)
      Pl=(2.d0-lr)*Pm1*miu-(1.d0-lr)*Pm2
      Pm2=Pm1
      Pm1=Pl
    end do

    Pprime=(dble(ym)*Pl)/sintheta**2
    cfwm=2.d0/(sintheta*Pprime)**2
    cfwm=cfwm*(2.d0*dpi)

  end subroutine cfweights_gauss_legendre

  subroutine coef(b_r0,xm,ym,miu,mius,cfwm)
    use mod_global_parameters

    integer, intent(in) :: xm,ym
    double precision, intent(in) :: b_r0(xm,ym),cfwm(ym),miu(ym),mius(ym)

    double complex :: Bm(0:xm-1,0:ym-1)
    double precision,dimension(xm) :: fftmr,fftmi
    double precision,dimension(0:lmax) :: N_mm
    double precision,dimension(ym) :: P_lm1,P_lm2,old_Pmm,P_l
    double precision :: mr,lr,c1,c2
    integer :: l,m,i,j,stat

    Bm=(0.d0,0.d0)
    do i=1,ym
      fftmr=b_r0(:,i)/dble(xm)
      fftmi=0.d0
      call fft_raw(fftmr,fftmi,xm,xm,xm,-1)
      Bm(:,i-1)=(fftmr+(0.d0,1.d0)*fftmi)
    end do
    N_mm(0)=1.d0/dsqrt(4.d0*dpi)
    do m=1,lmax
      N_mm(m)=-N_mm(m-1)*dsqrt(1.d0+1.d0/dble(2*m))
    end do
    !first do m=0
    P_lm2=N_mm(0)
    P_lm1=P_lm2*miu*dsqrt(3.d0)
    !set l=0 m=0 term
    flm(0,0)=sum(Bm(0,:)*P_lm2*cfwm)
    !set l=1 m=0 term
    flm(1,0)=sum(Bm(0,:)*P_lm1*cfwm)
    do l=2,lmax
      lr=dble(l)
      c1=dsqrt(4.d0-1.d0/lr**2) 
      c2=-(1.d0-1.d0/lr)*dsqrt((2.d0*lr+1.d0)/(2.d0*lr-3.d0))
      P_l=c1*miu*P_lm1+c2*P_lm2
      !set m=0 term for all other l's
      flm(l,0)=sum(Bm(0,:)*P_l*cfwm)
      P_lm2=P_lm1
      P_lm1=P_l
    end do

    !since only l modes from 0 to lmax are used
    Bm=2.d0*Bm

    !now the rest of the m's
    old_Pmm=N_mm(0)
    do m=1,lmax
      P_lm2=old_Pmm*mius*N_mm(m)/N_mm(m-1)
      P_lm1=P_lm2*miu*dsqrt(dble(2*m+3))
    !ACCURATE UP TO HERE
      old_Pmm=P_lm2
      !set l=m mode
      flm(m,m)=sum(Bm(m,:)*P_lm2*cfwm)
      !set l=m+1 mode
      if(m<lmax) flm(m+1,m)=sum(Bm(m,:)*P_lm1*cfwm)
      mr=dble(m)
      do l=m+2,lmax
        lr=dble(l)
        c1=dsqrt((4.d0*lr**2-1.d0)/(lr**2-mr**2))
        c2=-dsqrt(((2.d0*lr+1.d0)*((lr-1.d0)**2-mr**2))/((2.d0*lr-3.d0)*(lr**2-&
                  mr**2)))
        P_l=c1*miu*P_lm1+c2*P_lm2
        flm(l,m)=sum(Bm(m,:)*P_l*cfwm)
        P_lm2=P_lm1
        P_lm1=P_l
      end do
    end do

  end subroutine coef
{^IFTHREED
  subroutine pfss(ixI^L,ixO^L,Bpf,x)
    use mod_global_parameters

    integer, intent(in)           :: ixI^L,ixO^L
    double precision, intent(in)  :: x(ixI^S,1:ndim)
    double precision, intent(out) :: Bpf(ixI^S,1:ndir)

    double complex :: Bt(0:lmax,0:lmax,ixOmin1:ixOmax1)
    double precision :: phase(ixI^S,1:ndir),Bpfiv(ixOmin3:ixOmax3,ixOmin2:ixOmax2)
    double precision :: miu(ixOmin2:ixOmax2),mius(ixOmin2:ixOmax2),xr
    double precision :: tmp(ixOmin2:ixOmax2)
    integer :: l,m,ix^D,j,l1,l2,ntheta,nphi,ir,qlmax

    Bt=(0.d0,0.d0)
    nphi=ixOmax3-ixOmin3+1
    ntheta=ixOmax2-ixOmin2+1
    tmp(ixOmin2:ixOmax2)=x(ixOmin1,ixOmax2:ixOmin2:-1,ixOmin3,2)
    miu(ixOmin2:ixOmax2)=dcos(tmp(ixOmin2:ixOmax2))
    mius(ixOmin2:ixOmax2)=dsin(tmp(ixOmin2:ixOmax2))
    do ix1=ixOmin1,ixOmax1
      xr=x(ix1,ixOmin2,ixOmin3,1)
      if(trunc) then
        do ir=1,size(lmaxarray)
          if(xrg(ir)>=xr) exit
        end do
        if(ir>size(lmaxarray)) ir=size(lmaxarray)
        qlmax=lmaxarray(ir)
      else
        qlmax=lmax
      endif
    !Calculate Br
      do l=0,lmax
        do m=0,l
          Bt(l,m,ix1)=Alm(l,m)*dble(l)*xr**(l-1)-Blm(l,m)*dble(l+1)*xr**(-l-2)
        end do
      enddo
      call inv_sph_transform(Bt(:,:,ix1),x(ixOmin1,ixOmin2,&
           ixOmin3:ixOmax3,3),miu,mius,nphi,ntheta,Bpfiv,qlmax)
      do ix3=ixOmin3,ixOmax3
        do ix2=ixOmin2,ixOmax2
          Bpf(ix1,ix2,ix3,1)=Bpfiv(ix3,ixOmax2-ix2+ixOmin2)
        enddo
      enddo
    !Calculate Btheta
      do l=0,lmax
        do m=0,l
          if (l==0) then
            Bt(l,m,ix1)=-Rlm(l+1,m)*dble(l+2)*&
             (Alm(l+1,m)*xr**l+Blm(l+1,m)*xr**(-l-3))
          else if (l>=1 .and. l<=lmax-1) then
            Bt(l,m,ix1)=Rlm(l,m)*&
             dble(l-1)*(Alm(l-1,m)*xr**(l-2)+Blm(l-1,m)*&
             xr**(-l-1))-Rlm(l+1,m)*dble(l+2)*&
             (Alm(l+1,m)*xr**l+Blm(l+1,m)*xr**(-l-3))
          else
            Bt(l,m,ix1)=Rlm(l,m)*&
              dble(l-1)*(Alm(l-1,m)*xr**(l-2)+Blm(l-1,m)*xr**(-l-1))
          end if
        end do
      enddo
      call inv_sph_transform(Bt(:,:,ix1),x(ixOmin1,ixOmin2,&
           ixOmin3:ixOmax3,3),miu,mius,nphi,ntheta,Bpfiv,qlmax)
      do ix3=ixOmin3,ixOmax3
        do ix2=ixOmin2,ixOmax2
          Bpf(ix1,ix2,ix3,2)=Bpfiv(ix3,ixOmax2-ix2+ixOmin2)/mius(&
             ixOmax2-ix2+ixOmin2)
        enddo
      enddo
    
    !Calculate Bphi
      do l=0,lmax
        do m=0,l
          Bt(l,m,ix1)=(0.d0,1.d0)*m*(Alm(l,m)*xr**(l-1)+Blm(l,m)*xr**(-l-2))
        end do
      enddo
      call inv_sph_transform(Bt(:,:,ix1),x(ixOmin1,ixOmin2,&
           ixOmin3:ixOmax3,3),miu,mius,nphi,ntheta,Bpfiv,qlmax)
      do ix3=ixOmin3,ixOmax3
        do ix2=ixOmin2,ixOmax2
          Bpf(ix1,ix2,ix3,3)=Bpfiv(ix3,ixOmax2-ix2+ixOmin2)/mius(&
              ixOmax2-ix2+ixOmin2)
        enddo
      enddo
    enddo
  !Scalar Potential
  !       Potlc(ix^D)=Alm(l,m)*x(ix^D,1)**l+Blm(l,m)*x(ix^D,1)**(-l-1)
  
  !do ix3=ixOmin3,ixOmax3
  !    do ix2=ixOmin2,ixOmax2
  !    print*,x(ix1,ixOmax2-ix2+ixOmin2,ix3,2)
  !    print*,'miu==',miu(ixOmax2-ix2+ixOmin2)
  !    enddo
  !enddo
  
  end subroutine pfss

  subroutine inv_sph_transform(Bt,phi,miu,mius,nphi,ntheta,Bpf,qlmax)
    use mod_global_parameters

    integer, intent(in) :: nphi,ntheta,qlmax
    double complex, intent(in)  :: Bt(0:lmax,0:lmax)
    double precision, intent(in) :: phi(nphi),miu(ntheta),mius(ntheta)
    double precision, intent(out) :: Bpf(nphi,ntheta)

    double precision,dimension(0:lmax,0:lmax) :: cp,phase,Bamp
    double precision,dimension(ntheta) :: cp_1_0,cp_l_0,cp_lm1_0,cp_lm2_0,cp_m_m
    double precision,dimension(ntheta) :: cp_1_m,cp_l_m,cp_lm1_m,cp_lm2_m,cp_mp1_m
    double precision :: angpart(nphi)
    double precision :: ld,md,c1,c2,cp_0_0
    integer :: l,m,iph,ith 

    Bamp=abs(Bt)

    phase=atan2(dimag(Bt),dble(Bt))

    Bpf=0.d0
    !take care of modes where m=0
    cp_0_0=dsqrt(1.d0/(4.d0*dpi))
    !start with l=m=0 mode
    Bpf=Bpf+Bamp(0,0)*dcos(phase(0,0))*cp_0_0

    !proceed with l=1 m=0 mode
    cp_1_0=dsqrt(3.d0)*miu*cp_0_0
    do iph=1,nphi
     Bpf(iph,:)=Bpf(iph,:)+Bamp(1,0)*dcos(phase(1,0))*cp_1_0
    enddo

    !proceed with l modes for which m=0
    cp_lm1_0=cp_0_0
    cp_l_0=cp_1_0
    do l=2,qlmax
      ld=dble(l)
      cp_lm2_0=cp_lm1_0
      cp_lm1_0=cp_l_0
      c1=dsqrt(4.d0*ld**2-1.d0)/ld
      c2=dsqrt((2.d0*ld+1.d0)/(2.d0*ld-3.d0))*(ld-1.d0)/ld
      cp_l_0=c1*miu*cp_lm1_0-c2*cp_lm2_0
      do iph=1,nphi
       Bpf(iph,:)=Bpf(iph,:)+Bamp(l,0)*dcos(phase(l,0))*cp_l_0
      enddo
    enddo

    !loop through m's for m>0 and then loop through l's for each m
    cp_m_m=cp_0_0
    do m=1,qlmax
      md=dble(m)
      !first do l=m modes
      cp_m_m=-dsqrt(1.d0+1.d0/(2.d0*md))*mius*cp_m_m
      do iph=1,nphi 
        angpart(iph)=dcos(md*phi(iph)+phase(m,m))
      end do
      do ith=1,ntheta
        do iph=1,nphi
          Bpf(iph,ith)=Bpf(iph,ith)+Bamp(m,m)*angpart(iph)*cp_m_m(ith)
        enddo
      enddo

      !proceed with l=m+1 modes
      if(qlmax>=m+1) then
        cp_mp1_m=dsqrt(2.d0*md+3.d0)*miu*cp_m_m
        angpart=dcos(md*phi+phase(m+1,m))
        do ith=1,ntheta
          do iph=1,nphi
            Bpf(iph,ith)=Bpf(iph,ith)+Bamp(m+1,m)*angpart(iph)*cp_mp1_m(ith)
          enddo
        enddo
      endif

      !finish with the rest l  
      if(qlmax>=m+2) then
        cp_lm1_m=cp_m_m
        cp_l_m=cp_mp1_m
        do l=m+2,qlmax
          ld=dble(l)
          cp_lm2_m=cp_lm1_m
          cp_lm1_m=cp_l_m
          c1=dsqrt((4.d0*ld**2-1.d0)/(ld**2-md**2))
          c2=dsqrt((2.d0*ld+1.d0)*((ld-1.d0)**2-md**2)/(2.d0*ld-3.d0)/(ld**2-md**2))
          cp_l_m=c1*miu*cp_lm1_m-c2*cp_lm2_m
          angpart=dcos(md*phi+phase(l,m))
          do ith=1,ntheta
            do iph=1,nphi
              Bpf(iph,ith)=Bpf(iph,ith)+Bamp(l,m)*angpart(iph)*cp_l_m(ith)
            enddo
          enddo
        enddo
      endif
    enddo

  end subroutine inv_sph_transform
}

end module mod_pfss
