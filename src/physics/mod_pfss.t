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
  implicit none
  private

  double complex, allocatable :: flm(:,:),Alm(:,:),Blm(:,:)
  double precision, allocatable :: Rlm(:,:), xrg(:)
  double precision, public :: R_s=2.5d0, R_0=1.d0
  integer, allocatable :: lmaxarray(:)
  integer, public :: lmax=0
  logical, public :: trunc=.false.
 
  public :: harm_coef
{^IFTHREED
  public :: pfss
}
  
contains

  subroutine harm_coef(mapname)
    use mod_global_parameters

    double precision, allocatable :: b_r0(:,:)
    double precision, allocatable :: theta(:),phi(:),cfwm(:)
    double precision :: rsl,xrl,dxr
    integer :: xm,ym,l,m,amode,file_handle,il,ir,nlarr,nsh
    integer, dimension(MPI_STATUS_SIZE) :: statuss
    character(len=*) :: mapname
    character(len=80) :: fharmcoef
    logical :: aexist

    fharmcoef=mapname//'.coef'
    inquire(file=fharmcoef, exist=aexist)
    if(aexist) then
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
        call cfweights(ym,dcos(theta),cfwm)
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

  subroutine cfweights(ym,miu,cfwm)
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

  end subroutine cfweights

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
      call fft(fftmr,fftmi,xm,xm,xm,-1)
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
    integer :: l,m,ix^D,j,l1,l2,ntheta,nphi,ir,qlmax

    Bt=(0.d0,0.d0)
    nphi=ixOmax3-ixOmin3+1
    ntheta=ixOmax2-ixOmin2+1
    miu(ixOmin2:ixOmax2)=dcos(x(ixOmin1,ixOmax2:ixOmin2:-1,ixOmin3,2))
    mius(ixOmin2:ixOmax2)=dsin(x(ixOmin1,ixOmax2:ixOmin2:-1,ixOmin3,2))
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
  subroutine fft(a,b,ntot,n,nspan,isn)
  !  multivariate complex fourier transform, computed in place
  !    using mixed-radix fast fourier transform algorithm.
  !  by r. c. singleton, stanford research institute, sept. 1968
  !  arrays a and b originally hold the real and imaginary
  !    components of the data, and return the real and
  !    imaginary components of the resulting fourier coefficients.
  !  multivariate data is indexed according to the fortran
  !    array element successor function, without limit
  !    on the number of implied multiple subscripts.
  !    the subroutine is called once for each variate.
  !    the calls for a multivariate transform may be in any order.
  !  ntot is the total number of complex data values.
  !  n is the dimension of the current variable.
  !  nspan/n is the spacing of consecutive data values
  !    while indexing the current variable.
  !  the sign of isn determines the sign of the complex
  !    exponential, and the magnitude of isn is normally one.
  !  a tri-variate transform with a(n1,n2,n3), b(n1,n2,n3)
  !    is computed by
  !      call fft(a,b,n1*n2*n3,n1,n1,1)
  !      call fft(a,b,n1*n2*n3,n2,n1*n2,1)
  !      call fft(a,b,n1*n2*n3,n3,n1*n2*n3,1)
  !  for a single-variate transform,
  !    ntot = n = nspan = (number of complex data values), e.g.
  !      call fft(a,b,n,n,n,1)
  !  the data can alternatively be stored in a single complex array c
  !    in standard fortran fashion, i.e. alternating real and imaginary
  !    parts. then with most fortran compilers, the complex array c can
  !    be equivalenced to a real array a, the magnitude of isn changed
  !    to two to give correct indexing increment, and a(1) and a(2) used
  !    to pass the initial addresses for the sequences of real and
  !    imaginary values, e.g.
  !       complex c(ntot)
  !       real    a(2*ntot)
  !       equivalence (c(1),a(1))
  !       call fft(a(1),a(2),ntot,n,nspan,2)
  !  arrays at(maxf), ck(maxf), bt(maxf), sk(maxf), and np(maxp)
  !    are used for temporary storage.  if the available storage
  !    is insufficient, the program is terminated by a stop.
  !    maxf must be .ge. the maximum prime factor of n.
  !    maxp must be .gt. the number of prime factors of n.
  !    in addition, if the square-free portion k of n has two or
  !    more prime factors, then maxp must be .ge. k-1.
        double precision :: a(:),b(:)
  !  array storage in nfac for a maximum of 15 prime factors of n.
  !  if n has more than one square-free factor, the product of the
  !    square-free factors must be .le. 210
        dimension nfac(11),np(209)
  !  array storage for maximum prime factor of 23
        dimension at(23),ck(23),bt(23),sk(23)
        integer :: i,ii,maxp,maxf,n,inc,isn,nt,ntot,ks,nspan,kspan,nn,jc,jf,m
        integer :: k,j,jj,nfac,kt,np,kk,k1,k2,k3,k4,kspnn
        double precision :: c72,s72,s120,rad,radf,sd,cd,ak,bk,c1
        double precision :: s1,aj,bj,akp,ajp,ajm,akm,bkp,bkm,bjp,bjm,aa
        double precision :: bb,sk,ck,at,bt,s3,c3,s2,c2
        equivalence (i,ii)
  
  !  the following two constants should agree with the array dimensions.
        maxp=209
  ! Date: Wed, 9 Aug 1995 09:38:49 -0400
  ! From: ldm@apollo.numis.nwu.edu
        maxf=23
        s3=0.d0
        s2=0.d0
        c3=0.d0
        c2=0.d0
        if(n .lt. 2) return
        inc=isn
        c72=0.30901699437494742d0
        s72=0.95105651629515357d0
        s120=0.86602540378443865d0
        rad=6.2831853071796d0
        if(isn .ge. 0) go to 10
        s72=-s72
        s120=-s120
        rad=-rad
        inc=-inc
     10 nt=inc*ntot
        ks=inc*nspan
        kspan=ks
        nn=nt-inc
        jc=ks/n
        radf=rad*dble(jc)*0.5d0
        i=0
        jf=0
  !  determine the factors of n
        m=0
        k=n
        go to 20
     15 m=m+1
        nfac(m)=4
        k=k/16
     20 if(k-(k/16)*16 .eq. 0) go to 15
        j=3
        jj=9
        go to 30
     25 m=m+1
        nfac(m)=j
        k=k/jj
     30 if(mod(k,jj) .eq. 0) go to 25
        j=j+2
        jj=j**2
        if(jj .le. k) go to 30
        if(k .gt. 4) go to 40
        kt=m
        nfac(m+1)=k
        if(k .ne. 1) m=m+1
        go to 80
     40 if(k-(k/4)*4 .ne. 0) go to 50
        m=m+1
        nfac(m)=2
        k=k/4
     50 kt=m
        j=2
     60 if(mod(k,j) .ne. 0) go to 70
        m=m+1
        nfac(m)=j
        k=k/j
     70 j=((j+1)/2)*2+1
        if(j .le. k) go to 60
     80 if(kt .eq. 0) go to 100
        j=kt
     90 m=m+1
        nfac(m)=nfac(j)
        j=j-1
        if(j .ne. 0) go to 90
  !  compute fourier transform
    100 sd=radf/dble(kspan)
        cd=2.d0*dsin(sd)**2
        sd=dsin(sd+sd)
        kk=1
        i=i+1
        if(nfac(i) .ne. 2) go to 400
  !  transform for factor of 2 (including rotation factor)
        kspan=kspan/2
        k1=kspan+2
    210 k2=kk+kspan
        ak=a(k2)
        bk=b(k2)
        a(k2)=a(kk)-ak
        b(k2)=b(kk)-bk
        a(kk)=a(kk)+ak
        b(kk)=b(kk)+bk
        kk=k2+kspan
        if(kk .le. nn) go to 210
        kk=kk-nn
        if(kk .le. jc) go to 210
        if(kk .gt. kspan) go to 800
    220 c1=1.d0-cd
        s1=sd
    230 k2=kk+kspan
        ak=a(kk)-a(k2)
        bk=b(kk)-b(k2)
        a(kk)=a(kk)+a(k2)
        b(kk)=b(kk)+b(k2)
        a(k2)=c1*ak-s1*bk
        b(k2)=s1*ak+c1*bk
        kk=k2+kspan
        if(kk .lt. nt) go to 230
        k2=kk-nt
        c1=-c1
        kk=k1-k2
        if(kk .gt. k2) go to 230
        ak=c1-(cd*c1+sd*s1)
        s1=(sd*c1-cd*s1)+s1
        c1=2.d0-(ak**2+s1**2)
        s1=c1*s1
        c1=c1*ak
        kk=kk+jc
        if(kk .lt. k2) go to 230
        k1=k1+inc+inc
        kk=(k1-kspan)/2+jc
        if(kk .le. jc+jc) go to 220
        go to 100
  !  transform for factor of 3 (optional code)
    320 k1=kk+kspan
        k2=k1+kspan
        ak=a(kk)
        bk=b(kk)
        aj=a(k1)+a(k2)
        bj=b(k1)+b(k2)
        a(kk)=ak+aj
        b(kk)=bk+bj
        ak=-0.5d0*aj+ak
        bk=-0.5d0*bj+bk
        aj=(a(k1)-a(k2))*s120
        bj=(b(k1)-b(k2))*s120
        a(k1)=ak-bj
        b(k1)=bk+aj
        a(k2)=ak+bj
        b(k2)=bk-aj
        kk=k2+kspan
        if(kk .lt. nn) go to 320
        kk=kk-nn
        if(kk .le. kspan) go to 320
        go to 700
  !  transform for factor of 4
    400 if(nfac(i) .ne. 4) go to 600
        kspnn=kspan
        kspan=kspan/4
    410 c1=1.d0
        s1=0.d0
    420 k1=kk+kspan
        k2=k1+kspan
        k3=k2+kspan
        akp=a(kk)+a(k2)
        akm=a(kk)-a(k2)
        ajp=a(k1)+a(k3)
        ajm=a(k1)-a(k3)
        a(kk)=akp+ajp
        ajp=akp-ajp
        bkp=b(kk)+b(k2)
        bkm=b(kk)-b(k2)
        bjp=b(k1)+b(k3)
        bjm=b(k1)-b(k3)
        b(kk)=bkp+bjp
        bjp=bkp-bjp
        if(isn .lt. 0) go to 450
        akp=akm-bjm
        akm=akm+bjm
        bkp=bkm+ajm
        bkm=bkm-ajm
        if(s1 .eq. 0.d0) go to 460
    430 a(k1)=akp*c1-bkp*s1
        b(k1)=akp*s1+bkp*c1
        a(k2)=ajp*c2-bjp*s2
        b(k2)=ajp*s2+bjp*c2
        a(k3)=akm*c3-bkm*s3
        b(k3)=akm*s3+bkm*c3
        kk=k3+kspan
        if(kk .le. nt) go to 420
    440 c2=c1-(cd*c1+sd*s1)
        s1=(sd*c1-cd*s1)+s1
        c1=2.d0-(c2**2+s1**2)
        s1=c1*s1
        c1=c1*c2
        c2=c1**2-s1**2
        s2=2.d0*c1*s1
        c3=c2*c1-s2*s1
        s3=c2*s1+s2*c1
        kk=kk-nt+jc
        if(kk .le. kspan) go to 420
        kk=kk-kspan+inc
        if(kk .le. jc) go to 410
        if(kspan .eq. jc) go to 800
        go to 100
    450 akp=akm+bjm
        akm=akm-bjm
        bkp=bkm-ajm
        bkm=bkm+ajm
        if(s1 .ne. 0) go to 430
    460 a(k1)=akp
        b(k1)=bkp
        a(k2)=ajp
        b(k2)=bjp
        a(k3)=akm
        b(k3)=bkm
        kk=k3+kspan
        if(kk .le. nt) go to 420
        go to 440
  !  transform for factor of 5 (optional code)
    510 c2=c72**2-s72**2
        s2=2.d0*c72*s72
    520 k1=kk+kspan
        k2=k1+kspan
        k3=k2+kspan
        k4=k3+kspan
        akp=a(k1)+a(k4)
        akm=a(k1)-a(k4)
        bkp=b(k1)+b(k4)
        bkm=b(k1)-b(k4)
        ajp=a(k2)+a(k3)
        ajm=a(k2)-a(k3)
        bjp=b(k2)+b(k3)
        bjm=b(k2)-b(k3)
        aa=a(kk)
        bb=b(kk)
        a(kk)=aa+akp+ajp
        b(kk)=bb+bkp+bjp
        ak=akp*c72+ajp*c2+aa
        bk=bkp*c72+bjp*c2+bb
        aj=akm*s72+ajm*s2
        bj=bkm*s72+bjm*s2
        a(k1)=ak-bj
        a(k4)=ak+bj
        b(k1)=bk+aj
        b(k4)=bk-aj
        ak=akp*c2+ajp*c72+aa
        bk=bkp*c2+bjp*c72+bb
        aj=akm*s2-ajm*s72
        bj=bkm*s2-bjm*s72
        a(k2)=ak-bj
        a(k3)=ak+bj
        b(k2)=bk+aj
        b(k3)=bk-aj
        kk=k4+kspan
        if(kk .lt. nn) go to 520
        kk=kk-nn
        if(kk .le. kspan) go to 520
        go to 700
  !  transform for odd factors
    600 k=nfac(i)
        kspnn=kspan
        kspan=kspan/k
        if(k .eq. 3) go to 320
        if(k .eq. 5) go to 510
        if(k .eq. jf) go to 640
        jf=k
        s1=rad/dble(k)
        c1=dcos(s1)
        s1=dsin(s1)
        if(jf .gt. maxf) go to 998
        ck(jf)=1.d0
        sk(jf)=0.d0
        j=1
    630 ck(j)=ck(k)*c1+sk(k)*s1
        sk(j)=ck(k)*s1-sk(k)*c1
        k=k-1
        ck(k)=ck(j)
        sk(k)=-sk(j)
        j=j+1
        if(j .lt. k) go to 630
    640 k1=kk
        k2=kk+kspnn
        aa=a(kk)
        bb=b(kk)
        ak=aa
        bk=bb
        j=1
        k1=k1+kspan
    650 k2=k2-kspan
        j=j+1
        at(j)=a(k1)+a(k2)
        ak=at(j)+ak
        bt(j)=b(k1)+b(k2)
        bk=bt(j)+bk
        j=j+1
        at(j)=a(k1)-a(k2)
        bt(j)=b(k1)-b(k2)
        k1=k1+kspan
        if(k1 .lt. k2) go to 650
        a(kk)=ak
        b(kk)=bk
        k1=kk
        k2=kk+kspnn
        j=1
    660 k1=k1+kspan
        k2=k2-kspan
        jj=j
        ak=aa
        bk=bb
        aj=0.d0
        bj=0.d0
        k=1
    670 k=k+1
        ak=at(k)*ck(jj)+ak
        bk=bt(k)*ck(jj)+bk
        k=k+1
        aj=at(k)*sk(jj)+aj
        bj=bt(k)*sk(jj)+bj
        jj=jj+j
        if(jj .gt. jf) jj=jj-jf
        if(k .lt. jf) go to 670
        k=jf-j
        a(k1)=ak-bj
        b(k1)=bk+aj
        a(k2)=ak+bj
        b(k2)=bk-aj
        j=j+1
        if(j .lt. k) go to 660
        kk=kk+kspnn
        if(kk .le. nn) go to 640
        kk=kk-nn
        if(kk .le. kspan) go to 640
  !  multiply by rotation factor (except for factors of 2 and 4)
    700 if(i .eq. m) go to 800
        kk=jc+1
    710 c2=1.d0-cd
        s1=sd
    720 c1=c2
        s2=s1
        kk=kk+kspan
    730 ak=a(kk)
        a(kk)=c2*ak-s2*b(kk)
        b(kk)=s2*ak+c2*b(kk)
        kk=kk+kspnn
        if(kk .le. nt) go to 730
        ak=s1*s2
        s2=s1*c2+c1*s2
        c2=c1*c2-ak
        kk=kk-nt+kspan
        if(kk .le. kspnn) go to 730
        c2=c1-(cd*c1+sd*s1)
        s1=s1+(sd*c1-cd*s1)
        c1=2.d0-(c2**2+s1**2)
        s1=c1*s1
        c2=c1*c2
        kk=kk-kspnn+jc
        if(kk .le. kspan) go to 720
        kk=kk-kspan+jc+inc
        if(kk .le. jc+jc) go to 710
        go to 100
  !  permute the results to normal order---done in two stages
  !  permutation for square factors of n
    800 np(1)=ks
        if(kt .eq. 0) go to 890
        k=kt+kt+1
        if(m .lt. k) k=k-1
        j=1
        np(k+1)=jc
    810 np(j+1)=np(j)/nfac(j)
        np(k)=np(k+1)*nfac(j)
        j=j+1
        k=k-1
        if(j .lt. k) go to 810
        k3=np(k+1)
        kspan=np(2)
        kk=jc+1
        k2=kspan+1
        j=1
        if(n .ne. ntot) go to 850
  !  permutation for single-variate transform (optional code)
    820 ak=a(kk)
        a(kk)=a(k2)
        a(k2)=ak
        bk=b(kk)
        b(kk)=b(k2)
        b(k2)=bk
        kk=kk+inc
        k2=kspan+k2
        if(k2 .lt. ks) go to 820
    830 k2=k2-np(j)
        j=j+1
        k2=np(j+1)+k2
        if(k2 .gt. np(j)) go to 830
        j=1
    840 if(kk .lt. k2) go to 820
        kk=kk+inc
        k2=kspan+k2
        if(k2 .lt. ks) go to 840
        if(kk .lt. ks) go to 830
        jc=k3
        go to 890
  !  permutation for multivariate transform
    850 k=kk+jc
    860 ak=a(kk)
        a(kk)=a(k2)
        a(k2)=ak
        bk=b(kk)
        b(kk)=b(k2)
        b(k2)=bk
        kk=kk+inc
        k2=k2+inc
        if(kk .lt. k) go to 860
        kk=kk+ks-jc
        k2=k2+ks-jc
        if(kk .lt. nt) go to 850
        k2=k2-nt+kspan
        kk=kk-nt+jc
        if(k2 .lt. ks) go to 850
    870 k2=k2-np(j)
        j=j+1
        k2=np(j+1)+k2
        if(k2 .gt. np(j)) go to 870
        j=1
    880 if(kk .lt. k2) go to 850
        kk=kk+jc
        k2=kspan+k2
        if(k2 .lt. ks) go to 880
        if(kk .lt. ks) go to 870
        jc=k3
    890 if(2*kt+1 .ge. m) return
        kspnn=np(kt+1)
  !  permutation for square-free factors of n
        j=m-kt
        nfac(j+1)=1
    900 nfac(j)=nfac(j)*nfac(j+1)
        j=j-1
        if(j .ne. kt) go to 900
        kt=kt+1
        nn=nfac(kt)-1
        if(nn .gt. maxp) go to 998
        jj=0
        j=0
        go to 906
    902 jj=jj-k2
        k2=kk
        k=k+1
        kk=nfac(k)
    904 jj=kk+jj
        if(jj .ge. k2) go to 902
        np(j)=jj
    906 k2=nfac(kt)
        k=kt+1
        kk=nfac(k)
        j=j+1
        if(j .le. nn) go to 904
  !  determine the permutation cycles of length greater than 1
        j=0
        go to 914
    910 k=kk
        kk=np(k)
        np(k)=-kk
        if(kk .ne. j) go to 910
        k3=kk
    914 j=j+1
        kk=np(j)
        if(kk .lt. 0) go to 914
        if(kk .ne. j) go to 910
        np(j)=-j
        if(j .ne. nn) go to 914
        maxf=inc*maxf
  !  reorder a and b, following the permutation cycles
        go to 950
    924 j=j-1
        if(np(j) .lt. 0) go to 924
        jj=jc
    926 kspan=jj
        if(jj .gt. maxf) kspan=maxf
        jj=jj-kspan
        k=np(j)
        kk=jc*k+ii+jj
        k1=kk+kspan
        k2=0
    928 k2=k2+1
        at(k2)=a(k1)
        bt(k2)=b(k1)
        k1=k1-inc
        if(k1 .ne. kk) go to 928
    932 k1=kk+kspan
        k2=k1-jc*(k+np(k))
        k=-np(k)
    936 a(k1)=a(k2)
        b(k1)=b(k2)
        k1=k1-inc
        k2=k2-inc
        if(k1 .ne. kk) go to 936
        kk=k2
        if(k .ne. j) go to 932
        k1=kk+kspan
        k2=0
    940 k2=k2+1
        a(k1)=at(k2)
        b(k1)=bt(k2)
        k1=k1-inc
        if(k1 .ne. kk) go to 940
        if(jj .ne. 0) go to 926
        if(j .ne. 1) go to 924
    950 j=k3+1
        nt=nt-kspnn
        ii=nt-inc+1
        if(nt .ge. 0) go to 924
        return
  !  error finish, insufficient array storage
    998 isn=0
        print 999
        stop
    999 format(44h0array bounds exceeded within subroutine fft)
  end subroutine fft

end module mod_pfss
