
!! Copyright (C) 2002-2007 BigDFT group 
!! This file is distributed under the terms of the
!! GNU General Public License, see ~/COPYING file
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the list of contributors, see ~/AUTHORS 


!!****h* BigDFT/createKernel
!! NAME
!!    createKernel
!!
!! FUNCTION
!!    Allocate a pointer which corresponds to the zero-padded FFT slice needed for
!!    calculating the convolution with the kernel expressed in the interpolating scaling
!!    function basis. The kernel pointer is unallocated on input, allocated on output.
!!
!! SYNOPSIS
!!    geocode  Indicates the boundary conditions (BC) of the problem:
!!            'F' free BC, isolated systems.
!!                The program calculates the solution as if the given density is
!!                "alone" in R^3 space.
!!            'S' surface BC, isolated in y direction, periodic in xz plane                
!!                The given density is supposed to be periodic in the xz plane,
!!                so the dimensions in these direction mus be compatible with the FFT
!!                Beware of the fact that the isolated direction is y!
!!            'P' periodic BC.
!!                The density is supposed to be periodic in all the three directions,
!!                then all the dimensions must be compatible with the FFT.
!!                No need for setting up the kernel.
!!    iproc,nproc number of process, number of processes
!!    n01,n02,n03 dimensions of the real space grid to be hit with the Poisson Solver
!!    itype_scf   order of the interpolating scaling functions used in the decomposition
!!    hx,hy,hz    grid spacings. For the isolated BC case for the moment they are supposed to 
!!                be equal in the three directions
!!    kernel      pointer for the kernel FFT. Unallocated on input, allocated on output.
!!                Its dimensions are equivalent to the region of the FFT space for which the
!!                kernel is injective. This will divide by two each direction, 
!!                since the kernel for the zero-padded convolution is real and symmetric.
!!
!! WARNING
!!    Due to the fact that the kernel dimensions are unknown before the calling, the kernel
!!    must be declared as pointer in input of this routine.
!!    To avoid that, one can properly define the kernel dimensions by adding 
!!    the nd1,nd2,nd3 arguments to the PS_dim4allocation routine, then eliminating the pointer
!!    declaration.
!!
!! AUTHOR
!!    Luigi Genovese
!! CREATION DATE
!!    February 2007
!!
!! SOURCE
!!
subroutine createKernel(geocode,n01,n02,n03,hx,hy,hz,itype_scf,iproc,nproc,kernel)
  use mpi
  implicit none
  character(len=1), intent(in) :: geocode
  integer, intent(in) :: n01,n02,n03,itype_scf,iproc,nproc
  real(kind=8), intent(in) :: hx,hy,hz
  real(kind=8), pointer :: kernel(:)
  !local variables
  integer :: m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,i_stat
  integer :: jproc,nlimd,nlimk,jfd,jhd,jzd,jfk,jhk,jzk,npd,npk
  real(kind=8) :: hgrid

  call timing(iproc,'PSolvKernel   ','ON')

  hgrid=max(hx,hy,hz)

  if (iproc==iproc_verbose) write(*,'(1x,a)')&
          '------------------------------------------------------------ Poisson Kernel Creation'


  if (geocode == 'P') then
     
     if (iproc==iproc_verbose) write(*,'(1x,a)',advance='no')&
          'Poisson solver for periodic BC, no kernel calculation...'
     
     call F_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc)

     allocate(kernel(1),stat=i_stat)
     call memocc(i_stat,product(shape(kernel))*kind(kernel),'kernel','createkernel')

     nlimd=n2
     nlimk=0

  else if (geocode == 'S') then
     
     if (iproc==iproc_verbose) write(*,'(1x,a)',advance='no')&
          'Calculating Poisson solver kernel, surfaces BC...'

     !Build the Kernel
     call S_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc)

     allocate(kernel(nd1*nd2*nd3/nproc),stat=i_stat)
     call memocc(i_stat,product(shape(kernel))*kind(kernel),'kernel','createkernel')

     !the kernel must be built and scattered to all the processes
     call Surfaces_Kernel(n1,n2,n3,m3,nd1,nd2,nd3,hx,hz,hy,itype_scf,kernel,iproc,nproc)

     !last plane calculated for the density and the kernel
     nlimd=n2
     nlimk=n3/2+1

  else if (geocode == 'F') then

     if (iproc==iproc_verbose) write(*,'(1x,a)',advance='no')&
          'Calculating Poisson solver kernel, free BC...'

     !Build the Kernel
     call F_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc)
     allocate(kernel(nd1*nd2*nd3/nproc),stat=i_stat)
     call memocc(i_stat,product(shape(kernel))*kind(kernel),'kernel','createkernel')

     !the kernel must be built and scattered to all the processes

     call Free_Kernel(n01,n02,n03,n1,n2,n3,nd1,nd2,nd3,hx,hy,hz,itype_scf,iproc,nproc,kernel)


     !last plane calculated for the density and the kernel
     nlimd=n2/2
     nlimk=n3/2+1

  else
     
     if (iproc==iproc_verbose) &
          write(*,'(1x,a,3a)')'createKernel, geocode not admitted',geocode

     stop
  end if

  if (iproc==iproc_verbose) then
     write(*,'(a)')'done.'
     write(*,'(1x,2(a,i0))')&
          'Memory occ. per proc. (Bytes):  Density=',md1*md3*md2/nproc*8,&
          '  Kernel=',nd1*nd2*nd3/nproc*8
     write(*,'(1x,a,i0)')&
          '                                Full Grid Arrays=',n01*n02*n03*8
     !print the load balancing of the different dimensions on screen
     if (nproc > 1) then
        write(*,'(1x,a)')&
             'Load Balancing for Poisson Solver related operations:'
        jhd=1000
        jzd=1000
        npd=0
        load_balancing: do jproc=0,nproc-1
           !print *,'jproc,jfull=',jproc,jproc*md2/nproc,(jproc+1)*md2/nproc
           if ((jproc+1)*md2/nproc <= nlimd) then
              jfd=jproc
           else if (jproc*md2/nproc <= nlimd) then
              jhd=jproc
              npd=nint(real(nlimd-(jproc)*md2/nproc,kind=8)/real(md2/nproc,kind=8)*100.d0)
           else
              jzd=jproc
              exit load_balancing
           end if
        end do load_balancing
        write(*,'(1x,a,i3,a)')&
             'LB_density        : processors   0  -',jfd,' work at 100%'
        if (jfd < nproc-1) write(*,'(1x,a,i3,a,i3,1a)')&
             '                    processor     ',jhd,&
             '   works at ',npd,'%'
        if (jhd < nproc-1) write(*,'(1x,a,i3,1a,i3,a)')&
             '                    processors ',&
             jzd,'  -',nproc-1,' work at   0%'
        jhk=1000
        jzk=1000
        npk=0
        if (geocode /= 'P') then
           load_balancingk: do jproc=0,nproc-1
              !print *,'jproc,jfull=',jproc,jproc*nd3/nproc,(jproc+1)*nd3/nproc
              if ((jproc+1)*nd3/nproc <= nlimk) then
                 jfk=jproc
              else if (jproc*nd3/nproc <= nlimk) then
                 jhk=jproc
                 npk=nint(real(nlimk-(jproc)*nd3/nproc,kind=8)/real(nd3/nproc,kind=8)*100.d0)
              else
                 jzk=jproc
                 exit load_balancingk
              end if
           end do load_balancingk
           write(*,'(1x,a,i3,a)')&
                ' LB_kernel        : processors   0  -',jfk,' work at 100%'
           if (jfk < nproc-1) write(*,'(1x,a,i3,a,i3,1a)')&
                '                    processor     ',jhk,&
                '   works at ',npk,'%'
           if (jhk < nproc-1) write(*,'(1x,a,i3,1a,i3,a)')&
                '                    processors ',jzk,'  -',nproc-1,&
                ' work at   0%'
        end if
        write(*,'(1x,a)')&
             'Complete LB per proc.= 1/3 LB_density + 2/3 LB_kernel'
     end if

  end if
  call timing(iproc,'PSolvKernel   ','OF')

end subroutine createKernel


!!****h* BigDFT/Surfaces_Kernel
!! NAME
!!   Surfaces_Kernel
!!
!! FUNCTION
!!    Build the kernel of the Poisson operator with
!!    surfaces Boundary conditions
!!    in an interpolating scaling functions basis.
!!    Beware of the fact that the nonperiodic direction is y!
!!
!! SYNOPSIS
!!    n1,n2,n3           Dimensions for the FFT
!!    m3                 Actual dimension in non-periodic direction
!!    nker1,nker2,nker3  Dimensions of the kernel (nker3=n3/2+1) nker(1,2)=n(1,2)/2+1
!!    h1,h2,h3           Mesh steps in the three dimensions
!!    itype_scf          Order of the scaling function
!!    iproc,nproc        Number of process, number of processes
!!    karray             output array
!!
!! AUTHOR
!!    L. Genovese
!! CREATION DATE
!!    October 2006
!!
!! SOURCE
!!
subroutine Surfaces_Kernel(n1,n2,n3,m3,nker1,nker2,nker3,h1,h2,h3,itype_scf,karray,iproc,nproc)
  
  use mpi
  implicit none
  include 'perfdata.inc'
  
  !Arguments
  integer, intent(in) :: n1,n2,n3,m3,nker1,nker2,nker3,itype_scf,iproc,nproc
  real(kind=8), intent(in) :: h1,h2,h3
  real(kind=8), dimension(nker1,nker2,nker3/nproc), intent(out) :: karray
  
  !Local variables 
  !Better if higher (1024 points are enough 10^{-14}: 2*itype_scf*n_points)
  integer, parameter :: n_points = 2**6
  !Maximum number of points for FFT (should be same number in fft3d routine)
  integer, parameter :: nfft_max=24000
  
  real(kind=8), dimension(:), allocatable :: kernel_scf
  real(kind=8), dimension(:), allocatable :: x_scf ,y_scf
  !FFT arrays
  real(kind=8), dimension(:,:,:), allocatable :: halfft_cache,kernel
  real(kind=8), dimension(:,:,:,:), allocatable :: kernel_mpi
  real(kind=8), dimension(:,:), allocatable :: cossinarr,btrig
  integer, dimension(:), allocatable :: after,now,before
  
  real(kind=8) :: pi,dx,mu1,ponx,pony
  real(kind=8) :: a,b,c,d,feR,feI,foR,foI,fR,cp,sp,pion,x,value,diff
  integer :: n_scf,ncache,imu,ierr
  integer :: n_range,n_cell,num_of_mus,shift,istart,iend,ireim,jreim,j2st,j2nd,nact2
  integer :: i,i1,i2,i3,i_stat,i_all
  integer :: j2,ind1,ind2,jnd1,ic,inzee,nfft,ipolyord,jp2

  !coefficients for the polynomial interpolation
  real(kind=8), dimension(9,8) :: cpol
  !assign the values of the coefficients  
  cpol(:,:)=1.d0

  cpol(1,2)=.25d0
  
  cpol(1,3)=1.d0/3.d0
  
  cpol(1,4)=7.d0/12.d0
  cpol(2,4)=8.d0/3.d0
  
  cpol(1,5)=19.d0/50.d0
  cpol(2,5)=3.d0/2.d0
  
  cpol(1,6)=41.d0/272.d0
  cpol(2,6)=27.d0/34.d0
  cpol(3,6)=27.d0/272.d0
  
  cpol(1,7)=751.d0/2989.d0
  cpol(2,7)=73.d0/61.d0
  cpol(3,7)=27.d0/61.d0
  
  cpol(1,8)=-989.d0/4540.d0
  cpol(2,8)=-1472.d0/1135.d0
  cpol(3,8)=232.d0/1135.d0
  cpol(4,8)=-2624.d0/1135.d0

  !renormalize values
  cpol(1,1)=.5d0*cpol(1,1)
  cpol(1:2,2)=2.d0/3.d0*cpol(1:2,2)
  cpol(1:2,3)=3.d0/8.d0*cpol(1:2,3)
  cpol(1:3,4)=2.d0/15.d0*cpol(1:3,4)
  cpol(1:3,5)=25.d0/144.d0*cpol(1:3,5)
  cpol(1:4,6)=34.d0/105.d0*cpol(1:4,6)
  cpol(1:4,7)=2989.d0/17280.d0*cpol(1:4,7)
  cpol(1:5,8)=-454.d0/2835.d0*cpol(1:5,8)
  
  !assign the complete values
  cpol(2,1)=cpol(1,1)

  cpol(3,2)=cpol(1,2)

  cpol(3,3)=cpol(2,3)
  cpol(4,3)=cpol(1,3)

  cpol(4,4)=cpol(2,4)
  cpol(5,4)=cpol(1,4)

  cpol(4,5)=cpol(3,5)
  cpol(5,5)=cpol(2,5)
  cpol(6,5)=cpol(1,5)

  cpol(5,6)=cpol(3,6)
  cpol(6,6)=cpol(2,6)
  cpol(7,6)=cpol(1,6)

  cpol(5,7)=cpol(4,7)
  cpol(6,7)=cpol(3,7)
  cpol(7,7)=cpol(2,7)
  cpol(8,7)=cpol(1,7)

  cpol(6,8)=cpol(4,8)
  cpol(7,8)=cpol(3,8)
  cpol(8,8)=cpol(2,8)
  cpol(9,8)=cpol(1,8)

  !Number of integration points : 2*itype_scf*n_points
  n_scf=2*itype_scf*n_points
  !Allocations
  allocate(x_scf(0:n_scf),stat=i_stat)
  call memocc(i_stat,product(shape(x_scf))*kind(x_scf),'x_scf','surfaces_kernel')
  allocate(y_scf(0:n_scf),stat=i_stat)
  call memocc(i_stat,product(shape(y_scf))*kind(y_scf),'y_scf','surfaces_kernel')
  !Build the scaling function
  call scaling_function(itype_scf,n_scf,n_range,x_scf,y_scf)
  !Step grid for the integration
  dx = real(n_range,kind=8)/real(n_scf,kind=8)
  !Extend the range (no more calculations because fill in by 0.d0)
  n_cell = m3
  n_range = max(n_cell,n_range)


  !Allocations
  ncache=ncache_optimal
  !the HalFFT must be performed only in the third dimension,
  !and nker3=n3/2+1, hence
  if (ncache <= (nker3-1)*4) ncache=nker3-1*4

  !enlarge the second dimension of the kernel to be compatible with nproc
  nact2=nker2
  enlarge_ydim: do
     if (nproc*(nact2/nproc) /= nact2) then
        nact2=nact2+1
     else
        exit enlarge_ydim
     end if
  end do enlarge_ydim

  !array for the MPI procedure
  allocate(kernel(nker1,nact2/nproc,nker3),stat=i_stat)
  call memocc(i_stat,product(shape(kernel))*kind(kernel),'kernel','surfaces_kernel')
  allocate(kernel_mpi(nker1,nact2/nproc,nker3/nproc,nproc),stat=i_stat)
  call memocc(i_stat,product(shape(kernel_mpi))*kind(kernel_mpi),'kernel_mpi','surfaces_kernel')
  allocate(kernel_scf(n_range),stat=i_stat)
  call memocc(i_stat,product(shape(kernel_scf))*kind(kernel_scf),'kernel_scf','surfaces_kernel')
  allocate(halfft_cache(2,ncache/4,2),stat=i_stat)
  call memocc(i_stat,product(shape(halfft_cache))*kind(halfft_cache),'halfft_cache','surfaces_kernel')
  allocate(cossinarr(2,n3/2-1),stat=i_stat)
  call memocc(i_stat,product(shape(cossinarr))*kind(cossinarr),'cossinarr','surfaces_kernel')
  allocate(btrig(2,nfft_max),stat=i_stat)
  call memocc(i_stat,product(shape(btrig))*kind(btrig),'btrig','surfaces_kernel')
  allocate(after(7),stat=i_stat)
  call memocc(i_stat,product(shape(after))*kind(after),'after','surfaces_kernel')
  allocate(now(7),stat=i_stat)
  call memocc(i_stat,product(shape(now))*kind(now),'now','surfaces_kernel')
  allocate(before(7),stat=i_stat)
  call memocc(i_stat,product(shape(before))*kind(before),'before','surfaces_kernel')


  !constants
  pi=4.d0*datan(1.d0)

  !arrays for the halFFT
  call ctrig(n3/2,btrig,after,before,now,1,ic)

 
  !build the phases for the HalFFT reconstruction 
  pion=2.d0*pi/real(n3,kind=8)
  do i3=2,n3/2
     x=real(i3-1,kind=8)*pion
     cossinarr(1,i3-1)= dcos(x)
     cossinarr(2,i3-1)=-dsin(x)
  end do
  !kernel=0.d0
  !kernel_mpi=0.d0

  !calculate the limits of the FFT calculations
  !that can be performed in a row remaining inside the cache
  num_of_mus=ncache/(2*n3)

  diff=0.d0
  !order of the polynomial to be used for integration (must be a power of two)
  ipolyord=8 !this part should be incorporated inside the numerical integration
  !here we have to choice the piece of the x-y grid to cover

  !let us now calculate the fraction of mu that will be considered 
  j2st=iproc*(nact2/nproc)
  j2nd=min((iproc+1)*(nact2/nproc),n2/2+1)

  do ind2=(n1/2+1)*j2st+1,(n1/2+1)*j2nd,num_of_mus
     istart=ind2
     iend=min(ind2+(num_of_mus-1),(n1/2+1)*j2nd)
     nfft=iend-istart+1
     shift=0

     !initialization of the interesting part of the cache array
     halfft_cache(:,:,:)=0.d0

     if (istart == 1) then
        !i2=1
        shift=1

        call calculates_green_opt_muzero(n_range,n_scf,ipolyord,x_scf,y_scf,&
             cpol(1,ipolyord),dx,kernel_scf)

        !copy of the first zero value
        halfft_cache(1,1,1)=0.d0

        do i3=1,m3

           value=0.5d0*h3*kernel_scf(i3)
           !index in where to copy the value of the kernel
           call indices(ireim,num_of_mus,n3/2+i3,1,ind1)
           !index in where to copy the symmetric value
           call indices(jreim,num_of_mus,n3/2+2-i3,1,jnd1)
           halfft_cache(ireim,ind1,1) = value
           halfft_cache(jreim,jnd1,1) = value

        end do

     end if

     loopimpulses : do imu=istart+shift,iend

        !here there is the value of mu associated to hgrid
        !note that we have multiplicated mu for hgrid to be comparable 
        !with mu0ref

        !calculate the proper value of mu taking into account the periodic dimensions
        !corresponding value of i1 and i2
        i1=mod(imu,n1/2+1)
        if (i1==0) i1=n1/2+1
        i2=(imu-i1)/(n1/2+1)+1
        ponx=real(i1-1,kind=8)/real(n1,kind=8)
        pony=real(i2-1,kind=8)/real(n2,kind=8)
        
        mu1=2.d0*pi*sqrt((ponx/h1)**2+(pony/h2)**2)*h3

        call calculates_green_opt(n_range,n_scf,itype_scf,ipolyord,x_scf,y_scf,&
             cpol(1,ipolyord),mu1,dx,kernel_scf)

        !readjust the coefficient and define the final kernel

        !copy of the first zero value
        halfft_cache(1,imu-istart+1,1) = 0.d0
        do i3=1,m3
           value=-0.5d0*h3/mu1*kernel_scf(i3)
           !write(80,*)mu1,i3,kernel_scf(i03)
           !index in where to copy the value of the kernel
           call indices(ireim,num_of_mus,n3/2+i3,imu-istart+1,ind1)
           !index in where to copy the symmetric value
           call indices(jreim,num_of_mus,n3/2+2-i3,imu-istart+1,jnd1)
           halfft_cache(ireim,ind1,1)=value
           halfft_cache(jreim,jnd1,1)=value
        end do

     end do loopimpulses

     !now perform the FFT of the array in cache
     inzee=1
     do i=1,ic
        call fftstp(num_of_mus,nfft,n3/2,num_of_mus,n3/2,&
             halfft_cache(1,1,inzee),halfft_cache(1,1,3-inzee),&
             btrig,after(i),now(i),before(i),1)
        inzee=3-inzee
     enddo
     !assign the values of the FFT array
     !and compare with the good results
     do imu=istart,iend

        !corresponding value of i1 and i2
        i1=mod(imu,n1/2+1)
        if (i1==0) i1=n1/2+1
        i2=(imu-i1)/(n1/2+1)+1

        j2=i2-j2st

        a=halfft_cache(1,imu-istart+1,inzee)
        b=halfft_cache(2,imu-istart+1,inzee)
        kernel(i1,j2,1)=a+b
        kernel(i1,j2,n3/2+1)=a-b

        do i3=2,n3/2
           ind1=imu-istart+1+num_of_mus*(i3-1)
           jnd1=imu-istart+1+num_of_mus*(n3/2+2-i3-1)
           cp=cossinarr(1,i3-1)
           sp=cossinarr(2,i3-1)
           a=halfft_cache(1,ind1,inzee)
           b=halfft_cache(2,ind1,inzee)
           c=halfft_cache(1,jnd1,inzee)
           d=halfft_cache(2,jnd1,inzee)
           feR=.5d0*(a+c)
           feI=.5d0*(b-d)
           foR=.5d0*(a-c)
           foI=.5d0*(b+d) 
           fR=feR+cp*foI-sp*foR
           kernel(i1,j2,i3)=fR
        end do
     end do

  end do

  i_all=-product(shape(cossinarr))*kind(cossinarr)
  deallocate(cossinarr,stat=i_stat)
  call memocc(i_stat,i_all,'cossinarr','surfaces_kernel')


  !give to each processor a slice of the third dimension
  if (nproc > 1) then
     call MPI_ALLTOALL(kernel,nker1*(nact2/nproc)*(nker3/nproc), &
          MPI_double_precision, &
          kernel_mpi,nker1*(nact2/nproc)*(nker3/nproc), &
          MPI_double_precision,MPI_COMM_WORLD,ierr)

!!$     !Maximum difference
!!$     max_diff = 0.d0
!!$     i1_max = 1
!!$     i2_max = 1
!!$     i3_max = 1
!!$     do i3=1,nker3/nproc
!!$        do i2=1,nact2/nproc
!!$           do i1=1,nker1
!!$              factor=abs(kernel(i1,i2,i3+iproc*(nker3/nproc))&
!!$                   -kernel_mpi(i1,i2,i3,iproc+1))
!!$              if (max_diff < factor) then
!!$                 max_diff = factor
!!$                 i1_max = i1
!!$                 i2_max = i2
!!$                 i3_max = i3
!!$              end if
!!$           end do
!!$        end do
!!$     end do
!!$     write(*,*) '------------------'
!!$     print *,'iproc=',iproc,'difference post-mpi, at',i1_max,i2_max,i3_max
!!$     write(unit=*,fmt="(1x,a,1pe12.4)") 'Max diff: ',max_diff,&
!!$          'calculated',kernel(i1_max,i2_max,i3_max+iproc*(nker3/nproc)),'post-mpi',kernel_mpi(i1_max,i2_max,i3_max,iproc+1)

     do jp2=1,nproc
        do i3=1,nker3/nproc
           do i2=1,nact2/nproc
              j2=i2+(jp2-1)*(nact2/nproc)
              if (j2 <= nker2) then
                 do i1=1,nker1
                    karray(i1,j2,i3)=&
                         kernel_mpi(i1,i2,i3,jp2)
                 end do
              end if
           end do
        end do
     end do

  else
     karray(1:nker1,1:nker2,1:nker3)=kernel(1:nker1,1:nker2,1:nker3)
  endif


  !De-allocations
  i_all=-product(shape(kernel))*kind(kernel)
  deallocate(kernel,stat=i_stat)
  call memocc(i_stat,i_all,'kernel','surfaces_kernel')
  i_all=-product(shape(kernel_mpi))*kind(kernel_mpi)
  deallocate(kernel_mpi,stat=i_stat)
  call memocc(i_stat,i_all,'kernel_mpi','surfaces_kernel')
  i_all=-product(shape(btrig))*kind(btrig)
  deallocate(btrig,stat=i_stat)
  call memocc(i_stat,i_all,'btrig','surfaces_kernel')
  i_all=-product(shape(after))*kind(after)
  deallocate(after,stat=i_stat)
  call memocc(i_stat,i_all,'after','surfaces_kernel')
  i_all=-product(shape(now))*kind(now)
  deallocate(now,stat=i_stat)
  call memocc(i_stat,i_all,'now','surfaces_kernel')
  i_all=-product(shape(before))*kind(before)
  deallocate(before,stat=i_stat)
  call memocc(i_stat,i_all,'before','surfaces_kernel')
  i_all=-product(shape(halfft_cache))*kind(halfft_cache)
  deallocate(halfft_cache,stat=i_stat)
  call memocc(i_stat,i_all,'halfft_cache','surfaces_kernel')
  i_all=-product(shape(kernel_scf))*kind(kernel_scf)
  deallocate(kernel_scf,stat=i_stat)
  call memocc(i_stat,i_all,'kernel_scf','surfaces_kernel')
  i_all=-product(shape(x_scf))*kind(x_scf)
  deallocate(x_scf,stat=i_stat)
  call memocc(i_stat,i_all,'x_scf','surfaces_kernel')
  i_all=-product(shape(y_scf))*kind(y_scf)
  deallocate(y_scf,stat=i_stat)
  call memocc(i_stat,i_all,'y_scf','surfaces_kernel')

end subroutine Surfaces_Kernel
!!***

subroutine calculates_green_opt(n,n_scf,itype_scf,intorder,xval,yval,c,mu,hres,g_mu)
  implicit none
  real(kind=8), parameter :: mu_max=0.2d0
  integer, intent(in) :: n,n_scf,intorder,itype_scf
  real(kind=8), intent(in) :: hres,mu
  real(kind=8), dimension(0:n_scf), intent(in) :: xval,yval
  real(kind=8), dimension(intorder+1), intent(in) :: c
  real(kind=8), dimension(n), intent(out) :: g_mu
  !local variables
  integer :: izero,ivalue,i,iend,ikern,n_iter,nrec,i_all,i_stat
  real(kind=8) :: f,x,filter,gleft,gright,gltmp,grtmp,fl,fr,x0,x1,ratio,mu0
  real(kind=8), dimension(:), allocatable :: green,green1


  !We calculate the number of iterations to go from mu0 to mu0_ref
  if (mu <= mu_max) then
     n_iter = 0
     mu0 = mu
  else
     n_iter=1
     loop_iter: do
        ratio=real(2**n_iter,kind=8)
        mu0=mu/ratio
        if (mu0 <= mu_max) then
           exit loop_iter
        end if
        n_iter=n_iter+1
     end do loop_iter
  end if

  !dimension needed for the correct calculation of the recursion
  nrec=2**n_iter*n

  allocate(green(-nrec:nrec),stat=i_stat)
  call memocc(i_stat,product(shape(green))*kind(green),'green','calculates_green_opt')


  !initialization of the branching value
  ikern=0
  izero=0
  initialization: do
     if (xval(izero)>=real(ikern,kind=8) .or. izero==n_scf) exit initialization
     izero=izero+1
  end do initialization
  green=0.d0
  !now perform the interpolation in right direction
  ivalue=izero
  gright=0.d0
  loop_right: do
     if(ivalue >= n_scf-intorder-1) exit loop_right
     do i=1,intorder+1
        x=xval(ivalue)-real(ikern,kind=8)
        f=yval(ivalue)*dexp(-mu0*x)
        filter=real(intorder,kind=8)*c(i)
        gright=gright+filter*f
        ivalue=ivalue+1
     end do
     ivalue=ivalue-1
  end do loop_right
  iend=n_scf-ivalue
  do i=1,iend
     x=xval(ivalue)-real(ikern,kind=8)
     f=yval(ivalue)*dexp(-mu0*x)
     filter=real(intorder,kind=8)*c(i)
     gright=gright+filter*f
     ivalue=ivalue+1
  end do
  gright=hres*gright

  !the scaling function is symmetric, so the same for the other direction
  gleft=gright

  green(ikern)=gleft+gright

  !now the loop until the last value
  do ikern=1,nrec
     gltmp=0.d0
     grtmp=0.d0
     ivalue=izero
     x0=xval(izero)
     loop_integration: do
        if (izero==n_scf)  exit loop_integration
        do i=1,intorder+1
           x=xval(ivalue)
           fl=yval(ivalue)*dexp(mu0*x)
           fr=yval(ivalue)*dexp(-mu0*x)
           filter=real(intorder,kind=8)*c(i)
           gltmp=gltmp+filter*fl
           grtmp=grtmp+filter*fr
           ivalue=ivalue+1
           if (xval(izero)>=real(ikern,kind=8) .or. izero==n_scf) then
              x1=xval(izero)
              exit loop_integration
           end if
           izero=izero+1
        end do
        ivalue=ivalue-1
        izero=izero-1
     end do loop_integration
     gleft=dexp(-mu0)*(gleft+hres*dexp(-mu0*real(ikern-1,kind=8))*gltmp)
     if (izero == n_scf) then
        gright=0.d0
     else
        gright=dexp(mu0)*(gright-hres*dexp(mu0*real(ikern-1,kind=8))*grtmp)
     end if
     green(ikern)=gleft+gright
     green(-ikern)=gleft+gright
     if (abs(green(ikern)) <= 1.d-20) then
        nrec=ikern
        exit
     end if
     !print *,ikern,izero,n_scf,gltmp,grtmp,gleft,gright,x0,x1,green(ikern)
  end do
  !now we must calculate the recursion
  allocate(green1(-nrec:nrec),stat=i_stat)
  call memocc(i_stat,product(shape(green1))*kind(green1),'green1','calculates_green_opt')
  !Start the iteration to go from mu0 to mu
  call scf_recursion(itype_scf,n_iter,nrec,green(-nrec),green1(-nrec))

  do i=1,min(n,nrec)
     g_mu(i)=green(i-1)
  end do
  do i=min(n,nrec)+1,n
     g_mu(i)=0.d0
  end do
  

  i_all=-product(shape(green))*kind(green)
  deallocate(green,stat=i_stat)
  call memocc(i_stat,i_all,'green','calculates_green_opt')
  i_all=-product(shape(green1))*kind(green1)
  deallocate(green1,stat=i_stat)
  call memocc(i_stat,i_all,'green1','calculates_green_opt')

end subroutine calculates_green_opt

subroutine calculates_green_opt_muzero(n,n_scf,intorder,xval,yval,c,hres,green)
  implicit none
  integer, intent(in) :: n,n_scf,intorder
  real(kind=8), intent(in) :: hres
  real(kind=8), dimension(0:n_scf), intent(in) :: xval,yval
  real(kind=8), dimension(intorder+1), intent(in) :: c
  real(kind=8), dimension(n), intent(out) :: green
  !local variables
  integer :: izero,ivalue,i,iend,ikern
  real(kind=8) :: x,y,filter,gl0,gl1,gr0,gr1,c0,c1

  !initialization of the branching value
  ikern=0
  izero=0
  initialization: do
     if (xval(izero)>=real(ikern,kind=8) .or. izero==n_scf) exit initialization
     izero=izero+1
  end do initialization
  green=0.d0
  !first case, ikern=0
  !now perform the interpolation in right direction
  ivalue=izero
  gr1=0.d0
  loop_right: do
     if(ivalue >= n_scf-intorder-1) exit loop_right
     do i=1,intorder+1
        x=xval(ivalue)
        y=yval(ivalue)
        filter=real(intorder,kind=8)*c(i)
        gr1=gr1+filter*x*y
        ivalue=ivalue+1
     end do
     ivalue=ivalue-1
  end do loop_right
  iend=n_scf-ivalue
  do i=1,iend
     x=xval(ivalue)
     y=yval(ivalue)
     filter=real(intorder,kind=8)*c(i)
     gr1=gr1+filter*x*y
     ivalue=ivalue+1
  end do
  gr1=hres*gr1
  !the scaling function is symmetric
  gl1=-gr1
  gl0=0.5d0
  gr0=0.5d0

  green(1)=2.d0*gr1

  !now the loop until the last value
  do ikern=1,n-1
     c0=0.d0
     c1=0.d0
     ivalue=izero
     loop_integration: do
        if (izero==n_scf)  exit loop_integration
        do i=1,intorder+1
           x=xval(ivalue)
           y=yval(ivalue)
           filter=real(intorder,kind=8)*c(i)
           c0=c0+filter*y
           c1=c1+filter*y*x
           ivalue=ivalue+1
           if (xval(izero)>=real(ikern,kind=8) .or. izero==n_scf) then
              exit loop_integration
           end if
           izero=izero+1
        end do
        ivalue=ivalue-1
        izero=izero-1
     end do loop_integration
     c0=hres*c0
     c1=hres*c1

     gl0=gl0+c0
     gl1=gl1+c1
     gr0=gr0-c0
     gr1=gr1-c1
     !general case
     green(ikern+1)=real(ikern,kind=8)*(gl0-gr0)+gr1-gl1
     !print *,ikern,izero,n_scf,gltmp,grtmp,gleft,gright,x0,x1,green(ikern)
  end do

end subroutine calculates_green_opt_muzero

subroutine indices(realimag,nelem,intrn,extrn,index)

  implicit none
  integer, intent(in) :: intrn,extrn,nelem
  integer, intent(out) :: realimag,index
  !local
  integer :: i
  !real or imaginary part
  realimag=2-mod(intrn,2)
  !actual index over half the length
  i=(intrn+1)/2
  !check
  if (2*(i-1)+realimag /= intrn) then
     print *,'error, index=',intrn,'realimag=',realimag,'i=',i
  end if
  !complete index to be assigned
  index=extrn+nelem*(i-1)

end subroutine indices


!!****h* BigDFT/Free_Kernel
!! NAME
!!   Free_Kernel
!!
!! FUNCTION
!!    Build the kernel of a gaussian function
!!    for interpolating scaling functions.
!!    Do the parallel HalFFT of the symmetrized function and stores into
!!    memory only 1/8 of the grid divided by the number of processes nproc
!!
!! SYNOPSIS
!!    Build the kernel (karray) of a gaussian function
!!    for interpolating scaling functions
!!    $$ K(j) = \sum_k \omega_k \int \int \phi(x) g_k(x'-x) \delta(x'- j) dx dx' $$
!!
!!    n01,n02,n03        Mesh dimensions of the density
!!    nfft1,nfft2,nfft3  Dimensions of the FFT grid (HalFFT in the third direction)
!!    n1k,n2k,n3k        Dimensions of the kernel FFT
!!    hgrid              Mesh step
!!    itype_scf          Order of the scaling function (8,14,16)
!!
!! AUTHORS
!!    T. Deutsch, L. Genovese
!! CREATION DATE
!!    February 2006
!!
!! SOURCE
!!
subroutine Free_Kernel(n01,n02,n03,nfft1,nfft2,nfft3,n1k,n2k,n3k,&
     hx,hy,hz,itype_scf,iproc,nproc,karray)

 implicit none

 !Arguments
 integer, intent(in) :: n01,n02,n03,nfft1,nfft2,nfft3,n1k,n2k,n3k,itype_scf,iproc,nproc
 real(kind=8), intent(in) :: hx,hy,hz
 real(kind=8), dimension(n1k,n2k,n3k/nproc), intent(out) :: karray

 !Local variables
 !Do not touch !!!!
 integer, parameter :: n_gauss = 89
 !Better if higher (1024 points are enough 10^{-14}: 2*itype_scf*n_points)
 integer, parameter :: n_points = 2**6

 !Better p_gauss for calculation
 !(the support of the exponential should be inside [-n_range/2,n_range/2])
 real(kind=8), parameter :: p0_ref = 1.d0
 real(kind=8), dimension(n_gauss) :: p_gauss,w_gauss

 real(kind=8), dimension(:), allocatable :: kern_1_scf,x_scf ,y_scf
 real(kind=8), dimension(:,:), allocatable :: kernel_scf
 real(kind=8), dimension(:,:,:), allocatable :: kp


 real(kind=8) :: ur_gauss,dr_gauss,acc_gauss,pgauss,kern,a_range
 real(kind=8) :: factor,factor2,dx,absci,p0gauss,p0_cell,u1,u2,u3
 real(kind=8) :: a1,a2,a3,hgrid,pref1,pref2,pref3,p01,p02,p03,kern1,kern2,kern3
 integer :: n_scf,nker1,nker2,nker3
 integer :: i_gauss,n_range,n_cell,istart,iend,istart1
 integer :: i,n_iter,i1,i2,i3,i_kern,i_stat,i_all
 integer :: i01,i02,i03,n1h,n2h,n3h,nit1,nit2,nit3

 !grid spacing
 hgrid=max(hx,hy,hz)
 !Number of integration points : 2*itype_scf*n_points
 n_scf=2*itype_scf*n_points
 !Set karray

 !here we must set the dimensions for the fft part, starting from the nfft
 !remember that actually nfft2 is associated to n03 and viceversa
 
 !dimensions that define the center of symmetry
 n1h=nfft1/2
 n2h=nfft2/2
 n3h=nfft3/2

 !Auxiliary dimensions only for building the FFT part
 nker1=nfft1
 nker2=nfft2
 nker3=nfft3/2+1

 !adjusting the last two dimensions to be multiples of nproc
 do
    if(modulo(nker2,nproc) == 0) exit
    nker2=nker2+1
 end do
 do
    if(modulo(nker3,nproc) == 0) exit
    nker3=nker3+1
 end do

 !this will be the array of the kernel in the real space
 allocate(kp(n1h+1,n3h+1,nker2/nproc),stat=i_stat)
 call memocc(i_stat,product(shape(kp))*kind(kp),'kp','free_kernel')

 !defining proper extremes for the calculation of the
 !local part of the kernel

 istart=iproc*nker2/nproc+1
 iend=min((iproc+1)*nker2/nproc,n2h+n03)

 istart1=max(istart,n2h-n03+2)
 if(iproc .eq. 0) istart1=n2h-n03+2

 !Allocations
 allocate(x_scf(0:n_scf),stat=i_stat)
 call memocc(i_stat,product(shape(x_scf))*kind(x_scf),'x_scf','free_kernel')
 allocate(y_scf(0:n_scf),stat=i_stat)
 call memocc(i_stat,product(shape(y_scf))*kind(y_scf),'y_scf','free_kernel')

 !Build the scaling function
 call scaling_function(itype_scf,n_scf,n_range,x_scf,y_scf)
 !Step grid for the integration
 dx = real(n_range,kind=8)/real(n_scf,kind=8)
 !Extend the range (no more calculations because fill in by 0.d0)
 n_cell = max(n01,n02,n03)
 n_range = max(n_cell,n_range)

 !Lengthes of the box (use FFT dimension)
 a1 = hx * real(n01,kind=8)
 a2 = hy * real(n02,kind=8)
 a3 = hz * real(n03,kind=8)


 !Initialization of the gaussian (Beylkin)
 call gequad(n_gauss,p_gauss,w_gauss,ur_gauss,dr_gauss,acc_gauss)
 !In order to have a range from a_range=sqrt(a1*a1+a2*a2+a3*a3)
 !(biggest length in the cube)
 !We divide the p_gauss by a_range**2 and a_gauss by a_range
 a_range = sqrt(a1*a1+a2*a2+a3*a3)
 factor = 1.d0/a_range
 !factor2 = factor*factor
 factor2 = 1.d0/(a1*a1+a2*a2+a3*a3)
 do i_gauss=1,n_gauss
    p_gauss(i_gauss) = factor2*p_gauss(i_gauss)
 end do
 do i_gauss=1,n_gauss
    w_gauss(i_gauss) = factor*w_gauss(i_gauss)
 end do

 kp(:,:,:)=0.d0

 !Allocations
 allocate(kern_1_scf(-n_range:n_range),stat=i_stat)
 call memocc(i_stat,product(shape(kern_1_scf))*kind(kern_1_scf),'kern_1_scf','free_kernel')
 !add the treatment for inhomogeneous hgrids
 if (hx == hy .and. hy == hz) then
    allocate(kernel_scf(-n_range:n_range,1),stat=i_stat)
    call memocc(i_stat,product(shape(kernel_scf))*kind(kernel_scf),'kernel_scf','free_kernel')

    hgrid=hx

    !To have a correct integration
    p0_cell = p0_ref/(hgrid*hgrid)

    !Use in this order (better for accuracy).
    loop_gauss1: do i_gauss=n_gauss,1,-1
       !Gaussian
       pgauss = p_gauss(i_gauss)
       !We calculate the number of iterations to go from pgauss to p0_ref
       n_iter = nint((log(pgauss) - log(p0_cell))/log(4.d0))
       if (n_iter <= 0)then
          n_iter = 0
          p0gauss = pgauss
       else
          p0gauss = pgauss/4.d0**n_iter
       end if

       !Stupid integration
       !Do the integration with the exponential centered in i_kern
       kernel_scf(:,1) = 0.d0
       do i_kern=0,n_range
          kern = 0.d0
          do i=0,n_scf
             absci = x_scf(i) - real(i_kern,kind=8)
             absci = absci*absci*hgrid**2
             kern = kern + y_scf(i)*dexp(-p0gauss*absci)
          end do
          kernel_scf(i_kern,1) = kern*dx
          kernel_scf(-i_kern,1) = kern*dx
          if (abs(kern) < 1.d-18) then
             !Too small not useful to calculate
             exit
          end if
       end do

       !Start the iteration to go from p0gauss to pgauss
       call scf_recursion(itype_scf,n_iter,n_range,kernel_scf,kern_1_scf)

       !Add to the kernel (only the local part)

       do i3=istart1,iend
          i03 = i3 - n2h -1
          do i2=1,n02
             i02 = i2-1
             do i1=1,n01
                i01 = i1-1
                kp(i1,i2,i3-istart+1) = kp(i1,i2,i3-istart+1) + w_gauss(i_gauss)* &
                     kernel_scf(i01,1)*kernel_scf(i02,1)*kernel_scf(i03,1)
             end do
          end do
       end do


    end do loop_gauss1

 else

    allocate(kernel_scf(-n_range:n_range,3),stat=i_stat)
    call memocc(i_stat,product(shape(kernel_scf))*kind(kernel_scf),'kernel_scf','free_kernel')

    !To have a correct integration
    pref1 = p0_ref/(hx*hx)
    pref2 = p0_ref/(hy*hy)
    pref3 = p0_ref/(hz*hz)

    !Use in this order (better for accuracy).
    loop_gauss: do i_gauss=n_gauss,1,-1
       !Gaussian
       pgauss = p_gauss(i_gauss)
       !We calculate the number of iterations to go from pgauss to p0_ref
       nit1 = max(nint((log(pgauss) - log(pref1))/log(4.d0)),0)
       nit2 = max(nint((log(pgauss) - log(pref2))/log(4.d0)),0)
       nit3 = max(nint((log(pgauss) - log(pref3))/log(4.d0)),0)
       p01=pgauss/4.d0**nit1
       p02=pgauss/4.d0**nit2
       p03=pgauss/4.d0**nit3

       !Stupid integration
       !Do the integration with the exponential centered in i_kern
       kernel_scf(:,:) = 0.d0
       do i_kern=0,n_range
          kern1=0.d0
          kern2=0.d0
          kern3=0.d0
          do i=0,n_scf
             absci = x_scf(i) - real(i_kern,kind=8)
             u1=-p01*absci*absci*hx**2
             u2=-p02*absci*absci*hy**2
             u3=-p03*absci*absci*hz**2
             u1=dexp(u1)
             u2=dexp(u2)
             u3=dexp(u3)
             kern1=kern1+y_scf(i)*u1
             kern2=kern2+y_scf(i)*u2
             kern3=kern3+y_scf(i)*u3
          end do
          kernel_scf(i_kern,1) = kern1*dx
          kernel_scf(-i_kern,1) = kern1*dx
          kernel_scf(i_kern,2) = kern2*dx
          kernel_scf(-i_kern,2) = kern2*dx
          kernel_scf(i_kern,3) = kern3*dx
          kernel_scf(-i_kern,3) = kern3*dx
          if (abs(kern1)+abs(kern2)+abs(kern3) < 3.d-18) then
             !Too small not useful to calculate
             exit
          end if
       end do

       !Start the iteration to go from p0gauss to pgauss
       call scf_recursion(itype_scf,nit1,n_range,kernel_scf(-n_range,1),kern_1_scf)
       call scf_recursion(itype_scf,nit2,n_range,kernel_scf(-n_range,2),kern_1_scf)
       call scf_recursion(itype_scf,nit3,n_range,kernel_scf(-n_range,3),kern_1_scf)

       !Add to the kernel (only the local part)

       do i3=istart1,iend
          i03 = i3 - n2h -1
          do i2=1,n02
             i02 = i2-1
             do i1=1,n01
                i01 = i1-1
                kp(i1,i2,i3-istart+1) = kp(i1,i2,i3-istart+1) + w_gauss(i_gauss)* &
                     kernel_scf(i01,1)*kernel_scf(i02,2)*kernel_scf(i03,3)
             end do
          end do
       end do

    end do loop_gauss

 end if

 !De-allocations
 i_all=-product(shape(kernel_scf))*kind(kernel_scf)
 deallocate(kernel_scf,stat=i_stat)
 call memocc(i_stat,i_all,'kernel_scf','free_kernel')
 i_all=-product(shape(kern_1_scf))*kind(kern_1_scf)
 deallocate(kern_1_scf,stat=i_stat)
 call memocc(i_stat,i_all,'kern_1_scf','free_kernel')
 i_all=-product(shape(x_scf))*kind(x_scf)
 deallocate(x_scf,stat=i_stat)
 call memocc(i_stat,i_all,'x_scf','free_kernel')
 i_all=-product(shape(y_scf))*kind(y_scf)
 deallocate(y_scf,stat=i_stat)
 call memocc(i_stat,i_all,'y_scf','free_kernel')

!!!!END KERNEL CONSTRUCTION

!!$ if(iproc .eq. 0) print *,"Do a 3D PHalFFT for the kernel"

 call kernelfft(nfft1,nfft2,nfft3,nker1,nker2,nker3,n1k,n2k,n3k,nproc,iproc,&
      kp,karray)

 !De-allocations
 i_all=-product(shape(kp))*kind(kp)
 deallocate(kp,stat=i_stat)
 call memocc(i_stat,i_all,'kp','free_kernel')

end subroutine Free_Kernel
!!***

subroutine inserthalf(n1,n3,lot,nfft,i1,zf,zw)
  implicit none
  integer, intent(in) :: n1,n3,lot,nfft,i1
  real(kind=8), dimension(n1/2+1,n3/2+1), intent(in) :: zf
  real(kind=8), dimension(2,lot,n3/2), intent(out) :: zw
  !local variables
  integer :: l1,l3,i01,i03r,i03i,i3

  i3=0
  do l3=1,n3,2
     i3=i3+1
     i03r=abs(l3-n3/2-1)+1
     i03i=abs(l3-n3/2)+1
     do l1=1,nfft
        i01=abs(l1-1+i1-n1/2-1)+1
        zw(1,l1,i3)=zf(i01,i03r)
        zw(2,l1,i3)=zf(i01,i03i)
     end do
  end do

end subroutine inserthalf


!!****h* BigDFT/kernelfft
!! NAME
!!   kernelfft
!!
!! FUNCTION
!!     (Based on suitable modifications of S.Goedecker routines)
!!     Calculates the FFT of the distributed kernel
!!
!! SYNOPSIS
!!     zf:          Real kernel (input)
!!                  zf(i1,i2,i3)
!!     zr:          Distributed Kernel FFT 
!!                  zr(2,i1,i2,i3)
!!     nproc:       number of processors used as returned by MPI_COMM_SIZE
!!     iproc:       [0:nproc-1] number of processor as returned by MPI_COMM_RANK
!!     n1,n2,n3:    logical dimension of the transform. As transform lengths 
!!                  most products of the prime factors 2,3,5 are allowed.
!!                  The detailed table with allowed transform lengths can 
!!                  be found in subroutine CTRIG
!!     nd1,nd2,nd3: Dimensions of work arrays
!!
!! RESTRICTIONS on USAGE
!!     Copyright (C) Stefan Goedecker, Cornell University, Ithaca, USA, 1994
!!     Copyright (C) Stefan Goedecker, MPI Stuttgart, Germany, 1999
!!     Copyright (C) 2002 Stefan Goedecker, CEA Grenoble
!!     This file is distributed under the terms of the
!!     GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .
!!
!! AUTHORS
!!    S. Goedecker, L. Genovese
!!
!! CREATION DATE
!!     February 2006
!!
!! SOURCE
!!
subroutine kernelfft(n1,n2,n3,nd1,nd2,nd3,nk1,nk2,nk3,nproc,iproc,zf,zr)
  use mpi
  implicit none
  include 'perfdata.inc'
  !Arguments
  integer, intent(in) :: n1,n2,n3,nd1,nd2,nd3,nk1,nk2,nk3,nproc,iproc
  real(kind=8), dimension(n1/2+1,n3/2+1,nd2/nproc), intent(in) :: zf
  real(kind=8), dimension(nk1,nk2,nk3/nproc), intent(inout) :: zr
  !Local variables
  !Maximum number of points for FFT (should be same number in fft3d routine)
  integer, parameter :: nfft_max=24000
  integer :: ncache,lzt,lot,ma,mb,nfft,ic1,ic2,ic3,Jp2st,J2st
  integer :: j2,j3,i1,i3,i,j,inzee,ierr,i_all,i_stat
  real(kind=8) :: twopion
  !work arrays for transpositions
  real(kind=8), dimension(:,:,:), allocatable :: zt
  !work arrays for MPI
  real(kind=8), dimension(:,:,:,:,:), allocatable :: zmpi1
  real(kind=8), dimension(:,:,:,:), allocatable :: zmpi2
  !cache work array
  real(kind=8), dimension(:,:,:), allocatable :: zw
  !FFT work arrays
  real(kind=8), dimension(:,:), allocatable :: trig1,trig2,trig3,cosinarr
  integer, dimension(:), allocatable :: after1,now1,before1, & 
       after2,now2,before2,after3,now3,before3

  !Body

  !check input
  if (nd1.lt.n1) stop 'ERROR:nd1'
  if (nd2.lt.n2) stop 'ERROR:nd2'
  if (nd3.lt.n3/2+1) stop 'ERROR:nd3'
  if (mod(nd3,nproc).ne.0) stop 'ERROR:nd3'
  if (mod(nd2,nproc).ne.0) stop 'ERROR:nd2'
  
  !defining work arrays dimensions
  ncache=ncache_optimal
  if (ncache <= max(n1,n2,n3/2)*4) ncache=max(n1,n2,n3/2)*4
  lzt=n2
  if (mod(n2,2).eq.0) lzt=lzt+1
  if (mod(n2,4).eq.0) lzt=lzt+1
  
  !Allocations
  allocate(trig1(2,nfft_max),stat=i_stat)
  call memocc(i_stat,product(shape(trig1))*kind(trig1),'trig1','kernelfft')
  allocate(after1(7),stat=i_stat)
  call memocc(i_stat,product(shape(after1))*kind(after1),'after1','kernelfft')
  allocate(now1(7),stat=i_stat)
  call memocc(i_stat,product(shape(now1))*kind(now1),'now1','kernelfft')
  allocate(before1(7),stat=i_stat)
  call memocc(i_stat,product(shape(before1))*kind(before1),'before1','kernelfft')
  allocate(trig2(2,nfft_max),stat=i_stat)
  call memocc(i_stat,product(shape(trig2))*kind(trig2),'trig2','kernelfft')
  allocate(after2(7),stat=i_stat)
  call memocc(i_stat,product(shape(after2))*kind(after2),'after2','kernelfft')
  allocate(now2(7),stat=i_stat)
  call memocc(i_stat,product(shape(now2))*kind(now2),'now2','kernelfft')
  allocate(before2(7),stat=i_stat)
  call memocc(i_stat,product(shape(before2))*kind(before2),'before2','kernelfft')
  allocate(trig3(2,nfft_max),stat=i_stat)
  call memocc(i_stat,product(shape(trig3))*kind(trig3),'trig3','kernelfft')
  allocate(after3(7),stat=i_stat)
  call memocc(i_stat,product(shape(after3))*kind(after3),'after3','kernelfft')
  allocate(now3(7),stat=i_stat)
  call memocc(i_stat,product(shape(now3))*kind(now3),'now3','kernelfft')
  allocate(before3(7),stat=i_stat)
  call memocc(i_stat,product(shape(before3))*kind(before3),'before3','kernelfft')
  allocate(zw(2,ncache/4,2),stat=i_stat)
  call memocc(i_stat,product(shape(zw))*kind(zw),'zw','kernelfft')
  allocate(zt(2,lzt,n1),stat=i_stat)
  call memocc(i_stat,product(shape(zt))*kind(zt),'zt','kernelfft')
  allocate(zmpi2(2,n1,nd2/nproc,nd3),stat=i_stat)
  call memocc(i_stat,product(shape(zmpi2))*kind(zmpi2),'zmpi2','kernelfft')
  allocate(cosinarr(2,n3/2),stat=i_stat)
  call memocc(i_stat,product(shape(cosinarr))*kind(cosinarr),'cosinarr','kernelfft')
  if (nproc.gt.1) then
     allocate(zmpi1(2,n1,nd2/nproc,nd3/nproc,nproc),stat=i_stat)
     call memocc(i_stat,product(shape(zmpi1))*kind(zmpi1),'zmpi1','kernelfft')
  end if

  
  !calculating the FFT work arrays (beware on the HalFFT in n3 dimension)
  call ctrig(n3/2,trig3,after3,before3,now3,1,ic3)
  call ctrig(n1,trig1,after1,before1,now1,1,ic1)
  call ctrig(n2,trig2,after2,before2,now2,1,ic2)
  
  !Calculating array of phases for HalFFT decoding
  twopion=8.d0*datan(1.d0)/real(n3,kind=8)
  do i3=1,n3/2
     cosinarr(1,i3)= dcos(twopion*real(i3-1,kind=8))
     cosinarr(2,i3)=-dsin(twopion*real(i3-1,kind=8))
  end do
  
  !transform along z axis

  lot=ncache/(2*n3)
  if (lot.lt.1) stop 'kernelfft:enlarge ncache for z'
  
  do j2=1,nd2/nproc
     !this condition ensures that we manage only the interesting part for the FFT
     if (iproc*(nd2/nproc)+j2.le.n2) then
        do i1=1,n1,lot
           ma=i1
           mb=min(i1+(lot-1),n1)
           nfft=mb-ma+1

           !inserting real data into complex array of half lenght
           !input: I1,I3,J2,(Jp2)

           call inserthalf(n1,n3,lot,nfft,i1,zf(1,1,j2),zw(1,1,1))

           !performing FFT
           inzee=1
           do i=1,ic3
              call fftstp(lot,nfft,n3/2,lot,n3/2,zw(1,1,inzee),zw(1,1,3-inzee), &
                   trig3,after3(i),now3(i),before3(i),1)
              inzee=3-inzee
           enddo
           !output: I1,i3,J2,(Jp2)

           !unpacking FFT in order to restore correct result, 
           !while exchanging components
           !input: I1,i3,J2,(Jp2)
           call scramble_unpack(i1,j2,lot,nfft,n1,n3,nd2,nproc,nd3,zw(1,1,inzee),zmpi2,cosinarr)
           !output: I1,J2,i3,(Jp2)
        end do
     endif
  end do

  !Interprocessor data transposition
  !input: I1,J2,j3,jp3,(Jp2)
  if (nproc.gt.1) then
     !communication scheduling
     call MPI_ALLTOALL(zmpi2,2*n1*(nd2/nproc)*(nd3/nproc), &
          MPI_double_precision, &
          zmpi1,2*n1*(nd2/nproc)*(nd3/nproc), &
          MPI_double_precision,MPI_COMM_WORLD,ierr)
     ! output: I1,J2,j3,Jp2,(jp3)
  endif


  do j3=1,nd3/nproc
     !this condition ensures that we manage only the interesting part for the FFT
     if (iproc*(nd3/nproc)+j3.le.n3/2+1) then
        Jp2st=1
        J2st=1
        
        !transform along x axis
        lot=ncache/(4*n1)
        if (lot.lt.1) stop 'kernelfft:enlarge ncache for x'
        
        do j=1,n2,lot
           ma=j
           mb=min(j+(lot-1),n2)
           nfft=mb-ma+1

           !reverse ordering
           !input: I1,J2,j3,Jp2,(jp3)
           if (nproc.eq.1) then
              call mpiswitch(j3,nfft,Jp2st,J2st,lot,n1,nd2,nd3,nproc,zmpi2,zw(1,1,1))
           else
              call mpiswitch(j3,nfft,Jp2st,J2st,lot,n1,nd2,nd3,nproc,zmpi1,zw(1,1,1))
           endif
           !output: J2,Jp2,I1,j3,(jp3)

           !performing FFT
           !input: I2,I1,j3,(jp3)          
           inzee=1
           do i=1,ic1-1
              call fftstp(lot,nfft,n1,lot,n1,zw(1,1,inzee),zw(1,1,3-inzee), &
                   trig1,after1(i),now1(i),before1(i),1)
              inzee=3-inzee
           enddo
           !storing the last step into zt
           i=ic1
           call fftstp(lot,nfft,n1,lzt,n1,zw(1,1,inzee),zt(1,j,1), & 
                trig1,after1(i),now1(i),before1(i),1)
           !output: I2,i1,j3,(jp3)
        end do

        !transform along y axis, and taking only the first half
        lot=ncache/(4*n2)
        if (lot.lt.1) stop 'kernelfft:enlarge ncache for y'

        do j=1,nk1,lot
           ma=j
           mb=min(j+(lot-1),nk1)
           nfft=mb-ma+1

           !reverse ordering
           !input: I2,i1,j3,(jp3)
           call switch(nfft,n2,lot,n1,lzt,zt(1,1,j),zw(1,1,1))
           !output: i1,I2,j3,(jp3)

           !performing FFT
           !input: i1,I2,j3,(jp3)
           inzee=1
           do i=1,ic2
              call fftstp(lot,nfft,n2,lot,n2,zw(1,1,inzee),zw(1,1,3-inzee), &
                   trig2,after2(i),now2(i),before2(i),1)
              inzee=3-inzee
           enddo

           call realcopy(lot,nfft,n2,nk1,nk2,zw(1,1,inzee),zr(j,1,j3))
          
        end do
        !output: i1,i2,j3,(jp3)
     endif
  end do

  !De-allocations
  i_all=-product(shape(trig1))*kind(trig1)
  deallocate(trig1,stat=i_stat)
  call memocc(i_stat,i_all,'trig1','kernelfft')
  i_all=-product(shape(after1))*kind(after1)
  deallocate(after1,stat=i_stat)
  call memocc(i_stat,i_all,'after1','kernelfft')
  i_all=-product(shape(now1))*kind(now1)
  deallocate(now1,stat=i_stat)
  call memocc(i_stat,i_all,'now1','kernelfft')
  i_all=-product(shape(before1))*kind(before1)
  deallocate(before1,stat=i_stat)
  call memocc(i_stat,i_all,'before1','kernelfft')
  i_all=-product(shape(trig2))*kind(trig2)
  deallocate(trig2,stat=i_stat)
  call memocc(i_stat,i_all,'trig2','kernelfft')
  i_all=-product(shape(after2))*kind(after2)
  deallocate(after2,stat=i_stat)
  call memocc(i_stat,i_all,'after2','kernelfft')
  i_all=-product(shape(now2))*kind(now2)
  deallocate(now2,stat=i_stat)
  call memocc(i_stat,i_all,'now2','kernelfft')
  i_all=-product(shape(before2))*kind(before2)
  deallocate(before2,stat=i_stat)
  call memocc(i_stat,i_all,'before2','kernelfft')
  i_all=-product(shape(trig3))*kind(trig3)
  deallocate(trig3,stat=i_stat)
  call memocc(i_stat,i_all,'trig3','kernelfft')
  i_all=-product(shape(after3))*kind(after3)
  deallocate(after3,stat=i_stat)
  call memocc(i_stat,i_all,'after3','kernelfft')
  i_all=-product(shape(now3))*kind(now3)
  deallocate(now3,stat=i_stat)
  call memocc(i_stat,i_all,'now3','kernelfft')
  i_all=-product(shape(before3))*kind(before3)
  deallocate(before3,stat=i_stat)
  call memocc(i_stat,i_all,'before3','kernelfft')
  i_all=-product(shape(zmpi2))*kind(zmpi2)
  deallocate(zmpi2,stat=i_stat)
  call memocc(i_stat,i_all,'zmpi2','kernelfft')
  i_all=-product(shape(zw))*kind(zw)
  deallocate(zw,stat=i_stat)
  call memocc(i_stat,i_all,'zw','kernelfft')
  i_all=-product(shape(zt))*kind(zt)
  deallocate(zt,stat=i_stat)
  call memocc(i_stat,i_all,'zt','kernelfft')
  i_all=-product(shape(cosinarr))*kind(cosinarr)
  deallocate(cosinarr,stat=i_stat)
  call memocc(i_stat,i_all,'cosinarr','kernelfft')
  if (nproc.gt.1) then
     i_all=-product(shape(zmpi1))*kind(zmpi1)
     deallocate(zmpi1,stat=i_stat)
     call memocc(i_stat,i_all,'zmpi1','kernelfft')
  end if

end subroutine kernelfft


subroutine realcopy(lot,nfft,n2,nk1,nk2,zin,zout)
  implicit none
  integer, intent(in) :: nfft,lot,n2,nk1,nk2
  real(kind=8), dimension(2,lot,n2), intent(in) :: zin
  real(kind=8), dimension(nk1,nk2), intent(out) :: zout
  !local variables
  integer :: i,j
 
  do i=1,nk2
     do j=1,nfft
        zout(j,i)=zin(1,j,i)
     end do
  end do

end subroutine realcopy


subroutine switch(nfft,n2,lot,n1,lzt,zt,zw)
  implicit real(kind=8) (a-h,o-z)
  dimension zw(2,lot,n2),zt(2,lzt,n1)

  do j=1,nfft
     do i=1,n2
        zw(1,j,i)=zt(1,i,j)
        zw(2,j,i)=zt(2,i,j)
     end do
  end do

end subroutine switch


subroutine mpiswitch(j3,nfft,Jp2st,J2st,lot,n1,nd2,nd3,nproc,zmpi1,zw)
  implicit real(kind=8) (a-h,o-z)
  dimension zmpi1(2,n1,nd2/nproc,nd3/nproc,nproc),zw(2,lot,n1)

  mfft=0
  do Jp2=Jp2st,nproc
     do J2=J2st,nd2/nproc
        mfft=mfft+1
        if (mfft.gt.nfft) then
           Jp2st=Jp2
           J2st=J2
           return
        endif
        do I1=1,n1
           zw(1,mfft,I1)=zmpi1(1,I1,J2,j3,Jp2)
           zw(2,mfft,I1)=zmpi1(2,I1,J2,j3,Jp2)
        end do
     end do
     J2st=1
  end do

end subroutine mpiswitch


subroutine gequad(nterms,p,w,urange,drange,acc)
! 
  implicit real(kind=8) (a-h,o-z)
  real(kind=8) :: p(*),w(*)
!
!       range [10^(-9),1] and accuracy ~10^(-8);
!
  p(1)=4.96142640560223544d19
  p(2)=1.37454269147978052d19
  p(3)=7.58610013441204679d18
  p(4)=4.42040691347806996d18
  p(5)=2.61986077948367892d18
  p(6)=1.56320138155496681d18
  p(7)=9.35645215863028402d17
  p(8)=5.60962910452691703d17
  p(9)=3.3666225119686761d17
  p(10)=2.0218253197947866d17
  p(11)=1.21477756091902017d17
  p(12)=7.3012982513608503d16
  p(13)=4.38951893556421099d16
  p(14)=2.63949482512262325d16
  p(15)=1.58742054072786174d16
  p(16)=9.54806587737665531d15
  p(17)=5.74353712364571709d15
  p(18)=3.455214877389445d15
  p(19)=2.07871658520326804d15
  p(20)=1.25064667315629928d15
  p(21)=7.52469429541933745d14
  p(22)=4.5274603337253175d14
  p(23)=2.72414006900059548d14
  p(24)=1.63912168349216752d14
  p(25)=9.86275802590865738d13
  p(26)=5.93457701624974985d13
  p(27)=3.5709554322296296d13
  p(28)=2.14872890367310454d13
  p(29)=1.29294719957726902d13
  p(30)=7.78003375426361016d12
  p(31)=4.68148199759876704d12
  p(32)=2.8169955024829868d12
  p(33)=1.69507790481958464d12
  p(34)=1.01998486064607581d12
  p(35)=6.13759486539856459d11
  p(36)=3.69320183828682544d11
  p(37)=2.22232783898905102d11
  p(38)=1.33725247623668682d11
  p(39)=8.0467192739036288d10
  p(40)=4.84199582415144143d10
  p(41)=2.91360091170559564d10
  p(42)=1.75321747475309216d10
  p(43)=1.0549735552210995d10
  p(44)=6.34815321079006586d9
  p(45)=3.81991113733594231d9
  p(46)=2.29857747533101109d9
  p(47)=1.38313653595483694d9
  p(48)=8.32282908580025358d8
  p(49)=5.00814519374587467d8
  p(50)=3.01358090773319025d8
  p(51)=1.81337994217503535d8
  p(52)=1.09117589961086823d8
  p(53)=6.56599771718640323d7
  p(54)=3.95099693638497164d7
  p(55)=2.37745694710665991d7
  p(56)=1.43060135285912813d7
  p(57)=8.60844290313506695d6
  p(58)=5.18000974075383424d6
  p(59)=3.116998193057466d6
  p(60)=1.87560993870024029d6
  p(61)=1.12862197183979562d6
  p(62)=679132.441326077231d0
  p(63)=408658.421279877969d0
  p(64)=245904.473450669789d0
  p(65)=147969.568088321005d0
  p(66)=89038.612357311147d0
  p(67)=53577.7362552358895d0
  p(68)=32239.6513926914668d0
  p(69)=19399.7580852362791d0
  p(70)=11673.5323603058634d0
  p(71)=7024.38438577707758d0
  p(72)=4226.82479307685999d0
  p(73)=2543.43254175354295d0
  p(74)=1530.47486269122675d0
  p(75)=920.941785160749482d0
  p(76)=554.163803906291646d0
  p(77)=333.46029740785694d0
  p(78)=200.6550575335041d0
  p(79)=120.741366914147284d0
  p(80)=72.6544243200329916d0
  p(81)=43.7187810415471025d0
  p(82)=26.3071631447061043d0
  p(83)=15.8299486353816329d0
  p(84)=9.52493152341244004d0
  p(85)=5.72200417067776041d0
  p(86)=3.36242234070940928d0
  p(87)=1.75371394604499472d0
  p(88)=0.64705932650658966d0
  p(89)=0.072765905943708247d0
!
  w(1)=47.67445484528304247d10
  w(2)=11.37485774750442175d9
  w(3)=78.64340976880190239d8
  w(4)=46.27335788759590498d8
  w(5)=24.7380464827152951d8
  w(6)=13.62904116438987719d8
  w(7)=92.79560029045882433d8
  w(8)=52.15931216254660251d8
  w(9)=31.67018011061666244d8
  w(10)=1.29291036801493046d8
  w(11)=1.00139319988015862d8
  w(12)=7.75892350510188341d7
  w(13)=6.01333567950731271d7
  w(14)=4.66141178654796875d7
  w(15)=3.61398903394911448d7
  w(16)=2.80225846672956389d7
  w(17)=2.1730509180930247d7
  w(18)=1.68524482625876965d7
  w(19)=1.30701489345870338d7
  w(20)=1.01371784832269282d7
  w(21)=7.86264116300379329d6
  w(22)=6.09861667912273717d6
  w(23)=4.73045784039455683d6
  w(24)=3.66928949951594161d6
  w(25)=2.8462050836230259d6
  w(26)=2.20777394798527011d6
  w(27)=1.71256191589205524d6
  w(28)=1.32843556197737076d6
  w(29)=1.0304731275955989d6
  w(30)=799345.206572271448d0
  w(31)=620059.354143595343d0
  w(32)=480986.704107449333d0
  w(33)=373107.167700228515d0
  w(34)=289424.08337412132d0
  w(35)=224510.248231581788d0
  w(36)=174155.825690028966d0
  w(37)=135095.256919654065d0
  w(38)=104795.442776800312d0
  w(39)=81291.4458222430418d0
  w(40)=63059.0493649328682d0
  w(41)=48915.9040455329689d0
  w(42)=37944.8484018048756d0
  w(43)=29434.4290473253969d0
  w(44)=22832.7622054490044d0
  w(45)=17711.743950151233d0
  w(46)=13739.287867104177d0
  w(47)=10657.7895710752585d0
  w(48)=8267.42141053961834d0
  w(49)=6413.17397520136448d0
  w(50)=4974.80402838654277d0
  w(51)=3859.03698188553047d0
  w(52)=2993.51824493299154d0
  w(53)=2322.1211966811754d0
  w(54)=1801.30750964719641d0
  w(55)=1397.30379659817038d0
  w(56)=1083.91149143250697d0
  w(57)=840.807939169209188d0
  w(58)=652.228524366749422d0
  w(59)=505.944376983506128d0
  w(60)=392.469362317941064d0
  w(61)=304.444930257324312d0
  w(62)=236.162932842453601d0
  w(63)=183.195466078603525d0
  w(64)=142.107732186551471d0
  w(65)=110.23530215723992d0
  w(66)=85.5113346705382257d0
  w(67)=66.3325469806696621d0
  w(68)=51.4552463353841373d0
  w(69)=39.9146798429449273d0
  w(70)=30.9624728409162095d0
  w(71)=24.018098812215013d0
  w(72)=18.6312338024296588d0
  w(73)=14.4525541233150501d0
  w(74)=11.2110836519105938d0
  w(75)=8.69662175848497178d0
  w(76)=6.74611236165731961d0
  w(77)=5.23307018057529994d0
  w(78)=4.05937850501539556d0
  w(79)=3.14892659076635714d0
  w(80)=2.44267408211071604d0
  w(81)=1.89482240522855261d0
  w(82)=1.46984505907050079d0
  w(83)=1.14019261330527007d0
  w(84)=0.884791217422925293d0
  w(85)=0.692686387080616483d0
  w(86)=0.585244576897023282d0
  w(87)=0.576182522545327589d0
  w(88)=0.596688817388997178d0
  w(89)=0.607879901151108771d0
!
  urange = 1.d0
  drange=1d-08
  acc   =1d-08
!
end subroutine gequad
