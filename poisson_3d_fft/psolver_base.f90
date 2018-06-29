
!! Copyright (C) 2002-2007 BigDFT group 
!! This file is distributed under the terms of the
!! GNU General Public License, see ~/COPYING file
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the list of contributors, see ~/AUTHORS 


subroutine P_PoissonSolver(n1,n2,n3,nd1,nd2,nd3,md1,md2,md3,nproc,iproc,zf&
             ,scal,hx,hy,hz)
  use mpi
  implicit none
  !to be preprocessed
  include 'perfdata.inc'
  !Arguments
  integer, intent(in) :: n1,n2,n3,nd1,nd2,nd3,md1,md2,md3,nproc,iproc
  real(kind=8), intent(in) :: scal,hx,hy,hz
  real(kind=8), dimension(md1,md3,md2/nproc), intent(inout) :: zf
  !Local variables
  !Maximum number of points for FFT (should be same number in fft3d routine)
  integer, parameter :: nfft_max=24000
  integer :: ncache,lzt,lot,ma,mb,nfft,ic1,ic2,ic3,Jp2stb,J2stb,Jp2stf,J2stf
  integer :: j2,j3,i1,i3,i,j,inzee,ierr,i_all,i_stat
  !work arrays for transpositions
  real(kind=8), dimension(:,:,:), allocatable :: zt
  !work arrays for MPI
  real(kind=8), dimension(:,:,:,:,:), allocatable :: zmpi1
  real(kind=8), dimension(:,:,:,:), allocatable :: zmpi2
  !cache work array
  real(kind=8), dimension(:,:,:), allocatable :: zw
  !FFT work arrays
  real(kind=8), dimension(:,:), allocatable :: btrig1,btrig2,btrig3, &
       ftrig1,ftrig2,ftrig3
  integer, dimension(:), allocatable :: after1,now1,before1, & 
       after2,now2,before2,after3,now3,before3
  
  !Body

  call timing(iproc,'PSolv_comput  ','ON')
  ! check input
!!$  if (mod(n1,2).ne.0) stop 'Parallel convolution:ERROR:n1' !this can be avoided
!!$  if (mod(n2,2).ne.0) stop 'Parallel convolution:ERROR:n2' !this can be avoided
!!$  if (mod(n3,2).ne.0) stop 'Parallel convolution:ERROR:n3' !this can be avoided
  if (nd1.lt.n1/2+1) stop 'Parallel convolution:ERROR:nd1' 
  if (nd2.lt.n2/2+1) stop 'Parallel convolution:ERROR:nd2' 
  if (nd3.lt.n3/2+1) stop 'Parallel convolution:ERROR:nd3' 
  if (md1.lt.n1) stop 'Parallel convolution:ERROR:md1'
  if (md2.lt.n2) stop 'Parallel convolution:ERROR:md2'
  if (md3.lt.n3) stop 'Parallel convolution:ERROR:md3'
  if (mod(nd3,nproc).ne.0) stop 'Parallel convolution:ERROR:nd3'
  if (mod(md2,nproc).ne.0) stop 'Parallel convolution:ERROR:md2'
  
  !defining work arrays dimensions
  ncache=ncache_optimal
  if (ncache <= max(n1,n2,n3)*4) ncache=max(n1,n2,n3)*4

  if (timing_flag == 1 .and. iproc ==0) print *,'parallel ncache=',ncache

  lzt=n2
  if (mod(n2,2).eq.0) lzt=lzt+1
  if (mod(n2,4).eq.0) lzt=lzt+1 !maybe this is useless
  
  !Allocations
  allocate(btrig1(2,nfft_max),stat=i_stat)
  call memocc(i_stat,product(shape(btrig1))*kind(btrig1),'btrig1','p_poissonsolver')
  allocate(ftrig1(2,nfft_max),stat=i_stat)
  call memocc(i_stat,product(shape(ftrig1))*kind(ftrig1),'ftrig1','p_poissonsolver')
  allocate(after1(7),stat=i_stat)
  call memocc(i_stat,product(shape(after1))*kind(after1),'after1','p_poissonsolver')
  allocate(now1(7),stat=i_stat)
  call memocc(i_stat,product(shape(now1))*kind(now1),'now1','p_poissonsolver')
  allocate(before1(7),stat=i_stat)
  call memocc(i_stat,product(shape(before1))*kind(before1),'before1','p_poissonsolver')
  allocate(btrig2(2,nfft_max),stat=i_stat)
  call memocc(i_stat,product(shape(btrig2))*kind(btrig2),'btrig2','p_poissonsolver')
  allocate(ftrig2(2,nfft_max),stat=i_stat)
  call memocc(i_stat,product(shape(ftrig2))*kind(ftrig2),'ftrig2','p_poissonsolver')
  allocate(after2(7),stat=i_stat)
  call memocc(i_stat,product(shape(after2))*kind(after2),'after2','p_poissonsolver')
  allocate(now2(7),stat=i_stat)
  call memocc(i_stat,product(shape(now2))*kind(now2),'now2','p_poissonsolver')
  allocate(before2(7),stat=i_stat)
  call memocc(i_stat,product(shape(before2))*kind(before2),'before2','p_poissonsolver')
  allocate(btrig3(2,nfft_max),stat=i_stat)
  call memocc(i_stat,product(shape(btrig3))*kind(btrig3),'btrig3','p_poissonsolver')
  allocate(ftrig3(2,nfft_max),stat=i_stat)
  call memocc(i_stat,product(shape(ftrig3))*kind(ftrig3),'ftrig3','p_poissonsolver')
  allocate(after3(7),stat=i_stat)
  call memocc(i_stat,product(shape(after3))*kind(after3),'after3','p_poissonsolver')
  allocate(now3(7),stat=i_stat)
  call memocc(i_stat,product(shape(now3))*kind(now3),'now3','p_poissonsolver')
  allocate(before3(7),stat=i_stat)
  call memocc(i_stat,product(shape(before3))*kind(before3),'before3','p_poissonsolver')
  allocate(zw(2,ncache/4,2),stat=i_stat)
  call memocc(i_stat,product(shape(zw))*kind(zw),'zw','p_poissonsolver')
  allocate(zt(2,lzt,n1),stat=i_stat)
  call memocc(i_stat,product(shape(zt))*kind(zt),'zt','p_poissonsolver')
  allocate(zmpi2(2,n1,md2/nproc,nd3),stat=i_stat)
  call memocc(i_stat,product(shape(zmpi2))*kind(zmpi2),'zmpi2','p_poissonsolver')
  if (nproc.gt.1) then
     allocate(zmpi1(2,n1,md2/nproc,nd3/nproc,nproc),stat=i_stat)
     call memocc(i_stat,product(shape(zmpi1))*kind(zmpi1),'zmpi1','p_poissonsolver')
  end if


  !calculating the FFT work arrays (beware on the HalFFT in n3 dimension)
  call ctrig(n3,btrig3,after3,before3,now3,1,ic3)
  call ctrig(n1,btrig1,after1,before1,now1,1,ic1)
  call ctrig(n2,btrig2,after2,before2,now2,1,ic2)
  do  j=1,n1
     ftrig1(1,j)= btrig1(1,j)
     ftrig1(2,j)=-btrig1(2,j)
  enddo
  do  j=1,n2
     ftrig2(1,j)= btrig2(1,j)
     ftrig2(2,j)=-btrig2(2,j)
  enddo
  do  j=1,n3
     ftrig3(1,j)= btrig3(1,j)
     ftrig3(2,j)=-btrig3(2,j)
  enddo

  ! transform along z axis
  lot=ncache/(4*n3)
  if (lot.lt.1) then  
     write(6,*) & 
          'convolxc_off:ncache has to be enlarged to be able to hold at' // &  
          'least one 1-d FFT of this size even though this will' // & 
          'reduce the performance for shorter transform lengths'
     stop
  endif
  do j2=1,md2/nproc
     !this condition ensures that we manage only the interesting part for the FFT
     if (iproc*(md2/nproc)+j2.le.n2) then
        do i1=1,n1,lot
           ma=i1
           mb=min(i1+(lot-1),n1)
           nfft=mb-ma+1
           !inserting real data into complex array of half lenght
           call P_fill_upcorn(md1,md3,lot,nfft,n3,zf(i1,1,j2),zw(1,1,1))

           !performing FFT
           !input: I1,I3,J2,(Jp2)
           inzee=1
           do i=1,ic3
              call fftstp(lot,nfft,n3,lot,n3,zw(1,1,inzee),zw(1,1,3-inzee), &
                   btrig3,after3(i),now3(i),before3(i),1)
              inzee=3-inzee
           enddo

           !output: I1,i3,J2,(Jp2)
           !exchanging components
           !input: I1,i3,J2,(Jp2)
           call scramble_P(i1,j2,lot,nfft,n1,n3,md2,nproc,nd3,zw(1,1,inzee),zmpi2)
           !output: I1,J2,i3,(Jp2)
        end do
     end if
  end do

  !Interprocessor data transposition
  !input: I1,J2,j3,jp3,(Jp2)
  if (nproc.gt.1) then
     call timing(iproc,'PSolv_comput  ','OF')

     call timing(iproc,'PSolv_commun  ','ON')

     !communication scheduling
     call MPI_ALLTOALL(zmpi2,2*n1*(md2/nproc)*(nd3/nproc), &
          MPI_double_precision, &
          zmpi1,2*n1*(md2/nproc)*(nd3/nproc), &
          MPI_double_precision,MPI_COMM_WORLD,ierr)

     call timing(iproc,'PSolv_commun  ','OF')

     call timing(iproc,'PSolv_comput  ','ON')
  endif
  !output: I1,J2,j3,Jp2,(jp3)

  !now each process perform complete convolution of its planes
  do j3=1,nd3/nproc
     !this condition ensures that we manage only the interesting part for the FFT
     if (iproc*(nd3/nproc)+j3.le.n3/2+1) then
      Jp2stb=1
      J2stb=1
      Jp2stf=1
      J2stf=1
        
        ! transform along x axis
        lot=ncache/(4*n1)
        if (lot.lt.1) then  
           write(6,*) & 
                'convolxc_off:ncache has to be enlarged to be able to hold at' // &  
                'least one 1-d FFT of this size even though this will' // & 
                'reduce the performance for shorter transform lengths'
           stop
        endif
        
        do j=1,n2,lot
           ma=j
           mb=min(j+(lot-1),n2)
           nfft=mb-ma+1

           !reverse index ordering, leaving the planes to be transformed at the end
           !input: I1,J2,j3,Jp2,(jp3)
           if (nproc.eq.1) then
              call P_mpiswitch_upcorn(j3,nfft,Jp2stb,J2stb,lot,n1,md2,nd3,nproc,zmpi2,zw(1,1,1))
           else
              call P_mpiswitch_upcorn(j3,nfft,Jp2stb,J2stb,lot,n1,md2,nd3,nproc,zmpi1,zw(1,1,1))
           endif
           !output: J2,Jp2,I1,j3,(jp3)

           !performing FFT
           !input: I2,I1,j3,(jp3)
           inzee=1
           do i=1,ic1-1
              call fftstp(lot,nfft,n1,lot,n1,zw(1,1,inzee),zw(1,1,3-inzee), &
                   btrig1,after1(i),now1(i),before1(i),1)
              inzee=3-inzee
           enddo

           !storing the last step into zt array
           i=ic1
           call fftstp(lot,nfft,n1,lzt,n1,zw(1,1,inzee),zt(1,j,1), & 
                btrig1,after1(i),now1(i),before1(i),1)           
           !output: I2,i1,j3,(jp3)
        end do

        !transform along y axis
        lot=ncache/(4*n2)
        if (lot.lt.1) then  
           write(6,*) & 
                'convolxc_off:ncache has to be enlarged to be able to hold at' // &  
                'least one 1-d FFT of this size even though this will' // & 
                'reduce the performance for shorter transform lengths'
           stop
        endif

        do j=1,n1,lot
           ma=j
           mb=min(j+(lot-1),n1)
           nfft=mb-ma+1

           !reverse ordering 
           !input: I2,i1,j3,(jp3)
           call P_switch_upcorn(nfft,n2,lot,n1,lzt,zt(1,1,j),zw(1,1,1))
           !output: i1,I2,j3,(jp3)
           
           !performing FFT
           !input: i1,I2,j3,(jp3)
           inzee=1
           do i=1,ic2
              call fftstp(lot,nfft,n2,lot,n2,zw(1,1,inzee),zw(1,1,3-inzee), &
                   btrig2,after2(i),now2(i),before2(i),1)
              inzee=3-inzee
           enddo
           !output: i1,i2,j3,(jp3)


           !Multiply with kernel in fourier space
           i3=iproc*(nd3/nproc)+j3
           call P_multkernel(n1,n2,n3,lot,nfft,j,i3,zw(1,1,inzee),hx,hy,hz)
           
           !TRANSFORM BACK IN REAL SPACE
           
           !transform along y axis
           !input: i1,i2,j3,(jp3)
           do i=1,ic2
              call fftstp(lot,nfft,n2,lot,n2,zw(1,1,inzee),zw(1,1,3-inzee), &
                   ftrig2,after2(i),now2(i),before2(i),-1)
              inzee=3-inzee
           enddo

           !reverse ordering
           !input: i1,I2,j3,(jp3)
           call P_unswitch_downcorn(nfft,n2,lot,n1,lzt,zw(1,1,inzee),zt(1,1,j))
           !output: I2,i1,j3,(jp3)
        end do
        
        !transform along x axis
        !input: I2,i1,j3,(jp3)
        lot=ncache/(4*n1)
        do j=1,n2,lot
           ma=j
           mb=min(j+(lot-1),n2)
           nfft=mb-ma+1

           !performing FFT
           i=1
           call fftstp(lzt,nfft,n1,lot,n1,zt(1,j,1),zw(1,1,1), &
                ftrig1,after1(i),now1(i),before1(i),-1)
           
           inzee=1
           do i=2,ic1
              call fftstp(lot,nfft,n1,lot,n1,zw(1,1,inzee),zw(1,1,3-inzee), &
                   ftrig1,after1(i),now1(i),before1(i),-1)
              inzee=3-inzee
           enddo
           !output: I2,I1,j3,(jp3)

           !reverse ordering
           !input: J2,Jp2,I1,j3,(jp3)
           if (nproc.eq.1) then
              call P_unmpiswitch_downcorn(j3,nfft,Jp2stf,J2stf,lot,n1,md2,nd3,nproc,zw(1,1,inzee),zmpi2)
           else
              call P_unmpiswitch_downcorn(j3,nfft,Jp2stf,J2stf,lot,n1,md2,nd3,nproc,zw(1,1,inzee),zmpi1)
           endif
           ! output: I1,J2,j3,Jp2,(jp3)
        end do
     endif
  end do

  !Interprocessor data transposition
  !input: I1,J2,j3,Jp2,(jp3)
  if (nproc.gt.1) then
     call timing(iproc,'PSolv_comput  ','OF')

     call timing(iproc,'PSolv_commun  ','ON')
     !communication scheduling
     call MPI_ALLTOALL(zmpi1,2*n1*(md2/nproc)*(nd3/nproc), &
          MPI_double_precision, &
          zmpi2,2*n1*(md2/nproc)*(nd3/nproc), &
          MPI_double_precision,MPI_COMM_WORLD,ierr)
     call timing(iproc,'PSolv_commun  ','OF')

     call timing(iproc,'PSolv_comput  ','ON')

  endif
     !output: I1,J2,j3,jp3,(Jp2)
  !transform along z axis
  !input: I1,J2,i3,(Jp2)
  lot=ncache/(4*n3)
  do j2=1,md2/nproc
     !this condition ensures that we manage only the interesting part for the FFT
     if (iproc*(md2/nproc)+j2.le.n2) then
        do i1=1,n1,lot
           ma=i1
           mb=min(i1+(lot-1),n1)
           nfft=mb-ma+1

           !reverse ordering
           !input: I1,J2,i3,(Jp2)
           call unscramble_P(i1,j2,lot,nfft,n1,n3,md2,nproc,nd3,zmpi2,zw(1,1,1))
           !output: I1,i3,J2,(Jp2)

           !performing FFT
           !input: I1,i3,J2,(Jp2)           
           inzee=1
           do i=1,ic3
              call fftstp(lot,nfft,n3,lot,n3,zw(1,1,inzee),zw(1,1,3-inzee), &
                   ftrig3,after3(i),now3(i),before3(i),-1)
              inzee=3-inzee
           enddo
           !output: I1,I3,J2,(Jp2)

           !rebuild the output array
           call P_unfill_downcorn(md1,md3,lot,nfft,n3,zw(1,1,inzee),zf(i1,1,j2),scal)

        end do
     endif
  end do

  !De-allocations  
  i_all=-product(shape(btrig1))*kind(btrig1)
  deallocate(btrig1,stat=i_stat)
  call memocc(i_stat,i_all,'btrig1','p_poissonsolver')
  i_all=-product(shape(ftrig1))*kind(ftrig1)
  deallocate(ftrig1,stat=i_stat)
  call memocc(i_stat,i_all,'ftrig1','p_poissonsolver')
  i_all=-product(shape(after1))*kind(after1)
  deallocate(after1,stat=i_stat)
  call memocc(i_stat,i_all,'after1','p_poissonsolver')
  i_all=-product(shape(now1))*kind(now1)
  deallocate(now1,stat=i_stat)
  call memocc(i_stat,i_all,'now1','p_poissonsolver')
  i_all=-product(shape(before1))*kind(before1)
  deallocate(before1,stat=i_stat)
  call memocc(i_stat,i_all,'before1','p_poissonsolver')
  i_all=-product(shape(btrig2))*kind(btrig2)
  deallocate(btrig2,stat=i_stat)
  call memocc(i_stat,i_all,'btrig2','p_poissonsolver')
  i_all=-product(shape(ftrig2))*kind(ftrig2)
  deallocate(ftrig2,stat=i_stat)
  call memocc(i_stat,i_all,'ftrig2','p_poissonsolver')
  i_all=-product(shape(after2))*kind(after2)
  deallocate(after2,stat=i_stat)
  call memocc(i_stat,i_all,'after2','p_poissonsolver')
  i_all=-product(shape(now2))*kind(now2)
  deallocate(now2,stat=i_stat)
  call memocc(i_stat,i_all,'now2','p_poissonsolver')
  i_all=-product(shape(before2))*kind(before2)
  deallocate(before2,stat=i_stat)
  call memocc(i_stat,i_all,'before2','p_poissonsolver')
  i_all=-product(shape(btrig3))*kind(btrig3)
  deallocate(btrig3,stat=i_stat)
  call memocc(i_stat,i_all,'btrig3','p_poissonsolver')
  i_all=-product(shape(ftrig3))*kind(ftrig3)
  deallocate(ftrig3,stat=i_stat)
  call memocc(i_stat,i_all,'ftrig3','p_poissonsolver')
  i_all=-product(shape(after3))*kind(after3)
  deallocate(after3,stat=i_stat)
  call memocc(i_stat,i_all,'after3','p_poissonsolver')
  i_all=-product(shape(now3))*kind(now3)
  deallocate(now3,stat=i_stat)
  call memocc(i_stat,i_all,'now3','p_poissonsolver')
  i_all=-product(shape(before3))*kind(before3)
  deallocate(before3,stat=i_stat)
  call memocc(i_stat,i_all,'before3','p_poissonsolver')
  i_all=-product(shape(zmpi2))*kind(zmpi2)
  deallocate(zmpi2,stat=i_stat)
  call memocc(i_stat,i_all,'zmpi2','p_poissonsolver')
  i_all=-product(shape(zw))*kind(zw)
  deallocate(zw,stat=i_stat)
  call memocc(i_stat,i_all,'zw','p_poissonsolver')
  i_all=-product(shape(zt))*kind(zt)
  deallocate(zt,stat=i_stat)
  call memocc(i_stat,i_all,'zt','p_poissonsolver')
  if (nproc.gt.1) then
     i_all=-product(shape(zmpi1))*kind(zmpi1)
     deallocate(zmpi1,stat=i_stat)
     call memocc(i_stat,i_all,'zmpi1','p_poissonsolver')
  end if
  call timing(iproc,'PSolv_comput  ','OF')
end subroutine P_PoissonSolver


subroutine P_mpiswitch_upcorn(j3,nfft,Jp2stb,J2stb,lot,n1,md2,nd3,nproc,zmpi1,zw)
  implicit none
  integer, intent(in) :: j3,nfft,lot,n1,md2,nd3,nproc
  integer, intent(inout) :: Jp2stb,J2stb
  real(kind=8), dimension(2,n1,md2/nproc,nd3/nproc,nproc), intent(in) :: zmpi1
  real(kind=8), dimension(2,lot,n1), intent(out) :: zw
  !local variables
  integer :: I1,J2,Jp2,mfft

  mfft=0
  do Jp2=Jp2stb,nproc
     do J2=J2stb,md2/nproc
      mfft=mfft+1
      if (mfft.gt.nfft) then
           Jp2stb=Jp2
           J2stb=J2
           return
      end if
        do I1=1,n1
           zw(1,mfft,I1)=zmpi1(1,I1,J2,j3,Jp2)
           zw(2,mfft,I1)=zmpi1(2,I1,J2,j3,Jp2)
        end do
     end do
     J2stb=1
  end do

end subroutine P_mpiswitch_upcorn


subroutine P_switch_upcorn(nfft,n2,lot,n1,lzt,zt,zw)
  implicit none
  integer, intent(in) :: nfft,n2,lot,n1,lzt
  real(kind=8), dimension(2,lzt,n1), intent(in) :: zt
  real(kind=8), dimension(2,lot,n2), intent(out) :: zw
  !local variables
  integer :: i,j

  do j=1,nfft
     do i=1,n2
      zw(1,j,i)=zt(1,i,j)
      zw(2,j,i)=zt(2,i,j)
     end do
  end do

end subroutine P_switch_upcorn


subroutine P_unswitch_downcorn(nfft,n2,lot,n1,lzt,zw,zt)
  implicit none
  integer, intent(in) :: nfft,n2,lot,n1,lzt
  real(kind=8), dimension(2,lot,n2), intent(in) :: zw
  real(kind=8), dimension(2,lzt,n1), intent(out) :: zt
  !local variables
  integer :: i,j

  do j=1,nfft
     do i=1,n2
      zt(1,i,j)=zw(1,j,i)
      zt(2,i,j)=zw(2,j,i)
     end do
  end do

end subroutine P_unswitch_downcorn


subroutine P_unmpiswitch_downcorn(j3,nfft,Jp2stf,J2stf,lot,n1,md2,nd3,nproc,zw,zmpi1)
  implicit none
  integer, intent(in) :: j3,nfft,lot,n1,md2,nd3,nproc
  integer, intent(inout) :: Jp2stf,J2stf
  real(kind=8), dimension(2,lot,n1), intent(in) :: zw
  real(kind=8), dimension(2,n1,md2/nproc,nd3/nproc,nproc), intent(out) :: zmpi1
  !local variables
  integer :: I1,J2,Jp2,mfft

  mfft=0
  do Jp2=Jp2stf,nproc
     do J2=J2stf,md2/nproc
      mfft=mfft+1
      if (mfft.gt.nfft) then
           Jp2stf=Jp2
           J2stf=J2
           return
      end if
        do I1=1,n1
           zmpi1(1,I1,J2,j3,Jp2)=zw(1,mfft,I1)
           zmpi1(2,I1,J2,j3,Jp2)=zw(2,mfft,I1)
        end do
     end do
     J2stf=1
  end do
end subroutine P_unmpiswitch_downcorn


!!****f* BigDFT/P_fill_downcorn
!! NAME
!!   P_fill_downcorn
!!
!! FUNCTION
!!     (Based on suitable modifications of S.Goedecker routines)
!!     Restore data into output array
!!
!! SYNOPSIS
!!     zf:          Original distributed density as well as
!!                  Distributed solution of the poisson equation (inout)
!!     zw:          FFT work array
!!     n3:          (twice the) dimension of the last FFTtransform.
!!     md1,md3:     Dimensions of the undistributed part of the real grid
!!     nfft:        number of planes
!!     scal:        Needed to achieve unitarity and correct dimensions
!!
!! WARNING
!!     Assuming that high frequencies are in the corners 
!!     and that n3 is multiple of 4   
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
subroutine P_unfill_downcorn(md1,md3,lot,nfft,n3,zw,zf,scal)
  implicit none
  !Arguments
  integer, intent(in) :: md1,md3,lot,nfft,n3
  real(kind=8), dimension(2,lot,n3), intent(in) :: zw
  real(kind=8), dimension(md1,md3),intent(inout) :: zf
  real(kind=8), intent(in) :: scal
  !Local variables
  integer :: i3,i1
  real(kind=8) :: pot1

  do i3=1,n3
     do i1=1,nfft
        pot1 = scal*zw(1,i1,i3)
      zf(i1,i3)= pot1 
     end do
  end do

end subroutine P_unfill_downcorn
!!**

subroutine P_fill_upcorn(md1,md3,lot,nfft,n3,zf,zw)
  implicit none
  integer, intent(in) :: md1,md3,lot,nfft,n3
  real(kind=8), dimension(md1,md3), intent(in) :: zf
  real(kind=8), dimension(2,lot,n3), intent(out) :: zw
  !local variables
  integer :: i1,i3

  do i3=1,n3
     do i1=1,nfft
      zw(1,i1,i3)=zf(i1,i3)
      zw(2,i1,i3)=0.d0
     end do
  end do

end subroutine P_fill_upcorn


!!****f* BigDFT/scramble_P
!! NAME
!!   scramble_P
!!
!! FUNCTION
!!     (Based on suitable modifications of S.Goedecker routines)
!!     Assign the correct planes to the work array zmpi2
!!     in order to prepare for interprocessor data transposition.
!!
!! SYNOPSIS
!!     zmpi2:          Work array for multiprocessor manipulation (output)
!!     zw:             Work array (input)
!!     n1,n3:          logical dimension of the FFT transform, reference for work arrays
!!     md2,nd3:        Dimensions of real grid and of the kernel, respectively
!!     i1,j2,lot,nfft: Starting points of the plane and number of remaining lines
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
subroutine scramble_P(i1,j2,lot,nfft,n1,n3,md2,nproc,nd3,zw,zmpi2)
  implicit none
  !Arguments
  integer, intent(in) :: i1,j2,lot,nfft,n1,n3,md2,nproc,nd3
  real(kind=8), dimension(2,lot,n3), intent(in) :: zw
  real(kind=8), dimension(2,n1,md2/nproc,nd3), intent(out) :: zmpi2
  !Local variables
  integer :: i3,i

  do i3=1,n3/2+1
     do i=0,nfft-1
        zmpi2(1,i1+i,j2,i3)=zw(1,i+1,i3)
        zmpi2(2,i1+i,j2,i3)=zw(2,i+1,i3)
     end do
  end do

end subroutine scramble_P
!!**
 
!!****f* BigDFT/unscramble_P
!! NAME
!!   unscramble_P
!!
!! FUNCTION
!!     (Based on suitable modifications of S.Goedecker routines)
!!     Insert the correct planes of the work array zmpi2
!!     in order to prepare for backward FFT transform
!!
!! SYNOPSIS
!!     zmpi2:          Work array for multiprocessor manipulation (input)
!!     zw:             Work array (output)
!!     n1,n3:          logical dimension of the FFT transform, reference for work arrays
!!     md2,nd3:        Dimensions of real grid and of the kernel, respectively
!!     i1,j2,lot,nfft: Starting points of the plane and number of remaining lines
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
subroutine unscramble_P(i1,j2,lot,nfft,n1,n3,md2,nproc,nd3,zmpi2,zw)
  implicit none
  !Arguments
  integer, intent(in) :: i1,j2,lot,nfft,n1,n3,md2,nproc,nd3
  real(kind=8), dimension(2,lot,n3), intent(out) :: zw
  real(kind=8), dimension(2,n1,md2/nproc,nd3), intent(in) :: zmpi2
  !Local variables
  integer :: i3,i,j3

  i3=1
  do i=0,nfft-1
     zw(1,i+1,i3)=zmpi2(1,i1+i,j2,i3)
     zw(2,i+1,i3)=zmpi2(2,i1+i,j2,i3)
  end do

  do i3=2,n3/2+1
     j3=n3+2-i3
     do i=0,nfft-1
        zw(1,i+1,j3)= zmpi2(1,i1+i,j2,i3)
        zw(2,i+1,j3)=-zmpi2(2,i1+i,j2,i3)
        zw(1,i+1,i3)= zmpi2(1,i1+i,j2,i3)
        zw(2,i+1,i3)= zmpi2(2,i1+i,j2,i3)
     end do
  end do

end subroutine unscramble_P
!!**


!!****f* BigDFT/P_multkernel
!! NAME
!!   P_multkernel
!!
!! FUNCTION
!!     (Based on suitable modifications of S.Goedecker routines)
!!     Multiply with the kernel taking into account its symmetry
!!     Conceived to be used into convolution loops
!!
!! SYNOPSIS
!!     zw:         Work array (input/output)
!!     n1,n2:      logical dimension of the FFT transform, reference for zw
!!     nd1,nd2:    Dimensions of POT
!!     jS,j3,nfft: starting point of the plane and number of remaining lines
!!
!! RESTRICTIONS on USAGE
!!     Copyright (C) Stefan Goedecker, Cornell University, Ithaca, USA, 1994
!!     Copyright (C) Stefan Goedecker, MPI Stuttgart, Germany, 1999
!!     Copyright (C) 2002 Stefan Goedecker, CEA Grenoble
!!     This file is distributed under the terms of the
!!      GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .
!!
!! AUTHORS
!!    S. Goedecker, L. Genovese
!!
!! CREATION DATE
!!     February 2006
!!
!! SOURCE
!!
subroutine P_multkernel(n1,n2,n3,lot,nfft,jS,i3,zw,hx,hy,hz)
  implicit none
  !Argments
  integer, intent(in) :: n1,n2,n3,lot,nfft,jS,i3
  real(kind=8), intent(in) :: hx,hy,hz
  real(kind=8), dimension(2,lot,n2), intent(inout) :: zw
  !Local variables
  integer :: i1,j1,i2,j2,j3
  real(kind=8) :: ker,pi,fourpi2,mu3,p1,p2

  pi=4.d0*datan(1.d0)
  fourpi2=4.d0*pi**2
  j3=i3!n3/2+1-abs(n3/2+2-i3)
  mu3=real(j3-1,kind=8)/real(n3,kind=8)
  mu3=(mu3/hy)**2 !beware of the exchanged dimension
  !Body
  !generic case
  do i2=1,n2
     do i1=1,nfft
        j1=i1+jS-1
        j1=j1-(j1/(n1/2+2))*n1 !n1/2+1-abs(n1/2+2-jS-i1)
        j2=i2-(i2/(n2/2+2))*n2 !n2/2+1-abs(n2/2+1-i2)
        p1=real(j1-1,kind=8)/real(n1,kind=8)
        p2=real(j2-1,kind=8)/real(n2,kind=8)
        ker=-fourpi2*((p1/hx)**2+(p2/hz)**2+mu3)!beware of the exchanged dimension
        if (ker/=0.d0) ker=1.d0/ker
        zw(1,i1,i2)=zw(1,i1,i2)*ker
        zw(2,i1,i2)=zw(2,i1,i2)*ker
     end do
  end do

end subroutine P_multkernel
!!**


!!****f* BigDFT/multkernel
!! NAME
!!   multkernel
!!
!! FUNCTION
!!     (Based on suitable modifications of S.Goedecker routines)
!!     Multiply with the kernel taking into account its symmetry
!!     Conceived to be used into convolution loops
!!
!! SYNOPSIS
!!     pot:      Kernel, symmetric and real, half the length
!!     zw:       Work array (input/output)
!!     n1,n2:    logical dimension of the FFT transform, reference for zw
!!     nd1,nd2:  Dimensions of POT
!!     jS, nfft: starting point of the plane and number of remaining lines
!!
!! RESTRICTIONS on USAGE
!!     Copyright (C) Stefan Goedecker, Cornell University, Ithaca, USA, 1994
!!     Copyright (C) Stefan Goedecker, MPI Stuttgart, Germany, 1999
!!     Copyright (C) 2002 Stefan Goedecker, CEA Grenoble
!!     This file is distributed under the terms of the
!!      GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .
!!
!! AUTHORS
!!    S. Goedecker, L. Genovese
!!
!! CREATION DATE
!!     February 2006
!!
!! SOURCE
!!
subroutine multkernel(nd1,nd2,n1,n2,lot,nfft,jS,pot,zw)
  implicit none
  !Argments
  integer, intent(in) :: nd1,nd2,n1,n2,lot,nfft,jS
  real(kind=8), dimension(nd1,nd2), intent(in) :: pot
  real(kind=8), dimension(2,lot,n2), intent(inout) :: zw
  !Local variables
  integer :: j,j1,i2,j2

  !Body
  
  !case i2=1
  do j=1,nfft
     j1=j+jS-1
     !isign=(j1/(n1/2+2))
     !j1=(1-2*isign)*j1+isign*(n1+2) !n1/2+1-abs(n1/2+2-jS-i1)
     j1=j1+(j1/(n1/2+2))*(n1+2-2*j1)
!!$     j1=n1/2+1-abs(n1/2+2-jS-j)!this stands for j1=min(jS-1+j,n1+3-jS-j)
     zw(1,j,1)=zw(1,j,1)*pot(j1,1)
     zw(2,j,1)=zw(2,j,1)*pot(j1,1)
  end do

  !generic case
  do i2=2,n2/2
     do j=1,nfft
        j1=j+jS-1
        j1=j1+(j1/(n1/2+2))*(n1+2-2*j1)
!!$        j1=n1/2+1-abs(n1/2+2-jS-j)
        j2=n2+2-i2
        zw(1,j,i2)=zw(1,j,i2)*pot(j1,i2)
        zw(2,j,i2)=zw(2,j,i2)*pot(j1,i2)
        zw(1,j,j2)=zw(1,j,j2)*pot(j1,i2)
        zw(2,j,j2)=zw(2,j,j2)*pot(j1,i2)
     end do
  end do
  
  !case i2=n2/2+1
  do j=1,nfft
     j1=j+jS-1
     j1=j1+(j1/(n1/2+2))*(n1+2-2*j1)
!!$     j1=n1/2+1-abs(n1/2+2-jS-j)
     j2=n2/2+1
     zw(1,j,j2)=zw(1,j,j2)*pot(j1,j2)
     zw(2,j,j2)=zw(2,j,j2)*pot(j1,j2)
  end do

end subroutine multkernel
!!**


!!!HERE POT MUST BE THE KERNEL (BEWARE THE HALF DIMENSION)


!!****f* BigDFT/S_PoissonSolver
!! NAME
!!   PoissonSolver
!!
!! FUNCTION
!!     (Based on suitable modifications of S.Goedecker routines)
!!     Applies the local FFT space Kernel to the density in Real space.
!!     Does NOT calculate the LDA exchange-correlation terms
!!
!! SYNOPSIS
!!     zf:          Density (input/output)
!!                  ZF(i1,i3,i2)
!!                  i1=1,md1 , i2=1,md2/nproc , i3=1,md3 
!!     pot:         Kernel, only the distributed part (REAL)
!!                  POT(i1,i2,i3)
!!                  i1=1,nd1 , i2=1,nd2 , i3=1,nd3/nproc 
!!     nproc:       number of processors used as returned by MPI_COMM_SIZE
!!     iproc:       [0:nproc-1] number of processor as returned by MPI_COMM_RANK
!!     n1,n2,n3/2:  logical dimension of the transform. As transform lengths 
!!                  most products of the prime factors 2,3,5 are allowed.
!!                  The detailed table with allowed transform lengths can 
!!                  be found in subroutine CTRIG
!!     md1,md2,md3: Dimension of ZF
!!     nd1,nd2,nd3/nproc: Dimension of POT
!!     scal:        factor of renormalization of the FFT in order to acheve unitarity 
!!                  and the correct dimension
!!     hx,hy,hz:    grid spacing, used for integrating eharthree
!!     ehartree:    hartree energy of the potential
!!
!! RESTRICTIONS on USAGE
!!     Copyright (C) Stefan Goedecker, Cornell University, Ithaca, USA, 1994
!!     Copyright (C) Stefan Goedecker, MPI Stuttgart, Germany, 1999
!!     Copyright (C) 2002 Stefan Goedecker, CEA Grenoble
!!     This file is distributed under the terms of the
!!      GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .
!!
!! AUTHORS
!!    S. Goedecker, L. Genovese
!!
!! CREATION DATE
!!    October 2006
!!
!! SOURCE
!!
subroutine S_PoissonSolver(n1,n2,n3,nd1,nd2,nd3,md1,md2,md3,nproc,iproc,pot,zf&
             ,scal,hx,hy,hz)!,ehartree)
  use mpi
  implicit none
  include 'perfdata.inc'
  !Arguments
  integer, intent(in) :: n1,n2,n3,nd1,nd2,nd3,md1,md2,md3,nproc,iproc
  real(kind=8), intent(in) :: scal,hx,hy,hz
  !real(kind=8), intent(out) :: ehartree
  real(kind=8), dimension(nd1,nd2,nd3/nproc), intent(in) :: pot
  real(kind=8), dimension(md1,md3,md2/nproc), intent(inout) :: zf
  !Local variables
  !Maximum number of points for FFT (should be same number in fft3d routine)
  integer, parameter :: nfft_max=24000
  integer :: ncache,lzt,lot,ma,mb,nfft,ic1,ic2,ic3,Jp2stb,J2stb,Jp2stf,J2stf
  integer :: j2,j3,i1,i3,i,j,inzee,ierr,i_all,i_stat
  real(kind=8) :: twopion!,ehartreetmp
  !work arrays for transpositions
  real(kind=8), dimension(:,:,:), allocatable :: zt
  !work arrays for MPI
  real(kind=8), dimension(:,:,:,:,:), allocatable :: zmpi1
  real(kind=8), dimension(:,:,:,:), allocatable :: zmpi2
  !cache work array
  real(kind=8), dimension(:,:,:), allocatable :: zw
  !FFT work arrays
  real(kind=8), dimension(:,:), allocatable :: btrig1,btrig2,btrig3, &
       ftrig1,ftrig2,ftrig3,cosinarr
  integer, dimension(:), allocatable :: after1,now1,before1, & 
       after2,now2,before2,after3,now3,before3
  
  !Body

  call timing(iproc,'PSolv_comput  ','ON')
  ! check input
  !if (mod(n1,2).ne.0) stop 'Parallel convolution:ERROR:n1' !this can be avoided
  !if (mod(n2,2).ne.0) stop 'Parallel convolution:ERROR:n2' !this can be avoided
  if (mod(n3,2).ne.0) stop 'Parallel convolution:ERROR:n3'
  if (nd1.lt.n1/2+1) stop 'Parallel convolution:ERROR:nd1' 
  if (nd2.lt.n2/2+1) stop 'Parallel convolution:ERROR:nd2' 
  if (nd3.lt.n3/2+1) stop 'Parallel convolution:ERROR:nd3'
  if (md1.lt.n1) stop 'Parallel convolution:ERROR:md1'
  if (md2.lt.n2) stop 'Parallel convolution:ERROR:md2'
  if (md3.lt.n3/2) stop 'Parallel convolution:ERROR:md3'
  if (mod(nd3,nproc).ne.0) stop 'Parallel convolution:ERROR:nd3'
  if (mod(md2,nproc).ne.0) stop 'Parallel convolution:ERROR:md2'
  
  !defining work arrays dimensions
  ncache=ncache_optimal
  if (ncache <= max(n1,n2,n3/2)*4) ncache=max(n1,n2,n3/2)*4

  if (timing_flag == 1 .and. iproc ==0) print *,'parallel ncache=',ncache

  lzt=n2
  if (mod(n2,2).eq.0) lzt=lzt+1
  if (mod(n2,4).eq.0) lzt=lzt+1 !maybe this is useless
  
  !Allocations
  allocate(btrig1(2,nfft_max),stat=i_stat)
  call memocc(i_stat,product(shape(btrig1))*kind(btrig1),'btrig1','s_poissonsolver')
  allocate(ftrig1(2,nfft_max),stat=i_stat)
  call memocc(i_stat,product(shape(ftrig1))*kind(ftrig1),'ftrig1','s_poissonsolver')
  allocate(after1(7),stat=i_stat)
  call memocc(i_stat,product(shape(after1))*kind(after1),'after1','s_poissonsolver')
  allocate(now1(7),stat=i_stat)
  call memocc(i_stat,product(shape(now1))*kind(now1),'now1','s_poissonsolver')
  allocate(before1(7),stat=i_stat)
  call memocc(i_stat,product(shape(before1))*kind(before1),'before1','s_poissonsolver')
  allocate(btrig2(2,nfft_max),stat=i_stat)
  call memocc(i_stat,product(shape(btrig2))*kind(btrig2),'btrig2','s_poissonsolver')
  allocate(ftrig2(2,nfft_max),stat=i_stat)
  call memocc(i_stat,product(shape(ftrig2))*kind(ftrig2),'ftrig2','s_poissonsolver')
  allocate(after2(7),stat=i_stat)
  call memocc(i_stat,product(shape(after2))*kind(after2),'after2','s_poissonsolver')
  allocate(now2(7),stat=i_stat)
  call memocc(i_stat,product(shape(now2))*kind(now2),'now2','s_poissonsolver')
  allocate(before2(7),stat=i_stat)
  call memocc(i_stat,product(shape(before2))*kind(before2),'before2','s_poissonsolver')
  allocate(btrig3(2,nfft_max),stat=i_stat)
  call memocc(i_stat,product(shape(btrig3))*kind(btrig3),'btrig3','s_poissonsolver')
  allocate(ftrig3(2,nfft_max),stat=i_stat)
  call memocc(i_stat,product(shape(ftrig3))*kind(ftrig3),'ftrig3','s_poissonsolver')
  allocate(after3(7),stat=i_stat)
  call memocc(i_stat,product(shape(after3))*kind(after3),'after3','s_poissonsolver')
  allocate(now3(7),stat=i_stat)
  call memocc(i_stat,product(shape(now3))*kind(now3),'now3','s_poissonsolver')
  allocate(before3(7),stat=i_stat)
  call memocc(i_stat,product(shape(before3))*kind(before3),'before3','s_poissonsolver')
  allocate(zw(2,ncache/4,2),stat=i_stat)
  call memocc(i_stat,product(shape(zw))*kind(zw),'zw','s_poissonsolver')
  allocate(zt(2,lzt,n1),stat=i_stat)
  call memocc(i_stat,product(shape(zt))*kind(zt),'zt','s_poissonsolver')
  allocate(zmpi2(2,n1,md2/nproc,nd3),stat=i_stat)
  call memocc(i_stat,product(shape(zmpi2))*kind(zmpi2),'zmpi2','s_poissonsolver')
  allocate(cosinarr(2,n3/2),stat=i_stat)
  call memocc(i_stat,product(shape(cosinarr))*kind(cosinarr),'cosinarr','s_poissonsolver')
  if (nproc.gt.1) then 
     allocate(zmpi1(2,n1,md2/nproc,nd3/nproc,nproc),stat=i_stat)
     call memocc(i_stat,product(shape(zmpi1))*kind(zmpi1),'zmpi1','s_poissonsolver')
  end if


  !calculating the FFT work arrays (beware on the HalFFT in n3 dimension)
  call ctrig(n3/2,btrig3,after3,before3,now3,1,ic3)
  call ctrig(n1,btrig1,after1,before1,now1,1,ic1)
  call ctrig(n2,btrig2,after2,before2,now2,1,ic2)
  do  j=1,n1
     ftrig1(1,j)= btrig1(1,j)
     ftrig1(2,j)=-btrig1(2,j)
  enddo
  do  j=1,n2
     ftrig2(1,j)= btrig2(1,j)
     ftrig2(2,j)=-btrig2(2,j)
  enddo
  do  j=1,n3/2
     ftrig3(1,j)= btrig3(1,j)
     ftrig3(2,j)=-btrig3(2,j)
  enddo

  !Calculating array of phases for HalFFT decoding
  twopion=8.d0*datan(1.d0)/real(n3,kind=8)
  do i3=1,n3/2
     cosinarr(1,i3)= dcos(twopion*real(i3-1,kind=8))
     cosinarr(2,i3)=-dsin(twopion*real(i3-1,kind=8))
  end do

  !initializing integral
  !ehartree=0.d0

  ! transform along z axis
  lot=ncache/(2*n3)
  if (lot.lt.1) then  
     write(6,*) & 
          'convolxc_off:ncache has to be enlarged to be able to hold at' // &  
          'least one 1-d FFT of this size even though this will' // & 
          'reduce the performance for shorter transform lengths'
     stop
  endif
  
  do j2=1,md2/nproc
     !this condition ensures that we manage only the interesting part for the FFT
     if (iproc*(md2/nproc)+j2.le.n2) then
        do i1=1,n1,lot
           ma=i1
           mb=min(i1+(lot-1),n1)
           nfft=mb-ma+1

           !inserting real data into complex array of half lenght
           call halfill_upcorn(md1,md3,lot,nfft,n3,zf(i1,1,j2),zw(1,1,1))

           !performing FFT
           !input: I1,I3,J2,(Jp2)
           inzee=1
           do i=1,ic3
              call fftstp(lot,nfft,n3/2,lot,n3/2,zw(1,1,inzee),zw(1,1,3-inzee), &
                   btrig3,after3(i),now3(i),before3(i),1)
              inzee=3-inzee
           enddo
           !output: I1,i3,J2,(Jp2)

           !unpacking FFT in order to restore correct result, 
           !while exchanging components
           !input: I1,i3,J2,(Jp2)
           call scramble_unpack(i1,j2,lot,nfft,n1,n3,md2,nproc,nd3,zw(1,1,inzee),zmpi2,cosinarr)
           !output: I1,J2,i3,(Jp2)
        end do
     end if
  end do

  !Interprocessor data transposition
  !input: I1,J2,j3,jp3,(Jp2)
  if (nproc.gt.1) then
     call timing(iproc,'PSolv_comput  ','OF')

     call timing(iproc,'PSolv_commun  ','ON')
     !communication scheduling
     call MPI_ALLTOALL(zmpi2,2*n1*(md2/nproc)*(nd3/nproc), &
          MPI_double_precision, &
          zmpi1,2*n1*(md2/nproc)*(nd3/nproc), &
          MPI_double_precision,MPI_COMM_WORLD,ierr)

     call timing(iproc,'PSolv_commun  ','OF')

     call timing(iproc,'PSolv_comput  ','ON')
  endif
  !output: I1,J2,j3,Jp2,(jp3)

  !now each process perform complete convolution of its planes
  do j3=1,nd3/nproc
     !this condition ensures that we manage only the interesting part for the FFT
     if (iproc*(nd3/nproc)+j3.le.n3/2+1) then
      Jp2stb=1
      J2stb=1
      Jp2stf=1
      J2stf=1
        
        ! transform along x axis
        lot=ncache/(4*n1)
        if (lot.lt.1) then  
           write(6,*) & 
                'convolxc_off:ncache has to be enlarged to be able to hold at' // &  
                'least one 1-d FFT of this size even though this will' // & 
                'reduce the performance for shorter transform lengths'
           stop
        endif
        
        do j=1,n2,lot
           ma=j
           mb=min(j+(lot-1),n2)
           nfft=mb-ma+1

           !reverse index ordering, leaving the planes to be transformed at the end
           !input: I1,J2,j3,Jp2,(jp3)
           if (nproc.eq.1) then
              call S_mpiswitch_upcorn(j3,nfft,Jp2stb,J2stb,lot,n1,md2,nd3,nproc,zmpi2,zw(1,1,1))
           else
              call S_mpiswitch_upcorn(j3,nfft,Jp2stb,J2stb,lot,n1,md2,nd3,nproc,zmpi1,zw(1,1,1))
           endif
           !output: J2,Jp2,I1,j3,(jp3)

           !performing FFT
           !input: I2,I1,j3,(jp3)
           inzee=1
           do i=1,ic1-1
              call fftstp(lot,nfft,n1,lot,n1,zw(1,1,inzee),zw(1,1,3-inzee), &
                   btrig1,after1(i),now1(i),before1(i),1)
              inzee=3-inzee
           enddo

           !storing the last step into zt array
           i=ic1
           call fftstp(lot,nfft,n1,lzt,n1,zw(1,1,inzee),zt(1,j,1), & 
                btrig1,after1(i),now1(i),before1(i),1)           
           !output: I2,i1,j3,(jp3)
        end do

        !transform along y axis
        lot=ncache/(4*n2)
        if (lot.lt.1) then  
           write(6,*) & 
                'convolxc_off:ncache has to be enlarged to be able to hold at' // &  
                'least one 1-d FFT of this size even though this will' // & 
                'reduce the performance for shorter transform lengths'
           stop
        endif

        do j=1,n1,lot
           ma=j
           mb=min(j+(lot-1),n1)
           nfft=mb-ma+1

           !reverse ordering 
           !input: I2,i1,j3,(jp3)
           call S_switch_upcorn(nfft,n2,lot,n1,lzt,zt(1,1,j),zw(1,1,1))
           !output: i1,I2,j3,(jp3)
           
           !performing FFT
           !input: i1,I2,j3,(jp3)
           inzee=1
           do i=1,ic2
              call fftstp(lot,nfft,n2,lot,n2,zw(1,1,inzee),zw(1,1,3-inzee), &
                   btrig2,after2(i),now2(i),before2(i),1)
              inzee=3-inzee
           enddo
           !output: i1,i2,j3,(jp3)
           
           !Multiply with kernel in fourier space
           call multkernel(nd1,nd2,n1,n2,lot,nfft,j,pot(1,1,j3),zw(1,1,inzee))
           
           !TRANSFORM BACK IN REAL SPACE
           
           !transform along y axis
           !input: i1,i2,j3,(jp3)
           do i=1,ic2
              call fftstp(lot,nfft,n2,lot,n2,zw(1,1,inzee),zw(1,1,3-inzee), &
                   ftrig2,after2(i),now2(i),before2(i),-1)
              inzee=3-inzee
           enddo

           !reverse ordering
           !input: i1,I2,j3,(jp3)
           call S_unswitch_downcorn(nfft,n2,lot,n1,lzt,zw(1,1,inzee),zt(1,1,j))
           !output: I2,i1,j3,(jp3)
        end do
        
        !transform along x axis
        !input: I2,i1,j3,(jp3)
        lot=ncache/(4*n1)
        do j=1,n2,lot
           ma=j
           mb=min(j+(lot-1),n2)
           nfft=mb-ma+1

           !performing FFT
           i=1
           call fftstp(lzt,nfft,n1,lot,n1,zt(1,j,1),zw(1,1,1), &
                ftrig1,after1(i),now1(i),before1(i),-1)
           
           inzee=1
           do i=2,ic1
              call fftstp(lot,nfft,n1,lot,n1,zw(1,1,inzee),zw(1,1,3-inzee), &
                   ftrig1,after1(i),now1(i),before1(i),-1)
              inzee=3-inzee
           enddo
           !output: I2,I1,j3,(jp3)

           !reverse ordering
           !input: J2,Jp2,I1,j3,(jp3)
           if (nproc.eq.1) then
              call S_unmpiswitch_downcorn(j3,nfft,Jp2stf,J2stf,lot,n1,md2,nd3,nproc,zw(1,1,inzee),zmpi2)
           else
              call S_unmpiswitch_downcorn(j3,nfft,Jp2stf,J2stf,lot,n1,md2,nd3,nproc,zw(1,1,inzee),zmpi1)
           endif
           ! output: I1,J2,j3,Jp2,(jp3)
        end do
     endif
  end do


  !Interprocessor data transposition
  !input: I1,J2,j3,Jp2,(jp3)
  if (nproc.gt.1) then
     call timing(iproc,'PSolv_comput  ','OF')

     call timing(iproc,'PSolv_commun  ','ON')
     !communication scheduling
     call MPI_ALLTOALL(zmpi1,2*n1*(md2/nproc)*(nd3/nproc), &
          MPI_double_precision, &
          zmpi2,2*n1*(md2/nproc)*(nd3/nproc), &
          MPI_double_precision,MPI_COMM_WORLD,ierr)
     call timing(iproc,'PSolv_commun  ','OF')

     call timing(iproc,'PSolv_comput  ','ON')
  endif

 !output: I1,J2,j3,jp3,(Jp2)

  !transform along z axis
  !input: I1,J2,i3,(Jp2)
  lot=ncache/(2*n3)
  do j2=1,md2/nproc
     !this condition ensures that we manage only the interesting part for the FFT
     if (iproc*(md2/nproc)+j2.le.n2) then
        do i1=1,n1,lot
           ma=i1
           mb=min(i1+(lot-1),n1)
           nfft=mb-ma+1

           !reverse ordering and repack the FFT data in order to be backward HalFFT transformed
           !input: I1,J2,i3,(Jp2)
           call unscramble_pack(i1,j2,lot,nfft,n1,n3,md2,nproc,nd3,zmpi2,zw(1,1,1),cosinarr)
           !output: I1,i3,J2,(Jp2)

           !performing FFT
           !input: I1,i3,J2,(Jp2)           
           inzee=1
           do i=1,ic3
              call fftstp(lot,nfft,n3/2,lot,n3/2,zw(1,1,inzee),zw(1,1,3-inzee), &
                   ftrig3,after3(i),now3(i),before3(i),-1)
              inzee=3-inzee
           enddo
           !output: I1,I3,J2,(Jp2)

           !rebuild the output array
           call unfill_downcorn(md1,md3,lot,nfft,n3,zw(1,1,inzee),zf(i1,1,j2)&
                ,scal)!,ehartreetmp)

           !integrate local pieces together
           !ehartree=ehartree+0.5d0*ehartreetmp*hx*hy*hz
        end do
     endif
  end do

  !De-allocations  
  i_all=-product(shape(btrig1))*kind(btrig1)
  deallocate(btrig1,stat=i_stat)
  call memocc(i_stat,i_all,'btrig1','s_poissonsolver')
  i_all=-product(shape(ftrig1))*kind(ftrig1)
  deallocate(ftrig1,stat=i_stat)
  call memocc(i_stat,i_all,'ftrig1','s_poissonsolver')
  i_all=-product(shape(after1))*kind(after1)
  deallocate(after1,stat=i_stat)
  call memocc(i_stat,i_all,'after1','s_poissonsolver')
  i_all=-product(shape(now1))*kind(now1)
  deallocate(now1,stat=i_stat)
  call memocc(i_stat,i_all,'now1','s_poissonsolver')
  i_all=-product(shape(before1))*kind(before1)
  deallocate(before1,stat=i_stat)
  call memocc(i_stat,i_all,'before1','s_poissonsolver')
  i_all=-product(shape(btrig2))*kind(btrig2)
  deallocate(btrig2,stat=i_stat)
  call memocc(i_stat,i_all,'btrig2','s_poissonsolver')
  i_all=-product(shape(ftrig2))*kind(ftrig2)
  deallocate(ftrig2,stat=i_stat)
  call memocc(i_stat,i_all,'ftrig2','s_poissonsolver')
  i_all=-product(shape(after2))*kind(after2)
  deallocate(after2,stat=i_stat)
  call memocc(i_stat,i_all,'after2','s_poissonsolver')
  i_all=-product(shape(now2))*kind(now2)
  deallocate(now2,stat=i_stat)
  call memocc(i_stat,i_all,'now2','s_poissonsolver')
  i_all=-product(shape(before2))*kind(before2)
  deallocate(before2,stat=i_stat)
  call memocc(i_stat,i_all,'before2','s_poissonsolver')
  i_all=-product(shape(btrig3))*kind(btrig3)
  deallocate(btrig3,stat=i_stat)
  call memocc(i_stat,i_all,'btrig3','s_poissonsolver')
  i_all=-product(shape(ftrig3))*kind(ftrig3)
  deallocate(ftrig3,stat=i_stat)
  call memocc(i_stat,i_all,'ftrig3','s_poissonsolver')
  i_all=-product(shape(after3))*kind(after3)
  deallocate(after3,stat=i_stat)
  call memocc(i_stat,i_all,'after3','s_poissonsolver')
  i_all=-product(shape(now3))*kind(now3)
  deallocate(now3,stat=i_stat)
  call memocc(i_stat,i_all,'now3','s_poissonsolver')
  i_all=-product(shape(before3))*kind(before3)
  deallocate(before3,stat=i_stat)
  call memocc(i_stat,i_all,'before3','s_poissonsolver')
  i_all=-product(shape(zmpi2))*kind(zmpi2)
  deallocate(zmpi2,stat=i_stat)
  call memocc(i_stat,i_all,'zmpi2','s_poissonsolver')
  i_all=-product(shape(zw))*kind(zw)
  deallocate(zw,stat=i_stat)
  call memocc(i_stat,i_all,'zw','s_poissonsolver')
  i_all=-product(shape(zt))*kind(zt)
  deallocate(zt,stat=i_stat)
  call memocc(i_stat,i_all,'zt','s_poissonsolver')
  i_all=-product(shape(cosinarr))*kind(cosinarr)
  deallocate(cosinarr,stat=i_stat)
  call memocc(i_stat,i_all,'cosinarr','s_poissonsolver')
  if (nproc.gt.1) then
     i_all=-product(shape(zmpi1))*kind(zmpi1)
     deallocate(zmpi1,stat=i_stat)
     call memocc(i_stat,i_all,'zmpi1','s_poissonsolver')
  end if

  call timing(iproc,'PSolv_comput  ','OF')

end subroutine S_PoissonSolver
!!**


subroutine S_mpiswitch_upcorn(j3,nfft,Jp2stb,J2stb,lot,n1,md2,nd3,nproc,zmpi1,zw)
  implicit none
  integer, intent(in) :: j3,nfft,lot,n1,md2,nd3,nproc
  integer, intent(inout) :: Jp2stb,J2stb
  real(kind=8), dimension(2,n1,md2/nproc,nd3/nproc,nproc), intent(in) :: zmpi1
  real(kind=8), dimension(2,lot,n1), intent(out) :: zw
  !local variables
  integer :: I1,J2,Jp2,mfft

  mfft=0
  do Jp2=Jp2stb,nproc
     do J2=J2stb,md2/nproc
      mfft=mfft+1
      if (mfft.gt.nfft) then
           Jp2stb=Jp2
           J2stb=J2
           return
      end if
        do I1=1,n1
           zw(1,mfft,I1)=zmpi1(1,I1,J2,j3,Jp2)
           zw(2,mfft,I1)=zmpi1(2,I1,J2,j3,Jp2)
        end do
     end do
     J2stb=1
  end do

end subroutine S_mpiswitch_upcorn

subroutine S_switch_upcorn(nfft,n2,lot,n1,lzt,zt,zw)
  implicit none
  integer, intent(in) :: nfft,n2,lot,n1,lzt
  real(kind=8), dimension(2,lzt,n1), intent(in) :: zt
  real(kind=8), dimension(2,lot,n2), intent(out) :: zw
  !local variables
  integer :: i,j

  do j=1,nfft
     do i=1,n2
      zw(1,j,i)=zt(1,i,j)
      zw(2,j,i)=zt(2,i,j)
     end do
  end do

end subroutine S_switch_upcorn

subroutine S_unswitch_downcorn(nfft,n2,lot,n1,lzt,zw,zt)
  implicit none
  integer, intent(in) :: nfft,n2,lot,n1,lzt
  real(kind=8), dimension(2,lot,n2), intent(in) :: zw
  real(kind=8), dimension(2,lzt,n1), intent(out) :: zt
  !local variables
  integer :: i,j

  do j=1,nfft
     do i=1,n2
      zt(1,i,j)=zw(1,j,i)
      zt(2,i,j)=zw(2,j,i)
     end do
  end do

end subroutine S_unswitch_downcorn

subroutine S_unmpiswitch_downcorn(j3,nfft,Jp2stf,J2stf,lot,n1,md2,nd3,nproc,zw,zmpi1)
  implicit none
  integer, intent(in) :: j3,nfft,lot,n1,md2,nd3,nproc
  integer, intent(inout) :: Jp2stf,J2stf
  real(kind=8), dimension(2,lot,n1), intent(in) :: zw
  real(kind=8), dimension(2,n1,md2/nproc,nd3/nproc,nproc), intent(out) :: zmpi1
  !local variables
  integer :: I1,J2,Jp2,mfft

  mfft=0
  do Jp2=Jp2stf,nproc
     do J2=J2stf,md2/nproc
      mfft=mfft+1
      if (mfft.gt.nfft) then
           Jp2stf=Jp2
           J2stf=J2
           return
      end if
        do I1=1,n1
           zmpi1(1,I1,J2,j3,Jp2)=zw(1,mfft,I1)
           zmpi1(2,I1,J2,j3,Jp2)=zw(2,mfft,I1)
        end do
     end do
     J2stf=1
  end do
end subroutine S_unmpiswitch_downcorn


!!****f* BigDFT/unfill_downcorn
!! NAME
!!   unfill_downcorn
!!
!! FUNCTION
!!     (Based on suitable modifications of S.Goedecker routines)
!!     Restore data into output array, calculating in the meanwhile
!!     Hartree energy of the potential 
!!
!! SYNOPSIS
!!     zf:          Original distributed density as well as
!!                  Distributed solution of the poisson equation (inout)
!!     zw:          FFT work array
!!     n3:          (twice the) dimension of the last FFTtransform.
!!     md1,md3:     Dimensions of the undistributed part of the real grid
!!     nfft:        number of planes
!!     scal:        Needed to achieve unitarity and correct dimensions
!!     ehartreetmp: Hartree energy
!!
!! WARNING
!!     Assuming that high frequencies are in the corners 
!!     and that n3 is multiple of 4   
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
subroutine unfill_downcorn(md1,md3,lot,nfft,n3,zw,zf&
     ,scal)!,ehartreetmp)
  implicit none
  !Arguments
  integer, intent(in) :: md1,md3,lot,nfft,n3
  real(kind=8), dimension(2,lot,n3/2), intent(in) :: zw
  real(kind=8), dimension(md1,md3),intent(inout) :: zf
  real(kind=8), intent(in) :: scal
  !real(kind=8), intent(out) :: ehartreetmp
  !Local variables
  integer :: i3,i1
  real(kind=8) :: pot1

  !Body
  !ehartreetmp=0.d0
  do i3=1,n3/4
     do i1=1,nfft
        pot1 = scal*zw(1,i1,i3)
        !ehartreetmp =ehartreetmp + pot1* zf(i1,2*i3-1)
      zf(i1,2*i3-1)= pot1 
      pot1 = scal*zw(2,i1,i3)
        !ehartreetmp =ehartreetmp + pot1* zf(i1,2*i3)
      zf(i1,2*i3)= pot1 
     enddo
  end do
  
end subroutine unfill_downcorn
!!**


subroutine halfill_upcorn(md1,md3,lot,nfft,n3,zf,zw)
  implicit real(kind=8) (a-h,o-z)
  dimension zw(2,lot,n3/2),zf(md1,md3)
! WARNING: Assuming that high frequencies are in the corners 
!          and that n3 is multiple of 4
!in principle we can relax this condition
      
  do i3=1,n3/4
     do i1=1,nfft
        zw(1,i1,i3)=0.d0
        zw(2,i1,i3)=0.d0
     end do
  end do
  do i3=n3/4+1,n3/2
     do i1=1,nfft
        zw(1,i1,i3)=zf(i1,2*i3-1-n3/2)
        zw(2,i1,i3)=zf(i1,2*i3-n3/2)
     end do
  end do
      
end subroutine halfill_upcorn


!!****f* BigDFT/scramble_unpack
!! NAME
!!   scramble_unpack
!!
!! FUNCTION
!!     (Based on suitable modifications of S.Goedecker routines)
!!     Assign the correct planes to the work array zmpi2
!!     in order to prepare for interprocessor data transposition.
!!     In the meanwhile, it unpacks the data of the HalFFT in order to prepare for
!!     multiplication with the kernel
!!
!! SYNOPSIS
!!     zmpi2:          Work array for multiprocessor manipulation (output)
!!     zw:             Work array (input)
!!     cosinarr:      Array of the phases needed for unpacking
!!     n1,n3:          logical dimension of the FFT transform, reference for work arrays
!!     md2,nd3:        Dimensions of real grid and of the kernel, respectively
!!     i1,j2,lot,nfft: Starting points of the plane and number of remaining lines
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
subroutine scramble_unpack(i1,j2,lot,nfft,n1,n3,md2,nproc,nd3,zw,zmpi2,cosinarr)
  implicit none
  !Arguments
  integer, intent(in) :: i1,j2,lot,nfft,n1,n3,md2,nproc,nd3
  real(kind=8), dimension(2,lot,n3/2), intent(in) :: zw
  real(kind=8), dimension(2,n3/2), intent(in) :: cosinarr
  real(kind=8), dimension(2,n1,md2/nproc,nd3), intent(out) :: zmpi2
  !Local variables
  integer :: i3,i,ind1,ind2
  real(kind=8) ::  a,b,c,d,cp,sp,feR,feI,foR,foI,fR,fI
  
  !Body

  !case i3=1 and i3=n3/2+1
  do i=0,nfft-1
     a=zw(1,i+1,1)
     b=zw(2,i+1,1)
     zmpi2(1,i1+i,j2,1)=a+b
     zmpi2(2,i1+i,j2,1)=0.d0
     zmpi2(1,i1+i,j2,n3/2+1)=a-b
     zmpi2(2,i1+i,j2,n3/2+1)=0.d0
  end do
  !case 2<=i3<=n3/2
  do i3=2,n3/2
     ind1=i3
     ind2=n3/2-i3+2
     cp=cosinarr(1,i3)
     sp=cosinarr(2,i3)
     do i=0,nfft-1
        a=zw(1,i+1,ind1)
        b=zw(2,i+1,ind1)
        c=zw(1,i+1,ind2)
        d=zw(2,i+1,ind2)
        feR=.5d0*(a+c)
        feI=.5d0*(b-d)
        foR=.5d0*(a-c)
        foI=.5d0*(b+d) 
        fR=feR+cp*foI-sp*foR
        fI=feI-cp*foR-sp*foI
        zmpi2(1,i1+i,j2,ind1)=fR
        zmpi2(2,i1+i,j2,ind1)=fI
     end do
  end do

end subroutine scramble_unpack
!!**

 
!!****f* BigDFT/unscramble_pack
!! NAME
!!   unscramble_pack
!!
!! FUNCTION
!!     (Based on suitable modifications of S.Goedecker routines)
!!     Insert the correct planes of the work array zmpi2
!!     in order to prepare for backward FFT transform
!!     In the meanwhile, it packs the data in order to be transformed with the HalFFT 
!!     procedure
!!
!! SYNOPSIS
!!     zmpi2:          Work array for multiprocessor manipulation (input)
!!     zw:             Work array (output)
!!     cosinarr:       Array of the phases needed for packing
!!     n1,n3:          logical dimension of the FFT transform, reference for work arrays
!!     md2,nd3:        Dimensions of real grid and of the kernel, respectively
!!     i1,j2,lot,nfft: Starting points of the plane and number of remaining lines
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
subroutine unscramble_pack(i1,j2,lot,nfft,n1,n3,md2,nproc,nd3,zmpi2,zw,cosinarr)
  implicit none
  !Arguments
  integer, intent(in) :: i1,j2,lot,nfft,n1,n3,md2,nproc,nd3
  real(kind=8), dimension(2,lot,n3/2), intent(out) :: zw
  real(kind=8), dimension(2,n3/2), intent(in) :: cosinarr
  real(kind=8), dimension(2,n1,md2/nproc,nd3), intent(in) :: zmpi2
  !Local variables
  integer :: i3,i,indA,indB
  real(kind=8) ::  a,b,c,d,cp,sp,re,ie,ro,io,rh,ih

  !Body

  do i3=1,n3/2
     indA=i3
     indB=n3/2+2-i3
     cp=cosinarr(1,i3)
     sp=cosinarr(2,i3)
     do i=0,nfft-1
        a=zmpi2(1,i1+i,j2,indA)
        b=zmpi2(2,i1+i,j2,indA)
        c= zmpi2(1,i1+i,j2,indB)
        d=-zmpi2(2,i1+i,j2,indB)
        re=(a+c)
        ie=(b+d)
        ro=(a-c)*cp-(b-d)*sp
        io=(a-c)*sp+(b-d)*cp
        rh=re-io 
        ih=ie+ro
        zw(1,i+1,indA)=rh
        zw(2,i+1,indA)=ih
     end do
  end do

end subroutine unscramble_pack
!!**


!!****f* BigDFT/F_PoissonSolver
!! NAME
!!   F_PoissonSolver
!!
!! FUNCTION
!!     (Based on suitable modifications of S.Goedecker routines)
!!     Applies the local FFT space Kernel to the density in Real space.
!!     Calculates also the LDA exchange-correlation terms
!!
!! SYNOPSIS
!!     zf:          Density (input/output)
!!                  ZF(i1,i3,i2)
!!                  i1=1,md1 , i2=1,md2/nproc , i3=1,md3 
!!     pot:         Kernel, only the distributed part (REAL)
!!                  POT(i1,i2,i3)
!!                  i1=1,nd1 , i2=1,nd2 , i3=1,nd3/nproc
!!     zfpot_ion:   Distributed array of the ionization potential 
!!     nproc:       number of processors used as returned by MPI_COMM_SIZE
!!     iproc:       [0:nproc-1] number of processor as returned by MPI_COMM_RANK
!!     n1,n2,n3:    logical dimension of the transform. As transform lengths 
!!                  most products of the prime factors 2,3,5 are allowed.
!!                  The detailed table with allowed transform lengths can 
!!                  be found in subroutine CTRIG
!!     md1,md2,md3: Dimension of ZF
!!     nd1,nd2,nd3: Dimension of POT
!!     scal:        factor of renormalization of the FFT in order to acheve unitarity 
!!                  and the correct dimension
!!     hgrid:       grid spacing, used for integrating quantities
!!     ehartree:    hartree energy of the potential
!!     eexcu,vexcu: LDA exchange correlation terms   
!!
!! RESTRICTIONS on USAGE
!!     Copyright (C) Stefan Goedecker, Cornell University, Ithaca, USA, 1994
!!     Copyright (C) Stefan Goedecker, MPI Stuttgart, Germany, 1999
!!     Copyright (C) 2002 Stefan Goedecker, CEA Grenoble
!!     This file is distributed under the terms of the
!!      GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .
!!
!! AUTHORS
!!    S. Goedecker, L. Genovese
!!
!! CREATION DATE
!!     February 2006
!!
!! SOURCE
!!
subroutine F_PoissonSolver(n1,n2,n3,nd1,nd2,nd3,md1,md2,md3,nproc,iproc,pot,zf&
             ,scal)!,hgrid)!,ehartree)
  use mpi
  implicit none
  include 'perfdata.inc'
  !Arguments
  integer, intent(in) :: n1,n2,n3,nd1,nd2,nd3,md1,md2,md3,nproc,iproc
  real(kind=8), intent(in) :: scal!,hgrid
  !real(kind=8), intent(out) :: ehartree
  real(kind=8), dimension(nd1,nd2,nd3/nproc), intent(in) :: pot
  real(kind=8), dimension(md1,md3,md2/nproc), intent(inout) :: zf
  !Local variables
  !Maximum number of points for FFT (should be same number in fft3d routine)
  integer, parameter :: nfft_max=24000
  integer :: ncache,lzt,lot,ma,mb,nfft,ic1,ic2,ic3,Jp2stb,J2stb,Jp2stf,J2stf
  integer :: j2,j3,i1,i3,i,j,inzee,ierr,i_all,i_stat
  real(kind=8) :: twopion!,ehartreetmp
  !work arrays for transpositions
  real(kind=8), dimension(:,:,:), allocatable :: zt
  !work arrays for MPI
  real(kind=8), dimension(:,:,:,:,:), allocatable :: zmpi1
  real(kind=8), dimension(:,:,:,:), allocatable :: zmpi2
  !cache work array
  real(kind=8), dimension(:,:,:), allocatable :: zw
  !FFT work arrays
  real(kind=8), dimension(:,:), allocatable :: btrig1,btrig2,btrig3,&
       ftrig1,ftrig2,ftrig3,cosinarr
  integer, dimension(:), allocatable :: after1,now1,before1, & 
       after2,now2,before2,after3,now3,before3

  call timing(iproc,'PSolv_comput  ','ON')
  
  !Body
  ! check input
  if (mod(n1,2).ne.0) stop 'Parallel convolution:ERROR:n1'
  if (mod(n2,2).ne.0) stop 'Parallel convolution:ERROR:n2'
  if (mod(n3,2).ne.0) stop 'Parallel convolution:ERROR:n3'
  if (nd1.lt.n1/2+1) stop 'Parallel convolution:ERROR:nd1'
  if (nd2.lt.n2/2+1) stop 'Parallel convolution:ERROR:nd2'
  if (nd3.lt.n3/2+1) stop 'Parallel convolution:ERROR:nd3'
  if (md1.lt.n1/2) stop 'Parallel convolution:ERROR:md1'
  if (md2.lt.n2/2) stop 'Parallel convolution:ERROR:md2'
  if (md3.lt.n3/2) stop 'Parallel convolution:ERROR:md3'
  if (mod(nd3,nproc).ne.0) stop 'Parallel convolution:ERROR:nd3'
  if (mod(md2,nproc).ne.0) stop 'Parallel convolution:ERROR:md2'
  
  !defining work arrays dimensions
  
  ncache=ncache_optimal
  if (ncache <= max(n1,n2,n3/2)*4) ncache=max(n1,n2,n3/2)*4
  lzt=n2/2
  if (mod(n2/2,2).eq.0) lzt=lzt+1
  if (mod(n2/2,4).eq.0) lzt=lzt+1
  
  !Allocations
  allocate(btrig1(2,nfft_max),stat=i_stat)
  call memocc(i_stat,product(shape(btrig1))*kind(btrig1),'btrig1','f_poissonsolver')
  allocate(ftrig1(2,nfft_max),stat=i_stat)
  call memocc(i_stat,product(shape(ftrig1))*kind(ftrig1),'ftrig1','f_poissonsolver')
  allocate(after1(7),stat=i_stat)
  call memocc(i_stat,product(shape(after1))*kind(after1),'after1','f_poissonsolver')
  allocate(now1(7),stat=i_stat)
  call memocc(i_stat,product(shape(now1))*kind(now1),'now1','f_poissonsolver')
  allocate(before1(7),stat=i_stat)
  call memocc(i_stat,product(shape(before1))*kind(before1),'before1','f_poissonsolver')
  allocate(btrig2(2,nfft_max),stat=i_stat)
  call memocc(i_stat,product(shape(btrig2))*kind(btrig2),'btrig2','f_poissonsolver')
  allocate(ftrig2(2,nfft_max),stat=i_stat)
  call memocc(i_stat,product(shape(ftrig2))*kind(ftrig2),'ftrig2','f_poissonsolver')
  allocate(after2(7),stat=i_stat)
  call memocc(i_stat,product(shape(after2))*kind(after2),'after2','f_poissonsolver')
  allocate(now2(7),stat=i_stat)
  call memocc(i_stat,product(shape(now2))*kind(now2),'now2','f_poissonsolver')
  allocate(before2(7),stat=i_stat)
  call memocc(i_stat,product(shape(before2))*kind(before2),'before2','f_poissonsolver')
  allocate(btrig3(2,nfft_max),stat=i_stat)
  call memocc(i_stat,product(shape(btrig3))*kind(btrig3),'btrig3','f_poissonsolver')
  allocate(ftrig3(2,nfft_max),stat=i_stat)
  call memocc(i_stat,product(shape(ftrig3))*kind(ftrig3),'ftrig3','f_poissonsolver')
  allocate(after3(7),stat=i_stat)
  call memocc(i_stat,product(shape(after3))*kind(after3),'after3','f_poissonsolver')
  allocate(now3(7),stat=i_stat)
  call memocc(i_stat,product(shape(now3))*kind(now3),'now3','f_poissonsolver')
  allocate(before3(7),stat=i_stat)
  call memocc(i_stat,product(shape(before3))*kind(before3),'before3','f_poissonsolver')
  allocate(zw(2,ncache/4,2),stat=i_stat)
  call memocc(i_stat,product(shape(zw))*kind(zw),'zw','f_poissonsolver')
  allocate(zt(2,lzt,n1),stat=i_stat)
  call memocc(i_stat,product(shape(zt))*kind(zt),'zt','f_poissonsolver')
  allocate(zmpi2(2,n1,md2/nproc,nd3),stat=i_stat)
  call memocc(i_stat,product(shape(zmpi2))*kind(zmpi2),'zmpi2','f_poissonsolver')
  allocate(cosinarr(2,n3/2),stat=i_stat)
  call memocc(i_stat,product(shape(cosinarr))*kind(cosinarr),'cosinarr','f_poissonsolver')
  if (nproc.gt.1) then
     allocate(zmpi1(2,n1,md2/nproc,nd3/nproc,nproc),stat=i_stat)
     call memocc(i_stat,product(shape(zmpi1))*kind(zmpi1),'zmpi1','f_poissonsolver')
  end if


  !calculating the FFT work arrays (beware on the HalFFT in n3 dimension)
  call ctrig(n3/2,btrig3,after3,before3,now3,1,ic3)
  call ctrig(n1,btrig1,after1,before1,now1,1,ic1)
  call ctrig(n2,btrig2,after2,before2,now2,1,ic2)
  do  j=1,n1
     ftrig1(1,j)= btrig1(1,j)
     ftrig1(2,j)=-btrig1(2,j)
  enddo
  do  j=1,n2
     ftrig2(1,j)= btrig2(1,j)
     ftrig2(2,j)=-btrig2(2,j)
  enddo
  do  j=1,n3
     ftrig3(1,j)= btrig3(1,j)
     ftrig3(2,j)=-btrig3(2,j)
  enddo

  !Calculating array of phases for HalFFT decoding
  twopion=8.d0*datan(1.d0)/real(n3,kind=8)
  do i3=1,n3/2
     cosinarr(1,i3)= dcos(twopion*real(i3-1,kind=8))
     cosinarr(2,i3)=-dsin(twopion*real(i3-1,kind=8))
  end do

  !initializing integrals
  !ehartree=0.d0

  ! transform along z axis
  lot=ncache/(2*n3)
  if (lot.lt.1) then  
     write(6,*) & 
          'convolxc_on:ncache has to be enlarged to be able to hold at' // &  
          'least one 1-d FFT of this size even though this will' // & 
          'reduce the performance for shorter transform lengths'
     stop
  endif
  
  do j2=1,md2/nproc
     !this condition ensures that we manage only the interesting part for the FFT
     if (iproc*(md2/nproc)+j2.le.n2/2) then
        do i1=1,(n1/2),lot
           ma=i1
           mb=min(i1+(lot-1),(n1/2))
           nfft=mb-ma+1

           !inserting real data into complex array of half lenght
           call halfill_upcorn(md1,md3,lot,nfft,n3,zf(i1,1,j2),zw(1,1,1))

           !performing FFT
           !input: I1,I3,J2,(Jp2)
           inzee=1
           do i=1,ic3
              call fftstp(lot,nfft,n3/2,lot,n3/2,zw(1,1,inzee),zw(1,1,3-inzee), &
                   btrig3,after3(i),now3(i),before3(i),1)
              inzee=3-inzee
           enddo
           !output: I1,i3,J2,(Jp2)

           !unpacking FFT in order to restore correct result, 
           !while exchanging components
           !input: I1,i3,J2,(Jp2)
           call scramble_unpack(i1,j2,lot,nfft,n1/2,n3,md2,nproc,nd3,zw(1,1,inzee),zmpi2,cosinarr)
           !output: I1,J2,i3,(Jp2)
        end do
     endif
  end do

  !Interprocessor data transposition
  !input: I1,J2,j3,jp3,(Jp2)
  if (nproc.gt.1) then
     call timing(iproc,'PSolv_comput  ','OF')

     call timing(iproc,'PSolv_commun  ','ON')
     !communication scheduling
     call MPI_ALLTOALL(zmpi2,n1*(md2/nproc)*(nd3/nproc), &
          MPI_double_precision, &
          zmpi1,n1*(md2/nproc)*(nd3/nproc), &
          MPI_double_precision,MPI_COMM_WORLD,ierr)
     call timing(iproc,'PSolv_commun  ','OF')

     call timing(iproc,'PSolv_comput  ','ON')
  endif
  !output: I1,J2,j3,Jp2,(jp3)

  !now each process perform complete convolution of its planes
  do j3=1,nd3/nproc
     !this condition ensures that we manage only the interesting part for the FFT
     if (iproc*(nd3/nproc)+j3.le.n3/2+1) then
      Jp2stb=1
      J2stb=1
      Jp2stf=1
      J2stf=1
        
        ! transform along x axis
        lot=ncache/(4*n1)
        if (lot.lt.1) then  
           write(6,*) & 
                'convolxc_on:ncache has to be enlarged to be able to hold at' // &  
                'least one 1-d FFT of this size even though this will' // & 
                'reduce the performance for shorter transform lengths'
           stop
        endif
        
        do j=1,n2/2,lot
           ma=j
           mb=min(j+(lot-1),n2/2)
           nfft=mb-ma+1

           !reverse index ordering, leaving the planes to be transformed at the end
           !input: I1,J2,j3,Jp2,(jp3)
           if (nproc.eq.1) then
              call mpiswitch_upcorn(j3,nfft,Jp2stb,J2stb,lot,n1,md2,nd3,nproc,zmpi2,zw(1,1,1))
           else
              call mpiswitch_upcorn(j3,nfft,Jp2stb,J2stb,lot,n1,md2,nd3,nproc,zmpi1,zw(1,1,1))
           endif
           !output: J2,Jp2,I1,j3,(jp3)
           
           !performing FFT
           !input: I2,I1,j3,(jp3)
           inzee=1
           do i=1,ic1-1
              call fftstp(lot,nfft,n1,lot,n1,zw(1,1,inzee),zw(1,1,3-inzee), &
                   btrig1,after1(i),now1(i),before1(i),1)
              inzee=3-inzee
           enddo

           !storing the last step into zt array
           i=ic1
           call fftstp(lot,nfft,n1,lzt,n1,zw(1,1,inzee),zt(1,j,1), & 
                btrig1,after1(i),now1(i),before1(i),1)           
           !output: I2,i1,j3,(jp3)
        end do

        !transform along y axis
        lot=ncache/(4*n2)
        if (lot.lt.1) then  
           write(6,*) & 
                'convolxc_on:ncache has to be enlarged to be able to hold at' // &  
                'least one 1-d FFT of this size even though this will' // & 
                'reduce the performance for shorter transform lengths'
           stop
        endif

        do j=1,n1,lot
           ma=j
           mb=min(j+(lot-1),n1)
           nfft=mb-ma+1

           !reverse ordering 
           !input: I2,i1,j3,(jp3)
           call switch_upcorn(nfft,n2,lot,n1,lzt,zt(1,1,j),zw(1,1,1))
           !output: i1,I2,j3,(jp3)
           
           !performing FFT
           !input: i1,I2,j3,(jp3)
           inzee=1
           do i=1,ic2
              call fftstp(lot,nfft,n2,lot,n2,zw(1,1,inzee),zw(1,1,3-inzee), &
                   btrig2,after2(i),now2(i),before2(i),1)
              inzee=3-inzee
           enddo
           !output: i1,i2,j3,(jp3)
           
           !Multiply with kernel in fourier space
           call multkernel(nd1,nd2,n1,n2,lot,nfft,j,pot(1,1,j3),zw(1,1,inzee))
           
           !TRANSFORM BACK IN REAL SPACE
           
           !transform along y axis
           !input: i1,i2,j3,(jp3)
           do i=1,ic2
              call fftstp(lot,nfft,n2,lot,n2,zw(1,1,inzee),zw(1,1,3-inzee), &
                   ftrig2,after2(i),now2(i),before2(i),-1)
              inzee=3-inzee
           enddo

           !reverse ordering
           !input: i1,I2,j3,(jp3)
           call unswitch_downcorn(nfft,n2,lot,n1,lzt,zw(1,1,inzee),zt(1,1,j))
           !output: I2,i1,j3,(jp3)
        end do
        
        !transform along x axis
        !input: I2,i1,j3,(jp3)
        lot=ncache/(4*n1)
        do j=1,n2/2,lot
           ma=j
           mb=min(j+(lot-1),n2/2)
           nfft=mb-ma+1

           !performing FFT
           i=1
           call fftstp(lzt,nfft,n1,lot,n1,zt(1,j,1),zw(1,1,1), &
                ftrig1,after1(i),now1(i),before1(i),-1)
           
           inzee=1
           do i=2,ic1
              call fftstp(lot,nfft,n1,lot,n1,zw(1,1,inzee),zw(1,1,3-inzee), &
                   ftrig1,after1(i),now1(i),before1(i),-1)
              inzee=3-inzee
           enddo
           !output: I2,I1,j3,(jp3)

           !reverse ordering
           !input: J2,Jp2,I1,j3,(jp3)
           if (nproc.eq.1) then
              call unmpiswitch_downcorn(j3,nfft,Jp2stf,J2stf,lot,n1,md2,nd3,nproc,zw(1,1,inzee),zmpi2)
           else
              call unmpiswitch_downcorn(j3,nfft,Jp2stf,J2stf,lot,n1,md2,nd3,nproc,zw(1,1,inzee),zmpi1)
           endif
           ! output: I1,J2,j3,Jp2,(jp3)
        end do
     endif
  end do
  
  
  !Interprocessor data transposition
  !input: I1,J2,j3,Jp2,(jp3)
  if (nproc.gt.1) then
     call timing(iproc,'PSolv_comput  ','OF')

     call timing(iproc,'PSolv_commun  ','ON')
     !communication scheduling
     call MPI_ALLTOALL(zmpi1,n1*(md2/nproc)*(nd3/nproc), &
          MPI_double_precision, &
          zmpi2,n1*(md2/nproc)*(nd3/nproc), &
          MPI_double_precision,MPI_COMM_WORLD,ierr)
     call timing(iproc,'PSolv_commun  ','OF')

     call timing(iproc,'PSolv_comput  ','ON')
     !output: I1,J2,j3,jp3,(Jp2)
  endif

  !transform along z axis
  !input: I1,J2,i3,(Jp2)
  lot=ncache/(2*n3)
  do j2=1,md2/nproc
     !this condition ensures that we manage only the interesting part for the FFT
     if (iproc*(md2/nproc)+j2.le.n2/2) then
        do i1=1,(n1/2),lot
           ma=i1
           mb=min(i1+(lot-1),(n1/2))
           nfft=mb-ma+1

           !reverse ordering and repack the FFT data in order to be backward HalFFT transformed
           !input: I1,J2,i3,(Jp2)
           call unscramble_pack(i1,j2,lot,nfft,n1/2,n3,md2,nproc,nd3,zmpi2,zw(1,1,1),cosinarr)
           !output: I1,i3,J2,(Jp2)

           !performing FFT
           !input: I1,i3,J2,(Jp2)           
           inzee=1
           do i=1,ic3
              call fftstp(lot,nfft,n3/2,lot,n3/2,zw(1,1,inzee),zw(1,1,3-inzee), &
                   ftrig3,after3(i),now3(i),before3(i),-1)
              inzee=3-inzee
           enddo
           !output: I1,I3,J2,(Jp2)

           !calculates the exchange correlation terms locally and rebuild the output array
           call unfill_downcorn(md1,md3,lot,nfft,n3,zw(1,1,inzee),zf(i1,1,j2)&
                ,scal)!,ehartreetmp)

           !integrate local pieces together
           !ehartree=ehartree+0.5d0*ehartreetmp*hgrid**3
        end do
     endif
  end do
  
  !De-allocations  
  i_all=-product(shape(btrig1))*kind(btrig1)
  deallocate(btrig1,stat=i_stat)
  call memocc(i_stat,i_all,'btrig1','f_poissonsolver')
  i_all=-product(shape(ftrig1))*kind(ftrig1)
  deallocate(ftrig1,stat=i_stat)
  call memocc(i_stat,i_all,'ftrig1','f_poissonsolver')
  i_all=-product(shape(after1))*kind(after1)
  deallocate(after1,stat=i_stat)
  call memocc(i_stat,i_all,'after1','f_poissonsolver')
  i_all=-product(shape(now1))*kind(now1)
  deallocate(now1,stat=i_stat)
  call memocc(i_stat,i_all,'now1','f_poissonsolver')
  i_all=-product(shape(before1))*kind(before1)
  deallocate(before1,stat=i_stat)
  call memocc(i_stat,i_all,'before1','f_poissonsolver')
  i_all=-product(shape(btrig2))*kind(btrig2)
  deallocate(btrig2,stat=i_stat)
  call memocc(i_stat,i_all,'btrig2','f_poissonsolver')
  i_all=-product(shape(ftrig2))*kind(ftrig2)
  deallocate(ftrig2,stat=i_stat)
  call memocc(i_stat,i_all,'ftrig2','f_poissonsolver')
  i_all=-product(shape(after2))*kind(after2)
  deallocate(after2,stat=i_stat)
  call memocc(i_stat,i_all,'after2','f_poissonsolver')
  i_all=-product(shape(now2))*kind(now2)
  deallocate(now2,stat=i_stat)
  call memocc(i_stat,i_all,'now2','f_poissonsolver')
  i_all=-product(shape(before2))*kind(before2)
  deallocate(before2,stat=i_stat)
  call memocc(i_stat,i_all,'before2','f_poissonsolver')
  i_all=-product(shape(btrig3))*kind(btrig3)
  deallocate(btrig3,stat=i_stat)
  call memocc(i_stat,i_all,'btrig3','f_poissonsolver')
  i_all=-product(shape(ftrig3))*kind(ftrig3)
  deallocate(ftrig3,stat=i_stat)
  call memocc(i_stat,i_all,'ftrig3','f_poissonsolver')
  i_all=-product(shape(after3))*kind(after3)
  deallocate(after3,stat=i_stat)
  call memocc(i_stat,i_all,'after3','f_poissonsolver')
  i_all=-product(shape(now3))*kind(now3)
  deallocate(now3,stat=i_stat)
  call memocc(i_stat,i_all,'now3','f_poissonsolver')
  i_all=-product(shape(before3))*kind(before3)
  deallocate(before3,stat=i_stat)
  call memocc(i_stat,i_all,'before3','f_poissonsolver')
  i_all=-product(shape(zmpi2))*kind(zmpi2)
  deallocate(zmpi2,stat=i_stat)
  call memocc(i_stat,i_all,'zmpi2','f_poissonsolver')
  i_all=-product(shape(zw))*kind(zw)
  deallocate(zw,stat=i_stat)
  call memocc(i_stat,i_all,'zw','f_poissonsolver')
  i_all=-product(shape(zt))*kind(zt)
  deallocate(zt,stat=i_stat)
  call memocc(i_stat,i_all,'zt','f_poissonsolver')
  i_all=-product(shape(cosinarr))*kind(cosinarr)
  deallocate(cosinarr,stat=i_stat)
  call memocc(i_stat,i_all,'cosinarr','f_poissonsolver')
  if (nproc.gt.1) then
     i_all=-product(shape(zmpi1))*kind(zmpi1)
     deallocate(zmpi1,stat=i_stat)
     call memocc(i_stat,i_all,'zmpi1','f_poissonsolver')
  end if

  call timing(iproc,'PSolv_comput  ','OF')
end subroutine F_PoissonSolver
!!**


subroutine switch_upcorn(nfft,n2,lot,n1,lzt,zt,zw)
  implicit real(kind=8) (a-h,o-z)
  dimension zw(2,lot,n2),zt(2,lzt,n1)
! WARNING: Assuming that high frequencies are in the corners 
!          and that n2 is multiple of 2

! Low frequencies 
  do j=1,nfft
     do i=n2/2+1,n2
        zw(1,j,i)=zt(1,i-n2/2,j)
        zw(2,j,i)=zt(2,i-n2/2,j)
     end do
  end do

! High frequencies 
  do i=1,n2/2
     do j=1,nfft
        zw(1,j,i)=0.d0
        zw(2,j,i)=0.d0
     end do
  end do

end subroutine switch_upcorn

        
subroutine mpiswitch_upcorn(j3,nfft,Jp2stb,J2stb,lot,n1,md2,nd3,nproc,zmpi1,zw)
  implicit real(kind=8) (a-h,o-z)
  dimension zmpi1(2,n1/2,md2/nproc,nd3/nproc,nproc),zw(2,lot,n1)
! WARNING: Assuming that high frequencies are in the corners 
!          and that n1 is multiple of 2

  mfft=0
  do Jp2=Jp2stb,nproc
     do J2=J2stb,md2/nproc
        mfft=mfft+1
        if (mfft.gt.nfft) then
        Jp2stb=Jp2
        J2stb=J2
        return
        endif
        do I1=1,n1/2
           zw(1,mfft,I1)=0.d0
           zw(2,mfft,I1)=0.d0
        end do
        do I1=n1/2+1,n1
           zw(1,mfft,I1)=zmpi1(1,I1-n1/2,J2,j3,Jp2)
           zw(2,mfft,I1)=zmpi1(2,I1-n1/2,J2,j3,Jp2)
        end do
     end do
     J2stb=1
  end do

end subroutine mpiswitch_upcorn


subroutine unswitch_downcorn(nfft,n2,lot,n1,lzt,zw,zt)
  implicit real(kind=8) (a-h,o-z)
  dimension zw(2,lot,n2),zt(2,lzt,n1)
! WARNING: Assuming that high frequencies are in the corners 
!          and that n2 is multiple of 2

! Low frequencies
  do j=1,nfft
     do i=1,n2/2
        zt(1,i,j)=zw(1,j,i)
        zt(2,i,j)=zw(2,j,i)
     end do
  end do

end subroutine unswitch_downcorn


subroutine unmpiswitch_downcorn(j3,nfft,Jp2stf,J2stf,lot,n1,md2,nd3,nproc,zw,zmpi1)
  implicit real(kind=8) (a-h,o-z)
  dimension zmpi1(2,n1/2,md2/nproc,nd3/nproc,nproc),zw(2,lot,n1)
! WARNING: Assuming that high frequencies are in the corners 
!          and that n1 is multiple of 2

  mfft=0
  do Jp2=Jp2stf,nproc
     do J2=J2stf,md2/nproc
        mfft=mfft+1
        if (mfft.gt.nfft) then
           Jp2stf=Jp2
           J2stf=J2
           return
        endif
        do I1=1,n1/2
           zmpi1(1,I1,J2,j3,Jp2)=zw(1,mfft,I1)
           zmpi1(2,I1,J2,j3,Jp2)=zw(2,mfft,I1)
        end do
     end do
     J2stf=1
  end do

end subroutine unmpiswitch_downcorn


!!****f* BigDFT/F_unfill_downcorn
!! NAME
!!   F_unfill_downcorn
!!
!! FUNCTION
!!     (Based on suitable modifications of S.Goedecker routines)
!!     Restore data into output array, calculating in the meanwhile
!!     Hartree energy of the potential 
!!
!! SYNOPSIS
!!     zf:          Original distributed density as well as
!!                  Distributed solution of the poisson equation (inout)
!!     zw:          FFT work array
!!     n3:          (twice the) dimension of the last FFTtransform.
!!     md1,md3:     Dimensions of the undistributed part of the real grid
!!     nfft:        number of planes
!!     scal:        Needed to achieve unitarity and correct dimensions
!!     ehartreetmp: Hartree energy
!!
!! WARNING
!!     Assuming that high frequencies are in the corners 
!!     and that n3 is multiple of 4   
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
subroutine F_unfill_downcorn(md1,md3,lot,nfft,n3,zw,zf&
     ,scal,ehartreetmp)
  implicit none
  !Arguments
  integer, intent(in) :: md1,md3,lot,nfft,n3
  real(kind=8), dimension(2,lot,n3/2), intent(in) :: zw
  real(kind=8), dimension(md1,md3),intent(inout) :: zf
  real(kind=8), intent(in) :: scal
  real(kind=8), intent(out) :: ehartreetmp
  !Local variables
  integer :: i3,i1
  real(kind=8) :: pot1

  !Body
  ehartreetmp=0.d0
  do i3=1,n3/4
     do i1=1,nfft
        pot1 = scal*zw(1,i1,i3)
        ehartreetmp =ehartreetmp + pot1* zf(i1,2*i3-1)
        zf(i1,2*i3-1)= pot1 
        pot1 = scal*zw(2,i1,i3)
        ehartreetmp =ehartreetmp + pot1* zf(i1,2*i3)
        zf(i1,2*i3)= pot1 
     enddo
  end do
  
end subroutine F_unfill_downcorn
!!**
