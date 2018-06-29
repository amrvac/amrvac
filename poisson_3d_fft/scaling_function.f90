
!! Copyright (C) 2002-2007 BigDFT group 
!! This file is distributed under the terms of the
!! GNU General Public License, see ~/COPYING file
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the list of contributors, see ~/AUTHORS 


!!****h* BigDFT/scaling_function
!! NAME
!!   scaling_function
!!
!! FUNCTION
!!   Calculate the values of a scaling function in real uniform grid
!!
!! SOURCE
!!
subroutine scaling_function(itype,nd,nrange,a,x)

  implicit none
  !Arguments
  !Type of interpolating functions
  integer, intent(in) :: itype
  !Number of points: must be 2**nex
  integer, intent(in) :: nd
  integer, intent(out) :: nrange
  real(kind=8), dimension(0:nd), intent(out) :: a,x
  !Local variables
  real(kind=8), dimension(:), allocatable :: y
  integer :: i,nt,ni,i_all,i_stat
  
  !Only itype=8,14,16,20,24,30,40,50,60,100
  select case(itype)
  case(8,14,16,20,24,30,40,50,60,100)
     !O.K.
  case default
     print *,"Only interpolating functions 8, 14, 16, 20, 24, 30, 40, 50, 60, 100"
     stop
  end select
!!$  write(unit=*,fmt="(1x,a,i0,a)") &
!!$       "Use interpolating scaling functions of ",itype," order"

  !Give the range of the scaling function
  !from -itype to itype
  ni=2*itype
  nrange = ni
  allocate(y(0:nd),stat=i_stat)
  call memocc(i_stat,product(shape(y))*kind(y),'y','scaling_function')
  
  ! plot scaling function
  call zero(nd+1,x)
  call zero(nd+1,y)
  nt=ni
  x(nt/2-1)=1.d0
  loop1: do
     nt=2*nt
     ! write(6,*) 'nd,nt',nd,nt
     select case(itype)
     case(8)
        call back_trans_8(nd,nt,x,y)
     case(14)
        call back_trans_14(nd,nt,x,y)
     case(16)
        call back_trans_16(nd,nt,x,y)
     case(20)
        call back_trans_20(nd,nt,x,y)
     case(24)
        call back_trans_24(nd,nt,x,y)
     case(30)
        call back_trans_30(nd,nt,x,y)
     case(40)
        call back_trans_40(nd,nt,x,y)
     case(50)
        call back_trans_50(nd,nt,x,y)
     case(60)
        call back_trans_60(nd,nt,x,y)
     case(100)
        call back_trans_100(nd,nt,x,y)
     end select
     call dcopy(nt,y,1,x,1)
     if (nt.eq.nd) then
        exit loop1
     end if
  end do loop1

  !open (unit=1,file='scfunction',status='unknown')
  do i=0,nd
     a(i) = real(i*ni,kind=8)/real(nd,kind=8)-(.5d0*real(ni,kind=8)-1.d0)
     !write(1,*) a(i),x(i)
  end do
  !close(1)

  i_all=-product(shape(y))*kind(y)
  deallocate(y,stat=i_stat)
  call memocc(i_stat,i_all,'y','scaling_function')
end subroutine scaling_function
!!***


!!****h* BigDFT/wavelet_function
!! NAME
!!   wavelet_function
!!
!! FUNCTION
!!   Calculate the values of the wavelet function in a real uniform mesh.
!!
!! SOURCE
!!
subroutine wavelet_function(itype,nd,a,x)

  implicit none
  !Arguments
  !Type of the interpolating scaling function
  integer, intent(in) :: itype
  !must be 2**nex
  integer, intent(in) :: nd
  real(kind=8), dimension(0:nd), intent(out) :: a,x
  !Local variables
  real(kind=8), dimension(:), allocatable :: y
  integer :: i,nt,ni,i_all,i_stat

  !Only itype=8,14,16,20,24,30,40,50,60,100
  Select case(itype)
  case(8,14,16,20,24,30,40,50,60,100)
     !O.K.
  case default
     print *,"Only interpolating functions 8, 14, 16, 20, 24, 30, 40, 50, 60, 100"
     stop
  end select

  !Give the range of the scaling function
  !from -itype to itype
  ni=2*itype

  allocate(y(0:nd),stat=i_stat)
  call memocc(i_stat,product(shape(y))*kind(y),'y','wavelet_function')
  
  ! plot wavelet 
  call zero(nd+1,x)
  call zero(nd+1,y)
  nt=ni
  x(nt+nt/2-1)=1.d0
  loop3: do
     nt=2*nt
     !write(6,*) 'nd,nt',nd,nt
     select case(itype)
     case(8)
        call back_trans_8(nd,nt,x,y)
     case(14)
        call back_trans_14(nd,nt,x,y)
     case(16)
        call back_trans_16(nd,nt,x,y)
     case(20)
        call back_trans_20(nd,nt,x,y)
     case(24)
        call back_trans_24(nd,nt,x,y)
     case(30)
        call back_trans_30(nd,nt,x,y)
     case(40)
        call back_trans_40(nd,nt,x,y)
     case(50)
        call back_trans_50(nd,nt,x,y)
     case(60)
        call back_trans_60(nd,nt,x,y)
     case(100)
        call back_trans_100(nd,nt,x,y)
     end select
     call dcopy(nd,y,1,x,1)
     if (nt.eq.nd) then
        exit loop3
     end if
  end do loop3

  !open (unit=1,file='wavelet',status='unknown')
  do i=0,nd-1
     a(i) = real(i*ni,kind=8)/real(nd,kind=8)-(.5d0*real(ni,kind=8)-.5d0)
     !write(1,*) a(i),x(i)
  end do
  !close(1)

  i_all=-product(shape(y))*kind(y)
  deallocate(y,stat=i_stat)
  call memocc(i_stat,i_all,'y','wavelet_function')
 
end subroutine wavelet_function
!!***


!!****h* BigDFT/scf_recursion
!! NAME
!!   scf_recursion
!!
!! FUNCTION
!!   Do iterations to go from p0gauss to pgauss
!!   order interpolating scaling function
!!
!! SOURCE
!!
subroutine scf_recursion(itype,n_iter,n_range,kernel_scf,kern_1_scf)
  implicit none
  !Arguments
  integer, intent(in) :: itype,n_iter,n_range
  real(kind=8), intent(inout) :: kernel_scf(-n_range:n_range)
  real(kind=8), intent(out) :: kern_1_scf(-n_range:n_range)
  !Local variables

  !Only itype=8,14,16,20,24,30,40,50,60,100
  select case(itype)
  case(8,14,16,20,24,30,40,50,60,100)
     !O.K.
  case default
     print *,"Only interpolating functions 8, 14, 16, 20, 24, 30, 40, 50, 60, 100"
     stop
  end select

  select case(itype)
  case(8)
     call scf_recursion_8(n_iter,n_range,kernel_scf,kern_1_scf)
  case(14)
     call scf_recursion_14(n_iter,n_range,kernel_scf,kern_1_scf)
  case(16)
     call scf_recursion_16(n_iter,n_range,kernel_scf,kern_1_scf)
  case(20)
     call scf_recursion_20(n_iter,n_range,kernel_scf,kern_1_scf)
  case(24)
     call scf_recursion_24(n_iter,n_range,kernel_scf,kern_1_scf)
  case(30)
     call scf_recursion_30(n_iter,n_range,kernel_scf,kern_1_scf)
  case(40)
     call scf_recursion_40(n_iter,n_range,kernel_scf,kern_1_scf)
  case(50)
     call scf_recursion_50(n_iter,n_range,kernel_scf,kern_1_scf)
  case(60)
     call scf_recursion_60(n_iter,n_range,kernel_scf,kern_1_scf)
  case(100)
     call scf_recursion_100(n_iter,n_range,kernel_scf,kern_1_scf)
  end select

end subroutine scf_recursion
!!***


!!****h* BigDFT/zero
!! NAME
!!   zero
!!
!! FUNCTION
!!   Set to zero an array x(n)
!!
!! SOURCE
!!
subroutine zero(n,x)
  implicit none
  !Arguments
  integer, intent(in) :: n
  real(kind=8), intent(out) :: x(n)
  !Local variables
  integer :: i
  do i=1,n
     x(i)=0.d0
  end do
end subroutine zero
!!***


!!****h* BigDFT/for_trans_8
!! NAME
!!   for_trans_8
!!
!! FUNCTION
!!   forward wavelet transform
!!   nd: length of data set
!!   nt length of data in data set to be transformed
!!   m filter length (m has to be even!)
!!   x input data, y output data
!!
!! SOURCE
!!
subroutine for_trans_8(nd,nt,x,y)
  implicit none
  !Arguments
  integer, intent(in) :: nd,nt
  real(kind=8), intent(in) :: x(0:nd-1)
  real(kind=8), intent(out) :: y(0:nd-1)
  !Local variables
  integer :: i,j,ind

  include 'lazy_8.inc'

  do i=0,nt/2-1
     y(     i)=0.d0
     y(nt/2+i)=0.d0
     
     do j=-m+1,m
        
        ! periodically wrap index if necessary
        ind=j+2*i
        loop99: do
           if (ind.lt.0) then 
              ind=ind+nt
              cycle loop99
           end if
           if (ind.ge.nt) then 
              ind=ind-nt
              cycle loop99
           end if
           exit loop99
        end do loop99

        y(     i)=y(     i)+cht(j)*x(ind)
        y(nt/2+i)=y(nt/2+i)+cgt(j)*x(ind)
     end do
     
  end do

end subroutine for_trans_8
!!***


!!****h* BigDFT/back_trans_8
!! NAME
!!   back_trans_8
!!
!! FUNCTION
!!
!! SOURCE
!!
subroutine back_trans_8(nd,nt,x,y)
  ! backward wavelet transform
  ! nd: length of data set
  ! nt length of data in data set to be transformed
  ! m filter length (m has to be even!)
  ! x input data, y output data
  implicit none
  !Arguments
  integer, intent(in) :: nd,nt
  real(kind=8), intent(in) :: x(0:nd-1)
  real(kind=8), intent(out) :: y(0:nd-1)
  !Local variables
  integer :: i,j,ind

  include 'lazy_8.inc'
  
  do i=0,nt/2-1
     y(2*i+0)=0.d0
     y(2*i+1)=0.d0
     
     do j=-m/2,m/2-1
        
        ! periodically wrap index if necessary
        ind=i-j
        loop99: do
           if (ind.lt.0) then 
              ind=ind+nt/2
              cycle loop99
           end if
           if (ind.ge.nt/2) then 
              ind=ind-nt/2
              cycle loop99
           end if
           exit loop99
        end do loop99

        y(2*i+0)=y(2*i+0) + ch(2*j-0)*x(ind)+cg(2*j-0)*x(ind+nt/2)
        y(2*i+1)=y(2*i+1) + ch(2*j+1)*x(ind)+cg(2*j+1)*x(ind+nt/2)
     end do
        
  end do
        
end subroutine back_trans_8
!!***


!!****h* BigDFT/ftest_8
!! NAME
!!   ftest_8
!!
!! FUNCTION
!!   Tests the 4 orthogonality relations of the filters
!!
!! SOURCE
!!
subroutine ftest_8
  implicit none
  !Arguments
  !Local variables
  character(len=*), parameter :: fmt22 = "(a,i3,i4,4(e17.10))"
  integer :: i,j,l
  real(kind=8) :: t1,t2,t3,t4,eps

  include 'lazy_8.inc'
  
  ! do i=-m,m
  ! write(6,*) i,ch(i),cg(i)
  ! end do
  
  do i=-m,m
     do j=-m,m
        t1=0.d0
        t2=0.d0
        t3=0.d0
        t4=0.d0
        do l=-3*m,3*m
           if ( l-2*i.ge.-m .and. l-2*i.le.m  .and. &
                l-2*j.ge.-m .and. l-2*j.le.m ) then
              t1=t1+ch(l-2*i)*cht(l-2*j)
              t2=t2+cg(l-2*i)*cgt(l-2*j)
              t3=t3+ch(l-2*i)*cgt(l-2*j)
              t4=t4+cht(l-2*i)*cg(l-2*j)
           end if
        end do
        eps=1.d-10
        if (i.eq.j) then
           if (abs(t1-1.d0).gt.eps .or. abs(t2-1.d0).gt.eps .or. &
             & abs(t3).gt.eps  .or. abs(t4).gt.eps ) then 
              write(6,fmt22) 'Orthogonality ERROR', i,j,t1,t2,t3,t4
           end if
        else
           if (abs(t1).gt.eps .or. abs(t2).gt.eps .or. &
             & abs(t3).gt.eps  .or. abs(t4).gt.eps ) then
              write(6,fmt22) 'Orthogonality ERROR', i,j,t1,t2,t3,t4
           end if
        end if
     end do
  end do
  
  write(6,*) 'FILTER TEST PASSED'
  
end subroutine ftest_8
!!***


!!****h* BigDFT/scf_recursion_8
!! NAME
!!   scf_recursion_8
!!
!! FUNCTION
!!   Do iterations to go from p0gauss to pgauss
!!   8th-order interpolating scaling function
!!
!! SOURCE
!!
subroutine scf_recursion_8(n_iter,n_range,kernel_scf,kern_1_scf)
  implicit none
  !Arguments
  integer, intent(in) :: n_iter,n_range
  real(kind=8), intent(inout) :: kernel_scf(-n_range:n_range)
  real(kind=8), intent(out) :: kern_1_scf(-n_range:n_range)
  !Local variables
  real(kind=8) :: kern,kern_tot
  integer :: i_iter,i,j,ind

  include "lazy_8.inc"

  !Start the iteration to go from p0gauss to pgauss
  loop_iter_scf: do i_iter=1,n_iter
     kern_1_scf(:) = kernel_scf(:)
     kernel_scf(:) = 0.d0
     loop_iter_i: do i=0,n_range
        kern_tot = 0.d0
        do j=-m,m
           ind = 2*i-j
           if (abs(ind) > n_range) then
              kern = 0.d0
           else
              kern = kern_1_scf(ind)
           end if
           kern_tot = kern_tot + ch(j)*kern
        end do
        if (kern_tot == 0.d0) then
           !zero after (be sure because strictly == 0.d0)
           exit loop_iter_i
        else
           kernel_scf( i) = 0.5d0*kern_tot
           kernel_scf(-i) = kernel_scf(i)
        end if
     end do loop_iter_i
  end do loop_iter_scf
end subroutine scf_recursion_8
!!***


!!****h* BigDFT/for_trans_14
!! NAME
!!   for_trans_14
!!
!! FUNCTION
!!   forward wavelet transform
!!   nd: length of data set
!!   nt length of data in data set to be transformed
!!   m filter length (m has to be even!)
!!   x input data, y output data
!!
!! SOURCE
!!
subroutine for_trans_14(nd,nt,x,y)
  implicit none
  !Arguments
  integer, intent(in) :: nd,nt
  real(kind=8), intent(in) :: x(0:nd-1)
  real(kind=8), intent(out) :: y(0:nd-1)
  !Local variables
  integer :: i,j,ind

  include 'lazy_14.inc'

  do i=0,nt/2-1
     y(     i)=0.d0
     y(nt/2+i)=0.d0
     
     do j=-m+1,m
        
        ! periodically wrap index if necessary
        ind=j+2*i
        loop99: do
           if (ind.lt.0) then 
              ind=ind+nt
              cycle loop99
           end if
           if (ind.ge.nt) then 
              ind=ind-nt
              cycle loop99
           end if
           exit loop99
        end do loop99

        y(     i)=y(     i)+cht(j)*x(ind)
        y(nt/2+i)=y(nt/2+i)+cgt(j)*x(ind)
     end do
  end do
  
end subroutine for_trans_14
!!***


!!****h* BigDFT/back_trans_14
!! NAME
!!   back_trans_14
!!
!! FUNCTION
!!   backward wavelet transform
!!   nd: length of data set
!!   nt length of data in data set to be transformed
!!   m filter length (m has to be even!)
!!   x input data, y output data
!!
!! SOURCE
!!
subroutine back_trans_14(nd,nt,x,y)
  implicit none
  !Arguments
  integer, intent(in) :: nd,nt
  real(kind=8), intent(in) :: x(0:nd-1)
  real(kind=8), intent(out) :: y(0:nd-1)
  !Local variables
  integer :: i,j,ind

  include 'lazy_14.inc'
  
  do i=0,nt/2-1
     y(2*i+0)=0.d0
     y(2*i+1)=0.d0
     
     do j=-m/2,m/2-1
        
        ! periodically wrap index if necessary
        ind=i-j
        loop99: do
           if (ind.lt.0) then 
              ind=ind+nt/2
              cycle loop99
           end if
           if (ind.ge.nt/2) then 
              ind=ind-nt/2
              cycle loop99
           end if
           exit loop99
        end do loop99

        y(2*i+0)=y(2*i+0) + ch(2*j-0)*x(ind)+cg(2*j-0)*x(ind+nt/2)
        y(2*i+1)=y(2*i+1) + ch(2*j+1)*x(ind)+cg(2*j+1)*x(ind+nt/2)

     end do
  end do
        
end subroutine back_trans_14
!!***


!!****h* BigDFT/ftest_14
!! NAME
!!   ftest_14
!!
!! FUNCTION
!!   Tests the 4 orthogonality relations of the filters
!!
!! SOURCE
!!
subroutine ftest_14
  implicit none
  !Arguments
  !Local variables
  character(len=*), parameter :: fmt22 = "(a,i3,i4,4(e17.10))"
  real(kind=8) :: t1,t2,t3,t4,eps
  integer :: i,j,l

  include 'lazy_14.inc'
  
  ! do i=-m,m
  ! write(6,*) i,ch(i),cg(i)
  ! end do
  
  do i=-m,m
     do j=-m,m
        t1=0.d0
        t2=0.d0
        t3=0.d0
        t4=0.d0
        do l=-3*m,3*m
           if ( l-2*i.ge.-m .and. l-2*i.le.m  .and. &
                l-2*j.ge.-m .and. l-2*j.le.m ) then
              t1=t1+ch(l-2*i)*cht(l-2*j)
              t2=t2+cg(l-2*i)*cgt(l-2*j)
              t3=t3+ch(l-2*i)*cgt(l-2*j)
              t4=t4+cht(l-2*i)*cg(l-2*j)
           end if
        end do
        eps=1.d-10
        if (i.eq.j) then
           if (abs(t1-1.d0).gt.eps .or. abs(t2-1.d0).gt.eps .or. &
             & abs(t3).gt.eps  .or. abs(t4).gt.eps ) then 
              write(6,fmt22) 'Orthogonality ERROR', i,j,t1,t2,t3,t4
           end if
        else
           if (abs(t1).gt.eps .or. abs(t2).gt.eps .or. &
             & abs(t3).gt.eps  .or. abs(t4).gt.eps ) then
              write(6,fmt22) 'Orthogonality ERROR', i,j,t1,t2,t3,t4
           end if
        end if
     end do
  end do
  
  write(6,*) 'FILTER TEST PASSED'
  
end subroutine ftest_14
!!***


!!****h* BigDFT/scf_recursion_14
!! NAME
!!   scf_recursion_14
!!
!! FUNCTION
!!   Do iterations to go from p0gauss to pgauss
!!   14th-order interpolating scaling function
!!
!! SOURCE
!!
subroutine scf_recursion_14(n_iter,n_range,kernel_scf,kern_1_scf)
  implicit none
  !Arguments
  integer, intent(in) :: n_iter,n_range
  real(kind=8), intent(inout) :: kernel_scf(-n_range:n_range)
  real(kind=8), intent(out) :: kern_1_scf(-n_range:n_range)
  !Local variables
  real(kind=8) :: kern,kern_tot
  integer :: i_iter,i,j,ind

  include "lazy_14.inc"

  !Start the iteration to go from p0gauss to pgauss
  loop_iter_scf: do i_iter=1,n_iter
     kern_1_scf(:) = kernel_scf(:)
     kernel_scf(:) = 0.d0
     loop_iter_i: do i=0,n_range
        kern_tot = 0.d0
        do j=-m,m
           ind = 2*i-j
           if (abs(ind) > n_range) then
              kern = 0.d0
           else
              kern = kern_1_scf(ind)
           end if
           kern_tot = kern_tot + ch(j)*kern
        end do
        if (kern_tot == 0.d0) then
           !zero after (be sure because strictly == 0.d0)
           exit loop_iter_i
        else
           kernel_scf( i) = 0.5d0*kern_tot
           kernel_scf(-i) = kernel_scf(i)
        end if
     end do loop_iter_i
  end do loop_iter_scf
end subroutine scf_recursion_14
!!***


!!****h* BigDFT/for_trans_16
!! NAME
!!   for_trans_16
!!
!! FUNCTION
!!   forward wavelet transform
!!   nd: length of data set
!!   nt length of data in data set to be transformed
!!   m filter length (m has to be even!)
!!   x input data, y output data
!!
!! SOURCE
!!
subroutine for_trans_16(nd,nt,x,y)
  implicit none
  !Arguments
  integer, intent(in) :: nd,nt
  real(kind=8), intent(in) :: x(0:nd-1)
  real(kind=8), intent(out) :: y(0:nd-1)
  !Local variables
  integer :: i,j,ind

  include 'lazy_16.inc'

  do i=0,nt/2-1
     y(     i)=0.d0
     y(nt/2+i)=0.d0
     
     do j=-m+1,m
        
        ! periodically wrap index if necessary
        ind=j+2*i
        loop99: do
           if (ind.lt.0) then 
              ind=ind+nt
              cycle loop99
           end if
           if (ind.ge.nt) then 
              ind=ind-nt
              cycle loop99
           end if
           exit loop99
        end do loop99

        y(     i)=y(     i)+cht(j)*x(ind)
        y(nt/2+i)=y(nt/2+i)+cgt(j)*x(ind)
     end do
  end do
  
end subroutine for_trans_16
!!***


!!****h* BigDFT/back_trans_16
!! NAME
!!   back_trans_16
!!
!! FUNCTION
!!   backward wavelet transform
!!   nd: length of data set
!!   nt length of data in data set to be transformed
!!   m filter length (m has to be even!)
!!   x input data, y output data
!!
!! SOURCE
!!
subroutine back_trans_16(nd,nt,x,y)
  implicit none
  !Arguments
  integer, intent(in) :: nd,nt
  real(kind=8), intent(in) :: x(0:nd-1)
  real(kind=8), intent(out) :: y(0:nd-1)
  !Local variables
  integer :: i,j,ind

  include 'lazy_16.inc'
  
  do i=0,nt/2-1
     y(2*i+0)=0.d0
     y(2*i+1)=0.d0
     
     do j=-m/2,m/2-1
        
        ! periodically wrap index if necessary
        ind=i-j
        loop99: do
           if (ind.lt.0) then 
              ind=ind+nt/2
              cycle loop99
           end if
           if (ind.ge.nt/2) then 
              ind=ind-nt/2
              cycle loop99
           end if
           exit loop99
        end do loop99

        y(2*i+0)=y(2*i+0) + ch(2*j-0)*x(ind)+cg(2*j-0)*x(ind+nt/2)
        y(2*i+1)=y(2*i+1) + ch(2*j+1)*x(ind)+cg(2*j+1)*x(ind+nt/2)
     end do
  end do
        
end subroutine back_trans_16
!!***


!!****h* BigDFT/ftest_16
!! NAME
!!   ftest_16
!!
!! FUNCTION
!!   Tests the 4 orthogonality relations of the filters
!!
!! SOURCE
!!
subroutine ftest_16
  implicit none
  !Arguments
  !Local variables
  character(len=*), parameter :: fmt22 = "(a,i3,i4,4(e17.10))"
  real(kind=8) :: t1,t2,t3,t4,eps
  integer :: i,j,l

  include 'lazy_16.inc'
  
  ! do i=-m,m
  ! write(6,*) i,ch(i),cg(i)
  ! end do
  
  do i=-m,m
     do j=-m,m
        t1=0.d0
        t2=0.d0
        t3=0.d0
        t4=0.d0
        do l=-3*m,3*m
           if ( l-2*i.ge.-m .and. l-2*i.le.m  .and. &
                l-2*j.ge.-m .and. l-2*j.le.m ) then
              t1=t1+ch(l-2*i)*cht(l-2*j)
              t2=t2+cg(l-2*i)*cgt(l-2*j)
              t3=t3+ch(l-2*i)*cgt(l-2*j)
              t4=t4+cht(l-2*i)*cg(l-2*j)
           end if
        end do
        eps=1.d-10
        if (i.eq.j) then
           if (abs(t1-1.d0).gt.eps .or. abs(t2-1.d0).gt.eps .or. &
             & abs(t3).gt.eps  .or. abs(t4).gt.eps ) then 
              write(6,fmt22) 'Orthogonality ERROR', i,j,t1,t2,t3,t4
           end if
        else
           if (abs(t1).gt.eps .or. abs(t2).gt.eps .or. &
             & abs(t3).gt.eps  .or. abs(t4).gt.eps ) then
              write(6,fmt22) 'Orthogonality ERROR', i,j,t1,t2,t3,t4
           end if
        end if
     end do
  end do
  
  write(6,*) 'FILTER TEST PASSED'
  
end subroutine ftest_16
!!***


!!****h* BigDFT/scf_recursion_16
!! NAME
!!   scf_recursion_16
!!
!! FUNCTION
!!   Do iterations to go from p0gauss to pgauss
!!   16th-order interpolating scaling function
!!
!! SOURCE
!!
subroutine scf_recursion_16(n_iter,n_range,kernel_scf,kern_1_scf)
  implicit none
  !Arguments
  integer, intent(in) :: n_iter,n_range
  real(kind=8), intent(inout) :: kernel_scf(-n_range:n_range)
  real(kind=8), intent(out) :: kern_1_scf(-n_range:n_range)
  !Local variables
  real(kind=8) :: kern,kern_tot
  integer :: i_iter,i,j,ind

  include "lazy_16.inc"

  !Start the iteration to go from p0gauss to pgauss
  loop_iter_scf: do i_iter=1,n_iter
     kern_1_scf(:) = kernel_scf(:)
     kernel_scf(:) = 0.d0
     loop_iter_i: do i=0,n_range
        kern_tot = 0.d0
        do j=-m,m
           ind = 2*i-j
           if (abs(ind) > n_range) then
              kern = 0.d0
           else
              kern = kern_1_scf(ind)
           end if
           kern_tot = kern_tot + ch(j)*kern
        end do
        if (kern_tot == 0.d0) then
           !zero after (be sure because strictly == 0.d0)
           exit loop_iter_i
        else
           kernel_scf( i) = 0.5d0*kern_tot
           kernel_scf(-i) = kernel_scf(i)
        end if
     end do loop_iter_i
  end do loop_iter_scf
end subroutine scf_recursion_16
!!***


!!****h* BigDFT/for_trans_20
!! NAME
!!   for_trans_20
!!
!! FUNCTION
!!   forward wavelet transform
!!   nd: length of data set
!!   nt length of data in data set to be transformed
!!   m filter length (m has to be even!)
!!   x input data, y output data
!!
!! SOURCE
!!
subroutine for_trans_20(nd,nt,x,y)
  implicit none
  !Arguments
  integer, intent(in) :: nd,nt
  real(kind=8), intent(in) :: x(0:nd-1)
  real(kind=8), intent(out) :: y(0:nd-1)
  !Local variables
  integer :: i,j,ind

  include 'lazy_20.inc'

  do i=0,nt/2-1
     y(     i)=0.d0
     y(nt/2+i)=0.d0
     
     do j=-m+1,m
        
        ! periodically wrap index if necessary
        ind=j+2*i
        loop99: do
           if (ind.lt.0) then 
              ind=ind+nt
              cycle loop99
           end if
           if (ind.ge.nt) then 
              ind=ind-nt
              cycle loop99
           end if
           exit loop99
        end do loop99

        y(     i)=y(     i)+cht(j)*x(ind)
        y(nt/2+i)=y(nt/2+i)+cgt(j)*x(ind)
     end do
  end do
  
end subroutine for_trans_20
!!***


!!****h* BigDFT/back_trans_20
!! NAME
!!   back_trans_20
!!
!! FUNCTION
!!   backward wavelet transform
!!   nd: length of data set
!!   nt length of data in data set to be transformed
!!   m filter length (m has to be even!)
!!   x input data, y output data
!!
!! SOURCE
!!
subroutine back_trans_20(nd,nt,x,y)
  implicit none
  !Arguments
  integer, intent(in) :: nd,nt
  real(kind=8), intent(in) :: x(0:nd-1)
  real(kind=8), intent(out) :: y(0:nd-1)
  !Local variables
  integer :: i,j,ind

  include 'lazy_20.inc'
  
  do i=0,nt/2-1
     y(2*i+0)=0.d0
     y(2*i+1)=0.d0
     
     do j=-m/2,m/2-1
        
        ! periodically wrap index if necessary
        ind=i-j
        loop99: do
           if (ind.lt.0) then 
              ind=ind+nt/2
              cycle loop99
           end if
           if (ind.ge.nt/2) then 
              ind=ind-nt/2
              cycle loop99
           end if
           exit loop99
        end do loop99

        y(2*i+0)=y(2*i+0) + ch(2*j-0)*x(ind)+cg(2*j-0)*x(ind+nt/2)
        y(2*i+1)=y(2*i+1) + ch(2*j+1)*x(ind)+cg(2*j+1)*x(ind+nt/2)
     end do
  end do
        
end subroutine back_trans_20
!!***


!!****h* BigDFT/ftest_20
!! NAME
!!   ftest_20
!!
!! FUNCTION
!!   Tests the 4 orthogonality relations of the filters
!!
!! SOURCE
!!
subroutine ftest_20
  implicit none
  !Arguments
  !Local variables
  character(len=*), parameter :: fmt22 = "(a,i3,i4,4(e17.10))"
  real(kind=8) :: t1,t2,t3,t4,eps
  integer :: i,j,l

  include 'lazy_20.inc'
  
  ! do i=-m,m
  ! write(6,*) i,ch(i),cg(i)
  ! end do
  
  do i=-m,m
     do j=-m,m
        t1=0.d0
        t2=0.d0
        t3=0.d0
        t4=0.d0
        do l=-3*m,3*m
           if ( l-2*i.ge.-m .and. l-2*i.le.m  .and. &
                l-2*j.ge.-m .and. l-2*j.le.m ) then
              t1=t1+ch(l-2*i)*cht(l-2*j)
              t2=t2+cg(l-2*i)*cgt(l-2*j)
              t3=t3+ch(l-2*i)*cgt(l-2*j)
              t4=t4+cht(l-2*i)*cg(l-2*j)
           end if
        end do
        eps=1.d-10
        if (i.eq.j) then
           if (abs(t1-1.d0).gt.eps .or. abs(t2-1.d0).gt.eps .or. &
             & abs(t3).gt.eps  .or. abs(t4).gt.eps ) then 
              write(6,fmt22) 'Orthogonality ERROR', i,j,t1,t2,t3,t4
           end if
        else
           if (abs(t1).gt.eps .or. abs(t2).gt.eps .or. &
             & abs(t3).gt.eps  .or. abs(t4).gt.eps ) then
              write(6,fmt22) 'Orthogonality ERROR', i,j,t1,t2,t3,t4
           end if
        end if
     end do
  end do
  
  write(6,*) 'FILTER TEST PASSED'
  
end subroutine ftest_20
!!***


!!****h* BigDFT/scf_recursion_20
!! NAME
!!   scf_recursion_20
!!
!! FUNCTION
!!   Do iterations to go from p0gauss to pgauss
!!   16th-order interpolating scaling function
!!
!! SOURCE
!!
subroutine scf_recursion_20(n_iter,n_range,kernel_scf,kern_1_scf)
  implicit none
  !Arguments
  integer, intent(in) :: n_iter,n_range
  real(kind=8), intent(inout) :: kernel_scf(-n_range:n_range)
  real(kind=8), intent(out) :: kern_1_scf(-n_range:n_range)
  !Local variables
  real(kind=8) :: kern,kern_tot
  integer :: i_iter,i,j,ind

  include "lazy_20.inc"

  !Start the iteration to go from p0gauss to pgauss
  loop_iter_scf: do i_iter=1,n_iter
     kern_1_scf(:) = kernel_scf(:)
     kernel_scf(:) = 0.d0
     loop_iter_i: do i=0,n_range
        kern_tot = 0.d0
        do j=-m,m
           ind = 2*i-j
           if (abs(ind) > n_range) then
              kern = 0.d0
           else
              kern = kern_1_scf(ind)
           end if
           kern_tot = kern_tot + ch(j)*kern
        end do
        if (kern_tot == 0.d0) then
           !zero after (be sure because strictly == 0.d0)
           exit loop_iter_i
        else
           kernel_scf( i) = 0.5d0*kern_tot
           kernel_scf(-i) = kernel_scf(i)
        end if
     end do loop_iter_i
  end do loop_iter_scf
end subroutine scf_recursion_20
!!***


!!****h* BigDFT/for_trans_24
!! NAME
!!   for_trans_24
!!
!! FUNCTION
!!   forward wavelet transform
!!   nd: length of data set
!!   nt length of data in data set to be transformed
!!   m filter length (m has to be even!)
!!   x input data, y output data
!!
!! SOURCE
!!
subroutine for_trans_24(nd,nt,x,y)
  implicit none
  !Arguments
  integer, intent(in) :: nd,nt
  real(kind=8), intent(in) :: x(0:nd-1)
  real(kind=8), intent(out) :: y(0:nd-1)
  !Local variables
  integer :: i,j,ind

  include 'lazy_24.inc'

  do i=0,nt/2-1
     y(     i)=0.d0
     y(nt/2+i)=0.d0
     
     do j=-m+1,m
        
        ! periodically wrap index if necessary
        ind=j+2*i
        loop99: do
           if (ind.lt.0) then 
              ind=ind+nt
              cycle loop99
           end if
           if (ind.ge.nt) then 
              ind=ind-nt
              cycle loop99
           end if
           exit loop99
        end do loop99

        y(     i)=y(     i)+cht(j)*x(ind)
        y(nt/2+i)=y(nt/2+i)+cgt(j)*x(ind)
     end do
  end do
  
end subroutine for_trans_24
!!***


!!****h* BigDFT/back_trans_24
!! NAME
!!   back_trans_24
!!
!! FUNCTION
!!   backward wavelet transform
!!   nd: length of data set
!!   nt length of data in data set to be transformed
!!   m filter length (m has to be even!)
!!   x input data, y output data
!!
!! SOURCE
!!
subroutine back_trans_24(nd,nt,x,y)
  implicit none
  !Arguments
  integer, intent(in) :: nd,nt
  real(kind=8), intent(in) :: x(0:nd-1)
  real(kind=8), intent(out) :: y(0:nd-1)
  !Local variables
  integer :: i,j,ind

  include 'lazy_24.inc'
  
  do i=0,nt/2-1
     y(2*i+0)=0.d0
     y(2*i+1)=0.d0
     
     do j=-m/2,m/2-1
        
        ! periodically wrap index if necessary
        ind=i-j
        loop99: do
           if (ind.lt.0) then 
              ind=ind+nt/2
              cycle loop99
           end if
           if (ind.ge.nt/2) then 
              ind=ind-nt/2
              cycle loop99
           end if
           exit loop99
        end do loop99

        y(2*i+0)=y(2*i+0) + ch(2*j-0)*x(ind)+cg(2*j-0)*x(ind+nt/2)
        y(2*i+1)=y(2*i+1) + ch(2*j+1)*x(ind)+cg(2*j+1)*x(ind+nt/2)
     end do
  end do
        
end subroutine back_trans_24
!!***


!!****h* BigDFT/ftest_24
!! NAME
!!   ftest_24
!!
!! FUNCTION
!!   Tests the 4 orthogonality relations of the filters
!!
!! SOURCE
!!
subroutine ftest_24
  implicit none
  !Arguments
  !Local variables
  character(len=*), parameter :: fmt22 = "(a,i3,i4,4(e17.10))"
  real(kind=8) :: t1,t2,t3,t4,eps
  integer :: i,j,l

  include 'lazy_24.inc'
  
  ! do i=-m,m
  ! write(6,*) i,ch(i),cg(i)
  ! end do
  
  do i=-m,m
     do j=-m,m
        t1=0.d0
        t2=0.d0
        t3=0.d0
        t4=0.d0
        do l=-3*m,3*m
           if ( l-2*i.ge.-m .and. l-2*i.le.m  .and. &
                l-2*j.ge.-m .and. l-2*j.le.m ) then
              t1=t1+ch(l-2*i)*cht(l-2*j)
              t2=t2+cg(l-2*i)*cgt(l-2*j)
              t3=t3+ch(l-2*i)*cgt(l-2*j)
              t4=t4+cht(l-2*i)*cg(l-2*j)
           end if
        end do
        eps=1.d-10
        if (i.eq.j) then
           if (abs(t1-1.d0).gt.eps .or. abs(t2-1.d0).gt.eps .or. &
             & abs(t3).gt.eps  .or. abs(t4).gt.eps ) then 
              write(6,fmt22) 'Orthogonality ERROR', i,j,t1,t2,t3,t4
           end if
        else
           if (abs(t1).gt.eps .or. abs(t2).gt.eps .or. &
             & abs(t3).gt.eps  .or. abs(t4).gt.eps ) then
              write(6,fmt22) 'Orthogonality ERROR', i,j,t1,t2,t3,t4
           end if
        end if
     end do
  end do
  
  write(6,*) 'FILTER TEST PASSED'
  
end subroutine ftest_24
!!***


!!****h* BigDFT/scf_recursion_24
!! NAME
!!   scf_recursion_24
!!
!! FUNCTION
!!   Do iterations to go from p0gauss to pgauss
!!   16th-order interpolating scaling function
!!
!! SOURCE
!!
subroutine scf_recursion_24(n_iter,n_range,kernel_scf,kern_1_scf)
  implicit none
  !Arguments
  integer, intent(in) :: n_iter,n_range
  real(kind=8), intent(inout) :: kernel_scf(-n_range:n_range)
  real(kind=8), intent(out) :: kern_1_scf(-n_range:n_range)
  !Local variables
  real(kind=8) :: kern,kern_tot
  integer :: i_iter,i,j,ind

  include "lazy_24.inc"

  !Start the iteration to go from p0gauss to pgauss
  loop_iter_scf: do i_iter=1,n_iter
     kern_1_scf(:) = kernel_scf(:)
     kernel_scf(:) = 0.d0
     loop_iter_i: do i=0,n_range
        kern_tot = 0.d0
        do j=-m,m
           ind = 2*i-j
           if (abs(ind) > n_range) then
              kern = 0.d0
           else
              kern = kern_1_scf(ind)
           end if
           kern_tot = kern_tot + ch(j)*kern
        end do
        if (kern_tot == 0.d0) then
           !zero after (be sure because strictly == 0.d0)
           exit loop_iter_i
        else
           kernel_scf( i) = 0.5d0*kern_tot
           kernel_scf(-i) = kernel_scf(i)
        end if
     end do loop_iter_i
  end do loop_iter_scf
end subroutine scf_recursion_24
!!***


!!****h* BigDFT/for_trans_30
!! NAME
!!   for_trans_30
!!
!! FUNCTION
!!   forward wavelet transform
!!   nd: length of data set
!!   nt length of data in data set to be transformed
!!   m filter length (m has to be even!)
!!   x input data, y output data
!!
!! SOURCE
!!
subroutine for_trans_30(nd,nt,x,y)
  implicit none
  !Arguments
  integer, intent(in) :: nd,nt
  real(kind=8), intent(in) :: x(0:nd-1)
  real(kind=8), intent(out) :: y(0:nd-1)
  !Local variables
  integer :: i,j,ind

  include 'lazy_30.inc'

  do i=0,nt/2-1
     y(     i)=0.d0
     y(nt/2+i)=0.d0
     
     do j=-m+1,m
        
        ! periodically wrap index if necessary
        ind=j+2*i
        loop99: do
           if (ind.lt.0) then 
              ind=ind+nt
              cycle loop99
           end if
           if (ind.ge.nt) then 
              ind=ind-nt
              cycle loop99
           end if
           exit loop99
        end do loop99

        y(     i)=y(     i)+cht(j)*x(ind)
        y(nt/2+i)=y(nt/2+i)+cgt(j)*x(ind)
     end do
  end do
  
end subroutine for_trans_30
!!***


!!****h* BigDFT/back_trans_30
!! NAME
!!   back_trans_30
!!
!! FUNCTION
!!   backward wavelet transform
!!   nd: length of data set
!!   nt length of data in data set to be transformed
!!   m filter length (m has to be even!)
!!   x input data, y output data
!!
!! SOURCE
!!
subroutine back_trans_30(nd,nt,x,y)
  implicit none
  !Arguments
  integer, intent(in) :: nd,nt
  real(kind=8), intent(in) :: x(0:nd-1)
  real(kind=8), intent(out) :: y(0:nd-1)
  !Local variables
  integer :: i,j,ind

  include 'lazy_30.inc'
  
  do i=0,nt/2-1
     y(2*i+0)=0.d0
     y(2*i+1)=0.d0
     
     do j=-m/2,m/2-1
        
        ! periodically wrap index if necessary
        ind=i-j
        loop99: do
           if (ind.lt.0) then 
              ind=ind+nt/2
              cycle loop99
           end if
           if (ind.ge.nt/2) then 
              ind=ind-nt/2
              cycle loop99
           end if
           exit loop99
        end do loop99

        y(2*i+0)=y(2*i+0) + ch(2*j-0)*x(ind)+cg(2*j-0)*x(ind+nt/2)
        y(2*i+1)=y(2*i+1) + ch(2*j+1)*x(ind)+cg(2*j+1)*x(ind+nt/2)
     end do
  end do
        
end subroutine back_trans_30
!!***


!!****h* BigDFT/ftest_30
!! NAME
!!   ftest_30
!!
!! FUNCTION
!!   Tests the 4 orthogonality relations of the filters
!!
!! SOURCE
!!
subroutine ftest_30
  implicit none
  !Arguments
  !Local variables
  character(len=*), parameter :: fmt22 = "(a,i3,i4,4(e17.10))"
  real(kind=8) :: t1,t2,t3,t4,eps
  integer :: i,j,l

  include 'lazy_30.inc'
  
  ! do i=-m,m
  ! write(6,*) i,ch(i),cg(i)
  ! end do
  
  do i=-m,m
     do j=-m,m
        t1=0.d0
        t2=0.d0
        t3=0.d0
        t4=0.d0
        do l=-3*m,3*m
           if ( l-2*i.ge.-m .and. l-2*i.le.m  .and. &
                l-2*j.ge.-m .and. l-2*j.le.m ) then
              t1=t1+ch(l-2*i)*cht(l-2*j)
              t2=t2+cg(l-2*i)*cgt(l-2*j)
              t3=t3+ch(l-2*i)*cgt(l-2*j)
              t4=t4+cht(l-2*i)*cg(l-2*j)
           end if
        end do
        eps=1.d-10
        if (i.eq.j) then
           if (abs(t1-1.d0).gt.eps .or. abs(t2-1.d0).gt.eps .or. &
             & abs(t3).gt.eps  .or. abs(t4).gt.eps ) then 
              write(6,fmt22) 'Orthogonality ERROR', i,j,t1,t2,t3,t4
           end if
        else
           if (abs(t1).gt.eps .or. abs(t2).gt.eps .or. &
             & abs(t3).gt.eps  .or. abs(t4).gt.eps ) then
              write(6,fmt22) 'Orthogonality ERROR', i,j,t1,t2,t3,t4
           end if
        end if
     end do
  end do
  
  write(6,*) 'FILTER TEST PASSED'
  
end subroutine ftest_30
!!***


!!****h* BigDFT/scf_recursion_30
!! NAME
!!   scf_recursion_30
!!
!! FUNCTION
!!   Do iterations to go from p0gauss to pgauss
!!   16th-order interpolating scaling function
!!
!! SOURCE
!!
subroutine scf_recursion_30(n_iter,n_range,kernel_scf,kern_1_scf)
  implicit none
  !Arguments
  integer, intent(in) :: n_iter,n_range
  real(kind=8), intent(inout) :: kernel_scf(-n_range:n_range)
  real(kind=8), intent(out) :: kern_1_scf(-n_range:n_range)
  !Local variables
  real(kind=8) :: kern,kern_tot
  integer :: i_iter,i,j,ind

  include "lazy_30.inc"

  !Start the iteration to go from p0gauss to pgauss
  loop_iter_scf: do i_iter=1,n_iter
     kern_1_scf(:) = kernel_scf(:)
     kernel_scf(:) = 0.d0
     loop_iter_i: do i=0,n_range
        kern_tot = 0.d0
        do j=-m,m
           ind = 2*i-j
           if (abs(ind) > n_range) then
              kern = 0.d0
           else
              kern = kern_1_scf(ind)
           end if
           kern_tot = kern_tot + ch(j)*kern
        end do
        if (kern_tot == 0.d0) then
           !zero after (be sure because strictly == 0.d0)
           exit loop_iter_i
        else
           kernel_scf( i) = 0.5d0*kern_tot
           kernel_scf(-i) = kernel_scf(i)
        end if
     end do loop_iter_i
  end do loop_iter_scf
end subroutine scf_recursion_30
!!***


!!****h* BigDFT/for_trans_40
!! NAME
!!   for_trans_40
!!
!! FUNCTION
!!   forward wavelet transform
!!   nd: length of data set
!!   nt length of data in data set to be transformed
!!   m filter length (m has to be even!)
!!   x input data, y output data
!!
!! SOURCE
!!
subroutine for_trans_40(nd,nt,x,y)
  implicit none
  !Arguments
  integer, intent(in) :: nd,nt
  real(kind=8), intent(in) :: x(0:nd-1)
  real(kind=8), intent(out) :: y(0:nd-1)
  !Local variables
  integer :: i,j,ind

  include 'lazy_40.inc'

  do i=0,nt/2-1
     y(     i)=0.d0
     y(nt/2+i)=0.d0
     
     do j=-m+1,m
        
        ! periodically wrap index if necessary
        ind=j+2*i
        loop99: do
           if (ind.lt.0) then 
              ind=ind+nt
              cycle loop99
           end if
           if (ind.ge.nt) then 
              ind=ind-nt
              cycle loop99
           end if
           exit loop99
        end do loop99

        y(     i)=y(     i)+cht(j)*x(ind)
        y(nt/2+i)=y(nt/2+i)+cgt(j)*x(ind)
     end do
  end do
  
end subroutine for_trans_40
!!***


!!****h* BigDFT/back_trans_40
!! NAME
!!   back_trans_40
!!
!! FUNCTION
!!   backward wavelet transform
!!   nd: length of data set
!!   nt length of data in data set to be transformed
!!   m filter length (m has to be even!)
!!   x input data, y output data
!!
!! SOURCE
!!
subroutine back_trans_40(nd,nt,x,y)
  implicit none
  !Arguments
  integer, intent(in) :: nd,nt
  real(kind=8), intent(in) :: x(0:nd-1)
  real(kind=8), intent(out) :: y(0:nd-1)
  !Local variables
  integer :: i,j,ind

  include 'lazy_40.inc'
  
  do i=0,nt/2-1
     y(2*i+0)=0.d0
     y(2*i+1)=0.d0
     
     do j=-m/2,m/2-1
        
        ! periodically wrap index if necessary
        ind=i-j
        loop99: do
           if (ind.lt.0) then 
              ind=ind+nt/2
              cycle loop99
           end if
           if (ind.ge.nt/2) then 
              ind=ind-nt/2
              cycle loop99
           end if
           exit loop99
        end do loop99

        y(2*i+0)=y(2*i+0) + ch(2*j-0)*x(ind)+cg(2*j-0)*x(ind+nt/2)
        y(2*i+1)=y(2*i+1) + ch(2*j+1)*x(ind)+cg(2*j+1)*x(ind+nt/2)
     end do
  end do
        
end subroutine back_trans_40
!!***


!!****h* BigDFT/ftest_40
!! NAME
!!   ftest_40
!!
!! FUNCTION
!!   Tests the 4 orthogonality relations of the filters
!!
!! SOURCE
!!
subroutine ftest_40
  implicit none
  !Arguments
  !Local variables
  character(len=*), parameter :: fmt22 = "(a,i3,i4,4(e17.10))"
  real(kind=8) :: t1,t2,t3,t4,eps
  integer :: i,j,l

  include 'lazy_40.inc'
  
  ! do i=-m,m
  ! write(6,*) i,ch(i),cg(i)
  ! end do
  
  do i=-m,m
     do j=-m,m
        t1=0.d0
        t2=0.d0
        t3=0.d0
        t4=0.d0
        do l=-3*m,3*m
           if ( l-2*i.ge.-m .and. l-2*i.le.m  .and. &
                l-2*j.ge.-m .and. l-2*j.le.m ) then
              t1=t1+ch(l-2*i)*cht(l-2*j)
              t2=t2+cg(l-2*i)*cgt(l-2*j)
              t3=t3+ch(l-2*i)*cgt(l-2*j)
              t4=t4+cht(l-2*i)*cg(l-2*j)
           end if
        end do
        eps=1.d-10
        if (i.eq.j) then
           if (abs(t1-1.d0).gt.eps .or. abs(t2-1.d0).gt.eps .or. &
             & abs(t3).gt.eps  .or. abs(t4).gt.eps ) then 
              write(6,fmt22) 'Orthogonality ERROR', i,j,t1,t2,t3,t4
           end if
        else
           if (abs(t1).gt.eps .or. abs(t2).gt.eps .or. &
             & abs(t3).gt.eps  .or. abs(t4).gt.eps ) then
              write(6,fmt22) 'Orthogonality ERROR', i,j,t1,t2,t3,t4
           end if
        end if
     end do
  end do
  
  write(6,*) 'FILTER TEST PASSED'
  
end subroutine ftest_40
!!***


!!****h* BigDFT/scf_recursion_40
!! NAME
!!   scf_recursion_40
!!
!! FUNCTION
!!   Do iterations to go from p0gauss to pgauss
!!   16th-order interpolating scaling function
!!
!! SOURCE
!!
subroutine scf_recursion_40(n_iter,n_range,kernel_scf,kern_1_scf)
  implicit none
  !Arguments
  integer, intent(in) :: n_iter,n_range
  real(kind=8), intent(inout) :: kernel_scf(-n_range:n_range)
  real(kind=8), intent(out) :: kern_1_scf(-n_range:n_range)
  !Local variables
  real(kind=8) :: kern,kern_tot
  integer :: i_iter,i,j,ind

  include "lazy_40.inc"

  !Start the iteration to go from p0gauss to pgauss
  loop_iter_scf: do i_iter=1,n_iter
     kern_1_scf(:) = kernel_scf(:)
     kernel_scf(:) = 0.d0
     loop_iter_i: do i=0,n_range
        kern_tot = 0.d0
        do j=-m,m
           ind = 2*i-j
           if (abs(ind) > n_range) then
              kern = 0.d0
           else
              kern = kern_1_scf(ind)
           end if
           kern_tot = kern_tot + ch(j)*kern
        end do
        if (kern_tot == 0.d0) then
           !zero after (be sure because strictly == 0.d0)
           exit loop_iter_i
        else
           kernel_scf( i) = 0.5d0*kern_tot
           kernel_scf(-i) = kernel_scf(i)
        end if
     end do loop_iter_i
  end do loop_iter_scf
end subroutine scf_recursion_40
!!***


!!****h* BigDFT/for_trans_50
!! NAME
!!   for_trans_50
!!
!! FUNCTION
!!   forward wavelet transform
!!   nd: length of data set
!!   nt length of data in data set to be transformed
!!   m filter length (m has to be even!)
!!   x input data, y output data
!!
!! SOURCE
!!
subroutine for_trans_50(nd,nt,x,y)
  implicit none
  !Arguments
  integer, intent(in) :: nd,nt
  real(kind=8), intent(in) :: x(0:nd-1)
  real(kind=8), intent(out) :: y(0:nd-1)
  !Local variables
  integer :: i,j,ind

  include 'lazy_50.inc'

  do i=0,nt/2-1
     y(     i)=0.d0
     y(nt/2+i)=0.d0
     
     do j=-m+1,m
        
        ! periodically wrap index if necessary
        ind=j+2*i
        loop99: do
           if (ind.lt.0) then 
              ind=ind+nt
              cycle loop99
           end if
           if (ind.ge.nt) then 
              ind=ind-nt
              cycle loop99
           end if
           exit loop99
        end do loop99

        y(     i)=y(     i)+cht(j)*x(ind)
        y(nt/2+i)=y(nt/2+i)+cgt(j)*x(ind)
     end do
  end do
  
end subroutine for_trans_50
!!***


!!****h* BigDFT/back_trans_50
!! NAME
!!   back_trans_50
!!
!! FUNCTION
!!   backward wavelet transform
!!   nd: length of data set
!!   nt length of data in data set to be transformed
!!   m filter length (m has to be even!)
!!   x input data, y output data
!!
!! SOURCE
!!
subroutine back_trans_50(nd,nt,x,y)
  implicit none
  !Arguments
  integer, intent(in) :: nd,nt
  real(kind=8), intent(in) :: x(0:nd-1)
  real(kind=8), intent(out) :: y(0:nd-1)
  !Local variables
  integer :: i,j,ind

  include 'lazy_50.inc'
  
  do i=0,nt/2-1
     y(2*i+0)=0.d0
     y(2*i+1)=0.d0
     
     do j=-m/2,m/2-1
        
        ! periodically wrap index if necessary
        ind=i-j
        loop99: do
           if (ind.lt.0) then 
              ind=ind+nt/2
              cycle loop99
           end if
           if (ind.ge.nt/2) then 
              ind=ind-nt/2
              cycle loop99
           end if
           exit loop99
        end do loop99

        y(2*i+0)=y(2*i+0) + ch(2*j-0)*x(ind)+cg(2*j-0)*x(ind+nt/2)
        y(2*i+1)=y(2*i+1) + ch(2*j+1)*x(ind)+cg(2*j+1)*x(ind+nt/2)
     end do
  end do
        
end subroutine back_trans_50
!!***


!!****h* BigDFT/ftest_50
!! NAME
!!   ftest_50
!!
!! FUNCTION
!!   Tests the 4 orthogonality relations of the filters
!!
!! SOURCE
!!
subroutine ftest_50
  implicit none
  !Arguments
  !Local variables
  character(len=*), parameter :: fmt22 = "(a,i3,i4,4(e17.10))"
  real(kind=8) :: t1,t2,t3,t4,eps
  integer :: i,j,l

  include 'lazy_50.inc'
  
  ! do i=-m,m
  ! write(6,*) i,ch(i),cg(i)
  ! end do
  
  do i=-m,m
     do j=-m,m
        t1=0.d0
        t2=0.d0
        t3=0.d0
        t4=0.d0
        do l=-3*m,3*m
           if ( l-2*i.ge.-m .and. l-2*i.le.m  .and. &
                l-2*j.ge.-m .and. l-2*j.le.m ) then
              t1=t1+ch(l-2*i)*cht(l-2*j)
              t2=t2+cg(l-2*i)*cgt(l-2*j)
              t3=t3+ch(l-2*i)*cgt(l-2*j)
              t4=t4+cht(l-2*i)*cg(l-2*j)
           end if
        end do
        eps=1.d-10
        if (i.eq.j) then
           if (abs(t1-1.d0).gt.eps .or. abs(t2-1.d0).gt.eps .or. &
             & abs(t3).gt.eps  .or. abs(t4).gt.eps ) then 
              write(6,fmt22) 'Orthogonality ERROR', i,j,t1,t2,t3,t4
           end if
        else
           if (abs(t1).gt.eps .or. abs(t2).gt.eps .or. &
             & abs(t3).gt.eps  .or. abs(t4).gt.eps ) then
              write(6,fmt22) 'Orthogonality ERROR', i,j,t1,t2,t3,t4
           end if
        end if
     end do
  end do
  
  write(6,*) 'FILTER TEST PASSED'
  
end subroutine ftest_50
!!***


!!****h* BigDFT/scf_recursion_50
!! NAME
!!   scf_recursion_50
!!
!! FUNCTION
!!   Do iterations to go from p0gauss to pgauss
!!   16th-order interpolating scaling function
!!
!! SOURCE
!!
subroutine scf_recursion_50(n_iter,n_range,kernel_scf,kern_1_scf)
  implicit none
  !Arguments
  integer, intent(in) :: n_iter,n_range
  real(kind=8), intent(inout) :: kernel_scf(-n_range:n_range)
  real(kind=8), intent(out) :: kern_1_scf(-n_range:n_range)
  !Local variables
  real(kind=8) :: kern,kern_tot
  integer :: i_iter,i,j,ind

  include "lazy_50.inc"

  !Start the iteration to go from p0gauss to pgauss
  loop_iter_scf: do i_iter=1,n_iter
     kern_1_scf(:) = kernel_scf(:)
     kernel_scf(:) = 0.d0
     loop_iter_i: do i=0,n_range
        kern_tot = 0.d0
        do j=-m,m
           ind = 2*i-j
           if (abs(ind) > n_range) then
              kern = 0.d0
           else
              kern = kern_1_scf(ind)
           end if
           kern_tot = kern_tot + ch(j)*kern
        end do
        if (kern_tot == 0.d0) then
           !zero after (be sure because strictly == 0.d0)
           exit loop_iter_i
        else
           kernel_scf( i) = 0.5d0*kern_tot
           kernel_scf(-i) = kernel_scf(i)
        end if
     end do loop_iter_i
  end do loop_iter_scf
end subroutine scf_recursion_50
!!***


!!****h* BigDFT/for_trans_60
!! NAME
!!   for_trans_60
!!
!! FUNCTION
!!   forward wavelet transform
!!   nd: length of data set
!!   nt length of data in data set to be transformed
!!   m filter length (m has to be even!)
!!   x input data, y output data
!!
!! SOURCE
!!
subroutine for_trans_60(nd,nt,x,y)
  implicit none
  !Arguments
  integer, intent(in) :: nd,nt
  real(kind=8), intent(in) :: x(0:nd-1)
  real(kind=8), intent(out) :: y(0:nd-1)
  !Local variables
  integer :: i,j,ind

  include 'lazy_60.inc'

  do i=0,nt/2-1
     y(     i)=0.d0
     y(nt/2+i)=0.d0
     
     do j=-m+1,m
        
        ! periodically wrap index if necessary
        ind=j+2*i
        loop99: do
           if (ind.lt.0) then 
              ind=ind+nt
              cycle loop99
           end if
           if (ind.ge.nt) then 
              ind=ind-nt
              cycle loop99
           end if
           exit loop99
        end do loop99

        y(     i)=y(     i)+cht(j)*x(ind)
        y(nt/2+i)=y(nt/2+i)+cgt(j)*x(ind)
     end do
  end do
  
end subroutine for_trans_60
!!***


!!****h* BigDFT/back_trans_60
!! NAME
!!   back_trans_60
!!
!! FUNCTION
!!   backward wavelet transform
!!   nd: length of data set
!!   nt length of data in data set to be transformed
!!   m filter length (m has to be even!)
!!   x input data, y output data
!!
!! SOURCE
!!
subroutine back_trans_60(nd,nt,x,y)
  implicit none
  !Arguments
  integer, intent(in) :: nd,nt
  real(kind=8), intent(in) :: x(0:nd-1)
  real(kind=8), intent(out) :: y(0:nd-1)
  !Local variables
  integer :: i,j,ind

  include 'lazy_60.inc'
  
  do i=0,nt/2-1
     y(2*i+0)=0.d0
     y(2*i+1)=0.d0
     
     do j=-m/2,m/2-1
        
        ! periodically wrap index if necessary
        ind=i-j
        loop99: do
           if (ind.lt.0) then 
              ind=ind+nt/2
              cycle loop99
           end if
           if (ind.ge.nt/2) then 
              ind=ind-nt/2
              cycle loop99
           end if
           exit loop99
        end do loop99

        y(2*i+0)=y(2*i+0) + ch(2*j-0)*x(ind)+cg(2*j-0)*x(ind+nt/2)
        y(2*i+1)=y(2*i+1) + ch(2*j+1)*x(ind)+cg(2*j+1)*x(ind+nt/2)
     end do
  end do
        
end subroutine back_trans_60
!!***


!!****h* BigDFT/ftest_60
!! NAME
!!   ftest_60
!!
!! FUNCTION
!!   Tests the 4 orthogonality relations of the filters
!!
!! SOURCE
!!
subroutine ftest_60
  implicit none
  !Arguments
  !Local variables
  character(len=*), parameter :: fmt22 = "(a,i3,i4,4(e17.10))"
  real(kind=8) :: t1,t2,t3,t4,eps
  integer :: i,j,l

  include 'lazy_60.inc'
  
  ! do i=-m,m
  ! write(6,*) i,ch(i),cg(i)
  ! end do
  
  do i=-m,m
     do j=-m,m
        t1=0.d0
        t2=0.d0
        t3=0.d0
        t4=0.d0
        do l=-3*m,3*m
           if ( l-2*i.ge.-m .and. l-2*i.le.m  .and. &
                l-2*j.ge.-m .and. l-2*j.le.m ) then
              t1=t1+ch(l-2*i)*cht(l-2*j)
              t2=t2+cg(l-2*i)*cgt(l-2*j)
              t3=t3+ch(l-2*i)*cgt(l-2*j)
              t4=t4+cht(l-2*i)*cg(l-2*j)
           end if
        end do
        eps=1.d-10
        if (i.eq.j) then
           if (abs(t1-1.d0).gt.eps .or. abs(t2-1.d0).gt.eps .or. &
             & abs(t3).gt.eps  .or. abs(t4).gt.eps ) then 
              write(6,fmt22) 'Orthogonality ERROR', i,j,t1,t2,t3,t4
           end if
        else
           if (abs(t1).gt.eps .or. abs(t2).gt.eps .or. &
             & abs(t3).gt.eps  .or. abs(t4).gt.eps ) then
              write(6,fmt22) 'Orthogonality ERROR', i,j,t1,t2,t3,t4
           end if
        end if
     end do
  end do
  
  write(6,*) 'FILTER TEST PASSED'
  
end subroutine ftest_60
!!***


!!****h* BigDFT/scf_recursion_60
!! NAME
!!   scf_recursion_60
!!
!! FUNCTION
!!   Do iterations to go from p0gauss to pgauss
!!   16th-order interpolating scaling function
!!
!! SOURCE
!!
subroutine scf_recursion_60(n_iter,n_range,kernel_scf,kern_1_scf)
  implicit none
  !Arguments
  integer, intent(in) :: n_iter,n_range
  real(kind=8), intent(inout) :: kernel_scf(-n_range:n_range)
  real(kind=8), intent(out) :: kern_1_scf(-n_range:n_range)
  !Local variables
  real(kind=8) :: kern,kern_tot
  integer :: i_iter,i,j,ind

  include "lazy_60.inc"

  !Start the iteration to go from p0gauss to pgauss
  loop_iter_scf: do i_iter=1,n_iter
     kern_1_scf(:) = kernel_scf(:)
     kernel_scf(:) = 0.d0
     loop_iter_i: do i=0,n_range
        kern_tot = 0.d0
        do j=-m,m
           ind = 2*i-j
           if (abs(ind) > n_range) then
              kern = 0.d0
           else
              kern = kern_1_scf(ind)
           end if
           kern_tot = kern_tot + ch(j)*kern
        end do
        if (kern_tot == 0.d0) then
           !zero after (be sure because strictly == 0.d0)
           exit loop_iter_i
        else
           kernel_scf( i) = 0.5d0*kern_tot
           kernel_scf(-i) = kernel_scf(i)
        end if
     end do loop_iter_i
  end do loop_iter_scf
end subroutine scf_recursion_60
!!***


!!****h* BigDFT/for_trans_100
!! NAME
!!   for_trans_100
!!
!! FUNCTION
!!   forward wavelet transform
!!   nd: length of data set
!!   nt length of data in data set to be transformed
!!   m filter length (m has to be even!)
!!   x input data, y output data
!!
!! SOURCE
!!
subroutine for_trans_100(nd,nt,x,y)
  implicit none
  !Arguments
  integer, intent(in) :: nd,nt
  real(kind=8), intent(in) :: x(0:nd-1)
  real(kind=8), intent(out) :: y(0:nd-1)
  !Local variables
  integer :: i,j,ind

  include 'lazy_100.inc'

  do i=0,nt/2-1
     y(     i)=0.d0
     y(nt/2+i)=0.d0
     
     do j=-m+1,m
        
        ! periodically wrap index if necessary
        ind=j+2*i
        loop99: do
           if (ind.lt.0) then 
              ind=ind+nt
              cycle loop99
           end if
           if (ind.ge.nt) then 
              ind=ind-nt
              cycle loop99
           end if
           exit loop99
        end do loop99

        y(     i)=y(     i)+cht(j)*x(ind)
        y(nt/2+i)=y(nt/2+i)+cgt(j)*x(ind)
     end do
  end do
  
end subroutine for_trans_100
!!***


!!****h* BigDFT/back_trans_100
!! NAME
!!   back_trans_100
!!
!! FUNCTION
!!   backward wavelet transform
!!   nd: length of data set
!!   nt length of data in data set to be transformed
!!   m filter length (m has to be even!)
!!   x input data, y output data
!!
!! SOURCE
!!
subroutine back_trans_100(nd,nt,x,y)
  implicit none
  !Arguments
  integer, intent(in) :: nd,nt
  real(kind=8), intent(in) :: x(0:nd-1)
  real(kind=8), intent(out) :: y(0:nd-1)
  !Local variables
  integer :: i,j,ind

  include 'lazy_100.inc'
  
  do i=0,nt/2-1
     y(2*i+0)=0.d0
     y(2*i+1)=0.d0
     
     do j=-m/2,m/2-1
        
        ! periodically wrap index if necessary
        ind=i-j
        loop99: do
           if (ind.lt.0) then 
              ind=ind+nt/2
              cycle loop99
           end if
           if (ind.ge.nt/2) then 
              ind=ind-nt/2
              cycle loop99
           end if
           exit loop99
        end do loop99

        y(2*i+0)=y(2*i+0) + ch(2*j-0)*x(ind)+cg(2*j-0)*x(ind+nt/2)
        y(2*i+1)=y(2*i+1) + ch(2*j+1)*x(ind)+cg(2*j+1)*x(ind+nt/2)
     end do
  end do
        
end subroutine back_trans_100
!!***


!!****h* BigDFT/ftest_100
!! NAME
!!   ftest_100
!!
!! FUNCTION
!!   Tests the 4 orthogonality relations of the filters
!!
!! SOURCE
!!
subroutine ftest_100
  implicit none
  !Arguments
  !Local variables
  character(len=*), parameter :: fmt22 = "(a,i3,i4,4(e17.10))"
  real(kind=8) :: t1,t2,t3,t4,eps
  integer :: i,j,l

  include 'lazy_100.inc'
  
  ! do i=-m,m
  ! write(6,*) i,ch(i),cg(i)
  ! end do
  
  do i=-m,m
     do j=-m,m
        t1=0.d0
        t2=0.d0
        t3=0.d0
        t4=0.d0
        do l=-3*m,3*m
           if ( l-2*i.ge.-m .and. l-2*i.le.m  .and. &
                l-2*j.ge.-m .and. l-2*j.le.m ) then
              t1=t1+ch(l-2*i)*cht(l-2*j)
              t2=t2+cg(l-2*i)*cgt(l-2*j)
              t3=t3+ch(l-2*i)*cgt(l-2*j)
              t4=t4+cht(l-2*i)*cg(l-2*j)
           end if
        end do
        eps=1.d-10
        if (i.eq.j) then
           if (abs(t1-1.d0).gt.eps .or. abs(t2-1.d0).gt.eps .or. &
             & abs(t3).gt.eps  .or. abs(t4).gt.eps ) then 
              write(6,fmt22) 'Orthogonality ERROR', i,j,t1,t2,t3,t4
           end if
        else
           if (abs(t1).gt.eps .or. abs(t2).gt.eps .or. &
             & abs(t3).gt.eps  .or. abs(t4).gt.eps ) then
              write(6,fmt22) 'Orthogonality ERROR', i,j,t1,t2,t3,t4
           end if
        end if
     end do
  end do
  
  write(6,*) 'FILTER TEST PASSED'
  
end subroutine ftest_100
!!***


!!****h* BigDFT/scf_recursion_100
!! NAME
!!   scf_recursion_100
!!
!! FUNCTION
!!   Do iterations to go from p0gauss to pgauss
!!   16th-order interpolating scaling function
!!
!! SOURCE
!!
subroutine scf_recursion_100(n_iter,n_range,kernel_scf,kern_1_scf)
  implicit none
  !Arguments
  integer, intent(in) :: n_iter,n_range
  real(kind=8), intent(inout) :: kernel_scf(-n_range:n_range)
  real(kind=8), intent(out) :: kern_1_scf(-n_range:n_range)
  !Local variables
  real(kind=8) :: kern,kern_tot
  integer :: i_iter,i,j,ind

  include "lazy_100.inc"

  !Start the iteration to go from p0gauss to pgauss
  loop_iter_scf: do i_iter=1,n_iter
     kern_1_scf(:) = kernel_scf(:)
     kernel_scf(:) = 0.d0
     loop_iter_i: do i=0,n_range
        kern_tot = 0.d0
        do j=-m,m
           ind = 2*i-j
           if (abs(ind) > n_range) then
              kern = 0.d0
           else
              kern = kern_1_scf(ind)
           end if
           kern_tot = kern_tot + ch(j)*kern
        end do
        if (kern_tot == 0.d0) then
           !zero after (be sure because strictly == 0.d0)
           exit loop_iter_i
        else
           kernel_scf( i) = 0.5d0*kern_tot
           kernel_scf(-i) = kernel_scf(i)
        end if
     end do loop_iter_i
  end do loop_iter_scf
end subroutine scf_recursion_100
!!***
