program test_fft_core
  use mod_fft
  implicit none

  integer, parameter :: n1=12,n2=15
  integer, parameter :: sizes(5)=[12,15,20,240,324]
  double complex :: input(n1,n2),transformed(n1,n2),reference(n1,n2),phase
  double complex, allocatable :: line(:),original(:)
  double precision :: real_part(n1,n2),imag_part(n1,n2),error,roundtrip_error,twopi
  integer :: i,j,p,q,n

  twopi=8.d0*datan(1.d0)
  call random_number(real_part)
  call random_number(imag_part)
  input=dcmplx(real_part,imag_part)
  transformed=input
  call fft_2d(transformed,.false.)

  reference=(0.d0,0.d0)
  do q=1,n2
    do p=1,n1
      do j=1,n2
        do i=1,n1
          phase=exp(dcmplx(0.d0,-twopi*(dble((p-1)*(i-1))/n1+&
             dble((q-1)*(j-1))/n2)))
          reference(p,q)=reference(p,q)+input(i,j)*phase
        end do
      end do
    end do
  end do
  error=maxval(abs(transformed-reference))
  if(error>1.d-11) error stop 'FFT does not match direct DFT'

  call fft_2d(transformed,.true.)
  error=maxval(abs(transformed-input))
  if(error>1.d-12) error stop '2D FFT round trip failed'
  roundtrip_error=error

  do i=1,size(sizes)
    n=sizes(i)
    allocate(line(n),original(n))
    line=[(dcmplx(dsin(0.13d0*j),dcos(0.17d0*j)),j=1,n)]
    original=line
    call fft_1d(line,.false.)
    call fft_1d(line,.true.)
    error=maxval(abs(line-original))
    if(error>1.d-12) error stop '1D FFT round trip failed'
    deallocate(line,original)
  end do

  if(fft_size_supported(29)) error stop '29 should be unsupported'
  if(fft_next_supported(29)/=30) error stop 'next supported size should be 30'
  print '(a,es12.4)', 'FFT core tests passed; 2D round-trip error=',roundtrip_error
end program test_fft_core
