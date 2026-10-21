program test_lfff_fft_modes
  use mod_lfff, only: lfff_fft_transfer,lfff_balance_bottom_flux
  implicit none

  double precision :: b,d,q,pi
  double precision :: bz(2,2),before,after,correction
  integer :: status

  pi=4.d0*datan(1.d0)

  call lfff_fft_transfer(4.d0,0.d0,0.3d0,1.d0,.false.,b,d,status)
  q=2.d0
  call require(status==0,'open potential mode returned an error')
  call require_close(b,dexp(-q*0.3d0),'open potential Bz transfer')
  call require_close(d,q*dexp(-q*0.3d0),'open potential derivative transfer')

  call lfff_fft_transfer(4.d0,0.d0,0.d0,1.d0,.true.,b,d,status)
  call require(status==0,'closed potential bottom returned an error')
  call require_close(b,1.d0,'closed potential bottom Bz')
  call lfff_fft_transfer(4.d0,0.d0,1.d0,1.d0,.true.,b,d,status)
  call require_close(b,0.d0,'closed potential top Bz')

  q=dsqrt(3.d0)
  call lfff_fft_transfer(4.d0,1.d0,0.3d0,1.d0,.false.,b,d,status)
  call require(status==0,'open LFFF mode returned an error')
  call require_close(b,dexp(-q*0.3d0),'open LFFF Bz transfer')
  call require_close(d,q*dexp(-q*0.3d0),'open LFFF derivative transfer')

  call lfff_fft_transfer(0.25d0,1.d0,0.3d0,1.d0,.false.,b,d,status)
  call require(status==1,'open oscillatory LFFF mode was not rejected')

  call lfff_fft_transfer(4.d0,1.d0,0.d0,1.d0,.true.,b,d,status)
  call require_close(b,1.d0,'closed evanescent LFFF bottom Bz')
  call lfff_fft_transfer(4.d0,1.d0,1.d0,1.d0,.true.,b,d,status)
  call require_close(b,0.d0,'closed evanescent LFFF top Bz')

  call lfff_fft_transfer(0.25d0,1.d0,0.d0,1.d0,.true.,b,d,status)
  call require(status==0,'closed oscillatory LFFF mode returned an error')
  call require_close(b,1.d0,'closed oscillatory LFFF bottom Bz')
  call lfff_fft_transfer(0.25d0,1.d0,1.d0,1.d0,.true.,b,d,status)
  call require_close(b,0.d0,'closed oscillatory LFFF top Bz')

  call lfff_fft_transfer(1.d0,1.d0,0.25d0,1.d0,.true.,b,d,status)
  call require(status==0,'critical closed LFFF mode returned an error')
  call require_close(b,0.75d0,'critical closed LFFF Bz transfer')
  call require_close(d,1.d0,'critical closed LFFF derivative transfer')

  call lfff_fft_transfer(1.d0,dsqrt(pi*pi+1.d0),0.5d0,1.d0,.true.,b,d,status)
  call require(status==2,'closed LFFF resonance was not rejected')

  bz=reshape((/2.d0,-2.d0,1.d0,-0.8d0/),shape(bz))
  call lfff_balance_bottom_flux(bz,'strict',0.1d0,before,after,correction,status)
  call require(status==4,'strict flux treatment accepted an imbalanced map')
  call require_close(before,0.2d0/5.8d0,'strict input flux imbalance')
  call require_close(sum(bz),0.2d0,'strict treatment changed the magnetogram')

  call lfff_balance_bottom_flux(bz,'subtract_mean',0.1d0,before,after,&
     correction,status)
  call require(status==1,'mean-subtraction treatment was not applied')
  call require_close(correction,0.05d0,'mean-subtraction correction')
  call require_close(sum(bz),0.d0,'mean-subtraction residual flux')
  call require(after<=1.d-14,'mean-subtraction residual imbalance')

  bz=reshape((/1.d0,1.d0,1.d0,-0.1d0/),shape(bz))
  call lfff_balance_bottom_flux(bz,'subtract_mean',0.1d0,before,after,&
     correction,status)
  call require(status==5,'nearly unipolar map bypassed the safety threshold')
  call require_close(sum(bz),2.9d0,'rejected treatment changed the magnetogram')

  print '(a)', 'FFT potential/LFFF mode tests passed'

contains

  subroutine require(condition,message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if(.not.condition) then
      write(*,*) trim(message)
      error stop
    end if
  end subroutine require

  subroutine require_close(actual,expected,message)
    double precision, intent(in) :: actual,expected
    character(len=*), intent(in) :: message
    call require(dabs(actual-expected)<=1.d-12*max(1.d0,dabs(expected)),message)
  end subroutine require_close

end program test_lfff_fft_modes
