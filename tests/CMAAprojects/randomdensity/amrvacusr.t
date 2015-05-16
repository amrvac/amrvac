!=============================================================================
! amrvacusr.t.ding
! setup.pl -d=22 -phi=0 -z=0 -g=14,14 -p=hd -eos=default -nf=0 -ndust=0 -u=nul -arch=default
!=============================================================================
INCLUDE:amrvacnul/speciallog.t
INCLUDE:amrvacnul/specialsource.t
INCLUDE:amrvacnul/specialimpl.t
INCLUDE:amrvacnul/specialbound.t
INCLUDE:amrvacnul/usrflags.t
!=============================================================================
!define user specific variables 
module usr_initial
implicit none
double precision,save:: delrho,A0,beta
double precision,allocatable,save:: xgrid(:),ygrid(:)
integer,save:: Nx,Ny,ntmp
double precision,allocatable,save:: rand(:,:)
double precision,save:: rand_avg,rand_std
double precision,save:: Lx,Ly,x0,y0,Sigx,Sigy

end module usr_initial

!==================================
subroutine initglobaldata_usr
use usr_initial
include 'amrvacdef.f'

integer,dimension(:),allocatable:: seed
integer::  seed_size,ix

eqpar(gamma_)=5.d0/3.d0
Ly=xprobmax2-xprobmin2
Lx=xprobmax1-xprobmin1
nx=nint(Lx/dx(1,1))
if (.not. allocated(xgrid)) allocate(xgrid(nx))
do ix=1,Nx
 xgrid(ix)=xprobmin1+ix*dx(1,1)
enddo
ny=nint(Ly/dx(2,1))
if (.not. allocated(ygrid)) allocate(ygrid(ny))
do ix=1,Ny
 ygrid(ix)=xprobmin2+ix*dx(2,1)
enddo
ntmp=Nx*Ny
if (.not. allocated(rand)) allocate(rand(Nx,Ny))



if(mype==0)then
    print *,'Number of modes (x,y)=',nx,ny
    call random_seed(SIZE=seed_size)
    allocate(seed(seed_size))
    call random_seed(GET=seed(1:seed_size))
   ! print *,'Random number:',ntmp
   ! print *,'Seeds: ',seed
    call random_number(rand) ! generate ntmp random numbers [0,1]
   
endif
call MPI_BARRIER(icomm,ierrmpi)
if(npe>1) call MPI_BCAST(rand,ntmp,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
 ! print *,'Communication done!'
end subroutine initglobaldata_usr
!=============================================================================
subroutine initonegrid_usr(ixG^L,ix^L,w,x)
use usr_initial
! initialize one grid

include 'amrvacdef.f'

integer, intent(in) :: ixG^L,ix^L
double precision, intent(in) :: x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)
integer:: i,j,ix^D
double precision :: rho(ixG^T),pulse(ixG^T)
double precision :: rand_vec(1:Ny)
double precision :: tmp
logical :: patchw(ixG^T)
logical, save :: first=.true.
!----------------------------------------------------------------------------

if(first)then
  if(mype==0) then
     write(*,*)'Simulate 2D acoustic wave in random inhomogenueous medium!'
  endif
  first=.false.
endif

delrho=0.2
rand_avg=SUM(rand)/DBLE(ntmp)
rand_std=SQRT(SUM((rand-rand_avg)**2)/DBLE(Ntmp-1))
rho(ixG^T)=one
do ix1=ixGmin1,ixGmax1  
do ix2=ixGmin2,ixGmax2
   do i=1,Ny ! interpolate at x1 along xgrid for each y-position
     call interp1d(Nx,xgrid,rand(1:Nx,i),x(ix1,ix2,1),tmp)
     rand_vec(i)=tmp 
   enddo
    call interp1d(Ny,ygrid,rand_vec(1:Ny),x(ix1,ix2,2),tmp)
   ! scale to delrho above average value
  rho(ix1,ix2)=max(rho(ix1,ix2)+(tmp-rand_avg)/rand_std*delrho,0.1)
enddo
enddo

w(ix^S,rho_)=rho(ix^S)
w(ix^S,v1_)=zero
w(ix^S,v2_)=zero
! Point-source pertubation 
A0=0.05d0
x0=0.d0
y0=0.d0 ! start position of the pulse
Sigx=5d0
Sigy=5d0 ! x-and y-width of the pulse
w(ix^S,p_)=one*(1+A0*DEXP(-(x(ix^S,1)-x0)**2/Sigx**2)*DEXP(-(x(ix^S,2)-y0)**2/Sigy**2))
! Uncomment the following two lines for plane wave study

!x0=xprobmin1+0.1*Lx // position to launch a wave
!w(ix^S,p_)=one*(1+A0*DEXP(-(x(ix^S,1)-x0)**2/sigx**2))


patchw(ixG^S)=.false.
call conserve(ixG^L,ix^L,w,x,patchw)

end subroutine initonegrid_usr
!=============================================================================
subroutine interp1d(nin,xin,yin,xout1,yout1)
! subroutine to interpolate a values from a 1D function 
! validated on 12 June 2014 DY
! robustness improvement is suggested to account all possible scenario
! extrapolation included on 12 june 2014 DY 
  integer,intent(in) :: nin 
  DOUBLE PRECISION,intent(in) :: xin(1:nin),yin(1:nin)
  DOUBLE PRECISION,intent(in) :: xout1
  DOUBLE PRECISION,intent(out) :: yout1
  integer*4 :: ixa 
  DOUBLE PRECISION :: xa1,xb1,ya1,yb1,xtmp
 ! if (mype==0) then  
 !   write(*,*) nin
 !   write(*,*) size(xin),size(yin)
 !   write(*,*) xout1 
 ! endif
  ixa=minloc(abs(xin-xout1),1)
  xtmp=xin(ixa)
  ! not in range of xin[0] and xin[nin]
  ! use extraplotion 
  if (ixa .eq. 1) then 
     xa1=xin(1)
     xb1=xin(2)
     ya1=yin(1)
     yb1=yin(2)
    yout1=ya1+(yb1-ya1)*(xout1-xa1)/(xb1-xa1)
   ! print*,xa1,xb1
  endif
  if (ixa .eq. nin) then
     xa1=xin(ixa-1)
     xb1=xin(ixa)
     ya1=yin(ixa-1)
     yb1=yin(ixa)
   yout1=ya1+(yb1-ya1)*(xout1-xa1)/(xb1-xa1)
  endif

  if (ixa .lt. nin .and. ixa .gt. 1) then
    if (xtmp .gt. xout1) then
      ixa=ixa-1
    endif
    xa1=xin(ixa)
    xb1=xin(ixa+1)
    ya1=yin(ixa)
    yb1=yin(ixa+1)
    yout1=ya1+(yb1-ya1)*(xout1-xa1)/(xb1-xa1)
  endif
Return
End


!=============================================================================
! amrvacusr.t.ding
!=============================================================================
