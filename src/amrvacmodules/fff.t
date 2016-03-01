!##############################################################################
! module fff
! PURPOSE:
! Program to extrapolate linear force-free fields in 3D Cartesian coordinates,
! based on exact Green function method (Chiu & Hilton 1977 ApJ 212,873).
! Usage:
!1 In amrvacusr.t, put a line:
!  INCLUDE:amrvacmodules/fff.t
!2 In the subroutine initglobaldata_usr of amrvacusr.t:
!  To extrapolate a linear force free field from an analytical magnetogram, 
!  add lines like:
!  
!  logical, save :: firstusrglobaldata=.true.
! 
!  if(firstusrglobaldata) then
!
!    call init_bc_fff(600,360)
!    firstusrglobaldata=.false.
!  endif
!
!  600 x 360 is the resolution of the magnetogram. Users can use subroutine
!  init_bc_fff as an example and make a similar subroutine creating any
!  analytical magnetograms.
!  To extrapolate a linear force free field from a observed magnetogram 
!  prepared in a data file, e.g., 'hmiM720sxxxx.dat' replace 
!  "call init_bc_fff(xx,xx)" used above by
!  call init_bc_fff_data('hmiM720sxxxx.dat',Lunit,Bunit)
!  Lunit and Bunit are dimensionless unit for length and magnetic field.
!  'hmiM720sxxxx.dat' must be a binary file containing nx1,nx2,xc1,xc2,dxm1,
!  dxm2, Bz0(nx1,nx2). Integers nx1 and nx2 give the resolution of the 
!  uniform-grid magentogram. Others are double-precision floats. xc1 and xc2
!  are coordinates of the central point of the magnetogram. dxm1 and dxm2 
!  are the cell sizes for each direction, Bz0 is the vertical conponent 
!  of magetic field on the solar surface from observations.
!3 In the subroutine initonegrid_usr of amrvacusr.t,
!  add lines like:
!
!  double precision :: Bf(ixG^S,1:ndir), alpha, zshift
!
!  alpha=0.d0     ! potential field
!  !alpha=0.08d0  ! non-potential linear force-free field
!  zshift=0.05d0  ! lift your box zshift heigher to the bottom magnetogram
!  call calc_lin_fff(ixG^L,ix^L,Bf,x,alpha,zshift) 
!
! Notice that the resolution of input magnetogram must be better than the best
!  resolution of your AMR grid to have a good behavior in the bottom layer
! when zshift=0.
!============================================================================= 
module fff_global
implicit none

integer, save :: nx1,nx2
double precision, save :: Bzmax,darea
double precision, allocatable, save :: Bz0(:,:)
double precision, allocatable, save :: xa1(:),xa2(:)

end module fff_global
!============================================================================= 
subroutine init_b_fff_data(magnetogramname,qLunit,qBunit)
use fff_global

include 'amrvacdef.f'
double precision, intent(in) :: qLunit,qBunit
double precision :: xc1,xc2,dxm1,dxm2
integer, dimension(MPI_STATUS_SIZE) :: statuss
integer :: file_handle,i
character(len=*), intent(in) :: magnetogramname
logical :: aexist
!-----------------------------------------------------------------------------
! nx1,nx2 are numbers of cells for each direction
! xc1,xc2 are coordinates of the central point of the magnetogram
! dxm1,dxm2 are cell sizes for each direction
! Bz0 is the 2D Bz magnetogram
inquire(file=magnetogramname,exist=aexist)
if(.not. aexist) then
  if(mype==0) write(*,'(2a)') "can not find file:",magnetogramname
  call mpistop("no input magnetogram----init_b_fff_data")
end if
call MPI_FILE_OPEN(icomm,magnetogramname,MPI_MODE_RDONLY,MPI_INFO_NULL,&
                   file_handle,ierrmpi)
call MPI_FILE_READ_ALL(file_handle,nx1,1,MPI_INTEGER,statuss,ierrmpi)
call MPI_FILE_READ_ALL(file_handle,nx2,1,MPI_INTEGER,statuss,ierrmpi)
allocate(Bz0(nx1,nx2))
call MPI_FILE_READ_ALL(file_handle,xc1,1,MPI_DOUBLE_PRECISION,statuss,ierrmpi)
call MPI_FILE_READ_ALL(file_handle,xc2,1,MPI_DOUBLE_PRECISION,statuss,ierrmpi)
call MPI_FILE_READ_ALL(file_handle,dxm1,1,MPI_DOUBLE_PRECISION,statuss,ierrmpi)
call MPI_FILE_READ_ALL(file_handle,dxm2,1,MPI_DOUBLE_PRECISION,statuss,ierrmpi)
call MPI_FILE_READ_ALL(file_handle,Bz0,nx1*nx2,MPI_DOUBLE_PRECISION,&
                       statuss,ierrmpi)
call MPI_FILE_CLOSE(file_handle,ierrmpi)
allocate(xa1(nx1))
allocate(xa2(nx2))
xa1(nx1/2)=xc1
xa2(nx2/2)=xc2
do i=nx1/2+1,nx1
  xa1(i)=xa1(nx1/2)+dble(i-nx1/2)*dxm1
enddo
do i=nx1/2-1,1,-1
  xa1(i)=xa1(nx1/2)+dble(i-nx1/2)*dxm1
enddo
do i=nx2/2+1,nx2
  xa2(i)=xa2(nx2/2)+dble(i-nx2/2)*dxm2
enddo
do i=nx2/2-1,1,-1
  xa2(i)=xa2(nx2/2)+dble(i-nx2/2)*dxm2
enddo
! declare and define global variables Lunit and Bunit to be your length unit in
! cm and magnetic strength unit in Gauss first
dxm1=dxm1/qLunit
dxm2=dxm2/qLunit
xa1=xa1/qLunit
xa2=xa2/qLunit
darea=dxm1*dxm2
Bz0=Bz0/qBunit
Bzmax=maxval(dabs(Bz0(:,:)))

! normalize b
Bz0=Bz0/Bzmax
if(mype==0) then
  print*,'magnetogram xrange:',minval(xa1),maxval(xa1)
  print*,'magnetogram yrange:',minval(xa2),maxval(xa2)
end if

if(mype==0) then
  print*,'extrapolating 3D force-free field from an observed Bz '
  print*,'magnetogram of',nx1,'by',nx2,'pixels. Bzmax=',Bzmax
endif

end subroutine init_b_fff_data
!============================================================================= 
subroutine init_b_fff(qnx1,qnx2)
use fff_global

include 'amrvacdef.f'

double precision :: dxm1,dxm2,delx1,delx2,xo1,xo2,yo1,yo2,coB,B0
integer :: i1,i2,qnx1,qnx2
!-----------------------------------------------------------------------------
nx1=qnx1
nx2=qnx2
allocate(Bz0(nx1,nx2))
allocate(xa1(nx1))
allocate(xa2(nx2))
dxm1=(xprobmax1-xprobmin1)/dble(nx1)
dxm2=(xprobmax2-xprobmin2)/dble(nx2)
darea=dxm1*dxm2
do i1=1,nx1
  xa1(i1)=xprobmin1+(dble(i1)-0.5d0)*dxm1
enddo
do i2=1,nx2
  xa2(i2)=xprobmin2+(dble(i2)-0.5d0)*dxm2
enddo

xo1=(xprobmax1+xprobmin1)/2.d0
xo2=xo1
yo1=(xprobmax2+xprobmin2)/2.d0-(xprobmax2-xprobmin2)*0.15d0
yo2=(xprobmax2+xprobmin2)/2.d0+(xprobmax2-xprobmin2)*0.15d0
delx1=(xprobmax1-xprobmin1)*0.3d0
delx2=(xprobmax2-xprobmin2)*0.12d0

! Bz are composed by two elliptic 2D Gaussian functions to mimic a bipole 
B0=10.d0
do i2=1,nx2
  do i1=1,nx1
    Bz0(i1,i2)=B0*(dexp(-((xa1(i1)-xo1)**2/delx1**2+&
                    (xa2(i2)-yo1)**2/delx2**2)/2.d0)-&
                    dexp(-((xa1(i1)-xo2)**2/delx1**2+&
                    (xa2(i2)-yo2)**2/delx2**2)/2.d0))
  enddo
enddo

Bzmax=maxval(dabs(Bz0(:,:)))
! normalize b
Bz0=Bz0/Bzmax
if(mype==0) then
  print*,'extrapolating 3D force-free field from an analytical Bz '
  print*,'magnetogram of',nx1,'by',nx2,'pixels. Bzmax=',Bzmax
endif

end subroutine init_b_fff
!============================================================================= 
subroutine calc_lin_fff(ixI^L,ixO^L,Bf,x,alpha,zshift)
! PURPOSE: 
! Calculation to determine linear FFF from the field on 
! the lower boundary (Chiu and Hilton 1977 ApJ 212,873). 
! NOTE: Only works for Cartesian coordinates 
! INPUT: Bf,x
! OUTPUT: updated b in w 
use fff_global

include 'amrvacdef.f'

integer, intent(in) :: ixI^L, ixO^L
double precision, intent(in) :: x(ixI^S,1:ndim),alpha,zshift
double precision, intent(inout) :: Bf(ixI^S,1:ndir)

double precision :: cos_az(ixO^S),sin_az(ixO^S),zk(ixO^S)
double precision :: r, r2, r3, bigr, cos_ar, sin_ar, xsum, ysum, zsum
double precision :: g, dgdz, gx, gy, gz, twopiinv
integer :: ix1,ix2,ix3,ixp1,ixp2
!-----------------------------------------------------------------------------

Bf=0.d0
twopiinv = 0.5d0/dpi
! get cos and sin arrays
zk(ixO^S)=x(ixO^S,3)-xprobmin3+zshift
cos_az(ixO^S)=dcos(alpha*zk(ixO^S))
sin_az(ixO^S)=dsin(alpha*zk(ixO^S))
! 5 loops, for each grid integrate over x and y of Bz0
do ix3=ixOmin3,ixOmax3
  do ix2=ixOmin2,ixOmax2
    do ix1=ixOmin1,ixOmax1
      xsum = 0.d0
      ysum = 0.d0
      zsum = 0.d0
      do ixp2=1,nx2
        do ixp1=1,nx1
          bigr=dsqrt((x(ix1,ix2,ix3,1)-xa1(ixp1))**2+&
                    (x(ix1,ix2,ix3,2)-xa2(ixp2))**2)
          if(bigr>smalldouble) then
            r2=bigr**2+zk(ix1,ix2,ix3)**2
            r=dsqrt(r2)
            r3=r**3
            cos_ar=dcos(alpha*r)
            sin_ar=dsin(alpha*r)
            bigr=1.d0/bigr
            g=(zk(ix1,ix2,ix3)*cos_ar/r-cos_az(ix1,ix2,ix3))*bigr
            dgdz=(cos_ar*(1.d0/r-zk(ix1,ix2,ix3)**2/r3)&
                 -alpha*zk(ix1,ix2,ix3)**2*sin_ar/r2&
                 +alpha*sin_az(ix1,ix2,ix3))*bigr
            gx=Bz0(ixp1,ixp2)*((x(ix1,ix2,ix3,1)-xa1(ixp1))*dgdz&
               +alpha*g*(x(ix1,ix2,ix3,2)-xa2(ixp2)))*bigr
            gy=Bz0(ixp1,ixp2)*((x(ix1,ix2,ix3,2)-xa2(ixp2))*dgdz&
               -alpha*g*(x(ix1,ix2,ix3,1)-xa1(ixp1)))*bigr
            gz=Bz0(ixp1,ixp2)*(zk(ix1,ix2,ix3)*cos_ar/r3+alpha*&
               zk(ix1,ix2,ix3)*sin_ar/r2)
            xsum=xsum+gx*darea
            ysum=ysum+gy*darea
            zsum=zsum+gz*darea
          end if
        end do
      end do
      Bf(ix1,ix2,ix3,1)=xsum*twopiinv
      Bf(ix1,ix2,ix3,2)=ysum*twopiinv
      Bf(ix1,ix2,ix3,3)=zsum*twopiinv
    end do
  end do
end do
Bf=Bf*Bzmax

end subroutine calc_lin_fff
!============================================================================= 
subroutine optimization_fff
! PURPOSE: 
! To improve linear FFF to be more force-free and divergence-free 
! INPUT: w,x
! OUTPUT: updated b in w 
use fff_global

include 'amrvacdef.f'

integer :: i,iigrid, igrid, iitmax
double precision :: dtfff,dtminfff,abs_frac_diff,l_value,lf_value,delta_l
double precision :: l_grid,ld_grid,l_pe,ld_pe,l_new,ld_new
logical :: patchw(ixG^T)
!-----------------------------------------------------------------------------

if(mype==0) print*, 'improve force-free field using &
                     optimization method'

! normalize B field
do iigrid=1,igridstail; igrid=igrids(iigrid);
{#IFDEF ENERGY
  call primitive(ixG^LL,ixM^LL,pw(igrid)%w,px(igrid)%x)
}
  pw(igrid)%w(:^D&,b1_:b3_)=pw(igrid)%w(:^D&,b1_:b3_)/Bzmax
end do
call getbc(t,ixG^LL,pw,pwCoarse,pgeo,pgeoCoarse,.false.,b0_,ndir)
dtfff=1.d-5
dtminfff=1.d-16
abs_frac_diff=1.d-5
delta_l=1.d0
l_value=bigdouble
! optimization cycling
i=0
iitmax=1000
do
  if(i>=iitmax .or. dtfff < dtminfff) then
    if(mype==0) print*, 'Fail to converge! L=',l_new,' L_divb=',&
      ld_new,' dL=',delta_l,' i=',i, 'dtfff=',dtfff
    exit
  endif
  l_pe=0.d0
  ld_pe=0.d0
  do iigrid=1,igridstail; igrid=igrids(iigrid);
     ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
     call fff_grid(dtfff,ixG^LL,ixM^LL,pw(igrid)%w,px(igrid)%x,&
                   l_grid,ld_grid,.false.)
     l_pe=l_pe+l_grid
     ld_pe=ld_pe+ld_grid
  end do
  call MPI_ALLREDUCE(l_pe,l_new,1,MPI_DOUBLE_PRECISION,MPI_SUM, &
                         icomm,ierrmpi)
  call MPI_ALLREDUCE(ld_pe,ld_new,1,MPI_DOUBLE_PRECISION,MPI_SUM, &
                         icomm,ierrmpi)
  if(mype==0 .and. i==0) print*,'i=',i,'L=',l_new,' L_force=',&
                     l_new-ld_new,'L_divb=',ld_new
  delta_l=(l_value-l_new)/l_value
  if(delta_l>0.d0) then ! Success to converge
    !if(mype==0) print*,'converging i=',i
    if(dabs(delta_l) < abs_frac_diff) then
      lf_value=l_value-ld_new
      if(mype==0) print*, 'Succeed to converge! L=',l_value,&
        ' L_force=',lf_value,' L_divb=',ld_new,' i=',i, 'dtfff=',dtfff
      exit
    endif
    l_value=l_new
    dtfff=1.01d0*dtfff
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
       pwold(igrid)%w=pw(igrid)%w
       call fff_grid(dtfff,ixG^LL,ixM^LL,pw(igrid)%w,px(igrid)%x,&
                     l_grid,ld_grid,.true.)
    end do
    call getbc(t,ixG^LL,pw,pwCoarse,pgeo,pgeoCoarse,.false.,b0_,ndir)
  else ! Reset b
    dtfff=0.5d0*dtfff ! Reset dtfff
    !if(mype==0) print*,'redo i=',i,'dtfff=',dtfff
    ! roll back to previous w to redo with smaller time step
    do iigrid=1,igridstail; igrid=igrids(iigrid);
      pw(igrid)%w=pwold(igrid)%w
    end do
    call getbc(t,ixG^LL,pw,pwCoarse,pgeo,pgeoCoarse,.false.,b0_,ndir)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
       call fff_grid(dtfff,ixG^LL,ixM^LL,pw(igrid)%w,px(igrid)%x,&
                     l_grid,ld_grid,.true.)
    end do
    call getbc(t,ixG^LL,pw,pwCoarse,pgeo,pgeoCoarse,.false.,b0_,ndir)
  endif
  i=i+1
enddo
! de-normalize b and update total energy
do iigrid=1,igridstail; igrid=igrids(iigrid);
  pw(igrid)%w(:^D&,b1_:b3_)=pw(igrid)%w(:^D&,b1_:b3_)*Bzmax
{#IFDEF ENERGY
  patchw(ixG^T)=.false.
  call conserve(ixG^LL,ixM^LL,pw(igrid)%w,px(igrid)%x,patchw)
}
end do

end subroutine optimization_fff
!============================================================================= 
subroutine fff_grid(dtfff,ixI^L,ixO^L,w,x,l_grid,ld_grid,evolve)
use fff_global

include 'amrvacdef.f'

integer, intent(in) :: ixI^L, ixO^L
double precision, intent(in) :: x(ixI^S,1:ndim),dtfff
double precision, intent(inout) :: w(ixI^S,1:nw)
double precision, intent(out) :: l_grid,ld_grid
logical, intent(in) :: evolve

integer :: idirmin,idir,ix^L
double precision :: b2(ixG^T),b(ixG^T,1:ndir),j(ixG^T,1:ndir),divb(ixG^T)
double precision :: omega(ixG^T,1:ndir),omega2(ixG^T),fun(ixG^T,1:ndir)
double precision :: tmp(ixG^T),tmpvector(ixG^T,1:ndir)
double precision, external :: obj_funct_int
!-----------------------------------------------------------------------------

ix^L=ixO^L^LADD4;
if (ixI^L^LTix^L|.or.|.or.) &
   call mpistop("Need four extra layers for optimization_fff_grid")
b(ixI^S,1:ndir)=w(ixI^S,b0_+1:b0_+ndir)
! b2 = |b|^2
b2(ixI^S)=(^C&w(ixI^S,b^C_)**2+)
ix^L=ixI^L^LSUB2;
! j = curl b
call curlvector(b,ixI^L,ix^L,j,idirmin,1,ndir)
! divb = div b
call divvectorS(b,ixI^L,ix^L,divb)
! omega = j x b
omega(ix^S,1)=j(ix^S,2)*b(ix^S,3)-j(ix^S,3)*b(ix^S,2)
omega(ix^S,2)=j(ix^S,3)*b(ix^S,1)-j(ix^S,1)*b(ix^S,3)
omega(ix^S,3)=j(ix^S,1)*b(ix^S,2)-j(ix^S,2)*b(ix^S,1)
! omega = ((curl b) x b - (div b) b) / |b|^2
where(b2(ix^S)/=0.d0)
  omega(ix^S,1)=(omega(ix^S,1)-divb(ix^S)*b(ix^S,1))/b2(ix^S)
  omega(ix^S,2)=(omega(ix^S,2)-divb(ix^S)*b(ix^S,2))/b2(ix^S)
  omega(ix^S,3)=(omega(ix^S,3)-divb(ix^S)*b(ix^S,3))/b2(ix^S)
elsewhere
  omega(ix^S,1)=0.d0
  omega(ix^S,2)=0.d0
  omega(ix^S,3)=0.d0
endwhere
! omega2 = |omega|^2
omega2(ixO^S)=(^C&omega(ixO^S,^C)**2+)

! l_grid: the integral of b2*omega2 over the volume of grid 
l_grid=obj_funct_int(ixI^L,ixO^L,b2,omega2)
! ld_grid: the integral of |div B|^2 over the volume of grid
ld_grid=obj_funct_int(ixI^L,ixO^L,divb,divb)
! update b
if(evolve) then
  ! fun = - omega x j
  fun(ixO^S,1)=-(omega(ixO^S,2)*j(ixO^S,3)-omega(ixO^S,3)*j(ixO^S,2))
  fun(ixO^S,2)=-(omega(ixO^S,3)*j(ixO^S,1)-omega(ixO^S,1)*j(ixO^S,3))
  fun(ixO^S,3)=-(omega(ixO^S,1)*j(ixO^S,2)-omega(ixO^S,2)*j(ixO^S,1))
  ! tmpvector = omega x b
  tmpvector(ix^S,1)=omega(ix^S,2)*b(ix^S,3)-omega(ix^S,3)*b(ix^S,2)
  tmpvector(ix^S,2)=omega(ix^S,3)*b(ix^S,1)-omega(ix^S,1)*b(ix^S,3)
  tmpvector(ix^S,3)=omega(ix^S,1)*b(ix^S,2)-omega(ix^S,2)*b(ix^S,1)
  ! b2 = omega dot b
  b2(ix^S)=(^C&omega(ix^S,^C)*b(ix^S,^C)+)
  ix^L=ixI^L^LSUB4;
  ! j = grad (omega dot b)
  do idir=1,ndir
    call gradientS(b2,ixI^L,ix^L,idir,tmp)
    j(ixO^S,idir)=tmp(ixO^S)
  enddo
  ! fun = - omega x j - grad (omega dot b)
  fun(ixO^S,1)=fun(ixO^S,1)-j(ixO^S,1)
  fun(ixO^S,2)=fun(ixO^S,2)-j(ixO^S,2)
  fun(ixO^S,3)=fun(ixO^S,3)-j(ixO^S,3)
  ! j = curl (omega x b)
  call curlvector(tmpvector,ixI^L,ix^L,j,idirmin,1,ndir)
  ! fun = curl (omega x b) - omega x j - grad (omega dot b) + omega divb + |omega|^2 b
  fun(ixO^S,1)=fun(ixO^S,1)+j(ixO^S,1)+omega(ixO^S,1)*divb(ixO^S)+omega2(ixO^S)*b(ixO^S,1)
  fun(ixO^S,2)=fun(ixO^S,2)+j(ixO^S,2)+omega(ixO^S,2)*divb(ixO^S)+omega2(ixO^S)*b(ixO^S,2)
  fun(ixO^S,3)=fun(ixO^S,3)+j(ixO^S,3)+omega(ixO^S,3)*divb(ixO^S)+omega2(ixO^S)*b(ixO^S,3)
  ! b = b + dtfff*fun
  w(ixO^S,b0_+1:b0_+ndir)=b(ixO^S,1:ndir)+dtfff*fun(ixO^S,1:ndir)
end if

end subroutine fff_grid
!============================================================================= 
function obj_funct_int(ixI^L,ixO^L,obj1,obj2)
! PURPOSE: 
! obj_funct_int: the integral of obj1*obj2 over the volume; 
use fff_global

include 'amrvacdef.f'

integer, intent(in) :: ixI^L, ixO^L
double precision, intent(in) :: obj1(ixG^T),obj2(ixG^T)

integer :: ix^D
double precision :: obj_funct_int,dvolume
!-----------------------------------------------------------------------------
obj_funct_int=0.d0
dvolume=dxlevel(1)*dxlevel(2)*dxlevel(3)
{do ix^DB=ixOmin^DB,ixOmax^DB\}
   obj_funct_int=obj_funct_int+obj1(ix^D)*obj2(ix^D)*dvolume
{end do\}
end function obj_funct_int
!============================================================================= 
