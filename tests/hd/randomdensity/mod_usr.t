!=============================================================================
! Simulate 2D acoustic wave in random inhomogenueous medium
!=============================================================================
module mod_usr
  use mod_hd

  double precision,allocatable, save:: xgrid(:),ygrid(:)
  integer, save:: Nx,Ny,ntmp
  double precision,allocatable, save:: rand(:,:)

contains

  subroutine usr_init()

    call set_coordinate_system("Cartesian_2D")

    usr_init_one_grid => initonegrid_usr
    usr_set_parameters => initglobaldata_usr
    usr_aux_output      => specialvar_output
    usr_add_aux_names   => specialvarnames_output
    usr_var_for_errest  => myvar_for_errest
    usr_refine_grid     => specialrefine_grid


    call hd_activate()
  end subroutine usr_init

  subroutine initglobaldata_usr

integer,dimension(:),allocatable:: seed
integer::  seed_size,ix

double precision :: Lx, Ly


   hd_gamma=5.0d0/3.0d0
   Ly=xprobmax2-xprobmin2
   Lx=xprobmax1-xprobmin1
   Nx=nint(Lx/dx(1,1))
   if (.not. allocated(xgrid)) allocate(xgrid(Nx))
   do ix=1,Nx
      xgrid(ix)=xprobmin1+ix*dx(1,1)
   enddo
   Ny=nint(Ly/dx(2,1))
   if (.not. allocated(ygrid)) allocate(ygrid(Ny))
   do ix=1,Ny
      ygrid(ix)=xprobmin2+ix*dx(2,1)
   enddo
   ntmp=Nx*Ny
   if (.not. allocated(rand)) allocate(rand(Nx,Ny))

   if(mype==0)then
      print *,'Number of modes (x,y)=',Nx,Ny
      call random_seed(SIZE=seed_size)
      allocate(seed(seed_size))
      call random_seed(GET=seed(1:seed_size))
      call random_number(rand) ! generate ntmp random numbers [0,1]
   endif
   call MPI_BARRIER(icomm,ierrmpi)
   if(npe>1) call MPI_BCAST(rand,ntmp,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)

  end subroutine initglobaldata_usr

  ! initialize one grid
  subroutine initonegrid_usr(ixG^L,ix^L,w,x)
    use mod_dust
    integer, intent(in)             :: ixG^L, ix^L
    double precision, intent(in)    :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

integer:: i,ix^D
double precision :: rho(ixG^S)
double precision :: rand_vec(1:Ny)
double precision :: tmp,delrho,rand_avg,rand_std,A0,x0,y0,Sigx,Sigy
logical, save :: first=.true.


if(first)then
  if(mype==0) then
     write(*,*)'Simulate 2D acoustic wave in random inhomogenueous medium!'
  endif
  first=.false.
endif

delrho=0.2
rand_avg=SUM(rand)/DBLE(ntmp)
rand_std=SQRT(SUM((rand-rand_avg)**2)/DBLE(Ntmp-1))
rho(ix^S)=one
do ix1=ixmin1,ixmax1
do ix2=ixmin2,ixmax2
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
w(ix^S,mom(1))=zero
w(ix^S,mom(2))=zero
! Point-source pertubation
A0=0.2d0
x0=0.d0
y0=0.d0 ! start position of the pulse
Sigx=5d0
Sigy=5d0 ! x-and y-width of the pulse
w(ix^S,p_)=one*(1+A0*DEXP(-(x(ix^S,1)-x0)**2/Sigx**2)*DEXP(-(x(ix^S,2)-y0)**2/Sigy**2))
! Uncomment the following two lines for plane wave study

!x0=xprobmin1+0.1*Lx // position to launch a wave
!w(ix^S,p_)=one*(1+A0*DEXP(-(x(ix^S,1)-x0)**2/sigx**2))


call hd_to_conserved(ixG^L,ix^L,w,x)

  end subroutine initonegrid_usr

  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)

    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)

    double precision :: pth(ixI^S),wlocal(ixI^S,1:nw)
    double precision :: rho(ixI^S),gradrho(ixI^S),drho(ixI^S)
    double precision :: kk,kk0,grhomax,kk1
    integer :: idims

! Example: assuming nwauxio=3 at convert stage 

! first store temperature
    wlocal(ixI^S,1:nw)=w(ixI^S,1:nw)
    call hd_get_pthermal(wlocal,x,ixI^L,ixO^L,pth)
    w(ixO^S,nw+1)=pth(ixO^S)/w(ixO^S,rho_)

! then compute schlieren density plot
rho(ixI^S)=wlocal(ixI^S,rho_)
gradrho(ixO^S)=zero
do idims=1,ndim
   select case(typegrad)
   case("central")
     call gradient(rho,ixI^L,ixO^L,idims,drho)
   case("limited")
     call gradientS(rho,ixI^L,ixO^L,idims,drho)
   end select
   gradrho(ixO^S)=gradrho(ixO^S)+drho(ixO^S)**2.0d0
enddo
gradrho(ixO^S)=dsqrt(gradrho(ixO^S))
kk=5.0d0
kk0=0.001d0
kk1=1.0d0
grhomax=2000.0d0

w(ixO^S,nw+2)=dexp(-kk*(gradrho(ixO^S)-kk0*grhomax)/(kk1*grhomax-kk0*grhomax))

w(ixO^S,nw+3)=dlog10(w(ixO^S,rho_))


  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
    character(len=*) :: varnames

    varnames='T Schlier logrho'
  end subroutine specialvarnames_output

subroutine interp1d(nin,xin,yin,xout1,yout1)
! subroutine to interpolate values from a 1D function
  integer,intent(in) :: nin
  DOUBLE PRECISION,intent(in) :: xin(1:nin),yin(1:nin)
  DOUBLE PRECISION,intent(in) :: xout1
  DOUBLE PRECISION,intent(out) :: yout1
  integer*4 :: ixa
  DOUBLE PRECISION :: xa1,xb1,ya1,yb1,xtmp

  ixa=minloc(abs(xin-xout1),1)
  xtmp=xin(ixa)
  ! not in range of xin[0] and xin[nin]
  ! use extraplotion
  if (ixa == 1) then
     xa1=xin(1)
     xb1=xin(2)
     ya1=yin(1)
     yb1=yin(2)
     yout1=ya1+(yb1-ya1)*(xout1-xa1)/(xb1-xa1)
  endif
  if (ixa == nin) then
     xa1=xin(ixa-1)
     xb1=xin(ixa)
     ya1=yin(ixa-1)
     yb1=yin(ixa)
     yout1=ya1+(yb1-ya1)*(xout1-xa1)/(xb1-xa1)
  endif

  if (ixa < nin .and. ixa > 1) then
    if (xtmp > xout1) then
      ixa=ixa-1
    endif
    xa1=xin(ixa)
    xb1=xin(ixa+1)
    ya1=yin(ixa)
    yb1=yin(ixa+1)
    yout1=ya1+(yb1-ya1)*(xout1-xa1)/(xb1-xa1)
  endif

end subroutine interp1d

  subroutine specialrefine_grid(igrid,level,ixG^L,ix^L,qt,w,x,refine,coarsen)
    ! Enforce additional refinement or coarsening
    ! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.
    ! you must set consistent values for integers refine/coarsen:
    ! refine = -1 enforce to not refine
    ! refine =  0 doesn't enforce anything
    ! refine =  1 enforce refinement
    ! coarsen = -1 enforce to not coarsen
    ! coarsen =  0 doesn't enforce anything
    ! coarsen =  1 enforce coarsen
    integer, intent(in) :: igrid, level, ixG^L, ix^L
    double precision, intent(in) :: qt, w(ixG^S,1:nw), x(ixG^S,1:ndim)
    integer, intent(inout) :: refine, coarsen

    double precision:: rad(ixG^T)

    rad(ix^S) = dsqrt(x(ix^S,1)**2 + x(ix^S,2)**2)
    if( qt<10.d0.and.any(rad(ix^S) <= 7.d0) ) then
      coarsen = -1
      refine  = 1
    endif
    if( qt<10.0d0.and.all(rad(ix^S) > 15.d0) ) then
      coarsen = 1
      refine  = -1
    endif

  end subroutine specialrefine_grid

  subroutine myvar_for_errest(ixI^L,ixO^L,iflag,w,x,var)
      use mod_global_parameters
      integer, intent(in)           :: ixI^L,ixO^L,iflag
      double precision, intent(in)  :: w(ixI^S,1:nw), x(ixI^S,1:ndim)
      double precision, intent(out) :: var(ixI^S)

      if (iflag >nw+1)call mpistop(' iflag error')
      call hd_get_pthermal(w,x,ixI^L,ixO^L,var)

  end subroutine myvar_for_errest


end module mod_usr
