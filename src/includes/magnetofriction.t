!> module magnetofriction.t
!> Purpose: use magnetofrictional method to relax 3D magnetic field to 
!>          force-free field
!> 01.04.2016 developed by Chun Xia and Yang Guo
!> Usage:
!> 1. in definitions.h:
!>    #define MAGNETOFRICTION
!> 2. in amrvac.par:
!>    &stoplist
!>     itmaxmf=60000 ! set the maximum iteration number
!>    &savelist
!>     ditsavemf=20000 ! set iteration interval for data output
!>    &methodlist
!>     typeadvance='onestep' ! time marching scheme, or 'twostep','threestep'
!>     typefull1=13*'cd4' ! or 'tvdlf', 'fd'
!>     typelimiter1= 13*'koren' ! or 'vanleer','cada3','mp5' so on
!>    &amrlist
!>     ditregrid=20 ! set iteration interval for adjusting AMR 
!>    &paramlist
!>     cmf_c=0.3    ! stability coefficient controls numerical stability
!>     cmf_y=0.2    ! frictional velocity coefficient
!>     cmf_divb=0.1 ! divb cleaning coefficient controls diffusion speed of divb
module mod_magnetofriction
  implicit none

  double precision :: cmf_c
  double precision :: cmf_y
  double precision :: cmf_divb
  double precision :: tmf
  integer :: ditsavemf
  integer :: itmaxmf
  logical :: mf_advance
subroutine magnetofriction

use mod_global_parameters
use mod_input_output

double precision :: dvolume(ixG^T),dsurface(ixG^T),dvone
double precision :: dtfff,dtfff_pe,dtnew,dx^D
double precision :: cwsin_theta_new,cwsin_theta_old
double precision :: sum_jbb,sum_jbb_ipe,sum_j,sum_j_ipe,sum_l_ipe,sum_l
double precision :: f_i_ipe,f_i,volumepe,volume,tmpt,time_in
double precision, external :: integral_grid
integer :: i,iigrid, igrid, idims,ix^D,hxM^LL,fhmf,tmpit,i^D
logical :: patchwi(ixG^T)
!-----------------------------------------------------------------------------
time_in=MPI_WTIME()
if(mype==0) write(*,*) 'Evolving to force-free field using magnetofricitonal method...'
if(prolongprimitive) call mpistop('use prolongprimitive=.false. in MF module')
mf_advance=.false.
dtfff=1.d-2
tmpt=t
tmpit=it
tmf=t
i=it
! update ghost cells
call getbc(tmf,0.d0,pw,0,nwflux)
if(snapshotini==-1 .and. i==0) then
  call saveamrfile(1)
  call saveamrfile(2)
end if
mf_advance=.true.
! convert conservative variables to primitive ones which are used during MF
do iigrid=1,igridstail; igrid=igrids(iigrid);
   call primitive(ixG^LL,ixM^LL,pw(igrid)%w,px(igrid)%x)
end do
! calculate magnetofrictional velocity
call mf_velocity_update(pw,dtfff)
! update velocity in ghost cells
bcphys=.false.
call getbc(tmf,0.d0,pw,v0_,ndir)
bcphys=.true.
! calculate initial values of metrics
if(i==0) then
  call metrics
  call printlog_mf 
end if
! magnetofrictional loops
do
  ! calculate time step based on Cmax= Alfven speed + abs(frictional speed)
  dtfff_pe=bigdouble
  cmax_mype=zero
  do iigrid=1,igridstail; igrid=igrids(iigrid);
    pwold(igrid)%w(ixG^T,b0_+1:b0_+ndir)=pw(igrid)%w(ixG^T,b0_+1:b0_+ndir)
    if (.not.slab) mygeo => pgeo(igrid)
    if (B0field) then
       myB0_cell => pB0_cell(igrid)
       {^D&myB0_face^D => pB0_face^D(igrid)\}
    end if
    typelimiter=typelimiter1(node(plevel_,igrid))
    typegradlimiter=typegradlimiter1(node(plevel_,igrid))
    ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
    call getdtfff_courant(pw(igrid)%w,px(igrid)%x,ixG^LL,ixM^LL,dtnew)
    dtfff_pe=min(dtfff_pe,dtnew)
  end do
  call MPI_ALLREDUCE(dtfff_pe,dtfff,1,MPI_DOUBLE_PRECISION,MPI_MIN, &
                     icomm,ierrmpi)
  call MPI_ALLREDUCE(cmax_mype,cmax_global,1,MPI_DOUBLE_PRECISION,MPI_MAX,&
                     icomm,ierrmpi)

  ! =======
  ! evolve
  ! =======
  call advectmf(1,ndim,tmf,dtfff)

  ! clean divergence of magnetic field 
  do iigrid=1,igridstail; igrid=igrids(iigrid);
    if (.not.slab) mygeo => pgeo(igrid)
    if (B0field) then
       myB0_cell => pB0_cell(igrid)
       {^D&myB0_face^D => pB0_face^D(igrid)\}
    end if
    typelimiter=typelimiter1(node(plevel_,igrid))
    typegradlimiter=typegradlimiter1(node(plevel_,igrid))
    ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
    {do i^DB=-1,1\}
       if (i^D==0|.and.) cycle
       if (neighbor_type(i^D,igrid)==2 .or. neighbor_type(i^D,igrid)==4) then
          leveljump(i^D)=.true.
       else
          leveljump(i^D)=.false.
       end if
    {end do\}
    call divbclean(ixG^LL,ixM^LL,pw(igrid)%w,px(igrid)%x,dtfff)
  end do
  ! update B in ghost cells
  call getbc(tmf+dtfff,dtfff,pw,b0_,ndir)
  ! calculate magnetofrictional velocity
  call mf_velocity_update(pw,dtfff)
  ! update velocity in ghost cells
  bcphys=.false.
  call getbc(tmf+dtfff,dtfff,pw,v0_,ndir)
  bcphys=.true.

  i=i+1
  tmf=tmf+dtfff
  if(mod(i,10)==0) then
    ! calculate metrics
    call metrics
    call printlog_mf
  end if
  if(mod(i,ditsavemf)==0) then
    it=i
    t=tmf
    do iigrid=1,igridstail; igrid=igrids(iigrid);
      call conserve(ixG^LL,ixG^LL,pw(igrid)%w,px(igrid)%x,patchfalse)
    end do
    mf_advance=.false.
    call saveamrfile(1)
    call saveamrfile(2)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       call primitive(ixG^LL,ixG^LL,pw(igrid)%w,px(igrid)%x)
    end do
    mf_advance=.true.
    if(mype==0) then
      write(*,*) "itmf=",i
      write(*,*) '<CW sin theta>:',cwsin_theta_new
      write(*,*) '<f_i>:',f_i
      write(*,*) '----------------------------------------------------------'
    end if
  end if
  ! reconstruct AMR grid every 10 step
  if(mod(i,ditregrid)==0 .and. mxnest>1) call resettree
  if (i>=itmaxmf) then
    if(mod(i,10)/=0) then
      ! calculate metrics
      call metrics
      call printlog_mf
    end if
    if(mype==0) then
      write (*,*) 'The magnetofrictional iteration has been terminated because &
                     it reaches the maximum iteration step!'
      write (*,*) 'The total iteration step is:', i   
    end if
    exit
  end if
enddo
! set velocity back to zero and convert primitive variables back to conservative ones
do iigrid=1,igridstail; igrid=igrids(iigrid);
   pw(igrid)%w(ixG^T,v0_+1:v0_+ndir)=zero
   call conserve(ixG^LL,ixM^LL,pw(igrid)%w,px(igrid)%x,patchfalse)
end do
t=tmpt
it=tmpit
if (mype==0) call MPI_FILE_CLOSE(fhmf,ierrmpi)
mf_advance=.false.
if(mype==0) write(*,*) 'Magnetofriction phase took : ',MPI_WTIME()-time_in,' sec'
contains
!=============================================================================
! internal procedures start
!=============================================================================
subroutine metrics

sum_jbb_ipe = 0.d0
sum_j_ipe = 0.d0
sum_l_ipe = 0.d0
f_i_ipe = 0.d0
volumepe=0.d0
do iigrid=1,igridstail; igrid=igrids(iigrid);
  if (.not.slab) mygeo => pgeo(igrid)
  if (B0field) then
     myB0_cell => pB0_cell(igrid)
     {^D&myB0_face^D => pB0_face^D(igrid)\}
  end if
  if(slab) then
    dvone={rnode(rpdx^D_,igrid)|*}
    dvolume(ixM^T)=dvone
    dsurface(ixM^T)=two*(^D&dvone/rnode(rpdx^D_,igrid)+)
  else
    dvolume(ixM^T)=pgeo(igrid)%dvolume(ixM^T)
    dsurface(ixM^T)= ^D&pgeo(igrid)%surfaceC^D(ixM^T)+
    do idims=1,ndim
      hxM^LL=ixM^LL-kr(idims,^D);
      select case(idims)
      {case(^D)
         dsurface(ixM^T)=dsurface(ixM^T)+pgeo(igrid)%surfaceC^D(hxM^T) \}
      end select
    end do
  end if
  ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
  call mask_inner(ixG^LL,ixM^LL,pw(igrid)%w,px(igrid)%x)
  sum_jbb_ipe = sum_jbb_ipe+integral_grid_mf(ixG^LL,ixM^LL,pw(igrid)%w,&
    px(igrid)%x,1,patchwi)
  sum_j_ipe   = sum_j_ipe+integral_grid_mf(ixG^LL,ixM^LL,pw(igrid)%w,&
    px(igrid)%x,2,patchwi)
  f_i_ipe=f_i_ipe+integral_grid_mf(ixG^LL,ixM^LL,pw(igrid)%w,&
    px(igrid)%x,3,patchwi)
  sum_l_ipe   = sum_l_ipe+integral_grid_mf(ixG^LL,ixM^LL,pw(igrid)%w,&
    px(igrid)%x,4,patchwi)
end do
call MPI_ALLREDUCE(sum_jbb_ipe,sum_jbb,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
   icomm,ierrmpi)
call MPI_ALLREDUCE(sum_j_ipe,sum_j,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
   icomm,ierrmpi)
call MPI_ALLREDUCE(f_i_ipe,f_i,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
   icomm,ierrmpi)
call MPI_ALLREDUCE(sum_l_ipe,sum_l,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
   icomm,ierrmpi)
call MPI_ALLREDUCE(volumepe,volume,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
   icomm,ierrmpi)
! sin of the j weighted mean angle between current and magnetic field
cwsin_theta_new = sum_jbb/sum_j
! mean normalized divergence of magnetic field
f_i = f_i/volume
sum_j=sum_j/volume
sum_l=sum_l/volume

end subroutine metrics
!=============================================================================
subroutine mask_inner(ixI^L,ixO^L,w,x)

integer, intent(in)         :: ixI^L,ixO^L
double precision, intent(in):: w(ixI^S,nw),x(ixI^S,1:ndim)
double precision            :: xO^L
integer                     :: ix^D
!-----------------------------------------------------------------------------
if(slab) then
  xOmin1 = xprobmin1 + 0.05d0*(xprobmax1-xprobmin1)
  xOmax1 = xprobmax1 - 0.05d0*(xprobmax1-xprobmin1)
  xOmin2 = xprobmin2 + 0.05d0*(xprobmax2-xprobmin2)
  xOmax2 = xprobmax2 - 0.05d0*(xprobmax2-xprobmin2)
  xOmin3 = xprobmin3
  xOmax3 = xprobmax3 - 0.05d0*(xprobmax3-xprobmin3)
else
  xOmin1 = xprobmin1
  xOmax1 = xprobmax1 - 0.05d0*(xprobmax1-xprobmin1)
  xOmin2 = xprobmin2 + 0.05d0*(xprobmax2-xprobmin2)
  xOmax2 = xprobmax2 - 0.05d0*(xprobmax2-xprobmin2)
  xOmin3 = xprobmin3 + 0.05d0*(xprobmax3-xprobmin3)
  xOmax3 = xprobmax3 - 0.05d0*(xprobmax3-xprobmin3)
end if

{do ix^DB=ixOmin^DB,ixOmax^DB\}
    if(x(ix^D,1) > xOmin1 .and. x(ix^D,1) < xOmax1 .and. &
       x(ix^D,2) > xOmin2 .and. x(ix^D,2) < xOmax2 .and. &
       x(ix^D,3) > xOmin3 .and. x(ix^D,3) < xOmax3) then
      patchwi(ix^D)=.true.
      volumepe=volumepe+dvolume(ix^D)
    else
      patchwi(ix^D)=.false.
    endif
{end do\}

end subroutine mask_inner
!============================================================================= 
subroutine printlog_mf
integer :: amode, status(MPI_STATUS_SIZE)
character(len=800) :: filename,filehead
character(len=2048) :: line,datastr
logical, save :: logmfopened=.false.
!-----------------------------------------------------------------------------
if(mype==0) then
  if(.not.logmfopened) then
    ! generate filename
    write(filename,"(a,a)") TRIM(base_filename), "_mflog.csv"

    amode=ior(MPI_MODE_CREATE,MPI_MODE_WRONLY)
    amode=ior(amode,MPI_MODE_APPEND)
    call MPI_FILE_OPEN(MPI_COMM_SELF,filename,amode,MPI_INFO_NULL,fhmf,ierrmpi)
    logmfopened=.true.
    filehead="itmf,<f_i>,<CW sin theta>,<Current>,<Lorenz force>"
    call MPI_FILE_WRITE(fhmf,filehead,len_trim(filehead), &
                        MPI_CHARACTER,status,ierrmpi)
    call MPI_FILE_WRITE(fhmf,achar(10),1,MPI_CHARACTER,status,ierrmpi)
  end if
  line=''
  write(datastr,'(i6,a)') i,','
  line=trim(line)//trim(datastr)
  write(datastr,'(es13.6,a)') f_i,','
  line=trim(line)//trim(datastr)
  write(datastr,'(es13.6,a)') cwsin_theta_new,','
  line=trim(line)//trim(datastr)
  write(datastr,'(es13.6,a)') sum_j,','
  line=trim(line)//trim(datastr)
  write(datastr,'(es13.6)') sum_l
  line=trim(line)//trim(datastr)//new_line('A')
  call MPI_FILE_WRITE(fhmf,line,len_trim(line),MPI_CHARACTER,status,ierrmpi)
end if

end subroutine printlog_mf
!=============================================================================
function integral_grid_mf(ixI^L,ixO^L,w,x,iw,patchwi)

integer, intent(in)                :: ixI^L,ixO^L,iw
double precision, intent(in)       :: x(ixI^S,1:ndim)
double precision                   :: w(ixI^S,nw+nwauxio)
logical, intent(in) :: patchwi(ixI^S)

double precision, dimension(ixI^S,1:ndir) :: bvec,qvec,current
double precision :: integral_grid_mf,tmp(ixI^S),b_mag(ixI^S)
integer :: ix^D,i,idirmin,idir,jdir,kdir
!-----------------------------------------------------------------------------
integral_grid_mf=0.d0
select case(iw)
 case(1)
  ! Sum(|JxB|/|B|)
  if(B0field) then
    ^C&bvec(ixI^S,^C)=w(ixI^S,b^C_)+myB0_cell%w(ixI^S,^C);
  else
    ^C&bvec(ixI^S,^C)=w(ixI^S,b^C_);
  endif
  call curlvector(bvec,ixI^L,ixO^L,current,idirmin,1,ndir)
  ! calculate Lorentz force
  qvec(ixO^S,1:ndir)=zero
  do idir=1,ndir; do jdir=1,ndir; do kdir=idirmin,3
     if(lvc(idir,jdir,kdir)/=0)then
        tmp(ixO^S)=current(ixO^S,jdir)*bvec(ixO^S,kdir)
        if(lvc(idir,jdir,kdir)==1)then
           qvec(ixO^S,idir)=qvec(ixO^S,idir)+tmp(ixO^S)
        else
           qvec(ixO^S,idir)=qvec(ixO^S,idir)-tmp(ixO^S)
        endif
     endif
  enddo; enddo; enddo
  
  {do ix^DB=ixOmin^DB,ixOmax^DB\}
     if(patchwi(ix^D)) integral_grid_mf=integral_grid_mf+dsqrt((^C&qvec(ix^D,^C)**2+)/&
                       (^C&bvec(ix^D,^C)**2+))*dvolume(ix^D)
  {end do\}
 case(2)
  ! Sum(|J|)
  if(B0field) then
    ^C&bvec(ixI^S,^C)=w(ixI^S,b^C_)+myB0_cell%w(ixI^S,^C);
  else
    ^C&bvec(ixI^S,^C)=w(ixI^S,b^C_);
  endif
  call curlvector(bvec,ixI^L,ixO^L,current,idirmin,1,ndir)
  {do ix^DB=ixOmin^DB,ixOmax^DB\}
     if(patchwi(ix^D)) integral_grid_mf=integral_grid_mf+dsqrt(^C&current(ix^D,^C)**2+)*&
                       dvolume(ix^D)
  {end do\}
 case(3)
  ! f_i solenoidal property of B: (dvolume |div B|)/(dsurface |B|)
  if(B0field) then
    ^C&bvec(ixI^S,^C)=w(ixI^S,b^C_)+myB0_cell%w(ixI^S,^C);
  else
    ^C&bvec(ixI^S,^C)=w(ixI^S,b^C_);
  endif
  call divvector(bvec,ixI^L,ixO^L,tmp)
  {do ix^DB=ixOmin^DB,ixOmax^DB\}
     if(patchwi(ix^D)) integral_grid_mf=integral_grid_mf+dabs(tmp(ix^D))*&
           dvolume(ix^D)**2/dsqrt(^C&bvec(ix^D,^C)**2+)/dsurface(ix^D)
  {end do\}
 case(4)
  ! Sum(|JxB|)
  if(B0field) then
    ^C&bvec(ixI^S,^C)=w(ixI^S,b^C_)+myB0_cell%w(ixI^S,^C);
  else
    ^C&bvec(ixI^S,^C)=w(ixI^S,b^C_);
  endif
  call curlvector(bvec,ixI^L,ixO^L,current,idirmin,1,ndir)
  ! calculate Lorentz force
  qvec(ixO^S,1:ndir)=zero
  do idir=1,ndir; do jdir=1,ndir; do kdir=idirmin,3
     if(lvc(idir,jdir,kdir)/=0)then
        tmp(ixO^S)=current(ixO^S,jdir)*bvec(ixO^S,kdir)
        if(lvc(idir,jdir,kdir)==1)then
           qvec(ixO^S,idir)=qvec(ixO^S,idir)+tmp(ixO^S)
        else
           qvec(ixO^S,idir)=qvec(ixO^S,idir)-tmp(ixO^S)
        endif
     endif
  enddo; enddo; enddo
  
  {do ix^DB=ixOmin^DB,ixOmax^DB\}
     if(patchwi(ix^D)) integral_grid_mf=integral_grid_mf+dsqrt(^C&qvec(ix^D,^C)**2+)*dvolume(ix^D)
  {end do\}
end select
return
end function integral_grid_mf
!=============================================================================
! internal procedures end
!=============================================================================
end subroutine magnetofriction
!=============================================================================
subroutine mf_velocity_update(pwa,dtfff)

use mod_global_parameters

double precision, intent(in) :: dtfff 
integer :: i,iigrid, igrid
type(walloc) :: pwa(ngridshi)
double precision :: vhatmax,vhatmax_pe,vhatmaxgrid
!-----------------------------------------------------------------------------
vhatmax_pe=smalldouble
do iigrid=1,igridstail; igrid=igrids(iigrid);
  if (.not.slab) mygeo => pgeo(igrid)
  if (B0field) then
     myB0_cell => pB0_cell(igrid)
     {^D&myB0_face^D => pB0_face^D(igrid)\}
  end if
  ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
  call vhat(pwa(igrid)%w,px(igrid)%x,ixG^LL,ixM^LL,vhatmaxgrid)
  vhatmax_pe=max(vhatmax_pe,vhatmaxgrid)
end do
call MPI_ALLREDUCE(vhatmax_pe,vhatmax,1,MPI_DOUBLE_PRECISION,MPI_MAX, &
                       icomm,ierrmpi)
do iigrid=1,igridstail; igrid=igrids(iigrid);
  if (.not.slab) mygeo => pgeo(igrid)
  if (B0field) then
     myB0_cell => pB0_cell(igrid)
     {^D&myB0_face^D => pB0_face^D(igrid)\}
  end if
  typelimiter=typelimiter1(node(plevel_,igrid))
  typegradlimiter=typegradlimiter1(node(plevel_,igrid))
  ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
  ! calculate frictional velocity
  call frictional_velocity(pwa(igrid)%w,px(igrid)%x,ixG^LL,ixM^LL,vhatmax,dtfff)
end do

end subroutine mf_velocity_update
!============================================================================= 
subroutine vhat(w,x,ixI^L,ixO^L,vhatmaxgrid)

! Calculate v_hat 

use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L
double precision, intent(inout)  :: w(ixI^S,nw)
double precision, intent(in)  :: x(ixI^S,1:ndim)
double precision, intent(out) :: vhatmaxgrid

double precision              :: current(ixI^S,7-2*ndir:3),tmp(ixI^S),dxhm
integer :: idirmin,idir,jdir,kdir
!-----------------------------------------------------------------------------

call getcurrent(w,ixI^L,ixO^L,idirmin,current)
w(ixI^S,v0_+1:v0_+ndir)=0.d0
! calculate Lorentz force
do idir=1,ndir; do jdir=1,ndir; do kdir=idirmin,3
   if(lvc(idir,jdir,kdir)/=0)then
      if(B0field) then
        tmp(ixO^S)=current(ixO^S,jdir)*(w(ixO^S,b0_+kdir)+myB0_cell%w(ixO^S,kdir))
      else
        tmp(ixO^S)=current(ixO^S,jdir)*w(ixO^S,b0_+kdir)
      endif
      if(lvc(idir,jdir,kdir)==1)then
         w(ixO^S,v0_+idir)=w(ixO^S,v0_+idir)+tmp(ixO^S)
      else
         w(ixO^S,v0_+idir)=w(ixO^S,v0_+idir)-tmp(ixO^S)
      endif
   endif
enddo; enddo; enddo

if(B0field) then
  tmp(ixO^S)=( ^C&(w(ixO^S,b^C_)+myB0_cell%w(ixO^S,^C))**2+ )         ! |B|**2
else
  tmp(ixO^S)=( ^C&w(ixO^S,b^C_)**2+ )         ! |B|**2
endif

dxhm=dble(ndim)/(^D&1.0d0/dxlevel(^D)+)
^C&w(ixO^S,v^C_)=dxhm*w(ixO^S,v^C_)/tmp(ixO^S);
! ^C&w(ixO^S,v^C_)=w(ixO^S,v^C_)/tmp(ixO^S);
vhatmaxgrid=maxval(dsqrt( ^C&w(ixO^S,v^C_)**2+ ))

end subroutine vhat
!============================================================================= 
subroutine frictional_velocity(w,x,ixI^L,ixO^L,qvmax,qdt)

use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L
double precision, intent(in) :: x(ixI^S,1:ndim),qdt,qvmax
double precision, intent(inout) :: w(ixI^S,1:nw)

double precision :: dxhm,disbd(5),bfzone^D
integer :: ix^D
logical :: buffer
!-----------------------------------------------------------------------------
dxhm=dble(ndim)/(^D&1.0d0/dxlevel(^D)+)
dxhm=cmf_c*cmf_y/qvmax*dxhm/qdt
! dxhm=cmf_c*cmf_y/qvmax
^C&w(ixO^S,v0_+^C)=w(ixO^S,v0_+^C)*dxhm;
buffer=.true.
if (buffer) then
  bfzone1=0.05d0*(xprobmax1-xprobmin1)
  bfzone2=0.05d0*(xprobmax2-xprobmin2)
  bfzone3=0.05d0*(xprobmax3-xprobmin3)
  if(slab) then
    {do ix^DB=ixOmin^DB,ixOmax^DB\}
       disbd(1)=x(ix^D,1)-(xprobmin1+0.5d0*dxlevel(1))
       disbd(2)=(xprobmax1-0.5d0*dxlevel(1))-x(ix^D,1)
       disbd(3)=x(ix^D,2)-(xprobmin2+0.5d0*dxlevel(2))
       disbd(4)=(xprobmax2-0.5d0*dxlevel(2))-x(ix^D,2)
       disbd(5)=(xprobmax3-0.5d0*dxlevel(3))-x(ix^D,3)

       if(disbd(1)<bfzone1) then
         w(ix^D,v1_:v3_)=(1.d0-((bfzone1-disbd(1))/bfzone1)**2)*w(ix^D,v1_:v3_)
       endif
       if(disbd(2)<bfzone1) then
         w(ix^D,v1_:v3_)=(1.d0-((bfzone1-disbd(2))/bfzone1)**2)*w(ix^D,v1_:v3_)
       endif
       if(disbd(3)<bfzone2) then
         w(ix^D,v1_:v3_)=(1.d0-((bfzone2-disbd(3))/bfzone2)**2)*w(ix^D,v1_:v3_)
       endif
       if(disbd(4)<bfzone2) then
         w(ix^D,v1_:v3_)=(1.d0-((bfzone2-disbd(4))/bfzone2)**2)*w(ix^D,v1_:v3_)
       endif
       if(disbd(5)<bfzone3) then
         w(ix^D,v1_:v3_)=(1.d0-((bfzone3-disbd(5))/bfzone3)**2)*w(ix^D,v1_:v3_)
       endif
    {end do\}
  else
    {do ix^DB=ixOmin^DB,ixOmax^DB\}
       disbd(2)=(xprobmax1-0.5d0*dxlevel(1))-x(ix^D,1)
       disbd(3)=x(ix^D,2)-(xprobmin2+0.5d0*dxlevel(2))
       disbd(4)=(xprobmax2-0.5d0*dxlevel(2))-x(ix^D,2)
       disbd(5)=(xprobmax3-0.5d0*dxlevel(3))-x(ix^D,3)
       disbd(1)=x(ix^D,3)-(xprobmin3+0.5d0*dxlevel(3))

       if(disbd(2)<bfzone1) then
         w(ix^D,v1_:v3_)=(1.d0-((bfzone1-disbd(2))/bfzone1)**2)*w(ix^D,v1_:v3_)
       endif
       if(disbd(3)<bfzone2) then
         w(ix^D,v1_:v3_)=(1.d0-((bfzone2-disbd(3))/bfzone2)**2)*w(ix^D,v1_:v3_)
       endif
       if(disbd(4)<bfzone2) then
         w(ix^D,v1_:v3_)=(1.d0-((bfzone2-disbd(4))/bfzone2)**2)*w(ix^D,v1_:v3_)
       endif
       if(disbd(5)<bfzone3) then
         w(ix^D,v1_:v3_)=(1.d0-((bfzone3-disbd(5))/bfzone3)**2)*w(ix^D,v1_:v3_)
       endif
       if(disbd(1)<bfzone3) then
         w(ix^D,v1_:v3_)=(1.d0-((bfzone3-disbd(1))/bfzone3)**2)*w(ix^D,v1_:v3_)
       endif
    {end do\}
  end if
end if
end subroutine frictional_velocity
!=============================================================================
subroutine advectmf(idim^LIM,qt,qdt)

!  integrate all grids by one step of its delta(t)

! This subroutine is in VAC terminology equivalent to
! `advect' (with the difference that it will `advect' all grids)

use mod_global_parameters

integer, intent(in) :: idim^LIM
double precision, intent(in) :: qt, qdt

integer :: iigrid, igrid
!-----------------------------------------------------------------------------
! copy w instead of wold because of potential use of dimsplit or sourcesplit
do iigrid=1,igridstail; igrid=igrids(iigrid);
   allocate (pw1(igrid)%w(ixG^T,1:nw))
   pw1(igrid)%w=pw(igrid)%w
end do

istep=0

select case (typeadvance)
 case ("onestep")
   call advect1mf(typefull1,qdt,one,    idim^LIM,qt,          pw1,qt,pw, pwold)
 case ("twostep")
   ! predictor step
   call advect1mf(typepred1,qdt,half,   idim^LIM,qt,          pw,qt,pw1,pwold)
   ! corrector step
   call advect1mf(typefull1,qdt,one,    idim^LIM,qt+half*qdt, pw1,qt,pw, pwold)
 case ("threestep")
   ! three step Runge-Kutta in accordance with Gottlieb & Shu 1998
   call advect1mf(typefull1,qdt,one,    idim^LIM,qt,          pw ,qt,pw1,pwold)

   do iigrid=1,igridstail; igrid=igrids(iigrid);
      allocate (pw2(igrid)%w(ixG^T,1:nw))
      pw2(igrid)%w(ixG^T,1:nwflux)=0.75d0*pw(igrid)%w(ixG^T,1:nwflux)+0.25d0*&
        pw1(igrid)%w(ixG^T,1:nwflux)
   end do

   call advect1mf(typefull1,qdt,0.25d0, idim^LIM,qt+qdt,pw1,qt+dt*0.25d0,pw2,pwold)

   do iigrid=1,igridstail; igrid=igrids(iigrid);
      pw(igrid)%w(ixG^T,1:nwflux)=1.0d0/3.0d0*pw(igrid)%w(ixG^T,1:nwflux)+&
        2.0d0/3.0d0*pw2(igrid)%w(ixG^T,1:nwflux)
   end do   
   call advect1mf(typefull1,qdt,2.0d0/3.0d0, idim^LIM,qt+qdt/2.0d0,pw2,&
          qt+qdt/3.0d0,pw,pwold)
 case default
   write(unitterm,*) "typeadvance=",typeadvance
   write(unitterm,*) "Error in advectmf: Unknown time integration method"
   call mpistop("Correct typeadvance")
end select

do iigrid=1,igridstail; igrid=igrids(iigrid);
   deallocate (pw1(igrid)%w)
   select case (typeadvance)
     case ("threestep")
       deallocate (pw2(igrid)%w)
   end select
end do

end subroutine advectmf
!=============================================================================
subroutine advect1mf(method,dtin,dtfactor,idim^LIM,qtC,pwa,qt,pwb,pwc)

! Integrate all grids by one partial step
! This subroutine is equivalent to VAC's `advect1', but does
! the advection for all grids
use mod_global_parameters

integer, intent(in) :: idim^LIM
double precision, intent(in) :: dtin,dtfactor, qtC, qt
character(len=*), intent(in) :: method(nlevelshi)
type(walloc) :: pwa(ngridshi), pwb(ngridshi), pwc(ngridshi)

double precision :: qdt
integer :: iigrid, igrid, level
logical :: setigrid
!-----------------------------------------------------------------------------
istep=istep+1

if (levmax>levmin) then
   if (istep==nstep.or.nstep>2) &
   call init_comm_fix_conserve(idim^LIM)
end if

! loop over all grids to arrive at equivalent

! opedit: Just advance the active grids: 
qdt=dtfactor*dtin
do iigrid=1,igridstail; igrid=igrids(iigrid);
   level=node(plevel_,igrid)
   call process1_gridmf(method(level),igrid,qdt,ixG^LL,idim^LIM,qtC,&
                   pwa(igrid)%w,qt,pwb(igrid)%w,pwc(igrid)%w)
end do

! opedit: Send flux for all grids, expects sends for all 
! nsend_fc(^D), set in connectivity.t.

if (levmax>levmin) then
  if (istep==nstep.or.nstep>2) then
     do iigrid=1,igridstail; igrid=igrids(iigrid);
        call sendflux(igrid,idim^LIM)
     end do
     call fix_conserve(pwb,idim^LIM)
  end if
end if
   
! update B in ghost cells
call getbc(qt+qdt,qdt,pwb,b0_,ndir)
! calculate magnetofrictional velocity
call mf_velocity_update(pwb,qdt)
! update magnetofrictional velocity in ghost cells
bcphys=.false.
call getbc(qt+qdt,qdt,pwb,v0_,ndir)
bcphys=.true.

end subroutine advect1mf
!=============================================================================
subroutine process1_gridmf(method,igrid,qdt,ixG^L,idim^LIM,qtC,wCT,qt,w,wold)

! This subroutine is equivalent to VAC's `advect1' for one grid
use mod_global_parameters

character(len=*), intent(in) :: method
integer, intent(in) :: igrid, ixG^L, idim^LIM
double precision, intent(in) :: qdt, qtC, qt
double precision :: wCT(ixG^S,1:nw), w(ixG^S,1:nw), wold(ixG^S,1:nw)
double precision :: dx^D, fC(ixG^S,1:nwflux,1:ndim)
integer :: ixO^L
!-----------------------------------------------------------------------------
dx^D=rnode(rpdx^D_,igrid);
^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
saveigrid=igrid
fC=0.d0

if (.not.slab) mygeo => pgeo(igrid)
if (B0field) then
   myB0_cell => pB0_cell(igrid)
   {^D&myB0_face^D => pB0_face^D(igrid)\}
end if
typelimiter=typelimiter1(node(plevel_,igrid))
typegradlimiter=typegradlimiter1(node(plevel_,igrid))

ixO^L=ixG^L^LSUBdixB;
select case (method)
 case ("cd4")
   !================================
   ! 4th order central difference
   !================================ 
   call centdiff4mf(qdt,ixG^L,ixO^L,idim^LIM,qtC,wCT,qt,w,wold,fC,dx^D,px(igrid)%x)
 case ("tvdlf")
   !================================
   ! TVDLF
   !================================ 
   call tvdlfmf(qdt,ixG^L,ixO^L,idim^LIM,qtC,wCT,qt,w,wold,fC,dx^D,px(igrid)%x)
 case ('hancock')
   ! hancock predict (first) step for twostep tvdlf and tvdmu scheme
   call hancockmf(qdt,ixG^L,ixO^L,idim^LIM,qtC,wCT,qt,w,dx^D,px(igrid)%x)
 case ("fd")
   !================================
   ! finite difference
   !================================ 
   call fdmf(qdt,ixG^L,ixO^L,idim^LIM,qtC,wCT,qt,w,wold,fC,dx^D,px(igrid)%x)
case default
   if(mype==0) write(unitterm,*)'Error in advect1_gridmf:',method,' is not there!'
   call mpistop("The scheme is not implemented.")
end select

if (levmax>levmin) then
  if (istep==nstep.or.nstep>2) &
    call storeflux(igrid,fC,idim^LIM)
end if

end subroutine process1_gridmf
!============================================================================
subroutine upwindLRmf(ixI^L,ixL^L,ixR^L,idims,w,wCT,wLC,wRC,x,dxdim)

! Determine the upwinded wLC(ixL) and wRC(ixR) from w. 
! the wCT is only used when PPM is exploited.

use mod_global_parameters

integer, intent(in) :: ixI^L, ixL^L, ixR^L, idims
double precision, intent(in) :: dxdim
double precision, dimension(ixI^S,1:nw) :: w, wCT
double precision, dimension(ixI^S,1:nw) :: wLC, wRC
double precision, dimension(ixI^S,1:ndim) :: x

integer :: jxR^L, ixC^L, jxC^L, iw, ixtest^L
double precision :: wLtmp(ixI^S,1:nw), wRtmp(ixI^S,1:nw)
double precision :: ldw(ixI^S), dwC(ixI^S)
logical, dimension(ixI^S) :: flagL, flagR
logical, dimension(ixI^S) :: patchw, patchwLC, patchwRC

character*79 :: savetypelimiter
!-----------------------------------------------------------------------------

if(typelimiter/='ppm' .and. typelimiter /= 'mp5')then
 jxR^L=ixR^L+kr(idims,^D);
 ixCmax^D=jxRmax^D; ixCmin^D=ixLmin^D-kr(idims,^D);
 jxC^L=ixC^L+kr(idims,^D);

 do iw=1,nwflux
   if (loglimit(iw)) then
      w(ixCmin^D:jxCmax^D,iw)=dlog10(w(ixCmin^D:jxCmax^D,iw))
      wLC(ixL^S,iw)=dlog10(wLC(ixL^S,iw))
      wRC(ixR^S,iw)=dlog10(wRC(ixR^S,iw))
   end if

   dwC(ixC^S)=w(jxC^S,iw)-w(ixC^S,iw)

   savetypelimiter=typelimiter
   if(savetypelimiter=='koren') typelimiter='korenL'
   if(savetypelimiter=='cada')  typelimiter='cadaL'
   if(savetypelimiter=='cada3') typelimiter='cada3L'
   call dwlimiter2(dwC,ixI^L,ixC^L,idims,ldw,dxdim)

   wLtmp(ixL^S,iw)=wLC(ixL^S,iw)+half*ldw(ixL^S)
   if(savetypelimiter=='koren')then
     typelimiter='korenR'
     call dwlimiter2(dwC,ixI^L,ixC^L,idims,ldw,dxdim)
   endif
   if(savetypelimiter=='cada')then
     typelimiter='cadaR'
     call dwlimiter2(dwC,ixI^L,ixC^L,idims,ldw,dxdim)
   endif
   if(savetypelimiter=='cada3')then
     typelimiter='cada3R'
     call dwlimiter2(dwC,ixI^L,ixC^L,idims,ldw,dxdim)
   endif
   wRtmp(ixR^S,iw)=wRC(ixR^S,iw)-half*ldw(jxR^S)
   typelimiter=savetypelimiter

   if (loglimit(iw)) then
      w(ixCmin^D:jxCmax^D,iw)=10.0d0**w(ixCmin^D:jxCmax^D,iw)
      wLtmp(ixL^S,iw)=10.0d0**wLtmp(ixL^S,iw)
      wRtmp(ixR^S,iw)=10.0d0**wRtmp(ixR^S,iw)
   end if
 end do

 call checkw(.true.,ixI^L,ixL^L,wLtmp,flagL)
 call checkw(.true.,ixI^L,ixR^L,wRtmp,flagR)

 do iw=1,nwflux
   where (flagL(ixL^S).and.flagR(ixR^S))
      wLC(ixL^S,iw)=wLtmp(ixL^S,iw)
      wRC(ixR^S,iw)=wRtmp(ixR^S,iw)
   end where

   if (loglimit(iw)) then
      where (.not.(flagL(ixL^S).and.flagR(ixR^S)))
         wLC(ixL^S,iw)=10.0d0**wLC(ixL^S,iw)
         wRC(ixR^S,iw)=10.0d0**wRC(ixR^S,iw)
      end where
   end if
 enddo
else if (typelimiter .eq. 'ppm') then
 call PPMlimiter(ixI^L,ixM^LL,idims,w,wCT,wLC,wRC)
else
 call MP5limiter(ixI^L,ixL^L,idims,w,wLC,wRC)
endif

end subroutine upwindLRmf
!=============================================================================
subroutine getfluxmf(w,x,ixI^L,ixO^L,iw,idims,f,transport)

! Calculate non-transport flux f_idim[iw] within ixO^L.

use mod_global_parameters

integer, intent(in)             :: ixI^L, ixO^L, iw, idims
double precision, intent(in)    :: w(ixI^S,nw)
double precision, intent(in)    :: x(ixI^S,1:ndim)
double precision,intent(out)    :: f(ixI^S)
!.. local ..
logical :: transport
integer :: idirmin, idir
!-----------------------------------------------------------------------------
transport=.true.

select case (iw)
   ! f_i[b_k]=v_i*b_k-m_k/rho*b_i
   {case (b^C_)
      if (idims==^C) then
         ! f_i[b_i] should be exactly 0, so we do not use the transport flux
         f(ixO^S)=zero
         transport=.false.
      else
         f(ixO^S)= -w(ixO^S,b0_+idims)*w(ixO^S,v0_+^C)
         if (B0field) then
            f(ixO^S)=f(ixO^S) &
                     +w(ixO^S,v0_+idims)*myB0%w(ixO^S,^C) &
                     -myB0%w(ixO^S,idims)*w(ixO^S,v0_+^C)
         end if
      end if\}
end select

end subroutine getfluxmf
!=============================================================================
subroutine tvdlfmf(qdt,ixI^L,ixO^L,idim^LIM,qtC,wCT,qt,wnew,wold,fC,dx^D,x)

! method=='tvdlf'  --> 2nd order TVD-Lax-Friedrich scheme.
! method=='tvdlf1' --> 1st order TVD-Lax-Friedrich scheme.

use mod_global_parameters

double precision, intent(in)                         :: qdt, qtC, qt, dx^D
integer, intent(in)                                  :: ixI^L, ixO^L, idim^LIM
double precision, dimension(ixI^S,1:ndim), intent(in) ::  x
double precision, dimension(ixI^S,1:ndim)             ::  xi
double precision, dimension(ixI^S,1:nw)               :: wCT, wnew, wold
double precision, dimension(ixI^S,1:nwflux,1:ndim)        :: fC

double precision, dimension(ixI^S,1:nw) :: wLC, wRC
double precision, dimension(ixI^S)      :: fLC, fRC
double precision, dimension(ixI^S)      :: cmaxC
double precision :: dxinv(1:ndim),dxdim(1:ndim)
integer :: idims, iw, ix^L, hxO^L, ixC^L, ixCR^L, jxC^L, kxC^L, kxR^L
logical :: transport, logiB
logical, dimension(ixI^S) :: patchw
!-----------------------------------------------------------------------------

logiB=(BnormLF.and.b0_>0)
! The flux calculation contracts by one in the idim direction it is applied.
! The limiter contracts the same directions by one more, so expand ixO by 2.
ix^L=ixO^L;
do idims= idim^LIM
   ix^L=ix^L^LADD2*kr(idims,^D);
end do
if (ixI^L^LTix^L|.or.|.or.) &
   call mpistop("Error in tvdlfmf: Nonconforming input limits")

^D&dxinv(^D)=-qdt/dx^D;
^D&dxdim(^D)=dx^D;
fC=0.d0
do idims= idim^LIM
   if (B0field) then
      select case (idims)
      {case (^D)
         myB0 => myB0_face^D\}
      end select
   end if

   hxO^L=ixO^L-kr(idims,^D);
   ! ixC is centered index in the idim direction from ixOmin-1/2 to ixOmax+1/2
!   ixCmax^D=ixOmax^D; ixCmin^D=hxOmin^D;

{#IFDEF FCT
! Flux-interpolated constrained transport needs one more layer:
   ixCmax^D=ixOmax^D+1; ixCmin^D=hxOmin^D-1;
}{#IFNDEF FCT
   ixCmax^D=ixOmax^D; ixCmin^D=hxOmin^D;
}


   ! Calculate wRC=uR_{j+1/2} and wLC=uL_j+1/2 
   jxC^L=ixC^L+kr(idims,^D);

   kxCmin^D=ixImin^D; kxCmax^D=ixImax^D-kr(idims,^D);
   kxR^L=kxC^L+kr(idims,^D);
   ixCR^L=ixC^L;
 
   wRC(kxC^S,1:nwflux)=wCT(kxR^S,1:nwflux)
   wLC(kxC^S,1:nwflux)=wCT(kxC^S,1:nwflux)

   ! Get interface positions:
   xi(kxC^S,1:ndim) = x(kxC^S,1:ndim)
   xi(kxC^S,idims) = half* ( x(kxR^S,idims)+x(kxC^S,idims) )

   call upwindLRmf(ixI^L,ixCR^L,ixCR^L,idims,wCT,wCT,wLC,wRC,x,dxdim(idims))

   ! For the high order Lax-Friedrich TVDLF scheme the limiter is based on
   ! the maximum eigenvalue, it is calculated in advance.
   ! determine mean state and store in wLC
   wLC(ixC^S,1:nwflux)= &
         half*(wLC(ixC^S,1:nwflux)+wRC(ixC^S,1:nwflux))
   call getcmaxfff(wLC,xi,ixG^LL,ixC^L,idims,cmaxC)

   ! We regain wLC for further use
   wLC(ixC^S,1:nwflux)=two*wLC(ixC^S,1:nwflux)-wRC(ixC^S,1:nwflux)

   ! Calculate fLC=f(uL_j+1/2) and fRC=f(uR_j+1/2) for each iw
   do iw=b0_+1,b0_+ndir
      call getfluxmf(wLC,xi,ixG^LL,ixC^L,iw,idims,fLC,transport)
      call getfluxmf(wRC,xi,ixG^LL,ixC^L,iw,idims,fRC,transport)
      if (transport) then
         fLC(ixC^S)=fLC(ixC^S)+wLC(ixC^S,v0_+idims)*wLC(ixC^S,iw)
         fRC(ixC^S)=fRC(ixC^S)+wRC(ixC^S,v0_+idims)*wRC(ixC^S,iw)
      end if
      ! To save memory we use fLC to store (F_L+F_R)/2=half*(fLC+fRC)
      fLC(ixC^S)=half*(fLC(ixC^S)+fRC(ixC^S))

      ! Add TVDLF dissipation to the flux
      ! To save memory we use fRC to store -cmax*half*(w_R-w_L)
      if (.not.logiB .and. iw==b0_+idims) then
        fRC(ixC^S)=0.d0
      else
        fRC(ixC^S)=-tvdlfeps*cmaxC(ixC^S)*half*(wRC(ixC^S,iw)-wLC(ixC^S,iw))
      end if
      ! fLC contains physical+dissipative fluxes
      fLC(ixC^S)=fLC(ixC^S)+fRC(ixC^S)

      if (slab) then
         fC(ixC^S,iw,idims)=fLC(ixC^S)
      else
         select case (idims)
         {case (^D)
            fC(ixC^S,iw,^D)=mygeo%surfaceC^D(ixC^S)*fLC(ixC^S)\}
         end select
      end if

   end do ! Next iw
end do ! Next idims

{#IFDEF FCT
call fct_average(ixI^L,ixO^L,fC)
}

!Now update the state:
do idims= idim^LIM
   hxO^L=ixO^L-kr(idims,^D);
   do iw=b0_+1,b0_+ndir

      ! Multiply the fluxes by -dt/dx since Flux fixing expects this
      if (slab) then
         fC(ixI^S,iw,idims)=dxinv(idims)*fC(ixI^S,iw,idims)
         wnew(ixO^S,iw)=wnew(ixO^S,iw) &
              + (fC(ixO^S,iw,idims)-fC(hxO^S,iw,idims))
      else
         select case (idims)
         {case (^D)
            fC(ixI^S,iw,^D)=-qdt*fC(ixI^S,iw,idims)
            wnew(ixO^S,iw)=wnew(ixO^S,iw) &
              + (fC(ixO^S,iw,^D)-fC(hxO^S,iw,^D))/mygeo%dvolume(ixO^S)\}
         end select
      end if

   end do ! Next iw

end do ! Next idims

if (.not.slab) call addgeometrymf(qdt,ixI^L,ixO^L,wCT,wnew,x)

end subroutine tvdlfmf
!=============================================================================
subroutine hancockmf(qdt,ixI^L,ixO^L,idim^LIM,qtC,wCT,qt,wnew,dx^D,x)

! The non-conservative Hancock predictor for TVDLFmf

! on entry:
! input available on ixI^L=ixG^L asks for output on ixO^L=ixG^L^LSUBdixB

! one entry: (predictor): wCT -- w_n        wnew -- w_n   qdt=dt/2

! on exit :  (predictor): wCT -- w_n        wnew -- w_n+1/2


! FCT not implemented here

use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L, idim^LIM
double precision, intent(in) :: qdt, qtC, qt, dx^D, x(ixI^S,1:ndim)
double precision, intent(inout) :: wCT(ixI^S,1:nw), wnew(ixI^S,1:nw)

double precision, dimension(ixI^S,1:nw) :: wLC, wRC
double precision, dimension(ixI^S) :: fLC, fRC
double precision :: dxinv(1:ndim),dxdim(1:ndim)
integer :: idims, iw, ix^L, hxO^L, ixtest^L
logical :: transport
logical, dimension(ixI^S) :: patchw
!-----------------------------------------------------------------------------

! Expand limits in each idims direction in which fluxes are added
ix^L=ixO^L;
do idims= idim^LIM
   ix^L=ix^L^LADDkr(idims,^D);
end do
if (ixI^L^LTix^L|.or.|.or.) &
   call mpistop("Error in Hancockmf: Nonconforming input limits")

^D&dxinv(^D)=-qdt/dx^D;
^D&dxdim(^D)=dx^D;
do idims= idim^LIM
   if (B0field) then
      select case (idims)
      {case (^D)
         myB0 => myB0_face^D\}
      end select
   end if

   ! Calculate w_j+g_j/2 and w_j-g_j/2
   ! First copy all variables, then upwind wLC and wRC.
   ! wLC is to the left of ixO, wRC is to the right of wCT.
   hxO^L=ixO^L-kr(idims,^D);

   wRC(hxO^S,1:nwflux)=wCT(ixO^S,1:nwflux)
   wLC(ixO^S,1:nwflux)=wCT(ixO^S,1:nwflux)

   call upwindLRmf(ixI^L,ixO^L,hxO^L,idims,wCT,wCT,wLC,wRC,x,dxdim(idims))

   ! Advect w(iw)
   do iw=b0_+1,b0_+ndir
      ! Calculate the fLC and fRC fluxes
      call getfluxmf(wRC,x,ixI^L,hxO^L,iw,idims,fRC,transport)
      call getfluxmf(wLC,x,ixI^L,ixO^L,iw,idims,fLC,transport)
      if (transport) then
         fRC(hxO^S)=fRC(hxO^S)+wRC(hxO^S,v0_+idims)*wRC(hxO^S,iw)
         fLC(ixO^S)=fLC(ixO^S)+wLC(ixO^S,v0_+idims)*wLC(ixO^S,iw)
      end if

      ! Advect w(iw)
      if (slab) then
         wnew(ixO^S,iw)=wnew(ixO^S,iw)+dxinv(idims)* &
                          (fLC(ixO^S)-fRC(hxO^S))
      else
         select case (idims)
         {case (^D)
            wnew(ixO^S,iw)=wnew(ixO^S,iw)-qdt/mygeo%dvolume(ixO^S) &
                  *(mygeo%surfaceC^D(ixO^S)*fLC(ixO^S) &
                   -mygeo%surfaceC^D(hxO^S)*fRC(hxO^S))\}
         end select
      end if
   end do
end do ! next idims

if (.not.slab) call addgeometrymf(qdt,ixI^L,ixO^L,wCT,wnew,x)

end subroutine hancockmf
!=============================================================================
subroutine fdmf(qdt,ixI^L,ixO^L,idim^LIM,qtC,wCT,qt,wnew,wold,fC,dx^D,x)

use mod_global_parameters

double precision, intent(in)                                     :: qdt, qtC, qt, dx^D
integer, intent(in)                                              :: ixI^L, ixO^L, idim^LIM
double precision, dimension(ixI^S,1:ndim), intent(in)            :: x

double precision, dimension(ixI^S,1:nw), intent(inout)           :: wCT, wnew, wold
double precision, dimension(ixI^S,1:nwflux,1:ndim), intent(out)  :: fC

double precision, dimension(ixI^S)                               :: fCT
double precision, dimension(ixI^S,1:nw)                          :: fm, fp, fmR, fpL
double precision, dimension(ixI^S)                               :: v
double precision                                                 :: dxinv(1:ndim), dxdims
logical                                                          :: transport
integer                                                          :: idims, iw, ixC^L, ix^L, hxO^L, ixCR^L
!-----------------------------------------------------------------------------


^D&dxinv(^D)=-qdt/dx^D;
do idims= idim^LIM

   select case (idims)
      {case (^D) 
      dxdims = dx^D\}
   end select
   if (B0field) then
      myB0 => myB0_cell
   end if

   ! Get fluxes for the whole grid (mesh+dixB)
   {^D& ixCmin^D = ixOmin^D - dixB * kr(idims,^D)\}
   {^D& ixCmax^D = ixOmax^D + dixB * kr(idims,^D)\}

   hxO^L=ixO^L-kr(idims,^D);
   ! ix is centered index in the idim direction from ixOmin-1/2 to ixOmax+1/2
   ixmax^D=ixOmax^D; ixmin^D=hxOmin^D;

   ixCR^L=ixC^L;

   do iw=b0_+1,b0_+ndir
      call getfluxmf(wCT,x,ixG^LL,ixCR^L,iw,idims,fCT,transport)
      if (transport) fCT(ixCR^S) = fCT(ixCR^S) + wCT(ixCR^S,v0_+idims) * wCT(ixCR^S,iw)
      ! Lax-Friedrich splitting:
      fp(ixCR^S,iw) = half * (fCT(ixCR^S) + tvdlfeps * cmax_global * wCT(ixCR^S,iw))
      fm(ixCR^S,iw) = half * (fCT(ixCR^S) - tvdlfeps * cmax_global * wCT(ixCR^S,iw))
   end do ! iw loop
  
   ! now do the reconstruction of fp and fm:
   call reconstructL(ixI^L,ix^L,idims,fp,fpL,dxdims)
   call reconstructR(ixI^L,ix^L,idims,fm,fmR,dxdims)

   do iw=b0_+1,b0_+ndir
      if (slab) then
         fC(ix^S,iw,idims) = dxinv(idims) * (fpL(ix^S,iw) + fmR(ix^S,iw))
         wnew(ixO^S,iw)=wnew(ixO^S,iw)+ &
              (fC(ixO^S,iw,idims)-fC(hxO^S,iw,idims))
      else
         select case (idims)
         {case (^D)
            fC(ix^S,iw,^D)=-qdt*mygeo%surfaceC^D(ix^S) * (fpL(ix^S,iw) + fmR(ix^S,iw))
            wnew(ixO^S,iw)=wnew(ixO^S,iw)+ &
              (fC(ixO^S,iw,^D)-fC(hxO^S,iw,^D))/mygeo%dvolume(ixO^S)\}
         end select
      end if
   end do ! iw loop

end do !idims loop

if (.not.slab) call addgeometry(qdt,ixI^L,ixO^L,wCT,wnew,x)

end subroutine fdmf
!=============================================================================
subroutine centdiff4mf(qdt,ixI^L,ixO^L,idim^LIM,qtC,wCT,qt,w,wold,fC,dx^D,x)

! Advance the flow variables from t to t+qdt within ixO^L by
! fourth order centered differencing in space 
! for the dw/dt+dF_i(w)/dx_i=S type equation.
! wCT contains the time centered variables at time qtC for flux and source.
! w is the old value at qt on input and the new value at qt+qdt on output.

use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L, idim^LIM
double precision, intent(in) :: qdt, qtC, qt, dx^D
double precision :: wCT(ixI^S,1:nw), w(ixI^S,1:nw), wold(ixI^S,1:nw)
double precision, intent(in) :: x(ixI^S,1:ndim)
double precision, dimension(ixI^S,1:ndim) ::  xi
double precision :: fC(ixI^S,1:nwflux,1:ndim)

double precision :: v(ixI^S,ndim), f(ixI^S)
double precision, dimension(ixI^S,1:nw) :: wLC, wRC
double precision, dimension(ixI^S)      :: vLC, vRC,cmaxLC,cmaxRC
double precision :: dxinv(1:ndim), dxdim(1:ndim)
integer :: idims, iw, idirmin,ix^D
integer :: ix^L, hxO^L, ixC^L, jxC^L, hxC^L, kxC^L, kkxC^L, kkxR^L
logical :: transport,patchw(ixI^S)
!-----------------------------------------------------------------------------
! two extra layers are needed in each direction for which fluxes are added.
ix^L=ixO^L;
do idims= idim^LIM
   ix^L=ix^L^LADD2*kr(idims,^D);
end do

if (ixI^L^LTix^L|.or.|.or.) then
   call mpistop("Error in evolve_CentDiff4: Non-conforming input limits")
end if
^D&dxinv(^D)=-qdt/dx^D;
^D&dxdim(^D)=dx^D;

! Add fluxes to w
do idims= idim^LIM
   if (B0field) then
      select case (idims)
      {case (^D)
         myB0 => myB0_face^D\}
      end select
   end if

   ix^L=ixO^L^LADD2*kr(idims,^D);
   hxO^L=ixO^L-kr(idims,^D);

   ixCmin^D=hxOmin^D; ixCmax^D=ixOmax^D;
   hxC^L=ixC^L-kr(idims,^D);
   jxC^L=ixC^L+kr(idims,^D);
   kxC^L=ixC^L+2*kr(idims,^D);

   kkxCmin^D=ixImin^D; kkxCmax^D=ixImax^D-kr(idims,^D);
   kkxR^L=kkxC^L+kr(idims,^D);
   wRC(kkxC^S,1:nwflux)=wCT(kkxR^S,1:nwflux)
   wLC(kkxC^S,1:nwflux)=wCT(kkxC^S,1:nwflux)

   ! Get interface positions:
   xi(kkxC^S,1:ndim) = x(kkxC^S,1:ndim)
   xi(kkxC^S,idims) = half* ( x(kkxR^S,idims)+x(kkxC^S,idims) )

   call upwindLRmf(ixI^L,ixC^L,ixC^L,idims,wCT,wCT,wLC,wRC,x,dxdim(idims))

   ! Calculate velocities from upwinded values
   call getcmaxfff(wLC,xi,ixG^LL,ixC^L,idims,cmaxLC)
   call getcmaxfff(wRC,xi,ixG^LL,ixC^L,idims,cmaxRC)
   ! now take the maximum of left and right states
   vLC(ixC^S)=max(cmaxRC(ixC^S),cmaxLC(ixC^S))

   do iw=b0_+1,b0_+ndir
      ! Get non-transported flux
      call getfluxmf(wCT,xi,ixI^L,ix^L,iw,idims,f,transport)
      ! Add transport flux
      if (transport) f(ix^S)=f(ix^S)+wCT(ix^S,v0_+idims)*wCT(ix^S,iw)
      ! Center flux to interface
      ! f_i+1/2= (-f_(i+2) +7 f_(i+1) + 7 f_i - f_(i-1))/12
      fC(ixC^S,iw,idims)=(-f(kxC^S)+7.0d0*(f(jxC^S)+f(ixC^S))-f(hxC^S))/12.0d0
      ! add rempel dissipative flux, only second order version for now
      ! one could gradually reduce the dissipative flux to improve solutions 
      ! for computing steady states (Keppens et al. 2012, JCP)
      fC(ixC^S,iw,idims)=fC(ixC^S,iw,idims)-tvdlfeps*half*vLC(ixC^S) &
                                     *(wRC(ixC^S,iw)-wLC(ixC^S,iw))

      if (slab) then
         fC(ixC^S,iw,idims)=dxinv(idims)*fC(ixC^S,iw,idims)
         ! result: f_(i+1/2)-f_(i-1/2) = [-f_(i+2)+8(f_(i+1)+f_(i-1))-f_(i-2)]/12
         w(ixO^S,iw)=w(ixO^S,iw)+(fC(ixO^S,iw,idims)-fC(hxO^S,iw,idims))
      else
         select case (idims)
         {case (^D)
            fC(ixC^S,iw,^D)=-qdt*mygeo%surfaceC^D(ixC^S)*fC(ixC^S,iw,^D)
            w(ixO^S,iw)=w(ixO^S,iw)+ &
                 (fC(ixO^S,iw,^D)-fC(hxO^S,iw,^D))/mygeo%dvolume(ixO^S)\}
         end select
      end if
   end do    !next iw
end do       !next idims
if (.not.slab) call addgeometrymf(qdt,ixI^L,ixO^L,wCT,w,x)

end subroutine centdiff4mf
!!=============================================================================
subroutine getdtfff_courant(w,x,ixI^L,ixO^L,dtnew)

! compute CFL limited dt (for variable time stepping)

use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L
double precision, intent(in) :: x(ixI^S,1:ndim)
double precision, intent(inout) :: w(ixI^S,1:nw), dtnew

double precision :: courantmax, dxinv(1:ndim)
double precision :: cmax(ixI^S),tmp(ixI^S),alfven(ixI^S)
integer :: idims
!-----------------------------------------------------------------------------
dtnew=bigdouble
courantmax=zero
^D&dxinv(^D)=one/dxlevel(^D);

do idims=1,ndim
   call getcmaxfff(w,x,ixI^L,ixO^L,idims,cmax)
   cmax_mype = max(cmax_mype,maxval(cmax(ixO^S)))
   if (.not.slab) then
      tmp(ixO^S)=cmax(ixO^S)/mygeo%dx(ixO^S,idims)
      courantmax=max(courantmax,maxval(tmp(ixO^S)))
   else
      tmp(ixO^S)=cmax(ixO^S)*dxinv(idims)
      courantmax=max(courantmax,maxval(tmp(ixO^S)))
   end if
end do
! courantmax='max( c/dx)'
if (courantmax>smalldouble)  dtnew=min(dtnew,cmf_c/courantmax)

end subroutine getdtfff_courant
!=============================================================================
subroutine getcmaxfff(w,x,ixI^L,ixO^L,idims,cmax)
use mod_global_parameters

logical :: new_cmax,needcmin
integer, intent(in) :: ixI^L, ixO^L, idims
double precision, intent(in)    :: x(ixI^S,1:ndim),w(ixI^S,1:nw)
double precision, intent(out) :: cmax(ixI^S)
!-----------------------------------------------------------------------------

! calculate alfven speed
if(B0field) then
  cmax(ixO^S)=dsqrt((^C&(w(ixO^S,b^C_)+myB0%w(ixO^S,^C))**2+ )/w(ixO^S,rho_))
else
  cmax(ixO^S)=dsqrt((^C&w(ixO^S,b^C_)**2+ )/w(ixO^S,rho_))
endif
cmax(ixO^S)=cmax(ixO^S)+dabs(w(ixO^S,v0_+idims))

end subroutine getcmaxfff
!=============================================================================
subroutine divbclean(ixI^L,ixO^L,w,x,qdt)

! Add Janhunen's and Linde's divB related sources to w
use mod_global_parameters

integer, intent(in)             :: ixI^L, ixO^L
double precision, intent(in)    :: x(ixI^S,1:ndim),qdt
double precision, intent(inout) :: w(ixI^S,1:nw)
integer :: idims, ix^L, ixp^L, i^D, iside
double precision :: divb(ixI^S),graddivb(ixI^S),bdivb(ixI^S,1:ndir)
!-----------------------------------------------------------------------------

! Calculate div B
ix^L=ixO^L^LADD1;
call getdivb(w,ixI^L,ix^L,divb)
! for AMR stability, retreat one cell layer from the boarders of level jump
ixp^L=ixO^L;
!do idims=1,ndim
!  select case(idims)
!   {case(^D)
!      do iside=1,2
!        i^DD=kr(^DD,^D)*(2*iside-3);
!        if(leveljump(i^DD)) then
!          if(iside==1) then
!            ixpmin^D=ixOmin^D-i^D
!          else
!            ixpmax^D=ixOmax^D-i^D
!          end if
!        end if
!      end do
!   \}
!  end select
!end do

! Add Linde's diffusive terms
do idims=1,ndim
   ! Calculate grad_idim(divb)
   call gradient(divb,ixI^L,ixp^L,idims,graddivb)

   ! Multiply by Linde's eta*dt = divbdiff*(c_max*dx)*dt = divbdiff*dx**2
   if (slab) then
      graddivb(ixp^S)=graddivb(ixp^S)*cmf_divb/(^D&1.0d0/dxlevel(^D)**2+)
   else
      graddivb(ixp^S)=graddivb(ixp^S)*cmf_divb &
                      /(^D&1.0d0/mygeo%dx(ixp^S,^D)**2+)
   end if
   ! B_idim += eta*grad_idim(divb)
   w(ixp^S,b0_+idims)=w(ixp^S,b0_+idims)+&
         graddivb(ixp^S)-qdt*w(ixp^S,v0_+idims)*divb(ixp^S)
end do

end subroutine divbclean
!============================================================================= 
subroutine addgeometrymf(qdt,ixI^L,ixO^L,wCT,w,x)

! Add geometrical source terms to w

use mod_global_parameters

integer, intent(in)             :: ixI^L, ixO^L
double precision, intent(in)    :: qdt, x(ixI^S,1:ndim)
double precision, intent(inout) :: wCT(ixI^S,1:nw), w(ixI^S,1:nw)
!.. local ..
double precision :: tmp(ixI^S)
integer          :: iw
!-----------------------------------------------------------------------------

select case (typeaxial)
case ('slab')
   ! No source terms in slab symmetry
case ('cylindrical')
{^IFPHI
     ! s[Bphi]=(Bphi*vr-Br*vphi)/radius
     tmp(ixO^S)=(wCT(ixO^S,bphi_)*wCT(ixO^S,v1_) &
                -wCT(ixO^S,br_)*wCT(ixO^S,v3_))
     w(ixO^S,bphi_)=w(ixO^S,bphi_)+qdt*tmp(ixO^S)/x(ixO^S,1)
}
case ('spherical')
   do iw=1,nwflux
      select case (iw)
{^NOONEC
      ! s[b2]=(vr*Btheta-vtheta*Br)/r
      !       + cot(theta)*psi/r
      case (b2_)
         tmp(ixO^S)= wCT(ixO^S,v1_)*wCT(ixO^S,b2_) &
                    -wCT(ixO^S,v2_)*wCT(ixO^S,b1_)
         if (B0field) then
            tmp(ixO^S)=tmp(ixO^S)+wCT(ixO^S,v1_)*myB0_cell%w(ixO^S,2) &
                       -wCT(ixO^S,v2_)*myB0_cell%w(ixO^S,1)
         end if
         ! Divide by radius and add to w
         w(ixO^S,iw)=w(ixO^S,iw)+qdt*tmp(ixO^S)/x(ixO^S,1)
}
{^IFTHREEC
      ! s[b3]=(vr*Bphi-vphi*Br)/r
      !       -cot(theta)*(vphi*Btheta-vtheta*Bphi)/r
      case (b3_)
         tmp(ixO^S)=wCT(ixO^S,v1_)*wCT(ixO^S,b3_) &
                 -wCT(ixO^S,v3_)*wCT(ixO^S,b1_){^NOONED &
                -(wCT(ixO^S,v3_)*wCT(ixO^S,b2_) &
                 -wCT(ixO^S,v2_)*wCT(ixO^S,b3_))*dcos(x(ixO^S,2)) &
                               /dsin(x(ixO^S,2)) }
         if (B0field) then
            tmp(ixO^S)=tmp(ixO^S)+wCT(ixO^S,v1_)*myB0_cell%w(ixO^S,3) &
               -wCT(ixO^S,v3_)*myB0_cell%w(ixO^S,1){^NOONED &
               -(wCT(ixO^S,v3_)*myB0_cell%w(ixO^S,2) &
                -wCT(ixO^S,v2_)*myB0_cell%w(ixO^S,3))*dcos(x(ixO^S,2)) &
                               /dsin(x(ixO^S,2)) }
         end if
         ! Divide by radius and add to w
         w(ixO^S,iw)=w(ixO^S,iw)+qdt*tmp(ixO^S)/x(ixO^S,1)
}
      end select
   end do
end select

end subroutine addgeometrymf
end mod_magnetofriction
!=============================================================================
