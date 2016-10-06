!=============================================================================
! amrvacusr.t.gemmhd
!=============================================================================
INCLUDE:amrvacnul/usrflags.t
!=============================================================================
subroutine initglobaldata_usr

use mod_global_parameters
!-----------------------------------------------------------------------------
eqpar(gamma_)=1.66666667d0
select case(iprob)
 case(1)
  eqpar(eta_)=5.0d-3
 case(2)
  eqpar(eta_)=1.0d-3
 case(3)
  eqpar(eta_)=1.0d-4
 case(4)
  eqpar(eta_)=1.0d-4
case default
  call mpistop('iprob not given')
end select
eqpar(sheetl_)=0.5d0
eqpar(rhorat_)=0.2d0
eqpar(T0_)=0.5d0
eqpar(psi0_)=0.1d0
eqpar(llx_)=xprobmax1-xprobmin1
eqpar(lly_)=xprobmax2-xprobmin2

end subroutine initglobaldata_usr
!=============================================================================
subroutine initonegrid_usr(ixG^L,ix^L,w,x)

! initialize one grid

use mod_global_parameters

integer, intent(in) :: ixG^L,ix^L
double precision, intent(in) :: x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)

integer,dimension(:),allocatable:: seed
real::             ranx1d(ixGlo1:ixGhi1)
{^IFTHREED
real::             ranx2d(ixGlo1:ixGhi1,ixGlo3:ixGhi3)
\}
double precision:: ranx(ixG^T),dv,sigma
integer :: ix2,seed_size

logical :: patchw(ixG^T)
logical, save :: first=.true.
!----------------------------------------------------------------------------

if(typephys/='mhd')call mpistop("this is an MHD problem: set typephys!")

{^IFONED   call mpistop("GEM problem is multiD")}

{^IFTWOD
if(mype==0.and.first)then
   write(*,*)'Doing 2D GEM challenge, resistive MHD'
   write(*,*)'iprob=',iprob
   write(*,*)'resistivity equal to=',eqpar(eta_)
   write(*,*)'Hall parameter equal to=',eqpar(etah_)
   if(iprob==4)write(*,*)'random velocity perturbation added'
   first=.false.
endif
\}

{^IFTHREED
if(mype==0.and.first)then
   write(*,*)'Doing 3D GEM challenge, resistive MHD'
   write(*,*)'iprob=',iprob
   write(*,*)'resistivity equal to=',eqpar(eta_)
   write(*,*)'Hall parameter equal to=',eqpar(etah_)
   if(iprob==4)write(*,*)'random velocity perturbation added'
   first=.false.
endif
\}

! no initial velocity for cases iprob=1,2,3
w(ixG^S,v1_)  =zero
w(ixG^S,v2_)  =zero
{^IFTHREED 
w(ixG^S,v3_)=zero
\}

if(iprob==4)then
 dv=0.01d0
 sigma=2.0d0*eqpar(sheetl_)
 if(first)then
      write(*,*)'seeding random number generator, on mype==',mype
      call random_seed(SIZE=seed_size)
      allocate(seed(seed_size))
      call random_seed(GET=seed(1:seed_size))
 endif
 {^IFTWOD
 call random_number(ranx1d(ixGmin1:ixGmax1))
 do ix2=ixGmin2,ixGmax2
   ranx(ixGmin1:ixGmax1,ix2)=ranx1d(ixGmin1:ixGmax1)-0.5d0
 enddo
 \}
 {^IFTHREED
 call random_number(ranx2d(ixGmin1:ixGmax1,ixGmin3:ixGmax3))
 do ix2=ixGmin2,ixGmax2
   ranx(ixGmin1:ixGmax1,ix2,ixGmin3:ixGmax3)=ranx2d(ixGmin1:ixGmax1,ixGmin3:ixGmax3)-0.5d0
 enddo
 \}
 w(ixG^S,v1_)=dv*ranx(ixG^S)*dexp(-(x(ixG^S,2)/sigma)**2)
 w(ixG^S,v2_)=dv*ranx(ixG^S)*dexp(-(x(ixG^S,2)/sigma)**2)
 {^IFTHREED
 w(ixG^S,v3_)=dv*ranx(ixG^S)*dexp(-(x(ixG^S,2)/sigma)**2)
 \}
endif

! set up the 1D equilibrium variation
w(ixG^S,b1_)  =dtanh(x(ixG^S,2)/eqpar(sheetl_))
w(ixG^S,b2_)  =zero
{^IFTHREED 
w(ixG^S,b3_)=zero
\}

! add the 2D island perturbation
w(ixG^S,b1_)= w(ixG^S,b1_)-eqpar(psi0_) &
              *dcos(two*dpi*x(ixG^S,1)/eqpar(llx_))              &
              *dsin(dpi*x(ixG^S,2)/eqpar(lly_))*dpi/eqpar(lly_)
w(ixG^S,b2_)= w(ixG^S,b2_)+eqpar(psi0_) &
              *dsin(two*dpi*x(ixG^S,1)/eqpar(llx_))              &
              *dcos(dpi*x(ixG^S,2)/eqpar(lly_))*two*dpi/eqpar(llx_)

w(ixG^S,rho_) =eqpar(rhorat_)+one/(dcosh(x(ixG^S,2)/eqpar(sheetl_))**2) 
w(ixG^S,p_)   =eqpar(T0_) &
              *(eqpar(rhorat_)+one/(dcosh(x(ixG^S,2)/eqpar(sheetl_))**2))

patchw(ixG^S)=.false.
call conserve(ixG^L,ixG^L,w,x,patchw)

end subroutine initonegrid_usr
!=============================================================================
subroutine specialbound_usr(qt,ixG^L,ixO^L,iw,iB,w,x)

! special boundary types, user defined

use mod_global_parameters

integer, intent(in) :: ixO^L, iw, iB, ixG^L
double precision, intent(in) :: qt, x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)

logical :: patchw(ixG^T)
integer :: ix2
!----------------------------------------------------------------------------

end subroutine specialbound_usr
!=============================================================================
subroutine specialsource(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L, iw^LIM
double precision, intent(in) :: qdt, qtC, qt, wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
double precision, intent(inout) :: w(ixI^S,1:nw)

! coordinate info can be used in region ixO

! integer :: iw
! double precision :: s(ixG^T)
!-----------------------------------------------------------------------------

end subroutine specialsource
!=============================================================================
subroutine getdt_special(w,ixG^L,ix^L,dtnew,dx^D,x)

use mod_global_parameters

integer, intent(in) :: ixG^L,ix^L
double precision, intent(in) :: dx^D, x(ixG^S,1:ndim), w(ixG^S,1:nw)
double precision, intent(inout) :: dtnew
!-----------------------------------------------------------------------------

dtnew=bigdouble

end subroutine getdt_special
!=============================================================================
subroutine specialeta(w,ixI^L,ix^L,idirmin,x,current,eta)

! Set the common "eta" array for resistive MHD based on w or the
! "current" variable which has components between idirmin and 3.

use mod_global_parameters

integer, intent(in) :: ixI^L, ix^L, idirmin
double precision, intent(in) :: w(ixI^S,nw), x(ixI^S,1:ndim)

double precision :: current(ixG^T,7-2*ndir:3), eta(ixG^T)
!-----------------------------------------------------------------------------

end subroutine specialeta
!=============================================================================
subroutine specialrefine_grid(igrid,level,ixG^L,ix^L,qt,w,x,refine,coarsen)

! Enforce additional refinement or coarsening
! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.

use mod_global_parameters

integer, intent(in) :: igrid, level, ix^L, ixG^L
double precision, intent(in) :: qt, w(ixG^S,1:nw), x(ixG^S,1:ndim)
integer, intent(inout) :: refine, coarsen
!-----------------------------------------------------------------------------

end subroutine specialrefine_grid
!=============================================================================
subroutine specialvarforerrest(ixI^L,ixO^L,iflag,w,var)

! this is the place to compute a local auxiliary variable to be used
! as refinement criterion for the Lohner error estimator only
!  -->it is then requiring and iflag>nw
! note that ixO=ixI=ixG, hence the term local (gradients need special attention!)

use mod_global_parameters

integer, intent(in)          :: ixI^L,ixO^L,iflag
double precision, intent(in) :: w(ixI^S,1:nw)
double precision, intent(out):: var(ixG^T)
!-----------------------------------------------------------------------------

if (iflag >nw)call mpistop(' iflag> nw, make change in parfile or in user file')

var(ixI^S) = zero

end subroutine specialvarforerrest

!=============================================================================
subroutine printlog_special

! printlog: calculates volume averaged mean values 

use mod_global_parameters

logical :: fileopen
integer :: iigrid, igrid, level, nleafs_level(1:nlevelshi), iw, i
double precision :: wmean(1:nw), volume(1:nlevelshi), volprob, voltotal
double precision :: dvolume(ixG^T), volumeflat(1:nlevelshi)
integer :: numlevels
integer, dimension(1:nlevelshi) :: isum_send, isum_recv
double precision, dimension(1:nw+1+nlevelshi) :: dsum_send, dsum_recv
character(len=80) :: filename
character(len=1024) :: line
logical, save :: opened=.false.
integer :: amode, status(MPI_STATUS_SIZE)
!-----------------------------------------------------------------------------
volume(1:mxnest)=zero
volumeflat(1:mxnest)=zero
wmean(1:nw)= zero
nleafs_level(1:mxnest)=0

do iigrid=1,igridstail; igrid=igrids(iigrid);
   level=node(plevel_,igrid)
   nleafs_level(level)=nleafs_level(level)+1
   volumeflat(level)=volumeflat(level)+ &
          {(rnode(rpxmax^D_,igrid)-rnode(rpxmin^D_,igrid))|*}
   if (slab) then
      dvolume(ixM^T)={rnode(rpdx^D_,igrid)|*}
   else
      dvolume(ixM^T)=pgeo(igrid)%dvolume(ixM^T)
      volume(level)=volume(level)+sum(dvolume(ixM^T))
   end if
   wmean(rho_)=wmean(rho_)+sum(dvolume(ixM^T)*pw(igrid)%w(ixM^T,rho_))
   wmean(m1_)=wmean(m1_)+half*sum(dvolume(ixM^T)*(pw(igrid)%w(ixM^T,m1_)**2)/pw(igrid)%w(ixM^T,rho_))
   wmean(m2_)=wmean(m2_)+half*sum(dvolume(ixM^T)*(pw(igrid)%w(ixM^T,m2_)**2)/pw(igrid)%w(ixM^T,rho_))
 {^IFTHREED
   wmean(m3_)=wmean(m3_)+half*sum(dvolume(ixM^T)*(pw(igrid)%w(ixM^T,m3_)**2)/pw(igrid)%w(ixM^T,rho_))
 }
   wmean(e_)=wmean(e_)+sum(dvolume(ixM^T)*pw(igrid)%w(ixM^T,e_))
   wmean(b1_)=wmean(b1_)+half*sum(dvolume(ixM^T)* &
       ((pw(igrid)%w(ixM^T,b1_)-dtanh(px(igrid)%x(ixM^T,2)/eqpar(sheetl_)))**2+pw(igrid)%w(ixM^T,b2_)**2))
   wmean(b2_)=wmean(b2_)+half*sum(dvolume(ixM^T)*(pw(igrid)%w(ixM^T,b2_)**2))
 {^IFTHREED
   wmean(b3_)=wmean(b3_)+half*sum(dvolume(ixM^T)*(pw(igrid)%w(ixM^T,b3_)**2))
 }
end do
if (slab) volume(levmin:levmax)=volumeflat(levmin:levmax)

voltotal=sum(volume(levmin:levmax))

numlevels=levmax-levmin+1
dsum_send(1:nw)=wmean(1:nw)
dsum_send(nw+1)=voltotal
dsum_send(nw+2:nw+1+numlevels)=volumeflat(levmin:levmax)
call MPI_REDUCE(dsum_send,dsum_recv,nw+1+numlevels,MPI_DOUBLE_PRECISION, &
                MPI_SUM,0,icomm,ierrmpi)
isum_send(1:numlevels)=nleafs_level(levmin:levmax)
call MPI_REDUCE(isum_send,isum_recv,numlevels,MPI_INTEGER, &
                MPI_SUM,0,icomm,ierrmpi)

if (mype==0) then

   wmean(1:nw)=dsum_recv(1:nw)
   voltotal=dsum_recv(nw+1)
   volumeflat(levmin:levmax)=dsum_recv(nw+2:nw+1+numlevels)
   nleafs_level(levmin:levmax)=isum_recv(1:numlevels)

   wmean=wmean/voltotal

   ! determine coverage in coordinate space
   volprob={(xprobmax^D-xprobmin^D)|*}
   volumeflat(levmin:levmax)=volumeflat(levmin:levmax)/volprob

   if (.not.opened) then
      ! generate filename
      write(filename,"(a,a)") TRIM(filenamelog),".log"

      amode=ior(MPI_MODE_CREATE,MPI_MODE_WRONLY)
      call MPI_FILE_OPEN(MPI_COMM_SELF,filename,amode, &
                         MPI_INFO_NULL,log_fh,ierrmpi)
      opened=.true.

      call MPI_FILE_WRITE(log_fh,fileheadout,len_trim(fileheadout), &
                          MPI_CHARACTER,status,ierrmpi)
      !!call MPI_FILE_WRITE(log_fh,new_line('a'),1,MPI_CHARACTER,status,ierrmpi)
      call MPI_FILE_WRITE(log_fh,achar(10),1,MPI_CHARACTER,status,ierrmpi)

      i=len_trim(wnames)-1
      do level=1,mxnest
          i=i+3
          write(wnames(i:i+1),"(a,i1)") "c",level
      end do
      do level=1,mxnest
          i=i+3
          write(wnames(i:i+1),"(a,i1)") "n",level
      end do
      if (time_accurate) then
         if(residmin>smalldouble) then
           write(line,'(a15,a79)')"it   t  dt res ",wnames
         else
           write(line,'(a15,a79)')"it   t   dt    ",wnames
         endif
      else
         if(residmin>smalldouble) then
           write(line,'(a7,a79)')"it res ",wnames
         else
           write(line,'(a7,a79)')"it     ",wnames
         endif
      end if
      call MPI_FILE_WRITE(log_fh,line,len_trim(line),MPI_CHARACTER, &
                          status,ierrmpi)
   end if
   !!call MPI_FILE_WRITE(log_fh,new_line('a'),1,MPI_CHARACTER,status,ierrmpi)
   call MPI_FILE_WRITE(log_fh,achar(10),1,MPI_CHARACTER,status,ierrmpi)

   if (time_accurate) then
      if(residmin>smalldouble) then
         write(line,'(i7,3(e13.5))')it,t,dt,residual
      else
         write(line,'(i7,2(e13.5))')it,t,dt
      endif
   else
      if(residmin>smalldouble) then
         write(line,'(i7,1(e13.5))')it,residual
      else
         write(line,'(i7)')it
      endif
   end if
   call MPI_FILE_WRITE(log_fh,line,len_trim(line), &
                       MPI_CHARACTER,status,ierrmpi)
   do iw=1,nw
      write(line,'(e13.5)')wmean(iw)
      call MPI_FILE_WRITE(log_fh,line,len_trim(line), &
                          MPI_CHARACTER,status,ierrmpi)
   end do
   do level=1,mxnest
      write(line,'(e13.5)')volumeflat(level)
      call MPI_FILE_WRITE(log_fh,line,len_trim(line), &
                          MPI_CHARACTER,status,ierrmpi)
   end do
   do level=1,mxnest
      write(line,'(i6)') nleafs_level(level)
      call MPI_FILE_WRITE(log_fh,line,len_trim(line), &
                          MPI_CHARACTER,status,ierrmpi)
   end do
end if

end subroutine printlog_special
!=============================================================================
subroutine process_grid_usr(igrid,level,ixI^L,ixO^L,qt,w,x)

! this subroutine is ONLY to be used for computing auxiliary variables
! which happen to be non-local (like div v), and are in no way used for
! flux computations. As auxiliaries, they are also not advanced

use mod_global_parameters

integer, intent(in):: igrid,level,ixI^L,ixO^L
double precision, intent(in):: qt,x(ixI^S,1:ndim)
double precision, intent(inout):: w(ixI^S,1:nw)
!-----------------------------------------------------------------------------

end subroutine process_grid_usr
!=============================================================================
subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)

! this subroutine can be used in convert, to add auxiliary variables to the
! converted output file, for further analysis using tecplot, paraview, ....
! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
!
! the array normconv can be filled in the (nw+1:nw+nwauxio) range with
! corresponding normalization values (default value 1)

use mod_global_parameters

integer, intent(in)                :: ixI^L,ixO^L
double precision, intent(in)       :: x(ixI^S,1:ndim)
double precision                   :: w(ixI^S,nw+nwauxio)
double precision                   :: normconv(0:nw+nwauxio)


double precision :: current(ixG^T,7-2*ndir:3)

double precision                   :: gradrho(ixG^T),rho(ixG^T),drho(ixG^T)
double precision                   :: kk,kk0,grhomax,kk1
integer                            :: idims


integer, parameter:: idirmin0=7-2*ndir
integer :: idirmin,idir
double precision :: curlv(ixG^T,7-2*ndir:3),vvec(ixG^T,1:ndir)
double precision:: divb(ixG^T)
!-----------------------------------------------------------------------------

call getcurrent(w,ixI^L,ixO^L,idirmin,current)
w(ixO^S,nw+1)=current(ixO^S,3)

rho(ixI^S)=w(ixI^S,rho_)
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
kk0=0.01d0
kk1=1.0d0
grhomax=1000.0d0

w(ixO^S,nw+2)=dexp(-kk*(gradrho(ixO^S)-kk0*grhomax)/(kk1*grhomax-kk0*grhomax))

end subroutine specialvar_output
!=============================================================================
subroutine specialvarnames_output

! newly added variables to be concatenated with the primnames/wnames string

use mod_global_parameters
!-----------------------------------------------------------------------------

primnames= TRIM(primnames)//' '//'jz schlier'
wnames=TRIM(wnames)//' '//'jz schlier'

end subroutine specialvarnames_output
!=============================================================================
subroutine specialset_B0(ixI^L,ixO^L,x,wB0)

! Here one can add a steady (time-independent) potential background field

use mod_global_parameters

integer, intent(in)           :: ixI^L,ixO^L
double precision, intent(in)  :: x(ixG^T,1:ndim)
double precision, intent(inout) :: wB0(ixI^S,1:ndir)
!-----------------------------------------------------------------------------
call mpistop(' abs(Busr)> 0, make change in parfile or in user file')

wB0(ixO^S,1:ndir)=wB0(ixO^S,1:ndir)
!!wB0(ixO^S,1:ndir)=wB0(ixO^S,1:ndir)+user defined steady potential field

end subroutine specialset_B0
!=============================================================================
subroutine specialsource_impl(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L, iw^LIM
double precision, intent(in) :: qdt, qtC, qt, x(ixI^S,1:ndim)
double precision, intent(inout) :: w(ixI^S,1:nw), wCT(ixI^S,1:nw)
!-----------------------------------------------------------------------------

end subroutine specialsource_impl
!=============================================================================
subroutine getdt_impl(w,ixG^L,ix^L,dtnew,dx^D,x)

use mod_global_parameters

integer, intent(in) :: ixG^L,ix^L
double precision, intent(in) :: dx^D, x(ixG^S,1:ndim)
! note that depending on strictsmall etc, w values may change
double precision, intent(inout) :: w(ixG^S,1:nw), dtnew
!-----------------------------------------------------------------------------
dtnew=bigdouble

end subroutine getdt_impl

!=============================================================================
subroutine bc_int(qt,ixG^L,ixO^L,w,x)

! internal boundary, user defined
!
! This subroutine can be used to artificially overwrite ALL conservative 
! variables in a user-selected region of the mesh, and thereby act as
! an internal boundary region. It is called just before external (ghost cell)
! boundary regions will be set by the BC selection. Here, you could e.g. 
! want to introduce an extra variable (nwextra, to be distinguished from nwaux)
! which can be used to identify the internal boundary region location.
! Its effect should always be local as it acts on the mesh.
!

use mod_global_parameters

integer, intent(in) :: ixG^L,ixO^L
double precision, intent(in) :: qt
double precision, intent(inout) :: w(ixG^S,1:nw)
double precision, intent(in) :: x(ixG^S,1:ndim)

! .. local ..
!logical :: patchw(ixG^T)
!----------------------------------------------------------------------------

call mpistop("bc_int not defined")

end subroutine bc_int
!=============================================================================
! amrvacusr.t.gemmhd
!=============================================================================
