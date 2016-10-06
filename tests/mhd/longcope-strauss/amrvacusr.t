!=============================================================================
! amrvacusr.t.islandmhd22

! INCLUDE:amrvacnul/specialini.t
!INCLUDE:amrvacnul/speciallog.t
INCLUDE:amrvacnul/specialbound.t
!INCLUDE:amrvacnul/specialsource.t
!=============================================================================
subroutine initglobaldata_usr

use mod_global_parameters
!-----------------------------------------------------------------------------
eqpar(gamma_)=5.0d0/3.0d0

eqpar(etah_)=0.0000d0

select case(iprob)
 case(1)
   eqpar(eta_)=zero
 case(2)
   eqpar(eta_)=0.001d0
   eqpar(etah_)=0.001d0
 case(3)
   eqpar(eta_)=0.0001d0
 case(4)
   eqpar(eta_)=0.0005d0
 case(5)
   eqpar(eta_)=0.0004d0
 case(6)
   eqpar(eta_)=0.0003d0
 case(7)
   eqpar(eta_)=0.0002d0
 case(8)
   eqpar(eta_)=3.5d-5
 case(9)
   eqpar(eta_)=6.94d-5
endselect 

end subroutine initglobaldata_usr
!=============================================================================
subroutine initonegrid_usr(ixG^L,ix^L,w,x)

! initialize one grid 

use mod_global_parameters

integer, intent(in) :: ixG^L,ix^L
double precision, intent(in) :: x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)

logical::          first
double precision:: a0val,phi0val,p0val,T0val,tpi
double precision:: minbb,maxbb,minbeta,maxbeta,minrat,maxrat
double precision:: tmp(ixG^T)

logical, dimension(ixG^T)           :: patchw(ixG^T)

data first/.true./
!----------------------------------------------------------------------------

if (typephys/='mhd') then
   call mpistop("test problem is MHD problem: set typephys!")
end if

   ! setamrvac -d=22 -phi=0 -z=0 -g=16,16 -p=mhd -u=islandmhd22
   {^IFONED   call mpistop("prob is 2D") }
   {^IFTHREED call mpistop("prob is 2D") }
   {^IFTWOD
   a0val=one/(dpi*dsqrt(two))
   phi0val=0.001d0
   p0val=0.1d0
   T0val=one
   tpi=two*dpi
   w(ix^S,p_)=p0val+two*(dsin(tpi*x(ix^S,1))*dsin(tpi*x(ix^S,2)))**2
   w(ix^S,rho_)=w(ix^S,p_)/T0val
   w(ix^S,b1_)=tpi*a0val*dsin(tpi*x(ix^S,1))*dcos(tpi*x(ix^S,2))
   w(ix^S,b2_)=-tpi*a0val*dcos(tpi*x(ix^S,1))*dsin(tpi*x(ix^S,2))
   w(ix^S,v1_)=tpi*phi0val*dsin(tpi*x(ix^S,2))
   w(ix^S,v2_)=tpi*phi0val*dsin(tpi*x(ix^S,1))
   if(first .and. mype==0)then
      write(*,*)'Doing 2D MHD, Longcope Strauss test'
      write(*,*)'a0val,phi0val,p0val:'
      write(*,*) a0val,phi0val,p0val
   tmp(ix^S)=0.5d0*(w(ix^S,b1_)**2+w(ix^S,b2_)**2)
   minbb=minval(tmp(ix^S))
   maxbb=maxval(tmp(ix^S))
      write(*,*) 'min B2=',minbb,' max B2=',maxbb
   tmp(ix^S)=two*w(ix^S,p_)/( w(ix^S,b1_)**2+ w(ix^S,b2_)**2)
   minbeta=minval(tmp(ix^S))
   maxbeta=maxval(tmp(ix^S))
      write(*,*) 'min beta=',minbeta,' max beta=',maxbeta
   tmp(ix^S)= w(ix^S,rho_)*( w(ix^S,v1_)**2+ w(ix^S,v2_)**2) &
             /( w(ix^S,b1_)**2+ w(ix^S,b2_)**2)
   minrat=minval(tmp(ix^S))
   maxrat=maxval(tmp(ix^S))
      write(*,*) 'min kin/mag=',minrat,' max kin/mag=',maxrat
   endif
   if(first)then
      first=.false.
   endif
   patchw(ix^S)=.false.
   call conserve(ixG^L,ix^L,w,x,patchw)
   \}

end subroutine initonegrid_usr
!=============================================================================
subroutine specialsource(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

! Calculate w(iw)=w(iw)+qdt*SOURCE[wCT,qtC,x] within ixO for all indices
! iw=iwmin...iwmax.  wCT is at time qCT

use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L, iw^LIM
double precision, intent(in) :: qdt, qtC, qt, x(ixI^S,1:ndim)
double precision, intent(inout) :: wCT(ixI^S,1:nw), w(ixI^S,1:nw)

! integer :: iw
! double precision :: s(ixG^T)
!-----------------------------------------------------------------------------

! do iw= iw^LIM
!    select case(iw)
!    case(m1_)
!       ! The source is based on the time centered wCT
!       call getmyforce(wCT,ixO^L,s)
!       w(ixO^S,m1_)=w(ixO^S,m1_) + qdt*s(ixO^S)
!    case(e_)
!       call getmyheating(wCT,ixO^L,s)
!       w(ixO^S,e_) =w(ixO^S,e_)  + qdt*s(ixO^S)
!    end select
! end do

end subroutine specialsource
!=============================================================================
subroutine getdt_special(w,ixG^L,ix^L,dtnew,dx^D,x)

! Limit "dt" further if necessary, e.g. due to the special source terms.
! The getdt_courant (CFL condition) and the getdt subroutine in the AMRVACPHYS
! module have already been called.

use mod_global_parameters

integer, intent(in) :: ixG^L, ix^L
double precision, intent(in) :: dx^D, x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw), dtnew
!-----------------------------------------------------------------------------

dtnew=bigdouble

end subroutine getdt_special
!=============================================================================
subroutine specialeta(w,ixI^L,ix^L,idirmin,x,current,eta)

! Set the "eta" array for resistive MHD based on w or the
! "current" variable which has components between idirmin and 3.

use mod_global_parameters

integer, intent(in) :: ixI^L, ix^L, idirmin
double precision, intent(in) :: w(ixI^S,nw), x(ixI^S,1:ndim)

double precision :: current(ixG^T,7-2*ndir:3), eta(ixG^T)
!-----------------------------------------------------------------------------

!  eta(ix^S)=...

call mpistop("specialeta is not defined")

end subroutine specialeta
!=============================================================================
subroutine specialrefine_grid(igrid,level,ixG^L,ix^L,qt,w,x,refine,coarsen)

! Enforce additional refinement or coarsening
! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.

use mod_global_parameters

integer, intent(in) :: igrid, level, ixG^L, ix^L
double precision, intent(in) :: qt, w(ixG^S,1:nw), x(ixG^S,1:ndim)
integer, intent(inout) :: refine, coarsen
double precision :: xmid
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

var(ixI^S) = {^C& w(ixI^S,v^C_)**2 +}

end subroutine specialvarforerrest
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
subroutine specialset_B0(ixI^L,ixO^L,x,wB0)

! Here one can add a steady (time-independent) potential background field

use mod_global_parameters

integer, intent(in)           :: ixI^L,ixO^L
double precision, intent(in)  :: x(ixG^T,1:ndim)
double precision, intent(inout) :: wB0(ixI^S,1:ndir)
!-----------------------------------------------------------------------------

end subroutine specialset_B0
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


integer:: idirmin
double precision :: current(ixG^T,7-2*ndir:3)
double precision :: maxjz
double precision :: maxjzglobal
double precision :: maxjz_send,maxjz_recv
!-----------------------------------------------------------------------------
volume(1:mxnest)=zero
volumeflat(1:mxnest)=zero
wmean(1:nw)= zero
nleafs_level(1:mxnest)=0
maxjz = 0
maxjzglobal = 0

do iigrid=1,igridstail; igrid=igrids(iigrid);
   ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
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
   do iw=1,nw
      wmean(iw)=wmean(iw)+sum(dvolume(ixM^T)*pw(igrid)%w(ixM^T,iw))
   end do
   call getcurrent(pw(igrid)%w,ixG^LL,ixM^LL,idirmin,current)
   maxjz = maxval(dabs(current(ixM^T,3)))
   maxjzglobal=max(maxjz,maxjzglobal)
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

call MPI_REDUCE(maxjzglobal,maxjz_recv,1,MPI_DOUBLE_PRECISION, &
		MPI_MAX,0,icomm,ierrmpi)



if (mype==0) then

   wmean(1:nw)=dsum_recv(1:nw)
   maxjz = maxjz_recv
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
      amode=ior(amode,MPI_MODE_APPEND)
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
          if(level<10) then
            if (i+1<79) write(wnames(i:i+1),"(a,i1)") "c",level
          else
            if (i+2<79) write(wnames(i:i+2),"(a,i2)") "c",level
          endif
      end do
      do level=1,mxnest
          i=i+3
          if(level<10) then
            if (i+1<79) write(wnames(i:i+1),"(a,i1)") "n",level
          else
            if (i+2<79) write(wnames(i:i+2),"(a,i2)") "n",level
          endif
      end do

      if(residmin>smalldouble) then
        write(line,'(a15,a79)')"it   t  dt res ",wnames
      else
        write(line,'(a15,a79)')"it   t  dt jz  ",wnames
      endif

      call MPI_FILE_WRITE(log_fh,line,len_trim(line),MPI_CHARACTER, &
                          status,ierrmpi)
   end if
   !!call MPI_FILE_WRITE(log_fh,new_line('a'),1,MPI_CHARACTER,status,ierrmpi)
   call MPI_FILE_WRITE(log_fh,achar(10),1,MPI_CHARACTER,status,ierrmpi)

   
   write(line,'(i7,3(e13.5))')it,t,dt,maxjz

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

integer:: idirmin

double precision :: current(ixG^T,7-2*ndir:3)
!-----------------------------------------------------------------------------

call getcurrent(w,ixI^L,ixO^L,idirmin,current)

! Example: assuming nwauxio=1 at convert stage and desire to see -w(1)
w(ixO^S,nw+1)=current(ixO^S,3)

end subroutine specialvar_output
!=============================================================================
subroutine specialvarnames_output

! newly added variables to be concatenated with the primnames/wnames string

use mod_global_parameters
!-----------------------------------------------------------------------------

primnames= TRIM(primnames)//' '//'jz'
wnames=TRIM(wnames)//' '//'jz'

end subroutine specialvarnames_output
!=============================================================================
      subroutine flag_grid_usr(qt,ixG^L,ixO^L,w,x,flag)
      
      use mod_global_parameters

      integer, intent(in)             :: ixG^L, ixO^L
      integer, intent(inout)          :: flag
      double precision, intent(in)    :: qt
      double precision, intent(inout) :: w(ixG^S,1:nw)
      double precision, intent(in)    :: x(ixG^S,1:ndim)
      
      ! flag=0 : Treat as normal domain
      ! flag=1 : Treat as passive, but reduce by safety belt
      ! flag=2 : Treat as passive and don't care about safety belt

!-----------------------------------------------------------------------------
      flag = 0
      
      end subroutine flag_grid_usr
!=============================================================================

!=============================================================================
! amrvacusr.t.islandmhd22
!=============================================================================
