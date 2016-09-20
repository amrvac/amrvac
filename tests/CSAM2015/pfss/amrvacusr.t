!=============================================================================
! amrvacusr.t.pfss
!=============================================================================
INCLUDE:amrvacnul/specialbound.t
!INCLUDE:amrvacnul/specialsource.t
INCLUDE:amrvacnul/specialimpl.t
INCLUDE:amrvacnul/usrflags.t
INCLUDE:amrvacmodules/pfss.t
!=============================================================================
subroutine initglobaldata_usr 
use harm_coef_data
use mod_global_parameters

logical, save :: firstglobalusr=.true.
!-----------------------------------------------------------------------------
! CGS Unit
mHunit=1.67262d-24 ! g
k_B=1.3806d-16     ! erg*K^-1
miu0=4.d0*dpi      ! Gauss^2 cm^2 dyne^-1

Lunit=6.955d10 !cm
nHunit=1.d9 !cm^-3
Teunit=1.d6 !K

runit=1.4d0*mHunit*nHunit     ! 2.341668e-15 g*cm^-3
punit=2.3d0*nHunit*k_B*Teunit ! 0.317538 erg*cm^-3
Bunit=dsqrt(miu0*punit)       ! 1.99757357615242 Gauss
!vunit=Bunit/dsqrt(miu0*runit) ! 11644884.6777562 cm/s 
vunit=dsqrt(punit/runit) ! 11644884.6777562 cm/s 
tunit=Lunit/vunit             ! 5972.57954240224 s about 100 min

eqpar(grav1_)=-2.74d4*Lunit/vunit**2
eqpar(grav2_)=0.d0
eqpar(grav3_)=0.d0

!R_s=5.d0
! base density and temperature
rhob=1.d0
Tiso=1.d6/Teunit

if(firstglobalusr) then
  trunc=.false.
  lmax=720
  call harm_coef('hmisynopticmap2111.dat')
  firstglobalusr=.false.
endif

end subroutine initglobaldata_usr 
!=============================================================================
subroutine initonegrid_usr(ixG^L,ix^L,w,x)

!initialize one grid within ixI^L
use harm_coef_data, only: R_s, R_0
use mod_global_parameters

integer, intent(in) :: ixG^L, ix^L
double precision, intent(in) :: x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)

double precision :: rss,bpf(ixG^T,1:ndir),bpfss(ixG^T,1:ndir),xrss(ixG^S,1:ndim)
double precision :: fmass,BB(ixG^T)
integer :: ix1
logical, save:: first=.true.
!-----------------------------------------------------------------------------
rss=R_s
if(mype==0 .and. first) then
  print*,'Global Sun initializing grids'
  first=.false.
endif

w(ix^S,rho_)=rhob*dexp(eqpar(grav1_)*R_0**2/Tiso*(1.d0/R_0-1.d0/(x(ix^S,1))))

if(B0field) then
  w(ix^S,b1_)=0.d0
  w(ix^S,b2_)=0.d0
  w(ix^S,b3_)=0.d0
else 
  xrss=x
  xrss(ixG^S,r_)=rss
  call pfss(ixG^L,ix^L,bpfss,xrss) 
  if(any(x(ix^S,r_)<rss)) then
    call pfss(ixG^L,ix^L,bpf,x) 
    w(ix^S,b1_:b3_)=bpf(ix^S,1:3)
    where(x(ix^S,r_)>=rss)
      w(ix^S,b1_)=bpfss(ix^S,1)*(rss/x(ix^S,r_))**2
      w(ix^S,b2_)=0.d0
      w(ix^S,b3_)=0.d0
    endwhere
  else
    w(ix^S,b1_)=bpfss(ix^S,1)*(rss/x(ix^S,r_))**2
    w(ix^S,b2_)=0.d0
    w(ix^S,b3_)=0.d0
  endif
endif
w(ix^S,b1_:b3_)=w(ix^S,b1_:b3_)/Bunit

w(ix^S,m1_)=0.d0
w(ix^S,m2_)=0.d0
w(ix^S,m3_)=0.d0

end subroutine initonegrid_usr
!=============================================================================
subroutine getggrav(ggrid,ixI^L,ixO^L,x)

use harm_coef_data, only: R_0
use mod_global_parameters

integer, intent(in) :: ixI^L,ixO^L
double precision, intent(in) :: x(ixI^S,1:ndim)
double precision, intent(out) :: ggrid(ixG^T)

integer :: ix1
!---------------------------------------------------------------------------
! calculate gravity
do ix1=ixOmin1,ixOmax1
  ggrid(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=eqpar(grav1_)*(R_0/&
      x(ix1,ixOmin2,ixOmin3,1))**2
enddo

end subroutine getggrav
!=============================================================================
subroutine addsource_gravSA(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

! w[iw]=w[iw]+qdt*S[wCT,qtC,x] where S is the source based on wCT within ixO
! gravity distribution along a magnetic loop (a circular arc)

use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L, iw^LIM
double precision, intent(in) :: qdt, qtC, qt, x(ixI^S,1:ndim)
double precision, intent(inout) :: wCT(ixI^S,1:nw), w(ixI^S,1:nw)

double precision :: ggrid(ixG^T)
integer :: iw, idims
!---------------------------------------------------------------------------

call getggrav(ggrid,ixI^L,ixO^L,x)

! add sources from gravity
do iw= iw^LIM
   select case (iw)
   case (m^D_)
     ! dm_i/dt= +rho*g_i
      idims=iw-m0_
      if (abs(eqpar(grav0_+idims))>smalldouble) &
          w(ixO^S,m0_+idims)=w(ixO^S,m0_+idims) &
              +qdt*ggrid(ixO^S)*wCT(ixO^S,rho_)
{#IFDEF ENERGY
   case (e_)
     ! de/dt= +g_i*m_i
      do idims=1,ndim
         if (abs(eqpar(grav0_+idims))>smalldouble) &
            w(ixO^S,e_)=w(ixO^S,e_) &
              +qdt*ggrid(ixO^S)*wCT(ixO^S,m0_+idims)
      end do
}
   end select
end do

end subroutine addsource_gravSA
!=============================================================================
subroutine getdt_grav(w,ixG^L,ix^L,dtnew,dx^D,x)

use mod_global_parameters

integer, intent(in) :: ixG^L, ix^L
double precision, intent(in) :: dx^D, x(ixG^S,1:ndim), w(ixG^S,1:nw)
double precision, intent(inout) :: dtnew

double precision:: dxinv(1:ndim), dtgrav
integer:: idims
!----------------------------------------------------------------------------

^D&dxinv(^D)=one/dx^D;
dtgrav=bigdouble
do idims=1,ndim
   if(abs(eqpar(grav0_+idims))>zero)&
   dtgrav=min(dtgrav,one/sqrt(abs(eqpar(grav0_+idims))*dxinv(idims)))
enddo

dtnew=dtgrav

end subroutine getdt_grav
!=============================================================================
subroutine specialsource(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L, iw^LIM
double precision, intent(in) :: qdt, qtC, qt, x(ixI^S,1:ndim)
double precision, intent(inout) :: w(ixI^S,1:nw), wCT(ixI^S,1:nw)
!-----------------------------------------------------------------------------
call addsource_gravSA(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

end subroutine specialsource
!=============================================================================
subroutine getdt_special(w,ixG^L,ix^L,dtnew,dx^D,x)

use mod_global_parameters

integer, intent(in) :: ixG^L,ix^L
double precision, intent(in) :: dx^D, x(ixG^S,1:ndim)
! note that depending on strictsmall etc, w values may change 
double precision, intent(inout) :: w(ixG^S,1:nw), dtnew
!-----------------------------------------------------------------------------
dtnew=bigdouble

call getdt_grav(w,ixG^L,ix^L,dtnew,dx^D,x)

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

! you must set consistent values for integers refine/coarsen:

! refine = -1 enforce to not refine
! refine =  0 doesn't enforce anything
! refine =  1 enforce refinement

! coarsen = -1 enforce to not coarsen
! coarsen =  0 doesn't enforce anything
! coarsen =  1 enforce coarsen

use mod_global_parameters

integer, intent(in) :: igrid, level, ixG^L, ix^L
double precision, intent(in) :: qt, w(ixG^S,1:nw), x(ixG^S,1:ndim)
integer, intent(inout) :: refine, coarsen
!-----------------------------------------------------------------------------

!if(any(x(ix^S,2)<1.1d0 .and. x(ix^S,2)>2.2d0 .and. x(ix^S,3)<3.4d0 .and. &
!   x(ix^S,3)>4.8d0)) then
!  if(level>2) then
!    refine=-1
!    coarsen=1
!  else if(level==2) then
!    refine=-1
!  endif
!endif

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

use mod_global_parameters
!-----------------------------------------------------------------------------

call printlog_default
call spatial_integral_w

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

double precision                   :: csound2(ixG^T)
!-----------------------------------------------------------------------------
if(B0field) then
  w(ixO^S,nw+1)=w(ixO^S,b1_)+myb0_cell%w(ixO^S,1)
  w(ixO^S,nw+2)=w(ixO^S,b2_)+myb0_cell%w(ixO^S,2)
  w(ixO^S,nw+3)=w(ixO^S,b3_)+myb0_cell%w(ixO^S,3)
else
  w(ixO^S,nw+1)=w(ixO^S,b1_)
  w(ixO^S,nw+2)=w(ixO^S,b2_)
  w(ixO^S,nw+3)=w(ixO^S,b3_)
endif

end subroutine specialvar_output
!=============================================================================
subroutine specialvarnames_output

! newly added variables to be concatenated with the primnames/wnames string

use mod_global_parameters
!-----------------------------------------------------------------------------

primnames= TRIM(primnames)//' '//'br bth bph'
wnames=TRIM(wnames)//' '//'br bth bph'

end subroutine specialvarnames_output
!=============================================================================
subroutine specialset_B0(ixI^L,ixO^L,x,wB0)

! Here one can add a steady (time-independent) potential background field

use mod_global_parameters

integer, intent(in)           :: ixI^L,ixO^L
double precision, intent(in)  :: x(ixG^T,1:ndim)
double precision, intent(inout) :: wB0(ixI^S,1:ndir)

double precision :: bpf(ixG^T,1:ndir)
integer :: ix^L
!-----------------------------------------------------------------------------

ix^L=ixI^L^LSUB2;
call pfss(ixI^L,ixO^L,bpf,x)

wB0(ixO^S,1:ndir)=wB0(ixO^S,1:ndir)+bpf(ixO^S,1:ndir)

end subroutine specialset_B0
!=============================================================================
subroutine spatial_integral_w

use mod_global_parameters

double precision :: dvolume(ixG^T), dsurface(ixG^T),timephy,dvone
double precision, allocatable :: integral_ipe(:), integral_w(:)
double precision, external :: integral_grid

integer           :: nregions,ireg,ncellpe,ncell,idims,hxM^LL,nx^D
integer           :: iigrid,igrid,status(MPI_STATUS_SIZE),ni
character(len=100):: filename,region
character(len=1024) :: line, datastr
logical           :: patchwi(ixG^T),alive
!-----------------------------------------------------------------------------

nregions=1
! number of integrals to perform
ni=3
allocate(integral_ipe(ni),integral_w(ni))
integral_ipe=0.d0
integral_w=0.d0
nx^D=ixMhi^D-ixMlo^D+1;
do ireg=1,nregions
 select case(ireg)
 case(1)
   region='fulldomain'
 case(2)
   region='cropped'
 end select
 ncellpe=0 
 do iigrid=1,igridstail; igrid=igrids(iigrid);
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
   if (.not.slab) mygeo => pgeo(igrid)
   if (B0field) then
      myB0_cell => pB0_cell(igrid)
     {^D&myB0_face^D => pB0_face^D(igrid)\}
   end if
   typelimiter=typelimiter1(node(plevel_,igrid))
   typegradlimiter=typegradlimiter1(node(plevel_,igrid))
   patchwi(ixG^T)=.false.
   select case(region)
   case('fulldomain')
      patchwi(ixM^T)=.true.
      ncellpe=ncellpe+{nx^D*}
   case('cropped')
      call mask_grid(ixG^LL,ixM^LL,pw(igrid)%w,px(igrid)%x,patchwi,ncellpe)
   case default
      call mpistop("region not defined")
   end select
   integral_ipe(1)=integral_ipe(1)+ &
             integral_grid(ixG^LL,ixM^LL,pw(igrid)%w,px(igrid)%x,dvolume,dsurface,1,patchwi)
   integral_ipe(2)=integral_ipe(2)+ &
             integral_grid(ixG^LL,ixM^LL,pw(igrid)%w,px(igrid)%x,dvolume,dsurface,2,patchwi)
   integral_ipe(3)=integral_ipe(3)+ &
             integral_grid(ixG^LL,ixM^LL,pw(igrid)%w,px(igrid)%x,dvolume,dsurface,3,patchwi)
 end do
 call MPI_ALLREDUCE(integral_ipe,integral_w,ni,MPI_DOUBLE_PRECISION,&
                      MPI_SUM,icomm,ierrmpi)
 call MPI_ALLREDUCE(ncellpe,ncell,1,MPI_INTEGER,MPI_SUM,icomm,ierrmpi)
 integral_w(3)=integral_w(3)/dble(ncell)
 timephy=t
 if(mype==0) then
   write(filename,"(a,a,a)") TRIM(filenamelog),TRIM(region),".csv"
   inquire(file=filename,exist=alive)
   if(alive) then
     open(unit=21,file=filename,form='formatted',status='old',access='append')
   else
     open(unit=21,file=filename,form='formatted',status='new')
     write(21,'(a)') 'time, f_i, CWsin, theta'
   endif
   write(datastr,'(es11.4, a)') timephy,','
   line=datastr
   write(datastr,"(es11.4, a)") integral_w(3),','
   line = trim(line)//trim(datastr)
   write(datastr,"(es11.4, a)") integral_w(1)/(integral_w(2)+smalldouble),','
   line = trim(line)//trim(datastr)
   write(datastr,"(es11.4)") dasin(integral_w(1)/(integral_w(2)+smalldouble))*180.d0/dpi
   line = trim(line)//trim(datastr)
   write(21,'(a)') trim(line)
   close(21)
 endif
enddo
deallocate(integral_ipe,integral_w)

end subroutine spatial_integral_w
!=============================================================================
subroutine mask_grid(ixI^L,ixO^L,w,x,patchwi,cellcount)

use mod_global_parameters

integer, intent(in)                :: ixI^L,ixO^L
double precision, intent(in)       :: x(ixI^S,1:ndim)
double precision                   :: w(ixI^S,nw+nwauxio)
logical, intent(inout)             :: patchwi(ixG^T)

double precision  ::  buff
integer                            :: ix^D,cellcount
!-----------------------------------------------------------------------------
buff=0.05d0*(xprobmax1-xprobmin1)
{do ix^DB=ixOmin^DB,ixOmax^DB\}
   if(x(ix^D,1)>xprobmin1+buff .and. x(ix^D,1)<xprobmax1-buff .and. &
      x(ix^D,2)>xprobmin2+buff .and. x(ix^D,2)<xprobmax2-buff) then
     patchwi(ix^D)=.true.
     cellcount=cellcount+1
   else
     patchwi(ix^D)=.false.
   endif
{end do\}
return
end subroutine mask_grid
!=============================================================================
function integral_grid(ixI^L,ixO^L,w,x,dvolume,dsurface,intval,patchwi)

use mod_global_parameters

integer, intent(in)                :: ixI^L,ixO^L,intval
double precision, intent(in)       :: x(ixI^S,1:ndim),dvolume(ixG^T),dsurface(ixG^T)
double precision, intent(in)       :: w(ixI^S,nw)
logical, intent(in) :: patchwi(ixG^T)

double precision, dimension(ixG^T,1:ndir) :: bvec,qvec
double precision :: current(ixG^T,7-2*ndir:3),tmp(ixG^T)
double precision :: integral_grid,mcurrent
integer :: ix^D,idirmin,idir,jdir,kdir
!-----------------------------------------------------------------------------

integral_grid=0.d0
select case(intval)
 case(1)
  ! current times sin theta (between J and B)
  if(B0field) then
    ^C&bvec(ixI^S,^C)=w(ixI^S,b^C_)+myB0_cell%w(ixI^S,^C);
  else
    ^C&bvec(ixI^S,^C)=w(ixI^S,b^C_);
  endif
  call getcurrent(w,ixI^L,ixO^L,idirmin,current)
  ! calculate Lorentz force
  qvec(ixO^S,1:ndir)=zero
  do idir=1,ndir; do jdir=idirmin,3; do kdir=1,ndir
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
     if(patchwi(ix^D)) integral_grid=integral_grid+dsqrt((^C&qvec(ix^D,^C)**2+)/&
      (^C&bvec(ix^D,^C)**2+))*dvolume(ix^D)
  {end do\}
 case(2)
  call getcurrent(w,ixI^L,ixO^L,idirmin,current)
  ! integral of current
  {do ix^DB=ixOmin^DB,ixOmax^DB\}
     mcurrent=0.d0
     do idir=7-2*ndir,3
        mcurrent=current(ix^D,idir)**2+mcurrent
     end do
     if(patchwi(ix^D))  integral_grid=integral_grid+dsqrt(mcurrent)*dvolume(ix^D)
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
     if(patchwi(ix^D)) integral_grid=integral_grid+dabs(tmp(ix^D))*&
           dvolume(ix^D)/(dsqrt(^C&bvec(ix^D,^C)**2+)+smalldouble)/dsurface(ix^D)
  {end do\}
 case default
     call mpistop("intval not defined")
end select

return
end function integral_grid
!=============================================================================
! amrvacusr.t.pfss
!=============================================================================
