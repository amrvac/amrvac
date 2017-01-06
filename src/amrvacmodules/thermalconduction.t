!##############################################################################
!> module thermalconduction.t -- thermal conduction for HD and MHD
!> 10.07.2011 developed by Chun Xia and Rony Keppens
!> 01.09.2012 moved to modules folder by Oliver Porth
!> 13.10.2013 optimized further by Chun Xia
!> 12.03.2014 implemented RKL2 super timestepping scheme to reduce iterations 
!> and improve stability and accuracy up to second order in time by Chun Xia.
!> 23.08.2014 implemented saturation and perpendicular TC by Chun Xia
!> 
!> PURPOSE: 
!> IN MHD ADD THE HEAT CONDUCTION SOURCE TO THE ENERGY EQUATION
!> S=DIV(KAPPA_i,j . GRAD_j T)
!> where KAPPA_i,j = kappa b_i b_j + kappe (I - b_i b_j)
!> b_i b_j = B_i B_j / B**2, I is the unit matrix, and i, j= 1, 2, 3 for 3D
!> IN HD ADD THE HEAT CONDUCTION SOURCE TO THE ENERGY EQUATION
!> S=DIV(kappa . GRAD T)
!> USAGE:
!> 1. in amrvacusr.t add: 
!>   a. INCLUDE:amrvacmodules/thermalconduction.t
!>   b. in subroutine initglobaldata_usr, add 
!>        unit_length=<your length unit>
!>        unit_density=<your density unit>
!>        unit_velocity=<your velocity unit>
!>        unit_temperature=<your temperature unit>
!>        unit_numberdensity=<your number density unit>
!>        unit_magneticfield=<your magnetic field unit>
!>        call init_thermalconduction
!>      to add conductivity for solar corona. Note that:
!>      kappa=kappa0*unit_temperature**3.5d0/unit_length/unit_density/unit_velocity**3 
!>      kappa0=8.d-7 erg/cm/s/K**3.5 or 8.d-12 J/m/s/K**3.5
!>      therml conductivity perpendicular to magnetic field: 
!>      kappe=kappe*unit_numberdensity**2/unit_magneticfield**2/unit_temperature**3*kappa
!>      kappe0=4.d-26 in cgs or 4.d-30 in SI
!> 2. in definitions.h :
!>    #define TCRKL2
!> 3. in the methodlist of amrvac.par add:
!>    conduction=.true.
!>    if consider saturation of thermal conduction flux, in the methodlist of amrvac.par add:
!>    TCsaturate=.true.
!>    phi coefficient of saturated flux
!>    TCphi=1.d0
!>    if consider thermal conduction perpendicular to magnetic field, in the methodlist of amrvac.par add:
!>    TCperpendicular=.true.
!=============================================================================
module mod_thermalconduction
implicit none
double precision :: kappa,kappe
integer, dimension(-1:1), save :: ixS_srl_^L, ixR_srl_^L, ixS_r_^L
integer, dimension(0:3), save :: ixR_r_^L, ixS_p_^L, ixR_p_^L
integer, dimension(-1:1^D&), save :: type_send_srl, type_recv_srl, type_send_r
integer, dimension(0:3^D&), save :: type_recv_r, type_send_p, type_recv_p
integer :: ixM^L, ixCoG^L, ixCoM^L
end module mod_thermalconduction
!=============================================================================
subroutine init_thermalconduction
use mod_thermalconduction
use mod_global_parameters

integer :: ixG^L, i^D, ic^D, inc^D
integer :: dixBCo, interpolation_order
integer :: nx^D, nxCo^D
!-----------------------------------------------------------------------------
! Spitzer thermal conductivity
kappa=8.d-7*unit_temperature**3.5d0/unit_length/unit_density/unit_velocity**3 
! thermal conductivity perpendicular to magnetic field
kappe=4.d-26*unit_numberdensity**2/unit_magneticfield**2/unit_temperature**3*kappa

ixG^L=ixG^LL;
ixM^L=ixG^L^LSUBdixB;
ixCoGmin^D=1;
ixCoGmax^D=ixGmax^D/2+dixB;
ixCoM^L=ixCoG^L^LSUBdixB;

nx^D=ixMmax^D-ixMmin^D+1;
nxCo^D=nx^D/2;

select case (typeghostfill)
case ("copy")
   interpolation_order=1
case ("linear","unlimit")
   interpolation_order=2
case default
   write (unitterm,*) "Undefined typeghostfill ",typeghostfill
   call mpistop("")
end select
dixBCo=int((dixB+1)/2)

if (dixBCo+interpolation_order-1>dixB) then
   call mpistop("interpolation order for prolongation in getbc to high")
end if

{
ixS_srl_min^D(-1)=ixMmin^D
ixS_srl_min^D(0) =ixMmin^D
ixS_srl_min^D(1) =ixMmax^D+1-dixB
ixS_srl_max^D(-1)=ixMmin^D-1+dixB
ixS_srl_max^D(0) =ixMmax^D
ixS_srl_max^D(1) =ixMmax^D

ixR_srl_min^D(-1)=1
ixR_srl_min^D(0) =ixMmin^D
ixR_srl_min^D(1) =ixMmax^D+1
ixR_srl_max^D(-1)=dixB
ixR_srl_max^D(0) =ixMmax^D
ixR_srl_max^D(1) =ixGmax^D

ixS_r_min^D(-1)=ixCoMmin^D
ixS_r_min^D(0) =ixCoMmin^D
ixS_r_min^D(1) =ixCoMmax^D+1-dixB
ixS_r_max^D(-1)=ixCoMmin^D-1+dixB
ixS_r_max^D(0) =ixCoMmax^D
ixS_r_max^D(1) =ixCoMmax^D

ixR_r_min^D(0)=1
ixR_r_min^D(1)=ixMmin^D
ixR_r_min^D(2)=ixMmin^D+nxCo^D
ixR_r_min^D(3)=ixMmax^D+1
ixR_r_max^D(0)=dixB
ixR_r_max^D(1)=ixMmin^D-1+nxCo^D
ixR_r_max^D(2)=ixMmax^D
ixR_r_max^D(3)=ixGmax^D

ixS_p_min^D(0)=ixMmin^D-(interpolation_order-1)
ixS_p_min^D(1)=ixMmin^D-(interpolation_order-1)
ixS_p_min^D(2)=ixMmin^D+nxCo^D-dixBCo-(interpolation_order-1)
ixS_p_min^D(3)=ixMmax^D+1-dixBCo-(interpolation_order-1)
ixS_p_max^D(0)=ixMmin^D-1+dixBCo+(interpolation_order-1)
ixS_p_max^D(1)=ixMmin^D-1+nxCo^D+dixBCo+(interpolation_order-1)
ixS_p_max^D(2)=ixMmax^D+(interpolation_order-1)
ixS_p_max^D(3)=ixMmax^D+(interpolation_order-1)

ixR_p_min^D(0)=ixCoMmin^D-dixBCo-(interpolation_order-1)
ixR_p_min^D(1)=ixCoMmin^D-(interpolation_order-1)
ixR_p_min^D(2)=ixCoMmin^D-dixBCo-(interpolation_order-1)
ixR_p_min^D(3)=ixCoMmax^D+1-(interpolation_order-1)
ixR_p_max^D(0)=dixB+(interpolation_order-1)
ixR_p_max^D(1)=ixCoMmax^D+dixBCo+(interpolation_order-1)
ixR_p_max^D(2)=ixCoMmax^D+(interpolation_order-1)
ixR_p_max^D(3)=ixCoMmax^D+dixBCo+(interpolation_order-1)
\}

{do i^DB=-1,1\}
   if (i^D==0|.and.) cycle
   call get_bc_comm_type(type_send_srl(i^D),ixS_srl_^L(i^D),ixG^L)
   call get_bc_comm_type(type_recv_srl(i^D),ixR_srl_^L(i^D),ixG^L)
   call get_bc_comm_type(type_send_r(i^D),ixS_r_^L(i^D),ixCoG^L)
   {do ic^DB=1+int((1-i^DB)/2),2-int((1+i^DB)/2)
      inc^DB=2*i^DB+ic^DB\}
      call get_bc_comm_type(type_recv_r(inc^D),ixR_r_^L(inc^D),ixG^L)
      call get_bc_comm_type(type_send_p(inc^D),ixS_p_^L(inc^D),ixG^L)
      call get_bc_comm_type(type_recv_p(inc^D),ixR_p_^L(inc^D),ixCoG^L)
   {end do\}
{end do\}
contains
!=============================================================================
subroutine get_bc_comm_type(comm_type,ix^L,ixG^L)

integer, intent(inout) :: comm_type
integer, intent(in) :: ix^L, ixG^L

integer, dimension(ndim+1) :: size, subsize, start
!-----------------------------------------------------------------------------
^D&size(^D)=ixGmax^D;
size(ndim+1)=nw
^D&subsize(^D)=ixmax^D-ixmin^D+1;
subsize(ndim+1)=1
^D&start(^D)=ixmin^D-1;
start(ndim+1)=e_-1

call MPI_TYPE_CREATE_SUBARRAY(ndim+1,size,subsize,start,MPI_ORDER_FORTRAN, &
                              MPI_DOUBLE_PRECISION,comm_type,ierrmpi)
call MPI_TYPE_COMMIT(comm_type,ierrmpi)

end subroutine get_bc_comm_type
!=============================================================================
end subroutine init_thermalconduction
!=============================================================================
subroutine thermalconduction
! Meyer 2012 MNRAS 422,2102
use mod_thermalconduction
use mod_global_parameters


double precision :: omega1,cmu,cmut,cnu,cnut
double precision, allocatable :: bj(:)
double precision, save :: qdt
integer:: iigrid, igrid,j
integer, save :: s
logical :: evenstep, first=.true.
!-----------------------------------------------------------------------------
if(e_<1) call mpistop("Thermal conduction requires by e_>0!")

if(first) then
  if(dt/dtimpl < 0.5d0) then
    s=1
  else if(dt/dtimpl < 2.d0) then
    s=2
  else
    s=ceiling((dsqrt(9.d0+8.d0*dt/dtimpl)-1.d0)/2.d0)
    ! only use odd s number
    s=s/2*2+1
  endif
  qdt=dt*0.5d0
  if(mype==0 .and. .false.) write(*,*) 'supertime steps:',s,' normal subcycles:',&
                              ceiling(dt/dtimpl/2.d0),'time step ratio:',dt/dtimpl
  first=.false.
else
  first=.true.
end if

do iigrid=1,igridstail; igrid=igrids(iigrid);
   allocate(pw1(igrid)%w(ixG^T,1:nw))
   allocate(pw2(igrid)%w(ixG^T,1:nw))
   allocate(pw3(igrid)%w(ixG^T,1:nw))
   pw1(igrid)%w=pw(igrid)%w
   pw2(igrid)%w=pw(igrid)%w
   pw3(igrid)%w=pwold(igrid)%w
end do

allocate(bj(0:s))
bj(0)=1.d0/3.d0
bj(1)=bj(0)
if(s>1) then
  omega1=4.d0/dble(s**2+s-2)
  cmut=omega1/3.d0
else
  cmut=1.d0
endif

!$OMP PARALLEL DO PRIVATE(igrid)
do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
  if (.not.slab) mygeo => pgeo(igrid)
  if (B0field) then
     myB0_cell => pB0_cell(igrid)
     {^D&myB0_face^D => pB0_face^D(igrid)\}
  end if
  typelimiter=typelimiter1(node(plevel_,igrid))
  typegradlimiter=typegradlimiter1(node(plevel_,igrid))
  ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
  call evolve_step1(cmut,qdt,ixG^LL,ixM^LL,pw1(igrid)%w,pw(igrid)%w,&
                    px(igrid)%x,pw3(igrid)%w)
end do
!$OMP END PARALLEL DO
call getbce(t,pw1)
if(s==1) then
  do iigrid=1,igridstail; igrid=igrids(iigrid);
    pw(igrid)%w(ixG^T,e_)=pw1(igrid)%w(ixG^T,e_)
    deallocate(pw1(igrid)%w)
    deallocate(pw2(igrid)%w)
    deallocate(pw3(igrid)%w)
  end do
  return
endif
evenstep=.true.
do j=2,s
  bj(j)=dble(j**2+j-2)/dble(2*j*(j+1))
  cmu=dble(2*j-1)/dble(j)*bj(j)/bj(j-1)
  cmut=omega1*cmu
  cnu=dble(1-j)/dble(j)*bj(j)/bj(j-2)
  cnut=(bj(j-1)-1.d0)*cmut
  if(evenstep) then
!$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
      if(.not.slab) mygeo => pgeo(igrid)
      if(B0field) then
        myB0_cell => pB0_cell(igrid)
        {^D&myB0_face^D => pB0_face^D(igrid)\}
      end if
      typelimiter=typelimiter1(node(plevel_,igrid))
      typegradlimiter=typegradlimiter1(node(plevel_,igrid))
      ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
      call evolve_stepj(cmu,cmut,cnu,cnut,qdt,ixG^LL,ixM^LL,pw1(igrid)%w,&
                        pw2(igrid)%w,pw(igrid)%w,px(igrid)%x,pw3(igrid)%w)
    end do
!$OMP END PARALLEL DO
    call getbce(t,pw2)
    evenstep=.false.
  else
!$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
      if(.not.slab) mygeo => pgeo(igrid)
      if(B0field) then
         myB0_cell => pB0_cell(igrid)
        {^D&myB0_face^D => pB0_face^D(igrid)\}
      end if
      typelimiter=typelimiter1(node(plevel_,igrid))
      typegradlimiter=typegradlimiter1(node(plevel_,igrid))
      ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
      call evolve_stepj(cmu,cmut,cnu,cnut,qdt,ixG^LL,ixM^LL,pw2(igrid)%w,&
                        pw1(igrid)%w,pw(igrid)%w,px(igrid)%x,pw3(igrid)%w)
    end do
!$OMP END PARALLEL DO
    call getbce(t,pw1)
    evenstep=.true.
  end if 
end do
if(evenstep) then
  do iigrid=1,igridstail; igrid=igrids(iigrid);
    pw(igrid)%w(ixG^T,e_)=pw1(igrid)%w(ixG^T,e_)
    deallocate(pw1(igrid)%w)
    deallocate(pw2(igrid)%w)
    deallocate(pw3(igrid)%w)
  end do 
else
  do iigrid=1,igridstail; igrid=igrids(iigrid);
    pw(igrid)%w(ixG^T,e_)=pw2(igrid)%w(ixG^T,e_)
    deallocate(pw1(igrid)%w)
    deallocate(pw2(igrid)%w)
    deallocate(pw3(igrid)%w)
  end do 
end if
deallocate(bj)

end subroutine thermalconduction
!=============================================================================
subroutine evolve_stepj(qcmu,qcmut,qcnu,qcnut,qdt,ixI^L,ixO^L,w1,w2,w,x,wold)

use mod_global_parameters

integer, intent(in) :: ixI^L,ixO^L
double precision, intent(in) :: qcmu,qcmut,qcnu,qcnut,qdt
double precision, intent(in) :: w1(ixI^S,1:nw),w(ixI^S,1:nw),wold(ixI^S,1:nw)
double precision, intent(in) :: x(ixI^S,1:ndim)
double precision, intent(inout) :: w2(ixI^S,1:nw)

double precision :: tmp(ixI^S),tmp1(ixI^S),tmp2(ixI^S)
!-----------------------------------------------------------------------------
{^IFMHDPHYS
call heatconduct_mhd(tmp,tmp1,tmp2,ixI^L,ixO^L,w1,x)
}
{^IFHDPHYS
call heatconduct_hd(tmp,tmp1,tmp2,ixI^L,ixO^L,w1,x)
}
w2(ixO^S,e_)=qcmu*w1(ixO^S,e_)+qcnu*w2(ixO^S,e_)+(1.d0-qcmu-qcnu)*w(ixO^S,e_)&
            +qcmut*qdt*tmp(ixO^S)+qcnut*wold(ixO^S,e_)

if(fixsmall) call smallvalues(w2,x,ixI^L,ixO^L,'evolve_stepj')
end subroutine evolve_stepj
!=============================================================================
subroutine evolve_step1(qcmut,qdt,ixI^L,ixO^L,w1,w,x,wold)

use mod_global_parameters

integer, intent(in) :: ixI^L,ixO^L
double precision, intent(in) :: qcmut, qdt, w(ixI^S,1:nw), x(ixI^S,1:ndim)
double precision, intent(out) ::w1(ixI^S,1:nw),wold(ixI^S,1:nw)

double precision :: tmp(ixI^S),tmp1(ixI^S),tmp2(ixI^S),Te(ixI^S)
integer :: lowindex(ndim), ix^D
!-----------------------------------------------------------------------------
{^IFMHDPHYS
call heatconduct_mhd(tmp,tmp1,tmp2,ixI^L,ixO^L,w,x)
}
{^IFHDPHYS
call heatconduct_hd(tmp,tmp1,tmp2,ixI^L,ixO^L,w,x)
}

wold(ixO^S,e_)=qdt*tmp(ixO^S)
! update internal energy
tmp1(ixO^S) = tmp1(ixO^S) + qcmut*wold(ixO^S,e_)

! ensure you never trigger negative pressure 
! hence code up energy change with respect to kinetic and magnetic
! part(nonthermal)
if(smallT>0.d0) then
  Te(ixO^S)=tmp1(ixO^S)*(eqpar(gamma_)-1.d0)/w(ixO^S,rho_)
endif
if(strictsmall) then
  if(smallT>0.d0 .and. any(Te(ixO^S)<smallT)) then
    lowindex=minloc(Te(ixO^S))
    ^D&lowindex(^D)=lowindex(^D)+ixOmin^D-1;
    write(*,*)'too small temperature = ',minval(Te(ixO^S)),'at x=',&
   x(^D&lowindex(^D),1:ndim),lowindex,' with limit=',smallT,&
   ' on time=',t,' step=',it,' where density=',w(^D&lowindex(^D),rho_),&
   ' velocity=',&
    dsqrt((^C&w(^D&lowindex(^D),m^C_)**2+)/w(^D&lowindex(^D),rho_)**2)
{^IFMHDPHYS
    write(*,*)'magnetic field=',dsqrt((^C&w(^D&lowindex(^D),b^C_)**2+))
}
    call mpistop("==evolve_step1: too small temperature==")
  end if
  if(any(tmp1(ixO^S)<smalle)) then
    lowindex=minloc(tmp1(ixO^S))
    ^D&lowindex(^D)=lowindex(^D)+ixOmin^D-1;
    write(*,*)'too small internal energy = ',minval(tmp1(ixO^S)),'at x=',&
   x(^D&lowindex(^D),1:ndim),lowindex,' with limit=',smalle,&
   ' on time=',t,' step=',it,' where density=',w(^D&lowindex(^D),rho_),&
   ' velocity=',&
    dsqrt((^C&w(^D&lowindex(^D),m^C_)**2+)/w(^D&lowindex(^D),rho_)**2)
{^IFMHDPHYS
    write(*,*)'magnetic field=',dsqrt((^C&w(^D&lowindex(^D),b^C_)**2+))
}
    call mpistop("==evolve_step1: too small internal energy==")
  end if
  w1(ixO^S,e_) = tmp2(ixO^S)+tmp1(ixO^S)
else
 {do ix^DB=ixOmin^DB,ixOmax^DB\}
    if(smallT>0.d0) then
      if(Te(ix^D)<smallT) then
        w1(ix^D,e_) = tmp2(ix^D)+smallT*w(ix^D,rho_)/(eqpar(gamma_)-1.d0)
      else
        w1(ix^D,e_) = tmp2(ix^D)+tmp1(ix^D)
      end if
    else
      if(tmp1(ix^D)<smalle) then
        w1(ix^D,e_) = tmp2(ix^D)+smalle
      else
        w1(ix^D,e_) = tmp2(ix^D)+tmp1(ix^D)
      end if
    end if
 {end do\}
end if
!if(fixsmall) call smallvalues(w1,x,ixI^L,ixO^L,'evolve_step1')

end subroutine evolve_step1
!=============================================================================
subroutine addsource_heatconduct_mhd(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

use mod_global_parameters

integer, intent(in) :: ixI^L,ixO^L,iw^LIM
double precision, intent(in) :: qdt,qtC,qt, x(ixI^S,1:ndim),wCT(ixI^S,1:nw)
double precision, intent(inout) ::w(ixI^S,1:nw)
double precision :: tmp(ixI^S),tmp1(ixI^S),tmp2(ixI^S)
integer :: ix^D
integer, dimension(ndim)       :: lowindex
!------------------------------------------------------------------------------

if(e_<1) call mpistop("Heat conduction requires by e_>0!")
!heat conduct rate tmp, old internal energy tmp1, old nonthermal energy tmp2
call heatconduct_mhd(tmp,tmp1,tmp2,ixI^L,ixO^L,w,x)
! update internal energy
tmp1(ixO^S) = tmp1(ixO^S) + qdt*tmp(ixO^S)

! code up energy change with respect to kinetic and magnetic part(nonthermal)
if(strictsmall) then
  if(any(tmp1(ixO^S)<smalle)) then
    lowindex=minloc(tmp1(ixO^S))
    ^D&lowindex(^D)=lowindex(^D)+ixOmin^D-1;
    write(*,*)'too small internal energy = ',minval(tmp1(ixO^S)),'at x=',&
  x(^D&lowindex(^D),1:ndim),lowindex,' with limit=',smalle,&
 ' on time=',t,' step=',it,' where density=',w(^D&lowindex(^D),rho_),&
 ' velocity=',&
 dsqrt((^C&w(^D&lowindex(^D),m^C_)**2+)/w(^D&lowindex(^D),rho_)**2),&
    ' magnetic field=',dsqrt((^C&w(^D&lowindex(^D),b^C_)**2+))
    call mpistop("=strictsmall in sourceheatcond_mhd: too small energy==")
  else
    w(ixO^S,e_) = tmp2(ixO^S)+tmp1(ixO^S)
  endif
else
 {do ix^DB=ixOmin^DB,ixOmax^DB\}
    if(tmp1(ix^D)>smalle) then
      w(ix^D,e_) = tmp2(ix^D)+tmp1(ix^D)
    else
      w(ix^D,e_) = tmp2(ix^D)+smalle
    end if
 {end do\}
endif

end subroutine addsource_heatconduct_mhd
!=============================================================================
subroutine heatconduct_mhd(tmp,tmp1,tmp2,ixI^L,ixO^L,w,x)
use mod_thermalconduction
use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L
double precision, intent(in) ::  x(ixI^S,1:ndim), w(ixI^S,1:nw)
!! tmp store the heat conduction energy changing rate
double precision, intent(out) :: tmp(ixI^S),tmp1(ixI^S),tmp2(ixI^S)

double precision, dimension(ixI^S,1:ndir) :: mf,qvec
double precision, dimension(ixI^S,1:ndim) :: gradT,qvecsat
double precision, dimension(:^D&,:), allocatable :: qvec_per, qvec_max
double precision, dimension(ixI^S) :: B2inv, BgradT,cs3,qflux,qsatflux
integer, dimension(ndim) :: lowindex
integer :: ix^L,idims,ix^D
logical :: Bnull(ixI^S)
!-----------------------------------------------------------------------------
ix^L=ixO^L^LADD2;
if (ixI^L^LTix^L|.or.|.or.) &
   call mpistop("Need two extra layers for thermal conduction")

! store kinetic+magnetic energy before addition of heat conduction source
tmp2(ixI^S)=half*( (^C&w(ixI^S,m^C_)**2+)/w(ixI^S,rho_)+ &
  (^C&w(ixI^S,b^C_)**2+) )

! Calculate einternal=e-0.5*(2ek+2eb)
tmp1(ixI^S)=w(ixI^S,e_)-tmp2(ixI^S)
! Clip off negative pressure if smallp is set
if(strictsmall) then
   if(any(tmp1(ixI^S)<smalle)) then
     lowindex=minloc(tmp1(ixI^S))
     ^D&lowindex(^D)=lowindex(^D)+ixImin^D-1;
     write(*,*)'too low internal energy = ',minval(tmp1(ixI^S)),' at x=',&
     x(^D&lowindex(^D),1:ndim),lowindex,' with limit=',smalle,' on time=',t,&
     ' step=',it,' where density=',w(^D&lowindex(^D),rho_),' velocity=',&
     dsqrt((^C&w(^D&lowindex(^D),m^C_)**2+)/w(^D&lowindex(^D),rho_)**2),&
     ' magnetic field=',dsqrt((^C&w(^D&lowindex(^D),b^C_)**2+))
     call mpistop("=== strictsmall in heatcond_mhd: low internal energy ===")
   end if
else
{do ix^DB=ixImin^DB,ixImax^DB\}
   if(tmp1(ix^D)<smalle) then
      tmp1(ix^D)=smalle
   end if
{end do\}
end if
! compute the temperature
tmp(ixI^S)=tmp1(ixI^S)*(eqpar(gamma_)-one)/w(ixI^S,rho_)
if(smallT>0.d0) then
  if(strictsmall) then
     if(any(tmp(ixI^S)<smallT)) then
       lowindex=minloc(tmp(ixI^S))
       ^D&lowindex(^D)=lowindex(^D)+ixImin^D-1;
       write(*,*)'too low temperature = ',minval(tmp(ixI^S)),' at x=',&
       x(^D&lowindex(^D),1:ndim),lowindex,' with limit=',smallT,' on time=',t,&
       ' step=',it,' where density=',w(^D&lowindex(^D),rho_),' velocity=',&
       dsqrt((^C&w(^D&lowindex(^D),m^C_)**2+)/w(^D&lowindex(^D),rho_)**2),&
       ' magnetic field=',dsqrt((^C&w(^D&lowindex(^D),b^C_)**2+))
       call mpistop("=== strictsmall in heatcond_mhd: low temperature ===")
     end if
  else
  {do ix^DB=ixImin^DB,ixImax^DB\}
     if(tmp(ix^D)<smallT) then
        tmp(ix^D)=smallT
     end if
  {end do\}
  end if
end if
do idims=1,ndim
  ! idirth component of gradient of temperature at cell center
  select case(typegrad)
  case("central")
    ix^L=ixI^L^LSUB1;
    call gradient(tmp,ixI^L,ix^L,idims,B2inv)
  case("limited")
    ix^L=ixI^L^LSUB2;
    call gradientS(tmp,ixI^L,ix^L,idims,B2inv)
  end select
  ! get grad T in all directions
  gradT(ix^S,idims)=B2inv(ix^S)
end do
! B
if(B0field) then
 ^C&mf(ix^S,^C)=w(ix^S,b^C_)+myB0_cell%w(ix^S,^C);
else
 ^C&mf(ix^S,^C)=w(ix^S,b^C_);
end if
! B^-2
B2inv(ix^S)= ^C&mf(ix^S,^C)**2+
Bnull(ixI^S)=.false.
where(B2inv(ix^S)/=0.d0)
  B2inv(ix^S)=1.d0/B2inv(ix^S)
elsewhere
  Bnull(ix^S)=.true.
end where

BgradT(ix^S)=(^D&mf(ix^S,^D)*gradT(ix^S,^D)+)*B2inv(ix^S)
cs3(ix^S)=kappa*dsqrt(tmp(ix^S))*tmp(ix^S)**2*BgradT(ix^S)

do idims=1,ndim
  qvec(ix^S,idims)=mf(ix^S,idims)*cs3(ix^S)
end do

if(TCsaturate) then
  ! Cowie and Mckee 1977 ApJ, 211, 135
  cs3(ix^S)=dsqrt(tmp(ix^S)**3)
  ! unsigned saturated TC flux = 5 phi rho c**3
  qsatflux(ix^S)=5.d0*TCphi*w(ix^S,rho_)*cs3(ix^S)
  ! strength of classic TC flux
  qflux(ix^S)=dsqrt(^D&qvec(ix^S,^D)**2+)
  ! sign(b * Grad Te) 5 phi rho c**3 Bi/B 
  {do ix^DB=ixmin^DB,ixmax^DB\}
    if(.false. .and. qflux(ix^D)>qsatflux(ix^D) .and. idims==1) write(*,*)& 
      'it',it,' ratio=',qflux(ix^D)/qsatflux(ix^D),' TC saturated at ',&
      x(ix^D,:),' rho',w(ix^D,rho_),' Te',tmp(ix^D)
    if(qflux(ix^D)>qsatflux(ix^D)) then
    ! saturated TC flux = sign(b * Grad Te) 5 phi rho c**3
      qsatflux(ix^D)=sign(1.d0,BgradT(ix^D))*qsatflux(ix^D)*dsqrt(B2inv(ix^D))
      do idims=1,ndim
        qvec(ix^D,idims)=qsatflux(ix^D)*mf(ix^D,idims)
      end do
    end if
  {end do\}
end if

if(TCperpendicular) then
! van der Linden and Goossens 1991 SoPh 131, 79; Orlando et al 2008 ApJ 678, 274
  allocate(qvec_per(ixI^S,1:ndim), qvec_max(ixI^S,1:ndim))
  do idims=1,ndim
    ! q_per = kappe n^2 B^-2 Te^-0.5 (Grad Te - (e_b . Grad Te )e_b) 
    qvec_per(ix^S,idims)=kappe*w(ix^S,rho_)**2*B2inv(ix^S)/dsqrt(tmp(ix^S))&
                         *(gradT(ix^S,idims)-BgradT(ix^S)*mf(ix^S,idims))
  end do
  ! maximal thermal conducting (saturated) flux
  do idims=1,ndim
    qvec_max(ix^S,idims)=kappa*dsqrt(tmp(ix^S)**5)*gradT(ix^S,idims)
  end do
  if(TCsaturate) then
    qsatflux(ix^S)=5.d0*TCphi*w(ix^S,rho_)*cs3(ix^S)
    qflux(ix^S)=dsqrt(^D&qvec_max(ix^S,^D)**2+)
    {do ix^DB=ixmin^DB,ixmax^DB\}
      if(.false. .and. qflux(ix^D)>qsatflux(ix^D) .and. idims==1) write(*,*) & 
        'it',it,' ratio=',qflux(ix^D)/qsatflux(ix^D),' TC_PER saturated at ',&
        x(ix^D,:),' rho',w(ix^D,rho_),' Te',tmp(ix^D)
      if(qflux(ix^D)>qsatflux(ix^D)) then
        qsatflux(ix^D)=qsatflux(ix^D)/qflux(ix^D)
        do idims=1,ndim
          qvec_max(ix^D,idims)=qsatflux(ix^D)*qvec_max(ix^D,idims)
        end do
      end if
    {end do\}
  end if
  ! maximal thermal conducting flux perpendicular to magnetic field
  qvec_max(ix^S,1:ndim)=qvec_max(ix^S,1:ndim)-qvec(ix^S,1:ndim)
  qsatflux(ix^S)= ^D&qvec_max(ix^S,^D)**2+
  qflux(ix^S)= ^D&qvec_per(ix^S,^D)**2+
  {do ix^DB=ixmin^DB,ixmax^DB\}
     if(qflux(ix^D)>qsatflux(ix^D)) then
       qvec(ix^D,1:ndim)=qvec(ix^D,1:ndim)+qvec_max(ix^D,1:ndim)
     else
       qvec(ix^D,1:ndim)=qvec(ix^D,1:ndim)+qvec_per(ix^D,1:ndim)
     end if   
  {end do\}
end if

{do ix^DB=ixmin^DB,ixmax^DB\}
  if(Bnull(ix^D)) then
    qvec(ix^D,1:ndim)=kappa*dsqrt(tmp(ix^D)**5)*gradT(ix^D,1:ndim)
    if(TCsaturate) then
      ! unsigned saturated TC flux = 5 phi rho c**3
      qsatflux(ix^D)=5.d0*TCphi*w(ix^D,rho_)*cs3(ix^D)
      qflux(ix^D)=dsqrt(sum(qvec(ix^D,1:ndim)**2))
      if(qflux(ix^D)>qsatflux(ix^D)) then
        qsatflux(ix^D)=qsatflux(ix^D)/qflux(ix^D)
        do idims=1,ndim
          qvec(ix^D,idims)=qsatflux(ix^D)*qvec(ix^D,idims)
        end do
      end if
    end if
  end if  
{end do\}

select case(typediv)
   case("central")
      call divvector(qvec,ixI^L,ixO^L,tmp)
   case("limited")
      call divvectorS(qvec,ixI^L,ixO^L,tmp)
end select

end subroutine heatconduct_mhd
!=============================================================================
subroutine getdt_heatconduct_mhd(w,ixG^L,ix^L,dtnew,dx^D,x)

!Check diffusion time limit dt < dtTCpar*dx_i**2/((gamma-1)*kappa_i/rho)
!where                      kappa_i=kappa*B_i**2/B**2
!and                        T=p/rho
use mod_thermalconduction
use mod_global_parameters

integer, intent(in) :: ixG^L, ix^L
double precision, intent(in) :: dx^D, x(ixG^S,1:ndim)
! note that depending on strictsmall etc, w values may change 
! through call to getpthermal
double precision, intent(inout) :: w(ixG^S,1:nw), dtnew

double precision :: dxinv(1:ndim),mf(ixG^S,1:ndir)
double precision :: tmp2(ixG^S),tmp(ixG^S),Te(ixG^S),B2inv(ixG^S)
double precision :: dtdiff_tcond, dtdiff_tsat
integer          :: idim,ix^D
integer, dimension(ndim)       :: lowindex
!-----------------------------------------------------------------------------
^D&dxinv(^D)=one/dx^D;

! Calculate pressure=(gamma-1)*(e-0.5*(2ek+2eb))
tmp(ix^S)=(eqpar(gamma_)-one)*(w(ix^S,e_)- &
       half*(({^C&w(ix^S,m^C_)**2+})/w(ix^S,rho_)&
       +{ ^C&w(ix^S,b^C_)**2+}))
! Clip off negative pressure if smallp is set
if(strictsmall) then
  if(any(tmp(ix^S)<minp)) then
    lowindex=minloc(tmp(ix^S))
    ^D&lowindex(^D)=lowindex(^D)+ixmin^D-1;
    write(*,*)'low pressure = ',minval(tmp(ix^S)),' at x=',&
    x(^D&lowindex(^D),1:ndim),lowindex,' with limit=',minp,' on time=',t,&
    ' step=',it
   call mpistop("=== strictsmall in getdt_heatconduct_mhd: low pressure ===")
  end if
else
{do ix^DB=ixmin^DB,ixmax^DB\}
   if(tmp(ix^D)<minp) then
      tmp(ix^D)=minp
   end if
{end do\}
end if
!temperature
Te(ix^S)=tmp(ix^S)/w(ix^S,rho_)
!kappa_i
tmp(ix^S)=kappa*dsqrt(Te(ix^S)**5)
!(gamma-1)*kappa_i/rho
tmp(ix^S)=(eqpar(gamma_)-one)*tmp(ix^S)/w(ix^S,rho_)

! B
if(B0field) then
 ^C&mf(ix^S,^C)=w(ix^S,b^C_)+myB0_cell%w(ix^S,^C);
else
 ^C&mf(ix^S,^C)=w(ix^S,b^C_);
end if
! B^-2
B2inv(ix^S)= ^C&mf(ix^S,^C)**2+
where(B2inv(ix^S)/=0.d0)
  B2inv(ix^S)=1.d0/B2inv(ix^S)
end where
do idim=1,ndim
   ! B_i**2/B**2
   where(B2inv(ix^S)/=0.d0)
     tmp2(ix^S)=mf(ix^S,idim)**2*B2inv(ix^S)
   elsewhere
     tmp2(ix^S)=1.d0
   end where
   ! dt< dtTCpar * dx_idim**2/((gamma-1)*kappa_i/rho*B_i**2/B**2)
   dtdiff_tcond=dtTcpar/maxval(tmp(ix^S)*tmp2(ix^S)*dxinv(idim)**2)
   if(TCsaturate) then
     ! dt< dtTCpar* dx_idim**2/((gamma-1)*sqrt(Te)*5*phi)
     ! with an empirical coefficient dx_idim
     dtdiff_tsat=dtTCpar/maxval((eqpar(gamma_)-1.d0)*dsqrt(Te(ix^S))*&
                 5.d0*TCphi*dxinv(idim)**2)
     ! choose the slower flux (bigger time step) between classic and saturated
     dtdiff_tcond=max(dtdiff_tcond,dtdiff_tsat)
   end if
   ! limit the time step
   dtnew=min(dtnew,dtdiff_tcond)
end do
dtnew=dtnew/dble(ndim)

end subroutine getdt_heatconduct_mhd
!=============================================================================
subroutine addsource_heatconduct_hd(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

use mod_global_parameters

integer, intent(in) :: ixI^L,ixO^L,iw^LIM
double precision, intent(in) :: qdt,qtC,qt, x(ixI^S,1:ndim),wCT(ixI^S,1:nw)
double precision, intent(inout) ::w(ixI^S,1:nw)

double precision :: tmp(ixI^S),tmp1(ixI^S),tmp2(ixI^S)
integer:: ix^L,idim,ix^D
integer, dimension(ndim)       :: lowindex
!------------------------------------------------------------------------------
if(e_<1) call mpistop("Thermal conduction requires e_>0!")
!get heat conduct rate in tmp, old thermal energy in tmp1, old nonthermal energy in tmp2
call heatconduct_hd(tmp,tmp1,tmp2,ixI^L,ixO^L,w,x)
! update internal energy
tmp1(ixO^S) = tmp1(ixO^S) + qdt*tmp(ixO^S)

if(strictsmall) then
  if(any(tmp1(ixO^S)<smalle)) then
    lowindex=minloc(tmp1(ixO^S))
    ^D&lowindex(^D)=lowindex(^D)+ixOmin^D-1;
    write(*,*)'too small internal energy = ',minval(tmp1(ixO^S)),' at x=',&
    x(^D&lowindex(^D),1:ndim),lowindex,' with limit=',smalle,&
    ' on time=',t
    call mpistop("=== strictsmall in heatconduct: too small internal energy ===")
  end if
endif

! ensure you never trigger negative pressure 
! hence code up energy change with respect to kinetic and magnetic part
{do ix^DB=ixOmin^DB,ixOmax^DB\}
 if(tmp1(ix^D)>smalle) then
   w(ix^D,e_) = tmp2(ix^D)+tmp1(ix^D)
 else
   w(ix^D,e_) = tmp2(ix^D)+smalle
 endif
{end do\}

return
end subroutine addsource_heatconduct_hd
!=============================================================================
subroutine heatconduct_hd(tmp,tmp1,tmp2,ixI^L,ixO^L,w,x)
use mod_thermalconduction
use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L
double precision, intent(in) ::  x(ixI^S,1:ndim), w(ixI^S,1:nw)
!! tmp store the heat conduction energy changing rate
double precision, intent(out) :: tmp(ixI^S),tmp1(ixI^S),tmp2(ixI^S)
double precision :: qvec(ixI^S,1:ndir),Te(ixI^S),qflux(ixI^S),qsatflux(ixI^S)
integer:: ix^L,idims,ix^D
integer, dimension(ndim)       :: lowindex
!-----------------------------------------------------------------------------

ix^L=ixO^L^LADD2;
if (ixI^L^LTix^L|.or.|.or.) &
   call mpistop("Need two extra layers for thermal conduction")

! store old kinetic energy
tmp2(ixI^S)=half*(^C&w(ixI^S,m^C_)**2+ )/w(ixI^S,rho_)
! Calculate pressure=(gamma-1)*(e-0.5*(2ek+2eb))
tmp(ixI^S)=(eqpar(gamma_)-one)*(w(ixI^S,e_)-tmp2(ixI^S))
! Clip off negative pressure if smallp is set
if(strictsmall) then
  if(any(tmp(ixI^S)<minp)) then
    lowindex=minloc(tmp(ixI^S))
    write(*,*)'low pressure = ',minval(tmp(ixI^S)),' at x=',&
    x(^D&lowindex(^D),1:ndim),lowindex,' with limit=',minp,' on time=',t
    call mpistop("=== strictsmall in heatconduct: low pressure ===")
  end if
else
   {do ix^DB=ixImin^DB,ixImax^DB\}
     if(tmp(ix^D)<minp) then
      tmp(ix^D)=minp
     end if
   {end do\}
end if
! store old internal energy
tmp1(ixI^S)=tmp(ixI^S)/(eqpar(gamma_)-one)

! compute temperature before source addition
Te(ixI^S)=tmp(ixI^S)/w(ixI^S,rho_)
if(smallT>0.d0) then
  if(strictsmall) then
     if(any(Te(ixI^S)<smallT)) then
       lowindex=minloc(Te(ixI^S))
       ^D&lowindex(^D)=lowindex(^D)+ixImin^D-1;
       write(*,*)'too low temperature = ',minval(Te(ixI^S)),' at x=',&
       x(^D&lowindex(^D),1:ndim),lowindex,' with limit=',smallT,' on time=',t,&
       ' step=',it,' where density=',w(^D&lowindex(^D),rho_),' velocity=',&
       dsqrt((^C&w(^D&lowindex(^D),m^C_)**2+)/w(^D&lowindex(^D),rho_)**2)
       call mpistop("=== strictsmall in heatcond_hd: low temperature ===")
     end if
  else
  {do ix^DB=ixImin^DB,ixImax^DB\}
     if(Te(ix^D)<smallT) then
        Te(ix^D)=smallT
     end if
  {end do\}
  end if
end if
! compute grad T and store grad T vector
do idims=1,ndim
  ! idirth component of gradient of temperature at cell center
   select case(typegrad)
   case("central")
     ix^L=ixI^L^LSUB1;
     call gradient(Te,ixI^L,ix^L,idims,tmp)
   case("limited")
     ix^L=ixI^L^LSUB2;
     call gradientS(Te,ixI^L,ix^L,idims,tmp)
   end select
   qvec(ix^S,idims)=tmp(ix^S)*kappa*dsqrt(Te(ix^S)**5)
end do

if(TCsaturate) then
  ! unsigned saturated TC flux = 5 phi rho c**3
  qsatflux(ix^S)=5.d0*TCphi*w(ix^S,rho_)*dsqrt(Te(ix^S)**3)
  qflux(ix^S)=dsqrt(^D&qvec(ix^S,^D)**2+)
  {do ix^DB=ixmin^DB,ixmax^DB\}
    if(qflux(ix^D)>qsatflux(ix^D)) then
      qsatflux(ix^D)=qsatflux(ix^D)/qflux(ix^D)
      do idims=1,ndim
        qvec(ix^D,idims)=qsatflux(ix^D)*qvec(ix^D,idims)
      end do
    end if
  {end do\}
end if
select case(typediv)
   case("central")
      call divvector(qvec,ixI^L,ixO^L,tmp)
   case("limited")
      call divvectorS(qvec,ixI^L,ixO^L,tmp)
end select
end subroutine heatconduct_hd
!=============================================================================
subroutine getdt_heatconduct_hd(w,ixG^L,ix^L,dtnew,dx^D,x)

! Check diffusion time limit dt < dtTCpar * dx_i**2 / ((gamma-1)*kappa_i/rho)
use mod_thermalconduction
use mod_global_parameters

integer, intent(in) :: ixG^L, ix^L
double precision, intent(in) :: dx^D, x(ixG^S,1:ndim)
! note that depending on strictsmall etc, w values may change 
! through call to getpthermal
double precision, intent(inout) :: w(ixG^S,1:nw), dtnew

double precision :: dxinv(1:ndim), tmp(ixG^S), Te(ixG^S)
double precision :: dtdiff_tcond,dtdiff_tsat
integer          :: idim,ix^D
integer, dimension(ndim)       :: lowindex
!-----------------------------------------------------------------------------
^D&dxinv(^D)=one/dx^D;

! Determine pressure and then the direction independent part of
tmp(ix^S)=( ^C&w(ix^S,m^C_)**2+ )/w(ix^S,rho_)
! Calculate pressure=(gamma-1)*(e-0.5*(2ek))
tmp(ix^S)=(eqpar(gamma_)-one)*(w(ix^S,e_)-half*tmp(ix^S))
! Clip off negative pressure if smallp is set
if(strictsmall) then
  if(any(tmp(ix^S)<minp)) then
    lowindex=minloc(tmp(ix^S))
    ^D&lowindex(^D)=lowindex(^D)+ixmin^D-1;
    write(*,*)'low pressure = ',minval(tmp(ix^S)),' at x=',&
    x(^D&lowindex(^D),1:ndim),lowindex,' with limit=',minp,' on time=',t
    call mpistop("=== strictsmall in getdt_heatconduct_hd: low pressure ===")
  end if
else
  {do ix^DB=ixmin^DB,ixmax^DB\}
    if(tmp(ix^D)<minp) then
     tmp(ix^D)=minp
    end if
  {end do\}
endif
Te(ix^S)=tmp(ix^S)/w(ix^S,rho_)
tmp(ix^S)=(eqpar(gamma_)-one)*kappa*dsqrt((Te(ix^S))**5)/w(ix^S,rho_)

do idim=1,ndim
   ! dt< dtTCpar * dx_idim**2/((gamma-1)*kappa_idim/rho)
   dtdiff_tcond=dtTCpar/maxval(tmp(ix^S)*dxinv(idim)**2)
   if(TCsaturate) then
     ! dt< dtTCpar* dx_idim**2/((gamma-1)*sqrt(Te)*5*phi)
     dtdiff_tsat=dtTCpar/maxval((eqpar(gamma_)-1.d0)*dsqrt(Te(ix^S))*&
                 5.d0*TCphi*dxinv(idim)**2)
     ! choose the slower flux (bigger time scale) between classic and saturated
     dtdiff_tcond=max(dtdiff_tcond,dtdiff_tsat)
   end if
   ! limit the time step
   dtnew=min(dtnew,dtdiff_tcond)
enddo
dtnew=dtnew/dble(ndim)

end subroutine getdt_heatconduct_hd
!=============================================================================
!> swap ghost cells e value to update temperature
subroutine getbce(time,pwuse)
use mod_thermalconduction
use mod_global_parameters

double precision, intent(in)               :: time
type(walloc), dimension(ngridshi)          :: pwuse

integer :: ixG^L, idims, iside, nwstart,nwbc
integer :: my_neighbor_type, ipole
integer :: iigrid, igrid, ineighbor, ipe_neighbor
integer :: nrecvs, nsends, isizes
integer :: ixR^L, ixS^L
integer :: ixB^L
integer :: k^L
integer :: i^D, n_i^D, ic^D, inc^D, n_inc^D
integer :: isend_buf(npwbuf), ipwbuf
type(walloc) :: pwbuf(npwbuf)
logical  :: isphysbound

double precision :: time_bcin
{#IFDEF STRETCHGRID
double precision :: logGl,qstl
}
!-----------------------------------------------------------------------------
time_bcin=MPI_WTIME()
ixG^L=ixG^LL;
nwstart=e_-1
nwbc=1
if (internalboundary) then 
   call getintbc(time,ixG^L,pwuse)
end if

! default : no singular axis
ipole=0

irecv=0
isend=0
isend_buf=0
ipwbuf=1
nrecvs=nrecv_bc_srl+nrecv_bc_r
nsends=nsend_bc_srl+nsend_bc_r

if (nrecvs>0) then
   allocate(recvstatus(MPI_STATUS_SIZE,nrecvs),recvrequest(nrecvs))
   recvrequest=MPI_REQUEST_NULL
end if
if (nsends>0) then
   allocate(sendstatus(MPI_STATUS_SIZE,nsends),sendrequest(nsends))
   sendrequest=MPI_REQUEST_NULL
end if

do iigrid=1,igridstail; igrid=igrids(iigrid);
      saveigrid=igrid
      
   {do i^DB=-1,1\}
      if (i^D==0|.and.) cycle

      my_neighbor_type=neighbor_type(i^D,igrid)
      select case (my_neighbor_type)
      case (3)
         call bc_recv_srl
      case (4)
         call bc_recv_restrict
      end select
   {end do\}
end do

do iigrid=1,igridstail; igrid=igrids(iigrid);
      saveigrid=igrid
   if (any(neighbor_type(:^D&,igrid)==2)) then
      ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);

      call coarsen_grid(pwuse(igrid)%w,px(igrid)%x,ixG^L,ixM^L,pwCoarse(igrid)%w,pxCoarse(igrid)%x, &
                        ixCoG^L,ixCoM^L,pgeo(igrid),pgeoCoarse(igrid), &
                        coarsenprimitive,.true.)
   end if
   {do i^DB=-1,1\}
      if (i^D==0|.and.) cycle

      {^IFPHI ipole=neighbor_pole(i^D,igrid)}
      my_neighbor_type=neighbor_type(i^D,igrid)
      select case (my_neighbor_type)
      case (2)
         call bc_send_restrict
      case (3)
         call bc_send_srl
      end select
   {end do\}
end do

if (irecv/=nrecvs) then
   call mpistop("number of recvs in phase1 in amr_ghostcells is incorrect")
end if
if (isend/=nsends) then
   call mpistop("number of sends in phase1 in amr_ghostcells is incorrect")
end if

if (irecv>0) then
   call MPI_WAITALL(irecv,recvrequest,recvstatus,ierrmpi)
   deallocate(recvstatus,recvrequest)
end if
if (isend>0) then
   call MPI_WAITALL(isend,sendrequest,sendstatus,ierrmpi)
   deallocate(sendstatus,sendrequest)
   do ipwbuf=1,npwbuf
      if (isend_buf(ipwbuf)/=0) deallocate(pwbuf(ipwbuf)%w)
   end do
end if


irecv=0
isend=0
isend_buf=0
ipwbuf=1
nrecvs=nrecv_bc_p
nsends=nsend_bc_p
if (nrecvs>0) then
   allocate(recvstatus(MPI_STATUS_SIZE,nrecvs),recvrequest(nrecvs))
   recvrequest=MPI_REQUEST_NULL
end if
if (nsends>0) then
   allocate(sendstatus(MPI_STATUS_SIZE,nsends),sendrequest(nsends))
   sendrequest=MPI_REQUEST_NULL
end if

do iigrid=1,igridstail; igrid=igrids(iigrid);
   saveigrid=igrid
   
   {do i^DB=-1,1\}
      if (i^D==0|.and.) cycle

      my_neighbor_type=neighbor_type(i^D,igrid)
      if (my_neighbor_type==2) call bc_recv_prolong
   {end do\}
end do
do iigrid=1,igridstail; igrid=igrids(iigrid);
   saveigrid=igrid
   ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
  if (any(neighbor_type(:^D&,igrid)==4)) then
   {do i^DB=-1,1\}
      if (i^D==0|.and.) cycle

      {^IFPHI ipole=neighbor_pole(i^D,igrid)}
      my_neighbor_type=neighbor_type(i^D,igrid)
      if (my_neighbor_type==4) call bc_send_prolong
   {end do\}
  end if
end do


if (irecv/=nrecvs) then
   call mpistop("number of recvs in phase2 in amr_ghostcells is incorrect")
end if
if (isend/=nsends) then
   call mpistop("number of sends in phase2 in amr_ghostcells is incorrect")
end if

if (irecv>0) then
   call MPI_WAITALL(irecv,recvrequest,recvstatus,ierrmpi)
   deallocate(recvstatus,recvrequest)
end if

do iigrid=1,igridstail; igrid=igrids(iigrid);
   saveigrid=igrid
   ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
   if (any(neighbor_type(:^D&,igrid)==2)) then
      {do i^DB=-1,1\}
         if (i^D==0|.and.) cycle
         my_neighbor_type=neighbor_type(i^D,igrid)
         if (my_neighbor_type==2) call bc_prolong
      {end do\}
   end if
end do

if (isend>0) then
   call MPI_WAITALL(isend,sendrequest,sendstatus,ierrmpi)
   deallocate(sendstatus,sendrequest)
   do ipwbuf=1,npwbuf
      if (isend_buf(ipwbuf)/=0) deallocate(pwbuf(ipwbuf)%w)
   end do
end if

time_bc=time_bc+(MPI_WTIME()-time_bcin)

contains
!=============================================================================
! internal procedures
!=============================================================================
subroutine bc_send_srl
!-----------------------------------------------------------------------------
ineighbor=neighbor(1,i^D,igrid)
ipe_neighbor=neighbor(2,i^D,igrid)


if (ipole==0) then
   n_i^D=-i^D;
   if (ipe_neighbor==mype) then
      ixS^L=ixS_srl_^L(i^D);
      ixR^L=ixR_srl_^L(n_i^D);
      pwuse(ineighbor)%w(ixR^S,e_)=pwuse(igrid)%w(ixS^S,e_)
   else
      isend=isend+1
      itag=(3**^ND+4**^ND)*(ineighbor-1)+{(n_i^D+1)*3**(^D-1)+}
      call MPI_ISEND(pwuse(igrid)%w,1,type_send_srl(i^D), &
                     ipe_neighbor,itag,icomm,sendrequest(isend),ierrmpi)
   end if
else
   ixS^L=ixS_srl_^L(i^D);
   select case (ipole)
   {case (^D)
      n_i^D=i^D^D%n_i^DD=-i^DD;\}
   end select
   if (ipe_neighbor==mype) then
      ixR^L=ixR_srl_^L(n_i^D);
      call pole_copy(pwuse(ineighbor),ixR^L,pwuse(igrid),ixS^L)
   else
      if (isend_buf(ipwbuf)/=0) then
         call MPI_WAIT(sendrequest(isend_buf(ipwbuf)), &
                       sendstatus(:,isend_buf(ipwbuf)),ierrmpi)
         deallocate(pwbuf(ipwbuf)%w)
      end if
      allocate(pwbuf(ipwbuf)%w(ixS^S,e_))
      call pole_copy(pwbuf(ipwbuf),ixS^L,pwuse(igrid),ixS^L)
      isend=isend+1
      isend_buf(ipwbuf)=isend
      itag=(3**^ND+4**^ND)*(ineighbor-1)+{(n_i^D+1)*3**(^D-1)+}
      isizes={(ixSmax^D-ixSmin^D+1)*}*nwbc
      call MPI_ISEND(pwbuf(ipwbuf)%w,isizes,MPI_DOUBLE_PRECISION, &
                     ipe_neighbor,itag,icomm,sendrequest(isend),ierrmpi)
      ipwbuf=1+modulo(ipwbuf,npwbuf)
   end if
end if

end subroutine bc_send_srl
!=============================================================================
subroutine bc_send_restrict
!-----------------------------------------------------------------------------
ic^D=1+modulo(node(pig^D_,igrid)-1,2);
if ({.not.(i^D==0.or.i^D==2*ic^D-3)|.or.}) return

ineighbor=neighbor(1,i^D,igrid)
ipe_neighbor=neighbor(2,i^D,igrid)

if (ipole==0) then
   n_inc^D=-2*i^D+ic^D;
   if (ipe_neighbor==mype) then
      ixS^L=ixS_r_^L(i^D);
      ixR^L=ixR_r_^L(n_inc^D);
      pwuse(ineighbor)%w(ixR^S,e_)=pwCoarse(igrid)%w(ixS^S,e_)
   else
      isend=isend+1
      itag=(3**^ND+4**^ND)*(ineighbor-1)+3**^ND+{n_inc^D*4**(^D-1)+}
      call MPI_ISEND(pwCoarse(igrid)%w,1,type_send_r(i^D), &
                     ipe_neighbor,itag,icomm,sendrequest(isend),ierrmpi)
   end if
else
   ixS^L=ixS_r_^L(i^D);
   select case (ipole)
   {case (^D)
      n_inc^D=2*i^D+(3-ic^D)^D%n_inc^DD=-2*i^DD+ic^DD;\}
   end select
   if (ipe_neighbor==mype) then
      ixR^L=ixR_r_^L(n_inc^D);
      call pole_copy(pwuse(ineighbor),ixR^L,pwCoarse(igrid),ixS^L)
   else
      if (isend_buf(ipwbuf)/=0) then
         call MPI_WAIT(sendrequest(isend_buf(ipwbuf)), &
                       sendstatus(:,isend_buf(ipwbuf)),ierrmpi)
         deallocate(pwbuf(ipwbuf)%w)
      end if
      allocate(pwbuf(ipwbuf)%w(ixS^S,e_))
      call pole_copy(pwbuf(ipwbuf),ixS^L,pwCoarse(igrid),ixS^L)
      isend=isend+1
      isend_buf(ipwbuf)=isend
      itag=(3**^ND+4**^ND)*(ineighbor-1)+3**^ND+{n_inc^D*4**(^D-1)+}
      isizes={(ixSmax^D-ixSmin^D+1)*}*nwbc
      call MPI_ISEND(pwbuf(ipwbuf)%w,isizes,MPI_DOUBLE_PRECISION, &
                     ipe_neighbor,itag,icomm,sendrequest(isend),ierrmpi)
      ipwbuf=1+modulo(ipwbuf,npwbuf)
   end if
end if

end subroutine bc_send_restrict
!=============================================================================
subroutine bc_send_prolong
integer :: ii^D
!-----------------------------------------------------------------------------
{do ic^DB=1+int((1-i^DB)/2),2-int((1+i^DB)/2)
   inc^DB=2*i^DB+ic^DB\}

   ixS^L=ixS_p_^L(inc^D);

   ineighbor=neighbor_child(1,inc^D,igrid)
   ipe_neighbor=neighbor_child(2,inc^D,igrid)

   if (ipole==0) then
      n_i^D=-i^D;
      n_inc^D=ic^D+n_i^D;
      if (ipe_neighbor==mype) then
         ixR^L=ixR_p_^L(n_inc^D);
         pwCoarse(ineighbor)%w(ixR^S,e_)=pwuse(igrid)%w(ixS^S,e_)
      else
         isend=isend+1
         itag=(3**^ND+4**^ND)*(ineighbor-1)+3**^ND+{n_inc^D*4**(^D-1)+}
         call MPI_ISEND(pwuse(igrid)%w,1,type_send_p(inc^D), &
                        ipe_neighbor,itag,icomm,sendrequest(isend),ierrmpi)
      end if
   else
      select case (ipole)
      {case (^D)
         n_inc^D=inc^D^D%n_inc^DD=ic^DD-i^DD;\}
      end select
      if (ipe_neighbor==mype) then
         ixR^L=ixR_p_^L(n_inc^D);
         call pole_copy(pwCoarse(ineighbor),ixR^L,pwuse(igrid),ixS^L)
      else
         if (isend_buf(ipwbuf)/=0) then
            call MPI_WAIT(sendrequest(isend_buf(ipwbuf)), &
                          sendstatus(:,isend_buf(ipwbuf)),ierrmpi)
            deallocate(pwbuf(ipwbuf)%w)
         end if
         allocate(pwbuf(ipwbuf)%w(ixS^S,e_))
         call pole_copy(pwbuf(ipwbuf),ixS^L,pwuse(igrid),ixS^L)
         isend=isend+1
         isend_buf(ipwbuf)=isend
         itag=(3**^ND+4**^ND)*(ineighbor-1)+3**^ND+{n_inc^D*4**(^D-1)+}
         isizes={(ixSmax^D-ixSmin^D+1)*}*nwbc
         call MPI_ISEND(pwbuf(ipwbuf)%w,isizes,MPI_DOUBLE_PRECISION, &
                        ipe_neighbor,itag,icomm,sendrequest(isend),ierrmpi)
         ipwbuf=1+modulo(ipwbuf,npwbuf)
      end if
   end if
{end do\}

end subroutine bc_send_prolong
!=============================================================================
subroutine bc_recv_srl
!-----------------------------------------------------------------------------
ipe_neighbor=neighbor(2,i^D,igrid)
if (ipe_neighbor/=mype) then
   irecv=irecv+1
   itag=(3**^ND+4**^ND)*(igrid-1)+{(i^D+1)*3**(^D-1)+}
   call MPI_IRECV(pwuse(igrid)%w,1,type_recv_srl(i^D), &
                  ipe_neighbor,itag,icomm,recvrequest(irecv),ierrmpi)
end if

end subroutine bc_recv_srl
!=============================================================================
subroutine bc_recv_restrict
!-----------------------------------------------------------------------------
{do ic^DB=1+int((1-i^DB)/2),2-int((1+i^DB)/2)
   inc^DB=2*i^DB+ic^DB\}
   ipe_neighbor=neighbor_child(2,inc^D,igrid)
   if (ipe_neighbor/=mype) then
      irecv=irecv+1
      itag=(3**^ND+4**^ND)*(igrid-1)+3**^ND+{inc^D*4**(^D-1)+}
      call MPI_IRECV(pwuse(igrid)%w,1,type_recv_r(inc^D), &
                     ipe_neighbor,itag,icomm,recvrequest(irecv),ierrmpi)
   end if
{end do\}

end subroutine bc_recv_restrict
!=============================================================================
subroutine bc_recv_prolong
!-----------------------------------------------------------------------------
ic^D=1+modulo(node(pig^D_,igrid)-1,2);
if ({.not.(i^D==0.or.i^D==2*ic^D-3)|.or.}) return

ipe_neighbor=neighbor(2,i^D,igrid)
if (ipe_neighbor/=mype) then
   irecv=irecv+1
   inc^D=ic^D+i^D;
   itag=(3**^ND+4**^ND)*(igrid-1)+3**^ND+{inc^D*4**(^D-1)+}
   call MPI_IRECV(pwCoarse(igrid)%w,1,type_recv_p(inc^D), &
                  ipe_neighbor,itag,icomm,recvrequest(irecv),ierrmpi)  
end if

end subroutine bc_recv_prolong
!=============================================================================
subroutine bc_prolong

integer :: ixFi^L,ixCo^L,ii^D
double precision :: dxFi^D, dxCo^D, xFimin^D, xComin^D, invdxCo^D
!-----------------------------------------------------------------------------
ixFi^L=ixR_srl_^L(i^D);

dxFi^D=rnode(rpdx^D_,igrid);
dxCo^D=two*dxFi^D;
invdxCo^D=1.d0/dxCo^D;

xFimin^D=rnode(rpxmin^D_,igrid)-dble(dixB)*dxFi^D;
xComin^D=rnode(rpxmin^D_,igrid)-dble(dixB)*dxCo^D;
{#IFDEF STRETCHGRID
qst=qsts(node(plevel_,igrid))
logG=logGs(node(plevel_,igrid))
qstl=qsts(node(plevel_,igrid)-1)
logGl=logGs(node(plevel_,igrid)-1)
xFimin1=rnode(rpxmin1_,igrid)*qst**(-dixB)
xComin1=rnode(rpxmin1_,igrid)*qstl**(-dixB)
}

! moved the physical boundary filling here, to only fill the
! part needed

ixComin^D=int((xFimin^D+(dble(ixFimin^D)-half)*dxFi^D-xComin^D)*invdxCo^D)+1-1;
ixComax^D=int((xFimin^D+(dble(ixFimax^D)-half)*dxFi^D-xComin^D)*invdxCo^D)+1+1;

call interpolation_linear(pwuse(igrid),ixFi^L,dxFi^D,xFimin^D, &
                           pwCoarse(igrid),dxCo^D,invdxCo^D,xComin^D)

end subroutine bc_prolong
!=============================================================================
subroutine interpolation_linear(pwFi,ixFi^L,dxFi^D,xFimin^D, &
                                pwCo,dxCo^D,invdxCo^D,xComin^D)

integer, intent(in) :: ixFi^L
double precision, intent(in) :: dxFi^D, xFimin^D,dxCo^D, invdxCo^D, xComin^D
type(walloc) :: pwCo, pwFi

integer :: ixCo^D, jxCo^D, hxCo^D, ixFi^D, ix^D, iw, idims
double precision :: xCo^D, xFi^D, eta^D
double precision :: slopeL, slopeR, slopeC, signC, signR
double precision :: slope(ndim)
!-----------------------------------------------------------------------------
{do ixFi^DB = ixFi^LIM^DB
   ! cell-centered coordinates of fine grid point
   xFi^DB=xFimin^DB+(dble(ixFi^DB)-half)*dxFi^DB

   ! indices of coarse cell which contains the fine cell
   ixCo^DB=int((xFi^DB-xComin^DB)*invdxCo^DB)+1

   ! cell-centered coordinate for coarse cell
   xCo^DB=xComin^DB+(dble(ixCo^DB)-half)*dxCo^DB\}
{#IFDEF STRETCHGRID
   xFi1=xFimin1/(one-half*logG)*qst**(ixFi1-1)
   do ixCo1=1,ixCoGmax1
     xCo1=xComin1/(one-half*logGl)*qstl**(ixCo1-1)
     if(dabs(xFi1-xCo1)<half*logGl*xCo1) exit
   end do
}
   ! normalized distance between fine/coarse cell center
   ! in coarse cell: ranges from -0.5 to 0.5 in each direction
   ! (origin is coarse cell center)
   if (slab) then
      eta^D=(xFi^D-xCo^D)*invdxCo^D;
   else
      ix^D=2*int((ixFi^D+ixMlo^D)/2)-ixMlo^D;
      {eta^D=(xFi^D-xCo^D)*invdxCo^D &
            *two*(one-pgeo(igrid)%dvolume(ixFi^DD) &
            /sum(pgeo(igrid)%dvolume(ix^D:ix^D+1^D%ixFi^DD))) \}
{#IFDEF STRETCHGRID
      eta1=(xFi1-xCo1)/(logGl*xCo1)*two*(one-pgeo(igrid)%dvolume(ixFi^D) &
            /sum(pgeo(igrid)%dvolume(ix1:ix1+1^%1ixFi^D))) 
}
   end if
   iw=e_
   do idims=1,ndim
      hxCo^D=ixCo^D-kr(^D,idims)\
      jxCo^D=ixCo^D+kr(^D,idims)\

      slopeL=pwCo%w(ixCo^D,iw)-pwCo%w(hxCo^D,iw)
      slopeR=pwCo%w(jxCo^D,iw)-pwCo%w(ixCo^D,iw)
      slopeC=half*(slopeR+slopeL)

      ! get limited slope
      signR=sign(one,slopeR)
      signC=sign(one,slopeC)
      select case(typeprolonglimit)
      case('minmod')
        slope(idims)=signR*max(zero,min(dabs(slopeR), &
                                          signR*slopeL))
      case('woodward')
        slope(idims)=two*signR*max(zero,min(dabs(slopeR), &
                           signR*slopeL,signR*half*slopeC))
      case('mcbeta')
        slope(idims)=signR*max(zero,min(mcbeta*dabs(slopeR), &
                           mcbeta*signR*slopeL,signR*slopeC))
      case('koren')
        slope(idims)=signR*max(zero,min(two*signR*slopeL, &
         (dabs(slopeR)+two*slopeL*signR)*third,two*dabs(slopeR)))
      case default
        slope(idims)=signC*max(zero,min(dabs(slopeC), &
                          signC*slopeL,signC*slopeR))
      end select
   end do

   ! Interpolate from coarse cell using limited slopes
   pwFi%w(ixFi^D,e_)=pwCo%w(ixCo^D,e_)+{(slope(^D)*eta^D)+}

{end do\}

end subroutine interpolation_linear
!=============================================================================
subroutine pole_copy(pwrecv,ixR^L,pwsend,ixS^L)

integer, intent(in) :: ixR^L, ixS^L
type(walloc) :: pwrecv, pwsend

integer :: iw, iB
!-----------------------------------------------------------------------------
select case (ipole)
{case (^D)
   iside=int((i^D+3)/2)
   iB=2*(^D-1)+iside
   do iw=nwstart+1,nwstart+nwbc
      select case (typeB(iw,iB))
      case ("symm")
         pwrecv%w(ixR^S,iw) = pwsend%w(ixSmax^D:ixSmin^D:-1^D%ixS^S,iw)
      case ("asymm")
         pwrecv%w(ixR^S,iw) =-pwsend%w(ixSmax^D:ixSmin^D:-1^D%ixS^S,iw)
      case default
         call mpistop("Boundary condition at pole should be symm or asymm")
      end select
   end do \}
end select

end subroutine pole_copy
!=============================================================================
end subroutine getbce
