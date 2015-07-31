!=============================================================================
subroutine set_pole

include 'amrvacdef.f'
!-----------------------------------------------------------------------------
select case (typeaxial)
case ("spherical") 
case ("cylindrical")
  if (1==-1.and.periodB(1)) then
    if(mod(ng1(1),2)/=0) &
      call mpistop("Number of meshes in phi-direction should be even!")
    if(abs(xprobmin1)<smalldouble) then
      if(mype==0) write(unitterm,*) "Will apply pi-periodic conditions at r=0"
      poleB(1,1)=.true.
    else
      if(mype==0) write(unitterm,*) "There is no cylindrical axis!"
    end if
  end if
end select

end subroutine set_pole
!=============================================================================
subroutine getgridgeo(igrid)

include 'amrvacdef.f'

integer, intent(in) :: igrid

integer :: ixmin1,ixmax1, ixCoGmin1,ixCoGmax1, ixCoMmin1,ixCoMmax1, ixComin1,&
   ixComax1, ixCoCoGmin1,ixCoCoGmax1, ixGextmin1,ixGextmax1
double precision :: xmin1, dx1
!-----------------------------------------------------------------------------
!ix^L=ixM^LL^LADD1;
ixmin1=ixGlo1+1;ixmax1=ixGhi1-1;
if (2*int(dixB/2)==dixB) then
   ixGextmin1=ixGlo1;ixGextmax1=ixGhi1;
else
   ixGextmin1=ixGlo1-1;ixGextmax1=ixGhi1+1;
end if


allocate(pgeo(igrid)%surfaceC1(ixmin1-1:ixmax1), pgeo(igrid)%surface1(ixmin1&
   -1:ixmax1), pgeo(igrid)%dvolume(ixGextmin1:ixGextmax1), &
   pgeo(igrid)%dx(ixGextmin1:ixGextmax1,1:ndim))

dx1=rnode(rpdx1_,igrid);
xmin1=rnode(rpxmin1_,igrid);
call fillgeo(pgeo(igrid),ixGlo1,ixGhi1,ixGextmin1,ixGextmax1,xmin1,dx1,&
   .false.)

if (errorestimate==1) then
   ixCoGmin1=1; ixCoGmax1=ixGhi1/2+dixB;
   if (2*int(dixB/2)==dixB) then
      ixGextmin1=ixCoGmin1;ixGextmax1=ixCoGmax1;
   else
      ixGextmin1=ixCoGmin1-1;ixGextmax1=ixCoGmax1+1;
   end if
   ixCoMmin1=ixCoGmin1+dixB;ixCoMmax1=ixCoGmax1-dixB;
   ixComin1=ixCoMmin1-1;ixComax1=ixCoMmax1+1;
   ixComin1=ixCoGmin1+1;ixComax1=ixCoGmax1-1;

   allocate(pgeoCoarse(igrid)%surfaceC1(ixComin1-1:ixComax1),&
       pgeoCoarse(igrid)%surface1(ixComin1-1:ixComax1), pgeoCoarse&
      (igrid)%dvolume(ixGextmin1:ixGextmax1), pgeoCoarse(igrid)%dx&
      (ixGextmin1:ixGextmax1,1:ndim))

   dx1=two*rnode(rpdx1_,igrid);
   call fillgeo(pgeoCoarse(igrid),ixCoGmin1,ixCoGmax1,ixGextmin1,ixGextmax1,&
      xmin1,dx1,.false.)

   ixCoCoGmin1=1; ixCoCoGmax1=ixCoGmax1/2+dixB;
   if (2*int(dixB/2)==dixB) then
      ixGextmin1=ixCoCoGmin1;ixGextmax1=ixCoCoGmax1;
   else
      ixGextmin1=ixCoCoGmin1-1;ixGextmax1=ixCoCoGmax1+1;
   end if


   allocate(pgeoCoCo(igrid)%dvolume(ixGextmin1:ixGextmax1))

   dx1=4.0d0*rnode(rpdx1_,igrid);
   call fillgeo(pgeoCoCo(igrid),ixCoCoGmin1,ixCoCoGmax1,ixGextmin1,ixGextmax1,&
      xmin1,dx1,.true.)
else
   ixCoGmin1=1; ixCoGmax1=ixGhi1/2+dixB;
   if (2*int(dixB/2)==dixB) then
      ixGextmin1=ixCoGmin1;ixGextmax1=ixCoGmax1;
   else
      ixGextmin1=ixCoGmin1-1;ixGextmax1=ixCoGmax1+1;
   end if


   allocate(pgeoCoarse(igrid)%dvolume(ixGextmin1:ixGextmax1))

   dx1=two*rnode(rpdx1_,igrid);
   call fillgeo(pgeoCoarse(igrid),ixCoGmin1,ixCoGmax1,ixGextmin1,ixGextmax1,&
      xmin1,dx1,.true.)
end if

end subroutine getgridgeo
!=============================================================================
subroutine putgridgeo(igrid)

  include 'amrvacdef.f'

integer, intent(in) :: igrid
!-----------------------------------------------------------------------------
if (errorestimate==1) then
   deallocate(pgeo(igrid)%surfaceC1,pgeo(igrid)%surface1,pgeo(igrid)%dvolume,&
      pgeo(igrid)%dx, pgeoCoarse(igrid)%surfaceC1,pgeoCoarse(igrid)%surface1,&
      pgeoCoarse(igrid)%dvolume,pgeoCoarse(igrid)%dx,pgeoCoCo(igrid)%dvolume)
   
else
   deallocate(pgeo(igrid)%surfaceC1,pgeo(igrid)%surface1,pgeo(igrid)%dvolume,&
      pgeo(igrid)%dx,pgeoCoarse(igrid)%dvolume)
end if

end subroutine putgridgeo
!=============================================================================
subroutine fillgeo(pgeogrid,ixGmin1,ixGmax1,ixGextmin1,ixGextmax1,xmin1,dx1,&
   need_only_volume)

include 'amrvacdef.f'

type(geoalloc) :: pgeogrid
integer, intent(in) :: ixGmin1,ixGmax1, ixGextmin1,ixGextmax1
double precision, intent(in) :: xmin1, dx1
logical, intent(in) :: need_only_volume

integer :: idims, ix, ixMmin1,ixMmax1, ixmin1,ixmax1, ixCmin1,ixCmax1
double precision :: x(ixGextmin1:ixGextmax1,ndim)
!-----------------------------------------------------------------------------
ixMmin1=ixGmin1+dixB;ixMmax1=ixGmax1-dixB;
!ix^L=ixM^L^LADD1;
ixmin1=ixGmin1+1;ixmax1=ixGmax1-1;

select case (typeaxial)
case ("slabtest")

   pgeogrid%dvolume(ixGextmin1:ixGextmax1) = dx1

   if (need_only_volume) return

   ixCmin1=ixmin1-kr(1,1); ixCmax1=ixmax1;
   pgeogrid%surfaceC1(ixCmin1:ixCmax1)= one
   pgeogrid%surface1(ixCmin1:ixCmax1) = one
   
   
   

   pgeogrid%dx(ixGextmin1:ixGextmax1,1)=dx1;


case ("spherical")

   do idims=1,min(ndim,2)
      select case(idims)
      case(1)
         do ix = ixGextmin1,ixGextmax1
            x(ix,1)=xmin1+(dble(ix-dixB)-half)*dx1
         end do
      end select
   end do

   if(typespherical==0) then
     pgeogrid%dvolume(ixGextmin1:ixGextmax1)=(x(ixGextmin1:ixGextmax1,1)**2&
        +dx1**2/12.0d0)*dx1 
   else
     pgeogrid%dvolume(ixGextmin1:ixGextmax1)=(x(ixGextmin1:ixGextmax1,1)&
        **2)*dx1 
   endif

   if (need_only_volume) return


   ixCmin1=ixmin1-kr(1,1); ixCmax1=ixmax1;
   if(typespherical==0) then
       pgeogrid%surfaceC1(ixCmin1:ixCmax1)=(x(ixCmin1:ixCmax1,1)+half*dx1)**2 
   else
       pgeogrid%surfaceC1(ixCmin1:ixCmax1)=(x(ixCmin1:ixCmax1,1)+half*dx1)**2 
   endif

   

   

   ixCmin1=ixmin1-kr(1,1); ixCmax1=ixmax1;
   if(typespherical==0) then
       pgeogrid%surface1(ixCmin1:ixCmax1)=x(ixCmin1:ixCmax1,1)**2 
   else
      pgeogrid%surface1(ixCmin1:ixCmax1)=x(ixCmin1:ixCmax1,1)**2 
   endif

   

   

   pgeogrid%dx(ixGextmin1:ixGextmax1,1)=dx1
   
   

case ("cylindrical")

   do ix = ixGextmin1,ixGextmax1
      x(ix,1)=xmin1+(dble(ix-dixB)-half)*dx1
   end do

   pgeogrid%dvolume(ixGextmin1:ixGextmax1)=dabs(x(ixGextmin1:ixGextmax1,1))&
      *dx1
   pgeogrid%dvolume(ixGextmin1:ixGextmax1)=dabs(half*((x&
      (ixGextmin1:ixGextmax1,1)+half*dx1)**2-(x(ixGextmin1:ixGextmax1,1)&
      -half*dx1)**2))

   if (need_only_volume) return


   ixCmin1=ixmin1-kr(1,1); ixCmax1=ixmax1;
   !!pgeogrid%surfaceC1(ixC^S)=(x(ixC^S,1)+half*dx1){^DE&*dx^DE }
   pgeogrid%surfaceC1(ixCmin1:ixCmax1)=dabs((x(ixCmin1:ixCmax1,1)+half*dx1))
   
   

   ixCmin1=ixmin1-kr(1,1); ixCmax1=ixmax1;
   !!pgeogrid%surface1(ixC^S)=x(ixC^S,1){^DE&*dx^DE }
   pgeogrid%surface1(ixCmin1:ixCmax1)=dabs(x(ixCmin1:ixCmax1,1))
   
   


   pgeogrid%dx(ixGextmin1:ixGextmax1,1)=dx1
   
   

case default

   call mpistop("Sorry, typeaxial unknown")
   
end select

end subroutine fillgeo
!=============================================================================
subroutine gradient(q,ixmin1,ixmax1,idir,gradq)

! Calculate gradient of a scalar q within ixL in direction idir

include 'amrvacdef.f'

integer :: ixmin1,ixmax1, idir
double precision :: q(ixGlo1:ixGhi1), gradq(ixGlo1:ixGhi1)

double precision :: qC(ixGlo1:ixGhi1),invdx
integer :: jxmin1,jxmax1, hxmin1,hxmax1, ixCmin1,ixCmax1, jxCmin1,jxCmax1 

!-----------------------------------------------------------------------------

invdx=1.d0/dxlevel(idir)
if (slab) then

   jxmin1=ixmin1+kr(idir,1);jxmax1=ixmax1+kr(idir,1);
   hxmin1=ixmin1-kr(idir,1);hxmax1=ixmax1-kr(idir,1);
   gradq(ixmin1:ixmax1) = half*(q(jxmin1:jxmax1)-q(hxmin1:hxmax1))*invdx

else
   hxmin1=ixmin1-kr(idir,1);hxmax1=ixmax1-kr(idir,1);
   ixCmin1=hxmin1;ixCmax1=ixmax1;
   jxCmin1=ixCmin1+kr(idir,1);jxCmax1=ixCmax1+kr(idir,1);
   select case(idir)
   case(1)
      qC(ixCmin1:ixCmax1)=mygeo%surfaceC1(ixCmin1:ixCmax1)*half&
         *(q(ixCmin1:ixCmax1)+q(jxCmin1:jxCmax1))
      gradq(ixmin1:ixmax1)=(qC(ixmin1:ixmax1)-qC(hxmin1:hxmax1))&
         /mygeo%dvolume(ixmin1:ixmax1)
      ! Substract difference divergence and gradient
      gradq(ixmin1:ixmax1)=gradq(ixmin1:ixmax1)-q(ixmin1:ixmax1) &
                     *(mygeo%surfaceC1(ixmin1:ixmax1)-mygeo%surfaceC1&
                        (hxmin1:hxmax1)) &
                    /mygeo%dvolume(ixmin1:ixmax1) 
   end select
end if

end subroutine gradient
!=============================================================================
subroutine gradientS(q,ixmin1,ixmax1,idir,gradq)

! Calculate gradient of a scalar q within ixL in direction idir
! first use limiter to go from cell center to edge

include 'amrvacdef.f'

integer :: ixmin1,ixmax1, idir
double precision :: q(ixGlo1:ixGhi1), gradq(ixGlo1:ixGhi1)
double precision :: dxdim

double precision :: qC(ixGlo1:ixGhi1)
double precision,dimension(ixGlo1:ixGhi1):: qL,qR,dqC,ldq,invdx
integer                          :: hxmin1,hxmax1,ixCmin1,ixCmax1,jxCmin1,&
   jxCmax1,gxCmin1,gxCmax1,hxCmin1,hxCmax1,idummy
character*79 :: savetypelimiter,savetypegradlimiter,save2typelimiter
!-----------------------------------------------------------------------------

invdx=1.d0/dxlevel(idir)
hxmin1=ixmin1-kr(idir,1);hxmax1=ixmax1-kr(idir,1);
ixCmin1=hxmin1;ixCmax1=ixmax1;
jxCmin1=ixCmin1+kr(idir,1);jxCmax1=ixCmax1+kr(idir,1);
gxCmin1=ixCmin1-kr(idir,1);gxCmax1=jxCmax1;
hxCmin1=gxCmin1+kr(idir,1);hxCmax1=gxCmax1+kr(idir,1);
idummy=0

savetypelimiter=typelimiter
savetypegradlimiter=typegradlimiter
! set the gradient limiter here
typelimiter=typegradlimiter
qR(gxCmin1:gxCmax1) = q(hxCmin1:hxCmax1)
qL(gxCmin1:gxCmax1) = q(gxCmin1:gxCmax1)
if (typelimiter/='ppm') then
   dqC(gxCmin1:gxCmax1)= qR(gxCmin1:gxCmax1)-qL(gxCmin1:gxCmax1)
   save2typelimiter=typelimiter
   if(save2typelimiter=='koren') typelimiter='korenL'
   if(save2typelimiter=='cada')  typelimiter='cadaL'
   if(save2typelimiter=='cada3') typelimiter='cada3L'
   dxdim=dxlevel(idir)
   call dwlimiter2(dqC,gxCmin1,gxCmax1,idummy,idir,ldq,dxdim)
   qL(ixCmin1:ixCmax1) = qL(ixCmin1:ixCmax1) + half*ldq(ixCmin1:ixCmax1)
   if(save2typelimiter=='koren')then
     typelimiter='korenR'
     call dwlimiter2(dqC,gxCmin1,gxCmax1,idummy,idir,ldq,dxdim)
   endif
   if(save2typelimiter=='cada')then
     typelimiter='cadaR'
     call dwlimiter2(dqC,gxCmin1,gxCmax1,idummy,idir,ldq,dxdim)
   endif
   if(save2typelimiter=='cada3')then
     typelimiter='cada3R'
     call dwlimiter2(dqC,gxCmin1,gxCmax1,idummy,idir,ldq,dxdim)
   endif
   typelimiter=save2typelimiter
   qR(ixCmin1:ixCmax1) = qR(ixCmin1:ixCmax1) - half*ldq(jxCmin1:jxCmax1)
else
   call PPMlimitervar(ixGlo1,ixGhi1,ixMlo1,ixMhi1,idir,q,q,qL,qR)
endif
! set the method limiter back
typelimiter=savetypelimiter
typegradlimiter=savetypegradlimiter

if (slab) then
   gradq(ixmin1:ixmax1)=half*(qR(ixmin1:ixmax1)-qL(hxmin1:hxmax1))*invdx
else
   select case(idir)
   case(1)
    gradq(ixmin1:ixmax1)=(qR(ixmin1:ixmax1)-qL(hxmin1:hxmax1))&
       /mygeo%dx(ixmin1:ixmax1,idir) 
   end select
end if

end subroutine gradientS
!=============================================================================
subroutine divvector(qvec,ixImin1,ixImax1,ixOmin1,ixOmax1,divq)

! Calculate divergence of a vector qvec within ixL

include 'amrvacdef.f'

integer :: ixImin1,ixImax1,ixOmin1,ixOmax1
double precision :: qvec(ixGlo1:ixGhi1,1:ndir), divq(ixGlo1:ixGhi1)

double precision :: qC(ixGlo1:ixGhi1), invdx(1:ndim)
integer :: jxOmin1,jxOmax1, hxOmin1,hxOmax1, ixCmin1,ixCmax1, jxCmin1,jxCmax1,&
    idims, ixmin1,ixmax1 
!-----------------------------------------------------------------------------

ixmin1=ixOmin1-1;ixmax1=ixOmax1+1;

if (ixImin1>ixmin1.or.ixImax1<ixmax1) call mpistop&
   ("Error in divvector: Non-conforming input limits")
invdx=1.d0/dxlevel
divq(ixOmin1:ixOmax1)=zero
if (slab) then
  do idims=1,ndim

     jxOmin1=ixOmin1+kr(idims,1);jxOmax1=ixOmax1+kr(idims,1);
     hxOmin1=ixOmin1-kr(idims,1);hxOmax1=ixOmax1-kr(idims,1);
     divq(ixOmin1:ixOmax1)=divq(ixOmin1:ixOmax1)+half*(qvec(jxOmin1:jxOmax1,&
        idims)-qvec(hxOmin1:hxOmax1,idims))*invdx(idims)

  end do
else
  do idims=1,ndim
     hxOmin1=ixOmin1-kr(idims,1);hxOmax1=ixOmax1-kr(idims,1);
     ixCmin1=hxOmin1;ixCmax1=ixOmax1;
     jxCmin1=ixCmin1+kr(idims,1);jxCmax1=ixCmax1+kr(idims,1);
     select case(idims)
     case(1)
        qC(ixCmin1:ixCmax1)=mygeo%surfaceC1(ixCmin1:ixCmax1)*half&
           *(qvec(ixCmin1:ixCmax1,idims)+qvec(jxCmin1:jxCmax1,idims))
        divq(ixOmin1:ixOmax1)=divq(ixOmin1:ixOmax1)+qC(ixOmin1:ixOmax1)&
           -qC(hxOmin1:hxOmax1) 
      end select
  end do
  divq(ixOmin1:ixOmax1)=divq(ixOmin1:ixOmax1)/mygeo%dvolume(ixOmin1:ixOmax1)
end if


end subroutine divvector 
!=============================================================================
subroutine curlvector(qvec,ixImin1,ixImax1,ixOmin1,ixOmax1,curlvec,idirmin,&
   idirmin0,ndir0)

! Calculate curl of a vector qvec within ixL

include 'amrvacdef.f'

integer :: ixImin1,ixImax1,ixOmin1,ixOmax1,idirmin,ixmin1,ixmax1,idir,jdir,&
   kdir,hxOmin1,hxOmax1,jxOmin1,jxOmax1,ndir0,idirmin0
double precision :: qvec(ixGlo1:ixGhi1,1:ndir0),curlvec(ixGlo1:ixGhi1,&
   idirmin0:3), invdx(1:ndim)
double precision :: tmp(ixGlo1:ixGhi1),tmp2(ixGlo1:ixGhi1),&
   surface(ixGlo1:ixGhi1),mydx(ixGlo1:ixGhi1)
!-----------------------------------------------------------------------------

ixmin1=ixOmin1-1;ixmax1=ixOmax1+1;

if (ixImin1>ixmin1.or.ixImax1<ixmax1) call mpistop&
   ("Error in curl: Non-conforming input limits")

! Calculate curl within ixL: CurlV_i=eps_ijk*d_j V_k
! Curl can have components (idirmin0:3)
! Determine exact value of idirmin while doing the loop.

invdx=1.d0/dxlevel
idirmin=4
curlvec(ixOmin1:ixOmax1,idirmin0:3)=zero

do idir=idirmin0,3; do jdir=1,ndim; do kdir=1,ndir0
   if(lvc(idir,jdir,kdir)/=0)then
      tmp(ixmin1:ixmax1)=qvec(ixmin1:ixmax1,kdir)
      hxOmin1=ixOmin1-kr(jdir,1);hxOmax1=ixOmax1-kr(jdir,1);
      jxOmin1=ixOmin1+kr(jdir,1);jxOmax1=ixOmax1+kr(jdir,1);
      if(slab)then


         tmp2(ixOmin1:ixOmax1)=half*(tmp(jxOmin1:jxOmax1)-tmp&
            (hxOmin1:hxOmax1))*invdx(jdir)

      else
         ! approximate formula, reduces to slab case
         ! and avoids staggering

         if (kdir .le. ndim) then 
            mydx(ixmin1:ixmax1)=mygeo%dx(ixmin1:ixmax1,kdir)
         else 
            mydx(ixmin1:ixmax1)=one
         end if

         select case(idir)
           case(1)
             surface(ixOmin1:ixOmax1)=mygeo%surface1(ixOmin1:ixOmax1)
             tmp2(ixOmin1:ixOmax1)=half*(mydx(jxOmin1:jxOmax1)&
                *tmp(jxOmin1:jxOmax1) &
                              -mydx(hxOmin1:hxOmax1)*tmp(hxOmin1:hxOmax1)) &
                     /surface(ixOmin1:ixOmax1) 
          end select
      endif
      if(lvc(idir,jdir,kdir)==1)then
         curlvec(ixOmin1:ixOmax1,idir)=curlvec(ixOmin1:ixOmax1,idir)&
            +tmp2(ixOmin1:ixOmax1)
      else
         curlvec(ixOmin1:ixOmax1,idir)=curlvec(ixOmin1:ixOmax1,idir)&
            -tmp2(ixOmin1:ixOmax1)
      endif
      if(idir<idirmin)idirmin=idir
   endif
enddo; enddo; enddo;

end subroutine curlvector 
!=============================================================================
subroutine divvectorS(qvec,ixImin1,ixImax1,ixOmin1,ixOmax1,divq)

! Calculate divergence of a vector qvec within ixL
! using limited extrapolation to cell edges

include 'amrvacdef.f'

integer :: ixImin1,ixImax1,ixOmin1,ixOmax1
double precision :: qvec(ixGlo1:ixGhi1,1:ndir), divq(ixGlo1:ixGhi1)

double precision,dimension(ixGlo1:ixGhi1):: qL,qR,dqC,ldq
double precision :: dxdim, invdx(1:ndim)

integer :: hxOmin1,hxOmax1,ixCmin1,ixCmax1,jxCmin1,jxCmax1,idims,ixmin1,&
   ixmax1,gxCmin1,gxCmax1,hxCmin1,hxCmax1,idummy
character*79, save :: savetypelimiter,savetypegradlimiter,save2typelimiter
!-----------------------------------------------------------------------------
ixmin1=ixOmin1-2;ixmax1=ixOmax1+2;

if (ixImin1>ixmin1.or.ixImax1<ixmax1) call mpistop&
   ("Error in divvectorS: Non-conforming input limits")

idummy=0
invdx=1.d0/dxlevel
divq(ixOmin1:ixOmax1)=zero
do idims=1,ndim
   hxOmin1=ixOmin1-kr(idims,1);hxOmax1=ixOmax1-kr(idims,1);
   ixCmin1=hxOmin1;ixCmax1=ixOmax1;
   jxCmin1=ixCmin1+kr(idims,1);jxCmax1=ixCmax1+kr(idims,1);
   gxCmin1=ixCmin1-kr(idims,1);gxCmax1=jxCmax1;
   hxCmin1=gxCmin1+kr(idims,1);hxCmax1=gxCmax1+kr(idims,1);
   savetypelimiter=typelimiter
   savetypegradlimiter=typegradlimiter
   ! set the gradient limiter here
   typelimiter=savetypegradlimiter
  !! {if(gxCmin^D<ixGlo1)then
  !!   gxCmin^D=ixGlo1
  !!   hxCmin^D=hxCmin^D+1
  !! endif \}
  !! {if(hxCmax^D>ixGhi1)then
  !!   hxCmax^D=ixGhi1
  !!   gxCmax^D=gxCmax^D-1
  !! endif \}
   qR(gxCmin1:gxCmax1) = qvec(hxCmin1:hxCmax1,idims)
   qL(gxCmin1:gxCmax1) = qvec(gxCmin1:gxCmax1,idims)
   if(typelimiter/='ppm') then
      dqC(gxCmin1:gxCmax1)= qR(gxCmin1:gxCmax1)-qL(gxCmin1:gxCmax1)
      save2typelimiter=typelimiter
      if(save2typelimiter=='koren') typelimiter='korenL'
      if(save2typelimiter=='cada')  typelimiter='cadaL'
      if(save2typelimiter=='cada3') typelimiter='cada3L'
      dxdim=dxlevel(idims)
      call dwlimiter2(dqC,gxCmin1,gxCmax1,idummy,idims,ldq,dxdim)
      qL(ixCmin1:ixCmax1) = qL(ixCmin1:ixCmax1) + half*ldq(ixCmin1:ixCmax1)
      if(save2typelimiter=='koren')then
         typelimiter='korenR'
         call dwlimiter2(dqC,gxCmin1,gxCmax1,idummy,idims,ldq,dxdim)
       endif
      if(save2typelimiter=='cada')then
         typelimiter='cadaR'
         call dwlimiter2(dqC,gxCmin1,gxCmax1,idummy,idims,ldq,dxdim)
       endif
      if(save2typelimiter=='cada3')then
         typelimiter='cada3R'
         call dwlimiter2(dqC,gxCmin1,gxCmax1,idummy,idims,ldq,dxdim)
       endif
       typelimiter=save2typelimiter
      qR(ixCmin1:ixCmax1) = qR(ixCmin1:ixCmax1) - half*ldq(jxCmin1:jxCmax1)
   else if (typelimiter .eq. 'ppm') then
      dqC(ixImin1:ixImax1)=qvec(ixImin1:ixImax1,idims)
      call PPMlimitervar(ixGlo1,ixGhi1,ixMlo1,ixMhi1,idims,dqC,dqC,qL,qR)
   else
      call mpistop('typelimiter unknown in divvectorS')
   endif
   ! set the method limiter back
   typegradlimiter=typelimiter
   typelimiter=savetypelimiter

   if (slab) then
     divq(ixOmin1:ixOmax1)=divq(ixOmin1:ixOmax1)+half*(qR(ixOmin1:ixOmax1)&
        -qL(hxOmin1:hxOmax1))*invdx(idims)
   else
     select case(idims)
     case(1)
        qR(ixCmin1:ixCmax1)=mygeo%surfaceC1(ixCmin1:ixCmax1)&
           *qR(ixCmin1:ixCmax1)
        qL(ixCmin1:ixCmax1)=mygeo%surfaceC1(ixCmin1:ixCmax1)&
           *qL(ixCmin1:ixCmax1)
        divq(ixOmin1:ixOmax1)=divq(ixOmin1:ixOmax1)+qR(ixOmin1:ixOmax1)&
           -qL(hxOmin1:hxOmax1) 
      end select
   end if
end do
if(.not.slab) divq(ixOmin1:ixOmax1)=divq(ixOmin1:ixOmax1)/mygeo%dvolume&
   (ixOmin1:ixOmax1)

end subroutine divvectorS
!=============================================================================
subroutine extremaq(ixImin1,ixImax1,ixOmin1,ixOmax1,q,nshift,qMax,qMin)

include 'amrvacdef.f'

integer,intent(in)           :: ixImin1,ixImax1,ixOmin1,ixOmax1
double precision, intent(in) :: q(ixImin1:ixImax1)
integer,intent(in)           :: nshift

double precision, intent(out) :: qMax(ixGlo1:ixGhi1),qMin(ixGlo1:ixGhi1)

integer           :: ixsmin1,ixsmax1,ixsRmin1,ixsRmax1,ixsLmin1,ixsLmax1,&
   idims,jdims,kdims,ishift,i,j 
!-------------------------------------------------------------------------
do ishift=1,nshift
 idims=1
 ixsRmin1=ixOmin1+ishift*kr(idims,1);ixsRmax1=ixOmax1+ishift*kr(idims,1);
 ixsLmin1=ixOmin1-ishift*kr(idims,1);ixsLmax1=ixOmax1-ishift*kr(idims,1);
 if (ishift==1) then
   qMax(ixOmin1:ixOmax1)=max(q(ixOmin1:ixOmax1),q(ixsRmin1:ixsRmax1),&
      q(ixsLmin1:ixsLmax1))
   qMin(ixOmin1:ixOmax1)=min(q(ixOmin1:ixOmax1),q(ixsRmin1:ixsRmax1),&
      q(ixsLmin1:ixsLmax1))
 else
   qMax(ixOmin1:ixOmax1)=max(qMax(ixOmin1:ixOmax1),q(ixsRmin1:ixsRmax1),&
      q(ixsLmin1:ixsLmax1))
   qMin(ixOmin1:ixOmax1)=min(qMin(ixOmin1:ixOmax1),q(ixsRmin1:ixsRmax1),&
      q(ixsLmin1:ixsLmax1))
 end if
 
 
enddo

end subroutine  extremaq
!=============================================================================
subroutine extremaw(ixImin1,ixImax1,ixOmin1,ixOmax1,w,nshift,wMax,wMin)

include 'amrvacdef.f'

integer,intent(in)            :: ixImin1,ixImax1,ixOmin1,ixOmax1
double precision, intent(in)  :: w(ixImin1:ixImax1,1:nw)
integer,intent(in)            :: nshift

double precision, intent(out) :: wMax(ixGlo1:ixGhi1,1:nwflux),&
   wMin(ixGlo1:ixGhi1,1:nwflux)

integer          :: ixsmin1,ixsmax1,ixsRmin1,ixsRmax1,ixsLmin1,ixsLmax1,idims,&
   jdims,kdims,ishift,i,j
!-------------------------------------------------------------------------
do ishift=1,nshift
 idims=1
 ixsRmin1=ixOmin1+ishift*kr(idims,1);ixsRmax1=ixOmax1+ishift*kr(idims,1);
 ixsLmin1=ixOmin1-ishift*kr(idims,1);ixsLmax1=ixOmax1-ishift*kr(idims,1);
 if (ishift==1) then
    wMax(ixOmin1:ixOmax1,1:nwflux)= max(w(ixOmin1:ixOmax1,1:nwflux),&
       w(ixsRmin1:ixsRmax1,1:nwflux),w(ixsLmin1:ixsLmax1,1:nwflux))
    wMin(ixOmin1:ixOmax1,1:nwflux)= min(w(ixOmin1:ixOmax1,1:nwflux),&
       w(ixsRmin1:ixsRmax1,1:nwflux),w(ixsLmin1:ixsLmax1,1:nwflux))
 else
    wMax(ixOmin1:ixOmax1,1:nwflux)= max(wMax(ixOmin1:ixOmax1,1:nwflux),&
       w(ixsRmin1:ixsRmax1,1:nwflux),w(ixsLmin1:ixsLmax1,1:nwflux))
    wMin(ixOmin1:ixOmax1,1:nwflux)= min(wMin(ixOmin1:ixOmax1,1:nwflux),&
       w(ixsRmin1:ixsRmax1,1:nwflux),w(ixsLmin1:ixsLmax1,1:nwflux))
 end if
 
 
enddo

end subroutine  extremaw
!=============================================================================
subroutine extremaa(ixImin1,ixImax1,ixOmin1,ixOmax1,a,nshift,aMin)

include 'amrvacdef.f'

integer,intent(in)           :: ixImin1,ixImax1,ixOmin1,ixOmax1
double precision, intent(in) :: a(ixGlo1:ixGhi1)
integer,intent(in)           :: nshift

double precision, intent(out) :: aMin(ixGlo1:ixGhi1)

integer          :: ixsmin1,ixsmax1,ixsRmin1,ixsRmax1,ixsLmin1,ixsLmax1,idims,&
   jdims,kdims,ishift,i,j
!-------------------------------------------------------------------------
do ishift=1,nshift
  idims=1
  ixsRmin1=ixOmin1+ishift*kr(idims,1);ixsRmax1=ixOmax1+ishift*kr(idims,1);
  ixsLmin1=ixOmin1-ishift*kr(idims,1);ixsLmax1=ixOmax1-ishift*kr(idims,1);
  aMin(ixOmin1:ixOmax1)=min(a(ixsRmin1:ixsRmax1),a(ixOmin1:ixOmax1),&
     a(ixsLmin1:ixsLmax1))
  
  
end do

end subroutine extremaa
!=============================================================================
subroutine locate_in_table(xpoint,table,nxtp,ixtp,res)
!  Fast search table to find xpoint's location in the table
!  record the distance between xpoint and the its cloest element table(ixtp)
! INPUT:
!  xpoint : the point you want to locate in the table
!  table : 1D table you want to find location in
!  nxtp : number of elements in the table
! OUTPUT:
!  ixtp : closest element's index in table
!  res : offset(distance) from the closest element in table
include 'amrvacdef.f'

integer, intent(in) :: nxtp
double precision,intent(in)   :: xpoint,table(nxtp)
double precision, intent(out) :: res
integer, intent(out) :: ixtp

integer :: jl,jc,jh
!-----------------------------------------------------------------------------
if(xpoint < table(1)) then
  ixtp=1
  res=(xpoint-table(1))/(table(2)-table(1))
else if (xpoint > table(nxtp)) then
  ixtp=nxtp
  res=(xpoint-table(nxtp))/(table(nxtp)-table(nxtp-1))
else
  jl=0
  jh=nxtp+1
  do
    if (jh-jl <= 1) exit
    jc=(jh+jl)/2
    if (xpoint >= table(jc)) then
        jl=jc
    else
        jh=jc
    end if
  end do
  res=(xpoint-table(jh-1))/(table(jh)-table(jh-1))
  if(res<=0.5d0) then
    ixtp=jh-1
  else
    ixtp=jh
    res=res-1.d0
  endif
end if
end subroutine locate_in_table
!=============================================================================
