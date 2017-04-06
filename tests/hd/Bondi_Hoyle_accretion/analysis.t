!=============================================================================
subroutine write_analysis
! This is an example file how to use the analysis capability.  You can schedule
! this routine using the slot 5 in itsave, dtsave and ditsave.
! To use, just copy this file to your working directory and make your modifications.
use constants

! 15BAAL14My/22_0008 : to be reported
! Brand new version for BHL axisymm bhl_25D

!!$! printlog: calculates volume averaged mean values
!!$use mod_timing
!!$use mod_forest,only:nleafs,nleafs_active,nleafs_level

include 'amrvacdef.f'
character(len=20):: userconvert_type
!-----------------------------------------------------------------------------
integer, parameter :: nradii=128

logical :: fileopen
integer :: iigrid, igrid, level, iw
double precision :: volume(1:nlevelshi),  voltotal(nradii)
double precision :: re, ke, me, te
double precision :: dvolume(ixG^T), inBubble(ixG^T), lfac_pw(ixG^T), volumeflat(1:nlevelshi)
integer :: numlevels
double precision, dimension(1:nradii) :: dsum_send, dsum_recv
! double precision, dimension(2*1230) :: dcat_send
! double precision, allocatable, dimension(:) :: dcat_recv

character(len=80) :: filename
character(len=1024) :: line
logical, save :: opened2!=.false.
logical, save :: file_exists=.false.
integer :: amode, status(MPI_STATUS_SIZE)
double precision :: trcut
integer :: ii, jj, kk
double precision :: mdot(nradii)
double precision :: dS1(ixG^T)
double precision :: r_ring(nradii)
! double precision :: mdot_BHL ! For normalization
! double precision :: hull(1230,2), cvx_hull(1230,2)
! double precision :: press(ixG^T)
! double precision :: kh, ah, bh ! from y-kh=sqrt(ah2+bh*x2)
! double precision :: ap, bp     ! from y   =ap*x2+bp
! double precision :: NN, x2, x4, y2, yx2, cc,rshock,x2y2, x, yp, yh ! For parabolic fit of the shock front
! double precision :: slope(1230)
!-----------------------------------------------------------------------------
mdot=zero
! allocate(dcat_recv(npe*2*1230))
!print*, 'O', igridstail
do iigrid=1,igridstail; igrid=igrids(iigrid);
   dS1(ixM^T)=pgeo(igrid)%surfaceC1(ixM^T)
   ! r_ring(1)=1.05d-1
   ! mdot(1)=mdot(1)+sum(pw(igrid)%w(ixM^T,m1_)*dS1(ixM^T)*(-one),&
   !    MASK=(px(igrid)%x(ixM^T,1)>r_ring(1)-half*5.d-3 &
   !    .and. px(igrid)%x(ixM^T,1)<r_ring(1)+half*5.d-3))
   do ii=1,nradii
      r_ring(ii)=(xprobmin1*(qst/(one+half*logG)))*qst**(dble(ii-1))
      mdot(ii)=mdot(ii)+sum(pw(igrid)%w(ixM^T,m1_)*dS1(ixM^T)*(-one),&
         MASK=(px(igrid)%x(ixM^T,1)>r_ring(ii)*(one-half*logG) &
         .and. px(igrid)%x(ixM^T,1)<r_ring(ii)*(one+half*logG)))
      ! r_ring(ii)=(xprobmin1/(one-half*logG))*&
      !      ((one+half*logG)/(one-half*logG))**(dble(ii-1))
      ! mdot(ii)=mdot(ii)+sum(pw(igrid)%w(ixM^T,m1_)*dS1(ixM^T)*(-one),&
      !    MASK=(px(igrid)%x(ixM^T,1)>r_ring(ii)*(one-half*logG) &
      !    .and. px(igrid)%x(ixM^T,1)<r_ring(ii)*(one+half*logG)))
   enddo
   ! - - - - - - Shock hull - - - - - -
   ! call getpthermal(pw(igrid)%w,px(igrid)%x,ixG^LL,ixM^LL,press)
   ! !print*, igrid
   ! do jj=ixMlo1,ixMhi1 ! jj_prev+1,jj_prev+count()
   !    do kk=ixMlo2,ixMhi2
   !       ! Within 10 accretion radii, identify the shock w/ the temperature border 5.d4 and the density border 2*eqpar(rho0_)
   !       if ((any((press(ixM^T)/pw(igrid)%w(ixM^T,rho_))*eqpar(Tscale_)<5.d4 .and. &
   !             dsqrt((px(igrid)%x(ixM^T,1)*dsin(px(igrid)%x(ixM^T,2))-px(igrid)%x(jj,kk,1)*dsin(px(igrid)%x(jj,kk,2)))**two+&
   !                   (px(igrid)%x(ixM^T,1)*dcos(px(igrid)%x(ixM^T,2))-px(igrid)%x(jj,kk,1)*dcos(px(igrid)%x(jj,kk,2)))**two)<&
   !                      0.9d0*dsqrt(two)*px(igrid)%x(jj,kk,1)*logG) .and. &
   !          any((press(ixM^T)/pw(igrid)%w(ixM^T,rho_))*eqpar(Tscale_)>5.d4 .and. &
   !             dsqrt((px(igrid)%x(ixM^T,1)*dsin(px(igrid)%x(ixM^T,2))-px(igrid)%x(jj,kk,1)*dsin(px(igrid)%x(jj,kk,2)))**two+&
   !                   (px(igrid)%x(ixM^T,1)*dcos(px(igrid)%x(ixM^T,2))-px(igrid)%x(jj,kk,1)*dcos(px(igrid)%x(jj,kk,2)))**two)<&
   !                      0.9d0*dsqrt(two)*px(igrid)%x(jj,kk,1)*logG)) .or. &
   !          (any(pw(igrid)%w(ixM^T,rho_)<two*eqpar(rho0_) .and.  &
   !             dsqrt((px(igrid)%x(ixM^T,1)*dsin(px(igrid)%x(ixM^T,2))-px(igrid)%x(jj,kk,1)*dsin(px(igrid)%x(jj,kk,2)))**two+&
   !                   (px(igrid)%x(ixM^T,1)*dcos(px(igrid)%x(ixM^T,2))-px(igrid)%x(jj,kk,1)*dcos(px(igrid)%x(jj,kk,2)))**two)<&
   !                      0.9d0*dsqrt(two)*px(igrid)%x(jj,kk,1)*logG) .and. &
   !          any(pw(igrid)%w(ixM^T,rho_)>two*eqpar(rho0_) .and.  &
   !             dsqrt((px(igrid)%x(ixM^T,1)*dsin(px(igrid)%x(ixM^T,2))-px(igrid)%x(jj,kk,1)*dsin(px(igrid)%x(jj,kk,2)))**two+&
   !                   (px(igrid)%x(ixM^T,1)*dcos(px(igrid)%x(ixM^T,2))-px(igrid)%x(jj,kk,1)*dcos(px(igrid)%x(jj,kk,2)))**two)<&
   !                      0.9d0*dsqrt(two)*px(igrid)%x(jj,kk,1)*logG))) then
   !       hull(ind,1)=px(igrid)%x(jj,kk,1)*dsin(px(igrid)%x(jj,kk,2))
   !       hull(ind,2)=px(igrid)%x(jj,kk,1)*dcos(px(igrid)%x(jj,kk,2))
   !       ind=ind+1
   !       endif
   !    enddo
   ! enddo
enddo

!if (ind>1230 .or. ind==0) call mpistop("Ne manqueriez-vous pas de points par hasard?")

dsum_send(1:nradii)=mdot(1:nradii)
! dcat_send(1:1230)=hull(1:1230,1)
! dcat_send(1231:2460)=hull(1:1230,2)

call MPI_REDUCE(dsum_send,dsum_recv,nradii,MPI_DOUBLE_PRECISION, &
                MPI_SUM,0,icomm,ierrmpi)
! call MPI_GATHER(dcat_send,2460,MPI_DOUBLE_PRECISION,dcat_recv,2460,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)

! hull=0

if (mype==0) then

!    NN=0
!    rshock=0
!    ind=0
!    kk=0
!    x2=0
!    x4=0
!    yp=0
!    yx2=0
!    do ii=1,npe
!       do jj=1,size(dcat_send)/2 ! always 2
!          if (dcat_recv((ii-1)*size(dcat_send)+jj)/=zero .and. dcat_recv((ii-1)*size(dcat_send)+size(dcat_send)/2+jj)/=zero) then
!             ! Use the points in the close forward angle (3  orthoradial cells on one  side) to compute rshock
!             if     (datan(dabs(dcat_recv((ii-1)*size(dcat_send)+jj)/dcat_recv((ii-1)*size(dcat_send)+size(dcat_send)/2+jj)))<&
!                 3.d0*((xprobmax2-xprobmin2)/dble(nxlone2)) .and. &
!                 dcat_recv((ii-1)*size(dcat_send)+size(dcat_send)/2+jj)>zero) then
!                rshock=rshock+dsqrt(dcat_recv((ii-1)*size(dcat_send)+jj)**two+dcat_recv((ii-1)*size(dcat_send)+size(dcat_send)/2+jj)**two)
!                ind=ind+1
!             endif
!             ! Use the points in the close forward angle (10 orthoradial cells on each side) to compute the curvature radius
!             if (datan(dabs(dcat_recv((ii-1)*size(dcat_send)+jj)/dcat_recv((ii-1)*size(dcat_send)+size(dcat_send)/2+jj)))<&
!                 10.d0*((xprobmax2-xprobmin2)/dble(nxlone2))) then
!                 kk=kk+1+1
!                 x2 = x2 +  (dcat_recv((ii-1)*size(dcat_send)+jj))**2.d0 + (-dcat_recv((ii-1)*size(dcat_send)+jj))**2.d0
!                 x4  = x4 + (dcat_recv((ii-1)*size(dcat_send)+jj))**4.d0 + (-dcat_recv((ii-1)*size(dcat_send)+jj))**4.d0
!                 yp  = yp  + dcat_recv((ii-1)*size(dcat_send)+size(dcat_send)/2+jj) + &
!                             dcat_recv((ii-1)*size(dcat_send)+size(dcat_send)/2+jj)
!                 yx2= yx2+ ( dcat_recv((ii-1)*size(dcat_send)+jj))**2.d0 * dcat_recv((ii-1)*size(dcat_send)+size(dcat_send)/2+jj) +&
!                           (-dcat_recv((ii-1)*size(dcat_send)+jj))**2.d0 * dcat_recv((ii-1)*size(dcat_send)+size(dcat_send)/2+jj)
!              endif
!             NN = NN + 2
!             hull(NN,1)=   dcat_recv((ii-1)*size(dcat_send)+jj)
!             hull(NN,2)=   dcat_recv((ii-1)*size(dcat_send)+size(dcat_send)/2+jj)
!             hull(NN+1,1)=-dcat_recv((ii-1)*size(dcat_send)+jj)
!             hull(NN+1,2)= dcat_recv((ii-1)*size(dcat_send)+size(dcat_send)/2+jj)
!          endif
!       enddo
!    enddo
!    deallocate(dcat_recv)
!    ap=(yp*x2-dble(kk)*yx2)/(x2*x2-dble(kk)*x4)
!    rshock=rshock/dble(ind)
!    hull=hull*normvar(0)
!    call convex_hull(size(hull(:,1)),hull(:,:),cvx_hull(:,:))
!    hull=hull/normvar(0)

   ! rshock=rshock/dble(ind)
   ! bp=rshock
   ! NN=0
   ! x=0
   ! x2=0
   ! yp=0
   ! jj=0
   ! do ii=1,10
   !    NN=maxloc(hull(1,:),DIM=1)
   !    print*, 'aásasd', hull(1,NN)*normvar(0)
   !    jj = jj + 1
   !    x  = x + hull(1,NN)
   !    x2 = x2+ hull(1,NN)**two
   !    yp  = yp + hull(2,NN)
   !    hull(1,NN)=zero
   ! enddo
   ! kh=(yp*x2-yp*x*x)/(dble(jj)*x2-x*x)
   ! ah=kh-rshock
   ! NN=0
   ! x=0
   ! x2=0
   ! x4=0
   ! yh=0
   ! yp=0
   ! y2=0
   ! x2y2=0
   ! yx2=0
   ! do ii=1,npe
   !    do jj=1,size(dcat_send)/2 ! always 2
   !       if (dcat_recv((ii-1)*size(dcat_send)+jj)/=zero .and. dcat_recv((ii-1)*size(dcat_send)+size(dcat_send)/2+jj)/=zero) then
   !          NN  = NN  + one + one
   !          x2 = x2 +  (dcat_recv((ii-1)*size(dcat_send)+jj))**2.d0 + (-dcat_recv((ii-1)*size(dcat_send)+jj))**2.d0
   !          x4  = x4 + (dcat_recv((ii-1)*size(dcat_send)+jj))**4.d0 + (-dcat_recv((ii-1)*size(dcat_send)+jj))**4.d0
   !          yh  = yh  + (dcat_recv((ii-1)*size(dcat_send)+size(dcat_send)/2+jj)-kh) + &
   !                    (dcat_recv((ii-1)*size(dcat_send)+size(dcat_send)/2+jj)-kh)
   !          yp  = yp  + dcat_recv((ii-1)*size(dcat_send)+size(dcat_send)/2+jj) + &
   !                      dcat_recv((ii-1)*size(dcat_send)+size(dcat_send)/2+jj)
   !          y2 = y2 + (dcat_recv((ii-1)*size(dcat_send)+size(dcat_send)/2+jj)-kh)**2.d0 + &
   !                    (dcat_recv((ii-1)*size(dcat_send)+size(dcat_send)/2+jj)-kh)**2.d0
   !          x2y2=x2y2+ (dcat_recv((ii-1)*size(dcat_send)+jj))**2.d0 * ((dcat_recv((ii-1)*size(dcat_send)+size(dcat_send)/2+jj)-kh)**2.d0) + &
   !                    (-dcat_recv((ii-1)*size(dcat_send)+jj))**2.d0 * ((dcat_recv((ii-1)*size(dcat_send)+size(dcat_send)/2+jj)-kh)**2.d0)
   !          yx2= yx2+ ( dcat_recv((ii-1)*size(dcat_send)+jj))**2.d0 * dcat_recv((ii-1)*size(dcat_send)+size(dcat_send)/2+jj) +&
   !                    (-dcat_recv((ii-1)*size(dcat_send)+jj))**2.d0 * dcat_recv((ii-1)*size(dcat_send)+size(dcat_send)/2+jj)
   !          ! Tangent of the angle relative to normal to the polar axis (the "horizontal")
   !          !slope(ind)=dcat_recv((ii-1)*size(dcat_send)+size(dcat_send)/2+jj)/(dcat_recv((ii-1)*size(dcat_send)+jj))
   !          ! slope(ind)=dabs(datan(dabs(-two*aa*dcat_recv((ii-1)*size(dcat_send)+jj)))-&
   !          ! datan(dabs(dcat_recv((ii-1)*size(dcat_send)+jj)/dcat_recv((ii-1)*size(dcat_send)+size(dcat_send)/2+jj))))
   !          ! print*, slope(ind), dabs(datan(dabs(-two*aa*dcat_recv((ii-1)*size(dcat_send)+jj)))-&
   !          ! datan(dabs(dcat_recv((ii-1)*size(dcat_send)+size(dcat_send)/2+jj)/dcat_recv((ii-1)*size(dcat_send)+jj))))
   !          ! ind=ind+1
   !       endif
   !    enddo
   ! enddo
   !
   ! close(123)
   ! deallocate(dcat_recv)
   ! ! Even 2nd order polynomial => bb, x, x3, yx = 0
   ! ! aa= (NN*yx2-y*x2)/(x4*NN-x2*x2)
   ! ! cc = (yx2*x2-y*x4)/(x2*x2-NN*x4)
   ! bh=(x2y2-ah**two*x2)/x4
   ! ap=(yp*x2-NN*yx2)/(x2*x2-NN*x4)
   ! !bp=(yp*x4-yx2*x2)/(NN*x4-x2*x2)
   ! !bb=(x2y2-x2*y2)/(x4-x2*x2)
   ! !cc=(x2*x2y2-x4*y2)/(x2*x2-x4*NN)
   ! print*, kh*normvar(0), (ah*normvar(0))**two, bh
   ! print*, ap/normvar(0), bp*normvar(0)
   ! ! - - - - - - Shock hull - - - - - -

   mdot(1:nradii)=dsum_recv(1:nradii)

   ! if(.not. saveprim) then
   !    call mpistop("Wrong normalization in analysis.t")
   ! endif

   ! Mdot in units of Mdot_BHL=pi*Racc^2*rho_inf*v_inf
   ! Adimensioned mdot_BHL (is equal to dpi*2 here w/ Racc, rho0/2 and vel as normalization units)
   ! dpi/two since only half of the theta space (ie only from 0 to pi => half area)
   ! mdot_BHL=(dpi/two)*((two*eqpar(gm_))/eqpar(vel_)**two)*eqpar(rho0_)*eqpar(vel_)
   ! Factor two*dpi to integrate over phi (2.5D --> 3D)
	mdot=mdot*(two*dpi) ! to account for the azimuthal phi third direction
	mdot=mdot/dpi ! to normalize w/ the dimensionless BHL Mdot

   if (.not.opened2) then
      ! generate filename
      write(filename,"(a,a)") "analysis.txt"
      inquire(file=filename,exist=file_exists)

      if (.not. file_exists) then
         open(unit=unitanalysis,file=filename,status='unknown',access='append')
         !!write(unitanalysis,"(a)",advance="no") trim('# t [years] rho u1 u2 u3 p b1 b2 b3 tr2 lfac xi evolve at (r,z)=(')
         !!write(unitanalysis,*) ' '
         !!write(unitanalysis,*) ' '
         !write(unitanalysis,'(2(es14.6))',advance="no") myx(1)*normvar(0),myx(2)*normvar(0)
         !write(unitanalysis,"(a)") trim(') [cm]')
      else
         open(unit=unitanalysis,file=filename,status='unknown',access='append')
      end if
      opened2=.true.
   end if

   opened2=.true. ! bricolage sans lequel opened2 revient toujours a .false.

   write(unitanalysis,*) t, mdot(1)

   !!write(unitanalysis,*) 'Time :', t!*normt
   ! rshock is computed from the avg of the forefront points
   ! le maxval est le point le plus en avant de l enveloppe convexe (generalement > rshock)
   ! write(unitanalysis,*) 'Shock position (in R_E)   :', rshock/eqpar(Racc_), maxval(hull(:,2))/eqpar(Racc_)
   ! write(unitanalysis,*) 'Avg & discrepancy         :', half*(rshock+maxval(hull(:,2)))/eqpar(Racc_), 1.d2*(dabs(rshock-maxval(hull(:,2)))/rshock), '%'
   ! write(unitanalysis,*) 'Curvature radius at front :', one/((two*dabs(ap))/eqpar(Racc_)) ! Factor two from the 1st order derivative of a second order polynome
   !!write(unitanalysis,*) '  Radius r  |  Mdot(r)'
   !!do ii=1,nradii
      ! r_ring in units of R_E
      ! Mdot in units of Mdot_BHL=pi*Racc^2*rho_inf*v_inf
      !!write(unitanalysis,'(7(es12.4))') r_ring(ii), mdot(ii)!r_ring(ii)*eqpar(Racc_), mdot(ii)
   !!enddo
   !!write(unitanalysis,*) ' '

	! Save only the last
   ! if(t>=tmax .or. it>=itmax) then
   ! ! Tricky way to write the convex hull points in order and
   ! ! so to be able to link them with a line (without quinconce effect)
   ! open(2,file=TRIM(filenameanalysis)//'convex_hull')
   ! do ii=1,size(hull(:,1))
   !    if (cvx_hull(ii,1)/=0.d0 .and. cvx_hull(ii,2)/=0.d0) then
   !       write(2,*)  cvx_hull(ii,1), cvx_hull(ii,2)
   !    endif
   ! enddo
   ! ! (avoid reverse loops as much as possible -
   ! ! depending on the optimization flag during the compilation -O0, -O1, -O2 or -O3, it can make nasty things)
   ! do jj=1,size(hull(:,1))
   !    ii=size(hull(:,1))-(jj-1)
   !    if (cvx_hull(ii,1)/=0.d0 .and. cvx_hull(ii,2)/=0.d0) then
   !       write(2,*) -cvx_hull(ii,1), cvx_hull(ii,2)
   !    endif
   ! enddo
   ! close(2)
   ! endif

end if ! mype == 0

call MPI_BARRIER(icomm,ierrmpi)

end subroutine write_analysis
!=============================================================================



!=============================================================================
subroutine convex_hull(ndots,x,y)
implicit none
double precision, dimension(ndots,2), intent(in) :: x
double precision, dimension(ndots,2), intent(out) :: y
integer, intent(in) :: ndots
integer :: i

y=0.d0

 ! open(1,file='x.dat')
 ! do i=1,ndots
 ! write(1,*)  x(i,1), x(i,2)
 ! write(1,*) -x(i,1), x(i,2)
 ! enddo
 ! close(1)

call algo(ndots,x,y)

! open(2,file='y.dat')
! do i=1,ndots
!    if (y(i,1)/=0.d0 .and. y(i,2)/=0.d0) then
!       write(2,*)  y(i,1), y(i,2)
!       write(2,*) -y(i,1), y(i,2)
!    endif
! enddo
! close(2)

end

!=============================================================================
! Cf Wikipedia article Quick hull for the algo
! y is too big but no easy way to work on unknown size objects in Fortran,
! contrary to C (with the vector objects eg, or append functions)
!=============================================================================
subroutine algo(ndots,x,y)
implicit none
integer, intent(in) :: ndots
double precision, dimension(ndots,2), intent(in) :: x
double precision, dimension(ndots,2), intent(out) :: y
double precision :: yline, yline0, ylinek, xc, yc, delx, dely, dist
double precision, dimension(ndots) :: distmax
integer, dimension(ndots) :: check
integer :: i,j,k1,k2
integer, dimension(ndots) :: keys, keystot
logical, dimension(ndots) :: dans

y=0.d0
distmax=0.d0
keys=0

! - - - - - - - - - - - - - - - - - - - - -
! First octogon
! - - - - - - - - - - - - - - - - - - - - -
! Upper segments
keys(1)=minloc(x(:,1),DIM=1)
keys(2)=minloc(x(:,1)-x(:,2),DIM=1)
keys(3)=maxloc(x(:,2),DIM=1)
keys(4)=maxloc(x(:,1)+x(:,2),DIM=1)
! Lower segments
keys(5)=maxloc(x(:,1),DIM=1)
keys(6)=maxloc(x(:,1)-x(:,2),DIM=1)
keys(7)=minloc(x(:,2),DIM=1)
keys(8)=minloc(x(:,1)+x(:,2),DIM=1)

! open(3,file='oct.dat')
! do i=1,8
! write(3,*) x(keys(i),1), x(keys(i),2)
! enddo
! close(3)

dans=.false.
do i=1,ndots
   if (any(keys(:)==i)) cycle

   do j=1,8
      k1=keys(j)
      if (j/=8) k2=keys(j+1)
      if (j==8) k2=keys(1)
      if ((x(i,1)-x(k1,1))*(x(i,1)-x(k2,1))>0.d0) cycle ! If point i not above/below segment (j,j+1)
      yline0=x(keys(1),2)+(x(keys(5),2)-x(keys(1),2)) * ( (x(i,1)-x(keys(1),1)) / (x(keys(5),1)-x(keys(1),1)) )
      if (j<5 .and. x(i,2)<yline0) cycle
      if (j>4 .and. x(i,2)>yline0) cycle
      delx=x(k2,1)-x(k1,1)
      dely=x(k2,2)-x(k1,2)
      yline = x(k1,2) + dely * ( (x(i,1)-x(k1,1)) / delx )
      if (yline<x(i,2) .and. x(i,2)<yline0) dans(i)=.true.
      if (yline>x(i,2) .and. x(i,2)>yline0) dans(i)=.true.
   enddo
enddo
! We keep only the points out of the octogon
! - - - - - - - - - - - - - - - - - - - - -

keystot=0

! - - - - - - - - - - -
! Suppress doublons
! - - - - - - - - - - -
j=1
keystot(1)=keys(1)
do i=2,8
   if (any(keystot(:)==keys(i))) cycle
   j=j+1
   keystot(j)=keys(i)
enddo
keys=keystot
! - - - - - - - - - - -

check=1
do while(count(check/=0)/=0)
print*, 'loop'
distmax=0.d0
check=0
! Loop on vertex
do j=1,count(keys/=0)
      k1=keys(j)
      ! Initial line which separates upper and lower half, evaluated for the vertex (upper or lower segment?)
      ylinek=x(minloc(x(:,1),DIM=1),2)+(x(maxloc(x(:,1),DIM=1),2)-x(minloc(x(:,1),DIM=1),2)) * &
              ( (x(k1,1)-x(minloc(x(:,1),DIM=1),1)) / (x(maxloc(x(:,1),DIM=1),1)-x(minloc(x(:,1),DIM=1),1)) )
      if (j/=count(keys/=0)) k2=keys(j+1)
      if (j==count(keys/=0)) k2=keys(1)   ! Congruence
      if (k1==k2) then
         print*, 'WTF1?'
         !!stop
      endif
      ! Loop on all points
      do i=1,ndots
         if(dans(i) .EQV. .true.) cycle ! not if inside
         if (any(keys(:)==i)) cycle     ! not if it is itself a vertex
         if ((x(i,1)-x(k1,1))*(x(i,1)-x(k2,1))>0.d0) cycle ! If point i not above/below segment (j,j+1)
         yline=x(minloc(x(:,1),DIM=1),2)+(x(maxloc(x(:,1),DIM=1),2)-x(minloc(x(:,1),DIM=1),2)) * &
              ( (x(i,1)-x(minloc(x(:,1),DIM=1),1)) / (x(maxloc(x(:,1),DIM=1),1)-x(minloc(x(:,1),DIM=1),1)) )
         if ((x(k1,2)>ylinek .or. k1==minloc(x(:,1),DIM=1)) .and. x(i,2)<yline) cycle ! Upper segment / Lower half
         if ((x(k1,2)<ylinek .or. k1==maxloc(x(:,1),DIM=1)) .and. x(i,2)>yline) cycle ! Lower segment / Upper half
         delx=x(k2,1)-x(k1,1)
         dely=x(k2,2)-x(k1,2)
         xc = (x(i,2)-x(k1,2)+x(k1,1)*(dely/delx)+x(i,1)*(delx/dely))/(dely/delx+delx/dely)
         yc = x(k1,2) + dely * ( (xc-x(k1,1)) / delx )
         dist=(x(i,1)-xc)**2.d0+(x(i,2)-yc)**2.d0 ! Distance au projete orthogonal
         ! Store the maximal distance for this segment
         if (dist>distmax(j)) then
            check(j)=i
            distmax(j)=dist
         endif
      enddo
      ! Reorder the vertex, inserting the new ones
      if (check(j)/=0) then
         if (j/=count(keys/=0)) then
            ! Shift all the ones below...
            do i=count(keystot/=0),j+1+count(check/=0)-1,-1!1,8-j
               keystot(i+1)=keystot(i)
            enddo
            ! ... and insert the new one
            keystot(j+1+count(check/=0)-1)=check(j)
         else
            keystot(count(keystot/=0)+1)=check(j)
         endif
      endif
   enddo
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - -
   ! CLEAN POINT INSIDE ADDING NEW POINTS TO " DANS"
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - -
   keys=keystot
   do i=1,ndots
      if(dans(i) .EQV. .true.) cycle ! not if inside
      if (any(keys(:)==i)) cycle
      do j=1,count(keys/=0)
         k1=keys(j)
         if (j/=count(keys/=0)) k2=keys(j+1)
         if (j==count(keys/=0)) k2=keys(1)
         if ((x(i,1)-x(k1,1))*(x(i,1)-x(k2,1))>0.d0) cycle ! If point i not above/below segment (j,j+1)
         ylinek=x(minloc(x(:,1),DIM=1),2)+(x(maxloc(x(:,1),DIM=1),2)-x(minloc(x(:,1),DIM=1),2)) * &
              ( (x(k1,1)-x(minloc(x(:,1),DIM=1),1)) / (x(maxloc(x(:,1),DIM=1),1)-x(minloc(x(:,1),DIM=1),1)) )
         yline=x(minloc(x(:,1),DIM=1),2)+(x(maxloc(x(:,1),DIM=1),2)-x(minloc(x(:,1),DIM=1),2)) * &
              ( (x(i,1)-x(minloc(x(:,1),DIM=1),1)) / (x(maxloc(x(:,1),DIM=1),1)-x(minloc(x(:,1),DIM=1),1)) )
         if ((x(k1,2)>ylinek .or. k1==minloc(x(:,1),DIM=1)) .and. x(i,2)<yline) cycle ! Upper segment / Lower half
         if ((x(k1,2)<ylinek .or. k1==maxloc(x(:,1),DIM=1)) .and. x(i,2)>yline) cycle ! Lower segment / Upper half
         delx=x(k2,1)-x(k1,1)
         dely=x(k2,2)-x(k1,2)
         yline = x(k1,2) + dely * ( (x(i,1)-x(k1,1)) / delx )
         if (yline>x(i,2) .and. (x(k1,2)>ylinek .or. k1==minloc(x(:,1),DIM=1))) dans(i)=.true.
         if (yline<x(i,2) .and. (x(k1,2)<ylinek .or. k1==maxloc(x(:,1),DIM=1))) dans(i)=.true.
      enddo
   enddo
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - -
enddo

do i=1,count(keys/=0)
   y(i,1)=x(keys(i),1)
   y(i,2)=x(keys(i),2)
enddo

end
!=============================================================================
