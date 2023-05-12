!> Program to extrapolate linear force-free fields in 3D Cartesian coordinates,
!> based on exact Green function method (Chiu & Hilton 1977 ApJ 212,873).
!>
!> Usage:
!> 1 In the subroutine usr_set_parameters of mod_usr.t:
!>  To extrapolate a linear force free field from a observed magnetogram 
!>  prepared in a data file, e.g., 'hmiM720sxxxx.dat' replace 
!>  call init_bc_fff_data('hmiM720sxxxx.dat',unit_length,unit_magneticfield)
!>  'hmiM720sxxxx.dat' must be a binary file containing nx1,nx2,xc1,xc2,dxm1,
!>  dxm2, Bz0(nx1,nx2). Integers nx1 and nx2 give the resolution of the 
!>  uniform-grid magentogram. Others are double-precision floats. xc1 and xc2
!>  are coordinates of the central point of the magnetogram. dxm1 and dxm2 
!>  are the cell sizes for each direction, Bz0 is the vertical conponent 
!>  of magetic field on the solar surface from observations.
!>2 In the subroutine usr_init_one_grid of mod_usr.t,
!>  add lines like:
!>
!>  double precision :: Bf(ixG^S,1:ndir), alpha, zshift
!>
!>  alpha=0.d0     ! potential field
!>  !alpha=0.08d0  ! non-potential linear force-free field
!>  zshift=0.05d0  ! lift your box zshift heigher to the bottom magnetogram
!>  call calc_lin_fff(ixG^L,ix^L,Bf,x,alpha,zshift) 
!>
!>3 Notice that the resolution of input magnetogram must be better than the best
!>  resolution of your AMR grid to have a good behavior close to the bottom layer
module mod_lfff
  implicit none
  
  integer, save :: nx1,nx2
  double precision, save :: Bzmax,darea
  double precision, allocatable, save :: Bz0(:,:)
  double precision, allocatable, save :: xa1(:),xa2(:)
  
contains

  subroutine init_b_fff_data(magnetogramname,qLunit,qBunit)
    use mod_global_parameters

    double precision, intent(in) :: qLunit,qBunit
    double precision :: xc1,xc2,dxm1,dxm2
    integer, dimension(MPI_STATUS_SIZE) :: statuss
    integer :: file_handle,i
    character(len=*), intent(in) :: magnetogramname
    logical :: aexist
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
    do i=1,nx1 
      xa1(i) = xc1 + (dble(i) - dble(nx1)/2.d0 - 0.5d0)*dxm1
    enddo
    do i=1,nx2
      xa2(i) = xc2 + (dble(i) - dble(nx2)/2.d0 - 0.5d0)*dxm2
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

{^IFTHREED
  subroutine calc_lin_fff(ixI^L,ixO^L,Bf,x,alpha,zshift,idir)
  ! PURPOSE: 
  ! Calculation to determine linear FFF from the field on 
  ! the lower boundary (Chiu and Hilton 1977 ApJ 212,873). 
  ! NOTE: Only works for Cartesian coordinates 
  ! INPUT: Bf,x
  ! OUTPUT: updated b in w 
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    integer, optional, intent(in) :: idir
    double precision, intent(in) :: x(ixI^S,1:ndim),alpha,zshift
    double precision, intent(inout) :: Bf(ixI^S,1:ndir)

    double precision, dimension(ixO^S) :: cos_az,sin_az,zk,bigr,r,r2,r3,cos_ar,sin_ar,g,dgdz
    double precision :: gx(ixO^S,1:ndim),twopiinv
    integer :: idim,ixp1,ixp2

    Bf=0.d0
    twopiinv = 0.5d0/dpi*Bzmax*darea
    ! get cos and sin arrays
    zk(ixO^S)=x(ixO^S,3)-xprobmin3+zshift
    cos_az(ixO^S)=dcos(alpha*zk(ixO^S))
    sin_az(ixO^S)=dsin(alpha*zk(ixO^S))
    ! looping Bz0 pixels
    do ixp2=1,nx2
      do ixp1=1,nx1
        bigr(ixO^S)=dsqrt((x(ixO^S,1)-xa1(ixp1))**2+&
                          (x(ixO^S,2)-xa2(ixp2))**2)
        r2=bigr**2+zk**2
        r=dsqrt(r2)
        r3=r**3
        cos_ar=dcos(alpha*r)
        sin_ar=dsin(alpha*r)
        where(bigr/=0.d0)
          bigr=1.d0/bigr
        end where
        where(r/=0.d0)
          r=1.d0/r
        end where
        where(r2/=0.d0)
          r2=1.d0/r2
        end where
        where(r3/=0.d0)
          r3=1.d0/r3
        end where
        g=(zk*cos_ar*r-cos_az)*bigr
        dgdz=(cos_ar*(r-zk**2*r3)-alpha*zk**2*sin_ar*r2+alpha*sin_az)*bigr
        do idim=1,ndim
          if(present(idir).and.idim/=idir) cycle
          select case(idim)
          case(1)
            Bf(ixO^S,1)=Bf(ixO^S,1)+Bz0(ixp1,ixp2)*((x(ixO^S,1)-xa1(ixp1))*dgdz(ixO^S)&
                     +alpha*g(ixO^S)*(x(ixO^S,2)-xa2(ixp2)))*bigr(ixO^S)
          case(2)
            Bf(ixO^S,2)=Bf(ixO^S,2)+Bz0(ixp1,ixp2)*((x(ixO^S,2)-xa2(ixp2))*dgdz(ixO^S)&
                     -alpha*g(ixO^S)*(x(ixO^S,1)-xa1(ixp1)))*bigr(ixO^S)
          case(3)
            Bf(ixO^S,3)=Bf(ixO^S,3)+Bz0(ixp1,ixp2)*(zk(ixO^S)*cos_ar(ixO^S)*r3(ixO^S)+alpha*&
                                        zk(ixO^S)*sin_ar(ixO^S)*r2(ixO^S))
          end select
        end do
      end do
    end do
    Bf(ixO^S,:)=Bf(ixO^S,:)*twopiinv

  end subroutine calc_lin_fff

  subroutine get_potential_field_potential(ixI^L,ixO^L,potential,x,zshift)
  ! PURPOSE: 
  ! Calculation scalar potential of potential field given
  ! Bz at photosphere (Schmidt 1964 NASSP). 
  ! NOTE: Only works for Cartesian coordinates 
  ! INPUT: x,zshift
  ! OUTPUT: potential
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim),zshift
    double precision, intent(inout) :: potential(ixI^S)

    double precision, dimension(ixO^S) :: zk,bigr
    integer :: ixp1,ixp2

    zk(ixO^S)=x(ixO^S,3)-xprobmin3+zshift
    potential=0.d0
    ! looping Bz0 pixels see equation (2)
    !$OMP PARALLEL DO PRIVATE(bigr) REDUCTION(+:potential)
    do ixp2=1,nx2
      do ixp1=1,nx1
        bigr(ixO^S)=dsqrt((x(ixO^S,1)-xa1(ixp1))**2+&
                          (x(ixO^S,2)-xa2(ixp2))**2+&
                          zk(ixO^S)**2)
        potential(ixO^S)=potential(ixO^S)+0.5d0*Bz0(ixp1,ixp2)/bigr*darea/dpi
      end do
    end do
    !$OMP END PARALLEL DO
  end subroutine get_potential_field_potential

  subroutine get_potential_field_potential_sphere(ixI^L,x,potential,nth,nph,magnetogram,theta,phi,r_sphere)
  ! PURPOSE: 
  ! Calculation scalar potential of potential field given
  ! Bz at photosphere (Schmidt 1964 NASSP). 
  ! NOTE: Only works for spherical coordinates 
  ! OUTPUT: potential
    use mod_global_parameters
    integer, intent(in) :: ixI^L,nth,nph
    real*8, intent(in) :: x(ixI^S,1:ndim)
    ! magnetogram Br on photosphere
    real*8, intent(in) :: magnetogram(nth,nph)
    ! theta and phi grid of the photospheric magnetogram
    real*8, intent(in) :: theta(nth),phi(nph)
    ! radius of photosphere
    real*8, intent(in) :: r_sphere
    real*8, intent(out) :: potential(ixI^S)

    real*8 :: area(nth),distance(ixI^S),dtheta_half,dphi,inv2pi
    integer :: ix1,ix2

    potential=0.d0
    ! assume uniformly discretized theta and phi
    dtheta_half=0.5d0*(theta(2)-theta(1))
    dphi=phi(2)-phi(1)
    area(1:nth)=2.d0*r_sphere**2*sin(theta(1:nth))*sin(dtheta_half)*sin(dphi)
    inv2pi=-1.d0/(2.d0*dpi)
    
    !$OMP PARALLEL DO PRIVATE(distance) REDUCTION(+:potential)
    do ix2=1,nph,2
      do ix1=1,nth,2
        distance(ixI^S)=sqrt(x(ixI^S,1)**2+r_sphere**2-2.d0*x(ixI^S,1)*r_sphere*&
        (sin(x(ixI^S,2))*sin(theta(ix1))*cos(phi(ix2)-x(ixI^S,3))+cos(x(ixI^S,2))*&
        cos(theta(ix1))))
        where(distance(ixI^S)/=0.d0)
          distance(ixI^S)=1.d0/distance(ixI^S)
        end where
        potential(ixI^S)=potential(ixI^S)+inv2pi*magnetogram(ix1,ix2)*distance(ixI^S)*area(ix1)
      end do
    end do
    !$OMP END PARALLEL DO

  end subroutine get_potential_field_potential_sphere

  !> get potential magnetic field energy given normal B on all boundaries
  subroutine potential_field_energy_mg(benergy)
  ! NOTE: Solve Poisson equation of scalar potential using multigrid solver
  ! Only works for 3D Cartesian coordinates 
    use mod_global_parameters
    use mod_multigrid_coupling
    use mod_forest
    use mod_geometry
    real*8, intent(out) :: benergy
    real*8  :: block_Benergy(ixM^T), mype_Benergy
    real*8, allocatable :: tmp(:,:,:)
    integer :: iigrid, igrid, idir
    integer :: id,nc, lvl, ixI^L

    call get_potential_field_potential_mg()
    ixImin^D=ixMlo^D-1;
    ixImax^D=ixMhi^D+1;
    allocate(tmp(ixI^S))
    mype_Benergy=0.d0
    do iigrid=1,igridstail; igrid=igrids(iigrid);
      block=>ps(igrid)
      ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
      id = igrid_to_node(igrid, mype)%node%id
      lvl= mg%boxes(id)%lvl
      nc = mg%box_size_lvl(lvl)
      block_Benergy=0.d0
      do idir=1,ndim
        call gradient(mg%boxes(id)%cc({0:nc+1},mg_iphi),ixI^L,ixM^LL,idir,tmp)
        block_Benergy(ixM^T)=block_Benergy(ixM^T)+tmp(ixM^T)**2
      end do
      mype_Benergy=mype_Benergy+sum(0.5d0*block_Benergy(ixM^T)*block%dvolume(ixM^T))
    end do

    call MPI_ALLREDUCE(mype_Benergy,benergy,1,MPI_DOUBLE_PRECISION,MPI_SUM,icomm,ierrmpi)

  end subroutine potential_field_energy_mg

  !>  Solve Poisson equation of scalar potential using multigrid solver
  subroutine get_potential_field_potential_mg()
    use mod_global_parameters
    use mod_multigrid_coupling
    double precision :: max_residual, res
    integer :: i, max_its

    mg%operator_type = mg_laplacian
    max_its=50
    max_residual=1.d-10

    ! Set boundary conditions
    mg%bc(:, mg_iphi)%bc_type = mg_bc_neumann

    do i=1,2*ndim
       mg%bc(i, mg_iphi)%boundary_cond => multigrid_bc
    end do

    ! Solve laplacian(phi) = 0
    do i=1,max_its
       !call mg_fas_vcycle(mg, max_res=res)
       call mg_fas_fmg(mg,.true.,max_res=res)
       if(mype==0) write(*,*) 'mg iteration',i,'residual:', res
       if(res < max_residual) exit
    end do
    if(res > max_residual) call mpistop("get potential multigrid: no convergence")

  end subroutine get_potential_field_potential_mg

  !> To set boundary condition on physical boundaries for mg Poisson solver
  subroutine multigrid_bc(box, nc, iv, nb, bc_type, bc)
    use mod_global_parameters
    use mod_multigrid_coupling
    type(mg_box_t), intent(in)    :: box
    integer, intent(in)           :: nc
    integer, intent(in)           :: iv      !< Index of variable
    integer, intent(in)           :: nb      !< number of boundary from 1 to 6 for 3D
    integer, intent(out)          :: bc_type !< Type of b.c.
    ! mg boundary values
    double precision, intent(out) :: bc(nc, nc)
    ! mg coordinates of boundary layer
    double precision              :: rr(nc, nc, 3)
    double precision :: rmina,rminb,rmaxa,rmaxb,xmina,xminb,xmaxa,xmaxb
    ! store normal magnetic field on the boundary
    double precision :: wbn(ixG^T)
    double precision, allocatable :: xcoarse(:,:,:)
    integer :: iigrid,igrid,ix^D,idir,ixbca,ixbcb,ixbcn,dlvl,wnc

    bc_type = mg_bc_neumann

    call mg_get_face_coords(box, nb, nc, rr)

    idir=(nb+1)/2
    bc=0.d0
    ! send normal B on boundaries from AMRVAC to mg Poisson solver
    ! for Neumann boundary condition
    select case(idir)
    case(1)
      if(mod(nb,2)==0) then
        ixbcn=ixMhi1
      else
        ixbcn=ixMlo1-1
      end if
      rmina=rr(1,1,2)-0.5d0*box%dr(2)
      rmaxa=rr(nc,1,2)+0.5d0*box%dr(2)
      rminb=rr(1,1,3)-0.5d0*box%dr(3)
      rmaxb=rr(1,nc,3)+0.5d0*box%dr(3)
      do iigrid=1,igridstail; igrid=igrids(iigrid);
        block=>ps(igrid)
        if(.not.block%is_physical_boundary(nb)) cycle
        if(stagger_grid) then
          !wbn(ixbcn^%1ixG^T)=block%ws(ixbcn^%1ixG^T,idir)
          wbn(ixbcn^%1ixG^T)=block%ws(ixbcn^%1ixG^T,idir)
        else
          wbn(ixbcn^%1ixG^T)=half*(block%w(ixbcn^%1ixG^T,iw_mag(idir))+block%w(ixbcn+1^%1ixG^T,iw_mag(idir)))
        end if
        xmina=block%x(1,1,1,2)-0.5d0*rnode(rpdx2_,igrid)
        xmaxa=block%x(1,ixGhi2,1,2)+0.5d0*rnode(rpdx2_,igrid)
        xminb=block%x(1,1,1,3)-0.5d0*rnode(rpdx3_,igrid)
        xmaxb=block%x(1,1,ixGhi3,3)+0.5d0*rnode(rpdx3_,igrid)
        if(xmina<rr(1,1,2) .and. xmaxa>rr(nc,1,2) .and. &
           xminb<rr(1,1,3) .and. xmaxb>rr(1,nc,3)) then
          do ix2=1,nc
            do ix1=1,nc
              ixbca=ceiling((rr(ix1,ix2,2)-xmina)/rnode(rpdx2_,igrid))
              ixbcb=ceiling((rr(ix1,ix2,3)-xminb)/rnode(rpdx3_,igrid))
              bc(ix1,ix2)=wbn(ixbcn,ixbca,ixbcb)
            end do
          end do
        else if(block%x(1,ixMlo2,1,2)>rmina .and. block%x(1,ixMhi2,1,2)<rmaxa .and. &
                block%x(1,1,ixMlo3,3)>rminb .and. block%x(1,1,ixMhi3,3)<rmaxb) then
          dlvl=node(plevel_,igrid)-box%lvl
          wnc=nc/2**dlvl
          allocate(xcoarse(wnc,wnc,2))
          do ix2=1,wnc
            do ix1=1,wnc
              xcoarse(ix1,ix2,1)=sum(block%x(1,(ix1-1)*2**dlvl+1+nghostcells:ix1*2**dlvl+nghostcells,1,2))/dble(2**dlvl)
              xcoarse(ix1,ix2,2)=sum(block%x(1,1,(ix2-1)*2**dlvl+1+nghostcells:ix2*2**dlvl+nghostcells,3))/dble(2**dlvl)
              ixbca=ceiling((xcoarse(ix1,ix2,1)-rmina)/box%dr(2))
              ixbcb=ceiling((xcoarse(ix1,ix2,2)-rminb)/box%dr(3))
              bc(ixbca,ixbcb)=sum(wbn(ixbcn,(ix1-1)*2**dlvl+1+nghostcells:ix1*2**dlvl+nghostcells,&
                                  (ix2-1)*2**dlvl+1+nghostcells:ix2*2**dlvl+nghostcells))/dble(2**(2*dlvl))
            end do
          end do
          deallocate(xcoarse)
        end if
      end do
    case(2)
      if(mod(nb,2)==0) then
        ixbcn=ixMhi2
      else
        ixbcn=ixMlo2-1
      end if
      rmina=rr(1,1,1)-0.5d0*box%dr(1)
      rmaxa=rr(nc,1,1)+0.5d0*box%dr(1)
      rminb=rr(1,1,3)-0.5d0*box%dr(3)
      rmaxb=rr(1,nc,3)+0.5d0*box%dr(3)
      do iigrid=1,igridstail; igrid=igrids(iigrid);
        block=>ps(igrid)
        if(.not.block%is_physical_boundary(nb)) cycle
        if(stagger_grid) then
          wbn(ixbcn^%2ixG^T)=block%ws(ixbcn^%2ixG^T,idir)
        else
          wbn(ixbcn^%2ixG^T)=half*(block%w(ixbcn^%2ixG^T,iw_mag(idir))+block%w(ixbcn+1^%2ixG^T,iw_mag(idir)))
        end if
        xmina=block%x(1,1,1,1)-0.5d0*rnode(rpdx1_,igrid)
        xmaxa=block%x(ixGhi1,1,1,1)+0.5d0*rnode(rpdx1_,igrid)
        xminb=block%x(1,1,1,3)-0.5d0*rnode(rpdx3_,igrid)
        xmaxb=block%x(1,1,ixGhi3,3)+0.5d0*rnode(rpdx3_,igrid)
        if(xmina<rr(1,1,1) .and. xmaxa>rr(nc,1,1) .and. &
           xminb<rr(1,1,3) .and. xmaxb>rr(1,nc,3)) then
          do ix2=1,nc
            do ix1=1,nc
              ixbca=ceiling((rr(ix1,ix2,1)-xmina)/rnode(rpdx1_,igrid))
              ixbcb=ceiling((rr(ix1,ix2,3)-xminb)/rnode(rpdx3_,igrid))
              bc(ix1,ix2)=wbn(ixbca,ixbcn,ixbcb)
            end do
          end do
        else if(block%x(ixMlo1,1,1,1)>rmina .and. block%x(ixMhi1,1,1,1)<rmaxa .and. &
                block%x(1,1,ixMlo3,3)>rminb .and. block%x(1,1,ixMhi3,3)<rmaxb) then
          dlvl=node(plevel_,igrid)-box%lvl
          wnc=nc/2**dlvl
          allocate(xcoarse(wnc,wnc,2))
          do ix2=1,wnc
            do ix1=1,wnc
              xcoarse(ix1,ix2,1)=sum(block%x((ix1-1)*2**dlvl+1+nghostcells:ix1*2**dlvl+nghostcells,1,1,1))/dble(2**dlvl)
              xcoarse(ix1,ix2,2)=sum(block%x(1,1,(ix2-1)*2**dlvl+1+nghostcells:ix2*2**dlvl+nghostcells,3))/dble(2**dlvl)
              ixbca=ceiling((xcoarse(ix1,ix2,1)-rmina)/box%dr(1))
              ixbcb=ceiling((xcoarse(ix1,ix2,2)-rminb)/box%dr(3))
              bc(ixbca,ixbcb)=sum(wbn((ix1-1)*2**dlvl+1+nghostcells:ix1*2**dlvl+nghostcells,ixbcn,&
                                  (ix2-1)*2**dlvl+1+nghostcells:ix2*2**dlvl+nghostcells))/dble(2**(2*dlvl))
            end do
          end do
          deallocate(xcoarse)
        end if
      end do
    case(3)
      if(mod(nb,2)==0) then
        ixbcn=ixMhi3
      else
        ixbcn=ixMlo3-1
      end if
      rmina=rr(1,1,1)-0.5d0*box%dr(1)
      rmaxa=rr(nc,1,1)+0.5d0*box%dr(1)
      rminb=rr(1,1,2)-0.5d0*box%dr(2)
      rmaxb=rr(1,nc,2)+0.5d0*box%dr(2)
      do iigrid=1,igridstail; igrid=igrids(iigrid);
        block=>ps(igrid)
        if(.not.block%is_physical_boundary(nb)) cycle
        if(stagger_grid) then
          wbn(ixbcn^%3ixG^T)=block%ws(ixbcn^%3ixG^T,idir)
        else
          wbn(ixbcn^%3ixG^T)=half*(block%w(ixbcn^%3ixG^T,iw_mag(idir))+block%w(ixbcn+1^%3ixG^T,iw_mag(idir)))
        end if
        xmina=block%x(1,1,1,1)-0.5d0*rnode(rpdx1_,igrid)
        xmaxa=block%x(ixGhi1,1,1,1)+0.5d0*rnode(rpdx1_,igrid)
        xminb=block%x(1,1,1,2)-0.5d0*rnode(rpdx2_,igrid)
        xmaxb=block%x(1,ixGhi2,1,2)+0.5d0*rnode(rpdx2_,igrid)
        if(xmina<rr(1,1,1) .and. xmaxa>rr(nc,1,1) .and. &
           xminb<rr(1,1,2) .and. xmaxb>rr(1,nc,2)) then
          do ix2=1,nc
            do ix1=1,nc
              ixbca=ceiling((rr(ix1,ix2,1)-xmina)/rnode(rpdx1_,igrid))
              ixbcb=ceiling((rr(ix1,ix2,2)-xminb)/rnode(rpdx2_,igrid))
              bc(ix1,ix2)=wbn(ixbca,ixbcb,ixbcn)
            end do
          end do
        else if(block%x(ixMlo1,1,1,1)>rmina .and. block%x(ixMhi1,1,1,1)<rmaxa .and. &
                block%x(1,ixMlo2,1,2)>rminb .and. block%x(1,ixMhi2,1,2)<rmaxb) then
          dlvl=node(plevel_,igrid)-box%lvl
          wnc=nc/2**dlvl
          allocate(xcoarse(wnc,wnc,2))
          do ix2=1,wnc
            do ix1=1,wnc
              xcoarse(ix1,ix2,1)=sum(block%x((ix1-1)*2**dlvl+1+nghostcells:ix1*2**dlvl+nghostcells,1,1,1))/dble(2**dlvl)
              xcoarse(ix1,ix2,2)=sum(block%x(1,(ix2-1)*2**dlvl+1+nghostcells:ix2*2**dlvl+nghostcells,1,2))/dble(2**dlvl)
              ixbca=ceiling((xcoarse(ix1,ix2,1)-rmina)/box%dr(1))
              ixbcb=ceiling((xcoarse(ix1,ix2,2)-rminb)/box%dr(2))
              bc(ixbca,ixbcb)=sum(wbn((ix1-1)*2**dlvl+1+nghostcells:ix1*2**dlvl+nghostcells,&
                                  (ix2-1)*2**dlvl+1+nghostcells:ix2*2**dlvl+nghostcells,ixbcn))/dble(2**(2*dlvl))
            end do
          end do
          deallocate(xcoarse)
        end if
      end do
    end select

  end subroutine multigrid_bc
}
end module mod_lfff
