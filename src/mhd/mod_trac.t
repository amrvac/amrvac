module mod_trac
  use mod_global_parameters
  use mod_mhd
  implicit none
  ! common
  integer :: numFL,numLP
  double precision :: dL,Tmax,trac_delta,T_bott
  double precision, allocatable :: xFi(:,:)
  ! TRACB
  integer :: numxT^D
  double precision :: dxT^D 
  double precision :: xTmin(ndim),xTmax(ndim)
  double precision, allocatable :: xT(:^D&,:)
  ! TRACL mask
  integer, allocatable :: trac_grid(:),ground_grid(:)
  integer :: ngrid_trac,ngrid_ground
  logical, allocatable :: trac_pe(:)
contains

  subroutine init_trac_line(mask)
    logical, intent(in) :: mask 
    integer :: refine_factor,ix^D,ix(ndim),j,iFL,numL(ndim),finegrid
    double precision :: lengthFL
    double precision :: xprobmin(ndim),xprobmax(ndim),domain_nx(ndim)
    integer :: memxFi

    trac_delta=0.25d0
    !----------------- init magnetic field -------------------!
    refine_factor=2**(refine_max_level-1)
    ^D&xprobmin(^D)=xprobmin^D\
    ^D&xprobmax(^D)=xprobmax^D\
    ^D&domain_nx(^D)=domain_nx^D\
    dL=(xprobmax(ndim)-xprobmin(ndim))/(domain_nx(ndim)*refine_factor)
    ! max length of a field line
    if (.not. mask) phys_trac_mask=xprobmax^ND
    {^IFTWOD
      lengthFL=(phys_trac_mask-xprobmin(2))*3.d0
    }
    {^IFTHREED
      lengthFL=(phys_trac_mask-xprobmin(3))*3.d0
    }
    
    numLP=floor(lengthFL/dL)
    numL=1
    numFL=1
    do j=1,ndim-1
      ! number of field lines, every 4 finest grid, every direction
      finegrid=mhd_trac_finegrid
      numL(j)=floor((xprobmax(j)-xprobmin(j))/dL/finegrid)
      numFL=numFL*numL(j)
    end do
    allocate(xFi(numFL,ndim))
    xFi(:,ndim)=xprobmin(ndim)+dL/50.d0
    {do ix^DB=1,numL(^DB)\}
      ^D&ix(^D)=ix^D\
      iFL=0
      do j=ndim-1,1,-1
        iFL=iFL+(ix(j)-(ndim-1-j))*(numL(j))**(ndim-1-j)
      end do
      xFi(iFL,1:ndim-1)=xprobmin(1:ndim-1)+finegrid*dL*ix(1:ndim-1)-finegrid*dL/2.d0
    {end do\}
    {^IFTWOD
    if(mype .eq. 0) write(*,*) 'NOTE:  2D TRAC method take the y-dir == grav-dir'
    }
    {^IFTHREED
    if(mype .eq. 0) write(*,*) 'NOTE:  3D TRAC method take the z-dir == grav-dir'
    }

    memxFi=floor(8*numFL*numLP*ndim/1e6)
    if (mype==0) write(*,*) 'Memory requirement for each processor in TRAC:'
    if (mype==0) write(*,*)  memxFi,' MB'

    allocate(trac_pe(0:npe-1),trac_grid(max_blocks),ground_grid(max_blocks))

  end subroutine init_trac_line

  subroutine init_trac_block(mask)
    logical, intent(in) :: mask 
    integer :: refine_factor,finegrid,iFL,j
    integer :: ix^D,ixT^L
    integer :: numL(ndim),ix(ndim)
    double precision :: lengthFL
    double precision :: ration,a0
    double precision :: xprobmin(ndim),xprobmax(ndim),dxT(ndim)

    refine_factor=2**(refine_max_level-1)
    ^D&xprobmin(^D)=xprobmin^D\
    ^D&xprobmax(^D)=xprobmax^D\
    ^D&dxT^D=(xprobmax^D-xprobmin^D)/(domain_nx^D*refine_factor/block_nx^D)\
    ^D&dxT(ndim)=dxT^D\
    finegrid=mhd_trac_finegrid
    {^IFONED
      dL=dxT1/finegrid
    }
    {^NOONED
      dL=min(dxT^D)/finegrid
    }
    ! table for interpolation
    ^D&xTmin(^D)=xprobmin^D\
    ^D&xTmax(^D)=xprobmax^D\
    if(mask) xTmax(ndim)=phys_trac_mask
    ! max length of a field line
    if(mask) then
      lengthFL=maxval(xTmax-xprobmin)*3.d0
    else
      lengthFL=maxval(xprobmax-xprobmin)*3.d0
    end if
    numLP=floor(lengthFL/dL)
    ^D&numxT^D=ceiling((xTmax(^D)-xTmin(^D)-smalldouble)/dxT^D)\
    allocate(xT(numxT^D,ndim))
    ixTmin^D=1\
    ixTmax^D=numxT^D\
    {do j=1,numxT^D
      xT(j^D%ixT^S,^D)=(j-0.5d0)*dxT^D+xTmin(^D)
    end do\}
    if(mask) xTmax(ndim)=maxval(xT(:^D&,ndim))+half*dxT(ndim)
    numL=1
    numFL=1
    do j=1,ndim-1
      ! number of field lines, every 4 finest grid, every direction
      numL(j)=floor((xprobmax(j)-xprobmin(j))/dL)
      numFL=numFL*numL(j)
    end do
    allocate(xFi(numFL,ndim))
    xFi(:,ndim)=xprobmin(ndim)+dL/50.d0
    {do ix^DB=1,numL(^DB)\}
      ^D&ix(^D)=ix^D\
      iFL=0
      do j=ndim-1,1,-1
        iFL=iFL+(ix(j)-(ndim-1-j))*(numL(j))**(ndim-1-j)
      end do
      xFi(iFL,1:ndim-1)=xprobmin(1:ndim-1)+dL*ix(1:ndim-1)-dL/2.d0
    {end do\}
    {^IFTWOD
    if(mype .eq. 0) write(*,*) 'NOTE:  2D TRAC method take the y-dir == grav-dir'
    }
    {^IFTHREED
    if(mype .eq. 0) write(*,*) 'NOTE:  3D TRAC method take the z-dir == grav-dir'
    }
  end subroutine init_trac_block

  subroutine TRAC_simple(tco_global,trac_alfa,T_peak)
    double precision, intent(in) :: tco_global, trac_alfa,T_peak
    integer :: iigrid, igrid

    do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
      {^IFONED
      ps(igrid)%special_values(1)=tco_global
      }
      if(ps(igrid)%special_values(1)<trac_alfa*ps(igrid)%special_values(2)) then
        ps(igrid)%special_values(1)=trac_alfa*ps(igrid)%special_values(2)
      end if
      if(ps(igrid)%special_values(1) .lt. T_bott) then
        ps(igrid)%special_values(1)=T_bott
      else if(ps(igrid)%special_values(1) .gt. 0.2d0*T_peak) then
        ps(igrid)%special_values(1)=0.2d0*T_peak
      end if
      ps(igrid)%wextra(ixG^T,iw_tcoff)=ps(igrid)%special_values(1)
      !> special values(2) to save old tcutoff
      ps(igrid)%special_values(2)=ps(igrid)%special_values(1)
    end do
  end subroutine TRAC_simple

  subroutine LTRAC(T_peak)
    double precision, intent(in) :: T_peak
    integer :: iigrid, igrid
    integer :: ixO^L,trac_tcoff

    ixO^L=ixM^LL;
    trac_tcoff=iw_tcoff
    do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
      where(ps(igrid)%wextra(ixO^S,trac_tcoff) .lt. T_bott)
        ps(igrid)%wextra(ixO^S,trac_tcoff)=T_bott
      else where(ps(igrid)%wextra(ixO^S,trac_tcoff) .gt. 0.2d0*T_peak)
        ps(igrid)%wextra(ixO^S,trac_tcoff)=0.2d0*T_peak
      end where
    end do
  end subroutine LTRAC

  subroutine TRACL(mask,T_peak)
    logical, intent(in) :: mask
    double precision, intent(in) :: T_peak
    integer :: ix^D,j
    integer :: iigrid, igrid
    double precision :: xF(numFL,numLP,ndim)
    integer :: numR(numFL),ix^L
    double precision :: Tlcoff(numFL)
    integer :: ipel(numFL,numLP),igridl(numFL,numLP)
    logical :: forwardl(numFL)

    Tmax=T_peak
    xF=zero
    
    do j=1,ndim
      xF(:,1,j)=xFi(:,j)
    end do

    call MPI_BARRIER(icomm,ierrmpi)
    ! record pe and grid for mask region
    call update_pegrid()
    ! get temperature for the calculation of Te gradient, 
    ! which will be stored in wextra(ixI^S,Tcoff_) temporarily
    call get_Te_grid()
    ! find out direction for tracing B field
    call get_Btracing_dir(ipel,igridl,forwardl)
    ! trace Bfield and get Tcoff for each field line
    call get_Tcoff_line(xF,numR,Tlcoff,ipel,igridl,forwardl,mask)
    ! init vairable Tcoff_ and Tweight_
    call init_trac_Tcoff()
    ! get cell Tcoff via interpolation
    call interp_Tcoff(xF,ipel,igridl,numR,Tlcoff)
    ! convert primitive data back to conserved data
    call MPI_BARRIER(icomm,ierrmpi)


  end subroutine TRACL

  subroutine TRACB(mask,T_peak)
    logical, intent(in) :: mask
    double precision, intent(in) :: T_peak
    integer :: peArr(numxT^D),gdArr(numxT^D),numR(numFL)
    double precision :: Tcoff(numxT^D),Tcmax(numxT^D),Bdir(numxT^D,ndim)
    double precision :: xF(numFL,numLP,ndim),Tcoff_line(numFL)
    integer :: xpe(numFL,numLP,2**ndim)
    integer :: xgd(numFL,numLP,2**ndim)

    Tmax=T_peak
    Tcoff=zero
    Tcmax=zero
    Bdir=zero
    peArr=-1
    gdArr=-1
    call block_estable(mask,Tcoff,Tcmax,Bdir,peArr,gdArr)
    xF=zero
    numR=0
    Tcoff_line=zero
    call block_trace_mfl(mask,Tcoff,Tcoff_line,Tcmax,Bdir,peArr,gdArr,xF,numR,xpe,xgd)
    call block_interp_grid(mask,xF,numR,xpe,xgd,Tcoff_line)
  end subroutine TRACB

  subroutine block_estable(mask,Tcoff,Tcmax,Bdir,peArr,gdArr)
    logical :: mask
    double precision :: Tcoff(numxT^D),Tcoff_recv(numxT^D)
    double precision :: Tcmax(numxT^D),Tcmax_recv(numxT^D)
    double precision :: Bdir(numxT^D,ndim),Bdir_recv(numxT^D,ndim)
    integer :: peArr(numxT^D),peArr_recv(numxT^D)
    integer :: gdArr(numxT^D),gdArr_recv(numxT^D)
    integer :: xc^L,xd^L,ix^D
    integer :: iigrid,igrid,numxT,intab
    double precision :: xb^L

    Tcoff_recv=zero
    Tcmax_recv=zero
    Bdir_recv=zero
    peArr_recv=-1
    gdArr_recv=-1
    !> combine table from different processors 
    xcmin^D=nghostcells+1\
    xcmax^D=block_nx^D+nghostcells\
    do iigrid=1,igridstail; igrid=igrids(iigrid);
      ps(igrid)%wextra(:^D&,Tweight_)=zero
      ps(igrid)%wextra(:^D&,Tcoff_)=zero
      ^D&xbmin^D=rnode(rpxmin^D_,igrid)-xTmin(^D)\
      ^D&xbmax^D=rnode(rpxmax^D_,igrid)-xTmin(^D)\
      xdmin^D=nint(xbmin^D/dxT^D)+1\
      xdmax^D=ceiling((xbmax^D-smalldouble)/dxT^D)\
      {do ix^D=xdmin^D,xdmax^D \}
        intab=0
        {if (ix^D .le. numxT^D) intab=intab+1 \}
        if(intab .eq. ndim) then
          !> in principle, no overlap will happen here
          Tcoff(ix^D)=max(Tcoff(ix^D),ps(igrid)%special_values(1))
          Tcmax(ix^D)=ps(igrid)%special_values(2)
          !> abs(Bdir) <= 1, so that Bdir+2 should always be positive
          Bdir(ix^D,1:ndim)=ps(igrid)%special_values(3:3+ndim-1)+2.d0
          peArr(ix^D)=mype
          gdArr(ix^D)=igrid
        end if
      {end do\}
    end do 
    call MPI_BARRIER(icomm,ierrmpi)
    numxT={numxT^D*}
    call MPI_ALLREDUCE(peArr,peArr_recv,numxT,MPI_INTEGER,&
                       MPI_MAX,icomm,ierrmpi)
    call MPI_ALLREDUCE(gdArr,gdArr_recv,numxT,MPI_INTEGER,&
                       MPI_MAX,icomm,ierrmpi)
    call MPI_ALLREDUCE(Tcoff,Tcoff_recv,numxT,MPI_DOUBLE_PRECISION,&
                       MPI_MAX,icomm,ierrmpi)
    call MPI_ALLREDUCE(Bdir,Bdir_recv,numxT*ndim,MPI_DOUBLE_PRECISION,&
                       MPI_MAX,icomm,ierrmpi)
    if(.not. mask) then
      call MPI_ALLREDUCE(Tcmax,Tcmax_recv,numxT,MPI_DOUBLE_PRECISION,&
                         MPI_MAX,icomm,ierrmpi)
    end if
    peArr=peArr_recv
    gdArr=gdArr_recv
    Tcoff=Tcoff_recv
    Bdir=Bdir_recv-2.d0
    if(.not. mask) Tcmax=Tcmax_recv
  end subroutine block_estable

  subroutine block_trace_mfl(mask,Tcoff,Tcoff_line,Tcmax,Bdir,peArr,gdArr,xF,numR,xpe,xgd)
    integer :: i,j,k,k^D,ix_next^D
    logical :: mask,flag,first
    double precision :: Tcoff(numxT^D),Tcoff_line(numFL)
    double precision :: Tcmax(numxT^D),Tcmax_line(numFL)
    double precision :: xF(numFL,numLP,ndim)
    integer :: ix_mod(ndim,2),numR(numFL)
    double precision :: alfa_mod(ndim,2)
    double precision :: nowpoint(ndim),nowgridc(ndim)
    double precision :: Bdir(numxT^D,ndim)
    double precision :: init_dir,now_dir1(ndim),now_dir2(ndim)
    integer :: peArr(numxT^D),xpe(numFL,numLP,2**ndim)
    integer :: gdArr(numxT^D),xgd(numFL,numLP,2**ndim)

    do i=1,numFL
      flag=.true.
      ^D&k^D=ceiling((xFi(i,^D)-xTmin(^D)-smalldouble)/dxT^D)\
      Tcoff_line(i)=Tcoff(k^D)
      if(.not. mask) Tcmax_line(i)=Tcmax(k^D)
      ix_next^D=k^D\
      j=1
      xF(i,j,:)=xFi(i,:)
      do while(flag)
        nowpoint(:)=xF(i,j,:)
        nowgridc(:)=xT(ix_next^D,:)
        first=.true.
        if(j .eq. 1) then 
          call RK_Bdir(nowgridc,nowpoint,ix_next^D,now_dir1,Bdir,&
                       ix_mod,first,init_dir)
        else
          call RK_Bdir(nowgridc,nowpoint,ix_next^D,now_dir1,Bdir,&
                       ix_mod,first)
        end if
        {^IFTWOD
          xgd(i,j,1)=gdArr(ix_mod(1,1),ix_mod(2,1))
          xgd(i,j,2)=gdArr(ix_mod(1,2),ix_mod(2,1))
          xgd(i,j,3)=gdArr(ix_mod(1,1),ix_mod(2,2))
          xgd(i,j,4)=gdArr(ix_mod(1,2),ix_mod(2,2))
          xpe(i,j,1)=peArr(ix_mod(1,1),ix_mod(2,1))
          xpe(i,j,2)=peArr(ix_mod(1,2),ix_mod(2,1))
          xpe(i,j,3)=peArr(ix_mod(1,1),ix_mod(2,2))
          xpe(i,j,4)=peArr(ix_mod(1,2),ix_mod(2,2))
        }
        {^IFTHREED
          xgd(i,j,1)=gdArr(ix_mod(1,1),ix_mod(2,1),ix_mod(3,1))
          xgd(i,j,2)=gdArr(ix_mod(1,2),ix_mod(2,1),ix_mod(3,1))
          xgd(i,j,3)=gdArr(ix_mod(1,1),ix_mod(2,2),ix_mod(3,1))
          xgd(i,j,4)=gdArr(ix_mod(1,2),ix_mod(2,2),ix_mod(3,1))
          xgd(i,j,5)=gdArr(ix_mod(1,1),ix_mod(2,1),ix_mod(3,2))
          xgd(i,j,6)=gdArr(ix_mod(1,2),ix_mod(2,1),ix_mod(3,2))
          xgd(i,j,7)=gdArr(ix_mod(1,1),ix_mod(2,2),ix_mod(3,2))
          xgd(i,j,8)=gdArr(ix_mod(1,2),ix_mod(2,2),ix_mod(3,2))
          xpe(i,j,1)=peArr(ix_mod(1,1),ix_mod(2,1),ix_mod(3,1))
          xpe(i,j,2)=peArr(ix_mod(1,2),ix_mod(2,1),ix_mod(3,1))
          xpe(i,j,3)=peArr(ix_mod(1,1),ix_mod(2,2),ix_mod(3,1))
          xpe(i,j,4)=peArr(ix_mod(1,2),ix_mod(2,2),ix_mod(3,1))
          xpe(i,j,5)=peArr(ix_mod(1,1),ix_mod(2,1),ix_mod(3,2))
          xpe(i,j,6)=peArr(ix_mod(1,2),ix_mod(2,1),ix_mod(3,2))
          xpe(i,j,7)=peArr(ix_mod(1,1),ix_mod(2,2),ix_mod(3,2))
          xpe(i,j,8)=peArr(ix_mod(1,2),ix_mod(2,2),ix_mod(3,2))
        }
        nowpoint(:)=nowpoint(:)+init_dir*now_dir1*dL
        {if(nowpoint(^D) .gt. xTmax(^D) .or. nowpoint(^D) .lt. xTmin(^D)) then
          flag=.false.
        end if\}
        if(mask .and. nowpoint(ndim) .gt. phys_trac_mask) then
          flag=.false.
        end if
        if(flag) then
          first=.false.
          ^D&ix_next^D=ceiling((nowpoint(^D)-xTmin(^D)-smalldouble)/dxT^D)\
          nowgridc(:)=xT(ix_next^D,:)
          call RK_Bdir(nowgridc,nowpoint,ix_next^D,now_dir2,Bdir,&
                       ix_mod,first)
          xF(i,j+1,:)=xF(i,j,:)+init_dir*dL*half*(now_dir1+now_dir2)
          {if(xF(i,j+1,^D) .gt. xTmax(^D) .or. xF(i,j+1,^D) .lt. xTmin(^D)) then
            flag=.false.
          end if\}
          if(mask .and. xF(i,j+1,ndim) .gt. phys_trac_mask) then
            flag=.false.
          end if
          if(flag) then
            ^D&ix_next^D=ceiling((xF(i,j+1,^D)-xTmin(^D)-smalldouble)/dxT^D)\
            j=j+1
            Tcoff_line(i)=max(Tcoff_line(i),Tcoff(ix_next^D))
            if(.not.mask) Tcmax_line(i)=max(Tcmax_line(i),Tcmax(ix_next^D))
          end if
        end if
      end do
      numR(i)=j
      if(mask) then
        if(Tcoff_line(i) .gt. Tmax*0.2d0) then
          Tcoff_line(i)=Tmax*0.2d0
        end if
      else
        if(Tcoff_line(i) .gt. Tcmax_line(i)*0.2d0) then
          Tcoff_line(i)=Tcmax_line(i)*0.2d0
        end if
      end if
    end do
  end subroutine block_trace_mfl

  subroutine RK_Bdir(nowgridc,nowpoint,ix_next^D,now_dir,Bdir,ix_mod,first,init_dir)
    double precision :: nowpoint(ndim),nowgridc(ndim)
    integer :: ix_mod(ndim,2)
    double precision :: alfa_mod(ndim,2)
    integer :: ix_next^D,k^D
    double precision :: now_dir(ndim)
    double precision :: Bdir(numxT^D,ndim)
    logical :: first
    double precision, optional :: init_dir

    {if(nowpoint(^D) .gt. xTmin(^D)+half*dxT^D .and. nowpoint(^D) .lt. xTmax(^D)-half*dxT^D) then
      if(nowpoint(^D) .le. nowgridc(^D)) then
        ix_mod(^D,1)=ix_next^D-1
        ix_mod(^D,2)=ix_next^D
        alfa_mod(^D,1)=abs(nowgridc(^D)-nowpoint(^D))/dxT^D
        alfa_mod(^D,2)=one-alfa_mod(^D,1)
      else
        ix_mod(^D,1)=ix_next^D
        ix_mod(^D,2)=ix_next^D+1
        alfa_mod(^D,2)=abs(nowgridc(^D)-nowpoint(^D))/dxT^D
        alfa_mod(^D,1)=one-alfa_mod(^D,2)
      end if
    else
      ix_mod(^D,:)=ix_next^D
      alfa_mod(^D,:)=half
    end if\}
    now_dir=zero
    {^IFTWOD
    do k1=1,2
    do k2=1,2
      now_dir=now_dir + Bdir(ix_mod(1,k1),ix_mod(2,k2),:)*alfa_mod(1,k1)*alfa_mod(2,k2)
    end do
    end do
    } 
    {^IFTHREED
    do k1=1,2
    do k2=1,2
    do k3=1,2
      now_dir=now_dir + Bdir(ix_mod(1,k1),ix_mod(2,k2),ix_mod(3,k3),:)&
             *alfa_mod(1,k1)*alfa_mod(2,k2)*alfa_mod(3,k3)
    end do
    end do
    end do
    }
    if(present(init_dir)) then
      init_dir=sign(one,now_dir(ndim))
    end if
  end subroutine RK_Bdir

  subroutine block_interp_grid(mask,xF,numR,xpe,xgd,Tcoff_line)
    logical :: mask
    double precision :: xF(numFL,numLP,ndim)
    integer :: numR(numFL)
    integer :: xpe(numFL,numLP,2**ndim)
    integer :: xgd(numFL,numLP,2**ndim)
    double precision :: Tcoff_line(numFL)
    double precision :: weightIndex,weight,ds
    integer :: i,j,k,igrid,iigrid,ixO^L,ixc^L,ixc^D
    double precision :: dxMax^D,dxb^D

    !> interpolate lines into grid cells
    weightIndex=2.d0
    dxMax^D=2.d0*dL;
    ixO^L=ixM^LL;
    do i=1,numFL
      do j=1,numR(i) 
        do k=1,2**ndim
          if(mype .eq. xpe(i,j,k)) then
            igrid=xgd(i,j,k)
            if(igrid .le. igrids(igridstail)) then
              ^D&dxb^D=rnode(rpdx^D_,igrid)\
              ^D&ixcmin^D=floor((xF(i,j,^D)-dxMax^D-ps(igrid)%x(ixOmin^DD,^D))/dxb^D)+ixOmin^D\
              ^D&ixcmax^D=floor((xF(i,j,^D)+dxMax^D-ps(igrid)%x(ixOmin^DD,^D))/dxb^D)+ixOmin^D\
              {if (ixcmin^D<ixOmin^D) ixcmin^D=ixOmin^D\}
              {if (ixcmax^D>ixOmax^D) ixcmax^D=ixOmax^D\}
              {do ixc^D=ixcmin^D,ixcmax^D\}
                ds=0.d0
                {ds=ds+(xF(i,j,^D)-ps(igrid)%x(ixc^DD,^D))**2\}
                ds=sqrt(ds)
                if(ds .le. 0.099d0*dL) then
                  weight=(1/(0.099d0*dL))**weightIndex
                else
                  weight=(1/ds)**weightIndex
                endif
                ps(igrid)%wextra(ixc^D,Tweight_)=ps(igrid)%wextra(ixc^D,Tweight_)+weight
                ps(igrid)%wextra(ixc^D,Tcoff_)=ps(igrid)%wextra(ixc^D,Tcoff_)+weight*Tcoff_line(i)
              {enddo\}
            else
              call mpistop("we need to check here 366Line in mod_trac.t")
            end if
          end if
        end do
      end do
    end do
    ! finish interpolation
    do iigrid=1,igridstail; igrid=igrids(iigrid);
      where (ps(igrid)%wextra(ixO^S,Tweight_)>0.d0)
        ps(igrid)%wextra(ixO^S,Tcoff_)=ps(igrid)%wextra(ixO^S,Tcoff_)/ps(igrid)%wextra(ixO^S,Tweight_)
      elsewhere
        ps(igrid)%wextra(ixO^S,Tcoff_)=0.2d0*Tmax
      endwhere
    enddo
  end subroutine block_interp_grid

   ! TODO remove, not used? 
!  subroutine get_Tmax_grid(x,w,ixI^L,ixO^L,Tmax_grid)
!    integer, intent(in) :: ixI^L,ixO^L
!    double precision, intent(in) :: x(ixI^S,1:ndim)
!    double precision, intent(out) :: w(ixI^S,1:nw)
!    double precision :: Tmax_grid
!    double precision :: pth(ixI^S),Te(ixI^S)
!
!    call mhd_get_pthermal(w,x,ixI^L,ixO^L,pth)
!    Te(ixO^S)=pth(ixO^S)/w(ixO^S,rho_)
!    Tmax_grid=maxval(Te(ixO^S))
!  end subroutine get_Tmax_grid  

  subroutine init_trac_Tcoff()
    integer :: ixI^L,ixO^L,igrid,iigrid

    ixI^L=ixG^LL;
    ixO^L=ixM^LL;

    do iigrid=1,igridstail; igrid=igrids(iigrid);
      ps(igrid)%wextra(ixI^S,Tcoff_)=0.d0
      ps(igrid)%wextra(ixI^S,Tweight_)=0.d0
    end do

  end subroutine init_trac_Tcoff

  subroutine update_pegrid()

    ngrid_trac=0
    ngrid_ground=0
    trac_pe(:)=.false.
    trac_grid(:)=0
    ground_grid(:)=0

    ! traverse gridtable to information of grid inside mask region
    call traverse_gridtable()

  end subroutine update_pegrid

  subroutine traverse_gridtable()
    use mod_global_parameters

    double precision :: dxb^D,xb^L
    integer :: iigrid,igrid,j
    logical, allocatable :: trac_pe_recv(:)
    double precision :: hcmax_bt

    allocate(trac_pe_recv(npe))

    hcmax_bt=xprobmin^ND+(xprobmax^ND-xprobmin^ND)/(domain_nx^ND*2**(refine_max_level-1))

    do iigrid=1,igridstail; igrid=igrids(iigrid);
      xbmin^ND=rnode(rpxmin^ND_,igrid)
      if (xbmin^ND<phys_trac_mask) then
        trac_pe(mype)=.true.
        ngrid_trac=ngrid_trac+1
        trac_grid(ngrid_trac)=igrid
        if (xbmin^ND<hcmax_bt) then
          ngrid_ground=ngrid_ground+1
          ground_grid(ngrid_ground)=igrid
        endif
      endif
    enddo

    call MPI_ALLREDUCE(trac_pe,trac_pe_recv,npe,MPI_LOGICAL,MPI_LOR,icomm,ierrmpi)
    trac_pe=trac_pe_recv

    deallocate(trac_pe_recv)

  end subroutine traverse_gridtable

  subroutine get_Te_grid()
    integer :: ixI^L,ixO^L,igrid,iigrid,j

    ixI^L=ixG^LL;
    ixO^L=ixM^LL;

    do j=1,ngrid_trac
      igrid=trac_grid(j)
      
      call mhd_get_pthermal(ps(igrid)%w,ps(igrid)%x,ixI^L,ixI^L,ps(igrid)%wextra(ixI^S,Tcoff_))
      !TODO move  outside loop for optimziation? 
      if(has_equi_rho0) then
        ps(igrid)%wextra(ixI^S,Tcoff_)=ps(igrid)%wextra(ixI^S,Tcoff_)/&
          (ps(igrid)%w(ixI^S,rho_) + ps(igrid)%equi_vars(ixI^S,equi_rho0_,0))
      else
        ps(igrid)%wextra(ixI^S,Tcoff_)=ps(igrid)%wextra(ixI^S,Tcoff_)/ps(igrid)%w(ixI^S,rho_)
      endif
    enddo

  end subroutine get_Te_grid

  subroutine get_Btracing_dir(ipel,igridl,forwardl)
    integer :: ipel(numFL,numLP),igridl(numFL,numLP)
    logical :: forwardl(numFL)

    integer :: igrid,ixO^L,iFL,j,ix^D,idir,ixb^D,ixbb^D
    double precision :: xb^L,dxb^D,xd^D,factor,Bh
    integer :: numL(ndim),ixmin(ndim),ixmax(ndim),ix(ndim)
    logical :: forwardRC(numFL)

    ipel=-1
    igridl=0
    forwardl=.TRUE.

    ^D&ixOmin^D=ixmlo^D;
    ^D&ixOmax^D=ixmhi^D;

    do j=1,ngrid_ground
      igrid=ground_grid(j)
      ^D&dxb^D=rnode(rpdx^D_,igrid);
      ^D&xbmin^D=rnode(rpxmin^D_,igrid);
      ^D&xbmax^D=rnode(rpxmax^D_,igrid);
      ^D&ixmin(^D)=floor((xbmin^D-xprobmin^D)/(mhd_trac_finegrid*dL))+1;
      ^D&ixmax(^D)=floor((xbmax^D-xprobmin^D)/(mhd_trac_finegrid*dL));
      ^D&numL(^D)=floor((xprobmax^D-xprobmin^D)/(mhd_trac_finegrid*dL));
      ixmin(^ND)=1
      ixmax(^ND)=1
      numL(^ND)=1

      {do ix^DB=ixmin(^DB),ixmax(^DB)\}
        ^D&ix(^D)=ix^D;
        iFL=0
        do idir=ndim-1,1,-1
          iFL=iFL+(ix(idir)-(ndim-1-idir))*(numL(idir))**(ndim-1-idir)
        enddo
        ipel(iFL,1)=mype
        igridl(iFL,1)=igrid

        ^D&ixb^D=floor((xFi(iFL,^D)-ps(igrid)%x(ixOmin^DD,^D))/dxb^D)+ixOmin^D;
        ^D&xd^D=(xFi(iFL,^D)-ps(igrid)%x(ixb^DD,^D))/dxb^D;
        Bh=0.d0
        {do ixbb^D=0,1\}
          factor={abs(1-ix^D-xd^D)*}
          if (B0field) then
            Bh=Bh+factor*(ps(igrid)%w(ixb^D+ixbb^D,mag(^ND))+ps(igrid)%B0(ixb^D+ixbb^D,^ND,0))
          else
            Bh=Bh+factor*ps(igrid)%w(ixb^D+ixbb^D,mag(^ND))
          endif
        {enddo\}
        if (Bh>0) then
          forwardl(iFL)=.TRUE.
        else
          forwardl(iFL)=.FALSE.
        endif
      {enddo\}
    enddo

    call MPI_ALLREDUCE(forwardl,forwardRC,numFL,MPI_LOGICAL,&
                         MPI_LAND,icomm,ierrmpi)
    forwardl=forwardRC

  end subroutine get_Btracing_dir

  subroutine get_Tcoff_line(xFL,numR,TcoffFL,ipeFL,igridFL,forwardFL,mask)
    use mod_trace_field

    double precision :: xFL(numFL,numLP,ndim)
    integer :: numR(numFL)
    double precision :: TcoffFL(numFL),TmaxFL(numFL)
    integer :: ipeFL(numFL,numLP),igridFL(numFL,numLP)
    logical :: forwardFL(numFL)
    logical, intent(in) :: mask

    integer :: nwP,nwL,iFL,iLP
    double precision :: wPm(numFL,numLP,2),wLm(numFL,1+2)
    character(len=std_len) :: ftype,tcondi

    nwP=2
    nwL=2
    ftype='Bfield'
    tcondi='TRAC'    
    call trace_field_multi(xFL,wPm,wLm,dL,numFL,numLP,nwP,nwL,forwardFL,ftype,tcondi)
    do iFL=1,numFL
      numR(iFL)=int(wLm(iFL,1))
      TcoffFL(iFL)=wLm(iFL,2)
      TmaxFL(iFL)=wLm(iFL,3)
      if(mask) then
        if(TcoffFL(iFL)>0.2d0*Tmax) TcoffFL(iFL)=0.2d0*Tmax
      else
        TmaxFL(iFL)=wLm(iFL,3)
        if(TcoffFL(iFL)>0.2d0*TmaxFL(iFL)) TcoffFL(iFL)=0.2d0*TmaxFL(iFL)
      end if

      if(TcoffFL(iFL)<T_bott) TcoffFL(iFL)=T_bott
    enddo

    do iFL=1,numFL
      if (numR(iFL)>0) then
        do iLP=1,numR(iFL)
          ipeFL(iFL,iLP)=int(wPm(iFL,iLP,1))
          igridFL(iFL,iLP)=int(wPm(iFL,iLP,2))
        enddo
      endif
    enddo


  end subroutine get_Tcoff_line

  subroutine interp_Tcoff(xF,ipel,igridl,numR,Tlcoff)
    double precision :: xF(numFL,numLP,ndim)
    integer :: numR(numFL),ipel(numFL,numLP),igridl(numFL,numLP)
    double precision :: Tlcoff(numFL)

    integer :: iFL,iLP,ixO^L,ixI^L,ixc^L,ixb^L,ixc^D
    integer :: igrid,j,ipmin,ipmax,igrid_nb
    double precision :: dxb^D,dxMax^D,xb^L,Tcnow
    double precision :: xFnow(ndim)
    integer :: weightIndex,idn^D,ixmax^ND
    double precision :: ds,weight

    weightIndex=2

    ^D&ixImin^D=ixglo^D;
    ^D&ixImax^D=ixghi^D;
    ^D&ixOmin^D=ixmlo^D;
    ^D&ixOmax^D=ixmhi^D;

    do iFL=1,numFL
      iLP=1
      Tcnow=Tlcoff(iFL)
      do while (iLP<=numR(iFL))      
        ! find points in the same grid
        ipmin=iLP
        do while (ipel(iFL,ipmin)/=mype .and. ipmin<=numR(iFL))
          ipmin=ipmin+1
        enddo
        igrid=igridl(iFL,ipmin)
        ipmax=ipmin
        do while (ipel(iFL,ipmax)==mype .and. igridl(iFL,ipmax+1)==igrid .and. ipmax<numR(iFL))
          ipmax=ipmax+1
        enddo

        ! parameters for the grid
        ^D&dxb^D=rnode(rpdx^D_,igrid);
        ^D&dxMax^D=4*dxb^D;

        do iLP=ipmin,ipmax
          xFnow(:)=xF(iFL,iLP,:)
          ^D&ixbmin^D=floor((xFnow(^D)-dxMax^D-ps(igrid)%x(ixOmin^DD,^D))/dxb^D)+ixOmin^D;
          ^D&ixbmax^D=floor((xFnow(^D)+dxMax^D-ps(igrid)%x(ixOmin^DD,^D))/dxb^D)+ixOmin^D;

          ! do interpolation inside the grid
          {ixcmin^D=max(ixbmin^D,ixOmin^D)\}
          {ixcmax^D=min(ixbmax^D,ixOmax^D)\}
          xbmin^ND=rnode(rpxmin^ND_,igrid)
          xbmax^ND=rnode(rpxmax^ND_,igrid)
          ixmax^ND=floor((phys_trac_mask-xbmin^ND)/dxb^ND)+ixOmin^ND
          if (xbmax^ND>phys_trac_mask) ixcmax^ND=min(ixmax^ND,ixcmax^ND)
          {do ixc^D=ixcmin^D,ixcmax^D\}
            ds=0.d0
            {ds=ds+(xFnow(^D)-ps(igrid)%x(ixc^DD,^D))**2\}
            ds=sqrt(ds)
            if(ds<1.0d-2*dxb1) then
              weight=(1/(1.0d-2*dxb1))**weightIndex
            else
              weight=(1/ds)**weightIndex
            endif
            ps(igrid)%wextra(ixc^D,Tweight_)=ps(igrid)%wextra(ixc^D,Tweight_)+weight
            ps(igrid)%wextra(ixc^D,Tcoff_)=ps(igrid)%wextra(ixc^D,Tcoff_)+weight*Tcnow
          {enddo\}
          
          ! do interpolation in neighbor grids that have the same level and pe
          {
            if (ixbmin^D<ixOmin^D) then
              idn^DD=0;
              idn^D=-1
              if (neighbor(2,idn^DD,igrid)==mype .and. neighbor_type(idn^DD,igrid)==neighbor_sibling) then
                igrid_nb=neighbor(1,idn^DD,igrid)
                ixcmin^DD=max(ixbmin^DD,ixOmin^DD);
                ixcmax^DD=min(ixbmax^DD,ixOmax^DD);
                ixcmin^D=ixOmax^D+(ixbmin^D-ixOmin^D)
                ixcmax^D=ixOmax^D
                xbmin^ND=rnode(rpxmin^ND_,igrid_nb)
                xbmax^ND=rnode(rpxmax^ND_,igrid_nb)
                ixmax^ND=floor((phys_trac_mask-xbmin^ND)/dxb^ND)+ixOmin^ND
                if (xbmax^ND>phys_trac_mask) ixcmax^ND=min(ixmax^ND,ixcmax^ND)

                {do ixc^DD=ixcmin^DD,ixcmax^DD;}
                  ds=0.d0
                  {ds=ds+(xFnow(^DD)-ps(igrid_nb)%x({ixc^DD},^DD))**2;}
                  ds=sqrt(ds)
                  if(ds<1.0d-2*dxb1) then
                    weight=(1/(1.0d-2*dxb1))**weightIndex
                  else
                    weight=(1/ds)**weightIndex
                  endif
                  ps(igrid_nb)%wextra(ixc^DD,Tweight_)=ps(igrid_nb)%wextra(ixc^DD,Tweight_)+weight
                  ps(igrid_nb)%wextra(ixc^DD,Tcoff_)=ps(igrid_nb)%wextra(ixc^DD,Tcoff_)+weight*Tcnow
                {enddo;}
              endif
            endif

            if (ixbmax^D>ixOmin^D) then
              idn^DD=0;
              idn^D=1
              if (neighbor(2,idn^DD,igrid)==mype .and. neighbor_type(idn^DD,igrid)==neighbor_sibling) then
                igrid_nb=neighbor(1,idn^DD,igrid)
                xbmin^ND=rnode(rpxmin^ND_,igrid_nb)
                if (xbmin^ND<phys_trac_mask) then
                  ixcmin^DD=max(ixbmin^DD,ixOmin^DD);
                  ixcmax^DD=min(ixbmax^DD,ixOmax^DD);
                  ixcmin^D=ixOmin^D
                  ixcmax^D=ixOmin^D+(ixbmax^D-ixOmax^D)
                  xbmax^ND=rnode(rpxmax^ND_,igrid_nb)
                  ixmax^ND=floor((phys_trac_mask-xbmin^ND)/dxb^ND)+ixOmin^ND
                  if (xbmax^ND>phys_trac_mask) ixcmax^ND=min(ixmax^ND,ixcmax^ND)

                  {do ixc^DD=ixcmin^DD,ixcmax^DD;}
                    ds=0.d0
                    {ds=ds+(xFnow(^DD)-ps(igrid_nb)%x({ixc^DD},^DD))**2;}
                    ds=sqrt(ds)
                    if(ds<1.0d-2*dxb1) then
                      weight=(1/(1.0d-2*dxb1))**weightIndex
                    else
                      weight=(1/ds)**weightIndex
                    endif
                    ps(igrid_nb)%wextra(ixc^DD,Tweight_)=ps(igrid_nb)%wextra(ixc^DD,Tweight_)+weight
                    ps(igrid_nb)%wextra(ixc^DD,Tcoff_)=ps(igrid_nb)%wextra(ixc^DD,Tcoff_)+weight*Tcnow
                  {enddo;}
                endif
              endif
            endif
          \}
        enddo

        iLP=ipmax+1
      enddo
    enddo


    do j=1,ngrid_trac
      igrid=trac_grid(j)
      where(ps(igrid)%wextra(ixO^S,Tweight_)>0.d0)
        ps(igrid)%wextra(ixO^S,Tcoff_)=ps(igrid)%wextra(ixO^S,Tcoff_)/ps(igrid)%wextra(ixO^S,Tweight_)
      elsewhere
        ps(igrid)%wextra(ixO^S,Tcoff_)=T_bott
      endwhere
    enddo

  end subroutine interp_Tcoff

end module mod_trac
