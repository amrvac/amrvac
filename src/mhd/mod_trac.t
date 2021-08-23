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
  integer :: group_world,group_trac,comm_trac,mype_trac,npe_trac
  logical, allocatable :: trac_pe(:)
contains

  subroutine init_trac_line(mask)
    logical, intent(in) :: mask 
    integer :: refine_factor,ix^D,ix(ndim),j,iFL,numL(ndim),finegrid
    double precision :: lengthFL
    double precision :: xprobmin(ndim),xprobmax(ndim),domain_nx(ndim)

    trac_delta=0.25d0
    !----------------- init magnetic field -------------------!
    refine_factor=2**(refine_max_level-1)
    ^D&xprobmin(^D)=xprobmin^D\
    ^D&xprobmax(^D)=xprobmax^D\
    ^D&domain_nx(^D)=domain_nx^D\
    dL=(xprobmax(ndim)-xprobmin(ndim))/(domain_nx(ndim)*refine_factor)
    ! max length of a field line
    if(mask) then
      {^IFTWOD
        lengthFL=(phys_trac_mask-xprobmin(2))*3.d0
      }
      {^IFTHREED
        lengthFL=(phys_trac_mask-xprobmin(3))*3.d0
      }
    else
      lengthFL=maxval(xprobmax-xprobmin)*3.d0
    end if
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

    if (mask) then
      call init_pegrid_mask()
      call create_trac_communicator()
      call init_trac_Tcoff()
    endif

  end subroutine init_trac_line

  subroutine update_mask_info()
    ! free current 
    call free_trac_communicator()
    ! record pe and grid for mask region
    call update_pegrid_mask()
    ! create a communicator for trac
    call update_trac_communicator()
    ! init vairable Tcoff_ and Tweight_
    call init_trac_Tcoff()

  end subroutine update_mask_info

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
    integer :: ixO^L,trac_tcoff

    ixO^L=ixM^LL;
    trac_tcoff=iw_tcoff
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
      ps(igrid)%w(ixO^S,trac_tcoff)=ps(igrid)%special_values(1)
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
      where(ps(igrid)%w(ixO^S,trac_tcoff) .lt. T_bott)
        ps(igrid)%w(ixO^S,trac_tcoff)=T_bott
      else where(ps(igrid)%w(ixO^S,trac_tcoff) .gt. 0.2d0*T_peak)
        ps(igrid)%w(ixO^S,trac_tcoff)=0.2d0*T_peak
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

    do j=1,ndim
      xF(:,1,j)=xFi(:,j)
    end do

    if (mask) then
      if (trac_pe(mype)) then
        ! convert conserved data to primitive for the calculation of temperature gradient
        call convert_data_trac()
        ! find out direction for tracing B field
        call get_Btracing_dir(ipel,igridl,forwardl)
        ! trace Bfield and get Tcoff for each field line
        call trace_Bfield_mask(xF,numR,Tlcoff,ipel,igridl,forwardl)
        ! get cell Tcoff via interpolation
        call interp_Tcoff_mask(xF,ipel,igridl,numR,Tlcoff)
        ! convert primitive data back to conserved data
        call convert_data_back_trac()
      endif
    else
      call init_trac_Tcoff()
      call get_Tlcoff(xF,numR,Tlcoff,ipel,igridl)
      call interp_Tcoff(xF,ipel,igridl,numR,Tlcoff)
    endif

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
      ps(igrid)%w(:^D&,Tweight_)=zero
      ps(igrid)%w(:^D&,Tcoff_)=zero
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
                ps(igrid)%w(ixc^D,Tweight_)=ps(igrid)%w(ixc^D,Tweight_)+weight
                ps(igrid)%w(ixc^D,Tcoff_)=ps(igrid)%w(ixc^D,Tcoff_)+weight*Tcoff_line(i)
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
      where (ps(igrid)%w(ixO^S,Tweight_)>0.d0)
        ps(igrid)%w(ixO^S,Tcoff_)=ps(igrid)%w(ixO^S,Tcoff_)/ps(igrid)%w(ixO^S,Tweight_)
      elsewhere
        ps(igrid)%w(ixO^S,Tcoff_)=0.2d0*Tmax
      endwhere
    enddo
  end subroutine block_interp_grid

  subroutine get_Tmax_grid(x,w,ixI^L,ixO^L,Tmax_grid)
    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(out) :: w(ixI^S,1:nw)
    double precision :: Tmax_grid
    double precision :: pth(ixI^S),Te(ixI^S)

    call mhd_get_pthermal(w,x,ixI^L,ixO^L,pth)
    Te(ixO^S)=pth(ixO^S)/w(ixO^S,rho_)
    Tmax_grid=maxval(Te(ixO^S))
  end subroutine get_Tmax_grid  

  subroutine get_Tlcoff(xF,numR,Tlcoff,ipel,igridl)
    use mod_point_searching
    double precision :: xF(numFL,numLP,ndim)
    integer :: numR(numFL),ix1
    double precision :: Tlcoff(numFL)
    integer :: ipel(numFL,numLP),igridl(numFL,numLP)
    integer :: numRT
    double precision :: xFt(numLP,ndim)
    integer :: ipef(numLP),igridf(numLP)
    double precision :: Tcofl
    logical :: forward
    double precision :: xp(1:ndim),wp(1:nw)

    do ix1=1,numFL
      xFt(1,:)=xF(ix1,1,:)
      xp=xFt(1,:)
      call get_point_w(xp,wp,'conserved')
      if (wp(mag(ndim))<=0) then
        forward=.false.
      else
        forward=.true.
      endif
      call trace_Bfield_trac(xFt,ipef,igridf,Tcofl,numLP,numRT,forward)
      xF(ix1,1:numRT,:)=xFt(1:numRT,:)
      numR(ix1)=numRT
      ipel(ix1,1:numRT)=ipef(1:numRT)
      igridl(ix1,1:numRT)=igridf(1:numRT)
      Tlcoff(ix1)=Tcofl
    enddo
  end subroutine get_Tlcoff

  subroutine interp_Tcoff(xF,ipel,igridl,numR,Tlcoff)
    double precision :: xF(numFL,numLP,ndim)
    integer :: numR(numFL),ipel(numFL,numLP),igridl(numFL,numLP)
    double precision :: Tlcoff(numFL)
    double precision :: dxMax^D,weight,ds
    integer :: ix1,ix2,ixT^L,ixT^D,weightIndex
    integer :: ixO^L,ixc^L,igrid,iigrid,ixc^D
    double precision :: dxb^D

    weightIndex=2
    ^D&ixOmin^D=ixmlo^D\
    ^D&ixOmax^D=ixmhi^D\
    do ix1=1,numFL
      do ix2=1,numR(ix1)
        if (ipel(ix1,ix2)==mype) then
          igrid=igridl(ix1,ix2)
          ^D&dxb^D=rnode(rpdx^D_,igrid)\
          ^D&dxMax^D=4*dxb^D\ 
          ^D&ixcmin^D=floor((xF(ix1,ix2,^D)-dxMax^D-ps(igrid)%x(ixOmin^DD,^D))/dxb^D)+ixOmin^D\
          ^D&ixcmax^D=floor((xF(ix1,ix2,^D)+dxMax^D-ps(igrid)%x(ixOmin^DD,^D))/dxb^D)+ixOmin^D\
          {if(ixcmin^D<ixOmin^D) ixcmin^D=ixOmin^D\}
          {if(ixcmax^D>ixOmax^D) ixcmax^D=ixOmax^D\}
          if (ps(igrid)%x(ixcmax^D,ndim)>phys_trac_mask) then
            ixcmax^ND=floor((phys_trac_mask-ps(igrid)%x(ixOmin^D,ndim))/dxb^ND)+ixOmin^ND
          endif
          {do ixc^D=ixcmin^D,ixcmax^D\}
            ds=0.d0
            {ds=ds+(xF(ix1,ix2,^D)-ps(igrid)%x(ixc^DD,^D))**2\}
            ds=sqrt(ds)
            if(ds<1.0d-2*dxb1) then
              weight=(1/(1.0d-2*dxb1))**weightIndex
            else
              weight=(1/ds)**weightIndex
            endif
            ps(igrid)%w(ixc^D,Tweight_)=ps(igrid)%w(ixc^D,Tweight_)+weight
            ps(igrid)%w(ixc^D,Tcoff_)=ps(igrid)%w(ixc^D,Tcoff_)+weight*Tlcoff(ix1)
          {enddo\}
        endif
      enddo
    enddo
    ! finish interpolation
    do iigrid=1,igridstail; igrid=igrids(iigrid);
      where(ps(igrid)%w(ixO^S,Tweight_)>0.d0)
        ps(igrid)%w(ixO^S,Tcoff_)=ps(igrid)%w(ixO^S,Tcoff_)/ps(igrid)%w(ixO^S,Tweight_)
      elsewhere
        ps(igrid)%w(ixO^S,Tcoff_)=0.d0
      endwhere
    enddo
  end subroutine interp_Tcoff

  subroutine trace_Bfield_trac(xf,ipef,igridf,Tcofl,numP,numRT,forward)
    ! trace a field line
    use mod_usr_methods
    use mod_global_parameters
    use mod_particle_base
    integer :: numP,numRT
    double precision :: xf(numP,ndim)
    integer :: ipef(numP),igridf(numP)
    double precision :: Tcofl,Tlmax
    logical :: forward
    double precision :: dxb^D,xb^L
    integer :: indomain
    integer :: igrid,igrid_now,igrid_next,j
    integer :: ipe_now,ipe_next
    double precision :: xp_in(ndim),xp_out(ndim),x3d(3)
    integer :: ipoint_in,ipoint_out
    double precision :: statusB(7+2*ndim)
    logical :: stopB
    double precision :: Te_info(3)

    xf(2:numP,:)=0
    ipef=-1
    igridf=0
    Te_info=0.d0  ! Te_pre,Tcofl,Tlmax
    ! check whether or the first point is inside simulation box. if yes, find
    ! the pe and igrid for the point
    indomain=0
    numRT=0
    {if (xf(1,^DB)>=xprobmin^DB .and. xf(1,^DB)<xprobmax^DB) indomain=indomain+1\}
    if(indomain==ndim) then
      numRT=1
       ! find pe and igrid
       x3d=0.d0
       do j=1,ndim
         x3d(j)=xf(1,j)
       enddo
      call find_particle_ipe(x3d,igrid_now,ipe_now)
      stopB=.FALSE.
      ipoint_in=1
    else
      if (mype==0) then
        call MPISTOP('magnetic field tracing error: given point is not in simulation box!')
      endif
    endif
    ! other points in field line
    do while(stopB .eqv. .FALSE.)
      if (mype==ipe_now) then
        igrid=igrid_now
        ! looking for points in one pe
        call trace_in_pe(igrid,ipoint_in,xf,igridf,numP,forward,statusB,Te_info)
      endif
      ! comunication
      call MPI_BCAST(statusB,7+2*ndim,MPI_DOUBLE_PRECISION,ipe_now,icomm,ierrmpi)
      ! prepare for next step
      ipoint_out=int(statusB(1))
      ipe_next=int(statusB(2))
      igrid_next=int(statusB(3))
      if (int(statusB(4))==1) then
        stopB=.TRUE.
        numRT=ipoint_out-1
      endif
      xf(ipoint_out-1,1:ndim)=statusB(4+1:4+ndim)
      xf(ipoint_out,1:ndim)=statusB(4+ndim+1:4+ndim+ndim)
      Te_info(1)=statusB(5+2*ndim)   ! temperature at ipoint_out-1
      Te_info(2)=statusB(6+2*ndim)   ! current Tcofl
      Te_info(3)=statusB(7+2*ndim)   ! current Tlmax
      ipef(ipoint_in:ipoint_out-1)=ipe_now
      ! pe and grid of next point
      ipe_now=ipe_next
      igrid_now=igrid_next
      ipoint_in=ipoint_out
    enddo

    ! get cutoff temperature
    Tcofl=Te_info(2)
    Tlmax=Te_info(3)
    if(Tcofl>0.2d0*Tlmax) Tcofl=0.2d0*Tlmax
    if(Tcofl<T_bott) Tcofl=T_bott

  end subroutine trace_Bfield_trac

  subroutine trace_in_pe(igrid,ipoint_in,xf,igridf,numP,forward,statusB,Te_info)
    integer :: igrid,ipoint_in,numP
    double precision :: xf(numP,ndim)
    integer :: igridf(numP)
    logical :: forward
    double precision :: statusB(7+2*ndim),Te_info(3)
    integer :: ipe_next,igrid_next,ip_in,ip_out,j
    logical :: newpe,stopB
    double precision :: xfout(ndim)
    integer :: indomain

    ip_in=ipoint_in
    newpe=.FALSE.
    do while(newpe .eqv. .FALSE.)
      ! looking for points in given grid    
      call find_points_trac(igrid,ip_in,ip_out,xf,Te_info,numP,forward)
      igridf(ip_in:ip_out-1)=igrid
      ip_in=ip_out

      indomain=0
      {if (xf(ip_out,^DB)>=xprobmin^DB .and. xf(ip_out,^DB)<xprobmax^DB) indomain=indomain+1\}

      ! when next point is out of domain, stop
      ! when next point is out of given grid, find next grid  
      if (indomain/=ndim) then
        newpe=.TRUE.
        stopB=.TRUE.
      else if(ip_out<numP) then
        stopB=.FALSE.
        xfout=xf(ip_out,:)
        call find_next_grid_trac(igrid,igrid_next,ipe_next,xfout,newpe,stopB)
      else
        newpe=.TRUE.
        stopB=.TRUE.
      endif
      if(newpe) then
        statusB(1)=ip_out
        statusB(2)=ipe_next
        statusB(3)=igrid_next
        statusB(4)=0
        if (stopB) statusB(4)=1
        statusB(4+1:4+ndim)=xf(ip_out-1,1:ndim)
        statusB(4+ndim+1:4+ndim+ndim)=xf(ip_out,1:ndim)
        statusB(5+2*ndim)=Te_info(1)
        statusB(6+2*ndim)=Te_info(2)
        statusB(7+2*ndim)=Te_info(3)
      endif
      if(newpe .eqv. .FALSE.) igrid=igrid_next
    enddo
  end subroutine trace_in_pe

  subroutine find_points_trac(igrid,ip_in,ip_out,xf,Te_info,numP,forward)
    integer :: igrid,ip_in,ip_out,numP
    double precision :: xf(numP,ndim),Te_info(3)
    logical :: forward

    double precision :: dxb^D,xb^L
    integer          :: ip,inblock,ixI^L,j
    double precision :: Bg(ixg^T,ndir),pth(ixg^T),Te(ixg^T)
    double precision :: Bnow(ndir),Bnext(ndir)
    double precision :: xfpre(ndim),xfnow(ndim),xfnext(ndim)
    double precision :: Tpre,Tnow,Tnext,dTds,Lt,Lr,ds

    ^D&ixImin^D=ixglo^D;
    ^D&ixImax^D=ixghi^D;
    ^D&dxb^D=rnode(rpdx^D_,igrid);
    ^D&xbmin^D=rnode(rpxmin^D_,igrid);
    ^D&xbmax^D=rnode(rpxmax^D_,igrid);

    ! get grid magnetic field and temperature
    if (B0field) then
      do j=1,ndir
        Bg(ixI^S,j)=ps(igrid)%w(ixI^S,mag(j))+ps(igrid)%B0(ixI^S,j,0)
      enddo
    else
      do j=1,ndir
        Bg(ixI^S,j)=ps(igrid)%w(ixI^S,mag(j))
      enddo
    endif
    call mhd_get_pthermal(ps(igrid)%w,ps(igrid)%x,ixI^L,ixI^L,pth)
    Te(ixI^S)=pth(ixI^S)/ps(igrid)%w(ixI^S,rho_)

    ! main loop
    MAINLOOP: do ip=ip_in,numP-1
      ! get temperature to calculate dT/ds
      xfnow(:)=xf(ip,:)
      call get_BTe_local(xfnow,Bnow,Tnow,ps(igrid)%x,Bg,Te,ixI^L,dxb^D)
      if (ip>1) then 
        xfpre(:)=xf(ip-1,:)
        Tpre=Te_info(1)   ! temperature of previous point
      else
        xfpre(:)=0.d0
        Tpre=0.d0
      endif
      ds=dxb1
      call get_x_next(ip,xfpre,xfnow,Bnow,xfnext,ds,forward)
      xf(ip+1,:)=xfnext(:)
      call get_BTe_local(xfnext,Bnext,Tnext,ps(igrid)%x,Bg,Te,ixI^L,dxb^D)

      ! calculate dT/ds and update Te_info
      call get_dTds(ip,xfpre,xfnow,xfnext,Tpre,Tnow,Tnext,dTds)
      Te_info(1)=Tnow   ! temperature of previous point
      if (ip==1) then
        Lt=0.d0
        Te_info(2)=T_bott   ! current Tcofl
        Te_info(3)=Tnow     ! current Tlmax
      else
        Lt=0.d0
        if (dTds>0.d0) then
          Lt=Tnow/dTds
          Lr=ds
          ! renew cutoff temperature
          if(Lr>trac_delta*Lt) then
            if (Tnow>Te_info(2)) Te_info(2)=Tnow
          endif
        endif
        if (Tnow>Te_info(3)) Te_info(3)=Tnow
      endif

      ip_out=ip+1
      ! whether or not next point is in this block/grid
      inblock=0
      {if(xf(ip+1,^DB)>=xbmin^DB .and. xf(ip+1,^DB)<xbmax^DB) inblock=inblock+1\}
      if(inblock/=ndim) then
        ! exit loop if next point is not in this block
        exit MAINLOOP
      endif

    enddo MAINLOOP

  end subroutine find_points_trac

  subroutine get_dTds(ip,xfpre,xfnow,xfnext,Tpre,Tnow,Tnext,dTds)

    integer :: ip,j
    double precision :: xfpre(ndim),xfnow(ndim),xfnext(ndim)
    double precision :: Tpre,Tnow,Tnext
    double precision :: ds1,ds2,dT1,dT2,dTds

    ds1=0.d0
    ds2=0.d0
    do j=1,ndim
      ds1=ds1+(xfnext(j)-xfnow(j))**2
      if (ip>1) ds2=ds2+(xfpre(j)-xfnow(j))**2
    enddo
    ds1=sqrt(ds1)
    ds2=sqrt(ds2)
    dT1=Tnext-Tnow
    dT2=Tpre-Tnow

    if (ip>1) then
      dTds=abs(dT1*ds2**2-dT2*ds1**2)/(ds1*ds2*(ds1+ds2))
    else
      dTds=abs(dT1/ds1)
    endif

  end subroutine get_dTds

  subroutine get_x_next(ip,xfpre,xfnow,Bnow,xfnext,ds,forward)

    integer :: ip,j
    double precision :: xfpre(1:ndim),xfnow(1:ndim),xfnext(1:ndim),Bnow(1:ndir)
    double precision :: ds,Btotal,dxf(1:ndim),dxftot
    logical :: forward

    Btotal=0.d0
    do j=1,ndim
      Btotal=Btotal+(Bnow(j))**2
    enddo
    Btotal=dsqrt(Btotal)

    ! if local magnetic field is 0
    if(Btotal==0) then
      if (ip>1) then
        dxftot=0.d0
        do j=1,ndim
          dxf(j)=xfnow(j)-xfpre(j)
          dxftot=dxftot+dxf(j)**2
        enddo
        dxftot=sqrt(dxftot)
        dxf(:)=dxf(:)*ds/dxftot
      else
        dxf(:)=0.d0
        dxf(ndim)=ds
      endif
    else
    ! find next point based on magnetic field direction
      if(forward .eqv. .TRUE.) then
        do j=1,ndim
          dxf(j)=ds*Bnow(j)/Btotal
        enddo
      else
        do j=1,ndim
          dxf(j)=-ds*Bnow(j)/Btotal
        enddo
      endif
    endif

    ! next point
    do j=1,ndim
      xfnext(j)=xfnow(j)+dxf(j)
    enddo

  end subroutine get_x_next

  subroutine get_BTe_local(xfloc,Bloc,Tloc,x,Bg,Te,ixI^L,dxb^D)

    integer :: ixI^L
    double precision :: dxb^D,Tloc
    double precision :: xfloc(ndim),Bloc(ndir)
    double precision :: x(ixI^S,ndim),Bg(ixI^S,ndir),Te(ixI^S)

    integer          :: ixb^D,ix^D,ixbl^D,j
    double precision :: xd^D
    double precision :: factor(0:1^D&)
    double precision :: Bnear(0:1^D&,ndir),Tnear(0:1^D&)

    ^D&ixbl^D=floor((xfloc(^D)-x(ixImin^DD,^D))/dxb^D)+ixImin^D\
    ^D&xd^D=(xfloc(^D)-x(ixbl^DD,^D))/dxb^D\

    {do ix^D=0,1\}
      factor(ix^D)={abs(1-ix^D-xd^D)*}
    {enddo\}

    ! Bfield and Te for interpolation
    {do ix^DB=0,1\}
      do j=1,ndir
        Bnear(ix^D,j)=Bg(ixbl^D+ix^D,j)
      enddo
      Tnear(ix^D)=Te(ixbl^D+ix^D)
    {enddo\}

    Bloc(:)=0.d0
    Tloc=0.d0
    ! interpolation
    {do ix^DB=0,1\}
      do j=1,ndir
        Bloc(j)=Bloc(j)+Bnear(ix^D,j)*factor(ix^D)
      enddo
      Tloc=Tloc+Tnear(ix^D)*factor(ix^D)
    {enddo\}

  end subroutine get_BTe_local

  subroutine find_next_grid_trac(igrid,igrid_next,ipe_next,xf1,newpe,stopB)
    ! check the grid and pe of next point
    use mod_usr_methods
    use mod_forest
    integer :: igrid,igrid_next,ipe_next
    double precision :: xf1(ndim)
    double precision :: dxb^D,xb^L,xbmid^D
    logical :: newpe,stopB
    integer :: idn^D,my_neighbor_type,inblock
    integer :: ic^D,inc^D,ipe_neighbor,igrid_neighbor
    double precision :: xbn^L

    ^D&xbmin^D=rnode(rpxmin^D_,igrid)\
    ^D&xbmax^D=rnode(rpxmax^D_,igrid)\
    inblock=0
    ! direction of next grid
    idn^D=0\
    {if (xf1(^D)<=xbmin^D) idn^D=-1\}
    {if (xf1(^D)>=xbmax^D) idn^D=1\}
    my_neighbor_type=neighbor_type(idn^D,igrid)
    igrid_neighbor=neighbor(1,idn^D,igrid)
    ipe_neighbor=neighbor(2,idn^D,igrid)
    ! ipe and igrid of next grid
    select case(my_neighbor_type)
    case(neighbor_boundary)
      ! next point is not in simulation box
      newpe=.TRUE.
      stopB=.TRUE.
    case(neighbor_coarse)
      ! neighbor grid has lower refinement level      
      igrid_next=igrid_neighbor
      ipe_next=ipe_neighbor
      if(mype==ipe_neighbor) then
        newpe=.FALSE.
      else
        newpe=.TRUE.
      endif
    case(neighbor_sibling)
      ! neighbor grid has lower refinement level 
      igrid_next=igrid_neighbor
      ipe_next=ipe_neighbor
      if(mype==ipe_neighbor) then
        newpe=.FALSE.
      else
        newpe=.TRUE.
      endif
    case(neighbor_fine)
      ! neighbor grid has higher refinement level 
      {xbmid^D=(xbmin^D+xbmax^D)/2.d0\}
      ^D&inc^D=1\
      {if(xf1(^D)<=xbmin^D) inc^D=0\}
      {if(xf1(^D)>xbmin^D .and. xf1(^D)<=xbmid^D) inc^D=1\}
      {if(xf1(^D)>xbmid^D .and. xf1(^D)<xbmax^D) inc^D=2\}
      {if(xf1(^D)>=xbmax^D) inc^D=3\}
      ipe_next=neighbor_child(2,inc^D,igrid)
      igrid_next=neighbor_child(1,inc^D,igrid)
      if(mype==ipe_next) then
        newpe=.FALSE.
      else
        newpe=.TRUE.
      endif
    end select
  end subroutine find_next_grid_trac

  subroutine init_trac_Tcoff()
    integer :: ixI^L,ixO^L,igrid,iigrid

    ^D&ixOmin^D=ixmlo^D;
    ^D&ixOmax^D=ixmhi^D;
    do iigrid=1,igridstail; igrid=igrids(iigrid);
      ps(igrid)%w(ixO^S,Tcoff_)=0.d0
      ps(igrid)%w(ixO^S,Tweight_)=0.d0
    end do

  end subroutine init_trac_Tcoff

  subroutine init_pegrid_mask()
    use mod_forest, only: tree_node_ptr, tree_root
    integer :: iglv1^L,igm^D
    type(tree_node_ptr) :: branch

    allocate(trac_pe(0:npe-1),trac_grid(max_blocks),ground_grid(max_blocks))

    npe_trac=0
    ngrid_trac=0
    ngrid_ground=0
    trac_pe(:)=.false.
    trac_grid(:)=0
    ground_grid(:)=0

    ! traverse tree to get information of grid inside mask region
    ^D&iglv1min^D=1;
    ^D&iglv1max^D=ceiling((xprobmax^D-xprobmin^D)/dg^D(1));
    iglv1max^ND=ceiling((phys_trac_mask-xprobmin^ND)/dg^ND(1))
    {do igm^D=iglv1min^D,iglv1max^D\}
      branch=tree_root(igm^D)
      call traverse_tree_mask(branch)
    {enddo\}

  end subroutine init_pegrid_mask

  subroutine update_pegrid_mask()
    use mod_forest, only: tree_node_ptr, tree_root
    integer :: iglv1^L,igm^D
    type(tree_node_ptr) :: branch

    npe_trac=0
    ngrid_trac=0
    ngrid_ground=0
    trac_pe(:)=.false.
    trac_grid(:)=0
    ground_grid(:)=0

    ! traverse tree to get information of grid inside mask region
    ^D&iglv1min^D=1;
    ^D&iglv1max^D=ceiling((xprobmax^D-xprobmin^D)/dg^D(1));
    iglv1max^ND=ceiling((phys_trac_mask-xprobmin^ND)/dg^ND(1))
    {do igm^D=iglv1min^D,iglv1max^D\}
      branch=tree_root(igm^D)
      call traverse_tree_mask(branch)
    {enddo\}

  end subroutine update_pegrid_mask

  subroutine create_trac_communicator()
    integer, allocatable :: ranks_trac(:)
    integer :: ipe,j

    do j=0,npe-1
      if (trac_pe(j)) then
        npe_trac=npe_trac+1
      endif
    enddo

    allocate(ranks_trac(0:npe_trac-1))
    ipe=0
    do j=0,npe-1
      if (trac_pe(j)) then
        ranks_trac(ipe)=j
        ipe=ipe+1
      endif
    enddo

    call MPI_COMM_GROUP(icomm,group_world,ierrmpi)
    call MPI_GROUP_INCL(group_world,npe_trac,ranks_trac,group_trac,ierrmpi)
    call MPI_COMM_CREATE(icomm,group_trac,comm_trac,ierrmpi)

    deallocate(ranks_trac)

  end subroutine create_trac_communicator

  subroutine free_trac_communicator()

    if (trac_pe(mype)) call MPI_COMM_FREE(comm_trac,ierrmpi)
    call MPI_GROUP_FREE(group_trac,ierrmpi)

  end subroutine free_trac_communicator

  subroutine update_trac_communicator()
    integer, allocatable :: ranks_trac(:)
    integer :: ipe,j,numProcs

    do j=0,npe-1
      if (trac_pe(j)) then
        npe_trac=npe_trac+1
      endif
    enddo

    allocate(ranks_trac(0:npe_trac-1))
    ipe=0
    do j=0,npe-1
      if (trac_pe(j)) then
        ranks_trac(ipe)=j
        ipe=ipe+1
      endif
    enddo

    call MPI_GROUP_INCL(group_world,npe_trac,ranks_trac,group_trac,ierrmpi)
    call MPI_COMM_CREATE(icomm,group_trac,comm_trac,ierrmpi)

  end subroutine update_trac_communicator

  recursive subroutine traverse_tree_mask(branch)
    use mod_forest, only: tree_node_ptr, tree_root

    type(tree_node_ptr) :: branch,sub_branch
    integer :: ic^D,ic^L,ilevel,igmax^ND

    if (.not. branch%node%leaf) then
      ^D&icmin^D=1;
      ^D&icmax^D=2;
      ilevel=branch%node%level+1
      igmax^ND=ceiling((phys_trac_mask-xprobmin^ND)/dg^ND(ilevel))
      icmax^ND=min(2,igmax^ND-(2*branch%node%ig^ND)+2)
      {do ic^D=icmin^D,icmax^D\}
        sub_branch=branch%node%child(ic^D)
        call traverse_tree_mask(sub_branch)
      {enddo\}
    else
      trac_pe(branch%node%ipe)=.true.
      if (mype==branch%node%ipe) then
        ngrid_trac=ngrid_trac+1
        trac_grid(ngrid_trac)=branch%node%igrid
        if (branch%node%ig^ND==1) then
          ngrid_ground=ngrid_ground+1
          ground_grid(ngrid_ground)=branch%node%igrid
        endif
      endif
    endif

  end subroutine traverse_tree_mask

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
          Bh=Bh+factor*ps(igrid)%w(ixb^D+ixbb^D,mag(^ND))
        {enddo\}
        if (Bh>0) then
          forwardl(iFL)=.TRUE.
        else
          forwardl(iFL)=.FALSE.
        endif
      {enddo\}
    enddo

    call MPI_ALLREDUCE(forwardl,forwardRC,numFL,MPI_LOGICAL,&
                         MPI_LAND,comm_trac,ierrmpi)
    forwardl=forwardRC

  end subroutine get_Btracing_dir

  subroutine convert_data_trac()
    integer :: ixI^L,ixO^L
    integer :: igrid,j,idir

    ^D&ixImin^D=ixglo^D;
    ^D&ixImax^D=ixghi^D;
    ^D&ixOmin^D=ixmlo^D;
    ^D&ixOmax^D=ixmhi^D;

    do j=1,ngrid_trac
      igrid=trac_grid(j)
      call mhd_to_primitive(ixI^L,ixI^L,ps(igrid)%w,ps(igrid)%x)
      if (B0field) then
        do idir=1,ndir
          ps(igrid)%w(ixI^S,mag(idir))=ps(igrid)%w(ixI^S,mag(idir))+ps(igrid)%B0(ixI^S,idir,0)
        enddo
      endif
      ! init Tcoff_ and Tweight_ for masked grids
      ps(igrid)%w(ixO^S,Tcoff_)=0.d0
      ps(igrid)%w(ixO^S,Tweight_)=0.d0
    enddo

  endsubroutine convert_data_trac

  subroutine output_trac_grid(ngrid_ground,ground_grid)
    integer :: ngrid_ground
    integer :: ground_grid(max_blocks)
    character (30) :: fname,tempc
    integer :: igrid,iigrid,j
    double precision :: xb^L

    write(fname, '(a,i2.2,a)') 'trac_grid_',mype,'.txt'
    open(1,file=fname)
    write(1,*) 'ngrid_trac'
    !write(1,*) ngrid_trac
    write(1,*) ngrid_ground
    write(1,*) 'igrid  xbmin1 xbmin2 xbmax1 xbmax2'
    !do j=1,ngrid_trac
    do j=1,ngrid_ground
      !igrid=trac_grid(j)
      igrid=ground_grid(j)
      ^D&xbmin^D=rnode(rpxmin^D_,igrid);
      ^D&xbmax^D=rnode(rpxmax^D_,igrid);
      write(1,'(i8, e15.7, e15.7, e15.7, e15.7)') igrid,xb^L
    enddo
    close(1)

  end subroutine output_trac_grid

  subroutine trace_Bfield_mask(xFL,numR,TcoffFL,ipeFL,igridFL,forwardFL)
    use mod_particle_base

    double precision :: xFL(numFL,numLP,ndim)
    integer :: numR(numFL)
    double precision :: TcoffFL(numFL)
    integer :: ipeFL(numFL,numLP),igridFL(numFL,numLP)
    logical :: forwardFL(numFL)

    integer :: iFL,iLP,j,numSend,ipe_next,igrid_next,step
    integer :: ipoint_in(numFL)
    logical :: stopB
    logical :: myFL(numFL),continueFL(numFL)
    double precision :: statusFL(numFL,5+ndim),statusRC(numFL,5+ndim)

    step=0
    numR=0
    stopB=.FALSE.
    myFL=.FALSE.
    continueFL=.TRUE.
    ipoint_in=1
    TcoffFL=0

    do iFL=1,numFL
      if (ipeFL(iFL,1)==mype) myFL(iFL)=.TRUE.
    enddo

    do while(stopB .eqv. .FALSE.)
      ! tracing multiple field lines inside pe
      statusFL=0
      do iFL=1,numFL
        if (myFL(iFL) .and. continueFL(iFL)) then
          call traceFL_mask(xFL,ipeFL,igridFL,ipoint_in,iFL,forwardFL,statusFL,TcoffFL)
          numR(iFL)=int(statusFL(iFL,1))-1
          !print *, iFL,mype,statusFL(iFL,:)
        endif
      enddo

      step=step+1

      ! comunicating tracing results
      numSend=numFL*(5+ndim)
      call MPI_ALLREDUCE(statusFL,statusRC,numSend,MPI_DOUBLE_PRECISION,&
                       MPI_SUM,comm_trac,ierrmpi)
      statusFL=statusRC

      ! preparing for next tracing
      ! statusFL
      ! 1: ipoint_in(iFL)
      ! 2: continueFL(iFL)
      ! 3: ipe_next
      ! 4: igrid_next
      ! 5: current Tcofl
      ! 5+1:5+ndim: xFL(iFL,ipoint_in(iFL)
      stopB=.TRUE.
      myFL=.FALSE.
      continueFL=.TRUE.
      do iFL=1,numFL
        ipoint_in(iFL)=int(statusFL(iFL,1))
        if (int(statusFL(iFL,2))==0) continueFL(iFL)=.FALSE.
        if (continueFL(iFL)) then
          stopB=.FALSE.
          ipe_next=int(statusFL(iFL,3))
          if (mype==ipe_next) then
            myFL(iFL)=.TRUE.
            igrid_next=int(statusFL(iFL,4))
            igridFL(iFL,ipoint_in(iFL))=igrid_next
            xFL(iFL,ipoint_in(iFL),1:ndim)=statusFL(iFL,5+1:5+ndim)
          endif
        endif
        TcoffFL(iFL)=max(statusFL(iFL,5),TcoffFL(iFL))
      enddo
    enddo

    ! update Tcoff base on tracing results
    do iFL=1,numFL
      if(TcoffFL(iFL)>0.2d0*Tmax) TcoffFL(iFL)=0.2d0*Tmax
      if(TcoffFL(iFL)<T_bott) TcoffFL(iFL)=T_bott
    enddo

  end subroutine trace_Bfield_mask

  subroutine traceFL_mask(xFL,ipeFL,igridFL,ipoint_in,iFL,forwardFL,statusFL,TcoffFL)
    double precision :: xFL(numFL,numLP,ndim)
    integer :: ipeFL(numFL,numLP),igridFL(numFL,numLP),ipoint_in(numFL)
    integer :: iFL
    logical :: forwardFL(numFL)
    double precision :: statusFL(numFL,5+ndim),TcoffFL(numFL)

    integer :: ip_in,ip_out,igrid,indomain,ipe_next,igrid_next
    logical :: forward,newpe,stopB
    double precision :: xfout(ndim),Tcoff

    ip_in=ipoint_in(iFL)
    forward=forwardFL(iFL)
    igrid=igridFL(iFL,ip_in)
    newpe=.FALSE.
    Tcoff=TcoffFL(iFL)

    do while(newpe .eqv. .FALSE.)
      ! looking for points in given grid    
      call find_points_mask(xFL,igrid,iFL,ip_in,ip_out,Tcoff,forward)
      igridFL(iFL,ip_in:ip_out-1)=igrid
      ipeFL(iFL,ip_in:ip_out-1)=mype
      ip_in=ip_out

      indomain=0
      {if (xFL(iFL,ip_out,^DB)>=xprobmin^DB .and. xFL(iFL,ip_out,^DB)<xprobmax^DB) indomain=indomain+1\}

      ! when next point is out of domain, stop
      ! when next point is out of given grid, find next grid  
      if (indomain/=ndim) then
        newpe=.TRUE.
        stopB=.TRUE.
      else if(ip_out<numLP .and. xFL(iFL,ip_out,ndim)<phys_trac_mask) then
        stopB=.FALSE.
        xfout=xFL(iFL,ip_out,:)
        call find_next_grid_trac(igrid,igrid_next,ipe_next,xfout,newpe,stopB)
      else
        newpe=.TRUE.
        stopB=.TRUE.
      endif
      if(newpe) then
        statusFL(iFL,1)=ip_out
        statusFL(iFL,2)=1
        if (stopB) statusFL(iFL,2)=0
        statusFL(iFL,3)=ipe_next
        statusFL(iFL,4)=igrid_next
        statusFL(iFL,5)=Tcoff
        statusFL(iFL,5+1:5+ndim)=xFL(iFL,ip_out,1:ndim)
      endif
      if(newpe .eqv. .FALSE.) igrid=igrid_next
    enddo
  end subroutine traceFL_mask

  subroutine find_points_mask(xFL,igrid,iFL,ip_in,ip_out,Tcoff,forward)
    double precision :: xFL(numFL,numLP,ndim)
    integer :: igrid,iFL,ip_in,ip_out
    double precision :: Tcoff
    logical :: forward

    double precision :: dxb^D,xb^L
    integer          :: ip,inblock,ixI^L,j
    double precision :: Bnow(ndir),Bnext(ndir),Bpre(ndir)
    double precision :: xfpre(ndim),xfnow(ndim),xfnext(ndim)
    double precision :: Tpre,Tnow,Tnext,dTds,Lt,Lr,ds

    ^D&ixImin^D=ixglo^D;
    ^D&ixImax^D=ixghi^D;
    ^D&dxb^D=rnode(rpdx^D_,igrid);
    ^D&xbmin^D=rnode(rpxmin^D_,igrid);
    ^D&xbmax^D=rnode(rpxmax^D_,igrid);

    ! main loop
    MAINLOOP: do ip=ip_in,numLP-1
      xfnow(:)=xFL(iFL,ip,:)
      ds=dxb^ND

      ! get temperature to calculate dT/ds
      if (ip==ip_in) then
        call get_localBT_mask(xfnow,Bnow,Tnow,ps(igrid)%x,ps(igrid)%w,ixI^L,dxb^D)
        ! location and temperature of next point
        call get_x_neighbor(ip,xfnow,Bnow,xfnext,ds,forward)
        xFL(iFL,ip+1,:)=xfnext(:)
        call get_localBT_mask(xfnext,Bnext,Tnext,ps(igrid)%x,ps(igrid)%w,ixI^L,dxb^D)
        ! location and temperature of previous point
        call get_x_neighbor(ip,xfnow,Bnow,xfpre,-ds,forward)
        call get_localBT_mask(xfpre,Bpre,Tpre,ps(igrid)%x,ps(igrid)%w,ixI^L,dxb^D)
      else
        xfpre(:)=xFL(iFL,ip-1,:)
        Tpre=Tnow
        Tnow=Tnext
        Bnow=Bnext
        call get_x_neighbor(ip,xfnow,Bnow,xfnext,ds,forward)
        xFL(iFL,ip+1,:)=xfnext(:)
        call get_localBT_mask(xfnext,Bnext,Tnext,ps(igrid)%x,ps(igrid)%w,ixI^L,dxb^D)
      endif

      ! calculate dT/ds and update Te_info
      call get_dTds(ip,xfpre,xfnow,xfnext,Tpre,Tnow,Tnext,dTds)
      if (ip==1) then
        Lt=0.d0
        Tcoff=T_bott 
      else
        Lt=0.d0
        if (dTds>0.d0) then
          Lt=Tnow/dTds
          Lr=ds
          ! renew cutoff temperature
          if(Lr>trac_delta*Lt) then
            if (Tnow>Tcoff) Tcoff=Tnow
          endif
        endif
      endif

      ip_out=ip+1
      ! whether or not next point is in this block/grid
      inblock=0
      {if(xFL(iFL,ip+1,^DB)>=xbmin^DB .and. xFL(iFL,ip+1,^DB)<xbmax^DB) inblock=inblock+1\}
      if(inblock/=ndim .or. xFL(iFL,ip+1,ndim)>phys_trac_mask) then
        ! exit loop if next point is not in this block
        exit MAINLOOP
      endif

    enddo MAINLOOP

  end subroutine find_points_mask

  subroutine get_localBT_mask(xfloc,Bloc,Tloc,x,w,ixI^L,dxb^D)

    integer :: ixI^L
    double precision :: dxb^D,Tloc
    double precision :: xfloc(ndim),Bloc(ndir)
    double precision :: x(ixI^S,ndim),w(ixI^S,nw)

    integer          :: ixb^D,ix^D,ixbl^D,j
    double precision :: xd^D
    double precision :: factor(0:1^D&)
    double precision :: Bnear(0:1^D&,ndir),Tnear(0:1^D&)

    ^D&ixbl^D=floor((xfloc(^D)-x(ixImin^DD,^D))/dxb^D)+ixImin^D\
    ^D&xd^D=(xfloc(^D)-x(ixbl^DD,^D))/dxb^D\

    {do ix^D=0,1\}
      factor(ix^D)={abs(1-ix^D-xd^D)*}
    {enddo\}

    ! Bfield and Te for interpolation
    {do ix^DB=0,1\}
      do j=1,ndir
        Bnear(ix^D,j)=w(ixbl^D+ix^D,mag(j))
      enddo
      Tnear(ix^D)=w(ixbl^D+ix^D,p_)/w(ixbl^D+ix^D,rho_)
    {enddo\}

    Bloc(:)=0.d0
    Tloc=0.d0
    ! interpolation
    {do ix^DB=0,1\}
      do j=1,ndir
        Bloc(j)=Bloc(j)+Bnear(ix^D,j)*factor(ix^D)
      enddo
      Tloc=Tloc+Tnear(ix^D)*factor(ix^D)
    {enddo\}

  end subroutine get_localBT_mask

  subroutine get_x_neighbor(ip,xfnow,Bnow,xfnb,ds,forward)
    integer :: ip,j
    double precision :: xfnow(1:ndim),xfnb(1:ndim),Bnow(1:ndir)
    double precision :: ds,Btotal,dxf(1:ndim),dxftot
    logical :: forward

    Btotal=0.d0
    do j=1,ndim
      Btotal=Btotal+(Bnow(j))**2
    enddo
    Btotal=dsqrt(Btotal)

    ! if local magnetic field is 0
    if(Btotal==0) then
      dxf=0.d0
      if(forward .eqv. .TRUE.) then
        dxf(^ND)=ds
      else
        dxf(^ND)=-ds
      endif
    else
    ! find next point based on magnetic field direction
      if(forward .eqv. .TRUE.) then
        do j=1,ndim
          dxf(j)=ds*Bnow(j)/Btotal
        enddo
      else
        do j=1,ndim
          dxf(j)=-ds*Bnow(j)/Btotal
        enddo
      endif
    endif

    ! next point
    do j=1,ndim
      xfnb(j)=xfnow(j)+dxf(j)
    enddo

  end subroutine get_x_neighbor

  subroutine interp_Tcoff_mask(xF,ipel,igridl,numR,Tlcoff)
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
            ps(igrid)%w(ixc^D,Tweight_)=ps(igrid)%w(ixc^D,Tweight_)+weight
            ps(igrid)%w(ixc^D,Tcoff_)=ps(igrid)%w(ixc^D,Tcoff_)+weight*Tcnow
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
                  ps(igrid_nb)%w(ixc^DD,Tweight_)=ps(igrid_nb)%w(ixc^DD,Tweight_)+weight
                  ps(igrid_nb)%w(ixc^DD,Tcoff_)=ps(igrid_nb)%w(ixc^DD,Tcoff_)+weight*Tcnow
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
                    ps(igrid_nb)%w(ixc^DD,Tweight_)=ps(igrid_nb)%w(ixc^DD,Tweight_)+weight
                    ps(igrid_nb)%w(ixc^DD,Tcoff_)=ps(igrid_nb)%w(ixc^DD,Tcoff_)+weight*Tcnow
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
      where(ps(igrid)%w(ixO^S,Tweight_)>0.d0)
        ps(igrid)%w(ixO^S,Tcoff_)=ps(igrid)%w(ixO^S,Tcoff_)/ps(igrid)%w(ixO^S,Tweight_)
      elsewhere
        ps(igrid)%w(ixO^S,Tcoff_)=T_bott
      endwhere
    enddo

  end subroutine interp_Tcoff_mask

  subroutine convert_data_back_trac()
    integer :: ixI^L,ixO^L
    integer :: igrid,j,idir

    ^D&ixImin^D=ixglo^D;
    ^D&ixImax^D=ixghi^D;
    ^D&ixOmin^D=ixmlo^D;
    ^D&ixOmax^D=ixmhi^D;

    do j=1,ngrid_trac
      igrid=trac_grid(j)
      if (B0field) then
        do idir=1,ndir
          ps(igrid)%w(ixI^S,mag(idir))=ps(igrid)%w(ixI^S,mag(idir))-ps(igrid)%B0(ixI^S,idir,0)
        enddo
      endif
      call mhd_to_conserved(ixI^L,ixI^L,ps(igrid)%w,ps(igrid)%x)
    enddo

  endsubroutine convert_data_back_trac

end module mod_trac
