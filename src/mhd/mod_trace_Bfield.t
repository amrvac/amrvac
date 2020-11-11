!> mod_trace_Bfield -- trace magnetic field line
!>
!> This module should be used in a global subroutine, since it will 
!> search blocks/grids controlled by different processors.
!>
!> We can use subroutine trace_Bfield to trace a magnetic field line.
!> The location of a point in this field line need to be provided at
!> xf(1,:), and then other points will be calculated. The locations and
!> the plasma parameter of the points will be return. 
!>
!> wB(:,rho_): density
!> wB(:,mom(:)): velocity
!> wB(:,p_): pressure
!> wB(:,mag(:)): magnetic field
!> wB(:,nw+1:nw+ndir): current
!> dL: distance between two points in the same field line
!> numP: the number of points we wants to return.
!> numRT: how many points are valid.
!> forward=true: tracing B field parallel
!> forward=false: tracing B field antiparallel
!> wRT(i)=true: calculate and return wB(:,i)
module mod_trace_Bfield
  use mod_mhd
  implicit none

contains

  subroutine trace_Bfield(xf,wB,dL,numP,numRT,forward,wRT)
    ! trace a field line
    use mod_usr_methods
    use mod_global_parameters
    use mod_particle_base

    integer :: numP,numRT
    double precision :: xf(numP,ndim),wB(numP,nw+ndir)
    double precision :: dL
    logical :: forward
    logical :: wRT(nw+ndir)

    double precision :: dxb^D,xb^L
    integer :: indomain
    integer :: igrid,igrid_now,igrid_next,j
    integer :: ipe_now,ipe_next
    double precision :: xp_in(ndim),xp_out(ndim),x3d(3)
    integer :: ipoint_in,ipoint_out
    double precision :: statusB(4+ndim)
    integer :: pe_record(numP)
    logical :: stopB
    !double precision :: data_share(numP,ndim+nw+ndir)
    double precision :: data_send(numP),data_recv(numP)
    integer :: nums1,nums2,ipoint,nwRT,iRT

    wB=0
    xf(2:numP,:)=0
    pe_record=-1


    ! check whether or the first point is inside simulation box. if yes, find
    ! the pe and igrid for the point
    indomain=0
    numRT=0
    {if (xf(1,^DB)>=xprobmin^DB .and. xf(1,^DB)<xprobmax^DB) indomain=indomain+1\}
    if (indomain==ndim) then
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
        call find_points_in_pe(igrid,ipoint_in,xf,wB,numP,dL,forward,wRT,statusB)
      endif

      ! comunication
      call MPI_BCAST(statusB,4+ndim,MPI_DOUBLE_PRECISION,ipe_now,icomm,ierrmpi)

      ! prepare for next step
      ipoint_out=int(statusB(1))
      ipe_next=int(statusB(2))
      igrid_next=int(statusB(3))
      if (int(statusB(4))==1) then
        stopB=.TRUE.
        numRT=ipoint_out-1
      endif
      do j=1,ndim
        xf(ipoint_out,j)=statusB(4+j)
      enddo

      ! pe and grid of next point
      do j=ipoint_in,ipoint_out-1
        pe_record(j)=ipe_now
      enddo
      ipe_now=ipe_next
      igrid_now=igrid_next
      ipoint_in=ipoint_out
    enddo

    

    ! comunications between processors
    do j=1,numRT
      if (mype/=pe_record(j)) then
        xf(j,:)=0.d0
        wB(j,:)=0.d0
      endif
    enddo

    do j=1,ndim
      data_send=xf(:,j)
      call MPI_ALLREDUCE(data_send,data_recv,numP,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,icomm,ierrmpi)
      xf(:,j)=data_recv
    enddo

    do j=1,nw
      if (wRT(j)) then
        data_send=wB(:,j)
        call MPI_ALLREDUCE(data_send,data_recv,numP,MPI_DOUBLE_PRECISION,&
                             MPI_SUM,icomm,ierrmpi)
        wB(:,j)=data_recv
      else
        wB(:,j)=0.d0
      endif
    enddo

  end subroutine trace_Bfield

  subroutine find_points_in_pe(igrid,ipoint_in,xf,wB,numP,dL,forward,wRT,statusB)
 
    integer :: igrid,ipoint_in,numP
    double precision :: xf(numP,ndim),wB(numP,nw+ndir)
    double precision :: dL
    logical :: forward
    logical :: wRT(nw+ndir)
    double precision :: statusB(4+ndim)    

    integer :: ipe_next,igrid_next,ip_in,ip_out,j
    logical :: newpe,stopB
    double precision :: xfout(ndim) 
 
    ip_in=ipoint_in
    newpe=.FALSE.


    do while(newpe .eqv. .FALSE.)

      ! looking for points in given grid    
      call find_points_interp(igrid,ip_in,ip_out,xf,wB,numP,dL,forward,wRT)

      ip_in=ip_out

      ! when next point is out of given grid, find next grid  
      if (ip_out<numP) then
        stopB=.FALSE.
        xfout=xf(ip_out,:)
        call find_next_grid(igrid,igrid_next,ipe_next,xfout,newpe,stopB)
      else
        newpe=.TRUE.
        stopB=.TRUE.
      endif

      if (newpe) then
        statusB(1)=ip_out
        statusB(2)=ipe_next
        statusB(3)=igrid_next
        statusB(4)=0
        if (stopB) statusB(4)=1
        do j=1,ndim
          statusB(4+j)=xf(ip_out,j)
        enddo
      endif
  
      if (newpe .eqv. .FALSE.) igrid=igrid_next
    enddo

  end subroutine find_points_in_pe

  subroutine find_points_interp(igrid,ip_in,ip_out,xf,wB,numP,dL,forward,wRT)
    use mod_global_parameters

    integer :: igrid,ip_in,ip_out,numP
    double precision :: xf(numP,ndim),wB(numP,nw+ndir)
    double precision :: dL
    logical :: forward
    logical :: wRT(nw+ndir)

    double precision :: dxb^D,xb^L
    integer          :: ip,inblock,ixI^L,j
    double precision :: Bg(ixg^T,ndir)
    double precision :: xfnow(ndim),wBnow(nw+ndir)
    double precision :: xs1(ndim),xs2(ndim),K1(ndim),K2(ndim)

    ^D&ixImin^D=ixglo^D\
    ^D&ixImax^D=ixghi^D\

    ^D&dxb^D=rnode(rpdx^D_,igrid)\
    ^D&xbmin^D=rnode(rpxmin^D_,igrid)\
    ^D&xbmax^D=rnode(rpxmax^D_,igrid)\

    if (B0field) then
      do j=1,ndir
        Bg(ixI^S,j)=ps(igrid)%w(ixI^S,mag(j))+ps(igrid)%B0(ixI^S,j,0)
      enddo
    else
      do j=1,ndir
        Bg(ixI^S,j)=ps(igrid)%w(ixI^S,mag(j))
      enddo
    endif


    ! main loop
    MAINLOOP: do ip=ip_in,numP-1

      ! get local values for variable via interpolation
      xfnow(:)=xf(ip,:)
      call get_w_local(xfnow,wBnow,ps(igrid)%x,ps(igrid)%w,Bg,ixI^L,dxb^D,wRT)
      wB(ip,:)=wBnow(:)

      ! integrate magnetic field with Runge-Kutta method
      xs1(:)=xf(ip,:)
      call get_K(xs1,ps(igrid)%x,Bg,K1,ixI^L,dxb^D)
      if (forward) then
        xs2(:)=xf(ip,:)+dL*K1(:)
      else
        xs2(:)=xf(ip,:)-dL*K1(:)
      endif
      call get_K(xs2,ps(igrid)%x,Bg,K2,ixI^L,dxb^D)
      if (forward) then
        xf(ip+1,:)=xf(ip,:)+dL*(0.5*K1(:)+0.5*K2(:))
      else
        xf(ip+1,:)=xf(ip,:)-dL*(0.5*K1(:)+0.5*K2(:))
      endif

      ip_out=ip+1

      ! whether or not next point is in this block/grid
      inblock=0
      {if (xf(ip+1,^DB)>=xbmin^DB .and. xf(ip+1,^DB)<xbmax^DB) inblock=inblock+1\}
      if (inblock/=ndim) then
        ! exit loop if next point is not in this block
        exit MAINLOOP
      endif

    enddo MAINLOOP

  end subroutine find_points_interp

  subroutine get_w_local(xfn,wBn,x,w,Bg,ixI^L,dxb^D,wRT)

    integer :: ixI^L
    double precision :: dxb^D
    double precision :: xfn(ndim),wBn(nw+ndir)
    double precision :: x(ixI^S,ndim),w(ixI^S,nw),Bg(ixI^S,ndir)
    logical :: wRT(nw+ndir)

    integer          :: ixb^D,ix^D,ixbl^D,j
    double precision :: xd^D
    double precision :: factor(0:1^D&)
    double precision :: wBnear(0:1^D&,nw+ndir)

    integer :: idirmin,idir,jdir,kdir
    integer           :: hxO^L,jxO^L,nxO^L
    double precision :: tmp(0:1^D&)
    double precision :: ek,eb
    double precision :: dxb(ndim)

    ^D&ixbl^D=floor((xfn(^D)-x(ixImin^DD,^D))/dxb^D)+ixImin^D\
    ^D&xd^D=(xfn(^D)-x(ixbl^DD,^D))/dxb^D\
    ^D&dxb(^D)=dxb^D\

    {do ix^D=0,1\}
      factor(ix^D)={abs(1-ix^D-xd^D)*}
    {enddo\}

    {do ix^DB=0,1\}
      ! Bfield for interpolation
      do j=1,ndir
        wBnear(ix^D,mag(j))=Bg(ixbl^D+ix^D,j)
      enddo

      ! density for interpolation
      if (wRT(rho_)) then
        wBnear(ix^D,rho_)=w(ixbl^D+ix^D,rho_)
      endif

      ! velocity for interpolation
      do j=1,ndir
        if (wRT(mom(j))) then
          wBnear(ix^D,mom(j))=w(ixbl^D+ix^D,mom(j))/w(ixbl^D+ix^D,rho_)
        endif
      enddo

      ! pressure for interpolation
      if (wRT(p_)) then
        wBnear(ix^D,p_)=w(ixbl^D+ix^D,e_)
        do j=1,ndir
          ek=0.5d0*(w(ixbl^D+ix^D,mom(j)))**2/w(ixbl^D+ix^D,rho_)
          eb=0.5d0*w(ixbl^D+ix^D,mag(j))**2
          wBnear(ix^D,p_)=wBnear(ix^D,p_)-ek-eb
        enddo
        wBnear(ix^D,p_)=wBnear(ix^D,p_)*(mhd_gamma-1)
      endif
    {enddo\}

    ! current for interpolation
    if (wRT(nw+1) .or. wRT(nw+2) .or. wRT(nw+3)) then
      ^D&nxOmin^D=0\
      ^D&nxOmax^D=1\
      idirmin=7-2*ndir
      do idir=idirmin,3; do jdir=1,ndim; do kdir=1,ndir
        if (lvc(idir,jdir,kdir)/=0) then
          hxO^L=ixbl^D+nxO^L-kr(jdir,^D);
          jxO^L=ixbl^D+nxO^L+kr(jdir,^D);
          tmp(nxO^S)=half*(Bg(jxO^S,kdir)-Bg(hxO^S,kdir))/dxb(jdir)

          if (lvc(idir,jdir,kdir)==1) then
            wBnear(nxO^S,nw+idir)=wBnear(nxO^S,nw+idir)+tmp(nxO^S)
          else
            wBnear(nxO^S,nw+idir)=wBnear(nxO^S,nw+idir)-tmp(nxO^S)
          endif
        endif
      enddo; enddo; enddo
    endif

    wBn=0.d0
    {do ix^DB=0,1\}
      do j=1,nw+ndir
        wBn(j)=wBn(j)+wBnear(ix^D,j)*factor(ix^D)
      enddo
    {enddo\}

  end subroutine get_w_local

  subroutine get_K(xfn,x,Bg,K,ixI^L,dxb^D)

    integer :: ixI^L
    double precision :: dxb^D
    double precision :: xfn(ndim),K(ndim)
    double precision :: x(ixI^S,ndim),Bg(ixI^S,ndir)

    integer          :: ixb^D,ix^D,ixbl^D,j
    double precision :: xd^D
    double precision :: Bx(ndim),factor(0:1^D&)
    double precision :: Btotal

    ^D&ixbl^D=floor((xfn(^D)-x(ixImin^DD,^D))/dxb^D)+ixImin^D\
    ^D&xd^D=(xfn(^D)-x(ixbl^DD,^D))/dxb^D\

    {do ix^D=0,1\}
      factor(ix^D)={abs(1-ix^D-xd^D)*}
    {enddo\}
    Bx=0.d0
    {do ix^DB=0,1\}
      do j=1,ndim
        Bx(j)=Bx(j)+Bg(ixbl^D+ix^D,j)*factor(ix^D)
      enddo
    {enddo\}

    do j=1,ndim
      Btotal=Btotal+(Bx(j))**2
    enddo
    Btotal=dsqrt(Btotal)

    if (Btotal==0.d0) then
      K=0.d0
      K(1)=1.d0
    else
      K(1:ndim)=Bx(1:ndim)/Btotal
    endif

  end subroutine get_K

  subroutine find_next_grid(igrid,igrid_next,ipe_next,xf1,newpe,stopB)
    ! check the grid and pe of next point
    use mod_usr_methods
    use mod_global_parameters
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
    case (neighbor_boundary)
      ! next point is not in simulation box
      newpe=.TRUE.
      stopB=.TRUE.

    case(neighbor_coarse)
      ! neighbor grid has lower refinement level      
      igrid_next=igrid_neighbor
      ipe_next=ipe_neighbor
      if (mype==ipe_neighbor) then
        newpe=.FALSE.
      else
        newpe=.TRUE.
      endif

    case(neighbor_sibling)
      ! neighbor grid has lower refinement level 
      igrid_next=igrid_neighbor
      ipe_next=ipe_neighbor
      if (mype==ipe_neighbor) then
        newpe=.FALSE.
      else
        newpe=.TRUE.
      endif

    case(neighbor_fine)
      ! neighbor grid has higher refinement level 
      {xbmid^D=(xbmin^D+xbmax^D)/2.d0\}
      ^D&inc^D=1\ 
      {if (xf1(^D)<=xbmin^D) inc^D=0\}
      {if (xf1(^D)>xbmin^D .and. xf1(^D)<=xbmid^D) inc^D=1\}
      {if (xf1(^D)>xbmid^D .and. xf1(^D)<xbmax^D) inc^D=2\}
      {if (xf1(^D)>=xbmax^D) inc^D=3\}
      ipe_next=neighbor_child(2,inc^D,igrid)
      igrid_next=neighbor_child(1,inc^D,igrid)
      if (mype==ipe_next) then
        newpe=.FALSE.
      else
        newpe=.TRUE.
      endif

    end select

  end subroutine find_next_grid
  
end module mod_trace_Bfield
