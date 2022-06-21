!> mod_trace_Vfield -- trace flow streamline
!>
!> This module should be used in a global subroutine, since it will 
!> search blocks/grids controlled by different processors.
!>
!> We can use subroutine trace_Vfield to trace a flow streamline.
!> The start location of a point on this streamline needs to be provided by
!> xf(1,1:ndim), and then other points on it will be calculated. The locations and
!> the plasma parameters at these points will be returned. 
!>
!> wpV(:,1:nw): stores the primitive variables at the streamline points
!>
!> use phys_to_conserved, phys_to_primitive and phys_get_v_idim 
!> to get local primitive values as well as velocity
!>
!> dL: distance between two points on the same streamline
!> numP: the number of points on a streamline we want to return.
!> numRT: how many points are valid.
!> forward=true: tracing V field parallel
!> forward=false: tracing V field antiparallel
!>
module mod_trace_Vfield
  implicit none

contains

  subroutine trace_Vfield(xf,wpV,dL,numP,numRT,forward)
    ! trace a streamline
    use mod_usr_methods
    use mod_global_parameters
    use mod_particle_base

    integer :: numP,numRT
    double precision :: xf(numP,ndim),wpV(numP,nw)
    double precision :: dL
    logical :: forward

    integer :: indomain
    integer :: igrid,igrid_now,igrid_next,j
    integer :: ipe_now,ipe_next
    double precision :: xp_in(ndim),xp_out(ndim),x3d(3)
    integer :: ipoint_in,ipoint_out
    double precision :: statusB(4+ndim)
    integer :: pe_record(numP)
    logical :: stopB
    double precision :: data_send(numP),data_recv(numP)

    wpV=0.0d0
    xf(2:numP,:)=0.0d0
    pe_record=-1

    ! check whether the given first point is inside simulation box. if yes, find
    ! the pe and igrid for this point
    indomain=0
    numRT=0
    {if (xf(1,^DB)>=xprobmin^DB .and. xf(1,^DB)<xprobmax^DB) indomain=indomain+1\}
    if (indomain==ndim) then
      numRT=1
      ! find pe and igrid
      x3d(1:ndim)=xf(1,1:ndim)
      call find_particle_ipe(x3d,igrid_now,ipe_now)
      stopB=.FALSE.
      ipoint_in=1
    else
      if (mype==0) then
        call MPISTOP('streamline tracing error: given point is not in simulation box!')
      endif
    endif

    ! other points on streamline
    do while(stopB .eqv. .FALSE.)

      if (mype==ipe_now) then
        igrid=igrid_now
        ! looking for points in one pe
        call find_points_in_pe(igrid,ipoint_in,xf,wpV,numP,dL,forward,statusB)
      endif

      ! communication
      call MPI_BCAST(statusB,4+ndim,MPI_DOUBLE_PRECISION,ipe_now,icomm,ierrmpi)

      ! prepare for next step
      ipoint_out=int(statusB(1))
      ipe_next  =int(statusB(2))
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

    ! communications between processors
    do j=1,numRT
      if (mype/=pe_record(j)) then
        xf(j,:)=0.d0
        wpV(j,:)=0.d0
      endif
    enddo

    do j=1,ndim
      data_send=xf(:,j)
      call MPI_ALLREDUCE(data_send,data_recv,numP,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,icomm,ierrmpi)
      xf(:,j)=data_recv
    enddo

    do j=1,nw
      data_send=wpV(:,j)
      call MPI_ALLREDUCE(data_send,data_recv,numP,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,icomm,ierrmpi)
      wpV(:,j)=data_recv
    enddo

  end subroutine trace_Vfield

  subroutine find_points_in_pe(igrid,ipoint_in,xf,wpV,numP,dL,forward,statusB)
    use mod_global_parameters
    integer :: igrid,ipoint_in,numP
    double precision :: xf(numP,ndim),wpV(numP,nw)
    double precision :: dL
    logical :: forward
    double precision :: statusB(4+ndim)    

    integer :: ipe_next,igrid_next,ip_in,ip_out,j
    logical :: newpe,stopB
    double precision :: xfout(ndim) 
 
    ip_in=ipoint_in
    newpe=.FALSE.

    do while(newpe .eqv. .FALSE.)

      ! looking for points in given grid    
      call find_points_interp(igrid,ip_in,ip_out,xf,wpV,numP,dL,forward)

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

  subroutine find_points_interp(igrid,ip_in,ip_out,xf,wpV,numP,dL,forward)
    use mod_global_parameters
    use mod_physics, only: phys_get_v

    integer :: igrid,ip_in,ip_out,numP
    double precision :: xf(numP,ndim),wpV(numP,nw)
    double precision :: dL
    logical :: forward

    double precision :: dxb^D,xb^L
    integer          :: ip,inblock,ixI^L,j
    double precision :: Vg(ixG^T,ndir)
    double precision :: xfnow(ndim),wpVnow(nw)
    double precision :: xs1(ndim),xs2(ndim),K1(ndim),K2(ndim)

    ^D&ixImin^D=ixGlo^D\
    ^D&ixImax^D=ixGhi^D\

    ^D&dxb^D=rnode(rpdx^D_,igrid)\
    ^D&xbmin^D=rnode(rpxmin^D_,igrid)\
    ^D&xbmax^D=rnode(rpxmax^D_,igrid)\

    ! get flow field
    call phys_get_v(ps(igrid)%w,ps(igrid)%x,ixI^L,ixI^L,Vg)

    ! main loop
    MAINLOOP: do ip=ip_in,numP-1

      ! get local values for variables via interpolation
      xfnow(:)=xf(ip,:)
      call get_wprim_local(xfnow,wpVnow,ps(igrid)%x,ps(igrid)%w,ixI^L,dxb^D)
      wpV(ip,:)=wpVnow(:)

      ! integrate flow field with Runge-Kutta method
      xs1(:)=xf(ip,:)
      call get_K(xs1,ps(igrid)%x,Vg,K1,ixI^L,dxb^D)
      if (forward) then
        xs2(:)=xf(ip,:)+dL*K1(:)
      else
        xs2(:)=xf(ip,:)-dL*K1(:)
      endif
      call get_K(xs2,ps(igrid)%x,Vg,K2,ixI^L,dxb^D)
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

  subroutine get_wprim_local(xfn,wpVn,x,w,ixI^L,dxb^D)
    use mod_global_parameters
    use mod_physics, only: phys_to_primitive, phys_to_conserved

    integer :: ixI^L
    double precision :: dxb^D
    double precision :: xfn(ndim),wpVn(nw)
    double precision :: x(ixI^S,ndim),w(ixI^S,nw)

    integer          :: ixb^D,ix^D,ixbl^D,j,iw,ixO^L
    double precision :: xd^D
    double precision :: factor(0:1^D&)
    double precision :: wpVnear(0:1^D&,nw)
    double precision :: dxb(ndim)

    ^D&ixbl^D=floor((xfn(^D)-x(ixImin^DD,^D))/dxb^D)+ixImin^D\
    ^D&xd^D=(xfn(^D)-x(ixbl^DD,^D))/dxb^D\
    ^D&dxb(^D)=dxb^D\

    {do ix^D=0,1\}
      factor(ix^D)={abs(1-ix^D-xd^D)*}
    {enddo\}

    ^D&ixOmin^D=ixbl^D;
    ^D&ixOmax^D=ixbl^D+1;
    call phys_to_primitive(ixI^L,ixO^L,w,x)
    {do ix^DB=0,1\}
        wpVnear(ix^D,1:nw)=w(ixbl^D+ix^D,1:nw)
    {enddo\}
    call phys_to_conserved(ixI^L,ixO^L,w,x)

    wpVn=0.d0
    {do ix^DB=0,1\}
      do j=1,nw
        wpVn(j)=wpVn(j)+wpVnear(ix^D,j)*factor(ix^D)
      enddo
    {enddo\}

  end subroutine get_wprim_local

  subroutine get_K(xfn,x,Vg,K,ixI^L,dxb^D)
    use mod_global_parameters
    ! Note: here only to ndim: what for ndir>ndim?

    integer :: ixI^L
    double precision :: dxb^D
    double precision :: xfn(ndim),K(ndim)
    double precision :: x(ixI^S,ndim),Vg(ixI^S,ndir)

    integer          :: ixb^D,ix^D,ixbl^D,j
    double precision :: xd^D
    double precision :: Vx(ndim),factor(0:1^D&)
    double precision :: Vtotal

    ^D&ixbl^D=floor((xfn(^D)-x(ixImin^DD,^D))/dxb^D)+ixImin^D\
    ^D&xd^D=(xfn(^D)-x(ixbl^DD,^D))/dxb^D\

    {do ix^D=0,1\}
      factor(ix^D)={abs(1-ix^D-xd^D)*}
    {enddo\}

    Vx=0.d0
    {do ix^DB=0,1\}
      do j=1,ndim
        Vx(j)=Vx(j)+Vg(ixbl^D+ix^D,j)*factor(ix^D)
      enddo
    {enddo\}

    Vtotal=0.d0
    do j=1,ndim
      Vtotal=Vtotal+(Vx(j))**2
    enddo
    Vtotal=dsqrt(Vtotal)

    if (Vtotal==0.d0) then
      K=0.d0
      K(1)=1.d0
    else
      K(1:ndim)=Vx(1:ndim)/Vtotal
    endif

  end subroutine get_K

  subroutine find_next_grid(igrid,igrid_next,ipe_next,xf1,newpe,stopB)
    ! check the grid and pe of next point
    use mod_usr_methods
    use mod_global_parameters
    use mod_forest

    integer :: igrid,igrid_next,ipe_next
    double precision :: xf1(ndim)
    logical :: newpe,stopB

    integer :: idn^D,my_neighbor_type
    integer :: inc^D,ipe_neighbor,igrid_neighbor
    double precision :: xb^L,xbmid^D

    ^D&xbmin^D=rnode(rpxmin^D_,igrid)\
    ^D&xbmax^D=rnode(rpxmax^D_,igrid)\

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
      ! neighbor grid has same refinement level 
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
  
end module mod_trace_Vfield
