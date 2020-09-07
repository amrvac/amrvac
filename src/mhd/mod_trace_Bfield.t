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
!> interp=true: using interpolation in the calculation
module mod_trace_Bfield
  use mod_mhd
  implicit none

contains

  subroutine trace_Bfield(xf,wB,dL,numP,numRT,forward,wRT,interp)
    ! trace a field line
    use mod_usr_methods
    use mod_global_parameters
    use mod_particle_base

    integer :: numP,numRT
    double precision :: xf(numP,ndim),wB(numP,nw+ndir)
    double precision :: dL
    logical :: forward,interp
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
        call find_points_in_pe(igrid,ipoint_in,xf,wB,numP,dL,forward,wRT,statusB,interp)
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

  subroutine find_points_in_pe(igrid,ipoint_in,xf,wB,numP,dL,forward,wRT,statusB,interp)
 
    integer :: igrid,ipoint_in,numP
    double precision :: xf(numP,ndim),wB(numP,nw+ndir)
    double precision :: dL
    logical :: forward,interp
    logical :: wRT(nw+ndir)
    double precision :: statusB(4+ndim)    

    integer :: ipe_next,igrid_next,ip_in,ip_out,j
    logical :: newpe,stopB
    double precision :: xfout(ndim) 
 
    ip_in=ipoint_in
    newpe=.FALSE.


    do while(newpe .eqv. .FALSE.)

      ! looking for points in given grid    
      if (interp) then
        call find_points_interp(igrid,ip_in,ip_out,xf,wB,numP,dL,forward,wRT)
      else
        call find_points_nointerp(igrid,ip_in,ip_out,xf,wB,numP,dL,forward,wRT)
      endif

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

    integer :: igrid,ip_in,ip_out,numP
    double precision :: xf(numP,ndim),wB(numP,nw+ndir)
    double precision :: dL
    logical :: forward
    logical :: wRT(nw+ndir)

    integer          :: ixO^L,ixO^D,j
    double precision :: dxf(ndim)
    double precision :: dxb^D,xb^L,xd^D
    integer          :: ixb^D,ix^D,ixbl^D,ip
    double precision :: wBnear(0:1^D&,nw+ndir)
    double precision :: Bx(ndim),factor(0:1^D&)
    double precision :: Btotal,maxft,Bp

    integer :: idirmin,ixI^L,idir,jdir,kdir
    double precision :: dxb(ndim)
    integer :: hxO^L,jxO^L,nxO^L
    double precision :: tmp(0:1^D&)
    double precision :: ek,eb

    integer :: inblock


    ^D&ixImin^D=ixglo^D\
    ^D&ixImax^D=ixghi^D\
    ^D&ixOmin^D=ixmlo^D\
    ^D&ixOmax^D=ixmhi^D\

    ^D&dxb^D=rnode(rpdx^D_,igrid)\
    ^D&xbmin^D=rnode(rpxmin^D_,igrid)\
    ^D&xbmax^D=rnode(rpxmax^D_,igrid)\
    ^D&dxb(^D)=dxb^D\

    B0field=.TRUE.

    ! main loop
    MAINLOOP: do ip=ip_in,numP-1

      ! do interpolation
      ^D&ixbl^D=floor((xf(ip,^D)-ps(igrid)%x(ixOmin^DD,^D))/dxb^D)+ixOmin^D\
      ^D&xd^D=(xf(ip,^D)-ps(igrid)%x(ixbl^DD,^D))/dxb^D\
      wBnear=0

      {do ix^DB=0,1\}
        ! Bfield for interpolation

        do j=1,ndir
          if (B0field) then
            wBnear(ix^D,mag(j))=ps(igrid)%w(ixbl^D+ix^D,mag(j))+&
                                     ps(igrid)%B0(ixbl^D+ix^D,j,0)
          else
            wBnear(ix^D,mag(j))=ps(igrid)%w(ixbl^D+ix^D,mag(j))
          endif
        enddo

        ! density for interpolation
        if (wRT(rho_)) then
          wBnear(ix^D,rho_)=ps(igrid)%w(ixbl^D+ix^D,rho_)
        endif

        ! velocity for interpolation
        do j=1,ndir
          if (wRT(mom(j))) then
            wBnear(ix^D,mom(j))=ps(igrid)%w(ixbl^D+ix^D,mom(j))/&
                                ps(igrid)%w(ixbl^D+ix^D,rho_)
          endif
        enddo

        ! pressure for interpolation
        if (wRT(p_)) then
          wBnear(ix^D,p_)=ps(igrid)%w(ixbl^D+ix^D,e_)
          do j=1,ndir
            ek=0.5d0*(ps(igrid)%w(ixbl^D+ix^D,mom(j)))**2/&
               ps(igrid)%w(ixbl^D+ix^D,rho_)
            eb=0.5d0*ps(igrid)%w(ixbl^D+ix^D,mag(j))**2
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
            tmp(nxO^S)=half*(ps(igrid)%w(jxO^S,mag(kdir))-&
                             ps(igrid)%w(hxO^S,mag(kdir)))/dxb(jdir)
            if (B0field) then
              tmp(nxO^S)=tmp(nxO^S)+half*(ps(igrid)%B0(jxO^S,kdir,0)-&
                                    ps(igrid)%B0(hxO^S,kdir,0))/dxb(jdir)
            endif

            if (lvc(idir,jdir,kdir)==1) then
              wBnear(nxO^S,nw+idir)=wBnear(nxO^S,nw+idir)+tmp(nxO^S)
            else
              wBnear(nxO^S,nw+idir)=wBnear(nxO^S,nw+idir)-tmp(nxO^S)
            endif
          endif
        enddo; enddo; enddo
      endif


      ! interpolation factor
      {do ix^D=0,1\}
        factor(ix^D)={abs(1-ix^D-xd^D)*}
      {enddo\}

      ! do interpolation to get local magnetic field
      Bx=0
      wB(ip,:)=0
      {do ix^DB=0,1\}
        {Bx(^DB)=Bx(^DB)+wBnear(ix^DD,mag(^DB))*factor(ix^DD)\}
        do j=1,nw+ndir
          wB(ip,j)=wB(ip,j)+wBnear(ix^D,j)*factor(ix^D)
        enddo
      {enddo\}


      ! local magnetic field strength
      Btotal=0.0d0
      do j=1,ndim
        Btotal=Btotal+(Bx(j))**2
      enddo
      Btotal=dsqrt(Btotal)

      ! if local magnetic field is 0
      if (Btotal==0) then
        maxft=factor(0^D&)
        {do ix^DB=0,1\}
          ! local B equls the B of the closest point
          Bp=0
          do j=1,ndim
            Bp=Bp+(wBnear(ix^D,mag(j)))**2
          enddo

          if (factor(ix^D)>=maxft .and. Bp/=0) then
            Bx(:)=wBnear(ix^D,mag(:))
          endif
          Btotal=Bp
        {enddo\}

        ! all the point near hear is 0
        if (Btotal==0) then
          Bx(:)=1
          Btotal=dsqrt(1.0d0*ndim)
        endif
      endif


      ! find next point based on magnetic field direction
      if (forward .eqv. .TRUE.) then
        do j=1,ndim
          dxf(j)=dL*Bx(j)/Btotal
        enddo
      else
        do j=1,ndim
          dxf(j)=-dL*Bx(j)/Btotal
        enddo
      endif

      ! next point
      do j=1,ndim
        xf(ip+1,j)=xf(ip,j)+dxf(j)
      enddo
      ip_out=ip+1

      !print *, ip,xf(ip,:),Bx

      ! whether or not next point is in this block/grid
      inblock=0
      {if (xf(ip+1,^DB)>=xbmin^DB .and. xf(ip+1,^DB)<xbmax^DB) inblock=inblock+1\}
      if (inblock/=ndim) then
        ! exit loop if next point is not in this block
        exit MAINLOOP
      endif

    enddo MAINLOOP

  end subroutine find_points_interp

  subroutine find_points_nointerp(igrid,ip_in,ip_out,xf,wB,numP,dL,forward,wRT)

    integer :: igrid,ip_in,ip_out,numP
    double precision :: xf(numP,ndim),wB(numP,nw+ndir)
    double precision :: dL
    logical :: forward
    logical :: wRT(nw+ndir)

    double precision :: Bg(ixg^T,ndir)
    integer          :: ixO^L,ixO^D,ixb^D,ix^D,j,ip
    double precision :: dxf(ndim),Bx(ndim),dxb(ndim)
    double precision :: dxb^D,xb^L,xd^D
    integer :: hxO^D,jxO^D
    integer :: idirmin,ixI^L,idir,jdir,kdir
    double precision :: ek,eb,Btotal,tmp
    integer :: inblock

    ^D&ixImin^D=ixglo^D\
    ^D&ixImax^D=ixghi^D\
    ^D&ixOmin^D=ixmlo^D\
    ^D&ixOmax^D=ixmhi^D\

    ^D&dxb^D=rnode(rpdx^D_,igrid)\
    ^D&xbmin^D=rnode(rpxmin^D_,igrid)\
    ^D&xbmax^D=rnode(rpxmax^D_,igrid)\
    ^D&dxb(^D)=dxb^D\


    ! main loop
    MAINLOOP: do ip=ip_in,numP-1

      ! looking for the cell index
      ^D&ixb^D=floor((xf(ip,^D)-xbmin^D)/dxb^D)+ixOmin^D\

      ! magnetic field      
      wB(ip,mag(1:ndir))=ps(igrid)%w(ixb^D,mag(1:ndir))
      if (B0field) then
        wB(ip,mag(1:ndir))=wB(ip,mag(1:ndir))+ps(igrid)%B0(ixb^D,1:ndir,0)
      endif

      ! density
      if (wRT(rho_)) then
        wB(ip,rho_)=ps(igrid)%w(ixb^D,rho_)
      endif

      ! velocity
      do j=1,ndir
        if (wRT(mom(j))) then
          wB(ip,mom(j))=ps(igrid)%w(ixb^D,mom(j))/&
                        ps(igrid)%w(ixb^D,rho_)
        endif
      enddo

      ! pressure
      if (wRT(p_)) then
        wB(ip,p_)=ps(igrid)%w(ixb^D,e_)
        do j=1,ndir
          ek=0.5d0*(ps(igrid)%w(ixb^D,mom(j)))**2/&
             ps(igrid)%w(ixb^D,rho_)
          eb=0.5d0*ps(igrid)%w(ixb^D,mag(j))**2
          wB(ip,p_)=wB(ip,p_)-ek-eb
        enddo
      endif

      ! current
      if (wRT(nw+1) .or. wRT(nw+2) .or. wRT(nw+3)) then
        idirmin=7-2*ndir
        do idir=idirmin,3; do jdir=1,ndim; do kdir=1,ndir
          if (lvc(idir,jdir,kdir)/=0) then
            hxO^D=ixb^D-kr(jdir,^D);
            jxO^D=ixb^D+kr(jdir,^D);
            tmp=half*(ps(igrid)%w(jxO^D,mag(kdir))-&
                      ps(igrid)%w(hxO^D,mag(kdir)))/dxb(jdir)
            if (B0field) then
              tmp=tmp+half*(ps(igrid)%B0(jxO^D,kdir,0)-&
                            ps(igrid)%B0(hxO^D,kdir,0))/dxb(jdir)
            endif

            if (lvc(idir,jdir,kdir)==1) then
              wB(ip,nw+idir)=wB(ip,nw+idir)+tmp
            else
              wB(ip,nw+idir)=wB(ip,nw+idir)-tmp
            endif
          endif
        enddo; enddo; enddo
      endif



      ! local magnetic field strength
      Bx(1:ndim)=wB(ip,mag(1:ndim))
      Btotal=0.0d0
      do j=1,ndim
        Btotal=Btotal+(Bx(j))**2
      enddo
      Btotal=dsqrt(Btotal)

      ! if local magnetic field is 0
      if (Btotal==0) then
        if (ip>1) then
          Bx(1:ndim)=wB(ip,mag(1:ndim))
        else
          Bx=1
        endif

        do j=1,ndim
          Btotal=Btotal+(Bx(j))**2
        enddo
        Btotal=dsqrt(Btotal)
      endif

      ! find next point based on magnetic field direction
      if (forward .eqv. .TRUE.) then
        do j=1,ndim
          dxf(j)=dL*Bx(j)/Btotal
        enddo
      else
        do j=1,ndim
          dxf(j)=-dL*Bx(j)/Btotal
        enddo
      endif

      ! next point
      do j=1,ndim
        xf(ip+1,j)=xf(ip,j)+dxf(j)
      enddo
      ip_out=ip+1

      ! whether or not next point is in this block/grid
      inblock=0
      {if (xf(ip+1,^DB)>=xbmin^DB .and. xf(ip+1,^DB)<xbmax^DB) inblock=inblock+1\}
      if (inblock/=ndim) then
        ! exit loop if next point is not in this block
        exit MAINLOOP
      endif

    enddo MAINLOOP


  end subroutine find_points_nointerp

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
