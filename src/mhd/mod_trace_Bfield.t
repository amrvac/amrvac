module mod_trace_Bfield
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

  use mod_mhd

  contains

  subroutine trace_Bfield(xf,wB,dL,numP,numRT,forward,wRT,interp,ix1)
    ! trace a field line
    use mod_usr_methods
    use mod_global_parameters

    integer :: numP,numRT,ix1
    double precision :: xf(numP,ndim),wB(numP,nw+ndir)
    double precision :: dL
    logical :: forward,interp
    logical :: wRT(nw+ndir)

    double precision :: dxb^D,xb^L
    integer :: inblock,indomain
    integer :: iigrid,igrid,igrid_now,igrid_next,j
    integer :: ipe_now,ipe_next,mainpe
    integer :: status(mpi_status_size)
    double precision :: xp_in(ndim),xp_out(ndim)
    integer :: ipoint_in,ipoint_out
    double precision :: statusB(4+ndim)
    integer :: pe_record(numP),pe_grid(2)
    logical :: stopB
    double precision :: xf_send(numP,ndim),wB_send(numP,nw+ndir)
    double precision :: data_share(numP,ndim+nw+ndir)
    integer :: nums1,nums2,ipoint,nwRT,iRT


    wB=0
    xf(2:numP,:)=0
    pe_record=-1

    ! set processor 0 as main processor
    mainpe=0

    ! check whether or the first point is inside simulation box. if yes, find
    ! the pe and igrid for the point
    indomain=0
    numRT=0
    {if (xf(1,^DB)>=xprobmin^DB .and. xf(1,^DB)<xprobmax^DB) indomain=indomain+1\}
    if (indomain==ndim) then
      numRT=1

      ! find the grid and pe of the first point
      LOOP1: do iigrid=1,igridstail; igrid=igrids(iigrid);
        ^D&xbmin^D=rnode(rpxmin^D_,igrid)\
        ^D&xbmax^D=rnode(rpxmax^D_,igrid)\
        inblock=0
        {if (xf(1,^DB)>=xbmin^DB .and. xf(1,^DB)<xbmax^DB) inblock=inblock+1\}
        if (inblock==ndim) then
          pe_grid(1)=mype
          pe_grid(2)=igrid
          call MPI_SEND(pe_grid,2,MPI_INTEGER,mainpe,0,icomm,ierrmpi)
          exit LOOP1
        endif
      enddo LOOP1

      if (mype==mainpe) then
        call MPI_RECV(pe_grid,2,MPI_INTEGER,MPI_ANY_SOURCE,0,icomm,status,ierrmpi)
      endif
      call MPI_BCAST(pe_grid,2,MPI_INTEGER,mainpe,icomm,ierrmpi)
    endif


    ! looking for points in the field line
    stopB=.FALSE.
    ipoint_in=1
    ipe_now=pe_grid(1)
    igrid_now=pe_grid(2)

    do while(stopB .eqv. .FALSE.)

      if (mype==ipe_now) then
        igrid=igrid_now

        ! looking for points in one pe
        call find_points_in_pe(igrid,ipoint_in,xf,wB,numP,dL,forward,wRT,statusB,interp,ix1)
        
        call MPI_SEND(statusB,4+ndim,MPI_DOUBLE_PRECISION,mainpe,0,icomm,ierrmpi)
      endif

      ! comunication
      if (mype==mainpe) then
        call MPI_RECV(statusB,4+ndim,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,0,icomm,status,ierrmpi)
      endif
      call MPI_BCAST(statusB,4+ndim,MPI_DOUBLE_PRECISION,mainpe,icomm,ierrmpi)

      ! preparing for next step
      ipe_now=int(statusB(1))
      ipoint_out=int(statusB(2))
      igrid_now=int(statusB(3))
      if (int(statusB(4))==1) then
        stopB=.TRUE.
        numRT=ipoint_out-1
      endif
      do j=1,ndim
        xf(ipoint_out,j)=statusB(4+j)
      enddo

      ! check whether or next point is inside simulation box.
      indomain=0
      {if (xf(ipoint_out,^DB)>xprobmin^DB .and. xf(ipoint_out,^DB)<xprobmax^DB) indomain=indomain+1\}
      if (indomain/=ndim) then
        numRT=ipoint_out-1
        stopB=.TRUE.
      endif


      ! next point is in another pe, find out grid and pe numbers
      if (indomain==ndim) then
        call find_grid(ipe_now,igrid_now,ipe_next,igrid_next,xf(ipoint_out,:),ipoint_out)
      endif

      ! pe and grid of next point
      do j=ipoint_in,ipoint_out-1
        pe_record(j)=ipe_now
      enddo
      ipe_now=ipe_next
      igrid_now=igrid_next
      ipoint_in=ipoint_out

    enddo


    ! comunications between processors
    nwRT=0
    do j=1,nw+ndir
      if (wRT(j)) nwRT=nwRT+1
    enddo

    ipoint_in=1
    ipoint_out=1
    do while (pe_record(ipoint_in)>=0 .and. pe_record(ipoint_in)<=npe)
      nums1=0
      nums2=ndim+nwRT
      ipe_now=pe_record(ipoint_in)
      LOOPS: do ipoint=ipoint_in,numP
        if (pe_record(ipoint)/=ipe_now) exit LOOPS
        nums1=nums1+1
      enddo LOOPS
      ipoint_out=ipoint_in+nums1

      if (mype==ipe_now) then
        do j=1,ndim
          data_share(ipoint_in:ipoint_out-1,j)=xf(ipoint_in:ipoint_out-1,j)
        enddo
        do j=1,nw+ndir
          if (wRT(j)) then
            data_share(ipoint_in:ipoint_out-1,j+ndim)=wB(ipoint_in:ipoint_out-1,j)
          endif
        enddo
      endif

      if (nums2>0) then
        call MPI_BCAST(data_share(ipoint_in:ipoint_out-1,1:nwRT+ndim),nums1*nums2,&
                       MPI_DOUBLE_PRECISION,ipe_now,icomm,ierrmpi)
      endif

      ipoint_in=ipoint_out
    enddo

    do j=1,ndim
      xf(:,j)=data_share(:,j)
    enddo
    iRT=0
    do j=1,nw+ndir
      if (wRT(j)) then
        iRT=iRT+1
        wB(:,j)=data_share(:,iRT+ndim)
      endif
    enddo

  end subroutine trace_Bfield

  subroutine find_points_in_pe(igrid,ipoint_in,xf,wB,numP,dL,forward,wRT,statusB,interp,ix1)
 
    integer :: igrid,ipoint_in,numP,ix1
    double precision :: xf(numP,ndim),wB(numP,nw+ndir)
    double precision :: dL
    logical :: forward,interp
    logical :: wRT(nw+ndir)
    double precision :: statusB(4+ndim)    

    integer :: ipe_next,igrid_next,stopB,ip_in,ip_out
    logical :: newpe
  
    ip_in=ipoint_in
    newpe=.FALSE.


    do while(newpe .eqv. .FALSE.)

      ! looking for points in given grid    
      if (interp) then
        call find_points_interp(igrid,ip_in,ip_out,xf,wB,numP,dL,forward,wRT,ix1)
      else
        call find_points_nointerp(igrid,ip_in,ip_out,xf,wB,numP,dL,forward,wRT)
      endif

      ip_in=ip_out

      !if (ix1==519) print *, ip_out,numP

      ! when next point is out of given grid, find next grid  
      call find_next_grid(igrid,igrid_next,ipe_next,xf(ip_out,:),newpe)

      if (ip_out>=numP) newpe=.TRUE.

      if (newpe) then
        statusB(1)=mype
        statusB(2)=ip_out
        statusB(3)=igrid
        statusB(4)=0
        if (ip_out>=numP) statusB(4)=1
        do j=1,ndim
          statusB(4+j)=xf(ip_out,j)
        enddo
      endif
  
      if (newpe .eqv. .FALSE.) igrid=igrid_next
    enddo

  end subroutine find_points_in_pe

  subroutine find_points_interp(igrid,ip_in,ip_out,xf,wB,numP,dL,forward,wRT,iFL)

    integer :: igrid,ip_in,ip_out,numP,iFL
    double precision :: xf(numP,ndim),wB(numP,nw+ndir)
    double precision :: dL
    logical :: forward
    logical :: wRT(nw+ndir)

    integer          :: ixO^L,ixO^D,j
    double precision :: dxf(ndim)
    double precision :: dxb^D,xb^L,xd^D
    integer          :: ixb^D,ix^D,ixbl^D
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
    integer          :: ixO^L,ixO^D,ixb^D,ix^D,j
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
  
  subroutine find_next_grid(igrid,igrid_next,ipe_next,xf1,newpe)
    ! check the grid and pe of next point
    use mod_usr_methods
    use mod_global_parameters
    use mod_forest

    integer :: igrid,igrid_next,ipe_next
    double precision :: xf1(ndim)
    double precision :: dxb^D,xb^L

    integer :: inblock,inblock_n,inblock_nc
    integer :: ix^D
    logical :: newpe

    integer :: igrid_nb,ipe_nb

    newpe=.TRUE.

    ^D&xbmin^D=rnode(rpxmin^D_,igrid)\
    ^D&xbmax^D=rnode(rpxmax^D_,igrid)\
    inblock=0
    {if (xf1(^D)>=xbmin^D .and. xf1(^D)<xbmax^D) inblock=inblock+1\}
    if (inblock==ndim) then
      ! in the same grid with previous point
      igrid_next=igrid
      ipe_next=mype
      newpe=.FALSE.
    else
      ! not in the same grid
      {do ix^D=-1,1,1\}
        ! check neighbor
        igrid_nb=neighbor(1,ix^D,igrid)
        ipe_nb=neighbor(2,ix^D,igrid)
        if (mype==ipe_nb .and. igrid_inuse(igrid_nb,ipe_nb) .and. &
            igrid_nb>=0 .and. igrid_nb<=max_blocks) then
          ^D&xbmin^D=rnode(rpxmin^D_,igrid_nb)\
          ^D&xbmax^D=rnode(rpxmax^D_,igrid_nb)\
          inblock_n=0
          {if (xf1(^D)>=xbmin^D .and. xf1(^D)<xbmax^D) inblock_n=inblock_n+1\}
          if (inblock_n==ndim) then
            ! in neighbor
            igrid_next=igrid_nb
            ipe_next=mype
            newpe=.FALSE.
          endif
        endif

        ! check neighbor_child
        igrid_nb=neighbor_child(1,ix^D,igrid)
        ipe_nb=neighbor_child(2,ix^D,igrid)
        if (mype==ipe_nb .and. igrid_inuse(igrid_nb,ipe_nb) .and. &
            igrid_nb>=0 .and. igrid_nb<=max_blocks) then
          ^D&xbmin^D=rnode(rpxmin^D_,igrid_nb)\
          ^D&xbmax^D=rnode(rpxmax^D_,igrid_nb)\
          inblock_nc=0
          {if (xf1(^D)>=xbmin^D .and. xf1(^D)<xbmax^D) inblock_nc=inblock_nc+1\}
          if (inblock_nc==ndim) then
            ! in neighbor child
            igrid_next=igrid_nb
            ipe_next=mype
            newpe=.FALSE.
          endif
        endif
      {enddo\}
    endif

  end subroutine find_next_grid

  subroutine find_grid(ipe_now,igrid_now,ipe_next,igrid_next,xf1,nj)
    ! find for grid and pe numbers
    use mod_usr_methods
    use mod_global_parameters
    use mod_forest

    integer :: igrid,iigrid,igrid_now,ipe_now,igrid_next,ipe_next
    double precision :: xf1(ndim)

    double precision :: dxb^D,xb^L
    integer :: ix^D,i,j,nj
    integer :: found,found_recv,mainpe,inblock
    integer :: status(mpi_status_size)
    integer :: pe_grid_nb(2*(3**ndim),2),pe_grid(2)
    integer :: numblock,flag

    double precision :: igrid_nb,ipe_nb

    mainpe=0
    found=0
    numblock=2*(3**ndim)
    ipe_next=-1
    igrid_next=-1
    pe_grid(1)=-1
    pe_grid(2)=-1

    ! record neighbor blocks
    if (mype==ipe_now) then
      j=1
      {do ix^D=-1,1,1\}
        pe_grid_nb(j,2)=neighbor(1,ix^D,igrid_now)
        pe_grid_nb(j+numblock/2,2)=neighbor_child(1,ix^D,igrid_now)
        pe_grid_nb(j,1)=neighbor(2,ix^D,igrid_now)
        pe_grid_nb(j+numblock/2,1)=neighbor_child(2,ix^D,igrid_now)
        j=j+1
      {enddo\}
    endif


    call MPI_BCAST(pe_grid_nb,2*numblock,MPI_INTEGER,ipe_now,icomm,ierrmpi)

    ! check neighbors
    LOOPNB: do j=1,numblock
      if (mype==pe_grid_nb(j,1) .and. igrid_inuse(pe_grid_nb(j,2),pe_grid_nb(j,1)) .and. &
          pe_grid_nb(j,2)>=0 .and. pe_grid_nb(j,2)<=max_blocks) then
        ^D&xbmin^D=rnode(rpxmin^D_,pe_grid_nb(j,2))\
        ^D&xbmax^D=rnode(rpxmax^D_,pe_grid_nb(j,2))\
        inblock=0
        {if (xf1(^D)>=xbmin^D .and. xf1(^D)<xbmax^D) inblock=inblock+1\}
        if (inblock==ndim) then
          pe_grid(2)=pe_grid_nb(j,2)
          pe_grid(1)=mype
          found=1
          call MPI_SEND(pe_grid,2,MPI_INTEGER,ipe_now,nj+10,icomm,ierrmpi)
          exit LOOPNB
        endif
      endif
    enddo LOOPNB

    call MPI_ALLREDUCE(found,found_recv,1,MPI_INTEGER,MPI_SUM,icomm,ierrmpi)
    found=found_recv

    ! point is not in neighbors
    if (found==0) then
      LOOP1: do iigrid=1,igridstail; igrid=igrids(iigrid);
        !call update_block_para(igrid,dxb^D,xb^L)
        ^D&xbmin^D=rnode(rpxmin^D_,igrid)\
        ^D&xbmax^D=rnode(rpxmax^D_,igrid)\
        inblock=0
        {if (xf1(^D)>=xbmin^D .and. xf1(^D)<xbmax^D) inblock=inblock+1\}
        if (inblock==ndim) then
          pe_grid(2)=igrid
          pe_grid(1)=mype
          call MPI_SEND(pe_grid,2,MPI_INTEGER,ipe_now,nj+10,icomm,ierrmpi)
          exit LOOP1
        endif
      enddo LOOP1
    endif


    if (mype==ipe_now) then
      call MPI_RECV(pe_grid,2,MPI_INTEGER,MPI_ANY_SOURCE,nj+10,icomm,status,ierrmpi)
    endif
    call MPI_BCAST(pe_grid,2,MPI_INTEGER,ipe_now,icomm,ierrmpi)
    ipe_next=pe_grid(1)
    igrid_next=pe_grid(2)

  end subroutine find_grid
  
end module mod_trace_Bfield
