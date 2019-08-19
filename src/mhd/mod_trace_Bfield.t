!> mod_trace_Bfield -- trace magnetic field line
!>
!> This module should be used in a global subroutine, since it will 
!> search blocks/grids controlled by different processors.
!>
!> We can use subroutine trace_Bfield to trace a magnetic field line.
!> The location of a point in this field line need to be provided at
!> xf(1,:), and then the subroutine will calculate and return the 
!> locations of some other points in this field line. The locations and
!> the plasma parameter of the points will be return. Parameters are 
!> recorded in wB (rho, mom, e, B, current)
!> 
!> numP is the number of points we wants to return. Sometimes field
!> lines can go out of the simulation box. It means that only some of 
!> the points are valid. numRT records that how many points are valid.
!>
!> The disctance between two adjacent points is determined by dL and dirct.
!> For example, if dirct=1, the distance in x direction (dx) is a constant. 
!> This works when 1<=dirct<=ndim. For the conditions that dirct<1 or
!> dirct>ndim, dL=sqrt(dx^2+dy^2+dz^2).
!>
!> Parameter forward determines the direction we trace a field line.
!> If forward is true, the tracing direction is parallel to the B field
!> direction. If forward is false, the tracing direction is antiparallel
!> to the B field direction.

module mod_trace_Bfield
  use mod_mhd

  contains

    subroutine trace_Bfield(xf,wB,dL,dirct,numP,numRT,forward)
      ! trace a field line
      use mod_usr_methods
      use mod_global_parameters

      integer :: numP,numRT,dirct
      double precision :: xf(numP,ndim),wB(numP,nw+ndir)
      double precision :: dL
      logical :: forward

      double precision :: dxb^D,xb^L
      integer :: inblock,indomain
      integer :: iigrid,igrid,j
      integer :: igrid_now,igrid_next,igrid_rc
      integer :: ipe_now,ipe_next,ipe_rc
      integer :: mainpe
      integer :: status(mpi_status_size)
      logical :: newpe
      double precision :: xtemp(ndim)

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
            igrid_now=igrid
            ipe_now=mype
            call MPI_SEND(igrid_now,1,MPI_INTEGER,mainpe,0,icomm,ierrmpi)
            call MPI_SEND(ipe_now,1,MPI_INTEGER,mainpe,1,icomm,ierrmpi)
            exit LOOP1
          endif
        enddo LOOP1

        if (mype==mainpe) then
          call MPI_RECV(igrid_rc,1,MPI_INTEGER,MPI_ANY_SOURCE,0,icomm,status,ierrmpi)
          call MPI_RECV(ipe_rc,1,MPI_INTEGER,MPI_ANY_SOURCE,1,icomm,status,ierrmpi)
          igrid_now=igrid_rc
          ipe_now=ipe_rc
        endif
        call MPI_BCAST(igrid_now,1,MPI_INTEGER,mainpe,icomm,ierrmpi)
        call MPI_BCAST(ipe_now,1,MPI_INTEGER,mainpe,icomm,ierrmpi)
      endif


      ! if the first point is not in simulation box
      ! other points in field line    
      if (indomain==ndim) then
        LOOPMAIN: do j=1,numP-1
          newpe=.TRUE.
          mainpe=ipe_now
          if (mype==ipe_now) then
            igrid=igrid_now
            ! find next point in this field line and get density of current point
            call find_next(igrid,xf(j,:),xf(j+1,:),wB(j,:),dL,dirct,forward)
            ! check the pe of next point
            call check_next(igrid,igrid_next,ipe_next,xf(j+1,:),newpe)
          endif

          call MPI_BCAST(newpe,1,MPI_LOGICAL,mainpe,icomm,ierrmpi)
          call MPI_BCAST(xf(j+1,:),ndim,MPI_DOUBLE_PRECISION,mainpe,icomm,ierrmpi)
          call MPI_BCAST(wB(j,:),nw+ndir,MPI_DOUBLE_PRECISION,mainpe,icomm,ierrmpi)

          ! check whether or next point is inside simulation box.
          indomain=0
          {if (xf(j+1,^DB)>xprobmin^DB .and. xf(j+1,^DB)<xprobmax^DB) indomain=indomain+1\}
          if (indomain==ndim) then
            numRT=j+1
          else
            exit LOOPMAIN
          endif

          ! next point is in another pe, find out grid and pe numbers
          if (newpe .eqv. .TRUE.) then
            call find_grid(ipe_now,igrid_now,ipe_next,igrid_next,xf(j+1,:),j)
          endif

          ! prepare for next point 
          if (newpe .eqv. .TRUE.) then
            ipe_now=ipe_next
            igrid_now=igrid_next
          else
            if (mype==ipe_now) then
              igrid_now=igrid_next
            endif
          endif

        enddo LOOPMAIN
      endif


      ! density and magnetic field for the last point
      if (numRT==numP) then
        mainpe=ipe_now
        if (mype==ipe_now) then
          call find_next(igrid_now,xf(numRT,:),xtemp,wB(numRT,:), &
                         dL,dirct,forward)
        endif
        call MPI_BCAST(wB(numRT,:),nw+ndir,MPI_DOUBLE_PRECISION,mainpe,icomm,ierrmpi)
      endif

    end subroutine trace_Bfield

    subroutine find_next(igrid,xf0,xf1,wB0,dL,dirct,forward)
      !find next point
      use mod_usr_methods
      use mod_global_parameters
      use mod_forest

      integer :: igrid,dirct
      double precision :: xf0(ndim),xf1(ndim),wB0(nw+ndir)
      double precision :: dL
      logical :: forward

      integer          :: ixO^L,ixO^D,j
      double precision :: dxf(ndim)
      double precision :: dxb^D,xb^L,xd^D
      integer          :: ixb^D,ix^D,ixbl^D
      double precision :: wBnear(0:1^D&,nw+ndir)
      double precision :: Bx(ndim),factor(0:1^D&)
      double precision :: Btotal,maxft,Bp

      double precision :: current(ixg^T,7-2*ndir:3)
      integer :: idirmin,ixI^L,idir,jdir,kdir
      double precision :: dB(ndir,ndim),dxb(ndim)
      integer :: hxO^L,jxO^L,nxO^L
      double precision :: tmp(0:1^D&)


      ^D&ixImin^D=ixglo^D\
      ^D&ixImax^D=ixghi^D\
      ^D&ixOmin^D=ixmlo^D\
      ^D&ixOmax^D=ixmhi^D\

      ^D&dxb^D=rnode(rpdx^D_,igrid)\
      ^D&xbmin^D=rnode(rpxmin^D_,igrid)\
      ^D&xbmax^D=rnode(rpxmax^D_,igrid)\

      ^D&ixbl^D=floor((xf0(^D)-ps(igrid)%x(ixOmin^DD,^D))/dxb^D)+ixOmin^D\
      ^D&xd^D=(xf0(^D)-ps(igrid)%x(ixbl^DD,^D))/dxb^D\

      ! get current
      !call get_current(ps(igrid)%w,ixI^L,ixO^L,idirmin,current)
      current=0
      idirmin=7-2*ndir
      ^D&dxb(^D)=dxb^D\
      wBnear=0


      ! nearby B field for interpolation
      if (B0field) then
        {do ix^DB=0,1\}
          do j=1,ndir
            wBnear(ix^D,mag(j))=ps(igrid)%w(ixbl^D+ix^D,mag(j))+&
                          ps(igrid)%B0(ixbl^D+ix^D,j,0)
          enddo
        {enddo\}
      else
        {do ix^DB=0,1\}
          do j=1,ndir
            wBnear(ix^D,mag(j))=ps(igrid)%w(ixbl^D+ix^D,mag(j))
          enddo
        {enddo\}
      endif

      ! nearby current for interpolation
      ^D&nxOmin^D=0\
      ^D&nxOmax^D=1\
      do idir=idirmin,3; do jdir=1,ndim; do kdir=1,ndir
        if (lvc(idir,jdir,kdir)/=0) then
          hxO^L=ixbl^D+nxO^L-kr(jdir,^D);
          jxO^L=ixbl^D+nxO^L+kr(jdir,^D);
          if (B0field) then
            tmp(nxO^S)=half*(ps(igrid)%w(jxO^S,mag(kdir))+&
                             ps(igrid)%B0(jxO^S,kdir,0)-&
                             ps(igrid)%w(hxO^S,mag(kdir))-&
                             ps(igrid)%B0(hxO^S,kdir,0))/dxb(jdir)
          else
            tmp(nxO^S)=half*(ps(igrid)%w(jxO^S,mag(kdir))-&
                            ps(igrid)%w(hxO^S,mag(kdir)))/dxb(jdir)
          endif  

          if (lvc(idir,jdir,kdir)==1) then
            wBnear(nxO^S,nw+idir)=wBnear(nxO^S,nw+idir)+tmp(nxO^S)
          else
            wBnear(nxO^S,nw+idir)=wBnear(nxO^S,nw+idir)-tmp(nxO^S)
          endif
        endif
      enddo; enddo; enddo


      ! other parameter for interpolation
      {do ix^DB=0,1\}
        do j=1,nw
          if (j<mag(1) .or. j>mag(ndir)) then
            wBnear(ix^D,j)=ps(igrid)%w(ixbl^D+ix^D,j)
          endif
        enddo
      {enddo\}


      ! interpolation factor
      {do ix^D=0,1\}
        factor(ix^D)={abs(1-ix^D-xd^D)*}
      {enddo\}

      ! do interpolation to get local magnetic field
      Bx=0
      wB0=0
      {do ix^DB=0,1\}
        {Bx(^DB)=Bx(^DB)+wBnear(ix^DD,mag(^DB))*factor(ix^DD)\}
        do j=1,nw+ndir
          wB0(j)=wB0(j)+wBnear(ix^D,j)*factor(ix^D)
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
          if (dirct>=1 .and. dirct<=ndim) then
            dxf(j)=dL*Bx(j)/abs(Bx(dirct))
          else
            dxf(j)=dL*Bx(j)/Btotal
          endif
        enddo
      else
        do j=1,ndim
          if (dirct>=1 .and. dirct<=ndim) then
            dxf(j)=-dL*Bx(j)/abs(Bx(dirct))
          else
            dxf(j)=-dL*Bx(j)/Btotal
          endif
        enddo
      endif

      do j=1,ndim
        xf1(j)=xf0(j)+dxf(j)
      enddo

    end subroutine find_next

    subroutine check_next(igrid,igrid_next,ipe_next,xf1,newpe)
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

    end subroutine check_next

    subroutine find_grid(ipe_now,igrid_now,ipe_next,igrid_next,xf1,nj)
      ! find for grid and pe numbers
      use mod_usr_methods
      use mod_global_parameters
      use mod_forest

      integer :: igrid,iigrid,igrid_now,ipe_now,igrid_next,ipe_next
      integer :: igrid_recv,ipe_recv
      double precision :: xf1(ndim)

      double precision :: dxb^D,xb^L
      integer :: inblock
      integer :: ix^D,i,j,nj
      integer :: found,found_recv,mainpe
      integer :: status(mpi_status_size)
      integer :: grid_nb(2*(3**ndim)),pe_nb(2*(3**ndim))
      integer :: numblock,flag

      double precision :: igrid_nb,ipe_nb

      mainpe=0
      found=0
      numblock=2*(3**ndim)

      ! record neighbor blocks
      if (mype==ipe_now) then
        j=1
        {do ix^D=-1,1,1\}
          grid_nb(j)=neighbor(1,ix^D,igrid_now)
          pe_nb(j)=neighbor(2,ix^D,igrid_now)
          grid_nb(j+numblock/2)=neighbor_child(1,ix^D,igrid_now)
          pe_nb(j+numblock/2)=neighbor_child(2,ix^D,igrid_now)
          j=j+1
        {enddo\}
      endif

      call MPI_BCAST(grid_nb,numblock,MPI_INTEGER,ipe_now,icomm,ierrmpi)
      call MPI_BCAST(pe_nb,numblock,MPI_INTEGER,ipe_now,icomm,ierrmpi)

      ! check neighbors
      LOOPNB: do j=1,numblock
        if (mype==pe_nb(j) .and. igrid_inuse(grid_nb(j),pe_nb(j)) .and. &
            grid_nb(j)>=0 .and. grid_nb(j)<=max_blocks) then
          ^D&xbmin^D=rnode(rpxmin^D_,grid_nb(j))\
          ^D&xbmax^D=rnode(rpxmax^D_,grid_nb(j))\
          inblock=0
          {if (xf1(^D)>=xbmin^D .and. xf1(^D)<xbmax^D) inblock=inblock+1\}
          if (inblock==ndim) then
            igrid_next=grid_nb(j)
            ipe_next=mype
            found=1
            call MPI_SEND(igrid_next,1,MPI_INTEGER,ipe_now,nj+10,icomm,ierrmpi)
            call MPI_SEND(ipe_next,1,MPI_INTEGER,ipe_now,nj+11,icomm,ierrmpi)
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
            igrid_next=igrid
            ipe_next=mype
            call MPI_SEND(igrid_next,1,MPI_INTEGER,ipe_now,nj+10,icomm,ierrmpi)
            call MPI_SEND(ipe_next,1,MPI_INTEGER,ipe_now,nj+11,icomm,ierrmpi)
            exit LOOP1
          endif
        enddo LOOP1
      endif


      if (mype==ipe_now) then
        call MPI_RECV(igrid_recv,1,MPI_INTEGER,MPI_ANY_SOURCE,nj+10,icomm,status,ierrmpi)
        call MPI_RECV(ipe_recv,1,MPI_INTEGER,MPI_ANY_SOURCE,nj+11,icomm,status,ierrmpi)
        igrid_next=igrid_recv
        ipe_next=ipe_recv
      endif
      call MPI_BCAST(igrid_next,1,MPI_INTEGER,ipe_now,icomm,ierrmpi)
      call MPI_BCAST(ipe_next,1,MPI_INTEGER,ipe_now,icomm,ierrmpi)

    end subroutine find_grid

end module mod_trace_Bfield
