module mod_point_searching
  use mod_global_parameters
  use mod_physics

  implicit none

contains

  subroutine get_point_w(xp,wp)
    ! for given point (xp), provide the plasma parameters (wp) at this point

    double precision :: xp(1:ndim),wp(1:nw)

    double precision :: dxb^D,xb^L
    integer :: inblock,indomain
    integer :: iigrid,igrid,igrid_now,igrid_next,j
    integer :: mainpe
    integer :: status(mpi_status_size)
    integer :: ixO^L,ixO^D,ixbl^D,ix^D
    double precision :: xd^D,factor(0:1^D&)

    mainpe=0

    ! check whether or the first point is inside simulation box. if yes, find
    ! the pe and igrid for the point
    indomain=0
    {if (xp(^DB)>=xprobmin^DB .and. xp(^DB)<xprobmax^DB) indomain=indomain+1\}
    if (indomain==ndim) then

      ! find the grid and pe of the first point
      loop1: do iigrid=1,igridstail; igrid=igrids(iigrid);
        ^D&xbmin^D=rnode(rpxmin^D_,igrid)\
        ^D&xbmax^D=rnode(rpxmax^D_,igrid)\
        inblock=0
        {if (xp(^DB)>=xbmin^DB .and. xp(^DB)<xbmax^DB) inblock=inblock+1\}
        if (inblock==ndim) then
          ^D&dxb^D=rnode(rpdx^D_,igrid)\
          ^D&ixOmin^D=ixmlo^D\
          ^D&ixbl^D=floor((xp(^D)-ps(igrid)%x(ixOmin^DD,^D))/dxb^D)+ixOmin^D\
          ^D&xd^D=(xp(^D)-ps(igrid)%x(ixbl^DD,^D))/dxb^D\
          {do ix^D=0,1\}
            factor(ix^D)={abs(1-ix^D-xd^D)*}
          {enddo\}

          ! do interpolation to get point values
          wp=0.d0
          {do ix^D=0,1\}
            do j=1,nw
              wp(j)=wp(j)+factor(ix^D)*ps(igrid)%w(ixbl^D+ix^D,j)
            enddo
          {enddo\}

          ! do interpolation to get point magnetic field
          if(physics_type=='mhd') then
            wp(iw_mag(1):iw_mag(ndir))=0.d0
            do j=1,ndir
              if (b0field) then
                wp(iw_mag(j))=ps(igrid)%w(ixbl^D+ix^D,iw_mag(j))+&
                                      ps(igrid)%B0(ixbl^D+ix^D,j,0)
              else
                wp(iw_mag(j))=ps(igrid)%w(ixbl^D+ix^D,iw_mag(j))
              endif
            enddo
          endif

          call MPI_SEND(wp,nw,MPI_DOUBLE_PRECISION,mainpe,0,icomm,ierrmpi)
          exit loop1
        endif
      enddo loop1

      if (mype==mainpe) then
        call MPI_RECV(wp,nw,MPI_DOUBLE_PRECISION,mpi_any_source,0,icomm,status,ierrmpi)
      endif
      call MPI_BCAST(wp,nw,MPI_DOUBLE_PRECISION,mainpe,icomm,ierrmpi)
    endif

  end subroutine get_point_w

  subroutine get_cell_w(xp,wc)
    ! for given point (xp), looking for corresponding cell and then provide
    ! the plasma parameters (wc) of this cell

    double precision :: xp(1:ndim),wc(1:nw)

    double precision :: dxb^D,xb^L
    integer :: inblock,indomain
    integer :: iigrid,igrid,igrid_now,igrid_next,j
    integer :: mainpe
    integer :: status(mpi_status_size)
    integer :: ixO^L,ixO^D,ixb^D

    mainpe=0

    ! check whether or the first point is inside simulation box. if yes, find
    ! the pe and igrid for the point
    indomain=0
    {if (xp(^DB)>=xprobmin^DB .and. xp(^DB)<xprobmax^DB) indomain=indomain+1\}
    if (indomain==ndim) then

      ! find the grid and pe of the first point
      loop1: do iigrid=1,igridstail; igrid=igrids(iigrid);
        ^D&xbmin^D=rnode(rpxmin^D_,igrid)\
        ^D&xbmax^D=rnode(rpxmax^D_,igrid)\
        inblock=0
        {if (xp(^DB)>=xbmin^DB .and. xp(^DB)<xbmax^DB) inblock=inblock+1\}
        if (inblock==ndim) then
          ^D&dxb^D=rnode(rpdx^D_,igrid)\
          ^D&ixOmin^D=ixmlo^D\

          ! looking for the cell index
          ^D&ixb^D=floor((xp(^D)-xbmin^D)/dxb^D)+ixomin^D\

          wc=0.d0
          do j=1,nw
            wc(j)=ps(igrid)%w(ixb^D,j)
          enddo

          ! for magnetic field
          if(physics_type=='mhd') then
            wc(iw_mag(1):iw_mag(ndir))=0.d0
            do j=1,ndir
              if (b0field) then
                wc(iw_mag(j))=ps(igrid)%w(ixb^D,iw_mag(j))+&
                              ps(igrid)%B0(ixb^D,j,0)
              else
                wc(iw_mag(j))=ps(igrid)%w(ixb^D,iw_mag(j))
              endif
            enddo
          endif

          call MPI_SEND(wc,nw,MPI_DOUBLE_PRECISION,mainpe,0,icomm,ierrmpi)
          exit loop1
        endif
      enddo loop1

      if (mype==mainpe) then
        call MPI_RECV(wc,nw,MPI_DOUBLE_PRECISION,mpi_any_source,0,icomm,status,ierrmpi)
      endif
      call MPI_BCAST(wc,nw,MPI_DOUBLE_PRECISION,mainpe,icomm,ierrmpi)
    endif

  end subroutine get_cell_w

  subroutine get_cell_index(xp,ipe,iblock,ixc^D)
    ! for given point (xp), provide the corrsponding processor (ipe), 
    ! grid number (iblock) and cell index (ixc^D)

    double precision :: xp(ndim)
    integer :: ipe,iblock,ixc^D

    double precision :: dxb^D,xb^L
    integer :: inblock,indomain
    integer :: iigrid,igrid,igrid_now,igrid_next,j
    integer :: mainpe,datas(ndim+2)
    integer :: status(mpi_status_size)
    integer :: ixO^L,ixO^D,ixb(ndim)

    mainpe=0
    datas=0

    ! check whether or the first point is inside simulation box. if yes, find
    ! the pe and igrid for the point
    indomain=0
    {if (xp(^DB)>=xprobmin^DB .and. xp(^DB)<xprobmax^DB) indomain=indomain+1\}
    if (indomain==ndim) then

      ! find the grid and pe of the first point
      loop1: do iigrid=1,igridstail; igrid=igrids(iigrid);
        ^D&xbmin^D=rnode(rpxmin^D_,igrid)\
        ^D&xbmax^D=rnode(rpxmax^D_,igrid)\
        inblock=0
        {if (xp(^DB)>=xbmin^DB .and. xp(^DB)<xbmax^DB) inblock=inblock+1\}
        if (inblock==ndim) then
          ^D&dxb^D=rnode(rpdx^D_,igrid)\
          ^D&ixOmin^D=ixmlo^D\

          ! looking for the cell index
          ^D&ixb(^D)=floor((xp(^D)-xbmin^D)/dxb^D)+ixomin^D\

          datas(1)=mype
          datas(2)=igrid
          do j=1,ndim
            datas(2+j)=ixb(j)
          enddo

          call MPI_SEND(datas,2+ndim,MPI_INTEGER,mainpe,0,icomm,ierrmpi)
          exit loop1
        endif
      enddo loop1

      if (mype==mainpe) then
        call MPI_RECV(datas,2+ndim,MPI_INTEGER,mpi_any_source,0,icomm,status,ierrmpi)
      endif
      call MPI_BCAST(datas,2+ndim,MPI_INTEGER,mainpe,icomm,ierrmpi)
    endif

    ipe=datas(1)
    iblock=datas(2)
    ^D&ixc^D=datas(2+^D)\

  end subroutine get_cell_index

end module mod_point_searching

