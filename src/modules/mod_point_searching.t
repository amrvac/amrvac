! for Cartesian coordinate system
module mod_point_searching
  use mod_global_parameters
  use mod_physics
  use mod_particle_base
  use mod_comm_lib, only: mpistop
  implicit none

contains

  subroutine get_point_w_ingrid(igrid,xp,wp,variable_type)
    double precision :: xp(1:ndim),wp(1:nw)
    character(*) :: variable_type
    integer, intent(in) :: igrid

    double precision :: dxb^D,xd^D,xb^L,temp
    double precision :: factor(0:1^D&)
    integer :: ixI^L,ixO^L
    integer :: ixbl^D,ix^D,ixA^L,j,ingrid

    ixI^L=ixG^LL;
    ixO^L=ixM^LL;
    ^D&xbmin^D=rnode(rpxmin^D_,igrid)-rnode(rpdx^D_,igrid);
    ^D&xbmax^D=rnode(rpxmax^D_,igrid)+rnode(rpdx^D_,igrid);
    ingrid=0
    {if (xp(^DB)>xbmin^DB .and. xp(^DB)<xbmax^DB) ingrid=ingrid+1\}

    if (ingrid==ndim) then
      ^D&dxb^D=rnode(rpdx^D_,igrid)\
      ^D&ixbl^D=floor((xp(^D)-ps(igrid)%x(ixOmin^DD,^D))/dxb^D)+ixOmin^D\
      ^D&xd^D=(xp(^D)-ps(igrid)%x(ixbl^DD,^D))/dxb^D\
      ^D&ixAmin^D=ixbl^D\
      ^D&ixAmax^D=ixbl^D+1\

      {do ix^D=0,1\}
        factor(ix^D)={abs(1-ix^D-xd^D)*}
      {enddo\}

      select case(variable_type)
        case('primitive')
          call phys_to_primitive(ixI^L,ixA^L,ps(igrid)%w,ps(igrid)%x)
        case('conserved')

        case default
          call mpistop("get_point_w_ingrid: Unknown variable type!")
      end select

      wp=0.d0
      {do ix^D=0,1\}
        do j=1,nw
          wp(j)=wp(j)+factor(ix^D)*ps(igrid)%w(ixbl^D+ix^D,j)
        enddo
      {enddo\}

      select case(variable_type)
        case('primitive')
        call phys_to_conserved(ixI^L,ixA^L,ps(igrid)%w,ps(igrid)%x)
        case('conserved')

        case default
          call mpistop("get_point_w_ingrid: Unknown variable type!")
      end select

      if(physics_type=='mhd') then
        wp(iw_mag(1):iw_mag(ndir))=0.d0
        {do ix^D=0,1\}
          do j=1,ndir
            if (b0field) then
              temp=ps(igrid)%w(ixbl^D+ix^D,iw_mag(j))+&
                                    ps(igrid)%B0(ixbl^D+ix^D,j,0)
              wp(iw_mag(j))=wp(iw_mag(j))+factor(ix^D)*temp
            else
              temp=ps(igrid)%w(ixbl^D+ix^D,iw_mag(j))
              wp(iw_mag(j))=wp(iw_mag(j))+factor(ix^D)*temp
            endif
          enddo
        {enddo\}
      endif
    else
      call mpistop("get_point_w_ingrid: The point is not in given grid!")
    endif

  end subroutine get_point_w_ingrid

  subroutine get_point_w(xp,wp,variable_type)
    ! for given point (xp), provide the plasma parameters (wp) at this point

    double precision :: xp(1:ndim),wp(1:nw)
    character(*) :: variable_type

    double precision :: x3d(3)
    integer :: indomain,ipe,igrid,j

    indomain=0
    {if (xp(^DB)>=xprobmin^DB .and. xp(^DB)<xprobmax^DB) indomain=indomain+1\}
    if (indomain==ndim) then

      ! find pe and igrid
      x3d=0.d0
      do j=1,ndim
        x3d(j)=xp(j)
      enddo
      call find_particle_ipe(x3d,igrid,ipe)

      !do interpolation to get point values
      if (mype==ipe) then
        call get_point_w_ingrid(igrid,xp,wp,variable_type)
      endif

      call MPI_BCAST(wp,nw,MPI_DOUBLE_PRECISION,ipe,icomm,ierrmpi)
    endif

  end subroutine get_point_w

  subroutine get_cell_w(xp,wc,variable_type)
    ! for given point (xp), looking for corresponding cell and then provide
    ! the plasma parameters (wc) of this cell

    double precision :: xp(1:ndim),wc(1:nw)
    character(*) :: variable_type

    double precision :: x3d(3)
    double precision :: dxb^D,xb^L
    integer :: indomain,ixO^L,ixb^D,ixA^L,ixI^L
    integer :: ipe,igrid,j

    indomain=0
    {if (xp(^DB)>=xprobmin^DB .and. xp(^DB)<xprobmax^DB) indomain=indomain+1\}
    if (indomain==ndim) then

      ! find pe and igrid
      x3d=0.d0
      do j=1,ndim
        x3d(j)=xp(j)
      enddo
      call find_particle_ipe(x3d,igrid,ipe)

      if (mype==ipe) then
        ! looking for the cell index
        ^D&xbmin^D=rnode(rpxmin^D_,igrid)\
        ^D&xbmax^D=rnode(rpxmax^D_,igrid)\
        ^D&dxb^D=rnode(rpdx^D_,igrid)\
        ^D&ixOmin^D=ixmlo^D\
        ^D&ixImin^D=ixglo^D\
        ^D&ixImax^D=ixghi^D\
        ^D&ixb^D=floor((xp(^D)-xbmin^D)/dxb^D)+ixOmin^D\
        ^D&ixAmin^D=ixb^D\
        ^D&ixAmax^D=ixb^D\

        wc=0.d0
        select case(variable_type)
          case('primitive')
            call phys_to_primitive(ixI^L,ixA^L,ps(igrid)%w,ps(igrid)%x)
            do j=1,nw
              wc(j)=ps(igrid)%w(ixb^D,j)
            enddo
            call phys_to_conserved(ixI^L,ixA^L,ps(igrid)%w,ps(igrid)%x)
          case('conserved')
            do j=1,nw
              wc(j)=ps(igrid)%w(ixb^D,j)
            enddo
          case default
            call mpistop("get_point_w: Unknown variable type!")
        end select

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
      endif

      call MPI_BCAST(wc,nw,MPI_DOUBLE_PRECISION,ipe,icomm,ierrmpi)
    endif

  end subroutine get_cell_w

  subroutine get_cell_index(xp,ipe,igrid,ixc^D)
    ! for given point (xp), provide the corrsponding processor (ipe), 
    ! grid number (igrid) and cell index (ixc^D)

    double precision :: xp(1:ndim)
    integer :: ipe,igrid,ixc^D

    double precision :: x3d(1:3)
    double precision :: dxb^D,xb^L
    integer :: indomain,ixO^L,j
    integer :: datas(1:ndim+2)

    indomain=0
    {if (xp(^DB)>=xprobmin^DB .and. xp(^DB)<xprobmax^DB) indomain=indomain+1\}

    if (indomain==ndim) then
      ! find pe and igrid
      x3d=0.d0
      do j=1,ndim
        x3d(j)=xp(j)
      enddo
      call find_particle_ipe(x3d,igrid,ipe)

      if (mype==ipe) then
        ! looking for the cell index
        ^D&xbmin^D=rnode(rpxmin^D_,igrid)\
        ^D&xbmax^D=rnode(rpxmax^D_,igrid)\
        ^D&dxb^D=rnode(rpdx^D_,igrid)\
        ^D&ixOmin^D=ixmlo^D\
        ^D&ixc^D=floor((xp(^D)-xbmin^D)/dxb^D)+ixOmin^D\
        ^D&datas(^D)=ixc^D\
        datas(ndim+1)=ipe
        datas(ndim+2)=igrid
      endif

      call MPI_BCAST(datas,ndim+2,MPI_INTEGER,ipe,icomm,ierrmpi)
    endif

    ^D&ixc^D=datas(^D)\
    ipe=datas(ndim+1)
    igrid=datas(ndim+2)

  end subroutine get_cell_index

end module mod_point_searching
