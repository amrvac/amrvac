!> Module for handling split source terms (split from the fluxes)
module mod_source

  implicit none
  private

  public :: addsource_all
  public :: addsource2

contains

  subroutine addsource_all(prior)
    use mod_global_parameters
    use mod_ghostcells_update

    logical, intent(in) :: prior

    double precision :: qdt, qt
    integer :: iigrid, igrid, i^D

    if ((.not.prior).and.&
         (typesourcesplit=='sf' .or. typesourcesplit=='ssf')) return

    if (prior) then
       qt=global_time
    else
       qt=global_time+dt
    end if
    !$OMP PARALLEL DO PRIVATE(igrid,qdt,i^D)
    do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
       qdt=dt_grid(igrid)
       block=>pw(igrid)
       {do i^DB=-1,1\}
       if (i^D==0|.and.) cycle
       if (neighbor_type(i^D,igrid)==2 .or. neighbor_type(i^D,igrid)==4) then
          leveljump(i^D)=.true.
       else
          leveljump(i^D)=.false.
       end if
       {end do\}
       call addsource1_grid(igrid,qdt,qt,pw(igrid)%w)
    end do
    !$OMP END PARALLEL DO

    call getbc(qt,0.d0,0,nwflux+nwaux)

  end subroutine addsource_all

  subroutine addsource1_grid(igrid,qdt,qt,w)

    use mod_global_parameters

    integer, intent(in) :: igrid
    double precision, intent(in) :: qdt, qt
    double precision, intent(inout) :: w(ixG^T,nw)

    double precision :: w1(ixG^T,nw)

    saveigrid=igrid
    typelimiter=type_limiter(node(plevel_,igrid))
    typegradlimiter=type_gradient_limiter(node(plevel_,igrid))

    ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);

    w1(ixG^T,1:nw)=w(ixG^T,1:nw)

    select case (typesourcesplit)
    case ('sf')
       call addsource2(qdt  ,ixG^LL,ixM^LL,1,nw,qt,w1,qt,w,pw(igrid)%x,.true.)
    case ('sfs')
       call addsource2(qdt/2,ixG^LL,ixM^LL,1,nw,qt,w1,qt,w,pw(igrid)%x,.true.)
    case ('ssf')
       call addsource2(qdt/2,ixG^LL,ixG^LL,1,nw,qt,w,qt,w1,pw(igrid)%x,.true.)
       call addsource2(qdt  ,ixG^LL,ixM^LL,1,nw,qt,w1,qt,w,pw(igrid)%x,.true.)
    case ('ssfss')
       call addsource2(qdt/4,ixG^LL,ixG^LL,1,nw,qt,w,qt,w1,pw(igrid)%x,.true.)
       call addsource2(qdt/2,ixG^LL,ixM^LL,1,nw,qt,w1,qt,w,pw(igrid)%x,.true.)
    case default
       write(unitterm,*)'No such typesourcesplit=',typesourcesplit
       call mpistop("Error: Unknown typesourcesplit!")
    end select

  end subroutine addsource1_grid

  subroutine addsource2(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,qsourcesplit)

    ! Add source within ixO for iws: w=w+qdt*S[wCT]

    use mod_global_parameters
    use mod_physics, only: phys_add_source
    use mod_usr_methods, only: usr_source
    ! differences with VAC is in iw^LIM and in declaration of ranges for wCT,w

    integer, intent(in) :: ixI^L, ixO^L, iw^LIM
    double precision, intent(in) :: qdt, qtC, qt
    double precision, intent(in) :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    logical, intent(in) :: qsourcesplit

    ! user defined sources, typically explicitly added
    if ((qsourcesplit .eqv. source_split_usr) .and. associated(usr_source)) then
       call usr_source(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)
    end if

    ! physics defined sources, typically explicitly added,
    ! along with geometrical source additions
    call phys_add_source(qdt,ixI^L,ixO^L,wCT,w,x,qsourcesplit)

  end subroutine addsource2

end module mod_source
