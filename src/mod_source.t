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
    !-----------------------------------------------------------------------------

    if ((.not.prior).and.&
         (typesourcesplit=='sf' .or. typesourcesplit=='ssf')) return

    if (prior) then
       qt=t
    else
       qt=t+dt
    end if
    !$OMP PARALLEL DO PRIVATE(igrid,qdt,i^D)
    do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
       qdt=dt_grid(igrid)
       if (B0field) then
          myB0_cell => pB0_cell(igrid)
          {^D&myB0_face^D => pB0_face^D(igrid)\}
       end if
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

    call getbc(qt,0.d0,pw,0,nwflux+nwaux)

  end subroutine addsource_all

  subroutine addsource1_grid(igrid,qdt,qt,w)

    use mod_global_parameters

    integer, intent(in) :: igrid
    double precision, intent(in) :: qdt, qt
    double precision, intent(inout) :: w(ixG^T,nw)

    double precision :: w1(ixG^T,nw)
    !-----------------------------------------------------------------------------

    if (.not.slab) mygeo => pgeo(igrid)

    saveigrid=igrid
    typelimiter=typelimiter1(node(plevel_,igrid))
    typegradlimiter=typegradlimiter1(node(plevel_,igrid))

    ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);

    w1(ixG^T,1:nw)=w(ixG^T,1:nw)

    select case (typesourcesplit)
    case ('sf')
       call addsource2(qdt  ,ixG^LL,ixM^LL,1,nw,qt,w1,qt,w,px(igrid)%x,.true.)
    case ('sfs')
       call addsource2(qdt/2,ixG^LL,ixM^LL,1,nw,qt,w1,qt,w,px(igrid)%x,.true.)
    case ('ssf')
       call addsource2(qdt/2,ixG^LL,ixG^LL,1,nw,qt,w,qt,w1,px(igrid)%x,.true.)
       call addsource2(qdt  ,ixG^LL,ixM^LL,1,nw,qt,w1,qt,w,px(igrid)%x,.true.)
    case ('ssfss')
       call addsource2(qdt/4,ixG^LL,ixG^LL,1,nw,qt,w,qt,w1,px(igrid)%x,.true.)
       call addsource2(qdt/2,ixG^LL,ixM^LL,1,nw,qt,w1,qt,w,px(igrid)%x,.true.)
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
    double precision, intent(in) :: qdt, qtC, qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: wCT(ixI^S,1:nw), w(ixI^S,1:nw)
    logical, intent(in) :: qsourcesplit
    !-----------------------------------------------------------------------------

    ! user defined sources, typically explicitly added
    if ((qsourcesplit .eqv. ssplituser) .and. associated(usr_source)) then
       call usr_source(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)
    end if

    ! physics defined sources, typically explicitly added,
    ! along with geometrical source additions
    call phys_add_source(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,qsourcesplit)

  end subroutine addsource2

end module mod_source
