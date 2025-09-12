!> Module for handling split source terms (split from the fluxes)
module mod_source
  implicit none
  public

  !> How to apply dimensional splitting to the source terms, see
  !> @ref discretization.md
  integer :: sourcesplit =-1
  integer, parameter :: sourcesplit_sfs    = 0
  integer, parameter :: sourcesplit_sf     = 1
  integer, parameter :: sourcesplit_ssf    = 2
  integer, parameter :: sourcesplit_ssfss  = 3

  public :: add_split_source
  public :: addsource2

contains

  subroutine add_split_source(prior)
    use mod_global_parameters
    use mod_ghostcells_update
    use mod_physics, only: phys_global_source_after, phys_to_primitive
    use mod_supertimestepping, only: is_sts_initialized, sts_add_source,sourcetype_sts,&
                                      sourcetype_sts_prior, sourcetype_sts_after, sourcetype_sts_split   
    use mod_comm_lib, only: mpistop

    logical, intent(in) :: prior
    ! This variable, later allocated on the thread stack, causes segmentation fault
    ! when openmp is used with intel. That could be solved otherwise, by increasing
    ! the thread stack size, but not using it at all could speed up. 
    double precision :: w1(ixG^T,1:nw)
    double precision :: qt
    integer :: iigrid, igrid
    logical :: src_active

    ! add stiff source terms via super time stepping
    if(is_sts_initialized()) then
        select case (sourcetype_sts)
          case (sourcetype_sts_prior)
            if(prior) then
              call sts_add_source(dt)
            end if
          case (sourcetype_sts_after)
            if(.not. prior) then
              call sts_add_source(dt)
            end if
          case (sourcetype_sts_split)
            call sts_add_source(0.5d0*dt)
          end select
    end if
    src_active = .false.

    if ((.not.prior).and.&
         (sourcesplit==sourcesplit_sf .or. sourcesplit==sourcesplit_ssf)) return

    if (prior) then
       qt=global_time
    else
       qt=global_time+dt
    end if

    if(any_source_split) then
      ! add normal split source terms
      select case (sourcesplit)
      case (sourcesplit_sfs)
        !$OMP PARALLEL DO PRIVATE(igrid,w1)
        do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
           block=>ps(igrid)
           ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
           w1=ps(igrid)%w
           call phys_to_primitive(ixG^LL,ixG^LL,w1,ps(igrid)%x)
           call addsource2(0.5d0*dt,0.5d0,ixG^LL,ixM^LL,1,nw,qt,ps(igrid)%w,w1,qt,ps(igrid)%w,&
                ps(igrid)%x,.true.,src_active)
        end do
        !$OMP END PARALLEL DO
      case (sourcesplit_sf)
        !$OMP PARALLEL DO PRIVATE(igrid,w1)
        do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
           block=>ps(igrid)
           ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
           w1=ps(igrid)%w
           call phys_to_primitive(ixG^LL,ixG^LL,w1,ps(igrid)%x)
           call addsource2(dt  ,1d0,ixG^LL,ixM^LL,1,nw,qt,ps(igrid)%w,w1,qt,ps(igrid)%w,&
                ps(igrid)%x,.true.,src_active)
        end do
        !$OMP END PARALLEL DO
      case (sourcesplit_ssf)
        !$OMP PARALLEL DO PRIVATE(igrid,w1)
        do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
           block=>ps(igrid)
           ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
           w1=ps(igrid)%w
           call phys_to_primitive(ixG^LL,ixG^LL,w1,ps(igrid)%x)
           call addsource2(0.5d0*dt,0.5d0,ixG^LL,ixG^LL,1,nw,qt,ps(igrid)%w,w1,qt,ps(igrid)%w,&
                ps(igrid)%x,.true.,src_active)
           call addsource2(dt  ,1d0,ixG^LL,ixM^LL,1,nw,qt,ps(igrid)%w,w1,qt,ps(igrid)%w,&
                ps(igrid)%x,.true.,src_active)
        end do
        !$OMP END PARALLEL DO
      case (sourcesplit_ssfss)
        !$OMP PARALLEL DO PRIVATE(igrid,w1)
        do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
           block=>ps(igrid)
           ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
           w1=ps(igrid)%w
           call phys_to_primitive(ixG^LL,ixG^LL,w1,ps(igrid)%x)
           call addsource2(0.25d0*dt,0.25d0,ixG^LL,ixG^LL,1,nw,qt,ps(igrid)%w,w1,qt,ps(igrid)%w,&
                ps(igrid)%x,.true.,src_active)
           call addsource2(0.5d0*dt,0.5d0,ixG^LL,ixM^LL,1,nw,qt,ps(igrid)%w,w1,qt,ps(igrid)%w,&
                ps(igrid)%x,.true.,src_active)
        end do
        !$OMP END PARALLEL DO
      case default
         write(unitterm,*)'No such type of sourcesplit=',sourcesplit
         call mpistop("Error: Unknown type of sourcesplit!")
      end select
    end if

    if (.not. prior .and. associated(phys_global_source_after)) then
       call phys_global_source_after(dt, qt, src_active)
    end if

    if (src_active) then
       call getbc(qt,0.d0,ps,iwstart,nwgc)
    end if

  end subroutine add_split_source

  !> Add source within ixO for iws: w=w+qdt*S[wCT]
  subroutine addsource2(qdt,dtfactor,ixI^L,ixO^L,iw^LIM,qtC,wCT,wCTprim,qt,&
       w,x,qsourcesplit,src_active)
    use mod_global_parameters
    use mod_physics, only: phys_add_source
    use mod_usr_methods, only: usr_source
    ! differences with VAC is in iw^LIM and in declaration of ranges for wCT,w

    integer, intent(in)              :: ixI^L, ixO^L, iw^LIM
    double precision, intent(in)     :: qdt, dtfactor, qtC, qt
    double precision, intent(in)     :: wCT(ixI^S,1:nw), wCTprim(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout)  :: w(ixI^S,1:nw)
    logical, intent(in)              :: qsourcesplit
    logical, intent(inout), optional :: src_active

    logical                          :: tmp_active

    tmp_active = .false.

    ! user defined sources, typically explicitly added
    if ((qsourcesplit .eqv. source_split_usr) .and. associated(usr_source)) then
       tmp_active = .true.
       call usr_source(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)
    end if

    ! physics defined sources, typically explicitly added,
    ! along with geometrical source additions
    call phys_add_source(qdt,dtfactor,ixI^L,ixO^L,wCT,wCTprim,w,x,qsourcesplit,tmp_active)

    if (present(src_active)) src_active = src_active .or. tmp_active

  end subroutine addsource2

end module mod_source
