module mod_usr
  use mod_rho
  implicit none
  double precision, allocatable, save :: myBC(:,:,:,:) ! coord - coord - time - level
  integer :: time_old_bc, time_next_bc ! integer time allowed between snapshots (& dt cst)
  character(len=50) :: filename_bc_old, filename_bc_next

contains

  subroutine usr_init()
    usr_set_parameters => initglobaldata_usr
    usr_init_one_grid  => initonegrid_usr
    usr_special_bc     => specialbound_usr
    usr_process_global => process_global_usr
    usr_refine_grid    => specialrefine_grid

    call rho_activate()

  end subroutine usr_init

  character(len=20) function str(k)
       !   "Convert an integer to string."
       integer, intent(in) :: k
       write (str, *) k
       str = adjustl(str)
  end function str

  subroutine initglobaldata_usr
    integer :: i, j, ilev
    ! inflow from x minimal boundary
    allocate(myBC(domain_nx2*2**(refine_max_level-1),domain_nx3*2**(refine_max_level-1),2,refine_max_level))
    myBC(:,:,:,:)=0.0d0
    time_old_bc =0
    time_next_bc=1
    write (filename_bc_old,'(I2.2)') time_old_bc
    write (filename_bc_next,'(I2.2)') time_next_bc
    filename_bc_old ='bc_t'//trim(filename_bc_old)
    filename_bc_next='bc_t'//trim(filename_bc_next)
    do ilev=1,refine_max_level 
       ! 1st time snapshot 
       open(42,file='bc/'//trim(filename_bc_old)//'_AMR'//trim(str(ilev))//'.dat')
        ! has to adjust according to dimension!!!
       do i=1,domain_nx2*(2**(ilev-1))
          do j=1,domain_nx3*(2**(ilev-1))
             read(42,*) myBC(i,j,1,ilev)
          enddo
       enddo
      close(42)
      ! 2nd time snapshot 
      open(42,file='bc/'//trim(filename_bc_next)//'_AMR'//trim(str(ilev))//'.dat')
      do i=1,domain_nx2*(2**(ilev-1))
          do j=1,domain_nx3*(2**(ilev-1))
             read(42,*) myBC(i,j,2,ilev)
          enddo
      enddo
      close(42)
    enddo
    !! can check for time-independent part here: for up to t=1 then fixed
    myBC(:,:,2,1:refine_max_level)=myBC(:,:,1,1:refine_max_level) 

  end subroutine initglobaldata_usr

  subroutine initonegrid_usr(ixG^L,ix^L,w,x)

    ! initialize one grid 

    integer, intent(in) :: ixG^L, ix^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    w(ix^S,rho_)     = 0.5d0

  end subroutine initonegrid_usr

  ! special boundary types, user defined
  ! user must assign conservative variables in bounderies
  subroutine specialbound_usr(qt,ixG^L,ixO^L,iB,w,x)
    integer, intent(in) :: ixG^L, ixO^L, iB
    double precision, intent(in) :: qt, x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    integer :: ig^D, lvl, ix1, ix^L
    double precision :: dx2ratio

    select case(iB)
    case(1)
       ! === Left inflow boundary ===!
       dx2ratio=dx(2,1)/(x(ixOmin1,ixOmin2+1,ixOmin3,2)-x(ixOmin1,ixOmin2,ixOmin3,2))
       lvl=1+nint(dlog(dx2ratio)/dlog(2.0d0)) 
       ^D&ig^D=1+floor((x(ixOmin1,ixOmin2,ixOmin3,^D)-xprobmin^D)/(dx(^D,lvl)*block_nx^D))\
       !print *,'**********************'
       !print *,lvl,(x(ixOmin1,ixOmin2+1,ixOmin3,2)-x(ixOmin1,ixOmin2,ixOmin3,2))
       !print *,lvl,dx(2,1),dx(2,1)/(x(ixOmin1,ixOmin2+1,ixOmin3,2)-x(ixOmin1,ixOmin2,ixOmin3,2)),dxlevel(2)
       ! note: can not use dxlevel here, may not conform due to coarse bc call
       ^D&ixmin^D=(ig^D-1)*block_nx^D+ixOmin^D-nghostcells\
       ^D&ixmax^D=(ig^D-1)*block_nx^D+ixOmax^D-nghostcells\
       !print *,ixmin2,ixmax2,ixmin3,ixmax3
       !print *,ixOmin1,ixOmax1,ixOmin2,ixOmax2,ixOmin3,ixOmax3
       ! note: may not be entire block size to be filled!!!
       !!if(ixOmax2-ixOmin2+1/=block_nx2.or.ixOmax3-ixOmin3+1/=block_nx3) call mpistop('mismatch')
       ! Linear interpolation between 2 successive time snapshots
       do ix1=ixOmin1,ixOmax1
          w(ix1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)=&
               myBC(ixmin2:ixmax2,ixmin3:ixmax3,1,lvl)+&
                     ((qt-time_old_bc)/(time_next_bc-time_old_bc))*&
              (myBC(ixmin2:ixmax2,ixmin3:ixmax3,2,lvl)-&
               myBC(ixmin2:ixmax2,ixmin3:ixmax3,1,lvl))
       enddo
    case default
       call mpistop('boundary not defined')
    end select

  end subroutine specialbound_usr


  subroutine specialrefine_grid(igrid,level,ixG^L,ix^L,qt,w,x,refine,coarsen)
    ! Enforce additional refinement or coarsening
    ! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.
    ! you must set consistent values for integers refine/coarsen:
    ! refine = -1 enforce to not refine
    ! refine =  0 doesn't enforce anything
    ! refine =  1 enforce refinement
    ! coarsen = -1 enforce to not coarsen
    ! coarsen =  0 doesn't enforce anything
    ! coarsen =  1 enforce coarsen
    integer, intent(in) :: igrid, level, ixG^L, ix^L
    double precision, intent(in) :: qt, w(ixG^S,1:nw), x(ixG^S,1:ndim)
    integer, intent(inout) :: refine, coarsen

    ! test with different levels of refinement enforced
    if (all(x(ix^S,2) < zero) .and. all(x(ix^S,3)>zero)) then
       refine=-1
    else
       refine=1
    endif

  end subroutine specialrefine_grid

  ! This subroutine is called at the beginning of each time step 
  ! by each processor. No communication is specified, so the user
  ! has to implement MPI routines if information has to be shared
  subroutine process_global_usr(iit,qt)
    
    integer, intent(in)          :: iit
    double precision, intent(in) :: qt
    integer :: dt_bc, i, j,ilev

    if (qt>time_next_bc) then
       ! The snapshot called "next" becomes "old" 
       myBC(:,:,1,1:refine_max_level)=myBC(:,:,2,1:refine_max_level)
       dt_bc=time_next_bc-time_old_bc ! to save dt, CONSTANT
       time_old_bc=time_next_bc
       filename_bc_old=filename_bc_next
       ! Read the new BC data file
       time_next_bc=time_old_bc+dt_bc
       write (filename_bc_next,'(I2.2)') time_next_bc
       filename_bc_next='bc_t'//trim(filename_bc_next)
       do ilev=1,refine_max_level 
          ! next time snapshot 
          open(42,file='bc/'//trim(filename_bc_next)//'_AMR'//trim(str(ilev))//'.dat')
          ! has to adjust according to dimension!!!
          do i=1,domain_nx2*2**(ilev-1)
             do j=1,domain_nx3*2**(ilev-1)
               read(42,*) myBC(i,j,2,ilev)
             enddo
          enddo
          close(42)
       enddo
    endif

  end subroutine process_global_usr

end module mod_usr

