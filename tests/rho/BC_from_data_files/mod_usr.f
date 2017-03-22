module mod_usr
  use mod_rho
  implicit none
  double precision, allocatable, save :: myBC1(:,:,:), myBC2(:,:,:)
  integer :: time_old_bc, time_next_bc !only integer time allowed between snapshots
  character(len=50) :: filename_bc_old, filename_bc_next

contains

  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods

    usr_set_parameters => initglobaldata_usr
    usr_init_one_grid  => initonegrid_usr
    usr_special_bc     => specialbound_usr
    usr_process_global => process_global_usr
    usr_refine_grid    => specialrefine_grid

    call rho_activate()

  end subroutine usr_init

  subroutine initglobaldata_usr
    use mod_global_parameters
    integer :: i, j
    character(len=2) :: t0, t1

    allocate(myBC1(domain_nx2,domain_nx3,2))
    allocate(myBC2(2*domain_nx2,2*domain_nx3,2))
    time_old_bc =0
    time_next_bc=1
    write (filename_bc_old,'(I2.2)') time_old_bc
    write (filename_bc_next,'(I2.2)') time_next_bc
    filename_bc_old ='bc_t'//trim(filename_bc_old)
    filename_bc_next='bc_t'//trim(filename_bc_next)
    ! 1st time snapshot / 1st lvl of AMR
    open(42,file='bc/'//trim(filename_bc_old)//'_AMR1.dat')
    do i=1,domain_nx2
       do j=1,domain_nx3
          read(42,*) myBC1(i,j,1)
       enddo
    enddo
    close(42)
    ! 1st time snapshot / 2nd lvl of AMR
    open(42,file='bc/'//trim(filename_bc_old)//'_AMR2.dat')
    do i=1,domain_nx2*2
       do j=1,domain_nx3*2
          read(42,*) myBC2(i,j,1)
       enddo
    enddo
    close(42)
    ! 2nd time snapshot / 1st lvl of AMR
    open(42,file='bc/'//trim(filename_bc_next)//'_AMR1.dat')
    do i=1,domain_nx2
       do j=1,domain_nx3
          read(42,*) myBC1(i,j,2)
       enddo
    enddo
    close(42)
    ! 2nd time snapshot / 2nd lvl of AMR
    open(42,file='bc/'//trim(filename_bc_next)//'_AMR2.dat')
    do i=1,domain_nx2*2
       do j=1,domain_nx3*2
          read(42,*) myBC2(i,j,2)
       enddo
    enddo
    close(42)

    ! to verify orientation of the indexes
!!$    print*, myBC1(1,2), myBC1(2,2)
!!$    print*, myBC1(1,1), myBC1(2,1)
!!$    print*, ' '
!!$    print*, myBC2(1,2), myBC2(2,2)
!!$    print*, myBC2(1,1), myBC2(2,1)

  end subroutine initglobaldata_usr

  subroutine initonegrid_usr(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
     ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w,x)

    ! initialize one grid 

    use mod_global_parameters

    integer, intent(in) :: ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
        ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
    double precision, intent(in) :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:nw)

    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rho_)     = 0.5d0

  end subroutine initonegrid_usr

 ! special boundary types, user defined
  ! user must assign conservative variables in bounderies
  subroutine specialbound_usr(qt,ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,&
     ixGmax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,iB,w,x)
    use mod_global_parameters

    integer, intent(in) :: ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, iB
    double precision, intent(in) :: qt, x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:nw)

    logical :: patchw(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3)
    double precision:: rinlet(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3)

    integer :: ig2tmp, ig3tmp, lvl

    select case(iB)
    case(1)
       ! === Left boundary ===!
       lvl=int(1.01*(((xprobmax2-xprobmin2)/real(domain_nx2))/(x(ixOmin1,&
          ixOmin2+1,ixOmin3,2)-x(ixOmin1,ixOmin2,ixOmin3,2))))
!       print*, lvl
!!$       print*, dxlevel(1), dxlevel(2), dxlevel(3)
!!$       print*, x(ixOmin1+1,ixOmin2,ixOmin3,1)-x(ixOmin1,ixOmin2,ixOmin3,1)
!!$       print*, x(ixOmin1,ixOmin2+1,ixOmin3,2)-x(ixOmin1,ixOmin2,ixOmin3,2)
!!$       print*, x(ixOmin1,ixOmin2,ixOmin3+1,3)-x(ixOmin1,ixOmin2,ixOmin3,3)
       ig2tmp=1+int((x(ixOmin1,ixOmin2,ixOmin3,&
          2)-xprobmin2)/((xprobmax2-xprobmin2)/(real(domain_nx2)/real(&
          block_nx2/(2**(lvl-1))))))
       ig3tmp=1+int((x(ixOmin1,ixOmin2,ixOmin3,&
          3)-xprobmin3)/((xprobmax3-xprobmin3)/(real(domain_nx3)/real(&
          block_nx3/(2**(lvl-1))))))
       ! Linear interpolation between 2 successive time snapshots
       if (lvl==1) then          
          w(ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             rho_)=myBC1(ixGmin2+(ig2tmp-1)*block_nx2:ixGmax2-2*nghostcells+&
             (ig2tmp-1)*block_nx2,ixGmin3+&
             (ig3tmp-1)*block_nx3:ixGmax3-2*nghostcells+(ig3tmp-1)*block_nx3,&
             1)+((qt-time_old_bc)/(time_next_bc-time_old_bc))*(myBC1(ixGmin2+&
             (ig2tmp-1)*block_nx2:ixGmax2-2*nghostcells+(ig2tmp-1)*block_nx2,&
             ixGmin3+(ig3tmp-1)*block_nx3:ixGmax3-2*nghostcells+&
             (ig3tmp-1)*block_nx3,2)-myBC1(ixGmin2+&
             (ig2tmp-1)*block_nx2:ixGmax2-2*nghostcells+(ig2tmp-1)*block_nx2,&
             ixGmin3+(ig3tmp-1)*block_nx3:ixGmax3-2*nghostcells+&
             (ig3tmp-1)*block_nx3,1))
       elseif (lvl==2) then
          w(ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             rho_)=myBC2(ixGmin2+(ig2tmp-1)*block_nx2:ixGmax2-2*nghostcells+&
             (ig2tmp-1)*block_nx2,ixGmin3+&
             (ig3tmp-1)*block_nx3:ixGmax3-2*nghostcells+(ig3tmp-1)*block_nx3,&
             1)+((qt-time_old_bc)/(time_next_bc-time_old_bc))*(myBC2(ixGmin2+&
             (ig2tmp-1)*block_nx2:ixGmax2-2*nghostcells+(ig2tmp-1)*block_nx2,&
             ixGmin3+(ig3tmp-1)*block_nx3:ixGmax3-2*nghostcells+&
             (ig3tmp-1)*block_nx3,2)-myBC2(ixGmin2+&
             (ig2tmp-1)*block_nx2:ixGmax2-2*nghostcells+(ig2tmp-1)*block_nx2,&
             ixGmin3+(ig3tmp-1)*block_nx3:ixGmax3-2*nghostcells+&
             (ig3tmp-1)*block_nx3,1))

!!$       elseif (lvl==3) then
!!$       w(ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)=myBC1(ixGmin2+(ig2tmp-1)*block_nx2:ixGmax2-2*nghostcells+(ig2tmp-1)*block_nx2,ixGmin3+(ig3tmp-1)*block_nx3:ixGmax3-2*nghostcells+(ig3tmp-1)*block_nx3)
       endif
       !print*, ixOmin2-nghostcells, ixOmin2, ixOmax2-nghostcells, ixOmax2
       !print*, ixGmin2-nghostcells, ixGmin2, ixGmax2-nghostcells, ixGmax2
          

       !w(ixOmin1,:,:,rho_)=w(ixOmax1,:,:,rho_)
!!$       where (x(ixO^S,2)>zero .and. x(ixO^S,3)>zero)
!!$          w(ixO^S,rho_) = 0.5d0
!!$       endwhere
!!$       where (x(ixO^S,2)<zero .and. x(ixO^S,3)>zero)
!!$          w(ixO^S,rho_) = one     
!!$       endwhere     
!!$       where (x(ixO^S,2)>zero .and. x(ixO^S,3)<zero)
!!$          w(ixO^S,rho_) = 1.5d0
!!$       endwhere
!!$       where (x(ixO^S,2)<zero .and. x(ixO^S,3)<zero)
!!$          w(ixO^S,rho_) = two
!!$       endwhere

    case default
       call mpistop('boundary not defined')
    end select

  end subroutine specialbound_usr


  subroutine specialrefine_grid(igrid,level,ixGmin1,ixGmin2,ixGmin3,ixGmax1,&
     ixGmax2,ixGmax3,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,qt,w,x,refine,&
     coarsen)
    ! Enforce additional refinement or coarsening
    ! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.
    ! you must set consistent values for integers refine/coarsen:
    ! refine = -1 enforce to not refine
    ! refine =  0 doesn't enforce anything
    ! refine =  1 enforce refinement
    ! coarsen = -1 enforce to not coarsen
    ! coarsen =  0 doesn't enforce anything
    ! coarsen =  1 enforce coarsen
    use mod_global_parameters

    integer, intent(in) :: igrid, level, ixGmin1,ixGmin2,ixGmin3,ixGmax1,&
       ixGmax2,ixGmax3, ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
    double precision, intent(in) :: qt, w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:nw), x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:ndim)
    integer, intent(inout) :: refine, coarsen

    ! So as we can test the code with different levels of refinment enforced
    if (all(x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
       2) < zero) .and. all(x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
       3)>zero)) then
       refine=-1
    else
       refine=1
    endif

  end subroutine specialrefine_grid

  ! This subroutine is called at the beginning of each time step 
  ! by each processor. No communication is specified, so the user
  ! has to implement MPI routines if information has to be shared
  subroutine process_global_usr(iit,qt)
    
    use mod_global_parameters
    
    integer, intent(in)          :: iit
    double precision, intent(in) :: qt
    integer :: dt_bc, i, j

    if (qt>time_next_bc) then
       print*, 'switch'
       ! The snapshot called "next" becomes "old" 
       myBC1(:,:,1)=myBC1(:,:,2)
       myBC2(:,:,1)=myBC2(:,:,2)
       dt_bc=time_next_bc-time_old_bc ! to save dt, CONSTANT
       time_old_bc=time_next_bc
       filename_bc_old=filename_bc_next
       ! Read the new BC data file
       time_next_bc=time_old_bc+dt_bc
       write (filename_bc_next,'(I2.2)') time_next_bc
       filename_bc_next='bc_t'//trim(filename_bc_next)
       ! 2nd time snapshot / 1st lvl of AMR
       open(42,file='bc/'//trim(filename_bc_next)//'_AMR1.dat')
       do i=1,domain_nx2
          do j=1,domain_nx3
             read(42,*) myBC1(i,j,2)
          enddo
       enddo
       close(42)
       ! 2nd time snapshot / 2nd lvl of AMR
       open(42,file='bc/'//trim(filename_bc_next)//'_AMR2.dat')
       do i=1,domain_nx2*2
          do j=1,domain_nx3*2
             read(42,*) myBC2(i,j,2)
          enddo
       enddo
       close(42)        

    endif

  end subroutine process_global_usr

end module mod_usr

