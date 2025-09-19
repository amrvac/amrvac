module mod_usr
  use mod_amrvac
  use mod_physics

  implicit none

  double precision :: rhodens  = 10.0d0
  double precision :: rholight = 1.0d0
  double precision :: kx       = 4.0d0 * acos(-1.0d0)
  double precision :: kz       = 32d0*atan(1.0d0)
  double precision :: pint     = 2.5d0
  double precision :: dlin     = 0.025d0
  double precision :: sigma    = 0.00125d0

contains

  subroutine usr_init()


    call set_coordinate_system("Cartesian_3D")

!todojesse this reference does not work for some reason ...
    usr_refine_grid => specialrefine_grid

    usr_init_one_grid => initonegrid_usr

    usr_before_main_loop => blockcheck_usr

    call phys_activate()

  end subroutine usr_init

  subroutine blockcheck_usr()
    integer:: igrid, iigrid, iw, ix1, ix2, ix3

    print *, 'checking blocks'
    
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       do iw = 1, nw
          do ix3 = ixGlo3, ixGhi3
             do ix2 = ixGlo2, ixGhi2
                do ix1 = ixGlo1, ixGhi1
!                   print *, ix1,ix2,ix3,iw,igrid,ps(igrid)%w(ix1,ix2,ix3,iw)
                   if (ps(igrid)%w(ix1,ix2,ix3,iw) /= ps(igrid)%w(ix1,ix2,ix3,iw)) then
                      print *, 'NaN in block: ', ix1,ix2,ix3,iw,igrid
                   end if
                end do
             end do
          end do
       end do
    end do

  end subroutine blockcheck_usr

  subroutine initonegrid_usr(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
     ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w,x)
    integer, intent(in)             :: ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,&
       ixGmax3, ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
    double precision, intent(in)    :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:nw)
    double precision                :: vextra

    select case(iprob)
    case(1)
      vextra=0.0d0
    case(2)
      vextra = 10.0d0
    case(3)
      vextra=100.0d0
    case default
      error stop "Unknown iprob"
    endselect

    where(x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,&
       2)<0.25d0.or.x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,&
       2)>0.75d0)
      w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,rho_)=rholight
    elsewhere
      w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,rho_)=rhodens
    endwhere

    w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,p_)=pint

    where((x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,&
       2)<=(0.25d0-dlin)).or.(x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,2)>=(0.75d0+dlin)))
      w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,mom(1))=vextra-0.5d0
    endwhere

    where((x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,&
       2)>=(0.25d0+dlin)).and.(x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,2)<=(0.75d0-dlin)))
      w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,mom(1))=vextra+0.5d0
    endwhere

    where((x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,&
       2)>(0.25d0-dlin)).and.(x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,2)<(0.25d0+dlin)))
      w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,&
         mom(1))=(x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,&
         2)-(0.25d0-dlin))* ((vextra+0.5d0)-(vextra-0.5d0))/(2.0d0*dlin)+&
         (vextra-0.5d0)
    endwhere

    where((x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,&
       2)>(0.75d0-dlin)).and.(x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,2)<(0.75d0+dlin)))
      w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,&
         mom(1))=(x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,&
         2)-(0.75d0-dlin))* ((vextra-0.5d0)-(vextra+0.5d0))/(2.0d0*dlin)+&
         (vextra+0.5d0)
    endwhere

    w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,&
       mom(2))=0.01d0*dsin(kx*x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1))* (dexp(-0.5d0*(x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,2)-0.25d0)**2/sigma)+dexp(-0.5d0*(x(ixGmin1:ixGmax1,&
       ixGmin2:ixGmax2,ixGmin3:ixGmax3,2)-0.75d0)**2/sigma))
    
    w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,&
       mom(3))=0.1d0*dsin(kz*x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,&
       3))* (dexp(-0.5d0*(x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,&
       2)-0.25d0)**2/sigma)+dexp(-0.5d0*(x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,2)-0.75d0)**2/sigma))
   

    call phys_to_conserved(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
       ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w,x)

  end subroutine initonegrid_usr

  subroutine specialrefine_grid(igrid,level,ixGmin1,ixGmin2,ixGmin3,&
    ixGmax1,ixGmax2,ixGmax3,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,&
    qt,w,x,refine,coarsen)

    use mod_global_parameters

    ! Enforce additional refinement or coarsening
    ! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.

    ! you must set consistent values for integers refine/coarsen:

    ! refine = -1 enforce to not refine
    ! refine =  0 doesn't enforce anything
    ! refine =  1 enforce refinement

    ! coarsen = -1 enforce to not coarsen
    ! coarsen =  0 doesn't enforce anything
    ! coarsen =  1 enforce coarsen

    integer, intent(in)             :: igrid, level, ixGmin1,ixGmin2,&
        ixGmin3,ixGmax1,ixGmax2,ixGmax3, ixmin1,ixmin2,ixmin3,ixmax1,&
        ixmax2,ixmax3
    double precision, intent(in)    :: qt, x(ixGmin1:ixGmax1,&
         ixGmin2:ixGmax2,ixGmin3:ixGmax3,1:ndim)
    double precision, intent(in)    :: w(ixGmin1:ixGmax1,&
         ixGmin2:ixGmax2,ixGmin3:ixGmax3,1:nw)
    integer, intent(inout) :: refine, coarsen

    if ( any(x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,2) & 
         <0.75d0) .and. any(x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,2) & 
         >0.75d0) &
         .or. &
         any(x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,2) & 
         <0.25d0) .and. any(x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,2) & 
         >0.25d0)) then
       coarsen = -1
       refine  = 1
    else
       coarsen = 0
       refine  = 0
    end if

  end subroutine specialrefine_grid

end module mod_usr
