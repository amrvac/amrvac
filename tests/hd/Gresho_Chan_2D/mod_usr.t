module mod_usr
  use mod_hd

  implicit none

  double precision :: v0^D, Mach

  integer, parameter :: dp = kind(0.0d0)
  integer :: i_err_p, i_err_r
  integer :: i_sol

contains

  !> Read this module's parameters from a file
  subroutine usr_params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /usr_list/ v0^D, Mach

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, usr_list, end=111)
111    close(unitpar)
    end do

  end subroutine usr_params_read

  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods
    use mod_variables

    call usr_params_read(par_files)

    usr_init_one_grid => GC_init_one_grid
    usr_process_grid => set_error
    usr_print_log => print_error
    usr_refine_grid    => specialrefine_grid

    call set_coordinate_system("Cartesian_2D")
    call hd_activate()

    i_err_r = var_set_extravar("rho_error", "rho_error")
    i_err_p = var_set_extravar("p_error", "p_error")
    i_sol = var_set_extravar("p_solution", "p_solution")

  end subroutine usr_init

  ! Initialize one grid
  subroutine GC_init_one_grid(ixG^L,ix^L,w,x)
    use mod_global_parameters

    integer, intent(in) :: ixG^L, ix^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)
    
    double precision :: rr(ixG^S),vphi(ixG^S),cosphi(ixG^S),sinphi(ixG^S)
    double precision :: xmid^D

    {^IFONED   call mpistop("This is a 2D HD problem") }
    {^IFTHREED call mpistop("This is a 2D HD problem") }

    ^D&xmid^D=xprobmin^D+0.5d0*(xprobmax^D-xprobmin^D);
    rr(ix^S)=dsqrt( {^D&(x(ix^S,^D)-xmid^D)**2|+} )
    w(ix^S,rho_) = 1.0d0
    {^IFTWOD
    w(ix^S,e_)   = p_solution(x(ix^S, 1), x(ix^S, 2), xmid1, xmid2)
    vphi(ix^S)   = vphi_solution(x(ix^S, 1), x(ix^S, 2), xmid1, xmid2)
    sinphi(ix^S)=(x(ix^S,2)-xmid2)/rr(ix^S)
    cosphi(ix^S)=(x(ix^S,1)-xmid1)/rr(ix^S)
    w(ix^S,mom(1)) = v01-vphi(ix^S)*sinphi(ix^S)
    w(ix^S,mom(2)) = v02+vphi(ix^S)*cosphi(ix^S)
    }

    call hd_to_conserved(ixG^L,ix^L,w,x)

  end subroutine GC_init_one_grid

  elemental function p_solution(x, y, xc, yc) result(val)
    use mod_global_parameters
    real(dp), intent(in) :: x, y, xc, yc
    real(dp)             :: val
    real(dp) :: p0,rad

    p0=1.0d0/(hd_gamma*Mach**2)
    rad=dsqrt((x-xc)**2+(y-yc)**2)
    if(rad<0.2d0)then
      val=p0+12.5d0*rad**2
    else
      if(rad<0.4d0)then
        val=p0+12.5d0*rad**2+4.0d0*(1.0d0-5.0d0*rad-dlog(0.2d0)+dlog(rad))
      else
        val=p0-2.0d0+4.0d0*dlog(2.0d0)
      endif
    endif
  end function p_solution

  elemental function vphi_solution(x, y, xc, yc) result(val)
    use mod_global_parameters
    real(dp), intent(in) :: x, y, xc, yc
    real(dp)             :: val
    real(dp) :: rad

    rad=dsqrt((x-xc)**2+(y-yc)**2)
    if(rad<0.2d0)then
      val=5.0d0*rad
    else
      if(rad<0.4d0)then
        val=2.0d0-5.0d0*rad
      else
        val=0.0d0
      endif
    endif
  end function vphi_solution

  subroutine set_error(igrid,level,ixI^L,ixO^L,qt,w,x)
    use mod_global_parameters
    integer, intent(in)             :: igrid,level,ixI^L,ixO^L
    double precision, intent(in)    :: qt,x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    integer :: ic1,ic2
    double precision:: xmid^D, rr(ixI^S), xshift^D

    w(ixO^S,i_sol) =1.0d0/(hd_gamma*Mach**2)-2.0d0+4.0d0*dlog(2.0d0)
    ! determine centers of all 9 potentially overlapping circles
    if(v01>=0.0d0)then
       xshift1=v01*qt-floor(qt*v01/(xprobmax1-xprobmin1))*(xprobmax1-xprobmin1)
    else
       xshift1=v01*qt+floor(qt*dabs(v01)/(xprobmax1-xprobmin1))*(xprobmax1-xprobmin1)
    endif
    if(v02>=0.0d0)then
       xshift2=v02*qt-floor(qt*v02/(xprobmax2-xprobmin2))*(xprobmax2-xprobmin2)
    else
       xshift2=v02*qt+floor(qt*dabs(v02)/(xprobmax2-xprobmin2))*(xprobmax2-xprobmin2)
    endif
    do ic1=-1,1
      xmid1=xprobmin1+(ic1+0.5d0)*(xprobmax1-xprobmin1)+xshift1;
      do ic2=-1,1
         xmid2=xprobmin2+(ic2+0.5d0)*(xprobmax2-xprobmin2)+xshift2;
         rr(ixO^S)=dsqrt( {^D&(x(ixO^S,^D)-xmid^D)**2|+} )
         where(rr(ixO^S)<0.4d0)
           w(ixO^S,i_sol) = p_solution(x(ixO^S, 1), x(ixO^S, 2),xmid1,xmid2)
         endwhere
      enddo
    enddo
    call hd_to_primitive(ixI^L,ixO^L,w,x)
    w(ixO^S,i_err_p) = w(ixO^S,e_) - w(ixO^S,i_sol)
    call hd_to_conserved(ixI^L,ixO^L,w,x)
    w(ixO^S,i_err_r) = w(ixO^S,rho_) - 1.0d0
  end subroutine set_error

  subroutine print_error()
    use mod_global_parameters
    use mod_input_output, only: get_volume_average
    double precision   :: modes(nw, 2), volume

    call get_volume_average(1, modes(:, 1), volume)
    call get_volume_average(2, modes(:, 2), volume)

    if (mype == 0) then
       write(*, "(A,3E16.8)") " time -- RMSE (rho/p):",  &
           global_time, sqrt(modes(i_err_r, 2)), sqrt(modes(i_err_p, 2))
    end if
  end subroutine print_error

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
    use mod_global_parameters

    integer, intent(in) :: igrid, level, ixG^L, ix^L
    double precision, intent(in) :: qt, w(ixG^S,1:nw), x(ixG^S,1:ndim)
    integer, intent(inout) :: refine, coarsen
    double precision:: xmid^D, rr(ixG^S)
    logical:: mask(ixG^S)

    ^D&xmid^D=xprobmin^D+0.5d0*(xprobmax^D-xprobmin^D);
    rr(ix^S)=dsqrt( {^D&(x(ix^S,^D)-xmid^D)**2|+} )
    ! test with different levels of refinement enforced
    if (qt<smalldouble.and.any((rr(ix^S)<0.3d0).and.(rr(ix^S)>0.1d0))) then
       refine=1
    endif
    mask(ix^S)=dabs(w(ix^S,rho_)-1.0d0)>1.0d-3
    if (any(mask(ix^S))) then
       refine=1
    endif

  end subroutine specialrefine_grid

end module mod_usr
