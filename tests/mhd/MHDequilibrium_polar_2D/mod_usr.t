!> Equilibrium setup for purely planar 2D polar
! implemented a balance between toroidal field tension and centrifugal force
! combined with constant total pressure
! embedded in a uniform, unmagnetized medium
! This tests our implementation of current evaluations (typecurl), among other things.
module mod_usr
  use mod_mhd

  implicit none

  ! Input values for 
  !    radius: Rjet
  !    density and B factors: rhojet, rhocloud, apar, pjet, Bazi
  double precision :: Rjet, rhojet, rhocloud, apar, pjet, Bazi


  integer, parameter :: dp = kind(0.0d0)
  integer :: i_err_p, i_err_r, i_err_b
  integer :: i_sol_p, i_sol_r, i_sol_b, i_totp

contains

  !> Read this module's parameters from a file
  subroutine usr_params_read(files)
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /usr_list/ rhojet, rhocloud, Rjet, apar, pjet, Bazi

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, usr_list, end=111)
111    close(unitpar)
    end do

  end subroutine usr_params_read

  subroutine usr_init()
    use mod_variables

    call usr_params_read(par_files)

    usr_set_parameters  => initglobaldata_usr
    usr_init_one_grid => CCC_init_one_grid
    usr_process_grid => set_error
    usr_print_log => print_error
    usr_refine_grid    => specialrefine_grid
    usr_aux_output      => specialvar_output
    usr_add_aux_names   => specialvarnames_output

    {^IFONED call mpistop("2D polar problem") }
    {^IFTWOD call set_coordinate_system("polar_2D") }
    {^IFTHREED call mpistop("2D polar problem") }

    call mhd_activate()

    i_err_r = var_set_extravar("rho_error", "rho_error")
    i_sol_r = var_set_extravar("rho_solution", "rho_solution")
    i_err_p = var_set_extravar("p_error", "p_error")
    i_sol_p = var_set_extravar("p_solution", "p_solution")
    i_err_b = var_set_extravar("b_error", "b_error")
    i_sol_b = var_set_extravar("b_solution", "b_solution")
    i_totp  = var_set_extravar("totp", "totp")

  end subroutine usr_init

  subroutine initglobaldata_usr
    character(len=20) :: printsettingformat

    printsettingformat='(1x,A50,ES15.7,A7)'
    
    if(mype==0) then
      write(*,*) "Jet setup:"
      write(*,printsettingformat) "density in jet ",rhojet," input"
      write(*,printsettingformat) "density in cloud ",rhocloud," input"
      write(*,printsettingformat) "pressure jet ",pjet," input"
      write(*,printsettingformat) "Jet azimuthal B strength ",Bazi," input"
    end if

    if(mype==0) then
      write(*,*) "Deduced dimensionless values:"
      write(*,printsettingformat) "density ratio jet/cloud ",rhojet/rhocloud," output"
    end if


  end subroutine initglobaldata_usr

  !> Initialize one grid
  subroutine CCC_init_one_grid(ixG^L,ix^L,w,x)
    integer, intent(in) :: ixG^L, ix^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)
    
    {^IFONED   call mpistop("This is a 2D polar MHD problem") }

    {^NOONED
    w(ix^S,rho_) = rho_solution(x(ix^S, 1))
    w(ix^S,mag(1)) = 0.0d0
    w(ix^S,mag(2)) = bphi_solution(x(ix^S, 1))
    w(ix^S,mom(1)) = 0.0d0
    w(ix^S,mom(2)) = w(ix^S,mag(2))/dsqrt(w(ix^S,rho_))
    w(ix^S,e_)   = pjet-0.5d0*w(ix^S,mag(2))**2.0d0
    }

    call mhd_to_conserved(ixG^L,ix^L,w,x)

  end subroutine CCC_init_one_grid

  elemental function rho_solution(rad) result(val)
    real(dp), intent(in) :: rad
    real(dp)             :: val

    if(rad<Rjet)then
      val=rhojet
    else
      val=rhocloud
    endif

  end function rho_solution

  elemental function bphi_solution(rad) result(val)
    real(dp), intent(in) :: rad
    real(dp)             :: val

    if(rad<Rjet)then
      val=Bazi*TANH(rad/apar)
    else
      val=0.0d0
    endif

  end function bphi_solution


  subroutine set_error(igrid,level,ixI^L,ixO^L,qt,w,x)
    integer, intent(in)             :: igrid,level,ixI^L,ixO^L
    double precision, intent(in)    :: qt,x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    w(ixO^S,i_sol_r) = rho_solution(x(ixO^S, 1))
    w(ixO^S,i_err_r) = w(ixO^S,rho_) - w(ixO^S,i_sol_r)
    w(ixO^S,i_sol_p) = pjet-0.5d0*bphi_solution(x(ixO^S, 1))**2.0d0
    w(ixO^S,i_sol_b) =dsqrt(bphi_solution(x(ixO^S, 1))**2)
    call mhd_to_primitive(ixI^L,ixO^L,w,x)
    w(ixO^S,i_err_p) = w(ixO^S,e_) - w(ixO^S,i_sol_p)
    {^NOONED
    w(ixO^S,i_totp) = w(ixO^S,e_) +0.5d0* &
       (w(ixO^S,mag(1))**2+w(ixO^S,mag(2))**2)
    }
    call mhd_to_conserved(ixI^L,ixO^L,w,x)
    {^NOONED
    w(ixO^S,i_err_b) = dsqrt( w(ixO^S,mag(1))**2+w(ixO^S,mag(2))**2 ) &
                          - w(ixO^S,i_sol_b)
    }
  end subroutine set_error

  subroutine print_error()
    use mod_input_output, only: get_volume_average
    double precision   :: modes(nw, 2), volume

    call get_volume_average(1, modes(:, 1), volume)
    call get_volume_average(2, modes(:, 2), volume)

    if (mype == 0) then
       write(*, "(A,4E16.8)") " CONVTEST (global_time,rho_2,p_2,b_2):",  &
           global_time, sqrt(modes(i_err_r, 2)), sqrt(modes(i_err_p, 2)), sqrt(modes(i_err_b, 2))
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
    integer, intent(in) :: igrid, level, ixG^L, ix^L
    double precision, intent(in) :: qt, w(ixG^S,1:nw), x(ixG^S,1:ndim)
    integer, intent(inout) :: refine, coarsen

    ! test with different levels of refinement enforced
    if (qt<smalldouble.and.any((x(ix^S,1)<1.1d0).and.(x(ix^S,1)>0.9d0))) then
       refine=1
    endif

  end subroutine specialrefine_grid

  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
  ! this subroutine can be used in convert, to add auxiliary variables to the
  ! converted output file, for further analysis using tecplot, paraview, ....
  ! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
  ! the array normconv can be filled in the (nw+1:nw+nwauxio) range with
  ! corresponding normalization values (default value 1)
    use mod_radiative_cooling
    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)

    double precision :: pth(ixI^S),B2(ixI^S)
    double precision :: divb(ixI^S),wlocal(ixI^S,1:nw)
    double precision :: Btotal(ixI^S,1:ndir),current(ixI^S,7-2*ndir:3)
    integer :: idirmin,idir

    wlocal(ixI^S,1:nw)=w(ixI^S,1:nw)
    ! store current
    call get_current(wlocal,ixI^L,ixO^L,idirmin,current)
    ! in 2D polar, only z-component of current exists
    w(ixO^S,nw+1)=current(ixO^S,idirmin)

  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
    ! newly added variables need to be concatenated with the w_names/primnames string
    character(len=*) :: varnames

    varnames='j3'

  end subroutine specialvarnames_output

end module mod_usr
