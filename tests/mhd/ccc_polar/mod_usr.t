!> Equilibrium setups of a current-carrying cylinder
! implemented a choice between three 1D cylindrical equilibria,
! embedded in a uniform, (un)magnetized medium
! This tests our implementation of current evaluations (typecurl), among other things.
module mod_usr
  use mod_mhd

  implicit none

  ! Input values for 
  !    choice of equilibrium: equilibrium_version
  !    dimensionless parameters: beta1, qfac1, Rvacs, drat, invbext
  double precision :: beta1, qfac1, Rvacs, drat, invbext
  double precision :: pr01, Bz0, rho0, Jfac0, Lz
  character(len=std_len) :: equilibrium_version

  integer, parameter :: dp = kind(0.0d0)
  integer :: i_err_p, i_err_r, i_err_b
  integer :: i_sol_p, i_sol_r, i_sol_b, i_totp

contains

  !> Read this module's parameters from a file
  subroutine usr_params_read(files)
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /usr_list/ beta1, qfac1,  &
                 Rvacs, drat, invbext, equilibrium_version

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, usr_list, end=111)
111    close(unitpar)
    end do

    equilibrium_version=trim(equilibrium_version)
    
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

    {^IFTWOD call set_coordinate_system("polar_2.5D") }
    {^IFTHREED call set_coordinate_system("polar_3D")  }

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
    double precision :: rootfac, alfa

    printsettingformat='(1x,A50,ES15.7,A7)'
    
    {^IFTHREED
    Lz=(xprobmax3-xprobmin3)
    }
    {^IFTWOD
    Lz=xprobmax1
    }
    if(mype==0) then
      write(*,*) "Equilibrium chosen:", equilibrium_version
      write(*,printsettingformat) "beta at r=1 ",beta1," input"
      write(*,printsettingformat) "q-factor at r=1 ",qfac1," input"
      write(*,printsettingformat) "density contrast in flux tube ",drat," input"
      write(*,printsettingformat) "ratio alfven speed at r=0 to external sound speed ",Rvacs," input"
      write(*,printsettingformat) "inverse beta external ",invbext," input"
      write(*,printsettingformat) "using Lz equal to ",Lz," input"
    end if

    select case(equilibrium_version)
       case('TokamakCurrent')
          pr01=beta1*(1.0d0+invbext)/((beta1+1.0d0)*mhd_gamma)
          Jfac0=6.0d0*dsqrt(2.0d0*(1.0d0+invbext)/((beta1+1.0d0)*mhd_gamma))/dsqrt(1.0d0+(Lz*qfac1/dpi)**2)
          rootfac=2.0d0*(1.0d0+invbext)/((beta1+1.0d0)*mhd_gamma)-Jfac0**2/36.0d0
          if (rootfac<smalldouble) then
             if(mype==0)then
                print *,'rootfac=',rootfac
                print *,2.0d0*(1.0d0+invbext)/((beta1+1.0d0)*mhd_gamma)
                print *,Jfac0**2/36.0d0
             endif
             call mpistop("inconsistent parameters for TC")
          endif
          Bz0=dsqrt(rootfac)
          rho0=Bz0**2/Rvacs**2
       case('Sakanaka')
          pr01=beta1*(1.0d0+invbext)/((beta1+1.0d0)*mhd_gamma)
          Jfac0=dsqrt(200.0d0*(1.0d0+invbext)/((beta1+1.0d0)*mhd_gamma))/dsqrt(1.0d0+(Lz*qfac1/dpi)**2)
          alfa=-12.0d0+40.0d0/3.0d0-125.0d0/8.0d0+63.0d0/5.0d0+18.0d0/7.0d0-9.0d0/16.0d0+1.0d0/18.0d0
          rootfac=Jfac0**2*(alfa/5.0d0+(Lz*qfac1/(10.0d0*dpi))**2)-2.0d0*pr01*(dexp(4.0d0)-1.0d0)
          if (rootfac<smalldouble) then
             if(mype==0)then
                print *,'rootfac=',rootfac
                print *,'alfa=',alfa
                print *,2.0d0*pr01*(dexp(4.0d0)-1.0d0)
                print *,Jfac0**2*(alfa/5.0d0+(Lz*qfac1/(10.0d0*dpi))**2)
             endif
             call mpistop("inconsistent parameters for Sakanaka")
          endif
          Bz0=dsqrt(rootfac)
          rho0=Bz0**2/Rvacs**2
       case('GoldHoyle')
          pr01=beta1*(1.0d0+invbext)/((beta1+1.0d0)*mhd_gamma)
          Jfac0=dpi/(Lz*qfac1)
          Bz0=dsqrt(2.0d0*(1.0d0+invbext)*(1.0d0+Jfac0**2)/((beta1+1.0d0)*mhd_gamma))
          rho0=Bz0**2/Rvacs**2
       case default
          call mpistop("Unknown equilibrium, choose TokamakCurrent, Sakanaka or GoldHoyle")
    end select

    if(mype==0) then
      write(*,*) "Deduced dimensionless values:"
      write(*,printsettingformat) "pressure at r=1 ",pr01," output"
      write(*,printsettingformat) "density at r=0 ",rho0," output"
      write(*,printsettingformat) "axial B field at r=0 ",Bz0," output"
      write(*,printsettingformat) "current value ",Jfac0," output"
    end if

  end subroutine initglobaldata_usr

  !> Initialize one grid
  subroutine CCC_init_one_grid(ixG^L,ix^L,w,x)
    integer, intent(in) :: ixG^L, ix^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)
    
    {^IFONED   call mpistop("This is a multi-D MHD problem, in 2.5D or 3D") }

    {^NOONED
    w(ix^S,rho_) = rho_solution(x(ix^S, 1))
    w(ix^S,e_)   = p_solution(x(ix^S, 1))
    w(ix^S,mag(1)) = 0.0d0
    w(ix^S,mag(2)) = bphi_solution(x(ix^S, 1))
    w(ix^S,mag(3)) = bz_solution(x(ix^S, 1))
    w(ix^S,mom(1)) = 0.0d0
    w(ix^S,mom(2)) = 0.0d0
    w(ix^S,mom(3)) = 0.0d0
    }

    call mhd_to_conserved(ixG^L,ix^L,w,x)

  end subroutine CCC_init_one_grid

  elemental function p_solution(rad) result(val)
    real(dp), intent(in) :: rad
    real(dp)             :: val

    if(rad<1.0d0)then
       select case(equilibrium_version)
         case('TokamakCurrent')
          val=pr01+Jfac0**2/6.0d0*(47.0d0/120.0d0    &
              -3.0d0*rad**2/2.0d0+9.0d0*rad**4/4.0d0 &
              -5.0d0*rad**6/3.0d0+5.0d0*rad**8/8.0d0-rad**10/10.d0)
         case('Sakanaka')
          val=pr01*dexp(4.0d0-4.0d0*rad**2)
         case('GoldHoyle')
          val=pr01
       end select
    else
      val=1.0d0/mhd_gamma
    endif

  end function p_solution

  elemental function rho_solution(rad) result(val)
    real(dp), intent(in) :: rad
    real(dp)             :: val

    if(rad<1.0d0)then
       val=rho0*(1.0d0-(1.0d0-drat)*rad**2)
    else
      val=1.0d0
    endif

  end function rho_solution

  elemental function bphi_solution(rad) result(val)
    real(dp), intent(in) :: rad
    real(dp)             :: val

    if(rad<1.0d0)then
       select case(equilibrium_version)
         case('TokamakCurrent')
          val=Jfac0*rad*(3.0d0-3.0d0*rad**2+rad**4)/6.0d0
         case('Sakanaka')
          val=Jfac0*(5.0d0*rad-10.d0*rad**3+10.0d0*rad**5-5.0d0*rad**7+rad**9)/10.0d0
         case('GoldHoyle')
          val=Jfac0*rad*Bz0/(1.0d0+(Jfac0*rad)**2)
       end select
    else
      val=0.0d0
    endif

  end function bphi_solution

  elemental function bz_solution(rad) result(val)
    real(dp), intent(in) :: rad
    real(dp)             :: val

    if(rad<1.0d0)then
       select case(equilibrium_version)
         case('TokamakCurrent')
          val=Bz0
         case('Sakanaka')
          val=dsqrt(Bz0**2-0.2d0*Jfac0**2*(5.0d0*rad**2/2.0d0 &
                  -15.0d0*rad**4/2.0d0+40.0d0*rad**6/3.0d0   &
                 -125.0d0*rad**8/8.0d0+63.0d0*rad**10/5.0d0 &
                 -7.0d0*rad**12+18.0d0*rad**14/7.0d0        &
                 -9.0d0*rad**16/16.0d0+rad**18/18.d0) &
               +2.0d0*pr01*dexp(4.0d0)*(1.0d0-dexp(-4.0d0*rad**2)))
         case('GoldHoyle')
          val=Bz0/(1.0d0+(Jfac0*rad)**2)
       end select
    else
      val=dsqrt(invbext*2.0d0/mhd_gamma)
    endif

  end function bz_solution


  subroutine set_error(igrid,level,ixI^L,ixO^L,qt,w,x)
    integer, intent(in)             :: igrid,level,ixI^L,ixO^L
    double precision, intent(in)    :: qt,x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    w(ixO^S,i_sol_r) =1.0d0
    w(ixO^S,i_sol_p) =1.0d0/mhd_gamma
    w(ixO^S,i_sol_b) =dsqrt(invbext*2.0d0/mhd_gamma)
    w(ixO^S,i_totp)  =1.0d0/mhd_gamma
    {^NOONED
    where(x(ixO^S,1)<1.0d0)
      w(ixO^S,i_sol_p) = p_solution(x(ixO^S, 1))
      w(ixO^S,i_sol_r) = rho_solution(x(ixO^S, 1))
      w(ixO^S,i_sol_b) = dsqrt(  bz_solution(x(ixO^S, 1))**2 &
                              +bphi_solution(x(ixO^S, 1))**2)
    endwhere
    }
    w(ixO^S,i_err_r) = w(ixO^S,rho_) - w(ixO^S,i_sol_r)
    call mhd_to_primitive(ixI^L,ixO^L,w,x)
    w(ixO^S,i_err_p) = w(ixO^S,e_) - w(ixO^S,i_sol_p)
    {^NOONED
    w(ixO^S,i_totp) = w(ixO^S,e_) +0.5d0* &
       (w(ixO^S,mag(1))**2+w(ixO^S,mag(2))**2+w(ixO^S,mag(3))**2)
    }
    call mhd_to_conserved(ixI^L,ixO^L,w,x)
    {^NOONED
    w(ixO^S,i_err_b) = dsqrt( w(ixO^S,mag(1))**2+w(ixO^S,mag(2))**2+w(ixO^S,mag(3))**2 ) &
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
    double precision :: Btotal(ixI^S,1:ndir),curlvec(ixI^S,1:ndir)
    integer :: idirmin,idir

    wlocal(ixI^S,1:nw)=w(ixI^S,1:nw)
    ! store current
    call get_current(wlocal,ixI^L,ixO^L,idirmin,curlvec)
    do idir=1,ndir
      w(ixO^S,nw+idir)=curlvec(ixO^S,idir)
    end do

  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
    ! newly added variables need to be concatenated with the w_names/primnames string
    character(len=*) :: varnames

    varnames='j1 j2 j3'

  end subroutine specialvarnames_output

end module mod_usr
