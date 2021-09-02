!> Advection of a current-carrying cylinder
! implemented a choice between three 1D cylindrical equilibria,
! embedded in a uniform, (un)magnetized medium, advected at constant speed
! This tests Galilean invariance, among other things.
module mod_usr
  use mod_mhd

  implicit none

  ! Input values for
  !    uniform advection velocity: Mach, theta0, phi0
  !    choice of equilibrium: equilibrium_version
  !    dimensionless parameters: beta1, qfac1, Rvacs, drat, invbext
  double precision :: Mach, theta0, phi0, beta1, qfac1, Rvacs, drat, invbext
  double precision :: pr01, Bz0, rho0, Jfac0, Lz
  character(len=std_len) :: equilibrium_version

  ! The initial center of the cylinder
  double precision :: xmid^D

  ! Velocity
  double precision :: v01, v02

  integer, parameter :: dp = kind(0.0d0)
  integer :: i_err_p, i_err_r, i_err_b
  integer :: i_sol_p, i_sol_r, i_sol_b, i_totp
  integer :: i_divb

  double precision, parameter :: mask_halfwidth = 0.1d0

contains

  !> Read this module's parameters from a file
  subroutine usr_params_read(files)
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /usr_list/ Mach, theta0, phi0, beta1, qfac1,  &
         Rvacs, drat, invbext, equilibrium_version

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, usr_list, end=111)
111    close(unitpar)
    end do

    ! convert degrees to radians
    theta0=theta0*2.0d0*dpi/360.0d0
    phi0=phi0*2.0d0*dpi/360.0d0
    equilibrium_version=trim(equilibrium_version)

  end subroutine usr_params_read

  subroutine usr_init()
    use mod_variables

    call usr_params_read(par_files)

    usr_set_parameters  => initglobaldata_usr
    usr_init_one_grid => CCC_init_one_grid
    usr_print_log => print_error
    usr_write_analysis => analyze_forces_on_grid
    usr_refine_grid    => specialrefine_grid
    usr_aux_output      => specialvar_output
    usr_add_aux_names   => specialvarnames_output

    {^IFTWOD call set_coordinate_system("Cartesian_2.5D") }
    {^IFTHREED call set_coordinate_system("Cartesian_3D")  }

    call mhd_activate()

    i_err_r = var_set_extravar("rho_error", "rho_error")
    i_sol_r = var_set_extravar("rho_solution", "rho_solution")
    i_err_p = var_set_extravar("p_error", "p_error")
    i_sol_p = var_set_extravar("p_solution", "p_solution")
    i_err_b = var_set_extravar("b_error", "b_error")
    i_sol_b = var_set_extravar("b_solution", "b_solution")
    i_totp  = var_set_extravar("totp", "totp")
    i_divb  = var_set_extravar("divb", "divb")

  end subroutine usr_init

  subroutine initglobaldata_usr
    character(len=20) :: printsettingformat
    double precision :: rootfac, alfa

    ^D&xmid^D=xprobmin^D+0.5d0*(xprobmax^D-xprobmin^D);

    v01= Mach*dsin(theta0)*dcos(phi0)
    v02= Mach*dsin(theta0)*dsin(phi0)

    printsettingformat='(1x,A50,ES15.7,A7)'

    {^IFTHREED
    Lz=(xprobmax3-xprobmin3)
    }
    {^IFTWOD
    Lz=(xprobmax1-xprobmin1)
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
          Jfac0=6.0d0*sqrt(2.0d0*(1.0d0+invbext)/((beta1+1.0d0)*mhd_gamma))/sqrt(1.0d0+(Lz*qfac1/dpi)**2)
          rootfac=2.0d0*(1.0d0+invbext)/((beta1+1.0d0)*mhd_gamma)-Jfac0**2/36.0d0
          if (rootfac<smalldouble) then
             if(mype==0)then
                print *,'rootfac=',rootfac
                print *,2.0d0*(1.0d0+invbext)/((beta1+1.0d0)*mhd_gamma)
                print *,Jfac0**2/36.0d0
             endif
             call mpistop("inconsistent parameters for TC")
          endif
          Bz0=sqrt(rootfac)
          rho0=Bz0**2/Rvacs**2
       case('Sakanaka')
          pr01=beta1*(1.0d0+invbext)/((beta1+1.0d0)*mhd_gamma)
          Jfac0=sqrt(200.0d0*(1.0d0+invbext)/((beta1+1.0d0)*mhd_gamma))/sqrt(1.0d0+(Lz*qfac1/dpi)**2)
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
          Bz0=sqrt(rootfac)
          rho0=Bz0**2/Rvacs**2
       case('GoldHoyle')
          pr01=beta1*(1.0d0+invbext)/((beta1+1.0d0)*mhd_gamma)
          Jfac0=dpi/(Lz*qfac1)
          Bz0=sqrt(2.0d0*(1.0d0+invbext)*(1.0d0+Jfac0**2)/((beta1+1.0d0)*mhd_gamma))
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
      write(*,printsettingformat) "v_x ", v01," output"
      write(*,printsettingformat) "v_y ", v02," output"
    end if

  end subroutine initglobaldata_usr

  !> Initialize one grid
  subroutine CCC_init_one_grid(ixG^L,ix^L,w,x)
    integer, intent(in) :: ixG^L, ix^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    double precision :: rr(ixG^S),bphi(ixG^S),cosphi(ixG^S),sinphi(ixG^S)

    {^IFONED   call mpistop("This is a multi-D MHD problem, in 2.5D or 3D") }

    {^NOONED
    rr(ix^S)=sqrt( (x(ix^S,1)-xmid1)**2+(x(ix^S,2)-xmid2)**2 )
    w(ix^S,rho_) = rho_solution(x(ix^S, 1)-xmid1, x(ix^S, 2)-xmid2)
    w(ix^S,e_)   = p_solution(x(ix^S, 1)-xmid1, x(ix^S, 2)-xmid2)
    bphi(ix^S)   = bphi_solution(x(ix^S, 1)-xmid1, x(ix^S, 2)-xmid2)
    sinphi(ix^S)=(x(ix^S,2)-xmid2)/rr(ix^S)
    cosphi(ix^S)=(x(ix^S,1)-xmid1)/rr(ix^S)
    w(ix^S,mag(1)) = -bphi(ix^S)*sinphi(ix^S)
    w(ix^S,mag(2)) = +bphi(ix^S)*cosphi(ix^S)
    w(ix^S,mag(3)) = bz_solution(x(ix^S, 1)-xmid1, x(ix^S, 2)-xmid2)
    w(ix^S,mom(1)) = Mach*dsin(theta0)*dcos(phi0)
    w(ix^S,mom(2)) = Mach*dsin(theta0)*dsin(phi0)
    w(ix^S,mom(3)) = Mach*dcos(theta0)
    w(ix^S,tracer(1))=0.0d0
    where(rr(ix^S)<0.99d0)
      w(ix^S,tracer(1))=1.0d0
    endwhere
    where((rr(ix^S)>0.99d0).and.(rr(ix^S)<1.01d0))
      w(ix^S,tracer(1))=0.5d0
    endwhere
    }

    call mhd_to_conserved(ixG^L,ix^L,w,x)

  end subroutine CCC_init_one_grid

  elemental function p_solution(x, y) result(val)
    real(dp), intent(in) :: x, y
    real(dp)             :: val
    real(dp)             :: rad

    rad=sqrt(x**2+y**2)
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

  elemental function rho_solution(x, y) result(val)
    real(dp), intent(in) :: x, y
    real(dp)             :: val
    real(dp)             :: rad

    rad=sqrt(x**2+y**2)
    if(rad<1.0d0)then
       val=rho0*(1.0d0-(1.0d0-drat)*rad**2)
    else
      val=1.0d0
    endif

  end function rho_solution

  elemental function bphi_solution(x, y) result(val)
    real(dp), intent(in) :: x, y
    real(dp)             :: val
    real(dp)             :: rad

    rad=sqrt(x**2+y**2)
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

  elemental function bx_solution(x, y) result(val)
    real(dp), intent(in) :: x, y
    real(dp)             :: val
    real(dp)             :: rad, sinphi

    rad=max(epsilon(1.0d0), sqrt(x**2+y**2))
    sinphi=y/rad

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
    val = -val * sinphi
  end function bx_solution

  elemental function by_solution(x, y) result(val)
    real(dp), intent(in) :: x, y
    real(dp)             :: val
    real(dp)             :: rad, cosphi

    rad=max(epsilon(1.0d0), sqrt(x**2+y**2))
    cosphi=x/rad

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
    val = val * cosphi
  end function by_solution

  elemental function bz_solution(x, y) result(val)
    real(dp), intent(in) :: x, y
    real(dp)             :: val
    real(dp) :: rad

    rad=sqrt(x**2+y**2)
    if(rad<1.0d0)then
       select case(equilibrium_version)
         case('TokamakCurrent')
          val=Bz0
         case('Sakanaka')
          val=sqrt(Bz0**2-0.2d0*Jfac0**2*(5.0d0*rad**2/2.0d0 &
                  -15.0d0*rad**4/2.0d0+40.0d0*rad**6/3.0d0   &
                 -125.0d0*rad**8/8.0d0+63.0d0*rad**10/5.0d0 &
                 -7.0d0*rad**12+18.0d0*rad**14/7.0d0        &
                 -9.0d0*rad**16/16.0d0+rad**18/18.d0) &
               +2.0d0*pr01*dexp(4.0d0)*(1.0d0-dexp(-4.0d0*rad**2)))
         case('GoldHoyle')
          val=Bz0/(1.0d0+(Jfac0*rad)**2)
       end select
    else
      val=sqrt(invbext*2.0d0/mhd_gamma)
    endif

  end function bz_solution

  subroutine set_error(ixI^L,ixO^L,qt,w,x)
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(in)    :: qt,x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    integer :: ic1,ic2
    double precision :: rr(ixI^S)
    double precision :: x1(ixO^S), x2(ixO^S)

    w(ixO^S,i_sol_r) =1.0d0
    w(ixO^S,i_sol_p) =1.0d0/mhd_gamma
    w(ixO^S,i_sol_b) =sqrt(invbext*2.0d0/mhd_gamma)
    w(ixO^S,i_totp)  =1.0d0/mhd_gamma

    {^NOONED
    x1 = x(ixO^S,1) - v01 * qt - xprobmin1
    x1 = xprobmin1 + modulo(x1, xprobmax1-xprobmin1) - xmid1
    x2 = x(ixO^S,2) - v02 * qt - xprobmin2
    x2 = xprobmin2 + modulo(x2, xprobmax2-xprobmin2) - xmid2

    rr(ixO^S) = sqrt(x1**2 + x2**2)
    where(rr(ixO^S)<1.0d0)
       w(ixO^S,i_sol_p) = p_solution(x1, x2)
       w(ixO^S,i_sol_r) = rho_solution(x1, x2)
       w(ixO^S,i_sol_b) = sqrt(bz_solution(x1, x2)**2 &
            + bphi_solution(x1, x2)**2)
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
    w(ixO^S,i_err_b) = sqrt( &
         (w(ixO^S,mag(1)) - bx_solution(x1, x2))**2 + &
         (w(ixO^S,mag(2)) - by_solution(x1, x2))**2 + &
         (w(ixO^S,mag(3)) - bz_solution(x1, x2))**2)
    }
  end subroutine set_error

  subroutine print_error()
    use mod_input_output, only: get_volume_average
    double precision   :: modes(nw, 2), volume
    double precision :: divb(ixG^T),sumdivb_mype,sumdivb
    double precision :: current(ixG^T,7-2*ndir:3),Jval(ixG^T)
    double precision :: jmax_mype,jmax,jmin_mype,jmin
    character(len=100):: filename
    character(len=1024) :: line, datastr
    integer :: iigrid, igrid, idirmin
    logical :: first_time = .true.

    ! get divb (not normalized since also unmagnetized regions present),
    ! first sum over all grids per processor
    sumdivb_mype=0.0d0
    ! Get maximal and mimimal current value
    jmax_mype = -bigdouble
    jmin_mype = bigdouble
    ! Loop over all the grids local to processor
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       block=>ps(igrid)
       ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
       call set_error(ixG^LL,ixM^LL,global_time,ps(igrid)%w, ps(igrid)%x)
       call get_divb(ps(igrid)%w,ixG^LL,ixM^LL,divb)
       sumdivb_mype=sumdivb_mype+sum(dabs(divb(ixM^T)))
       call get_current(ps(igrid)%w,ixG^LL,ixM^LL,idirmin,current)
       Jval(ixM^T)=sqrt( current(ixM^T,1)**2+current(ixM^T,2)**2+current(ixM^T,3)**2 )
       !jmax_mype=max(jmax_mype,maxval(Jval(ixM^T)))
       jmax_mype=max(jmax_mype,maxval(current(ixM^T,3)))
       jmin_mype=min(jmin_mype,minval(Jval(ixM^T)))
    end do

    call get_volume_average(1, modes(:, 1), volume)
    call get_volume_average(2, modes(:, 2), volume)

    ! now add up all contributions from all processors
    ! only reduce to mype=0 processor, used for output
    call MPI_REDUCE(sumdivb_mype,sumdivb,1,MPI_DOUBLE_PRECISION, &
                 MPI_SUM,0,icomm,ierrmpi)
    call MPI_REDUCE(jmax_mype, jmax, 1, MPI_DOUBLE_PRECISION, &
         MPI_MAX, 0, icomm, ierrmpi)
    call MPI_REDUCE(jmin_mype, jmin, 1, MPI_DOUBLE_PRECISION, &
         MPI_MIN, 0, icomm, ierrmpi)
    ! divide by total domain volume
    {^IFTWOD
    sumdivb=sumdivb/((xprobmax1-xprobmin1)*(xprobmax2-xprobmin2)) }
    {^IFTHREED
    sumdivb=sumdivb/((xprobmax1-xprobmin1)*(xprobmax2-xprobmin2)*(xprobmax3-xprobmin3)) }

    if(mype==0) then
      filename = trim(base_filename) // '_errors.csv'

      if(first_time) then
        open(unit=21,file=filename,form='formatted')
        write(21,'(a)') 'time, rho-error, p-error, b-error, divbsum, jmax, jmin'
        first_time = .false.
      else
        open(unit=21,file=filename,form='formatted',status='old',access='append')
      endif
      write(datastr,'(es11.4, 2a)') global_time,', '
      line=datastr
      write(datastr,"(es12.5, 2a)") sqrt(modes(i_err_r, 2)),', '
      line = trim(line)//trim(datastr)
      write(datastr,"(es12.5, 2a)") sqrt(modes(i_err_p, 2)),', '
      line = trim(line)//trim(datastr)
      write(datastr,"(es12.5, 2a)") sqrt(modes(i_err_b, 2)),', '
      line = trim(line)//trim(datastr)
      write(datastr,"(es12.5, 2a)") sumdivb,', '
      line = trim(line)//trim(datastr)
      write(datastr,"(es12.5, 2a)") jmax,', '
      line = trim(line)//trim(datastr)
      write(datastr,"(es12.5)") jmin
      line = trim(line)//trim(datastr)
      write(21,'(a)') trim(line)
      close(21)
    endif

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
    double precision:: rr(ixG^S)

    {^NOONED
    rr(ix^S)=sqrt( (x(ix^S,1)-xmid1)**2+(x(ix^S,2)-xmid2)**2 )
    }
    ! test with different levels of refinement enforced
    if (qt<smalldouble.and.any((rr(ix^S)<1.1d0).and.(rr(ix^S)>0.9d0))) then
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

    double precision :: qt
    double precision :: wlocal(ixI^S,1:nw)
    double precision :: curlvec(ixI^S,1:ndir)
    double precision :: divb(ixI^S)
    integer :: idirmin,idir
    integer :: ic1,ic2
    double precision:: rr(ixI^S)
    double precision :: x1(ixO^S), x2(ixO^S)

    wlocal(ixI^S,1:nw)=w(ixI^S,1:nw)
    ! store current
    call get_current(wlocal,ixI^L,ixO^L,idirmin,curlvec)
    do idir=1,ndir
      w(ixO^S,nw+idir)=curlvec(ixO^S,idir)
    end do
    ! test the mask: should be following the fluxtube
    w(ixO^S,nw+4)=0.0d0
    qt=global_time

    {^NOONED
    x1 = x(ixO^S,1) - v01 * qt - xprobmin1
    x1 = xprobmin1 + modulo(x1, xprobmax1-xprobmin1) - xmid1
    x2 = x(ixO^S,2) - v02 * qt - xprobmin2
    x2 = xprobmin2 + modulo(x2, xprobmax2-xprobmin2) - xmid2
    rr(ixO^S)=sqrt(x1**2 + x2**2)
    where(rr(ixO^S)< 1.0d0 - mask_halfwidth)
       w(ixO^S,nw+4)=1.0d0
    endwhere
    where(rr(ixO^S) > 1.0d0 - mask_halfwidth .and. &
         rr(ixO^S) < 1.0d0 + mask_halfwidth)
       w(ixO^S,nw+4)=0.5d0
    endwhere
    }

    ! output divB
    call get_divb(w,ixI^L,ixO^L,divb)
    w(ixO^S,i_divb)=divb(ixO^S)

    call set_error(ixI^L,ixO^L,global_time,w,x)

  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
    ! newly added variables need to be concatenated with the w_names/primnames string
    character(len=*) :: varnames

    varnames='j1 j2 j3 mask'

  end subroutine specialvarnames_output

  subroutine analyze_forces_on_grid()
  ! all analysis is handled in compute_integrated_quantities
  ! here we ensure that the data is in primitive variables
    integer :: iigrid, igrid

    ! first switch all grids (with ghost cells) to primitive variables
    do iigrid=1,igridstail; igrid=igrids(iigrid)
       call mhd_to_primitive(ixG^LL,ixG^LL,ps(igrid)%w,ps(igrid)%x)
    end do

    call compute_integrated_quantities

    ! end with switching all back to conservative variables
    do iigrid=1,igridstail; igrid=igrids(iigrid)
       call mhd_to_conserved(ixG^LL,ixG^LL,ps(igrid)%w,ps(igrid)%x)
    end do

  end subroutine analyze_forces_on_grid

  subroutine compute_integrated_quantities
  ! here we do analysis on the whole grid tree:
  ! output is written to the ***analysis.int file
  integer :: nregions, nintegrals, ireg, intval, iigrid, igrid, idim, idirmin
  character(len=std_len), allocatable :: region_name(:)
  character(len=std_len) :: fstring
  double precision, allocatable :: integral_ipe(:,:),integral_w(:,:)
  double precision :: pgrad(ixG^T,1:3), tmp(ixG^T), dvolume(ixG^T), curlvec(ixG^T,1:ndir)
  logical :: first_time = .true.
  character(len=100):: filename


   ! identify region within, at edge, or outside fluxtube, or full domain
   nregions=4
   allocate(region_name(nregions))
   region_name(1)='inside'
   region_name(2)='edge'
   region_name(3)='outside'
   region_name(4)='domain'
   ! number of integrals to compute
   nintegrals=4
   ! format for 4 integrals to report, with time (so 5 DP)
   fstring='(5(es12.4))'
   allocate(integral_ipe(nregions,nintegrals),integral_w(nregions,nintegrals))
   integral_ipe=0.d0
   integral_w=0.d0

   ! Loop over all the grids
   do iigrid = 1, igridstail
      igrid = igrids(iigrid)
      block=>ps(igrid)
      ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
      ! Determine the volume of the grid cells
      if (slab) then
         dvolume(ixM^T) = {rnode(rpdx^D_,igrid)|*}
      else
         dvolume(ixM^T) = ps(igrid)%dvolume(ixM^T)
      end if

      call get_current(ps(igrid)%w,ixG^LL,ixM^LL,idirmin,curlvec)

      ! pressure gradient
      do idim=1,ndim
        select case(typegrad)
        case("central")
          call gradient(ps(igrid)%w(ixG^T,p_),ixG^LL,ixM^LL,idim,tmp)
        case("limited")
          call gradientS(ps(igrid)%w(ixG^T,p_),ixG^LL,ixM^LL,idim,tmp)
        end select
        pgrad(ixM^T,idim) = tmp(ixM^T)
      end do
      if(ndim==2)pgrad(ixM^T,3) =0.0d0

      do intval=1,nintegrals
         do ireg=1,nregions
              integral_ipe(ireg,intval)=integral_ipe(ireg,intval)+ &
                integral_grid(ixG^LL,ixM^LL,ps(igrid)%w,curlvec,pgrad,ps(igrid)%x,dvolume,ireg,intval)
         enddo
      enddo
   enddo
   call MPI_ALLREDUCE(integral_ipe,integral_w,nregions*nintegrals,MPI_DOUBLE_PRECISION, &
       MPI_SUM,icomm,ierrmpi)

   if(mype==0) then
       do ireg=1,nregions
          write(filename,"(a,a,a)") TRIM(base_filename),'_'//TRIM(region_name(ireg)), &
                                  "_analysis.int"
          if(first_time) then
            open(unit=21,file=filename,form='formatted')
            write(unit=21, fmt=*) "time volume lorentz lorentz-pressure inv-beta"
          else
            open(unit=21,file=filename,form='formatted',status='old',access='append')
          end if

          write(21,fstring) global_time,integral_w(ireg,1:nintegrals)
          close(21)
       enddo
       first_time = .false.
   endif

   deallocate(integral_ipe,integral_w)
   deallocate(region_name)

 end subroutine compute_integrated_quantities

 function integral_grid(ixI^L,ixO^L,w,Jvec,pgrad,x,dvolume,ireg,intval)

  ! note: on entry, we have primitives in w vector,
  !     and current in Jvec, and pressure gradient in pgrad
  !  we will quantify
  !  1)total volume of region
  !  2)magnitude of JxB force in regions identified by mask
  !  3)magnitude of (JxB-grad p) force in regions identified by mask
  !  4)inverse plasma beta in regions identified by mask

  integer, intent(in)                :: ixI^L,ixO^L,ireg,intval
  double precision, intent(in)       :: x(ixI^S,1:ndim),dvolume(ixI^S)
  double precision, intent(in)       :: w(ixI^S,1:nw),Jvec(ixI^S,1:ndir),pgrad(ixI^S,1:3)

  logical :: patchw(ixI^S)
  double precision :: mask(ixI^S)
  double precision :: qt, rr(ixI^S)
  double precision :: integral_grid,Bf(ixI^S,1:ndir),JcrosB(ixI^S,1:3),totforce(ixI^S)
  integer :: ix^D,ic1,ic2
  double precision :: x1(ixO^S), x2(ixO^S)

  mask(ixO^S)=0.0d0
  qt=global_time

  {^NOONED
  x1 = x(ixO^S,1) - v01 * qt - xprobmin1
  x1 = xprobmin1 + modulo(x1, xprobmax1-xprobmin1) - xmid1
  x2 = x(ixO^S,2) - v02 * qt - xprobmin2
  x2 = xprobmin2 + modulo(x2, xprobmax2-xprobmin2) - xmid2
  rr(ixO^S)=sqrt(x1**2 + x2**2)
  where(rr(ixO^S) < 1.0d0 - mask_halfwidth)
     mask(ixO^S)=1.0d0
  endwhere
  where(rr(ixO^S) > 1.0d0 - mask_halfwidth .and. &
       rr(ixO^S) < 1.0d0 + mask_halfwidth)
     mask(ixO^S)=0.5d0
  endwhere
  }
  ! mask the region of interest
  select case(ireg)
    case (1) ! inside
      patchw(ixO^S)=(mask(ixO^S)>0.95d0)
    case (2) ! edge
      patchw(ixO^S)=((mask(ixO^S)<=0.95d0).and.(mask(ixO^S)>=0.05d0))
    case (3) ! outside
      patchw(ixO^S)=(mask(ixO^S)<0.05d0)
    case (4) ! full domain
      patchw(ixO^S)=.true.
  end select

  integral_grid=0.d0
  select case(intval)
     case(1)
       ! volume integration
       {do ix^DB=ixOmin^DB,ixOmax^DB\}
          if(patchw(ix^D)) integral_grid=integral_grid+dvolume(ix^D)
       {end do\}
     case(2)
       ! Lorentz force
       ! magnetic field
       Bf(ixI^S,1:ndir)=w(ixI^S,mag(1:ndir))

       JcrosB(ixO^S,1)=Jvec(ixO^S,2)*Bf(ixO^S,3)-Jvec(ixO^S,3)*Bf(ixO^S,2)
       JcrosB(ixO^S,2)=Jvec(ixO^S,3)*Bf(ixO^S,1)-Jvec(ixO^S,1)*Bf(ixO^S,3)
       JcrosB(ixO^S,3)=Jvec(ixO^S,1)*Bf(ixO^S,2)-Jvec(ixO^S,2)*Bf(ixO^S,1)

       totforce(ixO^S)=sqrt( JcrosB(ixO^S,1)**2+JcrosB(ixO^S,2)**2+JcrosB(ixO^S,3)**2 )
       {do ix^DB=ixOmin^DB,ixOmax^DB\}
         if(patchw(ix^D)) integral_grid=integral_grid+totforce(ix^D)*dvolume(ix^D)
       {end do\}
     case(3)
       ! Lorentz force and pressure gradient
       ! magnetic field
       Bf(ixO^S,1:ndir)=w(ixO^S,mag(1:ndir))

       JcrosB(ixO^S,1)=Jvec(ixO^S,2)*Bf(ixO^S,3)-Jvec(ixO^S,3)*Bf(ixO^S,2)
       JcrosB(ixO^S,2)=Jvec(ixO^S,3)*Bf(ixO^S,1)-Jvec(ixO^S,1)*Bf(ixO^S,3)
       JcrosB(ixO^S,3)=Jvec(ixO^S,1)*Bf(ixO^S,2)-Jvec(ixO^S,2)*Bf(ixO^S,1)

       totforce(ixO^S)=sqrt(  (JcrosB(ixO^S,1)-pgrad(ixO^S,1))**2  &
                              +(JcrosB(ixO^S,2)-pgrad(ixO^S,2))**2  &
                              +(JcrosB(ixO^S,3)-pgrad(ixO^S,3))**2 )
       {do ix^DB=ixOmin^DB,ixOmax^DB\}
         if(patchw(ix^D)) integral_grid=integral_grid+totforce(ix^D)*dvolume(ix^D)
       {end do\}
     case(4)
       ! inverse plasma beta integration
       Bf(ixO^S,1:ndir)=w(ixO^S,mag(1:ndir))
       {do ix^DB=ixOmin^DB,ixOmax^DB\}
          if(patchw(ix^D)) integral_grid=integral_grid+ &
                (Bf(ix^D,1)**2+Bf(ix^D,2)**2+Bf(ix^D,3)**2)/(two*w(ix^D,p_))*dvolume(ix^D)
       {end do\}
      case default
        call mpistop("intval not defined")
    end select
    return
  end function integral_grid

end module mod_usr
