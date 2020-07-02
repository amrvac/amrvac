module mod_usr
 use mod_hd
 implicit none
 double precision :: beta, p0, T0, v0

contains

 subroutine usr_init()

   unit_length        = 1.d9
   unit_temperature   = 1.d6
   unit_numberdensity = 1.d9

   !usr_set_parameters  => initglobaldata_usr
   usr_init_one_grid   => initonegrid_usr
   usr_special_bc      => specialbound_usr
   usr_set_surface     => special_surface
   usr_aux_output      => extra_var_output
   usr_add_aux_names   => extra_var_names_output
   call usr_params_read(par_files)
   call set_coordinate_system("Cartesian_1D_expansion")
   call hd_activate()
 end subroutine usr_init

 subroutine usr_params_read(files)
    use mod_global_parameters, only:unitpar
    use mod_constants
    character(len=*), intent(in) :: files(:)
    integer                      :: n
    namelist /InitalConst_list/ p0, v0, T0, beta
    do n=1, size(files)
      open(unitpar, file=trim(files(n)), status='old')
      read(unitpar, InitalConst_list, end=111)
      111 close(unitpar)
    end do
 end subroutine usr_params_read


  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    w(ixO^S,rho_)=p0   !*beta/T0
    w(ixO^S,mom(1))=v0

    call hd_to_conserved(ixI^L, ixO^L, w, x)

  end subroutine initonegrid_usr


  subroutine specialbound_usr(qt,ixI^L,ixO^L,iB,w,x)
    integer, intent(in) :: ixO^L, iB, ixI^L
    double precision, intent(in) :: qt,x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    !Set here the left boundary, velocity source
    select case(iB)
    case(1)
      w(ixO^S,rho_)= p0    !/T0
      w(ixOmin1,mom(1))= v0/w(ixOmin1,rho_)
      w(ixOmax1,mom(1))= v0/w(ixOmax1,rho_)
      call hd_to_conserved(ixI^L,ixO^L,w,x)

      case default
       call mpistop("Special boundary is not defined for this region")
     end select

  end subroutine specialbound_usr

  subroutine special_surface(ixI^L,x,dx,exp_factor,del_exp_factor,exp_factor_primitive)
    integer, intent(in)             :: ixI^L
    double precision, intent(in)    :: dx(ixI^S,1), x(ixI^S,1)
    double precision, intent(out)   :: exp_factor(ixI^S), del_exp_factor(ixI^S)
    double precision, intent(out)   :: exp_factor_primitive(ixI^S)

    !> quadratic expansion factor
    !where (x(ixI^S,1) .ge. 4.50 .and. x(ixI^S,1) .le. 5.5d0)
      exp_factor(ixI^S) = (x(ixI^S,1)-5.0d0)**2 + 0.5d0
      del_exp_factor(ixI^S) = 2.0d0*(x(ixI^S,1)-5.0d0)
      exp_factor_primitive(ixI^S) = (x(ixI^S,1)**2+dx(ixI^S,1)**2/12.0d0-10.0d0*x(ixI^S,1)+26.0d0)*dx(ixI^S,1)
    !elsewhere
      !exp_factor(ixI^S) = 0.75d0
      !del_exp_factor(ixI^S) = 0.0d0
      !exp_factor_primitive(ixI^S) = 0.75d0*dx(ixI^S,1)
    !endwhere

    !> function from Mikic et al. (2013)
    !exp_factor(ixI^S) = 1.0d0/(1.0d0 + 10.d0*(exp(-x(ixI^S,1)/14.d0) + exp(-(xprobmax1-x(ixI^S,1))/14.d0)))
    !del_exp_factor(ixI^S) = -0.7143*(exp(-(xprobmax1-x(ixI^S,1))/14.d0)-exp(-x(ixI^S,1)/14.d0))/exp_factor(ixI^S)**2
    !exp_factor_primitive(ixI^S) = exp_factor(ixI^S)*dx(ixI^S,1)
  end subroutine special_surface


  subroutine extra_var_output(ixI^L, ixO^L, w, x, normconv)
    use mod_physics
    integer, intent(in)              :: ixI^L,ixO^L
    double precision, intent(in)     :: x(ixI^S,1)
    double precision                 :: w(ixI^S,nw+nwauxio), wlocal(ixI^S,1:nw)
    double precision                 :: normconv(0:nw+nwauxio)
    double precision                 :: exp_factor(ixI^S),del_exp_factor(ixI^S)
    double precision                 :: exp_factor_primitive(ixI^S),dx(ixI^S,1)

    ! store expansion factor
    call special_surface(ixI^L,x,dx,exp_factor,del_exp_factor,exp_factor_primitive)
    w(ixO^S, nw+1) = exp_factor
  end subroutine extra_var_output

  subroutine extra_var_names_output(varnames)
    character(len=*)  :: varnames

    varnames = "Ax"

  end subroutine extra_var_names_output

end module mod_usr
