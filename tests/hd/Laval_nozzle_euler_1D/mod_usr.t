! test for 1D Laval nozzle flow
module mod_usr
 use mod_hd
 implicit none
 double precision :: mach0

contains

 subroutine usr_init()
   call usr_params_read(par_files)

   usr_init_one_grid   => initonegrid_usr
   usr_special_bc      => specialbound_usr
   usr_set_surface     => special_surface
   usr_aux_output      => extra_var_output
   usr_add_aux_names   => extra_var_names_output

   call set_coordinate_system("Cartesian_1D_expansion")
   call hd_activate()
 end subroutine usr_init

 subroutine usr_params_read(files)
    use mod_global_parameters, only:unitpar
    use mod_constants
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /usr_list/ mach0
    do n=1, size(files)
      open(unitpar, file=trim(files(n)), status='old')
      read(unitpar, usr_list, end=111)
      111 close(unitpar)
    end do

 end subroutine usr_params_read


  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    ! reference density at inlet
    w(ixO^S,rho_)=1.0d0
    ! velocity as inlet mach value
    w(ixO^S,mom(1))=mach0
    if(hd_energy)then
       ! pressure to make inlet sound speed unity
       w(ixO^S,e_)=1.0d0/hd_gamma
    endif

    call hd_to_conserved(ixI^L, ixO^L, w, x)

  end subroutine initonegrid_usr


  subroutine specialbound_usr(qt,ixI^L,ixO^L,iB,w,x)
    integer, intent(in) :: ixO^L, iB, ixI^L
    double precision, intent(in) :: qt,x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    integer:: ixOInt^L,ix1

    select case(iB)
    case(1)
      ! fix the left inlet boundary
      w(ixO^S,rho_)= 1.0d0
      w(ixO^S,mom(1))= mach0
      if(hd_energy)then
         w(ixO^S,e_)=1.0d0/hd_gamma
      endif
      call hd_to_conserved(ixI^L,ixO^L,w,x)
    case(2)
      ixOInt^L=ixO^L;
      ixOIntmin1=ixOmin1-1 
      ixOIntmax1=ixOmin1-1 
      call hd_to_primitive(ixI^L,ixOInt^L,w,x)
      if(hd_energy)then
      ! extrapolate density, velocity, pressure
         do ix1=ixOmin1,ixOmax1
            w(ix1,rho_)=w(ixOmin1-1,rho_)
            w(ix1,mom(1))=w(ixOmin1-1,mom(1))
            w(ix1,e_)=w(ixOmin1-1,e_)
         enddo
      else
      ! extrapolate velocity, density 
         do ix1=ixOmin1,ixOmax1
            w(ix1,rho_)=w(ixOmin1-1,rho_)
            w(ix1,mom(1))=w(ixOmin1-1,mom(1))
         enddo
      endif
      call hd_to_conserved(ixI^L,ixOInt^L,w,x)
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

    where(dabs(x(ixI^S,1))>=2.0d0)
       exp_factor(ixI^S)=5.0d0
       del_exp_factor(ixI^S)=0.0d0
       exp_factor_primitive(ixI^S)=5.0d0*dx(ixI^S,1)
    elsewhere
       exp_factor(ixI^S)=x(ixI^S,1)**2+1.0d0
       del_exp_factor(ixI^S)=2.0d0*x(ixI^S,1)
       exp_factor_primitive(ixI^S)=dx(ixI^S,1)*(x(ixI^S,1)**2+1.0d0+dx(ixI^S,1)**2/12.0d0)
    endwhere

  end subroutine special_surface

  subroutine extra_var_output(ixI^L, ixO^L, w, x, normconv)
    use mod_physics
    integer, intent(in)              :: ixI^L,ixO^L
    double precision, intent(in)     :: x(ixI^S,1)
    double precision                 :: w(ixI^S,nw+nwauxio), wlocal(ixI^S,1:nw)
    double precision                 :: normconv(0:nw+nwauxio)
    double precision                 :: pp(ixI^S),vv(ixI^S),cs2(ixI^S)
    double precision                 :: exp_factor(ixI^S),del_exp_factor(ixI^S),tmp(ixI^S)
    double precision                 :: exp_factor_primitive(ixI^S),dx(ixI^S,1)

    ! store expansion factor, but note that dx is unknown here
    dx(ixI^S,1)=0.0d0
    call special_surface(ixI^L,x,dx,exp_factor,del_exp_factor,exp_factor_primitive)
    w(ixO^S, nw+1) = exp_factor(ixO^S)
    
    call hd_get_pthermal(w,x,ixI^L,ixO^L,pp)
    w(ixO^S,nw+2) = pp(ixO^S)/w(ixO^S,rho_)
    vv(ixO^S)=w(ixO^S,mom(1))/w(ixO^S,rho_)
    call hd_get_csound2(w,x,ixI^L,ixO^L,cs2)
    w(ixO^S,nw+3) =vv(ixO^S)/dsqrt(cs2(ixO^S))
    w(ixO^S,nw+4) =vv(ixO^S)*exp_factor(ixO^S)*w(ixO^S,rho_)
    if(hd_energy.or.hd_gamma/=1.0d0)then
       w(ixO^S,nw+5)=0.5d0*vv(ixO^S)**2+(pp(ixO^S)/w(ixO^S,rho_))*(hd_gamma/(hd_gamma-1.0d0))
    else
      print *,'using lnrho: energy,gamma,adiab'
      print *,hd_energy,hd_gamma,hd_adiab
       w(ixO^S,nw+5)=0.5d0*vv(ixO^S)**2+hd_adiab*dlog(w(ixO^S,rho_))
    endif

  end subroutine extra_var_output

  subroutine extra_var_names_output(varnames)
    character(len=*)  :: varnames

    varnames = "A T M Am E"

  end subroutine extra_var_names_output

end module mod_usr
