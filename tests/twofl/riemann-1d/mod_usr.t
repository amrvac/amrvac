module mod_usr
  use mod_twofl

  implicit none

contains



  !> Read this module's parameters from a file
  subroutine usr_params_read(files)
    character(len=*), intent(in) :: files(:)
  end subroutine usr_params_read

  subroutine usr_init()
    call usr_params_read(par_files)
    usr_init_one_grid => initonegrid_usr
    usr_var_for_errest => special_errest
    usr_special_convert => usrspecial_convert
    call set_coordinate_system("Cartesian_1.5D")
    call twofl_activate()
  end subroutine usr_init
  
  subroutine usrspecial_convert(qunitconvert)
    integer, intent(in) :: qunitconvert
    character(len=20):: userconvert_type
  end subroutine usrspecial_convert

  !refine taking into account only vertical velocity decoupling
  subroutine special_errest(ixI^L,ixO^L,iflag,w,x,var)
    integer, intent(in)           :: ixI^L,ixO^L,iflag
    double precision, intent(in)  :: w(ixI^S,1:nw),x(ixI^S,1:ndim)
    double precision, intent(out) :: var(ixI^S)

    var(ixO^S) = abs(w(ixO^S,mom_c(1))/w(ixO^S,rho_c_)-w(ixO^S,mom_n(1))/w(ixO^S,rho_n_))

  end subroutine special_errest

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    ! initialize one grid
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision,parameter :: rhotot0 = 1d0
    double precision,parameter :: bx0 = 1d-1
    double precision,parameter :: by0 = -1d0
    double precision,parameter :: H_ion_fr  =  1d-1
    double precision,parameter :: He_abundance  =  0d0
    double precision,parameter :: He_ion_fr  =  0d0
    double precision,parameter :: He_ion_fr2  =  0d0


    double precision :: petot0 
    double precision :: B0
    double precision :: beta 
    double precision :: pe_nc_fr, rho_nc_fr

    rho_nc_fr = (1d0 - H_ion_fr)/H_ion_fr
    pe_nc_fr = (1d0 - H_ion_fr + He_abundance*(1d0 - He_ion_fr))/((2d0 + He_ion_fr2) * He_abundance * He_ion_fr  + 2d0 * H_ion_fr)
    B0 = sqrt(bx0**2 + by0**2)
    petot0 = B0**2/2d0
    beta = 2d0/twofl_gamma

    w(ixO^S,mom_c(1:2))=0d0
    w(ixO^S,mom_n(1:2))=0d0
    
    w(ixO^S,mag(1)) = bx0
    w(ixO^S,mag(2)) = by0
    
    w(ixO^S,rho_c_) = rhotot0/(rho_nc_fr + 1d0)
    w(ixO^S,rho_n_) = rho_nc_fr * w(ixO^S,rho_c_)

    w(ixO^S,e_c_) = petot0/(pe_nc_fr + 1d0)
    w(ixO^S,e_n_) = pe_nc_fr * w(ixO^S,e_c_)


    call twofl_to_conserved(ixI^L,ixO^L,w,x)


  end subroutine initonegrid_usr








end module mod_usr
