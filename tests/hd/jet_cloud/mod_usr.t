module mod_usr
  use mod_hd
  implicit none
  double precision  :: beta, eta_jet, ca, mach, rc

contains

  subroutine usr_init()
  
    call set_coordinate_system("Cartesian_2D")
    
    usr_set_parameters  => initglobaldata_usr
    usr_init_one_grid   => initonegrid_usr
    usr_special_bc      => specialbound_usr
    usr_aux_output      => extra_var_output
    usr_add_aux_names   => extra_var_names_output

    call hd_activate()

  end subroutine usr_init
  
  subroutine initglobaldata_usr
  
    ! jet to cloud density ratio parameter
    beta    = 0.04d0
    ! jet to ambient density
    eta_jet = 3.00d0
    ! ambient sound speed
    ca      = 1.00d0
    ! jet Mach number
    mach    = 10.0d0
    ! cloud to jet radii ratio
    rc      = 1.5d0
    
  end subroutine initglobaldata_usr

  subroutine initonegrid_usr(ixI^L, ixO^L, w, x)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    double precision, intent(inout) :: w(ixI^S, 1:nw)

    double precision                :: rinlet(ixI^S), rcloud(ixI^S)
    double precision                :: x1, x2, x3, sigma


    ! jet comes in from left edge
    rinlet(ixO^S) = abs(x(ixO^S, 2))
    {^IFTHREED
    rinlet(ixO^S) = dsqrt(x(ixO^S, 2)**2+x(ixO^S,3)**2)
    }
    where(rinlet(ixO^S) <= 1.0d0 .and. abs(x(ixO^S, 1) - xprobmin1) <= 2.5d0)
       ! configure jet
       w(ixO^S, rho_)   = 1.0d0
       w(ixO^S, mom(1)) = mach * ca
       w(ixO^S, mom(2)) = 0.0d0
       {^IFTHREED
       w(ixO^S, mom(3)) = 0.0d0
       }
       w(ixO^S, e_)     = ca**2 / (hd_gamma * eta_jet)
    elsewhere
       ! configure ambient medium
       w(ixO^S, rho_)   = 1.0d0 / eta_jet
       w(ixO^S, mom(1)) = 0.0d0
       w(ixO^S, mom(2)) = 0.0d0
       {^IFTHREED
       w(ixO^S, mom(3)) = 0.0d0
       }
       w(ixO^S, e_)     = ca**2 / (hd_gamma * eta_jet)
    endwhere

    ! configure cloud, center coordinates
    x1    = 0.0d0
    x2    = 1.2d0
    x3    = 0.0d0
    sigma = 0.75d0 * rc ! Gaussian width

    rcloud(ixO^S) = (^D&(x(ixO^S,^D)-x^D)**2|+) 

    where(sqrt(rcloud(ixO^S)) <= rc)
       w(ixO^S, rho_) = 1.0d0/eta_jet + (1.0d0/(beta**2)) * exp(-rcloud(ixO^S) / (sigma*sigma))
    endwhere

    call hd_to_conserved(ixI^L, ixO^L, w, x)

  end subroutine initonegrid_usr

  subroutine specialbound_usr(qt, ixI^L, ixO^L, iB, w, x)
    integer, intent(in)             :: ixI^L, ixO^L, iB
    double precision, intent(in)    :: qt, x(ixI^S, 1:ndim)
    double precision, intent(inout) :: w(ixI^S, 1:nw)

    double precision                :: rinlet(ixI^S)

    select case(iB)
    case(1)
       ! fixed left boundary
       rinlet(ixO^S) = abs(x(ixO^S, 2))
       {^IFTHREED
       rinlet(ixO^S) = dsqrt(x(ixO^S, 2)**2+x(ixO^S,3)**2)
       }
       where(rinlet(ixO^S) < 1.0d0)
          w(ixO^S, rho_)   = 1.0d0
          w(ixO^S, mom(1)) = mach * ca
          w(ixO^S, mom(2)) = 0.0d0
          {^IFTHREED
          w(ixO^S, mom(3)) = 0.0d0
          }
          w(ixO^S, e_)     = ca**2 / (hd_gamma * eta_jet)
       elsewhere
          w(ixO^S, rho_)   = 1.0d0/eta_jet
          w(ixO^S, mom(1)) = 0.0d0
          w(ixO^S, mom(2)) = 0.0d0
          {^IFTHREED
          w(ixO^S, mom(3)) = 0.0d0
          }
          w(ixO^S, e_)     = ca**2 / (hd_gamma * eta_jet)
       endwhere
       call hd_to_conserved(ixI^L, ixO^L, w, x)
    case default
       call mpistop('boundary not defined')
    end select

  end subroutine specialbound_usr
  
  subroutine extra_var_output(ixI^L, ixO^L, w, x, normconv)
    use mod_physics
    use mod_global_parameters
    use mod_radiative_cooling
    integer, intent(in)              :: ixI^L, ixO^L
    double precision, intent(in)     :: x(ixI^S, 1:ndim)
    double precision                 :: w(ixI^S, nw+nwauxio), wlocal(ixI^S, 1:nw)
    double precision                 :: normconv(0:nw+nwauxio)

    double precision                 :: pth(ixI^S)
    
    wlocal(ixI^S, 1:nw) = w(ixI^S, 1:nw)

    ! store temperature
    call phys_get_pthermal(wlocal, x, ixI^L, ixO^L, pth)
    w(ixO^S, nw+1) = pth(ixO^S) / w(ixO^S, rho_)
    
  end subroutine extra_var_output
  
  subroutine extra_var_names_output(varnames)
    character(len=*)  :: varnames
    
    varnames = "T"
    
  end subroutine extra_var_names_output
    
end module mod_usr
