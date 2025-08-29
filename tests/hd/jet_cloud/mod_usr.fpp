! inlining to_conservative for performance:
module mod_usr
  use mod_amrvac
  use mod_physics
  
  implicit none
  
  double precision  :: beta, eta_jet, ca, mach, rc
  !$acc declare create( beta, eta_jet, ca, mach, rc )

contains
  
  subroutine usr_init()
  
    call set_coordinate_system("Cartesian_3D")
    
    usr_set_parameters  => initglobaldata_usr
    usr_init_one_grid   => initonegrid_usr
    usr_aux_output      => extra_var_output
    usr_add_aux_names   => extra_var_names_output

    call phys_activate()

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

    !$acc update device( beta, eta_jet, ca, mach, rc )
    
  end subroutine initglobaldata_usr

  subroutine initonegrid_usr(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, x)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nw)

    double precision                :: rinlet(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3), rcloud(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision                :: x1, x2, x3, sigma


    ! jet comes in from left edge
    rinlet(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = abs(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3, 2))
    
    rinlet(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = dsqrt(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3, 2)**2+x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,3)**2)
   
    where(rinlet(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) <= 1.0d0 .and. abs(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3, 1) - xprobmin1) <= 2.5d0)
       ! configure jet
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3, rho_)   = 1.0d0
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3, mom(1)) = mach * ca
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3, mom(2)) = 0.0d0
       
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3, mom(3)) = 0.0d0
      
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_)     = ca**2 / (hd_gamma * eta_jet)
    elsewhere
       ! configure ambient medium
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           rho_)   = 1.0d0 / eta_jet
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3, mom(1)) = 0.0d0
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3, mom(2)) = 0.0d0
       
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3, mom(3)) = 0.0d0
      
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_)     = ca**2 / (hd_gamma * eta_jet)
    endwhere

    ! configure cloud, center coordinates
    x1    = 0.0d0
    x2    = 1.2d0
    x3    = 0.0d0
    sigma = 0.75d0 * rc ! Gaussian width

    rcloud(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = ((x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1)-x1)**2+(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       2)-x2)**2+(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       3)-x3)**2) 

    where(sqrt(rcloud(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)) <= rc)
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           rho_) = 1.0d0/eta_jet + (1.0d0/(beta**2)) * &
          exp(-rcloud(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) / (sigma*sigma))
    endwhere

    call phys_to_conserved(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, x)

  end subroutine initonegrid_usr

  subroutine specialbound_usr(qt, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, iB, w, x)
    !$acc routine vector

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, iB
    double precision, intent(in)    :: qt, x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nw)

    
     double precision                :: rinlet2, inv_gamma_m1
     integer                         :: ix1,ix2,ix3


     select case(iB)

     ! fixed left boundary
     case(1)

        inv_gamma_m1 = 1.0d0/(hd_gamma - 1.0d0)
        
        !$acc loop collapse(3) vector
        do ix3=ixOmin3,ixOmax3
           do ix2=ixOmin2,ixOmax2
              do ix1=ixOmin1,ixOmax1

                 rinlet2 = x(ix1,ix2,ix3, 2)**2 + x(ix1,ix2,ix3, 3)**2

                 if (rinlet2 < 1.0d0) then
                    w(ix1,ix2,ix3, rho_)   = 1.0d0
                    w(ix1,ix2,ix3, mom(1)) = mach * ca
                    w(ix1,ix2,ix3, mom(2)) = 0.0d0
                    w(ix1,ix2,ix3, mom(3)) = 0.0d0
                    w(ix1,ix2,ix3, e_)     = ca**2 / (hd_gamma * eta_jet)
                 else
                    w(ix1,ix2,ix3, rho_)   = 1.0d0 / eta_jet
                    w(ix1,ix2,ix3, mom(1)) = 0.0d0
                    w(ix1,ix2,ix3, mom(2)) = 0.0d0
                    w(ix1,ix2,ix3, mom(3)) = 0.0d0
                    w(ix1,ix2,ix3, e_)     = ca**2 / (hd_gamma * eta_jet)
                 end if

                 ! inline to_conservative() for performance:
                 w(ix1,ix2,ix3, e_) = w(ix1,ix2,ix3, e_) * inv_gamma_m1 + 0.5_dp * w(ix1,ix2,ix3, rho_) * &
                      sum(w(ix1,ix2,ix3,iw_mom(1:ndim))**2)
                 
                 w(ix1,ix2,ix3, iw_mom(1)) = w(ix1,ix2,ix3, iw_rho) * w(ix1,ix2,ix3, iw_mom(1))


                 w(ix1,ix2,ix3, iw_mom(2)) = w(ix1,ix2,ix3, iw_rho) * w(ix1,ix2,ix3, iw_mom(2))


                 w(ix1,ix2,ix3, iw_mom(3)) = w(ix1,ix2,ix3, iw_rho) * w(ix1,ix2,ix3, iw_mom(3))
              end do
           end do
        end do

     end select

   end subroutine specialbound_usr
  
  subroutine extra_var_output(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, x, normconv)
    use mod_physics
    use mod_global_parameters

    integer, intent(in)              :: ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)     :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)
    double precision                 :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw+nwauxio), wlocal(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nw)
    double precision                 :: normconv(0:nw+nwauxio)

    double precision                 :: pth(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    
    wlocal(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
        1:nw) = w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3, 1:nw)

    ! AGILE: tbd
    ! store temperature
!    call phys_get_pthermal(wlocal, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
!       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, pth)
!    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
!        nw+1) = pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
!       ixOmin3:ixOmax3) / w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
!        rho_)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
        nw+1) = 1.0d0
    
  end subroutine extra_var_output
  
  subroutine extra_var_names_output(varnames)
    character(len=*)  :: varnames
    
    varnames = "T"
    
  end subroutine extra_var_names_output
    
end module mod_usr
