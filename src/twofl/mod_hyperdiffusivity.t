!> Subroutines for Roe-type Riemann solver for HD
module mod_hyperdiffusivity


  implicit none

  !all public
  !private
  !public :: hyperdiffusivity_init


contains

  subroutine hyperdiffusivity_init()
    use mod_global_parameters
    use mod_geometry
    use mod_physics, only: phys_req_diagonal

    !print*, slab_uniform !!THIS IS FALSE, why??

    !if(coordinate .ne. Cartesian .or. .not. slab_uniform) then
    if(coordinate .ne. Cartesian) then
      call mpistop("Hyperdiffusivity only implemented for Cartesian uniform grid")
    endif

    nghostcells = max(nghostcells,3)
    phys_req_diagonal = .true.

  end subroutine hyperdiffusivity_init


  subroutine hyp_coeff(ixI^L, ixO^L, var, idimm, nu_hyp)
    use mod_global_parameters
    integer, intent(in)                :: ixI^L, idimm
    integer, intent(out)                :: ixO^L
    double precision, intent(in)       :: var(ixI^S)
    double precision, intent(out)       :: nu_hyp(ixI^S)

    double precision        :: tmp(ixI^S), tmp2(ixI^S)
    integer :: hx^L, hx1f^L, hx1b^L, hx2b^L

    double precision, parameter :: eps=1e-8


   hxmin^D=ixImin^D+2*kr(idimm,^D); 
   hxmax^D=ixImax^D-kr(idimm,^D); 

   hx1f^L=hx^L+kr(idimm,^D);  
   hx1b^L=hx^L-kr(idimm,^D);  
   hx2b^L=hx^L-2*kr(idimm,^D);  

  !store d3 in tmp
  tmp(hx^S)= abs(3d0 * (var(hx^S) - var(hx1b^S)) - (var(hx1f^S)-var(hx2b^S)))
  !store d1 in tmp2
  tmp2(hx^S)= abs((var(hx^S) - var(hx1b^S)))

  ixOmin^D=hxmin^D+kr(idimm,^D);
  ixOmax^D=hxmax^D-kr(idimm,^D);
  hx1f^L=ixO^L+kr(idimm,^D);  
  hx1b^L=ixO^L-kr(idimm,^D);  

  nu_hyp(ixO^S) = max(tmp(hx1b^S),tmp(ixO^S),tmp(hx1f^S))/max(tmp2(hx1b^S),tmp2(ixO^S),tmp2(hx1f^S),eps)  

  !print*, "HYP IXO ", ixO^L 

  end subroutine hyp_coeff



  subroutine div_vel_coeff(ixI^L, ixO^L, vel, idimm, nu_vel)
    use mod_global_parameters
    integer, intent(in)                :: ixI^L, idimm
    integer, intent(out)                :: ixO^L
    double precision, intent(in)       :: vel(ixI^S,1:ndir)
    double precision, intent(out)       :: nu_vel(ixI^S)
    integer :: hx1f^L, hx1b^L

    integer :: ii

   ixOmin^D=ixImin^D+1; 
   ixOmax^D=ixImax^D-1; 

   ixOmin^D=ixOmin^D+kr(idimm,^D); 

   nu_vel(ixO^S) = 0d0 

  ! ndim should always be smaller or equal to ndir !
   do ii = 1, ndim
    hx1b^L=ixO^L-kr(ii,^D);
    if(ii==idimm) then
      nu_vel(ixO^S) = nu_vel(ixO^S) + (vel(ixO^S,ii) - vel(hx1b^S,ii))/dxlevel(ii) 
    else
      hx1f^L=ixO^L+kr(ii,^D);  
      nu_vel(ixO^S) = nu_vel(ixO^S) + (vel(hx1f^S,ii) - vel(hx1b^S,ii))/(2d0*dxlevel(ii)) 
    endif

    where(nu_vel(ixO^S) < 0d0)
      nu_vel(ixO^S) = -nu_vel(ixO^S)
    elsewhere
      nu_vel(ixO^S) = 0d0
    endwhere

   enddo 

  !print*, "DIV_VEL IXO ", ixO^L 

  end subroutine div_vel_coeff



  !var has cell centered values
  !nu_hyper is defined at the interfaces
  subroutine second_same_deriv(ixI^L, ixO^L, nu_hyper, var, idimm, res)
    use mod_global_parameters
    integer, intent(in)                :: ixI^L, idimm
    integer, intent(out)                :: ixO^L
    double precision, intent(in)       :: nu_hyper(ixI^S),var(ixI^S)
    double precision, intent(out)       :: res(ixI^S)

    integer :: hxf^L, hxb^L

    ixOmin^D=ixImin^D+3;
    ixOmax^D=ixImax^D-3;

    hxf^L=ixO^L+kr(idimm,^D);
    hxb^L=ixO^L-kr(idimm,^D);

    res(ixO^S) = 1d0/(dxlevel(idimm)**2)*&
          (nu_hyper(hxf^S) * (var(hxf^S)-var(ixO^S))-&
          nu_hyper(ixO^S) * (var(ixO^S)-var(hxb^S)))
   !print*, "SECOND SAME DERIV IXO ", ixO^L 

  end subroutine second_same_deriv

  !var has cell centered values
  !var2 has cell centered values
  !nu_hyper is defined at the interfaces
  subroutine second_same_deriv2(ixI^L, ixO^L, nu_hyper, var2, var, idimm, res)
    use mod_global_parameters
    integer, intent(in)                :: ixI^L, idimm
    integer, intent(out)                :: ixO^L
    double precision, intent(in)       :: nu_hyper(ixI^S),var(ixI^S), var2(ixI^S)
    double precision, intent(out)       :: res(ixI^S)

    integer :: hxf^L, hxb^L

    ixOmin^D=ixImin^D+3;
    ixOmax^D=ixImax^D-3;

    hxf^L=ixO^L+kr(idimm,^D);
    hxb^L=ixO^L-kr(idimm,^D);

    res(ixO^S) = 1d0/(2d0*dxlevel(idimm)**2)*&
          (nu_hyper(hxf^S) * (var2(hxf^S)+var2(ixO^S)) *  (var(hxf^S)-var(ixO^S))-&
          nu_hyper(ixO^S) * (var2(hxb^S)+var2(ixO^S)) * (var(ixO^S)-var(hxb^S)))
    !print*, "SECOND SAME DERIV2 IXO ", ixO^L 

  end subroutine second_same_deriv2

  !idimm inner derivative, idimm2 outer
  ! deriv_idimm2(nu * deriv_idimm (u) )
  !var has cell centered values
  !nu_hyper is defined at the interfaces
  subroutine second_cross_deriv(ixI^L, ixO^L, nu_hyper, var, idimm, idimm2, res)
    use mod_global_parameters
    integer, intent(in)                :: ixI^L, idimm,idimm2
    integer, intent(out)                :: ixO^L
    double precision, intent(in)       :: nu_hyper(ixI^S),var(ixI^S)
    double precision, intent(out)       :: res(ixI^S)

    integer :: hxfi^L, hxbi^L, hxfifj^L, hxbibj^L, hxfibj^L, hxbifj^L

    ixOmin^D=ixImin^D+3;
    ixOmax^D=ixImax^D-3;

    hxfi^L=ixO^L+kr(idimm2,^D);
    hxbi^L=ixO^L-kr(idimm2,^D);

    hxfifj^L=hxfi^L+kr(idimm,^D);
    hxfibj^L=hxfi^L-kr(idimm,^D);

    hxbifj^L=hxbi^L+kr(idimm,^D);
    hxbibj^L=hxbi^L-kr(idimm,^D);

    res(ixO^S) = 1d0/(8d0*dxlevel(idimm) * dxlevel(idimm2))*&
          ((nu_hyper(hxfifj^S) + nu_hyper(hxfi^S)) * (var(hxfifj^S)-var(hxfibj^S))-&
          (nu_hyper(hxbifj^S) + nu_hyper(hxbi^S)) * (var(hxbifj^S)-var(hxbibj^S)))

    !print*, "SECOND CROSS DERIV IXO ", ixO^L 

  end subroutine second_cross_deriv


  !idimm inner derivative, idimm2 outer
  ! deriv_idimm2(nu * var2 * deriv_idimm (u) )
  !var has cell centered values
  !var2 has cell centered values
  !nu_hyper is defined at the interfaces (with numbering index left center): center_i-1 interface_i center_i
  subroutine second_cross_deriv2(ixI^L, ixO^L, nu_hyper, var2, var, idimm, idimm2, res)
    use mod_global_parameters
    integer, intent(in)                :: ixI^L, idimm,idimm2
    integer, intent(out)                :: ixO^L
    double precision, intent(in)       :: nu_hyper(ixI^S),var(ixI^S),var2(ixI^S)
    double precision, intent(out)       :: res(ixI^S)

    integer :: hxfi^L, hxbi^L, hxfifj^L, hxbibj^L, hxfibj^L, hxbifj^L

    ixOmin^D=ixImin^D+3;
    ixOmax^D=ixImax^D-3;

    hxfi^L=ixO^L+kr(idimm2,^D);
    hxbi^L=ixO^L-kr(idimm2,^D);

    hxfifj^L=hxfi^L+kr(idimm,^D);
    hxfibj^L=hxfi^L-kr(idimm,^D);

    hxbifj^L=hxbi^L+kr(idimm,^D);
    hxbibj^L=hxbi^L-kr(idimm,^D);

    res(ixO^S) = 1d0/(8d0*dxlevel(idimm) * dxlevel(idimm2))*&
          ((nu_hyper(hxfifj^S) + nu_hyper(hxfi^S)) * var2(hxfi^S) * (var(hxfifj^S)-var(hxfibj^S))-&
          (nu_hyper(hxbifj^S) + nu_hyper(hxbi^S)) * var2(hxbi^S) * (var(hxbifj^S)-var(hxbibj^S)))

    !print*, "SECOND CROSS DERIV2 IXO ", ixO^L 

  end subroutine second_cross_deriv2

end module mod_hyperdiffusivity
