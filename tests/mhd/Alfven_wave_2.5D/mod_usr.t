module mod_usr
  use mod_mhd
  implicit none

contains

  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods

    usr_init_one_grid => initonegrid_usr

    call set_coordinate_system("Cartesian_2.5D")

    call mhd_activate()

  end subroutine usr_init

  subroutine initonegrid_usr(ixG^L,ix^L,w,x)
    ! initialize one grid
    use mod_global_parameters

    integer, intent(in) :: ixG^L,ix^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    double precision       :: A, phi(ixG^S), ca, alpha, beta
    double precision       :: v(ixG^S,1:3), b(ixG^S,1:3), vrot(ixG^S,1:3), brot(ixG^S,1:3), k(1:^ND)
    logical, save :: first=.true.

    A     = 0.1d0
    ca    = 1.0d0
    k(1)  = 2.0d0 * dpi

    {^IFONED
    alpha = 0.0d0
    beta  = 0.0d0
    phi(ix^S) = k(1) * x(ix^S,1)
    }

    {^IFTWOD 
    alpha = atan(2.0d0)
    beta  = 0.0d0
    k(2)  = k(1) * tan(alpha)
    phi(ix^S) = k(1) *x(ix^S,1) + k(2) * x(ix^S,2)
    }

    {^IFTHREED
    alpha = atan(2.0d0)
    beta  = atan(2.0d0)
    k(2)  = k(1) * tan(alpha)
    k(3)  = k(1) * tan(beta)
    phi(ix^S) = k(1) *x(ix^S,1) + k(2) * x(ix^S,2) + k(3) * x(ix^S,3)
    }

    w(ix^S,rho_) = 1.0d0
    w(ix^S,p_)   = 0.1d0
    
    v(ix^S,1)  = 0.0d0
    v(ix^S,2)  = A * sin(phi(ix^S))
    v(ix^S,3)  = A * cos(phi(ix^S))
    
    b(ix^S,1)  =   sqrt(w(ix^S,rho_)) * ca
    b(ix^S,2)  = - sqrt(w(ix^S,rho_)) * ca * A * sin(phi(ix^S))
    b(ix^S,3)  = - sqrt(w(ix^S,rho_)) * ca * A * cos(phi(ix^S))

    call rotate(ixG^L,ix^L,v,vrot,alpha,beta)
    call rotate(ixG^L,ix^L,b,brot,alpha,beta)

    w(ix^S,mom(:)) = vrot(ix^S,:)
    w(ix^S,mag(:)) = brot(ix^S,:)

    call mhd_to_conserved(ixG^L,ix^L,w,x)

  end subroutine initonegrid_usr

  subroutine rotate(ixG^L,ix^L,v,vrot,alpha,beta)
    use mod_global_parameters

    integer, intent(in) :: ixG^L,ix^L
    double precision, intent(in) :: v(ixG^S,1:3)
    double precision, intent(out) :: vrot(ixG^S,1:3)
    double precision, intent(in) :: alpha, beta

    double precision       :: gamma

    gamma = atan(cos(alpha)*tan(beta))

    vrot(ix^S,1) = cos(gamma) * cos(alpha) * v(ix^S,1) - sin(alpha) * v(ix^S,2) - sin(gamma) * cos(alpha) * v(ix^S,3)
    vrot(ix^S,2) = cos(gamma) * sin(alpha) * v(ix^S,1) + cos(alpha) * v(ix^S,2) - sin(gamma) * sin(alpha) * v(ix^S,3)
    vrot(ix^S,3) = sin(gamma)              * v(ix^S,1)                          + cos(gamma) * v(ix^S,3)

  end subroutine rotate

end module mod_usr
