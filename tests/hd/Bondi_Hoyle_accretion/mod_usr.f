module mod_usr
  use mod_hd
  implicit none
  double precision :: mach, rho0, vel, Racc, gm, pini, pw

contains

  subroutine usr_init()
    use mod_usr_methods

    call set_coordinate_system('spherical_2D')

    usr_set_parameters=> initglobaldata_usr
    usr_init_one_grid => initonegrid_usr
    usr_source        => pt_grav_source
    usr_get_dt        => get_dt_pt_grav
    usr_special_bc    => specialbound_usr
    usr_internal_bc   => intern_bc

    call hd_activate()

  end subroutine usr_init

  subroutine initglobaldata_usr
    
    use mod_global_parameters
    
    hd_gamma=5.d0/3.d0
    mach=4.d0
    vel=one
    gm=half
    rho0=one
    Racc=one
    pini=((1/hd_gamma)*rho0*vel**two)/mach**two

  end subroutine initglobaldata_usr

  subroutine initonegrid_usr(ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixmin1,ixmin2,&
     ixmax1,ixmax2,w,x)

    ! initialize one grid 

    use mod_global_parameters

    integer, intent(in) :: ixGmin1,ixGmin2,ixGmax1,ixGmax2, ixmin1,ixmin2,&
       ixmax1,ixmax2
    double precision, intent(in) :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:nw)

    w(ixmin1:ixmax1,ixmin2:ixmax2,rho_)   =  rho0
    w(ixmin1:ixmax1,ixmin2:ixmax2,mom(1)) = -rho0*vel*dcos(x(ixmin1:ixmax1,&
       ixmin2:ixmax2,2))
    w(ixmin1:ixmax1,ixmin2:ixmax2,mom(2)) =  rho0*vel*dsin(x(ixmin1:ixmax1,&
       ixmin2:ixmax2,2))
    w(ixmin1:ixmax1,ixmin2:ixmax2,e_)     =  &
       pini/(hd_gamma-one)+half*rho0*vel**two

  end subroutine initonegrid_usr

  subroutine pt_grav_source(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,iwmin,iwmax,qtC,wCT,qt,w,x)

    use mod_global_parameters

    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2, iwmin,iwmax
    double precision, intent(in) :: qdt, qtC, qt
    double precision, intent(in) :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,1:nw),&
        x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(1))=w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,mom(1)) - qdt * wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       rho_)    * (gm/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)**two)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_ )=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       e_ ) - qdt * wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       mom(1))  * (gm/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)**two)

  end subroutine pt_grav_source

  subroutine get_dt_pt_grav(w,ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixmin1,ixmin2,&
     ixmax1,ixmax2,dtnew,dx1,dx2,x)

    use mod_global_parameters

    integer, intent(in) :: ixGmin1,ixGmin2,ixGmax1,ixGmax2, ixmin1,ixmin2,&
       ixmax1,ixmax2
    double precision, intent(in) :: dx1,dx2, x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:nw),&
        dtnew

    dtnew=bigdouble
    dtnew=min(dtnew,minval(dsqrt(block%dx(ixmin1:ixmax1,ixmin2:ixmax2,&
       1)/(gm/x(ixmin1:ixmax1,ixmin2:ixmax2,1)**two))))

  end subroutine get_dt_pt_grav

  subroutine specialbound_usr(qt,ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixBmin1,&
     ixBmin2,ixBmax1,ixBmax2,iB,w,x)

    ! special boundary types, user defined

    use mod_global_parameters

    integer, intent(in) :: ixGmin1,ixGmin2,ixGmax1,ixGmax2, ixBmin1,ixBmin2,&
       ixBmax1,ixBmax2, iB
    double precision, intent(in) :: qt, x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:nw)
    double precision :: rho(ixGmin1:ixGmax1,ixGmin2:ixGmax2)
    
    ! Cont grad
      rho(ixBmax1,ixBmin2:ixBmax2)=w(ixBmax1+1,ixBmin2:ixBmax2,&
         rho_)+(w(ixBmax1+2,ixBmin2:ixBmax2,rho_)-w(ixBmax1+1,ixBmin2:ixBmax2,&
         rho_))*((x(ixBmax1      ,ixBmin2:ixBmax2,1)-x(ixBmax1+1,&
         ixBmin2:ixBmax2,1))/(x(ixBmax1+2,ixBmin2:ixBmax2,1)-x(ixBmax1+1,&
         ixBmin2:ixBmax2,1)))
      w(ixBmax1,ixBmin2:ixBmax2,rho_)=rho(ixBmax1,ixBmin2:ixBmax2)
      ! rho*vr*r^2
      w(ixBmax1,ixBmin2:ixBmax2,mom(1)) = min(zero,w(ixBmax1+1,ixBmin2:ixBmax2,&
         mom(1))*( x(ixBmax1+1,ixBmin2:ixBmax2,1)/x(ixBmax1,ixBmin2:ixBmax2,&
         1) )**two)
      ! rho*vr*vth*r^3 cont (ie r*vth cont)
      w(ixBmax1,ixBmin2:ixBmax2,mom(2)) = rho(ixBmax1,&
         ixBmin2:ixBmax2) * (w(ixBmax1+1,ixBmin2:ixBmax2,mom(2))/w(ixBmax1+1,&
         ixBmin2:ixBmax2,rho_))*(x(ixBmax1+1,ixBmin2:ixBmax2,1)/x(ixBmax1,&
         ixBmin2:ixBmax2,1))
      ! (e+P+rho*phi)*v_r*r^2    cont
      w(ixBmax1,ixBmin2:ixBmax2,e_) = (one/hd_gamma)*( (hd_gamma-one)*half*( ( &
         (min(zero,w(ixBmax1+1,ixBmin2:ixBmax2,mom(1))*( x(ixBmax1+1,&
         ixBmin2:ixBmax2,1)/x(ixBmax1,ixBmin2:ixBmax2,&
         1) )**two))**two+(rho(ixBmax1,ixBmin2:ixBmax2) * (w(ixBmax1+1,&
         ixBmin2:ixBmax2,mom(2))/w(ixBmax1+1,ixBmin2:ixBmax2,&
         rho_))*(x(ixBmax1+1,ixBmin2:ixBmax2,1)/x(ixBmax1,ixBmin2:ixBmax2,&
         1)))**two )/rho(ixBmax1,ixBmin2:ixBmax2))+ rho(ixBmax1,&
         ixBmin2:ixBmax2)*(gm/dabs(x(ixBmax1,ixBmin2:ixBmax2,&
         1))) +(rho(ixBmax1,ixBmin2:ixBmax2)/w(ixBmax1+1,ixBmin2:ixBmax2,&
         rho_))*(hd_gamma*w(ixBmax1+1,ixBmin2:ixBmax2,&
         e_)-(hd_gamma-one)*half*((w(ixBmax1+1,ixBmin2:ixBmax2,&
         mom(1))**two+w(ixBmax1+1,ixBmin2:ixBmax2,mom(2))**two)/w(ixBmax1+1,&
         ixBmin2:ixBmax2,rho_)) -w(ixBmax1+1,ixBmin2:ixBmax2,&
         rho_)*(gm/dabs(x(ixBmax1+1,ixBmin2:ixBmax2,1))) ) )

   w(ixBmin1,ixBmin2:ixBmax2,:)=w(ixBmax1,ixBmin2:ixBmax2,:)

  end subroutine specialbound_usr

  subroutine intern_bc(level,qt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,w,x)

    use mod_global_parameters

    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,level
    double precision, intent(in) :: qt
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision :: R1, R2, rho(ixImin1:ixImax1,ixImin2:ixImax2),&
        zeta(ixImin1:ixImax1,ixImin2:ixImax2)

    R1=half*Racc ! xprobmax1/(335./20.)
    R2=one *Racc ! xprobmax1/(335./35.)
    zeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=half*x(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,1)*dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       2))*( one + dsqrt(one+4.d0*(gm/vel**two)*((one-dcos(x(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,2)))/(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       1)*dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2))**two))) )
    rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=rho0*(zeta(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)**two/(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       1)*dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2))*(two*zeta(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)-x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       1)*dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2)))))
    ! Is it ahead the shock? THEN FIX! No need to compare to B-K solution
    where ( x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       2)<dpi/two .AND. x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       1)>R2/(1+((R2-R1)/R1)*dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2))) )
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)= rho(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(1)) = -rho0*(zeta(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)**two/(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          1)*dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          2))*(two*zeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2)-x(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,1)*dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          2)))))*dsqrt(vel**two+(two*gm)/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          1)-((zeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2)*vel)/x(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,1))**two)
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(2)) = rho0*(zeta(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)**two/(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          1)*dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          2))*(two*zeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2)-x(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,1)*dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          2)))))*((zeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2)*vel)/x(&
          ixOmin1:ixOmax1,ixOmin2:ixOmax2,1))
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          e_)  = ((rho0*vel**two)/hd_gamma)*(one/(mach**two*(hd_gamma-one))+&
          half)*(zeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2)**two/(x(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,1)*dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          2))*(two*zeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2)-x(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,1)*dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          2)))))+(rho0/hd_gamma)*(zeta(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)**two/(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          1)*dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          2))*(two*zeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2)-x(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,1)*dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          2)))))*(half*(hd_gamma-one)* (dsqrt(vel**two+&
          (two*gm)/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)-((zeta(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)*vel)/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          1))**two)**two+((zeta(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)*vel)/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          1))**two) + gm/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1))
    endwhere

  end subroutine intern_bc

end module mod_usr
