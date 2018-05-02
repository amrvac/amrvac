module mod_usr
  use mod_hd
  implicit none
  double precision :: rhoBC, velBC, gm, MdotIni, BernIni, rs

contains

  subroutine usr_init()

    call set_coordinate_system('spherical')

    usr_set_parameters=> initglobaldata_usr
    usr_init_one_grid => initonegrid_usr
    usr_source        => pt_grav_source
    usr_get_dt        => get_dt_pt_grav
    usr_special_bc    => specialbound_usr
    usr_refine_grid   => specialrefine_grid
    usr_aux_output    => specialvar_output
    usr_add_aux_names => specialvarnames_output

    call hd_activate()

  end subroutine usr_init

  subroutine initglobaldata_usr

    ! hd_gamma=5.d0/3.d0
    ! hd_adiab=one/hd_gamma
    ! mach=4.d0
   !  vel=one
    gm=half
   !  rho0=one
   !  Racc=one
    ! pini=((1/hd_gamma)*rho0*vel**two)/mach**two
    if (hd_gamma==one) then
      rhoBC=(dpi*dexp(1.5d0)/4.d0) / (4.d0*dpi)
    else
      rhoBC=( (dpi/4.d0) * (2.d0/(5.d0-3.d0*hd_gamma))**((5.d0-3.d0*hd_gamma)/(2.d0*(hd_gamma-1.d0))) ) / (4.d0*dpi)
    endif
    velBC=one ! = Mdot / ( 4pi*r2 * rhoBC)
    MdotIni=(dpi/4.d0)*(2.d0/(5.d0-3.d0*hd_gamma))**((5.d0-3.d0*hd_gamma)/(2.d0*hd_gamma-2.d0))
    BernIni=one/(hd_gamma-one)
    rs=(5.d0-3.d0*hd_gamma)/8.d0

  end subroutine initglobaldata_usr

  subroutine initonegrid_usr(ixG^L,ix^L,w,x)

    ! initialize one grid


    integer, intent(in) :: ixG^L, ix^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    double precision :: rin
    integer :: ix, iy

    if (hd_gamma==one) then
      w(ix^S,rho_)   =  (rhoBC/(xprobmax1**1.5d0))/10000.d0
      w(ix^S,mom(1)) = -((rhoBC/(xprobmax1**1.5d0))/10000.d0)*(velBC/dsqrt(xprobmin1))
    else
      rin=1.0277248004393015E-002 ! 1.0304870605468750d-2 ! ** to modify **
      open(123,file='compute_analytic_initial_state/iniState.dat')
      print*, x(3,1)
      if (stretched_grid) then
        do ix=1,domain_nx1/block_nx1
           ! ** to modify **
           print*, ix, rin*(xprobmax1/xprobmin1)**((dble((ix-1)*block_nx1))/dble(domain_nx1))
           if (dabs(rin*(xprobmax1/xprobmin1)**((dble((ix-1)*block_nx1))/dble(domain_nx1))-x(3,1))/x(3,1)<0.0001d0) then
              exit
           endif
        enddo
      else
        do ix=1,domain_nx1/block_nx1
          !  print*, ix, rin+(xprobmax1-xprobmin1)*((dble((ix-1)*block_nx1))/dble(domain_nx1))
           if (dabs(rin+(xprobmax1-xprobmin1)*((dble((ix-1)*block_nx1))/dble(domain_nx1))-x(3,1))/x(3,1)<0.0001d0) then
              exit
           endif
        enddo
      endif
      if (ix>1) then
      do iy=1,(ix-1)*block_nx1
          read(123,*)
      enddo
      endif
      do ix=ixmin1,ixmax1
         read(123,*) rin, w(ix,rho_), w(ix,mom(1))
      enddo
      close(123)
    endif

  end subroutine initonegrid_usr

  subroutine pt_grav_source(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

    integer, intent(in) :: ixI^L, ixO^L, iw^LIM
    double precision, intent(in) :: qdt, qtC, qt
    double precision, intent(in) :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    w(ixO^S,mom(1))=w(ixO^S,mom(1)) - &
         qdt * wCT(ixO^S,rho_) * (gm/x(ixO^S,1)**two)

  end subroutine pt_grav_source

  subroutine get_dt_pt_grav(w,ixG^L,ix^L,dtnew,dx^D,x)

    integer, intent(in) :: ixG^L, ix^L
    double precision, intent(in) :: dx^D, x(ixG^S,1:ndim)
    double precision, intent(in) :: w(ixG^S,1:nw)
    double precision, intent(inout) :: dtnew

    dtnew=bigdouble
    dtnew=min(dtnew,minval(dsqrt(block%dx(ix^S,1)/(gm/x(ix^S,1)**two))))

  end subroutine get_dt_pt_grav

  subroutine specialbound_usr(qt,ixG^L,ixB^L,iB,w,x)

    ! special boundary types, user defined

    integer, intent(in) :: ixG^L, ixB^L, iB
    double precision, intent(in) :: qt, x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)
    double precision :: rho(ixG^S)

    double precision :: rho_tmp(ixG^S)
    double precision :: eps=1.d-10, df, v, f
    integer :: ix

    select case(iB)
    case(1)
       do ix=ixBmin1,ixBmax1
          ! Cont grad
          rho_tmp(ix)=w(ixBmax1+1,rho_)+&
             (w(ixBmax1+2,rho_)-w(ixBmax1+1,rho_))*&
             ((x(ix      ,1)-x(ixBmax1+1,1))/&
             (x(ixBmax1+2,1)-x(ixBmax1+1,1)))
          w(ix,rho_)=rho_tmp(ix)
          ! rho*vr*r^2
          w(ix,mom(1)) = min(zero,w(ixBmax1+1,mom(1))*&
               ( x(ixBmax1+1,1)/x(ix,1) )**two)
       enddo

    case(2)

    if (hd_gamma==one) then

      do ix=ixBmin1,ixBmax1
         if (x(ix,1)<0.25d0) then
            v  = 5.*velBC/dsqrt(x(ix,1)*4.d0) ! first guess
         elseif (x(ix,1)>0.25d0) then
            v  =   velBC/dsqrt(x(ix,1)*4.d0) ! first guess
         endif
         f  = half*v**2.-dlog(v)-2.*dlog(x(ix,1))-1./(2.d0*x(ix,1)) + (3.d0/2.d0-2.d0*dlog(4.d0))
         do while(abs(f)>eps)
            df = v-1./v
            if (v-f/df>0.d0) then
               v  = v - f/df
            else
               v = v/2.
            endif
            f  = half*v**2.-dlog(v)-2.*dlog(x(ix,1))-1./(2.d0*x(ix,1)) + (3.d0/2.d0-2.d0*dlog(4.d0))
         enddo
         w(ix,mom(1) )=-rhoBC/x(ix,1)**two
         w(ix,rho_)=(rhoBC/x(ix,1)**two) / dabs(v)
      enddo

    else

       ! no str, gm=1.05, 16384
       !1.0496525209920142d0
       !1.0496311096571693d0
       ! no str, gm=1.05, 16384
       !1.0496525209920142d0
       !1.0496311096571693d0
       ! code by hand the values output interactively by compute_solution_for_any_gamma
       ! ** to modify **
       w(ixBmin1,rho_)= 1.0498505790046455d0 !1.0498383583500150d0
       w(ixBmax1,rho_)= 1.0471712764685581d0 !1.0471712764685581d0 !1.0471603272669170d0
       ! Slope set and used for extrapolation
       ! w(ixBmin1,rho_)= w(ixBmin1-1,rho_)+(x(ixBmin1,1)-x(ixBmin1-1,1))*((1.0471161492737657d0-1.0497890536683687d0)/(x(ixBmax1,1)-x(ixBmin1,1)))
       ! w(ixBmax1,rho_)= w(ixBmin1-1,rho_)+(x(ixBmax1,1)-x(ixBmin1-1,1))*((1.0471161492737657d0-1.0497890536683687d0)/(x(ixBmax1,1)-x(ixBmin1,1)))
       w(ixBmin1,mom(1) )=-MdotIni/(4.d0*dpi*x(ixBmin1,1)**two) !-2.6221506513911630d-3
       w(ixBmax1,mom(1) )=-MdotIni/(4.d0*dpi*x(ixBmax1,1)**two) !-2.3538708990921897d-3

       ! do ix=ixBmin1,ixBmax1
       !    w(ix,mom(1) )= - ( (one/16.d0) * ((two/(5.d0-3.d0*gamma))**((5.d0-3.d0*gamma)/(two*(gamma-one)))) ) / x(ix,1)**two
       !    w(ix,rho_)=w(ix-1,rho_)
       ! enddo

    endif

    case default
      call mpistop("Special boundary is not defined for this region")
    end select

  end subroutine specialbound_usr

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

    if (rs>x(ixmin1,1) .and. rs<x(ixmax1,1)) then
       refine=1
    else
       refine=-1
    endif

  end subroutine specialrefine_grid

  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
  integer, intent(in)                :: ixI^L,ixO^L
  double precision, intent(in)       :: x(ixI^S,1:ndim)
  double precision                   :: w(ixI^S,nw+nwauxio)
  double precision                   :: normconv(0:nw+nwauxio)

   w(ixO^S,nw+1) = dabs(w(ixO^S,mom(1))/w(ixO^S,rho_))/dsqrt(hd_gamma*hd_adiab*w(ixO^S,rho_)**(hd_gamma-one))
   w(ixO^S,nw+2) = dabs(dabs(w(ixO^S,mom(1))*x(ixO^S,1)**two)*4.d0*dpi-MdotIni)/MdotIni
   w(ixO^S,nw+3) = dabs(half*w(ixO^S,mom(1))**two/w(ixO^S,rho_)**two-gm/x(ixO^S,1)+(one/(hd_gamma-one))*w(ixO^S,rho_)**(hd_gamma-one)-BernIni)/BernIni
   
  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
    character(len=*) :: varnames
    varnames='Mach Mdot Berny'

  end subroutine specialvarnames_output

end module mod_usr
