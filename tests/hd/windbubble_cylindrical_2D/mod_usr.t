module mod_usr
  use mod_hd

  implicit none

  double precision :: Mdot,vwind,Twind,rhoISM,vISM,TISM,Rstar,Rwind,Tscale,Lscale
  integer :: icase

contains

  !> Read this module's parameters from a file
  subroutine usr_params_read(files)
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /usr_list/ icase

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, usr_list, end=111)
111    close(unitpar)
    end do

  end subroutine usr_params_read

  subroutine usr_init()
    use mod_global_parameters

    call usr_params_read(par_files)

    unit_length        = 3.0857D18
    unit_temperature   = 1.0d7**2.0d0/ (kb_cgs/mp_cgs)
    unit_numberdensity = 10.0**(-25)/mp_cgs

    usr_set_parameters  => initglobaldata_usr
    usr_init_one_grid   => wind_init_one_grid
    usr_special_bc      => specialbound_usr
    usr_refine_grid     => specialrefine_grid
    usr_source          => special_source
    usr_aux_output      => specialvar_output
    usr_add_aux_names   => specialvarnames_output
    usr_var_for_errest  => myvar_for_errest

    call set_coordinate_system("cylindrical_2D")
    call hd_activate()


  end subroutine usr_init

  subroutine initglobaldata_usr()
    use mod_global_parameters

    hd_gamma=5.0d0/3.0d0

select case( icase )
 case(1) ! O-star in cold medium
   Mdot  = 1.0d-6*const_msun/const_years
   vwind = 2.0d8
   Twind = 1.0d4
   rhoISM= (10.0d0)**(-22.5)
   vISM  = 5.0d6
   TISM  = 1.0d2
   Rstar = 5.0d13
   ! this last one is dimensionless
   Rwind = 2.0d-1
 case(2)  ! RSG in cold medium
   Mdot  = 1.0d-4*const_msun/const_years
   vwind = 2.0d6
   Twind = 1.0d3
   rhoISM= (10.0d0)**(-23.5)
   vISM  = 5.0d6
   TISM  = 1.0d2
   Rstar = 5.0d13
   ! this last one is dimensionless
   Rwind = 2.0d-1
 case(3)  ! WR in cold medium
   Mdot  = 1.0d-5*const_msun/const_years
   vwind = 2.5d8
   Twind = 2.0d4
   rhoISM= (10.0d0)**(-22.5)
   vISM  = 5.0d6
   TISM  = 1.0d2
   Rstar = 5.0d13
   ! this last one is dimensionless
   Rwind = 2.0d-1
 case default
     call mpistop("This problem has not been defined")
end select


    length_convert_factor    = unit_length
    w_convert_factor(rho_)   = 10.0**(-25)
    w_convert_factor(mom(1))  = 1.0d7
    w_convert_factor(mom(2))  = w_convert_factor(mom(1))
    w_convert_factor(p_)      = w_convert_factor(rho_)*w_convert_factor(mom(1))*w_convert_factor(mom(1))
    time_convert_factor       = length_convert_factor/w_convert_factor(mom(1))

    Tscale = (1.0D0/(w_convert_factor(mom(1))**2.0d0)) * kb_cgs/mp_cgs
    Lscale =  w_convert_factor(rho_)*time_convert_factor/((mp_cgs*w_convert_factor(mom(1)))**2.0)

    if(mype == 0) then
       write(*,1004) 'time_convert_factor:     ', time_convert_factor
       write(*,1004) 'length_convert_factor:   ', length_convert_factor
       write(*,1004) 'w_convert_factor(mom(1)):', w_convert_factor(mom(1))
       write(*,1004) 'w_convert_factor(rho_):  ', w_convert_factor(rho_)
       write(*,1004) 'w_convert_factor(p_):    ', w_convert_factor(p_)
       write(*,*)
       write(*,1004) 'accel                    ',  &
           w_convert_factor(mom(1))*w_convert_factor(mom(1))/length_convert_factor
       write(*,*)
       write(*,1002) 1.0d0/Tscale
       write(*,1003) Lscale
       write(*,*)
    endif

Rstar = Rstar / length_convert_factor

  if(mype==0) then
      print *, 'unit_density = ', unit_density
      print *, 'unit_pressure = ', unit_pressure
      print *, 'unit_velocity = ', unit_velocity
      print *, 'unit_time = ', unit_time
  end if

1002 format('Temperature unit: ', 1x1pe12.5)
1003 format('Luminosity scale: ', 1x1pe12.5)
1004 format(a25,1x,1pe12.5)
  end subroutine initglobaldata_usr

  ! Initialize one grid
  subroutine wind_init_one_grid(ixG^L,ix^L,w,x)
    use mod_global_parameters
    integer, intent(in) :: ixG^L, ix^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)
    
    double precision :: rad(ixG^S),cosTh(ixG^S),sinTh(ixG^S)

    rad(ix^S) = dsqrt(x(ix^S,1)**2 + x(ix^S,2)**2)
    cosTh(ix^S) = x(ix^S,2)/rad(ix^S)
    sinTh(ix^S) = x(ix^S,1)/rad(ix^S)

    where ( rad(ix^S)>= Rwind )
      w(ix^S,rho_) = rhoISM/w_convert_factor(rho_)
      w(ix^S,mom(1))  = zero
      w(ix^S,mom(2))  = -vISM  /  w_convert_factor(mom(2))
      w(ix^S,p_)   = w(ix^S,rho_)*TISM*Tscale
    elsewhere
      w(ix^S,rho_) = Mdot/(4.0D0*dpi*vwind * (rad(ix^S)*length_convert_factor)**2 ) &
                   / w_convert_factor(rho_)
      w(ix^S,mom(1))  = (vwind /w_convert_factor(mom(1))) * sinTh(ix^S)
      w(ix^S,mom(2))  = (vwind /w_convert_factor(mom(2))) * cosTh(ix^S)
      w(ix^S,p_)     =  w(ix^S,rho_)*Twind*Tscale
    end where

    call hd_to_conserved(ixG^L,ix^L,w,x)

  end subroutine wind_init_one_grid

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

    double precision:: rad(ixG^T)

    rad(ix^S) = dsqrt(x(ix^S,1)**2 + x(ix^S,2)**2)
    if( any(rad(ix^S) <= 1.5*Rwind) ) then
      coarsen = -1
      refine  = 1
    endif

  end subroutine specialrefine_grid

  subroutine specialbound_usr(qt,ixG^L,ixO^L,iB,w,x)
    use mod_global_parameters
    integer, intent(in) :: ixG^L, ixO^L, iB
    double precision, intent(in) :: qt, x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    select case(iB)
    case(4)
      w(ixO^S,rho_)   = rhoISM/w_convert_factor(rho_)
      w(ixO^S,mom(1)) = zero
      w(ixO^S,mom(2)) = -vISM/w_convert_factor(mom(2))
      w(ixO^S,p_)     = w(ixO^S,rho_)*TISM*Tscale
      call hd_to_conserved(ixG^L,ixO^L,w,x)
    case default
      call mpistop("This boundary is not supposed to be special")
    end select

  end subroutine specialbound_usr

  subroutine special_source(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)
    use mod_global_parameters
    integer, intent(in) :: ixI^L, ixO^L, iw^LIM
    double precision, intent(in) :: qdt, qtC, qt
    double precision, intent(in) :: x(ixI^S,1:ndim), wCT(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: rad(ixI^S), cosTh(ixI^S), sinTh(ixI^S)

    ! use of special source as an internal boundary....

    rad(ixO^S) = dsqrt(x(ixO^S,1)**2 + x(ixO^S,2)**2)
    cosTh(ixO^S) = x(ixO^S,2)/rad(ixO^S)
    sinTh(ixO^S) = x(ixO^S,1)/rad(ixO^S)

    where ( rad(ixO^S)< Rwind )
      w(ixO^S,rho_)  = Mdot/(4.0D0*dpi*vwind* (rad(ixO^S)*length_convert_factor)**2 ) &
                      / w_convert_factor(rho_)
      w(ixO^S,mom(1)) = (vwind /w_convert_factor(mom(1))) *sinTh(ixO^S)*w(ixO^S,rho_)
      w(ixO^S,mom(2)) = (vwind /w_convert_factor(mom(1))) *cosTh(ixO^S)*w(ixO^S,rho_)
      w(ixO^S,e_)    = w(ixO^S,rho_)*Twind*Tscale/(hd_gamma-one)+ &
           half*(w(ixO^S,mom(1))**2.0d0+w(ixO^S,mom(2))**2.0d0)/w(ixO^S,rho_)
    end where

  end subroutine special_source

  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)

    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)

    double precision :: pth(ixI^S),wlocal(ixI^S,1:nw)

    wlocal(ixI^S,1:nw)=w(ixI^S,1:nw)
    call hd_get_pthermal(wlocal,x,ixI^L,ixO^L,pth)
    w(ixO^S,nw+1)=pth(ixO^S)/w(ixO^S,rho_)
  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
    character(len=*) :: varnames

    varnames='Te'
  end subroutine specialvarnames_output

  subroutine myvar_for_errest(ixI^L,ixO^L,iflag,w,x,var)
      use mod_global_parameters
      integer, intent(in)           :: ixI^L,ixO^L,iflag
      double precision, intent(in)  :: w(ixI^S,1:nw), x(ixI^S,1:ndim)
      double precision, intent(out) :: var(ixI^S)

      if (iflag >nw+1)call mpistop(' iflag error')
      var(ixO^S) = dsqrt(w(ixO^S,mom(1))**2 + w(ixO^S,mom(2))**2 )/w(ixO^S,rho_)

  end subroutine myvar_for_errest

end module mod_usr
