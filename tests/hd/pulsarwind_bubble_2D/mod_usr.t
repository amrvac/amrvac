! pulsar wind bubble through SNR and ISM, as in Swaluw et al 2003
module mod_usr
  use mod_hd

  implicit none

  ! input quantities (with dimensions, in cgs)
  !     L0 input luminosity of pulsar
  !     rhoISM input density in cgs of ISM
  !     Lscale sets unit of length
  !     vinf is measure for terminal speed of pulsar wind
  !     Vrel is relative speed of pulsar to shocked SNR medium 
  double precision :: L0,rhoISM,Lscale,vinf,Vrel
  ! derived dimensional quantities
  !     Note: vinf and L0 actually sets Mdotpw
  !           Vsnr derived from Vrel
  !           Vpsr derived from Vrel
  !           Rts termination shock standoff distance
  !           Tcr crossing time for pulsar to meet up with SNR shell
  !           Psh from rhoISM and Vsnr
  !           Soundsh from rhoISM and Vsnr (Psh)
  double precision :: Vsnr,Vpsr,Mdotpw,Psh,Soundsh,Rts,Tcr
  ! scaled or dimensionless quantities
  !     dimensionless size of source region
  double precision :: Rsource

contains

  !> Read this module's parameters from a file
  subroutine usr_params_read(files)
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /usr_list/ L0,rhoISM,Lscale,vinf,Vrel

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, usr_list, end=111)
111    close(unitpar)
    end do

  end subroutine usr_params_read

  subroutine usr_init()
    use mod_global_parameters
    double precision :: eta,Esnr,secinyr

    call usr_params_read(par_files)

    Mdotpw=(2.0d0*L0/vinf**2)
    Vsnr=4.0d0*Vrel/7.0d0
    Vpsr=5.0d0*Vsnr/2.0d0
    Psh=3.0d0*rhoISM*Vsnr**2/4.0d0
    Soundsh=dsqrt(hd_gamma*Psh/(4.0d0*rhoISM))
    eta=0.83d0
    Rts=eta*dsqrt(L0/(2.0d0*dpi*rhoISM*(Vpsr**2)*vinf))
    Esnr=1.0d51
    Tcr=1.27d0*(Esnr/(rhoISM*Vpsr**5))**(1.0d0/3.0d0)
    secinyr=60.0d0*60.0d0*24.0d0*365.0d0

    ! unit length is input (0.1 parsec, in cgs)
    unit_length        = Lscale
    ! unit velocity from sound speed in shocked shell
    unit_velocity      = Soundsh
    ! unit numberdensity from rhoISM (in cgs)
    unit_numberdensity = rhoISM/mp_cgs

    if (mype==0)then
        print *,'unit length        from INPUT Lscale       in cgs=',Lscale
        print *,'unit velocity      from INPUT Vrel,rhoISM  in cgs=',unit_velocity
        print *,'unit numberdensity from INPUT rhoISM       in cgs=',unit_numberdensity
        print *,'unit time should be in cgs=',Lscale/Soundsh
        print *,'unit p    should be in cgs=',rhoISM*Soundsh**2
        print *,'unit rho  should be in cgs=',rhoISM
        print *,'INPUT Mdotpw in cgs=',Mdotpw
        print *,'INPUT L0     in cgs=',L0
        print *,'INPUT Vrel   in cgs=',Vrel
        print *,'INPUT vinf   in cgs=',vinf
        print *,'DERIVE Vsnr      in cgs=',Vsnr
        print *,'DERIVE Vpsr      in cgs=',Vpsr
        print *,'DERIVE Psh       in cgs=',Psh
        print *,'DERIVE Soundsh   in cgs=',Soundsh
        print *,'DERIVE Rts       in cgs=',Rts
        print *,'DERIVE Tcr       in cgs=',Tcr
        print *,'unit time   in yr=',(Lscale/Soundsh)/secinyr
        print *,'DERIVE Tcr  in yr=',Tcr/secinyr
        print *,'DERIVE Tcr  in code units=',Tcr/(Lscale/Soundsh)
    endif

    usr_set_parameters  => initglobaldata_usr
    usr_init_one_grid   => wind_init_one_grid
    usr_special_bc      => specialbound_usr
    usr_refine_grid     => specialrefine_grid
    usr_source          => special_source
    usr_get_dt          => special_getdt
    usr_aux_output      => specialvar_output
    usr_add_aux_names   => specialvarnames_output

    call set_coordinate_system("cylindrical_2D")
    call hd_activate()


  end subroutine usr_init

  subroutine initglobaldata_usr()
    use mod_global_parameters

    hd_gamma=5.0d0/3.0d0

    ! turn source amplitudes dimensionless
    Mdotpw=Mdotpw/(rhoISM*Lscale**2*Soundsh)
    L0=L0/(rhoISM*Lscale**2*(Soundsh)**3)
    Vrel=Vrel/Soundsh
    vinf=vinf/Soundsh
    Vsnr=Vsnr/Soundsh
    Psh=Psh/(rhoISM*Soundsh**2)
    Rts=Rts/Lscale
    if(Rts<dx(1,refine_max_level).or.Rts<dx(2,refine_max_level))then
       if(mype==0) print *,'WARNING:unresolved Rts: may need to increase levels'
    endif
    Rsource=0.05*Rts
    !if(Rsource<dx(1,refine_max_level).or.Rsource<dx(2,refine_max_level))then
    !     print *,Rsource,dx(1,refine_max_level),dx(2,refine_max_level)
    !     print *,'Resetting Rsource to grid size instead of 0.05Rts'
    !     Rsource=2.0d0*max(dx(1,refine_max_level),dx(2,refine_max_level))
    !     print *,'reset Rsource  =',Rsource
    !     if(Rsource>0.1d0*Rts)then
    !        call mpistop('unresolved Rts: increase levels')
    !    endif
    !endif
    if (mype==0)then
        print *,'-------------------'
        print *,'using hd_gamma    =',hd_gamma
        print *,'using He_abundance=',He_abundance
        print *,'-------------------'
        print *,'scaled Mdotpw=',Mdotpw
        print *,'scaled L0    =',L0
        print *,'scaled Vrel  =',Vrel,' SHOULD BE 3.13'
        print *,'scaled vinf  =',vinf
        print *,'scaled Vsnr  =',Vsnr
        print *,'scaled Psh   =',Psh
        print *,'scaled Rts   =',Rts
        print *,'derived Rsource  =',Rsource
        print *,'-------------------'
    endif
    
  if(mype==0) then
      print *, 'unit_density  = ', unit_density
      print *, 'unit_pressure = ', unit_pressure
      print *, 'unit_velocity = ', unit_velocity
      print *, 'unit_time     = ', unit_time
      print *,'--------------------------------'
  end if

  end subroutine initglobaldata_usr

  ! Initialize one grid
  subroutine wind_init_one_grid(ixG^L,ix^L,w,x)
    use mod_global_parameters
    integer, intent(in) :: ixG^L, ix^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)
    
    ! PWN into the (uniform) shocked SNR shell
    !  this shocked SNR shell has rho=4*rhoISM (compression ratio 4)
    w(ix^S,rho_)    = 4.0d0
    w(ix^S,mom(1))  = zero
    !  we set here the relative (dimensionless) velocity
    w(ix^S,mom(2))  = Vrel
    ! the shocked SNR shell has pressure Psh
    w(ix^S,p_)      = Psh
    
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
    if( any(rad(ix^S) <= 4.0d0*Rsource) ) then
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
    case(3)
      ! special bottom boundary: inflowing matter 
      w(ixO^S,rho_)   = 4.0d0
      w(ixO^S,mom(1)) = zero
      w(ixO^S,mom(2)) = Vrel
      w(ixO^S,p_)     = Psh
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

    double precision :: rad(ixI^S),VolSource

    ! use of special source as an internal boundary....
    !   constantly feed mass and (thermal) energy into source region at specified rates
    rad(ixO^S) = dsqrt(x(ixO^S,1)**2 + x(ixO^S,2)**2)
    VolSource=4.0d0*dpi*(Rsource**3)/3.0d0
    where ( rad(ixO^S)< Rsource )
      w(ixO^S,rho_)  = w(ixO^S,rho_)+qdt*Mdotpw/VolSource
      w(ixO^S,e_)    = w(ixO^S,e_)  +qdt*L0/VolSource
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
    w(ixO^S,nw+2)=dlog10(w(ixO^S,rho_))
  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
    character(len=*) :: varnames

    varnames='Te logrho'
  end subroutine specialvarnames_output

  subroutine special_getdt(w,ixG^L,ix^L,dtnew,dx^D,x)

    integer, intent(in) :: ixG^L, ix^L
    double precision, intent(in) :: dx^D, x(ixG^S,1:ndim)
    double precision, intent(in) :: w(ixG^S,1:nw)
    double precision, intent(inout) :: dtnew
    integer :: idims

    dtnew=bigdouble
    do idims=1,ndim 
       dtnew=min(dtnew,0.4d0*minval(block%dx(ix^S,idims)/vinf))
    enddo

  end subroutine special_getdt

end module mod_usr
