module mod_usr
  use mod_mhd
  use mod_eos, only: eos
  use mod_fld
  use mod_multigrid_coupling
    
  implicit none

    double precision:: rholeft,rhoright,slocx1
    double precision:: vleft(3),pleft,vright(3),pright
    double precision:: bx,byleft,bzleft,byright,bzright
    double precision:: Tleft,Tright

  ! Storing additional var in the dat file
  integer :: Tgas_,Trad_,pres_,velx_,vely_,amr_

contains

  subroutine usr_init()

    ! Note how we here must set three values that in turn define M-L-T
    unit_length        =1.001d10   ! cm
    unit_temperature   =1.1d3      ! K
    unit_numberdensity =1.001d10   ! cm^-3

    usr_set_parameters=> initglobaldata_usr

    usr_init_one_grid => initonegrid_usr

    ! Special Boundary conditions for MHD part
    usr_special_bc => boundary_conditions

   ! to add selected variables to the .dat file
    usr_modify_output => set_output_vars

    usr_aux_output      => specialvar_output
    usr_add_aux_names   => specialvarnames_output

    call set_coordinate_system('Cartesian_1.75D')

    call mhd_activate()

    ! to add selected variables to the .dat file
    Tgas_ = var_set_extravar("Tgas", "Tgas")
    Trad_ = var_set_extravar("Trad", "Trad")
    pres_ = var_set_extravar("pres", "pres")
    velx_ = var_set_extravar("velx", "velx")
    vely_ = var_set_extravar("vely", "vely")
    amr_ = var_set_extravar("level", "level")

  end subroutine usr_init

  subroutine initglobaldata_usr
   use mod_global_parameters

    bx=0.75d0 
    rholeft=1.0d0
    pleft=1.0d0 
    Tleft=pleft/(rholeft*RR)
    vleft(1:3)=0.0d0
    byleft=1.0d0
    bzleft=0.0d0
    
    rhoright=0.125d0
    pright=0.1d0
    Tright=pright/(rhoright*RR)
    vright(1:3)=0.0d0
    byright= -1.0d0
    bzright=0.0d0 

    ! here we set the dirichlet (or other) conditions for the MG solver
    mg%bc(1, mg_iphi)%bc_type = mg_bc_dirichlet
    mg%bc(1, mg_iphi)%bc_value = arad_norm*(Tleft**4.d0)
    mg%bc(2, mg_iphi)%bc_type = mg_bc_dirichlet
    mg%bc(2, mg_iphi)%bc_value = arad_norm*(Tright**4.d0)

  end subroutine initglobaldata_usr

  subroutine initonegrid_usr(ixG^L,ix^L,w,x)
    ! initialize one grid 
    integer, intent(in) :: ixG^L,ix^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    logical,save :: first=.true.

    if(first.and.mype==0) then
      print *,'BrioWu test'
      print *,'bx=',bx,' gamma=',eos%gamma
      print *,'LEFT:  by =',byleft,' bz=',bzleft
      print *,'LEFT:  rho=',rholeft,' p=',pleft
      print *,'LEFT:  T  =',Tleft,' Erad=',arad_norm*(Tleft)**4
      print *,'RIGHT: by =',byright,' bz=',bzright
      print *,'RIGHT: rho=',rhoright,' p=',pright
      print *,'RIGHT:  T =',Tright,' Erad=',arad_norm*(Tright)**4
      first=.false.
    endif

    slocx1=half*(xprobmax1+xprobmin1)
    where(x(ixG^S,1)<=slocx1)
       w(ixG^S,rho_)     = rholeft
       w(ixG^S,mom(1))   = rholeft*vleft(1)
       w(ixG^S,mom(2))   = rholeft*vleft(2)
       w(ixG^S,mom(3))   = rholeft*vleft(3)
       w(ixG^S,e_ )      = pleft/(eos%gamma-1.0d0)+half*(rholeft* &
           (vleft(1)**2+vleft(2)**2+vleft(3)**2)+(bx**2+byleft**2+bzleft**2))
       w(ixG^S,mag(1))  = bx
       w(ixG^S,mag(2))  = byleft
       w(ixG^S,mag(3))  = bzleft
       w(ixG^S,r_e)     = arad_norm*(Tleft**4.d0)
    elsewhere
       w(ixG^S,rho_)     = rhoright
       w(ixG^S,mom(1))   = rhoright*vright(1)
       w(ixG^S,mom(2))   = rhoright*vright(2)
       w(ixG^S,mom(3))   = rhoright*vright(3)
       w(ixG^S,e_ )      = pright/(eos%gamma-1.0d0)+half*(rhoright* &
               (vright(1)**2+vright(2)**2+vright(3)**2)+(bx**2+byright**2+bzright**2))
       w(ixG^S,mag(1))  = bx
       w(ixG^S,mag(2))  = byright
       w(ixG^S,mag(3))  = bzright
       w(ixG^S,r_e)     = arad_norm*(Tright**4.d0)
    endwhere

  end subroutine initonegrid_usr

  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)

    double precision :: wlocal(ixI^S,nw)
    double precision :: lamb(ixI^S), R(ixI^S)

    wlocal(ixI^S,1:nw)=w(ixI^S,1:nw)
    call fld_get_fluxlimiter(wlocal,x,ixI^L,ixO^L,lamb,R,2,fld_fl)
    w(ixO^S,nw+1)=lamb(ixO^S)
    w(ixO^S,nw+2)=R(ixO^S)
  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
    character(len=*) :: varnames
    varnames='Lambda R' 
  end subroutine specialvarnames_output


  subroutine boundary_conditions(qdt,qt,ixI^L,ixB^L,iB,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixB^L, iB
    double precision, intent(in)    :: qdt, qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    select case (iB)
    case(1)
       w(ixB^S,rho_)     = rholeft
       w(ixB^S,mom(1))   = rholeft*vleft(1)
       w(ixB^S,mom(2))   = rholeft*vleft(2)
       w(ixB^S,mom(3))   = rholeft*vleft(3)
       w(ixB^S,e_ )      = pleft/(eos%gamma-1.0d0)+half*(rholeft* &
           (vleft(1)**2+vleft(2)**2+vleft(3)**2)+(bx**2+byleft**2+bzleft**2))
       w(ixB^S,mag(1))  = bx
       w(ixB^S,mag(2))  = byleft
       w(ixB^S,mag(3))  = bzleft
       w(ixB^S,r_e)     = arad_norm*(Tleft**4.d0)
    case(2)  
       w(ixB^S,rho_)     = rhoright
       w(ixB^S,mom(1))   = rhoright*vright(1)
       w(ixB^S,mom(2))   = rhoright*vright(2)
       w(ixB^S,mom(3))   = rhoright*vright(3)
       w(ixB^S,e_ )      = pright/(eos%gamma-1.0d0)+half*(rhoright* &
               (vright(1)**2+vright(2)**2+vright(3)**2)+(bx**2+byright**2+bzright**2))
       w(ixB^S,mag(1))  = bx
       w(ixB^S,mag(2))  = byright
       w(ixB^S,mag(3))  = bzright
       w(ixB^S,r_e)     = arad_norm*(Tright**4.d0)
    case default
      call mpistop('boundary not known')
    end select
  end subroutine boundary_conditions

  subroutine set_output_vars(ixI^L,ixO^L,qt,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(in)    :: qt,x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: Trad(ixI^S),Tgas(ixI^S),pth(ixI^S)

    call eos%get_thermal_pressure(w,x,ixI^L,ixO^L,pth)
    call mhd_get_trad(w,x,ixI^L,ixO^L,Trad)
    call eos%get_temperature_from_etot(w,x,ixI^L,ixO^L,Tgas)
    w(ixO^S,Tgas_)=Tgas(ixO^S)
    w(ixO^S,Trad_)=Trad(ixO^S)
    w(ixO^S,pres_)=pth(ixO^S)
    w(ixO^S,velx_)=w(ixO^S,mom(1))/w(ixO^S,rho_)
    w(ixO^S,vely_)=w(ixO^S,mom(2))/w(ixO^S,rho_)
    ! output the AMR level (assuming uniform grid blocks)
    w(ixO^S,amr_)=dlog(((xprobmax1-xprobmin1)/domain_nx1)/dxlevel(1))/dlog(2.0d0)+1.0d0

  end subroutine set_output_vars


end module mod_usr
