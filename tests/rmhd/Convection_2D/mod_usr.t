module mod_usr
  use mod_rmhd
  use mod_viscosity, only: vc_mu
  implicit none
  double precision :: usr_grav,bstr,temptop

contains

  subroutine usr_init()
    call set_coordinate_system("Cartesian")

    usr_set_parameters  => initglobaldata_usr
    usr_init_one_grid   => initonegrid_usr
    usr_special_bc      => specialbound_usr
    usr_special_mg_bc   => mg_boundary_conditions 
    usr_gravity         => gravity
    usr_refine_grid     => specialrefine_grid
    usr_aux_output      => specialvar_output
    usr_add_aux_names   => specialvarnames_output
    usr_print_log       => energies_log

    unit_length        = 1.001d8 ! cm
    unit_temperature   = 1.001d3 ! K
    unit_numberdensity = 1.001d9 ! cm^-3

    call rmhd_activate()
  end subroutine usr_init

  subroutine initglobaldata_usr
    integer:: mpoly
    double precision:: zeta0,rhat,sigma,zz0,qchand
    double precision:: gamma, qchi, qmpoly, eta2
    logical, save:: first=.true.

    if (mype == 0.and.first)then
     print *,"unit_numberdensity = ",unit_numberdensity
     print *,"unit_temperature = ",unit_temperature
     print *,"unit_length = ",unit_length
     print *,"He_abundance = ",He_abundance
     print *,"He_ion_fr = ",He_ion_fr
     print *,"kB_cgs = ",kB_cgs
     print *,"mp_cgs = ",mp_cgs
     print *,"const_c = ",const_c
     print *,"c_norm = ",c_norm
     print *,"unit_density = ",unit_density
     print *,"unit_pressure = ",unit_pressure
     print *,"unit_velocity = ",unit_velocity
     print *,"unit_time = ",unit_time
     print *,"unit_magneticfield = ",unit_magneticfield
     print *,"unit_radflux = ",unit_radflux
     print *,"unit_opacity = ",unit_opacity
    endif

    ! This setup relates to Hurlburt and Toomre
    ! and calculates the eqpar array from problem-specific input parameters
    ! ---------------------------------------------------------------------
    !
    ! problem specific input parameters for
    ! Hurlburt and Toomre magnetoconvection problem
    !
    ! we need (apart from the aspect ratio which is set in par-file):
    ! i)equilibrium parameters: (1) ratio of specific heats (gamma)
    ! ------------------------- (2) the polytropic index mpoly (m)
    !                            --> sets gravity g
    !                           (3) the density contrast chi
    !                            --> sets top temperature zz0 (z0)
    !                           (4) the Chandrasekhar number qchand (Q)
    !                            --> initial field strength B (need nu and eta)
    ! ii)parameters setting the dissipative coefficients:
    ! ---------------------------------------------------
    !    (1) Prandtl number sigma (viscous/thermal conduction)
    !    (2) the Rayleigh number at half-depth rhat (degree of instability)
    !    (3) the magnetic Prandtl number zeta0
    !      (magnetic diffusion/thermal conduction)
    !
    ! hence fixing   i) gamma, mpoly, chi, qchand
    !               ii) sigma, rhat, zeta0
    !
    ! allows calculation of all equation parameters, namely
    ! the general purpose ones:
    ! ------------------------
    !  rmhd_gamma mhd_eta tc_k_para eqpar(grav1.2_) vc_mu
    !
    ! and the problem specific one:
    ! ----------------------------
    !  temptop for setting the top temperature
    !  bstr for setting the B field
    !
    
    ! equilibrium parameters
    gamma=1.66666667d0
    mpoly=1
    qchi=11.0d0
    qchand=72.0d0
    
    qmpoly=dble(mpoly)
    zz0=one/(qchi**(one/qmpoly)-one)
    
    if (mype==0.and.first) then
      write(*,*)'gamma, mpoly, qchi,qchand'
      write(*,*)gamma, mpoly, qchi, qchand
      write(*,*)'Deduced top dimensionless temperature z0=',zz0
    endif
    
    temptop=zz0
    rmhd_gamma=gamma
    usr_grav=-(qmpoly+one)
    
    ! dissipative parameters
    sigma=1.0d0
    rhat=1.0d5
    zeta0=0.25d0
    
    if (mype==0.and.first) then
      write(*,*) 'sigma,rhat,zeta0:'
      write(*,*) sigma,rhat,zeta0
    endif
    
    eta2=(qmpoly+one)*(gamma/(gamma-1)-(qmpoly+one))*((gamma-one)/gamma)*&
         ((zz0+half)**(two*qmpoly-one)/zz0**(two*qmpoly))*&
         (zeta0**2/(sigma*rhat))
    
    if(eta2<=smalldouble)then
       if(mype==0.and.first) write(*,*) 'eta2=',eta2
       call mpistop("Negative or too small value for eta**2")
    endif
    
    rmhd_eta=dsqrt(eta2)
    vc_mu=rmhd_eta*sigma/zeta0
    tc_fl%tc_k_para=(gamma/(gamma-one))*rmhd_eta/zeta0
    bstr=dsqrt(qchand*vc_mu*rmhd_eta)
    
    if (mype==0.and.first) then
      write(*,*)'dimensionless values for dissipative coefficients:'
      write(*,*)'resistivity          eta=',rmhd_eta
      write(*,*)'viscosity             mu=',vc_mu
      write(*,*)'thermal conduction tc_k_para=',tc_fl%tc_k_para
      write(*,*)'dimensionless magnetic field strength:',bstr
      first=.false.
    endif

  end subroutine initglobaldata_usr

  subroutine initonegrid_usr(ixG^L,ix^L,w,x)
    use mod_global_parameters 
    integer, intent(in) :: ixG^L,ix^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    double precision :: dvx,dvy,dvz,nkx,nky,nkz,zz0,qmpoly
    logical, save :: first=.true.

    ! velocity perturbation parameters
    dvx=1.0d-4
    dvy=1.0d-4
    
    nkx=two*dpi/3.0d0
    nky=two*dpi
    {^IFTHREED
    dvz=1.0d-4
    nkz=two*two*dpi
    }
    if(first)then
      if(mype==0) then
         write(*,*)'Simulating rmhd convection: Hurlburt-Toomre'
         write(*,*) 'velocity perturbation: dvx,dvy,nkx,nky'
         write(*,*) dvx,dvy,nkx,nky
         {^IFTHREED
         write(*,*) '3D setup!'
         write(*,*) dvz,nkz
         }
      endif
    endif
    zz0=temptop
    qmpoly=-one-usr_grav
    w(ix^S,rho_)=((zz0+one-x(ix^S,2))/zz0)**qmpoly
    w(ix^S,mag(1))=zero
    w(ix^S,mag(2))=bstr
    {^IFTHREED w(ix^S,mag(3))=zero }
    w(ix^S,p_)= zz0*(((zz0+one-x(ix^S,2))/zz0)**(qmpoly+one))
    w(ix^S,mom(1))=dvx*dsin(x(ix^S,1)*nkx)*dsin(x(ix^S,2)*nky){^IFTHREED *dsin(x(ix^S,3)*nkz)}
    w(ix^S,mom(2))=dvy*dsin(x(ix^S,1)*nkx)*dsin(x(ix^S,2)*nky){^IFTHREED *dsin(x(ix^S,3)*nkz)}
    {^IFTHREED w(ix^S,mom(3))=zero }
    w(ix^S,r_e) = const_rad_a*((zz0+one-x(ix^S,2))*unit_temperature)**4.d0/unit_pressure 
    call rmhd_to_conserved(ixG^L,ix^L,w,x)

    if(mype == 0.and.first)then
       print *,"radn energy magnitude ~ ", const_rad_a*(unit_temperature)**4.d0/unit_pressure 
       print *,"opacity magnitude ~ ", 0.4/unit_opacity 
       print *,"unit_numberdensity = ",unit_numberdensity
       print *,"unit_temperature = ",unit_temperature
       print *,"unit_length = ",unit_length
       print *,"unit_density = ",unit_density
       print *,"unit_pressure = ",unit_pressure
       print *,"unit_velocity = ",unit_velocity
       print *,"unit_time = ",unit_time
       print *,"unit_magneticfield = ",unit_magneticfield
       print *,"unit_radflux = ",unit_radflux
       print *,"unit_opacity = ",unit_opacity
       print *,"He_abundance = ",He_abundance
       print *,"eq_state_units = ",eq_state_units
       print *,"H_ion_fr = ",H_ion_fr
       print *,"He_ion_fr2 = ",He_ion_fr2
       first=.false.
    endif
  end subroutine initonegrid_usr

  subroutine specialbound_usr(qt,ixG^L,ixO^L,iB,w,x)
    use mod_global_parameters 
    ! special boundary types, user defined
    integer, intent(in) :: ixO^L, iB, ixG^L
    double precision, intent(in) :: qt, x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    double precision :: tempb(ixG^S)
    double precision :: delydelx,delydelz
    integer :: ix^D,ixIM^L
    
    select case(iB)
    case(3)
      ! special bottom boundary
      ! ensure fixed temperature gradient in ghost layers
      ! use asymm/symm for B-components (ensuring vertical field)
      ! use divB correction from central difference for vertical field
      ! use symm/asymm for v-components (ensuring no flow-through)
      ! density profile: use symmetry
      {^IFTWOD
      ! in nghostcells rows above the bottom boundary: switch to primitive
      ixIMmin2=ixOmax2+1;ixIMmax2=ixOmax2+nghostcells;
      ixIMmin1=ixOmin1;ixIMmax1=ixOmax1;
      call rmhd_to_primitive(ixG^L,ixIM^L,w,x)
      do ix2=ixOmax2,ixOmin2,-1
         w(ixOmin1:ixOmax1,ix2,rho_)= w(ixOmin1:ixOmax1,2*ixOmax2+1-ix2,rho_)
         w(ixOmin1:ixOmax1,ix2,mom(1)) = w(ixOmin1:ixOmax1,2*ixOmax2+1-ix2,mom(1))
         w(ixOmin1:ixOmax1,ix2,mom(2)) =-w(ixOmin1:ixOmax1,2*ixOmax2+1-ix2,mom(2))
         w(ixOmin1:ixOmax1,ix2,mag(1)) =-w(ixOmin1:ixOmax1,2*ixOmax2+1-ix2,mag(1))
         w(ixOmin1:ixOmax1,ix2,mag(2)) = w(ixOmin1:ixOmax1,2*ixOmax2+1-ix2,mag(2))
         !case1
         w(ixOmin1:ixOmax1,ix2,r_e)= w(ixOmin1:ixOmax1,2*ixOmax2+1-ix2,r_e) 
      enddo
      ! fill temperature array: extrapolate linearly with fixed dT/dy=-1, from bottom row
      do ix2=ixOmin2,ixOmax2
        tempb(ixOmin1:ixOmax1,ix2)=w(ixOmin1:ixOmax1,ixOmax2+1,p_)/w(ixOmin1:ixOmax1,ixOmax2+1,rho_) &
                                 -(x(ixOmin1:ixOmax1,ix2,2)-x(ixOmin1:ixOmax1,ixOmax2+1,2))
      enddo
      ! determine pressure
      w(ixO^S,p_)=w(ixO^S,rho_)*tempb(ixO^S)
      }
      {^IFTHREED
      ! in nghostcells rows above the bottom boundary: switch to primitive
      ixIMmin2=ixOmax2+1;ixIMmax2=ixOmax2+nghostcells;
      ixIMmin1=ixOmin1;ixIMmax1=ixOmax1;
      ixIMmin3=ixOmin3;ixIMmax3=ixOmax3;
      call rmhd_to_primitive(ixG^L,ixIM^L,w,x)
      do ix2=ixOmax2,ixOmin2,-1
         w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,rho_)= w(ixOmin1:ixOmax1,2*ixOmax2+1-ix2,ixOmin3:ixOmax3,rho_)
         w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,mom(1)) = w(ixOmin1:ixOmax1,2*ixOmax2+1-ix2,ixOmin3:ixOmax3,mom(1))
         w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,mom(2)) =-w(ixOmin1:ixOmax1,2*ixOmax2+1-ix2,ixOmin3:ixOmax3,mom(2))
         w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,mom(3)) = w(ixOmin1:ixOmax1,2*ixOmax2+1-ix2,ixOmin3:ixOmax3,mom(3))
         w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,mag(1)) =-w(ixOmin1:ixOmax1,2*ixOmax2+1-ix2,ixOmin3:ixOmax3,mag(1))
         w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,mag(2)) = w(ixOmin1:ixOmax1,2*ixOmax2+1-ix2,ixOmin3:ixOmax3,mag(2))
         w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,mag(3)) =-w(ixOmin1:ixOmax1,2*ixOmax2+1-ix2,ixOmin3:ixOmax3,mag(3))
      enddo
      ! fill temperature array: extrapolate linearly with fixed dT/dy=-1, from bottom row
      do ix2=ixOmin2,ixOmax2
        tempb(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3)=w(ixOmin1:ixOmax1,ixOmax2+1,ixOmin3:ixOmax3,p_) &
                                                  /w(ixOmin1:ixOmax1,ixOmax2+1,ixOmin3:ixOmax3,rho_) &
               -(x(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,2)-x(ixOmin1:ixOmax1,ixOmax2+1,ixOmin3:ixOmax3,2))
      enddo
      ! determine pressure
      w(ixO^S,p_)=w(ixO^S,rho_)*tempb(ixO^S)
      }
      ! now reset the inner mesh values to conservative
      call rmhd_to_conserved(ixG^L,ixIM^L,w,x)
      ! now switch to conservative in full bottom ghost layer
      call rmhd_to_conserved(ixG^L,ixO^L,w,x)
    case(4)
      ! special top boundary
      ! ensure the fixed temperature in ghost layers
      ! use asymm/symm for B-components (ensuring vertical field)
      ! use divB correction from central difference for vertical field
      ! use symm/asymm for v-components (ensuring no flow-through)
      ! density profile: use symmetry
    
      {^IFTWOD
      ! in nghostcells rows below top boundary: switch to primitive
      ixIMmin2=ixOmin2-nghostcells;ixIMmax2=ixOmin2-1;
      ixIMmin1=ixOmin1;ixIMmax1=ixOmax1;
      call rmhd_to_primitive(ixG^L,ixIM^L,w,x)
      do ix2=ixOmin2,ixOmax2,+1
          w(ixOmin1:ixOmax1,ix2,rho_)= w(ixOmin1:ixOmax1,2*ixOmin2-ix2-1,rho_)
          w(ixOmin1:ixOmax1,ix2,mom(1)) = w(ixOmin1:ixOmax1,2*ixOmin2-ix2-1,mom(1))
          w(ixOmin1:ixOmax1,ix2,mom(2)) =-w(ixOmin1:ixOmax1,2*ixOmin2-ix2-1,mom(2))
          w(ixOmin1:ixOmax1,ix2,mag(1)) =-w(ixOmin1:ixOmax1,2*ixOmin2-ix2-1,mag(1))
          w(ixOmin1:ixOmax1,ix2,mag(2)) = w(ixOmin1:ixOmax1,2*ixOmin2-ix2-1,mag(2))
      enddo
      w(ixO^S,p_)=w(ixO^S,rho_)*temptop
      !case1
      w(ixO^S,r_e)= const_rad_a*(temptop*unit_temperature)**4.d0/unit_pressure 
    
      }
      {^IFTHREED
      ! in nghostcells rows below top boundary: switch to primitive
      ixIMmin2=ixOmin2-nghostcells;ixIMmax2=ixOmin2-1;
      ixIMmin1=ixOmin1;ixIMmax1=ixOmax1;
      ixIMmin3=ixOmin3;ixIMmax3=ixOmax3;
      call rmhd_to_primitive(ixG^L,ixIM^L,w,x)
      do ix2=ixOmin2,ixOmax2,+1
          w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,rho_)= w(ixOmin1:ixOmax1,2*ixOmin2-ix2-1,ixOmin3:ixOmax3,rho_)
          w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,mom(1)) = w(ixOmin1:ixOmax1,2*ixOmin2-ix2-1,ixOmin3:ixOmax3,mom(1))
          w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,mom(2)) =-w(ixOmin1:ixOmax1,2*ixOmin2-ix2-1,ixOmin3:ixOmax3,mom(2))
          w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,mom(3)) = w(ixOmin1:ixOmax1,2*ixOmin2-ix2-1,ixOmin3:ixOmax3,mom(3))
          w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,mag(1)) =-w(ixOmin1:ixOmax1,2*ixOmin2-ix2-1,ixOmin3:ixOmax3,mag(1))
          w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,mag(2)) = w(ixOmin1:ixOmax1,2*ixOmin2-ix2-1,ixOmin3:ixOmax3,mag(2))
          w(ixOmin1:ixOmax1,ix2,ixOmin3:ixOmax3,mag(3)) =-w(ixOmin1:ixOmax1,2*ixOmin2-ix2-1,ixOmin3:ixOmax3,mag(3))
      enddo
      w(ixO^S,p_)=w(ixO^S,rho_)*temptop
      }
      ! now reset the inner mesh values to conservative
      call rmhd_to_conserved(ixG^L,ixIM^L,w,x)
      ! now switch to conservative in full bottom ghost layer
      call rmhd_to_conserved(ixG^L,ixO^L,w,x)
    case default
       call mpistop("Special boundary is not defined for this region")
    end select

  end subroutine specialbound_usr

  subroutine mg_boundary_conditions(iB)
    use mod_global_parameters
    use mod_multigrid_coupling

    integer, intent(in)             :: iB

    !case1
    select case (iB)
    case (3)
      mg%bc(iB, mg_iphi)%bc_type = mg_bc_neumann
      mg%bc(iB, mg_iphi)%bc_value = 0.0 !gradE
      
    case (4)
      mg%bc(iB, mg_iphi)%bc_type = mg_bc_dirichlet
      mg%bc(iB, mg_iphi)%bc_value = (const_rad_a*(temptop*unit_temperature)**4)/unit_pressure

    case default
      print *, "Not a standard: ", typeboundary(r_e, iB)
      error stop "Set special bound for this Boundary "
    end select
  end subroutine mg_boundary_conditions
    

  ! Calculate gravitational acceleration in each dimension
  subroutine gravity(ixI^L,ixO^L,wCT,x,gravity_field)
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(in)    :: wCT(ixI^S,1:nw)
    double precision, intent(out)   :: gravity_field(ixI^S,ndim)

    gravity_field(ixO^S,:)=0.d0
    gravity_field(ixO^S,2)=usr_grav

  end subroutine gravity

  subroutine specialrefine_grid(igrid,level,ixG^L,ix^L,qt,w,x,refine,coarsen)
    ! Enforce additional refinement or coarsening
    ! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.
    integer, intent(in) :: igrid, level, ix^L, ixG^L
    double precision, intent(in) :: qt, w(ixG^S,1:nw), x(ixG^S,1:ndim)
    integer, intent(inout) :: refine, coarsen

    if (any(x(ix^S,2)<=xprobmin2+0.05d0)) then
      refine=1
      coarsen=-1
    endif
  end subroutine specialrefine_grid

  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
    ! this subroutine can be used in convert, to add auxiliary variables to the
    ! converted output file, for further analysis using tecplot, paraview, ....
    ! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
    !
    ! the array normconv can be filled in the (nw+1:nw+nwauxio) range with
    ! corresponding normalization values (default value 1)
    use mod_fld
    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)

    double precision :: tmp(ixI^S),divb(ixI^S),current(ixI^S,7-2*ndir:3), lamb(ixO^S), R(ixO^S)
    double precision :: radflux(ixO^S,1:ndim)
    integer          :: idirmin

    ! output Te
    call rmhd_get_pthermal(w,x,ixI^L,ixO^L,tmp)
    w(ixO^S,nw+1)=tmp(ixO^S)/w(ixO^S,rho_)
    ! output B 
    w(ixO^S,nw+2)=dsqrt(sum(w(ixO^S,mag(:))**2,dim=ndim+1))
    ! output divB1
    call get_divb(w,ixI^L,ixO^L,divb)
    w(ixO^S,nw+3)=divb(ixO^S)
    ! output the plasma beta p*2/B**2
    w(ixO^S,nw+4)=tmp(ixO^S)*two/sum(w(ixO^S,mag(:))**2,dim=ndim+1)
    ! store current
    call get_current(w,ixI^L,ixO^L,idirmin,current)
    w(ixO^S,nw+5)=current(ixO^S,3)

    call fld_get_fluxlimiter(w,x,ixI^L,ixO^L,lamb,R,1)
    w(ixO^S,nw+6)=lamb(ixO^S)
    w(ixO^S,nw+7)=R(ixO^S)
    call fld_get_radflux(w,x,ixI^L,ixO^L,radflux,1)
    w(ixO^S,nw+8)=radflux(ixO^S,1)
    w(ixO^S,nw+9)=radflux(ixO^S,2)

  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
    ! newly added variables need to be concatenated with the w_names/primnames string
    character(len=*) :: varnames

    varnames='Te B divB beta j3 Lambda R F1 F2'
  end subroutine specialvarnames_output


  subroutine printlog_energy_averages
  
    use mod_input_output
    use mod_timing
    use mod_forest, only: nleafs, nleafs_active, nleafs_level
    use mod_global_parameters

    logical              :: fileopen
    integer              :: i, iw, level
    double precision     :: wmean(1:nw), total_volume
    double precision     :: volume_coverage(refine_max_level)
    integer              :: nx^D, nc, ncells, dit
    double precision     :: dtTimeLast, now, cellupdatesPerSecond
    double precision     :: activeBlocksPerCore, wctPerCodeTime, timeToFinish
    character(len=40)    :: fmt_string
    
    !character(len=*), parameter :: fmt_i  = 'i8'     ! Integer format
    !character(len=*), parameter :: fmt_r  = 'es16.8' ! Default precision
    !character(len=*), parameter :: fmt_r2 = 'es10.2' ! Two digits
    
    character(len=80)    :: filename
    character(len=2048)  :: line
    logical, save        :: opened  = .false.
    integer              :: amode, istatus(MPI_STATUS_SIZE)
    integer, parameter   :: my_unit = 20

    double precision :: KineticEnergy_avg, MagneticEnergy_avg, RadiationEnergy_avg, InternalEnergy_avg, Current_avg

    ! Compute the volume-average of w**1 = w
    call get_volume_average(1, wmean, total_volume)

    call get_volume_average_func(kinetic, KineticEnergy_avg, total_volume)
    call get_volume_average_func(magnetic, MagneticEnergy_avg, total_volume)
    call get_volume_average_func(radiation, RadiationEnergy_avg, total_volume)
    call get_volume_average_func(internal, InternalEnergy_avg, total_volume)

    ! Compute the volume coverage
    call get_volume_coverage(volume_coverage)

    if (mype == 0) then

       ! To compute cell updates per second, we do the following:
       nx^D=ixMhi^D-ixMlo^D+1;
       nc={nx^D*}
       ncells = nc * nleafs_active

       ! assumes the number of active leafs haven't changed since last compute.
       now        = MPI_WTIME()
       dit        = it - itTimeLast
       dtTimeLast = now - timeLast
       itTimeLast = it
       timeLast   = now
       cellupdatesPerSecond = dble(ncells) * dble(nstep) * &
            dble(dit) / (dtTimeLast * dble(npe))

       ! blocks per core:
       activeBlocksPerCore = dble(nleafs_active) / dble(npe)

       ! Wall clock time per code time unit in seconds:
       wctPerCodeTime = dtTimeLast / max(dit * dt, epsilon(1.0d0))

       ! Wall clock time to finish in hours:
       timeToFinish = (time_max - global_time) * wctPerCodeTime / 3600.0d0

       ! On first entry, open the file and generate the header
       if (.not. opened) then

          filename = trim(base_filename) // "_energy_avgs.log"

          ! Delete the log when not doing a restart run
          if (restart_from_file == undefined) then
             open(unit=my_unit,file=trim(filename),status='replace')
             close(my_unit, status='delete')
          end if

          amode    = ior(MPI_MODE_CREATE,MPI_MODE_WRONLY)
          amode    = ior(amode,MPI_MODE_APPEND)

          call MPI_FILE_OPEN(MPI_COMM_SELF, filename, amode, &
               MPI_INFO_NULL, log_fh, ierrmpi)

          opened = .true.

          ! Start of file headern
          line = "it global_time dt"
          do level=1,nw
             i = len_trim(line) + 2
             write(line(i:),"(a,a)") trim(cons_wnames(level)), " "
          end do

          ! Volume coverage per level
          do level = 1, refine_max_level
             i = len_trim(line) + 2
             write(line(i:), "(a,i0)") "c", level
          end do

          ! Cell counts per level
          do level=1,refine_max_level
             i = len_trim(line) + 2
             write(line(i:), "(a,i0)") "n", level
          end do

          ! Rest of file header
          line = trim(line) // " | Xload Xmemory 'Cell_Updates /second/core'"
          line = trim(line) // " 'Active_Blocks/Core' 'Wct Per Code Time [s]'"
          line = trim(line) // " 'TimeToFinish [hrs]'"

          ! Only write header if not restarting
          if (restart_from_file == undefined .or. reset_time) then
            call MPI_FILE_WRITE(log_fh, trim(line) // new_line('a'), &
                 len_trim(line)+1, MPI_CHARACTER, istatus, ierrmpi)
          end if
       end if

       ! Construct the line to be added to the log

       fmt_string = '(' // fmt_i // ',2' // fmt_r // ')'
       write(line, fmt_string) it, global_time, dt
       i = len_trim(line) + 2

       write(fmt_string, '(a,i0,a)') '(', nw, fmt_r // ')'
       write(line(i:), fmt_string) wmean(1:nw)
       i = len_trim(line) + 2

       write(fmt_string, '(a,i0,a)') '(', refine_max_level, fmt_r // ')'
       write(line(i:), fmt_string) volume_coverage(1:refine_max_level)
       i = len_trim(line) + 2

       write(fmt_string, '(a,i0,a)') '(', refine_max_level, fmt_i // ')'
       write(line(i:), fmt_string) nleafs_level(1:refine_max_level)
       i = len_trim(line) + 2

       fmt_string = '(a,6' // fmt_r2 // ')'
       write(line(i:), fmt_string) '| ', Xload, Xmemory, cellupdatesPerSecond, &
            activeBlocksPerCore, wctPerCodeTime, timeToFinish

       call MPI_FILE_WRITE(log_fh, trim(line) // new_line('a') , &
            len_trim(line)+1, MPI_CHARACTER, istatus, ierrmpi)
    end if

  end subroutine printlog_energy_averages


  subroutine energies_log
   use mod_global_parameters
   use mod_input_output, only: printlog_default, get_volume_average_func

   integer, parameter :: my_unit = 123
   character(len=80)  :: fmt_string = '(6(es12.4))', filename
   logical, save      :: alive, visited=.false.
   double precision   :: total_volume

   double precision :: KineticEnergy_avg, MagneticEnergy_avg, RadiationEnergy_avg, InternalEnergy_avg, J_max

   ! First output the standard log file:
   call printlog_default

   ! Now make the custom _c.log file:
   ! Calculate average temperature:
   call get_volume_average_func(kinetic, KineticEnergy_avg, total_volume)
   call get_volume_average_func(magnetic, MagneticEnergy_avg, total_volume)
   call get_volume_average_func(radiation, RadiationEnergy_avg, total_volume)
   call get_volume_average_func(internal, InternalEnergy_avg, total_volume)
   call get_max_current_magnitude(J_max)

   filename = trim(base_filename) // "_c.log"

   if (.not. visited) then
     ! Delete the log when not doing a restart run
     if (restart_from_file == undefined) then
        open(unit=my_unit,file=trim(filename),form='formatted',status='replace')
        write(my_unit,'(a)') '#Global_time KE_avg Mag_avg Rad_evg IE_avg J_max'
     end if
     visited = .true.
  end if

  if (mype == 0) then
     write(filename,"(a)") filename
     inquire(file=filename,exist=alive)
     if(alive) then
       open(unit=my_unit,file=filename,form='formatted',status='old',access='append')
     else
       open(unit=my_unit,file=filename,form='formatted',status='new')
     endif

     ! if number of output doubles is increased, don't forget to change the fmt_string above
      write(my_unit, fmt_string) global_time, KineticEnergy_avg, MagneticEnergy_avg, RadiationEnergy_avg, InternalEnergy_avg, J_max
      close(my_unit)
    end if


  end subroutine energies_log




  ! Function that calculates kinetic energy, to be used in get_volume_average_func
  pure function kinetic(w_vec, w_size) result(kin_energy)
    integer, intent(in)          :: w_size
    double precision, intent(in) :: w_vec(w_size)
    double precision             :: kin_energy

    kin_energy = 0.5*(w_vec(mom(1))**2 + w_vec(mom(2))**2) / w_vec(rho_) ! divide by rho for conserved
  end function kinetic

  ! Function that calculates magnetic energy, to be used in get_volume_average_func
  pure function magnetic(w_vec, w_size) result(mag_energy)
    integer, intent(in)          :: w_size
    double precision, intent(in) :: w_vec(w_size)
    double precision             :: mag_energy

    mag_energy = 0.5*(w_vec(mag(1))**2 + w_vec(mag(2))**2)
  end function magnetic

  ! Function that calculates radiation energy, to be used in get_volume_average_func
  pure function radiation(w_vec, w_size) result(rad_energy)
    integer, intent(in)          :: w_size
    double precision, intent(in) :: w_vec(w_size)
    double precision             :: rad_energy

    rad_energy = w_vec(r_e)
  end function radiation

  ! Function that calculates internal energy, to be used in get_volume_average_func
  pure function internal(w_vec, w_size) result(int_energy)
    integer, intent(in)          :: w_size
    double precision, intent(in) :: w_vec(w_size)
    double precision             :: int_energy

    int_energy = w_vec(e_)-0.5*(w_vec(mom(1))**2 + w_vec(mom(2))**2) / w_vec(rho_)-0.5*(w_vec(mag(1))**2 + w_vec(mag(2))**2)
  end function internal


  subroutine get_max_current_magnitude(J_max)

   use mod_global_parameters
   double precision, intent(out) :: J_max

   integer                       :: iigrid, igrid
   double precision              :: J_max_mype, J_max_recv

   double precision :: current(ixG^T,7-2*ndir:3)
   double precision :: tmp(ixG^T)
   integer :: idirmin

   J_max_mype = -bigdouble

  ! Loop over all the grids and keep track of the temporary min/max
   do iigrid = 1, igridstail
      igrid = igrids(iigrid)

      call get_current(ps(igrid)%w,ixG^LL,ixM^LL,idirmin,current)
      tmp(ixM^T) = sum(current(ixM^T,:)**2,dim=ndim+1)
      J_max_mype = max(J_max_mype,maxval(tmp(ixM^T)))
    end do


   call mpi_allreduce(J_max_mype, J_max_recv, 1, mpi_double_precision, &
        mpi_max, icomm, ierrmpi)

   J_max = J_max_recv

  end subroutine get_max_current_magnitude


end module mod_usr
