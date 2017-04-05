!##############################################################################
!> module thermal conduction for HD and MHD
!> 10.07.2011 developed by Chun Xia and Rony Keppens
!> 01.09.2012 moved to modules folder by Oliver Porth
!> 13.10.2013 optimized further by Chun Xia
!> 12.03.2014 implemented RKL2 super timestepping scheme to reduce iterations 
!> and improve stability and accuracy up to second order in time by Chun Xia.
!> 23.08.2014 implemented saturation and perpendicular TC by Chun Xia
!> 12.01.2017 modulized by Chun Xia
!> 
!> PURPOSE: 
!> IN MHD ADD THE HEAT CONDUCTION SOURCE TO THE ENERGY EQUATION
!> S=DIV(KAPPA_i,j . GRAD_j T)
!> where KAPPA_i,j = kappa b_i b_j + kappe (I - b_i b_j)
!> b_i b_j = B_i B_j / B**2, I is the unit matrix, and i, j= 1, 2, 3 for 3D
!> IN HD ADD THE HEAT CONDUCTION SOURCE TO THE ENERGY EQUATION
!> S=DIV(kappa . GRAD T)
!> USAGE:
!> 1. in mod_usr.t -> subroutine usr_init(), add 
!>        unit_length=<your length unit>
!>        unit_numberdensity=<your number density unit>
!>        unit_velocity=<your velocity unit>
!>        unit_temperature=<your temperature unit>
!>    before call (m)hd_activate()
!> 2. to switch on thermal conduction in the (m)hd_list of amrvac.par add:
!>    (m)hd_thermal_conduction=.true.
!> 3. in the tc_list of amrvac.par :
!>    tc_perpendicular=.true.  ! (default .false.) turn on thermal conduction perpendicular to magnetic field 
!>    tc_saturate=.false.  ! (default .true. ) turn off thermal conduction saturate effect
!>    tc_dtpar=0.9/0.45/0.3 ! stable time step coefficient for 1D/2D/3D, decrease it for more stable run

module mod_thermal_conduction

  implicit none
  !> Coefficient of thermal conductivity
  double precision, public :: kappa

  !> Coefficient of thermal conductivity perpendicular to magnetic field
  double precision, public :: kappe

  !> Time step of thermal conduction
  double precision :: dt_tc

  !> Number of sub-steps of supertime stepping
  integer, public :: s

  !> Index of the density (in the w array)
  integer, private, parameter              :: rho_ = 1

  !> Indices of the momentum density
  integer, allocatable, private, protected :: mom(:)

  !> Index of the energy density (-1 if not present)
  integer, private, protected              :: e_

  !> Indices of the magnetic field
  integer, allocatable, private, protected :: mag(:)

  !> The adiabatic index
  double precision, private :: tc_gamma

  !> The smallest allowed energy
  double precision, private :: smalle

  !> The smallest allowed density
  double precision, private :: minrho

  !> The smallest allowed pressure
  double precision, private :: minp

  !> Time step coefficient
  double precision, private :: tc_dtpar=0.9d0

  !> Maximal number of sub-cycles within one fuild time step
  integer :: ncyclemax=1000

  !> Calculate thermal conduction perpendicular to magnetic field (.true.) or not (.false.)
  logical, private :: tc_perpendicular=.false.

  !> Consider thermal conduction saturation effect (.true.) or not (.false.)
  logical, private :: tc_saturate=.true.

  !> Logical switch for prepare mpi datatype only once
  logical, private :: first=.true.

  procedure(thermal_conduction), pointer   :: phys_thermal_conduction => null()
  procedure(get_heatconduct), pointer   :: phys_get_heatconduct => null()
  procedure(getdt_heatconduct), pointer :: phys_getdt_heatconduct => null()

  abstract interface
    subroutine thermal_conduction()
    ! Meyer 2012 MNRAS 422,2102
      use mod_global_parameters
      use mod_ghostcells_update
    end subroutine thermal_conduction

    subroutine get_heatconduct(tmp,tmp1,tmp2,ixI^L,ixO^L,w,x)
      use mod_global_parameters
      
      integer, intent(in) :: ixI^L, ixO^L
      double precision, intent(in) ::  x(ixI^S,1:ndim), w(ixI^S,1:nw)
      !! tmp store the heat conduction energy changing rate
      double precision, intent(out) :: tmp(ixI^S),tmp1(ixI^S),tmp2(ixI^S)
    end subroutine get_heatconduct

    subroutine getdt_heatconduct(w,ixG^L,ix^L,dtnew,dx^D,x)
      use mod_global_parameters
      
      integer, intent(in) :: ixG^L, ix^L
      double precision, intent(in) :: dx^D, x(ixG^S,1:ndim)
      ! note that depending on strictsmall etc, w values may change 
      ! through call to getpthermal
      double precision, intent(inout) :: w(ixG^S,1:nw), dtnew
    end subroutine getdt_heatconduct
  end interface

contains
  !> Read this module"s parameters from a file
  subroutine tc_params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /tc_list/ tc_perpendicular, tc_saturate, tc_dtpar

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, tc_list, end=111)
111    close(unitpar)
    end do

  end subroutine tc_params_read

  !> Initialize the module
  subroutine thermal_conduction_init(phys_gamma)
    use mod_global_parameters
    use mod_physics

    double precision, intent(in) :: phys_gamma

    integer :: nwx,idir

    tc_gamma=phys_gamma

    tc_dtpar=tc_dtpar/dble(ndim)

    call tc_params_read(par_files)
  
    phys_thermal_conduction => do_thermal_conduction
    if(physics_type=='hd') then
      phys_get_heatconduct   => hd_get_heatconduct
      phys_getdt_heatconduct => hd_getdt_heatconduct
    else if(physics_type=='mhd') then
      phys_get_heatconduct   => mhd_get_heatconduct
      phys_getdt_heatconduct => mhd_getdt_heatconduct
    end if

    ! Determine flux variables
    nwx = 1                  ! rho (density)

    allocate(mom(ndir))
    do idir = 1, ndir
       nwx    = nwx + 1
       mom(idir) = nwx       ! momentum density
    end do

    nwx = nwx + 1
    e_     = nwx          ! energy density

    allocate(mag(ndir))
    do idir = 1, ndir
       nwx    = nwx + 1
       mag(idir) = nwx       ! magnetic field
    end do

    minp   = max(0.0d0, small_pressure)
    minrho = max(0.0d0, small_density)
    smalle = minp/(tc_gamma - 1.0d0)


    if(SI_unit) then
      ! Spitzer thermal conductivity with SI units
      kappa=8.d-12*unit_temperature**3.5d0/unit_length/unit_density/unit_velocity**3 
      ! thermal conductivity perpendicular to magnetic field
      kappe=4.d-30*unit_numberdensity**2/unit_magneticfield**2/unit_temperature**3*kappa
    else
      ! Spitzer thermal conductivity with cgs units
      kappa=8.d-7*unit_temperature**3.5d0/unit_length/unit_density/unit_velocity**3 
      ! thermal conductivity perpendicular to magnetic field
      kappe=4.d-26*unit_numberdensity**2/unit_magneticfield**2/unit_temperature**3*kappa
    end if

  end subroutine thermal_conduction_init

  subroutine do_thermal_conduction()
  ! Meyer 2012 MNRAS 422,2102
    use mod_global_parameters
    use mod_ghostcells_update
    
    double precision :: omega1,cmu,cmut,cnu,cnut
    double precision, allocatable :: bj(:)
    integer:: iigrid, igrid,j
    logical :: evenstep

    ! point bc mpi datatype to partial type for thermalconduction
    type_send_srl=>type_send_srl_p
    type_recv_srl=>type_recv_srl_p
    type_send_r=>type_send_r_p
    type_recv_r=>type_recv_r_p
    type_send_p=>type_send_p_p
    type_recv_p=>type_recv_p_p 
    ! create bc mpi datatype for ghostcells update
    if(first) then
      call create_bc_mpi_datatype(e_-1,1)
      first=.false.
    end if
    
    do iigrid=1,igridstail; igrid=igrids(iigrid);
      if(.not. allocated(pw(igrid)%w2)) allocate(pw(igrid)%w2(ixG^T,1:nw))
      if(.not. allocated(pw(igrid)%w3)) allocate(pw(igrid)%w3(ixG^T,1:nw))
      pw(igrid)%w1=pw(igrid)%w
      pw(igrid)%w2=pw(igrid)%w
      pw(igrid)%w3=pw(igrid)%w
    end do
    
    allocate(bj(0:s))
    bj(0)=1.d0/3.d0
    bj(1)=bj(0)
    if(s>1) then
      omega1=4.d0/dble(s**2+s-2)
      cmut=omega1/3.d0
    else
      omega1=0.d0
      cmut=1.d0
    endif
    
    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
      block=>pw(igrid)
      typelimiter=limiter(node(plevel_,igrid))
      typegradlimiter=gradient_limiter(node(plevel_,igrid))
      ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
      call evolve_step1(cmut,dt_tc,ixG^LL,ixM^LL,pw(igrid)%w1,pw(igrid)%w,&
                        pw(igrid)%x,pw(igrid)%w3)
      pw(igrid)%wb=>pw(igrid)%w1
    end do
    !$OMP END PARALLEL DO
    bcphys=.false.
    call getbc(global_time,0.d0,e_-1,1)
    if(s==1) then
      do iigrid=1,igridstail; igrid=igrids(iigrid);
        pw(igrid)%w(ixG^T,e_)=pw(igrid)%w1(ixG^T,e_)
        pw(igrid)%wb=>pw(igrid)%w
      end do
      ! point bc mpi data type back to full type for (M)HD
      type_send_srl=>type_send_srl_f
      type_recv_srl=>type_recv_srl_f
      type_send_r=>type_send_r_f
      type_recv_r=>type_recv_r_f
      type_send_p=>type_send_p_f
      type_recv_p=>type_recv_p_f
      bcphys=.true.
      deallocate(bj)
      return
    endif
    evenstep=.true.
    do j=2,s
      bj(j)=dble(j**2+j-2)/dble(2*j*(j+1))
      cmu=dble(2*j-1)/dble(j)*bj(j)/bj(j-1)
      cmut=omega1*cmu
      cnu=dble(1-j)/dble(j)*bj(j)/bj(j-2)
      cnut=(bj(j-1)-1.d0)*cmut
      if(evenstep) then
    !$OMP PARALLEL DO PRIVATE(igrid)
        do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
          block=>pw(igrid)
          typelimiter=limiter(node(plevel_,igrid))
          typegradlimiter=gradient_limiter(node(plevel_,igrid))
          ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
          call evolve_stepj(cmu,cmut,cnu,cnut,dt_tc,ixG^LL,ixM^LL,pw(igrid)%w1,&
                            pw(igrid)%w2,pw(igrid)%w,pw(igrid)%x,pw(igrid)%w3)
          pw(igrid)%wb=>pw(igrid)%w2
        end do
    !$OMP END PARALLEL DO
        call getbc(global_time,0.d0,e_-1,1)
        evenstep=.false.
      else
    !$OMP PARALLEL DO PRIVATE(igrid)
        do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
          block=>pw(igrid)
          typelimiter=limiter(node(plevel_,igrid))
          typegradlimiter=gradient_limiter(node(plevel_,igrid))
          ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
          call evolve_stepj(cmu,cmut,cnu,cnut,dt_tc,ixG^LL,ixM^LL,pw(igrid)%w2,&
                            pw(igrid)%w1,pw(igrid)%w,pw(igrid)%x,pw(igrid)%w3)
          pw(igrid)%wb=>pw(igrid)%w1
        end do
    !$OMP END PARALLEL DO
        call getbc(global_time,0.d0,e_-1,1)
        evenstep=.true.
      end if 
    end do
    if(evenstep) then
      do iigrid=1,igridstail; igrid=igrids(iigrid);
        pw(igrid)%w(ixG^T,e_)=pw(igrid)%w1(ixG^T,e_)
        pw(igrid)%wb=>pw(igrid)%w
      end do 
    else
      do iigrid=1,igridstail; igrid=igrids(iigrid);
        pw(igrid)%w(ixG^T,e_)=pw(igrid)%w2(ixG^T,e_)
        pw(igrid)%wb=>pw(igrid)%w
      end do 
    end if
    deallocate(bj)
    ! point bc mpi data type back to full type for (M)HD
    type_send_srl=>type_send_srl_f
    type_recv_srl=>type_recv_srl_f
    type_send_r=>type_send_r_f
    type_recv_r=>type_recv_r_f
    type_send_p=>type_send_p_f
    type_recv_p=>type_recv_p_f
    bcphys=.true.
  
  end subroutine do_thermal_conduction

  subroutine evolve_stepj(qcmu,qcmut,qcnu,qcnut,qdt,ixI^L,ixO^L,w1,w2,w,x,wold)
    
    use mod_global_parameters
    
    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: qcmu,qcmut,qcnu,qcnut,qdt
    double precision, intent(in) :: w1(ixI^S,1:nw),w(ixI^S,1:nw),wold(ixI^S,1:nw)
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w2(ixI^S,1:nw)
    
    double precision :: tmp(ixI^S),tmp1(ixI^S),tmp2(ixI^S)

    call phys_get_heatconduct(tmp,tmp1,tmp2,ixI^L,ixO^L,w1,x)

    w2(ixO^S,e_)=qcmu*w1(ixO^S,e_)+qcnu*w2(ixO^S,e_)+(1.d0-qcmu-qcnu)*w(ixO^S,e_)&
                +qcmut*qdt*tmp(ixO^S)+qcnut*wold(ixO^S,e_)
    
  end subroutine evolve_stepj

  subroutine evolve_step1(qcmut,qdt,ixI^L,ixO^L,w1,w,x,wold)
    use mod_global_parameters
    
    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: qcmut, qdt, w(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(out) ::w1(ixI^S,1:nw),wold(ixI^S,1:nw)
    
    double precision :: tmp(ixI^S),tmp1(ixI^S),tmp2(ixI^S),Te(ixI^S)
    integer :: lowindex(ndim), ix^D

    call phys_get_heatconduct(tmp,tmp1,tmp2,ixI^L,ixO^L,w,x)

    wold(ixO^S,e_)=qdt*tmp(ixO^S)
    ! update internal energy
    tmp1(ixO^S) = tmp1(ixO^S) + qcmut*wold(ixO^S,e_)
    
    ! ensure you never trigger negative pressure 
    ! hence code up energy change with respect to kinetic and magnetic
    ! part(nonthermal)
    if(small_temperature>0.d0) then
      Te(ixO^S)=tmp1(ixO^S)*(tc_gamma-1.d0)/w(ixO^S,rho_)
    endif
    if(strictsmall) then
      if(small_temperature>0.d0 .and. any(Te(ixO^S)<small_temperature)) then
        lowindex=minloc(Te(ixO^S))
        ^D&lowindex(^D)=lowindex(^D)+ixOmin^D-1;
        write(*,*)'too small temperature = ',minval(Te(ixO^S)),'at x=',&
       x(^D&lowindex(^D),1:ndim),lowindex,' with limit=',small_temperature,&
       ' on time=',global_time,' step=',it, 'where w(1:nwflux)=',w(^D&lowindex(^D),1:nwflux)
        call mpistop("==evolve_step1: too small temperature==")
      end if
      if(any(tmp1(ixO^S)<smalle)) then
        lowindex=minloc(tmp1(ixO^S))
        ^D&lowindex(^D)=lowindex(^D)+ixOmin^D-1;
        write(*,*)'too small internal energy = ',minval(tmp1(ixO^S)),'at x=',&
       x(^D&lowindex(^D),1:ndim),lowindex,' with limit=',smalle,&
       ' on time=',global_time,' step=',it, 'where w(1:nwflux)=',w(^D&lowindex(^D),1:nwflux)
        call mpistop("==evolve_step1: too small internal energy==")
      end if
      w1(ixO^S,e_) = tmp2(ixO^S)+tmp1(ixO^S)
    else
     {do ix^DB=ixOmin^DB,ixOmax^DB\}
        if(small_temperature>0.d0) then
          if(Te(ix^D)<small_temperature) then
            w1(ix^D,e_) = tmp2(ix^D)+small_temperature*w(ix^D,rho_)/(tc_gamma-1.d0)
          else
            w1(ix^D,e_) = tmp2(ix^D)+tmp1(ix^D)
          end if
        else
          if(tmp1(ix^D)<smalle) then
            w1(ix^D,e_) = tmp2(ix^D)+smalle
          else
            w1(ix^D,e_) = tmp2(ix^D)+tmp1(ix^D)
          end if
        end if
     {end do\}
    end if
  
  end subroutine evolve_step1

  subroutine mhd_get_heatconduct(tmp,tmp1,tmp2,ixI^L,ixO^L,w,x)
    use mod_global_parameters
    
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) ::  x(ixI^S,1:ndim), w(ixI^S,1:nw)
    !! tmp store the heat conduction energy changing rate
    double precision, intent(out) :: tmp(ixI^S),tmp1(ixI^S),tmp2(ixI^S)
    
    double precision, dimension(ixI^S,1:ndir) :: mf,qvec
    double precision, dimension(ixI^S,1:ndim) :: gradT,qvecsat
    double precision, dimension(:^D&,:), allocatable :: qvec_per, qvec_max
    double precision, dimension(ixI^S) :: B2inv, BgradT,cs3,qflux,qsatflux
    integer, dimension(ndim) :: lowindex
    integer :: ix^L,idims,idir,ix^D
    logical :: Bnull(ixI^S)

    ix^L=ixO^L^LADD2;
    if (ixI^L^LTix^L|.or.|.or.) &
       call mpistop("Need two extra layers for thermal conduction")
    
    ! tmp2 store kinetic+magnetic energy before addition of heat conduction source
    tmp2(ixI^S) = 0.5d0 * (sum(w(ixI^S,mom(:))**2,dim=ndim+1)/w(ixI^S,rho_) + &
         sum(w(ixI^S,mag(:))**2,dim=ndim+1))

    ! tmp1 store internal energy
    tmp1(ixI^S)=w(ixI^S,e_)-tmp2(ixI^S)

    ! Clip off negative pressure if small_pressure is set
    if(strictsmall) then
       if (any(tmp1(ixI^S)<smalle)) then
         lowindex=minloc(tmp1(ixI^S))
         ^D&lowindex(^D)=lowindex(^D)+ixImin^D-1;
         write(*,*)'too low internal energy = ',minval(tmp1(ixI^S)),' at x=',&
         x(^D&lowindex(^D),1:ndim),lowindex,' with limit=',smalle,' on time=',global_time
       end if
    else
    {do ix^DB=ixImin^DB,ixImax^DB\}
       if(tmp1(ix^D)<smalle) then
          tmp1(ix^D)=smalle
       end if
    {end do\}
    end if
    ! compute the temperature
    tmp(ixI^S)=tmp1(ixI^S)*(tc_gamma-one)/w(ixI^S,rho_)
    if(small_temperature>0.d0) then
      if(strictsmall) then
         if(any(tmp(ixI^S)<small_temperature)) then
           lowindex=minloc(tmp(ixI^S))
           ^D&lowindex(^D)=lowindex(^D)+ixImin^D-1;
           write(*,*)'too low temperature = ',minval(tmp(ixI^S)),' at x=',&
           x(^D&lowindex(^D),1:ndim),lowindex,' with limit=',small_temperature,' on time=',global_time
         end if
      else
      {do ix^DB=ixImin^DB,ixImax^DB\}
         if(tmp(ix^D)<small_temperature) then
            tmp(ix^D)=small_temperature
         end if
      {end do\}
      end if
    end if
    do idims=1,ndim
      ! idirth component of gradient of temperature at cell center
      select case(typegrad)
      case("central")
        ix^L=ixI^L^LSUB1;
        call gradient(tmp,ixI^L,ix^L,idims,B2inv)
      case("limited")
        ix^L=ixI^L^LSUB2;
        call gradientS(tmp,ixI^L,ix^L,idims,B2inv)
      end select
      ! get grad T in all directions
      gradT(ix^S,idims)=B2inv(ix^S)
    end do
    ! B
    if(B0field) then
      mf(ix^S,1:ndir)=w(ix^S,mag(1):mag(ndir))+block%B0(ix^S,1:ndir,0)
    else
      mf(ix^S,1:ndir)=w(ix^S,mag(1):mag(ndir));
    end if
    ! B^-2
    B2inv=0.d0
    do idir=1,ndir
      B2inv(ix^S)=B2inv(ix^S)+mf(ix^S,idir)**2
    end do
    Bnull(ixI^S)=.false.
    where(B2inv(ix^S)/=0.d0)
      B2inv(ix^S)=1.d0/B2inv(ix^S)
    elsewhere
      Bnull(ix^S)=.true.
    end where
    
    BgradT(ix^S)=(^D&mf(ix^S,^D)*gradT(ix^S,^D)+)*B2inv(ix^S)
    cs3(ix^S)=kappa*dsqrt(tmp(ix^S))*tmp(ix^S)**2*BgradT(ix^S)
    
    do idims=1,ndim
      qvec(ix^S,idims)=mf(ix^S,idims)*cs3(ix^S)
    end do

    if(tc_saturate) then
      ! consider saturation (Cowie and Mckee 1977 ApJ, 211, 135)
      cs3(ix^S)=dsqrt(tmp(ix^S)**3)
      ! unsigned saturated TC flux = 5 phi rho c**3
      qsatflux(ix^S)=5.d0*w(ix^S,rho_)*cs3(ix^S)
      ! strength of classic TC flux
      qflux(ix^S)=dsqrt(^D&qvec(ix^S,^D)**2+)
      ! sign(b * Grad Te) 5 phi rho c**3 Bi/B 
      {do ix^DB=ixmin^DB,ixmax^DB\}
        if(.false. .and. qflux(ix^D)>qsatflux(ix^D) .and. idims==1) write(*,*)& 
          'it',it,' ratio=',qflux(ix^D)/qsatflux(ix^D),' TC saturated at ',&
          x(ix^D,:),' rho',w(ix^D,rho_),' Te',tmp(ix^D)
        if(qflux(ix^D)>qsatflux(ix^D)) then
        ! saturated TC flux = sign(b * Grad Te) 5 phi rho c**3
          qsatflux(ix^D)=sign(1.d0,BgradT(ix^D))*qsatflux(ix^D)*dsqrt(B2inv(ix^D))
          do idims=1,ndim
            qvec(ix^D,idims)=qsatflux(ix^D)*mf(ix^D,idims)
          end do
        end if
      {end do\}
    end if
    
    if(tc_perpendicular) then
    ! consider thermal conduction perpendicular to magnetic field
    ! van der Linden and Goossens 1991 SoPh 131, 79; Orlando et al 2008 ApJ 678, 274
      allocate(qvec_per(ixI^S,1:ndim), qvec_max(ixI^S,1:ndim))
      do idims=1,ndim
        ! q_per = kappe n^2 B^-2 Te^-0.5 (Grad Te - (e_b . Grad Te )e_b) 
        qvec_per(ix^S,idims)=kappe*w(ix^S,rho_)**2*B2inv(ix^S)/dsqrt(tmp(ix^S))&
                             *(gradT(ix^S,idims)-BgradT(ix^S)*mf(ix^S,idims))
      end do
      ! maximal thermal conducting (saturated) flux
      do idims=1,ndim
        qvec_max(ix^S,idims)=kappa*dsqrt(tmp(ix^S)**5)*gradT(ix^S,idims)
      end do

      if(tc_saturate) then
        ! consider saturation (Cowie and Mckee 1977 ApJ, 211, 135)
        qsatflux(ix^S)=5.d0*w(ix^S,rho_)*cs3(ix^S)
        qflux(ix^S)=dsqrt(^D&qvec_max(ix^S,^D)**2+)
        {do ix^DB=ixmin^DB,ixmax^DB\}
          if(.false. .and. qflux(ix^D)>qsatflux(ix^D) .and. idims==1) write(*,*) & 
            'it',it,' ratio=',qflux(ix^D)/qsatflux(ix^D),' TC_PER saturated at ',&
            x(ix^D,:),' rho',w(ix^D,rho_),' Te',tmp(ix^D)
          if(qflux(ix^D)>qsatflux(ix^D)) then
            qsatflux(ix^D)=qsatflux(ix^D)/qflux(ix^D)
            do idims=1,ndim
              qvec_max(ix^D,idims)=qsatflux(ix^D)*qvec_max(ix^D,idims)
            end do
          end if
        {end do\}
      end if

      ! maximal thermal conducting flux perpendicular to magnetic field
      qvec_max(ix^S,1:ndim)=qvec_max(ix^S,1:ndim)-qvec(ix^S,1:ndim)
      qsatflux(ix^S)= ^D&qvec_max(ix^S,^D)**2+
      qflux(ix^S)= ^D&qvec_per(ix^S,^D)**2+
      {do ix^DB=ixmin^DB,ixmax^DB\}
         if(qflux(ix^D)>qsatflux(ix^D)) then
           qvec(ix^D,1:ndim)=qvec(ix^D,1:ndim)+qvec_max(ix^D,1:ndim)
         else
           qvec(ix^D,1:ndim)=qvec(ix^D,1:ndim)+qvec_per(ix^D,1:ndim)
         end if   
      {end do\}
    end if
    
    {do ix^DB=ixmin^DB,ixmax^DB\}
      if(Bnull(ix^D)) then
        ! consider magnetic null point
        qvec(ix^D,1:ndim)=kappa*dsqrt(tmp(ix^D)**5)*gradT(ix^D,1:ndim)
        if(tc_saturate) then
          ! unsigned saturated TC flux = 5 phi rho c**3
          qsatflux(ix^D)=5.d0*w(ix^D,rho_)*cs3(ix^D)
          qflux(ix^D)=dsqrt(sum(qvec(ix^D,1:ndim)**2))
          if(qflux(ix^D)>qsatflux(ix^D)) then
            qsatflux(ix^D)=qsatflux(ix^D)/qflux(ix^D)
            do idims=1,ndim
              qvec(ix^D,idims)=qsatflux(ix^D)*qvec(ix^D,idims)
            end do
          end if
        end if
      end if  
    {end do\}
    
    ! store thermal conduction source term in tmp
    select case(typediv)
       case("central")
          call divvector(qvec,ixI^L,ixO^L,tmp)
       case("limited")
          call divvectorS(qvec,ixI^L,ixO^L,tmp)
    end select
  
  end subroutine mhd_get_heatconduct

  subroutine mhd_getdt_heatconduct(w,ixI^L,ixO^L,dtnew,dx^D,x)
    !Check diffusion time limit dt < tc_dtpar*dx_i**2/((gamma-1)*kappa_i/rho)
    !where                      kappa_i=kappa*B_i**2/B**2
    !and                        T=p/rho
    use mod_global_parameters
    use mod_physics
    
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: dx^D, x(ixI^S,1:ndim)
    ! note that depending on strictsmall etc, w values may change 
    ! through call to getpthermal
    double precision, intent(inout) :: w(ixI^S,1:nw), dtnew
    
    double precision :: dxinv(1:ndim),mf(ixI^S,1:ndir)
    double precision :: tmp2(ixI^S),tmp(ixI^S),Te(ixI^S),B2inv(ixI^S)
    double precision :: dtdiff_tcond, dtdiff_tsat
    integer          :: idim,idir,ix^D
    integer, dimension(ndim)       :: lowindex

    ^D&dxinv(^D)=one/dx^D;
    
    call phys_get_pthermal(w,x,ixI^L,ixO^L,tmp)
    ! Clip off negative pressure if small_pressure is set
    if(strictsmall) then
      if(any(tmp(ixO^S)<minp)) then
        lowindex=minloc(tmp(ixO^S))
        ^D&lowindex(^D)=lowindex(^D)+ixOmin^D-1;
        write(*,*)'low pressure = ',minval(tmp(ixO^S)),' at x=',&
        x(^D&lowindex(^D),1:ndim),lowindex,' with limit=',minp,' on time=',global_time,&
        ' step=',it
       call mpistop("=== strictsmall in getdt_heatconduct_mhd: low pressure ===")
      end if
    else
    {do ix^DB=ixOmin^DB,ixOmax^DB\}
       if(tmp(ix^D)<minp) then
          tmp(ix^D)=minp
       end if
    {end do\}
    end if
    !temperature
    Te(ixO^S)=tmp(ixO^S)/w(ixO^S,rho_)
    !kappa_i
    tmp(ixO^S)=kappa*dsqrt(Te(ixO^S)**5)
    !(gamma-1)*kappa_i/rho
    tmp(ixO^S)=(tc_gamma-one)*tmp(ixO^S)/w(ixO^S,rho_)
    
    ! B
    if(B0field) then
      mf(ixO^S,1:ndir)=w(ixO^S,mag(1):mag(ndir))+block%B0(ixO^S,1:ndir,0)
    else
      mf(ixO^S,1:ndir)=w(ixO^S,mag(1):mag(ndir))
    end if
    ! B^-2
    B2inv=0.d0
    do idir=1,ndir
      B2inv(ixO^S)=B2inv(ixO^S)+mf(ixO^S,idir)**2
    end do
    where(B2inv(ixO^S)/=0.d0)
      B2inv(ixO^S)=1.d0/B2inv(ixO^S)
    end where

    do idim=1,ndim
       ! B_i**2/B**2
       where(B2inv(ixO^S)/=0.d0)
         tmp2(ixO^S)=mf(ixO^S,idim)**2*B2inv(ixO^S)
       elsewhere
         tmp2(ixO^S)=1.d0
       end where
       ! dt< tc_dtpar * dx_idim**2/((gamma-1)*kappa_i/rho*B_i**2/B**2)
       dtdiff_tcond=tc_dtpar/maxval(tmp(ixO^S)*tmp2(ixO^S)*dxinv(idim)**2)
       if(tc_saturate) then
         ! dt< tc_dtpar* dx_idim**2/((gamma-1)*sqrt(Te)*5*phi)
         ! with an empirical coefficient dx_idim
         dtdiff_tsat=tc_dtpar/maxval((tc_gamma-1.d0)*dsqrt(Te(ixO^S))*&
                     5.d0*dxinv(idim)**2)
         ! choose the slower flux (bigger time step) between classic and saturated
         dtdiff_tcond=max(dtdiff_tcond,dtdiff_tsat)
       end if
       ! limit the time step
       dtnew=min(dtnew,dtdiff_tcond)
    end do
    dtnew=dtnew/dble(ndim)
  
  end subroutine mhd_getdt_heatconduct

  subroutine hd_get_heatconduct(tmp,tmp1,tmp2,ixI^L,ixO^L,w,x)
    use mod_global_parameters
    
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) ::  x(ixI^S,1:ndim), w(ixI^S,1:nw)
    !! tmp store the heat conduction energy changing rate
    double precision, intent(out) :: tmp(ixI^S),tmp1(ixI^S),tmp2(ixI^S)
    double precision :: qvec(ixI^S,1:ndir),Te(ixI^S),qflux(ixI^S),qsatflux(ixI^S)
    integer:: ix^L,idims,ix^D
    integer, dimension(ndim)       :: lowindex
    
    ix^L=ixO^L^LADD2;
    if (ixI^L^LTix^L|.or.|.or.) &
       call mpistop("Need two extra layers for thermal conduction")
    
    ! store old kinetic energy
    tmp2(ixI^S)=half*sum(w(ixI^S,mom(:))**2,dim=ndim+1)/w(ixI^S,rho_)
    ! Calculate pressure=(gamma-1)*(e-0.5*(2ek+2eb))
    tmp(ixI^S)=(tc_gamma-one)*(w(ixI^S,e_)-tmp2(ixI^S))
    ! Clip off negative pressure if small_pressure is set
    if(strictsmall) then
      if(any(tmp(ixI^S)<minp)) then
        lowindex=minloc(tmp(ixI^S))
        write(*,*)'low pressure = ',minval(tmp(ixI^S)),' at x=',&
        x(^D&lowindex(^D),1:ndim),lowindex,' with limit=',minp,' on time=',global_time
        call mpistop("=== strictsmall in heatconduct: low pressure ===")
      end if
    else
       {do ix^DB=ixImin^DB,ixImax^DB\}
         if(tmp(ix^D)<minp) then
          tmp(ix^D)=minp
         end if
       {end do\}
    end if
    ! store old internal energy
    tmp1(ixI^S)=tmp(ixI^S)/(tc_gamma-one)
    
    ! compute temperature before source addition
    Te(ixI^S)=tmp(ixI^S)/w(ixI^S,rho_)
    if(small_temperature>0.d0) then
      if(strictsmall) then
         if(any(Te(ixI^S)<small_temperature)) then
           lowindex=minloc(Te(ixI^S))
           ^D&lowindex(^D)=lowindex(^D)+ixImin^D-1;
           write(*,*)'too low temperature = ',minval(Te(ixI^S)),' at x=',&
           x(^D&lowindex(^D),1:ndim),lowindex,' with limit=',small_temperature,' on time=',global_time,&
           ' step=',it,' where density=',w(^D&lowindex(^D),rho_),' velocity=',&
           dsqrt(sum(w(^D&lowindex(^D),mom(:))**2,dim=1)/w(^D&lowindex(^D),rho_)**2)
           call mpistop("=== strictsmall in heatcond_hd: low temperature ===")
         end if
      else
      {do ix^DB=ixImin^DB,ixImax^DB\}
         if(Te(ix^D)<small_temperature) then
            Te(ix^D)=small_temperature
         end if
      {end do\}
      end if
    end if
    ! compute grad T and store grad T vector
    do idims=1,ndim
      ! idirth component of gradient of temperature at cell center
       select case(typegrad)
       case("central")
         ix^L=ixI^L^LSUB1;
         call gradient(Te,ixI^L,ix^L,idims,tmp)
       case("limited")
         ix^L=ixI^L^LSUB2;
         call gradientS(Te,ixI^L,ix^L,idims,tmp)
       end select
       qvec(ix^S,idims)=tmp(ix^S)*kappa*dsqrt(Te(ix^S)**5)
    end do

    if(tc_saturate) then
      ! consider saturation with unsigned saturated TC flux = 5 phi rho c**3
      qsatflux(ix^S)=5.d0*w(ix^S,rho_)*dsqrt(Te(ix^S)**3)
      qflux(ix^S)=dsqrt(^D&qvec(ix^S,^D)**2+)
      {do ix^DB=ixmin^DB,ixmax^DB\}
        if(qflux(ix^D)>qsatflux(ix^D)) then
          qsatflux(ix^D)=qsatflux(ix^D)/qflux(ix^D)
          do idims=1,ndim
            qvec(ix^D,idims)=qsatflux(ix^D)*qvec(ix^D,idims)
          end do
        end if
      {end do\}
    end if

    select case(typediv)
      case("central")
        call divvector(qvec,ixI^L,ixO^L,tmp)
      case("limited")
        call divvectorS(qvec,ixI^L,ixO^L,tmp)
    end select

  end subroutine hd_get_heatconduct

  subroutine hd_getdt_heatconduct(w,ixI^L,ixO^L,dtnew,dx^D,x)
    ! Check diffusion time limit dt < tc_dtpar * dx_i**2 / ((gamma-1)*kappa_i/rho)
    use mod_global_parameters
    use mod_physics
    
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: dx^D, x(ixI^S,1:ndim)
    ! note that depending on strictsmall etc, w values may change 
    ! through call to getpthermal
    double precision, intent(inout) :: w(ixI^S,1:nw), dtnew
    
    double precision :: dxinv(1:ndim), tmp(ixI^S), Te(ixI^S)
    double precision :: dtdiff_tcond,dtdiff_tsat
    integer          :: idim,ix^D
    integer, dimension(ndim)       :: lowindex

    ^D&dxinv(^D)=one/dx^D;

    call phys_get_pthermal(w,x,ixI^L,ixO^L,tmp)
    ! Clip off negative pressure if small_pressure is set
    if(strictsmall) then
      if(any(tmp(ixO^S)<minp)) then
        lowindex=minloc(tmp(ixO^S))
        ^D&lowindex(^D)=lowindex(^D)+ixOmin^D-1;
        write(*,*)'low pressure = ',minval(tmp(ixO^S)),' at x=',&
        x(^D&lowindex(^D),1:ndim),lowindex,' with limit=',minp,' on time=',global_time
        call mpistop("=== strictsmall in getdt_heatconduct_hd: low pressure ===")
      end if
    else
      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        if(tmp(ix^D)<minp) then
         tmp(ix^D)=minp
        end if
      {end do\}
    end if
    Te(ixO^S)=tmp(ixO^S)/w(ixO^S,rho_)
    tmp(ixO^S)=(tc_gamma-one)*kappa*dsqrt((Te(ixO^S))**5)/w(ixO^S,rho_)
    
    do idim=1,ndim
       ! dt< tc_dtpar * dx_idim**2/((gamma-1)*kappa_idim/rho)
       dtdiff_tcond=tc_dtpar/maxval(tmp(ixO^S)*dxinv(idim)**2)
       if(tc_saturate) then
         ! dt< tc_dtpar* dx_idim**2/((gamma-1)*sqrt(Te)*5*phi)
         dtdiff_tsat=tc_dtpar/maxval((tc_gamma-1.d0)*dsqrt(Te(ixO^S))*&
                     5.d0*dxinv(idim)**2)
         ! choose the slower flux (bigger time scale) between classic and saturated
         dtdiff_tcond=max(dtdiff_tcond,dtdiff_tsat)
       end if
       ! limit the time step
       dtnew=min(dtnew,dtdiff_tcond)
    end do
    dtnew=dtnew/dble(ndim)
  
  end subroutine hd_getdt_heatconduct

end module mod_thermal_conduction
