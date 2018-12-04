module mod_srhd_phys

  implicit none
  private

  !> The adiabatic index
  double precision, public :: srhd_gamma = 5.d0/3.0d0

  !> Whether ideal eos is used
  logical, public, protected              :: srhd_ideal = .false.

  !> Whether Mathews eos is used
  logical, public, protected              :: srhd_mathews = .false.

! Probably not necessary
  !> Index of the energy density (-1 if not present)
  integer, public, protected              :: e_

  !> Index of the density (lab frame)
  integer, public, protected              :: d_

  !> Index of the  density (primitive?)
  integer, public, protected              :: s0_

!  !> Index of the density
!  integer, public, protected              :: rho_

  !> Indices of the momentum density (this is the old s1,s2,s3)
  integer, allocatable, public, protected :: rmom(:)

  !> Index of the total energy density 
  integer, public, protected              :: tau_

  !> Index of the pressure
  integer, public, protected              :: p_

  !> Number of tracer species
  integer, public, protected              :: hd_n_tracer = 0

  !> Indices of the tracers
  integer, allocatable, public, protected :: tracer(:)

  !Auxiliary variables
  !> Index of Lorentz factor
  integer, public, protected :: lfac_
  
  !> Index of pressure
  integer, public, protected :: p_

  integer, parameter :: nvector=1                             ! No. vector vars
  integer, dimension(nvector), parameter :: iw_vector=(/ s0_ /)

  integer, parameter:: gamma_=1,neqpar=1                     ! equation parameters

  double precision :: smalltau,smallxi,minrho,minp

! Public methods ??
!  public :: srhd_phys_init
!  public :: srhd_kin_en
!  public :: srhd_get_pthermal
!  public :: srhd_to_conserved
!  public :: srhd_to_primitive

contains

!---------------------------------------------------------------------
  subroutine srhd_read_params(files)
    use mod_global_parameters
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /srhd_list/srhd_energy, srhd_n_tracer,srhd_gamma
    
    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, hd_list, end=111)
111    close(unitpar)
    end do

  end subroutine srhd_read_params
!---------------------------------------------------------------------
!
!add the srhd write info subroutine ?
!
!----------------------------------------- ---------------------------

!> Initialize the module
  subroutine srhd_phys_init()
    use mod_global_parameters
    use mod_physics

    integer :: itr, idir

    call srhd_read_params(par_files)

    physics_type = "srhd"
    phys_energy  = srhd_energy

    ! Determine flux variables
    rho_ = var_set_rho()
    d_ = rho_

    allocate(rmom(ndir))
    rmom(:) = var_set_momentum(ndir)

    ! Set index of energy variable
    e_ = var_set_energy()
!   related to the choice of eos ??
!    if (hd_energy) then
!       e_ = var_set_energy()
!       p_ = e_
!    else
!       e_ = -1
!       p_ = -1
!    end if

    lfac_ = var_set_extravar("lfac", "lfac")
    p_ = var_set_extravar("pressure", "pressure")

    phys_get_dt              => srhd_get_dt
    phys_get_cmax            => srhd_get_cmax
    phys_get_cbounds         => srhd_get_cbounds
    phys_get_flux            => srhd_get_flux
    phys_add_source_geom     => srhd_add_source_geom
    phys_add_source          => srhd_add_source
    phys_to_conserved        => srhd_to_conserved
    phys_to_primitive        => srhd_to_primitive
    phys_check_params        => srhd_check_params
    phys_check_w             => srhd_check_w
    phys_get_pthermal        => srhd_get_pthermal
    phys_write_info          => srhd_write_info
    phys_handle_small_values => srhd_handle_small_values

    ! Whether diagonal ghost cells are required for the physics (TODO: true?)
    phys_req_diagonal = .true.

    ! derive units from basic units (TO BE ADDED FOR SRHD)
!    call srhd_physical_units()

    allocate(tracer(srhd_n_tracer))

    ! Set starting index of tracers
    do itr = 1, srhd_n_tracer
       tracer(itr) = var_set_fluxvar("trc", "trp", itr, need_bc=.false.)
    end do


    ! Check whether custom flux types have been defined
    if (.not. allocated(flux_type)) then
       allocate(flux_type(ndir, nw))
       flux_type = flux_default
    else if (any(shape(flux_type) /= [ndir, nw])) then
       call mpistop("phys_check error: flux_type has wrong shape")
    end if

    nvector      = 1 ! No. vector vars
    allocate(iw_vector(nvector))
    iw_vector(1) = rmom(1) - 1

  end subroutine srhd_phys_init
!==============================================================================
!  subroutine srhd_physical_units
!    use mod_global_parameters
!    double precision :: mp,kB
!    ! Derive scaling units
!    if(SI_unit) then
!      mp=mp_SI
!      kB=kB_SI
!    else
!      mp=mp_cgs
!      kB=kB_cgs
!    end if
!    if(unit_velocity==0) then
!      unit_density=(1.d0+4.d0*He_abundance)*mp*unit_numberdensity
!      unit_pressure=(2.d0+3.d0*He_abundance)*unit_numberdensity*kB*unit_temperature
!      unit_velocity=dsqrt(unit_pressure/unit_density)
!      unit_time=unit_length/unit_velocity
!    else
!      unit_density=(1.d0+4.d0*He_abundance)*mp*unit_numberdensity
!      unit_pressure=unit_density*unit_velocity**2
!      unit_temperature=unit_pressure/((2.d0+3.d0*He_abundance)*unit_numberdensity*kB)
!      unit_time=unit_length/unit_velocity
!    end if

!  end subroutine srhd_physical_units
!==============================================================================
  subroutine checkglobaldata

    use mod_global_parameters

    !-----------------------------------------------------------------------------
    minp   = max(zero,smallp)
    minrho = max(zero,smallrho)
    call smallvaluesEOS

  end subroutine checkglobaldata
!=============================================================================
  subroutine srhd_check_w(checkprimitive,ixI^L,ixO^L,w,flag)
    use mod_global_parameters

    logical, intent(in)          :: checkprimitive
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S,nw)
    logical, intent(out)         :: flag(ixG^T)
!    double precision             :: tmp(ixI^S)
    !-----------------------------------------------------------------------------

    flag(iO^S)= 0

    if (checkprimitive) then
       if(useprimitiveRel)then
          ! check   rho>=0, p>=smallp
          flag(ixO^S) = (w(ixO^S,rho_) >= minrho).and. &
               (w(ixO^S,pp_)  >= minp)
       else
          ! check  v^2 < 1, rho>=0, p>=smallp
          flag(ixO^S) = (sum(w(ixO^S,rmom(:))**2.0d0) < one).and. &
               (w(ixO^S,rho_) >= minrho).and. &
               (w(ixO^S,pp_)  >= minp)
       endif
    else
       ! Check D>=0 and lower limit for tau
       flag(ixO^S) = (w(ixO^S,d_)    >= minrho).and. &
            (w(ixO^S,tau_)  >= smalltau)
    endif

  end subroutine srhd_check_w
  !=============================================================================
  subroutine srhd_to_conserved(ixI^L,ixO^L,w,x,patchw)

    ! Transform primitive variables into conservative ones
    ! (rho,v,p) ---> (D,S,tau)
    ! call to smallvalues
    ! --> latter only used for correcting procedure in correctaux
    ! --> input array patchw for spatially selective transformation

    use mod_global_parameters

    integer, intent(in)               :: ixI^L, ixO^L
    double precision, intent(inout)   :: w(ixI^S,nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    logical, intent(in)               :: patchw(ixG^T)
    integer                         :: idir, itr

    double precision,dimension(ixG^T) :: xi
    !-----------------------------------------------------------------------------

    if(useprimitiveRel)then
       ! assumes four velocity computed in primitive (rho u p) with u=lfac*v
       xi(ixO^S)=one+(sum(w(ixO^S,uvel(:))**2.0d0))
       ! fill the auxiliary variable lfac_ (Lorentz Factor) and p_ (pressure)
       w(ixO^S,lfac_)=dsqrt(xi(ixO^S))
       w(ixO^S,p_)=w(ixO^S,pp_)
    else
       ! assumes velocity in primitive (rho v p) 
       xi(ixO^S)=one-(sum(w(ixO^S,vvel(:))**2.0d0))
       ! fill the auxiliary variable lfac_ (Lorentz Factor) and p_ (pressure)
       w(ixO^S,lfac_)=one/dsqrt(xi(ixO^S))
       w(ixO^S,p_)=w(ixO^S,pp_)
    endif

    call Enthalpy(w,ixI^L,ixO^L,patchw,xi)

    if(useprimitiveRel)then
       ! compute xi=Lfac w  (enthalphy w)
       xi(ixO^S)=w(ixO^S,lfac_)*xi(ixO^S)
    else
       ! compute xi=Lfac^2 w  (enthalphy w)
       xi(ixO^S)=w(ixO^S,lfac_)*w(ixO^S,lfac_)*xi(ixO^S)     
    endif

    if(useprimitiveRel)then
       w(ixO^S,d_)=w(ixO^S,rho_)*w(ixO^S,lfac_) 
       do idir = 1, ndir
       w(ixO^S, rmom(idir)) = xi(ixO^S)*w(ixO^S,uvel(idir))
       end do
!       ^C&w(ixO^S,s^C_)=xi(ixO^S)*w(ixO^S,u^C_);
       w(ixO^S,tau_)=xi(ixO^S)*w(ixO^S,lfac_) - w(ixO^S,p_) - w(ixO^S,d_)
    else
       w(ixO^S,d_)=w(ixO^S,rho_)*w(ixO^S,lfac_) 
       do idir = 1, ndir
       w(ixO^S, rmom(idir)) = xi(ixO^S)*w(ixO^S,vvel(idir))
       end do
!       ^C&w(ixO^S,s^C_)=xi(ixO^S)*w(ixO^S,v^C_);
       w(ixO^S,tau_)=xi(ixO^S) - w(ixO^S,p_) - w(ixO^S,d_)
    endif

!    if (srhd_n_tracer)
    ! We got D, now we can get the conserved tracers:
!    {w(ixO^S,tr^FL_) = w(ixO^S,d_)*w(ixO^S,tr^FL_)\}
!    end if

    if(check_small_values) call srhd_handle_smallvalues(.false.,w,x,ixI^L,ixO^L,"srhd_to_conserved")

  end subroutine srhd_to_conserved
  !=============================================================================
  subroutine srhd_to_primitive(ixI^L,ixO^L,w,x)

    ! Transform conservative variables into primitive ones
    ! (D,S,tau) ---> (rho,v,p)

    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S,nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    integer                         :: itr, idir
    !-----------------------------------------------------------------------------
    ! Calculate pressure and Lorentz factor from conservative variables only
    ! these are put in lfac_ and p_ auxiliaries
    call getaux(.true.,w,x,ixI^L,ixO^L,'primitive')
    ! note: on exit from getaux: gauranteed 
    !    xi=(d+tau+p)>smallp*gamma/(gamma-1), lfac>=1, p>smallp

    ! replace conservative with primitive variables
    ! compute velocity
    if(useprimitiveRel)then
       do idir=1, ndir
       w(ixO^S,uvel(idur))=w(ixO^S,lfac_)*w(ixO^S,rmom(idir))/(w(ixO^S,d_)+w(ixO^S,tau_)+w(ixO^S,p_));
       end do
    else
       do idir=1, ndir
       w(ixO^S,vvel(idir))=w(ixO^S,rmom(idir))/(w(ixO^S,d_)+w(ixO^S,tau_)+w(ixO^S,p_));
       end do
    endif

    ! compute density
    w(ixO^S,rho_)=w(ixO^S,d_)/w(ixO^S,lfac_)
    ! fill pressure
    w(ixO^S,pp_)=w(ixO^S,p_)

 !   {#IFDEF TRACER
 !   ! We got lor, rho, Dtr, now we can get the tracers:
 !   {^FL&w(ixO^S,tr^FL_) = w(ixO^S,Dtr^FL_)/w(ixO^S,lfac_)/w(ixO^S,rho_)\}
 !   }

!    if (tlow>zero) call fixp_usr(ixI^L,ixO^L,w,x)

  end subroutine srhd_to_primitive(ixI^L,ixO^L,w,x) 
  !=============================================================================
  subroutine e_to_rhos(ixI^L,ixO^L,w,x)

    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    double precision :: w(ixI^S,nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    !-----------------------------------------------------------------------------

 !   ADD CORRECT FLAG FOR EACH EOS
 !   if (hd_energy) then
 !      w(ixO^S, e_) = (hd_gamma - 1.0d0) * w(ixO^S, rho_)**(1.0d0 - hd_gamma) * &
 !           (w(ixO^S, e_) - hd_kin_en(w, ixI^L, ixO^L))
 !   else
 !      call mpistop("energy from entropy can not be used with -eos = iso !")
 !   end if

    call mpistop("e to rhos for SRHDEOS unavailable")

  end subroutine e_to_rhos
  !=============================================================================
  subroutine rhos_to_e(ixI^L,ixO^L,w,x)

    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    double precision :: w(ixI^S,nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    !-----------------------------------------------------------------------------

 !   ADD CORRECT FLAG FOR EACH EOS
!    if (hd_energy) then
!       w(ixO^S, e_) = w(ixO^S, rho_)**(hd_gamma - 1.0d0) * w(ixO^S, e_) &
!            / (hd_gamma - 1.0d0) + hd_kin_en(w, ixI^L, ixO^L)
!    else
!       call mpistop("entropy from energy can not be used with -eos = iso !")
!    end if

    call mpistop("rhos to e for SRHDEOS unavailable")

  end subroutine rhos_to_e
  !=============================================================================
 ! subroutine ppmflatcd(ixI^L,ixO^L,ixL^L,ixR^L,w,d2w,drho,dp)

 !   use mod_global_parameters

 !   integer, intent(in)           :: ixI^L,ixO^L,ixL^L,ixR^L
 !   double precision, intent(in)  :: w(ixI^S,nw),d2w(ixG^T,1:nwflux)

 !   double precision, intent(out) :: drho(ixG^T),dp(ixG^T)
    !-----------------------------------------------------------------------------

 !   if(useprimitive)then
 !      call Calcule_Geffect(w,ixI^L,ixO^L,.false.,drho)
 !      drho(ixO^S) =drho(ixO^S)*dabs(d2w(ixO^S,rho_))&
 !           /min(w(ixL^S,rho_),w(ixR^S,rho_))
 !      dp(ixO^S) = dabs(d2w(ixO^S,pp_))/min(w(ixL^S,pp_),w(ixR^S,pp_))
 !   end if

 ! end subroutine ppmflatcd
  !=============================================================================
 ! subroutine ppmflatsh(ixI^L,ixO^L,ixLL^L,ixL^L,ixR^L,ixRR^L,idims,w,drho,dp,dv)

 !   use mod_global_parameters

 !   integer, intent(in)           :: ixI^L,ixO^L,ixLL^L,ixL^L,ixR^L,ixRR^L
 !   integer, intent(in)           :: idims
 !   double precision, intent(in)  :: w(ixI^S,nw)

 !    double precision, intent(out) :: drho(ixG^T),dp(ixG^T),dv(ixG^T)

 !   double precision :: v(ixG^T)
    !-----------------------------------------------------------------------------

 !   if(useprimitive)then
       ! eq. B15, page 218, Mignone and Bodo 2005, ApJS (beta1)
 !      where (dabs(w(ixRR^S,pp_)-w(ixLL^S,pp_))>smalldouble)
 !         drho(ixO^S) = dabs((w(ixR^S,pp_)-w(ixL^S,pp_))&
 !              /(w(ixRR^S,pp_)-w(ixLL^S,pp_)))
 !      elsewhere
 !         drho(ixO^S) = zero
 !      end where

       !  eq. B76, page 48, Miller and Collela 2002, JCP 183, 26 
       !  use "dp" to save sound speed
 !      call getcsound2(w,ixI^L,ixO^L,.false.,dv,dp)

 !      dp(ixO^S) = dabs(w(ixR^S,pp_)-w(ixL^S,pp_))&
 !           /(w(ixO^S,rho_)*dp(ixO^S))
 !      if (useprimitiveRel) then
 !         v(ixI^S)  = w(ixI^S,u0_+idims)/w(ixI^S,lfac_)
 !      else
 !         v(ixI^S)  = w(ixI^S,v0_+idims)
 !      end if
 !      call gradient(v,ixI^L,ixO^L,idims,dv)
 !   end if

 ! end subroutine ppmflatsh
  !=============================================================================
  subroutine srhd_get_v(w,x,ixI^L,ixO^L,idim,v)

    ! Calculate v_idim=m_idim/rho within ixO^L

    use mod_global_parameters

    integer, intent(in)           :: ixI^L, ixO^L, idim
    double precision, intent(in)  :: w(ixI^S,nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out) :: v(ixG^T)
    !-----------------------------------------------------------------------------
    ! assuming getv FOLLOWS a getaux call for updated p_ entry and
    ! well-behaved xi=d+tau+p

    ! first case only happens when d=0 and p/tau are at enforced minimal
    ! values, namely p=smallp and tau=smallp/(gamma-1) (no velocities)
    ! This will typically NEVER be the case, but still...
    where(w(ixO^S,d_)+w(ixO^S,tau_)+w(ixO^S,p_)==smallxi)
       v(ixO^S)=zero
    elsewhere
       v(ixO^S)=w(ixO^S,s0_+idim)/ &
            (w(ixO^S,d_)+w(ixO^S,tau_)+w(ixO^S,p_))
    endwhere

  end subroutine srhd_get_v
  !=============================================================================
  subroutine srhd_get_cmax(new_cmax,w,x,ixI^L,ixO^L,idims,cmax,cmin,needcmin)

    ! Calculate cmax_idim=csound+abs(v_idim) within ixO^L

    use mod_global_parameters

    logical, intent(in)                               :: new_cmax,needcmin
    integer, intent(in)                               :: ixI^L, ixO^L, idims
    double precision, dimension(ixI^S,nw), intent(in) :: w
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, dimension(ixG^T), intent(out)   :: cmax,cmin

    double precision, dimension(ixG^T)                :: csound2,rhoh,vidim2,v2,vidim
    !-----------------------------------------------------------------------------
    !== ZM calculation of sound speed using the EOS ==!
    call getcsound2(w,ixI^L,ixO^L,.true.,rhoh,csound2)

    if(.not.needcmin)then
       v2(ixO^S)=(sum(w(ixO^S,rmom(:))**2.0d0)/ &
            ((rhoh(ixO^S)*w(ixO^S,lfac_)*w(ixO^S,lfac_))**2.0d0)
 !      v2(ixO^S)=({^C&w(ixO^S,s^C_)**2.0d0+})/ &
 !           ((rhoh(ixO^S)*w(ixO^S,lfac_)*w(ixO^S,lfac_))**2.0d0)
    else
       vidim(ixO^S)=(w(ixO^S,s0_+idims)/ &
            (rhoh(ixO^S)*w(ixO^S,lfac_)*w(ixO^S,lfac_)))
    endif

    if(.not.needcmin)then
       if (ndir==1) then
          cmax(ixO^S)= (dsqrt(v2(ixO^S))+dsqrt(csound2(ixO^S)))/ &
               (one+dsqrt(csound2(ixO^S)*v2(ixO^S)))
       else
          vidim2(ixO^S)=(w(ixO^S,s0_+idims)/ &
               (rhoh(ixO^S)*w(ixO^S,lfac_)*w(ixO^S,lfac_)))**2.0d0

          cmax(ixO^S)=( dsqrt(vidim2(ixO^S))*(one-csound2(ixO^S)) + &
               dsqrt(csound2(ixO^S)*(one-v2(ixO^S))*( &
               one-v2(ixO^S)*csound2(ixO^S)-vidim2(ixO^S)*(one-csound2(ixO^S))) &
               ) ) / (one-v2(ixO^S)*csound2(ixO^S))
       endif
    else
       if (ndir==1) then
          cmax(ixO^S)= min(one,max(zero,(vidim(ixO^S)+dsqrt(csound2(ixO^S)))/ &
               (one+dsqrt(csound2(ixO^S))*vidim(ixO^S))))
          cmin(ixO^S)= max(-one,min(zero,(vidim(ixO^S)-dsqrt(csound2(ixO^S)))/ &
               (one-dsqrt(csound2(ixO^S))*vidim(ixO^S))))
       else
          v2(ixO^S)=(sum(w(ixO^S,mom(:))**2.0d0)/ &
               (rhoh(ixO^S)*w(ixO^S,lfac_)*w(ixO^S,lfac_))**2.0d0
 !         v2(ixO^S)=({^C&w(ixO^S,s^C_)**2.0d0+})/ &
 !              (rhoh(ixO^S)*w(ixO^S,lfac_)*w(ixO^S,lfac_))**2.0d0

          cmax(ixO^S)=min(one,max(zero,(vidim(ixO^S)*(one-csound2(ixO^S)) + &
               dsqrt(csound2(ixO^S)*(one-v2(ixO^S))*( &
               one-v2(ixO^S)*csound2(ixO^S)-vidim(ixO^S)**2.0d0*(one-csound2(ixO^S))) &
               ) ) / (one-v2(ixO^S)*csound2(ixO^S))))

          cmin(ixO^S)=max(-one,min(zero,(vidim(ixO^S)*(one-csound2(ixO^S)) - &
               dsqrt(csound2(ixO^S)*(one-v2(ixO^S))*( &
               one-v2(ixO^S)*csound2(ixO^S)-vidim(ixO^S)**2.0d0*(one-csound2(ixO^S))) &
               ) ) / (one-v2(ixO^S)*csound2(ixO^S))))
       endif
    endif

  end subroutine srhd_get_cmax
  !=============================================================================
  subroutine srhd_get_flux(w,x,ixI^L,ixO^L,iw,idim,f,transport)

    ! Calculate non-transport flux f_idim[iw] within ixO^L.

    use mod_global_parameters

    integer, intent(in)           :: ixI^L, ixO^L, iw, idim
    double precision, intent(in)  :: w(ixI^S,nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out) ::  f(ixG^T)
    logical, intent(out)          :: transport
    !-----------------------------------------------------------------------------
    if(iw==s0_+idim)then
       ! f_i[s_i]=v_i*s_i + p
       f(ixO^S)=w(ixO^S,p_) 
!       {#IFDEF TRACER
!       {else if (iw==tr^FL_) then 
!       f(ixO^S)=zero\}
!       }

    else if(iw==tau_)then
       ! f_i[tau]=v_i*tau + v_i*p
       ! first case only happens when d=0 and p/tau are at enforced minimal
       ! values, namely p=smallp and tau=smallp/(gamma-1) (no velocities)
       ! This will typically NEVER be the case, but still...
       where(w(ixO^S,d_)+w(ixO^S,tau_)+w(ixO^S,p_)==smallxi)
          f(ixO^S)=zero
       elsewhere
          f(ixO^S)=w(ixO^S,s0_+idim)*w(ixO^S,p_)/ &
               (w(ixO^S,d_)+w(ixO^S,tau_)+w(ixO^S,p_))
       endwhere
    else
       f(ixO^S)=zero
    endif

    transport=.true.

  end subroutine srhd_get_flux
  !=============================================================================
  subroutine srhd_get_flux_forhllc(w,x,ixI^L,ixO^L,iw,idims,f,transport)

    ! Calculate non-transport flux f_idim[iw] within ixO^L.

    use mod_global_parameters

    integer, intent(in)           :: ixI^L,ixO^L,iw,idims
    double precision, intent(in)  :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out) :: f(ixG^T,1:nwflux)
    logical, intent(out)          :: transport
    !----------------------------------------------

    ! assuming getflux FOLLOWS a getaux call for updated p_ entry and
    ! well-behaved xi=d+tau+p, write:

    if(iw==s0_+idims)then
       ! f_i[s_i]=v_i*s_i + p
       f(ixO^S,iw)=w(ixO^S,p_)
!       {#IFDEF TRACER
!       {else if (iw==tr^FL_) then 
!       f(ixO^S,iw)=zero\}
!       }

    else if(iw==tau_)then
       ! f_i[tau]=v_i*tau + v_i*p
       ! first case only happens when d=0 and p/tau are at enforced minimal
       ! values, namely p=smallp and tau=smallp/(gamma-1) (no velocities)
       ! This will typically NEVER be the case, but still...
       where(w(ixO^S,d_)+w(ixO^S,tau_)+w(ixO^S,p_)==smallxi)
          f(ixO^S,iw)=zero
       elsewhere
          f(ixO^S,iw)=w(ixO^S,s0_+idims)*w(ixO^S,p_)/ &
               (w(ixO^S,d_)+w(ixO^S,tau_)+w(ixO^S,p_))
       endwhere
    else
       f(ixO^S,iw)=zero
    endif

    transport=.true.

  end subroutine srhd_get_flux_forhllc
  !=============================================================================
  subroutine srhd_con2prim(pressure,lfac,d,s^C,tau,ierror)
    !use ieee_arithmetic
    use mod_global_parameters

    double precision:: pressure,lfac
    double precision:: d,s^C,tau
    integer:: ierror

    integer:: ni,niiter
    double precision:: pcurrent,pnew
    double precision:: er,er1,ff,df,dp,v^C
    double precision:: pmin,lfac2inv,pLabs,pRabs,pprev
    double precision:: s2overcubeG2rh,sqrs
    double precision:: xicurrent
    double precision:: oldff1,oldff2
    double precision:: Nff
    double precision:: pleft,pright,pnewi
    integer::nit,n2it,ni2,ni3
    double precision:: h,dhdp
    !-----------------------------------------------------------------------------

    ierror=0
    ! ierror=0 : ok
    !
    ! ierror<>0
    !
    ! ierror=1 : error on entry: must have D>=minrho, tau>=smalltau
    ! ierror=2 : maxitnr reached without convergence
    ! ierror=3 : final pressure value < smallp or xi<smallxi during iteration
    ! ierror=4 : final v^2=1 hence problem as lfac=1/0
    ! ierror=5 : nonmonotonic function f?
    ! ierror=7 : stop due to strictnr violation

    if(d<minrho .or. tau<smalltau) then
       ierror=1
       return
    endif

    sqrs={s^C**2.0d0+}

    if(sqrs==zero)then
       call pressureNoFlow(pressure,tau,d)
       lfac=one
       return
    endif

    ! left and right brackets for p-range
    pmin=dsqrt(sqrs)/(one-dmaxvel)-tau-d
    pLabs=max(minp,pmin)
    pRabs=1.0d99
    ! start value from input
    !opedit: try to use previous value:
    if (pLabs .le. pressure .and. pressure .le. pRabs) then
       pcurrent = pressure
    else
       pcurrent = pLabs
    end if

    er1=one
    pprev=pcurrent

    ! Fudge Parameters
    oldff1=1.0d7  ! High number
    oldff2=1.0d9  ! High number bigger then oldff1
    n2it = 0
    nit  = 0

    LoopNR:  do ni=1,maxitnr
       nit = nit + 1
       !============= Controle ~1~=============!
       if(nit>maxitnr/4)then
          !print *,'ni,er,p',ni,er,pcurrent
          ! mix pressure value for convergence
          pcurrent=half*(pcurrent+pprev)
          ! relax accuracy requirement
          er1=10.*er1
          nit = nit - maxitnr/10
       endif
       !=======================================!

       niiter=ni  
       xicurrent=tau+d+pcurrent

       if(xicurrent<smallxi) then
          if(strictgetaux)then
             print*,'!--- amrvacphys/t.srhd-- con2prim ---!'
             print *,'stop: too small xi iterate:',xicurrent
             print *,'for pressure iterate p',pcurrent
             print *,'pressure bracket pLabs pRabs',pLabs,pRabs
             print *,'iteration number:',ni
             print *,'values for d,s,tau,s^2:',d,s^C,tau,sqrs
          endif
          ierror=3 
          return
       endif

       {v^C=s^C/xicurrent\}
       lfac2inv=one - ({v^C**2.0d0+})
       if(lfac2inv>zero) then
          lfac=one/dsqrt(lfac2inv)
       else
          if(strictgetaux)then
             print*,'!--- amrvacphys/t.srhd-- con2prim ---!'
             print *,'stop: negative or zero factor 1-v^2:',lfac2inv
             print *,'for pressure iterate p',pcurrent
             print *,'absolute pressure bracket pLabs pRabs',pLabs,pRabs
             print *,'iteration number:',ni
             print *,'values for d,s,tau,s^2:',d,s^C,tau,sqrs
             print *,'values for v,xi:',v^C,xicurrent
          endif
          ierror=4
          return
       endif


       s2overcubeG2rh=sqrs/(xicurrent**3.0d0)
       !== ZM calculation done using the EOS ==!
       call FuncEnthalpy(pcurrent,lfac2inv,d,s^C,tau,sqrs,xicurrent,&
            s2overcubeG2rh,h,dhdp,ierror)
       !=======================================!   
       ff=-xicurrent*lfac2inv + h 
       df=- two*sqrs/(xicurrent)**2.0d0  + dhdp - lfac2inv

       if (ff*df==zero) then
          if (ff==zero) then
             exit ! zero found
          else
             if(strictgetaux)print *,'stop: df becomes zero, non-monotonic f(p)!!'
             ierror=5
             return
          endif
       else 
          pnew=pcurrent-ff/df
          if (ff*df>zero) then
             ! pressure iterate has decreased
             ! restrict to left 
             pnew=max(pnew,pLabs)
          else  ! ff*df<0
             ! pressure iterate has increased
             ! restrict to right 
             pnew=min(pnew,pRabs)
          endif
       endif


       ! handle special case where NR incorrectly believes in convergence
       if(pnew == pLabs .and. pcurrent==pnew .and. &
            abs(ff)> absaccnr .and. sqrs > zero)then
          pnewi=pnew
          ! try 2 higher pressure values to locate a sign change for f(p)
          LoopCor:  do ni2=1,2
             !=====================!
             pcurrent=pnewi*500.0d0
             xicurrent=tau+d+pcurrent
             {v^C=s^C/xicurrent\}
             lfac2inv=one - ({v^C**2.0d0+})
             !=====================!

             !=====================!
             if(lfac2inv>zero)then
                lfac=one/dsqrt(lfac2inv)
             else
                ierror=4
                return
             endif
             !=====================!

             s2overcubeG2rh=-sqrs/(xicurrent**3.0d0)
             !==== Calculate enthalpy and derivative ====!
             call Bisection_Enthalpy(pcurrent,lfac2inv,d,s^C,&
                  tau,sqrs,xicurrent,h,ierror)
             Nff=-xicurrent*lfac2inv + h

             !== Save old value of pressure ==!
             pnewi=pcurrent
             !================================!

             !== find the interval where is the root ==!
             if(Nff * ff <=zero)then
                pnew=pcurrent
                exit LoopCor
             endif
             !=========================================!
          enddo LoopCor

          !== No possible solution, correct all including the conservatives ==!
          if( Nff*ff>zero)then

             ! following is in accord with trick done in smallvalues
             d   = 2.0d0*(one + 10.0d0 * minrho) * minrho
             tau = 2.0d0*(one + 10.0d0 * smalltau) * smalltau
             {^C&s^C =zero;}
             pressure     = (eqpar(gamma_)-one)*tau
             lfac = one

             if(strictnr)ierror=7
             ! leave the do loop here
             return
          endif
       endif
       !===============================================!
       dp=pcurrent-pnew
       er=two*dabs(dp)/(pnew+pcurrent)
       if(((er<tolernr*er1).or.(dabs(dp)<absaccnr))) exit LoopNR
       !===============================================!

       ! For very small values of pressure, NR algorithm is not efficient to
       ! find root, use Euler algorithm to find precise value of pressure
       if((dabs(oldff2-ff) < 1.0d-8 .or. niiter >= maxitnr-maxitnr/20).and.&
            ff * oldff1 < zero    .and.  dabs(ff)>absaccnr)then

          n2it=n2it+1
          if(n2it <= 3) pcurrent=half*(pnew+pcurrent)
          if(n2it >3)then
             pright =pcurrent
             pleft=pprev
             pcurrent=half*(pleft+pright)
             Dicho:  do ni3=1,maxitnr
                !===================!
                xicurrent=tau+d+pcurrent
                {v^C=s^C/xicurrent\}
                lfac2inv=one - ({v^C**2.0d0+})
                if(lfac2inv>zero)then
                   lfac=one/dsqrt(lfac2inv)
                else
                   ierror=4
                   return
                endif
                !===================!


                !== ZM calculation done using the EOS ==!
                call Bisection_Enthalpy(pnew,lfac2inv,d,s^C,&
                     tau,sqrs,xicurrent,h,ierror)
                Nff=-xicurrent*lfac2inv + h 
                !=======================================!
                !==== Iterate ====!
                if(ff * Nff < zero)then
                   pleft=pcurrent
                else
                   pright=pcurrent
                endif

                pcurrent=half*(pleft+pright)
                !==================!

                !=== The iteration converge ===!
                if(2.0d0*dabs(pleft-pright)/(pleft+pright) < absaccnr &
                     .or. dabs(ff)<absaccnr)then
                   pnew=pcurrent
                   exit LoopNR
                endif
                !==============================!

                !=== conserve the last value of Nff ===!
                ff=Nff
                !======================================!
             enddo    Dicho
          endif

       else
          !====== There is no problems, continue the NR iteration ======!
          pprev=pcurrent
          pcurrent=pnew
          !=============================================================!
       endif


       !=== keep the values of the 2 last ff ===!
       oldff2=oldff1
       oldff1=ff
       !========================================!
    enddo LoopNR

    if(niiter==maxitnr)then
       !print*,' ff = ',ff,' df = ',df
       !print*,'reachs maxitnr = ', niiter
       ierror=2
       return
    endif

    if(pcurrent<minp) then
       ierror=3
       return
    endif

    !--end result for pressure and lorentz factor------!
    pressure=pcurrent
    xicurrent=tau+d+pressure
    {v^C = s^C/xicurrent\}
    lfac2inv=one - ({v^C**2.0d0+})
    if(lfac2inv>zero) then
       lfac=one/dsqrt(lfac2inv)
    else
       ierror=4
       return
    endif
    !------------------------------!

  end subroutine srhd_con2prim
  !=============================================================================
  subroutine srhd_add_geometry(qdt,ixI^L,ixO^L,wCT,w,x)

    ! Add geometrical source terms to w

    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(in)    :: wCT(ixI^S,1:nw)
    double precision, intent(inout) ::  w(ixI^S,1:nw)

    integer          :: iw,idir, h1x^L{^NOONED, h2x^L}
    double precision :: tmp(ixG^T)
    logical          :: angmomfix=.false.
    !-----------------------------------------------------------------------------
!!! now handled in tvdmusclf/getaux/hll variants
!!!if(typeaxial /= 'slab') call getaux(.true.,wCT,ixI^L,ixO^L,'addgeometry')

    select case (typeaxial)
    case ('slab')
       ! No source terms in slab symmetry
    case ('spherical')
       h1x^L=ixO^L-kr(1,^D); {^NOONED h2x^L=ixO^L-kr(2,^D);}
       do iw=1,nwflux
          select case(iw)
             ! s[s1]=((Stheta**2+Sphi**2)/xi+2*p)/r
          case(s1_)
             tmp(ixO^S)=wCT(ixO^S,p_)*x(ixO^S,1) &
                  *(mygeo%surfaceC1(ixO^S)-mygeo%surfaceC1(h1x^S)) &
                  /mygeo%dvolume(ixO^S){&^CE&
                  +wCT(ixO^S,s^CE_)**2.0d0&
                  /(wCT(ixO^S,tau_)+wCT(ixO^S,d_)+wCT(ixO^S,p_)) }
             {^NOONEC
             ! s[s2]=-(Sr*Stheta/xi)/r
             !       + cot(theta)*(Sphi**2/xi+p)/r
          case(s2_)
             }
             {^NOONED
             tmp(ixO^S) = +wCT(ixO^S,p_)
             w(ixO^S,iw)=w(ixO^S,iw) &
                  +qdt*tmp(ixO^S)*(mygeo%surfaceC2(ixO^S)-mygeo%surfaceC2(h2x^S)) &
                  /mygeo%dvolume(ixO^S)
             }
             {^NOONEC
             tmp(ixO^S)=-(wCT(ixO^S,s1_)*wCT(ixO^S,s2_)&
                  /(wCT(ixO^S,tau_)+wCT(ixO^S,d_)+wCT(ixO^S,p_)))
             }
             {^IFTHREEC
             {^NOONED
             tmp(ixO^S)=tmp(ixO^S)+(wCT(ixO^S,s3_)**2.0&
                  /(wCT(ixO^S,tau_)+wCT(ixO^S,d_)+wCT(ixO^S,p_))) &
                  *dcos(x(ixO^S,2))/dsin(x(ixO^S,2))
             }
             ! s[s3]=-(sphi*sr/xi)/r
             !       -cot(theta)*(stheta*sphi/xi)/r
          case(s3_)
             if(.not.angmomfix) &
                  tmp(ixO^S)=-(wCT(ixO^S,s3_)*wCT(ixO^S,s1_)&
                  /(wCT(ixO^S,tau_)+wCT(ixO^S,d_)+wCT(ixO^S,p_))){^NOONED &
                  -(wCT(ixO^S,s2_)*wCT(ixO^S,s3_)&
                  /(wCT(ixO^S,tau_)+wCT(ixO^S,d_)+wCT(ixO^S,p_))) &
                  *dcos(x(ixO^S,2))/dsin(x(ixO^S,2)) }
             }
          end select
          ! Divide by radius and add to w
          if(iw==s1_{^NOONEC.or.iw==s2_}{^IFTHREEC  &
               .or.(iw==s3_.and..not.angmomfix)}) &
               w(ixO^S,iw)=w(ixO^S,iw)+qdt*tmp(ixO^S)/x(ixO^S,1)
       end do
    case ('cylindrical')
       do iw=1,nwflux
          select case (iw)
             ! source[sr]=sphi*vphi/radius + p/radius
          case (sr_)
             w(ixO^S,sr_)=w(ixO^S,sr_)+qdt*wCT(ixO^S,p_)/x(ixO^S,1)
             tmp(ixO^S)=zero
             {^IFPHI
             tmp(ixO^S)=tmp(ixO^S)+wCT(ixO^S,sphi_)**2/ &
                  (wCT(ixO^S,tau_)+wCT(ixO^S,d_)+wCT(ixO^S,p_))
             ! source[sphi]=(-sphi*vr)/radius
          case (sphi_)
             tmp(ixO^S)=-wCT(ixO^S,sphi_)*wCT(ixO^S,sr_)/ &
                  (wCT(ixO^S,tau_)+wCT(ixO^S,d_)+wCT(ixO^S,p_))
             }
          end select
          ! Divide by radius and add to w
          if (iw==sr_.or.iw==sphi_) then
             w(ixO^S,iw)=w(ixO^S,iw)+qdt*tmp(ixO^S)/x(ixO^S,1)
          end if
       end do
    end select

  end subroutine srhd_add_geometry
!=============================================================================
  subroutine srhd_correctaux(ixI^L,ixO^L,w,x,patchierror,subname)

    use mod_global_parameters

    integer, intent(in)         :: ixI^L, ixO^L
    integer, intent(in)         :: patchierror(ixG^T)
    character(len=*), intent(in)   :: subname
    double precision, intent(inout):: w(ixI^S,1:nw)
    double precision, intent(in)      :: x(ixI^S,1:ndim)

    integer        :: iw, kxO^L, ix^D, i
    logical        :: patchw(ixG^T)
    !-----------------------------------------------------------------------------

    {do ix^D= ixO^LIM^D\}
    ! point with local failure identified by patchierror
    ! patchierror=-1 then from smallvalues call
    ! patchierror=1  then from getaux      call
    if (patchierror(ix^D)/=0) then
       ! verify in cube with border width nflatgetaux the presence
       ! of cells where all went ok
       do i=1,nflatgetaux
          {kxOmin^D= max(ix^D-i,ixOmin^D);
          kxOmax^D= min(ix^D+i,ixOmax^D);\}
          ! in case cells are fine within smaller cube than 
          ! the userset nflatgetaux: use that smaller cube
          if (any(patchierror(kxO^S)==0)) exit
       end do
       if (any(patchierror(kxO^S)==0))then
          ! within surrounding cube, cells without problem were found
          if (useprimitive) then
             patchw(kxO^S)=(patchierror(kxO^S)/=0)
             call primitiven(ixI^L,kxO^L,w,patchw)
          end if
          ! faulty cells are corrected by averaging here
          ! only average those which were ok and replace faulty cells
          do iw = 1,nw
             w(ix^D,iw)=sum(w(kxO^S,iw),patchierror(kxO^S)==0)&
                  /count(patchierror(kxO^S)==0)
          end do
          if (useprimitive) then
             ! in addition to those switched to primitive variables
             ! above, also switch the corrected variables
             patchw(ix^D)=.false.
             call conserven(ixI^L,kxO^L,w,patchw)
          end if
       else
          ! no cells without error were found in cube of size nflatgetaux
          ! --> point of no recovery
          write(*,*)'Getaux error:',patchierror(ix^D),'ix^D=',ix^D
          !write(*,*)'New ','d=',w(ix^D,d_),'s=', &
          !        {^C&w(ix^D,s^C_)},'tau=',w(ix^D,tau_)
          !write(*,*)'position  ', px(saveigrid)%x(ix^D, 1:ndim)
          write(*,*)'Called from: ',subname
          if (patchierror(ix^D)<0) then
             call mpistop("---------correctaux from smallvalues-----")
          else
             call mpistop("---------correctaux from getaux----------")
          end if
       end if
    end if
    {enddo^D&\}

  end subroutine srhd_correctaux
  !=============================================================================
  subroutine srhd_handle_small_values(w,x,ixI^L,ixO^L,subname)

    use mod_global_parameters

    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    character(len=*), intent(in)    ::subname

    !!integer                         :: posvec(ndim)
    integer, dimension(ixG^T)       :: patchierror
    !-----------------------------------------------------------------------------

    if(any(w(ixO^S,d_) < minrho) .or. any(w(ixO^S,tau_) < smalltau))then
       if(strictsmall)then
          write(*,*)'smallvalues :: for tau = ',minval(w(ixO^S,tau_)), &
               'With limit=',smalltau,' For D =', minval(w(ixO^S,d_)),&
               ' With limit=',smallrho
          write(*,*)'From::  ', subname,' iteration ', it
          call mpistop("Smallvalues with strictsmall=T failed")
       else
          if(strictgetaux)then
             where(w(ixO^S,d_) < minrho .or. w(ixO^S,tau_) < smalltau)
                w(ixO^S,d_)  = 2.0*(1.0d0 + 10.0d0 * minrho)*minrho
                w(ixO^S,tau_)= 2.0*(1.0d0 + 10.0d0 * minp)*smalltau
                {^C&w(ixO^S,s^C_) =zero;}
                w(ixO^S,lfac_)=one
             end where
          else
             where(w(ixO^S,d_) < minrho .or. w(ixO^S,tau_) < smalltau)
                patchierror(ixO^S)=-1
             else where
                patchierror(ixO^S)=0
             end where
             call correctaux(ixI^L,ixO^L,w,x,patchierror,subname)
          end if
       end if ! strict
    end if

    {#IFDEF TRACER
    if(any(w(ixO^S,Dtr1_) < minrho) .or. &
         any(w(ixO^S,Dtr1_) > 10.d10 ))then
       where(w(ixO^S,Dtr1_) < minrho .or. &
            w(ixO^S,Dtr1_) > 10.d10) 
          w(ixO^S,Dtr1_) = 0.d0
       end where
    end if
    if(any(w(ixO^S,Dtr1_) .NE. w(ixO^S,Dtr1_)))then
       where(w(ixO^S,Dtr1_) .NE. w(ixO^S,Dtr1_)) 
          w(ixO^S,Dtr1_) = 0.d0
       end where
    end if\}

    end subroutine srhd_handle_small_values
  !=============================================================================
  subroutine srhd_getaux(clipping,w,x,ixI^L,ixO^L,subname)

    ! Calculate auxilary variables ixO^L from non-auxiliary entries in w
    ! clipping can be set to .true. to e.g. correct unphysical pressures,
    ! densities, v>c,  etc.

    use mod_global_parameters
    use mod_small_values
    logical, intent(in)             :: clipping
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S,nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    character(len=*), intent(in)    :: subname

    integer          :: err,ix^D
    double precision :: dold,tauold,{^C&sold^C_},pold,lfacold
    integer          :: patchierror(ixG^T)
    !-----------------------------------------------------------------------------

    ! artificial replacement of small D and tau values
    ! with corresponding smallrho/smallp settings,
    ! together with nullifying momenta
    call srhd_handle_small_values(w,x,ixI^L,ixO^L,subname)

    ! we compute auxiliaries p, lfac from D,S,tau
    ! put the p and lfac in the auxiliary fields lfac_ and p_
    ! on entry: p_ field may contain previous value for iteration
    ! however, for filling ghost cells, this does not exist, so we no longer
    ! use it

    {do ix^D= ixO^LIM^D\}
    dold=w(ix^D,d_)
    tauold=w(ix^D,tau_)
    pold=w(ix^D,p_)
    lfacold=w(ix^D,lfac_)
    { ^C&sold^C_=w(ix^D,s^C_);}

    call srhd_con2prim(w(ix^D,p_),w(ix^D,lfac_), &
         w(ix^D,d_),{^C&w(ix^D,s^C_)},w(ix^D,tau_),err)
    patchierror(ix^D)=err
    if (err/=0.and.strictgetaux) then
       write(*,*)'Getaux error:',err,'ix^D=',ix^D,' global timestep it=',it
       write(*,*)'Postion: ',x(ix^D,1:ndim)
       write(*,*)'start value for p=',pold
       write(*,*)'start value for lfac=',lfacold
       write(*,*)'end value for p=',w(ix^D,p_)
       write(*,*)'end value for lfac=',w(ix^D,lfac_)
       write(*,*)'exit  d=',w(ix^D,d_),'s=',{^C&w(ix^D,s^C_)},'tau=',w(ix^D,tau_)
       write(*,*)'input d=',dold,'s=',{^C&sold^C_},'tau=',tauold
       write(*,*)'Called from: ',subname
       call mpistop("problem in getaux")
    endif
    {enddo^D&\}

    if(.not.strictgetaux.and.any(patchierror(ixO^S)/=0)) &
         call srhd_correctaux(ixI^L,ixO^L,w,x,patchierror,subname)

  end subroutine srhd_getaux
!=============================================================================
end module mod_srhd_phys
