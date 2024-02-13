!> Jet setup, SRHD 2D cylindrical version
!> Using ApJ 920, 144 (2021) as guideline
module mod_usr
  use mod_srhd

  implicit none

  double precision :: rjet,rhob,etarho,rhojet,pb,zetap,pjet,lfacjet,vjet
  double precision :: Qjet, Mdotjet, p0val, t0val, tcross, vhead

contains

  !> Read this module's parameters from a file
  subroutine usr_params_read(files)
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /usr_list/ rjet, rhob, etarho, pb, zetap, lfacjet

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, usr_list, end=111)
111    close(unitpar)
    end do

  end subroutine usr_params_read

  subroutine usr_init()
    use mod_variables

    call usr_params_read(par_files)

    ! units IN CGS
    unit_numberdensity=0.001d0  ! 10^-3 per cubic cm
    unit_length=3.086d+21       ! 1 kpc in cm

    usr_set_parameters  => initglobaldata_usr
    usr_init_one_grid   => Jet_init_one_grid
    usr_special_bc      => specialbound_usr
    usr_refine_grid     => specialrefine_grid
    usr_aux_output      => specialvar_output
    usr_add_aux_names   => specialvarnames_output

    call set_coordinate_system("cylindrical_2D")

    call srhd_activate()

  end subroutine usr_init

  subroutine initglobaldata_usr
    character(len=20) :: printsettingformat

    printsettingformat='(1x,A50,ES15.7,A7)'

    if(mype==0) then
      write(*,*) "Jet (Seo et al ApJ 2021, 920, 144) setup:"
      write(*,printsettingformat) "dimensionless jet radius ",rjet," input"
      write(*,printsettingformat) "dimensionless density in medium ",rhob," input"
      write(*,printsettingformat) "density in jet: via  ",etarho," input"
      write(*,printsettingformat) "dimensionless pressure background ",pb," input"
      write(*,printsettingformat) "pressure in jet: via ",zetap," input"
      write(*,printsettingformat) "jet lorentz factor ",lfacjet," input"
    end if

    rhojet=etarho*rhob
    pjet=zetap*pb
    vjet=dsqrt(1.0d0-1.0d0/lfacjet**2)

    if(mype==0) then
      write(*,*) "Deduced dimensionless values:"
      write(*,printsettingformat) "jet density  ",rhojet," output"
      write(*,printsettingformat) "jet pressure ",pjet," output"
      write(*,printsettingformat) "jet velocity ",vjet," output"
    end if

  end subroutine initglobaldata_usr

  !> Initialize one grid
  subroutine Jet_init_one_grid(ixG^L,ix^L,w,x)
    integer, intent(in) :: ixG^L, ix^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)
    
    double precision :: R(ixG^S),Z(ixG^S)
    double precision,dimension(^D&1) :: rhoL,pL,rhohL
    double precision :: hjet,bracket,hb,etar,vhead,tcross
    integer :: ixL^L,itr
    logical :: first
    data first/.true./
    character(len=20) :: printsettingformat

    printsettingformat='(1x,A50,ES15.7,A7)'

    if(mype==0.and.first)then
      ixL^L=1;
      rhoL(^D&1)=rhojet
      pL(^D&1)=pjet
      call srhd_get_enthalpy_eos(ixL^L,rhoL,pL,rhohL)
      hjet=rhohL(^D&1)/rhojet
      bracket=(lfacjet**2*rhojet*hjet-lfacjet*rhojet)*unit_density*unit_velocity**2
      Qjet=dpi*(rjet*unit_length)**2*vjet*unit_velocity*bracket 
      bracket=(lfacjet**2*(rhojet*unit_density)*hjet*(vjet*unit_velocity)**2+pjet*unit_pressure)
      Mdotjet=dpi*(rjet*unit_length)**2*bracket
      p0val=rhob*unit_density*unit_velocity**2
      t0val=rjet*unit_length/unit_velocity
      ! switch to years
      t0val=t0val/(365.0d0*24.0d0*60.0d0*60.0d0)
      
      rhoL(^D&1)=rhob
      pL(^D&1)=pb
      call srhd_get_enthalpy_eos(ixL^L,rhoL,pL,rhohL)
      hb=rhohL(^D&1)/rhob
      etar=rhojet*hjet*lfacjet**2/(rhob*hb)
      vhead=(vjet*unit_velocity)*dsqrt(etar)/(dsqrt(etar)+one)
      tcross=rjet*unit_length/vhead
      ! switch to years
      tcross=tcross/(365.0d0*24.0d0*60.0d0*60.0d0)

      write(*,*) "Deduced values (with dimensions):"
      write(*,printsettingformat) "p0 value  ",p0val," output"
      write(*,printsettingformat) "t0 value in years ",t0val," output"
      write(*,printsettingformat) "tcross in years ",tcross," output"
      write(*,printsettingformat) "one crossing time is in code units",rjet*unit_length/vhead/unit_time," output"
      write(*,printsettingformat) "jet power  ",Qjet," output"
      write(*,printsettingformat) "jet thrust ",Mdotjet," output"
      write(*,printsettingformat) "pb in cgs",pb*unit_pressure," output"
      write(*,printsettingformat) "Tb in K",pb*unit_pressure &
                 /((2.0d0+3.0d0*He_abundance)*unit_numberdensity*kB_cgs)," output"
      write(*,*) "Deduced values (dimensionless):"
      write(*,printsettingformat) "jet enthalpy ",hjet," output"
      write(*,printsettingformat) "background enthalpy ",hb," output"
      write(*,printsettingformat) "vhead/c ",vhead/unit_velocity," output"
      write(*,printsettingformat) "pb/p0",pb*unit_pressure/p0val," in-output"

      write(*,*) "units from"
      write(*,*) "SI_unit ",SI_unit," should be false!!"
      write(*,printsettingformat) "unit_velocity ",unit_velocity," fixed"
      write(*,printsettingformat) "unit n ",unit_numberdensity," input"
      write(*,printsettingformat) "unit length ",unit_length," input"
      write(*,printsettingformat) "He-abundance ",He_abundance," input"
      write(*,printsettingformat) "unit rho ",unit_density," output"
      write(*,printsettingformat) "unit p ",unit_pressure," output"
      write(*,printsettingformat) "unit time ",unit_time," output"
      first=.false.
    endif

    ! radial direction R
    R(ix^S)=x(ix^S,1)
    ! axial direction Z
    Z(ix^S)=x(ix^S,2)

    ! radial velocity vR
    w(ix^S,mom(1))  = 0.0d0
    where(R(ix^S)<rjet .and. Z(ix^S)<0.0d0)
       w(ix^S,rho_) = rhojet
       w(ix^S,p_)   = pjet
       ! beware to fill with four-velocity
       w(ix^S,mom(2))= vjet*lfacjet
    elsewhere
       w(ix^S,rho_) = rhob
       w(ix^S,p_)   = pb
       w(ix^S,mom(2))= 0.0d0
    end where
    if(srhd_n_tracer>0)then
      do itr=1,srhd_n_tracer 
          where(R(ix^S)<rjet .and. Z(ix^S)<0.0d0)
             w(ix^S,tracer(itr))=1.0d0
          elsewhere
             w(ix^S,tracer(itr))=0.0d0
          end where
      enddo
    endif

    ! computes auxiliary lfac and xi from primitives
    ! and replaces rho-p-v*lfac with conservatives
    call srhd_to_conserved(ixG^L,ix^L,w,x)

  end subroutine Jet_init_one_grid

  subroutine specialbound_usr(qt,ixI^L,ixO^L,iB,w,x)
    ! special boundary types, user defined
    integer, intent(in) :: ixO^L, iB, ixI^L
    double precision, intent(in) :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: R(ixI^S)
    integer :: ix2, ixOInt^L, itr

    R(ixO^S)=x(ixO^S,1)

    select case(iB)
    ! special bottom (Z=0) boundary
    ! fixing all profiles within Rjet
    ! (a)symmetry beyond Rjet
    case(3)
      ! switch internal zone above boundary zone to primitive variables
      ixOInt^L=ixO^L;
      ixOIntmin2=ixOmax2+1
      ixOIntmax2=ixOmax2+nghostcells
      call srhd_to_primitive(ixI^L,ixOInt^L,w,x)
      ! prescribe solution at jet inlet
      where(dabs(R(ixO^S)) <= rjet)
         w(ixO^S,rho_)   = rhojet
         w(ixO^S,mom(1)) = zero
         w(ixO^S,mom(2)) = vjet*lfacjet
         w(ixO^S,p_)     = pjet
      endwhere

      ! extrapolate beyond jet radius
      do ix2 = ixOmin2,ixOmax2
         where(dabs(R(ixOmin1:ixOmax1,ix2)) > rjet)
           w(ixOmin1:ixOmax1,ix2,rho_)  = w(ixOmin1:ixOmax1,2*ixOmax2-ix2+1,rho_)
           w(ixOmin1:ixOmax1,ix2,mom(1))= w(ixOmin1:ixOmax1,2*ixOmax2-ix2+1,mom(1))
           w(ixOmin1:ixOmax1,ix2,mom(2))= w(ixOmin1:ixOmax1,2*ixOmax2-ix2+1,mom(2))
           w(ixOmin1:ixOmax1,ix2,p_)    = w(ixOmin1:ixOmax1,2*ixOmax2-ix2+1,p_)
         endwhere
      enddo

      if(srhd_n_tracer>0)then
        do itr=1,srhd_n_tracer 
          where(dabs(R(ixO^S))<=rjet)
             w(ixO^S,tracer(itr))=1.0d0
          elsewhere
             w(ixO^S,tracer(itr))=0.0d0
          end where
        enddo
      endif

      ! switch to conservative variables in internal zone
      call srhd_to_conserved(ixI^L,ixOInt^L,w,x)
      ! switch to conservative variables in ghost cells
      call srhd_to_conserved(ixI^L,ixO^L,w,x)

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
    double precision:: R(ixG^S), Z(ixG^S)

    R(ix^S)=x(ix^S,1)
    Z(ix^S)=x(ix^S,2)

    if (any((R(ix^S) <= 1.1d0*rjet) .and. (Z(ix^S)<= 0.5d0))) refine=1

  end subroutine specialrefine_grid

  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
  ! this subroutine can be used in convert, to add auxiliary variables to the
  ! converted output file, for further analysis using tecplot, paraview, ....
  ! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
  ! the array normconv can be filled in the (nw+1:nw+nwauxio) range with
  ! corresponding normalization values (default value 1)
    use mod_radiative_cooling
    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)

    double precision :: csound2(ixO^S)
    double precision :: pth(ixI^S),geff(ixI^S)
    double precision :: wlocal(ixI^S,1:nw)

    wlocal(ixI^S,1:nw)=w(ixI^S,1:nw)
    ! output pressure
    call srhd_get_pthermal(wlocal,x,ixI^L,ixO^L,pth)
    w(ixO^S,nw+1)=pth(ixO^S)

    ! output mach number v_z/c_sound
    call srhd_get_csound2(wlocal,x,ixI^L,ixO^L,csound2)
    w(ixO^S,nw+2)=wlocal(ixO^S,mom(2))/(wlocal(ixO^S,xi_)*dsqrt(csound2(ixO^S)))

    ! output effective gamma
    call srhd_get_Geff_eos(wlocal,ixI^L,ixO^L,.true.,geff)
    w(ixO^S,nw+3)=geff(ixO^S)

    ! output log10(rho)
    w(ixO^S,nw+4)=dlog10(w(ixO^S,d_)/w(ixO^S,lfac_))

  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
    ! newly added variables need to be concatenated with the w_names/primnames string
    character(len=*) :: varnames

    varnames='pthermal vzmach geff log10rho'

  end subroutine specialvarnames_output

end module mod_usr
