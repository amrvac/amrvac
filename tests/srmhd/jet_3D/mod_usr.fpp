#:mute
#:include "mod_gpu_directives.fpp"
#:endmute

module mod_usr

  use mod_amrvac
  use mod_physics

  implicit none

! User parameters.
  double precision :: rjet,zjet,rhob,etarho,rhojet,pb,zetap,pjet,lfacjet,vjet
  double precision :: Qjet, Mdotjet, p0val, t0val, tcross, vhead

  ${GPU_DECLARE_CREATE('rjet,zjet,rhob,etarho,rhojet,pb,zetap,pjet,lfacjet,vjet')}$
  ${GPU_DECLARE_CREATE('Qjet,Mdotjet,p0val,t0val,tcross,vhead')}$

contains

  subroutine usr_init()
    use mod_global_parameters
    implicit none

    call usr_params_read(par_files)

    ! units IN CGS
    unit_numberdensity=0.001d0  ! 10^-3 per cubic cm
    unit_length=3.086d+21       ! 1 kpc in cm

    call set_coordinate_system("Cartesian_3D")
    usr_set_parameters => set_parameters_usr
    usr_init_one_grid  => initonegrid_usr
    usr_aux_output     => specialvar_output
    usr_add_aux_names  => specialvarnames_output
    usr_process_adv_global => apply_floor

    call phys_activate()

  end subroutine usr_init

  subroutine usr_params_read(files)
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /usr_list/  rjet, zjet, rhob, etarho, pb, zetap, lfacjet

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, usr_list, end=111)
111    close(unitpar)
    end do

    ${GPU_UPDATE_DEVICE('rjet,zjet,rhob,etarho,pb,zetap,lfacjet')}$

    rhojet=etarho*rhob
    pjet=zetap*pb
    vjet=dsqrt(1.0d0-1.0d0/lfacjet**2)

    ${GPU_UPDATE_DEVICE('rhojet,pjet,vjet')}$

  end subroutine usr_params_read


  subroutine initonegrid_usr(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
     ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w,x)
    integer, intent(in)             :: ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,&
       ixGmax3, ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
    double precision, intent(in)    :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:nw)
    integer :: itr

    associate(&
      w_ => w(ixGmin1:ixGmax1, ixGmin2:ixGmax2, ixGmin3:ixGmax3, :),&
      x_ => x(ixGmin1:ixGmax1, ixGmin2:ixGmax2, ixGmin3:ixGmax3, :))

     where(x_(:,:,:,3)<zjet.and.(x_(:,:,:,1)**2+x_(:,:,:,2)**2)<rjet**2)
       w_(:,:,:,iw_rho)    = rhojet
       w_(:,:,:,iw_e)      = pjet
       w_(:,:,:,iw_mom(1)) = 0.0d0
       w_(:,:,:,iw_mom(2)) = 0.0d0
       w_(:,:,:,iw_mom(3)) = vjet*lfacjet
     elsewhere
       w_(:,:,:,iw_rho)    = rhob
       w_(:,:,:,iw_e)      = pb
       w_(:,:,:,iw_mom(1)) = 0.0d0
       w_(:,:,:,iw_mom(2)) = 0.0d0
       w_(:,:,:,iw_mom(3)) = 0.0d0
     endwhere

     w_(:,:,:,iw_mag(1)) = 0.0d0
     w_(:,:,:,iw_mag(2)) = 0.0d0
     w_(:,:,:,iw_mag(3)) = 0.0d0

     #:if defined(N_TRACER)
     if(srmhd_n_tracer>0)then
       do itr=1,srmhd_n_tracer
          where(x_(:,:,:,3)<zjet.and.(x_(:,:,:,1)**2+x_(:,:,:,2)**2)<rjet**2)
            w_(:,:,:,tracer(itr)) = w_(:,:,:,iw_rho)*lfacjet
          elsewhere
            w_(:,:,:,tracer(itr)) = 0.0d0
          endwhere
       enddo
     endif
     #:endif
    end associate

    call phys_to_conserved(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
       ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w,x)

  end subroutine initonegrid_usr

  subroutine print_params()
    implicit none
    character(len=*), parameter :: fmt = '(A20,1X,SP,ES12.2)'
    character(len=20) :: printsettingformat
    double precision :: rhoL,pL,rhohL,E_th,E
    double precision :: hjet,bracket,hb,etar,vhead,tcross


   printsettingformat='(1x,A50,ES15.7,A7)'

    print *, '------------------ PARAMETERS ------------------'
      write(*,*) "Jet (Seo et al ApJ 2021, 920, 144) setup:"
      write(*,printsettingformat) "dimensionless jet radius ",rjet," input"
      write(*,printsettingformat) "dimensionless jet length ",zjet," input"
      write(*,printsettingformat) "dimensionless density in medium ",rhob," input"
      write(*,printsettingformat) "density in jet: via  ",etarho," input"
      write(*,printsettingformat) "dimensionless pressure background ",pb," input"
      write(*,printsettingformat) "pressure in jet: via ",zetap," input"
      write(*,printsettingformat) "jet lorentz factor ",lfacjet," input"
      print fmt,  'srmhd_gamma',       srmhd_gamma
      write(*,*) "Deduced dimensionless values:"
      write(*,printsettingformat) "jet density  ",rhojet," output"
      write(*,printsettingformat) "jet pressure ",pjet," output"
      write(*,printsettingformat) "jet velocity ",vjet," output"
    print *, '------------------------------------------------'
    write(*,*) "units from"
      write(*,*) "SI_unit ",SI_unit," should be false!!"
      write(*,printsettingformat) "unit_velocity ",unit_velocity," fixed"
      write(*,printsettingformat) "unit n ",unit_numberdensity," input"
      write(*,printsettingformat) "unit length ",unit_length," input"
      write(*,printsettingformat) "He-abundance ",He_abundance," input"
      write(*,printsettingformat) "unit rho ",unit_density," output"
      write(*,printsettingformat) "unit p ",unit_pressure," output"
      write(*,printsettingformat) "unit time ",unit_time," output"
    print *, '------------------------------------------------'

      rhoL=rhojet
      pL=pjet
      !!call srmhd_get_enthalpy_eos(rhoL,pL,rhohL)
      if(srmhd_eos) then
        E_th = pL*1/(srmhd_gamma-1)
        E    = E_th+dsqrt(E_th**2+rhoL**2)
        rhohL = 0.5_dp*((srmhd_gamma+1.0_dp)*E &
                   - (srmhd_gamma-1)*rhoL*(rhoL/E))

      else
        rhohL = rhoL+srmhd_gamma/(srmhd_gamma-1)*pL
      end if
      hjet=rhohL/rhojet
      bracket=(lfacjet**2*rhojet*hjet-lfacjet*rhojet)*unit_density*unit_velocity**2
      Qjet=dpi*(rjet*unit_length)**2*vjet*unit_velocity*bracket
      bracket=(lfacjet**2*(rhojet*unit_density)*hjet*(vjet*unit_velocity)**2+pjet*unit_pressure)
      Mdotjet=dpi*(rjet*unit_length)**2*bracket
      p0val=rhob*unit_density*unit_velocity**2
      t0val=rjet*unit_length/unit_velocity
      ! switch to years
      t0val=t0val/(365.0d0*24.0d0*60.0d0*60.0d0)
      
      rhoL=rhob
      pL=pb
      !!call srmhd_get_enthalpy_eos(rhoL,pL,rhohL)
      if(srmhd_eos) then
        E_th = pL*1/(srmhd_gamma-1)
        E    = E_th+dsqrt(E_th**2+rhoL**2)
        rhohL = 0.5_dp*((srmhd_gamma+1.0_dp)*E &
                   - (srmhd_gamma-1)*rhoL*(rhoL/E))
      else
        rhohL = rhoL+srmhd_gamma/(srmhd_gamma-1)*pL
      end if
      hb=rhohL/rhob
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
    print *, '------------------------------------------------'


  end subroutine

  subroutine set_parameters_usr()
    use mod_global_parameters
    use mod_physics
    implicit none

    if (mype == 0) call print_params()
  end subroutine set_parameters_usr


! NOTE: This must be named exactly `specialbound_usr` in AGILE since during
!   compilation it deals with this function.
  ! NOTE: This is ran on the GPU.
  subroutine specialbound_usr(&
       qt,&
       ixImin1, ixImin2, ixImin3, ixImax1, ixImax2, ixImax3,&
       ixOmin1, ixOmin2, ixOmin3, ixOmax1, ixOmax2, ixOmax3,&
       iB, w, x)
    use mod_global_parameters
    use mod_physics_vars
    !use srmhd_con2prim, only: ideal_eos, synge_eos
    ${GPU_ROUTINE_VECTOR()}$
    implicit none
    integer, intent(in) :: ixImin1, ixImin2, ixImin3, ixImax1, ixImax2, ixImax3
    integer, intent(in) :: ixOmin1, ixOmin2, ixOmin3, ixOmax1, ixOmax2, ixOmax3
    integer, intent(in) :: iB
    double precision, intent(in) :: qt
    double precision, intent(inout),&
         dimension(ixImin1:ixImax1, ixImin2:ixImax2, ixImin3:ixImax3, 1:nw_phys) :: w
    double precision, intent(in),&
         dimension(ixImin1:ixImax1, ixImin2:ixImax2, ixImin3:ixImax3, 1:ndim) :: x
    integer              :: ix1, ix2, ix3, iw
    double precision     :: rhoh, E_th, E, b_sqr

    select case(iB)
    case(5)
       ! select the first grid-internal layer above the boundary
#ifndef _OPENMP
      ! Vector-level parallelization within a subroutine does not work with OpenMP.
      ! TBD if this routine can be made seq (i.e. cell-based).
      ${GPU_LOOP_VECTOR("collapse(3)")}$
#endif
      do ix3 = ixOmin3, ixOmax3
         do ix2 = ixOmin2, ixOmax2
            do ix1 = ixOmin1, ixOmax1
               if (x(ix1,ix2,ix3,1)**2+x(ix1,ix2,ix3,2)**2 < rjet**2) then
                                 
                  if(srmhd_eos) then
                     E_th = pjet * 1/(srmhd_gamma-1)
                     E    = E_th + dsqrt(E_th**2+rhojet**2)
                     ! writing rho/E on purpose, for numerics 
                     rhoh = 0.5_dp * ((srmhd_gamma+1.0_dp) * E &
                          - (srmhd_gamma-1) * rhojet * (rhojet/E))
                  else
                     rhoh = rhojet + srmhd_gamma/(srmhd_gamma-1)*pjet
                  end if

                  b_sqr = w(ix1,ix2,ix3,iw_mag(1))**2 &
                        + w(ix1,ix2,ix3,iw_mag(2))**2 &
                        + w(ix1,ix2,ix3,iw_mag(3))**2
                  
                  w(ix1,ix2,ix3,iw_rho)    = rhojet*lfacjet
                  w(ix1,ix2,ix3,iw_mom(1)) = 0.0d0
                  w(ix1,ix2,ix3,iw_mom(2)) = 0.0d0
                  w(ix1,ix2,ix3,iw_mom(3)) = rhoh * lfacjet**2 * vjet + b_sqr*vjet - w(ix1,ix2,ix3,iw_mag(3)) * vjet**2 
                  w(ix1,ix2,ix3,iw_e)      = &
                       rhoh * lfacjet**2 - pjet - rhojet * lfacjet &
                       + 0.5_dp*(b_sqr*(1 + vjet**2) - (w(ix1,ix2,ix3,iw_mag(3)) * vjet)**2)

                  w(ix1,ix2,ix3,iw_mag(1)) =   0
                  w(ix1,ix2,ix3,iw_mag(2)) =   0
                  w(ix1,ix2,ix3,iw_mag(3)) = 0

               else
! Calling primitive here is somehow broken.  To fix!!!!                  
!                  call to_primitive(u)
                  w(ix1,ix2,ix3,iw_rho)    =   w(ix1,ix2,2*ixOmax3-ix3+1,iw_rho)
                  w(ix1,ix2,ix3,iw_mom(1)) =   w(ix1,ix2,2*ixOmax3-ix3+1,iw_mom(1))
                  w(ix1,ix2,ix3,iw_mom(2)) =   w(ix1,ix2,2*ixOmax3-ix3+1,iw_mom(2))
                  w(ix1,ix2,ix3,iw_mom(3)) = - w(ix1,ix2,2*ixOmax3-ix3+1,iw_mom(3))
                  w(ix1,ix2,ix3,iw_e)      =   w(ix1,ix2,2*ixOmax3-ix3+1,iw_e)
                  w(ix1,ix2,ix3,iw_mag(1)) =   w(ix1,ix2,2*ixOmax3-ix3+1,iw_mag(1))
                  w(ix1,ix2,ix3,iw_mag(2)) =   w(ix1,ix2,2*ixOmax3-ix3+1,iw_mag(2))
                  w(ix1,ix2,ix3,iw_mag(3)) = - w(ix1,ix2,2*ixOmax3-ix3+1,iw_mag(3))

               endif
            end do
         end do
      end do

#:if defined(N_TRACER)
      if(srmhd_n_tracer>0)then
         do iw=1,srmhd_n_tracer
#ifndef _OPENMP
            ${GPU_LOOP_VECTOR("collapse(3)")}$
#endif
            do ix3 = ixOmin3, ixOmax3
               do ix2 = ixOmin2, ixOmax2
                  do ix1 = ixOmin1, ixOmax1
                     if (x(ix1,ix2,ix3,1)**2+x(ix1,ix2,ix3,2)**2<rjet**2) then
                        w(ix1,ix2,ix3,tracer(iw))    = w(ix1,ix2,ix3,iw_rho)
                     else
                        w(ix1,ix2,ix3,tracer(iw))    = 0.0d0
                     endif
                  end do
               end do
            end do
         enddo
      endif

  #:endif

      ! Curently broken, cannot call anything here.      
      ! ! switch to conserved in ghost cells
      ! ${GPU_LOOP_VECTOR("collapse(3)")}$
      ! do ix3 = ixOmin3, ixOmax3
      !    do ix2 = ixOmin2, ixOmax2
      !       do ix1 = ixOmin1, ixOmax1
      !          do iw = 1, nw_phys
      !             u(iw) = w(ix1,ix2,ix3,iw)
      !          end do
      !          call to_conservative(u)
      !          do iw = 1, nw_phys
      !             w(ix1,ix2,ix3,iw) = u(iw)
      !          end do
      !       end do
      !    end do
      ! end do

   end select

  end subroutine specialbound_usr


   subroutine specialvar_output( &
      ixImin1, ixImin2, ixImin3, ixImax1, ixImax2, ixImax3,&
      ixOmin1, ixOmin2, ixOmin3, ixOmax1, ixOmax2, ixOmax3,&
      w, x, normconv)
   ! this subroutine can be used in convert, to add auxiliary variables to the
   ! converted output file, for further analysis using tecplot, paraview, ....
   ! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
   !
   ! the array normconv can be filled in the (nw+1:nw+nwauxio) range with
   ! corresponding normalization values (default value 1)
    use mod_global_parameters
    use mod_physics
    implicit none
    integer, intent(in) :: ixImin1, ixImin2, ixImin3, ixImax1, ixImax2, ixImax3
    integer, intent(in) :: ixOmin1, ixOmin2, ixOmin3, ixOmax1, ixOmax2, ixOmax3
    double precision,&
      dimension(ixImin1:ixImax1, ixImin2:ixImax2, ixImin3:ixImax3, 1:nw+nwauxio) :: w
    double precision, intent(in),&
      dimension(ixImin1:ixImax1, ixImin2:ixImax2, ixImin3:ixImax3, 1:ndim) :: x
    double precision                   :: normconv(0:nw+nwauxio)

    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,nw+1)= &
       dlog10(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,d_) &
             /w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,lfac_))

   end subroutine specialvar_output

   subroutine specialvarnames_output(varnames)

   ! newly added variables
   character(len=*) :: varnames
   varnames='log10rho'

   end subroutine specialvarnames_output


  subroutine usr_refine_grid(level,qt,w,x,refineflag,coarsenflag,norefineflag,nocoarsenflag)
    use mod_global_parameters
    ${GPU_ROUTINE_SEQ()}$

    ! Enforce additional refinement or coarsening
    ! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.

    integer, intent(in)             :: level
    double precision, intent(in)    :: qt
    double precision, intent(in)    :: x(1:3)
    double precision, intent(in)    :: w(1:nw_phys)
    logical, intent(inout) :: refineflag, coarsenflag, norefineflag, nocoarsenflag

    if (x(3) < 1.1d0*zjet .and. (x(1)**2 + x(2)**2) < rjet**2) then
      refineflag = .true.
      nocoarsenflag = .true.
    end if

    ! Never force coarsen or block automatic refine
    coarsenflag = .false.
    norefineflag = .false.

  end subroutine usr_refine_grid


  !> Density and pressure floor, called every timestep via
  !> usr_process_adv_global, straight on the device: no host/device traffic,
  !> and rho and p are floored through a full to_primitive/to_conservative
  !> round trip - reusing mod_physics's own copies (public via the plain
  !> `use mod_physics` above, since mod_physics instantiates the srmhd macros
  !> itself) rather than duplicating or re-instantiating them here - so xi_
  !> and lfac_ stay consistent with the floored state rather than silently
  !> disagreeing with it.
  !>
  !> The thresholds are srmhd's own srmhd_small_density and srmhd_small_pressure
  !> from &srmhd_list, not case-local parameters, so that this hook, the
  !> pressure bracket in con2prim and the floor fix_prim_state applies to
  !> reconstructed face states all agree on what "too small" means. A no-op
  !> unless at least one of them is set positive.
  !>
  !> Covers the whole block (ixG), ghost cells included, rather than just
  !> the mesh. getbc only refreshes a *destination* state at the end of each
  !> RK substage, never the source state a substage starts from, and nothing
  !> in the main loop calls getbc between this hook and the next timestep's
  !> first substage - so a mesh-only floor would leave stale, unfloored
  !> ghost values feeding straight into the next flux computation at every
  !> block boundary, defeating the point of flooring. Since the floor is a
  !> pointwise max(value, floor), flooring a block's own ghost cell and
  !> flooring the neighbour's interior cell it mirrors give bit-identical
  !> results, so covering ixG keeps every ghost consistent with its owning
  !> interior without any ghost exchange at all.
  subroutine apply_floor(iit, qt)
    use mod_global_parameters
    implicit none
    integer, intent(in)          :: iit
    double precision, intent(in) :: qt

    integer :: iigrid, igrid, ix1, ix2, ix3
    double precision :: u(1:nw_phys)

    if (srmhd_small_density <= 0.0d0 .and. srmhd_small_pressure <= 0.0d0) return

    ${GPU_PARALLEL_LOOP_GANG("private(igrid)")}$
    do iigrid = 1, igridstail_active
       igrid = igrids_active(iigrid)

       ${GPU_LOOP_VECTOR("collapse(3) private(u)")}$
       do ix3 = ixGlo3, ixGhi3
          do ix2 = ixGlo2, ixGhi2
             do ix1 = ixGlo1, ixGhi1

                u = bg(1)%w(ix1, ix2, ix3, 1:nw_phys, igrid)
                call to_primitive(u)

                u(iw_rho) = max(u(iw_rho), srmhd_small_density)
                u(iw_e)   = max(u(iw_e),   srmhd_small_pressure)

                call to_conservative(u)
                bg(1)%w(ix1, ix2, ix3, 1:nw_phys, igrid) = u

             end do
          end do
       end do
    end do

  end subroutine apply_floor


end module mod_usr
