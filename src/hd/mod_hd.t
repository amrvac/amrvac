!> Hydrodynamics physics module
module mod_hd
  use mod_physics
  use mod_hd_parameters

  implicit none
  private

  ! Public methods
  public :: hd_activate

contains

  subroutine hd_activate()
    use mod_hd_ppm
    call hd_init()
    call hd_ppm_init() ! fixme: do we really need the ppm module?
  end subroutine hd_activate

  !> Read this module's parameters from a file
  ! fixme: update the hd_list
  subroutine hd_read_params(files)
    use mod_global_parameters
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /hd_list/ evolve_hydro, hd_gamma, hd_gamma_ad,&
                       fluid_cmax, atmo_gamma, atmo_adiab

    do n = 1, size(files)
       rewind(unitpar)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, hd_list, end=111)
111    close(unitpar)
    end do

    if (mype == 0) then
       if ( .not. evolve_hydro ) write(*,*) "Hydro is not evolving"
    end if

  end subroutine hd_read_params

  ! fixme: remove these subroutine
  !> Write this module's parameters to a snapsoht
  subroutine hd_write_info(fh)
    use mod_global_parameters
    integer, intent(in)                 :: fh
    integer, parameter                  :: n_par = 1
    double precision                    :: values(n_par)
    character(len=name_len)             :: names(n_par)
    integer, dimension(MPI_STATUS_SIZE) :: st
    integer                             :: er
    call MPI_FILE_WRITE(fh, n_par, 1, MPI_INTEGER, st, er)
    names(1) = "gamma"
    values(1) = 1.0d0!hd_gamma
    call MPI_FILE_WRITE(fh, values, n_par, MPI_DOUBLE_PRECISION, st, er)
    call MPI_FILE_WRITE(fh, names, n_par * name_len, MPI_CHARACTER, st, er)
  end subroutine hd_write_info

  !> Initialize the module
  subroutine hd_init()
    use mod_global_parameters
    use mod_hd_flux
    use mod_hd_add_source
    use mod_hd_convert

    integer :: idir

    call hd_read_params(par_files)

    physics_type = "hd"
  
    ! Determine primitive variables that need to be reconstructed
    rho_ = var_set_fluxvar('D','rho', var_type=hydro_var)
    
    ! W_vel corresponds here only to the velocity (maybe create a new variable for this)
    allocate(W_vel(ndir))
    do idir = 1, ndir
       W_vel(idir) = var_set_fluxvar('mom','v',idir, var_type=hydro_var)
    end do
    
    ! e is here the INTERNAL energy ( e = p/(gamma-1) ) and tau the total energy
    e_ = var_set_fluxvar('tau','e', var_type=hydro_var)

    ! the corresponding cons
    D_ = rho_
    allocate(mom(ndir))
    do idir = 1, ndir
       mom(idir) = W_vel(idir)
    end do
    tau_ = e_

    phys_get_dt              => hd_get_dt
    phys_get_cmax            => hd_get_cmax
    phys_check_params        => hd_check_params
    phys_check_prim          => hd_check_prim
    phys_write_info          => hd_write_info
    phys_handle_small_values => hd_handle_small_values

    call hd_convert_init()
    call hd_flux_init()
    call hd_add_source_init()

    ! Whether diagonal ghost cells are required for the physics
    phys_req_diagonal = .false.

    ! Check whether custom flux types have been defined
    if (.not. allocated(flux_type)) then
       allocate(flux_type(ndir, nwflux))
    else if (any(shape(flux_type) /= [ndir, nwflux])) then
       call mpistop("phys_check error: flux_type has wrong shape")
    end if

    flux_type = flux_default

    if ( .not.evolve_hydro ) then
       flux_type(:, D_)          = flux_nul
       flux_type(:, mom(1:ndir)) = flux_nul
       flux_type(:, tau_)        = flux_nul
    end if

    allocate(start_indices(1:number_species))
    allocate(stop_indices(1:number_species))

    start_indices(1) = D_
    stop_indices(1)  = max( mom(ndir), tau_)

  end subroutine hd_init

  subroutine hd_check_params
    use mod_global_parameters
    use mod_geometry, only: coordinate, spherical
    use mod_limiter

    if ( any(type_limiter(:) == limiter_ppm) ) phys_req_diagonal = .true.

  end subroutine hd_check_params

  !> Calculate cmax_idim within ixO^L (cmax = csound + fluid velocity)
  subroutine hd_get_cmax(ps_in, ixI^L, ixO^L, idim, cmax)
    use mod_global_parameters
    use mod_geometry
    type(state), intent(in)         :: ps_in
    integer, intent(in)             :: ixI^L, ixO^L, idim
    double precision, intent(inout) :: cmax(ixI^S)
    integer                         :: iw, ix^D
    double precision                :: csound(ixI^S)

    if (.not.ps_in%is_prim) call mpistop('hd_get_cmax expect prim input')
   
    call hd_get_intermediate_variables(ps_in,ixI^L,ixO^L,cs2=csound)
    csound(ixO^S) = dsqrt(csound(ixO^S))

    cmax(ixO^S) = dabs(ps_in%w(ixO^S, W_vel(idim)))+csound(ixO^S)

  end subroutine hd_get_cmax
  
  ! Unchanged from GRHD
  subroutine hd_get_dt(ps_in, ixI^L, ixO^L, dtnew)
    use mod_global_parameters

    type(state), intent(in)         :: ps_in
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: dtnew

    dtnew = bigdouble
  end subroutine hd_get_dt

  ! fixme: maybe we dont need this subroutine anymore
  ! Unchanged from GRHD
  !> Returns 0 in argument flag where values are ok
  subroutine hd_check_prim(ps_in, ixI^L, ixO^L, flag)
    use mod_global_parameters
    type(state), intent(in)      :: ps_in
    integer, intent(in)          :: ixI^L, ixO^L
    integer, intent(inout)       :: flag(ixI^S)
    associate( w => ps_in%w)
    flag(ixO^S) = 0
    ! this subroutine is used only in MP5, so we care about rho and eps only.
    where( (w(ixO^S, rho_) < smalldouble) .or. isnan(w(ixO^S, rho_)) ) 
       flag(ixO^S) = rho_
    end where
    end associate
  end subroutine hd_check_prim

end module mod_hd
