module mod_hd_parameters
  use mod_physics
  use mod_constants
  implicit none
  public

  !-------------------------------------------------------------------!
  ! Parameters for global settings
  !-------------------------------------------------------------------!
  logical                      :: evolve_hydro = .True.

  !> if used fluid cmax to calculate dt, otherwise use speed of light
  logical                      :: fluid_cmax   = .True.

  
  !-------------------------------------------------------------------!
  ! Parameters from AMRVAC
  !-------------------------------------------------------------------!  
  
  integer                      :: boundspeed = 1

  !-------------------------------------------------------------------!
  ! Parameters for ideal gas relations
  !-------------------------------------------------------------------!
  
  double precision             :: hd_gamma = 5.0d0/3.0d0
  double precision             :: hd_gamma_ad = 1.0d0

  !> The smallest allowed density
  double precision, public     :: small_rho_fac = 0.1d0
  double precision, public     :: small_rho_thr = 0.0d0
  double precision, public     :: small_rho     = smalldouble
  !> The smallest allowed e = rho * eps
  double precision, public     :: small_e       = smalldouble 
  !> The smallest allowed press
  double precision, public     :: small_press   = smalldouble 
  double precision             :: atmo_gamma    = 2.0d0
  double precision             :: atmo_adiab    = 1.0d2

  contains

  !> This subroutine fix the abnormal values in primitive/conserved variables
  subroutine hd_handle_small_values(ps_in, ixI^L, ixO^L, subname)
    use mod_global_parameters
    use mod_small_values
    use mod_geometry

    type(state), intent(inout)      :: ps_in
    integer, intent(in)             :: ixI^L,ixO^L
    character(len=*), intent(in)    :: subname

    integer                         :: ix^D
    double precision                :: radii_pt, x_tmp(1:^ND)
    double precision                :: eps_hat

    associate( w => ps_in%w, wextra => ps_in%wextra, x => ps_in%mesh%x, &
               m => ps_in%metric%vars )

    if (ps_in%is_prim) then
       ! if the input is primitive

       ! avoid coordinate singularities
       ! fixme: this might depends on different bc, but in general this should work
       if ( coordinate /= cartesian ) then
          where ( dabs(x(ixO^S,1)) < smalldouble ) 
             w(ixO^S, W_vel(1)) = 0.0d0
          end where
       end if
   
       select case (small_values_method)
       case ("replace")
          ! check the prim variables one by one
          {do ix^D = ixO^LIM^D \}
             x_tmp = x(ix^D,:)
             radii_pt = get_radii_pt(x_tmp)
             if ( w(ix^D, rho_) <= small_rho_thr) then
                ! atmosphere handling
                w(ix^D, rho_) = small_rho
                w(ix^D, W_vel(:)) = 0.0d0
                w(ix^D, e_) = atmo_adiab * w(ix^D, rho_)**atmo_gamma / ( atmo_gamma - 1.0d0 )
             else
                w(ix^D, e_) = max(0.0d0, w(ix^D, e_))
             end if
          {end do^D&\}
       case default
          ! nothing to do here
          !call small_values_error(w, x, ixI^L, ixO^L, flag, subname)
          return
       end select

    else
       ! if the input is conserved

       ! avoid coordinate singularities
       ! fixme: this might depends on different bc, but in general this should work
       if ( coordinate /= cartesian ) then
          where ( dabs(x(ixO^S,1)) < smalldouble ) 
             w(ixO^S, mom(1)) = 0.0d0
          end where
       end if
   
       select case (small_values_method)
       case ("replace")
          ! check the prim variables one by one
          {do ix^D = ixO^LIM^D \}
             x_tmp = x(ix^D,:)
             radii_pt = get_radii_pt(x_tmp)
             if ( w(ix^D, D_) <= small_rho_thr) then
                ! atmosphere handling
                w(ix^D, D_) = small_rho
                w(ix^D, mom(:)) = 0.0d0
                w(ix^D, tau_) = atmo_adiab * w(ix^D, D_)**atmo_gamma / ( atmo_gamma - 1.0d0 )
             end if
          {end do^D&\}
       case default
          ! nothing to do here
          !call small_values_error(w, x, ixI^L, ixO^L, flag, subname)
          return
       end select

    end if

    end associate
  end subroutine hd_handle_small_values
  
  !> get some useful variables from primitive
  subroutine hd_get_intermediate_variables(ps_in, ixI^L, ixO^L, press, cs2 )
     use mod_global_parameters
     use mod_geometry
     type(state), intent(in)                 :: ps_in
     integer, intent(in)                     :: ixI^L, ixO^L

     double precision, intent(out), optional :: cs2(ixI^S)   ! speed of sound squared
     double precision, intent(out), optional :: press(ixI^S) ! pressure

     integer                                 :: ix^D
     double precision                        :: tmp_prs(ixI^S), tmp_cs2(ixI^S)

     if (.not.ps_in%is_prim) call mpistop("get_inter err: the input has to be prim here")

     associate( prim => ps_in%w, wextra => ps_in%wextra, x => ps_in%mesh%x )

     {do ix^D = ixO^LIM^D \}
        tmp_prs(ix^D) = ( hd_gamma - 1.0d0 ) * prim(ix^D,e_)
        tmp_cs2(ix^D) = hd_gamma * tmp_prs(ix^D) / ( prim(ix^D,rho_) )
     {end do^D&\}
     if (present(press)) press(ixO^S) = tmp_prs(ixO^S)
     if (present(cs2)) cs2(ixO^S) = tmp_cs2(ixO^S)

     end associate
  end subroutine hd_get_intermediate_variables

end module mod_hd_parameters
