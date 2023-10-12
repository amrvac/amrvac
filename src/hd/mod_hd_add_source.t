module mod_hd_add_source
  use mod_physics
  use mod_hd_parameters

  implicit none
  private

  ! Public methods
  public :: hd_add_source_init

contains

  !> Initialize the module
  subroutine hd_add_source_init()
    use mod_global_parameters
    use mod_geometry
    if (coordinate /= cartesian) then
       ! we need geom source terms
       phys_add_source_geom     => hd_add_source_geom
    end if
    phys_add_source          => hd_add_source
  end subroutine hd_add_source_init

  !> Add gravitational source terms to w
  ! fixme: add these source terms
  subroutine hd_add_source(qdt, sCT, sold, ixI^L, ixO^L, cons, qsourcesplit, active)
    use mod_global_parameters
    use mod_geometry
    double precision, intent(in)    :: qdt
    type(state), intent(in)         :: sCT, sold
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: cons(ixI^S, 1:nwflux)
    logical, intent(in)             :: qsourcesplit
    logical, intent(inout)          :: active !< Needs to be set to true when active

    double precision                :: Delta_t(ixI^S) ! qdt * alp * psi6
    double precision                :: add_source(ixI^S,1:nwflux)
    double precision                :: vec_tmp1(ixI^S,1:3)
    double precision                :: rhoh(ixI^S)  ! rho * enthalpy
    double precision                :: press(ixI^S) ! pressure
    double precision                :: lfac(ixI^S)  ! Lorentz factor


    
  end subroutine hd_add_source

  ! Unchanged from GRHD because the terms should be the same
  !> Add geometrical source terms to w
  subroutine hd_add_source_geom(qdt, sCT, sCTp, ixI^L, ixO^L, cons)
    use mod_global_parameters
    use mod_geometry
    double precision, intent(in)    :: qdt
    type(state), intent(in)         :: sCT  ! physical state at preivous substep
    type(state), intent(in)         :: sCTp ! physical state at preivous substep, in primitive
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: cons(ixI^S, 1:nwflux)

    double precision                :: source(ixI^S)
    double precision                :: fluxCT(ixI^S,1:3,1:3)
    integer                         :: idir,jdir,kdir
    double precision                :: press(ixI^S)

    if ( .not. evolve_hydro ) return

    associate( x => sCT%mesh%x, chris_hat => sCT%mesh%christoffel, &
               m => sCT%metric%vars,                               &
               consCT => sCT%w, primCT => sCTp%w )
    
    call hd_get_intermediate_variables(sCTp, ixI^L, ixO^L, press=press)

    ! Only Momentum flux f^i_j is needed only
    fluxCT = 0.0d0
    do jdir = 1, ndir; do idir = 1, ndir
      fluxCT(ixO^S, idir, jdir) = primCT(ixO^S, W_vel(idir)) * consCT(ixO^S, mom(jdir))
    end do; end do
    ! and the pressure terms
    do jdir = 1, 3
       fluxCT(ixO^S, jdir,jdir) = ( fluxCT(ixO^S, jdir, jdir) + press(ixO^S) )
    end do

    select case (coordinate)
    case (cylindrical)
       ! geo_source[r] = Gamma^phi_{r phi}*f^phi_phi
       source(ixO^S) = chris_hat(ixO^S,2) * fluxCT(ixO^S,phi_,phi_)
       cons(ixO^S, mom(r_)) = cons(ixO^S, mom(r_)) + qdt * source(ixO^S)

    case (spherical)
       ! geo_source[r] = Gamma^2_{12}*f^2_2 + Gamma^3_{13}*f^3_3
       source(ixO^S) = chris_hat(ixO^S,1) * fluxCT(ixO^S,2,2) & 
                      +chris_hat(ixO^S,1) * fluxCT(ixO^S,3,3) 

       cons(ixO^S, mom(r_)) = cons(ixO^S, mom(r_)) + qdt * source(ixO^S)

       {^NOONED
       ! geo_source[theta] = Gamma^1_{22}*f^2_1 + Gamma^2_{21}*f^1_2 + Gamma^3_{23}*f^3_3
       ! where 
       ! Gamma^1_{22} = -r
       ! Gamma^2_{21} = 1/r
       ! Gamma^3_{23} = cot theta , only this term is related to theta
       source(ixO^S) = chris_hat(ixO^S,2) * fluxCT(ixO^S,2,1) & 
                      +chris_hat(ixO^S,1) * fluxCT(ixO^S,1,2) &
                      +chris_hat(ixO^S,5) * fluxCT(ixO^S,3,3) 

       cons(ixO^S, mom(theta_)) = cons(ixO^S, mom(theta_)) + qdt * source(ixO^S)
       ! geo_source[phi] = 0
       }
    end select
    end associate
  end subroutine hd_add_source_geom

end module mod_hd_add_source
