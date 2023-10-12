module mod_hd_convert
  use mod_physics
  use mod_hd_parameters

  implicit none
  private

  ! Public methods
  public :: hd_convert_init

contains

  !> Initialize the module
  subroutine hd_convert_init()

    phys_to_conserved    => hd_to_conserved
    phys_to_primitive    => hd_to_primitive
 
  end subroutine hd_convert_init

  !> Transform primitive variables into conservative ones
  subroutine hd_to_conserved(ixI^L, ixO^L, ps_in)
    use mod_global_parameters

    integer, intent(in)        :: ixI^L, ixO^L
    type(state), intent(inout) :: ps_in

    integer                    :: idir, iw

    associate( w => ps_in%w, wextra => ps_in%wextra )

    if (.not.ps_in%is_prim) call mpistop("The input is not prim")

    if ( fix_small_values ) then
       call hd_handle_small_values(ps_in, ixI^L, ixO^L, 'hd_to_conserved')
    end if

    !w(ixO^S, D_) is unchanged
    
    ! tau = e + kinetic energy
    w(ixO^S, tau_) = w(ixO^S, e_) + 0.5d0 * sum(w(ixO^S, W_vel(:))**2, dim=ndim+1) * w(ixO^S, rho_)
    
    do idir = 1, ndir
       w(ixO^S, mom(idir)) = w(ixO^S, W_vel(idir)) * w(ixO^S, rho_)
    end do
 
    ps_in%is_prim = .false.
    end associate
  end subroutine hd_to_conserved
  
  
  !> Transform conservative variables into primitve ones
  subroutine hd_to_primitive(ixI^L, ixO^L, ps_in)
    use mod_global_parameters

    integer, intent(in)        :: ixI^L, ixO^L
    type(state), intent(inout) :: ps_in

    integer                    :: idir, iw

    associate( w => ps_in%w, wextra => ps_in%wextra )

    if (ps_in%is_prim) call mpistop("The input is not cons")

    if ( fix_small_values ) then
       call hd_handle_small_values(ps_in, ixI^L, ixO^L, 'hd_to_primitive')
    end if

    !w(ixO^S, rho_) is unchanged
    
    do idir = 1, ndir
       w(ixO^S, W_vel(idir)) = w(ixO^S, mom(idir)) / w(ixO^S, rho_)
    end do
    
    ! e = tau - kinetic energy
    w(ixO^S, e_) = w(ixO^S, tau_) - 0.5d0 * sum(w(ixO^S, W_vel(:))**2, dim=ndim+1) * w(ixO^S, rho_)
 
    ps_in%is_prim = .true.
    end associate
  end subroutine hd_to_primitive


  
end module mod_hd_convert
