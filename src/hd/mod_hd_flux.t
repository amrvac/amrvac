module mod_hd_flux
  use mod_physics
  use mod_hd_parameters

  implicit none
  private

  ! Public methods
  public :: hd_flux_init

contains
  !> Initialize the module
  subroutine hd_flux_init()
    phys_get_cbounds         => hd_get_cbounds
    phys_get_flux            => hd_get_flux
  end subroutine hd_flux_init
  
  ! Taken and adapted from MPI-AMRVAC
  !> Calculate cmax_idim = csound + abs(v_idim) within ixO^L
  subroutine hd_get_cbounds(sL, sR, sLp, sRp, ixI^L, ixO^L, idim, cmax, cmin)
    use mod_global_parameters
    use mod_geometry

    ! left and right status
    type(state), intent(in)                   :: sL, sR, sLp, sRp
    integer, intent(in)                       :: ixI^L, ixO^L, idim
    double precision, intent(inout)           :: cmax(ixI^S, 1:number_species)
    double precision, intent(inout), optional :: cmin(ixI^S, 1:number_species)

    double precision :: wmean(ixI^S,nw)
    double precision, dimension(ixI^S) :: umean, dmean, csoundL, csoundR, tmp1,tmp2,tmp3
    integer :: ix^D
    
    if (.not. sLp%is_prim .or. .not. sRp%is_prim) call mpistop('hd_get_cbounds expect prim input')
    
    associate( wLp => sLp%w, wRp => sRp%w, &
               wLC => sL%w, wRC => sR%w)

    select case(boundspeed)
    case (1)
      ! This implements formula (10.52) from "Riemann Solvers and Numerical
      ! Methods for Fluid Dynamics" by Toro.

      tmp1(ixO^S)=dsqrt(wLp(ixO^S,rho_))
      tmp2(ixO^S)=dsqrt(wRp(ixO^S,rho_))
      tmp3(ixO^S)=1.d0/(dsqrt(wLp(ixO^S,rho_))+dsqrt(wRp(ixO^S,rho_)))
      umean(ixO^S)=(wLp(ixO^S,mom(idim))*tmp1(ixO^S)+wRp(ixO^S,mom(idim))*tmp2(ixO^S))*tmp3(ixO^S)

      call hd_get_intermediate_variables(sLp,ixI^L,ixO^L,cs2=csoundL)
      call hd_get_intermediate_variables(sRp,ixI^L,ixO^L,cs2=csoundR)

      dmean(ixO^S) = (tmp1(ixO^S)*csoundL(ixO^S)+tmp2(ixO^S)*csoundR(ixO^S)) * &
           tmp3(ixO^S) + 0.5d0*tmp1(ixO^S)*tmp2(ixO^S)*tmp3(ixO^S)**2 * &
           (wRp(ixO^S,mom(idim))-wLp(ixO^S,mom(idim)))**2

      dmean(ixO^S)=dsqrt(dmean(ixO^S))
      if(present(cmin)) then
        cmin(ixO^S,1)=umean(ixO^S)-dmean(ixO^S)
        cmax(ixO^S,1)=umean(ixO^S)+dmean(ixO^S)
        
      else
        cmax(ixO^S,1)=dabs(umean(ixO^S))+dmean(ixO^S)
      end if

    case (2)

      call mpistop('Not implemented yet')

    case (3)
      ! Miyoshi 2005 JCP 208, 315 equation (67)
      call hd_get_intermediate_variables(sLp,ixI^L,ixO^L,cs2=csoundL)
      call hd_get_intermediate_variables(sRp,ixI^L,ixO^L,cs2=csoundR)
      
      csoundL(ixO^S)=max(dsqrt(csoundL(ixO^S)),dsqrt(csoundR(ixO^S)))
      if(present(cmin)) then
        cmin(ixO^S,1)=min(wLp(ixO^S,mom(idim)),wRp(ixO^S,mom(idim)))-csoundL(ixO^S)
        cmax(ixO^S,1)=max(wLp(ixO^S,mom(idim)),wRp(ixO^S,mom(idim)))+csoundL(ixO^S)
        
      else
        cmax(ixO^S,1)=max(wLp(ixO^S,mom(idim)),wRp(ixO^S,mom(idim)))+csoundL(ixO^S)
      end if
      
    end select
    
    end associate

  end subroutine hd_get_cbounds

  ! Calculate flux f_idim[iw]
  subroutine hd_get_flux(s_cons, s_prim, ixI^L, ixO^L, idim, f)
    use mod_global_parameters
    use mod_geometry

    type(state), intent(in)         :: s_cons, s_prim
    integer, intent(in)             :: ixI^L, ixO^L, idim
    double precision, intent(out)   :: f(ixI^S, 1:nwflux)

    double precision                :: press(ixI^S) 
    integer                         :: idir

    associate( cons => s_cons%w, prim => s_prim%w, &
               x => s_cons%mesh%x, m => s_cons%metric%vars )

    if ( evolve_hydro ) then 
       call hd_get_intermediate_variables(s_prim, ixI^L, ixO^L, press=press)
       
       ! Density flux
       f(ixO^S, D_) = prim(ixO^S, W_vel(idim)) * cons(ixO^S, D_)
   
       ! Momentum flux
       do idir = 1, ndir
         f(ixO^S, mom(idir)) = prim(ixO^S, W_vel(idir))  * cons(ixO^S, mom(idir))
       end do
       ! add pressure term
       f(ixO^S, mom(idim)) = f(ixO^S, mom(idim)) +  press(ixO^S)
   
       ! Energy flux 
       f(ixO^S, tau_) = prim(ixO^S, W_vel(idim)) * cons(ixO^S, tau_) &
                         + press(ixO^S)  * prim(ixO^S,W_vel(idim))
    else
       f(ixO^S, D_)          = 0.0d0
       f(ixO^S, mom(1:ndir)) = 0.0d0
       f(ixO^S, tau_)        = 0.0d0

    end if

    end associate
  end subroutine hd_get_flux

end module mod_hd_flux
