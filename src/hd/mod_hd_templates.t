#:def addsource_local()
subroutine addsource_local(qdt, dtfactor, qtC, wCT, wCTprim, qt, wnew, x, qsourcesplit)
  !$acc routine seq
#:if defined('GRAVITY')
  use mod_usr, only: gravity_field
#:endif    
  real(dp), intent(in)     :: qdt, dtfactor, qtC, qt
  real(dp), intent(in)     :: wCT(nw_euler), wCTprim(nw_euler)
  real(dp), intent(in)     :: x(1:ndim)
  real(dp), intent(inout)  :: wnew(nw_euler)
  logical, intent(in)      :: qsourcesplit
  ! .. local ..
  integer                  :: idim
  real(dp)                 :: field

#:if defined('GRAVITY')
  do idim = 1, ndim
     field = gravity_field(wCT, x, idim)
     wnew(iw_mom(idim)) = wnew(iw_mom(idim)) + qdt * field * wCT(iw_rho)
     wnew(iw_e)         = wnew(iw_e) + qdt * field * wCT(iw_mom(idim))
  end do
#:endif  

end subroutine addsource_local
#:enddef

#:def to_primitive()
pure subroutine to_primitive(u)
  !$acc routine seq
  real(dp), intent(inout) :: u(nw_euler)

  {^D&
       u(iw_mom(^D)) = u(iw_mom(^D))/u(iw_rho)
  \}

  u(iw_e) = (hd_gamma-1.0_dp) * (u(iw_e) - &
       0.5_dp * u(iw_rho) * sum(u(iw_mom(1:ndim))**2) )

end subroutine to_primitive
#:enddef

#:def to_conservative()  
pure subroutine to_conservative(u)
  !$acc routine seq
  real(dp), intent(inout) :: u(nw_euler)
  real(dp)                :: inv_gamma_m1

  inv_gamma_m1 = 1.0d0/(hd_gamma - 1.0_dp)

  ! Compute energy from pressure and kinetic energy
  u(iw_e) = u(iw_e) * inv_gamma_m1 + &
       0.5_dp * u(iw_rho) * sum(u(iw_mom(1:ndim))**2)

  ! Compute momentum from density and velocity components
  {^D&
       u(iw_mom(^D)) = u(iw_rho) * u(iw_mom(^D))
  \}

end subroutine to_conservative
#:enddef

#:def get_flux()
subroutine get_flux(u, flux_dim, flux)
  !$acc routine seq
  real(dp), intent(in)  :: u(nw_euler)
  integer, intent(in)   :: flux_dim
  real(dp), intent(out) :: flux(nw_euler)
  real(dp)              :: inv_gamma_m1

  inv_gamma_m1 = 1.0d0/(hd_gamma - 1.0_dp)

  ! Density flux
  flux(iw_rho) = u(iw_rho) * u(iw_mom(flux_dim))

  ! Momentum flux with pressure term
  {^D&
       flux(iw_mom(^D)) = u(iw_rho) * u(iw_mom(^D)) * u(iw_mom(flux_dim))
  \}
  flux(iw_mom(flux_dim)) = flux(iw_mom(flux_dim)) + u(iw_e)

  ! Energy flux
  flux(iw_e) = u(iw_mom(flux_dim)) * (u(iw_e) * inv_gamma_m1 + &
       0.5_dp * u(iw_rho) * sum(u(iw_mom(1:ndim))**2) + u(iw_e))

end subroutine get_flux
#:enddef

#:def get_cmax()  
pure real(dp) function get_cmax(u, flux_dim) result(wC)
  !$acc routine seq
  real(dp), intent(in)  :: u(nw_euler)
  integer, intent(in)   :: flux_dim

  wC = sqrt(hd_gamma * u(iw_e) / u(iw_rho)) + abs(u(iw_mom(flux_dim)))

end function get_cmax
#:enddef  
