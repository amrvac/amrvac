module mod_rho

  implicit none
  public

contains

  subroutine rho_activate(par_files)
    character(len=*), intent(in) :: par_files(:)

    call rho_phys_init(par_files)
    call rho_roe_init()

  end subroutine rho_activate

end module mod_rho
