module mod_rho

  implicit none
  private

  ! Public methods
  public :: rho_activate

contains

  subroutine rho_activate(par_files)
    use mod_rho_phys, only: rho_phys_init
    use mod_rho_roe, only: rho_roe_init

    character(len=*), intent(in) :: par_files(:)

    call rho_phys_init(par_files)
    call rho_roe_init()

  end subroutine rho_activate

end module mod_rho
