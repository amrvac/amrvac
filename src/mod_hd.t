module mod_hd

  implicit none
  private

  ! Public methods
  public :: hd_activate

contains

  subroutine hd_activate(par_files)
    use mod_hd_phys, only: hd_phys_init
    use mod_hd_hllc, only: hd_hllc_init
    use mod_hd_roe, only: hd_roe_init

    character(len=*), intent(in) :: par_files(:)

    call hd_phys_init(par_files)
    call hd_hllc_init()
    call hd_roe_init()
  end subroutine hd_activate

end module mod_hd
