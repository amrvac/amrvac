module mod_hd

  implicit none
  public

contains

  subroutine hd_activate(par_files)
    character(len=*), intent(in) :: par_files(:)

    call hd_phys_init(par_files)
    call hd_hllc_init()
    call hd_roe_init()
  end subroutine hd_activate

end module mod_hd
