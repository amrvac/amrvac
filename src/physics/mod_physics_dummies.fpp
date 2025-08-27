! Dummy routines which can be overwritten by a physics-dependent implementation

#:def get_cs2()
!> obtain the squared sound speed
pure real(dp) function get_cs2(u) result(cs2)
  !$acc routine seq
  real(dp), intent(in)  :: u(nw_phys)
  
  
end function get_cs2
#:enddef
