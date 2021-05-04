module mod_usr

  use mod_rd
  implicit none
    
  integer, parameter :: dp = kind(0.0d0)

contains

  subroutine usr_init()
     call set_coordinate_system('Cartesian')
     call rd_activate()
        
     usr_init_one_grid => lor_init
  end subroutine usr_init

  subroutine lor_init(ixG^L,ix^L,w,x)
     integer, intent(in)             :: ixG^L, ix^L
     double precision, intent(in)    :: x(ixG^S,1:ndim)
     double precision, intent(inout) :: w(ixG^S,1:nw)
     double precision                :: urand(ix^S)

     call random_number(urand) ! Adding just a tiny bit of noise
     w(ix^S, u_) = -10.0d0 + 1.0d-3 * urand
     w(ix^S, v_) =  30.0d0 - 1.0d-3 * urand
     w(ix^S, w_) =  20.0d0 + 1.0d-3 * urand
  end subroutine lor_init

end module mod_usr
