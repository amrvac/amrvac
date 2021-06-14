module mod_usr

  use mod_rd
  implicit none

contains

  subroutine usr_init()
     call set_coordinate_system('Cartesian')
     call rd_activate()
        
     usr_init_one_grid => ebr_init
  end subroutine usr_init

  subroutine ebr_init(ixG^L,ix^L,w,x)
     integer, intent(in)             :: ixG^L, ix^L
     double precision, intent(in)    :: x(ixG^S,1:ndim)
     double precision, intent(inout) :: w(ixG^S,1:nw)
     double precision                :: urand(ix^S)
     double precision                :: ssu, ssv, ssw

     ssu = br_A
     ssv = br_B/br_A
     ssw = (br_A * br_C) / (br_D) 

     call random_number(urand)
     w(ix^S, u_) = ssu + 1.0d-1*urand
     call random_number(urand)
     w(ix^S, v_) = ssv + 1.0d-1*urand
     call random_number(urand)
     w(ix^S, w_) = ssw + 1.0d-1*urand
  end subroutine ebr_init

end module mod_usr
