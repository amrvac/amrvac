# 1 "mod_usr.f"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "mod_usr.f"
module mod_usr

  use mod_ard
  implicit none

contains

  subroutine usr_init()
     call set_coordinate_system('Cartesian')
     call ard_activate()
        
     usr_init_one_grid => ebr_init
  end subroutine usr_init

  subroutine ebr_init(ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixmin1,ixmin2,ixmax1,&
     ixmax2,w,x)
     integer, intent(in)             :: ixGmin1,ixGmin2,ixGmax1,ixGmax2,&
         ixmin1,ixmin2,ixmax1,ixmax2
     double precision, intent(in)    :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
        1:ndim)
     double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
        1:nw)
     double precision                :: urand(ixmin1:ixmax1,ixmin2:ixmax2)
     double precision                :: ssu, ssv, ssw

     ssu = br_A
     ssv = br_B/br_A
     ssw = (br_A * br_C) / (br_D) 

     call random_number(urand)
     w(ixmin1:ixmax1,ixmin2:ixmax2, u_) = ssu + 1.0d-1*urand
     call random_number(urand)
     w(ixmin1:ixmax1,ixmin2:ixmax2, v_) = ssv + 1.0d-1*urand
     call random_number(urand)
     w(ixmin1:ixmax1,ixmin2:ixmax2, w_) = ssw + 1.0d-1*urand
  end subroutine ebr_init

end module mod_usr
