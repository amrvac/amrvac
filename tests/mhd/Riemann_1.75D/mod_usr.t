module mod_usr
  use mod_mhd
  implicit none

contains

  subroutine usr_init()
    usr_set_parameters=> initglobaldata_usr
    usr_init_one_grid => initonegrid_usr

    call set_coordinate_system('Cartesian_1.75D')
    call mhd_activate()

  end subroutine usr_init

  subroutine initglobaldata_usr
    select case(iprob)
     case(21)
       mhd_gamma=2.0d0
    endselect 

  end subroutine initglobaldata_usr

  subroutine initonegrid_usr(ixG^L,ix^L,w,x)
    ! initialize one grid 
    integer, intent(in) :: ixG^L,ix^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    double precision:: rholeft,rhoright,slocx^D
    double precision:: vleft(3),pleft,vright(3),pright
    double precision:: bx,byleft,bzleft,byright,bzright
    logical,save :: first=.true.

    select case(iprob)
     case(22)
       ! Torrilhon test for 1.75D MHD setup.pl -d=1, Cartesian_1.75D

       bx=one
       rholeft=one
       pleft=one
       vleft=zero
       byleft=one
       bzleft=zero
    
       rhoright=0.2d0
       pright=0.2d0
       vright=zero
       byright=cos(3.0d0)
       bzright=sin(3.0d0)

       if(first.and.mype==0) then
         print *,'Torrilhon test'
         print *,'by=',byright,' bz=',bzright
         first=.false.
       endif

     case(21)
       ! Brio-Wu 1.75D MHD test: setup.pl -d=1, Cartesian_1.75D 

       bx       = 0.75d0
       rholeft  = one
       vleft    = zero
       byleft   = one
       bzleft   = zero
       pleft    = (mhd_gamma-one)*(1.78125d0-half*(bx**2+byleft**2+bzleft**2))

       rhoright = 0.125d0
       vright   = zero
       byright  = -1.0d0
       bzright  = zero
       pright   = (mhd_gamma-one)*(0.88125d0-half*(bx**2+byright**2+bzright**2))

     case(0)
       ! stationary contact : setup.pl -d=1, Cartesian_1.75D

       bx       = zero
       rholeft  = 10.0d0
       vleft    = zero
       byleft   = zero
       bzleft   = zero
       pleft    = 20.0d0

       rhoright = 1.0d0
       vright   = zero
       byright  = zero
       bzright  = zero
       pright   = 20.0d0

     case(1)
       ! S. Li/ JCP, 203, 2005, 344-357, test 1
       !setup.pl -d=1, Cartesian_1.75D

       bx       = 2.0d0/dsqrt(4.0d0*dpi)
       rholeft  = 1.08d0
       vleft(1) = 1.2d0
       vleft(2) = 0.01d0
       vleft(3) = 0.5d0
       pleft    = 0.95d0
       byleft   = 3.6d0/dsqrt(4.0d0*dpi)
       bzleft   = 2.0d0/dsqrt(2.0d0*dpi)
    
       rhoright = 1.0d0
       vright   = zero
       pright   = one
       byright  = 4.0d0/dsqrt(4.0d0*dpi)
       bzright  = 2.0d0/dsqrt(2.0d0*dpi)
    
     case(2)
       ! Ryu, Dongsu; Jones, T. W., 1995, ApJ, 442, 228, test 1, fig 1a
       !setup.pl -d=1, Cartesian_1.75D

       bx       = 5.0d0/dsqrt(4.0d0*dpi)
       rholeft  = 1.0d0
       vleft(1) = 10.0d0
       vleft(2:)= zero
       byleft   = 5.0d0/dsqrt(4.0d0*dpi)
       bzleft   = zero
       pleft    = 20.0d0
    
       rhoright = 1.0d0
       vright(1)= -10.0d0
       vright(2:)  = zero
       byright  = 5.0d0/dsqrt(4.0d0*dpi)
       bzright  = zero
       pright   = 1.0d0

     case(3)
       ! variant of Ryu, Dongsu; Jones, T. W., 1995, ApJ, 442, 228, test 1

       bx       = 5.0d0/dsqrt(4.0d0*dpi)
       rholeft  = 1.0d0
       vleft(1) = 10.0d0
       vleft(2:)= zero
       byleft   = 5.0d0/dsqrt(4.0d0*dpi)
       bzleft   = zero
       pleft    = 20.0d0
    
       rhoright = 1.0d0
       vright(1)= -10.0d0
       vright(2:)= zero
       byright  = 5.0d0/dsqrt(4.0d0*dpi)
       bzright  = zero
       pright   = 20.0d0
    
     case(4)
       !  Ryu, Dongsu; Jones, T. W., 1995, ApJ, 442, 228, test 2, fig 1b

       bx       = 3.0d0/dsqrt(4.0d0*dpi)
       rholeft  = 1.0d0
       vleft    = zero
       byleft   = 5.0d0/dsqrt(4.0d0*dpi)
       bzleft   = zero
       pleft    = 1.0d0
    
       rhoright = 0.1d0
       vright   = zero
       byright  = 2.0d0/dsqrt(4.0d0*dpi)
       bzright  = zero
       pright   = 10.0d0

     case(5)
       !  Ryu, Dongsu; Jones, T. W., 1995, ApJ, 442, 228, test 3  (fig 2a)
       !  setup.pl -d=1

       bx       = 2.0d0/dsqrt(4.0d0*dpi)
       rholeft  = 1.08d0
       vleft(1) = 1.2d0
       vleft(2) = 0.01d0
       vleft(3) = 0.5d0
       byleft   = 3.6d0/dsqrt(4.0d0*dpi)
       bzleft   = 2.0d0/dsqrt(4.0d0*dpi)
       pleft    = 0.95d0
    
       rhoright = 1.0d0
       vright   = zero
       byright  = 4.0d0/dsqrt(4.0d0*dpi)
       bzright  = 2.0d0/dsqrt(4.0d0*dpi)
       pright   = 1.0d0

     case(6)
       !  Ryu, Dongsu; Jones, T. W., 1995, ApJ, 442, 228, test 4 (fig 2b)
       !  setup.pl -d=1

       bx       = 3.0d0/dsqrt(4.0d0*dpi)
       rholeft  = 1.0d0
       vleft    = 0.0d0
       byleft   = 6.0d0/dsqrt(4.0d0*dpi)
       bzleft   = 0.0d0
       pleft    = 1.0d0

       rhoright = 0.1d0
       vright(1)= 0.0d0
       vright(2)= 2.0d0
       vright(3)= 1.0d0
       byright  = 1.0d0/dsqrt(4.0d0*dpi)
       bzright  = 0.0d0
       pright   = 10.0d0

     case(7)
       !  Ryu, Dongsu; Jones, T. W., 1995, ApJ, 442, 228, test 5 (fig 3a)
       !  setup.pl -d=1

       bx       = zero
       rholeft  = 0.1d0
       vleft(1) = 50.0d0
       vleft(2:)= 0.0d0
       byleft   = -1.0d0/dsqrt(4.0d0*dpi)
       bzleft   = -2.0d0/dsqrt(4.0d0*dpi)
       pleft    = 0.4d0
    
       rhoright = 0.1d0
       vright   = 0.0d0
       byright  = 1.0d0/dsqrt(4.0d0*dpi)
       bzright  = 2.0d0/dsqrt(4.0d0*dpi)
       pright   = 0.2d0
    
     case(8)
       !  Ryu, Dongsu; Jones, T. W., 1995, ApJ, 442, 228, test 6 (fig 3b)
       !  setup.pl -d=1

       bx       = zero
       rholeft  = 1.0d0
       vleft(1) = -1.0d0
       vleft(2:)= 0.0d0
       byleft   = 1.0d0
       bzleft   = 0.0d0
       pleft    = 1.0d0
    
       rhoright = 1.0d0
       vright(1)= 1.0d0
       vright(2:)= 0.0d0
       byright  = 1.0d0
       bzright  = 0.0d0
       pright   = 1.0d0
    
     case(9)
       !  Ryu, Dongsu; Jones, T. W., 1995, ApJ, 442, 228, test 7 (fig 4a)

       bx       = 1.0d0
       rholeft  = 1.0d0
       vleft    = zero
       byleft   = one
       bzleft   = zero
       pleft    = 1.0d0
    
       rhoright = 0.2d0
       vright   = zero
       byright  = zero
       bzright  = zero
       pright   = 0.1d0
    
     case(10)
       !  Ryu, Dongsu; Jones, T. W., 1995, ApJ, 442, 228, test 8 (fig 4b)

       bx       = 1.3d0
       rholeft  = 0.4d0
       vleft(1) =-0.66991d0
       vleft(2) = 0.98263d0
       vleft(3) = zero
       byleft   = 0.0025293d0
       bzleft   = zero
       pleft    = 0.52467d0
    
       rhoright = 1.0d0
       vright   = zero
       byright  = one
       bzright  = zero
       pright   = 1.0d0
     case(11)
       !  Ryu, Dongsu; Jones, T. W., 1995, ApJ, 442, 228, test 9 (fig 4c)

       bx       = 0.75d0
       rholeft  = 0.65d0
       vleft(1) = 0.667d0
       vleft(2) =-0.257d0
       vleft(3) = zero
       byleft   = 0.55d0
       bzleft   = zero
       pleft    = 0.5d0
    
       rhoright = 1.0d0
       vright(1)= 0.4d0
       vright(2)=-0.94d0
       vright(3)= zero
       byright  = zero
       bzright  = zero
       pright   = 0.75d0
    
     case(12)
       !  Ryu, Dongsu; Jones, T. W., 1995, ApJ, 442, 228, test 10 (fig 4d)

       bx       = 0.7d0
       rholeft  = 1.0d0
       vleft    = 0.0d0
       byleft   = 0.0d0
       bzleft   = zero
       pleft    = 1.0d0
    
       rhoright = 0.3d0
       vright(1)= 0.0d0
       vright(2)= 0.0d0
       vright(3)= 1.0d0
       byright  = 1.0d0
       bzright  = zero
       pright   = 0.2d0
    
     case(13)
       !  Ryu, Dongsu; Jones, T. W., 1995, ApJ, 442, 228, test 11 (fig 5a)

       bx       = 0.75d0
       rholeft  = 1.0d0
       vleft    = 0.0d0
       byleft   = 1.0d0
       bzleft   = zero
       pleft    = 1.0d0
    
       rhoright = 0.125d0
       vright   = 0.0d0
       byright  =-1.0d0
       bzright  = zero
       pright   = 0.1d0

     case(14)
       !  Ryu, Dongsu; Jones, T. W., 1995, ApJ, 442, 228, test 12 (fig 5b)

       bx       = 1.3d0
       rholeft  = 1.0d0
       vleft    = 0.0d0
       byleft   = 1.0d0
       bzleft   = zero
       pleft    = 1.0d0
    
       rhoright = 0.4d0
       vright   = 0.0d0
       byright  =-1.0d0
       bzright  = zero
       pright   = 0.4d0

    case default
       write(unitterm,*)'Undefined Iprob in Userfile ',iprob
       Call mpistop(' --- initonegrid_usr ---')
    end  select 

    slocx1=half*(xprobmax1+xprobmin1)
    where({^D&x(ixG^S,^D)<=slocx^D|.or.})
       w(ixG^S,rho_)     = rholeft
       w(ixG^S,mom(1))   = vleft(1)
       w(ixG^S,mom(2))   = vleft(2)
       w(ixG^S,mom(3))   = vleft(3)
       w(ixG^S,p_ )      = pleft
       w(ixG^S,mag(1) )  = bx
       w(ixG^S,mag(2) )  = byleft
       w(ixG^S,mag(3) )  = bzleft
    elsewhere
       w(ixG^S,rho_)     = rhoright
       w(ixG^S,mom(1))   = vright(1)
       w(ixG^S,mom(2))   = vright(2)
       w(ixG^S,mom(3))   = vright(3)
       w(ixG^S,p_  )     = pright
       w(ixG^S,mag(1) )  = bx
       w(ixG^S,mag(2) )  = byright
       w(ixG^S,mag(3) )  = bzright
    endwhere
  
    call mhd_to_conserved(ixG^L,ix^L,w,x)
  
  end subroutine initonegrid_usr

end module mod_usr
