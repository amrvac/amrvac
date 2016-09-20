!=============================================================================
! amrvacusr.t.testmhd

! INCLUDE:amrvacnul/specialini.t
INCLUDE:amrvacnul/speciallog.t
INCLUDE:amrvacnul/specialbound.t
INCLUDE:amrvacnul/specialimpl.t
!INCLUDE:amrvacnul/specialsource.t
INCLUDE:amrvacnul/usrflags.t
!=============================================================================
subroutine initglobaldata_usr

use mod_global_parameters
!-----------------------------------------------------------------------------
eqpar(gamma_)=5.0d0/3.0d0
eqpar(eta_)=zero

select case(iprob)
 case(16)
   eqpar(eta_)=0.00001d0
 case(18)
   eqpar(eta_)=0.0001d0
 case(20)
   eqpar(gamma_)=1.4d0
 case(21)
   eqpar(gamma_)=2.0d0
endselect 

end subroutine initglobaldata_usr
!=============================================================================
subroutine initonegrid_usr(ixG^L,ix^L,w,x)

! initialize one grid 

use mod_global_parameters

integer, intent(in) :: ixG^L,ix^L
double precision, intent(in) :: x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)

logical::          first
integer::          n1
double precision:: rholeft,rhoright,eleft,eright,m^Cleft,m^Cright,slocx^D
double precision:: v^Cleft,pleft,v^Cright,pright,gamm
double precision:: bx,byleft,bzleft,byright,bzright
double precision:: qv,width,dv,k1,sigma,rho0,p0,dA,dM,alfa,Lx,yjet
double precision:: rhodisk,rr0,rr1,x1disk,x2disk,v0
double precision:: cs,machs,macha,ff,Rjet
double precision, dimension(ixG^T) :: r,fslope,rinv 

integer,dimension(:),allocatable:: seed
integer::  ix2,seed_size
real::             ranx1d(ixGlo1:ixGhi1)
double precision:: ranx(ixG^T)

logical, dimension(ixG^T)           :: patchw(ixG^T)

data first/.true./
!----------------------------------------------------------------------------
oktest = index(teststr,'initonegrid_usr')>=1
if (oktest) write(unitterm,*) ' === initonegrid_usr  (in ) : ', &
                'ixG^L : ',ixG^L

if (typephys/='mhd') then
   call mpistop("test problems are all MHD problems: set typephys!")
end if

select case(iprob)
 case(22)
   ! Torrilhon test for 1.75D MHD setamrvac -d=13
   gamm=eqpar(gamma_)

   bx=one

   rholeft=one
   pleft=one
   m^Cleft=zero;
   byleft=one
   bzleft=zero

   rhoright=0.2d0
   pright=0.2d0
   m^Cright=zero;
   byright=cos(3.0d0)
   bzright=sin(3.0d0)

   if(first.and.mype==0) then
     print *,'Torrilhon test'
     print *,'by=',byright,' bz=',bzright
     first=.false.
   endif

   eleft=pleft/(gamm-one)+half*(bx**2+byleft**2+bzleft**2)
   eright=pright/(gamm-one)+half*(bx**2+byright**2+bzright**2)

   slocx1=zero
   where({^D&x(ixG^S,^D)<=slocx^D|.or.})
      w(ixG^S,rho_)     = rholeft
      {^C&w(ixG^S,m^C_) = m^Cleft \}
      w(ixG^S,e_ )      = eleft
      w(ixG^S,b1_ )     = bx
      w(ixG^S,b2_ )     = byleft
      {^IFTHREEC w(ixG^S,b3_ )     = bzleft}
   elsewhere
      w(ixG^S,rho_)     = rhoright
      {^C&w(ixG^S,m^C_) = m^Cright \}
      w(ixG^S,e_  )     = eright
      w(ixG^S,b1_ )     = bx
      w(ixG^S,b2_ )     = byright
      {^IFTHREEC w(ixG^S,b3_ )     = bzright}
   endwhere
 case(21)
   !setamrvac -d=13 -z=0 -phi=0 -g=16 -p=mhd -u=testmhd
   ! Brio-Wu
   {^NOONED call mpistop("iprob=-1 is 1.5D MHD brio-wu") }

   bx       = 0.75d0

   rholeft  = one
   {^C&v^Cleft=zero\}
   byleft   = one
   bzleft   = zero
   pleft    = (eqpar(gamma_)-one)*(1.78125d0-half*(bx**2+byleft**2+bzleft**2))

   rhoright = 0.125d0
   {v^Cright=zero \}
   byright  = -1.0d0
   bzright  = zero
   pright    = (eqpar(gamma_)-one)*(0.88125d0-half*(bx**2+byright**2+bzright**2))

   slocx1=half*(xprobmax1-xprobmin1)
   where({^D&x(ixG^S,^D)<=slocx^D|.or.})
      w(ixG^S,rho_)     = rholeft
      {^C&w(ixG^S,v^C_) = v^Cleft \}
      w(ixG^S,p_ )      = pleft
      w(ixG^S,b1_ )     = bx
      w(ixG^S,b2_ )     = byleft
      {^IFTHREEC w(ixG^S,b3_ )     = bzleft}
   elsewhere
      w(ixG^S,rho_)     = rhoright
      {^C&w(ixG^S,v^C_) = v^Cright \}
      w(ixG^S,p_  )     = pright
      w(ixG^S,b1_ )     = bx
      w(ixG^S,b2_ )     = byright
      {^IFTHREEC w(ixG^S,b3_ )     = bzright}
   endwhere

   patchw(ix^S)=.false.
   call conserve(ixG^L,ix^L,w,x,patchw)
 case(0)
   !setamrvac -d=13 -z=0 -phi=0 -g=16 -p=mhd -u=testmhd
   ! stationary contact
   {^NOONED call mpistop("iprob=0 is 1.5D (M)HD stationary contact") }

   bx       = zero

   rholeft  = 10.0d0
   {^C&v^Cleft=zero\}
   byleft   = zero
   bzleft   = zero
   pleft    = 20.0d0

   rhoright = 1.0d0
   {v^Cright=zero \}
   byright  = zero
   bzright  = zero
   pright   = 20.0d0

   slocx1=half*(xprobmax1-xprobmin1)
   where({^D&x(ixG^S,^D)<=slocx^D|.or.})
      w(ixG^S,rho_)     = rholeft
      {^C&w(ixG^S,v^C_) = v^Cleft \}
      w(ixG^S,p_ )      = pleft
      w(ixG^S,b1_ )     = bx
      w(ixG^S,b2_ )     = byleft
      {^IFTHREEC w(ixG^S,b3_ )     = bzleft}
   elsewhere
      w(ixG^S,rho_)     = rhoright
      {^C&w(ixG^S,v^C_) = v^Cright \}
      w(ixG^S,p_  )     = pright
      w(ixG^S,b1_ )     = bx
      w(ixG^S,b2_ )     = byright
      {^IFTHREEC w(ixG^S,b3_ )     = bzright}
   endwhere

   patchw(ix^S)=.false.
   call conserve(ixG^L,ix^L,w,x,patchw)
 case(1)
   !setamrvac -d=13 -z=0 -phi=0 -g=16 -p=mhd -u=testmhd
   ! S. Li/ JCP, 203, 2005, 344-357, test 1
   {^NOONED call mpistop("iprob=1 is 1.5D MHD") }

   bx       = 2.0d0/dsqrt(4.0d0*dpi)

   rholeft  = 1.08d0
   v1left   = 1.2d0
   v2left   = 0.01d0
   {^IFTHREEC v3left   = 0.5d0}
   pleft    = 0.95d0
   byleft   = 3.6d0/dsqrt(4.0d0*dpi)
   bzleft   = 2.0d0/dsqrt(2.0d0*dpi)

   rhoright = 1.0d0
   v1right  = zero
   v2right  = zero
   {^IFTHREEC v3right  = zero}
   pright   = one
   byright  = 4.0d0/dsqrt(4.0d0*dpi)
   bzright  = 2.0d0/dsqrt(2.0d0*dpi)

   slocx1=half*(xprobmax1-xprobmin1)
   where({^D&x(ixG^S,^D)<=slocx^D|.or.})
      w(ixG^S,rho_)     = rholeft
      {^C&w(ixG^S,v^C_) = v^Cleft \}
      w(ixG^S,p_ )      = pleft
      w(ixG^S,b1_ )     = bx
      w(ixG^S,b2_ )     = byleft
      {^IFTHREEC w(ixG^S,b3_ )     = bzleft}
   elsewhere
      w(ixG^S,rho_)     = rhoright
      {^C&w(ixG^S,v^C_) = v^Cright \}
      w(ixG^S,p_  )     = pright
      w(ixG^S,b1_ )     = bx
      w(ixG^S,b2_ )     = byright
      {^IFTHREEC w(ixG^S,b3_ )     = bzright}
   endwhere

   patchw(ix^S)=.false.
   call conserve(ixG^L,ix^L,w,x,patchw)
 case(2)
   !setamrvac -d=13 -z=0 -phi=0 -g=16 -p=mhd -u=testmhd
   ! Ryu, Dongsu; Jones, T. W., 1995, ApJ, 442, 228, test 1, fig 1a
   {^NOONED call mpistop("iprob=2 is 1.5D MHD test") }

   bx       = 5.0d0/dsqrt(4.0d0*dpi)

   rholeft  = 1.0d0
   v1left   = 10.0d0
   v2left   = zero
   {^IFTHREEC  v3left   = zero}
   byleft   = 5.0d0/dsqrt(4.0d0*dpi)
   bzleft   = zero
   pleft    = 20.0d0

   rhoright = 1.0d0
   v1right  = -10.0d0
   v2right  = zero
   {^IFTHREEC v3right  = zero}
   byright  = 5.0d0/dsqrt(4.0d0*dpi)
   bzright  = zero
   pright   = 1.0d0

   slocx1=half*(xprobmax1-xprobmin1)
   where({^D&x(ixG^S,^D)<=slocx^D|.or.})
      w(ixG^S,rho_)     = rholeft
      {^C&w(ixG^S,v^C_) = v^Cleft \}
      w(ixG^S,p_ )      = pleft
      w(ixG^S,b1_ )     = bx
      w(ixG^S,b2_ )     = byleft
      {^IFTHREEC w(ixG^S,b3_ )     = bzleft}
   elsewhere
      w(ixG^S,rho_)     = rhoright
      {^C&w(ixG^S,v^C_) = v^Cright \}
      w(ixG^S,p_  )     = pright
      w(ixG^S,b1_ )     = bx
      w(ixG^S,b2_ )     = byright
      {^IFTHREEC w(ixG^S,b3_ )     = bzright}
   endwhere


   patchw(ix^S)=.false.
   call conserve(ixG^L,ix^L,w,x,patchw)
 case(3)
   ! variant of Ryu, Dongsu; Jones, T. W., 1995, ApJ, 442, 228, test 1
   {^NOONED call mpistop("iprob=3 is 1.5D MHD test") }

   bx       = 5.0d0/dsqrt(4.0d0*dpi)

   rholeft  = 1.0d0
   v1left   = 10.0d0
   v2left   = zero
   {^IFTHREEC v3left   = zero}
   byleft   = 5.0d0/dsqrt(4.0d0*dpi)
   bzleft   = zero
   pleft    = 20.0d0

   rhoright = 1.0d0
   v1right  = -10.0d0
   v2right  = zero
   {^IFTHREEC v3right  = zero}
   byright  = 5.0d0/dsqrt(4.0d0*dpi)
   bzright  = zero
   pright   = 20.0d0

   slocx1=half*(xprobmax1-xprobmin1)
   where({^D&x(ixG^S,^D)<=slocx^D|.or.})
      w(ixG^S,rho_)     = rholeft
      {^C&w(ixG^S,v^C_) = v^Cleft \}
      w(ixG^S,p_ )      = pleft
      w(ixG^S,b1_ )     = bx
      w(ixG^S,b2_ )     = byleft
      {^IFTHREEC w(ixG^S,b3_ )     = bzleft}
   elsewhere
      w(ixG^S,rho_)     = rhoright
      {^C&w(ixG^S,v^C_) = v^Cright \}
      w(ixG^S,p_  )     = pright
      w(ixG^S,b1_ )     = bx
      w(ixG^S,b2_ )     = byright
      {^IFTHREEC w(ixG^S,b3_ )     = bzright}
   endwhere


   patchw(ix^S)=.false.
   call conserve(ixG^L,ix^L,w,x,patchw)
 case(4)
   !  Ryu, Dongsu; Jones, T. W., 1995, ApJ, 442, 228, test 2, fig 1b
   {^NOONED call mpistop("iprob=4 is 1.5D MHD test") }

   bx       = 3.0d0/dsqrt(4.0d0*dpi)

   rholeft  = 1.0d0
   v1left   = zero
   v2left   = zero
   {^IFTHREEC v3left   = zero}
   byleft   = 5.0d0/dsqrt(4.0d0*dpi)
   bzleft   = zero
   pleft    = 1.0d0

   rhoright = 0.1d0
   v1right  = zero
   v2right  = zero
   {^IFTHREEC v3right  = zero}
   byright  = 2.0d0/dsqrt(4.0d0*dpi)
   bzright  = zero
   pright   = 10.0d0

   slocx1=half*(xprobmax1-xprobmin1)
   where({^D&x(ixG^S,^D)<=slocx^D|.or.})
      w(ixG^S,rho_)     = rholeft
      {^C&w(ixG^S,v^C_) = v^Cleft \}
      w(ixG^S,p_ )      = pleft
      w(ixG^S,b1_ )     = bx
      w(ixG^S,b2_ )     = byleft
      {^IFTHREEC w(ixG^S,b3_ )     = bzleft}
   elsewhere
      w(ixG^S,rho_)     = rhoright
      {^C&w(ixG^S,v^C_) = v^Cright \}
      w(ixG^S,p_  )     = pright
      w(ixG^S,b1_ )     = bx
      w(ixG^S,b2_ )     = byright
      {^IFTHREEC w(ixG^S,b3_ )     = bzright}
   endwhere


   patchw(ix^S)=.false.
   call conserve(ixG^L,ix^L,w,x,patchw)
 case(5)
   !  Ryu, Dongsu; Jones, T. W., 1995, ApJ, 442, 228, test 3  (fig 2a)
   !  setamrvac -d=13 -z=0 -phi=0 -g=16 -p=mhd -u=testmhd
   {^NOONED call mpistop("iprob=5 is 1.5D MHD test") }

   bx       = 2.0d0/dsqrt(4.0d0*dpi)

   rholeft  = 1.08d0
   v1left   = 1.2d0
   v2left   = 0.01d0
   {^IFTHREEC v3left   = 0.5d0}
   byleft   = 3.6d0/dsqrt(4.0d0*dpi)
   bzleft   = 2.0d0/dsqrt(4.0d0*dpi)
   pleft    = 0.95d0

   rhoright = 1.0d0
   v1right  = zero
   v2right  = zero
   {^IFTHREEC    v3right  = zero}
   byright  = 4.0d0/dsqrt(4.0d0*dpi)
   bzright  = 2.0d0/dsqrt(4.0d0*dpi)
   pright   = 1.0d0

   slocx1=half*(xprobmax1-xprobmin1)
   where({^D&x(ixG^S,^D)<=slocx^D|.or.})
      w(ixG^S,rho_)     = rholeft
      {^C&w(ixG^S,v^C_) = v^Cleft \}
      w(ixG^S,p_ )      = pleft
      w(ixG^S,b1_ )     = bx
      w(ixG^S,b2_ )     = byleft
      {^IFTHREEC w(ixG^S,b3_ )     = bzleft}
   elsewhere
      w(ixG^S,rho_)     = rhoright
      {^C&w(ixG^S,v^C_) = v^Cright \}
      w(ixG^S,p_  )     = pright
      w(ixG^S,b1_ )     = bx
      w(ixG^S,b2_ )     = byright
      {^IFTHREEC w(ixG^S,b3_ )     = bzright}
   endwhere

   patchw(ix^S)=.false.
   call conserve(ixG^L,ix^L,w,x,patchw)
 case(6)
   !  Ryu, Dongsu; Jones, T. W., 1995, ApJ, 442, 228, test 4 (fig 2b)
   !  setamrvac -d=13 -z=0 -phi=0 -g=16 -p=mhd -u=testmhd
   {^NOONED call mpistop("iprob=6 is 1.5D MHD test") }

   bx       = 3.0d0/dsqrt(4.0d0*dpi)

   rholeft  = 1.0d0
   v1left   = 0.0d0
   v2left   = 0.0d0
   {^IFTHREEC v3left   = 0.0d0}
   byleft   = 6.0d0/dsqrt(4.0d0*dpi)
   bzleft   = 0.0d0
   pleft    = 1.0d0

   rhoright = 0.1d0
   v1right  = 0.0d0
   v2right  = 2.0d0
   {^IFTHREEC v3right  = 1.0d0}
   byright  = 1.0d0/dsqrt(4.0d0*dpi)
   bzright  = 0.0d0
   pright   = 10.0d0

   slocx1=half*(xprobmax1-xprobmin1)
   where({^D&x(ixG^S,^D)<=slocx^D|.or.})
      w(ixG^S,rho_)     = rholeft
      {^C&w(ixG^S,v^C_) = v^Cleft \}
      w(ixG^S,p_ )      = pleft
      w(ixG^S,b1_ )     = bx
      w(ixG^S,b2_ )     = byleft
      {^IFTHREEC w(ixG^S,b3_ )     = bzleft}
   elsewhere
      w(ixG^S,rho_)     = rhoright
      {^C&w(ixG^S,v^C_) = v^Cright \}
      w(ixG^S,p_  )     = pright
      w(ixG^S,b1_ )     = bx
      w(ixG^S,b2_ )     = byright
      {^IFTHREEC w(ixG^S,b3_ )     = bzright}
   endwhere

   patchw(ix^S)=.false.
   call conserve(ixG^L,ix^L,w,x,patchw)
 case(7)
   !  Ryu, Dongsu; Jones, T. W., 1995, ApJ, 442, 228, test 5 (fig 3a)
   !  setamrvac -d=13 -z=0 -phi=0 -g=16 -p=mhd -u=testmhd
   {^NOONED call mpistop("iprob=7 is 1.5D MHD test") }

   bx       = zero

   rholeft  = 0.1d0
   v1left   = 50.0d0
   v2left   = 0.0d0
   {^IFTHREEC v3left   = 0.0d0}
   byleft   = -1.0d0/dsqrt(4.0d0*dpi)
   bzleft   = -2.0d0/dsqrt(4.0d0*dpi)
   pleft    = 0.4d0

   rhoright = 0.1d0
   v1right  = 0.0d0
   v2right  = 0.0d0
   {^IFTHREEC v3right  = 0.0d0}
   byright  = 1.0d0/dsqrt(4.0d0*dpi)
   bzright  = 2.0d0/dsqrt(4.0d0*dpi)
   pright   = 0.2d0

   slocx1=half*(xprobmax1-xprobmin1)
   where({^D&x(ixG^S,^D)<=slocx^D|.or.})
      w(ixG^S,rho_)     = rholeft
      {^C&w(ixG^S,v^C_) = v^Cleft \}
      w(ixG^S,p_ )      = pleft
      w(ixG^S,b1_ )     = bx
      w(ixG^S,b2_ )     = byleft
      {^IFTHREEC w(ixG^S,b3_ )     = bzleft}
   elsewhere
      w(ixG^S,rho_)     = rhoright
      {^C&w(ixG^S,v^C_) = v^Cright \}
      w(ixG^S,p_  )     = pright
      w(ixG^S,b1_ )     = bx
      w(ixG^S,b2_ )     = byright
      {^IFTHREEC w(ixG^S,b3_ )     = bzright}
   endwhere

   patchw(ix^S)=.false.
   call conserve(ixG^L,ix^L,w,x,patchw)
 case(8)
   !  Ryu, Dongsu; Jones, T. W., 1995, ApJ, 442, 228, test 6 (fig 3b)
   !  setamrvac -d=13 -z=0 -phi=0 -g=16 -p=mhd -u=testmhd
   {^NOONED call mpistop("iprob=8 is 1.5D MHD test") }

   bx       = zero

   rholeft  = 1.0d0
   v1left   = -1.0d0
   v2left   = 0.0d0
   {^IFTHREEC v3left   = 0.0d0}
   byleft   = 1.0d0
   bzleft   = 0.0d0
   pleft    = 1.0d0

   rhoright = 1.0d0
   v1right  = 1.0d0
   v2right  = 0.0d0
   {^IFTHREEC v3right  = 0.0d0}
   byright  = 1.0d0
   bzright  = 0.0d0
   pright   = 1.0d0

   slocx1=half*(xprobmax1-xprobmin1)
   where({^D&x(ixG^S,^D)<=slocx^D|.or.})
      w(ixG^S,rho_)     = rholeft
      {^C&w(ixG^S,v^C_) = v^Cleft \}
      w(ixG^S,p_ )      = pleft
      w(ixG^S,b1_ )     = bx
      w(ixG^S,b2_ )     = byleft
      {^IFTHREEC w(ixG^S,b3_ )     = bzleft}
   elsewhere
      w(ixG^S,rho_)     = rhoright
      {^C&w(ixG^S,v^C_) = v^Cright \}
      w(ixG^S,p_  )     = pright
      w(ixG^S,b1_ )     = bx
      w(ixG^S,b2_ )     = byright
      {^IFTHREEC w(ixG^S,b3_ )     = bzright}
   endwhere

   patchw(ix^S)=.false.
   call conserve(ixG^L,ix^L,w,x,patchw)
 case(9)
   !  Ryu, Dongsu; Jones, T. W., 1995, ApJ, 442, 228, test 7 (fig 4a)
   {^NOONED call mpistop("iprob=9 is 1.5D MHD test") }

   bx       = 1.0d0

   rholeft  = 1.0d0
   v1left   = zero
   v2left   = zero
   {^IFTHREEC v3left   = zero}
   byleft   = one
   bzleft   = zero
   pleft    = 1.0d0

   rhoright = 0.2d0
   v1right  = zero
   v2right  = zero
   {^IFTHREEC v3right  = zero}
   byright  = zero
   bzright  = zero
   pright   = 0.1d0

   slocx1=half*(xprobmax1-xprobmin1)
   where({^D&x(ixG^S,^D)<=slocx^D|.or.})
      w(ixG^S,rho_)     = rholeft
      {^C&w(ixG^S,v^C_) = v^Cleft \}
      w(ixG^S,p_ )      = pleft
      w(ixG^S,b1_ )     = bx
      w(ixG^S,b2_ )     = byleft
      {^IFTHREEC w(ixG^S,b3_ )     = bzleft}
   elsewhere
      w(ixG^S,rho_)     = rhoright
      {^C&w(ixG^S,v^C_) = v^Cright \}
      w(ixG^S,p_  )     = pright
      w(ixG^S,b1_ )     = bx
      w(ixG^S,b2_ )     = byright
      {^IFTHREEC w(ixG^S,b3_ )     = bzright}
   endwhere

   patchw(ix^S)=.false.
   call conserve(ixG^L,ix^L,w,x,patchw)
 case(10)
   !  Ryu, Dongsu; Jones, T. W., 1995, ApJ, 442, 228, test 8 (fig 4b)
   {^NOONED call mpistop("iprob=10 is 1.5D MHD test") }

   bx       = 1.3d0

   rholeft  = 0.4d0
   v1left   =-0.66991d0
   v2left   = 0.98263d0
   {^IFTHREEC v3left   = zero}
   byleft   = 0.0025293d0
   bzleft   = zero
   pleft    = 0.52467d0

   rhoright = 1.0d0
   v1right  = zero
   v2right  = zero
   {^IFTHREEC v3right  = zero}
   byright  = one
   bzright  = zero
   pright   = 1.0d0

   slocx1=half*(xprobmax1-xprobmin1)
   where({^D&x(ixG^S,^D)<=slocx^D|.or.})
      w(ixG^S,rho_)     = rholeft
      {^C&w(ixG^S,v^C_) = v^Cleft \}
      w(ixG^S,p_ )      = pleft
      w(ixG^S,b1_ )     = bx
      w(ixG^S,b2_ )     = byleft
      {^IFTHREEC w(ixG^S,b3_ )     = bzleft}
   elsewhere
      w(ixG^S,rho_)     = rhoright
      {^C&w(ixG^S,v^C_) = v^Cright \}
      w(ixG^S,p_  )     = pright
      w(ixG^S,b1_ )     = bx
      w(ixG^S,b2_ )     = byright
      {^IFTHREEC w(ixG^S,b3_ )     = bzright}
   endwhere

   patchw(ix^S)=.false.
   call conserve(ixG^L,ix^L,w,x,patchw)
 case(11)
   !  Ryu, Dongsu; Jones, T. W., 1995, ApJ, 442, 228, test 9 (fig 4c)
   {^NOONED call mpistop("iprob=11 is 1.5D MHD test") }

   bx       = 0.75d0

   rholeft  = 0.65d0
   v1left   = 0.667d0
   v2left   =-0.257d0
   {^IFTHREEC  v3left   = zero}
   byleft   = 0.55d0
   bzleft   = zero
   pleft    = 0.5d0

   rhoright = 1.0d0
   v1right  = 0.4d0
   v2right  =-0.94d0
   {^IFTHREEC v3right  = zero}
   byright  = zero
   bzright  = zero
   pright   = 0.75d0

   slocx1=half*(xprobmax1-xprobmin1)
   where({^D&x(ixG^S,^D)<=slocx^D|.or.})
      w(ixG^S,rho_)     = rholeft
      {^C&w(ixG^S,v^C_) = v^Cleft \}
      w(ixG^S,p_ )      = pleft
      w(ixG^S,b1_ )     = bx
      w(ixG^S,b2_ )     = byleft
      {^IFTHREEC w(ixG^S,b3_ )     = bzleft}
   elsewhere
      w(ixG^S,rho_)     = rhoright
      {^C&w(ixG^S,v^C_) = v^Cright \}
      w(ixG^S,p_  )     = pright
      w(ixG^S,b1_ )     = bx
      w(ixG^S,b2_ )     = byright
      {^IFTHREEC w(ixG^S,b3_ )     = bzright}
   endwhere

   patchw(ix^S)=.false.
   call conserve(ixG^L,ix^L,w,x,patchw)
 case(12)
   !  Ryu, Dongsu; Jones, T. W., 1995, ApJ, 442, 228, test 10 (fig 4d)
   {^NOONED call mpistop("iprob=12 is 1.5D MHD test") }

   bx       = 0.7d0

   rholeft  = 1.0d0
   v1left   = 0.0d0
   v2left   = 0.0d0
   {^IFTHREEC v3left   = zero}
   byleft   = 0.0d0
   bzleft   = zero
   pleft    = 1.0d0

   rhoright = 0.3d0
   v1right  = 0.0d0
   v2right  = 0.0d0
   {^IFTHREEC v3right  = 1.0d0}
   byright  = 1.0d0
   bzright  = zero
   pright   = 0.2d0

   slocx1=half*(xprobmax1-xprobmin1)+xprobmin1
   where({^D&x(ixG^S,^D)<=slocx^D|.or.})
      w(ixG^S,rho_)     = rholeft
      {^C&w(ixG^S,v^C_) = v^Cleft \}
      w(ixG^S,p_ )      = pleft
      w(ixG^S,b1_ )     = bx
      w(ixG^S,b2_ )     = byleft
      {^IFTHREEC w(ixG^S,b3_ )     = bzleft}
   elsewhere
      w(ixG^S,rho_)     = rhoright
      {^C&w(ixG^S,v^C_) = v^Cright \}
      w(ixG^S,p_  )     = pright
      w(ixG^S,b1_ )     = bx
      w(ixG^S,b2_ )     = byright
      {^IFTHREEC w(ixG^S,b3_ )     = bzright}
   endwhere

   patchw(ix^S)=.false.
   call conserve(ixG^L,ix^L,w,x,patchw)
 case(13)
   !  Ryu, Dongsu; Jones, T. W., 1995, ApJ, 442, 228, test 11 (fig 5a)
   {^NOONED call mpistop("iprob=13 is 1.5D MHD test") }

   bx       = 0.75d0

   rholeft  = 1.0d0
   v1left   = 0.0d0
   v2left   = 0.0d0
   {^IFTHREEC v3left   = zero}
   byleft   = 1.0d0
   bzleft   = zero
   pleft    = 1.0d0

   rhoright = 0.125d0
   v1right  = 0.0d0
   v2right  = 0.0d0
   {^IFTHREEC v3right  = 0.0d0}
   byright  =-1.0d0
   bzright  = zero
   pright   = 0.1d0

   slocx1=half*(xprobmax1-xprobmin1)+xprobmin1
   where({^D&x(ixG^S,^D)<=slocx^D|.or.})
      w(ixG^S,rho_)     = rholeft
      {^C&w(ixG^S,v^C_) = v^Cleft \}
      w(ixG^S,p_ )      = pleft
      w(ixG^S,b1_ )     = bx
      w(ixG^S,b2_ )     = byleft
      {^IFTHREEC w(ixG^S,b3_ )     = bzleft}
   elsewhere
      w(ixG^S,rho_)     = rhoright
      {^C&w(ixG^S,v^C_) = v^Cright \}
      w(ixG^S,p_  )     = pright
      w(ixG^S,b1_ )     = bx
      w(ixG^S,b2_ )     = byright
      {^IFTHREEC w(ixG^S,b3_ )     = bzright}
   endwhere

   patchw(ix^S)=.false.
   call conserve(ixG^L,ix^L,w,x,patchw)
 case(14)
   !  Ryu, Dongsu; Jones, T. W., 1995, ApJ, 442, 228, test 12 (fig 5b)
   {^NOONED call mpistop("iprob=14 is 1.5D MHD test") }

   bx       = 1.3d0

   rholeft  = 1.0d0
   v1left   = 0.0d0
   v2left   = 0.0d0
   {^IFTHREEC v3left   = zero}
   byleft   = 1.0d0
   bzleft   = zero
   pleft    = 1.0d0

   rhoright = 0.4d0
   v1right  = 0.0d0
   v2right  = 0.0d0
   {^IFTHREEC v3right  = 0.0d0}
   byright  =-1.0d0
   bzright  = zero
   pright   = 0.4d0

   slocx1=half*(xprobmax1-xprobmin1)
   where({^D&x(ixG^S,^D)<=slocx^D|.or.})
      w(ixG^S,rho_)     = rholeft
      {^C&w(ixG^S,v^C_) = v^Cleft \}
      w(ixG^S,p_ )      = pleft
      w(ixG^S,b1_ )     = bx
      w(ixG^S,b2_ )     = byleft
      {^IFTHREEC w(ixG^S,b3_ )     = bzleft}
   elsewhere
      w(ixG^S,rho_)     = rhoright
      {^C&w(ixG^S,v^C_) = v^Cright \}
      w(ixG^S,p_  )     = pright
      w(ixG^S,b1_ )     = bx
      w(ixG^S,b2_ )     = byright
      {^IFTHREEC w(ixG^S,b3_ )     = bzright}
   endwhere

   patchw(ix^S)=.false.
   call conserve(ixG^L,ix^L,w,x,patchw)
 case(15)
   ! setamrvac -d=22 -phi=0 -z=0 -g=16,16 -p=mhd -u=testmhd
   {^IFONED   call mpistop("prob 15 is 2D") }
   {^IFTHREED call mpistop("prob 15 is 2D") }
   qv=0.645d0
   width=0.05d0
   dv=0.01d0
   k1=two*dpi
   sigma=0.20d0
   {^IFTWOD
   w(ixG^S,rho_)=one
   w(ixG^S,e_)=one
   w(ixG^S,b1_)=0.129d0
   w(ixG^S,b2_)=zero
   w(ixG^S,m1_)=qv*tanh(x(ixG^S,2)/width)
   w(ixG^S,m2_)=dv*sin(k1*x(ixG^S,1))*exp(-(x(ixG^S,2)/sigma)**2)
   w(ixG^S,e_)=w(ixG^S,e_)/(eqpar(gamma_)-1)+&
     half*(w(ixG^S,rho_)*(^C&w(ixG^S,m^C_)**2+)+(^C&w(ixG^S,b^C_)**2+))
   if(first .and. mype==0)then
      write(*,*)'Doing 2D MHD, Kelvin-Helmholtz problem'
      write(*,*)'qv, width, dv, k1, sigma:'
      write(*,*)qv,width,dv,k1,sigma
      first=.false.
   endif
   \}
 case(16)
   ! setamrvac -d=22 -phi=0 -z=0 -g=16,16 -p=mhd -u=testmhd
   {^IFONED   call mpistop("prob 16 is 2D") }
   {^IFTHREED call mpistop("prob 16 is 2D") }
   {^IFTWOD
   w(ixG^S,rho_)=one
   w(ixG^S,e_)=one
   where(x(ixG^S,2)>=zero)
      w(ixG^S,b1_)=0.129d0
   elsewhere
      w(ixG^S,b1_)=-0.129d0
   endwhere
   w(ixG^S,b2_)=zero
   qv=0.645d0
   width=0.05d0
   dv=0.01d0
   k1=two*dpi
   sigma=0.20d0
   if(first .and. mype==0)then
      write(*,*)'Doing 2D MHD, Reversed Kelvin-Helmholtz problem'
      write(*,*)'qv, width, dv, k1, sigma:'
      write(*,*)qv,width,dv,k1,sigma
      write(*,*)'Assuming eta set, using value:',eqpar(eta_)
      first=.false.
   endif
   w(ixG^S,m1_)=qv*tanh(x(ixG^S,2)/width)
   w(ixG^S,m2_)=dv*sin(k1*x(ixG^S,1))*exp(-(x(ixG^S,2)/sigma)**2)
   w(ixG^S,e_)=w(ixG^S,e_)/(eqpar(gamma_)-1)+&
     half*(w(ixG^S,rho_)*(^C&w(ixG^S,m^C_)**2+)+(^C&w(ixG^S,b^C_)**2+))
   \}
 case(17)
   ! setamrvac -d=22 -phi=0 -z=0 -g=16,16 -p=mhd -u=testmhd
   {^IFONED   call mpistop("prob 17 is 2D") }
   {^IFTHREED call mpistop("prob 17 is 2D") }
   {^IFTWOD
   rho0=25.0d0/9.0d0
   p0=5.0d0/3.0d0
   w(ixG^S,rho_)=rho0
   if(first .and. mype==0)then
      write(*,*)'Doing 2D ideal MHD, Orszag Tang problem'
      write(*,*)'rho - p - gamma:',rho0,p0,eqpar(gamma_)
      first=.false.
   endif
   w(ixG^S,m1_)=-rho0*sin(x(ixG^S,2))
   w(ixG^S,m2_)= rho0*sin(x(ixG^S,1))
   w(ixG^S,b1_)=-sin(x(ixG^S,2))
   w(ixG^S,b2_)= sin(two*x(ixG^S,1))
   w(ixG^S,e_)=p0/(eqpar(gamma_)-one)+&
     half*((^C&w(ixG^S,m^C_)**2+)/rho0+(^C&w(ixG^S,b^C_)**2+))
   \}
 case(18)
   ! setamrvac -d=23 -phi=0 -z=0 -g=16,16 -p=mhd -u=testmhd
   {^IFONED   call mpistop("prob 18 is 2.5D") }
   {^IFTHREED call mpistop("prob 18 is 2.5D") }
   {^IFTWOD
   dA=0.2d0
   dM=3.0d0
   width=one
   dv=0.01d0
   n1=1
   sigma=one
   alfa=0.35d0
   Lx=two*dpi/alfa
   rho0=one
   w(ixG^S,rho_)=rho0
   if(first .and. mype==0)then
      write(*,*)'Doing 2.5D resistive MHD, Wake problem'
      write(*,*)'A  - M - gamma:',dA,dM,eqpar(gamma_)
      write(*,*)'dv n1 sigma:',dv,n1,sigma
      write(*,*)'alfa Lx:',alfa,Lx
      write(*,*)'Assuming eta set, using value:',eqpar(eta_)
      first=.false.
   endif
   w(ixG^S,m1_)=rho0*(one-one/cosh(x(ixG^S,2)))
   w(ixG^S,m2_)=rho0*(dv*(sin(n1*two*dpi*x(ixG^S,1)/Lx)*tanh(x(ixG^S,2))&
               +cos(n1*two*dpi*x(ixG^S,1)/Lx)/cosh(x(ixG^S,2)))&
               *exp(-(x(ixG^S,2)/sigma)**2))
   w(ixG^S,b1_)=dA*tanh(x(ixG^S,2)/width)
   w(ixG^S,b2_)=zero
   if (ndir/=3) call mpistop("this is a mhd23 problem!")
   w(ixG^S,m^NC_)=zero
   w(ixG^S,b^NC_)=dA/cosh(x(ixG^S,2)/width)
   p0=one/(dM**2)/eqpar(gamma_)
   w(ixG^S,e_)=p0/(eqpar(gamma_)-1)+&
     half*((^C&w(ixG^S,m^C_)**2+)/rho0+(^C&w(ixG^S,b^C_)**2+))
   \}
 case(19)
   ! setamrvac -d=22 -phi=0 -z=0 -g=16,16 -p=mhd -u=testmhd
   ! Baty and Keppens, A&A 447, 9, 2006, case Fig 13
   {^IFONED   call mpistop("prob 19 is 2D") }
   {^IFTHREED call mpistop("prob 19 is 2D") }
   {^IFTWOD
   w(ixG^S,rho_)=one
   w(ixG^S,e_)=one
   cs=dsqrt(eqpar(gamma_))
   machs=3.0d0
   macha=7.0d0
   qv=cs*machs
   w(ixG^S,b1_)=qv/macha
   w(ixG^S,b2_)=zero
   ff=1.25d0
   Rjet=0.125d0
   dv=0.01d0
   sigma=0.2d0
   if(first)then
      write(*,*)'seeding random number generator, on mype==',mype
      call random_seed(SIZE=seed_size)
      allocate(seed(seed_size))
      call random_seed(GET=seed(1:seed_size))
   endif
   call random_number(ranx1d(ixGmin1:ixGmax1))
   do ix2=ixGmin2,ixGmax2
      ranx(ixGmin1:ixGmax1,ix2)=ranx1d(ixGmin1:ixGmax1)-0.5d0
   enddo
   if(first .and. mype==0)then
      write(*,*)'Doing 2D ideal MHD, Double Kelvin-Helmholtz problem'
      write(*,*)'cs, B0, Ma, Ms, V, Rj, F, dv, sigma:'
      write(*,*) cs,qv/macha,macha,machs,qv,Rjet,ff,dv,sigma
   endif
   if(first)then
      first=.false.
   endif
   w(ixG^S,m1_)=half*qv*(one-dtanh(ff*dabs(x(ixG^S,2))/Rjet-ff*Rjet/dabs(x(ixG^S,2))))
   w(ixG^S,m2_)=dv*ranx(ixG^S)*dexp(-((dabs(x(ixG^S,2))-Rjet)/sigma)**2)
   w(ixG^S,e_)=w(ixG^S,e_)/(eqpar(gamma_)-one)+&
     half*(w(ixG^S,rho_)*(^C&w(ixG^S,m^C_)**2+)+(^C&w(ixG^S,b^C_)**2+))
   \}
 case(20)
   ! setamrvac -d=22 -phi=0 -z=0 -g=16,16 -p=mhd -u=testmhd
   {^IFONED   stop'problem is 2D'}
   {^IFTHREED stop'problem is 2D'}
   rhodisk=10.0d0
   rr0=0.1d0
   rr1=0.115d0
   x1disk=xprobmin1+half*(xprobmax1-xprobmin1)
   {^NOONED
   x2disk=xprobmin2+half*(xprobmax2-xprobmin2)
   }
   v0=2.0d0
   if(first .and. mype==0)then
     write(*,*)'Doing 2D ideal MHD, rotor problem'
     write(*,*)'rr0,rr1:'
     write(*,*)rr0,rr1
     write(*,*)'gamma,v0:'
     write(*,*)eqpar(gamma_),v0
     first=.false.
   endif
  {^IFTWOD
  ! pressure
  w(ixG^S,e_)=one
  w(ixG^S,b1_)=5.0d0/dsqrt(4.0d0*dpi)
  w(ixG^S,b2_)=zero
  r(ixG^S)=dsqrt((x(ixG^S,1)-x1disk)**2+(x(ixG^S,2)-x2disk)**2)
  where(r(ixG^S)>rr0)
    fslope(ixG^S)=(rr1-r(ixG^S))/(rr1-rr0)
    rinv(ixG^S)=one/r(ixG^S)
  elsewhere
    fslope(ixG^S)=one
    rinv(ixG^S)=one/rr0
  endwhere
  where(r(ixG^S)<rr1)
    w(ixG^S,rho_)=one+(rhodisk-one)*fslope(ixG^S)
    w(ixG^S,m1_)=-v0*fslope(ixG^S)*(x(ixG^S,2)-x2disk)*rinv(ixG^S)
    w(ixG^S,m2_)=v0*fslope(ixG^S)*(x(ixG^S,1)-x1disk)*rinv(ixG^S)
  elsewhere
    w(ixG^S,rho_)=one
    w(ixG^S,m1_)=zero
    w(ixG^S,m2_)=zero
  endwhere

  w(ixG^S,e_)=w(ixG^S,e_)/(eqpar(gamma_)-1)+&
      half*(w(ixG^S,rho_)*(^C&w(ixG^S,m^C_)**2+)+(^C&w(ixG^S,b^C_)**2+))
  w(ixG^S,m1_)=w(ixG^S,rho_)*w(ixG^S,m1_)
  w(ixG^S,m2_)=w(ixG^S,rho_)*w(ixG^S,m2_)
  \}

case default
            write(unitterm,*)'Undefined Iprob in Userfile ',iprob
   Call mpistop(' --- initonegrid_usr ---')
end  select 

end subroutine initonegrid_usr
!=============================================================================
subroutine specialsource(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

! Calculate w(iw)=w(iw)+qdt*SOURCE[wCT,qtC,x] within ixO for all indices
! iw=iwmin...iwmax.  wCT is at time qCT

use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L, iw^LIM
double precision, intent(in) :: qdt, qtC, qt, x(ixI^S,1:ndim)
double precision, intent(inout) :: wCT(ixI^S,1:nw), w(ixI^S,1:nw)

! integer :: iw
! double precision :: s(ixG^T)
!-----------------------------------------------------------------------------

! do iw= iw^LIM
!    select case(iw)
!    case(m1_)
!       ! The source is based on the time centered wCT
!       call getmyforce(wCT,ixO^L,s)
!       w(ixO^S,m1_)=w(ixO^S,m1_) + qdt*s(ixO^S)
!    case(e_)
!       call getmyheating(wCT,ixO^L,s)
!       w(ixO^S,e_) =w(ixO^S,e_)  + qdt*s(ixO^S)
!    end select
! end do

end subroutine specialsource
!=============================================================================
subroutine getdt_special(w,ixG^L,ix^L,dtnew,dx^D,x)

! Limit "dt" further if necessary, e.g. due to the special source terms.
! The getdt_courant (CFL condition) and the getdt subroutine in the AMRVACPHYS
! module have already been called.

use mod_global_parameters

integer, intent(in) :: ixG^L, ix^L
double precision, intent(in) :: dx^D, x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw), dtnew
!-----------------------------------------------------------------------------

dtnew=bigdouble

end subroutine getdt_special
!=============================================================================
subroutine specialeta(w,ixI^L,ix^L,idirmin,x,current,eta)

! Set the "eta" array for resistive MHD based on w or the
! "current" variable which has components between idirmin and 3.

use mod_global_parameters

integer, intent(in) :: ixI^L, ix^L, idirmin
double precision, intent(in) :: w(ixI^S,nw), x(ixI^S,1:ndim)

double precision :: current(ixG^T,7-2*ndir:3), eta(ixG^T)
!-----------------------------------------------------------------------------

!  eta(ix^S)=...

call mpistop("specialeta is not defined")

end subroutine specialeta
!=============================================================================
subroutine specialrefine_grid(igrid,level,ixG^L,ix^L,qt,w,x,refine,coarsen)

! Enforce additional refinement or coarsening
! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.

use mod_global_parameters

integer, intent(in) :: igrid, level, ixG^L, ix^L
double precision, intent(in) :: qt, w(ixG^S,1:nw), x(ixG^S,1:ndim)
integer, intent(inout) :: refine, coarsen
!-----------------------------------------------------------------------------

{^IFTWOD if (iprob==16.and. (any(dabs(x(ix^S,2)) < 0.1d0))) refine=1 \}

end subroutine specialrefine_grid
!=============================================================================
subroutine specialvarforerrest(ixI^L,ixO^L,iflag,w,var)

! this is the place to compute a local auxiliary variable to be used
! as refinement criterion for the Lohner error estimator only
!  -->it is then requiring and iflag>nw
! note that ixO=ixI=ixG, hence the term local (gradients need special attention!)

use mod_global_parameters

integer, intent(in)          :: ixI^L,ixO^L,iflag
double precision, intent(in) :: w(ixI^S,1:nw)
double precision, intent(out):: var(ixG^T)
!-----------------------------------------------------------------------------

if (iflag >nw)call mpistop(' iflag> nw, make change in parfile or in user file')

var(ixI^S) = zero

end subroutine specialvarforerrest
!=============================================================================
subroutine specialset_B0(ixI^L,ixO^L,x,wB0)

! Here one can add a steady (time-independent) potential background field

use mod_global_parameters

integer, intent(in)           :: ixI^L,ixO^L
double precision, intent(in)  :: x(ixI^S,1:ndim)
double precision, intent(inout) :: wB0(ixI^S,1:ndir)
!-----------------------------------------------------------------------------
call mpistop(' abs(Busr)> 0, make change in parfile or in user file')

wB0(ixO^S,1:ndir)=wB0(ixO^S,1:ndir)
!!wB0(ixO^S,1:ndir)=wB0(ixO^S,1:ndir)+user defined steady potential field

end subroutine specialset_B0
!=============================================================================
! amrvacusr.t.testmhd
!=============================================================================
