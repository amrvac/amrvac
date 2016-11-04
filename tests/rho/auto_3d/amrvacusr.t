!=============================================================================
! amrvacusr.t.testrho

! INCLUDE:amrvacnul.specialini.t
INCLUDE:amrvacnul/speciallog.t
INCLUDE:amrvacnul/specialbound.t
INCLUDE:amrvacnul/specialsource.t
INCLUDE:amrvacnul/specialimpl.t
INCLUDE:amrvacnul/usrflags.t
!=============================================================================
subroutine initglobaldata_usr

use mod_global_parameters
!----------------------------------------------------------------------------

! {^IFONED   eqpar(v1_)=one }
! {^IFTWOD   eqpar(v1_)=one; eqpar(v2_)=one }
! {^IFTHREED eqpar(v1_)=one; eqpar(v2_)=one; eqpar(v3_)=one }

end subroutine initglobaldata_usr
!=============================================================================
subroutine initonegrid_usr(ixG^L,ix^L,w,x)

! initialize one grid 

use mod_global_parameters
use mod_rho, only: rho_

integer, intent(in) :: ixG^L, ix^L
double precision, intent(in) :: x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)

double precision:: rr2(ix^S)
double precision:: rhoflat,rhosquare,slocx^D,slocxmid^D,widthhalf

logical::          maskv(ix^S),maska(ix^S),maskc(ix^S)
double precision :: radius, xcircle^D
double precision:: xc1,yc1,xa1,xa2,ya1,ya2,xb1,xb2,yb1,yb2,xc2,yc2, &
                   rad,rad2,alp,nsig
!----------------------------------------------------------------------------

rhoflat  = 0.5d0 
rhosquare= 2.0d0 

! iprob=1 is a pure 1D Riemann problem, solvable in 1D, 2D, 3D
if (iprob==1) then
    slocx^D=0.2d0;
    where({^D&x(ix^S,^D)<=slocx^D|.and.})
       w(ix^S,rho_)     = rhosquare
    elsewhere
       w(ix^S,rho_)     = rhoflat
   endwhere
else if (iprob==2) then
    slocxmid^D=xprobmin^D+half*(xprobmax^D-xprobmin^D);
    widthhalf=0.1d0
    where({^D&((x(ix^S,^D)>slocxmid^D-widthhalf).and.&
               (x(ix^S,^D)<slocxmid^D+widthhalf))|.and.})
       w(ix^S,rho_)     = rhosquare
    elsewhere
       w(ix^S,rho_)     = rhoflat
   endwhere
else if (iprob==3) then
   {^IFONED call mpistop("iprob=3 is 2D!") }
   {^IFTHREED call mpistop("iprob=3 is 2D!") }
   {^IFTWOD
   xc1=0.25d0
   yc1=0.50d0
   rad=0.23d0
   rad2=0.13d0
   alp=dpi/3.0d0
   xa1=xc1
   ya1=yc1-rad
   xa2=xc1-rad*cos(alp)
   ya2=yc1+rad*sin(alp)
   xb1=xa1
   yb1=ya1
   xb2=xc1+rad*cos(alp)
   yb2=yc1+rad*sin(alp)
   xc2=xc1
   yc2=ya2+sqrt(rad2**2-(xa2-xc2)**2)
   maskv(ix^S)= ((x(ix^S,1)-xc1)**2+(x(ix^S,2)-yc1)**2 <= rad**2) &
            .and.(x(ix^S,2)>= (ya2-ya1)*(x(ix^S,1)-xa1)/(xa2-xa1)+ya1) & 
            .and.(x(ix^S,2)>= (yb2-yb1)*(x(ix^S,1)-xb1)/(xb2-xb1)+yb1) & 
            .and.((x(ix^S,1)-xc2)**2+(x(ix^S,2)-yc2)**2 > rad2**2) 
   xc1=0.45d0
   yc1=0.475d0
   xa1=xc1
   ya1=yc1+rad
   xa2=xc1-rad*cos(alp)
   ya2=yc1-rad*sin(alp)
   xb1=xa1
   yb1=ya1
   xb2=xc1+rad*cos(alp)
   yb2=yc1-rad*sin(alp)
   xc2=xc1
   yc2=ya2-sqrt(rad2**2-(xa2-xc2)**2)
   maska(ix^S)= ((x(ix^S,1)-xc1)**2+(x(ix^S,2)-yc1)**2 <= rad**2) &
            .and.(x(ix^S,2)<= (ya2-ya1)*(x(ix^S,1)-xa1)/(xa2-xa1)+ya1) & 
            .and.(x(ix^S,2)<= (yb2-yb1)*(x(ix^S,1)-xb1)/(xb2-xb1)+yb1) & 
            .and.((x(ix^S,1)-xc2)**2+(x(ix^S,2)-yc2)**2 > rad2**2) 
   xc1=0.75d0
   yc1=0.50d0
   alp=half*dpi-alp
   xa1=xc1-rad
   ya1=yc1
   xa2=xc1+rad*cos(alp)
   ya2=yc1+rad*sin(alp)
   xb1=xa1
   yb1=ya1
   xb2=xc1+rad*cos(alp)
   yb2=yc1-rad*sin(alp)
   yc2=yc1
   xc2=xa2+sqrt(rad2**2-(ya2-yc2)**2)
   maskc(ix^S)= ((x(ix^S,1)-xc1)**2+(x(ix^S,2)-yc1)**2 <= rad**2) &
            .and.(x(ix^S,2)<= (ya2-ya1)*(x(ix^S,1)-xa1)/(xa2-xa1)+ya1) & 
            .and.(x(ix^S,2)>= (yb2-yb1)*(x(ix^S,1)-xb1)/(xb2-xb1)+yb1) & 
            .and.((x(ix^S,1)-xc2)**2+(x(ix^S,2)-yc2)**2 > rad2**2) 
   where(maskv(ix^S).or.maska(ix^S).or.maskc(ix^S))
       w(ix^S,rho_)     = rhosquare
   elsewhere
       w(ix^S,rho_)     = rhoflat
   endwhere 
   }
else if (iprob==4) then
   {^IFONED call mpistop("iprob=4 is 2D!") }
   {^IFTHREED call mpistop("iprob=4 is 2D!") }
   {^IFTWOD
   xc1=0.5d0
   yc1=0.5d0
   rad=0.1d0
   rhosquare=0.6d0
   rhoflat=0.5d0
   maska(ix^S)=((x(ix^S,1)-xc1)**2+(x(ix^S,2)-yc1)**2 <= rad**2)
   where(maska(ix^S))
       w(ix^S,rho_)     = rhosquare+ &
          (rhoflat-rhosquare)*(x(ix^S,1)-xc1)**2/rad**2 + &
          (rhoflat-rhosquare)*(x(ix^S,2)-yc1)**2/rad**2  
   elsewhere
       w(ix^S,rho_)     = rhoflat
   endwhere 
   }
else if (iprob==5) then
   {^IFONED call mpistop("iprob=5 is 2D!") }
   {^IFTHREED call mpistop("iprob=5 is 2D!") }
   {^IFTWOD
   xc1=0.5d0
   yc1=0.5d0
   rad=0.1d0
   nsig=two
   rhosquare=2.0d0
   maska(ix^S)=((x(ix^S,1)-xc1)**2+(x(ix^S,2)-yc1)**2<=(nsig*rad)**2)
   rr2(ix^S)=(x(ix^S,1)-xc1)**2+(x(ix^S,2)-yc1)**2
   where(maska(ix^S))
       w(ix^S,rho_)     = rhosquare*exp(-rr2(ix^S)/rad**2) 
   elsewhere
       w(ix^S,rho_)     = rhosquare*exp(-nsig**2)
   endwhere 
   }
else if (iprob==6) then
   radius = 0.2d0
   xcircle^D=zero;
   where(radius**2> ^D&(x(ix^S,^D)-xcircle^D)**2+ )
      w(ix^S,rho_)     = rhosquare
   elsewhere
      w(ix^S,rho_)     = rhoflat
   endwhere
else
    call mpistop("iprob not available!")
end if

end subroutine initonegrid_usr
!=============================================================================
! amrvacusr.t.testrho
!=============================================================================
