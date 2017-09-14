program nr
implicit none
integer :: i, N=128, iloop ! 16384 **
double precision :: rho_prev1=0., rho_prev2=0., eps=1.d-11, df, rho, f, gm, rs, g, r, dpi=3.14159265359d0, rin, rout, r1
logical :: stretched=.true.
rin =0.01d0
rout=10.d0
r1  = 1.0277248004393015d-2 ! 1.0304870605468750d-2
gm=1.0001d0 !5./3.d0 ! must be > 1
! g is Mdot / (rhoinf*(2GM/cinf2)**2*cinf) (ie dimensionless Mdot)
g=(dpi/4.d0)*(2.d0/(5.d0-3.d0*gm))**((5.d0-3.d0*gm)/(2.d0*gm-2.d0))
! rs in units of 2GM/cinf2
rs=(5.d0-3.d0*gm)/8.d0
! open(1,file='func.dat')
! r=0.1d0
! do i=1,100
!    rho=0.1*1.07**float(i)
!    write(1,*) rho, (gm/(gm-1.d0))*(rho**(gm-1.d0)-1.d0)-0.5d0          /r+(g**2.d0/(32.d0*dpi*dpi))/((r**4.d0)*(rho**2.d0))
! enddo
! close(1)
! stop
print*, 'rs', rs
print*, 'rho@rs', &
((1.d0/16.d0)*(8.d0/(5.d0-3.d0*gm))**2.d0*(2.d0/(5.d0-3.d0*gm))**((5.d0-3.d0*gm)/(2.d0*gm-2.d0)))**(2.d0/(gm+1.d0))
! - - -
r=rs
rho=1.d0
do i=-2,-1
  if (stretched) then
    r=r1*(rout/rin)**(dble(i-1)/dble(N))
  else
    r=r1+(rout-rin)*(dble(i-1)/dble(N))
  endif
! first guess
! f=(gm/(gm-1.d0))*(rho**(gm-1.d0)-1.d0)-0.5d0    /r+(g**2.d0/(32.d0*dpi*dpi))/((r**4.d0)*(rho**2.d0))
f=(1.d0/(gm-1.d0))*(rho**(gm-1.d0)-1.d0)-0.5d0/r+(g**2.d0/(32.d0*dpi*dpi))/((r**4.d0)*(rho**2.d0))
! f=(1.d0/(gm-1.d0))*(rho**(gm+1.d0)-rho )-0.5d0*rho**2.d0/r+(g**2.d0/(32.d0*dpi*dpi))/(r**4.d0)
do while(abs(f)>eps)
! df=gm*rho**(gm-2.d0)-(2.d0*((g**2.d0/(32.d0*dpi*dpi))/((r**4.d0))))/rho**3.d0
df=rho**(gm-2.d0)-(2.d0*((g**2.d0/(32.d0*dpi*dpi))/((r**4.d0))))/rho**3.d0
! df=(gm/(gm-1.d0))*rho**gm+(1.d0/(gm-1.d0))*(rho**gm-1.d0)-1.d0/r
if (rho-f/df>0.d0) then
  rho  = rho - f/df
else
  rho = rho/2.
endif
! f=(gm/(gm-1.d0))*(rho**(gm-1.d0)-1.d0)-0.5d0    /r+(g**2.d0/(32.d0*dpi*dpi))/((r**4.d0)*(rho**2.d0))
f=(1.d0/(gm-1.d0))*(rho**(gm-1.d0)-1.d0)-0.5d0/r+(g**2.d0/(32.d0*dpi*dpi))/((r**4.d0)*(rho**2.d0))
! f=(1.d0/(gm-1.d0))*(rho**(gm+1.d0)-rho )-0.5d0*rho**2.d0/r+(g**2.d0/(32.d0*dpi*dpi))/(r**4.d0)
enddo
print*, 'inner ghost cell :', r, rho, -g/(4.d0*dpi*r**2.d0)
enddo
! - - -
r=rs
open(1,file='iniState.dat')
rho=1.d0
do i=1,N
  if (stretched) then
    r=r1*(rout/rin)**(dble(i-1)/dble(N))
  else
    r=r1+(rout-rin)*(dble(i-1)/dble(N))
  endif
print*, i, r/rs
! first guess
! f=(gm/(gm-1.d0))*(rho**(gm-1.d0)-1.d0)-0.5d0    /r+(g**2.d0/(32.d0*dpi*dpi))/((r**4.d0)*(rho**2.d0))
f=(1.d0/(gm-1.d0))*(rho**(gm-1.d0)-1.d0)-0.5d0    /r+(g**2.d0/(32.d0*dpi*dpi))/((r**4.d0)*(rho**2.d0))
! f=(1.d0/(gm-1.d0))*(rho**(gm+1.d0)-rho )-0.5d0*rho**2.d0/r+(g**2.d0/(32.d0*dpi*dpi))/(r**4.d0)
iloop=0
do while(abs(f)>eps)
! df=gm*rho**(gm-2.d0)-(2.d0*((g**2.d0/(32.d0*dpi*dpi))/((r**4.d0))))/rho**3.d0
df=rho**(gm-2.d0)-(2.d0*((g**2.d0/(32.d0*dpi*dpi))/((r**4.d0))))/rho**3.d0
! df=(gm/(gm-1.d0))*rho**gm+(1.d0/(gm-1.d0))*(rho**gm-1.d0)-1.d0/r
if (rho-f/df>0.d0) then
  rho  = rho - f/df
else
  rho = rho/2.
endif
! f=(gm/(gm-1.d0))*(rho**(gm-1.d0)-1.d0)-0.5d0    /r+(g**2.d0/(32.d0*dpi*dpi))/((r**4.d0)*(rho**2.d0))
f=(1.d0/(gm-1.d0))*(rho**(gm-1.d0)-1.d0)-0.5d0/r+(g**2.d0/(32.d0*dpi*dpi))/((r**4.d0)*(rho**2.d0))
! f=(1.d0/(gm-1.d0))*(rho**(gm+1.d0)-rho )-0.5d0*rho**2.d0/r+(g**2.d0/(32.d0*dpi*dpi))/(r**4.d0)
iloop=iloop+1
if (iloop>100) then
   print*, 'extrapolation'
   rho=rho_prev2+(rho_prev2-rho_prev1)*&
   (r-r1*(rout/rin)**(dble(i-2)/dble(N)))/(r1*(rout/rin)**(dble(i-2)/dble(N))-r1*(rout/rin)**(dble(i-3)/dble(N)))
   exit
endif
enddo
! save r, rho, rho*v, mach
write(1,*), r, rho, -g/(4.d0*dpi*r**2.d0), (g/(4.d0*dpi*r**2.d0*rho))/(rho**((gm-1.d0)/2.d0))
if (i==1) then
   rho_prev2=rho
endif
rho_prev1=rho_prev2
rho_prev2=rho
enddo
close(1)
! - - -
! before leaving, compute the values to hand code in the ghost cells
do i=N+1,N+2
  if (stretched) then
    r=r1*(rout/rin)**(dble(i-1)/dble(N))
  else
    r=r1+(rout-rin)*(dble(i-1)/dble(N))
  endif
! first guess
! f=(gm/(gm-1.d0))*(rho**(gm-1.d0)-1.d0)-0.5d0    /r+(g**2.d0/(32.d0*dpi*dpi))/((r**4.d0)*(rho**2.d0))
f=(1.d0/(gm-1.d0))*(rho**(gm-1.d0)-1.d0)-0.5d0/r+(g**2.d0/(32.d0*dpi*dpi))/((r**4.d0)*(rho**2.d0))
! f=(1.d0/(gm-1.d0))*(rho**(gm+1.d0)-rho )-0.5d0*rho**2.d0/r+(g**2.d0/(32.d0*dpi*dpi))/(r**4.d0)
do while(abs(f)>eps)
! df=gm*rho**(gm-2.d0)-(2.d0*((g**2.d0/(32.d0*dpi*dpi))/((r**4.d0))))/rho**3.d0
df=rho**(gm-2.d0)-(2.d0*((g**2.d0/(32.d0*dpi*dpi))/((r**4.d0))))/rho**3.d0
! df=(gm/(gm-1.d0))*rho**gm+(1.d0/(gm-1.d0))*(rho**gm-1.d0)-1.d0/r
if (rho-f/df>0.d0) then
  rho  = rho - f/df
else
  rho = rho/2.
endif
! f=(gm/(gm-1.d0))*(rho**(gm-1.d0)-1.d0)-0.5d0    /r+(g**2.d0/(32.d0*dpi*dpi))/((r**4.d0)*(rho**2.d0))
f=(1.d0/(gm-1.d0))*(rho**(gm-1.d0)-1.d0)-0.5d0/r+(g**2.d0/(32.d0*dpi*dpi))/((r**4.d0)*(rho**2.d0))
! f=(1.d0/(gm-1.d0))*(rho**(gm+1.d0)-rho )-0.5d0*rho**2.d0/r+(g**2.d0/(32.d0*dpi*dpi))/(r**4.d0)
enddo
print*, 'outer ghost cell :', r, rho, -g/(4.d0*dpi*r**2.d0)
enddo
end
