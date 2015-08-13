!#############################################################################
! Module amrvacphys/- srrmhd
! Special relativistic resistive MHD as in Komissarov 2007
! 2015-30-07 by Oliver Porth

INCLUDE:amrvacnul/roe.t
INCLUDE:amrvacnul/hllc.t
!=============================================================================
subroutine checkglobaldata
include 'amrvacdef.f'
!-----------------------------------------------------------------------------
minrho = max(zero,smallrho)
{#IFDEF ISO
minp=eqpar(adiab_)*minrho**eqpar(gamma_)
govergminone =eqpar(gamma_)/(eqpar(gamma_)-one)
}
{#IFDEF GAMMA
govergminone =eqpar(gamma_)/(eqpar(gamma_)-one)
minp  = max(zero,smallp)
smallxi=minrho+minp*govergminone
smalltau = minp/(eqpar(gamma_)-one)
}
{#IFDEF SYNGE
call smallvaluesEOS
}

! Check if ssplitresis = .true.
if (ssplitresis .neqv. .true.) &
     call mpistop('Please run with ssplitresis = .true.')

! We require three vector components
if (^NC .ne. 3) &
     call mpistop('SRRMHD needs three vector components.')

end subroutine checkglobaldata
!=============================================================================
subroutine initglobaldata

! place to set entropy fixes etc, absent for now

include 'amrvacdef.f'
!-----------------------------------------------------------------------------
eqpar(adiab_)=1.0d0

{#IFNDEF SYNGE
eqpar(gamma_)=4.0d0/3.0d0
}
{#IFDEF SYNGE
eqpar(gamma_)=5.d0/3.d0
}

eqpar(kappa_)= 0.5d0

if(strictzero)limitvalue=zero
if(.not.strictzero)limitvalue=smalldouble**2.0d0

end subroutine initglobaldata
!=============================================================================
subroutine checkw(checkprimitive,ixI^L,ixO^L,w,flag)

include 'amrvacdef.f'
  
logical, intent(in)          :: checkprimitive
integer, intent(in)          :: ixI^L,ixO^L
double precision, intent(in) :: w(ixI^S,1:nw)
logical, intent(out)         :: flag(ixG^T)
!-----------------------------------------------------------------------------

flag(ixG^T)=.true.

{#IFDEF ENERGY
if (checkprimitive) then
  if(useprimitiveRel)then
     ! check   rho>=0, p>=smallp
     flag(ixO^S) = (w(ixO^S,rho_) > minrho).and. &
                   (w(ixO^S,pp_)  >=minp &
{#IFDEF EPSINF         .and. w(ixO^S,rho1_) > minrho&
                    .and. w(ixO^S,n_) > minrho})
  else
    ! check  v^2 < 1, rho>=0, p>=minp
     ! v^2 < 1, rho>0, p>0
     flag(ixO^S) = (({^C&w(ixO^S,v^C_)**2.0d0+})< one).and. &
                   (w(ixO^S,rho_) > minrho).and. &
                   (w(ixO^S,pp_) >=minp&
{#IFDEF EPSINF         .and. w(ixO^S,rho1_) > minrho&
                    .and. w(ixO^S,n_) > minrho})
  endif
else
  ! checks on the conservative variables
  flag(ixO^S)= (w(ixO^S,d_)   > minrho).and. &
               (w(ixO^S,tau_) > smalltau&
{#IFDEF EPSINF         .and. w(ixO^S,Drho1_) > minrho&
                    .and. w(ixO^S,Dn_) > minrho})
end if
}
{#IFNDEF ENERGY
if (checkprimitive) then
  if(useprimitiveRel)then
     ! check   rho>=0
     flag(ixO^S) = (w(ixO^S,rho_) > minrho &
{#IFDEF EPSINF         ).and. w(ixO^S,rho1_) > minrho&
                    .and. w(ixO^S,n_) > minrho})
  else
    ! check  v^2 < 1, rho>=0
     ! v^2 < 1, rho>0, p>0
     flag(ixO^S) = (({^C&w(ixO^S,v^C_)**2.0d0+})< one).and. &
                   (w(ixO^S,rho_) > minrho &
{#IFDEF EPSINF          .and. w(ixO^S,rho1_) > minrho&
                    .and. w(ixO^S,n_) > minrho})
  endif
else
  ! checks on the conservative variables
  flag(ixO^S)= (w(ixO^S,d_)   > minrho &
{#IFDEF EPSINF         .and. w(ixO^S,Drho1_) > minrho&
                    .and. w(ixO^S,Dn_) > minrho})
end if
}

end subroutine checkw
!=============================================================================
subroutine conserve(ixI^L,ixO^L,w,x,patchw)

! Transform primitive variables into conservative ones
! (rho,v,p,B,E) ---> (D,S,tau,B,E,lfac,xi)
! call to smallvalues
! --> latter only used for correcting procedure in correctaux
! --> input array patchw for spatially selective transformation

include 'amrvacdef.f'

integer, intent(in)               :: ixI^L,ixO^L
double precision, intent(inout)   :: w(ixI^S,nw)
logical, intent(in)               :: patchw(ixG^T)
double precision, intent(in)      :: x(ixI^S,1:ndim)

!-----------------------------------------------------------------------------

call conserven(ixI^L,ixO^L,w,patchw)

if(fixsmall) call smallvalues(w,x,ixI^L,ixO^L,patchw(ixO^S),'conserve')

end subroutine conserve
!=============================================================================
subroutine conserven(ixI^L,ixO^L,w,patchw)

! Transform primitive variables into conservative ones
! (rho,v,p,B) ---> (D,S,tau,B,lfac,xi)
! no call to smallvalues
! --> latter only used for correcting procedure in correctaux
! --> input array patchw for spatially selective transformation

include 'amrvacdef.f'

integer, intent(in)               :: ixI^L,ixO^L
double precision, intent(inout)   :: w(ixI^S,nw)
logical, intent(in)               :: patchw(ixG^T)
! .. local ..
double precision, dimension(ixG^T):: ecrossbi, sqrU, sqrE, sqrB, rhoh
!-----------------------------------------------------------------------------

! fill the auxilary variable lfac (lorentz factor)
where(.not.patchw(ixO^S))
   sqrU(ixO^S)    = {^C&w(ixO^S,u^C_)**2|+}
   w(ixO^S,lfac_) = dsqrt(one+sqrU(ixO^S)) 
endwhere

! fill the auxilary variable xi and density D
! with enthalpy w: xi= lfac^2 rhoh
! density: d = lfac * rho
call Enthalpy(w,ixI^L,ixO^L,patchw,rhoh)
where(.not.patchw(ixO^S))
   w(ixO^S,xi_)  = w(ixO^S,lfac_)*w(ixO^S,lfac_)* rhoh(ixO^S)
   w(ixO^S,d_)   = w(ixO^S,rho_)*w(ixO^S,lfac_)
endwhere

{#IFDEF TRACER
! We got D, now we can get the conserved tracers:
where(.not.patchw(ixO^S))
   {^FL&w(ixO^S,tr^FL_) = w(ixO^S,d_)*w(ixO^S,tr^FL_)
\}
endwhere
}{#IFDEF EPSINF
where(.not.patchw(ixO^S))
   w(ixO^S,Drho1_)   = w(ixO^S,rho1_)*w(ixO^S,lfac_)
   w(ixO^S,Drho0_)   = w(ixO^S,Drho1_)*w(ixO^S,rho0_)
   w(ixO^S,Dn_)      = w(ixO^S,n_)*w(ixO^S,lfac_)
   w(ixO^S,Dn0_)     = w(ixO^S,Dn_)*w(ixO^S,n0_)
   w(ixO^S,Depsinf_) = w(ixO^S,epsinf_)*w(ixO^S,Drho1_)**(2.0d0/3.0d0) &
        *w(ixO^S,lfac_)**(1.0d0/3.0d0)
endwhere
}

! fill the vector S
! s = E x B + xi v
{^C& 
call ecrossb(ixI^L,ixO^L,^C,w,patchw,ecrossbi)
where(.not.patchw(ixO^S))
   w(ixO^S,s^C_) = ecrossbi(ixO^S) + w(ixO^S,xi_) * w(ixO^S,u^C_)/w(ixO^S,lfac_)
end where
\}

! tau = 1/2 (E^2+B^2) + xi - p - d
! instead of E use tau= E - D
where(.not.patchw(ixO^S))
   sqrB(ixO^S)    = {^C&w(ixO^S,B^C_)**2|+}
   sqrE(ixO^S)    = {^C&w(ixO^S,e^C_)**2|+}
   w(ixO^S,tau_)  = half * (sqrE(ixO^S)+sqrB(ixO^S)) + w(ixO^S,xi_) - w(ixO^S,pp_) - w(ixO^S,d_)
endwhere

end subroutine conserven
!=============================================================================
subroutine ecrossb(ixI^L,ixO^L,idir,w,patchw,res)

include 'amrvacdef.f'

integer, intent(in)                :: ixI^L,ixO^L
double precision, intent(in)       :: w(ixI^S,1:nw)
integer         , intent(in)       :: idir
logical, intent(in)                :: patchw(ixG^T)
double precision, intent(out)      :: res(ixG^T)
! .. local ..
integer                            :: j,k
!-----------------------------------------------------------------------------

res(ixO^S) = zero
do j=1,^NC
   do k=1,^NC
      if (j .eq. k .or. j .eq. idir .or. k .eq. idir) cycle
      where(.not.patchw(ixO^S))
         res(ixO^S) = res(ixO^S) + lvc(idir,j,k)*w(ixO^S,e0_+j)*w(ixO^S,b0_+k)
      endwhere
   end do
end do

end subroutine ecrossb
!=============================================================================
subroutine vcrossb(ixI^L,ixO^L,idir,w,x,patchw,res)

! Used with conserved variables, so call to getv needed.
! assuming auxilaries are set.  
include 'amrvacdef.f'

integer, intent(in)                :: ixI^L,ixO^L
double precision, intent(in)       :: w(ixI^S,1:nw)
double precision, intent(in)       :: x(ixI^S,1:ndim)
integer         , intent(in)       :: idir
logical, intent(in)                :: patchw(ixG^T)
double precision, intent(out)      :: res(ixG^T)
! .. local ..
integer                            :: j,k
double precision, dimension(ixG^T) :: vj
!-----------------------------------------------------------------------------

res(ixO^S) = zero
do j=1,^NC
   do k=1,^NC
      if (j .eq. k .or. j .eq. idir .or. k .eq. idir) cycle
      call getv(w,x,ixI^L,ixO^L,j,vj)
      where(.not.patchw(ixO^S))
         res(ixO^S) = res(ixO^S) + lvc(idir,j,k)*vj(ixO^S)*w(ixO^S,b0_+k)
      endwhere
   end do
end do

end subroutine vcrossb
!=============================================================================
subroutine primitive(ixI^L,ixO^L,w,x)

! Transform conservative variables into primitive ones
! (D,S,tau,B)-->(rho,v,p,B,lfac,xi)

include 'amrvacdef.f'

integer, intent(in)                :: ixI^L,ixO^L
double precision, intent(inout)    :: w(ixI^S,1:nw)
double precision, intent(in)       :: x(ixI^S,1:ndim)
! .. local ..
logical, dimension(ixG^T)          :: patchw
!-----------------------------------------------------------------------------

patchw(ixO^S) = .false.
! calculate lorentz factor and xi from conservatives only
call getaux(.true.,w,x,ixI^L,ixO^L,'primitive')
call primitiven(ixI^L,ixO^L,w,patchw)


{#IFDEF ENERGY
if (tlow>zero) call fixp_usr(ixI^L,ixO^L,w,x)
}

end subroutine primitive
!==============================================================================
subroutine primitiven(ixI^L,ixO^L,w,patchw)

! same as primitive, but with padding by patchw
!  --> needed in correctaux, avoiding recursive calling by smallvalues
! assumes updated xi and lfac, and does not use smallxi cutoff
! Transform conservative variables into primitive ones
! (D,S,tau,B)-->(rho,v,p,B,lfac,xi)

include 'amrvacdef.f'

integer, intent(in)                    :: ixI^L,ixO^L
double precision, intent(inout)        :: w(ixI^S,1:nw)
logical, intent(in),dimension(ixG^T)   :: patchw
! .. local ..
double precision, dimension(ixG^T)     :: sqrB, tmpP, ecrossbi
!-----------------------------------------------------------------------------

call Pressuren(w,ixI^L,ixO^L,.true.,tmpP,patchw)
where(.not.patchw(ixO^S))
   w(ixO^S,pp_)=tmpP(ixO^S)
end where

where(.not.patchw(ixO^S))
   w(ixO^S,rho_)=w(ixO^S,d_)/w(ixO^S,lfac_)
!   sqrB(ixO^S) = {^C&w(ixO^S,b^C_)**2.0d0+}
end where


! u = lfac/xi * (S - E x B)
{^C&
call ecrossb(ixI^L,ixO^L,^C,w,patchw,ecrossbi)
where(.not.patchw(ixO^S))
   w(ixO^S,u^C_) = w(ixO^S,lfac_)/w(ixO^S,xi_) * (w(ixO^S,s^C_)-ecrossbi(ixO^S))
end where
\}

{#IFDEF TRACER
! We got lor, rho, Dtr, now we can get the tracers:
where(.not.patchw(ixO^S))
   {^FL&w(ixO^S,tr^FL_) = w(ixO^S,Dtr^FL_)/w(ixO^S,lfac_)/w(ixO^S,rho_)\}
endwhere
}
{#IFDEF EPSINF
where(.not.patchw(ixO^S))
w(ixO^S,rho1_)   = w(ixO^S,Drho1_)/w(ixO^S,lfac_)
w(ixO^S,rho0_)   = w(ixO^S,Drho0_)/w(ixO^S,lfac_)/w(ixO^S,rho1_)
w(ixO^S,n_)      = w(ixO^S,Dn_)/w(ixO^S,lfac_)
w(ixO^S,n0_)     = w(ixO^S,Dn0_)/w(ixO^S,lfac_)/w(ixO^S,n_)
w(ixO^S,epsinf_) = w(ixO^S,Depsinf_)/(w(ixO^S,lfac_) &
     *w(ixO^S,rho1_)**(2.0d0/3.0d0))
end where
}
end subroutine primitiven
!==============================================================================
subroutine e_to_rhos(ixI^L,ixO^L,w,x)

include 'amrvacdef.f'
integer:: ixI^L,ixO^L
double precision:: w(ixI^S,nw)
double precision, intent(in)      :: x(ixI^S,1:ndim)
!-------------------------------------------------------------------------

call mpistop("e to rhos for SRMHDEOS unavailable")

end subroutine e_to_rhos
!=============================================================================
subroutine rhos_to_e(ixI^L,ixO^L,w,x)

include 'amrvacdef.f'

integer:: ixI^L,ixO^L
double precision:: w(ixI^S,nw)
double precision, intent(in)      :: x(ixI^S,1:ndim)
!-----------------------------------------------------------------------------

call mpistop("e to rhos for SRMHDEOS unavailable")

end subroutine rhos_to_e
!=============================================================================
subroutine ppmflatcd(ixI^L,ixO^L,ixL^L,ixR^L,w,d2w,drho,dp)

include 'amrvacdef.f'

integer, intent(in)           :: ixI^L,ixO^L,ixL^L,ixR^L
double precision, intent(in)  :: w(ixI^S,nw),d2w(ixG^T,1:nwflux)

double precision, intent(inout) :: drho(ixG^T),dp(ixG^T)

!-----------------------------------------------------------------------------

call mpistop("PPM with flatsh=.true. can not be used with srrmhd !")


end subroutine ppmflatcd
!=============================================================================
subroutine ppmflatsh(ixI^L,ixO^L,ixLL^L,ixL^L,ixR^L,ixRR^L,idims,w,drho,dp,dv)

include 'amrvacdef.f'

integer, intent(in)           :: ixI^L,ixO^L,ixLL^L,ixL^L,ixR^L,ixRR^L
integer, intent(in)           :: idims
double precision, intent(in)  :: w(ixI^S,nw)

double precision, intent(inout) :: drho(ixG^T),dp(ixG^T),dv(ixG^T)

!-----------------------------------------------------------------------------

call mpistop("PPM with flatsh=.true. can not be used with srrmhd !")

end subroutine ppmflatsh
!=============================================================================
subroutine getv(w,x,ixI^L,ixO^L,idims,v)

! v = (s - E x B)/xi

include 'amrvacdef.f'
  
integer, intent(in)                              :: ixI^L,ixO^L,idims
double precision, intent(in)                     :: w(ixI^S,1:nw)
double precision, intent(in)                     :: x(ixI^S,1:ndim)
double precision, dimension(ixG^T), intent(out)  :: v
! .. local ..
double precision, dimension(ixG^T)               :: ecrossbi
logical, dimension(ixG^T)                        :: patchw
!-----------------------------------------------------------------------------

patchw(ixO^S) = .false.
call ecrossb(ixI^L,ixO^L,idims,w,patchw,ecrossbi)
where(w(ixO^S,xi_) <= smallxi)
   v(ixO^S) = zero
elsewhere
   v(ixO^S) = ( w(ixO^S,s0_+idims) - ecrossbi(ixO^S) ) / w(ixO^S,xi_)
endwhere

end subroutine getv
!=============================================================================
subroutine getcmax(new_cmax,w,x,ixI^L,ixO^L,idims,cmax,cmin,needcmin)
  
! Calculate cmax_idim within ixO^L
! new_cmax is not used
  
include 'amrvacdef.f'
  
logical, intent(in)           :: new_cmax,needcmin
integer, intent(in)           :: ixI^L,ixO^L,idims
double precision, intent(in)  :: w(ixI^S,1:nw)
double precision, intent(in)  :: x(ixI^S,1:ndim)
double precision, intent(out) :: cmax(ixG^T),cmin(ixG^T)
!-----------------------------------------------------------------------------

cmax(ixO^S) =   one
cmin(ixO^S) = - one

end subroutine getcmax
!=============================================================================
subroutine getflux(w,x,ixI^L,ixO^L,iw,idims,f,transport)

! Calculate non-transport flux f_idim[iw] within ixO^L.
include 'amrvacdef.f'

integer, intent(in)                :: ixI^L,ixO^L,iw,idims
double precision, intent(in)       :: w(ixI^S,nw)
double precision, intent(in)       :: x(ixI^S,1:ndim)
double precision, intent(out)      :: f(ixG^T)
logical, intent(out)               :: transport
! .. local ..
integer                            :: k
double precision, dimension(ixG^T) :: vidims, vc, p, sqrE, sqrB, ecrossbi
!-----------------------------------------------------------------------------
transport=.true.

if (iw==d_) then
   f(ixO^S)=zero
{#IFDEF TRACER
{else if (iw==tr^FL_) then 
      f(ixO^S)=zero\}
}
{#IFDEF EPSINF
else if (iw==epsinf_) then 
      f(ixO^S)=zero
else if (iw==rho0_) then 
      f(ixO^S)=zero
else if (iw==rho1_) then 
      f(ixO^S)=zero
else if (iw==n0_) then 
      f(ixO^S)=zero
else if (iw==n_) then 
      f(ixO^S)=zero
}
! f_i[phib_]=b_{i}
else if(iw==phib_) then
   transport=.false. 
{#IFNDEF FCT
       f(ixO^S)=w(ixO^S,b0_+idims)
}{#IFDEF FCT
       f(ixO^S)=zero
}
! f_i[psi_]=e_{i}
else if(iw==psi_) then
   transport=.false. 
   f(ixO^S)=w(ixO^S,e0_+idims)

! F^i[q] = q vi (transport term) + lfac*sigma*(Ei+ v x B - (E dot v) vi)
else if (iw==q_) then
   transport = .true.
   call vcrossb(ixI^L,ixO^L,idims,w,x,patchfalse,ecrossbi) ! storing vxB in tmpvar for ecrossb
   sqrE(ixO^S) = ({^C& w(ixO^S,s^C_)*w(ixO^S,e^C_)|+}) / w(ixO^S,xi_)  ! storing EdotV in tmpvar for sqrE
   call getv(w,x,ixI^L,ixO^L,idims,vidims)
   f(ixO^S) = w(ixO^S,lfac_)/eqpar(eta_) * (w(ixO^S,e0_+idims) + ecrossbi(ixO^S) - sqrE(ixO^S) * vidims(ixO^S))
else

!f^i[B_c] = lvc(c,i,k) Ek + phi delta(i,c)
{if (iw==b^C_) then
   transport=.false. 
   if (idims .ne. ^C) then
      do k=1,^NC
         if (idims .eq. k .or. ^C .eq. k) cycle
         f(ixO^S) = lvc(^C,idims,k) * w(ixO^S,e0_+k)
      end do
   else
      f(ixO^S)=w(ixO^S,phib_)
   end if
end if\}

!f^i[E_c] = - lvc(c,i,k) Bk + psi delta(i,c)
{if (iw==e^C_) then
   transport=.false. 
   if (idims .ne. ^C) then 
      do k=1,^NC
         if (idims .eq. k .or. ^C .eq. k) cycle
         f(ixO^S) = - lvc(^C,idims,k) * w(ixO^S,b0_+k)
      end do
   else
      f(ixO^S) = w(ixO^S,psi_)
   end if
end if\}


sqrE(ixO^S) = {^C& w(ixO^S,e0_+^C)**2 |+}
sqrB(ixO^S) = {^C& w(ixO^S,b0_+^C)**2 |+}

! f^i[S_c] = - E_c E_i - B_c B_i + xi v_c v_i + (1/2(E**2+B**2)+p) delta(c,i)
{if (iw==s^C_) then
   transport=.false. 
   call getv(w,x,ixI^L,ixO^L,idims,vidims)
   call getv(w,x,ixI^L,ixO^L,^C,vc)
   f(ixO^S) = - w(ixO^S,e^C_)*w(ixO^S,e0_+idims) - w(ixO^S,b^C_)*w(ixO^S,b0_+idims) &
        + w(ixO^S,xi_)*vidims(ixO^S)*vc(ixO^S)
   if (idims==^C) then !!! add (1/2(E**2+B**2)+p)
      call Pressure(w,ixI^L,ixO^L,.true.,p)
      f(ixO^S)    = f(ixO^S) + half*(sqrE(ixO^S) + sqrB(ixO^S)) + p(ixO^S)
   end if
end if\}

! f^i[tau] = tau vi (transport term) + vi*(p -1/2 E^2 -1/2 B^2) + E x B
if (iw==e_) then
   transport = .true.
   call Pressure(w,ixI^L,ixO^L,.true.,p)
   call getv(w,x,ixI^L,ixO^L,idims,vidims)
   call ecrossb(ixI^L,ixO^L,idims,w,patchfalse,ecrossbi)
   f(ixO^S) = vidims(ixO^S) * (p(ixO^S) - half*(sqrE(ixO^S)+sqrB(ixO^S)) ) + ecrossbi(ixO^S)
end if

end if

end subroutine getflux
!=============================================================================
subroutine getfluxforhllc(w,x,ixI^L,ixO^L,iw,idims,f,transport)

! Calculate non-transport flux f_idim[iw] within ixO^L.

include 'amrvacdef.f'

integer, intent(in)::          ixI^L,ixO^L,iw,idims
double precision, intent(in)  :: w(ixI^S,nw)
double precision, intent(in)      :: x(ixI^S,1:ndim)
double precision, intent(out) :: f(ixG^T,1:nwflux)
logical, intent(out)          :: transport

!-----------------------------------------------------------------------------

call getflux(w,x,ixI^L,ixO^L,iw,idims,f(ixG^T,iw),transport)

end subroutine getfluxforhllc
!=============================================================================
subroutine addgeometry(qdt,ixI^L,ixO^L,wCT,w,x)

! Add geometrical source terms to w

include 'amrvacdef.f'

integer, intent(in)                :: ixI^L, ixO^L
double precision, intent(in)       :: qdt
double precision, intent(in)       :: x(ixI^S,1:ndim)
double precision, intent(inout)    :: wCT(ixI^S,1:nw), w(ixI^S,1:nw)

integer                            :: iw,ix,idir, h1x^L{^NOONED, h2x^L}
double precision, dimension(ixG^T) :: tmp, sqrVdotB, sqrB, P
logical                            :: angmomfix=.false.
!-----------------------------------------------------------------------------

if(typeaxial /= 'slab')then
  call mpistop('Only slab geometry implemented so far')
endif

end subroutine addgeometry
!=============================================================================
subroutine addsource(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,qsourcesplit)

! w[iws]=w[iws]+qdt*S[iws,wCT] where S is the source based on wCT within ixO

include 'amrvacdef.f'

integer, intent(in)             :: ixI^L, ixO^L, iw^LIM
double precision, intent(in)    :: qdt, qtC, qt
double precision, intent(in)    :: x(ixI^S,1:ndim)
double precision, intent(in)    :: wCT(ixI^S,1:nw)
double precision, intent(inout) :: w(ixI^S,1:nw)
logical, intent(in)             :: qsourcesplit
! .. local ..
double precision                :: dx^D
!-----------------------------------------------------------------------------
dx^D=dxlevel(^D);

! Two sources: Sb is added via Strang-splitting, Sa is unsplit.
! Using ssplitresis to distinguish (needs to be set true!!!)
if(qsourcesplit .eqv. ssplitresis) then
   call addsource_b(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,dx^D)
else
   call addsource_a(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,dx^D)
end if


! now update the auxiliaries to the new state
call getaux(.true.,w,x,ixI^L,ixO^L,'addsource')

end subroutine addsource
!=============================================================================
subroutine addsource_a(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,dx^D)
include 'amrvacdef.f'

integer, intent(in) :: ixI^L, ixO^L, iw^LIM
double precision, intent(in) :: qdt, qtC, qt, wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
double precision, intent(in) :: dx^D
double precision, intent(inout) :: w(ixI^S,1:nw)
! .. local ..
double precision, dimension(ixG^T)   :: vidir
!-----------------------------------------------------------------------------

! S[psi_] = q
w(ixO^S,psi_) = w(ixO^S,psi_) + qdt * wCT(ixO^S,q_)

! S[Ei_] = - q vi
{^C&
call getv(wCT,x,ixI^L,ixO^L,^C,vidir)
w(ixO^S,e^C_) = w(ixO^S,e^C_) - qdt * wCT(ixO^S,q_) * vidir(ixO^S)
\}

end subroutine addsource_a
!=============================================================================
subroutine addsource_b(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x,dx^D)
include 'amrvacdef.f'

integer, intent(in) :: ixI^L, ixO^L, iw^LIM
double precision, intent(in) :: qdt, qtC, qt, wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
double precision, intent(in) :: dx^D
double precision, intent(inout) :: w(ixI^S,1:nw)
! .. local ..
double precision, dimension(ixG^T,1:ndir)  :: Epar, Eperp, Eperpstar, v
double precision, dimension(ixG^T)         :: v2, EdotV
!-----------------------------------------------------------------------------

! S[phib_] = -kappa phib
w(ixO^S,phib_) = w(ixO^S,phib_)*exp(-qdt*eqpar(kappa_))

! S[psi_] = - kappa psi
w(ixO^S,psi_) = w(ixO^S,psi_)*exp(-qdt*eqpar(kappa_))

! S[E] = - sigma lfac (E + v x B - (Edotv) v)
! Split the source in parallel and perpendicular components
! Then solve the linear evolution equation (assuming constant background state)
! Then add again both components.

{^C& 
call getv(wCT,x,ixI^L,ixO^L,^C,v(ixG^T,^C))
\}
v2(ixO^S) = ({^C& v(ixO^S,^C)**2|+})
EdotV(ixO^S) = {^C& w(ixO^S,e^C_)*v(ixO^S,^C) |+}

{^C&
Epar(ixO^S,^C)  = EdotV(ixO^S)*v(ixO^S,^C)/v2(ixO^S)
Eperp(ixO^S,^C) = w(ixO^S,e^C_) - Epar(ixO^S,^C)
call vcrossb(ixI^L,ixO^L,^C,wCT,x,patchfalse,Eperpstar(ixG^T,^C))
Eperpstar(ixO^S,^C) = - Eperpstar(ixO^S,^C)
\}

! Handle zero velocity separate:
where (v2(ixO^S) .lt. smalldouble**2)
{^C&
     w(ixO^S,e^C_) = w(ixO^S,e^C_) * exp(-qdt/eqpar(eta_))
\}
elsewhere
{^C&
     Epar(ixO^S,^C)  = Epar(ixO^S,^C) * exp(-qdt/(eqpar(eta_)*wCT(ixO^S,lfac_)))
     Eperp(ixO^S,^C) = Eperpstar(ixO^S,^C) + (Eperp(ixO^S,^C) - Eperpstar(ixO^S,^C))&
          * exp(-qdt*wCT(ixO^S,lfac_)/eqpar(eta_))
     w(ixO^S,e^C_)   = Epar(ixO^S,^C) + Eperp(ixO^S,^C)
\}
endwhere


end subroutine addsource_b
!=============================================================================
subroutine getdt(w,ixG^L,ix^L,dtnew,dx^D,x)

! Limit the timestp due to the resistive source addition
  
include 'amrvacdef.f'

integer, intent(in)              :: ixG^L, ix^L
double precision, intent(in)     :: dx^D, x(ixG^S,1:ndim)
double precision, intent(inout)  :: w(ixG^S,1:nw), dtnew
!-----------------------------------------------------------------------------

dtnew = dtdiffpar * eqpar(eta_)/maxval(w(ix^S,lfac_))

end subroutine getdt
!=============================================================================
subroutine getcurrent(ixI^L,ixO^L,w,x,primvar,current)

! \bm{J} = \frac{\Gamma}{\eta} [\bm{E+v\times B-(E\cdot v)v}] + q \bm{v}
! Can be used with both primitive and conserved variables.
include 'amrvacdef.f'

integer, intent(in)                       :: ixO^L, ixI^L
double precision, intent(in)              :: w(ixI^S,1:nw), x(ixI^S,1:ndim)
logical,intent(in)                        :: primvar
double precision, intent(out)             :: current(ixG^T,1:ndir)
! .. local ..
double precision, dimension(ixG^T)        :: Eperpstar, vi, EdotV
logical                                   :: primvar_sub
integer                                   :: j,k
!-----------------------------------------------------------------------------

if (primvar) then
   EdotV(ixO^S) = ({^C& w(ixO^S,u^C_)*w(ixO^S,e^C_)|+}) / w(ixO^S,lfac_)

   {^C&
   Eperpstar(ixO^S) = zero
   do j=1,ndir
      do k=1,ndir
         if (^C .eq. j .or. ^C .eq. k .or. j .eq. k) cycle
         Eperpstar(ixO^S) = Eperpstar(ixO^S) &
              + lvc(^C,j,k) * w(ixO^S,u0_+j)*w(ixO^S,b0_+k)
      end do
   end do
   Eperpstar(ixO^S) = Eperpstar(ixO^S)/w(ixO^S,lfac_)
   
   vi(ixO^S) = w(ixO^S,u^C_)/w(ixO^S,lfac_)
   
   current(ixO^S,^C) = w(ixO^S,e^C_) + Eperpstar(ixO^S) &
        - EdotV(ixO^S)*vi(ixO^S)
   current(ixO^S,^C) = w(ixO^S,lfac_)/eqpar(eta_) * current(ixO^S,^C) &
         + w(ixO^S,q_)*vi(ixO^S)
   \}

else
   EdotV(ixO^S) = ({^C& w(ixO^S,s^C_)*w(ixO^S,e^C_)|+}) / w(ixO^S,xi_)

   {^C&
   call vcrossb(ixI^L,ixO^L,^C,w,x,patchfalse,Eperpstar)
   call getv(w,x,ixI^L,ixO^L,^C,vi)

   current(ixO^S,^C) = w(ixO^S,e^C_) + Eperpstar(ixO^S) &
        - EdotV(ixO^S)*vi(ixO^S)
   current(ixO^S,^C) = w(ixO^S,lfac_)/eqpar(eta_) * current(ixO^S,^C) &
         + w(ixO^S,q_)*vi(ixO^S)
   \}
end if


end subroutine getcurrent
!=============================================================================
!=============================================================================
! end module amrvacphys/- srrmhd
!#############################################################################
