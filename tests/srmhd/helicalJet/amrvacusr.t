!=============================================================================
! amrvacusr.t.GCsrmhdJet
! setamrvac -d=23 -phi=3 -z=2 -g=16,16 -p=srmhd -u=GCsrmhdJet
!=============================================================================
!INCLUDE:amrvacnul/speciallog.t
INCLUDE:amrvacnul/specialimpl.t
INCLUDE:amrvacnul/correctaux_usr.t
INCLUDE:amrvacnul/usrflags.t
!=============================================================================
subroutine initglobaldata_usr

use mod_global_parameters
!-----------------------------------------------------------------------------

eqpar(Rj_)=1.5d0
eqpar(Zj_)=3.0d0

eqpar(nn_)=4.0d0
eqpar(aa_)=5.0d0

{#IFDEF GLM
eqpar(Cr_)=0.18d0
}

select case(iprob)
 case(1)
    ! Reference case, both tor/pol field, rot, beta near unity, low sigma, g=5/3
    eqpar(gamma_) =5.0d0/3.0d0
    ! azimuthal B field
    eqpar(Bazi_)=one
    ! eqpar(Bo_) fixes the axial magnetic field amplitude 
    eqpar(Bo_)=one
    ! Bc is a small constant magnetic field = magnetic field of the cloud
    eqpar(Bc_)=0.01d0
    eqpar(rhojet_)=100.0d0
    eqpar(rhocloud_)=1000.0d0
    ! velocity magnitude
    eqpar(alpha_)=9.99d0
    ! pressure
    eqpar(Pj_)=one
 case(2)
    !  both tor/pol field, rot, g=4/3
    eqpar(gamma_) =4.0d0/3.0d0
    ! azimuthal B field
    eqpar(Bazi_)=one
    ! eqpar(Bo_) fixes the axial magnetic field amplitude 
    eqpar(Bo_)=one
    ! Bc is a small constant magnetic field = magnetic field of the cloud
    eqpar(Bc_)=0.01d0
    eqpar(rhojet_)=100.0d0
    eqpar(rhocloud_)=1000.0d0
    ! velocity magnitude
    eqpar(alpha_)=9.99d0
    ! pressure
    eqpar(Pj_)=one
 case(3)
    !  weak azimuthal, almost pure pol field, no rot, g=5/3
    eqpar(gamma_) =5.0d0/3.0d0
    ! azimuthal B field
    eqpar(Bazi_)=0.01d0
    ! eqpar(Bo_) fixes the axial magnetic field amplitude 
    eqpar(Bo_)=one
    ! Bc is a small constant magnetic field = magnetic field of the cloud
    eqpar(Bc_)=0.01d0
    eqpar(rhojet_)=100.0d0
    eqpar(rhocloud_)=1000.0d0
    ! velocity magnitude
    eqpar(alpha_)=999.9d0
    ! pressure
    eqpar(Pj_)=one
 case(4)
    !  pure tor field, rot, g=5/3
    eqpar(gamma_) =5.0d0/3.0d0
    ! azimuthal B field
    eqpar(Bazi_)=one
    ! eqpar(Bo_) fixes the axial magnetic field amplitude 
    eqpar(Bo_)=zero
    ! Bc is a small constant magnetic field = magnetic field of the cloud
    eqpar(Bc_)=zero
    eqpar(rhojet_)=100.0d0
    eqpar(rhocloud_)=1000.0d0
    ! velocity magnitude
    eqpar(alpha_)=9.99d0
    ! pressure
    eqpar(Pj_)=one
 case(5)
    !  tor/pol field, rot, increased density of medium (1/50), g=5/3
    eqpar(gamma_) =5.0d0/3.0d0
    ! azimuthal B field
    eqpar(Bazi_)=one
    ! eqpar(Bo_) fixes the axial magnetic field amplitude 
    eqpar(Bo_)=one
    ! Bc is a small constant magnetic field = magnetic field of the cloud
    eqpar(Bc_)=0.01d0
    eqpar(rhojet_)=100.0d0
    eqpar(rhocloud_)=500.0d0
    ! velocity magnitude
    eqpar(alpha_)=9.99d0
    ! pressure
    eqpar(Pj_)=one
 case(6)
    !  tor/pol field, rot, higher p, g=5/3
    eqpar(gamma_) =5.0d0/3.0d0
    ! azimuthal B field
    eqpar(Bazi_)=one
    ! eqpar(Bo_) fixes the axial magnetic field amplitude 
    eqpar(Bo_)=one
    ! Bc is a small constant magnetic field = magnetic field of the cloud
    eqpar(Bc_)=0.01d0
    eqpar(rhojet_)=100.0d0
    eqpar(rhocloud_)=1000.0d0
    ! velocity magnitude
    eqpar(alpha_)=9.99d0
    ! pressure
    eqpar(Pj_)=two
 case(7)
    !  stronger tor field
    eqpar(gamma_) =5.0d0/3.0d0
    ! azimuthal B field
    eqpar(Bazi_)=two
    ! eqpar(Bo_) fixes the axial magnetic field amplitude 
    eqpar(Bo_)=one
    ! Bc is a small constant magnetic field = magnetic field of the cloud
    eqpar(Bc_)=0.01d0
    eqpar(rhojet_)=100.0d0
    eqpar(rhocloud_)=1000.0d0
    ! velocity magnitude
    eqpar(alpha_)=4.99d0
    ! pressure
    eqpar(Pj_)=one
 case(8)
    !  lower pressure
    eqpar(gamma_) =5.0d0/3.0d0
    ! azimuthal B field
    eqpar(Bazi_)=one
    ! eqpar(Bo_) fixes the axial magnetic field amplitude 
    eqpar(Bo_)=one
    ! Bc is a small constant magnetic field = magnetic field of the cloud
    eqpar(Bc_)=0.01d0
    eqpar(rhojet_)=100.0d0
    eqpar(rhocloud_)=1000.0d0
    ! velocity magnitude
    eqpar(alpha_)=9.99d0
    ! pressure
    eqpar(Pj_)=0.5d0
 case(9)
    !  lower velocity
    eqpar(gamma_) =5.0d0/3.0d0
    ! azimuthal B field
    eqpar(Bazi_)=one
    ! eqpar(Bo_) fixes the axial magnetic field amplitude 
    eqpar(Bo_)=one
    ! Bc is a small constant magnetic field = magnetic field of the cloud
    eqpar(Bc_)=0.01d0
    eqpar(rhojet_)=100.0d0
    eqpar(rhocloud_)=1000.0d0
    ! velocity magnitude
    eqpar(alpha_)=9.9d0
    ! pressure
    eqpar(Pj_)=one
 case(10)
    !  stronger pol field
    eqpar(gamma_) =5.0d0/3.0d0
    ! azimuthal B field
    eqpar(Bazi_)=one
    ! eqpar(Bo_) fixes the axial magnetic field amplitude 
    eqpar(Bo_)=two
    ! Bc is a small constant magnetic field = magnetic field of the cloud
    eqpar(Bc_)=0.01d0
    eqpar(rhojet_)=100.0d0
    eqpar(rhocloud_)=1000.0d0
    ! velocity magnitude
    eqpar(alpha_)=9.99d0
    ! pressure
    eqpar(Pj_)=one
 case(11)
    !  lower v
    eqpar(gamma_) =5.0d0/3.0d0
    ! azimuthal B field
    eqpar(Bazi_)=one
    ! eqpar(Bo_) fixes the axial magnetic field amplitude 
    eqpar(Bo_)=one
    ! Bc is a small constant magnetic field = magnetic field of the cloud
    eqpar(Bc_)=0.01d0
    eqpar(rhojet_)=100.0d0
    eqpar(rhocloud_)=1000.0d0
    ! velocity magnitude
    eqpar(alpha_)=4.99d0
    ! pressure
    eqpar(Pj_)=one
end select

end subroutine initglobaldata_usr
!=============================================================================
subroutine initonegrid_usr(ixG^L,ix^L,w,x)

! initialize one grid

use mod_global_parameters

integer, intent(in) :: ixG^L,ix^L
double precision, intent(in) :: x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)

double precision:: n
DOUBLE PRECISION:: R(ixG^T),Z(ixG^T),hlpR(ixG^T),hlpphi(ixG^T)

logical:: patchw(ixG^T)
!----------------------------------------------------------------------------

n= eqpar(nn_)  ! Parameter for the magnetic field 

! eqpar(Rj_) radial extension of the jet
R(ixG^S)=x(ixG^S,1)/eqpar(Rj_) 
! eqpar(Zj_) vertical size of the jet at time=0
Z(ixG^S)=x(ixG^S,2)/eqpar(Zj_) 

where(R(ixG^S)<one .and. Z(ixG^S)<one)
     w(ixG^S,rho_) = eqpar(rhojet_)
     w(ixG^S,bphi_)= eqpar(Bazi_)*TANH(R(ixG^S)*eqpar(Rj_)/eqpar(aa_))
elsewhere
     w(ixG^S,rho_) = eqpar(rhocloud_)
     w(ixG^S,bphi_)= zero
end where
   
w(ixG^S,vr_)  = zero

hlpphi(ixG^S)=cosh(Z(ixG^S)**n)
hlpR(ixG^S)=(cosh(R(ixG^S)**two))**2.0d0
w(ixG^S,bz_)= eqpar(Bo_)/(hlpphi(ixG^S)*hlpR(ixG^S)) +  eqpar(Bc_)

hlpR(ixG^S)=(n*(Z(ixG^S)**(n-1.))*TANH(Z(ixG^S)**n)*TANH(R(ixG^S)**2.))/&
            (two*hlpphi(ixG^S)*R(ixG^S))
w(ixG^S,br_)= hlpR(ixG^S)*eqpar(Bo_)*eqpar(Rj_)/eqpar(Zj_)

w(ixG^S,vz_)=eqpar(alpha_)*w(ixG^S,bphi_) &
             /(dsqrt(w(ixG^S,rho_))*R(ixG^S)*eqpar(Rj_)/eqpar(aa_))
w(ixG^S,vphi_)=w(ixG^S,bphi_)/dsqrt(w(ixG^S,rho_))
w(ixG^S,pp_)= eqpar(Pj_)+(half*(eqpar(Bo_))**2.0d0-&
    ( w(ixG^S,bphi_)**2.0d0 + w(ixG^S,bz_)**2.0d0)*half)


w(ixG^S,lfac_) = one/DSQRT(one - ({^C&w(ixG^S,v^C_)**2+}))

if(useprimitiveRel)then
    {^C&w(ixG^S,u^C_)=w(ixG^S,lfac_)*w(ixG^S,v^C_)\}
endif

patchw(ixG^S)=.false.
call conserve(ixG^L,ixG^L,w,x,patchw)


return
end subroutine initonegrid_usr
!=============================================================================
subroutine specialbound_usr(qt,ixG^L,ixO^L,iw,iB,w,x)

! special boundary types, user defined

use mod_global_parameters

integer, intent(in) :: ixO^L, iw, iB, ixG^L
double precision, intent(in) :: qt, x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw)

double precision::R(ixG^T),Z(ixG^T),n
double precision::hlpR(ixG^T),hlpphi(ixG^T)
logical ::patchw(ixG^T)
integer ::ix1,ix2
!----------------------------------------------------------------------------

n= eqpar(nn_)

R(ixO^S)=x(ixO^S,1)/eqpar(Rj_)
Z(ixO^S)=x(ixO^S,2)/eqpar(Zj_)

select case(iB)
case(3)
select case(iw)
case(rho_)

hlpphi(ixO^S)=cosh(Z(ixO^S)**n)
hlpR(ixO^S)=(cosh(R(ixO^S)**two))**2.0d0
where(ABS(R(ixO^S)) <= one)
     w(ixO^S,rho_) = eqpar(rhojet_)
     w(ixO^S,bphi_)= eqpar(Bazi_)*TANH(R(ixO^S)*eqpar(Rj_)/eqpar(aa_))
     w(ixO^S,bz_)  = eqpar(Bo_)/(hlpphi(ixO^S)*hlpR(ixO^S)) +  eqpar(Bc_)
     w(ixO^S,vr_)  = zero
endwhere

hlpR(ixO^S)=(n*(Z(ixO^S)**(n-1.))*TANH(Z(ixO^S)**n)*TANH(R(ixO^S)**2.))/&
            (two*hlpphi(ixO^S)*R(ixO^S))
where(ABS(R(ixO^S)) <= one)
    w(ixO^S,br_)= hlpR(ixO^S)*eqpar(Bo_)*eqpar(Rj_)/eqpar(Zj_)
    w(ixO^S,vz_)=eqpar(alpha_)*w(ixO^S,bphi_) &
             /(dsqrt(w(ixO^S,rho_))*R(ixO^S)*eqpar(Rj_)/eqpar(aa_))
    w(ixO^S,vphi_)=w(ixO^S,bphi_)/dsqrt(w(ixO^S,rho_))
    w(ixO^S,pp_)= eqpar(Pj_)+(half*(eqpar(Bo_))**2.0d0-&
       ( w(ixO^S,bphi_)**2.0d0 + w(ixO^S,bz_)**2.0d0)*half)
endwhere


where(ABS(R(ixO^S)) <= one)
   w(ixO^S,lfac_) = one/DSQRT(one - ({^C&w(ixO^S,v^C_)**2+}))
endwhere

 
do ix2 = ixOmin2,ixOmax2
where(abs(R(ixOmin1:ixOmax1,ix2)) > one)
   w(ixOmin1:ixOmax1,ix2,d_)    = w(ixOmin1:ixOmax1,2*ixOmax2-ix2+1,d_)
   w(ixOmin1:ixOmax1,ix2,sr_)   = w(ixOmin1:ixOmax1,2*ixOmax2-ix2+1,sr_)
   w(ixOmin1:ixOmax1,ix2,sz_)   =-w(ixOmin1:ixOmax1,2*ixOmax2-ix2+1,sz_)
   w(ixOmin1:ixOmax1,ix2,sphi_) = w(ixOmin1:ixOmax1,2*ixOmax2-ix2+1,sphi_)
   w(ixOmin1:ixOmax1,ix2,e_)    = w(ixOmin1:ixOmax1,2*ixOmax2-ix2+1,e_)
   w(ixOmin1:ixOmax1,ix2,br_)   =-w(ixOmin1:ixOmax1,2*ixOmax2-ix2+1,br_)
   w(ixOmin1:ixOmax1,ix2,bz_)   = w(ixOmin1:ixOmax1,2*ixOmax2-ix2+1,bz_)
   w(ixOmin1:ixOmax1,ix2,bphi_) =-w(ixOmin1:ixOmax1,2*ixOmax2-ix2+1,bphi_)
   w(ixOmin1:ixOmax1,ix2,lfac_) = w(ixOmin1:ixOmax1,2*ixOmax2-ix2+1,lfac_)
   w(ixOmin1:ixOmax1,ix2,xi_)   = w(ixOmin1:ixOmax1,2*ixOmax2-ix2+1,xi_)
endwhere
enddo

{#IFDEF GLM
w(ixO^S,psi_) = zero
}

if(useprimitiveRel)then
  where(abs(R(ixO^S)) <= one)
     {^C&w(ixO^S,u^C_)=w(ixO^S,lfac_)*w(ixO^S,v^C_)\}
  endwhere
endif

where(abs(R(ixO^S))<= one)
   patchw(ixO^S) = .false.
elsewhere
   patchw(ixO^S) = .true.
endwhere
call conserve(ixG^L,ixO^L,w,x,patchw)

case(vr_)
    ! dummy
case(vz_)
    ! dummy
case(vphi_)
    ! dummy
case(e_)
    ! dummy
case(br_)
    ! dummy
case(bz_)
    ! dummy
case(bphi_)
    ! dummy
case(lfac_)
    ! dummy
case(xi_)
    ! dummy
end select
end select

return
end subroutine specialbound_usr

!=============================================================================
subroutine bc_int(level,qt,ixG^L,ixO^L,w,x)

! internal boundary, user defined
!
! This subroutine can be used to artificially overwrite ALL conservative 
! variables in a user-selected region of the mesh, and thereby act as
! an internal boundary region. It is called just before external (ghost cell)
! boundary regions will be set by the BC selection. Here, you could e.g. 
! want to introduce an extra variable (nwextra, to be distinguished from nwaux)
! which can be used to identify the internal boundary region location.
! Its effect should always be local as it acts on the mesh.
!

use mod_global_parameters

integer, intent(in) :: ixG^L,ixO^L,level
double precision, intent(in) :: qt
double precision, intent(inout) :: w(ixG^S,1:nw)
double precision, intent(in) :: x(ixG^S,1:ndim)

! .. local ..
!logical :: patchw(ixG^T)
!----------------------------------------------------------------------------

call mpistop("bc_int not defined")

! just to give an example for relativistic MHD
!  -----------------------------------------
!patchw(ixO^S)=.true.
!where (({^D&x(ixO^S,^D)**2+})<half**2.0d0) 
!    patchw(ixO^S) = .false.
!  ^C&w(ixO^S,v^C_)=zero;
!  ^C&w(ixO^S,b^C_)=zero;
!    w(ixO^S,b3_) = one
!    w(ixO^S,v1_) = 0.99
!    w(ixO^S,rho_) = 1.d0
!    w(ixO^S,pp_)  = 2.0d0
!    w(ixO^S,lfac_)=one/dsqrt(one-({^C&w(ixO^S,v^C_)**2.0d0+}))
!end where
!!if (useprimitiveRel) then
!!  where (({^D&x(ixO^S,^D)**2+})<half**2.0d0) 
!!  {^C&w(ixO^S,u^C_)=w(ixO^S,lfac_)*w(ixO^S,v^C_);\}
!!  end where
!!endif
!call conserve(ixG^L,ixO^L,w,x,patchw)

end subroutine bc_int
!=============================================================================
subroutine specialsource(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

use mod_global_parameters

integer, intent(in) :: ixI^L, ixO^L, iw^LIM
double precision, intent(in) :: qdt, qtC, qt, x(ixI^S,1:ndim)
double precision, intent(inout) :: w(ixI^S,1:nw), wCT(ixI^S,1:nw)

! coordinate info can be used in region ixO

! integer :: iw
! double precision :: s(ixG^T)
!-----------------------------------------------------------------------------

end subroutine specialsource
!=============================================================================
subroutine getdt_special(w,ixG^L,ix^L,dtnew,dx^D,x)

use mod_global_parameters

integer, intent(in) :: ixG^L, ix^L
double precision, intent(in) :: dx^D, x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw), dtnew
!-----------------------------------------------------------------------------
dtnew=bigdouble

end subroutine getdt_special
!=============================================================================
subroutine specialeta(w,ixI^L,ix^L,idirmin,x,current,eta)

! Set the common "eta" array for resistive MHD based on w or the
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

integer, intent(in) :: igrid, level, ix^L, ixG^L
double precision, intent(in) :: qt, w(ixG^S,nw), x(ixG^S,ndim)
integer, intent(inout) :: refine, coarsen

double precision::R(ixG^T),Z(ixG^T)
!-----------------------------------------------------------------------------
R(ix^S)=x(ix^S,1)/eqpar(Rj_)
Z(ix^S)=x(ix^S,2)/eqpar(Zj_)

if (any((R(ix^S) <= 3.0d0) .and. (Z(ix^S)<= 3.0d0))) refine=1

end subroutine specialrefine_grid
!=============================================================================
subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)

! this subroutine can be used in convert, to add auxiliary variables to the
! converted output file, for further analysis using tecplot, paraview, ....
! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
!
! the array normconv can be filled in the (nw+1:nw+nwauxio) range with
! corresponding normalization values (default value 1)

use mod_global_parameters

integer, intent(in)                :: ixI^L,ixO^L
double precision, intent(in)       :: x(ixI^S,1:ndim)
double precision                   :: w(ixI^S,nw+nwauxio)
double precision                   :: normconv(0:nw+nwauxio)
! .. local ..
double precision:: divb(ixG^T)
double precision :: current(ixG^T,7-2*ndir:3)
integer          :: idirmin
!-----------------------------------------------------------------------------

!call mpistop("special output file undefined")
call getdivb(w,ixI^L,ixO^L,divb)
w(ixO^S,nw+1)=divb(ixO^S)

end subroutine specialvar_output
!=============================================================================
subroutine specialvarnames_output

! newly added variables need to be concatenated with the varnames/primnames string

use mod_global_parameters
!-----------------------------------------------------------------------------

!call mpistop("special varnames and primnames undefined")
primnames= TRIM(primnames)//' '//'divb'
wnames=TRIM(wnames)//' '//'divb'

end subroutine specialvarnames_output
!=============================================================================



















!=============================================================================
subroutine printlog_special

use mod_global_parameters
!-----------------------------------------------------------------------------
call mpistop("special log file undefined")

end subroutine printlog_special
!=============================================================================
subroutine process_grid_usr(igrid,level,ixI^L,ixO^L,qt,w,x)

! this subroutine is ONLY to be used for computing auxiliary variables
! which happen to be non-local (like div v), and are in no way used for
! flux computations. As auxiliaries, they are also not advanced

use mod_global_parameters

integer, intent(in):: igrid,level,ixI^L,ixO^L
double precision, intent(in):: qt,x(ixI^S,1:ndim)
double precision, intent(inout):: w(ixI^S,1:nw)
!-----------------------------------------------------------------------------

end subroutine process_grid_usr
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
double precision, intent(in)  :: x(ixG^T,1:ndim)
double precision, intent(inout) :: wB0(ixI^S,1:ndir)
!-----------------------------------------------------------------------------
call mpistop(' abs(Busr)> 0, make change in parfile or in user file')

wB0(ixO^S,1:ndir)=wB0(ixO^S,1:ndir)
!!wB0(ixO^S,1:ndir)=wB0(ixO^S,1:ndir)+user defined steady potential field

end subroutine specialset_B0
!=============================================================================
! amrvacusr.t.GCsrmhdJet
!=============================================================================
