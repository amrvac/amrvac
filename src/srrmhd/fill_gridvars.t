{#IFDEF PARTICLES
!=============================================================================
subroutine fill_gridvars(mygridvars,mypw)

! Not jet adjusted for SRRMHD case!!!

use mod_gridvars
use constants
include 'amrvacdef.f'

type(walloc), dimension(ngridshi), intent(in)   :: mypw
type(walloc), dimension(ngridshi), intent(out)  :: mygridvars
integer                                   :: igrid, iigrid, idir
double precision, dimension(ixG^T,1:ndir) :: beta
double precision, dimension(ixG^T,1:nw)   :: w
{#IFDEF PARTICLES_GCA
double precision, dimension(ixG^T,1:ndir) :: ue, bhat
double precision, dimension(ixG^T)        :: kappa, kappa_B, absB, tmp
}
!-----------------------------------------------------------------------------

do iigrid=1,igridstail; igrid=igrids(iigrid);

   mygridvars(igrid)%w(ixG^T,1:ngridvars) = 0.0d0
   ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
   w(ixG^T,1:nw) = mypw(igrid)%w(ixG^T,1:nw)
   call primitive(ixG^LL,ixG^LL,w,px(igrid)%x)

{#IFDEF PARTICLES_ADVECT
   ! fill with velocity:
   mygridvars(igrid)%w(ixG^T,vp1_:vp^NC_) = w(ixG^T,v1_:v^NC_) * UNIT_VELOCITY
   if (useprimitiveRel) then
      do idir=vp1_,vp^NC_
         mygridvars(igrid)%w(ixG^T,idir) = mygridvars(igrid)%w(ixG^T,idir)/w(ixG^T,lfac_) 
      end do
   end if
}

{#IFDEF PARTICLES_LORENTZ

! fill with magnetic field:
   mygridvars(igrid)%w(ixG^T,bp1_:bp^NC_) = w(ixG^T,b1_:b^NC_) &
        * sqrt(4.0d0*dpi*UNIT_VELOCITY**2.0d0 * UNIT_DENSITY)


! get the four-velocity (we get the electric field later from e = b x beta):
   mygridvars(igrid)%w(ixG^T,up1_:up^NC_) = w(ixG^T,u1_:u^NC_)
   if (.not.useprimitiveRel) then
      do idir=1,ndir
         mygridvars(igrid)%w(ixG^T,up1_+idir-1) = mygridvars(igrid)%w(ixG^T,up1_+idir-1) * w(ixG^T,lfac_) 
      end do
   end if
   

!   select case (typeaxial)
!   case ('slab')
!      mygridvars(igrid)%w(ixG^T,ep1_) = mygridvars(igrid)%w(ixG^T,2) * beta(ixG^T,3) &
!           - mygridvars(igrid)%w(ixG^T,3) * beta(ixG^T,2)
!      
!      mygridvars(igrid)%w(ixG^T,ep2_) = mygridvars(igrid)%w(ixG^T,3) * beta(ixG^T,1) &
!           - mygridvars(igrid)%w(ixG^T,1) * beta(ixG^T,3)
!      
!      mygridvars(igrid)%w(ixG^T,ep3_) = mygridvars(igrid)%w(ixG^T,1) * beta(ixG^T,2) &
!           - mygridvars(igrid)%w(ixG^T,2) * beta(ixG^T,1)
!{^IFPHI
!   case ('cylindrical')
!      mygridvars(igrid)%w(ixG^T,ep1_) = mygridvars(igrid)%w(ixG^T,^PHI) * beta(ixG^T,^Z) &
!           - mygridvars(igrid)%w(ixG^T,^Z) * beta(ixG^T,^PHI)
!      
!      mygridvars(igrid)%w(ixG^T,ep^PHI_) = mygridvars(igrid)%w(ixG^T,^Z) * beta(ixG^T,1) &
!           - mygridvars(igrid)%w(ixG^T,1) * beta(ixG^T,^Z)
!      
!      mygridvars(igrid)%w(ixG^T,ep^Z_) = mygridvars(igrid)%w(ixG^T,1) * beta(ixG^T,^PHI) &
!           - mygridvars(igrid)%w(ixG^T,^PHI) * beta(ixG^T,1)
!}
!   case default
!      call mpistop('geometry not implemented in fill_gridvars, or phi not properly set')
!   end select

}

{#IFDEF PARTICLES_GCA

   ! fill with magnetic field:
   mygridvars(igrid)%w(ixG^T,bp1_:bp^NC_) = w(ixG^T,b1_:b^NC_) &
        * sqrt(4.0d0*dpi*UNIT_VELOCITY**2.0d0 * UNIT_DENSITY)

   ! get the electric field
   beta(ixG^T,1:ndir) = w(ixG^T,u1_:u^NC_)
   if (useprimitiveRel) then
      do idir=1,ndir
         beta(ixG^T,idir) = beta(ixG^T,idir)/w(ixG^T,lfac_) 
      end do
   end if

   select case (typeaxial)
   case ('slab')
      mygridvars(igrid)%w(ixG^T,ep1_) = mygridvars(igrid)%w(ixG^T,2) * beta(ixG^T,3) &
           - mygridvars(igrid)%w(ixG^T,3) * beta(ixG^T,2)
      
      mygridvars(igrid)%w(ixG^T,ep2_) = mygridvars(igrid)%w(ixG^T,3) * beta(ixG^T,1) &
           - mygridvars(igrid)%w(ixG^T,1) * beta(ixG^T,3)
      
      mygridvars(igrid)%w(ixG^T,ep3_) = mygridvars(igrid)%w(ixG^T,1) * beta(ixG^T,2) &
           - mygridvars(igrid)%w(ixG^T,2) * beta(ixG^T,1)
   case default
      call mpistop('sorry, PARTICLES_GCA currently only works for cartesian coordinates.')
   end select

   ! grad(kappa B)
   absB(ixG^T) = sqrt({^C& mygridvars(igrid)%w(ixG^T,bp^C_)**2 |+})
   ue(ixG^T,1) = mygridvars(igrid)%w(ixG^T,ep2_) * mygridvars(igrid)%w(ixG^T,bp3_) &
        - mygridvars(igrid)%w(ixG^T,ep3_) * mygridvars(igrid)%w(ixG^T,bp2_)
   ue(ixG^T,2) = mygridvars(igrid)%w(ixG^T,ep3_) * mygridvars(igrid)%w(ixG^T,bp1_) &
        - mygridvars(igrid)%w(ixG^T,ep1_) * mygridvars(igrid)%w(ixG^T,bp3_)
   ue(ixG^T,3) = mygridvars(igrid)%w(ixG^T,ep1_) * mygridvars(igrid)%w(ixG^T,bp2_) &
        - mygridvars(igrid)%w(ixG^T,ep2_) * mygridvars(igrid)%w(ixG^T,bp1_)
   do idir=1,ndir
      ue(ixG^T,idir) = ue(ixG^T,idir) * CONST_c / absB(ixG^T)**2
   end do

   kappa(ixG^T) = sqrt(1.0d0 - ({^C& ue(ixG^T,^C)**2 |+})/CONST_c**2)
   kappa_B(ixG^T) = kappa(ixG^T) * absB(ixG^T)

   do idir=1,ndim
      call gradient(kappa_B,ixG^LL^LSUB1,idir,tmp)
      mygridvars(igrid)%w(ixG^T,grad_kappa_B1_-1+idir) = tmp(ixG^T)/UNIT_LENGTH
   end do

   {^C& bhat(ixG^T,^C) = mygridvars(igrid)%w(ixG^T,bp^C_) / absB(ixG^T)\}

   ! (b dot grad) b and the other directional derivatives
   
   {^C&
   do idir=1,ndim

      call gradient(bhat(ixG^T,^C),ixG^LL^LSUB1,idir,tmp)
      mygridvars(igrid)%w(ixG^T,b_dot_grad_b^C_) = mygridvars(igrid)%w(ixG^T,b_dot_grad_b^C_) &
           + bhat(ixG^T,idir) * tmp(ixG^T)/UNIT_LENGTH
      mygridvars(igrid)%w(ixG^T,ue_dot_grad_b^C_) = mygridvars(igrid)%w(ixG^T,ue_dot_grad_b^C_) &
           + ue(ixG^T,idir) * tmp(ixG^T)/UNIT_LENGTH

      call gradient(ue(ixG^T,^C),ixG^LL^LSUB1,idir,tmp)
      mygridvars(igrid)%w(ixG^T,b_dot_grad_ue^C_) = mygridvars(igrid)%w(ixG^T,b_dot_grad_ue^C_) &
           + bhat(ixG^T,idir) * tmp(ixG^T)/UNIT_LENGTH
      mygridvars(igrid)%w(ixG^T,ue_dot_grad_ue^C_) = mygridvars(igrid)%w(ixG^T,ue_dot_grad_b^C_) &
           + ue(ixG^T,idir) * tmp(ixG^T)/UNIT_LENGTH

   end do
   }

}

end do

end subroutine fill_gridvars
!=============================================================================
}
