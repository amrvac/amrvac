{#IFDEF PARTICLES
!=============================================================================
subroutine fill_gridvars(mygridvars,mypw)

use mod_gridvars
use constants
include 'amrvacdef.f'

type(walloc), dimension(ngridshi), intent(in)   :: mypw
type(walloc), dimension(ngridshi), intent(out)  :: mygridvars
integer                                   :: igrid, iigrid, idir
double precision, dimension(ixG^T,1:ndir) :: beta
double precision, dimension(ixG^T,1:nw)   :: w
double precision                          :: current(ixG^T,7-2*ndir:3)
integer                                   :: idirmin
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
}

{#IFDEF PARTICLES_SPACECRAFT
   ! fill with velocity:
   mygridvars(igrid)%w(ixG^T,vp1_:vp^NC_) = w(ixG^T,v1_:v^NC_) * UNIT_VELOCITY
}

{#IFDEF PARTICLES_LORENTZ

   ! fill with magnetic field:
   mygridvars(igrid)%w(ixG^T,bp1_:bp^NC_) = w(ixG^T,b1_:b^NC_)

! fill with electric field
   current = zero
   call getcurrent(w,ixG^LL,ixG^LLIM^D^LSUB1,idirmin,current)
   select case (typeaxial)
   case ('slab')
      mygridvars(igrid)%w(ixG^T,ep1_) = mygridvars(igrid)%w(ixG^T,bp2_) * w(ixG^T,v3_) &
           - mygridvars(igrid)%w(ixG^T,bp3_) * w(ixG^T,v2_) + eqpar(eta_) * current(ixG^T,1)

      mygridvars(igrid)%w(ixG^T,ep2_) = mygridvars(igrid)%w(ixG^T,bp3_) * w(ixG^T,v1_) &
           - mygridvars(igrid)%w(ixG^T,bp1_) * w(ixG^T,v3_) + eqpar(eta_) * current(ixG^T,2)

      mygridvars(igrid)%w(ixG^T,ep3_) = mygridvars(igrid)%w(ixG^T,bp1_) * w(ixG^T,v2_) &
           - mygridvars(igrid)%w(ixG^T,bp2_) * w(ixG^T,v1_) + eqpar(eta_) * current(ixG^T,3)
   case default
      call mpistop('typeaxial not implemented')
   end select

! scale to cgs units:
   mygridvars(igrid)%w(ixG^T,bp1_:bp^NC_) = &
        mygridvars(igrid)%w(ixG^T,bp1_:bp^NC_) * sqrt(4.0d0*dpi*UNIT_VELOCITY**2.0d0 * UNIT_DENSITY)
   mygridvars(igrid)%w(ixG^T,ep1_:ep^NC_) = &
        mygridvars(igrid)%w(ixG^T,ep1_:ep^NC_) * sqrt(4.0d0*dpi*UNIT_VELOCITY**2.0d0 * UNIT_DENSITY) * UNIT_VELOCITY / CONST_c

}

{#IFDEF PARTICLES_GCA
   
   ! fill with magnetic field:
   mygridvars(igrid)%w(ixG^T,bp1_:bp^NC_) = w(ixG^T,b1_:b^NC_)

! fill with electric field
   current = zero
   call getcurrent(w,ixG^LL,ixG^LLIM^D^LSUB1,idirmin,current)
   select case (typeaxial)
   case ('slab')
      mygridvars(igrid)%w(ixG^T,ep1_) = mygridvars(igrid)%w(ixG^T,bp2_) * w(ixG^T,v3_) &
           - mygridvars(igrid)%w(ixG^T,bp3_) * w(ixG^T,v2_) + eqpar(eta_) * current(ixG^T,1)

      mygridvars(igrid)%w(ixG^T,ep2_) = mygridvars(igrid)%w(ixG^T,bp3_) * w(ixG^T,v1_) &
           - mygridvars(igrid)%w(ixG^T,bp1_) * w(ixG^T,v3_) + eqpar(eta_) * current(ixG^T,2)

      mygridvars(igrid)%w(ixG^T,ep3_) = mygridvars(igrid)%w(ixG^T,bp1_) * w(ixG^T,v2_) &
           - mygridvars(igrid)%w(ixG^T,bp2_) * w(ixG^T,v1_) + eqpar(eta_) * current(ixG^T,3)
   case ('cylindrical')
      call mpistop('cylindrical not implemented')
   end select

! scale to cgs units:
   mygridvars(igrid)%w(ixG^T,bp1_:bp^NC_) = &
        mygridvars(igrid)%w(ixG^T,bp1_:bp^NC_) * sqrt(4.0d0*dpi*UNIT_VELOCITY**2.0d0 * UNIT_DENSITY)
   mygridvars(igrid)%w(ixG^T,ep1_:ep^NC_) = &
        mygridvars(igrid)%w(ixG^T,ep1_:ep^NC_) * sqrt(4.0d0*dpi*UNIT_VELOCITY**2.0d0 * UNIT_DENSITY) * UNIT_VELOCITY / CONST_c




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
      call gradient(kappa_B,ixG^LL,ixG^LL^LSUB1,idir,tmp)
      mygridvars(igrid)%w(ixG^T,grad_kappa_B1_-1+idir) = tmp(ixG^T)/UNIT_LENGTH
   end do

   {^C& bhat(ixG^T,^C) = mygridvars(igrid)%w(ixG^T,bp^C_) / absB(ixG^T)\}

   ! (b dot grad) b and the other directional derivatives
   
   {^C&
   do idir=1,ndim

      call gradient(bhat(ixG^T,^C),ixG^LL,ixG^LL^LSUB1,idir,tmp)
      mygridvars(igrid)%w(ixG^T,b_dot_grad_b^C_) = mygridvars(igrid)%w(ixG^T,b_dot_grad_b^C_) &
           + bhat(ixG^T,idir) * tmp(ixG^T)/UNIT_LENGTH
      mygridvars(igrid)%w(ixG^T,ue_dot_grad_b^C_) = mygridvars(igrid)%w(ixG^T,ue_dot_grad_b^C_) &
           + ue(ixG^T,idir) * tmp(ixG^T)/UNIT_LENGTH

      call gradient(ue(ixG^T,^C),ixG^LL,ixG^LL^LSUB1,idir,tmp)
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
