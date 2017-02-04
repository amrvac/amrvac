!============================================================================
! Author: Allard Jan van Marle, version march 2011
!============================================================================
subroutine pointgrav(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

!
!  Calculate velocity change due to point source gravity
!  As it stands now, the points source has to be in the same
!  'dimensionality' as the grid. E.g. It is not possible to 
!  have an xy grid, with the point source localized in the 
!  z-direction.
!  On the other hand, it IS possible to have the point source 
!  in a location that is not part of the grid. This will in 
!  fact often be the case in spherical symmetries with the 
!  gravitational source at r=0
!
!  the scaled value for the central mass: ptmass has to be set in
!  the parfile, together with the coordinates x1ptms,x2ptms,x3ptms
!


use mod_global_parameters


integer, intent(in)             :: ixI^L, ixO^L, iw^LIM
double precision, intent(in)    :: qdt, qtC, qt, x(ixI^S,1:ndim)
double precision, intent(inout) :: wCT(ixI^S,1:nw), w(ixI^S,1:nw)

double precision                :: r2i(ixG^T,1:ndim)
!-----------------------------------------------------------------------------

if (ptmass == zero) return


call getgvector(ixI^L,ixO^L,x,r2i)

w(ixO^S,e_) =  w(ixO^S,e_)    &
            - qdt *ptmass*wCT(ixO^S,s1_)*r2i(ixO^S,1)   
{^NOONED
w(ixO^S,e_) =  w(ixO^S,e_)    &
            - qdt *ptmass*wCT(ixO^S,s2_)*r2i(ixO^S,2)
}
{^IFTHREED
w(ixO^S,e_) =  w(ixO^S,e_)  &
             - qdt *ptmass*wCT(ixO^S,s3_)*r2i(ixO^S,3)
}

!
!  Update momentum
!

w(ixO^S,s1_) =  w(ixO^S,s1_)  &
             - qdt *ptmass*wCT(ixO^S,rho_)*r2i(ixO^S,1)
{^NOONED
w(ixO^S,s2_) =  w(ixO^S,s2_)  &
             - qdt *ptmass*wCT(ixO^S,rho_)*r2i(ixO^S,2)
}
{^IFTHREED
w(ixO^S,s3_) =  w(ixO^S,s3_)  &
             - qdt *ptmass*wCT(ixO^S,rho_)*r2i(ixO^S,3)
}

return 

end subroutine pointgrav
!==============================================================================
subroutine getdt_pointgrav(w,ixG^L,ix^L,dtnew,dx^D,x)

!
! Limits timestep for gravitational pointsource
!


use mod_global_parameters

integer, intent(in)             :: ixG^L, ix^L
double precision, intent(in)    :: dx^D, x(ixG^S,1:ndim)
double precision, intent(inout) :: w(ixG^S,1:nw), dtnew

double precision                :: r2i(ixG^T,1:ndim)

double precision                :: dtgrav(1:ndim), tmp(ixG^T), dxinv(1:ndim)
integer                         :: idims
!-----------------------------------------------------------------------------

dtgrav = bigdouble
if( ptmass == zero ) return

^D&dxinv(^D)=one/dx^D;

call getgvector(ixG^L,ix^L,x,r2i)

r2i(ix^S,1:ndim)=dabs(r2i(ix^S,1:ndim))

do idims=1,ndim
   if (.not.slab) then
      tmp(ix^S)     = mygeo%dx(ix^S,idims)/(min(smalldouble, ptmass*r2i(ix^S,idims)))
      dtgrav(idims) = minval(tmp(ix^S))
   else
      dtgrav(idims) = minval(one/(ptmass*r2i(ix^S,idims)*dxinv(idims)))
   end if
enddo

dtnew = sqrt(minval(dtgrav(1:ndim)))

return
end subroutine getdt_pointgrav
!===========================================================================
subroutine getgvector(ixI^L,ixO^L,x,r2i)

!
!  returns the vector components of the 1/r^2 vector:
!  dx1/r^3,    
!  dx2/r^3,  
!  dx3/r^3,
!  between a pre-determined pointsource and a point in the grid.
!


use mod_global_parameters

integer, intent(in)           :: ixO^L,ixI^L
double precision, intent(in)  :: x(ixI^S,1:ndim)
double precision, intent(out) :: r2i(ixG^T,1:ndim)

double precision              :: rstptms, rctptms

double precision, allocatable :: r2(:^D&)
{^NOONED
double precision, allocatable :: rst(:^D&), rct(:^D&)
}
!-----------------------------------------------------------------------------

allocate(r2(ixO^S))


select case (typeaxial)
case('slab')
{^IFONED   
   r2(ixO^S) = (x(ixO^S,1)-x1ptms)**2
}      
{^IFTWOD
   r2(ixO^S) = (x(ixO^S,1)-x1ptms)**2  &
             + (x(ixO^S,2)-x2ptms)**2
}
{^IFTHREED
   r2(ixO^S) = (x(ixO^S,1)-x1ptms)**2  &
             + (x(ixO^S,2)-x2ptms)**2  &
             + (x(ixO^S,3)-x3ptms)**2
} 


   r2i(ixO^S,1) = (x(ixO^S,1)-x1ptms) &
                / (r2(ixO^S) * sqrt(r2(ixO^S)) + smalldouble)
{^NOONED
   r2i(ixO^S,2) = (x(ixO^S,2)-x2ptms) &
                / (r2(ixO^S) * sqrt(r2(ixO^S)) + smalldouble)
}
{^IFTHREED		 
   r2i(ixO^S,3) = (x(ixO^S,3)-x3ptms) &
                / (r2(ixO^S) * sqrt(r2(ixO^S)) + smalldouble)
}		 

case('cylindrical')
{^IFONED   
   r2(ixO^S)    = (x(ixO^S,r_)-x1ptms)**2   
   r2i(ixO^S,1) = one / (r2(ixO^S) + smalldouble)
}
{^IFTWOD
   if(^PHI == 2) then
      r2(ixO^S)    = (x(ixO^S,r_)-x1ptms)**2 &
                   + two*x(ixO^S,r_)*x1ptms  &
	           * (one -cos(x(ixO^S,2)-x2ptms))

      r2i(ixO^S,1) = (x(ixO^S,r_)-x1ptms*cos(x(ixO^S,2)-x2ptms)) &
                   / (r2(ixO^S) * sqrt(r2(ixO^S)) + smalldouble)   
 
      r2i(ixO^S,2) = x1ptms * sin(x(ixO^S,2)-x2ptms)             &
                   / (r2(ixO^S) * sqrt(r2(ixO^S)) + smalldouble)      		       
   else
      r2(ixO^S)    = (x(ixO^S,r_)-x1ptms)**2 &
                   + (x(ixO^S,2)-x2ptms)**2   

      r2i(ixO^S,1) = (x(ixO^S,r_)-x1ptms)    &
                   / (r2(ixO^S) * sqrt(r2(ixO^S)) + smalldouble)
      r2i(ixO^S,2) = (x(ixO^S,2)-x2ptms)    &
                   / (r2(ixO^S) * sqrt(r2(ixO^S)) + smalldouble)      
   endif 
}
{^IFTHREED
   r2(ixO^S) = (x(ixO^S,r_)-x1ptms)**2  &
             + (x(ixO^S,z_)-x2ptms)**2  &  
             + two*x(ixO^S,r_)*x1ptms   &
	     * (one -cos(x(ixO^S,phi_)-x3ptms))

   r2i(ixO^S,1) = (x(ixO^S,r_)-x1ptms*cos(x(ixO^S,phi_)-x3ptms)) &
                / (r2(ixO^S) * sqrt(r2(ixO^S)) + smalldouble)   
   r2i(ixO^S,2) = (x(ixO^S,z_)-x2ptms)                           &
                / (r2(ixO^S) * sqrt(r2(ixO^S)) + smalldouble)   
   r2i(ixO^S,3) = x1ptms * sin(x(ixO^S,phi_)-x3ptms)             &
                / (r2(ixO^S) * sqrt(r2(ixO^S)) + smalldouble)      		       
}


case('spherical')
   rstptms = x1ptms * sin(x2ptms)
   rctptms = x1ptms * cos(x2ptms)
   
{^IFONED   
   r2(ixO^S)   = (x(ixO^S,r_)-x1ptms)**2   
   r2i(ixO^S,1) = one / (r2(ixO^S) + smalldouble)
}
{^IFTWOD
   allocate(rct(ixO^S))
   allocate(rst(ixO^S))
      
   rct(ixO^S) = x(ixO^S,r_)*cos(x(ixO^S,2))
   rst(ixO^S) = x(ixO^S,r_)*sin(x(ixO^S,2))

   r2(ixO^S) = (rct(ixO^S)-rctptms)**2  &
             + (rst(ixO^S)-rstptms)**2

 
   r2i(ixO^S,1) = ((rct(ixO^S)-rctptms) * cos(x(ixO^S,2))     &
                + (rst(ixO^S)-rstptms)  * sin(x(ixO^S,2)))    &
                / (r2(ixO^S) * sqrt(r2(ixO^S)) + smalldouble)
   r2i(ixO^S,2) = (-(rct(ixO^S)-rctptms) * sin(x(ixO^S,2))    & 
                + (rst(ixO^S) - rstptms) * cos(x(ixO^S,2)))   & 
                / (r2(ixO^S) * sqrt(r2(ixO^S)) + smalldouble)
  
   deallocate(rst)
   deallocate(rct) 
}   
{^IFTHREED
   allocate(rct(ixO^S))
   allocate(rst(ixO^S))

   rct(ixO^S) = x(ixO^S,r_)*cos(x(ixO^S,2))
   rst(ixO^S) = x(ixO^S,r_)*sin(x(ixO^S,2))

   r2(ixO^S) = (rct(ixO^S)-rctptms)**2    &
             + (rst(ixO^S)-rstptms)**2    &
	     + two * rst(ixO^S) * rstptms &
	     * (one-cos(x(ixO^S,phi_)-x3ptms))

   r2i(ixO^S,1) = ((rct(ixO^S)-rctptms) * cos(x(ixO^S,2))        &
                + (rst(ixO^S)-rstptms*cos(x(ixO^S,phi_)-x3ptms)) &
     	        * sin(x(ixO^S,2)))    &
                / (r2(ixO^S) * sqrt(r2(ixO^S)) + smalldouble)
   r2i(ixO^S,2) = (-(rct(ixO^S)-rctptms) * sin(x(ixO^S,2))       & 
                + (rst(ixO^S)-rstptms*cos(x(ixO^S,phi_)-x3ptms)) &
	        * cos(x(ixO^S,2)))  & 
                / (r2(ixO^S) * sqrt(r2(ixO^S)) + smalldouble)    
   r2i(ixO^S,3) = rstptms * sin(x(ixO^S,phi_)-x3ptms)            &
                / (r2(ixO^S) * sqrt(r2(ixO^S)) + smalldouble)
 

   deallocate(rct)
   deallocate(rst)
}
   end select 

deallocate(r2)

return

end subroutine getgvector
!===========================================================================
