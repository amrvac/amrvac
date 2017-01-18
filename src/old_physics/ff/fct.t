{#IFDEF FCT
!=============================================================================
subroutine fct_average(ixI^L, ixO^L, fC)

! Performs the average for the flux constrained transport from Toth 2000. 
! His equation (25)

use mod_global_parameters

integer, intent(in)                :: ixI^L, ixO^L
double precision, intent(inout)    :: fC(ixI^S,1:nwflux,1:ndim)

double precision                   :: fB(ixG^T,1:ndir,1:ndim)
integer                            :: iwdim, iwdir, idim, idir
integer                            :: ixJp^L, ixKp^L, ixKm^L, ixC^L, ixJpKm^L
!-----------------------------------------------------------------------------


do idim = 1, ndim ! idim is the direction of the flux interface
   do idir = 1, ndim ! idir is the component of the field to consider
      if (idim.eq.idir) cycle
      iwdir = b0_+idir; iwdim = b0_+idim

! Assemble indices:
      {ixCmax^D=ixOmax^D;}
      {ixCmin^D=ixOmin^D-kr(idim,^D);} ! Extend range in flux direction
 
      {ixJp^L=ixC^L+kr(idim,^D);}
      {ixKp^L=ixC^L+kr(idir,^D);}
      {ixKm^L=ixC^L-kr(idir,^D);}
      {ixJpKm^L=ixJp^L-kr(idir,^D);}


! Perform flux average:
         fB(ixC^S,idir,idim) = &
              0.125d0 * ( 2.0d0* fC(ixC^S,iwdir,idim) + fC(ixKp^S,iwdir,idim) &
              + fC(ixKm^S,iwdir,idim) &
              - fC(ixC^S,iwdim,idir)  &
              - fC(ixJp^S,iwdim,idir)  &
              - fC(ixKm^S,iwdim,idir)  &
              - fC(ixJpKm^S,iwdim,idir)  &
              )

   end do ! idir
end do ! idim

! Overwrite the new flux entries:
do idim = 1, ndim ! idim is the direction of the flux interface
   do idir = 1, ndim ! idir is the component of the field to consider
      iwdir = b0_+idir
      if (idim.eq.idir) then
         ! To conserve divb to machine precision, this needs to be zero.
         ! However, since restriction and prolongation can introduce divb
         ! when using AMR, additional measures for divb-control must be taken. 
         ! These rely on the normal field flux component (e.g. as in GLM).
         ! Thus when more than one level is used, we dont set this zero 
         ! to allow additional divb-control to work.
         if (mxnest .eq. 1) fC(ixI^S,iwdir,idim) = zero
         cycle
      end if

      fC(ixI^S,iwdir,idim) = fB(ixI^S,idir,idim)

   end do ! idir
end do ! idim

end subroutine fct_average
!=============================================================================
subroutine b_from_vectorpotential(ixI^L, ixO^L, idirmin, idirmax, w, x)

! Implemented for Cartesian coordinates

use mod_global_parameters

integer, intent(in)                :: ixI^L, ixO^L, idirmin, idirmax
double precision, intent(inout)    :: w(ixI^S,1:nw)
double precision, intent(in)       :: x(ixI^S,1:ndim)

integer                            :: ixC^L, ixCp^L, ixCm^L, ixOm^L, idim, idir, j, k
double precision                   :: xC(ixG^T,1:ndim), A(ixG^T,1:ndir), dxC(ixG^T,1:ndim)
double precision                   :: B(ixG^T,idirmin:idirmax)
!-----------------------------------------------------------------------------

! get the corner-coordinates:

{ixCmax^D=ixOmax^D;}
{ixCmin^D=ixOmin^D-1;} ! Extend range by one

do idim = 1, ndim
   {ixCp^L=ixC^L+kr(idim,^D);}
   {ixCm^L=ixC^L-kr(idim,^D);}
   xC(ixC^S,idim) = half * (x(ixCp^S,idim) + x(ixC^S,idim))
   dxC(ixC^S,idim) = xC(ixC^S,idim) - xC(ixCm^S,idim)
end do

! initialize the vectorpotential:

call initvecpot_usr(ixG^LL, ixC^L, xC, A)


! take the curl of the vectorpotential: 

B(:^D&,:) = zero
do idir = idirmin, idirmax
   do j = 1, ndim
      if (idir .eq. j) cycle
      {ixCmin^D=ixOmin^D-kr(idir,^D);}
      {ixCm^L=ixC^L-kr(j,^D);}
      do k = 1,ndir
         B(ixC^S,idir) = B(ixC^S,idir) &
              + lvc(idir,j,k) * (A(ixC^S,k)-A(ixCm^S,k))/dxC(ixC^S,j)
      end do
   end do
end do



! Average to the cell centers and fill solution array:

do idir = idirmin, idirmax
   {ixOm^L=ixO^L-kr(idir,^D);}
   w(ixO^S,b0_+idir) = half * (B(ixO^S,idir) + B(ixOm^S,idir))
end do

end subroutine b_from_vectorpotential
!=============================================================================
}
