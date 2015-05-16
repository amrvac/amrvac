subroutine hlld(method,qdt,ixI^L,ixO^L,idim^LIM, &
                     qtC,wCT,qt,wnew,wold,fC,dx^D,x)

! method=='hlld'  --> 2nd order HLLD scheme.
! method=='hlld1' --> 1st order HLLD scheme.
! method=='hlldd' --> 2nd order HLLD+tvdlf scheme.
! method=='hlldd1'--> 1st order HLLD+tvdlf scheme.

include 'amrvacdef.f'

character(len=*), intent(in)                         :: method
double precision, intent(in)                         :: qdt, qtC, qt, dx^D
integer, intent(in)                                  :: ixI^L, ixO^L, idim^LIM
double precision, dimension(ixI^S,1:ndim), intent(in) ::  x
double precision, dimension(ixI^S,1:ndim)             ::  xi
double precision, dimension(ixI^S,1:nw)               :: wCT, wnew, wold
double precision, dimension(ixI^S,1:nwflux,1:ndim)  :: fC

double precision, dimension(ixG^T,1:nw)            :: wLC, wRC
double precision, dimension(ixG^T)                 :: vLC, vRC
double precision, dimension(ixG^T)                 :: cmaxC,cminC

double precision, dimension(1:ndim)                :: dxinv
double precision                                   :: dxdim

integer, dimension(ixG^T)                          :: patchf
integer :: idims, iw, ix^L, hxO^L, ixC^L, jxC^L, kxC^L, kxR^L
logical :: transport, new_cmax, CmaxMeanState, logiB
logical, dimension(ixG^T) :: patchw

!-----------------------------------------------------------------------------

call mpistop("hlld not yet implemented")

end subroutine hlld
